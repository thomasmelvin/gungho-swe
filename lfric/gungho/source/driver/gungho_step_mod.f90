!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Steps the gungho app through one timestep

!> @details Handles the stepping (for a single timestep) of the
!>          gungho app

module gungho_step_mod

  use clock_mod,                      only : clock_type
  use conservation_algorithm_mod,     only : conservation_algorithm
  use constants_mod,                  only : i_def, r_def, l_def
  use field_collection_mod,           only : field_collection_type
  use field_mod,                      only : field_type
  use formulation_config_mod,         only : use_physics,             &
                                             moisture_formulation,    &
                                             moisture_formulation_dry
  use geometric_constants_mod,        only : get_da_at_w2
  use gungho_model_data_mod,          only : model_data_type
  use io_config_mod,                  only : write_conservation_diag, &
                                             write_minmax_tseries
  use log_mod,                        only : LOG_LEVEL_INFO
  use log_mod,                        only : log_event, &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_TRACE

  use mesh_mod,                       only : mesh_type
  use minmax_tseries_mod,             only : minmax_tseries
  use mr_indices_mod,                 only : nummr
  use semi_implicit_timestep_alg_mod, only : semi_implicit_alg_step
  use moist_dyn_mod,                  only : num_moist_factors
  use moisture_conservation_alg_mod,  only : moisture_conservation_alg
  use moisture_fluxes_alg_mod,        only : moisture_fluxes_alg
  use rk_alg_timestep_mod,            only : rk_alg_step
  use timestepping_config_mod,        only : method, &
                                             method_semi_implicit, &
                                             method_rk,            &
                                             method_no_timestepping

  implicit none

  private
  public gungho_step

  contains

  !> @brief Steps the gungho app through one timestep
  !> @param[in] mesh      The primary mesh
  !> @param[in] twod_mesh The two-dimensional mesh
  !> @param[inout] model_data The working data set for the model run
  !> @param[in] timestep number of current timestep
  subroutine gungho_step( mesh,       &
                          twod_mesh,  &
                          model_data, &
                          clock )

    implicit none

    type(mesh_type), intent(in), pointer    :: mesh
    type(mesh_type), intent(in), pointer    :: twod_mesh
    type( model_data_type ), target, intent(inout) :: model_data
    class(clock_type),               intent(in)    :: clock

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_type ),            pointer :: moist_dyn(:) => null()
    type( field_collection_type ), pointer :: adv_fields_all_outer => null()
    type( field_collection_type ), pointer :: adv_fields_last_outer => null()
    type( field_collection_type ), pointer :: derived_fields => null()
    type( field_collection_type ), pointer :: radiation_fields => null()
    type( field_collection_type ), pointer :: microphysics_fields => null()
    type( field_collection_type ), pointer :: orography_fields => null()
    type( field_collection_type ), pointer :: turbulence_fields => null()
    type( field_collection_type ), pointer :: convection_fields => null()
    type( field_collection_type ), pointer :: cloud_fields => null()
    type( field_collection_type ), pointer :: surface_fields => null()
    type( field_collection_type ), pointer :: soil_fields => null()
    type( field_collection_type ), pointer :: snow_fields => null()
    type( field_collection_type ), pointer :: chemistry_fields => null()
    type( field_collection_type ), pointer :: aerosol_fields => null()
    type( field_collection_type ), pointer :: lbc_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: dA => null()  ! Areas of faces

    real(r_def) :: dt

    logical(l_def) :: use_moisture

    write( log_scratch_space, '("/", A, "\ ")' ) repeat( "*", 76 )
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, &
           '(A,I0)' ) 'Start of timestep ', clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    use_moisture = ( moisture_formulation /= moisture_formulation_dry )

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    mr => model_data%mr
    moist_dyn => model_data%moist_dyn
    adv_fields_all_outer => model_data%adv_fields_all_outer
    adv_fields_last_outer => model_data%adv_fields_last_outer
    derived_fields => model_data%derived_fields
    radiation_fields => model_data%radiation_fields
    microphysics_fields => model_data%microphysics_fields
    orography_fields => model_data%orography_fields
    turbulence_fields => model_data%turbulence_fields
    convection_fields => model_data%convection_fields
    cloud_fields => model_data%cloud_fields
    surface_fields => model_data%surface_fields
    soil_fields => model_data%soil_fields
    snow_fields => model_data%snow_fields
    chemistry_fields => model_data%chemistry_fields
    aerosol_fields => model_data%aerosol_fields
    lbc_fields => model_data%lbc_fields

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')
    dA => get_da_at_w2(mesh%get_id())

    ! Get timestep parameters from clock
    dt = real(clock%get_seconds_per_step(), r_def)

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        call semi_implicit_alg_step(u, rho, theta, exner, mr, moist_dyn,       &
                                    adv_fields_all_outer,                      &
                                    adv_fields_last_outer,                     &
                                    derived_fields, radiation_fields,          &
                                    microphysics_fields, orography_fields,     &
                                    turbulence_fields, convection_fields,      &
                                    cloud_fields, surface_fields,              &
                                    soil_fields, snow_fields,                  &
                                    chemistry_fields, aerosol_fields,          &
                                    lbc_fields, clock, mesh,                   &
                                    twod_mesh)
      case( method_rk )             ! RK
        call rk_alg_step(u, rho, theta, moist_dyn, exner, mr, dt )
      case( method_no_timestepping )
        write( log_scratch_space, &
           '(A, A)' ) 'CAUTION: Running with no timestepping. ' // &
                      ' Prognostic fields not evolved'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

    end select

    if ( write_conservation_diag ) then
      call conservation_algorithm( rho,              &
                                   u,                &
                                   theta,            &
                                   exner )
      if ( use_moisture ) then
        call moisture_conservation_alg( rho,              &
                                        mr,               &
                                        'After timestep' )
        if ( use_physics ) call moisture_fluxes_alg( microphysics_fields, &
                                                     convection_fields,   &
                                                     turbulence_fields,   &
                                                     dA,                  &
                                                     dt )
      end if

      if (write_minmax_tseries) call minmax_tseries(u, 'u', mesh)

      call u%log_minmax(LOG_LEVEL_INFO, ' u')
      call theta%log_minmax(LOG_LEVEL_INFO, 'theta')

    end if

    write( log_scratch_space, &
           '(A,I0)' ) 'End of timestep ', clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '("\", A, "/ ")' ) repeat( "*", 76 )
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine gungho_step

end module gungho_step_mod
