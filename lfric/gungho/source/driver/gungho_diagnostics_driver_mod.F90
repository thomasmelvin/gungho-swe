!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs diagnostics from gungho/lfric_atm

!> @details Calls the routine that generates diagnostic output for
!>          gungho/lfric_atm. This is only a temporary
!>          hard-coded solution in lieu of a proper dianostic system

module gungho_diagnostics_driver_mod

  use clock_mod,                 only : clock_type
  use constants_mod,             only : i_def, str_def
  use boundaries_config_mod,     only : limited_area, output_lbcs
  use diagnostics_io_mod,        only : write_scalar_diagnostic, &
                                        write_vector_diagnostic
  use diagnostics_calc_mod,      only : write_divergence_diagnostic, &
                                        write_hydbal_diagnostic, &
                                        write_vorticity_diagnostic
  use field_collection_iterator_mod, &
                                 only : field_collection_iterator_type
  use field_collection_mod,      only : field_collection_type
  use diagnostic_alg_mod,        only : column_total_diagnostics_alg, &
                                        calc_wbig_diagnostic_alg, &
                                        pressure_diag_alg
  use gungho_model_data_mod,     only : model_data_type
  use field_mod,                 only : field_type
  use field_parent_mod,          only : field_parent_type, write_interface
  use lfric_xios_write_mod,      only : write_field_edge
  use formulation_config_mod,    only : use_physics,             &
                                        moisture_formulation,    &
                                        moisture_formulation_dry
  use fs_continuity_mod,         only : W3, Wtheta
  use integer_field_mod,         only : integer_field_type
  use initialization_config_mod, only : ls_option,          &
                                        ls_option_analytic, &
                                        ls_option_file
  use moist_dyn_mod,             only : num_moist_factors
  use mr_indices_mod,            only : nummr, mr_names
  use log_mod,                   only : log_event, &
                                        LOG_LEVEL_INFO
  use mesh_mod,                  only : mesh_type
  use geometric_constants_mod,   only : get_panel_id, get_height
  use io_config_mod,             only: subroutine_timers, use_xios_io, write_fluxes
  use timer_mod,                 only: timer


#ifdef UM_PHYSICS
  use pmsl_alg_mod,              only : pmsl_alg
#endif

  implicit none

  private
  public gungho_diagnostics_driver

contains

  !> @brief Outputs simple diagnostics from Gungho/LFRic
  !> @param[in] mesh       The primary mesh
  !> @param[in] twod_mesh  The 2d mesh
  !> @param[in] model_data The working data set for the model run
  !> @param[in] timestep The timestep at which the fields are valid
  !> @param[in] nodal_output_on_w3 Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  subroutine gungho_diagnostics_driver( mesh,       &
                                        twod_mesh,  &
                                        model_data, &
                                        clock,      &
                                        nodal_output_on_w3 )

    implicit none

    type(mesh_type),       intent(in), pointer :: mesh
    type(mesh_type),       intent(in), pointer :: twod_mesh
    type(model_data_type), intent(in), target  :: model_data
    class(clock_type),     intent(in)          :: clock
    logical,               intent(in)          :: nodal_output_on_w3

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_collection_type ), pointer :: lbc_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_type ),            pointer :: moist_dyn(:) => null()
    type( field_collection_type ), pointer :: derived_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: h_u => null()
    type( field_type), pointer :: v_u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: panel_id => null()
    type( field_type), pointer :: height_w3 => null()
    type( field_type), pointer :: height_wth => null()
    type( field_type), pointer :: exner_in_wth => null()
    type( field_type), pointer :: lbc_u => null()
    type( field_type), pointer :: lbc_theta => null()
    type( field_type), pointer :: lbc_rho => null()
    type( field_type), pointer :: lbc_exner => null()
    type( field_type), pointer :: lbc_m_v=> null()
    type( field_type), pointer :: lbc_q=> null()
    type( field_type), pointer :: theta_in_w3 => null()
    type( field_type), pointer :: u_in_w2h => null()
    type( field_type), pointer :: v_in_w2h => null()
    type( field_type), pointer :: w_in_wth => null()

    ! Iterator for field collection
    type(field_collection_iterator_type)  :: iterator

    ! A pointer used for retrieving fields from collections
    ! when iterating over them
    class( field_parent_type ), pointer :: field_ptr  => null()

    procedure(write_interface), pointer  :: tmp_write_ptr => null()

    character(str_def) :: name

    integer :: i, fs

    if ( subroutine_timers ) call timer('gungho_diagnostics_driver')

    call log_event("Gungho: writing diagnostic output", LOG_LEVEL_INFO)

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    lbc_fields => model_data%lbc_fields
    mr => model_data%mr
    moist_dyn => model_data%moist_dyn
    derived_fields => model_data%derived_fields
    panel_id => get_panel_id(mesh%get_id())
    height_w3 => get_height(W3, mesh%get_id())
    height_wth => get_height(Wtheta, mesh%get_id())

    ! Can't just iterate through the prognostic/diagnostic collections as
    ! some fields are scalars and some fields are vectors, so explicitly
    ! extract all fields from the collections and output each of them in turn
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    ! Scalar fields
    call write_scalar_diagnostic('rho', rho, &
                                 clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('theta', theta, &
                                 clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('exner', exner, &
                                 clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('height_w3', height_w3, &
                                 clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('height_wth', height_wth, &
                                 clock, mesh, nodal_output_on_w3)

    ! Vector fields
    if (use_physics .and. use_xios_io .and. .not. write_fluxes) then
      ! These have already been calculated, so no need to recalculate them
      u_in_w2h => derived_fields%get_field('u_in_w2h')
      v_in_w2h => derived_fields%get_field('v_in_w2h')
      w_in_wth => derived_fields%get_field('w_in_wth')
      tmp_write_ptr => write_field_edge
      call u_in_w2h%set_write_behaviour(tmp_write_ptr)
      call v_in_w2h%set_write_behaviour(tmp_write_ptr)
      if (clock%is_initialisation()) then
        call u_in_w2h%write_field("init_u_in_w2h")
        call v_in_w2h%write_field("init_v_in_w2h")
        call w_in_wth%write_field("init_w_in_wth")
      else
        call u_in_w2h%write_field("u_in_w2h")
        call v_in_w2h%write_field("v_in_w2h")
        call w_in_wth%write_field("w_in_wth")
      end if
    else
      call write_vector_diagnostic('u', u, &
                                 clock, mesh, nodal_output_on_w3)
    end if
    call write_vorticity_diagnostic( u, clock )

    ! Moisture fields
    if ( moisture_formulation /= moisture_formulation_dry ) then
      do i=1,nummr
        call write_scalar_diagnostic( trim(mr_names(i)), mr(i), &
                                      clock, mesh, nodal_output_on_w3 )
      end do
    end if

    if (limited_area) then
      if (output_lbcs) then
        lbc_theta => lbc_fields%get_field('lbc_theta')
        lbc_u => lbc_fields%get_field('lbc_u')
        lbc_rho => lbc_fields%get_field('lbc_rho')
        lbc_exner => lbc_fields%get_field('lbc_exner')

        h_u => lbc_fields%get_field('lbc_h_u')
        v_u => lbc_fields%get_field('lbc_v_u')

        ! Scalar fields
        call write_scalar_diagnostic('lbc_rho', lbc_rho, &
                                 clock, mesh, nodal_output_on_w3)
        call write_scalar_diagnostic('lbc_theta', lbc_theta, &
                                 clock, mesh, nodal_output_on_w3)
        call write_scalar_diagnostic('lbc_exner', lbc_exner, &
                                 clock, mesh, nodal_output_on_w3)
        call write_scalar_diagnostic('readlbc_v_u', v_u, &
                                 clock, mesh, nodal_output_on_w3)

        if ( moisture_formulation /= moisture_formulation_dry ) then
          lbc_m_v => lbc_fields%get_field('lbc_m_v')
          call write_scalar_diagnostic('lbc_m_v', lbc_m_v, &
                                   clock, mesh, nodal_output_on_w3)
          lbc_q => lbc_fields%get_field('lbc_q')
          call write_scalar_diagnostic('lbc_q', lbc_q, &
                                   clock, mesh, nodal_output_on_w3)
          lbc_rho => lbc_fields%get_field('lbc_rho_r2')
          call write_scalar_diagnostic('lbc_rho_r2', lbc_rho, &
                                   clock, mesh, nodal_output_on_w3)
        end if

        ! Vector fields
        call write_vector_diagnostic('lbc_u', lbc_u, &
                                 clock, mesh, nodal_output_on_w3)
        call write_vector_diagnostic('readlbc_h_u', h_u, &
                                 clock, mesh, nodal_output_on_w3)
      endif
    endif

    ! Derived physics fields (only those on W3 or Wtheta)
    if (use_physics .and. .not. clock%is_initialisation()) then

      call iterator%initialise(derived_fields)
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        select type(field_ptr)
          type is (field_type)
            fs = field_ptr%which_function_space()
            if ( fs == W3 .or. fs == Wtheta ) then
              name = trim(adjustl( field_ptr%get_name() ))
              call write_scalar_diagnostic( trim(name), field_ptr, &
                                            clock,                 &
                                            mesh, nodal_output_on_w3 )
            end if
        end select
      end do
      field_ptr => null()

      ! Get w_in_wth for WBig calculation
      w_in_wth => derived_fields%get_field('w_in_wth')
      call calc_wbig_diagnostic_alg(w_in_wth, mesh)

      ! Pressure diagnostics
      exner => prognostic_fields%get_field('exner')
      call pressure_diag_alg(exner)

      exner_in_wth => derived_fields%get_field('exner_in_wth')
      call pressure_diag_alg(exner_in_wth)

#ifdef UM_PHYSICS
      ! Call PMSL algorithm
      theta => prognostic_fields%get_field('theta')
      call pmsl_alg(exner, derived_fields, theta, twod_mesh)
#endif
      theta_in_w3 => derived_fields%get_field('theta_in_w3')
      call column_total_diagnostics_alg(rho, mr, theta_in_w3, exner, mesh, twod_mesh)

    end if

    if (ls_option /= ls_option_file .and. ls_option /= ls_option_analytic) then
      ! Other derived diagnostics with special pre-processing
      ! Don't output for the tangent linear model
      call write_divergence_diagnostic( u, clock, mesh )
      call write_hydbal_diagnostic( theta, moist_dyn, exner, mesh )
    end if

    if ( subroutine_timers ) call timer('gungho_diagnostics_driver')
  end subroutine gungho_diagnostics_driver

end module gungho_diagnostics_driver_mod
