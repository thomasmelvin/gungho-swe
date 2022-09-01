!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Steps the shallow water miniapp through one timestep.
!!
!> @details Handles the stepping (for a single timestep) of the
!!          shallow water app, using either a semi-implicit or an
!!          explicit time stepping scheme.

module shallow_water_step_mod

  use clock_mod,                      only: clock_type
  use constants_mod,                  only: i_def
  use field_mod,                      only: field_type
  use field_collection_mod,           only: field_collection_type
  use log_mod,                        only: log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_INFO,    &
                                            LOG_LEVEL_ERROR
  use swe_timestep_alg_mod,           only: swe_timestep_alg_si, &
                                            swe_timestep_alg_rk
  use shallow_water_model_data_mod,   only: model_data_type
  use shallow_water_settings_config_mod, &
                                      only: time_scheme,               &
                                            time_scheme_semi_implicit, &
                                            time_scheme_explicit

  implicit none

  private
  public shallow_water_step

  contains

  !> @brief Steps the shallow water miniapp through one timestep.
  !> @details Extracts prognostic fields from model_data then calls
  !!          swe_timestep_alg to step the shallow water miniapp
  !!          through one timestep. An iterated semi-implicit or an
  !!          explicit Runge-Kutta timestepping scheme is used.
  !> @param [in,out] model_data   The working data set for the model run
  !> @param [in]     clock        Model time
  subroutine shallow_water_step( model_data, &
                                 clock )

    implicit none

    type( model_data_type ), target, intent(inout) :: model_data
    class(clock_type),               intent(in)    :: clock

    type( field_collection_type ), pointer :: prognostic_fields => null()

    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: buoyancy => null()
    type(field_type), pointer :: q => null()

    prognostic_fields => model_data%prognostic_fields

    wind     => prognostic_fields%get_field('wind')
    geopot   => prognostic_fields%get_field('geopot')
    buoyancy => prognostic_fields%get_field('buoyancy')
    q        => prognostic_fields%get_field('q')

    write( log_scratch_space, &
           '(A,I0)' ) 'Start of timestep ', clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    select case( time_scheme )

    case( time_scheme_semi_implicit )
      call swe_timestep_alg_si(wind, geopot, buoyancy, q, model_data%s_geopot)
    case ( time_scheme_explicit )
      call swe_timestep_alg_rk(wind, geopot, buoyancy, q, model_data%s_geopot)
    case default
      call log_event("No valid time stepping scheme selected", LOG_LEVEL_ERROR)

    end select

  end subroutine shallow_water_step

end module shallow_water_step_mod
