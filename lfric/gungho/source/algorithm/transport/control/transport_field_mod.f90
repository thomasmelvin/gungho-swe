!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains central routine for transporting fields.
!> @details Contains routine to transport a (multidata) field pointing to
!!          particular routines based on the specified transport options.

module transport_field_mod

  use constants_mod,                    only: r_def
  use field_mod,                        only: field_type
  use log_mod,                          only: log_event, LOG_LEVEL_ERROR
  use ffsl_control_alg_mod,             only: ffsl_conservative_control, &
                                              ffsl_advective_control
  use split_transport_mod,              only: split_transport_control
  use transport_metadata_mod,           only: transport_metadata_type
  use transport_runtime_alg_mod,        only: transport_runtime_type
  use transport_runtime_collection_mod, only: get_transport_runtime
  use transport_enumerated_types_mod,   only: scheme_mol_3d,         &
                                              scheme_ffsl_3d,        &
                                              scheme_split,          &
                                              direction_3d,          &
                                              equation_conservative, &
                                              equation_advective,    &
                                              equation_consistent
  use mol_conservative_alg_mod,         only: mol_conservative_alg
  use mol_advective_alg_mod,            only: mol_advective_alg
  use mol_consistent_alg_mod,           only: mol_consistent_alg

  implicit none

  private

  public :: transport_field

contains

  !> @brief Central routine for transporting fields.
  !> @details Performs a whole time step, solving the transport equation for
  !!          a (multidata) field.
  !> @param[in,out] field_np1          Field to return at end of transport step
  !> @param[in]     field_n            Field at the start of the transport step
  !> @param[in]     model_dt           Time difference across time step
  !> @param[in]     transport_metadata Contains the configuration options for
  !!                                   transporting these fields
  subroutine transport_field(field_np1, field_n, model_dt, transport_metadata)

    implicit none

    ! Arguments
    type(field_type),              intent(inout) :: field_np1
    type(field_type),              intent(in)    :: field_n
    real(kind=r_def),              intent(in)    :: model_dt
    type(transport_metadata_type), intent(in)    :: transport_metadata
    type(transport_runtime_type),  pointer       :: transport_runtime => null()

    ! Reset the counter for tracer transport steps
    transport_runtime => get_transport_runtime(field_n%get_mesh())
    call transport_runtime%reset_tracer_step_ctr()
    nullify( transport_runtime )

    ! First choose scheme, and for full 3D schemes then choose equation
    select case ( transport_metadata%get_scheme() )

    ! -------------------------------------------------------------------------!
    ! Full 3D Method of Lines scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_mol_3d )
      ! Choose form of transport equation
      select case ( transport_metadata%get_equation() )
      case ( equation_conservative )
         call mol_conservative_alg(field_np1, field_n, &
                                   direction_3d, transport_metadata)

      case ( equation_advective )
         call mol_advective_alg(field_np1, field_n, &
                                direction_3d, transport_metadata)

      case ( equation_consistent )
         call mol_consistent_alg(field_np1, field_n, &
                                 direction_3d, transport_metadata)

      case default
        call log_event('Trying to solve unrecognised form of transport equation', &
                        LOG_LEVEL_ERROR)

      end select

    ! -------------------------------------------------------------------------!
    ! Full 3D Flux-Form Semi-Lagrangian scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_ffsl_3d )
      ! Choose form of transport equation
      select case ( transport_metadata%get_equation() )
      case ( equation_conservative )
        call ffsl_conservative_control(field_np1, field_n, direction_3d, &
                                       model_dt, transport_metadata)

      case ( equation_advective )
        call ffsl_advective_control(field_np1, field_n, direction_3d, &
                                    model_dt, transport_metadata)

      case default
        call log_event('Trying to solve unrecognised form of transport equation', &
                        LOG_LEVEL_ERROR)

      end select

    ! -------------------------------------------------------------------------!
    ! Some split horizontal/vertical transport scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_split )
      call split_transport_control(field_np1, field_n, model_dt, &
                                   transport_metadata)

    case default
      call log_event('Trying to transport with unrecognised scheme', &
                      LOG_LEVEL_ERROR)

    end select

  end subroutine transport_field

end module transport_field_mod
