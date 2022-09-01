!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs the intial conditions output for Gungho
!>
!>
module initial_output_mod

  use clock_mod,                  only : clock_type
  use constants_mod,              only : i_def, l_def
  use io_config_mod,              only : write_diag
  use gungho_model_data_mod,      only : model_data_type
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use lfric_xios_clock_mod,       only : lfric_xios_clock_type
  use mesh_mod,                   only : mesh_type
  use io_context_mod,             only : io_context_type
  use lfric_xios_context_mod,     only : lfric_xios_context_type
  use limited_area_constants_mod, only : write_masks
  use boundaries_config_mod,      only : limited_area, lbc_method, &
                                         lbc_method_onion_layer
  implicit none

  private
  public :: write_initial_output

contains

  !> @brief Outputs simple diagnostics from Gungho/LFRic
  !> @param[in] mesh       The primary mesh
  !> @param[in] twod_mesh  The 2d mesh
  !> @param[in] model_data The working data set for the model run
  !> @param[in] io_context The model IO context
  !> @param[in] nodal_output_on_w3 Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  subroutine write_initial_output( mesh, twod_mesh, model_data, &
                                   io_context, nodal_output_on_w3 )

    implicit none

    type(mesh_type), pointer, intent(in)    :: mesh
    type(mesh_type), pointer, intent(in)    :: twod_mesh
    type(model_data_type),    intent(in)    :: model_data
    class(io_context_type),   intent(inout) :: io_context
    logical(l_def),           intent(in)    :: nodal_output_on_w3

    class(clock_type), pointer :: clock => null()

    ! Call clock initial step before initial conditions output
    clock => io_context%get_clock()
    select type( clock )
    type is (lfric_xios_clock_type)
        call clock%initial_step()
    end select

    if (clock%is_initialisation() .and. write_diag) then
      ! Calculation and output of initial conditions
      call gungho_diagnostics_driver( mesh,       &
                                      twod_mesh,  &
                                      model_data, &
                                      clock,      &
                                      nodal_output_on_w3 )

      if ( limited_area )then
        if ( lbc_method == lbc_method_onion_layer ) call write_masks()
      end if
    end if

  end subroutine write_initial_output

end module initial_output_mod
