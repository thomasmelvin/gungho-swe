!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Divides the horizontal Piola wind values by Det(J) at W2 dofs.
!> @details Calculates the departure wind in the horizontal by dividing the
!!          Piola wind values by Det(J) at horizontal W2 dofs. The departure
!!          wind is used by the flux form semi-Lagrangian (FFSL) scheme.
module ffsl_hori_dep_wind_kernel_mod

  use argument_mod,       only : arg_type, func_type,       &
                                 GH_FIELD, GH_READ, GH_INC, &
                                 GH_REAL, ANY_SPACE_9,      &
                                 GH_DIFF_BASIS, GH_BASIS,   &
                                 CELL_COLUMN, GH_EVALUATOR
  use constants_mod,      only : r_def, i_def
  use flux_direction_mod, only : x_direction, y_direction
  use fs_continuity_mod,  only : W2
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: ffsl_hori_dep_wind_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                      &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),         &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),         &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9) &
         /)
    type(func_type) :: meta_funcs(2) = (/                    &
         func_type(W2,          GH_BASIS),                   &
         func_type(ANY_SPACE_9, GH_DIFF_BASIS)               &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: ffsl_hori_dep_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: ffsl_hori_dep_wind_code

contains

  !> @brief Divides the horizontal Piola wind values by Det(J) at W2 dofs.
  !> @param[in]     nlayers          Number of layers
  !> @param[in,out] u_departure_wind Output field containing the departure
  !!                                 wind used to calculate departure points
  !> @param[in]     u_piola          Field for the Piola wind
  !> @param[in]     detj_at_w2       Input field containing the detj values at W2 locations
  !> @param[in]     ndf              Number of degrees of freedom per cell for the output field
  !> @param[in]     undf             Number of unique degrees of freedom for the output field
  !> @param[in]     map              Dofmap for the cell at the base of the column for the output field
  !> @param[in]     direction        The direction in which the winds are calculated
  subroutine ffsl_hori_dep_wind_code(nlayers,          &
                                     u_departure_wind, &
                                     u_piola,          &
                                     detj_at_w2,       &
                                     ndf, undf, map,   &
                                     direction         &
                                     )

    implicit none

    ! Arguments
    integer(kind=i_def),                        intent(in)    :: nlayers
    integer(kind=i_def),                        intent(in)    :: ndf, undf
    integer(kind=i_def), dimension(ndf),        intent(in)    :: map
    real(kind=r_def), dimension(undf),          intent(in)    :: u_piola
    real(kind=r_def), dimension(undf),          intent(in)    :: detj_at_w2
    real(kind=r_def), dimension(undf),          intent(inout) :: u_departure_wind
    integer(kind=i_def),                        intent(in)    :: direction

    ! Internal variables
    integer(kind=i_def)                  :: df, k
    real(kind=r_def)                     :: mult_factor

    ! Change the sign of the output winds depending on whether it is an x or y
    ! direction update
    if (direction == x_direction ) then
      mult_factor = 1.0_r_def
    else if (direction == y_direction ) then
      mult_factor = -1.0_r_def
    end if

    ! Divide the input Piola wind values by the corresponding detJ value.
    do k = 0, nlayers-1
      do df = 1,4
        u_departure_wind(map(df)+k) = mult_factor*u_piola(map(df)+k)/detj_at_w2(map(df)+k)
      end do
      u_departure_wind(map(5)+k) = -9999.0
      u_departure_wind(map(6)+k) = -9999.0
    end do

  end subroutine ffsl_hori_dep_wind_code

end module ffsl_hori_dep_wind_kernel_mod
