!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the average of nearest W3 space to W1 space.
!> @details Kernel to average a W3 lower-level field to the horizontal parts
!!          of a W1 field. The method is valid for the top-most DoFs, horizontal
!!          and lower level of the lowest-order finite elements on a cubed-sphere mesh.

module w3_to_w1_average_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_INC, GH_READ,   &
                                CELL_COLUMN
  use constants_mod,     only : i_def, r_def
  use fs_continuity_mod, only : W1, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: w3_to_w1_average_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/            &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3)  &
         /)
    integer :: iterates_over = CELL_COLUMN
  contains
    procedure, nopass :: w3_to_w1_average_code
  end type

  ! -----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: w3_to_w1_average_code

  contains

  !> @brief Computes a 1-2-1 Filter from W3 to W1 space.
  !> @param[in]     nlayers   Number of layers
  !> @param[in,out] field_w1  Output field from filter on W0 space
  !> @param[in]     field_w3  Input field for filter on W3 space
  !> @param[in]     ndf_w1    Number of degrees of freedom per cell for W1
  !> @param[in]     undf_w1   Number of unique degrees of freedom for W1
  !> @param[in]     map_w1    Dofmap for the cell at the base of the column for W1
  !> @param[in]     ndf_w3    Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3   Number of unique degrees of freedom for W3
  !> @param[in]     map_w3  ss  Dofmap for the cell at the base of the column for W3
  subroutine w3_to_w1_average_code(nlayers,                 &
                                   field_w1, field_w3,      &
                                   ndf_w1, undf_w1, map_w1, &
                                   ndf_w3, undf_w3, map_w3)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w1, undf_w1
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3

    integer(kind=i_def), intent(in), dimension(ndf_w1) :: map_w1
    integer(kind=i_def), intent(in), dimension(ndf_w3) :: map_w3

    real(kind=r_def), intent(inout), dimension(undf_w1) :: field_w1
    real(kind=r_def), intent(in),    dimension(undf_w3) :: field_w3

    ! Internal variables
    integer(kind=i_def) :: df, k

    do k = 0, nlayers-1
      do df = 5,8 ! Loop at the top
        field_w1(map_w1(df) + k) = field_w1(map_w1(df) + k) + field_w3(map_w3(1) + k)/4.0_r_def
      end do
    end do

  end subroutine w3_to_w1_average_code

end module w3_to_w1_average_kernel_mod
