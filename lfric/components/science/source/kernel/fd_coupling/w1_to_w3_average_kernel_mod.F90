!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the average of nearest W1 space to W3 space.
!> @details Kernel to average a W1 lower-level field to the horizontal parts
!!          of a W3 field. The method is valid for the top-most DoFs, horizontal
!!          and lower level of the lowest-order finite elements on a cubed-sphere mesh.

module w1_to_w3_average_kernel_mod

  use argument_mod,      only : arg_type,               &
                                GH_FIELD, GH_REAL,      &
                                GH_READWRITE, GH_READ,  &
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
  type, public, extends(kernel_type) :: w1_to_w3_average_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W1)  &
         /)
    integer :: iterates_over = CELL_COLUMN
  contains
    procedure, nopass :: w1_to_w3_average_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: w1_to_w3_average_code

  contains

  !> @brief Computes a 1-2-1 Filter from W1 to W3 space.
  !> @param[in]     nlayers       Number of layers
  !> @param[in,out] field_w3      Output field from filter on W3 space
  !> @param[in]     field_w1      Input field for filter on W1 space
  !> @param[in]     rmultiplicity Reciprocal of how many times the dof has been visited in total
  !> @param[in]     ndf_w3        Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3       Number of unique degrees of freedom for W3
  !> @param[in]     map_w3        Dofmap for the cell at the base of the column for W3
  !> @param[in]     ndf_w1        Number of degrees of freedom per cell for W1
  !> @param[in]     undf_w1       Number of unique degrees of freedom for W1
  !> @param[in]     map_w1        Dofmap for the cell at the base of the column for W1
  subroutine w1_to_w3_average_code(nlayers,                              &
                                   field_w3, field_w1, rmultiplicity_w1, &
                                   ndf_w3, undf_w3, map_w3,              &
                                   ndf_w1, undf_w1, map_w1)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_w1, undf_w1

    integer(kind=i_def), intent(in), dimension(ndf_w3) :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_w1) :: map_w1

    real(kind=r_def), intent(inout), dimension(undf_w3) :: field_w3
    real(kind=r_def), intent(in),    dimension(undf_w1) :: field_w1, rmultiplicity_w1

    ! Internal variables
    integer(kind=i_def) :: df, k

    do k = 0, nlayers-1
      field_w3(map_w3(1) + k) = 0.0_r_def
      do df = 5,8 ! Loop at the top
        field_w3(map_w3(1) + k) = field_w3(map_w3(1) + k) +                 &
                    field_w1(map_w1(df) + k)*rmultiplicity_w1(map_w1(df) + k)
      end do
    end do

  end subroutine w1_to_w3_average_code

end module w1_to_w3_average_kernel_mod
