!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Kernel to get the vertical grid spacing of the W3 cells.

module get_dz_w3_kernel_mod

  use argument_mod,          only : arg_type,          &
                                    GH_FIELD, GH_REAL, &
                                    GH_READ, GH_WRITE, &
                                    CELL_COLUMN
  use fs_continuity_mod,     only : W3, W2
  use constants_mod,         only : r_def, i_def
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : T, B

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: get_dz_w3_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/             &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: get_dz_w3_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: get_dz_w3_code

  contains

  !> @details Calculate the vertical grid spacing of the W3 cells (dz) using the
  !!          physical height at W2 DoFs. This kernel is only designed for
  !!          lowest order finite elements.
  !> @param[in]     nlayers Number of layers
  !> @param[in,out] dz      Vertical grid spacing for a W3 cell
  !> @param[in]     height  The height at W2 points
  !> @param[in]     ndf_w3  Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3 Number of unique degrees of freedom for W3
  !> @param[in]     map_w3  The dofmap for the cell at the base of the column for W3
  !> @param[in]     ndf_w2  Number of degrees of freedom for W2 per cell
  !> @param[in]     undf_w2 Number of unique degrees of freedom for W2
  !> @param[in]     map_w2  The dofmap for the cell at the base of the column for W2
  subroutine get_dz_w3_code( nlayers, &
                             dz,      &
                             height,  &
                             ndf_w3,  &
                             undf_w3, &
                             map_w3,  &
                             ndf_w2,  &
                             undf_w2, &
                             map_w2 )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)   :: nlayers
    integer(kind=i_def), intent(in)   :: undf_w3
    integer(kind=i_def), intent(in)   :: ndf_w3
    integer(kind=i_def), intent(in)   :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in)   :: undf_w2
    integer(kind=i_def), intent(in)   :: ndf_w2
    integer(kind=i_def), intent(in)   :: map_w2(ndf_w2)
    real(kind=r_def), intent(inout)   :: dz(undf_w3)
    real(kind=r_def), intent(in)      :: height(undf_w2)

    ! Vertical cell index
    integer(kind=i_def) :: k

    ! dz is the difference of height at the top and bottom of a cell
    do k = 0, nlayers-1
      dz(map_w3(1) + k) = height(map_w2(T)+k) - height(map_w2(B)+k)
    end do

  end subroutine get_dz_w3_code

end module get_dz_w3_kernel_mod
