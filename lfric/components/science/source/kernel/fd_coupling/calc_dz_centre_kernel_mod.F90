!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the difference between cell face layer heights,
!!        with the result evaluated at cell centre layer heights.
module calc_dz_centre_kernel_mod

use argument_mod,  only: arg_type,                              &
                         GH_FIELD, GH_SCALAR, GH_READ, GH_INC,  &
                         CELL_COLUMN, ANY_SPACE_1, ANY_SPACE_2, &
                         GH_INTEGER, GH_REAL
use constants_mod, only: r_def, i_def
use kernel_mod,    only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_dz_centre_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                         &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,  ANY_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_SPACE_2), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)               &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: calc_dz_centre_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: calc_dz_centre_code

contains

!> @details For lowest-order fields only, calculates the finite-difference
!!          vertical distance between adjacent cell top/bottom faces. As such,
!!          the result is evaluated at cell centre in the vertical. The code
!!          works for any horizontally located position within the cell
!!          (centre or horizontal faces), by use of the input variable
!!          "df_to_do" which describes how many horizontal locations to
!!          calculate dz for.
!> @param[in]     nlayers     Number of layers
!> @param[in,out] dz          Distance between faces (evaluated at cell centre)
!> @param[in]     height      Height field of cell faces
!> @param[in]     df_to_do    Number of dofs to calculate dz for
!> @param[in]     ndf_dz      Number of degrees of freedom per cell for dz
!> @param[in]     undf_dz     Number of unique degrees of freedom for dz
!> @param[in]     map_dz      Dofmap for the cell at the base of the column for dz
!> @param[in]     ndf_height  Number of degrees of freedom per cell for height
!> @param[in]     undf_height Number of unique degrees of freedom for height
!> @param[in]     map_height  Dofmap for the cell at the base of the column for height
subroutine calc_dz_centre_code(nlayers, dz, height, df_to_do, &
                               ndf_dz, undf_dz, map_dz,       &
                               ndf_height, undf_height, map_height)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, df_to_do
  integer(kind=i_def), intent(in) :: ndf_height, ndf_dz, undf_height, undf_dz

  integer(kind=i_def), dimension(ndf_dz), intent(in)     :: map_dz
  integer(kind=i_def), dimension(ndf_height), intent(in) :: map_height

  real(kind=r_def), dimension(undf_dz), intent(inout)  :: dz
  real(kind=r_def), dimension(undf_height), intent(in) :: height

  ! Internal variables
  integer(kind=i_def) :: df, k

  do df = 1, df_to_do
    do k = 0, nlayers-1
      dz(map_dz(df)+k) = height(map_height(df)+k+1) &
                       - height(map_height(df)+k)
    end do
  end do

end subroutine calc_dz_centre_code

end module calc_dz_centre_kernel_mod
