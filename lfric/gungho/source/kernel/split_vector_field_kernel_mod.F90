!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Split a the wind field (W2) into horizontal and vertical components.
module split_vector_field_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type, &
                                    GH_FIELD, GH_REAL,   &
                                    GH_READ, GH_INC,     &
                                    GH_BASIS,            &
                                    CELL_COLUMN, GH_EVALUATOR
use constants_mod,           only : r_def, i_def, EPS
use fs_continuity_mod,       only : W2

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: split_vector_field_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/            &
       arg_type(GH_FIELD, GH_REAL, GH_INC,  W2), &
       arg_type(GH_FIELD, GH_REAL, GH_INC,  W2), &
       arg_type(GH_FIELD, GH_REAL, GH_READ, W2)  &
       /)
  type(func_type) :: meta_funcs(1) = (/          &
       func_type(W2, GH_BASIS)                   &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: split_vector_field_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: split_vector_field_code
contains

!> @details Split a wind field (uvw) into horizontal (uv) and
!>          vertical (w) components.
!>          It uses a vertical basis vector (0,0,1) to identify which
!>          components are vertical or horizontal using the dot-product.
!> @param[in]     nlayers Number of layers
!> @param[in,out] uv Horizontal wind
!> @param[in,out] w Vertical wind
!> @param[in]     uvw 3D wind
!> @param[in]     ndf Number of degrees of freedom per cell
!> @param[in]     undf Number of unique degrees of freedom
!> @param[in]     map Dofmap for the cell at the base of the column
!> @param[in]     basis Basis functions evaluated at nodal points

subroutine split_vector_field_code(nlayers,              &
                                   uv, w, uvw,           &
                                   ndf, undf, map, basis )

  implicit none

  ! Arguments
  integer(kind=i_def),                 intent(in) :: nlayers
  integer(kind=i_def),                 intent(in) :: ndf, undf
  integer(kind=i_def), dimension(ndf), intent(in) :: map

  real(kind=r_def), dimension(3,ndf,ndf), intent(in)    :: basis
  real(kind=r_def), dimension(undf),      intent(inout) :: uv, w
  real(kind=r_def), dimension(undf),      intent(in)    :: uvw

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def), dimension(3) :: z_hat

  z_hat = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)
  do k = 0, nlayers-1
    do df = 1, ndf
      if ( abs(dot_product(basis(:,df,df),z_hat)) > EPS ) then
        w(map(df)+k) = uvw(map(df)+k)
        uv(map(df)+k) = 0.0_r_def
      else
        uv(map(df)+k) = uvw(map(df)+k)
        w(map(df)+k)  = 0.0_r_def
      end if
    end do
  end do

end subroutine split_vector_field_code

end module split_vector_field_kernel_mod
