!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief The kernel computes edge lengths at theta points.
!> @details Kernel to compute the edge lengths at theta points by finding the
!> minimum edge length of the four w2 edges surrounding each theta point
module calc_delta_at_wtheta_kernel_mod

use argument_mod,            only : arg_type,          &
                                    GH_FIELD, GH_REAL, &
                                    GH_READ, GH_WRITE, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def, LARGE_REAL_POSITIVE
use fs_continuity_mod,       only : W2, Wtheta
use kernel_mod,              only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_delta_at_wtheta_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                 &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W2)      &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: calc_delta_at_wtheta_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: calc_delta_at_wtheta_code
contains

!> @brief compute edge lengths based on cell volume and face area.
!! @param[in] nlayers Number of layers
!! @param[in,out] delta_at_wtheta Minimum edge length at theta points
!! @param[in] dx_at_w2 edge length on w2 points
!! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
!! @param[in] undf_wtheta  Number of unique degrees of freedom  for wtheta
!! @param[in] map_wtheta Dofmap for the cell at the base of the column for wtheta
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number of unique degrees of freedom for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
subroutine calc_delta_at_wtheta_code(nlayers,                 &
                       delta_at_wtheta, dx_at_w2,             &
                       ndf_wtheta, undf_wtheta, map_wtheta,   &
                       ndf_w2, undf_w2, map_w2)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_wtheta, undf_w2, undf_wtheta

  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(undf_wtheta), intent(inout) :: delta_at_wtheta
  real(kind=r_def), dimension(undf_w2), intent(in)        :: dx_at_w2

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def)    :: mindx

  delta_at_wtheta(map_wtheta(1)) = LARGE_REAL_POSITIVE

  do k = 0, nlayers-1
    mindx = LARGE_REAL_POSITIVE
    do df = 1, 4
      mindx = min(mindx,dx_at_w2(map_w2(df) + k))
    end do
    delta_at_wtheta(map_wtheta(2) + k) = mindx
    delta_at_wtheta(map_wtheta(1) + k) = min(mindx, delta_at_wtheta(map_wtheta(1) + k))

  end do

end subroutine calc_delta_at_wtheta_code

end module calc_delta_at_wtheta_kernel_mod
