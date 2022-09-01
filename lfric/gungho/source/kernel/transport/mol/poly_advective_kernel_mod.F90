!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes horizontal advective update through fitting a high order
!!        upwind reconstruction.
!> @details Compute the advective update (u.grad) for a tracer field using a high order
!!          polynomial fit to the integrated tracer values. The stencil used for the
!!          polynomial is centred on the upwind cell.
!!          This method is only valid for lowest order elements.

module poly_advective_kernel_mod

use argument_mod,      only : arg_type,          &
                              GH_FIELD, GH_REAL, &
                              GH_WRITE, GH_READ, &
                              CELL_COLUMN
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W1, W2, Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: poly_advective_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                 &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),     &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W1)      &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: poly_advective_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly_advective_code
contains

!> @brief Computes the horizontal advective update for a tracer
!> @param[in]     nlayers   Number of layers
!> @param[in,out] advective Advective update field to compute
!> @param[in]     wind      Wind field
!> @param[in]     tracer    Pointwise tracer field to advect stored on edge centres
!> @param[in]     ndf_wt    Number of degrees of freedom per cell
!> @param[in]     undf_wt   Number of unique degrees of freedom for the
!!                          advective_update field
!> @param[in]     map_wt    Dofmap for the cell at the base of the column
!> @param[in]     ndf_w2    Number of degrees of freedom per cell for the wind fields
!> @param[in]     undf_w2   Number of unique degrees of freedom for the wind fields
!> @param[in]     map_w2    Dofmap for the cell at the base of the column for the wind fields
!> @param[in]     ndf_w1    Number of degrees of freedom per cell for the
!!                          reconstructed tracer
!> @param[in]     undf_w1   Number of unique degrees of freedom for the
!!                          reconstructed tracer
!> @param[in]     map_w1    Dofmap for the cell at the base of the column for the
!!                          reconstructed tracer
subroutine poly_advective_code( nlayers,              &
                                advective,            &
                                wind,                 &
                                tracer,               &
                                ndf_wt,               &
                                undf_wt,              &
                                map_wt,               &
                                ndf_w2,               &
                                undf_w2,              &
                                map_w2,               &
                                ndf_w1,               &
                                undf_w1,              &
                                map_w1                &
                              )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt, ndf_w2, ndf_w1
  integer(kind=i_def), intent(in)                    :: undf_wt, undf_w2, undf_w1
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w1), intent(in) :: map_w1

  real(kind=r_def), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_w1), intent(in)    :: tracer

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def)    :: u, v, dtdx, dtdy

  ! Horizontal advective update computation
  ! Bottom point
  u =  0.25_r_def*( wind(map_w2(1)) + wind(map_w2(3)) )
  v = -0.25_r_def*( wind(map_w2(2)) + wind(map_w2(4)) )
  dtdx = tracer(map_w1(3)) - tracer(map_w1(1))
  dtdy = tracer(map_w1(4)) - tracer(map_w1(2))
  advective(map_wt(1)) = u*dtdx + v*dtdy

  do k = 1, nlayers - 1
    u =  0.25_r_def*( wind(map_w2(1) + k - 1 ) + wind(map_w2(1) + k ) &
                    + wind(map_w2(3) + k - 1 ) + wind(map_w2(3) + k ) )
    v = -0.25_r_def*( wind(map_w2(2) + k - 1 ) + wind(map_w2(2) + k ) &
                    + wind(map_w2(4) + k - 1 ) + wind(map_w2(4) + k ) )
    dtdx = tracer(map_w1(3) + k) - tracer(map_w1(1) + k)
    dtdy = tracer(map_w1(4) + k) - tracer(map_w1(2) + k)
    advective(map_wt(1)+k) = u*dtdx + v*dtdy
  end do
  ! Final point on the top
  k = nlayers-1
  u =  0.25_r_def*( wind(map_w2(1) + k) + wind(map_w2(3) + k) )
  v = -0.25_r_def*( wind(map_w2(2) + k) + wind(map_w2(4) + k) )
  dtdx = tracer(map_w1(11) + k) - tracer(map_w1(9) + k)
  dtdy = tracer(map_w1(12) + k) - tracer(map_w1(10) + k)
  advective(map_wt(2)+k) = u*dtdx + v*dtdy

end subroutine poly_advective_code

end module poly_advective_kernel_mod
