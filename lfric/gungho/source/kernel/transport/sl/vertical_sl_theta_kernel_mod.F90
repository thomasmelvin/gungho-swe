!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a semi-Lagragian advection of theta-type tracer using
!!        the vertical wind (w) only.
!> @Details The 1D vertical advective transport equation for a Wtheta variable
!!          is solved using a semi-Lagragian advection scheme.

module vertical_sl_theta_kernel_mod

use argument_mod,          only : arg_type,              &
                                  GH_FIELD, GH_SCALAR,   &
                                  GH_READWRITE, GH_READ, &
                                  GH_REAL, CELL_COLUMN
use fs_continuity_mod,     only : W2, Wtheta
use constants_mod,         only : r_def, i_def
use kernel_mod,            only : kernel_type
! TODO #3011: these config options should be passed through as arguments
use transport_config_mod,  only : vertical_sl_order,       &
                                  vertical_sl_order_cubic, &
                                  vertical_sl_order_quintic

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed
!>                                      by the PSy layer.
type, public, extends(kernel_type) :: vertical_sl_theta_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                      &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,      W2    ), & ! departure points
       arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, Wtheta)  & ! theta
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_sl_theta_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vertical_sl_theta_code

contains

!-------------------------------------------------------------------------------
!> @details This kernel calculates the departure point of w/theta-points using
!!          only w (i.e., vertical motion only), then interpolates theta at the
!!          departure point using 1d-Cubic-Lagrange interpolation.
!> @param[in]     nlayers      The number of layers
!> @param[in]     dep_pts_z    The vertical departure distance used for SL advection
!> @param[in,out] theta        The theta field at time level n --> theta_d (SL)
!> @param[in]     ndf_w2       The number of degrees of freedom per cell
!!                             on W2 space
!> @param[in]     undf_w2      The number of unique degrees of freedom
!!                             on W2 space
!> @param[in]     map_w2       The dofmap for the cell at the base of the column
!!                             on W2 space
!> @param[in]     ndf_wtheta   The number of degrees of freedom per cell
!!                             on Wtheta space
!> @param[in]     undf_wtheta  The number of unique degrees of freedom
!!                             on Wtheta space
!> @param[in]     map_wtheta   The dofmap for the cell at the base of the column
!!                             on Wtheta space
!-------------------------------------------------------------------------------

subroutine vertical_sl_theta_code( nlayers,                             &
                                   dep_pts_z,                           &
                                   theta,                               &
                                   ndf_w2, undf_w2, map_w2,             &
                                   ndf_wtheta, undf_wtheta, map_wtheta  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                         :: nlayers
  integer(kind=i_def), intent(in)                         :: ndf_w2
  integer(kind=i_def), intent(in)                         :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)      :: map_w2
  integer(kind=i_def), intent(in)                         :: ndf_wtheta
  integer(kind=i_def), intent(in)                         :: undf_wtheta
  integer(kind=i_def), dimension(ndf_wtheta), intent(in)  :: map_wtheta

  real(kind=r_def), dimension(undf_w2), intent(in)        :: dep_pts_z
  real(kind=r_def), dimension(undf_wtheta), intent(inout) :: theta

  integer(kind=i_def) :: ncelledges
  real(kind=r_def)    :: zd, z_surf, z_top
  real(kind=r_def),allocatable :: dist(:)
  real(kind=r_def),allocatable :: theta_local(:)
  real(kind=r_def),allocatable :: theta_d_local(:)
  integer(kind=i_def) :: k,k0,k1,k2,k3,k4,k5
  real(kind=r_def)    :: c0,c1,c2,c3,c4,c5
  real(kind=r_def)    :: sm3,sm2,sm1,ss,sp1,sp2

  ncelledges = nlayers + 1_i_def
  allocate(dist(1:ncelledges))
  allocate(theta_local(1:ncelledges))
  allocate(theta_d_local(1:ncelledges))

  ! Extract and fill local column from global variables
  z_surf = 1.0_r_def
  z_top  = real(ncelledges,r_def)

  do k=0,nlayers
       dist(k+1)        = dep_pts_z(map_w2(5)+k)
       theta_local(k+1) = theta(map_wtheta(1)+k)
  end do

  if ( vertical_sl_order == vertical_sl_order_cubic ) then

    do k=1,ncelledges
      zd = real(k,r_def) - dist(k)
      zd = min(z_top,max(z_surf,zd))
      k2 = floor(zd)
      ss = zd - real(k2, r_def)
      k1 = max(1_i_def, k2 - 1_i_def)
      k3 = min(ncelledges, k2 + 1_i_def)
      k4 = min(ncelledges, k2 + 2_i_def)

      sm1 = ss - 1.0_r_def
      sm2 = ss - 2.0_r_def
      sp1 = ss + 1.0_r_def

      c1 = -(1.0_r_def/6.0_r_def) * ss  * sm1 * sm2
      c2 =  (1.0_r_def/2.0_r_def) * sp1 * sm1 * sm2
      c3 = -(1.0_r_def/2.0_r_def) * ss  * sp1 * sm2
      c4 =  (1.0_r_def/6.0_r_def) * ss  * sp1 * sm1

      ! Do linear intepolation if you are next to the boundary.
      ! This if could be removed but this is equivalent to imposing
      ! zero-gradient assumption near the top and bottom boundaries

      if( k1 == k2 .or. k3 == k4) then
         c1 = 0.0_r_def
         c4 = 0.0_r_def
         c2 = 1.0_r_def - ss
         c3 = ss
      end if

      theta_d_local(k) = c1*theta_local(k1) + c2*theta_local(k2) + &
                         c3*theta_local(k3) + c4*theta_local(k4)
    end do

  else if ( vertical_sl_order == vertical_sl_order_quintic) then

    do k=1,ncelledges
      zd = real(k,r_def) - dist(k)
      zd = min(z_top,max(z_surf,zd))
      k2 = floor(zd)
      ss = zd - real(k2, r_def)
      k1 = max(1_i_def, k2 - 1_i_def)
      k0 = max(1_i_def, k2 - 2_i_def)
      k3 = min(ncelledges, k2 + 1_i_def)
      k4 = min(ncelledges, k2 + 2_i_def)
      k5 = min(ncelledges, k2 + 3_i_def)

      sm1 = ss - 1.0_r_def
      sm2 = ss - 2.0_r_def
      sm3 = ss - 3.0_r_def
      sp1 = ss + 1.0_r_def
      sp2 = ss + 2.0_r_def

      c0 = -(1.0_r_def/120.0_r_def) * sp1 * ss  * sm1 * sm2 * sm3
      c1 =  (1.0_r_def/24.0_r_def ) * sp2 * ss  * sm1 * sm2 * sm3
      c2 = -(1.0_r_def/12.0_r_def ) * sp2 * sp1 * sm1 * sm2 * sm3
      c3 =  (1.0_r_def/12.0_r_def ) * sp2 * sp1 * ss  * sm2 * sm3
      c4 = -(1.0_r_def/24.0_r_def ) * sp2 * sp1 * ss  * sm1 * sm3
      c5 =  (1.0_r_def/120.0_r_def) * sp2 * sp1 * ss  * sm1 * sm2

      if( k5 == k4 .or. k1 == k0) then
         c0 = 0.0_r_def
         c5 = 0.0_r_def
         c1 = -(1.0_r_def/6.0_r_def) * ss  * sm1 * sm2
         c2 =  (1.0_r_def/2.0_r_def) * sp1 * sm1 * sm2
         c3 = -(1.0_r_def/2.0_r_def) * ss  * sp1 * sm2
         c4 =  (1.0_r_def/6.0_r_def) * ss  * sp1 * sm1
      end if
      if( k1 == k2 .or. k3 == k4) then
         c0 = 0.0_r_def
         c1 = 0.0_r_def
         c4 = 0.0_r_def
         c5 = 0.0_r_def
         c2 = 1.0_r_def - ss
         c3 = ss
      end if

      theta_d_local(k) = c0*theta_local(k0) + c1*theta_local(k1) + c2*theta_local(k2) + &
                         c3*theta_local(k3) + c4*theta_local(k4) + c5*theta_local(k5)
    end do

  end if

  do k=0,nlayers
    theta(map_wtheta(1)+k) = theta_d_local(k+1)
  end do

  deallocate(dist      )
  deallocate(theta_local  )
  deallocate(theta_d_local)

end subroutine vertical_sl_theta_code

end module vertical_sl_theta_kernel_mod
