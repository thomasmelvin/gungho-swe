!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a vertical semi-Lagragian advection of rho-type field
!!        using the vertical wind (w) only.
!> @Details The 1D vertical advective transport equation for a W3 variable
!!          is solved using a semi-Lagragian advection scheme.

module vertical_sl_rho_kernel_mod

use argument_mod,          only : arg_type,              &
                                  GH_FIELD, GH_SCALAR,   &
                                  GH_REAL, GH_INTEGER,   &
                                  GH_READWRITE, GH_READ, &
                                  CELL_COLUMN
use fs_continuity_mod,     only : W2, W3
use constants_mod,         only : r_def, i_def, tiny_eps
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
type, public, extends(kernel_type) :: vertical_sl_rho_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                     &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2), & ! departure points
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3)  & ! rho
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_sl_rho_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vertical_sl_rho_code

contains

!-------------------------------------------------------------------------------
!> @details This kernel calculates the departure point of w/theta-points using
!!          only w (i.e., vertical motion only), then interpolate theta at the
!!          departure point using 1d-Cubic-Lagrange interpolation.
!> @param[in]     nlayers      The number of layers
!> @param[in]     dep_pts_z    The vertical departure distance used for SL advection
!> @param[in,out] rho          The rho field at time level n --> rho after SL-advection
!> @param[in]     ndf_w2       The number of degrees of freedom per cell
!!                             on W2 space
!> @param[in]     undf_w2      The number of unique degrees of freedom
!!                             on W2 space
!> @param[in]     map_w2       The dofmap for the cell at the base of the column
!!                             on W2 space
!> @param[in]     ndf_w3       The number of degrees of freedom per cell
!!                             on w3 space
!> @param[in]     undf_w3      The number of unique degrees of freedom
!!                             on w3 space
!> @param[in]     map_w3       The dofmap for the cell at the base of the column
!!                             on w3 space
!-------------------------------------------------------------------------------

subroutine vertical_sl_rho_code( nlayers,                 &
                                 dep_pts_z,               &
                                 rho,                     &
                                 ndf_w2, undf_w2, map_w2, &
                                 ndf_w3, undf_w3, map_w3 )


  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                         :: nlayers
  integer(kind=i_def), intent(in)                         :: ndf_w2
  integer(kind=i_def), intent(in)                         :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)      :: map_w2
  integer(kind=i_def), intent(in)                         :: ndf_w3
  integer(kind=i_def), intent(in)                         :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)      :: map_w3

  real(kind=r_def), dimension(undf_w2), intent(in)        :: dep_pts_z
  real(kind=r_def), dimension(undf_w3), intent(inout)     :: rho

  real(kind=r_def), allocatable :: dist(:)
  real(kind=r_def), allocatable :: dist_w2(:)
  real(kind=r_def), allocatable :: zm(:)
  real(kind=r_def), allocatable :: zd(:)
  real(kind=r_def), allocatable :: f0(:)
  real(kind=r_def), allocatable :: fd(:)

  integer(kind=i_def) :: nz,nzl
  integer(kind=i_def) :: k,k0,k1,k2,k3,k4,k5
  real(kind=r_def)    :: c0,c1,c2,c3,c4,c5
  real(kind=r_def)    :: sm3,sm2,sm1,ss,sp1,sp2

  nz = nlayers
  nzl = nz + 1
  allocate( dist_w2(1:nzl) )
  allocate( dist(1:nz) )
  allocate( zm(1:nz) )
  allocate( zd(1:nz) )
  allocate( f0(1:nz) )
  allocate( fd(1:nz) )

  ! Extract and fill local column from global data
  ! Map departure points into 1d-array dist_w2
  do k=0,nlayers
    dist_w2(k+1) = dep_pts_z(map_w2(5)+k)
  end do
  ! Map global field into 1d-array fo
  do k=0,nlayers - 1
    f0(k+1) = rho(map_w3(1)+k)
  end do

  ! Create a local 1-d SL problem
  ! dist = departure distance at centre (zm) of cells
  do k=1,nz
    dist(k) = 0.5_r_def * ( dist_w2(k) + dist_w2(k+1) )
    zm(k) = real(k-1,r_def)
  end do
  ! zd = departure of zm
  do k=1,nz
     zd(k) = zm(k) - dist(k)
     zd(k) = min(zm(nz),max(zm(1),zd(k)))
  end do

  if ( vertical_sl_order == vertical_sl_order_cubic) then

    do k=1,nz
      k2 = floor(zd(k)) + 1
      ss = zd(k) - zm(k2)
      k1 = max(1 , k2 - 1 )
      k3 = min(nz, k2 + 1 )
      k4 = min(nz, k2 + 2 )

      sm1 = ss - 1.0_r_def
      sm2 = ss - 2.0_r_def
      sp1 = ss + 1.0_r_def

      c1 = -(1.0_r_def/6.0_r_def) * ss  * sm1 * sm2
      c2 =  (1.0_r_def/2.0_r_def) * sp1 * sm1 * sm2
      c3 = -(1.0_r_def/2.0_r_def) * ss  * sp1 * sm2
      c4 =  (1.0_r_def/6.0_r_def) * ss  * sp1 * sm1

      ! Do linear intepolation if you are next to the boundary.
      ! This could be removed but this is equivalent to imposing
      ! zero-gradient assumption near the top and bottom boundaries

      if( k1 == k2 .or. k3 == k4) then
         c1 = 0.0_r_def
         c4 = 0.0_r_def
         c2 = 1.0_r_def - ss
         c3 = ss
      end if
      fd(k) = c1*f0(k1) + c2*f0(k2) + c3*f0(k3) + c4*f0(k4)
    end do

  else if ( vertical_sl_order == vertical_sl_order_quintic) then

    do k=1,nz
      k2 = floor(zd(k)) + 1
      ss = zd(k) - zm(k2)
      k1 = max(1 , k2 - 1 )
      k0 = max(1 , k2 - 2 )
      k3 = min(nz, k2 + 1 )
      k4 = min(nz, k2 + 2 )
      k5 = min(nz, k2 + 3 )

      sm1 = ss - 1.0_r_def
      sm2 = ss - 2.0_r_def
      sm3 = ss - 3.0_r_def
      sp1 = ss + 1.0_r_def
      sp2 = ss + 2.0_r_def

      ! If the stencil extend beyond the boundaries
      !    reduces the order of interpolation

      if( k0 == k1 .or. k4 == k5 ) then
        if( k1 == k2 .or. k3 == k4 ) then
        ! Revert to linear interpolation
          c0 = 0.0_r_def
          c1 = 0.0_r_def
          c2 = 1.0_r_def - ss
          c3 = ss
          c4 = 0.0_r_def
          c5 = 0.0_r_def
        else
        ! Revert to cubic interpolation
          c0 = 0.0_r_def
          c1 = -(1.0_r_def/6.0_r_def) * ss  * sm1 * sm2
          c2 =  (1.0_r_def/2.0_r_def) * sp1 * sm1 * sm2
          c3 = -(1.0_r_def/2.0_r_def) * ss  * sp1 * sm2
          c4 =  (1.0_r_def/6.0_r_def) * ss  * sp1 * sm1
          c5 = 0.0_r_def
        end if

      else
       ! Perform quintic interpolation
        c0 = -(1.0_r_def/120.0_r_def) * sp1 * ss  * sm1 * sm2 * sm3
        c1 =  (1.0_r_def/24.0_r_def ) * sp2 * ss  * sm1 * sm2 * sm3
        c2 = -(1.0_r_def/12.0_r_def ) * sp2 * sp1 * sm1 * sm2 * sm3
        c3 =  (1.0_r_def/12.0_r_def ) * sp2 * sp1 * ss  * sm2 * sm3
        c4 = -(1.0_r_def/24.0_r_def ) * sp2 * sp1 * ss  * sm1 * sm3
        c5 =  (1.0_r_def/120.0_r_def) * sp2 * sp1 * ss  * sm1 * sm2
      end if

      fd(k) = c0*f0(k0) + c1*f0(k1) + c2*f0(k2) + &
              c3*f0(k3) + c4*f0(k4) + c5*f0(k5)
    end do

  end if

  ! Map the local 1d-array back to the global field

  do k=0,nlayers - 1
    rho(map_w3(1)+k) = fd(k+1)
  end do

  deallocate( zm, zd, dist, dist_w2, f0, fd )

end subroutine vertical_sl_rho_code

end module vertical_sl_rho_kernel_mod
