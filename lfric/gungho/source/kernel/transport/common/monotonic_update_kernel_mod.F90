!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which enforces monotonicity on an advective update.
!> @details Modifies the advective update A(theta) such that
!!          \f[ M_\theta (theta^{n+1} - theta^{n}) + \Delta t A(\theta) \f]
!!          Returns a monotonic theta^{n+1}. The monotonticity is enforced by
!!          ensuring that the the implied theta^{n+1} lies within the range
!!          [min(theta_s),max(theta_s)] for all s in the stencil of values used
!!          to compute theta.
module monotonic_update_kernel_mod

use argument_mod,      only : arg_type,              &
                              GH_FIELD, GH_SCALAR,   &
                              GH_REAL, GH_INTEGER,   &
                              GH_READWRITE, GH_READ, &
                              STENCIL, CROSS, CELL_COLUMN
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: monotonic_update_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                         &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                 &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta, STENCIL(CROSS)), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                 &
       arg_type(GH_SCALAR, GH_REAL,    GH_READ),                              &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                               &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: monotonic_update_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: monotonic_update_code

contains

!> @brief Computes the horizontal fluxes for a tracer density.
!> @param[in]     nlayers      Number of layers
!> @param[in,out] adv          Advective update field to apply monotonicity to
!> @param[in]     theta        Field to be advected, used to compute the advective output
!> @param[in]     stencil_size Size of the stencil (number of cells)
!> @param[in]     stencil_map  Dofmaps for the stencil
!> @param[in]     inv_mt       Lumped inverse of the Wtheta mass matrix
!> @param[in]     dt           Timestep
!> @param[in]     order        Desired polynomial order used to reconstruct the high order
!!                             theta approximation
!> @param[in]     ndf_wt       Number of degrees of freedom per cell
!> @param[in]     undf_wt      Number of unique degrees of freedom for the tracer field
!> @param[in]     map_wt       Dofmap for the tracer field
subroutine monotonic_update_code( nlayers,              &
                                  adv,                  &
                                  theta,                &
                                  stencil_size,         &
                                  stencil_map,          &
                                  inv_mt,               &
                                  dt,                   &
                                  global_order,         &
                                  ndf_wt,               &
                                  undf_wt,              &
                                  map_wt)

  implicit none

  ! Arguments
  integer(kind=i_def),                                 intent(in) :: nlayers
  integer(kind=i_def),                                 intent(in) :: ndf_wt
  integer(kind=i_def),                                 intent(in) :: undf_wt
  integer(kind=i_def),                                 intent(in) :: stencil_size
  integer(kind=i_def), dimension(ndf_wt,stencil_size), intent(in) :: stencil_map
  integer(kind=i_def), dimension(ndf_wt),              intent(in) :: map_wt

  real(kind=r_def), dimension(undf_wt), intent(inout) :: adv
  real(kind=r_def), dimension(undf_wt), intent(in)    :: theta
  real(kind=r_def), dimension(undf_wt), intent(in)    :: inv_mt
  real(kind=r_def),                     intent(in)    :: dt
  integer(kind=i_def),                  intent(in)    :: global_order

  ! Internal variables
  integer(kind=i_def) :: k, m, stencil, cell, p, ijkp, order, pmin, pmax
  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  real(kind=r_def) :: theta0, theta_max, theta_min, a_min, a_max

  ! Compute the offset map for all even orders up to order
  allocate( smap(global_order+3,0:global_order/2) )
  smap(:,:) = 0
  do m = 0,global_order,2
    do stencil = 1,m+3
      smap(stencil,m/2) = - 1 - m/2 + (stencil-1)
    end do
  end do

  ! Enforce monotonicity of the advective update across the stencil
  do k = 0, nlayers
    theta0    = theta(stencil_map(1,1)+k)
    theta_max = theta(stencil_map(1,1)+k)
    theta_min = theta(stencil_map(1,1)+k)
    ! Max min across the horzontal stencil (this assumes a cross stencil which
    ! is incorrect for 2D interpolation but should still be a good guess)
    do cell = 2,stencil_size
      theta_max = max(theta_max,theta(stencil_map(1,cell)+k))
      theta_min = min(theta_min,theta(stencil_map(1,cell)+k))
    end do
    ! Max min across the vertical stencil
    order = min(global_order, min(2*(k-1), 2*(nlayers-1-k)))

    ! These limits should be p=1,order+3k
    ! and ijkp = stencil_map(1,1) + k + smap(p,order/2)
    ! but this caused a segmentation fault so for now use simpler
    ! k-1,k+1 limits
    pmin = max(0,       k-1)
    pmax = min(nlayers, k+1)
    do p = pmin,pmax
      ijkp = stencil_map(1,1) + p
      theta_max = max(theta_max,theta(ijkp))
      theta_min = min(theta_min,theta(ijkp))
    end do
    a_min = (theta0 - theta_max)/(dt*inv_mt(stencil_map(1,1)+k))
    a_max = (theta0 - theta_min)/(dt*inv_mt(stencil_map(1,1)+k))
    adv(map_wt(1)+k) = min(a_max,max(adv(stencil_map(1,1)+k),a_min))
  end do

  deallocate( smap )

end subroutine monotonic_update_code

end module monotonic_update_kernel_mod
