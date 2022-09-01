!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes exner from equation of state and hydrostatic balance
!> Only used for lowest order/FD fields.
!> Once we are able to loop over an LBC mesh, then this kernel
!> can be replaced with standard balancing kernels that don't
!> require masking.

module lbc_balance_kernel_mod

use argument_mod,               only : arg_type, func_type,      &
                                       GH_FIELD, GH_REAL,        &
                                       GH_SCALAR,                &
                                       GH_READ, GH_WRITE, GH_READWRITE,        &
                                       GH_BASIS, CELL_COLUMN, GH_EVALUATOR
use constants_mod,              only : r_def, i_def
use fs_continuity_mod,          only : Wtheta, W3
use kernel_mod,                 only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: lbc_balance_kernel_type
  private
  type(arg_type) :: meta_args(11) = (/                  &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),     & ! exner
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! rho
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta), & ! theta
       arg_type(GH_FIELD*3, GH_REAL, GH_READ,  Wtheta), & ! moist_dyn_factors
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! height_w3
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! w3_lbc_mask
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),          & ! gravity
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),          & ! p_zero
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),          & ! kappa
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),          & ! rd
       arg_type(GH_SCALAR,  GH_REAL, GH_READ)           & ! cp
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: lbc_balance_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: lbc_balance_code

contains

!> @brief Computes density from equation of state
!! @param[in] nlayers Number of layers
!! @param[in,out] exner Exner pressure field
!! @param[in] rho Density field
!! @param[in] theta Potential temperature field
!! @param[in] height_w3 Height coordinate in w3
!! @param[in] w3_lbc_mask LBC mask for w3 space
!! @param[in] gravity The planet gravity
!! @param[in] p_zero Reference surface pressure
!! @param[in] kappa Ratio of Rd and cp
!! @param[in] rd Gas constant for dry air
!! @param[in] cp Specific heat of dry air at constant pressure
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] ndf_wt Number of degrees of freedom per cell for wtheta
!! @param[in] undf_wt Number of unique degrees of freedom  for wtheta
!! @param[in] map_wt Dofmap for the cell at the base of the column for wt
subroutine lbc_balance_code(nlayers, exner, rho, theta,                           &
                                     moist_dyn_gas, moist_dyn_tot, moist_dyn_fac, &
                                     height_w3, w3_lbc_mask,                      &
                                     gravity, p_zero, kappa, rd, cp,              &
                                     ndf_w3, undf_w3, map_w3,                     &
                                     ndf_wt, undf_wt, map_wt )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf_w3, undf_w3,  ndf_wt, undf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(undf_w3),  intent(inout) :: exner
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: rho
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: height_w3, w3_lbc_mask
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: moist_dyn_gas, moist_dyn_tot, moist_dyn_fac
  real(kind=r_def), dimension(undf_wt),  intent(in) :: theta
  real(kind=r_def),                      intent(in) :: gravity
  real(kind=r_def),                      intent(in) :: p_zero
  real(kind=r_def),                      intent(in) :: kappa
  real(kind=r_def),                      intent(in) :: rd
  real(kind=r_def),                      intent(in) :: cp

  ! Internal variables
  integer(kind=i_def)                  :: k, dft
  real(kind=r_def)                     :: theta_moist, dz

  if (w3_lbc_mask(map_w3(1)) > 0.0_r_def)then

    ! Compute exner from eqn of state in lowest level
    theta_moist = 0.0_r_def
    do dft = 1, 2
      theta_moist = theta_moist + &
         0.5_r_def*theta( map_wt(dft)) * moist_dyn_gas(map_wt(dft))
    end do

    exner(map_w3(1)) = &
       (rd*rho(map_w3(1))*theta_moist/p_zero)**(kappa/(1.0_r_def-kappa))

    ! Exner on other levels from hydrostatic balance
    do k = 1, nlayers-1
      dz = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
      theta_moist = moist_dyn_gas(map_wt(1)+k) * theta(map_wt(1)+k) /   &
         moist_dyn_tot(map_wt(1)+k)
      exner(map_w3(1)+k) = exner(map_w3(1)+k-1) - gravity * dz / (cp * theta_moist)
    end do

  end if

end subroutine lbc_balance_code

end module lbc_balance_kernel_mod
