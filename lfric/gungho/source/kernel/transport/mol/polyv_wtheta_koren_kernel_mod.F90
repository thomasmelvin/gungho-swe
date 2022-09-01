!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel which computes vertical cell egdes value based on the monotone
!!        Koren scheme.
!> @details The kernel computes advective increment at Wtheta-points, using
!!          edge-values at w3-points computed with the Koren scheme.
!!  Ref.   Koren, B. (1993), "A robust upwind discretisation method for advection,
!!         diffusion and source terms", in Vreugdenhil, C.B.; Koren, B. (eds.),
!!         Numerical Methods for Advection/Diffusion Problems, Braunschweig:
!!         Vieweg, p. 117, ISBN 3-528-07645-3.

module polyv_wtheta_koren_kernel_mod

use argument_mod,         only : arg_type, GH_FIELD,     &
                                 GH_SCALAR, GH_REAL,     &
                                 GH_INTEGER, GH_LOGICAL, &
                                 GH_READWRITE,           &
                                 GH_READ, CELL_COLUMN
use constants_mod,        only : r_def, i_def, l_def, tiny_eps, EPS
use fs_continuity_mod,    only : W2, Wtheta
use kernel_mod,           only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: polyv_wtheta_koren_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                            &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                                  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: polyv_wtheta_koren_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: polyv_wtheta_koren_code

contains

!> @brief Computes the vertical fluxes for a tracer density.
!> @param[in]     nlayers      Number of layers
!> @param[in,out] advective    Advective update to increment
!> @param[in]     wind         Wind field
!> @param[in]     tracer       Tracer field to advect
!> @param[in]     ndata        Number of data points per dof location
!> @param[in]     logspace     Perform interpolation in log space
!> @param[in]     ndf_wt       Number of degrees of freedom per cell
!> @param[in]     undf_wt      Number of unique degrees of freedom for the tracer field
!> @param[in]     map_wt       Cell dofmaps for the tracer space
!> @param[in]     ndf_w2       Number of degrees of freedom per cell
!> @param[in]     undf_w2      Number of unique degrees of freedom for the flux &
!!                             wind fields
!> @param[in]     map_w2       Dofmap for the cell at the base of the column

subroutine polyv_wtheta_koren_code( nlayers,              &
                                    advective,            &
                                    wind,                 &
                                    tracer,               &
                                    ndata,                &
                                    logspace,             &
                                    ndf_wt,               &
                                    undf_wt,              &
                                    map_wt,               &
                                    ndf_w2,               &
                                    undf_w2,              &
                                    map_w2                )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: ndata

  real(kind=r_def), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: tracer
  logical(kind=l_def),                  intent(in)    :: logspace

  !Internal variables
  real(kind=r_def), dimension(nlayers+1)   :: wind_1d, dtracerdz
  real(kind=r_def), dimension(nlayers)     :: tracer_w3, wind_w3
  real(kind=r_def), dimension(0:nlayers+2) :: tracer_1d
  integer(kind=i_def) :: k, k1, k2, k3
  real(kind=r_def)    :: x, y, r, r1, r2, phi

  !Extract vertical 1d-arrays from global data
  do k=0,nlayers
        wind_1d(k+1) = wind(map_w2(5)+k)
      tracer_1d(k+1) = tracer(map_wt(1)+k)
  end do
  do k = 1,nlayers
     wind_w3(k) = 0.5_r_def *( wind_1d(k) + wind_1d(k+1))
  end do

  ! Add 2 extra points for tracer_1d array outside the boundaries
  ! to avoid treating the boundaries differently from the rest
  ! of the column.
  !
  ! Ideally these extra points can be constructed by
  ! preserving the slope near the boundaries, but this
  ! approach can cause the data to go outside the bounds
  ! of the initial (unextended) data. Therefore we use zero
  ! slope assumption.

  tracer_1d(0) = tracer_1d(1)
  tracer_1d(nlayers+2) = tracer_1d(nlayers+1)

  ! Apply log to tracer_1d if required
  ! If using the logspace option, the tracer is forced to be positive
  if (logspace) then
      do k=0,nlayers+2
         tracer_1d(k) = log(max(EPS,abs(tracer_1d(k))))
      end do
  end if

  ! Compute tracers at w3-points (edges of cells centred around theta-points)
  do k = 1, nlayers
     if ( wind_w3(k) > 0.0_r_def ) then
        k3 = k + 1_i_def
        k2 = k3 - 1_i_def
        k1 = k3 - 2_i_def
     else
        k3 = k
        k2 = k3 + 1_i_def
        k1 = k3 + 2_i_def
     end if
     x = tracer_1d(k2) - tracer_1d(k1)
     y = tracer_1d(k3) - tracer_1d(k2)
     r = (y + tiny_eps)/(x + tiny_eps)
     r1 = 2.0_r_def*r
     r2 = (1.0_r_def + r1)/3.0_r_def
     phi = max(0.0_r_def, min(r1,r2,2.0_r_def))
     tracer_w3(k) = tracer_1d(k2) + 0.5_r_def*phi*x
  end do

  ! Revert from log of tracer_w3 if required
  if (logspace) then
      do k = 1, nlayers
         tracer_w3(k) = exp(tracer_w3(k))
      end do
  end if

  do k = 2, nlayers
     dtracerdz(k) = tracer_w3(k) - tracer_w3(k-1)
  end do

  do k = 1, nlayers - 1
     advective(map_wt(1) + k ) = advective(map_wt(1) + k ) &
                               + wind(map_w2(5) + k )      &
                               * dtracerdz(k+1)
  end do

end subroutine polyv_wtheta_koren_code

end module polyv_wtheta_koren_kernel_mod
