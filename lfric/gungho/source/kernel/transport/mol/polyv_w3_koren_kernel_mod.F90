!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel which computes vertical cell egdes value based on the monotone
!!        Koren scheme.
!> @details The kernel computes vertical cell egdes at W2-points from
!!          W3-field using the Koren scheme.
!!  Ref.   Koren, B. (1993), "A robust upwind discretisation method for advection,
!!         diffusion and source terms", in Vreugdenhil, C.B.; Koren, B. (eds.),
!!         Numerical Methods for Advection/Diffusion Problems, Braunschweig:
!!         Vieweg, p. 117, ISBN 3-528-07645-3.

module polyv_w3_koren_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_LOGICAL,                  &
                              GH_INC, GH_READ, GH_BASIS,   &
                              CELL_COLUMN, GH_EVALUATOR
use constants_mod,     only : r_def, i_def, l_def, tiny_eps, EPS
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: polyv_w3_koren_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,   W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              &
       /)
  type(func_type) :: meta_funcs(1) = (/                                 &
       func_type(W2, GH_BASIS)                                          &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: polyv_w3_koren_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: polyv_w3_koren_code

contains

!> @brief Computes the vertical reconstructions for a tracer.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] reconstruction Mass reconstruction field to compute
!> @param[in]     wind           Wind field
!> @param[in]     tracer         Tracer field
!> @param[in]     ndata          Number of data points per dof location
!> @param[in]     logspace       Perform interpolation in log space;
!> @param[in]     ndf_w2         Number of degrees of freedom per cell
!> @param[in]     undf_w2        Number of unique degrees of freedom for the reconstruction &
!!                               wind fields
!> @param[in]     map_w2         Dofmap for the cell at the base of the column
!> @param[in]     basis_w2       Basis function array evaluated at w2 nodes
!> @param[in]     ndf_w3         Number of degrees of freedom per cell
!> @param[in]     undf_w3        Number of unique degrees of freedom for the tracer field
!> @param[in]     map_w3         Cell dofmaps for the tracer space

subroutine polyv_w3_koren_code( nlayers,                           &
                                reconstruction,                    &
                                wind,                              &
                                tracer,                            &
                                ndata,                             &
                                logspace,                          &
                                ndf_w2,                            &
                                undf_w2,                           &
                                map_w2,                            &
                                basis_w2,                          &
                                ndf_w3,                            &
                                undf_w3,                           &
                                map_w3                             )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(undf_w2), intent(inout) :: reconstruction
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)    :: tracer
  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2

  logical(kind=l_def), intent(in) :: logspace

  ! Local variables
  integer(kind=i_def)                             :: k, k1, k2, k3
  real(kind=r_def)                                :: x, y, r, r1, r2, phi
  integer(kind=i_def), parameter                  :: ext = 2
  real(kind=r_def), dimension(nlayers+1)          :: wind_1d, tracer_edge
  real(kind=r_def), dimension(1-ext:nlayers+ext)  :: tracer_1d

  ! Extract the global data into 1d-array
  do k = 0,nlayers
    wind_1d(k+1) = wind(map_w2(5)+k)
  end do
  do k = 0,nlayers - 1
    tracer_1d(k+1) = tracer(map_w3(1)+k)
  end do

  ! Extend the tracer_1d array to go beyond the boundaries
  ! to avoid treating the boundaries differently from the rest
  ! of the column.
  !
  ! Ideally these extra points can be constructed by
  ! preserving the slope near the boundaries, but this
  ! approach can cause the data to go outside the initial bounds.
  ! Therefore we use zero slope assumption.

  tracer_1d(1-ext:0) = tracer_1d(1)
  tracer_1d(nlayers+1:nlayers+ext) = tracer_1d(nlayers)

  ! Apply log to tracer_1d if required
  ! If using the logspace option, the tracer is forced to be positive
  if (logspace) then
      do k = 1-ext,nlayers+ext
         tracer_1d(k) = log(max(EPS,abs(tracer_1d(k))))
      end do
  end if

  do k = 1, nlayers + 1
     if ( wind_1d(k) > 0.0_r_def ) then
        k3 = k
        k2 = k3 - 1_i_def
        k1 = k3 - 2_i_def
     else
        k3 = k  - 1_i_def
        k2 = k3 + 1_i_def
        k1 = k3 + 2_i_def
     end if

     x = tracer_1d(k2) - tracer_1d(k1)
     y = tracer_1d(k3) - tracer_1d(k2)
     r = (y + tiny_eps)/(x + tiny_eps)
     r1 = 2.0_r_def*r
     r2 = (1.0_r_def + r1)/3.0_r_def
     phi = max(0.0_r_def, min(r1,r2,2.0_r_def))
     tracer_edge(k) = tracer_1d(k2) + 0.5_r_def*phi*x
  end do

  if (logspace) then
      do k=1,nlayers+1
         tracer_edge(k) = exp(tracer_edge(k))
      end do
  end if

  do k = 0, nlayers
     reconstruction(map_w2(5) + k ) =  tracer_edge(k+1)
  end do

end subroutine polyv_w3_koren_code

end module polyv_w3_koren_kernel_mod
