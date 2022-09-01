!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!> @brief Kernel which computes horizontal cell edges value based on the monotone
!!        Koren scheme.
!> @details The kernel computes the edge values at W2-points from a W3-tracer field.
!!  Ref.   Koren, B. (1993), "A robust upwind discretisation method for advection,
!!         diffusion and source terms", in Vreugdenhil, C.B.; Koren, B. (eds.),
!!         Numerical Methods for Advection/Diffusion Problems, Braunschweig:
!!         Vieweg, p. 117, ISBN 3-528-07645-3.
module polyh_w3_koren_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_INC, GH_READ,             &
                              STENCIL, CROSS, GH_BASIS,    &
                              CELL_COLUMN, GH_EVALUATOR,   &
                              outward_normals_to_horizontal_faces
use constants_mod,     only : r_def, i_def, tiny_eps
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: polyh_w3_koren_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,   W2),                        & ! Flux
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! Wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(CROSS)),        & ! Density
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              & ! ndata
       /)
  type(func_type) :: meta_funcs(1) = (/                                   &
       func_type(W2, GH_BASIS)                                            &
       /)
  type(reference_element_data_type) :: meta_reference_element(1) = (/     &
       reference_element_data_type( outward_normals_to_horizontal_faces ) &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: polyh_w3_koren_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: polyh_w3_koren_code

contains

!> @brief Computes the horizontal reconstruction for a tracer field.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] reconstruction Reconstructed W2 field to compute
!> @param[in]     wind           Wind field
!> @param[in]     tracer         Tracer field
!> @param[in]     stencil_size   Size of the stencil (number of cells)
!> @param[in]     stencil_map    Dofmaps for the stencil
!> @param[in]     ndata          Number of data points per dof location
!> @param[in]     ndf_w2         Number of degrees of freedom per cell
!> @param[in]     undf_w2        Number of unique degrees of freedom for the
!!                               reconstruction & wind fields
!> @param[in]     map_w2         Dofmap for the cell at the base of the column
!> @param[in]     basis_w2       Basis function array evaluated at w2 nodes
!> @param[in]     ndf_w3         Number of degrees of freedom per cell
!> @param[in]     undf_w3        Number of unique degrees of freedom for the tracer field
!> @param[in]     map_w3         Dofmap for the cell at the base of the column for the tracer field
!> @param[in]     nfaces_re_h    Number of horizontal neighbours
!> @param[in]     outward_normals_to_horizontal_faces Vector of normals to the
!!                                                    reference element horizontal
!!                                                    "outward faces"
subroutine polyh_w3_koren_code( nlayers,              &
                                reconstruction,       &
                                wind,                 &
                                tracer,               &
                                stencil_size,         &
                                stencil_map,          &
                                ndata,                &
                                ndf_w2,               &
                                undf_w2,              &
                                map_w2,               &
                                basis_w2,             &
                                ndf_w3,               &
                                undf_w3,              &
                                map_w3,               &
                                nfaces_re_h,          &
                                outward_normals_to_horizontal_faces )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: stencil_size
  integer(kind=i_def), intent(in)                    :: nfaces_re_h

  real(kind=r_def), dimension(undf_w2), intent(inout) :: reconstruction
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)    :: tracer
  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2
  integer(kind=i_def), dimension(ndf_w3,stencil_size), intent(in) :: stencil_map
  real(kind=r_def), intent(in) :: outward_normals_to_horizontal_faces(3,nfaces_re_h)

  ! Internal variables
  integer(kind=i_def)                      :: k, df
  real(kind=r_def)                         :: direction
  real(kind=r_def), dimension(nfaces_re_h) :: v_dot_n
  real(kind=r_def)                         :: edge_tracer
  integer(kind=i_def), dimension(3,4) :: point=reshape((/4,1,2,5,1,3,2,1,4,3,1,5/),shape(point))
  real(kind=r_def)                    :: x, y, r, phi, r1, r2

  ! for order = 2 the cross stencil map is
  !      | 5 |
  !  | 2 | 1 | 4 |
  !      | 3 |

  do df = 1,nfaces_re_h
    v_dot_n(df) =  dot_product(basis_w2(:,df,df),outward_normals_to_horizontal_faces(:,df))
  end do

  do k = 0, nlayers - 1
    do df = 1,nfaces_re_h
      ! Check if this is the upwind cell
      direction = wind(map_w2(df) + k )*v_dot_n(df)
      if ( direction > 0.0_r_def ) then
        x = tracer(stencil_map(1,point(2,df))+k) - tracer(stencil_map(1,point(1,df))+k)
        y = tracer(stencil_map(1,point(3,df))+k) - tracer(stencil_map(1,point(2,df))+k)
        r = (y + tiny_eps)/(x + tiny_eps)
        r1 = 2.0_r_def*r
        r2 = ( 1.0_r_def + r1 )/ 3.0_r_def
        phi = max (0.0_r_def, min(r1,r2,2.0_r_def))
        edge_tracer = tracer(stencil_map(1,point(2,df))+k) + 0.5_r_def*phi*x
        reconstruction(map_w2(df) + k ) = edge_tracer
      end if
    end do
  end do

end subroutine polyh_w3_koren_code

end module polyh_w3_koren_kernel_mod
