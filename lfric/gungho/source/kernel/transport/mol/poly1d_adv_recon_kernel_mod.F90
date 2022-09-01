!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes horizontal tracer values through fitting a high
!!        order 1D upwind reconstruction.
!> @details Computes the reconstruction for a tracer field using a high order
!!          polynomial fit to the integrated tracer values. The stencil used
!!          for the polynomial is centred on the upwind cell.
!!          A 1D cross stencil centred around the cell is used to build the
!!          polynomial representation of the tracer field.
!!          This method is only valid for lowest order elements.
module poly1d_adv_recon_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_INC, GH_READ,             &
                              STENCIL, CROSS, GH_BASIS,    &
                              CELL_COLUMN, GH_EVALUATOR,   &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              outward_normals_to_horizontal_faces
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W1, W2, Wtheta
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_adv_recon_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,   W1),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  Wtheta, STENCIL(CROSS)),    &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              &
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
  procedure, nopass :: poly1d_adv_recon_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly1d_adv_recon_code

contains

!> @brief Computes the horizontal polynomial interpolation of a tracer.
!> @param[in]  nlayers Number of layers
!> @param[in,out] reconstruction Reconstructed tracer field to compute
!> @param[in]  wind Wind field
!> @param[in]  tracer Pointwise tracer field to reconstruct
!> @param[in]  stencil_size Size of the stencil (number of cells)
!> @param[in]  stencil_map Dofmaps for the stencil
!> @param[in]  coeff Array of polynomial coefficients for interpolation
!> @param[in]  ndata Number of data points per dof location
!> @param[in]  order Desired order of polynomial reconstruction
!> @param[in]  ndf_w1 Number of degrees of freedom per cell
!> @param[in]  undf_w1 Number of unique degrees of freedom for the
!!                     reconstructed field
!> @param[in]  map_w1 Dofmap for the cell at the base of the column
!> @param[in]  ndf_w2 Number of degrees of freedom per cell
!> @param[in]  undf_w2 Number of unique degrees of freedom for the wind field
!> @param[in]  map_w2 Dofmap for the cell at the base of the column
!> @param[in]  basis_w2 Basis function array evaluated at w1 nodes
!> @param[in]  ndf_wt Number of degrees of freedom per cell
!> @param[in]  undf_wt Number of unique degrees of freedom for the tracer field
!> @param[in]  map_wt Dofmap for the cell at the base of the column for the tracer field
!> @param[in]  ndf_c Number of degrees of freedom per cell for the coeff space
!> @param[in]  undf_c Total number of degrees of freedom for the coeff space
!> @param[in]  map_c Dofmap for the coeff space
!> @param[in]  nfaces_re_h Number of horizontal neighbours
!> @param[in]  outward_normals_to_horizontal_faces Vector of normals to the
!!                                                 reference element horizontal
!!                                                 "outward faces"
subroutine poly1d_adv_recon_code( nlayers,              &
                                  reconstruction,       &
                                  wind,                 &
                                  tracer,               &
                                  stencil_size,         &
                                  stencil_map,          &
                                  coeff,                &
                                  ndata,                &
                                  order,                &
                                  ndf_w1,               &
                                  undf_w1,              &
                                  map_w1,               &
                                  ndf_w2,               &
                                  undf_w2,              &
                                  map_w2,               &
                                  basis_w2,             &
                                  ndf_wt,               &
                                  undf_wt,              &
                                  map_wt,               &
                                  ndf_c,                &
                                  undf_c,               &
                                  map_c,                &
                                  nfaces_re_h,          &
                                  outward_normals_to_horizontal_faces )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: ndf_w1
  integer(kind=i_def), intent(in)                    :: undf_w1
  integer(kind=i_def), dimension(ndf_w1), intent(in) :: map_w1
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), intent(in)                    :: ndf_c
  integer(kind=i_def), intent(in)                    :: undf_c
  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: order
  integer(kind=i_def), intent(in)                    :: stencil_size
  integer(kind=i_def), intent(in)                    :: nfaces_re_h

  real(kind=r_def), dimension(undf_w1), intent(inout) :: reconstruction
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: tracer
  real(kind=r_def), dimension(undf_c),  intent(in)    :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w1), intent(in) :: basis_w2

  integer(kind=i_def), dimension(ndf_wt,stencil_size), intent(in) :: stencil_map

  real(kind=r_def), intent(in) :: outward_normals_to_horizontal_faces(:,:)

  ! Internal variables
  integer(kind=i_def)                      :: k, df, p, face, stencil,        &
                                              stencil_depth, depth, face_mod, &
                                              ijp, kp, km
  real(kind=r_def)                         :: direction, v_dot_n
  real(kind=r_def), dimension(0:nlayers)   :: polynomial_tracer

  integer(kind=i_def), dimension(order+1,nfaces_re_h) :: map1d

  ! Compute 1d map from the cross stencil
  ! i.e for order = 2 the stencil map is
  !      | 5 |
  !  | 2 | 1 | 4 |
  !      | 3 |
  ! so map1d is
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! First cell is always the centre cell
  stencil_depth = order/2
  map1d(1,:) = 1
  do face = 1,nfaces_re_h
    depth=1
    face_mod = mod(face+1,2) * stencil_depth
    do stencil = 2,stencil_depth+1
      map1d(stencil+depth-1, face) = stencil + face_mod
      map1d(stencil+depth, face) = stencil + order + face_mod
      depth=depth+1
    end do
  end do

  do df = 1,nfaces_re_h
    polynomial_tracer(:) = 0.0_r_def
    do p = 1, order+1
      ijp = (p - 1 + (df-1)*(order+1)) + map_c(1)
      stencil = map1d(p,df)
      do k = 0, nlayers
        polynomial_tracer(k) = polynomial_tracer(k) &
                             + tracer( stencil_map(1,stencil) + k )*coeff( ijp )
      end do
    end do
    v_dot_n = dot_product(basis_w2(:,df,df),outward_normals_to_horizontal_faces(:,df))
    do k = 0, nlayers
      km = max(0,k-1)
      kp = min(k,nlayers-1)
      direction = (wind(map_w2(df) + kp ) + wind(map_w2(df) + km ))*v_dot_n
      if ( direction >  0.0_r_def ) reconstruction(map_w1(df) + k ) = polynomial_tracer(k)
    end do
  end do

end subroutine poly1d_adv_recon_code

end module poly1d_adv_recon_kernel_mod
