!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Transforms a field on the original mesh to the shifted mesh.
!> @details Performs the matrix-vector calculation to transform from a field on the
!> original mesh to its counterpart on the shifted mesh.
!> It multiplies the vector representing the original field by two matrices:
!> T_L, which represents the transform integrals over the lower-half levels of
!> the mesh, and T_U which represents the transform integrals over the upper-half
!> levels of the mesh.
!> With x the vector for the original field, y which represents the new shifted
!> field is calculated from
!> y[i] = T_L[i,j]*x[j] + T_U[i,j-1]*x[j-1]
!>

module matrix_vector_shifted_kernel_mod

use argument_mod,            only : arg_type,                 &
                                    GH_FIELD, GH_OPERATOR,    &
                                    GH_REAL, GH_READ, GH_INC, &
                                    ANY_SPACE_1, ANY_SPACE_2, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def
use kernel_mod,              only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: matrix_vector_shifted_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                     &
       arg_type(GH_FIELD,    GH_REAL, GH_INC,  ANY_SPACE_1),              & ! y_sh
       arg_type(GH_FIELD,    GH_REAL, GH_READ, ANY_SPACE_2),              & ! x
       arg_type(GH_OPERATOR, GH_REAL, GH_READ, ANY_SPACE_1, ANY_SPACE_2), & ! T_L
       arg_type(GH_OPERATOR, GH_REAL, GH_READ, ANY_SPACE_1, ANY_SPACE_2)  & ! T_U
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: matrix_vector_shifted_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: matrix_vector_shifted_code

contains

!> @brief Computes y_sh[i] = T_L[i,j]*x[j] + T_U[i,j-1]*x[j-1]
!> @param[in] cell Horizontal cell index
!> @param[in] nlayers Number of layers
!> @param[in] x Input data from original space
!> @param[in,out] y_sh Output lhs in shifted space
!> @param[in] ncell_3d_L Total number of cells in original mesh
!> @param[in] T_L Local matrix assembly of operator T_L which goes over lower half cells.
!> @param[in] ncell_3d_U Total number of cells in original mesh
!> @param[in] T_U Local matrix assembly of operator T_U which goes over upper half cells
!> @param[in] ndf_sh Number of degrees of freedom per cell for the shifted field
!> @param[in] undf_sh Unique number of degrees of freedom for the shifted field
!> @param[in] map_sh Dofmap for the cell at the base of the column for the shifted field.
!> @param[in] ndf_orig Number of degrees of freedom per cell for the original field
!> @param[in] undf_orig Unique number of degrees of freedom  for the original field
!> @param[in] map_orig Dofmap for the cell at the base of the column for the original field
subroutine matrix_vector_shifted_code(                     &
                                        cell,              &
                                        nlayers_sh,        &
                                        y_sh,              &
                                        x,                 &
                                        ncell_3d_L,        &
                                        T_L,               &
                                        ncell_3d_U,        &
                                        T_U,               &
                                        ndf_sh,            &
                                        undf_sh,           &
                                        map_sh,            &
                                        ndf_orig,          &
                                        undf_orig,         &
                                        map_orig           &
                                      )

  implicit none

  ! Arguments
  integer(kind=i_def),                      intent(in) :: cell, nlayers_sh
  integer(kind=i_def),                      intent(in) :: ncell_3d_L, ncell_3d_U
  integer(kind=i_def),                      intent(in) :: undf_orig, ndf_orig
  integer(kind=i_def),                      intent(in) :: undf_sh, ndf_sh
  integer(kind=i_def), dimension(ndf_orig), intent(in) :: map_orig
  integer(kind=i_def), dimension(ndf_sh),   intent(in) :: map_sh

  real(kind=r_def), dimension(undf_orig),                 intent(in)  :: x
  real(kind=r_def), dimension(undf_sh),                 intent(inout) :: y_sh
  real(kind=r_def), dimension(ndf_sh,ndf_orig,ncell_3d_L), intent(in) :: T_L
  real(kind=r_def), dimension(ndf_sh,ndf_orig,ncell_3d_U), intent(in) :: T_U

  ! Internal variables
  integer(kind=i_def)                   :: df, k, ik
  real(kind=r_def), dimension(ndf_orig) :: x_e
  real(kind=r_def), dimension(ndf_sh)   :: y_sh_e

  ! Bottom layer first, k = 0 with only contribution from T_L
  do df = 1, ndf_orig
    x_e(df) = x(map_orig(df))
  end do

  ik = (cell-1)*(nlayers_sh-1) + 1
  y_sh_e = matmul(T_L(:,:,ik), x_e)

  do df = 1, ndf_sh
    y_sh(map_sh(df)) = y_sh(map_sh(df)) + y_sh_e(df)
  end do

  ! Do values for generic internal level
  do k = 1, nlayers_sh-2

    ! Since x_e is now correct for the upper level calculation,
    ! we reuse it from last iteration
    y_sh_e = matmul(T_U(:,:,ik), x_e)

    do df = 1, ndf_sh
      y_sh(map_sh(df)+k) = y_sh(map_sh(df)+k) + y_sh_e(df)
    end do

    ! Get x values
    do df = 1, ndf_orig
      x_e(df) = x(map_orig(df)+k)
    end do

    ik = (cell-1)*(nlayers_sh-1) + k + 1
    y_sh_e = matmul(T_L(:,:,ik), x_e)
    do df = 1, ndf_sh
       y_sh(map_sh(df)+k) = y_sh(map_sh(df)+k) + y_sh_e(df)
    end do
  end do

  ! Do top level, with only contribution from T_U. Use previous x_e again.
  y_sh_e = matmul(T_U(:,:,ik), x_e)

  do df = 1, ndf_sh
    y_sh(map_sh(df)+nlayers_sh-1) = y_sh(map_sh(df)+nlayers_sh-1) + y_sh_e(df)
  end do


end subroutine matrix_vector_shifted_code

end module matrix_vector_shifted_kernel_mod
