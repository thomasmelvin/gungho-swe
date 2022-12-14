!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Copies the accumulated diagonal elements of a mass matrix, as stored in
!>        a field, back to a LMA representation

!> @details The kernel in mm_diagonal_mod can be used to extract the diagonal
!>          elements of a mass matrix and stores them in a field. For discontinuous
!>          fields, this field stores \f$D_i=\sum_e M^{(e)}_{ii}\f$ where
!>          \f$e\f$ are all elements that share the unknown (global) i.
!>          In other words, \f$D_i\f$ is the diagonal of the assembled mass matrix,
!>          whereas \f$M^{(e)}_{ii}\f$ is the diagonal of the local mass matrix. Using the
!>          data in the field \f$D\f$, kernel constructs an assembled diagonal mass matrix
!>          on every element, i.e. \f$\tilde{M}^{(e)}_{ij}= D_i\delta_{ij}\f$.


module mm_diagonal_assembled_kernel_mod

use argument_mod,            only : arg_type,              &
                                    GH_FIELD, GH_OPERATOR, &
                                    GH_READ, GH_WRITE,     &
                                    GH_REAL, ANY_SPACE_1,  &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def
use kernel_mod,              only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: mm_diagonal_assembled_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                     &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_SPACE_1),             &
       arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_1, ANY_SPACE_1) &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: mm_diagonal_assembled_kernel_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: mm_diagonal_assembled_kernel_code
contains

!> @brief Given a field, stores the assembled diagonal of a mass matrix in a LMA
!> @param[in]  cell Horizontal cell index
!> @param[in]  nlayers Number of vertical layers
!> @param[in]  mm_diag Field array containing the assembled diagonal entries
!>                     of the mass matrix
!> @param[in]  ncell_3d Total number of cells
!> @param[in,out] mass_matrix Array holding mass matrix values after kernel execution
!> @param[in]  ndf Number of degrees of freedom per cell
!> @param[in]  undf Unique number of degrees of freedom
!> @param[in]  map Dofmap for the cell at the base of the column

subroutine mm_diagonal_assembled_kernel_code(cell,        &
                                             nlayers,     &
                                             mm_diag,     &
                                             ncell_3d,    &
                                             mass_matrix, &
                                             ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def),                              intent(in)    :: cell, nlayers
  integer(kind=i_def),                              intent(in)    :: ncell_3d
  integer(kind=i_def),                              intent(in)    :: ndf, undf
  real   (kind=r_def), dimension(undf),             intent(in)    :: mm_diag
  real   (kind=r_def), dimension(ndf,ndf,ncell_3d), intent(inout) :: mass_matrix
  integer(kind=i_def), dimension(ndf),              intent(in)    :: map

  ! Internal variables
  integer(kind=i_def) :: df, k, ik

  do k = 0, nlayers-1
    ik = (cell-1)*nlayers + k + 1
    mass_matrix(:,:,ik) = 0.0_r_def
    do df = 1,ndf
      ! Set diagonal values of matrix
      mass_matrix(df,df,ik) = mm_diag(map(df)+k)
    end do
  end do

end subroutine mm_diagonal_assembled_kernel_code

end module mm_diagonal_assembled_kernel_mod
