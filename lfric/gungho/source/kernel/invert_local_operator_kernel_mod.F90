!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Provides access to the members of the w0_kernel class.
module invert_local_operator_kernel_mod

  use argument_mod,      only: arg_type,             &
                               GH_OPERATOR, GH_REAL, &
                               GH_READ, GH_WRITE,    &
                               ANY_SPACE_1, CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W3
  use kernel_mod,        only: kernel_type
  use matrix_invert_mod, only: matrix_invert

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: invert_local_operator_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                      &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_1, ANY_SPACE_1), &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  ANY_SPACE_1, ANY_SPACE_1)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: invert_local_operator_code
  end type invert_local_operator_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: invert_local_operator_code

contains

!> @brief This subroutine computes the mass matrix for the W0 space
!! @param[in] cell Cell Number
!! @param[in] nlayers Number of layers
!! @param[in] ncell3d_inv ncell*nlayers
!! @param[in,out] matrix_inv Inverse matrix
!! @param[in] ncell3d ncell*nlayers
!! @param[in] matrix Input matrix
!! @param[in] ndf Number of degrees of freedom per cell
subroutine invert_local_operator_code(cell, nlayers, ncell3d_inv,  &
                                      matrix_inv, ncell3d, matrix, &
                                      ndf)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell
  integer(kind=i_def), intent(in) :: nlayers, ndf
  integer(kind=i_def), intent(in) :: ncell3d, ncell3d_inv

  real(kind=r_def), dimension(ndf,ndf,ncell3d),     intent(in)     :: matrix
  real(kind=r_def), dimension(ndf,ndf,ncell3d_inv), intent(inout)  :: matrix_inv

  ! Internal variables
  integer(kind=i_def) :: k, ik

  ! Loop over layers: Start from 1 as in this loop k is not an offset
  do k = 1, nlayers
    ik = k + (cell-1)*nlayers
    call matrix_invert(matrix(:,:,ik),matrix_inv(:,:,ik),ndf)
  end do

end subroutine invert_local_operator_code

end module invert_local_operator_kernel_mod
