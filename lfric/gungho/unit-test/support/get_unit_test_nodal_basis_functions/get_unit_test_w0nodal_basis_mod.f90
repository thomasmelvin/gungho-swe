!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_w0nodal_basis_mod
! A module containing a collection of helper routines that provide canned
! basis functions (and differential basis functions) evaluated on w0 nodal
! points for use when writing unit tests of kernels.
!
! Not all the combinations of (diff) basis functions have been added. There are
! many of them and only a few will ever be used, so they are being added as
! necessary.
!
! If you require a (diff) basis function that you think should be in this file
! but isn't, then feel free to add it.

  use constants_mod, only : r_def

  implicit none

  private

  public :: get_wchi_w0nodal_basis

  contains

!---------------------------------------------------------------------

  subroutine get_wchi_w0nodal_basis(basis_w0_on_wchi)
    ! Return the basis function for a field on a w0 function space
    ! evaluated on wchi nodal points

    implicit none

    real(r_def), allocatable :: basis_w0_on_wchi(:,:,:)
    integer :: i

    ! Lowest order scalar. 2 W0 basis functions. 1 WCHI nodal points.
    allocate(basis_w0_on_wchi(1,8,8))

    basis_w0_on_wchi(:,:,:) = 0.0_r_def
    do i=1,8
      basis_w0_on_wchi(1,i,i) = 1.0_r_def
    end do

  end subroutine get_wchi_w0nodal_basis

!---------------------------------------------------------------------

end module get_unit_test_w0nodal_basis_mod
