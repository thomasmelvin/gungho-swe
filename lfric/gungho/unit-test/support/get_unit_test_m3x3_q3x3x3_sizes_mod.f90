!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_m3x3_q3x3x3_sizes_mod
! A module containing a collection of helper routines that provide the sizes of
! canned  dofmaps on a simple 3x3 mesh

  use constants_mod, only : i_def

  implicit none

  private

  public :: get_w0_m3x3_q3x3x3_size,       &
            get_w1_m3x3_q3x3x3_size,       &
            get_w2_m3x3_q3x3x3_size,       &
            get_w2v_m3x3_q3x3x3_size,      &
            get_w2h_m3x3_q3x3x3_size,      &
            get_w2broken_m3x3_q3x3x3_size, &
            get_w2trace_m3x3_q3x3x3_size,  &
            get_w3_m3x3_q3x3x3_size,       &
            get_wtheta_m3x3_q3x3x3_size,   &
            get_wchi_m3x3_q3x3x3_size

  contains

!---------------------------------------------------------------------

  subroutine get_w0_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                      dim_space, dim_space_diff, &
                                      nqp_h, nqp_v, &
                                      nlayers )
    ! For a field on a lowest-order W0 function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature

    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 8
    undf = 9*(nlay+1)
    ncells = 9
    dim_space = 1
    dim_space_diff = 3
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w0_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_w1_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                      dim_space, dim_space_diff, &
                                      nqp_h, nqp_v, &
                                      nlayers )
    ! For a field on a lowest-order W1 function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 12
    undf= 9*nlay + 18*(nlay+1)
    ncells=9
    dim_space = 3
    dim_space_diff = 3
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w1_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_w2_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                      dim_space, dim_space_diff, &
                                      nqp_h, nqp_v, &
                                      nlayers )
    ! For a field on a lowest-order W2 function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 6
    undf= 9*(nlay+1) + 18*nlay
    ncells=9
    dim_space = 3
    dim_space_diff = 1
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w2_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_w2v_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                       dim_space, dim_space_diff, &
                                       nqp_h, nqp_v, &
                                       nlayers )
    ! For a field on a lowest-order W2V function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 2
    undf= 9*(nlay+1)
    ncells=9
    dim_space = 3
    dim_space_diff = 1
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w2v_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_w2h_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                       dim_space, dim_space_diff, &
                                       nqp_h, nqp_v, &
                                       nlayers )
    ! For a field on a lowest-order W2V function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 4
    undf= 18*nlay
    ncells=9
    dim_space = 3
    dim_space_diff = 1
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w2h_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_w2broken_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                            dim_space, dim_space_diff, &
                                            nqp_h, nqp_v, &
                                            nlayers )
    ! For a field on a lowest-order W2broken function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 6
    undf= 54*nlay
    ncells=9
    dim_space = 3
    dim_space_diff = 1
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w2broken_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_w2trace_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                           dim_space, dim_space_diff, &
                                           nqp_h, nqp_v, &
                                           nlayers )
    ! For a field on a lowest-order W2trace function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 6
    undf= 9*(nlay+1) + 18*nlay
    ncells=9
    dim_space = 1
    dim_space_diff = 3
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w2trace_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_w3_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                      dim_space, dim_space_diff, &
                                      nqp_h, nqp_v, &
                                      nlayers )
    ! For a field on a lowest-order W3 function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 1
    undf= 9*nlay
    ncells=9
    dim_space = 1
    dim_space_diff = 3
    nqp_h = 9
    nqp_v = 3

  end subroutine get_w3_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_wtheta_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                          dim_space, dim_space_diff, &
                                          nqp_h, nqp_v, &
                                          nlayers )
    ! For a field on a lowest-order Wtheta function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 2
    undf= 9*(nlay+1)
    ncells=9
    dim_space = 1
    dim_space_diff = 3
    nqp_h = 9
    nqp_v = 3

  end subroutine get_wtheta_m3x3_q3x3x3_size

!---------------------------------------------------------------------

  subroutine get_wchi_m3x3_q3x3x3_size( ndf, undf, ncells, &
                                        dim_space, dim_space_diff, &
                                        nqp_h, nqp_v, &
                                        nlayers )
    ! For a coordinate field on Wchi function space, return:
    !  * ndf: no. of dofs in a cell
    !  * undf: no. of unique dofs in a 3x3 domain
    !  * ncells: no. of cells in a 3x3 domain
    !  * dim_space: no of dimensions in this function space
    !  * dim_space_diff: no. of dims in this funct space when differentiated
    !  * nqp_h: no. of quadrature points in the horiz for a 3x3x3 quadrature
    !  * nqp_v: no. of quadrature points in the vertical for a 3x3x3 quadrature
    ! Note that the coordinate field in GungHo uses next-to-lowest order Wchi,
    ! and this routine returns the numbers it requires rather than lowest-order
    ! numbers.
    implicit none
    integer(i_def), intent(out) :: ndf
    integer(i_def), intent(out) :: undf
    integer(i_def), intent(out) :: ncells
    integer(i_def), intent(out) :: dim_space, dim_space_diff
    integer(i_def), intent(out) :: nqp_h, nqp_v
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: nlay

    if (present(nlayers))then
      nlay=nlayers
    else
      nlay=3
    end if

    ndf = 8
    undf= 72*nlay
    ncells=9
    dim_space = 1
    dim_space_diff = 3
    nqp_h = 9
    nqp_v = 3

  end subroutine get_wchi_m3x3_q3x3x3_size

!---------------------------------------------------------------------

end module get_unit_test_m3x3_q3x3x3_sizes_mod
