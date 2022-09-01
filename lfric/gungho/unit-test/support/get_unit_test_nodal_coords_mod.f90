!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Helper routines that provide canned nodal coordinates on a simple 3x3 mesh
! for use when writing unit tests of kernels.

module get_unit_test_nodal_coords_mod

  use constants_mod, only : r_def

  implicit none

  private

  public :: get_w0_nodal_coords,       &
            get_w1_nodal_coords,       &
            get_w2_nodal_coords,       &
            get_w2h_nodal_coords,      &
            get_w2v_nodal_coords,      &
            get_w2broken_nodal_coords, &
            get_w2trace_nodal_coords,  &
            get_w3_nodal_coords,       &
            get_wtheta_nodal_coords,   &
            get_wchi_nodal_coords

  contains

    ! Calculates the coordinates of nodes on vertices.
    subroutine get_w0_nodal_coords( coords_w0 )

      implicit none

      real(r_def), allocatable, intent(out) :: coords_w0(:,:)

      allocate(coords_w0(3,8))
      coords_w0 = reshape( [0.00_r_def, 0.00_r_def, 0.00_r_def, &
                            1.00_r_def, 0.00_r_def, 0.00_r_def, &
                            1.00_r_def, 1.00_r_def, 0.00_r_def, &
                            0.00_r_def, 1.00_r_def, 0.00_r_def, &
                            0.00_r_def, 0.00_r_def, 1.00_r_def, &
                            1.00_r_def, 0.00_r_def, 1.00_r_def, &
                            1.00_r_def, 1.00_r_def, 1.00_r_def, &
                            0.00_r_def, 1.00_r_def, 1.00_r_def], [3,8] )

    end subroutine get_w0_nodal_coords

    ! Calculates the coordinates of nodes on edges.
    subroutine get_w1_nodal_coords( coords_w1 )

      implicit none

      real(r_def), allocatable, intent(out) :: coords_w1(:,:)

      allocate(coords_w1(3,12))
      coords_w1 = reshape( [0.00_r_def, 0.50_r_def, 0.00_r_def, &
                            0.50_r_def, 0.00_r_def, 0.00_r_def, &
                            1.00_r_def, 0.50_r_def, 0.00_r_def, &
                            0.50_r_def, 1.00_r_def, 0.00_r_def, &
                            0.00_r_def, 0.00_r_def, 0.50_r_def, &
                            1.00_r_def, 0.00_r_def, 0.50_r_def, &
                            1.00_r_def, 1.00_r_def, 0.50_r_def, &
                            0.00_r_def, 1.00_r_def, 0.50_r_def, &
                            0.00_r_def, 0.50_r_def, 1.00_r_def, &
                            0.50_r_def, 0.00_r_def, 1.00_r_def, &
                            1.00_r_def, 0.50_r_def, 1.00_r_def, &
                            0.50_r_def, 1.00_r_def, 1.00_r_def], [3,12] )

    end subroutine get_w1_nodal_coords

  ! Calculates the coordinates of nodes on faces.
  subroutine get_w2_nodal_coords( coords_w2 )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_w2(:,:)

    allocate(coords_w2(3,6))
    coords_w2 = reshape( [0.00_r_def, 0.50_r_def, 0.50_r_def, &
                          0.50_r_def, 0.00_r_def, 0.50_r_def, &
                          1.00_r_def, 0.50_r_def, 0.50_r_def, &
                          0.50_r_def, 1.00_r_def, 0.50_r_def, &
                          0.50_r_def, 0.50_r_def, 0.00_r_def, &
                          0.50_r_def, 0.50_r_def, 1.00_r_def], [3,6] )

  end subroutine get_w2_nodal_coords

  ! Calculates the coordinates of nodes on vertical faces.
  subroutine get_w2v_nodal_coords( coords_w2v )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_w2v(:,:)

    allocate(coords_w2v(3,2))
    coords_w2v = reshape( [0.50_r_def, 0.50_r_def, 0.00_r_def, &
                           0.50_r_def, 0.50_r_def, 1.00_r_def], [3,2] )

  end subroutine get_w2v_nodal_coords

  ! Calculates the coordinates of nodes on horizontal faces.
  subroutine get_w2h_nodal_coords( coords_w2h )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_w2h(:,:)

    allocate(coords_w2h(3,4))
    coords_w2h = reshape( [0.00_r_def, 0.50_r_def, 0.50_r_def, &
                           0.50_r_def, 0.00_r_def, 0.50_r_def, &
                           1.00_r_def, 0.50_r_def, 0.50_r_def, &
                           0.50_r_def, 1.00_r_def, 0.50_r_def], [3,4] )

  end subroutine get_w2h_nodal_coords

  ! Calculates the coordinates of nodes on broken faces.
  subroutine get_w2broken_nodal_coords( coords_w2broken )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_w2broken(:,:)

    allocate(coords_w2broken(3,6))
    coords_w2broken = reshape( [0.00_r_def, 0.50_r_def, 0.50_r_def, &
                                0.50_r_def, 0.00_r_def, 0.50_r_def, &
                                1.00_r_def, 0.50_r_def, 0.50_r_def, &
                                0.50_r_def, 1.00_r_def, 0.50_r_def, &
                                0.50_r_def, 0.50_r_def, 0.00_r_def, &
                                0.50_r_def, 0.50_r_def, 1.00_r_def], [3,6] )

  end subroutine get_w2broken_nodal_coords

  ! Calculates the coordinates of nodes on trace faces.
  subroutine get_w2trace_nodal_coords( coords_w2trace )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_w2trace(:,:)

    allocate(coords_w2trace(3,6))
    coords_w2trace = reshape( [0.00_r_def, 0.50_r_def, 0.50_r_def, &
                               0.50_r_def, 0.00_r_def, 0.50_r_def, &
                               1.00_r_def, 0.50_r_def, 0.50_r_def, &
                               0.50_r_def, 1.00_r_def, 0.50_r_def, &
                               0.50_r_def, 0.50_r_def, 0.00_r_def, &
                               0.50_r_def, 0.50_r_def, 1.00_r_def], [3,6] )

  end subroutine get_w2trace_nodal_coords

  ! Calculates the coordinates of nodes on centres.
  subroutine get_w3_nodal_coords( coords_w3 )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_w3(:,:)

    allocate(coords_w3(3,1))
    coords_w3 = reshape( [0.50_r_def, 0.50_r_def, 0.50_r_def], [3,1] )

  end subroutine get_w3_nodal_coords

  ! Calculates the coordinates of nodes on theta..
  subroutine get_wtheta_nodal_coords( coords_wtheta )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_wtheta(:,:)

    allocate(coords_wtheta(3,2))
    coords_wtheta = reshape( [0.50_r_def, 0.50_r_def, 0.00_r_def, &
                              0.50_r_def, 0.50_r_def, 1.00_r_def], [3,2] )

  end subroutine get_wtheta_nodal_coords

  ! Calculates the coordinates of nodes on chi.
  subroutine get_wchi_nodal_coords( coords_wchi )

    implicit none

    real(r_def), allocatable, intent(out) :: coords_wchi(:,:)

    allocate(coords_wchi(3,1))
    coords_wchi = reshape( [0.50_r_def, 0.50_r_def, 0.50_r_def], [3,1] )

  end subroutine get_wchi_nodal_coords

end module get_unit_test_nodal_coords_mod
