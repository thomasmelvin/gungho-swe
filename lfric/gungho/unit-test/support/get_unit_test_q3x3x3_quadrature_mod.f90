!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_q3x3x3_quadrature_mod
! A module containing a collection of helper routines that provide canned
! quadrature information on a simple 3x3 mesh for use when writing
! unit tests of kernels

  use constants_mod, only : r_def

  implicit none

  private

  public :: get_gaussian_q3x3x3_quadrature_points_xy,  &
            get_gaussian_q3x3x3_quadrature_points_z,   &
            get_gaussian_q3x3x3_quadrature_weights_xy, &
            get_gaussian_q3x3x3_quadrature_weights_z
  contains

!---------------------------------------------------------------------

  subroutine get_gaussian_q3x3x3_quadrature_points_xy(points_xy)
    ! Calculate the x-y co-ords for the points in a 3x3 gaussian quadrature
    implicit none
    real(r_def), allocatable, intent(out) :: points_xy(:,:)

    allocate(points_xy(9,2))
    points_xy = reshape( [0.112701665379258_r_def, 0.112701665379258_r_def, 0.112701665379258_r_def,   &
                            0.500000000000000_r_def, 0.500000000000000_r_def, 0.500000000000000_r_def, &
                            0.887298334620742_r_def, 0.887298334620742_r_def, 0.887298334620742_r_def, &
                          0.112701665379258_r_def, 0.500000000000000_r_def, 0.887298334620742_r_def,   &
                            0.112701665379258_r_def, 0.500000000000000_r_def, 0.887298334620742_r_def, &
                            0.112701665379258_r_def, 0.500000000000000_r_def, 0.887298334620742_r_def], [9,2] )

  end subroutine get_gaussian_q3x3x3_quadrature_points_xy

!---------------------------------------------------------------------

  subroutine get_gaussian_q3x3x3_quadrature_points_z(points_z)
    ! Calculate the z co-ords for the points in a 3x3 gaussian quadrature
    implicit none
    real(r_def), allocatable, intent(out) :: points_z(:)

    allocate(points_z(3))
    points_z = [0.112701665379258_r_def, 0.500000000000000_r_def, 0.887298334620742_r_def]

  end subroutine get_gaussian_q3x3x3_quadrature_points_z

!---------------------------------------------------------------------

  subroutine get_gaussian_q3x3x3_quadrature_weights_xy(weights_xy)
    ! Calculate the x-y weights for the points in a 3x3 gaussian quadrature
    implicit none
    real(r_def), allocatable, intent(out) :: weights_xy(:)

    allocate(weights_xy(9))
    weights_xy = [0.077160493827160_r_def, 0.123456790123456_r_def, 0.077160493827160_r_def, &
                  0.123456790123456_r_def, 0.197530864197531_r_def, 0.123456790123456_r_def, &
                  0.077160493827160_r_def, 0.123456790123456_r_def, 0.077160493827160_r_def]

  end subroutine get_gaussian_q3x3x3_quadrature_weights_xy

!---------------------------------------------------------------------

  subroutine get_gaussian_q3x3x3_quadrature_weights_z(weights_z)
    ! Calculate the z weights for the points in a 3x3 gaussian quadrature
    implicit none
    real(r_def), allocatable, intent(out) :: weights_z(:)

    allocate(weights_z(3))
    weights_z = [0.277777777777776_r_def, 0.444444444444444_r_def, 0.277777777777776_r_def]

  end subroutine get_gaussian_q3x3x3_quadrature_weights_z

!---------------------------------------------------------------------

end module get_unit_test_q3x3x3_quadrature_mod
