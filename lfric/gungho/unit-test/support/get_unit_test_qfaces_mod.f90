!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_qfaces_mod

! A module containing a collection of helper routines that provide
! canned basis functions for lowest order function spaces with 3x3
! quadrature points on faces. These are for use when writing unit
! tests of kernels.
! Typically, data can be obtained for the side (horizontal) faces,
! the top and botton (vertical) faces, or both.
! The canned data was computed by the LFRic infrastructure code.
!
  use constants_mod, only : r_def, i_def, l_def

  implicit none

  private

  ! Number of quadrature points per face
  integer(i_def) :: nqp = 9

! Function for sizes

  public :: get_number_quadrature_points_per_face

! Functions for weights
  public :: get_quadrature_faces_h_weights,          &
            get_quadrature_faces_v_weights,          &
            get_quadrature_faces_hv_weights

! Functions for basis functions
! Basis function data exists only for existing tests implied by the
! list of public subroutines. If new options are required the
! necessary data will need to be computed and added to the appropriate
! private subroutine

! Getters for basis functions
! Key: get_[fn space]_[cell shape]_[horizontal and/or vertical]_basis
  public :: get_w2_qfaces_cube_h_basis,              &
            get_w3_qfaces_cube_h_basis,              &
            get_wtheta_qfaces_cube_h_basis,          &
            get_w2trace_qfaces_cube_hv_basis,        &
            get_w2trace_qfaces_cube_h_basis,         &
            get_w2trace_qfaces_cube_v_basis,         &
            get_w2broken_qfaces_cube_hv_basis,       &
            get_w2broken_qfaces_cube_h_basis,        &
            get_w2broken_qfaces_cube_v_basis

contains

  ! Number of quadrature points per face
  subroutine get_number_quadrature_points_per_face(points_per_face)

    implicit none

    integer :: points_per_face

    points_per_face = nqp

  end subroutine get_number_quadrature_points_per_face

!------------------------------------------------------------------
  ! Return the number of quadrature points per face and the weights
  ! of each point per face

  ! Separate routines are provided for horizontal, vertical or both
  ! sets of faces.

  ! Horizontal
  subroutine get_quadrature_faces_h_weights(wqp)

    implicit none

    real(r_def), allocatable, intent(out)  :: wqp(:,:)

    integer :: nfaces = 4

    allocate( wqp( nqp, nfaces) )

    call get_quadrature_faces_weights( wqp, .true., .false.)

  end subroutine get_quadrature_faces_h_weights

  ! Vertical
  subroutine get_quadrature_faces_v_weights(wqp)

    implicit none

    real(r_def), allocatable, intent(out)  :: wqp(:,:)

    integer :: nfaces = 2

    allocate( wqp( nqp, nfaces) )

    call get_quadrature_faces_weights( wqp, .false., .true.)

  end subroutine get_quadrature_faces_v_weights

  ! Both horizontal and vertical
  subroutine get_quadrature_faces_hv_weights(wqp)

    implicit none

    real(r_def), allocatable, intent(out)  :: wqp(:,:)

    integer :: nfaces = 6

    allocate( wqp( nqp, nfaces) )

    call get_quadrature_faces_weights( wqp, .true., .true.)

  end subroutine get_quadrature_faces_hv_weights

!---------------------------------------------------------------------
  ! Return the basis functions for a field on a W2 function space.
  ! The functions are computed on a 3x3 grid on each face of the gungho
  ! cube reference cell.

  ! Separate routines are provided for horizontal, vertical or both
  ! sets of faces where available

  ! Horizontal faces
  subroutine get_w2_qfaces_cube_h_basis(basis_w2_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2_face(:,:,:,:)

    call get_w2_qfaces_cube_basis(basis_w2_face, .true., .false.)
  end subroutine get_w2_qfaces_cube_h_basis

!---------------------------------------------------------------------
  ! Return the basis functions for a field on a W3 function space.
  ! The functions are computed on a 3x3 grid on each face of the gungho
  ! cube reference cell.

  ! Separate routines are provided for horizontal, vertical or both
  ! sets of faces where available

  ! Horizontal faces
  subroutine get_w3_qfaces_cube_h_basis(basis_w3_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w3_face(:,:,:,:)

    call get_w3_qfaces_cube_basis(basis_w3_face, .true., .false.)
  end subroutine get_w3_qfaces_cube_h_basis

!---------------------------------------------------------------------
  ! Return the basis functions for a field on a Wtheta function space.
  ! The functions are computed on a 3x3 grid on each face of the gungho
  ! cube reference cell.

  ! Separate routines are provided for horizontal, vertical or both
  ! sets of faces where available

  ! Horizontal faces
  subroutine get_wtheta_qfaces_cube_h_basis(basis_wtheta_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_wtheta_face(:,:,:,:)

    call get_wtheta_qfaces_cube_basis(basis_wtheta_face, .true., .false.)
  end subroutine get_wtheta_qfaces_cube_h_basis


!---------------------------------------------------------------------
  ! Return the basis functions for a field on a W2trace function space.
  ! The functions are computed on a 3x3 grid on each face of the gungho
  ! cube reference cell.

  ! Separate routines are provided for horizontal, vertical or both
  ! sets of faces.

  ! Horizontal faces
  subroutine get_w2trace_qfaces_cube_h_basis(basis_w2trace_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2trace_face(:,:,:,:)

    call get_w2trace_qfaces_cube_basis(basis_w2trace_face, .true., .false.)
  end subroutine get_w2trace_qfaces_cube_h_basis

  ! Vertical faces
  subroutine get_w2trace_qfaces_cube_v_basis(basis_w2trace_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2trace_face(:,:,:,:)

    call get_w2trace_qfaces_cube_basis(basis_w2trace_face, .false., .true.)
  end subroutine get_w2trace_qfaces_cube_v_basis

  ! Both horizontal and vertical faces
  subroutine get_w2trace_qfaces_cube_hv_basis(basis_w2trace_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2trace_face(:,:,:,:)

    call get_w2trace_qfaces_cube_basis(basis_w2trace_face, .true., .true.)
  end subroutine get_w2trace_qfaces_cube_hv_basis

!---------------------------------------------------------------------
  ! Return the basis functions for a field on a W2broken function space.
  ! The functions are computed on a 3x3 grid on each face of the gungho
  ! cube reference cell.

  ! Separate routines are provided for horizontal, vertical or both
  ! sets of faces.

  ! Horizontal faces
  subroutine get_w2broken_qfaces_cube_h_basis(basis_w2broken_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2broken_face(:,:,:,:)

    call get_w2broken_qfaces_cube_basis(basis_w2broken_face, .true., .false.)
  end subroutine get_w2broken_qfaces_cube_h_basis

  ! Vertical faces
  subroutine get_w2broken_qfaces_cube_v_basis(basis_w2broken_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2broken_face(:,:,:,:)

    call get_w2broken_qfaces_cube_basis(basis_w2broken_face, .false., .true.)
  end subroutine get_w2broken_qfaces_cube_v_basis

  ! Both horizontal and vertical faces
  subroutine get_w2broken_qfaces_cube_hv_basis(basis_w2broken_face)

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2broken_face(:,:,:,:)

    call get_w2broken_qfaces_cube_basis(basis_w2broken_face, .true., .true.)
  end subroutine get_w2broken_qfaces_cube_hv_basis

! PUBLIC subroutines above
!---------------------------------------------------
! PRIVATE subroutines below

! Private routines with the actual data. Only data for the public
! getter routines is provided. If the required public getter routine
! is unavailable then it implies you will need to add the public
! function above and also add data below into an existing or new
! private subroutine.

  subroutine get_quadrature_faces_weights(wqp, horizontal, vertical)
    ! Return the number of quadrature points per face and the weights
    ! of each point per face
    implicit none

    real(r_def),allocatable, intent(out)  :: wqp(:,:)
    logical(l_def),           intent(in)   :: horizontal
    logical(l_def),           intent(in)   :: vertical

    integer(i_def) :: nfaces
    integer(i_def) :: face

    nfaces = 0

    if (horizontal) nfaces = nfaces + 4
    if (vertical) nfaces = nfaces + 2

    allocate(wqp(nqp,nfaces))

    do face = 1, nfaces
      wqp(:,face) = (/   &
        0.077160493827160_r_def, 0.123456790123456_r_def, 0.077160493827160_r_def, &
        0.123456790123456_r_def, 0.197530864197531_r_def, 0.123456790123456_r_def, &
        0.077160493827160_r_def, 0.123456790123456_r_def, 0.077160493827160_r_def /)
    end do
  end subroutine get_quadrature_faces_weights

  subroutine get_w2_qfaces_cube_basis(basis_w2_face, horizontal, vertical)
    ! Return the basis functions for a field on a W2 function space.
    ! The functions are computed on a 3x3 grid on each face of the gungho
    ! cube reference cell.

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2_face(:,:,:,:)
    logical(l_def),           intent(in)   :: horizontal
    logical(l_def),           intent(in)   :: vertical

    integer(i_def) :: nfaces
    integer(i_def) :: face

    nfaces = 0

    if (horizontal) nfaces = nfaces + 4
    if (vertical) nfaces = nfaces + 2

    ! 3-dimension vector space. 6 basis functions, 9 quadrature points.
    allocate(basis_w2_face(3, 6, 9, nfaces))

    face=1

    if (horizontal) then
      basis_w2_face(:,:,:,:) = reshape( [   &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def, &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def   &
       ], [   3,   6,   9,   4] )
      face = face + 4
    end if
    if (vertical) then
      ! No data required yet
    end if
  end subroutine get_w2_qfaces_cube_basis

!---------------------------------------------------------------------

  subroutine get_w3_qfaces_cube_basis(basis_w3_face, horizontal, vertical)
    ! Return the basis functions for a field on a W3 function space.
    ! The functions are computed on a 3x3 grid on each face of the gungho
    ! cube reference cell.

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w3_face(:,:,:,:)
    logical(l_def),           intent(in)   :: horizontal
    logical(l_def),           intent(in)   :: vertical

    integer(i_def) :: nfaces
    integer(i_def) :: face

    nfaces = 0

    if (horizontal) nfaces = nfaces + 4
    if (vertical) nfaces = nfaces + 2

    ! 3-dimension vector space. 1 basis functions, 9 quadrature points.
    allocate(basis_w3_face(1, 1, 9, nfaces))

    face=1
    if (horizontal) then
      basis_w3_face(:,:,:,:) = reshape( [   &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def,  &
        1.000000000000000_r_def,  1.000000000000000_r_def,  1.000000000000000_r_def   &
      ], [   1,   1,   9,   4] )
      face = face + 4
    end if
    if (vertical) then
      ! No data required yet
    end if
  end subroutine get_w3_qfaces_cube_basis

!---------------------------------------------------------------------

  subroutine get_wtheta_qfaces_cube_basis(basis_wtheta_face, horizontal, vertical)
    ! Return the basis functions for a field on a Wtheta function space.
    ! The functions are computed on a 3x3 grid on each face of the gungho
    ! cube reference cell.

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_wtheta_face(:,:,:,:)
    logical(l_def),           intent(in)   :: horizontal
    logical(l_def),           intent(in)   :: vertical

    integer(i_def) :: nfaces
    integer(i_def) :: face

    nfaces = 0

    if (horizontal) nfaces = nfaces + 4
    if (vertical) nfaces = nfaces + 2

    ! 1-dimension scalar space. 2 basis functions, 9 quadrature points.
    allocate(basis_wtheta_face(1, 2, 9, nfaces))

    face=1
    if (horizontal) then
      basis_wtheta_face(:,:,:,:)  = reshape( [   &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
        0.887298334620742_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
        0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def   &
        ], [   1,   2,   9,   4] )

      face = face + 4
    end if
    if (vertical) then
      ! No data required yet
    end if
  end subroutine get_wtheta_qfaces_cube_basis


  subroutine get_w2trace_qfaces_cube_basis(basis_w2trace_face, horizontal, vertical)
    ! Return the basis functions for a field on a W2 function space.
    ! The functions are computed on a 3x3 grid on each face of the gungho
    ! cube reference cell.

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2trace_face(:,:,:,:)
    logical(l_def),           intent(in)   :: horizontal
    logical(l_def),           intent(in)   :: vertical

    integer(i_def) :: nfaces
    integer(i_def) :: face

    nfaces = 0

    if (horizontal) nfaces = nfaces + 4
    if (vertical) nfaces = nfaces + 2

    ! 1-dimension scalar space. 6 basis functions, 9 quadrature points.
    allocate(basis_w2trace_face(1, 6, 9, nfaces))

    face=1

    if (horizontal) then
      basis_w2trace_face(:,:,:,1:4) = reshape( [   &
       1.000000000000000_r_def,  0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def,  1.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def,  1.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def,  1.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def,  1.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def,  1.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def,  1.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def,  1.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def,  1.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def,  1.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
      -0.000000000000000_r_def,  0.887298334620742_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
      -0.000000000000000_r_def,  0.887298334620742_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
      -0.000000000000000_r_def,  0.887298334620742_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
      -0.000000000000000_r_def,  0.500000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
      -0.000000000000000_r_def,  0.500000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
      -0.000000000000000_r_def,  0.500000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
      -0.000000000000000_r_def,  0.112701665379258_r_def,  1.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
      -0.000000000000000_r_def,  0.112701665379258_r_def,  1.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
      -0.000000000000000_r_def,  0.112701665379258_r_def,  1.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def, -0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def, -0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def, -0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def, -0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def, -0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def, -0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def, -0.000000000000000_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def, -0.000000000000000_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def, -0.000000000000000_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def   &
       ], [   1,   6,   9,   4] )

      face = face + 4
    end if

    if (vertical) then
      basis_w2trace_face(:,:,:,face:face+1) = reshape( [   &
       0.887298334620742_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.500000000000000_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.112701665379258_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.887298334620742_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.887298334620742_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.500000000000000_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def,  1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.887298334620742_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.500000000000000_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.112701665379258_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.887298334620742_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.500000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.112701665379258_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.887298334620742_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.500000000000000_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.112701665379258_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def, -0.000000000000000_r_def,  1.000000000000000_r_def   &
       ], [   1,   6,   9,   2] )
    end if
  end subroutine get_w2trace_qfaces_cube_basis


  subroutine get_w2broken_qfaces_cube_basis(basis_w2broken_face, horizontal, vertical)
    ! Return the basis functions for a field on a W2 function space.
    ! The functions are computed on a 3x3 grid on each face of the gungho
    ! cube reference cell.

    implicit none

    real(r_def), allocatable, intent(out)  :: basis_w2broken_face(:,:,:,:)
    logical(l_def),           intent(in)   :: horizontal
    logical(l_def),           intent(in)   :: vertical

    integer(i_def) :: nfaces
    integer(i_def) :: face

    nfaces = 0

    if (horizontal) nfaces = nfaces + 4
    if (vertical) nfaces = nfaces + 2

    ! 3-dimension vector space. 6 basis functions, 9 quadrature points.
    allocate(basis_w2broken_face(3, 6, 9, nfaces))

    face=1

    if (horizontal) then
      basis_w2broken_face(:,:,:,1:4) = reshape( [   &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       1.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.500000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def,  0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -1.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.112701665379258_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.887298334620742_r_def   &
       ], [   3,   6,   9,   4] )

      face = face + 4
    end if

    if (vertical) then
      basis_w2broken_face(:,:,:,face:face+1) = reshape( [   &
      0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.500000000000000_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.500000000000000_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def,  &
       0.112701665379258_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.112701665379258_r_def,  0.000000000000000_r_def,  &
       0.887298334620742_r_def,  0.000000000000000_r_def,  0.000000000000000_r_def,  &
       0.000000000000000_r_def, -0.887298334620742_r_def,  0.000000000000000_r_def,  &
      -0.000000000000000_r_def, -0.000000000000000_r_def, -0.000000000000000_r_def,  &
       0.000000000000000_r_def,  0.000000000000000_r_def,  1.000000000000000_r_def   &
       ], [   3,   6,   9,   2] )
    end if
  end subroutine get_w2broken_qfaces_cube_basis

end module get_unit_test_qfaces_mod
