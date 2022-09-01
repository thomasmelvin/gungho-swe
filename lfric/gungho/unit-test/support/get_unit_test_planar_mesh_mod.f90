!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_planar_mesh_mod

! A module containing a collection of helper routines that provide
! canned data for a simple 3x3 biperiodic planar mesh or a reference
! cube within it. These are for use when writing unit tests of
! kernels.

  use constants_mod, only : r_def, i_def

  implicit none

  private

  integer(i_def), parameter :: ncells = 9

  public :: get_m3x3_adjacent_face,                  &
            get_normals_to_faces,                    &
            get_normals_to_horizontal_faces,         &
            get_normals_to_vertical_faces,           &
            get_outward_normals_to_faces,            &
            get_outward_normals_to_horizontal_faces, &
            get_outward_normals_to_vertical_faces

contains

  subroutine get_m3x3_adjacent_face(adjacent_face)
    ! The adjacent_face array answers the question: if a cell has
    ! edges numbered 1,2,3 and 4, then what does the neighbouring cell
    ! number the shared edge?  For a regular planar mesh the answer is
    ! simple, and the same for all cells.
    ! This routine provides the answer for a 3x3 mesh used by many tests.
    implicit none

    integer(i_def), intent(out), allocatable :: adjacent_face(:,:)

    integer(i_def) :: cell
    integer(i_def) :: nfaces = 4

    allocate ( adjacent_face( nfaces, ncells) )

    do cell = 1, ncells
      ! As there is no relative rotation between cells, all cells
      ! in the mesh have the same information.
      adjacent_face(:,cell) = (/ 3, 4, 1, 2/)
    end do

  end subroutine get_m3x3_adjacent_face

  subroutine get_normals_to_faces(normals_to_faces)
    ! Return the coordinates of normal vector from each face
    ! as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: normals_to_faces(:,:)
    real(r_def),              allocatable :: normals_h(:,:), normals_v(:,:)

    allocate ( normals_to_faces(3, 6) )
    allocate ( normals_h(3, 4), normals_v(3, 2) )

    ! Populate normals to horizontal faces
    call get_normals_to_horizontal_faces( normals_h )
    normals_to_faces(:,1:4) = normals_h
    ! Populate normals to vertical faces
    call get_normals_to_vertical_faces( normals_v )
    normals_to_faces(:,5:6) = normals_v

  end subroutine get_normals_to_faces

  subroutine get_normals_to_horizontal_faces(normals_to_horizontal_faces)
    ! Return the coordinates of normal vector from each horizontal face
    ! as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: normals_to_horizontal_faces(:,:)

    allocate ( normals_to_horizontal_faces(3, 4) )

    normals_to_horizontal_faces(:,1) = (/ 1.0_r_def , 0.0_r_def , 0.0_r_def /)
    normals_to_horizontal_faces(:,2) = (/ 0.0_r_def ,-1.0_r_def , 0.0_r_def /)
    normals_to_horizontal_faces(:,3) = (/ 1.0_r_def , 0.0_r_def , 0.0_r_def /)
    normals_to_horizontal_faces(:,4) = (/ 0.0_r_def ,-1.0_r_def , 0.0_r_def /)

  end subroutine get_normals_to_horizontal_faces

  subroutine get_normals_to_vertical_faces(normals_to_vertical_faces)
    ! Return the coordinates of normal vector from each vertical face
    ! as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: normals_to_vertical_faces(:,:)

    allocate ( normals_to_vertical_faces(3, 2) )

    normals_to_vertical_faces(:,1) = (/ 0.0_r_def , 0.0_r_def , 1.0_r_def /)
    normals_to_vertical_faces(:,2) = (/ 0.0_r_def , 0.0_r_def , 1.0_r_def /)

  end subroutine get_normals_to_vertical_faces

  subroutine get_outward_normals_to_faces(outward_normals_to_faces)
    ! Return the coordinates of a vector pointing outward from a
    ! cell for each face as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: outward_normals_to_faces(:,:)
    real(r_def),              allocatable :: normals_h(:,:), normals_v(:,:)

    allocate ( outward_normals_to_faces(3, 6) )
    allocate ( normals_h(3, 4), normals_v(3, 2) )

    ! Populate outward normals to horizontal faces
    call get_outward_normals_to_horizontal_faces( normals_h )
    outward_normals_to_faces(:,1:4) = normals_h
    ! Populate outward normals to vertical faces
    call get_outward_normals_to_vertical_faces( normals_v )
    outward_normals_to_faces(:,5:6) = normals_v

  end subroutine get_outward_normals_to_faces

  subroutine get_outward_normals_to_horizontal_faces(outward_normals_to_horizontal_faces)
    ! Return the coordinates of a vector pointing outward from a
    ! cell for each horizontal face as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: outward_normals_to_horizontal_faces(:,:)

    allocate ( outward_normals_to_horizontal_faces(3, 4) )

    outward_normals_to_horizontal_faces(:,1) = (/ -1.0_r_def,  0.0_r_def,  0.0_r_def /)
    outward_normals_to_horizontal_faces(:,2) = (/  0.0_r_def, -1.0_r_def,  0.0_r_def /)
    outward_normals_to_horizontal_faces(:,3) = (/  1.0_r_def,  0.0_r_def,  0.0_r_def /)
    outward_normals_to_horizontal_faces(:,4) = (/  0.0_r_def,  1.0_r_def,  0.0_r_def /)

  end subroutine get_outward_normals_to_horizontal_faces

  subroutine get_outward_normals_to_vertical_faces(outward_normals_to_vertical_faces)
    ! Return the coordinates of a vector pointing outward from a
    ! cell for each vertical face as copied from the reference cube
    implicit none

    real(r_def), intent(out), allocatable :: outward_normals_to_vertical_faces(:,:)

    allocate ( outward_normals_to_vertical_faces(3, 2) )

    outward_normals_to_vertical_faces(:,1) = (/  0.0_r_def,  0.0_r_def, -1.0_r_def /)
    outward_normals_to_vertical_faces(:,2) = (/  0.0_r_def,  0.0_r_def,  1.0_r_def /)

  end subroutine get_outward_normals_to_vertical_faces

end module get_unit_test_planar_mesh_mod
