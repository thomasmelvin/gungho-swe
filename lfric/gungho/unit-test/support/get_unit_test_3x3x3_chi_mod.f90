!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Simple module/kernel to generate coordinate fields and return them for use in
!> unit tests

module get_unit_test_3x3x3_chi_mod

  use constants_mod,                       only : i_def, r_def

  implicit none

  private

  public :: get_w0_3x3x3_field
  public :: get_wchi_3x3x3_field
  public :: get_wchi_3x3x3_latlon_field
  contains

!---------------------------------------------------------------------

  subroutine get_w0_3x3x3_field(chi1, chi2, chi3, dx, dy, dz, &
                                dofmap, nlayers)
    !external imports

    implicit none

    ! For a field on a lowest-order W0 function space:
    !  - the dimensional scaling factors (dx, dy, dz) are passed as input
    !  - the dofmap is passed as input
    !  - the number of layers in the z-direction is passed as input
    !  - a 3-dimensional (chi1,chi2,chi3) coordinate field is returned
    real(r_def), intent(out) :: chi1(:), chi2(:), chi3(:)
    real(r_def), intent(in) :: dx, dy, dz

    integer(i_def), intent(in) :: dofmap(:,:)
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: cell, i, j, k, nlay

    nlay = 3
    if (present(nlayers))then
      nlay=nlayers
    end if

    ! Compute coordinates
    cell = 1
    do j = 1,3
      do i = 1,3
        do k = 0,nlay
          chi1(dofmap(1,cell)+k) = real(i-1)*dx
          chi2(dofmap(1,cell)+k) = real(j-1)*dy
          chi3(dofmap(1,cell)+k) = real(k)*dz
        end do
        cell = cell + 1
      end do
    end do

  end subroutine get_w0_3x3x3_field

  subroutine get_wchi_3x3x3_field(chi1, chi2, chi3, dx, dy, dz, &
                                  dzdx, dofmap, nlayers)
    !external imports

    implicit none

    ! For a field on WChi function space:
    !  - the dimensional scaling factors (dx, dy, dz) are passed as input
    !  - a gradient in the x-direction (dxdz) is passed as input
    !  - the dofmap is passed as input
    !  - the number of layers in the z-direction is passed as input
    !  - a 3-dimensional (chi1,chi2,chi3) coordinate field is returned
    real(r_def), intent(out) :: chi1(:), chi2(:), chi3(:)
    real(r_def), intent(in) :: dx, dy, dz, dzdx

    integer(i_def), intent(in) :: dofmap(:,:)
    integer(i_def), optional, intent(in) :: nlayers

    integer(i_def) :: cell, i, j, k, nlay

    nlay = 2
    if (present(nlayers))then
      nlay=nlayers-1
    end if

    ! Compute coordinates, for P1DG case
    ! Uniform mesh of size (dx,dy,dz) with constant slope (dxdz) in the x-direction
    ! Note, since this is a DG1 field the dof locations do not correspond to the
    ! vertex locations and so the reference element indices cannot be used
    cell = 1
    do j = 1,3
      do i = 1,3
        do k = 0,nlay
          chi1(dofmap(1:7:2, cell) + k) = real(i-1, r_def)*dx
          chi1(dofmap(2:8:2, cell) + k) = real(i, r_def)*dx

          chi2(dofmap(1:2, cell) + k) = real(j-1, r_def)*dy
          chi2(dofmap(5:6, cell) + k) = real(j-1, r_def)*dy
          chi2(dofmap(3:4, cell) + k) = real(j, r_def)*dy
          chi2(dofmap(7:8, cell) + k) = real(j, r_def)*dy

          chi3(dofmap(1, cell) + k) = real(k, r_def)*dz + real(i-1,r_def)*dzdx
          chi3(dofmap(2, cell) + k) = real(k, r_def)*dz + real(i  ,r_def)*dzdx
          chi3(dofmap(3, cell) + k) = real(k, r_def)*dz + real(i-1,r_def)*dzdx
          chi3(dofmap(4, cell) + k) = real(k, r_def)*dz + real(i  ,r_def)*dzdx
          chi3(dofmap(5, cell) + k) = real(k+1, r_def)*dz + real(i-1,r_def)*dzdx
          chi3(dofmap(6, cell) + k) = real(k+1, r_def)*dz + real(i  ,r_def)*dzdx
          chi3(dofmap(7, cell) + k) = real(k+1, r_def)*dz + real(i-1,r_def)*dzdx
          chi3(dofmap(8, cell) + k) = real(k+1, r_def)*dz + real(i  ,r_def)*dzdx
        end do
        cell = cell + 1
      end do
    end do

  end subroutine get_wchi_3x3x3_field

  subroutine get_wchi_3x3x3_latlon_field(chi1, chi2, chi3,       &
                                            dlat, dlon, dr, radius, &
                                            dofmap, nlayers)
    !external imports

    use coord_transform_mod, only : llr2xyz

    implicit none

    ! For a field on WChi function space:
    !  - the dimensional scaling factors (da, db, dr) are
    !    the lat-lon-radius increments passed as input (lat, lon) in radians
    !    first dof will be on the equator/meridian
    !  - radius is the surface radius
    !  - the dofmap is passed as input
    !  - the number of layers in the r-direction is passed as input
    !  - a 3-dimensional (chi1,chi2,chi3) coordinate field is returned

    real(r_def), intent(out) :: chi1(:), chi2(:), chi3(:)
    real(r_def), intent(in)  :: dlat, dlon, dr, radius

    integer(i_def), intent(in) :: dofmap(:,:)
    integer(i_def), optional, intent(in) :: nlayers

    real(r_def) :: lat_i, lat_im1, lon_j, lon_jm1, r_k, r_kp1
    real(r_def) :: x, y, z
    integer(i_def) :: cell, i, j, k, nlay, df

    nlay = 2
    if (present(nlayers))then
      nlay=nlayers-1
    end if

    ! Compute coordinates, for P1DG case
    ! Note, since this is a DG1 field the dof locations do not correspond to the
    ! vertex locations and so the reference element indices cannot be used
    cell = 1
    do j = 1,3
      lon_j = real(j, r_def)*dlon
      lon_jm1 = real(j-1, r_def)*dlon
      do i = 1,3
        lat_i = real(i, r_def)*dlat
        lat_im1 = real(i-1, r_def)*dlat
        do k = 0,nlay
          r_k=radius + real(k, r_def)*dr
          r_kp1=radius + real(k+1, r_def)*dr
          df=1
          call llr2xyz(lon_jm1, lat_im1, r_k, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
          df=2
          call llr2xyz(lon_jm1, lat_i, r_k, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
          df=3
          call llr2xyz(lon_j, lat_i, r_k, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
          df=4
          call llr2xyz(lon_j, lat_im1, r_k, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
          df=5
          call llr2xyz(lon_jm1, lat_im1, r_kp1, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
          df=6
          call llr2xyz(lon_jm1, lat_i, r_kp1, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
          df=7
          call llr2xyz(lon_j, lat_i, r_kp1, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
          df=8
          call llr2xyz(lon_j, lat_im1, r_kp1, x, y, z)
          chi1(dofmap(df, cell) + k) = x
          chi2(dofmap(df, cell) + k) = y
          chi3(dofmap(df, cell) + k) = z
        end do
        cell = cell + 1
      end do
    end do

  end subroutine get_wchi_3x3x3_latlon_field

end module get_unit_test_3x3x3_chi_mod
