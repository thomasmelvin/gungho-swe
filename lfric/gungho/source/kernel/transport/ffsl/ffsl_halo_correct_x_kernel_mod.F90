!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Kernel which performs halo exchange dependent on halo cell orientation
!!        when the exchange is in the x direction relative to the cubed-sphere
!!        panel's orientation.

module ffsl_halo_correct_x_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_WRITE, GH_READ, &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: ffsl_halo_correct_x_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/             &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_halo_correct_x_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: ffsl_halo_correct_x_code

contains

!> @brief Kernel which performs halo exchange in the x direction.
!> @param[in] nlayers           The number of layers
!> @param[in,out] rho_out       Density field after halo exchange
!> @param[in] rho_x_in          Density field created by the x-direction FFSL update
!> @param[in] rho_y_in          Density field created by the y-direction FFSL update
!> @param[in] cell_orientation  Orientation of cells, in particular halo cells
!> @param[in] ndf_w3            The number of degrees of freedom per cell
!> @param[in] undf_w3           The number of unique degrees of freedom
!> @param[in] map_w3            Array holding the W3 dofmap
subroutine ffsl_halo_correct_x_code( nlayers,              &
                                     rho_out,              &
                                     rho_x_in,             &
                                     rho_y_in,             &
                                     cell_orientation,     &
                                     ndf_w3,               &
                                     undf_w3,              &
                                     map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                       :: nlayers
  integer(kind=i_def), intent(in)                       :: ndf_w3
  integer(kind=i_def), intent(in)                       :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)    :: map_w3
  real(kind=r_def), dimension(undf_w3), intent(in)      :: cell_orientation
  real(kind=r_def), dimension(undf_w3), intent(in)      :: rho_x_in
  real(kind=r_def), dimension(undf_w3), intent(in)      :: rho_y_in
  real(kind=r_def), dimension(undf_w3), intent(inout)   :: rho_out

  integer(kind=i_def) :: k
  integer(kind=i_def) :: orientation

  orientation = int(cell_orientation(map_w3(1)))

  do k=0,nlayers-1

    if (orientation == 1 .OR. orientation == 3 ) then
      rho_out(map_w3(1)+k) = rho_x_in(map_w3(1)+k)
    else if (orientation == 2 .OR. orientation == 4 ) then
      rho_out(map_w3(1)+k) = rho_y_in(map_w3(1)+k)
    end if

  end do

end subroutine ffsl_halo_correct_x_code

end module ffsl_halo_correct_x_kernel_mod
