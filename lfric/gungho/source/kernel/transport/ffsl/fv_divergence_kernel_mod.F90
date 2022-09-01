!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the finite-volume divergence in one direction.
!> @details The flux form semi-Lagrangian (FFSL) scheme updates density in
!!          the x, y and z directions separately. This code calculates the
!!          divergence for either the x, y, or z direction. The scheme is a simple
!!          finite difference of the fluxes at opposite cell edges and is designed to
!!          work only with lowest order W2 and W3 spaces.

module fv_divergence_kernel_mod

  use argument_mod,       only : arg_type,            &
                                 GH_FIELD, GH_SCALAR, &
                                 GH_REAL, GH_INTEGER, &
                                 GH_WRITE, GH_READ,   &
                                 CELL_COLUMN
  use constants_mod,      only : r_def, i_def
  use flux_direction_mod, only : z_direction
  use fs_continuity_mod,  only : W2, W3
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: fv_divergence_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                 &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ )      &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: fv_divergence_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: fv_divergence_code

contains

!> @brief Computes the finite-volume divergence in either the x, y or z direction.
!> @param[in]     nlayers             The number of layers
!> @param[in,out] mass_divergence     The mass_divergence values in W3 space
!> @param[in]     cell_orientation    The orientation of the cells
!> @param[in]     mass_flux           The flux values which are calculated
!> @param[in]     direction           The direction in which to calculate the fluxes
!> @param[in]     ndf_w3              The number of degrees of freedom per cell
!> @param[in]     undf_w3             The number of unique degrees of freedom
!> @param[in]     map_w3              The dofmap for the cell at the base of the column
!> @param[in]     ndf_w2              The number of degrees of freedom per cell
!> @param[in]     undf_w2             The number of unique degrees of freedom
!> @param[in]     map_w2              The dofmap for the cell at the base of the column
subroutine fv_divergence_code( nlayers,              &
                               mass_divergence,      &
                               cell_orientation,     &
                               mass_flux,            &
                               direction,            &
                               ndf_w3,               &
                               undf_w3,              &
                               map_w3,               &
                               ndf_w2,               &
                               undf_w2,              &
                               map_w2 )

  use cosmic_flux_mod, only : dof_to_update

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                       :: nlayers
  integer(kind=i_def), intent(in)                       :: ndf_w3
  integer(kind=i_def), intent(in)                       :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)    :: map_w3
  real(kind=r_def), dimension(undf_w3), intent(inout)   :: mass_divergence
  real(kind=r_def), dimension(undf_w3), intent(in)      :: cell_orientation
  integer(kind=i_def), intent(in)                       :: ndf_w2
  integer(kind=i_def), intent(in)                       :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)    :: map_w2
  real(kind=r_def), dimension(undf_w2), intent(in)      :: mass_flux
  integer(kind=i_def), intent(in)                       :: direction

  integer(kind=i_def) :: k
  integer(kind=i_def) :: local_dofs(1:2)

  if (direction == z_direction) then
    local_dofs(1) = 5_i_def
    local_dofs(2) = 6_i_def
  else
    local_dofs = dof_to_update(int(cell_orientation(map_w3(1)),i_def),direction)
  end if

  ! This kernel has been designed to work with lowest order W2 and W3 spaces.
  ! As is the case for all of the code associated with the Cosmic transport
  ! scheme.

  do k=0,nlayers-1

    mass_divergence( map_w3(1)+k ) = mass_flux(map_w2(local_dofs(2))+k) -     &
                                            mass_flux(map_w2(local_dofs(1))+k)

  end do

end subroutine fv_divergence_code

end module fv_divergence_kernel_mod
