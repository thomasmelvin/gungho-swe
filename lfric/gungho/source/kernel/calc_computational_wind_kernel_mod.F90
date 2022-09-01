!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the wind used for calculating the departure points in the
!>        dimensionally split advection scheme.
!>
!> Input is the Piola wind and output is the departure wind. The
!> u_departure_wind variable is being used to store the departure winds and we
!> wish to have a positive u_departure_wind representing a postive physical
!> wind. Since u_departure_wind should not be multiplied by the basis
!> functions themselves in order to evaluate them, we remove any dependency on
!> the direction of the basis functions by multiplying the u_piola wind values
!> by the nodal basis functions.
!>
module calc_computational_wind_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_REAL, GH_INC, &
                              GH_READ, ANY_SPACE_9,      &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              GH_DIFF_BASIS, GH_BASIS,   &
                              CELL_COLUMN, GH_EVALUATOR
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2
use kernel_mod,        only : kernel_type

implicit none

private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_computational_wind_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                       &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(2) = (/                                  &
         func_type(W2,          GH_BASIS),                                 &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: calc_computational_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: calc_computational_wind_code

  contains

  !> @param[in] nlayers Number of layers
  !> @param[in,out] u_departure_wind Output field containing the departure wind
  !!                                 used to calculate departure points
  !> @param[in] u_piola Field for the Piola wind
  !> @param[in] chi1 Coordinates in the first direction
  !> @param[in] chi2 Coordinates in the second direction
  !> @param[in] chi3 Coordinates in the third direction
  !> @param[in] panel_id Field giving the ID for mesh panels.
  !> @param[in] ndf Number of degrees of freedom per cell for the output field
  !> @param[in] undf Number of unique degrees of freedom for the output field
  !> @param[in] map Dofmap for the cell at the base of the column for the output field
  !> @param[in] nodal_basis_u Basis functions evaluated at the nodal points for the W2 field
  !> @param[in] ndf_chi Number of degrees of freedom per cell for the coordinate field
  !> @param[in] undf_chi Number of unique degrees of freedom for the coordinate field
  !> @param[in] map_chi Dofmap for the cell at the base of the column for the coordinate field
  !> @param[in] basis_chi Basis functions of the coordinate space evaluated at the nodal points
  !> @param[in] diff_basis_chi Differential basis functions of the coordinate space
  !!                           evaluated at the nodal points
  !> @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
  !> @param[in] undf_pid Number of unique degrees of freedom for panel_id
  !> @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
  subroutine calc_computational_wind_code(nlayers,                                  &
                                          u_departure_wind,                         &
                                          u_piola,                                  &
                                          chi1, chi2, chi3, panel_id,               &
                                          ndf, undf, map, nodal_basis_u,            &
                                          ndf_chi, undf_chi, map_chi,               &
                                          basis_chi, diff_basis_chi,                &
                                          ndf_pid, undf_pid, map_pid                &
                                          )

  use coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                        intent(in)    :: nlayers
  integer(kind=i_def),                        intent(in)    :: ndf, undf
  integer(kind=i_def),                        intent(in)    :: ndf_chi, undf_chi
  integer(kind=i_def),                        intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf),        intent(in)    :: map
  real(kind=r_def), dimension(3,ndf,ndf),     intent(in)    :: nodal_basis_u
  integer(kind=i_def), dimension(ndf_chi),    intent(in)    :: map_chi
  integer(kind=i_def), dimension(ndf_pid),    intent(in)    :: map_pid
  real(kind=r_def), dimension(undf),          intent(in)    :: u_piola
  real(kind=r_def), dimension(undf_chi),      intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid),      intent(in)    :: panel_id
  real(kind=r_def), dimension(undf),          intent(inout) :: u_departure_wind
  real(kind=r_def), dimension(1,ndf_chi,ndf), intent(in)    :: basis_chi
  real(kind=r_def), dimension(3,ndf_chi,ndf), intent(in)    :: diff_basis_chi

  ! Internal variables
  integer(kind=i_def) :: df, k, ipanel
  real(kind=r_def) :: jacobian(3,3,ndf), dj(ndf)
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, ndf, chi1_e, chi2_e, chi3_e,  &
                             ipanel, basis_chi, diff_basis_chi, jacobian, dj)
    do df = 1,ndf
      u_departure_wind(map(df)+k) = u_departure_wind(map(df)+k) + 0.5_r_def* &
          dot_product(nodal_basis_u(:,df,df),abs(nodal_basis_u(:,df,df)))*    &
          u_piola(map(df)+k)/dj(df)
    end do
  end do

  end subroutine calc_computational_wind_code

end module calc_computational_wind_kernel_mod
