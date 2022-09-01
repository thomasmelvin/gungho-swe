!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Injects a Wtheta field into the shifted W3 space.
!> @details Converts a Wtheta field to shifted W3 by identifying the values
!!          at the degrees of freedom as being the same. Some interpolation is
!!          performed in the top and bottom half-cells. This kernel is only
!!          designed to work for the lowest order elements.

module inject_wt_to_sh_w3_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_WRITE,         &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      CELL_COLUMN
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : Wtheta

  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! Psy layer.
  !!
  type, public, extends(kernel_type) :: inject_wt_to_sh_w3_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta)                     &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: inject_wt_to_sh_w3_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: inject_wt_to_sh_w3_code

contains

!> @brief Inject a Wtheta field to a shifted W3 field.
!> @param[in] nlayers_sh  Number of layers in the vertically-shifted mesh.
!> @param[in,out] field_w3_sh Output field in shifted W3.
!> @param[in] field_wt The original field in Wtheta.
!> @param[in] ndf_w3_sh Number of degrees of freedom per cell for W3 shifted
!> @param[in] undf_w3_sh Number of unique degrees of freedom for W3 shifted
!> @param[in] map_w3_sh Dofmap for W3 shifted
!> @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!> @param[in] undf_wt Number of unique degrees of freedom for Wtheta
!> @param[in] map_wt Dofmap for Wtheta
subroutine inject_wt_to_sh_w3_code(                &
                                    nlayers_sh,    &
                                    field_w3_sh,   &
                                    field_wt,      &
                                    ndf_w3_sh,     &
                                    undf_w3_sh,    &
                                    map_w3_sh,     &
                                    ndf_wt,        &
                                    undf_wt,       &
                                    map_wt         &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w3_sh, ndf_wt
  integer(kind=i_def),                           intent(in) :: undf_w3_sh, undf_wt
  integer(kind=i_def), dimension(ndf_w3_sh),     intent(in) :: map_w3_sh
  integer(kind=i_def), dimension(ndf_wt),        intent(in) :: map_wt

  real(kind=r_def),    dimension(undf_w3_sh), intent(inout) :: field_w3_sh
  real(kind=r_def),    dimension(undf_wt),       intent(in) :: field_wt

  ! Internal variables
  integer(kind=i_def) :: k

  ! Assume lowest order so only a single DoF per cell
  ! Bottom boundary value
  field_w3_sh(map_w3_sh(1)) = 0.75_r_def * field_wt(map_wt(1)) &
                              + 0.25_r_def * field_wt(map_wt(1)+1)

  ! Top boundary value
  k = nlayers_sh - 1
  field_w3_sh(map_w3_sh(1)+k) = 0.25_r_def * field_wt(map_wt(1)+k-1) &
                                + 0.75_r_def * field_wt(map_wt(1)+k)

  ! All interior levels
  do k = 1, nlayers_sh - 2
    field_w3_sh(map_w3_sh(1)+k) = field_wt(map_wt(1)+k)
  end do

end subroutine inject_wt_to_sh_w3_code

end module inject_wt_to_sh_w3_kernel_mod
