!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Injects a shifted W3 field to the Wtheta space.
!> @details Converts a shifted W3 field to Wtheta by identifying the values
!!          at the degrees of freedom as being the same. Some extrapolation is
!!          performed at the top and bottom boundaries. This kernel is only
!!          designed to work for the lowest order elements.

module inject_sh_w3_to_wt_kernel_mod

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
  type, public, extends(kernel_type) :: inject_sh_w3_to_wt_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: inject_sh_w3_to_wt_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: inject_sh_w3_to_wt_code

contains

!> @brief Inject a shifted W3 field to a Wtheta field.
!> @param[in] nlayers  Number of layers in the primary mesh.
!> @param[in,out] field_wt Output field in Wtheta.
!> @param[in] field_w3_sh The original field in shifted W3.
!> @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!> @param[in] undf_wt Number of unique degrees of freedom for Wtheta
!> @param[in] map_wt Dofmap for Wtheta
!> @param[in] ndf_w3_sh Number of degrees of freedom per cell for W3 shifted
!> @param[in] undf_w3_sh Number of unique degrees of freedom for W3 shifted
!> @param[in] map_w3_sh Dofmap for W3 shifted
subroutine inject_sh_w3_to_wt_code(                &
                                    nlayers,       &
                                    field_wt,      &
                                    field_w3_sh,   &
                                    ndf_wt,        &
                                    undf_wt,       &
                                    map_wt,        &
                                    ndf_w3_sh,     &
                                    undf_w3_sh,    &
                                    map_w3_sh      &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers
  integer(kind=i_def),                           intent(in) :: ndf_w3_sh, ndf_wt
  integer(kind=i_def),                           intent(in) :: undf_w3_sh, undf_wt
  integer(kind=i_def), dimension(ndf_w3_sh),     intent(in) :: map_w3_sh
  integer(kind=i_def), dimension(ndf_wt),        intent(in) :: map_wt

  real(kind=r_def),    dimension(undf_w3_sh),    intent(in) :: field_w3_sh
  real(kind=r_def),    dimension(undf_wt),    intent(inout) :: field_wt

  ! Internal variables
  integer(kind=i_def) :: k

  ! Assume lowest order so only a single DoF per cell
  ! Bottom boundary value
  field_wt(map_wt(1)) = (4.0_r_def / 3.0_r_def) * field_w3_sh(map_w3_sh(1)) &
                        - (1.0_r_def / 3.0_r_def) * field_w3_sh(map_w3_sh(1)+1)

  ! Top boundary value
  k = nlayers
  field_wt(map_wt(1)+k) = - (1.0_r_def / 3.0_r_def) * field_w3_sh(map_w3_sh(1)+k-1) &
                          + (4.0_r_def / 3.0_r_def) * field_w3_sh(map_w3_sh(1)+k)

  ! All interior levels
  do k = 1, nlayers - 1
    field_wt(map_wt(1)+k) = field_w3_sh(map_w3_sh(1)+k)
  end do

end subroutine inject_sh_w3_to_wt_code

end module inject_sh_w3_to_wt_kernel_mod
