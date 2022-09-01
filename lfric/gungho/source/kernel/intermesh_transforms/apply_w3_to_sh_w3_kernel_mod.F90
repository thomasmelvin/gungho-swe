!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Applies an operator to take a field from W3 to shifted W3.
!> @details Applies a pre-computed operator to map a W3 field to the shifted W3
!!          function space. The operator is stored in two W3 fields, which
!!          represent contributions from the lower and upper halves of cells.
!!          This is only designed to worked for the lowest-order elements.
module apply_w3_to_sh_w3_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_WRITE,         &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      CELL_COLUMN
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W3

  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: apply_w3_to_sh_w3_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                      &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), & ! rhs_w3_sh
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        & ! field_w3
         arg_type(GH_FIELD*2, GH_REAL, GH_READ,  W3)                         & ! T_ip1, T_i
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: apply_w3_to_sh_w3_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: apply_w3_to_sh_w3_code

contains

!> @brief Applies an operator to take a field from W3 to shifted W3.
!> @param[in] nlayers_sh  Number of layers in the vertically-shifted mesh.
!> @param[in,out] rhs_w3_sh Right hand side vector to compute.
!> It is a field in w3 shifted.
!> @param[in] field_w3 The original wind in w3 from which to compute the RHS.
!> @param[in] T_ip1 The below diagonal components of the transform matrix. Is a field in w3.
!> @param[in] T_i The above diagonal components of the transform matrix. Is a field in w3.
!> @param[in] ndf_w3_sh Number of degrees of freedom per cell for w3 shifted
!> @param[in] undf_w3_sh Number of (local) unique degrees of freedom for w3 shifted
!> @param[in] map_w3_sh Dofmap for the cell at the base of the column for w3 shifted
!> @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!> @param[in] undf_w3 Number of (local) unique degrees of freedom for w3
!> @param[in] map_w3 Dofmap for the cell at the base of the column for w3
subroutine apply_w3_to_sh_w3_code(                &
                                   nlayers_sh,    &
                                   rhs_w3_sh,     &
                                   field_w3,      &
                                   T_ip1,         &
                                   T_i,           &
                                   ndf_w3_sh,     &
                                   undf_w3_sh,    &
                                   map_w3_sh,     &
                                   ndf_w3,        &
                                   undf_w3,       &
                                   map_w3         &
                                 )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w3_sh, ndf_w3
  integer(kind=i_def),                           intent(in) :: undf_w3_sh, undf_w3
  integer(kind=i_def), dimension(ndf_w3_sh),     intent(in) :: map_w3_sh
  integer(kind=i_def), dimension(ndf_w3),        intent(in) :: map_w3

  real(kind=r_def),    dimension(undf_w3_sh), intent(inout) :: rhs_w3_sh
  real(kind=r_def),    dimension(undf_w3),       intent(in) :: field_w3
  real(kind=r_def),    dimension(undf_w3),       intent(in) :: T_ip1, T_i

  ! Internal variables
  integer(kind=i_def) :: k

  ! Assume lowest order so only a single DoF per cell
  ! Bottom boundary value
  rhs_w3_sh(map_w3_sh(1)) = T_i(map_w3(1)) * field_w3(map_w3(1))

  ! Top boundary value
  rhs_w3_sh(map_w3_sh(1)+nlayers_sh-1) = &
    T_ip1(map_w3(1)+nlayers_sh-2) * field_w3(map_w3(1)+nlayers_sh-2)

  ! All interior levels
  do k = 1, nlayers_sh - 2
    rhs_w3_sh(map_w3_sh(1)+k) = T_i(map_w3(1)+k) * field_w3(map_w3(1)+k) + &
                                T_ip1(map_w3(1)+k-1) * field_w3(map_w3(1)+k-1)
  end do

end subroutine apply_w3_to_sh_w3_code

end module apply_w3_to_sh_w3_kernel_mod
