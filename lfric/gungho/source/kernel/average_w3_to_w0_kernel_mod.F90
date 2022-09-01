!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Constructs a w0 field from the w3 input field by averaging the
!>        corner values.
!>        If the weight is multiplicity_w0, then this will be a pointwise
!>        mean value around the corners.
!>
!>        Only intended for use with lowest order W3 elements.
!>
module average_w3_to_w0_kernel_mod

  use argument_mod,            only : arg_type,         &
                                      GH_FIELD, GH_INC, &
                                      GH_REAL, GH_READ, &
                                      CELL_COLUMN
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W0, W3

  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: average_w3_to_w0_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/            &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W0), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W0)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: average_w3_to_w0_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: average_w3_to_w0_code

contains

!> @brief
!! @param[in] nlayers Number of layers
!! @param[in,out] field_w0 Surface altitude on lowest order w0
!! @param[in] field_w3 Surface altitude on lowest order w3
!! @param[in] weight_w0 weight used in accumulation of w3 data
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0 Total number of degrees of freedom for w0
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] ndf_w3 Number of degrees of freedom per cell for  for w3
!! @param[in] undf_w3 Total number of degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
subroutine average_w3_to_w0_code(nlayers,                 &
                              field_w0,                   &
                              field_w3,                   &
                              weight_w0,                  &
                              ndf_w0, undf_w0, map_w0,    &
                              ndf_w3, undf_w3, map_w3     &
                              )

  implicit none

  ! Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: undf_w0, ndf_w0
  integer, intent(in) :: undf_w3, ndf_w3

  integer, dimension(ndf_w3), intent(in) :: map_W3
  integer, dimension(ndf_w0), intent(in) :: map_w0

  real(kind=r_def), dimension(undf_w0), intent(inout) :: field_w0
  real(kind=r_def), dimension(undf_w3), intent(in)    :: field_w3
  real(kind=r_def), dimension(undf_w0), intent(in)    :: weight_w0

  ! Internal variables
  integer(kind=i_def) :: df

  ! Simple local average around each cell corner
  do df = 1, ndf_w0
    field_w0(map_w0(df)) = field_w0(map_w0(df)) &
       + field_w3(map_w3(1))/weight_w0(map_w0(df))
  end do

end subroutine average_w3_to_w0_code

end module average_w3_to_w0_kernel_mod
