!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Assign field to a value only at a single column
!> @details Sets field to zero everwhere except for the first column, where they are
!>          set to one. The purpose of this kernel is to have a way of setting the field
!>          to zero almost everwhere to test the kernel in
!>          columnwise_op_asm_diag_hmht_kernel_mod.F90
!>          Note that in parallel this might not work, i.e. the values in
!>          several columns will be set to one as a module variable is used to determine
!>          whether the field has already been set.
!>

module assign_field_single_column_kernel_mod

use argument_mod,            only : arg_type,            &
                                    GH_FIELD, GH_REAL,   &
                                    GH_INC, ANY_SPACE_1, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def
use kernel_mod,              only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: assign_field_single_column_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                   &
       arg_type(GH_FIELD, GH_REAL, GH_INC, ANY_SPACE_1) &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: assign_field_single_column_code
end type

! Logical which checks if the non-zero column has already been set.
logical :: set_already = .false.

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: assign_field_single_column_code
contains

!> @brief Sets all field entries to 0, except for the values in the first encountered
!>        column.
!> @param[in]     nlayers Number of vertical layers
!> @param[in,out] x Output data
!> @param[in]     ndf Number of degrees of freedom per cell for the output field
!> @param[in]     undf Unique number of degrees of freedom  for the output field
!> @param[in]     map Dofmap for the cell at the base of the column for the output field
subroutine assign_field_single_column_code(nlayers,     &
                                           x,           &
                                           ndf, undf, map)
  implicit none

  ! Arguments
  integer(kind=i_def),                   intent(in)    :: nlayers
  integer(kind=i_def),                   intent(in)    :: undf, ndf
  real   (kind=r_def), dimension(undf),  intent(inout) :: x
  integer(kind=i_def), dimension(ndf),   intent(in)    :: map

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def)    :: val

  ! Check if one column has already been set. If yet, set value to 1, otherwise
  ! set it to 0.
  if (set_already) then
     val = 0.0_r_def
  else
     val = 1.0_r_def
     set_already = .true.
  end if

  ! Assign values
  do k = 0, nlayers-1
    do df = 1, ndf
       x(map(df)+k) = val
    end do
  end do

end subroutine assign_field_single_column_code

end module assign_field_single_column_kernel_mod
