!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Assign random values to a field
!> @details Sets all field values to uniformly distributed random numbers
!> in the range [0,scale] for some specified scale.

module assign_field_random_kernel_mod

use argument_mod,            only : arg_type,            &
                                    GH_FIELD, GH_REAL,   &
                                    GH_INC, ANY_SPACE_1, &
                                    CELL_COLUMN, GH_READ, &
                                    GH_SCALAR
use constants_mod,           only : r_def, i_def
use kernel_mod,              only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: assign_field_random_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                  &
    arg_type(GH_FIELD,  GH_REAL, GH_INC, ANY_SPACE_1), &
    arg_type(GH_SCALAR, GH_REAL, GH_READ )             &
  /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: assign_field_random_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: assign_field_random_code
contains

!> @brief Sets all field entries to random values
!> @param[in] nlayers Number of layers
!> @param[in,out] x Output data
!> @param[in] scale Output values are in range [0, scale]
!> @param[in] ndf Number of degrees of freedom per cell for the output field
!> @param[in] undf Unique number of degrees of freedom  for the output field
!> @param[in] map Dofmap for the cell at the base of the column for the output field
subroutine assign_field_random_code(nlayers,       &
                                    x,             &
                                    scale,         &
                                    ndf, undf, map)
  implicit none

  ! Arguments
  integer(kind=i_def),                   intent(in)    :: nlayers
  integer(kind=i_def),                   intent(in)    :: undf, ndf
  integer(kind=i_def), dimension(ndf),   intent(in)    :: map
  real   (kind=r_def), dimension(undf),  intent(inout) :: x
  real   (kind=r_def),                   intent(in)    :: scale

  ! Internal variables
  integer(kind=i_def)              :: df, k
  real(kind=r_def), dimension(ndf) :: random_values

  do k = 0, nlayers-1
    call random_number(random_values(:))
    do df = 1, ndf
      x(map(df)+k) = random_values(df) * scale
    end do
  end do

end subroutine assign_field_random_code

end module assign_field_random_kernel_mod
