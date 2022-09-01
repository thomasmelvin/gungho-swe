!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Enforces a lower bound on a field. If any element of the field array
!>        is below this value, then it is set to the lower bound.
!>        This can be used to remove negative values or clip any field to
!>        any desired minimum value.
module enforce_lower_bound_kernel_mod

  use argument_mod,  only : arg_type,            &
                            GH_FIELD, GH_SCALAR, &
                            GH_REAL, GH_READ,    &
                            GH_INC, ANY_SPACE_1, &
                            CELL_COLUMN
  use constants_mod, only : i_def, r_def
  use kernel_mod,    only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  ! The type declaration for the kernel. Contains the metadata needed by the
  ! Psy layer.
  !
  type, public, extends(kernel_type) :: enforce_lower_bound_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                     &
         arg_type(GH_FIELD,  GH_REAL, GH_INC, ANY_SPACE_1), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ )             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: enforce_lower_bound_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: enforce_lower_bound_code

contains

!> @brief Returns field=max(field,lower bound)
!! @param[in] nlayers Number of layers
!! @param[in,out] field Field
!! @param[in] lower_bound The lower bound
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
subroutine enforce_lower_bound_code(nlayers, field, lower_bound, &
                                    ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, undf
  integer(kind=i_def), dimension(ndf), intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: field
  real(kind=r_def), intent(in) :: lower_bound

  ! Internal variables
  integer(kind=i_def) :: df, k

  do k = 0, nlayers-1
    do df = 1, ndf

      ! Clip field
      if (field(map(df)+k) < lower_bound) field(map(df)+k) = lower_bound

    end do
  end do

end subroutine enforce_lower_bound_code

end module enforce_lower_bound_kernel_mod
