!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to sample a field at nodal points of another field
module sample_field_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,      &
                                    GH_FIELD, GH_REAL,        &
                                    GH_READ, GH_INC,          &
                                    ANY_SPACE_1, ANY_SPACE_2, &
                                    GH_BASIS, CELL_COLUMN,    &
                                    GH_EVALUATOR
use constants_mod,           only : r_def, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: sample_field_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                     &
       arg_type(GH_FIELD, GH_REAL, GH_INC,  ANY_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_SPACE_2)  &
       /)
  type(func_type) :: meta_funcs(1) = (/                   &
       func_type(ANY_SPACE_2, GH_BASIS)                   &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: sample_field_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: sample_field_code
contains

!> @brief Sample a field at nodal points of another field
!! @param[in]     nlayers       Number of layers
!! @param[in,out] field_1       Field to hold sampled values
!! @param[in]     rmultiplicity Inverse of how many times the dof has been visited in total
!! @param[in]     field_2       Field to take values from
!! @param[in]     ndf_1         Number of degrees of freedom per cell for output field
!! @param[in]     undf_1        Number of unique degrees of freedom for output space
!! @param[in]     map_1         Dofmap for the cell at the base of the column for output space
!! @param[in]     ndf_2         Number of degrees of freedom per cell for the field to be advected
!! @param[in]     undf_2        Number of unique degrees of freedom for the advected field
!! @param[in]     map_2         Dofmap for the cell at the base of the column for the field to be advected
!! @param[in]     basis_2       Basis functions evaluated at Gaussian quadrature points
subroutine sample_field_code(nlayers,                         &
                             field_1, rmultiplicity, field_2, &
                             ndf_1, undf_1, map_1,            &
                             ndf_2, undf_2, map_2, basis_2    &
                            )

  implicit none

  ! Arguments
  integer(kind=i_def),                   intent(in) :: nlayers
  integer(kind=i_def),                   intent(in) :: ndf_1, ndf_2, undf_1, undf_2
  integer(kind=i_def), dimension(ndf_1), intent(in) :: map_1
  integer(kind=i_def), dimension(ndf_2), intent(in) :: map_2

  real(kind=r_def), dimension(1,ndf_2,ndf_1), intent(in)    :: basis_2
  real(kind=r_def), dimension(undf_1),        intent(inout) :: field_1
  real(kind=r_def), dimension(undf_1),        intent(in)    :: rmultiplicity
  real(kind=r_def), dimension(undf_2),        intent(in)    :: field_2

  ! Internal variables
  integer(kind=i_def) :: df, df_2, k, ijk
  real(kind=r_def)    :: f_at_node

  do k = 0, nlayers-1
    do df = 1, ndf_1
      f_at_node = 0.0_r_def
      do df_2 = 1,ndf_2
        f_at_node = f_at_node + field_2(map_2(df_2)+k)*basis_2(1,df_2,df)
      end do
      ijk = map_1(df) + k
      field_1( ijk ) = field_1( ijk ) + f_at_node*rmultiplicity( ijk )
    end do
  end do

end subroutine sample_field_code

end module sample_field_kernel_mod
