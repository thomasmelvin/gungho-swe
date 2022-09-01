!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the sum of values of a W3 field for a column.
!> @details The kernel computes the sum of values in a W3 field (corresponding
!> to masses) for each column. Thus this can be used to calculate the mass of each
!> column, which is returned as a W3 field with the total mass stored in the bottom
!> layer.
!> \f[ \sum_k M_k \f]
!>
module compute_w3_column_total_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_WRITE, GH_READ, &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_w3_column_total_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/             &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: compute_w3_column_total_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_w3_column_total_code

contains

!> @brief Compute the sum of values in a column
!! @param[in] nlayers The number of layers
!! @param[in,out] column_total A field in W3 whose values in the bottom layer
!!                             will be the column totals
!! @param[in] mass A field in W3 space to be summed up
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
subroutine compute_w3_column_total_code(                                    &
                                         nlayers,                           &
                                         column_total,                      &
                                         mass,                              &
                                         ndf_w3, undf_w3, map_w3            &
                                       )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_w3, undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(undf_w3),  intent(inout) :: column_total
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: mass

  ! Internal variables
  integer(kind=i_def) :: df, k

  column_total(:) = 0.0_r_def

  do k = 0, nlayers-1
    do df = 1, ndf_w3
      column_total( map_w3(df) ) = column_total( map_w3(df) ) + mass( map_w3(df) + k )
    end do
  end do

end subroutine compute_w3_column_total_code

end module compute_w3_column_total_kernel_mod
