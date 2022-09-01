!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Place a field into a section of a multi-data field
!> @details Insert data from a smaller (multi-data) field into the required
!>          section of a larger multi-data field. It is assumed that all
!>          data from the input field will be written from (but not all
!>          elements of the output field will be written to).

module multi_insert_kernel_mod

  use argument_mod,  only: arg_type, CELL_COLUMN,     &
                           GH_FIELD, GH_SCALAR,       &
                           GH_REAL, GH_INTEGER,       &
                           GH_READ, GH_WRITE,         &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: multi_insert_kernel_type
      private
      type(arg_type) :: meta_args(4) = (/                                        &
           arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
           arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
           arg_type(GH_SCALAR, GH_INTEGER, GH_READ                            ), &
           arg_type(GH_SCALAR, GH_INTEGER, GH_READ                            )  &
           /)
      integer :: operates_on = CELL_COLUMN
  contains
      procedure, nopass :: multi_insert_code
  end type

  public :: multi_insert_code

contains

  !> @param[in]     nlayers      The number of layers
  !> @param[in,out] output_field Output field to write to
  !> @param[in]     input_field  Input field to write from
  !> @param[in]     start_ind    First index of output field to update
  !> @param[in]     ndata        Number of multi-data elements to update
  !> @param[in]     ndf_out      Number of DOFs per cell for output field
  !> @param[in]     undf_out     Number of total DOFs for output field
  !> @param[in]     map_out      Dofmap for cell for output fields
  !> @param[in]     ndf_in       Number of DOFs per cell for input field
  !> @param[in]     undf_in      Number of total DOFs for input field
  !> @param[in]     map_in       Dofmap for cell for input fields
  subroutine multi_insert_code(nlayers,                    &
                               output_field,               &
                               input_field,                &
                               start_ind,                  &
                               ndata,                      &
                               ndf_out, undf_out, map_out, &
                               ndf_in, undf_in, map_in)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndata, start_ind
    integer(kind=i_def), intent(in) :: ndf_in, undf_in
    integer(kind=i_def), intent(in) :: map_in(ndf_in)
    integer(kind=i_def), intent(in) :: ndf_out, undf_out
    integer(kind=i_def), intent(in) :: map_out(ndf_out)

    real(kind=r_def), intent(in)    :: input_field(undf_in)
    real(kind=r_def), intent(inout) :: output_field(undf_out)

    integer(kind=i_def) :: i

    ! Convert input field to output field
    do i = 0, ndata-1
      output_field(map_out(1)+i+start_ind-1) = input_field(map_in(1)+i)
    end do

  end subroutine multi_insert_code

end module multi_insert_kernel_mod
