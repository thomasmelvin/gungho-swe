!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Calculates the weighted average of a data field given some weights

module weighted_ave_kernel_mod

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
  type, public, extends(kernel_type) :: weighted_ave_kernel_type
      private
      type(arg_type) :: meta_args(5) = (/                                        &
           arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
           arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
           arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
           arg_type(GH_SCALAR, GH_INTEGER, GH_READ                            ), &
           arg_type(GH_SCALAR, GH_INTEGER, GH_READ                            )  &
           /)
      integer :: operates_on = CELL_COLUMN
  contains
      procedure, nopass :: weighted_ave_code
  end type

  public :: weighted_ave_code

contains

  !> @param[in]     nlayers      The number of layers
  !> @param[in,out] weighted_ave Output field of weighted sum
  !> @param[in]     data_field   Input field of data to sum
  !> @param[in]     weight_field Input field of weights for data
  !> @param[in]     start_ind    First index of output field to update
  !> @param[in]     num_data     Number of multi-data elements to update
  !> @param[in]     ndf_out      Number of DOFs per cell for output field
  !> @param[in]     undf_out     Number of total DOFs for output field
  !> @param[in]     map_out      Dofmap for cell for output fields
  !> @param[in]     ndf_in       Number of DOFs per cell for input fields
  !> @param[in]     undf_in      Number of total DOFs for input fields
  !> @param[in]     map_in       Dofmap for cell for input fieldss
  subroutine weighted_ave_code(nlayers,                    &
                               weighted_ave,               &
                               data_field,                 &
                               weight_field,               &
                               start_ind,                  &
                               num_data,                   &
                               ndf_out, undf_out, map_out, &
                               ndf_in, undf_in, map_in)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, num_data, start_ind
    integer(kind=i_def), intent(in) :: ndf_in, undf_in
    integer(kind=i_def), intent(in) :: map_in(ndf_in)
    integer(kind=i_def), intent(in) :: ndf_out, undf_out
    integer(kind=i_def), intent(in) :: map_out(ndf_out)

    real(kind=r_def), intent(in)    :: data_field(undf_in)
    real(kind=r_def), intent(in)    :: weight_field(undf_in)
    real(kind=r_def), intent(inout) :: weighted_ave(undf_out)

    integer(kind=i_def) :: i, i_start, i_end
    real(kind=r_def)    :: total_weight

    ! Calculate weighted average over data field
    i_start = start_ind - 1_i_def
    i_end = start_ind + num_data - 2_i_def

    weighted_ave(map_out(1)) = 0.0_r_def
    total_weight = sum(weight_field(map_in(1)+i_start:map_in(1)+i_end))

    do i = i_start, i_end
      if (weight_field(map_in(1)+i) > 0.0_r_def) then
        weighted_ave(map_out(1)) = weighted_ave(map_out(1)) + &
           data_field(map_in(1)+i) * weight_field(map_in(1)+i) / total_weight
      end if
    end do

  end subroutine weighted_ave_code

end module weighted_ave_kernel_mod
