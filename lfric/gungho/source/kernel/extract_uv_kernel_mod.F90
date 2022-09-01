!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Extracts the horizontal dof components from a W2 wind field

!> @details Extracts horizontal dof components from a 3D wind field
!>          on W2 and places them in a W2H field

module extract_uv_kernel_mod

use kernel_mod,               only: kernel_type
use argument_mod,             only: arg_type,          &
                                    GH_FIELD, GH_REAL, &
                                    GH_READ, GH_INC,   &
                                    CELL_COLUMN
use constants_mod,            only: r_def, i_def
use fs_continuity_mod,        only: W2, W2H
use kernel_mod,               only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: extract_uv_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/             &
       arg_type(GH_FIELD, GH_REAL, GH_INC,  W2H), &
       arg_type(GH_FIELD, GH_REAL, GH_READ, W2)   &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: extract_uv_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: extract_uv_code
contains

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in,out] h_wind Real array, horizontal components of wind
!! @param[in] u_wind Real array, 3d wind field
!! @param[in] ndf_w2h The number of degrees of freedom per cell for w2h
!! @param[in] undf_w2h The number of unique degrees of freedom for w2h
!! @param[in] map_w2h Integer array holding the dofmap for the cell at the
!>            base of the column for w2h
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number of unique degrees of freedom for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the
!>            base of the column for w2
subroutine extract_uv_code(nlayers,                    &
                           h_wind,                     &
                           u_wind,                     &
                           ndf_w2h, undf_w2h, map_w2h, &
                           ndf_w2, undf_w2, map_w2     &
                           )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2

  real(kind=r_def), dimension(undf_w2h), intent(inout) :: h_wind
  real(kind=r_def), dimension(undf_w2), intent(in)     :: u_wind
  integer(kind=i_def), dimension(ndf_w2h),  intent(in) :: map_w2h
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2

  ! Internal variables
  integer(kind=i_def) :: k, df

  do k = 0, nlayers - 1

    do df=1,4
      h_wind(map_w2h(df) + k) = u_wind(map_w2(df) + k )
    end do

  end do

end subroutine extract_uv_code

end module extract_uv_kernel_mod
