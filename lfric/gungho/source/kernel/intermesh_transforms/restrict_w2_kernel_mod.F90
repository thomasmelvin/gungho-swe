!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Perform the restriction operation from a fine grid W2 field to a
!!        coarse grid W2 field
!> @details Restrict the W2 fine grid field over a number of cells into a
!!          W2 coarse grid field. The fine grid cells are must be exactly
!!          nested in a coarse grid cell. The coarse field is obtained by
!!          summing contributions from the fine field multiplied by weights.
!!          This method is only designed for the lowest order W2 spaces.
module restrict_w2_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_INC,           &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_SPACE_2, CELL_COLUMN
use constants_mod,           only: i_def, r_def
use fs_continuity_mod,       only: W2
use kernel_mod,              only: kernel_type
use reference_element_mod,   only: W, S, E, N, B

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: restrict_w2_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                          &
       arg_type(GH_FIELD, GH_REAL, GH_INC,   W2,          mesh_arg=GH_COARSE), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  ), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  )  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: restrict_w2_kernel_code
end type restrict_w2_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: restrict_w2_kernel_code

contains

  !> @brief Subroutine to perform the W2 restriction operation
  !> @param[in]     nlayers         Number of layers in a model column
  !> @param[in]     cell_map        Map of which fine grid cells lie in the
  !!                                coarse grid cell
  !> @param[in]     ncell_f_per_c_x Number of fine cells per coarse cell in
  !!                                the x-direction
  !> @param[in]     ncell_f_per_c_y Number of fine cells per coarse cell in
  !!                                the y-direction
  !> @param[in]     ncell_f         Number of cells in the fine grid
  !> @param[in,out] coarse_field    Coarse grid W2 field to compute
  !> @param[in]     fine_field      Fine grid  W2 field to restrict
  !> @param[in]     rmultiplicity   A fine grid W2 field containing weights
  !!                                the reciprocal multiplicity of nodes
  !> @param[in]     undf_c          Total number of unique degrees of freedom
  !!                                for the coarse grid field
  !> @param[in]     dfmap_c         Cell dof-map for the coarse grid field
  !> @param[in]     ndf_f           Number of degrees of freedom per cell for
  !>                                the fine grid field
  !> @param[in]     undf_f          Total number of unique degrees of freedom
  !!                                for the fine grid field
  !> @param[in]     dfmap_f         Cell dof-map for the fine grid field
  subroutine restrict_w2_kernel_code(nlayers,         &
                                     cell_map,        &
                                     ncell_f_per_c_x, &
                                     ncell_f_per_c_y, &
                                     ncell_f,         &
                                     coarse_field,    &
                                     fine_field,      &
                                     rmultiplicity,   &
                                     undf_c,          &
                                     dfmap_c,         &
                                     ndf_f,           &
                                     undf_f,          &
                                     dfmap_f          )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_f_per_c_x
    integer(kind=i_def), intent(in) :: ncell_f_per_c_y
    integer(kind=i_def), intent(in) :: ncell_f
    integer(kind=i_def), intent(in) :: ndf_f
    integer(kind=i_def), intent(in) :: undf_f, undf_c

    ! Fields
    real(kind=r_def),    dimension(undf_c), intent(inout) :: coarse_field
    real(kind=r_def),    dimension(undf_f), intent(in)    :: fine_field
    real(kind=r_def),    dimension(undf_f), intent(in)    :: rmultiplicity

    ! Maps
    integer(kind=i_def), dimension(ndf_f, ncell_f), intent(in) :: dfmap_f
    integer(kind=i_def), dimension(ndf_f),          intent(in) :: dfmap_c
    integer(kind=i_def), dimension(ncell_f_per_c_x, ncell_f_per_c_y), intent(in) :: cell_map

    ! Internal variables
    integer(kind=i_def) :: df, k, lp_x, lp_y
    real(kind=r_def)    :: coarse_N, coarse_S, coarse_E, coarse_W, coarse_B

    !---------------------------------------------------------------------------
    ! Horizontal components
    !---------------------------------------------------------------------------

    do k = 0, nlayers-1

      ! The rows and columns forming the cell map match the arrangment
      ! of fine cells within the coarse cell

      ! These are aligned as follows with the LFRic directions:
      !         N
      !   |--------------|
      !   |    row 1     |
      !   |c            c|
      !   |o            o|
      ! W |l            l| E
      !   |              |
      !   |1           nx|
      !   |    row ny    |
      !   |--------------|
      !          S

      coarse_N = 0.0_r_def
      coarse_S = 0.0_r_def

      do lp_x = 1, ncell_f_per_c_x
        ! N edge is first row of cell map
        lp_y = 1
        df = N
        coarse_N = coarse_N + ( rmultiplicity(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                               * fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) )

        ! S edge is last row of cell map
        lp_y = ncell_f_per_c_y
        df = S
        coarse_S = coarse_S + ( rmultiplicity(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                               * fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) )
      end do

      coarse_field(dfmap_c(N)+k) = coarse_field(dfmap_c(N)+k) + coarse_N
      coarse_field(dfmap_c(S)+k) = coarse_field(dfmap_c(S)+k) + coarse_S

      coarse_E = 0.0_r_def
      coarse_W = 0.0_r_def

      do lp_y = 1, ncell_f_per_c_y
        ! W edge is first column of cell map
        lp_x = 1
        df = W
        coarse_W = coarse_W + ( rmultiplicity(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                               * fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) )

        ! E edge is last column of cell map
        lp_x = ncell_f_per_c_x
        df = E
        coarse_E = coarse_E + ( rmultiplicity(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                               * fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) )
      end do

      coarse_field(dfmap_c(E)+k) = coarse_field(dfmap_c(E)+k) + coarse_E
      coarse_field(dfmap_c(W)+k) = coarse_field(dfmap_c(W)+k) + coarse_W

    end do

    !---------------------------------------------------------------------------
    ! Vertical components
    !---------------------------------------------------------------------------

    ! Only do bottom value of cell
    ! Loop over an extra layer to get the very top
    df = B

    do k = 0, nlayers
      coarse_B = 0.0_r_def
      do lp_y = 1, ncell_f_per_c_y
        do lp_x = 1, ncell_f_per_c_x
          coarse_B = coarse_B + fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k)
        end do
      end do
      coarse_field(dfmap_c(df)+k) = coarse_B
    end do

  end subroutine restrict_w2_kernel_code

end module restrict_w2_kernel_mod
