!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Perform the prolongation operation from a coarse grid W2 field to a
!!        fine grid W2 field
!> @details Prolong the W2 coarse grid correction into all cells in a W2 field
!!          on a fine grid. The fine grid cells are contained in that coarse
!!          grid cell. Values are found by multiplying by weights. DoFs internal
!!          to the coarse cell are treated differently to those DoFs on the
!!          edges of the coarse cell:
!!          - internal DoFs: values interpolated from coarse cell values
!!          - edge DoFs: divide coarse cell flux by fine cell area fraction
!!          - vertical components: divide coarse cell flux by fine cell area fraction
!!          This method is only designed for the lowest order W2 spaces.

module prolong_w2_kernel_mod

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

type, public, extends(kernel_type) :: prolong_w2_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                           &
       arg_type(GH_FIELD, GH_REAL, GH_INC,   W2,          mesh_arg=GH_FINE ),   &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_COARSE ), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W2,          mesh_arg=GH_FINE )    &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: prolong_w2_kernel_code
end type prolong_w2_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: prolong_w2_kernel_code

contains

  !> @brief Subroutine to perform the W2 prolongation operation
  !> @param[in]     nlayers         Number of layers in a model column
  !> @param[in]     cell_map        Map of which fine grid cells lie in the
  !!                                coarse grid cell
  !> @param[in]     ncell_f_per_c_x Number of fine cells per coarse cell in
  !!                                the x-direction
  !> @param[in]     ncell_f_per_c_y Number of fine cells per coarse cell in
  !!                                the y-direction
  !> @param[in]     ncell_f         Number of cells in the fine grid
  !> @param[in,out] fine_field      Fine grid  W2 field to restrict
  !> @param[in]     coarse_field    Coarse grid W2 field to compute
  !> @param[in]     weights         A fine grid W2 field containing weights
  !!                                for the mapping process
  !> @param[in]     undf_f          Total number of unique degrees of freedom
  !!                                for the fine grid field
  !> @param[in]     ndf_f           Number of degrees of freedom per cell
  !!                                for the fine grid field
  !> @param[in]     dfmap_f         Cell dof-map for the fine grid field
  !> @param[in]     undf_c          Total number of unique degrees of freedom
  !!                                for the coarse grid field
  !> @param[in]     dfmap_c         Cell dof-map for the coarse grid field
  subroutine prolong_w2_kernel_code(nlayers,         &
                                    cell_map,        &
                                    ncell_f_per_c_x, &
                                    ncell_f_per_c_y, &
                                    ncell_f,         &
                                    fine_field,      &
                                    coarse_field,    &
                                    weights,         &
                                    ndf_f,           &
                                    undf_f,          &
                                    dfmap_f,         &
                                    undf_c,          &
                                    dfmap_c          )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_f_per_c_x
    integer(kind=i_def), intent(in) :: ncell_f_per_c_y
    integer(kind=i_def), intent(in) :: ncell_f
    integer(kind=i_def), intent(in) :: ndf_f
    integer(kind=i_def), intent(in) :: undf_f, undf_c

    ! Fields
    real(kind=r_def),    dimension(undf_f), intent(inout) :: fine_field
    real(kind=r_def),    dimension(undf_c), intent(in)    :: coarse_field
    real(kind=r_def),    dimension(undf_f), intent(in)    :: weights

    ! Maps
    integer(kind=i_def), dimension(ndf_f, ncell_f), intent(in) :: dfmap_f
    integer(kind=i_def), dimension(ndf_f),          intent(in) :: dfmap_c
    integer(kind=i_def), dimension(ncell_f_per_c_x, ncell_f_per_c_y), intent(in) :: cell_map

    ! Internal arguments
    integer(kind=i_def) :: df, opp_df, k, lp_x, lp_y
    real(kind=r_def)    :: edge_weight, internal_weight

    !===========================================================================
    ! Horizontal components
    !===========================================================================

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

      !-------------------------------------------------------------------------
      ! Edges of coarse cell
      !-------------------------------------------------------------------------

      ! The weights here describe the flux through the edge of the fine cell
      ! This should be the fractional area of the coarse cell edge taken up
      ! by the fine cell. If fine cells have equal area then this will be the
      ! reciprocal of the number of fine cells along that edge

      do lp_x = 1, ncell_f_per_c_x
        ! N edge is first row of cell map
        lp_y = 1
        df = N
        fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                        &
                                    weights(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                                    * coarse_field(dfmap_c(df)+k)

        ! S edge is last row of cell map
        lp_y = ncell_f_per_c_y
        df = S
        fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                        &
                                    weights(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                                    * coarse_field(dfmap_c(df)+k)
      end do

      do lp_y = 1, ncell_f_per_c_y
        ! W edge is first column of cell map
        lp_x = 1
        df = W
        fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                        &
                                    weights(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                                    * coarse_field(dfmap_c(df)+k)

        ! E edge is last column of cell map
        lp_x = ncell_f_per_c_x
        df = E
        fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                        &
                                    weights(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                                    * coarse_field(dfmap_c(df)+k)
      end do

      !-------------------------------------------------------------------------
      ! DoFs internal to coarse cell
      !-------------------------------------------------------------------------

      ! For internal cells, we interpolate. The weights at these DoFs give the
      ! relative weighting of the WEST and NORTH components. The contributions
      ! from the east and south components are given by 1 minus those weights.
      ! All must be combined with the area fraction weights to give the true
      ! values. As the fine cells are internal to a coarse cell, they will all
      ! have the same orientation as one another.

      ! Do West/East DoFs
      ! These loops will only execute if there are internal DoFs
      df = W
      opp_df = E
      do lp_y = 1, ncell_f_per_c_y
        edge_weight = weights(dfmap_f(df,cell_map(1,lp_y))+k)

        do lp_x = 2, ncell_f_per_c_x
          internal_weight = weights(dfmap_f(df,cell_map(lp_x,lp_y))+k)

          fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                      &
               edge_weight * ( internal_weight * coarse_field(dfmap_c(df)+k) + &
               (1.0_r_def - internal_weight) * coarse_field(dfmap_c(opp_df)+k) )
        end do
      end do

      ! Do North/South DoFs
      ! These loops will only execute if there are internal DoFs
      df = N
      opp_df = S
      do lp_x = 1, ncell_f_per_c_x
        edge_weight = weights(dfmap_f(df,cell_map(lp_x,1))+k)

        do lp_y = 2, ncell_f_per_c_y
          internal_weight = weights(dfmap_f(df,cell_map(lp_x,lp_y))+k)

          fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                      &
               edge_weight * ( internal_weight * coarse_field(dfmap_c(df)+k) + &
               (1.0_r_def - internal_weight) * coarse_field(dfmap_c(opp_df)+k) )

        end do
      end do

    end do

    !===========================================================================
    ! Vertical components
    !===========================================================================

    ! Only do bottom value of cell
    ! Loop over an extra layer to get the very top
    df = B

    do k = 0, nlayers
      do lp_y = 1, ncell_f_per_c_y
        do lp_x = 1, ncell_f_per_c_x
          fine_field(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                      &
                                    weights(dfmap_f(df,cell_map(lp_x,lp_y))+k) &
                                    * coarse_field(dfmap_c(df)+k)
        end do
      end do
    end do

  end subroutine prolong_w2_kernel_code

end module prolong_w2_kernel_mod
