!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes weights for the prologation operation on a W2 field.
!> @details Computes weights for performing the prolongation operation to take a
!!          W2 field from a coarse mesh to a fine mesh. The weights are stored
!!          in a W2 field on the fine mesh. For horizontal W2 components, the
!!          fine cell DoFs on the edges of coarse cells have weight given by
!!          the reciprocal of the number of fine cells touching that edge of
!!          the coarse cell. At fine cell DoFs on the interior of the coarse
!!          cell, weights correspond to interpolation between the edges of the
!!          coarse cell. The DoFs for all vertical W2 components of the fine
!!          cells should lie on the edge of the coarse cell, so are treated like
!!          these horizontal components.
!!          This is only designed to work with the lowest order W2 fields and
!!          for fine cells that are exactly nested within a coarse cell.

module weights_prolong_w2_kernel_mod

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

type, public, extends(kernel_type) :: weights_prolong_w2_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                           &
       arg_type(GH_FIELD, GH_REAL, GH_INC,   W2,          mesh_arg=GH_FINE   ), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_COARSE )  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: weights_prolong_w2_kernel_code
end type weights_prolong_w2_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: weights_prolong_w2_kernel_code

contains

  !> @brief Computes weights for the prologation operation on a W2 field.
  !> @param[in]     nlayers         Number of layers in a model column
  !> @param[in]     cell_map        Map of which fine grid cells lie in the
  !!                                coarse grid cell
  !> @param[in]     ncell_f_per_c_x Number of fine cells per coarse cell in
  !!                                the x-direction
  !> @param[in]     ncell_f_per_c_y Number of fine cells per coarse cell in
  !!                                the y-direction
  !> @param[in]     ncell_f         Number of cells in the fine grid
  !> @param[in,out] weights_fine    The fine grid W2 field that contains the
  !!                                weights to be computed here for prolongation
  !> @param[in]     dummy_coarse    An unused coarse W2 field supplied to ensure
  !!                                that appropriate coarse field properties are
  !!                                provided.
  !> @param[in]     undf_f          Total number of unique degrees of freedom
  !!                                for the fine grid field
  !> @param[in]     ndf_f           Number of degrees of freedom per cell
  !!                                for the fine grid field
  !> @param[in]     dfmap_f         Cell dof-map for the fine grid field
  !> @param[in]     undf_c          Total number of unique degrees of freedom
  !!                                for the coarse grid field
  !> @param[in]     dfmap_c         Cell dof-map for the coarse grid field
  subroutine weights_prolong_w2_kernel_code(nlayers,             &
                                            cell_map,            &
                                            ncell_f_per_c_x,     &
                                            ncell_f_per_c_y,     &
                                            ncell_f,             &
                                            weights_fine,        &
                                            dummy_coarse,        &
                                            ndf_f,               &
                                            undf_f,              &
                                            dfmap_f,             &
                                            undf_c,              &
                                            dfmap_c              )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_f_per_c_x
    integer(kind=i_def), intent(in) :: ncell_f_per_c_y
    integer(kind=i_def), intent(in) :: ncell_f
    integer(kind=i_def), intent(in) :: ndf_f
    integer(kind=i_def), intent(in) :: undf_f, undf_c

    ! Fields
    real(kind=r_def),    dimension(undf_f), intent(inout) :: weights_fine
    real(kind=r_def),    dimension(undf_c), intent(in)    :: dummy_coarse

    ! Maps
    integer(kind=i_def), dimension(ndf_f, ncell_f), intent(in) :: dfmap_f
    integer(kind=i_def), dimension(ndf_f),          intent(in) :: dfmap_c
    integer(kind=i_def), dimension(ncell_f_per_c_x, ncell_f_per_c_y), intent(in) :: cell_map

    ! Internal arguments
    integer(kind=i_def) :: df, k, lp_x, lp_y

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
        weights_fine(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                      &
                                        1.0_r_def / real(ncell_f_per_c_x, r_def)

        ! S edge is last row of cell map
        lp_y = ncell_f_per_c_y
        df = S
        weights_fine(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                      &
                                        1.0_r_def / real(ncell_f_per_c_x, r_def)
      end do

      do lp_y = 1, ncell_f_per_c_y
        ! W edge is first column of cell map
        lp_x = 1
        df = W
        weights_fine(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                      &
                                        1.0_r_def / real(ncell_f_per_c_y, r_def)

        ! E edge is last column of cell map
        lp_x = ncell_f_per_c_x
        df = E
        weights_fine(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                      &
                                        1.0_r_def / real(ncell_f_per_c_y, r_def)
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
      do lp_y = 1, ncell_f_per_c_y
        do lp_x = 2, ncell_f_per_c_x

          weights_fine(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                    &
              real(ncell_f_per_c_x + 1 - lp_x, r_def) / real(ncell_f_per_c_x, r_def)
        end do
      end do

      ! Do North/South DoFs
      ! These loops will only execute if there are internal DoFs
      df = N
      do lp_x = 1, ncell_f_per_c_x
        do lp_y = 2, ncell_f_per_c_y

          weights_fine(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                    &
              real(ncell_f_per_c_y + 1 - lp_y, r_def) / real(ncell_f_per_c_y, r_def)

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
          ! Weight is reciprocal of total number of fine cells per coarse
          weights_fine(dfmap_f(df,cell_map(lp_x,lp_y))+k) =                    &
                      1.0_r_def / real(ncell_f_per_c_x * ncell_f_per_c_y, r_def)
        end do
      end do
    end do

  end subroutine weights_prolong_w2_kernel_code

end module weights_prolong_w2_kernel_mod
