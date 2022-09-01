!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Perform the restriction operation from a fine grid W2 field to a
!!        coarse grid W2 field, using a mask e.g. a limited area mask, for
!!        an intensive field.
!> @details Restrict the W2 fine grid field over a number of cells into a
!!          W2 coarse grid field. The fine grid cells must be exactly
!!          nested in a coarse grid cell. The coarse field is obtained by
!!          summing contributions from the fine field multiplied by weights.
!!          This method is only designed for the lowest order W2 spaces.
!!          Deal with the N,S,E,W dofs separately. For example, for the W dofs
!!          For a coarse mesh cell j and W dof, that contains Nf fine mesh cells
!!          coarse(j,W) = sum_{i=1,Nf} mask(i,W) * fine(i,W) /  sum_{i=1,Nf} mask(i,W)
module restrict_masked_w2_kernel_mod

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
public :: restrict_masked_w2_kernel_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: restrict_masked_w2_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                       &
    arg_type(GH_FIELD, GH_REAL, GH_INC,   W2,          mesh_arg=GH_COARSE), & ! coarse field
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  ), & ! fine field
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  ), & ! rmultiplicity
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  )  & ! fine mesh
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: restrict_masked_w2_kernel_code
end type restrict_masked_w2_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

contains

  !> @brief Subroutine to perform the W2 restriction operation
  !> @param[in]     nlayers         Number of layers in a model column
  !> @param[in]     cell_map        Map of which fine grid cells lie in the
  !!                                coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x Number of fine cells per coarse cell in
  !!                                the x-direction
  !> @param[in]     ncell_fine_per_coarse_y Number of fine cells per coarse cell in
  !!                                the y-direction
  !> @param[in]     ncell_fine      Number of cells in the fine grid
  !> @param[in,out] coarse_field    Coarse grid W2 field to compute
  !> @param[in]     fine_field      Fine grid  W2 field to restrict
  !> @param[in]     rmultiplicity   A fine grid W2 field containing weights
  !!                                the reciprocal multiplicity of nodes
  !> @param[in]     mask_fine       W2 mask on the fine grid
  !> @param[in]     undf_coarse     Total number of unique degrees of freedom
  !!                                for the coarse grid field
  !> @param[in]     dfmap_coarse    Cell dof-map for the coarse grid field
  !> @param[in]     ndf_fine        Number of degrees of freedom per cell for
  !>                                the fine grid field
  !> @param[in]     undf_fine       Total number of unique degrees of freedom
  !!                                for the fine grid field
  !> @param[in]     dfmap_fine      Cell dof-map for the fine grid field
  subroutine restrict_masked_w2_kernel_code(                  &
                                     nlayers,                 &
                                     cell_map,                &
                                     ncell_fine_per_coarse_x, &
                                     ncell_fine_per_coarse_y, &
                                     ncell_fine,              &
                                     coarse_field,            &
                                     fine_field,              &
                                     rmultiplicity,           &
                                     mask_fine,               &
                                     undf_coarse,             &
                                     dfmap_coarse,            &
                                     ndf_fine,                &
                                     undf_fine,               &
                                     dfmap_fine          )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf_fine
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse

    ! Fields
    real(kind=r_def),    dimension(undf_coarse), intent(inout) :: coarse_field
    real(kind=r_def),    dimension(undf_fine), intent(in)    :: fine_field
    real(kind=r_def),    dimension(undf_fine), intent(in)    :: rmultiplicity
    real(kind=r_def),    dimension(undf_fine), intent(in)    :: mask_fine

    ! Maps
    integer(kind=i_def), dimension(ndf_fine, ncell_fine), intent(in) :: dfmap_fine
    integer(kind=i_def), dimension(ndf_fine),             intent(in) :: dfmap_coarse
    integer(kind=i_def), &
      dimension(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y), intent(in) :: cell_map

    ! Internal variables
    integer(kind=i_def) :: df, k, lp_x, lp_y, face
    real(kind=r_def)    :: new_coarse
    real(kind=r_def)    :: non_zero_cells

    integer(kind=i_def), parameter                     :: n_faces = 5
    integer(kind=i_def), dimension(n_faces), parameter :: face_order = [W,S,E,N,B]
    integer(kind=i_def), dimension(n_faces)            :: lp_x_start, lp_x_end, lp_y_start, lp_y_end
    real(kind=r_def),    dimension(n_faces)            :: scaling(n_faces)

    !---------------------------------------------------------------------------
    ! Define cells to average over for each df
    !---------------------------------------------------------------------------

    ! The rows and columns forming the cell map match the arrangment
    ! of fine cells within the coarse cell
    !
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

    do face = 1, size(face_order)
      df = face_order(face)

      select case(df)
      case(N)
        ! N edge is first row of cell map
        lp_x_start(df) = 1
        lp_x_end(df)   = ncell_fine_per_coarse_x
        lp_y_start(df) = 1
        lp_y_end(df)   = 1

      case(S)
        ! S edge is last row of cell map
        lp_x_start(df) = 1
        lp_x_end(df)   = ncell_fine_per_coarse_x
        lp_y_start(df) = ncell_fine_per_coarse_y
        lp_y_end(df)   = ncell_fine_per_coarse_y

      case(W)
        ! W edge is first column of cell map
        lp_x_start(df) = 1
        lp_x_end(df)   = 1
        lp_y_start(df) = 1
        lp_y_end(df)   = ncell_fine_per_coarse_y

      case(E)
        ! E edge is last column of cell map
        lp_x_start(df) = ncell_fine_per_coarse_x
        lp_x_end(df)   = ncell_fine_per_coarse_x
        lp_y_start(df) = 1
        lp_y_end(df)   = ncell_fine_per_coarse_y

      case default
        lp_x_start(df) = 1
        lp_x_end(df)   = ncell_fine_per_coarse_x
        lp_y_start(df) = 1
        lp_y_end(df)   = ncell_fine_per_coarse_y

      end select
    end do

    !---------------------------------------------------------------------------
    ! Calculate scaling based on bottom layer of mask
    !---------------------------------------------------------------------------

    scaling(:) = 0.0_r_def
    do face = 1, size(face_order)
      df = face_order(face)
      non_zero_cells = 0.0_r_def

      do lp_y = lp_y_start(df), lp_y_end(df)
        do lp_x = lp_x_start(df), lp_x_end(df)
          non_zero_cells = non_zero_cells + &
                           mask_fine( dfmap_fine(df,cell_map(lp_x,lp_y)) )
        end do
      end do

      if (non_zero_cells > 0.1_r_def) then
        scaling(df) = 1.0_r_def / non_zero_cells
      else
        scaling(df) = 1.0_r_def
      end if
    end do

    !---------------------------------------------------------------------------
    ! Horizontal components
    !---------------------------------------------------------------------------
    ! Shared dofs, so use rmultiplicity and increment the coarse field.

    do k = 0, nlayers-1

      do face = 1, 4
        df = face_order(face)
        new_coarse = 0.0_r_def

        do lp_y = lp_y_start(df), lp_y_end(df)
          do lp_x = lp_x_start(df), lp_x_end(df)
            new_coarse = new_coarse + ( rmultiplicity(dfmap_fine(df,cell_map(lp_x,lp_y))+k) &
                                        * mask_fine(dfmap_fine(df,cell_map(lp_x,lp_y)))     &
                                        * fine_field(dfmap_fine(df,cell_map(lp_x,lp_y))+k) )
          end do
        end do
        coarse_field(dfmap_coarse(df)+k) = coarse_field(dfmap_coarse(df)+k) + new_coarse * scaling(df)
      end do
    end do

    !---------------------------------------------------------------------------
    ! Vertical components
    !---------------------------------------------------------------------------
    ! Only visit dofs once, so don't use rmultiplicity and set coarse field equal
    ! rather than incrementing.
    !
    ! Only do bottom value of cell
    ! Loop over an extra layer to get the very top
    df = B
    do k = 0, nlayers
      new_coarse = 0.0_r_def

      do lp_y=lp_y_start(df), lp_y_end(df)
        do lp_x=lp_x_start(df), lp_x_end(df)
          new_coarse = new_coarse + ( mask_fine(dfmap_fine(df,cell_map(lp_x,lp_y)))    &
                                  * fine_field(dfmap_fine(df,cell_map(lp_x,lp_y))+k) )
        end do
      end do
      coarse_field(dfmap_coarse(df)+k) = new_coarse * scaling(df)
    end do

  end subroutine restrict_masked_w2_kernel_code

end module restrict_masked_w2_kernel_mod
