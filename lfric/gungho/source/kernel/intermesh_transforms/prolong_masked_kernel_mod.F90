!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Perform the prolongation operation from a coarse grid field to a fine
!!        grid field, using a mask, for an intensive field.
!> @details Prolong the coarse grid correction into a subset of corresponding
!!          sub-cells on the fine grid, as defined by the masking field on the
!!          fine mesh. e.g. using a limited area mask field.
!!          Can only be used with finite elements at lowest-order.
!!          For a coarse mesh cell j, that contains Nf fine mesh cells
!!          fine(i) = mask(i) * coarse(j), for i=1,Nf
module prolong_masked_kernel_mod

use constants_mod, only: i_def, r_def
use kernel_mod,    only: kernel_type
use argument_mod,  only: arg_type,                  &
                         GH_FIELD, GH_REAL,         &
                         GH_READ, GH_READWRITE,     &
                         ANY_DISCONTINUOUS_SPACE_1, &
                         ANY_DISCONTINUOUS_SPACE_2, &
                         GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private
public :: prolong_masked_kernel_code

type, public, extends(kernel_type) :: prolong_masked_kernel_type
   private
   type(arg_type) :: meta_args(3) = (/                                       &
        arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1, & ! field_fine
                                                  mesh_arg=GH_FINE),         &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2, & ! field_coarse
                                                  mesh_arg=GH_COARSE ),      &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1, & ! mask_fine
                                                  mesh_arg=GH_FINE )         &
        /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: prolong_masked_kernel_code
end type prolong_masked_kernel_type

contains

  !> @brief Subroutine to perform the masked prolongation operation
  !> @param[in] nlayers         Number of layers in a model column
  !> @param[in] cell_map        Map of which fine grid cells lie in the coarse
  !!                            grid cell
  !> @param[in] ncell_fine_per_coarse_x Number of fine cells per coarse cell in
  !!                            x-direction
  !> @param[in] ncell_fine_per_coarse_y Number of fine cells per coarse cell in
  !!                            y-direction
  !> @param[in] ncell_fine      Number of cells in the fine grid
  !> @param[in,out] field_fine  Fine grid field to update
  !> @param[in] field_coarse    Coarse grid field to prolong
  !> @param[in] mask_fine       Mask on fine grid
  !> @param[in] ndf             Number of degrees of freedom per cell on both
  !!                            the coarse and fine grid
  !> @param[in] undf_fine       Total number of degrees of freedom on the fine
  !!                            grid
  !> @param[in] dofmap_fine     Cell dofmap on the fine grid
  !> @param[in] undf_coarse     Total number of degrees of freedom on the
  !!                            coarse grid
  !> @param[in] dofmap_coarse   Cell dofmap on the coarse grid
  subroutine prolong_masked_kernel_code( nlayers,                 &
                                         cell_map,                &
                                         ncell_fine_per_coarse_x, &
                                         ncell_fine_per_coarse_y, &
                                         ncell_fine,              &
                                         field_fine,              &
                                         field_coarse,            &
                                         mask_fine,               &
                                         ndf,                     &
                                         undf_fine,               &
                                         dofmap_fine,             &
                                         undf_coarse,             &
                                         dofmap_coarse )

    implicit none

    integer(kind=i_def),                             intent(in) :: nlayers
    integer(kind=i_def),                             intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def),                             intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def),                             intent(in) :: ncell_fine
    integer(kind=i_def),                             intent(in) :: ndf
    integer(kind=i_def), dimension(ndf, ncell_fine), intent(in) :: dofmap_fine
    integer(kind=i_def), dimension(ndf),             intent(in) :: dofmap_coarse
    integer(kind=i_def),                             intent(in) :: undf_fine, undf_coarse
    integer(kind=i_def), &
      dimension(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y), intent(in) :: cell_map

    real(kind=r_def), dimension(undf_coarse),        intent(in)    :: field_coarse
    real(kind=r_def), dimension(undf_fine),          intent(inout) :: field_fine
    real(kind=r_def), dimension(undf_fine),          intent(in)    :: mask_fine

    integer(kind=i_def) :: df, k, lp_x, lp_y

    ! Loop over vertical layers from bottom to top
    ! Wtheta field with ndf=2: update the lower df and loops from k=0 to nlayers
    ! This uses the fact that for Wtheta, map(upper df)+ nlayers-1 = map(lower df) + nlayers
    ! W3 field with ndf=1: update the only (cell centre) df and loop from k=0 to nlayers-1
    ! Use the mask_fine from the bottom-layer (at present mask_fine is stored
    ! as a 3D field - but it may be stored as a 2D field in the future).
    df = 1
    do k = 0, nlayers-1 + (ndf-1)
      do lp_y = 1, ncell_fine_per_coarse_y
        do lp_x = 1, ncell_fine_per_coarse_x
          field_fine( dofmap_fine(df,cell_map(lp_x,lp_y)) + k ) =         &
                  field_fine( dofmap_fine(df,cell_map(lp_x,lp_y)) + k ) + &
                  mask_fine( dofmap_fine(df,cell_map(lp_x,lp_y)) ) *      &
                  field_coarse( dofmap_coarse(df)+k )
        end do
      end do
    end do

  end subroutine prolong_masked_kernel_code

end module prolong_masked_kernel_mod
