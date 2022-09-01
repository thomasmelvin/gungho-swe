!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Samples a field from W2 to shifted W2.
!> @details Calculates a shifted W2 field from a W2 by sampling the velocity
!!          and using area weightings. This preserves a horizontal velocity that
!!          varies linearly in the vertical.
!!          This kernel only works with the lowest-order finite element spaces
!!          on quadrilateral cells.
module sample_w2_to_sh_w2_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_INC,           &
                                    ANY_SPACE_2, ANY_SPACE_5,  &
                                    ANY_DISCONTINUOUS_SPACE_7, &
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2, Wtheta, W3
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : N, E, S, W, T, B

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: sample_w2_to_sh_w2_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  ANY_SPACE_2),               & ! field_w2_sh
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                        & ! field_w2
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_SPACE_2),               & ! area_w2_sh
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                        & ! area_w2
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_SPACE_5),               & ! area_w2_dl
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_7), & ! height_wt_sh
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta)                     & ! height_wt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: sample_w2_to_sh_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: sample_w2_to_sh_w2_code

contains

!> @brief Samples a W2 field in the W2 shifted space.
!>
!> @param[in] nlayers_sh Number of layers in the shifted mesh
!> @param[in,out] field_w2_sh Field in the shifted W2 space to be returned.
!> @param[in] field_w2 Original field in W2 to be used.
!> @param[in] area_w2_sh The areas of cell faces of the shifted mesh. In W2 shifted.
!> @param[in] area_w2 The areas of cell faces of the primal mesh. In W2.
!> @param[in] area_w2_dl The areas of cell faces on the double mesh. In W2 on the double mesh.
!> @param[in] height_wt_sh Heights of the DoFs in the wt shifted space.
!> @param[in] height_wt The heights of the wt DoFs.
!> @param[in] ndf_w2_sh Number of degrees of freedom per cell for W2 shifted
!> @param[in] undf_w2_sh Number of (local) unique degrees of freedom for W2 shifted
!> @param[in] map_w2_sh Dofmap for the cell at the base of the column for W2 shifted
!> @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!> @param[in] undf_w2 Number of (local) unique degrees of freedom for W2
!> @param[in] map_w2 Dofmap for the cell at the base of the column for W2
!> @param[in] ndf_wt_sh Number of degrees of freedom per cell for Wtheta shifted
!> @param[in] undf_wt_sh Number of (local) unique degrees of freedom for Wtheta shifted
!> @param[in] map_wt_sh Dofmap for the cell at the base of the column for Wtheta shifted
!> @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!> @param[in] undf_w3 Number of (local) unique degrees of freedom for W3
!> @param[in] map_w3 Dofmap for the cell at the base of the column for W3
!> @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!> @param[in] undf_wt Number of (local) unique degrees of freedom for Wtheta
!> @param[in] map_wt Dofmap for the cell at the base of the column for Wtheta
subroutine sample_w2_to_sh_w2_code( nlayers_sh,       &
                                    field_w2_sh,      &
                                    field_w2,         &
                                    area_w2_sh,       &
                                    area_w2,          &
                                    area_w2_dl,       &
                                    height_wt_sh,     &
                                    height_wt,        &
                                    ndf_w2_sh,        &
                                    undf_w2_sh,       &
                                    map_w2_sh,        &
                                    ndf_w2,           &
                                    undf_w2,          &
                                    map_w2,           &
                                    ndf_w2_dl,        &
                                    undf_w2_dl,       &
                                    map_w2_dl,        &
                                    ndf_wt_sh,        &
                                    undf_wt_sh,       &
                                    map_wt_sh,        &
                                    ndf_wt,           &
                                    undf_wt,          &
                                    map_wt            &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w2_sh, ndf_w2
  integer(kind=i_def),                           intent(in) :: ndf_w2_dl
  integer(kind=i_def),                           intent(in) :: ndf_wt_sh, ndf_wt
  integer(kind=i_def),                           intent(in) :: undf_w2_sh, undf_w2
  integer(kind=i_def),                           intent(in) :: undf_w2_dl
  integer(kind=i_def),                           intent(in) :: undf_wt_sh, undf_wt
  integer(kind=i_def), dimension(ndf_w2_sh),     intent(in) :: map_w2_sh
  integer(kind=i_def), dimension(ndf_w2),        intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w2_dl),     intent(in) :: map_w2_dl
  integer(kind=i_def), dimension(ndf_wt_sh),     intent(in) :: map_wt_sh
  integer(kind=i_def), dimension(ndf_wt),        intent(in) :: map_wt

  real(kind=r_def),    dimension(undf_w2_sh), intent(inout) :: field_w2_sh
  real(kind=r_def),    dimension(undf_w2),       intent(in) :: field_w2
  real(kind=r_def),    dimension(undf_w2_sh),    intent(in) :: area_w2_sh
  real(kind=r_def),    dimension(undf_w2),       intent(in) :: area_w2
  real(kind=r_def),    dimension(undf_w2_dl),    intent(in) :: area_w2_dl
  real(kind=r_def),    dimension(undf_wt_sh),    intent(in) :: height_wt_sh
  real(kind=r_def),    dimension(undf_wt),       intent(in) :: height_wt

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  integer(kind=i_def) :: horizontal_dofs(4)


  ! We don't want to do operations twice for DoFs that are shared between cells
  ! It would be good to find a way to avoid duplicating this calculation!
  horizontal_dofs = (/ N, E, S, W /)

  ! Loop over layers (but not top and bottom)
  do k = 1, nlayers_sh - 2

    ! Loop over horizontal W2 DoFs
    ! Values are weighted by the shared area between the meshes
    do j = 1, size(horizontal_dofs)
      df = horizontal_dofs(j)
      field_w2_sh(map_w2_sh(df)+k) =                                  &
        (area_w2_dl(map_w2_dl(df)+2*k) / area_w2(map_w2(df)+k))       &
        * field_w2(map_w2(df)+k)                                      &
        + (area_w2_dl(map_w2_dl(df)+2*k-1) / area_w2(map_w2(df)+k-1)) &
        * field_w2(map_w2(df)+k-1)
    end do
  end do
  do k = 1, nlayers_sh - 1
    ! Loop over vertical W2 DoFs. Only need to do bottom DoF of each cell.
    ! Values are obtained from linear interpolation.
    df = B
    field_w2_sh(map_w2_sh(df)+k) = area_w2_sh(map_w2_sh(df)+k) * &
      ( field_w2(map_w2(df)+k-1) / area_w2(map_w2(df)+k-1) + &
        (field_w2(map_w2(df+1)+k-1) / area_w2(map_w2(df+1)+k-1) - &
          field_w2(map_w2(df)+k-1) / area_w2(map_w2(df)+k-1)) * &
        (height_wt_sh(map_wt_sh(1)+k) - height_wt(map_wt(1)+k-1)) / &
        (height_wt(map_wt(2)+k-1) - height_wt(map_wt(1)+k-1)))
  end do

  ! Top and bottom values are the same as the original space
  field_w2_sh(map_w2_sh(B)) = field_w2(map_w2(B))
  field_w2_sh(map_w2_sh(T)+nlayers_sh-1) = field_w2(map_w2(T)+nlayers_sh-2)

  ! Finally, adjust horizontal DoFs in the top and bottom layers
  do j = 1, size(horizontal_dofs)
    df = horizontal_dofs(j)
    ! Bottom layer
    field_w2_sh(map_w2_sh(df)) =                                  &
      (area_w2_dl(map_w2_dl(df)) / area_w2(map_w2(df)))     &
      * field_w2(map_w2(df))

    ! Top layer
    k = nlayers_sh-1
    field_w2_sh(map_w2_sh(df)+k) =                                  &
      (area_w2_dl(map_w2_dl(df)+2*k-1) / area_w2(map_w2(df)+k-1)) &
      * field_w2(map_w2(df)+k-1)
  end do

end subroutine sample_w2_to_sh_w2_code

end module sample_w2_to_sh_w2_kernel_mod
