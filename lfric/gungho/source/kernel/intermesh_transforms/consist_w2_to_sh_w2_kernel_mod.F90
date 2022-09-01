!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Maps a field from W2 to shifted W2.
!> @details Calculates a shifted W2 field from a W2 field. The new horizontal
!!          and vertical fluxes take the average from each of the half-levels
!!          from the original field. In conjunction with the mapping from W3 to
!!          shifted W3 this provides a consistent mapping of fluxes and
!!          densities to the shifted mesh.
!!          This kernel only works for the lowest-order elements.
module consist_w2_to_sh_w2_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_INC,           &
                                    ANY_SPACE_2, CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : N, E, S, W, T, B

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! Psy layer.
  !!
  type, public, extends(kernel_type) :: consist_w2_to_sh_w2_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  ANY_SPACE_2),               &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2),                        &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2)                         &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: consist_w2_to_sh_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: consist_w2_to_sh_w2_code

contains

!> @brief Maps a field from W2 to the W2 shifted space.
!> @param[in] nlayers_sh Number of layers in the shifted mesh
!> @param[in,out] field_w2_sh Field in the shifted W2 space to be returned.
!> @param[in] field_w2 Original field in W2 to be used.
!> @param[in] rmultiplicity_w2 Reciprocal of nodal multiplicity field for W2
!> @param[in] ndf_w2_sh Number of degrees of freedom per cell for W2 shifted
!> @param[in] undf_w2_sh Number of (local) unique degrees of freedom for W2 shifted
!> @param[in] map_w2_sh Dofmap for the cell at the base of the column for W2 shifted
!> @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!> @param[in] undf_w2 Number of (local) unique degrees of freedom for W2
!> @param[in] map_w2 Dofmap for the cell at the base of the column for W2
subroutine consist_w2_to_sh_w2_code(  nlayers_sh,       &
                                      field_w2_sh,      &
                                      field_w2,         &
                                      rmultiplicity_w2, &
                                      ndf_w2_sh,        &
                                      undf_w2_sh,       &
                                      map_w2_sh,        &
                                      ndf_w2,           &
                                      undf_w2,          &
                                      map_w2            &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w2_sh, ndf_w2
  integer(kind=i_def),                           intent(in) :: undf_w2_sh, undf_w2
  integer(kind=i_def), dimension(ndf_w2_sh),     intent(in) :: map_w2_sh
  integer(kind=i_def), dimension(ndf_w2),        intent(in) :: map_w2

  real(kind=r_def),    dimension(undf_w2_sh), intent(inout) :: field_w2_sh
  real(kind=r_def),    dimension(undf_w2),       intent(in) :: field_w2
  real(kind=r_def),    dimension(undf_w2),       intent(in) :: rmultiplicity_w2

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  integer(kind=i_def) :: horizontal_dofs(4)


  ! We don't want to do operations twice for DoFs that are shared between cells
  ! It would be good to find a way to avoid duplicating this calculation!
  horizontal_dofs = (/ N, E, S, W /)

  ! Loop over layers of original mesh
  do k = 0, nlayers_sh - 2

    ! Loop over horizontal W2 DoFs
    ! Fluxes from the original mesh cells each contribute 1/2 to shifted fluxes
    do j = 1, size(horizontal_dofs)
      df = horizontal_dofs(j)

      field_w2_sh(map_w2_sh(df)+k) = field_w2_sh(map_w2_sh(df)+k)              &
        + rmultiplicity_w2(map_w2(df)+k) * 0.5_r_def * field_w2(map_w2(df)+k)

      field_w2_sh(map_w2_sh(df)+k+1) = field_w2_sh(map_w2_sh(df)+k+1)          &
        + rmultiplicity_w2(map_w2(df)+k) * 0.5_r_def * field_w2(map_w2(df)+k)

    end do
  end do

  do k = 1, nlayers_sh - 1
    ! Loop over vertical W2 DoFs. Only need to do bottom DoF of each cell.
    ! Values are the average from the overlapping cells on the original mesh.
    df = B
    field_w2_sh(map_w2_sh(df)+k) = &
      0.5_r_def * (field_w2(map_w2(df)+k-1) + field_w2(map_w2(df)+k) )
  end do

  ! Top and bottom values are the same as the original space
  field_w2_sh(map_w2_sh(B)) = field_w2(map_w2(B))
  field_w2_sh(map_w2_sh(T)+nlayers_sh-1) = field_w2(map_w2(T)+nlayers_sh-2)

end subroutine consist_w2_to_sh_w2_code

end module consist_w2_to_sh_w2_kernel_mod
