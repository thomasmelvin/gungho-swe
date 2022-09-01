!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Forms the integrals for projecting from W3 to shifted W3.
!> @details Calculates the integrals in the transform matrix for the RHS of the
!> projection from W3 to shifted W3.
!> As these lie on different vertical extrusions, the integrals must be
!> performed over the double level mesh.
!> The resulting (non-square) matrix is bi-diagonal.
!> We store the components in W3 fields.
!> This is only designed to worked for the lowest-order elements.
module proj_w3_to_sh_w3_rhs_op_kernel_mod

  use argument_mod,          only : arg_type, func_type,       &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_WRITE,         &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    ANY_DISCONTINUOUS_SPACE_9, &
                                    GH_BASIS, GH_DIFF_BASIS,   &
                                    ANY_SPACE_9, CELL_COLUMN,  &
                                    GH_QUADRATURE_XYoZ
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W3
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: proj_w3_to_sh_w3_rhs_op_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_FIELD*2, GH_REAL, GH_WRITE, W3),                        & ! T_ip1, T_i
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! dummy_w3_sh
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               & ! chi_dl
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_9)  & ! panel_id
         /)
    type(func_type) :: meta_funcs(1) = (/                                    &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: proj_w3_to_sh_w3_rhs_op_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: proj_w3_to_sh_w3_rhs_op_code

contains

!> @brief Forms the RHS for projecting a field from W3 to the W3 shifted space.
!>
!> @param[in] nlayers Number of layers in the original mesh
!> @param[in,out] T_ip1 The below diagonal components of the transform matrix. Is a field in W3.
!> @param[in,out] T_i The above diagonal components of the transform matrix. Is a field in W3.
!> @param[in] dummy_w3_sh A dummy field in the shifted W3 space.
!> @param[in] chi_dl_1 The 1st coordinate field in Wchi for double level mesh.
!> @param[in] chi_dl_2 The 2nd coordinate field in Wchi for double level mesh.
!> @param[in] chi_dl_3 The 3rd coordinate field in Wchi for double level mesh.
!> @param[in] panel_id A field giving the ID for the mesh panels.
!> @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!> @param[in] undf_w3 Number of (local) unique degrees of freedom for W3
!> @param[in] map_w3 Dofmap for the cell at the base of the column for W3
!> @param[in] ndf_w3_sh Number of degrees of freedom per cell for W3 shifted
!> @param[in] undf_w3_sh Number of (local) unique degrees of freedom for W3 shifted
!> @param[in] map_w3_sh Dofmap for the cell at the base of the column for W3 shifted
!> @param[in] ndf_wchi_dl Number of degrees of freedom per cell for double level WChi.
!> @param[in] undf_wchi_dl Number of (local) unique degrees of freedom for double level WChi.
!> @param[in] map_wch_dl Dofmap for the cell at the base of the column for double level WChi.
!> @param[in] chi_basis 4-dim array for holding the WChi basis functions
!>                      evaluated at quadrature points.
!> @param[in] chi_diff_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!> @param[in] nqp_h Number of quadrature points in the horizontal
!> @param[in] nqp_v Number of quadrature points in the vertical
!> @param[in] wqp_h Horizontal quadrature weights
!> @param[in] wqp_v Vertical quadrature weights
subroutine proj_w3_to_sh_w3_rhs_op_code( nlayers,        &
                                         T_ip1,          &
                                         T_i,            &
                                         dummy_w3_sh,    &
                                         chi_dl_1,       &
                                         chi_dl_2,       &
                                         chi_dl_3,       &
                                         panel_id,       &
                                         ndf_w3,         &
                                         undf_w3,        &
                                         map_w3,         &
                                         ndf_w3_sh,      &
                                         undf_w3_sh,     &
                                         map_w3_sh,      &
                                         ndf_wchi_dl,    &
                                         undf_wchi_dl,   &
                                         map_wchi_dl,    &
                                         chi_basis,      &
                                         chi_diff_basis, &
                                         ndf_pid,        &
                                         undf_pid,       &
                                         map_pid,        &
                                         nqp_h,          &
                                         nqp_v,          &
                                         wqp_h,          &
                                         wqp_v           &
                                       )

  use coordinate_jacobian_mod,  only: coordinate_jacobian


  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers
  integer(kind=i_def),                           intent(in) :: nqp_h, nqp_v
  integer(kind=i_def),                           intent(in) :: ndf_w3_sh, ndf_w3
  integer(kind=i_def),                           intent(in) :: ndf_wchi_dl, ndf_pid
  integer(kind=i_def),                           intent(in) :: undf_w3_sh, undf_w3
  integer(kind=i_def),                           intent(in) :: undf_wchi_dl, undf_pid
  integer(kind=i_def), dimension(ndf_w3_sh),     intent(in) :: map_w3_sh
  integer(kind=i_def), dimension(ndf_w3),        intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wchi_dl),   intent(in) :: map_wchi_dl
  integer(kind=i_def), dimension(ndf_pid),       intent(in) :: map_pid

  real(kind=r_def),    dimension(undf_w3),    intent(inout) :: T_ip1, T_i
  real(kind=r_def),    dimension(undf_w3_sh),    intent(in) :: dummy_w3_sh
  real(kind=r_def),    dimension(undf_wchi_dl),  intent(in) :: chi_dl_1, chi_dl_2, chi_dl_3
  real(kind=r_def),    dimension(undf_pid),      intent(in) :: panel_id
  real(kind=r_def),    dimension(nqp_h),         intent(in) :: wqp_h
  real(kind=r_def),    dimension(nqp_v),         intent(in) :: wqp_v

  real(kind=r_def), dimension(3, ndf_wchi_dl, nqp_h, nqp_v), intent(in) :: chi_diff_basis
  real(kind=r_def), dimension(3, ndf_wchi_dl, nqp_h, nqp_v), intent(in) :: chi_basis

  ! Internal variables
  integer(kind=i_def)                          :: df, k
  integer(kind=i_def)                          :: qp1, qp2
  real(kind=r_def)                             :: T_ip1_e, T_i_e
  real(kind=r_def), dimension(ndf_wchi_dl)     :: lower_chi_1_e, &
                                                  lower_chi_2_e, &
                                                  lower_chi_3_e, &
                                                  upper_chi_1_e, &
                                                  upper_chi_2_e, &
                                                  upper_chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: lower_dj, upper_dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: lower_jac, upper_jac

  integer(kind=i_def) :: ipanel
  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    ! Extract coordinates for lower and upper half cells
    ! Lower is lower half of orig layer and upper is upper half of orig layer
    do df = 1, ndf_wchi_dl
      lower_chi_1_e(df) = chi_dl_1( map_wchi_dl(df) + 2*k )
      lower_chi_2_e(df) = chi_dl_2( map_wchi_dl(df) + 2*k )
      lower_chi_3_e(df) = chi_dl_3( map_wchi_dl(df) + 2*k )
      upper_chi_1_e(df) = chi_dl_1( map_wchi_dl(df) + 2*k + 1 )
      upper_chi_2_e(df) = chi_dl_2( map_wchi_dl(df) + 2*k + 1 )
      upper_chi_3_e(df) = chi_dl_3( map_wchi_dl(df) + 2*k + 1 )
    end do

    ! Get dj for lower and upper half cells
    call coordinate_jacobian(ndf_wchi_dl, nqp_h, nqp_v, lower_chi_1_e, lower_chi_2_e, &
                             lower_chi_3_e, ipanel, chi_basis, chi_diff_basis,        &
                             lower_jac, lower_dj)
    call coordinate_jacobian(ndf_wchi_dl, nqp_h, nqp_v, upper_chi_1_e, upper_chi_2_e, &
                             upper_chi_3_e, ipanel, chi_basis, chi_diff_basis,        &
                             upper_jac, upper_dj)

    ! Initialise values to zero
    T_ip1_e = 0.0_r_def
    T_i_e = 0.0_r_def

    do qp1 = 1, nqp_h
      do qp2 = 1, nqp_v
        ! Assume lowest order cells, with only one DoF per cell in W3

        ! T_ip1 is volume of upper half cell
        T_ip1_e = T_ip1_e + wqp_h(qp1) * wqp_v(qp2) * lower_dj(qp1,qp2)

        ! T_i is volume of lower half cell
        T_i_e = T_i_e + wqp_h(qp1) * wqp_v(qp2) * upper_dj(qp1,qp2)
      end do
    end do

    T_ip1( map_w3(1) + k) = T_ip1_e
    T_i( map_w3(1) + k) = T_i_e
  end do

end subroutine proj_w3_to_sh_w3_rhs_op_code

end module proj_w3_to_sh_w3_rhs_op_kernel_mod
