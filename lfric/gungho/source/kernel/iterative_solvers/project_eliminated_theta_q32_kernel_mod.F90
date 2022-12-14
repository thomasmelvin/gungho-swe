!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the q32 matrix for analytic elimination of theta. The family
!!        of qXY matrices come from elimination of theta from the mixed solver.
!> @details Operator to map the velocity into the left hand side of the
!!          equation of state: q32 = q32 +  const*<sigma, (dtheta/dchi3)/theta * k.v>
!!          where v is a basis function in the W2 space,
!!          sigma is a test function in the W3 space and k is unit vector
!!          in the vertical direction of the reference cell.
!!          For more details, see the solver section of
!!          https://code.metoffice.gov.uk/trac/lfric/wiki/GhaspSupport/Documentation
module project_eliminated_theta_q32_kernel_mod

  use argument_mod,      only: arg_type, func_type,     &
                               GH_OPERATOR, GH_FIELD,   &
                               GH_REAL, GH_SCALAR,      &
                               GH_READ, GH_WRITE,       &
                               GH_BASIS, GH_DIFF_BASIS, &
                               CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only: i_def, r_def
  use fs_continuity_mod, only: W3, W2, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: project_eliminated_theta_q32_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                   &
        arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W2), &
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta), &
        arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3), &
        arg_type(GH_SCALAR,   GH_REAL, GH_READ)           &
        /)
    type(func_type) :: meta_funcs(3) = (/                 &
        func_type(W3,     GH_BASIS),                      &
        func_type(W2,     GH_BASIS),                      &
        func_type(Wtheta, GH_BASIS, GH_DIFF_BASIS)        &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: project_eliminated_theta_q32_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public project_eliminated_theta_q32_code

contains

!> @brief Compute the q32 matrix that arises from analytic elimination of theta
!!        in the equation of state:
!!        q32 = q32 +  const*<sigma, (dtheta/dchi3)/theta * k.v>.
!> @param[in]     cell           Horizontal cell index.
!> @param[in]     nlayers        Number of layers
!> @param[in]     ncell_3d       Number of cells in the 3D mesh
!> @param[in,out] q32_op         Projection matrix
!> @param[in]     theta          Potential temperature field
!> @param[in]     ncell_3d1      Number of cells in the 3D mesh
!> @param[in]     inv_m3         Inverse mass matrix for the W3 space
!> @param[in]     const          Constant scalar to multiply operator by
!> @param[in]     ndf_w3         Degrees of freedom per cell for the pressure space
!> @param[in]     basis_w3       Basis function for the pressure space
!!                               evaluated on quadrature points
!> @param[in]     ndf_w2         Degrees of freedom per cell for the velocity space
!> @param[in]     basis_w2       Vector basis function for the velocity space
!!                               evaluated on quadrature points
!> @param[in]     ndf_wt         Degrees of freedom per cell for the theta space
!> @param[in]     undf_wt        Total degrees of freedom for the theta space
!> @param[in]     map_wt         Cell dofmap for the theta space
!> @param[in]     basis_wt       Basis function for the theta space
!!                               evaluated on quadrature points
!> @param[in]     diff_basis_wt  Differential basis function for the theta space
!!                               evaluated on quadrature points
!> @param[in]     nqp_h          Number of horizontal quadrature points
!> @param[in]     nqp_v          Number of vertical quadrature points
!> @param[in]     wqp_h          Horizontal quadrature weights
!> @param[in]     wqp_v          Vertical quadrature weights
subroutine project_eliminated_theta_q32_code(cell, nlayers, ncell_3d, &
                                             q32_op,                  &
                                             theta,                   &
                                             ncell_3d1, inv_m3,       &
                                             const,                   &
                                             ndf_w3, basis_w3,        &
                                             ndf_w2, basis_w2,        &
                                             ndf_wt, undf_wt, map_wt, &
                                             basis_wt, diff_basis_wt, &
                                             nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ncell_3d, cell, ncell_3d1
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt

  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(3, ndf_w2, nqp_h, nqp_v), intent(in) :: basis_w2
  real(kind=r_def), dimension(1, ndf_w3, nqp_h, nqp_v), intent(in) :: basis_w3
  real(kind=r_def), dimension(1, ndf_wt, nqp_h, nqp_v), intent(in) :: basis_wt
  real(kind=r_def), dimension(3, ndf_wt, nqp_h, nqp_v), intent(in) :: diff_basis_wt

  real(kind=r_def), dimension(ndf_w3, ndf_w2, ncell_3d),  intent(inout) :: q32_op
  real(kind=r_def), dimension(ndf_w3, ndf_w3, ncell_3d1), intent(in)    :: inv_m3

  real(kind=r_def), dimension(undf_wt), intent(in) :: theta
  real(kind=r_def),                     intent(in) :: const
  real(kind=r_def), dimension(nqp_h),   intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v),   intent(in) :: wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, df2, k, ik
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def)                            :: dthetadz_q, theta_q, integrand
  real(kind=r_def), dimension(ndf_w3, ndf_w2) :: proj

  do k = 0, nlayers-1
     ik = 1 + k + (cell-1)*nlayers

    proj(:, :) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h

        theta_q = 0.0_r_def
        dthetadz_q = 0.0_r_def
        do df = 1, ndf_wt
          dthetadz_q = dthetadz_q + theta(map_wt(df)+k)*diff_basis_wt(3, df, qp1, qp2)
          theta_q    = theta_q    + theta(map_wt(df)+k)*basis_wt(1, df, qp1, qp2)
        end do
        ! Ensure that dtheta/dz (and hence the static stability, N^2 =
        ! g/theta*dtheta/dz) is positive
        dthetadz_q = max(1.0_r_def, dthetadz_q)
        do df2 = 1, ndf_w2
          integrand = wqp_h(qp1) * wqp_v(qp2) * const &
                    * dthetadz_q/theta_q * basis_w2(3,df2,qp1,qp2)
          do df = 1, ndf_w3
            proj(df,df2) = proj(df,df2) + basis_w3(1,df,qp1,qp2)*integrand
          end do
        end do
       end do
    end do
    q32_op(:,:,ik) = q32_op(:,:,ik) + matmul(inv_m3(:,:,ik), proj(:,:))
  end do

end subroutine project_eliminated_theta_q32_code

end module project_eliminated_theta_q32_kernel_mod
