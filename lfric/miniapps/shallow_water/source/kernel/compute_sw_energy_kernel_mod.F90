!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes cell integrated energy and potential enstrophy.
!!
!> @details The kernel computes the cell integrated energy,
!!          \f[ E = \frac{1}{2} h ( (u^2+v^2) + gh + 2 g h_s), \f]
!!          potential enstrophy
!!          \f[ Z = 1/2 phi q^2, \f]
!!          and potential vorticity q.
!!
module compute_sw_energy_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_WRITE, GH_READ, &
                                GH_REAL, ANY_SPACE_9,        &
                                ANY_DISCONTINUOUS_SPACE_3,   &
                                GH_BASIS, GH_DIFF_BASIS,     &
                                GH_SCALAR,                   &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

!---------------------------------------------------------------------------
! Public types
!---------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>
type, public, extends(kernel_type) :: compute_sw_energy_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                      &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        & ! energy
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        & ! enstrophy
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        & ! pv
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        & ! geopot
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        & ! q
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W2),                        & ! wind
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        & ! suface geopot
       arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               & ! chi
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  & ! panel ID
       /)
  type(func_type) :: meta_funcs(3) = (/                                   &
       func_type(W2,          GH_BASIS),                                  &
       func_type(W3,          GH_BASIS),                                  &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                    &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_sw_energy_code
end type

!---------------------------------------------------------------------------
! Contained functions/subroutines
!---------------------------------------------------------------------------
public :: compute_sw_energy_code
contains

!> @brief Compute the cell integrated total energy
!! @param[in] nlayers The number of layers
!! @param[in,out] energy    The cell integrated energy
!! @param[in,out] enstrophy The cell integrated potential enstrophy
!! @param[in,out] pv        The cell integrated potential vorticity
!! @param[in] geopot The geopotential
!! @param[in] q The potential vorticity
!! @param[in] u The velocity array
!! @param[in] s_geopot The surface geopotential
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id A field giving the ID for mesh panels
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis 4-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number of unique degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis 4-dim array holding basis functions evaluated at quadrature points
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Wchi basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Differential Wchi basis functions evaluated at gaussian quadrature point
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_sw_energy_code( nlayers,                           &
                                   energy, enstrophy, pv,             &
                                   geopot, q, u, s_geopot,            &
                                   chi_1, chi_2, chi_3, panel_id,     &
                                   ndf_w3, undf_w3, map_w3, w3_basis, &
                                   ndf_w2, undf_w2, map_w2, w2_basis, &
                                   ndf_chi, undf_chi, map_chi,        &
                                   chi_basis, chi_diff_basis,         &
                                   ndf_pid, undf_pid, map_pid,        &
                                   nqp_h, nqp_v, wqp_h, wqp_v         &
                                   )

  use coordinate_jacobian_mod,  only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                        intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def),                        intent(in) :: ndf_w2, ndf_w3
  integer(kind=i_def),                        intent(in) :: ndf_pid, ndf_chi
  integer(kind=i_def),                        intent(in) :: undf_w2, undf_w3
  integer(kind=i_def),                        intent(in) :: undf_pid, undf_chi
  integer(kind=i_def), dimension(ndf_w2),     intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3),     intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_chi),    intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),    intent(in) :: map_pid

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),     intent(in) :: w3_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),     intent(in) :: w2_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v),    intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v),    intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w3),     intent(inout) :: energy
  real(kind=r_def), dimension(undf_w3),     intent(inout) :: enstrophy
  real(kind=r_def), dimension(undf_w3),     intent(inout) :: pv
  real(kind=r_def), dimension(undf_w2),     intent(in)    :: u
  real(kind=r_def), dimension(undf_w3),     intent(in)    :: q, geopot, s_geopot
  real(kind=r_def), dimension(undf_chi),    intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid),    intent(in)    :: panel_id

  real(kind=r_def), dimension(nqp_h),       intent(in)  ::  wqp_h
  real(kind=r_def), dimension(nqp_v),       intent(in)  ::  wqp_v

  ! Internal variables
  integer(kind=i_def)   :: df, k, loc
  integer(kind=i_def)   :: qp1, qp2
  integer(kind=i_def)   :: ipanel

  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w3)          :: energy_e, enstrophy_e, pv_e
  real(kind=r_def), dimension(ndf_w2)          :: u_e
  real(kind=r_def), dimension(ndf_w3)          :: geopot_e, s_geopot_e
  real(kind=r_def), dimension(ndf_w3)          :: h_e, s_h_e, q_e

  real(kind=r_def) :: u_at_quad(3), q_at_quad
  real(kind=r_def) :: geopot_at_quad, s_geopot_at_quad, ke_term, pe_term
  real(kind=r_def) :: h_at_quad, s_h_at_quad

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             ipanel, chi_basis, chi_diff_basis, jac, dj)

    do df = 1, ndf_w3
      geopot_e(df)   = geopot( map_w3(df) + k )
      s_geopot_e(df) = s_geopot( map_w3(df) + k )
      h_e(df)   = geopot( map_w3(df) + k ) / 9.80616_r_def
      s_h_e(df) = s_geopot( map_w3(df) + k ) / 9.80616_r_def
      q_e(df)   = q( map_w3(df) + k )
      energy_e(df) = 0.0_r_def
      enstrophy_e(df) = 0.0_r_def
    end do
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
    end do
    ! Compute the energy integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        geopot_at_quad = 0.0_r_def
        s_geopot_at_quad = 0.0_r_def
        h_at_quad = 0.0_r_def
        s_h_at_quad = 0.0_r_def
        q_at_quad = 0.0_r_def
        do df = 1, ndf_w3
          geopot_at_quad   = geopot_at_quad   + geopot_e(df)  *w3_basis(1,df,qp1,qp2)
          s_geopot_at_quad = s_geopot_at_quad + s_geopot_e(df)*w3_basis(1,df,qp1,qp2)
          h_at_quad   = h_at_quad   + h_e(df)  *w3_basis(1,df,qp1,qp2)
          s_h_at_quad = s_h_at_quad + s_h_e(df)*w3_basis(1,df,qp1,qp2)
          q_at_quad  = q_at_quad + q_e(df)*w3_basis(1,df,qp1,qp2)
        end do
        ! PE term
        pe_term = geopot_at_quad + 2.0_r_def*s_geopot_at_quad
        ! k.e term
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
                       + u_e(df)*w2_basis(:,df,qp1,qp2)
        end do
        ke_term = dot_product(matmul(jac(:,:,qp1,qp2),u_at_quad), &
                                        matmul(jac(:,:,qp1,qp2),u_at_quad))/(dj(qp1,qp2)**2.0_r_def)
        do df = 1, ndf_w3
          energy_e(df) = energy_e(df) + 0.5_r_def*wqp_h(qp1)*wqp_v(qp2)*h_at_quad &
                  *(ke_term + pe_term)*dj(qp1,qp2)

          enstrophy_e(df) = enstrophy_e(df) + 0.5_r_def*wqp_h(qp1)*wqp_v(qp2)*q_at_quad &
                  *q_at_quad*geopot_at_quad*dj(qp1,qp2)

          pv_e(df) = pv_e(df) + wqp_h(qp1)*wqp_v(qp2)*q_at_quad
        end do
      end do
    end do
    do df = 1, ndf_w3
      energy( map_w3(df) + k ) = energy_e(df)
      enstrophy( map_w3(df) + k ) = enstrophy_e(df)
      pv( map_w3(df) + k ) = pv_e(df)
    end do
  end do

end subroutine compute_sw_energy_code

end module compute_sw_energy_kernel_mod
