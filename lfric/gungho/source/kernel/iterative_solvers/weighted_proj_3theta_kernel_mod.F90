!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Compute the locally assembled projection operator from the
!!        potential temperature space to the W3 space weighted by the inverse
!!        of the reference potential temperature.
!!
!> @details Compute \f[ \left<\sigma,\frac{\gamma}{\theta^{*}}\right>\f] where
!!          sigma is in W3 and gamma is in the potential temperature space.
!!
module weighted_proj_3theta_kernel_mod

  use argument_mod,            only: arg_type, func_type,       &
                                     GH_OPERATOR, GH_FIELD,     &
                                     GH_READ, GH_WRITE,         &
                                     GH_REAL, ANY_SPACE_1,      &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     GH_BASIS, GH_DIFF_BASIS,   &
                                     CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only: r_def, i_def
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use fs_continuity_mod,       only: W3, Wtheta
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: weighted_proj_3theta_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, Wtheta),               &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_1),              &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(3) = (/                                    &
         func_type(W3,          GH_BASIS),                                   &
         func_type(Wtheta,      GH_BASIS),                                   &
         func_type(ANY_SPACE_1, GH_BASIS, GH_DIFF_BASIS)                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: weighted_proj_3theta_code
  end type weighted_proj_3theta_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: weighted_proj_3theta_code

contains

!> @brief Compute the weighted projection operator from Wtheta to W3
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers
!! @param[in] ncell_3d Total number of cells in the 3d mesh
!! @param[in,out] projection Projection operator
!! @param[in] theta Potential temperature array
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!! @param[in] undf_w3 Total number of degrees of freedom for W3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for W3
!! @param[in] basis_w3 Basis functions evaluated at quadrature points
!! @param[in] ndf_wtheta Number of degrees of freedom per cell for the theta space
!! @param[in] undf_wtheta Number of unique degrees of freedom for theta field
!! @param[in] map_wtheta Dofmap for the cell at the base of the column
!! @param[in] basis_wtheta Basis functions evaluated at quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for the coordinate field
!! @param[in] undf_chi Number of unique degrees of freedom for chi field
!! @param[in] map_chi Dofmap for the cell at the base of the column
!! @param[in] basis_chi Wchi basis functions evaluated at Gaussian quadrature points
!! @param[in] diff_basis_chi Wchi derivatives of basis functions
!!                                evaluated at Gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine weighted_proj_3theta_code(cell, nlayers, ncell_3d,             &
                                     projection,                          &
                                     theta,                               &
                                     chi1, chi2, chi3, panel_id,          &
                                     ndf_w3, undf_w3, map_w3,             &
                                     basis_w3,                            &
                                     ndf_wtheta, undf_wtheta, map_wtheta, &
                                     basis_wtheta,                        &
                                     ndf_chi, undf_chi, map_chi,          &
                                     basis_chi, diff_basis_chi,           &
                                     ndf_pid, undf_pid, map_pid,          &
                                     nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell, nlayers, ncell_3d, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_wtheta, ndf_chi, ndf_pid
  integer(kind=i_def), intent(in) :: undf_w3, undf_wtheta, undf_chi, undf_pid

  integer(kind=i_def), dimension(ndf_wtheta),  intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_chi),     intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w3),      intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

  real(kind=r_def), dimension(ndf_w3,ndf_wtheta,ncell_3d),  intent(inout)  :: projection

  real(kind=r_def), dimension(1,ndf_wtheta, nqp_h,nqp_v), intent(in) :: basis_wtheta
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v),     intent(in) :: basis_chi
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v),     intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),      intent(in) :: basis_w3

  real(kind=r_def), dimension(undf_wtheta),  intent(in) :: theta
  real(kind=r_def), dimension(undf_chi),     intent(in) :: chi1
  real(kind=r_def), dimension(undf_chi),     intent(in) :: chi2
  real(kind=r_def), dimension(undf_chi),     intent(in) :: chi3
  real(kind=r_def), dimension(undf_pid),     intent(in) :: panel_id

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  ! Internal variables
  integer(kind=i_def)                          :: df, df1, df2, k, ik, loc, ipanel
  integer(kind=i_def)                          :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: theta_quad
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers - 1

    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi1_e(df) = chi1(loc)
      chi2_e(df) = chi2(loc)
      chi3_e(df) = chi3(loc)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             ipanel, basis_chi, diff_basis_chi, jac, dj)

    ik = k + 1 + (cell-1)*nlayers
    projection(:,:,ik) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        theta_quad = 0.0_r_def
        do df = 1,ndf_wtheta
          theta_quad = theta_quad + theta(map_wtheta(df)+k) &
                                   *basis_wtheta(1,df,qp1,qp2)
        end do
        integrand = wqp_h(qp1)*wqp_v(qp2)/theta_quad*dj(qp1,qp2)

        do df2 = 1, ndf_wtheta
          do df1 = 1, ndf_w3
            projection(df1,df2,ik) = projection(df1,df2,ik)               &
                                   + integrand*basis_w3(1,df1,qp1,qp2)    &
                                              *basis_wtheta(1,df2,qp1,qp2)
          end do
        end do
      end do
    end do
  end do

end subroutine weighted_proj_3theta_code

end module weighted_proj_3theta_kernel_mod
