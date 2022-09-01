!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes RHS of the continuity equation for the nonlinear equations.
!>
!> The kernel computes the RHS of the continuity equation for the nonlinear
!> equations, that is: rrho = -div(F) where F is the mass flux.
!>
module rrho_kernel_mod

  use argument_mod,      only : arg_type, func_type,     &
                                GH_FIELD, GH_REAL,       &
                                GH_READ, GH_WRITE,       &
                                GH_BASIS, GH_DIFF_BASIS, &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W0, W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: rrho_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/             &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2)  &
         /)
    type(func_type) :: meta_funcs(2) = (/           &
         func_type(W3, GH_BASIS),                   &
         func_type(W2, GH_DIFF_BASIS)               &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: rrho_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: rrho_code

contains

!> @brief Compute the right hand side of the continuity equation
!! @param[in] nlayers Number of layers
!! @param[in,out] r_rho Right hand side of the continuity equation
!! @param[in] u Velocity
!! @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!! @param[in] undf_w3 Number of (local) unique degrees of freedom
!! @param[in] map_w3 Dofmap for the cell at the base of the column for W3
!! @param[in] w3_basis Basis functions evaluated at quadrature points
!! @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!! @param[in] undf_w2 Number of (local) unique degrees of freedom
!! @param[in] map_w2 Dofmap for the cell at the base of the column for W2
!! @param[in] w2_diff_basis Differential basis functions evaluated at quadrature points
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine rrho_code(nlayers,                                &
                     r_rho, u,                               &
                     ndf_w3, undf_w3, map_w3, w3_basis,      &
                     ndf_w2, undf_w2, map_w2, w2_diff_basis, &
                     nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3
  integer(kind=i_def), intent(in) :: undf_w2, undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_diff_basis

  real(kind=r_def), dimension(undf_w3), intent(inout) :: r_rho
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, k
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_w2) :: u_e
  real(kind=r_def), dimension(ndf_w3) :: rrho_e
  real(kind=r_def) :: div_u_at_quad

  do k = 0, nlayers-1
    do df = 1, ndf_w3
      rrho_e(df) = 0.0_r_def
    end do
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        div_u_at_quad = 0.0_r_def
        do df = 1, ndf_w2
          div_u_at_quad = div_u_at_quad + u_e(df)*w2_diff_basis(1,df,qp1,qp2)
        end do
        do df = 1, ndf_w3
          rrho_e(df) = rrho_e(df) &
                     - wqp_h(qp1)*wqp_v(qp2)*w3_basis(1,df,qp1,qp2)*div_u_at_quad
        end do
      end do
    end do
    do df = 1, ndf_w3
      r_rho( map_w3(df) + k ) =  rrho_e(df)
    end do
  end do

end subroutine rrho_code

end module rrho_kernel_mod
