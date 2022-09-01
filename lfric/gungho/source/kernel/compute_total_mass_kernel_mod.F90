!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the cell integrated mass.
!>
module compute_total_mass_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_WRITE, GH_READ, &
                                GH_REAL, ANY_SPACE_9,        &
                                ANY_DISCONTINUOUS_SPACE_3,   &
                                GH_BASIS, GH_DIFF_BASIS,     &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: compute_total_mass_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                       &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),              &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(2) = (/                                   &
         func_type(W3,          GH_BASIS),                                  &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_total_mass_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_total_mass_code

contains

!> @brief Computes the cell integrated mass for a single model column
!! @param[in] nlayers Number of layers
!! @param[in,out] mass Cell integrated mass
!! @param[in] rho Density
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id A field giving the ID for mesh panels
!! @param[in] ndf_w3 Number of degrees of freedom per cell
!! @param[in] undf_w3 Total number of degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Total number of degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Differential basis functions evaluated at
!!                           gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of horizontal quadrature points
!! @param[in] wqp_v Weights of vertical quadrature points
subroutine compute_total_mass_code(                                   &
                                   nlayers,                           &
                                   mass, rho,                         &
                                   chi_1, chi_2, chi_3, panel_id,     &
                                   ndf_w3, undf_w3, map_w3, w3_basis, &
                                   ndf_chi, undf_chi, map_chi,        &
                                   chi_basis, chi_diff_basis,         &
                                   ndf_pid, undf_pid, map_pid,        &
                                   nqp_h, nqp_v, wqp_h, wqp_v         &
                                  )

  use coordinate_jacobian_mod, only : coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                      intent(in)  :: nlayers, nqp_h, nqp_v
  integer(kind=i_def),                      intent(in)  :: ndf_w3, undf_w3
  integer(kind=i_def),                      intent(in)  :: ndf_chi, undf_chi
  integer(kind=i_def),                      intent(in)  :: ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf_w3),   intent(in)  :: map_w3
  integer(kind=i_def), dimension(ndf_chi),  intent(in)  :: map_chi
  integer(kind=i_def), dimension(ndf_pid),  intent(in)  :: map_pid

  real(kind=r_def),    dimension(undf_w3),  intent(inout) :: mass
  real(kind=r_def),    dimension(undf_w3),  intent(in)    :: rho
  real(kind=r_def),    dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def),    dimension(undf_pid), intent(in)    :: panel_id

  real(kind=r_def),    dimension(nqp_h),    intent(in)  :: wqp_h
  real(kind=r_def),    dimension(nqp_v),    intent(in)  :: wqp_v

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),  intent(in) :: w3_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis


  ! Internal variables
  integer(kind=i_def)                          :: df, k, ipanel
  integer(kind=i_def)                          :: qp1, qp2

  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_w3)          :: rho_e, mass_e
  real(kind=r_def)                             :: integrand

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Compute the LHS integrated over one cell and solve
  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_e(df) = chi_1(map_chi(df) + k)
      chi_2_e(df) = chi_2(map_chi(df) + k)
      chi_3_e(df) = chi_3(map_chi(df) + k)
    end do

    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do

    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e, &
                             ipanel, chi_basis, chi_diff_basis, jac, dj)

    do df = 1, ndf_w3
      mass_e(df) = 0.0_r_def

      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          integrand =  w3_basis(1,df,qp1,qp2) * rho_e(df) * dj(qp1,qp2)
          mass_e(df) = mass_e(df) + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do

    mass( map_w3(df) + k ) = mass_e(df)

    end do
  end do

end subroutine compute_total_mass_code

end module compute_total_mass_kernel_mod
