!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------


!> @brief Compute the initial surface geopotential field.
!> @details Computes initial surface geopotential perturbation field from a given
!!          analytic expression.
module initial_surface_geopot_kernel_mod

  use argument_mod,         only: arg_type, func_type,                     &
                                  GH_FIELD, GH_WRITE, GH_READ, GH_INTEGER, &
                                  ANY_SPACE_9, GH_BASIS, GH_REAL,          &
                                  GH_DIFF_BASIS, CELL_COLUMN, GH_EVALUATOR
  use fs_continuity_mod,    only: W3
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type
  use shallow_water_settings_config_mod, &
                            only: swe_test

  implicit none

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: initial_surface_geopot_kernel_type
      private
      type(arg_type) :: meta_args(2) = (/                      &
          arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),         &
          arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9) &
          /)
      type(func_type) :: meta_funcs(1) = (/                    &
          func_type(ANY_SPACE_9, GH_BASIS)                     &
          /)
      integer :: iterates_over = CELL_COLUMN
      integer :: gh_shape = GH_EVALUATOR
  contains
      procedure, nopass :: initial_surface_geopot_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public initial_surface_geopot_code

contains

  !> @brief Compute the initial surface geopotential field.
  !> @param[in]     nlayers      The number of model layers
  !> @param[in,out] geopotential The surface geopotential to compute
  !> @param[in]     chi_1        Real array, the x component of the coordinate field
  !> @param[in]     chi_2        Real array, the y component of the coordinate field
  !> @param[in]     chi_3        Real array, the z component of the coordinate field
  !> @param[in]     ndf_w3       The number of degrees of freedom per cell for w3
  !> @param[in]     undf_w3      The number of total degrees of freedom for w3
  !> @param[in]     map_w3       Integer array holding the dofmap for the cell at the base of the column
  !> @param[in]     ndf_wx       The number of degrees of freedom per cell
  !> @param[in]     undf_wx      The total number of degrees of freedom
  !> @param[in]     map_wx       Integer array holding the dofmap for the cell at the base of the column
  !> @param[in]     wx_basis     Real 5-dim array holding basis functions evaluated at quadrature points
  subroutine initial_surface_geopot_code(nlayers, geopot,               &
                                         chi_1, chi_2, chi_3,           &
                                         ndf_w3, undf_w3, map_w3,       &
                                         ndf_wx, undf_wx, map_wx, wx_basis)

    use analytic_geopot_profiles_mod, only : analytic_surface_geopot

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_w3, ndf_wx, undf_w3, undf_wx
    integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_wx), intent(in) :: map_wx
    real(kind=r_def), dimension(undf_w3),          intent(inout) :: geopot
    real(kind=r_def), dimension(undf_wx),          intent(in)    :: chi_1, chi_2, chi_3
    real(kind=r_def), dimension(1,ndf_wx,ndf_w3),  intent(in)    :: wx_basis

    ! Internal variables
    integer(kind=i_def)                 :: df, df0
    real(kind=r_def), dimension(ndf_wx) :: chi_1_e, chi_2_e, chi_3_e
    real(kind=r_def)                    :: x(3)

    ! compute the pointwise surface geopotential profile
    do df0 = 1, ndf_wx
      chi_1_e(df0) = chi_1( map_wx(df0) )
      chi_2_e(df0) = chi_2( map_wx(df0) )
      chi_3_e(df0) = chi_3( map_wx(df0) )
    end do

    do df = 1, ndf_w3
      x(:) = 0.0_r_def
      do df0 = 1, ndf_wx
        x(1) = x(1) + chi_1_e(df0)*wx_basis(1,df0,df)
        x(2) = x(2) + chi_2_e(df0)*wx_basis(1,df0,df)
        x(3) = x(3) + chi_3_e(df0)*wx_basis(1,df0,df)
      end do

      geopot(map_w3(df)) = analytic_surface_geopot(x, swe_test)

    end do

  end subroutine initial_surface_geopot_code

end module initial_surface_geopot_kernel_mod
