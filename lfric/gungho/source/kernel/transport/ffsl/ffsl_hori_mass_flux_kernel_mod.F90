!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the horizontal fluxes for the FFSL transport scheme.
!> @details This kernel computes the mass fluxes used within the flux form
!!          semi-Lagrangian (FFSL) scheme. It calculates the mass fluxes using
!!          finite volume values rather than using finite-element coefficients
!!          and evaluating the mass fluxes by multiplying by the basis functions.
!!          The mass fluxes are calculated using the PPM + FFSL method
!!          of integrating the density over the departure distance for
!!          each cell face. This kernel is designed to calculate the fluxes in one
!!          direction only and so calculates the fluxes for two of the cell faces only.

module ffsl_hori_mass_flux_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_INC, GH_READ,   &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: ffsl_hori_mass_flux_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/            &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_hori_mass_flux_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: ffsl_hori_mass_flux_code

contains

!> @brief Computes the horizontal fluxes for the FFSL advection scheme.
!> @param[in]     nlayers     Number of layers
!> @param[in]     undf_w3     Number of unique degrees of freedom
!> @param[in]     ndf_w3      Number of degrees of freedom per cell
!> @param[in]     map_w3      Dofmap for the cell at the base of the column
!> @param[in]     rho         Density values in W3
!> @param[in]     a0_coeffs   Coefficients for the subgrid approximation of density
!> @param[in]     a1_coeffs   Coefficients for the subgrid approximation of density
!> @param[in]     a2_coeffs   Coefficients for the subgrid approximation of density
!> @param[in]     undf_w2     Number of unique degrees of freedom
!> @param[in]     ndf_w2      Number of degrees of freedom per cell
!> @param[in]     map_w2      Dofmap for the cell at the base of the column
!> @param[in,out] flux        Flux values which are calculated
!> @param[in]     dep_pts     Departure points
!> @param[in]     stencil_length Length of the 1D stencil
!> @param[in]     stencil_map Dofmaps for the stencil
!> @param[in]     direction   Direction in which to calculate the fluxes
!> @param[in]     deltaT      Length of a timestep
subroutine ffsl_hori_mass_flux_code( nlayers,              &
                                     undf_w3,              &
                                     ndf_w3,               &
                                     map_w3,               &
                                     rho,                  &
                                     a0_coeffs,            &
                                     a1_coeffs,            &
                                     a2_coeffs,            &
                                     undf_w2,              &
                                     ndf_w2,               &
                                     map_w2,               &
                                     flux,                 &
                                     dep_pts,              &
                                     stencil_length,       &
                                     stencil_map,          &
                                     direction,            &
                                     deltaT )

  use cosmic_flux_mod,    only : calc_stencil_ordering,                &
                                 frac_and_int_part,                    &
                                 calc_integration_limits,              &
                                 populate_array,                       &
                                 map_cell_index,                       &
                                 return_part_mass
  use flux_direction_mod, only : x_direction

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                       :: nlayers
  integer(kind=i_def), intent(in)                       :: ndf_w3
  integer(kind=i_def), intent(in)                       :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)    :: map_w3
  real(kind=r_def), dimension(undf_w3), intent(in)      :: rho
  real(kind=r_def), dimension(undf_w3), intent(in)      :: a0_coeffs
  real(kind=r_def), dimension(undf_w3), intent(in)      :: a1_coeffs
  real(kind=r_def), dimension(undf_w3), intent(in)      :: a2_coeffs
  integer(kind=i_def), intent(in)                       :: ndf_w2
  integer(kind=i_def), intent(in)                       :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)    :: map_w2
  real(kind=r_def), dimension(undf_w2), intent(inout)   :: flux
  real(kind=r_def), dimension(undf_w2), intent(in)      :: dep_pts
  integer(kind=i_def), intent(in)                       :: stencil_length
  integer(kind=i_def), intent(in)                       :: stencil_map(1:stencil_length)
  integer(kind=i_def), intent(in)                       :: direction
  real(kind=r_def),    intent(in)                       :: deltaT

  ! Internal variables
  real(kind=r_def) :: mass_total
  real(kind=r_def) :: departure_dist
  real(kind=r_def) :: rho_local(1:stencil_length)
  real(kind=r_def) :: a0_local(1:stencil_length)
  real(kind=r_def) :: a1_local(1:stencil_length)
  real(kind=r_def) :: a2_local(1:stencil_length)
  real(kind=r_def) :: fractional_distance
  real(kind=r_def) :: mass_frac
  real(kind=r_def) :: mass_from_whole_cells
  real(kind=r_def) :: left_integration_limit
  real(kind=r_def) :: right_integration_limit
  real(kind=r_def) :: subgrid_coeffs(3)

  integer(kind=i_def), allocatable :: index_array(:)
  integer(kind=i_def), allocatable :: local_density_index(:)

  integer(kind=i_def) :: stencil_ordering(1:stencil_length)
  integer(kind=i_def) :: k
  integer(kind=i_def) :: dof_iterator
  integer(kind=i_def) :: ii
  integer(kind=i_def) :: edge_options(1:2), local_dofs(1:2)
  integer(kind=i_def) :: n_cells_to_sum


  call calc_stencil_ordering(stencil_length,stencil_ordering)

  edge_options = (/ 0, 1 /)
  if (direction == x_direction) then
    local_dofs = (/ 1, 3 /)
  else
    local_dofs = (/ 2, 4 /)
  end if

  do k=0,nlayers-1

    do dof_iterator=1,2

      departure_dist = dep_pts( map_w2(local_dofs(dof_iterator)) + k )

      ! Rearrange data such that it is in the order 1 | 2 | 3 | 4 | 5 | 6 | 7 etc

      do ii=1,stencil_length
        rho_local(ii) = rho( stencil_map(stencil_ordering(ii)) + k )
        a0_local(ii)  = a0_coeffs( stencil_map(stencil_ordering(ii)) + k )
        a1_local(ii)  = a1_coeffs( stencil_map(stencil_ordering(ii)) + k )
        a2_local(ii)  = a2_coeffs( stencil_map(stencil_ordering(ii)) + k )
      end do

      ! Calculates number of cells of interest and fraction of a cell to add.
      call frac_and_int_part(departure_dist,n_cells_to_sum,fractional_distance)

      ! Calculates the left and right integration limits for the fractional cell.
      call calc_integration_limits( departure_dist,             &
                                    fractional_distance,        &
                                    left_integration_limit,     &
                                    right_integration_limit )

      allocate(index_array(n_cells_to_sum))
      allocate(local_density_index(n_cells_to_sum))

      call populate_array(n_cells_to_sum,index_array,departure_dist,edge_options(dof_iterator))

      do ii=1,n_cells_to_sum
        local_density_index(ii) = map_cell_index(index_array(ii),stencil_length)
      end do

      mass_from_whole_cells = sum(rho_local(local_density_index(1:n_cells_to_sum-1)))

      subgrid_coeffs = (/ a0_local(local_density_index(n_cells_to_sum)), &
                          a1_local(local_density_index(n_cells_to_sum)), &
                          a2_local(local_density_index(n_cells_to_sum)) /)

      mass_frac = return_part_mass(3,subgrid_coeffs,left_integration_limit,right_integration_limit)

      mass_total = mass_from_whole_cells + mass_frac

      flux( map_w2(local_dofs(dof_iterator)) + k ) = &
         sign(1.0_r_def,dep_pts( map_w2(local_dofs(dof_iterator)) + k )) &
         * mass_total / deltaT

      if (allocated(index_array)) deallocate(index_array)
      if (allocated(local_density_index)) deallocate(local_density_index)

    end do

  end do

end subroutine ffsl_hori_mass_flux_code

end module ffsl_hori_mass_flux_kernel_mod
