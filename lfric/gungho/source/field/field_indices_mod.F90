!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief Define indices for the prognostic field vectors
!>
!> @details Define and set indices for the prognostic field vectors
module field_indices_mod

  implicit none

  ! For gungho model
  integer, parameter :: igh_u = 1  ! wind
  integer, parameter :: igh_t = 2  ! potential temperature
  integer, parameter :: igh_d = 3  ! density
  integer, parameter :: igh_p = 4  ! Exner pressure
  integer, parameter :: igh_uv = igh_u ! uv wind if split - must be the same as igh_u
  integer, parameter :: igh_w = 5  ! w wind if split

  ! For gravity wave mini app
  integer, parameter :: igw_u = 1  ! wind
  integer, parameter :: igw_p = 2  ! pressure
  integer, parameter :: igw_b = 3  ! buoyancy
  integer, parameter :: igw_uv = igw_u ! uv wind if split - must be the same as igw_u
  integer, parameter :: igw_w = 4  ! w wind if split

  ! For shallow water mini app
  integer, parameter :: isw_u = 1  ! wind
  integer, parameter :: isw_g = 2  ! geopotential
  integer, parameter :: isw_b = 3  ! buoyancy
  integer, parameter :: isw_q = 4  ! vorticity

  ! For Semi-Implicit Solver
  integer, parameter :: isol_u = 1       ! wind
  integer, parameter :: isol_uv = isol_u ! uv wind if split - must be the same as isol_u
  integer, protected :: isol_p           ! Exner pressure
  integer, protected :: isol_d           ! density
  integer, protected :: isol_t           ! potential temperature
  integer, protected :: isol_w           ! w wind if split

contains

  !> @brief Set the solver indices to be a contiguous block of indices.
  subroutine set_solver_field_indices()

    use constants_mod,           only: IMDI
    use mixed_solver_config_mod, only: eliminate_variables, &
                                       eliminate_variables_none

    implicit none

    if ( eliminate_variables == eliminate_variables_none ) then
      ! Set isol_ to match igh_
      isol_t = igh_t
      isol_d = igh_d
      isol_p = igh_p
      isol_w = igh_w
    else
      ! isol_t and isol_d are not needed and we want the remaining values to be
      ! contiguous
      isol_t = IMDI
      isol_d = IMDI
      isol_p = 2
      isol_w = 3
    end if

  end subroutine set_solver_field_indices

end module field_indices_mod
