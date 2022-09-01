!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing operator related classes.
!>
!> @details Implements the locally assembled operator (i.e. the stencil is
!>          assembled in each cell of the 3d grid)

module operator_mod

  ! Eventually the precision of the operator will be set in a module held
  ! within the model (as it is model information). For now, PSyclone is
  ! expecting to "use" the definitions from operator_mod, so it is set here
#if (RDEF_PRECISION == 32)
  use operator_r32_mod, only: operator_type       => operator_r32_type, &
                              operator_proxy_type => operator_r32_proxy_type
#else
  use operator_r64_mod, only: operator_type       => operator_r64_type, &
                              operator_proxy_type => operator_r64_proxy_type
#endif

#if (R_SOLVER_PRECISION == 32)
  use operator_r32_mod, only: r_solver_operator_type => operator_r32_type, &
                        r_solver_operator_proxy_type => operator_r32_proxy_type
#else
  use operator_r64_mod, only: r_solver_operator_type => operator_r64_type, &
                        r_solver_operator_proxy_type => operator_r64_proxy_type
#endif

  implicit none
! Removing the following "private" statement is a workaround for a bug that
! appeared in Intel v19. Every item in the module has an explicit access set,
! so not setting the default has no effect. See ticket #3326 for details
!  private
  public :: operator_type, operator_proxy_type
  public :: r_solver_operator_type, r_solver_operator_proxy_type

end module operator_mod
