!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module lfric_xios_constants_mod

  use, intrinsic :: iso_fortran_env, only: int16
  use               constants_mod,   only: dp_native, i_def

  implicit none

  private
  public :: dp_xios, xios_max_int

  !< XIOS kind for double precision fields
  integer, parameter :: dp_xios = dp_native

  !< The largest integer that can be output by XIOS
  integer(kind=i_def), parameter :: xios_max_int = huge(0_int16)

end module lfric_xios_constants_mod