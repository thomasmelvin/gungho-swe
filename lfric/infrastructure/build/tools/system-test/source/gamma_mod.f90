module gamma_mod
  implicit none
  private
  public :: gamma_print
contains
  subroutine gamma_print()
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    write(output_unit, '("In Gamma")')
  end subroutine gamma_print
end module gamma_mod

module delta_mod
  implicit none
  private
  public :: delta_print
contains
  subroutine delta_print()
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    write(output_unit, '("In Delta")')
  end subroutine delta_print
end module delta_mod
