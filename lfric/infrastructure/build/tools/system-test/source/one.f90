program one
  use iso_fortran_env, only : output_unit
  use alpha_mod, only: alpha_sub
  use gamma_mod, only: gamma_print
  use epsilon_mod, only : cart_type
  implicit none
  real :: something
  type(cart_type) :: first
  type(cart_type) :: second
  call alpha_sub(something)
  write(output_unit, '("From alpha: ", F5.2)') something
  call gamma_print()
  first = cart_type(1,2)
  second = cart_type(3,4)
  write(output_unit, '("From epsilon: ", F5.2)') first%distance(second)
end program one
