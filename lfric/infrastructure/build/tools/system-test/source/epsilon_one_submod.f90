submodule (epsilon_mod) epsilon_one_submod
  use eta_mod, only : pythagorus
  implicit none
contains
  real pure module function distance(this, other)
    implicit none
    class(cart_type), intent(in) :: this
    type(cart_type),  intent(in) :: other
    distance = pythagorus(other%x - this%x, other%y - this%y)
  end function distance
end submodule epsilon_one_submod
