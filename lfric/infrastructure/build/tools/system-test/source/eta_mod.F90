module eta_mod
  implicit none
  private
  public :: pythagorus
contains
  real pure function pythagorus(s1, s2)
    implicit none
    integer, intent(in) :: s1
    integer, intent(in) :: s2
    pythagorus = sqrt(real(s1)**2 + real(s2)**2)
  end function pythagorus
end module eta_mod
