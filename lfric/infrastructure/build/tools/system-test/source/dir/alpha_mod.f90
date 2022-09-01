module alpha_mod
  use eta_mod, only: pythagorus
  implicit none
  private
  public :: alpha_sub
contains
  subroutine alpha_sub(fred)
    implicit none
    real,intent(out) :: fred
    fred = pythagorus(2,4)
  end subroutine alpha_sub
end module alpha_mod
