module mismatched_mod
  use delta_mod, only: delta_print
  implicit none
contains
  subroutine beta_beater()
    implicit none
    integer :: i
    do i=0, 4
      call delta_print()
    end do
  end subroutine beta_beater
end module mismatched_mod
