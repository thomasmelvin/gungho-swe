module epsilon_mod
  implicit none
  private
  type, public :: cart_type
    private
    integer :: x
    integer :: y
  contains
    private
    procedure, public :: distance
  end type cart_type
  interface cart_type
    procedure new_cart_type
  end interface cart_type
  interface
    module function new_cart_type(x,y)
      implicit none
      integer, intent(in) :: x
      integer, intent(in) :: y
      type(cart_type) :: new_cart_type
    end function new_cart_type
    real pure module function distance(this, other)
      implicit none
      class(cart_type), intent(in) :: this
      type(cart_type),  intent(in) :: other
    end function distance
  end interface
end module epsilon_mod
