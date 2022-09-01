!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Single level operator for the pressure equation
!>
module abstract_field_operator_mod
  use field_mod,                    only : field_type
  use function_space_mod,           only: function_space_type
  use constants_mod,                only : i_def, r_def
  implicit none

  !>@brief Abstract single level operator type
  !>
  !>@details Abstract type which allows the application of an operator
  !> \f$ A \f$ to a field and division by the diagonal of this operator, i.e
  !> it implements the operations \f$ y = Ax\f$ and
  !> \f$ y = D^{-1}x\f$ (given \f$x\f$) where \f$D\f$ is the diagonal of
  !> \f$A\f$.
  type, public, abstract :: abstract_field_operator_type
     private
   contains
     procedure (apply_interface),          deferred :: apply
     procedure (apply_inv_diag_interface), deferred :: apply_inv_diag
  end type abstract_field_operator_type

  abstract interface
     !> @brief Abstract interface for the apply method.
     !>
     !> @param[in] self a single level operator
     !> @param[in] x a field the linear operator is applied to
     !> @param[inout] y a field, the result.
     subroutine apply_interface(self, x, y)
       import :: abstract_field_operator_type
       import :: field_type
       class(abstract_field_operator_type), intent(inout) :: self
       class(field_type),          intent(in)    :: x
       class(field_type),          intent(inout) :: y
     end subroutine apply_interface
  end interface

  abstract interface
     !> @brief Abstract interface defined for application of the inverse diagonal
     !>
     !> @param[in] self a single level operator
     !> @param[in] x a field the linear operator is applied to
     !> @param[inout] y a field, the result.

     subroutine apply_inv_diag_interface(self, x, y)
       import :: abstract_field_operator_type
       import :: field_type
       class(abstract_field_operator_type), intent(inout) :: self
       class(field_type),          intent(in)    :: x
       class(field_type),          intent(inout) :: y
     end subroutine apply_inv_diag_interface
  end interface

end module abstract_field_operator_mod
