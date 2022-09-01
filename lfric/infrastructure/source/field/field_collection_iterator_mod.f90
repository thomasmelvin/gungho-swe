!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides functionality for iterating over all members of a field
!>        collection
!>
!> @details Provides functionality for iteratively returning every member
!>          of a field collection. The order of the fields returned is
!>          not defined and can change if the implementation of the field
!>          collection is changed
!
module field_collection_iterator_mod

  use constants_mod,           only: i_def, l_def
  use field_collection_mod,    only: field_collection_type
  use field_mod,               only: field_type, &
                                     field_pointer_type
  use field_r32_mod,           only: field_r32_type, &
                                     field_r32_pointer_type
  use field_r64_mod,           only: field_r64_type, &
                                     field_r64_pointer_type
  use field_parent_mod,        only: field_parent_type
  use integer_field_mod,       only: integer_field_type, &
                                     integer_field_pointer_type
  use log_mod,                 only: log_event, log_scratch_space, &
                                     LOG_LEVEL_ERROR
  use linked_list_data_mod,    only: linked_list_data_type
  use linked_list_mod,         only: linked_list_type, &
                                     linked_list_item_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type that iterates through a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_iterator_type
    private
    !> Dummy allocatable - workaround for gcc bug (ref:61767)
    integer(i_def), allocatable :: dummy_for_gnu
    !> A pointer to the field collection being iterated over
    type(field_collection_type), pointer :: collection
    !> A pointer to the linked list item within the collection that will
    !> contain the next field to be returned
    type(linked_list_item_type), pointer :: current
  contains
    procedure, public :: initialise => initialise_iter
    procedure, public :: next
    procedure, public :: has_next
  end type field_collection_iterator_type

  !-----------------------------------------------------------------------------
  ! Type that iterates through the real fields in a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_real_iterator_type
    private
    !> Dummy allocatable - workaround for gcc bug (ref:61767)
    integer(i_def), allocatable :: dummy_for_gnu
    !> A pointer to the field collection being iterated over
    type(field_collection_type), pointer :: collection
    !> A pointer to the linked list item within the collection that will
    !> contain the next real field to be returned
    type(linked_list_item_type), pointer :: current
  contains
    procedure, public :: initialise => initialise_real_iter
    procedure, public :: next => next_real
    procedure, public :: has_next => has_next_real
  end type field_collection_real_iterator_type

  !-----------------------------------------------------------------------------
  ! Type that iterates through the integer fields in a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_integer_iterator_type
    private
    !> Dummy allocatable - workaround for gcc bug (ref:61767)
    integer(i_def), allocatable :: dummy_for_gnu
    !> A pointer to the field collection being iterated over
    type(field_collection_type), pointer :: collection
    !> A pointer to the linked list item within the collection that will
    !> contain the next integer field to be returned
    type(linked_list_item_type), pointer :: current
  contains
    procedure, public :: initialise => initialise_integer_iter
    procedure, public :: next => next_integer
    procedure, public :: has_next => has_next_integer
  end type field_collection_integer_iterator_type

contains

!> Initialise a field collection iterator
!> @param [in] collection The collection to iterate over
subroutine initialise_iter(self, collection)

  implicit none

  class(field_collection_iterator_type) :: self
  type(field_collection_type), target :: collection

  ! Store a pointer to the collection being iterated over
  self%collection => collection

  ! Start the iterator at the beginning of the field list.
  nullify(self%current)
  self%current => self%collection%get_next_item(self%current)

  if(.not.associated(self%current))then
    write(log_scratch_space, '(2A)') &
       'Cannot create an iterator on an empty field collection: ', &
        trim(self%collection%get_name())
    call log_event( log_scratch_space, LOG_LEVEL_ERROR)
  end if

end subroutine initialise_iter

!> Initialise an iterator over the real fields in a field collection
!> @param [in] collection The collection to iterate over
subroutine initialise_real_iter(self, collection)

  implicit none

  class(field_collection_real_iterator_type) :: self
  type(field_collection_type), target :: collection

  ! Store a pointer to the collection being iterated over
  self%collection => collection

  ! Start the iterator at the beginning of the field list.
  nullify(self%current)
  self%current => self%collection%get_next_item(self%current)

  do
    if(.not.associated(self%current))then
      write(log_scratch_space, '(2A)') &
         'Cannot create a r_def real field iterator on field collection: ', &
          trim(self%collection%get_name())
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Make sure first field pointed to in list is a real field
    select type(listfield => self%current%payload)
      type is (field_type)
        exit
      type is (field_pointer_type)
        exit
    end select
    self%current => self%collection%get_next_item(self%current)
  end do

end subroutine initialise_real_iter


!> Initialise an iterator over the integer fields in a field collection
!> @param [in] collection The collection to iterate over
subroutine initialise_integer_iter(self, collection)

  implicit none

  class(field_collection_integer_iterator_type) :: self
  type(field_collection_type), target :: collection

  ! Store a pointer to the collection being iterated over
  self%collection => collection

  ! Start the iterator at the beginning of the field list.
  nullify(self%current)
  self%current => self%collection%get_next_item(self%current)

  do
    if(.not.associated(self%current))then
      write(log_scratch_space, '(2A)') &
         'Cannot create an integer field iterator on field collection: ', &
          trim(self%collection%get_name())
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Make sure first field pointed to in list is an integer field
    select type(listfield => self%current%payload)
      type is (integer_field_type)
        exit
      type is (integer_field_pointer_type)
        exit
    end select
    self%current => self%collection%get_next_item(self%current)
  end do

end subroutine initialise_integer_iter


!> Returns the next field from the collection
!> @return A polymorphic field pointer to either a field or field pointer
!>         that is next in the collection
function next(self) result (field)

  implicit none

  class(field_collection_iterator_type), intent(inout), target :: self
  class(field_parent_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (field_r32_type)
      field => listfield
    type is (field_r64_type)
      field => listfield
    type is (integer_field_type)
      field => listfield
    type is (field_r32_pointer_type)
      field => listfield%field_ptr
    type is (field_r64_pointer_type)
      field => listfield%field_ptr
    type is (integer_field_pointer_type)
      field => listfield%field_ptr
  end select
  ! Move the current item pointer onto the next field in the collection
  self%current => self%collection%get_next_item(self%current)

end function next

!> Returns the next real field from the collection
!> @return A field pointer to the next (r_def real) field in the collection
function next_real(self) result (field)

  implicit none

  class(field_collection_real_iterator_type), intent(inout), target :: self
  type(field_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (field_type)
      field => listfield
    type is (field_pointer_type)
      field => listfield%field_ptr
  end select

  ! Move the current item pointer onto the next real field in the collection
  self%current => self%collection%get_next_item(self%current)
  do
    if(.not.associated(self%current))then
      exit
    end if
    ! Make sure field pointed to in list is a real field
    select type(listfield => self%current%payload)
      type is (field_type)
        exit
      type is (field_pointer_type)
        exit
    end select
    self%current => self%collection%get_next_item(self%current)
  end do

end function next_real

!> Returns the next integer field from the collection
!> @return A field pointer to the next (integer) field in the collection
function next_integer(self) result (field)

  implicit none

  class(field_collection_integer_iterator_type), intent(inout), target :: self
  type(integer_field_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (integer_field_type)
      field => listfield
    type is (integer_field_pointer_type)
      field => listfield%field_ptr
  end select

  ! Move the current item pointer onto the next integer field in the collection
  self%current => self%collection%get_next_item(self%current)
  do
    if(.not.associated(self%current))then
      exit
    end if
    ! Make sure field pointed to in list is an integer field
    select type(listfield => self%current%payload)
      type is (integer_field_type)
        exit
      type is (integer_field_pointer_type)
        exit
    end select
    self%current => self%collection%get_next_item(self%current)
  end do

end function next_integer


!> Checks if there are any further fields in the collection being iterated over
!> @return next true if there is another field in the collection, and false if
!> there isn't.
function has_next(self) result(next)
  implicit none
  class(field_collection_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next

!> Checks if there are any further real fields in the collection being
!> iterated over.
!> @return next true if there is another real field in the collection, and
!> false if there isn't.
function has_next_real(self) result(next)
  implicit none
  class(field_collection_real_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next_real

!> Checks if there are any further integer fields in the collection being
!> iterated over.
!> @return next true if there is another integer field in the collection, and
!> false if there isn't.
function has_next_integer(self) result(next)
  implicit none
  class(field_collection_integer_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next_integer

end module field_collection_iterator_mod
