!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Holds and manages a collection of fieldspec objects
!>
!> @details A singleton container which holds a collection of fieldspec
!>          objects that are used to manage information about fields that
!>          will be created during a run

module fieldspec_collection_mod

  use constants_mod,             only: i_def, l_def
  use fieldspec_mod,             only: fieldspec_type
  use linked_list_mod,           only: linked_list_type, linked_list_item_type
  use log_mod,                   only: log_event, log_scratch_space, &
                                       LOG_LEVEL_ERROR

  implicit none

  private

  !===========================================================================
  !> @brief Holds a collection of fieldspecs in a linked list
  type, public :: fieldspec_collection_type

    private

    ! List of the global_mesh_type objects in this collection.
    type(linked_list_type) :: fieldspec_list

  contains
    procedure :: check_unique_id_in_use
    procedure :: generate_and_add_fieldspec
    procedure :: add_fieldspec
    procedure :: get_fieldspec
    procedure :: get_iterator
    procedure :: get_length
    procedure :: clear
    final     :: fieldspec_collection_destructor

  end type fieldspec_collection_type

  interface fieldspec_collection_type
    module procedure fieldspec_collection_constructor
  end interface

  !> @brief Singleton instance of a fieldspec_collection_type object.
  !>
  type(fieldspec_collection_type), public, allocatable, target :: &
      fieldspec_collection

  !===========================================================================
  !> @brief Iterates through a fieldspec collection
  type, public :: fieldspec_collection_iterator_type

    private

    !> A pointer to the fieldspec within the collection that will be
    !> the next to be returned
    type(linked_list_item_type), pointer :: current

  contains
    procedure, public :: next
    procedure, public :: has_next
  end type fieldspec_collection_iterator_type

  interface fieldspec_collection_iterator_type
    module procedure fieldspec_collection_iterator_constructor
  end interface
  !===========================================================================

contains

  !===========================================================================
  !> @brief Constructs an empty fieldspec collection object.
  !> @return The constructed fieldspec collection object.
  !>
  function fieldspec_collection_constructor() result( self )

    implicit none

    type(fieldspec_collection_type), pointer :: self

    if(.not. allocated(fieldspec_collection)) then
      allocate(fieldspec_collection)
      fieldspec_collection%fieldspec_list = linked_list_type()
    end if

    self => fieldspec_collection

  end function fieldspec_collection_constructor

  !===========================================================================
  !> @brief Tears down object prior to being freed.
  subroutine fieldspec_collection_destructor( self )

    implicit none

    type (fieldspec_collection_type), intent(inout) :: self

    call self%clear()

    return
  end subroutine fieldspec_collection_destructor


  !===========================================================================
  !> @brief Returns true if a given unique_id is already used by a
  !>        fieldspec in the collection
  !> @param [in] unique_id The unique_id to check
  !>
  function check_unique_id_in_use( self, unique_id ) result (fieldspec_exists)

    implicit none

    class(fieldspec_collection_type),   intent(inout) :: self
    character(len=*),                   intent(in)    :: unique_id

    logical(l_def)                                    :: fieldspec_exists

    fieldspec_exists = ( associated( self%get_fieldspec( trim(unique_id) ) ) )

  end function check_unique_id_in_use

  !===========================================================================
  !> @brief Creates a new fieldspec object and adds to the collection
  !> @param[in] unique_id The unique id of the field
  !> @param[in] field_group_id The id of the field's group
  !> @param[in] mesh_id The mesh id of the field
  !> @param[in] function_space The function space to create the field with
  !> @param[in] order The element element order of the function space
  !> @param[in] field_kind The kind of the field
  !> @param[in] field_type The data type of the field
  !> @param[in] io_driver The io driver used for the field
  !> @param[in] checksum Whether a checksum will be written for the field
  subroutine generate_and_add_fieldspec( self,            &
                                         unique_id,       &
                                         field_group_id,  &
                                         mesh_id,         &
                                         function_space,  &
                                         order,           &
                                         field_kind,      &
                                         field_type,      &
                                         io_driver,       &
                                         checksum )

    implicit none

    class(fieldspec_collection_type),  intent(inout) :: self
    character(*),                      intent(in)    :: unique_id
    character(*),                      intent(in)    :: field_group_id
    integer(i_def),                    intent(in)    :: mesh_id
    integer(i_def),                    intent(in)    :: function_space
    integer(i_def),                    intent(in)    :: order
    integer(i_def),                    intent(in)    :: field_kind
    integer(i_def),                    intent(in)    :: field_type
    integer(i_def),                    intent(in)    :: io_driver
    logical(l_def),                    intent(in)    :: checksum

    type(fieldspec_type) :: new_fieldspec

    ! Create new fieldspec object
    new_fieldspec = fieldspec_type( unique_id,       &
                                    field_group_id,  &
                                    mesh_id,         &
                                    function_space,  &
                                    order,           &
                                    field_kind,      &
                                    field_type,      &
                                    io_driver,       &
                                    checksum )

    call self%add_fieldspec( new_fieldspec )

    return
  end subroutine generate_and_add_fieldspec

  !===========================================================================
  !> @brief Adds the given fieldspec object to the collection
  !> @param[in] fieldspec The fieldspec to be added to the collection
  !>
  subroutine add_fieldspec( self, fieldspec )

    implicit none

    class(fieldspec_collection_type), intent(inout) :: self
    type(fieldspec_type),             intent(in)    :: fieldspec

    if ( self%check_unique_id_in_use( fieldspec%get_unique_id() ) ) then
      write(log_scratch_space, '(3A)') &
            'The field unique id "', trim( fieldspec%get_unique_id() ), &
            '" is already in use'

      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! Add it to the list
    call self%fieldspec_list%insert_item( fieldspec )

    return
  end subroutine add_fieldspec

  !===========================================================================
  !> @brief Requests a fieldspec object by unique id
  !>
  !> @param[in] unique_id The unique_id of fieldspec object
  !>
  !> @return A pointer to the fieldspec object
  !>
  function get_fieldspec( self, unique_id ) result( fieldspec )

    implicit none

    class(fieldspec_collection_type)    :: self
    character(*),   intent(in)          :: unique_id

    type(fieldspec_type), pointer :: fieldspec

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: loop => null()

    ! nullify this here as you can't do it on the declaration because it's also on the function def
    fieldspec => null()

    ! start at the head of the mesh collection linked list
    loop => self%fieldspec_list%get_head()

    do
      ! If list is empty or we're at the end of list and we didn't find the
      ! unique id, return a null pointer
      if ( .not. associated(loop) ) then
        nullify(fieldspec)
        exit
      end if

      ! 'cast' to the fieldspec_type and then we can get at the
      ! information in the payload
      select type( m => loop%payload )
        type is (fieldspec_type)
          fieldspec => m
      end select

      ! Now check for the id we want and keep looping until we
      ! find it or get to the end of the list
      if ( trim(unique_id) /= trim(fieldspec%get_unique_id()) ) then

        nullify(fieldspec)
        loop => loop%next

      else

        exit

      end if

    end do


    nullify(loop)

  end function get_fieldspec

  !===========================================================================
  !> @brief Returns an iterator on the fieldspec collection
  !> @return An iterator on the fieldspec collection
  function get_iterator(self) result(iterator)

    implicit none

    class(fieldspec_collection_type), intent(in) :: self
    type(fieldspec_collection_iterator_type) :: iterator

    iterator=fieldspec_collection_iterator_type(self)

  end function get_iterator

  !===========================================================================
  !> @brief Returns the number of fieldspec objects stored in this collection
  !> @return The number of fieldspec objects stored in this collection
  !>
  function get_length( self ) result( count )

    implicit none

    class(fieldspec_collection_type)    :: self
    integer(i_def)                      :: count

    count = self%fieldspec_list%get_length()

  end function get_length

  !===========================================================================
  !> @brief Forced clear of all the fieldspec objects in the collection.
  !>
  !> This routine should not need to be called manually except (possibly) in
  !> pfunit tests
  !>
  subroutine clear( self )

    ! Clear all items from the linked list in the collection
    implicit none

    class(fieldspec_collection_type), intent(inout) :: self

    call self%fieldspec_list%clear()

    return
  end subroutine clear

  !============================================================================
  !> @brief Constructor for a fieldspec collection iterator
  !> @param [in] collection The collection to iterate over
  function fieldspec_collection_iterator_constructor(collection) result(self)

    implicit none

    type(fieldspec_collection_type) :: collection
    type(fieldspec_collection_iterator_type) :: self

    ! Start the iterator at the beginning of the fieldspec list.
    self%current => collection%fieldspec_list%get_head()
    if (.not.associated(self%current)) then
      write(log_scratch_space, '(A)') &
         'Cannot create an iterator on an empty fieldspec collection'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end function fieldspec_collection_iterator_constructor

  !============================================================================
  !> @brief Returns the next fieldspec from the collection
  !> @return fieldspec The next fieldspec from the collection
  function next(self) result (fieldspec)

    implicit none

    class(fieldspec_collection_iterator_type), intent(inout), target :: self
    class(fieldspec_type), pointer :: fieldspec

    ! 'cast' to the fieldspec_type
    select type(listfieldspec => self%current%payload)
      type is (fieldspec_type)
        fieldspec => listfieldspec
    end select

    ! Move the current fieldspec pointer onto the next fieldspec in the collection
    self%current => self%current%next

  end function next

  !============================================================================
  !> Checks if there are any further fieldspecs in the collection being
  !> iterated over
  !> @return next Logical showing if there is another fieldspec in the collection
  function has_next(self) result(next)

    implicit none

    class(fieldspec_collection_iterator_type), intent(in) :: self
    logical(l_def) :: next

    next = .true.
    if (.not.associated(self%current)) next = .false.

  end function has_next

end module fieldspec_collection_mod
