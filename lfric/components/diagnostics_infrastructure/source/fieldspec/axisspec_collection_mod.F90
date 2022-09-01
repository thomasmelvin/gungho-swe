!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Holds and manages a collection of axisspec objects
!>
!> @details A singleton container which holds a collection of axisspec
!>          objects that are used to manage information about axes that
!>          will be created during a run

module axisspec_collection_mod

  use constants_mod,             only: i_def, l_def, r_def, str_def, i_native
  use axisspec_mod,              only: axisspec_type
  use linked_list_mod,           only: linked_list_type, linked_list_item_type
  use log_mod,                   only: log_event, log_scratch_space, &
                                       LOG_LEVEL_ERROR
  use non_spatial_dimension_mod, only: NUMERICAL, CATEGORICAL

  implicit none

  private

  !===========================================================================
  !> @brief Holds a collection of axisspecs in a linked list
  type, public :: axisspec_collection_type

    private

    ! List of the global_mesh_type objects in this collection.
    type(linked_list_type) :: axisspec_list

  contains
    procedure :: unique_id_in_use
    procedure :: generate_and_add_axisspec
    procedure :: add_axisspec
    procedure :: get_axisspec
    procedure :: get_length
    procedure :: clear
    final     :: axisspec_collection_destructor

  end type axisspec_collection_type

  interface axisspec_collection_type
    module procedure axisspec_collection_constructor
  end interface

  !> @brief Singleton instance of an axisspec_collection_type object.
  !>
  type(axisspec_collection_type), public, allocatable, target :: &
      axisspec_collection

contains

  !===========================================================================
  !> @brief Constructs an empty axisspec collection object.
  !> @return The constructed axisspec collection object.
  !>
  function axisspec_collection_constructor() result( self )

    implicit none

    type(axisspec_collection_type), pointer :: self

    if(.not. allocated(axisspec_collection)) then
      allocate(axisspec_collection)
      axisspec_collection%axisspec_list = linked_list_type()
    end if

    self => axisspec_collection

  end function axisspec_collection_constructor

  !===========================================================================
  !> @brief Tears down object prior to being freed.
  subroutine axisspec_collection_destructor( self )

    implicit none

    type (axisspec_collection_type), intent(inout) :: self

    call self%clear()

    return
  end subroutine axisspec_collection_destructor


  !===========================================================================
  !> @brief Throws an error if the given unique_id is already used by an
  !>        axisspec in the collection
  !> @param [in] unique_id The unique_id to check
  !>
  function unique_id_in_use( self, unique_id ) result(id_in_use)

    implicit none

    class(axisspec_collection_type),        intent(inout) :: self
    character(len=*),                       intent(in)    :: unique_id

    logical(l_def)                                        :: id_in_use

    id_in_use = associated( self%get_axisspec( trim(unique_id) ) )

    return
  end function unique_id_in_use

  !===========================================================================
  !> @brief Creates a new axisspec object and adds to the collection
  !>
  !> Should contain only one axis definition, either:
  !>   numeric_axis_def for NUMERICAL axes (a real array of axis values)
  !>   label_axis_def for CATEGORICAL axes (a character array of labels)
  !>
  !> @param[in] unique_id The unique id of the axis
  !> @param[in] axis_category Type of the axis definition. This should
  !>                                be either NUMERICAL or CATEGORICAL
  !> @param[in] numeric_axis_def Real array for defining NUMERICAL axes
  !> @param[in] label_axis_def Character array for defining CATEGORICAL axes
  subroutine generate_and_add_axisspec(self,               &
                                       unique_id,          &
                                       axis_category,      &
                                       numeric_axis_def,   &
                                       label_axis_def)

    implicit none

    class(axisspec_collection_type),      intent(inout) :: self
    character(*),                         intent(in)    :: unique_id
    integer(i_native),                    intent(in)    :: axis_category
    real(r_def),        optional,         intent(in)    :: numeric_axis_def(:)
    character(str_def), optional,         intent(in)    :: label_axis_def(:)

    type(axisspec_type) :: new_axisspec

    if (present(numeric_axis_def) .and. present(label_axis_def)) then
      call log_event("Axis cannot have both a 'numeric' and a 'label' axis_def", &
                     LOG_LEVEL_ERROR)
    end if

    ! Create new axisspec object
    if (axis_category == NUMERICAL) then
      new_axisspec = axisspec_type(unique_id, &
                                   NUMERICAL, &
                                   numeric_axis_def=numeric_axis_def)
    else if (axis_category == CATEGORICAL) then
      new_axisspec = axisspec_type(unique_id, &
                                   CATEGORICAL, &
                                   label_axis_def=label_axis_def)
    else
      call log_event("Unrecognised axis category", LOG_LEVEL_ERROR)
    end if

    call self%add_axisspec( new_axisspec )

    return
  end subroutine generate_and_add_axisspec

  !===========================================================================
  !> @brief Adds the given axisspec object to the collection
  !> @param[in] axisspec The axisspec to be added to the collection
  !>
  subroutine add_axisspec( self, axisspec )

    implicit none

    class(axisspec_collection_type), intent(inout) :: self
    type(axisspec_type),             intent(in)    :: axisspec

    character(str_def)                                  :: unique_id

    unique_id = axisspec%get_unique_id()

    if (.not. self%unique_id_in_use( unique_id )) then
      ! Add it to the list
      call self%axisspec_list%insert_item( axisspec )
    else
      call log_event("The axis id '" // unique_id // "' is already in use", &
                     LOG_LEVEL_ERROR)
    end if

    return
  end subroutine add_axisspec

  !===========================================================================
  !> @brief Requests an axisspec object by unique id
  !>
  !> @param[in] unique_id The unique_id of axisspec object
  !>
  !> @return A pointer to the axisspec object
  !>
  function get_axisspec( self, unique_id ) result( axisspec )

    implicit none

    class(axisspec_collection_type)     :: self
    character(*),   intent(in)          :: unique_id

    type(axisspec_type), pointer :: axisspec

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: loop => null()

    ! nullify this here as you can't do it on the declaration because it's also on the function def
    axisspec => null()

    ! start at the head of the mesh collection linked list
    loop => self%axisspec_list%get_head()

    do
      ! If list is empty or we're at the end of list and we didn't find the
      ! unique id, return a null pointer
      if ( .not. associated(loop) ) then
        nullify(axisspec)
        exit
      end if

      ! 'cast' to the axisspec_type and then we can get at the
      ! information in the payload
      select type( m => loop%payload )
        type is (axisspec_type)
          axisspec => m
      end select

      ! Now check for the id we want and keep looping until we
      ! find it or get to the end of the list
      if ( trim(unique_id) /= trim(axisspec%get_unique_id()) ) then

        nullify(axisspec)
        loop => loop%next

      else

        exit

      end if

    end do


    nullify(loop)

  end function get_axisspec

  !===========================================================================
  !> @brief Returns the number of axisspec objects stored in this collection
  !> @return The number of axisspec objects stored in this collection
  !>
  function get_length( self ) result( count )

    implicit none

    class(axisspec_collection_type)    :: self
    integer(i_def)                      :: count

    count = self%axisspec_list%get_length()

  end function get_length

  !===========================================================================
  !> @brief Forced clear of all the axisspec objects in the collection.
  !>
  !> This routine should not need to be called manually except (possibly) in
  !> pfunit tests
  !>
  subroutine clear( self )

    ! Clear all items from the linked list in the collection
    implicit none

    class(axisspec_collection_type), intent(inout) :: self

    call self%axisspec_list%clear()

    return
  end subroutine clear

end module axisspec_collection_mod
