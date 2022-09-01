!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the fieldspec class
!>
!> @details A fieldspec contains all relevant information to create a field


module fieldspec_mod

  use constants_mod,        only: i_def, str_def, l_def
  use linked_list_mod,      only: linked_list_type, &
                                  linked_list_item_type
  use linked_list_data_mod, only: linked_list_data_type
  use axisspec_mod,         only: axisspec_type

  implicit none

  private

  type, extends(linked_list_data_type), public :: fieldspec_type

    private

    !> Unique id used by the diagnostic system to identify the field
    character(str_def)                :: unique_id
    !> Other information used to create a field
    character(str_def)                :: field_group_id
    integer(i_def)                    :: mesh_id
    integer(i_def)                    :: function_space
    integer(i_def)                    :: order
    integer(i_def)                    :: field_kind
    integer(i_def)                    :: field_type
    integer(i_def)                    :: io_driver
    logical(l_def)                    :: checksum
    type(axisspec_type), pointer      :: vertical_axis => null()
    type(linked_list_type)            :: non_spatial_dimension_list

  contains

    !> Getter to return the unique_id
    procedure, public :: get_unique_id

    !> Getter to return the field_group_id
    procedure, public :: get_field_group_id

    !> Getter to return the mesh_id
    procedure, public :: get_mesh_id

    !> Getter to return the function_space
    procedure, public :: get_function_space

    !> Getter to return the order
    procedure, public :: get_order

    !> Getter to return the kind
    procedure, public :: get_kind

    !> Getter to return the type
    procedure, public :: get_type

    !> Getter to return the io_driver
     procedure, public :: get_io_driver

    !> Getter to return checksum
     procedure, public :: get_checksum

    !> Getter to return specification for vertical axis
     procedure, public :: get_vertical_axis

    !> Getter to return specification for non-spatial dimension
     procedure, public :: get_non_spatial_dimension

    !> Getter to return iterator for non-spatial dimensions
     procedure, public :: get_non_spatial_dimension_iterator

    !> Getter to return specification for vertical axis
     procedure, public :: set_vertical_axis

    !> Getter to return specification for non-spatial dimension
     procedure, public :: add_non_spatial_dimension

  end type fieldspec_type

  interface fieldspec_type
    module procedure fieldspec_constructor
  end interface

  !---------------------------------------------------------------------------
  !> @brief Iterates through the non-spatial dimensions of the field
  type, public :: non_spatial_dimension_iterator_type

    private

    !> A pointer to the non-spatial dimension axisspec that will be
    !> the next to be returned
    type(linked_list_item_type), pointer :: current

  contains
    procedure, public :: next
    procedure, public :: has_next
  end type non_spatial_dimension_iterator_type

  interface non_spatial_dimension_iterator_type
    module procedure non_spatial_dimension_iterator_constructor
  end interface


contains


  !> Construct a <code>fieldspec_type</code> object.
  !>
  !> @param [in] unique_id A unique identifer for the field
  !> @param [in] mesh_id The mesh to create field with
  !> @param [in] function_space The function space to create field with
  !> @param [in] order The order of the field
  !> @param [in] field_kind The kind of the field
  !> @param [in] field_type The type of the field
  !> @param [in] io_driver The name of the io_driver used to write this field
  !> @param [in] checksum Logical showing if a checksum should be written for &
  !>                      this field
  !> @return self the fieldspec object
  !>
  function fieldspec_constructor( unique_id,                   &
                                  field_group_id,              &
                                  mesh_id,                     &
                                  function_space,              &
                                  order,                       &
                                  field_kind,                  &
                                  field_type,                  &
                                  io_driver,                   &
                                  checksum )                   &
                                 result(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    character(*),                      intent(in) :: unique_id
    character(*),                      intent(in) :: field_group_id
    integer(i_def),                    intent(in) :: mesh_id
    integer(i_def),                    intent(in) :: function_space
    integer(i_def),                    intent(in) :: order
    integer(i_def),                    intent(in) :: field_kind
    integer(i_def),                    intent(in) :: field_type
    integer(i_def),                    intent(in) :: io_driver
    logical(l_def),                    intent(in) :: checksum

    type(fieldspec_type), target :: self

    self%unique_id                  = trim(unique_id)
    self%field_group_id             = trim(field_group_id)
    self%mesh_id                    = mesh_id
    self%function_space             = function_space
    self%order                      = order
    self%field_kind                 = field_kind
    self%field_type                 = field_type
    self%io_driver                  = io_driver
    self%checksum                   = checksum
    self%vertical_axis             => null()
    self%non_spatial_dimension_list = linked_list_type()

  end function fieldspec_constructor


  !> Getter for unique_id
  !> @param[in]  self  fieldspec_type
  !> @return unique_id
  function get_unique_id(self) result(unique_id)

    implicit none

    class(fieldspec_type), intent(in) :: self
    character(str_def) :: unique_id

    unique_id = trim(self%unique_id)

  end function get_unique_id

  !> Getter for field_group_id
  !> @param[in]  self  fieldspec_type
  !> @return field_group_id
  function get_field_group_id(self) result(field_group_id)

    implicit none

    class(fieldspec_type), intent(in) :: self
    character(str_def) :: field_group_id

    field_group_id = trim(self%field_group_id)

  end function get_field_group_id

  !> Getter for mesh_id
  !> @param[in]  self  fieldspec_type
  !> @return mesh_id
  function get_mesh_id(self) result(mesh_id)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: mesh_id

    mesh_id = self%mesh_id

  end function get_mesh_id

  !> Getter for function_space
  !> @param[in]  self  fieldspec_type
  !> @return function_space
  function get_function_space(self) result(function_space)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: function_space

    function_space = self%function_space

  end function get_function_space

  !> Getter for order
  !> @param[in]  self  fieldspec_type
  !> @return order
  function get_order(self) result(order)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: order

    order = self%order

  end function get_order

  !> Getter for field_kind
  !> @param[in]  self  fieldspec_type
  !> @return field_kind
  function get_kind(self) result(field_kind)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: field_kind

    field_kind = self%field_kind

  end function get_kind

  !> Getter for field_type
  !> @param[in]  self  fieldspec_type
  !> @return field_type
  function get_type(self) result(field_type)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: field_type

    field_type = self%field_type

  end function get_type

  !> Getter for io_driver
  !> @param[in]  self  fieldspec_type
  !> @return io_driver
  function get_io_driver(self) result(io_driver)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: io_driver

    io_driver = self%io_driver

  end function get_io_driver

  !> Getter for checksum
  !> @param[in]  self  fieldspec_type
  !> @return checksum
  function get_checksum(self) result(checksum)

    implicit none

    class(fieldspec_type), intent(in) :: self
    logical(l_def) :: checksum

    checksum = self%checksum

  end function get_checksum

  !> Getter for vertical_axis
  !> @param[in]  self  fieldspec_type
  !> @return vertical_axis
  function get_vertical_axis(self) result(vertical_axis)

    implicit none

    class(fieldspec_type), intent(in) :: self
    type(axisspec_type),      pointer :: vertical_axis

    vertical_axis => null()

    vertical_axis => self%vertical_axis

  end function get_vertical_axis


  !> @brief Requests a non-spatial dimension axisspec object by unique id
  !>
  !> @param[in] unique_id The unique_id of axisspec object
  !>
  !> @return A pointer to the axisspec object
  !>
  function get_non_spatial_dimension( self, unique_id ) result( axisspec )

    implicit none

    class(fieldspec_type)               :: self
    character(*),   intent(in)          :: unique_id

    type(axisspec_type), pointer :: axisspec

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: loop => null()

    ! nullify this here as you can't do it on the declaration because it's also on the function def
    axisspec => null()

    ! start at the head of the mesh collection linked list
    loop => self%non_spatial_dimension_list%get_head()

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

  end function get_non_spatial_dimension

  !> Setter for vertical_axis
  !> @param[in] vertical_axis
  subroutine set_vertical_axis( self, vertical_axis )

    implicit none

    class(fieldspec_type),     intent(inout)   :: self
    type(axisspec_type), pointer, intent(in)   :: vertical_axis

    self%vertical_axis => vertical_axis

    return
  end subroutine set_vertical_axis

  !> Adds a non-spatial dimension to the linked list
  !> @param[in] non_spatial_dimension
  subroutine add_non_spatial_dimension( self, non_spatial_dimension )

    implicit none

    class(fieldspec_type),     intent(inout)   :: self
    type(axisspec_type), pointer, intent(in)   :: non_spatial_dimension

    call self%non_spatial_dimension_list%insert_item(non_spatial_dimension)

    return
  end subroutine add_non_spatial_dimension

  !> @brief Returns an iterator for the non-spatial dimensions
  !> @return An iterator on the non-spatial dimensions
  function get_non_spatial_dimension_iterator(self) result(iterator)

    implicit none

    class(fieldspec_type), intent(in) :: self
    type(non_spatial_dimension_iterator_type) :: iterator

    iterator=non_spatial_dimension_iterator_type(self)

  end function get_non_spatial_dimension_iterator

  !============================================================================
  !> @brief Constructor for a non-spatial-dimension iterator
  !> @param [in] fieldspec The fieldspec whose non-spatial dimensions are to &
  !>                       be iterated over
  function non_spatial_dimension_iterator_constructor(fieldspec) result(self)

    implicit none

    type(fieldspec_type) :: fieldspec
    type(non_spatial_dimension_iterator_type) :: self

    ! Start the iterator at the beginning of the fieldspec list.
    self%current => fieldspec%non_spatial_dimension_list%get_head()

  end function non_spatial_dimension_iterator_constructor

  !> @brief Returns the next axisspec from the non-spatial dimensions
  !> @return axisspec The next axisspec from the non-spatial dimensions
  function next(self) result (axisspec)

    implicit none

    class(non_spatial_dimension_iterator_type), intent(inout), target :: self
    class(axisspec_type), pointer :: axisspec

    ! 'cast' to the axisspec_type
    select type(listaxisspec => self%current%payload)
      type is (axisspec_type)
        axisspec => listaxisspec
    end select

    ! Move the current axisspec pointer onto the next axisspec
    self%current => self%current%next

  end function next

  !> Checks if there are any further axisspecs in the collection being
  !> iterated over
  !> @return next Logical showing if there is another axisspec in the collection
  function has_next(self) result(next)

    implicit none

    class(non_spatial_dimension_iterator_type), intent(in) :: self
    logical(l_def) :: next

    next = .true.
    if (.not.associated(self%current)) next = .false.

  end function has_next

end module fieldspec_mod
