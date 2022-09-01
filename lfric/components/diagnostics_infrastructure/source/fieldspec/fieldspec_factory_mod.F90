!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A factory for generating a fieldspec object
!>
!> @details A fieldspec contains all relevant information to create a field


module fieldspec_factory_mod

  use constants_mod,        only: i_def, str_def, l_def
  use fieldspec_mod,        only: fieldspec_type
  use io_driver_enum_mod,   only: WRITE_FIELD_FACE

  implicit none

  private

  !> @brief Factory used to build up fieldspec objects while reading XML
  !>        before adding them to the collection
  type, public :: fieldspec_factory_type

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

  contains

    !> clear field properties
    procedure, public :: initialise

    !> create and return fieldspec object
    procedure, public :: finalise

    !> setter to return the unique_id
    procedure, public :: set_unique_id

    !> setter to return the field_group_id
    procedure, public :: set_field_group_id

    !> setter to return the mesh_id
    procedure, public :: set_mesh_id

    !> setter to return the function_space
    procedure, public :: set_function_space

    !> setter to return the order
    procedure, public :: set_order

    !> setter to return the kind
    procedure, public :: set_kind

    !> setter to return the type
    procedure, public :: set_type

    !> setter to return the io_driver
    procedure, public :: set_io_driver

    !> setter to return checksum
    procedure, public :: set_checksum

  end type fieldspec_factory_type

  interface fieldspec_factory_type
    module procedure fieldspec_factory_constructor
  end interface


contains


  !> @brief Constructs a <code>fieldspec_factory_type</code> object.
  !>
  !> @return self The fieldspec_factory object
  !>
  function fieldspec_factory_constructor() result(self)

    implicit none

    type(fieldspec_factory_type), target :: self

  end function fieldspec_factory_constructor


  !> @brief Sets all field properties to null values, readying the
  !> factory to be populated
  !>
  subroutine initialise( self )

    implicit none

    class(fieldspec_factory_type), intent(inout) :: self

    self%unique_id           = ""
    self%field_group_id      = ""
    self%mesh_id             = 0_i_def
    self%function_space      = 0_i_def
    self%order               = 0_i_def
    self%field_kind          = 0_i_def
    self%field_type          = 0_i_def
    self%io_driver           = WRITE_FIELD_FACE  ! default to write_field_face
    self%checksum            = .false.

    return
  end subroutine initialise

  !> @brief Returns a new fieldspec object using the newly populated field properties
  !> @return Fieldspec object created from the fieldspec_factory values
  !>
  function finalise( self ) result( new_fieldspec )

    implicit none

    class(fieldspec_factory_type),    intent(in)    :: self
    type(fieldspec_type)                            :: new_fieldspec

    new_fieldspec = fieldspec_type( self%unique_id,      &
                                    self%field_group_id, &
                                    self%mesh_id,        &
                                    self%function_space, &
                                    self%order,          &
                                    self%field_kind,     &
                                    self%field_type,     &
                                    self%io_driver,      &
                                    self%checksum)
  end function finalise


  !> Setter for the unique_id
  !> @param[in] unique_id
  subroutine set_unique_id( self, unique_id )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    character(len=*),                 intent(in)       :: unique_id

    self%unique_id = unique_id

    return
  end subroutine set_unique_id

  !> Setter for the field_group_id
  !> @param[in] field_group_id
  subroutine set_field_group_id( self, field_group_id )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    character(len=*),                 intent(in)       :: field_group_id

    self%field_group_id = field_group_id

    return
  end subroutine set_field_group_id

  !> Setter for the mesh_id
  !> @param[in] mesh_id
  subroutine set_mesh_id( self, mesh_id )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    integer(i_def),                   intent(in)       :: mesh_id

    self%mesh_id = mesh_id

    return
  end subroutine set_mesh_id

  !> Setter for the function_space
  !> @param[in] function_space
  subroutine set_function_space( self, function_space )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    integer(i_def),                   intent(in)       :: function_space

    self%function_space = function_space

    return
  end subroutine set_function_space

  !> Setter for the order
  !> @param[in] order
  subroutine set_order( self, order )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    integer(i_def),                   intent(in)       :: order

    self%order = order

    return
  end subroutine set_order

  !> Setter for the field_kind
  !> @param[in] field_kind
  subroutine set_kind( self, field_kind )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    integer(i_def),                   intent(in)       :: field_kind

    self%field_kind = field_kind

    return
  end subroutine set_kind

  !> Setter for the field_type
  !> @param[in] field_type
  subroutine set_type( self, field_type )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    integer(i_def),                   intent(in)       :: field_type

    self%field_type = field_type

    return
  end subroutine set_type

  !> Setter for the io_driver
  !> @param[in] io_driver
  subroutine set_io_driver( self, io_driver )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    integer(i_def),                   intent(in)       :: io_driver

    self%io_driver = io_driver

    return
  end subroutine set_io_driver

  !> Setter for checksum
  !> @param[in] checksum
  subroutine set_checksum( self, checksum )

    implicit none

    class(fieldspec_factory_type),    intent(inout)    :: self
    logical(l_def),                   intent(in)       :: checksum

    self%checksum = checksum

    return
  end subroutine set_checksum

end module fieldspec_factory_mod
