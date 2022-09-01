!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Define enumerator variables that describe the different types of field type.
!>
!> @details Enumerator variables that describe the different types of field types
!>          that can be used to construct fieldspec objects

module field_type_enum_mod

  use constants_mod,      only: i_native, str_short
  use log_mod,            only: log_event, log_scratch_space, log_level_error

  implicit none

  private
  public :: name_from_field_type, field_type_from_name

  character(*), private, parameter :: module_name = 'field_type_enum_mod'
  !-------------------------------------------------------------------------------
  ! Module parameters
  !-------------------------------------------------------------------------------
  integer(i_native), public, parameter :: integer_type = 475
  integer(i_native), public, parameter :: real_type = 721
  integer(i_native), public, parameter :: field_type_enumerator(2) &
          = [ integer_type, real_type ]
  character(str_short), public, parameter :: field_type_name(2) &
          = [character(str_short) :: 'INTEGER_TYPE', 'REAL_TYPE']

contains

  !> Gets the name corresponding to a particular field_type identifier.
  !>
  !> @param[in] field_type One of the field_type enumerations.
  !>
  !> @return String holding the field_type name.
  !>
  function name_from_field_type( field_type )
    implicit none

    integer(i_native), intent(in) :: field_type
    character(:), allocatable     :: name_from_field_type
    integer(i_native)             :: field_type_index

    field_type_index = 1
    do
      if (field_type_enumerator(field_type_index) == field_type) then
        name_from_field_type = trim(field_type_name(field_type_index))
        return
      end if
      field_type_index = field_type_index + 1
      if (field_type_index > ubound(field_type_enumerator, 1)) then
        write(log_scratch_space, &
                '(A, ": Unrecognised field_type: ",I0)') module_name, field_type
        call log_event(log_scratch_space, log_level_error)
      end if
    end do

  end function name_from_field_type

  !> Gets the field_type identifier corresponding to a particular name.
  !>
  !> @param[in] name String holding the field_type name.
  !>
  !> @return One of the field_type enumerations.
  !>
  function field_type_from_name( name )
    implicit none

    character(*), intent(in) :: name
    integer(i_native)        :: field_type_from_name
    integer(i_native)        :: field_type_index

    field_type_index = 1
    do
      if (field_type_name(field_type_index) == name) then
        field_type_from_name = field_type_enumerator(field_type_index)
        return
      end if
      field_type_index = field_type_index + 1
      if (field_type_index > ubound(field_type_enumerator, 1)) then
        call log_event("Unknown field type " // name, log_level_error)
      end if
    end do

  end function field_type_from_name

end module field_type_enum_mod