!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Define enumerator variables that describe the different types of io_driver.
!>
!> @details Enumerator variables that describe the different types of io_driver
!>          that can be used to construct fieldspec objects


module io_driver_enum_mod

  use constants_mod, only : i_native, str_def
  use log_mod, only : log_event, log_scratch_space, log_level_error

  implicit none
  private
  public :: name_from_io_driver, io_driver_from_name

  character(*), private, parameter :: module_name = 'io_driver_enum_mod'
  !-------------------------------------------------------------------------------
  ! Module parameters
  !-------------------------------------------------------------------------------
  integer(i_native), public, parameter :: WRITE_FIELD_FACE = 1
  integer(i_native), public, parameter :: io_driver_enumerator(1) = [WRITE_FIELD_FACE]
  character(str_def), public, parameter :: io_driver_name(1) &
          = [character(str_def) :: 'WRITE_FIELD_FACE']

contains

  !> Gets the name corresponding to a particular io_driver identifier.
  !>
  !> @param[in] io_driver One of the io_driver enumerations.
  !>
  !> @return String holding the io_driver name.
  !>
  function name_from_io_driver(io_driver)
    implicit none

    integer(i_native), intent(in) :: io_driver
    character(str_def) :: name_from_io_driver
    integer(i_native) :: io_driver_index

    io_driver_index = 1
    do
      if (io_driver_enumerator(io_driver_index) == io_driver) then
        name_from_io_driver = io_driver_name(io_driver_index)
        return
      end if

      io_driver_index = io_driver_index + 1
      if (io_driver_index > ubound(io_driver_enumerator, 1)) then
        write(log_scratch_space, &
                '(A, ": Unrecognised io_driver: ",I0)') module_name, io_driver
        call log_event(log_scratch_space, log_level_error)
      end if
    end do
  end function name_from_io_driver

  !> Gets the io_driver identifier corresponding to a particular name.
  !>
  !> @param[in] name String holding the io_driver name.
  !>
  !> @return One of the io_driver enumerations.
  !>
  function io_driver_from_name(name)
    implicit none

    character(*), intent(in) :: name
    integer(i_native) :: io_driver_from_name
    integer(i_native) :: io_driver_index

    io_driver_index = 1
    do
      if (io_driver_name(io_driver_index) == name) then
        io_driver_from_name = io_driver_enumerator(io_driver_index)
        return
      end if
      io_driver_index = io_driver_index + 1
      if (io_driver_index > ubound(io_driver_name, 1)) then
        call log_event("Unknown io driver " // name, log_level_error)
      end if
    end do
  end function io_driver_from_name

end module io_driver_enum_mod