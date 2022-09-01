!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Simple context supporting only a simple time-step clock.
!>
module simple_io_context_mod

  use clock_mod,         only : clock_type
  use constants_mod,     only : i_native, r_second
  use field_mod,         only : field_type
  use file_mod,          only : file_type
  use io_context_mod,    only : io_context_type
  use log_mod,           only : log_event, log_level_error
  use model_clock_mod,   only : model_clock_type
  use step_calendar_mod, only : step_calendar_type

  implicit none

  private

  !> @brief I/O context with no external dependencies.
  !>
  !> This context may be used where simplified I/O which does not require
  !> additional libraries is desirable.
  !>
  type, public, extends(io_context_type) :: simple_io_context_type
    private
    type(model_clock_type), allocatable :: clock
    class(file_type),       allocatable :: filelist(:)
  contains
    private
    procedure, public :: initialise
    procedure, public :: get_clock
    procedure, public :: get_filelist
  end type simple_io_context_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Set up a declared simple_io_context object.
  !>
  !> @param [in]     id                Unique identifying string.
  !> @param [in]     communicator      MPI communicator used by context.
  !> @param [in]     chi               Array of coordinate fields
  !> @param [in]     panel_id          Panel ID field
  !> @param [in]     start_time        Time of first step.
  !> @param [in]     finish_time       Time of last step.
  !> @param [in]     spinup_period     Number of seconds in spinup period.
  !> @param [in]     seconds_per_step  Number of seconds in a time step.
  !> @param [in]     calendar_start    Start date for calendar
  !> @param [in]     calendar_type     Type of calendar.
  !> @param [in]     list_of_files     List of file objects attached to the
  !!                                   context
  !>
  subroutine initialise( this,                    &
                         id, communicator,        &
                         chi, panel_id,           &
                         start_time, finish_time, &
                         spinup_period,           &
                         seconds_per_step,        &
                         calendar_start,          &
                         calendar_type,           &
                         list_of_files )

    implicit none

    class(simple_io_context_type), intent(inout) :: this
    character(*),                  intent(in)    :: id
    integer(i_native),             intent(in)    :: communicator
    class(field_type),             intent(in)    :: chi(:)
    class(field_type),             intent(in)    :: panel_id
    character(*),                  intent(in)    :: start_time
    character(*),                  intent(in)    :: finish_time
    real(r_second),                intent(in)    :: spinup_period
    real(r_second),                intent(in)    :: seconds_per_step
    character(*),                  intent(in)    :: calendar_start
    character(*),                  intent(in)    :: calendar_type
    class(file_type),    optional, intent(in)    :: list_of_files(:)

    type(step_calendar_type), allocatable :: calendar
    integer(i_native)                     :: rc

    allocate( calendar, stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate calendar", log_level_error )
    end if

    allocate( this%clock, stat=rc )
    if (rc /= 0) then
      call log_event( "Failed to allocate clock", log_level_error )
    end if
    this%clock = model_clock_type( calendar, start_time, finish_time, &
                                   seconds_per_step, spinup_period )

    ! Attach optional file list
    if (present(list_of_files)) then
        allocate( this%filelist, source=list_of_files )
    end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Gets the clock associated with this context.
  !>
  !> @return Context's clock object.
  !>
  function get_clock( this ) result(clock)

    implicit none

    class(simple_io_context_type), intent(in), target :: this
    class(clock_type), pointer :: clock

    clock => this%clock

  end function get_clock

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets the context's file list
  !>
  !> @return  The list of file objects used by the model
  !>
  function get_filelist( this ) result(filelist)

    implicit none

    class(simple_io_context_type), target, intent(in) :: this
    class(file_type), pointer :: filelist(:)

    filelist => this%filelist

  end function get_filelist

end module simple_io_context_mod
