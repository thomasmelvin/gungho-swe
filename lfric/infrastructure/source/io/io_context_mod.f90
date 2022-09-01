!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Specifies the interface for all I/O context classes.
!>
module io_context_mod

  use clock_mod,     only : clock_type
  use constants_mod, only : i_native, r_second
  use field_mod,     only : field_type
  use file_mod,      only : file_type

  implicit none

  private

  !> @brief All context classes inherit this interface.
  !>
  type, public, abstract :: io_context_type
    private
  contains
    private
    procedure(initialise_if),   public, deferred :: initialise
    procedure(get_clock_if),    public, deferred :: get_clock
    procedure(get_filelist_if), public, deferred :: get_filelist
  end type io_context_type

  abstract interface
    !> Initialises the context object
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
    subroutine initialise_if( this, id, communicator,  &
                              chi, panel_id,           &
                              start_time, finish_time, &
                              spinup_period,           &
                              seconds_per_step,        &
                              calendar_start,          &
                              calendar_type,           &
                              list_of_files )
      import io_context_type, field_type, file_type, r_second, i_native
      implicit none
      class(io_context_type),     intent(inout) :: this
      character(*),               intent(in)    :: id
      integer(i_native),          intent(in)    :: communicator
      class(field_type),          intent(in)    :: chi(:), panel_id
      character(*),               intent(in)    :: start_time, finish_time
      real(r_second),             intent(in)    :: spinup_period, seconds_per_step
      character(*),               intent(in)    :: calendar_start, calendar_type
      class(file_type), optional, intent(in)    :: list_of_files(:)
    end subroutine initialise_if
  end interface

  abstract interface
    !> Gets the clock associated with this context.
    !>
    !> @return Clock object.
    !>
    function get_clock_if( this ) result(clock)
      import clock_type, io_context_type
      implicit none
      class(io_context_type), intent(in), target :: this
      class(clock_type), pointer :: clock
    end function get_clock_if
  end interface

  abstract interface
    !> Gets the context's file list
    !>
    !> @return The list of file objects used by the model
    !>
    function get_filelist_if( this ) result(filelist)
      import file_type, io_context_type
      implicit none
      class(io_context_type), intent(in), target :: this
      class(file_type), pointer :: filelist(:)
    end function get_filelist_if
  end interface

contains

end module io_context_mod
