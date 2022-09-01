!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Wrap the XIOS context in an object for easier management and cleaner coding.
!>
module lfric_xios_context_mod

  use clock_mod,            only : clock_type
  use constants_mod,        only : i_native, &
                                   r_second, &
                                   l_def
  use field_mod,            only : field_type
  use file_mod,             only : file_type
  use io_context_mod,       only : io_context_type
  use lfric_xios_file_mod,  only : lfric_xios_file_type
  use log_mod,              only : log_event,       &
                                   log_level_error, &
                                   log_level_info
  use step_calendar_mod,    only : step_calendar_type
  use lfric_xios_clock_mod, only : lfric_xios_clock_type
  use lfric_xios_setup_mod, only : init_xios_dimensions, &
                                   setup_xios_files
  use lfric_xios_file_mod,  only : lfric_xios_file_type
  use lfric_xios_utils_mod, only : parse_date_as_xios
  use xios,                 only : xios_context,                  &
                                   xios_context_initialize,       &
                                   xios_close_context_definition, &
                                   xios_context_finalize,         &
                                   xios_date,                     &
                                   xios_define_calendar,          &
                                   xios_get_handle,               &
                                   xios_set_current_context,      &
                                   xios_context_finalize
  use mod_wait,             only : init_wait

  implicit none

  private

  !> Manages interactions with XIOS.
  !>
  type, public, extends(io_context_type) :: lfric_xios_context_type
    private
    character(:),                 allocatable :: id
    type(xios_context)                        :: handle
    class(lfric_xios_clock_type), allocatable :: clock
    type(lfric_xios_file_type),   allocatable :: filelist(:)

  contains
    private
    procedure, public :: initialise
    procedure, public :: get_clock
    procedure, public :: get_filelist
    final :: finalise
  end type lfric_xios_context_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Set up an XIOS context object.
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
  subroutine initialise( this, id, communicator,          &
                         chi, panel_id,                   &
                         start_time, finish_time,         &
                         spinup_period, seconds_per_step, &
                         calendar_start, calendar_type,   &
                         list_of_files)

    implicit none

    class(lfric_xios_context_type), intent(inout) :: this
    character(*),                   intent(in)    :: id
    integer(i_native),              intent(in)    :: communicator
    class(field_type),              intent(in)    :: chi(:)
    class(field_type),              intent(in)    :: panel_id
    character(*),                   intent(in)    :: start_time
    character(*),                   intent(in)    :: finish_time
    real(r_second),                 intent(in)    :: spinup_period
    real(r_second),                 intent(in)    :: seconds_per_step
    character(*),                   intent(in)    :: calendar_start
    character(*),                   intent(in)    :: calendar_type
    class(file_type),     optional, intent(in)    :: list_of_files(:)

    type(step_calendar_type), allocatable :: calendar
    integer(i_native)                     :: rc
    type(xios_date)                       :: calendar_start_xios

    call xios_context_initialize( id, communicator )
    call xios_get_handle( id, this%handle )
    call xios_set_current_context( this%handle )

    ! Calendar start is adjusted when clock is initialised
    calendar_start_xios = parse_date_as_xios(trim(adjustl(calendar_start)))
    call xios_define_calendar( type=calendar_type,              &
                               time_origin=calendar_start_xios, &
                               start_date=calendar_start_xios )

    allocate( calendar, stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate calendar", log_level_error )
    end if

    allocate( this%clock,                                            &
              source=lfric_xios_clock_type( calendar,                &
                                            start_time, finish_time, &
                                            seconds_per_step,        &
                                            spinup_period ),         &
              stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate clock", log_level_error )
    end if

    ! Attach optional file list
    if (present(list_of_files)) then
      select type(list_of_files)
      type is (lfric_xios_file_type)
        this%filelist = list_of_files

      class default
        call log_event("Files attached to lfric_xios_context_type must be "// &
                       "of lfric_xios_file_type", log_level_error)
      end select
    else
      allocate(this%filelist(0))
    end if

    ! Run XIOS setup routines
    call init_xios_dimensions(chi, panel_id)
    call setup_xios_files(this%filelist, this%clock)

    call xios_close_context_definition()

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Clean up an XIOS context object.
  !>
  subroutine finalise( this )

    implicit none

    type(lfric_xios_context_type), intent(inout) :: this

    integer(i_native) :: i

    call xios_context_finalize()

    ! We have closed the context on our end, but we need to make sure that XIOS
    ! has closed the files for all servers before we process them.
    call init_wait()

    do i = 1, size(this%filelist)
      call this%filelist(i)%file_close()
    end do

  end subroutine finalise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets the clock associated with this context.
  !>
  !> @return Clock object.
  !>
  function get_clock( this ) result(clock)

    implicit none

    class(lfric_xios_context_type), intent(in), target :: this
    class(clock_type), pointer :: clock

    clock => this%clock

  end function get_clock

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets the context's file list
  !>
  !> @return  The list of XIOS file type objects used by the model
  !>
  function get_filelist( this ) result( filelist )

    implicit none

    class(lfric_xios_context_type), target, intent(in) :: this
    class(file_type), pointer :: filelist(:)

    filelist => this%filelist

  end function get_filelist

end module lfric_xios_context_mod
