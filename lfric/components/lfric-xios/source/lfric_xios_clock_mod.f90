!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Integrate the model clock with the XIOS clock.
!>
module lfric_xios_clock_mod

  use calendar_mod,    only : calendar_type
  use model_clock_mod, only : model_clock_type
  use constants_mod,   only : i_timestep, r_second, &
                              l_def
  use timer_mod,       only : timer
  use xios,            only : operator(+),         &
                              xios_date,           &
                              xios_duration,       &
                              xios_get_start_date, &
                              xios_set_start_date, &
                              xios_set_timestep,   &
                              xios_update_calendar

  implicit none

  private

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Integrates the model's clock with XIOS.
  !>
  type, public, extends(model_clock_type) :: lfric_xios_clock_type
    private
    integer :: step_offset
    logical :: uses_timer = .false.
  contains
    private
    procedure, public :: initial_step
    procedure, public :: tick
  end type lfric_xios_clock_type

  interface lfric_xios_clock_type
    procedure lfric_xios_clock_constructor
  end interface lfric_xios_clock_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up an XIOS clock object.
  !>
  !> @param [in] calendar          Interprets human-readable times.
  !> @param [in] first             Time of first step.
  !> @param [in] last              Time of last step.
  !> @param [in] seconds_per_step  Length of a time step in seconds.
  !> @param [in] spinup_period     Number of seconds in spinup period.
  !> @param [in] timer_flag        Flag for use of subroutine timers.
  !>
  function lfric_xios_clock_constructor( calendar,         &
                                         first,            &
                                         last,             &
                                         seconds_per_step, &
                                         spinup_period,    &
                                         timer_flag ) result(new_clock)

    implicit none

    class(calendar_type),         intent(in)    :: calendar
    character(*),                 intent(in)    :: first
    character(*),                 intent(in)    :: last
    real(r_second),               intent(in)    :: seconds_per_step
    real(r_second),               intent(in)    :: spinup_period
    logical(l_def), optional,     intent(in)    :: timer_flag
    type(lfric_xios_clock_type) :: new_clock

    type(xios_duration) :: xios_since_timestep_zero, &
                           timestep_length_for_xios
    type(xios_date)     :: xios_start_date

    if ( present(timer_flag) ) then
      new_clock%uses_timer = timer_flag
    end if

    new_clock%model_clock_type = model_clock_type( calendar,         &
                                                   first,            &
                                                   last,             &
                                                   seconds_per_step, &
                                                   spinup_period )
    new_clock%step_offset = new_clock%get_first_step() - 1

    ! Set the current date by adding the run length so far to the run start date
    ! obtained from XIOS
    call xios_get_start_date(xios_start_date)
    xios_since_timestep_zero%second &
        = new_clock%seconds_from_steps(new_clock%step_offset)
    xios_start_date = xios_start_date + xios_since_timestep_zero
    call xios_set_start_date(xios_start_date)

    ! Set the XIOS time-step from the model clock
    timestep_length_for_xios%second = new_clock%get_seconds_per_step()
    call xios_set_timestep( timestep_length_for_xios )

  end function lfric_xios_clock_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs the first clock step during the intialisation phase - updates
  !> XIOS calendar without ticking model clock
  !>
  subroutine initial_step( this )

    implicit none

    class(lfric_xios_clock_type), intent(inout) :: this

    if ( this%is_initialisation() ) then
      if ( this%uses_timer ) call timer('xios_update_calendar')
      call xios_update_calendar( this%get_step() - this%get_first_step() + 1 )
      if ( this%uses_timer ) call timer('xios_update_calendar')
    end if

  end subroutine initial_step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Advances the clock by one step.
  !>
  function tick( this )

    implicit none

    class(lfric_xios_clock_type), intent(inout) :: this
    logical                                     :: tick

    tick = this%model_clock_type%tick()

    if ( .not. this%is_initialisation() ) then
      if ( this%uses_timer ) call timer('xios_update_calendar')
      call xios_update_calendar( this%get_step() - this%get_first_step() + 1 )
      if ( this%uses_timer ) call timer('xios_update_calendar')
    end if

  end function tick

end module lfric_xios_clock_mod
