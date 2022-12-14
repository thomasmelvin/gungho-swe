!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> A timestepping clock implementation.
!>
module model_clock_mod

  use calendar_mod,  only : calendar_type
  use clock_mod,     only : clock_type
  use constants_mod, only : i_timestep, r_def, r_second, l_def
  use log_mod,       only : log_event, log_level_error, log_scratch_space, &
                            log_set_timestep, log_forget_timestep

  implicit none

  private

  !> A clock which ticks in timesteps.
  !>
  !> The overall run (from step 1 to step @f$n@f$) may be split into a
  !> number of partial runs (subset of overall run). The clock ticks over a
  !> partial run from timestep @f$a@f$ to step @f$b@f$ where @f$a \geq 1@f$,
  !> @f$a < b@f$ and @f$b \leq n@f$.
  !>
  !> The clock exists in a number of states accessed by various "is_*" methods.
  !>
  type, public, extends(clock_type) :: model_clock_type
    private
    class(calendar_type), allocatable :: calendar
    integer(i_timestep) :: current_step
    integer(i_timestep) :: first_step
    integer(i_timestep) :: last_step
    integer(i_timestep) :: last_spinup_step
    real(r_second)      :: seconds_per_step
    real(r_def)         :: spinup_fraction
    logical             :: initialisation_phase
    logical             :: starting
  contains
    private
    procedure, public :: tick
    procedure, public :: get_first_step
    procedure, public :: get_step
    procedure, public :: get_last_step
    procedure, public :: get_seconds_per_step
    procedure, public :: seconds_from_steps
    procedure, public :: get_calendar
    procedure, public :: get_spinup_fraction
    procedure, public :: is_initialisation
    procedure, public :: is_running
    procedure, public :: is_spinning_up
    procedure :: calculate_spinup_fraction
  end type model_clock_type

  interface model_clock_type
    procedure model_clock_constructor
  end interface model_clock_type

contains

  !> @brief Sets up the clock object before use.
  !>
  !> @param[in] calendar Means to interpret human dates.
  !> @param[in] first First date in the current run.
  !> @param[in] last Last date in the current run.
  !> @param[in] seconds_per_step Length of a timestep in seconds.
  !> @param[in] spinup_period Length of spinup period in seconds. May be zero.
  !>
  function model_clock_constructor( calendar,         &
                                    first,            &
                                    last,             &
                                    seconds_per_step, &
                                    spinup_period ) result(new_clock)

    implicit none

    class(calendar_type),     intent(in)    :: calendar
    character(*),             intent(in)    :: first
    character(*),             intent(in)    :: last
    real(r_second),           intent(in)    :: seconds_per_step
    real(r_second),           intent(in)    :: spinup_period
    type(model_clock_type) :: new_clock

    allocate( new_clock%calendar, source=calendar )
    new_clock%first_step = new_clock%calendar%parse_instance( first )
    if (new_clock%first_step < 1) then
      write(log_scratch_space, '("First clock step must be positive")')
      call log_event(log_scratch_space, log_level_error)
    end if
    new_clock%last_step = new_clock%calendar%parse_instance( last )

    if (new_clock%last_step < new_clock%first_step) then
      write(log_scratch_space, '("Last clock step must be after first")')
      call log_event(log_scratch_space, log_level_error)
    end if

    if (seconds_per_step <= 0.0_r_second) then
      call log_event( 'Delta T must be greater than zero.', log_level_error )
    else
      new_clock%seconds_per_step = seconds_per_step
    end if

    new_clock%current_step = new_clock%first_step
    new_clock%last_spinup_step = ceiling( spinup_period / seconds_per_step )

    new_clock%spinup_fraction = new_clock%calculate_spinup_fraction()

    new_clock%initialisation_phase = (new_clock%current_step == 1_i_timestep)
    new_clock%starting = .true.

  end function model_clock_constructor


  !> Gets the calendar the clock is working on.
  !>
  !> @return Calendar pointer should never be unassociated.
  !>
  function get_calendar( this )

    implicit none

    class(model_clock_type), intent(in), target :: this
    class(calendar_type), pointer         :: get_calendar

    get_calendar => this%calendar

  end function get_calendar


  !> Gets the first step in the current run.
  !>
  !> @return Timestep, always greater than zero.
  !>
  function get_first_step( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    integer(i_timestep) :: get_first_step

    get_first_step = this%first_step

  end function get_first_step


  !> Gets the last step in the current run.
  !>
  !> @return Timestep, may be the same as the first step.
  !>
  function get_last_step( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    integer(i_timestep) :: get_last_step

    get_last_step = this%last_step

  end function get_last_step


  !> Gets the length of a timestep.
  !>
  !> @return Timestep length in seconds. Always greater than zero.
  !>
  function get_seconds_per_step( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    real(r_second) :: get_seconds_per_step

    get_seconds_per_step = this%seconds_per_step

  end function get_seconds_per_step


  !> Gets the current spinup period fraction.
  !>
  !> This function returns the value of a ramp from 0.0 at the start of the
  !> first timestep to 1.0 at the end of the last step in the spinup period.
  !>
  !> The actual value returned is the average across the timestep. In the case
  !> of a linear ramp such as we have this is equivalent to taking the value
  !> at the mid-point of the step.
  !>
  !> If this function is called after the spinup period is complete it will
  !> return 1.0.
  !>
  !> @returns Value between 0.0 and 1.0.
  !>
  function get_spinup_fraction( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    real(r_def) :: get_spinup_fraction

    get_spinup_fraction = this%spinup_fraction

  end function get_spinup_fraction


  !> Gets the current timestep.
  !>
  !> @return Timestep between first and last.
  !>
  function get_step( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    integer(i_timestep) :: get_step

    get_step = min(this%current_step, this%last_step)

  end function get_step


  !> Indicates whether the clock is in the "initialisation" state.
  !>
  !> The clock is in the initialisation state for only step 1, the first step
  !> of the overall run. The real value of this state comes from the fact that
  !> the clock is in it from the moment it is initialised and thus can be used
  !> before the first calls to is_running() and tick().
  !>
  !> @return True if clock is in initialisation state.
  !>
  function is_initialisation( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    logical :: is_initialisation

    is_initialisation = this%initialisation_phase

  end function is_initialisation


  !> Indicates whether the clock is in the "running" state.
  !>
  !> The "running" state exists while the clock's current step is between its
  !> first and last step. i.e. it indicates the time period over which the
  !> model should run.
  !>
  !> @return True if clock is in running state.
  !>
  function is_running( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    logical :: is_running

    is_running = (this%current_step <= this%last_step)

  end function is_running


  !> Indicates whether the clock is in the "spinning up" state.
  !>
  !> There is a "spinup" period which starts at step 1 and runs for a user
  !> defined length of time. (possibly no time) This period may span several
  !> partial run periods.
  !>
  !> The duration of the spinup phase is specified in seconds during
  !> initialisation and is rounded up to the nearest whole timestep.
  !>
  !> While in the spin up period calculate_spinup_fraction() may be used.
  !>
  !> @return True if clock is in spinup state.
  !>
  function is_spinning_up( this )

    implicit none

    class(model_clock_type), intent(in) :: this
    logical :: is_spinning_up

    is_spinning_up = (this%current_step <= this%last_spinup_step)

  end function is_spinning_up


  !> Converts a number of timesteps to a number of seconds.
  !>
  !> @return Seconds will be positive or negative depending on period.
  !>
  function seconds_from_steps( this, period )

    implicit none

    class(model_clock_type),   intent(in) :: this
    integer(i_timestep), intent(in) :: period
    real(r_second) :: seconds_from_steps

    seconds_from_steps = period * this%seconds_per_step

  end function seconds_from_steps


  !> Advances the clock by one timestep.
  !>
  !> @return True if clock is still running.
  !>
  function tick( this )

    implicit none

    class(model_clock_type), intent(inout) :: this
    logical :: tick

    if (this%starting) then
      this%starting = .false.
      this%initialisation_phase = .false.
    else
      this%current_step = this%current_step + 1
    end if
    if (this%is_running()) then
      call log_set_timestep( this%current_step )
    else
      call log_forget_timestep()
    end if
    this%spinup_fraction = this%calculate_spinup_fraction()
    tick = this%is_running()

  end function tick


  ! Gets the fraction of the elapsed spinup period.
  !
  ! Returns fraction between 0.0 and 1.0.
  !
  function calculate_spinup_fraction( this ) result(frac)

    implicit none

    class(model_clock_type), intent(inout) :: this
    real(r_def) :: frac

    if (this%last_spinup_step > 0.0) then ! avoids divide by zero
      frac = min( (real(this%current_step, r_def) - 0.5_r_def) &
                  / real(this%last_spinup_step, r_def), &
                  1.0_r_def )
    else
      frac = 1.0_r_def
    end if

  end function calculate_spinup_fraction

end module model_clock_mod
