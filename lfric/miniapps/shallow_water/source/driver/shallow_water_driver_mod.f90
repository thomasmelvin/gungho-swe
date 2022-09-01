!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Drives the execution of the shallow_water miniapp.
!!
module shallow_water_driver_mod

  use clock_mod,                     only: clock_type
  use constants_mod,                 only: i_def, i_native, imdi, r_def
  use driver_io_mod,                 only: get_clock, get_io_context
  use field_mod,                     only: field_type
  use field_collection_mod,          only: field_collection_type
  use io_config_mod,                 only: write_diag,           &
                                           diagnostic_frequency, &
                                           nodal_output_on_w3
  use io_context_mod,                only: io_context_type
  use log_mod,                       only: log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_INFO
  use mesh_collection_mod,           only: mesh_collection, &
                                           mesh_collection_type
  use mesh_mod,                      only: mesh_type
  use shallow_water_mod,             only: program_name
  use shallow_water_diagnostics_mod, only: shallow_water_diagnostics
  use shallow_water_model_mod,       only: initialise_infrastructure, &
                                           initialise_model,          &
                                           finalise_infrastructure,   &
                                           finalise_model
  use shallow_water_model_data_mod,  only: model_data_type,       &
                                           create_model_data,     &
                                           initialise_model_data, &
                                           output_model_data,     &
                                           finalise_model_data
  use shallow_water_step_mod,        only: shallow_water_step

  implicit none

  private

  public :: initialise, &
            run,        &
            finalise

  ! Model run working data set
  type (model_data_type) :: model_data

  ! Coordinate field
  type(field_type), target :: chi(3)

  ! Mesh
  integer(i_def) :: mesh_id = imdi
  type(mesh_type), pointer :: mesh => null()

  class(io_context_type), pointer :: io_context => null()

contains

  !=============================================================================
  !> @brief   Sets up required state in preparation for run.
  !> @details Initialises the infrastructure and the fields stored in
  !!          model_data, then sets the initial conditions for the run.
  !> @param[in] filename Configuration namelist file
  !!
  subroutine initialise()

    implicit none

    class(clock_type), pointer :: clock => null()

    call log_event( 'Initialising Infrastructure ...', LOG_LEVEL_INFO )
    ! Initialise infrastructure (from shallow_water_model_mod.F90) and setup constants
    call initialise_infrastructure( program_name,       &
                                    mesh,               &
                                    chi )

    io_context => get_io_context()
    clock => get_clock()

    call log_event( 'Creating model data ...', LOG_LEVEL_INFO )
    ! Instantiate the fields stored in model_data
    call create_model_data( model_data,   &
                            mesh )

    call log_event( 'Initialising model data ...', LOG_LEVEL_INFO )
    ! Initialise the fields stored in the model_data
    call initialise_model_data( model_data, mesh, clock )

    ! Initial output: we only want these once at the beginning of a run
    if (clock%is_initialisation() .and. write_diag) then
      call log_event( 'Output of initial diagnostics ...', LOG_LEVEL_INFO )
      ! Calculation and output of diagnostics
      call shallow_water_diagnostics( mesh,               &
                                      model_data,         &
                                      clock,              &
                                      nodal_output_on_w3, &
                                      1_i_def )
    end if

    ! Model configuration initialisation
    call initialise_model( mesh,    &
                           model_data )

  end subroutine initialise

  !=============================================================================
  !> @brief Performs time stepping for the shallow water miniapp.
  !!
  subroutine run()

    implicit none

    class(clock_type), pointer :: clock => null()

    clock => get_clock()

    write(log_scratch_space,'(A)') 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    !--------------------------------------------------------------------------
    ! Model step
    !--------------------------------------------------------------------------
    do while (clock%tick())

      call shallow_water_step( model_data, &
                               clock  )

      ! Use diagnostic output frequency to determine whether to write
      ! diagnostics on this timestep

      if ( ( mod(clock%get_step(), diagnostic_frequency) == 0 ) &
           .and. ( write_diag ) ) then

        call shallow_water_diagnostics(mesh,               &
                                       model_data,         &
                                       clock,              &
                                       nodal_output_on_w3, &
                                       0_i_def)


      end if

    end do ! while clock%is_running()

  end subroutine run

  !=============================================================================
  !> @brief Tidies up after a run.
  !!
  subroutine finalise()

    implicit none

    class(clock_type), pointer :: clock => null()

    clock => get_clock()

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_INFO )

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( model_data, clock )

    ! Model configuration finalisation
    call finalise_model( mesh_id,    &
                         model_data, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module shallow_water_driver_mod
