!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the GungHo model.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gungho_driver_mod

  use cli_mod,                    only : get_initial_filename
  use clock_mod,                  only : clock_type
  use constants_mod,              only : i_def, i_native, imdi
  use driver_io_mod,              only : get_clock, get_io_context
  use field_mod,                  only : field_type
  use gungho_mod,                 only : program_name
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model, &
                                         finalise_infrastructure, &
                                         finalise_model
  use gungho_model_data_mod,      only : model_data_type, &
                                         create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use gungho_step_mod,            only : gungho_step
  use initial_output_mod,         only : write_initial_output
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use io_context_mod,             only : io_context_type
  use initialization_config_mod,  only : lbc_option,               &
                                         lbc_option_gungho_file,   &
                                         lbc_option_um2lfric_file, &
                                         ancil_option,             &
                                         ancil_option_updating
  use init_gungho_lbcs_alg_mod,   only : update_lbcs_file_alg
  use log_mod,                    only : log_event,           &
                                         log_level_always,    &
                                         log_level_info,      &
                                         log_scratch_space
  use mesh_mod,                   only : mesh_type
  use derived_config_mod,         only : l_esm_couple
#ifdef UM_PHYSICS
  use variable_fields_mod,        only : update_variable_fields
  use update_ancils_alg_mod,      only : update_ancils_alg
#endif
#ifdef COUPLED
  use esm_couple_config_mod,      only : l_esm_couple_test
  use coupler_mod,                only : cpl_snd, cpl_rcv, cpl_fld_update
  use coupler_diagnostics_mod,    only : save_sea_ice_frac_previous
#endif

  implicit none

  private
  public initialise, run, finalise

  ! Model run working data set
  type (model_data_type) :: model_data

  type(mesh_type), pointer :: mesh              => null()
  type(mesh_type), pointer :: twod_mesh         => null()
  type(mesh_type), pointer :: shifted_mesh      => null()
  type(mesh_type), pointer :: double_level_mesh => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  subroutine initialise()

    implicit none

    class(io_context_type), pointer :: io_context => null()

    character(:), allocatable :: filename

    call get_initial_filename( filename )

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( filename,             &
                                    program_name,         &
                                    mesh,                 &
                                    twod_mesh,            &
                                    shifted_mesh,         &
                                    double_level_mesh,    &
                                    model_data            )

    ! Instantiate the fields stored in model_data
    call create_model_data( model_data, mesh, twod_mesh )

    ! Initialise the fields stored in the model_data
    call initialise_model_data( model_data )

    ! Initial output
    io_context => get_io_context()
    call write_initial_output( mesh, twod_mesh, model_data, &
                               io_context, nodal_output_on_w3 )

    ! Model configuration initialisation
    call initialise_model( mesh, model_data )

#ifdef COUPLED
    ! Placeholder for ESM coupling initialisation code.
    ! Check we have a value for related namelist control variable
    write(log_scratch_space,'(A,L1)') program_name//': Couple flag l_esm_couple_test: ', &
                                     l_esm_couple_test
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
#endif


  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm
  !>       based upon the configuration.
  !>
  subroutine run()

    implicit none

    class(clock_type), pointer :: clock => null()

    write(log_scratch_space,'(A)') 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    if(l_esm_couple) then
      write(log_scratch_space,'(A)') 'Configuration is coupled to ocean'
      call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )
    else
      write(log_scratch_space,'(A)') 'Configuration is not coupled to ocean'
      call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )
    endif

    clock => get_clock()
    do while (clock%tick())

#ifdef COUPLED
      if(l_esm_couple) then

         write(log_scratch_space, *) 'Coupling timestep: ', clock%get_step() - 1
         call log_event( log_scratch_space, LOG_LEVEL_INFO )

         call save_sea_ice_frac_previous(model_data%depository)

         ! Receive all incoming (ocean/seaice fields) from the coupler
         call cpl_rcv(model_data%cpl_rcv, model_data%depository, clock)

         ! Send all outgoing (ocean/seaice driving fields) to the coupler
         call cpl_snd(model_data%cpl_snd, model_data%depository, clock)

      endif
#endif

      if ( lbc_option == lbc_option_gungho_file .or. &
           lbc_option == lbc_option_um2lfric_file) then

        call update_lbcs_file_alg( model_data%lbc_times_list, &
                                   clock, model_data%lbc_fields )
      endif

      ! Perform a timestep
      call gungho_step( mesh, twod_mesh, model_data, &
                        clock )

      ! Use diagnostic output frequency to determine whether to write
      ! diagnostics on this timestep

      if ( ( mod(clock%get_step(), diagnostic_frequency) == 0 ) &
           .and. ( write_diag ) ) then

        ! Calculation and output diagnostics
        call gungho_diagnostics_driver( mesh,       &
                                        twod_mesh,  &
                                        model_data, &
                                        clock,      &
                                        nodal_output_on_w3 )
      end if

#ifdef UM_PHYSICS
      ! Update time-varying ancillaries
      ! This is done last in the timestep, because the time data of the
      ! ancillaries needs to be that of the start of timestep, but the
      ! first thing that happens in a timestep is that the clock ticks to the
      ! end of timestep date.
      if (ancil_option == ancil_option_updating) then
        call update_variable_fields( model_data%ancil_times_list, &
                                     clock, model_data%ancil_fields )
        call update_ancils_alg( model_data%ancil_times_list, &
                                clock, model_data%ancil_fields, &
                                model_data%surface_fields)
      end if
#endif

    end do ! end ts loop

#ifdef COUPLED
    if (l_esm_couple) then
       ! Ensure coupling fields are updated at the end of a cycle to ensure the values
       ! stored in and recovered from checkpoint dumps are correct and reproducible
       ! when (re)starting subsequent runs!
       call cpl_fld_update(model_data%cpl_snd, model_data%depository, clock)
    endif
#endif

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Output the fields stored in the model_data (checkpoint and dump)
    call output_model_data( model_data )

    ! Model configuration finalisation
    call finalise_model( model_data, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module gungho_driver_mod
