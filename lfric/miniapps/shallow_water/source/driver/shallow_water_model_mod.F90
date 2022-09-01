!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!!        and the shallow_water model simulation.
module shallow_water_model_mod

  use assign_orography_field_mod,     only: assign_orography_field
  use checksum_alg_mod,               only: checksum_alg
  use cli_mod,                        only: get_initial_filename
  use clock_mod,                      only: clock_type
  use configuration_mod,              only: final_configuration
  use conservation_algorithm_mod,     only: conservation_algorithm
  use constants_mod,                  only: i_def, i_native, &
                                            PRECISION_REAL
  use convert_to_upper_mod,           only: convert_to_upper
  use count_mod,                      only: count_type, halo_calls
  use derived_config_mod,             only: set_derived_config
  use driver_log_mod,                 only: init_logger, final_logger
  use driver_comm_mod,                only: init_comm, final_comm
  use driver_fem_mod,                 only: init_fem, final_fem
  use driver_io_mod,                  only: init_io, final_io, get_clock, &
                                            filelist_populator
  use driver_mesh_mod,                only: init_mesh, final_mesh
  use field_mod,                      only: field_type
  use field_parent_mod,               only: write_interface
  use field_collection_mod,           only: field_collection_type
  use global_mesh_collection_mod,     only: global_mesh_collection, &
                                            global_mesh_collection_type
  use io_config_mod,                  only: subroutine_timers,       &
                                            subroutine_counters,     &
                                            use_xios_io,             &
                                            write_conservation_diag, &
                                            write_dump,              &
                                            write_minmax_tseries
  use lfric_xios_file_mod,            only: lfric_xios_file_type
  use linked_list_mod,                only: linked_list_type
  use local_mesh_collection_mod,      only: local_mesh_collection, &
                                            local_mesh_collection_type
  use log_mod,                        only: log_event,          &
                                            log_set_level,      &
                                            log_scratch_space,  &
                                            LOG_LEVEL_INFO
  use mesh_collection_mod,            only: mesh_collection, &
                                            mesh_collection_type
  use mesh_mod,                       only: mesh_type
  use minmax_tseries_mod,             only: minmax_tseries,      &
                                            minmax_tseries_init, &
                                            minmax_tseries_final
  use mpi_mod,                        only: get_comm_size, &
                                            get_comm_rank
  use runtime_constants_mod,          only: create_runtime_constants, &
                                            final_runtime_constants
  use shallow_water_mod,              only: load_configuration
  use shallow_water_model_data_mod,   only: model_data_type
  use shallow_water_setup_io_mod,     only: init_shallow_water_files
  use timer_mod,                      only: timer, output_timer, init_timer
  use timestepping_config_mod,        only: dt, &
                                            spinup_period
  use xios,                           only: xios_update_calendar

  implicit none

  private
  public :: initialise_infrastructure, &
            initialise_model,          &
            finalise_infrastructure,   &
            finalise_model

  contains

  !=============================================================================
  !> @brief Initialises the infrastructure and sets up constants used by the model.
  !> @param[in]     program_name An identifier given to the model begin run
  !> @param[in,out] mesh         The 3D mesh
  !> @param[in,out] chi          A size 3 array of fields holding the coordinates of the mesh
  subroutine initialise_infrastructure(program_name, &
                                       mesh,         &
                                       chi           )

    implicit none

    character(*),      intent(in)                    :: program_name
    type(mesh_type),   intent(inout), pointer        :: mesh
    type(field_type),  intent(inout)                 :: chi(3)

    type(field_type), target :: panel_id
    type(mesh_type), pointer :: twod_mesh => null()

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    character(len=*), parameter :: io_context_name = "shallow_water"
    character(:),   allocatable :: filename

    integer(i_native) :: communicator

    class(clock_type), pointer :: clock

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------

    ! Set up the MPI communicator for later use
    call init_comm( program_name, communicator )

    call get_initial_filename( filename )
    call load_configuration( filename )

    call init_logger( get_comm_size(), get_comm_rank(), program_name)

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Initialise timers and counters
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) then
      call init_timer()
      call timer(program_name)
    end if

    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_INFO )

    if ( subroutine_counters ) then
      allocate(halo_calls, source=count_type('halo_calls'))
      call halo_calls%counter(program_name)
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------

    ! TODO Stencil depth needs to be taken from configuration options
    ! Create the mesh
    call init_mesh( get_comm_rank(),  get_comm_size(), mesh, twod_mesh, &
                    input_stencil_depth = 2_i_def )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem(mesh, chi, panel_id)

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
    ! domain and context

    files_init_ptr => init_shallow_water_files
    call init_io( io_context_name, communicator, &
                  chi, panel_id,                 &
                  populate_filelist=files_init_ptr )

    clock => get_clock()

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh, twod_mesh, chi, panel_id, dt)

  end subroutine initialise_infrastructure

  !=============================================================================
  !> @brief Initialises the shallow_water application.
  !> @param[in]     mesh       The primary mesh
  !> @param[in,out] model_data The working data set for the model run
  subroutine initialise_model( mesh,    &
                               model_data )

    use swe_timestep_alg_mod, only: swe_timestep_alg_init

    implicit none

    type(mesh_type),       pointer, intent(in)    :: mesh
    type(model_data_type), target,  intent(inout) :: model_data

    type(field_collection_type), pointer :: prognostic_fields => null()
    type(field_collection_type), pointer :: diagnostic_fields => null()

    ! Prognostic fields
    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: buoyancy => null()
    type(field_type), pointer :: q => null()

    call log_event( 'shallow_water: Initialising miniapp ...', LOG_LEVEL_INFO )

    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields

    wind     => prognostic_fields%get_field("wind")
    buoyancy => prognostic_fields%get_field("buoyancy")
    geopot   => prognostic_fields%get_field("geopot")
    q        => prognostic_fields%get_field("q")

    ! Initialise transport and shallow water model
    call swe_timestep_alg_init( mesh, wind, geopot, buoyancy, q )

    call log_event( 'shallow_water: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine initialise_model


  !=============================================================================
  !> @brief Finalises infrastructure and constants used by the model.
  !> @param[in] program_name      The program name
  subroutine finalise_infrastructure(program_name)

    implicit none

    character(*), intent(in) :: program_name

    !-------------------------------------------------------------------------
    ! Finalise constants
    !-------------------------------------------------------------------------

    call final_runtime_constants()

    !-------------------------------------------------------------------------
    ! Finalise timers and counters
    !-------------------------------------------------------------------------

    if ( subroutine_timers ) then
      call timer(program_name)
      call output_timer()
    end if

    if ( subroutine_counters ) then
      call halo_calls%counter(program_name)
      call halo_calls%output_counters()
    end if

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------

    call final_mesh()
    call final_fem()

    !-------------------------------------------------------------------------
    ! Final logging before infrastructure is destroyed
    !-------------------------------------------------------------------------
    call final_logger( program_name )

    !-------------------------------------------------------------------------
    ! Finalise infrastructure
    !-------------------------------------------------------------------------

    ! Finalise namelist configurations
    call final_configuration()

    call final_comm()

  end subroutine finalise_infrastructure

  !=============================================================================
  !> @brief Finalise the shallow_water application.
  !> @param[in]     mesh_id      The identifier of the primary mesh
  !> @param[in,out] model_data   The working data set for the model run
  !> @param[in]     program_name An identifier given to the model begin run
  subroutine finalise_model( mesh_id,    &
                             model_data, &
                             program_name )

    use swe_timestep_alg_mod, only: swe_timestep_alg_final

    implicit none

    integer(i_def),                intent(in)    :: mesh_id
    type(model_data_type), target, intent(inout) :: model_data
    character(*),                  intent(in)    :: program_name

    type( field_collection_type ), pointer :: prognostic_fields => null()

    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: q => null()

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    prognostic_fields => model_data%prognostic_fields

    wind   => prognostic_fields%get_field('wind')
    geopot => prognostic_fields%get_field('geopot')
    q      => prognostic_fields%get_field('q')

    ! Checksums
    call checksum_alg( program_name, wind, 'wind', geopot, 'geopot', q, 'q' )

    ! Finalise transport
    call swe_timestep_alg_final()

  end subroutine finalise_model

end module shallow_water_model_mod
