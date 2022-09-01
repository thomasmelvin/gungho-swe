!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Create the prognostic fields and place them in the depository.
!!
!> @details Create the prognostic fields and place them both in the
!!          depository field collection and put pointers to them in the
!!          prognostic field collection. Also create surface geopotential
!!          and put it into depository field collection.

module create_shallow_water_prognostics_mod

  use constants_mod,                      only: i_def
  use field_mod,                          only: field_type
  use field_parent_mod,                   only: write_interface,            &
                                                checkpoint_write_interface, &
                                                checkpoint_read_interface
  use field_collection_mod,               only: field_collection_type
  use finite_element_config_mod,          only: element_order
  use function_space_collection_mod,      only: function_space_collection
  use fs_continuity_mod,                  only: W1, W2, W3
  use io_config_mod,                      only: write_diag,       &
                                                use_xios_io,      &
                                                checkpoint_write, &
                                                checkpoint_read
  use pure_abstract_field_mod,            only: pure_abstract_field_type
  use lfric_xios_write_mod,               only: write_field_face, &
                                                checkpoint_write_xios
  use lfric_xios_read_mod,                only: checkpoint_read_xios
  use io_mod,                             only: checkpoint_read_netcdf, &
                                                checkpoint_write_netcdf
  use log_mod,                            only: log_event,      &
                                                LOG_LEVEL_INFO, &
                                                LOG_LEVEL_ERROR
  use mesh_mod,                           only: mesh_type

  implicit none

  private

  public :: create_shallow_water_prognostics

  contains

  !> @brief Create the prognostic fields and place them in the depository.
  !> @details Creates the prognostic fields for the shallow water miniapp.
  !!          The prognostics are all the fields in the depository - so the
  !!          two collections contain the same fields.
  !> @param [in]     mesh        Mesh to initialise variables on
  !> @param [in,out] depository  A collection of all fields used in the miniapp
  !> @param [in,out] prognostics A collection of all the prognostic fields.
  !> @param [in,out] s_geopot    Surface geopotential.
  subroutine create_shallow_water_prognostics( mesh,        &
                                               depository,  &
                                               prognostics, &
                                               s_geopot)

    implicit none

    type(mesh_type), pointer,    intent(in)    :: mesh
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostics
    type(field_type),            intent(inout) :: s_geopot

    type(field_type) :: wind
    type(field_type) :: buoyancy
    type(field_type) :: geopot
    type(field_type) :: q

    class(pure_abstract_field_type), pointer  :: tmp_ptr => null()

    integer(i_def) :: buoyancy_space
    integer(i_def) :: vorticity_space

    procedure(write_interface),            pointer :: tmp_write_ptr
    procedure(checkpoint_write_interface), pointer :: tmp_checkpoint_write_ptr => null()
    procedure(checkpoint_read_interface),  pointer :: tmp_checkpoint_read_ptr => null()



    call depository%initialise(name='depository', table_len=100)
    call prognostics%initialise(name="prognostics", table_len=100)

    ! Create prognostic fields
    buoyancy_space = W3
    call log_event( 'shallow_water: Using V2 for buoyancy', LOG_LEVEL_INFO )
    vorticity_space = W3
    call log_event( 'shallow_water: Using V2 for vorticity', LOG_LEVEL_INFO )

    call wind%initialise( vector_space =                                                &
                          function_space_collection%get_fs(mesh, element_order, W2), &
                          name="wind" )
    call buoyancy%initialise( vector_space =                                                            &
                              function_space_collection%get_fs(mesh, element_order, buoyancy_space), &
                              name="buoyancy" )
    call geopot%initialise( vector_space =                                                &
                            function_space_collection%get_fs(mesh, element_order, W3), &
                            name="geopot" )
    call q%initialise( vector_space =                                                             &
                       function_space_collection%get_fs(mesh, element_order, vorticity_space), &
                       name="q")
    call s_geopot%initialise( vector_space =                                                &
                              function_space_collection%get_fs(mesh, element_order, W3), &
                              name="s_geopot" )

    ! Set I/O behaviours for diagnostic output
    if (write_diag .and. use_xios_io) then
       ! Fields that are output on the XIOS face domain
       tmp_write_ptr => write_field_face
       call wind%set_write_behaviour(tmp_write_ptr)
       call geopot%set_write_behaviour(tmp_write_ptr)
       call buoyancy%set_write_behaviour(tmp_write_ptr)
       call q%set_write_behaviour(tmp_write_ptr)
    end if

    ! Set I/O behaviours for checkpoint / restart
    if ( checkpoint_write .or. checkpoint_read) then
      if (use_xios_io) then
        ! Use XIOS for checkpoint / restart
        tmp_checkpoint_write_ptr => checkpoint_write_xios
        tmp_checkpoint_read_ptr => checkpoint_read_xios
        call log_event( 'GungHo: Using XIOS for checkpointing...', LOG_LEVEL_INFO )
      else
        ! Use old checkpoint and restart methods
        tmp_checkpoint_write_ptr => checkpoint_write_netcdf
        tmp_checkpoint_read_ptr => checkpoint_read_netcdf
        call log_event( 'GungHo: Using NetCDF for checkpointing...', LOG_LEVEL_INFO )
      end if

      call wind%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call geopot%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call buoyancy%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call q%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call s_geopot%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)

      call wind%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call geopot%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call buoyancy%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call q%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call s_geopot%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)

    end if

    ! Put the prognostic fields into the depository
    call depository%add_field(wind)
    call depository%add_field(geopot)
    call depository%add_field(buoyancy)
    call depository%add_field(q)
    call depository%add_field(s_geopot)

    ! Put pointers to the prognostic fields into the prognostic collection
    tmp_ptr => depository%get_field('wind')
    call prognostics%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('geopot')
    call prognostics%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('buoyancy')
    call prognostics%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('q')
    call prognostics%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('s_geopot')
    call prognostics%add_reference_to_field(tmp_ptr)

  end subroutine create_shallow_water_prognostics

end module create_shallow_water_prognostics_mod
