!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Creates finite difference (fd) prognostic fields
!> @details Handles creation of finite difference fields

module create_fd_prognostics_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use field_parent_mod,               only : read_interface, &
                                             write_interface
  use lfric_xios_read_mod,            only : read_field_face, &
                                             read_field_edge
  use lfric_xios_write_mod,           only : write_field_face, &
                                             write_field_edge
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W3, Wtheta, W2H
  use log_mod,                        only : log_event,      &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use init_ancils_mod,                only : setup_ancil_field
  use initialization_config_mod,      only : ancil_option,                &
                                             ancil_option_start_dump,     &
                                             ancil_option_fixed,          &
                                             ancil_option_updating,       &
                                             read_w2h_wind
  use nlsizes_namelist_mod,           only : sm_levels
  use jules_control_init_mod,         only : n_land_tile, n_sea_ice_tile
  use jules_physics_init_mod,         only : snow_lev_tile
  use derived_config_mod,             only : l_esm_couple

  implicit none
  private
  public :: create_fd_prognostics

contains
  !>@brief Routine to create a field collection for finite difference
  !>       prognostic fields
  !> @param[in] mesh        The mesh
  !> @param[in] twod_mesh   The current 2d mesh
  !> @param[in,out] fd_field_collection The collection object to store fields in
  !> @param[in,out] depository The depository field collection
  subroutine create_fd_prognostics( mesh, twod_mesh, fd_field_collection, &
                                    depository)

    implicit none

    type( mesh_type ), intent(in), pointer     :: mesh
    type( mesh_type ), intent(in), pointer     :: twod_mesh

    type(field_collection_type), intent(inout) :: fd_field_collection
    type(field_collection_type), intent(inout) :: depository

    procedure(read_interface), pointer  :: tmp_read_ptr => null()
    procedure(write_interface), pointer  :: tmp_write_ptr => null()

    ! FD field declarations
    type( field_type ) :: ew_wind_in_w3 ! U wind
    type( field_type ) :: ns_wind_in_w3 ! V wind
    type( field_type ) :: h_wind_in_w2h ! Horizontal wind (i.e. on W2h dofs)
    type( field_type ) :: dry_rho_in_w3 ! Dry rho
    ! Vertical theta levels
    type( field_type ) :: upward_wind_in_wtheta ! W wind
    type( field_type ) :: theta_in_wtheta ! Potential temp
    type( field_type ) :: mv_in_wtheta    ! Vapour mix ratio
    type( field_type ) :: mcl_in_wtheta   ! Cloud liquid mix ratio
    type( field_type ) :: mcf_in_wtheta   ! Clould ice mix ratio
    type( field_type ) :: mr_in_wtheta    ! Rain mix ratio

    call log_event( 'Physics: Creating Finite Difference prognostics...', LOG_LEVEL_INFO )

    if (element_order > 0) then
      call log_event( 'Finite diff prognostics: requires lowest order elements'&
           , LOG_LEVEL_ERROR )
    end if

    ! Create the field collection
    call fd_field_collection%initialise(name="fd_prognostics", table_len=100)

    if ( read_w2h_wind )then
       ! In this case we read in directly onto the W2H dofs
       tmp_read_ptr => read_field_edge
       tmp_write_ptr => write_field_edge
       call h_wind_in_w2h%initialise( vector_space = &
         function_space_collection%get_fs(mesh, element_order, W2H), &
         name='h_wind')
       call h_wind_in_w2h%set_read_behaviour(tmp_read_ptr)
       call h_wind_in_w2h%set_write_behaviour(tmp_write_ptr)
       call fd_field_collection%add_field(h_wind_in_w2h)

    else

      ! Setup I/O behaviour handler. In the case of FD prognostic fields these
      ! are currently read from a UM2LFRic dump
      tmp_read_ptr => read_field_face
      tmp_write_ptr => write_field_face

      ! Create the fields, set the I/O behaviour and add to
      ! the field collection
      !========================================================================
      ! W3 fields - rho levels
      !========================================================================

      call ew_wind_in_w3%initialise( vector_space = &
        function_space_collection%get_fs(mesh, element_order, W3), &
        name='ew_wind_in_w3')

      call ew_wind_in_w3%set_read_behaviour(tmp_read_ptr)
      call ew_wind_in_w3%set_write_behaviour(tmp_write_ptr)

      call fd_field_collection%add_field(ew_wind_in_w3)

      call ns_wind_in_w3%initialise( vector_space = &
        function_space_collection%get_fs(mesh, element_order, W3), &
        name='ns_wind_in_w3')

      call ns_wind_in_w3%set_read_behaviour(tmp_read_ptr)
      call ns_wind_in_w3%set_write_behaviour(tmp_write_ptr)

      call fd_field_collection%add_field(ns_wind_in_w3)

    end if

    tmp_read_ptr => read_field_face
    tmp_write_ptr => write_field_face

    call dry_rho_in_w3%initialise( vector_space = &
         function_space_collection%get_fs(mesh, element_order, W3), &
         name='dry_rho_in_w3')

    call dry_rho_in_w3%set_read_behaviour(tmp_read_ptr)
    call dry_rho_in_w3%set_write_behaviour(tmp_write_ptr)

    call fd_field_collection%add_field(dry_rho_in_w3)

    !========================================================================
    ! Wtheta fields - theta levels
    !========================================================================

    call upward_wind_in_wtheta%initialise( vector_space = &
         function_space_collection%get_fs(mesh, element_order, Wtheta), &
         name='upward_wind_in_wtheta')

    call upward_wind_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call upward_wind_in_wtheta%set_write_behaviour(tmp_write_ptr)

    call fd_field_collection%add_field(upward_wind_in_wtheta)

    call theta_in_wtheta%initialise( vector_space =        &
         function_space_collection%get_fs(mesh, element_order, Wtheta), &
         name='theta_in_wtheta')

    call theta_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call theta_in_wtheta%set_write_behaviour(tmp_write_ptr)

    call fd_field_collection%add_field(theta_in_wtheta)

    call mv_in_wtheta%initialise( vector_space =           &
         function_space_collection%get_fs(mesh, element_order, Wtheta), &
         name='mv_in_wtheta')

    call mv_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mv_in_wtheta%set_write_behaviour(tmp_write_ptr)

    call fd_field_collection%add_field(mv_in_wtheta)

    call mcl_in_wtheta%initialise( vector_space =          &
         function_space_collection%get_fs(mesh, element_order, Wtheta), &
         name='mcl_in_wtheta')

    call mcl_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mcl_in_wtheta%set_write_behaviour(tmp_write_ptr)

    call fd_field_collection%add_field(mcl_in_wtheta)

    call mcf_in_wtheta%initialise( vector_space =          &
         function_space_collection%get_fs(mesh, element_order, Wtheta), &
         name='mcf_in_wtheta')

    call mcf_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mcf_in_wtheta%set_write_behaviour(tmp_write_ptr)

    call fd_field_collection%add_field(mcf_in_wtheta)

    call mr_in_wtheta%initialise( vector_space =          &
         function_space_collection%get_fs(mesh, element_order, Wtheta), &
         name='mr_in_wtheta')

    call mr_in_wtheta%set_read_behaviour(tmp_read_ptr)
    call mr_in_wtheta%set_write_behaviour(tmp_write_ptr)

    call fd_field_collection%add_field(mr_in_wtheta)

    !========================================================================
    ! Physics fields
    !========================================================================
    ! turbulence fields
    call setup_ancil_field("zh", depository, fd_field_collection, mesh, &
                           twod_mesh, twod=.true.)
    ! cloud fields
    call setup_ancil_field("area_fraction", depository, &
                           fd_field_collection, mesh, twod_mesh)
    call setup_ancil_field("bulk_fraction", depository, &
                           fd_field_collection, mesh, twod_mesh)
    call setup_ancil_field("liquid_fraction", depository, &
                           fd_field_collection, mesh, twod_mesh)
    call setup_ancil_field("frozen_fraction", depository, &
                           fd_field_collection, mesh, twod_mesh)
    ! surface fields
    call setup_ancil_field("z0msea", depository, fd_field_collection, &
                           mesh, twod_mesh, twod=.true.)

    if (ancil_option == ancil_option_start_dump) then
      call setup_ancil_field("tstar", depository, fd_field_collection, &
                             mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("ozone", depository, fd_field_collection, &
                             mesh, twod_mesh)
    else if (ancil_option == ancil_option_fixed .or. &
             ancil_option == ancil_option_updating) then
      ! convection fields
      call setup_ancil_field("dd_mf_cb", depository, fd_field_collection, &
                             mesh, twod_mesh, twod=.true.)
      ! soil fields
      call setup_ancil_field("soil_sat_frac", depository, fd_field_collection, &
                             mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("water_table", depository, fd_field_collection, &
                             mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("wetness_under_soil", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true.)
      call setup_ancil_field("soil_temperature", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=sm_levels)
      call setup_ancil_field("soil_moisture", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=sm_levels)
      ! surface fields
      call setup_ancil_field("surface_conductance", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true.)
      call setup_ancil_field("can_water_in", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("land_tile_temp", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("tstar_sea_ice", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_sea_ice_tile)
      call setup_ancil_field("sea_ice_temperature", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_sea_ice_tile)

      ! For coupled models get the sea ice fraction and thickness from the
      ! dump
      if (l_esm_couple) then
         call setup_ancil_field("sea_ice_fraction", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_sea_ice_tile)
         call setup_ancil_field("sea_ice_thickness", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_sea_ice_tile)
      endif

      ! snow fields
      call setup_ancil_field("tile_snow_mass_in", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("n_snow_layers_in", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("snow_depth_in", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("tile_snow_rgrain_in", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("snow_soot", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true.)
      call setup_ancil_field("snow_under_canopy_in", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("snowpack_density_in", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=n_land_tile)
      call setup_ancil_field("snow_layer_thickness", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=snow_lev_tile)
      call setup_ancil_field("snow_layer_ice_mass", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=snow_lev_tile)
      call setup_ancil_field("snow_layer_liq_mass", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=snow_lev_tile)
      call setup_ancil_field("snow_layer_temp", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=snow_lev_tile)
      call setup_ancil_field("snow_layer_rgrain", depository, &
                             fd_field_collection, mesh, twod_mesh, &
                             twod=.true., ndata=snow_lev_tile)
    end if

    call log_event( 'Physics: Finite diff prognostics created', LOG_LEVEL_INFO )

  end subroutine create_fd_prognostics

end module create_fd_prognostics_mod
