!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module init_ancils_mod

  use constants_mod,                  only : i_def, l_def, str_def,   &
                                             r_def
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use field_mod,                      only : field_type
  use field_parent_mod,               only : read_interface, &
                                             write_interface
  use io_config_mod,                  only : use_xios_io
  use linked_list_mod,                only : linked_list_type
  use lfric_xios_read_mod,            only : read_field_face, &
                                             read_field_single_face, &
                                             read_field_time_var
  use lfric_xios_write_mod,           only : write_field_face, &
                                             write_field_single_face
  use field_collection_mod,           only : field_collection_type
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W3, WTheta
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use lfric_xios_time_axis_mod,       only : time_axis_type, update_interface
  use jules_control_init_mod,         only : n_land_tile
  use jules_surface_types_mod,        only : npft
  use dust_parameters_mod,            only : ndiv
  use initialization_config_mod,      only : ancil_option,ancil_option_updating
  use aerosol_config_mod,             only : glomap_mode, glomap_mode_ukca, &
                                             glomap_mode_climatology
  use jules_surface_config_mod,       only : l_vary_z0m_soil
  use surface_config_mod,             only : sea_alb_var_chl, albedo_obs
  use radiation_config_mod,           only : topography, topography_slope, &
                                             topography_horizon, &
                                             n_horiz_ang, n_horiz_layer
  use derived_config_mod,             only : l_esm_couple

  implicit none

  public   :: create_fd_ancils,         &
              setup_ancil_field

contains

  !> @details Organises fields to be read from ancils into ancil_fields
  !           collection then reads them.
  !> @param[in,out] depository The depository field collection
  !> @param[out] ancil_fields Collection for ancillary fields
  !> @param[in] mesh      The current 3d mesh
  !> @param[in] twod_mesh The current 2d mesh
  subroutine create_fd_ancils( depository, ancil_fields, mesh, &
                               twod_mesh, ancil_times_list )

    implicit none

    type( field_collection_type ), intent( inout ) :: depository
    type( field_collection_type ), intent( out )   :: ancil_fields

    type( mesh_type ), intent(in), pointer :: mesh
    type( mesh_type ), intent(in), pointer :: twod_mesh

    type(linked_list_type), intent(out) :: ancil_times_list

    ! Pointer to time-axis update procedure
    procedure(update_interface), pointer :: tmp_update_ptr => null()

    ! Time axis objects for different ancil groups - must be saved to be
    ! available after function call
    type(time_axis_type), save :: sea_time_axis
    type(time_axis_type), save :: sst_time_axis
    type(time_axis_type), save :: sea_ice_time_axis
    type(time_axis_type), save :: aerosol_time_axis
    type(time_axis_type), save :: albedo_vis_time_axis
    type(time_axis_type), save :: albedo_nir_time_axis
    type(time_axis_type), save :: pft_time_axis
    type(time_axis_type), save :: ozone_time_axis
    type(time_axis_type), save :: em_bc_bf_time_axis
    type(time_axis_type), save :: em_bc_ff_time_axis
    type(time_axis_type), save :: em_bc_bb_time_axis
    type(time_axis_type), save :: em_dms_lnd_time_axis
    type(time_axis_type), save :: dms_ocn_time_axis
    type(time_axis_type), save :: em_mterp_time_axis
    type(time_axis_type), save :: em_om_bf_time_axis
    type(time_axis_type), save :: em_om_ff_time_axis
    type(time_axis_type), save :: em_om_bb_time_axis
    type(time_axis_type), save :: em_so2_lo_time_axis
    type(time_axis_type), save :: em_so2_hi_time_axis
    type(time_axis_type), save :: h2o2_limit_time_axis
    type(time_axis_type), save :: ho2_time_axis
    type(time_axis_type), save :: no3_time_axis
    type(time_axis_type), save :: o3_time_axis
    type(time_axis_type), save :: oh_time_axis

    ! Time axis options
    logical(l_def),   parameter :: interp_flag=.true.

    ! Set pointer to time axis read behaviour
    tmp_update_ptr => read_field_time_var

    ! Set up ancil_fields collection
    write(log_scratch_space,'(A,A)') "Create ancil fields: "// &
          "Setting up ancil field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    call ancil_fields%initialise(name='ancil_fields', table_len=100)

    ! Here ancil fields are set up with a call to setup_ancil_field. For ancils
    ! that are time-varying, the time-axis is passed to the setup_ancil_field
    ! subroutine.

    !=====  LAND ANCILS  =====
    call setup_ancil_field("land_area_fraction", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("land_tile_fraction", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.,          &
                              ndata=n_land_tile)
    call pft_time_axis%initialise("plant_func_time",          &
                                  file_id="plant_func_ancil", &
                                  interp_flag=interp_flag, pop_freq="five_days")
    call setup_ancil_field("canopy_height", depository, ancil_fields,         &
                              mesh, twod_mesh, twod=.true., ndata=npft, &
                              time_axis=pft_time_axis)
    call setup_ancil_field("leaf_area_index", depository, ancil_fields,       &
                              mesh, twod_mesh, twod=.true., ndata=npft, &
                              time_axis=pft_time_axis)
    call pft_time_axis%set_update_behaviour(tmp_update_ptr)
    call ancil_times_list%insert_item(pft_time_axis)

    !=====  SEA ANCILS  =====
    if ( sea_alb_var_chl ) then
      call sea_time_axis%initialise("sea_time", file_id="sea_ancil", &
                                    interp_flag=interp_flag, pop_freq="five_days")
      call setup_ancil_field("chloro_sea", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.,                      &
                              time_axis=sea_time_axis)
      call sea_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(sea_time_axis)
    end if

    call sst_time_axis%initialise("sst_time", file_id="sst_ancil", &
                                  interp_flag=interp_flag, pop_freq="daily")
    call setup_ancil_field("tstar_sea", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.,                   &
                              time_axis=sst_time_axis)
    call sst_time_axis%set_update_behaviour(tmp_update_ptr)
    call ancil_times_list%insert_item(sst_time_axis)

    !=====  SEA ICE ANCILS  =====
    if (.not. l_esm_couple) then
      call sea_ice_time_axis%initialise("sea_ice_time", file_id="sea_ice_ancil", &
                                      interp_flag=interp_flag, pop_freq="daily")
      call setup_ancil_field("sea_ice_thickness", depository, ancil_fields, &
                mesh, twod_mesh, twod=.true., time_axis=sea_ice_time_axis)
      call setup_ancil_field("sea_ice_fraction", depository, ancil_fields, &
                mesh, twod_mesh, twod=.true., time_axis=sea_ice_time_axis)
      call sea_ice_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(sea_ice_time_axis)
    endif

    !=====  RADIATION ANCILS  =====
    if ( albedo_obs ) then
      call albedo_vis_time_axis%initialise("albedo_vis_time",          &
                                           file_id="albedo_vis_ancil", &
                                           interp_flag=interp_flag,    &
                                           pop_freq="five_days")
      call setup_ancil_field("albedo_obs_vis", depository, ancil_fields, &
                             mesh, twod_mesh, twod=.true.,         &
                             time_axis=albedo_vis_time_axis)
      call albedo_vis_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(albedo_vis_time_axis)

      call albedo_nir_time_axis%initialise("albedo_nir_time",          &
                                           file_id="albedo_nir_ancil", &
                                           interp_flag=interp_flag,    &
                                           pop_freq="five_days")
      call setup_ancil_field("albedo_obs_nir", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.,        &
                              time_axis=albedo_nir_time_axis)
      call albedo_nir_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(albedo_nir_time_axis)
    end if

    !=====  SOIL ANCILS  =====
    call setup_ancil_field("soil_albedo", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.)
    if ( l_vary_z0m_soil ) then
      call setup_ancil_field("soil_roughness", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
    endif
    call setup_ancil_field("soil_thermal_cond", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("soil_moist_wilt", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("soil_moist_crit", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("soil_moist_sat", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("soil_cond_sat", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("soil_thermal_cap", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("soil_suction_sat", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("clapp_horn_b", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("mean_topog_index", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("stdev_topog_index", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)

    !=====  OROGRAPHY ANCILS  =====
    call setup_ancil_field("sd_orog", depository, ancil_fields, mesh, &
                              twod_mesh, twod=.true.)
    call setup_ancil_field("grad_xx_orog", depository, ancil_fields, mesh,&
                              twod_mesh, twod=.true.)
    call setup_ancil_field("grad_xy_orog", depository, ancil_fields, mesh,&
                              twod_mesh, twod=.true.)
    call setup_ancil_field("grad_yy_orog", depository, ancil_fields, mesh,&
                              twod_mesh, twod=.true.)
    call setup_ancil_field("peak_to_trough_orog", depository, ancil_fields,  &
                              mesh, twod_mesh, twod=.true.)
    call setup_ancil_field("silhouette_area_orog", depository, ancil_fields, &
                              mesh, twod_mesh, twod=.true.)
    if (topography == topography_slope .or. &
        topography == topography_horizon) then
      call setup_ancil_field("grad_x_orog", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("grad_y_orog", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true.)
    end if
    if (topography == topography_horizon) then
      call setup_ancil_field("horizon_angle", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true., &
                                ndata=n_horiz_ang*n_horiz_layer)
      call setup_ancil_field("horizon_aspect", depository, ancil_fields, &
                                mesh, twod_mesh, twod=.true., &
                                ndata=n_horiz_ang)
    end if

    !=====  OZONE ANCIL  =====
    call ozone_time_axis%initialise("ozone_time", file_id="ozone_ancil", &
                                    interp_flag=interp_flag, pop_freq="monthly")
    call setup_ancil_field("ozone", depository, ancil_fields, mesh, &
                             twod_mesh, time_axis=ozone_time_axis)
    call ozone_time_axis%set_update_behaviour(tmp_update_ptr)
    call ancil_times_list%insert_item(ozone_time_axis)

    !=====  AEROSOL ANCILS  =====
    if (glomap_mode == glomap_mode_climatology) then
      call aerosol_time_axis%initialise("aerosols_time",          &
                                      file_id="aerosols_ancil", &
                                      interp_flag=interp_flag,  &
                                      pop_freq="five_days")
      call setup_ancil_field("acc_sol_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("acc_sol_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("acc_sol_su", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("acc_sol_ss", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("n_acc_sol",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_sol_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_sol_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_sol_su", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("n_ait_sol",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_ins_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("ait_ins_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("n_ait_ins",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_bc", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_om", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_su", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("cor_sol_ss", depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      call setup_ancil_field("n_cor_sol",  depository, ancil_fields, mesh,  &
                             twod_mesh, time_axis=aerosol_time_axis)
      ! The following fields will need adding when dust is available in the
      ! ancillary file:
      !   acc_sol_du, cor_sol_du, n_acc_ins, acc_ins_du, n_cor_ins, cor_ins_du
      call aerosol_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(aerosol_time_axis)
    end if

    !=====  EMISSION ANCILS  =====
    if ( glomap_mode == glomap_mode_ukca   .and.                         &
         ancil_option == ancil_option_updating )  then
      ! -- Single level ancils
      call em_bc_bf_time_axis%initialise("em_bc_bf_time",                &
                                       file_id="emiss_bc_biofuel_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
      call setup_ancil_field("emiss_bc_biofuel", depository, ancil_fields,   &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_bc_bf_time_axis)
      call em_bc_bf_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_bc_bf_time_axis)

      call em_bc_ff_time_axis%initialise("em_bc_ff_time",                &
                                       file_id="emiss_bc_fossil_ancil",  &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
      call setup_ancil_field("emiss_bc_fossil", depository, ancil_fields,    &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_bc_ff_time_axis)
      call em_bc_ff_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_bc_ff_time_axis)

      call em_dms_lnd_time_axis%initialise("em_dms_lnd_time",              &
                                         file_id="emiss_dms_land_ancil",   &
                                         interp_flag=interp_flag,          &
                                         pop_freq="five_days")
      call setup_ancil_field("emiss_dms_land", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_dms_lnd_time_axis)
      call em_dms_lnd_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_dms_lnd_time_axis)

      call dms_ocn_time_axis%initialise("dms_ocn_time",                 &
                                      file_id="dms_conc_ocean_ancil",   &
                                      interp_flag=interp_flag,          &
                                      pop_freq="five_days")
      call setup_ancil_field("dms_conc_ocean", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=dms_ocn_time_axis)
      call dms_ocn_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(dms_ocn_time_axis)

      call em_mterp_time_axis%initialise("em_mterp_time",               &
                                       file_id="emiss_monoterp_ancil",  &
                                       interp_flag=interp_flag,         &
                                       pop_freq="five_days")
      call setup_ancil_field("emiss_monoterp", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_mterp_time_axis)
      call em_mterp_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_mterp_time_axis)

      call em_om_bf_time_axis%initialise("em_om_bf_time",                &
                                       file_id="emiss_om_biofuel_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
      call setup_ancil_field("emiss_om_biofuel", depository, ancil_fields,   &
                           mesh, twod_mesh, twod=.true.,               &
                           time_axis=em_om_bf_time_axis)
      call em_om_bf_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_om_bf_time_axis)

      call em_om_ff_time_axis%initialise("em_om_ff_time",                &
                                       file_id="emiss_om_fossil_ancil",  &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
      call setup_ancil_field("emiss_om_fossil", depository, ancil_fields,   &
                           mesh, twod_mesh, twod=.true.,              &
                           time_axis=em_om_ff_time_axis)
      call em_om_ff_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_om_ff_time_axis)

      call em_so2_lo_time_axis%initialise("em_so2_lo_time",              &
                                        file_id="emiss_so2_low_ancil",   &
                                        interp_flag=interp_flag,         &
                                        pop_freq="five_days")
      call setup_ancil_field("emiss_so2_low", depository, ancil_fields,     &
                           mesh, twod_mesh, twod=.true.,              &
                           time_axis=em_so2_lo_time_axis)
      call em_so2_lo_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_so2_lo_time_axis)

      call em_so2_hi_time_axis%initialise("em_so2_hi_time",              &
                                        file_id="emiss_so2_high_ancil",  &
                                        interp_flag=interp_flag,         &
                                        pop_freq="five_days")
      call setup_ancil_field("emiss_so2_high", depository, ancil_fields,    &
                           mesh, twod_mesh, twod=.true.,              &
                           time_axis=em_so2_hi_time_axis)
      call em_so2_hi_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_so2_hi_time_axis)

      call setup_ancil_field("soil_clay", depository, ancil_fields,         &
                           mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("soil_sand", depository, ancil_fields,         &
                           mesh, twod_mesh, twod=.true.)
      call setup_ancil_field("dust_mrel", depository, ancil_fields,         &
                           mesh, twod_mesh, twod=.true., ndata=ndiv)

      ! -- 3-D ancils
      !-- natural SO2 emissions, currently single-time
      call setup_ancil_field("emiss_so2_nat", depository, ancil_fields,     &
                             mesh, twod_mesh)

      call em_bc_bb_time_axis%initialise("em_bc_bb_time",                &
                                       file_id="emiss_bc_biomass_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
      call setup_ancil_field("emiss_bc_biomass", depository, ancil_fields,  &
                           mesh, twod_mesh,                           &
                           time_axis=em_bc_bb_time_axis)   ! 3-D
      call em_bc_bb_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_bc_bb_time_axis)

      call em_om_bb_time_axis%initialise("em_om_bb_time",                &
                                       file_id="emiss_om_biomass_ancil", &
                                       interp_flag=interp_flag,          &
                                       pop_freq="five_days")
      call setup_ancil_field("emiss_om_biomass", depository, ancil_fields,  &
                           mesh, twod_mesh,                           &
                           time_axis=em_om_bb_time_axis)
      call em_om_bb_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(em_om_bb_time_axis)

      !=====  OFFLINE OXIDANT ANCILS  =====
      call h2o2_limit_time_axis%initialise("h2o2_limit_time",         &
                                         file_id="h2o2_limit_ancil",  &
                                         interp_flag=interp_flag,     &
                                         pop_freq="five_days")
      call setup_ancil_field("h2o2_limit", depository, ancil_fields,        &
                           mesh, twod_mesh,                           &
                           time_axis=h2o2_limit_time_axis)
      call h2o2_limit_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(h2o2_limit_time_axis)

      call ho2_time_axis%initialise("ho2_time", file_id="ho2_ancil", &
                                    interp_flag=interp_flag,         &
                                    pop_freq="five_days")
      call setup_ancil_field("ho2", depository, ancil_fields,               &
                           mesh, twod_mesh, time_axis=ho2_time_axis)
      call ho2_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(ho2_time_axis)

      call no3_time_axis%initialise("no3_time", file_id="no3_ancil", &
                                    interp_flag=interp_flag,         &
                                    pop_freq="five_days")
      call setup_ancil_field("no3", depository, ancil_fields,               &
                           mesh, twod_mesh, time_axis=no3_time_axis)
      call no3_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(no3_time_axis)

      call o3_time_axis%initialise("o3_time", file_id="o3_ancil",    &
                                   interp_flag=interp_flag,          &
                                   pop_freq="five_days")
      call setup_ancil_field("o3", depository, ancil_fields,                &
                           mesh, twod_mesh, time_axis=o3_time_axis)
      call o3_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(o3_time_axis)

      call oh_time_axis%initialise("oh_time", file_id="oh_ancil",    &
                                   interp_flag=interp_flag,          &
                                   pop_freq="five_days")
      call setup_ancil_field("oh", depository, ancil_fields,                &
                           mesh, twod_mesh, time_axis=oh_time_axis)
      call oh_time_axis%set_update_behaviour(tmp_update_ptr)
      call ancil_times_list%insert_item(oh_time_axis)

    endif  ! ancil_updating, glomap_ukca

    ! Now the field collection is set up, the fields will be initialised in
    ! gungho_model_data_mod

  end subroutine create_fd_ancils

  !> @details Adds fields to the ancil collection, sets up their read and write
  !>      behaviour and creates them in the depository if they do not yet exist
  !> @param[in] name The field name
  !> @param[in, out] depository The depository field collection
  !> @param[in, out] ancil_fields The ancil field collection
  !> @param[in] mesh                The current 3d mesh
  !> @param[in, optional] twod_mesh The current 2d mesh
  !> @param[in, optional] ndata Number of non-spatial dimensions for multi-data
  !>                            field
  !> @param[in, out, optional] time_axis Time axis associated with ancil field
  subroutine setup_ancil_field( name, depository, ancil_fields, mesh, &
                                twod_mesh, twod, ndata, time_axis )

    implicit none

    character(*),                   intent(in)    :: name
    type( field_collection_type ),  intent(inout) :: depository
    type( field_collection_type ),  intent(inout) :: ancil_fields
    type( mesh_type ),    pointer,  intent(in)    :: mesh
    type( mesh_type ),    pointer,  intent(in)    :: twod_mesh
    logical(l_def),       optional, intent(in)    :: twod
    integer(i_def),       optional, intent(in)    :: ndata
    type(time_axis_type), optional, intent(inout) :: time_axis

    ! Local variables
    type(field_type)          :: new_field
    integer(i_def)            :: ndat, time_ndat
    logical(l_def)            :: twod_field
    integer(i_def), parameter :: fs_order = 0

    ! Pointers
    type(function_space_type),       pointer :: vec_space => null()
    procedure(read_interface),       pointer :: tmp_read_ptr => null()
    procedure(write_interface),      pointer :: tmp_write_ptr => null()
    class(field_type),               pointer :: fld_ptr => null()
    class(pure_abstract_field_type), pointer :: abs_fld_ptr => null()

    ! Set field ndata if argument is present, else leave as default value
    if (present(ndata)) then
      ndat = ndata
    else
      ndat = 1
    end if
    if (present(twod)) then
      twod_field = twod
    else
      twod_field = .false.
    end if

    ! If field does not yet exist, then create it
    if ( .not. depository%field_exists( name ) ) then
      write(log_scratch_space,'(3A,I6)') &
           "Creating new field for ", trim(name)
      call log_event(log_scratch_space,LOG_LEVEL_INFO)
      if (twod_field) then
        vec_space => function_space_collection%get_fs( twod_mesh, fs_order, &
                                                       W3, ndat )
        tmp_write_ptr => write_field_single_face
      else
        vec_space => function_space_collection%get_fs( mesh, fs_order, &
                                                       WTheta, ndat )
        tmp_write_ptr => write_field_face
       end if
      call new_field%initialise( vec_space, name=trim(name) )
      call new_field%set_write_behaviour(tmp_write_ptr)
      ! Add the new field to the field depository
      call depository%add_field(new_field)
    end if

    ! If field is time-varying, also create field storing raw data to be
    ! interpolated
    if (present(time_axis)) then
      write(log_scratch_space,'(3A,I6)') &
           "Creating time axis field for ", trim(name)
      call log_event(log_scratch_space,LOG_LEVEL_INFO)

      ! Multiply ndat by the number of time windows
      time_ndat = ndat * time_axis%get_window_size()
      if (twod_field) then
        vec_space => function_space_collection%get_fs( twod_mesh, fs_order, &
                                                       W3, time_ndat )
      else
        vec_space => function_space_collection%get_fs( mesh, fs_order, &
                                                       WTheta, time_ndat )
      end if
      call new_field%initialise( vec_space, name=trim(name) )
      call time_axis%add_field(new_field)
    end if

    ! Get a field pointer from the depository
    fld_ptr => depository%get_field(name)

    if (.not. present(time_axis)) then
      !Set up field read behaviour for 2D and 3D fields
      if (twod_field) then
        tmp_read_ptr => read_field_single_face
      else
        tmp_read_ptr => read_field_face
      end if
      ! Set field read behaviour for target field
      call fld_ptr%set_read_behaviour(tmp_read_ptr)
    end if

    ! Add the field pointer to the target field collection
    abs_fld_ptr => depository%get_field(name)
    call ancil_fields%add_reference_to_field(abs_fld_ptr)

    ! Nullify pointers
    nullify(vec_space)
    nullify(tmp_read_ptr)
    nullify(tmp_write_ptr)
    nullify(fld_ptr)

  end subroutine setup_ancil_field

end module init_ancils_mod
