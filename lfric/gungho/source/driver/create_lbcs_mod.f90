!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create LBC fields.
!> @details Create LBC field collection and add fields.
module create_lbcs_mod

  use constants_mod,              only : i_def, l_def, str_def
  use log_mod,                    only : log_event,             &
                                         log_scratch_space,     &
                                         LOG_LEVEL_INFO,        &
                                         LOG_LEVEL_ERROR
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use fs_continuity_mod,          only : W0, W2, W3, Wtheta, W2h
  use function_space_mod,         only : function_space_type
  use mesh_mod,                   only : mesh_type
  use mr_indices_mod,             only : nummr,                      &
                                         mr_names
  use linked_list_mod,            only : linked_list_type
  use lfric_xios_time_axis_mod,   only : time_axis_type,        &
                                         update_interface
  use lfric_xios_read_mod,        only : read_field_time_var
  use init_time_axis_mod,         only : setup_field
  use initialization_config_mod,  only : lbc_option,             &
                                         lbc_option_analytic,    &
                                         lbc_option_gungho_file, &
                                         lbc_option_um2lfric_file

  implicit none

  public  :: create_lbc_fields

  contains

  !> @brief   Create and add LBC fields.
  !> @details Create the lateral boundary condition field collection.
  !!          On every timestep these fields will be updated and used by the
  !!          limited area model.
  !> @param[in]     mesh              The current 3d mesh.
  !> @param[in,out] depository        Main collection of all fields in memory.
  !> @param[in,out] prognostic_fields The prognostic variables in the model.
  !> @param[in,out] lbc_fields        The lbc_fields used on every timestep.
  subroutine create_lbc_fields( mesh, depository, prognostic_fields, &
                                lbc_fields, lbc_times_list )

    implicit none

    type(mesh_type), intent(in), pointer :: mesh

    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(inout) :: lbc_fields
    type(linked_list_type),      intent(out)   :: lbc_times_list

    logical(l_def)                             :: checkpoint_restart_flag
    procedure(update_interface), pointer       :: tmp_update_ptr => null()

    type(time_axis_type), save                 :: lbc_time_axis
    logical(l_def),   parameter                :: cyclic=.false.
    logical(l_def),   parameter                :: interp_flag=.true.
    character(len=*), parameter                :: axis_id="lbc_axis"
    character(str_def)                         :: name
    integer(i_def)                             :: imr

    call log_event( 'Create LBC fields', LOG_LEVEL_INFO )

    call lbc_fields%initialise( name='lbc_fields', table_len=100 )

    select case( lbc_option )

      case ( lbc_option_analytic )

        checkpoint_restart_flag = .true.

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_theta", Wtheta, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_u", W2, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_rho", W3, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_exner", W3, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "boundary_u_diff", W2, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "boundary_u_driving", W2, mesh, checkpoint_restart_flag )

        do imr = 1, nummr
          name = trim('lbc_' // adjustl(mr_names(imr)) )
          call setup_field( lbc_fields, depository, prognostic_fields, &
             name, wtheta, mesh, checkpoint_restart_flag )
        enddo

      case ( lbc_option_gungho_file )

        checkpoint_restart_flag = .false.
        ! Set pointer to time axis read behaviour
        tmp_update_ptr => read_field_time_var

        call lbc_time_axis%initialise( "lbc_time", file_id="lbc", yearly=cyclic, &
                                       interp_flag = interp_flag )

        !------ Fields updated directly from LBC file-----------------

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_theta", Wtheta, mesh, checkpoint_restart_flag,    &
            time_axis=lbc_time_axis )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_rho", W3, mesh, checkpoint_restart_flag,          &
            time_axis=lbc_time_axis )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_exner", W3, mesh, checkpoint_restart_flag,        &
            time_axis=lbc_time_axis )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_h_u", W2h, mesh, checkpoint_restart_flag,         &
            time_axis=lbc_time_axis )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_v_u", Wtheta, mesh, checkpoint_restart_flag,      &
            time_axis=lbc_time_axis )

        call lbc_time_axis%set_update_behaviour(tmp_update_ptr)
        call lbc_times_list%insert_item(lbc_time_axis)

        !----- Fields derived from the fields in the LBC file---------

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_u", W2, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "boundary_u_diff", W2, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "boundary_u_driving", W2, mesh, checkpoint_restart_flag )

      case ( lbc_option_um2lfric_file )

        checkpoint_restart_flag = .false.

        ! Set pointer to time axis read behaviour
        tmp_update_ptr => read_field_time_var

        call lbc_time_axis%initialise( "lbc_time", file_id="lbc", yearly=cyclic, &
                                       interp_flag = interp_flag )

        !------ Fields updated directly from LBC file-----------------

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_theta", Wtheta, mesh, checkpoint_restart_flag,    &
            time_axis=lbc_time_axis )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_rho_r2", W3, mesh, checkpoint_restart_flag,       &
            time_axis=lbc_time_axis )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_h_u", W2h, mesh, checkpoint_restart_flag,         &
            time_axis=lbc_time_axis  )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_v_u", Wtheta, mesh, checkpoint_restart_flag,      &
            time_axis=lbc_time_axis  )

        ! Specific humidities (Initially these are read in, but
        ! in the long run we should move to just use the mixing
        ! ratios)
        call setup_field( lbc_fields, depository, prognostic_fields, &
           'lbc_q', wtheta, mesh, checkpoint_restart_flag,        &
           time_axis=lbc_time_axis )
        call setup_field( lbc_fields, depository, prognostic_fields, &
           'lbc_qcl', wtheta, mesh, checkpoint_restart_flag,      &
           time_axis=lbc_time_axis )
        call setup_field( lbc_fields, depository, prognostic_fields, &
           'lbc_qcf', wtheta, mesh, checkpoint_restart_flag,      &
           time_axis=lbc_time_axis )
        call setup_field( lbc_fields, depository, prognostic_fields, &
           'lbc_qrain', wtheta, mesh, checkpoint_restart_flag,    &
           time_axis=lbc_time_axis )

        call lbc_time_axis%set_update_behaviour(tmp_update_ptr)
        call lbc_times_list%insert_item(lbc_time_axis)

        !----- Fields derived from the fields in the LBC file---------
        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_rho", W3, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_exner", W3, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "lbc_u", W2, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "boundary_u_diff", W2, mesh, checkpoint_restart_flag )

        call setup_field( lbc_fields, depository, prognostic_fields, &
           "boundary_u_driving", W2, mesh, checkpoint_restart_flag )

        ! Mixing ratios
        name = trim('lbc_m_v')
        call setup_field( lbc_fields, depository, prognostic_fields, &
             name, wtheta, mesh, checkpoint_restart_flag )
        name = trim('lbc_m_cl')
        call setup_field( lbc_fields, depository, prognostic_fields, &
             name, wtheta, mesh, checkpoint_restart_flag )
        name = trim('lbc_m_ci')
        call setup_field( lbc_fields, depository, prognostic_fields, &
             name, wtheta, mesh, checkpoint_restart_flag )
        name = trim('lbc_m_r')
        call setup_field( lbc_fields, depository, prognostic_fields, &
             name, wtheta, mesh, checkpoint_restart_flag )

      case default
        call log_event( 'This lbc_option not available', LOG_LEVEL_ERROR )
    end select

  end subroutine create_lbc_fields

end module create_lbcs_mod
