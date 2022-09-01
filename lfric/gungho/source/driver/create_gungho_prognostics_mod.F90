!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Create empty fields for use by the gungho model
!> @details Creates the empty prognostic fields that will be later
!>          initialise and used by the gungho model. Field creation
!>          consists of three stages: constructing the field object,
!>          placing it in the depository (so it doesn't go out of scope) and
!>          putting a pointer to the depository version of the field into
!>          a 'prognostic_fields' field collection

module create_gungho_prognostics_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use field_parent_mod,               only : write_interface, &
                                             checkpoint_write_interface, &
                                             checkpoint_read_interface
  use field_collection_mod,           only : field_collection_type
  use finite_element_config_mod,      only : element_order
  use fs_continuity_mod,              only : W0, W2, W3, Wtheta
  use function_space_collection_mod , only : function_space_collection
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO
  use mesh_mod,                       only : mesh_type
  use mr_indices_mod,                 only : nummr, &
                                             mr_names
  use moist_dyn_mod,                  only : num_moist_factors
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use lfric_xios_read_mod,            only : checkpoint_read_xios
  use lfric_xios_write_mod,           only : write_field_node, &
                                             write_field_face, &
                                             checkpoint_write_xios
  use io_mod,                         only : checkpoint_write_netcdf, &
                                             checkpoint_read_netcdf
  use io_config_mod,                  only : use_xios_io,     &
                                             write_diag,      &
                                             checkpoint_read, &
                                             checkpoint_write
  use derived_config_mod,             only : l_esm_couple
  implicit none

  private
  public :: create_gungho_prognostics

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Create empty fields to be used as prognostics by the gungho model
  !> @param[in]    mesh       The current 3d mesh
  !> @param[inout] depository A collection of all fields that need to be
  !>                          kept in scope
  !> @param[inout] prognostic_fields A collection of the fields that make up the
  !>                                 prognostic variables in the model
  !> @param[inout] diagnostic_fields A collection of the fields that make up the
  !>                                 diagnostic variables in the model
  !> @param[inout] adv_fields_all_outer A collection of fields to be advected every outer iteration
  !> @param[inout] adv_fields_last_outer A collection of all fields to be advected on last outer iteration
  !> @param[inout] mr An array of fields that hold the moisture mixing ratios
  !> @param[inout] moist_dyn An array of the moist dynamics fields
  subroutine create_gungho_prognostics( mesh, depository, &
                                        prognostic_fields, diagnostic_fields, &
                                        adv_fields_all_outer, &
                                        adv_fields_last_outer, &
                                        mr, moist_dyn )
    implicit none

    type(mesh_type), intent(in), pointer :: mesh

    type(field_collection_type), intent(inout):: depository
    type(field_collection_type), intent(inout):: prognostic_fields
    type(field_collection_type), intent(out)  :: diagnostic_fields
    type(field_collection_type), intent(inout):: adv_fields_all_outer
    type(field_collection_type), intent(inout):: adv_fields_last_outer

    type( field_type ), intent(inout), target :: mr(nummr)
    type( field_type ), intent(inout)         :: moist_dyn(num_moist_factors)

    class(pure_abstract_field_type), pointer  :: tmp_ptr => null()

    integer(i_def)                            :: imr

    procedure(write_interface),            pointer :: tmp_write_ptr => null()
    procedure(checkpoint_write_interface), pointer :: tmp_checkpoint_write_ptr => null()
    procedure(checkpoint_read_interface),  pointer :: tmp_checkpoint_read_ptr => null()

    ! Temp fields to create prognostics
    type( field_type )                         :: u, rho, theta, exner

    logical                                    :: create_depository

    call log_event( 'GungHo: Creating prognostics...', LOG_LEVEL_INFO )

    ! Create the depository, prognostics and diagnostics field collections.
    ! if coupled configuration depository and prognostics created in cpl_fields
!> @todo this is a temporary solution related to the limitation of the XIOS
!>       currently used. This routine should return to original version
!>       when cpl_define can be called in any place in the code.
!>       See #2710 test branch for details.

    if(l_esm_couple) then
      create_depository = .false.
    else
      create_depository = .true.
    endif

    if (create_depository) then
       call depository%initialise(name='depository', table_len=100)
       call prognostic_fields%initialise(name="prognostics", table_len=100)
    endif
    call diagnostic_fields%initialise(name="diagnostics", table_len=100)
    ! Create collection of fields to be advected
    call adv_fields_last_outer%initialise(name='adv_fields_last_outer', table_len=100)
    call adv_fields_all_outer%initialise(name='adv_fields_all_outer', table_len=100)

    ! Create prognostic fields
    call theta%initialise( vector_space = &
                        function_space_collection%get_fs(mesh, element_order, Wtheta), &
                        name= "theta" )
    call u%initialise( vector_space = &
                        function_space_collection%get_fs(mesh, element_order, W2), &
                        name = "u" )
    call rho%initialise( vector_space = &
                        function_space_collection%get_fs(mesh, element_order, W3), &
                        name = "rho" )
    call exner%initialise( vector_space = &
                        function_space_collection%get_fs(mesh, element_order, W3), &
                        name = "exner" )

    ! The moisture mixing ratio fields (mr) and moist dynamics fields
    ! (moist_dyn) are always passed into the timestep algorithm, so are
    ! always created here, even when moisture_formulation = 'dry'
    do imr = 1,nummr
      call mr(imr)%initialise( vector_space = &
      function_space_collection%get_fs(mesh, element_order, theta%which_function_space()), &
                            name = trim(mr_names(imr)) )
    end do

    ! Auxilliary fields holding moisture-dependent factors for dynamics
    do imr = 1, num_moist_factors
      call moist_dyn(imr)%initialise( vector_space = &
      function_space_collection%get_fs(mesh, element_order, theta%which_function_space()) )
    end do

    ! Set I/O behaviours for diagnostic output

    if (write_diag .and. use_xios_io) then

       ! Set diagnostic output handlers

       ! Face domain

       tmp_write_ptr => write_field_face

       ! Vector fields that are projected to scalar components
       call u%set_write_behaviour(tmp_write_ptr)

       ! Scalar fields
       call rho%set_write_behaviour(tmp_write_ptr)
       call exner%set_write_behaviour(tmp_write_ptr)

       ! Theta is a special case as it can be on face (if function space is WTheta)
       ! or node (if function space is W0)
       if (theta%which_function_space() == Wtheta) then

         call theta%set_write_behaviour(tmp_write_ptr)

       else

        tmp_write_ptr => write_field_node

        call theta%set_write_behaviour(tmp_write_ptr)

       end if

       ! Moisture uses the same type of field write as Theta

       call theta%get_write_behaviour(tmp_write_ptr)

       do imr = 1,nummr
         call mr(imr)%set_write_behaviour(tmp_write_ptr)
       end do

    end if

    if ( checkpoint_write .or. checkpoint_read) then

      if ( use_xios_io ) then

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

      call u%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call rho%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call theta%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
      call exner%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)

      call u%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call rho%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call theta%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      call exner%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)

      do imr = 1,nummr
        call mr(imr)%set_checkpoint_write_behaviour(tmp_checkpoint_write_ptr)
        call mr(imr)%set_checkpoint_read_behaviour(tmp_checkpoint_read_ptr)
      end do

    end if

    ! Populate the depository
    call depository%add_field( theta )
    call depository%add_field( rho )
    call depository%add_field( u )
    call depository%add_field( exner )

    ! Populate the prognostic field collection
    tmp_ptr => depository%get_field('theta')
    call prognostic_fields%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('rho')
    call prognostic_fields%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('u')
    call prognostic_fields%add_reference_to_field(tmp_ptr)
    tmp_ptr => depository%get_field('exner')
    call prognostic_fields%add_reference_to_field(tmp_ptr)
    ! The moisture mixing ratios always need checkpointing otherwise
    ! they are uninitialised on a restart
    do imr = 1,nummr
      tmp_ptr => mr(imr)
      call prognostic_fields%add_reference_to_field(tmp_ptr)
    end do

    nullify( tmp_write_ptr, tmp_checkpoint_write_ptr, tmp_checkpoint_read_ptr )

  end subroutine create_gungho_prognostics


end module create_gungho_prognostics_mod
