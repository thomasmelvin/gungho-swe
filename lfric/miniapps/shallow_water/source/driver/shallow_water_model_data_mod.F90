!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for shallow water model run working data set inclduing methods
!!        to initialise, copy and finalise the data set.
!!
!> @details This module provides a type to hold all the model fields and methods to
!!          initialise (create and read), copy and finalise (write and destroy) the
!!          data contained within the type.
!!
module shallow_water_model_data_mod

  use clock_mod,                            only: clock_type
  use field_mod,                            only: field_type
  use field_parent_mod,                     only: write_interface
  use field_collection_mod,                 only: field_collection_type
  use files_config_mod,                     only: checkpoint_stem_name
  use constants_mod,                        only: i_def, l_def
  use log_mod,                              only: log_event,      &
                                                  LOG_LEVEL_INFO, &
                                                  LOG_LEVEL_ERROR
  use io_config_mod,                        only: checkpoint_read,  &
                                                  checkpoint_write, &
                                                  write_dump
  use lfric_xios_read_mod,                  only: read_checkpoint,  &
                                                  read_state
  use lfric_xios_write_mod,                 only: write_checkpoint, &
                                                  write_state
  use mesh_mod,                             only: mesh_type
  use create_shallow_water_prognostics_mod, only: create_shallow_water_prognostics
  use swe_init_fields_alg_mod,              only: swe_init_fields_alg
  use linked_list_mod,                      only: linked_list_type
  use variable_fields_mod,                  only: init_variable_fields

  implicit none

  private

  !> Holds the working data set for a model run and other working state.
  !>
  type :: model_data_type

    private

    !> Stores all the fields used by the model - the remaining collections store
    !> pointers to the data in the depository.
    type(field_collection_type), public :: depository
    !> All the prognostic fields
    type(field_collection_type), public   :: prognostic_fields
    !> All the diagnostic fields
    type(field_collection_type), public   :: diagnostic_fields
    !> Surface geopotential field
    type(field_type), public  :: s_geopot

  end type model_data_type

  public :: model_data_type,       &
            create_model_data,     &
            finalise_model_data,   &
            initialise_model_data, &
            output_model_data

contains

  !=============================================================================
  !> @brief Create the fields contained in model_data.
  !> @param[in,out] model_data   The working data set for a model run
  !> @param[in]     mesh        Mesh to initialise variables on
  subroutine create_model_data( model_data,   &
                                mesh )

    implicit none

    type(model_data_type), intent(inout) :: model_data
    type(mesh_type), pointer, intent(in) :: mesh

    !-------------------------------------------------------------------------
    ! Instantiate the fields
    !-------------------------------------------------------------------------

    ! Create prognostics

    call create_shallow_water_prognostics( mesh,                         &
                                           model_data%depository,        &
                                           model_data%prognostic_fields, &
                                           model_data%s_geopot )

  end subroutine create_model_data

  !=============================================================================
  !> @brief Initialises the working data set dependent of namelist configuration.
  !> @param[in,out] model_data The working data set for a model run
  !> @param[in]     mesh       Mesh to initialise variables on
  !> @param[in]     clock      Model time
  subroutine initialise_model_data( model_data, &
                                    mesh,       &
                                    clock )

    implicit none

    type(model_data_type), intent(inout) :: model_data
    type(mesh_type), pointer, intent(in) :: mesh
    class(clock_type),        intent(in) :: clock

    ! Initialise prognostic fields
    if (checkpoint_read) then                 ! Recorded check point to start from
      call read_checkpoint(model_data%depository, &
                           clock%get_first_step() - 1, checkpoint_stem_name)

    else                                      ! No check point to start from
      call swe_init_fields_alg(mesh,                         &
                               model_data%prognostic_fields, &
                               model_data%s_geopot)
    end if

  end subroutine initialise_model_data

  !=============================================================================
  !> @brief Writes out a checkpoint and dump file dependent on namelist
  !!        options.
  !> @param[in,out] model_data The working data set for the model run
  !> @param[in]     clock      Model time.
  subroutine output_model_data( model_data, &
                                clock )

    implicit none

    type(model_data_type), intent(inout), target :: model_data
    class(clock_type),     intent(in)            :: clock

    type(field_collection_type), pointer :: prognostic_fields => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields

    !=================== Write fields to checkpoint files ====================!
    if( checkpoint_write ) then
       call write_checkpoint( prognostic_fields, clock, checkpoint_stem_name )
    end if

  end subroutine output_model_data

  !=============================================================================
  !> @brief Routine to destroy all the field collections in the working data set.
  !> @param[in,out] model_data The working data set for a model run
  subroutine finalise_model_data( model_data )

    implicit none

    type(model_data_type), intent(inout) :: model_data

    ! Clear all the fields in each field collection
    call model_data%depository%clear()
    call model_data%prognostic_fields%clear()
    call model_data%diagnostic_fields%clear()

    call log_event( 'finalise_model_data: all fields have been cleared', &
                     LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module shallow_water_model_data_mod
