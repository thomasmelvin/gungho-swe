!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation and finalisation of IO within the
!! driver layer. This module also contains the model's io_context object and
!! associated getter routines
!>
module driver_io_mod

  use clock_mod,               only: clock_type
  use constants_mod,           only: i_native
  use field_mod,               only: field_type
  use file_mod,                only: file_type
  use io_context_mod,          only: io_context_type
  use io_config_mod,           only: use_xios_io
  use log_mod,                 only: log_event, log_level_error
  use time_config_mod,         only: timestep_end, timestep_start,  &
                                     calendar_start, calendar_type, &
                                     key_from_calendar_type
  use timestepping_config_mod, only: dt, spinup_period
  use simple_io_context_mod,   only: simple_io_context_type
#ifdef USE_XIOS
  use lfric_xios_context_mod,  only: lfric_xios_context_type
  use lfric_xios_file_mod,     only: lfric_xios_file_type
#endif

  implicit none

  public :: init_io, final_io,  &
            get_io_context,     &
            get_clock,          &
            filelist_populator, &
            append_file_to_list
  private

  class(io_context_type), allocatable, target :: context

  abstract interface
    subroutine filelist_populator(files_list)
      import file_type
      class(file_type), allocatable, intent(out) :: files_list(:)
    end subroutine filelist_populator
  end interface

contains

  !> @brief  Initialises the model I/O
  !>
  !> @param[in] id                A string identifier for the model
  !> @param[in] communicator      The ID for the model MPI communicator
  !> @param[in] chi               The model coordinate field
  !> @param[in] panel_id          Field containing the panel ID for each mesh
  !!                              vertex
  !> @param[in] populate_filelist Optional procedure for creating a list of
  !!                              file descriptions used by the model I/O
  subroutine init_io( id, communicator, &
                      chi, panel_id,    &
                      populate_filelist )

    implicit none

    character(*),                    intent(in) :: id
    integer(i_native),               intent(in) :: communicator
    class(field_type),               intent(in) :: chi(:)
    class(field_type),               intent(in) :: panel_id
    procedure(filelist_populator), &
                  optional, pointer, intent(in) :: populate_filelist

    class(file_type), allocatable :: file_list(:)
    integer(i_native) :: rc

    ! Allocate IO context type based on model configuration
    if ( use_xios_io ) then
#ifdef USE_XIOS
      allocate( lfric_xios_context_type::context, stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate LFRic-XIOS context object", &
                        log_level_error )
      end if
#endif

    else
      allocate( simple_io_context_type::context, stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate simple I/O context object", &
                        log_level_error )
      end if

    end if

    ! Populate list of files if present and intialise context object
    if (present(populate_filelist)) then
      call populate_filelist(file_list)
      call context%initialise( id, communicator,                      &
                               chi, panel_id,                         &
                               timestep_start,                        &
                               timestep_end,                          &
                               spinup_period, dt,                     &
                               calendar_start,                        &
                               key_from_calendar_type(calendar_type), &
                               file_list )

    else
      call context%initialise( id, communicator,                      &
                               chi, panel_id,                         &
                               timestep_start,                        &
                               timestep_end,                          &
                               spinup_period, dt,                     &
                               calendar_start,                        &
                               key_from_calendar_type(calendar_type) )

    end if

  end subroutine init_io

  !> @brief  Finalises the model I/O
  subroutine final_io()

    implicit none

    deallocate(context)

  end subroutine final_io

  !> @brief  Returns the model io context.
  function get_io_context() result(context_ptr)

    implicit none

    class(io_context_type), pointer :: context_ptr

    context_ptr => context

  end function get_io_context

  !> @brief  Returns the model clock.
  function get_clock() result(clock_ptr)

    implicit none

    class(clock_type), pointer :: clock_ptr

    clock_ptr => context%get_clock()

  end function get_clock

  !> @brief  Appends a file to the end of a list of files.
  !>
  !> @param[in]      file      File to be added
  !> @param[in,out]  filelist  List of files to be added to
  subroutine append_file_to_list(file, filelist)

    implicit none

    class(file_type),              intent(in)    :: file
    class(file_type), allocatable, intent(inout) :: filelist(:)

    class(file_type), allocatable :: new_filelist(:)

    if (.not. allocated(filelist)) then
      select type(file)
#ifdef USE_XIOS
      type is (lfric_xios_file_type)
        allocate(lfric_xios_file_type::new_filelist(1))

        ! We just allocated this but it is still polymorphic, so we need to cast
        select type(new_filelist)
        type is (lfric_xios_file_type)
          new_filelist(1) = file
        end select ! type(new_filelist)
#endif
      end select ! type(file)

    else
      select type(filelist)
#ifdef USE_XIOS
      type is (lfric_xios_file_type)
        allocate(lfric_xios_file_type::new_filelist(size(filelist)+1))

        ! We just allocated this but it is still polymorphic, so we need to cast
        select type(new_filelist)
        type is (lfric_xios_file_type)
          new_filelist(1:size(filelist)) = filelist

          select type(file)
          type is (lfric_xios_file_type)
            new_filelist(size(filelist)+1) = file
          end select ! type(file)

        end select ! type(new_filelist)
#endif

      end select ! type(filelist)

    end if

    call move_alloc(new_filelist, filelist)

  end subroutine append_file_to_list

end module driver_io_mod