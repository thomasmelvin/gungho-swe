!-------------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Sets up I/O configuration from within shallow_water.
!> @details Collects configuration information relevant for the I/O subsystem
!!          and formats it so that it can be passed to the infrastructure.
module shallow_water_setup_io_mod

  use constants_mod,                 only: i_def, i_native, &
                                           str_def, str_max_filename
  use driver_io_mod,                 only: append_file_to_list
  use file_mod,                      only: file_type
  use files_config_mod,              only: checkpoint_stem_name
  use io_config_mod,                 only: diagnostic_frequency, &
                                           checkpoint_write,     &
                                           checkpoint_read,      &
                                           write_dump
  use lfric_xios_file_mod,           only: lfric_xios_file_type
  use linked_list_mod,               only: linked_list_type, &
                                           linked_list_item_type
  use time_config_mod,               only: timestep_start, &
                                           timestep_end

  implicit none

  private
  public :: init_shallow_water_files

  contains

  !> @brief  Sets up I/O configuration.
  !> @param[out] files_list The list of I/O files
  subroutine init_shallow_water_files(files_list)

    implicit none

    class(file_type), allocatable, intent(out) :: files_list(:)

    type(lfric_xios_file_type)      :: tmp_file
    character(len=str_max_filename) :: checkpoint_write_fname, &
                                       checkpoint_read_fname
    integer(i_def)                  :: ts_start, ts_end
    integer(i_native)               :: rc

    ! Get time configuration in integer form
    read(timestep_start,*,iostat=rc)  ts_start
    read(timestep_end,*,iostat=rc)    ts_end

    ! Setup diagnostic output file
    call tmp_file%file_new("lfric_diag")
    call tmp_file%configure(xios_id="lfric_diag", freq=diagnostic_frequency)
    call append_file_to_list(tmp_file, files_list)

    ! Setup checkpoint writing context information
    if ( checkpoint_write ) then
      ! Create checkpoint filename from stem and end timestep
      write(checkpoint_write_fname,'(A,A,I6.6)') &
                           trim(checkpoint_stem_name),"_", ts_end
      call tmp_file%file_new(checkpoint_write_fname)
      call tmp_file%configure(xios_id="lfric_checkpoint_write", &
                              freq=ts_end,                      &
                              field_group_id="checkpoint_fields")
      call append_file_to_list(tmp_file, files_list)
    end if

    ! Setup checkpoint reading context information
    if ( checkpoint_read ) then
      ! Create checkpoint filename from stem and (start - 1) timestep
      write(checkpoint_read_fname,'(A,A,I6.6)') &
                   trim(checkpoint_stem_name),"_", (ts_start - 1)
      call tmp_file%file_new(checkpoint_read_fname)
      call tmp_file%configure(xios_id="lfric_checkpoint_read", &
                              freq=ts_start - 1,               &
                              field_group_id="checkpoint_fields")
      call append_file_to_list(tmp_file, files_list)
    end if

  end subroutine init_shallow_water_files

end module shallow_water_setup_io_mod
