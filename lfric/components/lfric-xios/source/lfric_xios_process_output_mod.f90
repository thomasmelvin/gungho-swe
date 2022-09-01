!-------------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief  Module to hold routines for processing XIOS output for compliance
!>         with downstream data requirements
!>
module lfric_xios_process_output_mod

  use io_config_mod,       only: file_convention,       &
                                 file_convention_ugrid, &
                                 file_convention_cf
  use lfric_ncdf_file_mod, only: lfric_ncdf_file_type, &
                                 LFRIC_NCDF_WRITE,     &
                                 LFRIC_NCDF_OPEN
  use mpi_mod,             only: get_comm_rank

  implicit none

  public :: process_output_file
  private

contains

!> @brief Processes a NetCDF file produced by XIOS to align with the LFRic
!!        UGRID file format
!>
!> @param[in] file_path  The path to the NetCDF file to be edited
subroutine process_output_file(file_path)

  implicit none

  character(len=*), intent(in) :: file_path

  type(lfric_ncdf_file_type) :: file_ncdf
  logical                    :: file_exists

  ! Output processing must be done in serial
  if (get_comm_rank() /= 0) return

  ! If file has not been written out, then don't attempt to process it
  inquire(file=trim(file_path), exist=file_exists)
  if (.not. file_exists) return

  ! Open output file
  file_ncdf = lfric_ncdf_file_type( trim(file_path),           &
                                    open_mode=LFRIC_NCDF_OPEN, &
                                    io_mode=LFRIC_NCDF_WRITE )

  call format_version(file_ncdf)

  call file_ncdf%close_file()

end subroutine process_output_file

!> @brief Tags output file with the current version number of the LFRic file
!!        format
!>
!> @param[in] file_ncdf  The netcdf file to be edited
subroutine format_version(file_ncdf)

  implicit none

  type(lfric_ncdf_file_type), intent(inout) :: file_ncdf

  call file_ncdf%set_attribute("description", "LFRic file format v0.1.0")

  select case(file_convention)
  case (file_convention_ugrid)
    call file_ncdf%set_attribute("Conventions", "UGRID-1.0")

  case (file_convention_cf)
    call file_ncdf%set_attribute("Conventions", "CF")

  end select

end subroutine format_version

end module lfric_xios_process_output_mod