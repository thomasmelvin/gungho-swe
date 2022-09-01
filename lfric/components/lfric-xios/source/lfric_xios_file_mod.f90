!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Module containing a file type for XIOS interface
!>
module lfric_xios_file_mod

  use constants_mod,                 only: i_def, str_def, str_max_filename
  use lfric_xios_process_output_mod, only: process_output_file
  use file_mod,       only: file_type
  use log_mod,        only: log_event, log_level_error

  implicit none

private

!> @brief Container for file properties need by XIOS
!>
type, public, extends(file_type)  :: lfric_xios_file_type
  private

  !> Unique identifier for XIOS file handle
  character(str_def)          :: xios_id
  !> Path to file
  character(str_max_filename) :: path
  !> Flag denoting if file is to be read in
  logical :: is_input_file
  !> File output frequency
  integer(i_def)              :: output_freq
  !> XIOS ID of associated field group
  character(str_def)          :: field_group = "unset"

contains
  procedure, public :: file_new
  procedure, public :: file_open
  procedure, public :: file_close
  procedure, public :: configure
  procedure, public :: get_xios_id
  procedure, public :: get_path
  procedure, public :: get_output_freq
  procedure, public :: get_field_group

end type lfric_xios_file_type

contains

!> @brief  Initialises an LFRic-XIOS file object
!> @param[in] file_name Path to the file
subroutine file_new(self, file_name)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self
  character(len=*),            intent(in)    :: file_name

  self%path = file_name

end subroutine file_new

!> @brief  Opens the LFRic definition of an output file after it has been
!!         opened by XIOS
!> @param[in] file_name Path to the file
subroutine file_open(self, file_name)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self
  character(len=*),            intent(in)    :: file_name

  ! TODO Write orography

end subroutine file_open

!> @brief  Performs final routines before closure of XIOS file
subroutine file_close(self)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self

  if (.not. self%is_input_file) then
    call process_output_file(trim(self%path)//".nc")
  end if

end subroutine file_close

!>  @brief  Configures a file object for use with XIOS
!>
!>  @param[in]           xios_id         The XIOS name for the file
!>  @param[in,optional]  io_mode_read    Logical flag denoting that the file is
!!                                       to be read into the model
!>  @param[in,optional]  freq            The file output frequency
!>  @param[in,optional]  field_group_id  The associated field group id
!>
subroutine configure(self, xios_id, freq, io_mode_read, field_group_id)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self
  character(len=*),            intent(in)    :: xios_id
  logical,          optional,  intent(in)    :: io_mode_read
  integer(i_def),   optional,  intent(in)    :: freq
  character(len=*), optional,  intent(in)    :: field_group_id

  self%xios_id = xios_id

  ! XIOS default for file mode is "write"
  if (present(io_mode_read)) then
    self%is_input_file = io_mode_read
  else
    self%is_input_file = .false.
  end if

  if (present(freq)) then
    if (freq < 1) then
      call log_event( "XIOS files cannot have negative frequency", &
                      log_level_error )
    end if
    self%output_freq = freq
  else
    ! -999 used as flag for frequency not to be set
    self%output_freq = -999
  end if

  if (present(field_group_id)) self%field_group = field_group_id

end subroutine configure

!> Getter for XIOS file ID
!> @return  output_xios_id  The XIOS ID for the file
function get_xios_id(self) result(output_xios_id)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self

  character(str_def) :: output_xios_id

  output_xios_id = self%xios_id

end function get_xios_id

!> Getter for file path
!> @return output_path The path to the file
function get_path(self) result(output_path)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self

  character(str_max_filename) :: output_path

  output_path = self%path

end function get_path

!> Getter for file output frequency
!> @return  file_freq  The file output frequency
function get_output_freq(self) result(file_freq)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self

  integer(i_def) :: file_freq

  file_freq = self%output_freq

end function get_output_freq

!> Getter for associated field group ID
!> @return  field_group_id  The associated field group id
function get_field_group(self) result(field_group_id)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self

  character(str_def) :: field_group_id

  field_group_id = self%field_group

end function get_field_group

end module lfric_xios_file_mod
