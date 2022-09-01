!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief  File handler for NetCDF files.
module lfric_ncdf_file_mod

  use constants_mod, only: i_def, str_max_filename
  use log_mod,       only: log_event, log_scratch_space, LOG_LEVEL_ERROR
  use netcdf,        only: nf90_open, nf90_create, nf90_close,           &
                           nf90_write, nf90_nowrite, nf90_clobber,       &
                           nf90_strerror, nf90_noerr, nf90_64bit_offset, &
                           nf90_put_att, nf90_global, nf90_enddef

  implicit none

  private

  !> @brief  NetCDF file type.
  type, public :: lfric_ncdf_file_type

    private

    integer(kind=i_def)             :: ncid
    character(len=str_max_filename) :: name
    integer(kind=i_def)             :: mode

  contains

    procedure, private :: create_file
    procedure, private :: open_file
    procedure, public  :: close_file
    procedure, public  :: get_id
    procedure, public  :: get_io_mode
    procedure, public  :: set_attribute
    procedure, public  :: close_definition

  end type

  interface lfric_ncdf_file_type
    module procedure lfric_ncdf_file_constructor
  end interface

  ! IO Mode enumerations
  integer(kind=i_def), public, parameter :: LFRIC_NCDF_OPEN   = 485
  integer(kind=i_def), public, parameter :: LFRIC_NCDF_CREATE = 653
  integer(kind=i_def), public, parameter :: LFRIC_NCDF_READ   = 971
  integer(kind=i_def), public, parameter :: LFRIC_NCDF_WRITE  = 248

contains

  !> @brief  Constructor for an LFRic NetCDF file object.
  !>
  !> @param[in]  name       The name of the NetCDF file
  !> @param[in]  open_mode  The opening mode of the NetCDF file
  !> @param[in]  io_mode    The IO mode of the NetCDF file
  function lfric_ncdf_file_constructor(name, open_mode, io_mode) result(self)

    implicit none

    type(lfric_ncdf_file_type)             :: self
    character(len=*),           intent(in) :: name
    integer(kind=i_def),        intent(in) :: open_mode
    integer(kind=i_def),        intent(in) :: io_mode

    self%name = trim(name)

    select case (io_mode)
      case (LFRIC_NCDF_READ)
        self%mode = nf90_nowrite

      case (LFRIC_NCDF_WRITE)
        self%mode = nf90_write

      case default
        call log_event( "[lfric_ncdf_file_mod] - Invalid option for ncdf " // &
                        "file IO mode in lfric_ncdf_file_type constructor",   &
                        LOG_LEVEL_ERROR )
    end select

    select case (open_mode)
    case (LFRIC_NCDF_OPEN)
      call self%open_file()

    case (LFRIC_NCDF_CREATE)
      call self%create_file()

    case default
      call log_event( "[lfric_ncdf_file_mod] - Invalid option for ncdf " //    &
                      "file opening mode in lfric_ncdf_file_type constructor", &
                      LOG_LEVEL_ERROR )
    end select

    return

  end function lfric_ncdf_file_constructor

  !> @brief  Create a new NetCDF file.
  !> @details Creates an opens a new, clean NetCDF file. If a file of the
  !!          same name already exists, this routine will clobber it.
  subroutine create_file(self)

    implicit none

    class(lfric_ncdf_file_type), intent(inout) :: self

    integer(kind=i_def)         :: ierr
    character(len=*), parameter :: routine = 'create_file'

    ierr = nf90_create( path=trim(self%name),                      &
                        cmode=ior(nf90_clobber,nf90_64bit_offset), &
                        ncid=self%ncid )

    call check_err(ierr, routine, self%name)

    return

  end subroutine create_file

  !> @brief  Open a NetCDF file for IO.
  subroutine open_file(self)

    implicit none

    class(lfric_ncdf_file_type), intent(inout) :: self

    integer(kind=i_def)         :: ierr
    character(len=*), parameter :: routine = 'open_file'

    ierr = nf90_open( trim(self%name), self%mode, self%ncid )

    call check_err(ierr, routine, self%name)

    return

  end subroutine open_file

  !> @brief  Closes the NetCDF file.
  subroutine close_file(self)

    implicit none

    class(lfric_ncdf_file_type), intent(inout) :: self

    integer(kind=i_def)         :: ierr
    character(len=*), parameter :: routine = 'close_file'

    ierr = nf90_close( self%ncid )

    call check_err(ierr, routine, self%name)

    return

  end subroutine close_file

  !> @brief  Returns the ncid of the file.
  !>
  !> @return  The NetCDF file ID
  function get_id(self) result(id)

    implicit none

    class(lfric_ncdf_file_type), intent(in) :: self

    integer(kind=i_def) :: id

    id = self%ncid

    return

  end function get_id

  !> @brief  Returns the file's nf90 io_mode enumerator.
  !>
  !> @return  NetCDF IO mode enumerator
  function get_io_mode(self) result(io_mode)

    implicit none

    class(lfric_ncdf_file_type), intent(in) :: self

    integer(kind=i_def) :: io_mode

    io_mode = self%mode

    return

  end function get_io_mode

  !> @brief    Assigns global attributes to the NetCDF file.
  !>
  !> @param[in]  attr_name   The name of attribute to be set
  !> @param[in]  attr_value  The value of the attribute
  subroutine set_attribute(self, attr_name, attr_value)

    implicit none

    class(lfric_ncdf_file_type), intent(in) :: self
    character(len=*),            intent(in) :: attr_name
    character(len=*),            intent(in) :: attr_value

    integer(kind=i_def)         :: ierr
    character(len=*), parameter :: routine = 'set_attribute'

    ierr = nf90_put_att(self%ncid, nf90_global, attr_name, &
                        trim(attr_value))

    call check_err(ierr, routine, self%name)

    return

  end subroutine set_attribute

  !> @brief  Closes the NetCDF file definition.
  subroutine close_definition(self)

    implicit none

    class(lfric_ncdf_file_type), intent(inout) :: self

    integer(kind=i_def)         :: ierr
    character(len=*), parameter :: routine = 'close_definition'

    ierr = nf90_enddef(self%ncid)

    call check_err(ierr, routine, self%name)

    return

  end subroutine close_definition

  !> @brief    Calls logger on error.
  !> @details  Checks the error code returned by the NetCDF file. If an error is
  !!           detected, the relevant error message is passed to the logger.
  !>
  !> @param[in] ierr     The error code to check
  !> @param[in] routine  The routine name that call the error check
  !> @param[in] filename The file name that raised the error
  subroutine check_err(ierr, routine, filename)

    implicit none

    integer(kind=i_def), intent(in) :: ierr
    character(len=*),    intent(in) :: routine, filename

    if ( ierr /= nf90_noerr ) then
      write(log_scratch_space,*) "Error in lfric_ncdf_file_mod ['"//routine// &
                                 "'] for file '"//trim(filename)//"': "//     &
                                 trim(nf90_strerror(ierr))
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
    end if

    return

  end subroutine check_err

end module lfric_ncdf_file_mod