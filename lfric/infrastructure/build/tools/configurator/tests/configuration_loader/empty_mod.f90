!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the loading of namelists.
!
module empty_mod

  use constants_mod, only : i_native, l_def, str_def, str_max_filename
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank, broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: read_configuration, ensure_configuration, final_configuration

contains

  ! Reads configuration namelists from a file.
  !
  ! [in] filename File holding the namelists.
  !
  ! TODO: Assumes namelist tags come at the start of lines.
  ! TODO: Support "namelist file" namelists which recursively call this
  !       procedure to load other namelist files.
  !
  subroutine read_configuration( filename )

    use io_utility_mod, only : open_file, close_file

    implicit none

    character(*), intent(in) :: filename

    integer(i_native) :: local_rank

    character(str_def), allocatable :: namelists(:)
    integer(i_native) :: unit = -1

    local_rank = get_comm_rank()

    if (local_rank == 0) unit = open_file( filename )

    call get_namelist_names( unit, local_rank, namelists )

    call read_configuration_namelists( unit, local_rank, &
                                       namelists, filename )

    if (local_rank == 0) call close_file( unit )

  end subroutine read_configuration

  ! Finds names of all namelists present in file.
  !
  ! [in] unit File holding namelists.
  ! [out] names of namelist in file (in order).
  !
  subroutine get_namelist_names( unit, local_rank, names )

    use io_utility_mod, only : read_line

    implicit none

    integer(i_native),  intent(in)                 :: unit
    integer(i_native),  intent(in)                 :: local_rank
    character(str_def), intent(inout), allocatable :: names(:)

    character(str_def), allocatable :: names_temp(:)
    ! TODO: Buffer is large enough for a fair sized string and a filename.
    !       Ideally it should be dynamically sized for the length of the
    !       incoming data but I'm not sure how best to achieve that at the
    !       moment. #1752
    character(str_def + str_max_filename) :: buffer
    logical(l_def)     :: continue_read
    ! Number of names - technically a scalar but must be defined as a
    ! single element array to be broadcast-able
    integer(i_native)  :: namecount(1)

    namecount = 0
    if (local_rank == 0) then
      text_line_loop: do

        continue_read = read_line( unit, buffer )
        if ( .not. continue_read ) exit text_line_loop

        ! TODO: Assumes namelist tags are at the start of lines. #1753
        !
        if (buffer(1:1) == '&') then
          namecount = namecount + 1
          allocate(names_temp(namecount(1)))
          if (namecount(1) > 1) then
            names_temp(1:namecount(1)-1) = names
          end if
          names_temp(namecount(1)) = trim(buffer(2:))
          call move_alloc(names_temp, names)
        end if
      end do text_line_loop
      rewind(unit)
    end if

    call broadcast( namecount, 1, 0 )

    if (local_rank /= 0) then
      allocate(names(namecount(1)))
    end if

    call broadcast( names, namecount(1)*str_def, 0 )

  end subroutine get_namelist_names

  ! Checks that the requested namelists have been loaded.
  !
  ! [in]  names List of namelists.
  ! [out] success_mask Marks corresponding namelists as having failed.
  !
  ! [return] Overall success.
  !
  function ensure_configuration( names, success_mask )

    implicit none

    character(*),             intent(in)  :: names(:)
    logical(l_def), optional, intent(out) :: success_mask(:)
    logical(l_def)                        :: ensure_configuration

    integer(i_native) :: i
    logical           :: configuration_found = .True.

    if (present(success_mask) &
        .and. (size(success_mask, 1) /= size(names, 1))) then
      call log_event( 'Arguments "names" and "success_mask" to function' &
                      // '"ensure_configuration" are different shapes',  &
                      LOG_LEVEL_ERROR )
    end if

    ensure_configuration = .True.

    name_loop: do i = 1, size(names)
      select case(trim( names(i) ))
        case default
          write( log_scratch_space, '(A, A, A)' )          &
               "Tried to ensure unrecognised namelist """, &
               trim(names(i)),                             &
               """ was loaded"
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

      ensure_configuration = ensure_configuration .and. configuration_found

      if (present(success_mask)) success_mask(i) = configuration_found

    end do name_loop

  end function ensure_configuration

  subroutine read_configuration_namelists( unit, local_rank, &
                                           namelists, filename )
    implicit none

    integer(i_native),  intent(in) :: unit
    integer(i_native),  intent(in) :: local_rank
    character(str_def), intent(in) :: namelists(:)
    character(*),       intent(in) :: filename

    integer(i_native) :: i

    ! Read the namelists
    do i = 1, size(namelists)
      select case (trim(namelists(i)))
        case default
          write( log_scratch_space, '(A)' )        &
                "Unrecognised namelist """//        &
                trim(namelists(i))//                &
                """ found in file "//               &
                trim(filename)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do

    ! Perform post load actions
    do i = 1, size(namelists)
      select case (trim(namelists(i)))
        case default
          write( log_scratch_space, '(A)' )        &
                "Unrecognised namelist """//        &
                trim(namelists(i))//                &
                """ found in file "//               &
                trim(filename)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do

  end subroutine read_configuration_namelists

  subroutine final_configuration()

    implicit none

    return
  end subroutine final_configuration

end module empty_mod
