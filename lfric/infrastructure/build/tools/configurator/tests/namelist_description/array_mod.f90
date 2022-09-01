!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the aerial namelist.
!>
module aerial_config_mod

  use constants_mod, only: i_def, &
                           i_native, &
                           r_def, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key
  use wibble_mod, only: esize

  implicit none

  private
  public :: read_aerial_namelist, postprocess_aerial_namelist, &
            aerial_is_loadable, aerial_is_loaded, aerial_final

  integer(i_native), parameter, public :: max_array_size = 100

  character(str_def), public, protected :: absolute(5) = cmdi
  integer(i_def), public, protected, allocatable :: inlist(:)
  integer(i_native), public, protected :: lsize = imdi
  real(r_def), public, protected, allocatable :: outlist(:)
  integer(i_def), public, protected, allocatable :: unknown(:)

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_aerial_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_aerial_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data

    integer(i_native) :: buffer_integer_i_native(1)

    namelist /aerial/ absolute, &
                      inlist, &
                      lsize, &
                      outlist, &
                      unknown

    integer(i_native) :: condition

    missing_data = 0

    allocate( inlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    allocate( outlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    allocate( unknown(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    absolute = cmdi
    inlist = imdi
    lsize = imdi
    outlist = rmdi
    unknown = imdi

    if (local_rank == 0) then

      read( file_unit, nml=aerial, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_integer_i_native(1) = lsize

    call broadcast( buffer_integer_i_native, 1, 0 )

    lsize = buffer_integer_i_native(1)


    call broadcast( absolute, size(absolute, 1)*str_def, 0 )
    call broadcast( inlist, size(inlist, 1), 0 )
    call broadcast( outlist, size(outlist, 1), 0 )
    call broadcast( unknown, size(unknown, 1), 0 )

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_aerial_namelist()

    implicit none

    integer(i_native) :: condition
    integer(i_native) :: array_size
    integer(i_def), allocatable :: new_inlist(:)
    real(r_def), allocatable :: new_outlist(:)
    integer(i_def), allocatable :: new_unknown(:)
    integer(i_native) :: index_unknown

    condition  = 0
    array_size = 0

    array_size = lsize
    allocate( new_inlist(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_inlist(:array_size) = inlist(:array_size)
    call move_alloc( new_inlist, inlist )
    if (allocated(new_inlist)) deallocate( new_inlist)
    array_size = esize
    allocate( new_outlist(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_outlist(:array_size) = outlist(:array_size)
    call move_alloc( new_outlist, outlist )
    if (allocated(new_outlist)) deallocate( new_outlist)
    do index_unknown=ubound(unknown, 1), 1, -1
      if (unknown(index_unknown) /= imdi) exit
    end do
    array_size = index_unknown
    allocate( new_unknown(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_unknown(:array_size) = unknown(:array_size)
    call move_alloc( new_unknown, unknown )
    if (allocated(new_unknown)) deallocate( new_unknown)

  end subroutine postprocess_aerial_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function aerial_is_loadable()

    implicit none

    logical :: aerial_is_loadable

    aerial_is_loadable = .not. namelist_loaded

  end function aerial_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function aerial_is_loaded()

    implicit none

    logical :: aerial_is_loaded

    aerial_is_loaded = namelist_loaded

  end function aerial_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine aerial_final()

    implicit none

    absolute = cmdi
    lsize = imdi

    if ( allocated(inlist) ) deallocate(inlist)
    if ( allocated(outlist) ) deallocate(outlist)
    if ( allocated(unknown) ) deallocate(unknown)

    return
  end subroutine aerial_final


end module aerial_config_mod
