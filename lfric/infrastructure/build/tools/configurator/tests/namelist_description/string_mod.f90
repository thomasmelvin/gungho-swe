!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the mirth namelist.
!>
module mirth_config_mod

  use constants_mod, only: i_native, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key
  use random_config_mod, only: biggles

  implicit none

  private
  public :: read_mirth_namelist, postprocess_mirth_namelist, &
            mirth_is_loadable, mirth_is_loaded, mirth_final

  integer(i_native), parameter, public :: max_array_size = 100

  character(str_def), public, protected, allocatable :: chortle(:)
  character(str_def), public, protected :: chuckle = cmdi
  character(str_def), public, protected :: guffaw(3) = cmdi
  character(str_def), public, protected, allocatable :: hysterics(:)

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_mirth_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_mirth_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data

    character(str_def) :: buffer_character_str_def(1)

    namelist /mirth/ chortle, &
                     chuckle, &
                     guffaw, &
                     hysterics

    integer(i_native) :: condition

    missing_data = 0

    allocate( hysterics(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "hysterics"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    allocate( chortle(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "chortle"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    chortle = cmdi
    chuckle = cmdi
    guffaw = cmdi
    hysterics = cmdi

    if (local_rank == 0) then

      read( file_unit, nml=mirth, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_character_str_def(1) = chuckle

    call broadcast( buffer_character_str_def, 1*str_def, 0 )

    chuckle = buffer_character_str_def(1)


    call broadcast( chortle, size(chortle, 1)*str_def, 0 )
    call broadcast( guffaw, size(guffaw, 1)*str_def, 0 )
    call broadcast( hysterics, size(hysterics, 1)*str_def, 0 )

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_mirth_namelist()

    implicit none

    integer(i_native) :: condition
    integer(i_native) :: array_size
    character(str_def), allocatable :: new_hysterics(:)
    integer(i_native) :: index_hysterics
    character(str_def), allocatable :: new_chortle(:)

    condition  = 0
    array_size = 0

    do index_hysterics=ubound(hysterics, 1), 1, -1
      if (hysterics(index_hysterics) /= cmdi) exit
    end do
    array_size = index_hysterics
    allocate( new_hysterics(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "hysterics"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_hysterics(:array_size) = hysterics(:array_size)
    call move_alloc( new_hysterics, hysterics )
    if (allocated(new_hysterics)) deallocate( new_hysterics)
    array_size = biggles
    allocate( new_chortle(array_size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "chortle"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_chortle(:array_size) = chortle(:array_size)
    call move_alloc( new_chortle, chortle )
    if (allocated(new_chortle)) deallocate( new_chortle)

  end subroutine postprocess_mirth_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function mirth_is_loadable()

    implicit none

    logical :: mirth_is_loadable

    mirth_is_loadable = .not. namelist_loaded

  end function mirth_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function mirth_is_loaded()

    implicit none

    logical :: mirth_is_loaded

    mirth_is_loaded = namelist_loaded

  end function mirth_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine mirth_final()

    implicit none

    chuckle = cmdi
    guffaw = cmdi

    if ( allocated(chortle) ) deallocate(chortle)
    if ( allocated(hysterics) ) deallocate(hysterics)

    return
  end subroutine mirth_final


end module mirth_config_mod
