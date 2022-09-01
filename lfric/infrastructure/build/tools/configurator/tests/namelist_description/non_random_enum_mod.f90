!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the telly namelist.
!>
module telly_config_mod

  use constants_mod, only: i_native, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: tubbies_from_key, key_from_tubbies, &
            read_telly_namelist, postprocess_telly_namelist, &
            telly_is_loadable, telly_is_loaded, telly_final

  integer(i_native), public, parameter :: tubbies_inky = 3
  integer(i_native), public, parameter :: tubbies_lala = 1
  integer(i_native), public, parameter :: tubbies_po = 2

  integer(i_native), public, protected :: tubbies = emdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: tubbies_key(3) &
          = [character(len=str_def) :: 'inky', &
                                       'lala', &
                                       'po']

  integer(i_native), parameter :: tubbies_value(3) &
          = [3_i_native, &
             1_i_native, &
             2_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function tubbies_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for tubbies enumeration in telly namelist.'
      tubbies_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      return
    end if

    key_index = 1
    do
      if (trim(tubbies_key(key_index)) == trim(key)) then
        tubbies_from_key = tubbies_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(tubbies_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for telly tubbies")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function tubbies_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_tubbies( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (tubbies_value(value_index) == emdi) then
        key_from_tubbies = unset_key
        return
      else if (tubbies_value(value_index) == value) then
        key_from_tubbies = tubbies_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(tubbies_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in telly tubbies")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_tubbies

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_telly_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        tubbies )

  end subroutine read_telly_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_tubbies )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data
    integer(i_native), intent(out) :: dummy_tubbies

    integer(i_native) :: buffer_integer_i_native(1)

    character(str_def) :: tubbies

    namelist /telly/ tubbies

    integer(i_native) :: condition

    missing_data = 0

    tubbies = unset_key

    if (local_rank == 0) then

      read( file_unit, nml=telly, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_tubbies = tubbies_from_key( tubbies )

    end if

    buffer_integer_i_native(1) = dummy_tubbies

    call broadcast( buffer_integer_i_native, 1, 0 )

    dummy_tubbies = buffer_integer_i_native(1)


    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_telly_namelist()

    implicit none


  end subroutine postprocess_telly_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function telly_is_loadable()

    implicit none

    logical :: telly_is_loadable

    telly_is_loadable = .not. namelist_loaded

  end function telly_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function telly_is_loaded()

    implicit none

    logical :: telly_is_loaded

    telly_is_loaded = namelist_loaded

  end function telly_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine telly_final()

    implicit none

    tubbies = emdi

    return
  end subroutine telly_final


end module telly_config_mod
