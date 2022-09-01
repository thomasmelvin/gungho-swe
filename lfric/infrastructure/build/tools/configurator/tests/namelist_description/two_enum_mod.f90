!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the twoenum namelist.
!>
module twoenum_config_mod

  use constants_mod, only: i_native, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: first_from_key, key_from_first, &
            second_from_key, key_from_second, &
            read_twoenum_namelist, postprocess_twoenum_namelist, &
            twoenum_is_loadable, twoenum_is_loaded, twoenum_final

  integer(i_native), public, parameter :: first_one = 1952457118
  integer(i_native), public, parameter :: first_three = 1813125082
  integer(i_native), public, parameter :: first_two = 533081353
  integer(i_native), public, parameter :: second_ay = 1248446338
  integer(i_native), public, parameter :: second_bee = 144118421
  integer(i_native), public, parameter :: second_see = 359914450

  integer(i_native), public, protected :: first = emdi
  integer(i_native), public, protected :: second = emdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: first_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']
  character(str_def), parameter :: second_key(3) &
          = [character(len=str_def) :: 'ay', &
                                       'bee', &
                                       'see']

  integer(i_native), parameter :: first_value(3) &
          = [1952457118_i_native, &
             1813125082_i_native, &
             533081353_i_native]
  integer(i_native), parameter :: second_value(3) &
          = [1248446338_i_native, &
             144118421_i_native, &
             359914450_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function first_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for first enumeration in twoenum namelist.'
      first_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      return
    end if

    key_index = 1
    do
      if (trim(first_key(key_index)) == trim(key)) then
        first_from_key = first_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(first_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for twoenum first")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function first_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_first( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (first_value(value_index) == emdi) then
        key_from_first = unset_key
        return
      else if (first_value(value_index) == value) then
        key_from_first = first_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(first_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in twoenum first")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_first

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function second_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for second enumeration in twoenum namelist.'
      second_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      return
    end if

    key_index = 1
    do
      if (trim(second_key(key_index)) == trim(key)) then
        second_from_key = second_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(second_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for twoenum second")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function second_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_second( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (second_value(value_index) == emdi) then
        key_from_second = unset_key
        return
      else if (second_value(value_index) == value) then
        key_from_second = second_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(second_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in twoenum second")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_second

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_twoenum_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        first, &
                        second )

  end subroutine read_twoenum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_first, &
                            dummy_second )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data
    integer(i_native), intent(out) :: dummy_first
    integer(i_native), intent(out) :: dummy_second

    integer(i_native) :: buffer_integer_i_native(2)

    character(str_def) :: first
    character(str_def) :: second

    namelist /twoenum/ first, &
                       second

    integer(i_native) :: condition

    missing_data = 0

    first = unset_key
    second = unset_key

    if (local_rank == 0) then

      read( file_unit, nml=twoenum, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_first = first_from_key( first )
      dummy_second = second_from_key( second )

    end if

    buffer_integer_i_native(1) = dummy_first
    buffer_integer_i_native(2) = dummy_second

    call broadcast( buffer_integer_i_native, 2, 0 )

    dummy_first = buffer_integer_i_native(1)
    dummy_second = buffer_integer_i_native(2)


    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_twoenum_namelist()

    implicit none


  end subroutine postprocess_twoenum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function twoenum_is_loadable()

    implicit none

    logical :: twoenum_is_loadable

    twoenum_is_loadable = .not. namelist_loaded

  end function twoenum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function twoenum_is_loaded()

    implicit none

    logical :: twoenum_is_loaded

    twoenum_is_loaded = namelist_loaded

  end function twoenum_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine twoenum_final()

    implicit none

    first = emdi
    second = emdi

    return
  end subroutine twoenum_final


end module twoenum_config_mod
