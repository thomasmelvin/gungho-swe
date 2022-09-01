!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the enum namelist.
!>
module enum_config_mod

  use constants_mod, only: i_native, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: value_from_key, key_from_value, &
            read_enum_namelist, postprocess_enum_namelist, &
            enum_is_loadable, enum_is_loaded, enum_final

  integer(i_native), public, parameter :: value_one = 1695414371
  integer(i_native), public, parameter :: value_three = 839906103
  integer(i_native), public, parameter :: value_two = 246150388

  integer(i_native), public, protected :: value = emdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: value_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_native), parameter :: value_value(3) &
          = [1695414371_i_native, &
             839906103_i_native, &
             246150388_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function value_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for value enumeration in enum namelist.'
      value_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      return
    end if

    key_index = 1
    do
      if (trim(value_key(key_index)) == trim(key)) then
        value_from_key = value_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(value_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for enum value")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function value_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_value( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (value_value(value_index) == emdi) then
        key_from_value = unset_key
        return
      else if (value_value(value_index) == value) then
        key_from_value = value_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(value_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in enum value")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_value

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_enum_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        value )

  end subroutine read_enum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_value )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data
    integer(i_native), intent(out) :: dummy_value

    integer(i_native) :: buffer_integer_i_native(1)

    character(str_def) :: value

    namelist /enum/ value

    integer(i_native) :: condition

    missing_data = 0

    value = unset_key

    if (local_rank == 0) then

      read( file_unit, nml=enum, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_value = value_from_key( value )

    end if

    buffer_integer_i_native(1) = dummy_value

    call broadcast( buffer_integer_i_native, 1, 0 )

    dummy_value = buffer_integer_i_native(1)


    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_enum_namelist()

    implicit none


  end subroutine postprocess_enum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function enum_is_loadable()

    implicit none

    logical :: enum_is_loadable

    enum_is_loadable = .not. namelist_loaded

  end function enum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function enum_is_loaded()

    implicit none

    logical :: enum_is_loaded

    enum_is_loaded = namelist_loaded

  end function enum_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine enum_final()

    implicit none

    value = emdi

    return
  end subroutine enum_final


end module enum_config_mod
