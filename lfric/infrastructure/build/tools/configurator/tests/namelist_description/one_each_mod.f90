!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only: i_def, &
                           i_long, &
                           i_native, &
                           i_short, &
                           l_def, &
                           r_def, &
                           r_double, &
                           r_second, &
                           r_single, &
                           str_def, &
                           str_max_filename
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: enum_from_key, key_from_enum, &
            read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, test_final

  integer(i_native), public, parameter :: enum_one = 189779348
  integer(i_native), public, parameter :: enum_three = 1061269036
  integer(i_native), public, parameter :: enum_two = 1625932035

  integer(i_def), public, protected :: dint = imdi
  logical(l_def), public, protected :: dlog = .false.
  real(r_def), public, protected :: dreal = rmdi
  character(str_def), public, protected :: dstr = cmdi
  integer(i_native), public, protected :: enum = emdi
  character(str_max_filename), public, protected :: fstr = cmdi
  integer(i_long), public, protected :: lint = imdi
  real(r_double), public, protected :: lreal = rmdi
  integer(i_short), public, protected :: sint = imdi
  real(r_single), public, protected :: sreal = rmdi
  real(r_second), public, protected :: treal = rmdi
  integer(i_def), public, protected :: vint = imdi
  real(r_def), public, protected :: vreal = rmdi
  character(str_def), public, protected :: vstr = cmdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: enum_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_native), parameter :: enum_value(3) &
          = [189779348_i_native, &
             1061269036_i_native, &
             1625932035_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function enum_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for enum enumeration in test namelist.'
      enum_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      return
    end if

    key_index = 1
    do
      if (trim(enum_key(key_index)) == trim(key)) then
        enum_from_key = enum_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(enum_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for test enum")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function enum_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_enum( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (enum_value(value_index) == emdi) then
        key_from_enum = unset_key
        return
      else if (enum_value(value_index) == value) then
        key_from_enum = enum_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(enum_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in test enum")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_enum

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_test_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        enum )

  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_enum )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data
    integer(i_native), intent(out) :: dummy_enum

    character(str_def) :: buffer_character_str_def(2)
    character(str_max_filename) :: buffer_character_str_max_filename(1)
    integer(i_def) :: buffer_integer_i_def(2)
    integer(i_long) :: buffer_integer_i_long(1)
    integer(i_native) :: buffer_integer_i_native(1)
    integer(i_short) :: buffer_integer_i_short(1)
    integer(i_native) :: buffer_logical_l_def(1)
    real(r_def) :: buffer_real_r_def(2)
    real(r_double) :: buffer_real_r_double(1)
    real(r_second) :: buffer_real_r_second(1)
    real(r_single) :: buffer_real_r_single(1)

    character(str_def) :: enum

    namelist /test/ dint, &
                    dlog, &
                    dreal, &
                    dstr, &
                    enum, &
                    fstr, &
                    lint, &
                    lreal, &
                    sint, &
                    sreal, &
                    treal, &
                    vint, &
                    vreal, &
                    vstr

    integer(i_native) :: condition

    missing_data = 0

    dint = imdi
    dlog = .false.
    dreal = rmdi
    dstr = cmdi
    enum = unset_key
    fstr = cmdi
    lint = imdi
    lreal = rmdi
    sint = imdi
    sreal = rmdi
    treal = rmdi
    vint = imdi
    vreal = rmdi
    vstr = cmdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_enum = enum_from_key( enum )

    end if

    buffer_integer_i_def(2) = dint
    buffer_logical_l_def(1) = merge( 1, 0, dlog )
    buffer_real_r_def(2) = dreal
    buffer_character_str_def(2) = dstr
    buffer_integer_i_native(1) = dummy_enum
    buffer_character_str_max_filename(1) = fstr
    buffer_integer_i_long(1) = lint
    buffer_real_r_double(1) = lreal
    buffer_integer_i_short(1) = sint
    buffer_real_r_single(1) = sreal
    buffer_real_r_second(1) = treal
    buffer_integer_i_def(1) = vint
    buffer_real_r_def(1) = vreal
    buffer_character_str_def(1) = vstr

    call broadcast( buffer_character_str_def, 2*str_def, 0 )
    call broadcast( buffer_character_str_max_filename, 1*str_max_filename, 0 )
    call broadcast( buffer_integer_i_def, 2, 0 )
    call broadcast( buffer_integer_i_long, 1, 0 )
    call broadcast( buffer_integer_i_native, 1, 0 )
    call broadcast( buffer_integer_i_short, 1, 0 )
    call broadcast( buffer_logical_l_def, 1, 0 )
    call broadcast( buffer_real_r_def, 2, 0 )
    call broadcast( buffer_real_r_double, 1, 0 )
    call broadcast( buffer_real_r_second, 1, 0 )
    call broadcast( buffer_real_r_single, 1, 0 )

    dint = buffer_integer_i_def(2)
    dlog = buffer_logical_l_def(1) /= 0
    dreal = buffer_real_r_def(2)
    dstr = buffer_character_str_def(2)
    dummy_enum = buffer_integer_i_native(1)
    fstr = buffer_character_str_max_filename(1)
    lint = buffer_integer_i_long(1)
    lreal = buffer_real_r_double(1)
    sint = buffer_integer_i_short(1)
    sreal = buffer_real_r_single(1)
    treal = buffer_real_r_second(1)
    vint = buffer_integer_i_def(1)
    vreal = buffer_real_r_def(1)
    vstr = buffer_character_str_def(1)


    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    implicit none


  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine test_final()

    implicit none

    dint = imdi
    dlog = .false.
    dreal = real(rmdi,r_def)
    dstr = cmdi
    enum = emdi
    fstr = cmdi
    lint = imdi
    lreal = real(rmdi,r_double)
    sint = imdi
    sreal = real(rmdi,r_single)
    treal = real(rmdi,r_second)
    vint = imdi
    vreal = real(rmdi,r_def)
    vstr = cmdi

    return
  end subroutine test_final


end module test_config_mod
