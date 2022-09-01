!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the teapot namelist.
!>
module teapot_config_mod

  use constants_mod, only: i_native, &
                           r_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: read_teapot_namelist, postprocess_teapot_namelist, &
            teapot_is_loadable, teapot_is_loaded, teapot_final

  real(r_def), public, protected :: bar = rmdi
  real(r_def), public, protected :: foo = rmdi
  real(r_def), public, protected :: fum = rmdi

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_teapot_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_teapot_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data

    real(r_def) :: buffer_real_r_def(2)

    namelist /teapot/ foo, &
                      fum

    integer(i_native) :: condition

    missing_data = 0

    bar = rmdi
    foo = rmdi
    fum = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=teapot, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = foo
    buffer_real_r_def(2) = fum

    call broadcast( buffer_real_r_def, 2, 0 )

    foo = buffer_real_r_def(1)
    fum = buffer_real_r_def(2)

   ! Parameter name bar: derived by computation
    bar = foo ** 2

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_teapot_namelist()

    implicit none


  end subroutine postprocess_teapot_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function teapot_is_loadable()

    implicit none

    logical :: teapot_is_loadable

    teapot_is_loadable = .not. namelist_loaded

  end function teapot_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function teapot_is_loaded()

    implicit none

    logical :: teapot_is_loaded

    teapot_is_loaded = namelist_loaded

  end function teapot_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine teapot_final()

    implicit none

    bar = real(rmdi,r_def)
    foo = real(rmdi,r_def)
    fum = real(rmdi,r_def)

    return
  end subroutine teapot_final


end module teapot_config_mod
