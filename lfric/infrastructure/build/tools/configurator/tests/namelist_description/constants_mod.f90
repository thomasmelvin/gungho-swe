!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the cheese namelist.
!>
module cheese_config_mod

  use constants_mod, only: i_native, &
                           r_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, FUDGE, imdi, rmdi, unset_key

  implicit none

  private
  public :: read_cheese_namelist, postprocess_cheese_namelist, &
            cheese_is_loadable, cheese_is_loaded, cheese_final

  real(r_def), public, protected :: fred = rmdi
  real(r_def), public, protected :: wilma = rmdi

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_cheese_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_cheese_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only: i_def

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_def)                :: missing_data

    real(r_def) :: buffer_real_r_def(1)

    namelist /cheese/ fred

    integer(i_native) :: condition

    missing_data = 0

    fred = rmdi
    wilma = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=cheese, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = fred

    call broadcast( buffer_real_r_def, 1, 0 )

    fred = buffer_real_r_def(1)

   ! Parameter name wilma: derived by computation
    wilma = fred * FUDGE

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_cheese_namelist()

    implicit none


  end subroutine postprocess_cheese_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function cheese_is_loadable()

    implicit none

    logical :: cheese_is_loadable

    cheese_is_loadable = .not. namelist_loaded

  end function cheese_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function cheese_is_loaded()

    implicit none

    logical :: cheese_is_loaded

    cheese_is_loaded = namelist_loaded

  end function cheese_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine cheese_final()

    implicit none

    fred = real(rmdi,r_def)
    wilma = real(rmdi,r_def)

    return
  end subroutine cheese_final


end module cheese_config_mod
