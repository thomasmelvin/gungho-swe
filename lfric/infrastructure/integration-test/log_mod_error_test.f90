!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

! A very simply program which just logs an error.
!
program log_mod_error_test

  use, intrinsic :: iso_fortran_env, only : error_unit
  use log_mod, only : initialise_logging, finalise_logging, log_event, &
                      LOG_LEVEL_ERROR
  use mpi_mod, only : initialise_comm, store_comm, &
                      finalise_comm, &
                      get_comm_size, get_comm_rank

  implicit none

  integer :: total_ranks, local_rank, comm

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm(comm)
  ! Store the communicator for later use
  call store_comm(comm)
  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()
  call initialise_logging(local_rank, total_ranks, 'log_mod_error_test')

  ! Testing an error occurring on rank 0 of a parallel
  ! application. The choice of rank needs to be consistent with the
  ! matching Python routine which should query the output from the
  ! first rank: the PET0* file.

  ! Note that we initiate the error on the first rank only. Initiating
  ! an error on all ranks will sometimes result in the first rank
  ! being killed off by the error report in another rank before it has
  ! processed the error. If the first rank is prevented from reporting
  ! an error, the test of the output will fail to find the expected
  ! error message.

  if (local_rank == 0) then
    call log_event( 'An error was logged.', LOG_LEVEL_ERROR )
  end if

  ! Finalise mpi and release the communicator
  call finalise_comm()

  ! Finalise the logging system
  call finalise_logging()

end program log_mod_error_test
