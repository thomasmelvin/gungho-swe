!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module empty_mod

  use constants_mod, only : i_native
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank

  implicit none

  private

  integer(i_native) :: local_rank = -1
  integer(i_native), parameter :: temporary_unit = 3

contains

end module empty_mod
