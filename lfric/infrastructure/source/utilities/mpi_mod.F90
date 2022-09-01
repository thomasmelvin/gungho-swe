!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides access to the MPI related functionality
!>
!> Provides access to global reduction functions, all_gather, broadcasts and
!. generation of the halo exchange redistribution object. In order for that
!> functionality to work, the subroutine store_comm must first be called to
!> store a valid MPI communicator.
!>
module mpi_mod

  use constants_mod, only : i_def, i_halo_index, i_native,              &
                            l_def, r_double, r_single, str_def,         &
                            real_type, integer_type, logical_type
#ifdef NO_MPI
  ! No "use mpi" in non-mpi build
#else
  use mpi,           only : mpi_comm_rank, mpi_comm_size, mpi_finalize, &
                            mpi_init, mpi_success, mpi_comm_world,      &
                            mpi_max, mpi_min, mpi_sum,                  &
                            mpi_character, mpi_double_precision,        &
                            mpi_integer, mpi_integer1, mpi_integer2,    &
                            mpi_integer8, mpi_logical, mpi_real4
#endif
  use log_mod,       only : log_event, LOG_LEVEL_ERROR

  implicit none

  private
  public initialise_comm, finalise_comm, store_comm, clear_comm, &
         get_comm, is_comm_set, &
         global_sum, global_min, global_max, &
         all_gather, broadcast, &
         get_mpi_datatype, &
         get_comm_size, get_comm_rank

  ! The mpi communicator
  integer(i_def), private :: comm=-999, comm_size=-999, comm_rank=-999
  ! Flag marks whether an MPI communicator has been stored
  logical(l_def), private :: comm_set = .false.

  ! Generic interface for specific broadcast functions
  interface broadcast
   module procedure broadcast_l_def,    &
                    broadcast_i_def,    &
                    broadcast_r_double, &
                    broadcast_r_single, &
                    broadcast_str
  end interface

  ! Generic interface for specific global_sum functions
  interface global_sum
   module procedure global_sum_i_def,    &
                    global_sum_r_double, &
                    global_sum_r_single
  end interface

  ! Generic interface for specific max functions
  interface global_max
   module procedure global_max_i_def,    &
                    global_max_r_double, &
                    global_max_r_single
  end interface

  ! Generic interface for specific min functions
  interface global_min
   module procedure global_min_i_def,    &
                    global_min_r_double, &
                    global_min_r_single
  end interface

contains

  !> Initialises MPI and returns mpi_comm_world as the communicator
  !>
  !> @param out_comm The MPI communicator.
  !>
  subroutine initialise_comm(out_comm)
    implicit none
    integer(i_native), intent(out) :: out_comm
    integer(i_native) :: ierr

#ifdef NO_MPI
    ! Don't initialise mpi in non-mpi build.
    out_comm = 0
#else
    call mpi_init(ierr)
    if (ierr /= mpi_success) &
          call log_event('Unable to initialise MPI', LOG_LEVEL_ERROR )
    out_comm = mpi_comm_world
#endif
  end subroutine initialise_comm

  !> Stores the MPI communicator in a private variable, ready for later use.
  !>
  !> @param in_comm The MPI communicator to be stored.
  !>
  subroutine store_comm(in_comm)
    implicit none
    integer(i_def), intent(in) :: in_comm
    integer(i_def) :: ierr

    comm = in_comm
#ifdef NO_MPI
    ! Set default values for number of ranks and local rank in non-mpi build
    comm_size = 1
    comm_rank = 0
#else
    call mpi_comm_size(comm,comm_size,ierr)
    call mpi_comm_rank(comm,comm_rank,ierr)
#endif
    comm_set = .true.
  end subroutine store_comm

  !> Finalises MPI
  !>
  subroutine finalise_comm()
    implicit none
    integer(i_def) :: ierr

#ifdef NO_MPI
    ! Don't finalise mpi in non-mpi build
#else
    call mpi_finalize(ierr)
    if (ierr /= mpi_success) &
          call log_event('Unable to finalise MPI', LOG_LEVEL_ERROR )
#endif
    comm = -999
    comm_size = -999
    comm_rank = -999
    comm_set = .false.
  end subroutine finalise_comm

  !> Clears the stored MPI communicator
  !>
  subroutine clear_comm()
    implicit none
    comm = -999
    comm_size = -999
    comm_rank = -999
    comm_set = .false.
  end subroutine clear_comm

  !> Returns the stored MPI communicator
  !> @return communicator The stored MPI communicator
  !>
  function get_comm() result(communicator)
    implicit none
    integer(i_def) :: communicator
    communicator = comm
  end function get_comm

  !> Returns whether the MPI communicator has been stored
  !> @return comm_state A flag indicating whether the MPI communicator is stored
  !>
  function is_comm_set() result(comm_state)
    implicit none
    logical(l_def) :: comm_state
    comm_state = comm_set
  end function is_comm_set

  !> Calculates the global sum of a collection of real local sums
  !>
  !> @param l_sum The sum of the reals on the local partition
  !> @param g_sum The calculated global sum
  !>
  subroutine global_sum_r_double(l_sum, g_sum)
    implicit none
    real(r_double), intent(in)  :: l_sum
    real(r_double), intent(out) :: g_sum

    integer(i_def) :: err

#ifdef NO_MPI
    ! Global sum and local sum are the same thing in a non-mpi build
    g_sum = l_sum
#else
    if(comm_set)then
      ! Generate global sum
      call mpi_allreduce( l_sum, g_sum, 1, get_mpi_datatype( real_type, r_double ), &
                          mpi_sum, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to real global_sum failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_sum failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_sum_r_double


  !> Calculates the global sum of a collection of real32 local sums
  !>
  !> @param l_sum The sum of the reals on the local partition
  !> @param g_sum The calculated global sum
  !>
  subroutine global_sum_r_single(l_sum, g_sum)
    implicit none
    real(r_single), intent(in)  :: l_sum
    real(r_single), intent(out) :: g_sum

    integer(i_def) :: err

#ifdef NO_MPI
    ! Global sum and local sum are the same thing in a non-mpi build
    g_sum = l_sum
#else
    if(comm_set)then
      ! Generate global sum
      call mpi_allreduce( l_sum, g_sum, 1, get_mpi_datatype( real_type, r_single), &
                          mpi_sum, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to real global_sum failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_sum failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_sum_r_single


  !> Calculates the global sum of a collection of integer local sums
  !>
  !> @param l_sum The sum of the integers on the local partition
  !> @param g_sum The calculated global sum
  !>
  subroutine global_sum_i_def(l_sum, g_sum)
    implicit none
    integer(i_def), intent(in)  :: l_sum
    integer(i_def), intent(out) :: g_sum

    integer(i_def):: err

#ifdef NO_MPI
    ! Global sum and local sum are the same thing in a non-mpi build
    g_sum = l_sum
#else
    if(comm_set)then
      ! Generate global sum
      call mpi_allreduce( l_sum, g_sum, 1, get_mpi_datatype( integer_type, i_def ), &
                          mpi_sum, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to integer global_sum failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_sum failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_sum_i_def


  !> Calculates the global minimum of a collection of local real minimums
  !>
  !> @param l_min The min on the local partition
  !> @param g_min The calculated global minimum
  !>
  subroutine global_min_r_double(l_min, g_min)
    implicit none
    real(r_double), intent(in)  :: l_min
    real(r_double), intent(out) :: g_min

    integer(i_def)  :: err

#ifdef NO_MPI
    ! Global minimum and local minimum are the same thing in a non-mpi build
    g_min = l_min
#else
    if(comm_set)then
      ! Generate global min
      call mpi_allreduce( l_min, g_min, 1, get_mpi_datatype( real_type, r_double ), &
                          mpi_min, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_min failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_min failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_min_r_double

  !> Calculates the global minimum of a collection of local real32 minimums
  !>
  !> @param l_min The min on the local partition
  !> @param g_min The calculated global minimum
  !>
  subroutine global_min_r_single(l_min, g_min)
    implicit none
    real(r_single), intent(in)  :: l_min
    real(r_single), intent(out) :: g_min

    integer(i_def)  :: err

#ifdef NO_MPI
    ! Global minimum and local minimum are the same thing in a non-mpi build
    g_min = l_min
#else
    if(comm_set)then
      ! Generate global min
      call mpi_allreduce( l_min, g_min, 1, get_mpi_datatype( real_type, r_single), &
                          mpi_min, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_min failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_min failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_min_r_single


  !> Calculates the global minimum of a collection of local integer minimums
  !>
  !> @param l_min The min on the local partition
  !> @param g_min The calculated global minimum
  !>
  subroutine global_min_i_def(l_min, g_min)
    implicit none
    integer(i_def), intent(in)  :: l_min
    integer(i_def), intent(out) :: g_min

    integer(i_def)  :: err

#ifdef NO_MPI
    ! Global minimum and local minimum are the same thing in a non-mpi build
    g_min = l_min
#else
    if(comm_set)then
      ! Generate global min
      call mpi_allreduce( l_min, g_min, 1, get_mpi_datatype( integer_type, i_def ), &
                          mpi_min, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_min failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_min failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_min_i_def


  !> Calculates the global maximum of a collection of local real maximums
  !>
  !> @param l_min The max on the local partition
  !> @param g_max The calculated global maximum
  !>
  subroutine global_max_r_double(l_max, g_max)
    implicit none
    real(r_double), intent(in)  :: l_max
    real(r_double), intent(out) :: g_max

    integer(i_def)  :: err

#ifdef NO_MPI
    ! Global maximum and local maximum are the same thing in a non-mpi build
    g_max = l_max
#else
    if(comm_set)then
      ! Generate global max
      call mpi_allreduce( l_max, g_max, 1, get_mpi_datatype( real_type, r_double ), &
                          mpi_max, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_max failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_max failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_max_r_double


  !> Calculates the global maximum of a collection of local real32 maximums
  !>
  !> @param l_min The max on the local partition
  !> @param g_max The calculated global maximum
  !>
  subroutine global_max_r_single(l_max, g_max)
    implicit none
    real(r_single), intent(in)  :: l_max
    real(r_single), intent(out) :: g_max

    integer(i_def)  :: err

#ifdef NO_MPI
    ! Global maximum and local maximum are the same thing in a non-mpi build
    g_max = l_max
#else
    if(comm_set)then
      ! Generate global max
      call mpi_allreduce( l_max, g_max, 1, get_mpi_datatype( real_type, r_single), &
                          mpi_max, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_max failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_max failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_max_r_single


  !> Calculates the global maximum of a collection of local integer maximums
  !>
  !> @param l_min The max on the local partition
  !> @param g_max The calculated global maximum
  !>
  subroutine global_max_i_def(l_max, g_max)
    implicit none
    integer(i_def), intent(in)  :: l_max
    integer(i_def), intent(out) :: g_max

    integer(i_def)  :: err

#ifdef NO_MPI
    ! Global maximum and local maximum are the same thing in a non-mpi build
    g_max = l_max
#else
    if(comm_set)then
      ! Generate global max
      call mpi_allreduce( l_max, g_max, 1, get_mpi_datatype( integer_type, i_def ), &
                          mpi_max, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to global_max failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to global_max failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif

  end subroutine global_max_i_def


  !> Gather integer data from all MPI tasks into a single array in all MPI tasks
  !> The data in send_buffer from the jth process is received by every
  !> process and placed in the jth block of the buffer recv_buffer.
  !>
  !> @param send_buffer The buffer of data to be sent to all MPI tasks
  !> @param recv_buffer The buffer into which the gathered data will be placed
  !> @param count The number of items in send_buffer
  subroutine all_gather(send_buffer, recv_buffer, count)
    implicit none
    integer(i_def), intent(in)  :: send_buffer(:)
    integer(i_def), intent(out) :: recv_buffer(:)
    integer(i_def), intent(in)  :: count

    integer(i_def) :: err

#ifdef NO_MPI
    ! Send and recv buffers in a gather are the same thing in a non-mpi build
    recv_buffer = send_buffer
#else
    if(comm_set)then
      call mpi_allgather(send_buffer, count, get_mpi_datatype( integer_type, i_def ), &
                         recv_buffer, count, get_mpi_datatype( integer_type, i_def ), &
                         comm, err)
      if (err /= mpi_success) &
        call log_event('Call to all_gather failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to all_gather failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end subroutine all_gather



  !> Broadcasts logical data from the root MPI task to all other MPI tasks
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_l_def(buffer, count, root)

    implicit none

    logical(l_def), intent(inout) :: buffer(:)
    integer(i_def), intent(in)    :: count
    integer(i_def), intent(in)    :: root

    integer(i_def) :: err

#ifdef NO_MPI
    ! In a non-mpi build there is nowhere to broadcast to - so do nothing
#else
    if(comm_set)then
      call mpi_bcast( buffer, count, MPI_LOGICAL, root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to logical broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end subroutine broadcast_l_def

  !> Broadcasts integer data from the root MPI task to all other MPI tasks
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_i_def(buffer, count, root)

    implicit none

    integer(i_def), intent(inout) :: buffer(:)
    integer(i_def), intent(in)    :: count
    integer(i_def), intent(in)    :: root

    integer(i_def) :: err

#ifdef NO_MPI
    ! In a non-mpi build there is nowhere to broadcast to - so do nothing
#else
    if(comm_set)then
      call mpi_bcast( buffer, count, get_mpi_datatype( integer_type, i_def ), root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to integer broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end subroutine broadcast_i_def

  !> Broadcasts double real data from the root MPI task to all other MPI tasks.
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_r_double(buffer, count, root)

    implicit none

    real(r_double), intent(inout) :: buffer(:)
    integer(i_def), intent(in)    :: count
    integer(i_def), intent(in)    :: root

    integer(i_def) :: err

#ifdef NO_MPI
    ! In a non-mpi build there is nowhere to broadcast to - so do nothing
#else
    if(comm_set)then
      call mpi_bcast( buffer, count, &
                      get_mpi_datatype( real_type, r_double ), &
                      root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to real broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end subroutine broadcast_r_double

  !> Broadcasts single real data from the root MPI task to all other MPI tasks.
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_r_single(buffer, count, root)

    implicit none

    real(r_single), intent(inout) :: buffer(:)
    integer(i_def), intent(in)    :: count
    integer(i_def), intent(in)    :: root

    integer(i_def) :: err

#ifdef NO_MPI
    ! In a non-mpi build there is nowhere to broadcast to - so do nothing
#else
    if(comm_set)then
      call mpi_bcast( buffer, count, &
                      get_mpi_datatype( real_type, r_single ), &
                      root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to real broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end subroutine broadcast_r_single

  !> Broadcasts character data from the root MPI task to all other MPI tasks
  !>
  !> @param buffer On the root MPI task, contains the data to broadcast,
  !>               on other tasks the data from root task will be writen to here
  !> @param count The number of items in buffer
  !> @param root The MPI task from which data will be broadcast
  subroutine broadcast_str(buffer, count, root)

    implicit none

    character(len=*), intent(inout) :: buffer(:)
    integer(i_def),   intent(in)    :: count
    integer(i_def),   intent(in)    :: root

    integer(i_def) :: err

#ifdef NO_MPI
    ! In a non-mpi build there is nowhere to broadcast to - so do nothing
#else
    if(comm_set)then
      call mpi_bcast( buffer, count, MPI_CHARACTER, root, comm, err )
      if (err /= mpi_success) &
        call log_event('Call to string broadcast failed with an MPI error.', &
                       LOG_LEVEL_ERROR )
    else
      call log_event( &
      'Call to broadcast failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end subroutine broadcast_str


  !> Returns the appropriate MPI datatype enumerator for all the Fortran
  !> kinds supported by the LFRic distributed memory code
  !>
  !> @param fortran_type An integer parameter indicating the Fortran data type
  !> @param fortran_kind A Fortran kind variable
  !> @return mpi_datatype The MPI datatype enumerator associated with the
  !>                      given Fortran type and kind
  function get_mpi_datatype( fortran_type, fortran_kind ) result(mpi_datatype)
    use, intrinsic :: iso_fortran_env, only : real128, real64, real32, &
                                              int64, int32, int16, int8
    implicit none
    integer(i_native), intent(in) :: fortran_type
    integer(i_native), intent(in) :: fortran_kind
    integer(i_native)             :: mpi_datatype

#ifdef NO_MPI
    ! In a non-mpi build the mpi datatype is meaningless - just return zero
    mpi_datatype = 0
#else
   ! Determine MPI datatype enumerator from a Fortran kind.
   ! (To support a new Fortran kind, just add a new case clause)
    select case (fortran_type)
    case (real_type)
      ! In the case where the data is real
      select case (fortran_kind)
      case (real32)
        mpi_datatype = MPI_REAL4
      case (real64)
        mpi_datatype = MPI_DOUBLE_PRECISION
      case (real128)
        call log_event( 'Attempt to use real128 Fortran kind used for MPI comms - &
           &NOT YET SUPPORTED', LOG_LEVEL_ERROR )
      case default
        call log_event( 'Unrecognised Fortran kind used for MPI comms', &
           LOG_LEVEL_ERROR )
      end select
    case (integer_type)
      ! In the case where the data is integer
      select case (fortran_kind)
      case (int8)
        mpi_datatype = MPI_INTEGER1
      case (int16)
        mpi_datatype = MPI_INTEGER2
      case (int32)
        mpi_datatype = MPI_INTEGER
      case (int64)
        mpi_datatype = MPI_INTEGER8
      case default
        call log_event( 'Unrecognised Fortran kind used for MPI comms', &
           LOG_LEVEL_ERROR )
      end select
    case (logical_type)
      mpi_datatype = MPI_LOGICAL
    end select
#endif

  end function get_mpi_datatype

  !> Returns the number of MPI ranks in the communicator
  !>
  !> @return c_size The number of MPI ranks in the communicator
  function get_comm_size() result(c_size)
    implicit none
    integer(i_def) :: c_size
#ifdef NO_MPI
    ! A non-mpi run is serial, therefore, number of ranks has to be one
    c_size = 1
#else
    if(comm_set)then
      c_size = comm_size
    else
      call log_event( &
      'Call to get_com_size failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end function get_comm_size

  !> Returns the number of the local MPI rank
  !>
  !> @return c_size The number of the local MPI rank
  function get_comm_rank() result(c_rank)
    implicit none
    integer(i_def) :: c_rank
#ifdef NO_MPI
    ! A non-mpi run is serial, therefore, local rank is always rank zero
    c_rank = 0
#else
    if(comm_set)then
      c_rank = comm_rank
    else
      call log_event( &
      'Call to get_com_rank failed. Must call store_comm first',&
      LOG_LEVEL_ERROR )
    end if
#endif
  end function get_comm_rank

end module mpi_mod
