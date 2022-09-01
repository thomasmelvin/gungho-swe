!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Shallow water program support functions.
!> @details Functions to support shallow water miniapp. Originally these
!!          were "block" constructs within the program but neither
!!          GNU or Intel Fortran where properly able to cope with that.

module shallow_water_mod

  use log_mod, only : log_event,         &
                      log_scratch_space, &
                      LOG_LEVEL_ERROR,   &
                      LOG_LEVEL_INFO,    &
                      LOG_LEVEL_TRACE,   &
                      LOG_LEVEL_DEBUG

  use configuration_mod, only : read_configuration, &
                                ensure_configuration

  implicit none

  private
  public :: load_configuration, program_name

  character(*), parameter :: program_name = 'shallow_water'

contains

  !> @brief Loads run-time configuration and ensures everything is ship-shape.
  !> @param[in] filename I/O unit for file holding namelists.
  subroutine load_configuration( filename )

    implicit none

    character(*), intent(in) :: filename

    character(*), parameter :: &
                            required_configuration(8) = ['base_mesh             ', &
                                                         'files                 ', &
                                                         'formulation           ', &
                                                         'planet                ', &
                                                         'shallow_water_settings', &
                                                         'time                  ', &
                                                         'timestepping          ', &
                                                         'transport             ']

    logical              :: okay
    logical, allocatable :: success_map(:)
    integer              :: i

    allocate( success_map(size(required_configuration)) )

    call log_event( 'Loading '//program_name//' configuration ...', &
                    LOG_LEVEL_INFO )

    call read_configuration( filename )

    okay = ensure_configuration( required_configuration, success_map )
    if (.not. okay) then
      write( log_scratch_space, '(A)' ) &
                             'The following required namelists were not loaded:'
      do i = 1,size(required_configuration)
        if (.not. success_map(i)) &
          log_scratch_space = trim(log_scratch_space) // ' ' &
                              // required_configuration(i)
      end do
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    deallocate( success_map )

  end subroutine load_configuration

end module shallow_water_mod
