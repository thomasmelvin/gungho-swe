!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> Dynamo program support functions.
!>
!> Originally these were "block" constructs within the program but neither
!> GNU or Intel Fortran where properly able to cope with that.
!>
module gungho_mod

  use log_mod, only : log_event,         &
                      log_scratch_space, &
                      LOG_LEVEL_ALWAYS,  &
                      LOG_LEVEL_ERROR,   &
                      LOG_LEVEL_TRACE,   &
                      LOG_LEVEL_DEBUG

  use configuration_mod, only : read_configuration,   &
                                ensure_configuration

  implicit none

  private
  public :: load_configuration, program_name

  character(*), parameter :: program_name = 'gungho'

contains

  !> Loads run-time configuration and ensures everything is ship-shape.
  !>
  !> @param file_unit I/O unit for file holding namelists.
  !>
  subroutine load_configuration( filename )

    use check_configuration_mod,  only : check_configuration

    implicit none

    character(*), intent(in) :: filename

    character(*), parameter :: &
                  required_configuration(12) = ['finite_element             ', &
                                                'formulation                ', &
                                                'base_mesh                  ', &
                                                'initial_wind               ', &
                                                'planet                     ', &
                                                'solver                     ', &
                                                'mixed_solver               ', &
                                                'helmholtz_solver           ', &
                                                'timestepping               ', &
                                                'extrusion                  ', &
                                                'transport                  ', &
                                                'orography                  ']

    logical              :: okay
    logical, allocatable :: success_map(:)
    integer              :: i

    allocate( success_map(size(required_configuration)) )

    call log_event( 'Loading '//program_name//' configuration ...', &
                    LOG_LEVEL_ALWAYS )

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

    call check_configuration()

  end subroutine load_configuration

end module gungho_mod
