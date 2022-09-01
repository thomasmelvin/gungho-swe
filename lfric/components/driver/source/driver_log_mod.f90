module driver_log_mod

use constants_mod,        only: i_native
use convert_to_upper_mod, only: convert_to_upper
use log_mod,              only: log_event,          &
                                log_set_level,      &
                                log_scratch_space,  &
                                initialise_logging, &
                                finalise_logging,   &
                                LOG_LEVEL_ALWAYS,   &
                                LOG_LEVEL_ERROR,    &
                                LOG_LEVEL_WARNING,  &
                                LOG_LEVEL_INFO,     &
                                LOG_LEVEL_DEBUG,    &
                                LOG_LEVEL_TRACE
use logging_config_mod,   only: run_log_level,          &
                                key_from_run_log_level, &
                                RUN_LOG_LEVEL_ERROR,    &
                                RUN_LOG_LEVEL_INFO,     &
                                RUN_LOG_LEVEL_DEBUG,    &
                                RUN_LOG_LEVEL_TRACE,    &
                                RUN_LOG_LEVEL_WARNING

implicit none

public :: init_logger, final_logger
private

contains

subroutine init_logger(total_ranks, local_rank, program_name)

  implicit none

  character(len=*),  intent(in) :: program_name
  integer(i_native), intent(in) :: total_ranks, local_rank

  integer(i_native) :: log_level

  call initialise_logging(local_rank, total_ranks, program_name)

  select case (run_log_level)
  case( RUN_LOG_LEVEL_ERROR )
    log_level = LOG_LEVEL_ERROR
  case( RUN_LOG_LEVEL_WARNING )
    log_level = LOG_LEVEL_WARNING
  case( RUN_LOG_LEVEL_INFO )
    log_level = LOG_LEVEL_INFO
  case( RUN_LOG_LEVEL_DEBUG )
    log_level = LOG_LEVEL_DEBUG
  case( RUN_LOG_LEVEL_TRACE )
    log_level = LOG_LEVEL_TRACE
  case default
    call log_event( "Invalid option for run_log_level", LOG_LEVEL_ERROR )
  end select

  call log_set_level( log_level )

  write(log_scratch_space,'(A)')                              &
      'Runtime message logging severity set to log level: '// &
      convert_to_upper(key_from_run_log_level(run_log_level))
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

end subroutine init_logger

subroutine final_logger(program_name)

  implicit none

  character(len=*), intent(in) :: program_name

  ! Final logging before infrastructure is destroyed
  call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

  ! Finalise the logging system
  call finalise_logging()

end subroutine final_logger

end module driver_log_mod
