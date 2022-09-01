!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page gung_ho GungHo Program
!> This is a code that uses the LFRic infrastructure to build a model that
!> just includes the GungHo dynamical core.

!> @brief Main program used to illustrate gungho functionality.

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program gungho

  use gungho_driver_mod, only : initialise, run, finalise

  implicit none

  call initialise()

  call run()

  call finalise()

end program gungho
