!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different time steps
!> @details This enumarator contains all the timesteps available for an LFRic
!> field. It is used by vertical_dimensions_mod.f90 and field meta data
!> definition files for specifying the timestep that the field is output on
!> This file is solely for unit testing purposes
module time_step_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY. DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)
    enumerator :: TEST_TIMESTEP_1 = 3001, &
                  TEST_TIMESTEP_2 = 3002
  end enum
end module time_step_enum_mod