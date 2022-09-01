!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different model levels
!> @details This enumerator contains all the possible model levels
!> for an LFRic field. It is used by vertical_dimensions_mod.f90 and
!> field meta data definition files for specifying model levels
!> This file is solely for unit testing purposes

module levels_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY.
  !> DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)
    enumerator :: BOTTOM_ATMOSPHERIC_LEVEL = 1001, &
                  TOP_ATMOSPHERIC_LEVEL = 1002, &
                  BOTTOM_SOIL_LEVEL = 1003, &
                  TOP_SOIL_LEVEL = 1004
  end enum
end module levels_enum_mod