!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different interpolation methods
!> @details This enumerator contains all the possible interpolation
!> methods for an LFRic field. It is used by diagnostics_mod.f90 for
!> specifying the recommended interpolation type.
!> This file is solely for unit testing purposes

module interpolation_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY.
  !> DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)

    enumerator :: TEST_BILINEAR = 4002
    enumerator :: TEST_TRILINEAR = 4003

  end enum

end module interpolation_enum_mod