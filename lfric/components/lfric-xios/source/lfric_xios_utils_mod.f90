!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Collection of small utility procedures for lfric-xios
!>
module lfric_xios_utils_mod

    use xios,           only : xios_date
    use constants_mod,  only : i_native

    implicit none
    private
    public :: parse_date_as_xios

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Interpret a string as an XIOS date object.
    !> Expected format is yyyy-mm-dd hh:mm:ss
    !>
    !> @param [in] date_str The string representation of the date
    !> @result     date_obj The xios_date object represented by the input string
    !>
    function parse_date_as_xios( date_str ) result( date_obj )
      implicit none
      character(*), intent(in)      :: date_str
      type(xios_date)               :: date_obj
      integer(i_native)             :: y, mo, d, h, mi, s, size

      size = len(date_str)

      ! Indexing from end to support arbitrarily long year
      read( date_str(1      :size-15), * ) y
      read( date_str(size-13:size-12), * ) mo
      read( date_str(size-10:size-9 ), * ) d
      read( date_str(size-7 :size-6 ), * ) h
      read( date_str(size-4 :size-3 ), * ) mi
      read( date_str(size-1 :size   ), * ) s

      date_obj = xios_date( y, mo, d, h, mi, s )

    end function parse_date_as_xios

  end module lfric_xios_utils_mod