!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!>  @brief Module containing an abstract file type
!!
!!  @details Provides an abstract file type, together with abstract
!!           procedure interfaces.
!-------------------------------------------------------------------------------
module file_mod

implicit none

private

!-------------------------------------------------------------------------------
!> @brief Abstract file type
!!
!! @details  Defines the interface for a family of IO strategies,
!!           which extend this abstract type.
!-------------------------------------------------------------------------------

!--------- File type ---------
type, abstract, public :: file_type
  private

contains
  procedure (new_open_interface ),      deferred :: file_open
  procedure (new_open_interface ),      deferred :: file_new
  procedure (close_interface),          deferred :: file_close

end type file_type

!-------------------------------------------------------------------------------
! Abstract interfaces
!-------------------------------------------------------------------------------
abstract interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Open an existing file, or create a new file.
  !!
  !! @param[in] self               The file strategy object.
  !! @param[in] file_name          Filename
  !-----------------------------------------------------------------------------

  subroutine new_open_interface(self, file_name)
    import :: file_type

    !Arguments
    class(file_type),       intent(inout) :: self
    character(len=*),       intent(in)    :: file_name

  end subroutine new_open_interface

  !-----------------------------------------------------------------------------
  !> @brief  Interface: Close a file
  !!
  !! @param[in] self               The file strategy object.
  !-----------------------------------------------------------------------------

  subroutine close_interface(self)
    import :: file_type

    !Arguments
    class(file_type), intent(inout) :: self

  end subroutine close_interface

end interface

contains

end module file_mod