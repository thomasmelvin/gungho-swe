!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_entity_map_mod
! A module containing a collection of helper routines that provide canned
! entity maps which (for each function space) list the entities on which each
! dof lives on a reference cube
! Separate functions are supplied for lowest order and "higher order".
! Essentially the list of entities is:
! 1: The single cell volume
! 1000+X, where X is 1 to 6 : The 6 faces
! 2000+X, where X is 1 to 12: The 12 edges
! 3000+X, where X is 1 to 8 : The 8 verticies

  use constants_mod, only: i_def

  implicit none

  private

  public :: get_w0_o0_entity_map,            &
            get_w1_o0_entity_map,            &
            get_w2_o0_entity_map,            &
            get_w3_o0_entity_map,            &
            get_w2v_o0_entity_map,           &
            get_w2h_o0_entity_map,           &
            get_w2broken_o0_entity_map,      &
            get_w2trace_o0_entity_map,       &
            get_wtheta_o0_entity_map,        &
            get_wchi_o0_entity_map,          &
            get_w0_o1_entity_map,            &
            get_w1_o1_entity_map,            &
            get_w2_o1_entity_map,            &
            get_w3_o1_entity_map,            &
            get_w2v_o1_entity_map,           &
            get_w2h_o1_entity_map,           &
            get_w2broken_o1_entity_map,      &
            get_w2trace_o1_entity_map,       &
            get_wtheta_o1_entity_map,        &
            get_wchi_o1_entity_map
contains

!------------------------------------------------------

! Subroutines for returning lowest-order entity maps

subroutine get_w0_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 8 ) )
  entity_map = (/ 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008 /)

end subroutine get_w0_o0_entity_map

!------------------------------------------------------

subroutine get_w1_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(12 ) )
  entity_map = (/ 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, &
                  2009, 2010, 2011, 2012 /)

 end subroutine get_w1_o0_entity_map

!------------------------------------------------------

subroutine get_w2_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 6 ))
  entity_map = (/ 1001, 1002, 1003, 1004, 1005, 1006 /)

end subroutine get_w2_o0_entity_map

!------------------------------------------------------

subroutine get_w3_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 1 ))
  entity_map = (/    1 /)

end subroutine get_w3_o0_entity_map

!------------------------------------------------------

subroutine get_w2trace_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 6 ))
  entity_map = (/ 1001, 1002, 1003, 1004, 1005, 1006 /)

end subroutine get_w2trace_o0_entity_map

!------------------------------------------------------

subroutine get_w2broken_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 6 ))
  entity_map = (/ 1001, 1002, 1003, 1004, 1005, 1006 /)

end subroutine get_w2broken_o0_entity_map

!------------------------------------------------------

subroutine get_wchi_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 1 ))
  entity_map = (/    1 /)

end subroutine get_wchi_o0_entity_map

!------------------------------------------------------

subroutine get_w2v_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 2 ))
  entity_map = (/ 1005, 1006 /)

end subroutine get_w2v_o0_entity_map

!------------------------------------------------------

subroutine get_w2h_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 4 ))
  entity_map = (/ 1001, 1002, 1003, 1004 /)

end subroutine get_w2h_o0_entity_map

!------------------------------------------------------

subroutine get_wtheta_o0_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 2 ))
  entity_map = (/ 1005, 1006 /)

end subroutine get_wtheta_o0_entity_map

!------------------------------------------------------

! Subroutines for returning order = 1 entity maps (higher order)

subroutine get_w0_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(27 ) )
  entity_map = (/    1, 1001, 1002, 1003, 1004, 1005, 1006, 2001, &
                  2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, &
                  2010, 2011, 2012, 3001, 3002, 3003, 3004, 3005, &
                  3006, 3007, 3008 /)

end subroutine get_w0_o1_entity_map

!------------------------------------------------------

subroutine get_w1_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(54 ) )
  entity_map = (/    1,    1,    1,    1,    1,    1, 1001, 1001, &
                  1001, 1001, 1002, 1002, 1002, 1002, 1003, 1003, &
                  1003, 1003, 1004, 1004, 1004, 1004, 1005, 1005, &
                  1005, 1005, 1006, 1006, 1006, 1006, 2001, 2001, &
                  2002, 2002, 2003, 2003, 2004, 2004, 2005, 2005, &
                  2006, 2006, 2007, 2007, 2008, 2008, 2009, 2009, &
                  2010, 2010, 2011, 2011, 2012, 2012 /)

end subroutine get_w1_o1_entity_map

!------------------------------------------------------

subroutine get_w2_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(36 ))
  entity_map = (/    1,    1,    1,    1,    1,    1,    1,    1, &
                     1,    1,    1,    1, 1001, 1001, 1001, 1001, &
                  1002, 1002, 1002, 1002, 1003, 1003, 1003, 1003, &
                  1004, 1004, 1004, 1004, 1005, 1005, 1005, 1005, &
                  1006, 1006, 1006, 1006 /)


end subroutine get_w2_o1_entity_map

!------------------------------------------------------

subroutine get_w3_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 8 ))
  entity_map = (/    1,    1,    1,    1,    1,    1,    1,    1 /)

end subroutine get_w3_o1_entity_map

!------------------------------------------------------

subroutine get_w2trace_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(24 ))
  entity_map = (/ 1001, 1001, 1001, 1001, 1002, 1002, 1002, 1002, &
                  1003, 1003, 1003, 1003, 1004, 1004, 1004, 1004, &
                  1005, 1005, 1005, 1005, 1006, 1006, 1006, 1006 /)

end subroutine get_w2trace_o1_entity_map

!------------------------------------------------------

subroutine get_w2broken_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(36 ))
  entity_map = (/    1,    1,    1,    1,    1,    1,    1,    1, &
                     1,    1,    1,    1, 1001, 1001, 1001, 1001, &
                  1002, 1002, 1002, 1002, 1003, 1003, 1003, 1003, &
                  1004, 1004, 1004, 1004, 1005, 1005, 1005, 1005, &
                  1006, 1006, 1006, 1006 /)


end subroutine get_w2broken_o1_entity_map

!------------------------------------------------------

subroutine get_wchi_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map( 8 ))
  entity_map = (/    1,    1,    1,    1,    1,    1,    1,    1 /)

end subroutine get_wchi_o1_entity_map

!------------------------------------------------------

subroutine get_w2v_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(12 ))
  entity_map = (/    1,    1,    1,    1, 1005, 1005, 1005, 1005, &
                  1006, 1006, 1006, 1006 /)

end subroutine get_w2v_o1_entity_map

!------------------------------------------------------

subroutine get_w2h_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(24 ))
  entity_map = (/    1,    1,    1,    1,    1,    1,    1,    1, &
                  1001, 1001, 1001, 1001, 1002, 1002, 1002, 1002, &
                  1003, 1003, 1003, 1003, 1004, 1004, 1004, 1004 /)

end subroutine get_w2h_o1_entity_map

!------------------------------------------------------

subroutine get_wtheta_o1_entity_map(entity_map)

  implicit none

  integer(i_def), intent(out),allocatable :: entity_map(:)

  allocate(entity_map(12 ))
  entity_map = (/    1,    1,    1,    1, 1005, 1005, 1005, 1005, &
                  1006, 1006, 1006, 1006 /)

end subroutine get_wtheta_o1_entity_map

end module get_unit_test_entity_map_mod
