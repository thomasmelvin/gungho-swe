!-------------------------------------------------------------------------------
!(c) Crown copyright 2020 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!>@brief Updates time-varying fields from linked list of time axis objects
module variable_fields_mod

  use clock_mod,                     only : clock_type
  use constants_mod,                 only : r_def, r_second, str_def
  use field_collection_mod,          only : field_collection_type
  use lfric_xios_time_axis_mod,      only : time_axis_type
  use linked_list_mod,               only : linked_list_type, &
                                            linked_list_item_type

  implicit none

  private
  public :: init_variable_fields, &
            update_variable_fields

contains

  subroutine init_variable_fields(time_axis_list, clock, fields)

    implicit none

    type(linked_list_type),      intent(in)    :: time_axis_list
    class(clock_type),           intent(in)    :: clock
    type(field_collection_type), intent(inout) :: fields

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: loop => null()
    type(time_axis_type),        pointer :: time_axis => null()

    ! start at the head of the time_axis linked list
    loop => time_axis_list%get_head()
    do
      ! If list is empty or we're at the end of list, exit
      if ( .not. associated(loop) ) then
        exit
      end if

      ! Select the time_axis_type to get at the information in the list payload
      select type( list_item => loop%payload )
        type is (time_axis_type)
          time_axis => list_item

          ! Align time window and populate model data
          call time_axis%align()
          call time_axis%update_fields()
          ! Only populate the fields on a new run
          if ( clock%is_initialisation() ) then
            call time_axis%populate_model_fields(fields)
          end if

      end select

      loop => loop%next

    end do

    nullify(loop)
    nullify(time_axis)

  end subroutine init_variable_fields

  subroutine update_variable_fields(time_axis_list, clock, fields)

    implicit none

    type(linked_list_type),      intent(in)    :: time_axis_list
    class(clock_type),           intent(in)    :: clock
    type(field_collection_type), intent(inout) :: fields

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: loop => null()
    type(time_axis_type),        pointer :: time_axis => null()

    ! Start at the head of the time_axis linked list
    loop => time_axis_list%get_head()
    do
      ! If list is empty or we're at the end of list, exit
      if ( .not. associated(loop) ) then
        exit
      end if

      ! Select the time_axis_type to get at the information in the list payload
      select type( list_item => loop%payload )
        type is (time_axis_type)
          time_axis => list_item

          ! Step time axis forward
          call time_axis%step()

          ! Populate the model fields from the time axis data
          if ( time_axis%populate_fields() ) then
            call time_axis%populate_model_fields(fields)
          end if

      end select

      loop => loop%next

    end do

    nullify(loop)
    nullify(time_axis)

  end subroutine update_variable_fields

end module variable_fields_mod
