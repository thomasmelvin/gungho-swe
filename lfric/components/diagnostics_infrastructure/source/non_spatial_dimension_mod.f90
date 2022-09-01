!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the non_spatial dimension types
!>
!> @details This allows for the creation of an object that represents the
!> non-spatial dimension meta data for a field object.
!> A meta_data_type can be supplied with, as an optional argument, a
!> non_spatial_dimension_type to define its non-spatial dimensional meta data.


module non_spatial_dimension_mod

  use constants_mod,                              only: str_def,   &
                                                        str_short, &
                                                        r_def,     &
                                                        str_long,  &
                                                        i_native
  use log_mod,                                    only: log_event, &
                                                        LOG_LEVEL_ERROR

  implicit none

  private
  public :: NUMERICAL, CATEGORICAL

  !> Enumerator listing the possible non-spatial dimension categories
  !>
  !> Numerical - a numerical axis of real numbers
  !> Categorical - a list of categories with strings for labels
  enum, bind(c)
    enumerator :: NUMERICAL,      &
                  CATEGORICAL
  end enum

  !> Defines the configurable non_spatial dimension
  type, public :: non_spatial_dimension_type

    !> @brief Name of the non-spatial dimension
    character(str_def)                              :: dimension_name

    !> @brief Optional. Category of the dimension: either numerical or categorical
    !>
    !> Specified using the enums NUMERICAL / CATEGORICAL within
    !> \link non_spatial_dimension_mod.f90 non_spatial_dimension_mod.f90
    !> \endlink
    integer(i_native)                               :: dimension_category

    !> Description to be displayed in the Rose GUI help text
    character(str_long)                             :: help_text

    !> @brief Optional. Definition of the axis if it is labelled using strings
    !>
    !> Array of strings to use as the labels
    !> e.g. `['fine','medium','coarse']`
    character(str_short), dimension(:), allocatable :: label_definition

    !> @brief Definition of the axis if it is described by real numbers
    !>
    !> Array of real numbers. These define the **boundaries** of the data
    !> points.
    real     (r_def),     dimension(:), allocatable :: axis_definition

    !> @brief Optional. Units of the non-spatial dimension
    !>
    !> If the dimension has a suitable unit this should specify it so that
    !> the resultant netCDF file has the unit included.
    character(str_def)                              :: non_spatial_units

  contains

    procedure, public :: get_dimension_name
    procedure, public :: get_axis_definition
    procedure, public :: get_label_definition

  end type non_spatial_dimension_type

  interface non_spatial_dimension_type
    module procedure non_spatial_dimension_type_constructor
  end interface non_spatial_dimension_type

  contains

  !> Construct a <code>non_spatial_dimension_type</code> object.
  !> @brief The constructor for a non_spatial_dimension_type object
  !> @param [in] dimension_name The name of the non-spatial dimension
  !> @return self the non_spatial_dimension object
  function non_spatial_dimension_type_constructor(dimension_name,     &
                                                  dimension_category, &
                                                  help_text,          &
                                                  axis_definition,    &
                                                  label_definition,   &
                                                  non_spatial_units) result(self)

    implicit none

    character(*),                                 intent(in) :: dimension_name
    integer  (i_native),                          intent(in) :: dimension_category
    character(*),                                 intent(in) :: help_text
    real     (r_def),     dimension(:), optional, intent(in) :: axis_definition
    character(str_short), dimension(:), optional, intent(in) :: label_definition
    character(*),                       optional, intent(in) :: non_spatial_units

    type(non_spatial_dimension_type) :: self

    self%dimension_name = dimension_name
    self%dimension_category = dimension_category
    self%help_text = help_text

    if (present(axis_definition) .and. present(label_definition)) then
      call log_event( "A non_spatial_dimension_type cannot have an "// &
                      "axis_definition and a label_definition", LOG_LEVEL_ERROR )

    else if (dimension_category == NUMERICAL .and. present(label_definition)) then
      call log_event( "Numerical non-spatial dimensions cannot have a "// &
                      "label_definition", LOG_LEVEL_ERROR )

    else if (dimension_category == CATEGORICAL .and. present(axis_definition)) then
      call log_event( "Categorical non-spatial dimensions cannot have "// &
                      "an axis_definition", LOG_LEVEL_ERROR )
    end if

    if(present(axis_definition)) then
      allocate(self%axis_definition, source=axis_definition)
    end if

    if (present(label_definition)) then
      allocate(self%label_definition, source=label_definition)
    end if

    if (present(non_spatial_units)) then
      self%non_spatial_units = non_spatial_units
    end if

  end function non_spatial_dimension_type_constructor


  !> Getter for dimension_name
  !> @param[in]  self non_spatial_dimension_type
  !> @return dimension name
  function get_dimension_name(self) result(dimension_name)

    implicit none

    class(non_spatial_dimension_type), intent(in) :: self

    character(str_def)                            :: dimension_name

    dimension_name = self%dimension_name

  end function get_dimension_name


  !> Getter for axis_definition
  !> @param[in]  self  non_spatial_dimension_type
  !> @return axis_definition
  function get_axis_definition(self) result(axis_definition)

    implicit none

    class(non_spatial_dimension_type), intent(in) :: self

    real(r_def), allocatable                      :: axis_definition(:)

    allocate(axis_definition, source=self%axis_definition)

  end function get_axis_definition


  !> Getter for a label_definition
  !> @param[in]  self  non_spatial_dimension_type
  !> @return label_definition
  function get_label_definition(self) result(label_definition)

    implicit none

    class(non_spatial_dimension_type),   intent(in) :: self

    character(str_short), dimension(:), allocatable :: label_definition

    label_definition = self%label_definition

  end function get_label_definition


end module non_spatial_dimension_mod