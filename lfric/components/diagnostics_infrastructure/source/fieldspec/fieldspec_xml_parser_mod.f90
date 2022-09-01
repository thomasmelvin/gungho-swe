!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Procedures to parse an iodef.xml file and add a fieldspec object to
!> the collection
!>
!> @details Uses the FoX library to parse an iodef.xml file and pass the values
!> to a fieldspec_factory object. This is then used to create a
!> fieldspec object , which is added to the fieldspec_collection
!>
!> The parser expects to find an iodef.xml file with a field_definition
!> section with the following structure (shown here with details
!> ignored by the parser removed for clarity):
!>
!>     <field_definition>
!>       <field_group id="example_field_group_1" enabled=".TRUE.">
!>         <field id="example_field_1">
!>           <variable name="variable_1_name">variable_1_value</variable>
!>         </field>
!>       </field_group>
!>     </field_definition>

module fieldspec_xml_parser_mod

  use constants_mod,                only: i_def, str_def, l_def, r_def
  use axisspec_collection_mod,      only: axisspec_collection_type
  use field_type_enum_mod,          only: field_type_from_name
  use fieldspec_mod,                only: fieldspec_type
  use fieldspec_collection_mod,     only: fieldspec_collection_type
  use fieldspec_factory_mod,        only: fieldspec_factory_type
  use fs_continuity_mod,            only: functionspace_from_name
  use fox_sax,                      only: XML_T, dictionary_t, getValue, hasKey, &
                                          parse, open_xml_file, close_xml_t
  use io_driver_enum_mod,           only: io_driver_from_name
  use non_spatial_dimension_mod,    only: NUMERICAL, CATEGORICAL
  use log_mod,                      only: log_event, log_scratch_space, &
                                          LOG_LEVEL_ERROR

  implicit none

  private
  public :: populate_fieldspec_collection

  type(fieldspec_collection_type),     pointer :: fieldspec_collection => null()
  type(axisspec_collection_type),      pointer :: axisspec_collection => null()

  !> @brief Factory used to build up fieldspec objects while reading XML
  !>        before adding them to the collection
  type(fieldspec_factory_type) :: fieldspec_factory
  logical :: fieldspec_factory_initialised = .false.


  ! Module-level variables used by XML handler methods called from populate()

  ! Tells the parser when to skip lines in the XML file
  logical :: ignore_element = .false.

  ! Flags used to check which element the XML parser is currently within
  logical :: in_axis_def = .false.
  logical :: in_field_def = .false.
  logical :: in_field_group = .false.
  logical :: in_field = .false.
  logical :: in_field_variable = .false.

  ! Stores the name of the current XML "field" element's
  ! "variable" while it is processed by the XML event handlers
  character(str_def) :: xml_field_variable
  character(str_def) :: field_group_id
  character(str_def) :: field_id
  character(str_def), allocatable :: axis_ids(:)

contains

  !===========================================================================
  !> @brief Populates the fieldspec collection from an XIOS iodef.xml file
  !> @param[in] iodef_filepath The path to the iodef.xml file
  !>
  subroutine populate_fieldspec_collection( iodef_filepath )

    implicit none

    character(len = *),               intent(in)       :: iodef_filepath
    type(XML_T)                                        :: parser
    integer                                            :: iostatus

    fieldspec_collection => fieldspec_collection_type()
    axisspec_collection => axisspec_collection_type()

    call open_xml_file(parser, iodef_filepath, iostatus)
    if (iostatus /= 0) then
      write(log_scratch_space, '(A)') 'Error opening file'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
      stop
    end if

    ! Call FoX's parse subroutine to begin reading through the XML file
    ! Tells the FoX API the mapping of interfaces to the subroutines containing
    ! instructions for event handling
    ! (e.g event when reaching an opening/closing tag of an XML element)
    call parse( parser, &
                startElement_handler = startElement_handler, &
                characters_handler = text_handler, &
                endElement_handler = endElement_handler )

    ! end parsing of XML file
    call close_xml_t(parser)

    return
  end subroutine populate_fieldspec_collection

  !===========================================================================
  !> @brief "Switches" on/off the logicals used to keep track of which
  !> XML elements the parser is currently within
  !>
  !> Intended for use by XML start/endElement_handler rather than being called
  !> manually
  !>
  !> @param[in] element_name The name of the flag whose value is to be set
  !> @param[in] new_value The value to set the flag to (.true. or .false.)
  subroutine switch_xml_element_flag( element_name, new_value )

    implicit none

    character(len=*), intent(in)    :: element_name
    logical, intent(in)             :: new_value

    select case (element_name)

      case ("axis_definition")
        in_axis_def = new_value

      case ("field_definition")
        in_field_def = new_value

      case ("field_group")
        in_field_group = new_value

      case ("field")
        in_field = new_value

      case("variable")
        if (in_field) then
          in_field_variable = new_value
        end if

    end select

    return
  end subroutine switch_xml_element_flag

  !===========================================================================
  !> @brief Handles the event where the FoX SAX parser reaches an XML element's
  !>        opening tag
  !>
  !> This subroutine is only called in the populate subroutine using the FoX API
  !> and is not to be called manually
  !>
  !> @param[in] namespaceURI
  !> @param[in] localname
  !> @param[in] name
  !> @param[in] attributes
  subroutine startElement_handler( namespaceURI, localname, name, attributes )

    implicit none

    character(len = *), intent(in)    :: namespaceURI
    character(len = *), intent(in)    :: localname
    character(len = *), intent(in)    :: name
    type(dictionary_t), intent(in)    :: attributes

    ! Strings parsed from xml
    character(str_def)                :: axis_size_str
    character(str_def)                :: numeric_axis_def_str
    character(str_def)                :: label_axis_def_str
    character(str_def)                :: grid_id

    ! Attributes for axisspec objects
    character(str_def)                :: axis_id
    integer(i_def)                    :: axis_size
    real(r_def),        allocatable   :: numeric_axis_def(:)
    character(str_def), allocatable   :: label_axis_def(:)

    ! Ignore element if parser is still in a field group that is not enabled
    if (ignore_element) then
      return
    end if

    ! Set the flag corresponding to this element to true to signify parser is
    ! inside it
    call switch_xml_element_flag(name, .true.)

    ! If in axis node, then create axisspec object for it
    if (in_axis_def .and. name == 'axis') then

      axis_id = getValue(attributes, 'id')
      axis_size_str = getValue(attributes, 'n_glo')
      numeric_axis_def_str = getValue(attributes, 'value')
      label_axis_def_str = getValue(attributes, 'label')

      read(axis_size_str, '(i3)') axis_size

      ! Create an axis definiton for numerical axes
      if (numeric_axis_def_str /= "") then
        numeric_axis_def = numeric_axis_from_string(numeric_axis_def_str, axis_size)
        call axisspec_collection%generate_and_add_axisspec( &
              axis_id, &
              NUMERICAL, &
              numeric_axis_def=numeric_axis_def)

      ! Create a label definition for categorical axes
      else if (label_axis_def_str /= "") then
        label_axis_def = label_axis_from_string(label_axis_def_str, axis_size)
        call axisspec_collection%generate_and_add_axisspec( &
              axis_id, &
              CATEGORICAL, &
              label_axis_def=label_axis_def)
      end if

    else if (in_field_def) then

      ! Set parser to ignore current group of fields if they are not enabled
      if (name == "field_group") then
        field_group_id = getValue(attributes, "id")
        if (getValue(attributes, "enabled") == ".FALSE.") then
          ignore_element = .true.
          return
        end if
      end if

      ! Store field's ID to be put in fieldspec object in endElement_handler
      if (in_field_group) then
        if (name == "field" .and.  hasKey(attributes, "id")) then

          ! Clear the fieldspec factory of any previous data
          call fieldspec_factory%initialise()
          fieldspec_factory_initialised = .true.
          field_id = getValue(attributes, "id")
          call fieldspec_factory%set_unique_id( field_id )
          call fieldspec_factory%set_field_group_id( field_group_id )

          ! Get axis ids from grid id
          grid_id = getValue(attributes, "grid_ref")
          axis_ids = axis_ids_from_grid_id(grid_id)

        end if

        ! Store name of variable while it is processed by text_handler
        if (name == "variable" .and. in_field) then
          xml_field_variable = getValue(attributes, "name")
        end if

      end if

    end if

    return
  end subroutine startElement_handler

  !===========================================================================
  !> @brief Handles the event where the FoX SAX parser
  !>        starts reading an XML element's text
  !>
  !> The text corresponds to the value of variable whose name is stored in
  !> xml_field_variable. It is stored in the fieldspec_factory to be to be put
  !> in a fieldspec object in endElement_handler
  !>
  !> This subroutine is only called in the populate subroutine using the FoX API
  !> and is not to be called manually
  !>
  !> @param[in] text The text to be handled
  subroutine text_handler( text )
    use fox_common, only : rts
    implicit none

    character(len = *), intent(in) :: text
    integer(i_def)                 :: int_from_char
    character(str_def)             :: trimmed_char
    logical(l_def)                 :: checksum

    if (in_field .and. in_field_variable) then

      select case(xml_field_variable)
        case ("mesh_id")
          ! rts function to convert XML text to integer
          call rts(text, int_from_char)
          call fieldspec_factory%set_mesh_id(int_from_char)

        case ("function_space")
          ! rts function used to trim whitespace on BOTH sides
          call rts(text, trimmed_char)
          call fieldspec_factory%set_function_space( functionspace_from_name(trimmed_char) )

        case ("order", "element_order")
          call rts(text, int_from_char)
          call fieldspec_factory%set_order(int_from_char)

        case ("field_kind")
          call rts(text, int_from_char)
          call fieldspec_factory%set_kind(int_from_char)

        case ("field_type")
          call rts(text, trimmed_char)
          call fieldspec_factory%set_type( field_type_from_name(trimmed_char) )

        case ("io_driver")
          call rts(text, trimmed_char)
          call fieldspec_factory%set_io_driver( io_driver_from_name(trimmed_char))

        case ("checksum")
          call rts(text, checksum)
          call fieldspec_factory%set_checksum(checksum)
      end select

    end if

    return
  end subroutine text_handler

  !===========================================================================
  !> @brief Handles the event where the FoX SAX parser reaches an XML element's
  !>        closing tag
  !>
  !> This subroutine is only called in the populate subroutine using the FoX API
  !> and is not to be called manually
  !>
  !> @param[in] namespaceURI
  !> @param[in] localname
  !> @param[in] name
  subroutine endElement_handler( namespaceURI, localname, name )
    implicit none

    character(len = *),              intent(in)  :: namespaceURI
    character(len = *),              intent(in)  :: localname
    character(len = *),              intent(in)  :: name

    type(fieldspec_type), pointer :: fieldspec
    integer(i_def) :: i

    ! Stop parser ignoring lines now that it has finished reading the disabled
    ! field_group
    if (ignore_element .and. name == "field_group") then
      ignore_element = .false.
    end if

    ! Create fieldspec object when the parser is done collating its
    ! properties from the XML field
    if (in_field_def .and. in_field_group .and. name == "field" &
            .and. fieldspec_factory_initialised) then

      call fieldspec_collection%add_fieldspec( fieldspec_factory%finalise() )
      fieldspec => fieldspec_collection%get_fieldspec(field_id)
      ! Set each axis in the fieldspec
      do i = 1, size(axis_ids)
        if (axis_ids(i)(7:15) == 'vert_axis') then
          call fieldspec%set_vertical_axis( &
                  axisspec_collection%get_axisspec(axis_ids(i)))
        else
          call fieldspec%add_non_spatial_dimension( &
                  axisspec_collection%get_axisspec(axis_ids(i)))
        end if
      end do
      fieldspec_factory_initialised = .false.
      deallocate(axis_ids)

    end if

    ! Switch off flag for current element now the parser has finished reading it
    call switch_xml_element_flag(name, .false.)

    return
  end subroutine endElement_handler


  !===========================================================================
  !> Get an array of the level values from the string in the iodef.xml file
  !> @param[in] axis_string The XIOS string defining the axis
  !> @param[in] array_size The integer size of the array
  !> @return The axis definition array
  function numeric_axis_from_string(axis_string, array_size) result(axis_array)

    implicit none

    real(r_def), allocatable          :: axis_array(:)

    character(str_def)                :: axis_string
    integer(i_def)                    :: array_size

    logical(l_def)                    :: in_array = .false.
    character(str_def)                :: current_value = ''
    integer(i_def)                    :: char_index
    integer(i_def)                    :: array_index
    real(r_def)                       :: array_value

    allocate(axis_array(array_size))
    array_index = 1
    current_value = ''

    ! Loop through string representing axis
    do char_index = 1, len(trim(axis_string))

      ! Entering array of axis definition
      if (axis_string(char_index:char_index) == '[') then
        in_array = .true.

      ! Exiting arrray of axis definition, add last value to axis_array
      else if (axis_string(char_index:char_index) == ']') then
        in_array = .false.
        read(current_value, '(f10.0)') array_value
        axis_array(array_index) = array_value

      else if (in_array) then

        ! Append character to current value of array element
        if (axis_string(char_index:char_index) /= ' ') then
          current_value = trim(current_value) // axis_string(char_index:char_index)

        ! Spaces separate array elements so add current value to axis_array
        else
          read(current_value, '(f10.0)') array_value
          axis_array(array_index) = array_value
          array_index =  array_index + 1
          current_value = ''

        end if
      end if
    end do

  end function numeric_axis_from_string


  !===========================================================================
  !> Get an array of the labels from the string in the iodef.xml file
  !> @param[in] label_string The XIOS string defining the axis labels
  !> @param[in] array_size The integer size of the array
  !> @return The label definition array
  function label_axis_from_string(label_string, array_size) result(label_array)

    implicit none

    character(str_def), allocatable   :: label_array(:)

    character(str_def)                :: label_string
    integer(i_def)                    :: array_size

    logical(l_def)                    :: in_array = .false.
    character(str_def)                :: current_value = ''
    integer(i_def)                    :: char_index
    integer(i_def)                    :: array_index

    allocate(label_array(array_size))
    array_index = 1
    current_value = ''

    ! Loop through string representing axis
    do char_index = 1, len(trim(label_string))

      ! Entering array of label definition
      if (label_string(char_index:char_index) == '[') then
        in_array = .true.

      ! Exiting arrray of label definition, add last value to label_array
      else if (label_string(char_index:char_index) == ']') then
        in_array = .false.
        label_array(array_index) = current_value

      else if (in_array) then

        ! Append character to current value of array element
        if (label_string(char_index:char_index) /= ' ') then
          current_value = trim(current_value) // label_string(char_index:char_index)

        ! Spaces separate array elements so add current value to label_array
        else
          label_array(array_index) = current_value
          array_index =  array_index + 1
          current_value = ''

        end if
      end if
    end do

  end function label_axis_from_string

  !===========================================================================
  !> Get an array of all ids of the axes that make up a grid from its id
  !>
  !> Assumes a grid id with n axes of the form:
  !>     <axis_1_id>__<axis_2_id>__ ... <axis_n_id>__<domain_id>_grid
  !>
  !> @param[in] grid_id The id of the grid
  !> @return Character array of the id of each axis in the grid
  function axis_ids_from_grid_id(grid_id) result(axis_ids)

    implicit none

    character(str_def)              :: grid_id
    character(str_def), allocatable :: axis_ids(:)

    integer(i_def)                  :: array_size
    integer(i_def)                  :: char_index
    integer(i_def)                  :: array_index
    character(str_def)              :: current_value

    array_size = 0
    array_index = 1
    current_value = ''

    ! Loop through string to count how many axes are in name
    do char_index = 1, len(trim(grid_id)) - 1
      if (grid_id(char_index:char_index+1) == '__') then
        array_size = array_size + 1
      end if
    end do

    allocate(axis_ids(array_size))

    ! Loop through string representing axis
    do char_index = 1, len(trim(grid_id)) - 1

      ! Append character to current value of axis id
      if (grid_id(char_index:char_index+1) /= '__') then
        ! Ignore underscore at start of axis id from '__' separator
        if (current_value == '_') then
          current_value = grid_id(char_index:char_index)
        else
          current_value = trim(current_value) // grid_id(char_index:char_index)
        end if

      ! Double underscores separate axis ids so add current value to array
      else
        axis_ids(array_index) = current_value
        array_index =  array_index + 1
        current_value = ''

        ! Exit loop when we have added all the axis ids
        if (array_index > array_size) exit

      end if
    end do

  end function axis_ids_from_grid_id

end module fieldspec_xml_parser_mod
