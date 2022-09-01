!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief A module for the meta_data type
!>
!> @details This type will hold all meta data information for a given field.
!> Objects of this type are created by the field meta data definition files

module diagnostics_mod

  use constants_mod,                only: i_def, str_def, l_def, str_long, &
                                          i_native
  use vertical_dimension_types_mod, only: range_vertical_dimension_meta_data_type, &
                                          list_vertical_dimension_meta_data_type, &
                                          abstract_vertical_meta_data_type
  use non_spatial_dimension_mod,    only: non_spatial_dimension_type
  use misc_meta_data_mod,           only: misc_meta_data_type
  use field_synonym_mod,            only: field_synonym_type

  implicit none
  private

  type, public :: field_meta_data_type

    private
    !> @brief Unique id used by the diagnostic system to identify the field
    !>
    !>     unique_id = "<section>__<field>"
    !>
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span> <span style="font-weight:bold;font-size: 0.8em;
    !> color:#0f79be">NetCDF Metadata</span> <span style="font-weight:bold;
    !> font-size: 0.8em; color:#CF4922">User Configuration Information</span>
    !>
    !> This is the unique identifier for this prog/diagnostic throughout the
    !> whole LFRic system and in output files.
    !>
    !> This **must** be globally unique, but given the inclusion of the
    !> `<section>` the field only needs to be unique within that section. The
    !> `<field> `aspect *should* avoid acronyms and abbreviations if possible.
    !> The `<field>` _might_ be the same as the name used in external
    !> standards, if that is specific and clear enough. It *should* be clear
    !> to those without domain expertise what this is, although the description
    !> can be used to further expand on this.
    !>
    !> This is in lower_snake_case format with the two parts separated with
    !> double underscore.
    !>
    !> e.g. `large_scale_precipitation__rainfall_rate`
    character(str_def)  ::  unique_id

    !> @brief Optional(ish). Long name of the field. Often specified by an
    !> external standard such as CF.
    !>
    !>     long_name = "..."
    !> <span style="font-weight:bold;font-size: 0.8em; color:#0f79be">NetCDF
    !> Metadata</span>
    !>
    !> Space for a longer name - particularly if acronyms or abbreviations have
    !> been included in the unique_id. If no expansion is needed then the
    !> `<field>` aspect of the unique_id can be reused.
    !> This is often used in lieu of a standard_name if the field has not yet
    !> been adopted by a standard but this may also be defined by external
    !> conventions.
    !>
    !> Convention states you __should__ have at least one of long_name OR
    !> standard_name.
    !>
    !> This is in lower_snake_case format.
    !>
    !> e.g. `rainfall_rate_from_large_scale_precipitation`
    character(str_def)  ::  long_name


    !> @brief The SI unit of measure for the field. Often specified by an
    !> external standard such as CF.
    !>
    !>     units = "K"
    !> <span style="font-weight:bold;font-size: 0.8em; color:#0f79be">NetCDF
    !> Metadata</span>
    !>
    !> The unit of measure for the field. Care must be taken to confirm
    !> character cases and spacing are accurate. IF CF/CMIP synonyms are
    !> specified then this can be validated against those standards.
    !>
    !> It will be very easy to end up with a multitude of mis-matched but
    !> identical units. For example the SI unit of speed could be:
    !> `m/s` OR `m / s` OR `m s-1` OR `m.s-1` OR `m s^-1` OR `m.s^-1`.
    !> The current preference is to avoid special characters (`/` and `^`)
    !> e.g. `m s-1` for speed.
    !>
    !> SI derived units are encouraged for ease of reading.
    character(str_def)  ::  units

    !> @brief Function space.
    !> For prognostics, which function space does this need to exist on.
    !> For diagnostics, what function space is it possible to request this on
    !>
    !>     function_space = W3
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span>
    !>
    !> This should be an enumeration from \link fs_continuity_mod.F90
    !> fs_continuity_mod \endlink. This specifies the function space the field
    !> operates on, and therefore (in combination with order) the size of the
    !> field and the distribution of the data points. For more information on
    !> function spaces see LFRic introductory documentation
    integer(i_native)   ::  function_space

    !> @brief The IO driver used by this field
    !>
    !>     io_driver = "write_field_face"
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span>
    !>
    !> This specifies the function by which LFRic communicates with XIOS when
    !> transmitting the field. It is coupled to Order and Function Space and
    !> will be automatically calculated in a future revision.
    character(str_def)  ::  io_driver

    !> @brief The order of the function space. Cannot be negative
    !>
    !>     order = 0
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span>
    !>
    !> This is an integer specifying the default order of the finite element
    !> analysis for the field. Normally this is 0th order (single value /
    !> mesh entity).
    integer(i_def)      ::  order

    !> @brief Contains Rose triggering syntax.
    !>
    !>     trigger = "__checksum: true;"
    !><span style="font-weight:bold;font-size: 0.8em; color:#CF4922">User
    !> Configuration Information</span>
    !>
    !> This is Rose triggering syntax and allows the field to block others. The
    !> typical starting value is "__checksum: true;" to automatically
    !> show/hide the optional checksum field when the field is turned on.
    !>
    !> The Rose triggering syntax requires that the restriction be defined on
    !> the blocking item, to release the blocked item. This will require
    !> changing metadata in sections beyond the field being edited if it has
    !> a prerequisite outside of the section. It is hoped that future work
    !> will expand this with a mechanism to define prerequisites.
    !>
    !> The fields are identified in Rose using their unique_id (with section
    !> removed), section and group in the format
    !>     field_config:<section>:<group>=<field>:true
    !> e.g.
    !>     trigger = "__checksum: true;" // &
    !>               " field_config:radiation:general=radiative_flux:true;"
    !> would only allow the field radiative_flux in the general group within
    !> the radiation section if the field this trigger is attached to is
    !> enabled.
    !>
    !> It is also possible to reference these fields from within the science
    !> configuration and vice versa - It will all be stitched together as one
    !> Rose app. e.g. the science option to specify the calculation method for
    !> a particular field is only available if that field is enabled.
    !>
    !> The triggering syntax should be used to ensure that user configurations
    !> are valid (within the limits of the Rose syntax).
    character(str_def)  ::  trigger

    !> @brief This string is displayed in the Rose gui under the field
    !>
    !>     description = "something helpful about what this field is / does"
    !> <span style="font-weight:bold;font-size: 0.8em; color:#CF4922">User
    !> Configuration Information</span>
    !>
    !> The description property is primarily for help text in Rose. The user
    !> will already have access to the field's Name, Unit of Measure, Function
    !> Space, Data Type, Time Step, Interpolation, Vertical Dimensions,
    !> Synonyms, and Non Spatial Dimensions so these don't need repeating.
    character(str_long) ::  description

    !> @brief Enumerators are used to represent the data type
    !>
    !>     data_type = REAL_TYPE
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span>
    !>
    !> This allows specification of the type of field to create - Real
    !> vs Integer. These are enumerations pulled from the
    !> \link constants_mod constants mod \endlink.
    integer(i_native)   ::  data_type

    !> @brief Timestep information.
    !> Specifies which timesteps diagnostics or prognostics are available on.
    !> This is set to an enumerated value, not the actual number of time steps
    !>
    !>     time_step = STANDARD_TIMESTEP
    !> <span style="font-weight:bold;font-size: 0.8em; color:#CF4922">User
    !> Configuration Information</span> <span style="font-weight:bold;
    !> font-size: 0.8em; color:#288172">LFRic Field Definition</span>
    !>
    !> This property is used to specify if fields are only available to be
    !> calculated on particular timestep patterns. These are registered using
    !> the model's local time_step_enum_mod.f90 file
    !>
    !> This feature is yet to be implemented fully.
    integer(i_native)   ::  time_step

    !> @brief Interpolation method. Yet to be implemented
    !>
    !>     recommended_interpolation = BILINEAR
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span> <span style="font-weight:bold;font-size: 0.8em;
    !> color:#0f79be">NetCDF Metadata</span>
    !>
    !> This property is used to specify the suggested interpolation method.
    !> These are registered using the model's local interpolation_enum_mod e.g.
    !> `um_physics/source/diagnostics_meta/meta_types/interpolation_enum_mod.f90`
    !>
    !> This feature is yet to be implemented fully.
    integer(i_native)   ::  recommended_interpolation

    !> @brief Level of packing for the field
    !>
    !>      packing = 0
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span> <span style="font-weight:bold;font-size: 0.8em;
    !> color:#0f79be">NetCDF Metadata</span>
    !>
    !> This property is a placeholder for specifying the packing approach to
    !> take for a particular field. This feature is yet to be implemented
    !> fully.
    integer(i_def)      ::  packing

    !> @brief Optional. Vertical dimension specification
    !>
    !>     vertical_dimension = model_height_dimension(...)
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span>
    !>
    !> A field_meta_data_type can be supplied with, as an optional argument, a
    !> vertical_dimension_type to define its vertical dimension meta data.
    !> Current types are:
    !>  - \link
    !> vertical_dimension_types_mod::list_vertical_dimension_meta_data_type
    !> list vertical dimension \endlink
    !>  - \link
    !> vertical_dimension_types_mod::range_vertical_dimension_meta_data_type
    !> range vertical dimension \endlink
    !>
    !> Most fields (including some 2d fields) will have a vertical dimension.
    !> This defines where the data is located vertically. This has a wide
    !> range of options but all are based on an \link
    !> vertical_dimension_types_mod::abstract_vertical_meta_data_type
    !> abstract_vertical_meta_data_type \endlink. This breaks down into two
    !> patterns:
    !> - When is it defined?
    !> - Model: definition is relative to a vertical axis defined as part of
    !>    the model configuration - e.g. standard vertical axis between 1st
    !>    atmospheric level and top level in boundary layer. This constrains
    !>    but defers the precise specification until user configuration.
    !> - Fixed: definition is hard coded by the developer - e.g. cloud amount
    !>    between 700-950m, 950-1100m and 1100-1300m.
    !>
    !> The vertical dimension also contains information about what standard
    !> name it is, the units and which direction positive is. These are
    !> wrapped in some standard configurations for common use cases - found
    !> in the local models vertical_dimensions_mod - e.g.
    !> `um_physics/source/diagnostics_meta/meta_types/
    !> vertical_dimensions_mod.f90`
    !> These lock down direction, units and standard name to allow just the
    !> axis values or bounds to be specified - e.g. model_depth_dimension
    !> specifies
    !>
    !>     standard_name='depth', units='m', positive=POSITIVE_DOWN
    !>
    !> For model dimensions the bounds are specified using standard level
    !> labels. These are specified in the local model levels_enum_mod file e.g.
    !> `um_physics/source/diagnostics_meta/meta_types/levels_enum_mod.f90`
    !>
    !> For fixed dimensions the bounds are specified as an array of values.
    class(abstract_vertical_meta_data_type), allocatable :: vertical_dimension

    !> @brief Optional. Non spatial dimension definitions
    !>
    !>     non_spatial_dimension = [non_spatial_dimension_type(...)]
    !> <span style="font-weight:bold;font-size: 0.8em; color:#288172">LFRic
    !> Field Definition</span>
    !>
    !> A field can have additional dimensions beyond the vertical. This
    !> allows for the field to contain additional data points such as
    !> wavelength or particle size.
    !> see \link non_spatial_dimension_mod::non_spatial_dimension_type non
    !> spatial dimension type \endlink for more details
    type(non_spatial_dimension_type), allocatable :: non_spatial_dimension(:)

    !> @brief Optional(ish). Standard name for the field. Often specified by an
    !> external standard such as CF.
    !>
    !>      standard_name = "..."
    !> <span style="font-weight:bold;font-size: 0.8em; color:#0f79be">NetCDF
    !> Metadata</span>
    !>
    !> If the prog/diagnostic is synonymous with a matching CF standard then
    !> standard_name can be used to populate the output metadata in the netcdf
    !> file.
    !>
    !> This is in lower_snake_case format but should conform with the standard.
    !>
    !> Convention states you __should__ have at least one of long_name OR
    !> standard_name.
    !>
    !> e.g. `rainfall_rate`
    character(str_def)  ::  standard_name

    !> @brief Optional. Positive direction of the field (vectors), optional
    !>
    !>     positive = "east"
    !> <span style="font-weight:bold;font-size: 0.8em; color:#0f79be">NetCDF
    !> Metadata</span>
    !>
    !> For fields that represent a vector to specify which direction a
    !> positive value is in (e.g. eastward wind: positive values are eastwards)
    integer(i_native)   ::  positive

    !> @brief Optional. Any synonyms for the field from list of known standards
    !>
    !>     synonyms = [field_synonym_type(...)]
    !> <span style="font-weight:bold;font-size: 0.8em; color:#CF4922">User
    !> Configuration Information</span> <span style="font-weight:bold;
    !> font-size: 0.8em; color:#0f79be">NetCDF Metadata</span>
    !>
    !> These are alternative identifiers for the field from specific standards.
    !> They take the form of (ENUM, Value) pairs within a field_synonym_type.
    !> more detail can be found in
    !> \link field_synonym_mod::field_synonym_type field synonym mod\endlink
    type(field_synonym_type), allocatable  :: synonyms(:)

    !> @brief Optional. Any misc data for the field
    !>
    !>     misc_meta_data = [&
    !>              misc_meta_data_type("source_satellite_data","EUMETSAT"),&
    !>              misc_meta_data_type(&
    !>                       "algorithm",&
    !>                       "https://doi.org/10.1145/3377713.3377808"&
    !>                       )
    !>              ]
    !> <span style="font-weight:bold;font-size: 0.8em; color:#0f79be">NetCDF
    !> Metadata</span>
    !>
    !> This is a list of misc_meta_data_type, which in turn is a pair of
    !> strings as key / value.
    !>
    !> This field is primarily to pass the key value pairs into the netcdf
    !> metadata. This allows you to add extra information to your output
    !> fields as needed. This might be references to standards that are not
    !> incorporated, extra information that's needed by a standard but not
    !> included by other properties, or any other helpful information as
    !> required. This information is also made available as part of the
    !> metadata json file for downstream system consumption.
    type(misc_meta_data_type), allocatable  :: misc_meta_data(:)

  contains

    !> @brief Getter for unique_id
    !> @param[in]  self  field_meta_data_type
    !> @return unique_id string
    procedure, public :: get_unique_id !> Returns the unique id string

  end type field_meta_data_type

  interface field_meta_data_type
    module procedure meta_data_constructor
  end interface

contains

  !> Construct a <code>meta_data_type</code> object.
  !> \fn meta_data_constructor
  !> @brief The constructor for a field_meta_data_type object
  !> @param [in] unique_id A unique identifier of the field
  !> @param [in] units SI Unit of measure for the field
  !> @param [in] function_space The function_space to create field with
  !> @param [in] order The order of the function space
  !> @param [in] io_driver The IO driver used by the field
  !> @param [in] trigger The triggering syntax for the field in Rose
  !> @param [in] description Description of the field shown in Rose
  !> @param [in] data_type The data type of the field
  !> @param [in] time_step The available time step of the field
  !> @param [in] recommended_interpolation The recommended interpolation method
  !> @param [in] packing Packing setting for the field
  !> @param [in,optional] standard_name The standard name of the field if it exists
  !> @param [in,optional] long_name The long name of the field
  !> @param [in,optional] positive The direction for positive numbers
  !> @param [in,optional] vertical_dimension The vertical dimension of the field
  !> @param [in,optional] non_spatial_dimension The non-spatial dimension(s) of the field
  !> @param [in,optional] synonyms Holds a key/value pair of field_synonyms_enum / string
  !> @param [in,optional] misc_meta_data Holds a key/value pair of strings
  !> used for any miscellaneous data that a field might need
  !> @return self the meta_data object
  !>
  function meta_data_constructor(unique_id,                     &
                                 units,                         &
                                 function_space,                &
                                 order,                         &
                                 io_driver,                     &
                                 trigger,                       &
                                 description,                   &
                                 data_type,                     &
                                 time_step,                     &
                                 recommended_interpolation,     &
                                 packing,                       &
                                 standard_name,                 &
                                 long_name,                     &
                                 positive,                      &
                                 vertical_dimension,            &
                                 non_spatial_dimension,         &
                                 synonyms,                      &
                                 misc_meta_data)                &
                                 result(self)

    implicit none

    character(*),                        intent(in) :: unique_id
    character(*),                        intent(in) :: units
    integer(i_native),                   intent(in) :: function_space
    integer(i_def),                      intent(in) :: order
    character(*),                        intent(in) :: io_driver
    character(*),                        intent(in) :: trigger
    character(*),                        intent(in) :: description
    integer(i_native),                   intent(in) :: data_type
    integer(i_native),                   intent(in) :: time_step
    integer(i_native),                   intent(in) :: recommended_interpolation
    integer(i_def),                      intent(in) :: packing
    character(*), optional,              intent(in) :: standard_name
    character(*), optional,              intent(in) :: long_name
    integer(i_native), optional,         intent(in) :: positive
    type(misc_meta_data_type), optional, intent(in) :: misc_meta_data(:)
    type(field_synonym_type), optional,  intent(in) :: synonyms(:)
    class(abstract_vertical_meta_data_type), optional, &
                                         intent(in) :: vertical_dimension
    type(non_spatial_dimension_type),        optional, &
                                         intent(in) :: non_spatial_dimension(:)

    type(field_meta_data_type)   :: self

    self%unique_id                 = unique_id
    self%units                     = units
    self%function_space            = function_space
    self%order                     = order
    self%io_driver                 = io_driver
    self%trigger                   = trigger
    self%description               = description
    self%data_type                 = data_type
    self%time_step                 = time_step
    self%recommended_interpolation = recommended_interpolation

    !> Handling optional arguments
    if(present(standard_name)) then
      self%standard_name = standard_name
    else
      self%standard_name = ""
    end if

    if(present(long_name)) then
      self%long_name = long_name
    else
      self%long_name = ""
    end if

    if(present(positive)) then
      self%positive = positive
    else
      self%positive = 0
    end if

    if(present(vertical_dimension)) then
      !> This is the way polymorphism is done in fortran 2003. You have to
      !> allocate the derived type instead of using the assignment operator (=)
      !> Fortran 2008 would support self%vertical_dimension = vertical_dimension
      allocate(self%vertical_dimension, source=vertical_dimension)
    end if

    if(present(non_spatial_dimension)) then
      allocate(self%non_spatial_dimension, source=non_spatial_dimension)
    end if

    if(present(misc_meta_data)) then
      allocate(self%misc_meta_data, source=misc_meta_data)
    end if

    if(present(synonyms)) then
      allocate(self%synonyms, source=synonyms)
    end if

  end function meta_data_constructor

  !> @brief Getter for unique_id
  !> @param[in]  self  field_meta_data_type
  !> @return unique_id string
  function get_unique_id(self) result(unique_id)

    implicit none

    class(field_meta_data_type), intent(in) :: self
    character(str_def) :: unique_id

    unique_id = trim(self%unique_id)

  end function get_unique_id

end module diagnostics_mod
