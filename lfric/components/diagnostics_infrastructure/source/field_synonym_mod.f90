!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the field_synonym class
!>
!> @details A key-value pair object to contain synonyms of a field.
!> Keys are enumerations (interger)

module field_synonym_mod

  use constants_mod,          only: str_def, i_native

  implicit none

  private

  type, public :: field_synonym_type

    !> @brief String containing the name of the type of data being stored
    !>
    !> This should come from a model field_synonyms_enum_mod, such as
    !> `um_physics/source/diagnostics_meta/meta_type/field_synonyms_enum_mod.f90`
    integer(i_native) :: type
    !> @brief Value of the data
    character(str_def) :: value

  end type field_synonym_type

  interface field_synonym_type
      module procedure field_synonym_constructor
  end interface field_synonym_type

contains
  !> Construct a <code>field_synonym_type</code> object.
  !>
  !> @param [in] type Describes what type of synonym is being stored here
  !> @param [in] value The actual data to be stored

  !> @return self the meta_data object
  !> Constructor named shortened for use in meta_mod
  function field_synonym_constructor(type, value) result(self)

    implicit none

    integer(i_native),           intent(in) :: type
    character(*),                intent(in) :: value

    type(field_synonym_type) :: self

    self%type = type
    self%value = value

  end function field_synonym_constructor

end module field_synonym_mod