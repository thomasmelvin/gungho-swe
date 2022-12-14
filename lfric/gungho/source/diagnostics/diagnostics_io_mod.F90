!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Module implementing a basic diagnostic system
!!
!!  @details Module implementing a basic diagnostic system
!-------------------------------------------------------------------------------
module diagnostics_io_mod

  use clock_mod,                     only: clock_type
  use constants_mod,                 only: i_def, i_timestep, &
                                           r_def, str_max_filename
  use physics_mappings_alg_mod,      only: split_wind_alg
  use diagnostic_alg_mod,            only: extract_w2h_diagnostic_alg,  &
                                           scalar_nodal_diagnostic_alg, &
                                           scalar_ugrid_diagnostic_alg, &
                                           vector_nodal_diagnostic_alg
  use io_config_mod,                 only: use_xios_io, write_fluxes
  use files_config_mod,              only: diag_stem_name
  use function_space_collection_mod, only: function_space_collection
  use finite_element_config_mod,     only: element_order
  use fs_continuity_mod,             only: Wtheta, W2H
  use project_output_mod,            only: project_output
  use io_mod,                        only: ts_fname, &
                                           nodal_write_field
  use lfric_xios_write_mod,          only: write_field_face, &
                                           write_field_edge
  use mesh_mod,                      only: mesh_type
  use field_mod,                     only: field_type
  use field_parent_mod,              only: write_interface
  use fs_continuity_mod,             only: W3

  implicit none
  private
  public :: write_scalar_diagnostic, &
            write_vector_diagnostic

contains

!-------------------------------------------------------------------------------
!>  @brief    Handles generic scalar diagnostic processing
!!
!!  @details  Handles generic scalar diagnostic processing
!!
!!> @param[in] field_name  Character string the field name
!!> @param[in] field       The field to output
!!> @param[in] ts          Timestep
!!> @param[in] mesh        Mesh
!!> @param[in] W3_project  Logical to allow projection to W3
!-------------------------------------------------------------------------------

subroutine write_scalar_diagnostic( field_name, field, &
                                    clock, mesh, W3_project )
  implicit none

  character(len=*),  intent(in)    :: field_name
  type(field_type),  intent(in)    :: field
  class(clock_type), intent(in)    :: clock
  type(mesh_type),   intent(in), pointer :: mesh
  logical,           intent(in)    :: W3_project

  integer(i_def), parameter       :: nodal_output_unit = 21

  ! Local Variables
  type(field_type)                :: nodal_coordinates(3)
  type(field_type)                :: output_field(3)
  type(field_type)                :: level
  character(len=str_max_filename) :: fname
  integer(i_timestep)             :: timestep

  procedure(write_interface), pointer  :: tmp_write_ptr => null()

  ! Nodal output
  if ( .not. (use_xios_io) )  then

    if (clock%is_initialisation()) then
      timestep = 0
    else
      timestep = clock%get_step()
    end if

    ! Always call straight nodal output

    ! Setup output filename

    fname=trim(ts_fname( trim(diag_stem_name), &
                         "nodal_",             &
                         field_name,           &
                         timestep,             &
                         ".m" ))

    ! Call diagnostic processing to create nodal field
    call scalar_nodal_diagnostic_alg( output_field, nodal_coordinates, &
                                      level, field_name, field,        &
                                      mesh, .false. )

    ! Call write routine
    call nodal_write_field(nodal_coordinates, level, output_field,    &
                           1, nodal_output_unit, fname)

    ! If projection to W3 was requested then output that as well

    if (W3_project .and. (field%which_function_space() /= W3)) then

      fname=trim(ts_fname( trim(diag_stem_name),  &
                           "nodal_w3projection_", &
                           field_name,            &
                           timestep,              &
                           ".m" ))


      ! Call diagnostic processing to create nodal field
      call scalar_nodal_diagnostic_alg( output_field, nodal_coordinates, &
                                        level, field_name, field,        &
                                        mesh, .true. )

      ! Call write routine
      call nodal_write_field( nodal_coordinates, level, output_field,    &
                              1, nodal_output_unit, fname)

    end if

  else ! (use_xios_io)

    ! XIOS UGRID output

    ! Call diagnostic processing to create ugrid field
    call scalar_ugrid_diagnostic_alg( output_field, field_name, field, &
                                      mesh, .false. )

    ! Set a field I/O method appropriately
    tmp_write_ptr => write_field_face
    call output_field(1)%set_write_behaviour(tmp_write_ptr)

    ! Call write on the output field

    ! Check if we need to write an initial field
    if (clock%is_initialisation()) then
       call output_field(1)%write_field(trim('init_'//field_name))
    else
       call output_field(1)%write_field(trim(field_name))
    end if

    nullify(tmp_write_ptr)

  end if


end subroutine write_scalar_diagnostic

!-------------------------------------------------------------------------------
!>  @brief    Handles generic vector diagnostic processing
!!
!!  @details  Handles generic vector diagnostic processing
!!
!!> @param[in] field_name  Character string the field name
!!> @param[in] field       The field to output
!!> @param[in] ts          Timestep
!!> @param[in] mesh        Mesh
!!> @param[in] W3_project  Logical to allow projection to W3
!-------------------------------------------------------------------------------

subroutine write_vector_diagnostic( field_name, field, &
                                    clock, mesh, W3_project )
  implicit none

  character(len=*),  intent(in)    :: field_name
  type(field_type),  intent(in)    :: field
  class(clock_type), intent(in)    :: clock
  type(mesh_type),   intent(in), pointer :: mesh
  logical,           intent(in)    :: W3_project

  integer(i_def), parameter       :: nodal_output_unit = 21

  ! Local Variables
  type(field_type)                :: nodal_coordinates(3)
  type(field_type)                :: output_field(3)
  type(field_type)                :: projected_field(3)
  type(field_type)                :: level
  type(field_type)                :: u1_wind, u2_wind, u3_wind
  type(field_type)                :: h_wind, v_wind
  character(len=str_max_filename) :: fname
  character(len=1)                :: uchar
  character(len=str_max_filename) :: field_name_new
  integer(i_def)                  :: i
  integer(i_def)                  :: output_dim
  integer(i_timestep)             :: timestep

  procedure(write_interface), pointer  :: tmp_write_ptr => null()

  ! Nodal output
  if ( .not. (use_xios_io) )  then

    if (clock%is_initialisation()) then
      timestep = 0
    else
      timestep = clock%get_step()
    end if

    ! Always call straight nodal output

    ! Setup output filename
    fname=trim(ts_fname(trim(diag_stem_name), &
                        "nodal_",             &
                        field_name,           &
                        timestep,             &
                        ".m"))

    call vector_nodal_diagnostic_alg( output_field, output_dim, &
                                      nodal_coordinates, level, &
                                      field_name, field )

    ! Call write routine
    call nodal_write_field(nodal_coordinates, level, output_field, &
                           output_dim, nodal_output_unit, fname)

    ! If projection to W3 was requested then output that as well

    if (W3_project) then

      output_dim = 3

      ! Project the field to the output field
      call project_output( field, projected_field, output_dim, W3 , mesh )

      do i =1,output_dim
         ! Write the component number into a new field name
         write(uchar,'(i1)') i
         field_name_new = trim("w3projection_"//field_name//uchar)

         ! Setup output filename
         fname=trim(ts_fname(trim(diag_stem_name), &
                             "nodal_",             &
                             field_name_new,       &
                             timestep,             &
                             ".m"))

         ! Call scalar output on each component
         call scalar_nodal_diagnostic_alg( output_field(i), nodal_coordinates,        &
                                           level, field_name_new, projected_field(i), &
                                           mesh, .false.)

         ! Call write routine
         call nodal_write_field(nodal_coordinates, level, output_field(i), &
                                1, nodal_output_unit, fname)
       end do


    end if

  else

  ! XIOS UGRID output

    ! Check for specific vector fields and applying appropriate processing

    ! Currently we need to force projection of vorticity (Xi) to W3 until
    ! we decide how it should be handled in UGRID
    if (field_name == 'xi') then

      output_dim = 3

      ! Project the field to the output field
      call project_output( field, projected_field, output_dim, W3 , mesh )

      ! Set up correct I/O handler for Xi projected to W3
      tmp_write_ptr => write_field_face

      do i =1,output_dim

        call projected_field(i)%set_write_behaviour(tmp_write_ptr)

        ! Write the component number into a new field name
        write(uchar,'(i1)') i
        field_name_new = trim(field_name//uchar)

        ! Check if we need to write an initial field
        if (clock%is_initialisation()) then
          call projected_field(i)%write_field(trim('init_'//field_name_new))
        else
          call projected_field(i)%write_field(trim(field_name_new))
        end if

     end do

   else ! wind fields in w2

      !---- Output wind as u1, u2, and u3 wind components ---

      call u1_wind%initialise( vector_space = &
               function_space_collection%get_fs(mesh, element_order, W2H) )

     call u2_wind%initialise( vector_space = &
              function_space_collection%get_fs(mesh,element_order, W2H) )

     call u3_wind%initialise( vector_space = &
              function_space_collection%get_fs(mesh,element_order, Wtheta) )

      call split_wind_alg( u1_wind, u2_wind, u3_wind, &
                           field, mesh )

      ! Set up I/O handler as these are derived fields
      tmp_write_ptr => write_field_face
      call u3_wind%set_write_behaviour(tmp_write_ptr)
      tmp_write_ptr => write_field_edge
      call u1_wind%set_write_behaviour(tmp_write_ptr)
      call u2_wind%set_write_behaviour(tmp_write_ptr)

      if (clock%is_initialisation()) then
        if (field_name == 'u') then
          call u1_wind%write_field("init_u_in_w2h")
          call u2_wind%write_field("init_v_in_w2h")
          call u3_wind%write_field("init_w_in_wth")
        else
          call u1_wind%write_field("init_"//trim(field_name)//"1")
          call u2_wind%write_field("init_"//trim(field_name)//"2")
          call u3_wind%write_field("init_"//trim(field_name)//"3")
        end if
      else
        if (field_name == 'u') then
          call u1_wind%write_field("u_in_w2h")
          call u2_wind%write_field("v_in_w2h")
          call u3_wind%write_field("w_in_wth")
        else
          call u1_wind%write_field(trim(field_name)//"1")
          call u2_wind%write_field(trim(field_name)//"2")
          call u3_wind%write_field(trim(field_name)//"3")
        end if
      end if

      if (write_fluxes) then

        !---- Output wind as w2h and wtheta fluxes ----

        ! Convert u to w2h (h_wind) and wtheta (v_wind)
        call extract_w2h_diagnostic_alg( h_wind, v_wind, field )

        tmp_write_ptr => write_field_face
        call v_wind%set_write_behaviour(tmp_write_ptr)
        tmp_write_ptr => write_field_edge
        call h_wind%set_write_behaviour(tmp_write_ptr)
        if (clock%is_initialisation()) then
          call h_wind%write_field( "init_h_"//trim(field_name) )
          call v_wind%write_field( "init_v_"//trim(field_name) )
        else
          call h_wind%write_field( "h_"//trim(field_name) )
          call v_wind%write_field( "v_"//trim(field_name) )
        end if

      end if ! Output fluxes

    end if ! Check for wind fields

    nullify(tmp_write_ptr)

  end if


end subroutine write_vector_diagnostic


end module diagnostics_io_mod
