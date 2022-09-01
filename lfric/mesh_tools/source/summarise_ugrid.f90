!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Summarise UGRID
!>
!> @brief   Utility to summarise of each mesh contained if NetCDF file which
!>          conforms to UGRID format convention.
!> @details Usage:
!>
!>          summarise_ugrid <filename>
!>          filename - UGRID file
!>
!-----------------------------------------------------------------------------
program summarise_ugrid

  use cli_mod,         only : get_initial_filename
  use constants_mod,   only : i_def, r_def, str_def, str_long, str_longlong, l_def
  use, intrinsic :: iso_fortran_env, only : output_unit
  use ncdf_quad_mod,   only : ncdf_quad_type
  use ugrid_2d_mod,    only : ugrid_2d_type
  use ugrid_file_mod,  only : ugrid_file_type
  use log_mod,         only : initialise_logging, finalise_logging, &
                              log_event, log_scratch_space,         &
                              LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use mpi_mod,         only : initialise_comm, store_comm, finalise_comm, &
                              get_comm_size, get_comm_rank

  implicit none

  character(:), allocatable :: filename

  class(ugrid_file_type), allocatable :: ugrid_file
  type(ugrid_2d_type)                 :: infile

  integer(i_def) :: n_meshes, i, j
  character(str_def), allocatable :: mesh_names(:)

  character(str_def) :: mesh_name
  character(str_def) :: geometry
  character(str_def) :: coord_sys
  character(str_def) :: topology

  logical(l_def)     :: periodic_x
  logical(l_def)     :: periodic_y

  integer(i_def)     :: max_stencil_depth
  character(str_longlong) :: constructor_inputs
  character(str_long)     :: target_mesh_names_str
  character(str_def), allocatable :: target_mesh_names(:)

  real(r_def)    :: north_pole(2)
  real(r_def)    :: null_island(2)

  character(str_def) :: fmt_str
  character(str_def) :: tmp_str

  integer(i_def) :: nodes, edges, faces
  integer(i_def) :: nodes_per_face, edges_per_face
  integer(i_def) :: nodes_per_edge, max_faces_per_node

  integer(i_def) :: comm, total_ranks, local_rank, nmaps

  ! Start up
  call initialise_comm(comm)
  call store_comm(comm)
  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()
  call initialise_logging(local_rank, total_ranks, "summarise")

  ! Get filename from command line
  call get_initial_filename( filename, 'UGRID mesh file' )

  ! Create object to manipulate UGRID conforming NetCDF file
  allocate(ncdf_quad_type::ugrid_file)
  call infile%set_file_handler(ugrid_file)

  ! Get the names of all mesh topologies in the UGRID file
  call infile%get_n_meshes( trim(filename), n_meshes )
  allocate( mesh_names(n_meshes) )
  call infile%get_mesh_names( trim(filename), mesh_names )

  call log_event(                                                         &
      '================================================================', &
      LOG_LEVEL_INFO )
  write (log_scratch_space,'(A)') &
      'File ('// trim(adjustl(filename))// ') contains ugrid mesh(es):'
  call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
  call log_event(                                                         &
      '================================================================', &
      LOG_LEVEL_INFO )

  do i=1, n_meshes

    ! Load in specified mesh from UGRID file into ugrid file object
    call infile%set_from_file_read( trim(mesh_names(i)), &
                                    trim(adjustl(filename)) )

    if (allocated(target_mesh_names)) deallocate(target_mesh_names)

    if (n_meshes > 1) then
      ! Extract data on the current mesh in the ugrid file object
      call infile%get_metadata(                                &
                      mesh_name          = mesh_name,          &
                      geometry           = geometry,           &
                      topology           = topology,           &
                      coord_sys          = coord_sys,          &
                      max_stencil_depth  = max_stencil_depth,  &
                      constructor_inputs = constructor_inputs, &
                      nmaps              = nmaps,              &
                      target_mesh_names  = target_mesh_names,  &
                      periodic_x         = periodic_x,         &
                      periodic_y         = periodic_y,         &
                      north_pole         = north_pole,         &
                      null_island        = null_island )
    else
      call infile%get_metadata(                                &
                      mesh_name          = mesh_name,          &
                      geometry           = geometry,           &
                      topology           = topology,           &
                      coord_sys          = coord_sys,          &
                      max_stencil_depth  = max_stencil_depth,  &
                      constructor_inputs = constructor_inputs, &
                      periodic_x         = periodic_x,         &
                      periodic_y         = periodic_y,         &
                      north_pole         = north_pole,         &
                      null_island        = null_island )
    end if

    call infile%get_dimensions( nodes, edges, faces,            &
                                nodes_per_face, edges_per_face, &
                                nodes_per_edge, max_faces_per_node )

    ! Write extracted data and log output
    write (log_scratch_space, '(A)') &
      '"'//trim(mesh_names(i))//'":'
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    fmt_str='(A,T24,A)'
    write ( log_scratch_space, fmt_str ) &
        '  Geometry: ', trim(geometry)
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) &
        '  Topology: ', trim(topology)
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) &
        '  Co-ordinate system: ', trim(coord_sys)
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    fmt_str='(A,T24,L1)'
    write ( log_scratch_space, fmt_str ) &
        '  Periodic X: ', periodic_x
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) &
        '  Periodic Y: ', periodic_y
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    fmt_str='(A,T24,A)'
    write ( log_scratch_space, fmt_str ) &
        '  Constructor inputs: ', trim(constructor_inputs)
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    fmt_str='(A,T24,I0)'
    write ( log_scratch_space, fmt_str ) '  Nodes:', nodes
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) '  Edges:', edges
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) '  Faces:', faces
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) &
        '  Nodes per face:', nodes_per_face
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) &
        '  Edges per face:', edges_per_face
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) &
        '  Nodes per edge:', nodes_per_edge
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    write ( log_scratch_space, fmt_str ) &
        '  Max. faces per node:', max_faces_per_node
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )


    if (allocated(target_mesh_names)) then
      if (size(target_mesh_names) > 0 .and. nmaps > 0) then
        target_mesh_names_str = '"'//trim(adjustl(target_mesh_names(1)))//'"'
        if (size(target_mesh_names) > 1) then
          do j=2, size(target_mesh_names)
            write(target_mesh_names_str,'(A)')                    &
                 trim(adjustl(target_mesh_names_str)) // ', "' // &
                 trim(target_mesh_names(j)) // '"'
          end do
        end if
      end if

      fmt_str='(A,T24,A)'
      write ( log_scratch_space, fmt_str ) &
           '  Maps to:', trim(adjustl(target_mesh_names_str))
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
      deallocate(target_mesh_names)

    end if

    if ( trim(geometry)  == 'spherical' .and. &
         trim(coord_sys) == 'll' ) then

      fmt_str = '("[",F10.2,",",F10.2,"]")'

      write(tmp_str, fmt_str) north_pole
      write( log_scratch_space, '(A)' ) &
         '  North pole [lon,lat]: '//trim(tmp_str)
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

      write(tmp_str, fmt_str) null_island
      write( log_scratch_space, '(A)' ) &
         '  Null Island [lon,lat]: '//trim(tmp_str)
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    end if

  end do

  call log_event(                                                         &
      '================================================================', &
      LOG_LEVEL_INFO )

  deallocate( mesh_names )

  ! Finalise the logging system
  call finalise_logging()

end program summarise_ugrid
