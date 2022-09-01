!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Planar mesh generator
!>
!> @brief   Utility to generate a planar surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          planar_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program planar_mesh_generator

  use cli_mod,           only: get_initial_filename
  use constants_mod,     only: i_def, l_def, r_def, str_def, i_native, &
                               cmdi, imdi, emdi
  use configuration_mod, only: read_configuration, final_configuration
  use gen_lbc_mod,       only: gen_lbc_type
  use gen_planar_mod,    only: gen_planar_type,          &
                               set_partition_parameters, &
                               NPANELS
  use rotation_mod,      only: get_target_north_pole,    &
                               get_target_null_island
  use global_mesh_mod,   only: global_mesh_type
  use halo_comms_mod,    only: initialise_halo_comms, &
                               finalise_halo_comms
  use io_utility_mod,    only: open_file, close_file
  use local_mesh_mod,    only: local_mesh_type
  use log_mod,           only: initialise_logging,       &
                               finalise_logging,         &
                               log_event, log_set_level, &
                               log_scratch_space,        &
                               LOG_LEVEL_INFO,           &
                               LOG_LEVEL_ERROR

  use mpi_mod,           only: initialise_comm, store_comm, finalise_comm, &
                               get_comm_size, get_comm_rank
  use ncdf_quad_mod,     only: ncdf_quad_type
  use partition_mod,     only: partition_type, partitioner_interface
  use reference_element_mod, only: reference_element_type, &
                                   reference_cube_type
  use remove_duplicates_mod, only: any_duplicates
  use ugrid_2d_mod,          only: ugrid_2d_type
  use ugrid_file_mod,        only: ugrid_file_type
  use ugrid_mesh_data_mod,   only: ugrid_mesh_data_type

  ! Configuration modules
  use mesh_config_mod,         only: mesh_filename, rotate_mesh, &
                                     n_meshes, mesh_names,       &
                                     mesh_maps, n_partitions,    &
                                     coord_sys, coord_sys_ll,    &
                                     coord_sys_xyz,              &
                                     key_from_coord_sys,         &
                                     topology,                   &
                                     topology_non_periodic,      &
                                     topology_periodic,          &
                                     topology_channel,           &
                                     key_from_topology,          &
                                     geometry, geometry_planar,  &
                                     geometry_spherical,         &
                                     key_from_geometry
  use partitioning_config_mod, only: max_stencil_depth
  use planar_mesh_config_mod,  only: edge_cells_x, edge_cells_y, &
                                     periodic_x, periodic_y,     &
                                     domain_x, domain_y,         &
                                     first_node,                 &
                                     create_lbc_mesh,            &
                                     lbc_rim_depth,              &
                                     lbc_parent_mesh
  use rotation_config_mod,     only: target_north_pole,           &
                                     target_null_island,          &
                                     rotation_target,             &
                                     ROTATION_TARGET_NULL_ISLAND, &
                                     ROTATION_TARGET_NORTH_POLE
  use coord_transform_mod,     only: rebase_longitude_range

  implicit none

  integer(i_def) :: communicator = -999
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  type(reference_cube_type) :: cube_element

  type(gen_planar_type),  allocatable :: mesh_gen(:)
  type(ugrid_2d_type),    allocatable :: ugrid_2d(:)
  class(ugrid_file_type), allocatable :: ugrid_file

  type(ugrid_2d_type) :: ugrid_2d_lbc

  integer(i_def) :: fsize
  integer(i_def) :: xproc, yproc
  integer(i_def) :: n_mesh_maps = 0
  integer(i_def) :: n_targets

  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: target_edge_cells_x(:)
  integer(i_def),     allocatable :: target_edge_cells_y(:)

  integer(i_def),     allocatable :: target_edge_cells_x_tmp(:)
  integer(i_def),     allocatable :: target_edge_cells_y_tmp(:)
  character(str_def), allocatable :: target_mesh_names_tmp(:)

  ! Partition variables
  procedure(partitioner_interface),  &
                          pointer     :: partitioner_ptr => null()
  type(global_mesh_type), target      :: global_mesh
  type(global_mesh_type), pointer     :: global_mesh_ptr => null()
  type(ugrid_mesh_data_type)          :: ugrid_mesh_data
  type(partition_type)                :: partition
  type(local_mesh_type)               :: local_mesh

  ! Switches
  logical(l_def) :: l_found = .false.
  logical(l_def) :: any_duplicate_names = .false.

  ! Temporary variables
  character(str_def), allocatable :: requested_mesh_maps(:)
  character(str_def) :: first_mesh
  character(str_def) :: second_mesh
  character(str_def) :: tmp_str
  character(str_def) :: check_mesh(2)
  integer(i_def)     :: first_mesh_edge_cells_x, first_mesh_edge_cells_y
  integer(i_def)     :: second_mesh_edge_cells_x,second_mesh_edge_cells_y
  real(kind=r_def)   :: set_north_pole(2)
  real(kind=r_def)   :: set_null_island(2)

  character(str_def) :: name
  logical(l_def)     :: lbc_generated
  type(gen_lbc_type) :: lbc_mesh_gen

  character(str_def)  :: lon_str
  character(str_def)  :: lat_str

  ! Counters
  integer(i_def)    :: i, j, k, l, n_voids
  integer(i_native) :: log_level

  !===================================================================
  ! 1.0 Set the logging level for the run, should really be able
  !     to set it from the command line as an option
  !===================================================================
  call log_set_level(LOG_LEVEL_INFO)


  !===================================================================
  ! 2.0 Start up
  !===================================================================
  cube_element = reference_cube_type()

  call initialise_comm(communicator)
  call store_comm(communicator)

  ! Initialise halo functionality
  call initialise_halo_comms( communicator )

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()
  call initialise_logging(local_rank, total_ranks, "planar")


  !===================================================================
  ! 3.0 Read in the control namelists from file
  !===================================================================
  call get_initial_filename( filename )
  call read_configuration( filename )
  deallocate( filename )

  ! The number of mesh maps in the namelist array is unbounded
  ! and so may contain unset/empty array elements. Remove
  ! these from the initial count of mesh-maps.
  n_voids = count(cmdi == mesh_maps) + count('' == mesh_maps)
  if ( n_voids == 0 ) then
    n_mesh_maps = size(mesh_maps)
  else
    n_mesh_maps = size(mesh_maps) - n_voids
  end if

  !===================================================================
  ! 4.0 Perform some error checks on the namelist inputs
  !===================================================================
  ! 4.1a Check the namelist file enumeration: geometry
  log_level = LOG_LEVEL_ERROR
  select case (geometry)

  case (geometry_spherical, geometry_planar)

  case (emdi)
    write(log_scratch_space,'(A)') &
        'Enumeration key for geometry has not been set.'
    call log_event(log_scratch_space, log_level)
  case default
    write(log_scratch_space,'(A)')                     &
        'Unrecognised enumeration key for geometry:'// &
        trim(key_from_geometry(geometry))
    call log_event(log_scratch_space, log_level)
  end select


  ! 4.1b Check the namelist file enumeration: topology
  log_level = LOG_LEVEL_ERROR
  select case (topology)

  case ( topology_periodic, &
         topology_channel,  &
         topology_non_periodic )

    select case (topology)

    case (topology_periodic)
      if (geometry == geometry_planar) then
        if (.not. (periodic_x .and. periodic_y)) then
          write(log_scratch_space,'(A)')                  &
              'A periodic planar regional mesh should '// &
              'have all boundaries as periodic.'
          call log_event(log_scratch_space, log_level)
        end if
      else
         write(log_scratch_space,'(A)')                 &
             'A periodic spherical regional mesh is '// &
             'unsupported by this generator.'
          call log_event(log_scratch_space, log_level)
      end if

    case (topology_channel)
      if ( periodic_x .eqv. periodic_y ) then
        write(log_scratch_space,'(A)')                    &
            'A channel regional mesh should only have '// &
            'a single periodic axis.'
        call log_event(log_scratch_space, log_level)
      end if

    case (topology_non_periodic)
      if ( periodic_x .or. periodic_y ) then
        write(log_scratch_space,'(A)')                   &
            'A non-periodic regional mesh should not '// &
            'have any periodic boundaries.'
        call log_event(log_scratch_space, log_level)
      end if

    end select

  case (emdi)
    write(log_scratch_space,'(A)') &
        'Enumeration key for topology has not been set.'
    call log_event(log_scratch_space, log_level)

  case default
    write(log_scratch_space,'(A)')                      &
        'Unrecognised enumeration key for topology: '// &
        key_from_topology(topology)
    call log_event(log_scratch_space, log_level)
  end select


  ! 4.1c Check the namelist file enumeration: coord_sys
  log_level = LOG_LEVEL_ERROR
  select case (coord_sys)

  case (coord_sys_ll, coord_sys_xyz)

  case (emdi)
    write(log_scratch_space,'(A)') &
        'Enumeration key for coord_sys has not been set.'
    call log_event(log_scratch_space, log_level)

  case default
    write(log_scratch_space,'(A)')                      &
        'Unrecognised enumeration key for coord_sys:'// &
        trim(key_from_coord_sys(coord_sys))
    call log_event(log_scratch_space, log_level)

  end select


  ! 4.2 Check the number of meshes requested.
  if (n_meshes < 1) then
    write(log_scratch_space,'(A,I0,A)') &
        'Invalid number of meshes requested, (',n_meshes,')'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.3 Check that there are enough entries of edge cells
  !     to match the number of meshes requested.
  if ( size(edge_cells_x) < n_meshes .or. &
       size(edge_cells_y) < n_meshes ) then
    write(log_scratch_space,'(A,I0,A)')                      &
        'Not enough data in edge_cells_x/edge_cells_y for ', &
        n_meshes,' meshe(s).'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.4 Check for missing data.
  if ( any(edge_cells_x == imdi) .or. &
       any(edge_cells_y == imdi) ) then
    write(log_scratch_space,'(A)') &
        'Missing data in namelist variable, edge_cells_x/edge_cells_y'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.5 Check that all meshes requested have unique names.
  any_duplicate_names = any_duplicates(mesh_names)
  if (any_duplicate_names)  then
    write(log_scratch_space,'(A)')       &
        'Duplicate mesh names found, '// &
        'all requested meshes must have unique names.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! 4.6 Check that all mesh map requests are unique.
  any_duplicate_names = any_duplicates(mesh_maps)
  if (any_duplicate_names)  then
    write(log_scratch_space,'(A)')          &
        'Duplicate mesh requests found, '// &
        'please remove duplicate requests.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! Perform a number of checks related to mesh map
  ! requests.
  if (n_mesh_maps > 0) then
    do i=1, n_mesh_maps

      tmp_str = mesh_maps(i)
      check_mesh(1) = tmp_str(:index(tmp_str,':')-1)
      check_mesh(2) = tmp_str(index(tmp_str,':')+1:)
      first_mesh    = check_mesh(1)
      second_mesh   = check_mesh(2)

      first_mesh_edge_cells_x  = imdi
      first_mesh_edge_cells_y  = imdi
      second_mesh_edge_cells_x = imdi
      second_mesh_edge_cells_y = imdi

      do j=1, n_meshes
        if (trim(mesh_names(j)) == trim(first_mesh)) then
          first_mesh_edge_cells_x = edge_cells_x(j)
          first_mesh_edge_cells_y = edge_cells_y(j)
        end if

        if (trim(mesh_names(j)) == trim(second_mesh)) then
          second_mesh_edge_cells_x = edge_cells_x(j)
          second_mesh_edge_cells_y = edge_cells_y(j)
        end if
      end do

      ! 4.7 Check that mesh names in the map request exist.
      do j=1, size(check_mesh)

        l_found = .false.
        do k=1, n_meshes
          if (trim(check_mesh(j)) == trim(mesh_names(k))) then
            l_found = .true.
          end if
        end do

        if ( .not. l_found ) then
          write(log_scratch_space,'(A)')      &
              'Mesh "'//trim(check_mesh(j))// &
              '" not configured for this file.'
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if
      end do

      ! 4.8 Check the map request is not mapping at mesh
      !     to itself.
      if (trim(first_mesh) == trim(second_mesh)) then
        write(log_scratch_space,'(A)')                &
            'Found identical adjacent mesh names "'// &
           trim(mesh_maps(i))//'", requested for mapping.'
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if

      ! 4.9 Check that the number of edge cells of the meshes
      !     are not the same.
      if ( (first_mesh_edge_cells_x == second_mesh_edge_cells_x ) .and. &
            first_mesh_edge_cells_y == second_mesh_edge_cells_y ) then
        write(log_scratch_space,'(A,I0,A)')                      &
            'Found identical adjacent mesh edge cells, (',       &
            first_mesh_edge_cells_x,',',first_mesh_edge_cells_y, &
            '), requested for mapping "'// &
            trim(first_mesh)//'"-"'//trim(second_mesh)//'".'
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if

    end do
  end if

  ! 4.10 Check the requested lbc parent lam exists
  if (create_lbc_mesh) then
    l_found = .false.
    do i=1, n_meshes
      if ( trim(mesh_names(i)) == trim(lbc_parent_mesh) ) then
        l_found=.true.
        exit
      end if
    end do
    if ( .not. l_found ) then
      write( log_scratch_space, '(A)')                  &
          'The parent mesh, '// trim(lbc_parent_mesh)// &
          ' specified for LBC mesh generation does not exist.'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
    end if
  end if


  !===================================================================
  ! 5.0 Create unique list of Requested Mesh maps.
  !     Each map request will create two maps, one in each direction
  !===================================================================
  if (n_mesh_maps > 0) then
    allocate(requested_mesh_maps(n_mesh_maps*2))
    j=1
    do i=1, n_mesh_maps
      tmp_str = mesh_maps(i)
      first_mesh  = tmp_str(:index(tmp_str,':')-1)
      second_mesh = tmp_str(index(tmp_str,':')+1:)
      write(requested_mesh_maps(j),   '(A)') &
          trim(first_mesh)//':'//trim(second_mesh)
      write(requested_mesh_maps(j+1), '(A)') &
          trim(second_mesh)//':'//trim(first_mesh)
      j=j+2
    end do
  end if


  !===================================================================
  ! 6.0 Report/Check what the code thinks is requested by user
  !===================================================================
  log_level=LOG_LEVEL_INFO

  write(log_scratch_space, '(A)')    &
      '===================================================================='
  call log_event( log_scratch_space, log_level )

  write(log_scratch_space, '(A)')    &
      'Mesh geometry: ' // trim(key_from_geometry(geometry))
  call log_event( log_scratch_space, log_level )

  write(log_scratch_space,'(A)')     &
      'Mesh topology: ' // trim(key_from_topology(topology))
  call log_event( log_scratch_space, log_level )

  write(log_scratch_space, '(A)')    &
      'Co-ordinate system: '// trim(key_from_coord_sys(coord_sys))
  call log_event( log_scratch_space, log_level )

  write(log_scratch_space, '(A,L1)') &
      'Periodic in x-axis: ', periodic_x
  call log_event( log_scratch_space, log_level )

  write(log_scratch_space, '(A,L1)') &
      'Periodic in y-axis: ', periodic_y
  call log_event( log_scratch_space, log_level )


  write(log_scratch_space, '(A)')    &
      '===================================================================='
  call log_event( log_scratch_space, log_level )
  call log_event( "Generating mesh(es):", log_level )

  ! 6.1 Generate objects which know how to generate each requested
  !     unique mesh.
  allocate( mesh_gen (n_meshes) )
  allocate( ugrid_2d (n_meshes) )

  ! 6.2 Assign temporary arrays for target meshes in requested maps
  if (n_mesh_maps > 0) then
    if (allocated(target_mesh_names_tmp))   deallocate(target_mesh_names_tmp)
    if (allocated(target_edge_cells_x_tmp)) deallocate(target_edge_cells_x_tmp)
    if (allocated(target_edge_cells_y_tmp)) deallocate(target_edge_cells_y_tmp)
    allocate( target_mesh_names_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_x_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_y_tmp(n_mesh_maps*2) )
  end if

  if ( geometry == geometry_spherical .and. &
       coord_sys == coord_sys_ll ) then
    if (rotate_mesh) then

      write(log_scratch_space, '(A)') &
         '  Rotation of mesh requested with: '
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

      select case ( rotation_target )
      case ( ROTATION_TARGET_NULL_ISLAND )
        ! Use the domain_centre (Null Island) rather than pole as input
        set_north_pole(:) = get_target_north_pole(target_null_island)
        set_null_island(:) = target_null_island
        write(log_scratch_space,'(A)')       &
           '    Target pole will be derived from Null Island.'
        call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

        write(lon_str,'(F10.2)') set_null_island(1)
        write(lat_str,'(F10.2)') set_null_island(2)
        write(log_scratch_space,'(A)')        &
           '    Null Island [lon,lat]: ['  // &
           trim(adjustl(lon_str)) // ',' //   &
           trim(adjustl(lat_str)) // ']'
        call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
      case ( ROTATION_TARGET_NORTH_POLE )
        set_north_pole(:)  = target_north_pole(:)
        set_null_island(:) = get_target_null_island(target_north_pole)
      end select

      ! Ensure the requested target longitudes are in the range -180,180
      set_null_island(1) = rebase_longitude_range( set_null_island(1), -180.0_r_def)
      set_north_pole(1)  = rebase_longitude_range( set_north_pole(1), -180.0_r_def)

      write(lon_str,'(F10.2)') set_north_pole(1)
      write(lat_str,'(F10.2)') set_north_pole(2)
      write(log_scratch_space,'(A)')        &
         '    Target pole [lon,lat]: ['  // &
         trim(adjustl(lon_str)) // ',' //   &
         trim(adjustl(lat_str)) // ']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

      write(lon_str,'(F10.2)') first_node(1)
      write(lat_str,'(F10.2)') first_node(2)
      write(log_scratch_space,'(A)')        &
         '    First node  [lon,lat]: ['  // &
         trim(adjustl(lon_str)) // ',' //   &
         trim(adjustl(lat_str)) // ']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    end if
  end if

  do i=1, n_meshes

    write(log_scratch_space,'(A,2(I0,A))')             &
       '  Creating Mesh: '// trim(mesh_names(i))//'(', &
                            edge_cells_x(i), ',',      &
                            edge_cells_y(i), ')'
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)


    ! 6.3 Get any target mappings requested for this mesh
    n_targets = 0
    if (n_mesh_maps > 0) then

      target_mesh_names_tmp = cmdi
      target_edge_cells_x_tmp = imdi
      target_edge_cells_y_tmp = imdi
      l=1
      do j=1, size(requested_mesh_maps)
        tmp_str= requested_mesh_maps(j)
        if (tmp_str( :index(tmp_str,':')-1) == trim(mesh_names(i))) then
          do k=1, n_meshes
            if ( trim(tmp_str( index(tmp_str,':')+1:)) ==  &
                 trim(mesh_names(k)) ) then
              target_mesh_names_tmp(l)   = trim(mesh_names(k))
              target_edge_cells_x_tmp(l) = edge_cells_x(k)
              target_edge_cells_y_tmp(l) = edge_cells_y(k)
              l=l+1
            end if
          end do
        end if
      end do
      n_targets=l-1
    end if  ! n_mesh_maps > 0

    ! 6.4 Call generation strategy
    if (n_targets == 0 .or. n_meshes == 1 ) then

      mesh_gen(i) = gen_planar_type(                          &
                        reference_element  = cube_element,    &
                        mesh_name          = mesh_names(i),   &
                        geometry           = geometry,        &
                        topology           = topology,        &
                        coord_sys          = coord_sys,       &
                        edge_cells_x       = edge_cells_x(i), &
                        edge_cells_y       = edge_cells_y(i), &
                        periodic_x         = periodic_x,      &
                        periodic_y         = periodic_y,      &
                        domain_x           = domain_x,        &
                        domain_y           = domain_y,        &
                        rotate_mesh        = rotate_mesh,     &
                        target_north_pole  = set_north_pole,  &
                        target_null_island = set_null_island, &
                        first_node         = first_node )

    else if (n_meshes > 1) then

      if (allocated(target_mesh_names))   deallocate(target_mesh_names)
      if (allocated(target_edge_cells_x)) deallocate(target_edge_cells_x)
      if (allocated(target_edge_cells_y)) deallocate(target_edge_cells_y)

      allocate( target_mesh_names(n_targets)   )
      allocate( target_edge_cells_x(n_targets) )
      allocate( target_edge_cells_y(n_targets) )
      target_mesh_names(:)   = target_mesh_names_tmp(:n_targets)
      target_edge_cells_x(:) = target_edge_cells_x_tmp(:n_targets)
      target_edge_cells_y(:) = target_edge_cells_y_tmp(:n_targets)

      write(log_scratch_space,'(A,I0)') '    Maps to:'
      do j=1, n_targets
        write(log_scratch_space,'(2(A,I0),A)') &
            trim(log_scratch_space)//' '//     &
            trim(target_mesh_names(j))//       &
            '(',target_edge_cells_x(j),',',    &
            target_edge_cells_y(j),')'
      end do

      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)

      mesh_gen(i) = gen_planar_type(                               &
                        cube_element, mesh_names(i),               &
                        geometry, topology, coord_sys,             &
                        edge_cells_x(i), edge_cells_y(i),          &
                        periodic_x, periodic_y,                    &
                        domain_x, domain_y,                        &
                        target_mesh_names   = target_mesh_names,   &
                        target_edge_cells_x = target_edge_cells_x, &
                        target_edge_cells_y = target_edge_cells_y, &
                        rotate_mesh         = rotate_mesh,         &
                        target_north_pole   = set_north_pole,      &
                        target_null_island  = set_null_island,     &
                        first_node          = first_node )

    else
      write(log_scratch_space, "(A,I0,A)") &
         '  Number of meshes is negative [', n_meshes,']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR)
    end if

    ! Pass the generation object to the ugrid file writer
    call ugrid_2d(i)%set_by_generator(mesh_gen(i))

    if (allocated(target_mesh_names))   deallocate(target_mesh_names)
    if (allocated(target_edge_cells_x)) deallocate(target_edge_cells_x)
    if (allocated(target_edge_cells_y)) deallocate(target_edge_cells_y)

  end do

  call log_event( "...generation complete.", LOG_LEVEL_INFO )
  write(log_scratch_space, '(A)') &
      '===================================================================='
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  !=================================================================
  ! 7.0 Partitioning
  !=================================================================
  if (n_partitions > 0 ) then

    ! 7.1 Set partitioning parameters.
    call set_partition_parameters( xproc, yproc, &
                                   partitioner_ptr )

    do i=1, n_meshes

      ! 7.2 Create global mesh object.
      call ugrid_mesh_data%set_by_ugrid_2d( ugrid_2d(i) )
      global_mesh = global_mesh_type( ugrid_mesh_data, NPANELS )
      call ugrid_mesh_data%clear()
      global_mesh_ptr => global_mesh

      ! 7.3 Create global mesh partitions
      do j=0,n_partitions-1
        partition = partition_type( global_mesh_ptr,   &
                                    partitioner_ptr,   &
                                    xproc, yproc,      &
                                    max_stencil_depth, &
                                    j, n_partitions )
        call local_mesh%initialise( global_mesh_ptr, partition )
      end do

      global_mesh_ptr => null()

    end do
  end if

  !===================================================================
  ! 8.0 Now the write out mesh to the NetCDF file
  !===================================================================
  do i=1, n_meshes

    if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

    call ugrid_2d(i)%set_file_handler(ugrid_file)

    if (i==1) then
      call ugrid_2d(i)%write_to_file( trim(mesh_filename) )
    else
      call ugrid_2d(i)%append_to_file( trim(mesh_filename) )
    end if

    inquire(file=mesh_filename, size=fsize)
    write( log_scratch_space, '(A,I0,A)')                 &
        'Adding mesh (' // trim(mesh_names(i)) //         &
        ') to ' // trim(adjustl(mesh_filename)) // ' - ', &
        fsize, ' bytes written.'

    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if (allocated(ugrid_file)) deallocate(ugrid_file)

  end do ! n_meshes


  !===================================================================
  ! 9.0 Now create/output LBC mesh
  !===================================================================
  ! A LBC mesh is created from a parent planar mesh strategy that has
  ! been generated. The name of the resulting LBC mesh will be:
  !
  !    <parent mesh name>-lbc
  !
  if (create_lbc_mesh) then

    lbc_generated = .false.

    do i=1, size(mesh_gen)

      if (lbc_generated) exit

      call mesh_gen(i)%get_metadata(mesh_name=name)
      if (trim(name) == trim(lbc_parent_mesh)) then
        lbc_mesh_gen = gen_lbc_type(mesh_gen(i), lbc_rim_depth)

        call ugrid_2d_lbc%set_by_generator(lbc_mesh_gen)
        if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

        call ugrid_2d_lbc%set_file_handler(ugrid_file)
        call ugrid_2d_lbc%append_to_file( trim(mesh_filename) )

        inquire(file=mesh_filename, size=fsize)
        write( log_scratch_space, '(A,I0,A)')                &
            'Adding lbc mesh for ' // trim(mesh_names(i)) // &
            ' to ' // trim(adjustl(mesh_filename)) // ' - ', &
            fsize, ' bytes written.'

        call log_event( log_scratch_space, LOG_LEVEL_INFO )
        if (allocated(ugrid_file)) deallocate(ugrid_file)

        lbc_generated = .true.

      end if
    end do
  end if

  !===================================================================
  ! 9.0 Clean up and Finalise
  !===================================================================

  if ( allocated( mesh_gen ) ) deallocate (mesh_gen)

  if ( allocated( requested_mesh_maps     ) ) deallocate (requested_mesh_maps)
  if ( allocated( target_mesh_names       ) ) deallocate (target_mesh_names)
  if ( allocated( target_edge_cells_x     ) ) deallocate (target_edge_cells_x)
  if ( allocated( target_edge_cells_y     ) ) deallocate (target_edge_cells_y)
  if ( allocated( target_mesh_names_tmp   ) ) deallocate (target_mesh_names_tmp)
  if ( allocated( target_edge_cells_x_tmp ) ) deallocate (target_edge_cells_x_tmp)
  if ( allocated( target_edge_cells_y_tmp ) ) deallocate (target_edge_cells_y_tmp)

  call finalise_halo_comms()

  call finalise_comm()

  call finalise_logging()

  call final_configuration()

end program planar_mesh_generator
