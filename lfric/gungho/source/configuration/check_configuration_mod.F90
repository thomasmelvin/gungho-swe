!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module check_configuration_mod

  use constants_mod,        only: i_def, l_def
  use mixing_config_mod,    only: smagorinsky,            &
                                  viscosity,              &
                                  viscosity_mu
  use subgrid_config_mod,   only: dep_pt_stencil_extent,  &
                                  rho_approximation_stencil_extent
  use transport_config_mod, only: operators,                   &
                                  operators_fv,                &
                                  consistent_metric,           &
                                  fv_horizontal_order,         &
                                  fv_vertical_order,           &
                                  profile_size,                &
                                  scheme,                      &
                                  horizontal_method,           &
                                  vertical_method,             &
                                  max_vert_cfl_calc,           &
                                  max_vert_cfl_calc_dep_point
  use transport_enumerated_types_mod,                     &
                            only: scheme_mol_3d,          &
                                  scheme_ffsl_3d,         &
                                  scheme_split,           &
                                  split_method_mol,       &
                                  split_method_ffsl

  implicit none

  private

  public :: check_configuration
  public :: check_any_scheme_mol
  public :: check_any_scheme_split
  public :: check_any_scheme_ffsl
  public :: check_horz_dep_pts
  public :: check_vert_dep_pts
  public :: get_required_stencil_depth

contains

  !> @brief Check the namelist configuration for unsupported combinations
  !>        of options and flag up errors and warnings
  subroutine check_configuration()
    use log_mod,                     only: log_event,                          &
                                           log_scratch_space,                  &
                                           LOG_LEVEL_ERROR,                    &
                                           LOG_LEVEL_WARNING,                  &
                                           LOG_LEVEL_INFO
    use constants_mod,               only: EPS,                                &
                                           r_def,                              &
                                           i_def
    use finite_element_config_mod,   only: cellshape,                          &
                                           cellshape_triangle,                 &
                                           element_order,                      &
                                           rehabilitate,                       &
                                           coord_order,                        &
                                           coord_system,                       &
                                           coord_system_alphabetaz,            &
                                           coord_system_lonlatz
    use formulation_config_mod,      only: use_physics,                        &
                                           use_wavedynamics,                   &
                                           dlayer_on
    use io_config_mod,               only: write_diag,                         &
                                           use_xios_io
    use planet_config_mod,           only: gravity,                            &
                                           radius,                             &
                                           omega,                              &
                                           rd,                                 &
                                           cp,                                 &
                                           p_zero,                             &
                                           scaling_factor
    use timestepping_config_mod,     only: method,                             &
                                           method_semi_implicit,               &
                                           dt,                                 &
                                           alpha,                              &
                                           outer_iterations,                   &
                                           inner_iterations
    use base_mesh_config_mod,        only: geometry,                           &
                                           geometry_spherical,                 &
                                           geometry_planar,                    &
                                           topology,                           &
                                           topology_fully_periodic,            &
                                           topology_non_periodic
    use damping_layer_config_mod,    only: dl_base,                            &
                                           dl_str
    use extrusion_config_mod,        only: domain_top
    use mixed_solver_config_mod,     only: reference_reset_freq
    use helmholtz_solver_config_mod, only:                                     &
                            helmholtz_solver_preconditioner => preconditioner, &
                            preconditioner_tridiagonal
    implicit none

      logical(kind=l_def) :: any_scheme_mol

      call log_event( 'Checking gungho configuration...', LOG_LEVEL_INFO )

      ! Check the options in the finite element namelist
      if ( cellshape == cellshape_triangle ) then
        write( log_scratch_space, '(A)' ) 'Triangular elements are unsupported'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( element_order < 0 ) then
        write( log_scratch_space, '(A,I4,A)' ) 'Invalid choice: element order ', &
          element_order, ' must be non-negative'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( .not. rehabilitate ) then
        write( log_scratch_space, '(A)' ) 'Only rehabilitated W3 function space is allowed'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_order < 0 ) then
        write( log_scratch_space, '(A,I4,A)' ) 'Invalid choice: coordinate order ', &
          coord_order, ' must be non-negative'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_order == 0 .and. &
           ( topology /= topology_non_periodic .or. geometry == geometry_planar ) ) then
        write( log_scratch_space, '(A)' ) 'For planar geometry or periodic meshes, coordinate order must be positive'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_system == coord_system_alphabetaz .and. geometry /= geometry_spherical ) then
        write( log_scratch_space, '(A)' ) '(alpha,beta) coordinate system is only valid with spherical geometry'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_system == coord_system_alphabetaz .and. topology /= topology_fully_periodic ) then
        ! This could change in future if we were to add meshes that were a single panel
        write( log_scratch_space, '(A)' ) '(alpha,beta) coordinate system is only valid with fully-periodic topology'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_system == coord_system_lonlatz .and. geometry /= geometry_spherical ) then
        write( log_scratch_space, '(A)' ) '(longitude,latitude) coordinate system is only valid with spherical geometry'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_system == coord_system_lonlatz .and. topology == topology_fully_periodic ) then
        write( log_scratch_space, '(A)' ) '(longitude,latitude) coordinate system is not valid with fully-periodic topology'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if


      ! Check the options in the formulation namelist
      if ( .not. use_physics .and. .not. use_wavedynamics ) then
        write( log_scratch_space, '(A)' ) 'Wave dynamics and physics turned off'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the io namelist
      if ( .not. write_diag ) then
        write( log_scratch_space, '(A)' ) 'Diagnostic output not enabled'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the planet namelist
      if ( gravity < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative gravity: ', &
          gravity
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( radius < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative radius: ', &
          radius
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( omega < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative omega: ', &
          omega
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( rd < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative Rd: ', &
          rd
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( cp < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative cp: ', &
          cp
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( p_zero < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative p_zero: ', &
          p_zero
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( scaling_factor < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Invalid choice: Zero or negative scaling factor: ', &
          scaling_factor
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Check the timestepping namelist
      if ( dt < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative dt: ', &
          dt
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( method == method_semi_implicit ) then
        if( alpha < 0.5_r_def ) then
          write( log_scratch_space, '(A,E16.8)' ) 'alpha < 1/2 likely to be unstable: ',&
            alpha
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
        if ( outer_iterations < 1 ) then
          write( log_scratch_space, '(A,I4)' ) 'Invalid Choice: outer_iterations must be at least 1:', &
          outer_iterations
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( inner_iterations < 1 ) then
          write( log_scratch_space, '(A,I4)' ) 'Invalid Choice: inner_iterations must be at least 1:', &
          inner_iterations
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if

      ! Check the transport namelist
      if ( geometry == geometry_spherical .and.  consistent_metric) then
        write( log_scratch_space, '(A)' ) 'Consistent metric option only valid for planar geometries'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      any_scheme_mol = check_any_scheme_mol()
      if (any_scheme_mol) then
        ! Check that flux orders are even
        if ( mod(fv_horizontal_order,2_i_def) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'fv_horizontal_order must be even'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( mod(fv_vertical_order,2_i_def) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'fv_vertical_order must be even'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if

      ! Check the mixing namelist
      if ( viscosity .and. geometry == geometry_spherical ) then
        write( log_scratch_space, '(A)' ) 'Viscosity might not work in spherical domains'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( viscosity .and. viscosity_mu < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Negative viscosity coefficient: Anti-diffusion: ', &
          viscosity_mu
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the damping layer namelist
      if ( dlayer_on .and. (dl_base < 0.0 .or. dl_base > domain_top) ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Damping layer base lies outside fo domain: ',&
          dl_base
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( dlayer_on .and. dl_str < 0.0 ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Damping layer strength is negative: ',&
          dl_str
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check for options that are invalid with higher order elements
      if ( element_order > 0 ) then
        if ( operators == operators_fv ) then
          write( log_scratch_space, '(A)' ) 'FV transport operators only valid for element_order = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( viscosity ) then
          write( log_scratch_space, '(A)' ) 'Viscosity only valid for element_order = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( helmholtz_solver_preconditioner == preconditioner_tridiagonal ) then
          write( log_scratch_space, '(A)' ) 'Tridiagonal helmholtz preconditioner only valid for  element_order = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( use_xios_io ) then
          write( log_scratch_space, '(A)' ) 'xios output may not work with element order > 0'
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
      end if

      if ( method == method_semi_implicit ) then
        ! Check the mixed solver namelist
        if ( reference_reset_freq > outer_iterations*inner_iterations ) then
          write( log_scratch_space, '(A)' ) 'reference_reset_freq greater than total number of si iterations'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( mod(outer_iterations*inner_iterations,reference_reset_freq) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'reference_reset_freq not a divsor of total number of si iterations'
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
      end if

      call log_event( '...Check gungho config done', LOG_LEVEL_INFO )

  end subroutine check_configuration


  !> @brief   Determine required stencil depth for the current configuration.
  !> @details Depending on the choice of science schemes the required local
  !>          mesh needs to support the anticipated stencils. This function
  !>          returns required stencil depth that needs to be supported.
  !> @return  stencil_depth
  !> @todo    stencil_depth will eventually become either a configuration
  !>          item, or will be autogenerated from kernel metadata.
  !>
  !===========================================================================
  function get_required_stencil_depth() result(stencil_depth)

    implicit none

    integer(kind=i_def) :: stencil_depth
    logical(kind=l_def) :: any_scheme_ffsl

    stencil_depth = 1

    ! Smagorinsky (or boundary layers) appears to need larger haloes
    if (smagorinsky) stencil_depth = max( stencil_depth, 2 )

    if (operators == operators_fv) then
      ! Need larger haloes for fv operators
      stencil_depth  = max( stencil_depth, fv_horizontal_order/2 )
    end if

    any_scheme_ffsl = check_any_scheme_ffsl()

    if ( any_scheme_ffsl ) then
        stencil_depth = max( stencil_depth,          &
                             dep_pt_stencil_extent + &
                             rho_approximation_stencil_extent )
    end if

  end function get_required_stencil_depth

  !> @brief   Determine whether any of the transport schemes are MoL
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using the Method of Lines
  !>          scheme
  !> @return  any_scheme_mol
  function check_any_scheme_mol() result(any_scheme_mol)

    implicit none

    logical(kind=l_def) :: any_scheme_mol
    integer(kind=i_def) :: i

    any_scheme_mol = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_mol_3d ) .or.                      &
           ( scheme(i) == scheme_split .and.                        &
             ( vertical_method(i) == split_method_mol .or.          &
               horizontal_method(i) == split_method_mol ) ) ) then
        any_scheme_mol = .true.
        exit
      end if
    end do

  end function check_any_scheme_mol

  !> @brief   Determine whether any of the transport schemes are split
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using the split vertical-
  !>          horizontal scheme
  !> @return  any_scheme_split
  function check_any_scheme_split() result(any_scheme_split)

    implicit none

    logical(kind=l_def) :: any_scheme_split
    integer(kind=i_def) :: i

    any_scheme_split = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split ) then
        any_scheme_split = .true.
        exit
      end if
    end do

  end function check_any_scheme_split

  !> @brief   Determine whether any of the transport schemes are FFSL
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using the Flux-Form
  !>          Semi-Lagrangian scheme
  !> @return  any_scheme_ffsl
  function check_any_scheme_ffsl() result(any_scheme_ffsl)

    implicit none

    logical(kind=l_def) :: any_scheme_ffsl
    integer(kind=i_def) :: i

    any_scheme_ffsl = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_ffsl_3d ) .or.                      &
           ( scheme(i) == scheme_split .and.                        &
             ( vertical_method(i) == split_method_ffsl .or.          &
               horizontal_method(i) == split_method_ffsl ) ) ) then
        any_scheme_ffsl = .true.
        exit
      end if
    end do

  end function check_any_scheme_ffsl

  !> @brief   Determine whether horizontal departure points need computing
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using a scheme that
  !>          requires horizontal departure points to be computed
  !> @return  any_horz_dep_pts
  function check_horz_dep_pts() result(any_horz_dep_pts)

    implicit none

    logical(kind=l_def) :: any_horz_dep_pts
    integer(kind=i_def) :: i

    any_horz_dep_pts = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_ffsl_3d ) .or.                     &
           ( scheme(i) == scheme_split .and.                        &
             horizontal_method(i) == split_method_ffsl ) ) then
        any_horz_dep_pts = .true.
        exit
      end if
    end do

  end function check_horz_dep_pts

  !> @brief   Determine whether vertical departure points need computing
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using a scheme that
  !>          requires vertical departure points to be computed
  !> @return  any_vert_dep_pts
  function check_vert_dep_pts() result(any_vert_dep_pts)

    implicit none

    logical(kind=l_def) :: any_vert_dep_pts
    integer(kind=i_def) :: i

    any_vert_dep_pts = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_ffsl_3d ) .or.                     &
           ( scheme(i) == scheme_split .and.                        &
             vertical_method(i) /= split_method_mol ) .or.          &
           ( max_vert_cfl_calc == max_vert_cfl_calc_dep_point ) ) then
        any_vert_dep_pts = .true.
        exit
      end if
    end do

  end function check_vert_dep_pts

end module check_configuration_mod
