!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Module for computing and outputting derived diagnostics
!!
!!  @details Computes various derived diagnostics that are written out
!!           by the diagnostic system. Also computes norms of some fields
!!           which are written to the output log.
!-------------------------------------------------------------------------------
module diagnostics_calc_mod

  use clock_mod,                     only: clock_type
  use constants_mod,                 only: i_def, r_def, str_max_filename
  use diagnostic_alg_mod,            only: divergence_diagnostic_alg,   &
                                           hydbal_diagnostic_alg,       &
                                           vorticity_diagnostic_alg
  use io_config_mod,                 only: use_xios_io,          &
                                           nodal_output_on_w3
  use files_config_mod,              only: diag_stem_name
  use project_output_mod,            only: project_output
  use io_mod,                        only: ts_fname, &
                                           nodal_write_field
  use lfric_xios_write_mod,          only: write_field_face, &
                                           write_field_edge
  use diagnostics_io_mod,            only: write_scalar_diagnostic,     &
                                           write_vector_diagnostic
  use field_mod,                     only: field_type
  use field_parent_mod,              only: write_interface
  use fs_continuity_mod,             only: W3
  use moist_dyn_mod,                 only: num_moist_factors
  use log_mod,                       only: log_event,         &
                                           log_set_level,     &
                                           log_scratch_space, &
                                           LOG_LEVEL_ERROR,   &
                                           LOG_LEVEL_INFO,    &
                                           LOG_LEVEL_DEBUG,   &
                                           LOG_LEVEL_TRACE
  use mesh_mod,                      only: mesh_type

  implicit none
  private
  public :: write_divergence_diagnostic, &
            write_hydbal_diagnostic,     &
            write_vorticity_diagnostic

contains

!-------------------------------------------------------------------------------
!>  @brief    Handles divergence diagnostic processing
!!
!!  @details  Handles divergence diagnostic processing
!!
!!> @param[in] u_field     The u field
!!> @param[in] ts          Timestep
!!> @param[in] mesh        Mesh
!-------------------------------------------------------------------------------

subroutine write_divergence_diagnostic(u_field, clock, mesh)
  implicit none

  type(field_type),  intent(in)    :: u_field
  class(clock_type), intent(in)    :: clock
  type(mesh_type),   intent(in), pointer :: mesh

  type(field_type)                :: div_field
  real(r_def)                     :: l2_norm


  procedure(write_interface), pointer  :: tmp_write_ptr

  ! Create the divergence diagnostic
  call divergence_diagnostic_alg( div_field, l2_norm, u_field, mesh )

  write( log_scratch_space, '(A,E16.8)' )  &
       'L2 of divergence =',l2_norm
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  if (use_xios_io) then
      !If using XIOS, we need to set a field I/O method appropriately
      tmp_write_ptr => write_field_face
      call div_field%set_write_behaviour(tmp_write_ptr)
  end if

  call write_scalar_diagnostic( 'divergence', div_field, &
                                clock, mesh, .false. )

  nullify(tmp_write_ptr)

end subroutine write_divergence_diagnostic

!-------------------------------------------------------------------------------
!>  @brief    Handles hydrostatic balance diagnostic processing
!!
!!  @details  Handles hydrostatic balance diagnostic processing
!!
!!> @param[in] theta_field   The theta field
!!> @param[in] exner_field   The exner field
!!> @param[in] mesh          Mesh
!-------------------------------------------------------------------------------

subroutine write_hydbal_diagnostic( theta_field, moist_dyn_field, exner_field,  &
                                    mesh )

  implicit none

  type(field_type), intent(in)    :: theta_field
  type(field_type), intent(in)    :: moist_dyn_field(num_moist_factors)
  type(field_type), intent(in)    :: exner_field
  type(mesh_type),  intent(in), pointer :: mesh

  real(r_def)                     :: l2_norm = 0.0_r_def

  call hydbal_diagnostic_alg(l2_norm, theta_field, moist_dyn_field,            &
                             exner_field, mesh)

  write( log_scratch_space, '(A,E16.8)' )  &
       'L2 of hydrostatic imbalance =', l2_norm
  call log_event( log_scratch_space, LOG_LEVEL_INFO )


end subroutine write_hydbal_diagnostic


!-------------------------------------------------------------------------------
!>  @brief    Handles vorticity diagnostic processing
!!
!!  @details  Handles vorticity diagnostic processing
!!
!!> @param[in] u_field   The wind field
!!> @param[in] timestep  Model timestep to index the output file
!-------------------------------------------------------------------------------

subroutine write_vorticity_diagnostic(u_field, clock)
  implicit none

  type(field_type),  intent(in) :: u_field
  class(clock_type), intent(in) :: clock

  type(field_type) :: vorticity

  call vorticity_diagnostic_alg(vorticity, u_field)

  call write_vector_diagnostic('xi', vorticity, clock, &
                               vorticity%get_mesh(), nodal_output_on_w3)

end subroutine write_vorticity_diagnostic

end module diagnostics_calc_mod
