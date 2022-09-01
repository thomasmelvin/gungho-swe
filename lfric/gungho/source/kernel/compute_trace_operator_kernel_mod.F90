!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes boundary integral part of rhs of the momentum equation for
!>        the nonlinear equations.
!>
!> @details
!> The kernel computes the boundary integral arising in the hybridized
!> equations. This equation enforces the continuity of the "broken" wind
!> variable across inter-elemental boundaries. The transpose of this
!> operator appears in the hybridized momentum equations, where the
!> traces approximate the pressure on mesh faces.
!>
!> This consists of: gamma*dot(u, normal_vector)
!>
!> where gamma is a test function in the trace space, and u is
!> a trial function in the broken W2 space. This operator maps
!> from the broken W2 space to the trace space.
!>
module compute_trace_operator_kernel_mod

  use argument_mod,              only: arg_type, func_type,             &
                                       reference_element_data_type,     &
                                       GH_OPERATOR, GH_FIELD,           &
                                       GH_REAL, GH_WRITE, GH_READ,      &
                                       GH_BASIS, GH_DIFF_BASIS,         &
                                       CELL_COLUMN,                     &
                                       GH_QUADRATURE_face,              &
                                       ANY_SPACE_9,                     &
                                       outward_normals_to_faces

  use reference_element_mod,     only: FACE_OFFSET
  use constants_mod,             only: r_def, i_def
  use fs_continuity_mod,         only: W2broken, W2trace
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer

  type, public, extends(kernel_type) :: compute_trace_operator_type
    private
    type(arg_type) :: meta_args(2) = (/                                 &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2trace, W2broken),   &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_9)          &
         /)
    type(func_type) :: meta_funcs(3) = (/                               &
         func_type(W2broken,    GH_BASIS),                              &
         func_type(W2trace,     GH_BASIS),                              &
         func_type(ANY_SPACE_9, GH_DIFF_BASIS)                          &
         /)
    type(reference_element_data_type) :: meta_reference_element(1) = (/ &
         reference_element_data_type( outward_normals_to_faces )        &
         /)
     integer :: gh_shape = GH_QUADRATURE_face
     integer :: operates_on = CELL_COLUMN
   contains
     procedure, nopass :: compute_trace_operator_code
   end type compute_trace_operator_type

   !-------------------------------------------------------------------------------
   ! Contained functions/subroutines
   !-------------------------------------------------------------------------------

   public :: compute_trace_operator_code

contains

  !> @brief Compute the boundary integral terms in the hybridized formulation of
  !!        the momentum and transmission equations (terms with Lagrange multipliers).
  !! @param[in] cell The cell id
  !! @param[in] nlayers Number of layers
  !! @param[in] ncell_3d ncell*nlayers
  !! @param[in,out] trace_op The operator coupling broken W2 and W2 trace functions
  !! @param[in] ndf_w2b Number of degrees of freedom per cell for W2broken space
  !! @param[in] ndf_w2t Number of degrees of freedom per cell for W2trace space
  !! @param[in] nqp Number of quadrature points on each face
  !! @param[in] wqp Quadrature weights for quadrature points on faces
  !! @param[in] w2b_basis Basis functions in W2broken space evaluated at gaussian
  !!                      quadrature points on horizontal and vertical faces
  !! @param[in] w2t_basis Basis functions in W2trace space evaluated at gaussian
  !!                      quadrature points on horizontal and vertical faces
  !! @param[in] face_entity_map Array mapping dof index to face entity
  !! @param[in] nfaces_re The number of faces (3D) or edges (2D) in each cell
  !! @param[in] outward_normals_to_faces Vector of normals to the reference
  !!                                     element "outward faces"
  !!
  subroutine compute_trace_operator_code( cell, nlayers, ncell_3d,   &
                                          trace_op,                  &
                                          ndf_w2b,                   &
                                          ndf_w2t,                   &
                                          nqp, wqp,                  &
                                          w2b_basis,                 &
                                          w2t_basis,                 &
                                          face_entity_map,           &
                                          nfaces_re,                 &
                                          outward_normals_to_faces )

    implicit none

    ! Argument declarations
    integer(kind=i_def),                                     intent(in) :: cell, nlayers, ncell_3d
    integer(kind=i_def),                                     intent(in) :: ndf_w2b, ndf_w2t
    integer(kind=i_def),                                     intent(in) :: nfaces_re, nqp

    real(kind=r_def), dimension(ndf_w2t, ndf_w2b, ncell_3d), intent(inout) :: trace_op
    real(kind=r_def), dimension(3, ndf_w2b, nqp, nfaces_re), intent(in)    :: w2b_basis
    real(kind=r_def), dimension(1, ndf_w2t, nqp, nfaces_re), intent(in)    :: w2t_basis

    integer(kind=i_def), dimension(ndf_w2t), intent(in) :: face_entity_map
    real(kind=r_def),                        intent(in) :: outward_normals_to_faces(:, :)

    ! Internal variables
    integer(kind=i_def) :: df, dfb, dft, k, ik, face, qp

    real(kind=r_def), dimension(nqp, nfaces_re), intent(in) :: wqp
    real(kind=r_def)                                        :: integrand

    ! Loop over layers
    do k = 0, nlayers - 1
      ik = k + 1 + (cell - 1) * nlayers

      do face = 1, nfaces_re
        do dfb = 1, ndf_w2b
          do dft = 1, ndf_w2t
            ! If trace dofs lie on the face, compute the integrand. Otherwise,
            ! don't accumulate.
            !
            ! The `face_entity_map` provides a unique label for each face of
            ! the reference element. Each face is labeled by the reference element
            ! using `face_idx + FACE_OFFSET`, where the `face_idx` is an enumerator
            ! for the number of each face. The ordering of faces (in order from 1 to 6)
            ! is: W, S, E, N, B, and T.
            if ( face_entity_map(dft) - FACE_OFFSET == face ) then
              trace_op(dft, dfb, ik) = 0.0_r_def
              do qp = 1, nqp
                ! NOTE: The surface integral of a Piola-mapped vector field (like
                ! a W2 or W2broken field) over the boundary of a physical cell
                ! is equivalent to surface integral over the reference cell. Hence,
                ! no need for Jacobians.
                integrand = wqp(qp, face) * w2t_basis(1, dft, qp, face)  &
                          * dot_product(w2b_basis(:, dfb, qp, face),     &
                                        outward_normals_to_faces(:, face))
                trace_op(dft, dfb, ik) = trace_op(dft, dfb, ik) + integrand
              end do
            end if ! End of trace dof / face check
          end do ! End of trace dof loop
        end do ! End of broken W2 dof loop
      end do ! End of face loop
    end do

  end subroutine compute_trace_operator_code

end module compute_trace_operator_kernel_mod
