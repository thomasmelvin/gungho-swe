!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the departure points in 1D.
!>
!> The kernel computes the departure points using the departure wind. The
!> kernel works in 1D only.
!>
module calc_departure_point_kernel_mod

  use argument_mod,                only : arg_type,                &
                                          GH_FIELD, GH_REAL,       &
                                          GH_READ, GH_WRITE,       &
                                          GH_SCALAR, GH_INTEGER,   &
                                          CELL_COLUMN
  use constants_mod,               only : r_def, i_def
  use fs_continuity_mod,           only : W2, W3
  use kernel_mod,                  only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_departure_point_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/              &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W3), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),   &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)       &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: calc_departure_point_code
  end type

  public :: calc_departure_point_code

contains

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief Subroutine to calculate the departure point
!! @param[in]    nlayers Number of model levels
!! @param[in,out] dep_pts Departure point values in W2 space
!! @param[in]    departure_pt_stencil_length  Length of stencil
!! @param[in]    n_iterations  Number of iterations
!! @param[in]    dt The model timestep length
!! @param[in]    undf_w2 Number of unique degrees of freedom for W2
!! @param[in]    ndf_w2 Number of degrees of freedom per cell in W2
!! @param[in]    stencil_map_w2  Stencil map in W2 space
!! @param[in]    undf_w3 Number of unique degrees of freedom for W3
!! @param[in]    ndf_w3 Number of degrees of freedom per cell in W3
!! @param[in]    stencil_map_w3  Stencil map in W3 space
!! @param[in]    cell_orientation Orientation of cells, inparticular halo cells
!! @param[in]    u_n Wind in W2 space at time n
!! @param[in]    u_np1 Wind in W2 space at time n+1
!! @param[in]    direction Direction in which to calculate departure points
!! @param[in]    dep_pt_method Enumeration of method to use
subroutine calc_departure_point_code( nlayers,                       &
                                      dep_pts,                       &
                                      departure_pt_stencil_length,   &
                                      n_iterations,                  &
                                      dt,                            &
                                      undf_w2,                       &
                                      ndf_w2,                        &
                                      stencil_map_w2,                &
                                      undf_w3,                       &
                                      ndf_w3,                        &
                                      stencil_map_w3,                &
                                      cell_orientation,              &
                                      u_n,                           &
                                      u_np1,                         &
                                      direction,                     &
                                      dep_pt_method )

  use departure_points_mod, only : calc_dep_point
  use cosmic_flux_mod,      only : calc_stencil_ordering, w2_dof, reorientate_w2field
  use flux_direction_mod,   only : x_direction, y_direction

  implicit none

  integer(kind=i_def), intent(in)         :: nlayers
  integer(kind=i_def), intent(in)         :: undf_w2
  real(kind=r_def), intent(inout)         :: dep_pts(1:undf_w2)
  integer(kind=i_def), intent(in)         :: departure_pt_stencil_length
  integer(kind=i_def), intent(in)         :: n_iterations
  integer(kind=i_def), intent(in)         :: ndf_w2
  integer(kind=i_def), intent(in)         :: undf_w3
  integer(kind=i_def), intent(in)         :: ndf_w3
  integer(kind=i_def), intent(in)         :: stencil_map_w2(1:ndf_w2,1:departure_pt_stencil_length)
  integer(kind=i_def), intent(in)         :: stencil_map_w3(1:ndf_w3,1:departure_pt_stencil_length)
  real(kind=r_def), intent(in)            :: u_n(1:undf_w2)
  real(kind=r_def), intent(in)            :: u_np1(1:undf_w2)
  real(kind=r_def), intent(in)            :: cell_orientation(1:undf_w3)
  integer(kind=i_def), intent(in)         :: direction
  integer(kind=i_def), intent(in)         :: dep_pt_method
  real(kind=r_def), intent(in)            :: dt

  real(kind=r_def)     :: xArrival
  real(kind=r_def)     :: u_n_local(1:departure_pt_stencil_length+1)
  real(kind=r_def)     :: u_np1_local(1:departure_pt_stencil_length+1)

  integer(kind=i_def)  :: nCellEdges
  integer(kind=i_def)  :: k, df1, df2, ii, jj
  integer(kind=i_def)  :: stencil_order_out(1:departure_pt_stencil_length)

  real(kind=r_def)     :: unordered_u_n(1:4,1:departure_pt_stencil_length)
  real(kind=r_def)     :: unordered_u_np1(1:4,1:departure_pt_stencil_length)
  real(kind=r_def)     :: orientated_unordered_u_n(1:4,1:departure_pt_stencil_length)
  real(kind=r_def)     :: orientated_unordered_u_np1(1:4,1:departure_pt_stencil_length)
  real(kind=r_def)     :: orientated_ordered_u_n(1:4,1:departure_pt_stencil_length)
  real(kind=r_def)     :: orientated_ordered_u_np1(1:4,1:departure_pt_stencil_length)

  nCellEdges = departure_pt_stencil_length+1

  call calc_stencil_ordering(departure_pt_stencil_length,stencil_order_out)

  if (direction == x_direction) then
    df1 = 1
    df2 = 3
  else if (direction == y_direction) then
    df1 = 2
    df2 = 4
  end if

  do k= 0, nlayers-1

    do ii = 1, departure_pt_stencil_length
      do jj = 1, 4
        unordered_u_n(jj,ii) = u_n(stencil_map_w2(jj,ii)+k)
        unordered_u_np1(jj,ii) = u_np1(stencil_map_w2(jj,ii)+k)
      end do
    end do

    do ii = 1, departure_pt_stencil_length
      orientated_unordered_u_n(:,ii) = reorientate_w2field(unordered_u_n(:,ii), int(cell_orientation(stencil_map_w3(1,ii))))
      orientated_unordered_u_np1(:,ii) = reorientate_w2field(unordered_u_np1(:,ii), int(cell_orientation(stencil_map_w3(1,ii))))
    end do

    do ii = 1, departure_pt_stencil_length
      orientated_ordered_u_n(:,ii) = orientated_unordered_u_n(:,stencil_order_out(ii))
      orientated_ordered_u_np1(:,ii) = orientated_unordered_u_np1(:,stencil_order_out(ii))
    end do

    nCellEdges = departure_pt_stencil_length+1
    do ii = 1, nCellEdges-1
      u_n_local(ii)      = orientated_ordered_u_n(df1,ii)
      u_np1_local(ii)    = orientated_ordered_u_np1(df1,ii)
    end do
    u_n_local(nCellEdges)      = orientated_ordered_u_n(df2,nCellEdges-1)
    u_np1_local(nCellEdges)    = orientated_ordered_u_np1(df2,nCellEdges-1)

    ! xArrival = 0.0 represents the arrival point at a flux edge in a finite
    ! element cell in the local coordinates
    xArrival = 0.0_r_def

    dep_pts(stencil_map_w2(df1,1)+k) = calc_dep_point(  xArrival,             &
                                                        nCellEdges,           &
                                                        u_n_local,            &
                                                        u_np1_local,          &
                                                        dt,                   &
                                                        dep_pt_method,        &
                                                        n_iterations )

    ! xArrival = 1.0 represents the opposite flux edge in a finite element cell
    ! in the local coordinates
    xArrival = 1.0_r_def

    dep_pts(stencil_map_w2(df2,1)+k) = calc_dep_point(  xArrival,             &
                                                        nCellEdges,           &
                                                        u_n_local,            &
                                                        u_np1_local,          &
                                                        dt,                   &
                                                        dep_pt_method,        &
                                                        n_iterations )

  end do

end subroutine calc_departure_point_code

end module calc_departure_point_kernel_mod
