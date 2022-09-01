!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the departure distances for cell faces in the
!!        vertical direction.
!> @details This code calculates the distance which is swept through a cell
!!          in the z direction during one timestep. The arrival point is the
!!          cell face and the departure point is calculated. Options for the
!!          calculation of the departure point are a single Euler timestep, the
!!          midpoint rule, the trapezoidal rule, or time-averaging. This kernel
!!          returns the distance between the arrival and departure point for
!!          each cell face in the vertical and this value is positive if the
!!          w wind (radial wind) is positive, i.e. increasing in height. The
!!          integer part of the distance is the number of complete cells and
!!          the fractional part is the fraction of the cell in which the
!!          departure point resides.

module vertical_deppt_kernel_mod

use argument_mod,                only : arg_type,              &
                                        GH_FIELD, GH_REAL,     &
                                        GH_INC, GH_READ,       &
                                        GH_SCALAR, GH_INTEGER, &
                                        CELL_COLUMN
use fs_continuity_mod,           only : W2
use constants_mod,               only : r_def, i_def
use kernel_mod,                  only : kernel_type
use departure_points_config_mod, only : vertical_limit,          &
                                        vertical_limit_boundary, &
                                        vertical_limit_exponential
implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vertical_deppt_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/             &
       arg_type(GH_FIELD,  GH_REAL, GH_INC,  W2), & ! dep_pts
       arg_type(GH_FIELD,  GH_REAL, GH_INC,  W2), & ! cfl
       arg_type(GH_FIELD,  GH_REAL, GH_READ, W2), & ! u_n
       arg_type(GH_FIELD,  GH_REAL, GH_READ, W2), & ! u_np1
       arg_type(GH_FIELD,  GH_REAL, GH_READ, W2), & ! heights_w2
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),  & ! iterations
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),  & ! method
       arg_type(GH_SCALAR, GH_REAL, GH_READ)      & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_deppt_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vertical_deppt_code

contains

!> @brief Kernel which computes the departure distances for cell faces in the
!!        vertical direction.
!> @param[in]     nlayers             The number of layers
!> @param[in,out] dep_pts_z           The departure distances in the vertical
!> @param[in,out] cfl                 The vertical CFL calculated via the
!!                                    vertical departure points
!> @param[in]     u_n                 The wind field at time level n
!> @param[in]     u_np1               The wind field at time level n+1
!> @param[in]     height_w2           Physical height of W2 dofs
!> @param[in]     n_dep_pt_iterations The number of departure point iterations
!> @param[in]     vertical_method     Enumerator for the vertical method to be
!!                                    used for computing the departure points
!> @param[in]     dt                  The model timestep length
!> @param[in]     ndf_w2              The number of degrees of freedom per cell
!> @param[in]     undf_w2             The number of unique degrees of freedom
!> @param[in]     map_w2              The dofmap for the cell at the base of the column
subroutine vertical_deppt_code(  nlayers,             &
                                 dep_pts_z,           &
                                 cfl,                 &
                                 u_n,                 &
                                 u_np1,               &
                                 height_w2,           &
                                 n_dep_pt_iterations, &
                                 vertical_method,     &
                                 dt,                  &
                                 ndf_w2,              &
                                 undf_w2,             &
                                 map_w2 )

  use departure_points_mod,        only : calc_vertical_dep_cfl, &
                                          vertical_increasing_check

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: ndf_w2
  integer(kind=i_def),                     intent(in)    :: undf_w2
  integer(kind=i_def),                     intent(in)    :: vertical_method
  integer(kind=i_def), dimension(ndf_w2),  intent(in)    :: map_w2
  real(kind=r_def),    dimension(undf_w2), intent(in)    :: u_n
  real(kind=r_def),    dimension(undf_w2), intent(in)    :: u_np1
  real(kind=r_def),    dimension(undf_w2), intent(in)    :: height_w2
  real(kind=r_def),                        intent(in)    :: dt
  real(kind=r_def),    dimension(undf_w2), intent(inout) :: dep_pts_z
  real(kind=r_def),    dimension(undf_w2), intent(inout) :: cfl

  integer(kind=i_def), intent(in) :: n_dep_pt_iterations

  integer(kind=i_def) :: k

  integer(kind=i_def) :: nCellEdges
  real(kind=r_def)    :: xArrival_comp
  real(kind=r_def)    :: xArrival_phys
  real(kind=r_def)    :: u_n_local(1:nlayers+1)
  real(kind=r_def)    :: u_np1_local(1:nlayers+1)
  real(kind=r_def)    :: u_dep(1:nlayers+1)
  real(kind=r_def)    :: height_local(1:nlayers+1)
  real(kind=r_def)    :: dz_local(1:nlayers)
  real(kind=r_def)    :: upwind_dz_n
  real(kind=r_def)    :: upwind_dz_np1
  real(kind=r_def)    :: dep_local(1:nlayers-1)
  real(kind=r_def)    :: cfl_local

  ! Number of cell edgs
  nCellEdges = nlayers+1

  ! Initialise local variables to zero
  u_dep        = 0.0_r_def
  u_n_local    = 0.0_r_def
  u_np1_local  = 0.0_r_def
  height_local = 0.0_r_def
  dz_local     = 0.0_r_def
  dep_local    = 0.0_r_def

  ! Set up dz (the spacing between vertical W2 dofs)
  do k=0,nlayers-1
    dz_local(k+1) = height_w2(map_w2(6)+k)-height_w2(map_w2(5)+k)
  end do

  ! Compute the physical velocity as departure_wind*upwind_dz
  do k=1,nlayers-1
    upwind_dz_n   = ( 0.5_r_def + sign( 0.5_r_def, u_n(map_w2(5)+k) ) ) * dz_local(k) + &
                    ( 0.5_r_def - sign( 0.5_r_def, u_n(map_w2(5)+k) ) ) * dz_local(k+1)
    upwind_dz_np1 = ( 0.5_r_def + sign( 0.5_r_def, u_np1(map_w2(5)+k) ) ) * dz_local(k) + &
                    ( 0.5_r_def - sign( 0.5_r_def, u_np1(map_w2(5)+k) ) ) * dz_local(k+1)
    u_dep(k+1)       = 0.5_r_def * ( u_n(map_w2(5)+k) + u_np1(map_w2(5)+k) )
    u_n_local(k+1)   = u_n(map_w2(5)+k) * upwind_dz_n
    u_np1_local(k+1) = u_np1(map_w2(5)+k) * upwind_dz_np1
    height_local(k+1) = height_w2(map_w2(5)+k)
  end do
  ! Apply vertical boundary conditions
  height_local(1) = height_w2(map_w2(5))
  u_dep(1)        = 0.0_r_def
  u_n_local(1)    = 0.0_r_def
  u_np1_local(1)  = 0.0_r_def
  height_local(nCellEdges) = height_w2(map_w2(6)+nlayers-1)
  u_dep(nCellEdges)        = 0.0_r_def
  u_n_local(nCellEdges)    = 0.0_r_def
  u_np1_local(nCellEdges)  = 0.0_r_def

  ! Apply vertical boundary conditions to the departure points.
  dep_pts_z( map_w2(5) ) =  0.0_r_def
  dep_pts_z( map_w2(6)+nlayers-1 ) =  0.0_r_def

  ! Loop over all layers except the bottom layer.
  ! This code is hard-wired to work with 6 W2 dofs per cell where dof=5 is the
  ! vertical dof at the bottom of the cell.
  do k=1,nlayers-1
    xArrival_comp = real(k,r_def)
    xArrival_phys = height_w2(map_w2(5)+k)
    call calc_vertical_dep_cfl( xArrival_comp,        &
                                xArrival_phys,        &
                                nCellEdges,           &
                                u_n_local,            &
                                u_np1_local,          &
                                u_dep,                &
                                height_local,         &
                                dt,                   &
                                vertical_method,      &
                                n_dep_pt_iterations,  &
                                vertical_limit,       &
                                dep_local(k),         &
                                cfl_local )
    cfl( map_w2(5) + k ) =  cfl_local
  end do

  if ( vertical_limit == vertical_limit_exponential ) then
    ! Ensure vertical departure points are monotonic
    call vertical_increasing_check(dep_local, nlayers)
  end if

  do k=1,nlayers-1
    xArrival_comp = real(k,r_def)
    dep_pts_z( map_w2(5) + k ) =  xArrival_comp - dep_local(k)
  end do

end subroutine vertical_deppt_code

end module vertical_deppt_kernel_mod
