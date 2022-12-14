!------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!------------------------------------------------------------------------------
!> @brief  Routines for calculating departure points in 1D used for the
!!         split advection scheme.
!!
!------------------------------------------------------------------------------
module departure_points_mod

use constants_mod, only : r_def, i_def
use log_mod,       only : log_event, LOG_LEVEL_ERROR, log_scratch_space
use departure_points_config_mod, only : method_euler,                &
                                        method_midpoint,             &
                                        method_trapezoidal,          &
                                        method_timeaverage,          &
                                        vertical_method_euler,       &
                                        vertical_method_midpoint,    &
                                        vertical_method_trapezoidal, &
                                        vertical_method_timeaverage, &
                                        vertical_limit,              &
                                        vertical_limit_boundary,     &
                                        vertical_limit_exponential

implicit none

private

public :: calc_dep_point
public :: calc_vertical_dep_cfl
public :: calc_uniform_vertical_dep_cfl
! The following subroutines are public in order to facilitate unit testing
public :: find_local_x_value
public :: find_local_vertical_value
public :: calc_u_at_x
public :: calc_u_in_vertical_comp
public :: calc_u_in_vertical_phys
public :: vertical_increasing_check

contains

  !----------------------------------------------------------------------------
  !> @brief  Calculates the distance between the arrival point and the
  !!         departure point in 1D. Note that the distance has sign (+/-) and
  !!         positive values represent the case when the wind is positive, such
  !!         that x_departure < x_arrival. The distance is negative if the wind
  !!         is negative.
  !!
  !! @param[in]   x_arrival    Arrival point in departure point calculation
  !! @param[in]   nCellEdges   Number of velocity values
  !! @param[in]   u_n          Velocity at cell edges at time n
  !! @param[in]   u_np1        Velocity at cell edges at time n+1
  !! @param[in]   deltaT       Time step length
  !! @param[in]   method       Integration method
  !! @param[in]   n_dep_pt_iterations Number of solver iterations
  !! @return      distance     Distance between arrival point and departure
  !!                           point
  !----------------------------------------------------------------------------
  function calc_dep_point(  x_arrival,           &
                            nCellEdges,          &
                            u_n,                 &
                            u_np1,               &
                            deltaT,              &
                            method,              &
                            n_dep_pt_iterations )  result(distance)

    implicit none

    real(kind=r_def), intent(in)    :: x_arrival
    integer(kind=i_def), intent(in) :: nCellEdges
    real(kind=r_def), intent(in)    :: u_n(1:nCellEdges)
    real(kind=r_def), intent(in)    :: u_np1(1:nCellEdges)
    real(kind=r_def), intent(in)    :: deltaT
    integer(kind=i_def), intent(in) :: method
    integer(kind=i_def), intent(in) :: n_dep_pt_iterations
    real(kind=r_def)                :: distance

    real(kind=r_def) :: u_arrival
    real(kind=r_def) :: u_at_midpoint
    real(kind=r_def) :: u_departure
    real(kind=r_def) :: x_at_mid_point
    real(kind=r_def) :: left_limit
    real(kind=r_def) :: right_limit
    real(kind=r_def) :: x_departure
    real(kind=r_def) :: u_arrival_np

    integer(kind=i_def) :: iLoop

    x_departure = x_arrival

    left_limit = real(-nCellEdges/2_i_def + 1_i_def, r_def)
    right_limit = real(nCellEdges/2_i_def, r_def)
    call test_value_in_limits(x_arrival,left_limit,right_limit)

    select case (method)

    case(method_euler) ! Euler's method

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

    case(method_trapezoidal) ! Trapezoidal

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

      do iLoop=1,n_dep_pt_iterations
        u_departure = calc_u_at_x(x_departure,nCellEdges,u_n)
        x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_departure)
        call test_value_in_limits(x_departure,left_limit,right_limit)
      end do

    case(method_midpoint) ! Mid-point

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*u_arrival
      call test_value_in_limits(x_departure,left_limit,right_limit)

      do iLoop=1,n_dep_pt_iterations
        x_at_mid_point = 0.5_r_def*(x_departure+x_arrival)
        u_at_midpoint = calc_u_at_x(x_at_mid_point,nCellEdges,u_n)
        x_departure = x_arrival - deltaT*u_at_midpoint
        call test_value_in_limits(x_departure,left_limit,right_limit)
      end do

    case(method_timeaverage) ! Time averaged velocity at arrival point

      u_arrival = calc_u_at_x(x_arrival,nCellEdges,u_n)
      u_arrival_np = calc_u_at_x(x_arrival,nCellEdges,u_np1)
      x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_arrival_np)
      call test_value_in_limits(x_departure,left_limit,right_limit)

    case default
      call log_event( " Departure point method undefined ", LOG_LEVEL_ERROR )
    end select

    distance = x_arrival - x_departure

  end function calc_dep_point


  !----------------------------------------------------------------------------
  !> @brief  Calculates the departure point in 1D in the vertical, and the
  !!         vertical CFL from the departure distance, using a physical vertical
  !!         grid. If the departure point exceeds the vertical domain it is
  !!         modified to lie within the domain but the CFL is unaffected.
  !!         The departure distance in computational space is calculated as
  !!         distance = arrival_point - departure_point.
  !!         Note that the distance has sign (+/-) and positive values represent
  !!         the case when the wind is positive, such that x_departure < x_arrival.
  !!         The distance is negative if the wind is negative.
  !!
  !! @param[in]     x_arrival_comp      Arrival point in computational space
  !! @param[in]     x_arrival_phys      Arrival point in physical space
  !! @param[in]     nCellEdges          Number of velocity values
  !! @param[in]     u_n                 Velocity*dz at cell edges at time n
  !! @param[in]     u_np1               Velocity*dz at cell edges at time n+1
  !! @param[in]     u_dep               Departure velocity at cell edges averaged
  !!                                    over time step
  !! @param[in]     height              Height of cell edges
  !! @param[in]     deltaT              Time step length
  !! @param[in]     vertical_method     Integration method
  !! @param[in]     n_dep_pt_iterations Number of solver iterations
  !! @param[in]     vertical_limit      Method to force departure points within the
  !!                                    vertical domain
  !! @param[in,out] departure_point     The departure point in reference space
  !!                                    (integer part is number of cells, fractional
  !!                                    part is fraction of departure cell)
  !! @param[in,out] cfl                 The CFL at the arrival point based on the
  !!                                    departure distance
  !----------------------------------------------------------------------------
  subroutine calc_vertical_dep_cfl( x_arrival_comp,         &
                                    x_arrival_phys,         &
                                    nCellEdges,             &
                                    u_n,                    &
                                    u_np1,                  &
                                    u_dep,                  &
                                    height,                 &
                                    deltaT,                 &
                                    vertical_method,        &
                                    n_dep_pt_iterations,    &
                                    vertical_limit,         &
                                    departure_point,        &
                                    cfl )

    implicit none

    real(kind=r_def), intent(in)        :: x_arrival_comp
    real(kind=r_def), intent(in)        :: x_arrival_phys
    integer(kind=i_def), intent(in)     :: nCellEdges
    real(kind=r_def), intent(in)        :: u_n(1:nCellEdges)
    real(kind=r_def), intent(in)        :: u_np1(1:nCellEdges)
    real(kind=r_def), intent(in)        :: u_dep(1:nCellEdges)
    real(kind=r_def), intent(in)        :: height(1:nCellEdges)
    real(kind=r_def), intent(in)        :: deltaT
    integer(kind=i_def), intent(in)     :: vertical_method
    integer(kind=i_def), intent(in)     :: n_dep_pt_iterations
    integer(kind=i_def), intent(in)     :: vertical_limit
    real(kind=r_def), intent(inout)     :: departure_point
    real(kind=r_def), intent(inout)     :: cfl

    real(kind=r_def) :: u_arrival
    real(kind=r_def) :: u_arrival_np
    real(kind=r_def) :: u_at_midpoint
    real(kind=r_def) :: u_departure
    real(kind=r_def) :: x_at_mid_point
    real(kind=r_def) :: left_limit
    real(kind=r_def) :: right_limit
    real(kind=r_def) :: x_departure
    real(kind=r_def) :: x_departure_comp
    real(kind=r_def) :: x_departure_corr
    real(kind=r_def) :: nn, dx, x0, x1

    integer(kind=i_def) :: iLoop

    ! Check computational arrival point is in range
    left_limit = 0.0_r_def
    right_limit = real(nCellEdges-1,r_def)
    call test_value_in_limits(x_arrival_comp,left_limit,right_limit)

    ! Compute physical departure point
    select case (vertical_method)

    case(vertical_method_euler) ! Euler's method

        u_arrival = calc_u_in_vertical_comp(x_arrival_comp,nCellEdges,u_np1)
        x_departure = x_arrival_phys - deltaT*u_arrival

    case(vertical_method_trapezoidal) ! Trapezoidal

        u_arrival = calc_u_in_vertical_comp(x_arrival_comp,nCellEdges,u_np1)
        x_departure = x_arrival_phys - deltaT*u_arrival

        do iLoop=1,n_dep_pt_iterations
          u_departure = calc_u_in_vertical_phys(x_departure,nCellEdges,u_n,height)
          x_departure = x_arrival_phys - deltaT*0.5_r_def*(u_arrival+u_departure)
        end do

    case(vertical_method_midpoint) ! Mid-point

        u_arrival = calc_u_in_vertical_comp(x_arrival_comp,nCellEdges,u_np1)
        x_departure = x_arrival_phys - deltaT*u_arrival

        do iLoop=1,n_dep_pt_iterations
          x_at_mid_point = 0.5_r_def*(x_departure+x_arrival_phys)
          u_at_midpoint = calc_u_in_vertical_phys(x_at_mid_point,nCellEdges,u_n,height)
          x_departure = x_arrival_phys - deltaT*u_at_midpoint
        end do

    case(vertical_method_timeaverage) ! Time averaged velocity at arrival point

        u_arrival = calc_u_in_vertical_comp(x_arrival_comp,nCellEdges,u_n)
        u_arrival_np = calc_u_in_vertical_comp(x_arrival_comp,nCellEdges,u_np1)
        x_departure = x_arrival_phys - deltaT*0.5_r_def*(u_arrival+u_arrival_np)

    case default
        call log_event( " Vertical departure point method undefined ", LOG_LEVEL_ERROR )
    end select

    ! Convert to computational departure distance ensuring it is within range

    nn = real(nCellEdges-1,r_def)
    x_departure_comp = 0.0_r_def

    if (x_departure < height(1)) then
      select case (vertical_limit)
      case(vertical_limit_boundary)
        ! Set the computational departure point to zero
        x_departure_comp = 0.0_r_def
      case(vertical_limit_exponential)
        ! Use exponential function to force departure point into the domain
        dx = height(2)-height(1)
        x_departure_corr = (x_arrival_phys-x_departure)/dx
        x_departure_comp = max( exp(-x_departure_corr/x_arrival_comp), 0.0_r_def )
      end select
      u_arrival = calc_u_in_vertical_comp(x_arrival_comp,nCellEdges,u_dep)
      cfl = deltaT*u_arrival
    else if (x_departure >= height(nCellEdges)) then
      select case (vertical_limit)
      case(vertical_limit_boundary)
        ! Set the computational departure point to be the top level
        x_departure_comp = real(nCellEdges-1,r_def)
      case(vertical_limit_exponential)
        ! Use exponential function to force departure point into the domain
        dx = height(nCellEdges)-height(nCellEdges-1)
        x_departure_corr = (x_arrival_phys-x_departure)/dx
        x_departure_comp = min(nn-exp( x_departure_corr/(nn-x_arrival_comp) ),nn)
      end select
      u_arrival = calc_u_in_vertical_comp(x_arrival_comp,nCellEdges,u_dep)
      cfl = deltaT*u_arrival
    else
      do iloop=1,nCellEdges-1
        if (x_departure .le. height(iloop+1) .AND. x_departure .gt. height(iloop)) then
            x0 = height(iloop)
            x1 = height(iloop+1)
            x_departure_comp = real(iloop-1, r_def)  + (x_departure-x0)/(x1-x0)
        end if
      end do
      cfl = x_arrival_comp - x_departure_comp
    end if

    departure_point = x_departure_comp

  end subroutine calc_vertical_dep_cfl


  !----------------------------------------------------------------------------
  !> @brief  Calculates the departure point in 1D in the vertical, and the
  !!         vertical CFL from the departure distance, using a uniform grid.
  !!         If the departure point exceeds the vertical domain it is modified
  !!         to lie within the domain but the CFL is unaffected.
  !!         The departure distance in computational space is calculated as
  !!         distance = arrival_point - departure_point.
  !!         Note that the distance has sign (+/-) and positive values represent
  !!         the case when the wind is positive, such that x_departure < x_arrival.
  !!         The distance is negative if the wind is negative.
  !!
  !! @param[in]     x_arrival           Arrival point in computational space
  !! @param[in]     nCellEdges          Number of velocity values
  !! @param[in]     u_n                 Velocity at cell edges at time n
  !! @param[in]     u_np1               Velocity at cell edges at time n+1
  !! @param[in]     deltaT              Time step length
  !! @param[in]     vertical_method     Integration method
  !! @param[in]     n_dep_pt_iterations Number of solver iterations
  !! @param[in]     vertical_limit      Method to force departure points within the
  !!                                    vertical domain
  !! @param[in,out] departure_point     The departure point in reference space
  !!                                    (integer part is number of cells, fractional
  !!                                    part is fraction of departure cell)
  !! @param[in,out] cfl                 The CFL at the arrival point based on the
  !!                                    departure distance
  !----------------------------------------------------------------------------
  subroutine calc_uniform_vertical_dep_cfl( x_arrival,              &
                                            nCellEdges,             &
                                            u_n,                    &
                                            u_np1,                  &
                                            deltaT,                 &
                                            vertical_method,        &
                                            n_dep_pt_iterations,    &
                                            vertical_limit,         &
                                            departure_point,        &
                                            cfl )

    implicit none

    real(kind=r_def), intent(in)        :: x_arrival
    integer(kind=i_def), intent(in)     :: nCellEdges
    real(kind=r_def), intent(in)        :: u_n(1:nCellEdges)
    real(kind=r_def), intent(in)        :: u_np1(1:nCellEdges)
    real(kind=r_def), intent(in)        :: deltaT
    integer(kind=i_def), intent(in)     :: vertical_method
    integer(kind=i_def), intent(in)     :: n_dep_pt_iterations
    integer(kind=i_def), intent(in)     :: vertical_limit
    real(kind=r_def), intent(inout)     :: departure_point
    real(kind=r_def), intent(inout)     :: cfl

    real(kind=r_def) :: u_arrival
    real(kind=r_def) :: u_arrival_np
    real(kind=r_def) :: u_at_midpoint
    real(kind=r_def) :: u_departure
    real(kind=r_def) :: x_at_mid_point
    real(kind=r_def) :: left_limit
    real(kind=r_def) :: right_limit
    real(kind=r_def) :: x_departure

    integer(kind=i_def) :: iLoop

    ! Check computational arrival point is in range
    left_limit = 0.0_r_def
    right_limit = real(nCellEdges-1,r_def)
    call test_value_in_limits(x_arrival,left_limit,right_limit)

    select case (vertical_method)

    case(vertical_method_euler) ! Euler's method

        u_arrival = calc_u_in_vertical_comp(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*u_arrival

    case(vertical_method_trapezoidal) ! Trapezoidal

        u_arrival = calc_u_in_vertical_comp(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*u_arrival

        do iLoop=1,n_dep_pt_iterations
          u_departure = calc_u_in_vertical_comp(x_departure,nCellEdges,u_n)
          x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_departure)
        end do

    case(vertical_method_midpoint) ! Mid-point

        u_arrival = calc_u_in_vertical_comp(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*u_arrival

        do iLoop=1,n_dep_pt_iterations
          x_at_mid_point = 0.5_r_def*(x_departure+x_arrival)
          u_at_midpoint = calc_u_in_vertical_comp(x_at_mid_point,nCellEdges,u_n)
          x_departure = x_arrival - deltaT*u_at_midpoint
        end do

    case(vertical_method_timeaverage) ! Time averaged velocity at arrival point

        u_arrival = calc_u_in_vertical_comp(x_arrival,nCellEdges,u_n)
        u_arrival_np = calc_u_in_vertical_comp(x_arrival,nCellEdges,u_np1)
        x_departure = x_arrival - deltaT*0.5_r_def*(u_arrival+u_arrival_np)

    case default
        call log_event( " Vertical departure point method undefined ", LOG_LEVEL_ERROR )
    end select

    if (x_departure < 0.0_r_def .OR. x_departure > real(nCellEdges-1,r_def)) then
      call test_value_in_vertical_limits(x_departure,left_limit,right_limit,vertical_limit)
      cfl = deltaT*u_arrival
    else
      cfl = x_arrival - x_departure
    end if

    departure_point = x_departure

  end subroutine calc_uniform_vertical_dep_cfl


  !----------------------------------------------------------------------------
  !> @brief  Calculates the location of a value x_in within a given stencil of
  !!         length nCellEdges.
  !!
  !! @param[in]    x_in  Arrival value, typically equal to 0.0.
  !! @param[in]    nCellEdges  Number of cell edges in a stencil
  !! @param[out]   iEdge       Index of cell edge to the left of x_in
  !! @param[out]   fractional_x_value  Fractional value of x_in
  !----------------------------------------------------------------------------
  subroutine find_local_x_value(x_in,nCellEdges,iEdge,fractional_x_value)

    implicit none

    real(kind=r_def), intent(in)     :: x_in
    integer(kind=i_def), intent(in)  :: nCellEdges
    integer(kind=i_def), intent(out) :: iEdge
    real(kind=r_def), intent(out)    :: fractional_x_value

    ! Check that the number of CellEdges is even
    if (modulo(nCellEdges,2_i_def) == 1_i_def) then
      call log_event( " Stencil length is incorrect ", LOG_LEVEL_ERROR )
    end if

    iEdge = floor(x_in)+nCellEdges/2_i_def

    ! Calculate distance from nearest lefthand cell edge
    fractional_x_value = abs(x_in - floor(x_in))

    if (iEdge < 1_i_def .OR. iEdge > nCellEdges) then
      call log_event( " Error in find_local_x_value routine ", LOG_LEVEL_ERROR )
    end if

  end subroutine find_local_x_value


  !----------------------------------------------------------------------------
  !> @brief  Calculates integer part and fractional part of x_in and the
  !!         subroutine is typically used in the vertical direction.
  !!
  !! @param[in]    x_in  Arrival value, typically equal to 0.0.
  !! @param[in]    nCellEdges  Number of cell edges in a stencil
  !! @param[out]   iEdge       Index of cell edge to the left of x_in
  !! @param[out]   fractional_x_value  Fractional value of x_in
  !----------------------------------------------------------------------------
  subroutine find_local_vertical_value(x_in,nCellEdges,iEdge,fractional_x_value)

    implicit none

    real(kind=r_def),    intent(in)       :: x_in
    integer(kind=i_def), intent(in)       :: nCellEdges
    integer(kind=i_def), intent(out)      :: iEdge
    real(kind=r_def),    intent(out)      :: fractional_x_value

    iEdge = floor(x_in) + 1_i_def

    ! Calculate distance from nearest lefthand cell edge
    fractional_x_value = abs(x_in - floor(x_in))

    if (iEdge < 1_i_def .OR. iEdge > nCellEdges) then
      call log_event( " Error in find_local_vertical_value routine ", LOG_LEVEL_ERROR )
    end if

  end subroutine find_local_vertical_value

  !----------------------------------------------------------------------------
  !> @brief  Subroutine which checks whether departure points are outside the
  !!         domain of interest defined by the stencil length.
  !!
  !! @param[in]   x_in         X value to be tested
  !! @param[in]   left_limit   Left hand bound
  !! @param[in]   right_limit  Right hand bound
  !----------------------------------------------------------------------------
  subroutine test_value_in_limits(x_in,left_limit,right_limit)

    implicit none

    real(kind=r_def), intent(in) ::   x_in
    real(kind=r_def), intent(in) ::   left_limit
    real(kind=r_def), intent(in) ::   right_limit

    if (x_in < left_limit .OR. x_in > right_limit) then
      write(log_scratch_space, '(A,E12.4E3,A,2E12.4E3)')           &
         'Departure distance ', x_in,                              &
         ' is out of bounds. Limits are ', left_limit, right_limit
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end subroutine test_value_in_limits

  !----------------------------------------------------------------------------
  !> @brief  Subroutine which checks whether departure points are outside the
  !!         vertical domain defined by the stencil length. We either 1) follow
  !!         Wood et al. 2009 to make the departure point fall within the domain at the
  !!         lowest level, or 2) set the departure point to be the boundary value.
  !!
  !! @param[in,out] x_dep           X value to be tested
  !! @param[in]     lower_limit     Lower bound
  !! @param[in]     upper_limit     Upper bound
  !! @param[in]     vertical_limit  Method to enforce the vertical limits
  !----------------------------------------------------------------------------
  subroutine test_value_in_vertical_limits(x_dep,lower_limit,upper_limit,vertical_limit)

    implicit none

    real(kind=r_def),    intent(inout) ::   x_dep
    real(kind=r_def),    intent(in)    ::   lower_limit
    real(kind=r_def),    intent(in)    ::   upper_limit
    integer(kind=i_def), intent(in)    ::   vertical_limit

    if (x_dep > upper_limit) then
      ! Set departure point value to upper limit
      x_dep = upper_limit
    end if

    if (x_dep < lower_limit) then
      select case (vertical_limit)
      case(vertical_limit_boundary)
        ! Set the departure point to be the lower limit
        x_dep = lower_limit
      case(vertical_limit_exponential)
        ! Following Wood et al. (2009) we use the formula
        ! x_departure = x_arrival exp(-dt * average_velocity/x_arrival)
        ! At the lowest level x_arrival = 1
        ! In our departure point calculation x_arrival - x_departure = dt * average_velocity
        x_dep = exp(x_dep - 1.0_r_def)
      end select
    end if

  end subroutine test_value_in_vertical_limits

  !----------------------------------------------------------------------------
  !> @brief  Returns an interpolated wind field value at x_in.
  !!
  !! @param[in]    x_in        Position at which to interpolate wind
  !! @param[in]    nCellEdges  Number of values in the local u field
  !! @param[in]    u_wind      Wind values
  !! @return       u_out       Interpolated wind value
  !----------------------------------------------------------------------------
  function calc_u_at_x(x_in,nCellEdges,u_wind) result(u_out)

    implicit none

    real(kind=r_def), intent(in)    ::  x_in
    integer(kind=i_def), intent(in) ::  nCellEdges
    real(kind=r_def), intent(in)    ::  u_wind(1:nCellEdges)
    real(kind=r_def)                ::  u_out

    real(kind=r_def)    :: fractional_x_value
    integer(kind=i_def) :: iEdge
    integer(kind=i_def) :: iCellRight

    call find_local_x_value(x_in,nCellEdges,iEdge,fractional_x_value)

    if (iEdge==nCellEdges) then
      iCellRight=iEdge
    else
      iCellRight = iEdge + 1_i_def
    end if

    u_out = (1.0_r_def-fractional_x_value)*u_wind(iEdge) +                &
                                  fractional_x_value*u_wind(iCellRight)

  end function calc_u_at_x

  !----------------------------------------------------------------------------
  !> @brief  Returns an interpolated wind field value in the vertical direction
  !!         at x_in assuming uniform integer values for arrival points.
  !!
  !! @param[in]    x_in        Position at which to interpolate wind
  !! @param[in]    nCellEdges  Number of values in the local u field
  !! @param[in]    u_wind      Wind values
  !! @return       u_out       Interpolated wind value
  !----------------------------------------------------------------------------
  function calc_u_in_vertical_comp(x_in,nCellEdges,u_wind) result(u_out)

    implicit none

    real(kind=r_def), intent(in)       ::  x_in
    integer(kind=i_def), intent(in)    ::  nCellEdges
    real(kind=r_def), intent(in)       ::  u_wind(1:nCellEdges)
    real(kind=r_def)                   ::  u_out

    real(kind=r_def)    :: fractional_x_value
    integer(kind=i_def) :: iEdge
    integer(kind=i_def) :: iCellRight

    call find_local_vertical_value(x_in,nCellEdges,iEdge,fractional_x_value)

    if (iEdge==nCellEdges) then
      iCellRight=iEdge
    else
      iCellRight = iEdge + 1_i_def
    end if

    if (x_in .le. 0.0_r_def .OR. x_in .ge. real(nCellEdges-1,r_def) ) then
      u_out=0.0_r_def
    else
      u_out = (1.0_r_def-fractional_x_value)*u_wind(iEdge) +                &
                                  fractional_x_value*u_wind(iCellRight)
    end if

  end function calc_u_in_vertical_comp

  !----------------------------------------------------------------------------
  !> @brief  Returns an interpolated wind field value in the vertical direction
  !!         at x_in using physical vertical values.
  !!
  !! @param[in]    x_in        Position at which to interpolate wind
  !! @param[in]    nCellEdges  Number of values in the local u field
  !! @param[in]    u_wind      Wind values
  !! @param[in]    height      Height of edges
  !! @return       u_out       Interpolated wind value
  !----------------------------------------------------------------------------
  function calc_u_in_vertical_phys(x_in,nCellEdges,u_wind,height) result(u_out)

    implicit none

    real(kind=r_def), intent(in)    ::  x_in
    integer(kind=i_def), intent(in) ::  nCellEdges
    real(kind=r_def), intent(in)    ::  u_wind(1:nCellEdges)
    real(kind=r_def), intent(in)    ::  height(1:nCellEdges)
    real(kind=r_def)                ::  u_out

    real(kind=r_def)    :: fractional_x_value
    integer(kind=i_def) :: iEdge
    integer(kind=i_def) :: iloop
    integer(kind=i_def) :: iCellRight
    real(kind=r_def)    :: x0
    real(kind=r_def)    :: x1

    if (x_in .le. height(1) .OR. x_in .ge. height(nCellEdges)) then
      u_out=0.0_r_def
    else
      do iloop=1,nCellEdges-1
        if (x_in .le. height(iloop+1) .AND. x_in .gt. height(iloop)) then

            iEdge = iloop

            x0 = height(iloop)
            x1 = height(iloop+1)

            fractional_x_value = (x_in-x0)/(x1-x0)

            if (iEdge==nCellEdges) then
              iCellRight=iEdge
            else
              iCellRight = iEdge+1
            end if

            u_out = (1.0_r_def-fractional_x_value)*u_wind(iEdge) +  &
                                  fractional_x_value*u_wind(iCellRight)

        end if

      end do

    end if

  end function calc_u_in_vertical_phys

  !----------------------------------------------------------------------------
  !> @brief  Subroutine which checks whether departure points are increasing with
  !!         model level in a column.
  !!
  !! @param[in,out] dep_pts          The departure points in a column
  !! @param[in]     nlayers          The number of layers in a column
  !----------------------------------------------------------------------------
  subroutine vertical_increasing_check(dep_pts, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    ::  nlayers
    real(kind=r_def),    intent(inout) ::  dep_pts(1:nlayers-1)

    integer(kind=i_def) :: k

    do k = 1, nlayers-2

      if (dep_pts(k+1) < dep_pts(k)) then
        ! If departure point above is less than departure point below, set
        ! the above point equal to the below point plus epsilon
        dep_pts(k+1) = dep_pts(k) + epsilon(1.0_r_def)

        ! check epsilon doesn't push point outside the domain
        if (dep_pts(k+1) > real(nlayers,r_def)) then
          dep_pts(k+1) = real(nlayers,r_def)
        end if
      end if

    end do

  end subroutine vertical_increasing_check

end module departure_points_mod
