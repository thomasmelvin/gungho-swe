!------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!------------------------------------------------------------------------------
!> @brief   Routines for calculating coefficients for subgrid rho representation.
!!
!! @details This module contains functions and subroutines which allow both
!!          linear and quadratic (PPM) representation of rho to be computed.
!------------------------------------------------------------------------------
module subgrid_rho_mod

use constants_mod, only: i_def, r_def, l_def, EPS

implicit none

private

public :: second_order_coeffs
public :: fourth_order_vertical_edge
public :: ppm_output

contains

  !----------------------------------------------------------------------------
  !> @brief  Minmod function which is a slope limiter for linear subgrid
  !!         representation of rho.
  !!
  !! @param[in]   a       Estimate of slope
  !! @param[in]   b       Estimate of slope
  !! @return      f       Output slope
  !----------------------------------------------------------------------------
  function minmod_function(a,b) result(f)

    implicit none

    real(kind=r_def), intent(in)    :: a
    real(kind=r_def), intent(in)    :: b
    real(kind=r_def)                :: f
    if ( (a*b) <= 0.0_r_def ) then
      f = 0.0_r_def
    else
      if (abs(a) < abs(b)) then
        f = a
      else
        f = b
      end if
    end if

  end function minmod_function

  !----------------------------------------------------------------------------
  !> @brief  Maxmod function which is a slope limiter for linear subgrid
  !!         representation of rho.
  !!
  !! @param[in]   a       Estimate of slope
  !! @param[in]   b       Estimate of slope
  !! @return      f       Output slope
  !----------------------------------------------------------------------------
  function maxmod_function(a,b) result(f)

    implicit none

    real(kind=r_def), intent(in)    :: a
    real(kind=r_def), intent(in)    :: b
    real(kind=r_def)                :: f
    if ( (a*b) <= 0.0_r_def ) then
      f = 0.0_r_def
    else
      if (abs(a) < abs(b)) then
        f = b
      else
        f = a
      end if
    end if

  end function maxmod_function

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a fourth-order interpolation.
  !> @details Uses a fourth-order interpolation to find the vertical cell edge
  !!          values of rho. The vertical grid spacing is used to compute the
  !!          mass, and a high-order polynomial is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge value.
  !!
  !> @param[in]   rho        Density values of four cells which have the ordering
  !!                         | 1 | 2 | 3 | 4 |
  !> @param[in]   dz         Height of each layer, with index the same as rho
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 | 3 | 4 |
  !!                         with edges  0   1   2   3   4
  !> @param[out]  edge_below The edge value located below layer k
  !!                         (layer k corresponds to cell 3 index above)
  !> @param[in]   positive   Ensures returned estimate of rho at the cell edge is
  !!                         positive (not yet implemented, see #1419)
  !> @param[in]   monotone   Ensures no over or undershoots are produced
  !!                         (not yet implemented, see #1419)
  !----------------------------------------------------------------------------
  subroutine fourth_order_vertical_edge(rho, dz, edge_to_do, positive, monotone, edge_below)

    implicit none

    real(kind=r_def),    intent(in)    :: rho(1:4)
    real(kind=r_def),    intent(in)    :: dz(1:4)
    integer(kind=i_def), intent(in)    :: edge_to_do
    logical(kind=l_def), intent(in)    :: positive
    logical(kind=l_def), intent(in)    :: monotone
    real(kind=r_def),    intent(out)   :: edge_below

    real(kind=r_def) :: z(0:4), dzs(1:4), dzsum, edge_height
    real(kind=r_def) :: dmass(1:4)
    real(kind=r_def) :: cmass(0:4)
    real(kind=r_def) :: poly_mass(1:4)
    real(kind=r_def) :: dl_dz(1:4)

    integer(kind=i_def) :: i

    ! Get scaling value
    dzsum = sum(dz)

    ! Get scaled dz
    dzs = dz / dzsum

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_def
    do i = 1, 4
      z(i) = z(i-1) + dzs(i)
    end do

    ! Get edge height to interpolate rho to
    edge_height = z(edge_to_do)

    ! Get mass scaled by height
    dmass = rho * dzs

    ! Get cumulative mass
    cmass(0) = 0.0_r_def
    do i = 1, 4
      cmass(i) = cmass(i-1) + dmass(i)
    end do

    ! Get cumulative mass divided by denominator of polynomial
    poly_mass(1) = cmass(1)/((z(1))*(z(1)-z(2))*(z(1)-z(3))*(z(1)-z(4)))
    poly_mass(2) = cmass(2)/((z(2))*(z(2)-z(1))*(z(2)-z(3))*(z(2)-z(4)))
    poly_mass(3) = cmass(3)/((z(3))*(z(3)-z(1))*(z(3)-z(2))*(z(3)-z(4)))
    poly_mass(4) = cmass(4)/((z(4))*(z(4)-z(1))*(z(4)-z(2))*(z(4)-z(3)))

    ! Calculate derivative of numerator of polynomial at edge height
    dl_dz    = 4.0_r_def*edge_height**3
    dl_dz(1) = dl_dz(1) - 3.0_r_def*(z(2)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_def*(z(3)*z(4) + z(2)*z(3) + z(2)*z(4))*edge_height - z(2)*z(3)*z(4)
    dl_dz(2) = dl_dz(2) - 3.0_r_def*(z(1)+z(3)+z(4))*edge_height**2 &
               + 2.0_r_def*(z(3)*z(4) + z(1)*z(3) + z(1)*z(4))*edge_height - z(1)*z(3)*z(4)
    dl_dz(3) = dl_dz(3) - 3.0_r_def*(z(1)+z(2)+z(4))*edge_height**2 &
               + 2.0_r_def*(z(2)*z(4) + z(1)*z(2) + z(1)*z(4))*edge_height - z(1)*z(2)*z(4)
    dl_dz(4) = dl_dz(4) - 3.0_r_def*(z(1)+z(2)+z(3))*edge_height**2 &
               + 2.0_r_def*(z(2)*z(3) + z(1)*z(2) + z(1)*z(3))*edge_height - z(1)*z(2)*z(3)

    ! Calculate value of edge below layer k
    edge_below = sum( poly_mass * dl_dz )

  end subroutine fourth_order_vertical_edge

  !----------------------------------------------------------------------------
  !> @brief  Returns the coefficients,a0,a1,a2 which are a quadratic
  !!         representation of rho within the cell, rho(x)=a0 + a1*x + a2*x^2
  !!         for 0<=x<=1. The dofmap for the density values is of the form
  !!         | 1 | 2 | 3 | 4 | 5 | where the subgrid coefficients are being
  !!         estimated for cell 3.
  !!
  !! @param[in]   density        Density values of five cells which have the
  !!                             ordering
  !!                             | 1 | 2 | 3 | 4 | 5 |
  !! @param[out]  coeffs         Coefficients for cell 3 with coeffs(1)=a0,
  !!                             coeffs(2)=a1, coeffs(3)=a2
  !! @param[in]   positive       Ensures returned estimate of rho at the cell
  !!                             edge is positive
  !! @param[in]   monotone       Ensures no over or undershoots are produced
  !----------------------------------------------------------------------------
  subroutine second_order_coeffs(density,coeffs,positive,monotone)

    implicit none

    real(kind=r_def), intent(in)    :: density(1:5)
    real(kind=r_def), intent(out)   :: coeffs(1:3)
    logical, intent(in)             :: positive
    logical, intent(in)             :: monotone

    real(kind=r_def)                :: cell_widths(1:5)
    real(kind=r_def)                :: density_cell_edge_left
    real(kind=r_def)                :: density_cell_edge_right

    coeffs(:) = 0.0_r_def
    cell_widths(:) = (/ 1.0_r_def,1.0_r_def,1.0_r_def,1.0_r_def,1.0_r_def /)

    density_cell_edge_left = calc_density_at_cell_edge(cell_widths(1:4),density(1:4),positive,monotone)
    density_cell_edge_right = calc_density_at_cell_edge(cell_widths(2:5),density(2:5),positive,monotone)
    call ppm_output(density_cell_edge_left,density_cell_edge_right,density(3),monotone,coeffs)

  end subroutine second_order_coeffs

  !----------------------------------------------------------------------------
  !> @brief  Calculates the estimated density at the edge of a cell required for
  !!         using PPM to estimate the quadratic subgrid representation of rho.
  !!         The function is passed four density values from consecutive cells
  !!         (which all lie in the same direction) with the dofmap
  !!         | 1 | 2 | 3 | 4 | and returns the estimated density value between
  !!         cells 2 and 3.
  !!         Positivity and monotonicity options are provided.
  !!
  !! @param[in]   cell_widths        All are assumed to be equal to 1.0 (computational domain)
  !! @param[in]   density            Has dof map of the form | 1 | 2 | 3 | 4 |
  !! @param[in]   positive           Ensures returned estimate of rho at the cell edge is positive
  !! @param[in]   monotone           Ensures no over or undershoots are produced
  !! @return      density_at_edge    Coefficients for cell 3 with coeffs(1)=a0,
  !!                                 coeffs(2)=a1, coeffs(3)=a2
  !----------------------------------------------------------------------------
  function calc_density_at_cell_edge(cell_widths,density,positive,monotone) result(density_at_edge)

    implicit none

    real(kind=r_def), intent(in)  :: cell_widths(1:4)
    real(kind=r_def), intent(in)  :: density(1:4)
    logical, intent(in)           :: positive
    logical, intent(in)           :: monotone

    real(kind=r_def) :: density_at_edge
    real(kind=r_def) :: mass(1:4)
    real(kind=r_def) :: y1,y2,y3,y4,y
    real(kind=r_def) :: m1,m2,m3,m4
    real(kind=r_def) :: t1,t2,t3


    density_at_edge= 0.0_r_def

    ! Convert density to mass although this is being done on the computational grid where
    ! cell widths are assumed equal to 1.0.
    mass = cell_widths*density

    y1 = cell_widths(1)
    y2 = y1+cell_widths(2)
    y3 = y2+cell_widths(3)
    y4 = y3+cell_widths(4)

    m1 = mass(1)
    m2 = mass(2) + m1
    m3 = mass(3) + m2
    m4 = mass(4) + m3

    m1 = m1 / (y1*(y1-y2)*(y1-y3)*(y1-y4))
    m2 = m2 / (y2*(y2-y1)*(y2-y3)*(y2-y4))
    m3 = m3 / (y3*(y3-y1)*(y3-y2)*(y3-y4))
    m4 = m4 / (y4*(y4-y1)*(y4-y2)*(y4-y3))

    y = cell_widths(1) + cell_widths(2)

    density_at_edge = &
        y2*(y2-y3)*(y2-y4)*m1 + &
        ((y2-y1)*(y2-y3)*(y2-y4) + y2*(y2-y3)*(y2-y4) + y2*(y2-y1)*(y2-y4) + y2*(y2-y1)*(y2-y3))*m2 + &
        y2*(y2-y1)*(y2-y4)*m3 + &
        y2*(y2-y1)*(y2-y3)*m4

!     If cell_widths=(/ 1.0,1.0,1.0,1.0,1.0/) then the above equation reduces to the equation below which agrees
!     with Colella and Woodward,JCP 54,1984, equation (1.9), see ticket #441 for more details.
!     density_at_edge = (7.0_r_def/12.0_r_def)*(density(2)+density(3))-(1.0_r_def/12.0_r_def)*(density(1)+density(4))

    if (positive) then
      density_at_edge = max(density_at_edge,0.0_r_def)
    end if

    if ( monotone ) then
      t1 = ( density_at_edge - density(2) )*( density(3) - density_at_edge )
      t2 = ( density(2) - density(1) )*( density(4) - density(3) )
      t3 = ( density_at_edge - density(2) )*( density(2) - density(1) )
      if ( t1 < 0.0_r_def .AND. ( t2 > 0.0_r_def .OR. t3 < 0.0_r_def ) ) then
         t1 = min(density(3),density(2))
         t2 = max(density(3),density(2))
         density_at_edge = min( t2, max(density_at_edge,t1) )
      end if
    end if

  end function calc_density_at_cell_edge

  !----------------------------------------------------------------------------
  !> @brief  Outputs the coefficients (a0,a1,a2) for the subgrid representation
  !!         rho(x) = a0 + a1*x + a2*x^2. Inputs are the density value for the
  !!         cell, and the left hand and right hand estimates of the density
  !!         for the cell. Given these three values a quadratic subgrid
  !!         approximation of rho can be made.
  !!
  !! @param[in]   density_cell_edge_left   Estimate of the density at x=0
  !! @param[in]   density_cell_edge_right  Estimate of the density at x=1
  !! @param[in]   density_of_cell          Average density of the cell
  !! @param[in]   monotone                 Ensures no over or undershoots
  !! @param[out]  coeffs                   coeffs(1)=a0, coeffs(2)=a1, coeffs(3)=a2
  !----------------------------------------------------------------------------
  subroutine ppm_output(density_cell_edge_left,density_cell_edge_right,density_of_cell,monotone,coeffs)

    implicit none

    real(kind=r_def), intent(in)    :: density_cell_edge_left
    real(kind=r_def), intent(in)    :: density_cell_edge_right
    real(kind=r_def), intent(in)    :: density_of_cell
    logical,          intent(in)    :: monotone
    real(kind=r_def), intent(out)   :: coeffs(1:3)

    real(kind=r_def) :: t1,t2,t3


    ! Calculate coefficients
    coeffs(1) = density_cell_edge_left
    coeffs(2) = -4.0_r_def*density_cell_edge_left - 2.0_r_def*density_cell_edge_right + 6.0_r_def*density_of_cell
    coeffs(3) =  3.0_r_def*density_cell_edge_left + 3.0_r_def*density_cell_edge_right - 6.0_r_def*density_of_cell

    !
    ! Check for subgrid monotonicity
    !
    if ( monotone ) then
      t1 = -0.5_r_def*coeffs(2)/(coeffs(3) + EPS)
      if ((t1+EPS)*(1.0_r_def+EPS-t1) > 0.0_r_def) then
        t2 = (density_cell_edge_right-density_of_cell) * (density_of_cell-density_cell_edge_left)
        t3 = abs(density_of_cell-density_cell_edge_left) - abs(density_cell_edge_right-density_of_cell)
        if ( t2 < EPS ) then
          coeffs(1) = density_of_cell
          coeffs(2) = 0.0_r_def
          coeffs(3) = 0.0_r_def
        else
          if ( t3 < 0.0_r_def ) then
            coeffs(1) = density_cell_edge_left
            coeffs(2) = 0.0_r_def
            coeffs(3) = 3.0_r_def*(density_of_cell - density_cell_edge_left)
          else
            coeffs(1) = -2.0_r_def*density_cell_edge_right + 3.0_r_def*density_of_cell
            coeffs(2) =  6.0_r_def*density_cell_edge_right - 6.0_r_def*density_of_cell
            coeffs(3) = -3.0_r_def*density_cell_edge_right + 3.0_r_def*density_of_cell
          end if
        end if
      end if
    end if

  end subroutine ppm_output

end module subgrid_rho_mod
