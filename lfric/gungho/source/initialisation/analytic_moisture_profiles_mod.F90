!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of moisture fields
!!          at a given point based upon a specified analytic formula
module analytic_moisture_profiles_mod

use constants_mod,                only : r_def, i_def, pi
use log_mod,                      only : log_event,                &
                                         log_scratch_space,        &
                                         LOG_LEVEL_ERROR
use idealised_config_mod,         only : test_grabowski_clark
use physics_common_mod,           only : qsaturation
use planet_config_mod,            only : recip_epsilon

implicit none

private

public :: analytic_moisture

contains

!> @brief Compute an analytic moisture field
!> @param[in] chi         Position in Cartesian coordinates
!> @param[in] temperature Air temperature in K
!> @param[in] pressure    Air pressure in Pa
!> @param[in] choice      Integer defining which specified formula to use
!> @result moisture The resulting moisture field
function analytic_moisture(chi, temperature, pressure, choice) result(moisture)

  implicit none
  real(kind=r_def),    intent(in) :: chi(3)
  real(kind=r_def),    intent(in) :: temperature
  real(kind=r_def),    intent(in) :: pressure
  integer(kind=i_def), intent(in) :: choice
  real(kind=r_def)                :: moisture
  real(kind=r_def)                :: r, r1, r2, xc, zc  ! Spatial distances
  real(kind=r_def)                :: h0, rel_hum        ! Relative humidities
  real(kind=r_def)                :: mr_sat             ! Saturation value

  ! We'll always need mr_sat so find it right away. Pressure is needed in mbar
  mr_sat = qsaturation(temperature, 0.01_r_def*pressure)

  select case( choice )

  ! Test from Grabowski and Clark (1991)
  ! Returns the relative humidity field
  case( test_grabowski_clark )
    ! Parameters
    xc = 0.0_r_def    ! central x position of bubble
    zc = 800.0_r_def  ! central z position of bubble
    r1 = 300.0_r_def  ! outer radius of bubble
    r2 = 200.0_r_def  ! inner radius of bubble
    h0 = 0.2_r_def    ! background relative humidity

    r = sqrt((chi(1)-xc)**2.0 + (chi(3)-zc)**2.0)
    if ( r <= r2 ) then
      rel_hum = 1.0_r_def
    else if ( r <= r1 ) then
      rel_hum = h0 + (1.0_r_def - h0) &
                      *( cos( pi*(r-r2) / (2.0_r_def*(r1-r2)) ) )**2.0_r_def
    else
      rel_hum = h0
    end if

    ! Now invert relative humidity expression to get water vapour mixing ratio
    moisture = rel_hum * mr_sat / &
               ( 1.0_r_def + (1.0_r_def - rel_hum)*mr_sat*recip_epsilon)

  case default
    ! In other cases, mixing ratio is set just under saturation value
    moisture = 0.99_r_def*mr_sat

  end select

end function analytic_moisture

end module analytic_moisture_profiles_mod
