!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_pressure_profiles_mod

use constants_mod,              only : r_def, pi
use log_mod,                    only : log_event,                &
                                       log_scratch_space,        &
                                       LOG_LEVEL_ERROR
use coord_transform_mod,        only : xyz2llr, central_angle
use idealised_config_mod,       only : test_cold_bubble_x,           &
                                       test_cold_bubble_y,           &
                                       test_warm_bubble,             &
                                       test_warm_bubble_3d,          &
                                       test_gaussian_hill,           &
                                       test_cosine_hill,             &
                                       test_cosine_bell,             &
                                       test_yz_cosine_hill,          &
                                       test_slotted_cylinder,        &
                                       test_constant_field,          &
                                       test_cosine_stripe,           &
                                       test_vortex_field,            &
                                       test_hadley_like_dcmip,       &
                                       test_gravity_wave,            &
                                       test_solid_body_rotation,     &
                                       test_solid_body_rotation_alt, &
                                       test_deep_baroclinic_wave,    &
                                       test_isentropic,              &
                                       test_isot_atm,                &
                                       test_isot_cold_atm,           &
                                       test_const_lapse_rate,        &
                                       test_dry_cbl,                 &
                                       test_snow,                    &
                                       test_shallow_conv,            &
                                       test_cos_phi,                 &
                                       test_cosine_bubble,           &
                                       test_div_free_reversible,     &
                                       test_eternal_fountain,        &
                                       test_curl_free_reversible,    &
                                       test_rotational,              &
                                       test_translational,           &
                                       test_vertical_cylinder,       &
                                       test_specified_profiles,      &
                                       test_bryan_fritsch,           &
                                       test_grabowski_clark
use initial_density_config_mod, only : r1, x1, y1, z1, r2, x2, y2, z2, &
                                       tracer_max, tracer_background
use base_mesh_config_mod,       only : geometry, &
                                       geometry_spherical
use planet_config_mod,          only : p_zero, Rd, kappa, scaled_radius
use reference_profile_mod,      only : reference_profile
use analytic_temperature_profiles_mod, only: analytic_temperature
use deep_baroclinic_wave_mod,   only : deep_baroclinic_wave

implicit none

private

public :: vortex_field
public :: analytic_pressure

contains

  !>@brief Compute the vortex field from Nair and Jablonowski 08
  !>@details  Equations below have been taken from Nair and Jablonowski, "Moving vortices
  !> on the sphere: A test case for horizontal advection problems", AMS 2008,
  !> equations (18)-(22)
  !> lat_pole and lon_pole denote the position of one of the vortices
  !> Parameter values have been taken from the paper and are currently
  !> hard-wired
  !>@param[in] lat Latitude of point
  !>@param[in] long Longitude of point
  !>@param[in] radius Distance from the centre of the planet of the point
  !>@param[in] time Time to compute the vortex at
  !>return pressure Value of the tracer field at this point
  function vortex_field(lat, long, radius, time) result(pressure)
    implicit none
    real(kind=r_def), intent(in) :: lat
    real(kind=r_def), intent(in) :: long
    real(kind=r_def), intent(in) :: radius
    real(kind=r_def), intent(in) :: time
    real(kind=r_def)             :: pressure

    real(kind=r_def) :: lat_pole, lon_pole
    real(kind=r_def) :: lat_dash, lon_dash
    real(kind=r_def) :: r0, v0, V
    real(kind=r_def) :: gamma_value, radial_distance, omega

    lat_pole = 0.0_r_def
    lon_pole = 0.0_r_def
    r0 = 3.0_r_def
    v0 = 38.61073731_r_def
    gamma_value = 5.0_r_def

    lat_dash = asin(sin(lat)*sin(lat_pole) &
                  + cos(lat)*cos(lat_pole)*cos(long-lon_pole))
    lon_dash = atan2(cos(lat)*sin(long-lon_pole),              &
                     cos(lat)*sin(lat_pole)*cos(long-lon_pole) &
                   - cos(lat_pole)*sin(lat) )

    radial_distance = r0*cos(lat_dash)
    V = v0*3.0_r_def*(sqrt(3.0_r_def)/2.0_r_def)*  &
                                (sinh(radial_distance)/cosh(radial_distance)**3)

    if (abs(radial_distance)<1E-10_r_def) then
      omega = 0.0_r_def
    else
      omega = V/(radius*radial_distance)
    end if

    pressure = 1.0_r_def - tanh((radial_distance/gamma_value)*     &
                                                      sin(lon_dash-omega*time))

  end function vortex_field

  !> @brief Compute an analytic exner pressure field
  !> @param[in] chi Position in physical coordinates
  !> @param[in] choice Integer defining which specified formula to use
  !> @result pressure The result pressure field
  function analytic_pressure(chi, choice, time) result(pressure)

    implicit none
    real(kind=r_def), intent(in) :: chi(3)
    integer,          intent(in) :: choice
    real(kind=r_def), intent(in) :: time
    real(kind=r_def)             :: pressure

    real(kind=r_def)             :: l, dt
    real(kind=r_def), parameter  :: XC = 0.0_r_def, &
                                    XR = 4000.0_r_def, &
                                    ZC_cold = 3000.0_r_def, &
                                    ZR = 2000.0_r_def
    real(kind=r_def)             :: long, lat, radius
    real(kind=r_def)             :: l1, l2
    real(kind=r_def)             :: h1, h2
    real(kind=r_def)             :: density, temperature
    real(kind=r_def)             :: t0
    real(kind=r_def)             :: u, v, w

    integer                      :: id

    if ( geometry == geometry_spherical ) then
      call xyz2llr(chi(1), chi(2), chi(3), long, lat, radius)
      call central_angle(long, lat, x1, y1, l1)
      call central_angle(long, lat, x2, y2, l2)
    else
      long = chi(1)
      lat  = chi(2)
      l1 = sqrt((long-x1)**2 + (lat-y1)**2)
      l2 = sqrt((long-x2)**2 + (lat-y2)**2)
    end if

    select case( choice )
    case (test_gravity_wave, test_isentropic, &
          test_isot_atm, test_isot_cold_atm,  &
          test_const_lapse_rate, test_specified_profiles)
      call reference_profile(pressure, density, temperature, chi, choice)

    case (test_cold_bubble_x, test_cold_bubble_y )
      if (choice == test_cold_bubble_x ) then
        id = 1
      else
        id = 2
      end if
      call reference_profile(pressure, density, temperature, chi, choice)
      l = sqrt( ((chi(id)-XC)/XR)**2 + ((chi(3)-ZC_cold)/ZR)**2 )
      if ( l <= 1.0_r_def ) then
        dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
        temperature = temperature - dt/pressure
        density = p_zero/(Rd*temperature) * pressure**( (1.0_r_def - kappa )/ kappa )
      end if

    !> No perturbation needed for warm bubble tests so just use background
    !> (isentropic) value
    case (test_warm_bubble, test_warm_bubble_3d, &
          test_bryan_fritsch, test_grabowski_clark )
      call reference_profile(pressure, density, temperature, chi, choice)

    case (test_gaussian_hill)
      h1 = tracer_max*exp( -(l1/r1)**2 )
      h2 = tracer_max*exp( -(l2/r2)**2 )
      pressure = h1 + h2

    case (test_cosine_hill)
      if ( l1 < r1 ) then
        h1 = tracer_background + (tracer_max/2.0_r_def)*(1.0_r_def+cos((l1/r1)*PI))
      else
        h1 = tracer_background
      end if
      if (l2 < r2) then
        h2 = tracer_background + (tracer_max/2.0_r_def)*(1.0_r_def+cos((l2/r2)*PI))
      else
        h2 = tracer_background
      end if
      pressure = h1+h2

    case (test_slotted_cylinder)
      ! Cylinder 1
      if ( l1 < r1 ) then
        if (abs(long-x1) > r1/6.0_r_def) then
          h1 = tracer_max
        else
          if (lat < y1-r1*5.0_r_def/12.0_r_def) then
            h1 = tracer_max
          else
            h1 = tracer_background
          end if
        end if
      else
        h1 = tracer_background
      end if
      ! Cylinder 2
      if ( l2 < r2 ) then
        if (abs(long-x2) > r2/6.0_r_def) then
          h2 = tracer_max
        else
          if (lat > y2+r2*5.0_r_def/12.0_r_def) then
            h2 = tracer_max
          else
            h2 = tracer_background
          end if
        end if
      else
        h2 = tracer_background
      end if
      pressure = h1 + h2

    case (test_constant_field)
      pressure = tracer_background

    case (test_cosine_stripe)
      l1 = sqrt((long-x1)**2)
      if ( l1 < r1 ) then
        pressure = tracer_background + (tracer_max/2.0_r_def)*(1.0_r_def+cos((l1/r1)*PI))
      else
        pressure = tracer_background
      end if

    case (test_vortex_field)
      pressure = vortex_field(lat,long,radius,time)

    case( test_yz_cosine_hill,       &
          test_cosine_bell,          &
          test_hadley_like_dcmip,    &
          test_eternal_fountain,     &
          test_curl_free_reversible, &
          test_div_free_reversible,  &
          test_rotational,           &
          test_translational,        &
          test_vertical_cylinder )
      ! This case is for transport of tracers and density only so it is not
      ! calculated here.

    case (test_solid_body_rotation, &
          test_solid_body_rotation_alt)
      t0 = 280.0_r_def
      temperature = analytic_temperature(chi, choice)
      pressure = t0/temperature
      density = p_zero/(Rd*temperature) * pressure**( (1.0_r_def - kappa )/ kappa )

    case (test_deep_baroclinic_wave)
      call deep_baroclinic_wave(long, lat, radius-scaled_radius, &
                                pressure, temperature, density, &
                                u, v, w)
    case(test_dry_cbl, test_shallow_conv, test_snow)
      call reference_profile(pressure, density, temperature, chi, choice)

    case( test_cos_phi )
      pressure = tracer_max*cos(lat)**4

    case( test_cosine_bubble )
      l1 = sqrt( ((chi(1) - x1)/r1)**2 + ((chi(3) - y1)/r2)**2 )
      if ( l1 < 1.0_r_def ) then
        pressure = tracer_background + tracer_max*cos(0.5_r_def*l1*PI)**2
      else
        pressure = tracer_background
      end if
    case default
      write( log_scratch_space, '(A)' )  'Invalid pressure profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

  end function analytic_pressure

end module analytic_pressure_profiles_mod
