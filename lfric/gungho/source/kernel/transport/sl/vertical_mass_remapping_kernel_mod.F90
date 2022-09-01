!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a vertical mass remapping of tracer-mass at W3.
!! @details The SLICE scheme (using different order reconstructions) is
!!          used to perform a vertical mass remapping of tracer mass.

module vertical_mass_remapping_kernel_mod

use argument_mod,         only : arg_type,              &
                                 GH_FIELD, GH_SCALAR,   &
                                 GH_REAL, GH_INTEGER,   &
                                 GH_READWRITE, GH_READ, &
                                 CELL_COLUMN
use fs_continuity_mod,    only : W2, W3
use constants_mod,        only : r_def, i_def
use kernel_mod,           only : kernel_type
! TODO #3011: these config options should be passed through as arguments
use transport_config_mod, only : slice_order_constant, slice_order_linear, &
                                 slice_order_parabola, slice_order_cubic

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed
!>                                      by the PSy layer.
type, public, extends(kernel_type) :: vertical_mass_remapping_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                     &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2), & ! departure points
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3), & ! mass
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)           &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_mass_remapping_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public ::  vertical_mass_remapping_code
private piecewise_reconstruction
private remapp_mass
private local_point_1d_array

contains

  !-------------------------------------------------------------------------------
  !> @details This kernel remaps mass from one grid to another, for
  !!          vertical motion only, then interpolates theta at the
  !!          departure point using 1D cubic-Lagrange interpolation.
  !> @param[in]     nlayers      The number of layers
  !> @param[in]     dep_pts_z    The vertical departure distance used for SL advection
  !> @param[in,out] mass         The mass field that needs to be remapped
  !> @param[in]     order        Order of the reconstruction of underlying function
  !> @param[in]     ndf_w2       The number of degrees of freedom per cell
  !!                             on W2 space
  !> @param[in]     undf_w2      The number of unique degrees of freedom
  !!                             on W2 space
  !> @param[in]     map_w2       The dofmap for the cell at the base of the column
  !!                             on W2 space
  !> @param[in]     ndf_w3       The number of degrees of freedom per cell
  !!                             on w3 space
  !> @param[in]     undf_w3      The number of unique degrees of freedom
  !!                             on w3 space
  !> @param[in]     map_w3       The dofmap for the cell at the base of the column
  !!                             on w3 space
  !-------------------------------------------------------------------------------

  subroutine vertical_mass_remapping_code( nlayers,                     &
                                           dep_pts_z,                   &
                                           mass,                        &
                                           order,                       &
                                           ndf_w2, undf_w2, map_w2,     &
                                           ndf_w3, undf_w3, map_w3      )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                         :: nlayers
  integer(kind=i_def), intent(in)                         :: ndf_w2
  integer(kind=i_def), intent(in)                         :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)      :: map_w2
  integer(kind=i_def), intent(in)                         :: ndf_w3
  integer(kind=i_def), intent(in)                         :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)      :: map_w3

  real(kind=r_def), dimension(undf_w2), intent(in)        :: dep_pts_z
  real(kind=r_def), dimension(undf_w3), intent(inout)     :: mass
  integer(kind=i_def), intent(in)                         :: order

  integer(kind=i_def) :: k, nz, nzl
  integer(kind=i_def) :: ext = 0

  real(kind=r_def), allocatable :: dist(:)
  real(kind=r_def), allocatable :: zl(:)
  real(kind=r_def), allocatable :: zd(:)
  real(kind=r_def), allocatable :: dz(:)
  real(kind=r_def), allocatable :: m0(:)
  real(kind=r_def), allocatable :: mn(:)
  real(kind=r_def), allocatable :: a0(:)
  real(kind=r_def), allocatable :: a1(:)
  real(kind=r_def), allocatable :: a2(:)
  real(kind=r_def), allocatable :: a3(:)

  nz = nlayers
  nzl = nz + 1_i_def
  allocate( dist(1-ext:nzl+ext))
  allocate( zl(1-ext:nzl+ext))
  allocate( zd(1-ext:nzl+ext))
  allocate( dz(1-ext:nz+ext))
  allocate( m0(1-ext:nz+ext))
  allocate( mn(1-ext:nz+ext))
  allocate( a0(1-ext:nz+ext))
  allocate( a1(1-ext:nz+ext))
  allocate( a2(1-ext:nz+ext))
  allocate( a3(1-ext:nz+ext))

  ! Extract and fill local column from global data
  do k=0,nlayers
    dist(k+1) = dep_pts_z(map_w2(5)+k)
  end do
  do k=0,nlayers - 1
    m0(k+1) = mass(map_w3(1)+k)
  end do

  ! Create a local 1-d vertical remapping problem
  if ( ext > 0 ) then
    do k=1-ext, 0
      m0(k) = m0(1)
    end do
    do k=nz+1,nz+ext
      m0(k) = m0(nz)
    end do
  end if

  do k=1-ext, nzl+ext
    zl(k) = real(k-1,r_def)
  end do
  do k=1-ext, nzl+ext
    zd(k) = zl(k) - dist(k)
    zd(k) = min( max(zl(1-ext),zd(k)), zl(nzl+ext) )
  end do
  do k=1-ext,nz+ext
    dz(k)  = 1.0_r_def
  end do

  call piecewise_reconstruction(zl,dz,m0,a0,a1,a2,a3,order,1-ext,nz+ext)
  call remapp_mass(zd,zl,dz,m0,a0,a1,a2,a3,1-ext,nz+ext,mn )

  ! Export the local new mass back to the global mass

  do k=0,nlayers - 1
    mass(map_w3(1)+k) = mn(k+1)
  end do

  deallocate(zl,zd,dz,m0,mn,a0,a1,a2,a3,dist)

  end subroutine vertical_mass_remapping_code

  !-------------------------------------------------------------------------------
  !> @brief   This routine returns the cubic coefficents (a0,a1,a2,a3) for each cell.
  !> @details These coefficents define a piece-wise function for each cell. This
  !>          function is defined such that the integral of the function over the
  !>          cell gives the mass of that cell.
  !> @param[in]  xl     The grid-edges of the given mass (mass(i) is the mass of
  !!                    cell [xl(i),xl(i+1)])
  !> @param[in]  dx     Size of the cells, dx(i)=xl(i+1)-xl(i)
  !> @param[in]  mass   The mass of grid cells
  !> @param[out] a0     The first coefficent of the cell reconstruction.
  !>                    Each cell has a cubic function f(x) = a0 + a1*x + a2*x^2 + a3*x^3
  !> @param[out] a1     The second coefficent of the cell reconstruction.
  !> @param[out] a2     The third  coefficent of the cell reconstruction.
  !> @param[out] a3     The fourth coefficent of the cell reconstruction.
  !> @param[in] order   Order the piecewise function (constant, linear, quadratic, or cubic)
  !> @param[in] ns      Index of the starting cell
  !> @param[in] nf      Index of the end cell
  !-------------------------------------------------------------------------------
  subroutine piecewise_reconstruction(xl,dx,mass,a0,a1,a2,a3,order,ns,nf)

   implicit none

   integer(kind=i_def), intent(in)                  :: ns, nf, order
   real(kind=r_def), dimension(ns:nf+1), intent(in) :: xl
   real(kind=r_def), dimension(ns:nf), intent(in)   :: dx, mass
   real(kind=r_def), dimension(ns:nf), intent(out)  :: a0, a1, a2, a3

   ! Local variables
   real(kind=r_def), dimension(ns:nf) :: rho_left_cv, rho_right_cv, slope_rho
   real(kind=r_def), dimension(ns:nf) :: y_centre_cv, rho_bar
   integer(kind=i_def), parameter :: min_cvs_required = 4_i_def
   real(kind=r_def), dimension(ns:nf+1) :: rho_left_cv_all,slope_im,slope_ip
   integer(kind=i_def) :: j, j0, j1, cv_start
   real(kind=r_def)    :: y0, y1, y2, y3, y4, y, m1, m2, m3, m4
   real(kind=r_def), parameter :: c2=2.0_r_def, c3=3.0_r_def, c4=4.0_r_def,  &
                                  c6=6.0_r_def, c9=9.0_r_def, c12=12.0_r_def

    j0 = int(min_cvs_required/2,i_def)
    j1 = min_cvs_required - 1_i_def

    do j = ns, nf
      y_centre_cv(j) = 0.5_r_def * ( xl(j) + xl(j+1) )
      a1(j) = 0.0_r_def
      a2(j) = 0.0_r_def
      a3(j) = 0.0_r_def
      rho_bar(j) = mass(j) / dx(j)
    end do

   ! Use piecewise-constant for less than 4 cells
    if ( ((nf - ns + 1 ) < min_cvs_required) .or.    &
                   (order == slice_order_constant) ) then
      do j = ns, nf
        a0(j) = rho_bar(j)
      end do

    else

      do j = ns, nf + 1
        cv_start = min( max( j - j0, ns ) + j1, nf ) - j1
        y0 = xl(cv_start)
        y1 = xl(cv_start+1) - y0
        y2 = xl(cv_start+2) - y0
        y3 = xl(cv_start+3) - y0
        y4 = xl(cv_start+4) - y0
        m1 =   mass(cv_start)
        m2 =   mass(cv_start+1) + m1
        m3 =   mass(cv_start+2) + m2
        m4 =   mass(cv_start+3) + m3
        m1 = m1 / ( y1*(y1-y2)*(y1-y3)*(y1-y4) )
        m2 = m2 / ( y2*(y2-y1)*(y2-y3)*(y2-y4) )
        m3 = m3 / ( y3*(y3-y1)*(y3-y2)*(y3-y4) )
        m4 = m4 / ( y4*(y4-y1)*(y4-y2)*(y4-y3) )

        y  = xl(j) - y0
        rho_left_cv_all(j) =                                                             &
             ( y*( y*(c4*y - c3*(y2+y3+y4)) + c2*(y3*y2+y2*y4+y3*y4) ) - y3*y2*y4 )*m1 + &
             ( y*( y*(c4*y - c3*(y1+y3+y4)) + c2*(y1*y3+y1*y4+y3*y4) ) - y1*y3*y4 )*m2 + &
             ( y*( y*(c4*y - c3*(y2+y4+y1)) + c2*(y2*y4+y1*y2+y1*y4) ) - y1*y2*y4 )*m3 + &
                  (y*(y*(c4*y-c3*(y2+y3+y1))+c2*(y3*y2+y1*y2+y1*y3))-y1*y2*y3)*m4
        y  = y_centre_cv( max( ns, j - 1 ) ) - y0
        slope_im(j) =                                                           &
                  ( y*(c12*y - c6*(y2+y3+y4)) + c2*(y3*y2+y2*y4+y3*y4) )*m1 +   &
                  ( y*(c12*y - c6*(y1+y3+y4)) + c2*(y1*y3+y1*y4+y3*y4) )*m2 +   &
                  ( y*(c12*y - c6*(y2+y4+y1)) + c2*(y2*y4+y1*y2+y1*y4) )*m3 +   &
                  ( y*(c12*y - c6*(y2+y3+y1)) + c2*(y3*y2+y1*y2+y1*y3) )*m4
        y  = y_centre_cv( min( j, nf    ) ) - y0
        slope_ip(j) =                                                           &
                  ( y*(c12*y - c6*(y2+y3+y4)) + c2*(y3*y2+y2*y4+y3*y4) )*m1 +   &
                  ( y*(c12*y - c6*(y1+y3+y4)) + c2*(y1*y3+y1*y4+y3*y4) )*m2 +   &
                  ( y*(c12*y - c6*(y2+y4+y1)) + c2*(y2*y4+y1*y2+y1*y4) )*m3 +   &
                  ( y*(c12*y - c6*(y2+y3+y1)) + c2*(y3*y2+y1*y2+y1*y3) )*m4
      end do

      do j = ns, nf
        slope_rho(j)   = 0.5_r_def * ( slope_im(j+1) + slope_ip(j) )
        slope_rho(j)   = slope_rho(j) * dx(j)
        rho_left_cv(j) = rho_left_cv_all(j)
        rho_right_cv(j)= rho_left_cv_all(j+1)
      end do

      if ( order == slice_order_linear ) then
         do j = ns, nf
           if ( (rho_bar(j)-rho_left_cv(j)) > (rho_bar(j)-rho_right_cv(j)) ) then
             a0(j) = rho_left_cv(j)
             a1(j) = c2*(rho_bar(j)-rho_left_cv(j))
           else
             a0(j) = c2*rho_bar(j)-rho_right_cv(j)
             a1(j) = c2*(rho_right_cv(j)-rho_bar(j))
           end if
         end do
      else if ( order == slice_order_parabola ) then
         do j = ns, nf
            a0(j) = rho_left_cv(j)
            a1(j) = -c4*rho_left_cv(j) - c2*rho_right_cv(j) + c6*rho_bar(j)
            a2(j) =  c3*rho_left_cv(j) + c3*rho_right_cv(j) - c6*rho_bar(j)
         end do
      else if ( order == slice_order_cubic ) then
         do j = ns, nf
            a0(j) = rho_left_cv(j)
            a1(j) = -c6*rho_left_cv(j) + c6*rho_bar(j)   - c2*slope_rho(j)
            a2(j) =  c9*rho_left_cv(j) - c3*rho_right_cv(j)  &
                     - c6*rho_bar(j)     + c6*slope_rho(j)
            a3(j) = -c4*rho_left_cv(j) + c4*rho_right_cv(j) - c4*slope_rho(j)
         end do
      end if
    end if

  end subroutine piecewise_reconstruction

  !-------------------------------------------------------------------------------
  !> @brief   This routine remaps conservatively the mass from a grid xl to a
  !>          superposed grid xld.
  !> @details The mass(i) is the mass of cell [xl(i),xl(i+1)] and
  !>           mass_d(i) is the mass of the cell [xld(i),xld(i+1)]
  !>           if xl and xld are bounded by the same boundaries then
  !>           sum{mass} = sum{mass_d}
  !> @param[in] xld  Departure points of grid xl
  !> @param[in] xl   The grid-edges of the given mass (mass(i) is the mass of
  !!                 cell [xl(i),xl(i+1)])
  !> @param[in] dx   Size of the cells, dx(i)=xl(i+1)-xl(i)
  !> @param[in] mass The mass of grid cells
  !> @param[in] a0   The first coefficent of the cell reconstruction.
  !>                 Each cell has a cubic function f(x) = a0 + a1*x + a2*x^2 + a3*x^3
  !> @param[in] a1   The second coefficent of the cell reconstruction.
  !> @param[in] a2   The third  coefficent of the cell reconstruction.
  !> @param[in] a3   The fourth coefficent of the cell reconstruction.
  !> @param[in] ns   Index of the starting cell
  !> @param[in] nf   Index of the end cell
  !> @param[out] mass_d  The remapped mass (mass_d(i) is the mass of cell [xld(i),xld(i+1)])
  !-------------------------------------------------------------------------------
  subroutine remapp_mass( xld, xl, dx, mass, a0, a1, a2, a3, ns, nf, mass_d )

    implicit none

    integer(kind=i_def), intent(in)                  :: ns,nf
    real(kind=r_def), dimension(ns:nf+1), intent(in) :: xl, xld
    real(kind=r_def), dimension(ns:nf), intent(in)   :: dx, mass, a0,a1,a2,a3
    real(kind=r_def), dimension(ns:nf), intent(out)  :: mass_d

    ! Local variables
    real(kind=r_def)    :: x1,x2,m1,m2,m3,h1,h2,s1,s2
    integer(kind=i_def) :: j, i, cv1, cv2
    real(kind=r_def)    :: c1, c2, c3, sig

    c1 = 0.5_r_def
    c2 = 1.0_r_def/3.0_r_def
    c3 = 0.25_r_def

    do j = ns, nf
       m1 = 0.0_r_def
       m2 = 0.0_r_def
       m3 = 0.0_r_def
       x1 = min(xld(j),xld(j+1))
       x2 = max(xld(j),xld(j+1))
       sig = sign(1.0_r_def, xld(j+1)-xld(j) )
       cv1 = local_point_1d_array(x1, xl, ns, nf)
       cv2 = local_point_1d_array(x2, xl, ns, nf)
       h1 = 1.0_r_def/dx(cv1)
       h2 = 1.0_r_def/dx(cv2)

       s1 = ( x1 - xl(cv1) ) * h1
       s2 = 1.0_r_def

       m1 = a0(cv1)*(s2-s1) + c1*a1(cv1)*(s2**2 - s1**2) + &
                              c2*a2(cv1)*(s2**3 - s1**3) + &
                              c3*a3(cv1)*(s2**4 - s1**4)
       m1 = m1 * dx(cv1)

       s1 = 0.0_r_def
       s2 = ( x2 - xl(cv2) ) * h2

       m2 = a0(cv2)*(s2-s1) + c1*a1(cv2)*(s2**2 - s1**2) + &
                              c2*a2(cv2)*(s2**3 - s1**3) + &
                              c3*a3(cv2)*(s2**4 - s1**4)
       m2 = m2 * dx(cv2)

       if ( (cv2 - cv1) == 0_i_def ) then
          m3 = -dx(cv1)*(a0(cv1) + c1*a1(cv1) + c2*a2(cv1) + c3*a3(cv1))
        else
          do i = cv1+1, cv2-1
            s1 = a0(i) + c1*a1(i) + c2*a2(i) + c3*a3(i)
            m3 = m3 + s1*dx(i)
          end do
       end if
       mass_d(j) = (m1 + m2 + m3)*sig
    end do

  end subroutine remapp_mass

  !-------------------------------------------------------------------------------
  !> @brief Compute location index of the point p in a 1d array x.
  !> @param[in] p Coordinate of the point
  !> @param[in] x The 1d array of grid-points
  !> @param[in] ns Start of index of grid-points
  !> @param[in] nf End of index of grid-points
  !> @return location index "l" such as x(l) < p < x(l+1)
  !-------------------------------------------------------------------------------
  function local_point_1d_array(p, x, ns, nf)

  implicit none

  integer(kind=i_def), intent(in) :: ns, nf
  real(kind=r_def), dimension(ns:nf), intent(in) :: x
  real(kind=r_def),    intent(in) :: p
  integer(kind=i_def)             :: jl,jm,ju
  integer(kind=i_def)             :: local_point_1d_array

  jl=ns-1
  ju=nf+1
  do while ((ju-jl) > 1 )
    jm=(ju+jl)/2
    if( (x(nf) > x(ns)) .eqv. (p > x(jm)) ) then
      jl=jm
    else
      ju=jm
    end if
  end do
  local_point_1d_array = max(ns,min(jl,nf))

  end function local_point_1d_array

end module vertical_mass_remapping_kernel_mod
