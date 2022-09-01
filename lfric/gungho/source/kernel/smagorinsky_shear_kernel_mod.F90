!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculate 3D rate of strain tensor at theta points needed for
!>        Smagorinsky diffusion coefficient.
!>
module smagorinsky_shear_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_READ, GH_WRITE,         &
                                STENCIL, CROSS, CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: smagorinsky_shear_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W2,     STENCIL(CROSS)),   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W2),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3)                        &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: smagorinsky_shear_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: smagorinsky_shear_code

contains

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Number of layers in the mesh
!! @param[in,out] shear 3D wind shear
!! @param[in] u_n Input wind field
!! @param[in] map_w2_stencil_size Number of cells in the stencil at the base
!!                                of the column for W2
!! @param[in] map_w2_stencil Array holding the stencil dofmap for the cell at
!!                           the base of the column for W2
!! @param[in] dx_at_w2 Grid length at the w2 dofs
!! @param[in] height_wth Height of wth space levels above surface
!! @param[in] height_w3 Height of w3 space levels above surface
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom for theta space
!! @param[in] map_wt Array holding the dofmap for the cell at the base of
!!                   the column for theta space
!! @param[in] ndf_w2 Number of degrees of freedom per cell for wind space
!! @param[in] undf_w2  Number of unique degrees of freedom for wind space
!! @param[in] map_w2 Array holding the dofmap for the cell at the base of
!!                   the column for wind space
!! @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!! @param[in] undf_w3  Number of unique degrees of freedom for W3
!! @param[in] map_w3 Array holding the dofmap for the cell at the base of
!!                   the column for W3
subroutine smagorinsky_shear_code( nlayers,                                 &
                                   shear,                                   &
                                   u_n,                                     &
                                   map_w2_stencil_size, map_w2_stencil,     &
                                   dx_at_w2,                                &
                                   height_wth,                              &
                                   height_w3,                               &
                                   ndf_wt, undf_wt, map_wt,                 &
                                   ndf_w2, undf_w2, map_w2,                 &
                                   ndf_w3, undf_w3, map_w3                  &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2, ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: map_w2_stencil_size
  integer(kind=i_def), dimension(ndf_w2,map_w2_stencil_size), intent(in)  :: map_w2_stencil
  integer(kind=i_def), dimension(ndf_w2),  intent(in)  :: map_w2
  integer(kind=i_def), dimension(ndf_w3),  intent(in)  :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in)  :: map_wt

  real(kind=r_def), dimension(undf_wt),  intent(inout) :: shear
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: u_n, dx_at_w2
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: height_wth
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: height_w3

  ! Internal variables
  integer(kind=i_def)                      :: k, km, kp, df
  real(kind=r_def)                         :: weight_pl_w3, weight_min_w3
  real(kind=r_def)                         :: weight_pl_wth, weight_min_wth
  real(kind=r_def)                         :: sum_sij, ssq12k, ssq11, ssq22, ssq33
  real(kind=r_def)                         :: ssq12up, ssq12, ssq13, ssq23
  real(kind=r_def), dimension(0:nlayers-1) :: idx2, idy2
  real(kind=r_def), dimension(0:nlayers-1,4) :: idx_w2
  real(kind=r_def), dimension(0:nlayers-1) :: dz_w3, idz_w3, idz_w3_2
  real(kind=r_def), dimension(1:nlayers-1) :: idz_wth
  real(kind=r_def), parameter              :: smallp=1.0e-14_r_def

  ! Assumed direction for derivatives in this kernel is:
  !  y
  !  ^
  !  |_> x
  !

  ! Layout of dofs for the stencil map
  ! dimensions of map are (ndf, ncell)
  ! Horizontally:
  !
  !   -- 4 --
  !   |     |
  !   1     3
  !   |     |
  !   -- 2 --
  !
  ! df = 5 is in the centre on the bottom face
  ! df = 6 is in the centre on the top face

  ! The layout of the cells in the stencil is:
  !
  !          -----
  !          |   |
  !          | 5 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !          |   |
  !          | 3 |
  !          -----
  !

  ! If the full stencil isn't available, we must be at the domain edge.
  ! Simply set the shear to 0 for now, and exit the routine.
  if (map_w2_stencil_size < 5_i_def) then
    do k = 0, nlayers
      shear(map_wt(1) + k) = 0.0_r_def
    end do
    return
  end if

  ! Compute horizontal grid spacings:
  ! N.B. not accounting for difference between w3 and wtheta heights
  do k = 0, nlayers - 1
    ! At cell centre, dx and dy are average of cell face values
    idx2(k) = (2.0_r_def/(dx_at_w2(map_w2(1)+k)+dx_at_w2(map_w2(3)+k)))**2
    idy2(k) = (2.0_r_def/(dx_at_w2(map_w2(2)+k)+dx_at_w2(map_w2(4)+k)))**2
    do df = 1,4
      ! Inverse dx at cell faces
      idx_w2(k,df) = 1.0_r_def/dx_at_w2(map_w2(df)+k)
    end do
  end do

  ! Compute vertical grid spacings:
  ! dz between adjacent wth/theta levels, held on w3/rho levels
  ! dz on 1st w3/rho level:
  k = 0
  kp = k + 1
  dz_w3(k) = height_wth(map_wt(1) + kp) - height_wth(map_wt(1) + k)
  idz_w3(k) = 1.0_r_def / dz_w3(k)
  idz_w3_2(k) = idz_w3(k)*idz_w3(k)

  do k = 1, nlayers - 1
    km = k - 1
    kp = k + 1

    ! dz between adjacent wth/theta levels, held on w3/rho levels
    dz_w3(k) = height_wth(map_wt(1) + kp) - height_wth(map_wt(1) + k)
    idz_w3(k) = 1.0_r_def / dz_w3(k)
    idz_w3_2(k) = idz_w3(k)*idz_w3(k)

    ! dz between adjacent w3/rho levels, held on wth/theta levels
    idz_wth(k) = 1.0_r_def / (height_w3(map_w3(1) + k) - height_w3(map_w3(1) + km))

  end do

  ! Calculate half-squared strain rate SSQ on wth points:
  ! shear(0) and shear(nlayers) are both initialised to zero and left unchanged

  k = 0
  ! ssq12: (du/dy + dv/dx)^2 on 1st w3 level
  ! To be averaged to wth/theta level later
  ssq12k = ( ( idx_w2(k,2) * (u_n(map_w2_stencil(1,1) + k) - u_n(map_w2_stencil(1,3) + k) ) +             &
            idx_w2(k,1) * (u_n(map_w2_stencil(2,1) + k) - u_n(map_w2_stencil(2,2) + k) ) )**2 +           &
            ( idx_w2(k,2) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(3,3) + k) ) +              &
            idx_w2(k,3) * (u_n(map_w2_stencil(2,4) + k) - u_n(map_w2_stencil(2,1) + k) ) )**2 +           &
            ( idx_w2(k,4) * (u_n(map_w2_stencil(1,5) + k) - u_n(map_w2_stencil(1,1) + k) ) +              &
            idx_w2(k,1) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(4,2) + k) ) )**2 +           &
            ( idx_w2(k,4) * (u_n(map_w2_stencil(3,5) + k) - u_n(map_w2_stencil(3,1) + k) ) +              &
            idx_w2(k,3) * (u_n(map_w2_stencil(4,4) + k) - u_n(map_w2_stencil(4,1) + k) ) )**2 ) / 4

  ! Set lowest level to 0
  shear(map_wt(1) + 0) = 0.0_r_def

  do k = 1, nlayers - 1
    km = k - 1
    kp = k + 1

    ! Vertical interpolation weights:
    ! Vertical interpolation of w3 to wth levels
    weight_pl_w3 = (height_wth(map_wt(1) + k) - height_w3(map_w3(1) + km)) /             &
                   (height_w3(map_w3(1) + k) - height_w3(map_w3(1) + km))
    weight_min_w3 = (height_w3(map_w3(1) + k) - height_wth(map_wt(1) + k)) /             &
                    (height_w3(map_w3(1) + k) - height_w3(map_w3(1) + km))

    ! Vertical interpolation of wth to wth levels
    weight_pl_wth = dz_w3(km) / (dz_w3(km) + dz_w3(k))
    weight_min_wth = dz_w3(k) / (dz_w3(km) + dz_w3(k))

    ! ssq11: 2 * backward difference (du/dx)^2 averaged to wth
    ssq11 = 2.0_r_def * (                                                                                     &
            weight_pl_w3 * idx2(k) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(1,1) + k) )**2 +      &
            weight_min_w3 * idx2(km) * (u_n(map_w2_stencil(3,1) + km) - u_n(map_w2_stencil(1,1) + km) )**2 )

    ! ssq22: 2 * backward difference (dv/dy)^2 averaged to wth
    ssq22 = 2.0_r_def * (                                                                                     &
            weight_pl_w3 * idy2(k) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(2,1) + k) )**2 +      &
            weight_min_w3 * idy2(km) * (u_n(map_w2_stencil(4,1) + km) - u_n(map_w2_stencil(2,1) + km) )**2 )

    ! ssq33: 2 * backward difference (dw/dz)^2 averaged to wth
    ssq33 = 2.0_r_def * (                                                                                         &
            weight_pl_wth * idz_w3_2(km) * (u_n(map_w2_stencil(6,1) + km) - u_n(map_w2_stencil(5,1) + km) )**2 +  &
            weight_min_wth * idz_w3_2(k) * (u_n(map_w2_stencil(6,1) + k) - u_n(map_w2_stencil(5,1) + k) )**2 )

    ! ssq13: (du/dz + dw/dx)^2 averaged to wth
    ! ssq31 = ssq13
    ssq13 = ( ( idz_wth(k) * (u_n(map_w2_stencil(1,1) + k) - u_n(map_w2_stencil(1,1) + km) ) +            &
            idx_w2(k,1) * (u_n(map_w2_stencil(5,1) + k) - u_n(map_w2_stencil(5,2) + k) ) )**2 +           &
            ( idz_wth(k) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(3,1) + km) ) +              &
            idx_w2(k,3) * (u_n(map_w2_stencil(5,4) + k) - u_n(map_w2_stencil(5,1) + k) ) )**2 ) / 2

    ! ssq23: (dw/dy + dv/dz)^2 averaged to wth
    ! ssq32 = ssq23
    ssq23 = ( ( idx_w2(k,4) * (u_n(map_w2_stencil(5,5) + k) - u_n(map_w2_stencil(5,1) + k) ) +            &
            idz_wth(k) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(4,1) + km) ) )**2 +           &
            ( idx_w2(k,2) * (u_n(map_w2_stencil(5,1) + k) - u_n(map_w2_stencil(5,3) + k) ) +              &
            idz_wth(k) * (u_n(map_w2_stencil(2,1) + k) - u_n(map_w2_stencil(2,1) + km) ) )**2 ) / 2

    ! ssq12: (du/dy + dv/dx)^2 on w3 level k
    ! ssq21 = ssq12
    ssq12up = ( ( idx_w2(k,2) * (u_n(map_w2_stencil(1,1) + k) - u_n(map_w2_stencil(1,3) + k) ) +          &
              idx_w2(k,1) * (u_n(map_w2_stencil(2,1) + k) - u_n(map_w2_stencil(2,2) + k) ) )**2 +         &
              ( idx_w2(k,2) * (u_n(map_w2_stencil(3,1) + k) - u_n(map_w2_stencil(3,3) + k) ) +            &
              idx_w2(k,3) * (u_n(map_w2_stencil(2,4) + k) - u_n(map_w2_stencil(2,1) + k) ) )**2 +         &
              ( idx_w2(k,4) * (u_n(map_w2_stencil(1,5) + k) - u_n(map_w2_stencil(1,1) + k) ) +            &
              idx_w2(k,1) * (u_n(map_w2_stencil(4,1) + k) - u_n(map_w2_stencil(4,2) + k) ) )**2 +         &
              ( idx_w2(k,4) * (u_n(map_w2_stencil(3,5) + k) - u_n(map_w2_stencil(3,1) + k) ) +            &
              idx_w2(k,3) * (u_n(map_w2_stencil(4,4) + k) - u_n(map_w2_stencil(4,1) + k) ) )**2 ) / 4

    ! average ssq21 to wth level k
    ssq12 = weight_pl_w3 * ssq12up + weight_min_w3 * ssq12k
    ! save the ssq12 value at this level for use at the next iteration of k
    ssq12k = ssq12up

    ! sum_ssqij would be 2*ssq11 + 2*ssq22 + 2*ssq33 + ssq13 +
    !                      ssq23 + ssq12 + ssq31 + ssq32 + ssq21
    ! hence by equivalence of the off-diagonal terms (ssq12 etc),
    ! ssq = 0.5 (sum_ssqij) is given by sum_sij below
    !
    ! the extra factor of 2 in the diagnonal terms comes from the fact
    ! that they should be (2*du/dx)^2, not 2*(du/dx)^2

    sum_sij = ssq11 + ssq22 + ssq33 + ssq13 + ssq23 + ssq12 + smallp

    shear(map_wt(1) + k) = SQRT(sum_sij)

  end do

  ! Set model top value
  shear(map_wt(1) + nlayers) = 0.0_r_def

end subroutine smagorinsky_shear_code

end module smagorinsky_shear_kernel_mod
