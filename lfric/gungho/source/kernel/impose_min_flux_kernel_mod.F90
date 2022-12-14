!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module impose_min_flux_kernel_mod
use argument_mod,            only : arg_type,                           &
                                    GH_FIELD, GH_OPERATOR, GH_READ,     &
                                    CELL_COLUMN, GH_REAL, GH_INC,       &
                                    GH_SCALAR
use fs_continuity_mod,       only : W3, W2
use constants_mod,           only : r_def, i_def, EPS
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: impose_min_flux_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                    &
       arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),     &
       arg_type(GH_FIELD,    GH_REAL, GH_INC,   W2),     &
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W2), &
       arg_type(GH_SCALAR,   GH_REAL, GH_READ),          &
       arg_type(GH_SCALAR,   GH_REAL, GH_READ)           &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: impose_min_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public impose_min_flux_code

contains

!> @brief Modify the flux so the updated field f(n+1) is greater than field_min
!!        i.e., f(n+1) = f(n) - dts*div(flux) >= field_min
!! @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!! @param[in] field Field at level time n
!! @param[inout] flux Input/output
!! @param[in] ncell_3d1 Total number of cells related to div
!! @param[in] div Divergence operator used in the update
!! @param[in] field_min The minimum value we want to enfore for the updated field
!! @param[in] dts  The time-step used in the update
!! @param[in] ndf1 Number of degrees of freedom per cell for the output field
!! @param[in] undf1 Unique number of degrees of freedom  for the output field
!! @param[in] map1 Dofmap for the cell at the base of the column for the output field
!! @param[in] ndf2 Number of degrees of freedom per cell for the input field
!! @param[in] undf2 Unique number of degrees of freedom for the input field
!! @param[in] map2 Dofmap for the cell at the base of the column for the input field
subroutine impose_min_flux_code(cell,              &
                                nlayers,           &
                                field, flux,       &
                                ncell_3d1,         &
                                div,               &
                                field_min,         &
                                dts,               &
                                ndf1, undf1, map1, &
                                ndf2, undf2, map2  )

  implicit none

  ! Arguments
  integer(kind=i_def),                  intent(in) :: cell, nlayers
  integer(kind=i_def),                  intent(in) :: ncell_3d1
  integer(kind=i_def),                  intent(in) :: undf1, ndf1
  integer(kind=i_def),                  intent(in) :: undf2, ndf2
  integer(kind=i_def), dimension(ndf1), intent(in) :: map1
  integer(kind=i_def), dimension(ndf2), intent(in) :: map2
  real(kind=r_def), dimension(undf2),               intent(inout) :: flux
  real(kind=r_def), dimension(undf1),               intent(in)    :: field
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d1), intent(in)    :: div
  real(kind=r_def), intent(in)                                    :: field_min
  real(kind=r_def), intent(in)                                    :: dts

  ! Internal variables
  integer(kind=i_def)                  :: k, ik, df1, df2
  real(kind=r_def), dimension(ndf2)    :: cell_fluxes
  real(kind=r_def)                     :: a, b, inc, inc_n, flux_scaler
  integer(kind=i_def), dimension(ndf2) :: flux_change_id

  do k = 0, nlayers-1

    do df2 = 1, ndf2
      cell_fluxes(df2) = flux(map2(df2)+k)
    end do
    ik = (cell-1)*nlayers + k + 1

    do df1 = 1, ndf1

       inc_n = 0.0_r_def
       flux_change_id = 0_i_def

       do df2 = 1, ndf2
         inc = - dts*div(df1,df2,ik)*cell_fluxes(df2)
         if ( inc < 0.0_r_def ) then
             inc_n = inc_n - inc
             flux_change_id(df2) = 1_i_def
         end if
       end do

       a = field(map1(df1)+k) - (field_min + EPS)
       b = a / max(inc_n, EPS)
       flux_scaler = min(max(0.0_r_def,b),1.0_r_def)

       do df2 = 1, ndf2
          if ( flux_change_id(df2) == 1_i_def ) then
              flux(map2(df2)+k) = flux(map2(df2)+k) * flux_scaler
          end if
       end do

    end do

  end do

end subroutine impose_min_flux_code

end module impose_min_flux_kernel_mod
