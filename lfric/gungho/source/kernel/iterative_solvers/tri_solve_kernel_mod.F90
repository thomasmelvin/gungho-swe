!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Tridiagonal solver using Thomas algorithm.
!>
!> @details Solve the triadiagonal system of equations
!>          \f[ A^+ y^+ + A^0 y + A^- y^- = x \]f
!>          for known x and matrix A using the Thomas algorithm
!>          Code is only valid for lowest order elements
!>
module tri_solve_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_READ, GH_WRITE, &
                               CELL_COLUMN
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: tri_solve_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/               &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3), &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3), &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  W3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tri_solve_code
  end type tri_solve_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tri_solve_code

contains

!> @brief Tridiagonal solver using Thomas algorithm
!> @param[in]  nlayers Number of levels to solve over
!> @param[in,out] y LHS field to solve for
!> @param[in]  x RHS field
!> @param[in]  tri_0 Centred part of tridiagonal matrix
!> @param[in]  tri_plus Upper diagonal part of tridiagonal matrix
!> @param[in]  tri_minus lower diagonal part of tridiagonal matrix
!> @param[in]  ndf Number of dofs per cell for all fields, should be = 1
!> @param[in]  undf Size of all field arrays
!> @param[in]  map Array containing the address of the first dof in the column
subroutine tri_solve_code(nlayers,                    &
                          y, x,                       &
                          tri_0, tri_plus, tri_minus, &
                          ndf, undf, map)

  implicit none

  integer(kind=i_def), intent(in) :: ndf, undf
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), dimension(ndf),  intent(in) :: map

  real(kind=r_def), dimension(undf), intent(inout) :: y
  real(kind=r_def), dimension(undf), intent(in)    :: x
  real(kind=r_def), dimension(undf), intent(in)    :: tri_0
  real(kind=r_def), dimension(undf), intent(in)    :: tri_plus
  real(kind=r_def), dimension(undf), intent(in)    :: tri_minus

  integer(kind=i_def)                  :: k, ij
  real(kind=r_def), dimension(nlayers) :: x_new, tri_plus_new
  real(kind=r_def)                     :: denom

  k  = 0
  ij = map(1)
  denom = 1.0_r_def/tri_0(ij+k)
  tri_plus_new(1) = tri_plus(ij+k)*denom
  x_new(1)        = x(ij+k)       *denom

  do k = 1,nlayers-1
    denom = 1.0_r_def/(tri_0(ij+k) - tri_minus(ij+k)*tri_plus_new(k))
    tri_plus_new(k+1) = tri_plus(ij+k)*denom
    x_new(k+1)        = (x(ij+k) - tri_minus(ij+k)*x_new(k))*denom
  end do

  k = nlayers-1
  y(ij+k) = x_new(k+1)
  do k = nlayers-2,0,-1
    y(ij+k) = x_new(k+1) - tri_plus_new(k+1)*y(ij+k+1)
  end do

end subroutine tri_solve_code

end module tri_solve_kernel_mod

