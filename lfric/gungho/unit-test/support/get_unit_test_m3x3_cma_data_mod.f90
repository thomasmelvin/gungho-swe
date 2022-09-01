!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_m3x3_cma_data_mod
! A module containing a collection of helper routines that provide the sizes of
! canned  data for CMA operators on a simple 3x3 mesh

  use constants_mod, only : i_def, r_def

  implicit none

  private
  public :: get_wtheta_w3_m3x3_cma_data
  public :: get_wtheta_m3x3_cma_data
  public :: get_cma_size
  public :: get_cma_prod_size
  public :: get_w3_w3_m3x3_cma_data
  public :: get_w2v_w2v_m3x3_cma_data
contains



  subroutine get_wtheta_w3_m3x3_cma_data(ndf_to, ndf_from, ncells, nlayers, map_to, map_from, &
                                      alpha, beta, g_m, g_p, bw,   &
                                      ncol, nrow, ind_map_to, col_bmap_to, ind_map_from, col_bmap_from)
    ! For a columnwise matrix operator on Wtheta-W3 space
    !! data needed (in)
    ! ndf number of dofs in a cell
    ! undf number unique dofs
    ! ncells number of cells
    ! nlayers number of layers
    ! map dofmap for the space
    !! data computed (out)
    !! scalars for the CMA
    ! alpha, beta, g_m, g_p, bw, ncol, nrow
    ! memory maps
    ! indirection_map, column_banded_map
    implicit none
    integer(kind=i_def), intent(in) :: ndf_to
    integer(kind=i_def), intent(in) :: ndf_from
    integer(kind=i_def), intent(in) :: ncells
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: map_to(ndf_to,ncells)
    integer(kind=i_def), intent(in) :: map_from(ndf_from,ncells)

    integer(kind=i_def), intent(out) :: alpha
    integer(kind=i_def), intent(out) :: beta
    integer(kind=i_def), intent(out) :: g_m
    integer(kind=i_def), intent(out) :: g_p
    integer(kind=i_def), intent(out) :: bw
    integer(kind=i_def), intent(out) :: ncol
    integer(kind=i_def), intent(out) :: nrow
    integer(kind=i_def), allocatable, dimension(:),   intent(out) :: ind_map_to, ind_map_from
    integer(kind=i_def), allocatable, dimension(:,:), intent(out) :: col_bmap_to, col_bmap_from

    integer(kind=i_def) :: nface_to, ninterior_to, nface_from, ninterior_from

    ! For Wtheta-W3 the values are
    nface_to = 1_i_def
    ninterior_to  = 0_i_def

    nface_from = 0_i_def
    ninterior_from  = 1_i_def

    call get_cma_size(nlayers, ndf_to, ndf_from, nface_to, ninterior_to, nface_from, ninterior_from, &
         alpha, beta, g_m, g_p, bw, ncol, nrow)

    call allocate_maps(col_bmap_to, ind_map_to, col_bmap_from, ind_map_from, &
                           nlayers, ndf_to, ndf_from, nrow, ncol)

    call make_maps(nlayers, ndf_to, ninterior_to, nface_to, col_bmap_to, ind_map_to, map_to)
    call make_maps(nlayers, ndf_from, ninterior_from, nface_from, col_bmap_from, ind_map_from, map_from)

  end subroutine get_wtheta_w3_m3x3_cma_data


  subroutine get_wtheta_m3x3_cma_data(ndf,  undf, ncells, nlayers, map, &
                                      alpha, beta, g_m, g_p, bw,   &
                                      ncol, nrow, indirection_map, column_banded_map)
    ! For a columnwise matrix operator on Wtheta-Wtheta space
    !! data needed (in)
    ! ndf number of dofs in a cell
    ! undf number unique dofs
    ! ncells number of cells
    ! nlayers number of layers
    ! map dofmap for the space
    !! data computed (out)
    !! scalars for the CMA
    ! alpha, beta, g_m, g_p, bw, ncol, nrow
    ! memory maps
    ! indirection_map, column_banded_map
    implicit none
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf
    integer(kind=i_def), intent(in) :: ncells
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: map(ndf,ncells)

    integer(kind=i_def), intent(out) :: alpha
    integer(kind=i_def), intent(out) :: beta
    integer(kind=i_def), intent(out) :: g_m
    integer(kind=i_def), intent(out) :: g_p
    integer(kind=i_def), intent(out) :: bw
    integer(kind=i_def), intent(out) :: ncol
    integer(kind=i_def), intent(out) :: nrow
    integer(kind=i_def), allocatable, intent(out) :: indirection_map(:)
    integer(kind=i_def), allocatable, intent(out) :: column_banded_map(:,:)

    integer(kind=i_def) :: nface, ninterior
    integer(kind=i_def) :: k, df
    integer(kind=i_def) :: temp1, temp2

    ! For Wtheta-Wtheta the values are
    nface = 1_i_def
    ninterior  = 0_i_def
    call get_cma_size(nlayers, ndf, ndf, nface, ninterior, nface, ninterior, &
                       alpha, beta, g_m, g_p, bw, ncol, nrow)

    allocate(indirection_map(ncol))
    allocate(column_banded_map(ndf, nlayers))

    call make_maps(nlayers, ndf, ninterior, nface, column_banded_map, indirection_map, map)

  end subroutine get_wtheta_m3x3_cma_data


  subroutine get_cma_size(nlayers, ndf_to, ndf_from, nface_to, ninterior_to, nface_from, ninterior_from, &
       alpha, beta, g_m, g_p, bw, ncol, nrow)
    ! Compute the sizes of the CMA scalars based on the to and from function spaces
    ! inputs
    ! nlayers the number of layers
    ! ndf_to   the number of degrees of freedom per cell on the "to" space
    ! ndf_from the number of degrees of freedom per cell on the "from" space
    ! nface_to the number of degrees of freedom per cell which are face valued on the "to" space, excluding top face for double counting
    ! ninterior_to the number of degrees of freedom per cell which cell-body valued on the "to" space
    ! nface_from the number of degrees of freedom per cell which are face valued on the "from" space, excluding top face for double counting
    ! ninterior_from the number of degrees of freedom per cell which cell-body valued on the "from" space
    !! Outputs - computed CMA scalars
    ! alpha, beta, g_m, g_p, bw, ncol, nrow
    implicit none
    integer(kind=i_def), intent(in)  :: nlayers, ndf_to, ndf_from, nface_to, ninterior_to, nface_from, ninterior_from
    integer(kind=i_def), intent(out) :: alpha, beta, g_m, g_p, bw, ncol, nrow

    nrow  = nlayers * ninterior_to + (nlayers+1) * nface_to
    ncol  = nlayers * ninterior_from + (nlayers+1) * nface_from
    alpha = nface_from + ninterior_from
    beta  = nface_to + ninterior_to
    g_m   = (ndf_from -1_i_def) * (ninterior_to + nface_to)
    g_p   = (ndf_to -1_i_def) * (ninterior_from + nface_from)


    ! should compute the gcd but for small spaces this is not necessary
    bw    = 1_i_def + ceiling( (real( (g_m+g_p), kind=r_def)/ &
         (real(beta, kind=r_def))), i_def)

  end subroutine get_cma_size

  subroutine allocate_maps(col_bmap_to, ind_map_to, col_bmap_from, ind_map_from, &
       nlayers, ndf_to, ndf_from, nrow, ncol)
    ! allocate the arrays necessary for CMA
    ! output are the column banded map (to and from space) and the indirection map (to and from space)
    ! Sizes are passed in as inputs
    ! nlayers, the number of layers
    ! ndf_to, the number of degrees of freedom per cell for the "to" space
    ! ndf_from, the number of degrees of freedom per cell for the "from" space
    implicit none
    integer(kind=i_def), allocatable, dimension(:,:), intent(out) :: col_bmap_to, col_bmap_from
    integer(kind=i_def), allocatable, dimension(:),   intent(out) :: ind_map_to, ind_map_from
    integer(kind=i_def),                              intent(in)  :: nlayers, ndf_to, ndf_from, nrow, ncol

    if(.not. allocated(col_bmap_to)) &
         allocate(col_bmap_to(ndf_to, nlayers))
    if(.not. allocated(col_bmap_from)) &
         allocate(col_bmap_from(ndf_from, nlayers))

    if (.not. allocated(ind_map_to)) &
         allocate(ind_map_to(nrow))
    if(.not. allocated(ind_map_from)) &
         allocate(ind_map_from(ncol))

  end subroutine allocate_maps

  subroutine make_maps(nlayers, ndf, ninterior, nface, col_bmap, ind_map, map)
    ! compute the CMA maps from the function space dofmap
    ! input nlayers the number of layers
    ! input ninterior, the number of degrees per cell with are cell-body valued
    ! input nface, the number of degrees of freedome per cell with are face values, excluding the upper face to avoid double counting
    ! output col_bmap, two dimensional map (CMA column banded map)
    ! output ind_map, one dimensional map (CMA indirection map)
    ! input map, the dofmap for the function space
    implicit none
    integer(kind=i_def), intent(in) :: nlayers, ndf, ninterior, nface
    integer(kind=i_def), dimension(:,:), intent(inout) :: col_bmap
    integer(kind=i_def), dimension(:),   intent(inout) :: ind_map
    integer(kind=i_def), dimension(:,:), intent(in)    :: map

    integer(kind=i_def) :: k, df

    do k=1_i_def, nlayers
       do df = 1_i_def, ndf
          col_bmap(df,k) = (k-1_i_def) * (ninterior + nface) + df
          ind_map(col_bmap(df,k)) = (map(df,1_i_def) - map(1_i_def,1_i_def) + 1_i_def) + k-1_i_def
       end do
    end do

  end subroutine make_maps

  subroutine get_cma_prod_size(alpha_A, beta_A, g_p_A, g_m_A, nrow_A, &
                               alpha_B, beta_B, g_p_B, g_m_B, ncol_B, &
                               alpha_C, beta_C, g_p_C, g_m_C, bw_C, nrow_C, ncol_C)
    ! C = B * A
    ! for CMA product, the scalars are not necessarily the same as to and from function spaces so they have to be recomputed
    ! The required inputs are different for the two input matrices
    ! input 5 scalars of CMA matrix A
    ! alpha, beta, g_p, g_m and nrow
    ! input 5 scalars of CMA matrix B
    ! alpha, beta, g_p, g_m, ncol
    ! Outputs 7 scalars of the CMA product matrix
    ! alpha, beta, g_p, g_m, bw, nrow, ncol
    implicit none
    integer(kind=i_def), intent(in)  :: alpha_A, beta_A, g_p_A, g_m_A, nrow_A
    integer(kind=i_def), intent(in)  :: alpha_B, beta_B, g_p_B, g_m_B, ncol_B
    integer(kind=i_def), intent(out) :: alpha_C, beta_C, g_p_C, g_m_C, bw_C, nrow_C, ncol_C

    alpha_C = alpha_A * alpha_B
    beta_C  = beta_A  * beta_B

    g_m_C   = alpha_B * g_m_A + beta_A * g_m_B
    g_p_C   = alpha_B * g_p_A + beta_A * g_p_B
    ! should compute the gcd but for small spaces this is not necessary
    bw_C    = 1_i_def + ceiling( (g_m_C+g_p_C)/(real(beta_C, kind=r_def)), i_def)

    nrow_C = nrow_A
    ncol_C = ncol_B

  end subroutine get_cma_prod_size

  subroutine get_w3_w3_m3x3_cma_data(ndf, ncells, nlayers, map, &
                                      alpha, beta, g_m, g_p, bw,   &
                                      ncol, nrow, ind_map, col_bmap)
    ! For a columnwise matrix operator on W3-W3 space
    !! data needed (in)
    ! ndf number of dofs in a cell
    ! undf number unique dofs
    ! ncells number of cells
    ! nlayers number of layers
    ! map dofmap for the space
    !! data computed (out)
    !! scalars for the CMA
    ! alpha, beta, g_m, g_p, bw, ncol, nrow
    ! memory maps
    ! indirection_map, column_banded_map
    implicit none
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: ncells
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), dimension(ndf, ncells), intent(in) :: map

    integer(kind=i_def), intent(out) :: alpha
    integer(kind=i_def), intent(out) :: beta
    integer(kind=i_def), intent(out) :: g_m
    integer(kind=i_def), intent(out) :: g_p
    integer(kind=i_def), intent(out) :: bw
    integer(kind=i_def), intent(out) :: ncol
    integer(kind=i_def), intent(out) :: nrow
    integer(kind=i_def), allocatable, dimension(:),   intent(out) :: ind_map
    integer(kind=i_def), allocatable, dimension(:,:), intent(out) :: col_bmap

    integer(kind=i_def) :: nface, ninterior

    ! For W3-W3 the values are
    nface = 0_i_def
    ninterior  = 1_i_def

    call get_cma_size(nlayers, ndf, ndf, nface, ninterior, nface, ninterior, &
         alpha, beta, g_m, g_p, bw, ncol, nrow)

    call allocate_maps(col_bmap, ind_map, col_bmap, ind_map, &
                           nlayers, ndf, ndf, nrow, ncol)

    call make_maps(nlayers, ndf, ninterior, nface, col_bmap, ind_map, map)

  end subroutine get_w3_w3_m3x3_cma_data

  subroutine get_w2v_w2v_m3x3_cma_data(ndf, ncells, nlayers, map, &
                                      alpha, beta, g_m, g_p, bw,   &
                                      ncol, nrow, ind_map, col_bmap)
    ! For a columnwise matrix operator on W2h-W2h space
    ! any resemblance to any other space, living or dead  such as wtheta is entirely
    ! coincidental, the value of assets may go up or down, your mileage may vary, etc.
    !! data needed (in)
    ! ndf number of dofs in a cell
    ! undf number unique dofs
    ! ncells number of cells
    ! nlayers number of layers
    ! map dofmap for the space
    !! data computed (out)
    !! scalars for the CMA
    ! alpha, beta, g_m, g_p, bw, ncol, nrow
    ! memory maps
    ! indirection_map, column_banded_map
    implicit none
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: ncells
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), dimension(ndf, ncells), intent(in) :: map

    integer(kind=i_def), intent(out) :: alpha
    integer(kind=i_def), intent(out) :: beta
    integer(kind=i_def), intent(out) :: g_m
    integer(kind=i_def), intent(out) :: g_p
    integer(kind=i_def), intent(out) :: bw
    integer(kind=i_def), intent(out) :: ncol
    integer(kind=i_def), intent(out) :: nrow
    integer(kind=i_def), allocatable, dimension(:),   intent(out) :: ind_map
    integer(kind=i_def), allocatable, dimension(:,:), intent(out) :: col_bmap

    integer(kind=i_def) :: nface, ninterior

    ! For W2v-W2v the values are
    nface = 1_i_def
    ninterior  = 0_i_def

    call get_cma_size(nlayers, ndf, ndf, nface, ninterior, nface, ninterior, &
         alpha, beta, g_m, g_p, bw, ncol, nrow)

    call allocate_maps(col_bmap, ind_map, col_bmap, ind_map, &
                           nlayers, ndf, ndf, nrow, ncol)

    call make_maps(nlayers, ndf, ninterior, nface, col_bmap, ind_map, map)

  end subroutine get_w2v_w2v_m3x3_cma_data

end module get_unit_test_m3x3_cma_data_mod
