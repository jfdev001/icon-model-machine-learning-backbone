! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Contains the implementation of the div,rot,recon mathematical
! operators employed by the shallow water prototype.
!
! @par To Do
! Boundary exchange, nblks in presence of halos and dummy edge

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_math_divrot
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp, vp
USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_intp_data_strc,      ONLY: t_int_state, t_lsq
USE mo_interpol_config,     ONLY: lsq_high_set
USE mo_model_domain,        ONLY: t_patch
USE mo_grid_config,         ONLY: l_limited_area
USE mo_parallel_config,     ONLY: nproma
USE mo_exception,           ONLY: finish
USE mo_fortran_tools,       ONLY: init
use mo_lib_divrot,          ONLY: recon_lsq_cell_l_lib, recon_lsq_cell_l_svd_lib, &
                                  div_lib, div_avg_lib, &
                                  rot_vertex_atmos_lib, rot_vertex_ri_lib, &
                                  recon_lsq_cell_q_lib, recon_lsq_cell_q_svd_lib, &
                                  recon_lsq_cell_c_lib, recon_lsq_cell_c_svd_lib

! USE mo_timer,              ONLY: timer_start, timer_stop, timer_div

IMPLICIT NONE

PRIVATE

PUBLIC :: recon_lsq_cell_l, recon_lsq_cell_l_svd
PUBLIC :: recon_lsq_cell_q, recon_lsq_cell_q_svd
PUBLIC :: recon_lsq_cell_c, recon_lsq_cell_c_svd
PUBLIC :: div, div_avg
PUBLIC :: rot_vertex, rot_vertex_ri
PUBLIC :: rot_vertex_atmos

INTERFACE rot_vertex

  MODULE PROCEDURE rot_vertex_atmos

END INTERFACE


INTERFACE div

  MODULE PROCEDURE div3d
  MODULE PROCEDURE div3d_2field
  MODULE PROCEDURE div4d

END INTERFACE


CONTAINS


!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Rx = Q^T d.
!! R: upper triangular matrix (2 x 2)
!! Q: orthogonal matrix (3 x 2)
!! d: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! works only on triangular grid yet
!!
SUBROUTINE recon_lsq_cell_l( p_cc, ptr_patch, ptr_int_lsq, p_coeff, lacc, &
  &                          opt_slev, opt_elev, opt_rlstart,       &
  &                          opt_rlend, opt_lconsv, opt_acc_async )

  TYPE(t_patch), INTENT(IN)     :: &    !< patch on which computation
    &  ptr_patch                        !<is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN)          ::  &   !<  cell centered variable
    &  p_cc(:,:,:)

  LOGICAL, INTENT(IN)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< if true, conservative reconstruction is used
    &  opt_lconsv

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< optional async OpenACC
    &  opt_acc_async

  INTEGER :: slev, elev               !< vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
  LOGICAL :: l_consv

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(p_cc,2)
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF
  IF ( PRESENT(opt_lconsv) ) THEN
    l_consv = opt_lconsv
  ELSE
    l_consv = .FALSE.
  END IF

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  CALL recon_lsq_cell_l_lib( p_cc, ptr_int_lsq%lsq_idx_c, ptr_int_lsq%lsq_blk_c, &
    &                        ptr_int_lsq%lsq_qtmat_c, ptr_int_lsq%lsq_rmat_rdiag_c, &
    &                        ptr_int_lsq%lsq_rmat_utri_c, ptr_int_lsq%lsq_moments, p_coeff, &
    &                        i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                        slev, elev, nproma, l_consv, lacc=lacc, acc_async=opt_acc_async )

END SUBROUTINE recon_lsq_cell_l



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered linear
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! l    : linear reconstruction
!!
!! The least squares approach is used. Solves Ax = b via Singular
!! Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape:
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (2 x 3)
!! b: input vector (3 x 1)
!! x: solution vector (2 x 1)
!! only works on triangular grid yet
!!
SUBROUTINE recon_lsq_cell_l_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, lacc, &
  &                              opt_slev, opt_elev, opt_rlstart, opt_rlend,  &
  &                              opt_lconsv, opt_acc_async )

  TYPE(t_patch), INTENT(IN)     :: &    !< patch on which computation
    &  ptr_patch                        !< is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN)          ::  &   !<  cell centered variable
    &  p_cc(:,:,:)

  LOGICAL, INTENT(IN)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  &   !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< if true, conservative reconstruction is used
    &  opt_lconsv

  REAL(wp), INTENT(INOUT) ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)          !< (constant and gradients in latitudinal and
                                 !< longitudinal direction)

  LOGICAL, INTENT(IN), OPTIONAL ::  &   !< optional async OpenACC
    &  opt_acc_async

  INTEGER :: slev, elev              !< vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
  LOGICAL :: l_consv

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(p_cc,2)
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF
  IF ( PRESENT(opt_lconsv) ) THEN
    l_consv = opt_lconsv
  ELSE
    l_consv = .FALSE.
  END IF

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  CALL recon_lsq_cell_l_svd_lib( p_cc, ptr_int_lsq%lsq_idx_c, ptr_int_lsq%lsq_blk_c, &
    &                            ptr_int_lsq%lsq_pseudoinv, ptr_int_lsq%lsq_moments, p_coeff, &
    &                            i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                            slev, elev, nproma, l_consv, lacc=lacc, acc_async=opt_acc_async )

END SUBROUTINE recon_lsq_cell_l_svd


!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes the coefficients (derivatives) for a quadratic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (5 x 5)
!! Q  : orthogonal matrix (9 x 5)
!! Q^T: transposed of Q (5 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_q( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                          lacc, opt_slev, opt_elev, opt_rlstart, &
  &                          opt_rlend )

  TYPE(t_patch), INTENT(IN) ::   & !< patch on which computation
    &  ptr_patch                   !< is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  LOGICAL, INTENT(IN)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< derivatives, one mixed derivative and a constant
                                  !< coefficient for zonal and meridional direction

  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(p_cc,2)
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_startblk    = ptr_patch%cells%start_block(rl_start)
  i_endblk      = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  CALL recon_lsq_cell_q_lib( p_cc, ptr_int_lsq%lsq_idx_c, ptr_int_lsq%lsq_blk_c, &
    &                        ptr_int_lsq%lsq_rmat_rdiag_c, ptr_int_lsq%lsq_rmat_utri_c, &
    &                        ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_qtmat_c, p_coeff, &
    &                        i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                        slev, elev, nproma, ptr_patch%id, lsq_high_set%dim_c, l_limited_area, lacc=lacc)

END SUBROUTINE recon_lsq_cell_q



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered quadratic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! q    : quadratic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a quadratic polynomial,
!! using the least-squares method. The coefficients are provided at cell
!! centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (5 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (5 x 1)
!!
!! Coefficients:
!! p_coeff(jc,jk,jb,1) : C0
!! p_coeff(jc,jk,jb,2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb,3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb,4) : C3 (0.5*ddPhi_ddx)
!! p_coeff(jc,jk,jb,5) : C4 (0.5*ddPhi_ddy)
!! p_coeff(jc,jk,jb,6) : C5 (ddPhi_dxdy)
!!
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_q_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                              lacc, opt_slev, opt_elev, opt_rlstart, &
  &                              opt_rlend )

  TYPE(t_patch), INTENT(IN) ::   & !< patch on which computation
    &  ptr_patch                   !< is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  LOGICAL, INTENT(IN)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< derivatives, one mixed derivative and a constant
                                  !< coefficient for zonal and meridional direction

  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(p_cc,2)
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  CALL recon_lsq_cell_q_svd_lib( p_cc, ptr_int_lsq%lsq_idx_c, ptr_int_lsq%lsq_blk_c, &
    &                            ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_pseudoinv, p_coeff, &
    &                            i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                            slev, elev, nproma, ptr_patch%id, lsq_high_set%dim_c, l_limited_area, lacc=lacc)

END SUBROUTINE recon_lsq_cell_q_svd



!-------------------------------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes the coefficients (derivatives) for a cubic reconstruction,
!! using the the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!! Solves linear system Rx = Q^T d.
!! The matrices have the following size and shape:
!! R  : upper triangular matrix (9 x 9)
!! Q  : orthogonal matrix (9 x 9)
!! Q^T: transposed of Q (9 x 9)
!! d  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_c( p_cc, ptr_patch, ptr_int_lsq, p_coeff,  &
  &                          lacc, opt_slev, opt_elev, opt_rlstart,  &
  &                          opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  LOGICAL, INTENT(IN)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(p_cc,2)
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  CALL recon_lsq_cell_c_lib( p_cc, ptr_int_lsq%lsq_idx_c, ptr_int_lsq%lsq_blk_c, &
    &                        ptr_int_lsq%lsq_rmat_rdiag_c, ptr_int_lsq%lsq_rmat_utri_c, &
    &                        ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_qtmat_c, p_coeff, &
    &                        i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                        slev, elev, nproma, ptr_patch%id, lsq_high_set%dim_c, l_limited_area, lacc=lacc)

END SUBROUTINE recon_lsq_cell_c




!--------------------------------------------------
!
!
!>
!! Computes coefficients (i.e. derivatives) for cell centered cubic
!! reconstruction.
!!
!! DESCRIPTION:
!! recon: reconstruction of subgrid distribution
!! lsq  : least-squares method
!! cell : solution coefficients defined at cell center
!! c    : cubic reconstruction
!!
!! Computes unknown coefficients (derivatives) of a cubic polynomial,
!! using the least-squares method. The coefficients are provided at
!! cell centers in a local 2D cartesian system (tangential plane).
!!
!! Mathematically we solve Ax = b via Singular Value Decomposition (SVD)
!! x = PINV(A) * b
!!
!! Matrices have the following size and shape (triangular grid) :
!! PINV(A): Pseudo or Moore-Penrose inverse of A (via SVD) (9 x 9)
!! b  : input vector (LHS) (9 x 1)
!! x  : solution vector (unknowns) (9 x 1)
!!
!! Coefficients
!! p_coeff(jc,jk,jb, 1) : C0
!! p_coeff(jc,jk,jb, 2) : C1 (dPhi_dx)
!! p_coeff(jc,jk,jb, 3) : C2 (dPhi_dy)
!! p_coeff(jc,jk,jb, 4) : C3 (1/2*ddPhi_ddx)
!! p_coeff(jc,jk,jb, 5) : C4 (1/2*ddPhi_ddy)
!! p_coeff(jc,jk,jb, 6) : C5 (ddPhi_dxdy)
!! p_coeff(jc,jk,jb, 7) : C6 (1/6*dddPhi_dddx)
!! p_coeff(jc,jk,jb, 8) : C7 (1/6*dddPhi_dddy)
!! p_coeff(jc,jk,jb, 9) : C8 (1/2*dddPhi_ddxdy)
!! p_coeff(jc,jk,jb,10) : C9 (1/2*dddPhi_dxddy)
!!
!! works only on triangular grid yet
!!
!! !LITERATURE
!! Ollivier-Gooch et al (2002): A High-Order-Accurate Unstructured Mesh
!! Finite-Volume Scheme for the Advection-Diffusion Equation, J. Comput. Phys.,
!! 181, 729-752
!!
SUBROUTINE recon_lsq_cell_c_svd( p_cc, ptr_patch, ptr_int_lsq, p_coeff, &
  &                              lacc, opt_slev, opt_elev, opt_rlstart, &
  &                              opt_rlend )

  TYPE(t_patch), INTENT(IN) :: & !< patch on which computation
    &  ptr_patch                 !< is performed

  TYPE(t_lsq), TARGET, INTENT(IN) :: &  !< data structure for interpolation
    &  ptr_int_lsq

  REAL(wp), INTENT(IN) ::           & !< cell centered variable
    &  p_cc(:,:,:)

  LOGICAL, INTENT(IN)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
    &  opt_slev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level
    &  opt_elev

  INTEGER, INTENT(IN), OPTIONAL ::  & !< start and end values of refin_ctrl flag
    &  opt_rlstart, opt_rlend

  REAL(wp), INTENT(INOUT)  ::  &  !< cell based coefficients (geographical components)
    &  p_coeff(:,:,:,:)           !< physically this vector contains gradients, second
                                  !< and third derivatives, one mixed derivative and a
                                  !< constant coefficient for zonal and meridional
                                  !< direction

  INTEGER :: slev, elev           !< vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
  INTEGER :: rl_start_init, rl_end_init, i_startblk_init, i_endblk_init
  INTEGER :: i_startidx_init, i_endidx_init

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(p_cc,2)
  END IF
  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 2
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF


  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  ! values for the blocking to zero-init the lateral boundary points
  rl_start_init   = 1
  rl_end_init     = MAX(1,rl_start-1)

  i_startblk_init = ptr_patch%cells%start_block(rl_start_init)
  i_endblk_init   = ptr_patch%cells%end_block(rl_end_init)

  i_startidx_init = ptr_patch%cells%start_index(rl_start_init)
  i_endidx_init   = ptr_patch%cells%end_index(rl_end_init)

  CALL recon_lsq_cell_c_svd_lib( p_cc, ptr_int_lsq%lsq_idx_c, ptr_int_lsq%lsq_blk_c, &
    &                            ptr_int_lsq%lsq_moments, ptr_int_lsq%lsq_pseudoinv, p_coeff, &
    &                            i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                            i_startblk_init, i_endblk_init, i_startidx_init, i_endidx_init, &
    &                            slev, elev, nproma, ptr_patch%id, lsq_high_set%dim_c, l_limited_area, lacc=lacc)

END SUBROUTINE recon_lsq_cell_c_svd


!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges.
!! The midpoint rule is used for quadrature.
!! input:  lives on edges (velocity points)
!! output: lives on centers of triangles
!!
SUBROUTINE div3d( vec_e, ptr_patch, ptr_int, div_vec_c, lacc, &
  &               opt_slev, opt_elev, opt_rlstart, opt_rlend  )
!

  !  patch on which computation is performed
  TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

  ! Interpolation state
  TYPE(t_int_state), INTENT(in)     :: ptr_int

  ! edge based variable of which divergence
  ! is computed
  REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

  LOGICAL, INTENT(in)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

  INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

  INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

  ! cell based variable in which divergence is stored
  !REAL(wp), INTENT(out) ::  &
  REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

  INTEGER :: slev, elev     ! vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(vec_e,2)
  END IF

  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 1
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  !IF(ltimer) CALL timer_start(timer_div)

  CALL div_lib( vec_e, ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, &
  &             ptr_int%geofac_div, div_vec_c, &
  &             i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &             slev, elev, nproma, lacc=lacc)

  !IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div3d

SUBROUTINE div3d_2field( vec_e, ptr_patch, ptr_int, div_vec_c, lacc,           &
  &                      opt_slev, opt_elev, in2, out2, opt_rlstart, opt_rlend )

  !
  !  patch on which computation is performed
  !
  TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

  ! Interpolation state
  TYPE(t_int_state), INTENT(in)     :: ptr_int
  !
  ! edge based variable of which divergence
  ! is computed
  !
  REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

  ! second input field for more efficient processing in NH core
  REAL(wp), INTENT(in) ::  &
  &  in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

  LOGICAL, INTENT(in)           ::  &   !< if true, use OpenACC
    &  lacc

  INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

  INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

  INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

  !
  ! cell based variable in which divergence is stored
  !
  !REAL(wp), INTENT(out) ::  &
  REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

  ! second output field
  REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

  INTEGER :: slev, elev     ! vertical start and end level
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

  !-----------------------------------------------------------------------

  ! check optional arguments
  IF ( PRESENT(opt_slev) ) THEN
    slev = opt_slev
  ELSE
    slev = 1
  END IF
  IF ( PRESENT(opt_elev) ) THEN
    elev = opt_elev
  ELSE
    elev = UBOUND(vec_e,2)
  END IF

  IF ( PRESENT(opt_rlstart) ) THEN
    rl_start = opt_rlstart
  ELSE
    rl_start = 1
  END IF
  IF ( PRESENT(opt_rlend) ) THEN
    rl_end = opt_rlend
  ELSE
    rl_end = min_rlcell
  END IF

  ! values for the blocking
  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  !IF(ltimer) CALL timer_start(timer_div)

  CALL div_lib( vec_e, ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, &
  &           ptr_int%geofac_div, div_vec_c, in2, out2, &
  &           i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &           slev, elev, nproma, lacc=lacc)

  !IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div3d_2field

!-------------------------------------------------------------------------
!
!
!>
!! Special version of div that processes 4D fields in one step
!!
!! See standard routine (div3d) for further description
!!
SUBROUTINE div4d( ptr_patch, ptr_int, f4din, f4dout, dim4d, lacc, &
  &               opt_slev, opt_elev, opt_rlstart, opt_rlend      )

!  patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
! edge based 4D input field of which divergence is computed
!
REAL(wp), INTENT(in) ::  &
  &  f4din(:,:,:,:) ! dim: (nproma,nlev,nblks_e,dim4d)

INTEGER, INTENT(in) :: dim4d ! Last dimension of the input/output fields

LOGICAL, INTENT(in)           ::  &   !< if true, use OpenACC
  &  lacc

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev(dim4d)    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev(dim4d)    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

! 4D cell based variable in which divergence is stored
REAL(vp), INTENT(inout) ::  &
  &  f4dout(:,:,:,:) ! dim: (nproma,nlev,nblks_c,dim4d)

INTEGER :: slev(dim4d), elev(dim4d)     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(f4din,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

!IF(ltimer) CALL timer_start(timer_div)

CALL div_lib( ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, &
  &           ptr_int%geofac_div, f4din, f4dout, dim4d, &
  &           i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &           slev, elev, nproma, lacc=lacc)

!IF(ltimer) CALL timer_stop(timer_div)

END SUBROUTINE div4d

!-------------------------------------------------------------------------
!
!
!>
!! Computes discrete divergence of a vector field.
!!
!! Computes discrete divergence of a vector field
!! given by its components in the directions normal to triangle edges,
!! followed by bilinear averaging to remove checkerboard noise
!! (Combines div_midpoint and cell_avg_varwgt to increase computing efficiency)
!!
SUBROUTINE div_avg( vec_e, ptr_patch, ptr_int, avg_coeff, div_vec_c, lacc,   &
  &                 opt_in2, opt_out2, opt_slev, opt_elev, opt_rlstart,      &
  &                 opt_rlend )

!  patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  averaging coefficients
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)
!
! edge based variable of which divergence
! is computed
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

LOGICAL, INTENT(in)           ::  &   !< if true, use OpenACC
  &  lacc

! optional second input field for more efficient processing in NH core
REAL(wp), OPTIONAL, INTENT(in) ::  &
  &  opt_in2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
! cell based variable in which divergence is stored
!
REAL(wp), INTENT(inout) ::  &
  &  div_vec_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

! optional second output field
REAL(wp), OPTIONAL, INTENT(inout) ::  &
  &  opt_out2(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end, rl_start_l2, rl_end_l1
INTEGER :: i_startblk_in(3), i_endblk_in(3), i_startidx_in(3), i_endidx_in(3)

LOGICAL :: l2fields


!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(vec_e,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 1
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell
END IF
IF ( PRESENT(opt_in2) .AND. PRESENT(opt_out2)) THEN
  l2fields = .TRUE.
ELSE
  l2fields = .FALSE.
ENDIF

rl_start_l2 = rl_start + 1

IF ( PRESENT(opt_rlend) .AND. rl_end < 0 .AND. rl_end > min_rlcell ) THEN
  rl_end_l1 = rl_end - 1
ELSE
  rl_end_l1 = rl_end
END IF

! values for the blocking
i_startblk_in(1) = ptr_patch%cells%start_block(rl_start)
i_endblk_in(1)   = ptr_patch%cells%end_block(rl_end_l1)

i_startblk_in(2) = ptr_patch%cells%start_block(rl_start)
i_endblk_in(2)   = ptr_patch%cells%end_block(rl_start_l2)

i_startblk_in(3) = ptr_patch%cells%start_block(rl_start_l2)
i_endblk_in(3)   = ptr_patch%cells%end_block(rl_end)

i_startidx_in(1) = ptr_patch%cells%start_index(rl_start)
i_endidx_in(1)   = ptr_patch%cells%end_index(rl_end_l1)

i_startidx_in(2) = ptr_patch%cells%start_index(rl_start)
i_endidx_in(2)   = ptr_patch%cells%end_index(rl_start_l2)

i_startidx_in(3) = ptr_patch%cells%start_index(rl_start_l2)
i_endidx_in(3)   = ptr_patch%cells%end_index(rl_end)

CALL div_avg_lib( vec_e, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, &
  &               ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, &
  &               ptr_int%geofac_div, avg_coeff, div_vec_c, opt_in2, opt_out2, &
  &               i_startblk_in, i_endblk_in, i_startidx_in, i_endidx_in, &
  &               ptr_patch%nlev, ptr_patch%nblks_c, ptr_patch%id, l_limited_area, slev, elev, nproma, l2fields, lacc=lacc)

END SUBROUTINE div_avg

!-------------------------------------------------------------------------
!
!>
!! Computes discrete rotation.
!!
!! Computes discrete rotation at
!! vertices of triangle cells (centers of dual hexagon cells)
!! from a vector field given by its components in the directions normal
!! to triangle edges.
!! input:  lives on edges (velocity points)
!! output: lives on dual of cells (vertices for triangular grid)
!!
SUBROUTINE rot_vertex_atmos( vec_e, ptr_patch, ptr_int, rot_vec, lacc,  &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )

!  patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  edge based variable of which rotation is computed
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

LOGICAL, INTENT(in)           ::  &   !< if true, use OpenACC
  &  lacc

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!  vertex based variable in which rotation is stored
REAL(wp), INTENT(inout) ::  &
  &  rot_vec(:,:,:) ! dim: (nproma,nlev,nblks_v or nblks_e)

INTEGER :: slev, elev     ! vertical start and end level

INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(vec_e,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:rot_vertex_atmos',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

  ! values for the blocking
  i_startblk = ptr_patch%verts%start_block(rl_start)
  i_endblk   = ptr_patch%verts%end_block(rl_end)

  i_startidx_in = ptr_patch%verts%start_index(rl_start)
  i_endidx_in   = ptr_patch%verts%end_index(rl_end)

  CALL rot_vertex_atmos_lib( vec_e, ptr_patch%verts%edge_idx, ptr_patch%verts%edge_blk, &
                       ptr_int%geofac_rot, rot_vec, &
    &                  i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                  slev, elev, nproma, lacc=lacc)

END SUBROUTINE rot_vertex_atmos

!>
!! Same as above routine, but expects reversed index order (vertical first)
!! of the output field if __LOOP_EXCHANGE is specified. In addition, the
!! output field (vorticity) has single precision if __MIXED_PRECISION is specified
!!
!!
SUBROUTINE rot_vertex_ri( vec_e, ptr_patch, ptr_int, rot_vec, lacc,    &
  &                       opt_slev, opt_elev, opt_rlend, opt_acc_async )

!  patch on which computation is performed
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  edge based variable of which rotation is computed
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

LOGICAL, INTENT(in)           ::  &   !< if true, use OpenACC
  &  lacc

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlend   ! end value of refin_ctrl flag

LOGICAL, INTENT(IN), OPTIONAL ::  &
  &  opt_acc_async ! optional async OpenACC

!  vertex based variable in which rotation is stored
REAL(vp), INTENT(inout) ::  &
  &  rot_vec(:,:,:) ! dim: (nproma,nlev,nblks_v) or (nlev,nproma,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level

INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(vec_e,2)
END IF

rl_start = 2

IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

  ! values for the blocking
  i_startblk = ptr_patch%verts%start_block(rl_start)
  i_endblk   = ptr_patch%verts%end_block(rl_end)

  i_startidx_in = ptr_patch%verts%start_index(rl_start)
  i_endidx_in   = ptr_patch%verts%end_index(rl_end)

  CALL rot_vertex_ri_lib( vec_e, ptr_patch%verts%edge_idx, ptr_patch%verts%edge_blk, &
                          ptr_int%geofac_rot, rot_vec, &
    &                     i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
    &                     slev, elev, nproma, lacc=lacc, acc_async=opt_acc_async )

END SUBROUTINE rot_vertex_ri

END MODULE mo_math_divrot
