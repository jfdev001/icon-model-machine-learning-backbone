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

! Contains the implementation of the nabla mathematical operators.
!
! Contains the implementation of the mathematical operators
! employed by the shallow water prototype.
!
! @par To Do
! Boundary exchange, nblks in presence of halos and dummy edge

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_math_laplace
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                     ONLY: wp
USE mo_impl_constants,           ONLY: min_rlcell, min_rledge, min_rlvert
USE mo_intp_data_strc,           ONLY: t_int_state
USE mo_model_domain,             ONLY: t_patch
USE mo_grid_config,              ONLY: l_limited_area
USE mo_parallel_config,          ONLY: nproma, p_test_run
USE mo_exception,                ONLY: finish
USE mo_sync,                     ONLY: SYNC_C, SYNC_E, sync_patch_array
USE mo_math_gradients,           ONLY: grad_fd_norm
USE mo_fortran_tools,            ONLY: copy
USE mo_lib_divrot,               ONLY: div_lib, rot_vertex_atmos_lib
use mo_lib_laplace,              ONLY: nabla2_vec_atmos_lib, nabla2_scalar_lib, nabla2_scalar_avg_lib
USE mo_lib_interpolation_scalar, ONLY: edges2verts_scalar_lib, verts2edges_scalar_lib

IMPLICIT NONE

PRIVATE

PUBLIC :: nabla2_vec
PUBLIC :: nabla2_scalar, nabla2_scalar_avg
PUBLIC :: nabla4_vec
PUBLIC :: nabla4_scalar

INTERFACE nabla2_vec

  MODULE PROCEDURE nabla2_vec_atmos

END INTERFACE


CONTAINS


!-------------------------------------------------------------------------
!
!>
!!  Computes  laplacian of a vector field.
!!
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
SUBROUTINE nabla2_vec_atmos( vec_e, ptr_patch, ptr_int, nabla2_vec_e, lacc, &
  &                          opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which laplacian is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

LOGICAL, INTENT(in) ::  &
&  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: je, jk, jb
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_c, rl_end_c, rl_start_v, rl_end_v
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v

REAL(wp) ::  &
  &  z_div_c(nproma,ptr_patch%nlev,ptr_patch%nblks_c),  &
  &  z_rot_v(nproma,ptr_patch%nlev,ptr_patch%nblks_v)

INTEGER,  DIMENSION(:,:,:),   POINTER :: icidx, icblk, ividx, ivblk

!-----------------------------------------------------------------------
IF (p_test_run) THEN
  z_div_c(:,:,:)=0.0_wp
  z_rot_v(:,:,:)=0.0_wp
ENDIF

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
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 2)) THEN
    CALL finish ('mo_math_operators:nabla2_vec_atmos',  &
          &      'opt_rlstart must not be between 0 and 2')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 3
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

rl_start_c = rl_start/2

IF (rl_start > 0) THEN
  rl_start_v = (rl_start+1)/2
ELSE
  rl_start_v = (rl_start-1)/2
ENDIF

IF (rl_end > 0) THEN
  rl_end_c = (rl_end+1)/2
  rl_end_v = rl_end/2+1
ELSE
  rl_end_c = (rl_end-1)/2
  rl_end_v = rl_end/2-1
ENDIF

rl_end_c = MAX(min_rlcell,rl_end_c)
rl_end_v = MAX(min_rlvert,rl_end_v)

icidx => ptr_patch%edges%cell_idx
icblk => ptr_patch%edges%cell_blk
ividx => ptr_patch%edges%vertex_idx
ivblk => ptr_patch%edges%vertex_blk

i_startblk = ptr_patch%edges%start_block(rl_start)
i_endblk   = ptr_patch%edges%end_block(rl_end)

i_startidx_in = ptr_patch%edges%start_index(rl_start)
i_endidx_in   = ptr_patch%edges%end_index(rl_end)

! This values will be needed to call div_lib
i_startblk_c = ptr_patch%cells%start_block(rl_start_c)
i_endblk_c   = ptr_patch%cells%end_block(rl_end_c)

i_startidx_c = ptr_patch%cells%start_index(rl_start_c)
i_endidx_c   = ptr_patch%cells%end_index(rl_end_c)

! This values will be needed to call rot_vertex_lib
i_startblk_v = ptr_patch%verts%start_block(rl_start_v)
i_endblk_v   = ptr_patch%verts%end_block(rl_end_v)

i_startidx_v = ptr_patch%verts%start_index(rl_start_v)
i_endidx_v   = ptr_patch%verts%end_index(rl_end_v)

SELECT CASE (ptr_patch%geometry_info%cell_type)

CASE (3) ! (cell_type == 3)

CALL nabla2_vec_atmos_lib( vec_e, &
  &                        ptr_patch%edges%cell_idx, ptr_patch%edges%cell_blk, & ! require to calculate nabla2
  &                        ptr_patch%edges%vertex_idx, ptr_patch%edges%vertex_blk, & ! required to calculate nabla2
  &                        ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, & ! required for div_lib
  &                        ptr_patch%verts%edge_idx, ptr_patch%verts%edge_blk, & ! required for rot_vertex_lib
  &                        ptr_patch%edges%tangent_orientation, ptr_patch%edges%inv_primal_edge_length, &
  &                        ptr_patch%edges%inv_dual_edge_length, ptr_int%geofac_div, ptr_int%geofac_rot, &
  &                        nabla2_vec_e, & ! main output vector
  &                        i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, & ! required for div_lib
  &                        i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v, & ! required for rot_vertex_lib
  &                        i_startblk, i_endblk, i_startidx_in, i_endidx_in, & ! this four are needed to call get_indices_e_lib
  &                        ptr_patch%nlev, ptr_patch%nblks_c, ptr_patch%nblks_v, slev, elev, nproma, lacc=lacc )

END SELECT


END SUBROUTINE nabla2_vec_atmos


!-------------------------------------------------------------------------
!

!>
!! Computes biharmonic laplacian @f$\nabla ^4@f$ of a vector field without boundaries as used in atmospheric model.
!!
!! Computes biharmonic laplacian @f$\nabla ^4@f$ of a vector field without boundaries as used in atmospheric model.
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
SUBROUTINE nabla4_vec( vec_e, ptr_patch, ptr_int, nabla4_vec_e, lacc, &
  &                    opt_nabla2, opt_slev, opt_elev, opt_rlstart, opt_rlend )

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int
!
!  edge based variable of which laplacian is computed
!
REAL(wp), INTENT(in) ::  &
  &  vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

LOGICAL, INTENT(in) ::  &
  &  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla4_vec_e(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Optional argument for passing nabla2 to the calling program
! (to avoid double computation for Smagorinsky diffusion and nest boundary diffusion)
REAL(wp), INTENT(inout), TARGET, OPTIONAL  ::  &
  &  opt_nabla2(:,:,:) ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: rl_start_s1, rl_end_s1

REAL(wp), ALLOCATABLE, TARGET :: z_nabla2_vec_e(:,:,:) ! dim: (nproma,nlev,ptr_patch%nblks_e)
REAL(wp), POINTER :: p_nabla2(:,:,:)

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
  IF ((opt_rlstart >= 0) .AND. (opt_rlstart <= 4)) THEN
    CALL finish ('mo_math_operators:nabla4_vec',  &
          &      'opt_rlstart must not be between 0 and 4')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 5
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

IF (rl_start > 0) THEN
  rl_start_s1 = rl_start - 2
ELSE
  rl_start_s1 = rl_start + 2
ENDIF

IF (rl_end > 0) THEN
  rl_end_s1 = rl_end + 2
ELSE
  rl_end_s1 = rl_end - 2
ENDIF

rl_end_s1 = MAX(min_rledge,rl_end_s1)

IF (PRESENT(opt_nabla2) ) THEN
  p_nabla2 => opt_nabla2
ELSE
  ALLOCATE (z_nabla2_vec_e(nproma,ptr_patch%nlev,ptr_patch%nblks_e))

  p_nabla2 => z_nabla2_vec_e
ENDIF

!$ACC DATA CREATE(z_nabla2_vec_e) PRESENT(vec_e) PRESENT(nabla4_vec_e) &
!$ACC   PRESENT(ptr_patch, ptr_int) IF(lacc)

!
! apply second order Laplacian twice
!
IF (p_test_run) THEN
  p_nabla2(:,:,:) = 0.0_wp
!   rl_start_s1 = 1
!   rl_end_s1 = min_rledge
ENDIF

CALL nabla2_vec( vec_e, ptr_patch, ptr_int, p_nabla2, lacc=lacc, &
  &              opt_slev=slev, opt_elev=elev, opt_rlstart=rl_start_s1, &
  &              opt_rlend=rl_end_s1 )

CALL sync_patch_array(SYNC_E, ptr_patch, p_nabla2, lacc=lacc)

CALL nabla2_vec( p_nabla2, ptr_patch, ptr_int, nabla4_vec_e, lacc=lacc, &
  &              opt_slev=slev, opt_elev=elev, opt_rlstart=rl_start, &
  &              opt_rlend=rl_end )

IF (.NOT. PRESENT(opt_nabla2) ) THEN
  DEALLOCATE (z_nabla2_vec_e)
ENDIF

!$ACC END DATA

END SUBROUTINE nabla4_vec
!-----------------------------------------------------------------------

!>
!!  Computes laplacian @f$\nabla ^2 @f$ of a scalar field.
!!
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
SUBROUTINE nabla2_scalar( psi_c, ptr_patch, ptr_int, nabla2_psi_c, &
  &                       lacc, slev, elev, rl_start, rl_end )

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch         !< patch on which computation is performed
TYPE(t_int_state),     INTENT(in) :: ptr_int           !< interpolation state

REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) !< cells based variable of which biharmonic laplacian is computed, dim: (nproma,nlev,nblks_c)

LOGICAL,               INTENT(in) :: lacc              !< if TRUE, use OpenACC
INTEGER,               INTENT(in) :: slev              !< vertical start level
INTEGER,               INTENT(in) :: elev              !< vertical end level
INTEGER,               INTENT(in) :: rl_start,rl_end   !< start and end values of refin_ctrl flag

! cell based variable in which biharmonic laplacian is stored
REAL(wp), INTENT(inout) ::  &
  &  nabla2_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e

!-----------------------------------------------------------------------

! values for the blocking
i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

i_startblk_e = ptr_patch%edges%start_block(rl_start)
i_endblk_e   = ptr_patch%edges%end_block(rl_end)

i_startidx_e = ptr_patch%edges%start_index(rl_start)
i_endidx_e   = ptr_patch%edges%end_index(rl_end)

CALL nabla2_scalar_lib( psi_c, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, &
  &                     ptr_patch%edges%cell_idx, ptr_patch%edges%cell_blk, ptr_patch%edges%inv_dual_edge_length, & ! required for grad_fd_norm_lib
  &                     ptr_patch%cells%edge_idx, ptr_patch%cells%edge_blk, & ! required for div_lib
  &                     ptr_int%geofac_n2s, ptr_int%geofac_div, nabla2_psi_c, &
  &                     i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                     i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
  &                     ptr_patch%nlev, slev, elev, nproma, ptr_patch%nblks_e, ptr_patch%geometry_info%cell_type, lacc=lacc )

END SUBROUTINE nabla2_scalar

!-------------------------------------------------------------------------
!

!>
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging.
!!
!!  Computes Laplacian @f$\nabla ^2 @f$ of a scalar field, followed by weighted averaging
!!  with the neighboring cells to increase computing efficiency.
!!  NOTE: This optimized routine works for triangular grids only.
!! input:  lives on cells (mass points)
!! output: lives on cells
!!
SUBROUTINE nabla2_scalar_avg( psi_c, ptr_patch, ptr_int, avg_coeff, nabla2_psi_c, &
  &                           lacc, opt_slev, opt_elev )
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation state
TYPE(t_int_state), INTENT(in)     :: ptr_int

!  averaging coefficients
REAL(wp), INTENT(in) :: avg_coeff(:,:,:) ! dim: (nproma,nlev,nblks_c)

!
!  cells based variable of which biharmonic laplacian is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

LOGICAL, INTENT(in) ::  &
  &  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level


!
!  cell based variable in which biharmonic laplacian is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  nabla2_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end, rl_start_l2
INTEGER :: i_startblk_in(3), i_endblk_in(3), i_startidx_in(3), i_endidx_in(3)

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
  elev = UBOUND(psi_c,2)
END IF

rl_start = 2
rl_start_l2 = rl_start + 1
rl_end = min_rlcell

i_startblk_in(1) = ptr_patch%cells%start_block(rl_start)
i_endblk_in(1)   = ptr_patch%cells%end_block(rl_end)

i_startidx_in(1) = ptr_patch%cells%start_index(rl_start)
i_endidx_in(1)   = ptr_patch%cells%end_index(rl_end)

i_startblk_in(2) = ptr_patch%cells%start_block(rl_start)
i_endblk_in(2)   = ptr_patch%cells%end_block(rl_start_l2)

i_startidx_in(2) = ptr_patch%cells%start_index(rl_start)
i_endidx_in(2)   = ptr_patch%cells%end_index(rl_start_l2)

i_startblk_in(3) = ptr_patch%cells%start_block(rl_start_l2)
i_endblk_in(3)   = ptr_patch%cells%end_block(rl_end)

i_startidx_in(3) = ptr_patch%cells%start_index(rl_start_l2)
i_endidx_in(3)   = ptr_patch%cells%end_index(rl_end)

CALL nabla2_scalar_avg_lib( psi_c, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, &
  &                         ptr_int%geofac_n2s, avg_coeff, nabla2_psi_c, &
  &                         i_startblk_in, i_endblk_in, i_startidx_in, i_endidx_in, &
  &                         ptr_patch%nblks_c, ptr_patch%geometry_info%cell_type, ptr_patch%id, &
  &                         ptr_patch%nlev, slev, elev, nproma, l_limited_area, lacc=lacc)

END SUBROUTINE nabla2_scalar_avg

!-------------------------------------------------------------------------
!
!

!>
!!  Computes biharmonic laplacian @f$\nabla ^4 @f$ of a scalar field.
!!
!! input:  lives on edges (velocity points)
!! output: lives on edges
!!
SUBROUTINE nabla4_scalar( psi_c, ptr_patch, ptr_int, nabla4_psi_c, &
  &                       lacc, slev, elev, rl_start, rl_end, p_nabla2  )

TYPE(t_patch), TARGET, INTENT(in)    :: ptr_patch           !< patch on which computation is performed
TYPE(t_int_state),     INTENT(in)    :: ptr_int             !< interpolation state

REAL(wp),              INTENT(in)    ::  &
  &  psi_c(:,:,:) !< cells based variable of which biharmonic laplacian is computed, dim: (nproma,nlev,nblks_c)

LOGICAL,               INTENT(in)    ::  lacc               !< if TRUE, use OpenACC
INTEGER,               INTENT(in)    ::  slev               !< vertical start level
INTEGER,               INTENT(in)    ::  elev               !< vertical end level
INTEGER,               INTENT(in)    ::  rl_start, rl_end   !< start and end values of refin_ctrl flag

! cell based variable in which biharmonic laplacian is stored
REAL(wp), INTENT(inout) ::  nabla4_psi_c(:,:,:) ! dim: (nproma,nlev,nblks_c)
! argument for passing nabla2 to the calling program
! (to avoid double computation for Smagorinsky diffusion and nest boundary diffusion)
REAL(wp), INTENT(inout)  ::  &
  &  p_nabla2(:,:,:) ! dim: (nproma,nlev,nblks_e)
INTEGER :: rl_start_s1, rl_end_s1

!-----------------------------------------------------------------------

IF (rl_start > 0) THEN
  rl_start_s1 = rl_start - 1
ELSE
  rl_start_s1 = rl_start + 1
ENDIF
IF (rl_end > 0) THEN
  rl_end_s1 = rl_end + 1
ELSE
  rl_end_s1 = rl_end - 1
ENDIF

rl_end_s1 = MAX(min_rlcell,rl_end_s1)

!$ACC DATA PRESENT(psi_c) PRESENT(nabla4_psi_c) &
!$ACC   PRESENT(ptr_patch, ptr_int) IF(lacc)

! apply second order Laplacian twice
IF (p_test_run) p_nabla2(:,:,:) = 0.0_wp

CALL nabla2_scalar( psi_c, ptr_patch, ptr_int, p_nabla2, lacc, &
                    slev, elev, rl_start=rl_start_s1, rl_end=rl_end_s1 )

CALL sync_patch_array(SYNC_C, ptr_patch, p_nabla2, lacc=lacc)

CALL nabla2_scalar( p_nabla2, ptr_patch, ptr_int, nabla4_psi_c, lacc, &
                    slev, elev, rl_start=rl_start, rl_end=rl_end )

!$ACC END DATA

END SUBROUTINE nabla4_scalar


END MODULE mo_math_laplace
