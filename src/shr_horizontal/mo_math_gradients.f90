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

! Contains the implementation of the mathematical grad operators
! employed by the shallow water prototype.
!
! @par To Do
! Boundary exchange, nblks in presence of halos and dummy edge

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_math_gradients
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                      ONLY: wp, vp
USE mo_impl_constants,            ONLY: min_rlcell, min_rledge
USE mo_intp_data_strc,            ONLY: t_int_state
USE mo_model_domain,              ONLY: t_patch
USE mo_parallel_config,           ONLY: nproma
USE mo_run_config,                ONLY: timers_level
USE mo_exception,                 ONLY: finish
USE mo_timer,                     ONLY: timer_start, timer_stop, timer_grad
USE mo_fortran_tools,             ONLY: init
USE mo_lib_interpolation_scalar,  ONLY: cells2edges_scalar_lib
USE mo_lib_gradients,             ONLY: grad_fd_norm_lib, grad_fd_tang_lib
USE mo_lib_gradients,             ONLY: grad_fe_cell_lib, grad_green_gauss_cell_lib
USE mo_grid_config,               ONLY: l_limited_area

IMPLICIT NONE

PRIVATE

PUBLIC :: grad_fd_norm, grad_fd_tang
PUBLIC :: grad_green_gauss_cell
PUBLIC :: grad_fe_cell

INTERFACE grad_green_gauss_cell
  MODULE PROCEDURE grad_green_gauss_cell_adv
  MODULE PROCEDURE grad_green_gauss_cell_dycore
END INTERFACE

INTERFACE grad_fe_cell
  MODULE PROCEDURE grad_fe_cell_3d
  MODULE PROCEDURE grad_fe_cell_2d
END INTERFACE

CONTAINS

!-------------------------------------------------------------------------
!

!-------------------------------------------------------------------------
!
!!  Computes directional  derivative of a cell centered variable.
!!
!!  Computes directional  derivative of a cell centered variable
!!  with respect to direction normal to triangle edge.
!! input: lives on centres of triangles
!! output:  lives on edges (velocity points)
!!
SUBROUTINE grad_fd_norm( psi_c, ptr_patch, grad_norm_psi_e, lacc,   &
  &                      opt_slev, opt_elev, opt_rlstart, opt_rlend )

!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
!  cell based variable of which normal derivative is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_c(:,:,:)       ! dim: (nproma,nlev,nblks_c)

LOGICAL, INTENT(in) ::  &
  &  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which normal derivative is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  grad_norm_psi_e(:,:,:)  ! dim: (nproma,nlev,nblks_e)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!
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

IF ( PRESENT(opt_rlstart) ) THEN
  ! rl_start=1 means edges located along a lateral boundary of a nested
  ! domain. For those, gradient computation is not possible
  IF (opt_rlstart == 1) THEN
    CALL finish ('mo_math_operators:grad_fd_norm',  &
          &      'opt_rlstart must not be equal to 1')
  ENDIF
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rledge
END IF

i_startblk = ptr_patch%edges%start_block(rl_start)
i_endblk   = ptr_patch%edges%end_block(rl_end)

i_startidx_in = ptr_patch%edges%start_index(rl_start)
i_endidx_in   = ptr_patch%edges%end_index(rl_end)

IF (timers_level > 10) CALL timer_start(timer_grad)

CALL grad_fd_norm_lib( psi_c, ptr_patch%edges%cell_idx, ptr_patch%edges%cell_blk, &
  &                    ptr_patch%edges%inv_dual_edge_length, grad_norm_psi_e, &
  &                    i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                    slev, elev, nproma, lacc=lacc )

IF (timers_level > 10) CALL timer_stop(timer_grad)

END SUBROUTINE grad_fd_norm

!-------------------------------------------------------------------------
!
! RESTRUCT: @Marco: please adjust calls to this routine to your needs.
!!
!! Computes directional derivative of a vertex centered variable with.
!!
!! Computes directional derivative of a vertex centered variable with
!! respect to direction tanget to triangle edge. Notice that the
!! tangential direction is defined by
!!   iorient*(vertex2 - vertex1)
!! input: lives on vertices of triangles
!! output: lives on edges (velocity points)
!!
SUBROUTINE grad_fd_tang( psi_v, ptr_patch, grad_tang_psi_e, lacc,   &
  &                      opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
!
! vertex based variable of which tangential derivative is computed
!
REAL(wp), INTENT(in) ::  &
  &  psi_v(:,:,:)

LOGICAL, INTENT(in) ::  &
  &  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

!
!  edge based variable in which tangential derivative is stored
!
!REAL(wp), INTENT(out) ::  &
REAL(wp), INTENT(inout) ::  &
  &  grad_tang_psi_e(:,:,:)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
!
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
  elev = UBOUND(psi_v,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  ! The possible domain extent depends on the reconstruction algorithm
  ! used to calculate psi_v; the following values are valid for a 6-point
  ! stencil. In the hexagon case (where prognostic variables are located
  ! at vertices), rl_start may be set to 1.
  IF ((opt_rlstart >= 1) .AND. (opt_rlstart <= 2)) THEN
    CALL finish ('mo_math_operators:grad_fd_tang',  &
          &      'opt_rlstart must not be equal to 1 or 2')
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

! values for the blocking
i_startblk = ptr_patch%edges%start_block(rl_start)
i_endblk   = ptr_patch%edges%end_block(rl_end)

i_startidx_in = ptr_patch%edges%start_index(rl_start)
i_endidx_in   = ptr_patch%edges%end_index(rl_end)

CALL grad_fd_tang_lib( psi_v, ptr_patch%edges%vertex_idx, ptr_patch%edges%vertex_blk, &
  &                    ptr_patch%edges%primal_edge_length,ptr_patch%edges%tangent_orientation, grad_tang_psi_e, &
  &                    i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                    slev, elev, nproma, lacc=lacc )

END SUBROUTINE grad_fd_tang

!-------------------------------------------------------------------------
!
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! The triangular element is spanned by the cell circumcenters of the three
!! direct neighbours. In contrast to the Green-Gauss approach, this
!! approach does not involve the cell center value of the central triangle.
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons, Sec. 7.2, 7.6
!!
!!
SUBROUTINE grad_fe_cell_3d( p_cc, ptr_patch, ptr_int, p_grad, &
  &                         lacc, opt_slev, opt_elev,         &
  &                         opt_rlstart, opt_rlend            )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:,:)

LOGICAL, INTENT(in) ::  &
  &  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(vp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)      ! dim:(2,nproma,nlev,nblks_c)

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

i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

CALL grad_fe_cell_lib( p_cc, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, &
  &                        ptr_int%gradc_bmat, p_grad, &
  &                        i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                        slev, elev, nproma, ptr_patch%id, lacc=lacc )

END SUBROUTINE grad_fe_cell_3d




!-------------------------------------------------------------------------
!
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The gradient is computed by taking the derivative of the shape functions
!! for a three-node triangular element (Finite Element thinking).
!! The triangular element is spanned by the cell circumcenters of the three
!! direct neighbours. In contrast to the Green-Gauss approach, this
!! approach does not involve the cell center value of the central triangle.
!!
!! 2D version, i.e. for a single vertical level
!!
!! LITERATURE:
!! Fish. J and T. Belytschko, 2007: A first course in finite elements,
!!                                  John Wiley and Sons, Sec. 7.2, 7.6
!!
!!
SUBROUTINE grad_fe_cell_2d( p_cc, ptr_patch, ptr_int, p_grad, &
  &                         lacc, opt_rlstart, opt_rlend      )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:)

LOGICAL, INTENT(in) ::  &
  &  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(wp), INTENT(inout) ::  &
  &  p_grad(:,:,:)      ! dim:(2,nproma,nblks_c)

INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in

!-----------------------------------------------------------------------

! check optional arguments
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

i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

CALL grad_fe_cell_lib( p_cc, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, &
  &                    ptr_int%gradc_bmat, p_grad, &
  &                    i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                    nproma, ptr_patch%id, lacc=lacc )

END SUBROUTINE grad_fe_cell_2d

!-------------------------------------------------------------------------
!
!! Computes the cell centered gradient in geographical coordinates.
!!
!! The Green-Gauss approach is used. See for example:
!! http://www.cfd-online.com/Wiki/Gradient_computation
!!
SUBROUTINE grad_green_gauss_cell_adv( p_cc, ptr_patch, ptr_int, p_grad, lacc, &
  &                                   opt_slev, opt_elev, opt_p_face,   &
  &                                   opt_rlstart, opt_rlend            )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
!
!  data structure for interpolation
!
TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

!
!  cell centered variable
!
REAL(wp), INTENT(in) ::  &
  &  p_cc(:,:,:)

LOGICAL, INTENT(in) ::  &
  &  lacc    ! if TRUE, use OpenACC

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_slev    ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_elev    ! optional vertical end level

INTEGER, INTENT(in), OPTIONAL ::  &
  &  opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag
!
! cell based Green-Gauss reconstructed geographical gradient vector
!
REAL(vp), INTENT(inout) ::  &
  &  p_grad(:,:,:,:)      ! dim:(2,nproma,nlev,nblks_c)

! optional: calculated face values of cell centered quantity
REAL(wp), INTENT(inout), OPTIONAL ::  &
  &  opt_p_face(:,:,:)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx_in, i_endidx_in
INTEGER :: i_startblk_opt(2), i_endblk_opt(2), i_startidx_opt(2), i_endidx_opt(2)

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

i_startblk = ptr_patch%cells%start_block(rl_start)
i_endblk   = ptr_patch%cells%end_block(rl_end)

i_startidx_in = ptr_patch%cells%start_index(rl_start)
i_endidx_in   = ptr_patch%cells%end_index(rl_end)

! save face values in optional output field
! (the cell-to-edge interpolation is no longer needed otherwise because
!  of using precomputed geometrical factors)
IF ( PRESENT(opt_p_face) ) THEN

  i_startblk_opt(1) = ptr_patch%edges%start_block(1)
  i_endblk_opt(1)   = ptr_patch%edges%end_block(1)

  i_startblk_opt(2) = ptr_patch%edges%start_block(rl_start)
  i_endblk_opt(2)   = ptr_patch%edges%end_block(rl_end)

  i_startidx_opt(1) = ptr_patch%edges%start_index(1)
  i_endidx_opt(1)   = ptr_patch%edges%end_index(1)

  i_startidx_opt(2) = ptr_patch%edges%start_index(rl_start)
  i_endidx_opt(2)   = ptr_patch%edges%end_index(rl_end)

  CALL cells2edges_scalar_lib( p_cc, ptr_patch%edges%cell_idx, ptr_patch%edges%cell_blk, ptr_int%c_lin_e, opt_p_face, &
    &                          i_startblk_opt, i_endblk_opt, i_startidx_opt, i_endidx_opt, &
    &                          slev, elev, nproma, ptr_patch%id, l_limited_area, lfill_latbc=.FALSE., lacc=lacc)

ENDIF

CALL grad_green_gauss_cell_lib( p_cc, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, &
  &                             ptr_int%geofac_grg, p_grad, &
  &                             i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
  &                             slev, elev, nproma, ptr_patch%id, lacc=lacc )

END SUBROUTINE grad_green_gauss_cell_adv

SUBROUTINE grad_green_gauss_cell_dycore(p_ccpr, ptr_patch, ptr_int, p_grad, lacc,   &
    &                                   opt_slev, opt_elev, opt_rlstart, opt_rlend, &
    &                                   opt_acc_async)
  !
  !  patch on which computation is performed
  !
  TYPE(t_patch), TARGET, INTENT(in)     :: ptr_patch
  !
  !  data structure for interpolation
  !
  TYPE(t_int_state), TARGET, INTENT(in) :: ptr_int

  !  cell centered I/O variables
  !
  REAL(vp), INTENT(in) :: p_ccpr(:,:,:,:) ! perturbation fields passed from dycore (2,nproma,nlev,nblks_c)

  LOGICAL, INTENT(in) ::  &
    &  lacc    ! if TRUE, use OpenACC

  INTEGER, INTENT(in), OPTIONAL :: opt_slev    ! optional vertical start level

  INTEGER, INTENT(in), OPTIONAL :: opt_elev    ! optional vertical end level

  INTEGER, INTENT(in), OPTIONAL :: opt_rlstart, opt_rlend   ! start and end values of refin_ctrl flag

  LOGICAL, INTENT(IN), OPTIONAL :: opt_acc_async
  !
  ! cell based Green-Gauss reconstructed geographical gradient vector
  !
  REAL(vp), INTENT(inout) :: p_grad(:,:,:,:)      ! dim:(4,nproma,nlev,nblks_c)

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
#ifdef __SWAPDIM
    elev = UBOUND(p_ccpr,2)
#else
    elev = UBOUND(p_ccpr,3)
#endif
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

  i_startblk = ptr_patch%cells%start_block(rl_start)
  i_endblk   = ptr_patch%cells%end_block(rl_end)

  i_startidx_in = ptr_patch%cells%start_index(rl_start)
  i_endidx_in   = ptr_patch%cells%end_index(rl_end)

  CALL grad_green_gauss_cell_lib(p_ccpr, ptr_patch%cells%neighbor_idx, ptr_patch%cells%neighbor_blk, &
      &                          ptr_int%geofac_grg, p_grad, &
      &                          i_startblk, i_endblk, i_startidx_in, i_endidx_in, &
      &                          slev, elev, nproma, lacc=lacc, acc_async=opt_acc_async )

END SUBROUTINE grad_green_gauss_cell_dycore

END MODULE mo_math_gradients
