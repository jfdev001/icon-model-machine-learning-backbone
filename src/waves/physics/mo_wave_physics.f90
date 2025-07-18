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

! Contains the subroutines with wave physics parametrisation

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_physics

  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_parallel_config,     ONLY: nproma
  USE mo_physical_constants,  ONLY: grav
  USE mo_math_constants,      ONLY: dbl_eps, pi, pi2
  USE mo_fortran_tools,       ONLY: init
  USE mo_wave_types,          ONLY: t_wave_diag
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_wave_constants,      ONLY: EPS1, EMIN

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: air_sea
  PUBLIC :: last_prog_freq_ind
  PUBLIC :: impose_high_freq_tail
  PUBLIC :: tm1_tm2_periods_and_wm1_wm2_wavenumber
  PUBLIC :: wave_stress
  PUBLIC :: mean_frequency_and_total_energy
  PUBLIC :: compute_wave_number
  PUBLIC :: compute_group_velocity, wave_group_velocity_nt
  PUBLIC :: set_energy2emin
  PUBLIC :: mask_energy
  PUBLIC :: sdepth_lim

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_physics'

CONTAINS


  !>
  !! Calculation of group velocity.
  !!
  !! Wrapper routine for computing group velocity absolute values at
  !! cell centers and edge midpoints.
  !!
  SUBROUTINE compute_group_velocity(p_patch, wave_config, depth_c, depth_e, &
    &                               wave_num_c, wave_num_e, gv_c, gv_e)

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: wave_config
    REAL(wp),            INTENT(IN)   :: depth_c(:,:)      !< water depth at cell centers ( m )
    REAL(wp),            INTENT(IN)   :: depth_e(:,:)      !< water depth at edge midpoints ( m )
    REAL(wp),            INTENT(IN)   :: wave_num_c(:,:,:) !< wave number at cell center (1/m)
    REAL(wp),            INTENT(IN)   :: wave_num_e(:,:,:) !< wave number at edge midpoints (1/m)
    REAL(wp),            INTENT(INOUT):: gv_c(:,:,:)       !< group velocity at cell center (absolute value)  ( m/s )
    REAL(wp),            INTENT(INOUT):: gv_e(:,:,:)       !< group velocity edge midpoint (absolute value)  ( m/s )

    ! compute absolute value of group velocity at cell centers
    !
    CALL wave_group_velocity_c(           & !in
      &  p_patch     = p_patch,           & !in
      &  p_config    = wave_config,       & !in
      &  wave_num_c  = wave_num_c(:,:,:), & !in
      &  depth_c     = depth_c(:,:),      & !in
      &  gv_c        = gv_c(:,:,:))         !out

    ! compute absolute value of group velocity at edge midpoints
    !
    CALL wave_group_velocity_e(           & !in
      &  p_patch     = p_patch,           & !in
      &  p_config    = wave_config,       & !in
      &  wave_num_e  = wave_num_e(:,:,:), & !in
      &  depth_e     = depth_e(:,:),      & !in
      &  gv_e        = gv_e(:,:,:))         !out

  END SUBROUTINE compute_group_velocity


  !>
  !! Calculation of wave group velocity
  !!
  !! Calculation of shallow water cell centered
  !! wave group velocity
  !!
  SUBROUTINE wave_group_velocity_c(p_patch, p_config, wave_num_c, depth_c, gv_c)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_c'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(IN)   :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),            INTENT(IN)   :: depth_c(:,:)      !< water depth at cell centers (nproma,nblks_c) ( m )
    REAL(wp),            INTENT(INOUT):: gv_c(:,:,:)       !< group velocity (nproma,nfreqs,nblks_c)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,jc
    INTEGER :: nfreqs

    REAL(wp) :: gh, ak, akd, gv

    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    nfreqs = p_config%nfreqs

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jc,i_startidx,i_endidx,ak,akd,gv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1, nfreqs

        DO jc = i_startidx, i_endidx
          ! shallow water group velocity
          ak = wave_num_c(jc,jf,jb)
          akd = ak * depth_c(jc,jb)

          IF (akd <= 10.0_wp) THEN
            gv = 0.5_wp * SQRT(grav * TANH(akd)/ak) * (1.0_wp + 2.0_wp*akd/SINH(2.0_wp*akd))
          ELSE
            gv = gh / p_config%freqs(jf)
          END IF
          gv_c(jc,jf,jb) = gv
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_group_velocity_c


  !>
  !! Calculation of wave group velocity
  !!
  !! Calculation of shallow water edges centered
  !! wave group velocities
  !!
  SUBROUTINE wave_group_velocity_e(p_patch, p_config, wave_num_e, depth_e, gv_e)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_e'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    REAL(wp),            INTENT(IN)   :: wave_num_e(:,:,:) !< wave number (1/m)
    REAL(wp),            INTENT(IN)   :: depth_e(:,:)      !< water depth at cell edges (nproma,nblks_e) ( m )
    REAL(wp),            INTENT(INOUT):: gv_e(:,:,:)       !< group velocity (nproma,nfreqs,nblks_c)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,je
    INTEGER :: nfreqs

    REAL(wp) :: gh, ak, akd, gv

    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rledge
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    nfreqs = p_config%nfreqs

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,je,i_startidx,i_endidx,ak,akd,gv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1, nfreqs

        DO je = i_startidx, i_endidx
          ! shallow water group velocity
          ak = wave_num_e(je,jf,jb)
          akd = ak * depth_e(je,jb)

          IF (akd <= 10.0_wp) THEN
            gv = 0.5_wp * SQRT(grav * TANH(akd)/ak) * (1.0_wp + 2.0_wp*akd/SINH(2.0_wp*akd))
          ELSE
            gv = gh / p_config%freqs(jf)
          END IF
          gv_e(je,jf,jb) = gv
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE wave_group_velocity_e


  !>
  !! Calculation of wave group velocity
  !!
  !! Calculation of shallow water
  !! edge-normal and -tangential projections of
  !! wave group velocities using spectral directions
  !!
  SUBROUTINE wave_group_velocity_nt(p_patch, p_config, jf, gv_e, gvn_e, gvt_e)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':wave_group_velocity_nt'

    TYPE(t_patch),       INTENT(IN)   :: p_patch
    TYPE(t_wave_config), INTENT(IN)   :: p_config
    INTEGER ,            INTENT(IN)   :: jf           !< frequency index
    REAL(wp),            INTENT(IN)   :: gv_e(:,:,:)  !< group velocity (nproma,nfreqs,nblks_e)  ( m/s )
    REAL(wp),            INTENT(INOUT):: gvn_e(:,:,:) !< normal group velocity (nproma,ndirs,nblks_e)  ( m/s )
    REAL(wp),            INTENT(INOUT):: gvt_e(:,:,:) !< tangential group velocity (nproma,ndirs,nblks_e)  ( m/s )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jd,je,ic,jje,jjb
    INTEGER :: nfreqs, ndirs
    LOGICAL :: is_towards_coastline

    REAL(wp) :: gvu, gvv, gv

    i_rlstart  = 1
    i_rlend    = min_rledge
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    ndirs  = p_config%ndirs
    gv = grav / (2.0_wp * pi2 * p_config%freqs(jf))

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jd,je,i_startidx,i_endidx,gvu,gvv) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jd = 1, ndirs
        DO je = i_startidx, i_endidx
          gvu = gv_e(je,jf,jb) * p_config%sin_dir(jd)
          gvv = gv_e(je,jf,jb) * p_config%cos_dir(jd)

          gvn_e(je,jd,jb) = &
               gvu * p_patch%edges%primal_normal(je,jb)%v1 + &
               gvv * p_patch%edges%primal_normal(je,jb)%v2

          gvt_e(je,jd,jb) = &
               gvu * p_patch%edges%dual_normal(je,jb)%v1 + &
               gvv * p_patch%edges%dual_normal(je,jb)%v2
        END DO
      END DO
    END DO
!$OMP ENDDO

    ! Correction of normal to edge group velocity,
    ! avoiding of wave energy propagation from land
    ! and insuring full "outflow" of wave energy towards land.

    ! Set the wave group velocity to zero at the boundary edge
    ! in case of wave energy propagation towards the ocean,
    ! and set gn = deep water group velocity otherwise.

    ! We make use of the fact that the edge-normal velocity vector points
    ! * towards the coast, if
    !   cells%edge_orientation > 0 .AND. vn > 0
    !   OR
    !   cells%edge_orientation < 0 .AND. vn < 0
    ! * towards the sea, if
    !   cells%edge_orientation > 0 .AND. vn < 0
    !   OR
    !   cells%edge_orientation < 0 .AND. vn > 0

    ! hence:
    ! * towards the coast, if (cells%edge_orientation * vn) > 0
    ! * towards the sea,   if (cells%edge_orientation * vn) < 0

    ! For optimization, the index list of the coastal edge points is precomputed,
    ! and p_config%orient_coastedges contains the values of cells%edge_orientation

!$OMP DO PRIVATE(jd,ic,jje,jjb,is_towards_coastline) ICON_OMP_DEFAULT_SCHEDULE
!NEC$ outerloop_unroll(4)
    DO jd = 1, ndirs
!$NEC ivdep
      DO ic = 1, p_config%n_coastedges
        jje = p_config%idx_coastedges(ic)
        jjb = p_config%blk_coastedges(ic)
        is_towards_coastline = (p_config%orient_coastedges(ic)*gvn_e(jje,jd,jjb)) > 0._wp
        gvn_e(jje,jd,jjb) = MERGE(gv, 0.0_wp, is_towards_coastline)
      ENDDO  !jc
    ENDDO  !jd
!$OMP ENDDO NOWAIT

!$OMP END PARALLEL
  END SUBROUTINE wave_group_velocity_nt


  !>
  !! Calculation of total stress and sea surface roughness
  !!
  !! Adaptation of WAM 4.5 algorithm and code for calculation of total stress
  !! and sea surface roughness for ICON-waves (P.A.E.M. Janssen, 1990).
  !!
  SUBROUTINE air_sea(p_patch, wave_config, wsp10m, tauw, ustar, z0)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         & routine =  modname//'air_sea'

    TYPE(t_patch),      INTENT(IN)  :: p_patch
    TYPE(t_wave_config),INTENT(IN)  :: wave_config
    REAL(wp),           INTENT(IN)  :: wsp10m(:,:)!10m wind speed (nproma,nblks_c) ( m/s )
    REAL(wp),           INTENT(IN)  :: tauw(:,:)  !wave stress (nproma,nblks_c) ( (m/s)^2 )
    REAL(wp),           INTENT(OUT) :: ustar(:,:) !friction velocity (nproma,nblks_c) ( m/s )
    REAL(wp),           INTENT(OUT) :: z0(:,:)    !roughness length (nproma,nblks_c) ( m )

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,iter

    REAL(wp), PARAMETER :: TWOXMP1 = 3.0_wp
    REAL(wp), PARAMETER :: EPSUS   = 1.0E-6_wp

    !     *ACD*       COEFFICIENTS FOR SIMPLE CD(U10) RELATION
    !     *BCD*       CD = ACD + BCD*U10
    REAL(wp), PARAMETER :: ACD = 8.0E-4_wp
    REAL(wp), PARAMETER :: BCD = 8.0E-5_wp

    INTEGER, PARAMETER :: NITER = 15

    REAL(wp):: xkutop(nproma), tauold(nproma), ustm1(nproma), z0ch(nproma)
    REAL(wp):: xlogxl, alphaog, xologz0
    REAL(wp):: ustold, taunew, x, f, delf
    REAL(wp):: z0tot, z0vis, zz
    LOGICAL :: l_converged(nproma)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    xlogxl  = LOG(wave_config%XNLEV)
    alphaog = wave_config%ALPHA_CH / grav

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,iter,i_startidx,i_endidx,xkutop,ustold,tauold,ustm1, &
!$OMP            l_converged,x,z0ch,z0vis,z0tot,xologz0,f,zz,delf,taunew) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! initialization
      DO jc = i_startidx, i_endidx
        xkutop(jc) = wave_config%XKAPPA * wsp10m(jc,jb)
        ustold     = wsp10m(jc,jb) * SQRT(ACD + BCD * wsp10m(jc,jb))
        tauold(jc) = MAX(ustold**2,tauw(jc,jb)+EPS1)
        ustar(jc,jb) = SQRT(tauold(jc))
        ustm1(jc) = 1.0_wp/MAX(ustar(jc,jb),EPSUS)

        l_converged(jc) = .FALSE.
      END DO

      DO iter = 1,NITER
        DO jc = i_startidx, i_endidx
          IF (l_converged(jc)) CYCLE
          x        = tauw(jc,jb) / tauold(jc)
          z0ch(jc) = alphaog * tauold(jc) / SQRT(MAX(1.0_wp-x,EPS1))
          z0vis    = wave_config%RNUAIRM * ustm1(jc)
          z0tot    = z0ch(jc) + z0vis

          xologz0 = 1.0 / (xlogxl - LOG(z0tot))
          f = ustar(jc,jb) - xkutop(jc) * xologz0
          zz = ustm1(jc) &
               * (z0ch(jc) * (2.0_wp-TWOXMP1*X) / (1.0_wp -x) - z0vis) &
               / z0tot
          delf = 1.0_wp - xkutop(jc) * xologz0**2 * zz

          ustar(jc,jb) = ustar(jc,jb) - f / delf

          taunew = MAX(ustar(jc,jb)**2,tauw(jc,jb) + EPS1)

          ustar(jc,jb) = SQRT(taunew)

          IF (ABS(taunew-tauold(jc))<= dbl_eps) l_converged(jc) = .TRUE.

          ustm1(jc)  = 1.0_wp/MAX(ustar(jc,jb),EPSUS)
          tauold(jc) = taunew
        END DO  !jc
        !
        IF ( ALL(l_converged(i_startidx:i_endidx)) ) EXIT

      END DO  !iter
      !
      DO jc = i_startidx, i_endidx
        z0(jc,jb) = z0ch(jc)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE air_sea


  !>
  !! Calculation of mean frequency and total energy
  !!
  !! Integration over frequencies for calculation of mean frequency
  !! energy and total energy (combined for efficiency). Adaptation of WAM 4.5 code
  !! of the subroutines FEMEAN and TOTAL_ENERGY developed by S.D. HASSELMANN,
  !! optimized by L. Zambresky and H. Guenther, GKSS, 2001                              !
  !!
  SUBROUTINE mean_frequency_and_total_energy(p_patch, wave_config, tracer, llws, emean, emeanws, femean, femeanws)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'mean_frequency_and_total_energy'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp), INTENT(IN)    :: tracer(:,:,:,:) !energy spectral bins (nproma,ndirs,nblks_c,nfreqs)
    INTEGER,  INTENT(IN)    :: llws(:,:,:,:)   !=1 where wind_input is positive (nproma,ndirs,nblks_c,nfreqs)
    REAL(wp), INTENT(INOUT) :: emean(:,:)      !total energy (nproma,nblks_c)
    REAL(wp), INTENT(INOUT) :: emeanws(:,:)    !total windsea energy (nproma,nblks_c)
    REAL(wp), INTENT(INOUT) :: femean(:,:)     !mean frequency energy (nproma,nblks_c)
    REAL(wp), INTENT(INOUT) :: femeanws(:,:)   !mean windsea frequency energy (nproma,nblks_c)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd

    REAL(wp) :: temp(nproma,wave_config%nfreqs), temp_1(nproma,wave_config%nfreqs)
    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jd,jf,i_startidx,i_endidx,temp,temp_1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          temp(jc,jf)   = 0._wp
          temp_1(jc,jf) = 0._wp
        ENDDO

!nec$ outerloop_unroll(4)
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            temp(jc,jf)   = temp(jc,jf)   + tracer(jc,jd,jb,jf)
            temp_1(jc,jf) = temp_1(jc,jf) + tracer(jc,jd,jb,jf)*REAL(llws(jc,jd,jb,jf),wp)
          ENDDO
        ENDDO  ! jd
      END DO  ! jf

      DO jc = i_startidx, i_endidx
        emean(jc,jb)    = wc%MO_TAIL  * temp(jc,wc%nfreqs)
        emeanws(jc,jb)  = wc%MO_TAIL  * temp_1(jc,wc%nfreqs)
        femean(jc,jb)   = wc%MM1_TAIL * temp(jc,wc%nfreqs)
        femeanws(jc,jb) = wc%MM1_TAIL * temp_1(jc,wc%nfreqs)
      END DO

!nec$ outerloop_unroll(4)
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          emean(jc,jb)    = emean(jc,jb)    + temp(jc,jf)   * wc%DFIM(jf)
          emeanws(jc,jb)  = emeanws(jc,jb)  + temp_1(jc,jf) * wc%DFIM(jf)
          femean(jc,jb)   = femean(jc,jb)   + temp(jc,jf)   * wc%DFIMOFR(jf)
          femeanws(jc,jb) = femeanws(jc,jb) + temp_1(jc,jf) * wc%DFIMOFR(jf)
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        emean(jc,jb)    = MAX(emean(jc,jb),EMIN)
        emeanws(jc,jb)  = MAX(emeanws(jc,jb),EMIN)
        femean(jc,jb)   = emean(jc,jb) / MAX(femean(jc,jb),EMIN)
        femeanws(jc,jb) = emeanws(jc,jb) / MAX(femeanws(jc,jb),EMIN)
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE mean_frequency_and_total_energy


  !>
  !! Calculation of wave stress.
  !!
  !! Compute normalized wave stress from input source function
  !!
  !! Adaptation of WAM 4.5 code.
  !! STRESSO
  !!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT
  !!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990
  !!     J. BIDLOT             ECMWF FEBRUARY  1996-97
  !!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90
  !!     J. BIDLOT             ECMWF           2007  ADD MIJ
  !!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS
  !!
  !! Reference
  !!       R SNYDER ET AL,1981.
  !!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
  !!       P. JANSSEN, JPO, 1985
  !!
  SUBROUTINE wave_stress(p_patch, wave_config, dir10m, sl, tracer, p_diag)
     CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
          &  routine = modname//'wave_stress'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: dir10m(:,:)
    REAL(vp),                    INTENT(IN)    :: sl(:,:,:,:)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)
    TYPE(t_wave_diag),           INTENT(INOUT) :: p_diag

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd
    INTEGER :: jf_lp     ! last prognostic frequency index

    REAL(wp) :: gm1, const, sinplus, cosw
    REAL(wp) :: cmrhowgdfth
    REAL(wp) :: const1(nproma), const2(nproma)
    REAL(wp) :: rhowgdfth(nproma,wave_config%nfreqs)
    REAL(wp) :: cm(nproma,wave_config%nfreqs)
    REAL(wp) :: xstress(nproma), xstress_tot
    REAL(wp) :: ystress(nproma), ystress_tot
    REAL(wp) :: temp1(nproma), temp2(nproma)
    REAL(wp) :: sumt(nproma), sumx(nproma), sumy(nproma)
    REAL(wp) :: roair    ! air density

    TYPE(t_wave_config), POINTER :: wc => NULL()

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config

    gm1   = 1.0_wp/grav
    const = wc%delth*(pi2)**4*gm1
    roair = MAX(wc%roair,1._wp)


!$OMP PARALLEL
    CALL init(p_diag%phiaw, lacc=.FALSE.)
!$OMP BARRIER
!$OMP DO PRIVATE(jb,jc,jf,jd,jf_lp,i_startidx,i_endidx,cm,const1,const2,       &
!$OMP            rhowgdfth,sinplus,sumt,sumx,sumy,cmrhowgdfth,xstress,ystress, &
!$OMP            xstress_tot,ystress_tot,cosw,temp1,temp2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          cm(jc,jf) = p_diag%wave_num_c(jc,jf,jb) * 1.0_wp/(pi2*wc%freqs(jf))
          rhowgdfth(jc,jf) = MERGE(wc%rhowg_dfim(jf), 0.0_wp, jf <= p_diag%last_prog_freq_ind(jc,jb))
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx

        jf = p_diag%last_prog_freq_ind(jc,jb)
        IF (jf /= wc%nfreqs) rhowgdfth(jc,jf) = 0.5_wp * rhowgdfth(jc,jf)

        !initialisation
        xstress(jc) = 0._wp
        ystress(jc) = 0._wp

      END DO

      !sum
      DO jf = 1, MAXVAL(p_diag%last_prog_freq_ind(i_startidx:i_endidx,jb))
        DO jc = i_startidx, i_endidx
          sumt(jc) = 0._wp
          sumx(jc) = 0._wp
          sumy(jc) = 0._wp
        END DO

        DO jd = 1, wc%ndirs
          DO jc = i_startidx, i_endidx
            sinplus = MAX(sl(jc,jd,jb,jf),0._wp)
            sumt(jc) = sumt(jc) + sinplus
            sumx(jc) = sumx(jc) + sinplus * wc%sin_dir(jd)
            sumy(jc) = sumy(jc) + sinplus * wc%cos_dir(jd)
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          p_diag%phiaw(jc,jb) =  p_diag%phiaw(jc,jb) + sumt(jc)*rhowgdfth(jc,jf)
          cmrhowgdfth = cm(jc,jf) * rhowgdfth(jc,jf)
          xstress(jc) = xstress(jc) + sumx(jc)*cmrhowgdfth
          ystress(jc) = ystress(jc) + sumy(jc)*cmrhowgdfth
        END DO

      END DO  ! jf



      ! calculate high-frequency contribution to stress
      !
      DO jc = i_startidx, i_endidx
        temp1(jc)  = 0._wp
        temp2(jc)  = 0._wp
        const1(jc) = const * wc%freqs(p_diag%last_prog_freq_ind(jc,jb))**5 * gm1
        const2(jc) = roair * const * wc%freqs(p_diag%last_prog_freq_ind(jc,jb))**5
      ENDDO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          jf_lp = p_diag%last_prog_freq_ind(jc,jb)

          cosw = MAX(COS(wc%dirs(jd)-dir10m(jc,jb)),0.0_wp)
          temp1(jc) = temp1(jc) + tracer(jc,jd,jb,jf_lp) * cosw**3
          temp2(jc) = temp2(jc) + tracer(jc,jd,jb,jf_lp) * cosw**2
        END DO
      END DO

      CALL high_frequency_stress(wave_config        = wave_config,                     & !IN
        &                        i_startidx         = i_startidx,                      & !IN
        &                        i_endidx           = i_endidx,                        & !IN
        &                        last_prog_freq_ind = p_diag%last_prog_freq_ind(:,jb), & !IN
        &                        ustar              = p_diag%ustar(:,jb),              & !IN
        &                        z0                 = p_diag%z0(:,jb),                 & !IN
        &                        xlevtail           = p_diag%xlevtail(:,jb),           & !IN
        &                        tauhf1             = p_diag%tauhf1(:,jb),             & !INOUT
        &                        phihf1             = p_diag%phihf1(:,jb) )              !INOUT

      DO jc = i_startidx, i_endidx
        p_diag%tauhf(jc,jb) = const1(jc)*temp1(jc)*p_diag%tauhf1(jc,jb)
        p_diag%phihf(jc,jb) = const2(jc)*temp2(jc)*p_diag%phihf1(jc,jb)

        p_diag%phiaw(jc,jb) = p_diag%phiaw(jc,jb) + p_diag%phihf(jc,jb)

        xstress_tot = xstress(jc)/roair + p_diag%tauhf(jc,jb)*SIN(dir10m(jc,jb))
        ystress_tot = ystress(jc)/roair + p_diag%tauhf(jc,jb)*COS(dir10m(jc,jb))

        p_diag%tauw(jc,jb) = SQRT(xstress_tot**2+ystress_tot**2)
        p_diag%tauw(jc,jb) = MIN(p_diag%tauw(jc,jb),p_diag%ustar(jc,jb)**2 - EPS1)
        p_diag%tauw(jc,jb) = MAX(p_diag%tauw(jc,jb),0.0_wp)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE wave_stress


  !>
  !! Calculation of high frequency stress.
  !!
  !! Adaptation of WAM 4.5 code.
  !!
  SUBROUTINE high_frequency_stress(wave_config, i_startidx, i_endidx, last_prog_freq_ind, &
    &                              ustar, z0, xlevtail, tauhf1, phihf1)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//'high_frequency_stress'

    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    INTEGER,                     INTENT(IN)    :: i_startidx, i_endidx
    INTEGER,                     INTENT(IN)    :: last_prog_freq_ind(:) !<last frequency index of the prognostic range
    REAL(wp),                    INTENT(IN)    :: ustar(:)              !< friction velocity
    REAL(wp),                    INTENT(IN)    :: z0(:)                 !< roughness length
    REAL(wp),                    INTENT(IN)    :: xlevtail(:)           !< tail level
    REAL(wp),                    INTENT(INOUT) :: tauhf1(:)             !< high-frequency stress
    REAL(wp),                    INTENT(INOUT) :: phihf1(:)             !< high-frequency energy flux into ocean

    INTEGER :: jc, j
    REAL(wp):: gm1, x0g
    REAL(wp):: OMEGA, OMEGAC, OMEGACC
    REAL(wp):: UST(nproma), UST0, TAUW0, tauw(nproma)
    REAL(wp):: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
    REAL(wp):: DELZ(nproma), ZINF(nproma)
    REAL(wp):: FNC2, SQRTZ0OG(nproma), GZ0, SQRTGZ0(nproma), XLOGGZ0(nproma)

    TYPE(t_wave_config), POINTER :: wc => NULL()
    REAL(wp), PARAMETER :: ZSUP = 0.0_wp  !  LOG(1.)

    wc => wave_config

    gm1 = 1.0_wp/grav
    x0g = wc%X0TAUHF * grav

    DO jc = i_startidx, i_endidx
      OMEGAC      = pi2 * wc%freqs(last_prog_freq_ind(jc))
      UST0        = ustar(jc)
      TAUW0       = UST0**2.0_wp
      GZ0         = grav * z0(jc)
      OMEGACC     = MAX(OMEGAC,X0G/UST0)

      XLOGGZ0(jc) = LOG(GZ0)
      SQRTZ0OG(jc)= SQRT(z0(jc)*GM1)
      SQRTGZ0(jc) = 1.0_wp / SQRTZ0OG(jc)
      YC          = OMEGACC * SQRTZ0OG(jc)
      ZINF(jc)    = LOG(YC)
      DELZ(jc)    = MAX((ZSUP-ZINF(jc))/REAL(wc%JTOT_TAUHF-1),0.0_wp)

      tauw(jc)    = TAUW0
      UST(jc)     = UST0

      tauhf1(jc)  = 0.0_wp
      phihf1(jc)  = 0.0_wp
    ENDDO


    ! Integrals are integrated following a change of variable : Z=LOG(Y)
    DO J = 1, wc%JTOT_TAUHF
      DO jc = i_startidx, i_endidx
        Y         = EXP(ZINF(jc)+REAL(J-1)*DELZ(jc))
        OMEGA     = Y * SQRTGZ0(jc)
        CM1       = OMEGA * GM1
        ZX        = UST(jc) * CM1 + wc%zalp
        ZARG      = wc%XKAPPA/ZX
        ZLOG      = XLOGGZ0(jc) + 2.0_wp * LOG(CM1)+ZARG
        ZLOG      = MIN(ZLOG,0.0_wp)
        ZBETA     = EXP(ZLOG) * ZLOG**4

        FNC2      = ZBETA * tauw(jc) * wc%WTAUHF(J) * DELZ(jc)
        tauw(jc)  = MAX(tauw(jc)-xlevtail(jc) * FNC2,0.0_wp)
        UST(jc)   = SQRT(tauw(jc))

        tauhf1(jc) = tauhf1(jc) + FNC2
        phihf1(jc) = phihf1(jc) + FNC2/Y
      END DO
    END DO

    DO jc = i_startidx, i_endidx
      phihf1(jc) = SQRTZ0OG(jc) * phihf1(jc)
    END DO

  END SUBROUTINE high_frequency_stress


  !>
  !! Returns the last frequency index of prognostic part of spectrum
  !!
  !! Compute last frequency index of prognostic part of spectrum.
  !! Frequencies le MAX(tailfactor*max(fmnws,fm),tailfactor_pm*fpm),
  !! where fpm is the Pierson-Moskowitz frequency based on friction
  !! velocity. (fpm=g/(fric*zpi*ustar))
  !!
  !! Adaptation of WAM 4.5 code.
  !! FRCUTINDEX
  !!
  !! Initial revision by Mikhail Dobrynin, DWD (2019-10-10)
  !!
  SUBROUTINE last_prog_freq_ind(p_patch, wave_config, femeanws, femean, ustar, lpfi)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//'last_prog_freq_ind'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: femeanws(:,:) !< mean frequency wind sea wave energy
    REAL(wp),                    INTENT(IN)    :: femean(:,:)   !< mean frequency wave energy
    REAL(wp),                    INTENT(IN)    :: ustar(:,:)    !< friction velocity
    INTEGER,                     INTENT(INOUT) :: lpfi(:,:)     !< last frequency index
                                                                !  of prognostic part of spectrum
    ! local
    TYPE(t_wave_config), POINTER :: wc => NULL()
    !
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc
    !
    REAL(wp) :: fpmh, fppm, fm2, fpm, fpm4, inv_log_co
    !
    REAL(wp), PARAMETER :: epsus = 1.0e-6_wp
    REAL(wp), PARAMETER :: fric = 28.0_wp
    REAL(wp), PARAMETER :: tailfactor = 2.5_wp
    REAL(wp), PARAMETER :: tailfactor_pm = 3.0_wp

    wc => wave_config

    fpmh = tailfactor / wc%freqs(1)
    fppm = tailfactor_pm * grav / (fric * pi2 * wc%freqs(1))
    inv_log_co = 1.0_wp / LOG10(wc%co)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,i_startidx,i_endidx,fm2,fpm,fpm4) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        fm2 = MAX(femeanws(jc,jb),femean(jc,jb)) * fpmh
        fpm = fppm / MAX(ustar(jc,jb),epsus)
        fpm4 = MAX(fm2,fpm)
        lpfi(jc,jb) = NINT(LOG10(fpm4)*inv_log_co)+1
        lpfi(jc,jb) = MIN(MAX(1,lpfi(jc,jb)),wc%nfreqs)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE last_prog_freq_ind


  !>
  !! Impose high frequency tail to the spectrum
  !!
  !! Adaptation of WAM 4.5 code.
  !! IMPHFTAIL
  !!
  SUBROUTINE impose_high_freq_tail(p_patch, wave_config, wave_num_c, depth, last_prog_freq_ind, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'impose_high_freq_tail'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:)  !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    INTEGER,                     INTENT(IN)    :: last_prog_freq_ind(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)


    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jf,jc,jd
    INTEGER :: jf_lp     ! last prognostic frequency index

    REAL(wp) :: gh, ak, akd, tcgond, akm1
    REAL(wp) :: temp(nproma, wave_config%nfreqs)
    REAL(wp) :: tfac(nproma)

    wc => wave_config
    gh = grav / (4.0_wp * pi)

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jc,jd,jf_lp,i_startidx,i_endidx,ak,akd,tcgond,akm1,temp,tfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jf = 1, wc%nfreqs
        DO jc = i_startidx, i_endidx
          ak = wave_num_c(jc,jf,jb)
          IF (jf >=last_prog_freq_ind(jc,jb)) THEN
            akd = ak * depth(jc,jb)
            IF (akd.le.10.0_wp) THEN
              tcgond = 0.5_wp * SQRT(grav * TANH(akd)/ak) * (1.0_wp + 2.0_wp*akd/SINH(2.0_wp*akd))
            ELSE
              tcgond = gh / wc%freqs(jf)
            END IF
            akm1 = 1.0_wp/ak
            temp(jc,jf) = akm1**3/tcgond
          END IF
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        DO jf = last_prog_freq_ind(jc,jb)+1, wc%nfreqs
          temp(jc,jf) = temp(jc,jf) / temp(jc,last_prog_freq_ind(jc,jb))
        END DO
      END DO

      DO jd = 1, wc%ndirs
        DO jc = i_startidx, i_endidx
          jf_lp = last_prog_freq_ind(jc,jb)
          tfac(jc) = tracer(jc,jd,jb,jf_lp)
        ENDDO
        !
        DO jf = 1, wc%nfreqs
          DO jc = i_startidx, i_endidx
            IF (jf >=last_prog_freq_ind(jc,jb)+1) THEN
              tracer(jc,jd,jb,jf) = temp(jc,jf) * tfac(jc)
            END IF
          END DO  !jc
        END DO  !jf
      END DO  !jd

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE impose_high_freq_tail


  !>
  !! Calculation of TM1 and TM2 periods and WM1 and WM2 wavenumbers
  !!
  !! Adaptation of WAM 4.5 code.
  !! TM1_TM2_PERIODS_B and WM1_WM2_WAVENUMBER_B (combined for efficiency)
  !!
  !! C.Schneggenburger 08/97.
  !!
  !! Integration of spectra and adding of tail factors.
  !!   WM1 IS SQRT(1/K)*F
  !!   WM2 IS SQRT(K)*F
  !!
  SUBROUTINE tm1_tm2_periods_and_wm1_wm2_wavenumber(p_patch, wave_config, wave_num_c, tracer, emean, &
    &                                               tm1, tm2, f1mean, akmean, xkmean)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         & routine =  modname//'tm1_tm2_periods_and_wm1_wm2_wavenumber'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: wave_num_c(:,:,:) !< wave number (1/m)
    REAL(wp),                    INTENT(IN)    :: tracer(:,:,:,:)
    REAL(wp),                    INTENT(IN)    :: emean(:,:)    !< total wave energy
    REAL(wp),                    INTENT(INOUT) :: tm1(:,:)        !< tm1 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: tm2(:,:)        !< tm2 period (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: f1mean(:,:)     !< tm1 frequency (nproma,nblks_c)
    REAL(wp),                    INTENT(INOUT) :: akmean(:,:)   !< mean wavenumber based on SQRT(1/K)-moment, wm1
    REAL(wp),                    INTENT(INOUT) :: xkmean(:,:)   !< mean wavenumber based on SQRT(K)-moment, wm2

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd

    REAL(wp) :: temp(nproma, wave_config%nfreqs)
    REAL(wp) :: temp2(nproma, wave_config%nfreqs)

    wc => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,jd,i_startidx,i_endidx,temp,temp2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs

        DO jc = i_startidx, i_endidx
          temp(jc,jf) = 0._wp
          temp2(jc,jf)= SQRT(wave_num_c(jc,jf,jb))
        ENDDO
        !
        ! sum
!nec$ outerloop_unroll(4)
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            temp(jc,jf) = temp(jc,jf) + tracer(jc,jd,jb,jf)
          END DO
        ENDDO
      END DO  !jf

      !initialisation / tail part
      DO jc = i_startidx, i_endidx
        tm1(jc,jb)    = wc%MP1_TAIL * temp(jc,wc%nfreqs)
        tm2(jc,jb)    = wc%MP2_TAIL * temp(jc,wc%nfreqs)
        akmean(jc,jb) = wc%MM1_TAIL * SQRT(grav)/pi2 * temp(jc,wc%nfreqs)
        xkmean(jc,jb) = wc%MM1_TAIL * SQRT(grav)/pi2 * temp(jc,wc%nfreqs)
      END DO

      ! add all other frequencies
      DO jf = 1,wc%nfreqs
        DO jc = i_startidx, i_endidx
          tm1(jc,jb)    = tm1(jc,jb) + temp(jc,jf) * wc%dfim_fr(jf)
          tm2(jc,jb)    = tm2(jc,jb) + temp(jc,jf) * wc%dfim_fr2(jf)
          akmean(jc,jb) = akmean(jc,jb) + temp(jc,jf) / temp2(jc,jf) * wc%DFIM(jf)
          xkmean(jc,jb) = xkmean(jc,jb) + temp(jc,jf) * temp2(jc,jf) * wc%DFIM(jf)
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        IF (emean(jc,jb) > EMIN) THEN
          tm1(jc,jb)    = emean(jc,jb) / tm1(jc,jb)
          tm2(jc,jb)    = SQRT(emean(jc,jb) / tm2(jc,jb))
          akmean(jc,jb) = ( emean(jc,jb) / akmean(jc,jb))**2
          xkmean(jc,jb) = ( xkmean(jc,jb) / emean(jc,jb))**2
        ELSE
          tm1(jc,jb)    = 1._wp
          tm2(jc,jb)    = 1._wp
          akmean(jc,jb) = 1._wp
          xkmean(jc,jb) = 1._wp
        END IF
        f1mean(jc,jb) = 1.0_wp / tm1(jc,jb)
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE tm1_tm2_periods_and_wm1_wm2_wavenumber


  !>
  !! Calculation of wave number at centers and edges.
  !!
  !! Wrapper routine for computing wave number at cell centers and
  !! edge midpoints.
  !!
  SUBROUTINE compute_wave_number(p_patch, wave_config, depth_c, depth_e, wave_num_c, wave_num_e)

    TYPE(t_patch),       INTENT(IN)    :: p_patch
    TYPE(t_wave_config), INTENT(IN)    :: wave_config
    REAL(wp),            INTENT(IN)    :: depth_c(:,:)      !< bathymetric height at cell centers
    REAL(wp),            INTENT(IN)    :: depth_e(:,:)      !< bathymetric height at edge midpoints
    REAL(wp),            INTENT(INOUT) :: wave_num_c(:,:,:) !< Wave number as a function of circular frequency
    REAL(wp),            INTENT(INOUT) :: wave_num_e(:,:,:) !< Wave number as a function of circular frequency

    ! get wave number as a function of circular frequency and water depth
    ! at cell center
    CALL wave_number_c(p_patch     = p_patch,              & !IN
      &                wave_config = wave_config,          & !IN
      &                depth       = depth_c(:,:),         & !IN
      &                wave_num_c  = wave_num_c(:,:,:))      !OUT

    ! get wave number as a function of circular frequency and water depth
    ! at edge midpoint
    CALL wave_number_e(p_patch     = p_patch,             & !IN
      &                wave_config = wave_config,         & !IN
      &                depth       = depth_e(:,:),        & !IN
      &                wave_num_e  = wave_num_e(:,:,:))     !OUT

  END SUBROUTINE compute_wave_number


  !>
  !! Calculation of wave number at cell centers.
  !!
  !! Wave number as a function of circular frequency and water depth.
  !! Newtons method to solve the dispersion relation in shallow water.
  !! G. KOMEN, P. JANSSEN   KNMI              01/06/1986
  !! Adaptation of WAM 4.5 code, function AKI
  !!
  SUBROUTINE wave_number_c(p_patch, wave_config, depth, wave_num_c)

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: depth(:,:)        !< bathymetric height at cell centers (nproma)
    REAL(wp),                    INTENT(INOUT) :: wave_num_c(:,:,:) !< Wave number as a function of circular frequency
    ! local
    REAL(wp), PARAMETER  :: EBS = 0.0001_wp     !< RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
    REAL(wp), PARAMETER  :: DKMAX = 40.0_wp     !< MAXIMUM VALUE OF DEPTH*WAVENUMBER.
    REAL(wp)             :: BO, TH, STH
    INTEGER              :: jb, jc, jf          !< loop indices
    REAL(wp)             :: OM                  !< CIRCULAR FREQUENCY 2*pi*freq (nfreqs)
    LOGICAL              :: l_converged(nproma), all_converged
    REAL(wp)             :: wave_num(nproma)
    REAL(wp)             :: AKP
    INTEGER              :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER              :: i_startidx, i_endidx
    TYPE(t_wave_config), POINTER :: wc => NULL()
    !  ---------------------------------------------------------------------------- !
    !     1. START WITH MAXIMUM FROM DEEP AND EXTREME SHALLOW WATER WAVE NUMBER.    !

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    wc  => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jf,i_startidx,i_endidx,OM,wave_num,l_converged,all_converged,BO,AKP,TH,STH) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        OM = pi2 * wc%freqs(jf)
        AKP = 10000.0_wp
        DO jc=i_startidx, i_endidx
          !
          ! initialization
          wave_num(jc) = MAX( OM**2/(4.0_wp* grav), OM/(2.0_wp*SQRT(grav*depth(jc,jb))) )

          l_converged(jc) = (ABS(AKP-wave_num(jc)) <= EBS*wave_num(jc))
        ENDDO  ! jc
        all_converged =  ALL(l_converged(i_startidx:i_endidx))

        ! ---------------------------------------------------------------------------- !
        !     2. ITERATION LOOP.                                                       !

        DO WHILE (.NOT. all_converged)
          DO jc = i_startidx, i_endidx
            IF (.NOT.l_converged(jc)) THEN
              BO = depth(jc,jb)*wave_num(jc)
              IF (BO > DKMAX) THEN
                wave_num(jc) = OM**2/grav
                l_converged(jc) = .TRUE.
              ELSE
                AKP = wave_num(jc)
                TH  = grav*wave_num(jc)*TANH(BO)
                STH = SQRT(TH)
                wave_num(jc) = wave_num(jc) &
                  &                + (OM-STH)*STH*2.0_wp / (TH/wave_num(jc) + grav*BO/COSH(BO)**2)
                ! check for converged solution
                l_converged(jc) = (ABS(AKP-wave_num(jc)) <= EBS*wave_num(jc))
              END IF
            END IF  ! l_converged
          ENDDO  !jc
          all_converged = ALL(l_converged(i_startidx:i_endidx))
        ENDDO !while
        !
        DO jc = i_startidx, i_endidx
          wave_num_c(jc,jf,jb) = wave_num(jc)
        ENDDO
        !
      ENDDO ! jf
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE wave_number_c


  !>
  !! Calculation of wave number at edge midpoints.
  !!
  !! Wave number as a function of circular frequency and water depth.
  !! Newtons method to solve the dispersion relation in shallow water.
  !! G. KOMEN, P. JANSSEN   KNMI              01/06/1986
  !! Adaptation of WAM 4.5 code, function AKI
  !!
  SUBROUTINE wave_number_e(p_patch, wave_config, depth, wave_num_e)

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: depth(:,:)        !< bathymetric height at cell centers (nproma)
    REAL(wp),                    INTENT(INOUT) :: wave_num_e(:,:,:) !< Wave number as a function of circular frequency
    ! local
    REAL(wp), PARAMETER  :: EBS = 0.0001_wp     !< RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
    REAL(wp), PARAMETER  :: DKMAX = 40.0_wp     !< MAXIMUM VALUE OF DEPTH*WAVENUMBER.
    REAL(wp)             :: BO, TH, STH
    INTEGER              :: jb, je, jf          !< loop indices
    REAL(wp)             :: OM                  !< CIRCULAR FREQUENCY 2*pi*freq (nfreqs)
    LOGICAL              :: l_converged(nproma), all_converged
    REAL(wp)             :: wave_num(nproma)
    REAL(wp)             :: AKP
    INTEGER              :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER              :: i_startidx, i_endidx
    TYPE(t_wave_config), POINTER :: wc => NULL()
    !  ---------------------------------------------------------------------------- !
    !     1. START WITH MAXIMUM FROM DEEP AND EXTREME SHALLOW WATER WAVE NUMBER.    !

    i_rlstart  = 1
    i_rlend    = min_rledge
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    wc  => wave_config

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jb,jf,i_startidx,i_endidx,OM,wave_num,l_converged,all_converged,BO,AKP,TH,STH) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        OM = pi2 * wc%freqs(jf)
        AKP = 10000.0_wp
        DO je = i_startidx, i_endidx
          !
          ! initialization
          wave_num(je) = MAX( OM**2/(4.0_wp* grav), OM/(2.0_wp*SQRT(grav*depth(je,jb))) )

          l_converged(je) = (ABS(AKP-wave_num(je)) <= EBS*wave_num(je))
        ENDDO  ! je
        all_converged =  ALL(l_converged(i_startidx:i_endidx))

        ! ---------------------------------------------------------------------------- !
        !     2. ITERATION LOOP.                                                       !

        DO WHILE (.NOT. all_converged)
          DO je = i_startidx, i_endidx
            IF (.NOT.l_converged(je)) THEN
              BO = depth(je,jb)*wave_num(je)
              IF (BO > DKMAX) THEN
                wave_num(je) = OM**2/grav
                l_converged(je) = .TRUE.
              ELSE
                AKP = wave_num(je)
                TH  = grav*wave_num(je)*TANH(BO)
                STH = SQRT(TH)
                wave_num(je) = wave_num(je) &
                  &                + (OM-STH)*STH*2.0_wp / (TH/wave_num(je) + grav*BO/COSH(BO)**2)
                ! check for converged solution
                l_converged(je) = (ABS(AKP-wave_num(je)) <= EBS*wave_num(je))
              END IF
            END IF  ! l_converged
          ENDDO  !je
          all_converged = ALL(l_converged(i_startidx:i_endidx))
        ENDDO !while
        !
        DO je = i_startidx, i_endidx
          wave_num_e(je,jf,jb) = wave_num(je)
        ENDDO
        !
      ENDDO ! jf
    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE wave_number_e

  ! Adaptation of WAM 4.5 code, SUBROUTINE SDEPTHLIM
  !*    PURPOSE.
  !     --------
  ! Limits the level of total energy such that the wave height
  ! does not exceed the maximum wave height allowed for a given depth
  !
  SUBROUTINE sdepth_lim(p_patch, wave_config, depth, emean, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'sdepth_lim'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: depth(:,:)
    REAL(wp),                    INTENT(INOUT) :: emean(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)

    INTEGER  :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb,jc,jf,jd
    REAL(wp) :: em(nproma)

    REAL, PARAMETER :: gamd  = 0.8_wp  !! Parameter of depth limited wave height

    TYPE(t_wave_config), POINTER :: wc => NULL()

    wc  => wave_config

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,i_startidx,i_endidx,em) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        em(jc) = MIN((0.25_wp*gamd*depth(jc,jb))**2/emean(jc,jb),1._wp)
        emean(jc,jb) = em(jc) * emean(jc,jb)
      END DO

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            tracer(jc,jd,jb,jf) = tracer(jc,jd,jb,jf) * em(jc)
          END DO
        END DO
      END DO

    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE sdepth_lim

  !>
  !! Set wave spectrum to absolute allowed minimum
  !! Adaptation of WAM 4.5 code.
  !!
  SUBROUTINE set_energy2emin(p_patch, wave_config, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'set_energy2emin'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            tracer(jc,jd,jb,jf) = MAX(tracer(jc,jd,jb,jf),EMIN)
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE set_energy2emin


  !>
  !! Set wave spectrum to zero according to 0,1 mask by
  !! multiplication of tracers and mask
  !!
  SUBROUTINE mask_energy(p_patch, wave_config, mask, tracer)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = 'mask_energy'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    INTEGER,                     INTENT(IN)    :: mask(:,:)
    REAL(wp),                    INTENT(INOUT) :: tracer(:,:,:,:)

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb,jc,jf,jd

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! save some paperwork
    wc => wave_config
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            tracer(jc,jd,jb,jf) = tracer(jc,jd,jb,jf) * REAL(mask(jc,jb),wp)
          END DO
        END DO
      END DO
    END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE mask_energy


END MODULE mo_wave_physics
