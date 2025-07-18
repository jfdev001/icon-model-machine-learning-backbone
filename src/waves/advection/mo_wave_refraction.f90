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

! Refraction of spectral surface wave energy
!
! Literature: The WAM model - A third Generation Ocean Wave Prediction Model, 1988

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_wave_refraction

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_model_domain,        ONLY: t_patch
  USE mo_wave_config,         ONLY: t_wave_config
  USE mo_grid_config,         ONLY: grid_sphere_radius
  USE mo_parallel_config,     ONLY: nproma
  USE mo_math_constants,      ONLY: pi, pi2
  USE mo_loopindices,         ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_refraction'

  PUBLIC :: wave_refraction

CONTAINS

  !>
  !! Calculate grid, depth and current refraction
  !!
  !! Based on WAM shallow water with depth and current refraction
  !! P_SPHER_SHALLOW_CURR
  !!
  SUBROUTINE wave_refraction(p_patch, wave_config, dtime, wave_num_c, gv_c, depth, &
    &                        depth_grad, tracer_now, tracer_new)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
         routine = modname//':wave_refraction'

    TYPE(t_patch),               INTENT(IN)   :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)   :: wave_config
    REAL(wp),                    INTENT(IN)   :: dtime                ! integration time step [s]
    REAL(wp),                    INTENT(IN)   :: wave_num_c(:,:,:)
    REAL(wp),                    INTENT(IN)   :: gv_c(:,:,:)          ! group velocity at cell centers
    REAL(wp),                    INTENT(IN)   :: depth(:,:)
    REAL(wp),                    INTENT(IN)   :: depth_grad(:,:,:)    ! bathymetry gradient (2,nproma,nblks_c)
    REAL(wp), TARGET,            INTENT(IN)   :: tracer_now(:,:,:,:)  ! energy before transport
    REAL(wp),                    INTENT(INOUT):: tracer_new(:,:,:,:)  ! (nproma,ndirs,nblks_c,nfreqs)


    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jf,jd,isub
    INTEGER :: jdm1,jdp1                       ! index of direction -1/+1

    REAL(wp) :: DELTHR, DELTH, DELTR, DELTH0, sm, sp, akd, DTP, DTM, dDTC, temp, dtime_sub
    REAL(wp) :: thdd(nproma,wave_config%ndirs)
    REAL(wp) :: delta_ref(nproma,wave_config%ndirs,wave_config%nfreqs)
    REAL(wp) :: tsihkd(nproma), tan_lat(nproma)
    REAL(wp), TARGET :: tracer_tmp(nproma,wave_config%ndirs,wave_config%nfreqs)
    REAL(wp), POINTER :: tracer_ptr(:,:,:)

    wc => wave_config

    dtime_sub = dtime/REAL(wc%nsubs_refrac,wp) ! substepping time step

    DELTH = 2.0_wp*pi/REAL(wc%ndirs,wp)
    DELTR = DELTH * grid_sphere_radius
    DELTH0 = 0.5_wp * dtime_sub / DELTR
    DELTHR = 0.5_wp * dtime_sub / DELTH

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,i_startidx,i_endidx,isub,temp,akd,tsihkd,thdd,tracer_ptr, &
!$OMP            tracer_tmp,sm,sp,jdm1,jdp1,dtp,dtm,dDTC,delta_ref,tan_lat) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
           &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        tan_lat(jc) = TAN(p_patch%cells%center(jc,jb)%lat)
      ENDDO

      DO isub = 1,wc%nsubs_refrac
        IF (isub == 1) THEN
          tracer_ptr => tracer_now(:,:,jb,:)
        ELSE
          tracer_ptr => tracer_tmp
        ENDIF

        DO jf = 1,wc%nfreqs

          DO jc = i_startidx, i_endidx
            akd = wave_num_c(jc,jf,jb) * depth(jc,jb)
            IF (akd <= 10.0_wp) THEN
              tsihkd(jc) = (pi2 * wc%freqs(jf))/SINH(2.0_wp*akd)
            ELSE
              tsihkd(jc) = 0.0_wp
            END IF
          ENDDO

          DO jd = 1,wc%ndirs
            DO jc = i_startidx, i_endidx

              temp = (wc%sin_dir(jd) + wc%sin_dir(wc%dir_neig_ind(2,jd))) * depth_grad(2,jc,jb) &
                   - (wc%cos_dir(jd) + wc%cos_dir(wc%dir_neig_ind(2,jd))) * depth_grad(1,jc,jb)

              thdd(jc,jd) = temp * tsihkd(jc)

            END DO !jc
          END DO !jd


          DO jd = 1,wc%ndirs

            jdm1 = wc%dir_neig_ind(1,jd)  !index of direction - 1
            jdp1 = wc%dir_neig_ind(2,jd)  !index of direction + 1

            sm = DELTH0 * (wc%sin_dir(jd) + wc%sin_dir(jdm1))
            sp = DELTH0 * (wc%sin_dir(jd) + wc%sin_dir(jdp1))

            DO jc = i_startidx, i_endidx

              DTP = tan_lat(jc) * gv_c(jc,jf,jb)

              DTM = DTP * SM + thdd(jc,jdm1) * DELTHR
              DTP = DTP * SP + thdd(jc,jd)   * DELTHR

              dDTC = -MAX(0._wp , DTP) + MIN(0._wp , DTM)
              DTP  = -MIN(0._wp , DTP)
              DTM  =  MAX(0._wp , DTM)

              delta_ref(jc,jd,jf) = dDTC*tracer_ptr(jc,jd,jf) + DTM*tracer_ptr(jc,jdm1,jf) + DTP*tracer_ptr(jc,jdp1,jf)
            END DO !jc
          END DO !jd
        END DO !jf

        IF (isub < wc%nsubs_refrac) THEN
          DO jf = 1,wc%nfreqs
            DO jd = 1,wc%ndirs
              DO jc = i_startidx, i_endidx
                tracer_tmp(jc,jd,jf) = tracer_ptr(jc,jd,jf) + delta_ref(jc,jd,jf)
              END DO !jc
            END DO !jd
          END DO !jf
        ELSE
          DO jf = 1,wc%nfreqs
            DO jd = 1,wc%ndirs
              DO jc = i_startidx, i_endidx
                tracer_new(jc,jd,jb,jf) = tracer_new(jc,jd,jb,jf) + delta_ref(jc,jd,jf) + &
                                          (tracer_ptr(jc,jd,jf)-tracer_now(jc,jd,jb,jf))
              END DO !jc
            END DO !jd
          END DO !jf
        ENDIF

      ENDDO ! isub

    END DO !jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE wave_refraction

END MODULE mo_wave_refraction
