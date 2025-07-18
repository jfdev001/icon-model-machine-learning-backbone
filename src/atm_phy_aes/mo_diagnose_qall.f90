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

! Module containing subroutine summing up the mass fractions in air of all
! hydrometeors: cloud liquid water, cloud ice, rain, snow, and graupel

MODULE mo_diagnose_qall

  USE mo_kind            ,ONLY: wp

  USE mo_run_config      ,ONLY: num_lev, iqv, iqc, iqi, iqr, iqs, iqg
  USE mo_dynamics_config ,ONLY: nnew, nnew_rcf

  USE mo_nonhydro_state  ,ONLY: p_nh_state
  USE mo_aes_phy_memory  ,ONLY: prm_field

  USE mo_timer           ,ONLY: ltimer, timer_start, timer_stop, timer_qall

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: diagnose_qall

CONTAINS

  !-------------------------------------------------------------------

  SUBROUTINE diagnose_qall(jg, jb, jcs, jce)

    INTEGER , INTENT(in)    :: jg, jb, jcs, jce
    INTEGER                 :: jc, jk, jtl_trc

    IF (ASSOCIATED(prm_field(jg)%qall)) THEN
      !
      IF (ltimer) CALL timer_start(timer_qall)
      !
      jtl_trc = nnew_rcf(jg)
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, num_lev(jg)
        DO jc = jcs, jce
          prm_field(jg)%qall(jc,jk,jb) = p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqc) &
                                     & + p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqi) &
                                     & + p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqr) &
                                     & + p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqs) &
                                     & + p_nh_state(jg)%prog(jtl_trc)%tracer(jc,jk,jb,iqg)
        END DO ! jc
      END DO ! jk
      !$ACC END PARALLEL
      !
      IF (ltimer) CALL timer_stop(timer_qall)
      !
    END IF

  END SUBROUTINE  diagnose_qall

  !-------------------------------------------------------------------

END MODULE mo_diagnose_qall
