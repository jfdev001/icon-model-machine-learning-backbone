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

! Module containing subroutine summing up the vertical integrals of water vapour, liquid water,
! cloud ice, rain, snow, and graupel to the total cloud water

MODULE mo_diagnose_tcw

  USE mo_kind            ,ONLY: wp

  USE mo_run_config      ,ONLY: num_lev, iqv, iqc, iqi, iqr, iqs, iqg
  USE mo_dynamics_config ,ONLY: nnew, nnew_rcf

  USE mo_nonhydro_state  ,ONLY: p_nh_state
  USE mo_aes_phy_memory  ,ONLY: prm_field

  USE mo_timer           ,ONLY: ltimer, timer_start, timer_stop, timer_tcw

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: diagnose_tcw

CONTAINS

  !-------------------------------------------------------------------

  SUBROUTINE diagnose_tcw(jg, jb, jcs, jce)

    INTEGER , INTENT(in)    :: jg, jb, jcs, jce
    INTEGER                 :: jc, jtl_trc

    IF (ltimer) CALL timer_start(timer_tcw)

    jtl_trc = nnew_rcf(jg)

    IF (ASSOCIATED(prm_field(jg)%tcw)) THEN
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        prm_field(jg)%tcw(jc,jb) =   prm_field(jg)%mtrcvi(jc,jb,iqv) &
                                 & + prm_field(jg)%mtrcvi(jc,jb,iqc) &
                                 & + prm_field(jg)%mtrcvi(jc,jb,iqi) &
                                 & + prm_field(jg)%mtrcvi(jc,jb,iqr) &
                                 & + prm_field(jg)%mtrcvi(jc,jb,iqs) &
                                 & + prm_field(jg)%mtrcvi(jc,jb,iqg)
      END DO ! jc
      !$ACC END PARALLEL
      !
    END IF

    IF (ltimer) CALL timer_stop(timer_tcw)

  END SUBROUTINE  diagnose_tcw

  !-------------------------------------------------------------------

END MODULE mo_diagnose_tcw
