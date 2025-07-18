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

! Module containing subroutine summing up the vertical integrals of cloud liquid water
! and cloud ice.

MODULE mo_diagnose_clwvi

  USE mo_kind            ,ONLY: wp

  USE mo_run_config      ,ONLY: iqc, iqi

  USE mo_aes_phy_memory  ,ONLY: prm_field

  USE mo_timer           ,ONLY: ltimer, timer_start, timer_stop, timer_clwvi

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: diagnose_clwvi

CONTAINS

  !-------------------------------------------------------------------

  SUBROUTINE diagnose_clwvi(jg, jb, jcs, jce)

    INTEGER , INTENT(in)    :: jg, jb, jcs, jce
    INTEGER                 :: jc

    IF (ASSOCIATED(prm_field(jg)%clwvi)) THEN
      !
      IF (ltimer) CALL timer_start(timer_clwvi)
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = jcs, jce
        prm_field(jg)%clwvi(jc,jb) = prm_field(jg)%mtrcvi(jc,jb,iqc) &
                                 & + prm_field(jg)%mtrcvi(jc,jb,iqi)
      END DO ! jc
      !$ACC END PARALLEL
      !
      IF (ltimer) CALL timer_stop(timer_clwvi)
      !
    END IF

  END SUBROUTINE  diagnose_clwvi

  !-------------------------------------------------------------------

END MODULE mo_diagnose_clwvi
