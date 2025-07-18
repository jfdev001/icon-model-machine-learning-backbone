!
! mo_art_linoz
! This module provides a simple linear ozone chemistry
! first introduced by McLinden 2000
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.md for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_art_linoz
    ! ICON
    USE mo_math_constants,       ONLY: pi, pi_2, rad2deg
    USE mo_kind,                 ONLY: wp

    USE mtime,                   ONLY: datetime

    ! ART
    USE mo_art_chem_types_param, ONLY: t_chem_meta_linoz
    USE mo_art_data,             ONLY: p_art_data
    USE mo_art_chem_data,        ONLY: t_art_chem,  &
                                   &   t_art_linoz
    USE mo_art_atmo_data,        ONLY: t_art_atmo



  IMPLICIT NONE


  PRIVATE

  PUBLIC   ::                        &
      &   art_calc_linoz,            &
      &   art_calc_linoz_ana,        &
      &   art_calc_linoz_polarchem,  &
      &   art_calc_linoz_polarchem_lt

CONTAINS



! ----------------------------------
! --- Linear Ozone (mclinden 2002)
! --- by CS/IMK-ASF 10.09.2014
! --- Polar Ozoneloss with constant lifetime
! --- Lifetime assumption based on Sinnhuber, 2003
! --- Changes done by Katerina Kusakova, 09/2024:
! --- + Vertical coordinate changed from z_mc to z_logp 
!       for consistency with look-up tables
! --- + Tropospheric ozone production and surface deposition is calculated via
!       relaxing within the 3 lowermost levels towards fixed volume mixing 
!       ratio (o3_lbc) with lifetime (lt_lbc). This values can be changed via 
!       xml file. Per default (following McLinden 2002): o3_lbc=25ppb, 
!       lt_lbc=2 days.
! --- + Latitude threshold for polar ozoneloss was 
!       added to both polar routines to avoid activation in the tropics
! --- + xsza threshold for polar ozoneloss in art_calc_linoz_polarchem_lt
!        was changed from 85 degrees to 90 degrees according to M.Braun, 2021
! ----------------------------------
SUBROUTINE art_calc_linoz(jg,jb,jcs,jce, tracer, current_date, p_dtime )

!<
! SUBROUTINE art_calc_linoz
! Routine to calculate the ozone concentration
! based on McLinden 2000
! Author: Christian Stassen, KIT
! Initial Release: 2015-03-24
! Changed by Katerina Kusakova, 2024-09-12 
! + added relaxation of tropospheric ozone for 3 lowermost levels
! + pressure altitude z_log is used as vertical coordinate                    
!>

  INTEGER, INTENT(in) :: &
    &  jg,               &               !< patch id
    &  jb,jcs,jce                        !< loop indices
  TYPE(t_chem_meta_linoz), INTENT(inout) ::  &
    &  tracer                            !< structure for linoz tracer
  TYPE(datetime), POINTER, INTENT(IN)  ::  &
    &  current_date                      !< current simulation date
  REAL(wp), INTENT(IN) ::  &
    &  p_dtime                           !< model time step
  ! local variables
  INTEGER :: jk, jc,       &
     &    nlev                           !< number of full levels



  !------------------------------------
  !---Initialize the linoz routine
  !---Local variables:
  !------------------------------------

  INTEGER             :: ik,il,ilp
  REAL(wp)            :: alf,alf1
  REAL(wp)            :: z_logp(p_art_data(jg)%atmo%nproma, & 
    &                           p_art_data(jg)%atmo%nlev, &
    &                           p_art_data(jg)%atmo%nblks)     ! Pressure altitude
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART atmo fields
  TYPE(t_art_chem),POINTER    :: &
    &  art_chem                     !< Pointer to ART chem fields
  TYPE(t_art_linoz),POINTER    :: &
    &  art_linoz                    !< Pointer to ART linoz fields
 

  !------------------------------------
  !---Allocate the linoz routine
  !---Local variables:
  !------------------------------------

  art_atmo => p_art_data(jg)%atmo
  art_chem => p_art_data(jg)%chem
  art_linoz => p_art_data(jg)%chem%param%linoz

  nlev = art_atmo%nlev

  !Calculate pressure altitude from ICON pressure
  z_logp = 16._wp*LOG10(100000._wp/art_atmo%pres)

  DO jc = jcs,jce
    !  ---------------------------------------------
    !  Linear interpolation of Coefficients in Latitude and
    !  averaging in altitude
    !  ---------------------------------------------
    DO jk=1,nlev
      ik = INT( (art_atmo%lat(jc,jb)+pi_2)/pi*18._wp ) +1
      il = INT((z_logp(jc,jk,jb)-8._wp)/2._wp)
      ilp = il+1
      IF (il >= 25) THEN
        il = 25
        ilp = 25
      END IF
      IF (il < 1)  THEN
        il = 1
        ilp = 1
      END IF

      alf = (z_logp(jc,jk,jb) - REAL(il*2+8))/2.0_wp

      IF ((alf > 1.0_wp) .OR. (alf <= 0.0_wp)) THEN
        alf = 0.0_wp
        ilp = il
      END IF

      alf1 = 1.0_wp - alf


      art_linoz%linoz_tab1(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,1) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,1)

      art_linoz%linoz_tab2(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,2) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,2)

      art_linoz%linoz_tab3(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,3) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,3)

      art_linoz%linoz_tab4(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,4) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,4)

      art_linoz%linoz_tab5(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,5) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,5)

      art_linoz%linoz_tab6(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,6) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,6)

      art_linoz%linoz_tab7(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,7) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,7)

    ENDDO




    DO jk = 1,nlev-3

      IF (z_logp(jc,jk,jb) >= 10._wp) THEN

        tracer%tend(jc,jk,jb)     = (art_linoz%linoz_tab4(jc,jk)                              &
                         & + art_linoz%linoz_tab5(jc,jk) * (tracer%tracer(jc,jk,jb)           &
                         &    / art_chem%vmr2Nconc(jc,jk,jb) - art_linoz%linoz_tab1(jc,jk))   &
                         & + art_linoz%linoz_tab6(jc,jk)                                      &
                         &    * (art_atmo%temp(jc,jk,jb) - art_linoz%linoz_tab2(jc,jk))       &
                         & + art_linoz%linoz_tab7(jc,jk)                                      &
                         &    * (tracer%column(jc,jk,jb) - art_linoz%linoz_tab3(jc,jk)) )     &
                         & * art_chem%vmr2Nconc(jc,jk,jb)

        tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) *p_dtime
      ENDIF
    ENDDO

    DO jk = nlev-2,nlev
      tracer%tend(jc,jk,jb) =  (tracer%o3_lbc*art_chem%vmr2Nconc(jc,jk,jb) - tracer%tracer(jc,jk,jb)) &
                       & * (1._wp - EXP(-1._wp*p_dtime/tracer%lt_lbc))
    ENDDO

  ENDDO


END SUBROUTINE art_calc_linoz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_calc_linoz_ana(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
! SUBROUTINE art_calc_linoz_ana
! Routine to calculate the ozone concentration
! based on McLinden 2000
! Author: Christian Stassen, KIT
! Initial Release: 2015-03-24
! Changed by Katerina Kusakova, 2024-09-12 
! + added relaxation of tropospheric ozone for 3 lowermost levels
! + pressure altitude z_log is used as vertical coordinate                    
!>

  IMPLICIT NONE

  INTEGER, INTENT(in) :: &
    &  jg,               &                  !< patch id
    &  jb,jcs,jce                           !< loop indices
  TYPE(t_chem_meta_linoz), INTENT(inout) ::  &
    &  tracer                               !< structure for linoz tracer
  TYPE(datetime), POINTER, INTENT(IN)    ::  &
    &  current_date                         !< current simulation date
  REAL(wp), INTENT(IN) ::  &
    &  p_dtime                              !< model time step

  INTEGER :: jk, jc,       &
     &    nlev                              !< number of full levels

  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo                            !< Pointer to ART atmo fields
  TYPE(t_art_chem),POINTER    :: &
     &  art_chem                            !< Pointer to ART chem fields
  TYPE(t_art_linoz),POINTER   :: &
     &  art_linoz                           !< Pointer to ART linoz fields
  REAL(wp)                    :: z_logp(p_art_data(jg)%atmo%nproma, & 
    &                                   p_art_data(jg)%atmo%nlev, &
    &                                   p_art_data(jg)%atmo%nblks)     ! Pressure altitude

  REAL(wp) :: o3ss           !< steady state mixing ratio
  REAL(wp) :: tau            !< time constant, steady state

  INTEGER  :: ik,il,ilp
  REAL     :: alf,alf1

  !------------------------------------
  !---Allocate the linoz routine
  !---Local variables:
  !------------------------------------

  art_atmo => p_art_data(jg)%atmo
  art_chem => p_art_data(jg)%chem
  art_linoz => p_art_data(jg)%chem%param%linoz

  nlev = art_atmo%nlev

  ! Calculate the pressure altitude
  z_logp = 16._wp*LOG10(100000._wp/art_atmo%pres)

  DO jc = jcs,jce
    !  ---------------------------------------------
    !  Linear interpolation of Coefficients in Latitude and
    !  averaging in altitude
    !  ---------------------------------------------
    DO jk=1,nlev
      ik = INT( (art_atmo%lat(jc,jb)+pi_2)/pi*18._wp ) +1
      il = INT((z_logp(jc,jk,jb)-8._wp)/2._wp)
      ilp = il+1
      IF (il >= 25) THEN
        il = 25
        ilp = 25
      ENDIF

      IF (il < 1) THEN
        il = 1
        ilp = 1
      ENDIF

      alf = (z_logp(jc,jk,jb) - REAL(il*2+8))/2.0_wp

      IF ((alf > 1.0_wp) .OR. (alf <= 0.0_wp)) THEN
        alf = 0.0_wp
        ilp = il
      ENDIF

      alf1 = 1.0_wp - alf

      art_linoz%linoz_tab1(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,1) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,1)

      art_linoz%linoz_tab2(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,2) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,2)

      art_linoz%linoz_tab3(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,3) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,3)

      art_linoz%linoz_tab4(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,4) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,4)

      art_linoz%linoz_tab5(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,5) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,5)

      art_linoz%linoz_tab6(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,6) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,6)

      art_linoz%linoz_tab7(jc,jk) = alf1*art_linoz%tparm(il,ik,current_date%date%month,7) &
           &                        +alf*art_linoz%tparm(ilp,ik,current_date%date%month,7)

    ENDDO




    DO jk = 1,nlev-3
      IF (z_logp(jc,jk,jb) >= 10._wp) THEN

        tau = -1._wp / (art_linoz%linoz_tab5(jc,jk))

        o3ss = (art_linoz%linoz_tab1(jc,jk) + ( art_linoz%linoz_tab4(jc,jk)             &
                     & + art_linoz%linoz_tab6(jc,jk)                                    &
                     &    * (art_atmo%temp(jc,jk,jb) - art_linoz%linoz_tab2(jc,jk))     &
                     & + art_linoz%linoz_tab7(jc,jk)                                    &
                     &    * (tracer%column(jc,jk,jb) - art_linoz%linoz_tab3(jc,jk)) )   &
                     & * tau) * art_chem%vmr2Nconc(jc,jk,jb)

        tracer%tend(jc,jk,jb) = (o3ss - tracer%tracer(jc,jk,jb))   &
              &                  * (1._wp - EXP(-1._wp * p_dtime / tau))
      ENDIF
    ENDDO

    DO jk = nlev-2,nlev
      tracer%tend(jc,jk,jb) =  (tracer%o3_lbc*art_chem%vmr2Nconc(jc,jk,jb) - tracer%tracer(jc,jk,jb)) &
                       & * (1._wp - EXP(-1_wp*p_dtime/tracer%lt_lbc))
    ENDDO
  END DO


END SUBROUTINE art_calc_linoz_ana

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_calc_linoz_polarchem_lt(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
!---------------------------------------------------------
!--- Calculate non-linear heterogenous ozone loss for polar areas
!--- Lifetime assumption based on Sinnhuber, 2003
!--- Author: Christopher Diekmann, KIT
!--- xsza changed to 90 degrees according to M.Braun (KIT, 2021)
!--- Katerina Kusakova, KIT: Latitude threshold zlat= 45 degree added
!---------------------------------------------------------
!>
  INTEGER, INTENT(IN) ::   &
    &  jg,jb,jcs,jce   !< patch on which computation is performed, loop indices
  TYPE(datetime), POINTER, INTENT(in)    ::   &
    &  current_date    !< current date and time
  TYPE(t_chem_meta_linoz), INTENT(inout) ::   &
    &  tracer          !< tracer structure
  REAL(wp),                INTENT(in)    ::   &
    &  p_dtime         !< time step

  ! local variables
  REAL(wp),PARAMETER   :: xsza = 90._wp !< Threshold for solar zenith angle
  REAL(wp)             :: sza           !< solar zenith angle in degrees
  REAL(wp)             :: zlat          !< latitude in degrees
  INTEGER              :: jk, jc
  TYPE(t_art_atmo),POINTER    :: &
  &  art_atmo                     !< Pointer to ART diagnostic fields

  art_atmo => p_art_data(jg)%atmo

  !--------------------------------------------------------------------

  DO jc=jcs,jce
    zlat = art_atmo%lat(jc,jb)*rad2deg
    IF (ABS(zlat) > 45._wp) THEN
      DO jk = 1,art_atmo%nlev
        sza = art_atmo%sza_deg(jc,jb)
      
        IF ( art_atmo%temp(jc,jk,jb) <= tracer%Thet .AND. sza <= xsza ) THEN
          tracer%tend(jc,jk,jb) = tracer%tracer(jc,jk,jb)  &
                &  * EXP(-1._wp*p_dtime / tracer%o3lt_het) - tracer%tracer(jc,jk,jb) 
        END IF
      ENDDO
    END IF
  END DO

END SUBROUTINE art_calc_linoz_polarchem_lt
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_calc_linoz_polarchem(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
!---------------------------------------------------------
!--- Calculate non-linear heterogenous ozone loss for polar areas
!--- Lifetime assumption based on Sinnhuber, 2003
!--- Author: Christopher Diekmann, KIT
!--- Katerina Kusakova, KIT: Latitude threshold zlat= 45 degree added
!---------------------------------------------------------
!>
  INTEGER, INTENT(IN) ::   &
    &  jg,jb,jcs,jce   !< patch on which computation is performed, loop indices
  TYPE(datetime), POINTER, INTENT(in)    ::   &
    &  current_date    !< current date and time
  TYPE(t_chem_meta_linoz), INTENT(inout) ::   &
    &  tracer          !< tracer structure
  REAL(wp),                INTENT(in)    ::   &
    &  p_dtime         !< time step

  ! local variables
  REAL(wp), PARAMETER   :: xsza = 85     !< Threshold for solar zenith angle
  REAL(wp)              :: sza           !< solar zenith angle in degrees
  REAL(wp)              :: rate_lt
  REAL(wp)              :: zlat          !< latitude in degrees
  INTEGER               :: jk, jc
  TYPE(t_art_atmo),POINTER    :: &
  &  art_atmo                     !< Pointer to ART diagnostic fields
  TYPE(t_art_chem),POINTER    :: &
  &  art_chem                     !< Pointer to ART chemistry fields

  art_chem => p_art_data(jg)%chem
  art_atmo => p_art_data(jg)%atmo

  !--------------------------------------------------------------------

  DO jc = jcs,jce
    zlat = art_atmo%lat(jc,jb)*rad2deg
    IF (ABS(zlat) > 45._wp) THEN
      DO jk = 1,art_atmo%nlev
        sza = art_atmo%sza_deg(jc,jb)

        IF ( sza <= xsza ) THEN
          rate_lt            = 1._wp / tracer%o3lt_het  * tracer%cold_tracer(jc,jk,jb) &
                        &    / art_chem%vmr2Nconc(jc,jk,jb)

          tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb)+ tracer%tracer(jc,jk,jb)   & 
                        &       * EXP(-1._wp*p_dtime * rate_lt) - tracer%tracer(jc,jk,jb)
        END IF
      ENDDO
    END IF
    DO jk=1,art_atmo%nlev
      IF ((tracer%tracer(jc,jk,jb) +tracer%tend(jc,jk,jb))<= 1.E-25_wp*art_chem%vmr2Nconc(jc,jk,jb)) THEN
        tracer%tracer(jc,jk,jb) =  0.5E-25_wp*art_chem%vmr2Nconc(jc,jk,jb)
        tracer%tend(jc,jk,jb)   =  0.5E-25_wp*art_chem%vmr2Nconc(jc,jk,jb)
      ENDIF
    ENDDO

  END DO

END SUBROUTINE art_calc_linoz_polarchem

END MODULE mo_art_linoz
