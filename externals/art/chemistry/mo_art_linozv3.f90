!
! mo_art_linozv3
! Adapation of mo_art_linoz to LinozV3
!
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

MODULE mo_art_linozv3
    ! ICON
    USE mo_math_constants,       ONLY: pi, &
                                   &  rad2deg
    USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH
    USE mo_kind,                 ONLY: wp, i8

    USE mtime,                   ONLY: datetime
    
    !USE mo_read_interface,            ONLY: openInputFile, closeFile, read_1D, read_1D_extdim_time
    USE mo_exception,                 ONLY: finish
    USE mtime,                        ONLY: datetime, newDatetime, deallocateDatetime, julianday, newJulianday &
       &                                , deallocateJulianday, getJulianDayFromDatetime &
       &                                , newDateTime, deallocateDateTime, no_of_ms_in_a_day &
       &                                , getNoOfDaysInYearDateTime, getDayOfYearFromDateTime

    ! ART
    USE mo_art_chem_types_param, ONLY: t_chem_meta_linozv3
    USE mo_art_data,             ONLY: p_art_data
    USE mo_art_chem_data,        ONLY: t_art_chem,  &
                                   &   t_art_linozv3
    USE mo_art_atmo_data,        ONLY: t_art_atmo
!    USE mtime_datetime,          ONLY: getJulianDayFromDatetime

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC   ::                        &
      &   art_calc_linozv3,            &
      &   art_calc_linoz_anav3,        &
      &   art_calc_linoz_polarchemv3,  &
      &   art_calc_linoz_polarchem_ltv3
         
CONTAINS

! ----------------------------------

SUBROUTINE read_f107index(valid_datetime, f107_int,rjd_now, tracer)
  TYPE(t_chem_meta_linozv3), INTENT(INOUT) ::  &
    &  tracer       !< structure for linoz tracer   
  TYPE(datetime), POINTER  :: valid_datetime
  REAL(wp)   :: f107_int

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
                routine='read_f107index'

  !---LOCAL VARIABLES
  TYPE(julianday), POINTER :: jd_now
  REAL(wp):: w2, w1
  INTEGER             :: ntime
  INTEGER :: errno, iref1, iref2, index
  !--- the reference julian date in the f107 file
  REAL(wp), PARAMETER :: rjd_ref = 2396759._wp ! 1850-01-01
  REAL(wp)            :: rjd_now, time_diff

  !--- get the julian date of current date
  !    prepare c-pointers according mtime lib
  jd_now => newJulianday(0_i8,0_i8,errno)
  CALL getJulianDayFromDatetime(valid_datetime, jd_now)

  rjd_now = jd_now%day + jd_now%ms/86400._wp/1000._wp
  ! days since refdate
  time_diff = rjd_now - rjd_ref
  ntime = SIZE( tracer%f107(:,1))
  IF ( time_diff < tracer%f107(1,1) .OR. time_diff > tracer%f107(ntime,1) ) THEN
    CALL finish( routine, 'Model time outside f107 time series.' )
  ENDIF
  !--- find neighbours
  index = 1
  DO WHILE ( time_diff > tracer%f107(index,1) )
      index = index + 1
  ENDDO
  iref1 = index - 1
  iref2 = index
  !--- calculate weights
  w1 = ( tracer%f107(iref2,1) - time_diff )/(tracer%f107( iref2,1 ) - tracer%f107( iref1,1 ) )
  w2 = 1._wp - w1
  IF ( w2 < 0._wp .OR. w2 > 1._wp ) THEN
      CALL finish( routine, 'Weights exceed 0,1.' )
  END IF
  f107_int = w1*tracer%f107(iref1,2) + w2*tracer%f107(iref2,2)
END SUBROUTINE read_f107index

SUBROUTINE art_calc_linozv3(jg,jb,jcs,jce, tracer, current_date, p_dtime )

!<
! Adapation of SUBROUTINE art_calc_linoz to LinozV3
!>

  INTEGER, INTENT(IN) :: &
    &  jg,               &  !< patch id
    &  jb,jcs,jce           !< loop indices
  TYPE(t_chem_meta_linozv3), INTENT(INOUT) ::  &
    &  tracer       !< structure for linoz traceri   
  TYPE(datetime), POINTER, INTENT(IN)  ::  &
    &  current_date !< current simulation date
  REAL(wp), INTENT(IN) ::  &
    &  p_dtime      !< model time step
  ! local variables
  INTEGER ::   &
     &  jk, jc &
     &, nlev        !< number of full levels
 
  !------------------------------------
  !---Initialize the linoz routine
  !---Local variables:
  !------------------------------------

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
                routine='art_calc_linozv3'
  INTEGER             :: ik,il,ilp
  REAL(wp)            :: alf,alf1
  REAL(wp)            :: o3ss_bl  ! Boundary layer O3 VMR,in ppb
  REAL(wp)            :: tau_bl   ! Relaxation time in boundary layer, in seconds
  REAL(wp)            :: z_logp(p_art_data(jg)%atmo%nproma, & 
    &                           p_art_data(jg)%atmo%nlev, &
    &                           p_art_data(jg)%atmo%nblks)     ! Pressure altitude

  REAL(wp)            :: rjd_now
  REAL(wp)            :: f107_int
  ! The minmax spectra are derived from ATLAS3 measurements provided by Markus Kunze
  ! sf107 is taken for the respective periods
  REAL(wp), PARAMETER :: sol10cm_max = 230.00_wp ! in Solar Radio Flux Unit sfu = 10^4 Jansky at 10.7 cm
  REAL(wp), PARAMETER :: sol10cm_min = 79.143_wp ! sfu
  REAL(wp):: wf2, wf1


  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART atmo fields
  TYPE(t_art_chem),POINTER    :: &
    &  art_chem                     !< Pointer to ART chem fields
  TYPE(t_art_linozv3),POINTER    :: &
    &  art_linozv3                    !< Pointer to ART linoz fields

  !Set values for tropospheric O3
  o3ss_bl = 25e-9_wp
  tau_bl = 2.*86400._wp ! 2 days

  !------------------------------------
  !---Allocate the linoz routine
  !---Local variables:
  !------------------------------------

  art_atmo => p_art_data(jg)%atmo
  art_chem => p_art_data(jg)%chem
  art_linozv3 => p_art_data(jg)%chem%param%linozv3
   
  nlev = art_atmo%nlev

  !Calculate pressure altitude from ICON pressure
  z_logp = 16._wp*LOG10(100000._wp/art_atmo%pres)
  
  !set wf2 = 0.0 for solarmin and = 1.0 for solarmax
  SELECT CASE (tracer%sol_ssi)
    CASE (0)     
      wf2 = 0.5_wp
    CASE (1)
      wf2 = 0.0_wp
    CASE (2)
      wf2 = 1.0_wp
    CASE (3)
      CAll  read_f107index(current_date, f107_int,rjd_now, tracer )
      !calculate weighting factor for solmin and solmax
      wf2 = (f107_int - sol10cm_min)/(sol10cm_max -sol10cm_min)
  END SELECT
  IF ( wf2 < 0._wp .OR. wf2 > 1._wp ) THEN
      print *, wf2
      CALL finish( routine, 'Weights exceed 0,1.' )
  END IF
  wf1 = 1.0_wp - wf2
 
  DO jc = jcs,jce
    !  ---------------------------------------------
    !  Linear interpolation of Coefficients in Latitude and
    !  averaging in altitude
    !  ---------------------------------------------
    DO jk=1,nlev
      ik = INT( (art_atmo%lat(jc,jb)+pi/2.)/pi*18._wp ) +1
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

      alf = (z_logp(jc,jk,jb)-8._wp - il*2._wp)/2.0_wp

      IF ((alf > 1.0_wp) .OR. (alf <= 0.0_wp)) THEN
        alf = 0.0_wp
        ilp = il
      END IF
  
      alf1 = 1.0_wp - alf

      art_linozv3%linozv3_tab1(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,1) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,1) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,1) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,1)

      art_linozv3%linozv3_tab2(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,2) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,2) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,2) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,2)

      art_linozv3%linozv3_tab3(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,3) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,3) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,3) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,3)

      art_linozv3%linozv3_tab4(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,4) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,4) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,4) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,4)

      art_linozv3%linozv3_tab5(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,5) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,5) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,5) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,5)

      art_linozv3%linozv3_tab6(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,6) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,6) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,6) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,6)

      art_linozv3%linozv3_tab7(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,7) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,7) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,7) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,7)

      art_linozv3%linozv3_tab8(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,8) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,8) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,8) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,8)

      art_linozv3%linozv3_tab9(jc,jk) = wf1*alf1*art_linozv3%tparm_min(il,ik,current_date%date%month,9) &
           &                           +wf1*alf*art_linozv3%tparm_min(ilp,ik,current_date%date%month,9) &
           &                           +wf2*alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,9) &
           &                           +wf2*alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,9)

    ENDDO
    DO jk = 1,nlev-3
      IF (z_logp(jc,jk,jb) >= 10._wp) THEN
        tracer%tend(jc,jk,jb)     =    ( art_linozv3%linozv3_tab5(jc,jk)                          &
                         & + art_linozv3%linozv3_tab6(jc,jk) * (tracer%tracer(jc,jk,jb)           &
                         &    / art_chem%vmr2Nconc(jc,jk,jb) - art_linozv3%linozv3_tab1(jc,jk))   &
                         & + art_linozv3%linozv3_tab7(jc,jk) * (tracer%NOy_tracer(jc,jk,jb)       &
                         &    / art_chem%vmr2Nconc(jc,jk,jb) - art_linozv3%linozv3_tab2(jc,jk))   &
                         & + art_linozv3%linozv3_tab8(jc,jk)                                      &
                         &    * (art_atmo%temp(jc,jk,jb) - art_linozv3%linozv3_tab3(jc,jk))       &    
                         & + art_linozv3%linozv3_tab9(jc,jk)                                      &
                         &    * (tracer%column(jc,jk,jb) - art_linozv3%linozv3_tab4(jc,jk))       &
                         &              )                                                         &
                         & * art_chem%vmr2Nconc(jc,jk,jb)

        tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb)*p_dtime
        IF ( tracer%tracer(jc,jk,jb) + tracer%tend(jc,jk,jb) < 0._wp ) tracer%tend(jc,jk,jb) = -tracer%tracer(jc,jk,jb)
      ENDIF
    ENDDO

    DO jk = nlev-2,nlev
      tracer%tend(jc,jk,jb) =  (o3ss_bl*art_chem%vmr2Nconc(jc,jk,jb) - tracer%tracer(jc,jk,jb)) &
                       & * (1._wp - EXP(-1_wp*p_dtime/tau_bl))
    ENDDO

  ENDDO

END SUBROUTINE art_calc_linozv3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_calc_linoz_anav3(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
  !<
  ! Adaption of SUBROUTINE art_calc_linoz_ana to LinozV3
  !>

  IMPLICIT NONE

  INTEGER, INTENT(in) :: &
    &  jg,               &  !< patch id
    &  jb,jcs,jce           !< loop indices
  TYPE(t_chem_meta_linozv3), INTENT(inout) ::  &
    &  tracer       !< structure for linoz tracer
  TYPE(datetime), POINTER, INTENT(IN)  ::  &
    &  current_date !< current simulation date
  REAL(wp), INTENT(IN) ::  &
    &  p_dtime      !< model time step

  INTEGER :: jk, jc,             &
     &    nlev                                   !< number of full levels


  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo                     !< Pointer to ART atmo fields 
  TYPE(t_art_chem),POINTER    :: &
     &  art_chem                     !< Pointer to ART chem fields 
  TYPE(t_art_linozv3),POINTER    :: &
     &  art_linozv3                     !< Pointer to ART linoz fields 

  REAL(wp)    :: o3ss, tau

  INTEGER             :: ik,il,ilp
  REAL(wp)                :: alf,alf1
  
  !Variable for pressure altitude
  REAL(wp), DIMENSION(:, :, :), ALLOCATABLE :: z_logp

  !------------------------------------
  !---Allocate the linoz routine
  !---Local variables:
  !------------------------------------

  art_atmo => p_art_data(jg)%atmo
  art_chem => p_art_data(jg)%chem
  art_linozv3 => p_art_data(jg)%chem%param%linozv3

  ALLOCATE(z_logp(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))

  nlev = art_atmo%nlev
  z_logp = 16._wp*LOG10(100000._wp/art_atmo%pres)


  DO jc = jcs,jce
    !  ---------------------------------------------
    !  Linear interpolation of Coefficients in Latitude and
    !  averaging in altitude
    !  ---------------------------------------------
    DO jk=1,nlev
      ik = INT( (art_atmo%lat(jc,jb)+pi/2._wp)/pi*18._wp ) +1
      il = INT((z_logp(jc,jk,jb) - 8._wp)/2._wp)
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

      art_linozv3%linozv3_tab1(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,1) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,1)

      art_linozv3%linozv3_tab2(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,2) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,2)

      art_linozv3%linozv3_tab3(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,3) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,3)

      art_linozv3%linozv3_tab4(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,4) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,4)

      art_linozv3%linozv3_tab5(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,5) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,5)

      art_linozv3%linozv3_tab6(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,6) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,6)

      art_linozv3%linozv3_tab7(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,7) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,7)

      art_linozv3%linozv3_tab8(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,8) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,8)

      art_linozv3%linozv3_tab9(jc,jk) = alf1*art_linozv3%tparm_max(il,ik,current_date%date%month,9) &
           &                        +alf*art_linozv3%tparm_max(ilp,ik,current_date%date%month,9)
    ENDDO


    DO jk = 1,nlev
      IF (z_logp(jc,jk,jb) >= 10._wp) THEN

        tau = -1._wp / (art_linozv3%linozv3_tab6(jc,jk))

        o3ss =  (art_linozv3%linozv3_tab1(jc,jk) + (art_linozv3%linozv3_tab5(jc,jk)    &
                & + art_linozv3%linozv3_tab8(jc,jk) * (art_atmo%temp(jc,jk,jb)     &
                & - art_linozv3%linozv3_tab3(jc,jk)) + art_linozv3%linozv3_tab9(jc,jk) &
                & * (tracer%column(jc,jk,jb) - art_linozv3%linozv3_tab4(jc,jk)) )  &
                & * tau) * art_chem%vmr2Nconc(jc,jk,jb)

        tracer%tend(jc,jk,jb) = (o3ss - tracer%tracer(jc,jk,jb))               &
                & *(1._wp - EXP(-1._wp * p_dtime / tau))

      ELSE
        tracer%tend(jc,jk,jb) = tracer%tracer(jc,jk,jb) * EXP(-1._wp*p_dtime*tracer%des)  &
               &                  - tracer%tracer(jc,jk,jb)
      ENDIF

    ENDDO

  END DO

  DEALLOCATE(z_logp)
END SUBROUTINE art_calc_linoz_anav3

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_calc_linoz_polarchem_ltv3(jg,jb,jcs,jce ,&
            &             tracer, p_dtime )
!< 
!---------------------------------------------------------
!--- Calculate non-linear heterogenous ozone loss for polar areas
!--- Adaption of SUBROUTINE art_calc_linoz_polarchem_lt to LinozV3
!---------------------------------------------------------
!>
  INTEGER, INTENT(IN) ::   &
    &  jg,jb,jcs,jce   !< patch on which computation is performed, loop indices
!  TYPE(datetime), POINTER, INTENT(in)    ::   &
!    &  current_date    !< current date and time
  TYPE(t_chem_meta_linozv3), INTENT(inout) ::   &
    &  tracer          !< tracer structure
  REAL(wp),                INTENT(in)    ::   &
    &  p_dtime         !< time step

  ! local variables
  REAL(wp)   :: XSZA = 90._wp  !< Threshold for solar zenith angle
  REAL(wp)   :: sza            !< solar zenith angle in degrees
  REAL(wp)   :: zlat           !< latitude (deg)


  INTEGER    :: jk, jc
  TYPE(t_art_atmo),POINTER    :: &
  &  art_atmo                     !< Pointer to ART diagnostic fields

  art_atmo => p_art_data(jg)%atmo

  !--------------------------------------------------------------------

  DO jc=jcs,jce
    zlat = art_atmo%lat(jc,jb)*rad2deg
    IF (ABS(zlat) > 45._wp) THEN
      DO jk = 1,art_atmo%nlev
        sza = art_atmo%sza_deg(jc,jb)      
        IF ( art_atmo%temp(jc,jk,jb) <= tracer%Thet .AND. sza <= XSZA ) THEN
          tracer%tend(jc,jk,jb) = tracer%tracer(jc,jk,jb)  &
                &  * EXP(-1._wp*p_dtime / tracer%o3lt_het) - tracer%tracer(jc,jk,jb) 
        END IF
      ENDDO
    END IF
  END DO
  
END SUBROUTINE art_calc_linoz_polarchem_ltv3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_calc_linoz_polarchemv3(jg,jb,jcs,jce ,&
            &             tracer, p_dtime )
!< 
!---------------------------------------------------------
! Adaptation of SUBROUTINE art_calc_linoz_polarchem to LinozV3
!---------------------------------------------------------
!>
  INTEGER, INTENT(IN) ::   &
    &  jg,jb,jcs,jce   !< patch on which computation is performed, loop indices
!  TYPE(datetime), POINTER, INTENT(in)    ::   &
!    &  current_date    !< current date and time
  TYPE(t_chem_meta_linozv3), INTENT(inout) ::   &
    &  tracer          !< tracer structure
  REAL(wp),                INTENT(in)    ::   &
    &  p_dtime         !< time step

  ! local variables
  REAL(wp)   :: XSZA = 90._wp !< Threshold for solar zenith angle
  REAL(wp)   :: sza           !< solar zenith angle in degrees
  REAL(wp)   :: rate_lt
  REAL(wp)   :: zlat          !< latitude (deg)

  INTEGER    :: jk, jc
  TYPE(t_art_atmo),POINTER    :: &
  &  art_atmo                     !< Pointer to ART diagnostic fields
  TYPE(t_art_chem),POINTER    :: &
  &  art_chem                     !< Pointer to ART diagnostic fields

  art_chem => p_art_data(jg)%chem
  art_atmo => p_art_data(jg)%atmo

  !--------------------------------------------------------------------

  DO jc = jcs,jce
    zlat = art_atmo%lat(jc,jb)*rad2deg
    IF (ABS(zlat) > 45._wp) THEN
      DO jk = 1,art_atmo%nlev
        sza = art_atmo%sza_deg(jc,jb)
      
        IF ( sza <= XSZA ) THEN 
          rate_lt            = 1._wp / tracer%o3lt_het  * tracer%cold_tracer(jc,jk,jb) &
                            &         / art_chem%vmr2Nconc(jc,jk,jb)

          tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) + tracer%tracer(jc,jk,jb)   &
                 &              * exp(-1._wp*p_dtime * rate_lt) - tracer%tracer(jc,jk,jb)
        END IF   
      ENDDO
    END IF
  END DO
  
END SUBROUTINE art_calc_linoz_polarchemv3

END MODULE mo_art_linozv3
