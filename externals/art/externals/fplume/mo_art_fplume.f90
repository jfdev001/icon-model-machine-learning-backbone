!
! mo_art_fplume
! This module constucts the atmospheric profiles, searches the actual
! volcanic phase for FPlume, and starts the plume calculations
! (based on FPLUME-1.1 by A.Folch, G.Macedonio, A.Costa, February 2016)
!
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

MODULE mo_art_fplume

  USE mo_art_fplume_bpt,                ONLY: initialize_plume_BPT,initialize_plume_wind,    &
                                          &   solve_plume_BPT
  USE mo_art_external_types,            ONLY: t_art_volc_fplume
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mtime,                            ONLY: datetime,newDatetime,datetimeToString,         &
                                          &   max_datetime_str_len,newEvent,timedelta,       &
                                          &   event, isCurrentEventActive
  USE mo_art_fplume_types,              ONLY: t_fplume_phases
  USE mo_art_config,                    ONLY: art_config

  PRIVATE
 
  PUBLIC:: art_fplume

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_fplume(fplume_init,volc_fplume,z,rho,pres,temp,u,v,sh,jg,jb,jc,          &
               &      profile_nz,current_date,plume_MER_icon,plume_H_icon,              &
               &      plume_MER_H2O,plume_MER_SO2,plume_zv_icon,fplume_on, &
               &      plume_fine_ash_fraction,exit_water_fraction)
 
  IMPLICIT NONE

  TYPE(t_fplume_phases),INTENT(INOUT) ::  &
    &  fplume_init                          !< FPlume data container for one volcano
  TYPE(t_art_volc_fplume),INTENT(INOUT):: &
    &  volc_fplume                          !< Container for ash transport properties
  TYPE(datetime), POINTER, INTENT(IN)  :: &
    &  current_date                           !< mtime object containing current date (ICON)
  REAL(wp), INTENT(IN)                 :: &
    &  z(:,:),                            & !< Geometric height
    &  rho(:,:),                          & !< Density of air
    &  pres(:,:),                         & !< Air pressure
    &  temp(:,:),                         & !< Air temperature
    &  u(:,:),                            & !< Zonal wind
    &  v(:,:),                            & !< Meridional wind
    &  sh(:,:)                          !< specific humidity
  INTEGER, INTENT(IN)                  :: profile_nz,jg,jb,jc
  LOGICAL, INTENT(INOUT)               :: &
    &  fplume_on                              !< was FPlume active?
  ! Variables necessary for transport in ICON
  REAL(wp), INTENT(INOUT)   ::            &
    &  plume_zv_icon,                       & !< vent altitude (m)
    &  plume_MER_icon,                      & !< MER (kg/s) of each eruption phase
    &  plume_H_icon,                        & !< Height (m agl)
    &  plume_fine_ash_fraction,             & !< phase dependent fine ash fraction
    &  plume_MER_H2O,                       & !< MER of H2O (read from .inp file
    &  plume_MER_SO2,                       & !< MER of SO2 (read from .inp file)
    &  exit_water_fraction

  REAL(wp)                          :: &
    &  angle
  INTEGER                   :: &
    &  plume_phase,                       & !< current eruption phase number as defined in .ipn file
    &  i, idt2
  CHARACTER(len=MAX_CHAR_LENGTH),PARAMETER  :: thisroutine='mo_art_fplume:art_fplume'
  INTEGER               :: plume_status               ! status code
  INTEGER,PARAMETER     :: plume_ns = 300             ! number of plume sources (plume+umbrella)
  INTEGER,PARAMETER     :: plume_np = 200             ! number of plume sources (up to NBL)
  REAL(wp)    :: plume_MER
  REAL(wp)    :: fplume_min_height          ! use FPlume above (below Mastin)
  !
  REAL(wp)    :: plume_z(plume_ns)    ! z(ns) (elevation in m above terrain)
  !
  REAL(wp)    :: profile_z(profile_nz)       ! Height (a.s.l.)
  REAL(wp)    :: profile_rho(profile_nz)     ! Air density
  REAL(wp)    :: profile_p(profile_nz)       ! Pressure
  REAL(wp)    :: profile_T(profile_nz)       ! Temperature
  REAL(wp)    :: profile_sh(profile_nz)      ! Specific humidity
  REAL(wp)    :: profile_u(profile_nz)       ! Wind speed
  REAL(wp)    :: profile_ux(profile_nz)      ! x-component of wind vector
  REAL(wp)    :: profile_uy(profile_nz)      ! y-component of wind vector
  REAL(wp)    :: profile_psia(profile_nz)    ! Wind direction (origin at E,anticlockwise, radians)


  plume_status = 0 ! initialize plume status to avoid error when plume not active at beginning
  fplume_min_height=fplume_init%fplume_min_height
  
!----------- PHASES --------------------------------------------------------
  DO idt2=1,fplume_init%nphases
    IF (fplume_init%phase(idt2)%isActive(current_date) .AND. &
    &   fplume_init%phase(idt2)%plume_udt /=0) THEN
      fplume_on = .TRUE.
      plume_phase = idt2
      EXIT
    ELSE
      fplume_on = .FALSE.
      plume_phase = idt2
    ENDIF
  ENDDO
  
  exit_water_fraction = fplume_init%phase(plume_phase)%plume_wvdt                 &
          &           + fplume_init%phase(plume_phase)%plume_wldt                 &
          &           + fplume_init%phase(plume_phase)%plume_wsdt

!----------- Reverse atmospheric profile -----------------------------------
  DO i = 1,profile_nz
    IF (fplume_init%ext_meteo_prof) THEN
      profile_z   (i) = fplume_init%prof_z(i)*1.0e3_wp    ! km in m
      profile_rho (i) = fplume_init%prof_rho(i)           ! in kg/m3
      profile_p   (i) = fplume_init%prof_pres(i)*1.0e2_wp ! hPa in Pa
      profile_T   (i) = fplume_init%prof_T(i)             ! in K
      profile_sh  (i) = fplume_init%prof_sh(i)*1e-3_wp    ! g/kg in kg/kg
      profile_ux  (i) = fplume_init%prof_u(i)             ! m/s
      profile_uy  (i) = fplume_init%prof_v(i)             ! m/s
    ELSE
      profile_z   (profile_nz-i+1) = z(jc,i)           ! in m
      profile_rho (profile_nz-i+1) = rho(jc,i)         ! in kg/m3                                           
      profile_p   (profile_nz-i+1) = pres(jc,i)        ! in Pa                                            
      profile_T   (profile_nz-i+1) = temp(jc,i)        ! in K
      profile_sh  (profile_nz-i+1) = sh(jc,i)          ! kg/kg
      profile_ux  (profile_nz-i+1) = u(jc,i)           ! m/s                                           
      profile_uy  (profile_nz-i+1) = v(jc,i)           ! m/s                                           
    ENDIF
  ENDDO

  DO i = 1,profile_nz
     profile_u(i)  = SQRT(profile_ux(i)*profile_ux(i)+profile_uy(i)*profile_uy(i))
       IF (ABS(profile_ux(i))>1.0e-8_wp) THEN
          angle = ATAN2(profile_uy(i),profile_ux(i))*180.0_wp/pi
       ELSE
          IF (profile_uy(i)>1.0e-8_wp) THEN
             angle = 90.0_wp
          ELSE IF (profile_uy(i)< -1.0e-8_wp) THEN
             angle = 270.0_wp
          ELSE
             angle = 0.0_wp
          ENDIF
       ENDIF
       profile_PSIa(i) = angle*pi/180.0_wp     ! angle in Rad
  ENDDO
  !
  IF (fplume_init%plume_zv<profile_z(1)) &
     & CALL finish(thisroutine,'Vent altitude below surface')
  !
!-------------- Calculate plume properties ----------------------------------------------------------
  IF (fplume_on) THEN
    IF (fplume_init%phase(plume_phase)%plume_Hdt>=fplume_min_height .OR.     & 
      & fplume_init%plume_solve_for =='HEIGHT') THEN 

      CALL initialize_plume_bpt(plume_ns, plume_np, fplume_init)
      !
      CALL initialize_plume_wind(profile_nz,profile_z,profile_rho,profile_p,profile_T,  &
                            &  profile_sh,profile_u,profile_PSIa) 
      !
      plume_z      = 0.0_wp
      CALL solve_plume_bpt(fplume_init%phase(plume_phase),plume_status,     &
                      &  plume_MER,plume_z)

!----------- Reverse z Profile again----------------------------------------

      plume_MER_icon         = plume_MER
      plume_zv_icon          = fplume_init%plume_zv 
      plume_H_icon           = plume_z(plume_ns) - fplume_init%plume_zv
      plume_MER_SO2          = fplume_init%phase(plume_phase)%MER_SO2
      plume_MER_H2O          = fplume_init%phase(plume_phase)%MER_H2O
      plume_fine_ash_fraction= fplume_init%phase(plume_phase)%plume_fine_ash_fraction
      
    ELSE IF ((fplume_init%phase(plume_phase)%plume_Hdt < fplume_min_height .AND.    &
      & fplume_init%plume_solve_for == 'MFR')) THEN


      plume_MER_icon  = (3.295_wp*fplume_init%phase(plume_phase)%plume_Hdt/1.e3_wp)  &
                  &     **(4.15_wp)

      DO i = 1, plume_ns
        plume_z(i)    = fplume_init%phase(plume_phase)%plume_Hdt/plume_ns *i
      ENDDO
      plume_zv_icon           = fplume_init%plume_zv
      plume_H_icon            = fplume_init%phase(plume_phase)%plume_Hdt
      plume_MER_SO2           = fplume_init%phase(plume_phase)%MER_SO2
      plume_MER_H2O           = fplume_init%phase(plume_phase)%MER_H2O
      plume_fine_ash_fraction = fplume_init%phase(plume_phase)%plume_fine_ash_fraction

    ENDIF 
  ELSE
    plume_MER_icon            = 0.0_wp
    plume_zv_icon             = 0.0_wp
    plume_H_icon              = 0.0_wp
    plume_MER_SO2             = 0.0_wp
  ENDIF !fplume on?
  !
  IF (plume_status == 1) THEN
    CALL finish(thisroutine,'COLLAPSE of volcanic plume')
  ELSE
    WRITE (message_text,*) 'Plume characteristics have been calculated successfully'
    CALL message (thisroutine, message_text)
  ENDIF

  RETURN
END SUBROUTINE art_fplume
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_fplume
