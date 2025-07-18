!
! mo_art_simnoy
! This module provides a simplified photochemistry for
! N2O and NOy
! first introduced by Olsen, 2001
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

MODULE mo_art_simnoy

  ! ICON
  USE mo_physical_constants,        ONLY: earth_radius
  USE mo_kind,                      ONLY: wp, i8
  USE mo_math_constants,            ONLY: pi, pi_2, rad2deg
  USE mo_exception,                 ONLY: finish, message, message_text, warning
  USE mo_impl_constants,            ONLY: MAX_CHAR_LENGTH
  USE mo_art_impl_constants,        ONLY: IART_SIMNOY_UBC
  USE mtime,                        ONLY: datetime, newDatetime, deallocateDatetime, julianday, newJulianday &
       &                                , deallocateJulianday, getJulianDayFromDatetime &
       &                                , newDateTime, deallocateDateTime, no_of_ms_in_a_day &
       &                                , getNoOfDaysInYearDateTime, getDayOfYearFromDateTime
  ! ART
  USE mo_art_read_simnoy,      ONLY:  tparm_noy,tparm_noy2, lparm, lparm2,  &
                                 &    kparm,               &
                                 &    lparm_min, lparm_dist
  USE mo_art_data,             ONLY: p_art_data
  USE mo_art_atmo_data,        ONLY: t_art_atmo
  USE mo_art_chem_types_param, ONLY: t_chem_meta_simnoy
  USE mo_art_chem_data,        ONLY: t_art_chem
!  USE mtime_datetime,          ONLY: getJulianDayFromDatetime


  IMPLICIT NONE
  PRIVATE

  PUBLIC  ::  art_noy_polarchem, art_noy_polarchem_sedi
  PUBLIC  ::  art_n2onoy_chemistry_wmo, art_n2onoy_chemistry_pres, art_n2onoy_chemistry_extp


CONTAINS


! -----------------------------------------
! --- simplified N2O / NOy chemistry
! --- Olsen et al., 2001
! --- by IMK-ASF DATE
! -----------------------------------------

SUBROUTINE art_simnoy_get_tab_values(current_date,lat,pres,nlev, jcs, jce, &
                      &              n2onoy_tab1, n2onoy_tab2,             &
                      &              n2onoy_tab3, n2onoy_tab4)
  ! inout variables
  TYPE(datetime), INTENT(in)    ::  &
    &  current_date           !< actual date
  REAL(wp),       INTENT(in)    ::  &
    &  lat(:)           !< [lat] = rad
  REAL(wp),       INTENT(in)    ::  &
    &  pres(:,:)          !< pressure of full levels
  INTEGER,        INTENT(in)    ::  &
    &  nlev, jcs, jce
  REAL(wp), INTENT(inout)         ::  &
    &  n2onoy_tab1(:,:),              &
    &  n2onoy_tab2(:,:),              &
    &  n2onoy_tab3(:,:),              &
    &  n2onoy_tab4(:,:)

  ! local variables
  INTEGER              :: jk, ik, il, ilp, il2, ilp2, jc
  REAL(wp)             :: alf, alf1, alf2, alf3, zlogp

  !----------------------------------------------------------
  ! From io/mo_art_read_simnoy
  ! lparm = 20        ! #heights
  ! lparm2 = 31       ! #heights for the extended table
  ! kparm = 18        ! #latitudes -85 .... +85
  ! mparm = 12        ! #month
  ! nparm = 5         ! #tables
  ! lparm_min = 14.   ! <lowest height with available
  !                   !  values in table, (km)
  ! lparm_dist = 2.   ! <constant distance between
  !                   !  table height levels
  !-----------------------------------------------------------

  ! Tables are given on pressure height
  ! Calculate the pressure altitude according MPrather code
  ! z_logp = 16._wp*LOG10(100000._wp/art_atmo%pres)

  DO jc = jcs,jce
    DO jk=1,nlev
      zlogp = 16._wp*LOG10(100000._wp/pres(jc,jk)) 
      ! nearest neighbour for lat
      ik = INT( (lat(jc) + pi_2) / pi * kparm ) + 1

      ! linear interpolation for lev
      il = INT( ( zlogp - lparm_min) / lparm_dist) + 1
      ilp = il + 1

      IF (il >= lparm) THEN
        il = lparm
        ilp = lparm
      ENDIF

      IF (il < 1) THEN
        il = 1
        ilp = 1
      ENDIF

      ! relative distance to il = weight of level ilp
      alf = ( zlogp - lparm_min - (il-1)*lparm_dist )/lparm_dist

      IF ((alf > 1.0_wp) .OR. (alf <= 0.0_wp)) THEN
        alf = 0.0_wp
        ilp = il
      ENDIF

      alf1 = 1.0_wp - alf

      ! linear interpolation for lev(for table2)
      il2 = INT( (zlogp - lparm_min) / lparm_dist) + 1
      ilp2 = il2 + 1

      IF (il2 >= lparm2) THEN
        il2 = lparm2
        ilp2 = lparm2
      ENDIF

      IF (il2 < 1) THEN
        il2 = 1
        ilp2 = 1
      ENDIF

      alf2 = ( zlogp - lparm_min - (il2-1)*lparm_dist )/lparm_dist

      IF ((alf2 > 1.0_wp) .OR. (alf2 <= 0.0_wp)) THEN
        alf2 = 0.0_wp
        ilp2 = il2
      ENDIF

      alf3 = 1.0_wp - alf2

      ! creating interpolated coefficient lists
      ! note extended tables 2 and 4
      n2onoy_tab1(jc,jk) = alf1 * tparm_noy(il,ik,current_date%date%month,1) &
                        & + alf * tparm_noy(ilp,ik,current_date%date%month,1)
      n2onoy_tab2(jc,jk) = alf3 * tparm_noy2(il2,ik,current_date%date%month,2) &
                        & + alf2 * tparm_noy2(ilp2,ik,current_date%date%month,2)
      n2onoy_tab3(jc,jk) = alf1  * tparm_noy(il,ik,current_date%date%month,3) &
                        & + alf * tparm_noy(ilp,ik,current_date%date%month,3)
      n2onoy_tab4(jc,jk) = alf3 * tparm_noy2(il2,ik,current_date%date%month,4) &
                        & + alf2 * tparm_noy2(ilp2,ik,current_date%date%month,4)
    ENDDO
  END DO

END SUBROUTINE art_simnoy_get_tab_values

SUBROUTINE art_n2onoy_chemistry_extp(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime)
!<
! SUBROUTINE n2onoy_chemistry
! based on Olsen 2001
! Author: Christopher Diekmann, KIT
!>
  LOGICAL, SAVE                           :: lfirst = .TRUE.
  !TYPE(datetime), POINTER                 :: valid_datetime
  REAL(wp)                                :: rjd_now
  REAL(wp), DIMENSION(250)                :: ap_int
  !LOGICAL, SAVE                           :: lread_ap = .FALSE.
  REAL(wp), SAVE, ALLOCATABLE             :: noy(:),bg(:),press(:)
  REAL(wp)                                :: zlat
  INTEGER                                 :: nlev_top
  ! inout variables
  INTEGER, INTENT(in)                     :: jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: tracer
  TYPE(datetime), POINTER, INTENT(in)     :: current_date
  REAL(wp), INTENT(in)                    :: p_dtime

  ! local variables
  REAL(wp)                                :: dn2odt, noyvmr, dnoydt, noyloss
  INTEGER                                 :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER               :: art_atmo

  art_atmo => p_art_data(jg)%atmo

  nlev = art_atmo%nlev

  !--- only for UBC-NOx
  !--- define levels where to set upper boundary condition for NOx
  nlev_top = 3
  if ( tracer%ubc == IART_SIMNOY_UBC ) then
    IF ( lfirst ) THEN
      ALLOCATE(noy(nlev_top))
      ALLOCATE(press(nlev_top))
      ALLOCATE(bg(nlev_top))
      lfirst = .FALSE.
    END IF
  endif

  ! -------------------------------------------------
  ! interpolation of coefficients onto ICON grid
  ! -------------------------------------------------
  !--- note table2 for NOy loss extended in altitude for UBC-NOx
  CALL art_simnoy_get_tab_values(current_date, art_atmo%lat(:,jb),   &
          &                      art_atmo%pres(:,:,jb), nlev,        &
          &                      jcs, jce, tracer%n2onoy_tab1, tracer%n2onoy_tab2, &
          &                      tracer%n2onoy_tab3, tracer%n2onoy_tab4)


  ! ------------------------------------------------------------
  ! start chemistry
  ! calculate new mixing ratios for NOy and N2O
  ! ------------------------------------------------------------
  DO jc = jcs,jce  
    !--- now use Simnoy tables to calculate tendenciies
    DO jk=1,nlev
      ! N2O loss:
      dn2odt = (1.0_wp - EXP(-tracer%n2onoy_tab1(jc,jk)*p_dtime)) * tracer%n2o_tracer(jc,jk,jb)

      tracer%tend_n2o(jc,jk,jb) = - dn2odt

      ! NOy production:
      tracer%tend(jc,jk,jb) = 2.0_wp * dn2odt * tracer%n2onoy_tab3(jc,jk)

      ! NOy loss:
      noyvmr = tracer%tracer(jc,jk,jb) / p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)

      dnoydt = (1.0_wp - EXP(-tracer%n2onoy_tab2(jc,jk)*p_dtime)) * tracer%tracer(jc,jk,jb)
      noyloss = 2.0_wp * noyvmr * tracer%n2onoy_tab4(jc,jk)       &
              & / (noyvmr * tracer%n2onoy_tab4(jc,jk) + 1.0)
      tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb)  - noyloss*dnoydt
    ENDDO


    !--- UBC-NOx
    if ( tracer%ubc == IART_SIMNOY_UBC ) then
    !--- UBC-NOx parameterized by ap-index
      CAll  read_apindex(current_date, ap_int,rjd_now )
    !--- get UBC-NOx
    !--- UBC-NOx depends on latitude
      zlat = art_atmo%lat(jc,jb)*rad2deg
      press(:) = art_atmo%pres(jc,:,jb)/ 100._wp
      noy = 0._wp
      CALL ubcnox_calc_current_ubc(current_date,noy,bg,nlev_top,ap_int,press(1:nlev_top),zlat)
      !--- add positive tendency from UBC-NOx in nlev-top layers
      DO jk=1,nlev_top
      !--- positive tendency from UBC-NOx only if UBC-NOx > NOy+dNOy
      !--- otherwise NOy would be destroyed by UBC-NOx
        if ( tracer%tracer(jc,jk,jb) + tracer%tend(jc,jk,jb) < noy(jk) ) then
          tracer%tend(jc,jk,jb) = noy(jk) - tracer%tracer(jc,jk,jb)
        endif
      ENDDO
    endif

  END DO
END SUBROUTINE art_n2onoy_chemistry_extp

SUBROUTINE art_n2onoy_chemistry_wmo(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
! SUBROUTINE n2onoy_chemistry
! based on Olsen 2001
! Author: Christopher Diekmann, KIT
!>
  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer
  TYPE(datetime), POINTER, INTENT(in) :: &
    &  current_date
  REAL(wp), INTENT(in) :: &
    &  p_dtime

  ! local variables
  REAL(wp)             :: dn2odt, noyvmr, &
                        & dnoydt, noyloss
  INTEGER              :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo

  nlev = art_atmo%nlev

  ! -------------------------------------------------
  ! interpolation of coefficients onto ICON grid
  ! -------------------------------------------------

  CALL art_simnoy_get_tab_values(current_date, art_atmo%lat(:,jb),   &
          &                      art_atmo%pres(:,:,jb), nlev,        &
          &                      jcs, jce, tracer%n2onoy_tab1, tracer%n2onoy_tab2, &
          &                      tracer%n2onoy_tab3, tracer%n2onoy_tab4)


  ! ------------------------------------------------------------
  ! start chemistry
  ! calculate new mixing ratios for NOy and N2O
  ! ------------------------------------------------------------
  DO jc = jcs,jce
    DO jk = 1, art_atmo%ktrpwmo(jc,jb)
      IF (jk <= art_atmo%ktrpwmo(jc,jb)) THEN
        ! N2O loss:
        dn2odt = (1.0_wp- EXP(-tracer%n2onoy_tab1(jc,jk)*p_dtime)) * tracer%n2o_tracer(jc,jk,jb)

        tracer%tend_n2o(jc,jk,jb) = - dn2odt

        ! NOy production:
        tracer%tend(jc,jk,jb) = 2.0_wp * dn2odt * tracer%n2onoy_tab3(jc,jk)

        ! NOy loss:
        noyvmr = tracer%tracer(jc,jk,jb) / p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)

        dnoydt = (1.0_wp - EXP(-tracer%n2onoy_tab2(jc,jk)*p_dtime)) * tracer%tracer(jc,jk,jb)

        noyloss = 2.0_wp * noyvmr * tracer%n2onoy_tab4(jc,jk)     &
                & / (noyvmr * tracer%n2onoy_tab4(jc,jk) + 1.0)

        tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) - noyloss*dnoydt
      ELSE
        tracer%tend_n2o(jc,jk,jb) = tracer%n2o_tracer(jc,jk,jb)           &
                &                    * EXP(-1._wp*p_dtime*tracer%des_N2O) &
                &                    - tracer%n2o_tracer(jc,jk,jb)

        tracer%tend(jc,jk,jb) = tracer%tracer(jc,jk,jb)           &
                &                * EXP(-1._wp*p_dtime*tracer%des) &
                &                - tracer%tracer(jc,jk,jb)
      ENDIF
    ENDDO
  END DO
END SUBROUTINE art_n2onoy_chemistry_wmo

SUBROUTINE art_n2onoy_chemistry_pres(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
! SUBROUTINE n2onoy_chemistry
! based on Olsen 2001
! Author: Christopher Diekmann, KIT
!>
  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer
  TYPE(datetime), POINTER, INTENT(in) :: &
    &  current_date
  REAL(wp), INTENT(in) :: &
    &  p_dtime

  ! local variables
  REAL(wp)             :: dn2odt, noyvmr, &
                        & dnoydt, noyloss
  INTEGER              :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo

  nlev = art_atmo%nlev

  ! -------------------------------------------------
  ! interpolation of coefficients onto ICON grid
  ! -------------------------------------------------

  CALL art_simnoy_get_tab_values(current_date, art_atmo%lat(:,jb),   &
          &                      art_atmo%pres(:,:,jb), nlev,        &
          &                      jcs, jce, tracer%n2onoy_tab1, tracer%n2onoy_tab2, &
          &                      tracer%n2onoy_tab3, tracer%n2onoy_tab4)


  ! ------------------------------------------------------------
  ! start chemistry
  ! calculate new mixing ratios for NOy and N2O
  ! ------------------------------------------------------------

  DO jc = jcs,jce
    DO jk = 1,nlev
      IF (art_atmo%pres(jc,jk,jb) <= 90000._wp) THEN
        ! N2O loss:
        dn2odt = (1.0_wp - EXP(-tracer%n2onoy_tab1(jc,jk)*p_dtime)) * tracer%n2o_tracer(jc,jk,jb)

        tracer%tend_n2o(jc,jk,jb) =  - dn2odt

        ! NOy production:
        tracer%tend(jc,jk,jb) =  2.0_wp * dn2odt * tracer%n2onoy_tab3(jc,jk)

        ! NOy loss:
        noyvmr = tracer%tracer(jc,jk,jb) / p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)

        dnoydt = (1.0_wp - EXP(-tracer%n2onoy_tab2(jc,jk)*p_dtime)) * tracer%tracer(jc,jk,jb)
        noyloss = 2.0_wp * noyvmr * tracer%n2onoy_tab4(jc,jk)      &
                & / (noyvmr * tracer%n2onoy_tab4(jc,jk) + 1.0)

        tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) - noyloss*dnoydt

      ELSE
        tracer%tend_n2o(jc,jk,jb) = tracer%n2o_tracer(jc,jk,jb)           &
                &                    * EXP(-1._wp*p_dtime*tracer%des_N2O) &
                &                    - tracer%n2o_tracer(jc,jk,jb)

        tracer%tend(jc,jk,jb) = tracer%tracer(jc,jk,jb)               &
                &                * EXP(-1._wp*p_dtime * tracer%des)   &
                &                - tracer%tracer(jc,jk,jb)

      ENDIF
    ENDDO
  END DO

END SUBROUTINE art_n2onoy_chemistry_pres

SUBROUTINE art_noy_polarchem_sedi(jg,jb,jcs,jce ,&
            &             tracer, p_dtime )
!<
!---------------------------------------------------------
!--- Calculate non-linear heterogenous ozone loss for polar areas
!--- Author: Christopher Diekmann, KIT
!---------------------------------------------------------
!>

  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer
!  TYPE(datetime), POINTER, INTENT(in) :: &
!    &  current_date
  REAL(wp), INTENT(in) :: &
    &  p_dtime


  ! Local variables
  REAL(wp) :: diff_n_noy, n_noy, n_noy_new
  REAL(wp) :: rate

  INTEGER         :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo
  nlev = art_atmo%nlev

  !--------------------------------------------------------------------

  DO jc = jcs, jce
    diff_n_noy = 0._wp

    DO jk=1,nlev
      n_noy =  tracer%tracer(jc,jk,jb) * 1.e6_wp

      n_noy = n_noy + diff_n_noy

      IF (jk /= nlev) THEN
          rate = tracer%des_noysed * tracer%cold_tracer(jc,jk,jb)  &
              &   /  p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)
          n_noy_new = n_noy * EXP(-1._wp*p_dtime * rate)

          diff_n_noy = n_noy - n_noy_new
      END IF

      tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) +  n_noy_new * 1.e-6_wp          &
           &                - tracer%tracer(jc,jk,jb)
    END DO
  END DO

END SUBROUTINE art_noy_polarchem_sedi

SUBROUTINE art_noy_polarchem(jg,jb,jcs,jce ,&
            &             nlev, tracer)
!<
!---------------------------------------------------------
!--- Calculate non-linear heterogenous ozone loss for polar areas
!--- Author: Christopher Diekmann, KIT
!---------------------------------------------------------
!>

  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce,    &
    &  nlev
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer

  ! Local variables
  REAL(wp) :: diff_n_noy, n_noy, n_noy_new
  REAL(wp) :: diff_vmr_cold

  INTEGER  :: jk, jc

  !--------------------------------------------------------------------
  DO jc = jcs, jce
    diff_n_noy  = 0._wp

    DO jk=1,nlev
      ! vertical redistribution of noy
      n_noy = tracer%tracer(jc,jk,jb) * 1.e6_wp

      n_noy = n_noy + diff_n_noy
  
      IF (jk /= nlev) THEN 
        !diff_vmr_cold = tracer%p_cold_sed(jc,jk,jb) * 1.e-6_wp       &
        !              &  /  p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb) 
        ! The upper part is not restartable, due to p_cold_sed not being
        ! properly recreated at the restart (needs some code restructuring)
        ! AS A PRELIMARY SOLUTION P_COLD_SED IS SET TO A CONSTANT 0.0
        ! THIS BASICALLY TURNS THIS PROCESS OFF. THIS NEEDS PROPER
        ! TREATMENT BEFORE IT CAN BE ACTIVATED AGAIN.
        ! IF YOU CAN READ THIS MESSAGE, THIS WAS NOT YET DONE.
        diff_vmr_cold = 0.0_wp * 1.e-6_wp       &
                      &  /  p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb) 
        n_noy_new     = n_noy * (1._wp - diff_vmr_cold)

        diff_n_noy    = n_noy - n_noy_new
      END IF

      tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) + n_noy_new * 1.e-6_wp   &
           &                  - tracer%tracer(jc,jk,jb)
    ENDDO
  END DO

END SUBROUTINE art_noy_polarchem


SUBROUTINE read_apindex( valid_datetime, ap_block, rjd_now )
!< driver for reading the daily apindex for calculaton of upper boundary NOx
! reads a block of values to reduce io-overhead
!>
  !---
  TYPE(datetime), POINTER  :: valid_datetime
  INTEGER, PARAMETER        :: nwidth=250 ! span in days of used ap index
  REAL(wp), INTENT(inout)   :: ap_block(nwidth) ! ap indices for use in parameterization

  !CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !                routine='read_apindex'

  !---LOCAL VARIABLES
  TYPE(julianday), POINTER :: jd_now
  INTEGER, PARAMETER             :: ntime=2000
  REAL(wp), SAVE                 :: time_ref(ntime), apref(ntime)
  INTEGER :: errno, iref1
  REAL(wp), SAVE                 :: rjd_start
  REAL(wp)                       :: rjd_now
  LOGICAL, SAVE :: lfirst = .TRUE.

  !--- get the julian date of current date
  !    prepare c-pointers according mtime lib
  jd_now => newJulianday(0_i8,0_i8,errno)
  CALL getJulianDayFromDatetime(valid_datetime, jd_now)
  !- julian day now
  rjd_now = jd_now%day + jd_now%ms/86400._wp/1000._wp
  CALL deallocateJulianday( jd_now )
  !> start julian date of current list
  IF ( lfirst ) THEN
    rjd_start = 0._wp
    lfirst = .FALSE.
  ENDIF

  IF ( rjd_now - rjd_start > ntime - nwidth ) THEN
      !> read time series ap
      CALL ubcnox_read_ap( rjd_now-nwidth, time_ref, apref, ntime )
      rjd_start = rjd_now
  END IF
  iref1 = INT(rjd_now-rjd_start) + 1
  
  ap_block = apref(iref1-nwidth:iref1-1)
 
END SUBROUTINE read_apindex

SUBROUTINE ubcnox_read_ap( rjd_0, ts, ap, ntime)
!> @brief reads a limited number of ap indices for ubcnox boundary condition
!!
!! The routine reads a specific number of ap indices for the use of a NOy upper boundary.
!! Data are from CMIP time series, starting from 1850-01-01. 
!! Julian dates and ap indices starting from rjd_now are returned.
!! 
!! @param[in,out]  ts list of days of type wp
!! @param[in,out]  ap corresponding list of ap indices of type wp
!! @param[in] ntime length of list type integer
!! @param[in] rjd_now of type wp
!! @par Revision History
!! First version by Maryam (<2023-02-01>)
!! 
!!
  INTEGER,  INTENT(IN)            :: ntime
  REAL(wp), INTENT(INOUT)         :: ts(ntime),ap(ntime)
  REAL(wp), INTENT(IN)            :: rjd_0

  !--- the reference julian date in the ap file
  REAL(wp), PARAMETER :: rjd_ref = 2396759._wp ! 1850-01-01
  CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
                routine='ubcnox_read_ap'
  CHARACTER(LEN=256)             :: cfname   ! file name containing variables

  INTEGER                        :: ird, i, istat
  REAL(wp)                       :: rdum

  ts = 0._wp
  ap = 0._wp

  cfname = 'ap_daily.txt'
  !< index where to start list
  ird = int( rjd_0 - rjd_ref )
  IF ( ird < 1 ) THEN
    CALL message( routine, 'Outside available data. ap set to 0.' )
    RETURN
  ENDIF

  OPEN( unit=45, file=cfname, status='old' )
  DO i=1,ird
    READ(45,*) rdum, rdum
  ENDDO

  READ(45,*, IOSTAT=istat) ( ts(i), ap(i), i=1,ntime )

  IF ( istat > 0 ) THEN
    CALL finish( routine, 'Error reading file '//cfname )
  ELSE IF ( istat < 0 ) THEN
    ts = 0._wp
    ap = 0._wp
    CALL message( routine, 'End of ap list reached. Use zero for ap.' )
  ENDIF

  CLOSE(45)

  CALL message( routine, 'Ap index list read.' )

END SUBROUTINE ubcnox_read_ap
!-----------------------------------------------------------

  SUBROUTINE ubcnox_linear_interp(interp,xq,x1,x2,y1,y2)
! linear interpolation between two points
   REAL(wp), INTENT(OUT) :: interp         ! interpolated value
   REAL(wp), INTENT(IN) :: xq              ! x value for interpolated value
   REAL(wp), INTENT(IN) :: x1,x2,y1,y2     ! two points for interpolation

   REAL(wp) :: dx                          ! total distance for x
   REAL(wp) :: w1, w2                      ! weights for point1 and point2

   dx=x2-x1
   IF (abs(dx)<1.e-20_wp) THEN
     interp=y1
   ELSE
     w1=1._wp-(xq-x1)/dx
     w2=1._wp-(x2-xq)/dx
     interp=w1*y1+w2*y2
   END IF

END SUBROUTINE


!-------------------------------------------------------------------------
SUBROUTINE time_span_d(dtd,valid_datetime1,valid_datetime2)
!< calculates time span in days between two dates
!>
  TYPE(datetime), POINTER :: &
&   valid_datetime1,valid_datetime2
  INTEGER :: errno 
  REAL(wp) :: rjd_now2, rjd_now1, dtd ! dime span [d]
  TYPE(julianday), POINTER :: jd_now1,jd_now2
  jd_now1 => newJulianday(0_i8,0_i8,errno)
  jd_now2 => newJulianday(0_i8,0_i8,errno)
  CALL getJulianDayFromDatetime(valid_datetime1, jd_now1)
  CALL getJulianDayFromDatetime(valid_datetime2, jd_now2)
  rjd_now1 = jd_now1%day + jd_now1%ms/86400._wp/1000._wp
  rjd_now2 = jd_now2%day + jd_now2%ms/86400._wp/1000._wp
  dtd = (rjd_now2 - rjd_now1)
  CALL deallocateJulianday( jd_now1 )
  CALL deallocateJulianday( jd_now2 )

END SUBROUTINE time_span_d
!-------------------------------------------------------------------------


SUBROUTINE ubcnox_calc_current_ubc(current_date,noy,bg,nlev_top,ap_int,press,zlat)
  TYPE(datetime), POINTER, INTENT(in)  :: &
  &  current_date
  INTEGER, INTENT(in)                   :: nlev_top
  REAL(wp), INTENT(in)                  :: zlat     !< latitude (deg)i
  REAL(wp), INTENT(out)                 :: noy(nlev_top)   ! EPP-noy field (complete)
  REAL(wp), INTENT(out)                 :: bg(nlev_top)    ! EPP-noy field (only background)
  REAL(wp), INTENT(in)                  :: press(nlev_top) ! Pressure in hPascal
  REAL(wp), DIMENSION(250)              :: ap_int          ! Ap Index
 

  ! local variables
  LOGICAL, SAVE                         :: lfirst = .TRUE.
  INTEGER  :: jk
  REAL(wp), PARAMETER                   :: N_A = 6.02214129E23_wp
  INTEGER, PARAMETER                    :: nlrefh = 9, npref = 12
  REAL(wp), SAVE, DIMENSION(250)        :: dl
  REAL(wp), SAVE                        :: pref(npref)
  REAL(wp), SAVE                        :: txnr(npref)
  REAL(wp), SAVE                        :: txsr(npref)
  REAL(wp), SAVE                        :: lref(nlrefh*2)
  REAL(wp), SAVE                        :: latbdlo(nlrefh*2)
  REAL(wp), SAVE                        :: latbdhi(nlrefh*2)
  REAL(wp), SAVE                        :: lrefs(nlrefh)
  REAL(wp), SAVE                        :: lrefn(nlrefh)
  REAL(wp), SAVE                        :: tmnr(npref)
  REAL(wp), SAVE                        :: tmsr(npref)
  REAL(wp), SAVE                        :: tfnr(npref)
  REAL(wp), SAVE                        :: tfsr(npref)
  REAL(wp), SAVE                        :: nmnr(npref)
  REAL(wp), SAVE                        :: nmsr(npref)
  REAL(wp), SAVE                        :: wmnr(npref)
  REAL(wp), SAVE                        :: wmsr(npref)
  REAL(wp), SAVE                        :: wfnr(npref)
  REAL(wp), SAVE                        :: wfsr(npref)
  REAL(wp), SAVE                        :: lds(npref,nlrefh)
  REAL(wp), SAVE                        :: ldn(npref,nlrefh)
  REAL(wp), SAVE                        :: lde(npref,nlrefh)
  REAL(wp), SAVE                        :: amr(npref,nlrefh*2)
  REAL(wp), SAVE                        :: a1r(npref,nlrefh*2)
  REAL(wp), SAVE                        :: a2r(npref,nlrefh*2)
  REAL(wp), SAVE                        :: a3r(npref,nlrefh*2)
  REAL(wp), SAVE                        :: p1r(npref,nlrefh*2)
  REAL(wp), SAVE                        :: p2r(npref,nlrefh*2)
  REAL(wp), SAVE                        :: p3r(npref,nlrefh*2)
  REAL(wp), SAVE                        :: latweight(nlrefh*2)
  REAL(wp), SAVE                        :: area(nlrefh*2)
  
  REAL(wp), DIMENSION(nlrefh)           :: hs, hn

  REAL(wp) :: dn,ds          ! days since Jan, 1st (ds) and Jul, 1st (dn)
  REAL(wp) :: tn, ts, nn, ns, wn, ws, seasn, seass, txn, txs, xn, xs, &
              facs, facn
  REAL(wp), DIMENSION(250) :: filtern, filters, filtern_flip, filters_flip
  !REAL(wp) :: lpl, tm, xu, wu, fm, wm, xb, nne, we, rfac, tl, sease, xe           ! variables for ESE
  REAL(wp)  ::  xld
  ! local temporary variables
  INTEGER  :: i, j,l,k
  REAL(wp) :: ep(nlev_top)    ! EPP-noy field (without background)
  REAL(wp) :: tmnr1
  REAL(wp) :: tmnr2
  REAL(wp) :: pref1
  REAL(wp) :: pref2
  REAL(wp) :: lat1
  REAL(wp) :: lat2
  REAL(wp) :: fac   ! Factor 
  LOGICAL, DIMENSION(npref) :: mask
  INTEGER :: ind1, ind2, indlat1, indlat2
  REAL(wp), DIMENSION(nlrefh*2) :: a1,a2,a3,p1,p2,p3,am,bgh
  TYPE(datetime), POINTER :: valid_datetime1,valid_datetime2,valid_datetime3 


  IF ( lfirst ) THEN

    CALL read_ubcnox( npref, nlrefh, pref, lref, lrefs, lrefn, latbdlo, latbdhi &
                      & , txnr, txsr, tmnr, tmsr, tfnr, tfsr &
                      & , nmnr, nmsr, wmnr, wmsr, wfnr, wfsr, lds, ldn, lde &
                      & , amr, a1r, a2r, a3r, p1r, p2r, p3r )
    DO i=2,250
      dl(i)=i-1
    END DO
    dl(1) = 0.5_wp
    WHERE (lds<1.e-20_wp) lds = 1.e-20_wp
    WHERE (ldn<1.e-20_wp) ldn = 1.e-20_wp
    WHERE (lde<1.e-20_wp) lde = 1.e-20_wp

    !--- scale ld to density
    latweight = SIN(pi/180._wp*latbdhi)-SIN(pi/180._wp*latbdlo)       ! associated area weights
    area=latweight*2*pi*earth_radius*earth_radius*1.e4_wp     ! associated area in cm2 (note R_E in Meter)
    !--- GM -> density
    DO k=1,npref
      facs = 0._wp
      facn = 0._wp
      DO l=1,nlrefh
          facs = facs + lds(k,l)*area(l)/N_A
          facn = facn + ldn(k,l)*area(l+nlrefh)/N_A
      ENDDO
      lds(k,:) = lds(k,:)/facs*1.e4_wp
      ldn(k,:) = ldn(k,:)/facn*1.e4_wp
    ENDDO
    lfirst = .FALSE.
  ENDIF

  noy = 0._wp
  bg = 0._wp

  valid_datetime1 => newDatetime(current_date%date%year, 1, 1 , 0, 0, 0, 0)
  valid_datetime2 => newDatetime(current_date%date%year, current_date%date%month, current_date%date%day , 0, 0, 0, 0)
  IF ( current_date%date%month <7 ) THEN
    CALL time_span_d (dn,valid_datetime1,valid_datetime2)
    dn= dn+185
  else
    valid_datetime3 => newDatetime(current_date%date%year, 7, 1 , 0, 0, 0, 0)
    CALL time_span_d(dn,valid_datetime3,valid_datetime2)
    dn=dn+1
  end if
  CALL time_span_d (ds,valid_datetime1,valid_datetime2)
  CALL deallocateDateTime( valid_datetime1 )
  CALL deallocateDateTime( valid_datetime2 )
  CALL deallocateDateTime( valid_datetime3 )
  ds= ds+1

  DO jk=1,nlev_top
    !--- determine nearest ref pressure indices
    !    valid for whole loop!
    mask=(/ .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE./)
    ind1=MINLOC(ABS(pref-press(jk)),DIM=1,MASK=mask)
    mask(ind1)=.FALSE.
    pref1=pref(ind1)
    ind2=MAX(ind1-1,1)
    IF (pref1 > press(jk)) THEN
    ind2=MIN(ind1+1,12)
    END IF
    pref2=pref(ind2)
    tmnr1=tmnr(ind1)
    tmnr2=tmnr(ind2)
    CALL ubcnox_linear_interp(tn,press(jk),pref1,pref2,tmnr1,tmnr2)
    CALL ubcnox_linear_interp(ts,press(jk),pref1,pref2,tmsr(ind1),tmsr(ind2))
    CALL ubcnox_linear_interp(nn,press(jk),pref1,pref2,nmnr(ind1),nmnr(ind2))
    CALL ubcnox_linear_interp(ns,press(jk),pref1,pref2,nmsr(ind1),nmsr(ind2))
    CALL ubcnox_linear_interp(wn,press(jk),pref1,pref2,wmnr(ind1),wmnr(ind2))
    CALL ubcnox_linear_interp(ws,press(jk),pref1,pref2,wmsr(ind1),wmsr(ind2))
    seasn=4._wp*nn*EXP(-wn*(dn-tn))/(1._wp+EXP(-wn*(dn-tn)))**2
    seass=4._wp*ns*EXP(-ws*(ds-ts))/(1._wp+EXP(-ws*(ds-ts)))**2

    ind1=MINLOC(ABS(pref-press(jk)),DIM=1)
    pref1=pref(ind1)
    ind2=MAX(ind1-1,1)
    IF (pref1 > press(jk)) THEN
    ind2=MIN(ind1+1,12)
    END IF
    pref2=pref(ind2)

    DO j=1,9
      CALL ubcnox_linear_interp(hs(j),press(jk),pref1,pref2,lds(ind1,j),lds(ind2,j))
      CALL ubcnox_linear_interp(hn(j),press(jk),pref1,pref2,ldn(ind1,j),ldn(ind2,j))
    END DO

    IF (zlat<0._wp) THEN
      indlat1=MINLOC(ABS(lrefs-zlat),DIM=1)
      IF (zlat<lrefs(indlat1)) THEN
          indlat2=MAX(indlat1-1,1)
      END IF
      IF (zlat>=lrefs(indlat1)) THEN
          indlat2=MIN(indlat1+1,9)
      END IF
      lat1=lrefs(indlat1)
      lat2=lrefs(indlat2)
      CALL ubcnox_linear_interp(xld,zlat,lat1,lat2,hs(indlat1),hs(indlat2))
      IF (zlat<lrefs(1)) THEN
          xld=hs(1)
      END IF
      IF (zlat>lrefs(9)) THEN
          xld=hs(9)
      END IF
    END IF
    IF (zlat>=0._wp) THEN
      indlat1=MINLOC(ABS(lrefn-zlat),DIM=1)
      IF (zlat<lrefn(indlat1)) THEN
          indlat2=MAX(indlat1-1,1)
      END IF
      IF (zlat>=lrefn(indlat1)) THEN
          indlat2=MIN(indlat1+1,9)
      END IF
      lat1=lrefn(indlat1)
      lat2=lrefn(indlat2)
      CALL ubcnox_linear_interp(xld,zlat,lat1,lat2,hn(indlat1),hn(indlat2))
      IF (zlat>lrefn(9)) THEN
          xld=hn(9)
      END IF
      IF (zlat<lrefn(1)) THEN
          xld=hn(1)
      END IF
    END IF
    fac = xld

    ! filter function for calculation of weighted Ap
    CALL ubcnox_linear_interp(txn,press(jk),pref1,pref2,txnr(ind1),txnr(ind2))
    CALL ubcnox_linear_interp(txs,press(jk),pref1,pref2,txsr(ind1),txsr(ind2))
    filtern(1:250)=SQRT(1._wp/(dl(1:250)*dl(1:250)*dl(1:250))) &
                & *EXP(-(dl(1:250)-txn) * (dl(1:250)-txn)/(2*(SQRT(0.7_wp*txn)+6._wp)**2 * dl(1:250)/txn))
    filters(1:250)=SQRT(1._wp/(dl(1:250)*dl(1:250)*dl(1:250))) &
                & *EXP(-(dl(1:250)-txs) * (dl(1:250)-txs)/(2*(SQRT(0.7_wp*txs)+6._wp)**2 * dl(1:250)/txs))
    filtern=filtern/sum(filtern(:))
    filters=filters/sum(filters(:))
    ! flip indices for filtern and filters
    filtern_flip(1:250)=filtern(250:1:-1)
    filters_flip(1:250)=filters(250:1:-1)
    xn=sum(ap_int*filtern_flip(1:250))*seasn
    xs=sum(ap_int*filters_flip(1:250))*seass
    IF (zlat<0._wp) THEN
      ep(jk)=xs*fac
    END IF
    IF (zlat>=0._wp) THEN
      ep(jk)=xn*fac
    END IF

    ! calculate background NOy
    indlat1=minloc(abs(lref-zlat),DIM=1)
    IF (zlat<lref(indlat1)) THEN
        indlat2=MAX(indlat1-1,1)
    END IF
    IF (zlat>=lref(indlat1)) THEN
        indlat2=MIN(indlat1+1,18)
    END IF
    lat1=lref(indlat1)
    lat2=lref(indlat2)

    DO j=1,18
      CALL ubcnox_linear_interp(am(j),press(jk),pref1,pref2,amr(ind1,j),amr(ind2,j))
      CALL ubcnox_linear_interp(a1(j),press(jk),pref1,pref2,a1r(ind1,j),a1r(ind2,j))
      CALL ubcnox_linear_interp(a2(j),press(jk),pref1,pref2,a2r(ind1,j),a2r(ind2,j))
      CALL ubcnox_linear_interp(a3(j),press(jk),pref1,pref2,a3r(ind1,j),a3r(ind2,j))
      CALL ubcnox_linear_interp(p1(j),press(jk),pref1,pref2,p1r(ind1,j),p1r(ind2,j))
      CALL ubcnox_linear_interp(p2(j),press(jk),pref1,pref2,p2r(ind1,j),p2r(ind2,j))
      CALL ubcnox_linear_interp(p3(j),press(jk),pref1,pref2,p3r(ind1,j),p3r(ind2,j))
      bgh(j)=am(j)*(1._wp+a1(j)*sin(ds/365._wp*2._wp*pi+p1(j))+a2(j)*sin(ds/365._wp*4._wp*pi+p2(j))+a3(j)*sin(ds/365._wp*6._wp*pi+p3(j)))
    END DO
    CALL ubcnox_linear_interp(bg(jk),zlat,lat1,lat2,bgh(indlat1),bgh(indlat2))
    IF (zlat>lref(18)) THEN
      bg(jk)=bgh(18)
    END IF
    IF (zlat<lref(1)) THEN
      bg(jk)=bgh(1)
    END IF
    IF ((zlat<=5._wp) .and. (zlat>=0._wp)) THEN
      bg(jk)=bgh(10)
    END IF
    IF ((zlat>=-5._wp) .and. (zlat<0._wp)) THEN
      bg(jk)=bgh(9)
    END IF

  END DO ! top_levels

  noy=ep+bg

END SUBROUTINE ubcnox_calc_current_ubc

SUBROUTINE read_ubcnox( npref, nlrefh, pref, lref, lrefs, lrefn, latbdlo, latbdhi &
                      & , txnr, txsr, tmnr, tmsr, tfnr, tfsr &
                      & , nmnr, nmsr, wmnr, wmsr, wfnr, wfsr, lds, ldn, lde &
                      & , amr, a1r, a2r, a3r, p1r, p2r, p3r )


  INTEGER  :: npref, nlrefh
  REAL(wp) :: pref(npref)
  REAL(wp) :: lref(nlrefh*2)
  REAL(wp) :: txnr(npref)
  REAL(wp) :: txsr(npref)
  REAL(wp) :: latbdlo(nlrefh*2)
  REAL(wp) :: latbdhi(nlrefh*2)
  REAL(wp) :: lrefs(nlrefh)
  REAL(wp) :: lrefn(nlrefh)
  REAL(wp) :: tmnr(npref)
  REAL(wp) :: tmsr(npref)
  REAL(wp) :: tfnr(npref)
  REAL(wp) :: tfsr(npref)
  REAL(wp) :: nmnr(npref)
  REAL(wp) :: nmsr(npref)
  REAL(wp) :: nfnr(npref)
  REAL(wp) :: nfsr(npref)
  REAL(wp) :: wmnr(npref)
  REAL(wp) :: wmsr(npref)
  REAL(wp) :: wfnr(npref)
  REAL(wp) :: wfsr(npref)
  REAL(wp) :: lds(npref,nlrefh)
  REAL(wp) :: ldn(npref,nlrefh)
  REAL(wp) :: lde(npref,nlrefh)
  REAL(wp) :: amr(npref,nlrefh*2)
  REAL(wp) :: a1r(npref,nlrefh*2)
  REAL(wp) :: a2r(npref,nlrefh*2)
  REAL(wp) :: a3r(npref,nlrefh*2)
  REAL(wp) :: p1r(npref,nlrefh*2)
  REAL(wp) :: p2r(npref,nlrefh*2)
  REAL(wp) :: p3r(npref,nlrefh*2)
  

  INTEGER :: iubcunit, ios, i

  !read the data 
  OPEN(NEWUNIT=iubcunit, file='ubcnox_data.dat', status='old', iostat=ios )
  IF ( ios /= 0 ) THEN
    CALL finish( 'read_ubc_nox', 'Error reading file ubcnox_data.dat' )
  ENDIF

  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) pref
  READ(iubcunit,*) lref
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) latbdlo
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) latbdhi
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) lrefs
  READ(iubcunit,*) lrefn
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) txnr
  READ(iubcunit,*) txsr
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) tmnr
  READ(iubcunit,*) tmsr
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) tfnr
  READ(iubcunit,*) tfsr
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) nmnr
  READ(iubcunit,*) nmsr
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) nfnr
  READ(iubcunit,*) nfsr
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) wmnr
  READ(iubcunit,*) wmsr
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) wfnr
  READ(iubcunit,*) wfsr

  READ(iubcunit,*) ! skip line
  DO i=1,9 
    READ(iubcunit,*) lds(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,9 
    READ(iubcunit,*) ldn(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,9 
    READ(iubcunit,*) lde(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  READ(iubcunit,*) ! skip line
  DO i=1,18 
    READ(iubcunit,*) amr(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,18 
    READ(iubcunit,*) a1r(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,18 
    READ(iubcunit,*) a2r(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,18 
    READ(iubcunit,*) a3r(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,18 
    READ(iubcunit,*) p1r(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,18 
    READ(iubcunit,*) p2r(:,i)
  ENDDO
  READ(iubcunit,*) ! skip line
  DO i=1,18 
    READ(iubcunit,*) p3r(:,i)
  ENDDO
  
  CLOSE(iubcunit)

END SUBROUTINE read_ubcnox

END MODULE mo_art_simnoy
