!
! mo_art_aero_optical_props
! This module provides routines to calculate optical properties
! of aerosol at AERNONET and Ceilometer wavelengths
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

MODULE mo_art_aero_optical_props
! ICON
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: SUCCESS, MAX_CHAR_LENGTH
! ART
  USE mo_art_config,                    ONLY: art_config, IART_PATH_LEN
  USE mo_art_data,                      ONLY: t_art_data
  USE mo_art_modes,                     ONLY: t_fields_2mom,             &
    &                                         t_diag_optprops
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  USE mo_art_diag_types,                ONLY: t_art_aeronet,             &
    &                                         t_art_ceilo
  USE mo_art_read_opt_props,            ONLY: art_read_aeronet_optprops,              &
    &                                         art_read_satellite_optprops,            &
    &                                         art_read_singlewave_optprops
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_calc_aod
  PUBLIC :: art_calc_aodvar
  PUBLIC :: art_calc_bsc

CONTAINS

!
! private helper subroutines
!

INTEGER FUNCTION art_get_tracer_index(p_art_data, tracer_name) RESULT(i_tracer)
  !arguments
  TYPE(t_art_data), INTENT(IN) :: p_art_data
  CHARACTER(LEN=*), INTENT(IN) :: tracer_name
  !local variables
  INTEGER :: ierror

  CALL p_art_data%dict_tracer%get(TRIM(tracer_name), i_tracer, ierror)
  IF (ierror /= SUCCESS) i_tracer = 0

END FUNCTION art_get_tracer_index

SUBROUTINE art_find_diagnostic_optprops_for_tracer(jg, mode_name, diag_opt_props, success)
!<
! SUBROUTINE art_find_diagnostic_optprops_for_tracer
! This routine searches through the modes for the one with the given tracer names and
! sets the corresponding optical properties to the output variable 'diag_opt_props'
! If the desired mode is not found ICON finishes with an error message.
! Part of module: mo_art_aero_optical_props
! Author: Enrico P. Metzner, KIT
! Initial Release: 2024-11-20
!>
  !arguments
  INTEGER, INTENT(IN)                :: jg !< patch id
  CHARACTER(LEN=*), INTENT(IN)       :: mode_name
  TYPE(t_diag_optprops), INTENT(OUT) :: diag_opt_props
  LOGICAL, INTENT(OUT)               :: success
  !local variables
  TYPE(t_mode), POINTER              :: current_mode !< pointer to loop through mode structure

  success = .FALSE. !< first assume correct mode is not present
  current_mode => p_mode_state(jg)%p_mode_list%p%first_mode
  DO WHILE(ASSOCIATED(current_mode))
    IF (TRIM(current_mode%fields%name) == TRIM(mode_name)) THEN
      SELECT TYPE (fields=>current_mode%fields)
        CLASS IS (t_fields_2mom)
          diag_opt_props = fields%diag_opt_props
          success = .TRUE. !< found the correct mode
        CLASS DEFAULT
          CALL message('mo_art_aero_optical_props:art_find_diagnostic_optprops_for_tracer', &
            &          'mode '//TRIM(mode_name)//' is not of kind "2mom". '//&
            &          'Cannot return any diagnostic optprops!')
      END SELECT
      EXIT
    ENDIF
    current_mode => current_mode%next_mode
  ENDDO

END SUBROUTINE art_find_diagnostic_optprops_for_tracer



SUBROUTINE art_calc_aod_single_diag(diag, istart, iend, nlev, jg, jb,         &
    &                               tracer, dz, rho, n_wavel,                 &
    &                               n_modes, tracer_idx, mode_names, success)
!<
! SUBROUTINE art_calc_aod_single_diag
! This subroutine calculates the aerosol optical depth at 9 different wavelengths used by AERONET or only at 550nm
! for a single diagnostic. It can combine the extinction properties of up to 4 different tracers/modes, which are
! given by `tracernameX` and their corresponding index `imodeX` for the `tracer` array.
! It's a rewrite of the code developed by Daniel Rieger, Philipp Gasch and Carolin Walter (2014 - 2015)
!
! Part of Module: mo_art_aero_optical_props
! Author: Enrico P. Metzner, KIT
! Initial Release: 2024-11-20
! Modifications:
! yyyy-mm-dd: name, KIT
! - description
!>
  !arguments
  TYPE(t_art_aeronet), INTENT(INOUT) :: diag(:)                        !< diagnostics container
  LOGICAL, INTENT(INOUT) ::             success                        !< success of this routine
  INTEGER, INTENT(IN) ::                istart, iend,                & !< Start and end index of nproma loop
    &                                   nlev, jb,                    & !< Number of verical levels, Block index
    &                                   jg,                          & !< Patch id
    &                                   n_wavel,                     & !< number of wavelengths
    &                                   n_modes                        !< number of modes
  INTEGER, INTENT(IN) ::                tracer_idx(:)                  !< corresponding tracer indices
  REAL(wp), INTENT(IN) ::               rho(:,:,:),                  & !< Air density
    &                                   tracer(:,:,:,:),             & !< Tracer mixing ratios [kg/kg]
    &                                   dz(:,:,:)                      !< Layer height
  CHARACTER(LEN=*), INTENT(IN) ::       mode_names(:)                  !< mode names
  !local variables
  LOGICAL ::                        available(n_modes)
  INTEGER ::                        i, jk, jc, i_wavel, n_successes !, i550nm
  REAL(wp) ::                       ext_value
  TYPE(t_diag_optprops) ::          modediag(n_modes)

  n_successes = 0
  DO i = 1, n_modes
    IF (tracer_idx(i) > 0) THEN
      CALL art_find_diagnostic_optprops_for_tracer(jg, mode_names(i), modediag(i), available(i))
      IF (available(i)) THEN
        n_successes = n_successes + 1
      END IF
    ELSE
      available(i) = .FALSE.
    END IF
  END DO

  IF (n_successes == 0) THEN
    !TODO: add info, which diagnostic is here unsufficient established.
    CALL message("mo_art_aero_opt_props:art_calc_aod_single_diag", &
      &          "Could not find any mode contributing to this diagnostic.")
    success = .FALSE.
    return
  END IF

  IF (n_wavel == 1) THEN
    DO i_wavel = 1, n_wavel
      IF (ASSOCIATED(diag(i_wavel)%tau)) THEN
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_value = 0.0_wp;
            DO i = 1,n_modes
              IF (available(i)) THEN
                ext_value = ext_value + modediag(i)%ext_snglwave(i_wavel) &
                  &                   * tracer(jc,jk,jb,tracer_idx(i))
              ENDIF
            END DO
            !factor 1.e-6_wp converts tracer from ug/kg to g/kg as extinction is given in nm2/g
            ext_value = ext_value * rho(jc,jk,jb) * 1.e-6_wp
            diag(i_wavel)%tau(jc,jk,jb) = ext_value * dz(jc,jk,jb)
          ENDDO !jc
        ENDDO !jk
      ENDIF
    ENDDO ! i_wavel
  END IF

  IF (n_wavel == 9) THEN
    DO i_wavel = 1, n_wavel
      IF (ASSOCIATED(diag(i_wavel)%tau)) THEN
        DO jk = 1, nlev
          DO jc = istart, iend
            ext_value = 0.0_wp;
            DO i = 1,n_modes
              IF (available(i)) THEN
                ext_value = ext_value + modediag(i)%ext_aeronet(i_wavel) &
                  &                   * tracer(jc,jk,jb,tracer_idx(i))
              ENDIF
            END DO
            !factor 1.e-6_wp converts tracer from ug/kg to g/kg as extinction is given in nm2/g
            ext_value = ext_value * rho(jc,jk,jb) * 1.e-6_wp
            diag(i_wavel)%tau(jc,jk,jb) = ext_value * dz(jc,jk,jb)
          ENDDO !jc
        ENDDO !jk
      ENDIF
    ENDDO ! i_wavel
  END IF

  success = .TRUE.

END SUBROUTINE art_calc_aod_single_diag

SUBROUTINE art_calc_aodvar_aeronet(aeronet, istart, iend, nlev, jb, jg, &
  &                                n_wavel,                             &
  &                                ini_cmd, cmd, rho, dz, tracer,       &
  &                                tracer_idx, mode_name)
!<
! SUBROUTINE art_calc_aodvar_aeronet
! This subroutine calculates the aerosol optical depth at 9 different wavelengths used by AERONET for one mode
! with polynomial variable diameter approximation of the optical properties.
! The name and index within the tracer array are given by tracer_name and imode respectively.
! It's a rewrite of the code developed by Ali Hoshyaripour (2018)
!
! Part of Module: mo_art_aero_optical_props
! Author: Enrico P. Metzner, KIT
! Initial Release: 2024-11-20
! Modifications:
! yyyy-mm-dd: name, KIT
! - description
!>
  !arguments
  TYPE(t_art_aeronet), INTENT(INOUT) :: aeronet(:)                     !< diagnostics container
  INTEGER, INTENT(IN)                :: istart, iend,                & !< Start and end index of nproma loop
    &                                   nlev, jb,                    & !< Number of verical levels, Block index
    &                                   jg,                          & !< Patch id
    &                                   n_wavel,                     & !< number of wavelengths
    &                                   tracer_idx                     !< tracer container index
  REAL(wp), INTENT(IN)               :: ini_cmd,                     & !< Initial Diameter number conc.
    &                                   cmd(:,:),                    & !< Diameter with respect to number conc.
    &                                   rho(:,:,:),                  & !< Air density
    &                                   dz(:,:,:),                   & !< Layer height
    &                                   tracer(:,:,:,:)                !< Tracer mixing ratios [kg/kg]
  CHARACTER(LEN=*), INTENT(IN)       :: mode_name
  !local variables
  LOGICAL                        :: available
  INTEGER                        :: jk, jc, i_wavel
  REAL(wp)                       :: fac_cmd
  TYPE(t_diag_optprops)          :: diag_opt_props

  IF (tracer_idx < 1) THEN
    RETURN
  ENDIF

  available = .FALSE.
  CALL art_find_diagnostic_optprops_for_tracer(jg, mode_name, diag_opt_props, available)
  IF (.NOT. available) THEN
    CALL message('mo_art_aero_optical_props:art_calc_oadvar_aeronet', &
      &          'Cannot find optical properties of required mode for aod diagnostic with variable diameter modes.')
    RETURN
  END IF

  DO i_wavel = 1, n_wavel
    IF (ASSOCIATED(aeronet(i_wavel)%tau)) THEN
      DO jk = 1, nlev
        DO jc = istart, iend
          ! for 550 nm only
          IF (i_wavel==5 .AND. cmd(jc,jk) > 0.25_wp * ini_cmd) THEN
            fac_cmd = ini_cmd / cmd(jc,jk)
          ELSE
            fac_cmd = 1.0_wp
          ENDIF
          !factor 1.e-6_wp converts tracer from ug/kg to g/kg as extinction is given in nm2/g
          aeronet(i_wavel)%tau(jc,jk,jb) = aeronet(i_wavel)%tau(jc,jk,jb)             &
            &                            + diag_opt_props%ext_aeronet(i_wavel)        &
            &                            * tracer(jc,jk,jb,tracer_idx)                &
            &                            * fac_cmd * rho(jc,jk,jb) * 1.e-6_wp * dz(jc,jk,jb)
        ENDDO !jc
      ENDDO !jk
    ENDIF
  ENDDO ! i_wavel

END SUBROUTINE art_calc_aodvar_aeronet

SUBROUTINE art_calc_single_backscatter(ceilometer, attenuation, satellite,            &
  &                                    istart, iend, nlev, jg, jb, tracer, dz, rho,   &
  &                                    n_modes, tracer_idx, mode_names )
!<
! SUBROUTINE art_calc_single_backscatter
! This subroutine calculates the aerosol backscatter and attenuated backscattter at 3 different
! wavelengths used by lidars from ground and from satellite
! It combines the extinction and backscatter properties of multiple tracers/modes, which are
! given by `mode_names` and their corresponding indices `tracer_idx` for the `tracer` array.
! It's a rewrite of the code developed by D. Rieger, P. Gasch and C. Walter (2014 - 2015)
!
! Part of Module: mo_art_aero_optical_props
! Author: Enrico P. Metzner, KIT
! Initial Release: 2025-02-11
! Modifications:
!>
  !arguments
  TYPE(t_art_ceilo), INTENT(INOUT) ::   ceilometer(:),               & !< Diagnostics container for ceilometer
    &                                   attenuation(:),              & !< Diagnostics container for attenuation
    &                                   satellite(:)                   !< Diagnostics container for satellite
  INTEGER, INTENT(IN) ::                istart, iend,                & !< Start and end index of nproma loop
    &                                   nlev, jb,                    & !< Number of verical levels, Block index
    &                                   jg,                          & !< Patch id
    &                                   n_modes,                     & !< number of modes
    &                                   tracer_idx(:)                  !< Tracer indices
  REAL(wp), INTENT(IN) ::               rho(:,:,:),                  & !< Air density
    &                                   tracer(:,:,:,:),             & !< Tracer mixing ratios [kg/kg]
    &                                   dz(:,:,:)                      !< Layer height
  CHARACTER(LEN=*), INTENT(IN) ::       mode_names(:)                  !< Mode names
  !local variables
  LOGICAL ::               available(n_modes)
  INTEGER ::               i, jk, jkp1, jkm1, jc, i_wavel, n_successes !, i550nm
  REAL(wp) ::              ext_value, bsc_value
  TYPE(t_diag_optprops) :: modediag(n_modes)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: att_arr

  n_successes = 0
  DO i = 1,n_modes
    IF (tracer_idx(i) > 0) THEN
      CALL art_find_diagnostic_optprops_for_tracer(jg, mode_names(i), modediag(i), available(i))
      IF (available(i)) THEN
        n_successes = n_successes + 1
      ENDIF
    ELSE
      available(i) = .FALSE.
    ENDIF
  END DO

  IF (n_successes == 0) THEN
    !TODO: error message with specification of the bsc diag
    CALL message("mo_art_aero_opt_props:art_calc_single_backscatter", &
      &          "Could not find all modes corresponding to this backscatter diagnostic.")
    RETURN
  END IF

  DO i_wavel = 1, 3
    IF (ASSOCIATED(ceilometer(i_wavel)%bsc)) THEN
      DO jk=1,nlev
        DO jc = istart, iend
           bsc_value = 0.0_wp
           DO i = 1,n_modes
             IF (available(i)) THEN
               bsc_value = bsc_value + modediag(i)%bsc_satellite(i_wavel) &
                 &                   * tracer(jc,jk,jb,tracer_idx(i))
             END IF
           ENDDO
           !factor 1.e-6_wp converts tracer from ug/kg to g/kg as extinction is given in nm2/g
           ceilometer(i_wavel)%bsc(jc,jk,jb) = bsc_value * rho(jc,jk,jb) * 1.e-6_wp
        ENDDO !jc
      ENDDO !jk
    ENDIF
  ENDDO ! i_wavel

  !calculate attenuated backscatter
  ALLOCATE(att_arr(iend,nlev))

  ! ... for ceilometer
  DO i_wavel = 1, 3
    IF (ASSOCIATED(ceilometer(i_wavel)%bsc) .AND.     &
      & ASSOCIATED(attenuation(i_wavel)%ceil_bsc)) THEN
      DO jk = nlev, 1, -1
        jkp1 = MIN(nlev, jk+1)
        !NEC$ ivdep
        DO jc = istart, iend
          IF (jk == nlev) att_arr(jc,jk) = 0.0_wp
          ext_value = 0.0_wp
          DO i = 1,n_modes
            IF (available(i)) THEN
              ext_value = ext_value + modediag(i)%ext_satellite(i_wavel) &
                &                   * tracer(jc,jk,jb,tracer_idx(i))
            END IF
          ENDDO
          !factor 1.e-6_wp converts tracer from ug/kg to g/kg as extinction is given in nm2/g
          att_arr(jc,jk) = att_arr(jc,jkp1) + &
            &              ext_value * rho(jc,jk,jb) * 1.e-6_wp * dz(jc,jk,jb)
          IF (ceilometer(i_wavel)%bsc(jc,jk,jb) > 0.0_wp) THEN
            attenuation(i_wavel)%ceil_bsc(jc,jk,jb) = &
              &  ceilometer(i_wavel)%bsc(jc,jk,jb) * EXP(-2.0_wp * att_arr(jc,jk))
          ELSE
            attenuation(i_wavel)%ceil_bsc(jc,jk,jb) = 0.0_wp
          END IF
        ENDDO !jc
      ENDDO ! jk
    ENDIF
  ENDDO ! i_wavel

  ! ... for satellite
  DO i_wavel = 1, 3
    IF (ASSOCIATED(ceilometer(i_wavel)%bsc) .AND.     &
      & ASSOCIATED(satellite(i_wavel)%sat_bsc)) THEN
      DO jk = 1, nlev
        jkm1 = MAX(1, jk-1)
        !NEC$ ivdep
        DO jc = istart, iend
          IF (jk == 1) att_arr(jc,jk) = 0.0_wp
          ext_value = 0.0_wp
          DO i = 1,n_modes
            IF (available(i)) THEN
              ext_value = ext_value + modediag(i)%ext_satellite(i_wavel) &
                &                   * tracer(jc,jk,jb,tracer_idx(i))
            END IF
          ENDDO
          !factor 1.e-6_wp converts tracer from ug/kg to g/kg as extinction is given in nm2/g
          att_arr(jc,jk) = att_arr(jc,jkm1) + &
            &              ext_value * rho(jc,jk,jb) * 1.e-6_wp * dz(jc,jk,jb)
          IF (ceilometer(i_wavel)%bsc(jc,jk,jb) == 0.0_wp) THEN
            satellite(i_wavel)%sat_bsc(jc,jk,jb) = 0.0_wp
          ELSE
            satellite(i_wavel)%sat_bsc(jc,jk,jb) = &
              &  ceilometer(i_wavel)%bsc(jc,jk,jb) * EXP(-2.0_wp * att_arr(jc,jk))
          END IF
        ENDDO !jc
      ENDDO !jk
    ENDIF
  ENDDO ! i_wavel

  DEALLOCATE( att_arr )

END SUBROUTINE art_calc_single_backscatter


!
! public subroutines
!

SUBROUTINE art_calc_aod(rho, tracer, dz, istart, iend, nlev, jb, jg, var_med_dia, p_art_data)

!<
! SUBROUTINE art_calc_aod
! This subroutine calculates the aerosol optical depth at 9 different
! Wavelengths used by AERONET
! Extinction coefficient values calculated by P. Gasch, 2015 (Dust)
! and C. Walter (2015) (Ash)
! Part of Module: mo_art_aero_optical_props
! Author: Daniel Rieger, KIT
! Initial Release: 2014-08-04
! Modifications:
! 2014-11-24: Daniel Rieger, KIT
! - Put block loop into interface and adapted values to 550 nm
! 2015-12-15: Philipp Gasch, Carolin Walter, KIT
! - Combined dust/seas AOD routines to a more
!   generic routine for every aerosol
! 2025-01-10: Enrico P. Metzner, KIT
! - restructure the code
! - move hardcoded optical properties into netcdf file
!>

  REAL(wp), INTENT(in)   :: &
    &  rho(:,:,:),          & !< Air density
    &  tracer(:,:,:,:),     & !< Tracer mixing ratios [kg/kg]
    &  dz(:,:,:)              !< Layer height
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of nproma loop
    &  nlev, jb,            & !< Number of verical levels, Block index
    &  jg,                  & !< Patch id
    &  var_med_dia            !< control variable for varying median parametrization
                              !  for dust (1=varying median dia,0=constant median dia)
  TYPE(t_art_data),INTENT(inout) :: &
    &  p_art_data             !< Data container for ART
  ! Local variables
  LOGICAL  ::  success
  INTEGER  ::                            &
    ! tracer container indices for:
    &  iso4_sol_ait, iso4_sol_acc,       & !< soluble sulfate (SO4)
    &  iash_insol_acc, iash_insol_coa,   & !< insoluble ash
    &  iash_mixed_acc, iash_mixed_coa,   & !< mixed ash
    &  iash_giant,                       & !< giant ash
    &  idusta, idustb, idustc,           & !< dust (generic)
    &  idust_insol_acc, idust_insol_coa, & !< insoluble dust
    &  idust_giant,                      & !< giant dust
    &  iseasa, iseasb, iseasc,           & !< seasalt (generic)
    &  inacl_sol_acc, inacl_sol_coa,     & !< soluble sodium chloride (NaCl)
    &  ina_sol_acc, ina_sol_coa,         & !< soluble sodium ions (Na+)
    &  icl_sol_acc, icl_sol_coa,         & !< soluble chloride ions (Cl-)
    &  iasha, iashb, iashc,              & !< volcanic ash (generic)
    &  isoot,                            & !< soot (generic)
    &  isoot_insol_ait, isoot_insol_acc    !< insoluble soot

  iso4_sol_ait   = art_get_tracer_index(p_art_data, 'so4_sol_ait')
  iso4_sol_acc   = art_get_tracer_index(p_art_data, 'so4_sol_acc')
  iash_insol_acc = art_get_tracer_index(p_art_data, 'ash_insol_acc')
  iash_insol_coa = art_get_tracer_index(p_art_data, 'ash_insol_coa')
  iash_mixed_acc = art_get_tracer_index(p_art_data, 'ash_mixed_acc')
  iash_mixed_coa = art_get_tracer_index(p_art_data, 'ash_mixed_coa')
  iash_giant     = art_get_tracer_index(p_art_data, 'ash_giant')

  ! Calculate AOD at specific wavelengths : SO4 - just for 550nm
  IF (iso4_sol_ait > 0 .OR. iso4_sol_acc > 0 ) THEN
    CALL art_calc_aod_single_diag(p_art_data%diag%so4_sol_aeronet,        &
      &                           istart, iend, nlev, jg, jb,             &
      &                           tracer, dz, rho, 1,                     &
      &                           2, (/iso4_sol_ait,iso4_sol_acc/),       &
      &                           (/'sol_ait','sol_acc'/),                &
      &                           success )
  ENDIF

 ! Calculate AOD at specific wavelengths : ash insol - just for 550nm
  IF (iash_insol_acc > 0 .OR. iash_insol_coa > 0 ) THEN
    CALL art_calc_aod_single_diag(p_art_data%diag%ash_insol_aeronet,      &
      &                           istart, iend, nlev, jg, jb,             &
      &                           tracer, dz, rho, 1,                     &
      &                           2, (/iash_insol_acc,iash_insol_coa/),   &
      &                           (/'insol_acc','insol_coa'/),            &
      &                           success )
  ENDIF

 ! Calculate AOD at specific wavelengths : ash mixed - just for 550nm
  IF (iash_mixed_acc > 0 .OR. iash_mixed_coa > 0 ) THEN
    CALL art_calc_aod_single_diag(p_art_data%diag%ash_mixed_aeronet,      &
      &                           istart, iend, nlev, jg, jb,             &
      &                           tracer, dz, rho, 1,                     &
      &                           2, (/iash_mixed_acc,iash_mixed_coa/),   &
      &                           (/'mixed_acc','mixed_coa'/),            &
      &                           success )
  ENDIF

 ! Calculate AOD at specific wavelengths : ash giant - just for 550nm
  IF (iash_giant > 0 ) THEN
    CALL art_calc_aod_single_diag(p_art_data%diag%ash_giant_aeronet,      &
      &                           istart, iend, nlev, jg, jb,             &
      &                           tracer, dz, rho, 1,                     &
      &                           1, (/iash_giant/), (/'giant'/),         &
      &                           success )
  ENDIF

  IF (art_config(jg)%lart_diag_out) THEN

    ! ----------------------------------
    ! --- Calculate AOD at specific wavelengths: dust
    ! ----------------------------------
    IF (art_config(jg)%iart_dust > 0 .AND. var_med_dia == 0) THEN
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_aod',  &
        &  'iart_dust > 0 is not (yet) supported on GPU')
#endif

      idusta          = art_get_tracer_index(p_art_data, 'dusta')
      idustb          = art_get_tracer_index(p_art_data, 'dustb')
      idustc          = art_get_tracer_index(p_art_data, 'dustc')
      idust_insol_acc = 0
      idust_insol_coa = 0
      idust_giant     = 0

      IF (idusta == 0) THEN
        idust_insol_acc = art_get_tracer_index(p_art_data, 'dust_insol_acc')
        idust_insol_coa = art_get_tracer_index(p_art_data, 'dust_insol_coa')
        idust_giant     = art_get_tracer_index(p_art_data, 'dust_giant')
      ENDIF

      IF (idusta > 0 .AND. idustb > 0 .AND. idustc > 0) THEN
        CALL art_calc_aod_single_diag(p_art_data%diag%dust_aeronet,        &
          &                           istart, iend, nlev, jg, jb,          &
          &                           tracer, dz, rho, 9,                  &
          &                           3, (/idusta,idustb,idustc/),         &
          &                           (/'dusta','dustb','dustc'/),         &
          &                           success )
        IF (.NOT. success) THEN
          CALL art_calc_aod_single_diag(p_art_data%diag%dust_aeronet,        &
            &                           istart, iend, nlev, jg, jb,          &
            &                           tracer, dz, rho, 9,                  &
            &                           3, (/idusta,idustb,idustc/),         &
            &                           (/'insol_acc','insol_coa','giant    '/), &
            &                           success )
        END IF
      ELSE IF (idust_insol_acc > 0 .AND. idust_insol_coa > 0 .AND. idust_giant > 0) THEN
        CALL art_calc_aod_single_diag(p_art_data%diag%dust_aeronet,           &
          &                           istart, iend, nlev, jg, jb,             &
          &                           tracer, dz, rho, 9,                     &
          &                           2, (/idust_insol_acc,idust_insol_coa,idust_giant/), &
          &                           (/'insol_acc','insol_coa','giant    '/),            &
          &                           success )
        IF (.NOT. success) THEN
          CALL art_calc_aod_single_diag(p_art_data%diag%dust_aeronet,           &
            &                           istart, iend, nlev, jg, jb,             &
            &                           tracer, dz, rho, 9,                     &
            &                           2, (/idust_insol_acc,idust_insol_coa,idust_giant/), &
            &                           (/'dusta','dustb','dustc'/),            &
            &                           success )
        END IF
      ENDIF

    ENDIF

    ! ----------------------------------
    ! --- Calculate AOD at specific wavelengths: sea salt
    ! ----------------------------------
    IF (art_config(jg)%iart_seasalt > 0) THEN
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_aod',  &
        &  'iart_seasalt > 0 is not (yet) supported on GPU')
#endif

      iseasa        = art_get_tracer_index(p_art_data, 'seasa')
      iseasb        = art_get_tracer_index(p_art_data, 'seasb')
      iseasc        = art_get_tracer_index(p_art_data, 'seasc')
      ina_sol_acc   = 0
      ina_sol_coa   = 0
      icl_sol_acc   = 0
      icl_sol_coa   = 0
      inacl_sol_acc = 0
      inacl_sol_coa = 0

      IF (iseasa == 0) THEN
        ina_sol_acc   = art_get_tracer_index(p_art_data, 'na_sol_acc')
        ina_sol_coa   = art_get_tracer_index(p_art_data, 'na_sol_coa')
        icl_sol_acc   = art_get_tracer_index(p_art_data, 'cl_sol_acc')
        icl_sol_coa   = art_get_tracer_index(p_art_data, 'cl_sol_coa')
        inacl_sol_acc = art_get_tracer_index(p_art_data, 'nacl_sol_acc')
        inacl_sol_coa = art_get_tracer_index(p_art_data, 'nacl_sol_coa')
      ENDIF

      IF (iseasa > 0 .AND. iseasb > 0 .AND. iseasc > 0) THEN
        CALL art_calc_aod_single_diag(p_art_data%diag%seas_aeronet,        &
          &                           istart, iend, nlev, jg, jb,          &
          &                           tracer, dz, rho, 9,                  &
          &                           3, (/iseasa,iseasb,iseasc/),         &
          &                           (/'seasa','seasb','seasc'/),         &
          &                           success )
      ELSE IF (ina_sol_acc > 0 .AND. ina_sol_coa > 0 .AND.  &
        &      icl_sol_acc > 0 .AND. icl_sol_coa > 0) THEN
        CALL art_calc_aod_single_diag(p_art_data%diag%seas_aeronet,        &
          &                           istart, iend, nlev, jg, jb,          &
          &                           tracer, dz, rho, 9,                  &
          &                           4, (/ina_sol_acc,ina_sol_coa,icl_sol_acc,icl_sol_coa/),&
          &                           (/'sol_acc','sol_coa','sol_acc','sol_coa'/),           &
          &                           success )
      ELSE IF (inacl_sol_acc > 0 .AND. inacl_sol_coa > 0 ) THEN
        CALL art_calc_aod_single_diag(p_art_data%diag%seas_aeronet,        &
          &                           istart, iend, nlev, jg, jb,          &
          &                           tracer, dz, rho, 9,                  &
          &                           2, (/inacl_sol_acc,inacl_sol_coa/),  &
          &                           (/'sol_acc','sol_coa'/),             &
          &                           success )
      ENDIF

    ENDIF

    ! ----------------------------------
    ! --- Calculate AOD at specific wavelengths: volcanic ash
    ! ----------------------------------
    IF (art_config(jg)%iart_volcano == 2) THEN
!current implementation of art_calc_aod is not compatible with GPU
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_aod',  &
        &  'iart_volcano > 0 is not (yet) supported on GPU')
#endif

      iasha = art_get_tracer_index(p_art_data, 'asha')
      iashb = art_get_tracer_index(p_art_data, 'ashb')
      iashc = art_get_tracer_index(p_art_data, 'ashc')

      IF (iasha > 0 .AND. iashb > 0 .AND. iashc > 0) THEN
        CALL art_calc_aod_single_diag(p_art_data%diag%volc_aeronet,        &
          &                           istart, iend, nlev, jg, jb,          &
          &                           tracer, dz, rho, 9,                  &
          &                           3, (/iasha,iashb,iashc/),            &
          &                           (/'asha','ashb','ashc'/),            &
          &                           success )
      ENDIF

    ENDIF

    ! ----------------------------------
    ! --- Calculate AOD at specific wavelengths: soot
    ! ----------------------------------
    IF (art_config(jg)%iart_fire > 0) THEN
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_aod',  &
        &  'iart_fire > 0 is not (yet) supported on GPU')
#endif

      isoot           = art_get_tracer_index(p_art_data, 'soot')
      isoot_insol_ait = 0
      isoot_insol_acc = 0

      IF (isoot == 0) THEN
        isoot_insol_ait = art_get_tracer_index(p_art_data, 'soot_insol_ait')
        isoot_insol_acc = art_get_tracer_index(p_art_data, 'soot_insol_acc')
      ENDIF

      IF (isoot > 0) THEN
        CALL art_calc_aod_single_diag(p_art_data%diag%soot_aeronet,        &
          &                           istart, iend, nlev, jg, jb,          &
          &                           tracer, dz, rho, 9,                  &
          &                           1, (/isoot/), (/'soot'/),            &
          &                           success )
      ELSE IF (isoot_insol_ait > 0 .AND. isoot_insol_acc > 0) THEN
        !do it with insoluble soot Aitken mode, or only accumulation mode?
        CALL art_calc_aod_single_diag(p_art_data%diag%soot_aeronet,          &
          &                           istart, iend, nlev, jg, jb,            &
          &                           tracer, dz, rho, 9,                    &
          &                           2, (/isoot_insol_ait,isoot_insol_acc/),&
          &                           (/'insol_ait','insol_acc'/),           &
          &                           success )
      ENDIF

    ENDIF

  ENDIF !lart_diag_out

END SUBROUTINE art_calc_aod

SUBROUTINE art_calc_aodvar(rho, tracer, dz, istart, iend, nlev, jb, jg, p_art_data, &
  &                        cmd, ini_cmd, mode_name, l_init_aod)

!<
! SUBROUTINE art_calc_aodvar
! This subroutine calculates the aerosol optical depth at 9 different
! Wavelengths considering the variable cmd at 550 nm
! Part of Module: mo_art_aero_optical_props
! Author: Ali Hoshyaripour, KIT
! Initial Release: 2018-08-14
! Modifications:
! 2025-01-10: Enrico P. Metzner, KIT
! - restructure the code
! - move hardcoded optical properties into netcdf file
!>

  REAL(wp), INTENT(in)   :: &
    &  rho(:,:,:),          & !< Air density
    &  tracer(:,:,:,:),     & !< Tracer mixing ratios [kg/kg]
    &  dz(:,:,:),           & !< Layer height
    &  cmd(:,:),            &
    &  ini_cmd
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of nproma loop
    &  nlev, jb,            & !< Number of verical levels, Block index
    &  jg                     !< Patch id
  CHARACTER(*), INTENT(in) :: &
    &  mode_name
  LOGICAL, INTENT(in) :: &
    &  l_init_aod
  TYPE(t_art_data),INTENT(inout) :: &
    &  p_art_data             !< Data container for ART
  ! Local variables
  INTEGER ::                        &
    &  jc, jk, i_wavel,             &
    &  idusta, idustb, idustc         !< Tracer container indices

  !CHARACTER(LEN=MAX_CHAR_LENGTH) :: &
  !  &  thisroutine = "mo_art_aero_optical_props:art_calc_aodvar"

  idusta = 0
  idustb = 0
  idustc = 0

  IF (art_config(jg)%iart_dust > 0) THEN
    idusta = art_get_tracer_index(p_art_data, 'dusta')
    idustb = art_get_tracer_index(p_art_data, 'dustb')
    idustc = art_get_tracer_index(p_art_data, 'dustc')
  ENDIF

  IF (art_config(jg)%lart_diag_out) THEN

    ! ----------------------------------
    ! --- Calculate AOD at specific wavelengths: dust
    ! ----------------------------------
    IF (art_config(jg)%iart_dust > 0) THEN
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_aodvar',  &
        &  'iart_dust > 0 is not (yet) supported on GPU')
#endif

      IF (mode_name(:4) == 'dust' .AND. l_init_aod) THEN
        DO i_wavel = 1, 9
          IF (ASSOCIATED(p_art_data%diag%dust_aeronet(i_wavel)%tau)) THEN
            DO jk = 1, nlev
              DO jc = istart, iend
                p_art_data%diag%dust_aeronet(i_wavel)%tau(jc,jk,jb) = 0.0_wp
              ENDDO !jc
            ENDDO !jk
          ENDIF
        ENDDO ! i_wavel
      END IF

      SELECT CASE (mode_name)
      CASE ('dusta')
        CALL art_calc_aodvar_aeronet(p_art_data%diag%dust_aeronet,  &
          &                          istart, iend, nlev, jb, jg, 9, &
          &                          ini_cmd, cmd, rho, dz, tracer, &
          &                          idusta, 'dusta')
      CASE ('dustb')
        CALL art_calc_aodvar_aeronet(p_art_data%diag%dust_aeronet,  &
          &                          istart, iend, nlev, jb, jg, 9, &
          &                          ini_cmd, cmd, rho, dz, tracer, &
          &                          idustb, 'dustb')
      CASE ('dustc')
        CALL art_calc_aodvar_aeronet(p_art_data%diag%dust_aeronet,  &
          &                          istart, iend, nlev, jb, jg, 9, &
          &                          ini_cmd, cmd, rho, dz, tracer, &
          &                          idustc, 'dustc')
      CASE DEFAULT
        ! nothing!
      END SELECT

    ENDIF

  ENDIF !lart_diag_out

END SUBROUTINE art_calc_aodvar

SUBROUTINE art_calc_bsc(rho,tracer, dz, istart, iend, nlev, jb, jg, p_art_data)
!<
! SUBROUTINE art_calc_bsc
! This subroutine calculates the aerosol backscatter and attenuated
! backscattter at 3 different Wavelengths used by lidars
! from ground and from satellite
! Backscatter coefficient values calculated by C. Walter (2015) (Ash)
! Part of Module: mo_art_aero_optical_props
! Author: Carolin Walter, KIT
! Initial Release: 2016-03-14
! Modifications:
! 2025-01-10: Enrico P. Metzner, KIT
! - restructure the code
! - move hardcoded optical properties into netcdf file
!>
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:,:),          & !< Air density
    &  tracer(:,:,:,:),     & !< Tracer mixing ratios [kg/kg]
    &  dz(:,:,:)              !< Layer height
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of nproma loop
    &  nlev, jb, jg           !< Number of verical levels, Block index
  TYPE(t_art_data),INTENT(inout) :: &
    &  p_art_data             !< Data container for ART
  ! Local variables
  INTEGER                :: &
    !v Tracer container indices...
    &  idusta, idustb, idustc,           & !< dust
    &  iseasa, iseasb, iseasc,           & !< seasalt
    &  iasha, iashb, iashc,              & !< volcanic ash
    &  isoot,                            & !< soot
    &  isoot_insol_ait, isoot_insol_acc   !< insoluble soot

!  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
!    &  thisroutine = "mo_art_aero_optical_props:art_calc_bsc"

  idusta = art_get_tracer_index(p_art_data, 'dusta')
  idustb = art_get_tracer_index(p_art_data, 'dustb')
  idustc = art_get_tracer_index(p_art_data, 'dustc')
  
  iseasa = art_get_tracer_index(p_art_data, 'seasa')
  iseasb = art_get_tracer_index(p_art_data, 'seasb')
  iseasc = art_get_tracer_index(p_art_data, 'seasc')

  iasha = art_get_tracer_index(p_art_data, 'asha')
  iashb = art_get_tracer_index(p_art_data, 'ashb')
  iashc = art_get_tracer_index(p_art_data, 'ashc')

  IF (art_config(jg)%iart_fire > 0) THEN
    isoot_insol_ait = 0
    isoot_insol_acc = 0
    isoot = art_get_tracer_index(p_art_data, 'soot')
    IF (isoot == 0) THEN
      isoot_insol_ait = art_get_tracer_index(p_art_data, 'soot_insol_ait')
      isoot_insol_acc = art_get_tracer_index(p_art_data, 'soot_insol_acc')
    ENDIF
  ENDIF

  IF (art_config(jg)%lart_diag_out) THEN

    ! ----------------------------------
    ! --- Calculate backscatter at specific wavelengths dust
    ! ----------------------------------
    IF (idusta > 0 .AND. idustb > 0 .AND. idustc > 0) THEN
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_bsc',  &
        &  'idustX > 0 is not (yet) supported on GPU')
#endif

      CALL art_calc_single_backscatter(   &
        &     p_art_data%diag%dust_ceilo, &
        &     p_art_data%diag%dust_att,   &
        &     p_art_data%diag%dust_sat,   &
        &     istart, iend, nlev, jg, jb, &
        &     tracer, dz, rho,            &
        &     3, (/idusta,idustb,idustc/),&
        &     (/'dusta','dustb','dustc'/) )

    ENDIF

!! MARKER: seasalt needs to be adjusted for aerodyn-treatment !!
!! ============================================================
! ----------------------------------
! --- Calculate backscatter at specific wavelengths sea salt
! ----------------------------------
!    IF (iseasa > 0 .AND. iseasb > 0 .AND. iseasc > 0) THEN
!      bsc_seas_m1 = (/ 0.37411_wp,  0.27673_wp,  0.13381_wp /)
!      bsc_seas_m2 = (/ 0.02118_wp,  0.02521_wp,  0.03019_wp /)
!      bsc_seas_m3 = (/ 0.00423_wp,  0.00453_wp,  0.00469_wp /)
!      ext_seas_m1 = (/ 5.60142_wp,  4.30244_wp,  1.67568_wp /)
!      ext_seas_m2 = (/ 0.39031_wp,  0.40685_wp,  0.43607_wp /)
!      ext_seas_m3 = (/ 0.10244_wp,  0.10235_wp,  0.10624_wp /)
!
!      ALLOCATE(bsc_seas(iend,nlev))
!      DO i_wavel = 1, 3
!        IF (ASSOCIATED(p_art_data%diag%seas_ceilo(i_wavel)%bsc)) THEN
!          DO jk = 1, nlev
!            DO jc = istart, iend
!              bsc_seas(jc,jk) = ( bsc_seas_m1(i_wavel) * p_trac(jc,jk,iseasa) &
!                            &   + bsc_seas_m2(i_wavel) * p_trac(jc,jk,iseasb) &
!                            &   + bsc_seas_m3(i_wavel) * p_trac(jc,jk,iseasc) &
!                            &   ) * rho(jc,jk) * 1.e-6_wp
!              p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) = bsc_seas(jc,jk)
!            ENDDO !jc
!          ENDDO !jk
!        ENDIF
!      ENDDO ! i_wavel
!      DEALLOCATE(bsc_seas)
!
!      !calculate attenuated backscatter
!      ALLOCATE(att_seas(iend,nlev))
!      ALLOCATE(ext_seas(iend,nlev))
!
!      ! ... for ceilometer
!      DO i_wavel = 1, 3
!        IF (ASSOCIATED(p_art_data%diag%seas_ceilo(i_wavel)%bsc) .AND.     &
!          & ASSOCIATED(p_art_data%diag%seas_att(i_wavel)%ceil_bsc)) THEN
!          DO jk = nlev, 1, -1
!            jkp1 = MIN(nlev, jk+1)
!!NEC$ ivdep
!            DO jc = istart, iend
!              IF (jk == nlev) att_seas(jc,jk) = 0.0_wp
!              ext_seas(jc,jk) = ( ext_seas_m1(i_wavel) * p_trac(jc,jk,iseasa) &
!                            &   + ext_seas_m2(i_wavel) * p_trac(jc,jk,iseasb) &
!                            &   + ext_seas_m3(i_wavel) * p_trac(jc,jk,iseasc) &
!                            &   ) * rho(jc,jk) * 1.e-6_wp
!              att_seas(jc,jk) = att_seas(jc,jkp1) + ext_seas(jc,jk) * dz(jc,jk)
!              IF (p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) > 0.0_wp) THEN
!                p_art_data%diag%seas_att(i_wavel)%ceil_bsc(jc,jk,jb) = &
!                  &  p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) &
!                  &  * EXP(-2.0_wp * att_seas(jc,jk))
!              ELSE
!                p_art_data%diag%seas_att(i_wavel)%ceil_bsc(jc,jk,jb) = 0.0_wp
!              END IF
!            ENDDO !jc
!          ENDDO ! jk
!        ENDIF
!      ENDDO ! i_wavel
!
!      ! ... for satellite
!      DO i_wavel = 1, 3
!        IF (ASSOCIATED(p_art_data%diag%seas_ceilo(i_wavel)%bsc) .AND.     &
!          & ASSOCIATED(p_art_data%diag%seas_sat(i_wavel)%sat_bsc)) THEN
!          DO jk = 1, nlev
!            jkm1 = MAX(1, jk-1)
!!NEC$ ivdep
!            DO jc = istart, iend
!              IF (jk == 1) att_seas(jc,jk) = 0.0_wp
!              ext_seas(jc,jk) = ( ext_seas_m1(i_wavel) * p_trac(jc,jk,iseasa) &
!                            &   + ext_seas_m2(i_wavel) * p_trac(jc,jk,iseasb) &
!                            &   + ext_seas_m3(i_wavel) * p_trac(jc,jk,iseasc) &
!                            &   ) * rho(jc,jk) * 1.e-6_wp
!              att_seas(jc,jk) = att_seas(jc,jkm1) + ext_seas(jc,jk) * dz(jc,jk)
!              IF (p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) == 0.0_wp) THEN
!                p_art_data%diag%seas_sat(i_wavel)%sat_bsc(jc,jk,jb) = 0.0_wp
!              ELSE
!                p_art_data%diag%seas_sat(i_wavel)%sat_bsc(jc,jk,jb) = &
!                  &  p_art_data%diag%seas_ceilo(i_wavel)%bsc(jc,jk,jb) &
!                  &  * EXP(-2.0_wp * att_seas(jc,jk))
!              END IF
!            ENDDO !jc
!          ENDDO !jk
!        ENDIF
!      ENDDO ! i_wavel
!
!      DEALLOCATE(ext_seas)
!      DEALLOCATE(att_seas)
!    ENDIF !iart_seasalt >0

! ----------------------------------
! --- Calculate backscatter at specific wavelengths volcanic ash
! ----------------------------------
    IF (iasha > 0 .AND. iashb > 0 .AND. iashc > 0) THEN
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_bsc',  &
        &  'iashX > 0 is not (yet) supported on GPU')
#endif
      
      CALL art_calc_single_backscatter(   &
        &     p_art_data%diag%volc_ceilo, &
        &     p_art_data%diag%volc_att,   &
        &     p_art_data%diag%volc_sat,   &
        &     istart, iend, nlev, jg, jb, &
        &     tracer, dz, rho,            &
        &     3, (/iasha,iashb,iashc/),   &
        &     (/'asha','ashb','ashc'/) )
      
    ENDIF

    ! ----------------------------------
    ! --- Calculate backscatter at specific wavelengths for soot
    ! ----------------------------------
    IF (art_config(jg)%iart_fire > 0) THEN
      ! coated soot with cmd 150 nm
#ifdef _OPENACC
      CALL finish('mo_art_aero_optical_props:art_calc_bsc',  &
        &  'iart_fire > 0 is not (yet) supported on GPU')
#endif

      !Can combine isoot ('soot') with isoot_insol_ait,isoot_insol_acc ('insol_ait','insol_acc')
      !as they do not interfere with each other. Either the first is set, or the second.
      CALL art_calc_single_backscatter(          &
        &     p_art_data%diag%soot_ceilo,        &
        &     p_art_data%diag%soot_att,          &
        &     p_art_data%diag%soot_sat,          &
        &     istart, iend, nlev, jg, jb,        &
        &     tracer, dz, rho,                   &
        &     3, (/isoot,isoot_insol_ait,isoot_insol_acc/),&
        &     (/'soot     ','insol_ait','insol_acc'/) )
      
    ENDIF !iart_fire >0

  ENDIF !lart_diag_out

END SUBROUTINE art_calc_bsc

END MODULE mo_art_aero_optical_props
