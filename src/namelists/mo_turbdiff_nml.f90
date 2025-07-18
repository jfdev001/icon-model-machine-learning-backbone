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

! @brief: Namelist for turbulent processes (turbdiff)
!
! Completing the 'turbdiff_config'-setup contained in module 'mo_turbdiff_config' by including
!  input of "NAMELIST/turbdiff_nml/" and loading the domain-specific configuration state.
! This is done through subroutine 'read_turbdiff_namelist' (called in 'read_atmo_namelists').

MODULE mo_turbdiff_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_turbdiff_config !global USE: - all individal quantities being associated to components of 'turbdiff_config(jg)'
                         !            - SUB 'load_turbdiff_config' so as to load these components by updated values
                         !            - config-state vector 'turbdiff_config' (only for "!$ACC UPDATE DEVICE"-directive)
                         !Note(MR): This indicates that 'mo_turbdiff_nml' should be a part of 'mo_turbdiff_config'!
!dom_spec<
! USE mo_turbdiff_config, pat_len_def => pat_len !only, if 'pat_len'-values shall be domain-specific
!dom_spec>

  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_turbdiff_namelist

  !--------------------------------------------------------------------------------------------------

  ! Private auxilary declarations for namelist variables:
  !
  ! - Each of these variables is associated to a component of 'turbdiff_config(jg)'
  !    and may be a vector in case of domain-specific settings.

!dom_spec<
! REAL(wp):: pat_len(0:max_dom) !only, if 'pat_len'-values shall be domain-specific
!dom_spec>
  !--------------------------------------------------------------------------------------------------

  NAMELIST/turbdiff_nml/ &
  ! parameters, switches and selectors from 'mo_turbdiff_config':
    & impl_s, impl_t, tkhmin, tkmmin, tkhmin_strat, tkmmin_strat, &
    & imode_frcsmot, &
     & frcsmot, tkesmot, &
    & rlam_heat, rlam_mom, rat_lam, rat_sea, rat_glac, &
    & imode_charpar, &
     & alpha0, alpha0_max, alpha1, &
    & lconst_z0, &
     & const_z0, tur_len, pat_len, &
    & c_diff, a_stab, a_hshr, &
    & q_crit, &
    & ltkesso, ltkecon, ltkeshs, ltmpcor, lcpfluc, lsflcnd, &
    & ldiff_qi, ldiff_qs, lfreeslip, &
    & imode_tran, imode_turb, icldm_tran, icldm_turb, itype_wcld, itype_sher, &
    & imode_shshear, imode_tkesso, imode_snowsmot, imode_tkemini

  ! Note:
  ! The individual variable names applied in the namelist are taken from 'mo_turbdiff_config'
  !  and have already been initialized there by default values at declaration.
  ! Some of them will be overwritten by namelist-settings below.
  ! If domain-specific values shall be given by the above namelist for any variable, say 'pat_len',
  ! - add ", pat_len_def => pat_len" to "USE mo_turbdiff_config' above
  ! - define a local vector 'pat_len(:)' in the module-header just before the namelist declaration
  ! (see out-commented lines parenthesized by "!dom_spec<" and "!dom-spec>")

CONTAINS

  !-------------------------------------------------------------------------
  !
  !! Read Namelist for turbulent diffusion.
  !!
  !! This subroutine
  !! - reads the Namelist for turbulent diffusion
  !! - sets default values (for domain-specific quantities)
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (almost fully)
  !
  SUBROUTINE read_turbdiff_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_turbdiff_nml: read_turbdiff_nml'

    ! Note:
    ! If domain-specific values shall be given by the above namelist for any variable, say 'pat_len',
    ! - copy default to vector 'pat_len' via the line "pat_len = pat_len_def" in section 1.
    ! - include the line "pat_len_def = pat_len(jg)" just before "CALL load_turbdiff_config(jg)"
    ! (see out-commented lines parenthesized by "!dom_spec<" and "!dom-spec>")

    !------------------------------------------------------------------
    ! 1. default settings of namelist variables are taken from initialization
    !    of 'turbdiff_config' in MODULE 'mo_turbdiff_config'
    !------------------------------------------------------------------

!dom_spec<
!   pat_len = pat_len_def !only, if 'pat_len'-values shall be domain-specific
!dom_spec>

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------

    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('turbdiff_nml')
      READ(funit,NML=turbdiff_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('turbdiff_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, turbdiff_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, turbdiff_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, turbdiff_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml


    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    IF (.NOT. ltkesso) imode_tkesso = 0

    ! The product rlam_heat*rat_sea should be on the order of 10; otherwise, the surface fluxes over oceans
    ! will be unrealistic. Comment the following lines if you want to do it anyway.
    !
    IF (rlam_heat*rat_sea < 2._wp .OR. rlam_heat*rat_sea > 15._wp) THEN
      CALL finish( TRIM(routine), 'The product of rlam_heat and rat_sea should be between 2 and 15')
    ENDIF


    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 1,max_dom

!dom_spec<
!     pat_len_def = pat_len(jg) !only, if 'pat_len'-values shall be domain-specific
!dom_spec>

      CALL load_turbdiff_config(jg) !filling the full configuration-state 'turbdiff_config' (except its part 6.)

    END DO
    !$ACC UPDATE DEVICE(turbdiff_config) ASYNC(1)

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=turbdiff_nml)
      CALL store_and_close_namelist(funit, 'turbdiff_nml')
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=turbdiff_nml)


  END SUBROUTINE read_turbdiff_namelist

END MODULE mo_turbdiff_nml
