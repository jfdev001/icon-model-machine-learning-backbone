!
! mo_art_init_radiation
! This module provides initialization of radiation
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

MODULE mo_art_init_radiation
  ! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish  ! TEMP !!!
  ! ART
  USE mo_art_config,                    ONLY: art_config, IART_PATH_LEN
  USE mo_art_read_opt_props,            ONLY: art_read_prognostic_optprops,   &
                                          &   art_read_diagnostic_optprops,   &
                                          &   art_read_dia_factors,           &
                                          &   art_read_vardia_optprops
  IMPLICIT NONE
    
  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_init_radiation'
  
  PUBLIC :: art_init_radiation
  PUBLIC :: art_init_radiation_diagnostics
  
  CONTAINS
  
  SUBROUTINE art_init_radiation(jg, ncid, lut_optics, use_var_dia, &
    &                           ext_coeff, ssa_coeff, asy_coeff,        &
    &                           ext_param, ssa_param, asy_param,        &
    &                           ext_default, ssa_default, asy_default,  &
    &                           dia_min_factor, dia_max_factor)
  !<
  ! SUBROUTINE art_init_radiation
  ! This subroutine initializes the radiation
  ! Part of Module: mo_art_init_radiation
  ! Author: Daniel Rieger, Philipp Gasch, Carolin Walter, KIT
    ! Initial Release: 2015-12-15
  ! Modification:
  ! 2024-10-18: Enrico P. Metzner, KIT
  !  - moved look-up tables into own xml-file
  !  - reading of polynomial data is only done, if modes*.xml contains attrib cvar_dia="true" in XML-element <lut_optics>
  ! 2025-01-10: Enrico P. Metzner, KIT
  !  - moved look-up tables again, now into a netcdf-file, which is more suitable for a database
  !>
    INTEGER, INTENT(IN)            :: &
      &  jg, ncid
    CHARACTER(LEN=*), INTENT(in)   :: &
      &  lut_optics
    LOGICAL, INTENT(in)            :: &
      &  use_var_dia
    REAL(wp),POINTER,INTENT(inout) :: &
      &  ext_coeff(:),                &
      &  ssa_coeff(:),                &
      &  asy_coeff(:),                &
      &  ext_param(:,:),              &
      &  ssa_param(:,:),              &
      &  asy_param(:,:),              &
      &  ext_default(:,:),            &
      &  ssa_default(:,:),            &
      &  asy_default(:,:),            &
      &  dia_min_factor(:),           &
      &  dia_max_factor(:)
   ! Local Variables
    INTEGER ::                        &
      &  jspec,                       & !< Loop indices
      &  npoly                          !< degree of polynom+1 used for diameter parameterization of coefficients

    jspec=30
    npoly=4

    ALLOCATE(ext_coeff(jspec))
    ALLOCATE(ssa_coeff(jspec))
    ALLOCATE(asy_coeff(jspec))
    
    ALLOCATE(ext_param(jspec,npoly))
    ALLOCATE(ssa_param(jspec,npoly))
    ALLOCATE(asy_param(jspec,npoly))

    ALLOCATE(ext_default(jspec,2))
    ALLOCATE(ssa_default(jspec,2))
    ALLOCATE(asy_default(jspec,2))

    ALLOCATE(dia_min_factor(1))
    ALLOCATE(dia_max_factor(1))

    CALL art_read_prognostic_optprops(ncid, TRIM(lut_optics), ext_coeff, ssa_coeff, asy_coeff)

    CALL art_read_dia_factors(ncid, dia_min_factor, dia_max_factor)

    IF (use_var_dia) THEN
      CALL art_read_vardia_optprops(ncid, TRIM(lut_optics),               &
        &                           ext_param, ssa_param, asy_param,      &
        &                           ext_default, ssa_default, asy_default)
    ENDIF !use_var_dia

    IF (art_config(jg)%iart_nonsph > 0) THEN
#ifdef _OPENACC
      CALL finish('mo_art_init_radiation:art_init_radiation',  &
        &  'iart_nonsph > 0 is not (yet) supported on GPU')
#endif
    ENDIF !iart_nonsph>0

  END SUBROUTINE art_init_radiation




  SUBROUTINE art_init_radiation_diagnostics(ncid, lut_optics,             &
    &                                       ext_aeronet, ext_snglwave,    &
    &                                       ext_satellite, bsc_satellite)
  !<
  ! SUBROUTINE art_init_radiation_diagnostics
  ! This subroutine initializes the optical properties
  ! used for radiation diagnostics
  ! Part of Module: mo_art_init_radiation
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-02-10
  ! Modification:
  !>
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: &
      &  ncid
    CHARACTER(LEN=*), INTENT(IN)   :: &
      &  lut_optics
    REAL(wp),POINTER,INTENT(INOUT) :: &
      &  ext_aeronet(:),              &
      &  ext_snglwave(:),             &
      &  ext_satellite(:),            &
      &  bsc_satellite(:)

    ALLOCATE(ext_aeronet(9))

    ALLOCATE(ext_snglwave(1))

    ALLOCATE(ext_satellite(3))
    ALLOCATE(bsc_satellite(3))

    CALL art_read_diagnostic_optprops(ncid, TRIM(lut_optics), &
      &                               ext_aeronet,            &
      &                               ext_snglwave,           &
      &                               ext_satellite,          &
      &                               bsc_satellite)

  END SUBROUTINE art_init_radiation_diagnostics

END MODULE mo_art_init_radiation
