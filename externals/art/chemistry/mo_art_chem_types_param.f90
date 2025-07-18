!
! mo_art_chem_types_param
! This module provides data storage structures and constants for parametrised
! chemical tracers (thus extensions of t_chem_meta_param)
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

MODULE mo_art_chem_types_param
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_physical_constants,            ONLY: amd
  USE mo_exception,                     ONLY: finish, message, message_text
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_art_config,                    ONLY: art_config
  USE mo_radiation_config,              ONLY: irad_o3
  USE mo_read_interface,                ONLY: openInputFile, closeFile, read_1D, read_1D_extdim_time


! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN,    &
                                          &   IART_LINOZ_ANA,     &
                                          &   IART_CHEM_NO,       &
                                          &   IART_POLARCHEM,     &
                                          &   IART_LINOZ_LT,      &
                                          &   IART_SIMNOY_SEDI,   &
                                          &   IART_SIMNOY_WMO,    &
                                          &   IART_SIMNOY_PRES,   &
                                          &   IART_SIMNOY_EXTP,   &
                                          &   IART_SIMNOY_UBC,    &
                                          &   IART_SIMNOY_NOUBC
  USE mo_art_OH_chem_types,             ONLY: t_art_OH_chem_meta
  USE mo_art_kinetic_constants,         ONLY: art_determine_bimolecular_kinetic_constant,  &
                                          &   art_determine_termolecular_kinetic_constant, &
                                          &   art_get_CO_des_1d,art_get_CH4_des_1d
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_chem_data,                 ONLY: t_art_chem, &
                                          &   t_art_chem_indices
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_read_simnoy,               ONLY: n2onoy_read
  USE mo_art_chem_types,                ONLY: t_chem_meta_param,  &
                                          &   t_prod_list
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string
  
 

  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_chem_meta_lt
  PUBLIC :: t_chem_meta_OH
  PUBLIC :: t_chem_meta_linoz
  PUBLIC :: t_chem_meta_linozv3
  PUBLIC :: t_chem_meta_simnoy
  PUBLIC :: t_chem_meta_cold
  PUBLIC :: OH_get_destruct, linoz_fill_init, linoz_integrate, linozv3_fill_init, linozv3_integrate
  PUBLIC :: simnoy_fill_init


  !##################################################################
  ! type of tracers depleted with a constant or parametrised lifetime
  !##################################################################

  TYPE, extends(t_chem_meta_param) :: t_chem_meta_lt
    REAL(wp)              ::   &
      &  des            !< constant global destruction rate (inverse of lifetime) (s-1)

    CONTAINS 
      PROCEDURE :: get_destruct => lt_get_destruct
      PROCEDURE :: init => lt_init_arrays
      PROCEDURE :: integrate  => art_implicit_methods_3d
      PROCEDURE :: art_implicit_methods_1d
      PROCEDURE :: art_implicit_methods_3d

  END TYPE t_chem_meta_lt

  !###########################################################################
  ! type of cold tracer (it has special features and interacts with Simnoy and
  !                      linoz tracers, this is why it has its own type)
  !###########################################################################

  TYPE, extends(t_chem_meta_param) :: t_chem_meta_cold

    REAL(wp)              ::   &
      &  des                !< constant global destruction rate (inverse of lifetime) (s-1)
    REAL(wp)              ::   &
      &  des_coldsed        !< sedimentation destruction rate of cold tracer (s-1)
    INTEGER               ::   &
      &  polarchem          !< switch which sort of polar chemistry should be applied
    REAL(wp), ALLOCATABLE ::   &
      &  p_cold_sed(:,:,:)  !< sedimented part of the cold tracer (used for simnoy) (#/m3)


    CONTAINS 
      PROCEDURE :: get_destruct => cold_get_destruct
      PROCEDURE :: init => cold_init_arrays
      PROCEDURE :: integr_and_prescribe  => cold_integrate_and_prescribe

  END TYPE t_chem_meta_cold

  !##################################################################
  ! type of tracers depleted with with interactively calculated  OH
  !##################################################################

  TYPE, extends(t_chem_meta_param) :: t_chem_meta_OH
    REAL(wp) :: &
      &  des_1d               !< 1-D destruction rate
    REAL(wp), ALLOCATABLE ::   &
      &  des_star(:,:,:)      !< star value after predictor step (s-1)
    REAL(wp), ALLOCATABLE ::   &
      &  prod_star(:,:,:)     !< star production after predictor step ([tracer_star]/s)
    REAL(wp), ALLOCATABLE ::   &
      &  tracer_star(:,:,:)   !< tracer mass mixing ratio after predictor step (kg/kg)

    CONTAINS
      PROCEDURE :: get_tracer_name => OH_get_tracer_name
      PROCEDURE :: get_destruct => OH_get_destruct
      PROCEDURE :: get_prod_star => OH_get_prod_star
      PROCEDURE :: init => OH_init_arrays

  END TYPE t_chem_meta_OH

  !##################################################################
  ! type of Linoz ozone tracers 
  !##################################################################

  TYPE, extends(t_chem_meta_param) :: t_chem_meta_linoz
      
    REAL(wp), ALLOCATABLE ::  &
      &  tend(:,:,:)    !< tendency of ozone that is finally added to the tracer (kg/kg)
    INTEGER               ::  &
      &  O3_paramet     !< switch which solver for O3 should be used (analytic or not)
    INTEGER               ::  &
      &  O3_feed        !< radiation feedback or not?
    REAL(wp), POINTER ::   &
      &  cold_tracer(:,:,:) => NULL()  !< pointer to the cold tracer  concentration
    INTEGER               ::  &
      &  polarchem      !< switch which polar chemistry should be applied to the tracer
    REAL(wp)              ::  &
      &  o3lt_het       !< heterogeneous lifetime of O3, only if cold tracer
                        !  is not available (s)
    REAL(wp)              ::  &
      &  Thet           !< temperature below which polarchem applies (K)
    REAL(wp)              ::  &
      &  o3_lbc        !< lower boundary ozone volume mixing ratio
    REAL(wp)              ::  &
      &  lt_lbc         !< ozone lifetime for lower boundary layer
    REAL(wp), ALLOCATABLE ::  &
      &  column(:,:,:)  !< column value of ozone (DU)
    
    CONTAINS
      PROCEDURE :: init => linoz_init_arrays
      PROCEDURE :: fill_init => linoz_fill_init
      PROCEDURE :: set_tracer_linoz
      PROCEDURE :: integrate => linoz_integrate

  END TYPE t_chem_meta_linoz

  !##################################################################
  ! type of Linozv3 ozone tracers
  !##################################################################

  TYPE, extends(t_chem_meta_param) :: t_chem_meta_linozv3

    REAL(wp), ALLOCATABLE ::  &
      &  tend(:,:,:)    !< tendency of ozone that is finally added to the tracer (kg/kg)
    REAL(wp), ALLOCATABLE ::   &
      &  tend_NOy(:,:,:) !< tendency of NOy that is finally added to the tracer (kg/kg)
    INTEGER               ::  &
      &  O3_paramet     !< switch which solver for O3 should be used (analytic or not)
    INTEGER               ::  &
      &  O3_feed        !< radiation feedback or not?
    INTEGER               ::  &
      &  sol_ssi        !< switch for solar ssi: 0 = variabel, 1 = solar min, 2 = solar max
    REAL(wp)              ::  &
      &  des            !< tropospheric destruction rate from lifetime (s-1)
    REAL(wp), POINTER ::   &
      &  cold_tracer(:,:,:) => NULL()  !< pointer to the cold tracer  concentration
    REAL(wp), POINTER ::   &
      &  NOy_tracer(:,:,:) => NULL()  !< pointer to the noy tracer  concentration
    INTEGER               ::  &
      &  polarchem      !< switch which polar chemistry should be applied to the tracer
    REAL(wp)              ::  &
      &  o3lt_het       !< heterogeneous lifetime of O3, only if cold tracer
                        !  is not available (s)
    REAL(wp)              ::  &
      &  Thet           !< temperature below which polarchem applies (K)
    REAL(wp), ALLOCATABLE ::  &
      &  column(:,:,:)  !< column value of ozone (DU)
    REAL(wp), ALLOCATABLE ::  &
      &  f107(:,:)      !<  time, f107 from external file

    CONTAINS
      PROCEDURE :: init => linozv3_init_arrays
      PROCEDURE :: fill_init => linozv3_fill_init
      PROCEDURE :: set_tracer_linozv3
      PROCEDURE :: integrate => linozv3_integrate
      PROCEDURE :: convert_mmr => convert_vmr_to_mmrv3
      PROCEDURE :: convert_vmr => convert_mmr_to_vmrv3

  END TYPE t_chem_meta_linozv3


  !##################################################################
  ! type of Simnoy tracers 
  !##################################################################

  TYPE, extends(t_chem_meta_param) :: t_chem_meta_simnoy
      
    REAL(wp), ALLOCATABLE ::   &
      &  tend(:,:,:)     !< tendency of NOy that is finally added to the tracer (kg/kg)
    REAL(wp), ALLOCATABLE ::   &
      &  tend_n2o(:,:,:) !< tendency of N2O that is finally added to the tracer (kg/kg)
    REAL(wp)              ::  &
      &  des             !< constant destruction rate from lifetime that applies
                         !  dependent on cnoyn2o_tropo (s-1)
    REAL(wp)              ::  &
      &  des_n2o         !< destruction rate of N2O
    REAL(wp)              ::  &
      &  des_noysed      !< sedimentation destruction rate of NOy
    INTEGER               ::  &
      &  cnoyn2o_tropo   !< switch how to treat the troposphere
    INTEGER               ::  &
      &  ubc             !< switch how to treat the upper boundary
    REAL(wp), POINTER ::                &
      &  cold_tracer(:,:,:) => NULL(),  & !< pointers to the concentrations of
      &  n2o_tracer(:,:,:) => NULL(),   & !  cold, n2o and sedimented part of the cold tracer
      &  p_cold_sed(:,:,:) => NULL()
    INTEGER               ::   &
      &  polarchem       !< switch which polar chemistry should be applied to the tracer
    REAL(wp), ALLOCATABLE ::   &
      &  n2onoy_tab1(:,:),     & !< tabular values of the destruction rates for Simnoy
      &  n2onoy_tab2(:,:),     &
      &  n2onoy_tab3(:,:),     &
      &  n2onoy_tab4(:,:)


    CONTAINS
      PROCEDURE :: init => simnoy_init_arrays
      PROCEDURE :: fill_init => simnoy_fill_init
      PROCEDURE :: integrate => simnoy_integrate
      PROCEDURE :: set_tracer_simnoy

  END TYPE t_chem_meta_simnoy



CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for t_chem_meta_lt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE lt_init_arrays(this,nproma,nlev,nblks)
!<
! SUBROUTINE lt_init_arrays
! Initialise arrays for tracer destroyed by constant lifetime
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_lt), INTENT(inout) :: &
    &  this                     !< container with fields
  INTEGER, INTENT(in) :: &
    &  nproma,nlev,nblks    !< dimensions

  CALL this%init_param(nproma,nlev,nblks)

  IF (.NOT. ALLOCATED(this%des_3d)) ALLOCATE(this%des_3d(nproma,nlev,nblks))

  this%des_3d(:,:,:) = 0.0_wp
END SUBROUTINE lt_init_arrays



SUBROUTINE lt_get_destruct(this)
!<
! SUBROUTINE lt_get_destruct
! Calculate the destruction rate
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter and Michael Weimer, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_lt),INTENT(inout) :: &
    &  this         !< Container with fields
  REAL(wp) ::   &
    &  lt_tr        !< lifetime of the tracer (s)
  INTEGER :: &
    &  ierror
  CHARACTER(LEN = IART_VARNAMELEN) :: &
    &  tracer_name
  CHARACTER(:),ALLOCATABLE         :: &
    &  c_tmp

  CALL this%opt_meta%get('lifetime',lt_tr, ierror)
  IF (ierror == SUCCESS) THEN
    this%des = 1.0_wp/lt_tr
  ELSE
    CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
    WRITE(tracer_name,'(A)') c_tmp
    CALL finish('mo_art_chem_types:lt_get_destruct',      &
            &   'lifetime for '//TRIM(tracer_name)//' missing.')
  END IF

  this%des = 1.0_wp/lt_tr
  this%des_3d = this%des
  this%prod = 0._wp
END SUBROUTINE lt_get_destruct


SUBROUTINE art_implicit_methods_3d(this, p_dtime)
!<
! SUBROUTINE art_implicit_methods
! This subroutine decides via the value of dest which formular to use:
! - analytical formular is exact but instable for low dest values (lower than about 1e-10 1/s)
! - simplified formular is always stable but not exact
! Based on: dissertation Roland Ruhnke (Year?)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: around 2015-09-15                
! Modifications:
!>
  CLASS(t_chem_meta_lt),INTENT(inout) :: &
    &  this        !< Container with fields
  REAL(wp), INTENT(IN)              ::  &
    &  p_dtime     !< time step (s)
    

  IF (ANY(this%des_3d < 1.e-10_wp) )THEN
    ! simplified formular:
    this%tracer   = (this%prod    * p_dtime + this%tracer )/ &
                  & (1._wp + this%des_3d * p_dtime)
  ELSE
    ! analytical formular

    this%tracer   = this%prod / this%des_3d                            &
                       &   + ( this%tracer - this%prod / this%des_3d ) &
                       &   * exp( -1.0_wp * this%des_3d * p_dtime )
  END IF
END SUBROUTINE art_implicit_methods_3d

SUBROUTINE art_implicit_methods_1d(this, p_dtime)
!<
! SUBROUTINE art_implicit_methods
! This subroutine decides via the value of dest which formular to use:
! - analytical formular is exact but instable for low dest values (lower than about 1e-10 1/s)
! - simplified formular is always stable but not exact
! Based on: dissertation Roland Ruhnke (Year?)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: around 2015-09-15                
! Modifications:
!>
  IMPLICIT NONE
  CLASS(t_chem_meta_lt),INTENT(inout) :: &
     &  this         !< Container with fields
   REAL(wp), INTENT(IN)              ::  &
     &  p_dtime      !< time step (s)
   
  IF (this%des > 1.e-10_wp) THEN
    this%tracer   = (this%prod    * p_dtime + this%tracer ) / &
         &          (1._wp + this%des * p_dtime)

  ELSE
    this%tracer   = this%prod / this%des                        &
        &           + ( this%tracer - this%prod / this%des )    &
        &           * exp( -1.0_wp * this%des * p_dtime )
  END IF
END SUBROUTINE art_implicit_methods_1d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for t_chem_meta_cold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE cold_init_arrays(this,nproma,nlev,nblks)
!<
! SUBROUTINE cold_init_arrays
! Initialise arrays for cold tracer
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_cold), INTENT(inout) :: &
    &  this                  !< container with fields
  INTEGER, INTENT(in) :: &
    &  nproma,nlev,nblks  !< dimensions

  CHARACTER(LEN = IART_VARNAMELEN) :: &
    &  polarchem
  REAL(wp) :: &
    &  lt_tr, lt_coldsed
  INTEGER :: &
    &  ierror
  CHARACTER(:),ALLOCATABLE  :: &
    &  c_tmp

  CALL this%init_param(nproma,nlev,nblks)

  IF (.NOT. ALLOCATED(this%des_3d)) ALLOCATE(this%des_3d(nproma,nlev,nblks))
  IF (.NOT. ALLOCATED(this%p_cold_sed)) THEN
     ALLOCATE(this%p_cold_sed(nproma,nlev,nblks))
     this%p_cold_sed = 0._wp
  END IF

  this%des_3d(:,:,:) = 0.0_wp
  

  CALL this%opt_meta%get('lifetime',lt_tr,ierror)

  IF (ierror == SUCCESS) THEN
    this%des = 1.0_wp/lt_tr
  ELSE
    CALL finish('mo_art_chem_types:cold_init_arrays',      &
            &   'lifetime for TR_cold missing.')
  END IF

  CALL this%opt_meta%get('lt_sed', lt_coldsed,ierror)
  IF (ierror /= SUCCESS) lt_coldsed = lt_tr

  this%des_coldsed = 1._wp / lt_coldsed

  CALL key_value_storage_as_string(this%opt_meta,'polarchem', c_tmp, ierror)
  IF (ierror /= SUCCESS) THEN
    polarchem = 'on'
  ELSE
    WRITE(polarchem,'(A)') c_tmp
  END IF

  SELECT CASE (TRIM(polarchem))
    CASE ('on')
      this%polarchem = IART_POLARCHEM
    CASE('off')
      this%polarchem = IART_CHEM_NO
    CASE DEFAULT
      CALL finish('mo_art_chem_types', &
        &   'polarchem for Cold tracer must be one of ''on'' or ''off''.')
  END SELECT
END SUBROUTINE cold_init_arrays


SUBROUTINE cold_get_destruct(this, jg)
!<
! SUBROUTINE cold_init_arrays
! Set destruction rate for cold tracer
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_cold),INTENT(inout) :: &
    &  this     !< Container with field
  INTEGER, INTENT(in) :: &
    &  jg       !< patch id

  ! local
  INTEGER :: &
    &  jb, jc, jk, i_startidx, i_endidx
  TYPE(t_art_atmo), POINTER ::  &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jk = 1, art_atmo%nlev
      DO jc = i_startidx, i_endidx
        IF (art_atmo%temp(jc,jk,jb) <= 195._wp .AND. art_atmo%z_mc(jc,jk,jb)/1000._wp <= 50._wp) THEN
          this%des_3d(jc,jk,jb) = 0.0_wp
        ELSE
          this%des_3d(jc,jk,jb) = this%des
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE cold_get_destruct


SUBROUTINE cold_integrate_and_prescribe(this, jg, p_dtime, vmr2Nconc)
!<
! SUBROUTINE cold_implicit_methods
! This subroutine decides via the value of dest which formular to use:
! - analytical formular is exact but instable for low dest values (lower than about 1e-10 1/s)
! - simplified formular is always stable but not exact
! In addition it calculates the polar chemistry for the cold tracer if it should
! be used.
! Based on: dissertation Roland Ruhnke (Year?)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE
  CLASS(t_chem_meta_cold),INTENT(inout) :: &
    &  this              !< Container with fields
  INTEGER, INTENT(in)                 :: &
    &  jg                !< patch id
  REAL(wp), INTENT(IN)              ::  &
    &  p_dtime           !< time step (s)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)  !< volume mixing ratio (mol/mol) to number concentration (# / cm3)

  ! local variables
  INTEGER :: &
    &  i_startidx, i_endidx, jc, jk, jb
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  REAL(wp) :: &
    &  n_cold, n_cold_new


  art_atmo => p_art_data(jg)%atmo

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jk = 1,art_atmo%nlev
      DO jc = i_startidx, i_endidx
        IF (this%des_3d(jc,jk,jb) > 0.0_wp) THEN

          IF (this%des_3d(jc,jk,jb) < 1.e-10_wp) THEN
            ! simplified formular:
            this%tracer(jc,jk,jb)   = (this%prod(jc,jk,jb)    * p_dtime + this%tracer(jc,jk,jb) )/ &
                                      & (1._wp + this%des_3d(jc,jk,jb) * p_dtime)
          ELSE
            ! analytical formular
            this%tracer(jc,jk,jb)   = this%prod(jc,jk,jb) / this%des_3d(jc,jk,jb)   &
              &                    + ( this%tracer(jc,jk,jb) - this%prod(jc,jk,jb)  &
              &                                        / this%des_3d(jc,jk,jb) )    &
              &                    * exp( -1.0_wp * this%des_3d(jc,jk,jb) * p_dtime )
        
          END IF
        ELSE
          this%tracer(jc,jk,jb) = 1._wp * vmr2Nconc(jc,jk,jb)
        END IF
      END DO
    END DO
  END DO

  IF (this%polarchem == IART_POLARCHEM) THEN
    DO jb = art_atmo%i_startblk, art_atmo%i_endblk
      CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
  
      DO jk = 1,art_atmo%nlev
        DO jc = i_startidx, i_endidx
          ! Sedimentation of Cold tracer
          n_cold = this%tracer(jc,jk,jb) * 1.e6_wp
  
          IF (jk >= 2) THEN
            n_cold = n_cold + this%p_cold_sed(jc,jk-1,jb)
          END IF
  
          IF (jk /= art_atmo%nlev) THEN
            n_cold_new  = n_cold * exp( -1._wp * p_dtime * this%des_coldsed)
  
            this%p_cold_sed(jc,jk,jb) = n_cold - n_cold_new
          ELSE
            this%p_cold_sed(jc,jk,jb) = 0._wp
          END IF
  
          this%tracer(jc,jk,jb) = n_cold_new / 1.e6_wp
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE cold_integrate_and_prescribe


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Procedures related to OH chemistry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE OH_get_tracer_name(this,tracer_name)
!<
! SUBROUTINE OH_get_tracer_name
! Get the name of the tracer using the OH chemistry
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_OH), INTENT(in) ::  &
    &  this          !< container with fields
  CHARACTER(:), ALLOCATABLE, INTENT(inout) ::  &
    &  tracer_name   !< name of the tracer

  CALL key_value_storage_as_string(this%opt_meta,'name',tracer_name)

END SUBROUTINE OH_get_tracer_name


SUBROUTINE OH_init_arrays(this,nproma,nlev,nblks)
!<
! SUBROUTINE OH_init_arrays
! Initialise the structures for tracer using the OH chemistry
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_OH), INTENT(inout) :: &
    &  this                   !< container with fields
  INTEGER, INTENT(in) :: &
    &  nproma,nlev,nblks   !< dimensionsS

  ! local variables
  INTEGER :: &
    &  ierror              !< check if element exists
  CHARACTER(:), ALLOCATABLE :: &
    &  tracer_name         !< tracer name
  REAL(wp) :: &
    &  const_lt            !< constant lifetime at altitude higher than CH4 1 ppm (s)

  CALL this%init_param(nproma,nlev,nblks)

  ALLOCATE(this%des_3d(nproma,nlev,nblks))
  ALLOCATE(this%des_star(nproma,nlev,nblks))
  ALLOCATE(this%prod_star(nproma,nlev,nblks))
  ALLOCATE(this%tracer_star(nproma,nlev,nblks))

  this%des_3d = 0._wp
  this%des_star = 0._wp
  this%tracer_star = 0._wp
  this%prod_star = 0._wp

  CALL this%opt_meta%get('lifetime', const_lt, ierror)
  IF (ierror == SUCCESS) THEN
    this%des_1d = 1._wp / const_lt
  ELSE
    CALL this%get_tracer_name(tracer_name)
    CALL finish('mo_art_chem_types:OH_init_arrays',        &
           &    'Lifetime missing for '//TRIM(tracer_name)//'.')
  END IF
END SUBROUTINE OH_init_arrays


SUBROUTINE OH_get_destruct(this, jg, OH_chem_meta, star_vals)
!<
! SUBROUTINE OH_get_destruct
! This subroutine calculates the kinetic constants
! For the user: if new tracer should be added to this mechanism, please include herein 
!               its reaction rate with OH
! Based on: 
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-16
! Modifications: Khompat Satitkovitchai, 2025-01-10
! - added the new tracers (CFC-10, CFC-11 and CFC-12) with O1D reactions 
!>
  IMPLICIT NONE
  CLASS(t_chem_meta_OH), INTENT(inout), TARGET :: &
    &  this              !< container with fields
  INTEGER, INTENT(in)   :: &
    &  jg                !< patch id
  TYPE(t_art_OH_chem_meta), INTENT(inout) ::   &
    &  OH_chem_meta      !< meta data for OH chemistry
  LOGICAL, INTENT(in) ::  &
    &  star_vals         !< flag if star values should be used or not

  !local variables
  INTEGER :: jc,jk,jb, i_startidx, i_endidx   !< loop indices
  REAL(wp), POINTER :: &
    &  des(:,:,:)        !< destruction rate (pointer to star values or not)  (s-1)
  REAL(wp) ::        &
    &  N2_O2_Nconc,  &   !< number concentration of N2 + O2 [# / cm3]
    &  k_tr_OH,      &   !< kinetic constants of this tracer with OH (cm3 / # /s)
    &  k_C5H8_O3,    &
    &  k_tr_O1D,     &   !< kinetic constants of this tracer with O1D (cm3 / # /s) 
    &  A_factor,     &   !< pre-exponential factor 
    &  Ea_value          !< activation energy
  CHARACTER(:), ALLOCATABLE :: &
    &  tracer_name       !< name of the tracer

  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo          !< Pointer to ART atmo fields
  TYPE(t_art_chem),POINTER    :: &
    &  art_chem          !< Pointer to ART chem fields

  CALL this%get_tracer_name(tracer_name)

  art_chem => p_art_data(jg)%chem
  art_atmo => p_art_data(jg)%atmo
  

  IF (star_vals) THEN
    des => this%des_star
  ELSE
    des => this%des_3d
  END IF

  SELECT CASE (TRIM(tracer_name))
    !
    ! CO: reaction rate with OH
    !
    CASE ('TRCO')
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
        DO jc = i_startidx,i_endidx
          DO jk = 1,OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            CALL art_get_CO_des_1d(des(jc,jk,jb), art_atmo%pres(jc,jk,jb))
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            des(jc,jk,jb)  = MAX(OH_chem_meta%k_CO_OH(jc,jk,jb) * OH_chem_meta%OH_Nconc(jc,jk,jb), &
                            &  1.e-30_wp)
          END DO
        END DO
      END DO
    !
    ! CH4: reaction rate with OH
    !
    CASE ('TRCH4') 
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
        DO jc = i_startidx,i_endidx
          DO jk = 1,OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            CALL art_get_CH4_des_1d(des(jc,jk,jb), art_atmo%pres(jc,jk,jb))
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            des(jc,jk,jb)  = MAX(OH_chem_meta%k_CH4_OH(jc,jk,jb) * OH_chem_meta%OH_Nconc(jc,jk,jb), &
                              &  1.e-30_wp)
          END DO
        END DO
      END DO
    !
    ! CH3COCH3: reaction rate with OH
    !
    CASE ('TRCH3COCH3') 
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
        DO jc = i_startidx,i_endidx
          DO jk = 1,OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            des(jc,jk,jb) = this%des_1d
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            CALL art_determine_bimolecular_kinetic_constant        &
              & (k_tr_OH,art_atmo%temp(jc,jk,jb),                  &
              &   3.82e-11_wp,2000._wp,1.33e-13_wp)

            des(jc,jk,jb) = MAX( k_tr_OH * OH_chem_meta%OH_Nconc(jc,jk,jb) &
                               &   +  art_chem%photo%rate(jc,jk,jb,68)     & 
                               &   +  art_chem%photo%rate(jc,jk,jb,69),    &
                               &     1.e-30_wp)
          END DO
        END DO
      END DO
    !
    ! C5H8: reaction rate with OH
    !
    CASE ('TRC5H8')
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
        DO jc = i_startidx,i_endidx
          DO jk = 1,OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            des(jc,jk,jb) = this%des_1d
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            CALL art_determine_bimolecular_kinetic_constant         &
              & (k_tr_OH,art_atmo%temp(jc,jk,jb),3.1e-11_wp,-350._wp)

            CALL art_determine_bimolecular_kinetic_constant  &
              & (k_C5H8_O3,art_atmo%temp(jc,jk,jb),1e-14_wp,1970._wp)

            des(jc,jk,jb)  = MAX(k_tr_OH * OH_chem_meta%OH_Nconc(jc,jk,jb)              &
                             & + k_C5H8_O3 * OH_chem_meta%ozone_Nconc(jc,jk,jb),        &
                             & 1.e-30_wp)
          END DO
        END DO
      END DO
    !
    ! SO2: reaction rate with OH
    !
    CASE ('TRSO2') 
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
        DO jc = i_startidx,i_endidx
          DO jk = 1,OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            des(jc,jk,jb) = this%des_1d
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            N2_O2_Nconc = 0.9903 * art_chem%vmr2Nconc(jc,jk,jb)

            ! exponent "m" is empty in reaction SO2 + OH + M--> HSO3 + M and therefore set to 0  
            ! (same as in MECCA)
            CALL art_determine_termolecular_kinetic_constant &
                  &  (k_tr_OH,art_atmo%temp(jc,jk,jb),  &
                  &   3.3e-31_wp,4.3_wp,1.6e-12_wp,0.0_wp,N2_O2_Nconc,.FALSE.)

            des(jc,jk,jb)  = MAX(k_tr_OH * OH_chem_meta%OH_Nconc(jc,jk,jb), &
                               &     1.e-30_wp)
          END DO
        END DO
      END DO
    !
    ! NO2: reaction rate with OH
    !
    CASE ('TRNO2') 
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

        DO jc = i_startidx,i_endidx
          DO jk = 1,OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            des(jc,jk,jb) = this%des_1d
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            N2_O2_Nconc = 0.9903 * art_chem%vmr2Nconc(jc,jk,jb)
            CALL art_determine_termolecular_kinetic_constant &
                  &  (k_tr_OH,art_atmo%temp(jc,jk,jb), &
                  &   1.8e-30_wp,3.0_wp,2.8e-11_wp,0.0_wp,N2_O2_Nconc,.FALSE.)

            des(jc,jk,jb)  = MAX(k_tr_OH * OH_chem_meta%OH_Nconc(jc,jk,jb), &
                               &     1.e-30_wp)
          END DO
        END DO
      END DO
    !
    ! C2H6,C3H8,OCS,DMS,NH3: reaction rate with OH
    !
    CASE ('TRC2H6', 'TRC3H8', 'TROCS', 'TRDMS', 'TRNH3')
      SELECT CASE (TRIM(tracer_name))
        CASE ('TRC2H6')
          A_factor = 7.66e-12_wp
          Ea_value = 1020._wp
        CASE ('TRC3H8')
          A_factor = 7.6e-12_wp
          Ea_value = 585._wp
        CASE ('TROCS')
          A_factor = 7.2e-14_wp
          Ea_value = 1070._wp
        CASE ('TRDMS')
          A_factor = 1.2e-11_wp
          Ea_value = 280._wp
        CASE ('TRNH3')
          A_factor = 1.7e-12_wp
          Ea_value = 710._wp
      END SELECT

      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

        DO jc = i_startidx,i_endidx
          DO jk = 1,OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            des(jc,jk,jb) = this%des_1d
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            CALL art_determine_bimolecular_kinetic_constant &
                  &  (k_tr_OH, art_atmo%temp(jc,jk,jb), A_factor, Ea_value)

            des(jc,jk,jb) = MAX(k_tr_OH * OH_chem_meta%OH_Nconc(jc,jk,jb), &
                             &     1.e-30_wp)
          END DO
        END DO
      END DO
    !
    ! CFCs: reaction rate with O1D
    ! based on: JPL Publication No. 19-5
    !
    CASE ('TRCCl4', 'TRCFCl3', 'TRCF2Cl2')
      SELECT CASE (TRIM(tracer_name))
        CASE ('TRCCl4')
          A_factor = 3.3e-10_wp
          Ea_value = 0.0_wp
        CASE ('TRCFCl3')
          A_factor = 2.3e-10_wp
          Ea_value = 0.0_wp
        CASE ('TRCF2Cl2')
          A_factor = 1.4e-10_wp
          Ea_value = -25._wp
      END SELECT
 
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

        DO jc = i_startidx, i_endidx
          DO jk = 1, OH_chem_meta%level_CH4_gt_1ppm(jc,jb) - 1
            des(jc,jk,jb) = this%des_1d
          END DO

          DO jk = OH_chem_meta%level_CH4_gt_1ppm(jc,jb), art_atmo%nlev
            CALL art_determine_bimolecular_kinetic_constant &
                 & (k_tr_O1D, art_atmo%temp(jc,jk,jb), A_factor, Ea_value)

            des(jc,jk,jb) = MAX(k_tr_O1D * OH_chem_meta%O1D_Nconc(jc,jk,jb), 1.e-30_wp)
          END DO
        END DO
      END DO
    !
    CASE DEFAULT
      WRITE(message_text,*) 'no reaction rate of tracer ',  &
         &  TRIM(tracer_name), ' with OH found. Please include it.'
      CALL finish('mo_art_chem_types:OH_get_destruct',TRIM(message_text))
  END SELECT

END SUBROUTINE OH_get_destruct


SUBROUTINE OH_get_prod_star(this)
!<
! SUBROUTINE OH_get_prod_star
! This subroutine calculates the star values for the production
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_OH), INTENT(inout) :: &
    &  this
  ! local variables
  TYPE(t_prod_list), POINTER :: &
    &  current_prod

  this%prod_star = 0._wp

  IF (ASSOCIATED(this%first_prod)) THEN
    current_prod => this%first_prod

    DO WHILE(ASSOCIATED(current_prod))
      SELECT TYPE(edu => current_prod%educt)
        TYPE IS (t_chem_meta_OH)
          this%prod_star = this%prod_star + current_prod%factor * edu%des_star * edu%tracer_star
        CLASS DEFAULT
          this%prod_star = this%prod_star + edu%prod
      END SELECT

      current_prod => current_prod%next_prod
    END DO
  END IF
END SUBROUTINE OH_get_prod_star

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Linoz tracer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE linoz_init_arrays(this,nproma,nlev,nblks)
!<
! SUBROUTINE linoz_init_arrays
! Initialise the strcutures for linoz tracer
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter, KIT
! Initial Release: aronund 2018-10
! Modifications:
!>
  CLASS(t_chem_meta_linoz), INTENT(inout) :: &
    &  this               !< container with fields
  INTEGER, INTENT(in) :: &
    &  nproma,nlev,nblks  !< dimensions

  CALL this%init_param(nproma,nlev,nblks)

  IF(.NOT. ALLOCATED(this%tend)) THEN
    ALLOCATE(this%tend(nproma,nlev,nblks))
    this%tend = 0._wp
  END IF

  IF(.NOT. ALLOCATED(this%column)) ALLOCATE(this%column(nproma,nlev,nblks))
END SUBROUTINE linoz_init_arrays


SUBROUTINE linoz_fill_init(this, jg, idx_tracer)
!<
! SUBROUTINE linoz_fill_init
! Fill the initialised structures with values
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter, KIT
! Initial Release: aronund 2018-10
! Modifications:
!>
  CLASS(t_chem_meta_linoz),INTENT(inout) :: &
    &  this                !< Container with fields
  INTEGER, INTENT(in) :: &
    &  jg,               & !< patch id
    &  idx_tracer          !< index of the tracer

  ! local variables
  INTEGER ::   &
    &  ierror  !< flag if meta information exists in storage
  CHARACTER(LEN = IART_VARNAMELEN) ::  &
    &  O3_paramet,                     & !< read element from xml for how to solve linoz tracer
    &  tracer_name,                    & !< name of the tracer
    &  polarchem                         !< read element from xml for how to treat polar chemistry
  TYPE(t_art_chem_indices), POINTER :: &
    &  art_indices         !< Pointer to ART chem indices
  CHARACTER(:), ALLOCATABLE         :: &
    &  c_tmp

  art_indices => p_art_data(jg)%chem%indices
  
  CALL key_value_storage_as_string(this%opt_meta,'parametrization',c_tmp,ierror)
  IF (ierror == SUCCESS) THEN
    WRITE(O3_paramet,'(A)') c_tmp
    SELECT CASE (TRIM(O3_paramet))
      CASE('analytic')
        this%O3_paramet = IART_LINOZ_ANA
      CASE DEFAULT
        this%O3_paramet = IART_CHEM_NO
    END SELECT
  ELSE
    this%O3_paramet = IART_CHEM_NO
  END IF

  CALL this%opt_meta%get('feedback',this%O3_feed,ierror)
  IF (ierror /= SUCCESS) this%O3_feed = 0
  IF (this%O3_feed == 1 ) THEN
    IF ( irad_o3 /= 10 ) THEN
      CALL finish('mo_art_chem_types:linozv3_fill_init',            &
              &     'ozone feedback set, but irado3 /= 10. ' &
              &   //' Change NAMELIST_ accordingly.')
    ELSE
      art_config(jg)%O3_feedback = this%O3_feed
      CALL message( 'linozv3_fill_init', 'O3 feedback activated.' )
    ENDIF
  ENDIF

  CALL key_value_storage_as_string(this%opt_meta,'polarchem', c_tmp, ierror)
  IF (ierror /= SUCCESS) THEN
    polarchem = 'on'
  ELSE
    WRITE(polarchem,'(A)') c_tmp
  END IF

  SELECT CASE(TRIM(polarchem))
    CASE('off')
      this%polarchem = IART_CHEM_NO
    CASE('on')
      IF (art_indices%iTR_cold /= 0) THEN
        this%polarchem = IART_POLARCHEM
      ELSE
        this%polarchem = IART_LINOZ_LT
      END IF
    CASE('lifetime')
      this%polarchem = IART_LINOZ_LT
    CASE('coldtracer')
      IF (art_indices%iTR_cold /= 0) THEN
        this%polarchem = IART_POLARCHEM
      ELSE
        CALL finish('mo_art_chem_types:linoz_fill_init',            &
              &     'No Cold tracer given. Please include it into ' &
              &   //TRIM(art_config(jg)%cart_chemtracer_xml)//'.')
      END IF
    CASE DEFAULT
      CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
      WRITE(tracer_name,'(A)') c_tmp
      CALL finish('mo_art_chem_types:linoz_fill_init',             &
             &    'Unknown polarchem mode for '//TRIM(tracer_name) &
             &  //':'//TRIM(polarchem)//'. Should be one of '      &
             &  //'on, off, lifetime or  coldtracer.')
  END SELECT

  IF (this%polarchem == IART_POLARCHEM) THEN
    IF (idx_tracer < art_indices%iTR_cold) THEN
      CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
      WRITE(tracer_name,'(A)') c_tmp
      CALL finish('mo_art_chem_types:linoz_fill_init',    &
              &   'Cold tracer has to be set before '     &
              &  //TRIM(tracer_name)//' in xml file.')
    END IF
  END IF

  CALL this%opt_meta%get('lt_het',this%o3lt_het,ierror)
  IF (ierror /= SUCCESS) this%o3lt_het = 864000._wp


  CALL this%opt_meta%get('Thet',this%Thet,ierror)
  IF (ierror /= SUCCESS) this%Thet = 195._wp

  CALL this%opt_meta%get('o3_lbc',this%o3_lbc,ierror)
  IF (ierror /= SUCCESS) this%o3_lbc = 25.E-9_wp


  CALL this%opt_meta%get('lt_lbc',this%lt_lbc,ierror)
  IF (ierror /= SUCCESS) this%lt_lbc = 172800._wp
END SUBROUTINE linoz_fill_init

SUBROUTINE set_tracer_linoz(this,p_tracer_now, p_tracer_cold)
!<
! SUBROUTINE set_tracer_linoz
! Set current tracer values to internal structures (nnow and nnew problem)
! For linoz, this subroutine has to be generated because in case of an existing
! cold tracer, it has to be included, too.
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_linoz),INTENT(inout) :: &
    &  this                           !< Container with fields
  REAL(wp), INTENT(IN), TARGET ::  &
    &  p_tracer_now(:,:,:),        &  !< Mass mixing ratio of tracer (kg/kg)
    &  p_tracer_cold(:,:,:)           !< concentration of the cold tracer
  
  CALL this%set_tracer(p_tracer_now)

  this%cold_tracer => p_tracer_cold
END SUBROUTINE set_tracer_linoz


SUBROUTINE linoz_integrate(this)
!<
! SUBROUTINE linoz_integrate
! Add the calculated tendency to the tracer value
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter, KIT
! Initial Release: aronund 2018-10
! Modifications:
!>
  CLASS(t_chem_meta_linoz),INTENT(inout) :: &
    &  this                        !< Container with fields

  this%tracer = this%tracer + this%tend

END SUBROUTINE linoz_integrate



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Linozv3 tracer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE linozv3_init_arrays(this,nproma,nlev,nblks)
!<
! Adapation of SUBROUTINE linoz_init_arrays to LinozV3
!>
  CLASS(t_chem_meta_linozv3), INTENT(inout) :: &
    &  this               !< container with fields
  INTEGER, INTENT(in) :: &
    &  nproma,nlev,nblks  !< dimensions

  CALL this%init_param(nproma,nlev,nblks)

  IF(.NOT. ALLOCATED(this%tend)) THEN
    ALLOCATE(this%tend(nproma,nlev,nblks))
    this%tend = 0._wp
  END IF

  IF(.NOT. ALLOCATED(this%column)) ALLOCATE(this%column(nproma,nlev,nblks))
  IF(.NOT. ALLOCATED(this%tend_NOy)) ALLOCATE(this%tend_NOy(nproma,nlev,nblks))
  this%tend_NOy = 0._wp
END SUBROUTINE linozv3_init_arrays


SUBROUTINE linozv3_fill_init(this, jg, idx_tracer)
!<
! Adaptation of SUBROUTINE linoz_fill_init to LinozV3
!>
  CLASS(t_chem_meta_linozv3),INTENT(INOUT) :: &
    &  this                !< Container with fields
  INTEGER, INTENT(IN) :: &
    &  jg,               & !< patch id
    &  idx_tracer          !< index of the tracer

  ! local variables
  INTEGER ::   &
    &  ierror, &!< flag if meta information exists in storage
    &  ntime    !< length of f107 series
  REAL(wp) ::   &
    &  lt_tr   !< lifetime of the tracer (s)
  CHARACTER(LEN = IART_VARNAMELEN) ::  &
    &  O3_paramet,                     & !< read element from xml for how to solve linozv3 tracer
    &  tracer_name,                    & !< name of the tracer
    &  polarchem,                      & !< read element from xml for how to treat polar chemistry
    &  solar_ssi                         !< read element from xml for how to treat solar variability
  TYPE(t_art_chem_indices), POINTER :: &
    &  art_indices         !< Pointer to ART chem indices
  CHARACTER(:), ALLOCATABLE         :: &
    &  c_tmp, c_dum

  !-------------------------------------------------------------------------------------------------
  REAL(wp),POINTER               :: zvar1d(:),avar1d(:)
  CHARACTER(LEN = 256)             :: cfname   ! file name containing variables
  INTEGER                        :: ifile_id


  art_indices => p_art_data(jg)%chem%indices
  
  CALL this%opt_meta%get('lifetime',lt_tr,ierror)
  IF (ierror == SUCCESS) THEN
    this%des = 1.0_wp/lt_tr
  ELSE
    CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
    WRITE(tracer_name,'(A)') c_tmp
    CALL finish('mo_art_chem_types:linozv3_fill_init',      &
            &   'lifetime for '//TRIM(tracer_name)//' missing.')
  END IF

  CALL key_value_storage_as_string(this%opt_meta,'parametrization',c_tmp,ierror)
  WRITE(O3_paramet,'(A)') c_tmp
  IF (ierror == SUCCESS) THEN
    SELECT CASE (TRIM(O3_paramet))
      CASE('analytic')
        this%O3_paramet = IART_LINOZ_ANA
      CASE DEFAULT
        this%O3_paramet = IART_CHEM_NO
    END SELECT
  ELSE
    this%O3_paramet = IART_CHEM_NO
  END IF

  CALL key_value_storage_as_string(this%opt_meta,'solar_ssi',c_dum,ierror)
  IF (ierror == SUCCESS) THEN
    WRITE(solar_ssi,'(A)') c_dum
    SELECT CASE (TRIM(solar_ssi))
      CASE('solar mean')
        this%sol_ssi = 0
      CASE('solar min')
        this%sol_ssi = 1
      CASE('solar max')
        this%sol_ssi = 2
      CASE('solar cmip')
        this%sol_ssi = 3
      CASE DEFAULT
        CALL finish('mo_art_chem_types:linozv3_fill_init',      &
            &   'No valid option for solar_ssi.')
    END SELECT
  ELSE
    this%sol_ssi = 1
  END IF

  CALL this%opt_meta%get('feedback',this%O3_feed,ierror)
  IF (ierror /= SUCCESS) this%O3_feed = 0
  IF (this%O3_feed == 1 ) THEN
    IF ( irad_o3 /= 10 ) THEN
      CALL finish('mo_art_chem_types:linozv3_fill_init',            &
              &     'ozone feedback set, but irado3 /= 10. ' &
              &   //' Change NAMELIST_ accordingly.')
    ELSE
      art_config(jg)%O3_feedback = this%O3_feed
      CALL message( 'linozv3_fill_init', 'O3 feedback activated.' )
    ENDIF
  ENDIF

  CALL key_value_storage_as_string(this%opt_meta,'polarchem', c_tmp, ierror)
  WRITE(polarchem,'(A)') c_tmp
  IF (ierror /= SUCCESS) polarchem = 'on'


  SELECT CASE(TRIM(polarchem))
    CASE('off')
      this%polarchem = IART_CHEM_NO
    CASE('on')
      IF (art_indices%iTR_cold /= 0) THEN
        this%polarchem = IART_POLARCHEM
      ELSE
        this%polarchem = IART_LINOZ_LT
      END IF
    CASE('lifetime')
      this%polarchem = IART_LINOZ_LT
    CASE('coldtracer')
      IF (art_indices%iTR_cold /= 0) THEN
        this%polarchem = IART_POLARCHEM
      ELSE
        CALL finish('mo_art_chem_types:linozv3_fill_init',            &
              &     'No Cold tracer given. Please include it into ' &
              &   //TRIM(art_config(jg)%cart_chemtracer_xml)//'.')
      END IF
    CASE DEFAULT
      CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
      WRITE(tracer_name,'(A)') c_tmp
      CALL finish('mo_art_chem_types:linozv3_fill_init',             &
             &    'Unknown polarchem mode for '//TRIM(tracer_name) &
             &  //':'//TRIM(polarchem)//'. Should be one of '      &
             &  //'on, off, lifetime or  coldtracer.')
  END SELECT

  IF (this%polarchem == IART_POLARCHEM) THEN
    IF (idx_tracer < art_indices%iTR_cold) THEN
      CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
      WRITE(tracer_name,'(A)') c_tmp
      CALL finish('mo_art_chem_types:linozv3_fill_init',    &
              &   'Cold tracer has to be set before '     &
              &  //TRIM(tracer_name)//' in xml file.')
    END IF
  END IF

  CALL this%opt_meta%get('lt_het',this%o3lt_het,ierror)
  IF (ierror /= SUCCESS) this%o3lt_het = 864000._wp

  CALL this%opt_meta%get('Thet',this%Thet,ierror)
  IF (ierror /= SUCCESS) this%Thet = 195._wp

  IF ( this%sol_ssi == 3) THEN ! variable solar ssi
    !--- read solar f107      
    cfname = 'solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_18500101-22991231.nc'
              
    CALL openInputFile(ifile_id, cfname)

    CALL read_1D(file_id=ifile_id, variable_name='time',&
    &                 return_pointer=zvar1d )
    ntime = SIZE(zvar1d)
    CALL read_1D(file_id=ifile_id, variable_name='f107',&
    &                 return_pointer=avar1d )
    CALL closeFile(ifile_id)

    ALLOCATE( this%f107(ntime,2))
    this%f107(:,1) = zvar1d
    this%f107(:,2) = avar1d
    CALL message( 'linozv3_fill_init', 'f107 file '// cfname// ' with  records read.' )
    DEALLOCATE( zvar1d, avar1d )
  ENDIF

END SUBROUTINE linozv3_fill_init

SUBROUTINE set_tracer_linozv3(this,p_tracer_now,p_tracer_NOy, p_tracer_cold)
!<
! Adaption of SUBROUTINE set_tracer_linoz to LinozV3
!>
  CLASS(t_chem_meta_linozv3),INTENT(INOUT) :: &
    &  this                           !< Container with fields
  REAL(wp), INTENT(IN), TARGET ::  &
    &  p_tracer_now(:,:,:),        &  !< Mass mixing ratio of tracer (kg/kg)
    &  p_tracer_NOy(:,:,:),        &  !< Mass mixing ratio of tracer (kg/kg)
    &  p_tracer_cold(:,:,:)           !< concentration of the cold tracer
  
  CALL this%set_tracer(p_tracer_now)

  this%cold_tracer => p_tracer_cold
  this%NOy_tracer => p_tracer_NOy
END SUBROUTINE set_tracer_linozv3


SUBROUTINE linozv3_integrate(this)
!<
! Adapation of SUBROUTINE linoz_integrate to LinozV3
!>
  CLASS(t_chem_meta_linozv3),INTENT(INOUT) :: &
    &  this                        !< Container with fields

  this%tracer = this%tracer + this%tend
  !in this version of linozv3 tend of NOy is calculated by simnoy
  !this%NOy_tracer = this%NOy_tracer + this%tend_NOy
  WHERE( this%tracer < 0._wp ) this%tracer = 0._wp

END SUBROUTINE linozv3_integrate

SUBROUTINE convert_vmr_to_mmrv3(this)
!<
! Adapation of SUBROUTINE convert_vmr_to_mmr to LinozV3
!>
  CLASS(t_chem_meta_linozv3),INTENT(INOUT) :: &
    &  this           !< Container with fields
  REAL(wp) ::  &
    &  mol_weight     !< molar weight of the tracer (kg / mol)
  INTEGER ::   &
    &  ierror         !< flag if element exists in storage

  CALL this%opt_meta%get('mol_weight',mol_weight,ierror)
   
  this%tracer = this%tracer/(mol_weight*1000._wp/ amd)
  !this%NOy_tracer = this%NOy_tracer/(mol_weight*1000._wp/ amd)
END SUBROUTINE  convert_vmr_to_mmrv3


SUBROUTINE convert_mmr_to_vmrv3(this)
!<
! SUBROUTINE convert_mmr_to_vmr
! Convert the mass mixing ratio to volume mixing ratio
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter, KIT
! Initial Release: aronund 2018-10
! Modifications:
!>
  CLASS(t_chem_meta_linozv3),INTENT(INOUT) :: &
    &  this           !< Container with fields
  REAL(wp) ::  &
    &  mol_weight     !< molar weight of the tracer (kg / mol)
  INTEGER ::   &
    &  ierror         !< flag if element exists in storage

  CALL this%opt_meta%get('mol_weight',mol_weight,ierror)
  this%tracer = this%tracer * (mol_weight*1000._wp/ amd)
  this%NOy_tracer = this%NOy_tracer * (mol_weight*1000._wp/ amd)
END SUBROUTINE  convert_mmr_to_vmrv3

!##############################################################
! Simnoy routines
!##############################################################

SUBROUTINE simnoy_init_arrays(this, nproma, nlev, nblks)
!<
! SUBROUTINE simnoy_init_arrays
! Initialise the strcutures for simnoy tracer
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_simnoy), INTENT(INOUT) :: &
    &  this                     !< container with fields
  INTEGER, INTENT(in) :: &
    &  nproma, nlev, nblks  !< dimensions

  CALL this%init_param(nproma,nlev,nblks)
  
  IF (.NOT. ALLOCATED(this%tend)) ALLOCATE(this%tend(nproma,nlev,nblks))
  IF (.NOT. ALLOCATED(this%tend_n2o)) ALLOCATE(this%tend_n2o(nproma,nlev,nblks))
  IF (.NOT. ALLOCATED(this%n2onoy_tab1)) ALLOCATE(this%n2onoy_tab1(nproma,nlev))
  IF (.NOT. ALLOCATED(this%n2onoy_tab2)) ALLOCATE(this%n2onoy_tab2(nproma,nlev))
  IF (.NOT. ALLOCATED(this%n2onoy_tab3)) ALLOCATE(this%n2onoy_tab3(nproma,nlev))
  IF (.NOT. ALLOCATED(this%n2onoy_tab4)) ALLOCATE(this%n2onoy_tab4(nproma,nlev))

  this%tend(:,:,:) = 0.0_wp
  this%tend_n2o(:,:,:) = 0.0_wp
END SUBROUTINE simnoy_init_arrays


SUBROUTINE simnoy_fill_init(this, jg, p_prog_list)
!<
! SUBROUTINE simnoy_init_arrays
! Fill the initialised structures for simnoy tracer with values
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_simnoy), INTENT(INOUT) :: &
    &  this              !< container with fields
  INTEGER, INTENT(IN) :: &
    &  jg                !< patch id
  TYPE(t_var_list_ptr), INTENT(IN) :: &
    &  p_prog_list       !< list of prognostic variables
  CHARACTER(LEN = IART_VARNAMELEN) ::  &
    &  tracer_name,                    &  !< name of the tracer
    &  cnoyn2o_tropo,                  &  !< element in xml for how to treat tropospheric param.
    &  polarchem,                      &  !< element in xml for how to treat the polar chemistry
    &  ubc                                !< element in xml controls upper boundary treatment
  INTEGER :: &
    &  ierror           !< flag if element exists in storage
  REAL(wp) ::      &
    &  lt_noysed, lt_tr !< sedimentation lifetime of NOy and lifetime of the tracer (s)
  TYPE(t_var_metadata_dynamic), TARGET :: &
    &  info_dyn_N2O,              &  !< dynamic meta information of N2O, NOy and
    &  info_dyn_NOy,              &  !  cold tracer. They are used to check if
    &  info_dyn_cold                 !  the types are correct
  TYPE(t_art_chem_indices), POINTER :: &
    &  art_indices       !< pointer to ART chem indices
  CHARACTER(:), ALLOCATABLE    :: &
    &  c_tmp
  INTEGER :: iv

  TYPE(t_chem_meta_cold), POINTER :: tracer_ptr

  art_indices => p_art_data(jg)%chem%indices

  CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
  WRITE(tracer_name,'(A)') c_tmp

  IF (TRIM(tracer_name) == 'TRNOy') THEN
    IF (art_indices%iTRN2O /= 0) THEN
      !CALL get_tracer_info_dyn_by_idx(p_prog_list,art_indices%iTRN2O,info_dyn_N2O)
      DO iv = 1, p_prog_list%p%nvars
        IF(p_prog_list%p%vl(iv)%p%info%ncontained /= art_indices%iTRN2O) CYCLE
        info_dyn_N2O = p_prog_list%p%vl(iv)%p%info_dyn
      END DO

      SELECT TYPE (tracer => info_dyn_N2O%tracer)
        TYPE IS (t_chem_meta_simnoy)
          CALL tracer%opt_meta%get('lifetime', lt_tr)
          this%des_n2o = 1._wp / lt_tr
        CLASS DEFAULT
          CALL finish('mo_art_chem_types:simnoy_fill_init',      &
            &   'c_solve of N2O is not simnoy. Please correct it.')
      END SELECT
    ELSE
      CALL finish('mo_art_chem_types:simnoy_fill_init',      &
              &   'N2O is not present in tracers, although ' &
              & //'simnoy is chosen as parametrisation of N2O and NOy.')
    END IF

    CALL n2onoy_read()

    CALL key_value_storage_as_string(this%opt_meta,'polarchem',c_tmp, ierror)
    IF (ierror == SUCCESS) THEN
      WRITE(polarchem,'(A)') c_tmp
    ELSE
      IF (art_indices%iTR_cold /= 0) THEN
        polarchem = 'on'
      ELSE
        polarchem = 'off'
      END IF
    END IF

    SELECT CASE (TRIM(polarchem))
      CASE('on')
        this%polarchem = IART_POLARCHEM
      CASE('off')
        this%polarchem = IART_CHEM_NO
      CASE('noy_sedimentation')
        this%polarchem = IART_SIMNOY_SEDI
      CASE DEFAULT
        CALL finish('mo_art_chem_types:simnoy_fill_init',  &
                &   'Unknown Parameter for polarchem: '//TRIM(polarchem)  &
                & //', must be one of ''on'', ''off'', or ''noy_sedimentation''.')
    END SELECT

    IF (this%polarchem /= IART_CHEM_NO) THEN
      IF (art_indices%iTR_cold == 0) THEN
        CALL finish('mo_art_chem_types:simnoy_fill_init',  &
                &   'Cold tracer not present although polarchem' &
                & //' is set to ''on'' or ''noy_sedimentation''.')
      ELSE IF (art_indices%iTR_cold > art_indices%iTRNOy) THEN
        CALL finish('mo_art_chem_types:simnoy_fill_init',  &
                &   'Cold tracer must be set before NOy in the xml file.')
      ELSE
        IF (this%polarchem == IART_POLARCHEM) THEN
          !CALL get_tracer_info_dyn_by_idx(p_prog_list,art_indices%iTR_cold,info_dyn_cold)
          DO iv = 1, p_prog_list%p%nvars
            IF(p_prog_list%p%vl(iv)%p%info%ncontained /= art_indices%iTR_cold) CYCLE
            info_dyn_cold = p_prog_list%p%vl(iv)%p%info_dyn
          END DO
          SELECT TYPE (tracer => info_dyn_cold%tracer)
            TYPE IS (t_chem_meta_cold)
              IF (tracer%polarchem == IART_POLARCHEM) THEN
                !this%p_cold_sed => tracer%p_cold_sed
                tracer_ptr => tracer
                this%p_cold_sed => tracer_ptr%p_cold_sed
              ELSE
                CALL finish('mo_art_chem_types:simnoy_fill_init',  &
                       &    'Cold tracer must have polchem ''on'' '&
                       &  //'when combined with simnoy polarchem ''on''.')
              END IF
          END SELECT
        END IF
      END IF
    END IF
  
    CALL this%opt_meta%get('lt_sed', lt_noysed,ierror)
    IF (ierror /= SUCCESS) THEN
      IF (this%polarchem == IART_SIMNOY_SEDI) THEN
        lt_noysed = 2._wp * 86400._wp
      ELSE
        lt_noysed = 1.e30_wp
      END IF
    END IF
    this%des_noysed = 1._wp / lt_noysed

    CALL key_value_storage_as_string(this%opt_meta,'cnoyn2o_tropo', c_tmp, ierror)
    IF (ierror == SUCCESS) THEN
      WRITE(cnoyn2o_tropo,'(A)') c_tmp
      SELECT CASE (TRIM(cnoyn2o_tropo))
        CASE('EXTP')
          this%cnoyn2o_tropo = IART_SIMNOY_EXTP
        CASE('WMO')
          this%cnoyn2o_tropo = IART_SIMNOY_WMO
        CASE('PRES')
          this%cnoyn2o_tropo = IART_SIMNOY_PRES
        CASE DEFAULT
          CALL finish('mo_art_chem_types_param:simnoy_fill_init',  &
                 &     'cnoyn2o_tropo must be one of EXTP, WMO, or PRES.')
      END SELECT
    ELSE
      this%cnoyn2o_tropo = IART_SIMNOY_EXTP
    END IF

    CALL key_value_storage_as_string(this%opt_meta,'ubc', c_tmp, ierror)
    IF (ierror == SUCCESS) THEN
      WRITE(ubc,'(A)') c_tmp
      SELECT CASE (TRIM(ubc))
        CASE('UBC')
          this%ubc = IART_SIMNOY_UBC
        CASE('NOUBC')
          this%ubc = IART_SIMNOY_NOUBC
        CASE DEFAULT
          CALL finish('mo_art_chem_types_param:simnoy_fill_init',  &
                 &     'ubc must be one of UBC, NOUBC.')
      END SELECT
    ELSE
      this%ubc = IART_SIMNOY_NOUBC
    END IF

    CALL this%opt_meta%get('lifetime', lt_tr, ierror)
    IF (ierror == SUCCESS) THEN
      this%des = 1.0_wp/lt_tr
    ELSE
      CALL finish('mo_art_chem_types:simnoy_fill_init',      &
              &   'lifetime for TRNOy missing.')
    END IF

    

  ELSE IF (TRIM(tracer_name) == 'TRN2O') THEN

    CALL this%opt_meta%get('lifetime', lt_tr, ierror)
    IF (ierror == SUCCESS) THEN
      this%des = 1.0_wp/lt_tr
    ELSE
      CALL finish('mo_art_chem_types:simnoy_fill_init',      &
              &   'lifetime for TRN2O missing.')
    END IF

    ! Check for NOy
    IF (art_indices%iTRNOy /= 0) THEN
      !CALL get_tracer_info_dyn_by_idx(p_prog_list,art_indices%iTRNOy,info_dyn_NOy)
      DO iv = 1, p_prog_list%p%nvars
        IF(p_prog_list%p%vl(iv)%p%info%ncontained /= art_indices%iTRNOy) CYCLE
        info_dyn_NOy = p_prog_list%p%vl(iv)%p%info_dyn
      END DO

      SELECT TYPE (tracer => info_dyn_NOy%tracer)
        TYPE IS (t_chem_meta_simnoy)
          CALL tracer%opt_meta%get('lifetime', lt_tr)
        CLASS DEFAULT
          CALL finish('mo_art_chem_types:simnoy_fill_init',      &
            &   'c_solve of NOy is not simnoy. Please correct it.')
      END SELECT
    ELSE
      CALL finish('mo_art_chem_types:simnoy_fill_init',      &
              &   'NOy is not present in tracers, although ' &
              & //'simnoy is chosen as parametrisation of N2O and NOy.')
    END IF
  ELSE
    CALL finish('mo_art_chem_types:simnoy_fill_init',  &
            &   'Only TRNOy and TRN2O can participate in Simnoy.')
  END IF             
END SUBROUTINE simnoy_fill_init


SUBROUTINE set_tracer_simnoy(this,p_tracer_now, p_tracer_n2o, p_tracer_cold, p_prog_list, jg)
!<
! SUBROUTINE set_tracer_simnoy
! Set current tracer values to internal structures (nnow and nnew problem)
! For simnoy, this includes also the cold tracer and N2O and p_cold_sed (if
! allocated and used)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE

  CLASS(t_chem_meta_simnoy),INTENT(inout) :: &
    &  this                        !< Container with fields
  REAL(wp), INTENT(IN), TARGET ::  &
    &  p_tracer_now(:,:,:),        &  !< NOy mass mixing ratio (kg/kg)
    &  p_tracer_n2o(:,:,:)            !< N2O mass mixing ratio (kg/kg)
  REAL(wp), INTENT(in), TARGET, OPTIONAL ::  &
    &  p_tracer_cold(:,:,:)           !< Cold tracer concentration
  TYPE(t_var_list_ptr), INTENT(in), OPTIONAL :: &
    &  p_prog_list                    !< List of prognostic variables
  INTEGER, INTENT(in), OPTIONAL :: &
    &  jg                             !< patch id

  ! local
  INTEGER :: &
    &  iv                             !< loop index (list-element)

  TYPE(t_chem_meta_cold),POINTER :: tracer_ptr
  
  CALL this%set_tracer(p_tracer_now)

  this%n2o_tracer => p_tracer_n2o

  IF (PRESENT(p_tracer_cold)) THEN
    this%cold_tracer => p_tracer_cold

    IF (this%polarchem == IART_POLARCHEM) THEN

      DO iv = 1, p_prog_list%p%nvars
        IF (p_prog_list%p%vl(iv)%p%info%ncontained == p_art_data(jg)%chem%indices%iTR_cold) THEN
          SELECT TYPE (tracer => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
            TYPE IS (t_chem_meta_cold)
              !this%p_cold_sed => tracer%p_cold_sed
              tracer_ptr => tracer
              this%p_cold_sed => tracer_ptr%p_cold_sed
              EXIT
          END SELECT
        END IF

      END DO
    END IF
  END IF
END SUBROUTINE set_tracer_simnoy


SUBROUTINE simnoy_integrate(this)
!<
! SUBROUTINE get_tracer_simnoy
! Add the tendencies of NOy and N2O
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  CLASS(t_chem_meta_simnoy), INTENT(inout) :: &
    &  this

  this%tracer = this%tracer + this%tend
  WHERE( this%tracer < 0._wp ) this%tracer = 0._wp
  this%n2o_tracer = this%n2o_tracer + this%tend_n2o
  WHERE( this%n2o_tracer < 0._wp ) this%n2o_tracer = 0._wp
END SUBROUTINE simnoy_integrate


END MODULE mo_art_chem_types_param
  
