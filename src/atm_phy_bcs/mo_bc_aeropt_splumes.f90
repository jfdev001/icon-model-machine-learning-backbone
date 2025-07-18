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

! Read and apply monthly aerosol optical properties of S. Kinne
! from yearly files.

! ---------------------------
#include "consistent_fma.inc"
! ---------------------------

MODULE mo_bc_aeropt_splumes

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
  USE mo_read_interface,       ONLY: openInputFile, read_1D, &
                                   & read_bcast_real_2D, read_bcast_real_3D, &
                                   & closeFile
  USE mo_model_domain,         ONLY: p_patch
  USE mo_fortran_tools,        ONLY: assert_acc_device_only
  USE mo_math_constants,       ONLY: rad2deg
  USE mtime,                   ONLY: datetime, getDayOfYearFromDateTime, &
       &                             getNoOfDaysInYearDateTime

  IMPLICIT NONE

  PRIVATE
  PUBLIC                  :: setup_bc_aeropt_splumes, add_bc_aeropt_splumes

  INTEGER, PARAMETER      ::     &
       nplumes   = 9            ,& !< Number of plumes
       nfeatures = 2            ,& !< Number of features per plume
       ntimes    = 52           ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251             !< Number of years of available forcing
  CHARACTER(LEN=*), PARAMETER :: cfname = 'MACv2.0-SP_v1.nc'

  REAL(wp), ALLOCATABLE ::                    &
       plume_lat   (:)  ,& !< (nplumes) latitude where plume maximizes
       plume_lon   (:)  ,& !< (nplumes) longitude where plume maximizes
       beta_a      (:)  ,& !< (nplumes) parameter a for beta function
                           !< vertical profile
       beta_b      (:)  ,& !< (nplumes) parameter b for beta function
                           !< vertical profile
       aod_spmx    (:)  ,& !< (nplumes) aod at 550 for simple plume (maximum)
       aod_fmbg    (:)  ,& !< (nplumes) aod at 550 for fine mode
                           !< natural background (for twomey effect)
       asy550      (:)  ,& !< (nplumes) asymmetry parameter for plume at 550nm
       ssa550      (:)  ,& !< (nplumes) single scattering albedo for
                           !< plume at 550nm
       angstrom    (:)  ,& !< (nplumes) angstrom parameter for plume
       sig_lon_E   (:,:),& !< (nfeatures,nplumes) Eastward extent of
                           !< plume feature
       sig_lon_W   (:,:),& !< (nfeatures,nplumes) Westward extent of
                           !< plume feature
       sig_lat_E   (:,:),& !< (nfeatures,nplumes) Southward extent of
                           !< plume feature
       sig_lat_W   (:,:),& !< (nfeatures,nplumes) Northward extent of
                           !< plume feature
       theta       (:,:),& !< (nfeatures,nplumes) Rotation angle of feature
       ftr_weight  (:,:),& !< (nfeatures,nplumes) Feature weights =
                           !< (nfeatures + 1) to account for BB background
       year_weight (:,:)    ,& !< (nyear,nplumes) Yearly weight for plume
       ann_cycle   (:,:,:)     !< (nfeatures,ntimes,nplumes) annual cycle for feature
  LOGICAL                  :: sp_initialized = .FALSE.

  !$ACC DECLARE CREATE(angstrom, asy550, ssa550)

  CONTAINS

  ! -----------------------------------------------------------------
  ! setup_bc_aeropt_splumes:  This subroutine should be called at initialization to
  !            read the netcdf data that describes the simple plume
  !            climatology.  The information needs to be either read
  !            by each processor or distributed to processors.
  !
  SUBROUTINE setup_bc_aeropt_splumes
    !
    ! ----------
    !
    INTEGER           :: ifile_id

    CALL openInputFile(ifile_id, cfname)

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lat',&
                       & alloc_array=plume_lat,    file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lon',&
                       & alloc_array=plume_lon,    file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_a',   &
                       & alloc_array=beta_a,    file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_b',   &
                       & alloc_array=beta_b,    file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_spmx', &
                       & alloc_array=aod_spmx,    file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_fmbg', &
                       & alloc_array=aod_fmbg,    file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='ssa550',   &
                       & alloc_array=ssa550,    file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='asy550',   &
                       & alloc_array=asy550,    file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='angstrom', &
                       & alloc_array=angstrom,    file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_W',&
                       & alloc_array=sig_lat_W,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_E',&
                       & alloc_array=sig_lat_E,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_W',&
                       & alloc_array=sig_lon_W,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_E',&
                       & alloc_array=sig_lon_E,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='theta',    &
                       & alloc_array=theta,    file_name=cfname,             &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='setup_bc_aeropt_splumes'             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='ftr_weight',&
                       & alloc_array=ftr_weight,    file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),                &
                       & module_name='mo_bc_aeropt_splumes',                  &
                       & sub_prog_name='setup_bc_aeropt_splumes'              )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='year_weight',&
                       & alloc_array=year_weight,    file_name=cfname,         &
                       & variable_dimls=(/nyears,nplumes/),                    &
                       & module_name='mo_bc_aeropt_splumes',                   &
                       & sub_prog_name='setup_bc_aeropt_splumes'               )
    CALL read_3d_wrapper(ifile_id=ifile_id,        variable_name='ann_cycle', &
                       & alloc_array=ann_cycle,    file_name=cfname,          &
                       & variable_dimls=(/nfeatures,ntimes,nplumes/),         &
                       & module_name='mo_bc_aeropt_splumes',                  &
                       & sub_prog_name='setup_bc_aeropt_splumes'               )
    CALL closeFile(ifile_id)

    sp_initialized = .TRUE.

    !$ACC ENTER DATA CREATE(plume_lat, plume_lon, beta_a, beta_b, aod_spmx, aod_fmbg, ssa550) &
    !$ACC   CREATE(asy550, angstrom, sig_lat_W, sig_lat_E, sig_lon_W, sig_lon_E, theta) &
    !$ACC   CREATE(ftr_weight, year_weight, ann_cycle)

    !$ACC UPDATE DEVICE(plume_lat, plume_lon, beta_a, beta_b, aod_spmx, aod_fmbg, ssa550) &
    !$ACC   DEVICE(asy550, angstrom, sig_lat_W, sig_lat_E, sig_lon_W, sig_lon_E, theta) &
    !$ACC   DEVICE(ftr_weight, year_weight, ann_cycle) ASYNC(1)

  END SUBROUTINE setup_bc_aeropt_splumes
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
  ! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
  ! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
  ! week in the year and superimposed on the yearly mean value of the weight.
  !
  SUBROUTINE set_time_weight(current_dt, time_weight, time_weight_bg)
    !
    ! ----------
    !
    TYPE(datetime), INTENT(IN) :: current_dt !< Current date and time
    REAL(wp), INTENT(OUT) :: time_weight(nfeatures,nplumes) !< Time-weights to account for BB background
    REAL(wp), INTENT(OUT) :: time_weight_bg(nfeatures,nplumes) !< as time_weight but for natural background in Twomey effect

    INTEGER          ::  &
         iyear          ,& !< Integer year values between 1 and 156 (1850-2100)
         iweek          ,& !< Integer index (between 1 and ntimes); for ntimes=52 this corresponds to weeks (roughly)
         iplume            ! plume number
    !
    ! ----------
    !

    iyear = INT(current_dt%date%year) - 1849
    iweek = 1 + FLOOR( &
        & REAL(getDayOfYearFromDateTime(current_dt) - 1, wp) / REAL(getNoOfDaysInYearDateTime(current_dt), wp) * ntimes &
      )

    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) THEN
      CALL finish('mo_bc_aeropt_splumes:set_time_weight', 'time index out of bounds')
    END IF

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
    DO iplume=1,nplumes
      time_weight(1,iplume) = year_weight(iyear,iplume) * ann_cycle(1,iweek,iplume)
      time_weight(2,iplume) = year_weight(iyear,iplume) * ann_cycle(2,iweek,iplume)
      time_weight_bg(1,iplume) = ann_cycle(1,iweek,iplume)
      time_weight_bg(2,iplume) = ann_cycle(2,iweek,iplume)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE set_time_weight
  !
  ! ---------------------------------------------------------------------------------------------
  ! sp_plume_profile_550:
  !                  This subroutine calculates the simple plume aerosol and cloud active optical
  !                  properites at 550nm based on the the simple plume fit to the MPI Aerosol
  !                  Climatology (Version 2).
  !
  SUBROUTINE sp_plume_profile_550( &
        & nlevels, jcs, jce, nproma, time_weight, time_weight_bg, z, dz, oro, lon, lat, aod_550, dNovrN &
      )

    INTEGER, INTENT(IN) :: nlevels !< number of levels
    INTEGER, INTENT(IN) :: jcs !< start index in block
    INTEGER, INTENT(IN) :: jce !< end index in block
    INTEGER, INTENT(IN) :: nproma !< first dimension of 2d-vars as declared in calling (sub)program

    REAL(wp), INTENT(IN) :: time_weight(nfeatures,nplumes) !< Time-weights to account for BB background.
    REAL(wp), INTENT(IN) :: time_weight_bg(nfeatures,nplumes) !< as time_weight but for natural background in Twomey effect.

    REAL(wp), INTENT(IN) :: z(nproma,nlevels) !< height above sea-level (m)
    REAL(wp), INTENT(IN) :: dz(nproma,nlevels) !< level thickness (difference between half levels)
    REAL(wp), INTENT(IN) :: oro(nproma) !< orographic height (m)
    REAL(wp), INTENT(IN) :: lon(nproma) !< longitude in degrees E
    REAL(wp), INTENT(IN) :: lat(nproma) !< latitude in degrees N

    REAL(wp), INTENT(OUT) :: aod_550(nproma,nlevels,nplumes) !< AOD at 550nm at each level for each plume.
    REAL(wp), INTENT(INOUT), OPTIONAL :: dNovrN(nproma) !< anthropogenic increment to cloud drop number concentration

    REAL(wp) :: caod_sp(jcs:jce) !< column simple plume (anthropogenic) aod at 550 nm
    REAL(wp) :: caod_bg(jcs:jce) !< column fine-mode natural background aod at 550 nm
    REAL(wp) :: cw_an(jcs:jce) !< column weight for simple plume (anthropogenic) aod at 550 nm
    REAL(wp) :: cw_bg(jcs:jce) !< column weight for fine-mode indurstrial background aod at 550 nm
    REAL(wp) :: prof(jcs:jce,nlevels) !< scaled profile (by beta function)
    REAL(wp) :: beta_sum(jcs:jce) !< vertical sum of beta function

    REAL(wp) :: eta !< normalized height (by 15 km)
    REAL(wp) :: delta_lat !< latitude offset
    REAL(wp) :: delta_lon !< longitude offset
    REAL(wp) :: delta_lon_t !< threshold for maximum longitudinal plume extent used in transition from 360 to 0 degrees
    REAL(wp) :: a_plume1 !< gaussian longitude factor for feature 1
    REAL(wp) :: b_plume1 !< gaussian latitude factor for feature 1
    REAL(wp) :: a_plume2 !< gaussian longitude factor for feature 2
    REAL(wp) :: b_plume2 !< gaussian latitude factor for feature 2
    REAL(wp) :: lon1 !< rotated longitude for feature 1
    REAL(wp) :: lat1 !< rotated latitude for feature 1
    REAL(wp) :: lon2 !< rotated longitude for feature 2
    REAL(wp) :: lat2 !< rotated latitude for feature 2
    REAL(wp) :: f1 !< contribution from feature 1
    REAL(wp) :: f2 !< contribution from feature 2
    REAL(wp) :: f3 !< contribution from feature 1 in natural background of Twomey effect
    REAL(wp) :: f4 !< contribution from feature 2 in natural background of Twomey effect

    INTEGER :: iplume, icol, k
    LOGICAL :: have_dNovrN

    have_dNovrN = PRESENT(dNovrN)

    !$ACC DATA CREATE(caod_sp, caod_bg, cw_an, cw_bg, prof, beta_sum) NO_CREATE(dNovrN)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO icol=jcs,jce
      caod_sp(icol)  = 0.00_wp
      caod_bg(icol)  = 0.02_wp
    END DO

    !
    ! sum contribution from plumes to construct composite profiles of aerosol otpical properties
    !
    !$ACC LOOP SEQ
!PREVENT_INCONSISTENT_IFORT_FMA
    DO iplume=1,nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO icol=jcs,jce
        beta_sum(icol) = 0._wp
      END DO

      !$ACC LOOP SEQ
      DO k=1,nlevels
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(eta)
        DO icol=jcs,jce
          eta = MIN(MAX(0.0_wp, z(icol,k)/15000._wp), 1.0_wp)
          prof(icol,k) = (eta**(beta_a(iplume)-1._wp) * (1._wp-eta)**(beta_b(iplume)-1._wp))*dz(icol,k)
          beta_sum(icol) = beta_sum(icol) + prof(icol,k)
        END DO
      END DO

      !$ACC LOOP SEQ
      DO k=1,nlevels
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO icol=jcs,jce
          IF (z(icol,k) >= oro(icol)) THEN
            prof(icol,k) = prof(icol,k) / beta_sum(icol)
          ELSE
            prof(icol,k) = 0._wp
          END IF
        END DO
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(delta_lat, delta_lon, delta_lon_t, a_plume1) &
      !$ACC   PRIVATE(b_plume1, a_plume2, b_plume2, lon1, lat1, lon2, lat2, f1, f2, f3, f4)
      DO icol=jcs,jce
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat = lat(icol) - plume_lat(iplume)
        delta_lon = lon(icol) - plume_lon(iplume)
        delta_lon_t = MERGE (260._wp, 180._wp, iplume == 1)
        delta_lon = MERGE ( delta_lon-SIGN(360._wp,delta_lon) , delta_lon , ABS(delta_lon) > delta_lon_t)

        a_plume1  = 0.5_wp / (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0.0_wp)**2)
        b_plume1  = 0.5_wp / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0.0_wp)**2)
        a_plume2  = 0.5_wp / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0.0_wp)**2)
        b_plume2  = 0.5_wp / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0.0_wp)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        !
        lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
        lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)
        lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
        lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode background aerosol (cw_bg)
        !
        f1 = time_weight(1,iplume) * ftr_weight(1,iplume) * EXP(-(a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2))))
        f2 = time_weight(2,iplume) * ftr_weight(2,iplume) * EXP(-(a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))
        f3 = time_weight_bg(1,iplume) * ftr_weight(1,iplume) * EXP(-(a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2))))
        f4 = time_weight_bg(2,iplume) * ftr_weight(2,iplume) * EXP(-(a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))


        cw_an(icol) = f1 * aod_spmx(iplume) + f2 * aod_spmx(iplume)
        cw_bg(icol) = f3 * aod_fmbg(iplume) + f4 * aod_fmbg(iplume)
      END DO

      !$ACC LOOP SEQ
      DO k=1,nlevels
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO icol = jcs,jce
          aod_550(icol,k,iplume) = prof(icol,k) * cw_an(icol)
          caod_sp(icol) = caod_sp(icol) + aod_550(icol,k,iplume)
          caod_bg(icol) = caod_bg(icol) + prof(icol,k) * cw_bg(icol)
        END DO
      END DO
    END DO ! iplume

    IF (have_dNovrN) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO icol=jcs,jce
        dNovrN(icol) = LOG((1000.0_wp * (caod_sp(icol) + caod_bg(icol))) + 1.0_wp)/LOG((1000.0_wp * caod_bg(icol)) + 1.0_wp)
      END DO
    END IF

    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE sp_plume_profile_550

  !
  ! ---------------------------------------------------------------------------------------------
  ! sp_aop_profile_wavelength: This routine computes the wavelength-adjusted simple plume aerosol
  !                  optical properties away from 500nm.  It sums over nplumes to provide a
  !                  profile of aerosol optical properties on a host models vertical grid.
  !
  SUBROUTINE sp_aop_profile_wavelength(nlevels, jcs, jce, nproma, lambda, aod_550, aod_prof, ssa_prof, asy_prof)

    !$ACC ROUTINE GANG

    INTEGER, INTENT(IN) :: nlevels !< number of levels
    INTEGER, INTENT(IN) :: jcs !< start index in block
    INTEGER, INTENT(IN) :: jce !< end index in block
    INTEGER, INTENT(IN) :: nproma !< first dimension of 2d-vars as declared in calling (sub)program

    REAL(wp), INTENT(IN) :: lambda !< Wavelength [nm].
    REAL(wp), INTENT(IN) :: aod_550(nproma,nlevels,nplumes) !< AOD at 550nm at each level for each plume.

    REAL(wp), INTENT(OUT) :: aod_prof(nproma,nlevels) !< profile of aerosol optical depth @ lambda
    REAL(wp), INTENT(OUT) :: ssa_prof(nproma,nlevels) !< profile of single_scattering albedo @ lambda
    REAL(wp), INTENT(OUT) :: asy_prof(nproma,nlevels) !< profile of asymmetry parameter @ lambda

    REAL(wp) :: lfactor
    REAL(wp) :: ssa
    REAL(wp) :: asy
    REAL(wp) :: aod_lmd

    INTEGER :: iplume, k, icol

    !$ACC LOOP SEQ
    DO k = 1, nlevels
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(aod_lmd)
      DO icol = jcs, jce
        aod_prof(icol,k) = 0._wp
        ssa_prof(icol,k) = 0._wp
        asy_prof(icol,k) = 0._wp
      END DO
    END DO

    !$ACC LOOP SEQ
!PREVENT_INCONSISTENT_IFORT_FMA
    DO iplume = 1, nplumes
      lfactor = MIN(1.0_wp,700.0_wp/lambda)
      ssa = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1-ssa550(iplume)) * lfactor))
      asy = asy550(iplume) * SQRT(lfactor)

      !$ACC LOOP SEQ
      DO k = 1, nlevels
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(aod_lmd)
        DO icol = jcs, jce
          aod_lmd = aod_550(icol,k,iplume) * EXP(-angstrom(iplume) * LOG(lambda/550.0_wp))
          asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa * asy
          ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa
          aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd
        END DO
      END DO
    END DO

    !
    ! complete optical depth weighting
    !
    !$ACC LOOP SEQ
    DO k = 1, nlevels
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO icol = jcs, jce
        IF (ssa_prof(icol,k) > TINY(1._wp)) THEN
          asy_prof(icol,k) = asy_prof(icol,k)/ssa_prof(icol,k)
        ELSE
          asy_prof(icol,k) = 0.0_wp
        END IF

        IF (aod_prof(icol,k) > TINY(1._wp)) THEN
          ssa_prof(icol,k) = ssa_prof(icol,k)/aod_prof(icol,k)
        ELSE
          ssa_prof(icol,k) = 1.0_wp
        END IF
      END DO
    END DO

  END SUBROUTINE sp_aop_profile_wavelength

  ! -----------------------------------------------------------------------------------------------
  ! add_bc_aeropt_splumes:  This subroutine provides the interface to simple plume (sp) fit to the
  !                         MPI Aerosol Climatology (Version 2). It does so by collecting or
  !                         deriving spatio-temporal information and calling the simple plume
  !                         aerosol subroutine and incrementing the background aerosol properties
  !                         (and effective radius) with the anthropogenic plumes.
  !
  SUBROUTINE add_bc_aeropt_splumes( jg                                             ,&
     & jcs            ,jce            ,nproma         ,klev           ,jb          ,&
     & nb_sw          ,this_datetime  ,zf             ,dz             ,z_sfc       ,&
     & sw_wv1         ,sw_wv2         ,aod_sw_vr      ,ssa_sw_vr      ,asy_sw_vr   ,&
     & x_cdnc         ,lacc                                                         )
    !
    ! --- 0.1 Variables passed through argument list
    INTEGER, INTENT(IN) ::            &
         jg                          ,& !< domain index
         jcs                         ,& !< start index in current block
         jce                         ,& !< end index in current block
         nproma                      ,& !< block dimension
         klev                        ,& !< number of full levels
         jb                          ,& !< index for current block
         nb_sw                          !< number of bands in short wave

    TYPE(datetime), POINTER      :: this_datetime

    REAL(wp), INTENT (IN)        :: &
         zf(nproma,klev),            & !< geometric height at full level [m]
         dz(nproma,klev),            & !< geometric height thickness     [m]
         z_sfc(nproma),              & !< geometric height of surface    [m]
         sw_wv1(nb_sw),              & !< smallest wave number in each of the sw bands
         sw_wv2(nb_sw)                !< largest  wave number in each of the sw bands

    REAL(wp), INTENT (INOUT) ::       &
         aod_sw_vr(nproma,klev,nb_sw) ,& !< Aerosol shortwave optical depth
         ssa_sw_vr(nproma,klev,nb_sw) ,& !< Aerosol single scattering albedo
         asy_sw_vr(nproma,klev,nb_sw)    !< Aerosol asymmetry parameter
    REAL(wp), INTENT(OUT), OPTIONAL:: &
         x_cdnc(nproma)                  !< Scale factor for Cloud Droplet Number Concentration

    LOGICAL, OPTIONAL, INTENT(IN) :: lacc !< OpenACC flag.

    !
    ! --- 0.2 Local variables
    !
    INTEGER ::                        &
         jk                          ,& !< index for looping over vertical dimension
         jki                         ,& !< index for looping over vertical dimension for reversing
         jl                          ,& !< index for looping over block
         jwl                            !< index for looping over wavelengths

    REAL(wp) :: time_weight(nfeatures,nplumes)
    REAL(wp) :: time_weight_bg(nfeatures,nplumes)

    REAL(wp) :: aod_550(nproma,klev,nplumes)

    REAL(wp) ::                       &
         lambda                      ,& !< wavelength at central band wavenumber [nm]
         lon_sp(nproma)              ,& !< longitude passed to sp
         lat_sp(nproma)              ,& !< latitude passed to sp
         z_fl_vr(nproma,klev)        ,& !< level height [m], vertically reversed indexing (1=lowest level)
         dz_vr(nproma,klev)          ,& !< level thickness [m], vertically reversed
         sp_aod_vr(nproma,klev)      ,& !< simple plume aerosol optical depth, vertically reversed
         sp_ssa_vr(nproma,klev)      ,& !< simple plume single scattering albedo, vertically reversed
         sp_asy_vr(nproma,klev)         !< simple plume asymmetry factor, vertically reversed indexing

    CALL assert_acc_device_only('add_bc_aeropt_splumes',lacc)
    !
    ! ----------
    !
    ! initialize input data (by calling setup at first instance)
    !
    IF (.NOT.sp_initialized) CALL setup_bc_aeropt_splumes

    IF (this_datetime%date%year > 1850) THEN

      !$ACC DATA CREATE(time_weight, time_weight_bg, aod_550, lon_sp, lat_sp, z_fl_vr, dz_vr) &
      !$ACC   CREATE(sp_aod_vr, sp_ssa_vr, sp_asy_vr)

      !
      ! get time weights
      !
      CALL set_time_weight(this_datetime, time_weight, time_weight_bg)

      !
      ! --- 1.1 geographic information
      !

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk=1,klev
        jki=klev-jk+1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jl=jcs,jce
          dz_vr  (jl,jk) = dz(jl,jki)
          z_fl_vr(jl,jk) = zf(jl,jki)
        END DO
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jl=jcs,jce
        lon_sp(jl) = p_patch(jg)%cells%center(jl,jb)%lon * rad2deg
        lat_sp(jl) = p_patch(jg)%cells%center(jl,jb)%lat * rad2deg
      END DO
      !$ACC END PARALLEL

      !
      ! --- 1.2 Aerosol Shortwave properties
      !
      ! get aerosol optical properties in each band, and adjust effective radius
      !

      CALL sp_plume_profile_550( &
          & nlevels=klev, &
          & jcs=jcs, &
          & jce=jce, &
          & nproma=nproma, &
          & time_weight=time_weight(:,:), &
          & time_weight_bg=time_weight_bg(:,:), &
          & z=z_fl_vr(:,:), &
          & dz=dz_vr(:,:), &
          & oro=z_sfc(:), &
          & lon=lon_sp(:), &
          & lat=lat_sp(:), &
          & aod_550=aod_550(:,:,:), &
          & dNovrN=x_cdnc(:) &
        )

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jwl = 1,nb_sw
        lambda = 1.e7_wp/ (0.5_wp * (sw_wv1(jwl) + sw_wv2(jwl)))

        CALL sp_aop_profile_wavelength ( &
            & nlevels=klev, &
            & jcs=jcs, &
            & jce=jce, &
            & nproma=nproma, &
            & lambda=lambda, &
            & aod_550=aod_550, &
            & aod_prof=sp_aod_vr, &
            & ssa_prof=sp_ssa_vr, &
            & asy_prof=sp_asy_vr &
          )

        !$ACC LOOP SEQ
        DO jk=1,klev
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jl=jcs,jce
            asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl) * ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)    &
                 + sp_asy_vr(jl,jk)   * sp_ssa_vr(jl,jk)    * sp_aod_vr(jl,jk)
            ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl) * aod_sw_vr(jl,jk,jwl)                           &
                 + sp_ssa_vr(jl,jk)   * sp_aod_vr(jl,jk)
            aod_sw_vr(jl,jk,jwl) = aod_sw_vr(jl,jk,jwl) + sp_aod_vr(jl,jk)
            IF (ssa_sw_vr(jl,jk,jwl) > TINY(1.0_wp)) THEN
              asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl)/ssa_sw_vr(jl,jk,jwl)
            ELSE
              asy_sw_vr(jl,jk,jwl) = asy_sw_vr(jl,jk,jwl)
            END IF

            IF (aod_sw_vr(jl,jk,jwl) > TINY(1.0_wp)) THEN
              ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl)/aod_sw_vr(jl,jk,jwl)
            ELSE
              ssa_sw_vr(jl,jk,jwl) = ssa_sw_vr(jl,jk,jwl)
            END IF

          END DO
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC WAIT(1)
      !$ACC END DATA
    END IF

  END SUBROUTINE add_bc_aeropt_splumes

  SUBROUTINE read_1d_wrapper(ifile_id,                 variable_name,        &
                           & alloc_array,              file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable
                                                    !< to be read
    REAL(wp), ALLOCATABLE ,INTENT(out)  :: alloc_array(:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(1)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length, cj_length
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_1D(file_id=ifile_id,         variable_name=variable_name,      &
                 alloc_array=alloc_array                                     )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(alloc_array,1)/=variable_dimls(1)) THEN
         WRITE(ci_length,*) SIZE(alloc_array,1)
         WRITE(cj_length,*) variable_dimls(1)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('// &
                & TRIM(ADJUSTL(cj_length))//') has wrong dimension length '// &
                & TRIM(ADJUSTL(ci_length))//' in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_1d_wrapper
  SUBROUTINE read_2d_wrapper(ifile_id,                 variable_name,        &
                           & alloc_array,              file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable
                                                    !< to be read
    REAL(wp), ALLOCATABLE, INTENT(out)  :: alloc_array(:,:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(2)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(2), cj_length(2)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_2D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  alloc_array=alloc_array                          )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(alloc_array,1)/=variable_dimls(1) .OR. &
         & SIZE(alloc_array,2)/=variable_dimls(2)) THEN
         WRITE(ci_length(1),*) SIZE(alloc_array,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(alloc_array,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_2d_wrapper
  SUBROUTINE read_3d_wrapper(ifile_id,                 variable_name,        &
                           & alloc_array,              file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable
                                                    !< to be read
    REAL(wp), ALLOCATABLE, INTENT(out)  :: alloc_array(:,:,:) !< values of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(3)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(3), cj_length(3)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_3D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  alloc_array=alloc_array                          )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(alloc_array,1)/=variable_dimls(1) .OR. &
         & SIZE(alloc_array,2)/=variable_dimls(2) .OR. &
         & SIZE(alloc_array,3)/=variable_dimls(3)) THEN
         WRITE(ci_length(1),*) SIZE(alloc_array,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(alloc_array,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         WRITE(ci_length(3),*) SIZE(alloc_array,3)
         WRITE(cj_length(3),*) variable_dimls(3)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//','//&
                & TRIM(ADJUSTL(cj_length(3)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//','//&
                & TRIM(ADJUSTL(ci_length(3)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
!!$         WRITE(0,*) TRIM(ADJUSTL(message1))
!!$         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_3d_wrapper

END MODULE mo_bc_aeropt_splumes
