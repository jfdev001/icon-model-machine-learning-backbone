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
!
!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"
!
! Description of *gscp_ice*:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, cloud ice, water vapor, rain and snow due to cloud microphysical
!   processes related to the formation of grid scale precipitation. This
!   includes the sedimentation of rain and snow. The precipitation fluxes at
!   the surface are also calculated here.
!
!   This is the extended version of gscp_cloudice with the option for
!   prognostic number density of cloud ice
!
!   The following parameterizations for heterogenous ice nucleation are
!   available depending of the integer parameter ice_nucleation.
!      0: Cooper formula (legacy option)
!      1: INAS ice nucleation with prescribed one-dimensional dust profile
!      2: INAS ice nucleation using prognostic ART dust (number and mean surface area)
!
! Method:
!   Prognostic bulk microphysical parameterization.
!   The sedimentation of ice, snow and rain is computed implicitly.
!
!
!------------------------------------------------------------------------------

MODULE gscp_ice

!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

!------------------------------------------------------------------------------
! Microphysical constants and variables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------

USE mo_kind, ONLY: wp, i4
USE mo_math_constants    , ONLY: pi, rad2deg
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                                 r_d   => rd    , & !> gas constant for dry air
                                 lh_v  => alv   , & !! latent heat of vapourization
                                 lh_s  => als   , & !! latent heat of sublimation
                                 cp_d  => cpd   , & !! spec. heat of dry air at constant press
                                 cpdr  => rcpd  , & !! (spec. heat of dry air at constant press)^-1
                                 cvdr  => rcvd  , & !! (spec. heat of dry air at const vol)^-1
                                 b3    => tmelt , & !! melting temperature of ice/snow
                                 t0    => tmelt , & !! melting temperature of ice/snow
                                 N_avo => avo,    & !! Avogadro number [1/mol]
                                 k_b   => ak,     & !! Boltzmann constant [J/K]
                                 grav               !! acceleration due to Earth's gravity

USE mo_lookup_tables_constants, ONLY:  &
                                 b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b4w   => c4les     !!               -- " --
USE mo_thdyn_functions,    ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
                                 sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
                                 latent_heat_vaporization, &
                                 latent_heat_sublimation
USE mo_exception,          ONLY: message, message_text, finish

!------------------------------------------------------------------------------

USE gscp_data, ONLY: &
    ccsrim,    ccsagg,    ccsdep,    ccsvel,    ccsvxp,    ccslam,       &
    ccslxp,    ccsaxp,    ccsdxp,    ccshi1,    ccdvtp,    ccidep,       &
    ccswxp,    zconst,    zcevxp,    zbevxp,                             &
    zvzxp,     zxstar,    zxcmin,    zami,                               &
    v0snow,                                                              &
    zvz0r => zvz0r0,      zbev => zbev0,        zcev => zcev0,           &
    x13o8,     x1o2,      x27o16,    x7o4,      x7o8,      x1o3,         &
    zbvi,      zcac,      zccau,     zciau,     zcicri,                  &
    zcrcri,    zcrfrz,    zcrfrz1,   zcrfrz2,   zeps,      zkcac,        &
    zkphi1,    zkphi2,    zkphi3,    zmi0,      zmimax,    zmsmin,       &
    zn0s0,     zn0s1,     zn0s2,     znimax_thom,          zqmin,        &
    zrho0,     zthet,     zthn,      ztmix,     ztrfrz,                  &
    zvz0i,     x2o3,      x5o24,     zams => zams_ci, zasmel,            &
    zbsmel,    zcsmel,    icesedi_exp,                                   &
    iautocon,  isnow_n0temp, dist_cldtop_ref,   reduce_dep_ref,          &
    tmin_iceautoconv,     zceff_fac, zceff_min,                          &
    mma, mmb, v_sedi_rain_min, v_sedi_snow_min, ice_coeffs, ice2mom

!==============================================================================

IMPLICIT NONE
PRIVATE

!------------------------------------------------------------------------------
!! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: cloudice2mom

LOGICAL, PARAMETER :: &
     ! removed lorig_icon because it was useless and confusing
     !
  lred_depgrowth = .TRUE., &  ! switch for reduced depositional growth near tops of stratus clouds
                              ! (combined with increased 'ztmix' parameter in order not to degrade T2M in Siberian winter)
  lsedi_ice    = .TRUE. ,  &  ! switch for sedimentation of cloud ice (Heymsfield & Donner 1990 *1/3)
  lstickeff    = .TRUE. ,  &  ! switch for sticking coeff. (work from Guenther Zaengl)
  licenum      = .TRUE. ,  &  ! switch for 2mom cloud ice (for .false. it should become the 1mom cloudice scheme)
  lice_hom     = .TRUE. ,  &  ! switch for homogeneous ice nucleation (for 2mom ice)
  lice_lat     = .TRUE. ,  &  ! switch for latitude dependency of vice and sticking efficiency
  lice_relax   = .FALSE. , &  ! switch for relaxation of depositional growth
  lice_qvel    = .FALSE. , &  ! switch for simple q-dependent ice fall speed in 2mom ice scheme
  lsuper_coolw = .TRUE.       ! switch for improved supercooled liquid water (work from Felix Rieper)

!------------------------------------------------------------------------------
!> Parameters and variables which are global in this module
!------------------------------------------------------------------------------

REAL(wp), PARAMETER, DIMENSION(1:5) :: &
     cdust = (/ 286.0_wp, 0.017_wp, 256.7_wp, 0.080_wp, 200.75_wp/)    ! dust of Ullrich et al. (2007)

REAL(wp), PARAMETER ::     & ! constant parameters for INAS-based ice nucleation scheme
     numdust = 1e+4_wp,    & ! number density of dust
     diadust = 5e-7_wp,    & ! diameter of dust 0.5 mu
     sigdust = 2.50_wp,    & ! standard deviation of lognormal dust distribution
     sfcdust = pi * EXP( 2.0_wp * LOG( sigdust )**2 ) * diadust**2  ! total surface area of dust

REAL(wp), PARAMETER ::     &
     zximin  = zami * 2.0e-6_wp**3, & ! Minimum mean mass of cloud ice ~1e-13
     zximax  = zami * 500e-6_wp**3, & ! Maximum mean mass of cloud ice ~1e-08
     zxstick = 1.0E-11_wp,          & ! minimal crystal mass for sticking efficiency (gscp3 only)
     zvnvq   = 0.70_wp,             & ! ratio of sedimentation velocities of mass and number of ice
     tau_ice = 7200.0_wp              ! relaxation timescale for ice nuclei (2 hours = 7200 s)

REAL(wp), PARAMETER ::             &  ! some constants needed for Kaercher and Lohmann parameterization
     ni_hom_max = 5000e3_wp      , &  ! number of liquid aerosols between 100-5000 per liter   [1/m3]
     rho_ice = 900.0             , &  ! density of solid ice
     r_0     = 0.25e-6           , &  ! aerosol particle radius prior to freezing
     alpha_d = 0.5               , &  ! deposition coefficient (KL02; Spichtinger & Gierens 2009)
     M_w     = 18.01528e-3       , &  ! molecular mass of water [kg/mol]
     M_a     = 28.96e-3          , &  ! molecular mass of air [kg/mol]
     ma_w    = M_w / N_avo       , &  ! mass of water molecule [kg]
     svol    = ma_w / rho_ice         ! specific volume of a water molecule in ice

REAL(wp), PARAMETER ::             &  ! some constants no longer provided by mo_math_constants
     pi4 = pi/4.0_wp

REAL(wp), PARAMETER ::             &
     tropics = 25.0_wp                ! for special treatment of tropics in case of lice_lat=.true.

!==============================================================================

CONTAINS

#ifdef _OPENACC
! GPU code can't flush to zero double precision denormals
! So to avoid CPU-GPU differences we'll do it manually
FUNCTION make_normalized(v)
  !$ACC ROUTINE SEQ
  REAL(wp) :: v, make_normalized

  IF (ABS(v) <= 2.225073858507201e-308_wp) THEN
    make_normalized = 0.0_wp
  ELSE
    make_normalized = v
  END IF
END FUNCTION
#else
FUNCTION make_normalized(v)
  REAL(wp) :: v, make_normalized
    make_normalized = v
END FUNCTION
#endif

!==============================================================================
!> Module procedure "cloudice2mom" in "gscp_ice" for computing effects of
!!  grid scale precipitation including cloud water, cloud ice, rain and snow
!------------------------------------------------------------------------------

SUBROUTINE cloudice2mom (            &
  nvec,ke,                           & !! array dimensions
  ivstart,ivend, kstart,             & !! optional start/end indicies
  idbg,                              & !! optional debug level
  zdt, dz,                           & !! time step and vertical grid spacing
  t,p,rho,qv,qc,qi,qr,qs,            & !! prognostic variables
  qnc,                               & !! diagnostic cloud droplet number
  dustnum,                           & !! prognostic dust concentration
  dustsfc,                           & !! prognostic dust surface area
  qni,ninact,                        & !! prognostic cloud ice number
  w,                                 & !! vertical velocity (for homogeneous nucleation)
  tropicsmask,                       & !! mask for tropics in [0,1]
  qi0,qc0,                           & !! cloud ice/water threshold for autoconversion
  prr_gsp,prs_gsp,pri_gsp,           & !! surface precipitation rates
  qrsflux,                           & !! total precipitation flux
  l_cv,                              & !! cv switch
  ldass_lhn,                         & !! lhn switch
  ithermo_water,                     & !! choice of water thermodynamics
  inucleation,                       & !! choice of ice nucleation
  ldiag_ttend,     ldiag_qtend     , & !! switches for optional diagnostic tendency output
  ddt_tend_t     , ddt_tend_qv     , & !! tendencies which are not
  ddt_tend_qc    , ddt_tend_qi     , & !! used anywhere
  ddt_tend_qr    , ddt_tend_qs       )

!------------------------------------------------------------------------------
! Description:
!   This module procedure calculates the rates of change of temperature, cloud
!   water, cloud ice, water vapor, rain and snow due to cloud
!   microphysical processes related to the formation of grid scale
!   precipitation.
!   The variables are updated in this subroutine. Rain and snow are
!   prognostic variables. The precipitation fluxes at the surface are stored
!   on the corresponding global fields.
!
! Method:
!   The sedimentation of ice, rain and snow is computed implicitly.
!
! Vectorization:
!   Most computations in this routine are grouped in IF-clauses. But the IFs
!   inside DO-loops often hinder or even disables vectorization.
!   For the big IF-chunks, the condition is now checked at the beginning of
!   the subroutine and the corresponding indices are stored in extra index
!   arrays. The following loops then are running only over these indices,
!   avoiding at least some IF-clauses inside the loops.
!
!------------------------------------------------------------------------------
!! Declarations:
!!
!------------------------------------------------------------------------------
!! Modules used: These are declared in the module declaration section
!! -------------

!! Subroutine arguments:
!! --------------------

  INTEGER, INTENT(IN) ::  &
    nvec          ,    & !> number of horizontal points
    ke                     !! number of grid points in vertical direction

  INTEGER, INTENT(IN), OPTIONAL ::  &
    ivstart   ,    & !> optional start index for horizontal direction
    ivend     ,    & !! optional end index   for horizontal direction
    kstart    ,    & !! optional start index for the vertical index
    idbg             !! optional debug level

  REAL(KIND=wp), INTENT(IN) :: &
    zdt             ,    & !> time step for integration of microphysics     (  s  )
    qi0,qc0                !> cloud ice/water threshold for autoconversion

  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::      &   ! (ie,ke)
    dz              ,    & !> layer thickness of full levels                (  m  )
    rho             ,    & !! density of moist air                          (kg/m3)
    p                      !! pressure                                      ( Pa  )

  REAL(KIND=wp), DIMENSION(:), INTENT(IN)   ::      &   ! dim (ie)
    tropicsmask            !! mask for tropics

  LOGICAL, INTENT(IN), OPTIONAL :: &
    l_cv, &                !! if true, cv is used instead of cp
    ldass_lhn

  INTEGER, INTENT(IN), OPTIONAL :: &
    inucleation,         & !! ice nucleation choice
    ithermo_water          !! water thermodynamics

  LOGICAL, INTENT(IN), OPTIONAL :: &
    ldiag_ttend,         & ! if true, temperature tendency shall be diagnosed
    ldiag_qtend            ! if true, moisture tendencies shall be diagnosed

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    t               ,    & !> temperature                                   (  K  )
    qv              ,    & !! specific water vapor content                  (kg/kg)
    qc              ,    & !! specific cloud water content                  (kg/kg)
    qi              ,    & !! specific cloud ice   content                  (kg/kg)
    qr              ,    & !! specific rain content                         (kg/kg)
    qs              ,    & !! specific snow content                         (kg/kg)
    qni             ,    & !! specific cloud ice number                     ( 1/kg)
    ninact          ,    & !! number of activated ice nuclei                ( 1/kg)
    w                      !! vertical velocity                             (m/s)

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL ::   &  ! dim (ie,ke)
    dustnum         ,    & !! dust concentration                            ( 1/kg)
    dustsfc                !! mean surface area of dust                     ( m2  )

  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
       qrsflux        ! total precipitation flux (nudg)

  REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) ::   &   ! dim (ie)
    prr_gsp,             & !> precipitation rate of rain, grid-scale        (kg/(m2*s))
    prs_gsp,             & !! precipitation rate of snow, grid-scale        (kg/(m2*s))
    pri_gsp,             & !! precipitation rate of cloud ice, grid-scale   (kg/(m2*s))
    qnc                    !! cloud number concentration

  REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), OPTIONAL ::   &     ! dim (ie,ke)
    ddt_tend_t      , & !> tendency T                                       ( 1/s )
    ddt_tend_qv     , & !! tendency qv                                      ( 1/s )
    ddt_tend_qc     , & !! tendency qc                                      ( 1/s )
    ddt_tend_qi     , & !! tendency qi                                      ( 1/s )
    ddt_tend_qr     , & !! tendency qr                                      ( 1/s )
    ddt_tend_qs         !! tendency qs                                      ( 1/s )

  !! Local parameters: None, parameters are in module header, gscp_data or data_constants
  !! ----------------

  !> Local scalars:
  !! -------------

  INTEGER (KIND=i4   ) ::  &
    iv, k             !> loop indices

  REAL    (KIND=wp   ) :: nnr

  REAL    (KIND=wp   ) :: z_heat_cap_r !! reciprocal of cpdr or cvdr (depending on l_cv)

  INTEGER ::  &
    iv_start     ,    & !> start index for horizontal direction
    iv_end       ,    & !! end index for horizontal direction
    k_start      ,    & !! model level where computations start
    izdebug             !! debug level

  REAL    (KIND=wp   ) ::  &
    fpvsw,             & ! name of statement function
    fxna ,             & ! statement function for ice crystal number
    fxna_cooper ,      & ! statement function for ice crystal number, Cooper(1986)
    ztx  ,             & ! dummy argument for statement functions
    znimax,            & ! maximum number of cloud ice crystals
    znimix,            & ! number of ice crystals at ztmix -> threshold temp for mixed-phase clouds
    zpvsw0,            & ! sat.vap. pressure at melting temperature
    zqvsw0,            & ! sat.specific humidity at melting temperature
    zdtr ,             & ! reciprocal of timestep for integration
    zscsum, zscmax, zcorr,  & ! terms for limiting  total cloud water depletion
    zsrsum, zsrmax,    & ! terms for limiting  total rain water depletion
    zssmax,            & ! terms for limiting snow depletion
    znin,              & ! number of cloud ice crystals at nucleation
    fnuc,              & !FR: coefficient needed for Forbes (2012) SLW layer parameterization
    znid,              & ! number of cloud ice crystals for deposition
    zmi ,              & ! mass of a cloud ice crystal
    zsvidep, zsvisub,  & ! deposition, sublimation of cloud ice
    zsimax , zsisum , zsvmax,& ! terms for limiting total cloud ice depletion
    zqvsw,             & ! sat. specitic humidity at ice and water saturation
    zqvsidiff,         & ! qv-zqvsi
    ztfrzdiff,         & ! ztrfrz-t
    zztau, zxfac, zx1, zx2, ztt, &   ! some help variables
    ztau, zphi, zhi, zdvtp, ztc, zeff, zlog_10

  REAL    (KIND=wp   ) ::  &
    zqct   ,& ! layer tendency of cloud water
    zqvt   ,& ! layer tendency of water vapour
    zqit   ,& ! layer tendency of cloud ice mass
    znit   ,& ! layer tendency of cloud ice number
    zqrt   ,& ! layer tendency of rain
    zqst      ! layer tendency of snow

  REAL    (KIND=wp   ) ::  &
    tg,ppg,rhog,qvg,qcg,qig,qrg,qsg,     & ! copies of the respective state variable for each cell
    nig,niactg,                          & ! copies of ice number variables
    temp_c,                              & ! temperature in deg. Cesius
    zlnqrk,zlnqsk,zlnlogmi,              &
    ccswxp_ln1o2,zvzxp_ln1o2,zbvi_ln1o2, &
    alf,bet,m2s,m3s,hlp,maxevap,zvi,taudepii

  LOGICAL :: &
    llqr,llqs,llqc,llqi  !   switch for existence of qr, qs, qc, qi

  LOGICAL :: lldiag_ttend, lldiag_qtend

  INTEGER :: ice_nucleation

  REAL(KIND=wp), DIMENSION(nvec,ke) ::   &
    t_in               ,    & !> temperature                                   (  K  )
    qv_in              ,    & !! specific water vapor content                  (kg/kg)
    qc_in              ,    & !! specific cloud water content                  (kg/kg)
    qi_in              ,    & !! specific cloud ice   content                  (kg/kg)
    qr_in              ,    & !! specific rain content                         (kg/kg)
    qs_in                     !! specific snow content                         (kg/kg)


!! Local (automatic) arrays:
!! -------------------------

  REAL    (KIND=wp   ) ::   &
    zqvsi             ,     & !> sat. specitic humidity at ice and water saturation
    zqvsw_up    (nvec),     & ! sat. specitic humidity at ice and water saturation in previous layer
    zvzr        (nvec),     & !
    zvzs        (nvec),     & !
    zvzi        (nvec),     & ! terminal fall velocity of ice mass
    zvzin       (nvec),     & ! terminal fall velocity of ice number
    zpkr        (nvec),     & ! precipitation flux of rain
    zpks        (nvec),     & ! precipitation flux of snow
    zpki        (nvec),     & ! precipitation flux of ice mass
    zpkin       (nvec),     & ! precipitation flux of ice number
    zprvr       (nvec),     & !
    zprvs       (nvec),     & !
    zprvi       (nvec),     & !
    zprvin      (nvec)


  REAL    (KIND=wp   ) ::   &
    zcsdep            ,     & !
    zcidep            ,     & !
    zvz0s             ,     & !
    zcrim             ,     & !
    zcagg             ,     & !
    zbsdep            ,     & !
    zcslam            ,     & !
    zn0s              ,     & !
    zimr              ,     & !
    zims              ,     & !
    zimi              ,     & ! as zimi but for ice number
    zini              ,     & !
    zzar              ,     & !
    zzas              ,     & !
    zzai              ,     & !
    zzni              ,     & ! as zzai but for ice number
    zqrk              ,     & ! rain water per volume of air (kg/m3)
    zqsk              ,     & ! snow water per volume of air (kg/m3)
    zqik              ,     & ! ice water per volume of air (kg/m3)
    znik              ,     & ! ice number per volume of air (1/m3)
    zdtdh             ,     & !
    zpsati            ,     & ! saturation pressure over ice
    zpsatw            ,     & ! saturation pressure over liquid water
    zssi              ,     & ! ice saturation ratio
    zinas             ,     & ! ice nuclating active site density
    zndust            ,     & ! number of dust particles
    zsdust            ,     & ! surface area of dust particles
    znhet             ,     & ! number of heterogeneous ice particles
    z1orhog           ,     & ! 1/rhog
    zrho1o2           ,     & ! (rho0/rhog)**1/2
    zrhofac_qi        ,     & ! (rho0/rhog)**icesedi_exp
    zeln7o8qrk        ,     & !
    zeln7o4qrk        ,     & ! FR new
    zeln27o16qrk      ,     & !
    zeln13o8qrk       ,     & !
    zeln5o24qsk       ,     & !
    zeln2o3qsk


  REAL    (KIND=wp   ) ::  &
    scau   , & ! transfer rate due to autoconversion of cloud water
    scac   , & ! transfer rate due to accretion of cloud water
    snuc   , & ! transfer rate due nucleation of cloud ice
    snucn  , & ! transfer rate due nucleation of cloud ice number
    shom   , & ! transfer rate due homogeneous nucleation of cloud ice
    shomn  , & ! transfer rate due homogeneous nucleation of cloud ice number
    scfrz  , & ! transfer rate due homogeneous freezing of cloud water
    scfrzn , & ! transfer rate due homogeneous freezing of cloud number
    simelt , & ! transfer rate due melting of cloud ice
    sidep  , & ! transfer rate due depositional growth of cloud ice
    ssdep  , & ! transfer rate due depositional growth of snow
    sdau   , & ! transfer rate due depositional cloud ice autoconversion
    srim   , & ! transfer rate due riming of snow
    sshed  , & ! transfer rate due shedding
    sicri  , & ! transfer rate due cloud ice collection by rain (sink qi)
    srcri  , & ! transfer rate due cloud ice collection by rain (sink qr)
    sagg   , & ! transfer rate due aggregation of snow and cloud ice
    siau   , & ! transfer rate due autoconversion of cloud ice (mass)
    siaun  , & ! transfer rate due autoconversion of cloud ice (number)
    ssmelt , & ! transfer rate due melting of snow
    sev    , & ! transfer rate due evaporation of rain
    srfrz  , & ! transfer rate due to rainwater freezing
    reduce_dep,&!FR: coefficient: reduce deposition at cloud top (Forbes 2012)
    dist_cldtop(nvec) !FR: distance from cloud top layer

  REAL(KIND=wp)  :: & ! for homogeneous ice nucleation
       v_th,n_sat,flux,phi,cool,tau,delta,scr,wcr,ctau,acoeff(3),bcoeff(2),ri_dot,  &
       kappa,sqrtkap,ren,R_imfc,R_im,R_ik,ri_0,zri,mi_hom,ni_hom,ri_hom,w_pre,zdi

#ifdef __LOOP_EXCHANGE
   REAL (KIND = wp )  ::  zlhv(ke), zlhs(ke)
#else
   REAL (KIND = wp )  ::  zlhv(nvec), zlhs(nvec) ! Latent heat if vaporization and sublimation
#endif

  LOGICAL :: lvariable_lh   ! Use constant latent heat (default .true.)

!------------ End of header ---------------------------------------------------

!------------------------------------------------------------------------------
!! Begin Subroutine cloudice
!------------------------------------------------------------------------------

!> Statement functions
! -------------------

! saturation vapour pressure over water (fpvsw), over ice (fpvsi)
! and specific humidity at vapour saturation (fqvs)
  fpvsw(ztx)     = b1*EXP( b2w*(ztx-b3)/(ztx-b4w) )

! Number of activate ice crystals;  ztx is temperature
  fxna(ztx)   = 1.0E2_wp * EXP(0.2_wp * (t0 - ztx))
  fxna_cooper(ztx) = 5.0E+0_wp * EXP(0.304_wp * (t0 - ztx))   ! FR: Cooper (1986) used by Greg Thompson(2008)

#ifdef _OPENACC
  CALL finish('mo_nwp_gscp_interface: ', 'subroutine cloudice2mom (gscp=3) not available on GPU') ! not tested
#endif

! Define reciprocal of heat capacity of dry air (at constant pressure vs at constant volume)

  IF (PRESENT(l_cv)) THEN
    IF (l_cv) THEN
      z_heat_cap_r = cvdr
    ELSE
      z_heat_cap_r = cpdr
    ENDIF
  ELSE
    z_heat_cap_r = cpdr
  ENDIF

  IF (PRESENT(ithermo_water)) THEN
     lvariable_lh = (ithermo_water .NE. 0)
  ELSE  ! Default themodynamic is constant latent heat
     lvariable_lh = .false.
  END IF

!------------------------------------------------------------------------------
!  Section 1: Initial setting of local and global variables
!------------------------------------------------------------------------------
  ! Input data
  !$ACC DATA &
  !$ACC   PRESENT(dz, t, p, rho, qv, qc, qi, qr, qs, qnc) &
  !$ACC   PRESENT(prr_gsp, prs_gsp, qrsflux, pri_gsp) &
  !$ACC   PRESENT(w, qni, ninact) &
  ! automatic arrays
  !$ACC   CREATE(zvzr, zvzs, zvzi, zvzin) &
  !$ACC   CREATE(zpkr, zpks, zpki, zpkin) &
  !$ACC   CREATE(zprvr, zprvs, zprvi, zprvin, zqvsw_up) &
  !$ACC   CREATE(dist_cldtop, zlhv, zlhs, acoeff, bcoeff)

! Some constant coefficients
  IF( lsuper_coolw) THEN
    znimax = znimax_Thom         !znimax_Thom = 250.E+3_wp,
    znimix = fxna_cooper(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  ELSE
    znimax = fxna(zthn) ! Maximum number of cloud ice crystals
    znimix = fxna(ztmix) ! number of ice crystals at temp threshold for mixed-phase clouds
  END IF

  zpvsw0 = fpvsw(t0)  ! sat. vap. pressure for t = t0
  zlog_10 = LOG(10._wp) ! logarithm of 10

  ! Precomputations for optimization
  ccswxp_ln1o2   = EXP (ccswxp * LOG (0.5_wp))
  zvzxp_ln1o2    = EXP (zvzxp * LOG (0.5_wp))
  zbvi_ln1o2     = EXP (zbvi * LOG (0.5_wp))

! Optional arguments

  IF (PRESENT(ivstart)) THEN
    iv_start = ivstart
  ELSE
    iv_start = 1
  END IF
  IF (PRESENT(ivend)) THEN
    iv_end = ivend
  ELSE
    iv_end = nvec
  END IF
  IF (PRESENT(kstart)) THEN
    k_start = kstart
  ELSE
    k_start = 1
  END IF
  IF (PRESENT(idbg)) THEN
    izdebug = idbg
  ELSE
    izdebug = 0
  END IF
  IF (PRESENT(ldiag_ttend)) THEN
    lldiag_ttend = ldiag_ttend
  ELSE
    lldiag_ttend = .FALSE.
  ENDIF
  IF (PRESENT(ldiag_qtend)) THEN
    lldiag_qtend = ldiag_qtend
  ELSE
    lldiag_qtend = .FALSE.
  ENDIF

  IF (.not.(PRESENT(inucleation).and.PRESENT(dustnum).and.PRESENT(dustsfc))) THEN
    ice_nucleation = 1
  ELSE
    ice_nucleation = inucleation
  ENDIF


  !$ACC DATA CREATE(t_in) IF(lldiag_ttend)
  !$ACC DATA CREATE(qv_in, qc_in, qi_in, qr_in, qs_in) IF(lldiag_qtend)

  ! save input arrays for final tendency calculation
  IF (lldiag_ttend) THEN
    !$ACC KERNELS DEFAULT(NONE) ASYNC(1)
    t_in  = t
    !$ACC END KERNELS
  ENDIF
  IF (lldiag_qtend) THEN
    !$ACC KERNELS DEFAULT(NONE) ASYNC(1)
    qv_in = qv
    qc_in = qc
    qi_in = qi
    qr_in = qr
    qs_in = qs
    !$ACC END KERNELS
  END IF

! timestep for calculations
  zdtr  = 1.0_wp / zdt

! output for various debug levels
  IF (izdebug > 10) THEN
    CALL message('gscp_ice: ','cloudice2mom')
  END IF
  IF (izdebug > 15) THEN
    WRITE (message_text,'(A,E10.3)') '      zams   = ',zams   ; CALL message('',message_text)
    WRITE (message_text,'(A,E10.3)') '      ccslam = ',ccslam ; CALL message('',message_text)
  END IF
  IF (izdebug > 20) THEN
    WRITE (message_text,*) '   nvec = ',nvec       ; CALL message('',message_text)
    WRITE (message_text,*) '   ke = ',ke           ; CALL message('',message_text)
    WRITE (message_text,*) '   ivstart = ',ivstart ; CALL message('',message_text)
    WRITE (message_text,*) '   ivend   = ',ivend   ; CALL message('',message_text)
  END IF
  IF (izdebug > 50) THEN
#if defined( _OPENACC )
    CALL message('cloudice2mom','GPU-info : update host before cloudice')
#endif
    !$ACC UPDATE HOST(dz, t, p, rho, qv, qc, qi, qr, qs) ASYNC(1)
    !$ACC WAIT(1)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN dz  = ',MAXVAL(dz),MINVAL(dz)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN T   = ',MAXVAL(t),MINVAL(t)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN p   = ',MAXVAL(p),MINVAL(p)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN rho = ',MAXVAL(rho),MINVAL(rho)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qv  = ',MAXVAL(qv),MINVAL(qv)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qc  = ',MAXVAL(qc),MINVAL(qc)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qr  = ',MAXVAL(qr),MINVAL(qr)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qi  = ',MAXVAL(qi),MINVAL(qi)
    CALL message('',message_text)
    WRITE (message_text,'(A,2E10.3)') '      MAX/MIN qs  = ',MAXVAL(qs),MINVAL(qs)
    CALL message('',message_text)
  ENDIF

  ! Delete precipitation fluxes from previous timestep
  !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(iv_start, iv_end)
  !$ACC LOOP GANG VECTOR
  DO iv = iv_start, iv_end
    prr_gsp (iv) = 0.0_wp
    prs_gsp (iv) = 0.0_wp
    zpkr(iv)     = 0.0_wp
    zpks(iv)     = 0.0_wp
    zpki(iv)     = 0.0_wp
    zprvr(iv)    = 0.0_wp
    zprvs(iv)    = 0.0_wp
    zprvi(iv)    = 0.0_wp
    zvzr(iv)     = 0.0_wp
    zvzs(iv)     = 0.0_wp
    zvzi(iv)     = 0.0_wp
    IF (licenum) THEN
      zpkin(iv)  = 0.0_wp
      zprvin(iv) = 0.0_wp
      zvzin(iv)  = 0.0_wp
    END IF
    dist_cldtop(iv) = 0.0_wp
    zqvsw_up(iv) = 0.0_wp
    pri_gsp (iv) = 0.0_wp
#ifndef __LOOP_EXCHANGE
    zlhv(iv)     = lh_v
    zlhs(iv)     = lh_s
#endif
  END DO
  !$ACC END PARALLEL

  ! Initialize latent heats to constant values
#ifdef __LOOP_EXCHANGE
  zlhv(:) = lh_v
  zlhs(:) = lh_s
#endif

! *********************************************************************
! Loop from the top of the model domain to the surface to calculate the
! transfer rates  and sedimentation terms
! *********************************************************************

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP SEQ
#ifdef __LOOP_EXCHANGE
  DO iv = iv_start, iv_end  !loop over horizontal domain

    ! Calculate latent heats if necessary
    IF ( lvariable_lh ) THEN
      DO  k = k_start, ke  ! loop over levels
        tg      = make_normalized(t(iv,k))
        zlhv(k) = latent_heat_vaporization(tg)
        zlhs(k) = latent_heat_sublimation(tg)
      END DO
    END IF

    DO  k = k_start, ke  ! loop over levels
#else
  DO  k = k_start, ke  ! loop over levels

    ! Calculate latent heats if necessary
    IF ( lvariable_lh ) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(tg)
      DO  iv = iv_start, iv_end  !loop over horizontal domain
        tg      = make_normalized(t(iv,k))
        zlhv(iv) = latent_heat_vaporization(tg)
        zlhs(iv) = latent_heat_sublimation(tg)
      END DO
    END IF

    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(alf, bet, fnuc, llqc, llqi, llqr) &
    !$ACC   PRIVATE(llqs, m2s, m3s, maxevap, nnr, ppg, qcg) &
    !$ACC   PRIVATE(rhog, qig, qrg, qsg, qvg, reduce_dep) &
    !$ACC   PRIVATE(sagg, scac, scau, scfrz, sdau, sev, siau, siaun) &
    !$ACC   PRIVATE(sicri, sidep, simelt, snuc, srcri, srfrz) &
    !$ACC   PRIVATE(srim, ssdep, sshed, ssmelt, temp_c) &
    !$ACC   PRIVATE(tg, z1orhog, zbsdep, zcagg, zcidep, zcorr) &
    !$ACC   PRIVATE(zcrim, zcsdep, zcslam, zdtdh, zdvtp, zeff) &
    !$ACC   PRIVATE(zeln13o8qrk, zeln27o16qrk, zeln5o24qsk) &
    !$ACC   PRIVATE(zeln2o3qsk, zeln7o4qrk, zeln7o8qrk) &
    !$ACC   PRIVATE(zhi, zimi, zimr, zims, hlp) &
    !$ACC   PRIVATE(zlnlogmi, zvi, zlnqrk, zlnqsk) &
    !$ACC   PRIVATE(zmi, zn0s, znid, znin, zphi, zqct) &
    !$ACC   PRIVATE(zqik, zqit, zqrk, zqrt, zqsk, zqst) &
    !$ACC   PRIVATE(zqvsi, zqvsidiff, zqvsw, zqvsw0) &
    !$ACC   PRIVATE(zqvt, zrho1o2, zrhofac_qi, zscmax, zscsum) &
    !$ACC   PRIVATE(zsimax, zsisum, zsrmax, zsrsum) &
    !$ACC   PRIVATE(zssmax, zsvidep, zsvisub, zsvmax) &
    !$ACC   PRIVATE(ztau, ztc, ztfrzdiff, ztt, zvz0s, zx1, zx2) &
    !$ACC   PRIVATE(zxfac, zzai, zzar, zzas, zztau) &
    !$ACC   PRIVATE(nig, niactg, znik, zzni, zini, znit, snucn, scfrzn)
    DO iv = iv_start, iv_end  !loop over horizontal domain
#endif

      ! add part of latent heating calculated in subroutine graupel to model latent
      ! heating field: subtract temperature from model latent heating field
      IF (ldass_lhn) THEN
        qrsflux(iv,k) = 0.0_wp
      ENDIF

      !----------------------------------------------------------------------------
      ! Section 2: Check for existence of rain and snow
      !            Initialize microphysics and sedimentation scheme
      !----------------------------------------------------------------------------

      zcrim  = 0.0_wp
      zcagg  = 0.0_wp
      zbsdep = 0.0_wp
      zvz0s  = 0.0_wp
      zn0s   = zn0s0
      reduce_dep = 1.0_wp  !FR: Reduction coeff. for dep. growth of rain and ice

      !----------------------------------------------------------------------------
      ! 2.1: Preparations for computations and to check the different conditions
      !----------------------------------------------------------------------------

      qrg  = make_normalized(qr(iv,k))
      qsg  = make_normalized(qs(iv,k))
      qvg  = make_normalized(qv(iv,k))
      qcg  = make_normalized(qc(iv,k))
      qig  = make_normalized(qi(iv,k))
      IF (licenum) THEN
        nig    = make_normalized(qni(iv,k))
        niactg = make_normalized(ninact(iv,k))
        zmi = MAX(MIN(qig/(nig+zeps),zximax),zximin)  ! mean crystal mass
      ELSE
        nig = 0.0_wp
        niactg = 0.0_wp
      END IF
      tg   = t(iv,k)
      ppg  = p(iv,k)
      rhog = rho(iv,k)

      !..for density correction of fall speeds
      z1orhog = 1.0_wp/rhog
      hlp     = LOG(zrho0*z1orhog)
      zrho1o2 = EXP(hlp*x1o2)
      zrhofac_qi = EXP(hlp*icesedi_exp)

      zqrk = qrg * rhog
      zqsk = qsg * rhog
      zqik = qig * rhog
      IF (licenum) THEN
        znik = nig * rhog
      ELSE
        znik = 0.0_wp
      END IF

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqi = zqik > zqmin

      zdtdh = 0.5_wp * zdt / dz(iv,k)

      zzar = zqrk/zdtdh + zprvr(iv) + zpkr(iv)
      zzas = zqsk/zdtdh + zprvs(iv) + zpks(iv)
      zzai = zqik/zdtdh + zprvi(iv) + zpki(iv)
      IF (licenum) THEN
        zzni = znik/zdtdh + zprvin(iv) + zpkin(iv)
      ELSE
        zzni = 0.0_wp
      END IF

      zpkr(iv) = 0.0_wp
      zpks(iv) = 0.0_wp
      zpki(iv) = 0.0_wp
      zpkin(iv) = 0.0_wp

      !-------------------------------------------------------------------------
      ! qs_prepare:
      !-------------------------------------------------------------------------
      IF (llqs) THEN
        IF (isnow_n0temp == 1) THEN
          ! Calculate n0s using the temperature-dependent
          ! formula of Field et al. (2005)
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
          zn0s = zn0s1*EXP(zn0s2*ztc)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSEIF (isnow_n0temp == 2) THEN
          ! Calculate n0s using the temperature-dependent moment
          ! relations of Field et al. (2005) who assume bms=2.0
          ztc = tg - t0
          ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
          nnr  = 3._wp
          hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
              + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
              + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
          alf = EXP(hlp*zlog_10) ! 10.0_wp**hlp
          bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
              + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
              + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3
          m2s = qsg * rhog / zams  ! assumes bms=2
          m3s = alf*EXP(bet*LOG(m2s))
          hlp  = zn0s1*EXP(zn0s2*ztc)
          zn0s = 13.50_wp * m2s * (m2s / m3s)**3
          zn0s = MAX(zn0s,hlp)
          zn0s = MIN(zn0s,1.0E2_wp*hlp)
          zn0s = MIN(zn0s,1.0E9_wp)
          zn0s = MAX(zn0s,1.0E6_wp)
        ELSE
          ! Old constant n0s
          zn0s = 8.0E5_wp
        ENDIF
        zcrim  = ccsrim*zn0s
        zcagg  = ccsagg*zn0s
        zbsdep = ccsdep*SQRT(v0snow)
        zvz0s  = ccsvel*EXP(ccsvxp * LOG(zn0s))
        zlnqsk = zvz0s * EXP (ccswxp * LOG (zqsk)) * zrho1o2
        ! Prevent terminal fall speed of snow from being zero at the surface level
        IF ( k == ke ) zlnqsk = MAX( zlnqsk, v_sedi_snow_min )
        zpks(iv) = zqsk * zlnqsk
        IF (zvzs(iv) == 0.0_wp) THEN
          zvzs(iv) = zlnqsk * ccswxp_ln1o2
        ENDIF
      ENDIF ! qs_prepare

      ! sedimentation fluxes

      !-------------------------------------------------------------------------
      ! qr_sedi:
      !-------------------------------------------------------------------------

      IF (llqr) THEN
        zlnqrk = zvz0r * EXP (zvzxp * LOG (zqrk)) * zrho1o2
        ! Prevent terminal fall speed of rain from being zero at the surface level
        IF ( k == ke ) zlnqrk = MAX( zlnqrk, v_sedi_rain_min )
        zpkr(iv) = zqrk * zlnqrk
        IF (zvzr(iv) == 0.0_wp) THEN
          zvzr(iv) = zlnqrk * zvzxp_ln1o2
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      ! qi_sedi:
      !-------------------------------------------------------------------------

      IF (llqi) THEN
        IF (.not.licenum .or. lice_qvel) THEN
          ! sedimentation of qi
          zvi = zvz0i * EXP (zbvi * LOG (zqik)) * zrhofac_qi
          zpki(iv) = zqik * zvi
          IF (zvzi(iv) == 0.0_wp) THEN
            zvzi(iv) = zvi * zbvi_ln1o2
          ENDIF
          IF (licenum) THEN
            ! sedimentation of ice number assuming vn=vq for cloud ice
            zpkin(iv) = znik * zvi * zvnvq
            IF (zvzin(iv) == 0.0_wp) THEN
              zvzin(iv) = zvi * zbvi_ln1o2 * zvnvq
            ENDIF
          ELSE
            zpkin(iv) = 0.0_wp
            zvzin(iv) = 0.0_wp
          ENDIF
        ELSE ! size-dependent fall speed of ice
          zvi = vice2mom(zqik,zmi,zrhofac_qi,tropicsmask(iv))
          zpki(iv)  = zqik * zvi
          zpkin(iv) = znik * zvi * zvnvq
          IF (zvzin(iv) == 0.0_wp) THEN
            zvzi(iv)  = zvi * zbvi_ln1o2
            zvzin(iv) = zvi * zbvi_ln1o2 * zvnvq
          ENDIF
        END IF
      ENDIF


      ! Prevent terminal fall speeds of precip hydrometeors from being zero at the surface level
      IF ( k == ke ) THEN
        zvzr(iv) = MAX( zvzr(iv), v_sedi_rain_min )
        zvzs(iv) = MAX( zvzs(iv), v_sedi_snow_min )
      ENDIF

      !--------------------------------------------------------------------------
      ! 2.3: Second part of preparations
      !--------------------------------------------------------------------------

      zeln7o8qrk    = 0.0_wp
      zeln7o4qrk    = 0.0_wp !FR
      zeln27o16qrk  = 0.0_wp
      zeln13o8qrk   = 0.0_wp
      zeln5o24qsk   = 0.0_wp
      zeln2o3qsk    = 0.0_wp
      zsrmax        = 0.0_wp
      zssmax        = 0.0_wp
      zndust        = 0.0_wp ! just to tell the compiler that
      zsdust        = 0.0_wp ! there is no loop dependency

      zcsdep        = 3.367E-2_wp
      zcidep        = 1.3E-5_wp
      zcslam        = 1e10_wp

      scau          = 0.0_wp
      scac          = 0.0_wp
      snuc          = 0.0_wp
      snucn         = 0.0_wp
      shom          = 0.0_wp
      shomn         = 0.0_wp
      scfrz         = 0.0_wp
      scfrzn        = 0.0_wp
      simelt        = 0.0_wp
      sidep         = 0.0_wp
      ssdep         = 0.0_wp
      sdau          = 0.0_wp
      srim          = 0.0_wp
      sshed         = 0.0_wp
      sicri         = 0.0_wp
      srcri         = 0.0_wp
      sagg          = 0.0_wp
      siau          = 0.0_wp
      siaun         = 0.0_wp
      ssmelt        = 0.0_wp
      sev           = 0.0_wp
      srfrz         = 0.0_wp

      zpsati = sat_pres_ice(tg)
      zpsatw = sat_pres_water(tg)

      hlp = 1._wp/(rhog * r_v * tg)
      zqvsi = zpsati * hlp
      zqvsw = zpsatw * hlp

      zssi = rhog * qvg * r_v * tg / zpsati

      zpkr(iv)   = MIN( zpkr(iv) , zzar )
      zpks(iv)   = MIN( zpks(iv) , zzas )
      zpki(iv)   = MIN( zpki(iv) , zzai )
      IF (licenum) THEN
        zpkin(iv) = MIN( zpkin(iv) , zzni )
      ELSE
        zpkin(iv) = 0.0_wp
      END IF

      zzar   = zdtdh * (zzar-zpkr(iv))
      zzas   = zdtdh * (zzas-zpks(iv))
      zzai   = zdtdh * (zzai-zpki(iv))
      IF (licenum) THEN
        zzni   = zdtdh * (zzni-zpkin(iv))
      ELSE
        zzni   = 0.0_wp
      END IF

      zimr   = 1.0_wp / (1.0_wp + zvzr(iv) * zdtdh)
      zims   = 1.0_wp / (1.0_wp + zvzs(iv) * zdtdh)
      zimi   = 1.0_wp / (1.0_wp + zvzi(iv) * zdtdh)
      IF (licenum) THEN
        zini   = 1.0_wp / (1.0_wp + zvzin(iv) * zdtdh)
      ELSE
        zini   = 0.0_wp
      END IF

      zqrk   = zzar*zimr
      zqsk   = zzas*zims
      zqik   = zzai*zimi
      IF (licenum) THEN
        znik   = zzni*zini
      ELSE
        znik   = 0.0_wp
      END IF

      llqr = zqrk > zqmin
      llqs = zqsk > zqmin
      llqi =  qig > zqmin
      llqc =  qcg > zqmin

      !!----------------------------------------------------------------------------
      !! 2.4: IF (llqr): ic1
      !!----------------------------------------------------------------------------

      IF (llqr) THEN
        zlnqrk   = LOG (zqrk)
        IF ( qig+qcg > zqmin ) THEN
          zeln7o8qrk   = EXP (x7o8   * zlnqrk)
        ENDIF
        IF ( tg < ztrfrz ) THEN
          zeln7o4qrk   = EXP (x7o4   * zlnqrk) !FR new
          zeln27o16qrk = EXP (x27o16 * zlnqrk)
        ENDIF
        IF (llqi) THEN
          zeln13o8qrk  = EXP (x13o8  * zlnqrk)
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.5: IF (llqs): ic2
      !!----------------------------------------------------------------------------

      IF (llqs) THEN

        zlnqsk   = LOG (zqsk)
        zeln5o24qsk  = EXP (x5o24  * zlnqsk)
        zeln2o3qsk   = EXP (x2o3   * zlnqsk)

      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.6: -- This section is only used in the graupel scheme
      !!----------------------------------------------------------------------------

      !!----------------------------------------------------------------------------
      !! 2.7:  slope of snow PSD and coefficients for depositional growth (llqi,llqs; ic3)
      !!----------------------------------------------------------------------------

      IF (llqi .OR. llqs) THEN
        zdvtp  = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg
        zhi    = ccshi1*zdvtp*rhog*zqvsi/(tg*tg)
        hlp    = zdvtp / (1.0_wp + zhi)
        zcidep = ccidep * hlp

        IF (llqs) THEN
          zcslam = EXP(ccslxp * LOG(ccslam * zn0s / zqsk ))
          zcslam = MIN(zcslam,1.0E15_wp)
          zcsdep = 4.0_wp * zn0s * hlp
        ENDIF
      ENDIF

      !!----------------------------------------------------------------------------
      !! 2.8: Deposition nucleation for low temperatures below a threshold (llqv)
      !!----------------------------------------------------------------------------

      IF ( licenum .and. ice_nucleation > 0 ) THEN
        ! INAS-based deposition nucleation with exponential vertical profile of dust
        zndust = numdust * MAX(MIN(exp(5e-3_wp*(ppg-300e2)),1e2_wp),1.0_wp)
        zsdust = sfcdust
        IF (ice_nucleation == 2 .and. dustnum(iv,k) > zndust) THEN
          zndust = dustnum(iv,k)
          zsdust = dustsfc(iv,k)
        END IF
        IF ( tg < zthet .and. qvg > 8.E-6_wp .and. qvg > 1.01_wp*zqvsi ) THEN
          zinas = het_icenuc_inas_depo(tg,zssi)
          znhet = zndust * (1.0_wp - EXP(-MAX(MIN(zinas*zsdust,30.0_wp),0.0_wp)))
          znin  = MAX(znhet - niactg, 0.0_wp)
          snucn = z1orhog * znin * zdtr
          snuc  = zmi0 * snucn
        ENDIF
      ELSE
        IF ( tg < zthet .AND. qvg >  8.E-6_wp .AND. qig <= 0.0_wp .AND. qvg > zqvsi ) THEN
          IF( lsuper_coolw ) THEN
            znin  = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin  = MIN( fxna(tg), znimax )
          END IF
          snucn = z1orhog * znin * zdtr
          snuc  = zmi0 * snucn
        ENDIF
      END IF

      !!----------------------------------------------------------------------------
      !! 2.8: Homogeneous ice nucleation for low temperatures
      !!----------------------------------------------------------------------------

      IF (licenum .and. lice_hom .and. tg < zthn) THEN

        ! critical supersaturation for homogeneous nucleation
        scr = 2.349 - tg * (1.0_wp/ 259.00_wp)

        !zssi = 1.2_wp * (zssi-1.0_wp) + 1.0_wp   ! increase supersaturation
        wcr = w(iv,k)  !* 4.0_wp                  ! and increase w

        IF (zssi > scr .AND. nig < ni_hom_max ) THEN

          zmi = MAX(MIN(qig/(nig+zeps),zximax),zximin)
          zri = EXP( (1./3.)*LOG(zmi/(4./3.*pi*rho_ice)) ) ! assumes spherical ice

          v_th  = SQRT( 8.*k_b*tg/(pi*ma_w) )
          flux  = alpha_d * v_th/4.
          n_sat = zpsati / (k_b*tg)
          zdvtp = ccdvtp * EXP(1.94_wp * LOG(tg)) / ppg

          ! coeffs of supersaturation equation
          acoeff(1) = (zlhs(iv) * grav) / (cp_d * r_v * tg**2) - grav/(r_d * tg)
          acoeff(2) = 1.0_wp/n_sat
          acoeff(3) = (zlhs(iv)**2 * M_w * ma_w)/(cp_d * ppg * tg * M_a)

          ! coeffs of depositional growth equation
          bcoeff(1) = flux * svol * n_sat * (zssi - 1.0_wp)
          bcoeff(2) = flux / zdvtp

          ! pre-existing ice crystals included as reduced updraft speed
          ri_dot = bcoeff(1) / (1. + bcoeff(2) * zri)
          R_ik   = (4 * pi) / svol * nig * zri**2 * ri_dot
          w_pre  = (acoeff(2) + acoeff(3) * zssi)/(acoeff(1) * zssi) * R_ik  ! KHL06 Eq. 19
          w_pre  = MAX(w_pre,0.0_wp)

          IF (wcr > w_pre) THEN   ! homogenous nucleation event

            zqvsidiff = qvg-zqvsi
            zsvmax    = zqvsidiff * zdtr

            ! timescales of freezing event (see KL02, RM05, KHL06)
            cool    = grav / cp_d * wcr
            ctau    = tg * ( 0.004_wp*tg - 2.0_wp ) + 304.4_wp
            tau     = 1.0_wp / (ctau * cool)                       ! freezing timescale, eq. (5)
            delta   = (bcoeff(2) * r_0)                            ! dimless aerosol radius, eq.(4)
            phi     = acoeff(1) * zssi / ( acoeff(2) + acoeff(3)*zssi) * (wcr - w_pre)

            ! monodisperse approximation following KHL06
            kappa   = 2.0_wp * bcoeff(1) * bcoeff(2) * tau / (1.+ delta)**2   ! kappa, Eq. 8 KHL06
            sqrtkap = SQRT(kappa)                                             ! root of kappa
            ren     = 3. * sqrtkap / ( 2. + SQRT(1.+9.*kappa/pi) )            ! analy. approx. of erfc by RM05
            R_imfc  = 4. * pi * bcoeff(1)/bcoeff(2)**2 / svol
            R_im    = R_imfc / (1.+ delta) * ( delta**2 - 1. &
                   & + (1.+0.5*kappa*(1.+ delta)**2) * ren/sqrtkap)           ! RIM Eq. 6 KHL06

            ! number concentration and radius of ice particles
            ni_hom  = phi / R_im                                              ! ni Eq.9 KHL06
            ri_0    = 1. + 0.5 * sqrtkap * ren                                ! for Eq. 3 KHL06
            ri_hom  = (ri_0 * (1. + delta) - 1. ) / bcoeff(2)                 ! Eq. 3 KHL06 * REN = Eq.23 KHL06
            mi_hom  = (4./3. * pi * rho_ice) * ri_hom**3
            mi_hom  = MIN(MAX(mi_hom,zximin),zximax)

            ! nucleation rate
            shomn = MIN(MAX(z1orhog*ni_hom, 0.d0),ni_hom_max) * zdtr
            shom  = MIN(mi_hom*shomn, zsvmax)

          END IF
        END IF
      END IF

      !!--------------------------------------------------------------------------
      !! Section 3: Search for cloudy grid points with cloud water and
      !!            calculation of the conversion rates involving qc (ic6)
      !!--------------------------------------------------------------------------

      IF (llqc) THEN

        zscmax = qcg*zdtr
        IF( tg > zthn ) THEN
          IF (iautocon == 0) THEN
            ! Kessler (1969) autoconversion rate
            scau = zccau * MAX( qcg - qc0, 0.0_wp )
            scac = zcac  * qcg * zeln7o8qrk
          ELSEIF (iautocon == 1) THEN
            ! Seifert and Beheng (2001) autoconversion rate
            ! with constant cloud droplet number concentration qnc
            IF (qcg > 1.0E-6_wp) THEN
              ztau = MIN(1.0_wp-qcg/(qcg+qrg),0.9_wp)
              ztau = MAX(ztau,1.E-30_wp)
              hlp  = EXP(zkphi2*LOG(ztau))
              zphi = zkphi1 * hlp * (1.0_wp - hlp)**3
              scau = zconst * qcg*qcg*qcg*qcg/(qnc(iv)*qnc(iv)) &
                   * (1.0_wp + zphi/(1.0_wp - ztau)**2)
              zphi = (ztau/(ztau+zkphi3))**4
              scac = zkcac * qcg * qrg * zphi
            ELSE
              scau = 0.0_wp
              scac = 0.0_wp
            ENDIF
          ENDIF
          IF (llqs) THEN
            srim = zcrim*  EXP(ccsaxp * LOG(zcslam)) * qcg
          ENDIF

          IF( tg >= t0 ) THEN
            sshed = srim
            srim  = 0.0_wp
          ENDIF
          ! Check for maximum depletion of cloud water and adjust the
          ! transfer rates accordingly
          zscsum = scau + scac + srim + sshed
          zcorr  = zscmax / MAX( zscmax, zscsum )
          scau   = zcorr*scau
          scac   = zcorr*scac
          srim   = zcorr*srim
          sshed  = zcorr*sshed

        ELSE !tg >= zthn: ! hom. freezing of cloud and rain water
          scfrz = zscmax
          IF (licenum) THEN
            scfrzn = scfrz / MAX(zxcmin,MIN(qcg/qnc(iv),zxstar))
          END IF
        ENDIF
        ! Calculation of heterogeneous nucleation of cloud ice.
        ! This is done in this section, because we require water saturation
        ! for this process (i.e. the existence of cloud water) to exist.
        ! Heterogeneous nucleation is assumed to occur only when no
        ! cloud ice is present and the temperature is below a nucleation
        ! threshold.
        IF( licenum .and. ice_nucleation > 0 ) THEN
          IF ( tg <= 267.15_wp .and. tg > 235.0_wp ) THEN
            zinas = EXP( 151.548_wp - 0.521_wp*tg )
            znhet = zndust * (1.0_wp - EXP(-MAX(MIN(zinas*zsdust,30.0_wp),0.0_wp)))
            znin  = MAX(znhet - niactg, 0.0_wp)
            snucn = z1orhog * znin * zdtr
            snuc  = zmi0 * snucn
          END IF
        ELSE IF ( tg <= 267.15_wp .AND. qig <= 0.0_wp ) THEN
          IF (lsuper_coolw ) THEN
            znin  = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin = MIN( fxna(tg), znimax )
          END IF
          snucn = z1orhog * znin * zdtr
          snuc  = zmi0 * snucn
        ENDIF

        ! Calculation of in-cloud rainwater freezing
        IF ( tg < ztrfrz .AND. qrg > 0.1_wp*qcg ) THEN
          IF (lsuper_coolw) THEN
            srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
          ELSE
            ztfrzdiff = ztrfrz-tg
            srfrz = zcrfrz*ztfrzdiff*SQRT(ztfrzdiff)* zeln27o16qrk
          ENDIF
        ENDIF

        ! Calculation of reduction of depositional growth at cloud top (Forbes 2012)
        IF( k>k_start .AND. k<ke .AND. lred_depgrowth ) THEN
          znin = MIN(fxna_cooper(tg), znimax )
          fnuc = MIN(znin/znimix, 1.0_wp)

          !! distance from cloud top
          IF( qv(iv,k-1) + qc(iv,k-1) < zqvsw_up(iv) .AND. qi(iv,k-1) + qs(iv,k-1) < zqmin ) THEN ! upper cloud layer
            dist_cldtop(iv) = 0.0_wp    ! reset distance to upper cloud layer
          ELSE
            dist_cldtop(iv) = dist_cldtop(iv) + dz(iv,k)
          END IF

          ! with asymptotic behaviour dz -> 0 (xxx)
          !        reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
          !                             dist_cldtop(iv)/dist_cldtop_ref + &
          !                             (1.0_wp-reduce_dep_ref)*(zdh/dist_cldtop_ref)**4), 1.0_wp)

          ! without asymptotic behaviour dz -> 0
          reduce_dep = MIN(fnuc + (1.0_wp-fnuc)*(reduce_dep_ref + &
                        dist_cldtop(iv)/dist_cldtop_ref), 1.0_wp)

        END IF ! Reduction of dep. growth of snow/ice

      ENDIF

      !------------------------------------------------------------------------
      ! Section 4: Search for cold grid points with cloud ice and/or snow and
      !            calculation of the conversion rates involving qi and qs
      !------------------------------------------------------------------------

      IF ( llqi .OR. llqs ) THEN
        llqs =  zqsk > zqmin !zqsk > zqmin
        llqi =   qig > zqmin

        IF (tg<=t0) THEN           ! cold case

          zqvsidiff = qvg-zqvsi
          zsvmax    = zqvsidiff * zdtr

          IF (licenum) THEN
            znin = MAX( MIN( rhog*nig, znimax ), 1e-6_wp )
          ELSEIF (lsuper_coolw) THEN
            znin = MIN( fxna_cooper(tg), znimax )
          ELSE
            znin = MIN( fxna(tg), znimax )
          END IF
          zmi  = MAX( MIN( rhog*qig/znin, zmimax ), zmi0 )
          IF (lstickeff.and.licenum) THEN
            zeff = effi2mom(tg,zmi,tropicsmask(iv))
          ELSEIF (lstickeff) THEN
            zeff = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
            zeff = MAX(zeff, zceff_min, zceff_fac*(tg-tmin_iceautoconv))  ! Guenther for gscp=1
          ELSE
            zeff = MIN(EXP(0.09_wp*(tg-t0)),1.0_wp)
            zeff = MAX(zeff,0.2_wp)                                       ! original version
          END IF
          IF (licenum) THEN
            IF (zqik > 1e-9_wp) THEN
              ! selfcollection of cloud ice
              zdi   = EXP( x1o3*LOG(zmi/zami) )
              zvi   = ice2mom%a_vel * EXP( ice2mom%b_vel*LOG(zmi) ) * zrhofac_qi
              hlp   = pi4 * nig * zdi * zdi * zeff * zdt
              siau  = hlp * ice_coeffs%sc_delta_q * qig                              &
                  & * SQRT( ice_coeffs%sc_theta_q * zvi * zvi + 2.0*ice2mom%s_vel**2 )
              siaun = hlp * ice_coeffs%sc_delta_n * nig                              &
                  & * SQRT( ice_coeffs%sc_theta_n * zvi * zvi + 2.0*ice2mom%s_vel**2 )
            END IF
          ELSE
            siau = zciau * MAX( qig - qi0, 0.0_wp ) * zeff
          END IF
          sagg = zcagg * EXP(ccsaxp*LOG(zcslam)) * qig * zeff
          znid = rhog * qig/zmi
          IF (llqi) THEN
            zlnlogmi  = LOG (zmi)
            IF (licenum .and. lice_relax) THEN
              taudepii  = zcidep * znid * EXP(0.33_wp * zlnlogmi)
              sidep     = zqvsidiff * (1.0_wp - EXP(-zdt*taudepii)) * zdtr
            ELSE
              sidep     = zcidep * znid * EXP(0.33_wp * zlnlogmi) * zqvsidiff
            END IF
          ELSE
            sidep = 0.0_wp
          ENDIF
          zsvidep   = 0.0_wp
          zsvisub   = 0.0_wp
          ! for sedimenting quantities the maximum
          ! allowed depletion is determined by the predictor value.
          IF (lsedi_ice) THEN
            zsimax  = zzai*z1orhog*zdtr
          ELSE
            zsimax  = qig*zdtr
          ENDIF
          IF( sidep > 0.0_wp ) THEN
            IF (lred_depgrowth ) THEN
              sidep = sidep * reduce_dep  !FR new: depositional growth reduction
            END IF
            zsvidep = MIN( sidep, zsvmax )
          ELSEIF ( sidep < 0.0_wp ) THEN
            IF (k < ke) THEN
              zsvisub  =   MAX (   sidep,  zsvmax)
              zsvisub  = - MAX ( zsvisub, -zsimax)
            ELSE
              zsvisub = - MAX(sidep, zsvmax )
              IF (zsvisub > zsimax) THEN
                ! Prefer reducing precipitation over reducing sublimation in
                ! lowest level to reduce time-step dependence. Factor 2 because
                ! zpki enters the final precipitation with a factor 1/2 (Crank-
                ! Nicholson).
                zpki(iv) = zpki(iv) - 2._wp * (zsvisub - zsimax)*rhog*dz(iv,k)
                zsvisub = zsvisub + 0.5 * MIN(0._wp, zpki(iv) / (rhog*dz(iv,k)))
                zpki(iv) = MAX(0._wp, zpki(iv))
                zzai = zsvisub/(z1orhog*zdtr)
                zqik = zzai*zimi
                zsimax   = zsvisub
              ENDIF
            ENDIF
          ENDIF

          IF (llqi) THEN
            zlnlogmi   = LOG  (zmsmin/zmi)
            zztau      = 1.5_wp*( EXP(0.66_wp*zlnlogmi) - 1.0_wp)
            sdau       = zsvidep/MAX(zztau,zeps)
          ELSE
            sdau    =  0.0_wp
          ENDIF

          sicri      = zcicri * qig * zeln7o8qrk
          ! Allow growth of snow only if the existing amount of snow is sufficiently large
          ! for a meaningful distiction between snow and cloud ice
          IF (qsg > 1.e-7_wp) srcri = zcrcri * (qig/zmi) * zeln13o8qrk

          zxfac = 1.0_wp + zbsdep * EXP(ccsdxp*LOG(zcslam))
          ssdep = zcsdep * zxfac * zqvsidiff / (zcslam+zeps)**2
          ! Check for maximal depletion of vapor by sdep
          IF (ssdep > 0.0_wp) THEN
            !FR new: depositional growth reduction
            IF (lred_depgrowth ) THEN
              ssdep = ssdep*reduce_dep
            END IF
            ! GZ: This limitation is crucial for numerical stability in the tropics!
            ssdep = MIN(ssdep, zsvmax-zsvidep)
          END IF


          ! Suppress depositional growth of snow if the existing amount is too small for a
          ! a meaningful distiction between cloud ice and snow
          IF (qsg <= 1.e-7_wp) ssdep = MIN(ssdep, 0.0_wp)

          ! Check for maximal depletion of cloud ice
          zsisum = siau + sdau + sagg + sicri + zsvisub
          zcorr  = 0.0_wp
          IF( zsimax > 0.0_wp ) zcorr  = zsimax / MAX( zsimax, zsisum )
          sidep  = zsvidep - zcorr*zsvisub
          sdau   = zcorr*sdau
          siau   = zcorr*siau
          sagg   = zcorr*sagg
          sicri  = zcorr*sicri
          IF (licenum) THEN
            siaun = zcorr*siaun
          END IF

        ELSE ! tg > 0 - warm case

          !------------------------------------------------------------------------
          ! Section 5: Search for warm grid points with cloud ice and/or snow and
          !            calculation of the melting rates of qi and ps
          !------------------------------------------------------------------------

          ! cloud ice melts instantaneously
          IF (lsedi_ice) THEN
            simelt = zzai*z1orhog*zdtr
          ELSE
            simelt = qig*zdtr
          ENDIF

          zqvsw0     = zpvsw0/(rhog * r_v * t0)
          zx1        = (tg - t0) + zasmel*(qvg - zqvsw0)
          zx2        = 1.0_wp + zbsmel * zeln5o24qsk
          ssmelt     = zcsmel * zx1 * zx2 * zeln2o3qsk
          ssmelt     = MAX( ssmelt, 0.0_wp )

        ENDIF ! tg

      ENDIF

      !--------------------------------------------------------------------------
      ! Section 6: Search for grid points with rain in subsaturated areas
      !            and calculation of the evaporation rate of rain
      !--------------------------------------------------------------------------

      IF( (llqr) .AND. (qvg+qcg <= zqvsw)) THEN

        zlnqrk   = LOG (zqrk)
        zx1      = 1.0_wp + zbev * EXP (zbevxp  * zlnqrk)
        ! Limit evaporation rate in order to avoid overshoots towards supersaturation
        ! the pre-factor approximates (esat(T_wb)-e)/(esat(T)-e) at temperatures between 0 degC and 30 degC
        temp_c = tg - t0
        maxevap     = (0.61_wp-0.0163_wp*temp_c+1.111e-4_wp*temp_c**2)*(zqvsw-qvg)/zdt
        sev    = MIN(zcev*zx1*(zqvsw - qvg) * EXP (zcevxp  * zlnqrk), maxevap)

        ! Calculation of below-cloud rainwater freezing
        IF ( tg < ztrfrz ) THEN
          IF (lsuper_coolw) THEN
            !FR new: reduced rain freezing rate
            srfrz = zcrfrz1*(EXP(zcrfrz2*(ztrfrz-tg))-1.0_wp ) * zeln7o4qrk
          ELSE
            srfrz = zcrfrz*SQRT( (ztrfrz-tg)**3 ) * zeln27o16qrk
          ENDIF
        ENDIF


      ENDIF

      !--------------------------------------------------------------------------
      ! Section 7: Calculate the total tendencies of the prognostic variables.
      !            Update the prognostic variables in the interior domain.
      !--------------------------------------------------------------------------
      zsrmax = zzar*z1orhog*zdtr

      zssmax   = zzas * z1orhog * zdtr

      zsrsum = sev + srfrz + srcri
      zcorr  = 1.0_wp
      IF(zsrsum > 0._wp) THEN
        zcorr  = zsrmax / MAX( zsrmax, zsrsum )
      ENDIF
      sev   = zcorr*sev
      srfrz = zcorr*srfrz
      srcri = zcorr*srcri
      ssmelt = MIN(ssmelt, zssmax)

      IF (ssdep < 0.0_wp ) THEN
        ssdep = MAX(ssdep, - zssmax)
      ENDIF

      zqvt =   sev    - sidep  - ssdep  - snuc   - shom
      zqct =   simelt - scau   - scfrz  - scac   - sshed  - srim
      zqit =   snuc   + shom   + scfrz  - simelt - sicri  + sidep  - sdau   - sagg   - siau
      zqrt =   scau   + sshed  + scac   + ssmelt - sev    - srcri  - srfrz
      zqst =   siau   + sdau   + sagg   - ssmelt + sicri  + srcri  + srim   + ssdep + srfrz
      IF (licenum) THEN
        znit = snucn + shomn + scfrzn - siaun - ( sdau + sagg + sicri + simelt ) / zmi
      ELSE
        znit = 0.0_wp
      END IF

#ifdef __LOOP_EXCHANGE
      ztt = z_heat_cap_r*( zlhv(k)*(zqct+zqrt) + zlhs(k)*(zqit+zqst) )
#else
      ztt = z_heat_cap_r*( zlhv(iv)*(zqct+zqrt) + zlhs(iv)*(zqit+zqst) )
#endif

      ! Update variables and add qi to qrs for water loading
      IF (lsedi_ice) THEN
        qig = MAX ( 0.0_wp, (zzai*z1orhog + zqit*zdt)*zimi)
      ELSE
        qig = MAX ( 0.0_wp, qig + zqit*zdt)
      END IF
      IF (licenum) THEN
        nig = MAX ( 0.0_wp, (zzni*z1orhog + znit*zdt)*zini)
      ELSE
        nig = 0.0_wp
      END IF
      qrg = MAX ( 0.0_wp, (zzar*z1orhog + zqrt*zdt)*zimr)
      qsg = MAX ( 0.0_wp, (zzas*z1orhog + zqst*zdt)*zims)

      ! time integration of ninact with relaxation term
      IF (licenum) THEN
        niactg = MAX ( 0.0_wp, niactg + ( snucn - niactg/tau_ice ) * zdt )
      ELSE
        niactg = 0.0_wp
      END IF

      !----------------------------------------------------------------------
      ! Section 8: Store satuaration specitic humidity at ice and water
      !            saturation of this layer
      !----------------------------------------------------------------------
      zqvsw_up(iv) = zqvsw ! to be available in the layer below (layer k-1)

      !----------------------------------------------------------------------
      ! Section 10: Complete time step
      !----------------------------------------------------------------------

      IF ( k /= ke) THEN
        ! Store precipitation fluxes and sedimentation velocities
        ! for the next level
        zprvr(iv) = qrg*rhog*zvzr(iv)
        zprvs(iv) = qsg*rhog*zvzs(iv)
        zprvi(iv) = qig*rhog*zvzi(iv)
        IF (licenum) THEN
          zprvin(iv) = nig*rhog*zvzin(iv)
        END IF

        IF (zprvr(iv) <= zqmin) zprvr(iv)=0.0_wp
        IF (zprvs(iv) <= zqmin) zprvs(iv)=0.0_wp
        IF (licenum) THEN
          IF (zprvi(iv) <= zqmin) THEN
            zprvi(iv)  = 0.0_wp
            zprvin(iv) = 0.0_wp
          END IF
        ELSE
          IF (zprvi(iv) <= zqmin) zprvi(iv)=0.0_wp
        END IF

        ! for the latent heat nudging
        IF (ldass_lhn) THEN
          IF (lsedi_ice) THEN
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)+zprvi(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv)+zpki(iv))
          ELSE
            qrsflux(iv,k) = zprvr(iv)+zprvs(iv)
            qrsflux(iv,k) = 0.5_wp*(qrsflux(iv,k)+zpkr(iv)+zpks(iv))
          END IF
        ENDIF

        IF (qrg+qr(iv,k+1) <= zqmin) THEN
          zvzr(iv)= 0.0_wp
        ELSE
          zvzr(iv)= zvz0r * EXP(zvzxp*LOG((qrg+qr(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (qsg+qs(iv,k+1) <= zqmin) THEN
          zvzs(iv)= 0.0_wp
        ELSE
          zvzs(iv)= zvz0s * EXP(ccswxp*LOG((qsg+qs(iv,k+1))*0.5_wp*rhog)) * zrho1o2
        ENDIF
        IF (licenum) THEN ! gscp=3
          IF (qig+qi(iv,k+1) <= zqmin ) THEN
            zvzi(iv)  = 0.0_wp
            zvzin(iv) = 0.0_wp
          ELSEIF (lice_qvel) THEN
            ! gscp=3 with simple q-dependent fall speed
            zvzi(iv)  = zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_wp*rhog)) * zrhofac_qi
            zvzin(iv) = zvzi(iv) * zvnvq
          ELSE
            ! gscp=3 with size-dependent fall speed
            zmi = MAX(MIN( (qig+qi(iv,k+1))/(nig+qni(iv,k+1)+zeps), zximax), zximin)
            zvi = vice2mom(0.5_wp*(qig+qi(iv,k+1))*rhog,zmi,zrhofac_qi,tropicsmask(iv))
            zvzi(iv)  = zvi
            zvzin(iv) = zvi * zvnvq
          END IF
        ELSE ! gscp=1
          IF (qig+qi(iv,k+1) <= zqmin ) THEN
            zvzi(iv)= 0.0_wp
          ELSE
            zvzi(iv)= zvz0i * EXP(zbvi*LOG((qig+qi(iv,k+1))*0.5_wp*rhog)) * zrhofac_qi
          ENDIF
        END IF

      ELSE
        ! Precipitation fluxes at the ground
        prr_gsp(iv) = 0.5_wp * (qrg*rhog*zvzr(iv) + zpkr(iv))
        IF (lsedi_ice) THEN
          prs_gsp(iv) = 0.5_wp * (rhog*qsg*zvzs(iv) + zpks(iv))
          pri_gsp(iv) = 0.5_wp * (rhog*qig*zvzi(iv) + zpki(iv))
        ELSE
          prs_gsp(iv) = 0.5_wp * (qsg*rhog*zvzs(iv) + zpks(iv))
        END IF

        ! for the latent heat nudging
        IF (ldass_lhn) THEN
          qrsflux(iv,k) = prr_gsp(iv)+prs_gsp(iv)
        ENDIF


      ENDIF

      ! Update of prognostic variables or tendencies
      qr (iv,k) = qrg
      qs (iv,k) = qsg
      qi (iv,k) = qig
      t  (iv,k) = t (iv,k) + ztt*zdt
      qv (iv,k) = MAX ( 0.0_wp, qv(iv,k) + zqvt*zdt )
      qc (iv,k) = MAX ( 0.0_wp, qc(iv,k) + zqct*zdt )
      IF (licenum) THEN
        qni (iv,k)   = MIN( nig, qig/zximin) ! the MIN is mostly for qi=0 -> qni=0
        ninact(iv,k) = niactg
      ELSE
        qni (iv,k)   = 0.0_wp
        ninact(iv,k) = 0.0_wp
      END IF

#if !defined (_OPENACC) && !defined (__SX__)
      IF (izdebug > 15) THEN
        ! Check for negative values
        IF (qr(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qr'
          CALL message('',message_text)
        ENDIF
        IF (qc(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qc'
          CALL message('',message_text)
        ENDIF
        IF (qi(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qi'
          CALL message('',message_text)
        ENDIF
        IF (qs(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qs'
          CALL message('',message_text)
        ENDIF
        IF (qv(iv,k) < 0.0_wp) THEN
          WRITE(message_text,'(A)') ' WARNING: cloudice, negative value in qv'
          CALL message('',message_text)
        ENDIF
      ENDIF
#endif

    END DO  !loop over iv

  END DO ! loop over levels
  !$ACC END PARALLEL

!------------------------------------------------------------------------------
! final tendency calculation for ICON
!
! Note: as soon as we have a new satad subroutine in ICON, this tendency
! calculation will be done in the k-loop and the original 3D variables wont
! be used to store the new values. Then we wont need the _in variables anymore.
!------------------------------------------------------------------------------

! calculated pseudo-tendencies

  IF ( lldiag_ttend ) THEN
    !$ACC DATA &
    !$ACC   PRESENT(ddt_tend_t, t, t_in)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, ke, iv_start, iv_end, zdtr)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_t (iv,k) = (t (iv,k) - t_in (iv,k))*zdtr
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA
  ENDIF

  IF ( lldiag_qtend ) THEN
    !$ACC DATA &
    !$ACC   PRESENT(ddt_tend_qv, ddt_tend_qc, ddt_tend_qr, ddt_tend_qs) &
    !$ACC   PRESENT(ddt_tend_qi, qv_in, qc_in, qr_in, qs_in, qi_in)

    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, ke, iv_start, iv_end, zdtr)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=k_start,ke
      DO iv=iv_start,iv_end
        ddt_tend_qv(iv,k) = MAX(-qv_in(iv,k)*zdtr,(qv(iv,k) - qv_in(iv,k))*zdtr)
        ddt_tend_qc(iv,k) = MAX(-qc_in(iv,k)*zdtr,(qc(iv,k) - qc_in(iv,k))*zdtr)
        ddt_tend_qi(iv,k) = MAX(-qi_in(iv,k)*zdtr,(qi(iv,k) - qi_in(iv,k))*zdtr)
        ddt_tend_qr(iv,k) = MAX(-qr_in(iv,k)*zdtr,(qr(iv,k) - qr_in(iv,k))*zdtr)
        ddt_tend_qs(iv,k) = MAX(-qs_in(iv,k)*zdtr,(qs(iv,k) - qs_in(iv,k))*zdtr)
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  ENDIF

  IF (izdebug > 15) THEN
#ifdef _OPENACC
   CALL message('cloudice2mom', 'GPU-info : update host after cloudice')
#endif
   !$ACC UPDATE HOST(t, qv, qc, qi, qr, qs) ASYNC(1)
   !$ACC WAIT(1)
   CALL message('cloudice2mom', 'UPDATED VARIABLES')
   WRITE(message_text,'(A,2E20.9)') 'cloudice  T= ',&
    MAXVAL( t(:,:)), MINVAL(t(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qv= ',&
    MAXVAL( qv(:,:)), MINVAL(qv(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qc= ',&
    MAXVAL( qc(:,:)), MINVAL(qc(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qi= ',&
    MAXVAL( qi(:,:)), MINVAL(qi(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qr= ',&
    MAXVAL( qr(:,:)), MINVAL(qr(:,:) )
    CALL message('', TRIM(message_text))
   WRITE(message_text,'(A,2E20.9)') 'cloudice qs= ',&
    MAXVAL( qs(:,:)), MINVAL(qs(:,:) )
   CALL message('', TRIM(message_text))
  ENDIF
  !$ACC WAIT(1)

  !$ACC END DATA ! IF(lldiag_qtend)
  !$ACC END DATA ! IF(lldiag_ttend)
  !$ACC END DATA ! general

!------------------------------------------------------------------------------
! End of subroutine cloudice
!------------------------------------------------------------------------------

END SUBROUTINE cloudice2mom

!==============================================================================

FUNCTION het_icenuc_inas_depo(tk,ssi) RESULT(inas)
  !$ACC ROUTINE SEQ
  !
  ! INAS-based deposition nucleation formula, see Ullrich et al. (2017, JAS)
  !
  REAL(wp), INTENT(in)  :: tk           !< temperature in K
  REAL(wp), INTENT(in)  :: ssi          !< supersaturation over ice
  REAL(wp)              :: inas         !< ice nucleating active site density in m^-2
  REAL(wp)              :: temp, acotan, tfunc, qfunc, pfunc

  LOGICAL,  PARAMETER :: loptimized = .false.   ! only very minor speedup
  REAL(wp), PARAMETER :: pcoeff(4) = (/-1.10099003e+07_wp, 1.36991120e+05_wp,-5.77301530e+02_wp, 8.60558248e-01_wp/)
  REAL(wp), PARAMETER :: qcoeff(3) = (/ 8.65974989e+06_wp,-8.45240144e+04_wp, 2.29723476e+02_wp/)

  temp   = MIN(MAX(tk,190.0_wp),260.0_wp)

  IF (loptimized) THEN
    ! rational function approximation with MAE smaller than 0.1 %
    ! but only valid for cdust(2:5) = (/ 0.017, 256.7, 0.080, 200.75/)
    pfunc = pcoeff(1) + pcoeff(2)*temp + pcoeff(3)*temp**2 + pcoeff(4)*temp**3
    qfunc = qcoeff(1) + qcoeff(2)*temp + qcoeff(3)*temp**2
    tfunc = pfunc/qfunc
  ELSE
    acotan = pi/2.0_wp - ATAN(cdust(4) * (temp - cdust(5)))
    tfunc  = COS( cdust(2)*(temp-cdust(3)) )**2 * acotan/pi
  END IF

  inas = EXP( cdust(1)*EXP(0.25_wp*LOG(ssi-1.0_wp)) * tfunc )
  inas = MAX(MIN(inas,1e15_wp),1e5_wp) ! paper recommends upper limit of 1e15

END FUNCTION het_icenuc_inas_depo

FUNCTION effi2mom(temp,zmi,ztropics) RESULT(zeff)
  !$ACC ROUTINE SEQ
  REAL(wp), INTENT(in)  :: temp, zmi, ztropics
  REAL(wp) :: zeff, xlat, zxi, zfac

  zeff = estick(temp) ! Connolly with zeff(0 C) = 0.14

  IF (lice_lat) THEN
    ! lower value of zceff_min outside of tropics
    xlat = 1.0_wp - ztropics ! [0,1]
    IF (xlat > 0.0_wp) THEN
      zeff = (1.0_wp - 0.5*xlat) * zeff
      zxi  = MIN(zmi/zximax,1.0_wp)
      zfac = MAX(EXP(0.2_wp*xlat * LOG(zxi)),0.1_wp)
      zeff = MERGE(zeff, zceff_min*zfac, temp > zthn)
    ELSE
      zeff = MERGE(zeff, zceff_min, temp > zthn)
    END IF
  ELSE
    ! apply the same zceff_min everywhere
    zeff = MERGE(zeff, zceff_min, temp > zthn)
  END IF

END FUNCTION effi2mom

FUNCTION estick (temp) RESULT(e_i)
  !$ACC ROUTINE SEQ
  REAL(wp), INTENT(in)  :: temp

  REAL(wp) :: e_i, T_c

  T_c = temp - t0

  ! piecewise linear sticking efficiency with maximum at -15 C,
  ! inspired by Figure 14 of Connolly et al. ACP 2012, doi:10.5194/acp-12-2055-2012
  ! Value at -40 C is based on Kajikawa and Heymsfield as cited by Philips et al. (2015, JAS)
  ! but not used in the cloudice2mom scheme.
  IF ( T_c >= 0_wp ) THEN
    e_i = 0.14_wp
  ELSEIF ( T_c >= -10_wp ) THEN
    e_i = -0.01_wp*(T_c+10_wp)+0.24_wp
  ELSEIF ( T_c >= -15_wp ) THEN
    e_i = -0.08_wp*(T_c+15_wp)+0.64_wp
  ELSEIF ( T_c >= -20_wp ) THEN
    e_i =  0.10_wp*(T_c+20_wp)+0.14_wp
  ELSEIF ( T_c >= -40_wp ) THEN
    e_i = 0.005_wp*(T_c+40_wp)+0.04_wp
  ELSE
    e_i = 0.04_wp
  END IF

END FUNCTION estick

FUNCTION vice2mom(zqi,zmi,zrhofac,ztropics) RESULT(zvi)
  !$ACC ROUTINE SEQ
  REAL(wp), INTENT(in)  :: zqi, zmi, zrhofac, ztropics
  REAL(wp) :: zvi, zxi, bvi, zvs, xlat
  LOGICAL,  PARAMETER :: l2mom_sedi = .false.
  REAL(wp), PARAMETER ::          &
     zxavi = 15.0_wp,             & ! prefactor in v=zxavi*xi**zbxvi
     zxbvi = 0.26_wp                ! exponent in v=zxavi*xi**zbxvi

  ! fall velocity of cloud ice uses a one-moment formula in tropics
  ! and an additional size dependency only outside of tropics. This is somewhat
  ! inconsistent, but using a two-moment formulation leads to temperature biases
  ! below the tropical tropopause (as large as 0.5 K bias at 200 hPa after 48 h).
  ! In the tropics the model is strongly constrained by the deep convection
  ! scheme and the microphysics has to do what the convection scheme requires.
  ! Introducing the full dependency in the tropics would require to reformulate
  ! or at least retune the deep convection.
  ! An alternative is to use the original one-moment sedimentation everywhere
  ! but this leads either to a large OLR bias in mid-latitudes or a large TEMP bias
  ! in the tropics, which is both not really acceptable.
  ! Currently it is not recommended to use the two-moment sedimentation.

  IF (l2mom_sedi) THEN

    ! two-moment sedimentation
    zvi = zxavi*EXP(zxbvi*LOG(zmi))

  ELSE

    ! one-moment sedimentation with option for size dependency outside of tropics
    zvi = zvz0i*EXP(zbvi*LOG(zqi))

    IF (lice_lat) THEN
      ! multiplicative size dependency only outside of tropics
      xlat = 1.0_wp - ztropics ! [0,1]
      IF (xlat > 0.0_wp) THEN
        zxi = MIN(zmi/zximax,1.0_wp)
        bvi = MIN(MAX(0.2*xlat,0.0_wp),0.2_wp)
        zvi = zvi * MAX(EXP(bvi * LOG(zxi)),0.2_wp)
      END IF
      ! Stokes asymptotic for small qi
      bvi = (1.0_wp - 0.9*xlat)
      zvs = (zqi*1e6_wp*bvi)**x2o3
      zvi = 1.0_wp/(1.0_wp/zvi+1.0_wp/zvs)   ! could be replaced by MIN(zvi,zvs) for efficiency
    ELSE
      zxi = MIN(zmi/zximax,1.0_wp)
      zvi = zvi * MAX(EXP(0.2_wp * LOG(zxi)),0.2_wp)
    END IF

  END IF

  ! density correction
  zvi = zvi * zrhofac

END FUNCTION vice2mom

!==============================================================================

END MODULE gscp_ice
