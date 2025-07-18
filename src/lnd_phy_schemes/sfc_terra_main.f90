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

! Soil Vegetation Atmosphere Transfer (SVAT) scheme TERRA
! "Nihil in TERRA sine causa fit." (Cicero)!!

MODULE sfc_terra_main

!===============================================================================
!
! @par Description:
!   The module "sfc_terra.f90" performs calculations related to the
!   parameterization of soil processes. It contains the subroutine terra, which
!   is the combination of the former parts terra1 and terra2 of the COSMO-Model.
!
!   All parametric scalar and array data for this soil model routine are
!   defined in the data module sfc_terra_data.f90.
!   All global fields that are used by the soil model are passed through the
!   argument list.
!   All global scalar variables of the model that are used by the soil model
!   routine terra are imported by USE statements below.
!
! @par reference   This is an adaptation of subroutine terra_multlay in file
!  src_soil_multlay.f90 of the lm_f90 library (COSMO code). Equation numbers refer to
!  Doms, Foerstner, Heise, Herzog, Raschendorfer, Schrodin, Reinhardt, Vogel
!    (September 2005): "A Description of the Nonhydrostatic Regional Model LM",

!------------------------------------------------------------------------------
!
! Modules used:
#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num
#endif

!------------------------------------------------------------------------------

USE sfc_terra_data  ! All variables from this data module are used by
                    ! this module. These variables start with letter "c"

!------------------------------------------------------------------------------


USE mo_mpi,                ONLY : get_my_global_mpi_id
!
USE mo_kind,               ONLY: wp
!
USE mo_physical_constants, ONLY: t0_melt => tmelt,& ! absolute zero for temperature
                                 r_d   => rd    , & ! gas constant for dry air
                                 rvd_m_o=>vtmpc1, & ! r_v/r_d - 1
                                 lh_v  => alv   , & ! latent heat of vapourization
                                 lh_s  => als   , & ! latent heat of sublimation
                                 lh_f  => alf   , & ! latent heat of fusion
                                 cp_d  => cpd   , & ! specific heat of dry air at constant press
                                 sigma => stbo  , & ! Boltzmann-constant
                                 rho_w => rhoh2o, & ! density of liquid water (kg/m^3)
                                 rdocp => rd_o_cpd  ! r_d / cp_d
!
USE mo_lnd_nwp_config,     ONLY: lmulti_snow, &
  &                              itype_trvg, &
  &                              itype_canopy, &
  &                              tau_skin, &
  &                              lterra_urb, &
  &                              itype_ahf
!
!
USE mo_run_config,         ONLY: msg_level
USE mo_fortran_tools,      ONLY: set_acc_host_or_device
USE mo_lnd_nwp_config,     ONLY: lcuda_graph_lnd

USE sfc_terra_evap,        ONLY: calc_evapotranspiration
USE sfc_terra_transport,   ONLY: calc_hydrology, calc_infiltration, calc_heat_conductivity, &
  &                              calc_heat_conduction, calc_soil_water_melt
USE sfc_terra_snow,        ONLY: snow_update_freshsnow_factor, &
  &                              snow_calc_precipitation_phase_change, snow_single_prepare, &
  &                              snow_single_soil_forcing, snow_single_melt, &
  &                              snow_single_calc_temperature, snow_single_update_new_state, &
  &                              snow_multi_prepare, snow_multi_handle_snowfall, &
  &                              snow_multi_soil_forcing, snow_multi_melt, &
  &                              snow_multi_update_new_state

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE

! Silence unused variable warnings on non-OpenACC builds
#ifdef _OPENACC
# define OPENACC_SUPPRESS_UNUSED_LZACC
#else
# define OPENACC_SUPPRESS_UNUSED_LZACC IF (lzacc .AND. acc_async_queue > 0) THEN; END IF
#endif

!------------------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------------------

PUBLIC :: terra

!------------------------------------------------------------------------------
! Public variables
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!+ Computation of the first part of the soil parameterization scheme
!------------------------------------------------------------------------------

  SUBROUTINE terra         (         &
                  nvec             , & ! array dimensions
                  ivstart          , & ! start index for computations in the parallel program
                  ivend            , & ! end index for computations in the parallel program
                  iblock           , & ! number of block
                  ke_soil, ke_snow , &
                  ke_soil_hy       , & ! number of active soil moisture layers
                  zmls             , & ! processing soil level structure
                  icant            , & ! canopy type
                  dt               , & ! time step
!
                  soiltyp_subs     , & ! type of the soil (keys 0-9)                     --
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m  )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  laifac           , & ! ratio between current LAI and laimax            --
                  eai              , & ! earth area (evaporative surface area) index     --
                  skinc            , & ! skin conductivity                        ( W/m**2/K )
! for TERRA_URB
                  urb_isa          , & ! urban impervious surface area                 (  -  )
                  urb_ai           , & ! surface area index of the urban canopy        (  -  )
                  urb_h_bld        , & ! building height                               (  m  )
                  urb_hcap         , & ! volumetric heat capacity of urban material (J/m**3/K)
                  urb_hcon         , & ! thermal conductivity of urban material        (W/m/K)
                  ahf              , & ! anthropogenic heat flux                      (W/m**2)
!
                  heatcond_fac     , & ! tuning factor for soil thermal conductivity
                  heatcap_fac      , & ! tuning factor for soil heat capacity
                  hydiffu_fac      , & ! tuning factor for hydraulic diffusivity
                  !
                  rsmin2d          , & ! minimum stomata resistance                    ( s/m )
                  r_bsmin          , & ! minimum bare soil evap resistance             ( s/m )
                  z0               , & ! vegetation roughness length                   ( m   )
!
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv               , & ! specific water vapor content                  (kg/kg)
                  qc               , & ! specific liquid-water content                 (kg/kg)
                  qi               , & ! specific frozen-water content                 (kg/kg)
                  ptot             , & ! full pressure                                 ( Pa  )
                  ps               , & ! surface pressure                              ( Pa  )
!
                  t_snow_now       , & ! temperature of the snow-surface               (  K  )
                  t_snow_new       , & ! temperature of the snow-surface               (  K  )
!
                  t_snow_mult_now  , & ! temperature of the snow-surface               (  K  )
                  t_snow_mult_new  , & ! temperature of the snow-surface               (  K  )
!
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_s_new          , & ! temperature of the ground surface             (  K  )
!
                  t_sk_now         , & ! skin temperature                              (  K  )
                  t_sk_new         , & ! skin temperature                              (  K  )
!
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s             , & ! specific humidity at the surface              (kg/kg)
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  w_snow_new       , & ! water content of snow                         (m H2O)
!
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  rho_snow_new     , & ! snow density                                  (kg/m**3)
!
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  rho_snow_mult_new, & ! snow density                                  (kg/m**3)
!
                  h_snow           , & ! snow depth                                   (  m  )
                  h_snow_gp        , & ! grid-point averaged snow depth               (  m  )
                  meltrate         , & ! snow melting rate                             (kg/(m**2*s))
                  tsnred           , & ! snow temperature offset for calculating evaporation  (K)
!
                  w_i_now          , & ! water content of interception water           (m H2O)
                  w_i_new          , & ! water content of interception water           (m H2O)
!
                  t_so_now         , & ! soil temperature (main level)                 (  K  )
                  t_so_new         , & ! soil temperature (main level)                 (  K  )
!
                  w_so_now         , & ! total soil water content (ice + liquid water) (m H20)
                  w_so_new         , & ! total soil water content (ice + liquid water) (m H20)
!
                  w_so_ice_now     , & ! soil ice content                              (m H20)
                  w_so_ice_new     , & ! soil ice content                              (m H20)
!
!                 t_2m             , & ! temperature in 2m                             (  K  )
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m            , & ! meridional wind in 10m                        ( m/s )
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow          , & ! snow-cover fraction                           (  -  )
!
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
!
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  wtot_snow_new    , & ! total (liquid + solid) water content of snow  (m H2O)
!
                  dzh_snow_now     , & ! layer thickness between half levels in snow   (  m  )
                  dzh_snow_new     , & ! layer thickness between half levels in snow   (  m  )
!
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  conv_frac        , & ! convective area fraction as assumed in convection scheme
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  pri_gsp          , & ! precipitation rate of ice, grid-scale        (kg/m2*s)
                  prg_gsp          , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
!
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  tfv              , & ! laminar reduction factor for evaporation      ( -- )
                  tfvsn            , & ! reduction factor for snow evaporation from model-DA coupling     ( -- )
!
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs             , & !!!! photosynthetic active radiation            ( W/m2)
!
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g         , & ! soil water runoff; sum over forecast          (kg/m2)
                  resid_wso        , & ! soil water budget, residuum                   (kg/m2)
!
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2)
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2)
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2)
                  zlhfl_snow       , & ! latent   heat flux snow/air interface         (W/m2)
                  lhfl_bs          , & ! latent heat flux from bare soil evap.         (W/m2)
                  lhfl_pl          , & ! latent heat flux from plants                  (W/m2)
                  plevap           , & ! function of accumulated plant evaporation     (kg/m2)
                  rstom            , & ! stomatal resistance                           ( s/m )
                  zshfl_sfc        , & ! sensible heat flux surface interface          (W/m2)
                  zlhfl_sfc        , & ! latent   heat flux surface interface          (W/m2)
                  zqhfl_sfc        , & ! moisture      flux surface interface          (kg/m2/s)
                  ldiff_qi         , & ! turbulent diffusion of frozen water is active
                  ldepo_qw         , & ! deposition of (frozen or liquid) cloud water required
                  lres_soilwatb    , & ! flag for computing the soil water budget
                  lacc             , & ! flag for activating OpenACC
                  opt_acc_async_queue) ! OpenACC stream


!-------------------------------------------------------------------------------
! Declarations
!-------------------------------------------------------------------------------


  INTEGER, INTENT(IN)  ::  &
                  icant,             & ! canopy type
                  nvec,              & ! array dimensions
                  ivstart,           & ! start index for computations in the parallel program
                  ivend,             & ! end index for computations in the parallel program
                  iblock,            & ! number of block
                  ke_soil, ke_snow,  &
                  ke_soil_hy           ! number of active soil moisture layers
  REAL    (KIND = wp), DIMENSION(ke_soil+1), INTENT(IN) :: &
                  zmls                 ! processing soil level structure
  REAL    (KIND = wp), INTENT(IN)  ::  &
                  dt                   ! time step

  INTEGER, DIMENSION(nvec), INTENT(IN) :: &
                  soiltyp_subs         ! type of the soil (keys 0-9)                     --

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
                  plcov            , & ! fraction of plant cover                         --
                  rootdp           , & ! depth of the roots                            ( m  )
                  sai              , & ! surface area index                              --
                  tai              , & ! transpiration area index                        --
                  laifac           , & ! ratio between current LAI and laimax
                  eai              , & ! earth area (evaporative surface area) index     --
                  skinc            , & ! skin conductivity                        ( W/m**2/K )
! for TERRA_URB
                  urb_isa          , & ! urban impervious surface area                 (  -  )
                  urb_ai           , & ! surface area index of the urban canopy        (  -  )
                  urb_h_bld        , & ! building height                               (  m  )
                  urb_hcap         , & ! volumetric heat capacity of urban material (J/m**3/K)
                  urb_hcon         , & ! thermal conductivity of urban material        (W/m/K)
                  ahf              , & ! anthropogenic heat flux from ISA        (W/m**2(isa))
!
                  heatcond_fac     , & ! tuning factor for soil thermal conductivity
                  heatcap_fac      , & ! tuning factor for soil heat capacity
                  hydiffu_fac      , & ! tuning factor for hydraulic diffusivity
                  rsmin2d          , & ! minimum stomata resistance                    ( s/m )
                  r_bsmin          , & ! minimum bare soil evap resistance             ( s/m )
                  u                , & ! zonal wind speed                              ( m/s )
                  v                , & ! meridional wind speed                         ( m/s )
                  t                , & ! temperature                                   (  k  )
                  qv               , & ! specific water vapor content                  (kg/kg)
                  qc               , & ! specific liquid-water content                 (kg/kg)
                  qi               , & ! specific frozen-water content                 (kg/kg)
                  ptot             , & ! full pressure                                 ( Pa )
                  ps               , & ! surface pressure                              ( pa  )
                  h_snow_gp        , & ! grid-point averaged snow depth
                  tsnred           , & ! snow temperature offset for calculating evaporation (K)
                  u_10m            , & ! zonal wind in 10m                             ( m/s )
                  v_10m            , & ! meridional wind in 10m                        ( m/s )
                  prr_con          , & ! precipitation rate of rain, convective        (kg/m2*s)
                  prs_con          , & ! precipitation rate of snow, convective        (kg/m2*s)
                  conv_frac        , & ! convective area fraction
                  prr_gsp          , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                  prs_gsp          , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                  pri_gsp          , & ! precipitation rate of ice, grid-scale         (kg/m2*s)
                  prg_gsp          , & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
                  sobs             , & ! solar radiation at the ground                 ( W/m2)
                  thbs             , & ! thermal radiation at the ground               ( W/m2)
                  pabs                 !!!! photosynthetic active radiation            ( W/m2)

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
                  z0                   ! vegetation roughness length                    ( m )

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(IN) :: &
                  tfv              , & ! laminar reduction factor for evaporation      ( -- )
                  tfvsn                ! reduction factor for snow evaporation from model-DA coupling     ( -- )

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
                  plevap               ! function of accumulated plant evaporation     (kg/m2)

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(INOUT) :: &
                  t_snow_now       , & ! temperature of the snow-surface (K)
                  t_s_now          , & ! temperature of the ground surface             (  K  )
                  t_sk_now         , & ! skin temperature                              (  K  )
                  t_g              , & ! weighted surface temperature                  (  K  )
                  qv_s             , & ! specific humidity at the surface              (kg/kg)
                  w_snow_now       , & ! water content of snow                         (m H2O)
                  rho_snow_now     , & ! snow density                                  (kg/m**3)
                  h_snow           , & ! snow depth
                  w_i_now          , & ! water content of interception store           (m H2O)
                  freshsnow        , & ! indicator for age of snow in top of snow layer(  -  )
                  zf_snow          , & ! snow-cover fraction
                  tch              , & ! turbulent transfer coefficient for heat       ( -- )
                  tcm              , & ! turbulent transfer coefficient for momentum   ( -- )
                  runoff_s         , & ! surface water runoff; sum over forecast       (kg/m2)
                  runoff_g         , & ! soil water runoff; sum over forecast          (kg/m2)
                  resid_wso            ! residuum of the budget of soil water content  (kg/m2)

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(OUT) :: &
                  t_snow_new       , & !
                  t_s_new          , & ! temperature of the ground surface             (  K  )
                  t_sk_new         , & ! skin temperature                              (  K  )
                  w_snow_new       , & ! water content of snow                         (m H2O)
                  rho_snow_new     , & ! snow density                                  (kg/m**3)
                  meltrate         , & ! snow melting rate
                  w_i_new          , & ! water content of interception store           (m H2O)
                  zshfl_s          , & ! sensible heat flux soil/air interface         (W/m2)
                  zlhfl_s          , & ! latent   heat flux soil/air interface         (W/m2)
                  zshfl_snow       , & ! sensible heat flux snow/air interface         (W/m2)
                  zlhfl_snow       , & ! latent   heat flux snow/air interface         (W/m2)
                  rstom            , & ! stomata resistance                            ( s/m )
                  lhfl_bs              ! latent heat flux from bare soil evap.         ( W/m2)


  REAL    (KIND = wp), DIMENSION(nvec,0:ke_snow), INTENT(INOUT) :: &
                  t_snow_mult_now      ! temperature of the snow-surface               (  K  )
  REAL    (KIND = wp), DIMENSION(nvec,0:ke_snow), INTENT(OUT) :: &
                  t_snow_mult_new      ! temperature of the snow-surface               (  K  )

  REAL    (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(INOUT) :: &
                  rho_snow_mult_now, & ! snow density                                  (kg/m**3)
                  wliq_snow_now    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_now    , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_now         ! layer thickness between half levels in snow   (  m  )

  REAL    (KIND = wp), DIMENSION(nvec,ke_snow), INTENT(OUT) :: &
                  rho_snow_mult_new, & ! snow density                                  (kg/m**3)
                  wliq_snow_new    , & ! liquid water content in the snow              (m H2O)
                  wtot_snow_new    , & ! total (liquid + solid) water content of snow  (m H2O)
                  dzh_snow_new         ! layer thickness between half levels in snow   (  m  )

  REAL    (KIND = wp), DIMENSION(nvec,0:ke_soil+1), INTENT(INOUT) :: &
                  t_so_now             ! soil temperature (main level)                 (  K  )
  REAL    (KIND = wp), DIMENSION(nvec,0:ke_soil+1), INTENT(OUT) :: &
                  t_so_new             ! soil temperature (main level)                 (  K  )

  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(INOUT) :: &
                  w_so_now         , & ! total soil water content (ice + liquid water) (m H20)
                  w_so_ice_now         ! soil ice content                              (m H20)
  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
                  w_so_new         , & ! total soil water content (ice + liquid water) (m H20)
                  w_so_ice_new         ! soil ice content                              (m H20)



!US why +1  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
!US is really +1 in ICON
  REAL    (KIND = wp), DIMENSION(nvec,ke_soil+1), INTENT(OUT) :: &
                  lhfl_pl          ! average latent heat flux from plants              ( W/m2)

  REAL    (KIND = wp), DIMENSION(nvec), INTENT(OUT) :: &
                  zshfl_sfc        , & ! sensible heat flux surface interface          (W/m2)
                  zlhfl_sfc        , & ! latent   heat flux surface interface          (W/m2)
                  zqhfl_sfc            ! latent   heat flux surface interface          (W/m2)

  LOGICAL, INTENT(IN) ::  &
                  ldiff_qi         , & ! turbulent diffusion of frozen water is active
                  ldepo_qw         , & ! deposition of (frozen or liquid) cloud water required
                  lres_soilwatb        ! calculation soil water budget desired

  LOGICAL, OPTIONAL, INTENT(IN) :: lacc
  INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

  LOGICAL :: lzacc
  INTEGER :: acc_async_queue

!--------------------------------------------------------------------------------
! TERRA Declarations
! ------------------

! Local parameters:
! ----------------

  REAL(wp), PARAMETER :: H_SNOW_GLAC_MIN = 1._wp !< Minimum snow height on glaciers [m].
  REAL(wp), PARAMETER :: H_SNOW_MAX = 40._wp !< Maximum snow height [m].

! Local scalars:
! -------------

  INTEGER        ::  &

    ! Indices
    kso            , & ! loop index for soil moisture layers
    ksn            , & ! loop index for snow layers
    k              , & ! loop index for snow layers
    i              , & ! loop index in x-direction
    mstyp              ! soil type index

  REAL(wp) :: tv_s !< virtual temperature at the surface [K].
  REAL(wp) :: radfl_th_dn !< downward longwave radiation [W/m^2].

  REAL(wp) :: organic_fraction !< utility variable

  REAL(wp) :: zalpha_uf !< height-dependent conductivity factor for buildings/street environment

! Local (automatic) arrays:
! -------------------------

  REAL(KIND = wp)                 :: &

    ! Two-time level variables exchange with interface
    h_snow_now     (nvec)          , & ! snow height  (m)
    h_snow_new     (nvec)          , & ! snow height  (m)

    ! Model geometry
    dz_snow_flx    (nvec)          , & ! snow depth for snow temperature gradient

    ! Multi-layer snow model
    zhh_snow       (nvec,  ke_snow), & ! depth of the half level snow layers
    zhm_snow       (nvec,  ke_snow), & ! depth of snow main levels
    zdzh_snow      (nvec,  ke_snow), & ! layer thickness between half levels
    zdzm_snow      (nvec,  ke_snow), & ! distance between main levels
    dt_t_snow_mult (nvec,0:ke_snow)    ! tendency of t_snow

  REAL(KIND=wp)                  ::  &

    ! Connection to the atmosphere
    dew_rate       (nvec)          , & ! dew formation rate
    rime_rate      (nvec)          , & ! rime formation rate
    rain_dew_rate  (nvec)          , & ! total rain rate including formation of dew and deposition of water droplets
    snow_rime_rate (nvec)          , & ! total snow rate including formation of rime and deposition of frozen particles
    eva_bs         (nvec)          , & ! evaporation from bare soil
    rho_ch         (nvec)          , & ! transfer coefficient*rho*g
    th_atm         (nvec)          , & ! potential temperature of lowest layer
    evapotrans_snfr(nvec)          , & ! total evapotranspiration
    evapo_snow     (nvec)          , & ! total evaporation from snow surface
    radfl_th_snow  (nvec)          , & ! thermal flux at snow surface
    forcing_soil   (nvec)          , & ! total forcing at soil surface
    hfl_snow_soil  (nvec)          , & ! heat-flux through snow
    zrnet_s        (nvec)          , & ! net radiation
    lhfl_precip    (nvec)          , & ! Latent heat flux due to melting/freezing precipitation [W/m^2].
    hfl_anthrop    (nvec)          , & ! TERRA_URB: anthropogenic heat flux [W/m^2(tile)]

    ! Tendencies
    dt_w_i         (nvec)          , & ! tendency of water content of interception store
    dt_w_snow      (nvec)          , & ! tendency of snow water content
    dt_t_s         (nvec)          , & ! tendency of t_s_now
    dt_t_snow      (nvec)          , & ! tendency of t_snow
    dt_w_so        (nvec,  ke_soil)    ! tendency of water content [kg/(m**3 s)]

  REAL    (KIND=wp) ::  &

    !   Interception variables
    w_i_max         (nvec)              ! maximum water content of interception store

  REAL    (KIND=wp) ::  &
    rad_flx         (nvec)          , & ! total radiative flux at surface             (W/m^2)
    rho_atm         (nvec)          , & ! air density of lowest atmospheric layer     (kg/m^3)

    !additional variables for root distribution
    zqhfl_s        (nvec)          , & ! moisture flux at soil/air interface
    zqhfl_snow     (nvec)              ! moisture flux at snow/air interface

  REAL    (KIND=wp) ::  &

    ! Soil and plant parameters
    hcap_ml     (nvec,ke_soil+1)   , & ! heat capacity

    ! Hydraulic variables
    transp_ml   (nvec,ke_soil)     , & ! transpiration contribution by the different layers
    transp_sum  (nvec)             , & ! total transpiration (transpiration from all
                                       !    soil layers)
    fr_w_ml     (nvec,ke_soil+1)   , & !fractional total water content of soil layers
    infil_rate  (nvec)             , & ! infiltration rate
    fr_liq_ml                      , & ! fractional liqu. water content of soil layer
    fr_ice_ml   (nvec,ke_soil+1)   , & ! fractional ice content of soil layer
    runoff_grav (nvec,ke_soil+1)       ! main level water gravitation
                                       !    flux (for runoff calc.)

  REAL    (KIND=wp) ::  &

    ! Thermal variables
    t_snow_top  (nvec)             , & ! snow surface temperaure

    ! HEATCOND (soil moisture dependent thermal conductivity of the soil)
    hcond_hl    (nvec,ke_soil)     , & ! thermal conductivity on half levels,
                                       ! meaning at the (lower) boundaries of
                                       ! soil layers 1 to ke_soil
    hcap_total                     , & ! total volumetric heat capacity of soil
    hcap_soil                      , & ! volumetric heat capacity of bare soil
    hcap_snow   (nvec)                 ! heat capacity of snow [J/m^2].

  REAL(KIND=wp)                  ::  &

    !   Multi-snow layer parameters
    fr_snow_lim   (nvec)           , &
    zalas_mult    (nvec,  ke_snow) , & ! heat conductivity of snow
    zextinct      (nvec,  ke_snow) , & ! solar radiation extinction coefficient in snow (1/m)
    zfor_snow_mult(nvec)               ! total forcing at snow surface

  REAL    (KIND=wp) ::  &
    ! Auxiliary variables
    hcond_ml     (nvec,ke_soil+1)  , & ! thermal conductivity on full levels, meaning at the centers of soil layers 1 to ke_soil+1
    dqvdt_snow   (nvec)            , & ! first derivative of saturation specific humidity
                                       !    with respect to t_snow
    rho_snow    (nvec)             , & ! snow density used for computing heat capacity and conductivity
    swe_correction                 , & ! SWE correction for ensuring height limits [m(H2O)].
    rho_snow_top                       ! Top-layer snow density [kg/m^3].

  REAL    (KIND=wp) ::  &
    budget_w_so_start(nvec)      ! utility variables for soil water budget

  REAL(wp) :: sp_10m(nvec) !< 10m mean wind speed [m/s].
  REAL(wp) :: rain_rate(nvec) !< Total rain rate (convective + grid-scale) [kg/(m^2 s)].
  REAL(wp) :: snow_rate(nvec) !< Total snow rate (convective + grid-scale) excluding ice [kg/(m^2 s)].
  REAL(wp) :: graupel_rate(nvec) !< Total graupel rate [kg/(m^2 s)].

!- End of header
!==============================================================================

  CALL set_acc_host_or_device(lzacc, lacc)

  IF(PRESENT(opt_acc_async_queue)) THEN
    acc_async_queue = opt_acc_async_queue
  ELSE
    acc_async_queue = 1
  ENDIF

!------------------------------------------------------------------------------
! Begin Subroutine terra
!------------------------------------------------------------------------------
!==============================================================================
!  Computation of the diagnostic part I of the soil parameterization scheme
!  In this part, evaporation from the surface and transpiration is calculated.
!  A multi-layer soil water distribution is calculated by simple bulk
!  parameterisation or a Penman-Monteith-type version of evaporation
!  and transpiration which can be used alternatively.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section I.1: Initializations
!------------------------------------------------------------------------------

  ! for the soil water budget
  !$ACC DATA PRESENT(resid_wso) CREATE(budget_w_so_start) ASYNC(acc_async_queue) IF(lres_soilwatb)

  ! Subroutine parameters IN
  !$ACC DATA ASYNC(acc_async_queue) &
  !$ACC   PRESENT(ivend) &
  !$ACC   PRESENT(zmls) &
  !$ACC   PRESENT(soiltyp_subs, plcov, rootdp, sai, eai, tai) &
  !$ACC   PRESENT(laifac) &
  !$ACC   PRESENT(skinc) &
  !$ACC   PRESENT(urb_isa, urb_ai, urb_h_bld) &
  !$ACC   PRESENT(urb_hcap, urb_hcon, ahf) &
  !$ACC   PRESENT(heatcond_fac, heatcap_fac, hydiffu_fac) &
  !$ACC   PRESENT(rsmin2d, r_bsmin, u, v, t, qv, qc, qi, ptot, ps, h_snow_gp, u_10m) &
  !$ACC   PRESENT(v_10m, prr_con, prs_con, conv_frac, prr_gsp, prs_gsp, pri_gsp) &
  !$ACC   PRESENT(prg_gsp, sobs, thbs, pabs, tsnred, z0) &

  ! Subroutine parameters INOUT
  !$ACC   PRESENT(t_snow_now, t_s_now, t_sk_now, t_g) &
  !$ACC   PRESENT(qv_s, w_snow_now) &
  !$ACC   PRESENT(rho_snow_now, h_snow, w_i_now) &
  !$ACC   PRESENT(freshsnow, zf_snow, tch, tcm, tfv, tfvsn, runoff_s) &
  !$ACC   PRESENT(runoff_g, t_snow_mult_now, rho_snow_mult_now) &
  !$ACC   PRESENT(wliq_snow_now, wtot_snow_now, dzh_snow_now) &
  !$ACC   PRESENT(t_so_now, w_so_now, w_so_ice_now, plevap) &

  ! Subroutine parameters OUT
  !$ACC   PRESENT(t_snow_new, t_s_new, t_sk_new) &
  !$ACC   PRESENT(w_snow_new, rho_snow_new) &
  !$ACC   PRESENT(meltrate, w_i_new, zshfl_s) &
  !$ACC   PRESENT(zlhfl_s, zshfl_snow, zlhfl_snow, rstom) &
  !$ACC   PRESENT(lhfl_bs, t_snow_mult_new, rho_snow_mult_new) &
  !$ACC   PRESENT(wliq_snow_new, wtot_snow_new, dzh_snow_new) &
  !$ACC   PRESENT(t_so_new, w_so_new, w_so_ice_new, lhfl_pl) &
  !$ACC   PRESENT(zshfl_sfc, zlhfl_sfc, zqhfl_sfc) &

  ! Local arrays
  !$ACC   CREATE(h_snow_now, h_snow_new) &
  !$ACC   CREATE(dz_snow_flx, zhh_snow, zhm_snow, zdzh_snow) &
  !$ACC   CREATE(zdzm_snow, dt_t_snow_mult) &
  !$ACC   CREATE(rain_dew_rate, snow_rime_rate, dew_rate, rime_rate, eva_bs) &
  !$ACC   CREATE(rho_ch, th_atm) &
  !$ACC   CREATE(evapotrans_snfr, evapo_snow, radfl_th_snow, forcing_soil, hfl_snow_soil) &
  !$ACC   CREATE(zrnet_s, lhfl_precip, hfl_anthrop, dt_w_i, dt_w_snow, dt_t_s) &
  !$ACC   CREATE(dt_t_snow, dt_w_so, w_i_max) &
  !$ACC   CREATE(rad_flx) &
  !$ACC   CREATE(rho_atm) &
  !$ACC   CREATE(zqhfl_s, zqhfl_snow, hcap_ml) &
  !$ACC   CREATE(transp_ml, transp_sum) &
  !$ACC   CREATE(fr_w_ml, infil_rate, fr_ice_ml) &
  !$ACC   CREATE(runoff_grav) &
  !$ACC   CREATE(t_snow_top) &
  !$ACC   CREATE(hcond_ml) &
  !$ACC   CREATE(hcap_snow) &
  !$ACC   CREATE(fr_snow_lim, zalas_mult) &
  !$ACC   CREATE(zextinct, zfor_snow_mult, hcond_hl, dqvdt_snow) &
  !$ACC   CREATE(rho_snow) &
  !$ACC   CREATE(sp_10m) &
  !$ACC   CREATE(rain_rate) &
  !$ACC   CREATE(snow_rate) &
  !$ACC   CREATE(graupel_rate) &

  !$ACC   NO_CREATE(budget_w_so_start) &

  ! Terra data module fields
  !$ACC   PRESENT(zzhls, zdzhs, zdzms)

  IF (lres_soilwatb) THEN
    ! Calculation of soil water budget: store recent runoff fluxes
    ! the final runoff fluxes are subtracted at the end
    CALL prepare_water_budget_diagnostic ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & runoff_s=runoff_s(:), &
        & runoff_g=runoff_g(:), &
        & budget_w_so_start=budget_w_so_start(:), & ! out
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  ENDIF

! Prepare basic surface properties (for land-points only)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp, tv_s, swe_correction, rho_snow_top)
  DO i = ivstart, ivend
    mstyp     = soiltyp_subs(i)        ! soil type

    IF (lmulti_snow) THEN
      rho_snow_top = rho_snow_mult_now(i,1)
    ELSE
      rho_snow_top = rho_snow_now(i)
    END IF

    ! ensure that glacier snow height is at between H_SNOW_GLAC_MIN and H_SNOW_MAX
    IF (mstyp == IST_ICE) THEN
      swe_correction = &
          ! snow water equivalent SWE (w_snow) increased to low snow height (h_snow).
          &   MAX(0._wp, H_SNOW_GLAC_MIN / rho_w * rho_snow_top - w_snow_now(i)) & !pos
          ! snow water equivalent SWE (w_snow) reduced   to top snow height (h_snow).
          & + MIN(0._wp, H_SNOW_MAX / rho_w * rho_snow_top - w_snow_now(i)) !neg
    ELSE
      ! snow water equivalent SWE (w_snow) reduced   to top snow height (h_snow).
      swe_correction = MIN(0._wp, H_SNOW_MAX / rho_w * rho_snow_top - w_snow_now(i)) !neg
    END IF

    ! Update SWE and height, corrected water goes to surface runoff
    w_snow_now(i) = w_snow_now(i) + swe_correction
    h_snow(i)     = h_snow(i)     + swe_correction * rho_w / rho_snow_top
    runoff_s(i)   = runoff_s(i)   - swe_correction * rho_w

    IF (lmulti_snow) THEN
      wtot_snow_now(i,1) = wtot_snow_now(i,1) + swe_correction
      dzh_snow_now(i,1)  = dzh_snow_now(i,1) + swe_correction * rho_w / rho_snow_top
    END IF

    IF (lres_soilwatb) THEN
      ! The budget diagnostic has the modified w_snow_now as reference point. Remove
      ! the created water from the initial sum to move that point to the actual initial value.
      budget_w_so_start(i) = budget_w_so_start(i) - swe_correction * rho_w
    END IF

#ifndef __SX__
    hcond_ml     (i,:) = cala0 (mstyp)              ! heat conductivity parameter
#endif

    meltrate(i) = 0.0_wp
    sp_10m(i) = SQRT(u_10m(i)**2 + v_10m(i)**2)

    ! Sum total rain, snow, and graupel rates.
    rain_rate(i)    = prr_con(i) + prr_gsp(i)
    snow_rate(i)    = prs_con(i) + prg_gsp(i) + prs_gsp(i)
    graupel_rate(i) = prs_con(i) + prg_gsp(i)

    rad_flx(i) = sobs(i)+thbs(i)
    tv_s = t_g (i)*(1.0_wp + rvd_m_o*qv_s(i))
    rho_atm(i) = ps(i)/(r_d*tv_s)

    ! moisture and potential temperature of lowest atmospheric layer
    th_atm(i) = t(i) * EXP(rdocp*LOG(ps(i)/ptot(i)))
  ENDDO

#ifdef __SX__
!$NEC outerloop_unroll(7)
  !$ACC LOOP SEQ
  DO kso = 1, ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp)
    DO i = ivstart, ivend
    mstyp       = soiltyp_subs(i)        ! soil type
    hcond_ml(i,kso)  = cala0(mstyp)              ! heat conductivity parameter
   ENDDO
  ENDDO
#endif

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    transp_sum(i)         = 0.0_wp
    fr_w_ml(i,ke_soil+1)  = w_so_now(i,ke_soil+1)/zdzhs(ke_soil+1)
    lhfl_bs(i)            = 0.0_wp
    lhfl_pl(i,ke_soil+1)  = 0.0_wp
    rstom  (i)            = 0.0_wp
    hfl_anthrop(i)        = 0.0_wp         ! TERRA_URB: Anthropogenic heat flux
  ENDDO

  ! REORDER
!$NEC outerloop_unroll(7)
  !$ACC LOOP SEQ
  DO kso   = 1, ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      fr_w_ml(i,kso) = w_so_now(i,kso)/zdzhs(kso)
      transp_ml(i,kso) = 0.0_wp
      lhfl_pl(i,kso) = 0.0_wp
    ENDDO
  ENDDO

!$NEC outerloop_unroll(8)
  !$ACC LOOP SEQ
  DO kso   = 1, ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ! No soil moisture for Ice and Rock
      IF (soiltyp_subs(i) <= IST_ROCK) THEN
        w_so_now(i,kso)         = 0._wp
        w_so_ice_now(i,kso)     = 0._wp
      END IF

      w_so_new(i,kso)         = w_so_now(i,kso)
      w_so_ice_new(i,kso)     = w_so_ice_now(i,kso)
    END DO
  END DO
  !$ACC END PARALLEL

  CALL calc_heat_conductivity ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & ke_soil=ke_soil, &
      & soiltyp_subs=soiltyp_subs(:), &
      & z_ml=zmls(:), &
      & z_hl=zzhls(:), &
      & dz_hl=zdzhs(:), &
      & t_so_now=t_so_now(:,0:), &
      & w_so_now=w_so_now(:,:), &
      & w_so_ice_now=w_so_ice_now(:,:), &
      & plcov=plcov(:), &
      & root_depth=rootdp(:), &
      & heatcond_fac=heatcond_fac(:), &
      & t_snred=tsnred(:), &
      & z0=z0(:), &
      & urb_h_bld=urb_h_bld(:), &
      & urb_ai=urb_ai(:), &
      & urb_isa=urb_isa(:), &
      & urb_hcon=urb_hcon(:), &
      & hcond_ml=hcond_ml(:,:), & ! out
      & hcond_hl=hcond_hl(:,:), & ! out
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

  CALL limit_transfer_coefficients ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & ke_soil=ke_soil, &
      & ke_snow=ke_snow, &
      & dt=dt, &
      & u_atm=u, &
      & v_atm=v, &
      & t_atm=t, &
      & th_atm=th_atm, &
      & rho_atm=rho_atm, &
      & t_g=t_g, &
      & qv_atm=qv, &
      & qv_s=qv_s, &
      & hcond_hl=hcond_hl, &
      & t_so_now=t_so_now, &
      & dz_ml=zdzms, &
      & dz_hl=zdzhs, &
      & rad_flx=rad_flx, &
      & w_snow_now=w_snow_now, &
      & rho_snow_now=rho_snow_now, &
      & t_snow_now=t_snow_now, &
      & t_snow_mult_now=t_snow_mult_now, &
      & tch=tch, & ! inout
      & rho_ch=rho_ch, & ! inout
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

  CALL snow_update_freshsnow_factor ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & dt=dt, &
      & w_snow=w_snow_now(:), &
      & t_snow=t_snow_now(:), &
      & h_snow_gp=h_snow_gp(:), &
      & t_atm=t(:), &
      & sp_10m=sp_10m(:), &
      & rain_rate=rain_rate(:), &
      & snow_rate=snow_rate(:), &
      & freshsnow=freshsnow(:), & ! inout
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

!------------------------------------------------------------------------------
! Section I.2: temperatures, water contents (in mH2O), surface pressure,
!------------------------------------------------------------------------------

  IF (lmulti_snow) THEN
    CALL snow_multi_prepare ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_snow=ke_snow, &
        & t_s_now=t_s_now(:), &
        & w_i_now=w_i_now(:), &
        & w_snow_now=w_snow_now(:), &
        & h_snow_new=h_snow_new(:), &
        & h_snow_now=h_snow_now(:), &
        & h_snow=h_snow(:), &
        & t_snow_mult_now=t_snow_mult_now(:,:), &
        & wtot_snow_now=wtot_snow_now(:,:), &
        & dzh_snow_now=dzh_snow_now(:,:), &
        & rho_snow_mult_now=rho_snow_mult_now(:,:), &
        & zhh_snow=zhh_snow(:,:), &
        & zhm_snow=zhm_snow(:,:), &
        & zdzh_snow=zdzh_snow(:,:), &
        & zextinct=zextinct(:,:), &
        & t_snow_top=t_snow_top(:), &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  ELSE
    CALL snow_single_prepare ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & w_snow_now=w_snow_now(:), &
        & t_snow_now=t_snow_now(:), &
        & rho_snow_now=rho_snow_now(:), &
        & rho_snow_mult_now_top=rho_snow_mult_now(:,1), &
        & t_s_now=t_s_now(:), &
        & freshsnow=freshsnow(:), &
        & fr_snow=zf_snow(:), &
        & h_snow_new=h_snow_new(:), &
        & h_snow_now=h_snow_now(:), &
        & dz_snow_flx=dz_snow_flx(:), &
        & fr_snow_lim=fr_snow_lim(:), &
        & rho_snow=rho_snow(:), &
        & hcap_snow=hcap_snow(:), &
        & t_snow_top=t_snow_top(:), &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  END IF

  ! INVARIANT: if itype_canopy == 1, then t_sk == t_s
  IF (itype_canopy == 1) THEN
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      t_sk_now(i) = t_s_now(i)
    ENDDO
    !$ACC END PARALLEL
  END IF

!------------------------------------------------------------------------------
! Section I.4: Hydrology, 1.Section
!------------------------------------------------------------------------------

  CALL calc_evapotranspiration ( &
      & dt=dt, &
      & icant=icant, &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & n_soil=ke_soil, &
      & n_soil_hy=ke_soil_hy, &
      & u_atm=u, &
      & v_atm=v, &
      & qv_atm=qv, &
      & t_atm=t, &
      & rho_atm=rho_atm, &
      & p_s=ps, &
      & t_s=t_s_now, &
      & t_sk=t_sk_now, &
      & t_snow_top=t_snow_top, &
      & fr_snow=zf_snow, &
      & w_snow=w_snow_now, &
      & t_snred=tsnred, &
      & w_i=w_i_now, &
      & w_so=w_so_now, &
      & w_so_ice=w_so_ice_now, &
      & fr_w_ml=fr_w_ml, &
      & tai=tai, &
      & eai=eai, &
      & sai=sai, &
      & plcov=plcov, &
      & laifac=laifac, &
      & rad_flx=rad_flx, &
      & par_absorbed=pabs, &
      & z_ml=zmls, &
      & z_hl=zzhls, &
      & dz_hl=zdzhs, &
      & soiltyp_subs=soiltyp_subs, &
      & root_depth=rootdp, &
      & r_bsmin=r_bsmin, &
      & r_stommin=rsmin2d, &
      & tcm=tcm, &
      & tch=tch, &
      & tfv=tfv, &
      & tfvsn=tfvsn, &
      & z0=z0, &
      & rho_ch=rho_ch, &
      & urb_isa=urb_isa, &
      & plevap=plevap, & ! inout
      & qv_s=qv_s, & ! out
      & eva_bs=eva_bs, & ! out
      & lhfl_bs=lhfl_bs, & ! out
      & transp_ml=transp_ml, & ! out
      & transp_sum=transp_sum, & ! out
      & lhfl_pl=lhfl_pl(:,1:ke_soil), & ! out
      & r_stom=rstom, & ! out
      & dqvdt_snow=dqvdt_snow, & ! out
      & eva_w_i=dt_w_i, & ! out
      & eva_w_sn=dt_w_snow, & ! out
      & dew_rate=dew_rate, & ! out
      & rime_rate=rime_rate, & ! out
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

!------------------------------------------------------------------------------
! End of former module procedure terra1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Former SUBROUTINE terra2
!------------------------------------------------------------------------------

!   In the prognostic part II the equation of heat conduction and water
!   transport is solved for a multi-layer soil using the same vertical grid
!   Freezing/melting of soil water/ice is accounted for (optionally). A
!   simple one-layer snow model provides the snow surface temperature and
!   the snow water equivalent.

!------------------------------------------------------------------------------
! Section II.1: Initializations
!------------------------------------------------------------------------------

  ! Initialisations
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  IF (lmulti_snow) THEN
    !$ACC LOOP SEQ
    DO ksn = 0, ke_snow
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        dt_t_snow_mult(i,ksn)  = 0.0_wp
        dt_t_s(i)     = 0.0_wp
      END DO
    END DO
  ELSE
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      dt_t_snow(i)  = 0.0_wp
      dt_t_s(i)     = 0.0_wp
    END DO
  END IF

!------------------------------------------------------------------------------
! Section II.2: Prepare basic surface properties and create some local
!               arrays of surface related quantities (for land-points only)
!               Initialise some fields
!------------------------------------------------------------------------------

  !$ACC LOOP SEQ
  DO   kso = 1,ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(organic_fraction, hcap_soil, hcap_total, fr_liq_ml)
    DO i = ivstart, ivend
      ! Scale soil heat capacity with organic fraction -> Chadburn et al., 2015
      organic_fraction = MAX(0._wp, plcov(i) * (rootdp(i) - zmls(kso)) / MAX(rootdp(i),eps_div))
      hcap_soil = crhoc(soiltyp_subs(i)) * heatcap_fac(i)
      hcap_total = (1.0_wp-organic_fraction) * hcap_soil + organic_fraction*0.58E+06_wp

      fr_ice_ml(i,kso) = w_so_ice_now(i,kso) / zdzhs(kso)   ! ice frac.
      fr_liq_ml        = fr_w_ml(i,kso) - fr_ice_ml(i,kso)  ! liquid water frac.
      hcap_ml(i,kso)   = hcap_total + rho_w * fr_liq_ml * chc_w + &
                              rho_w * fr_ice_ml(i,kso) * chc_i
      IF (kso <= ke_soil) THEN
        dt_w_so(i,kso) = 0.0_wp
      END IF
    END DO
  END DO      !soil layers

  IF (lterra_urb) THEN
    ! HW: modification of soil heat capacity in urban areas:
    !     interpolation between buildings and non-buildings,
    !     modification decreases with respect to the natural soil below.
    !     Below urb_h_bld, everything is natural soil.

    !$ACC LOOP SEQ
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zalpha_uf)
      DO i = ivstart, ivend

        zalpha_uf   = MAX(0.0_wp, MIN(zmls(kso)/urb_h_bld(i), 1.0_wp))

        hcap_ml(i,kso) =        urb_isa(i)  * ( (1.0_wp - zalpha_uf) * urb_hcap(i)*urb_ai(i)   &
                                                +         zalpha_uf  * hcap_ml(i,kso)           ) &
                    + (1.0_wp - urb_isa(i)) * hcap_ml(i,kso)

      ENDDO
    ENDDO
  END IF
  !$ACC END PARALLEL

!------------------------------------------------------------------------------
! Section II.3: Estimate thermal surface fluxes
!------------------------------------------------------------------------------
  !$NEC inline_complete
  CALL calc_infiltration ( &
      & dt=dt, &
      & nvec=nvec, &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & soiltyp_subs=soiltyp_subs, &
      & w_i_now=w_i_now, &
      & sp_10m=sp_10m, &
      & rain_rate=rain_rate, &
      & rain_rate_con=prr_con, &
      & snow_rate=snow_rate, &
      & conv_frac=conv_frac, &
      & ice_rate=pri_gsp, &
      & qc_atm=qc, &
      & qi_atm=qi, &
      & evapotrans_snfr=evapotrans_snfr, & ! out
      & evapo_snow=evapo_snow, & ! out
      & infil_rate=infil_rate, & ! out
      & runoff_s=runoff_s, & ! inout
      & rho_ch=rho_ch, &
      & fr_snow=zf_snow, &
      & w_snow_now=w_snow_now, &
      & eva_bs=eva_bs, &
      & transp_sum=transp_sum, &
      & plcov=plcov, &
      & tai=tai, &
      & urb_isa=urb_isa, &
      & dew_rate=dew_rate, &
      & rime_rate=rime_rate, &
      & rain_dew_rate=rain_dew_rate, & ! out
      & snow_rime_rate=snow_rime_rate, & ! out
      & dt_w_i=dt_w_i, & ! inout
      & dt_w_snow=dt_w_snow, & ! inout
      & t_sk_now=t_sk_now, &
      & t_s_now=t_s_now, &
      & w_i_max=w_i_max, & ! out
      & ldiff_qi=ldiff_qi, &
      & ldepo_qw=ldepo_qw, &
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

!------------------------------------------------------------------------------
! Section II.4: Soil water transport and runoff from soil layers
!------------------------------------------------------------------------------

  !$NEC inline_complete
  CALL calc_hydrology ( &
      & dt=dt, &
      & ke_soil=ke_soil, &
      & ke_soil_hy=ke_soil_hy, &
      & nvec=nvec, &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & soiltyp_subs=soiltyp_subs, &
      & dz_hl=zdzhs, &
      & dz_ml=zdzms, &
      & infil_rate=infil_rate, &
      & runoff_grav=runoff_grav, & ! out
      & w_so_new=w_so_new, & ! out
      & w_so_now=w_so_now, &
      & w_so_ice_now=w_so_ice_now, &
      & runoff_s=runoff_s, & ! inout
      & runoff_g=runoff_g, & ! inout
      & eva_bs=eva_bs, &
      & transp_ml=transp_ml, &
      & dt_w_so=dt_w_so, & ! out
      & root_depth=rootdp, &
      & z_ml=zmls, &
      & z_hl=zzhls, &
      & hydiffu_fac=hydiffu_fac, &
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

!------------------------------------------------------------------------------
! Section II.5: Soil surface heat flux (thermal forcing)
!------------------------------------------------------------------------------
  ! Estimate thermal surface fluxes:
  ! Estimate thermal surface fluxes over snow covered and snow free
  ! part of surface based on area mean values calculated in radiation
  ! code (positive = downward)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  ! TERRA_URB: Set anthropogenic heat flux
  IF (lterra_urb .AND. itype_ahf >= 1) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      hfl_anthrop(i) = urb_isa(i)*ahf(i)
    END DO
  END IF

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(radfl_th_dn)
  DO i = ivstart, ivend

    radfl_th_dn =   sigma*(1._wp - ctalb) * ( (1._wp - zf_snow(i))* &
                        t_sk_now(i) + zf_snow(i)*t_snow_top(i) )**4 + thbs(i)
    radfl_th_snow(i) = - sigma*(1._wp - ctalb)*t_snow_top(i)**4 + radfl_th_dn

    ! the estimation of the solar component would require the availability
    ! of the diffuse and direct components of the solar flux
    !
    ! Forcing for snow-free soil:
    ! (evaporation, transpiration, formation of dew and rime are already
    !  weighted by correspondind surface fraction)
    ! net radiation, sensible and latent heat flux

    zrnet_s(i) = sobs(i) - sigma*(1._wp - ctalb)*t_sk_now(i)**4 + radfl_th_dn
    zshfl_s(i) = cp_d*rho_ch(i) * (th_atm(i) - t_sk_now(i))
    zlhfl_s(i) = MERGE(lh_v, lh_s, t_sk_now(i) >= t0_melt) * evapotrans_snfr(i) &
                  / MAX(eps_div,(1._wp - zf_snow(i)))  ! take out (1-f) scaling
    zqhfl_s(i) = evapotrans_snfr(i)/ MAX(eps_div,(1._wp - zf_snow(i)))  ! take out (1-f) scaling
  END DO
  !$ACC END PARALLEL

  CALL snow_calc_precipitation_phase_change ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & dt=dt, &
      & rain_dew_rate=rain_dew_rate(:), &
      & snow_rime_rate=snow_rime_rate(:), &
      & w_snow_now=w_snow_now(:), &
      & t_snow_top=t_snow_top(:), &
      & soiltyp_subs=soiltyp_subs(:), &
      & t_so_now=t_so_now(:,:), &
      & hcap_ml=hcap_ml(:,:), &
      & dz_hl=zdzhs(:), &
      & w_i_now=w_i_now(:), &
      & w_i_max=w_i_max(:), &
      & fr_w_ml_top=fr_w_ml(:,1), &
      & dt_w_i=dt_w_i(:), & ! inout
      & dt_w_snow=dt_w_snow(:), & ! inout
      & dt_w_so_top=dt_w_so(:,1), & ! inout
      & runoff_s=runoff_s(:), & ! inout
      & lhfl_precip=lhfl_precip(:), & ! out
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

  IF (lmulti_snow) THEN
    CALL snow_multi_handle_snowfall ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_snow=ke_snow, &
        & dt=dt, &
        & t_s_now=t_s_now(:), &
        & th_atm=th_atm(:), &
        & rain_dew_rate=rain_dew_rate(:), &
        & snow_rime_rate=snow_rime_rate(:), &
        & w_snow_now=w_snow_now(:), &
        & dt_w_snow=dt_w_snow(:), &
        & h_snow_now=h_snow_now(:), & ! inout
        & wtot_snow_now=wtot_snow_now(:,:), & ! inout
        & wliq_snow_now=wliq_snow_now(:,:), & ! inout
        & rho_snow_mult_now=rho_snow_mult_now(:,:), & ! inout
        & t_snow_mult_now=t_snow_mult_now(:,:), & ! inout
        & zhm_snow=zhm_snow(:,:), & ! inout
        & zhh_snow=zhh_snow(:,:), & ! inout
        & zdzh_snow=zdzh_snow(:,:), & ! inout
        & zdzm_snow=zdzm_snow(:,:), & ! inout
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )

    CALL snow_multi_soil_forcing (&
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_snow=ke_snow, &
        & ke_soil=ke_soil, &
        & dt=dt, &
        & t_so_now=t_so_now(:,:), &
        & th_atm=th_atm(:), &
        & fr_snow=zf_snow(:), &
        & hcond_ml=hcond_ml(:,:), &
        & dz_ml=zdzms(:), &
        & rho_ch=rho_ch(:), &
        & sobs=sobs(:), &
        & radfl_th_snow=radfl_th_snow(:), &
        & radfl_net_snfr=zrnet_s(:), &
        & shfl_snfr=zshfl_s(:), &
        & lhfl_snfr=zlhfl_s(:), &
        & lhfl_precip=lhfl_precip(:), &
        & rain_dew_rate=rain_dew_rate(:), &
        & evapo_snow=evapo_snow(:), &
        & w_snow_now=w_snow_now(:), &
        & dt_w_snow=dt_w_snow(:), &
        & rho_snow_mult_now=rho_snow_mult_now(:,:), &
        & zalas_mult=zalas_mult(:,:), & ! out
        & t_snow_mult_now=t_snow_mult_now(:,:), &
        & zhm_snow=zhm_snow(:,:), &
        & zextinct=zextinct(:,:), &
        & hfl_snow_soil=hfl_snow_soil(:), & ! out
        & forcing_soil=forcing_soil(:), & ! out
        & shfl_snow=zshfl_snow(:), & ! out
        & lhfl_snow=zlhfl_snow(:), & ! out
        & qhfl_snow=zqhfl_snow(:), & ! out
        & zfor_snow_mult=zfor_snow_mult(:), & ! out
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  ELSE
    !NEC$ inline_complete
    CALL snow_single_soil_forcing ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & dt=dt, &
        & t_s_now=t_s_now(:), &
        & fr_snow=zf_snow(:), &
        & t_snow_top=t_snow_top(:), &
        & w_snow_now=w_snow_now(:), &
        & rho_snow=rho_snow(:), &
        & dz_snow_flx=dz_snow_flx(:), &
        & dt_w_snow=dt_w_snow(:), &
        & radfl_net_snfr=zrnet_s(:), &
        & shfl_snfr=zshfl_s(:), &
        & lhfl_snfr=zlhfl_s(:), &
        & lhfl_precip=lhfl_precip(:), &
        & hfl_anthrop=hfl_anthrop(:), &
        & hfl_snow_soil=hfl_snow_soil(:), & ! out
        & forcing_soil=forcing_soil(:), & ! out
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  END IF

!------------------------------------------------------------------------------
! Section II.6: Solution of the heat conduction equation, freezing/melting
!               of soil water/ice (optionally)
!------------------------------------------------------------------------------

  ! EM: If the single-layer snow model is used, nothing changes;
  ! if the multi-layer snow model is used and zf_snow(i) == 1,
  ! then the heat conduction equation for the whole column "soil + snow" is solved
  ! (see below, after the statement IF (lmulti_snow) THEN);
  ! if the multi-layer snow model is used, but zf_snow(i) < 1._wp,
  ! then the two partial temperature updates are computed: first, the heat conduction equation
  ! for the snow-free surface is solved, second, the heat conduction equation
  ! for the whole column "soil + snow" for snow-covered surface is solved.
  ! Then, the two updates are merged together.

  !$NEC inline_complete
  CALL calc_heat_conduction ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & ke_soil=ke_soil, &
      & ke_snow=ke_snow, &
      & dt=dt, &
      & z_ml=zmls(:), &
      & dz_ml=zdzms(:), &
      & dz_hl=zdzhs(:), &
      & fr_snow=zf_snow(:), &
      & t_s_now=t_s_now(:), &
      & hcond_ml=hcond_ml(:,:), &
      & hcond_hl=hcond_hl(:,:), &
      & hcap_ml=hcap_ml(:,:), &
      & forcing_soil=forcing_soil(:), &
      & w_snow_now=w_snow_now(:), &
      & w_snow_new=w_snow_new(:), & ! out
      & wliq_snow_now=wliq_snow_now(:,:), &
      & wtot_snow_now=wtot_snow_now(:,:), &
      & rho_snow_mult_now=rho_snow_mult_now(:,:), &
      & zalas_mult=zalas_mult(:,:), &
      & zhm_snow=zhm_snow(:,:), &
      & zdzh_snow=zdzh_snow(:,:), &
      & zdzm_snow=zdzm_snow(:,:), &
      & zfor_snow_mult=zfor_snow_mult(:), &
      & hfl_snow_soil=hfl_snow_soil(:), & ! out (ml snow)
      & t_snow_mult_now=t_snow_mult_now(:,0:), &
      & t_snow_mult_new=t_snow_mult_new(:,0:), & ! out (ml snow)
      & dt_t_snow_mult=dt_t_snow_mult(:,0:), & ! out (ml snow)
      & t_so_now=t_so_now(:,:), &
      & t_so_new=t_so_new(:,:), & ! out
      & dt_w_snow=dt_w_snow(:), & ! inout (ml snow)
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

  !NEC$ inline_complete
  CALL calc_soil_water_melt ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & ke_soil=ke_soil, &
      & dt=dt, &
      & soiltyp_subs=soiltyp_subs(:), &
      & z_ml=zmls(:), &
      & dz_hl=zdzhs(:), &
      & root_depth=rootdp(:), &
      & plcov=plcov(:), &
      & dt_w_so=dt_w_so(:,:), &
      & hcap_ml=hcap_ml(:,:), &
      & w_so_now=w_so_now(:,:), &
      & w_so_ice_now=w_so_ice_now(:,:), &
      & w_so_ice_new=w_so_ice_new(:,:), & ! out
      & t_so_now=t_so_now(:,0:), &
      & t_so_new=t_so_new(:,0:), & ! out
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
    )

!------------------------------------------------------------------------------
! Section II.7: Energy budget and temperature prediction at snow-surface
!------------------------------------------------------------------------------

  IF (.NOT. lmulti_snow) THEN
    !NEC$ inline_complete
    CALL snow_single_calc_temperature ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & dt=dt, &
        & t_so_new_top=t_so_new(:,1), &
        & t_s_now=t_s_now(:), &
        & w_snow_now=w_snow_now(:), &
        & dt_w_snow=dt_w_snow(:), &
        & dz_snow_flx=dz_snow_flx(:), &
        & fr_snow=zf_snow(:), &
        & rho_snow=rho_snow(:), &
        & sobs=sobs(:), &
        & radfl_th_snow=radfl_th_snow(:), &
        & rho_ch=rho_ch(:), &
        & hcap_snow=hcap_snow(:), &
        & th_atm=th_atm(:), &
        & evapo_snow=evapo_snow(:), &
        & hfl_snow_soil=hfl_snow_soil(:), &
        & dqvdt_snow=dqvdt_snow(:), &
        & t_snow_top=t_snow_top(:), &
        & t_snow_new=t_snow_new(:), & ! out
        & dt_t_snow=dt_t_snow(:), & ! out
        & shfl_snow=zshfl_snow(:), & ! out
        & lhfl_snow=zlhfl_snow(:), & ! out
        & qhfl_snow=zqhfl_snow(:), & ! out
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
        )
  END IF

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    ! next line has to be changed if a soil surface temperature is
    ! predicted by the heat conduction equation
    dt_t_s (i) = (t_so_new(i,1) - t_s_now(i)) / dt
  END DO
  !$ACC END PARALLEL

!------------------------------------------------------------------------------
! Section II.8: Melting of snow ,infiltration and surface runoff of snow water
!------------------------------------------------------------------------------

! If the soil surface temperature predicted by the equation of heat conduction
! is used instead of using T_s = T_so(1), the following section has to be
! adjusted accordingly.

  IF (lmulti_snow) THEN
    !NEC$ inline_complete
    CALL snow_multi_melt ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_snow=ke_snow, &
        & dt=dt, &
        & soiltyp_subs=soiltyp_subs(:), &
        & dz_top=zdzhs(1), &
        & w_snow_now=w_snow_now(:), &
        & zf_snow=zf_snow(:), &
        & dt_w_snow=dt_w_snow(:), & ! inout
        & dt_w_so_top=dt_w_so(:,1), & ! inout
        & fr_w_top=fr_w_ml(:,1), &
        & runoff_s=runoff_s(:), & ! inout
        & sobs=sobs(:), &
        & zextinct=zextinct(:,:), &
        & t_snow_mult_new=t_snow_mult_new(:,0:), &
        & dt_t_snow_mult=dt_t_snow_mult(:,0:), & ! inout
        & zdzh_snow=zdzh_snow(:,:), & ! inout
        & wtot_snow_now=wtot_snow_now(:,:), & ! inout
        & wliq_snow_now=wliq_snow_now(:,:), & ! inout
        & rho_snow_mult_now=rho_snow_mult_now(:,:), & ! inout
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  ELSE
    !NEC$ inline_complete
    CALL snow_single_melt ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & dt=dt, &
        & soiltyp_subs=soiltyp_subs(:), &
        & dz_top=zdzhs(1), &
        & hcap_ml_top=hcap_ml(:,1), &
        & hcap_snow=hcap_snow(:), &
        & fr_snow_lim=fr_snow_lim(:), &
        & w_snow_now=w_snow_now(:), &
        & dt_w_snow=dt_w_snow(:), & ! inout
        & t_snow=t_snow_new(:), &
        & dt_t_snow=dt_t_snow(:), & ! inout
        & rho_snow_now=rho_snow_now(:), &
        & t_s_now=t_s_now(:), &
        & t_so_top=t_so_new(:,1), &
        & dt_w_so_top=dt_w_so(:,1), & ! inout
        & dt_t_s=dt_t_s(:), & ! inout
        & fr_w_top=fr_w_ml(:,1), &
        & fr_ice_top=fr_ice_ml(:,1), &
        & meltrate=meltrate(:), & ! out
        & runoff_s=runoff_s(:), & ! inout
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  END IF

!------------------------------------------------------------------------------
! Section II.9: Final updating of prognostic values
!------------------------------------------------------------------------------

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    ! t_snow is computed above
    ! t_snow(i,nnew)  = t_snow(i,nx) + dt*dt_t_snow(i)
    t_so_new(i,1)   = t_so_now(i,1) + dt*dt_t_s(i)         ! (*)

    ! Next line has to be changed, if the soil surface temperature
    ! t_so(i,0,nnew) predicted by the heat conduction equation is used
    t_s_new   (i)   = t_so_new(i,1)
    t_so_new  (i,0) = t_so_new(i,1)
    w_i_new   (i)   = w_i_now(i) + dt*dt_w_i(i)/rho_w
  END DO
  !$ACC END PARALLEL

  IF (lmulti_snow) THEN
    !NEC$ inline_complete
    CALL snow_multi_update_new_state ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_snow=ke_snow, &
        & dt=dt, &
        & t_snow_new=t_snow_new(:), & ! out
        & w_snow_new=w_snow_new(:), & ! out
        & w_snow_now=w_snow_now(:), &
        & dt_w_snow=dt_w_snow(:), &
        & rho_snow_new=rho_snow_new(:), & ! out
        & h_snow_new=h_snow_new(:), & ! out
        & t_snow_mult_new=t_snow_mult_new(:,0:), & ! out
        & t_snow_mult_now=t_snow_mult_now(:,0:), &
        & dt_t_snow_mult=dt_t_snow_mult(:,0:), &
        & dzh_snow_new=dzh_snow_new(:,:), & ! out
        & zdzh_snow=zdzh_snow(:,:), &
        & wtot_snow_new=wtot_snow_new(:,:), & ! out
        & wtot_snow_now=wtot_snow_now(:,:), &
        & rho_snow_mult_new=rho_snow_mult_new(:,:), & ! out
        & rho_snow_mult_now=rho_snow_mult_now(:,:), &
        & wliq_snow_new=wliq_snow_new(:,:), & ! out
        & wliq_snow_now=wliq_snow_now(:,:), &
        & t_so_new_top=t_so_new(:,0), &
        & w_i_new=w_i_new(:), & ! inout
        & zhh_snow=zhh_snow(:,:), & ! out
        & zhm_snow=zhm_snow(:,:), & ! out
        & zdzm_snow=zdzm_snow(:,:), & ! out
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  ELSE
    CALL snow_single_update_new_state ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_snow=ke_snow, &
        & dt=dt, &
        & dz_hl_top=zdzhs(1), &
        & t_snow_new=t_snow_new(:), & ! out
        & t_snow_now=t_snow_now(:), &
        & dt_t_snow=dt_t_snow(:), &
        & w_snow_new=w_snow_new(:), & ! out
        & w_snow_now=w_snow_now(:), &
        & dt_w_snow=dt_w_snow(:), &
        & rho_snow_new=rho_snow_new(:), & ! out
        & rho_snow_now=rho_snow_now(:), &
        & rho_snow_mult_new=rho_snow_mult_new(:,:), & ! out
        & rho_snow_mult_now=rho_snow_mult_now(:,:), &
        & h_snow_new=h_snow_new(:), & ! out
        & h_snow_gp=h_snow_gp(:), &
        & t_so_new_top=t_so_new(:,0), &
        & w_i_new=w_i_new(:), & ! inout
        & w_i_max=w_i_max(:), &
        & sp_10m=sp_10m(:), &
        & th_atm=th_atm(:), &
        & fr_w_top=fr_w_ml(:,1), &
        & dt_w_so_top=dt_w_so(:,1), & ! inout
        & runoff_s=runoff_s(:), & ! inout
        & soiltyp_subs=soiltyp_subs(:), &
        & snow_rate=snow_rate(:), &
        & ice_rate=pri_gsp(:), &
        & graupel_rate=graupel_rate(:), &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  END IF

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  ! Update skin temperature
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (itype_canopy == 1) THEN
      t_sk_new(i)    = t_s_new(i)
    ELSE IF (itype_canopy == 2) THEN

      ! Calculation of the skin temperature (snow free area)
      ! by Schulz and Vogel (2020), based on Viterbo and Beljaars (1995).
      ! A Newtonian relaxation approach is used to ensure numerical stability.

      IF (w_snow_now(i) > eps_soil .OR. w_snow_new(i) > eps_soil) THEN
        t_sk_new(i) = t_s_new(i) ! needs to be t_s rather than t_snow in order to obtain correct t_g afterwards
      ELSE
        t_sk_new(i) = t_sk_now(i) + 0.5_wp*(t_s_new(i) - t_s_now(i)) + dt/tau_skin * &
          ( (zrnet_s(i) + zshfl_s(i) + zlhfl_s(i) + lhfl_precip(i) + hfl_anthrop(i)) / &
            MAX(5.0_wp,skinc(i)) - (t_sk_now(i) - t_s_now(i)) )
      ENDIF

    END IF
  END DO

  !$ACC LOOP SEQ
  DO kso = 1,ke_soil
  !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      w_so_new(i,kso) = w_so_now(i,kso) + dt*dt_w_so(i,kso)/rho_w
    END DO
  END DO        ! soil layers

  ! Update of two-time level interface variables
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    h_snow(i) = h_snow_new(i)
    IF (w_i_new(i) <= 1.0E-4_wp*eps_soil) THEN
      runoff_s(i) = runoff_s(i) + w_i_new(i) * rho_w
      w_i_new(i) = 0.0_wp
    END IF
  END DO

  ! computation of the weighted turbulent fluxes at the boundary surface-atmosphere
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    zshfl_sfc(i) = zshfl_s(i)*(1._wp - zf_snow(i)) + zshfl_snow(i)*zf_snow(i)

    ! Undo earlier division by MAX(eps_div, 1._wp - zf_snow(i)) or MAX(eps_div, zf_snow(i))
    zlhfl_sfc(i) = zlhfl_s(i)*MAX(eps_div, 1._wp - zf_snow(i)) + zlhfl_snow(i)*MAX(eps_div, zf_snow(i))
    zqhfl_sfc(i) = zqhfl_s(i)*MAX(eps_div, 1._wp - zf_snow(i)) + zqhfl_snow(i)*MAX(eps_div, zf_snow(i))
  END DO
  !$ACC END PARALLEL

  ! This block tests the residuum of water mass content in soil
  IF (lres_soilwatb) THEN
    CALL update_water_budget_diagnostic ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_soil_hy=ke_soil_hy, &
        & dt=dt, &
        & budget_w_so_start=budget_w_so_start(:), & ! inout
        & rain_rate=rain_rate(:), &
        & snow_rate=snow_rate(:), &
        & ice_rate=pri_gsp(:), &
        & qhfl_sfc=zqhfl_sfc(:), &
        & runoff_s=runoff_s(:), &
        & runoff_g=runoff_g(:), &
        & w_snow_new=w_snow_new(:), &
        & w_snow_now=w_snow_now(:), &
        & w_i_new=w_i_new(:), &
        & w_i_now=w_i_now(:), &
        & w_so_new=w_so_new(:,:), &
        & w_so_now=w_so_now(:,:), &
        & resid_wso=resid_wso(:), & ! inout
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )
  ENDIF

  IF (.NOT. lcuda_graph_lnd) THEN
    !$ACC WAIT(acc_async_queue)
  END IF

! for general fields
!$ACC END DATA

! for optional fields related to soil water budget
!$ACC END DATA

  IF (msg_level >= 19) THEN
    !$ACC UPDATE HOST(ivend) ASYNC(acc_async_queue)
    IF (.NOT. lcuda_graph_lnd) THEN
      !$ACC WAIT(acc_async_queue)
    END IF
    DO i = ivstart, ivend

     IF (ABS(t_s_now(i)-t_s_new(i)) > 15.0_wp .or. ABS(t_sk_now(i)-t_sk_new(i)) > 15.0_wp) THEN
        WRITE(*,'(A        )') '                                '
        WRITE(*,'(A,2I5)'  ) 'SFC-DIAGNOSIS terra output:  iblock = ', iblock, i

        WRITE(*,'(A        )') ' Temperatures and Humidities: '
        WRITE(*,'(A,2F28.16)') '   t_s      now/new :  ', t_s_now(i),       t_s_new(i)
        WRITE(*,'(A,2F28.16)') '   t_snow   now/new :  ', t_snow_now(i),    t_snow_new(i)
        WRITE(*,'(A, F28.16)') '   t_g              :  ', t_g(i)
        WRITE(*,'(A, F28.16)') '   qv_s (out)       :  ', qv_s(i)
        WRITE(*,'(A,2F28.16)') '   w_snow   now/new :  ', w_snow_now(i),    w_snow_new(i)
        WRITE(*,'(A,2F28.16)') '   rho_snow now/new :  ', rho_snow_now(i),  rho_snow_new(i)
        WRITE(*,'(A,2F28.16)') '   h_snow   now/new :  ', h_snow_now(i),    h_snow_new(i)
        WRITE(*,'(A, F28.16)') '   fresh_snow       :  ', freshsnow(i)
        WRITE(*,'(A, F28.16)') '   zf_snow (out)    :  ', zf_snow(i)
        WRITE(*,'(A,2F28.16)') '   w_i      now/new :  ', w_i_now(i),       w_i_new(i)
        DO k = 0, ke_soil+1
          WRITE(*,'(A,I1,A,2F28.16)') '   t_so    (',k,')      :  ', t_so_now (i,k), t_so_new (i,k)
        ENDDO
        DO k = 1, ke_soil+1
          WRITE(*,'(A,I1,A,2F28.16)') '   w_so    (',k,')      :  ', w_so_now (i,k), w_so_new (i,k)
        ENDDO
        DO k = 1, ke_soil+1
          WRITE(*,'(A,I1,A,2F28.16)') '   w_so_ice(',k,')      :  ', w_so_ice_now (i,k), w_so_ice_new (i,k)
        ENDDO
        WRITE(*,'(A        )') '                                '
        WRITE(*,'(A        )') ' Fluxes etc.            :  '
        WRITE(*,'(A, F28.16)') '   tcm              :  ', tcm(i)
        WRITE(*,'(A, F28.16)') '   tch              :  ', tch(i)
        WRITE(*,'(A, F28.16)') '   runoff_s (out)       :  ', runoff_s(i)
        WRITE(*,'(A, F28.16)') '   runoff_g (out)       :  ', runoff_g(i)
        WRITE(*,'(A, F28.16)') '   rstom            : ', rstom(i)
        IF (itype_trvg == 3) THEN
          WRITE(*,'(A, F28.16)') '   plevap               :  ', plevap(i)
        ENDIF
        WRITE(*,'(A,2F28.16)') '   zshfl/zlhfl (surface):  ', zshfl_sfc(i), zlhfl_sfc(i)
        WRITE(*,'(A, F28.16)') '   zqhfl (surface)      :  ', zqhfl_sfc(i)
        WRITE(*,'(A,2F28.16)') '   zshfl/zlhfl (soil)   :  ', zshfl_s(i), zlhfl_s(i)
        WRITE(*,'(A, F28.16)') '   zqhfl (soil)         :  ', zqhfl_s(i)
        WRITE(*,'(A,2F28.16)') '   zshfl/lhfl  (snow)   :  ', zshfl_snow(i), zlhfl_snow(i)
        WRITE(*,'(A, F28.16)') '   zqhfl (snow)         :  ', zqhfl_snow(i)
        WRITE(*,'(A, F28.16)') '   hfl_snow_soil (heat flux throw snow) :  ', hfl_snow_soil(i)
        WRITE(*,'(A, F28.16)') '   lhfl_bs              :  ', lhfl_bs(i)
        DO k = 1, ke_soil
          WRITE(*,'(A,I1,A, F28.16)') '   lhfl_pl (',k,')      :  ', lhfl_pl(i,k)
        ENDDO
        DO k = 1, ke_soil+1
          WRITE(*,'(A,I1,A, F28.16)') '   runoff_grav [kg/m2/s] (',k,')      :  ', runoff_grav(i,k)
        ENDDO

      ENDIF
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
! End of module procedure terra
!------------------------------------------------------------------------------

END SUBROUTINE terra



SUBROUTINE limit_transfer_coefficients ( &
      & ivstart, ivend, nvec, ke_soil, ke_snow, dt, u_atm, v_atm, t_atm, th_atm, rho_atm, t_g, &
      & qv_atm, qv_s, hcond_hl, t_so_now, dz_ml, dz_hl, rad_flx, w_snow_now, rho_snow_now, &
      & t_snow_now, t_snow_mult_now, tch, rho_ch, lzacc, acc_async_queue &
    )

  ! Initialisations and conversion of tch to tmch

  INTEGER, INTENT(IN) :: nvec !< Array dimensions.
  INTEGER, INTENT(IN) :: ivstart !< Start index for computations in the parallel program.
  INTEGER, INTENT(IN) :: ivend !< End index for computations in the parallel program.
  INTEGER, INTENT(IN) :: ke_soil !< Number of active soil layers.
  INTEGER, INTENT(IN) :: ke_snow !< Number of snow layers.

  REAL(wp), INTENT(IN) :: dt !< Integration time-step [s].

  REAL(wp), INTENT(IN) :: u_atm(nvec) !< Zonal wind speed in lowest level [m/s].
  REAL(wp), INTENT(IN) :: v_atm(nvec) !< Meridional wind speed in lowest level [m/s].
  REAL(wp), INTENT(IN) :: t_atm(nvec) !< Temperature in lowest level [K].
  REAL(wp), INTENT(IN) :: th_atm(nvec) !< Potential temperature in lowest level [K].
  REAL(wp), INTENT(IN) :: rho_atm(nvec) !< Air density in lowest level [kg/m^3].
  REAL(wp), INTENT(IN) :: qv_atm(nvec) !< Specific humidity in lowest level [kg/kg].
  REAL(wp), INTENT(IN) :: t_g(nvec) !< Weighted surface temperature [K].
  REAL(wp), INTENT(IN) :: qv_s(nvec) !< Specific humidity at the surface [kg/kg].
  REAL(wp), INTENT(IN) :: rad_flx(nvec) !< Net solar + thermal radiation at the ground [W/m2].
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Snow water equivalent [m H2O].
  REAL(wp), INTENT(IN) :: rho_snow_now(nvec) !< Snow density [kg/m^3]
  REAL(wp), INTENT(IN) :: t_snow_now(nvec) !< Temperature of the snow-surface [K].

  REAL(wp), INTENT(IN) :: t_snow_mult_now(nvec,0:ke_snow) !< Snow-layer temperature [K].

  REAL(wp), INTENT(INOUT) :: tch(nvec) !< Turbulent transfer coefficient for heat [1].
  REAL(wp), INTENT(INOUT) :: rho_ch(nvec) !< Density times transfer velocity for heat [kg/(m^2 s)].

  REAL(wp), INTENT(IN) :: t_so_now(nvec,0:ke_soil+1) !< Soil temperature (main level) [K].

  REAL(wp), INTENT(IN) :: hcond_hl(nvec,ke_soil) !< Thermal conductivity across half levels [W/(K m)].

  REAL(wp), INTENT(IN) :: dz_ml(ke_soil+1) !< Distance beween main levels [m].
  REAL(wp), INTENT(IN) :: dz_hl(ke_soil+1) !< Distance between half levels [m].

  LOGICAL, INTENT(IN) :: lzacc
  INTEGER, INTENT(IN) :: acc_async_queue

  ! local variables
  INTEGER :: i, m_limit

  REAL(wp) :: zuv !< Wind velocity in lowest atmospheric layer [m/s].
  REAL(wp) :: zdt_atm !< Surface-atmosphere potential temperature difference [K].
  REAL(wp) :: zdq_atm !< Surface-atmosphere specific humidity difference [kg/kg].
  REAL(wp) :: zg1 !< Heat flux between first and second soil layer [W/m^2].
  REAL(wp) :: hfl_limit !< Limit to latent or sensible heat flux [W/m^2].
  REAL(wp) :: zthfl !< Total heat flux (latent + sensible) [W/m^2].
  REAL(wp) :: zeb1 !< Unconstrained energy budget of first soil layer [W/m^2].
  REAL(wp) :: ze_melt !< Energy required to fully melt snow [J/m^2].
  REAL(wp) :: zch_snow !< Snow heat capacity [J/(m^2 K)].
  REAL(wp) :: ztchv_max !< Maximum permissible transfer velocity [m/s].
  REAL(wp) :: zdT_snow !< Maximum temperature increment in snow layer before melting sets in [K].

  REAL(wp) :: zshfl(nvec) !< Sensible heat flux estimate [W/m^2].
  REAL(wp) :: zlhfl(nvec) !< Latent heat flux estimate [W/m^2].
  REAL(wp) :: ztchv(nvec) !< Unconstrained transfer velocity [m/s].

  LOGICAL :: limit_tch(nvec) !< Indicator for flux limitation problem.

  ! Parameters
  !> approximate average soil heat capacity [J/(m^3 K)].
  REAL(wp), PARAMETER :: zch_soil = 1.4E06_wp
  !> maximum allowed temperature increment per time step in uppermost soil layer [K].
  REAL(wp), PARAMETER :: zlim_dtdt = 2.5_wp

  OPENACC_SUPPRESS_UNUSED_LZACC

  m_limit = 0

  !$ACC DATA PRESENT(ivend) CREATE(limit_tch, zshfl, zlhfl, ztchv) ASYNC(acc_async_queue)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    limit_tch(i) = .false.  !  preset problem indicator
  END DO

#ifdef __INTEL_COMPILER
!DIR$ NOFMA
#endif

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zuv, hfl_limit, zdT_snow, zg1) &
  !$ACC   PRIVATE(zthfl, ze_melt, zch_snow, zeb1, ztchv_max, zdt_atm, zdq_atm)
  DO i = ivstart, ivend
    zuv = SQRT ( u_atm(i)**2 + v_atm(i)**2 )

    zdt_atm = th_atm(i)-t_g(i)
    zdq_atm = qv_atm(i)-qv_s(i)

    ! introduce an artificical upper boundary on transfer coefficients in cases
    ! where an extreme cooling/heating of topmost soil layer may be caused due
    ! to excessive sensible/latent heat fluxes (e.g. after creation of unbalanced
    ! structures in the data assimilation or following massive changes in radiative
    ! forcing due to infrequent radiation computations)

    ! estimate current energy budget of topmost soil layer

    ! heat flux between layers 1&2 based on current temperature profile
    zg1 = hcond_hl(i,1)*(t_so_now(i,1)-t_so_now(i,2)) / dz_ml(2)

    ! estimates of sensible and latent heat flux
    zshfl(i) = tch(i)*zuv*rho_atm(i)*cp_d*zdt_atm
    zlhfl(i) = tch(i)*zuv*rho_atm(i)*lh_v*zdq_atm


    hfl_limit = MIN(500.0_wp,200.0_wp+0.5_wp*ABS(rad_flx(i)))

    IF (zshfl(i)*zlhfl(i) >= 0._wp) THEN
      zthfl = zshfl(i) + zlhfl(i)
    ELSE IF (ABS(zshfl(i)) > ABS(zlhfl(i))) THEN
      zthfl = zshfl(i) + SIGN(MIN(hfl_limit,ABS(zlhfl(i))),zlhfl(i))
    ELSE
      zthfl = zlhfl(i) + SIGN(MIN(hfl_limit,ABS(zshfl(i))),zshfl(i))
    ENDIF

    IF (ABS(zthfl) <= eps_soil) zthfl=SIGN(eps_soil,zthfl)

    ! unconstrained estimated energy budget of topmost soil layer
    zeb1 = zthfl + rad_flx(i) - zg1

    ! energy required to melt existing snow
    ze_melt = w_snow_now(i)*rho_w*lh_f     ! (J/m**2)

    ! heat capacity of snow layer, limited to a snow depth of 1.5 m
    ! for consistency with subsequent calculations
    zch_snow = MIN(w_snow_now(i),1.5_wp*rho_snow_now(i)/rho_w)*rho_w*chc_i   ! (J/(m**2 K))

    ! constrain transfer coefficient, if energy budget  of topmost soil layer is:
    ! a) negative & surface layer is unstable (i.e   upward directed turbulent heat flux)
    ! b) positive & surface layer is stable   (i.e downward directed turbulent heat flux)

    IF (zeb1<0.0_wp .AND. zthfl<0.0_wp) THEN
      ! cooling of 1st soil layer&upward SHF+LHF

      ztchv_max = ( zlim_dtdt*(zch_soil*dz_hl(1)+zch_snow)/dt    &
                   + ABS(zg1-rad_flx(i)) ) / ABS(zthfl) * tch(i)*zuv
    ELSEIF (zeb1>0.0_wp .AND. zthfl>0.0_wp) THEN
      ! heating of 1st soil layer & downward SHF+LHF
      !   Note: The heat capacity of snow is only relevant for the difference
      !         between the actual temperature and the melting point. The mean
      !         snow temperature is set to the average of t_snow & t_so(1)
      IF (lmulti_snow) THEN
        zdT_snow=MIN(0._wp, t0_melt-0.5_wp*(t_snow_mult_now(i,1)+t_so_now(i,1)))
      ELSE
        zdT_snow=MIN(0._wp,t0_melt-0.5_wp*(t_snow_now(i)+t_so_now(i,1)))
      ENDIF
      ztchv_max = ( (zlim_dtdt*zch_soil*dz_hl(1)+zdT_snow*zch_snow+ze_melt)/dt  &
                   + ABS(zg1-rad_flx(i)) ) / zthfl * tch(i)*zuv
    ELSE
      ! unlimited transfer coefficient
      ztchv_max = HUGE(1._wp)
    ENDIF
                                                    ! required constraint as non-turbulent
                                                    ! energy budget components alone may
    ztchv_max = MAX( ztchv_max, eps_soil)      ! exceed the upper limit in the energy
                                                    ! budget of the 1st soil layer

    ! Additional limitation for better numerical stability at long time steps
    ztchv_max = MIN(ztchv_max,(4._wp*zlim_dtdt*(zch_soil*dz_hl(1)+zch_snow)/dt &
                  +ABS(rad_flx(i)))/ABS(zthfl)*tch(i)*zuv)

    ztchv(i) = tch(i)*zuv  ! transfer coefficient * velocity

    LIM: IF (ztchv(i) > ztchv_max) THEN
      tch(i)=ztchv_max/MAX(zuv,1.E-06_wp)
      limit_tch(i) = .true.          ! set switch for later use
    END IF LIM

    rho_ch(i) = tch(i)*zuv*rho_atm(i) + eps_soil
  ENDDO
  !$ACC END PARALLEL

  IF (msg_level >= 20) THEN
    !$ACC UPDATE HOST(limit_tch, u_atm, v_atm, ztchv, tch, zshfl, zlhfl, t_g, t_atm) ASYNC(acc_async_queue)
    !$ACC WAIT(acc_async_queue)

    ! counter for limitation of transfer coefficients
    m_limit = COUNT( limit_tch(:) )

    ! In debugging mode and if transfer coefficient occured for at least one grid point
    IF (m_limit > 0) THEN
      WRITE(*,'(1X,A,/,2(1X,A,F10.2,A,/),1X,A,F10.2,/,1X,A,F10.3,/)')                  &
             'terra1: transfer coefficient had to be constrained',                     &
             'model time step                                 :', dt      ,' seconds', &
             'max. temperature increment allowed per time step:',zlim_dtdt,' K',       &
             'upper soil model layer thickness                :', dz_hl(1)

      DO i = ivstart, ivend
        IF (limit_tch(i)) THEN
          zuv = SQRT (u_atm(i)**2 + v_atm(i)**2 )
          WRITE(*,*) 'TERRA flux limiter: TCH before and after, zshfl, zlhf, Tsfc, Tatm , seb', &
                      ztchv(i)/zuv, tch(i), zshfl(i), zlhfl(i), t_g(i), t_atm(i)

        END IF
      END DO
    ENDIF
  ENDIF

  !$ACC END DATA

END SUBROUTINE limit_transfer_coefficients


SUBROUTINE prepare_water_budget_diagnostic ( &
      & ivstart, ivend, nvec, runoff_s, runoff_g, budget_w_so_start, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart !< Array start index.
  INTEGER, INTENT(IN) :: ivend !< Array end index.
  INTEGER, INTENT(IN) :: nvec !< Array dimension.

  REAL(wp), INTENT(IN) :: runoff_s(nvec) !< Initial surface runoff [kg/m^2].
  REAL(wp), INTENT(IN) :: runoff_g(nvec) !< Initial subsurface runoff [kg/m^2].

  REAL(wp), INTENT(OUT) :: budget_w_so_start(nvec) !< Saved runoffs [kg/m^2].

  LOGICAL, INTENT(IN) :: lzacc
  INTEGER, INTENT(IN) :: acc_async_queue

  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    budget_w_so_start(i) = runoff_s(i) + runoff_g(i)
    ! in COSMO another procedure is used than in ICON to
    !   initialize the values of runoff_s and runoff_g (COSMO:
    !   actual values, ICON: zeros) -> for ICON the following init is possible
    !   budget_w_so_start(i) = 0.0_wp
  ENDDO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE


SUBROUTINE update_water_budget_diagnostic ( &
      & ivstart, ivend, nvec, ke_soil_hy, dt, budget_w_so_start, rain_rate, snow_rate, &
      & ice_rate, &
      & qhfl_sfc, runoff_s, runoff_g, w_snow_new, w_snow_now, w_i_new, w_i_now, w_so_new, &
      & w_so_now, resid_wso, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart !< Array start index.
  INTEGER, INTENT(IN) :: ivend !< Array end index.
  INTEGER, INTENT(IN) :: nvec !< Array dimension.
  INTEGER, INTENT(IN) :: ke_soil_hy !< Number of hydrologically active soil layers/

  REAL(wp), INTENT(IN) :: dt !< Time step [s].

  REAL(wp), INTENT(INOUT) :: budget_w_so_start(nvec) !< Saved runoffs [kg/m^2].
  REAL(wp), INTENT(IN) :: rain_rate(nvec) !< Total rain rate (convective + grid-scale) [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: snow_rate(nvec) !< Total snow rate (convective + grid-scale), excluding ice [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: ice_rate(nvec) !< Ice precipitation rate [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: qhfl_sfc(nvec) !< Vapor flux at surface [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: runoff_s(nvec) !< Final surface runoff [kg/m^2].
  REAL(wp), INTENT(IN) :: runoff_g(nvec) !< Final subsurface runoff [kg/m^2].
  REAL(wp), INTENT(IN) :: w_snow_new(nvec) !< Final snow-water equivalent [m H2O].
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Initial snow-water equivalent [m H2O].
  REAL(wp), INTENT(IN) :: w_i_new(nvec) !< Final interception water amount [m H2O].
  REAL(wp), INTENT(IN) :: w_i_now(nvec) !< Initial interception water amount [m H2O].
  REAL(wp), INTENT(IN) :: w_so_new(nvec, ke_soil_hy) !< Final total water in soil layer [m H2O].
  REAL(wp), INTENT(IN) :: w_so_now(nvec, ke_soil_hy) !< Initial total water in soil layer [m H2O].

  REAL(wp), INTENT(INOUT) :: resid_wso(nvec) !< Total water residual [kg/m^2].

  LOGICAL, INTENT(IN) :: lzacc
  INTEGER, INTENT(IN) :: acc_async_queue

  INTEGER :: i
  INTEGER :: k

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    ! the rain rates
    budget_w_so_start(i) = budget_w_so_start(i) + (rain_rate(i) + snow_rate(i) + ice_rate(i))*dt
    ! the evapotranspiration
    budget_w_so_start(i) = budget_w_so_start(i) + qhfl_sfc(i)*dt
    ! surface + subsurface runoff (subtraction because water is lost)
    ! The initial values were saved on entry to terra.
    budget_w_so_start(i) = budget_w_so_start(i) - runoff_s(i) - runoff_g(i)
    ! snow + interception storage (convert m H2O to kg/m^2)
    budget_w_so_start(i) = budget_w_so_start(i) - ( w_snow_new(i) - w_snow_now(i) ) * rho_w
    budget_w_so_start(i) = budget_w_so_start(i) - ( w_i_new(i)    - w_i_now(i)    ) * rho_w
    ! residuum = recharge - budget (RHS)
    resid_wso(i) = resid_wso(i) - budget_w_so_start(i)
  END DO

  !$ACC LOOP SEQ
  DO k = 1, ke_soil_hy
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ! Sum up water soil water content difference.
      resid_wso(i) = resid_wso(i) + (w_so_new(i,k) - w_so_now(i,k)) * rho_w
    END DO
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE

END MODULE sfc_terra_main
