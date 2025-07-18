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

! Data module for all parametric data in the soil model "terra"!

MODULE sfc_terra_data

!------------------------------------------------------------------------------
!
! Modules used:

USE mo_kind, ONLY: wp

!==============================================================================

IMPLICIT NONE

PUBLIC           ! All constants and variables in this module are public

!==============================================================================

! Global (i.e. public) Declarations:

  ! BATS parameters
  REAL  (KIND=wp) , PARAMETER ::  &
    crhowm     =    0.8_wp    , & !  BATS (1)
    cdmin      =    0.25E-9_wp, & !  BATS (m**2/s)
    cfinull    =    0.2_wp    , & !  BATS (m)
    ckrdi      =    1.0E-5_wp , & !  BATS (m/s)
    cdash      =    0.05_wp   , & !  BATS ((m/s)**1/2)
    clai       =    3.0_wp    , & !  BATS
    cparcrit   =  100.0_wp    , & !  BATS (W/m**2)
    ctend      =  313.15_wp   , & !  BATS (K)
    csatdef    = 4000.0_wp        !  BATS (Pa)

! 1. Data arrays for properties of different soil types (array index)
! -------------------------------------------------------------------
  ! Soil type indexes
  INTEGER, PARAMETER :: IST_ICE = 1 !< Soil type: ice.
  INTEGER, PARAMETER :: IST_ROCK = 2 !< Soil type: rock.
  INTEGER, PARAMETER :: IST_SAND = 3 !< Soil type: sand.
  INTEGER, PARAMETER :: IST_SLOAM = 4 !< Soil type: sandy loam.
  INTEGER, PARAMETER :: IST_LOAM = 5 !< Soil type: loam.
  INTEGER, PARAMETER :: IST_CLOAM = 6 !< Soil type: clay loam.
  INTEGER, PARAMETER :: IST_CLAY = 7 !< Soil type: clay.
  INTEGER, PARAMETER :: IST_PEAT = 8 !< Soil type: peat.
  INTEGER, PARAMETER :: IST_SEAWTR = 9 !< Soil type: sea water.
  INTEGER, PARAMETER :: IST_SEAICE = 10 !< Soil type: sea ice.
  INTEGER, PARAMETER :: IST_NUM = 10 !< Number of soil types.

  ! Scheme IDs
  !> BATS bare-soil evaporation.
  INTEGER, PARAMETER :: EVSL_BATS = 2
  !> Noilhan and Planton (1989) bare-soil evaporation.
  INTEGER, PARAMETER :: EVSL_NP89 = 3
  !> Resistance-based formulation [Schulz & Vogel (2020)].
  INTEGER, PARAMETER :: EVSL_RESIST = 4
  !> Resistance-based formulation with `c_soil=2` and `cr_bsmin` as tuning parameter.
  INTEGER, PARAMETER :: EVSL_RESIST_RBS = 5

  !> Heat conductivity based on assumption of a soil water content which is equal to the average
  !! between wilting point and field capacity.
  INTEGER, PARAMETER :: HCOND_AVG = 1
  !> Heat conductivity based on Peters-Lidard et al. (1998)
  INTEGER, PARAMETER :: HCOND_PL98 = 2
  !> Heat conductivity based on Peters-Lidard et al. (1998) with modified conductivity in upper soil
  !! due to vegetation.
  INTEGER, PARAMETER :: HCOND_PL98_VEG = 3

  INTEGER, PARAMETER :: TRVG_BATS = 2 !< BATS transpiration.
  INTEGER, PARAMETER :: TRVG_BATS_EXT = 3 !< Extended BATS transpiration scheme.

  INTEGER, PARAMETER :: ROOT_CONSTANT = 1 !< Constant root density.
  INTEGER, PARAMETER :: ROOT_EXPONENTIAL = 2 !< Exponentially decaying root density.

  INTEGER, PARAMETER :: itype_mire = 0 !< Switch for mire parameterization

  ! Initialization of soil type parameters
  INTEGER, PRIVATE :: i

  ! soil type:   ice    rock    sand    sandy   loam   clay      clay    peat    sea     sea
  ! (by index)                          loam           loam                     water    ice

  REAL(wp), PARAMETER, DIMENSION(IST_NUM) :: &
  ! a) parameters describing the soil water budget
    !> pore volume (fraction of volume)
    cporv = [ 1.E-10_wp, 1.E-10_wp, 0.364_wp  , 0.445_wp  , 0.455_wp  , 0.475_wp  , 0.507_wp  , 0.863_wp  , 1.E-10_wp, 1.E-10_wp], &
    !> field capacity (fraction of volume)
    cfcap = [ 1.E-10_wp, 1.E-10_wp, 0.196_wp  , 0.260_wp  , 0.340_wp  , 0.370_wp  , 0.463_wp  , 0.763_wp  , 1.E-10_wp, 1.E-10_wp], &
    !> plant wilting point (fraction of volume)
    cpwp  = [ 0.0_wp   , 0.0_wp   , 0.042_wp  , 0.100_wp  , 0.110_wp  , 0.185_wp  , 0.257_wp  , 0.265_wp  , 0.0_wp   ,  0.0_wp  ], &
    !> air dryness point (fraction of volume)
    cadp  = [ 0.0_wp   , 0.0_wp   , 0.012_wp  , 0.030_wp  , 0.035_wp  , 0.060_wp  , 0.065_wp  , 0.098_wp  , 0.0_wp   ,  0.0_wp  ], &
    !> minimum infiltration rate (kg/s*m**2)
    cik2  = [ 0.0_wp   , 0.0_wp   , 0.0035_wp , 0.0023_wp , 0.0010_wp , 0.0006_wp , 0.0001_wp , 0.0002_wp , 0.0_wp   ,  0.0_wp  ], &
    !> parameter for determination of hydr. conductivity (m/s)
    ckw0  = [ 0.0_wp   , 0.0_wp   , 479.E-7_wp, 943.E-8_wp, 531.E-8_wp, 764.E-9_wp, 85.E-9_wp , 58.E-9_wp , 0.0_wp   ,  0.0_wp  ], &
    !> parameter for determination of hydr. conductivity (1)
    ckw1  = [ 0.0_wp   , 0.0_wp   , -19.27_wp , -20.86_wp , -19.66_wp , -18.52_wp , -16.32_wp , -16.48_wp , 0.0_wp   ,  0.0_wp  ], &
    !> parameter for determination of hydr. diffusivity (m**2/s)
    cdw0  = [ 0.0_wp   , 0.0_wp   , 184.E-7_wp, 346.E-8_wp, 357.E-8_wp, 118.E-8_wp, 442.E-9_wp, 106.E-9_wp, 0.0_wp   ,  0.0_wp  ], &
    !> parameter for determination of hydr. diffusivity (1)
    cdw1  = [ 0.0_wp   , 0.0_wp   , -8.45_wp  , -9.47_wp  , -7.44_wp  , -7.76_wp  , -6.74_wp  , -5.97_wp  , 0.0_wp   ,  0.0_wp  ], &
    !> rock/ice/water indicator (hydrological calculations only for crock=1)
    crock = [ 0.0_wp   , 0.0_wp   , 1.0_wp    , 1.0_wp    , 1.0_wp    , 1.0_wp    , 1.0_wp    , 1.0_wp    , 0.0_wp   ,  0.0_wp  ], &
  ! b) parameters describing the soil heat budget
    !> soil heat capacity  (J/K*m**3)
    crhoc = [ 1.92E6_wp, 2.10E6_wp, 1.28E6_wp , 1.35E6_wp , 1.42E6_wp , 1.50E6_wp , 1.63E6_wp , 0.58E6_wp , 4.18E6_wp, 1.92E6_wp], &
    !> parameter for the determination of the soil heat conductivity (W/(K*m))
    cala0 = [ 2.26_wp  , 2.41_wp  , 0.30_wp   , 0.28_wp   , 0.25_wp   , 0.21_wp   , 0.18_wp   , 0.06_wp   , 1.0_wp   ,  2.26_wp ], &
    !> parameter for the determination of the soil heat conductivity (W/(K*m))
    cala1 = [ 2.26_wp  , 2.41_wp  , 2.40_wp   , 2.40_wp   , 1.58_wp   , 1.55_wp   , 1.50_wp   , 0.50_wp   , 1.0_wp   ,  2.26_wp ], &
    !> slope of solar albedo with respect to soil water content
    csalbw = [0.00_wp  , 0.00_wp  , 0.44_wp   , 0.27_wp   , 0.24_wp   , 0.23_wp   , 0.22_wp   , 0.10_wp   , 0.00_wp  ,  0.00_wp ], &
  ! c) additional parameters for soil water content dependent freezing/melting
    !> mean fraction of sand (weight percent)
    csandf = [0.0_wp   , 0.0_wp   , 90._wp    , 65._wp    , 40._wp    , 35._wp    , 15._wp    , 90._wp    , 0.00_wp  ,  0.00_wp ], &
    !> mean fraction of clay (weight percent)
    cclayf = [0.0_wp   , 0.0_wp   , 5.0_wp    , 10._wp    , 20._wp    , 35._wp    , 70._wp    , 5.0_wp    , 0.00_wp  ,  0.00_wp ], &
  ! d) additional parameters for the BATS scheme (Dickinson)
    !>  (m/s)
    ck0di = [ 1.E-4_wp , 1.E-4_wp , 2.E-4_wp  , 2.E-5_wp  , 6.E-6_wp  , 2.E-6_wp  , 1.E-6_wp  , 1.5E-6_wp , 0.00_wp  ,  0.00_wp ], &
    !>  (1)
    cbedi = [ 1.00_wp  , 1.00_wp  , 3.5_wp    , 4.8_wp    , 6.1_wp    , 8.6_wp    , 10.0_wp   , 9.0_wp    , 0.00_wp  ,  0.00_wp ], &
    !>  auxiliary variable
    clgk0 = [ (LOG10(MAX(1.0E-6_wp,ck0di(i)/ckrdi)), i = 1, IST_NUM) ]
  ! options for diffuse solar albedo; selection is made in mo_radiation_nml.
  ! These cannot be PARAMETERs because the pointer below points to one of them.
  REAL(wp), TARGET, DIMENSION(IST_NUM) :: &
    csalb1 = [0.70_wp  , 0.30_wp  , 0.30_wp   , 0.25_wp   , 0.25_wp   , 0.25_wp   , 0.25_wp   , 0.20_wp   , 0.07_wp  ,  0.70_wp ], &
    csalb2 = [0.70_wp  , 0.30_wp  , 0.30_wp   , 0.25_wp   , 0.25_wp   , 0.25_wp   , 0.25_wp   , 0.20_wp   , 0.06_wp  ,  0.70_wp ]
  !$ACC DECLARE COPYIN(csalb1, csalb2)

  !> Diffuse solar albedo (IST_NUM) [1].
  REAL(wp), POINTER :: csalb(:)
  !$ACC DECLARE CREATE(csalb)

! 2. Additional parameters for the soil model
! -------------------------------------------------------------------

  REAL  (KIND=wp), PARAMETER :: &
!==============================================================================

    csalb_p        = 0.15_wp  , & !  solar albedo of ground covered by plants
    csalb_snow     = 0.70_wp  , & !  solar albedo of ground covered by snow

    ! T.R. 2011-09-21 csalb_snow_min/max set to values used in GME
    csalb_snow_min = 0.500_wp , &
                           ! min. solar albedo of snow for forest free surfaces
    csalb_snow_max = 0.800_wp , &
                           ! max. solar albedo of snow for forest free surfaces
    ! T.R. 2011-09-21 snow albedos for forests set to values used in GME
    csalb_snow_fe  = 0.270_wp , &  ! solar albedo of snow for surfaces with evergreen forest
    csalb_snow_fd  = 0.320_wp , &  ! solar albedo of snow for surfaces with deciduous forest

    ctalb          = 0.004_wp , & !  thermal albedo ( of all soil types )
    cf_snow        = 0.0150_wp, & !  parameter for the calculation of the
                                  !  fractional snow coverage
  ! for the multi-layer soil model
    cwhc       = 0.04_wp      , & !  water holding capacity of snow ()
    chcond     = 0.01_wp      , & !  saturation hydraulic conductivity of snow ()
    ca2        = 6.6E-07_wp   , & !  activation energy (for snow metamorphosis) (J)
    csigma     = 75._wp       , & !  snow metamorphosis, Pa

  ! cf_w changed from 0.0004 to 0.0010 (in agreement with GME)
    cf_w       = 0.0010_wp    , & !  parameter for the calculation of the
                                  !  fractional water coverage

    csvoro     = 1.0000_wp    , & !  parameter to estimate the subgrid-scale
                                  !  variation of orography
    cik1       = 0.0020_wp    , & !  parameter for the determination of the
                                  !  maximum infiltaration
    cakw       = 0.8000_wp    , & !  parameter for averaging the water contents
                                  !  of the top and middle soil water layers to
                                  !  calculate the hydraulic diffusivity and
                                  !  conductiviy

    ctau1      = 1.0000_wp    , & !  first adjustment time period in EFR-method
    ctau2      = 5.0000_wp    , & !  second adjustment time period in EFR-method
    chc_i      = 2100.0_wp    , & !  heat capacity of ice
    chc_w      = 4180.0_wp    , & !  heat capacity of water

    cdzw12     = 0.1000_wp    , & !  thickness of upper soil water layer in
                                  !  two-layer model
    cdzw22     = 0.9000_wp    , & !  thickness of lower soil water layer in
                                  !  two-layer model
    cdzw13     = 0.0200_wp    , & !  thickness of upper soil water layer in
                                  !  three-layer model
    cdzw23     = 0.0800_wp    , & !  thickness of middle soil water layer in
                                  !  three-layer model
    cdzw33     = 0.9000_wp        !  thickness of lower soil water layer in
                                  !  three-layer model

  ! Monolithic TERRA modifies this
  REAL (KIND=wp) :: ctau_i = 7200.0_wp !< time constant for the drainage from the interception storage

  REAL (KIND=wp), PARAMETER ::  &
    cdsmin     = 0.0100_wp    , & !  minimum snow depth
    crhosmin   = 500.00_wp    , & !  minimum density of snow
    crhosmax   = 800.00_wp    , & !  maximum density of snow
    crhosmin_ml=  50.00_wp    , & !  minimum density of snow
    crhosmax_ml= 400.00_wp    , & !  maximum density of snow
    crhosminf  =  50.00_wp    , & !  minimum density of fresh snow
    crhosmaxf  = 150.00_wp    , & !  maximum density of fresh snow
    crhogminf  = 100.00_wp    , & !  minimum density of fresh graupel / convective snow
    crhogmaxf  = 200.00_wp    , & !  maximum density of fresh graupel / convective snow
    crhosmint  =  0.125_wp    , & !  value of time constant for ageing
                                  !  of snow at csnow_tmin (8 days)
    crhosmaxt  =   0.40_wp    , & !  maximum value of time constant for ageing
                                  !  of snow
    crhosmax_tmin = 200.00_wp , & ! maximum density of snow at csnow_tmin
    csnow_tmin = 258.15_wp    , & !  lower threshold temperature of snow for
                                  !  ageing and fresh snow density computation
                                  !  ( = 273.15-15.0)
    crhos_dw   = 300.00_wp    , & !  change of snow density with water content
    calasmin   = 0.2000_wp    , & !  minimum heat conductivity of snow (W/m K)
    calasmax   = 1.5000_wp    , & !  maximum heat conductivity of snow (W/m K)
    calas_dw   = 1.3000_wp    , & !  change of snow heat conductivity with
                                  !  water content                (W/(m**2) K)

    !Minimum and maximum value of stomatal resistance (s/m)
    !used by the Pen.-Mont. method for vegetation transpiration
    !(itype_trvg=2):
    crsmin     = 150.0_wp     , & !  BATS (s/m)
    crsmax     = 4000.0_wp        !  BATS (s/m)
    ! crsmax increased from 1000 to 4000 s/m (to reduce latent heat flux).


! 3. Additional variables for the soil geometry
! ---------------------------------------------

  ! these are allocated and computed in sfc_init
  REAL  (KIND=wp), ALLOCATABLE ::    &
    zzhls          (:)             , & ! depth of the half level soil layers in m
    zdzhs          (:)             , & ! layer thickness between half levels
    zdzms          (:)                 ! distance between main levels


! 4. Variables for TERRA_URB
! --------------------------

  ! Default urban fabric parameters are derived according to literature.
  ! They are obtained/tested for
  !      Toulouse, Basel (Wouters et al., 2015),
  !      Paris (De Ridder et al., 2013;
  !      Sarkar and De Ridder, 2010;
  !      Demuzere et al., 2008),

  REAL  (KIND=wp), PARAMETER ::       &
    ctalb_bm = 0.08_wp   , & !  default effective thermal albedo of building/road environment

    ! csalb_eff_uf = 0.80_wp,   & ! correction factor for effective albedo induced by the
    !                             ! urban fabric (street canyons).
    ! ! This value is based on observations and monte-carlo simulations, taking H/W ratio of 1.0
    ! ! and roof fraction of 0.5, see Pawlak et al.
    ! ! ????: http://nargeo.geo.uni.lodz.pl/~icuc5/text/P_4_6.pdf. csalb_eff_uf = (20% + 15%)/2

    csalb_bm = 0.213_wp  , & ! default short-wave albedo of building/road materials
                             ! in the urban fabric
      ! csalb_bm is chosen in such a way that the effective albedo for a dense urban
      ! environment csalb_bm * csalb_eff_uf is equal to 0.17

  ! Default surface area index of buildings/street environment. This is esimated from the
  ! squared thermal inertia = 3.8E6 (for Paris, see De Ridder et al.,2013) estimated from
  ! model simulations divided by the building material parameters below estimates
  ! c_rhoc_bm * c_ala_bm

    ! c_ai_uf = 2.0_wp,    &
    c_uf_h  = 15._wp     , & ! default height of building elements in urban fabric (metre)
    c_lnd_h = 0.01_wp    , & ! default height of natural soil elements in natural land (metre)
    c_htw   = 1.5_wp     , &
    c_roof  = 0.667_wp   , &

  ! default building/road 'material' properties...

    ! Value of 'concrete' is taken (see engineering-toolbox.com)
    c_rhoc_bm = 1.74E6_wp, & ! default specific heat times density of buildings,

    !  Value of 'medium concrete' is taken (higher boundary, as average value for
    !  concrete in,  see engineerintoolbox.com)
    c_ala_bm  = 0.87_wp  , & ! default building material heat conductivity:
    c_isa_runoff = 1.0_wp, & ! default fraction of water exiting from the impervious
                               ! surface leading to runoff. The remainder fraction is (potentially)
                               ! infiltrated in the neighbouring natural soil (switched off by default).
                               ! Some addaptation strategies (such as infiltration of roof water)
                               ! will lead to values less than 1
    c_isa_delt = 0.12_wp , & ! The maximum wet-surface fraction delta_max, see Wouters et al., 2015
    cwisamax = 1.31E-3_wp, & ! Maximum amount of water that can be stored by impervious surfaces  (mH2O)
                             ! (estimate for urban areas, Toulouse centre), see Wouters et al., 2015

    ! constants for calculation of anthropogenic heat, see Flanner, 2009
    cb1        = 0.451_wp, & !
    cb2        = 0.8_wp  , & !
    csig       = 0.3_wp  , & ! #sigma =0.18
    cmu        = 0.5_wp  , & !
    cA1        = -0.5_wp , & ! #A1 = -0.3
    cff        = 2.0_wp  , & ! #ff = 1.9
    calph      = 10.0_wp , & !
    ceps       = 0.25_wp     !



! 5. Additional control variables
! -------------------------------

  LOGICAL                   ::  &

    lsoilinit_dfi = .FALSE.         ! initialize soil after dfi forward launching

! 6. Epsilons (security constants)
! --------------------------------

  REAL  (KIND=wp), PARAMETER :: &

    ! Avoid division by zero, e.g. x = y / MAX(z,eps_div).
    eps_div  = 1.0E-6_wp      , &
!!    eps_div  = repsilon       , &
!! RUS
!! eps_div is used in divisions to avoid division by zero, repsilon would
!! be appropriate therefore. However, the testsuite fails with repsilon (about 1E-30)
!! as it is 'used to' the (in this context) huge epsilon of 1E-6.

    ! Multi-purpose epsilon in soil model (former zepsi).
    eps_soil = 1.0E-6_wp      , &

    ! Small value to check if temperatures have exceeded a fixed threshold
    ! such as the freezing point.  In double precision (15 decimal digits)
    ! a value as small value such as 1.0E-6 can be used. In single
    ! precision (6-7 decimal digits), however, the value has to be larger
    ! in order not to vanish. The current formulation is save for
    ! temperatures up to 500K.
    eps_temp = MAX(1.0E-6_wp,500.0_wp*EPSILON(1.0_wp)), &

    ! Extremely small value in order to prevent a floating point underflow
    ! in double precision.
    eps_nounderflow = 1.0E-5_wp * EPSILON(1.0_wp)


! 8. Soil ice parameterization
! ----------------------------

  REAL  (KIND=wp), PARAMETER ::  &
    T_ref_ice  =  0.1_wp,        & !degC Soil ice parameterization
    T_star_ice =  0.01_wp,       & !degC according to K. Schaefer and Jafarov, E.,2016, doi:10.5194/bg-13-1991-2016
    b_clay     = -0.3_wp,        &
    b_silt     = -0.5_wp,        &
    b_sand     = -0.9_wp,        &
    b_org      = -1.0_wp

  REAL(wp), PARAMETER :: rho_i = 910.0_wp !< Density of solid ice (soil model) [kg/m**3].

END MODULE sfc_terra_data
