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

! @brief: Configuration setup for turbulent processes (turbdiff)
!
! Configuration setup for the entire turblence parameterization (by M. Raschendorfer)
!  based on the modules 'turb_[diffusion|transfer|vertdiff|utilities]:
!
! This module 'mo_turbdiff_config' paricularly includes the whole content of the now obsolete module 'turb_data',
!  that is the declaration and initialization of all required parameters, switches and selectors.
!
! The module also contains a mapping of at least all included variable quantities that define the configuration-state
!  of the scheme to associated components of the domain-specific config-state vector 'turbdiff_config(jg)'.
! Some of these quantities are included into the namelist setting in module 'mo_turbdiff_nml'.

MODULE mo_turbdiff_config

! ---------------------------------------------------------------
! Declarations:
!
! Modules used:
! ---------------------------------------------------------------

USE mo_kind,               ONLY: wp           ! KIND-type parameter for real variables
USE mo_impl_constants,     ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

!===================================================================================

IMPLICIT NONE

PUBLIC

!!--------------------------------------------------------------------------------------------------
!! Basic configuration setup for atmospheric turbulence (turbdiff) and turbulent transfer (turbtran)
!!  including vertical diffusion (vertdff):
!!--------------------------------------------------------------------------------------------------

!#ifdef new_vers
!===================================================================================
! Fixed configuration parameters:
! -------------------------------
INTEGER, PARAMETER :: &
   ntmax=3,          & !maxmal number of time-levels for TKE
!
!  Indices for the two discriminated variable-types:
!
   mom=1,       & !momentum variables
   sca=mom+1,   & !scalar   variables
   ntyp=sca             !related number of variable-types ('mom' and 'sca')
!
! Indices associated to paricular model variables:
! ------------------------------------------------
!
INTEGER, PARAMETER :: &
   u_m=1,       & !zonal velocity-component at the mass center
   v_m=u_m+1,   & !meridional ,,      ,,    ,, ,,   ,,    ,,
   nvel=v_m,          & !number of velocity-components active for turbulece ('u_m', 'v_m')

   tem=nvel+1,  & !                       temperature
   tem_l=tem,   & !liquid-water           temperature
   tet=tem,     & !             potential temperature
   tet_l=tet,   & !liquid-water potential temperature

   vap=tem+1,   & !water vapor (specific humidity 'qv')
   h2o_g=vap,   & !total water ('qv+qc')
   nred=h2o_g,        & !number of progn. variables being active for turbulence after reduction by local sat.adj.
   ninv=nred-nvel,    & !number of scalar variables being conserved during 'vap'<->'liq'-transistions ('tet_l', 'h2o_g')

   liq=vap+1,   & !liquid water (mass fraction 'qc')
   nmvar=liq,         & !number of progn. variables beding active for turbulence, which is equal to the
                        !number of variables being included into the single-column turbulence statistics
   nscal=liq-nvel,    & !number of scalar variables being active for turbulence ('tem', 'vap', 'liq')

   w_m=nmvar,   & !vertical velocity-component (only used for optional signle-column turbulence statistics)

!  Notice that: 'u_m, v_m'        are within [1, nvel];
!  but:         'tet_l, h2o_g'    are within [nvel+1, nred]
!  and:         'tem (tet), vap, liq' within [nvel+1, nmvar],
!  while:       'w_m'            is equal to 'nred+1=nmvar=liq'.

   naux=5,            & !(positive) number of auxilary variables

   ndim=MAX(nmvar,naux) !(positive) limit of last dimension used for 'zaux' and 'zvari'

!  Note:
!  The index "0" for the last dimension of 'zaux' and 'zvari' is also used, and it refers to
!   'zaux' : saturation fraction (cl_cv)
!   'zvari': potential available energy per volume (pressure) on half levels or
!            vertical gradient of effectively availabel kinetic energy (per mass)
!             due to near surface thermal inhomogeneity, being called here Circulation Kinetic Energy (CKE)

LOGICAL, PARAMETER :: &
   ldynimp=.FALSE., &   !dynamical calculation of implicit weights for semi-implicit vertical diffusion
   lprecnd=.FALSE., &   !preconditioning of tridiagonal matrix      ,,     ,,          ,,        ,,

   lporous=.FALSE., &   !Vertically Resolved Roughness Layer (VRRL) representing a porous atmospheric medium
   !Note: The VRRL-treatment is not yet complete!

   ltst2ml =.FALSE., &   !test required, whether  2m-level is above the lowest main-level
   ltst10ml=.FALSE.      !test required, whether 10m-level is above the lowest half-level

   !Attention:
   !So far, the  2m-level is assumed to be always below the lowest main-level,
   !    and the 10m-level is assumed to be always below the lowest half-level!!
!#endif

!===================================================================================
! Parameters that may be used for tuning, special configurations and namelist-input:
! ----------------------------------------------------------------------------------

! Attention:
! The given initializations are default settings of the boundary layer parameters.
! Some of these initial parameter values may be changed afterwards by model input NAMELISTs.

! 1. Numerical parameters:
!-------------------------

  REAL(wp):: impl_s        =  1.20_wp   ! implicit weight near the surface (maximal value)
  REAL(wp):: impl_t        =  0.75_wp   ! implicit weight near top of the atmosphere (maximal value)

  INTEGER :: imode_tkvmini = 2          ! mode of calculating the minimal turbulent diff. coeffecients
                           ! 1: with a constant value 'tk[h|m]min*'
                           ! 2: with a stability dependent correction
  ! Minimal diffusion coefficients in [m^2/s] for vertical
  REAL(wp):: tkhmin        =  0.75_wp   ! scalar (heat) transport
  REAL(wp):: tkmmin        =  0.75_wp   ! momentum transport
  REAL(wp):: tkhmin_strat  =  0.75_wp   ! scalar (heat) transport, enhanced value for stratosphere
  REAL(wp):: tkmmin_strat  =  4.00_wp   ! momentum transport,      enhanced value for stratosphere

  REAL(wp):: ditsmot       =  0.00_wp   ! smoothing factor for direct time-step iterations

  INTEGER :: imode_frcsmot = 2          ! if "frcsmot>0", apply smoothing of TKE source terms
                           ! 1: globally or
                           ! 2: in the tropics only (if 'trop_mask' is present)
  REAL(wp):: frcsmot       =  0.00_wp   ! vertical smoothing factor for TKE forcing
  REAL(wp):: tkesmot       =  0.15_wp   ! time smoothing factor for TKE and diffusion coefficients

  REAL(wp):: frcsecu       =  1.00_wp   ! security factor for TKE-forcing       (<=1)
  REAL(wp):: tkesecu       =  1.00_wp   ! security factor in  TKE equation      (out of [0; 1])
  REAL(wp):: stbsecu       =  0.01_wp   ! security factor in stability function (out of ]0; 1])
  REAL(wp):: prfsecu       =  0.50_wp   ! relat. secur. fact. for prof. funct.  (out of ]0; 1[)

  REAL(wp):: epsi          = 1.0E-6_wp  ! relative limit of accuracy for comparison of numbers

  INTEGER :: it_end        = 1          ! number of initialization iterations (>=0)

! 2. Parameters describing physical properties of the lower boundary of the atmosphere:
!--------------------------------------------------------------------------------------

  REAL(wp):: rlam_heat     = 10.0_wp    ! scaling factor of the laminar boundary layer for heat
  REAL(wp):: rlam_mom      =  0.0_wp    ! scaling factor of the laminar boundary layer for momentum
                                        ! (should remain at "0" with the current formulation!)

  REAL(wp):: rat_lam       =  1.0_wp    ! vapour/heat ratio of laminar scaling factors (over land)
  REAL(wp):: rat_sea       =  0.8_wp    ! sea/land ratio of laminar scaling factors for heat (and vapor)
  REAL(wp):: rat_glac      =  3.0_wp    ! glacier/land ratio of laminar scaling factors for heat (and vapor)

  REAL(wp):: rat_can       =  1.0_wp    ! ratio of canopy height over sai*z0m

  INTEGER :: imode_nsf_wind= 1          ! mode of local wind-definition at near-surface levels (related to 'rsur_sher')
                                        !  (applied to 10m wind-diagnostics as well as to calculation of sea-surface roughness)
                           ! 1: ordinary wind speed (magnitude of grid-scale averaged wind-vector)
                           ! 2: local wind speed related to additional surface-shear by NTCs|LLDCs (at "rsur_sher>0")

  ! Scaling factor (out of [0, 1]) representing the considered fraction of additional shear-forcing by
  !  Non-Turbulent subgrid Circulations (NTCs) or via Lower Limits of Diffusion-Coefficients (LLDCs)
  !  being transmitted from level "P" (k=ke) to level "0" (k=ke1): top of R-layer by Land-Use (LU):
  REAL(wp):: rsur_sher     =  0.0_wp    ! (so far deactivated)
  ! Notice:
  ! Through "rsur_sher>0", a related amplification of shear-forcing at the surface is active.
  ! The LLDC-part is only included, if it is considered to substitue missing shear-forcing,
  !  which is expressed by "imode_tkemini=2"!

  INTEGER :: imode_charpar = 2          ! mode of estimating the Charnock-Parameter
                           ! 1: use a constant value
                           ! 2: use a wind-dependent value with a constant upper bound
                           ! 3: as "2", but with reduction at wind speeds above 25 m/s for more realistic TC wind speeds
  REAL(wp):: alpha0        =  0.0123_wp ! Charnock-parameter
  REAL(wp):: alpha0_max    =  0.0335_wp ! upper limit of velocity-dependent Charnock-parameter
  REAL(wp):: alpha0_pert   =  0.0_wp    ! additive ensemble perturbation of Charnock-parameter

  REAL(wp):: alpha1        =  0.7500_wp ! parameter scaling the molecular roughness of water waves

! 3. Parameters that should be external parameter fields being not yet available:
! -------------------------------------------------------------------------------

  REAL(wp):: c_lnd         = 2.0_wp      ! surface area density of the roughness elements over land
  REAL(wp):: c_sea         = 1.5_wp      ! surface area density of the waves over sea
  REAL(wp):: c_soil        = 1.0_wp      ! surface area density of the (evaporative) soil surface
  REAL(wp):: c_stm         = 0.0_wp      ! (so far deactivated)
! REAL(wp):: c_stm         = 2.5_wp      ! surface area density of stems and branches at the plant-covered part of the surface
  !Note:
  !"c_stm=0" matches with the previous not consistant scaling, considering "c_stm=2.5" only as an additional surface
  ! part that can hold interception water.
  REAL(wp):: e_surf        = 1.0_wp      ! exponent to get the effective surface area

  LOGICAL :: lconst_z0     = .FALSE.     ! horizontally homogeneous roughness length (for idealized testcases!) applied

  REAL(wp):: const_z0      =  0.001_wp   ! horizontally homogeneous roughness length (for idealized testcases!)

! 4. Parameters that should be dynamical fields being not yet available:
! ----------------------------------------------------------------------

  REAL(wp):: z0m_dia       =  0.2_wp     ! roughness length of a typical synoptic station [m]
  REAL(wp):: z0_ice        =  0.001_wp   ! roughness length of sea ice

! 5. Parameters for modelling turbulent diffusion:
! ------------------------------------------------

  REAL(wp):: tur_len       =  500.0_wp   ! asymptotic maximal turbulent distance [m]
  REAL(wp):: pat_len       =  100.0_wp   ! effective global length scale of subscale surface patterns over land [m]
                                         ! (should be dependent on location)
  REAL(wp):: len_min       =  1.0E-6_wp  ! minimal turbulent length scale [m]

  INTEGER :: imode_vel_min = 2           ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                           ! 1: with a constant value "tkesecu*vel_min"
                           ! 2: with a stability dependent correction of "tkesecu*vel_min"
  REAL(wp):: vel_min       =  0.01_wp    ! minimal velocity scale [m/s]
  REAL(wp):: vel_max       =  30.0_wp    ! maximal velocity scale [m/s]

  REAL(wp):: akt           =  0.4_wp     ! von-Karman constant

  ! Length-scale factors for pressure destruction of turbulent
  REAL(wp):: a_heat        =  0.74_wp    ! scalar (heat) transport
  REAL(wp):: a_mom         =  0.92_wp    ! momentum transport

  ! Length-scale factors for dissipation of turbulent
  REAL(wp):: d_heat        =  10.1_wp    ! scalar (temperature) variance
  REAL(wp):: d_mom         =  16.6_wp    ! momentum variance

  ! Length-scale factor for turbulent transport (vertical diffusion) of TKE
  REAL(wp):: c_diff        =  0.20_wp    ! (including turb. pressure-transport)

  ! Length-scale factor for the stability correction
  REAL(wp):: a_stab        =  0.00_wp    ! applied to integral turbulent length-scale

  ! Length-scale factor for separate horizontal shear circulations (related 'ltkeshs')
  REAL(wp):: a_hshr        =  1.00_wp    ! contributing to shear-production of TKE

  ! Dimensionless parameters used in the sub grid scale condensation scheme
  ! (statistical cloud scheme):

  REAL(wp):: clc_diag      =  0.5_wp     ! cloud cover at saturation
  REAL(wp):: q_crit        =  1.6_wp     ! critical value for normalized super-saturation

  REAL(wp):: c_scld        =  1.0_wp     ! shape-factor (0<=c_scld) applied to pure 'cl_cv' at the moist correct.
                                         !  by turbulent phase-transit., providing an eff. 'cl_cv', which scales
                                         !  the implicit liquid-water flux under turbulent sat.-adj.:
                        !  <1: small eff. 'cl_cv' even at large pure 'cl_cv'
                        !  =1:       eff. 'cl_cv' just equals   pure 'cl_cv'
                        !  >1: large eff. 'cl_cv' even at small pure 'cl_cv'

!===================================================================================
! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

  LOGICAL :: ltkesso       = .TRUE.      ! consider mechanical SSO-wake production in TKE-equation
  LOGICAL :: ltkecon       = .FALSE.     ! consider convective buoyancy production in TKE-equation
  LOGICAL :: ltkeshs       = .TRUE.      ! consider separ. horiz. shear production in TKE-equation
  LOGICAL :: ltkenst       = .TRUE.      ! consider produc. by near-surf. thermals in TKE-equation

  LOGICAL :: loutsso       = .TRUE.      ! consider mechanical SSO-wake production of TKE for output
  LOGICAL :: loutshs       = .TRUE.      ! consider separ. horiz. shear production of TKE for output
  LOGICAL :: loutnst       = .FALSE.     ! consider produc. by near-surf. thermals of TKE for output
  LOGICAL :: loutbms       = .FALSE.     ! consider TKE-production by turbulent buoyancy, total mechanical shear
                                         !  or grid-scale mechanical shear for additional output

  LOGICAL :: ltmpcor       = .FALSE.     ! consideration minor turbulent sources in the enthalpy budget
  LOGICAL :: lcpfluc       = .FALSE.     ! consideration of fluctuations of the heat capacity of air

  LOGICAL :: lexpcor       = .FALSE.     ! explicit warm-cloud correct. of implicitly calculated turbul. diff.
  LOGICAL :: lsflcnd       = .TRUE.      ! lower flux condition for vertical diffusion calculation
  LOGICAL :: lcirflx       = .FALSE.     ! consideration of non-turbulent fluxes related to near-surface circulations
  LOGICAL :: ldiff_qi      = .FALSE.     ! turbulent diffusion of cloud ice QI acitve
  LOGICAL :: ldiff_qs      = .FALSE.     ! turbulent diffusion of snow      QS acitve
  LOGICAL :: lfreeslip     = .FALSE.     ! free-slip lower boundary condition (use for idealized runs only!)

! Notice that the following switches are provided by the parameter-list of
! SUB 'turbdiff' or 'turbtran':

! lnsfdia                   :calculation of (synoptical) near-surface variables required
! lsfluse                   :use explicit heat flux densities at the suface
! ltkeinp                   :TKE present as input (at level k=ke1 for current time level 'ntur')
! lgz0inp                   :gz0 present as input

!===================================================================================
! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! -----------------------------------------------------------------------------

  INTEGER :: imode_tran    = 0           ! mode of TKE-equation in transfer   scheme            (compare 'imode_turb')
  INTEGER :: imode_turb    = 1           ! mode of TKE-equation in turbulence scheme
                           ! 0: diagnostic equation
                           ! 1: prognostic equation (default)
                           ! 2: prognostic equation (implicitly positive definit)
  INTEGER :: icldm_tran    = 2           ! mode of cloud representation in transfer   parametr. (compare 'icldm_turb')
  INTEGER :: icldm_turb    = 2           ! mode of cloud representation in turbulence parametr.
                           !-1: ignoring cloud water completely (pure dry scheme)
                           ! 0: no clouds considered (all cloud water is evaporated)
                           ! 1: only grid scale condensation possible
                           ! 2: also sub grid (turbulent) condensation considered
  INTEGER :: itype_wcld    = 2           ! type of water cloud diagnosis within the turbulence scheme:
                           ! 1: employing a scheme based on relative humitidy
                           ! 2: employing a statistical saturation adjustment
  INTEGER :: itype_sher    = 0    ! type of mean shear-production for TKE
                           ! 0: only vertical shear of horizontal wind
                           ! 1: previous plus horizontal shear correction
                           ! 2: previous plus shear from vertical velocity

! These are the settings for the ICON-like setup of the physics:
  INTEGER :: imode_stbcalc = 1           ! mode of calculating the stability function (related to 'stbsecu')
                        ! (-)1: always for unstable strat. using a restr. 'gama' in terms of prev. forc.
                        ! (-)2: only to avoid non-physic. solution or if current 'gama' is too large
                        ! negative values for additional preconditioning
  INTEGER :: ilow_def_cond = 2           ! type of the default condition at the lower boundary
                           ! 1: zero surface gradient
                           ! 2: zero surface value

  INTEGER :: imode_pat_len = 2           ! mode of determining a length scale of surface patterns used for the "circulation-term"
                                         !  and additional roughness by tile-variation of land-use:
                           ! 1: employing the constant value 'pat_len' only
                           !    - raw "circulation term" considered as to be due to thermal surface-patterns.
                           ! 2: using the standard deviat. of SGS orography as a lower limit
                           !    - raw "circulation term" considered as to be due to thermal SSO effect

  !Note:
  !The theoretical background of the "circulation term" has meanwhile been fundamentally revised.
  !Accordingly, it is going to be substituted by two complementary approaches:
  !i) a thermal SSO parameterization and ii) a new "circulation term" due to thermal surface patterns.
  !Hence in the following, the still active raw parameterization is referred to as raw "circulation term".

  INTEGER :: imode_shshear = 2           ! mode of calculat. the separated horizontal shear (related to 'ltkeshs', 'a_hshr')
                           ! 0: with a constant lenght scale and based on 3D-shear and incompressibility
                           ! 1: with a constant lenght scale and considering the trace constraint for the 2D-strain tensor
                           ! 2: with a Ri-number depend. length sclale correct. and the trace constraint for the 2D-strain tensor
  INTEGER :: imode_tkesso  = 1           ! mode of calculat. the SSO source term for TKE production (related to 'ltkesso')
                           ! 1: original implementation
                           ! 2: with a Ri-dependent reduction factor for Ri>1
                           ! 3: as "2", but additional reduction for mesh sizes < 2 km
  INTEGER :: imode_snowsmot= 1           ! mode to treating the aerodynamic surface-smoothing by snow
                           ! 0: no smoothing active at all
                           ! 1: no impact on SAI, but full smoothing of z0 (G. Zaengl's approach)
                           ! 2: "1", but with full smoothing of SAI: full smoothing of z0 and SAI
                           ! 3: dynamical smoothing of z0 and SAI dependent on snow- and roughness height

  INTEGER :: itype_2m_diag = 1           ! type of 2m-diagnostics for temperature and -dewpoint (related to 'z0m_dia')
                           ! 1: Considering a fictive surface roughness of a SYNOP lawn
                           ! 2: Considering the mean surface roughness of a grid box
                           !    and using an exponential roughness layer profile
  INTEGER :: imode_stadlim = 2           ! mode of limitting statist. saturation adjustment in SUB 'turb_cloud'
                                         ! (related to 'q_crit, clc_diag')
                           ! 1: only absolut upper limit of stand. dev. of local super-satur. (sdsd)
                           ! 2: relative limit of sdsd and upper limit of cloud-water
  INTEGER :: imode_trancnf = 2           ! mode of configuring the transfer-scheme (SUB 'turbtran')
                           ![1: eliminated (old version)]
                           ! 2: 1-st ConSAT: start. with estim. Ustar, without a laminar correct. for prof.-funct.;
                           !    interpol. Tet_l onto "0"-level; calcul. Tet_l-gradients directly;
                           !    without an upper bound for TKE-forcing; with transmit. skin-layer depth to turbul.
                           ! 3: 2-nd ConSAT: as "2", but with a hyperbolic interpol. of profile function
                           !    for stable stratification
                           ! 4: 3-rd ConSAT: as "3", but without using an upper interpolation node
  INTEGER :: imode_lamdiff = 1           ! mode of considering laminar diffusion within the surface layer
                           ! 1: only a limitation of "0"-level DCs at calculating resistnance-lenght values for the R-layer
                           ! 2: laminar limit of full transfer-resistance values
  INTEGER :: imode_tkemini = 1           ! mode of adapting q=SQRT(2*TKE) and the TMod. to Lower Limits for Diff. Coeffs. (LLDCs)
                                         ! (related to 'tk[h|m]min*')
                           ! 1: LLDC treated as corrections of stability length without any further adaptation
                           ! 2: TKE adapted to that part of LLDC representing so far missing shear forcing, while the
                           !     assumed part of LLDC representing missing drag-forces has no feedback to the TMod.
                           ! 3: Tuned variant of "2" that is suitable for operational forecasts
  INTEGER :: imode_suradap = 0           ! mode of adapting the Diff. Coeff. (DC) at Half-Level (HL) "P" (k=ke) or "0" (k=ke1) as input
                                         !  of the surface-layer profile-function between levels "0" and "P" (related to 'tk[h|m]min*')
                           ! 0: no adaptations at all
                           ! 1: removing the artific. drag  contribut. by Lower Limits of DCs (LLDCs) for momentum at "P"-level
                           ! 2: "1" and including not transmitted shear contributions by LLDCs and NTCs to DCs at level "0"
  INTEGER :: imode_tkediff = 2           ! mode of implicit TKE-Diffusion (related to 'c_diff')
                           ! 1: in terms of q=SQRT(2*TKE)
                           ! 2: in terms of TKE=0.5*q**2
  INTEGER :: imode_adshear = 2           ! mode of considering addit. shear by scale interaction
                                         !  (related to 'ltkesso, ltkeshs, ltkecon, ltkenst')
                           ! 1: not  considered for stability functions
                           ! 2: also considered for stability functions


!#ifdef new_vers
REAL (KIND=wp)     ::        &
  ! do we need it as TARGET?
  ! these variables are set in SUB 'turb_setup'
  c_tke,tet_g,rim, &
  c_m,c_h, b_m,b_h,  sm_0, sh_0, &
  d_0,d_1,d_2,d_3,d_4,d_5,d_6,   &
  a_3,a_5,a_6,                   &

  ! these parameters are physical constants which are either taken as
  ! they are or set to 0.0 in the turbulence for special applications
  tur_rcpv, & ! cp_v/cp_d - 1
  tur_rcpl    ! cp_l/cp_d - 1 (where cp_l=cv_l)

! Definition of used data types
! -----------------------------

TYPE modvar !model variable
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
             av(:,:) => NULL(), & !atmospheric values
             sv(:)   => NULL(), & !surface     values (concentration or flux density)
             at(:,:) => NULL()    !atmospheric time tendencies
     LOGICAL                                 ::         &
             fc                   !surface values are flux densities
     INTEGER                                 ::         &
             kstart  = 1          !start level for vertical diffusion
END TYPE modvar

TYPE turvar !turbulence variables
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
             tkv(:,:) => NULL(), & !turbulent coefficient for vert. diff.
             tsv(:)   => NULL()    !turbulent velocity at the surface
END TYPE turvar

TYPE varprf !variable profile
     REAL (KIND=wp), POINTER, CONTIGUOUS     ::         &
             bl(:,:), & !variable at boundary model levels
             ml(:,:)    !variable at main     model levels
END TYPE varprf
!#endif


TYPE :: t_turbdiff_config !configuration state of turbulence mode

! 1. Numerical parameters:
! ------------------------

  REAL(wp):: impl_s        ! implicit weight near the surface (maximal value)
  REAL(wp):: impl_t        ! implicit weight near top of the atmosphere (maximal value)

  INTEGER :: imode_tkvmini ! mode of calculating the minimal turbulent diff. coeffecients
                                           ! Minimal diffusion coefficients in [m^2/s]:
  REAL(wp):: tkhmin        !  for vertical scalar (heat) transport
  REAL(wp):: tkmmin        !  for vertical momentum transport
  REAL(wp):: tkhmin_strat  ! enhanced 'tkhmin' for stratosphere
  REAL(wp):: tkmmin_strat  ! enhanced 'tkmmin' for stratosphere

  REAL(wp):: ditsmot       ! smoothing factor for direct time-step iterations

  INTEGER :: imode_frcsmot ! mode of applying vertical smoothing of TKE forcing
  REAL(wp):: frcsmot       ! vertical smoothing factor for TKE forcing
  REAL(wp):: tkesmot       ! time smoothing factor for TKE and diffusion coefficients

  REAL(wp):: frcsecu       ! security factor for TKE-forcing       (<=1)
  REAL(wp):: tkesecu       ! security factor in  TKE equation      (out of [0; 1])
  REAL(wp):: stbsecu       ! security factor in stability function (out of ]0; 1])
  REAL(wp):: prfsecu       ! relat. secur. fact. for prof. funct.  (out of ]0; 1[)

  REAL(wp):: epsi          ! relative limit of accuracy for comparison of numbers

  INTEGER :: it_end        ! number of initialization iterations (>=0)

! 2. Parameters describing physical properties of the lower boundary of the atmosphere:
! -------------------------------------------------------------------------------------

  REAL(wp):: rlam_heat     ! scaling factor of the laminar boundary layer for heat
  REAL(wp):: rlam_mom      ! scaling factor of the laminar boundary layer for momentum

  REAL(wp):: rat_lam       ! vapour/heat ratio of laminar scaling factors (over land)
  REAL(wp):: rat_sea       ! sea/land ratio of laminar scaling factors for heat (and vapor)
  REAL(wp):: rat_glac      ! glacier/land ratio of laminar scaling factors for heat (and vapor)

  REAL(wp):: rat_can       ! ratio of canopy height over sai*z0m

  INTEGER :: imode_nsf_wind! mode of local wind-definition at near-surface levels (related to 'rsur_sher')
  REAL(wp):: rsur_sher     ! scaling factor for the considered fraction of additional surface shear-forcing

  INTEGER :: imode_charpar ! mode of estimating the Charnock-Parameter
  REAL(wp):: alpha0        ! Charnock-parameter
  REAL(wp):: alpha0_max    ! upper limit of velocity-dependent Charnock-parameter
  REAL(wp):: alpha0_pert   ! additive ensemble perturbation of Charnock-parameter

  REAL(wp):: alpha1        ! parameter scaling the molecular roughness of water waves

! 3. Parameters that should be external parameter fields being not yet available:
! -------------------------------------------------------------------------------

                                           ! Surface area density:
  REAL(wp):: c_lnd         !  of the roughness elements over land
  REAL(wp):: c_sea         !  of the waves over sea
  REAL(wp):: c_soil        !  of the (evaporative) soil surface
  REAL(wp):: c_stm         !  of stems and branches at the plant-covered part of the surface
  REAL(wp):: e_surf        ! exponent to get the effective surface area

  LOGICAL :: lconst_z0     ! horizontally homogeneous roughness length (for idealized testcases!) applied
  REAL(wp)::  const_z0     ! horizontally homogeneous roughness length (for idealized testcases!)

! 4. Parameters that should be dynamical fields being not yet available:
! ----------------------------------------------------------------------

  REAL(wp):: z0m_dia       ! roughness length of a typical synoptic station [m]
  REAL(wp):: z0_ice        ! roughness length of sea ice

! 5. Parameters for modelling turbulent diffusion:
! ------------------------------------------------

  REAL(wp):: tur_len       ! asymptotic maximal turbulent distance [m]
  REAL(wp):: pat_len       ! effective global length scale of subscale surface patterns over land [m]
  REAL(wp):: len_min       ! minimal turbulent length scale [m]

  INTEGER :: imode_vel_min ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
  REAL(wp):: vel_min       ! minimal velocity scale [m/s]
  REAL(wp):: vel_max       ! maximal velocity scale [m/s]

  REAL(wp):: akt           ! von-Karman constant

                           ! Length scale factors:
                           !  for pressure destruction of turbulent:
  REAL(wp):: a_heat        !   scalar (heat) transport
  REAL(wp):: a_mom         !   momentum transport
                           !  for dissipation of turbulent:
  REAL(wp):: d_heat        !   scalar (temperature) variance
  REAL(wp):: d_mom         !   momentum variance
  REAL(wp):: c_diff        !  for turbulent transport (vertical diffusion) of TKE
  REAL(wp):: a_stab        !  for the stability correction applied to integral turbulent length-scale
  REAL(wp):: a_hshr        !  for separate horizontal shear circulations (related to 'ltkeshs')


  REAL(wp):: clc_diag      ! cloud cover at saturation
  REAL(wp):: q_crit        ! critical value for normalized super-saturation
  REAL(wp):: c_scld        ! shape-factor (0<=c_scld) applied to pure 'cl_cv' at the moist correct.

  LOGICAL :: ltkesso       ! consider mechanical SSO-wake production in TKE-equation
  LOGICAL :: ltkecon       ! consider convective buoyancy production in TKE-equation
  LOGICAL :: ltkeshs       ! consider separ. horiz. shear production in TKE-equation
  LOGICAL :: ltkenst       ! consider produc. by near-surf. thermals in TKE-equation

  LOGICAL :: loutsso       ! consider mechanical SSO-wake production of TKE for output
  LOGICAL :: loutshs       ! consider separ. horiz. shear production of TKE for output
  LOGICAL :: loutnst       ! consider produc. by near-surf. thermals of TKE for output
  LOGICAL :: loutbms       ! consider TKE-production by turbulent buoyancy, total mechanical shear
  LOGICAL :: ltmpcor       ! consideration minor turbulent sources in the enthalpy budget
  LOGICAL :: lcpfluc       ! consideration of fluctuations of the heat capacity of air

  LOGICAL :: lexpcor       ! explicit warm-cloud correct. of implicitly calculated turbul. diff.
  LOGICAL :: lsflcnd       ! lower flux condition for vertical diffusion calculation
  LOGICAL :: lcirflx       ! consideration of non-turbulent fluxes related to near-surface circulations
  LOGICAL :: ldiff_qi      ! turbulent diffusion of cloud ice QI acitve
  LOGICAL :: ldiff_qs      ! turbulent diffusion of snow      QS acitve
  LOGICAL :: lfreeslip     ! free-slip lower boundary condition (use for idealized runs only!)

  INTEGER :: imode_tran    ! mode of TKE-equation in transfer   scheme            (compare 'imode_turb')
  INTEGER :: imode_turb    ! mode of TKE-equation in turbulence scheme
  INTEGER :: icldm_tran    ! mode of cloud representation in transfer   parametr. (compare 'icldm_turb')
  INTEGER :: icldm_turb    ! mode of cloud representation in turbulence parametr.
  INTEGER :: itype_wcld    ! type of water cloud diagnosis within the turbulence scheme:
  INTEGER :: itype_sher    ! type of mean shear-production for TKE

  INTEGER :: imode_stbcalc ! mode of calculating the stability function (related to 'stbsecu')
  INTEGER :: ilow_def_cond ! type of the default condition at the lower boundary
  INTEGER :: imode_pat_len ! mode of determining a length scale of surface patterns used for the "circulation-term"
  INTEGER :: imode_shshear ! mode of calculat. the separated horizontal shear (related to 'ltkeshs', 'a_hshr')
  INTEGER :: imode_tkesso  ! mode of calculat. the SSO source term for TKE production (related to 'ltkesso')
  INTEGER :: imode_snowsmot! mode to treating the aerodynamic surface-smoothing by snow
  INTEGER :: itype_2m_diag ! type of 2m-diagnostics for temperature and -dewpoint (related to 'z0m_dia')
  INTEGER :: imode_stadlim ! mode of limitting statist. saturation adjustment in SUB 'turb_cloud'
  INTEGER :: imode_trancnf ! mode of configuring the transfer-scheme (SUB 'turbtran')
  INTEGER :: imode_lamdiff ! mode of considering laminar diffusion within the surface layer
  INTEGER :: imode_tkemini ! mode of adapting q=SQRT(2*TKE) and the TMod. to Lower Limits for Diff. Coeffs. (LLDCs)
  INTEGER :: imode_suradap ! mode of adapting the Diff. Coeff. (DC) at Half-Level (HL) "P" (k=ke) or "0" (k=ke1) as input
  INTEGER :: imode_tkediff ! mode of implicit TKE-Diffusion (related to 'c_diff')
  INTEGER :: imode_adshear ! mode of considering addit. shear by scale interaction

  ! Note:
  ! This very list of configure-variables appears 3 times, since each variable is
  ! - declared and initialized by defaults individually
  ! - declared as associated component of this configure-state data structre
  ! - included into the loading of the config-vector 'turbdiff_config(jg)' via SUB 'load_turbdiff_config'.
  ! Both, the individual variables and SUB 'load_turbdiff_config', are also USEed in 'mo_turbdiff_nml',
  !  which saves another 2 listings of these variables there.
  ! For that reason, the individual variables are still required, and pure indirect initialization is not sufficient.
  ! The configure-state allows to USE just all these configure-variables into the modelling routines without
  !  the need to either mention each variable individually or to make the whole configure module accessible.
  ! As 'turbdiff_config' is a configure-state vector, also domain-specific values can be applied.

! 6. control parameters being set at interface level (and neither by pure initialization nor by namelist):
! --------------------------------------------------------------------------------------------------------

  INTEGER :: iinit         ! indicator for initialization level

  REAL(wp), DIMENSION(:), POINTER :: impl_weight ! implicit weights for tridiagonal solver

END TYPE t_turbdiff_config

!===================================================================================
! Declarations of utility variables:

! Turbulence parameters which are computed during model run
! ---------------------------------------------------------

TYPE(t_turbdiff_config), TARGET  :: turbdiff_config(max_dom) ! 'turbdiff' configuration state-vector (for each domain)
!$ACC DECLARE CREATE(turbdiff_config)

!===================================================================================

CONTAINS

!===================================================================================

SUBROUTINE load_turbdiff_config (jg)

! Loads the full configuration-state except part 6. in the above TYPE-specificatin of 't_turbdiff_config'
!  (regarding "control parameters being set at interface level (and neither by pure initialization nor by namelist)")
! Thus, each component of 'turbdiff_config(jg)' except those of part 6. needs to be included here; and these are all
!  not strictly constant configuration-quantities (that means without a "PARAMETER"-attribute) being initialized in this
!  module and not set at interface level!

    INTEGER, INTENT(IN) :: jg !patch index

      turbdiff_config(jg)%impl_s         = impl_s
      turbdiff_config(jg)%impl_t         = impl_t

      turbdiff_config(jg)%imode_tkvmini  = imode_tkvmini
      turbdiff_config(jg)%tkhmin         = tkhmin
      turbdiff_config(jg)%tkmmin         = tkmmin
      turbdiff_config(jg)%tkhmin_strat   = tkhmin_strat
      turbdiff_config(jg)%tkmmin_strat   = tkmmin_strat

      turbdiff_config(jg)%ditsmot        = ditsmot

      turbdiff_config(jg)%imode_frcsmot  = imode_frcsmot
      turbdiff_config(jg)%frcsmot        = frcsmot
      turbdiff_config(jg)%tkesmot        = tkesmot

      turbdiff_config(jg)%frcsecu        = frcsecu
      turbdiff_config(jg)%tkesecu        = tkesecu
      turbdiff_config(jg)%stbsecu        = stbsecu
      turbdiff_config(jg)%prfsecu        = prfsecu

      turbdiff_config(jg)%epsi           = epsi

      turbdiff_config(jg)%it_end         = it_end

      turbdiff_config(jg)%rlam_heat      = rlam_heat
      turbdiff_config(jg)%rlam_mom       = rlam_mom

      turbdiff_config(jg)%rat_lam        = rat_lam
      turbdiff_config(jg)%rat_sea        = rat_sea
      turbdiff_config(jg)%rat_glac       = rat_glac

      turbdiff_config(jg)%rat_can        = rat_can

      turbdiff_config(jg)%imode_nsf_wind = imode_nsf_wind
      turbdiff_config(jg)%rsur_sher      = rsur_sher

      turbdiff_config(jg)%imode_charpar  = imode_charpar
      turbdiff_config(jg)%alpha0         = alpha0
      turbdiff_config(jg)%alpha0_max     = alpha0_max
      turbdiff_config(jg)%alpha0_pert    = alpha0_pert
      turbdiff_config(jg)%alpha1         = alpha1

      turbdiff_config(jg)%c_lnd          = c_lnd
      turbdiff_config(jg)%c_sea          = c_sea
      turbdiff_config(jg)%c_soil         = c_soil
      turbdiff_config(jg)%c_stm          = c_stm
      turbdiff_config(jg)%e_surf         = e_surf

      turbdiff_config(jg)%lconst_z0      = lconst_z0
      turbdiff_config(jg)%const_z0       = const_z0

      turbdiff_config(jg)%z0m_dia        = z0m_dia
      turbdiff_config(jg)%z0_ice         = z0_ice

      turbdiff_config(jg)%tur_len        = tur_len
      turbdiff_config(jg)%pat_len        = pat_len
      turbdiff_config(jg)%len_min        = len_min

      turbdiff_config(jg)%imode_vel_min  = imode_vel_min
      turbdiff_config(jg)%vel_min        = vel_min
      turbdiff_config(jg)%vel_max        = vel_max

      turbdiff_config(jg)%akt            = akt

      turbdiff_config(jg)%a_heat         = a_heat
      turbdiff_config(jg)%a_mom          = a_mom
      turbdiff_config(jg)%d_heat         = d_heat
      turbdiff_config(jg)%d_mom          = d_mom
      turbdiff_config(jg)%c_diff         = c_diff
      turbdiff_config(jg)%a_stab         = a_stab
      turbdiff_config(jg)%a_hshr         = a_hshr


      turbdiff_config(jg)%clc_diag       = clc_diag
      turbdiff_config(jg)%q_crit         = q_crit
      turbdiff_config(jg)%c_scld         = c_scld

      turbdiff_config(jg)%ltkesso        = ltkesso
      turbdiff_config(jg)%ltkecon        = ltkecon
      turbdiff_config(jg)%ltkeshs        = ltkeshs
      turbdiff_config(jg)%ltkenst        = ltkenst

      turbdiff_config(jg)%loutsso        = loutsso
      turbdiff_config(jg)%loutshs        = loutshs
      turbdiff_config(jg)%loutnst        = loutnst
      turbdiff_config(jg)%loutbms        = loutbms
      turbdiff_config(jg)%ltmpcor        = ltmpcor
      turbdiff_config(jg)%lcpfluc        = lcpfluc

      turbdiff_config(jg)%lexpcor        = lexpcor
      turbdiff_config(jg)%lsflcnd        = lsflcnd
      turbdiff_config(jg)%lcirflx        = lcirflx
      turbdiff_config(jg)%ldiff_qi       = ldiff_qi
      turbdiff_config(jg)%ldiff_qs       = ldiff_qs
      turbdiff_config(jg)%lfreeslip      = lfreeslip

      turbdiff_config(jg)%imode_tran     = imode_tran
      turbdiff_config(jg)%imode_turb     = imode_turb
      turbdiff_config(jg)%icldm_tran     = icldm_tran
      turbdiff_config(jg)%icldm_turb     = icldm_turb
      turbdiff_config(jg)%itype_wcld     = itype_wcld
      turbdiff_config(jg)%itype_sher     = itype_sher

      turbdiff_config(jg)%imode_stbcalc  = imode_stbcalc
      turbdiff_config(jg)%ilow_def_cond  = ilow_def_cond
      turbdiff_config(jg)%imode_pat_len  = imode_pat_len
      turbdiff_config(jg)%imode_shshear  = imode_shshear
      turbdiff_config(jg)%imode_tkesso   = imode_tkesso
      turbdiff_config(jg)%imode_snowsmot = imode_snowsmot
      turbdiff_config(jg)%itype_2m_diag  = itype_2m_diag
      turbdiff_config(jg)%imode_stadlim  = imode_stadlim
      turbdiff_config(jg)%imode_trancnf  = imode_trancnf
      turbdiff_config(jg)%imode_lamdiff  = imode_lamdiff
      turbdiff_config(jg)%imode_tkemini  = imode_tkemini
      turbdiff_config(jg)%imode_suradap  = imode_suradap
      turbdiff_config(jg)%imode_tkediff  = imode_tkediff
      turbdiff_config(jg)%imode_adshear  = imode_adshear

END SUBROUTINE load_turbdiff_config

!===================================================================================

END MODULE mo_turbdiff_config
