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
! Source module for computing diffusion coefficients, Turbulent Kinetic Energy (TKE),
! Eddy Dissipation Rate (EDR), some particular scale-interaction source terms
! for TKE and some optional turbulent tendencies apart from turbulent diffusion:
!
! Description of *turb_diffusion*:
!   This  module calculates the tendencies for turbulent
!   vertical transport of momentum and heat and the coefficients
!   for turbulent diffusion as well.
!
!   The clousure is made on lever 2.5 (Mellor/Yamada) using a prognostic
!   TKE-equation and includes the formulation of a flow through a porous
!   medium (roughness layer)
!
!   The turbulence model (with some Prandtl-layer approximations is used
!   for the calculation of turbulent transfer between atmosphere and the
!   lower boundary too.
!
! The module contains the public subroutines :
!
!   turbdiff
!
! called from the turbulence interface routine of the model.
!
!-------------------------------------------------------------------------------

MODULE turb_diffusion

!-------------------------------------------------------------------------------
!
! COSMO-Documentation History:
!
!  The history of all these modifications is as follows, where those belonging to the fomal
!   reorganization of the whole package (atmospheric turbulence and surface-to-atmpsphere transfer)
!   are now in the header of MODULE 'turb_utilities', containing various common SUBs for 'turbdiff'
!   and 'turtran' (related moist thermodynamicds and the treatment of turbulent budget equations)
!   and also the blocked code for semi-implicit vertical diffusion. The new blocked version of SUB 'turbtran'
!   is now in MODULE 'turb_transfer':
!
!              2010/12/17 Matthias Raschendorfer
!  Introduction of a TKE-source term due to scale interaction with sub grid scale convection
!   using 'ltkecon' and the convective buoyant heat flux density in 'tket_conv'.
!              2011/02/18 Matthias Raschendorfer
!  Introduction of some minor formal modif. of some parts of the code mainly to achiev better vectorization
!   (results will be modyfied only because of numerical effects).
!  Introduction of OPTIONAL extra output fields 'tket_sso' and 'tket_hshr'.
!              2011/03/23 Matthias Raschendorfer
!  Substitution of run time allocations because of bad performance on some computers.
!              2011/08/26 Matthias Raschendorfer
!  Discriminating the variable for standard deviation of super-saturation (input) and cloud fraction (output)
!   in the CALL of 'turb_cloud'.
!  Changing the the definition of LOGICAL 'lcircterm'.
!  Introduction of a preconditioning of the tridiagonal system for "imode_turb=3".
!              2011/09/23 Matthias Raschendorfer
!  Changing 'k' to 'ke' in the initialization section of 'turbdiff' in case of "itype_tran = 3".
!  Formulating temp-grad for 'ltmpcor' for "k = ke1" similar to the other levels.
!  Removing a wrong multiplication with Exner-factor in case of "imode_turb = 4".
!  Adding a missing multiplication by 'dicke()' in case of ".NOT.limpltkediff".
!  Dirscriminating between two effective values of effective Prandtl-layer depth,
!   by introducing the arrays 'dz0(:,mom|sca)' and 'vh0(:,:)
!   removing wrong surface conditions in case "imode_turb >= 3 .AND. itype_tran = 2".
!              2011/12/08 Matthias Raschendorfer
!   Diffusion tendency (rather than flux-density) of pot. temp. is converted in that of ordinary temp.
!              2012/01/10 Matthias Raschendorfer
!  Correction of some bugs:
!   In case of "limpltkediff" 'ko' needs to be added by "1".
!   Coding mistake in n-loop bounds for surface layer gradients in 'turbdiff'.
!  Introduction of array 'rhoh' (density on main levels) to avoid a back interpolation from half levels.
!  Using a mass weighted density interpolation onto half levels as for all other variables.
!              2012/01/11 Matthias Raschendorfer
!  Some rearrangement and simplification:
!   Removing explicit vertical diffusion and expressing explicit (eg. moist) corrections
!    by corrected vertical profiles still using implicit formulation of vertical diffusion.
!   Now "imode_turb=2" is for semi-implicit diffusion within 'turbdiff' using a concentration condition
!    at the surface and "imode_turb=3" is for the same with a flux condition at the surface.
!    "imode_turb_turb=4" is obsolet now.
!   One loop for diffusion of all first order model variables.
!              2012/01/20 Matthias Raschendorfer
!  Reformulation of roughness layer drag by an implicit equation (keeping sign of wind components)
!  Expressing explicit circulation term (as far as it has divergence form) by a respective
!   q**2-profile correction within the implicit formulation of diffusion similar to the treatment
!   of the former explicit (moist) corrections
!  Introduction of effective flux profile calculation into SUB 'calc_impl_vert_diff'.
!  Setup of a positiv definit solution avoiding some of the former limitations.
!              2012/03/20 Matthias Raschendorfer
!  Rearrangement, modularization and revision of numerical treatment, such as:
!  Intoduction of the SUBs 'solve_turb_budgets' and 'adjust_satur_equil' in order to use the same code
!  for the same purpose in 'turbtran' and 'turbdiff.
!  Removing the explicit diffusion option.
!  Introducing the new driving SUB 'vert_grad_diff' organizing vertical diffusion when called in a loop of
!  several full level variables and managing options how to treat non-gradient-fluxes and additional
!  explicit time tendencies numerically, including smoothing options and preconditioning the linear system.
!  Parameters 'imode_tran' and 'imode_turb' define whether the TKE-equation is solved in a dianostic (1)
!  or a prognostic (2) form.
!  New parameter 'lsflcnd' is a flag for using a lower flux condition for vertical diffusion.
!  Introducing the flags 'lturatm', 'ltursrf', 'lmomdif', 'lscadif' arranging what tasks of
!  'organize_turbdiff' are required.
!  Introduction of 'ltkeshs' and removing the case "itype_sher = 3".
!              2014/07/28 Matthias Raschendorfer
!  Introduction of precalculated 'hdef2', 'hdiv', 'dwdx', 'dwdy' used for calculation of horizontal shear
!   including the scale separated non-turbulent part controlled by 'itype_sher'.
!  Simpler (physically identical) formulation of moist flux conversion
!   -> only numerical differences in the case of 'lexplcor'
!  Numerically more efficient formulation of Blackadar 'len_scale'
!   -> only numerical differences
!  Eliminating array 'wind', as 'u' and 'v' are already defined at mass positions.
!              2015/08/25 Matthias Raschendorfer
! Adopting other development within the ICON-version (namely by Guenther Zaengl) as switchable options
!  related to the following new selectors and switches:
!   imode_pat_len, imode_frcsmot, imode_shshear, imode_tkvmini, imode_charpar.
! Rearranging the development by Matthias Raschendorfer that had not yet been transferred to COSMO as
!  switchable options related to the following switches:
!   lsflcnd, ldynimp, lprecnd, ltkeshs, loutshs
!  and selectors:
!   imode_tkediff, imode_adshear, imode_tkemini
!  and a partly new (more consistent) interpretation of:
!   imode_turb, icldm_turb and  itype_sher
! Controlling numerical restrictions gradually namely by the parameters:
!  tndsmot, frcsmot
! Correction of some bugs (of the latest ICON-version) in particular related to the optional lower
!  concentration condition.
! Using the arrays 'tvm', 'tvh' and 'tkm', allowing an easier formulation of transfer-resistances.
!              2016-05-10 Ulrich Schaettler
! Splitting this module from the original module 'organize_turbdiff' as it was used by ICON before.
! Moving declarations, allocation and deallocations of ausxilary arrays into MODULE 'turb_data'.
!-------------------------------------------------------------------------------------------------------


! Modules used:

#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num
#endif

!-------------------------------------------------------------------------------
! Parameter for precision
!-------------------------------------------------------------------------------

USE mo_kind,         ONLY : wp, vp
USE mo_exception,    ONLY: message_text, message

!-------------------------------------------------------------------------------
! Mathematical and physical constants
!-------------------------------------------------------------------------------

USE mo_mpi,                ONLY : get_my_global_mpi_id
USE mo_exception,          ONLY : finish
USE mo_physical_constants, ONLY : &
!
! Physical constants and related variables:
! -------------------------------------------
!
    r_d      => rd,       & ! gas constant for dry air
    rvd_m_o  => vtmpc1,   & ! r_v/r_d - 1
    cp_d     => cpd,      & ! specific heat for dry air
    lh_v     => alv,      & ! evaporation heat
    lhocp    => alvdcp,   & ! lh_v / cp_d
    con_m,                & ! kinematic vsicosity of dry air (m2/s)
con_h,                & ! scalar conductivity of dry air (m2/s)
grav                    ! acceleration due to gravity

!-------------------------------------------------------------------------------
! From Flake model
!-------------------------------------------------------------------------------

USE sfc_flake_data, ONLY: &
h_Ice_min_flk      ! Minimum ice thickness [m]

!-------------------------------------------------------------------------------
! Turbulence data (should be the same in ICON and COSMO)
!-------------------------------------------------------------------------------

USE mo_turbdiff_config, ONLY :          &
    t_turbdiff_config, &

! Numerical constants and parameters:
! -----------------------------------
    ! derived parameters calculated in 'turb_setup'
    tet_g, rim, b_m, b_h, sm_0, sh_0,   &
    a_3, a_5 ,a_6,                      &
    tur_rcpv, tur_rcpl,                 &

    ! used derived types
    modvar, turvar, varprf, & !

    ldynimp,      & ! dynamical calculation of implicit weights
    lprecnd,      & ! preconditioning of tridiagonal matrix

    ! numbers and indices

    nvel    ,     & ! number of velocity-components active for turbulece ('u_m', 'v_m')
    naux    ,     & ! number of auxilary variables
    ndim    ,     & ! (positive) limit of last dimension used for 'zaux' and 'zvari'
    nmvar   ,     & ! number of included dynamically active prognostic model variables

    ntyp    ,     & ! number of variable types ('mom' and 'sca')

    mom     ,     & ! index for a momentum variable
    sca     ,     & ! index for a scalar   variable
    u_m     ,     & ! index for mass centered zonal      velocity-compont
    v_m     ,     & ! index for mass centered meridional  ,,         ,,
    tet_l   ,     & ! index for liquid-water potential temperature
    tet     ,     & ! index for potential temperature
    tem     ,     & ! index for temperature
    h2o_g   ,     & ! index for toatal water
    vap     ,     & ! index for water vapor
    liq             ! index for liquid water

!Note: It always holds: "tem=tet=tet_l=tem_l" and "vap=h2o_g" (respective usage of equal indices)!

!-------------------------------------------------------------------------------
! Control parameters for the run
!-------------------------------------------------------------------------------

USE turb_utilities,          ONLY:   &
    turb_setup,                      &
    adjust_satur_equil,              &
    solve_turb_budgets,              &
    prep_impl_vert_diff,             &
    calc_impl_vert_diff,             &
    vert_smooth,                     &
    zbnd_val, bound_level_interp

!-------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
    lsclm, latmflu, i_cal, i_mod, imb, &
    SHF, LHF
#endif
!SCLM---------------------------------------------------------------------------

USE mo_fortran_tools, ONLY: set_acc_host_or_device

!===============================================================================

IMPLICIT NONE

PUBLIC  :: turbdiff

!===============================================================================

REAL (KIND=wp), PARAMETER :: &
    z0 = 0.0_wp,    &
    z1 = 1.0_wp,    &
    z2 = 2.0_wp,    &
    z3 = 3.0_wp,    &
    z4 = 4.0_wp,    &
    z5 = 5.0_wp,    &
    z6 = 6.0_wp,    &
    z7 = 7.0_wp,    &
    z8 = 8.0_wp,    &
    z9 = 9.0_wp,    &
    z10=10.0_wp,    &

    z1d2=z1/z2     ,&
    z1d3=z1/z3     ,&
    z2d3=z2/z3     ,&
    z3d2=z3/z2

!===============================================================================

CONTAINS

!===============================================================================

#  define err_args

SUBROUTINE turbdiff ( &
!
          tdc,                                                       &
!
          iini, ltkeinp, l3dturb,                                    &
                lrunsso, lruncnv, lrunscm,                           &
                ldoexpcor, ldocirflx,                                &
!
          dt_var,dt_tke, nprv, ntur, ntim,                           &
!
          nvec, ke, ke1, kcm, iblock, ivstart, ivend,                &
!
!
          l_hori, hhl, dp0, trop_mask, innertrop_mask,               &
!
          gz0, l_pat, c_big, c_sml, r_air,                           &
!
          t_g, qv_s, ps,                                             &
          u, v, w, t, qv, qc, prs, rhoh, rhon, epr,                  &
!
          tvm, tvh, tfm, tfh, tfv, tkred_sfc, tkred_sfc_h,           &
          tke, tkvm, tkvh, tprn, rcld, tkhm, tkhh,                   &
          hdef2, hdiv, dwdx, dwdy,                                   &
!
          edr, tket_sso, tket_nstc, tket_conv, tket_hshr,            &
          tket_buoy, tket_fshr, tket_gshr,                           &
          u_tens, v_tens, t_tens,                                    &
          tketens, tketadv,                                          &
          ut_sso, vt_sso,                                            &
!
          zvari                                                      &
!
          err_args)

!-------------------------------------------------------------------------------
!
! Notes:
!
! All tendency parameters are OPTIONAL (except 'tketens'. If they are missing,
!  calculated tendencies of SUB 'turbdiff' are automatically added to the related
!  prognostic variables.
!
! It is also possible to use only one time level for TKE using "ntim=1" and thus "nprv=1=ntur".
!
! Description:
!
!     Es werden die Diffusionskoeffizienten berechnet und ggf. Anteile
!     der zeitlichen Tendenzen der turbulenten Diffusion bestimmt
!     und zu den Tendenzfeldern hinzuaddiert.
!     Optional wird eine explizite oder (teil-)implizite Berechnung der
!     Diffusionstendenzen oder aber nur eine Berechnung der Diffusions-
!     koeffizienten durchgefuehrt. Im letzten Fall wird dann ein
!     implizit zu berechnender Anteil der Diffusionstendenzen an
!     anderer Stelle (slow_tendencies) bestimmt.
!     Allerdings koennen dann zusaetzliche explizite Korrekturtendenzen
!     hier in tubdiff bestimmt werden.
!
! Method:
!
!     Die Berechnung basiert auf einer Schliessung 2-ter Ordnung auf
!     dem level 2.5 (nach Mellor/Yamada). Demnach wird also eine
!     prognostische Gleichung fuer die TKE geloest.
!     Ausser der TKE-Advektion, die zusammen mit den Advektionstendenzen
!     der anderen prognostischen Variablen an anderer Stelle berechnet
!     wird, geschieht die gesamte TKE-Prognose in diesem Unterprogramm.

!     Die Formulierung des Schemas erfolgt mit thermodynamischen
!     Variablen, die bei feuchtadiabatischen vertikalen Verrueckungen
!     erhalten bleiben (pot. Fluessigw.temp. und Gesamtwassergehalt),
!     so dass der Kondesationseffekt auf subskalige Vertikalbewegungen
!     beruecksichtigt wird.
!     Die turbulenten Flussdichten der Erhaltungsgroessen werden in
!     solche der Modellvariablen konvertiert, so dass die thermodyn.
!     Kopplung der Flussdichten richtig erhalten bleibt.

!     Angeschlossen ist auch ein optionales statistisches Wolkenschema
!     (nach Sommeria und Deardorff), SUB 'turb_cloud', welches auch
!     subskalige Bewoelkung mit Hilfe der ueber das Feld 'rcld' ausge-
!     gebenen Standardabweichung des Saettigungsdefizites (SDSS) berechnet.

!     Das Turbulenzschema wurde so verallgemeinert, dass es auch bei
!     einer vertikal aufgeloesten Bestandesschicht gueltig ist, indem
!     idealisierend von der Durchstroemung eines poroesen Mediums
!     ausgegangen wird. Die Bilanzgleichungen 1-ter und 2-ter Ordnung
!     enthalten dann zusaetzliche Terme, welche die Wechselwirkungen
!     mit der Bestandes-Matrix beschreiben. Dies wirkt sich zum einen
!     auf die Stabilitaetsfunktionen und zum anderen vor allem auf die
!     TKE-Gleichung aus, welche einen auf den Formwiderstand der
!     Bestandeselemente zurueckzufuehrenden zusaetzlichen Quellterm
!     (Nachlaufturbulenz) enthaelt. Ausserdem werden die turbulenten
!     Flussdichtedivergenzen noch um einen Zusatzterm, welcher der
!     Reduktion des lufterfuellten Volumens im Gitterelement Rechnung
!     traegt, erweitert. Der Effekt des Formwiderstandes in der
!     Impulsgleichung ist ebenfalls beruecksichtigt. Die zusaetzlichen
!     Tendenzterme, die auf die Flussdichten zwischen Bestandes-Matrix
!     und umgebender Luft zurueckzufuehren sind (Bestandesquellen),
!     muessen noch in einem separaten Bestandesmodell parametrisiert
!     werden und sind nicht Gegenstand des Turbulenzmodells.
!     Schliesslich wird auch der Effekt der Transformation von Turbulenz
!     auf der dominierenden Skala in kleinskalige dissipative Turbulenz
!     durch Wirbelbrechen an Koerpern mit
!          Laengenskalen der Abmessungen << turbulente Laengenskala
!     beruecksichtigt; was sich durch eine (von der Laengenskala und
!     Volumendichte jener sehr kleinen Bestandeselemente aubhaengige)
!     Modifikation der Modellkonstanten ausdruecken laesst.

!     Es wird auch versucht, den Effekt thermisch induzierter
!     Zirkulationen auf die TKE-Produktion zu beruecksichtigen.
!     Hierdurch wird (vor allem) der Austauch in der naechtlichen
!     Grenzschicht erhoeht, was der Tendenz des alten Schemas,
!     in Bodennaehe zu kalte und nicht schnell genug anwachsende
!     Inversionen zu produzieren, entgegenwirkt.

!     Optional kann die Berechnung der vertikalen Gradienten durch eine
!     nicht-lokale Variante erfolgen. Hierbei werden die Gradienten
!     mit Profilen gebildet, die mit einem ueber die stabilitaets-
!     abhaengige Laengenskala gebildeten gleitenden Mittel behandelt
!     wurden.
!
!     Die Bildung der Anteile der Diffusionstendenzen, die numerisch
!     durch die Multiplikation mit einer Tridiagonalmatrix ausdrueckbar
!     sind, kann (neben der expliziten Variante) auch implizit erfolgen
!     (im Falle der Berechnung von nich-lokalen Gradienten und fuer die
!     TKE allerdings nur explizit).
!     Bei expliziter Rechnung ist, um auch bei Zeitschritten von
!     mehreren Minuten numerisch stabil zu bleiben, eine Limitierung
!     der Groesse der Diffusionskoeffezienten und eine im vertikalen Integral
!     quellenfreie numerische Glaettung der Vertikalprofile der Diffusions-
!     tendenzen erforderlich, sowie eine teilimplizite Behandlung der
!     Diffusionstendenzen in der untersten Modellschicht, erforderlich.
!
!     Die unteren Randwerte der turbulenten Flussdichten werden ueber
!     die Transferkoeffizienten zwischen Erdboden und unterster
!     Modellschicht (tcm und tch) bestimmt.
!     Optional koennen die Transferkoeffizienten auch mit diesem
!     Unterprogramm bestimmt werden, indem das Turbulenzmodell auch auf
!     das Niveau d+z0 angewandt wird, wobei vertikale Gradienten in
!     diesem Niveau mit Hilfe der Prandtl-Schicht-Hypothese berechnet
!     werden.
!     In diesem Zusammenhang wird auch die Wirkung der laminaren
!     Grenzschicht behandelt.
!
!     Turbulente Horizontaldiffusion (um ein 3-d-Schema zu erhalten)
!     ist noch nicht enthalten, kann aber integriert werden.
!     Uebergabevariablen:

!-------------------------------------------------------------------------------

! Declarations
!-------------------------------------------------------------------------------

!Formal Parameters:
!-------------------------------------------------------------------------------

! Parameters controlling the call of 'turbtran':
! ----------------------------------------------

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

LOGICAL, INTENT(IN) :: &
  l3dturb,      & !a model run with 3D-(turbulent)-diffusion

  ltkeinp,      & !TKE present as input for current time level 'ntur'

  lruncnv,      & !convection scheme is active
  lrunsso,      & !SSO-Scheme is active
  lrunscm         !a Single Column run (default: FALSE)

LOGICAL, OPTIONAL, INTENT(OUT) :: & !not necessary for initialization
  ldoexpcor,    & !consider explicit warm-cloud correct. for turb. scalar fluxes
  ldocirflx       !consider circulation heat-flux

REAL (KIND=wp), INTENT(IN) :: &
  dt_var,       & !time step for ordinary prognostic variables            ( s )
  dt_tke          !time step for the 2-nd order porgnostic variable 'tke' ( s )

INTEGER,        INTENT(IN) :: &
                  !                              -1: no, but precalculation of constant parameters
  iini,         & !type of initialization (0: no, 1: separate before the time loop
                  !                             , 2: within the first time step)
  ntur,         & !current  time level of 'tke' valid after  prognostic incrementation
  nprv,         & !previous time level of 'tke valid before prognostic incrementation
  ntim            !number of 'tke' time levels

! Horizontal and vertical sizes of the fields and related variables:
! ------------------------------------------------------------------

INTEGER,        INTENT(IN) :: &
  nvec,         & ! number of horizontal grid points in the nproma-vector
  ke,           & ! index of the lowest main level
  ke1,          & ! index of the lowest half level (=ke+1)
  !Note:
  !Level 'ke1' represents the top of that part of the grid-scale R-layer that is generated by land-use obstacles.
  kcm,          & ! level index of the upper canopy bound
  iblock

! Start- and end-indices for the computations in the horizontal layers:
! ---------------------------------------------------------------------

INTEGER,        INTENT(IN) :: &
  ivstart,      & ! start index in the nproma vector
  ivend           ! end index in the nproma vector

! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! -----------------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: &
  hhl,          & ! height of model half(=boundary)-levels        ( m )
  dp0             ! pressure thickness of layer                   (pa )

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
  l_pat,        & ! effektive Laengenskala der therm. Inhomogenitaeten
                  ! der Erdbodenoberflaeche
  l_hori          ! horizontal grid spacing (m)

REAL (KIND=wp), DIMENSION(:,kcm:), TARGET, OPTIONAL, INTENT(IN) :: &
  c_big,        & ! effective drag coefficient of canopy elements
                  ! larger than or equal to the turbulent length scale (1/m)
  c_sml           ! effective drag coefficient of canopy elements
                  ! smaller than the turbulent length scale            (1/m)

REAL (KIND=wp), DIMENSION(:,kcm-1:), TARGET, OPTIONAL, INTENT(IN) :: &
  r_air           ! log of air containing fraction of a gridbox inside
                  ! the canopy                                          (1)

! Fields for surface values and soil|canopy model variables:
! ------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
  ps,           & ! surface pressure                              ( pa  )
  qv_s,         & ! specific water vapor content on the surface   (kg/kg)
  t_g             ! weighted surface temperature                  (  k  )

! Atmospheric model variables:
! ---------------------------------

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
                  ! main-level values of:
  u,            & ! zonal wind speed       (at mass positions)    ( m/s )
  v,            & ! meridional wind speed  (at mass positions)    ( m/s )
  t               ! temperature                                   (  k  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
                  ! main-level values of:
  qv,           & ! specific water vapor content                  (kg/kg)
  qc              ! specific cloud water content                  (kg/kg)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
  prs             ! atmospheric pressure (at main levels)         ( pa  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
                  ! main-level values of:
  rhoh,         & ! total density of air (at main levels)         (kg/m3)
  epr             ! exner pressure                                ( --- )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(OUT) :: &
  rhon            ! total density of air (at half levels)         (kg/m3)

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
  w               ! vertical wind speed (defined at half levels)  ( m/s )

! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
  gz0,           & ! roughness length * g of the vertically not resovled R-layer
                   ! valid for the grid-scale surface              (m2/s2)
  !Achtung: Der g-Faktor ist ueberfluessig!

  ! turbulent (transfer) velocity scales at the surface
  tvm,           & ! for momentum                                  ( m/s)
  tvh              ! for heat and moisture                         ( m/s)

  !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
  !vecolities  'tvm' and 'tvh' makes live much easier!!

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
  tfm,           & !
  ! Prandtl-layer fraction of total transfer-layer resistance      ( --- ) for momentum as INP
  ! factor for removing the pure drag-contribution of 'tkmmin'     ( --- ) for momentum as OUT
  tfh,           & !
  ! Prandtl-layer fraction of total transfer-layer resistance      ( --- ) for scalars  as INP
  ! additional shear-forcing corresponding to the impact of LLDCs  ( 1/s2) at "P"-level as OUT
  tfv              ! additional shear-forcing by NTCs              ( 1/s2) at "P"-level as OUT

  !Attention: "INTENT(OUT)" might cause not-intended default-settings for 'tfv' in case of "rsur_sher > 0"!

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
  tkred_sfc, tkred_sfc_h   ! reduction factors for minimum diffusion coefficients near the surface


! Atmospheric variables of the turbulence model:
! ------------------------------------------------

 REAL (KIND=wp), DIMENSION(nvec,ke1,ntim), TARGET, INTENT(INOUT) :: &
                   ! half-level values of:
  tke              ! q:=SQRT(2*TKE) with TKE='turbul. kin. energy' ( m/s )
  !Note:
  !'tke' is the "turbulent velocity" (in m/s) and NOT the (mass-density) of turb. kin. energy,
  ! which has the dimension m2/s2!
  !In case of "ntim=1", the actual parameter for 'tke' may be a 2-dim. array for a fix time level.

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
                   ! half-level values of:
  tkvm,          & ! turbulent diffusion coefficient for momentum  (m2/s )
  tkvh             ! turb. diff. coeff. for heat and other scalars (m2/s )
                   ! (both defined at half levels)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT)    :: &
  tprn             ! turbulent Prandtl-number (at half-levels)     ( --- )

REAL (KIND=wp), DIMENSION(:,:,0:), TARGET, INTENT(OUT) :: &
  zvari            ! 3-rd dim. > 0: (quasi-conserved) model variables at main levels (including the lower boundary)
                   !                or their (possibly non-local) vertical gradients (at half levels)
                   !                final output: converted effective vertical gradients of regular model variables
                   ! 3-rd dim. = 0: potentially available energy per volume (pressure) at half levels or
                   !                as output: effective vertical gradient of availabel Circulation Kinetic Energy
                   !                (per mass) due to near surface thermal inhomogeneity (CKE),
                   !                which is a related acceleration

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
  rcld             ! standard deviation of local super-saturation (SDSS)
                   !  at MAIN levels including the lower boundary  (---)
                   ! AUX: cloud-cover at main levels (as output of SUB 'adjust_satur_equil'
                   !        and later at half levels (as output of SUB 'bound_level_interp'
                   !                                  and input of SUB 'solve_turb_budgets')

! Variables used for 3D-shear calculations:
! -----------------------------------------------

REAL (KIND=vp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
                   ! half-level values of:
  hdef2,         & ! horizontal deformation square at half levels  ( 1/s2 )
  hdiv,          & ! horizontal divergence                   ,,    ( 1/s )

  dwdx,          & ! zonal      derivative of vertical wind  ,,    ( 1/s )
  dwdy             ! meridional derivative of vertical wind  ,,    ( 1/s )

!Note: These variables are the result of horizontal gradient operations
!       and need to be precalculated dependent on the respective horizontal grid!

! Tendency fields for the prognostic variables:
! -----------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
                   ! main-level values of:
  u_tens,        & ! u-tendency                                    ( m/s2)
  v_tens,        & ! v-tendency                                    ( m/s2)
  t_tens           ! t-tendency                                    ( K/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
                   ! half-level values of:
  tketens          ! diffusion tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
                   ! half-level values of:
  tketadv          ! advection tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=wp), DIMENSION(:),   TARGET, INTENT(IN)    :: &
  trop_mask,     & ! mask-factor (1: within tropics; 0: within extra-tropics)
                   ! used for vertical smoothing of TKE forcing terms
  innertrop_mask

REAL (KIND=wp), DIMENSION(:,:),         OPTIONAL, INTENT(IN)    :: &
                   ! main-level values of:
  ut_sso,        & ! u-tendency due to the SSO-Scheme              ( 1/s )
  vt_sso           ! v-tendency due to the SSO-Scheme              ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(OUT)    :: &
                   ! half-level values of:
  edr,           & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
  tket_sso,      & ! TKE-tendency due to SSO wake production       (m2/s3)
  tket_nstc,     & ! TKE-tendency due to near-surf. therm. circul. (m2/s3)
  tket_hshr,     & ! TKE-tendency due to separ. horiz. shear       (m2/s3)
  tket_buoy,     & ! TKE tendency due to turbulent buoyancy        (m2/s3)
  tket_fshr,     & ! TKE tendency due to full mech. shear forcing  (m2/s3)
  tket_gshr,     & ! TKE tendency due to grid-scale shear forcing  (m2/s3)
  tkhm,          & ! horizontal diffusion coefficient for momentum ( m2/s )
  tkhh             ! horizontal diffusion coefficient for scalars  ( m2/s )

REAL (KIND=wp), DIMENSION(:,:),         OPTIONAL, INTENT(IN)    :: &
                   ! half-level values of:
  tket_conv        ! TKE-tendency due to convective buoyancy       (m2/s3)


!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

INTEGER ::        &
  i, k,           & !horizontaler und vertikaler Laufindex
  kem,            & !ke oder ke1
  nvor,           & !laufende Zeittstufe des bisherigen TKE-Feldes
  it_start,       & !Startindex der Iterationen
  it_durch,       & !Durchgangsindex der Iterationen
  ndiff             !number of 1-st order variables

LOGICAL ::        &
  lini              !initialization required

REAL (KIND=wp) :: &
  fr_tke              ! z1/dt_tke

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
  ediss(:,:)    !pointer for density and eddy dissipation rate

! Lokale logical Variablen:

LOGICAL ::          &
  ldotkedif,        & !berechne (teil-)implizite Vert.diff von TKE
  lcircterm,        & !TKE-source according to raw "circul.-term"  to be considered
  loutthcrc,        & !                                            to be calc. as output
  ltkemcsso,        & !TKE-source due to ordinary mech. SSO circ.  to be considered
  loutmcsso,        & !                                            to be calc. as output
  ltkeshshr,        & !TKE-source due to separ. horiz. shear-circ. to be considered
  loutshshr,        & !                                            to be calc. as output
  lssintact,        & !trenne Skalen-Interaktionsterme vom mech. Forcing ab
  ltkeadapt,        & !full TKE-adaptation to shear-related part of LLDCs
  lsrfshear,        & !specific preparations for surface-layer adaptations to additional "P"-level shear
  luse_mask,        & !use tropical mask for vertical smoothing of TKE forcing terms

  lcond                !particular condition fulfilled

! Lokale Integer-Hilfsvariablen:

INTEGER ::          &
  ii,kk, n,m,       & !Indices fuer diverse Schleifen
  ku,k1,k2,         & !Schicht-Indices
  itndcon             !Index fuer Modus der Tendenzberuecksichtigung

! Lokale real Variablen:

REAL (KIND=wp) ::   &

! Hilfsvariablen:
  wert, val1, val2,   & ! Platzhalter fuer beliebige Zwischenergebnisse
  fakt, fakt1, fakt2, & !  ,,          ,,     ,,     Faktoren

! Platzh. fuer thermodynamische Hilfsgreossen
  flw_h2o_g,        & !                 rc/(1+(lh_v/cp_d)*d_qsat/d_T)
  flw_tet_l           !epr*d_qsat/d_T*rc/(1+(lh_v/cp_d)*d_qsat/d_T)

REAL (KIND=wp) ::   &

! Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2
  fh2,fm2,          &

! Platzh. fuer horiz. Geschw.-Komponenten und bel. Geschw.:
  vel1,vel2,velo,   &

! Platzh. fuer den Kehrwert von 'grav' und 'akt':
  edgrav,           &

! Platzh. fuer verschiedene Laengenmasse:
  com_len, hk,hu,   & ! allgem. Laengenskala, Hoehe ueber Grund  und untere Hoehenbegrenzung
  lh,lm,   & ! allgem. und stab.abh. turb. Laengenskalen fuer Skalare und Impuls
  edh                 ! Kehrwert von Schichtdicken

REAL (KIND=wp) ::   &

! Zwischenspeicher fuer
  phasdif,          & !Temperaturtendenz durch Phasendiffusion

! For empirical tuning of scale-interaction terms and minimal diffusion coefficient:
  x4, x4i             !

! Local arrays:

INTEGER ::          &
  ivtp(nmvar)         ! index of variable type

! Note:
! The following buffers wouldn't be necessary, if the related pointers above
! were allowed to be allocated at run time:

REAL (KIND=wp), DIMENSION(:,:), POINTER, CONTIGUOUS :: &
  prss, & ! near-surface pressure (Pa)
  tmps, & ! near-surface temperature-varible (K)
  vaps, & ! near-surface humidity-variable
  liqs, & ! near-surface liquid water content

  cur_prof, upd_prof, sav_prof, &
  expl_mom, impl_mom, invs_mom, &
  eff_flux

TYPE (varprf) :: pvar(0:naux+2) !vector of vertical variable profiles at main- and boundary levels

! these fields are still taken as local arrays, because the CRAY compiler cannot do the
! same optimizations with OMP threadprivate variables

REAL (KIND=wp), TARGET ::  &
  ! targets of used pointers
  diss_tar   (nvec,ke1)      ! target for eddy dissipation rate (m2/s3)

REAL (KIND=wp), TARGET ::  &
  ! internal atmospheric variables
  len_scale(nvec,ke1),     & ! turbulent length-scale (m)
  hor_scale(nvec,ke),      & ! effective hoprizontal length-scale used for sep. horiz. shear calc. (m)

  l_scal   (nvec),         & ! reduced maximal turbulent length scale due to horizontal grid spacing (m)

  fc_min   (nvec),         & ! minimal value for TKE-forcing (1/s2)

  shv      (nvec,ke1),     & ! shelf of any variable (related to additional shear of NTCs)
  frh      (nvec,ke1),     & ! thermal forcing (1/s2) or thermal acceleration (m/s2)
  frm      (nvec,ke1),     & ! mechan. forcing (1/s2) or mechan. accelaration (m/s2)
  ftm      (nvec,ke1),     & ! mechan. forcing (1/s2) by traditional (pure mean) shear

  dicke    (nvec,ke1),     & ! any (effective) depth of model layers (m) or other auxilary variables
  hlp      (nvec,ke1),     & ! any 'help' variable

  zaux     (nvec,ke1,ndim),& ! auxilary array containing thermodynamical properties on boundary levels:
                             ! (1:ex_fakt, 2:cp_fakt, 3:dQs/dT, 4:g_tet l, 5:g_vap)
                             ! or various auxilary variables for calculation of implicit vertical diffusion

  can      (nvec,kcm:ke1), & ! auxilary array valid for the vertically resolved canopy
  layr     (nvec),         & ! any variable at a specific layer
  lays     (nvec,2)          ! any (2-D) vector of variables at a specific layer

REAL (KIND=wp)         ::  &
  grad     (nvec,nmvar),   &  ! any vertical gradient
  hig      (nvec,2),       &  ! obere und untere Referenzhoehe bei der Bildung
                              ! nicht-lokaler Gradienten
  xri      (nvec,ke)          ! for empirical tuning of scale-interaction terms and minimal diffusion coefficient

! Hoehenvieaus, wie etwa bei der Berechnung der gemittelten Profile bei nicht-lokaler Gradient-Berechnung:
INTEGER                ::  &
  levs     (nvec,2)

LOGICAL, PARAMETER :: ldebug=.FALSE.

INTEGER :: my_cart_id, my_thrd_id

LOGICAL :: lzacc

!---- End of header ------------------------------------------------------------

!===============================================================================

! All variables and their tendencies are defined at horizontal mass positions.
! This routine does not contain any horizontal operations!

 lzacc = (iini <= 0) !not for ordinary initialization

 !Note: 'lzacc' is equal to ".NOT.lini", as defined by calling SUB 'turb_setup' further below.

 ndiff=nmvar      !number of 1-st order variables used in the turbulence model
                  !without additional tracer: these are treated in vertdiff

! According to module 'mo_turbdiff_config' it holds:
!      nmvar = nscal+nvel: number of model variables being dynamically active for turbulence
!              nvel  = 2    active horizontal wind components:  'u_m', 'v_m'
!                           u_m = 1:     zonal      wind
!                           v_m = 2:     meridional wind
!              nscal = 3    active scalar variables 1st order:  'tem', 'vap', 'liq'
!                           tem   = 3:   temperature
!                           vap   = 4:   water vapor mixing ration
!                           liq   = 5:   liquid water
!
!                 but also: tem_l = 3:   liquid water temperature
!                           tet   = 3:   potential temperature
!                           tet_l = 3:   moist (liquid water?) potential temperature
!                           h2o_g = 4:   total water content

 kem=ke !lowest model-layer, SUB 'turbdiff' is applied to

  ! Pointer assignments:

  IF (PRESENT(edr)) THEN
     ediss => edr
  ELSE
     ediss => diss_tar
  END IF

  prss => zvari(:,:,0)   ! half-level pressure (Pa)
  tmps => zvari(:,:,tet) ! half-level temperature-variable (K)
  vaps => zvari(:,:,vap) ! half-level humidity-variable
  liqs => zvari(:,:,liq) ! half-level liquid water content

!     Fuer die Turb.par. benutzter Variablensatz auf Hauptflaechen:
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)

!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     und hat einen der Werte u_m,v_m,tet_l,h2o_g,liq;
!                        bzw. tem,vap,liq.
!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)
!     Am Ende des Programms wird zvari mit den (Co-Varianzen) der
!     Geschwindigkeitskomponenten ueberschrieben, die fuer die
!     Berechung der Horizontaldiffusion benoetigt werden.
!     zvari() enthaelt spaeter auch die nmvar (nichtlokalen) vertikalen
!     Gradienten und auch die durch Wirkung der subskaligen Kondensation
!     veraenderten (effektiven) vertikalen Gradienten.
!     Zum Schluss enthaelt zvari() fuer die turbulente Horizontaldiff.
!     benoetigte Komponenten des turbulenten Spannungstensors.


  !Begin of GPU data region
  !Input
  !$ACC DATA &
  !Working arrays
  !$ACC   CREATE(diss_tar, ivtp, hig) &
  !$ACC   CREATE(len_scale, hor_scale, l_scal, fc_min) &
  !$ACC   CREATE(shv, frh, frm, ftm, dicke, hlp) &
  !$ACC   CREATE(zaux, can, layr, lays, grad, hig, xri, levs) &
  !$ACC   PRESENT(c_big, c_sml, r_air) &
  !$ACC   IF(lzacc)


  !Note:
  !If a tendency field of an ordinary prognostic variable is not present,
  !the related time step increment due to turbulent diffusion will be
  !added to the prognostic variable directly.

  edgrav=z1/grav

  DO n=1,ndiff
    IF (n <= nvel) THEN
      ivtp(n)=mom
    ELSE
      ivtp(n)=sca
    END IF
  END DO
  !$ACC UPDATE DEVICE(ivtp) ASYNC(1) IF(lzacc)

  IF (l3dturb .AND..NOT. (PRESENT(tkhm) .AND. PRESENT(tkhh))) THEN
    CALL finish("", 'ERROR *** 3D-diffusion with not present horiz. diff.coeffs. ***')
  END IF

!-------------------------------------------------------------------------------
  CALL turb_setup (tdc=tdc, i_st=ivstart, i_en=ivend, k_st=1, k_en=ke1, &
                   iini=iini, dt_tke=dt_tke, nprv=nprv, l_hori=l_hori, &
                   ps=ps, t_g=t_g, qv_s=qv_s, qc_a=qc(:,ke), &
                   lini=lini, it_start=it_start, nvor=nvor, fr_tke=fr_tke, &
                   l_scal=l_scal, fc_min=fc_min, &
                   prss=prss(:,ke1), tmps=tmps(:,ke1), vaps=vaps(:,ke1), liqs=liqs(:,ke1), rcld=rcld, &
                   lacc=lzacc)

  !$ACC WAIT(1)
!-------------------------------------------------------------------------------


my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
my_thrd_id = omp_get_thread_num()
#endif

!------------------------------------------------------------------------------------
! 0)  Berechnung der Erhaltungsvariablen (auf 'zvari') samt des Bedeckungsgrades
!     und thermodynamischer Hilfgroessen, sowie der turbulenten Laengenskalen:
!------------------------------------------------------------------------------------

  !Additional TKE-sources due to following processes can be calculated:
  lcircterm=(tdc%pat_len > z0) !according to raw "circulation term"

  ltkeshshr=(tdc%a_hshr > z0)  !separ. horiz. shear-circ.

                           !ordinary mech. SSO-circ.:
  ltkemcsso=(lrunsso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso))
            !SSO-scheme is running and related tendencies are present

  !Output required for addit. TKE-sources (STIC-terms) due to shear by NTCs related to:
  loutthcrc=(lcircterm .AND. tdc%loutnst .AND. PRESENT(tket_nstc)) !raw "circulation-term"
  loutshshr=(ltkeshshr .AND. tdc%loutshs .AND. PRESENT(tket_hshr)) !separated horiz. shear production
  loutmcsso=(ltkemcsso .AND. tdc%loutsso .AND. PRESENT(tket_sso))  !mechanical SSO production

  !Consideration of additional TKE-sources due to:
  lcircterm=(lcircterm .AND. tdc%ltkenst)            !raw "circulation-term"
  ltkeshshr=(ltkeshshr .AND. tdc%ltkeshs)            !separ. horiz. shear-circ.
  ltkemcsso=(ltkemcsso .AND. tdc%ltkesso)            !ordinary mech. SSO-circ.

  ldotkedif=(tdc%c_diff > z0)

  lssintact=((ltkemcsso.OR.ltkeshshr.OR.tdc%ltkecon) .AND. tdc%imode_adshear == 1)

  ltkeadapt=(tdc%imode_tkemini == 2) !full TKE-adaptation to shear-related part of LLDCs,
  lsrfshear=(tdc%rsur_sher > 0 .OR. (tdc%imode_trancnf < 4 .AND. tdc%imode_suradap >= 1))
            !specific preparations for surface-layer adaptations to additional "P"-level shear

  ! Thermodynamische Hilfsvariablen auf Hauptflaechen:
  CALL adjust_satur_equil( tdc=tdc, i1dim=nvec, khi=1, ktp=1,       & !in

           i_st=ivstart, i_en=ivend, k_st=1, k_en=ke,               & !in

           lcalrho=.FALSE., lcalepr=.FALSE.,                        & !in
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        & !in

           icldmod=tdc%icldm_turb,                                  & !in

           zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          & !in

           prs=prs, t=t,     qv=qv,    qc=qc,                       & !in

           psf=ps,                                                  & !in

           rcld=rcld,  & !inp: std. deviat. of local super-saturat.
                         !out: saturation fraction (cloud-cover)

           dens=rhoh,         exner=epr,                            & !out
           r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3),                    & !out
           g_tet=zaux(:,:,4), g_h2o=zaux(:,:,5),                    & !out

           tet_liq=zvari(:,:,tet_l), q_h2o=zvari(:,:,h2o_g),        & !out
                                     q_liq=zvari(:,:,liq),          & !out

           lacc=lzacc )

  ! Thermodynamische Hilfsvariablen auf Unterrand der Prandtl-Schicht:
  CALL adjust_satur_equil( tdc=tdc, i1dim=nvec, khi=1, ktp=1,       & !in

           i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1,            & !in

           lcalrho=.TRUE., lcalepr=.TRUE.,                          & !in
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        & !in

           icldmod=tdc%icldm_turb,                                  & !in

           zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          & !in

           !Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
           !         und an COSMO-Version angepasste Interpretation von "icldmod=-1":

           prs=prss, t=tmps, qv=vaps, qc=liqs,                      & !in (surface values at level 'ke1')

           psf=ps, fip=tfh,                                         & !in

           rcld=rcld,  & !inp: std. deviat. of local super-saturat.
                         !out: saturation fraction (cloud-cover)

           dens=rhon,         exner=zaux(:,:,1),                    & !out
           r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3),                    & !out
           g_tet=zaux(:,:,4), g_h2o=zaux(:,:,5),                    & !out

           tet_liq=zvari(:,:,tet_l), q_h2o=zvari(:,:,h2o_g),        & !inout (inp as target of 'tmps, vaps')
                                     q_liq=zvari(:,:,liq),          & !out

           lacc=lzacc )

  ! Beachte:
  !     'zvari(:,ke1,tet_l)' und 'zvari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen
  !      am Unterrand der Prandtl-Schicht (zero-level). Die Werte im Niveau 'ke' wurden dabei
  !      zur Interpolation benutzt.
  !     'zaux(:,ke1,1)' enthaelt den Exner-Faktor im zero-level.
  !     Das Feld 'zaux(:,:,1) wird im Folgenden mit dem Exner-Faktor auf NF belegt.

  ! Kommentar:
  !     Der 2-te Aufruf von 'adjust_satur_equil' stellt die unteren Randwerte der thermodyn. Variablen
  !      zur Verfuegung. Dies koennte in den 1-ten Aufruf integriert werden, wenn alle thermodyn.
  !      Modell-Variablen bis "k=ke1" allociert waeren. Dies wuerde Vieles vereinfachen!

  ! Berechnung der horizontalen Windgeschwindigkeiten und Schichtdicken:

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=1,ke
!DIR$ IVDEP
    DO i=ivstart, ivend
      dicke(i,k)=hhl(i,k)-hhl(i,k+1)  ! Berechnung der Schichtdicken und der Dichte auf Nebenflaechen
      zvari(i,k,u_m)=u(i,k)
      zvari(i,k,v_m)=v(i,k)
    END DO
  END DO
  !$ACC END PARALLEL

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
     zvari(i,ke1,u_m)=zvari(i,ke,u_m)*(z1-tfm(i))
     zvari(i,ke1,v_m)=zvari(i,ke,v_m)*(z1-tfm(i))
  END DO
  !$ACC END PARALLEL

  ! Noetige Interpolationen auf Nebenflaechen:

  n=0; pvar(n)%bl => rcld       ; pvar(n)%ml => rcld        !cl_covr
  n=1; pvar(n)%bl => zaux(:,:,n); pvar(n)%ml => epr         !ex_fakt
  m=MERGE( n, n+1, tdc%lcpfluc ) !interpolation of "r_cpd=zaux(:,:,2)" only at "lcpfluc=T"
  DO WHILE (m < naux)                                       !2:cp_fakt, 3:dQs/dT, 4:g_tet l, 5:g_vap
    n=n+1; m=m+1; pvar(n)%bl => zaux(:,:,m) ; pvar(n)%ml => zaux(:,:,m)
  END DO
  IF (lcircterm .OR. loutthcrc) THEN !Der bisherige "Zirkulationsterm" muss berechnet werden
    n=n+1; pvar(n)%bl => prss; pvar(n)%ml => prs            !air_pres
  END IF
  !___________________________________________________________________________
  !test: different interpolation weights for bl-interpolation of rho:<
  ! Volume-weighted interpolation:
  !CALL bound_level_interp( ivstart, ivend, 2, ke, &
  !                         nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=dicke, lacc=lzacc)
  ! Mass-weighted interpolation (included into mulit-variable interpolation:)
  n=n+1; pvar(n)%bl => rhon; pvar(n)%ml => rhoh             !air_dens
 !___________________________________________________________________________
  CALL bound_level_interp( ivstart, ivend, 2, ke, &
                           nvars=n+1, pvar=pvar, depth=dp0, auxil=hlp, lacc=lzacc)

  ! Berechnung der turbulenten Laengenscalen:

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
    len_scale(i,ke1)=gz0(i)*edgrav
  END DO
  !$ACC END PARALLEL

  IF (PRESENT(c_big) .AND. PRESENT(r_air)) THEN

    !US: Up to now it is kcm = ke+1 and the next vertical loop will not be executed!!
    !    If a canopy layer is implemented, kcm will be <= ke.

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=ke,kcm,-1 !Innerhalb des Bestandesmodells
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        IF (c_big(i,k) > z0) THEN
          ! Die turbulente Laengenskala wird durch die Laengenskala
          ! der lufterfuellten Zwischenraeume limitiert:
          len_scale(i,k)=MIN( dicke(i,k)+len_scale(i,k+1), z1/(c_big(i,k)*SQRT(z1/EXP(r_air(i,k))-z1)) )
        ELSE
          len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
  ENDIF

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP SEQ
  DO k=kcm-1,1,-1
!DIR$ IVDEP
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
    END DO
  END DO
  !$ACC END PARALLEL

  ! Uebergang von der maximalen turbulenten Laengenskala zur
  ! effektiven turbulenten Laengenskala:

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=ke1,1,-1
!DIR$ IVDEP
    DO i=ivstart, ivend
      len_scale(i,k)=tdc%akt*MAX( tdc%len_min, l_scal(i)*len_scale(i,k)/(l_scal(i)+len_scale(i,k)) )
    END DO
  END DO
  !$ACC END PARALLEL

  ! Initialisierung der Felder fuer tke,tkvh,tkvm:

!------------------------------------------------------------------------------------------------
  IF (lini) THEN  !nur beim allerersten Durchgang
!------------------------------------------------------------------------------------------------

    ! Erste Schaetzwerte aus vereinfachtem TKE-Gleichgewicht:

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=2,kem
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR PRIVATE(com_len, edh, fh2, fm2, fakt, lm, lh, val1, val2)
      DO i=ivstart, ivend

        ! Der Einfachheit halber erfolgt nur eine lokale Berechnung der vertikalen Gradienten:

        com_len=len_scale(i,k)
        edh=z2/(hhl(i,k+1)-hhl(i,k-1))

        grad(i,u_m  )=(zvari(i,k,u_m  )-zvari(i,k-1,u_m  ))*edh
        grad(i,v_m  )=(zvari(i,k,v_m  )-zvari(i,k-1,v_m  ))*edh
        grad(i,tet_l)=(zvari(i,k,tet_l)-zvari(i,k-1,tet_l))*edh
        grad(i,h2o_g)=(zvari(i,k,h2o_g)-zvari(i,k-1,h2o_g))*edh

        fh2=zaux(i,k,4)*grad(i,tet_l)+zaux(i,k,5)*grad(i,h2o_g)
        fm2=MAX( grad(i,u_m)**2+grad(i,v_m)**2, fc_min(i) )

        ! Vereinfachte Loesung mit "Rf=Ri":
        IF (fh2 >= (z1-rim)*fm2) THEN !die krit. Ri-Zahl wird ueberschritten
          !'lm' sowie 'lh' werden durch 'lm' bei der krit. Ri-Zahl angenaehert:
          fakt=z1/rim-z1
          lm=com_len*(sm_0-(a_6+a_3)*fakt)
          lh=lm
        ELSE
          fakt=fh2/(fm2-fh2)
          lm=com_len*(sm_0-(a_6+a_3)*fakt)
          lh=com_len*(sh_0-a_5*fakt)
        END IF

        IF (ltkeinp) THEN
          tke(i,k,1)=tke(i,k,ntur)
        ELSE
          val1=lm*fm2; val2=lh*fh2
          wert=MAX( val1-val2, rim*val1 )
          tke(i,k,1)=MAX( tdc%vel_min, SQRT(tdc%d_mom*com_len*wert) ) !Initialwert fuer SQRT(2TKE)
        END IF

        val1=con_m; tkvm(i,k)=lm*tke(i,k,1)
        val2=con_h; tkvh(i,k)=lh*tke(i,k,1)

        !Note:
        !'tk[h|m]min' are, fist of all, foreseen as lower limits for 'vertdiff'-calculations; hence,
        ! they are not required for initialization.
        !Nevertheless, since the 'tkv[m|h]' from the previous time-step are required as input of the
        ! Turbulence Model (TMod) in SUB 'solve_turb_budgets' (dependent on 'imode_stbcalc'), at least
        ! the laminar limit is used for securing a reasonable start of this kind of time-step iteration.

        IF (tdc%imode_tkemini >= 2) THEN !any adaptation of TKE and the TMod. to lower limits
          tke(i,k,1)=tke(i,k,1)*MAX( z1, val2/tkvh(i,k) ) !adapted 'tke'
        ENDIF
        IF (ltkeadapt) THEN !full adaptation of TKE and the TMod. to lower limits
          tprn(i,k)=tkvm(i,k)/tkvh(i,k) !turbulent Prandtl-number as calcuated by the simplified TMod.
                                        ! used for initialization
          !Note:
          !At this simplified turbulence diagnostics (being applied for initialization only) additional
          ! shear-forcing by NTCs is not considered.
          !Accordingly, the indirect shear-impact related to 'tk[h|m]min' is missing as well, not at least,
          ! because it's application should be connected with further modulation (e.g. dependent on Ri-number
          ! or on the distance from the surface).
          !Moreover, due to time-step smoothing of 'tke' (through 'tkesmot'), any not realistic and large deviation
          ! of 'tke' from the quilibrium-solution of the TMod, might have a quite long-standing detrimental impact.
          !Thus, for initialization, only the laminar limit (including the related TKE-adaptation) is applied,
          ! which secures a reasonable start of time-step interation.
          !See also notes related to 'ltkeadapt' further below!
        END IF

        tkvm(i,k)=MAX( val1, tkvm(i,k) ) !'tkvm' with lower limit
        tkvh(i,k)=MAX( val2, tkvh(i,k) ) !'tkvh' with lower limit

        ! Am Anfang konnte noch keine Diffusion von q=SQRT(2*TKE) berechnet werden:
        tketens(i,k)=z0

      END DO
    END DO
    !$ACC END PARALLEL

!DIR$ IVDEP
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      tke(i,1,1)=tke(i,2,1)
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO n=2,ntim
      DO k=1,kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          tke(i,k,n)=tke(i,k,1)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT

  END IF ! (lini)

  IF (ltkeadapt) THEN !full adaptation of TKE and the TMod. to lower limits
    !Note:
    !'tkv[m|h]' may have been artificially increased by the application of lower limits.
    !However, after this manipulation, only 'tkvh' can be assumed to represent a turbulent diff. coeff.,
    ! while 'tkvm' has possibly been further increased in order to compensate missing SGS drag in the
    ! momentum budgets. In contrast to real momentum fluxes, the virtual downward momentum flux, associated
    ! to this compensation of missing drag-forces, is not directly driven by additional shear.
    !As the amplification of 'tkvh' should purely be due to not yet considered shear by NTCs, the shear-induced
    ! part of the 'tkvm'-amplific. is equal to the full 'tkvh'-amplific. (provided that turbulent Prandtl-number
    ! would not be altered, if the missing shear was applied).

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=1,kem
!DIR$ IVDEP
      DO i=ivstart, ivend
        tkvm(i,k)=tprn(i,k)*tkvh(i,k) !current, true turbulent, diff. coeff. for momentum
      END DO
    END DO
    !$ACC END PARALLEL

    !Note:
    !Now, 'tkvm' is a true turbulent diff. coeff. relieved from artif. drag contrib. again,
    ! and both "tkvm/tke" as well as "tkvh/tke" can be used as an estimate for the turbulent
    ! stability length, required as input for the TMod. in SUB 'solve_turb_budgets'.
    !'tprn' contains the turbulent Prandtl-number as an output of the before applied TMod.
    ! (as it was called at the previous time step).
  END IF !TKE-source due to new version of thermal SSO-circ. required

!------------------------------------------------------------------------------------
! 1)  Berechnung der benoetigten vertikalen Gradienten und Abspeichern auf 'zvari':
!------------------------------------------------------------------------------------

  ! Am unteren Modellrand:

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
    lays(i,mom)=tvm(i)/(tkvm(i,ke1)*tfm(i))
    lays(i,sca)=tvh(i)/(tkvh(i,ke1)*tfh(i))
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO n=1,nmvar
!DIR$ IVDEP
    DO i=ivstart, ivend
      zvari(i,ke1,n)=(zvari(i,ke,n)-zvari(i,ke1,n))*lays(i,ivtp(n))
    END DO
  END DO
  !$ACC END PARALLEL

  ! An den darueberliegenden Nebenflaechen:

  ! Berechnung lokaler Gradienten:

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(com_len)
  DO k=ke,2,-1
!DIR$ IVDEP
    DO i=ivstart, ivend
      com_len=(hhl(i,k-1)-hhl(i,k+1))*z1d2
      hlp(i,k)=z1/com_len
      dicke(i,k)=rhon(i,k)*com_len*fr_tke
    END DO
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP SEQ
  DO n=1,nmvar
#ifdef __INTEL_COMPILER
    FORALL(k=2:ke,i=ivstart:ivend)                        &
             zvari(i,k,n)=(zvari(i,k-1,n)-zvari(i,k,n))*hlp(i,k)
#else
    !$ACC LOOP SEQ
    DO k=ke,2,-1
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        zvari(i,k,n)=(zvari(i,k-1,n)-zvari(i,k,n))*hlp(i,k)
      END DO
    END DO
#endif
  END DO
  !$ACC END PARALLEL

!------------------------------------------------------------------------------------
! 2)  Berechnung der verallgemeinerten Antriebsfunktionen einschliesslich der
!     Korrekturen innerhalb der Rauhigkeitsschicht (samt der Windtendenz durch Formreibung)
!     und der Scherung durch nicht-turbulente subskalige Stroemungen:
!------------------------------------------------------------------------------------

  ! Thermal forcing:
  !-------------------------

  ! Achtung:
  !'frh'(ke1) wird fuer Zirkulationsterm und Temperaturkorrektur benoetigt

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=2,ke1
!DIR$ IVDEP
    DO i=ivstart, ivend
      frh(i,k)=zaux(i,k,4)*zvari(i,k,tet_l) + zaux(i,k,5)*zvari(i,k,h2o_g)
    END DO
  END DO
  !$ACC END PARALLEL

  !Note: 'zaux(:,:,5)' is free now.

  ! Total mechanical forcing:
  !--------------------------

  !hdef2 = (d1v2+d2v1)^2 + (d1v1-d2v2)^2 !horizontal deformation square     (at half levels)
  !hdiv  = (d1v1+d2v2)                   !horizontal wind-divergence            ,,

  !dwdx !zonal      derivation of vertical wind                                 ,,
  !dwdy !meridional derivation of vertical wind                                 ,,

  !vel_div=hdiv+dzdw=0 !Incomressibility

  !itype_sher = 0 : only single column vertical shear
  !             1 : previous and additional 3D horiz. shear correction
  !             2 : previous and additional 3D vertc. shear correction

  !ltkeshshr: consider separated non-turbulent horizontal shear mode for TKE forcing

  ! Mechanical forcing by vertical shear:

  IF (tdc%itype_sher == 2 .AND. (PRESENT(dwdx) .AND. PRESENT(dwdy) .AND. PRESENT(hdiv))) THEN

    !Include 3D-shear correction by the vertical wind (employing incomressibility):

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=2,kem
!DIR$ IVDEP
      DO i=ivstart, ivend
        frm(i,k)=MAX( (zvari(i,k,u_m)+dwdx(i,k))**2 + (zvari(i,k,v_m)+dwdy(i,k))**2 &
                    + z3*hdiv(i,k)**2, fc_min(i) )
      END DO
    END DO
    !$ACC END PARALLEL

  ELSE

    !Load pure single column shear:

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=2,kem
!DIR$ IVDEP
      DO i=ivstart, ivend
        frm(i,k)=MAX( zvari(i,k,u_m)**2+zvari(i,k,v_m)**2, fc_min(i))
      END DO
    END DO
    !$ACC END PARALLEL

  END IF

  ! Mechanical forcing by horizontal shear:

  IF (PRESENT(hdef2)) THEN
    IF (tdc%itype_sher >= 1) THEN   !Apply horizontal 3D-shear correction:
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          frm(i,k)=frm(i,k)+hdef2(i,k) !extended shear
        END DO
      END DO
      !$ACC END PARALLEL
    END IF
  END IF

  IF (lssintact .OR. & !full shear (including NTC impact) and mean shear to be separated
      tdc%loutbms ) THEN !output of (traditional) purely grid-scale shear requested
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=2,kem
!DIR$ IVDEP
      DO i=ivstart, ivend
        ftm(i,k)=frm(i,k) !save traditional (pure mean) shear
      END DO
    END DO
    !$ACC END PARALLEL
  ELSEIF (tdc%rsur_sher > z0) THEN !addit. shear by NTCs|LLDCs at "P"-level to be transmitted to "0"-level
!DIR$ IVDEP
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
       ftm(i,ke)=frm(i,ke) !save traditional (pure mean) shear at lowest half level
    END DO
    !$ACC END PARALLEL
  END IF

  ! Preparation for Richardson-number-dependent factor used for correcting
  !  the additional TKE-production by horizontal shear and SSO as well as
  !  the minimum diffusion coefficients.
!GZ: For tuning. <
!---------------------------------------------------------------------------
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=2,ke
!DIR$ IVDEP
    DO i=ivstart, ivend
      xri(i,k)=EXP( z2d3*LOG( MAX( 1.e-6_wp, frm(i,k) ) / & !1/Ri**(2/3)
                              MAX( 1.e-5_wp, frh(i,k) ) ) )
    END DO
  END DO
  !$ACC END PARALLEL
!---------------------------------------------------------------------------
!GZ>

 !Additional impact by separated horizontal shear:

  IF (PRESENT(hdef2)) THEN

!   IF ((tdc%ltkeshs .OR. (tdc%loutshs .AND. PRESENT(tket_hshr))) .AND. PRESENT(hdiv)) THEN
    IF ((ltkeshshr .OR. loutshshr) .AND. PRESENT(hdiv)) THEN
      !Include separated horizontal shear mode:

      fakt=z1/(z2*sm_0)**2; wert=tdc%a_hshr*tdc%akt*z1d2

!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        layr(i)=wert*l_hori(i) !uncorrected effective horizontal length scale
      END DO
      !$ACC END PARALLEL

      IF (tdc%imode_shshear == 2) THEN !Ri-dependent length sclale correction
!GZ: For tuning. <
!---------------------------------------------------------------------------
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(x4i, x4)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            ! Factor for variable 3D horizontal-vertical length scale proportional to 1/SQRT(Ri),
            ! decreasing to zero in the lowest two kilometer above ground
            ! from ICON 180206
            x4i = MIN( 1._wp, 0.5e-3_wp*(hhl(i,k)-hhl(i,ke1)) )
            x4 = (3._wp - 2._wp*x4i)*x4i**2  !low-level reduction factor
            hor_scale(i,k) = layr(i)*MIN( 5.0_wp, MAX( 0.01_wp, x4*xri(i,k) ) ) &
                                    /MAX(1._wp,0.2_wp*tke(i,k,nvor))
          END DO
        END DO
        !$ACC END PARALLEL
!---------------------------------------------------------------------------
!GZ>
      ELSE

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            hor_scale(i,k) = layr(i)
          END DO
        END DO
        !$ACC END PARALLEL

      ENDIF

      !strain velocity of the separated horizontal shear mode saved in 'hlp':
      IF (tdc%imode_shshear == 0) THEN !former variant based on 3D-shear and incompressibility
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            hlp(i,k)=hor_scale(i,k)*SQRT(hdef2(i,k)+hdiv(i,k)**2) !not equal to trace of 2D-strain tensor
          END DO
        END DO
        !$ACC END PARALLEL
      ELSE !new variant in accordance with the trace constraint for the separated horizontal strain tensor
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(wert)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            wert=fakt*hdiv(i,k)
            hlp(i,k)=hor_scale(i,k)*(SQRT(wert**2+hdef2(i,k))-wert) !equal to trace of 2D-strain tensor
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      IF (ltkeshshr .AND. l3dturb) THEN
        ! Load related horizontal diffusion coefficients:
        fakt=sh_0/sm_0
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            tkhm(i,k)=hor_scale(i,k)*hlp(i,k) !for momentum related to the sep. shear mode
            tkhh(i,k)=fakt*tkhm(i,k)          !for scalars    ,,       ,,            ,,
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          hlp(i,k)=(hlp(i,k))**3/hor_scale(i,k) !additional TKE-source by related shear
        END DO
      END DO
      !$ACC END PARALLEL

      IF (loutshshr) THEN
        !Load output variable for the TKE-source by separated horiz. shear:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            tket_hshr(i,k)=hlp(i,k)
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      IF (ltkeshshr) THEN
        !Consider separated horizontal shear mode in mechanical forcing:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            frm(i,k)=frm(i,k)+hlp(i,k)/tkvm(i,k) !extended shear
          END DO
        END DO
        !$ACC END PARALLEL


      END IF    ! ltkeshshr

    END IF      ! (ltkeshshr .OR. loutshshr) .AND. PRESENT(hdiv)

  END IF !(PRESENT (hdef2))

  !Erweiterung der Vertikalscherung mit verallgemeinerten Scher-Beitraegen durch die nicht-turbulente
  ! subskalige Stroemung (STIC-Terme):

!------------------------------------------------------------------------------------------------
  IF (.NOT.lini) THEN !not for initialization since 'ut_sso' and 'vt_sso' may be calculated later
!------------------------------------------------------------------------------------------------

    IF (ltkemcsso .OR. loutmcsso) THEN !clac. or output of addit. TKE-sources due to ordinary mech. SSO-circ.
      !SSO-Schema ist aktiv, die SSO-Tendenzen des Windes sind vorhanden
      ! und die Berechnung der TKE-Quellen durch mech. SSO-Wirbel wird gewuenscht:

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=1,kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          hlp(i,k)=ut_sso(i,k)*u(i,k)+vt_sso(i,k)*v(i,k)
        END DO

          !Note:
          !Horizontal wind components and SSO-tendencies refer to horizontal mass centeres here.
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(wert)
      DO k=2,kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          wert=MAX( z0, -zbnd_val(hlp(i,k), hlp(i,k-1), dp0(i,k), dp0(i,k-1)) )

          !Note:
          !Although the SSO-tendencies 'ut_sso' and 'vt_sso' should never become positive,
          ! the 'MAX'-function is used for security.

          IF (loutmcsso) THEN !output of addit. TKE-sources due to ordinary mech. SSO-circ.
            tket_sso(i,k)=wert
          END IF

          IF (ltkemcsso) THEN !consideration of addit. TKE-sources due to ordinary mech. SSO-circulation
            IF (tdc%imode_tkesso == 1) THEN !without further manipulation
              frm(i,k)=frm(i,k) + wert/tkvm(i,k)
            ELSE IF (tdc%imode_tkesso == 2) THEN ! reduction in the presence of large Richardson numbers
              frm(i,k)=frm(i,k) + wert/tkvm(i,k)*MIN(1.0_wp,MAX(0.01_wp,xri(i,k)))
            ELSE IF (tdc%imode_tkesso == 3) THEN ! Reduce TKE production in the presence of large Richardson numbers
              frm(i,k)=frm(i,k) + wert/tkvm(i,k)*MIN(1.0_wp,MAX(0.01_wp,xri(i,k)))*MIN(1.0_wp,l_hori(i)/2000._wp)
            END IF
          END IF

        END DO
      END DO
      !$ACC END PARALLEL

    END IF !clac. or output of addit. TKE-sources due to ordinary mech. SSO-circ.

    !Adding shear due to sub-grid convective circulation:
    IF (lruncnv .AND. tdc%ltkecon .AND. PRESENT(tket_conv)) THEN
      !Convection scheme is active, it is desired for impacting turbulence and 'tket_conv' is present:

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2, kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          frm(i,k) = frm(i,k) + MAX( z0, tket_conv(i,k)/tkvm(i,k) )
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

!------------------------------------------------------------------------------------------------
  END IF  ! IF (.NOT.lini)
!------------------------------------------------------------------------------------------------

  IF (PRESENT(c_big) .AND. PRESENT(c_sml) .AND. kcm <= kem .AND. iini.NE.1) THEN
    ! Berechnung von Korrekturtermen innerhalb der Rauhigkeitsschicht
    ! (ausser Volumenterme, die zur Diffusion gehoeren):

    !US: at the moment kcm = kem+1, so this block is never executed.

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=kcm-1, kem !alle Hauptfl. innerh. d. R-Schicht (einschl. der Flaeche darueber)
      !$ACC LOOP GANG VECTOR PRIVATE(velo)
!DIR$ IVDEP
      DO i=ivstart, ivend
         velo=z1d2*(w(i,k)+w(i,k+1)) !Vertikalwind auf Hauptflaeche
         hlp(i,k)=SQRT(u(i,k)**2+v(i,k)**2+velo**2) !Windbetrag auf Hauptflaeche
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=kcm, kem !von oben nach unten durch Rauhigkeitsschicht

!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR PRIVATE(velo, wert)
      DO i=ivstart, ivend
        ! Formreibungskoeff. auf Hauptflaechen:
        wert=z1d2*(c_big(i,k)+c_sml(i,k)+c_big(i,k+1)+c_sml(i,k+1))

        ! Formreibungsfrequenz auf Hauptflaeche:
        wert=wert*hlp(i,k)

        ! Aufaddieren der Windtendenzen durch Fromreibung:
        u_tens(i,k)=u_tens(i,k)-wert*u(i,k)/(z1+dt_var*wert)
        v_tens(i,k)=v_tens(i,k)-wert*v(i,k)/(z1+dt_var*wert)

        ! Windbetrag auf Nebenflaechen:
!       can(i,k)=(hlp(i,k)*dp0(i,k-1)+hlp(i,k-1)*dp0(i,k))/(dp0(i,k-1)+dp0(i,k))
        can(i,k)=zbnd_val(hlp(i,k), hlp(i,k-1), dp0(i,k), dp0(i,k-1))

        ! Windbetrag unter Wirkung der Formreibung:
        velo=can(i,k)/(z1+can(i,k)*(c_big(i,k)+c_sml(i,k))*dt_var)

        ! Addition des Scherterms durch Nachlaufwirbel an grossen Rauhigkeitselementen:
        frm(i,k)=frm(i,k)+c_big(i,k)*velo**3/tkvm(i,k)

        ! Frequenz der kleinskaligen Rauhigkeitselemente:
        can(i,k)=c_sml(i,k)*can(i,k)

        !ToDo: Ev. neue Behandlung der Rauhigkeitsschicht einfuehren!
      END DO
    END DO !k=kcm,kem !von oben nach unten durch Rauhigkeitsschicht
    !$ACC END PARALLEL

  ENDIF !Berechnung von Korrekturtermen innerhalb der Rauhigkeitsschicht

  IF (tdc%rsur_sher > 0) THEN !shear factor by NTCs at lowest half-level required for surface layer
!DIR$ IVDEP
     !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=ivstart, ivend
        tfv(i)=frm(i,ke)-ftm(i,ke) !saved shear-forcing due to NTCs (without subsequent vertical smoothing)
     END DO
     !$ACC END PARALLEL
  END IF

  ! Optional output of source terms:

  IF (tdc%loutbms) THEN !output of additional TKE-sources required
     !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
     !$ACC LOOP SEQ
     DO k=2, ke1
!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
           IF (PRESENT(tket_buoy)) tket_buoy(i,k) = frh(i,k)*tkvh(i,k)
           IF (PRESENT(tket_fshr)) tket_fshr(i,k) = frm(i,k)*tkvm(i,k)
           IF (PRESENT(tket_gshr)) tket_gshr(i,k) = ftm(i,k)*tkvm(i,k)
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

  ! Check if vertical smoothing of TKE forcing terms is needed:
  IF (tdc%frcsmot > z0) THEN
    luse_mask=(tdc%imode_frcsmot == 2 .AND. .NOT.lini)
    IF (luse_mask) THEN
      lcond=(ANY(trop_mask(ivstart:ivend) > z0)) !not all points are masked out
    ELSE
      lcond=.TRUE.
    END IF

    IF (lcond) THEN !vertical smoothing required
      ! Optionale vertikale Glaettung des mechanischen Antriebs:
      CALL vert_smooth (i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
                        disc_mom=dicke, cur_tend=frm, vertsmot=tdc%frcsmot, &
                        smotfac=trop_mask, luse_mask=luse_mask, lacc=lzacc )

      ! Optionale vertikale Glaettung des thermischen Antriebs:
      CALL vert_smooth (i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
                        disc_mom=dicke, cur_tend=frh, vertsmot=tdc%frcsmot, &
                        smotfac=trop_mask, luse_mask=luse_mask, lacc=lzacc )
    END IF
  END IF

  ! Belegung von tkvh und tkvm mit den stabilitaetsabhaengigen Laengenmassen:

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(wert)
  DO k=2,kem !stability-dependent length-scales in 'tkvh/m' only for the here included atmospheric levels
!DIR$ IVDEP
    DO i=ivstart, ivend
      wert=z1/tke(i,k,nvor)
      tkvh(i,k)=tkvh(i,k)*wert
      tkvm(i,k)=tkvm(i,k)*wert
    END DO
  END DO
  !$ACC END PARALLEL

  !Note:
  !'tke' is always the turbulent velocity scale 'tvs' in [m/s] and NOT "TKE = 1/2*q^2"!
  !'tkv[m|h]' are the turbulent stability length-scales 'ls[m|h]:=S[m|h]*len_scale' in [m] now.
  !It has been modified by the impact of artificial lower limits compared to the original output of
  ! SUB 'solve_turb_budgets' (at the previous time step), unless "ltkeadapt=T" (which is "ltkeadapt=T") is applied!
  !At "imode_tkemini=1", 'tke' has not been enlarged according to applied LLDCs, and their full impact
  ! causes respectively enlarged values of 'S[m|h]' as input to the Tmod. called by SUB 'solve_turb_budgets.
  !AT "imode_tkemini=2", 'tke' has been fully adapted to applied LLDCs, and 'tkvm'(as input to SUB 'turbdiff')
  ! has also been relieved from additional artificial drag contributions by means of 'tprn'.
  !At "imode_tkemini=3", only a reduced (tuned) 'tke-adaptation has been applied and the additional drag-impact
  ! of 'tkvm' has not been removed (by means of 'tprn').

  IF (tdc%ltmpcor .AND. tdc%lcpfluc) THEN !consideration of temperature tendency due to phase-diffusion
    !  Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=2, ke1
!DIR$ IVDEP
      DO i=ivstart, ivend
        hlp(i,k)=zaux(i,k,1)*zvari(i,k,tet_l)-tet_g  !vertical temperature gradient
      END DO
    END DO
    !$ACC END PARALLEL

    IF (tdc%icldm_turb.NE.-1) THEN !water phase changes are possible
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          hlp(i,k)=hlp(i,k)+lhocp*zvari(i,k,liq) !liquid water correction
        END DO
      END DO
      !$ACC END PARALLEL
    END IF
  END IF

!------------------------------------------------------------------------------------
! 3)  Loesung der turbulenten Bilanzgleichungen (Bestimmung von TKE und der Stabilitaetsfunktionen
!      sowie der Standard-Abweichung der lokalen Uebersaettigung)
!      und Berechnung der turbulenten Diffusionskoeffizienten:
!------------------------------------------------------------------------------------

  DO it_durch=it_start, tdc%it_end !iteration

    !Die Schleife wird nur bei der Initialisierung (d.h. beim ersten Aufruf) wiederholt,
    !um TKE-Gleichgewicht anzunaehern. Die resultierenden TKE-Werte der Zeitstufe 'ntur'
    !gehoeren in diesem Fall dann zur Zeitstufe der uebrigen prognostischen Variablen
    !('nold' bei "leap-frog" oder 'nnow' bei 2-Zeitebenen).
    !Fuer die folgenden Aufrufe wird die Schleife nur einmal durchlaufen und liefert TKE-Werte
    !die gegenueber den Vorgaengerwerten um einen Zeitschritt weiter in der Zukunft liegen,
    !also wieder zur Zeitstufe der uebrigen prognostischen Variablen gehoeren.

    CALL solve_turb_budgets( tdc=tdc, it_s=it_durch, it_start=it_start,               & !in

                             i1dim=nvec, i_st=ivstart, i_en=ivend,                    & !in
                             khi=1, ktp=1, kcm=kcm, k_st=2, k_en=kem, k_sf=ke1,       & !in
                             !Note: In case of "kem=ke", all calculation for level 'ke1' should rather
                             !       be called by SUB 'turbtran' exclusively!

                             ntur=ntur, nvor=nvor,                                    & !in

                             lssintact=lssintact,    lupfrclim=.FALSE.,               & !in
                             lpres_edr=(lsrfshear .OR. PRESENT(edr)),                 & !in
                             ltkeinp=ltkeinp,                                         & !in
                             imode_stke=tdc%imode_turb,  imode_vel_min=1,             & !in

                             dt_tke=dt_tke, fr_tke=fr_tke,                            & !in

                             fm2=frm,  fh2=frh,  ft2=ftm,                             & !in
                             lsm=tkvm, lsh=tkvh, tls=len_scale,                       & !in(out)

                             fcd=can, tvt=tketens, avt=tketadv,                       & !in
                             tke=tke, ediss=ediss,                                    & !inout, out

                             lactcnv=(tdc%icldm_turb.NE.-1),                          & !in (activ. flux-conversion)
                             laddcnv=(tdc%ltmpcor .AND. tdc%lcpfluc),                    & !in (addit. flux-conversion)

                             exner=zaux(:,:,1), r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3), & !in

                             rcld=rcld,          & !inp: effective saturation fraction (cloud-cover)
                                                   !out: std. deviat. of local super-saturat.
                                                   !(only for last iteration step)

                             lcircterm=(lcircterm.OR.loutthcrc),                      & !in
                             dens=rhon,         & !air density at boundary levels     as in
                             l_pat=l_pat, l_hori=l_hori,                              & !in

                             grd=zvari,         & !inp: vert. grads. (incl. those of tet_l, h20_g and liq)
                                                  !out: vert. grads. (incl. those of tet,   vap   and liq
                                                  !     resuting from flux conversion, provided it is executed)
                             !'prss=>zvari(:,:,0) (as inp): half-level pressure, providing (if devided by density and a
                             !                               length of coherence) a related circulation acceleration
                             !                    (as out): circulation acceleration due to thermal near surface heterogeneity,
                             !                               provided that "lcircterm=T"
                             !'zvari'-output is only calculated at the last iteration step; and it's still equal to the input,
                             ! if these calculations are not executed.

                             lacc=lzacc)                                                !in

    IF (it_durch < tdc%it_end .AND. .NOT.ltkeinp) THEN
      nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
    END IF

  END DO !iteration

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=2, kem
!DIR$ IVDEP
     DO i=ivstart, ivend
        tkvh(i,k)=MAX(tkvh(i,k)*tke(i,k,ntur),con_h)
        tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
     END DO
  END DO
  !$ACC END PARALLEL

  !Note:
  !'tke' is always the turbulent velocity scale in [m/s] and NOT "TKE = 1/2*q^2"!
  !'tkv[m|h]' are no longer the stability length-scales 'S[m|h]*len_scale' now,
  ! but the turb. diff. coeffs. in [m^2/s] again.

  IF (ltkeadapt .OR. l3dturb) THEN !correction of 'tkvm' from artificial drag impact required
     !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
     DO k=2, kem
!DIR$ IVDEP
        DO i=ivstart, ivend
           tprn(i,k)=tkvm(i,k)/tkvh(i,k) !turbulent Prandtl-number as given by the above called TMod.
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
     tke(i,1,ntur)=tke(i,2,ntur) !kein TKE-Gradient am Oberrand
  END DO
  !$ACC END PARALLEL

  IF (lsrfshear) THEN !specific preparations for surface-layer adaptations to additional "P"-level shear
!DIR$ IVDEP
     !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=ivstart, ivend
        ! Save turb. diff. coeff. for mom. at the lowest half level "P" without any lower limit:
        tfm(i)=tkvm(i,ke)
        tfh(i)=tkvh(i,ke)
     END DO
     !$ACC END PARALLEL
  END IF
!------------------------------------------------------------------------------------
!  4) Berechnung der effectiven Diffusionskoeffizienten:
!------------------------------------------------------------------------------------

  ! Beschraenkung der Diffusionskoeffizienten nach unten:

  vel1=MAX(con_m, tdc%tkmmin); vel2=MAX(con_h, tdc%tkhmin)

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(fakt, fakt1, fakt2, wert, val1, val2, x4, x4i)
  DO k=2, ke
!DIR$ IVDEP
  DO i=ivstart, ivend

    val1=vel1; val2=vel2 !default minimum diffusion values

!GZ: For tuning. <
!---------------------------------------------------------------------------
    IF (tdc%imode_tkvmini == 2 .AND. iini.NE.1) THEN !using addit. empirical modificat. by means of 'xri'

      ! Factor for variable minimum diffusion coefficient proportional to 1/SQRT(Ri):

      fakt1=4.e-3_wp*(hhl(i,k)-hhl(i,ke1)) !precalculation including the tuning-factor
      fakt2=0.25_wp+7.5e-3_wp*(hhl(i,k)-hhl(i,ke1))
      fakt1=MERGE( fakt1, fakt2, (gz0(i) < 0.01_wp .AND. l_pat(i) > 0._wp) ) !special treatment of glaciers
      fakt1=tkred_sfc(i)*fakt1; fakt2=tkred_sfc_h(i)*fakt1
      val1=val1*MIN( 2.5_wp, MAX( 0.01_wp, MIN( z1, fakt1 )*xri(i,k) ) )
      val2=val2*MIN( 2.5_wp, MAX( 0.01_wp, MIN( z1, fakt2 )*xri(i,k) ) )

      IF (tdc%imode_tkemini == 3) THEN !restricted (tuned) adaptation of TKE and the TMod. to lower limits
        tke(i,k,ntur)=tke(i,k,ntur)*SQRT(MAX( z1, val2/tkvh(i,k)) ) !adapted 'tke'

        !Note:
        !In contrast to the full 'tke'-adaptation at "ltkeadapt=T" (just below), the enhanced stratospheric LLDCs
        ! are not included into the 'tke'-adpatation, and only the "SQRT" of the derived amplification factor is applied.
        !So far, both modifications appear to be necessary with the current state of the entire model.
      END IF

      IF (tdc%tkhmin_strat > z0 .OR. tdc%tkmmin_strat > z0) THEN
        ! Enhanced diffusion in the stratosphere - needed primarily for momentum because
        ! there is otherwise too little dynamic coupling between adjacent model levels

        fakt = MIN( z1, 2.e-4_wp*MAX( z0, hhl(i,k) - 12500._wp ) ) ! transition zone between 12.5 and 17.5 km
        ! Wider transition zone in the tropics in order to avoid too strong diffusion in the tropopause region
        x4  =z1-z1d3*trop_mask(i)*MIN(z1, 2.e-4_wp*MAX(z0, 22500._wp-hhl(i,k)) )
        x4i =z1-z2d3*innertrop_mask(i)*MIN(z1, 2.e-4_wp*MAX(z0, 27500._wp-hhl(i,k)) )
        wert=SQRT(xri(i,k))
        fakt=fakt*MIN( x4*1.5_wp, MAX( 0.25_wp, wert ) )
        val1=MAX( val1, tdc%tkmmin_strat*MIN(x4,x4i)*fakt ); val2=MAX( val2, tdc%tkhmin_strat*x4*fakt )
      END IF

      ! Remark (GZ): The enhanced stratospheric diffusion seems to parameterize a missing
      !              process outside the turbulence scheme, maybe momentum transports due
      !              to non-stationary gravity waves. This may also explain why we need a
      !              much larger minimum diffusion coefficient for momentum than for heat.

      !Note:
      !This kind of correction by LLDCs compensates shortcomings of other components of the model
      ! and should be substituted by a less ad-hoc approach by means of respective process descriptions!
      !So far, this treatment degenerates the turbulence model in its principal conception according to STIC!
      !The now implemented "adapted treatment of lower limits" by means of the turb. Prandtl.-number 'tprn'
      ! (activated at "imode_tkemini=2") may circumvent this problem, so as to avoid an unintended feed-back
      ! of an artificial 'tkvm'-contribution onto the turbulence scheme.
      !IF the overall effect of LLDCs appears too large now, the LLDC for scalars should be reduced accordingly!
    END IF !using addit. empirical modificat. by means of 'xri'
!---------------------------------------------------------------------------
!GZ> For tuning

    IF (ltkeadapt) THEN !full adaptation of TKE and the TMod. to lower limits
       tke(i,k,ntur)=tke(i,k,ntur)*MAX( z1, val2/tkvh(i,k) ) !adapted 'tke'

       !Note:
       !In contrast to the amplification of 'tkvm' by lower limits, this LLDC-amplification of 'tkvh' is assumed
       ! to be always only due to otherwise not considered shear by NTCs.
    END IF

    tkvm(i,k)=MAX( val1, tkvm(i,k) ) !'tkvm' with total lower limit
    tkvh(i,k)=MAX( val2, tkvh(i,k) ) !'tkvh' with total lower limit

!---------------------------------------------------
  END DO  ! i
  END DO  ! k
  !$ACC END PARALLEL
!---------------------------------------------------

  IF (lsrfshear) THEN !specific preparations for surface-layer adaptations to additional "P"-level shear
!DIR$ IVDEP
     !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
     !$ACC LOOP GANG VECTOR PRIVATE(fakt)
     DO i=ivstart, ivend
        fakt=tkvh(i,ke)/tfh(i) !amplification-factor for the additional shear-forcing introduced to 'tkv[m|h](:,ke)' by LLDCs
        tfh(i)=ediss(i,ke)*(fakt**2-z1)/tfm(i) !additional shear-forcing corresponding to the impact of LLDCs
        tfm(i)=fakt*(tfm(i)/tkvm(i,ke)) !reduction-factor for the pure drag-contribution that has been introduced
                                        ! to 'tkvm(:,ke)' by LLDCs

        !Note:
        !So far, 'tf[m|h]' have contained the DCs at the lowest half level "P" without any lower limit.
        !In order to express the partial derivation of the DCs as a function of shear-forcing, the equilibrium solution
        ! is employed, at which the implicit deviation of stability-functions is neglected.
     END DO
     !$ACC END PARALLEL
  END IF

  IF (l3dturb) THEN !horizontal diff. coeffs. required
     !Consider horizontal diffusion coefficients:

     !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
     !$ACC LOOP SEQ
     DO k=2,kem
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR PRIVATE(wert)
        DO i=ivstart, ivend
           wert=tprn(i,k)*tkvh(i,k) !true turb. diff. coeff. for mom. relieved from artif. drag contrib.
           IF (PRESENT(hdef2) .AND. PRESENT(hdiv) .AND. ltkeshshr) THEN
              !Add isotropic turbulent part to that part due to the sep. horiz. shear mode:
              tkhh(i,k)=tkhh(i,k)+tkvh(i,k)
              tkhm(i,k)=tkhm(i,k)+wert
           ELSE !no treatment of sep. horiz. shear mode has taken place
              !Load only the isotropic turbulent part:
              tkhh(i,k)=tkvh(i,k)
              tkhm(i,k)=wert
           END IF
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

!------------------------------------------------------------------------------------------------
  IF (iini.NE.1) THEN !not for the separate initialization before the time loop
    !(this block lasts until the end of SUB 'turbdiff')
!------------------------------------------------------------------------------------------------

    ldoexpcor=(tdc%lexpcor .AND. tdc%icldm_turb.NE.-1) !consider explicit warm-cloud correct. for turb. scalar fluxes
    ldocirflx=(tdc%lcirflx .AND. lcircterm) !consider circulation heat-flux

!------------------------------------------------------------------------------------
! 5)  Berechnung der zu TKE-Quellen gehoerigen Temperaturtendenzen
!     (ausser der Divergenz des Zirkulationsflusses):
!------------------------------------------------------------------------------------

    IF (tdc%ltmpcor) THEN
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(phasdif)
      DO k=2, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          ! Beachte:
          ! In 'hlp' steht an dieser Stelle der Temp.-Gradient (auf Neben-Flaechen)!
          IF (tdc%lcpfluc) THEN
             phasdif=tkvh(i,k)*hlp(i,k)*(tur_rcpv*zvari(i,k,vap)+tur_rcpl*zvari(i,k,liq))
             hlp(i,k)=ediss(i,k)/cp_d + phasdif
          ELSE
             hlp(i,k)=ediss(i,k)/cp_d
          END IF
          hlp(i,k)=len_scale(i,k)/zaux(i,k,2)*hlp(i,k)
        END DO

        ! Achtung:
        ! Wegen der nachfolgenden Interpolation auf Hauptflaechen muss 'hlp' auch fuer "k=ke1" ausgewertet
        !  werden und ist jetzt eine mit der turbulenten Laengenskala multiplizierte T-Tendenz.
        ! Allein schon wegen dieser Interpolation ist wohl keine numerische Energie-Erhaltung zu erwarten!
      END DO
      !$ACC END PARALLEL

!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        t_tens(i,1)=t_tens(i,1)+hlp(i,2)/(len_scale(i,1)+len_scale(i,2))
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,ke
!DIR$ IVDEP
        DO i=ivstart, ivend
          t_tens(i,k)=t_tens(i,k)+(hlp(i,k)+hlp(i,k+1)) &
                                         /(len_scale(i,k)+len_scale(i,k+1))
        END DO
      END DO
      !$ACC END PARALLEL

    END IF !ltmpcor

    IF (ldocirflx) THEN !consider heat-flux by thermal SSO circulations
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,ke
!DIR$ IVDEP
        DO i=ivstart, ivend
          hlp(i,k)=ABS(frh(i,k)*tkvh(i,k)) !save magnitude of updated buoyant TKE-source by turbulence
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

!------------------------------------------------------------------------------------
! 6) Berechnung der Diffusionstendenz von q=SQRT(2*TKE) einschliesslich der
!    q-Tendenz durch den Zirkulationsterm:
!------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! Vorbereitung zur Bestimmung der zugehoerigen Incremente von TKE=(q**2)/2:
    !--------------------------------------------------------------------------------

    upd_prof => zaux(:,:,1)
    sav_prof => zaux(:,:,2)

    IF (ldotkedif .OR. lcircterm) THEN
       !ldotkedif=(tdc%c_diff > 0.0): partly implicit vertical diffusion for TKE
       !lcircterm: raw "circulation term" computed together with TKE diffusion

      expl_mom => zaux(:,:,3)

      ! Diffusions-Koeffizienten auf NF:
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          IF (tdc%imode_tkediff == 2) THEN !Diffusion in terms of TKE
!___________________________________________________________________________
!test: TKE-Diffusion mit Stab.fnkt. fuer Skalare: <
            sav_prof(i,k)=tdc%c_diff*len_scale(i,k)*tke(i,k,ntur) !diff.-coeff. for TKE
! sav_prof(i,k)=tdc%c_diff*tkvh(i,k)
!test>
!___________________________________________________________________________

          ELSE !Diffusion in terms of q=SQRT(2*TKE)
            sav_prof(i,k)=tdc%c_diff*len_scale(i,k)*tke(i,k,ntur)**2 !diff.coeff. for q
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=3, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          expl_mom(i,k)=rhoh(i,k-1)*z1d2*(sav_prof(i,k-1)+sav_prof(i,k)) &
                                             /(hhl(i,k-1)-hhl(i,k))
        END DO
        !Beachte:
        !'expl_mom' bezieht sich auf HF, also die Fluss-Niveaus fuer die TKE (bzw. q-)-Diffusion.
        !Wegen der spaeteren Nutzung der SUBs 'prep_impl_vert_diff' und 'calc_impl_vert_diff'
        ! muss ein Fluss-Niveau (hier HF) ueber dem Variabl.-Niveau (hier NF) mit gleichem Index liegen.
      END DO
      !$ACC END PARALLEL

      IF (tdc%imode_tkediff == 2) THEN !Diffusion in terms of TKE
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2, ke1
!DIR$ IVDEP
          DO i=ivstart, ivend
            sav_prof(i,k)=z1d2*tke(i,k,ntur)**2 !TKE
          END DO
        END DO
        !$ACC END PARALLEL
      ELSE !Diffusion in terms of q=SQRT(2*TKE)
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2, ke
!DIR$ IVDEP
          DO i=ivstart, ivend
            sav_prof(i,k)=tke(i,k,ntur)         !q=SQRT(2*TKE)
            dicke(i,k)=dicke(i,k)*tke(i,k,ntur) !related effective discretization momentum
          END DO
        END DO
        !$ACC END PARALLEL

        !Das Feld 'dicke' wird bei "k=ke1" nicht benoetigt und war zuvor dort auch nicht belegt!

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
        DO i=ivstart, ivend
          sav_prof(i,ke1)=tke(i,ke1,ntur) !q at the surface layer
        END DO
        !$ACC END PARALLEL
      END IF

      !Beachte:
      !Das Feld 'tke' enthaelt nicht TKE sondern "q=SQRT(2TKE)"!

    END IF !(ldotkedif .OR. lcircterm)

    !--------------------------------------------------------------------------------
    ! Aufnahme des Zirkulationstermes mit Interpolation auf HF:
    !--------------------------------------------------------------------------------

    IF (lcircterm .OR. loutthcrc) THEN !Der bisherige "Zirkulationsterm" muss berechnet werden

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          ! Belegung von 'frh' mit der CKE-Flussdichte durch nicht-turbulente Zirkulationen, die durch thermische
          !  Inhomogenitaet an der Oberflaeche verursacht wird:

          frh(i,k)=rhon(i,k)*tkvh(i,k)*prss(i,k)*len_scale(i,k) !skalierte Flussdichte auf NF

          !Bemerkung:
          !'frh/len_scale' ist eine TKE-Flussdichte in [Kg/s3], deren Vertikalprofil im wesentlichen durch
          ! d_z(tet_v)**2 bestimmt ist, was zumindest in der Prandtl-Schicht prop. zu 1/len_scale ist.
          !Die nachfolgende lineare Interpolation auf Hauptflaechen erfolgt daher mit 'frh'!
        END DO
      END DO
      !$ACC END PARALLEL

      ! Interpolation der skalierten CKE-Flussdichte auf Hauptflaechen:

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=3,ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          frm(i,k)=(frh(i,k)+frh(i,k-1))/(len_scale(i,k)+len_scale(i,k-1)) !interpolierte Flussdichte auf HF
        END DO
      END DO
      !$ACC END PARALLEL

      !----------------------------------------------------------------------------------------------
      ! Explizite Berechnung des Zirkulationsterms und der Konvergenz des zugehoerigen Theta-Flusses:
      !----------------------------------------------------------------------------------------------

      IF (tdc%lcirflx .OR. loutthcrc) THEN !explicit TKE-source of raw "circulation-term" required

!DIR$ IVDEP
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          upd_prof(i,2)=frm(i,3)/dicke(i,2)
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=3,ke
!DIR$ IVDEP
          DO i=ivstart, ivend
            upd_prof(i,k)=-(frm(i,k)-frm(i,k+1))/dicke(i,k)
          END DO
        END DO
        !$ACC END PARALLEL

        !Note:
        !At "imode_tkediff == 2", 'upd_prof' is a profile of TKE-increments. In contrast, at
        !   "imode_tkediff == 1", 'dicke' contains the additional factor 'tke', and
        !                         'upd_prof' is a profile or q-increments.

        ! Explizite Berechnung der zum alten "Zirkulationsterm" gehoerigen CKE-Quelle:

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP SEQ
        DO k=2,ke !no Theta-flux due to thermal circulations at the surface
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR
          DO i=ivstart, ivend

            shv(i,k)=upd_prof(i,k)*fr_tke
            IF (tdc%imode_tkediff == 1) THEN !'upd_prf' conatins q-increments
               shv(i,k)=shv(i,k)*tke(i,k,ntur)
               !Note that 'upd_prf' contains the factor "1/dicke", where 'dicke' includes the factor 'tke(:,:,ntur)'.
            END IF

            IF (loutthcrc) THEN !TKE-source by raw "circulation-term" to be stored as STIC-term
               tket_nstc(i,k) = shv(i,k) !load TKE-source considered to be due to thermal SSO production
            END IF

          END DO
        END DO
        !$ACC END PARALLEL

      END IF !explicit TKE-source of raw "circulation-term" required

    END IF !Der bisherige "Zirkulationsterm" muss berechnet werden

!------------------------------------------------------------------------------------------------------
! 7) Addition des Theta-Gradienten, welcher zur nicht-turbulenten Theta-Flussdichte durch den
!     alten "Zirkulationsterm" oder durch thermische SSO  bzgl. des neuen "Zirkulationsterms" gehoert:
!------------------------------------------------------------------------------------------------------

    IF (ldocirflx) THEN !consider circulation heat-flux

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=2,ke !no Theta-flux due to thermal circulations at the surface
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR PRIVATE(wert)
        DO i=ivstart, ivend
           wert=shv(i,k) !(possibly reduced) buoyant heat-flux of NTCs
           wert=SIGN(z1,wert)*MIN(hlp(i,k),ABS(wert)) !limitation of buoyant heat-flux by NTCs through
                                                      ! the magnitude of buoyant heat-flux by turbulence
           wert=-wert/(zaux(i,k,4)*tkvh(i,k)) !virtual gradient of the related circul. Theta-flux

           If (ldoexpcor) THEN !converted explicit turbulent fluxes required for vert. diff.
              zvari(i,k,tet) = zvari(i,k,tet) + wert !store increased effective Theta-gradient
           ELSE !converted explicit turbulent fluxes are not being used for vertical diffusion
              zvari(i,k,tet) = wert !store only the circulation-contribution as Theta-gradient!!
           END IF
        END DO
      END DO
      !$ACC END PARALLEL

    END IF !consider circulation heat-flux

!------------------------------------------------------------------------------------
! 8)  Bestimmung des Zirkulationstermes als zusaetzliche TKE-Flussdichte:
!------------------------------------------------------------------------------------

    IF (lcircterm) THEN !Der bisherige "Zirkulationsterm" geht in die TKE-Gleichung ein

      ! Berechnung der TKE-Flussdichte-Konvergenz (einchliesslich der von CKE) und anderer TKE-Quellen
      !  fuer den naechsten Zeitschritt:

      ! Quasi-implizite Berechnung der Zirkulatinstendenz:

      cur_prof => hlp

!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        cur_prof(i,2)=sav_prof(i,2)
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=3,ke1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR !to be activated
        DO i=ivstart, ivend
          cur_prof(i,k)=(cur_prof(i,k-1)-sav_prof(i,k-1)+frm(i,k)/expl_mom(i,k))+sav_prof(i,k) !to be activated
        END DO
      END DO
      !$ACC END PARALLEL

      !Beachte:
      !'cur_prof' enthaelt ein virtuelles TKE-Profil (oder q-Profile bei "imode_tkediff=1"),
      ! dessen Diffusions-Tendenz die Zirkulations-Tendenz einschliesst.
      !Fuer die expliziten Diff.-Anteile wird ebenfalls 'cur_prof' benutzt.

      !Bereucksichtige Zirkulations-Tendenz:
      itndcon=0 !indem 'cur_prof' auf der rechten Seite der impliz. Diff.-Gl. benutzt wird.

    ELSEIF (ldotkedif) THEN

      cur_prof => sav_prof

      itndcon=0 !'cur_prof' wird auf der rechter Seite der impliz. Diff.-Gl.
                ! und fuer explizite Diff.-Anteile benutzt.

    END IF !Der bisherige "Zirkulationsterm" geht in die TKE-Gleichung ein

    IF (ldotkedif .OR. lcircterm) THEN

!----- --------------------------------------------------------------------
! 9)  Aufdatieren des TKE-Profils durch die (erweiterte) Diffusions-Tendenz
!--------------------------------------------------------------------------

      impl_mom => zaux(:,:,4)
      invs_mom => zaux(:,:,5)

      !'frm', 'frh' und 'len_scale' sind frei.
      !In den Diffusionroutinen wird vorausgesetzt, dass ein Flussniveau mit gleichem
      !Vertikalindex wie ein Konzentrationsniveau gerade ueber letzterem liegt.
      !Die bisherige Hauptflaechenindizierung musste daher fuer die Uebergabefelder
      !der Routinen 'prep_impl_vert_diff' und 'calc_impl_vert_diff' um "1" verschoben werden.

      CALL prep_impl_vert_diff( tdc=tdc,                                                   &
                                lsflucond=.FALSE., ldynimpwt=ldynimp, lprecondi=lprecnd,   &
                                i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1,                &
!Achtung: q_Diff:
!disc_mom=sav_prof*dicke,    expl_mom=expl_mom,                      &
                                disc_mom=dicke,    expl_mom=expl_mom,                      &
                                impl_mom=impl_mom, invs_mom=invs_mom,                      &
                                invs_fac=frh, scal_fac=frm,                                &
                                lacc=lzacc )

      ! Berechnung der vertikalen Diffusionstendenzen von TKE=z1d2*q**2:

      eff_flux => len_scale

      CALL calc_impl_vert_diff ( lsflucond=.FALSE.,lprecondi=lprecnd,                      &
                                 leff_flux=(kcm <= ke), itndcon=-itndcon,                  &
                                 i_st=ivstart, i_en=ivend,k_tp=1, k_sf=ke1,                &
!Achtung: q_Diff:
!disc_mom=sav_prof*dicke,    expl_mom=expl_mom,                      &
                                 disc_mom=dicke,    expl_mom=expl_mom,                     &
                                 impl_mom=impl_mom, invs_mom=invs_mom,                     &
                                 invs_fac=frh, scal_fac=frm, cur_prof=cur_prof,            &
                                 upd_prof=upd_prof, eff_flux=eff_flux, lacc=lzacc )

      !Beachte:
      !Bei "imode_tkediff=1" erfolgt q-Diffusion, so dass 'disc_mom' und 'expl_mom' den zusaetzlichen
      ! Faktor 'sav_prof'='tke'='q' (turb. Geschwindigkeit) enthalten!
      !'upd_prof' enthaelt jetzt die mit der Diffusionstendenz aufdatierten (modifizierte) TKE-Werte.
      !Weil "itndcon<=0", bleiben auf 'cur_prof' die Eingangsprofile erhalten.
      !'eff_flux' enthaelt die effektiven Flussdichten (positiv abwaerts) der (semi-)impliziten
      ! Vertikaldiffusion.

      IF (lcircterm) THEN !es wurden virtuelle Effektiv-Profile benutzt
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,ke
!DIR$ IVDEP
          DO i=ivstart, ivend
            upd_prof(i,k)=sav_prof(i,k)+upd_prof(i,k)-cur_prof(i,k) !aufdatierte echte Profile
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      IF (PRESENT(r_air)) THEN
        ! Zuschlag durch Volumenterm innerhalb der Rauhigkeitsschicht:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP SEQ
        DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR PRIVATE(wert)
          DO i=ivstart, ivend
!           wert=(eff_flux(i,k)*dp0(i,k)+eff_flux(i,k+1)*dp0(i,k-1))/(dp0(i,k)+dp0(i,k-1))
            wert=zbnd_val(eff_flux(i,k+1), eff_flux(i,k), dp0(i,k), dp0(i,k-1))
                 !effektive TKE-Flussdichte interpoliert auf die k-te Nebenflaeche,
                 ! wobei sich 'eff_flux(:,k)' auf die (k-1)-te Hauptflaeche bezieht!
            upd_prof(i,k)=upd_prof(i,k)+wert*z1d2*(r_air(i,k-1)-r_air(i,k+1))/dicke(i,k)
          END DO
        END DO
        !$ACC END PARALLEL
        !'upd_prof' enthaelt das mit dem Volunmenterm beaufschlagte aufdatierte Profil.
      ENDIF

!------------------------------------------------------------------------------------
!10)  Speichern der zugehoerigen q-Tendenzen:
!------------------------------------------------------------------------------------

      IF (tdc%imode_tkediff == 2) THEN !Diffusion in terms of TKE
        !'upd_prof' ist ein TKE-Profil:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,ke
!DIR$ IVDEP
          DO i=ivstart, ivend
!___________________________________________________________________________
!test: different versions of effective prognositc variable tke=q=SQRT(2*TKE): <
                  ! direct tendency of q:
            tketens(i,k)=( SQRT( 2*MAX( upd_prof(i,k), z0 ) ) - tke(i,k,ntur) )*fr_tke
            ! Note: This way appears to be numerically stable!!

! ! according to Dt(TKE)=Dt(0.5*q**2)=q*Dt(q) => dD(q)=Dt(TKE)/q:
! tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - sav_prof(i,k) )*fr_tke/tke(i,k,ntur)
! !Attention: this way appears to be numerically unstable!!
!test>
!___________________________________________________________________________

          END DO
        END DO
        !$ACC END PARALLEL

      ELSE !Diffusion in terms of q=SQRT(2*TKE)
        !'upd_prof' ist ein q-Profil:

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,ke
!DIR$ IVDEP
          DO i=ivstart, ivend
            !Achtung:
            tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - tke(i,k,ntur) )*fr_tke
          END DO
        END DO
        !$ACC END PARALLEL
      END IF
      !'tketens' enthaelt jetzt immer eine q-Tendenz!

      !Am Unterrand gibt es keine q-Tendenz durch Diffusionsterme:
!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        tketens(i,ke1)=z0
      END DO
      !$ACC END PARALLEL

    ELSE !keine q-Tendenzen, weder durch TKE-Diffusion noch durch den Zirkulationsterm

      ! Zuruecksetzen der q-Tendenzen:
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          tketens(i,k)=z0
        END DO
      END DO
      !$ACC END PARALLEL

    END IF   ! ldotkedif .OR. lcircterm

!------------------------------------------------------------------------------------
! 11) Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!     des Saettigungsdefizites:
!------------------------------------------------------------------------------------

!DIR$ IVDEP
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      rcld(i,1)=rcld(i,2)
    END DO
    !$ACC END PARALLEL

#ifdef __INTEL_COMPILER
    FORALL(k=2:kem-1, i=ivstart:ivend) &
        rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
#else
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=2,kem-1
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
      END DO
    END DO
    !$ACC END PARALLEL
#endif

    ! Fuer die unterste Hauptflaeche (k=ke) wird bei kem=ke
    ! der Wert auf der entspr. Nebenflaeche beibehalten.

!------------------------------------------------------------------------------------------------
  END IF !not for the separate initialization before the time loop
!------------------------------------------------------------------------------------------------

  !$ACC WAIT
  !$ACC END DATA

END SUBROUTINE turbdiff

!==============================================================================

END MODULE turb_diffusion
