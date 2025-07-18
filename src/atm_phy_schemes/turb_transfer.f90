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

! Source module for computing the coefficients for turbulent transfer
!
! Description of *turb_transfer*:
!   This  module calculates the coefficients for turbulent transfer.
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
!   turbtran
!
! called from the turbulence interface routine of the model.
!
!-------------------------------------------------------------------------------

MODULE turb_transfer

!-------------------------------------------------------------------------------
!
! Documentation History:
!
!  The history of all these modifications is as follows, where those belonging to the formal
!   reorganization of the whole package (atmospheric turbulence and surface-to-atmosphere transfer)
!   are now in the header of MODULE 'turb_utilities', containing various common SUBs for 'turbdiff'
!   and 'turtran' (related moist thermodynamicds and the treatment of turbulent budget equations)
!   and also the blocked code for semi-implicit vertical diffusion. The new blocked version of SUB 'turbtran'
!   is now in MODULE 'turb_transfer':
!
!              2010/09/30 Matthias Raschendorfer
!  Substitution of 'itype_diag_t2m' by 'itype_synd' being already used for that purpose
!   "1": alternative SYNOP-digansostics related to previous Lewis scheme (not included here)
!   "2": SYNOP-diagnostics according to current transfer scheme using SYNOP-z0 'z0d'
!   "3": like "2" but using 'z0d' only for 10m-wind and a specific roughness layer profile for
!        2m-temperature and -humidity.
!  Including the adiabatic lapse rate correction for 't_2m' for "itype_synd = 3" as well.
!              2011/03/23 Matthias Raschendorfer
!  Correction of two bugs introduced together with the last modifications in SUB 'turbtran'
!   (related to SUB 'diag_level' and the 'tet_l'-gradient used for the flux output).
!  Substitution of run time allocations because of bad performance on some computers.
!              2014/07/28 Matthias Raschendorfer
!  Removing a bug in formular for 'rcld(:,ke1)' in SUB 'turbtran'
!   -> influence on near-surface temperature and - humidity
!              2015/08/25 Matthias Raschendorfer
! Adopting other development within the ICON-version (namely by Guenther Zaengl) as switchable options
!  related to the following new selectors and switches:
!   imode_rat_sea, imode_vel_min, imode_charpar and lfreeslip.
! Rearranging the development by Matthias Raschendorfer that had not yet been transferred to COSMO as switchable
!  options related to the following switches:
!  and selectors:
!   itype_diag_t2m, imode_syndiag, imode_trancnf, imode_tkemini, imode_lamdiff
!  and a partly new (more consistent) interpretation of:
!   imode_tran, icldm_tran, itype_sher, and itype_diag_t2m
! Controlling numerical restrictions gradually namely by the parameter:
!  ditsmot
! Using the arrays 'tvm', 'tvh' and 'tkm', allowing an easier formulation of transfer-resistances.
!              2016-05-10 Ulrich Schaettler
! Splitting this module from the original module 'organize_turbdiff' as it was used by ICON before.
! Moving declarations, allocation and deallocations of ausxilary arrays into MODULE 'turb_data'.
!
!-------------------------------------------------------------------------------

! Modules used:

#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num
#endif
USE mo_exception,         ONLY: message_text, message
!-------------------------------------------------------------------------------
! Parameter for precision
!-------------------------------------------------------------------------------

USE mo_kind,         ONLY :   &
    wp              ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------
! Mathematical and physical constants
!-------------------------------------------------------------------------------

USE mo_mpi,                ONLY : get_my_global_mpi_id

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
    t0_melt  => tmelt,    & ! absolute zero for temperature (K)
    b3       => tmelt,    & !          -- " --

    rdv, con_m, con_h, grav

USE mo_lookup_tables_constants, ONLY : &
!
! Parameters for auxilary parametrizations:
! ------------------------------------------
!
    b1       => c1es,     & ! variables for computing the saturation steam pressure
    b2w      => c3les,    & ! over water (w) and ice (e)
    b4w      => c4les       !               -- " --

!-------------------------------------------------------------------------------
! Turbulence data (should be the same in ICON and COSMO)
!-------------------------------------------------------------------------------

USE mo_turbdiff_config, ONLY : &

    t_turbdiff_config, &

! Numerical constants and parameters:
! -----------------------------------

    ! derived parameters calculated in 'turb_setup'
    tet_g, rim, b_m, b_h, sm_0, sh_0,   &
    a_3, a_5 ,a_6,                      &
    tur_rcpv, tur_rcpl,                 &

    ! used derived types
    modvar,       & !

    ltst2ml,  &     ! test required, whether  2m-level is  above the lowest main-level
                    ! F:  2m-level is assumed to be always below the lowest main-level
    ltst10ml, &     ! test required, whether 10m-level is  above the lowest half-level
                    ! F: 10m-level is assumed to be always below the lowest half-level

    ! numbers and indices

    ninv    ,     & ! number of invariant scalar variables being conserved during 'vap'<->'liq'-transistions ('tet_l', 'h2o_g')

    nvel    ,     & ! number of velocity-components active for turbulece ('u_m', 'v_m')
    nmvar   ,     & ! number of progn. variables beding active for turbulence
    ndim    ,     & ! (positive) limit of last dimension used for 'zaux' and 'zvari'
    mom, sca,     & ! indices for momentum- and scalar-variables
    u_m     ,     & ! zonal velocity-component at the mass center
    v_m     ,     & ! meridional ,,      ,,    ,, ,,   ,,    ,,
    tet_l   ,     & ! liquid-water potential temperature
    h2o_g   ,     & ! total water ('qv+qc')
    tet     ,     & ! potential temperature
    vap     ,     & ! water vapor (specific humidity 'qv')
    liq             ! liquid water (mass fraction 'qc')

    !Note: It always holds: tem=tet=tet_l=tem_l and vap=h2o_g (respective usage of equal indices)!


!-------------------------------------------------------------------------------
! Control parameters for the run
!-------------------------------------------------------------------------------

! Switches controlling other physical parameterizations:
USE mo_lnd_nwp_config,       ONLY: lseaice, llake, lterra_urb, itype_kbmo
USE mo_atm_phy_nwp_config,   ONLY: lcuda_graph_turb_tran
!
USE turb_utilities,          ONLY:   &
    turb_setup,                      &
    adjust_satur_equil,              &
    solve_turb_budgets,              &
    zexner, zpsat_w


!-------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
    lsclm, lsurflu, i_cal, i_upd, i_mod, imb, &
    SHF, LHF
#endif
!SCLM---------------------------------------------------------------------------

USE mo_fortran_tools, ONLY: set_acc_host_or_device
!===============================================================================

IMPLICIT NONE

PUBLIC  :: turbtran

!===============================================================================

!-------------------------------------------------------------------------------

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

SUBROUTINE turbtran (                                                         &
!
          tdc,                                                                &
          iini, ltkeinp, igz0inp, lsrflux, lnsfdia, lrunscm,                  &
          ladsshr,                                                            &
!
          dt_tke, nprv, ntur, ntim,                                           &
!
          nvec, ke, ke1, kcm, iblock, ivstart, ivend,                         &
!
          l_pat, l_hori, hhl, fr_land, l_lake, l_sice, gz0, z0_waves,         &
          rlamh_fac, sai, urb_isa,                                            &
          t_g, qv_s, ps, u, v, t, qv, qc, epr,                                &
!
          tcm, tch, tvm, tvh, tfm, tfh, tfv, tkr,                             &

          tke, tkvm, tkvh, rcld,                                              &
          hdef2, dwdx, dwdy,           & ! optional for "itype_sher=2"
!
          edr, tketens,                                                       &
!
          t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m,                            &
          shfl_s, qvfl_s, umfl_s, vmfl_s,                                     &
!
          lacc, opt_acc_async_queue)
!-------------------------------------------------------------------------------
!
! Note:
! It is also possible to use only one time level for TKE using "ntim=1" and thus "nprv=1=ntur".
!
! Description:
!
!     Es werden die Transfer-Geschwindigkeiten fuer den Austausch von Impuls,
!     sowie fuehlbarer und latenter Waerme bestimmt und die Modellwerte
!     fuer die bodennahen Messwerte (in 2m und 10m) berechnet.
!
! Method:
!
!     Hierzu wird der gesamte Bereich von den festen Oberflachen am
!     Unterrand des Modells bis hin zur untersten Hauptflaeche in
!     die drei Teilbereiche:
!
!     - laminare Grenzschicht (L-Schicht)
!     - turbulente Bestandesschicht (B-Schicht)
!     - turbulente Prandtl-Schicht (P-Schicht)
!
!     aufgeteilt. Fuer jeden dieser Teilbereiche wird (getrennt nach
!     skalaren Eigenschaften und Impuls) ein zugehoeriger Transport-
!     widerstand berechnet, der gleich einer effektiven Widerstands-
!     laenge ( dz_(sg, g0, 0a)_(h,m) ) dividiert durch den Diffusions-
!     koeffizienten am Unterrand der P-Schicht (Niveau '0') ist.
!     Die Konzentrationen am Unterrand der B-Schicht, also im
!     Abstand der L-Schicht-Dicke entlang der festen Oberflaechen,
!     haben den Index 'g' (ground) und die Oberflaechenkonzentrationen
!     den Index 's'. Groessen fuer den Imoulst haben den Index 'm' (momentum)
!     und solche fuer skalare Eigenschaften 'h' (heat).
!     Der Widerstand der P-Schicht vom Niveau '0' bis zum Niveau 'a'
!     (atmospheric) der untersten Hauptflaeche wird durch vertikale
!     Integration der Modellgleichungen in P-Schicht-Approximation
!     (vertikal konstante Flussdichten, turbulente Laengenskala lin. Funkt.
!      von der Hoehe) gewonnen.
!     Dabei wird das atmosphaerische Turbulenzschema aus der Subroutine
!     'turbdiff' benutzt, so dass alo keine empirischen Profilfunktionen
!     benutzt werden. Zur Vereinfachung der Integration  wird das Produkt
!     aus turbulenter Geschwindigkeitsskala 'q' und der Stabilitaetsfunktion
!     s(h,m), alo die stabilitaetsabhaengige turb. Geschwindigkeitsskala
!     'v' innerhalb der P-Schicht als linear angesehen.
!     Die turb. Laengenskala im Niveau '0' wird mit der Rauhigkeitslaenge 'z0'
!     (multipliziert mit der v.Kaman-Konstanten) gleichgesetzt. Formal werden
!     dann fuer das Nieveau '0' Vertikalgradienten und auch Diffusions-
!     koeffizienten abgeleitet.
!     Unter der Annahme, dass 'v' innerhalb der B-Schicht konstant bleibt,
!     ergibt sich die laminare Widerstandslaenge dz_sg als prop. zu 'z0'
!     und die Widerstandsstrecke durch die B-Schicht als prop. zu
!     'z0*ln(delta/z0)', wobei 'delta' die Dicke der L-Schicht ist, die der
!     Abstand von einer ebenen Wand sein soll in dem der turbulente
!     Diffusionskoeffizient fuer Impuls gleich dem molekularen ist.
!     Ferner wird angenommen, dass die Widerstaende durch die L- und
!     B-Schicht prop. zur effektiven Quellflaech der Bestandeselemente
!     zuzueglich der Grundflaeche des Erdbodens sind. Die Bestandesoberflaechen
!     werden durch den Wert 'sai' (surface area index) ausgedrueckt und setzt
!     sich aus dem Flaechenindex der transpirierenden Oberflaechen 'lai'
!     (leaf area index) und dem fuer die nicht transpirierenden Flaechen
!     zusammen. Im Falle nicht benetzter Oberlfaechen hat die latente Waerme
!     i.a. eine kleinere Quellflaeche als die fuehlbare Waerme, so dass die
!     Wiederstaende fuer beide Groessen unterschieden werden muessten.
!     Um dies zu vermeiden, wird nur der Widerstand fuer die fuehlbare Waerme
!     berechnet. Dafuer wird aber bei der Berechnung der effektiven
!     Oberflaechenkonzentration 'qv_s' der spez. Feuchtigkeit in Subroutine
!     'terra1' dieser Effekt beruecksichtigt.
!     Beim vertikalen Impulstransport ist aber noch die zusaetzliche
!     Impulssenke innerhalb der B-Schicht durch die Wirkung der Formreibungs-
!     kraft zu beruecksichtigen, was durch einen zusaetzlichen Flaechenindex
!     'dai' (drag area index) bewerkstelligt wird.
!
!     Die Vertikalprofile aller Eigenschaften innerhalb der P-Schicht ergeben
!     sich aus dem vertikal integrierten Turbulenzmodell in P-Schicht-
!     Approximation zu logarithmischen Funktionen, welche durch die
!     thermische Schichtung modifiziert sind. Wie bereits erwaehnt, ist die
!     Stabilitaetsfunktion nur noch von Konstanten des atmosphaerischen
!     Turbulenzmodells abhaengig. Das Transferschema ist somit auch automatisch
!     konsistent zum oben anschliessenden Turbulenzmodell formuliert.
!
!     Die Profilfunktionen innerhalb der B-Schicht ergeben sich aus der
!     Annahme eines Gleichgewichtes zwischen vertikalen Flussdichtedivergenzen
!     und Quellstaerken durch die laminaren Grenzschichten der Rauhigkeits-
!     elemente bei vertikal konstanten Bestandeseigenschaften zu exponentiellen
!     Funktionen. Durch die Bedingung eines glatten Ueberganges zwischen beiden
!     Profiltypen im Niveau '0' und der Bedingung, dass im Abstand einer
!     effektiven Bestandesdicke 'Hb' unterhalb des Nieveaus '0' die Bestandes-
!     profile in die Konzentration am Unterrand der B-Schicht (Niveau mit
!     Index 'g') uebergehen, ist das gesamte Transferschema geschlossen und
!     es kann auch der "drag area index" 'dai', sowie die Bestandeshoehe
!     'Hb' selbst eliminiert werden.
!
!     Zur Charakterisierung des Oberflaechentransfers werden dann nur die
!     externen Parameter 'z0', 'sai', 'lai' und je ein globaler Parameter
!     fuer den laminaren Grenzschichtwiderstand des skalaren - und des
!     Impulstransportes benoetigt '(lam_(h,m)'. Hieraus koennte auch eine
!     aequivalente Rauhigkeitslaenge fuer Skalare 'z0h' berechnet werden.
!     Die Oberfalaechenkonzentrationen (Niveau mit Index 's') fuer die skalaren
!     Groessen werden im Modul 'terra' berechnet. Fuer den Impuls gilt die
!     Haftbedingung. Im Grundniveau des atmosphaerischen Modells
!     (Niveau '0') verschwindet also der Wind i.a. nicht; dies ist erst
!     entlang der festen Oberfalechen der Fall. Die bodennahen synoptischen
!     Niveaus werden nun vom Niveau 'z=-Hb', also von der effektiven Bestandes-
!     grundflaeche (Umsatzniveau) aus gezaehlt. Ist z.B. 'Hb>2m', werden
!     die 2m-Werte entlang der exponentiellen Bestandesprofile ausgeweret.
!     Ist 'Hb<2m', wird das logarithmische Profil in der Hoehe '2m-Hb' entlang
!     dinnerhalb der P-Schicht ausgewertet.
!     Die resultierenden Transfer-Geschwindigkeiten 'tv(h,m)' sind die Kehrwerte
!     des Gesamtwiderstandes von den festen Oberflaechen (Neviau 's') bis
!     zur untersten Modellhauptflaeche (Niveau 'a').
!     Die turbulenten Diffusionskoeffizienten 'tkv(h,m)' fuer den vertikalen
!     Index 'ke1', beziehen sich aber auf den Unterrand des atmosphaerischen
!     Modells (Niveau '0').
!     Mit Hilfe der Felder 'tf(mh)' werden noch Reduktionsfaktoren der
!     Transfer-Geschwindigkeiten durch die Wirkung der L-Schicht uebertragen.
!     Diese koennen im Modul 'terra' benutzt werden, um ev. das effektive
!     'qv_s' so zu bestimmen, als gaebe es fuer fuehlbare und latente Waerme
!     unterschiedliche Parameter fuer den laminaren Transportwiderstand.
!     Zu beachten ist, dass im Falle eines vertikal vom atmosphaerischen Modell
!     aufgeloesten 'Makrobestandes' (z.B. Bebauung, Wald oder subskalige
!     Orographie) das Transferschema genauso wie im Falle eines nicht
!     aufgeloesten Bestandes angewendet wird. Allerdings beziehen sich die
!     den Bestand des Transferschemas charakterisierenden externen Parameter
!     dann auf den nicht vertikal aufgeloesten verbleibenden 'Mikrobestand',
!     der ev. allein durch niedrigen Bewuchs gebildet wird.
!     Im Transferschema eingearbeitet ist auch dei iterative Bestimmmung der
!     Rauhigkeitslaenge der Meeresoberflaeche gemaess einer modifizierten
!     Charnock-Formel, bei der die Wellenerzeugung bei verschwindenden
!     mittleren Wind mit hilfe der zur TKE ausgedrueckt wird.
!
!-------------------------------------------------------------------------------

! Declarations
!-------------------------------------------------------------------------------

!Formal Parameters:
!-------------------------------------------------------------------------------

! Parameters controlling the call of 'organize_turbdiff':

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

LOGICAL, INTENT(IN) :: &

   lnsfdia,      & !calculation of (synoptical) near-surface variables required
   lsrflux,      & !calculation of surface flux densities in 'trubtran'


   ltkeinp,      & !TKE present as input for time level 'ntur' at level "0" (k=ke1)
                   ! and also at level "P" (k=ke) in case of "lini=T"

   lrunscm,      & !a Single Column run (default: FALSE)
   ladsshr         !treatment of additional surface-shear by NTCs or LLDCs active: "rsur_sher>0"

REAL (KIND=wp), INTENT(IN) :: &

   dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER,        INTENT(IN) :: &

   igz0inp,      & !1: gz0 is present as input, 2: z0_waves is provided as input
   iini,         & !type of initialization (0: no, 1: separate before the time loop
                   !                             , 2: within the first time step)
   ntur,         & !current  time level of 'tke' valid after  prognostic incrementation
   nprv,         & !previous time level of 'tke valid before prognostic incrementation
   ntim            !number of 'tke' time levels

INTEGER,        INTENT(IN) :: &

! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------

    nvec,         & ! number of grid points in the nproma-vector
    ke,           & ! index of the lowest model level: either Half-Level (HL) "P" or Main-Level (ML) "A"
    ke1,          & ! index of the lowest boundary/half-level: either surface- or "0"-level (k=ke+1)
    kcm,          & ! level index of the upper vertically-resolved canopy bound
    iblock

INTEGER,        INTENT(IN) :: &

! Start- and end-indices for the computations in the horizontal layers:
! -----------------------------------------------------------------------

    ivstart,      & ! start index in the nproma vector
    ivend           ! end index in the nproma vector

! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! --------------------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: &
!
    hhl             ! height of model half levels                   ( m )

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
    l_pat,        & ! effective length scale of near-surface circulation patterns [m]
                    !  (scaling the near-surface circulation acceleration)
    l_hori,       & ! horizontal grid spacing (m)
    rlamh_fac,    & ! scaling factor for rlam_heat
!
! External parameter fields:
! ----------------------------
    fr_land,      & ! land portion of a grid point area             ( 1 )
    sai,          & ! surface area index                            ( 1 )
    urb_isa         ! urban impervious surface area                 ( 1 )

! Fields for surface values and soil|canopy model variables:
! ------------------------------------------------------------

LOGICAL, DIMENSION(:), INTENT(IN) :: &
!
    l_lake,       & ! a lake surface
    l_sice          ! an ice surface

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
!
    ps,           & ! surface pressure                              ( pa  )
    qv_s,         & ! specific water vapor content on the surface   (kg/kg)
    t_g             ! weighted surface temperature                  (  k  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
! Atmospheric model variables:
! ---------------------------------
                    ! main-level values of:
     u,           & ! zonal wind speed       (at mass positions)    ( m/s )
     v,           & ! meridional wind speed  (at mass positions)    ( m/s )
     t,           & ! temperature                                   (  k  )
     qv,          & ! specific water vapor content                  (kg/kg)
     qc             ! specific cloud water content                  (kg/kg)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
     epr            ! exner pressure (at main levels)                (1)


REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
!
! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------
!
     gz0             ! roughness length * g of the vertically not
                     ! resolved canopy                               (m2/s2)
!Achtung: Der g-Faktor ist ueberfluessig!

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
     z0_waves        ! roughness length from wave model  (m)

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(OUT) :: &
     !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
     ! velocities 'tvm' and 'tvh' makes live much easier!!

     !turbulent transfer coefficients at the surface
     tcm,          & ! OUT: ... for momentum                         ( -- )
                     ! AUX: specific length-scale fraction           ( -- )
     tch             ! OUT: ... for scalars (heat and moisture)      ( -- )
                     ! AUX: specific length-scale fraction           ( -- )

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
     !turbulent (transfer) velocity scales at the surface
     tvm,          & ! ... for momentum                              ( m/s)
     tvh,          & ! ... for heat and moisture                     ( m/s)

     tfm,          & !
     ! factor for removing the pure drag-contribution of 'tkmmin'     ( --- ) for momentum as INP
     ! Prandtl-layer fraction of total transfer-layer resistance      ( --- ) for momentum as OUT
     tfh,          & !
     ! additional shear-forcing corresponding to the impact of LLDCs  ( 1/s2) at "P"-level as INP
     ! Prandtl-layer fraction of total transfer-layer resistance      ( --- ) for scalars  as OUT
     tfv             ! additional shear-forcing by NTCs               ( 1/s2) at "P"-level as INP

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(INOUT) :: &
     !reciprocal dimensionless diffusion coefficient at top of RL:
     tkr             ! Ustar/(q*Sm)_0                                ( -- )
                     ! INOUT: only, if "imode_trancnf >= 4"
                     ! AUX: also at initialization, if "imode_trancnf == 2 .OR. imode_trancnf == 3"

! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=wp), DIMENSION(nvec,ke1,ntim), TARGET, INTENT(INOUT) :: &
                     ! half-level values of:
     tke             ! q:=SQRT(2*TKE) with TKE='turbul. kin. energy' ( m/s )
                     ! (defined on half levels)
     !Note:
     !'tke' is the "turbulent velocity" (in m/s) and NOT the (mass-density) of turb. kin. energy,
     ! which has the dimension m2/s2!
     !In case of "ntim=1", the actual parameter for 'tke' may be a 2-dim. array for a fix time level.

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
                     ! half-level values of:
     tkvm,         & ! turbulent diffusion coefficient for momentum  (m2/s )
     tkvh            ! turbulent diffusion coefficient for heat      (m2/s )
                     ! (and other scalars)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
     rcld            ! standard deviation of local super-saturation (SDSS)
                     !  at MAIN levels including the lower boundary  (---)
                     ! AUX: cloud-cover at "0"-level (as output of SUB 'adjust_satur_equil'
                     !                                and input of SUB 'solve_turb_budgets')

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
                     ! half-level values of:
     hdef2,        & ! horizontal deformation square at half levels  ( 1/s2 )
     dwdx,         & ! zonal      derivative of vertical wind  ,,    ( 1/s )
     dwdy            ! meridional derivative of vertical wind  ,,    ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
                     ! half-level values of:
     tketens         ! diffusion tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(OUT) :: &
                     ! half-level values of:
     edr             ! eddy dissipation rate of TKE (EDR)            (m2/s3)

REAL (KIND=wp), DIMENSION(:), INTENT(OUT) :: &
!
! Diagnostic near surface variables:
! -----------------------------------------------
!
     t_2m,         & ! temperature in 2m                             (  K  )
     qv_2m,        & ! specific water vapor content in 2m            (kg/kg)
     td_2m,        & ! dew-point in 2m                               (  K  )
     rh_2m,        & ! relative humidity in 2m                       (  %  )
     u_10m,        & ! zonal wind in 10m                             ( m/s )
     v_10m           ! meridional wind in 10m                        ( m/s )

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
!
     shfl_s,       & ! sensible heat flux at the surface             (W/m2)    (positive downward)
     qvfl_s          ! water vapor   flux at the surface             (kg/m2/s) (positive downward)

REAL (KIND=wp), DIMENSION(:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
     umfl_s,       & ! u-momentum flux at the surface                (N/m2)    (positive downward)
     vmfl_s          ! v-momentum flux at the surface                (N/m2)    (positive downward)

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

LOGICAL :: lzacc
INTEGER :: acc_async_queue

INTEGER            :: my_cart_id, my_thrd_id

!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

INTEGER ::      &
    i, k,       & !horizontaler und vertikaler Laufindex
    k1,k2,ks,   & !spezifische Level-Indices
    n,          & !Index fuer diverse Schleifen
!
    nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes (wird zwischen Iterationen zu 'ntur')
    it_durch,   & !Durchgangsindex der Iterationen
    it_start      !Startindex der Iterationen

REAL (KIND=wp) :: &
    fr_tke,             & ! z1/dt_tke
    wert, val1, val2,   & ! Platzhalter fuer beliebige Zwischenergebnisse
    fakt, fac_m, fac_h, & !  ,,         ,,     ,,      Faktoren

!   Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 (fh2,fm2):
    fh2,fm2, &

!   Platzh. fuer horiz. Geschw.(-Komponenten), Geschw.-Quadrat und Druck:
    vel1,vel2, patm, &

!   Platzh. fuer Hoehendifferenzen und Laengaenskalen:
    dh,l_turb,lh,lm,z_surf,len1,len2, &
    fr_sd_h, &     !dimensionless resistance for scalars between the surface and the synoptic "0"-level
    h_2m, h_10m, & !level heights (equal 2m and 10m)
    a_2m, a_10m, & !turbulent distance of 2m- and 10m-level with respect to diag. roughness
    a_atm, &       !turbulent distance of the atmosp. level with respect to land-use roughness
    h_atm, &       !"A"-level heigth of transfer layer (atm. level)
    edgrav,  &     ! 1/grav

!   for TERRA_URB:
    zkbmo_dia, zkbmo_urb, zustar, & ! inverse Stanton number

!   Sonstiges:
    ren_m,ren_h, & !specific Re-numbers at top of land-use roughness
    g_z0_ice, g_len_min, g_alpha1_con_m, & !used for calculating sea-surface roughness
    edprfsecu,   & !"1/prfsecu" (out of ]0; 1[)
    z0wave_threshold, & ! threshold for using z0 provided by wave model
    xf             !scaling factor representing the volume-height ratio of the lowest atmospheric
                   ! full and half level

LOGICAL :: &
    lini,        & !initialization required
    lgz0ini        !initialization of roughness length over water and ice

! Local arrays:

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
!   pointer for variable layers:
    vel1_2d  (:),      &
    vel2_2d  (:),      &
    ta_2d    (:),      &
    qda_2d   (:),      &
!
    g_tet    (:),      &
    g_vap    (:),      &
    qsat_dT  (:),      &
!
    epr_2d   (:),      &
!
    l_tur_z0 (:),      &
!
    z_mom_tot(:),      & ! total roughness length
    a_atm_tot(:),      & ! total turb. dist. (of the adapted virtual atm. level "A")
    a_atm_mod(:),      & ! modified total turb. dist.
    h_atm_mod(:),      & ! modified heigth of the atm. level "A"
    prf_ren_m(:),      & ! log. profile Re-number for momentum
    prf_ren_h(:),      & ! log. profile Re-number for scalars
!
    prss     (:,:),    & ! near-surface pressure (Pa)
    tmps     (:,:),    & ! near-surface temperature-varible (K)
    vaps     (:,:),    & ! near-surface humidity-variable
    liqs     (:,:),    & ! near-surface liquid water content
!
    ediss    (:,:)       ! surface eddy-dissipation rate

REAL (KIND=wp), TARGET    :: &
  ! targets of used pointers
  diss_tar    (nvec,ke1:ke1), & ! eddy dissipation rate (m2/s3)

  ! internal atmospheric variables
  len_scale   (nvec,ke1:ke1), & ! turbulent length-scale (m)
  l_scal      (nvec),         & ! reduced maximal turbulent length scale due to horizontal grid spacing (m)

  fc_min      (nvec),         & ! minimal value for TKE-forcing (1/s2)

  rhon        (nvec,ke1:ke1), & ! boundary level air density (at "0"-level)     (Kg/m3)
  frh         (nvec,ke1:ke1), & ! thermal forcing (1/s2) or thermal acceleration (m/s2)
  frm         (nvec,ke1:ke1), & ! mechan. forcing (1/s2) or mechan. accelaration (m/s2)

  zaux        (nvec,ke1:ke1,ndim), &
                                ! auxilary array containing thermodynamical properties on boundary levels:
                                ! (1:ex_fakt, 2:cp_fakt, 3:dQs/dT, 4:g_tet l, 5:g_vap)
  zvari       (nvec,ke-1:ke1,0:ndim), &
                                ! set of variables used in the turbulent 2-nd order equations
                                ! and later their effective vertical gradients
                                ! 'zvari(:,:,0)' is reserved for half-level pressure or "circulation acceleration"
                                ! (finally, 'zvari(:,ke1,tet|vap)' also contains T2m and qv_2m)

  rcls        (nvec,ke1),     & ! double for standard deviation of the local super-saturation (SDSS)
                                ! ('rcls(:,ke1)' provides 2m-SDSS to, and receives 2m-cl_cov from, SUB 'adjust_satur_equil')

  tl_s_2d     (nvec),         & ! surface-level value of conserved temperature
                                !   (liquid water temperature) (K)
  qt_s_2d     (nvec),         & ! surface-level value of conserved humidity
                                !   (total water)
  vel_2d      (nvec,ke:ke),   & ! wind speed (m/s) at the lowest full model level

  velmin      (nvec),         & ! modified 'vel_min' used for tuning corrections (m/s)
                                ! (hyper-parameterizations)

  ! internal variables for the resistance model
  hk_2d       (nvec),         & ! mid/full-level height above ground belonging to 'k_2d' (m)
  hk1_2d      (nvec),         & ! mid/full-level height above ground of the previous layer (below) (m)
  hk2_2d      (nvec),         & ! 'hk1_2d, hk1_2d and hk2_2d' are also used a shelf for any level height

  a_atm_2d    (nvec),         & ! turbulent distance of the lowermost full model level

  h_top_2d    (nvec),         & ! boundary-level height of transfer layer (top  level "P")
  h_atm_2d    (nvec),         & ! mid     -level heigth of transfer layer (atm. level "A")
  h_can_2d    (nvec),         & ! effective canopy height (m) used as depth of the R-layer, applied
                                !  for calculation of the roughness layer resistance of momentum or
                                !  for determination of the synoptic 2m-level in case of a diagnostic
                                !  exponential R-layer profile (at "itype_2m_diag == 2")

  edh         (nvec),         & ! reciprocal of any layer depth
  z0m_2d      (nvec),         & ! mean  roughness length
  z0d_2d      (nvec),         & ! diag. roughness length (for the SYNOP lawn)
  z2m_2d      (nvec),         & ! height of 2m  level (above the surface) or total roughn. length for momentum
  z10m_2d     (nvec),         & ! height of 10m level (above the surface)

  rat_m_2d    (nvec),         & ! any surface layer ratio for momentum
                                !   (like Re-number or Stability factor)
  rat_h_2d    (nvec),         & ! any surface layer ratio for scalars
                                !   (like Re-number or Stability factor)
  fac_h_2d    (nvec),         & ! surface layer profile factor for scalars
  fac_m_2d    (nvec),         & ! surface layer profile factor for momentum

  frc_2d      (nvec),         & ! saved amplfification factor for wind-shear forcing

  val_m       (nvec),         & ! any specific value for momentum (Re-numer or resistance length)
  val_h       (nvec),         & ! any specific value for scalars  (Re-numer or resistance length)

  dz_sa_h     (nvec),         & ! total transfer resistance length for scalars (m)
  dz_0a_h     (nvec),         & ! turbulent Prandtl-layer resistance length for scalars (m)
  dz_s0_h     (nvec),         & ! total roughness layer resistance length for scalars (m)
  dz_g0_h     (nvec),         & ! turbulent roughness layer resistance length for scalars (m)
  dz_sg_h     (nvec),         & ! laminar resistance length for scalars (m)

  dz_s0_m     (nvec),         & ! total roughness layer resistance length for momentum (m)
  dz_0a_m     (nvec)            ! turbulent Prandtl-layer resistance length for momentum (m)

REAL (KIND=wp) ::             &
  dz_sa_m,                    & ! total transfer resistance length for momentum (m)
  dz_sg_m                       ! laminar resistance length for momentum (m)

REAL (KIND=wp) ::             &
  grad        (nvec,nmvar)      ! any vertical gradient

INTEGER        ::             &
  k_2d        (nvec)            ! index field of the upper level index to be used
                                !   for near surface diagn.
LOGICAL        ::   ldebug = .FALSE.

!---- End of header -----------------------------------------------------------

  CALL set_acc_host_or_device(lzacc, lacc)

  IF(PRESENT(opt_acc_async_queue)) THEN
    acc_async_queue = opt_acc_async_queue
  ELSE
    acc_async_queue = 1
  ENDIF

!==============================================================================
! Begin subroutine turbtran
!------------------------------------------------------------------------------

! 1)  Vorbereitungen:

  !GPU data region of all local variables except pointers which are set later on
  !$ACC DATA &
! local variables
  !$ACC   CREATE(diss_tar, len_scale, l_scal, fc_min) &
  !$ACC   CREATE(rhon, frh, frm, zaux, zvari, rcls) &
  !$ACC   CREATE(tl_s_2d, qt_s_2d, vel_2d, velmin) &
  !$ACC   CREATE(hk_2d, hk1_2d, hk2_2d, h_top_2d, h_atm_2d, a_atm_2d, h_can_2d) &
  !$ACC   CREATE(edh, val_m, val_h, z0m_2d, z0d_2d, z2m_2d) &
  !$ACC   CREATE(z10m_2d, rat_m_2d, rat_h_2d, fac_h_2d, fac_m_2d) &
  !$ACC   CREATE(frc_2d, dz_s0_m, dz_sg_h, dz_g0_h) &
  !$ACC   CREATE(dz_0a_m, dz_0a_h, dz_sa_h, dz_s0_h, grad) &
  !$ACC   CREATE(k_2d) &
  !$ACC   COPYIN(ivend) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)


  ! Pointer assignments:

  IF (PRESENT(edr)) THEN
     ediss => edr
  ELSE
     ediss => diss_tar
  END IF

  prss(1:,ke-1:) => zvari(:,:,0)    ! near-surface pressure (Pa)
  tmps(1:,ke-1:) => zvari(:,:,tet)  ! near-surface temperature-variable (K)
  vaps(1:,ke-1:) => zvari(:,:,vap)  ! near-surface humidity-variable
  liqs(1:,ke-1:) => zvari(:,:,liq)  ! near-surface liquid water content

  vel1_2d => zvari(:,ke ,u_m)
  vel2_2d => zvari(:,ke ,v_m)

  ta_2d   => zvari(:,ke ,tet_l)
  qda_2d  => zvari(:,ke ,h2o_g)

  epr_2d  => zaux(:,ke1,1)
  qsat_dT => zaux(:,ke1,3)
  g_tet   => zaux(:,ke1,4)
  g_vap   => zaux(:,ke1,5)

  l_tur_z0 => len_scale(:,ke1)

  prf_ren_h => val_h
  z_mom_tot => z0m_2d   !total roughness length                     equals local roughness length
  a_atm_tot => a_atm_2d !total turbulent distance                   equals local turbulent dist.
  a_atm_mod => a_atm_2d !modified total turbulent distance          equals local turbulent dist.
  h_atm_mod => h_atm_2d !modified heigth of the atmospheric level   equals local one
  prf_ren_m => val_h    !logarithmic profile Re-number for momentum equals that one for scalars

  !SUB-arguments and pointers explicitly operating in this code:
  !$ACC DATA PRESENT(fr_land, t_g, l_lake, l_sice) &
  !$ACC   PRESENT(prss, tmps, vaps, liqs) &
  !$ACC   PRESENT(z_mom_tot, a_atm_tot, a_atm_mod, h_atm_mod, prf_ren_m, prf_ren_h) &
  !$ACC   PRESENT(hhl, epr_2d, u, v, t) &
  !$ACC   PRESENT(ediss, l_pat) &
  !$ACC   PRESENT(epr, tke, tkvm, tkvh, gz0, z0_waves, tkr) &
  !$ACC   PRESENT(l_tur_z0, ps, sai, urb_isa, rlamh_fac) &
  !$ACC   PRESENT(tfm, tfh, tfv, qv_s, qv, qc) &
  !$ACC   PRESENT(dwdx, dwdy, hdef2, g_tet) &
  !$ACC   PRESENT(g_vap, tcm, tch, z0d_2d, shfl_s) &
  !$ACC   PRESENT(rcld, qsat_dT, qvfl_s, umfl_s, vmfl_s) &
  !$ACC   PRESENT(edr, tvh, tvm, t_2m) &
  !$ACC   PRESENT(qv_2m, qda_2d, ta_2d) &
  !$ACC   PRESENT(v_10m, u_10m, vel1_2d, vel2_2d, rh_2m, td_2m) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)

!-------------------------------------------------------------------------------
  CALL turb_setup (tdc=tdc, i_st=ivstart, i_en=ivend, k_st=ke, k_en=ke1, &
                   iini=iini, dt_tke=dt_tke, nprv=nprv, l_hori=l_hori, &
                   ps=ps, t_g=t_g, qv_s=qv_s, qc_a=qc(:,ke), &
                   lini=lini, it_start=it_start, nvor=nvor, fr_tke=fr_tke, &
                   l_scal=l_scal, fc_min=fc_min, &
                   prss=prss(:,ke1), tmps=tmps(:,ke1), vaps=vaps(:,ke1), liqs=liqs(:,ke1), rcld=rcld, &
                   lacc=lzacc, opt_acc_async_queue=acc_async_queue)

#ifdef ICON_USE_CUDA_GRAPH
   IF (lzacc .AND. lini .AND. lcuda_graph_turb_tran ) THEN
      CALL finish ('turbtran', 'initialization is not supported when capturing a graph with OpenACC')
   END IF
#endif
!-------------------------------------------------------------------------------

  my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
  my_thrd_id = omp_get_thread_num()
#endif

! 2)  Initialisierung der z0-Werte ueber Meer
!     und der laminaren Transferfaktoren:

      ! Berechnung einiger Hilfsgroessen und Initialisierung der Diffusionskoeff.:

      ! Unterste Hauptflaeche halbiert den Abstand zur untersten Nebenflaeche:
      xf=z2
      xf=z1/xf

      ! Hoehe des 2m- und 10m-Niveaus:
      h_2m  = z2
      h_10m = z10

      edprfsecu = z1/tdc%prfsecu !used for limitation of the profile funtion
      edgrav = z1/grav

      ! Fixed parameter used for calculating sea-surface roughness:
      g_z0_ice=grav*tdc%z0_ice; g_len_min=grav*tdc%len_min; g_alpha1_con_m=grav*tdc%alpha1*con_m

      ! Provisional tuning for minimum roughness length used from wave model
      z0wave_threshold = 1.e-4_wp

      ks=ke

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         ! Dicke der Modell-Prandtl-Schicht
         h_top_2d(i) = hhl(i,ke)-hhl(i,ke1)
         h_atm_2d(i) = h_top_2d(i)*xf

         ! Surface-Exner-pressure:
         epr_2d(i) = zexner(ps(i))

      END DO

      !$ACC LOOP SEQ
      DO k=ks, ke
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            vel_2d(i,k) = MAX( tdc%vel_min, SQRT( u(i,k)**2+v(i,k)**2 ) ) !wind speed
         END DO
      END DO

!<Tuning
!---------------------------------------------------------------------------------------
      IF (tdc%imode_vel_min == 2) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            ! stability-dependent minimum velocity serving as lower limit on surface TKE
            ! (parameterizes small-scale circulations developing over a strongly heated surface;
            ! tuned to get 0.75 m/s when the land surface is at least 7.5 K warmer than the air in the
            ! lowest model level; nothing is set over water because this turned out to induce
            ! detrimental effects in NH winter)

            velmin(i) = MAX( tdc%vel_min, MIN(0.75_wp, fr_land(i)*(t_g(i)/epr_2d(i) - t(i,ke)/epr(i,ke))/ &
                        LOG(2.e3_wp*h_atm_2d(i))) )
         END DO
      END IF
!---------------------------------------------------------------------------------------
!>Tuning: This kind of correction should be substituded by a less ad-hoc approach.

      IF (lini) THEN !only for initialization

         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(lgz0ini, l_turb, dh, vel1, vel2, fm2, fh2, fakt, lm, lh, wert, val1, val2)
         DO i=ivstart, ivend

            lgz0ini=(igz0inp == 0 .AND. fr_land(i) <= z1d2)
            !Note: This definition of a non-land surface is now in line with the ICON-definition
            !       using "frlnd_thrhld=z1d2"

            !Einfachste Schaetzung der Schubspannung als Impusls-Flussdichte durch die Nebenflaeche "P" (k=ke)
            ! mit Hilfe einer diagnostischen TKE ohne Beruecksichtigung von Feuchte-Effekten und unter Anwendung
            ! von "Rf=Ri" bei der Berechnung der Stabilitaetsfunktionen:

            l_turb=h_top_2d(i) !approx. turb. length scale at "P"-level (k=ke)

            l_turb=tdc%akt*MAX( tdc%len_min, l_turb/( z1+l_turb/l_scal(i) ) )

            dh=z1d2*(hhl(i,ke-1)-hhl(i,ke1))

            vel1=u(i,ke-1)
            vel2=u(i,ke  )
            grad(i,u_m)=(vel1-vel2)/dh

            vel1=v(i,ke-1)
            vel2=v(i,ke  )
            grad(i,v_m)=(vel1-vel2)/dh

            grad(i,tet_l)=(t(i,ke-1)-t(i,ke))/dh + tet_g

            fm2=MAX( grad(i,u_m)**2+grad(i,v_m)**2, fc_min(i) )
            fh2=grav*grad(i,tet_l)/t(i,ke)

            ! Vereinfachte Loesung mit Rf=Ri:
            IF (fh2 >= (z1-rim)*fm2) THEN !die krit. Ri-Zahl wird ueberschritten
               !'lm' sowie 'lh' werden durch 'lm' bei der krit. Ri-Zahl angenaehert:
               fakt=z1/rim-z1
               lm=l_turb*(sm_0-(a_6+a_3)*fakt)
               lh=lm
            ELSE
               fakt=fh2/(fm2-fh2)
               lm=l_turb*(sm_0-(a_6+a_3)*fakt)
               lh=l_turb*(sh_0-a_5*fakt)
            END IF

            val1=lm*fm2; val2=lh*fh2
            wert=MAX( val1-val2, rim*val1 )

            IF (ltkeinp) THEN
               tke(i,ke,nvor)=tke(i,ke,ntur)
            ELSE
               tke(i,ke,nvor)=SQRT(tdc%d_mom*l_turb*wert)
            END IF

            val1=con_m; tkvm(i,ke)=lm*tke(i,ke,nvor)
            val2=con_h; tkvh(i,ke)=lh*tke(i,ke,nvor)

            !Note:
            !'tk[h|m]min' are, first of all, foreseen as lower limits for 'vertdiff'-calculations;
            ! hence, they are not required for this particular initialization of "P"-level DCs.
            !Nevertheless, positive-definite initial 'tkv[m|h](:,ke)' are required for initialization
            ! of both, friction-velocity (and with it 'gz0') at sea-surfaces, and 'tkv[m|h](:,ke1)'.
            !While, at "imode_trancnf<4", these initial DC-values at "P"-and "0"-level determine the
            ! surface-layer profile-functions, the one at "0"-level needs also to be present as input
            ! of the Turbulence Model (TMod) in SUB 'solve_turb_budgets'.
            !For that purposes, at least the laminar limit appears to by reasonable at "P"-level.
            !'tke(:,ke,nvor)' is only required for the final 'gz0'-calculation for sea-surfaces
            ! (at the end of SUB 'turbtran').

            IF (tdc%imode_tkemini >= 2) THEN !adaptation of TKE and TMod. to lower limits
               tke(i,ke,nvor)=tke(i,ke,nvor)*MAX( z1, val2/tkvh(i,ke) ) !adapted 'tke'
               !Note:
               !Compare the respective 'tke'-adaptation in SUB 'turbdiff'.
            END IF

            tkvh(i,ke)=MAX(val2, tkvh(i,ke)) !'tkvh(:,ke)' with lower limit
            tkvm(i,ke)=MAX(val1, tkvm(i,ke)) !'tkvm(:,ke)' with full lower limit

            !Note:
            !Laminar LLDCs do not contain any artificial drag contribution. As, moreover, 'tkv[m|h](:,ke1)'
            ! are being estimated by means of 'tkv[m|h](:,ke)', 'imode_suradap' needs not to be considered
            ! within this initialization at all.

            vel2=MAX( tdc%epsi, tkvm(i,ke)*SQRT(fm2) ) !estimated Ustar**2
            vel1=SQRT(vel2) !Ustar

            IF (lgz0ini) THEN !initialization of roughness length for water- or ice-covered surface:
               IF ( l_sice(i) ) THEN !ice-covered surface
                  gz0(i)=g_z0_ice
               ELSE !water-covered surface
                  ! Basic Charnock-parameter; use enhanced value of "0.1" over lakes if "imode_charpar>1":
                  fakt=MERGE( tdc%alpha0, 0.1_wp, tdc%imode_charpar == 1 )
                  ! use velocity-dependent Charnock paramter over sea if "imode_charpar>1":
                  fakt=MERGE( alpha0_char( vel_2d(i,ke), tdc ), fakt, tdc%imode_charpar > 1 .AND. .NOT.l_lake(i) )

                  !Final diagnosed 'gz0' (shear-related dynamic contribution and laminar correction):
                  gz0(i)=MAX( g_len_min, fakt*vel2 + g_alpha1_con_m/vel1 )
               END IF
            END IF

            tkr(i)=l_turb*vel1                          !l_0*Ustar

            rat_m_2d(i)= tkr(i)/tkvm(i,ke)              !Ustar/(q*Sm)_P
            rat_h_2d(i)=(tkr(i)*sh_0)/(tkvh(i,ke)*sm_0) !Ustar/(q*Sh)_P*Sh(0)/Sm(0)

         END DO

      END IF !only for initialization

!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         z0m_2d(i)  = gz0(i)*edgrav !mean roughness-length
         l_tur_z0(i)= tdc%akt*z0m_2d(i) !turbulent length scale
         tcm(i) = z0m_2d(i)/(h_top_2d(i)+z0m_2d(i)) !unspecific length scale fraction
         tch(i) = tcm(i)       !default setting with unspecific length scale fraction

         a_atm_2d(i)=h_atm_2d(i)+z0m_2d(i) !turbulent distance of the lowermost full model level "A"

         !Note: The additional 'grav'-factor in 'gz0' is obsolete and should be removed!
      END DO

      IF (lini) THEN !only for initialization

         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            tkvm(i,ke1)=tkvm(i,ke)*tcm(i)*(tcm(i)+(z1-tcm(i))*rat_m_2d(i))
            tkvh(i,ke1)=tkvh(i,ke)*tch(i)*(tch(i)+(z1-tch(i))*rat_h_2d(i))

            tkr(i)=tkr(i)/tkvm(i,ke1) !Ustar/(q*Sm)_0
         END DO

      ELSEIF (ladsshr .OR. (tdc%imode_trancnf < 4 .AND. tdc%imode_suradap >= 1)) THEN
             !for surface-layer adaptations to additional "P"-level shear (only apart from initialization)

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt, wert, val1, val2)
         DO i=ivstart, ivend
            fakt=(tfm(i)*tkvm(i,ke))/tkvm(i,ke1)**2 !amplification factor for transmitted shear-forcing
            wert=MERGE( z1/ediss(i,ke1), tdc%d_mom*l_tur_z0(i)/tke(i,ke1,nvor)**3, PRESENT(edr) ) !reciprocal EDR
            wert=fakt*wert*tkvm(i,ke1) !related squared time-scale

            val1=MERGE( tfh(i), z0, (tdc%imode_tkemini > 1) ) + tfv(i) !total addit. shear-forc. potentially transmitted from level "P"
            val2=MERGE( z0, tfh(i), (tdc%imode_tkemini > 1) ) + MERGE( z1, z1-tdc%rsur_sher, iini <= 0 )*val1 !not yet transmitted
                 !total additional shear at level "P", which is "tfh+tfv" at the first time step.

            tfv(i)=z2/( z1+SQRT( z1+z4*wert*val2 ) ) !according reduction factor

            IF (ladsshr) frc_2d(i)=tdc%rsur_sher*fakt*tfv(i)*val1 !additional shear-forcing directly transmitted from level "P"

            !Note:
            !It is assumed that the momentum-flux related to the transmitted part of additional "P"-level shear
            ! is constant wihtin the transfer-layer between levels "0" and "P".
            !For estimation of the transmitted contribution of additional shear by NTCs and LLDCs at level "P" to level "0",
            ! the used "P"-level DC 'tkvm(:,ke)' is relieved from artificial drag-related contributions by means of 'tfm'.
            !'tfv' expresses a correction applied to the "0"-level DC 'tkvm(:,ke1)', which is due to additional "P-level shear
            ! being not transmitted to "0-level.
            !At the first time-step after initialization, the "0"-level DC has not yet been affected by any transmitted
            ! additional shear. Thus, in this case, the total additional shear at level "P" has not yet been transmitted!
         END DO

         IF (tdc%imode_trancnf < 4 .AND. tdc%imode_suradap >= 1) THEN
            !profile-factors (being calculated by means of an upper node) require LLDC-corrention:
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
            DO i=ivstart, ivend
               fakt=MERGE( SQRT( tfv(i) ), z1, (tdc%imode_suradap > 1) ) !reduction-factor due to not-transmitted
                                                                         ! additonal shear from level "P"

               !length-scale fractions including required reduction factors:
               tcm(i)=tcm(i)*fakt*tfm(i)
               tch(i)=tch(i)*fakt

               !Note:
               !These reductions avoid a distortion of the profile-function due to not consistent DC-profiles.
            END DO
         END IF
      END IF !for surface-layer adaptations to additional "P"-level shear (only apart from initialization)

      IF (tdc%imode_trancnf >= 4) THEN
         !Calculation the profile-factors without an upper node of diffusion coefficients (above the RL)
         !but based on previous values of Ustar and diffusion-coefficients (at the top of the RL):
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            rat_m_2d(i)=tkr(i)                                       !Ustar/(q*Sm)_0
            rat_h_2d(i)=tkr(i)*(tkvm(i,ke1)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(Sh(0)/Sm(0))
         END DO
      ELSE
         !Profile-factors by using the previous diffusion coefficients
         !without a laminar correction, but still based on the upper node:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            rat_m_2d(i)=tcm(i)*tkvm(i,ke)/tkvm(i,ke1) !(q*Sm)_P/(q*Sm)_0
            rat_h_2d(i)=tch(i)*tkvh(i,ke)/tkvh(i,ke1) !(q*Sh)_P/(q*Sh)_0
         END DO
      END IF

      !$ACC END PARALLEL

! 4)  Berechnung der Transfer-Geschwindigkeiten:

!----------------------------------------------
      DO it_durch=it_start, tdc%it_end !Iterationen
!----------------------------------------------

         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(z_surf, val1, val2, ren_m, ren_h, dz_sg_m, wert)
         DO i=ivstart, ivend
            z_surf= z0m_2d(i)/sai(i) !effective R-length
            val1=tdc%rlam_heat    & !Maximal DimensionLess Laminar Resistance-Complement (MDLLRC) for a rigid suface
                *rlamh_fac(i)     & !related tuning factor (depending on skin conductivity and analyzed T2M|RH2M bias)
                *(z1+(z1-REAL(NINT(fr_land(i)),wp))*(tdc%rat_sea-z1)) & !correction-factor for a (liquid) water surface
                *MERGE( tdc%rat_glac, 1._wp, gz0(i)<0.01_wp .AND. fr_land(i)>z1d2 ) !and for a glacier surface
            !Notes:
            !The applied identific. of a land surface is now in line with the ICON-definition using "frlnd_thrhld=z1d2"
            !The standard-corrections of MDLLRC may be due to possible motions of the surface by the action of wind drag.

            ! Effective resistance-length values of local R-layer for scalars:
            IF (tdc%imode_lamdiff == 1) THEN !laminar DC-limitation only
               ! Based on a laminar limit of turbulent Diffusion Coefficients (CD):

               tkvm(i,ke1)=MAX( con_m, tkvm(i,ke1) )
               tkvh(i,ke1)=MAX( con_h, tkvh(i,ke1) )

               ren_m=tkvm(i,ke1)/con_m !Re-number for momentum
               ren_h=tkvh(i,ke1)/con_h !Re-number for scalars

               dz_sg_h(i)=z_surf*val1*(ren_h/ren_m) !resistance-length related to MDLLRC for scalars

               dz_g0_h(i)=z_surf*LOG(ren_m) !resistance-length related to general dim.less resistance
                                            ! of R-layer due to pure turbulence (without laminar restriction)
                                            ! at neutral stratification but limited by the pure laminar resistance
            ELSE !permanent inclusion of laminar transport
               ! Based on laminar correction of R-layer resistances (by means of current DCs):

               ren_h=tkvh(i,ke1)/con_h !Re-number for scalars
               val2=LOG(ren_h+z1) !dim.less resistance of R-layer for scalars due to pure turbulence
                                  ! (without laminar reduction) but with the effect of parallel laminar transport

               dz_sg_h(i)=z_surf*MAX( z0, MIN( ren_h-val2, val1 ) ) !resistance-length related to reduced MDLLRC
                                                                    ! so as to keep the laminar-resistance limit
               dz_g0_h(i)=z_surf*val2 !resistance-length related to non-restricted turbulence at neutral stratif.

               !Note:
               !The used constant MDLLRC is only applicable in case of large Re-numbers at the "0"-level (k=ke1).
               !Due to the laminar resistance limit, "val1+val2<=ren_h" needs always to be fulfilled.
               !In the laminar limit, the dim.-less turb.-lam. R-layer resistance 'val2' requires no further DLLRC,
               ! and thus, the real LRC [in s/m] (and even more the related real DLLRC 'val1') has to vanish.
               !This form approximates a more-complicate general solution of the resistance integral for turbulent-laminar
               ! transport throughout an internal BL, being modified by R-elements and with arbitrarily stratification.
            END IF

            dz_s0_h(i)=dz_sg_h(i)+dz_g0_h(i) !through full R-layer

            ! Effective height of the R-layer:
            IF (tdc%rlam_mom > z0 .OR. tdc%itype_2m_diag == 2 ) THEN !R-height required
               h_can_2d(i)=tdc%rat_can*MERGE( sai(i)*z0m_2d(i), &
                                              dz_s0_h(i)*LOG(dz_s0_h(i)/dz_sg_h(i)), &
                                              dz_sg_h(i) == z0 )
            END IF

            ! Effective resistance-length values of local R-layer for momentum:
            IF (tdc%rlam_mom > z0) THEN !including R-layer resistance for momentum
               dz_sg_m=tdc%rlam_mom*z_surf !through laminar layer
               wert=z1d2*dz_sg_m
               dz_s0_m(i)=wert+SQRT(wert**2+h_can_2d(i)*dz_sg_m) !through full R-layer
            ELSE
               dz_s0_m(i)=z0 !no R-layer resistance for momentum
            END IF
         END DO

!--------------------------------------------------------------------------
         IF (lterra_urb .AND. (.NOT. itype_kbmo == 1)) THEN
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zkbmo_dia, zustar, zkbmo_urb)
           DO i=ivstart, ivend
              IF (urb_isa(i) > 0.0_wp) THEN
                zkbmo_dia    = dz_s0_h(i)/z0m_2d(i)
                zustar       = SQRT(tvm(i)*vel_2d(i,ke))
                IF (itype_kbmo == 2) THEN
                  ! Brutsaert Kanda parameterisation for bluff-body elements
                  ! for some reason, it doesn't withstand zero values of ustar
                  zkbmo_urb = MAX(0.1_wp, 1.29_wp * (z0m_2d(i)*zustar/con_m)**0.25_wp - z2)
                ELSE !(itype_kbmo == 3)
                  ! Zilitinkevich
                  zkbmo_urb = MAX(0.1_wp, 0.13_wp * (z0m_2d(i)*zustar/con_m)**0.45_wp)
                END IF

                dz_s0_h(i) = (urb_isa(i)*zkbmo_urb + (1.0_wp - urb_isa(i))*zkbmo_dia) * z0m_2d(i)
              ENDIF
           END DO
         END IF
!--------------------------------------------------------------------------

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
         DO i=ivstart, ivend

!           Profilfakoren der turbulenten Prandtl-Schicht:

            fakt=z0m_2d(i)/h_top_2d(i)

            IF (tdc%imode_trancnf < 4) THEN
               !Profile-factors by employing previous values of the diffusion-coefficients
               !at the top fo the roughness-layer (0) and also at the upper bound of the
               !lowest atm. model layer (P) as an upper node:

               rat_m_2d(i)=MIN( edprfsecu, MAX( tdc%prfsecu, rat_m_2d(i) ) ) !limited (q*Sm)_P/(q*Sm)_0
               rat_h_2d(i)=MIN( edprfsecu, MAX( tdc%prfsecu, rat_h_2d(i) ) ) !limited (q*Sh)_P/(q*Sh)_0

               !Note:
               !This quite arbitrary relative limitation of 'rat_m|h' avoids not realistic outlayers
               ! (mainly) due to numerical problems that may arise through the applied time-step iteration.
               !It is "0 < prfsecu < 1"; and 'prfsecu' should be chosen as small as possible in order
               ! to avoid artificial degenerations of the below solution of the P-layer resistance integral.
               !'rat_[m|h]' should be based on the ratio of pure turbulent DCs for the levels ("P" and "0"),
               ! thus, any (e.g. laminar) Lower Limitation of DCs (LLDC) artificially degenerates simulated
               ! stratification, unless it is implemented in accordance with the general BL-approximation.

               fac_m_2d(i)=(rat_m_2d(i)-z1)*fakt !non-stab. profile-factor for momentum
               fac_h_2d(i)=(rat_h_2d(i)-z1)*fakt !non-stab. profile-factor for scalars

            ELSE !Profile-factors without using the upper node

               fac_m_2d(i)=z1-rat_m_2d(i) !profile-factor for momentum
               fac_h_2d(i)=z1-rat_h_2d(i) !profile-factor for scalars

               !Note:
               !In this case, no restriction of 'fac_m|h' is necessary in order to avoid singularities
               ! in the below formula for the resistance integral.
            END IF

         END DO

!        Effektive Widerstandslaengen der turb. Prandtl-Schicht:

        ! Preparations with regard to the hyperbolic interpolation of profile-function:

         IF (tdc%imode_trancnf >= 3) THEN !hyperbolic interpolation of profile-function for stable strat.
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               prf_ren_h(i)=LOG(a_atm_2d(i)/z0m_2d(i)) !log. profile Re-number for scalars
            END DO
         END IF

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(a_atm, z_surf, h_atm, fac_h, fac_m)
         DO i=ivstart, ivend
            a_atm=a_atm_2d(i)
            z_surf=z0m_2d(i)
            h_atm=h_atm_2d(i)

            fac_m=fac_m_2d(i)

            IF (fac_m >= z0 .OR. tdc%imode_trancnf < 3) THEN
               !non-stable stratfication or based on linear interpolation of profile-function
               !for the velocity scale (q*Sm):
               dz_0a_m(i)=z_surf*MERGE( (z_surf*h_atm)/(a_atm_tot(i)*z_mom_tot(i)), &
                                        LOG( a_atm_mod(i)/(z_surf+fac_m*h_atm_mod(i)) )/(z1-fac_m), &
                                        fac_m == z1 )
            ELSEIF (tdc%imode_trancnf == 3) THEN
               !based on hyperbolic interpolation of profile-function (q*Sm) for stable stratification,
               ! at which the upper node for diffusion-coefficients is used:
               fac_m_2d(i)=-fac_m/rat_m_2d(i) !transformed profile-factor for stable stratification
               dz_0a_m(i)=z_surf*(z1-fac_m)*prf_ren_m(i)+fac_m*h_atm
            ELSE !(tdc%imode_trancnf == 4): without using the upper node for diff.-coefs.
               dz_0a_m(i)=(z_surf*prf_ren_m(i)-fac_m*h_atm)/(z1-fac_m)
            END IF

            fac_h=fac_h_2d(i)

            IF (fac_h >= z0 .OR. tdc%imode_trancnf < 3) THEN
               dz_0a_h(i)=z_surf*MERGE( h_atm/a_atm, &
                                        LOG( a_atm/(z_surf+fac_h*h_atm) )/(z1-fac_h), &
                                        fac_h == z1 )
            ELSEIF (tdc%imode_trancnf == 3) THEN !the upper node for diffusion-coefficients is used
               fac_h_2d(i)=-fac_h/rat_h_2d(i) !transformed profile-factor for stable stratification
               dz_0a_h(i)=z_surf*(z1-fac_h)*prf_ren_h(i)+fac_h*h_atm
            ELSE !(tdc%imode_trancnf == 4): without using the upper node for diff.-coefs.
               dz_0a_h(i)=(z_surf*prf_ren_h(i)-fac_h*h_atm)/(z1-fac_h)
            END IF

            !Note:
            !This solution corresponds to a quite general solution of the resistance integral for pure
            ! turbulent transport throughout an arbitrarily stratified internal BL above the R-layer
            ! at high Re-numbers (that means at negligible laminar effects).
         END DO

         IF (tdc%imode_lamdiff == 2) THEN !permanent inclusion of laminar transport
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               !Securing at least laminar transport based on current DCs:
               dz_0a_m(i) = MIN( h_atm_2d(i)*tkvm(i,ke1)/con_m, dz_0a_m(i) )
               dz_0a_h(i) = MIN( h_atm_2d(i)*tkvh(i,ke1)/con_h, dz_0a_h(i) )
            END DO
         END IF

         ! Combination of resistance-length values for momentum and scalars:

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(dz_sa_m)
         DO i=ivstart, ivend

!           Effektive Widerstandslaengen von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):

            ! total resistance-path of the transfer-layer for momentum:
            dz_sa_m    = dz_0a_m(i) + dz_s0_m(i)

            ! total resistance-path of the transfer-layer for scalars:
            dz_sa_h(i) = dz_s0_h(i) + dz_0a_h(i)

!           Reduktionsfaktoren fuer die Bestandesschicht incl. lam. Grenzschicht:

            tfm(i)=dz_0a_m(i)/dz_sa_m    !for momentum
            tfh(i)=dz_0a_h(i)/dz_sa_h(i) !for scalars

!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den Faktor 'rat_lam'
!           gegenueber fuehlbarer Waerme vergroesserten laminaren Transportwiderstandes:

            tfv(i)=z1/(z1+(tdc%rat_lam-z1)*dz_sg_h(i)/dz_sa_h(i))

            !Note:
            !So far, this reduction factor for evaporation is only applied in 'terra',
            ! hence, it is not affecting evaporation of the open sea or sea-ice!
         END DO

!        Berechnung der Erhaltungsgroessen in der Prandtl-Schicht:

         IF (tdc%icldm_tran == -1 .OR. tdc%ilow_def_cond == 2) THEN
            !conserved values at the rigid surface are temperature and humidity
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tl_s_2d(i)=t_g(i); qt_s_2d(i)=qv_s(i)
            END DO
         ELSE !conserved variables at the rigid surface depend on liquid water
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tl_s_2d(i)= t_g(i) - lhocp*liqs(i,ke1)
               qt_s_2d(i)=qv_s(i) +       liqs(i,ke1)
            END DO
         END IF

         !$ACC LOOP SEQ
         DO k=ks, ke
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               zvari(i,k,u_m)=u(i,k)
               zvari(i,k,v_m)=v(i,k)
            END DO
         END DO

         IF (tdc%icldm_tran == -1) THEN !no water phase change possible
            !$ACC LOOP SEQ
            DO k=ks, ke
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=ivstart, ivend
                  zvari(i,k,tet_l)=t(i,k)/epr(i,k)
                  zvari(i,k,h2o_g)=qv(i,k)
               END DO
            END DO
         ELSE !water phase changes are possible
            !$ACC LOOP SEQ
            DO k=ks, ke
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=ivstart, ivend
                  zvari(i,k,tet_l)=(t(i,k) - lhocp*qc(i,k))/epr(i,k)
                  zvari(i,k,h2o_g)=qv(i,k) +       qc(i,k)
               END DO
            END DO
         END IF

!        Thermodynamische Hilfsvariablen auf dem Unterrand der Prandtl-Schicht:

         !$ACC END PARALLEL

         CALL adjust_satur_equil( tdc=tdc, i1dim=nvec, khi=ke1, ktp=ke-1, & !in
!
              i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1,            & !in
!
              lcalrho=.TRUE.,  lcalepr=.FALSE., lcaltdv=.TRUE.,        & !in
              lpotinp=.FALSE., ladjout=.FALSE.,                        & !in
!
              icldmod=tdc%icldm_tran,                                  & !in
!
              zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          & !in
!
              prs=prss, t=tmps, qv=vaps, qc=liqs,                      & !in (surface values at "0"-level 'ke1')
!
              psf=ps, fip=tfh,                                         & !in
!
              rcld=rcld,  & !inp: std. deviat. of local super-saturat.
                            !out: saturation fraction (cloud-cover)
!
              dens=rhon,         exner=zaux(:,:,1),                    & !out
              r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3),                    & !out
              g_tet=zaux(:,:,4), g_h2o=zaux(:,:,5),                    & !out
!
              tet_liq=zvari(:,:,tet_l), q_h2o=zvari(:,:,h2o_g),        & !inout (inp as target of 'tmps, vaps')
                                        q_liq=zvari(:,:,liq),          & !out
!
              lacc=lzacc, opt_acc_async_queue=acc_async_queue )

         !Beachte:
         !'zvari(:,ke1,tet_l)' und 'zvari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen am Unterrand der
         ! Prandtl-Schicht, waehrend  ta_2d' => 'zvari(:,ke,tet_l) und 'qda_2d  => 'zvari(:,ke,h2o_g)' auf diese
         ! Groessen bzgl. der untersten Hauptflaeche zeigen, welche zur Interpolation der Groessen an der Oberflaeche
         ! auf jenen Unterrand der Prandtl-Schicht (Oberrand der Rauhigkeitsschicht) benutzt werden.

         ! Berechnung der benoetigten Vertikalgradienten und der TKE-Antriebe:

         !Vertikalgradienten des Horizontalwindes:

         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            val_m(i)=tfm(i)/dz_0a_m(i) !reciprocal resistance length of roughness- and laminar-layer
         END DO

         !$ACC LOOP SEQ
         DO n=1, nvel
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               zvari(i,ke1,n)=zvari(i,ke,n)*val_m(i) !vertical gradient of wind-component
            END DO
            !Beachte: Dies ist die Darstellung ohne Nutzung der unteren Randwerte der Prandtl-Schicht
         END DO

         !Scherungs-Antrieb der TKE:
         IF (tdc%itype_sher == 2 .AND. PRESENT(dwdx) .AND. PRESENT(dwdy)) THEN
            !Einschliesslich der 3D-Korrektur durch den Vertikalwind bzgl. der mittleren Hangneigung:
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=MAX( (zvari(i,ke1,u_m)+dwdx(i,ke1)*val_m(i))**2 &
                              +(zvari(i,ke1,v_m)+dwdy(i,ke1)*val_m(i))**2 &
                              +hdef2(i,ke1)*val_m(i)**2, fc_min(i) )
            END DO
            !Beachte:
            !'dwdx(ke1)', 'dwdy(ke1)' und 'hdef2(ke1)' beziehen sich auf die vorlaeufige Schichtdicke "1m".
            !Diese Felder sind in ICON fuer den Level 'ke1' nicht vorhanden!
         ELSE
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=MAX( zvari(i,ke1,u_m)**2+zvari(i,ke1,v_m)**2, fc_min(i) )
            END DO
         END IF
         IF (ladsshr) THEN !treatment of additional surface-shear by NTCs or LLDCs active
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=frm(i,ke1)+frc_2d(i) !full shear-forcing included transmitted shear from "P"-level
            END DO
         END IF

         ! Vertikalgradienten der dynamisch wirksamen Skalare:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            edh(i)=z1/dz_0a_h(i)   !reciprocal resistance length for momentum of turb. Prandt.-layer
            val_h(i)=tfh(i)*edh(i) !reciprocal resistance length of roughness- and laminar-layer
         END DO

        !$ACC LOOP SEQ
         DO n=tet_l, h2o_g
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               zvari(i,ke1,n)=(zvari(i,ke,n)-zvari(i,ke1,n))*edh(i)
            END DO
         END DO
         !'zvari(:,ke1,n)' enthaelt jetzt die Vertikalgradienten der Erhaltungsvariablen

         !Auftriebs-Antrieb der TKE:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            frh(i,ke1)=g_tet(i)*zvari(i,ke1,tet_l)+g_vap(i)*zvari(i,ke1,h2o_g)
         END DO

         ! Berechnung der Stabilitaetslaengen:

         IF (it_durch == it_start .AND. lini) THEN !Startinitialisierung

            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt, val1, val2, wert)
            DO i=ivstart, ivend
               IF (frh(i,ke1) >= (z1-rim)*frm(i,ke1)) THEN !die krit. Ri-Zahl wird ueberschritten
                  !'lm' sowie 'lh' werden durch 'lm' bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  tkvm(i,ke1)=l_tur_z0(i)*(sm_0-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=tkvm(i,ke1)
               ELSE
                  fakt=frh(i,ke1)/(frm(i,ke1)-frh(i,ke1))
                  tkvm(i,ke1)=l_tur_z0(i)*(sm_0-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=l_tur_z0(i)*(sh_0-a_5*fakt)
               END IF

               val1=tkvm(i,ke1)*frm(i,ke1)
               val2=tkvh(i,ke1)*frh(i,ke1)
               wert=MAX( val1-val2, rim*val1 )

               IF (.NOT.ltkeinp) THEN !TKE not present as input
                  tke(i,ke1,nvor)=MAX( SQRT(tdc%d_mom*l_tur_z0(i)*wert), tdc%vel_min )
               END IF

!                 Retrieving this peace of out-commented code:
               !Note:
               !'tkvm|h' are stability-dependent length-scales here, which should not be dependent on LLDCs.
            END DO

         ELSE ! mit Hilfe der vorhergehenden TKE-Werte

!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
            DO i=ivstart, ivend
               wert=z1/tke(i,ke1,nvor)
               tkvm(i,ke1)=tkvm(i,ke1)*wert
               tkvh(i,ke1)=tkvh(i,ke1)*wert
            END DO
         END IF

         !$ACC END PARALLEL

! 4f)    Bestimmung des neuen SQRT(2*TKE)-Wertes, der Stabilitaetsfuntionen und des SDSS:

         CALL solve_turb_budgets( tdc=tdc, it_s=it_durch, it_start=it_start,                    &

                                  i1dim=nvec, i_st=ivstart, i_en=ivend,                         & !in

                                  khi=ke1, ktp=ke-1, kcm=kcm, k_st=ke1, k_en=ke1, k_sf=ke1,     & !in

                                  ntur=ntur, nvor=nvor,                                         & !in

                                  lssintact=.FALSE.,      lupfrclim=.FALSE.,                    & !in
                                  lpres_edr=PRESENT(edr),                                       & !in
                                  ltkeinp=ltkeinp,                                              & !in

                                  imode_stke=tdc%imode_tran,  imode_vel_min=tdc%imode_vel_min,  & !in

                                  dt_tke=dt_tke, fr_tke=fr_tke,                                 & !in

                                  fm2=frm, fh2=frh, ft2=frm,                                    & !in
                                  lsm=tkvm, lsh=tkvh, tls=len_scale,                            & !in(out)

                                  tvt=tketens, velmin=velmin,                                   & !in
                                  tke=tke, ediss=ediss,                                         & !inout, out

                                  lactcnv=(tdc%icldm_tran.NE.-1 .AND. lsrflux),                 & !in (act. flux conversion)
                                  laddcnv=lsrflux,                                              & !in (add. flux-conversion)
                                  exner=zaux(:,:,1), r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3),      & !in

                                  rcld=rcld, & !inp: effective saturation fraction (cloud-cover)
                                               !out: std. deviat. of local super-saturat.
                                               !     (only for last iteration step)

                                  lcircterm=.FALSE.,                                            & !in
                                  dens=rhon, l_pat=l_pat, l_hori=l_hori,                        & !in

                                  grd=zvari,  & !inp: vert. grads. (incl. those of tet_l, h2o_g and liq)
                                                !out: vert. grads. (incl. those of tet,   vap   and liq
                                                !                   resulting from flux conversion)
                                  !'zvari'-output is only calculated at the last iteration step; and it's still equal
                                  ! to the input, if these calculations are not executed.
                                  !As "lcircterm=F", 'prss=>zvari(:,:,0)' is still near-surface pressure.

                                  lacc=lzacc, opt_acc_async_queue=acc_async_queue               ) !in

         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

 ! 4h)   Bestimmung der rein turbulenten Diffusionskoeffizienten und weiterer Turbulenzgroessen der Transferschicht:

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend

            tkvm(i,ke1)=tke(i,ke1,ntur)*tkvm(i,ke1)
            tkvh(i,ke1)=tke(i,ke1,ntur)*tkvh(i,ke1)
         END DO
         IF (tdc%imode_trancnf >= 4 .OR. it_durch < tdc%it_end) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
            DO i=ivstart, ivend
               wert=l_tur_z0(i)*SQRT(SQRT(frm(i,ke1))/tkvm(i,ke1)) !updated tkr=Ustar/(q*Sm)_0
               tkr(i)=MERGE( tdc%ditsmot*tkr(i) + (z1-tdc%ditsmot)*wert, wert, tdc%ditsmot > z0 )
            END DO
         END IF

         IF (it_durch < tdc%it_end) THEN !at least one additional iteration will take place

            IF (tdc%imode_trancnf == 2 .OR. tdc%imode_trancnf == 3) THEN
               !new version of initializing the profile-factors using Ustar,
               !but still epressing this factor in terms of "(q*Sx)_P/(q*Sx)_0":
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
               DO i=ivstart, ivend
                  fakt=h_top_2d(i)/z0m_2d(i) !(l_P-l_0)/l_0; l_0=akt*z0m
                  rat_m_2d(i)=z1+fakt*(z1-tkr(i))                                       !(q*Sm)_P/(q*Sm)_0
                  rat_h_2d(i)=z1+fakt*(z1-tkr(i)*(tkvm(i,ke1)*sh_0)/(tkvh(i,ke1)*sm_0)) !(q*Sh)_P/(q*Sh)_0
               END DO

            ELSEIF (tdc%imode_trancnf >= 4) THEN
               !new version of initializing the profile-factors and already expressing
               !them in terms of "Ustar/(q*Sh)_0*(Sh(0)/Sm(0))":

               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=ivstart, ivend
                  rat_m_2d(i)=tkr(i)                                       !Ustar/(q*Sm)_0
                  rat_h_2d(i)=tkr(i)*(tkvm(i,ke1)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(sh(0)/sm(0))
               END DO
            END IF

         END IF

         !$ACC END PARALLEL

         ! This should happen outside the OpenACC paralle region
         IF ( it_durch < tdc%it_end .AND. .NOT.ltkeinp) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF

!----------------------------------------------
      END DO !Iterationen
!----------------------------------------------

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

! 4i) Belegung der Felder fuer die Transfer-Geschwindigkeiten:

      IF (.NOT.lini .AND. tdc%ditsmot > 0) THEN
         !previous values of 'tvm|h' need to be saved
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            tcm(i)=tvm(i) !previous tvm
            tch(i)=tvh(i) !previous tvh
         END DO
      END IF
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         ! Transfer-Velocities:
         tvm(i)=tkvm(i,ke1)*val_m(i) !to be used
         tvh(i)=tkvh(i,ke1)*val_h(i) !to be used
      END DO
      IF (tdc%imode_lamdiff == 2) THEN !permanent inclusion of laminar transport
         !Securing at least laminar transport after the update of DCs:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
         DO i=ivstart, ivend
            wert=z1/(z0m_2d(i)/sai(i)+h_atm_2d(i))
            tvm(i)=MAX(con_m*wert, tvm(i))
            tvh(i)=MAX(con_h*wert, tvh(i))
         END DO
      END IF
      IF (.NOT.lini .AND. tdc%ditsmot > 0) THEN !smoothing of transfer velocity required
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            tvm(i)=tdc%ditsmot*tcm(i) + (z1-tdc%ditsmot)*tvm(i) !smoothed new tvm
            tvh(i)=tdc%ditsmot*tch(i) + (z1-tdc%ditsmot)*tvh(i) !smoothed new tvh
         END DO
      END IF
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
      DO i=ivstart, ivend
         ! Transfer-Coefficients:
         fakt=z1/vel_2d(i,ke)
         tcm(i)=tvm(i)*fakt
         tch(i)=tvh(i)*fakt
      END DO

!-----------------------------------------------
!-----------------------------------------------

! 4j) Berechnung der Enthalpie- und Impulsflussdichten sowie der EDR am Unterrand:

      IF (lsrflux.OR.lrunscm) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
         DO i=ivstart, ivend
            wert=rhon(i,ke1)*tkvh(i,ke1)

            shfl_s(i)=cp_d*wert*zvari(i,ke1,tet)*epr_2d(i)
            qvfl_s(i)=wert*zvari(i,ke1,vap)
            !Note: 'shfl_s' and 'qvfl_s' are positive downward and 'shfl_s' belogns to the T-equation!
        END DO
      END IF

      IF ((lsrflux.OR.lrunscm) .AND. (PRESENT(umfl_s).OR.PRESENT(vmfl_s))) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
         DO i=ivstart, ivend
            wert=rhon(i,ke1)*tkvm(i,ke1)
            IF (PRESENT(umfl_s)) umfl_s(i)=wert*zvari(i,ke1,u_m)
            IF (PRESENT(vmfl_s)) vmfl_s(i)=wert*zvari(i,ke1,v_m)
            !Note: 'umfl_s' and 'vmfl_s' are positive downward!
         END DO
      END IF

      !$ACC END PARALLEL

!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN
         IF (SHF%mod(0)%vst > i_cal .AND. SHF%mod(0)%ist == i_mod) THEN
            !measured SHF has to be used for forcing:
            shfl_s(imb)=SHF%mod(0)%val
         ELSEIF (lsurflu) THEN !SHF defined by explicit surface flux density
            SHF%mod(0)%val=shfl_s(imb)
            SHF%mod(0)%vst=MAX(i_upd, SHF%mod(0)%vst) !SHF is at least updated
         END IF
         IF (LHF%mod(0)%vst > i_cal .AND. LHF%mod(0)%ist == i_mod) THEN
            !measured LHF has to be used for forcing:
            qvfl_s(imb)=LHF%mod(0)%val / lh_v
         ELSEIF (lsurflu) THEN !LHF defined by explicit surface flux density
            LHF%mod(0)%val=qvfl_s(imb) * lh_v
            LHF%mod(0)%vst=MAX(i_upd, LHF%mod(0)%vst) !LHF is at least updated
         END IF
         !Note: LHF always is the latent heat flux connected with evaporation by definition,
         !      independent whether the surface is frozen or not!
      END IF
#endif
!SCLM --------------------------------------------------------------------------------

! 5)  Diagnose der meteorologischen Groessen im 2m- und 10m-Niveau:

      IF (lnsfdia) THEN !diagnostics at near surface levels required at this place

!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend

         !Einschraenkung von z0m_dia ueber Land:

         IF (fr_land(i) <= z1d2) THEN
            !Ueber See gibt es keinen synoptischen Garten
            z0d_2d(i)=z0m_2d(i)
         ELSE
            !Die Rauhigkeitslaenge einer SYNOP Station soll immer
            !kleiner als 10m bleiben:
            z0d_2d(i)=MIN( h_10m, tdc%z0m_dia )
         END IF

         !Festlegung der synoptischen Niveaus:

         IF (tdc%itype_2m_diag == 2) THEN !using an exponetial rougness layer profile
           z2m_2d (i) = h_2m -h_can_2d(i) !2m ueber dem Bodenniveau des Bestandes
           z10m_2d(i) = h_10m-z0m_2d(i)   !Hoehe, in der die turbulente Distanz 10m betraegt
         ELSE !using only a logarithmic profile above a SYNOP lawn
           z2m_2d (i) = h_2m
           z10m_2d(i) = h_10m
         END IF

         !Erste Belegung zweier benachbarter Modellniveaus:

         hk_2d(i)=h_atm_2d(i)
         hk1_2d(i)=z0
         k_2d(i)=ke

      END DO
      !$ACC END PARALLEL

!     Diagnose der 2m-Groessen:

      IF (ltst2ml) THEN !test required, whether 2m-level is above the lowest main-level
#ifdef _OPENACC
        CALL diag_level_gpu(ivstart, ivend, ke1, z2m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                            lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#else
        CALL diag_level(ivstart, ivend, ke1, z2m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                        lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#endif
      END IF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

      IF (tdc%itype_2m_diag == 2) THEN !using an exponential rougness layer profile

         val2=z1/tdc%epsi
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(val1, fakt, wert, z_surf, fac_h)
         DO i=ivstart, ivend
            IF (k_2d(i) == ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               IF (z2m_2d(i) < z0) THEN
!                 2m-Niveau liegt innerhalb der Bestandesschicht
!                 mit exponentiellen Vertikalprofilen:

                  val1=z2m_2d(i)/dz_s0_h(i)
                  IF (-val1 <= val2) THEN
                    fakt=dz_s0_h(i)/dz_sa_h(i)
                    fakt=MIN( z1, MAX( z0, fakt*EXP(val1) ) )
                  ELSE
                    fakt=z0
                  ENDIF
               ELSE
!                 2m-Niveau liegt innerhalb der Modell_Prandtl-Schicht
!                 mit logarithmischen Vertikalprofilen:

                  z_surf=z0m_2d(i)
                  fac_h=fac_h_2d(i)

                  IF (ABS(z1-fac_h) < tdc%epsi ) THEN
                     wert=z_surf*z2m_2d(i)/(z2m_2d(i)+z_surf)
                  ELSEIF (fac_h >= z0 .OR. tdc%imode_trancnf < 3) THEN
                     !non-stable strat. or using only linear interpolation of profile-function
                     !for the velocity scale (q*Sh):
                     wert=z_surf*LOG((z2m_2d(i)+z_surf)/(z_surf+fac_h*z2m_2d(i)))/(z1-fac_h)
                  ELSE !hyperbolic interpolation of profile-function (q*Sh) for stable stratification
                     wert=z2m_2d(i)/z_surf
                     IF (tdc%imode_trancnf == 3) THEN !only if the upper node for diffusion coefficients is used
                        wert=z_surf*((z1-fac_h)*LOG(wert+z1)+fac_h*wert)/(z1-fac_h)
                     ELSE !(tdc%imode_trancnf >= 4): without using the upper node for diff.-coefs.
                        wert=z_surf*(LOG(wert+z1)-fac_h*wert)/(z1-fac_h)
                     END IF
                  END IF
                  fakt=(dz_s0_h(i)+wert)/dz_sa_h(i)
               END IF

               tmps(i,ke1) = fakt*ta_2d(i) + (z1-fakt)*tl_s_2d(i)/epr_2d(i)
               vaps(i,ke1) = qt_s_2d(i) + fakt*(qda_2d(i)-qt_s_2d(i))

            END IF
         END DO

      ELSE !using only a logarithmic profile above a SYNOP lawn

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(a_atm, a_2m, fr_sd_h, val1, val2, fakt, wert, z_surf, fac_h)
         DO i=ivstart, ivend
            IF (k_2d(i) == ke .OR. .NOT.ltst2ml) THEN
               !2m-Niveau unterhalb der untersten Modell-Hauptflaeche
               !in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z_surf=z0d_2d(i)

               a_atm=h_atm_2d(i)+z_surf
               a_2m=h_2m+z_surf

               !Dimensionsloser Widerstand des Rauhigkeits-Schicht der SYNOP-Wiese (zwischen dem Unterrand
               !an der Erdoberflaeche und Oberrand bzgl. der Rauheigkeitslaenge 'z0d' einer SYNOP-Wiese
               !(wobei die turbulente Geschwindigkeitsskala und der Oberflaechenindex die fuer das gesamte
               !Gitterelemnt gueltigen Werte behalten):
               fr_sd_h=MAX( z0, dz_s0_h(i)/z0m_2d(i)+LOG(z_surf/z0m_2d(i)) )

               IF (tdc%imode_trancnf < 4) THEN !only if the  upper node for diffusion coefficients is used
                  fac_h_2d(i)=(rat_h_2d(i)-z1)*z_surf/h_top_2d(i) !re-defined profile-factor employing 'z0d'
                  !Attention:
                  !Possibly, the original profile-factor 'fac_h_2d' should be considered as a constant
                  ! of the vertical profile, rather than the re-defined one.
               END IF

               !Verhaeltnis der dimensionslosen Widerstaende:

               fac_h=fac_h_2d(i)

               IF (fac_h >= z0 .OR. tdc%imode_trancnf < 3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sh):
                  IF (fac_h == z1) THEN
                     val1=fr_sd_h+h_2m/a_2m
                     val2=fr_sd_h+h_atm_2d(i)/a_atm
                  ELSE
                     fakt=z1/(z1-fac_h)
                     val1=fr_sd_h+LOG(a_2m /(z_surf+fac_h*h_2m       ))*fakt
                     val2=fr_sd_h+LOG(a_atm/(z_surf+fac_h*h_atm_2d(i)))*fakt
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sh) for stable stratification
                  wert=z1/z_surf
                  IF (tdc%imode_trancnf == 3) THEN !only if the upper node for diffusion coefficients is used
                     fac_h_2d(i)=-fac_h/rat_h_2d(i) !transformed profile-factor for stable strat.
                     val1=fr_sd_h+(z1-fac_h)*LOG(a_2m *wert)+fac_h*h_2m       *wert
                     val2=fr_sd_h+(z1-fac_h)*LOG(a_atm*wert)+fac_h*h_atm_2d(i)*wert
                  ELSE !(tdc%imode_trancnf >= 4): without using the upper node for diff.-coefs.
                     fakt=z1/(z1-fac_h)
                     val1=fr_sd_h+(LOG(a_2m *wert)-fac_h*h_2m       *wert)*fakt
                     val2=fr_sd_h+(LOG(a_atm*wert)-fac_h*h_atm_2d(i)*wert)*fakt
                  END IF
               END IF

               fakt=val1/val2

               !Interpolationswerte fuer das synoptische 2m-Niveau:

               tmps(i,ke1) = fakt*ta_2d(i) + (z1-fakt)*tl_s_2d(i)/epr_2d(i)
               vaps(i,ke1) = qt_s_2d(i) + fakt*(qda_2d(i)-qt_s_2d(i))
               IF (tdc%icldm_tran > -1) THEN !water phase change is possible
                  fakt=h_2m/h_atm_2d(i)
                  rcls(i,ke1)=rcld(i,ke1)+fakt*(rcld(i,ke)-rcld(i,ke1))
               END IF
            END IF
         END DO
      END IF !using only a logarithmic profile above a SYNOP lawn

      !$ACC END PARALLEL

      IF (ltst2ml) THEN !test required, whether 2m-level is above the lowest main-level
      !$ACC KERNELS ASYNC(acc_async_queue) DEFAULT(PRESENT) IF(lzacc)
      k=MINVAL(k_2d(ivstart:ivend))
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      IF (k < ke) THEN !2m-level is above the lowest main-level at least for one grid point
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(k2, k1, fakt, wert, val1, val2)
      DO i=ivstart, ivend
         IF (k_2d(i) < ke) THEN
!           2m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!           trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommenen
!           linearen Profile der progn. Modellvariablen eine logarith. Interpolation:

            k2=k_2d(i); k1=k2+1

            fakt=z1/(hk1_2d(i)+z0d_2d(i))
            wert=(h_2m    +z0d_2d(i))*fakt
            fakt=(hk_2d(i)+z0d_2d(i))*fakt
            fakt=LOG(wert)/LOG(fakt)

            IF (tdc%icldm_tran == -1) THEN !no water phase change possible
               val2=qv(i,k2)          ; val1= qv(i,k1)        ; vaps(i,ke1)=val1+fakt*(val2-val1)
               val2= t(i,k2)/epr(i,k2); val1=t(i,k1)/epr(i,k1); tmps(i,ke1)=val1+fakt*(val2-val1)
            ELSE !water phase changes are possible
               val2=qv(i,k2)+qc(i,k2) ; val1=qv(i,k1)+qc(i,k1); vaps(i,ke1)=val1+fakt*(val2-val1)
               val2=(t(i,k2)-lhocp*qc(i,k2))/epr(i,k2); val1=(t(i,k1)-lhocp*qc(i,k1))/epr(i,k1)
                                                               tmps(i,ke1)=val1+fakt*(val2-val1)
               rcls(i,ke1)=rcld(i,k1)+fakt*(rcld(i,k2)-rcld(i,k1))
            END IF
         END IF
         !Note: In the case "k_2d(i) == ke" 'tmps' and 'vaps' have already been
         !      calculated above.
      END DO
      END IF !At least for one grid-point, the 2m-level is above the lowest main-level.
      !$ACC END PARALLEL
      END IF !test required, whether 2m-level is above the lowest main-level

      !Druck im 2m-Niveau:

!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
      DO i=ivstart, ivend
!test<
         wert=tmps(i,ke1)*(z1+rvd_m_o*vaps(i,ke1)) !angenaeherte virt. Temp.
         prss(i,ke1)=prss(i,ke1)                &  !Druck
                    *EXP(-(z2m_2d(i)-hk1_2d(i))*grav/(r_d*wert))
!prss(i,ke1)=prss(i,ke1)-(z2m_2d(i)-hk1_2d(i))*grav*rhon(i,ke1)
!test>
      END DO
      !$ACC END PARALLEL

      IF (tdc%icldm_tran == -1) THEN
!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
             t_2m(i)=tmps(i,ke1)
            qv_2m(i)=vaps(i,ke1)
         END DO
         !$ACC END PARALLEL
      ELSE
!        Berechnung der zugehoerigen Modell- und Feuchtevariablen im 2m-Niveau
!        aus den Erhalturngsvariablen.

         CALL adjust_satur_equil( tdc=tdc, i1dim=nvec, khi=ke1, ktp=ke-1,      & !in
!
              i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1,                    & !in
!
              lcalrho=.FALSE., lcalepr=.TRUE., lcaltdv=.FALSE.,                & !in
              lpotinp=.TRUE. , ladjout=.TRUE.,                                 & !in
!
              icldmod=tdc%icldm_tran,                                          & !in
!
              zrcpv=tur_rcpv, zrcpl=tur_rcpl,                                  & !in
!
              prs=prss, t=tmps, qv=vaps,                                       & !in (pres. and conserved variabs. at 2m-level)
!
              psf=ps,                                                          & !in
!
              ! Note: The follwing REAL variables will be overwritten by values valid for the 2m-level:
!
              rcld=rcls,  & !inp: std. deviat. of local super-saturat. at 2m-level
                            !out: saturation fraction (cloud-cover)    at 2m-level (not yet used)
!
              exner=zaux(:,:,1),                                               & !out (not usd)
                                                                                 !aux (internally)
!
              tet_liq=zvari(:,:,tet), q_h2o=zvari(:,:,vap),                    & !inp: as target of 'tmps, vaps'
                                                                                 !out: adjusted  variables at 2m-level
                                      q_liq=zvari(:,:,liq),                    & !out: cloud-water of 2m-fog (not yet used)
!
              lacc=lzacc, opt_acc_async_queue=acc_async_queue )

!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
             t_2m(i)=zvari(i,ke1,tet)
            qv_2m(i)=zvari(i,ke1,vap)
         END DO
         !$ACC END PARALLEL
      END IF

      IF (.NOT.tdc%lfreeslip) THEN !not for idealized dry runs with free-slip condition

!        Diagnose der 10m-Groessen:

#ifdef _OPENACC
         CALL diag_level_gpu(ivstart, ivend, ke1, z10m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                             lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#else
         CALL diag_level(ivstart, ivend, ke1, z10m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                         lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#endif


         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR PRIVATE(a_atm, a_10m, val1, val2, fakt, wert, z_surf, fac_m)
         DO i=ivstart, ivend

            IF (k_2d(i) == ke) THEN

!              10m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z_surf=z0d_2d(i)

               a_atm=h_atm_2d(i)+z_surf
               a_10m=h_10m+z_surf

               IF (tdc%imode_trancnf < 4) THEN !further corrected profile-factor using an upper node
                  fac_m_2d(i)=(rat_m_2d(i)-z1)*z_surf/h_top_2d(i)
                  !Attention:
                  !Possibly, the original profile-factor 'fac_m_2d' should be considered as a constant
                  ! of the vertical profile, rather than the re-defined one.
               END IF

               !Verhaeltnis der dimensionslosen Widerstaende:

               fac_m=fac_m_2d(i)

               IF (fac_m >= z0 .OR. tdc%imode_trancnf < 3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sm):
                  IF (fac_m == z1) THEN
                     val1=h_10m/a_10m
                     val2=h_atm_2d(i)/a_atm
                  ELSE
                     val1=LOG(a_10m/(z_surf+fac_m*h_10m))
                     val2=LOG(a_atm/(z_surf+fac_m*h_atm_2d(i)))
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sm) for stable stratification
                  wert=z1/z_surf
                  IF (tdc%imode_trancnf == 3) THEN !only if the upper node for diffusion coefficients is used
                     fac_m_2d(i)=-fac_m/rat_m_2d(i) !transformed profile-factor for stable stratification
                     val1=(z1-fac_m)*LOG(a_10m*wert)+fac_m*h_10m      *wert
                     val2=(z1-fac_m)*LOG(a_atm*wert)+fac_m*h_atm_2d(i)*wert
                  ELSE !(tdc%imode_trancnf >= 4): without using the upper node for diff.-coefs.
                     val1=LOG(a_10m*wert)-fac_m*h_10m      *wert
                     val2=LOG(a_atm*wert)-fac_m*h_atm_2d(i)*wert
                  END IF
               END IF

               fakt=val1/val2

               !Interpolationswerte fuer das synoptische 10m-Niveau:

               u_10m(i)=vel1_2d(i)*fakt; v_10m(i)=vel2_2d(i)*fakt

            END IF
         END DO

!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR PRIVATE(fakt, wert, k1, k2)
         DO i=ivstart, ivend
            IF (k_2d(i) < ke) THEN
!              10m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!              trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommen
!              linearen Profile der progn. Modellvariablen eine logarithm. Interpolation:

               IF (ltst10ml) THEN
                  k2=k_2d(i); k1=k2+1
               ELSE
                  k2=ke-1; k1=ke
               END IF

               fakt=z1/(hk1_2d(i)+z0d_2d(i))
               wert=(h_10m   +z0d_2d(i))*fakt
               fakt=(hk_2d(i)+z0d_2d(i))*fakt
               fakt=LOG(wert)/LOG(fakt)
               u_10m(i)=u(i,k1)+fakt*(u(i,k2)-u(i,k1))
               v_10m(i)=v(i,k1)+fakt*(v(i,k2)-v(i,k1))
            END IF
         END DO

         !$ACC END PARALLEL

      END IF !not for idealized dry runs with free-slip condition

      END IF !diagnostics at near surface levels required at this place

      !Notes:
      !'[u|v]_10m' are ialways required for extended calculation of sea-surface roughness-length.

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

      IF (tdc%lfreeslip) THEN !only for idealized dry runs with free-slip condition
!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
            qv_2m(i)=z0
            rh_2m(i)=z0
            td_2m(i)=z0
            u_10m(i)=u(i,ke)
            v_10m(i)=v(i,ke)
         END DO
      ELSEIF (.NOT.lnsfdia) THEN !neither "lfreeslip=T" nor a particular near-surface-diagnostics
!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
             t_2m(i)= t(i,ke)
            qv_2m(i)=qv(i,ke)
            u_10m(i)=u(i,ke)
            v_10m(i)=v(i,ke)
         END DO
      ELSE !usual near-surface-diagnostics active

!        Finale 2m-Diagnose:
!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR PRIVATE(patm, fakt, wert)
         DO i=ivstart, ivend
            patm=prss(i,ke1)*qv_2m(i) &
                /(rdv+(z1-rdv)*qv_2m(i))          !Wasserdampfdruck

            fakt=patm/zpsat_w( t_2m(i) )
            rh_2m(i)=100.0_wp*MIN( fakt, z1 )     !relative Feuchte

            !UB: old formulation
            wert=LOG(patm/b1)
            !UB: For dry atmosphere, the Teten's formula is not defined at vapor
            !    pressure patm=0 because of log(patm/b1). However, it converges
            !    to a dew point of the value of parameter b4w in case we impose
            !    a very small positive lower bound on  patm:
            ! erst mal wieder weggenommen, da man den Absturz will:
            !wert=LOG(MAX(patm,1.0E-16_wp)/b1)

            td_2m(i)=MIN( (b2w*b3-b4w*wert) &
                         /(b2w-wert), t_2m(i) )   !Taupunktstemperatur
         END DO
      END IF

      IF (ladsshr) THEN !treatment of additional surface-shear by NTCs or LLDCs active
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
         DO i=ivstart, ivend
            fakt=SQRT( frm(i,ke1)/(frm(i,ke1)-frc_2d(i)) ) !amplification-factor for near-surface wind-speed
                                                           ! due to shear by NTCs|LLDCs
            vel_2d(i,ke)=fakt*vel_2d(i,ke) !amplified wind-speed at level "A"
            IF (tdc%imode_nsf_wind == 2) THEN
               !Amplification of near-surface wind-diagnostics by NTCs or LLDCs:
               u_10m(i)=fakt*u_10m(i); v_10m(i)=fakt*v_10m(i)

               !Note:
               !These amplified values are also applied for extended calculation of sea-surface roughnes-length.
            END IF
         END DO
      END IF

      IF (igz0inp /= 1 .OR. lini) THEN
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(vel2, wert, fakt)
      DO i=ivstart, ivend

!        Diagnose von 'gz0' (fuer den naechsten Zeitschritt)
!        ueber Wasserflaechen mit der (anepassten) Charnock-Formel

         IF (fr_land(i) <= z1d2) THEN !a non-land surface
            !Note: This definition of a non-land surface is now in line with the ICON-definition
            !       using "frlnd_thrhld=z1d2"

            IF ( l_sice(i) ) THEN !ice-covered surface
               gz0(i)=g_z0_ice
            ELSE !water-covered surface
               vel2=(tke(i,ke,nvor)*z1d2)**2 !squared sub-grid wind complement
               !Note:
               !This refers to linear interpolation of turbulent velocity at the lowest atmospheric boundary-level
               ! onto the lowest main-level, with vanishing magnitude at "0"-level.
               !The length-scales contributing to 'tke(i,ke1,ntur)' are assumed to be too small as to generate
               ! surface waves that can consturct the roughness layer below. Hence, 'vel2' is not affected by them.
               !However, through 'tke(i,ke,nvor)', Ustar may also develop without vertical shear of mean wind,
               ! (e.g. at the convective limit).
               !In case of "ladsshr=T", 'vel_2d' and 'tvm' already contain some contribution by NTCs|LLDCs due
               ! to the impact of related additional shear at level "P" that has been transmitted to level "0".
               !Since 'vel2' is used here as a substitute for the transmission of additional shear at level "P",
               ! the impact of 'vel2' on effective Ustar**2 has to be restricted to the fraction "1-rsur_sher"
               ! of not transmitted shear-forcing.
               !At "ladsshr=T" and "imode_nsf_wind=2", even '[u|v]_10m' are affected by NTCs|LLDCs.

               !Effective squared friction velocity Ustar**2:
               vel2=tvm(i)*MAX( tdc%vel_min, SQRT( vel_2d(i,ke)**2 + (z1-tdc%rsur_sher)*vel2 ) )

               ! Charnock-parameter:

               ! Basic Charnock-parameter; use enhanced value of "0.1" over lakes if "imode_charpar>1":
               fakt=MERGE( tdc%alpha0, 0.1_wp, tdc%imode_charpar == 1 )

               ! use velocity-dependent Charnock paramter over sea if "imode_charpar>1":
               fakt=MERGE( alpha0_char( SQRT( u_10m(i)**2+v_10m(i)**2 ), tdc ), &
                           fakt, tdc%imode_charpar > 1 .AND. .NOT.l_lake(i) )

               !Substituting dynamic (shear-related) part of diagnosed R-length by 'z0_waves', if present and singnificant:
               wert=MERGE( grav*z0_waves(i), fakt*vel2, (igz0inp == 2 .AND. z0_waves(i) >= z0wave_threshold) )

               !Diagnosed sea-surface R-length including a laminar correction:
               wert=MAX( g_len_min, wert + g_alpha1_con_m/SQRT(vel2) )
               !Final 'gz0' with optional time-step smoothing:
               gz0(i)=MERGE( tdc%ditsmot*gz0(i)+(z1-tdc%ditsmot)*wert, wert, tdc%ditsmot > z0 )
            END IF
         END IF
      END DO
      ENDIF  !igz0inp

      !$ACC END PARALLEL

      !$ACC END DATA ! from acc data present
      !$ACC END DATA ! from acc data create

END SUBROUTINE turbtran

!==============================================================================

!+ Module procedure 'diag_level' for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level (i_st, i_en, ke1, zdia_2d, hhl, k_2d, hk_2d, hk1_2d, lacc, opt_acc_async_queue)
   INTEGER, INTENT(IN) :: &
!
      ke1, &
      i_st, i_en  !start end end indices of horizontal domain

   REAL (KIND=wp), INTENT(IN) :: &
!
      zdia_2d(:)  !diagnostic height

   INTEGER, INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                  !to be used for near surface diagnostics

   REAL (KIND=wp), INTENT(IN) :: &
!
     hhl(:,:)

   REAL (KIND=wp), INTENT(INOUT) :: &
!
     hk_2d(:), & !mid/full-level height above ground belonging to 'k_2d'
     hk1_2d(:)   !mid/full-level height above ground of the previous layer (below)

   LOGICAL, OPTIONAL, INTENT(IN) :: lacc ! If true, use openacc
   INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

   LOGICAL :: lzacc ! non-optional version of lacc
   INTEGER :: acc_async_queue

   INTEGER :: i

   LOGICAL :: lcheck

   CALL set_acc_host_or_device(lzacc, lacc)

   IF(PRESENT(opt_acc_async_queue)) THEN
       acc_async_queue = opt_acc_async_queue
   ELSE
       acc_async_queue = 1
   ENDIF

   lcheck=.TRUE. !check whether a diagnostic level is above the current layer

   !XL_ACCTMP:this could be implemented with an explcit K loop, may be faster on GPU (no atomic)

   DO WHILE (lcheck) !loop while previous layer had to be checked
      lcheck=.FALSE. !check next layer ony, if diagnostic level is at least once
                     !above the current layer

      !$ACC PARALLEL ASYNC(acc_async_queue) PRESENT(hhl, zdia_2d, k_2d, hk_2d, hk1_2d) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=i_st,i_en
         IF (hk_2d(i)<zdia_2d(i) .AND. k_2d(i)>1) THEN !diagnostic level is above current layer
            !$ACC ATOMIC WRITE
            lcheck=.TRUE. !for this point or any previous one, the diagnostic level is above the current layer
            !$ACC END ATOMIC
            k_2d(i)=k_2d(i)-1
            hk1_2d(i)=hk_2d(i)
            hk_2d(i)=(hhl(i,k_2d(i))+hhl(i,k_2d(i)+1))*z1d2-hhl(i,ke1)
          END IF
       END DO
       !$ACC END PARALLEL
   END DO

END SUBROUTINE diag_level

!==============================================================================

!+ Module procedure 'diag_level_gpu' for computing the upper level index
!+ used for near surface diganostics (GPU-version)

SUBROUTINE diag_level_gpu (i_st, i_en, ke1, zdia_2d, hhl, k_2d, hk_2d, hk1_2d, lacc, opt_acc_async_queue)

   INTEGER, INTENT(IN) :: &
!
      ke1, &
      i_st, i_en  !start end end indices of horizontal domain

   REAL (KIND=wp), INTENT(IN) :: &
!
      zdia_2d(:)  !diagnostic height

   INTEGER, INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                  !to be used for near surface diagnostics

   REAL (KIND=wp), INTENT(IN) :: &
!
     hhl(:,:)

   REAL (KIND=wp), INTENT(INOUT) :: &
!
     hk_2d(:), & !mid/full-level height above ground belonging to 'k_2d'
     hk1_2d(:)   !mid/full-level height above ground of the previous layer (below)

   LOGICAL, OPTIONAL, INTENT(IN) :: lacc
   INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

   LOGICAL :: lzacc
   INTEGER :: acc_async_queue

   INTEGER :: i, k, ke1_war

   CALL set_acc_host_or_device(lzacc, lacc)

   IF(PRESENT(opt_acc_async_queue)) THEN
       acc_async_queue = opt_acc_async_queue
   ELSE
       acc_async_queue = 1
   ENDIF

   ! NV HPC 23.3 workaround
   ke1_war = ke1
   ! Need to keep data section separate for now
   !$ACC DATA PRESENT(i_en) IF(lzacc)
   !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
   DO i=i_st,i_en
      IF (i>i_en) CYCLE ! NVHPC compiler WAR
      !$ACC LOOP SEQ
      DO k=k_2d(i)-1, 0, -1
         IF (hk_2d(i)<zdia_2d(i) .AND. k_2d(i)>1) THEN !diagnostic level is above current layer
            k_2d(i)=k
            hk1_2d(i)=hk_2d(i)
            hk_2d(i)=(hhl(i,k_2d(i))+hhl(i,k_2d(i)+1))*z1d2-hhl(i,ke1_war)
         ELSE
            EXIT
         END IF
      END DO
   END DO
   !$ACC END DATA

END SUBROUTINE diag_level_gpu

!==============================================================================

ELEMENTAL FUNCTION alpha0_char(u10, tdc)

   !$ACC ROUTINE SEQ

   ! Wind-speed dependent specification of the Charnock parameter based on suggestions by
   ! Jean Bidlot and Peter Janssen, ECMWF
   REAL (KIND=wp), INTENT(IN) :: u10 ! 10 m wind speed
   TYPE(t_turbdiff_config), INTENT(IN) :: tdc

   REAL (KIND=wp), PARAMETER  :: a=6.e-3_wp, b=5.5e-4_wp, &
                                 c=4.e-5_wp, d=6.e-5_wp,  &
                                 u2=17.5_wp, umax=40.0_wp
   REAL (KIND=wp) :: ulim, ured, alpha0_char

   ulim = MIN( u10, umax )
   ured = MAX( 0._wp, ulim - u2 )
   alpha0_char = MIN( tdc%alpha0_max, MAX( tdc%alpha0, a + tdc%alpha0_pert + ulim*(b + c*ulim - d*ured) ) )
   alpha0_char = MERGE( MIN( alpha0_char, 0.8_wp/MAX( 1._wp, u10 ) ), alpha0_char, tdc%imode_charpar==3 )

 END FUNCTION alpha0_char

!==============================================================================

END MODULE turb_transfer
