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

! Source module for turbulence utility routines

MODULE  turb_utilities

!==============================================================================
!
! Description:
!
!   Routines (module procedures) currently contained:
!     - turb_setup           : setting up the turbulence model
!     - adjust_satur_equil   : sub-grid scale moist physics in terms of a
!                              statistical saturation adjustment
!     - solve_turb_budgets   : solution of prognostic TKE-equation and the reduced
!                              linear system of the other diagnostic 2-nd-order equations
!                              including the calculation of Standard Deviatin of of Super-Saturation (SDSS)
!                              and of the near-surface circulation acceleration
!     - turb_cloud           : statistical cloud scheme called by 'adjust_satur_equil'
!     - turb_stat            : calculation of all 2-nd order moments (only for SC-diagnostics)
!     - vert_grad_diff       : main routine organizing the general calculation of
!                              semi-implicit vertical diffusion
!     - prep_impl_vert_diff  : setting up the tridiagonal system and performing some
!                              common calculations being independent on the diffused variable
!     - calc_impl_vert_diff  : final solutioin of the tridiagonal system
!     - vert_smooth          : vertical smoothing
!     - bound_level_interp   : interpolation of full level variables onto half-levels
!
!   and FUNCTIONs
!     - alpha0_char(u10)     : derived Charnock-parameter dependent on near surface wind speed
!     - zexner (zpres)       : Exner-factor
!     - zqvap (zpvap, zpdry) : satur. specif. humid. (new version)
!     - zpsat_w (ztemp)      : satur. vapor pressure over water
!     - zdqsdt (ztemp, zqsat): d_qsat/d_tem (new version)
!
! Documentation of changes to former versions of these subroutines:
!
!  Initial Release, based on SUB's and FUCNTIONS that had been originally contained in a single
!   module 'organize_turbdiff' of the common blocked schemes for atmospheric turbulence, vertical
!   diffusion and surface-to-atmosphere transfer after some major reorganization.
!   These routines are mainly called by SUB 'turbdiff' and SUB 'turbtran' contained in the MODULEs
!   'turb_diffusion' and 'turb_transfer' respectively.  In particular the common blocked routines for
!   semi-implicit vertical diffusion are now conatained in this MODULE.
!   Partly new (more consistent) interpretation of already existing selectors and some
!   further parameters gradually controlling numerical restrictions have been introduced.
!
!-------------------------------------------------------------------------------------------------------!
!  The history of all these modifications is as follows, where those belonging to the fomal
!   reorganization of the whole package (atmospheric turbulence and surface-to-atmpsphere transfer)
!   are now in the header of MODULE 'turb_utilities', containing various common SUBs for 'turbdiff'
!   and 'turtran' (related moist thermodynamicds and the treatment of turbulent 2-nd order budget equations)
!   and also the blocked code for semi-implicit vertical diffusion. The new blocked version of SUB 'turbtran'
!   is now in MODULE 'turb_transfer':
!
!              2010/12/30 Matthias Raschendorfer
!  Reorganization of the code for use in both models COSMO and ICON with various formal modifications:
!   MODULE 'src_turbdiff' CONTAINS now SUBs 'init_canopy', 'organize_turbdiff' and 'turb_cloud'.
!   The latter was before in 'meteo_utilities' with name 'cloud_diag' and contains a saturation
!   adjustment with due regard on turbulent fluctuations of thermodynamic model variables.
!   SUB 'organize_turbdiff' CONTAINS SUBs 'turbtran', 'turbdiff', 'stab_funct', 'diag_level'
!   and ('canopy_source'). Exept 'diag_level' they had been present before as well, partly in INCLUDE-files.
!   In accordance with ICON rules, SUB 'organize_turbdiff' is CALLed with the most needed parameters
!   in the SUB header. Some parameters are OPTIONAL, allowing to run the scheme in differen modes.
!              2011/06/17 Matthias Raschendorfer
!  Removing a bug related to the saving of the saturation fraction that had been introduced by
!  the ICON modifications.
!              2011/08/26 Matthias Raschendorfer
!  Introduction of a preconditioning of the tridiagonal system for "imode_turb=3".
!              2011/09/28 Matthias Raschendorfer
!  Correction in formula for specific humidity at saturation by using partial pressure of dry air.
!              2011/12/08 Matthias Raschendorfer
!  Introduction of SUBs 'prep_impl_vert_diff' and 'calc_impl_vert_diff' substituting the code segment for
!   implicit vertical diffusion calcul. for TKE or all other variables, if 'imode_turb' is one of 3 or 4.
!   These SUBs allow the calculation of semi implicit diffusion for arbitrary variables defined on
!   full- or  half-levels.
!              2011/12/27 Matthias Raschendorfer
!  Correction of some bugs:
!   Including missing definition of 'km1' for the surface level in SUB calc_impl_vert_diff
!   Moving a wrong ')' in interpolation of 't' onto the lower boundary of model atmosphere.
!              2012/03/20 Matthias Raschendorfer
!  Rearrangement, modularization and revision of numerical treatment, such as:
!  Intoduction of the SUBs 'solve_turb_budgets' and 'adjust_satur_equil' in order to use the same code
!  for the same purpose in 'turbtran' and 'turbdiff.
!  Introducing the new driving SUB 'vert_grad_diff' oranizing vertical diffusion.
!              2014/07/28 Matthias Raschendorfer
!  Some formal cleaning mainly in order to avoid code duplication
!   -> no influence of results, except due to non-associative multiplication using 'rprs'
!      and non-reversible calculation of the reciprocal: "z1/(temp/exner).NE.exner/temp"
!      both in SUB 'adjust_satur_equil'
!              2015/08/25 Matthias Raschendorfer
! Adopting other development within the ICON-version (namely by Guenther Zaengl) as switchable options
!  related to the following new selectors and switches:
!   imode_pat_len, imode_frcsmot and  lfreeslip.
! Rearranging the development by Matthias Raschendorfer that had not yet been transferred to COSMO as switchable
!  options related to the following switches:
!  and selectors:
!   imode_stbcorr, imode_qvsatur, imode_stadlim,
!  and a partly new (more consistent) interpretation of:
! Controlling numerical restrictions gradually namely by the parameters:
!  tndsmot, tkesmot, stbsmot, frcsecu, tkesecu, stbsecu
! Adopting the 3D-options from the COSMO-Version and correcting the horizontal limit of the turbulent length scale.
!              2016-05-10 Ulrich Schaettler
! Splitting this module from the original module 'organize_turbdiff' as it was used by ICON before.
! Moving declarations, allocation and deallocations of ausxilary arrays into MODULE 'turb_data'.
!
!==============================================================================
!
! Declarations:
!
! Modules used:
!------------------------------------------------------------------------------

USE mo_kind, ONLY: wp ! KIND-type parameter for real variables

USE mo_turbdiff_config, ONLY : &

    t_turbdiff_config, &

#ifdef SCLM
    nmvar,        & ! number of dynamically active model variables

    u_m     ,     & ! zonal velocity-component at the mass center
    v_m     ,     & ! meridional ,,      ,,    ,, ,,   ,,    ,,
    w_m             ! vertical   ,,      ,,    ,, ,,   ,,    ,, (only used for optional single-column diagnostics)
#endif

    ntyp    ,     & ! number of variable-types ('mom' and 'sca')
    mom     ,     & ! index for a momentum variable
    sca     ,     & ! index for a scalar   variable
    vap     ,     & ! water vapor
    liq     ,     & ! liquid water
    tet     ,     & ! potential temperature
    tet_l   ,     & ! liquid-water potentional temperture
    h2o_g   ,     & ! total water content

    !Note: It always holds: "tem=tet=tet_l=tem_l" and "vap=h2o_g" (respective usage of equal indices)!

    lporous ,     & ! vertically resolved roughness layer representing a porous atmospheric medium

    ! derived parameters calculated in 'turb_setup'
    c_tke,tet_g,rim,                  &
    c_m,c_h, b_m,b_h,  sm_0, sh_0,    &
    d_1,d_2,d_3,d_4,d_5,d_6,          &
    a_3,a_5,a_6,                      &
    tur_rcpv, tur_rcpl,               &

    ! used data-types
    varprf !variable profile

USE mo_physical_constants, ONLY : &
!
! Physical constants and related variables:
! -------------------------------------------
!
    r_d      => rd,       & ! gas constant for dry air
    rdv,                  & ! r_d / r_v
    o_m_rdv,              & ! 1 - r_d/r_v
    rvd_m_o  => vtmpc1,   & ! r_v/r_d - 1
    cp_d     => cpd,      & ! specific heat for dry air
    rcpv,                 & ! cp_v/cp_d - 1
    rcpl,                 & ! cp_l/cp_d - 1 (where cp_l=cv_l)
    rdocp    => rd_o_cpd, & ! r_d / cp_d
    lhocp    => alvdcp,   & ! lh_v / cp_d
    con_m,                & ! kinematic vsicosity of dry air (m2/s)
!
    grav,                 & ! acceleration due to gravity
    p0ref,                & ! reference pressure for Exner-function
    b3       => tmelt,    & ! temperature at melting point
!
    uc1, ucl                ! params. used in cloud-cover diagnostics based on rel. humidity

USE mo_math_constants, ONLY : &

    uc2      => sqrt3       ! SQRT(3) used in cloud-cover diagnostics based on rel. humid.

USE mo_lookup_tables_constants, ONLY : &
!
! Parameters for auxilary parametrizations:
! ------------------------------------------
!
    b1       => c1es,     & ! variables for computing the saturation steam pressure
    b2w      => c3les,    & ! over water (w) and ice (e)
    b2i      => c3ies,    & !               -- " --
    b4w      => c4les,    & !               -- " --
    b4i      => c4ies,    & !               -- " --
    b234w    => c5les       ! b2w * (b3 - b4w)
!
USE mo_lnd_nwp_config, ONLY: lterra_urb, itype_eisa
!
!------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global,    ONLY : &
!
    lsclm, i_cal, imb, &
!
    UUA, VVA, UWA, VWA, WWA, UST, TWA, QWA, TTA, TQA, QQA, &
    TKE_SCLM=>TKE, BOYPR, SHRPR, DISSI, TRANP
#endif
!SCLM---------------------------------------------------------------------------
USE mo_fortran_tools, ONLY: set_acc_host_or_device

#if defined(_OPENACC) && (__NVCOMPILER_MAJOR__ <= 21)

! Legacy Nvidia compilers are buggy and fail with attach in bound_level_interp
#define _PGI_LEGACY_WAR 1
USE openacc
#endif

USE mo_atm_phy_nwp_config, ONLY: lcuda_graph_turb_tran
!==============================================================================

IMPLICIT NONE

PUBLIC adjust_satur_equil, solve_turb_budgets, vert_grad_diff,   &
       prep_impl_vert_diff, calc_impl_vert_diff,                 &
       vert_smooth, bound_level_interp,                          &
       zbnd_val, zexner, zpsat_w,                                &
       turb_setup

REAL (KIND=wp), PARAMETER :: &
!
    z0   = 0.0_wp, &
    z1   = 1.0_wp, &
    z2   = 2.0_wp, &
    z3   = 3.0_wp, &
    z1d2 = z1/z2 , &
    z1d3 = z1/z3

REAL (KIND=wp), POINTER :: &
    a_h, a_m, d_h, d_m

!==============================================================================

CONTAINS

!==============================================================================

SUBROUTINE turb_setup ( tdc, i_st, i_en, k_st, k_en, &
                        iini, dt_tke, nprv, l_hori, &
                        ps, t_g, qv_s, qc_a, &
                        lini, it_start, nvor, fr_tke, l_scal, fc_min, &
                        prss, tmps, vaps, liqs, rcld, &
                        lacc, opt_acc_async_queue )

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

INTEGER, INTENT(IN) :: &
!
    k_st, k_en, & ! vertical   start- and end-index
    i_st, i_en, & ! horizontal start- and end-index
!
    iini,    & ! type of initialization (0: no, 1: separate before the time loop
               !                                2: within the first time step)
    nprv       ! previous time step index of 'tke'

REAL (KIND=wp) :: dt_tke !time step for TKE

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
    l_hori, &  ! horizontal grid spacing (m)
    ps,     &  ! surface pressure [Pa]
    t_g,    &  ! mean surface temperature [K]
    qv_s,   &  ! specific humidity at the surface
    qc_a       ! liquid water content of the lowermost model layer

LOGICAL, INTENT(OUT) :: lini    !initialization required

INTEGER, INTENT(OUT) ::  &
!
     it_start, & ! start index of iteration
     nvor        ! running time step index for previous 'tke' (will be updated in case of iterations)

REAL (KIND=wp), INTENT(OUT) :: fr_tke !1/dt_tke

REAL (KIND=wp), DIMENSION(:), INTENT(OUT) :: & !always evaluated
!
     l_scal, & ! reduced maximal turbulent length scale due to horizontal grid spacing (m)
     fc_min, & ! minimal value for TKE-forcing (s-2)
     ! effective surface-value within profiles of thermodynamic variables including the surface level:
     prss,   & ! pressure
     tmps,   & ! temperature
     vaps,   & ! specific humidity
     liqs      ! liquid water content at the surface

REAL (KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: & !modified only for initialization
!
     rcld      ! standard deviation of local super-saturation

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

LOGICAL :: lzacc
INTEGER :: acc_async_queue

INTEGER :: i,k

!-------------------------------------------------------------------------------
  CALL set_acc_host_or_device(lzacc, lacc)

  IF(PRESENT(opt_acc_async_queue)) THEN
      acc_async_queue = opt_acc_async_queue
  ELSE
      acc_async_queue = 1
  ENDIF

  fr_tke=z1/dt_tke

  a_h => tdc%a_heat ! length-scale factor for turbulent heat transport
  a_m => tdc%a_mom  ! length-scale factor for turbulent momentum transport
  d_h => tdc%d_heat ! length-scale factor for turbulent heat dissipation
  d_m => tdc%d_mom  ! length-scale factor for turbulent momentum dissipation

  !$ACC DATA PRESENT(l_hori, ps, t_g, qv_s, qc_a) &
  !$ACC   PRESENT(l_scal, fc_min, rcld) &
  !$ACC   PRESENT(prss, tmps, vaps, liqs) &
  !$ACC   COPYIN(i_en) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i=i_st, i_en
     l_scal(i)=MIN( z1d2*l_hori(i), tdc%tur_len )
 !__________________________________________________________________________
 !test: frm ohne fc_min-Beschraenkung: Bewirkt Unterschiede!
     fc_min(i)=(tdc%vel_min/MAX( l_hori(i), tdc%tur_len ))**2
 !   fc_min(i)=z0
 !__________________________________________________________________________
     prss(i)=ps(i)
     tmps(i)=t_g(i)
     vaps(i)=qv_s(i)
  END DO

  IF (tdc%ilow_def_cond.EQ.2) THEN !zero surface value of liquid water
!DIR$ IVDEP
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
        liqs(i)=z0
     END DO
  ELSE !constant liquid water within the transfer-layer
!DIR$ IVDEP
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
        liqs(i)=qc_a(i)
     END DO

     !Note: This setting belongs to a zero-flux condition at the surface.
  END IF

  !$ACC END PARALLEL

  IF (iini.GT.0) THEN !an initialization run
     lini=.TRUE.
     IF (iini.EQ.1) THEN !separate initialization before the time loop
        it_start=1 !only 'it_end' iteration steps for initialization
                   !and an additional one at the first time loop
     ELSE !initialization within the first time step
        it_start=0 !"it_end+1" iteration steps for initialization
     END IF

 !   Initializing the standard-deviation of super-saturation (SDSS):
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
     !$ACC LOOP SEQ
     DO k=k_st, k_en
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=i_st, i_en
           rcld(i,k)=z0 !no standard-deviat. of local super-saturation for all provided levels
        END DO
     END DO
     !$ACC END PARALLEL

  ELSE !not an initialization run
     lini=.FALSE.
     it_start=tdc%it_end !only a single iteration step
  END IF

  nvor=nprv !Eingangsbelegung von 'nvor' (wird bei Iterationen auf 'ntur' gesetzt)

 !Note:
 !A call with "iini=2" (at the first time step) provides the same result as a
 !  call with "iini=1" (before the   time loop) followed by a second
 !  call with "iini=0" (at the first time step).

 !US and these parameters are also needed after restart!!!!

 !   Calculation of derived parameters of the turbulence model:

     tet_g=grav/cp_d !adiabatic T-gradient

     c_tke=d_m**z1d3 !=exp( log(d_m)/3.0_wp )

     c_m=1.0_wp-1.0_wp/(a_m*c_tke)-6.0_wp*a_m/d_m !=3*0.08
     c_h=0.0_wp !kann auch als unabhaengiger Parameter aufgefasst werden

     b_m=1.0_wp-c_m
     b_h=1.0_wp-c_h

     d_1=1.0_wp/a_h
     d_2=1.0_wp/a_m
     d_3=9.0_wp*a_h
     d_4=6.0_wp*a_m
     d_5=3.0_wp*(d_h+d_4)
     d_6=d_3+3.0_wp*d_4

     rim=1.0_wp/(1.0_wp+(d_m-d_4)/d_5) !1-Rf_c

     a_3=d_3/(d_2*d_m)
     a_5=d_5/(d_1*d_m)
     a_6=d_6/(d_2*d_m)

     sh_0=(b_h-d_4/d_m)/d_1 !stability-function for scalars  at neutr. strat.
     sm_0=(b_m-d_4/d_m)/d_2 !stability-function for momentum at neutr. strat.

     IF (tdc%lcpfluc) THEN !only if cp is not treated as a constant
       tur_rcpv=rcpv
       tur_rcpl=rcpl
     ELSE
       tur_rcpv=0.0_wp
       tur_rcpl=0.0_wp
     END IF

     IF (.NOT. lcuda_graph_turb_tran) THEN
       !$ACC WAIT(acc_async_queue)
     END IF
     !$ACC END DATA

END SUBROUTINE turb_setup

!==============================================================================
!==============================================================================

SUBROUTINE adjust_satur_equil ( tdc, khi, ktp, i1dim, &
!
   i_st, i_en, k_st, k_en,            &
!
   lcalrho, lcalepr,                  &
   lcaltdv,                           &
   lpotinp, ladjout,                  &
!
   icldmod,                           &
!
   zrcpv, zrcpl,                      &
!
   prs, t, qv, qc,                    &
   psf, fip,                          &
!
   rcld,                              &
!
   exner, dens, r_cpd,                &
!
   qst_t,   g_tet, g_h2o,             &
   tet_liq, q_h2o, q_liq,             &
!
   lacc, opt_acc_async_queue )

#ifdef _OPENACC
!Issue with Cray compiler (tested with 8.4.4), rutime error not present virt
!DIR$ INLINENEVER adjust_satur_equil
#endif

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

INTEGER, INTENT(IN) :: &
  i1dim,      & !length of blocks
  khi,        & !extra start index of vertical dimension for auxilary arrays (highest level of operations)
  ktp,        & !extra start index of vertical dimension for auxilary arrays (topping level of operations)
!
  i_st, i_en, & !horizontal start- and end-indices
  k_st, k_en    !vertical   start- and end-indices

  !Note(GPU):
  ! 'i1dim' is required here because of the CUDA graphs
  ! When this routine is called, 'i_en' is only set on the GPU. The CPU value is just "0" (or garbage)
  ! So we need to use 'i1dim' array size (basically just nproma) instead of "i_st:i_en"

LOGICAL, INTENT(IN) :: &
  lcalrho, &   !density calculation required
  lcalepr, &   !exner pressure calculation required
  lcaltdv, &   !calculation of special thermodynamic variables required
  lpotinp, &   !input temperature is a potential one
  ladjout      !output of adjusted variables insted of conserved ones

INTEGER, INTENT(IN) :: &
  icldmod      !mode of water cloud representation in transfer parametr.
               !-1: ignoring cloud water completely (pure dry scheme)
               ! 0: no clouds considered (all cloud water is evaporated)
               ! 1: only grid scale condensation possible
               ! 2: also sub grid (turbulent) condensation considered

REAL (KIND=wp), DIMENSION(:,ktp:), INTENT(IN) :: &
  prs,     &   !current pressure [Pa]
  t, qv        !current temperature [K] and water vapor content (spec. humidity) [1]

REAL (KIND=wp), DIMENSION(:,ktp:), OPTIONAL, INTENT(IN) :: &
  qc           !current cloud water content [1]

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
  psf          !surface pressure [Pa]

REAL (KIND=wp), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
  fip          !interpolation factor with respect to an auxilary level [1]

REAL (KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: &
  rcld         !inp: standard deviation of local super-saturation (SDSS) at mid-levels [1]
               !out: saturation fraction (cloud-cover) [1]

REAL (KIND=wp), DIMENSION(:,khi:), INTENT(INOUT) :: &
  exner        !current Exner-factor [1]

REAL (KIND=wp), DIMENSION(:,ktp:), TARGET, INTENT(INOUT), CONTIGUOUS :: &
  tet_liq, &   !inp: liquid water potent. temp. (only if 'fip' is present)      [K]
               !out: liquid water potent. temp. (or adjust. 't' , if "ladjout")
    q_h2o, &   !inp: total  water content       (only if 'fip' is present)      [1]
               !out: total  water content       (or adjust. 'qv', if "ladjout")
    q_liq      !out: liquid water content after adjustment                      [1]


REAL (KIND=wp), DIMENSION(:,khi:), OPTIONAL, INTENT(INOUT) ::  &
  dens,  &     !current air density [Kg/m3]
  r_cpd        !cp/cpd [1]

REAL (KIND=wp), DIMENSION(:,khi:), TARGET, OPTIONAL, INTENT(INOUT) :: &
               !Attention: "INTENT(OUT)" might cause some not-intended default-setting for not-treated levels!
  qst_t, &     !out: d_qsat/d_T [1/K]
               !aux:                                    TARGET for 'qvap', if ".NOT.ladjout"
  g_tet, &     !out: g/tet-factor of tet_l-gradient in buoyancy term [m/(s2K)]
               !aux: total  water content of ref. lev.; TARGET for 'temp', if ".NOT.ladjout"
  g_h2o        !out: g    -factpr of q_h2o-gradient in buoyancy term [m/s2]
               !aux: liq. wat. pot. temp. of ref. lev.; TARGET for 'virt'

REAL (KIND=wp), INTENT(IN) :: &
  zrcpv,  &    !0-rcpv  switch for cp_v/cp_d - 1 [1]
  zrcpl        !0-rcpl  switch for cp_l/cp_d - 1 [1]

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

LOGICAL :: lzacc
INTEGER :: acc_async_queue

REAL (KIND=wp) :: &
  pdry,  &     !partial pressure of dry air [Pa]
  mcor         !moist correction [1]

REAL (KIND=wp), DIMENSION(i1dim,k_st:k_en) :: &
  rprs         !reduced pressure [Pa]

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
  temp(:,:), &  !corrected temperature          [K]
  qvap(:,:), &  !corrected water vapour content [1]
  virt(:,:)     !reciprocal virtual factor      [1]

INTEGER :: &
  i, k, &
  icldtyp  !index for type of cloud diagnostics in SUB 'turb_cloud'
!-------------------------------------------------------------------------------

   !Note:
   !If 'qc' is not present, it is assumed that 't' and 'qv' already contain liquid water temperature
   ! and total water content.
   !In case of "k_st=khi=k_en", the vertical loop is only for one level.
   !Since a surface variable for 'qc' is not used in general, it is not forseen here as well, assuming
   ! a zero value at the surface. This implicats the surface variables to be identical with
   ! liquid water temperature and total water content (conserved variables for moist water conversions).
   !If 'fip' is resent, the conseved variables 'tet_liq' and q_h2o' at level 'k' are an interpolation between
   ! these values and those of the level "k-1" using the weight 'fip'. This is used  in order to interpolate
   ! rigid surface values at level "k=ke1" and those valid for the lowermost atmospheric full level "k=ke"
   ! to an atmopspheric surface level at the top of the roughness layer formed by land use. In this case it is
   ! always "k_st=ke1=k_en".
   !In this version, atmospheric ice is not included to the adjustment process.

   CALL set_acc_host_or_device(lzacc, lacc)

   IF(PRESENT(opt_acc_async_queue)) THEN
      acc_async_queue = opt_acc_async_queue
   ELSE
      acc_async_queue = 1
   ENDIF

   !$ACC DATA PRESENT(prs, t, qv, exner, rcld) &
   !$ACC   PRESENT(tet_liq, q_h2o, q_liq) &
   !$ACC   PRESENT(qc, dens, r_cpd) &
   !$ACC   PRESENT(qst_t, g_tet, g_h2o, fip) &
   !$ACC   CREATE(rprs) &
   !$ACC   COPYIN(i_en) &
   !$ACC   ASYNC(acc_async_queue) IF(lzacc)

   !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

   !Calculation of Exner-pressure:
   IF (lcalepr) THEN
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            exner(i,k)=zexner(prs(i,k))
         END DO
      END DO
   END IF

   !Conserved variables (with respect to phase change):
   IF (icldmod.EQ.-1 .OR. .NOT.PRESENT(qc)) THEN
     !$ACC LOOP SEQ
     DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
              q_h2o(i,k)=qv(i,k)
            tet_liq(i,k)= t(i,k)
         END DO
      END DO
   ELSE !water phase changes are possible and 'qc' is present
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
              q_h2o(i,k)=qv(i,k) +       qc(i,k) !tot. wat. cont.
            tet_liq(i,k)= t(i,k) - lhocp*qc(i,k) !liq. wat. temp.
         END DO
      END DO
   END IF

   !Transformation in real liquid water temperature:
   IF (lpotinp) THEN
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            tet_liq(i,k)=exner(i,k)*tet_liq(i,k)
         END DO
      END DO
   END IF

   !$ACC END PARALLEL

   !Interpolation of conserved variables (with respect to phase change)
   !onto the zero-level of interest at the top of the roughness layer:
   IF (PRESENT(fip)) THEN
      k=k_en !only for the lowest level
!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=i_st,i_en
           q_h2o(i,k)=             q_h2o(i,k-1)*(1.0_wp-fip(i))+  q_h2o(i,k)*fip(i)
         tet_liq(i,k)=exner(i,k)*tet_liq(i,k-1)*(1.0_wp-fip(i))+tet_liq(i,k)*fip(i)
      END DO
      !$ACC END PARALLEL
      !Note:
      !'tet_liq' at level "k-1" needs to be present and it is assumed to be already
      ! a potential (liquid water) temperature there, where it is still a pure
      ! liquid water temperature at level 'k'!
      !The roughness layer between the current level (at the rigid surface)
      ! and the desired level (at the top of the roughness layer) of index "k=ke1"
      ! is assumend to be a pure CFL layer without any mass. Consequentley
      ! the Exner pressure values at both levels are treated like being equal!
   END IF

   !Calculating effective saturated volume-fraction with respect to liquid cloud-water (cloud cover):

   IF (icldmod.EQ.0 .OR. (icldmod.EQ.-1 .AND. .NOT.PRESENT(qc))) THEN
      !Entire cloud-water evaporates or is ignored:

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            rcld(i,k)=0.0_wp !vanishing cloud-cover
           q_liq(i,k)=0.0_wp ! and liquid-water content
         END DO
      END DO
      !$ACC END PARALLEL

   ELSEIF (icldmod.EQ.-1) THEN
      !Clouds are present, but do not participate in turbulent phase-transitions:

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            rcld(i,k)=0.0_wp  !vanishing cloud-cover
           q_liq(i,k)=qc(i,k) ! and present liquid-water content
         END DO
      END DO
      !$ACC END PARALLEL
   ELSEIF (icldmod.EQ.1 .AND. PRESENT(qc)) THEN
      !Present cloud-water (due to pure grid-scale phase-transitions) is applied:

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
             rcld(i,k)=MERGE( 1.0_wp, 0.0_wp, (qc(i,k) .GT. 0.0_wp) ) !cloud-cover dependent on present liquid-water
            q_liq(i,k)=qc(i,k) !liquid-water content
         END DO
      END DO
      !$ACC END PARALLEL

   ELSE !a special cloud diagnostics needs to be employed

      IF (icldmod.EQ.1) THEN !only grid scale clouds possible
         icldtyp=0
      ELSE !sub grid scale clouds possible
         icldtyp=tdc%itype_wcld !use specified type of cloud diagnostics
      END IF

      CALL turb_cloud( tdc, i1dim=i1dim, ktp=ktp,              &
           istart=i_st, iend=i_en, kstart=k_st, kend=k_en,     &
           icldtyp=icldtyp,                                    &
           prs=prs, t=tet_liq, qv=q_h2o,                       &
           psf=psf,                                            &
           rcld=rcld,   & !inp: standard deviation of local super-saturation
                          !out: saturation fraction (cloud-cover)
           clwc=q_liq,  & !out: liquid water content
           lacc=lzacc, opt_acc_async_queue=opt_acc_async_queue )
   END IF

   !Calculating thermodynamic auxilary quantities:

   IF (ladjout) THEN !output of adjusted non conservative variables
      qvap =>   q_h2o
      temp => tet_liq
   ELSE !output of conserved variables
      qvap =>   qst_t
      temp =>   g_tet
   END IF

   IF (lcaltdv .OR. lcalrho) THEN
      virt => g_h2o
   END IF

   !$ACC DATA PRESENT(qvap, temp, virt) ASYNC(acc_async_queue) IF(lzacc)

   !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
   IF (ladjout .OR. lcaltdv .OR. lcalrho .OR. PRESENT(r_cpd)) THEN
      IF (.NOT.ladjout .AND. icldmod.LE.0) THEN !'temp' and 'vap' equal conserv. vars.
        !$ACC LOOP SEQ
         DO k=k_st, k_en
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=i_st,i_en
               temp(i,k)=tet_liq(i,k)
               qvap(i,k)=  q_h2o(i,k)
            END DO
         END DO
      ELSEIF (icldmod.GT.0) THEN !'temp' and 'qvap' my be different form conserv. vars.
         !$ACC LOOP SEQ
         DO k=k_st, k_en
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=i_st,i_en
               temp(i,k)=tet_liq(i,k)+lhocp*q_liq(i,k) !corrected temperature
               qvap(i,k)=  q_h2o(i,k)-      q_liq(i,k) !corrected water vapor
            END DO
         END DO
      END IF
      !Note: In the remaining case "ladjout .AND. icldmod.LE.0" 'temp' and 'qvap'
      !      already point to the conserved variables.
   END IF

   IF (lcaltdv .OR. lcalrho) THEN
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            virt(i,k)=1.0_wp/(1.0_wp+rvd_m_o*qvap(i,k)-q_liq(i,k)) !rezipr. virtual factor
            rprs(i,k)=virt(i,k)*prs(i,k)                   !reduced pressure profile
         END DO
      END DO
   END IF

   IF (lcalrho .AND. PRESENT(dens)) THEN
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            dens(i,k)=rprs(i,k)/(r_d*temp(i,k))
         END DO
      END DO
   END IF

   IF (PRESENT(r_cpd)) THEN
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            IF (tdc%lcpfluc) THEN
               r_cpd(i,k)=1.0_wp+zrcpv*qvap(i,k)+zrcpl*q_liq(i,k) !Cp/Cpd
            ELSE
               r_cpd(i,k)=1.0_wp
            END IF
         END DO
      END DO
   END IF

   IF (.NOT.ladjout) THEN !the potential (liquid water) temperature values are requested
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st,i_en
            tet_liq(i,k)=tet_liq(i,k)/exner(i,k) !liquid water pot. temp.
         END DO
      END DO
   END IF

   IF (lcaltdv) THEN
      !$ACC LOOP SEQ
      DO k=k_st, k_en
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(pdry)
         DO i=i_st,i_en
            pdry=(1.0_wp-qvap(i,k))*rprs(i,k)                !partial pressure of dry air
            qst_t(i,k)=zdqsdt( temp(i,k), zqvap( zpsat_w( temp(i,k) ), pdry ) )
                                                             !d_qsat/d_T (new version)
         END DO
      END DO

      IF (icldmod.EQ.-1) THEN !no consideration of water phase changes
         !$ACC LOOP SEQ
         DO k=k_st, k_en
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=i_st,i_en
               g_h2o(i,k)=grav*(rvd_m_o*virt(i,k))              !g    -factor of q_h2o-gradient
               g_tet(i,k)=grav*(exner(i,k)/temp(i,k))           !g/tet-factor of tet_l-gradient
            END DO
         END DO
      ELSE !water phase changes are possible
         !$ACC LOOP SEQ
         DO k=k_st, k_en
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mcor)
            DO i=i_st,i_en
               rcld(i,k)=tdc%c_scld*rcld(i,k)/(1.0_wp+rcld(i,k)*(tdc%c_scld-1.0_wp)) !effective cloud-cover
               mcor=rcld(i,k)*(lhocp/temp(i,k)-(1.0_wp+rvd_m_o)*virt(i,k)) &
                        /(1.0_wp+qst_t(i,k)*lhocp)              !moist correction by turb. phase-transitions

               g_h2o(i,k)=grav*(rvd_m_o*virt(i,k)+mcor)         !g    -factor of q_h2o-gradient

               g_tet(i,k)=grav*(exner(i,k)/temp(i,k) &          !g/tet-factor of tet_l-gradient
                               -mcor*exner(i,k)*qst_t(i,k))
            END DO
         END DO
      END IF

   END IF
   !$ACC END PARALLEL
   IF (.NOT. lcuda_graph_turb_tran) THEN
      !$ACC WAIT(acc_async_queue)
   END IF
   !$ACC END DATA
   !$ACC END DATA

   !Die thermodynamischen Hilfsgroessen wurden hier unter Beruecksichtigung der diagnostizierten
   !Kondensationskorrektur gebildet, indem die entspr. korrigierten Werte fuer t, qv und qc benutzt wurden.
   !Beachte, dass es zu einer gewissen Inkosistenz kommt, wenn die Dichte nicht angepasst wird!

END SUBROUTINE adjust_satur_equil

!==============================================================================
!==============================================================================

SUBROUTINE solve_turb_budgets ( tdc, it_s, it_start, &
!
   i1dim, i_st, i_en,                &
   khi, ktp, kcm, k_st, k_en, k_sf,  &
!
   ntur, nvor,                       &
!
   lssintact, lupfrclim, lpres_edr,  &
!
   ltkeinp,                          &
!
   imode_stke, imode_vel_min,        &
!
   dt_tke, fr_tke,                   &
!
   fm2, fh2, ft2, lsm, lsh, tls,     &
!
   fcd, tvt, avt, velmin,            &
   tke, ediss,                       &
!
   lactcnv, laddcnv,                 &
!
   exner, r_cpd, qst_t,              &
!
   rcld, & !inp: saturation fraction (cloud-cover)
           !out: std. deviat. of local super-saturat.
!
   lcircterm,                        &
   dens,                             &
   l_pat, l_hori,                    &
!
   grd,                              &
!
   lacc, opt_acc_async_queue         )

!------------------------------------------------------------------------------

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

INTEGER, INTENT(IN) :: &  !
!
  it_s,                   & !iteration step
  it_start,               & ! start index of iteration
!
  i1dim,                  & !length of blocks
  i_st, i_en,             & !horizontal start- and end-indices
  khi,                    & !extra start index of vertical dimension for auxilary arrays (highest level of operations)
  ktp,                    & !extra start index of vertical dimension for auxilary arrays (topping level of operations)
  k_st, k_en,             & !vertical   start- and end-indices
  k_sf,                   & !vertical index for the surface level (top of grid-scale roughness layer)
  kcm,                    & !level index of the upper canopy bound
!
  ntur,                   & !current new time step index of tke
  nvor,                   & !current     time step index of tke
!
  imode_stke,             & !mode of solving the TKE equation
  imode_vel_min             !mode of determination of minimum turbulent velocity scale

LOGICAL, INTENT(IN)  :: &
!
  lssintact, & !seperate treatment of non-turbulent shear (by scale interaction) requested
  lupfrclim, & !enabling an upper limit for TKE-forcing
  lpres_edr, & !if edr is present in calling routine
  ltkeinp      !TKE present as input for current time level 'ntur'

REAL (KIND=wp), INTENT(IN) :: &
!
  dt_tke,    & !time step for the 2-nd order porgnostic variable 'tke'
  fr_tke       !1.0_wp/dt_tke

REAL (KIND=wp), DIMENSION(:,khi:), TARGET, INTENT(IN) :: &
!
  fm2,  &  !squared frequency of mechanical forcing   [1/s2]
  fh2,  &  !squared frequency of thermal    forcing   [1/s2]
  ft2      !squared frequency of pure turbulent shear [1/s2]

REAL (KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: &
!
  lsm, &   !turbulent length scale times value of stability function for momentum       [m]
  lsh      !turbulent length scale times value of stability function for scalars (heat) [m]

REAL (KIND=wp), DIMENSION(:,khi:), INTENT(INOUT) :: &
!
  tls      !turbulent master scale [m]

REAL (KIND=wp), DIMENSION(:,kcm:), OPTIONAL, INTENT(IN) :: &
!
  fcd      !frequency of small scale canopy drag [1/s]

REAL (KIND=wp), DIMENSION(:,khi:), TARGET, INTENT(IN) :: &
!
  tvt      !turbulent transport of turbulent velocity scale [m/s2]

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
  avt      !advective transport of turbulent velocity scale [m/s2]

REAL (KIND=wp), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
!
  velmin   !location-dependent minimum velocity [m/s]

REAL (KIND=wp), DIMENSION(:,:,:), TARGET, INTENT(INOUT) :: &
!
  tke      !q:=SQRT(2*TKE); TKE='turbul. kin. energy'       [ m/s ]
           ! (defined on half levels)

REAL (KIND=wp), DIMENSION(:,ktp:), INTENT(INOUT) :: &
           !Attention: "INTENT(OUT)" might cause some not-intended default-setting for not-treated levels!
!
  ediss    !eddy dissipation rate                           [m2/s3]

LOGICAL, INTENT(IN)  :: &
!
  lactcnv, & !flux conversion is activated at all (first of all, because water clouds are active)
  laddcnv    !fluxes of conserved variables need to be converted into those of non-conserved ones
             ! apart from an explicit moisture correction, which is being performed at "lexpcor=T"

REAL (KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: &
!
  rcld     !inp: saturation fraction [1]
           !out: standard deviation of local super-saturation (SDSS) [1]

REAL (KIND=wp), DIMENSION(:,khi:), INTENT(IN) :: &
!
  exner, & !Exner-pressure [1]
  r_cpd    !Cp/Cpd         [1]

REAL (KIND=wp), DIMENSION(:,khi:), INTENT(IN) :: &
  qst_t    !d_qsat/d_T [1/K]

LOGICAL, INTENT(IN)  :: &
!
  lcircterm  !calculation of the old "circulation term" required

REAL (KIND=wp), DIMENSION(:,khi:), INTENT(IN) :: &
!
  dens     !air density at half levels [Kg/m3]

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
  l_pat, & !effective length scale of near-surface circulation patterns [m]
           ! (scaling the near-surface circulation acceleration)
  l_hori   !horizontal grid spacing                                     [m]

REAL (KIND=wp), DIMENSION(:,ktp:,0:), INTENT(INOUT) :: &
!
  grd !inp (3-rd dim. > 0): gradients of 1-st order model-variables including those of conserved variables
      !                      for temperature and moisture needed for SUB 'turb_stat' during a SCM-run
      !                      and  for calculation of SDSS
      !    (3-rd dim. = 0): half-level pressure, providing (if devided by density and a length of coherence)
      !                      a related circulation acceleration
      !out (3-rd dim. > 0): effectiv vertical gradients (for Tet, qv, qc) after linear flux-conversion
      !    (3-rd dim. = 0): effective vertical gradient of available Circulation Kinetic Energy (per mass)
      !                      due to near surface thermal inhomogeneity (CKE), which is a related acceleration
      !                      being used to express the old "circulation term"

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

LOGICAL :: lzacc
INTEGER :: acc_async_queue

!------------------------------------------------------------------------------

! Local variables

INTEGER :: &
!
  i, k, k_tvs        !running indices for horizontal location loops or vertical levels

REAL (KIND=wp) :: &
!
  val1, val2, wert, fakt, & !auxilary values
  w1, w2, &                 !weighting factors
  q1, q2, q3, &             !turbulent velocity scales [m/s]
!
  a11, a12, a22, a21, &
  a3, a5, a6, be1, be2, r_b_m, &
  bb1, bb2, &
  a_1, a_2, &
!
  sm,sh, & !dimensionless stability functions
  det,   & !determinant of 2X2 matrix
!
  gm,gh, & !dimensionless forcing due to shear and buoyancy
  gama,  & !parameter of TKE-equilibrium
  gam0,  & !parameter of TKE-equilibrium (maximal value)
  tim2     !square of turbulent time scale [s2]

REAL (KIND=wp) :: &
!
  flw_h2o_g, & !weight of h2o_g-flux
  flw_tet_l    !weight of tet_l-flux

REAL (KIND=wp) :: &
!
  grav2, & !"grav**2"
  l_coh    !length-scale of coherence directly employed for expressing the near surface circulation acceleration

REAL (KIND=wp), DIMENSION(i1dim) :: &
!
  l_dis, & !dissipation length scale [m]
  l_frc, & !forcing     length scale [m]
  frc,   & !effective TKE-forcing (acceleration of turbulent motion) [m/s2]
  tvs_m    !minimal turbulent velocity scale [m/s]

REAL (KIND=wp), DIMENSION(i1dim,1), TARGET :: &
!
! tvs      !turbulent velocity scale          [m/s]
  tvs_u    !updated turbulent velocity scale  [m/s]

REAL (KIND=wp), DIMENSION(i1dim,-ntyp:10), TARGET :: &
!
  dd       !local derived turbulence parameter

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
!
  tvs_0(:,:),   & ! pointer for intermediate turbulent velocity scale [m/s]
  fm2_e(:,:)      ! pointer for the effective mechanical forcing [1/s2]

LOGICAL :: lpres_avt, lpres_fcd, lrogh_lay, alt_gama, lcorr, lstbsecu

INTEGER :: imode_stbcorr !mode of correcting the stability function (=ABS(imode_stbcalc)
!-------------------------------------------------------------------------------

  LOGICAL, PARAMETER :: lstfnct=.TRUE. !calculate stability functions
                                       !(otherwise the former values remain unchainged)
  lpres_avt=PRESENT(avt) !array for advection-tendency of TKE is present
  lpres_fcd=PRESENT(fcd) !array for small-scale canpy drag is present

  imode_stbcorr=ABS(tdc%imode_stbcalc) !mode of correcting the stability function
  alt_gama=(imode_stbcorr.EQ.1 .AND. .NOT.ltkeinp) !alternative gama-Berechnung

  IF (lssintact) THEN !seperate treatment of shear by scale interaction
     fm2_e => ft2 !effective shear is pure turbulent shear
  ELSE
     fm2_e => fm2 !effective shear is total mechanical shear
  END IF

  CALL set_acc_host_or_device(lzacc, lacc)

  IF(PRESENT(opt_acc_async_queue)) THEN
     acc_async_queue = opt_acc_async_queue
  ELSE
     acc_async_queue = 1
  ENDIF

  !Note: 'fm2_e' is being used in the expressions for thermal stability.

  ! Pointer assignment was moved out of the k loop to improve GPU performance
  IF (lpres_avt) THEN !explizite Addition der Advektions-Inkremente
     tvs_0 => tvs_u(:,:) !die effektiven TKE-Vorgaengerwerte werden auf 'tvs_u(:,1)' geschrieben
                         ! und anschliessend durch die neuen Werte ueberschrieben
     k_tvs = 1 !always used as k-index of 'tvs_0' within the k-loop below
  ELSE !benutze die Vorgaengerwerte
     tvs_0 => tke(:,:,nvor)
  END IF

  r_b_m=1.0_wp/b_m

  grav2=grav**2

  a_1=a_h; a_2=a_m

  IF (ltkeinp) THEN
     gam0=1.0_wp/d_m !TKE-equilibrium
  ELSE
     !Obere Schranke fuer die Abweichung 'gama' vom TKE-Gleichgewicht:
     gam0=tdc%stbsecu/d_m+(1.0_wp-tdc%stbsecu)*b_m/d_4
  END IF

#ifdef TST_CODE_VERS
  ! Precalculations for the case ".NOT.(alt_gama .OR. lrogh_lay)":
  a3=d_3*gam0*a_2; a5=d_5*gam0*a_1; a6=d_6*gam0*a_2
  wert=d_4*gam0; bb1=(b_h-wert)*a_1; bb2=(b_m-wert)*a_2
#endif

  !$ACC DATA PRESENT(fm2, fh2, ft2, lsm, lsh, tls) &
  !$ACC   PRESENT(fcd, tvt, avt, velmin, tke, ediss) &
  !$ACC   PRESENT(fm2_e, tvs_0) &
  !$ACC   CREATE(l_dis, l_frc, frc, tvs_m, tvs_u, dd) &
  !$ACC   COPYIN(i_en) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  ! Correction of turbulent master lenght-scale due to stable thermal stratification:
  IF (tdc%a_stab.GT.0.0_wp .AND. it_s==it_start) THEN
     !$ACC LOOP SEQ
     DO k=k_st,k_en !von oben nach unten
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
        DO i=i_st, i_en
           wert=tdc%a_stab*SQRT(MAX( 0.0_wp, fh2(i,k)) )
           tls(i,k)=tke(i,k,nvor)*tls(i,k)/(tke(i,k,nvor)+wert*tls(i,k))
        END DO
     END DO

     !Attention:
     !Even if "a_stab=0." has been set initially (in 'mo_turbdiff_config' or by 'turbdiff_nml'),
     ! this block may be executed due to PERTURBATIONS applied to 'a_stab'!
  END IF

  ! Pre-filling the fields for (derived) turbulence parameters with values valid for the free atmosphere:
!DIR$ IVDEP
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i=i_st, i_en
     dd(i,0)=d_m
     dd(i,1)=d_1; dd(i,2)=d_2; dd(i,3)=d_3
     dd(i,4)=d_4; dd(i,5)=d_5; dd(i,6)=d_6
     dd(i,7)=rim
     dd(i,8)=gam0
     dd(i,9)=a_1; dd(i,10)=a_2
  END DO

  IF (PRESENT(velmin) .AND. imode_vel_min.EQ.2) THEN !nutze variierendes 'tvs_m'
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
        tvs_m(i)=tdc%tkesecu*velmin(i) !effektiver variiernder Minimalwert fuer 'tvs_u'
     END DO
  ELSE
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
        tvs_m(i)=tdc%tkesecu*tdc%vel_min !effektiver konstanter Minimalwert fuer 'tvs_u'
     END DO
  END IF

  !$ACC END PARALLEL

  !XL_GPU_OPT : need to make k and i purely nested
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP SEQ PRIVATE(lrogh_lay, lstbsecu, w1, w2)
  DO k=k_st, k_en !ueber alle Schichten beginnend mit der freien Atm.

     lrogh_lay=(lporous .AND. k.GE.kcm .AND. lpres_fcd) !innerhalb der Rauhigkeitsschicht

     IF (lrogh_lay) THEN !innerhalb der Rauhigkeitsschicht

        ! Calculating roughness-layer corrections, except "volume terms" (treated as a part of diffusion):

!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
        DO i=i_st, i_en
!          Berechnung der modifizierten Modellparameter:

           wert=3.0_wp*tls(i,k)*fcd(i,k)/tke(i,k,nvor)

           dd(i,0)=d_m/(1.0_wp+d_m*wert)
           dd(i,1)=a_h/(1.0_wp+a_h*wert)
           dd(i,2)=a_m/(1.0_wp+2.0_wp*a_m*wert)
           dd(i,3)=9.0_wp*dd(i,1)
           dd(i,4)=6.0_wp*dd(i,2)
           dd(i,5)=3.0_wp*(d_h+dd(i,4))
           dd(i,6)=dd(i,3)+3.0_wp*dd(i,4)
           dd(i,1)=1.0_wp/dd(i,1)
           dd(i,2)=1.0_wp/dd(i,2)
           dd(i,7)=1.0_wp/(1.0_wp+(dd(i,0)-dd(i,4))/dd(i,5))
            IF (ltkeinp) THEN
               dd(i,8)=1.0_wp/dd(i,0) !TKE-equilibrium
            ELSE
              !Obere Schranke fuer die Abweichung 'gama' vom TKE-Gleichgewicht:
              dd(i,8)=tdc%stbsecu/dd(i,0)+(1.0_wp-tdc%stbsecu)*b_m/dd(i,4)
            END IF
            dd(i,9)=1.0_wp/dd(i,1); dd(i,10)=1.0_wp/dd(i,2)
        END DO
     END IF

     ! Total TKE-forcing in [m/s2]:

!DIR$ IVDEP
     !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(val1, val2)
     DO i=i_st, i_en
        l_dis(i)=tls(i,k)*dd(i,0)          !length scale of dissipation
        l_frc(i)=tls(i,k)*dd(i,4)*r_b_m    !length scale of forcing

        val1=lsm(i,k)*fm2_e(i,k); val2=lsh(i,k)*fh2(i,k)
!_______________________________________________________________
!test: Keine Beschraenkung von frc nach unten
        frc(i)=MAX( val1-val2, tdc%frcsecu*dd(i,7)*val1 )
! frc(i)=val1-val2
!_______________________________________________________________
     END DO

     !Beachte:
     !Bei "frcsecu=1" wird 'frc' so nach unten beschraenkt, dass die krit. Rf-Zahl (Rf=1-'dd(i,7)')
     ! im TKE-Quellgleichgewicht (ohne Transport und Zeittendenz) nicht ueberschritten wird.
     !Bei "frcsecu<1" wird 'frc' entsrechend weniger eingeschraenkt.
     !Bei "frcsecu=0" wird Rf nicht > 1 und die Summe der TKE-Quellterme wird nicht negativ.
     !Bei "frcsecu<0" kann "Rf>1" und die Summe der TKE-Quellterme auch negativ werden.

!________________________________________________________________
!test: keine Beschraenkung von frc nach oben (Auskommentierung)
!Achtung: Korrektur: Ermoeglicht obere 'frc'-Schranke im Transferschema (wie in COSMO-Version)
!    IF (tdc%frcsecu.GT.0) THEN
     IF (tdc%frcsecu.GT.0 .AND. lupfrclim) THEN
       !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=i_st, i_en
           frc(i)=MIN( frc(i), tdc%frcsecu*tke(i,k,nvor)**2/l_frc(i)+(1.0_wp-tdc%frcsecu)*frc(i) )
        END DO
     END IF
!________________________________________________________________
     !Beachte:
     !'tke(i,k,nvor)**2/l_frc(i)' ist der Maximal-Wert fuer 'frc', mit dem die Abweichung
     ! vom TKE-Gleichgewicht den Wert besitzt, der mit der gegebenen 'tke' bei neutraler Schichtung
     ! nicht ueberschritten werden kann.
     !Bei "frcsecu<=0" wird 'frc' nicht nach oben beschraenkt.

     IF (lssintact) THEN !shear forcing by scale interaction needs to be added in order to act as a
                         ! general positive-definete source-term for TKE
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=i_st, i_en
           frc(i)=frc(i)+lsm(i,k)*(fm2(i,k)-ft2(i,k)) !complete forcing
        END DO
     END IF

     ! Calculating updated turbulent velocity scale, which is SQRT(2*TKE):

     IF (.NOT.ltkeinp) THEN

        !Bestimmung der Zwischenwerte:

        IF (lpres_avt) THEN !explizite Addition der Advektions-Inkremente
           !Notice that "k_tvs=1" has been set for this case before; and it will not be changed during this k-loop!
!DIR$ IVDEP
          !$ACC LOOP GANG(STATIC: 1) VECTOR
           DO i=i_st, i_en
              tvs_0(i,k_tvs)=MAX( tvs_m(i), tke(i,k,nvor)+avt(i,k)*dt_tke )
           END DO

           !Beachte:
           !Die Advektions-Inkremente werden explizit in 'tvs_0' aufgenommen,
           ! damit die spaetere zeitliche Glaettung (tkesmot) nicht die Transportgeschindigkeit
           ! der Advektion beeinflusst (verlangsamt).
           !Da die Advektions-Inkremente auch von benachbarten 'tke'-Werten in horiz. Richtung abhaengen,
           ! werden sie auch nicht in die optionale Iteration gegen einen Gleichgewichtswert einbezogen
           ! und werden somit nur beim letzten Iterations-Schritt addiert.
        ELSE !benutze die Vorgaengerwerte
           k_tvs = k !k-index for tvs_0 only
        END IF

        !Integration von SQRT(2TKE):
        IF (imode_stke.EQ.1) THEN !1-st (former) type of prognostic solution
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q1, q2)
           DO i=i_st, i_en
              q1=l_dis(i)*fr_tke
              q2=MAX( 0.0_wp, tvs_0(i,k_tvs)+tvt(i,k)*dt_tke )+frc(i)*dt_tke
              tvs_u(i,1)=q1*(SQRT(1.0_wp+4.0_wp*q2/q1)-1.0_wp)*0.5_wp
           END DO
        ELSEIF (imode_stke.GE.2) THEN !2-nd (new) type of prognostic solution
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q1, fakt)
           DO i=i_st, i_en
              fakt=1.0_wp/(1.0_wp+2.0_wp*dt_tke*tvs_0(i,k_tvs)/l_dis(i))
              q1=fakt*(tvt(i,k)+frc(i))*dt_tke
              tvs_u(i,1)=q1+SQRT( q1**2+fakt*( tvs_0(i,k_tvs)**2+2.0_wp*dt_tke*con_m*fm2(i,k) ) )
           END DO
        ELSEIF (imode_stke.EQ.-1) THEN !diagn. solution of station. TKE-equation
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q2, q3)
           DO i=i_st, i_en
              q2=l_dis(i)*(tvt(i,k)+frc(i))
              q3=l_dis(i)*con_m*fm2(i,k)
              IF (q2.GE.0.0_wp) THEN
               ! tvs_u(i,1)=SQRT(q2+q3/tvs_0(i,k_tvs))
                 tvs_u(i,1)=EXP( z1d3*LOG( q2*tvs_0(i,k_tvs)+q3 ) )
              ELSEIF (q2.LT.0.0_wp) THEN
                 tvs_u(i,1)=SQRT( q3/( tvs_0(i,k_tvs)-q2/tvs_0(i,k_tvs) ))
              ELSE
                 tvs_u(i,1)=EXP( z1d3*LOG(q3) )
              END IF
           END DO
        ELSE !standard diagnostic solution
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR
           DO i=i_st, i_en
              tvs_u(i,1)=SQRT( l_dis(i)*MAX( frc(i), 0.0_wp ) )
           END DO
        END IF

        w1=tdc%tkesmot; w2=1.0_wp-tdc%tkesmot

!DIR$ IVDEP
       !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q2)
       DO i=i_st, i_en
          q2=SQRT( l_frc(i)*MAX( frc(i), 0.0_wp ) )
!________________________________________________________________
!test< ohne tkesecu*vel_min als untere Schranke
          tke(i,k,ntur)=MAX( tvs_m(i), tdc%tkesecu*q2, w1*tvs_0(i,k_tvs)+w2*tvs_u(i,1) )
!tke(i,k,ntur)=MAX(           tdc%tkesecu*q2, w1*tke(i,k,nvor) +w2*tvs_u(i,1) )
!test>
!________________________________________________________________

          !'q2' ist ein Minimalwert fuer 'tke', mit dem die Abweichung vom TKE-Gleichgewicht
          ! den Wert besitzt, der mit dem gegebenen 'frc' bei neutraler Schichtung nicht
          ! ueberschritten werden kann.
          !'tke(:,:,ntur)' ist der Wert, der im naechsten Prognoseschritt benutzt wird!
       END DO

     END IF ! (.NOT.ltkeinp)

     ! Calculating updated value for stability lenght (stability-factor times turbulent master length-scale):

     IF (lstfnct) THEN !calculation of stability function required

        lstbsecu=(tdc%imode_stbcalc.LT.0) !apply preconditioning

!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR &
        !$ACC   PRIVATE(d_1, d_2, d_3, d_4, d_5, d_6, a_1, a_2) &
        !$ACC   PRIVATE(gam0, gama, lcorr, sm, sh, det, tim2, gh, gm) &
        !$ACC   PRIVATE(be1, be2, a11, a12, a21, a22, wert) &
        !$ACC   PRIVATE(bb1, bb2, a3, a5, a6, val1, val2, fakt)
        DO i=i_st, i_en

           IF (lrogh_lay .OR. lzacc) THEN
              d_1=dd(i,1); d_2=dd(i,2); d_3=dd(i,3)
              d_4=dd(i,4); d_5=dd(i,5); d_6=dd(i,6)
              gam0=dd(i,8)

              a_1=dd(i,9); a_2=dd(i,10)
           END IF

           tim2=(tls(i,k)/tke(i,k,ntur))**2 !Quadrat der turbulenten Zeitskala

!-------------------------------------------------
           sh=1.0_wp; sm=1.0_wp
!-------------------------------------------------
           !Note:
           !Although 'sh/m' are being evaluated for each of the follwoing cases,
           ! default-values needs to be set for vectorization!

           lcorr=.TRUE.

           IF (imode_stbcorr.EQ.2 .OR. fh2(i,k).GE.0.0_wp) THEN !standard solution at given 'tke'

              ! Solution based on given 'tke':

              !Dimensionslose Antriebe der Turbulenz
              gh=fh2  (i,k)*tim2 !durch thermischen Auftrieb
              gm=fm2_e(i,k)*tim2 !durch mechanische Scherung

              IF (.NOT.lstbsecu) THEN !Koeffizientenbelegung ohne Praeconditionierung

                 be1=b_h
                 be2=b_m

                 a11=d_1+(d_5-d_4)*gh
                 a12=d_4*gm
                 a21=(d_6-d_4)*gh
                 a22=d_2+d_3*gh+d_4*gm

              ELSE !Koeffizientenbelegung mit Praeconditionierung:

                 wert=1.0_wp/tim2

                 a11=d_1*wert+(d_5-d_4)*fh2(i,k)
                 a12=d_4*fm2_e(i,k)
                 a21=(d_6-d_4)*fh2(i,k)
                 a22=d_2*wert+d_3*fh2(i,k)+d_4*fm2_e(i,k)

                 be1=1.0_wp/MAX( ABS(a11), ABS(a12) )
                 be2=1.0_wp/MAX( ABS(a21), ABS(a22) )

                 a11=be1*a11; a12=be1*a12
                 a21=be2*a21; a22=be2*a22

                 be1=be1*wert*b_h; be2=be2*wert*b_m

              END IF

              det=a11*a22-a12*a21
              sh=be1*a22-be2*a12
              sm=be2*a11-be1*a21

              IF (det.GT.0.0_wp .AND. sh.GT.0.0_wp .AND. sm.GT.0.0_wp) THEN !Loesung moeglich
                 !solution possible, which always holds at "fh2>=0":
                 det=1.0_wp/det
                 sh=sh*det
                 sm=sm*det

                 lcorr=MERGE( (sm*gm-sh*gh.GT.gam0), .FALSE., imode_stbcorr.EQ.2 )
              END IF

           END IF !standard solution at given 'tke'

           !Note:
           !The pure solution for the stability functions 'sh' and 'sm' through the above linear system, inserting
           ! given 'tke'-values from the just before solved TKE-equation, may become non-realizable at strongly
           ! unstable stratification, where this situation is connected with an infinite positive deviation
           ! 'gama:=sm*gm-sh*gh' from TKE-equilibrium "gama=dd(i,0)=:d_m".
           !This problem is circumvented by employing a modified solution based on a predescribed deviation 'gama',
           ! which is expressed by 'frc*tim2/tls' and an upper limit 'gam0'.
           !This modified solution is being executed at "lcorr=T" dependent on 'imode_stbcorr':
           ! "imode_stbcorr=1": Just for strictly non-stable stratification, that means at "fh2<0".
           ! "imode_stbcorr=2": Only, if 'gama' exceeds the threshold 'gam0' or if the standard solution is not possible.
           !  At "stbsecu=0", this can only happen for "fh2<0", while it may happen also for "fh2>=0" at "stbsecu>0".

           ! Correction with restricted 'gama':

           IF (lcorr) THEN !Korrektur noetig
              !Loesung bei der die TKE-Gleichung in Gleichgewichtsform eingesetzt ist,
              !wobei 'gama' die Abweichung vom Gleichgewicht darstellt,
              !die so beschraenkt wird, dass immer eine Loesung moeglich ist:

#ifdef TST_CODE_VERS
              !Attention:
              !With the below "IF"-statement, NEC does not vectorize this loop!
              IF (alt_gama .OR. lrogh_lay) THEN !below parameters are not constant
#endif
                 gama=MERGE( MIN( gam0, frc(i)*tim2/tls(i,k) ), gam0, alt_gama )
                 !Note:
                 !At "fh2<0", "gama<=gama0" always secures a realizable solution for any "0<stbsecu<=1".
                 !At "fh2>=0" (which implies "imode_stbcorr=2" and thus "alt_gama=F"), it is "0<gama=gam0<sm*gm-sh*gh",
                 ! which always secures a realizable solution due to the properties of governing relations.
                 !For "0<fh2->0" and "frcsecu=0", 'gam0' can get arbitrarily close to the critical value,
                 ! which may even be hit or exceeded due to rounding errors!
                 !Hence, in order to surely avoid non-positive 'sm'- or 'sh-values, "stbsecu>0" should always be applied.
                 wert=d_4*gama

                 bb1=(b_h-wert)*a_1
                 bb2=(b_m-wert)*a_2

                 a3=d_3*gama*a_2
                 a5=d_5*gama*a_1
                 a6=d_6*gama*a_2
#ifdef TST_CODE_VERS
              END IF
#endif

              val1=(fm2_e(i,k)*bb2+(a5-a3+bb1)*fh2(i,k))/(2.0_wp*bb1)
              val2=val1+SQRT(val1**2-(a6+bb2)*fh2(i,k)*fm2_e(i,k)/bb1)

              fakt=fh2(i,k)/(val2-fh2(i,k))
              sh=bb1-a5*fakt
              sm=sh*(bb2-a6*fakt)/(bb1-(a5-a3)*fakt)
           END IF !lcorr

           lsh(i,k)=tls(i,k)*sh
           lsm(i,k)=tls(i,k)*sm

        END DO !i=i_st, i_en

     END IF !(lstfnct)

!------------------------------------------------------------------------------------
#ifdef SCLM
     IF (lsclm .AND. it_s.EQ.tdc%it_end) THEN
        CALL turb_stat(k=k, &
                       lsm=  lsm(imb,k), lsh=lsh(imb,k),      &
                       fm2=fm2_e(imb,k), fh2=fh2(imb,k),      &
                       tls=  tls(imb,k), tvs=tke(imb,k,ntur), &
                       tvt=  tvt(imb,k), grd=grd(imb,k,:),    &
                       d_m=dd(imb,0), d_h=d_h, a_m=dd(imb,10))
     END IF
#endif
!SCLM--------------------------------------------------------------------------------

     IF ((lpres_edr .OR. tdc%ltmpcor).AND.lrogh_lay) THEN
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=i_st, i_en
           ediss(i,k)=dd(i,0) !save location dependent dissipation parameter
        END DO
     END IF

  END DO !k
  !$ACC END PARALLEL

  !$ACC DATA PRESENT(rcld) &
  !$ACC   PRESENT(exner, r_cpd, qst_t, dens, l_pat, l_hori, grd) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)

  IF (it_s.EQ.tdc%it_end) THEN !only for the last iteration step

     ! Calculating vertical acceleration (CKE-gradient) used for the raw "circulation term":

     IF (lcircterm) THEN !calculation of the old "circulation term" required

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lacc)
        !$ACC LOOP SEQ
        DO k=k_st, k_sf !for all boundary levels including the extra lowermost boundary
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt, l_coh)
           DO i=i_st, i_en
              fakt=z1-z2*ABS(rcld(i,k)-z1d2) !coherence factor due to cloud-cover
              l_coh=MAX( l_pat(i), SQRT(fakt*tls(i,k)*l_hori(i)) ) !combintion of pure land-use pattern scale
                                                                   ! with acontribution by fractional cloud-cover
              fakt=fh2(i,k)*grd(i,k,0)/(dens(i,k)*grav2)     != Rd/g*exnr*grad(tet_v) (dimensionless Tet_v-gradient)
              l_coh=l_coh*SIGN(z1,fakt)*MIN( ABS(fakt), z1 ) !vert. coeherence-scale of near-surace circ. patterns

              grd(i,k,0)=l_coh*fh2(i,k) !resulting circulation acceleration in [m/s2]
           END DO
        END DO
        !$ACC END PARALLEL
     END IF

     ! Calculating effective vertical gradients of liquid-water content:

     IF (lactcnv .AND. tdc%lexpcor) THEN !flux conversion with an adjusted liquid-water flux required
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lacc)
        !$ACC LOOP SEQ
        DO k=k_st, k_sf !for all boundary levels including the extra lowermost boundary
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(flw_h2o_g, flw_tet_l)
           DO i=i_st, i_en
              flw_h2o_g=rcld(i,k)/(z1+lhocp*qst_t(i,k))    !weight of h2o_g-flux dependent on cl-cov.
!             flw_tet_l=-flw_h2o_g*epr_qst_t               !weight of tet_l-flux
              flw_tet_l=-flw_h2o_g*(exner(i,k)*qst_t(i,k)) !weight of tet_l-flux

              ! Effective vertical gradient of liquid water content:
              grd(i,k,liq) = flw_tet_l*grd(i,k,tet_l) & !eff_grad(liq)
                           + flw_h2o_g*grd(i,k,h2o_g)
           END DO
        END DO
        !$ACC END PARALLEL
        !Note: Otherwise, 'grd(:,:,liq)' remains unchanged.
     END IF

     ! Calculating Standard Deviation of local Super-Saturation (SDSS):

     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lacc)
     !$ACC LOOP SEQ
     DO k=k_st, k_en !for all boundary levels included for the turbulence model (correct coding)
        !Note:
        !If SUB 'solve_turb_budgets' is called by SUB 'turbdiff' with "k_en=ke<k_sf", 'rcld' (SDSS) at level 'k_sf'
        ! has already been calculated through the call of SUB 'solve_turb_budgets' by SUB 'turbtran'.
        !Hence, in this case, 'rcld' at level 'k_sf' is already present (probably in form of a tile-aggregation),
        ! and the 'rcld'-calculation should not be executed once again (based on somehow aggregated input variables).
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=i_st, i_en
           ! Std. deviation of local super-saturation (SDSS):
           rcld(i,k)=SQRT( tls(i,k)*lsh(i,k)*d_h )* &                   !related effective length scale "times"
!                     ABS( epr_qst_t*grd(i,k,tet_l) - grd(i,k,h2o_g) )  !related effective gradient of super-sat. water
                      ABS( exner(i,k)*qst_t(i,k)*grd(i,k,tet_l) - grd(i,k,h2o_g) )
                                                                        !related effective gradient of super-sat. water
        END DO
     END DO
     !$ACC END PARALLEL

     ! Calculating effective vertical gradients of the usual, not quasi-conservative, scalar variables:

     IF (lactcnv .AND. (tdc%lexpcor.OR.laddcnv)) THEN !flux conversion towards non-conserved variables required
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lacc)
        !$ACC LOOP SEQ
        DO k=k_st, k_sf !for all boundary levels including the extra lowermost boundary
!DIR$ IVDEP
          !$ACC LOOP GANG(STATIC: 1) VECTOR
           DO i=i_st, i_en
              ! Effective vertical gradient of water vapor content and pot. temper.:
              grd(i,k,vap) =  grd(i,k,h2o_g)-grd(i,k,liq)                     !eff_grad(vap)
              grd(i,k,tet) = (grd(i,k,tet_l)+grd(i,k,liq)*lhocp/exner(i,k))   !eff_grad(tet)

           END DO
        END DO
        !$ACC END PARALLEL

     END IF

     ! Correcting the effective vertical gradient under consideration of a fluctuating mixed-phase heat-capacity:

     IF (tdc%lcpfluc) THEN !the effect by fluctuating heat capacity has to be considered
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lacc)
        !$ACC LOOP SEQ
        DO k=k_st, k_sf !for all boundary levels including the extra lowermost boundary
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR
           DO i=i_st, i_en
              grd(i,k,tet)=grd(i,k,tet_l)*r_cpd(i,k) ! grad(tet)*(Cp/Cpd)
           END DO
        END DO
        !$ACC END PARALLEL
     END IF

    ! Calculating Eddy-Dissipation Rate (EDR):

     IF (lpres_edr .OR. tdc%ltmpcor) THEN
        ! Sichern der TKE-Dissipation "edr=q**3/l_dis" u.a. als thermische Quelle:
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lacc)
        !$ACC LOOP SEQ
        DO k=k_st, k_sf !for all boundary levels including the extra lowermost boundary
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
           DO i=i_st, i_en
              IF (lrogh_lay) THEN
                 wert=ediss(i,k)
              ELSE
                 wert=d_m
              END IF
!test<
!ediss(i,k)=tke(i,k,ntur)**3/(wert*tls(i,k))

              ediss(i,k)=MIN( tke(i,k,nvor), tke(i,k,ntur), tdc%vel_max )**3/(wert*tls(i,k))
!test>
           END DO
        END DO
        !$ACC END PARALLEL
     END IF

  END IF !only for the last iteration step

  !$ACC END DATA
   IF (.NOT. lcuda_graph_turb_tran) THEN
     !$ACC WAIT(acc_async_queue)
   END IF
  !$ACC END DATA

END SUBROUTINE solve_turb_budgets

!==============================================================================
!==============================================================================

!SCLM--------------------------------------------------------------------------
#ifdef SCLM
SUBROUTINE turb_stat (lsm, lsh, fm2, fh2, d_m, d_h, a_m, tls, tvs, tvt, grd, k)

INTiEGER, INTENT(IN) :: &
!
    k         !vertical level index

REAL (KIND=wp), INTENT(IN) :: &
!
    lsm,    & !turbulent length scale times value of stability function for momentum       [m]
    lsh,    & !turbulent length scale times value of stability function for scalars (heat) [m]
    fm2,    & !squared forcing frequency for momentum       [1/s2]
    fh2,    & !squared forcing frequency for scalars (heat) [1/s2]
    d_m,    & !dissipation parameter for momentum       [1]
    d_h,    & !dissipation parameter for scalars (heat) [1]
    a_m,    & !return-to-isotropy parameter for momentum [1]
    tls,    & !turbulent length scale   [m]
    tvs,    & !turbulent velocity scale [m/s]
    tvt,    & !turbulent transport of turbulent velocity scale [m/s2]
    grd(0:)    !vector of vertical gradients of diffused (conserved) variables [{unit}/m]
               !Note: Index "0" is reserved for the vertical gradient of CKE
               !       due to near-surface thermal inhomogeneity (being an related acceleration)

REAL (KIND=wp) ::  &
!
   tts,                  & !turbulent master time scale
   ts_d,                 & !dissipation time scale [s]
   ts_m,                 & !return-to-isotropy time scale for momentum [s]
   tvs2,                 & !tvs**2 [m2/s2]
   tkm, tkh,             & !turbulent diffusion coefficient for momentum and scalars (heat) [m2/s]
   x1, x2, x3,           & !auxilary TKE source terme values [m2/s3]
   cvar(nmvar,nmvar)       !covariance matrix  [{unit1}*{unit2}]
   !Notice that (according to module 'mo_turbdiff_config') 'u_m', 'v_m', 'tet_l', 'h2o_g' and 'w_m' are all
   ! equal to "1", "2", "3, "4" and "5=nmvar" respectively.

   cvar(tet_l,tet_l)=d_h*tls*lsh*grd(tet_l)**2
   cvar(h2o_g,h2o_g)=d_h*tls*lsh*grd(h2o_g)**2
   cvar(tet_l,h2o_g)=d_h*tls*lsh*grd(tet_l)*grd(h2o_g)

   tkm=lsm*tvs; tkh=lsh*tvs

   cvar(u_m  ,w_m)=-tkm*grd(u_m)
   cvar(v_m  ,w_m)=-tkm*grd(v_m)
   cvar(tet_l,w_m)=-tkh*grd(tet_l)
   cvar(h2o_g,w_m)=-tkh*grd(h2o_g)

   TKE_SCLM%mod(k)%val=tvs        ; TKE_SCLM%mod(k)%vst=i_cal

   !Achtung: TKE%mod(k)%val zeigt z.Z. noch auf die alte TKE-Zeitstufe.
   !Somit wird also der alte TKE-Wert mit dem neuen ueberschrieben,
   !was aber ohne Bedeutung ist, weil ab jetzt der alte tke-Wert nicht
   !mehr benoetigt wird. Beachte, dass die Modelleinheit hier [m/s] ist.

   tvs2=tvs**2
   tts=tls/tvs
   ts_d=d_m*tts
   ts_m=a_m*tts

   x1=cvar(u_m,w_m)*grd(u_m)
   x2=cvar(v_m,w_m)*grd(v_m)
   x3=-tkh*fh2

   BOYPR%mod(k)%val=x3             ; BOYPR%mod(k)%vst=i_cal
   SHRPR%mod(k)%val=-(x1+x2)       ; SHRPR%mod(k)%vst=i_cal
   DISSI%mod(k)%val=tvs2/ts_d      ; DISSI%mod(k)%vst=i_cal
   TRANP%mod(k)%val=tvs*tvt        ; TRANP%mod(k)%vst=i_cal

   cvar(u_m,u_m)=z1d3*tvs2+ts_m*(-4.0_wp*x1+2.0_wp*x2-2.0_wp*x3)
   cvar(v_m,v_m)=z1d3*tvs2+ts_m*(+2.0_wp*x1-4.0_wp*x2-2.0_wp*x3)
   cvar(w_m,w_m)=tvs2-cvar(u_m,u_m)-cvar(v_m,v_m)

   UWA%mod(k)%val=cvar(u_m  ,w_m)  ; UWA%mod(k)%vst=i_cal
   VWA%mod(k)%val=cvar(v_m  ,w_m)  ; VWA%mod(k)%vst=i_cal
   TWA%mod(k)%val=cvar(tet_l,w_m)  ; TWA%mod(k)%vst=i_cal
   UST%mod(k)%val=SQRT(tkm*SQRT(fm2))
                                     UST%mod(k)%vst=i_cal
 ! LMO%mod(0)%val=-UST%mod(k)%val**3/x3
 !                                   LMO%mod(0)%vst=i_cal

   IF (cvar(u_m,u_m).GE.0.0_wp .AND. cvar(v_m,v_m).GE.0.0_wp .AND. &
       cvar(u_m,u_m)+cvar(v_m,v_m).LE.tvs2) THEN

      UUA%mod(k)%val=cvar(u_m,u_m)  ; UUA%mod(k)%vst=i_cal
      VVA%mod(k)%val=cvar(v_m,v_m)  ; VVA%mod(k)%vst=i_cal
      WWA%mod(k)%val=tvs2-UUA%mod(k)%val-VVA%mod(k)%val
                                      WWA%mod(k)%vst=i_cal
   END IF

   TTA%mod(k)%val=cvar(tet_l,tet_l); TTA%mod(k)%vst=i_cal
   QQA%mod(k)%val=cvar(h2o_g,h2o_g); QQA%mod(k)%vst=i_cal
   TQA%mod(k)%val=cvar(tet_l,h2o_g); TQA%mod(k)%vst=i_cal

END SUBROUTINE turb_stat
#endif
!SCLM--------------------------------------------------------------------------

!==============================================================================
!==============================================================================

SUBROUTINE turb_cloud ( tdc, i1dim, ktp, &
!
   istart, iend, kstart, kend,           &
!
   icldtyp,                              &
!
   prs, t, qv, qc,                       &
   psf,                                  &
!
   rcld, & !inp: standard deviation of local super-saturation (SDSS)
           !out: saturation fraction (cloud-cover)
   clwc, & !out: liquid water content
   lacc, opt_acc_async_queue )

!------------------------------------------------------------------------------
!
! Description:
!
!     This routine calculates the area fraction of a grid box covered
!     by stratiform (non-convective) clouds.
!     If subgrid-scale condensation is required, an additional
!     saturation adjustment is done.
!
! Method:
!
!     "icldtyp = 0" : grid scale diagnosis of local super-saturation
!     Cloud water is estimated by grid scale saturation adjustment,
!     which equals the case "icldtyp = 2" in calse of no variance.
!     So cloud-cover is either "0" or "1".
!
!     "icldtyp = 1" : empirical diagnosis of local super-saturation
!     Fractional cloud-cover is determined empirically from
!     relative humidity. Also, an in-cloud water content of sugrid-scale
!     clouds is determined as a fraction of the saturation specific
!     humidity. Both quantities define the grid-volume mean cloud water
!     content.
!
!     "icldtyp = 2": statistical diagnosis of local super-saturation
!     A Gaussion distribution is assumed for the local super-saturation
!     'dq := qt - qs' where 'qt := qv + ql' is the total water content and
!     'qs' is the saturation specific humidity. Using the standard deviation
!     of this distribution SDSS (as input) and the quasi-conservative
!     quantities 'qt' and 'tl' (liquid water temperature), a corrected liquid
!     water content is determined, which contains also the contributions by
!     subgrid-scale cloud processes. A corresponding cloudiness is calculated as well.
!
!------------------------------------------------------------------------------

! Subroutine arguments
!----------------------

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

! Scalar arguments with intent(in):

INTEGER, INTENT(IN) :: &  ! indices used for allocation of arrays
  ktp,               & ! extra start index of vertical dimension for auxilary arrays (topping level of operations)
  istart, iend,      & ! zonal      start and end index
  kstart, kend,      & ! vertical   start and end index
  i1dim,             & ! length of blocks
!
  icldtyp              ! type of cloud diagnostics

! Array arguments with intent(in):

REAL (KIND=wp), DIMENSION(:,ktp:), TARGET, INTENT(IN) :: &
  prs,   & !atmospheric pressure [Pa]
  t, qv    !temperature [K] and water vapor content (spec. humidity) [1]

REAL (KIND=wp), DIMENSION(:,ktp:), TARGET, OPTIONAL, INTENT(IN) :: &
  qc       !cloud water content (if not present, 't' and 'qv' are conserved variables) [1]

REAL (KIND=wp), DIMENSION(:,:), INTENT(INOUT) :: &
  rcld     !inp: standard deviation of local super-saturation (SDSS) [1]
           !out: saturation fraction (cloud-cover) [1]

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
  psf      !surface pressure [Pa]

! Array arguments with intent(out):

REAL (KIND=wp), DIMENSION(:,ktp:), INTENT(INOUT) :: &
           !Attention: "INTENT(OUT)" might cause some not-intended default-setting for not-treated levels!
  clwc     ! liquid water content [1]

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

LOGICAL :: lzacc
INTEGER :: acc_async_queue

! Local variables and constants
! -----------------------------

INTEGER :: &
  i,k ! loop indices

REAL (KIND=wp), PARAMETER :: &
!Achtung: Erweiterung
  asig_max = 0.001_wp,   & ! max. absoute  standard deviation of local super-saturation
  rsig_max = 0.05_wp,    & ! max. relative standard deviation of local super-saturation
  zclwfak  = 0.005_wp,   & ! fraction of saturation specific humidity
  zuc      = 0.95_wp       ! constant for critical relative humidity

REAL (KIND=wp) :: &
  qs, dq, gam, sig

REAL (KIND=wp), DIMENSION(:,:), POINTER :: &
  qt, tl !total water and liquid water temperature

REAL (KIND=wp), DIMENSION(i1dim,kstart:kend), TARGET :: &
  qt_tar, tl_tar

REAL (KIND=wp) :: &
  pdry, q, & !
  zsigma, zclc0, zq_max, zq_inv !

!------------------------------------------------------------------
!Festlegung der Formelfunktionen fuer die Turbulenzparametrisierung:

!------------ End of header ----------------------------------------

!-------------------------------------------------------------------
! Note:
! If 'qc' is not present, 't' and 'qv' are assumed to be already liquid water temperature (tl)
!  and total water content (qt), which are conserved for moist conversions
! In case of "kstart=khi=kend", the vertical loop is only for one level.
! Since a surface variable for 'qc' is not used in general, it is not forseen here as well, assuming
!  't' and 'qv' to already contain the conserved variables.
! In this version, atmospheric ice is not included to the adjustment process.

! Begin Subroutine turb_cloud
! ---------------------------

  CALL set_acc_host_or_device(lzacc, lacc)

  IF(PRESENT(opt_acc_async_queue)) THEN
      acc_async_queue = opt_acc_async_queue
  ELSE
      acc_async_queue = 1
  ENDIF

!Achtung: Korrektur: zsig_max -> epsi
!Achtung: test
! zclc0=MIN( tdc%clc_diag, 1.0_wp-rsig_max )
  zclc0=MIN( tdc%clc_diag, 1.0_wp-tdc%epsi )

  zq_max = tdc%q_crit*(1.0_wp/zclc0 - 1.0_wp)
  zq_inv = 1.0_wp/(zq_max + tdc%q_crit)

  ! Local array
  !XL_ACCTMP : replace with allocatables wk array
  !$ACC DATA COPYIN(iend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC DATA CREATE(qt_tar, tl_tar) ASYNC(acc_async_queue) IF(lzacc .AND. PRESENT(qc))

  IF (PRESENT(qc)) THEN
     qt => qt_tar
     tl => tl_tar

     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
     !$ACC LOOP SEQ
     DO k = kstart, kend
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i = istart, iend
           qt(i,k) = qv(i,k) +       qc(i,k) ! total water content
           tl(i,k) =  t(i,k) - lhocp*qc(i,k) ! liquid water temperature
        END DO
     END DO
     !$ACC END PARALLEL

  ELSE !'qv' and 't' already contain conserved variablesi
     qt => qv
     tl => t
  END IF

  !Note: 'qt' and 'tl' are not being changed in the following!

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP SEQ
  DO k = kstart, kend
    !Calculation of saturation properties with respect to "t=tl" and "qv=qt":
!DIR$ IVDEP
    !$ACC LOOP GANG VECTOR PRIVATE(pdry, zsigma, q, qs, dq, gam, sig)
    DO i = istart, iend

!mod_2011/09/28: zpres=patm -> zpres=pdry {
      pdry=( 1.0_wp-qt(i,k) )/( 1.0_wp+rvd_m_o*qt(i,k) )*prs(i,k) ! part. pressure of dry air
      qs = zqvap( zpsat_w( tl(i,k) ), pdry )                      ! saturation mixing ratio (new version)
      gam = 1.0_wp/( 1.0_wp + lhocp*zdqsdt( tl(i,k), qs ) )       ! slope factor (from new vers. of d_qsat/d_T)
!mod_2011/09/28: zpres=patm -> zpres=pdry }

      dq = qt(i,k) - qs                                             ! local super-saturation

      IF (icldtyp.EQ.1) THEN
        ! Calculation of cloud-cover and cloud water content
        ! using an empirical relative humidity criterion

        zsigma = prs(i,k)/psf(i)

        ! critical relative humidity
        sig = zuc - uc1 * zsigma * ( 1.0_wp - zsigma )*( 1.0_wp + uc2*(zsigma-0.5_wp) )

        ! cloud-cover:
        rcld(i,k) = MAX( 0.0_wp, MIN( 1.0_wp, zclc0*( (ucl*qt(i,k)/qs-sig)/(ucl-sig) )**2 ) )

        ! grid-volume water content:
        clwc(i,k) = rcld(i,k)*qs*zclwfak

        IF (dq.GT.0.0_wp) THEN
          ! adapted grid-volume water content:
          clwc(i,k) = clwc(i,k) + (gam*dq-clwc(i,k))*(rcld(i,k)-zclc0)/(1.0_wp-zclc0)
        END IF

      ELSE !cloud water diagnosis based on saturation adjustment

        IF ( icldtyp.EQ.2 ) THEN
          ! Statistical calculation of cloud-cover and cloud water content
          ! using the standard deviation of the local super-saturation

!Achtung: Erweiterung, um COSMO-Version abzubilden
          IF (tdc%imode_stadlim.EQ.1) THEN !absolute upper SDSS-limit
            sig = MIN( asig_max, rcld(i,k) )
            qs = sig*zq_max
          ELSE !relative upper SDSS-limit
            sig = MIN( rsig_max*qs, rcld(i,k) )
            qs = MIN( qt(i,k), sig*zq_max )
          ENDIF
        ELSE !grid scale adjustment wihtout any variance
          sig=0.0_wp; qs=dq
        END IF

        q = MERGE( MERGE( -tdc%q_crit, zq_max, dq.LE.0.0_wp ), dq/sig, sig.LE.0.0_wp )
        ! In case of "sig=0", the method is similar to grid-scale saturation adjustment.
        ! Otherwise, a fractional cloud-cover (saturation fraction) is diagnosed.

        !cloud-water 'rcld' and liquid-water content 'clwc':
        sig = MIN( 1.0_wp, MAX ( 0.0_wp, (q + tdc%q_crit)*zq_inv ) )
        clwc(i,k) = gam*MERGE( dq, qs, q.GE.zq_max )*sig**2
        rcld(i,k) = sig

      END IF

    END DO
  END DO
  !$ACC END PARALLEL
  IF (.NOT. lcuda_graph_turb_tran) THEN
    !$ACC WAIT(acc_async_queue)
  END IF
  !$ACC END DATA
  !$ACC END DATA

END SUBROUTINE turb_cloud

!==============================================================================
!==============================================================================

SUBROUTINE vert_grad_diff ( tdc, kcm,                  &
!
          i_st, i_en, k_tp, k_sf,                      &
!
          dt_var, ivtype, igrdcon, itndcon,            &
          linisetup, lnewvtype,                        &
!++++
          lsflucond, lsfgrduse,                        &
          ldynimpwt, lprecondi, leff_flux,             &
!
          rho, rho_s, rho_n, hhl, r_air, tkv, tsv,     &
!
          disc_mom, expl_mom, impl_mom, invs_mom,      &
          diff_dep, diff_mom, invs_fac, scal_fac,      &
!
          dif_tend, cur_prof, eff_flux, lacc )
!++++

!------------------------------------------------------------------------------

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

INTEGER, INTENT(IN) :: &
!
! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------
!
    kcm,          &  ! level index of the upper roughness layer bound
!
! Start- and end-indices for the computations in the horizontal layers:
! -----------------------------------------------------------------------
    i_st, i_en, &    ! start and end index for meridional direction
    k_tp, k_sf, &    ! vertical level indices for top and surface
!
    ivtype,     &    ! index for variable type
!
    igrdcon,    &    ! mode index for effective gradient consideration
                     ! 0: no consider., 1: use only the related profile correction
                     !                  2: use complete effective profile
                     !                  3: !use additional effective gradients of an extra non-turbulent flux-contribution
    itndcon          ! mode index of current tendency  consideration
                     ! 0: no consider., 1: consider curr. tend. in implicit equation only
                     !                  2: add curr. tend. increment to current profile
                     !                  3: add related diffusion correction to current profile

REAL (KIND=wp), INTENT(IN) :: &
!
    dt_var           ! time increment for vertical diffusion         ( s )

LOGICAL, INTENT(IN) :: &
!
    linisetup,  &    ! calculate initial setup
    lnewvtype,  &    ! calculate setup for a new variable type
    lsflucond,  &    ! surface flux condition
!++++
    lsfgrduse,  &    ! use explicit surface gradient
    ldynimpwt,  &    ! dynamical calculatin of implicit weights
    lprecondi        ! preconditioning of tridiagonal matrix
!++++


LOGICAL, INTENT(INOUT) :: &
!
    leff_flux        ! calculation of effective flux density required

                  ! DIMENSION(ie,ke)
REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
    rho              ! air density at full levels                   (Kg/m3)

                  ! DIMENSION(ie,ke1)
REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: &
!
    hhl              ! half level height                              (m)

!Attention: Notice the start index of vertical dimension!
                  ! DIMENSION(ie,ke1)
REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: &
!
    tkv              ! turbulent diffusion coefficient               (m/s2)

                  ! DIMENSION(ie)
REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
    tsv              ! turbulent velocity at the surface              (m/s)

                  ! DIMENSION(ie,kcm-1:ke1)
REAL (KIND=wp), DIMENSION(:,kcm-1:), OPTIONAL, TARGET, INTENT(IN) :: &
!
    r_air            ! log of air containing volume fraction

                  ! DIMENSION(ie,ke1)
REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
!
    rho_n            ! air density at half levels                   (Kg/m3)

                  ! DIMENSION(ie)
REAL (KIND=wp), DIMENSION(:), OPTIONAL, INTENT(IN) :: &
!
    rho_s            ! air density at the surface                   (Kg/m3)

REAL (KIND=wp), TARGET, INTENT(INOUT), CONTIGUOUS :: &
!
!   Auxilary arrays:
!
!++++
    disc_mom(:,:), &      ! prep_inp calc_inp: discretis. momentum (rho*dz/dt) on var. levels
    expl_mom(:,:), &      ! prep_inp         : diffusion  momentum (rho*K/dz))
                     ! prep_out calc_inp: explicit part of diffusion momentum
    impl_mom(:,:), &      ! prep_out calc_inp: implicit  part of diffusion momentum
    invs_mom(:,:), &      ! prep_out calc_inp: inverted momentum
    invs_fac(:,:), &      ! prep_out calc_inp: inversion factor
    scal_fac(:,:), &      ! prep_out calc_inp: scaling factor due to preconditioning
    diff_dep(:,:)         ! diffusion depth

                  ! DIMENSION(ie,ke1)
REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: &
!
    diff_mom         ! aux: saved comlete diffusion momentum (only in case of "itndcon.EQ.3")
!++++

REAL (KIND=wp), INTENT(INOUT), CONTIGUOUS :: &
!
!   Inp-out-variable: DIMENSION(ie,ke1)
!
    cur_prof(:,:), &      ! inp     : current   variable profile (will be overwritten!)
                          ! calc_inp: corrected variable profile
                          ! calc_out: current   variable profile including tendency increment
    eff_flux(:,:), &      ! inp     : effective gradient
                          ! out     : effective flux density (if "leff_flux=T")
                          ! aux     : downward flux density and related vertical profile increment
    dif_tend(:,:)         ! inp     : current time tendency of variable
                          ! calc_inp: updated   variable profile including tendency increment
                          ! calc_out: updated   variable profile including increments by tendency and diffusion
                          ! out     : pure (vertically smoothed) diffusion tendency

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
LOGICAL :: lzacc

INTEGER :: &
!
    i,k, &
    k_lw, k_hi, &       ! vertical index of the lowest and highest level used for diffusion
          k_gc          ! vertical index of the uppermost level with gradient correction


REAL (KIND=wp) :: &
!
    fr_var, &           ! 1/dt_var
    tkmin               ! effective minimal diffusion coefficient

REAL (KIND=wp), DIMENSION(:,:), POINTER, CONTIGUOUS :: &
!
    rhon, rhoh

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------

  CALL set_acc_host_or_device(lzacc, lacc)

  k_lw=k_sf-1; k_hi=k_tp+1
               k_gc=k_hi+1

  fr_var=1.0_wp/dt_var

!++++
!  IF (ivtype.EQ.mom) THEN
!     tkmin=MAX( con_m, tkmmin )
!  ELSE
!     tkmin=MAX( con_h, tkhmin )
!  END IF
!++++

! Initial setup and adoptions for new variable type:

  IF (linisetup .OR. lnewvtype) THEN

     IF (linisetup .OR. .NOT.PRESENT(rho_n)) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=k_hi,k_lw
!DIR$ IVDEP
           DO i=i_st,i_en
              expl_mom(i,k)=hhl(i,k)-hhl(i,k+1)
           END DO
        END DO
        !$ACC END PARALLEL
     END IF

     IF (PRESENT(rho_n)) THEN
        rhon => rho_n
     ELSE
!++++
        rhon => invs_mom ; rhoh => rho
!++++

        CALL bound_level_interp( i_st,i_en, k_hi+1,k_lw, &
                                 nvars=1, pvar=(/varprf(rhon,rhoh)/), depth=expl_mom)

!DIR$ IVDEP
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           rhon(i,k_sf)=rho_s(i)
        END DO
        !$ACC END PARALLEL
     END IF

     IF (linisetup) THEN
!DIR$ IVDEP
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           disc_mom(i,k_hi)=rho(i,k_hi)*expl_mom(i,k_hi)*fr_var
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=k_hi+1,k_lw
!DIR$ IVDEP
           DO i=i_st,i_en
              disc_mom(i,k)=rho(i,k)*expl_mom(i,k)*fr_var
              diff_dep(i,k)=0.5_wp*(expl_mom(i,k-1)+expl_mom(i,k))
           END DO
        END DO
        !$ACC END PARALLEL
      END IF


     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
     DO k=k_hi+1,k_lw
!DIR$ IVDEP
        DO i=i_st,i_en
!Achtung: Eventuell tkmin-Beschraenkung nur bei VDiff
         ! expl_mom(i,k)=MAX( tkmin, tkv(i,k) ) &
           expl_mom(i,k)=tkv(i,k) &
                        *rhon(i,k)/diff_dep(i,k)
        END DO
     END DO
     !$ACC END PARALLEL
!DIR$ IVDEP
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st,i_en
!Achtung: Einfuehrung von 'tsv': macht Unterschiede
        expl_mom(i,k_sf)=rhon(i,k_sf)*tsv(i)
        diff_dep(i,k_sf)=tkv(i,k_sf)/tsv(i)
     !  expl_mom(i,k_sf)=rhon(i,k_sf)*tkv(i,k_sf)/diff_dep(i,k_sf)
     !  diff_dep(i,k_sf)=dzs(i)
     END DO
     !$ACC END PARALLEL
     !Attention: 'tkmin' should be excluded for the surface level 'k_sf'!

     IF (itndcon.EQ.3) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=k_hi+1,k_sf
!DIR$ IVDEP
           DO i=i_st,i_en
              diff_mom(i,k)=expl_mom(i,k)
           END DO
        END DO
        !$ACC END PARALLEL
     END IF

!    This manipulation enforces always a complete decoupling from the surface:
     IF (tdc%lfreeslip) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           expl_mom(i,k_sf)=0.0_wp
        END DO
        !$ACC END PARALLEL
     END IF

     CALL prep_impl_vert_diff( tdc, lsflucond, ldynimpwt, lprecondi, &
          i_st, i_en, k_tp=k_tp, k_sf=k_sf, &
          disc_mom=disc_mom, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
          invs_fac=invs_fac, scal_fac=scal_fac, &
          lacc=lzacc )

  END IF

! Optional correction of vertical profiles by vertical integration of correction gradients:

  IF (lsfgrduse) THEN !use effective surface value from effective surface gradient
!DIR$ IVDEP
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st,i_en
        cur_prof(i,k_sf)=cur_prof(i,k_sf-1)-diff_dep(i,k_sf)*eff_flux(i,k_sf)
     END DO
     !$ACC END PARALLEL
  END IF

  IF (igrdcon.EQ.3) THEN !use additional effective gradients of an extra non-turbulent flux-contribution
!DIR$ IVDEP
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st,i_en
        eff_flux(i,k_sf)=0._wp !non-turbulent circulation fluxes always vanish at the surface
     END DO
     !$ACC END PARALLEL
  END IF

  IF (igrdcon.EQ.1 .OR. igrdcon.EQ.3) THEN !applying correction profiles or adding profile corrections

     !Increments of profile values between adjacent levels:

     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO k=k_sf,k_gc,-1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           IF (igrdcon.EQ.1) THEN !only a diffusion correction at given and complete vertical gradients
              cur_prof(i,k)=cur_prof(i,k-1)-cur_prof(i,k) !positive increment between adjacent levels
           ELSE !full vertical diffusion using effective gradients of an extra non-local flux-contribution
              cur_prof(i,k)=cur_prof(i,k)-cur_prof(i,k-1) !negative increment between adjacent levels
           END IF
        END DO
     END DO
     !$ACC END PARALLEL

     !Top-down integration of effective (partial) gradients, starting with given values at level "k=k_gc-1":

     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO k=k_gc,k_sf
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
            cur_prof(i,k)=cur_prof(i,k-1)+(cur_prof(i,k)-diff_dep(i,k)*eff_flux(i,k))
        END DO
     END DO
     !$ACC END PARALLEL
  ELSEIF (igrdcon.EQ.2) THEN !effektive total profile
     !Bottom-up integration of full effective gradients starting with given values at level "k=k_lw"=k_sf-1=ke:

     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO k=k_lw,k_gc,-1 !bottom-up loop
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           cur_prof(i,k-1)=cur_prof(i,k)+diff_dep(i,k)*eff_flux(i,k) !bottom-up integration of full gradient
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

  IF (itndcon.EQ.3) THEN !add diffusion correction from tendency to current profile
     !Related downward flux densities:
!DIR$ IVDEP
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st,i_en
        eff_flux(i,2)=-dif_tend(i,1)*disc_mom(i,1)*dt_var
     END DO
     !$ACC END PARALLEL

     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO k=k_hi+1,k_lw
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           eff_flux(i,k+1)=eff_flux(i,k)-dif_tend(i,k)*disc_mom(i,k)*dt_var
        END DO
     END DO
     !$ACC END PARALLEL

     !Virtual total vertical increment:
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO k=k_hi+1,k_sf
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           eff_flux(i,k)=eff_flux(i,k)/diff_mom(i,k)+(cur_prof(i,k-1)-cur_prof(i,k))
        END DO
     END DO
     !$ACC END PARALLEL

     !Related corrected profile:
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO k=k_hi+1,k_sf
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
           cur_prof(i,k)=cur_prof(i,k-1)-eff_flux(i,k)
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

  IF (itndcon.GE.1) THEN !calculate updated profile by adding tendency increment to current profile
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO k=k_hi,k_lw
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=i_st,i_en
            dif_tend(i,k)=cur_prof(i,k)+dif_tend(i,k)*dt_var
        END DO
     END DO
     !$ACC END PARALLEL
!DIR$ IVDEP
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st,i_en
        dif_tend(i,k_sf)=cur_prof(i,k_sf)
     END DO
     !$ACC END PARALLEL
  END IF

  IF (kcm.LE.k_lw) THEN
     leff_flux = .TRUE.
  END IF

  !Final solution of the semi-implicit diffusion equation:

  !Note:
  !'cur_prof' is the current profile including the gradient correction (if "igrdcon>0")
  !           including the virtual gradient correction of an explicit time tendency (if "itndcon=3").
  !'dif_tend' is only used, if "itndcon>0" and contains 'cur_prof' updated by the curr. tend. incr..

  CALL calc_impl_vert_diff ( lsflucond, lprecondi, leff_flux, itndcon, &
       i_st, i_en ,k_tp=k_tp, k_sf=k_sf, &
       disc_mom=disc_mom, expl_mom=expl_mom, impl_mom=impl_mom, invs_mom=invs_mom, &
       invs_fac=invs_fac, scal_fac=scal_fac, cur_prof=cur_prof, upd_prof=dif_tend, eff_flux=eff_flux, &
       lacc=lzacc )

  !Note:
  !'cur_prof' now contains the current profile updated by the current tendency increment (if "itndcon>0").
  !'dif_tend' now contains the final updated profile including vertical diffusion.

  !Calculation of time tendencies for pure vertical diffusion:
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=k_hi,k_lw
!DIR$ IVDEP
     DO i=i_st,i_en
        dif_tend(i,k)=(dif_tend(i,k)-cur_prof(i,k))*fr_var
     END DO
  END DO
  !$ACC END PARALLEL

! Volume correction within the roughness layer:

  IF (PRESENT(r_air)) THEN
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
     DO k=kcm,k_lw  !r_air-gradient within the roughness layer
!DIR$ IVDEP
        DO i=i_st,i_en
           cur_prof(i,k)=(r_air(i,k)-r_air(i,k+1))/(dt_var*disc_mom(i,k))
           dif_tend(i,k)=dif_tend(i,k)                        &
                        +0.5_wp*(eff_flux(i,k+1)+eff_flux(i,k)) & !flux on full level
                             *cur_prof(i,k)                     !r_air-gradient
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

END SUBROUTINE vert_grad_diff

!==============================================================================
!==============================================================================

SUBROUTINE prep_impl_vert_diff ( tdc, lsflucond, ldynimpwt, lprecondi, &
!
   i_st,i_en, k_tp, k_sf, &
!
!  disc_mom, expl_mom, impl_mom, invs_mom, invs_fac, scal_fac, impl_weight, &
   disc_mom, expl_mom, impl_mom, invs_mom, invs_fac, scal_fac, &
   lacc )

!Achtung: l-Schleifen -> k-Schleifen
!Achtung: Vorzeichenwechsel fur impl. momentum ist uebersichtlicher

TYPE(t_turbdiff_config), POINTER, INTENT(IN) :: tdc ! 'turbdiff' configuration state for a single patch (domain)

INTEGER, INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for first  horizontal dimension
   k_tp,k_sf    !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1

LOGICAL, INTENT(IN) :: &
!
   lsflucond, & !flux condition at the surface
   ldynimpwt, & !dynamical calculatin of implicit weights
   lprecondi    !preconditioning of tridiagonal matrix

REAL (KIND=wp), INTENT(IN) :: &
!
   disc_mom(:,:)    !discretis momentum (rho*dz/dt) on variable levels

REAL (KIND=wp), INTENT(INOUT) :: &
!
   expl_mom(:,:)    !inp: diffusion momentum (rho*K/dz)
                    !out: explicit part of diffusion momentum

REAL (KIND=wp), INTENT(OUT) :: &
!
   impl_mom(:,:), & !scaled implicit part of diffusion momentum
   invs_mom(:,:), & !inverted momentum
   invs_fac(:,:), & !inversion vactor
   scal_fac(:,:)    !scaling factor due to preconditioning

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
LOGICAL :: lzacc

INTEGER :: &
!
   i,k, &           !horizontal and vertical coordiante indices
   m                !level increment dependent on type of surface condition

!-------------------------------------------------------------------------

   CALL set_acc_host_or_device(lzacc, lacc)

   IF (lsflucond) THEN !a surface-flux condition
      m=2
   ELSE
      m=1
   END IF

!  Implicit and explicit weights:

   IF (ldynimpwt) THEN !dynamical determination of implicit weights
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=k_tp+2, k_sf+1-m
!DIR$ IVDEP
         DO i=i_st, i_en
            impl_mom(i,k)=expl_mom(i,k) &
                            ! *MAX(MIN(expl_mom(i,k)/impl_mom(i,k), tdc%impl_s), tdc%impl_t)
                              *MAX(tdc%impl_s-0.5_wp*impl_mom(i,k)/expl_mom(i,k), tdc%impl_t)
         END DO
      END DO
      !$ACC END PARALLEL
   ELSE !use precalculated implicit weights
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=k_tp+2, k_sf+1-m
!DIR$ IVDEP
         DO i=i_st, i_en
!Achtung:
            impl_mom(i,k)=expl_mom(i,k)*tdc%impl_weight(k)
!impl_mom(i,k)=expl_mom(i,k)*1.00_wp
         END DO
      END DO
      !$ACC END PARALLEL
   END IF

   !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
   !$ACC LOOP GANG VECTOR COLLAPSE(2)
   DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
      DO i=i_st, i_en
         expl_mom(i,k)=expl_mom(i,k)-impl_mom(i,k)
      END DO
   END DO
   !$ACC END PARALLEL
   !Notice that 'expl_mom' still contains the whole diffusion momentum at level 'k_sf'!

!  Inverse momentum vector:

   IF (lprecondi) THEN !apply symmetric preconditioning of tridiagonal matrix
!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=i_st, i_en
         scal_fac(i,k_tp+1)=1.0_wp/SQRT(disc_mom(i,k_tp+1)+impl_mom(i,k_tp+2))
         invs_mom(i,k_tp+1)=1.0_wp
      END DO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=k_tp+2, k_sf-m
!DIR$ IVDEP
         DO i=i_st, i_en
            scal_fac(i,k)=1.0_wp/SQRT(disc_mom(i,k)+impl_mom(i,k)+impl_mom(i,k+1))
         END DO
      END DO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=k_sf+1-m, k_sf-1 !only for a surface-flux condition and at level k_sf-1
!DIR$ IVDEP
         DO i=i_st, i_en
            scal_fac(i,k)=1.0_wp/SQRT(disc_mom(i,k)+impl_mom(i,k))
         END DO
      END DO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR
         DO i=i_st, i_en
            impl_mom(i,k)=scal_fac(i,k-1)*scal_fac(i,k)*impl_mom(i,k)
            invs_fac(i,k)=invs_mom(i,k-1)*impl_mom(i,k)
            invs_mom(i,k)=1.0_wp/( 1.0_wp-invs_fac(i,k)*impl_mom(i,k) )
         END DO
      END DO
      !$ACC END PARALLEL
   ELSE !without preconditioning
!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=i_st, i_en
         invs_mom(i,k_tp+1)=1.0_wp/(disc_mom(i,k_tp+1)+impl_mom(i,k_tp+2))
      END DO
      !$ACC LOOP SEQ
      DO k=k_tp+2, k_sf-m
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=i_st, i_en
            invs_fac(i,k)=invs_mom(i,k-1)*impl_mom(i,k)
            invs_mom(i,k)=1.0_wp/( disc_mom(i,k)+impl_mom(i,k+1) &
                              +impl_mom(i,k)*(1.0_wp-invs_fac(i,k)) )
         END DO
      END DO
      !$ACC LOOP SEQ
      DO k=k_sf+1-m, k_sf-1 !only for a surface-flux condition and at level k_sf-1
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=i_st, i_en
            invs_fac(i,k)=invs_mom(i,k-1)*impl_mom(i,k)
            invs_mom(i,k)=1.0_wp/( disc_mom(i,k) &
                              +impl_mom(i,k)*(1.0_wp-invs_fac(i,k)) )
         END DO
      END DO
      !$ACC END PARALLEL
   END IF

END SUBROUTINE prep_impl_vert_diff

!==============================================================================
!==============================================================================

SUBROUTINE calc_impl_vert_diff ( lsflucond, lprecondi, leff_flux, itndcon, &
!
   i_st,i_en, k_tp, k_sf, &
!
   disc_mom, expl_mom, impl_mom, invs_mom, invs_fac, scal_fac, cur_prof, upd_prof, eff_flux, &
   lacc )

INTEGER, INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for first  horizontal dimension
   k_tp,k_sf, & !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1
!
   itndcon      ! mode index of current tendency  consideration
                ! 0: no consider., (+/-)1: consider curr. tend. in implicit equation only
                !                  (+/-)2: add curr. tend. increment to current profile
                !                  (+/-)3: add related diffusion correction to current profile
                ! for "itndcon>0", 'cur_prof' will be overwritten by input of 'upd_prof'

LOGICAL, INTENT(IN) :: &
!
   lsflucond, & !flux condition at the surface
   lprecondi, & !preconditioning of tridiagonal matrix
   leff_flux    !calculation of effective flux densities required

REAL (KIND=wp), INTENT(IN) :: &
!
   expl_mom(:,:), & !(explicit part of) diffusion momentum (rho*K/dz)
   impl_mom(:,:), & !(scaled) implicit part of diffusion momentum
   disc_mom(:,:), & !discretis momentum (rho*dz/dt)
!
   invs_mom(:,:), & !inverted momentum
   invs_fac(:,:), & !inversion factor
   scal_fac(:,:)    !scaling factor due to preconditioning

REAL (KIND=wp), TARGET, CONTIGUOUS, INTENT(INOUT) :: &
!
   cur_prof(:,:), & !inp: current vertical variable profile (including gradient corrections)
                    !out: current vertical variable profile (including tendency increment, if "itndcon>=1")
   upd_prof(:,:), & !inp: if "|itndcon|>=1": updated vertical variable profile (includ. tendency increment)
                    !aux: interim solution of diffusion equation
                    !out: updated vertical variable profile by vertical diffusion
   !"itndcon>= 1": 'cur_prof' as output equals 'upd_prof' as input
   !"itndcon<=-1": 'cur_prof' keeps input profile
   !"itndcon = 0": 'upd_prof' is not used as input
!
   eff_flux(:,:)    !out: effective flux density belonging to (semi-)implicit diffusion
                    !     (positiv downward and only, if "leff_flux=.TRUE.")
                    !aux: explicit flux density (positiv upward) and
                    !     full right-hand side of the semi-implicit diffusion equation

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
LOGICAL :: lzacc

INTEGER :: &
!
   i,k              !horizontal and vertical coordiante indices

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
!
   old_prof(:,:), & !old variable profile value>
   rhs_prof(:,:)    !profile used at RHS of tridiag. system

!-----------------------------------------------------------------------

!  Preparation:

   CALL set_acc_host_or_device(lzacc, lacc)

  SELECT CASE ( ABS(itndcon) ) !discriminate modes of current tendency consideration
  CASE (0) ! no consideration of current tendency
    old_prof => cur_prof
    rhs_prof => cur_prof
  CASE (1) ! consideration of current tendency in implicit part only
    old_prof => upd_prof
    rhs_prof => cur_prof
  CASE (2) ! consideration of current tendency in all parts
    old_prof => upd_prof
    rhs_prof => upd_prof
  CASE (3) ! current profile already contains related gradient correction
    old_prof => cur_prof
    rhs_prof => cur_prof
  END SELECT


  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=k_tp+2, k_sf
!DIR$ IVDEP
     DO i=i_st, i_en
        eff_flux(i,k) = expl_mom(i,k) * ( rhs_prof(i,k  ) - rhs_prof(i,k-1) )
     END DO
  END DO
  !$ACC END PARALLEL
  !Notice that 'expl_mom(i,k_sf)' still contains the whole diffusion momentum!

!Achtung: Korrektur: Richtige Behandlung der unteren Konzentrations-Randbedingung
  IF (.NOT.lsflucond) THEN !only for a surface-concentration condition
!DIR$ IVDEP
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st, i_en
        eff_flux(i,k_sf) = eff_flux(i,k_sf) + impl_mom(i,k_sf) * rhs_prof(i,k_sf-1)
     END DO
     !$ACC END PARALLEL
     !Note:
     !At level 'k_sf' 'impl_mom' still contains the implicit part without scaling,
     ! and it vanishes at all in case of "lsflucond=T" (surface-flux condition)!
     !In all, level "k_sf-1" is treated (semi-)implicitly by this.
  END IF

! Resultant right-hand side flux:

  k=k_tp+1
!DIR$ IVDEP
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i=i_st, i_en
     eff_flux(i,k  ) = disc_mom(i,k  ) * old_prof(i,k ) + eff_flux(i,k+1)
  END DO
  !Note: Zero flux condition just below top level.

  !$ACC LOOP SEQ
  DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
        eff_flux(i,k  ) = disc_mom(i,k  ) * old_prof(i,k  ) + eff_flux(i,k+1) - eff_flux(i,k  )
     END DO
  END DO
  !$ACC END PARALLEL

  IF (lprecondi) THEN !preconditioning is active
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
     DO k=k_tp+1, k_sf-1
!DIR$ IVDEP
        DO i=i_st, i_en
           eff_flux(i,k  ) = scal_fac(i,k  ) * eff_flux(i,k  )
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

!  Save updated profiles (including explicit increments of current tendencies):

  IF (itndcon.GT.0) THEN !consideration of explicit tendencies
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
     DO k=k_tp+1, k_sf-1
!DIR$ IVDEP
        DO i=i_st, i_en
           cur_prof(i,k) = upd_prof(i,k)
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

!  Forward substitution:

   k=k_tp+1
!DIR$ IVDEP
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i=i_st, i_en
     upd_prof(i,k  ) = eff_flux(i,k  ) * invs_mom(i,k  )
  END DO

  !$ACC LOOP SEQ
  DO k=k_tp+2, k_sf-1
!DIR$ IVDEP
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
        upd_prof(i,k  ) = ( eff_flux(i,k  ) + impl_mom(i,k  ) * upd_prof(i,k-1) ) * invs_mom(i,k  )
     END DO
  END DO

!  Backward substitution:
  !$ACC LOOP SEQ
  DO k=k_sf-2, k_tp+1, -1
!DIR$ IVDEP
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
        upd_prof(i,k) = upd_prof(i,k  ) + invs_fac(i,k+1) * upd_prof(i,k+1)
     END DO
  END DO
  !$ACC END PARALLEL

  IF (lprecondi) THEN !preconditioning is active
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
     DO k=k_tp+1, k_sf-1
!DIR$ IVDEP
        DO i=i_st, i_en
           upd_prof(i,k  ) = scal_fac(i,k  ) * upd_prof(i,k  )
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

   !Note:
   !'cur_prof(i,k  )' is now the current profile including optional explicit tendencies
   !'upd_prof(i,k  )' contains the new profile value after diffusion

!  Effective flux density by vertical integration of diffusion tendencies:

  IF (leff_flux) THEN
     k=k_tp+1
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     !$ACC LOOP GANG(STATIC: 1) VECTOR
     DO i=i_st, i_en
       eff_flux(i,k) = 0.0_wp !upper zero-flux condition
     END DO

     !$ACC LOOP SEQ
     DO k=k_tp+2, k_sf
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO i=i_st, i_en
           eff_flux(i,k  ) = eff_flux(i,k-1) + (cur_prof(i,k-1) - upd_prof(i,k-1)) * disc_mom(i,k-1)
        END DO
     END DO
     !$ACC END PARALLEL
     !Note:
     !'eff_flux' is the vertical flux density of pure diffusion now (positive downward)
  END IF

END SUBROUTINE calc_impl_vert_diff

!==============================================================================
!==============================================================================

SUBROUTINE vert_smooth( i_st, i_en, k_tp, k_sf, &
                        cur_tend, disc_mom, vertsmot, smotfac, &
                        luse_mask, lacc )

INTEGER, INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for horizontal dimension
   k_tp,k_sf    !vertical level indices for top and surface
                !full level vars: k_tp=0; k_sf=ke+1
                !half level vars: k_tp=1; k_sf=ke+1

REAL (KIND=wp), INTENT(IN) :: &
!
   disc_mom(:,:), & !discretised momentum (rho*dz/dt) on variable levels
!
   vertsmot         !vertical smoothing factor

REAL (KIND=wp), INTENT(IN), OPTIONAL :: smotfac(:)

REAL (KIND=wp), INTENT(INOUT) :: &
!
   cur_tend(:,:)    !inp: current  vertical tendency profile
                    !out: smoothed vertical tendency profile
   !modifications takes place at levels k_tp+1 until k_sf-1.

LOGICAL, OPTIONAL, INTENT(IN) :: luse_mask

LOGICAL, INTENT(IN) :: lacc

LOGICAL :: lcond

LOGICAL :: lzacc

INTEGER :: &
!
   i,k, &           !horizontal and vertical coordiante indices
   j0,j1,j2         !altering indices "1" and "2"

REAL (KIND=wp) :: &
!
   sav_tend(SIZE(cur_tend,1),2), & !saved tendencies of current level
!
   versmot(SIZE(cur_tend,1)), & !smoothing factor of current level
   remfact(SIZE(cur_tend,1))    !remaining factor of current level

!----------------------------------------------------------------------------

!locals XL_GPUOPT: make this variables allocatables
   !$ACC DATA CREATE(sav_tend, versmot, remfact) PRESENT(cur_tend, disc_mom) IF(lacc)
   ! This conditional-create + no_create combo is to omit the ACC present checks.
   !$ACC DATA NO_CREATE(sav_tend, versmot, remfact, cur_tend, disc_mom, smotfac)

   IF (PRESENT(luse_mask)) THEN
     lcond=(PRESENT(smotfac) .AND. luse_mask)
   ELSE
     lcond=PRESENT(smotfac)
   END IF
   IF (lcond) THEN
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st,i_en
        versmot(i) = vertsmot*smotfac(i)
     END DO
     !$ACC END PARALLEL
   ELSE
     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lacc)
     !$ACC LOOP GANG VECTOR
     DO i=i_st,i_en
        versmot(i) = vertsmot
     END DO
     !$ACC END PARALLEL
   ENDIF

   !$ACC PARALLEL DEFAULT(PRESENT) PRIVATE(j0, j1, j2, k) ASYNC(1) IF(lacc)
   !$ACC LOOP GANG(STATIC: 1) VECTOR
   DO i=i_st,i_en
      remfact(i)=1.0_wp-versmot(i)
   END DO

   k=k_tp+1
   j1=1; j2=2
!DIR$ IVDEP
   !$ACC LOOP GANG(STATIC: 1) VECTOR
   DO i=i_st,i_en
      sav_tend(i,j1)=cur_tend(i,k)
      cur_tend(i,k) =remfact(i)* cur_tend(i,k)                   &
                    +versmot(i)* cur_tend(i,k+1)*disc_mom(i,k+1) &
                                                /disc_mom(i,k)
   END DO

   !$ACC LOOP GANG(STATIC: 1) VECTOR
   DO i=i_st,i_en
      remfact(i)=1.0_wp-2.0_wp*versmot(i)
   END DO

   !$ACC LOOP SEQ
   DO k=k_tp+2, k_sf-2
      j0=j1; j1=j2; j2=j0
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=i_st,i_en
         sav_tend(i,j1)=cur_tend(i,k)
         cur_tend(i,k) =remfact(i)* cur_tend(i,k)                    &
                       +versmot(i)*(sav_tend(i,j2) *disc_mom(i,k-1)  &
                                   +cur_tend(i,k+1)*disc_mom(i,k+1)) &
                                                   /disc_mom(i,k)
      END DO
   END DO

   !$ACC LOOP GANG(STATIC: 1) VECTOR
   DO i=i_st,i_en
      remfact(i)=1.0_wp-versmot(i)
   END DO

   k=k_sf-1
   j2=j1
!DIR$ IVDEP
   !$ACC LOOP GANG(STATIC: 1) VECTOR
   DO i=i_st,i_en
      cur_tend(i,k) =remfact(i)* cur_tend(i,k)                   &
                    +versmot(i)* sav_tend(i,j2) *disc_mom(i,k-1) &
                                                /disc_mom(i,k)
   END DO
   !$ACC END PARALLEL

   !$ACC WAIT IF(lacc)
   !$ACC END DATA
   !$ACC END DATA

END SUBROUTINE vert_smooth

!==============================================================================
!==============================================================================

SUBROUTINE bound_level_interp (i_st,i_en, k_st, k_en, &
                               nvars, pvar, depth, auxil, lacc)

INTEGER, INTENT(IN) :: &
!
   i_st,i_en, & !start and end index for horizontal dimension
   k_st,k_en, & !start and end index for vertical dimension
!
   nvars        !number of variables to be interpolated onto bound. levels

TYPE (varprf) :: pvar(1:) !'%ml': present variable profiles on ml as input
                          !'%bl': interpol. varib. profiles on bl as output

REAL (KIND=wp), TARGET, OPTIONAL, INTENT(IN) :: &
!
   depth(:,:)    !layer depth used for interpolation

REAL (KIND=wp), TARGET, OPTIONAL, INTENT(INOUT) :: &
!
   auxil(:,:)    !target for layer depth, when 'depth' contains boundary level height

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
LOGICAL :: lzacc

INTEGER :: i,k,n

LOGICAL :: ldepth, lauxil

REAL (KIND=wp), POINTER, CONTIGUOUS :: blvar(:,:), mlvar(:,:) !facilitates loop unrolling

   ldepth=PRESENT(depth) !'depth' has to be used
   lauxil=PRESENT(auxil) !'depth' contains boundary level height

!-------------------------------------------------------------------------

   ! WARNING: pvar%bl and pvar%ml MAY ALIAS!!
   ! OpenACC loops can't be collapsed!

   CALL set_acc_host_or_device(lzacc, lacc)

   ! OpenACC attachment
   !$ACC DATA CREATE(pvar) ASYNC(1) IF(lzacc)
   !$ACC DATA PRESENT(blvar, mlvar) IF(lzacc)

   DO n=1,nvars
#ifdef _PGI_LEGACY_WAR
      IF(lzacc) THEN
         !$ACC WAIT(1)
         CALL acc_attach(pvar(n)%bl)
         CALL acc_attach(pvar(n)%ml)
      END IF
#else
      !$ACC ENTER DATA ATTACH(pvar(n)%bl, pvar(n)%ml) ASYNC(1) IF(lzacc)
#endif
   END DO

   IF (ldepth) THEN !depth weighted interpolation

      IF (lauxil) THEN !precalculation of the reciprocal layer depth
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
         !$ACC LOOP SEQ
         DO k=k_en, k_st, -1
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=i_st, i_en
               auxil(i,k)=depth(i,k-1)/(depth(i,k-1)+depth(i,k))
            END DO
         END DO
         !$ACC END PARALLEL

         DO n=1, nvars
            blvar => pvar(n)%bl; mlvar => pvar(n)%ml !facilitates loop-unrolling
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP SEQ
            DO k=k_en, k_st, -1
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=i_st, i_en
                  blvar(i,k)=mlvar(i,k)*auxil(i,k)+mlvar(i,k-1)*(1._wp-auxil(i,k))
               END DO
            END DO
            !$ACC END PARALLEL
         END DO

      ELSE !no precalculation

         DO n=1, nvars
            blvar => pvar(n)%bl; mlvar => pvar(n)%ml !facilitates loop-unrolling
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP SEQ
            DO k=k_en, k_st, -1
!DIR$ IVDEP
               !$ACC LOOP GANG VECTOR
               DO i=i_st, i_en
                  blvar(i,k)=zbnd_val(mlvar(i,k), mlvar(i,k-1), depth(i,k), depth(i,k-1))
               END DO
            END DO
            !$ACC END PARALLEL
         END DO
      END IF

   ELSE !inverse of main level interpolation

      DO n=1, nvars
         blvar => pvar(n)%bl; mlvar => pvar(n)%ml !facilitates loop-unrolling
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
         !$ACC LOOP SEQ
         DO k=k_st, k_en
!DIR$ IVDEP
            !$ACC LOOP GANG VECTOR
            DO i=i_st, i_en
               blvar(i,k)=2.0_wp*mlvar(i,k)-mlvar(i,k+1)
            END DO
         END DO
         !$ACC END PARALLEL
      END DO

   END IF

   ! See comment above
   DO n=1,nvars
#ifdef _PGI_LEGACY_WAR
      IF(lzacc) THEN
         !$ACC WAIT(1)
         CALL acc_detach(pvar(n)%bl)
         CALL acc_detach(pvar(n)%ml)
      END IF
#else
      !$ACC EXIT DATA DETACH(pvar(n)%bl, pvar(n)%ml) ASYNC(1) IF(lzacc)
#endif
   END DO
   !$ACC WAIT(1)
   !$ACC END DATA
   !$ACC END DATA

END SUBROUTINE bound_level_interp

!==============================================================================
!==============================================================================

PURE FUNCTION  zbnd_val (val1, val2, dep1, dep2)
  !$ACC ROUTINE SEQ

  REAL (KIND=wp), INTENT(IN) :: val1, val2, dep1, dep2
  REAL (KIND=wp) :: zbnd_val

  zbnd_val=( val1*dep2+val2*dep1 )/( dep1+dep2 )

END FUNCTION zbnd_val

PURE FUNCTION zexner (zpres) !Exner-factor
  !$ACC ROUTINE SEQ

  REAL (KIND=wp), INTENT(IN) :: zpres
  REAL (KIND=wp) :: zexner

  zexner = EXP(rdocp*LOG(zpres/p0ref))

END FUNCTION zexner

PURE FUNCTION zpsat_w (ztemp) !satur. vapor pressure over water
  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN) :: ztemp
  REAL (KIND=wp) :: zpsat_w

  zpsat_w = b1*EXP(b2w*(ztemp-b3)/(ztemp-b4w))

END FUNCTION zpsat_w

PURE FUNCTION zqvap_old (zpvap, zpres) !satur. specif. humid. (old version)
  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN) :: zpvap, &
                                    zpres !pressure
  REAL (KIND=wp) :: zqvap_old

  zqvap_old = rdv*zpvap/(zpres-o_m_rdv*zpvap) !old form

END FUNCTION zqvap_old

PURE FUNCTION zqvap (zpvap, zpdry) !satur. specif. humid. (new version)
  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN) :: zpvap, &
                                zpdry !part. pressure of dry air
  REAL (KIND=wp) :: zqvap

!mod_2011/09/28: zpres=patm -> zpres=pdry {
  zqvap = rdv*zpvap
  zqvap = zqvap/(zpdry+zqvap)
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zqvap

PURE FUNCTION zdqsdt_old (ztemp, zqsat) !d_qsat/d_temp (old version)
  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN) :: ztemp, zqsat
  REAL (KIND=wp) :: zdqsdt_old

  zdqsdt_old=b234w*(1.0_wp+rvd_m_o*zqsat) & !old form
                  *zqsat/(ztemp-b4w)**2         !old form

END FUNCTION zdqsdt_old

PURE FUNCTION zdqsdt (ztemp, zqsat) !d_qsat/d_tem (new version)
  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN) :: ztemp, zqsat
  REAL (KIND=wp) :: zdqsdt

!mod_2011/09/28: zpres=patm -> zpres=pdry {
  zdqsdt=b234w*(1.0_wp-zqsat)*zqsat/(ztemp-b4w)**2
!mod_2011/09/28: zpres=patm -> zpres=pdry }

END FUNCTION zdqsdt

!==============================================================================

END MODULE turb_utilities
