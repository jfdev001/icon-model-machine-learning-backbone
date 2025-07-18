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

!NEC$ options "-finline-max-depth=3 -finline-max-function-size=2000"

!
! Two-moment bulk microphysics after Seifert, Beheng and Blahak
!
! Description:
! Provides various modules and subroutines for two-moment bulk microphysics
!

MODULE mo_2mom_mcrph_setup

  USE mo_kind,               ONLY: sp, wp
  USE mo_exception,          ONLY: finish, message, txt => message_text
  USE mo_math_constants,     ONLY: pi, pi4 => pi_4
  USE mo_physical_constants, ONLY: &
       & R_l   => rd,     & ! gas constant of dry air (luft)
       & R_d   => rv,     & ! gas constant of water vapor (dampf)
       & cp    => cpd,    & ! specific heat capacity of air at constant pressure
       & c_w   => clw,    & ! specific heat capacity of water
       & L_wd  => alv,    & ! specific heat of vaporization (wd: wasser->dampf)
       & L_ed  => als,    & ! specific heat of sublimation (ed: eis->dampf)
       & L_ew  => alf,    & ! specific heat of fusion (ew: eis->wasser)
       & T_3   => tmelt,  & ! melting temperature of ice
       & rho_w => rhoh2o, & ! density of liquid water
       & rho_ice => rhoice,&! density of pure ice
       & nu_l  => con_m,  & ! kinematic viscosity of air
       & D_v   => dv0,    & ! diffusivity of water vapor in air at 0 C
       & K_t   => con0_h, & ! heat conductivity of air
       & N_avo => avo,    & ! Avogadro number [1/mol]
       & k_b   => ak,     & ! Boltzmann constant [J/K]
       & grav               ! acceleration due to Earth's gravity
  USE mo_thdyn_functions, ONLY:     &
       & e_ws  => sat_pres_water,  & ! saturation pressure over liquid water
       & sat_pres_ice                ! saturation pressure over ice
  USE mo_2mom_mcrph_types, ONLY: &
       & particle, particle_frozen, particle_lwf, atmosphere, &
       & particle_sphere, particle_rain_coeffs, particle_cloud_coeffs, aerosol_ccn, &
       & particle_ice_coeffs, particle_snow_coeffs, particle_graupel_coeffs, &
       & particle_coeffs, collection_coeffs, rain_riming_coeffs, dep_imm_coeffs, &
       & coll_coeffs_ir_pm ! , lookupt_1D, lookupt_4D
  USE mo_fortran_tools,      ONLY: init
  USE mo_2mom_mcrph_config,  ONLY: t_cfg_2mom


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_2mom_mcrph_setup'

  TYPE(t_cfg_2mom) :: cfg_params !.. Container to hold some config params for the actual 2-mom call

  ! .. some cloud physics parameters
  REAL(wp), PARAMETER :: N_sc = 0.710_wp        !..Schmidt-Zahl (PK, S.541)
  REAL(wp), PARAMETER :: n_f  = 0.333_wp        !..Exponent von N_sc im Vent-koeff. (PK, S.541)

  INTEGER, PARAMETER :: cloud_type_default_gscp4 = 2603, ccn_type_gscp4 = 7
  INTEGER, PARAMETER :: cloud_type_default_gscp5 = 2603, ccn_type_gscp5 = 8

  ! AS: For gscp=4 use 2103 with ccn_type = 1 (HDCP2 IN and CCN schemes)
  !     For gscp=5 use 2603 with ccn_type = 8 (PDA ice nucleation and Segal&Khain CCN activation)

  ! AS: Runs without hail, e.g, 1503 are buggy and give a segmentation fault.
  !     So far I was not able to identify the problem, needs more detailed debugging.

  ! debug switches
  LOGICAL, PARAMETER     :: isdebug = .false.   ! use only when really desperate
  LOGICAL, PARAMETER     :: isprint = .true.    ! print-out initialization values

  REAL(wp), PARAMETER    :: pi6 = pi/6.0_wp, pi8 = pi/8.0_wp ! more pieces of pi

  ! Functions
  PUBLIC :: particle_mass, particle_meanmass, particle_diameter, particle_normdiameter
  PUBLIC :: particle_assign, particle_frozen_assign, particle_lwf_assign
  PUBLIC :: particle_velocity, particle_lwf_idx
  PUBLIC :: rain_mue_dm_relation
  PUBLIC :: moment_gamma
  PUBLIC :: coll_delta_11, coll_delta_12
  PUBLIC :: coll_theta_11, coll_theta_12
  ! Setup of processes
  PUBLIC :: setup_particle_coeffs, setup_cloud_autoconversion_sb
  PUBLIC :: setup_ice_selfcollection, setup_snow_selfcollection, setup_graupel_selfcollection
  PUBLIC :: setup_particle_collection_type1, setup_particle_collection_type2
  PUBLIC :: setup_particle_coll_pm_type1, setup_particle_coll_pm_type1_bfull
  PUBLIC :: set_ccn_cloud_type
  ! Constants
  PUBLIC :: n_f, N_sc
  PUBLIC :: cloud_type_default_gscp4, ccn_type_gscp4, cloud_type_default_gscp5, ccn_type_gscp5
  ! Types
  PUBLIC :: cfg_params

CONTAINS

  !*******************************************************************************
  ! Functions and subroutines working on particle class
  !*******************************************************************************
  ! (1) CLASS procedures for particle class
  !*******************************************************************************

  subroutine particle_assign(that,this)
    CLASS(particle), INTENT(in)   :: this
    TYPE(particle), INTENT(inout) :: that

    that%name = this%name
    that%nu = this%nu
    that%mu = this%mu
    that%x_max = this%x_max
    that%x_min = this%x_min
    that%a_geo = this%a_geo
    that%b_geo = this%b_geo
    that%a_vel = this%a_vel
    that%b_vel = this%b_vel
    that%a_ven = this%a_ven
    that%b_ven = this%b_ven
    that%cap   = this%cap
    that%vsedi_max = this%vsedi_max
    that%vsedi_min = this%vsedi_min
  END subroutine particle_assign

  subroutine particle_frozen_assign(that,this)
    TYPE(particle_frozen), INTENT(in)    :: this
    TYPE(particle_frozen), INTENT(inout) :: that

    that%name = this%name
    that%nu = this%nu
    that%mu = this%mu
    that%x_max = this%x_max
    that%x_min = this%x_min
    that%a_geo = this%a_geo
    that%b_geo = this%b_geo
    that%a_vel = this%a_vel
    that%b_vel = this%b_vel
    that%a_ven = this%a_ven
    that%b_ven = this%b_ven
    that%cap   = this%cap
    that%vsedi_max = this%vsedi_max
    that%vsedi_min = this%vsedi_min
    that%ecoll_c   = this%ecoll_c
    that%D_crit_c  = this%D_crit_c
    that%q_crit_c  = this%q_crit_c
    that%s_vel     = this%s_vel
  END subroutine particle_frozen_assign

  subroutine particle_lwf_assign(that,this)
    TYPE(particle_lwf), INTENT(in)    :: this
    TYPE(particle_lwf), INTENT(inout) :: that

    that%name = this%name
    that%nu = this%nu
    that%mu = this%mu
    that%x_max = this%x_max
    that%x_min = this%x_min
    that%a_geo = this%a_geo
    that%b_geo = this%b_geo
    that%a_vel = this%a_vel
    that%b_vel = this%b_vel
    that%a_ven = this%a_ven
    that%b_ven = this%b_ven
    that%cap   = this%cap
    that%vsedi_max = this%vsedi_max
    that%vsedi_min = this%vsedi_min

    that%ecoll_c   = this%ecoll_c
    that%D_crit_c  = this%D_crit_c
    that%q_crit_c  = this%q_crit_c
    that%s_vel     = this%s_vel

    that%lwf_cnorm1 = this%lwf_cnorm1
    that%lwf_cnorm2 = this%lwf_cnorm2
    that%lwf_cnorm3 = this%lwf_cnorm3
    that%lwf_cmelt1 = this%lwf_cmelt1
    that%lwf_cmelt2 = this%lwf_cmelt2
  END subroutine particle_lwf_assign

  ! mean mass with limiters, Eq. (94) of SB2006
  ELEMENTAL FUNCTION particle_meanmass(this,q,n) RESULT(xmean)

    !$ACC ROUTINE SEQ

    CLASS(particle), INTENT(in) :: this
    REAL(wp),        INTENT(in) :: q, n
    REAL(wp)                    :: xmean
    REAL(wp), PARAMETER         :: eps = 1e-20_wp

    xmean = MIN(MAX(q/(n+eps),this%x_min),this%x_max)
  END FUNCTION particle_meanmass

  ! mass-diameter relation, power law, Eq. (32) of SB2006
  ELEMENTAL FUNCTION particle_diameter(this,x) RESULT(D)

    !$ACC ROUTINE SEQ

    CLASS(particle), INTENT(in) :: this
    REAL(wp),        INTENT(in) :: x
    REAL(wp)                    :: D

    D = this%a_geo * EXP(this%b_geo*LOG(x))    ! D = a_geo * x**b_geo
  END FUNCTION particle_diameter

  ! inverse mass-diameter relation, power law, Eq. (32) of SB2006
  ELEMENTAL FUNCTION particle_mass(this,D) RESULT(x)

    !$ACC ROUTINE SEQ

    CLASS(particle), INTENT(in) :: this
    REAL(wp),        INTENT(in) :: D
    REAL(wp)                    :: x

    x = EXP((1.0_wp/this%b_geo)*LOG(D/this%a_geo))    ! x = (D/a_geo)**(1/b_geo)
  END FUNCTION particle_mass

  ! normalized diameter for rational function approx.
  ! in lwf-melting scheme
  PURE FUNCTION particle_normdiameter(this,D_m) RESULT(dnorm)
    CLASS(particle_lwf), intent(in) :: this
    REAL(wp), INTENT(in) :: D_m
    REAL(wp)             :: dnorm

    dnorm = MIN(MAX((LOG10(D_m*this%lwf_cnorm1)+this%lwf_cnorm2)*this%lwf_cnorm3,0.0_wp),1.0_wp)
    RETURN
  END FUNCTION particle_normdiameter

  ! lwf of mixed particle
  PURE FUNCTION particle_lwf_idx(this,i,j) RESULT(lwf)
    CLASS(particle_lwf), INTENT(in) :: this
    INTEGER,         INTENT(in) :: i,j
    REAL(wp)                    :: lwf
    REAL(wp), PARAMETER         :: eps = 1e-20_wp

    lwf = MAX(MIN(this%l(i,j)/(this%q(i,j)+eps),1.0_wp),0.0_wp)
    RETURN
  END FUNCTION particle_lwf_idx

  ! terminal fall velocity of particles, cf. Eq. (33) of SB2006
! Cray compiler does not support OpenACC in elemental or pure
#ifdef _CRAYFTN
  FUNCTION particle_velocity(this,x) RESULT(v)
#else
  ELEMENTAL FUNCTION particle_velocity(this,x) RESULT(v)
#endif
    !$ACC ROUTINE SEQ

    CLASS(particle), INTENT(in) :: this
    REAL(wp),        INTENT(in) :: x
    REAL(wp)                    :: v

    v = this%a_vel * EXP(this%b_vel * LOG(x))  ! v = a_vel * x**b_vel
  END FUNCTION particle_velocity

  ! mue-Dm relation of raindrops
! Cray compiler does not support OpenACC in elemental or pure
#ifdef _CRAYFTN
  FUNCTION rain_mue_dm_relation(this,D_m) result(mue)
#else
  PURE FUNCTION rain_mue_dm_relation(this,D_m) result(mue)
#endif

    !$ACC ROUTINE SEQ

    TYPE(particle_rain_coeffs), INTENT(in) :: this
    REAL(wp), INTENT(in) :: D_m
    REAL(wp)             :: mue, delta

    delta = this%cmu2*(D_m-this%cmu3)
    IF (D_m.LE.this%cmu3) THEN
      mue = this%cmu0*TANH((4.0_wp*delta)**2) + this%cmu4
    ELSE
      mue = this%cmu1*TANH(delta**2) + this%cmu4
    ENDIF
  END FUNCTION rain_mue_dm_relation

  !*******************************************************************************
  ! (2) More functions working on particle class, these are not CLASS procedures
  !*******************************************************************************

  ! bulk ventilation coefficient, Eq. (88) of SB2006
  REAL(wp) FUNCTION vent_coeff_a(parti,n)
    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: n
    CLASS(particle), INTENT(IN) :: parti

    vent_coeff_a = parti%a_ven * GAMMA((parti%nu+n+parti%b_geo)/parti%mu)                 &
         &                     / GAMMA((parti%nu+1.0_wp)/parti%mu)                        &
         &                   * ( GAMMA((parti%nu+1.0_wp)/parti%mu)                        &
         &                     / GAMMA((parti%nu+2.0_wp)/parti%mu) )**(parti%b_geo+n-1.0_wp)
  END FUNCTION vent_coeff_a

  ! bulk ventilation coefficient, Eq. (89) of SB2006
  REAL(wp) FUNCTION vent_coeff_b(parti,n)
    IMPLICIT NONE
    INTEGER, INTENT(in)         :: n
    CLASS(particle), INTENT(in) :: parti

    REAL(wp), PARAMETER :: m_f = 0.500 ! see PK, S.541. Do not change.

    vent_coeff_b = parti%b_ven                                                  &
         & * GAMMA((parti%nu+n+(m_f+1.0_wp)*parti%b_geo+m_f*parti%b_vel)/parti%mu)  &
         &             / GAMMA((parti%nu+1.0_wp)/parti%mu)                          &
         &           * ( GAMMA((parti%nu+1.0_wp)/parti%mu)                          &
         &             / GAMMA((parti%nu+2.0_wp)/parti%mu)                          &
         &             )**((m_f+1.0_wp)*parti%b_geo+m_f*parti%b_vel+n-1.0_wp)
  END FUNCTION vent_coeff_b

  ! complete mass moment of particle size distribution, Eq (82) of SB2006
  REAL(wp) FUNCTION moment_gamma(p,n)
    IMPLICIT NONE
    INTEGER, INTENT(in)           :: n
    CLASS(particle), INTENT(in)   :: p

    moment_gamma  = GAMMA((n+p%nu+1.0_wp)/p%mu) / GAMMA((p%nu+1.0_wp)/p%mu)        &
         &      * ( GAMMA((  p%nu+1.0_wp)/p%mu) / GAMMA((p%nu+2.0_wp)/p%mu) )**n
  END FUNCTION moment_gamma

  ! fractional mass moment of particle size distribution, i.e., this is
  ! Eq (82) of SB2006 with a non-integer exponent fexp
  REAL(wp) FUNCTION fracmoment_gamma(p,fexp)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: fexp
    CLASS(particle), INTENT(in)   :: p

    fracmoment_gamma  = GAMMA((fexp+p%nu+1.0_wp)/p%mu) / GAMMA((p%nu+1.0_wp)/p%mu)        &
         &          * ( GAMMA((     p%nu+1.0_wp)/p%mu) / GAMMA((p%nu+2.0_wp)/p%mu) )**fexp
  END FUNCTION fracmoment_gamma

  ! coefficient for slope of PSD, i.e., for lambda in Eq. (80) of SB2006
  REAL(wp) FUNCTION lambda_gamma(p,x)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: x
    CLASS(particle), INTENT(in)   :: p

    lambda_gamma  = ( GAMMA((p%nu+1.0_wp)/p%mu) / GAMMA((p%nu+2.0_wp)/p%mu) * x)**(-p%mu)
  END FUNCTION lambda_gamma

  ! coefficient for general collision integral.
  ! This function contains a generalized replacement of Eq. (90) of SB2006,
  ! because the latter is only correct for n=0 and 1.
  ! Due to different numeric calculation order, there are slight
  ! changes of the resulting coefficients in the 10th significant digit compared to
  ! the original function.
  REAL(wp) FUNCTION coll_delta(p1,n)
    IMPLICIT NONE
    CLASS(particle), INTENT(in) :: p1
    INTEGER, INTENT(in)         :: n

    coll_delta = GAMMA((2.0_wp*p1%b_geo+p1%nu+1.0_wp+n)/p1%mu)    &
         &                     / GAMMA((p1%nu+1.0_wp+n)/p1%mu)    &
         &     * GAMMA((p1%nu+1.0_wp)/p1%mu)**(2.0_wp*p1%b_geo)   &
         &     / GAMMA((p1%nu+2.0_wp)/p1%mu)**(2.0_wp*p1%b_geo)
    RETURN
  END FUNCTION coll_delta


  ! wrapper for coll_delta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_delta_11(p1,p2,n)
    CLASS(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n
    coll_delta_11 = coll_delta(p1,n)
    RETURN
  END FUNCTION coll_delta_11

  ! wrapper for coll_delta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_delta_22(p1,p2,n)
    CLASS(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n
    coll_delta_22 = coll_delta(p2,n)
    RETURN
  END FUNCTION coll_delta_22

  ! coefficient for general collision integral.
  ! This function contains a generalized replacement of Eq. (91) of SB2006,
  ! because the latter is only correct for n=0 and 1.
  ! Due to different numeric calculation order, there are slight
  ! changes of the resulting coefficients in the 10th significant digit compared to
  ! the original function.
  REAL(wp) FUNCTION coll_delta_12(p1,p2,n)
    CLASS(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_delta_12 = 2.0_wp * GAMMA((p1%b_geo+p1%nu+1.0_wp)/p1%mu)     &
         &                 / GAMMA((p1%nu+1.0_wp)/p1%mu)              &
         &                 * GAMMA((p1%nu+1.0_wp)/p1%mu)**(p1%b_geo)  &
         &                 / GAMMA((p1%nu+2.0_wp)/p1%mu)**(p1%b_geo)  &
         &                 * GAMMA((p2%b_geo+p2%nu+1.0_wp+n)/p2%mu)   &
         &                 / GAMMA((p2%nu+1.0_wp+n)/p2%mu)            &
         &                 * GAMMA((p2%nu+1.0_wp)/p2%mu)**(p2%b_geo)  &
         &                 / GAMMA((p2%nu+2.0_wp)/p2%mu)**(p2%b_geo)
    RETURN
  END FUNCTION coll_delta_12

  ! coefficient for general collision integral, Eq. (92) of SB2006
  REAL(wp) FUNCTION coll_theta(p1,n)
    CLASS(particle), INTENT(in) :: p1
    INTEGER, INTENT(in)         :: n

    coll_theta = GAMMA((2.0_wp*p1%b_vel+2.0_wp*p1%b_geo+p1%nu+1.0_wp+n)/p1%mu)    &
         &                     / GAMMA((2.0_wp*p1%b_geo+p1%nu+1.0_wp+n)/p1%mu)    &
         &                     * GAMMA((p1%nu+1.0_wp)/p1%mu)**(2.0_wp*p1%b_vel)   &
         &                     / GAMMA((p1%nu+2.0_wp)/p1%mu)**(2.0_wp*p1%b_vel)
    RETURN
  END FUNCTION coll_theta

  ! wrapper for coll_theta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_theta_11(p1,p2,n)
    CLASS(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_theta_11 = coll_theta(p1,n)
    RETURN
  END FUNCTION coll_theta_11

  ! wrapper for coll_theta (unnecessary and unused argument p2, but do not remove this)
  REAL(wp) FUNCTION coll_theta_22(p1,p2,n)
    CLASS(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_theta_22 = coll_theta(p2,n)
    RETURN
  END FUNCTION coll_theta_22

  ! coefficient for general collision integral, Eq. (93) of SB2006
  REAL(wp) FUNCTION coll_theta_12(p1,p2,n)
    CLASS(particle), INTENT(in) :: p1,p2
    INTEGER, INTENT(in)         :: n

    coll_theta_12 = 2.0_wp * GAMMA((p1%b_vel+2.0_wp*p1%b_geo+p1%nu+1.0_wp)/p1%mu)   &
         &                 / GAMMA((2.0_wp*p1%b_geo+p1%nu+1.0_wp)/p1%mu)            &
         &                 * GAMMA((p1%nu+1.0_wp)/p1%mu)**(p1%b_vel)                &
         &                 / GAMMA((p1%nu+2.0_wp)/p1%mu)**(p1%b_vel)                &
         &                 * GAMMA((p2%b_vel+2.0_wp*p2%b_geo+p2%nu+1.0_wp+n)/p2%mu) &
         &                 / GAMMA((2.0_wp*p2%b_geo+p2%nu+1.0_wp+n)/p2%mu)          &
         &                 * GAMMA((p2%nu+1.0_wp)/p2%mu)**(p2%b_vel)                &
         &                 / GAMMA((p2%nu+2.0_wp)/p2%mu)**(p2%b_vel)
    RETURN
  END FUNCTION coll_theta_12

  !********************************************************************************
  !
  ! Subroutines for setting up constant coeffs for computing generalized
  ! partial collision integrals T_ab^(n,m) based on spherical geometric collection kernel K_ab
  ! for spectral moments of order n and m of the gen-gamma PSDs, which are defined as
  !
  ! T_ab^(n,m) = int_xua^xoa int_xub^xob xa^n xb^m K_ab fa fb dxa dxb
  !
  ! These are for the parameterzation of collisions between two
  ! hydrometeor classes, where only a spectral part of the first class collides with a spectral part of the
  ! second class. The variable part of the coefficients has to be calculated
  ! for each time step in the corresponding collision subroutine (e.g., hail_rain_riming).

  ! usable for the constant (fixed) part of \delta_{aa}^{n,m} and \delta_{bb}^{m,n},
  ! regardless of lower or upper truncation:
  REAL(wp) FUNCTION coll_delta_aa_pm_fix(pa,pb,n,m)

    CLASS(PARTICLE), INTENT(in) :: pa,pb
    INTEGER, INTENT(in)         :: n,m

    coll_delta_aa_pm_fix =  ( GAMMA((pa%nu+1.0_wp)/pa%mu) / &
                              GAMMA((pa%nu+2.0_wp)/pa%mu) )**(2.0_wp*pa%b_geo) / &
               ( GAMMA((pa%nu+n+1.0_wp)/pa%mu) * GAMMA((pb%nu+m+1.0_wp)/pb%mu) )

    RETURN
  END FUNCTION coll_delta_aa_pm_fix

  ! usable for the constant (fixed) part of \delta_{ab}^{n,m},
  ! regardless of lower or upper truncation:
  REAL(wp) FUNCTION coll_delta_ab_pm_fix(pa,pb,n,m)

    CLASS(PARTICLE), INTENT(in) :: pa,pb
    INTEGER, INTENT(in)         :: n,m

    coll_delta_ab_pm_fix =  2.0_wp * &
         ( GAMMA((pa%nu+1.0_wp)/pa%mu) / GAMMA((pa%nu+2.0_wp)/pa%mu) )**(pa%b_geo) * &
         ( GAMMA((pb%nu+1.0_wp)/pb%mu) / GAMMA((pb%nu+2.0_wp)/pb%mu) )**(pb%b_geo) / &
         ( GAMMA((pa%nu+n+1.0_wp)/pa%mu) * GAMMA((pb%nu+m+1.0_wp)/pb%mu) )

    RETURN
  END FUNCTION coll_delta_ab_pm_fix

  ! usable for the constant (fixed) part of \theta_{aa}^{n,m} and \theta_{bb}^{m,n}
  ! for the approximation of the charact. velocity difference, regardless of lower or upper truncation:
  REAL(wp) FUNCTION coll_theta_aa_pm_fix(pa)

    CLASS(PARTICLE), INTENT(in) :: pa

    coll_theta_aa_pm_fix =  ( GAMMA((pa%nu+1.0_wp)/pa%mu) / &
                              GAMMA((pa%nu+2.0_wp)/pa%mu) )**(2.0_wp*pa%b_vel)

    RETURN
  END FUNCTION coll_theta_aa_pm_fix

  ! usable for the constant (fixed) part of \delta_{ab}^{n,m}
  ! for the approximation of the charact. velocity difference, regardless of lower or upper truncation:
  REAL(wp) FUNCTION coll_theta_ab_pm_fix(pa,pb)

    CLASS(PARTICLE), INTENT(in) :: pa,pb

    coll_theta_ab_pm_fix =  2.0_wp * &
         ( GAMMA((pa%nu+1.0_wp)/pa%mu) / GAMMA((pa%nu+2.0_wp)/pa%mu) )**(pa%b_vel) * &
         ( GAMMA((pb%nu+1.0_wp)/pb%mu) / GAMMA((pb%nu+2.0_wp)/pb%mu) )**(pb%b_vel)

    RETURN
  END FUNCTION coll_theta_ab_pm_fix


  ! usable for the constant (fixed) part of \delta_{aa}^{n,m},
  ! if integration over b is from 0 to infinity (full moment):
  REAL(wp) FUNCTION coll_delta_aa_pm_bfull_fix(pa,pb,n,m)

    CLASS(PARTICLE), INTENT(in) :: pa,pb
    INTEGER, INTENT(in)         :: n,m

    coll_delta_aa_pm_bfull_fix =  &
         ( GAMMA((pa%nu+1.0_wp)/pa%mu) / &
           GAMMA((pa%nu+2.0_wp)/pa%mu) )**(2.0_wp*pa%b_geo) / &
         GAMMA((pa%nu+n+1.0_wp)/pa%mu)

    RETURN
  END FUNCTION coll_delta_aa_pm_bfull_fix

  ! usable for the constant (fixed) part of \delta_{bb}^{n,m},
  ! if integration over b is from 0 to infinity (full moment):
  REAL(wp) FUNCTION coll_delta_bb_pm_bfull_fix(pa,pb,n,m)

    CLASS(PARTICLE), INTENT(in) :: pa,pb
    INTEGER, INTENT(in)         :: n,m

    coll_delta_bb_pm_bfull_fix =  &
         ( GAMMA((pa%nu+1.0_wp)/pa%mu) / &
           GAMMA((pa%nu+2.0_wp)/pa%mu) )**(2.0_wp*pa%b_geo) / &
         ( GAMMA((pa%nu+n+1.0_wp)/pa%mu) * GAMMA((pb%nu+m+1.0_wp)/pb%mu) ) * &
         GAMMA((pb%nu+2.0_wp*pb%b_geo+m+1.0_wp)/pb%mu)

    RETURN
  END FUNCTION coll_delta_bb_pm_bfull_fix

  ! usable for the constant (fixed) part of \delta_{ar}^{n,m},
  ! if integration over b is from 0 to infinity (full moment):
  REAL(wp) FUNCTION coll_delta_ab_pm_bfull_fix(pa,pb,n,m)

    CLASS(PARTICLE), INTENT(in) :: pa,pb
    INTEGER, INTENT(in)         :: n,m

    coll_delta_ab_pm_bfull_fix =  2.0_wp * &
         ( GAMMA((pa%nu+1.0_wp)/pa%mu) / GAMMA((pa%nu+2.0_wp)/pa%mu) )**(pa%b_geo) * &
         ( GAMMA((pb%nu+1.0_wp)/pb%mu) / GAMMA((pb%nu+2.0_wp)/pb%mu) )**(pb%b_geo) / &
         ( GAMMA((pa%nu+n+1.0_wp)/pa%mu) * GAMMA((pb%nu+m+1.0_wp)/pb%mu) ) * &
         GAMMA((pb%nu+pb%b_geo+m+1.0_wp)/pb%mu)

    RETURN
  END FUNCTION coll_delta_ab_pm_bfull_fix

  ! usable for the constant (fixed) part of \theta_{aa}^{n,m}
  ! for the approximation of the charact. velocity difference,
  ! if integration over b is from 0 to infinity (full moment):
  REAL(wp) FUNCTION coll_theta_aa_pm_bfull_fix(pa)

    CLASS(PARTICLE), INTENT(in) :: pa

    coll_theta_aa_pm_bfull_fix =  ( GAMMA((pa%nu+1.0_wp)/pa%mu) / &
                                    GAMMA((pa%nu+2.0_wp)/pa%mu) )**(2.0_wp*pa%b_vel)

    RETURN
  END FUNCTION coll_theta_aa_pm_bfull_fix

  ! usable for the constant (fixed) part of \theta_{bb}^{n,m}
  ! for the approximation of the charact. velocity difference,
  ! if integration over b is from 0 to infinity (full moment):
  REAL(wp) FUNCTION coll_theta_bb_pm_bfull_fix(pb,m)

    CLASS(PARTICLE), INTENT(in) :: pb
    INTEGER, INTENT(in)         :: m

    coll_theta_bb_pm_bfull_fix =  &
         ( GAMMA((pb%nu+1.0_wp)/pb%mu) / GAMMA((pb%nu+2.0_wp)/pb%mu) )**(2.0_wp*pb%b_vel) * &
           GAMMA((pb%nu+2.0_wp*pb%b_geo+2.0_wp*pb%b_vel+m+1.0_wp)/pb%mu) / &
           GAMMA((pb%nu+2.0_wp*pb%b_geo                +m+1.0_wp)/pb%mu)


    RETURN
  END FUNCTION coll_theta_bb_pm_bfull_fix

  ! usable for the constant (fixed) part of \theta_{ab}^{n,m}
  ! for the approximation of the charact. velocity difference,
  ! if integration over b is from 0 to infinity (full moment):
  REAL(wp) FUNCTION coll_theta_ab_pm_bfull_fix(pa,pb,m)

    CLASS(PARTICLE), INTENT(in) :: pa,pb
    INTEGER, INTENT(in)         :: m

    coll_theta_ab_pm_bfull_fix =  2.0_wp * &
         ( GAMMA((pa%nu+1.0_wp)/pa%mu) / GAMMA((pa%nu+2.0_wp)/pa%mu) )**(pa%b_vel) * &
         ( GAMMA((pb%nu+1.0_wp)/pb%mu) / GAMMA((pb%nu+2.0_wp)/pb%mu) )**(pb%b_vel) * &
           GAMMA((pb%nu+2.0_wp*pb%b_geo+pb%b_vel+m+1.0_wp)/pb%mu) / &
           GAMMA((pb%nu+2.0_wp*pb%b_geo         +m+1.0_wp)/pb%mu)

    RETURN
  END FUNCTION coll_theta_ab_pm_bfull_fix

  !==================================================================================================
  !
  ! Generic function for the often appearing form of the argument a of gamma functions in collision
  ! parameterizations
  !
  ! a = (c1*nu + c2*bgeo + c3*bvel + c4*n + 1) / mu
  !
  FUNCTION momarg_coll (p,n,c1,c2,c3,c4) RESULT (a)
    IMPLICIT NONE

    CLASS(PARTICLE), INTENT(in) :: p
    INTEGER, INTENT(in)         :: n
    REAL(wp), INTENT(in)        :: c1, c2, c3, c4
    REAL(wp)                    :: a

    a = ( c1*p%nu + c2*p%b_geo + c3*p%b_vel + c4*n + 1.0_wp) / p%mu

    RETURN
  END FUNCTION momarg_coll

  !==================================================================================================
  !
  ! arguments for the parameter a of gamma functions for use in the below collision parameterizations
  !
  ! a(1) = (nu + n + 1 ) / mu
  ! a(2) = (nu + bgeo + n + 1) / mu
  ! a(3) = (nu + 2*bgeo + n + 1) / mu
  ! a(4) = (nu + 2*bgeo + bvel + n + 1) / mu
  ! a(5) = (nu + 2*bgeo + 2*bvel + n + 1) / mu
  !
  FUNCTION momargs_coll_gam (p,n) RESULT (a)
    IMPLICIT NONE

    CLASS(PARTICLE), INTENT(in) :: p
    INTEGER, INTENT(in)         :: n
    REAL(wp)                    :: a(5)

    a(1) = momarg_coll (p,n,1.0_wp,0.0_wp,0.0_wp,1.0_wp)
    a(2) = momarg_coll (p,n,1.0_wp,1.0_wp,0.0_wp,1.0_wp)
    a(3) = momarg_coll (p,n,1.0_wp,2.0_wp,0.0_wp,1.0_wp)
    a(4) = momarg_coll (p,n,1.0_wp,2.0_wp,1.0_wp,1.0_wp)
    a(5) = momarg_coll (p,n,1.0_wp,2.0_wp,2.0_wp,1.0_wp)

    RETURN
  END FUNCTION momargs_coll_gam

  ! currently not used
  FUNCTION D_average_factor (parti)
    ! UB: Faktor zur Berechnung des mittleren Durchmessers von verallg. gammaverteilten Hydrometeoren:
    !     gueltig fuer D = a_geo * x^b_geo
    !     Berechnung des mittleren Durchmessers: D_average = parti%b_geo * D_average_factor * (q/qn)**parti%b_geo
    REAL(wp) :: D_average_factor
    CLASS(particle), INTENT(in) :: parti

    D_average_factor = &
         ( GAMMA( (parti%b_geo+parti%nu+1.0_wp)/parti%mu ) / &
           GAMMA( (parti%nu+1.0_wp)/parti%mu ) ) * &
         ( GAMMA( (parti%nu+1.0_wp)/parti%mu ) / GAMMA( (parti%nu+2.0_wp)/parti%mu ) ) ** parti%b_geo
  END FUNCTION D_average_factor

  !*******************************************************************************
  ! Fundamental physical relations
  !*******************************************************************************

  ! Molecular diffusivity of water vapor
  ELEMENTAL FUNCTION diffusivity(T,p) result(D_v)

    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN) :: T,p
    REAL(wp) :: D_v
    ! This is D_v = 8.7602e-5_wp * T_a**(1.81_wp) / p_a
    D_v = 8.7602e-5_wp * EXP(1.81_wp*LOG(T)) / p
    RETURN
  END FUNCTION diffusivity

  !*******************************************************************************
  ! saturation pressure over ice and liquid water                                *
  !*******************************************************************************

  ELEMENTAL REAL(wp) FUNCTION e_es (ta)
    !$ACC ROUTINE SEQ
    REAL(wp), INTENT(IN) :: ta
    e_es  = sat_pres_ice(ta)
    !e_es  = sat_pres_ice(ta,0)
  END FUNCTION e_es

  !  ELEMENTAL REAL(wp) FUNCTION e_ws (ta)
  !    REAL(wp), INTENT (IN) :: ta
  !    e_ws  = e_3 * EXP (A_w * (ta - T_3) / (ta - B_w))
  !  END FUNCTION e_ws_old

  ! FUNCTION e_ws_vec (ta,idim,jdim)
  !   INTEGER :: idim, jdim
  !   REAL(wp)               :: e_ws_vec(idim,jdim)
  !   REAL(wp), INTENT (IN)  :: ta(idim,jdim)
  !   e_ws_vec  = e_3 * EXP (A_w * (ta - T_3) / (ta - B_w))
  ! END FUNCTION e_ws_vec

  ! FUNCTION e_es_vec (ta,idim,jdim)
  !   INTEGER :: idim, jdim
  !   REAL(wp)               :: e_es_vec(idim,jdim)
  !   REAL(wp), INTENT (IN)  :: ta(idim,jdim)
  !   e_es_vec  = e_3 * EXP (A_e * (ta - T_3) / (ta - B_e))
  ! END FUNCTION e_es_vec


  !********************************************************************************
  ! setup subroutines
  !********************************************************************************

  SUBROUTINE setup_particle_coeffs(ptype,pcoeffs)
    CLASS(particle),        INTENT(in)    :: ptype
    CLASS(particle_coeffs), INTENT(inout) :: pcoeffs

    pcoeffs%c_i = 1.0 / ptype%cap
    pcoeffs%a_f = vent_coeff_a(ptype,1)
    pcoeffs%b_f = vent_coeff_b(ptype,1) * N_sc**n_f / sqrt(nu_l)
    pcoeffs%c_z = moment_gamma(ptype,2)

  END SUBROUTINE setup_particle_coeffs

  SUBROUTINE setup_cloud_autoconversion_sb(cloud,cloud_coeffs)
    CLASS(particle), INTENT(in) :: cloud
    TYPE(particle_cloud_coeffs), INTENT(inout) :: cloud_coeffs
    REAL(wp) :: nu, mu
    REAL(wp), PARAMETER :: kc_autocon  = 9.44e+9_wp  !..Long-Kernel

    nu = cloud%nu
    mu = cloud%mu
    IF (mu == 1.0) THEN
      !.. see SB2001
      cloud_coeffs%k_au  = kc_autocon / cloud%x_max * (1.0_wp / 20.0_wp) &
           &               * (nu+2.0_wp)*(nu+4.0_wp)/(nu+1.0_wp)**2
      cloud_coeffs%k_sc  = kc_autocon * (nu+2.0_wp)/(nu+1.0_wp)
    ELSE
      !.. see Eq. (3.44) of Seifert (2002)
      cloud_coeffs%k_au = kc_autocon / cloud%x_max * (1.0_wp / 20.0_wp)  &
           & * ( 2.0_wp * GAMMA((nu+4.0_wp)/mu)**1                           &
           &            * GAMMA((nu+2.0_wp)/mu)**1 * GAMMA((nu+1.0_wp)/mu)**2    &
           &   - 1.0_wp * GAMMA((nu+3.0_wp)/mu)**2 * GAMMA((nu+1.0_wp)/mu)**2 )  &
           &   / GAMMA((nu+2.0_wp)/mu)**4
      cloud_coeffs%k_sc = kc_autocon * cloud_coeffs%c_z
    ENDIF

  END SUBROUTINE setup_cloud_autoconversion_sb

  SUBROUTINE setup_ice_selfcollection(ice,ice_coeffs)
    CLASS(particle), INTENT(in)              :: ice
    TYPE(particle_ice_coeffs), INTENT(inout) :: ice_coeffs

    ! local variables
    REAL(wp) :: delta_n_11,delta_n_12,delta_n_22
    REAL(wp) :: delta_q_11,delta_q_12,delta_q_22
    REAL(wp) :: theta_n_11,theta_n_12,theta_n_22
    REAL(wp) :: theta_q_11,theta_q_12,theta_q_22

    CHARACTER(len=*), PARAMETER :: sroutine = 'setup_ice_selfcollection'

    delta_n_11 = coll_delta_11(ice,ice,0)
    delta_n_12 = coll_delta_12(ice,ice,0)
    delta_n_22 = coll_delta_22(ice,ice,0)
    delta_q_11 = coll_delta_11(ice,ice,0)
    delta_q_12 = coll_delta_12(ice,ice,1)
    delta_q_22 = coll_delta_22(ice,ice,1)

    theta_n_11 = coll_theta_11(ice,ice,0)
    theta_n_12 = coll_theta_12(ice,ice,0)
    theta_n_22 = coll_theta_22(ice,ice,0)
    theta_q_11 = coll_theta_11(ice,ice,0)
    theta_q_12 = coll_theta_12(ice,ice,1)
    theta_q_22 = coll_theta_22(ice,ice,1)

    ice_coeffs%sc_delta_n = delta_n_11 + delta_n_12 + delta_n_22
    ice_coeffs%sc_delta_q = delta_q_11 + delta_q_12 + delta_q_22
    ice_coeffs%sc_theta_n = theta_n_11 - theta_n_12 + theta_n_22
    ice_coeffs%sc_theta_q = theta_q_11 - theta_q_12 + theta_q_22

    IF (isdebug) THEN
      WRITE(txt,'(A,ES14.7)') "    a_ice      = ",ice%a_geo ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    b_ice      = ",ice%b_geo ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    alf_ice    = ",ice%a_vel ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    bet_ice    = ",ice%b_vel ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_n_11 = ",delta_n_11 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_n_12 = ",delta_n_12 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_n_22 = ",delta_n_22 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n_11 = ",theta_n_11 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n_12 = ",theta_n_12 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n_22 = ",theta_n_22 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_q_11 = ",delta_q_11 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_q_12 = ",delta_q_12 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_q_22 = ",delta_q_22 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_q_11 = ",theta_q_11 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_q_12 = ",theta_q_12 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_q_22 = ",theta_q_22 ; CALL message(sroutine,TRIM(txt))
    END IF
    IF (isprint) THEN
      WRITE(txt,'(A,ES14.7)') "    delta_n    = ",ice_coeffs%sc_delta_n ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n    = ",ice_coeffs%sc_theta_n ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_q    = ",ice_coeffs%sc_delta_q ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_q    = ",ice_coeffs%sc_theta_q ; CALL message(sroutine,TRIM(txt))
    END IF
  END SUBROUTINE setup_ice_selfcollection

  SUBROUTINE setup_snow_selfcollection(snow, snow_coeffs)
    CLASS(particle), INTENT(in)               :: snow
    TYPE(particle_snow_coeffs), INTENT(inout) :: snow_coeffs

    REAL(wp) :: delta_n_11,delta_n_12
    REAL(wp) :: theta_n_11,theta_n_12

    CHARACTER(len=*), PARAMETER :: sroutine = 'setup_snow_selfcollection'

    delta_n_11 = coll_delta_11(snow,snow,0)
    delta_n_12 = coll_delta_12(snow,snow,0)
    theta_n_11 = coll_theta_11(snow,snow,0)
    theta_n_12 = coll_theta_12(snow,snow,0)

    snow_coeffs%sc_delta_n = (2.0*delta_n_11 + delta_n_12)
    snow_coeffs%sc_theta_n = (2.0*theta_n_11 - theta_n_12)

    IF (isdebug) THEN
      WRITE(txt,'(A,ES14.7)') "    a_snow     = ",snow%a_geo ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    b_snow     = ",snow%b_geo ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    alf_snow   = ",snow%a_vel ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    bet_snow   = ",snow%b_vel ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_n_11 = ",delta_n_11 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_n_12 = ",delta_n_12 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n_11 = ",theta_n_11 ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n_12 = ",theta_n_12 ; CALL message(sroutine,TRIM(txt))
    END IF
    IF (isprint) THEN
      WRITE(txt,'(A,ES14.7)') "    delta_n    = ",snow_coeffs%sc_delta_n ; CALL message(sroutine,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n    = ",snow_coeffs%sc_theta_n ; CALL message(sroutine,TRIM(txt))
    END IF
  END SUBROUTINE setup_snow_selfcollection

  SUBROUTINE setup_graupel_selfcollection(graupel,graupel_coeffs)
    CLASS(particle), INTENT(in) :: graupel
    TYPE(particle_graupel_coeffs)  :: graupel_coeffs
    REAL(wp) :: delta_n_11,delta_n_12
    REAL(wp) :: theta_n_11,theta_n_12
    REAL(wp) :: delta_n, theta_n

    CHARACTER(len=*), PARAMETER :: routi = 'setup_graupel_selfcollection'

    delta_n_11 = coll_delta_11(graupel,graupel,0)
    delta_n_12 = coll_delta_12(graupel,graupel,0)
    theta_n_11 = coll_theta_11(graupel,graupel,0)
    theta_n_12 = coll_theta_12(graupel,graupel,0)

    delta_n = (2.0*delta_n_11 + delta_n_12)
    theta_n = (2.0*theta_n_11 - theta_n_12)**0.5

    graupel_coeffs%sc_coll_n  = pi8 * delta_n * theta_n

    IF (isprint) THEN
      WRITE(txt,'(A,ES14.7)') "    delta_n_11 = ",delta_n_11 ; CALL message(routi,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_n_12 = ",delta_n_12 ; CALL message(routi,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    delta_n    = ",delta_n    ; CALL message(routi,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n_11 = ",theta_n_11 ; CALL message(routi,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n_12 = ",theta_n_12 ; CALL message(routi,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    theta_n    = ",theta_n    ; CALL message(routi,TRIM(txt))
      WRITE(txt,'(A,ES14.7)') "    coll_n     = ",graupel_coeffs%sc_coll_n; CALL message(routi,TRIM(txt))
    END IF
  END SUBROUTINE setup_graupel_selfcollection

  SUBROUTINE setup_particle_collection_type1(ptype,qtype,coll_coeffs)
    CLASS(particle), INTENT(in) :: ptype, qtype
    TYPE(collection_coeffs)     :: coll_coeffs
    CHARACTER(len=*), PARAMETER :: routi = 'setup_particle_collection_type1'

    coll_coeffs%delta_n_aa = coll_delta_11(ptype,qtype,0)
    coll_coeffs%delta_n_ab = coll_delta_12(ptype,qtype,0)
    coll_coeffs%delta_n_bb = coll_delta_22(ptype,qtype,0)
    coll_coeffs%delta_q_aa = coll_delta_11(ptype,qtype,0)
    coll_coeffs%delta_q_ab = coll_delta_12(ptype,qtype,1)
    coll_coeffs%delta_q_bb = coll_delta_22(ptype,qtype,1)

    coll_coeffs%theta_n_aa = coll_theta_11(ptype,qtype,0)
    coll_coeffs%theta_n_ab = coll_theta_12(ptype,qtype,0)
    coll_coeffs%theta_n_bb = coll_theta_22(ptype,qtype,0)
    coll_coeffs%theta_q_aa = coll_theta_11(ptype,qtype,0)
    coll_coeffs%theta_q_ab = coll_theta_12(ptype,qtype,1)
    coll_coeffs%theta_q_bb = coll_theta_22(ptype,qtype,1)

  END SUBROUTINE setup_particle_collection_type1

  SUBROUTINE setup_particle_collection_type2(ptype,qtype,coll_coeffs)
    CLASS(particle), INTENT(in) :: ptype, qtype
    TYPE(rain_riming_coeffs)    :: coll_coeffs
    CHARACTER(len=*), PARAMETER :: routi = 'setup_particle_collection_type2'

    coll_coeffs%delta_n_aa = coll_delta_11(ptype,qtype,0)
    coll_coeffs%delta_n_ab = coll_delta_12(ptype,qtype,0)
    coll_coeffs%delta_n_bb = coll_delta_22(ptype,qtype,0)
    coll_coeffs%delta_q_aa = coll_delta_11(ptype,qtype,1) ! mass weighted
    coll_coeffs%delta_q_ab = coll_delta_12(ptype,qtype,1)
    coll_coeffs%delta_q_ba = coll_delta_12(qtype,ptype,1)
    coll_coeffs%delta_q_bb = coll_delta_22(ptype,qtype,1)

    coll_coeffs%theta_n_aa = coll_theta_11(ptype,qtype,0)
    coll_coeffs%theta_n_ab = coll_theta_12(ptype,qtype,0)
    coll_coeffs%theta_n_bb = coll_theta_22(ptype,qtype,0)
    coll_coeffs%theta_q_aa = coll_theta_11(ptype,qtype,1) ! mass weighted
    coll_coeffs%theta_q_ab = coll_theta_12(ptype,qtype,1)
    coll_coeffs%theta_q_ba = coll_theta_12(qtype,ptype,1)
    coll_coeffs%theta_q_bb = coll_theta_22(ptype,qtype,1)

  END SUBROUTINE setup_particle_collection_type2

  SUBROUTINE setup_particle_coll_pm_type1(pa,pb,coeffs)
    CLASS(particle), INTENT(in)          :: pa, pb
    TYPE(coll_coeffs_ir_pm), INTENT(out) :: coeffs
    CHARACTER(len=*), PARAMETER :: routi = 'setup_particle_coll_pm_type1'

    INTEGER              :: i, j

    ! coll a+b->a with partial integration range for both species

    ! prepare 0 and first partial moments collision terms:
    DO i=0,1
      coeffs%moma(i,:) = momargs_coll_gam (pa,i)
      coeffs%momb(i,:) = momargs_coll_gam (pb,i)
    END DO

    DO i=0,1
      DO j=0,1
        coeffs%delta_aa(i,j) = coll_delta_aa_pm_fix(pa,pb,i,j)
        coeffs%delta_bb(i,j) = coll_delta_aa_pm_fix(pb,pa,j,i)
        coeffs%delta_ab(i,j) = coll_delta_ab_pm_fix(pa,pb,i,j)
      END DO
    END DO

    coeffs%theta_aa(:,:) = coll_theta_aa_pm_fix(pa)
    coeffs%theta_bb(:,:) = coll_theta_aa_pm_fix(pb)
    coeffs%theta_ab(:,:) = coll_theta_ab_pm_fix(pa,pb)

    coeffs%lamfakt_a = ( GAMMA(coeffs%moma(1,1)) / GAMMA(coeffs%moma(0,1)) )**(pa%mu)
    coeffs%lamfakt_b = ( GAMMA(coeffs%momb(1,1)) / GAMMA(coeffs%momb(0,1)) )**(pb%mu)

  END SUBROUTINE setup_particle_coll_pm_type1

  SUBROUTINE setup_particle_coll_pm_type1_bfull(pa,pb,coeffs)
    CLASS(particle), INTENT(in)          :: pa, pb
    TYPE(coll_coeffs_ir_pm), INTENT(out) :: coeffs
    CHARACTER(len=*), PARAMETER :: routi = 'setup_particle_coll_pm_type1'

    INTEGER              :: i, n, m

    ! coll a+b->a with partner b beeing integrated from 0 to infty (full moment)
    !             and a beeing a partial moment

    ! prepare 0 and first partial moments collision terms:
    DO i=0,1
      coeffs%moma(i,:) = momargs_coll_gam (pa,i)
      coeffs%momb(i,:) = momargs_coll_gam (pb,i)
    END DO

    DO n=0,1
      DO m=0,1
        coeffs%delta_aa(n,m) = coll_delta_aa_pm_bfull_fix(pa,pb,n,m)
        coeffs%delta_bb(n,m) = coll_delta_bb_pm_bfull_fix(pa,pb,n,m)
        coeffs%delta_ab(n,m) = coll_delta_ab_pm_bfull_fix(pa,pb,n,m)
      END DO
    END DO

    coeffs%theta_aa(:,:) = coll_theta_aa_pm_bfull_fix(pa)
    DO m=0,1
      coeffs%theta_bb(:,m) = coll_theta_bb_pm_bfull_fix(pb,m)
      coeffs%theta_ab(:,m) = coll_theta_ab_pm_bfull_fix(pa,pb,m)
    END DO

    coeffs%lamfakt_a = ( GAMMA(coeffs%moma(1,1)) / GAMMA(coeffs%moma(0,1)) )**(pa%mu)
    coeffs%lamfakt_b = ( GAMMA(coeffs%momb(1,1)) / GAMMA(coeffs%momb(0,1)) )**(pb%mu)

  END SUBROUTINE setup_particle_coll_pm_type1_bfull

  SUBROUTINE set_ccn_cloud_type (zccn_type, zcloud_type, zccn_coeffs, zN_cn0, zcfg_2mom)

    INTEGER, INTENT(out)                   :: zccn_type, zcloud_type
    TYPE(aerosol_ccn), INTENT(out)         :: zccn_coeffs
    REAL(wp), OPTIONAL, INTENT(in)         :: zN_cn0
    TYPE(t_cfg_2mom), OPTIONAL, INTENT(in) :: zcfg_2mom

    CHARACTER(len=*), PARAMETER :: routine = 'set_ccn_cloud_type'

    IF (PRESENT(zN_cn0)) THEN
      IF (PRESENT(zcfg_2mom)) THEN
        IF (zcfg_2mom%ccn_type > 0) THEN
          zccn_type   = zcfg_2mom%ccn_type
        ELSE
          zccn_type   = ccn_type_gscp5
        END IF
      ELSE
        zccn_type   = ccn_type_gscp5
      END IF
      zcloud_type = cloud_type_default_gscp5 + 10 * zccn_type
    ELSE
      IF (PRESENT(zcfg_2mom)) THEN
        IF (zcfg_2mom%ccn_type > 0) THEN
          zccn_type   = zcfg_2mom%ccn_type
        ELSE
          zccn_type   = ccn_type_gscp4
        END IF
      ELSE
        zccn_type   = ccn_type_gscp4
      END IF
      zcloud_type = cloud_type_default_gscp4 + 10 * zccn_type
    END IF

    !..parameters for exponential decrease of N_ccn with height
    !  z0:  up to this height (m) constant unchanged value
    !  z1e: height interval at which N_ccn decreases by factor 1/e above z0_nccn

    zccn_coeffs%z0  = 4000.0_wp
    zccn_coeffs%z1e = 2000.0_wp

    ! min updraft speed for Segal&Khain activation
    zccn_coeffs%wcb_min = cfg_params%ccn_wcb_min

    ! characteristics of different kinds of CN
    ! (copied from COSMO 5.0 Segal & Khain nucleation subroutine)

    SELECT CASE(zccn_type)
    CASE(6)
      !... maritime case
      zccn_coeffs%Ncn0 = 100.0e6_wp   ! CN concentration at ground
      zccn_coeffs%Nmin =  35.0e6_wp   ! NOT relevant at the moment
      zccn_coeffs%lsigs = 0.4_wp      ! log(sigma_s)
      zccn_coeffs%R2    = 0.03_wp     ! in mum
      zccn_coeffs%etas  = 0.9_wp      ! soluble fraction
    CASE(7)
      !... intermediate case
      zccn_coeffs%Ncn0 = 250.0e6_wp
      zccn_coeffs%Nmin =  35.0e6_wp
      zccn_coeffs%lsigs = 0.4_wp
      zccn_coeffs%R2    = 0.03_wp       ! in mum
      zccn_coeffs%etas  = 0.8_wp        ! soluble fraction
    CASE(8)
      !... continental case
      zccn_coeffs%Ncn0 = 1700.0e6_wp
      zccn_coeffs%Nmin =   35.0e6_wp  ! NOT relevant at the moment
      zccn_coeffs%lsigs = 0.2_wp
      zccn_coeffs%R2    = 0.03_wp     ! in mum
      zccn_coeffs%etas  = 0.7_wp      ! soluble fraction
    CASE(9)
      !... "polluted" continental
      zccn_coeffs%Ncn0 = 3200.0e6_wp
      zccn_coeffs%Nmin =   35.0e6_wp    ! NOT relevant at the moment
      zccn_coeffs%lsigs = 0.2_wp
      zccn_coeffs%R2    = 0.03_wp       ! in mum
      zccn_coeffs%etas  = 0.7_wp        ! soluble fraction
    CASE(1)
       !... dummy values
       zccn_coeffs%Ncn0  =  200.0e6_wp
       zccn_coeffs%Nmin  =   10.0e6_wp  ! NOT relevant at the moment
       zccn_coeffs%lsigs = 0.0_wp
       zccn_coeffs%R2    = 0.0_wp
       zccn_coeffs%etas  = 0.0_wp
    CASE DEFAULT
       CALL finish(TRIM(routine),'Error: Invalid value for ccn_type')
    END SELECT

  END SUBROUTINE set_ccn_cloud_type

END MODULE mo_2mom_mcrph_setup
