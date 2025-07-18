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

!> TERRA heat and moisture transfer process descriptions.
MODULE sfc_terra_transport

  USE mo_kind, ONLY: wp

  USE mo_exception, ONLY: finish
  USE mo_physical_constants, ONLY: &
      & g => grav, &
      & lh_f => alf, &
      & rho_w => rhoh2o, &
      & t0_melt => tmelt

  USE mo_lnd_nwp_config, ONLY: &
      & cwimax_ml, &
      & itype_eisa, &
      & itype_evsl, &
      & itype_heatcond, &
      & itype_hydbound, &
      & lmulti_snow, &
      & lterra_urb

  USE sfc_terra_data, ONLY: &
      & b_clay, &
      & b_org, &
      & b_sand, &
      & b_silt, &
      & cadp, &
      & cala0, &
      & cala1, &
      & cclayf, &
      & cdw0, &
      & cdw1, &
      & cfcap, &
      & ckw0, &
      & ckw1, &
      & cporv, &
      & cpwp, &
      & csandf, &
      & ctau_i, &
      & cwisamax, &
      & eps_div, &
      & eps_temp, &
      & eps_soil, &
      & EVSL_RESIST, &
      & EVSL_RESIST_RBS, &
      & HCOND_AVG, &
      & HCOND_PL98, &
      & HCOND_PL98_VEG, &
      & IST_PEAT, &
      & IST_ROCK, &
      & IST_SAND, &
      & IST_SLOAM, &
      & itype_mire, &
      & T_ref_ice, &
      & T_star_ice

  USE sfc_terra_snow, ONLY: snow_multi_calc_heat_conduction
  USE sfc_terra_util, ONLY: solve_tridiag, zalfa

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_hydrology
  PUBLIC :: calc_infiltration
  PUBLIC :: calc_heat_conductivity
  PUBLIC :: calc_heat_conduction
  PUBLIC :: calc_soil_water_melt

  !> Number of layers contributing to surface runoff.
  INTEGER, PARAMETER :: SURFACE_RUNOFF_LAYERS = 1

  !> Number of hydrologically active soil layers in mires.
  INTEGER, PARAMETER :: ke_soil_hy_m = 4

  ! Silence unused variable warnings on non-OpenACC builds
#ifdef _OPENACC
# define OPENACC_SUPPRESS_UNUSED_LZACC
#else
# define OPENACC_SUPPRESS_UNUSED_LZACC IF (lzacc .AND. acc_async_queue > 0) THEN; END IF
#endif

CONTAINS

SUBROUTINE calc_hydrology ( &
      & dt, ke_soil, ke_soil_hy, nvec, ivstart, &
      & ivend, soiltyp_subs, dz_hl, dz_ml, infil_rate, runoff_grav, w_so_new, w_so_now, &
      & w_so_ice_now, runoff_s, runoff_g, eva_bs, transp_ml, dt_w_so, root_depth, z_ml, z_hl, hydiffu_fac, &
      & lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN)  :: ke_soil !< Number of soil moisture layers.
  INTEGER, INTENT(IN)  :: ke_soil_hy !< Number of active soil moisture layers.
  INTEGER, INTENT(IN)  :: nvec !< array dimensions
  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program

  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< type of the soil (keys 0-9)
  REAL(wp), INTENT(IN) :: dz_hl(ke_soil+1) !< Layer interface distance (n_soil) [m].
  REAL(wp), INTENT(IN) :: dz_ml(ke_soil+1) !< Main level distance (n_soil) [m].

  REAL(wp), INTENT(IN) :: dt !< integration time-step [s].

  REAL(wp), INTENT(IN) :: infil_rate(nvec) !< infiltration rate [kg/(m^2 s)].

  REAL(wp), INTENT(INOUT) :: runoff_grav(nvec,ke_soil+1) !< main level water gravitation flux (for runoff calc.)
  REAL(wp), INTENT(INOUT) :: w_so_new(nvec,ke_soil+1) !< new total water content (ice + liquid water) [m H2O]

  REAL(wp), INTENT(IN) :: w_so_now(nvec,ke_soil+1) !< current total water content (ice + liquid water) [m H2O]
  REAL(wp), INTENT(IN) :: w_so_ice_now(nvec,ke_soil+1) !< current ice content [m H2O]

  REAL(wp), INTENT(INOUT) :: runoff_s(nvec) !< surface water runoff; sum over forecast [kg/m2]
  REAL(wp), INTENT(INOUT) :: runoff_g(nvec) !< soil water runoff; sum over forecast [kg/m2]

  REAL(wp), INTENT(IN) :: eva_bs(nvec) !< evaporation from bare soil
  REAL(wp), INTENT(IN) :: transp_ml(nvec,ke_soil) !< transpiration contribution by the different layers

  REAL(wp), INTENT(INOUT) :: dt_w_so(nvec,ke_soil) !< tendency of water content [kg/(m**3 s)]

  REAL(wp), INTENT(IN) :: root_depth(nvec) !< Root depth [m].
  REAL(wp), INTENT(IN) :: z_ml(ke_soil+1) !< Layer center depth [m].
  REAL(wp), INTENT(IN) :: z_hl(ke_soil+1) !< Half-level depth [m]
  REAL(wp), INTENT(IN) :: hydiffu_fac(nvec) ! tuning factor for hydraulic diffusivity

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  ! local variables
  REAL(wp) :: dt_o_rho_w !< timestep / density of water [m^3 s/kg].
  REAL(wp) :: rho_w_o_dt !< density of water / timestep [kg/(m^3 s)].

    ! Water transport
  REAL(wp) :: zfr_ice !< reduction factor for water transport
  REAL(wp) :: fr_ice_ksom1 !< fractional ice content of actual layer - 1
  REAL(wp) :: fr_ice_kso !< fractional ice content of actual layer
  REAL(wp) :: fr_ice_ksop1 !< fractional ice content of actual layer + 1
  REAL(wp) :: fr_liq_ksom1 !< fractional liquid water content of actual layer - 1
  REAL(wp) :: fr_liq_ksom1_new !< fractional liquid water content of actual layer -1
  REAL(wp) :: fr_liq_kso !< fractional liquid water content of actual layer
  REAL(wp) :: fr_liq_kso_new !< fractional liquid water content of actual layer
  REAL(wp) :: fr_liq_ksop1 !< fractional liquid water content of actual layer + 1
  REAL(wp) :: fr_liq_ksom05 !< fractional liquid water content of actual layer - 1/2
  REAL(wp) :: fr_liq_ksop05 !< fractional liquid water content of actual layer + 1/2
  REAL(wp) :: zdlw_fr_ksom05 !< hydraulic diffusivity coefficient at half level above
  REAL(wp) :: zklw_fr_ksom05 !< hydraulic conductivity coefficient at half level above
  REAL(wp) :: zklw_fr_kso_new !< hydraulic conductivity coefficient at main level for actual runoff_g
  REAL(wp) :: zdlw_fr_kso !< hydraulic diff coefficient at main level
  REAL(wp) :: zklw_fr_kso !< hydraulic conductivity coefficient at main level
  REAL(wp) :: zklw_fr_ksom1 !< hydraulic conductivity coefficient at main level above
  REAL(wp) :: zdlw_fr_ksop05 !< hydraulic diffusivity coefficient at half level below
  REAL(wp) :: zklw_fr_ksop05 !< hydraulic conductivity coefficient at half level below
  REAL(wp) :: zcou_roffg !< indicator to sum up runoffg in hydrological active layers

  ! Implicit solution of hydraulic equation
  REAL(wp) :: z1dgam1 !< utility variable
  REAL(wp) :: zredm !< utility variable
  REAL(wp) :: zredm05 !< utility variable
  REAL(wp) :: zredp05 !< utility variable
  REAL(wp) :: zgam2m05 !< utility variable
  REAL(wp) :: zgam2p05 !< utility variable

  ! Hydraulic parameters
  REAL(wp) :: zro_sfak !< utility variable
  REAL(wp) :: zro_gfak !< utility variable
  REAL(wp) :: zfmb_fak !< utility variable
  REAL(wp) :: zdwg !< preliminary change of soil water content
  REAL(wp) :: zwgn !< preliminary change of soil water content
  REAL(wp) :: zredfu !< utility variable for runoff determination
  REAL(wp) :: zro !< utility variable for runoff determination
  REAL(wp) :: zro2 !< utility variable for runoff determination
  REAL(wp) :: zkorr !< utility variable for runoff determination

  ! Utility variables
  REAL(wp) :: lse_a(nvec,ke_soil_hy) !< LSE subdiagonal
  REAL(wp) :: lse_b(nvec,ke_soil_hy) !< LSE diagonal
  REAL(wp) :: lse_c(nvec,ke_soil_hy) !< LSE superdiagonal
  REAL(wp) :: lse_rhs(nvec,ke_soil_hy) !< LSE right-hand side
  REAL(wp) :: lse_sol(nvec,ke_soil_hy) !< LSE solution

  LOGICAL :: has_soil(nvec) !< Solver mask.
  LOGICAL :: mire_mask(nvec) !< Solver mask for mires.

  ! ground water as lower boundary of soil column
  REAL(wp) :: zdelta_sm
  REAL(wp) :: zdhydcond_dlwfr

  INTEGER :: mstyp
  INTEGER :: kso_end
  ! Indices
  INTEGER :: kso !< loop index for soil moisture layers
  INTEGER :: i !< loop index in x-direction

  REAL(wp) :: fr_liq_ml(nvec,ke_soil+1) !< fractional liqu. water content of soil layer
  REAL(wp) :: fr_ice_ml(nvec,ke_soil+1) !< fractional ice content of soil layer
  REAL(wp) :: fr_w_ml(nvec,ke_soil+1) !< fractional total water content of soil layers

  REAL(wp) :: pore_vol(nvec)
  REAL(wp) :: air_dryness_point(nvec)
  REAL(wp) :: zdw(nvec)
  REAL(wp) :: zdw1(nvec)
  REAL(wp) :: zkw(nvec)
  REAL(wp) :: zkw1(nvec)

  REAL(wp) :: zdw_mod !< Diffusivity modified by adaptive tuning [m^2/s]

  REAL(wp) :: zkw_m05 !< Conductivity at saturation across upper interface accounting for root depth [m/s]
  REAL(wp) :: zkw_p05 !< Conductivity at saturation across lower interface accounting for root depth [m/s]
  REAL(wp) :: zkw_ml !< Conductivity at saturation across layer accounting for root depth [m/s]

  ! end of definitions

  OPENACC_SUPPRESS_UNUSED_LZACC

  dt_o_rho_w = dt / rho_w
  rho_w_o_dt = rho_w / dt

  !$ACC DATA CREATE(mire_mask) ASYNC(acc_async_queue) IF(lzacc .AND. itype_mire == 1)

  !$ACC DATA CREATE(fr_w_ml, fr_liq_ml, fr_ice_ml, pore_vol, air_dryness_point, zdw, zdw1, zkw, zkw1) &
  !$ACC   CREATE(lse_a, lse_b, lse_c, lse_rhs, lse_sol, has_soil) NO_CREATE(mire_mask) PRESENT(ivend) &
  !$ACC   ASYNC(acc_async_queue)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP SEQ
  DO   kso = 1,ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      fr_w_ml (i,kso) = w_so_now(i,kso)/dz_hl(kso)
      fr_ice_ml(i,kso) = w_so_ice_now(i,kso)/dz_hl(kso)   ! ice frac.
      fr_liq_ml(i,kso) = fr_w_ml(i,kso) - fr_ice_ml(i,kso)  ! liquid water frac.

      runoff_grav(i,kso) = 0._wp
    END DO
  END DO      !soil layers

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp)
  DO i = ivstart, ivend
    mstyp = soiltyp_subs(i)
    pore_vol(i) = cporv(mstyp)
    air_dryness_point(i) = cadp(mstyp)
    zdw(i) = cdw0(mstyp)
    zdw1(i) = cdw1(mstyp)
    zkw(i) = ckw0(mstyp)
    zkw1(i) = ckw1(mstyp)

    has_soil(i) = (mstyp > IST_ROCK)
  END DO

  ! uppermost layer, kso = 1

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fr_ice_kso, fr_ice_ksop1, fr_liq_kso) &
  !$ACC   PRIVATE(fr_liq_ksop1, zdw_mod, zfr_ice, zredp05, fr_liq_ksop05) &
  !$ACC   PRIVATE(zdlw_fr_ksop05, zklw_fr_ksop05, z1dgam1, zgam2p05, zkw_p05)
  DO i = ivstart, ivend
    IF (has_soil(i)) THEN
      ! sedimentation and capillary transport in soil
      ! Note: The fractional liquid water content (concentration)  of each layer
      !       is normalized by the ice free fraction of each layer in order to
      !       obtain a representative concentration of liquid water in the
      !       'active' part of each soil layer
      !       Hydraulic diffusivity and conductivity coefficients are multiplied
      !       by a reduction factor depending on the maximum ice fraction of the
      !       adjacent layers in order to avoid the transport of liquid water
      !       in to the frozen part of the adjacent layer
      !
      ! GZ, 2017-11-07:
      ! I wonder if this is sufficient to describe the reduction of water transport in frozen soils.
      ! Certainly, the hydraulic conductivity does not drop abruptly to zero once first ice crystals form in the soil,
      ! but if rain falls after a longer-term frost period without snow cover, ponds form even if the soil surface
      ! has already thawed. Should the hydraulic conductivity include a term like ...*(1-tuning_factor*fr_ice_ml/zporv)?
      !
      fr_ice_kso   = fr_ice_ml(i,1)
      fr_ice_ksop1 = fr_ice_ml(i,2)
      fr_liq_kso   = fr_liq_ml(i,1)
      fr_liq_ksop1 = fr_liq_ml(i,2)

      ! adaptive parameter tuning for near-surface hydraulic diffusivity
      zdw_mod = zdw(i)*hydiffu_fac(i)**2

      !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
      zkw_p05 = zkw(i) * EXP(-2._wp*(z_hl(1)-root_depth(i)))

      ! compute reduction factor for transport coefficients
      zfr_ice  = max (fr_ice_kso,fr_ice_ksop1)
      zredp05  = 1._wp-zfr_ice/MAX(fr_liq_kso+fr_ice_kso,fr_liq_ksop1+fr_ice_ksop1)

      ! interpolated scaled liquid water fraction at layer interface
      fr_liq_ksop05  = 0.5_wp*(dz_hl(2)*fr_liq_kso+dz_hl(1)*fr_liq_ksop1) / dz_ml(2)

      zdlw_fr_ksop05= zredp05 * watrdiff_RT(zdw_mod, fr_liq_ksop05, zdw1(i), pore_vol(i), air_dryness_point(i))
      zklw_fr_ksop05= zredp05 * watrcon_RT (zkw_p05, fr_liq_ksop05, zkw1(i), pore_vol(i), air_dryness_point(i))

      ! coefficients for implicit flux computation
      z1dgam1     = dt/dz_hl(1)
      zgam2p05    = zdlw_fr_ksop05/dz_ml(2)
      lse_a(i,1) = 0._wp
      lse_b(i,1) = 1._wp+zalfa*zgam2p05*z1dgam1
      lse_c(i,1) = -zalfa * zgam2p05*z1dgam1
      lse_rhs(i,1) = fr_liq_ml(i,1) + infil_rate(i)*z1dgam1/rho_w  &
                         -zklw_fr_ksop05*z1dgam1                     &
                         +(1._wp - zalfa)* zgam2p05*z1dgam1*(fr_liq_ksop1 - fr_liq_kso)  &
                         +                 zgam2p05*z1dgam1*(fr_ice_ksop1-fr_ice_kso)

    END IF
  END DO

! inner layers 2 <=kso<=ke_soil_hy-1

  !$ACC LOOP SEQ
!DIR$ ivdep, prefervector
  DO kso =2,ke_soil_hy-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fr_ice_ksom1, fr_ice_kso) &
    !$ACC   PRIVATE(fr_ice_ksop1, fr_liq_ksom1, fr_liq_kso, fr_liq_ksop1) &
    !$ACC   PRIVATE(fr_liq_ksom05, fr_liq_ksop05, zfr_ice, zredm05, zredp05) &
    !$ACC   PRIVATE(zdlw_fr_ksom05, zdlw_fr_ksop05, zklw_fr_ksom05, zdw_mod) &
    !$ACC   PRIVATE(zklw_fr_ksop05, z1dgam1, zgam2m05, zgam2p05, zkw_m05, zkw_p05)
    DO i = ivstart, ivend
      IF (has_soil(i)) THEN
        ! sedimentation and capillary transport in soil
        fr_ice_ksom1 = fr_ice_ml(i,kso-1)
        fr_ice_kso   = fr_ice_ml(i,kso  )
        fr_ice_ksop1 = fr_ice_ml(i,kso+1)
        fr_liq_ksom1 = fr_liq_ml(i,kso-1)
        fr_liq_kso   = fr_liq_ml(i,kso  )
        fr_liq_ksop1 = fr_liq_ml(i,kso+1)

        ! adaptive parameter tuning for near-surface hydraulic diffusivity
        zdw_mod = MERGE(zdw(i)*hydiffu_fac(i), zdw(i), kso == 2)

        !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
        zkw_m05 = zkw(i) * EXP(-2._wp*(z_hl(kso-1)-root_depth(i)))
        zkw_p05 = zkw(i) * EXP(-2._wp*(z_hl(kso)-root_depth(i)))

        ! interpolated scaled liquid water content at interface to layer
        ! above and below
        fr_liq_ksom05 = 0.5_wp*(dz_hl(kso-1)*fr_liq_kso+   &
                               dz_hl(kso)*fr_liq_ksom1)/dz_ml(kso)
        fr_liq_ksop05 = 0.5_wp*(dz_hl(kso+1)*fr_liq_kso+   &
                               dz_hl(kso)*fr_liq_ksop1)/dz_ml(kso+1)

        ! compute reduction factor for coefficients
        zfr_ice          = max (fr_ice_kso,fr_ice_ksom1)
        zredm05 = 1._wp-zfr_ice/max (fr_liq_kso+fr_ice_kso,fr_liq_ksom1+fr_ice_ksom1)
        zfr_ice          = max (fr_ice_kso,fr_ice_ksop1)
        zredp05 = 1._wp-zfr_ice/max (fr_liq_kso+fr_ice_kso,fr_liq_ksop1+fr_ice_ksop1)


        zdlw_fr_ksom05= zredm05*watrdiff_RT(zdw_mod,fr_liq_ksom05,zdw1(i),pore_vol(i),air_dryness_point(i))
        zdlw_fr_ksop05= zredp05*watrdiff_RT(zdw_mod,fr_liq_ksop05,zdw1(i),pore_vol(i),air_dryness_point(i))
        zklw_fr_ksom05= zredm05*watrcon_RT(zkw_m05,fr_liq_ksom05,zkw1(i),pore_vol(i),air_dryness_point(i))
        zklw_fr_ksop05= zredp05*watrcon_RT(zkw_p05,fr_liq_ksop05,zkw1(i),pore_vol(i),air_dryness_point(i))


        ! coefficients for implicit flux computation
        z1dgam1 = dt/dz_hl(kso)
        zgam2m05  = zdlw_fr_ksom05/dz_ml(kso)
        zgam2p05  = zdlw_fr_ksop05/dz_ml(kso+1)
        lse_a (i,kso) = -zalfa*zgam2m05*z1dgam1
        lse_c (i,kso) = -zalfa*zgam2p05*z1dgam1
        lse_b (i,kso) = 1._wp +zalfa*(zgam2m05+zgam2p05)*z1dgam1
        lse_rhs (i,kso) = fr_liq_ml(i,kso)+                               &
                              z1dgam1*(-zklw_fr_ksop05+zklw_fr_ksom05)+ &
                              (1._wp-zalfa)*z1dgam1*                &
                              (zgam2p05*(fr_liq_ksop1-fr_liq_kso  )     &
                              -zgam2m05*(fr_liq_kso  -fr_liq_ksom1)   ) &
                             +z1dgam1*                                  &
                              (zgam2p05*(fr_ice_ksop1-fr_ice_kso  )   &
                              -zgam2m05*(fr_ice_kso-fr_ice_ksom1)   )

      END IF
    END DO
  END DO
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fr_ice_ksom1, fr_ice_kso, fr_liq_ksom1) &
  !$ACC   PRIVATE(fr_liq_kso, fr_liq_ksom05, zfr_ice, zredm05, zdlw_fr_ksom05) &
  !$ACC   PRIVATE(z1dgam1, zgam2m05, zklw_fr_ksom05, zkw_m05, kso, zdw_mod)
  DO i = ivstart, ivend
    IF (has_soil(i)) THEN
      IF (soiltyp_subs(i) == IST_PEAT .AND. itype_mire == 1) THEN
        kso = ke_soil_hy_m
      ELSE
        kso = ke_soil_hy
      END IF

      ! lowest active hydrological layer ke_soil_hy{,_m}-1
      fr_ice_ksom1  = fr_ice_ml(i,kso-1)
      fr_ice_kso    = fr_ice_ml(i,kso  )
      fr_liq_ksom1  = fr_liq_ml(i,kso-1)
      fr_liq_kso    = fr_liq_ml(i,kso  )
      fr_liq_ksom05 = 0.5_wp*(dz_hl(kso-1)*fr_liq_kso+ &
                            dz_hl(kso)*fr_liq_ksom1)/dz_ml(kso)

      ! adaptive parameter tuning for near-surface hydraulic diffusivity
      zdw_mod = MERGE(zdw(i)*hydiffu_fac(i), zdw(i), kso == 2)

      !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
      zkw_m05 = zkw(i) * EXP(-2._wp*(z_hl(kso-1)-root_depth(i)))

      zfr_ice          = max (fr_ice_kso,fr_ice_ksom1)
      zredm05 = 1._wp-zfr_ice/max (fr_liq_kso+fr_ice_kso,fr_liq_ksom1+fr_ice_ksom1)

      zdlw_fr_ksom05 = zredm05*watrdiff_RT(zdw_mod,fr_liq_ksom05,zdw1(i),pore_vol(i),air_dryness_point(i))
      zklw_fr_ksom05 = zredm05*watrcon_RT(zkw_m05,fr_liq_ksom05,zkw1(i),pore_vol(i),air_dryness_point(i))

      z1dgam1 = dt/dz_hl(kso)
      zgam2m05  = zdlw_fr_ksom05/dz_ml(kso)

      lse_a(i,kso) = -zalfa* zgam2m05*z1dgam1
      lse_b(i,kso) = 1._wp+ zalfa*zgam2m05*z1dgam1
      lse_c(i,kso) = 0.0_wp
      lse_rhs(i,kso) = fr_liq_kso+z1dgam1*zklw_fr_ksom05               &
                              +(1._wp-zalfa)*z1dgam1*               &
                              zgam2m05*(fr_liq_ksom1  - fr_liq_kso) &
                              +z1dgam1*                             &
                              zgam2m05*(fr_ice_ksom1-fr_ice_kso )
    END IF
  END DO

  !solve tridiagonal matrix
  IF (itype_mire == 1) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      mire_mask(i) = (soiltyp_subs(i) == IST_PEAT)
    END DO

    ! Solve peat cells down to ke_soil_hy_m.
    CALL solve_tridiag ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & nlev=ke_soil_hy_m, &
        & a=lse_a, &
        & b=lse_b, &
        & c=lse_c, &
        & d=lse_rhs, &
        & out=lse_sol, &
        & mask=mire_mask &
      )

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      mire_mask(i) = (has_soil(i) .AND. (soiltyp_subs(i) /= IST_PEAT))
    END DO

    ! Solve the full column for everything else.
    CALL solve_tridiag ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & nlev=ke_soil_hy, &
        & a=lse_a, &
        & b=lse_b, &
        & c=lse_c, &
        & d=lse_rhs, &
        & out=lse_sol, &
        & mask=mire_mask &
      )

    !$ACC LOOP SEQ
    DO kso = 1, ke_soil_hy
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF ((has_soil(i) .AND. soiltyp_subs(i) /= IST_PEAT) &
            & .OR. (soiltyp_subs(i) == IST_PEAT .AND. kso <= ke_soil_hy_m)) THEN
          w_so_new(i,kso) = lse_sol(i,kso)*dz_hl(kso) + w_so_ice_now(i,kso)
        END IF
      END DO
    END DO

  ELSE

    CALL solve_tridiag ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & nlev=ke_soil_hy, &
        & a=lse_a, &
        & b=lse_b, &
        & c=lse_c, &
        & d=lse_rhs, &
        & out=lse_sol, &
        & mask=has_soil &
      )

    !$ACC LOOP SEQ
    DO kso = 1, ke_soil_hy
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (has_soil(i)) THEN
          w_so_new(i,kso) = lse_sol(i,kso)*dz_hl(kso) + w_so_ice_now(i,kso)
        END IF
      END DO
    END DO
  END IF

  ! to ensure vertical constant water concentration profile beginning at
  ! layer ke_soil_hy{,_m} for energetic treatment only
  ! soil water climate layer(s)

  ! J. Helmert: Activate the ground water table with water diffusion from layers beyond ke_soil_hy
  !             for mire points.
  !             This should avoid a dry falling of the upper mire layers due to strong bare soil evaporation!

  ! ground water as lower boundary of soil column

  IF (itype_hydbound == 3) THEN
    !$ACC LOOP SEQ
    DO kso = ke_soil_hy+1,ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (has_soil(i)) THEN
          w_so_new(i,kso) = pore_vol(i)*dz_hl(kso)
        END IF
      END DO
    END DO
  ELSE
    !$ACC LOOP SEQ
    DO kso = ke_soil_hy+1,ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (has_soil(i)) THEN
          w_so_new(i,kso) = w_so_new(i,kso-1)*dz_hl(kso)/dz_hl(kso-1)
        END IF
      END DO
    END DO
  ENDIF

  IF (itype_mire == 1) THEN
    !$ACC LOOP SEQ
    DO kso = ke_soil_hy_m+1,ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        IF (soiltyp_subs(i) == IST_PEAT) THEN
          w_so_new(i,kso) = pore_vol(i)*dz_hl(kso)
        END IF
      END DO
    END DO
  END IF

  ! combine implicit part of sedimentation and capillary flux with explicit part
  ! (for soil water flux investigations only)
  !$ACC LOOP SEQ
  DO kso = 2,ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR &
    !$ACC   PRIVATE(fr_ice_ksom1, fr_ice_kso, fr_liq_ksom1_new) &
    !$ACC   PRIVATE(fr_liq_kso_new, fr_liq_ksom1, fr_liq_kso, zfr_ice, zredm05) &
    !$ACC   PRIVATE(fr_liq_ksom05, zdlw_fr_ksom05, zklw_fr_ksom05, zredm) &
    !$ACC   PRIVATE(zklw_fr_kso_new, zdelta_sm, zdlw_fr_kso, zklw_fr_kso) &
    !$ACC   PRIVATE(zklw_fr_ksom1, zdhydcond_dlwfr, zkw_ml, zkw_m05, zkw_p05, zdw_mod)
    DO i = ivstart, ivend
      IF (has_soil(i)) THEN
        fr_ice_ksom1 = fr_ice_ml(i,kso-1)
        fr_ice_kso   = fr_ice_ml(i,kso)
        fr_liq_ksom1 = fr_liq_ml(i,kso-1)
        fr_liq_kso   = fr_liq_ml(i,kso  )
        fr_liq_ksom1_new= w_so_new(i,kso-1)/dz_hl(kso-1) - fr_ice_ksom1
        fr_liq_kso_new  = w_so_new(i,kso  )/dz_hl(kso  ) - fr_ice_kso

        ! adaptive parameter tuning for near-surface hydraulic diffusivity
        zdw_mod = MERGE(zdw(i)*hydiffu_fac(i), zdw(i), kso == 2)

        !fc=2 1/m Exponential Ksat-profile decay parameter,see Decharme et al. (2006)
        zkw_m05 = zkw(i) * EXP(-2._wp*(z_hl(kso-1)-root_depth(i)))
        zkw_ml = zkw(i) * EXP(-2._wp*(z_ml(kso) - root_depth(i)))
        zkw_p05 = zkw(i) * EXP(-2._wp*(z_hl(kso)-root_depth(i)))

        !... additionally for runoff_g at lower level of lowest active water
        ! layer calculated with (upstream) main level soil water content
        ! compute reduction factor for transport coefficients
        zfr_ice          = max (fr_ice_kso,fr_ice_ksom1)
        zredm05 = 1._wp-zfr_ice/max (fr_liq_kso+fr_ice_kso,fr_liq_ksom1+fr_ice_ksom1)

        ! interpolated liquid water content at interface to layer above
        fr_liq_ksom05 =0.5_wp*(dz_hl(kso)*fr_liq_ksom1+dz_hl(kso-1)*fr_liq_kso) /dz_ml(kso)

        zdlw_fr_ksom05 = zredm05*watrdiff_RT(zdw_mod,fr_liq_ksom05,zdw1(i),pore_vol(i),air_dryness_point(i))
        zklw_fr_ksom05 = zredm05*watrcon_RT(zkw_m05,fr_liq_ksom05,zkw1(i),pore_vol(i),air_dryness_point(i))

        IF ((soiltyp_subs(i) == IST_PEAT .AND. itype_mire == 1 .AND. kso > ke_soil_hy_m) .OR. &
            (.NOT. (soiltyp_subs(i) == IST_PEAT .AND. itype_mire == 1) .AND. kso > ke_soil_hy)) THEN
          zdlw_fr_ksom05=0.0_wp   ! no flux gradient contribution below 2.5m
          zklw_fr_ksom05=0.0_wp   ! no gravitation flux below 2.5m
        END IF

        zredm = 1._wp-fr_ice_kso/(fr_liq_kso+fr_ice_kso)

        zklw_fr_kso_new = zredm*watrcon_RT(zkw_ml,fr_liq_kso_new,zkw1(i),pore_vol(i),air_dryness_point(i))

        ! actual gravitation water flux
        IF (w_so_new(i,kso).LT.1.01_wp*air_dryness_point(i)*dz_hl(kso)) THEN
          zklw_fr_kso_new=0._wp
        ENDIF
        runoff_grav(i,kso) =  - rho_w * zklw_fr_kso_new

        ! ground water as lower boundary of soil column
        IF (((itype_mire == 1 .AND. soiltyp_subs(i) == IST_PEAT) .AND. (kso == ke_soil_hy_m+1)) .OR. &
            (kso == ke_soil_hy+1 .AND. itype_hydbound == 3)) THEN
          zdelta_sm=( fr_liq_kso_new - fr_liq_ksom1_new )

          ! GZ, 2017-11-07: The variable staggering in the following expressions needs to be checked
          !   zredm05 is used for interface levels kso and kso-1; logically, "zredp05" would be
          !   needed in the first two expressions in addition, the indexing of zkw1 would be
          !   incorrect if this parameter is made level-dependent
          zdlw_fr_kso = zredm05*watrdiff_RT(zdw_mod,fr_liq_kso_new,    &
                          zdw1(i),pore_vol(i),air_dryness_point(i))
          zklw_fr_kso = zredm05*watrcon_RT(zkw_p05,                    &
                          fr_liq_kso_new,zkw1(i),pore_vol(i),air_dryness_point(i))
          zklw_fr_ksom1 = zredm05*watrcon_RT(zkw_m05,&
                          fr_liq_ksom1_new,zkw1(i),pore_vol(i),air_dryness_point(i))

          zdhydcond_dlwfr=( zklw_fr_kso - zklw_fr_ksom1 ) / zdelta_sm
          runoff_grav(i,kso-1)=runoff_grav(i,kso-1)+ zdhydcond_dlwfr / &
             (1.0_wp-EXP(-zdhydcond_dlwfr/zdlw_fr_kso*0.5_wp*dz_ml(kso)))* zdelta_sm
        ENDIF
      END IF
    END DO
  END DO

#ifdef __INTEL_COMPILER
!DIR$ NOFMA
#endif
  !$ACC LOOP SEQ PRIVATE(zro_sfak, zro_gfak)
!$NEC novector
  DO  kso = 1,ke_soil
    ! utility variables used to avoid if-constructs in following loops
    zro_sfak = MERGE(1._wp, 0._wp, SURFACE_RUNOFF_LAYERS > kso) ! 1.0 for 'surface runoff'
    zro_gfak = 1._wp - zro_sfak                                 ! 1.0 for 'ground runoff'

    ! - Compute subsoil runoff (runoff_g) as drainage flux through bottom
    !   of layer ke_soil_hy (as suggested by the Rhone-Aggregation
    !   Experiment)
    ! - soil moisture gradient related flux is switched off below
    !   (i.e. only sedimentation flux allowed between ke_soil_hy and ke_soil_hy+1)

    ! sedimentation and capillary transport in soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdwg, zredfu, zro, zwgn, zro2, zkorr, mstyp) &
    !$ACC   PRIVATE(kso_end, zfmb_fak, zcou_roffg)
    DO i = ivstart, ivend
      IF (has_soil(i)) THEN
        mstyp = soiltyp_subs(i)

        kso_end = MERGE(ke_soil_hy_m, ke_soil_hy, itype_mire == 1 .AND. mstyp == IST_PEAT)
        zfmb_fak = MERGE(1.0_wp, 0.0_wp, kso==kso_end)
        zcou_roffg = MERGE(1.0_wp, 0.0_wp, kso<=kso_end)

        ! first runoff calculation without consideration of
        ! evapotranspiration
        !zdwg calculated above by flux divergence has to be aequivalent with
        zdwg =  (w_so_new(i,kso)/dz_hl(kso)-fr_w_ml(i,kso))*dz_hl(kso) / dt_o_rho_w
        zdwg =  zdwg + runoff_grav(i,kso)*zfmb_fak
        zredfu =  MAX( 0.0_wp, MIN( 1.0_wp,(fr_w_ml(i,kso) -     &
                          cfcap(mstyp))/MAX(pore_vol(i) - cfcap(mstyp),eps_div)) )
        zredfu = MERGE(zredfu, 0._wp, zdwg >= 0._wp)
        zro    = zdwg*zredfu
        zdwg   = zdwg*(1._wp - zredfu)

        ! add evaporation (first layer only)
        ! and transpiration (for each layer)
        zdwg   = zdwg + MERGE(eva_bs(i), 0._wp, kso == 1) + transp_ml (i,kso)
        zwgn   = fr_w_ml(i,kso) + dt_o_rho_w*zdwg/dz_hl(kso)
        zro2   = rho_w_o_dt*dz_hl(kso)*MAX(0.0_wp, zwgn - pore_vol(i))
        zkorr  = rho_w_o_dt*dz_hl(kso)*MAX(0.0_wp, air_dryness_point(i) - zwgn )
        dt_w_so(i,kso)= zdwg + zkorr - zro2
        zro    = zro      + zro2
        runoff_s(i) = runoff_s(i) + zro*zro_sfak*dt
        ! only count runoff_g in the hydrological active layers
        runoff_g(i) = runoff_g(i) + zcou_roffg*zro*zro_gfak*dt
        ! runoff_g reformulation due to drainage flux through bottom of layer ke_soil_hy
        !  zfmb_fak is only 1 at the last active horizont and 0 elsewhere
        !  only at this level the gravitational settling has to be counted for runoff (bulk view)
        runoff_g(i) = runoff_g(i) - (runoff_grav(i,kso) * zfmb_fak &
                                          + zkorr) * dt
      END IF
    END DO
  END DO         ! end loop over soil layers
  !$ACC END PARALLEL
  !$ACC END DATA
  !$ACC END DATA

END SUBROUTINE calc_hydrology


SUBROUTINE calc_infiltration ( &
      & nvec, ivstart, ivend, dt, soiltyp_subs, w_i_now, sp_10m, rain_rate, rain_rate_con, snow_rate, &
      & conv_frac, ice_rate, qc_atm, qi_atm, evapotrans_snfr, evapo_snow, infil_rate, &
      & runoff_s, fr_snow, w_snow_now, rho_ch, eva_bs, transp_sum, plcov, tai, urb_isa, dew_rate, &
      & rime_rate, rain_dew_rate, snow_rime_rate, dt_w_i, dt_w_snow, t_sk_now, t_s_now, w_i_max, &
      & ldiff_qi, ldepo_qw, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN)  :: nvec !< array dimensions
  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program


  ! Timestep parameters
  REAL(wp), INTENT(IN) :: dt !< timestep [s]

  INTEGER, INTENT(IN) :: soiltyp_subs(nvec)

  REAL(wp), INTENT(IN) :: transp_sum(nvec) !< total transpiration (transpiration from all soil layers) [kg/(m^2 s)]
  REAL(wp), INTENT(IN) :: w_i_now(nvec) !< water cont. of interception store  [m H20]
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< snow water equivalent [m H20]

  ! Connection to the atmosphere
  REAL(wp), INTENT(IN) :: eva_bs(nvec) !< evaporation from bare soil [kg/(m^2 s)]
  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< snow-cover fraction [1]
  REAL(wp), INTENT(IN) :: rho_ch(nvec) !< transfer coefficient*rho*g

  ! Thermal variables
  REAL(wp), INTENT(IN) :: t_sk_now(nvec) !< Skin temperature [K]
  REAL(wp), INTENT(IN) :: t_s_now(nvec) !< soil surface temperature [K]

  REAL(wp), INTENT(IN) :: urb_isa(nvec) !< Urban impervious surface area index [m^2/m^2].

    !   Interception variables
  REAL(wp), INTENT(IN) :: dew_rate(nvec) !< dew-formation rate [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: rime_rate(nvec) !< rime-formation rate [kg/(m^2 s)].
  REAL(wp), INTENT(INOUT) :: rain_dew_rate(nvec) !< total rain rate including formation of dew [kg/(m^2 s)]
  REAL(wp), INTENT(INOUT) :: snow_rime_rate(nvec) !< total snow rate including formation of rime [kg/(m^2 s)]
  REAL(wp), INTENT(INOUT) :: infil_rate(nvec) !< infiltration rate [kg/(m^2 s)]
  REAL(wp), INTENT(INOUT) :: w_i_max(nvec) !< maximum interception store [m H2O]
  REAL(wp), INTENT(INOUT) :: evapotrans_snfr(nvec) !< total evapotranspiration [kg/(m^2 s)]
  REAL(wp), INTENT(INOUT) :: evapo_snow(nvec) !< total evaporation from snow surface [kg/(m^2 s)]

  ! Tendencies
  REAL(wp), INTENT(INOUT) :: dt_w_i(nvec) !< tendency of water content of interception store [kg/(m^2 s)]
  REAL(wp), INTENT(INOUT) :: dt_w_snow(nvec) !< tendency of snow water content [kg/(m^2 s)]

  ! Soil and plant parameters
  REAL(wp), INTENT(IN) :: plcov(nvec) !< fraction of plant cover [m^2/m^2]
  REAL(wp), INTENT(IN) :: tai(nvec) !< transpiration area index [m^2/m^2]

  REAL(wp), INTENT(INOUT) :: runoff_s(nvec) !< surface water runoff; sum over forecast [kg/m^2]

  REAL(wp), INTENT(IN) :: sp_10m(nvec) !< wind speed in 10m [m/s]
  REAL(wp), INTENT(IN) :: rain_rate(nvec) !< Total precipitation rate of rain [kg/(m^2 s)]
  REAL(wp), INTENT(IN) :: rain_rate_con(nvec) !< precipitation rate of rain, convective [kg/(m^2 s)]
  REAL(wp), INTENT(IN) :: snow_rate(nvec) !< Total precipitation rate of snow, excluding ice [kg/(m^2 s)]
  REAL(wp), INTENT(IN) :: conv_frac(nvec) !< convective area fraction as assumed in convection scheme
  REAL(wp), INTENT(IN) :: ice_rate(nvec) !< precipitation rate of ice, grid-scale [kg/(m^2 s)]
  REAL(wp), INTENT(IN) :: qc_atm(nvec) !< specific liquid-water content [kg/kg]
  REAL(wp), INTENT(IN) :: qi_atm(nvec) !< specific frozen-water content [kg/kg]

  LOGICAL, INTENT(IN) :: ldiff_qi !< Turbulent diffusion of frozen water is active.
  LOGICAL, INTENT(IN) :: ldepo_qw !< Deposition of (frozen or liquid) cloud water required.
  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  ! local variables

  INTEGER :: i !< loop index in x-direction

  REAL(wp) :: zalf !< utility variable
  REAL(wp) :: zro_inf !< surface runoff
  REAL(wp) :: dt_w_snow_prov !< provisional time rate of change of snow store
  REAL(wp) :: w_i_prov !< preliminary value of interception store [m H2O]
  REAL(wp) :: w_snow_prov !< preliminary value of snow store
  REAL(wp) :: dt_w_i_prov !< provisional time rate of change of interception store
  REAL(wp) :: zdwieps !< artificial change of small amount of interception store
  REAL(wp) :: zdwseps !< artificial change of small amount of snow store
  REAL(wp) :: zro_wi !< surface runoff due to limited infiltration capacity
  REAL(wp) :: zinf !< infiltration
  REAL(wp) :: zrime !< ground riming rate
  REAL(wp) :: zzz !< utility variable
  REAL(wp) :: ztau_i
  REAL(wp) :: ztfunc

  REAL(wp) :: dt_o_rho_w !< timestep / density of water [s m^3/kg]
  REAL(wp) :: rho_w_o_dt !< density of water / timestep [kg/(m^3 s)]

  REAL(wp) :: zvers(nvec) !< water supply for infiltration

  ! end of definitions

  OPENACC_SUPPRESS_UNUSED_LZACC

  ! time constant for infiltration of water from interception store
  ! must not be less than 2*time step
  ztau_i = MAX( ctau_i, dt )
  dt_o_rho_w = dt / rho_w
  rho_w_o_dt = rho_w / dt

  !$ACC DATA PRESENT(ivend) CREATE(zvers) ASYNC(acc_async_queue)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrime, zalf, zzz, zinf, dt_w_i_prov, zdwieps, zro_wi) &
  !$ACC   PRIVATE(zro_inf, dt_w_snow_prov, w_snow_prov, zdwseps, ztfunc)
  DO i = ivstart, ivend
    ! store forcing terms due to evapotranspiration, formation of dew
    ! and rime for later use
    evapotrans_snfr(i) = dt_w_i(i) + eva_bs(i) + transp_sum(i) + &
                    (1._wp-fr_snow(i))*(dew_rate(i) + rime_rate(i))
    evapo_snow(i) = dt_w_snow(i) + fr_snow(i)*(dew_rate(i) + rime_rate(i))

    ! add grid scale and convective precipitation (and graupel, if present)
    ! to dew and rime
    rain_dew_rate(i) = dew_rate(i) + rain_rate(i)
    zrime = rime_rate(i)
    snow_rime_rate(i) = rime_rate(i) + snow_rate(i) + ice_rate(i)

    ! add possible deposition of water droplets or ice particles:
    IF (ldepo_qw) THEN !cloud-particle deposition related to vertical diffsuion with zero-concentration condition
      rain_dew_rate(i) = rain_dew_rate(i) + rho_ch(i)*qc_atm(i) !dew -fall + liquid precipitation including deposition
      IF (ldiff_qi) snow_rime_rate(i) = snow_rime_rate(i) + rho_ch(i)*qi_atm(i) !rime-fall + frozn. precipitation including deposition
    END IF

    ! Decide whether riming is added to interception store or to snow cover
    IF (snow_rime_rate(i) >= 1.05_wp*zrime .OR. fr_snow(i) >= 0.9_wp) zrime = 0._wp
    snow_rime_rate(i) = snow_rime_rate(i) - zrime

    ! infiltration and surface run-off

    ! subtract evaporation from interception store to avoid negative
    ! values due to sum of evaporation+infiltration
    w_i_prov = w_i_now(i) + dt_w_i(i)*dt_o_rho_w
    w_i_prov = MAX(0.0_wp,w_i_prov)

    ! Linear transition between 0C (ztfunc=1) and 2C (ztfunc=0)
    ztfunc = MAX(0.0_wp,1.0_wp - MAX(0.0_wp,0.5_wp*(t_sk_now(i)-t0_melt)))
    w_i_max(i) = cwimax_ml*(1.0_wp+ztfunc)*MAX(ztfunc, eps_soil, MAX(2.5_wp*plcov(i),tai(i)))

    ! TERRA_URB: Puddles on the impervious surface area
    IF (lterra_urb .AND. (itype_eisa == 3)) THEN
      w_i_max(i) = urb_isa(i) * cwisamax + (1.0_wp - urb_isa(i)) * w_i_max(i)
    END IF

    zalf   = SQRT(MAX(0.0_wp,1.0_wp - w_i_prov/w_i_max(i)))

    ! water supply from interception store (if Ts above freezing)
    zzz    = MAX(0.1_wp, 0.4_wp - 0.05_wp*sp_10m(i))
    zinf   = MAX(0._wp,w_i_prov-zzz*w_i_max(i))*rho_w/ztau_i*        &
             (1._wp+0.75_wp*MAX(0._wp,sp_10m(i)-1._wp))*(1._wp-ztfunc)**2

    ! possible contribution of rain to infiltration
    !    IF (rain_dew_rate(i)-eps_soil > 0.0_wp) THEN
    !      zalf = MAX( zalf,                                                   &
    !             (rho_w_o_dt*MAX(0.0_wp, w_i_max(i)-zwinstr) + zinf)/rain_dew_rate(i) )
    !      zalf = MAX( 0.01_wp, MIN(1.0_wp, zalf) )

    ! if rain falls onto snow, all rain is considered for infiltration
    ! as no liquid water store is considered in the snowpack
    IF (w_snow_now(i) > 0.0_wp) zalf = 0.0_wp
      ! Increase infiltration and reduce surface runoff (bugfix)
      ! this deactivates filling the interception store!!

    ! rain freezes on the snow surface
    IF (lmulti_snow .AND. w_snow_now(i) > 0.0_wp) zalf = 1.0_wp

    ! interception store; convective precip is taken into account with a
    ! fractional area passed from the convection scheme
    dt_w_i_prov = zalf*(rain_dew_rate(i)+(conv_frac(i)-1._wp)*rain_rate_con(i)) + zrime + dt_w_i(i)-zinf
    w_i_prov = w_i_now(i) + dt_w_i_prov*dt_o_rho_w
    w_i_prov = MAX(0.0_wp, w_i_prov) !avoid negative values (security)
    zdwieps = 0.0_wp
    IF (w_i_prov > 0.0_wp .AND. w_i_prov < 1.0E-4_wp*eps_soil) THEN
      zdwieps    = w_i_prov*rho_w_o_dt
      runoff_s(i)= runoff_s(i) + zdwieps*dt
      dt_w_i_prov    = dt_w_i_prov - zdwieps
      w_i_prov    = 0.0_wp
    END IF

    ! add excess over w_i_max(i) to infiltration
    zro_wi       = rho_w_o_dt*MAX( 0.0_wp, w_i_prov-w_i_max(i) )
    dt_w_i_prov  = dt_w_i_prov - zro_wi
    dt_w_i(i)    = dt_w_i_prov
    zinf         = zinf + zro_wi
    IF (t_s_now(i) <= t0_melt) THEN
      ! add excess rime to snow
      snow_rime_rate(i) = snow_rime_rate(i) + zinf
      zinf = 0.0_wp
    ENDIF

    ! add rain contribution to water supply for infiltration
    ! surface runoff is evaluated after the calculation of infiltration
    zvers(i) = zinf + (1._wp - zalf)*rain_dew_rate(i) + (1._wp-conv_frac(i))*zalf*rain_rate_con(i)

    ! Avoid infiltration for rock, ice and snow-covered surfaces
    infil_rate(i) = zvers(i)*MERGE(1._wp, 0._wp, soiltyp_subs(i) > IST_ROCK)*(1.0_wp - fr_snow(i))

    ! Avoid infiltration for urban impervious surface area fraction
    IF (lterra_urb .AND. ((itype_eisa == 2) .OR. (itype_eisa == 3))) THEN
      infil_rate(i) = infil_rate(i) * (1.0_wp - urb_isa(i))
    END IF

    ! Add difference to surface runoff
    zro_inf = zvers(i) - infil_rate(i)

    runoff_s(i) = runoff_s(i) + zro_inf*dt

    ! change of snow water and interception water store
    ! (negligible residuals are added to the run-off)

    ! snow store
    dt_w_snow_prov = snow_rime_rate(i) + dt_w_snow(i)
    w_snow_prov  = w_snow_now(i) + dt_w_snow_prov*dt_o_rho_w
    w_snow_prov  = MAX(0.0_wp, w_snow_prov) ! avoid negative values (security)
    zdwseps  = 0.0_wp
    IF (w_snow_prov > 0.0_wp .AND. w_snow_prov < eps_soil) THEN
      ! shift marginal snow amounts to interception storage
!      IF (ztsnow_pm(i) > 0.0_wp) THEN
        zdwseps    = w_snow_prov*rho_w_o_dt
!        runoff_s(i) = runoff_s(i) + zdwseps*dt   ! previous implementation
        dt_w_snow_prov   = dt_w_snow_prov - zdwseps
        dt_w_i(i)  = dt_w_i(i) + zdwseps
!      END IF
    END IF
    dt_w_snow(i) = dt_w_snow_prov
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE calc_infiltration


SUBROUTINE calc_heat_conductivity ( &
      & ivstart, ivend, nvec, ke_soil, soiltyp_subs, z_ml, z_hl, dz_hl, t_so_now, w_so_now, &
      & w_so_ice_now, plcov, root_depth, heatcond_fac, t_snred, z0, urb_h_bld, urb_ai, urb_isa, &
      & urb_hcon, hcond_ml, hcond_hl, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart !< Start index for computations in the parallel program.
  INTEGER, INTENT(IN) :: ivend !< End index for computations in the parallel program.
  INTEGER, INTENT(IN) :: nvec !< Block length.
  INTEGER, INTENT(IN) :: ke_soil !< Number of active soil layers.

  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type.

  REAL(wp), INTENT(IN) :: z_ml(ke_soil+1) !< Soil layer depth [m].
  REAL(wp), INTENT(IN) :: z_hl(ke_soil) !< Lower soil interface depth [m].
  REAL(wp), INTENT(IN) :: dz_hl(ke_soil+1) !< Soil layer thickness [m].
  REAL(wp), INTENT(IN) :: t_so_now(nvec,0:ke_soil+1) !< Soil temperature [K].
  REAL(wp), INTENT(IN) :: w_so_now(nvec,ke_soil+1) !< Liquid soil water equivalent [m].
  REAL(wp), INTENT(IN) :: w_so_ice_now(nvec,ke_soil+1) !< Frozen soil water equivalent [m].
  REAL(wp), INTENT(IN) :: plcov(nvec) !< Plant-cover fraction [m**2/m**2].
  REAL(wp), INTENT(IN) :: root_depth(nvec) !< Depth of root zone [m].
  REAL(wp), INTENT(IN) :: heatcond_fac(nvec) !< Heat conductivity tuning factor [1].
  REAL(wp), INTENT(IN) :: t_snred(nvec) !< Snow evaporation temperature offset [K].
  REAL(wp), INTENT(IN) :: z0(nvec) !< Roughness length [m].
  REAL(wp), INTENT(IN) :: urb_h_bld(nvec) !< Urban tile building height [m].
  REAL(wp), INTENT(IN) :: urb_ai(nvec) !< Urban tile building area index [m**2/m**2].
  REAL(wp), INTENT(IN) :: urb_isa(nvec) !< Urban tile fraction [m**2/m**2].
  REAL(wp), INTENT(IN) :: urb_hcon(nvec) !< Urban tile building heat conductivity [W/(K m)].
  REAL(wp), INTENT(INOUT) :: hcond_ml(nvec,ke_soil+1) !< Heat conductivity in layers [W/(K m)].
  REAL(wp), INTENT(INOUT) :: hcond_hl(nvec,ke_soil) !< Heat coductivity across half levels [W/(K m)]

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp), PARAMETER :: rho_solids = 2700.0_wp !< Density of solid soil parts [kg/m**3].
  REAL(wp), PARAMETER :: lam_water_ln = LOG(0.57_wp) !< LOG(thermal conductivity of water)
  REAL(wp), PARAMETER :: lam_ice_ln = LOG(2.2_wp) !< LOG(thermal conductivity of ice)
  REAL(wp), PARAMETER :: lam_quartz_ln = LOG(7.7_wp) !< LOG(thermal conductivity of quartz)
  REAL(wp), PARAMETER :: ln_2 = LOG(2.0_wp)
  REAL(wp), PARAMETER :: ln_3 = LOG(3.0_wp)
  REAL(wp), PARAMETER :: ln_10 = LOG(10.0_wp)
  REAL(wp), PARAMETER :: ln_006 = LOG(0.06_wp)
  REAL(wp), PARAMETER :: zlamspeat = LOG(0.25_wp) !< AYu lambdas for peat (Lawrence & Slater, 2008)
  REAL(wp), PARAMETER :: zlamdrypeat = 0.05_wp !< AYu lambdas for peat (Lawrence & Slater, 2008)

  ! AYu 4 below for the new Kersten number parametrization
  ! for the peat (Russian construction database), ice and liquid water
  REAL(wp), PARAMETER :: zkeai = 0.6116_wp !< Kersten number scale factor for frozen peat.
  REAL(wp), PARAMETER :: zkebi = 1.4123_wp !< Kersten number exponent of saturation for frozen peat.

  INTEGER :: kso
  INTEGER :: i
  INTEGER :: mstyp

  !> Water volume fraction half-way between field capacity and wilting point [m**3/m**3].
  REAL(wp) :: zwqg
  REAL(wp) :: z4wdpv
  REAL(wp) :: zdlam
  REAL(wp) :: zalamtmp

  !> tuning constant for dry thermal conductivity formula, depends on EVSL choice.
  REAL(wp) :: lamdry_c1
  !> tuning constant for dry thermal conductivity formula, depends on EVSL choice.
  REAL(wp) :: lamdry_c2

  REAL(wp) :: fr_pore !< Pore volume fraction [m**3/m**3].
  REAL(wp) :: fr_open !< Pore volume not occupied by ice [m**3/m**3].
  REAL(wp) :: fr_quartz !< Quartz-sand fraction [1].
  REAL(wp) :: lam_minerals_ln !< LOG of conductivity of non-quartz fractions of the soil.
  REAL(wp) :: lam_solid_ln !< LOG of conductivity of solid soil (no ice, no liquid water).
  REAL(wp) :: lam_sat !< Conductivity of saturated soil [W/(K m)].
  REAL(wp) :: rho_dry !< Density of dry soil [kg/m**3].
  REAL(wp) :: lam_dry_soil !< Dry-soil conductivity [W/(K m)].
  REAL(wp) :: lam_dry_total !< Dry-soil conductivity including roots [W/(K m)].
  REAL(wp) :: fr_rootvol !< Root volume fraction [m**3/m**3].
  REAL(wp) :: s_r !< Saturation w.r.t. pore volume.
  REAL(wp) :: k_e !< Kersten number.
  REAL(wp) :: alpha

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)

  SELECT CASE (itype_heatcond)
  CASE (HCOND_AVG)
    ! heat conductivity based on assumption of a soil water content which is equal to the
    ! average between wilting point and field capacity

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp, zwqg, z4wdpv, zdlam, zalamtmp)
    DO i = ivstart, ivend
      mstyp = soiltyp_subs(i)
      zwqg = 0.5_wp*(cfcap(mstyp) + cpwp(mstyp))
      z4wdpv = 4._wp*zwqg/cporv(mstyp)
      ! heat conductivity
      zdlam = cala1(mstyp) - cala0(mstyp)
      zalamtmp = zdlam * (0.25_wp + 0.30_wp*zdlam / (1._wp+0.75_wp*zdlam)) &
          & * MIN( &
          &     z4wdpv, &
          &     1.0_wp + (z4wdpv-1.0_wp) * (1.0_wp+0.35_wp*zdlam) / (1.0_wp+1.95_wp*zdlam))
      hcond_ml(i,1) = cala0(mstyp) + zalamtmp
    ENDDO

    ! Conductivity is constant throughout the column. Copy values from first layer.
    ! As a consequence, the conductivity between layer centers is the same as in each layer.
    !NEC$ outerloop_unroll(7)
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        hcond_ml(i,kso) = hcond_ml(i,1)
        hcond_hl(i,kso) = hcond_ml(i,1)
      END DO
    END DO
    !$ACC END PARALLEL

  CASE (HCOND_PL98, HCOND_PL98_VEG)
    ! heat conductivity dependent on actual soil water content
    ! based on Peters-Lidard et al. (1998) and Johansen (1975),
    ! see also Block, Alexander (2007), Dissertation BTU Cottbus

    IF (itype_evsl == EVSL_RESIST .OR. itype_evsl == EVSL_RESIST_RBS) THEN
      lamdry_c1 = 437.0_wp
      lamdry_c2 = 0.901_wp
    ELSE
      lamdry_c1 = 64.7_wp
      lamdry_c2 = 0.947_wp
    ENDIF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !NEC$ outerloop_unroll(8)
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp, fr_pore, fr_open, fr_quartz) &
      !$ACC   PRIVATE(lam_minerals_ln, lam_solid_ln, lam_sat, rho_dry, lam_dry_soil, fr_rootvol) &
      !$ACC   PRIVATE(lam_dry_total, s_r, k_e, alpha)
      DO i = ivstart, ivend
        mstyp = soiltyp_subs(i)
        fr_pore = cporv(mstyp)
        fr_open = fr_pore - w_so_ice_now(i,kso)/dz_hl(kso)

        ! quartz content
        fr_quartz = csandf(mstyp)/100._wp

        ! LOG(thermal conductivity non-quartz)
        IF (fr_quartz >= 0.2_wp)  THEN
          lam_minerals_ln = ln_2
        ELSE
          lam_minerals_ln = ln_3
        END IF

        ! saturated thermal conductivity
        IF (mstyp == IST_PEAT .AND. itype_mire == 1) THEN
          lam_solid_ln = zlamspeat ! AYu for peat
        ELSE
          lam_solid_ln = lam_quartz_ln*fr_quartz + lam_minerals_ln*(1._wp-fr_quartz)
        ENDIF

        lam_sat = EXP(lam_solid_ln*(1.0_wp-fr_pore) + lam_ice_ln*(fr_pore-fr_open) + fr_open*lam_water_ln)

        ! dry thermal conductivity

        rho_dry = rho_solids*(1.0_wp-fr_pore)
        lam_dry_soil = (0.135_wp*rho_dry + lamdry_c1) / (rho_solids - lamdry_c2*rho_dry)
        ! missing: crushed rock formulation for dry thermal conductivity (see PL98)
        ! Scale lam_dry_total with organic fraction
        IF(z_ml(kso) < root_depth(i)) THEN
          fr_rootvol = plcov(i)*(root_depth(i)-z_ml(kso))/root_depth(i)
          ! Chadburn et al.,2015, Dankers et al., 2011
          lam_dry_total = EXP(LOG(lam_dry_soil)*(1._wp-fr_rootvol)+ln_006*fr_rootvol)
        ELSE
          fr_rootvol = 0._wp
          lam_dry_total = lam_dry_soil
        END IF

        IF (mstyp == IST_PEAT .AND. itype_mire == 1) lam_dry_total = zlamdrypeat !AYu for peat

        ! Kersten number

        s_r = MIN(1.0_wp, w_so_now(i,kso)/dz_hl(kso) / fr_pore)

        IF (t_so_now(i,kso) < t0_melt) THEN                         ! frozen
          k_e = s_r
          IF (mstyp == IST_PEAT .AND. itype_mire == 1) THEN
            k_e = zkeai * s_r**zkebi !AYu for peat
          END IF
        ELSE                                                         ! unfrozen
          k_e = 0.0_wp
          SELECT CASE (mstyp)
          CASE (IST_SAND, IST_SLOAM)
            ! coarse soil
            IF (s_r >= 0.05_wp) THEN
              k_e = 0.7_wp * LOG(s_r) / ln_10 + 1.0_wp
            ENDIF
          CASE DEFAULT
            ! fine soil (other)
            IF (s_r >= 0.1_wp) THEN
              k_e = LOG(s_r) / ln_10 + 1.0_wp
            ENDIF
          END SELECT
        ENDIF
        k_e = MAX(0.0_wp, MIN(1.0_wp, k_e))

        ! thermal conductivity

        ! tuning factor to indirectly account for the impact of vegetation, which does not depend on soil moisture
        IF(itype_heatcond == HCOND_PL98_VEG .AND. z_ml(kso) < 0.075_wp) THEN
          alpha = 12.5_wp*(0.075_wp-z_ml(kso))*fr_rootvol
        ELSE
          alpha = 0.0_wp
        ENDIF
        hcond_ml(i,kso) = heatcond_fac(i)*(k_e*(lam_sat - lam_dry_total) + lam_dry_total)*(1._wp-alpha) + alpha*0.06_wp

        ! heat conductivity is also artificially reduced on snow-free forest-covered tiles generated
        ! by the melting-rate parameterization
        IF (t_snred(i) < -1.0_wp .AND. z0(i) >= 0.4_wp) THEN
          alpha = MAX(0.0_wp,2.0_wp - ABS(t_snred(i)))
          hcond_ml(i,kso) = alpha*hcond_ml(i,kso) + (1.0_wp-alpha)*0.06_wp
        ENDIF
      ENDDO
    ENDDO

    !NEC$ nofuse
    !NEC$ outerloop_unroll(7)
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        ! mean heat conductivity
        hcond_hl(i,kso) = hcond_ml(i,kso)*hcond_ml(i,kso+1)*(z_ml(kso+1)-z_ml(kso)) &
                        / ( hcond_ml(i,kso)*(z_ml(kso+1)-z_hl(kso)) &
                        +   hcond_ml(i,kso+1)*(z_hl(kso)-z_ml(kso)) )
      ENDDO
    ENDDO
    !$ACC END PARALLEL

  CASE DEFAULT
    CALL finish('calc_heat_conductivity', 'unknown heat-conductivity scheme (itype_heatcond)')

  END SELECT

  IF (lterra_urb) THEN
    ! HW: Modification of the surface thermal conductivity:
    ! - according to the building materials,
    ! - area index of buildings (urb_ai), height of building elements (urb_h_bld)
    !   and building fraction (urb_isa)
    ! - area index of natural surfaces (c_lnd) and height of natural soil elements (c_lnd_h)
    !
    ! Because of the curvature of the surface, the uppermost soil layer heat
    ! transfer is larger compared to the thermal conductivity of a plan area.
    ! As a result, the effective thermal conductivity of the upper surface is increased.
    ! This is also the case in the layers beneath, in which the effective thermal conductivity
    ! decreases with depth.
    !
    ! This modification decreases with depth with respect to the
    ! natural soil below the buildings.

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(alpha)
      DO i = ivstart, ivend
        alpha = MAX(0.0_wp, MIN(z_ml(kso)/urb_h_bld(i), 1.0_wp))
        hcond_hl(i,kso) = urb_isa(i) * ( &
            &   (1.0_wp - alpha) * urb_hcon(i)*urb_ai(i) &
            &   + alpha * hcond_hl(i,kso)) &
            & + (1.0_wp - urb_isa(i)) * hcond_hl(i,kso)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
  END IF

  !$ACC END DATA

END SUBROUTINE calc_heat_conductivity

SUBROUTINE calc_heat_conduction ( &
      & ivstart, ivend, nvec, ke_soil, ke_snow, dt, z_ml, dz_ml, dz_hl, fr_snow, t_s_now, &
      & hcond_ml, hcond_hl, hcap_ml, forcing_soil, w_snow_now, w_snow_new, wliq_snow_now, wtot_snow_now, &
      & rho_snow_mult_now, zalas_mult, zhm_snow, zdzh_snow, zdzm_snow, zfor_snow_mult, hfl_snow_soil, &
      & t_snow_mult_now, t_snow_mult_new, dt_t_snow_mult, t_so_now, t_so_new, dt_w_snow, lzacc, &
      & acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart !< Start index for computations in the parallel program.
  INTEGER, INTENT(IN) :: ivend !< End index for computations in the parallel program.
  INTEGER, INTENT(IN) :: nvec !< Block length.
  INTEGER, INTENT(IN) :: ke_soil !< Number of active soil layers.
  INTEGER, INTENT(IN) :: ke_snow !< Number of snow layers.
  REAL(wp), INTENT(IN) :: dt !< Time step [s].

  REAL(wp), INTENT(IN) :: z_ml(ke_soil+1) !< Soil main level depth [m].
  REAL(wp), INTENT(IN) :: dz_ml(ke_soil+1) !< Distance between main levels [m].
  REAL(wp), INTENT(IN) :: dz_hl(ke_soil+1) !< Distance between half levels [m].

  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< Snow fraction [m**2(snow) / m**2(tile)].
  REAL(wp), INTENT(IN) :: t_s_now(nvec) !< Surface temperature [K].
  REAL(wp), INTENT(IN) :: hcond_ml(nvec,ke_soil+1) !< Heat conductivity in soil layers [W/(K m)]
  REAL(wp), INTENT(IN) :: hcond_hl(nvec,ke_soil) !< Heat conductivity across half levels [W/(K m)]
  REAL(wp), INTENT(IN) :: hcap_ml(nvec,ke_soil+1) !< Soil heat capacity [J/(m**3 K)].
  REAL(wp), INTENT(IN) :: forcing_soil(nvec) !< Total soil forcing [W/m**2].

  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Total snow water equivalent (current) [m].
  REAL(wp), INTENT(INOUT) :: w_snow_new(nvec) !< Total snow water equivalent (new) [m].
  REAL(wp), INTENT(IN) :: wliq_snow_now(nvec,ke_snow) !< Liquid snow water equivalent in layers [m].
  REAL(wp), INTENT(IN) :: wtot_snow_now(nvec,ke_snow) !< Total snow water equivalent in layers [m].
  REAL(wp), INTENT(IN) :: rho_snow_mult_now(nvec,ke_snow) !< Snow density in layers [kg/m**3].
  REAL(wp), INTENT(IN) :: zalas_mult(nvec,ke_snow) !< Snow heat conductivity [?].
  REAL(wp), INTENT(IN) :: zhm_snow(nvec,ke_snow) !< Depth of snow main levels [m].
  REAL(wp), INTENT(IN) :: zdzh_snow(nvec,ke_snow) !< Snow layer thickness [m].
  REAL(wp), INTENT(IN) :: zdzm_snow(nvec,ke_snow) !< Distance between main snow layers [m].
  REAL(wp), INTENT(IN) :: zfor_snow_mult(nvec) !< Total forcing on snow (multilayer model) [W/m**2].
  REAL(wp), INTENT(INOUT) :: hfl_snow_soil(nvec) !< Heat flux through snow [W/m**2].

  REAL(wp), INTENT(IN) :: t_snow_mult_now(nvec,0:ke_snow) !< Snow temperature (current) [K].
  REAL(wp), INTENT(INOUT) :: t_snow_mult_new(nvec,0:ke_snow) !< Snow temperature (current) [K].
  REAL(wp), INTENT(INOUT) :: dt_t_snow_mult(nvec,0:ke_snow) !< Snow temperature tendency [K/s].

  REAL(wp), INTENT(IN) :: t_so_now(nvec,0:ke_soil+1) !< Soil temperature (current) [K].
  REAL(wp), INTENT(INOUT) :: t_so_new(nvec,0:ke_soil+1) !< Soil temperature (new) [K].
  REAL(wp), INTENT(INOUT) :: dt_w_snow(nvec) !< Rate of change of SWE [m/s].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp) :: t_so_free_new(nvec,0:ke_soil+1) !< Snow-free soil temperature [K].

  !> Adjusted snow fraction for newly created or deleted snow [m**2(snow) / m**2(tile)].
  REAL(wp) :: sn_frac(nvec)

  LOGICAL :: partial_snow_flag(nvec) !< Partial snow cover flag (or use of single-layer model).

  REAL(wp) :: lse_a(nvec,0:ke_soil+1) !< Subdiagonal matrix entry.
  REAL(wp) :: lse_b(nvec,0:ke_soil+1) !< Diagonal matrix entry.
  REAL(wp) :: lse_c(nvec,0:ke_soil+1) !< Superdiagonal matrix entry.
  REAL(wp) :: lse_rhs(nvec,0:ke_soil+1) !< Right-hand side.

  REAL(wp) :: zakb1, zakb2

  INTEGER :: i, kso

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA CREATE(partial_snow_flag, lse_a, lse_b, lse_c, lse_rhs, sn_frac, t_so_free_new) &
  !$ACC   PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    sn_frac(i) = fr_snow(i)
    w_snow_new(i) = w_snow_now(i) + dt_w_snow(i) * (dt / rho_w)
    IF (w_snow_now(i) < eps_soil .AND. w_snow_new(i) >= eps_soil) THEN
      sn_frac(i) = 0.01_wp
    ELSEIF (w_snow_now(i) >= eps_soil .AND. w_snow_new(i) < eps_soil) THEN
      sn_frac(i) = 0._wp
    ENDIF
    partial_snow_flag(i) = .NOT. lmulti_snow .OR. sn_frac(i) < 1._wp
  END DO

  !$ACC LOOP SEQ
  DO kso = 1, ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      t_so_free_new(i,kso) = t_so_now(i,kso)
    END DO
  END DO        ! soil layers

  !$ACC LOOP SEQ
  DO kso = 2, ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zakb1, zakb2)
    DO i = ivstart, ivend
      IF (partial_snow_flag(i)) THEN
        ! for heat conductivity: hcond_ml is now 3D
        zakb1 = hcond_hl(i,kso-1)/hcap_ml(i,kso)
        zakb2 = hcond_hl(i,kso  )/hcap_ml(i,kso)
        lse_a(i,kso) = -zalfa*dt*zakb1/(dz_hl(kso)*dz_ml(kso))
        lse_c(i,kso) = -zalfa*dt*zakb2/(dz_hl(kso)*dz_ml(kso+1))
        lse_b(i,kso) = 1._wp - lse_a(i,kso) - lse_c(i,kso)
        lse_rhs(i,kso) = t_so_now(i,kso) +                                     &
               (1._wp - zalfa)*( - lse_a(i,kso)/zalfa*t_so_now(i,kso-1)+ &
               (lse_a(i,kso)/zalfa + lse_c(i,kso)/zalfa)*t_so_now(i,kso) -  &
                lse_c(i,kso)/zalfa*t_so_now(i,kso+1)  )
      END IF
    END DO
  END DO        ! soil layers

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zakb1, zakb2)
  DO i = ivstart, ivend
    IF (partial_snow_flag(i)) THEN
      ! for heat conductivity: hcond_ml is now 3D: here we need layer 1
      zakb1 = hcond_ml(i,1)/hcap_ml(i,1)
      zakb2 = hcond_hl(i,1)/hcap_ml(i,1)
      lse_a(i,1) = -zalfa*dt*zakb1/(dz_hl(1)*dz_ml(1))
      lse_c(i,1) = -zalfa*dt*zakb2/(dz_hl(1)*dz_ml(2))
      lse_b(i,1) = 1._wp - lse_a(i,1) - lse_c(i,1)
      lse_rhs(i,1) = t_so_now(i,1) + (1._wp - zalfa)* (                 &
                      - lse_a(i,1)/zalfa * t_s_now(i) +                     &
                      (lse_a(i,1) + lse_c(i,1))/zalfa * t_so_now(i,1) -    &
                       lse_c(i,1)/zalfa * t_so_now(i,2)   )
      lse_a(i,0) = 0.0_wp
      lse_b(i,0) = zalfa
      lse_c(i,0) = -zalfa
      ! EM: In the case of multi-layer snow model, forcing_soil(i) does not include the heat conductivity flux
      ! between soil and snow (hfl_snow_soil). It will be accounted for at the next semi-step, see below.
      lse_rhs(i,0)    = dz_ml(1) * forcing_soil(i)/hcond_ml(i,1)+(1._wp-zalfa)* &
                      (t_so_now(i,1) - t_s_now(i))
      lse_a(i,ke_soil+1) = 0.0_wp
      lse_b(i,ke_soil+1) = 1.0_wp
      lse_c(i,ke_soil+1) = 0.0_wp
      lse_rhs(i,ke_soil+1) = t_so_now(i,ke_soil+1)
    END IF
  END DO

  IF (.NOT. lmulti_snow) THEN
    ! Solve on every cell (partial_snow_flag is always TRUE).
    CALL solve_tridiag ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & nlev=ke_soil+2, & ! 0:ke_soil+1
        & a=lse_a(:,0:), &
        & b=lse_b(:,0:), &
        & c=lse_c(:,0:), &
        & d=lse_rhs(:,0:), &
        & out=t_so_free_new(:,0:) &
      )
    !$ACC LOOP SEQ
    DO kso = 1, ke_soil+1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
        t_so_new(i,kso) = t_so_free_new(i,kso)
      END DO
    END DO
  ELSE
    ! Use partial_snow_flag as mask.
    CALL solve_tridiag ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & nlev=ke_soil+2, & ! 0:ke_soil+1
      & a=lse_a(:,0:), &
      & b=lse_b(:,0:), &
      & c=lse_c(:,0:), &
      & d=lse_rhs(:,0:), &
      & out=t_so_free_new(:,0:), &
      & mask=partial_snow_flag(:) &
    )
  END IF
  !$ACC END PARALLEL

  IF (lmulti_snow) THEN

    ! If there is snow, the solution of the heat conduction equation
    ! goes through the whole column "soil+snow"

    CALL snow_multi_calc_heat_conduction ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & ke_snow=ke_snow, &
        & ke_soil=ke_soil, &
        & dt=dt, &
        & t_s_now=t_s_now(:), &
        & t_so_now=t_so_now(:,0:), &
        & t_so_new=t_so_new(:,0:), &
        & t_so_free_new=t_so_free_new(:,0:), &
        & hcond_ml=hcond_ml(:,:), &
        & hcond_hl=hcond_hl(:,:), &
        & hcap_ml=hcap_ml(:,:), &
        & z_ml=z_ml(:), &
        & dz_ml=dz_ml(:), &
        & dz_hl=dz_hl(:), &
        & t_snow_mult_now=t_snow_mult_now(:,0:), &
        & t_snow_mult_new=t_snow_mult_new(:,0:), &
        & dt_t_snow_mult=dt_t_snow_mult(:,0:), &
        & wliq_snow_now=wliq_snow_now(:,:), &
        & wtot_snow_now=wtot_snow_now(:,:), &
        & rho_snow_mult_now=rho_snow_mult_now(:,:), &
        & zalas_mult=zalas_mult(:,:), &
        & zhm_snow=zhm_snow(:,:), &
        & zdzh_snow=zdzh_snow(:,:), &
        & zdzm_snow=zdzm_snow(:,:), &
        & zfor_snow_mult=zfor_snow_mult(:), &
        & sn_frac=sn_frac(:), &
        & w_snow_now=w_snow_now(:), &
        & w_snow_new=w_snow_new(:), &
        & dt_w_snow=dt_w_snow(:), &
        & hfl_snow_soil=hfl_snow_soil(:), &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )

  END IF
  !$ACC END DATA

END SUBROUTINE calc_heat_conduction


SUBROUTINE calc_soil_water_melt ( &
      & ivstart, ivend, nvec, ke_soil, dt, soiltyp_subs, z_ml, dz_hl, root_depth, plcov, &
      & dt_w_so, hcap_ml, w_so_now, w_so_ice_now, w_so_ice_new, t_so_now, t_so_new, lzacc, &
      & acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart !< Start index for computations in the parallel program.
  INTEGER, INTENT(IN) :: ivend !< End index for computations in the parallel program.
  INTEGER, INTENT(IN) :: nvec !< Block length.
  INTEGER, INTENT(IN) :: ke_soil !< Number of active soil layers.
  REAL(wp), INTENT(IN) :: dt !< Time step [s]
  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type.
  REAL(wp), INTENT(IN) :: z_ml(ke_soil+1) !< Layer depth [m].
  REAL(wp), INTENT(IN) :: dz_hl(ke_soil+1) !< Layer thickness [m].
  REAL(wp), INTENT(IN) :: root_depth(nvec) !< Depth of root zone [m].
  REAL(wp), INTENT(IN) :: plcov(nvec) !< Plant-covered area fraction [m**2/m**2].
  REAL(wp), INTENT(IN) :: dt_w_so(nvec, ke_soil) !< Tendency of water content [m/s].
  REAL(wp), INTENT(IN) :: hcap_ml(nvec, ke_soil) !< Volumetric heat capacity [J/(m**3 K)].
  REAL(wp), INTENT(IN) :: w_so_now(nvec, ke_soil) !< Total layer water (liq+ice) content [m].
  REAL(wp), INTENT(IN) :: w_so_ice_now(nvec, ke_soil) !< Current layer ice content [m].
  REAL(wp), INTENT(INOUT) :: w_so_ice_new(nvec, ke_soil) !< New layer ice content [m].
  REAL(wp), INTENT(IN) :: t_so_now(nvec, 0:ke_soil) !< Current soil temperature [K].
  REAL(wp), INTENT(INOUT) :: t_so_new(nvec, 0:ke_soil) !< New soil temperature [K].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  !> Upper limit for logarithmic interpolation [K].
  REAL(wp), PARAMETER :: T_ZW_UP = t0_melt - 3._wp
  !> Lower limit for logarithmic interpolation [K].
  REAL(wp), PARAMETER :: T_ZW_LO = t0_melt - 40._wp
  !> Melting timescale [s].
  REAL(wp), PARAMETER :: TAU_MELT = 1800._wp

  REAL(wp), PARAMETER :: zpsi0 = 0.01_wp !< air entry potential at water saturation [m]

  ! Soil ice parameterization according to K. Schaefer and Jafarov, E.,2016,
  ! doi:10.5194/bg-13-1991-2016, Exponents and temperature constants can be found in sfc_terra_data.

  REAL(wp), PARAMETER :: zd = LOG((T_ref_ice-(T_ZW_LO-t0_melt))/T_star_ice)
  REAL(wp), PARAMETER :: zd1 = EXP(b_sand*zd)
  REAL(wp), PARAMETER :: zd2 = EXP(b_clay*zd)
  REAL(wp), PARAMETER :: zd3 = EXP(b_silt*zd)
  REAL(wp), PARAMETER :: zd4 = EXP(b_org*zd)

  REAL(wp) :: zpsis
  REAL(wp) :: fr_rootvol !< Volume fraction inside root zone [m**3/m**3].

  REAL(wp) :: b_por(nvec)
  REAL(wp) :: zaa(nvec)
  REAL(wp) :: w_m_up(nvec) !< Minimum liquid water in layer at T_ZW_UP [m].
  REAL(wp) :: w_m_low !< Minimum liquid water in layer at T_ZW_LO [m].
  !> Minimum liquid fraction of the non-organic part at T_ZW_DOWN [m**3(H2O)/m**3(pore)].
  REAL(wp) :: fr_frozen_lo(nvec)

  REAL(wp) :: zw_m !< Minimum liquid water amount [m].
  REAL(wp) :: w_liq !< Current liquid water amount [m].
  REAL(wp) :: znen !< ??? [1]
  REAL(wp) :: ztx !< ??? [K].
  REAL(wp) :: icemelt_in_m_per_k !< Melted ice per temperature change [m/K].
  REAL(wp) :: delta_w_ice !< Ice change during timestep [m].
  REAL(wp) :: w_so_new !< New total water content [m].
  REAL(wp) :: w_avail !< Liquid water available for freezing [m].
  REAL(wp) :: melting_time !< Timestep over melting timescale [1].

  INTEGER :: mstyp
  INTEGER :: i
  INTEGER :: kso

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) CREATE(b_por, zaa, w_m_up, fr_frozen_lo) ASYNC(acc_async_queue)

  melting_time = dt / TAU_MELT

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp, zpsis)
  DO i = ivstart, ivend
    mstyp = soiltyp_subs(i)        ! soil type
    zpsis = -zpsi0 * EXP(LOG(10._wp)*(1.88_wp - 0.013_wp*csandf(mstyp)))
    b_por(i) = 2.91_wp + 0.159_wp*cclayf(mstyp)
    zaa(i)   = g*zpsis/lh_f

    ! Liq. water content at -3 degC (without cporv * dz_hl)
    w_m_up(i) = EXP(-1.0_wp/b_por(i)*LOG((T_ZW_UP - t0_melt)/(T_ZW_UP*zaa(i))) )

    ! Determine liq. water content at -40 degC
    fr_frozen_lo(i) = 0.01_wp * ( &
        & csandf(mstyp)*zd1 + cclayf(mstyp)*zd2 + (100.0_wp-csandf(mstyp)-cclayf(mstyp))*zd3)
  ENDDO

  !$ACC LOOP SEQ
  DO kso = 1,ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ztx, zw_m, fr_rootvol, w_m_low, w_liq, znen) &
    !$ACC   PRIVATE(icemelt_in_m_per_k, delta_w_ice, w_so_new, w_avail)
    DO i = ivstart, ivend
      mstyp = soiltyp_subs(i)

      IF (mstyp > IST_ROCK) THEN
        ztx = t0_melt
        zw_m = cporv(mstyp) * dz_hl(kso)

        IF(t_so_new(i,kso) < t0_melt-eps_temp) THEN
          IF (z_ml(kso) < root_depth(i)) THEN
            ! Add in organic fraction.
            fr_rootvol = plcov(i)*(root_depth(i)-z_ml(kso))/root_depth(i)
            w_m_low = cporv(mstyp) * dz_hl(kso) * &
                & (fr_rootvol*zd4 + (1.0_wp-fr_rootvol)*fr_frozen_lo(i))
          ELSE
            w_m_low = cporv(mstyp) * dz_hl(kso) * fr_frozen_lo(i)
          END IF

          IF (t_so_new(i,kso) < T_ZW_LO) THEN
            zw_m = w_m_low
          ELSE IF (t_so_new(i,kso) < T_ZW_UP) THEN ! Logarithmic Interpolation between -3 degC and -40 degC
            zw_m = w_m_low * EXP((t_so_new(i,kso) - T_ZW_LO) * &
              (LOG(cporv(mstyp) * dz_hl(kso) * w_m_up(i)) - LOG(w_m_low))/(T_ZW_UP-T_ZW_LO))
          ELSE
            zw_m = zw_m * EXP(-1._wp/b_por(i) * &
              LOG((t_so_new(i,kso) - t0_melt)/(t_so_new(i,kso)*zaa(i))))
          END IF

          w_liq = MAX(eps_div, w_so_now(i,kso) -  w_so_ice_now(i,kso))
          znen = 1._wp - zaa(i) * EXP(b_por(i) * LOG(cporv(mstyp)*dz_hl(kso)/w_liq))
          ztx = t0_melt/znen
        ENDIF

        ztx = MIN(t0_melt, ztx)
        icemelt_in_m_per_k = hcap_ml(i,kso) * dz_hl(kso) / (lh_f * rho_w)
        delta_w_ice = -icemelt_in_m_per_k * (t_so_new(i,kso) - ztx)
        w_so_new = w_so_now(i,kso) + dt*dt_w_so(i,kso)/rho_w
        w_avail = w_so_new - zw_m - w_so_ice_now(i,kso)

        IF (t_so_new(i,kso) > t0_melt .AND. w_so_ice_now(i,kso) > 0.0_wp) THEN
          ! melting point adjustment (time scale 30 min)
          delta_w_ice = - MIN( &
              & w_so_ice_now(i,kso), &
              & melting_time * (t_so_new(i,kso)-t0_melt) * icemelt_in_m_per_k &
            )
        ELSE IF (delta_w_ice < 0.0_wp) THEN ! this branch contains cases of melting and freezing
          delta_w_ice = - MIN(-delta_w_ice, -w_avail, w_so_ice_now(i,kso))
          ! limit latent heat consumption due to melting to half the temperature increase since
          ! last time step or 2.5 K within 30 min; the freezing rate is limited below.
          delta_w_ice = - MIN( &
              & -delta_w_ice, &
              & MAX(2.5_wp * melting_time, 0.5_wp*(t_so_new(i,kso)-t_so_now(i,kso))) &
              &   * icemelt_in_m_per_k &
            )
        ELSE
          delta_w_ice = MIN(delta_w_ice, MAX(w_avail, 0.0_wp))
        ENDIF

        IF (delta_w_ice > 0.0_wp) THEN
          ! limit latent heat release due to freezing to half the difference from the melting point.
          delta_w_ice = MIN( &
              & delta_w_ice, &
              & 0.5_wp * MAX(0.0_wp, (t0_melt - t_so_new(i,kso))) * icemelt_in_m_per_k &
            )
        ENDIF
        w_so_ice_new(i,kso) = w_so_ice_now(i,kso) + delta_w_ice
        t_so_new    (i,kso) = t_so_new    (i,kso) + delta_w_ice / icemelt_in_m_per_k
      END IF
    END DO
  ENDDO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE calc_soil_water_melt


REAL(wp) FUNCTION watrcon_RT(ks,lw,kw1,pv,adp)

  !$ACC ROUTINE SEQ
  REAL(wp), INTENT(IN) :: ks,lw,kw1,pv,adp
  watrcon_RT = ks*EXP(kw1*MAX(0.0_wp,pv-lw)/(pv-adp))

END FUNCTION watrcon_RT

!==============================================================================

REAL(wp) FUNCTION watrdiff_RT(ds,lw,dw1,pv,adp)

  !$ACC ROUTINE SEQ
  REAL(wp), INTENT(IN) :: ds,lw,dw1,pv,adp
  watrdiff_RT = ds*EXP(dw1*MAX(0.0_wp,pv-lw)/(pv-adp))

END FUNCTION watrdiff_RT

END MODULE sfc_terra_transport
