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

!> TERRA evapotranspiration routines.
MODULE sfc_terra_evap

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish

  USE mo_math_constants, ONLY: pi
  USE mo_physical_constants, ONLY: &
      & b3 => t3, lh_v  => alv, o_m_rdv, rdv, rho_w => rhoh2o, rvd_m_o => vtmpc1, t0_melt => tmelt
  USE mo_lookup_tables_constants, ONLY: &
      & b1 => c1es, b2w => c3les, b2i => c3ies , b4w => c4les, b4i => c4ies

  USE mo_lnd_nwp_config, ONLY: &
      & cwimax_ml, itype_canopy, itype_evsl, itype_eisa, itype_trvg, itype_root, lstomata, &
      & lterra_urb

  USE sfc_terra_data, ONLY: &
      & cadp, cbedi, cdash, cdmin, cfcap, cfinull, cf_w, ck0di, clgk0, cparcrit, cporv, cpwp, &
      & crhowm, crsmax, crsmin, ctend, cwisamax, eps_div, eps_soil, itype_mire, &
      & EVSL_BATS, EVSL_NP89, EVSL_RESIST, EVSL_RESIST_RBS, ROOT_EXPONENTIAL, &
      & IST_ICE, IST_NUM, IST_PEAT, IST_ROCK, TRVG_BATS, TRVG_BATS_EXT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_evapotranspiration

  REAL(wp), PARAMETER :: BATS_ROOT_DENSITY = 3._wp !< Root density scaling parameter.


! Silence unused variable warnings on non-OpenACC builds
#ifdef _OPENACC
# define OPENACC_SUPPRESS_UNUSED_LZACC
#else
# define OPENACC_SUPPRESS_UNUSED_LZACC IF (lzacc .AND. acc_async_queue > 0) THEN; END IF
#endif

CONTAINS

!>
!! Compute evapotranspiration from soil, snow, interception, and plants.
!!
!!
SUBROUTINE calc_evapotranspiration ( &
      & dt, icant, ivstart, ivend, nvec, n_soil, n_soil_hy, u_atm, v_atm, qv_atm, t_atm, rho_atm, &
      & p_s, t_s, t_sk, t_snow_top, fr_snow, w_snow, t_snred, w_i, w_so, w_so_ice, fr_w_ml, tai, &
      & eai, sai, plcov, laifac, rad_flx, par_absorbed, z_ml, z_hl, dz_hl, soiltyp_subs, root_depth, &
      & r_bsmin, r_stommin, tcm, tch, tfv, tfvsn, z0, rho_ch, urb_isa, &
      & plevap, &
      & qv_s, eva_bs, lhfl_bs, transp_ml, transp_sum, lhfl_pl, r_stom, dqvdt_snow, eva_w_i, eva_w_sn, &
      & dew_rate, rime_rate, lzacc, acc_async_queue &
    )

  REAL(wp), INTENT(IN) :: dt !< Time step [s].
  INTEGER, INTENT(IN)  :: icant !< Canopy transfer scheme (1=Louis, 2=Raschendorfer, 3=EDMF).
  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec
  INTEGER, INTENT(IN)  :: n_soil !< Number of soil moisture layers.
  INTEGER, INTENT(IN)  :: n_soil_hy !< Number of active soil moisture layers.
  REAL(wp), INTENT(IN) :: u_atm(nvec) !< Zonal wind in lowest layer (n_cell) [m/s].
  REAL(wp), INTENT(IN) :: v_atm(nvec) !< Meridional wind in lowest layer (n_cell) [m/s].
  REAL(wp), INTENT(IN) :: qv_atm(nvec) !< Specific humidity in lowest level (n_cell) [kg(vap)/kg(air)].
  REAL(wp), INTENT(IN) :: t_atm(nvec) !< Temperature in lowest level (n_cell) [K].
  REAL(wp), INTENT(IN) :: rho_atm(nvec) !< Surface air density (n_cell) [kg(air)/m**3]
  REAL(wp), INTENT(IN) :: p_s(nvec) !< Surface pressure (n_cell) [Pa].
  REAL(wp), INTENT(IN) :: t_s(nvec) !< Surface temperature (n_cell) [K].
  REAL(wp), INTENT(IN) :: t_sk(nvec) !< Skin temperature (n_cell) [K].
  REAL(wp), INTENT(IN) :: t_snow_top(nvec) !< Snow surface temperature (n_cell) [K].
  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< Snow-covered fraction (n_cell) [m**2(snow)/m**2(tile)].
  REAL(wp), INTENT(IN) :: w_snow(nvec) !< Amount of snow on tile (n_cell) [m(H2O)].
  REAL(wp), INTENT(IN) :: t_snred(nvec) !< Snow temperature offset for calculating evaporation (n_cell) [K].
  REAL(wp), INTENT(IN) :: w_i(nvec) !< Water level in interception reservoir (n_cell) [m(H2O)].
  REAL(wp), INTENT(IN) :: w_so(nvec,n_soil+1) !< Water in soil layers (n_cell, n_soil) [m(H2O)].
  REAL(wp), INTENT(IN) :: w_so_ice(nvec,n_soil+1) !< Soil ice content (n_cell, n_soil) [m(H2O)].
  REAL(wp), INTENT(IN) :: fr_w_ml(nvec,n_soil+1) !< Water in soil layers relative to layer thickness (n_cell, n_soil) [m(H2O)/m(dz)].
  REAL(wp), INTENT(IN) :: tai(nvec) !< Transpiration area index (n_cell) [m**2(transp)/m**2(tile)].
  REAL(wp), INTENT(IN) :: eai(nvec) !< Evaporation area index (n_cell) [m**2(evap)/m**2(tile)].
  REAL(wp), INTENT(IN) :: sai(nvec) !< Surface area index (n_cell) [m**2(surf)/m**2(tile)].
  REAL(wp), INTENT(IN) :: plcov(nvec) !< Plant-covered fraction (n_cell) [m**2(plants)/m**2(tile)].
  REAL(wp), INTENT(IN) :: laifac(nvec) !< Ratio of current LAI to laimax (n_cell) [1]
  REAL(wp), INTENT(IN) :: rad_flx(nvec) !< Net radiation flux (n_cell) [W/m**2(tile)].
  REAL(wp), INTENT(IN) :: par_absorbed(nvec) !< Amount of absorbed PAR radiation (n_cell) [W/m**2].
  REAL(wp), INTENT(IN) :: z_ml(n_soil+1) !< Depth of main layers (n_soil) [m].
  REAL(wp), INTENT(IN) :: z_hl(n_soil+1) !< Depth of layer interfaces (n_soil) [m].
  REAL(wp), INTENT(IN) :: dz_hl(n_soil+1) !< Layer thickness (n_soil) [m].
  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type (n_cell).
  REAL(wp), INTENT(IN) :: root_depth(nvec) !< Root depth (n_cell) [m].
  REAL(wp), INTENT(IN) :: r_bsmin(nvec) !< Minimum bare-soil evaporation resistance (n_cell) [s/m].
  REAL(wp), INTENT(IN) :: r_stommin(nvec) !< Minimum stomata resistance (n_cell) [s/m].
  REAL(wp), INTENT(IN) :: tcm(nvec) !< Turbulent transfer coefficient for momentum (n_cell) [1].
  REAL(wp), INTENT(IN) :: tch(nvec) !< Turbulent transfer coefficient for heat (n_cell) [1].
  REAL(wp), INTENT(IN) :: tfv(nvec) !< laminar reduction factor for evaporation [1].
  REAL(wp), INTENT(IN) :: tfvsn(nvec) !< reduction factor for snow evaporation from model-DA coupling [1].
  REAL(wp), INTENT(IN) :: z0(nvec) !< Roughness length (n_cell) [m].
  REAL(wp), INTENT(IN) :: rho_ch(nvec) !< Surface air density times transfer velocity (n_cell) [kg(air)/(m**2 s)].
  REAL(wp), INTENT(IN) :: urb_isa(nvec) !< Urban impervious surface area fraction [m**2(urb)/m**2(tile)].

  REAL(wp), INTENT(INOUT) :: plevap(nvec) !< Accumulated plant evaporation since start of day (n_cell) [kg/m**2].

  REAL(wp), INTENT(OUT) :: qv_s(nvec) !< Ficticious surface humidity (n_cell) [kg(vap)/kg(air)].
  REAL(wp), INTENT(OUT) :: eva_bs(nvec) !< Bare-soil evaporation (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: lhfl_bs(nvec) !< Bare-soil latent heat flux (n_cell) [W/m**2].
  REAL(wp), INTENT(OUT) :: transp_ml(nvec,n_soil) !< Transpiration from soil layers (n_cell,n_soil) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: transp_sum(nvec) !< Total transpiration (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: lhfl_pl(nvec,n_soil) !< Latent heat flux through plants from soil layer (n_cell,n_soil) [W/m**2].
  REAL(wp), INTENT(OUT) :: r_stom(nvec) !< Stomata resistance (n_cell) [s/m].
  REAL(wp), INTENT(OUT) :: dqvdt_snow(nvec) !< dqsat/dT at snow temperature [kg(vap)/kg(air)/K].
  REAL(wp), INTENT(OUT) :: eva_w_i(nvec) !< Evaporation from interception water [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: eva_w_sn(nvec) !< Evaporation from snow [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: dew_rate(nvec) !< Rate of dew formation [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: rime_rate(nvec) !< Rate of rime formation [kg/(m**2 s)].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  REAL(wp) :: evapot_s(nvec) !< Potential evaporation for water [kg/(m**2 s)].
  REAL(wp) :: evapot_snow(nvec) !< Potential evaporation for snow [kg/(m**2 s)].
  REAL(wp) :: fr_w_i(nvec) !< Fraction of surface covered by interception water [m**2/m**2(tile)].

  REAL(wp) :: eva_sum !< Sum of evapotranspiration contributions [kg/(m**2 s)].
  REAL(wp) :: w_i_scale
  REAL(wp) :: b2iw, b4iw, b234iw, q_s, dq_s, q_snow, dq_snow
  REAL(wp) :: smth_heav, area_fac

  INTEGER :: i

  !----------------------------------------------------------------------------
  ! Section I.4.1: Evaporation from interception store and from snow cover,
  !----------------------------------------------------------------------------
  ! Evaporation and transpiration are negative, dew and rime
  ! positive quantities, since positive sign indicates a flux
  ! directed towards the earth's surface!

  !$ACC DATA PRESENT(ivend) CREATE(evapot_s, evapot_snow, fr_w_i) ASYNC(acc_async_queue)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(w_i_scale, b2iw, b4iw, b234iw, q_s, dq_s, q_snow) &
  !$ACC   PRIVATE(dq_snow, smth_heav, area_fac)
  DO i = ivstart, ivend
    ! Compute fraction covered by interception water.
    w_i_scale  = MAX( 0.25_wp * cf_w, 0.4_wp * cwimax_ml * MAX(2.5_wp * plcov(i), tai(i)))

    IF (lterra_urb .AND. itype_eisa == 3) THEN
      w_i_scale = urb_isa(i) * cwisamax + (1.0_wp - urb_isa(i)) * w_i_scale
    END IF

    fr_w_i(i) = MERGE( &
        & MAX(0.01_wp, 1.0_wp - EXP(MAX( -5.0_wp, -w_i(i)/w_i_scale))), &
        & 0._wp, &
        & w_i(i) >= 1.0E-4_wp*eps_soil &
      )

    ! Compute potential evaporation.
    b2iw = MERGE(b2w, b2i, t_sk(i) >= t0_melt)
    b4iw = MERGE(b4w, b4i, t_sk(i) >= t0_melt)
    q_s = zsf_qsat(zsf_psat_iw(t_sk(i), b2iw, b4iw), p_s(i))
    dq_s = qv_atm(i) - q_s
    IF (ABS(dq_s) < 0.01_wp * eps_soil) dq_s = 0.0_wp

    b2iw = MERGE(b2w, b2i, t_snow_top(i) >= t0_melt)
    b4iw = MERGE(b4w, b4i, t_snow_top(i) >= t0_melt)
    b234iw = b2iw*(b3 - b4iw)
    q_snow = zsf_qsat(zsf_psat_iw(t_snow_top(i) - MAX(0.0_wp,t_snred(i)), b2iw, b4iw), p_s(i))
    dq_snow = qv_atm(i) - q_snow
    IF (ABS(dq_snow) < 0.01_wp*eps_soil) dq_snow = 0.0_wp

    dqvdt_snow(i) = zsf_dqvdt_iw(t_snow_top(i), q_snow, b4iw, b234iw)

    evapot_s(i) = tfv(i) * rho_ch(i) * dq_s
    evapot_snow(i) = MERGE(tfv(i) * rho_ch(i) * dq_snow, 0._wp, t_snow_top(i) < t0_melt) &
        & * MERGE(tfvsn(i), 1._wp, dq_snow<0._wp)

    ! Evaporation from interception store if it contains water (w_i>0) and
    ! if evapot_s<0 indicates potential evaporation for temperature Ts
    ! amount of water evaporated is limited to total content of store

    ! Only 1/3 to 1/2 of the wet area actually participates, depending on skin temperature. Linear
    ! transition between 0C and 2C.
    smth_heav = MAX(0._wp, 1._wp - MAX(0._wp, 0.5_wp*(t_sk(i) - t0_melt)))
    area_fac = (1._wp + 0.5_wp*smth_heav)/3._wp

    eva_w_i(i) = MERGE(MAX( &
        & & ! Evaporate freely, ...
        &   area_fac * (1.0_wp - fr_snow(i)) * fr_w_i(i) * evapot_s(i), &
        & & ! ... but no more than the available water, ...
        &   -rho_w * w_i(i) / dt, &
        & & ! ... and no more than what 75% of net radiation or 300 W/m**2 can support.
        &   -MAX(300.0_wp,0.75_wp*rad_flx(i))/lh_v), &
        & 0._wp, &
        & evapot_s(i) < 0._wp &
      )

    ! Evaporation of snow, if snow exists (w_snow>0) and if evapot_snow<0
    ! indicates potential evaporation for temperature t_snow
    eva_w_sn(i) = MERGE( &
        & MAX(-rho_w * w_snow(i) / dt, fr_snow(i) * evapot_snow(i)), &
        & 0._wp, &
        & evapot_snow(i) < 0._wp &
      )

    ! Formation of dew or rime, if evapot_s > 0. Distinction between
    ! dew or rime is only controlled by sign of surface temperature
    ! and not affected by presence of snow !
    dew_rate(i) = MERGE(evapot_s(i), 0._wp, t_sk(i) >= t0_melt .AND. evapot_s(i) >= 0.)
    rime_rate(i) = MERGE(evapot_snow(i), 0._wp, t_sk(i) < t0_melt .AND. evapot_snow(i) >= 0.)
  END DO
  !$ACC END PARALLEL


  !----------------------------------------------------------------------------
  ! Section I.4.2: Bare soil evaporation
  !----------------------------------------------------------------------------
  SELECT CASE (itype_evsl)
  CASE (EVSL_BATS)
    CALL calc_evsl_bats ( &
        & dt=dt, &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & n_soil=n_soil, &
        & z_ml=z_ml, &
        & z_hl=z_hl, &
        & dz_hl=dz_hl, &
        & soiltyp_subs=soiltyp_subs, &
        & fr_w_i=fr_w_i, &
        & fr_snow=fr_snow, &
        & eai=eai, &
        & sai=sai, &
        & evapot_s=evapot_s, &
        & w_so=w_so(:,:), &
        & fr_w_top=fr_w_ml(:,1), &
        & & ! out
        & eva_bs=eva_bs, &
        & lhfl_bs=lhfl_bs, &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )

  CASE (EVSL_NP89)
    CALL calc_evsl_np89 ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & soiltyp_subs=soiltyp_subs, &
        & fr_snow=fr_snow, &
        & eai=eai, &
        & sai=sai, &
        & rho_ch=rho_ch, &
        & qv_atm=qv_atm, &
        & evapot_s=evapot_s, &
        & t_snow_top=t_snow_top, &
        & t_sk=t_sk, &
        & t_s=t_s, &
        & p_s=p_s, &
        & fr_w_top=fr_w_ml(:,1), &
        & & ! out
        & eva_bs=eva_bs, &
        & lhfl_bs=lhfl_bs, &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )

  CASE (EVSL_RESIST, EVSL_RESIST_RBS)
    CALL calc_evsl_resistance ( &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & soiltyp_subs=soiltyp_subs, &
        & fr_snow=fr_snow, &
        & eai=eai, &
        & sai=sai, &
        & rho_ch=rho_ch, &
        & rho_atm=rho_atm, &
        & r_bsmin=r_bsmin, &
        & evapot_s=evapot_s, &
        & fr_w_top=fr_w_ml(:,1), &
        & & ! out
        & eva_bs=eva_bs, &
        & lhfl_bs=lhfl_bs, &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )

  CASE DEFAULT
    CALL finish('terra:calc_evapotranspiration', 'Unknown bare-soil evaporation scheme')
  END SELECT

  !----------------------------------------------------------------------------
  ! Section I.4.3b: transpiration by plants, BATS version
  !----------------------------------------------------------------------------

  SELECT CASE (itype_trvg)
  CASE (TRVG_BATS, TRVG_BATS_EXT)
    !NEC$ inline_complete
    CALL calc_trvg_bats ( &
        & dt=dt, &
        & icant=icant, &
        & ivstart=ivstart, &
        & ivend=ivend, &
        & nvec=nvec, &
        & n_soil=n_soil, &
        & n_soil_hy=n_soil_hy, &
        & soiltyp_subs=soiltyp_subs, &
        & z_ml=z_ml, &
        & dz_hl=dz_hl, &
        & fr_w_ml=fr_w_ml, &
        & w_so_ice=w_so_ice, &
        & u_atm=u_atm, &
        & v_atm=v_atm, &
        & t_atm=t_atm, &
        & tcm=tcm, &
        & tch=tch, &
        & z0=z0, &
        & r_stommin=r_stommin, &
        & sai=sai, &
        & tai=tai, &
        & plcov=plcov, &
        & laifac=laifac, &
        & fr_snow=fr_snow, &
        & evapot_s=evapot_s, &
        & par_absorbed=par_absorbed, &
        & root_depth=root_depth, &
        & & ! inout
        & plevap=plevap, &
        & & ! out
        & transp_ml=transp_ml, &
        & transp_sum=transp_sum, &
        & lhfl_pl=lhfl_pl, &
        & r_stom=r_stom, &
        & lzacc=lzacc, &
        & acc_async_queue=acc_async_queue &
      )

  CASE DEFAULT
    CALL finish('terra:calc_evapotranspiration', 'Unknown plant transpiration scheme')
  END SELECT

  !----------------------------------------------------------------------------
  ! Section I.4.4: total evapotranspiration and
  !              associated fictitious soil humidity qv_s
  !----------------------------------------------------------------------------

  !NEC$ inline_complete
  CALL limit_evaporation ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & n_soil=n_soil, &
      & evapot_s=evapot_s, &
      & evapot_snow=evapot_snow, &
      & t_snred=t_snred, &
      & fr_snow=fr_snow, &
      & & ! inout
      & eva_w_i=eva_w_i, &
      & eva_w_sn=eva_w_sn, &
      & eva_bs=eva_bs, &
      & lhfl_bs=lhfl_bs, &
      & transp_sum=transp_sum, &
      & transp_ml=transp_ml(:,:), &
      & lhfl_pl=lhfl_pl(:,:), &
      & lzacc=lzacc, &
      & acc_async_queue=acc_async_queue &
  )

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(eva_sum)
  DO i = ivstart, ivend

    eva_sum = eva_w_sn(i)  & ! evaporation of snow
            + eva_w_i(i)  & ! evaporation from interception store
            + eva_bs(i)   & ! evaporation from bare soil
            + transp_sum(i)   & ! transpiration from all soil layers
            + dew_rate(i) & ! formation of dew
            + rime_rate(i)  ! formation of rime

    qv_s(i) = qv_atm(i) - eva_sum / (rho_ch(i) + eps_div)
  END DO
  !$ACC END PARALLEL

  !$ACC END DATA

END SUBROUTINE calc_evapotranspiration


!>
!! Bare-soil evaporation, BATS version.
!!
!! Calculation of bare soil evaporation after Dickinson (1984).
!!
!!
SUBROUTINE calc_evsl_bats ( &
      & dt, ivstart, ivend, nvec, n_soil, z_ml, z_hl, dz_hl, soiltyp_subs, fr_w_i, fr_snow, eai, sai, &
      & evapot_s, w_so, fr_w_top, eva_bs, lhfl_bs, lzacc, acc_async_queue &
    )

  REAL(wp), INTENT(IN) :: dt !< Time step [s].
  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec               !< array dimensions
  INTEGER, INTENT(IN)  :: n_soil            !< Number of soil moisture layers.
  REAL(wp), INTENT(IN) :: z_ml(n_soil+1) !< Depth of main level (n_soil) [m].
  REAL(wp), INTENT(IN) :: z_hl(n_soil+1) !< Depth of layer interface (n_soil) [m].
  REAL(wp), INTENT(IN) :: dz_hl(n_soil+1) !< Layer thickness (n_soil) [m].

  !> Soil type (n_cell).
  INTEGER, INTENT(IN) :: soiltyp_subs(nvec)
  !> Surface fraction covered by interception water (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: fr_w_i(nvec)
  !> Surface fraction covered by snow (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: fr_snow(nvec)
  !> Evaporating surface area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: eai(nvec)
  !> Surface area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: sai(nvec)
  !> Potential water evaporation (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: evapot_s(nvec)
  !> Soil water content (n_cell,n_soil) [m(H2O)].
  REAL(wp), INTENT(IN) :: w_so(nvec,n_soil+1)
  !> Relative top soil water content (n_cell) [m(H2O)/m(dz)].
  REAL(wp), INTENT(IN) :: fr_w_top(nvec)
  !> Bare-soil evaporation flux (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: eva_bs(nvec)
  !> Bare-soil latent heat flux (n_cell) [W/m**2].
  REAL(wp), INTENT(OUT) :: lhfl_bs(nvec)

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  REAL(wp) :: z10cm !< Actual depth limit for 10cm depth.
  REAL(wp) :: z100cm !< Actual depth limit for 100cm depth.

  REAL(wp) :: w_so_10cm(nvec) !< Water content up to ~10cm depth.
  REAL(wp) :: w_so_100cm(nvec) !< Water content up to ~100cm depth.
  REAL(wp) :: rs_so_10cm !< Saturation relative to pore volume up to ~10cm depth.
  REAL(wp) :: rs_so_100cm !< Saturation relative to pore volume up to ~100cm depth.

  REAL(wp) :: d
  REAL(wp) :: ck
  REAL(wp) :: fqmax !< Maximum sustainable evaporation flux [kg/(m**2 s)].
  REAL(wp) :: evapor !< Tentative evaporation flux [kg/(m**2 s)].
  REAL(wp) :: beta !< Fraction of potential evaporation that gets evaporated.
  REAL(wp) :: fr_w_top_new !< Relative water content of surface layer [m(H2O)/m(dz)].
  REAL(wp) :: sfc_frac_bs !< Fraction of evaporating surface to total surface [m**2/m**2].

  INTEGER :: i, mstyp

  ! Auxiliary parameters
  REAL(wp), PARAMETER :: bf1(*) =  [( &
      & 5.5_wp - 0.8_wp * cbedi(mstyp) * &
      &   (1.0_wp + 0.1_wp * (cbedi(mstyp) - 4.0_wp) * clgk0(mstyp)), &
      & mstyp = 1, IST_NUM &
    )]
  REAL(wp), PARAMETER :: bf2(*) = [( &
      & (cbedi(mstyp) - 3.7_wp + 5.0_wp/MAX(cbedi(mstyp), eps_div)) / (5.0_wp + cbedi(mstyp)), &
      & mstyp = 1, IST_NUM &
    )]
  REAL(wp), PARAMETER :: dmax(*) = [( &
      & cbedi(mstyp) * cfinull * ck0di(mstyp) / crhowm, &
      & mstyp = 1, IST_NUM &
    )]

  !$ACC DATA PRESENT(ivend) CREATE(w_so_10cm, w_so_100cm) ASYNC(acc_async_queue)

  CALL get_water_in_100cm_10cm ( &
      & ivstart, ivend, nvec, n_soil, z_ml, z_hl, w_so, z100cm, w_so_100cm, z10cm, w_so_10cm, &
      & lzacc, acc_async_queue &
    )

  ! Calculation of bare soil evaporation after Dickinson (1984)
  ! Determination of mean water content relative to volume of voids

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp, beta, rs_so_100cm, rs_so_10cm, d, ck, fqmax) &
  !$ACC   PRIVATE(evapor, fr_w_top_new, sfc_frac_bs)
  DO i = ivstart, ivend
    IF (evapot_s(i) < 0._wp) THEN   ! upwards directed potential evaporation
      mstyp = soiltyp_subs(i)

      beta = 0._wp

      ! Treatment of ice and rocks
      SELECT CASE (mstyp)
      CASE (IST_ICE, IST_ROCK)
        beta = 0._wp
        eva_bs(i) = 0._wp
        lhfl_bs(i) = 0._wp

      CASE DEFAULT ! Computations not for ice and rocks
        ! auxiliary quantities
        rs_so_100cm = w_so_100cm(i) / (z100cm * cporv(mstyp))
        rs_so_10cm  = w_so_10cm(i) / (z10cm * cporv(mstyp))
        d = 1.02_wp * dmax(mstyp) * EXP( (cbedi(mstyp)+2._wp) * LOG(rs_so_10cm) ) * &
            & EXP( bf1(mstyp) * LOG(rs_so_100cm/rs_so_10cm) )
        ck = (1.0_wp + 1550.0_wp * cdmin / dmax(mstyp)) * bf2(mstyp)

        ! maximum sustainable moisture flux in the uppermost surface
        ! layer in kg/(s*m**2)
        fqmax = - rho_w * ck * d * rs_so_100cm / SQRT(z100cm * z10cm)
        evapor = MAX(evapot_s(i), fqmax)

        ! Stop evaporating if new relative water content would be below air dryness point.
        fr_w_top_new = fr_w_top(i) + &
            & evapor * (1.0_wp - fr_w_i(i)) * (1.0_wp - fr_snow(i)) * &
            & eai(i) / sai(i) * dt / rho_w / dz_hl(1)
        IF (fr_w_top_new <= cadp(mstyp)) evapor = 0._wp

        beta = evapor / MIN(evapot_s(i), -eps_div)
        sfc_frac_bs = eai(i) / sai(i)

        IF (mstyp == IST_PEAT .AND. itype_mire == 1) THEN ! AYu mire block
          beta = 0.6_wp
          sfc_frac_bs = 1.0_wp
        ENDIF

        eva_bs(i) = beta * evapot_s(i)         & ! evaporation
        !!!              *(1.0_wp - fr_w_i(i)) & ! not water covered
                        *(1.0_wp - fr_snow(i)) & ! not snow covered
                        * sfc_frac_bs            ! relative source surface
                                                 ! of the bare soil
        lhfl_bs(i) = lh_v * eva_bs(i)
      END SELECT
    ELSE
      eva_bs(i) = 0._wp
      lhfl_bs(i) = 0._wp
    END IF  ! upwards directed potential evaporation
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE calc_evsl_bats


!>
!! Bare-soil evaporation after Noilhan and Planton (1989).
!!
!!
SUBROUTINE calc_evsl_np89 ( &
      & ivstart, ivend, nvec, soiltyp_subs, fr_snow, eai, sai, rho_ch, qv_atm, evapot_s, &
      & t_snow_top, t_sk, t_s, p_s, fr_w_top, eva_bs, lhfl_bs, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec               !< array dimensions
  INTEGER, INTENT(IN)  :: soiltyp_subs(nvec)
  !> Surface fraction covered by snow (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: fr_snow(nvec)
  !> Evaporating surface area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: eai(nvec)
  !> Surface area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: sai(nvec)
  !> Surface air density times transfer velocity (n_cell) [kg(air)/(m**2 s)].
  REAL(wp), INTENT(IN) :: rho_ch(nvec)
  !> Specific humidity in lowest atm. layer (n_cell) [kg(vap)/kg(air)].
  REAL(wp), INTENT(IN) :: qv_atm(nvec)
  !> Potential water evaporation (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: evapot_s(nvec)
  !> Snow surface temperature (n_cell) [K].
  REAL(wp), INTENT(IN) :: t_snow_top(nvec)
  !> Skin temperature (n_cell) [K].
  REAL(wp), INTENT(IN) :: t_sk(nvec)
  !> Surface temperature (n_cell) [K].
  REAL(wp), INTENT(IN) :: t_s(nvec)
  !> Surface pressure (n_cell) [Pa].
  REAL(wp), INTENT(IN) :: p_s(nvec)
  !> Relative top soil water content (n_cell) [m(H2O)/m(dz)].
  REAL(wp), INTENT(IN) :: fr_w_top(nvec)
  !> Bare-soil evaporation flux (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: eva_bs(nvec)
  !> Bare-soil latent heat flux (n_cell) [W/m**2].
  REAL(wp), INTENT(OUT) :: lhfl_bs(nvec)

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  REAL(wp) :: b2iw !< Vapor pressure scale coefficient over ice or water.
  REAL(wp) :: b4iw !< Vapor pressure offset coefficient over ice or water.
  REAL(wp) :: qs !< Saturation specific humidity.

  REAL(wp) :: evapor !< Tentative evaporation flux [kg/(m**2 s)].
  REAL(wp) :: alpha !< Relative surface humidity.
  REAL(wp) :: beta !< Fraction of potential evaporation that gets evaporated.

  INTEGER :: i, mstyp

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(alpha, beta, b2iw, b4iw, qs, evapor)
  DO i = ivstart, ivend
    IF (evapot_s(i) < 0._wp) THEN   ! upwards directed potential evaporation
      mstyp = soiltyp_subs(i)

      beta = 0._wp

      ! Treatment of ice and rocks
      SELECT CASE (mstyp)
      CASE (IST_ICE, IST_ROCK)
        beta = 0._wp
        eva_bs(i) = 0._wp
        lhfl_bs(i) = 0._wp

      CASE DEFAULT ! Computations not for ice and rocks
        IF (fr_w_top(i)> cfcap(mstyp)) THEN
          alpha = 1.0_wp
        ELSE
          alpha = 0.5_wp * (1.0_wp - COS( &
              & pi * (fr_w_top(i) - cadp(mstyp)) / ( cfcap(mstyp) - cadp(mstyp))))
        ENDIF

        IF (itype_canopy == 2) THEN
          b2iw = MERGE(b2w, b2i, t_sk(i) >= t0_melt)
          b4iw = MERGE(b4w, b4i, t_sk(i) >= t0_melt)
          qs = zsf_qsat(zsf_psat_iw(t_sk(i), b2iw, b4iw), p_s(i))
        ELSE ! IF (itype_canopy == 1) THEN
          b2iw = MERGE(b2w, b2i, t_snow_top(i) >= t0_melt)
          b4iw = MERGE(b4w, b4i, t_snow_top(i) >= t0_melt)
          qs = zsf_qsat(zsf_psat_iw(t_s(i), b2iw, b4iw), p_s(i))
        END IF
        evapor = MIN(0._wp, rho_ch(i) * (qv_atm(i) - alpha * qs))

        beta = evapor / MIN(evapot_s(i),-eps_div)
        eva_bs(i) = beta * evapot_s(i)         & ! evaporation
        !!!            *(1.0_wp - fr_w_i (i)) & ! not water covered
                       *(1.0_wp - fr_snow(i)) & ! not snow covered
                       * eai(i)/sai(i)          ! relative source surface
                                                ! of the bare soil
        lhfl_bs(i) = lh_v * eva_bs(i)
      END SELECT ! Computations not for ice and rocks
    ELSE
      eva_bs(i) = 0._wp
      lhfl_bs(i) = 0._wp
    END IF  ! upwards directed potential evaporation
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE calc_evsl_np89


!>
!! Bare-soil evaporation, resistance formulation.
!!
!! For a review see Schulz et al. (1998) and Schulz and Vogel (2020).
!!
!!
SUBROUTINE calc_evsl_resistance ( &
      & ivstart, ivend, nvec, soiltyp_subs, fr_snow, eai, sai, rho_ch, rho_atm, r_bsmin, evapot_s, fr_w_top, &
      & eva_bs, lhfl_bs, lzacc, acc_async_queue &
    )

  !> Soil type (n_cell).
  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec               !< array dimensions
  INTEGER, INTENT(IN)  :: soiltyp_subs(nvec)
  !> Surface fraction covered by snow (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: fr_snow(nvec)
  !> Evaporating surface area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: eai(nvec)
  !> Surface area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: sai(nvec)
  !> Surface air density times transfer velocity (n_cell) [kg(air)/(m**2 s)].
  REAL(wp), INTENT(IN) :: rho_ch(nvec)
  !> Air density in lowest atm. layer (n_cell) [kg(air)/m**3].
  REAL(wp), INTENT(IN) :: rho_atm(nvec)
  !> Minimum bare-soil evaporation resistance (n_cell) [s/m].
  REAL(wp), INTENT(IN) :: r_bsmin(nvec)
  !> Potential water evaporation (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: evapot_s(nvec)
  !> Relative top soil water content (n_cell) [m(H2O)/m(dz)].
  REAL(wp), INTENT(IN) :: fr_w_top(nvec)
  !> Bare-soil evaporation flux (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: eva_bs(nvec)
  !> Bare-soil latent heat flux (n_cell) [W/m**2].
  REAL(wp), INTENT(OUT) :: lhfl_bs(nvec)

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  REAL(wp) :: alpha !< Relative surface humidity.
  REAL(wp) :: r_bs !< Resistance for bare-soil evaporation [s/m].
  REAL(wp) :: beta !< Fraction of potential evaporation that gets evaporated.
  REAL(wp) :: sfc_frac_bs !< Surface fraction of bare soil [m**2(soil)/m**2(surf)].

  INTEGER :: i, mstyp

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(alpha, beta, sfc_frac_bs)
  DO i = ivstart, ivend
    IF (evapot_s(i) < 0._wp) THEN   ! upwards directed potential evaporation
      mstyp = soiltyp_subs(i)

      beta = 0._wp

      SELECT CASE (mstyp)
      CASE (IST_ICE, IST_ROCK)
        beta = 0._wp
        eva_bs(i) = 0._wp
        lhfl_bs(i) = 0._wp

      CASE DEFAULT ! Computations not for ice and rocks
        alpha = MAX( 0.0_wp, MIN( &
            & 1.0_wp, &
            & (fr_w_top(i) - cadp(mstyp)) / (cfcap(mstyp) - cadp(mstyp))))
        r_bs = r_bsmin(i) / (alpha + eps_soil)
        beta = 1.0_wp / (1.0_wp + rho_ch(i) * r_bs / rho_atm(i))
        sfc_frac_bs = eai(i) / sai(i)

        IF (mstyp == IST_PEAT .AND. itype_mire == 1) THEN      ! AYu mire block
          beta = 0.6_wp
          sfc_frac_bs = 1.0_wp
        ENDIF

        eva_bs(i) = beta * evapot_s(i)        & ! evaporation
        !!!           * (1.0_wp - fr_w_i (i)) & ! not water covered
                      * (1.0_wp - fr_snow(i)) & ! not snow covered
                      * sfc_frac_bs             ! relative source surface
                                                ! of the bare soil
        lhfl_bs(i) = lh_v * eva_bs(i)
      END SELECT ! Computations not for ice and rocks
    ELSE
      eva_bs(i) = 0._wp
      lhfl_bs(i) = 0._wp
    END IF  ! upwards directed potential evaporation
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE calc_evsl_resistance


!>
!! Plant transpiration, BATS version.
!!
!! This version is based on Dickinson's (1984) BATS scheme, simplified by
!! neglecting the water and energy transports between the soil and the plant
!! canopy. This leads to a Monteith combination formula for the computation
!! of plant transpiration.
!! Option itype_trvg=TRVG_BATS_EXT (=3) is an extended variant with an
!! additional diagnostic variable for accumulated plant evaporation since
!! sunrise, allowing for a better representation of the diurnal cycle of plant
!! evaporation (in particular trees).
!!
!!
SUBROUTINE calc_trvg_bats ( &
    & dt, icant, ivstart, ivend, nvec, n_soil, n_soil_hy, soiltyp_subs, z_ml, dz_hl, fr_w_ml, &
    & w_so_ice, u_atm, v_atm, t_atm, tcm, tch, z0, r_stommin, sai, tai, plcov, laifac, fr_snow, &
    & evapot_s, par_absorbed, root_depth, plevap, transp_ml, transp_sum, lhfl_pl, r_stom, lzacc, &
    & acc_async_queue &
  )

  REAL(wp), INTENT(IN) :: dt !< Time step [s].
  INTEGER, INTENT(IN)  :: icant !< Canopy transfer scheme (1=Louis, 2=Raschendorfer, 3=EDMF).
  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec               !< array dimensions
  INTEGER, INTENT(IN) :: n_soil_hy !< Number of active soil moisture layers.
  INTEGER, INTENT(IN) :: n_soil !< Number of soil moisture layers.
  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type (n_cell).
  REAL(wp), INTENT(IN) :: z_ml(n_soil+1) !< Depth of main level (n_soil) [m].
  REAL(wp), INTENT(IN) :: dz_hl(n_soil+1) !< Layer thickness (n_soil) [m].
  !> Relative soil water content (n_cell, n_soil) [m(H2O)/m(dz)].
  REAL(wp), INTENT(IN) :: fr_w_ml(nvec,n_soil+1)
  !> Soil ice content (n_cell, n_soil) [m(H2O)].
  REAL(wp), INTENT(IN) :: w_so_ice(nvec,n_soil+1)
  REAL(wp), INTENT(IN) :: u_atm(nvec) !< Zonal wind in lowest layer (n_cell) [m/s].
  REAL(wp), INTENT(IN) :: v_atm(nvec) !< Meridional wind in lowest layer (n_cell) [m/s].
  REAL(wp), INTENT(IN) :: t_atm(nvec) !< Temperature in lowest layer (n_cell) [K].
  REAL(wp), INTENT(IN) :: tcm(nvec) !< Turbulent transfer coefficient for momentum (n_cell) [1].
  REAL(wp), INTENT(IN) :: tch(nvec) !< Turbulent transfer coefficient for heat (n_cell) [1].
  REAL(wp), INTENT(IN) :: z0(nvec) !< Roughness length (n_cell) [m].
  REAL(wp), INTENT(IN) :: r_stommin(nvec) !< Minimum stomata resistance (n_cell) [s/m].
  !> Surface area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: sai(nvec)
  REAL(wp), INTENT(IN) :: tai(nvec) !< Transpiration area index (n_cell) [m**2/m**2].
  REAL(wp), INTENT(IN) :: plcov(nvec) !< Plant-covered fraction (n_cell) [m**2(plants)/m**2(tile)].
  REAL(wp), INTENT(IN) :: laifac(nvec) !< Ratio between current LAI and laimax (n_cell) [1].
  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< Snow-covered fraction (n_cell) [m**2(snow)/m**2(tile)].
  !> Potential water evaporation (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: evapot_s(nvec)
  REAL(wp), INTENT(IN) :: par_absorbed(nvec) !< Absorbed PAR radiation (n_cell) [W/m**2].
  REAL(wp), INTENT(IN) :: root_depth(nvec) !< Root depth (n_cell) [m]

  REAL(wp), INTENT(INOUT) :: plevap(nvec) !< Accumulated plant evaporation since sunrise [kg/m**2].

  REAL(wp), INTENT(OUT) :: transp_ml(nvec,n_soil) !< Transpiration from soil layers (n_cell, n_soil) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: transp_sum(nvec) !< Total Transpiration (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(OUT) :: lhfl_pl(nvec,n_soil) !< Latent heat flux through plants from soil layer (n_cell,n_soil) [W/m**2].
  REAL(wp), INTENT(OUT) :: r_stom(nvec) !< Stomata resistance (n_cell) [s/m].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  REAL(wp) :: root_depth_reg !< Root depth, regularized [m].
  REAL(wp) :: root_density !< Root density decay scale [1/m].
  REAL(wp) :: fr_root !< Root fraction in layer [1].
  REAL(wp) :: rootdz !< Effective root depth in layer [m].
  REAL(wp) :: rootdz_int(nvec) !< Integrated root density, effective root depth [m].
  REAL(wp) :: wrootdz(nvec,n_soil) !< Root water content in layer [m(H2O)].
  REAL(wp) :: wrootdz_int(nvec) !< Integrated root water content, normalized by `rootdz_int` [m(H2O)/m].

  REAL(wp) :: tlpmwp !< Turgor-loss point minus plant wilting point.
  REAL(wp) :: f_rad !< Radiation factor.
  REAL(wp) :: f_tem !< Temperature factor.
  REAL(wp) :: f_wat !< Water availability factor.
  REAL(wp) :: f_sat !< Relative humidity factor.
  REAL(wp) :: f_tot !< Product of factors.
  REAL(wp) :: uv !< Wind speed in lowest level [m/s].
  REAL(wp) :: rla !< Laminar canopy resistance [s/m].
  REAL(wp) :: ustar !< Friction velocity [m/s].
  REAL(wp) :: cond_stom !< Stomatal conductance [m/s].
  REAL(wp) :: catm !< Atmospheric conductance [m/s].
  REAL(wp) :: rveg !< Vegetation resistance [s/m].

  REAL(wp) :: rf_plevap !< Reduction factor for limiting transpiration in the evening [1].
  REAL(wp) :: rf_plevap_rad !< Multiplicative factor increasing amount of PAR necessary for `f_rad = 1` [1].
  REAL(wp) :: rf_plevap_stom !< Multiplicative factor increasing minimum stomatal resistance [1].

  REAL(wp) :: traleav(nvec) !< Leaf transpiration [kg/(m**2 s)].
  REAL(wp) :: tr_frac !< Fraction of transpiration that gets its water from this layer [1].

  INTEGER :: mstyp
  INTEGER :: i
  INTEGER :: kso

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA ASYNC(acc_async_queue) IF(lzacc) &
  !$ACC   PRESENT(ivend) CREATE(rootdz_int, wrootdz, wrootdz_int, traleav)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  ! Root distribution

  IF (itype_root == ROOT_EXPONENTIAL) THEN
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      rootdz_int(i) = 0.0_wp
      wrootdz_int(i) = 0.0_wp
    END DO

    !$ACC LOOP SEQ
    DO kso = 1, n_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fr_root, rootdz, mstyp, root_depth_reg, root_density)
      DO i = ivstart, ivend
        mstyp = soiltyp_subs(i)
        IF (mstyp > IST_ROCK .AND. evapot_s(i) < 0.0_wp) THEN
          root_depth_reg = MAX(0.001_wp,root_depth(i))
          root_density = BATS_ROOT_DENSITY / root_depth_reg

          ! consider the effect of root depth & root density
          fr_root = EXP (-root_density*z_ml(kso)) ! root density
          rootdz = fr_root * MIN( &
              & dz_hl(kso), &
              & MAX(0.0_wp, root_depth_reg - (z_ml(kso) - 0.5_wp * dz_hl(kso))))
          rootdz_int(i) = rootdz_int(i) + rootdz

          ! The factor of 10 ensures that plants do not extract notable amounts of water from partly frozen soil
          wrootdz(i,kso) = rootdz * MAX( &
              & cpwp(mstyp), &
              & fr_w_ml(i,kso) - 10.0_wp * w_so_ice(i,kso) / dz_hl(kso))
          wrootdz_int(i) = wrootdz_int(i) + wrootdz(i,kso)
        END IF  ! negative potential evaporation only
      END DO
    END DO

    ! Compute root zone integrated average of liquid water content
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      wrootdz_int(i) = wrootdz_int(i) / MAX(rootdz_int(i), eps_div)
    END DO

  ELSE   ! itype_root

    !$ACC LOOP SEQ
    DO kso = 1, n_soil
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(rootdz, mstyp, root_depth_reg)
      DO i = ivstart, ivend
        mstyp = soiltyp_subs(i)
        IF (mstyp > IST_ROCK .AND. evapot_s(i) < 0.0_wp) THEN
          root_depth_reg = MAX(0.001_wp,root_depth(i))

          rootdz = MIN( &
              & dz_hl(kso), &
              & MAX(0.0_wp, root_depth_reg - (z_ml(kso) - 0.5_wp*dz_hl(kso))))
          wrootdz(i,kso) = rootdz * (fr_w_ml(i,kso) - w_so_ice(i,kso) / dz_hl(kso))
          wrootdz_int(i) = wrootdz_int(i) + wrootdz(i,kso)
        END IF  ! negative potential evaporation only
      END DO
    END DO

    ! Normalize root-zone integrated water
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      rootdz_int(i) = MAX(0.001_wp, root_depth(i))
      wrootdz_int(i) = wrootdz_int(i) / rootdz_int(i)
    END DO

  ENDIF  ! itype_root


  ! Determination of the transfer functions CA, CF, and CV
  !$ACC LOOP GANG(STATIC: 1) VECTOR &
  !$ACC   PRIVATE(mstyp, uv, catm, ustar, rla, f_rad, f_wat, f_tem, f_sat, f_tot, rf_plevap) &
  !$ACC   PRIVATE(rf_plevap_rad, rf_plevap_stom, cond_stom, rveg)
  DO i = ivstart, ivend
    mstyp = soiltyp_subs(i)

    ! Zero out layer-accumulated transpiration.
    transp_sum(i) = 0._wp

    IF (mstyp > IST_ROCK .AND. evapot_s(i) < 0.0_wp) THEN
      ! upwards directed potential evaporation
      uv = SQRT (u_atm(i)**2 + v_atm(i)**2 )
      catm = tch(i)*uv           ! Function CA

      SELECT CASE ( icant )
        CASE (1)   ! Louis-transfer-scheme: additional laminar canopy resistance
          ustar = uv*SQRT(tcm(i))
          rla = 1.0_wp/MAX(cdash*SQRT(ustar),eps_div)
        CASE (2)   ! Raschendorfer transfer scheme: laminar canopy resistance already considered
          rla = 0.0_wp
        CASE DEFAULT ! additional laminar canopy resistance for other schemes
          rla = 1.0_wp/MAX(catm,eps_div)
      END SELECT

      ! to compute CV, first the stomatal resistance has to be determined
      ! this requires the determination of the F-functions:
      ! Radiation function
      IF (itype_trvg == TRVG_BATS_EXT) THEN
        ! modification depending on accumulated plant evaporation in order to reduce evaporation in the evening
        rf_plevap = 0.75_wp * MAX(eps_soil, ABS(plevap(i))) / MAX(0.2_wp, plcov(i))
        rf_plevap_rad = MIN(3._wp, MAX(1._wp, rf_plevap))
        ! stronger limitation for non-forest vegetation classes
        IF (z0(i) <= 0.4_wp) rf_plevap_rad = MIN(2._wp, rf_plevap_rad)
      ELSE
        rf_plevap_rad = 1._wp
        rf_plevap = 1._wp
      ENDIF
      f_rad = MAX(0._wp, MIN(1._wp, par_absorbed(i) / (cparcrit*rf_plevap_rad)))
      tlpmwp = (cfcap(mstyp) - cpwp(mstyp)) * &
          & (0.81_wp + 0.121_wp * ATAN(-86400._wp * evapot_s(i) - 4.75_wp))

      ! Soil water function
      f_wat = MAX(0._wp, MIN(1._wp,(wrootdz_int(i) - cpwp(mstyp))/tlpmwp))

      ! Temperature function
      ! T at lowest model level used (approximation of leaf height)
      f_tem = MAX(0._wp, MIN(1._wp, &
          & 4._wp * (t_atm(i)-t0_melt)*(ctend-t_atm(i))/(ctend-t0_melt)**2))

      ! f_sat fixed.
      f_sat = 1._wp

      f_tot = f_rad * f_wat * f_tem * f_sat

      IF (lstomata) THEN
        IF (itype_trvg == TRVG_BATS_EXT) THEN
          ! Modification of rsmin depending on accumulated plant evaporation; the z0 dependency
          ! is used to get a stronger effect for trees than for low vegetation

          rf_plevap_stom = MAX(0.5_wp+MIN(0.5_wp,1.0_wp-laifac(i)), EXP(SQRT(z0(i))*LOG(rf_plevap)) )
        ELSE
          rf_plevap_stom = 1.0_wp
        ENDIF
        cond_stom = (1._wp - f_tot)/crsmax + f_tot/MAX(40._wp,rf_plevap_stom*r_stommin(i))
      ELSE
        cond_stom = (1._wp - f_tot)/crsmax + f_tot/crsmin
      END IF

      r_stom(i) = 1._wp/cond_stom
      rveg = rla + r_stom(i)

      ! Transpiration rate of dry leaves:
      traleav(i) = evapot_s(i) * tai(i) / (sai(i) + rveg * catm)
    ELSE
      r_stom(i) = 0._wp
    END IF  ! upwards directed potential evaporation only
  END DO

  ! Consideration of water and snow coverage, distribution to the different
  ! soil layers

  !$ACC LOOP SEQ
  DO kso = 1, n_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp, tr_frac)
    DO i = ivstart, ivend
      mstyp = soiltyp_subs(i)

      IF (mstyp > IST_ROCK .AND. evapot_s(i) < 0.0_wp) THEN
        ! upwards potential evaporation

        tr_frac = wrootdz(i,kso) / (rootdz_int(i) * wrootdz_int(i))
        transp_ml(i,kso) = traleav(i) * tr_frac * & ! plant covered part
        !!!            (1.0_wp - fr_w_i(i)) * & ! not water covered
                       (1.0_wp - fr_snow(i))    ! not snow covered

        ! Limit evaporation such that the soil water content does not fall beyond the wilting point
        IF(fr_w_ml(i,kso) + transp_ml(i,kso) * dt / rho_w / dz_hl(kso) < cpwp(mstyp)) THEN
          transp_ml(i,kso) = MIN(0._wp, (cpwp(mstyp)-fr_w_ml(i,kso)) * dz_hl(kso) * rho_w / dt)
        END IF

        IF (soiltyp_subs(i) == IST_PEAT .AND. itype_mire == 1) THEN    ! AYu mire block
          transp_ml(i,kso) = 0.0_wp
        ENDIF

        lhfl_pl(i,kso) = lh_v * transp_ml(i,kso)
        transp_sum(i) = transp_sum(i) + transp_ml(i,kso)
      ELSE
        lhfl_pl(i,kso) = 0._wp
        transp_ml(i,kso) = 0._wp
        transp_sum(i) = 0._wp
      END IF  ! upwards directed potential evaporation only
    END DO
  END DO          ! loop over soil layers
  !$ACC END PARALLEL
  !$ACC END DATA

  IF (itype_trvg == TRVG_BATS_EXT) THEN
    CALL update_plevap ( &
        & dt, ivstart, ivend, nvec, n_soil, n_soil_hy, plcov, par_absorbed, lhfl_pl, plevap, &
        & lzacc, acc_async_queue &
      )
  END IF

END SUBROUTINE calc_trvg_bats


!>
!! Update accumulated plant evaporation since sunrise.
!!
!! The latent heat flux is offset by 75 W/m**2(plants) parametrizing the amount of water that can
!! be continuously supplied through the stem. That offset also drives the reset for the next day.
!! To ensure the value is reset, the recovery rate is increased at night
!! (absorbed PAR < 40 W/m**2).
!!
!!
SUBROUTINE update_plevap ( &
      & dt, ivstart, ivend, nvec, n_soil, n_soil_hy, plcov, par_absorbed, lhfl_pl, plevap, &
      & lzacc, acc_async_queue &
    )

  REAL(wp), INTENT(IN) :: dt !< Time step [s].
  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec               !< array dimensions
  INTEGER, INTENT(IN) :: n_soil !< Number of soil moisture layers.
  INTEGER, INTENT(IN) :: n_soil_hy !< Number of active soil moisture layers.
  REAL(wp), INTENT(IN) :: plcov(nvec) !< Plant-covered fraction (n_cell) [m**2(plants)/m**2(tile)].
  REAL(wp), INTENT(IN) :: par_absorbed(nvec) !< Absorbed PAR radiation (n_cell) [W/m**2].
  REAL(wp), INTENT(IN) :: lhfl_pl(nvec,n_soil) !< Latent heat flux through plants from soil layer (n_cell,n_soil) [W/m**2].
  REAL(wp), INTENT(INOUT) :: plevap(nvec) !< Accumulated plant evaporation since sunrise (n_cell) [kg/m**2].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  INTEGER :: i
  REAL(wp) :: recovery_rate

#ifdef __SX__
  REAL(wp) :: lhfl_pl_int(SIZE(plevap))
  INTEGER :: kso
#endif

  OPENACC_SUPPRESS_UNUSED_LZACC

  IF (itype_trvg /= TRVG_BATS_EXT) RETURN

#ifdef __SX__
  lhfl_pl_int(:) = 0._wp
  DO kso = 1, n_soil_hy
    DO i = ivstart, ivend
      lhfl_pl_int(i) = lhfl_pl_int(i) + lhfl_pl(i, kso)
    ENDDO
  ENDDO
#endif

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(recovery_rate)
  DO i = ivstart, ivend
    recovery_rate = MAX(1._wp, 0.1_wp*(50._wp - par_absorbed(i)))
    plevap(i) = MAX(-6._wp, MIN(0._wp, plevap(i) + recovery_rate*dt/lh_v * &
#ifdef __SX__
                   (lhfl_pl_int(i)+MAX(0.2_wp,plcov(i))*75._wp) ))
#else
                   (SUM(lhfl_pl(i,1:n_soil_hy))+MAX(0.2_wp,plcov(i))*75._wp) ))
#endif
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE update_plevap

!>
!! Limit evaporation from all sources to potential evaporation.
!!
!! In case the sum of evaporations exceeds potential evaporation, a suppression factor is applied that
!! uniformly reduces evaporation from all sources. Also stops bare-soil evaporation and plant transpiration
!! for snow-covered tiles.
!!
!!
SUBROUTINE limit_evaporation ( &
      & ivstart, ivend, nvec, n_soil, evapot_s, evapot_snow, t_snred, fr_snow, eva_w_i, eva_w_sn, &
      & eva_bs, lhfl_bs, transp_sum, transp_ml, lhfl_pl, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec
  INTEGER, INTENT(IN)  :: n_soil
  !> Potential water evaporation (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: evapot_s(nvec)
  !> Potential evaporation for snow (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: evapot_snow(nvec)
  !> Snow temperature offset for calculating evaporation (n_cell) [K].
  REAL(wp), INTENT(IN) :: t_snred(nvec)
  !> Snow-covered fraction (n_cell) [m**2(snow)/m**2(tile)].
  REAL(wp), INTENT(IN) :: fr_snow(nvec)

  !> Evaporation from interception water (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(INOUT) :: eva_w_i(nvec)
  !> Evaporation from snow (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(INOUT) :: eva_w_sn(nvec)
  !> Bare-soil evaporation (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(INOUT) :: eva_bs(nvec)
  !> Bare-soil latent heat flux (n_cell) [W/m**2].
  REAL(wp), INTENT(INOUT) :: lhfl_bs(nvec)
  !> Total transpiration (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(INOUT) :: transp_sum(nvec)
  !> Transpiration from soil layers (n_cell, n_soil) [kg/(m**2 s)].
  REAL(wp), INTENT(INOUT) :: transp_ml(nvec,n_soil)
  !> Latent heat flux through plants from soil layer (n_cell, n_soil) [W/m**2].
  REAL(wp), INTENT(INOUT) :: lhfl_pl(nvec,n_soil)

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  REAL(wp) :: eva_sum !< Sum of evapotranspiration contributions [kg/(m**2 s)].
  REAL(wp) :: evapot_wgt !< Snow-cover weighted potential evaporation [kg/(m**2 s)].
  REAL(wp) :: eva_red !< Reduction factor for evaporation.
  REAL(wp) :: tran_red !< Reduction factor for transpiration.

  INTEGER :: i

#ifdef __SX__
  !> Accumulated reduction factors for transpiration (NEC only).
  REAL(wp) :: tran_fac(SIZE(evapot_s))
  INTEGER :: kso
#endif

  OPENACC_SUPPRESS_UNUSED_LZACC

#ifdef __SX__
  tran_fac(:) = 1._wp
#endif

  ! Ensure that the sum of the evaporation terms does not exceed the potential evaporation
  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(eva_sum, evapot_wgt, eva_red, tran_red)
  DO i = ivstart, ivend
    eva_sum = eva_w_sn(i) + eva_w_i(i) + eva_bs(i) + transp_sum(i)
    ! snow-weighted potential evaporation
    evapot_wgt = fr_snow(i) * evapot_snow(i) + (1._wp - fr_snow(i)) * evapot_s(i)

    IF (evapot_wgt < 0._wp .AND. eva_sum < evapot_wgt) THEN
      eva_red = evapot_wgt / eva_sum
      eva_w_sn(i) = eva_w_sn(i) * eva_red
      eva_w_i(i) = eva_w_i(i) * eva_red
      eva_bs(i) = eva_bs(i) * eva_red
      transp_sum(i) = transp_sum(i) * eva_red
      lhfl_bs(i) = lhfl_bs(i) * eva_red
#ifdef __SX__
      tran_fac(i) = eva_red
#else
      transp_ml(i,:) = transp_ml(i,:) * eva_red
      lhfl_pl(i,:) = lhfl_pl(i,:) * eva_red
#endif
    ENDIF

    ! Negative values of t_snred indicate that snow is present on the corresponding snow tile
    ! and that the snow-free tile has been artificially generated by the melting-rate parameterization
    ! in this case, bare soil evaporation and, in the case of a long-lasting snow cover, plant evaporation, are turned off.
    IF (t_snred(i) < 0.0_wp .AND. evapot_s(i) < 0.0_wp) THEN
      eva_red = MAX(0.0_wp,1.0_wp-ABS(t_snred(i)))
      tran_red = MIN(1.0_wp,MAX(0.0_wp,2.0_wp-ABS(t_snred(i))))
      eva_bs(i) = eva_bs(i) * eva_red
      lhfl_bs(i) = lhfl_bs(i) * eva_red
      transp_sum(i) = transp_sum(i) * tran_red
#ifdef __SX__
      tran_fac(i) = tran_fac(i) * tran_red
#else
      transp_ml(i,:) = transp_ml(i,:) * tran_red
      lhfl_pl(i,:) = lhfl_pl(i,:) * tran_red
#endif
    ENDIF
  ENDDO
  !$ACC END PARALLEL
  !$ACC END DATA

#ifdef __SX__
  !$NEC unroll(7)
  DO kso = 1, n_soil
    transp_ml(:,kso) = transp_ml(:,kso) * tran_fac(:)
    lhfl_pl(:,kso) = lhfl_pl(:,kso) * tran_fac(:)
  ENDDO
#endif

END SUBROUTINE limit_evaporation


!>
!! Retrieve the amount of water in the first 10cm and 100cm of soil.
!!
!! Looks for the first layer deeper than 10cm (100cm) and adds up water contents in that layer and
!! above. Returns the actual depths of the lower boundaries of the selected layers.
!!
!!
SUBROUTINE get_water_in_100cm_10cm ( &
      & ivstart, ivend, nvec, n_soil, z_ml, z_hl, w_so, z100cm, w_so_100cm, z10cm, w_so_10cm, &
      & lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN)  :: ivstart !< start index for computations in the parallel program
  INTEGER, INTENT(IN)  :: ivend !< end index for computations in the parallel program
  INTEGER, INTENT(IN)  :: nvec
  INTEGER, INTENT(IN)  :: n_soil
  REAL(wp), INTENT(IN) :: z_ml(n_soil+1) !< Depth of main level (n_soil) [m].
  REAL(wp), INTENT(IN) :: z_hl(n_soil+1) !< Depth of layer interface (n_soil) [m].
  REAL(wp), INTENT(IN) :: w_so(nvec,n_soil+1) !< Water content in layer (n_cell,n_soil) [m(H2O)].
  REAL(wp), INTENT(OUT) :: z100cm !< Actual depth of 100cm layer [m].
  REAL(wp), INTENT(OUT) :: w_so_100cm(nvec) !< Water in 100cm layer (n_cell) [m(H2O)].
  REAL(wp), INTENT(OUT), OPTIONAL :: z10cm !< Actual depth of 10cm layer [m].
  REAL(wp), INTENT(OUT), OPTIONAL :: w_so_10cm(nvec) !< Water in 10cm layer (n_cell) [m(H2O)].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  INTEGER :: k10cm, k100cm, kso
  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)

  ! Find soil levels closest to 10cm and 100cm depth
  IF (PRESENT(z10cm)) THEN
    k10cm = FINDLOC(z_ml(:) <= 0.1_wp, .TRUE., DIM=1, BACK=.TRUE.)
    z10cm = z_hl(k10cm)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      w_so_10cm(i) = 0._wp
    END DO
    !$ACC END PARALLEL
  ELSE
    k10cm = 0
  END IF

  k100cm = FINDLOC(z_ml(:) <= 1._wp, .TRUE., DIM=1, BACK=.TRUE.)
  z100cm = z_hl(k100cm)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    w_so_100cm(i) = 0._wp
  END DO

  ! Determine total water content in 10cm and 100cm.
  !$ACC LOOP SEQ
  DO kso = 1, k100cm
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      ! If w_so_10cm is not present, k10cm = 0 and the condition is never true.
      IF (kso <= k10cm) w_so_10cm(i) = w_so_10cm(i) + w_so(i,kso)
      w_so_100cm(i) = w_so_100cm(i) + w_so(i,kso)
    ENDDO
  ENDDO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE get_water_in_100cm_10cm


REAL (KIND=wp) PURE FUNCTION zsf_psat_iw  (zstx, z2iw, z4iw)
  ! Saturation water vapour pressure over ice or water depending on temperature "zstx"

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: zstx, z2iw, z4iw
  zsf_psat_iw   = b1*EXP(z2iw*(zstx - b3)/(zstx - z4iw))

END FUNCTION zsf_psat_iw


REAL (KIND=wp) PURE FUNCTION zsf_qsat  (zspsatx, zspx)
  ! Specific humidity at saturation pressure (depending on the saturation water
  !  vapour pressure zspsatx" and the air pressure "zspx")

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: zspsatx, zspx
  zsf_qsat      = rdv*zspsatx/(zspx-o_m_rdv*zspsatx)

END FUNCTION zsf_qsat


REAL (KIND=wp) PURE FUNCTION zsf_dqvdt_iw (zstx, zsqsatx, z4iw, z234iw)
  ! First derivative of specific saturation humidity with respect to temperature
  ! (depending on temperature "zstx" and saturation specific humidity pressure
  !  "zsqsatx")

  !$ACC ROUTINE SEQ
  REAL (KIND=wp), INTENT(IN)  :: zstx, zsqsatx, z4iw, z234iw
  zsf_dqvdt_iw  = z234iw*(1.0_wp+rvd_m_o*zsqsatx)*zsqsatx/(zstx-z4iw)**2

END FUNCTION zsf_dqvdt_iw

END MODULE sfc_terra_evap
