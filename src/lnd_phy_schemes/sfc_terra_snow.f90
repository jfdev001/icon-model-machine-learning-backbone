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

!> Implementation of the TERRA single- and multilayer snow models.
MODULE sfc_terra_snow

  USE mo_kind, ONLY: wp

  USE mo_physical_constants, ONLY: &
      & cp_d => cpd, &
      & g => grav, &
      & lh_f => alf, &
      & lh_s => als, &
      & rho_w => rhoh2o, &
      & t0_melt => tmelt

  USE mo_lnd_nwp_config, ONLY: &
      & l2lay_rho_snow, &
      & max_toplaydepth

  USE sfc_terra_data, ONLY: &
      & ca2, &
      & cdsmin, &
      & cfcap, &
      & chc_i, &
      & chc_w, &
      & chcond, &
      & cporv, &
      & crhosmax_ml, &
      & crhosmax_tmin, &
      & crhosmaxf, &
      & crhosmaxt, &
      & crhosmin_ml, &
      & crhosminf, &
      & crhosmint, &
      & crhogmaxf, &
      & crhogminf, &
      & csigma, &
      & csnow_tmin, &
      & cwhc, &
      & eps_div, &
      & eps_nounderflow, &
      & eps_temp, &
      & eps_soil, &
      & IST_ICE, &
      & IST_ROCK, &
      & rho_i

  USE sfc_terra_util, ONLY: &
      & solve_tridiag, &
      & zalfa

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: snow_update_freshsnow_factor
  PUBLIC :: snow_calc_precipitation_phase_change

  PUBLIC :: snow_single_prepare
  PUBLIC :: snow_single_soil_forcing
  PUBLIC :: snow_single_melt
  PUBLIC :: snow_single_calc_temperature
  PUBLIC :: snow_single_update_new_state

  PUBLIC :: snow_multi_prepare
  PUBLIC :: snow_multi_handle_snowfall
  PUBLIC :: snow_multi_soil_forcing
  PUBLIC :: snow_multi_calc_heat_conduction
  PUBLIC :: snow_multi_melt
  PUBLIC :: snow_multi_update_new_state

  ! Silence unused variable warnings on non-OpenACC builds
#ifdef _OPENACC
# define OPENACC_SUPPRESS_UNUSED_LZACC
#else
# define OPENACC_SUPPRESS_UNUSED_LZACC IF (lzacc .AND. acc_async_queue > 0) THEN; END IF
#endif

  !> Maximum snow depth for heat transfer calculations (single-layer scheme) [m].
  REAL(wp), PARAMETER :: hmax_single_heattransfer = 1.5_wp

CONTAINS

!>
!! Update indicator for age of snow in top of snow layer.
!!
!! Note that cloud ice is deliberately disregarded here in order to avoid counting drifting snow as
!! fresh snow.
!!
!!
SUBROUTINE snow_update_freshsnow_factor ( &
      & ivstart, ivend, nvec, dt, w_snow, t_snow, h_snow_gp, t_atm, sp_10m, rain_rate, snow_rate, &
      & freshsnow, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart !< Start cell index.
  INTEGER, INTENT(IN) :: ivend !< End cell index
  INTEGER, INTENT(IN) :: nvec !< Cell dimension.
  REAL(wp), INTENT(IN) :: dt !< Time step [s].
  REAL(wp), INTENT(IN) :: w_snow(nvec) !< Snow amount (n_cell) [m(H2O)].
  REAL(wp), INTENT(IN) :: t_snow(nvec) !< Snow temperature (n_cell) [K].
  REAL(wp), INTENT(IN) :: h_snow_gp(nvec) !< Grid-point averaged snow height (n_class) [m].
  REAL(wp), INTENT(IN) :: t_atm(nvec) !< Temperature in the lowest atmospheric level (n_cell) [K].
  REAL(wp), INTENT(IN) :: sp_10m(nvec) !< 10m wind speed (n_cell) [m/s].
  REAL(wp), INTENT(IN) :: rain_rate(nvec) !< Total rain rate (n_cell) [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: snow_rate(nvec) !< Total snow rate excluding ice (n_cell) [kg/(m**2 s)].

  REAL(wp), INTENT(INOUT) :: freshsnow(nvec) !< Freshsnow factor (n_cell) [kg/(m**2 s)].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  REAL(wp) :: tau !< Freshsnow decay time scale [s].
  REAL(wp) :: speed_factor !< Wind-speed factor [m**2/s**2].
  REAL(wp) :: delta_decay
  REAL(wp) :: delta_newsnow

  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR &
  !$ACC   PRIVATE(tau, speed_factor, delta_decay, delta_newsnow)
  DO i = ivstart, ivend
    IF (w_snow(i) <= 0._wp) THEN
      ! if no snow exists, reinitialize age indicator
      freshsnow(i) = 1._wp
    ELSE

      ! temperature-dependent aging timescale: 7 days at freezing point, 28 days below -15 deg C
      tau = 86400._wp * MIN(28.0_wp, 7._wp + 1.4_wp * (t0_melt - MIN(t0_melt, t_snow(i))))

      ! wind-dependent snow aging: a thin snow cover tends to get broken under strong winds, which reduces the albedo
      ! an offset is added in order to ensure moderate aging for low snow depths
      speed_factor = MIN(300._wp, sp_10m(i)**2 + 12._wp )
      tau = MIN(tau, MAX(86400._wp, 2.e8_wp * MAX(0.05_wp, h_snow_gp(i)) / speed_factor))

      ! decay rate for fresh snow including contribution by rain (full aging after 10 mm of rain)
      delta_decay = dt/tau + dt * rain_rate(i) * 0.1_wp

      ! linear growth rate equals 1.0 in 1 day for a temperature-dependent snow rate between
      ! 10 mmH2O (kg/m**2) per day (0.1) and 5 mmH2O (kg/m**2) per day (0.2)
      delta_newsnow = dt * snow_rate(i) * (0.1_wp + MIN(0.1_wp, 0.02_wp * (t0_melt - t_atm(i))))

      ! reduce decay rate, if new snow is falling and as function of snow
      ! age itself
      delta_decay = (delta_decay - delta_newsnow) * freshsnow(i)
      delta_decay = MAX(delta_decay, 0._wp)

      freshsnow(i) = freshsnow(i) + delta_newsnow - delta_decay
      freshsnow(i) = MIN(1._wp, MAX(0._wp, freshsnow(i)))
    END IF
  ENDDO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_update_freshsnow_factor

!>
!! Compute immediate melting / freezing of precipitation.
!!
!! Handles snow falling on warm soil or rain on cold soil. The resulting latent heat flux is stored
!! and enters into the soil forcing later. Melting snow is moved from snow to the interception
!! reservoir. Freezing rain was added to the interception reservoir and infiltration by
!! calc_infiltration.
!!
SUBROUTINE snow_calc_precipitation_phase_change ( &
      & ivstart, ivend, nvec, dt, rain_dew_rate, snow_rime_rate, w_snow_now, t_snow_top, &
      & soiltyp_subs, t_so_now, hcap_ml, dz_hl, w_i_now, w_i_max, fr_w_ml_top, dt_w_i, dt_w_snow, &
      & dt_w_so_top, runoff_s, lhfl_precip, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart !< Start cell index.
  INTEGER, INTENT(IN) :: ivend !< End cell index
  INTEGER, INTENT(IN) :: nvec !< Cell dimension.

  REAL(wp), INTENT(IN) :: dt !< Time step [s].
  REAL(wp), INTENT(IN) :: rain_dew_rate(nvec) !< Combined rainfall and dew-formation rate [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: snow_rime_rate(nvec) !< Combined snowfall and rime-formation rate [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current SWE [m(H2O)].
  REAL(wp), INTENT(IN) :: t_snow_top(nvec) !< Snow (soil) top temperature if snow is (not) present [K].
  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type.
  REAL(wp), INTENT(IN) :: t_so_now(nvec, 2) !< Soil layer temperature (top 2 layers) [K].
  REAL(wp), INTENT(IN) :: hcap_ml(nvec, 2) !< Volumetric soil heat capacity (top 2 layers) [J/(m^2 K)].
  REAL(wp), INTENT(IN) :: dz_hl(2) !< Layer thickness (top 2 layers) [m].
  REAL(wp), INTENT(IN) :: w_i_now(nvec) !< Interception water [m(H2O)].
  REAL(wp), INTENT(IN) :: w_i_max(nvec) !< Interception reservoir capacity [m(H2O)].
  REAL(wp), INTENT(IN) :: fr_w_ml_top(nvec) !< Fractional total water content in top layer [m(H2O)/m].

  REAL(wp), INTENT(INOUT) :: dt_w_i(nvec) !< Tendency of interception water [kg/(m^2 s)].
  REAL(wp), INTENT(INOUT) :: dt_w_snow(nvec) !< Tendency of total snow [kg/(m^2 s)].
  REAL(wp), INTENT(INOUT) :: dt_w_so_top(nvec) !< Tendency of top-layer soil water [kg/(m^2 s)].
  REAL(wp), INTENT(INOUT) :: runoff_s(nvec) !< Surface runoff [kg/m^2].

  REAL(wp), INTENT(INOUT) :: lhfl_precip(nvec) !< Heat flux due to precipitation phase change [W/m^2].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.


  INTEGER :: i
  REAL(wp) :: heat_melt_max !< Maximum heat available for immediately melting falling snow [J/m^2].
  REAL(wp) :: fr_melt !< Fraction of precipitation that melts immediately.
  REAL(wp) :: w_i_prov !< Provisional updated interception water [m(H2O)].
  REAL(wp) :: overflow_rate !< Interception reservoir overflow rate [kg/(m^2 s)].

  REAL(wp), PARAMETER :: t_so_min = t0_melt-0.25_wp !< Minimum soil temperature after melting [K].

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

#ifdef __INTEL_COMPILER
  !DIR$ NOFMA
#endif

  !NEC$ sparse
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(heat_melt_max, fr_melt, w_i_prov, overflow_rate)
  DO i = ivstart, ivend
    lhfl_precip(i) = 0.0_wp

    ! thawing of snow falling on soil with Ts > T0
    !
    ! it turned out that ensuring (snow_rime_rate(i) > 0.0_wp) is not sufficient, as in rare cases
    ! a floating point underflow might occur for snow_rime_rate. By itself this is not a huge problem, however computing
    ! 1/snow_rime_rate(i) (as it is done below) results in a floating point overflow in these cases.
    IF (t_snow_top(i) >= t0_melt .AND. snow_rime_rate(i) > eps_nounderflow) THEN
      ! snow fall on soil with T>T0, snow water content increases interception store water content
      ! melting rate is limited such that the two upper soil levels are not cooled significantly below the freezing point
      heat_melt_max = &
          & (t_so_now(i,1) - t_so_min) * hcap_ml(i,1) * dz_hl(1) + &
          & MAX(0._wp, t_so_now(i,2) - t_so_min) * hcap_ml(i,2) * dz_hl(2)
      fr_melt = heat_melt_max / MAX(dt * lh_f * snow_rime_rate(i), heat_melt_max)
      lhfl_precip(i) = - lh_f * snow_rime_rate(i) * fr_melt
      dt_w_i(i) = dt_w_i (i) + snow_rime_rate(i) * fr_melt
      dt_w_snow(i) = dt_w_snow(i) - snow_rime_rate(i) * fr_melt

      ! avoid overflow of interception store, add possible excess to
      ! surface run-off
      w_i_prov = w_i_now(i) + dt_w_i(i)*dt/rho_w
      IF (w_i_prov > w_i_max(i)) THEN  ! overflow of interception store
        overflow_rate = (w_i_prov - w_i_max(i))*rho_w/dt
        dt_w_i(i)   = dt_w_i(i) - overflow_rate
        dt_w_so_top(i) = dt_w_so_top(i) + overflow_rate
        ! check for pore volume overflow and add the remainder to surface runoff
        overflow_rate = MAX( &
            & 0.0_wp, &
            & (fr_w_ml_top(i) - cporv(soiltyp_subs(i))) * dz_hl(1) * rho_w / dt + dt_w_so_top(i) &
          )
        dt_w_so_top(i) = dt_w_so_top(i) - overflow_rate
        runoff_s(i) = runoff_s(i) + overflow_rate * dt
      ENDIF                       ! overflow of interception store

    ! freezing of rain falling on soil with Ts < T0  (black-ice !!!)
    ELSEIF (w_snow_now(i) == 0.0_wp .AND. t_snow_top(i) < t0_melt .AND. rain_dew_rate(i) > 0.0_wp) THEN
      lhfl_precip  (i) = MIN( &
          & lh_f*rain_dew_rate(i), &
          & (t0_melt - t_snow_top(i)) * hcap_ml(i,1) * dz_hl(1) / dt &
        )
      ! keep freezing rain in interception storage rather than shifting it to snow
      ! dt_w_i (i) = dt_w_i (i) - rain_dew_rate(i)
      ! dt_w_snow(i) = dt_w_snow(i) + rain_dew_rate(i)
    END IF
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_calc_precipitation_phase_change

!--------------------------------------------------------------------------------------------------
! Single-Layer Snow Model
!--------------------------------------------------------------------------------------------------

!>
!! Perform sanity checks on prognostic snow variables and compute effective snow height and fraction.
!!
SUBROUTINE snow_single_prepare ( &
      & ivstart, ivend, nvec, w_snow_now, t_snow_now, rho_snow_now, rho_snow_mult_now_top, &
      & t_s_now, freshsnow, fr_snow, h_snow_new, h_snow_now, dz_snow_flx, fr_snow_lim, rho_snow, &
      & hcap_snow, t_snow_top, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec

  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current SWE [m(H2O)].
  REAL(wp), INTENT(INOUT) :: t_snow_now(nvec) !< Current snow top temperature [K].
  REAL(wp), INTENT(IN) :: rho_snow_now(nvec) !< Current snow density [kg/m^3].
  REAL(wp), INTENT(IN) :: rho_snow_mult_now_top(nvec) !< Current snow density in top layer [kg/m^3].

  REAL(wp), INTENT(INOUT) :: t_s_now(nvec) !< Current soil top temperature [K].

  REAL(wp), INTENT(IN) :: freshsnow(nvec) !< Freshsnow factor [1].
  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< Snow fraction [m^2(snow)/m^2(tile)].

  REAL(wp), INTENT(INOUT) :: h_snow_new(nvec) !< New geometric snow height [m].
  REAL(wp), INTENT(INOUT) :: h_snow_now(nvec) !< Current geometric snow height [m].
  REAL(wp), INTENT(INOUT) :: dz_snow_flx(nvec) !< Geometric snow height for heat-flux computation [m].
  REAL(wp), INTENT(INOUT) :: fr_snow_lim(nvec) !< Snow fraction with freshsnow-dependent lower limit [m].
  REAL(wp), INTENT(INOUT) :: rho_snow(nvec) !< Effective snow density [kg/m^3].
  REAL(wp), INTENT(INOUT) :: hcap_snow(nvec) !< Effective snow heat capacity [J/(m^2 K)].
  REAL(wp), INTENT(INOUT) :: t_snow_top(nvec) !< Effective snow top temperature [K].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number.

  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (w_snow_now(i) > 0.0_wp) THEN
      ! existence of snow
      ! --> no water in interception store and t_snow < t0_melt
      ! GZ: this effectively suppresses rime formation because deposition rates per time
      ! step are usually less than 1e-6 m (eps_temp, eps_soil)
 !!!  w_snow_now(i) = w_snow_now(i) + w_i_now(i)
 !!!  w_i_now(i) = 0.0_wp
      t_snow_now(i) = MIN (t0_melt - eps_temp, t_snow_now(i) )
    ELSE IF (t_snow_now(i) >= t0_melt) THEN
      ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
      t_s_now   (i) = MAX (t0_melt + eps_temp, t_s_now(i) )
      t_snow_now(i) = t_s_now(i)
    ELSE
      ! no snow and  t_snow < t0_melt
      ! --> t_snow = t_s and no water w_i in interception store
      t_s_now   (i) = MIN (t0_melt - eps_temp, t_s_now(i) )
      t_snow_now(i) = t_s_now(i)
 !!!  w_i_now(i) = 0.0_wp
    END IF
  ENDDO

  ! set t_snow_top to t_snow
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    t_snow_top(i) = t_snow_now(i)
    h_snow_now(i) = w_snow_now(i)*rho_w/rho_snow_now(i)
    h_snow_new(i) = h_snow_now(i)

    ! Decide which snow density is used for computing the heat capacity
    IF (l2lay_rho_snow) THEN
      rho_snow(i) = rho_snow_mult_now_top(i)
    ELSE
      rho_snow(i) = rho_snow_now(i)
    ENDIF

    ! constrain snow depth and consider this constraint for the computation
    ! of average snow density of snow layer
    fr_snow_lim(i) = MAX(0.01_wp,0.1_wp*freshsnow(i),fr_snow(i))
    dz_snow_flx(i) = MIN(MAX(cdsmin,h_snow_now(i)/fr_snow_lim(i)), hmax_single_heattransfer)

    hcap_snow(i) = chc_i * dz_snow_flx(i) * rho_snow(i)
  ENDDO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_single_prepare

!>
!! Compute heat flux into the soil top from snow and atmospheric fluxes.
!!
SUBROUTINE snow_single_soil_forcing ( &
      & ivstart, ivend, nvec, dt, t_s_now, fr_snow, t_snow_top, w_snow_now, rho_snow, dz_snow_flx, &
      & dt_w_snow, radfl_net_snfr, shfl_snfr, lhfl_snfr, lhfl_precip, hfl_anthrop, hfl_snow_soil, &
      & forcing_soil, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec

  REAL(wp), INTENT(IN) :: dt !< Time step [s].

  REAL(wp), INTENT(IN) :: t_s_now(nvec) !< Surface temperature (current) [K].
  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< Snow-covered fraction [1].

  REAL(wp), INTENT(IN) :: t_snow_top(nvec) !< Snow surface temperature [K]
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current snow water equivalent [m(H2O)].
  REAL(wp), INTENT(IN) :: rho_snow(nvec) !< Snow density [kg/m^3].
  REAL(wp), INTENT(IN) :: dz_snow_flx(nvec) !< Snow height used for computing heat flux [m].
  REAL(wp), INTENT(IN) :: dt_w_snow(nvec) !< SWE rate of change [m(H2O)/s].

  REAL(wp), INTENT(IN) :: radfl_net_snfr(nvec) !< Net radiation at soil surface [W/m^2].
  REAL(wp), INTENT(IN) :: shfl_snfr(nvec) !< Sensible heat flux at soil surface [W/m^2].
  REAL(wp), INTENT(IN) :: lhfl_snfr(nvec) !< Latent heat flux at soil surface [W/m^2].
  REAL(wp), INTENT(IN) :: lhfl_precip(nvec) !< Precipitation latent heat flux [W/m^2].
  REAL(wp), INTENT(IN) :: hfl_anthrop(nvec) !< Anthropogenic heat flux [W/m^2].

  REAL(wp), INTENT(INOUT) :: hfl_snow_soil(nvec) !< Heat flux between snow and soil [W/m^2].
  REAL(wp), INTENT(INOUT) :: forcing_soil(nvec) !< Soil surface forcing [W/m^2].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  INTEGER :: i

  REAL(wp) :: w_snow_prov
  REAL(wp) :: hcond_snow

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG VECTOR PRIVATE(w_snow_prov, hcond_snow)
  !$NEC sparse
  DO i = ivstart, ivend
    ! Influence of heatflux through snow on total forcing:
    w_snow_prov = w_snow_now(i) + dt_w_snow(i)*dt / rho_w
    hfl_snow_soil(i) = 0._wp

    IF (w_snow_prov > eps_soil) THEN
      !         heat conductivity of snow as funtion of water content
      ! BR      hcond_snow  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*w_snow_now(i)))
      !
      ! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
      !
      hcond_snow  = 2.22_wp*EXP(1.88_wp*LOG(rho_snow(i)/rho_i))

      ! BR 11/2005 Use alternative formulation for heat conductivity by Sun et al., 1999
      !            The water vapour transport associated conductivity is not included.

      !        hcond_snow   = 0.023_wp+(2.290_wp-0.023_wp)* &
      !                               (7.750E-05_wp*rho_snow(i,nx) + &
      !                                1.105E-06_wp*prho_snow(i,nx)**2)

      !          hfl_snow_soil(i) = hcond_snow*(ztsnow(i) - t_s_now(i))/zdz_snow_fl(i)
      hfl_snow_soil(i) = hcond_snow*(t_snow_top(i) - t_s_now(i))/dz_snow_flx(i)
    END IF

    ! Calculation of the surface energy balance

    ! total forcing for uppermost soil layer
    forcing_soil(i) = ( radfl_net_snfr(i) + shfl_snfr(i) + lhfl_snfr(i) )       &
                    * (1._wp - fr_snow(i)) + lhfl_precip(i) + hfl_anthrop(i) &
                + MERGE(fr_snow(i) * hfl_snow_soil(i), 0._wp, t_snow_top(i) < t0_melt)

  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_single_soil_forcing

!>
!! Compute new snow temperature and atmospheric fluxes into the snow surface.
!!
SUBROUTINE snow_single_calc_temperature ( &
      & ivstart, ivend, nvec, dt, t_so_new_top, t_s_now, w_snow_now, dt_w_snow, dz_snow_flx, &
      & fr_snow, rho_snow, sobs, radfl_th_snow, rho_ch, hcap_snow, th_atm, evapo_snow, &
      & hfl_snow_soil, dqvdt_snow, t_snow_top, t_snow_new, dt_t_snow, shfl_snow, lhfl_snow, &
      & qhfl_snow, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec

  REAL(wp), INTENT(IN) :: dt !< Time step [s].

  REAL(wp), INTENT(IN) :: t_so_new_top(nvec) !< New top-layer soil temperature [K].
  REAL(wp), INTENT(IN) :: t_s_now(nvec) !< Current soil surface temperature [K].

  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current SWE [m(H2O)].
  REAL(wp), INTENT(IN) :: dt_w_snow(nvec) !< Snow mass tendency [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: dz_snow_flx(nvec) !< Effective snow height for flux computation [m].
  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< Snow fraction [m^2(snow)/m^2(tile)].

  REAL(wp), INTENT(IN) :: rho_snow(nvec) !< Effective snow density [kg/m^3].

  REAL(wp), INTENT(IN) :: sobs(nvec) !< Net short-wave radiation at surface [W/m^2].
  REAL(wp), INTENT(IN) :: radfl_th_snow(nvec) !< Net thermal radiation at snow surface [W/m^2].
  REAL(wp), INTENT(IN) :: rho_ch(nvec) !< Density * heat transfer velocity [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: hcap_snow(nvec) !< Heat capacity of snow [J/(m^2 K)].
  REAL(wp), INTENT(IN) :: th_atm(nvec) !< Potential temperature of lowest atmospheric layer [K].
  REAL(wp), INTENT(IN) :: evapo_snow(nvec) !< Evaporation rate from snow [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: hfl_snow_soil(nvec) !< Heat flux from snow to soil [W/m^2].
  REAL(wp), INTENT(IN) :: dqvdt_snow(nvec) !< Derivative of sat. spec. humidity w.r.t. snow temperature [(kg/kg)/K].

  REAL(wp), INTENT(IN) :: t_snow_top(nvec) !< Snow surface temperature [K].
  REAL(wp), INTENT(INOUT) :: t_snow_new(nvec) !< New snow temperature [K].
  REAL(wp), INTENT(INOUT) :: dt_t_snow(nvec) !< Tendency of snow temperature [K/s].

  REAL(wp), INTENT(INOUT) :: shfl_snow(nvec) !< Sensible heat flux into snow surface [W/m^2(snow)].
  REAL(wp), INTENT(INOUT) :: lhfl_snow(nvec) !< Latent heat flux into snow surface [W/m^2(snow)].
  REAL(wp), INTENT(INOUT) :: qhfl_snow(nvec) !< Humidity flux into snow surface [kg/(m^2(snow) s)].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp) :: ztsn
  REAL(wp) :: w_snow_prov
  REAL(wp) :: zrnet_snow
  REAL(wp) :: zfor_snow
  REAL(wp) :: hcond_snow
  REAL(wp) :: ztsnow_im
  REAL(wp) :: zfak

  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(ztsn, w_snow_prov, zrnet_snow, zfor_snow, hcond_snow, ztsnow_im, zfak)
  DO i = ivstart, ivend
    ! next line has to be changed if a soil surface temperature is
    ! predicted by the heat conduction equation
    ztsn = t_so_new_top(i)
    t_snow_new(i) = ztsn ! default setting
    w_snow_prov = w_snow_now(i) + dt_w_snow(i) * dt / rho_w

    ! forcing contributions for snow formation of dew and rime are
    ! contained in ze_ges, heat fluxes must not be multiplied by
    ! snow covered fraction

    zrnet_snow    = sobs(i) + radfl_th_snow(i)
    shfl_snow(i) = rho_ch(i)*cp_d*(th_atm(i) - t_snow_top(i))
    lhfl_snow(i) = lh_s*evapo_snow(i) / MAX(eps_div, fr_snow(i))
    qhfl_snow(i) = evapo_snow(i) / MAX(eps_div, fr_snow(i))
    zfor_snow     = zrnet_snow + shfl_snow(i) + lhfl_snow(i)

    ! forecast of snow temperature Tsnow
    IF (t_snow_top(i) < t0_melt .AND. w_snow_prov > eps_soil) THEN
      t_snow_new(i) = t_snow_top(i) + dt * 2._wp * (zfor_snow - hfl_snow_soil(i))  &
                      /hcap_snow(i) - ( ztsn - t_s_now(i) )

      ! implicit formulation
      ! BR        hcond_snow  = MAX(calasmin,MIN(calasmax, calasmin + calas_dw*w_snow_prov(i)))
      ! BR 7/2005 Introduce new dependency of snow heat conductivity on snow density
      !
      hcond_snow  = 2.22_wp*EXP(1.88_wp*LOG(rho_snow(i)/rho_i))

      ztsnow_im    = - rho_ch(i) * (cp_d + dqvdt_snow(i) * lh_s)       &
                                        - hcond_snow/dz_snow_flx(i)
      zfak  = MAX(eps_div,1.0_wp - dt*zalfa*ztsnow_im/hcap_snow(i))
      t_snow_new(i) = t_snow_top(i) + (t_snow_new(i)-t_snow_top(i))/zfak
    END IF

    dt_t_snow(i) = (t_snow_new(i) - t_snow_top(i)) / dt
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_single_calc_temperature

!>
!! Handle melting of snow when the snow-top or soil temperature go above freezing.
!!
!! Basically this snow model uses heat fluxes to either heat the uppermost soil
!! layer, if the snow surface temperature exceeds t0_melt, or to heat the snow, if
!! the temperature of the uppermost soil layer exceeds t0_melt. Melting is considered
!! after this process, if the snow temperature equals t0_melt AND the temperature of
!! the uppermost soil layer exceeds t0_melt. The excess heat (t_so(1)-t0_melt) is used for
!! melting snow. In cases this melting may be postponed to the next time step.
!!
SUBROUTINE snow_single_melt ( &
      & ivstart, ivend, nvec, dt, soiltyp_subs, dz_top, hcap_ml_top, hcap_snow, fr_snow_lim, &
      & w_snow_now, dt_w_snow, t_snow, dt_t_snow, rho_snow_now, t_s_now, t_so_top, dt_w_so_top, &
      & dt_t_s, fr_w_top, fr_ice_top, meltrate, runoff_s, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec

  REAL(wp), INTENT(IN) :: dt

  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type.

  REAL(wp), INTENT(IN) :: dz_top !< Top layer thickness [m].
  REAL(wp), INTENT(IN) :: hcap_ml_top(nvec) !< Soil volumentric heat capacity [J/(m^3 K)]
  REAL(wp), INTENT(IN) :: hcap_snow(nvec) !< Snow heat capacity [J/(m^2 K)].

  REAL(wp), INTENT(IN) :: fr_snow_lim(nvec) !< Snow fraction, limited [m^2/m^2].
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current SWE [m].
  REAL(wp), INTENT(INOUT) :: dt_w_snow(nvec) !< Snow mass tendency [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: t_snow(nvec) !< Snow surface temperature [K].
  REAL(wp), INTENT(INOUT) :: dt_t_snow(nvec) !< Snow surface temperature tendency [K].
  REAL(wp), INTENT(IN) :: rho_snow_now(nvec) !< Snow density [kg/m^3].
  REAL(wp), INTENT(IN) :: t_s_now(nvec) !< Soil top temperature [K].
  REAL(wp), INTENT(IN) :: t_so_top(nvec) !< Top soil-layer temperature [K].
  REAL(wp), INTENT(INOUT) :: dt_w_so_top(nvec) !< Top soil water tendency [kg/(m^2 s)].
  REAL(wp), INTENT(INOUT) :: dt_t_s(nvec) !< Top soil-layer temperature tendency [K/s].
  REAL(wp), INTENT(IN) :: fr_w_top(nvec) !< Total water fraction in top soil layer [m^3/m^3].
  REAL(wp), INTENT(IN) :: fr_ice_top(nvec) !< Ice fraction in top soil layer [m^3/m^3].
  REAL(wp), INTENT(INOUT) :: meltrate(nvec) !< Snow meltrate [kg/(m^2 s)].
  REAL(wp), INTENT(INOUT) :: runoff_s(nvec) !< Surface runoff (accumulated) [kg/m^2].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp) :: w_snow_prov !< Provisional updated SWE [m].

  REAL(wp) :: field_cap
  REAL(wp) :: pore_vol

  REAL(wp) :: zdwsnm
  REAL(wp) :: zro
  REAL(wp) :: ze_avail
  REAL(wp) :: ze_total
  REAL(wp) :: zfr_melt
  REAL(wp) :: zfr_ice_free
  REAL(wp) :: ztsnew
  REAL(wp) :: ztsnownew
  REAL(wp) :: zdwgme
  REAL(wp) :: zredfu
  REAL(wp) :: zw_ovpv
  REAL(wp) :: zdelt_s

  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  ! If the soil surface temperature predicted by the equation of heat conduction
  ! is used instead of using T_s = T_so(1), the following section has to be
  ! adjusted accordingly.


  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(w_snow_prov, field_cap, pore_vol, zdwsnm, zro, ze_avail, ze_total, zfr_melt) &
  !$ACC   PRIVATE(zfr_ice_free, ztsnew, ztsnownew, zdwgme, zredfu, zw_ovpv, zdelt_s)
  DO i = ivstart, ivend

    w_snow_prov   = w_snow_now(i) + dt_w_snow(i) * dt / rho_w
    ze_avail      = 0.0_wp
    ze_total      = 0.0_wp
    zfr_melt      = 0.0_wp

    field_cap = cfcap(soiltyp_subs(i))
    pore_vol = cporv(soiltyp_subs(i))

    IF (w_snow_prov > eps_soil) THEN        ! points with snow cover only
      ! first case: T_snow > t0_melt: melting from above
      ! ----------
      IF (t_snow(i) > t0_melt .AND. t_so_top(i) < t0_melt ) THEN
        ! Limit max w_snow in melting conditions for consistency with heat capacity calculation
        zdwsnm       = MIN(1.5_wp*rho_snow_now(i)/rho_w,w_snow_prov)* &
              .5_wp*(t_snow(i) - (t0_melt - eps_temp))/ &
             (.5_wp* (t_s_now(i) - (t0_melt - eps_temp)) - lh_f/chc_i)
        zdwsnm       = zdwsnm * rho_w/dt
        dt_w_snow(i) = dt_w_snow(i) + zdwsnm
        meltrate(i)  = - zdwsnm
        ztsnownew    = t0_melt - eps_temp
        dt_t_snow(i) = dt_t_snow(i) + (ztsnownew - t_snow(i))/dt

        ! decide which parts of the meltwater are passed to w_so and runoff, respectively
        zro = meltrate(i)
        IF (soiltyp_subs(i) > IST_ROCK) THEN
          zfr_ice_free = 1._wp-fr_ice_top(i)/pore_vol
          zdwgme       = zfr_ice_free*zro       ! contribution to w_so
          zredfu       = MAX( &
              & 0.0_wp, &
              & MIN(1.0_wp, &
              &   (fr_w_top(i) - field_cap)/MAX(pore_vol-field_cap,eps_div) &
              & ) &
            )
          dt_w_so_top(i) = dt_w_so_top(i) + zdwgme*(1._wp - zredfu)
          zro = zro - zdwgme*(1._wp - zredfu)
        END IF

        ! zro-, zdw_so_dt-correction in case of pore volume overshooting
        zw_ovpv = MAX(0._wp, (fr_w_top(i)-pore_vol) * dz_top * rho_w / dt + dt_w_so_top(i) )
        zro = zro + zw_ovpv
        dt_w_so_top(i) = dt_w_so_top(i) - zw_ovpv
        runoff_s(i) = runoff_s(i) + zro * dt
      ENDIF ! melting from above

      IF (t_so_top(i) >= t0_melt) THEN
        !second case:  temperature of uppermost soil layer > t0_melt. First a
        !-----------   heat redistribution is performed. As a second step,
        !              melting of snow is considered.
        ! a) Heat redistribution
        ztsnew = t0_melt + eps_temp
        ztsnownew      = t_snow(i) + fr_snow_lim(i)*(t_so_top(i) - ztsnew) +  &
             2._wp*(t_so_top(i) - ztsnew)*hcap_ml_top(i)*dz_top/hcap_snow(i)
        dt_t_s(i)    = dt_t_s(i) + fr_snow_lim(i)*(ztsnew - t_so_top(i)) / dt
        dt_t_snow(i) = dt_t_snow(i) + (ztsnownew - t_snow(i)) / dt
        ! b) Melting of snow (if possible)
        IF (ztsnownew > t0_melt) THEN
          ze_avail     = 0.5_wp*(ztsnownew - t0_melt)*hcap_snow(i)*fr_snow_lim(i)
          ze_total     = lh_f*w_snow_prov*rho_w
          zfr_melt     = MIN(1.0_wp,ze_avail/ze_total)
          dt_t_snow(i) = dt_t_snow(i) + (t0_melt - ztsnownew) / dt
          zdelt_s      = MAX(0.0_wp,(ze_avail - ze_total)/(hcap_ml_top(i)*dz_top))
          dt_t_s(i)    = dt_t_s(i) + fr_snow_lim(i)*zdelt_s / dt

          IF (zfr_melt > 0.9999_wp) zfr_melt = 1._wp

          ! melted snow is allowed to penetrate the soil (up to field
          ! capacity), if the soil type is neither ice nor rock;
          ! else it contributes to surface run-off;
          ! fractional water content of the first soil layer determines
          ! a reduction factor which controls additional run-off
          zdwsnm      = zfr_melt * w_snow_prov * rho_w / dt  ! available water
          dt_w_snow(i) = dt_w_snow(i) - zdwsnm
          meltrate(i) = meltrate(i) + zdwsnm

          IF (soiltyp_subs(i) > IST_ROCK) THEN
            zredfu = MAX( &
                & 0.0_wp, &
                & MIN(1.0_wp, &
                &   (fr_w_top(i) - field_cap)/MAX(pore_vol-field_cap, eps_div) &
                & ) &
              )
            dt_w_so_top(i) = dt_w_so_top(i) + zdwsnm*(1._wp - zredfu)
            zro = zdwsnm*zredfu    ! Infiltration not possible
                                   ! for this fraction
          ELSE
            zro = zdwsnm      ! surface runoff
          END IF

          ! zro-, zdw_so_dt-correction in case of pore volume overshooting
          zw_ovpv = MAX(0._wp, (fr_w_top(i) - pore_vol) * dz_top * rho_w / dt + &
                     dt_w_so_top(i))
          zro = zro + zw_ovpv
          dt_w_so_top(i) = dt_w_so_top(i) - zw_ovpv

          runoff_s(i) = runoff_s(i) + zro * dt

        END IF   ! snow melting
      END IF     ! snow and/or soil temperatures
    END IF       ! points with snow cover only
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_single_melt


!>
!! Final update of new snow variables.
!!
!! Updates snow temperature, density and height as well as melting a fraction of falling
!! snow when air temperature is above 0.5 C.
!!
SUBROUTINE snow_single_update_new_state ( &
      & ivstart, ivend, nvec, ke_snow, dt, dz_hl_top, t_snow_new, t_snow_now, dt_t_snow, &
      & w_snow_new, w_snow_now, dt_w_snow, rho_snow_new, rho_snow_now, rho_snow_mult_new, &
      & rho_snow_mult_now, h_snow_new, h_snow_gp, t_so_new_top, w_i_new, w_i_max, sp_10m, &
      & th_atm, fr_w_top, dt_w_so_top, runoff_s, soiltyp_subs, snow_rate, ice_rate, &
      & graupel_rate, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: ke_snow

  REAL(wp), INTENT(IN) :: dt !< Time step [s].
  REAL(wp), INTENT(IN) :: dz_hl_top !< Top soil layer thickness [m]

  REAL(wp), INTENT(INOUT) :: t_snow_new(nvec) !< New snow temperature [K].
  REAL(wp), INTENT(IN) :: t_snow_now(nvec) !< Current snow temperature [K].
  REAL(wp), INTENT(IN) :: dt_t_snow(nvec) !< Tendency of snow temperature [K/s].

  REAL(wp), INTENT(INOUT) :: w_snow_new(nvec) !< New SWE [m(H2O)].
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current SWE [m(H2O)].
  REAL(wp), INTENT(IN) :: dt_w_snow(nvec) !< Snow mass tendency [kg/(m^2 s)].

  REAL(wp), INTENT(INOUT) :: rho_snow_new(nvec) !< New snow density [kg/m^3].
  REAL(wp), INTENT(IN) :: rho_snow_now(nvec) !< Current snow density [kg/m^3].

  REAL(wp), INTENT(INOUT) :: rho_snow_mult_new(nvec, ke_snow) !< New snow density (2-layer density) [kg/m^3].
  REAL(wp), INTENT(IN) :: rho_snow_mult_now(nvec, ke_snow) !< Current snow density (2-layer density) [kg/m^3].

  REAL(wp), INTENT(INOUT) :: h_snow_new(nvec) !< New geometric snow height [m].
  REAL(wp), INTENT(IN) :: h_snow_gp(nvec) !< Average geometric snow height on grid point [m].

  REAL(wp), INTENT(IN) :: t_so_new_top(nvec) !< Top soil layer temperature [K].

  REAL(wp), INTENT(INOUT) :: w_i_new(nvec) !< Content of interception reservoir [m(H2O)].
  REAL(wp), INTENT(IN) :: w_i_max(nvec) !< Interception reservoir capacity [m(H2O)].

  REAL(wp), INTENT(IN) :: sp_10m(nvec) !< 1m wind speed [m/s].
  REAL(wp), INTENT(IN) :: th_atm(nvec) !< Potential temperature of lowest atmospheric level [K].

  REAL(wp), INTENT(IN) :: fr_w_top(nvec) !< Fractional water content (top soil layer) [m^3(H2O)/m^3(pore)].

  REAL(wp), INTENT(INOUT) :: dt_w_so_top(nvec) !< Tendency of top-layer soil water [kg/(m^2 s)].

  REAL(wp), INTENT(INOUT) :: runoff_s(nvec) !< Surface runoff [kg/m^2].

  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type.
  REAL(wp), INTENT(IN) :: snow_rate(nvec) !< Total snow rate excluding ice [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: ice_rate(nvec) !< Ice precipitation rate [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: graupel_rate(nvec) !< Graupel rate [kg/(m^2 s)].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp) :: graupel_frac
  REAL(wp) :: runoff_frac

  REAL(wp) :: rho_graupel
  REAL(wp) :: rho_snow_existing
  REAL(wp) :: rho_snow_fresh
  REAL(wp) :: rho_snow_min
  REAL(wp) :: rho_snow_max

  REAL(wp) :: t_snow_rel
  REAL(wp) :: th_low_rel

  REAL(wp) :: tau_snow_days

  REAL(wp) :: total_snow_rate
  REAL(wp) :: total_snow_m
  REAL(wp) :: infiltration_m
  REAL(wp) :: runoff_rate
  REAL(wp) :: snow_acc_m
  REAL(wp) :: zw_ovpv

  REAL(wp) :: w_snow_gp
  REAL(wp) :: dw_g_melt

  INTEGER :: mstyp

  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(mstyp, t_snow_rel, tau_snow_days, rho_snow_max) &
  !$ACC   PRIVATE(rho_snow_existing, th_low_rel, rho_snow_min, rho_snow_fresh, rho_graupel) &
  !$ACC   PRIVATE(total_snow_rate, graupel_frac, total_snow_m, snow_acc_m, infiltration_m) &
  !$ACC   PRIVATE(dw_g_melt, runoff_rate, runoff_frac, zw_ovpv, w_snow_gp)
  !$NEC sparse
  DO i = ivstart, ivend
    w_snow_new(i) = w_snow_now(i) + dt*dt_w_snow(i)/rho_w

    t_snow_new(i) = t_snow_now(i) + dt*dt_t_snow(i)

    ! Reset t_snow_new to t_so(0) if no snow was present at the beginning of the time step
    ! The heat balance calculation is incomplete in this case and sometimes yields unreasonable results
    IF (w_snow_now(i) < eps_soil .AND. w_snow_new(i) >= eps_soil) THEN
      t_snow_new(i) = MIN(t0_melt,t_so_new_top(i))
    ENDIF

    ! Put small snow amounts to interception.
    IF (w_snow_new(i) <= eps_soil) THEN
      w_i_new(i)    = w_i_new(i) + w_snow_new(i)
      w_snow_new(i) = 0.0_wp
      t_snow_new(i) = t_so_new_top(i)
    ENDIF

    IF (w_i_new(i) <= 1.0E-4_wp*eps_soil) w_i_new(i) = 0.0_wp

    mstyp = soiltyp_subs(i)

    !BR 7/2005 Update snow density
    ! a) aging of existing snow
    !    temperature dependence of relaxation/ageing constant
    t_snow_rel = (t_snow_new(i) - csnow_tmin) / (t0_melt - csnow_tmin)
    tau_snow_days = MIN(MAX( &
        & 0.05_wp, &
        & crhosmint + (crhosmaxt - crhosmint) * t_snow_rel), &
        & crhosmaxt &
      )

    ! use 20 days in combination with temperature-dependent equilibrium density
    rho_snow_max = crhosmax_tmin + MAX(-0.25_wp, t_snow_rel) * (crhosmax_ml - crhosmax_tmin)
    rho_snow_existing = MAX( &
        & rho_snow_now(i), &
        & rho_snow_max+(rho_snow_now(i)-rho_snow_max) * EXP(-tau_snow_days*dt/86400._wp) &
      )

    ! b) density of fresh snow; the minimum density depends on wind speed to account for wind compression of powder snow
    th_low_rel = (th_atm(i) - csnow_tmin) / (t0_melt - csnow_tmin)
    rho_snow_min = crhosminf * MIN(MAX( &
        & 1._wp, &
        & 0.125_wp * sp_10m(i) + 0.5_wp), &
        & 2._wp &
      )
    rho_snow_fresh = MIN(MAX( &
        & rho_snow_min, &
        & rho_snow_min + (crhosmaxf - rho_snow_min) * (th_low_rel)**2), &
        & crhosmaxf &
      )
    rho_graupel = MIN(MAX( &
        & crhogminf, &
        & crhogminf + (crhogmaxf - crhogminf) * (th_low_rel)**2), &
        & crhogmaxf &
      )

    ! graupel fraction
    total_snow_rate = snow_rate(i) + ice_rate(i)
    graupel_frac = graupel_rate(i) / MAX(eps_soil, total_snow_rate)
    total_snow_m = total_snow_rate * (dt / rho_w)

    ! prevent accumulation of new snow if the air temperature is above 1 deg C with
    ! linear transition between 0.5 and 1 deg C
    IF (total_snow_m > 0.5_wp*eps_soil .AND. th_atm(i) > t0_melt + 0.5_wp) THEN
      ! Some of the snow may have melted or sublimated already.
      total_snow_m = MIN(w_snow_new(i), total_snow_m)

      ! part of the new snow that accumulates on the ground
      snow_acc_m = MAX(0._wp, total_snow_m*(t0_melt + 1._wp - th_atm(i))*2._wp)
      !
      ! the rest is transferred to the interception storage, soil moisture or runoff:
      w_i_new(i) = w_i_new(i) + total_snow_m - snow_acc_m
      IF (w_i_new(i) > w_i_max(i)) THEN  ! overflow of interception store
        infiltration_m = w_i_new(i) - w_i_max(i)
        w_i_new(i) = w_i_max(i)
      ELSE
        infiltration_m = 0.0_wp
      ENDIF

      IF (soiltyp_subs(i) > IST_ROCK) THEN
        dw_g_melt = infiltration_m * (rho_w / dt)
        runoff_rate = 0._wp
      ELSE
        dw_g_melt = 0._wp
        runoff_rate = infiltration_m * (rho_w / dt)
      END IF

      runoff_frac = MIN(MAX( &
          & 0._wp, &
          & (fr_w_top(i) - cfcap(mstyp)) / MAX(cporv(mstyp) - cfcap(mstyp), eps_div)), &
          & 1._wp &
        )
      dt_w_so_top(i) = dt_w_so_top(i) + dw_g_melt*(1._wp - runoff_frac)
      runoff_rate = runoff_rate + dw_g_melt * runoff_frac

      ! runoff_rate-, zdw_so_dt-correction in case of pore volume overshooting
      zw_ovpv = MAX(0._wp, &
          & (fr_w_top(i) - cporv(mstyp)) * dz_hl_top * (rho_w / dt) + dt_w_so_top(i) &
        )
      runoff_rate = runoff_rate + zw_ovpv
      dt_w_so_top(i) = dt_w_so_top(i) - zw_ovpv

      runoff_s(i) = runoff_s(i) + runoff_rate * dt

      ! correct SWE for immediately melted new snow
      w_snow_new(i) = w_snow_new(i) - (total_snow_m - snow_acc_m)
    ELSE
      snow_acc_m = total_snow_m
    END IF

    ! c) new snow density is computed by adding depths of existing and new snow
    IF (mstyp /= IST_ICE) THEN
      rho_snow_new(i) = (w_snow_now(i) + snow_acc_m) / &
          & ( MAX(w_snow_now(i), eps_soil) / rho_snow_existing &
          &   + snow_acc_m / ( (1.0_wp-graupel_frac)*rho_snow_fresh + graupel_frac*rho_graupel) &
          & )
    ELSE ! constant snow density over glaciers
      rho_snow_new(i) = rho_snow_now(i)
    ENDIF
    rho_snow_new(i) = MIN(MAX(crhosmin_ml, rho_snow_new(i)), crhosmax_ml)

    ! New calculation of snow height for single layer snow model
    h_snow_new(i) = w_snow_new(i) * rho_w / rho_snow_new(i)

    ! Calculation of top-layer snow density for two-layer snow density scheme
    IF (l2lay_rho_snow) THEN
      rho_snow_existing = MAX(rho_snow_mult_now(i,1), &
          & rho_snow_max + (rho_snow_mult_now(i,1) - rho_snow_max) * EXP(-tau_snow_days * dt/86400._wp))
      w_snow_gp = MIN(max_toplaydepth, h_snow_gp(i)) * rho_snow_mult_now(i,1) / rho_w
      rho_snow_mult_new(i,1) = (w_snow_gp + snow_acc_m) / &
          & (MAX(w_snow_gp, eps_div) / rho_snow_existing + snow_acc_m / rho_snow_fresh)
      rho_snow_mult_new(i,1) = MIN(MAX(crhosmin_ml, rho_snow_mult_new(i,1)), crhosmax_ml)
      rho_snow_mult_new(i,2) = rho_snow_new(i)
    ENDIF

  END DO
  !$ACC END PARALLEL

  !$ACC END DATA

END SUBROUTINE snow_single_update_new_state


!--------------------------------------------------------------------------------------------------
! Multi-Layer Snow Model
!--------------------------------------------------------------------------------------------------


SUBROUTINE snow_multi_prepare ( &
      & ivstart, ivend, nvec, ke_snow, t_s_now, w_i_now, w_snow_now, h_snow_new, h_snow_now, &
      & h_snow, t_snow_mult_now, wtot_snow_now, dzh_snow_now, rho_snow_mult_now, zhh_snow, &
      & zhm_snow, zdzh_snow, zextinct, t_snow_top, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: ke_snow

  REAL(wp), INTENT(INOUT) :: t_s_now(nvec)
  REAL(wp), INTENT(INOUT) :: w_i_now(nvec)

  REAL(wp), INTENT(INOUT) :: w_snow_now(nvec)
  REAL(wp), INTENT(INOUT) :: h_snow_new(nvec)
  REAL(wp), INTENT(INOUT) :: h_snow_now(nvec)
  REAL(wp), INTENT(IN) :: h_snow(nvec)

  REAL(wp), INTENT(INOUT) :: t_snow_mult_now(nvec,0:ke_snow)
  REAL(wp), INTENT(INOUT) :: wtot_snow_now(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: dzh_snow_now(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: rho_snow_mult_now(nvec,ke_snow)

  REAL(wp), INTENT(INOUT) :: zhh_snow(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: zhm_snow(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: zdzh_snow(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: zextinct(nvec,ke_snow)

  REAL(wp), INTENT(INOUT) :: t_snow_top(nvec)

  LOGICAL, INTENT(IN) :: lzacc
  INTEGER, INTENT(IN) :: acc_async_queue

  INTEGER :: ksn
  INTEGER :: i

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP SEQ
  DO ksn = 0,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_now(i) > 0.0_wp) THEN
        ! existence of snow
        t_snow_mult_now(i,ksn) = MIN (t0_melt - eps_temp, t_snow_mult_now(i,ksn) )
      ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
        ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
        t_snow_mult_now(i,ksn) = MAX (t0_melt + eps_temp, t_s_now(i) )
      ELSE
        ! no snow and  t_snow < t0_melt
        ! --> t_snow = t_s
        t_snow_mult_now(i,ksn) = MIN (t0_melt - eps_temp, t_s_now(i) )
      END IF
    ENDDO
  ENDDO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (w_snow_now(i) > 0.0_wp) THEN
      ! existence of snow
      ! --> no water in interception store and t_snow < t0_melt
      w_snow_now(i) = w_snow_now(i) + w_i_now(i)
      wtot_snow_now(i,1) = wtot_snow_now(i,1) + w_i_now(i)
      dzh_snow_now(i,1)  = dzh_snow_now(i,1)  + w_i_now(i) / rho_snow_mult_now(i,1)*rho_w
      w_i_now(i)         = 0.0_wp
      h_snow_now(i)      = h_snow(i)
    ELSE IF (t_snow_mult_now(i,ke_snow) >= t0_melt) THEN
      ! no snow and t_snow >= t0_melt --> t_s > t0_melt and t_snow = t_s
      t_s_now   (i) = t_snow_mult_now(i,ke_snow)
      h_snow_now(i) = 0.0_wp
    ELSE
      ! no snow and  t_snow < t0_melt
      ! --> t_snow = t_s and no water w_i in interception store
      t_s_now   (i) = t_snow_mult_now(i,ke_snow)
      w_i_now(i)    = 0.0_wp
      h_snow_now(i) = 0.0_wp
    END IF
  ENDDO

  ! some preparations for ksn==0 and ksn==1
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    zhh_snow(i,1) =  -h_snow_now(i) + dzh_snow_now(i,1)
    zhm_snow(i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._wp

    zdzh_snow(i,1) = dzh_snow_now(i,1)
    zextinct (i,1) = 0.13_wp*rho_snow_mult_now(i,1)+3.4_wp

    ! set t_snow_top to t_snow_mult_now(ksn=1)
    t_snow_top   (i) = t_snow_mult_now(i,1)

    ! reinitialize and recompute h_snow_now
    h_snow_now(i) = dzh_snow_now(i,1) ! zdzh_snow
  ENDDO

  !$ACC LOOP SEQ
  DO  ksn = 2,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      zhh_snow   (i,ksn) = zhh_snow(i,ksn-1) + dzh_snow_now(i,ksn)
      zhm_snow   (i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._wp

      zdzh_snow  (i,ksn) = dzh_snow_now(i,ksn)
      zextinct   (i,ksn) = 0.13_wp*rho_snow_mult_now(i,ksn)+3.4_wp

      ! build sum over all layers for zdzh_snow in h_snow_now
      h_snow_now (i)     = h_snow_now(i) + zdzh_snow(i,ksn)
    ENDDO
  ENDDO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    h_snow_new(i) = h_snow_now(i)
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_multi_prepare


!> Compute changes of snow layer heights, densities, temperature, and water content due to snowfall
!! and rain.
SUBROUTINE snow_multi_handle_snowfall ( &
      & ivstart, ivend, nvec, ke_snow, dt, t_s_now, th_atm, rain_dew_rate, snow_rime_rate, w_snow_now, dt_w_snow, &
      & h_snow_now, wtot_snow_now, wliq_snow_now, rho_snow_mult_now, t_snow_mult_now, zhm_snow, &
      & zhh_snow, zdzh_snow, zdzm_snow, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: ke_snow

  REAL(wp), INTENT(IN) :: dt

  REAL(wp), INTENT(IN) :: t_s_now(nvec)
  REAL(wp), INTENT(IN) :: th_atm(nvec)

  REAL(wp), INTENT(IN) :: rain_dew_rate(nvec)
  REAL(wp), INTENT(IN) :: snow_rime_rate(nvec)

  REAL(wp), INTENT(IN) :: w_snow_now(nvec)
  REAL(wp), INTENT(IN) :: dt_w_snow(nvec)

  REAL(wp), INTENT(INOUT) :: h_snow_now(nvec)

  REAL(wp), INTENT(INOUT) :: wtot_snow_now(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: wliq_snow_now(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: rho_snow_mult_now(nvec,ke_snow)

  REAL(wp), INTENT(INOUT) :: t_snow_mult_now(nvec, 0:ke_snow)
  REAL(wp), INTENT(INOUT) :: zhm_snow(nvec, ke_snow)
  REAL(wp), INTENT(INOUT) :: zhh_snow(nvec, ke_snow)
  REAL(wp), INTENT(INOUT) :: zdzh_snow(nvec, ke_snow)
  REAL(wp), INTENT(INOUT) :: zdzm_snow(nvec, ke_snow)

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp) :: w_snow_prov(nvec)
  REAL(wp) :: sum_weight(nvec)
  REAL(wp) :: dz_old(nvec, ke_snow)
  REAL(wp) :: z_old(nvec, ke_snow)
  REAL(wp) :: t_new(nvec, ke_snow)
  REAL(wp) :: rho_new(nvec, ke_snow)
  REAL(wp) :: wl_new(nvec, ke_snow)

  INTEGER :: i
  INTEGER :: k
  INTEGER :: ksn

  REAL(wp) :: zrho_snowf
  REAL(wp) :: weight
  REAL(wp) :: dt_o_rho_w

  OPENACC_SUPPRESS_UNUSED_LZACC

  dt_o_rho_w = dt / rho_w

  !$ACC DATA CREATE(w_snow_prov, sum_weight, dz_old, z_old, t_new, rho_new, wl_new) &
  !$ACC   PRESENT(ivend) ASYNC(acc_async_queue)
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrho_snowf)
  DO i = ivstart, ivend

    w_snow_prov(i) = w_snow_now(i) + dt_w_snow(i)*dt_o_rho_w

    IF (w_snow_prov(i).GT.eps_soil) THEN

      zrho_snowf = crhosminf+(crhosmaxf-crhosminf)* (th_atm(i)-csnow_tmin) &
                                              /(t0_melt          -csnow_tmin)
      zrho_snowf = MAX(crhosminf,MIN(crhosmaxf,zrho_snowf))

      IF(dt_w_snow(i)-snow_rime_rate(i)-rain_dew_rate(i).GT.0.0_wp) THEN

        wtot_snow_now(i,1) = MAX(wtot_snow_now(i,1) + dt_w_snow(i)*dt_o_rho_w, 0.0_wp)

        zhm_snow(i,1) = zhm_snow(i,1) - (dt_w_snow(i)-snow_rime_rate(i)-                             &
        rain_dew_rate(i))*dt/rho_i/2._wp- snow_rime_rate(i)*dt/zrho_snowf/2._wp- rain_dew_rate(i)*dt/rho_i/2._wp
        zdzh_snow(i,1) = zdzh_snow(i,1) + (dt_w_snow(i)-snow_rime_rate(i)-rain_dew_rate(i))*dt/rho_i +        &
        snow_rime_rate(i)*dt/zrho_snowf + rain_dew_rate(i)*dt/rho_i

        rho_snow_mult_now(i,1) = MAX(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), 0.0_wp)
      ELSE

        wtot_snow_now(i,1) = MAX(wtot_snow_now(i,1) + (snow_rime_rate(i)+rain_dew_rate(i))*dt_o_rho_w, 0.0_wp)

        zhm_snow(i,1)  = zhm_snow(i,1) - snow_rime_rate(i)*dt/zrho_snowf/2._wp- &
                           rain_dew_rate(i)*dt/rho_i/2._wp
        zdzh_snow(i,1) = zdzh_snow(i,1) + snow_rime_rate(i)*dt/zrho_snowf + rain_dew_rate(i)*dt/rho_i

        IF (wtot_snow_now(i,1) .GT. 0._wp) THEN
          rho_snow_mult_now(i,1) = MAX(wtot_snow_now(i,1)*rho_w/zdzh_snow(i,1), 0.0_wp)

          wtot_snow_now(i,1) = MAX(wtot_snow_now(i,1) &
                                + (dt_w_snow(i)-snow_rime_rate(i)-rain_dew_rate(i))*dt_o_rho_w,0.0_wp)

          zhm_snow(i,1)  = zhm_snow(i,1) - (dt_w_snow(i)-snow_rime_rate(i)-rain_dew_rate(i)) &
                               *dt/rho_snow_mult_now(i,1)/2._wp
          zdzh_snow(i,1) = zdzh_snow(i,1) + (dt_w_snow(i)-snow_rime_rate(i)-rain_dew_rate(i)) &
                               *dt/rho_snow_mult_now(i,1)
        ELSE
          rho_snow_mult_now(i,1) = 0.0_wp
          zdzh_snow(i,1) = 0.0_wp
        END IF

      END IF
    END IF
    h_snow_now(i) = 0.0_wp
    sum_weight(i) = 0.0_wp
  END DO

  !$ACC LOOP SEQ
  DO ksn = 1,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_prov(i).GT.eps_soil) THEN
        h_snow_now(i) = h_snow_now(i) + zdzh_snow(i,ksn)
      END IF
    END DO
  END DO

  k = MIN(2,ke_snow-1)
  !$ACC LOOP SEQ
  DO ksn = 1,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_prov(i) .GT. eps_soil) THEN
        IF (w_snow_now(i) .GT. eps_soil) THEN
          IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
            zhh_snow(i,ksn) = -MAX( h_snow_now(i)-max_toplaydepth, h_snow_now(i)/ke_snow*(ke_snow-ksn) )
          ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
            zhh_snow(i,ksn) = MIN( 8._wp*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
          ELSE ! distribute the remaining snow equally among the layers
            zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
          ENDIF
        ELSE ! a newly generated snow cover will not exceed max_toplaydepth
          zhh_snow(i,ksn) = -h_snow_now(i)/ke_snow*(ke_snow-ksn)
        END IF
      END IF
    END DO
  END DO


  !$ACC LOOP SEQ
  DO ksn = ke_snow,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(w_snow_prov(i) .GT. eps_soil .AND. w_snow_now(i) .GT. eps_soil) THEN
        dz_old(i,ksn) = zdzh_snow(i,ksn)
        z_old(i,ksn) = -sum_weight(i) - zdzh_snow(i,ksn)/2._wp
        sum_weight(i) = sum_weight(i) + zdzh_snow(i,ksn)
      END IF
    END DO
  END DO


  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF(w_snow_prov(i) .GT. eps_soil) THEN
      zhm_snow (i,1) = (-h_snow_now(i) + zhh_snow(i,1))/2._wp
      zdzh_snow(i,1) = zhh_snow(i,1) + h_snow_now(i)   !layer thickness betw. half levels of uppermost snow layer
      zdzm_snow(i,1) = zhm_snow(i,1) + h_snow_now(i)   !layer thickness between snow surface and
                                                       !  main level of uppermost layer
      IF(w_snow_now(i) .GT. eps_soil) THEN
        IF(dz_old(i,1).ne.0..and.rho_snow_mult_now(i,1).ne.0.) THEN
          wliq_snow_now(i,1) = wliq_snow_now(i,1)/dz_old(i,1)
        END IF
      END IF
    END IF
  END DO


  !$ACC LOOP SEQ
  DO ksn = 2,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(w_snow_prov(i) .GT. eps_soil) THEN
        zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._wp
        zdzh_snow(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
        zdzm_snow(i,ksn) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
        IF(w_snow_now(i) .GT. eps_soil) THEN
          IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_now(i,ksn).ne.0.) THEN
            wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn)/dz_old(i,ksn)
          END IF
        END IF
      END IF
    END DO
  END DO

  !$ACC LOOP SEQ
  DO ksn = ke_snow,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      t_new  (i,ksn) = 0.0_wp
      rho_new(i,ksn) = 0.0_wp
      wl_new (i,ksn) = 0.0_wp
    END DO

    !$ACC LOOP SEQ
    DO k = ke_snow,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(weight)
      DO i = ivstart, ivend
        IF(w_snow_prov(i) .GT. eps_soil .AND. w_snow_now(i) .GT. eps_soil) THEN

          weight = MIN(dz_old(i, k), &
               &       z_old(i, k) + dz_old(i, k)*0.5_wp &
               &         - zhm_snow(i, ksn) + zdzh_snow(i, ksn)*0.5_wp,&
               &       zhm_snow(i, ksn) + zdzh_snow(i,ksn)*0.5_wp &
               &         - z_old(i, k) + dz_old(i, k)*0.5_wp, &
               &       zdzh_snow(i,ksn))

          weight = (weight + ABS(weight)) * 0.5_wp
          weight = weight / zdzh_snow(i,ksn)
          t_new  (i,ksn) = t_new  (i,ksn) + t_snow_mult_now(i,k)*weight
          rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_now(i,k)*weight
          wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_now(i,k)*weight
        END IF
      END DO
    END DO
  END DO

  !$ACC LOOP SEQ
  DO ksn = ke_snow,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(w_snow_prov(i) .GT. eps_soil) THEN
        IF(w_snow_now(i) .GT. eps_soil) THEN
          t_snow_mult_now  (i,ksn      ) = t_new  (i,ksn)
          rho_snow_mult_now(i,ksn) = rho_new(i,ksn)
          wtot_snow_now    (i,ksn) = rho_new(i,ksn)*zdzh_snow(i,ksn)/rho_w
          wliq_snow_now    (i,ksn) = wl_new (i,ksn)*zdzh_snow(i,ksn)
        ELSE ! Remark: if there was now snow in the previous time step, snow depth
             ! will not exceed the limit for equipartitioning
          t_snow_mult_now  (i,ksn      ) = t_s_now(i)
          rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,1)
          wtot_snow_now    (i,ksn) = w_snow_prov(i)/ke_snow
        END IF
      END IF
    END DO
  END DO
  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_multi_handle_snowfall


!> Compute soil and snow forcing as well as snow heat conductivity.
SUBROUTINE snow_multi_soil_forcing ( &
      & ivstart, ivend, nvec, ke_snow, ke_soil, dt, t_so_now, th_atm, fr_snow, hcond_ml, &
      & dz_ml, rho_ch, sobs, radfl_th_snow, radfl_net_snfr, shfl_snfr, lhfl_snfr, lhfl_precip, &
      & rain_dew_rate, evapo_snow, w_snow_now, dt_w_snow, rho_snow_mult_now, zalas_mult, &
      & t_snow_mult_now, zhm_snow, zextinct, hfl_snow_soil, forcing_soil, shfl_snow, lhfl_snow, &
      & qhfl_snow, zfor_snow_mult, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: ke_snow
  INTEGER, INTENT(IN) :: ke_soil

  REAL(wp), INTENT(IN) :: dt !< Time step [s].

  REAL(wp), INTENT(IN) :: t_so_now(nvec, 0:ke_soil+1) !< Soil temperature (current) [K].
  REAL(wp), INTENT(IN) :: th_atm(nvec) !< Atmospheric potential temperature [K].
  REAL(wp), INTENT(IN) :: fr_snow(nvec) !< Snow-covered fraction [1].
  REAL(wp), INTENT(IN) :: hcond_ml(nvec, ke_soil+1) !< Heat conductivity in soil layers [W/(K m)].

  REAL(wp), INTENT(IN) :: dz_ml(ke_soil+1) !< Layer main level distance (n_soil) [m].

  REAL(wp), INTENT(IN) :: rho_ch(nvec) !< Surface air density times transfer velocity [kg(air)/(m^2 s)].
  REAL(wp), INTENT(IN) :: sobs(nvec) !< Net shortwave radiation at surface [W/m^2].
  REAL(wp), INTENT(IN) :: radfl_th_snow(nvec) !< Net longwave radiation at snow surface [W/m^2].
  REAL(wp), INTENT(IN) :: radfl_net_snfr(nvec) !< Net radiation at soil surface [W/m^2].
  REAL(wp), INTENT(IN) :: shfl_snfr(nvec) !< Sensible heat flux at soil surface [W/m^2].
  REAL(wp), INTENT(IN) :: lhfl_snfr(nvec) !< Latent heat flux at soil surface [W/m^2].
  REAL(wp), INTENT(IN) :: lhfl_precip(nvec) !< Precipitation latent heat flux [W/m^2].
  REAL(wp), INTENT(IN) :: rain_dew_rate(nvec) !< Rain and dew formation rate [kg/(m^2 s)].
  REAL(wp), INTENT(IN) :: evapo_snow(nvec) !< Total evaporation from snow surface [kg/(m^2 s)].

  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current snow water equivalent [m(H2O)].
  REAL(wp), INTENT(IN) :: dt_w_snow(nvec) !< SWE rate of change [m(H2O)/s].

  REAL(wp), INTENT(IN) :: rho_snow_mult_now(nvec,ke_snow) !< Snow density in layer [kg/m^3].
  REAL(wp), INTENT(INOUT) :: zalas_mult(nvec, ke_snow) !< Heat conductivity in snow layers [W/(K m)].

  REAL(wp), INTENT(IN) :: t_snow_mult_now(nvec, 0:ke_snow) !< Temperature in snow layers [K].
  REAL(wp), INTENT(IN) :: zhm_snow(nvec, ke_snow) !< Height of snow ML [m].
  REAL(wp), INTENT(IN) :: zextinct(nvec, ke_snow) !< Extinction coefficient in snow layers [1].

  REAL(wp), INTENT(INOUT) :: hfl_snow_soil(nvec) !< Heat flux between snow and soil [W/m^2].
  REAL(wp), INTENT(INOUT) :: forcing_soil(nvec) !< Soil surface forcing [W/m^2].

  REAL(wp), INTENT(INOUT) :: shfl_snow(nvec) !< Sensible heat flux at snow surface [W/m^2].
  REAL(wp), INTENT(INOUT) :: lhfl_snow(nvec) !< Latent heat flux at snow surface [W/m^2].
  REAL(wp), INTENT(INOUT) :: qhfl_snow(nvec) !< Water vapor flux at snow surface [W/m^2].
  REAL(wp), INTENT(INOUT) :: zfor_snow_mult(nvec) !< Snow surface forcing [W/m^2].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  INTEGER :: i
  INTEGER :: ksn

  LOGICAL :: snow_present_flag(nvec)

  REAL(wp) :: zrnet_snow

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA PRESENT(ivend) CREATE(snow_present_flag) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    snow_present_flag(i) = (w_snow_now(i) * rho_w + dt_w_snow(i)*dt > rho_w * eps_soil)
  END DO

  ! heat conductivity of snow as funtion of water content
  !$ACC LOOP SEQ
  DO ksn = 1, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (snow_present_flag(i)) THEN
        zalas_mult(i,ksn) = 2.22_wp*EXP(1.88_wp*LOG(rho_snow_mult_now(i,ksn)/rho_i))
      END IF
    END DO
  END DO


  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrnet_snow)
  DO i = ivstart, ivend
    hfl_snow_soil(i) = 0._wp
    IF (snow_present_flag(i)) THEN
      hfl_snow_soil(i) = ((zalas_mult(i,ke_snow)*(-zhm_snow(i,ke_snow))+hcond_ml(i,1)*dz_ml(1))/ &
                 (-zhm_snow(i,ke_snow)+dz_ml(1)) * &
                 (t_snow_mult_now(i,ke_snow) - t_so_now(i,1))/(-zhm_snow(i,ke_snow) &
                 +dz_ml(1)))*fr_snow(i)

      !  ! GZ: use formulation of single-layer snow model, which is numerically more stable
      !  hfl_snow_soil(i) = zalas_mult(i,ke_snow)*(t_snow_mult_now(i,ke_snow) - t_so_now(i,1))/ &
      !  MAX(-zhm_snow(i,ke_snow),cdsmin)

    END IF

    ! total forcing for uppermost soil layer
    !<em new solution
    !  forcing_soil(i) = ( radfl_net_snfr(i) + shfl_snfr(i) + lhfl_snfr(i) + lhfl_precip(i) ) * (1._wp - fr_snow(i)) &
    !                    + (1._wp-ztsnow_pm(i)) * hfl_snow_soil(i)
    forcing_soil(i) = ( radfl_net_snfr(i) + shfl_snfr(i) + lhfl_snfr(i) + lhfl_precip(i) ) * (1._wp - fr_snow(i))
    !em>

    IF(snow_present_flag(i)) THEN
      IF(zextinct(i,1) > 0.0_wp) THEN
        zrnet_snow = radfl_th_snow(i)
      ELSE
        zrnet_snow = sobs(i) + radfl_th_snow(i)
      END IF
    ELSE
      zrnet_snow = sobs(i) + radfl_th_snow(i)
    END IF
    shfl_snow(i) = rho_ch(i)*cp_d*(th_atm(i) - t_snow_mult_now(i,1))
    lhfl_snow(i) = lh_s*evapo_snow(i) / MAX(eps_div, fr_snow(i))
    qhfl_snow(i) = evapo_snow(i) / MAX(eps_div, fr_snow(i))
    zfor_snow_mult(i)  = (zrnet_snow + shfl_snow(i) + lhfl_snow(i) + lh_f*rain_dew_rate(i))*fr_snow(i)
  END DO

  !$ACC END PARALLEL
  !$ACC END DATA
END SUBROUTINE snow_multi_soil_forcing


!> Calculate heat conduction through multilayer snow.
SUBROUTINE snow_multi_calc_heat_conduction ( &
      & ivstart, ivend, nvec, ke_snow, ke_soil, dt, t_s_now, t_so_now, t_so_new, t_so_free_new, &
      & hcond_ml, hcond_hl, hcap_ml, z_ml, dz_ml, dz_hl, t_snow_mult_now, t_snow_mult_new, &
      & dt_t_snow_mult, wliq_snow_now, wtot_snow_now, rho_snow_mult_now, zalas_mult, zhm_snow, &
      & zdzh_snow, zdzm_snow, zfor_snow_mult, sn_frac, w_snow_now, w_snow_new, dt_w_snow, &
      & hfl_snow_soil, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: ke_snow
  INTEGER, INTENT(IN) :: ke_soil

  REAL(wp), INTENT(IN) :: dt !< Time step [s].

  REAL(wp), INTENT(IN) :: t_s_now(nvec) !< Soil surface temperature (current) [K].
  REAL(wp), INTENT(IN) :: t_so_now(nvec, 0:ke_soil+1) !< Soil temperature (current) [K].
  REAL(wp), INTENT(INOUT) :: t_so_new(nvec,0:ke_soil+1) !< Soil temperature (new) [K].
  REAL(wp), INTENT(IN) :: t_so_free_new(nvec,0:ke_soil+1) !< Snow-free soil temperature [K].
  REAL(wp), INTENT(IN) :: hcond_ml(nvec,ke_soil+1) !< Heat conductivity in soil layers [W/(K m)].
  REAL(wp), INTENT(IN) :: hcond_hl(nvec,ke_soil) !< Heat conductivity across half levels [W/(K m)]
  REAL(wp), INTENT(IN) :: hcap_ml(nvec,ke_soil+1) !< Soil heat capacity [J/(K m^3)].
  REAL(wp), INTENT(IN) :: z_ml(ke_soil+1) !< Soil layer depth [m].
  REAL(wp), INTENT(IN) :: dz_ml(ke_soil+1) !< Layer main level distance (n_soil) [m].
  REAL(wp), INTENT(IN) :: dz_hl(ke_soil+1) !< Layer interface distance (n_soil) [m].

  REAL(wp), INTENT(IN) :: t_snow_mult_now(nvec,0:ke_snow) !< Snow temperature (current) [K].
  REAL(wp), INTENT(INOUT) :: t_snow_mult_new(nvec,0:ke_snow) !< Snow temperature (new) [K].
  REAL(wp), INTENT(INOUT) :: dt_t_snow_mult(nvec,0:ke_snow) !< Snow temperature tendency [K/s].

  REAL(wp), INTENT(IN) :: wliq_snow_now(nvec,ke_snow) !< Liquid snow water equivalent in layers [m].
  REAL(wp), INTENT(IN) :: wtot_snow_now(nvec,ke_snow) !< Total snow water equivalent in layers [m].
  REAL(wp), INTENT(IN) :: rho_snow_mult_now(nvec,ke_snow) !< Snow density in layers [kg/m**3].
  REAL(wp), INTENT(IN) :: zalas_mult(nvec,ke_snow) !< Snow heat conductivity [W/(K m)].
  REAL(wp), INTENT(IN) :: zhm_snow(nvec,ke_snow) !< Depth of snow main levels [m].
  REAL(wp), INTENT(IN) :: zdzh_snow(nvec,ke_snow) !< Snow layer thickness [m].
  REAL(wp), INTENT(IN) :: zdzm_snow(nvec,ke_snow) !< Distance between main snow layers [m].
  REAL(wp), INTENT(IN) :: zfor_snow_mult(nvec) !< Total forcing on snow (multilayer model) [W/m**2].

  REAL(wp), INTENT(IN) :: sn_frac(nvec) !< Snow fraction [1].
  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Total snow water equivalent (current) [m].
  REAL(wp), INTENT(INOUT) :: w_snow_new(nvec) !< Total snow water equivalent (new) [m].
  REAL(wp), INTENT(INOUT) :: dt_w_snow(nvec) !< Rate of change of SWE [m/s].

  REAL(wp), INTENT(INOUT) :: hfl_snow_soil(nvec) !< Heat flux through snow [W/m**2].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag.
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue.

  LOGICAL :: snow_present_flag(nvec)

  REAL(wp) :: zswitch(nvec) !< Threshold for thin snow [mH2O].

  REAL(wp) :: lse_a(nvec,0:ke_soil+ke_snow+1) !< Subdiagonal matrix entry.
  REAL(wp) :: lse_b(nvec,0:ke_soil+ke_snow+1) !< Diagonal matrix entry.
  REAL(wp) :: lse_c(nvec,0:ke_soil+ke_snow+1) !< Superdiagonal matrix entry.
  REAL(wp) :: lse_rhs(nvec,0:ke_soil+ke_snow+1) !< Right-hand side.
  REAL(wp) :: lse_sol(nvec,0:ke_soil+ke_snow+1) !< Solution of the linear system.

  INTEGER :: i
  INTEGER :: ksn
  INTEGER :: kso

  REAL(wp) :: hcap_sl !< Volumetric heat capacity of snow layer [J/(K m^3)].
  REAL(wp) :: zakb !< Heat diffusion coefficient [m^2/s].
  REAL(wp) :: zakb1 !< Heat diffusion coefficient to upper layer [m^2/s].
  REAL(wp) :: zakb2 !< Heat diffusion coefficient to lower layer [m^2/s].

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA CREATE(snow_present_flag, zswitch, lse_a, lse_b, lse_c, lse_rhs, lse_sol) &
  !$ACC   PRESENT(ivend) ASYNC(acc_async_queue) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(hcap_sl, zakb)
  DO i = ivstart, ivend
    snow_present_flag(i) = sn_frac(i) > 0._wp

    IF (snow_present_flag(i)) THEN
      ! Uppermost snow layer, Neumann boundary condition
      hcap_sl = (wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w + &
            (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i)*rho_snow_mult_now(i,1)
      zakb      = zalas_mult(i,1)/hcap_sl

      lse_a(i,1) = -zalfa*dt*zakb/(zdzh_snow(i,1)*zdzm_snow(i,1))
      lse_c(i,1) = -zalfa*dt*zakb/(zdzh_snow(i,1)*zdzm_snow(i,2))
      lse_b(i,1) = 1._wp - lse_a(i,1) - lse_c(i,1)
      lse_rhs(i,1) = t_snow_mult_now(i,1) + (1._wp - zalfa)* (-lse_a(i,1)/zalfa * t_s_now(i) + &
                  (lse_a(i,1) + lse_c(i,1))/zalfa * t_so_now(i,1) - lse_c(i,1)/zalfa * t_so_now(i,2) )
      lse_a(i,0) = 0.0_wp
      lse_b(i,0) = zalfa
      lse_c(i,0) = -zalfa
      lse_rhs(i,0) = zdzm_snow(i,1) * zfor_snow_mult(i)/zalas_mult(i,1)+(1._wp-zalfa)* &
                 (t_so_now(i,1) - t_s_now(i))

      lse_c(i,0) = lse_c(i,0)/lse_b(i,0)
      lse_rhs(i,0) = lse_rhs(i,0)/lse_b(i,0)

      ! Lowermost soil layer, Dirichlet boundary condition
      lse_a(i,ke_snow + ke_soil+1) = 0.0_wp
      lse_b(i,ke_snow + ke_soil+1) = 1.0_wp
      lse_c(i,ke_snow + ke_soil+1) = 0.0_wp
      lse_rhs(i,ke_snow + ke_soil+1) = t_so_now(i,ke_soil+1)

      ! Lowermost snow layer, special treatment
      hcap_sl = (wliq_snow_now(i,ke_snow)/wtot_snow_now(i,ke_snow)*chc_w +  &
        & (wtot_snow_now(i,ke_snow)-wliq_snow_now(i,ke_snow))/                    &
        & wtot_snow_now(i,ke_snow)*chc_i)*rho_snow_mult_now(i,ke_snow)
      zakb = zalas_mult(i,ke_snow)/hcap_sl
      lse_a(i,ke_snow) = -zalfa*dt*zakb/(zdzh_snow(i,ke_snow)*zdzm_snow(i,ke_snow))
      zakb = (zalas_mult(i,ke_snow)/hcap_sl*(-zhm_snow(i,ke_snow))+hcond_ml(i,1)/ &
              hcap_ml(i,1)*z_ml(1))/(z_ml(1)-zhm_snow(i,ke_snow))
      lse_c(i,ke_snow) = -zalfa*dt*zakb/(zdzh_snow(i,ke_snow)*(z_ml(1)-zhm_snow(i,ke_snow)))
      lse_b(i,ke_snow) = 1.0_wp - lse_a(i,ke_snow) - lse_c(i,ke_snow)
      lse_rhs(i,ke_snow) = t_snow_mult_now(i,ke_snow) + &
        &     (1._wp - zalfa)*( - lse_a(i,ke_snow)/zalfa*t_snow_mult_now(i,ke_snow-1) + &
        &     (lse_a(i,ke_snow)/zalfa + lse_c(i,ke_snow)/zalfa)*t_snow_mult_now(i,ke_snow) - &
        &     lse_c(i,ke_snow)/zalfa*t_so_now(i,1)  )

      ! Uppermost soil layer, special treatment
      zakb = (zalas_mult(i,ke_snow)/hcap_sl*(-zhm_snow(i,ke_snow))+     &
              hcond_ml(i,1)/hcap_ml(i,1)*z_ml(1))/(z_ml(1)-zhm_snow(i,ke_snow))
      lse_a(i,ke_snow+1) = -zalfa*dt*zakb/(dz_hl(1)*(z_ml(1)-zhm_snow(i,ke_snow)))
      zakb = hcond_ml(i,1)/hcap_ml(i,1)
      lse_c(i,ke_snow+1) = -zalfa*dt*zakb/(dz_hl(1)*dz_ml(2))
      lse_b(i,ke_snow+1) = 1._wp - lse_a(i,ke_snow+1) - lse_c(i,ke_snow+1)
      lse_rhs(i,ke_snow+1) = t_so_now(i,1) + &
        &    (1._wp - zalfa)*( - lse_a(i,ke_snow+1)/zalfa*t_snow_mult_now(i,ke_snow) + &
        &    (lse_a(i,ke_snow+1)/zalfa + lse_c(i,ke_snow+1)/zalfa)*t_so_now(i,1)      - &
        &    lse_c(i,ke_snow+1)/zalfa*t_so_now(i,2)  )
    END IF
  END DO

  ! Snow layers
  !$ACC LOOP SEQ
  DO ksn = 2, ke_snow-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(hcap_sl, zakb)
    DO i = ivstart, ivend
      IF (snow_present_flag(i)) THEN
        hcap_sl = (wliq_snow_now(i,ksn)/wtot_snow_now(i,ksn)*chc_w + &
          & (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn))/ &
          & wtot_snow_now(i,ksn)*chc_i)*rho_snow_mult_now(i,ksn)
        zakb = zalas_mult(i,ksn)/hcap_sl
        lse_a(i,ksn) = -zalfa*dt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn))
        zakb = (zalas_mult(i,ksn)/hcap_sl*zdzm_snow(i,ksn)+zalas_mult(i,ksn+1)/hcap_sl*zdzm_snow(i,ksn+1))/ &
               (zdzm_snow(i,ksn)+zdzm_snow(i,ksn+1))
        lse_c(i,ksn) = -zalfa*dt*zakb/(zdzh_snow(i,ksn)*zdzm_snow(i,ksn+1))
        lse_b(i,ksn) = 1._wp - lse_a(i,ksn) - lse_c(i,ksn)
        lse_rhs(i,ksn) = t_snow_mult_now(i,ksn) + &
          &     (1._wp - zalfa)*( - lse_a(i,ksn)/zalfa*t_snow_mult_now(i,ksn-1) + &
          &     (lse_a(i,ksn)/zalfa + lse_c(i,ksn)/zalfa)*t_snow_mult_now(i,ksn)     - &
          &     lse_c(i,ksn)/zalfa*t_snow_mult_now(i,ksn+1)  )
      END IF
    END DO
  END DO                ! snow layers

  ! Soil layers
  !$ACC LOOP SEQ
  DO ksn = ke_snow+2, ke_snow + ke_soil
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(kso, zakb1, zakb2)
    DO i = ivstart, ivend
      IF (snow_present_flag(i)) THEN
        kso = ksn - ke_snow
        zakb1 = hcond_hl(i,kso-1)/hcap_ml(i,kso)
        zakb2 = hcond_hl(i,kso  )/hcap_ml(i,kso)
        lse_a(i,ksn) = -zalfa*dt*zakb1/(dz_hl(kso)*dz_ml(kso))
        lse_c(i,ksn) = -zalfa*dt*zakb2/(dz_hl(kso)*dz_ml(kso+1))
        lse_b(i,ksn) = 1._wp - lse_a(i,ksn) - lse_c(i,ksn)
        lse_rhs(i,ksn) = t_so_now(i,kso) + &
          &    (1._wp - zalfa)*( - lse_a(i,ksn)/zalfa*t_so_now(i,kso-1) +  &
          &    (lse_a(i,ksn)/zalfa + lse_c(i,ksn)/zalfa)*t_so_now(i,kso)     -  &
          &    lse_c(i,ksn)/zalfa*t_so_now(i,kso+1)  )
      END IF
    END DO
  END DO                ! soil layers

  CALL solve_tridiag ( &
      & ivstart=ivstart, &
      & ivend=ivend, &
      & nvec=nvec, &
      & nlev=ke_snow + ke_soil + 2, & ! 0:ke_snow+ke_soil+1
      & a=lse_a(:,0:), &
      & b=lse_b(:,0:), &
      & c=lse_c(:,0:), &
      & d=lse_rhs(:,0:), &
      & out=lse_sol(:,0:), &
      & mask=snow_present_flag(:) &
    )

  !$ACC LOOP SEQ
  DO ksn = 1, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (snow_present_flag(i)) THEN
        t_snow_mult_new(i,ksn) = lse_sol(i,ksn)
      END IF
    END DO
  END DO

  !$ACC LOOP SEQ
  DO kso = 1, ke_soil+1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (snow_present_flag(i)) THEN
        t_so_new(i,kso) = lse_sol(i,ke_snow + kso)*sn_frac(i) + t_so_free_new(i,kso)*(1._wp - sn_frac(i))
      ELSE
        t_so_new(i,kso) = t_so_free_new(i,kso)
      END IF
    END DO
  END DO

  !in case of thin snowpack (less than zswitch), apply single-layer snow model
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(hcap_sl)
  DO i = ivstart, ivend
    IF (snow_present_flag(i)) THEN
      t_snow_mult_new(i,0) = t_snow_mult_new(i,1)

      hfl_snow_soil(i) = ((zalas_mult(i,ke_snow)*(-zhm_snow(i,ke_snow))+hcond_ml(i,1)*dz_ml(1))/ &
                (-zhm_snow(i,ke_snow)+dz_ml(1)) * &
                (t_snow_mult_new(i,ke_snow) - t_so_now(i,1))/(-zhm_snow(i,ke_snow) &
                +dz_ml(1)))*sn_frac(i)

      hcap_sl = (wliq_snow_now(i,1)/wtot_snow_now(i,1)*chc_w + &
        (wtot_snow_now(i,1)-wliq_snow_now(i,1))/wtot_snow_now(i,1)*chc_i)*rho_snow_mult_now(i,1)
      zswitch(i) = (-zfor_snow_mult(i)+hfl_snow_soil(i))/50./hcap_sl*dt*ke_snow
      zswitch(i) = MAX(zswitch(i),1.E-03_wp)
    ELSE
      zswitch(i) = 0.0_wp
    END IF
  END DO

  !$ACC LOOP SEQ
  DO ksn = 1, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(w_snow_new(i) > eps_soil .AND. w_snow_new(i) < zswitch(i)) THEN

        IF((zfor_snow_mult(i)-hfl_snow_soil(i))*dt > w_snow_new(i)*rho_w*lh_f) THEN
          t_snow_mult_new(i,ksn) = t_so_now(i,0)
        ELSE IF(zfor_snow_mult(i)-hfl_snow_soil(i) > 0._wp) THEN
          t_snow_mult_new(i,ksn) = t_snow_mult_now(i,ksn) + &
            (zfor_snow_mult(i)-hfl_snow_soil(i))*dt/(chc_i*wtot_snow_now(i,ksn))/rho_w/ke_snow
        END IF
      END IF
    END DO
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF(w_snow_new(i) > eps_soil .AND. w_snow_new(i) < zswitch(i)) THEN

      IF((zfor_snow_mult(i)-hfl_snow_soil(i))*dt > w_snow_new(i)*rho_w*lh_f) THEN
        dt_w_snow(i) = dt_w_snow(i) - w_snow_new(i)*rho_w/dt
        w_snow_new(i)  = 0._wp
        t_snow_mult_new(i,0) = t_so_now(i,0)
      ELSE IF((zfor_snow_mult(i)-hfl_snow_soil(i)) > 0._wp) THEN
        t_snow_mult_new(i,0) = t_snow_mult_new(i,1)
      END IF

    END IF
  END DO

  !$ACC LOOP SEQ
  DO ksn = 0, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_now(i) + dt_w_snow(i) * (dt / rho_w) > eps_soil) THEN
        dt_t_snow_mult(i,ksn) = (t_snow_mult_new(i,ksn) - t_snow_mult_now(i,ksn)) / dt
      END IF
    END DO
  ENDDO

  !$ACC END PARALLEL
  !$ACC END DATA

END SUBROUTINE snow_multi_calc_heat_conduction


SUBROUTINE snow_multi_melt ( &
      & ivstart, ivend, nvec, ke_snow, dt, soiltyp_subs, dz_top, w_snow_now, zf_snow, dt_w_snow, &
      & dt_w_so_top, fr_w_top, runoff_s, sobs, zextinct, t_snow_mult_new, dt_t_snow_mult, &
      & zdzh_snow, wtot_snow_now, wliq_snow_now, rho_snow_mult_now, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: ke_snow

  REAL(wp), INTENT(IN) :: dt

  INTEGER, INTENT(IN) :: soiltyp_subs(nvec) !< Soil type.

  REAL(wp), INTENT(IN) :: dz_top !< Top layer thickness [m].

  REAL(wp), INTENT(IN) :: w_snow_now(nvec) !< Current SWE [m].
  REAL(wp), INTENT(IN) :: zf_snow(nvec) !< Snow fraction [m**2/m**2].
  REAL(wp), INTENT(INOUT) :: dt_w_snow(nvec) !< snow tendency [kg/(m**2 s)].
  REAL(wp), INTENT(INOUT) :: dt_w_so_top(nvec) !< Top soil water tendency [kg/(m**2 s)].
  REAL(wp), INTENT(IN) :: fr_w_top(nvec) !< Total water fraction in top soil layer [m**3/m**3].
  REAL(wp), INTENT(INOUT) :: runoff_s(nvec) !< Surface runoff (accumulated) [kg/m**2].

  REAL(wp), INTENT(IN) :: sobs(nvec) !< Net surface short-wave radiation [W/m**2].
  REAL(wp), INTENT(IN) :: zextinct(nvec,ke_snow) !< Extinction length scale [m**-1].
  REAL(wp), INTENT(IN) :: t_snow_mult_new(nvec,0:ke_snow) !< Snow temperature [K]
  REAL(wp), INTENT(INOUT) :: dt_t_snow_mult(nvec,0:ke_snow) !< Snow layer temperature tendency [K/s].
  REAL(wp), INTENT(INOUT) :: zdzh_snow(nvec,ke_snow) !< Snow layer thickness [m].
  REAL(wp), INTENT(INOUT) :: wtot_snow_now(nvec,ke_snow) !< Total SWE [m].
  REAL(wp), INTENT(INOUT) :: wliq_snow_now(nvec,ke_snow) !< Liquid water content [m].
  REAL(wp), INTENT(INOUT) :: rho_snow_mult_now(nvec,ke_snow) !< Snow density [kg/m**3].

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp) :: w_snow_prov(nvec) !< Provisional updated SWE [m].

  REAL(wp) :: field_cap
  REAL(wp) :: pore_vol

  REAL(wp) :: zdwsnm
  REAL(wp) :: zro
  REAL(wp) :: zredfu
  REAL(wp) :: zw_ovpv

  REAL(wp) :: ztsnownew_mult(nvec,0:ke_snow)
  REAL(wp) :: ze_in
  REAL(wp) :: ze_out(nvec)
  REAL(wp) :: ze_rad(nvec)
  REAL(wp) :: zqbase(nvec)
  REAL(wp) :: zp(nvec,ke_snow)
  REAL(wp) :: zcounter
  REAL(wp) :: zrefr
  REAL(wp) :: zmelt
  REAL(wp) :: zrho_dry_old
  REAL(wp) :: zsn_porosity
  REAL(wp) :: zp1
  REAL(wp) :: zfukt
  REAL(wp) :: zadd_dz
  REAL(wp) :: zdens_old
  REAL(wp) :: zeta
  REAL(wp) :: zq0

  REAL(wp) :: dt_recip

  INTEGER :: i, k, ksn

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA ASYNC(acc_async_queue) IF(lzacc) &
  !$ACC   PRESENT(ivend) CREATE(w_snow_prov, ztsnownew_mult, ze_out, ze_rad, zqbase, zp)

  dt_recip = 1._wp / dt

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend

    w_snow_prov(i) = w_snow_now(i) + dt_w_snow(i) * dt / rho_w

    ze_out  (i) = 0.0_wp
    zqbase  (i) = 0.0_wp

    ze_rad(i) = 0.0_wp
    IF(zextinct(i,1).gt.0.0_wp) ze_rad(i) = zf_snow(i) * sobs(i)

    ztsnownew_mult(i,0) = t_snow_mult_new(i,0)
  END DO

  !$ACC LOOP SEQ
  DO ksn = 1,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zrefr, zmelt, zrho_dry_old, ze_in, zcounter) &
    !$ACC   PRIVATE(zadd_dz, zsn_porosity, zp1, zfukt, zq0)
    DO i = ivstart, ivend
      IF (w_snow_prov(i) > eps_soil) THEN        ! points with snow cover only

        zrefr = 0.0_wp
        zmelt = 0.0_wp
        ztsnownew_mult(i,ksn) = t_snow_mult_new(i,ksn)

        IF(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn).GT.eps_soil .OR. &
          wtot_snow_now(i,ksn) - wliq_snow_now(i,ksn).GT.eps_soil) THEN
          zrho_dry_old = MAX(wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn), &
            &                     eps_soil)                                               &
            &                 *rho_w/(zdzh_snow(i,ksn) - wliq_snow_now(i,ksn))
        ELSE
          zrho_dry_old = rho_w
        END IF

        ztsnownew_mult(i,ksn) = (t_snow_mult_new(i,ksn)*wtot_snow_now(i,ksn) &
          &                       + t0_melt*zqbase(i)*dt)/(zqbase(i)*dt     &
          &                       + wtot_snow_now(i,ksn))

        IF(zextinct(i,ksn).eq.0.0_wp) THEN
          ze_in = ze_out(i)
        ELSE
          IF(ksn.eq.ke_snow) THEN     ! all the rest of radiation is absorbed by the lowermost snow layer
            ze_in = ze_out(i) + ze_rad(i)
          ELSE
            zcounter = EXP (-zextinct(i,ksn)*zdzh_snow(i,ksn))
            ze_in = ze_out(i) + ze_rad(i) * (1._wp - zcounter)
            ze_rad(i) = ze_rad(i) * zcounter
          END IF
        END IF

        ztsnownew_mult(i,ksn) = ztsnownew_mult(i,ksn) &
          &                       + ze_in*dt/(chc_i*wtot_snow_now(i,ksn))/rho_w
        wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) + zqbase(i)*dt
        wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zqbase(i)*dt

        zdzh_snow(i,ksn) = zdzh_snow(i,ksn) + zqbase(i)*dt

        rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*&
          &                                rho_w/zdzh_snow(i,ksn), &
          &                                0.0_wp)

        IF(ztsnownew_mult(i,ksn) .GT. t0_melt) THEN

          IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
            ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
              &      *dt_recip*rho_w
            zmelt = 0.0_wp
          ELSEIF(chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt)/lh_f <= &
            wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn)) THEN
            zmelt = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
              &          *dt_recip/lh_f
            ze_out(i) = 0.0_wp
            wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt*dt
          ELSE
            zmelt = (wtot_snow_now(i,ksn)-wliq_snow_now(i,ksn)) * dt_recip
            ze_out(i) = chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn)-t0_melt) &
              &      *dt_recip*rho_w - zmelt*lh_f*rho_w
            wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) + zmelt*dt
          END IF
          ztsnownew_mult(i,ksn) = t0_melt

        ELSE
          ! T<0
          IF(wliq_snow_now(i,ksn) .GT. -chc_i*wtot_snow_now(i,ksn) &
            & *(ztsnownew_mult(i,ksn) - t0_melt)/lh_f) THEN
            zrefr = -chc_i*wtot_snow_now(i,ksn)*(ztsnownew_mult(i,ksn) &
              &          - t0_melt)*dt_recip/lh_f
            ztsnownew_mult(i,ksn)   = t0_melt
            wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zrefr*dt
          ELSE
            zrefr = wliq_snow_now(i,ksn) * dt_recip
            wliq_snow_now(i,ksn) = 0.0_wp
            ztsnownew_mult(i,ksn)   = ztsnownew_mult(i,ksn) + zrefr*dt*lh_f &
              &                         /(chc_i*wtot_snow_now(i,ksn))
          END IF
          ze_out(i) = 0.0_wp

        END IF

        dt_t_snow_mult(i,ksn) = dt_t_snow_mult(i,ksn) + &
                                  (ztsnownew_mult(i,ksn) - t_snow_mult_new(i,ksn))*dt_recip
        IF(wtot_snow_now(i,ksn) .LE. wliq_snow_now(i,ksn)) THEN
          zqbase(i)           = wliq_snow_now(i,ksn)*dt_recip
          wliq_snow_now(i,ksn) = 0.0_wp
          wtot_snow_now(i,ksn) = 0.0_wp
          zdzh_snow(i,ksn)    = 0.0_wp
          rho_snow_mult_now(i,ksn)  = 0.0_wp
        ELSE
          IF(zrefr .GT. 0.0_wp .OR. zmelt .GT. 0.0_wp) THEN
            zadd_dz = 0.0_wp
            zadd_dz = MAX(zrefr,0._wp)*(-1.0_wp + 1.0_wp/rho_i*rho_w)*dt
            zadd_dz = MAX(zmelt,0._wp)*(-1.0_wp/zrho_dry_old*rho_w &
              &       + 1.0_wp)*dt
            zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) + zadd_dz
            rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
              &                            /zdzh_snow(i,ksn),0.0_wp)
            IF(wtot_snow_now(i,ksn) .LE. 0.0_wp) zdzh_snow(i,ksn) = 0.0_wp
            IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
              zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
              rho_snow_mult_now(i,ksn) = rho_w
            END IF
          END IF

          zsn_porosity = 1._wp - (rho_snow_mult_now(i,ksn)/rho_w -  &
                         wliq_snow_now(i,ksn)/zdzh_snow(i,ksn))/rho_i*rho_w - &
                         wliq_snow_now(i,ksn)/zdzh_snow(i,ksn)
          zsn_porosity = MAX(zsn_porosity,cwhc + 0.1_wp)
          zp1 = zsn_porosity - cwhc

          IF (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) .GT. cwhc) THEN
            zfukt             = (wliq_snow_now(i,ksn)/zdzh_snow(i,ksn) - cwhc)/zp1
            zq0               = chcond * zfukt**3
            zqbase(i)       = MIN(zq0*dt,wliq_snow_now(i,ksn))
            wliq_snow_now(i,ksn) = wliq_snow_now(i,ksn) - zqbase(i)
            wtot_snow_now(i,ksn) = wtot_snow_now(i,ksn) - zqbase(i)

            zdzh_snow(i,ksn) = zdzh_snow(i,ksn) - zqbase(i)
            zqbase(i)        = zqbase(i) * dt_recip

            IF(zdzh_snow(i,ksn) .LT. eps_soil*0.01_wp) THEN
              wliq_snow_now(i,ksn) = 0.0_wp
              wtot_snow_now(i,ksn) = 0.0_wp
              zdzh_snow(i,ksn)     = 0.0_wp
              rho_snow_mult_now(i,ksn)  = 0.0_wp
            ELSE
              rho_snow_mult_now(i,ksn) = MAX(wtot_snow_now(i,ksn)*rho_w &
                &                            /zdzh_snow(i,ksn),0.0_wp)
              IF(wtot_snow_now(i,ksn) .LE. 0.0_wp) zdzh_snow(i,ksn) = 0.0_wp
              IF(rho_snow_mult_now(i,ksn) .GT. rho_w) THEN
                zdzh_snow(i,ksn)   = zdzh_snow(i,ksn)*rho_snow_mult_now(i,ksn)/rho_w
                rho_snow_mult_now(i,ksn) = rho_w
              END IF
            END IF
          ELSE
            zqbase(i) = 0.0_wp
          END IF
        END IF

      END IF       ! points with snow cover only
    END DO
  END DO        ! snow layers

  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdwsnm, pore_vol, field_cap, zredfu, zro, zw_ovpv)
  DO i = ivstart, ivend
    IF (w_snow_prov(i) > eps_soil) THEN        ! points with snow cover only
      zdwsnm = zqbase(i)*rho_w       ! ksn == ke_snow
      dt_w_snow(i)  = dt_w_snow(i) - zdwsnm

      pore_vol = cporv(soiltyp_subs(i))
      field_cap = cfcap(soiltyp_subs(i))

      ! melted snow is allowed to penetrate the soil (up to field
      ! capacity), if the soil type is neither ice nor rock;
      ! else it contributes to surface run-off;
      ! fractional water content of the first soil layer determines
      ! a reduction factor which controls additional run-off

      IF (soiltyp_subs(i) > IST_ROCK) THEN
        zredfu = MAX( &
            & 0.0_wp, &
            & MIN(1.0_wp, &
            &   (fr_w_top(i) - field_cap)/MAX(pore_vol-field_cap, eps_soil) &
            & ) &
          )
        dt_w_so_top(i) = dt_w_so_top(i) + zdwsnm*(1._wp - zredfu)
        zro = zdwsnm * zredfu  ! Infiltration not possible
                               ! for this fraction
      ELSE
        zro = zdwsnm ! surface runoff
      END IF

      ! zro-, zdw_so_dt-correction in case of pore volume overshooting
      zw_ovpv = MAX(0._wp, (fr_w_top(i) - pore_vol) * dz_top * rho_w / dt +  &
                dt_w_so_top(i))
      zro = zro + zw_ovpv
      dt_w_so_top(i)= dt_w_so_top(i) - zw_ovpv

      runoff_s(i) = runoff_s(i) + zro * dt
    END IF       ! points with snow cover only
  END DO

  ! snow densification due to gravity and metamorphism
  !$ACC LOOP SEQ
  DO ksn = 2, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      zp(i,ksn) = 0.0_wp                         ! gravity, Pa
    END DO

    !$ACC LOOP SEQ
    DO k = ksn,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i = ivstart, ivend
         zp(i,ksn) = zp(i,ksn) + rho_snow_mult_now(i,k)*g*zdzh_snow(i,ksn)
      END DO
    END DO
  END DO

  !$ACC LOOP SEQ
  DO ksn = 2, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zdens_old, zeta)
    DO i = ivstart, ivend
      IF (w_snow_prov(i) > eps_soil) THEN        ! points with snow cover only
        IF(rho_snow_mult_now(i,ksn) .LT. 600._wp .AND. &
          rho_snow_mult_now(i,ksn) .NE. 0.0_wp) THEN
          zdens_old = rho_snow_mult_now(i,ksn)
          zeta =         &! compactive viscosity of snow
            ca2*EXP(19.3_wp*rho_snow_mult_now(i,ksn)/rho_i)* &
            EXP(67300._wp/8.31_wp/ztsnownew_mult(i,ksn))
          rho_snow_mult_now(i,ksn) = rho_snow_mult_now(i,ksn) + &
            dt*rho_snow_mult_now(i,ksn)*(csigma+zp(i,ksn))/zeta
          rho_snow_mult_now(i,ksn) = MIN(rho_snow_mult_now(i,ksn),rho_i)
          zdzh_snow(i,ksn)   = zdzh_snow(i,ksn) * zdens_old/rho_snow_mult_now(i,ksn)
        END IF
      END IF       ! points with snow cover only
    END DO
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (w_snow_prov(i) > eps_soil) THEN        ! points with snow cover only

      IF(ztsnownew_mult(i,0) .GT. t0_melt) THEN
        ztsnownew_mult(i,0) = t0_melt
        dt_t_snow_mult(i,0) = dt_t_snow_mult(i,0) +     &
                                (ztsnownew_mult(i,0) - t_snow_mult_new(i,0)) / dt
      END IF
    END IF       ! points with snow cover only
  END DO
  !$ACC END PARALLEL

  !$ACC END DATA

END SUBROUTINE snow_multi_melt


SUBROUTINE snow_multi_update_new_state ( &
      & ivstart, ivend, nvec, ke_snow, dt, t_snow_new, w_snow_new, w_snow_now, &
      & dt_w_snow, rho_snow_new, h_snow_new, t_snow_mult_new, t_snow_mult_now, dt_t_snow_mult, &
      & dzh_snow_new, zdzh_snow, wtot_snow_new, wtot_snow_now, rho_snow_mult_new, &
      & rho_snow_mult_now, wliq_snow_new, wliq_snow_now, t_so_new_top, w_i_new, &
      & zhh_snow, zhm_snow, zdzm_snow, lzacc, acc_async_queue &
    )

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: ke_snow

  REAL(wp), INTENT(IN) :: dt

  REAL(wp), INTENT(INOUT) :: t_snow_new(nvec)

  REAL(wp), INTENT(INOUT) :: w_snow_new(nvec)
  REAL(wp), INTENT(IN) :: w_snow_now(nvec)
  REAL(wp), INTENT(IN) :: dt_w_snow(nvec)

  REAL(wp), INTENT(INOUT) :: rho_snow_new(nvec)

  REAL(wp), INTENT(INOUT) :: h_snow_new(nvec)

  REAL(wp), INTENT(INOUT) :: t_snow_mult_new(nvec, 0:ke_snow)
  REAL(wp), INTENT(IN) :: t_snow_mult_now(nvec, 0:ke_snow)
  REAL(wp), INTENT(IN) :: dt_t_snow_mult(nvec, 0:ke_snow)

  REAL(wp), INTENT(INOUT) :: dzh_snow_new(nvec, ke_snow)
  REAL(wp), INTENT(IN) :: zdzh_snow(nvec, ke_snow)

  REAL(wp), INTENT(INOUT) :: wtot_snow_new(nvec, ke_snow)
  REAL(wp), INTENT(IN) :: wtot_snow_now(nvec, ke_snow)

  REAL(wp), INTENT(INOUT) :: rho_snow_mult_new(nvec, ke_snow)
  REAL(wp), INTENT(IN) :: rho_snow_mult_now(nvec, ke_snow)

  REAL(wp), INTENT(INOUT) :: wliq_snow_new(nvec, ke_snow)
  REAL(wp), INTENT(IN) :: wliq_snow_now(nvec, ke_snow)

  REAL(wp), INTENT(IN) :: t_so_new_top(nvec)

  REAL(wp), INTENT(INOUT) :: w_i_new(nvec)

  REAL(wp), INTENT(INOUT) :: zhh_snow(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: zhm_snow(nvec,ke_snow)
  REAL(wp), INTENT(INOUT) :: zdzm_snow(nvec,ke_snow)

  LOGICAL, INTENT(IN) :: lzacc !< OpenACC flag
  INTEGER, INTENT(IN) :: acc_async_queue !< OpenACC queue number

  REAL(wp) :: weight

  REAL(wp) :: sum_weight(nvec)
  REAL(wp) :: dz_old(nvec, ke_snow)
  REAL(wp) :: z_old(nvec, ke_snow)
  REAL(wp) :: t_new(nvec, ke_snow)
  REAL(wp) :: rho_new(nvec, ke_snow)
  REAL(wp) :: wl_new(nvec, ke_snow)

  INTEGER :: i
  INTEGER :: k
  INTEGER :: ksn

  OPENACC_SUPPRESS_UNUSED_LZACC

  !$ACC DATA ASYNC(acc_async_queue) IF(lzacc) &
  !$ACC   PRESENT(ivend) CREATE(sum_weight, dz_old, z_old, t_new, rho_new, wl_new)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

  ! First for ksn == 0
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    t_snow_mult_new  (i,0) = t_snow_mult_now(i,0) + dt*dt_t_snow_mult(i,0)
    t_snow_new(i) = t_snow_mult_new (i,0)
  ENDDO

  !$ACC LOOP SEQ
  DO ksn = 1, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      t_snow_mult_new  (i,ksn) = t_snow_mult_now(i,ksn) + &
        &                              dt*dt_t_snow_mult(i,ksn)
      dzh_snow_new     (i,ksn) = zdzh_snow(i,ksn)
      wtot_snow_new    (i,ksn) = wtot_snow_now(i,ksn)
      rho_snow_mult_new(i,ksn) = rho_snow_mult_now(i,ksn)
      wliq_snow_new    (i,ksn) = wliq_snow_now(i,ksn)
    ENDDO
  ENDDO

  ! Reset t_snow_new to t_so(0) if no snow was present at the beginning of the time step
  ! The heat balance calculation is incomplete in this case and sometimes yields unreasonable results
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    w_snow_new(i) = w_snow_now(i) + dt*dt_w_snow(i)/rho_w

    IF (w_snow_now(i) < eps_soil .AND. w_snow_new(i) >= eps_soil) THEN
      t_snow_new(i) = MIN(t0_melt,t_so_new_top(i))
    ENDIF
  ENDDO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (w_snow_now(i) < eps_soil .AND. w_snow_new(i) >= eps_soil) THEN
      t_snow_mult_new(i,:) = t_snow_new(i)
    ENDIF
  ENDDO

  ! Eliminate snow for multi-layer snow model, if w_snow = 0
  !$ACC LOOP SEQ
  DO ksn = 1, ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_new(i) <= eps_soil) THEN
        t_snow_mult_new(i,ksn) = t_so_new_top(i)
        wliq_snow_new(i,ksn) = 0.0_wp
        wtot_snow_new(i,ksn) = 0.0_wp
        rho_snow_mult_new(i,ksn) = 0.0_wp
        dzh_snow_new(i,ksn) = 0.0_wp
      ENDIF
    END DO
  END DO

  ! Put small snow amounts to interception.
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (w_snow_new(i) <= eps_soil) THEN
      w_i_new(i)    = w_i_new(i) + w_snow_new(i)
      w_snow_new(i) = 0.0_wp
      t_snow_new(i) = t_so_new_top(i)
    ENDIF

    IF (w_i_new(i) <= 1.0E-4_wp*eps_soil) w_i_new(i) = 0.0_wp
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    h_snow_new(i) = 0.0_wp
    sum_weight(i) = 0.0_wp
  END DO

  !$ACC LOOP SEQ
  DO ksn = 1,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(w_snow_new(i) .GT. eps_soil) THEN
        h_snow_new(i) = h_snow_new(i) + zdzh_snow(i,ksn)
      END IF
    END DO
  END DO

  k = MIN(2,ke_snow-1)
  !$ACC LOOP SEQ
  DO ksn = 1,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_new(i) .GT. eps_soil) THEN
        IF (ksn == 1) THEN ! Limit top layer to max_toplaydepth
          zhh_snow(i,ksn) = -MAX( h_snow_new(i)-max_toplaydepth, h_snow_new(i)/ke_snow*(ke_snow-ksn) )
        ELSE IF (ksn == 2 .AND. ke_snow > 2) THEN ! Limit second layer to 8*max_toplaydepth
          zhh_snow(i,ksn) = MIN( 8._wp*max_toplaydepth+zhh_snow(i,1), zhh_snow(i,1)/(ke_snow-1)*(ke_snow-ksn) )
        ELSE ! distribute the remaining snow equally among the layers
          zhh_snow(i,ksn) = zhh_snow(i,k)/(ke_snow-k)*(ke_snow-ksn)
        ENDIF
      ENDIF
    END DO
  END DO

  !$ACC LOOP SEQ
  DO ksn = ke_snow,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (w_snow_new(i) .GT. eps_soil) THEN
        dz_old(i,ksn) = dzh_snow_new(i,ksn)
        z_old(i,ksn) = -sum_weight(i) - dzh_snow_new(i,ksn)/2._wp
        sum_weight(i) = sum_weight(i) + dzh_snow_new(i,ksn)
      END IF
    END DO
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
      IF(w_snow_new(i) .GT. eps_soil) THEN
        zhm_snow(i,1) = (-h_snow_new(i) + zhh_snow(i,1))/2._wp

        ! layer thickness betw. half levels of uppermost snow layer
        dzh_snow_new(i,1) = zhh_snow(i,1) + h_snow_new(i)

        ! layer thickness between snow surface and main level of uppermost layer
        zdzm_snow(i,1     ) = zhm_snow(i,1) + h_snow_new(i)

        IF(dz_old(i,1).ne.0..and.rho_snow_mult_new(i,1).ne.0.) THEN
          wliq_snow_new(i,1) = wliq_snow_new(i,1)/dz_old(i,1)
        END IF
      END IF
  END DO

  !$ACC LOOP SEQ
  DO ksn = 2,ke_snow
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(w_snow_new(i) .GT. eps_soil) THEN
        zhm_snow(i,ksn) = (zhh_snow(i,ksn) + zhh_snow(i,ksn-1))/2._wp
        dzh_snow_new(i,ksn) = zhh_snow(i,ksn) - zhh_snow(i,ksn-1) ! layer thickness betw. half levels
        zdzm_snow(i,ksn     ) = zhm_snow(i,ksn) - zhm_snow(i,ksn-1) ! layer thickness betw. main levels
        IF(dz_old(i,ksn).ne.0..and.rho_snow_mult_new(i,ksn).ne.0.) THEN
          wliq_snow_new(i,ksn) = wliq_snow_new(i,ksn)/dz_old(i,ksn)
        END IF
      END IF
    END DO
  END DO

  !$ACC LOOP SEQ
  DO ksn = ke_snow,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      t_new  (i,ksn) = 0.0_wp
      rho_new(i,ksn) = 0.0_wp
      wl_new (i,ksn) = 0.0_wp
    END DO

    !$ACC LOOP SEQ
    DO k = ke_snow,1,-1
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(weight)
      DO i = ivstart, ivend
        IF(w_snow_new(i) .GT. eps_soil) THEN

          weight = MIN(&
                dz_old(i,k), &
                z_old(i,k) + dz_old(i,k)/2._wp &
                - zhm_snow(i,ksn) + dzh_snow_new(i,ksn)/2._wp , &
                zhm_snow(i,ksn) + dzh_snow_new(i,ksn)/2._wp &
                - z_old(i,k) + dz_old(i,k)/2._wp, &
                dzh_snow_new(i,ksn))

          weight = (weight + ABS(weight)) * 0.5_wp / dzh_snow_new(i,ksn)

          t_new  (i,ksn) = t_new  (i,ksn) + t_snow_mult_new  (i,k)*weight
          rho_new(i,ksn) = rho_new(i,ksn) + rho_snow_mult_new(i,k)*weight
          wl_new (i,ksn) = wl_new (i,ksn) + wliq_snow_new(i,k)*weight
        END IF
      END DO
    END DO
  END DO

  !$ACC LOOP SEQ
  DO ksn = ke_snow,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF(w_snow_new(i) > eps_soil) THEN
        t_snow_mult_new  (i,ksn) = t_new  (i,ksn)
        rho_snow_mult_new(i,ksn) = rho_new(i,ksn)
        wtot_snow_new    (i,ksn) = rho_new(i,ksn)*dzh_snow_new(i,ksn)/rho_w
        wliq_snow_new    (i,ksn) = wl_new (i,ksn)*dzh_snow_new(i,ksn)
      END IF
    END DO
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF(w_snow_new(i) > eps_soil) THEN
      rho_snow_new(i) = w_snow_new(i)/h_snow_new(i)*rho_w
    ELSE !JH
      rho_snow_new(i) = 250._wp ! workaround need to be inspected!!
    END IF
    IF(w_snow_new(i) > eps_soil) THEN
      ! linear extrapolation from t_snow_mult_new(i,2) and t_snow_mult_new(i,1) to t_snow_mult_new(i,0)
      t_snow_mult_new(i,0) = (t_snow_mult_new(i,1)*(2._wp*dzh_snow_new(i,1)+dzh_snow_new(i,2))- &
                            t_snow_mult_new(i,2)*dzh_snow_new(i,1))/(dzh_snow_new(i,1)+dzh_snow_new(i,2))
      ! limiter to prevent unphysical values and/or numerical instabilities
      t_snow_mult_new(i,0) = MIN(273.15_wp,t_snow_mult_new(i,0),t_snow_mult_new(i,1)+5.0_wp)
      t_snow_mult_new(i,0) = MAX(t_snow_mult_new(i,0),t_snow_mult_new(i,1)-5.0_wp)
    END IF
    t_snow_new(i) = t_snow_mult_new (i,0)
  END DO
  !$ACC END PARALLEL

  !$ACC END DATA

END SUBROUTINE snow_multi_update_new_state


END MODULE sfc_terra_snow
