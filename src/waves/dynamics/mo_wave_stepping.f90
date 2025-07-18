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

! Initializes and controls the time stepping in the wave model.

MODULE mo_wave_stepping
  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_impl_constants,           ONLY: SUCCESS
  USE mo_run_config,               ONLY: output_mode, ltestcase, ltransport, msg_level
  USE mo_name_list_output,         ONLY: write_name_list_output, istime4name_list_output, istime4name_list_output_dom
  USE mo_name_list_output_init,    ONLY: output_file
  USE mo_output_event_handler,     ONLY: get_current_jfile
  USE mo_parallel_config,          ONLY: proc0_offloading
  USE mo_time_config,              ONLY: t_time_config
  USE mo_runtime_diag,             ONLY: print_timestep_info, print_wave_stats
  USE mtime,                       ONLY: datetime, timedelta, &
       &                                 OPERATOR(+), OPERATOR(>=), OPERATOR(==)
  USE mo_util_mtime,               ONLY: is_event_active
  USE mo_model_domain,             ONLY: p_patch
  USE mo_grid_config,              ONLY: n_dom, nroot
  USE mo_io_units,                 ONLY: filename_max
  USE mo_master_config,            ONLY: isRestart, isInitFromRestart, getModelBaseDir
  USE mo_dynamics_config,          ONLY: nnow, nnew
  USE mo_fortran_tools,            ONLY: swap, copy
  USE mo_intp_data_strc,           ONLY: p_int_state
  USE mo_pp_scheduler,             ONLY: new_simulation_status, pp_scheduler_process
  USE mo_pp_tasks,                 ONLY: t_simulation_status
  USE mo_wave_adv_exp,             ONLY: init_wind_adv_test
  USE mo_init_wave_physics,        ONLY: init_wave_spectrum, init_wave_nonlinear, fetch_law, jonswap, min_energy
  USE mo_wave_state,               ONLY: p_wave_state
  USE mo_wave_ext_data_state,      ONLY: wave_ext_data
  USE mo_wave_forcing_state,       ONLY: wave_forcing_state
  USE mo_wave_diagnostics,         ONLY: calculate_output_diagnostics
  USE mo_wave_source,              ONLY: src_wind_input, src_dissipation, src_bottom_friction, &
    &                                    src_nonlinear_transfer, integrate_in_time_src, &
    &                                    src_wave_breaking
  USE mo_wave_physics,             ONLY: tm1_tm2_periods_and_wm1_wm2_wavenumber, &
       &                                 mean_frequency_and_total_energy, air_sea, last_prog_freq_ind, &
       &                                 impose_high_freq_tail, wave_stress, &
       &                                 mask_energy, compute_wave_number, compute_group_velocity, sdepth_lim
  USE mo_wave_config,              ONLY: wave_config, generate_filename
  USE mo_energy_propagation_config,ONLY: energy_propagation_config
  USE mo_wave_forcing_state,       ONLY: wave_forcing_state
  USE mo_wave_forcing,             ONLY: t_read_wave_forcing
  USE mo_wave_events,              ONLY: waveCheckpointEvent, waveRestartEvent
  USE mo_wave_td_update,           ONLY: update_speed_and_direction, update_ice_free_mask, &
    &                                    update_water_depth_and_grad
  USE mo_wave_advection_stepping,  ONLY: wave_step_advection
  USE mo_coupling_config,          ONLY: is_coupled_to_atmo
  USE mo_timer,                    ONLY: ltimer, timer_start, timer_stop, timer_coupling, timers_level
  USE mo_wave_timer,               ONLY: timer_wave_total, timer_wave_reader, timer_wave_time_integration, &
    &                                    timer_wave_src, timer_wave_src_wind_input, &
    &                                    timer_wave_src_dissipation, timer_wave_src_nonlinear, &
    &                                    timer_wave_diagnostics
  USE mo_wave_atmo_coupling,       ONLY: couple_wave_to_atmo
  ! restart
  USE mo_restart,                  ONLY: t_RestartDescriptor
  USE mo_restart_nml_and_att,      ONLY: getAttributesForRestarting
  USE mo_key_value_store,          ONLY: t_key_value_store


  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_stepping'

  PUBLIC :: perform_wave_stepping

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Organizes wave time stepping
  !!
  SUBROUTINE perform_wave_stepping (time_config, restartDescriptor)

    CHARACTER(len=*), PARAMETER :: routine = modname//':perform_wave_stepping'

    TYPE(t_time_config),                 INTENT(IN) :: time_config  !< information for time control
    CLASS(t_RestartDescriptor), POINTER, INTENT(IN) :: restartDescriptor

    TYPE(datetime),  POINTER :: mtime_current     => NULL() !< current datetime
    TYPE(timedelta), POINTER :: model_time_step   => NULL()

    !
    ! note that the following TARGET attribute is essential! Otherwise the pointer to the
    ! specific reader inside the time interpolator object (this%reader in time_intp_intp)
    ! will lose its association status.
    TYPE(t_read_wave_forcing), ALLOCATABLE, TARGET :: reader_wave_forcing(:)
    INTEGER                  :: jstep                       !< time step number
    INTEGER                  :: jstep_shift                 !< number of time steps for backward shifting
    LOGICAL                  :: lprint_timestep             !< print current datetime information
    LOGICAL                  :: lprint_wave_stats           !< print wave height information
    INTEGER                  :: jg, jlev
    INTEGER                  :: ierrstat
    REAL(wp)                 :: dtime                       !< model time step in seconds
    TYPE(t_simulation_status):: simulation_status
    LOGICAL                  :: l_nml_output                !< TRUE, if output is due at current timestep
    LOGICAL                  :: lwrite_checkpoint

    ! Time levels
    INTEGER :: n_new, n_now

    CHARACTER(LEN=filename_max) :: wave_forc_wind_fn(n_dom) ! forc_file_prefix+'_wind' for U and V 10 meter wind (m/s)
    CHARACTER(LEN=filename_max) :: wave_forc_ice_fn(n_dom)  ! forc_file_prefix+'_ice'  for sea ice concentration (fraction of 1)
    CHARACTER(LEN=filename_max) :: wave_forc_slh_fn(n_dom)  ! forc_file_prefix+'_slh'  for sea level height (m)
    CHARACTER(LEN=filename_max) :: wave_forc_osc_fn(n_dom)  ! forc_file_prefix+'_osc'  for U and V ocean surface currents (m/s)
    ! Restarting
    TYPE(t_key_value_store), POINTER :: restartAttributes
    INTEGER, ALLOCATABLE :: output_jfile(:)
    LOGICAL :: l_isStartdate, l_isExpStopdate, l_isRestart, l_isCheckpoint, l_doWriteRestart
    INTEGER :: i

    IF (timers_level >= 1) CALL timer_start(timer_wave_total)

    lprint_wave_stats = msg_level > 9

    ! convenience pointer
    mtime_current   => time_config%tc_current_date  ! current datetime
    model_time_step => time_config%tc_dt_model      ! model time step

    ! allocate temporary variable for restarting purposes
    IF (output_mode%l_nml) THEN
      ALLOCATE(output_jfile(SIZE(output_file)), STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed for output_jfile!')
    ENDIF

    IF (ltestcase) THEN
      !-----------------------------------------------------------------------
      ! advection experiment
      CALL message(routine,'test case run: advection experiment')

      DO jg = 1, n_dom
        ! Initialisation of 10 meter wind and sea ice
        CALL init_wind_adv_test(p_patch(jg), wave_config(jg), wave_forcing_state(jg))
        CALL update_ice_free_mask(p_patch    = p_patch(jg),                          & ! IN
          &                    sea_ice_c     = wave_forcing_state(jg)%sea_ice_c,     & ! IN
          &                    ice_free_mask = wave_forcing_state(jg)%ice_free_mask_c) ! OUT
      END DO
    ENDIF


    IF (is_coupled_to_atmo()) THEN
      CALL message(routine,'coupled run: forcing data are received from the atmo model...')
    ELSE

      IF (timers_level >= 5) CALL timer_start(timer_wave_reader)

      CALL message(routine,'standalone run: forcing data are read from file...')

      ALLOCATE(reader_wave_forcing(n_dom), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, 'Allocation failed for reader_wave_forcing')

      DO jg = 1, n_dom
        IF (wave_config(jg)%lread_forcing) THEN

          jlev = p_patch(jg)%level

          wave_forc_wind_fn(jg) = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_wind.nc",&
            &                 getModelBaseDir(), nroot, jlev, jg)
          wave_forc_ice_fn(jg)  = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_ice.nc", &
            &                 getModelBaseDir(), nroot, jlev, jg)
          wave_forc_slh_fn(jg)  = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_slh.nc", &
            &                 getModelBaseDir(), nroot, jlev, jg)
          wave_forc_osc_fn(jg)  = generate_filename(TRIM(wave_config(jg)%forc_file_prefix)//"_osc.nc", &
            &                 getModelBaseDir(), nroot, jlev, jg)

          ! initialize reader of external forcing data
          CALL reader_wave_forcing(jg)%init(p_patch             = p_patch(jg),           & !in
            &                               destination_time    = mtime_current,         & !in
            &                               wave_forc_wind_file = wave_forc_wind_fn(jg), & !in
            &                               wave_forc_ice_file  = wave_forc_ice_fn(jg),  & !in
            &                               wave_forc_slh_file  = wave_forc_slh_fn(jg),  & !in
            &                               wave_forc_osc_file  = wave_forc_osc_fn(jg) )   !in

          ! get initial forcing data set (read from file and copy to forcing state vector)
          CALL reader_wave_forcing(jg)%update_forcing(                                &
            &                destination_time = mtime_current,                        & !in
            &                u10m             = wave_forcing_state(jg)%u10m,          & !out
            &                v10m             = wave_forcing_state(jg)%v10m,          & !out
            &                sp10m            = wave_forcing_state(jg)%sp10m,         & !out
            &                dir10m           = wave_forcing_state(jg)%dir10m,        & !out
            &                sic              = wave_forcing_state(jg)%sea_ice_c,     & !out
            &                slh              = wave_forcing_state(jg)%sea_level_c,   & !out
            &                uosc             = wave_forcing_state(jg)%usoce_c,       & !out
            &                vosc             = wave_forcing_state(jg)%vsoce_c,       & !out
            &                sp_osc           = wave_forcing_state(jg)%sp_soce_c,     & !out
            &                dir_osc          = wave_forcing_state(jg)%dir_soce_c,    & !out
            &                ice_free_mask_c  = wave_forcing_state(jg)%ice_free_mask_c) !out

          ! update depth and gradient
          CALL update_water_depth_and_grad(p_patch = p_patch(jg),                        & !in
            &                     p_int_state      = p_int_state(jg),                    & !in
            &                     bathymetry_c     = wave_ext_data(jg)%bathymetry_c,     & !in
            &                     sea_level_c      = wave_forcing_state(jg)%sea_level_c, & !in
            &                     depth_c          = wave_ext_data(jg)%depth_c,          & !out
            &                     depth_e          = wave_ext_data(jg)%depth_e,          & !out
            &                     geo_depth_grad_c = wave_ext_data(jg)%geo_depth_grad_c)   !out
        ELSE

          WRITE(message_text,'(a,a,a)') 'No forcing files specified, testcase run is assumed.'
          CALL message(routine, message_text)

        END IF
      END DO

      IF (timers_level >= 5) CALL timer_stop(timer_wave_reader)

    END IF  ! is_coupled_to_atmo


    DO jg = 1, n_dom
      n_now  = nnow(jg)
      n_new  = nnew(jg)

      ! Calculate the minimum values of energy allowed
      ! for each frequency for a given wind speed bin
      ! from 1 to wave_config%jmax, and up to wave_config%umax
      CALL min_energy(wave_config(jg), p_wave_state(jg)%diag%flminfr_tab)

      IF (isRestart() .OR. isInitFromRestart()) THEN
        ! do nothing
      ELSE  ! coldstart

        ! Initialisation of the wave spectrum
        CALL fetch_law(                                       &
          &  p_patch     = p_patch(jg),                       & !in
          &  fetch       = wave_config(jg)%fetch,             & !in
          &  fpmax       = wave_config(jg)%fm,                & !in
          &  sp10m       = wave_forcing_state(jg)%sp10m(:,:), & !in
          &  fp          = p_wave_state(jg)%diag%fp(:,:),     & !out
          &  alphaj      = p_wave_state(jg)%diag%alphaj(:,:))   !out

        ! Initialisation of the wave spectrum
        CALL init_wave_spectrum(                                          &
          &  p_patch     = p_patch(jg),                                   & !in
          &  wave_config = wave_config(jg),                               & !in
          &  dir10m      = wave_forcing_state(jg)%dir10m(:,:),            & !in
          &  fp          = p_wave_state(jg)%diag%fp(:,:),                 & !in
          &  alphaj      = p_wave_state(jg)%diag%alphaj(:,:),             & !in
          &  et          = p_wave_state(jg)%diag%et(:,:,:),               & !out  ! purely diagnostic
          &  tracer      = p_wave_state(jg)%prog(n_now)%tracer(:,:,:,:))    !out
      END IF


      ! initialisation of the nonlinear transfer computations
      ! computes time-constant index arrays and weights
      CALL init_wave_nonlinear(wave_config = wave_config(jg),     & !in
        &                      p_diag      = p_wave_state(jg)%diag) !inout


      ! compute wave number at centers and edges
      CALL compute_wave_number(                            &
        &  p_patch     = p_patch(jg),                      & !in
        &  wave_config = wave_config(jg),                  & !in
        &  depth_c     = wave_ext_data(jg)%depth_c,        & !in
        &  depth_e     = wave_ext_data(jg)%depth_e,        & !in
        &  wave_num_c  = p_wave_state(jg)%diag%wave_num_c, & !out
        &  wave_num_e  = p_wave_state(jg)%diag%wave_num_e  ) !out


      ! compute group velocity at centers and edges
      CALL compute_group_velocity(                          &
        &  p_patch      = p_patch(jg),                      & !in
        &  wave_config  = wave_config(jg),                  & !in
        &  wave_num_c   = p_wave_state(jg)%diag%wave_num_c, & !in
        &  wave_num_e   = p_wave_state(jg)%diag%wave_num_e, & !in
        &  depth_c      = wave_ext_data(jg)%depth_c,        & !in
        &  depth_e      = wave_ext_data(jg)%depth_e,        & !in
        &  gv_c         = p_wave_state(jg)%diag%gv_c,       & !out
        &  gv_e         = p_wave_state(jg)%diag%gv_e)         !out

    END DO


    IF (isRestart()) THEN
      CALL getAttributesForRestarting(restartAttributes)
      ! get start counter for time loop from restart file:
      CALL restartAttributes%get("jstep", jstep)

    ELSE  ! no restart, or isInitFromRestart

      ! initialize time step counter
      !
      IF (time_config%timeshift%dt_shift < 0._wp) THEN
        ! get model time step in seconds for base domain
        dtime = time_config%get_model_timestep_sec(p_patch(1)%nest_level)
        jstep_shift = NINT(time_config%timeshift%dt_shift/dtime)
        jstep = jstep_shift
        !
        WRITE(message_text,'(a,i6,a)') 'Model start shifted backwards by ', ABS(jstep_shift),' time steps'
        CALL message(routine, message_text)
      ELSE
        jstep = 0
      ENDIF

      DO jg = 1,n_dom
        IF (.NOT. p_patch(jg)%ldom_active) CYCLE

        ! Calculate total and mean frequency energy
        CALL mean_frequency_and_total_energy(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%prog(n_now)%tracer, &
             p_wave_state(jg)%source%llws, &
             p_wave_state(jg)%diag%emean, & ! OUT
             p_wave_state(jg)%diag%emeanws, & ! OUT
             p_wave_state(jg)%diag%femean, & ! OUT
             p_wave_state(jg)%diag%femeanws) ! OUT

        ! Calculate roughness length and friction velocities
        CALL air_sea(p_patch(jg), wave_config(jg), &
             wave_forcing_state(jg)%sp10m, &
             p_wave_state(jg)%diag%tauw, &
             p_wave_state(jg)%diag%ustar, & ! OUT
             p_wave_state(jg)%diag%z0)      ! OUT

        ! Calculate tm1 period and f1 frequency and wavenumbers
        CALL tm1_tm2_periods_and_wm1_wm2_wavenumber(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%diag%wave_num_c, &
             p_wave_state(jg)%prog(n_now)%tracer, &
             p_wave_state(jg)%diag%emean, &
             p_wave_state(jg)%diag%tm1, &  ! OUT
             p_wave_state(jg)%diag%tm2, &  ! OUT
             p_wave_state(jg)%diag%f1mean, & !OUT
             p_wave_state(jg)%diag%akmean, & !OUT
             p_wave_state(jg)%diag%xkmean) !OUT


        IF (istime4name_list_output_dom(jg=jg, jstep=jstep)) THEN

          ! Calculation of diagnostic output parameters
          CALL calculate_output_diagnostics(p_patch = p_patch(jg),                    & ! IN
            &                      wave_config = wave_config(jg),                     & ! IN
            &                            sp10m = wave_forcing_state(jg)%sp10m,        & ! IN
            &                           dir10m = wave_forcing_state(jg)%dir10m,       & ! IN
            &                           depth  = wave_ext_data(jg)%depth_c,           & ! IN
            &                           tracer = p_wave_state(jg)%prog(n_now)%tracer, & ! IN
            &                           p_diag = p_wave_state(jg)%diag)                 ! INOUT
        ENDIF
      ENDDO

      !--------------------------------------------------------------------------
      ! loop over the list of internal post-processing tasks, e.g.
      ! interpolate selected fields to lat-lon
      simulation_status = new_simulation_status(l_first_step   = .TRUE.,                  &
        &                                       l_output_step  = .TRUE.,                  &
        &                                       l_dom_active   = p_patch(1:)%ldom_active, &
        &                                       i_timelevel_dyn= nnow, i_timelevel_phy= nnow)
      CALL pp_scheduler_process(simulation_status, lacc=.TRUE.)

      ! output at initial time
      !
      ! Ensure that the initial output is skipped for model runs with shifted start date
      IF (output_mode%l_nml .AND. (mtime_current >= time_config%tc_exp_startdate) ) THEN
        CALL write_name_list_output(jstep=0)
      END IF

    ENDIF  ! isRestart


    TIME_LOOP: DO

      ! update model date and time
      mtime_current = mtime_current + model_time_step
      jstep = jstep + 1

      ! store state of output files for restarting purposes
      IF (output_mode%l_nml .AND. jstep>=0 ) THEN
        DO i=1,SIZE(output_file)
          output_jfile(i) = get_current_jfile(output_file(i)%out_event)
        END DO
      ENDIF


      lprint_timestep   = msg_level > 2 .OR. MOD(jstep,25) == 0
      !
      IF (lprint_timestep) THEN
        CALL print_timestep_info(time_config, jstep)
      ENDIF


      DO jg = 1, n_dom

        IF (.NOT. p_patch(jg)%ldom_active) CYCLE

        n_now  = nnow(jg)
        n_new  = nnew(jg)

        IF (is_coupled_to_atmo()) THEN
          ! send and receive coupling fields
          !
          IF (ltimer) CALL timer_start(timer_coupling)
          CALL couple_wave_to_atmo(p_patch   = p_patch(jg),                     & ! IN
            &                      z0        = p_wave_state(jg)%diag%z0,        & ! IN
            &                      u10m      = wave_forcing_state(jg)%u10m,     & ! OUT
            &                      v10m      = wave_forcing_state(jg)%v10m,     & ! OUT
            &                      sea_ice_c = wave_forcing_state(jg)%sea_ice_c ) ! OUT
          IF (ltimer) CALL timer_stop(timer_coupling)

          ! update forcing state
          ! update wind speed and direction
          CALL update_speed_and_direction(p_patch = p_patch(jg),                   & ! IN
            &                               u     = wave_forcing_state(jg)%u10m,   & ! IN
            &                               v     = wave_forcing_state(jg)%v10m,   & ! IN
            &                              sp     = wave_forcing_state(jg)%sp10m,  & ! OUT
            &                              dir    = wave_forcing_state(jg)%dir10m)   ! OUT

          ! update ice-free mask
          CALL update_ice_free_mask(p_patch    = p_patch(jg),                          & ! IN
            &                    sea_ice_c     = wave_forcing_state(jg)%sea_ice_c,     & ! IN
            &                    ice_free_mask = wave_forcing_state(jg)%ice_free_mask_c) ! OUT

        ELSE
          ! get new forcing data (read from file and copy to forcing state vector)
          IF (wave_config(jg)%lread_forcing) THEN
            CALL reader_wave_forcing(jg)%update_forcing(                                &
              &                destination_time = mtime_current,                        & !in
              &                u10m             = wave_forcing_state(jg)%u10m,          & !out
              &                v10m             = wave_forcing_state(jg)%v10m,          & !out
              &                sp10m            = wave_forcing_state(jg)%sp10m,         & !out
              &                dir10m           = wave_forcing_state(jg)%dir10m,        & !out
              &                sic              = wave_forcing_state(jg)%sea_ice_c,     & !out
              &                slh              = wave_forcing_state(jg)%sea_level_c,   & !out
              &                uosc             = wave_forcing_state(jg)%usoce_c,       & !out
              &                vosc             = wave_forcing_state(jg)%vsoce_c,       & !out
              &                sp_osc           = wave_forcing_state(jg)%sp_soce_c,     & !out
              &                dir_osc          = wave_forcing_state(jg)%dir_soce_c,    & !out
              &                ice_free_mask_c  = wave_forcing_state(jg)%ice_free_mask_c) !out

            ! update depth and gradient
            CALL update_water_depth_and_grad(p_patch = p_patch(jg),                        & !in
              &                     p_int_state      = p_int_state(jg),                    & !in
              &                     bathymetry_c     = wave_ext_data(jg)%bathymetry_c,     & !in
              &                     sea_level_c      = wave_forcing_state(jg)%sea_level_c, & !in
              &                     depth_c          = wave_ext_data(jg)%depth_c,          & !out
              &                     depth_e          = wave_ext_data(jg)%depth_e,          & !out
              &                     geo_depth_grad_c = wave_ext_data(jg)%geo_depth_grad_c)   !out

          END IF
        END IF ! is_coupled_to_atmo()

        ! horizontal propagation of binned wave energy
        ! Here, we integrate the spectral energy equation in time without sources and sinks,
        ! only taking into account advection and refraction.
        ! If the horizontal propagation is deactivated, a simple copy is performed from
        ! prog(n_now)%tracer to prog(n_new)%tracer
        !
        IF (ltransport) THEN
          ! get model time step in seconds
          dtime = time_config%get_model_timestep_sec(p_patch(jg)%nest_level)
          !
          CALL wave_step_advection(p_patch                   = p_patch(jg),                           & !in
            &                      p_int_state               = p_int_state(jg),                       & !in
            &                      wave_config               = wave_config(jg),                       & !in
            &                      energy_propagation_config = energy_propagation_config(jg),         & !in
            &                      p_dtime                   = dtime,                                 & !in
            &                      wave_num_c                = p_wave_state(jg)%diag%wave_num_c,      & !in
            &                      gv_c                      = p_wave_state(jg)%diag%gv_c,            & !in
            &                      gv_e                      = p_wave_state(jg)%diag%gv_e,            & !in
            &                      depth_c                   = wave_ext_data(jg)%depth_c,             & !in
            &                      geo_depth_grad_c          = wave_ext_data(jg)%geo_depth_grad_c,    & !in
            &                      p_tracer_now              = p_wave_state(jg)%prog(n_now)%tracer,   & !in
            &                      p_tracer_new              = p_wave_state(jg)%prog(n_new)%tracer    ) !out
        ELSE
!$OMP PARALLEL
          CALL copy(src  = p_wave_state(jg)%prog(n_now)%tracer, &
            &       dest = p_wave_state(jg)%prog(n_new)%tracer, lacc=.FALSE.)
!$OMP END PARALLEL
        ENDIF

        IF (timers_level >= 5) CALL timer_start(timer_wave_src)
        !
        IF (timers_level >= 8) CALL timer_start(timer_wave_src_wind_input)
        ! Calculate total and mean frequency energy
        CALL mean_frequency_and_total_energy(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%prog(n_new)%tracer, &
             p_wave_state(jg)%source%llws,&
             p_wave_state(jg)%diag%emean, & ! OUT
             p_wave_state(jg)%diag%emeanws, & ! OUT
             p_wave_state(jg)%diag%femean, & ! OUT
             p_wave_state(jg)%diag%femeanws) ! OUT

        ! The advection of wave energy is not directly limited by depth and has no control over
        ! whether the amount of energy transported may lead to an unphysically high level
        ! for a given depth (only transit from deep to shallow regions are affected).
        ! Therefore, the level of energy is cutted by a depth limiter just after transport.
        ! This allows the subsequent calls of source function terms to start from
        ! the 'correct' level of energy.
        ! Reduce wave energy if larger than depth limited energy level
        CALL sdepth_lim(p_patch     = p_patch(jg),                       & ! IN
                        wave_config = wave_config(jg),                   & ! IN
                        depth       = wave_ext_data(jg)%depth_c,         & ! IN
                        emean       = p_wave_state(jg)%diag%emean,       & ! INOUT
                        tracer      = p_wave_state(jg)%prog(n_new)%tracer) ! INOUT

        ! Calculate tm1 period and f1 frequency and wavenumbers
        CALL tm1_tm2_periods_and_wm1_wm2_wavenumber(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%diag%wave_num_c, &
             p_wave_state(jg)%prog(n_new)%tracer, &
             p_wave_state(jg)%diag%emean, &
             p_wave_state(jg)%diag%tm1, &  ! OUT
             p_wave_state(jg)%diag%tm2, &  ! OUT
             p_wave_state(jg)%diag%f1mean, & !OUT
             p_wave_state(jg)%diag%akmean, & !OUT
             p_wave_state(jg)%diag%xkmean) !OUT

        ! Calculate roughness length and friction velocities
        CALL air_sea(p_patch(jg), wave_config(jg), &
             wave_forcing_state(jg)%sp10m, &
             p_wave_state(jg)%diag%tauw, &
             p_wave_state(jg)%diag%ustar, & ! OUT
             p_wave_state(jg)%diag%z0)      ! OUT

        ! Calculate wind input source function
        IF (wave_config(jg)%linput_sf1) THEN
          CALL src_wind_input(                                    &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  dir10m      = wave_forcing_state(jg)%dir10m,       & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_diag      = p_wave_state(jg)%diag,               & !in: ustar,z0,wave_num_c
            &  p_source    = p_wave_state(jg)%source)               !inout: llws,fl,sl
        END IF

        ! Update total and mean frequency energy
        CALL mean_frequency_and_total_energy(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%prog(n_new)%tracer, &
             p_wave_state(jg)%source%llws,&
             p_wave_state(jg)%diag%emean, & ! OUT
             p_wave_state(jg)%diag%emeanws, & ! OUT
             p_wave_state(jg)%diag%femean, & ! OUT
             p_wave_state(jg)%diag%femeanws) ! OUT

        ! Calculate last frequency index of prognostic part of spectrum
        CALL last_prog_freq_ind(                           &
          &  p_patch     = p_patch(jg),                    & !IN
          &  wave_config = wave_config(jg),                & !IN
          &  femeanws    = p_wave_state(jg)%diag%femeanws, & !IN
          &  femean      = p_wave_state(jg)%diag%femean,   & !IN
          &  ustar       = p_wave_state(jg)%diag%ustar,    & !IN
          &  lpfi        = p_wave_state(jg)%diag%last_prog_freq_ind) !OUT

        ! Calculate wave stress
        IF (wave_config(jg)%lwave_stress1) THEN
          CALL wave_stress(                                       &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  dir10m      = wave_forcing_state(jg)%dir10m,       & !in
            &  sl          = p_wave_state(jg)%source%sl,          & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_diag      = p_wave_state(jg)%diag                ) !IN : last_prog_freq_ind,ustar,z0
                                                                    !OUT: phiaw,tauw,tauhf,phihf
        END IF

        ! Update roughness length and friction velocities
        CALL air_sea(p_patch(jg), wave_config(jg), &
             wave_forcing_state(jg)%sp10m, &
             p_wave_state(jg)%diag%tauw, &
             p_wave_state(jg)%diag%ustar, & ! OUT
             p_wave_state(jg)%diag%z0)      ! OUT

        ! Impose high frequency tail to the spectrum
        CALL impose_high_freq_tail(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%diag%wave_num_c,         & !IN
             wave_ext_data(jg)%depth_c,                & !IN
             p_wave_state(jg)%diag%last_prog_freq_ind, & !IN
             p_wave_state(jg)%prog(n_new)%tracer)        !INOUT

        ! Recompute wind input source function
        IF (wave_config(jg)%linput_sf2) THEN
          CALL src_wind_input(                                    &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  dir10m      = wave_forcing_state(jg)%dir10m,       & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_diag      = p_wave_state(jg)%diag,               & !in: ustar,z0,wave_num_c
            &  p_source    = p_wave_state(jg)%source)               !inout: llws,fl,sl
        END IF

        ! Update wave stress
        IF (wave_config(jg)%lwave_stress2) THEN
          CALL wave_stress(                                       &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  dir10m      = wave_forcing_state(jg)%dir10m,       & !in
            &  sl          = p_wave_state(jg)%source%sl,          & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_diag      = p_wave_state(jg)%diag                ) !IN : last_prog_freq_ind,ustar,z0
                                                                    !OUT: phiaw,tauw,tauhf,phihf
        END IF
        IF (timers_level >= 8) CALL timer_stop(timer_wave_src_wind_input)

        ! Calculate dissipation source function
        IF (wave_config(jg)%ldissip_sf) THEN
          IF (timers_level >= 8) CALL timer_start(timer_wave_src_dissipation)

          CALL src_dissipation(                                   &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  wave_num_c  = p_wave_state(jg)%diag%wave_num_c,    & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_diag      = p_wave_state(jg)%diag,               & !in: f1mean,emean,xkmean
            &  p_source    = p_wave_state(jg)%source)               !inout: fl,sl

          IF (timers_level >= 8) CALL timer_stop(timer_wave_src_dissipation)
        END IF

        ! Calculate source function due to nonlinear transfer
        IF (wave_config(jg)%lnon_linear_sf) THEN
          IF (timers_level >= 8) CALL timer_start(timer_wave_src_nonlinear)

          CALL src_nonlinear_transfer(                            &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  depth       = wave_ext_data(jg)%depth_c,           & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_diag      = p_wave_state(jg)%diag,               & !in
            &  p_source    = p_wave_state(jg)%source)               !inout: fl,sl

          IF (timers_level >= 8) CALL timer_stop(timer_wave_src_nonlinear)
        END IF

        IF (timers_level >= 8) CALL timer_start(timer_wave_src_dissipation)
        ! Calculate dissipation due to bottom friction
        IF (wave_config(jg)%lbottom_fric_sf) THEN
          CALL src_bottom_friction(                               &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  wave_num_c  = p_wave_state(jg)%diag%wave_num_c,    & !in
            &  depth       = wave_ext_data(jg)%depth_c,           & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_source    = p_wave_state(jg)%source)               !inout: fl, sl
        END IF

        ! Calculate dissipation due to depth-induced wave breaking
        IF (wave_config(jg)%lwave_brk_sf) THEN
          CALL src_wave_breaking(                                 &
            &  p_patch     = p_patch(jg),                         & !in
            &  wave_config = wave_config(jg),                     & !in
            &  depth_c     = wave_ext_data(jg)%depth_c,           & !in
            &  tracer      = p_wave_state(jg)%prog(n_new)%tracer, & !in
            &  p_diag      = p_wave_state(jg)%diag,               & !inout, in: emean, f1mean out: hrms_frac, wbr_frac
            &  p_source    = p_wave_state(jg)%source)               !inout: fl, sl
        END IF
        IF (timers_level >= 8) CALL timer_stop(timer_wave_src_dissipation)
        !
        IF (timers_level >= 5) CALL timer_stop(timer_wave_src)

        IF (timers_level >= 5) CALL timer_start(timer_wave_time_integration)
        ! Calculate new spectrum
        CALL integrate_in_time_src(                           &
          &  p_patch     = p_patch(jg),                       & !in
          &  wave_config = wave_config(jg),                   & !in
          &  p_diag      = p_wave_state(jg)%diag,             & !in ustar, femeanws, femean
          &  p_source    = p_wave_state(jg)%source,           & !in sl, fl
          &  sp10m       = wave_forcing_state(jg)%sp10m,      & !in
          &  dir10m      = wave_forcing_state(jg)%dir10m,     & !in
          &  tracer      = p_wave_state(jg)%prog(n_new)%tracer) !inout

        ! Set energy to zero under the sea ice
        CALL mask_energy(p_patch(jg), wave_config(jg), &
             wave_forcing_state(jg)%ice_free_mask_c, & !IN
             p_wave_state(jg)%prog(n_new)%tracer) ! INOUT

        ! Update total and mean frequency energy
        CALL mean_frequency_and_total_energy(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%prog(n_new)%tracer, &
             p_wave_state(jg)%source%llws,&
             p_wave_state(jg)%diag%emean, & ! OUT
             p_wave_state(jg)%diag%emeanws, & ! OUT
             p_wave_state(jg)%diag%femean, & ! OUT
             p_wave_state(jg)%diag%femeanws) ! OUT

        ! Update high frequency tail
        CALL last_prog_freq_ind(                           &
          &  p_patch     = p_patch(jg),                    & !IN
          &  wave_config = wave_config(jg),                & !IN
          &  femeanws    = p_wave_state(jg)%diag%femeanws, & !IN
          &  femean      = p_wave_state(jg)%diag%femean,   & !IN
          &  ustar       = p_wave_state(jg)%diag%ustar,    & !IN
          &  lpfi        = p_wave_state(jg)%diag%last_prog_freq_ind) !OUT

        CALL impose_high_freq_tail(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%diag%wave_num_c,         & !IN
             wave_ext_data(jg)%depth_c,                & !IN
             p_wave_state(jg)%diag%last_prog_freq_ind, & !IN
             p_wave_state(jg)%prog(n_new)%tracer)        !INOUT

        ! Update total and mean frequency energy
        CALL mean_frequency_and_total_energy(p_patch(jg), wave_config(jg), &
             p_wave_state(jg)%prog(n_new)%tracer, &
             p_wave_state(jg)%source%llws,&
             p_wave_state(jg)%diag%emean, & ! OUT
             p_wave_state(jg)%diag%emeanws, & ! OUT
             p_wave_state(jg)%diag%femean, & ! OUT
             p_wave_state(jg)%diag%femeanws) ! OUT

        ! switch between time levels now and new for next time step
        CALL swap(nnow(jg), nnew(jg))
        !
        IF (timers_level >= 5) CALL timer_stop(timer_wave_time_integration)


        !--------------------------------------------------------------------------
        ! Output section
        !--------------------------------------------------------------------------

        IF (timers_level >= 5) CALL timer_start(timer_wave_diagnostics)

        IF (istime4name_list_output_dom(jg=jg, jstep=jstep)) THEN
          ! Calculation of diagnostic output parameters
          ! Calculation is performed only at output times
          CALL calculate_output_diagnostics(p_patch = p_patch(jg),                       & ! IN
            &                      wave_config = wave_config(jg),                        & ! IN
            &                            sp10m = wave_forcing_state(jg)%sp10m,           & ! IN
            &                           dir10m = wave_forcing_state(jg)%dir10m,          & ! IN
            &                           depth  = wave_ext_data(jg)%depth_c,              & ! IN
            &                           tracer = p_wave_state(jg)%prog(nnow(jg))%tracer, & ! IN
            &                           p_diag = p_wave_state(jg)%diag)                    ! INOUT
        ENDIF

        IF (lprint_wave_stats) THEN
          ! Print information on global maxima and minima to stdout
          CALL print_wave_stats(p_patch = p_patch(jg),                      & !IN
            &                   emean   = p_wave_state(jg)%diag%emean(:,:), & !IN
            &                   femean  = p_wave_state(jg)%diag%femean(:,:) ) !IN
        ENDIF

        IF (timers_level >= 5) CALL timer_stop(timer_wave_diagnostics)
      ENDDO ! jg

      l_nml_output = output_mode%l_nml .AND. jstep >= 0 .AND. istime4name_list_output(jstep)
      simulation_status = new_simulation_status(l_output_step  = l_nml_output,             &
        &                                       l_last_step    = (mtime_current >= time_config%tc_stopdate), &
        &                                       l_accumulation_step = .FALSE.,             &
        &                                       l_dom_active   = p_patch(1:)%ldom_active,  &
        &                                       i_timelevel_dyn= nnow,                     &
        &                                       i_timelevel_phy= nnow)
      CALL pp_scheduler_process(simulation_status, lacc=.TRUE.)


      IF (l_nml_output) THEN
        CALL write_name_list_output(jstep=jstep)
      END IF


      !--------------------------------------------------------------------------
      ! Write restart file
      !--------------------------------------------------------------------------
      ! check whether time has come for writing restart file

      ! default is to assume we do not write a checkpoint/restart file
      lwrite_checkpoint = .FALSE.

      ! if the model is not supposed to write output, do not write checkpoints
      IF (.NOT. output_mode%l_none) THEN
        ! to clarify the decision tree we use shorter and more expressive names:

        l_isStartdate    = (time_config%tc_startdate == mtime_current)
        l_isExpStopdate  = (time_config%tc_exp_stopdate == mtime_current)
        l_isRestart      = is_event_active(waveRestartEvent, mtime_current, proc0_offloading)
        l_isCheckpoint   = is_event_active(waveCheckpointEvent, mtime_current, proc0_offloading)
        l_doWriteRestart = time_config%tc_write_restart

        IF ( &
             !  if normal checkpoint or restart cycle has been reached, i.e. checkpoint+model stop
             &         (l_isRestart .OR. l_isCheckpoint)                     &
             &  .AND.                                                        &
             !  and the current date differs from the start date
             &        .NOT. l_isStartdate                                    &
             &  .AND.                                                        &
             !  and end of run has not been reached or restart writing has been disabled
             &        (.NOT. l_isExpStopdate .OR. l_doWriteRestart)          &
             & ) THEN
          lwrite_checkpoint = .TRUE.
        END IF
      END IF


      IF (lwrite_checkpoint) THEN
        DO jg = 1,n_dom
          CALL restartDescriptor%updatePatch(p_patch(jg), opt_ndom=n_dom)
        ENDDO

        ! trigger writing of restart files.
        CALL restartDescriptor%writeRestart(mtime_current, jstep, opt_output_jfile = output_jfile)
      END IF  ! lwrite_checkpoint


      IF (mtime_current >= time_config%tc_stopdate) THEN
        ! leave time loop
        EXIT TIME_LOOP
      END IF

    ENDDO TIME_LOOP


    ! cleanup
    IF (ALLOCATED(output_jfile)) THEN
      DEALLOCATE(output_jfile, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed for output_jfile!')
    ENDIF
    !
    IF (ALLOCATED(reader_wave_forcing)) THEN
      DO jg=1,n_dom

        CALL reader_wave_forcing(jg)%deinit()
        !
      ENDDO
      DEALLOCATE(reader_wave_forcing, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, 'Deallocation failed for reader_wave_forcing')
    ENDIF

    CALL message(routine,'finished')

    IF (timers_level >= 1) CALL timer_stop(timer_wave_total)

  END SUBROUTINE perform_wave_stepping

END MODULE mo_wave_stepping
