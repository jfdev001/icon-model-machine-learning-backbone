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
! wave-specifc runtime diagnostics which can be written to stdout
!
MODULE mo_runtime_diag

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_time_config,         ONLY: t_time_config
  USE mo_exception,           ONLY: message, message_text
  USE mo_util_mtime,          ONLY: mtime_utils, FMT_DDHHMMSS_DAYSEP
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: global_max
  USE mo_mpi,                 ONLY: process_mpi_stdio_id, p_comm_work, p_min

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: print_timestep_info
  PUBLIC :: print_wave_stats

CONTAINS

  !
  ! Write timestep information to stdout
  !
  SUBROUTINE print_timestep_info (tc, jstep)
    !
    TYPE(t_time_config), INTENT(IN) :: tc      !< time_config object
    INTEGER,             INTENT(IN) :: jstep   !< timestep counter

    CALL message('','')

    WRITE(message_text,'(a,i8,a,i0,a,5(i2.2,a),i3.3,a,a)') &
      &             'Time step waves: ', jstep, ', model time: ',                              &
      &             tc%tc_current_date%date%year,   '-', tc%tc_current_date%date%month,    '-',    &
      &             tc%tc_current_date%date%day,    ' ', tc%tc_current_date%time%hour,     ':',    &
      &             tc%tc_current_date%time%minute, ':', tc%tc_current_date%time%second,   '.',    &
      &             tc%tc_current_date%time%ms, ' forecast time ',                            &
      &             TRIM(mtime_utils%ddhhmmss(tc%tc_exp_startdate, &
      &                                       tc%tc_current_date, FMT_DDHHMMSS_DAYSEP))

    CALL message('',message_text)

  END SUBROUTINE print_timestep_info


  !
  ! Print detailed information on global maxima and minima to stdout
  ! for significant wave height hs and mean wave period tmp.
  !
  SUBROUTINE print_wave_stats (p_patch, emean, femean)
    !
    TYPE(t_patch),   INTENT(IN) :: p_patch
    REAL(wp),        INTENT(IN) :: emean(:,:)     !< total energy [m^2]
    REAL(wp),        INTENT(IN) :: femean(:,:)    !< mean frequency energy [m^2]

    ! local
    REAL(wp):: hs                                 !< significant wave height [m]
    REAL(wp):: tmp                                !< total mean wave period [s]
    REAL(wp):: hs_max, hs_min, tmp_max, tmp_min                    ! mpi process-specific max/min
    REAL(wp):: hs_max_glb, hs_min_glb, tmp_max_glb, tmp_min_glb    ! global max/min
    INTEGER :: i_rlstart, i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc

    i_rlstart  = 1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    hs_max  = -HUGE(1._wp)
    hs_min  =  HUGE(1._wp)
    tmp_max = -HUGE(1._wp)
    tmp_min =  HUGE(1._wp)
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,hs,tmp) REDUCTION(max:hs_max,tmp_max) &
!$OMP    REDUCTION(min:hs_min,tmp_min)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        hs  = 4.0_wp * SQRT(emean(jc,jb))  ! significant wave height [m]
        tmp = 1.0_wp / femean(jc,jb)       ! total mean wave period [s]
        hs_max  = MAX(hs_max, hs)
        hs_min  = MIN(hs_min, hs)
        tmp_max = MAX(tmp_max, tmp)
        tmp_min = MIN(tmp_min, tmp)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    hs_max_glb  = global_max(hs_max,  iroot=process_mpi_stdio_id)
    tmp_max_glb = global_max(tmp_max, iroot=process_mpi_stdio_id)
    hs_min_glb  = p_min(hs_min,  comm=p_comm_work, root=process_mpi_stdio_id)
    tmp_min_glb = p_min(tmp_min, comm=p_comm_work, root=process_mpi_stdio_id)


    WRITE(message_text,'(a,i3,a,e18.10,a,e18.10)') &
      &  'MIN/MAX hs  [m] in domain', p_patch%id, ':', hs_min_glb,' /',hs_max_glb
    CALL message('',message_text)
    WRITE(message_text,'(a,i3,a,e18.10,a,e18.10)') &
      &  'MIN/MAX tmp [s] in domain', p_patch%id, ':',  tmp_min_glb,' /',tmp_max_glb
    CALL message('',message_text)

  END SUBROUTINE print_wave_stats

END MODULE mo_runtime_diag
