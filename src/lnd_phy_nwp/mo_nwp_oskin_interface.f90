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
! This module is the interface between nwp_sfc_interface and the
! ocean skin parameterisation.  It calls the ocean skin parameterization
! from ECMWF: mo_voskin.f90.
!
! ---------------------------------------------------------------


!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_oskin_interface

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level
  USE mo_lnd_nwp_config,      ONLY: isub_water, itype_oskin_warm, itype_oskin_cold
  USE mo_voskin,              ONLY: voskin

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_oskin

CONTAINS

!-------------------------------------------------------------------------

  SUBROUTINE nwp_oskin  (dtime, p_patch, prm_diag, ext_data, lnd_prog_now, lnd_diag)

    REAL(wp),              INTENT(in)   :: dtime          !< time interval for surface
    TYPE(t_patch),         INTENT(in)   :: p_patch        !< grid/patch info
    TYPE(t_nwp_phy_diag),  INTENT(in)   :: prm_diag       !< atm phys vars
    TYPE(t_external_data), INTENT(in)   :: ext_data       !< external data
    TYPE(t_lnd_prog),      INTENT(in)   :: lnd_prog_now   !< prog vars for sfc
    TYPE(t_lnd_diag),      INTENT(inout):: lnd_diag       !< diag vars for sfc

    ! Local arrays  (local copies)
    !
    REAL(wp) :: shfl_s   (nproma)      ! sensible heat flux at the surface            [W/m^2]
    REAL(wp) :: lhfl_s   (nproma)      ! latent heat flux at the surface              [W/m^2]
    REAL(wp) :: lwflxsfc (nproma)      ! net long-wave radiation flux at the surface  [W/m^2]
    REAL(wp) :: swflxsfc (nproma)      ! net solar radiation flux at the surface      [W/m^2]
    REAL(wp) :: umfl_s   (nproma)      ! U momentum flux at sfc                       [N/m2]
    REAL(wp) :: vmfl_s   (nproma)      ! V momentum flux at sfc                       [N/m2]
    REAL(wp) :: u_10m_t  (nproma)      ! 10 U wind                                    [m/s]
    REAL(wp) :: v_10m_t  (nproma)      ! 10 V wind                                    [m/s]
    REAL(wp) :: t_g_old_t(nproma)      ! surface skin temperature old                 [K]
    REAL(wp) :: t_seasfc (nproma)      ! sea surface temperature (below warm layer)   [K]
    REAL(wp) :: sst_warm_layer(nproma) ! SST warm layer                               [K]
    REAL(wp) :: sst_cold_skin (nproma) ! SST cold skin                                [K]

    ! Local array bounds:
    !
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    ! blocks

    ! Local scalars:
    !
    INTEGER :: jc, jb, ic              ! loop indices
    INTEGER :: i_count

    CHARACTER(len=*), PARAMETER :: routine = 'mo_nwp_sfc_interface:nwp_oskin'

!-------------------------------------------------------------------------

    ! exclude nest boundary and halo points
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    IF (msg_level >= 15) THEN
      CALL message(routine, 'call nwp_oskin scheme')
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_count,ic,jc,shfl_s,lhfl_s,lwflxsfc,swflxsfc,        &
!$OMP            umfl_s,vmfl_s,u_10m_t,v_10m_t,t_g_old_t,t_seasfc,        &
!$OMP            sst_warm_layer,sst_cold_skin) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      ! Copy input fields

      i_count = ext_data%atm%list_seawtr%ncount(jb)

      IF (i_count == 0) CYCLE ! skip loop if the index list for the given block is empty

      DO ic = 1, i_count

        jc = ext_data%atm%list_seawtr%idx(ic,jb)

        umfl_s   (ic) = prm_diag%umfl_s_t  (jc,jb,isub_water)   ! U momentum flux at sfc          [N/m2]
        vmfl_s   (ic) = prm_diag%vmfl_s_t  (jc,jb,isub_water)   ! V momentum flux at sfc          [N/m2]
        u_10m_t  (ic) = prm_diag%u_10m_t   (jc,jb,isub_water)
        v_10m_t  (ic) = prm_diag%v_10m_t   (jc,jb,isub_water)
        shfl_s   (ic) = prm_diag%shfl_s_t  (jc,jb,isub_water)   ! sensible heat flux at sfc       [W/m^2]
        lhfl_s   (ic) = prm_diag%lhfl_s_t  (jc,jb,isub_water)   ! latent heat flux at sfc         [W/m^2]
        lwflxsfc (ic) = prm_diag%lwflxsfc_t(jc,jb,isub_water)   ! net lw radiation flux at sfc    [W/m^2]
        swflxsfc (ic) = prm_diag%swflxsfc_t(jc,jb,isub_water)   ! net solar radiation flux at sfc [W/m^2]
        t_g_old_t(ic) = lnd_prog_now%t_g_t (jc,jb,isub_water)
        t_seasfc (ic) = lnd_diag%t_seasfc  (jc,jb)              ! SST foundation temperature      [K]
        IF (itype_oskin_warm > 0) THEN
          sst_warm_layer(ic) = lnd_diag%sst_warm_layer(jc,jb)   ! warm layer (prognostic)         [K]
        ELSE
          sst_warm_layer(ic) = 0._wp
        ENDIF
      ENDDO  ! ic

      ! call ocean skin and warm layer scheme from ECMWF

      CALL voskin (                            &
                 & kidia  = 1                , & !start index
                 & kfdia  = i_count          , & !end index
                 & klon   = i_count          , & !length of index
                 & ptmst  = dtime            , & !in
                 & pssrfl = swflxsfc (:)     , & !in
                 & pslrfl = lwflxsfc (:)     , & !in
                 & pahfs  = shfl_s   (:)     , & !in
                 & pahfl  = lhfl_s   (:)     , & !in
                 & pustr  = umfl_s   (:)     , & !in
                 & pvstr  = vmfl_s   (:)     , & !in
                 & pu10   = u_10m_t  (:)     , & !in
                 & pv10   = v_10m_t  (:)     , & !in
                 & ptskm1m= t_g_old_t(:)     , & !in
                 & psst   = t_seasfc (:)     , & !in
                 & pdwarm = sst_warm_layer(:), & !inout
                 & pdcool = sst_cold_skin (:))   !out

      !  Recover fields from index list

      DO ic = 1, i_count
        jc = ext_data%atm%list_seawtr%idx(ic,jb)

        IF (itype_oskin_warm > 0) lnd_diag%sst_warm_layer(jc,jb) = sst_warm_layer(ic)
        IF (itype_oskin_cold > 0) lnd_diag%sst_cold_skin (jc,jb) = sst_cold_skin (ic)

      ENDDO  ! ic

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_oskin

END MODULE mo_nwp_oskin_interface
