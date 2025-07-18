! This module provides the functionality for a radiation multicall within the nwp physics.
!
! The radiation multicall calls the radiation routine multiple times instead of once
! with different components of aerosols included in each call. This allows calculating
! the direct radiation effect (DRE) for each aerosol species/mode explicitly.
! Further output variables for the DRE are included.
!
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
MODULE mo_art_radiation_multicall

  ! USE <module>, ONLY:
  USE mo_exception,           ONLY: finish, message_text
  USE mo_kind,                ONLY: wp, vp
  USE mo_model_domain,         ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_nonhydro_types,          ONLY: t_nh_diag
  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_impl_constants,          ONLY: success, min_rlcell_int
  USE mo_run_config,              ONLY: msg_level
  ! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_modes_linked_list,         ONLY: p_mode_state
  USE mo_art_diag_types,                ONLY: t_art_diag

  IMPLICIT NONE

  PRIVATE

  INTEGER :: ncallsrad, nccrad  !< number of total calls and current call

  PUBLIC  :: rad_multicall_dre, &
    &        rad_multicall_init, &
    &        rad_multicall_alloc, &
    &        rad_multicall_dealloc, &
    &        rad_multicall_finalize, &
    &        ncallsrad, nccrad, t_dre_flx_ref


  ! !>
  ! !! Array for flux references
  ! !!
  ! !! The type just provides a handle for the four flux
  ! !! references required by using the same type definition.
  ! !!
  TYPE :: t_dre_flx_ref
    REAL(wp), ALLOCATABLE ::   &
      swflxsfc(:,:),   & !< Flux reference: short wave, surface
      swflxtoa(:,:),   & !< Flux reference: short wave, top of atmosphere
      lwflxsfc(:,:),   & !< Flux reference: long  wave, surface
      lwflxtoa(:,:)      !< Flux reference: long  wave, top of atmosphere
  END TYPE t_dre_flx_ref

CONTAINS

! This subroutine initializes the multicall.
!
! Set the number of calls used in the multicall scheme.
!
SUBROUTINE rad_multicall_init (jg)

  INTEGER, INTENT(IN)    :: jg          !<  Domain index
  INTEGER           :: icall_type  !<  Type of the multicall scheme
  
  icall_type = art_config(jg)%irad_multicall
  IF ( icall_type > 0 .AND. art_config(jg)%iart_ari == 0 ) THEN
    CALL finish('art_radiation_multicall:', "multicall is only supported with online aerosol radiation interaction (iart_ari=1)")
  ENDIF
  
  IF ( icall_type == 0 ) THEN ! no multiple call
    ncallsrad = 1
  ELSEIF ( icall_type == 1 ) THEN ! double call
    ncallsrad = 2
  ELSEIF ( icall_type == 2 .OR. icall_type == 3 ) THEN ! multicall
    ncallsrad = 2 + p_mode_state(jg)%p_mode_list%p%nmodes
  ELSE
    WRITE (message_text, '(a,i2,a)') 'irad_multicall = ', icall_type, ' is not supported!'
    CALL finish('art_radiation_multicall:', message_text)
  ENDIF
  nccrad = 1

END SUBROUTINE rad_multicall_init



! This subroutine allocates memory for the multicall.
!
! Arrays for the flux references have to be allocated.
! These references save the fluxes of the first call
! (without aerosols) of the current multicall loop.
!
SUBROUTINE rad_multicall_alloc (dre_ref, kblks)

  INTEGER, INTENT(IN)    :: kblks !< dimension size
  TYPE(t_dre_flx_ref), INTENT(INOUT) :: dre_ref
  INTEGER :: ist

  ! Allocate memory
  ALLOCATE(dre_ref%swflxsfc(nproma, kblks), &
    &      dre_ref%swflxtoa(nproma, kblks), &
    &      dre_ref%lwflxsfc(nproma, kblks), &
    &      dre_ref%lwflxtoa(nproma, kblks), &
    &      STAT=ist)
  IF(ist/=success) CALL finish('mo_art_radiation_multicall:rad_multicall_alloc', 'allocation of DRE flux reference arrays failed')
  
  ! Set references to zero.
  dre_ref%swflxsfc (:,:) = 0._wp
  dre_ref%swflxtoa (:,:) = 0._wp
  dre_ref%lwflxsfc (:,:) = 0._wp
  dre_ref%lwflxtoa (:,:) = 0._wp

END SUBROUTINE rad_multicall_alloc



! This subroutine deallocates memory for the multicall.
!
! Deallocate the memory for the flux references.
!
SUBROUTINE rad_multicall_dealloc (dre_ref)

  TYPE(t_dre_flx_ref), INTENT(INOUT) :: dre_ref
  INTEGER :: ist

  DEALLOCATE(dre_ref%swflxsfc, &
    &        dre_ref%swflxtoa, &
    &        dre_ref%lwflxsfc, &
    &        dre_ref%lwflxtoa, &
    &        STAT=ist)
  IF(ist/=success) CALL finish('mo_art_radiation_multicall:rad_multicall_alloc', 'allocation of DRE flux reference arrays failed')

END SUBROUTINE rad_multicall_dealloc



! This subroutine calculates the direct radiative effect (DRE).
!
! The direct radiative effect is calculated by taking the difference
! of the current fluxes to the reference saved in the first call.
! If it is the first call these references are saved instead.
!
SUBROUTINE rad_multicall_dre (prm_diag, dre_flx_ref, pt_patch)

  TYPE(t_patch), TARGET, INTENT(IN)    :: pt_patch     !<grid/patch info.

  TYPE(t_nwp_phy_diag),  INTENT(IN)    :: prm_diag
  TYPE(t_dre_flx_ref),   INTENT(INOUT) :: dre_flx_ref

  TYPE(t_art_diag),      POINTER       :: art_diag !< Pointer to ART diagnostic fields

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: jc,jb                   !< loop indices
  
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = pt_patch%cells%start_block(rl_start)
  i_endblk   = pt_patch%cells%end_block(rl_end)
  art_diag => p_art_data(pt_patch%id)%diag

  IF ( nccrad == 1 ) THEN  ! First call - no aerosol
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            &            i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        dre_flx_ref%swflxsfc (jc,jb) = prm_diag%swflxsfc (jc,jb)
        dre_flx_ref%swflxtoa (jc,jb) = prm_diag%swflxtoa (jc,jb)
        dre_flx_ref%lwflxsfc (jc,jb) = prm_diag%lwflxsfc (jc,jb)
        dre_flx_ref%lwflxtoa (jc,jb) = prm_diag%lwflxtoa (jc,jb)
      ENDDO
    ENDDO ! blocks
  ELSEIF ( nccrad > 1 ) THEN   ! Later calls - with aerosols
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            &            i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        art_diag%dre_sw_sfc(jc, jb, nccrad - 1) = prm_diag%swflxsfc (jc,jb) - dre_flx_ref%swflxsfc (jc,jb)
        art_diag%dre_sw_toa(jc, jb, nccrad - 1) = prm_diag%swflxtoa (jc,jb) - dre_flx_ref%swflxtoa (jc,jb)
        art_diag%dre_lw_sfc(jc, jb, nccrad - 1) = prm_diag%lwflxsfc (jc,jb) - dre_flx_ref%lwflxsfc (jc,jb)
        art_diag%dre_lw_toa(jc, jb, nccrad - 1) = prm_diag%lwflxtoa (jc,jb) - dre_flx_ref%lwflxtoa (jc,jb)
      ENDDO
    ENDDO ! blocks
  ENDIF

END SUBROUTINE rad_multicall_dre



! This subroutine performs necessary calculations after all calls.
!
! It calculates (allaero DRE - calculated DRE) for the call type 'irad_multicall = 3'
! Thus, one obtains the DRE for the specific mode (not all except the mode).
!
SUBROUTINE rad_multicall_finalize (pt_patch)

  TYPE(t_patch), TARGET,   INTENT(IN)    :: pt_patch     !<grid/patch info.
  TYPE(t_art_diag), POINTER :: art_diag !< Pointer to ART diagnostic fields

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: jc,jb,nc                !< loop indices
  
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = pt_patch%cells%start_block(rl_start)
  i_endblk   = pt_patch%cells%end_block(rl_end)
  art_diag => p_art_data(pt_patch%id)%diag

  IF ( art_config(pt_patch%id)%irad_multicall == 3 ) THEN   ! Last call requires
    ! Loop over all calls except the last one (for allaero nothing) has to be done.
    ! For all other calls, subtract the calculated DRE from the allaero DRE to
    ! obtain the DRE for the specific mode (not all except the mode).
    DO nc = 1, ncallsrad - 2
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
              &            i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx
          art_diag%dre_sw_sfc(jc, jb, nc) = art_diag%dre_sw_sfc(jc, jb, ncallsrad - 1) - art_diag%dre_sw_sfc(jc, jb, nc)
          art_diag%dre_sw_toa(jc, jb, nc) = art_diag%dre_sw_toa(jc, jb, ncallsrad - 1) - art_diag%dre_sw_toa(jc, jb, nc)
          art_diag%dre_lw_sfc(jc, jb, nc) = art_diag%dre_lw_sfc(jc, jb, ncallsrad - 1) - art_diag%dre_lw_sfc(jc, jb, nc)
          art_diag%dre_lw_toa(jc, jb, nc) = art_diag%dre_lw_toa(jc, jb, ncallsrad - 1) - art_diag%dre_lw_toa(jc, jb, nc)
        ENDDO
      ENDDO ! blocks
    ENDDO
  ENDIF

  ! Loop over all calls and add the current DRE to the accumulated quantity.
  DO nc = 1, ncallsrad - 1
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
            &            i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        art_diag%acc_dre_sw_sfc(jc, jb, nc) = art_diag%acc_dre_sw_sfc(jc, jb, nc) + art_diag%dre_sw_sfc(jc, jb, nc)
        art_diag%acc_dre_sw_toa(jc, jb, nc) = art_diag%acc_dre_sw_toa(jc, jb, nc) + art_diag%dre_sw_toa(jc, jb, nc)
        art_diag%acc_dre_lw_sfc(jc, jb, nc) = art_diag%acc_dre_lw_sfc(jc, jb, nc) + art_diag%dre_lw_sfc(jc, jb, nc)
        art_diag%acc_dre_lw_toa(jc, jb, nc) = art_diag%acc_dre_lw_toa(jc, jb, nc) + art_diag%dre_lw_toa(jc, jb, nc)
      ENDDO
    ENDDO ! blocks
  ENDDO

END SUBROUTINE rad_multicall_finalize


END MODULE mo_art_radiation_multicall

