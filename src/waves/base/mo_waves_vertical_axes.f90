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

! Specification of vertical axes for the wave model

MODULE mo_waves_vertical_axes

  USE mo_kind,                              ONLY: dp
  USE mo_zaxis_type,                        ONLY: ZA_SURFACE, ZA_HEIGHT_10M,  &
    &                                             ZA_FREQ_GENERIC, ZA_DIR_GENERIC
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxisList
  USE mo_name_list_output_zaxes,            ONLY: single_level_axis, vertical_axis
  USE mo_level_selection_types,             ONLY: t_level_selection
  USE mo_wave_config,                       ONLY: t_wave_config, wave_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setup_zaxes_waves

CONTAINS


  SUBROUTINE setup_zaxes_waves(verticalAxisList, level_selection, log_patch_id)
    TYPE(t_verticalAxisList), INTENT(INOUT)       :: verticalAxisList
    TYPE(t_level_selection),  INTENT(IN), POINTER :: level_selection  ! in general non-associated for waves
    INTEGER,                  INTENT(IN)          :: log_patch_id

    ! local
    INTEGER :: k
    TYPE(t_wave_config), POINTER :: wc

    ! convenience pointer
    wc => wave_config(log_patch_id)

    ! --------------------------------------------------------------------------------------
    ! Definitions for single levels --------------------------------------------------------
    ! --------------------------------------------------------------------------------------

    ! surface level
    CALL verticalAxisList%append(single_level_axis(ZA_surface))

    ! Specified height level above ground: 10m
    CALL verticalAxisList%append(single_level_axis(ZA_height_10m, opt_level_value=10._dp))

    ! ZA_FREQ_GENERIC
    ! vertical axis for 3D fields which are a function of frequency
    CALL verticalAxisList%append(vertical_axis(                             &
      &                           za_type          = ZA_FREQ_GENERIC,       &
      &                           in_nlevs         = wc%nfreqs,             &
      &                           levels           = REAL(wc%freqs(:),dp),  &
      &                           level_selection  = level_selection,       &
      &                           opt_name         = "freq",                &
      &                           opt_unit         = "s-1"))

    ! ZA_DIR_GENERIC
    ! vertical axis for 3D fields which are a function of direction
    CALL verticalAxisList%append(vertical_axis(                             &
      &                           za_type          = ZA_DIR_GENERIC,        &
      &                           in_nlevs         = wc%ndirs,              &
      &                           levels           = REAL(wc%dirs(:),dp),   &
      &                           level_selection  = level_selection,       &
      &                           opt_name         = "dir",                 &
      &                           opt_unit         = "rad"))

  END SUBROUTINE setup_zaxes_waves

END MODULE mo_waves_vertical_axes
