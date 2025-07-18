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

! wave specific runtime timer

MODULE mo_wave_timer

  USE mo_real_timer,  ONLY: new_timer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_wave_timer

  PUBLIC :: timer_wave_total
  PUBLIC :: timer_wave_reader
  PUBLIC :: timer_wave_src
  PUBLIC :: timer_wave_src_wind_input
  PUBLIC :: timer_wave_src_dissipation
  PUBLIC :: timer_wave_src_nonlinear
  PUBLIC :: timer_wave_propagation
  PUBLIC :: timer_wave_energy_propagation
  PUBLIC :: timer_wave_grid_refraction
  PUBLIC :: timer_wave_time_integration
  PUBLIC :: timer_wave_diagnostics


  INTEGER :: timer_wave_total, &
    &        timer_wave_reader, &
    &        timer_wave_src, &
    &        timer_wave_src_wind_input, timer_wave_src_dissipation,  timer_wave_src_nonlinear, &
    &        timer_wave_propagation, &
    &        timer_wave_energy_propagation, timer_wave_grid_refraction, &
    &        timer_wave_time_integration, &
    &        timer_wave_diagnostics

CONTAINS

  !
  ! Initialization routine for wave-specific timer
  !
  SUBROUTINE init_wave_timer (ltimer)

    LOGICAL, INTENT(IN) :: ltimer  ! main switch

    IF (.NOT. ltimer)  RETURN

    timer_wave_total              = new_timer("wave_total")
    timer_wave_reader             = new_timer("wave_forcing_reader")
    timer_wave_propagation        = new_timer("wave_propagation")
    timer_wave_energy_propagation = new_timer("wave_energy_propagation")
    timer_wave_grid_refraction    = new_timer("wave_grid_refraction")
    timer_wave_src                = new_timer("wave_source")
    timer_wave_src_wind_input     = new_timer("wave_source_wind_input")
    timer_wave_src_dissipation    = new_timer("wave_source_dissipation")
    timer_wave_src_nonlinear      = new_timer("wave_source_non_linear")
    timer_wave_time_integration   = new_timer("wave_time_integration")
    timer_wave_diagnostics        = new_timer("wave_diagnostics")

  END SUBROUTINE init_wave_timer

END MODULE mo_wave_timer
