!
! mo_art_chemical_constants
! This module defines chemical constants for ICON-ART
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

MODULE mo_art_chemical_constants

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE

  PUBLIC

  !> Atmospheric gas fractions
  REAL(wp), PARAMETER :: n2_frac = 0.78084_wp  ! Fraction of nitrogen in dry air
  REAL(wp), PARAMETER :: o2_frac = 0.20946_wp  ! Fraction of oxygen in dry air

END MODULE mo_art_chemical_constants
