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
! Subroutine interface_cleo calls the coupling to the CLEO model
!

MODULE mo_interface_cleo

  USE mo_aes_phy_memory,      ONLY: prm_field
  USE mo_atmo_cleo_coupling,  ONLY: couple_atmo_to_cleo

  PUBLIC :: interface_cleo

  CONTAINS

  SUBROUTINE interface_cleo(jg)
      INTEGER, INTENT(IN) :: jg

      CALL couple_atmo_to_cleo(temperature      = prm_field(jg)%ta, &
      &                        pressure         = prm_field(jg)%pfull, &
      &                        tracers_data     = prm_field(jg)%qtrc_phy, &
      &                        eastward_wind    = prm_field(jg)%ua, &
      &                        northward_wind   = prm_field(jg)%va, &
      &                        vertical_wind    = prm_field(jg)%wa)
  END SUBROUTINE interface_cleo

END MODULE mo_interface_cleo
