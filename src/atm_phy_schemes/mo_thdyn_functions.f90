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

!  Description:
!  This module provides service utilities for meteorological calculations.
!
! Routines (module procedure)
!
!    - qsat_rho
!      Specific humidity at water saturation (with respect to flat surface)
!      depending on the temperature "temp" and the total density "rhotot")
!
!    - dqsatdT_rho
!       Partial derivative of the specific humidity at water saturation with
!       respect to the temperature at constant total density.
!
!    the following functions should  later be replaced by lookup tables
!     from mo_convect_tables!
!     - pres_sat_water
!       Saturation water vapour pressure
!
!     - pres_sat_ice
!       Saturation water vapour pressure
!
!     - spec_humi
!       Specific humidity at saturation pressure

MODULE mo_thdyn_functions


USE mo_kind, ONLY: wp, i4
USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                               rvd_m_o => vtmpc1 , & !! rv/rd-1._wp
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 cvd            , & !!
                                 cl    => clw   , & !! specific heat of water
                                 lwd   => alv   , & !! latent heat of vaporization
                                 led   => als   , & !! latent heat of sublimation
                                 b3    => tmelt , & !!
                                 tmelt

USE mo_lookup_tables_constants, ONLY:  &
                                 b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les , & !!               -- " --
                                 b234i => c5ies     !!               -- " --

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: spec_humi
  PUBLIC  :: sat_pres_water
  PUBLIC  :: sat_pres_ice
  PUBLIC  :: dqsatdT
  PUBLIC  :: dqsatdT_ice
  PUBLIC  :: qsat_rho
  PUBLIC  :: dqsatdT_rho
  PUBLIC  :: latent_heat_vaporization
  PUBLIC  :: latent_heat_sublimation
  PUBLIC  :: latent_heat_melting

  INTEGER, PARAMETER :: ipsat = 1    ! (1) Tetens (1930)
                                     ! (2) Murphy-Koop for liq and ice

  REAL (KIND=wp), PARAMETER :: cp_v = 1850._wp ! specific heat of water vapor J
                                                       ! at constant pressure
                                                       ! (Landolt-Bornstein)
  REAL (KIND=wp), PARAMETER :: ci = 2108.0_wp  ! specific heat of ice

  real(wp), parameter  ::      &
       &  c1i_mk = 9.550426,    &  ! coefficients in Murphy and Koop saturation vapor pressure
       &  c2i_mk = 5723.265,    &  ! over ice and over  liquid water
       &  c3i_mk = 3.53068,     &
       &  c4i_mk = 0.00728332,  &
       &  c1w_mk = 54.842763,   &
       &  c2w_mk = 6763.22,     &
       &  c3w_mk = 4.210,       &
       &  c4w_mk = 0.000367,    &
       &  c5w_mk = 53.878,      &
       &  c6w_mk = 1331.22,     &
       &  c7w_mk = 9.44523,     &
       &  c8w_mk = 0.014025,    &
       &  xi_mk  = 0.0415,      &
       &  t0_mk  = 218.8

CONTAINS

ELEMENTAL FUNCTION sat_pres_water(temp)
  IMPLICIT NONE

  REAL (KIND=wp)              :: sat_pres_water
  REAL (KIND=wp), INTENT(IN)  :: temp

  !$ACC ROUTINE SEQ

  IF (ipsat <= 1) THEN
    sat_pres_water = b1*EXP( b2w*(temp-b3)/(temp-b4w) )
  ELSEIF (ipsat == 2) THEN
    sat_pres_water = psatw_murphykoop(temp)
  END IF

END FUNCTION sat_pres_water

!==============================================================================

ELEMENTAL FUNCTION sat_pres_ice(temp)
  IMPLICIT NONE

  REAL (KIND=wp)              :: sat_pres_ice
  REAL (KIND=wp), INTENT(IN)  :: temp

  !$ACC ROUTINE SEQ

  IF (ipsat <= 1) THEN
    sat_pres_ice = b1*EXP( b2i*(temp-b3)/(temp-b4i) )
  ELSEIF (ipsat==2 .OR. ipsat==3) THEN
    sat_pres_ice = psati_murphykoop(temp)
  ENDIF

END FUNCTION sat_pres_ice

!==============================================================================

ELEMENTAL FUNCTION spec_humi(pvap,pres)
  IMPLICIT NONE

  REAL (KIND=wp)              :: spec_humi
  REAL (KIND=wp), INTENT(IN)  :: pres,pvap

  !$ACC ROUTINE SEQ

  spec_humi = rdv*pvap/( pres - o_m_rdv*pvap )

END FUNCTION spec_humi

! Routine with ACC ROUTINE SEQ cannot be elemental
#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION qsat_rho(temp, rhotot)
  !$ACC ROUTINE SEQ

    !-------------------------------------------------------------------------------
    !
    ! Description:
    !   Specific humidity at water saturation (with respect to flat surface)
    !   depending on the temperature "temp" and the total density "rhotot")
    !
    !-------------------------------------------------------------------------------

  ! input variables: temperature [K], total density [kg/m^2]:
  REAL (KIND=wp)             :: qsat_rho
  REAL (KIND=wp), INTENT(IN) :: temp, rhotot

  qsat_rho   = sat_pres_water(temp) / (rhotot * r_v * temp)

END FUNCTION qsat_rho

! Routine with ACC ROUTINE SEQ cannot be elemental
#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION dqsatdT_rho(zqsat, temp, rho)
  !$ACC ROUTINE SEQ

    !-------------------------------------------------------------------------------
    !
    ! Description:
    !   Partial derivative of the specific humidity at water saturation with
    !   respect to the temperature at constant total density.
    !   Depends on the temperature "temp" and the
    !   saturation specific humidity "zqsat".
    !
    !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input variables: specific saturation humidity [-], temperature [K]:
  REAL (KIND=wp)            :: dqsatdT_rho
  REAL (KIND=wp), INTENT(IN):: zqsat, temp, rho

  ! local variables:
  REAL(kind=wp) :: beta

  IF (ipsat == 1) THEN
    beta        = b234w/(temp-b4w)**2_i4 - 1.0_wp / temp
    dqsatdT_rho = beta * zqsat
  ELSEIF (ipsat == 2) THEN
    beta        = psatw_murphykoop_derivative(temp)
    dqsatdT_rho = beta/(r_v*rho*temp)- zqsat/temp
  END IF

END FUNCTION dqsatdT_rho

! UB_20100525<<

ELEMENTAL  FUNCTION dqsatdT (zqsat, temp)

    !-------------------------------------------------------------------------------
    !>
    !! Description:
    !!   Partial derivative of the specific humidity at water saturation with
    !!   respect to the temperature
    !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input variables: specific saturation humidity [-], temperature [K]:
  REAL (KIND=wp)            :: dqsatdT
  REAL (KIND=wp), INTENT(IN):: zqsat, temp

  dqsatdT = b234w * ( 1.0_wp + rvd_m_o*zqsat ) * zqsat &
                             / (temp-b4w)**2

END FUNCTION dqsatdT

ELEMENTAL  FUNCTION dqsatdT_ice (zqsat, temp)

    !-------------------------------------------------------------------------------
    !>
    !! Description:
    !!   Partial derivative of the specific humidity at ice saturation with
    !!   respect to the temperature
    !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input variables: specific saturation humidity [-], temperature [K]:
  REAL (KIND=wp)            :: dqsatdT_ice
  REAL (KIND=wp), INTENT(IN):: zqsat, temp
  !$ACC ROUTINE SEQ

  dqsatdT_ice = b234i * ( 1.0_wp + rvd_m_o*zqsat ) * zqsat &
                             / (temp-b4i)**2

END FUNCTION dqsatdT_ice

! Routine with ACC ROUTINE SEQ cannot be elemental
#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION latent_heat_vaporization(temp)

    !-------------------------------------------------------------------------------
    !>
    !! Description:
    !!   Latent heat of vaporization as internal energy and taking into account
    !!   Kirchoff's relations
    !-------------------------------------------------------------------------------

  IMPLICIT NONE
  REAL(KIND=wp)             :: latent_heat_vaporization
  REAL(KIND=wp), INTENT(IN) :: temp
  !$ACC ROUTINE SEQ

  latent_heat_vaporization = lwd + (cp_v - cl)*(temp-tmelt) - r_v*temp

END FUNCTION latent_heat_vaporization

! Routine with ACC ROUTINE SEQ cannot be elemental
#ifndef _OPENACC
ELEMENTAL &
#endif
FUNCTION latent_heat_sublimation(temp)

    !-------------------------------------------------------------------------------
    !>
    !! Description:
    !!   Latent heat of sublimation as internal energy and taking into account
    !!   Kirchoff's relations
    !-------------------------------------------------------------------------------

  IMPLICIT NONE
  REAL(KIND=wp)             :: latent_heat_sublimation
  REAL(KIND=wp), INTENT(IN) :: temp
  !$ACC ROUTINE SEQ

  latent_heat_sublimation = led + (cp_v - ci)*(temp-tmelt) - r_v*temp

END FUNCTION latent_heat_sublimation

ELEMENTAL FUNCTION latent_heat_melting(temp)

    !-------------------------------------------------------------------------------
    !>
    !! Description:
    !!   Latent heat of sublimation as internal energy and taking into account
    !!   Kirchoff's relations
    !-------------------------------------------------------------------------------

  !$ACC ROUTINE SEQ

  IMPLICIT NONE
  REAL(KIND=wp)             :: latent_heat_melting
  REAL(KIND=wp), INTENT(IN) :: temp

  latent_heat_melting = lwd - led + (ci - cl)*(temp-tmelt)

END FUNCTION latent_heat_melting

ELEMENTAL FUNCTION psatw_murphykoop(tk)
  IMPLICIT NONE
  REAL(KIND=wp)             :: psatw_murphykoop
  REAL(KIND=wp), intent(IN) :: tk

  ! Eq. (10) of Murphy and Koop (2005)
  psatw_murphykoop = exp( c1w_mk - c2w_mk/tk - c3w_mk*log(tk) + c4w_mk*tk &
       & + tanh(xi_mk*(tk-t0_mk))*(c5w_mk-c6w_mk/tk-c7w_mk*log(tk)+c8w_mk*tk) )

END FUNCTION psatw_murphykoop

ELEMENTAL FUNCTION psati_murphykoop(tk)
  IMPLICIT NONE
  REAL(KIND=wp)             :: psati_murphykoop
  REAL(KIND=wp), intent(IN) :: tk

  ! Eq. (7) of Murphy and Koop (2005)
  psati_murphykoop = exp(c1i_mk - c2i_mk/tk + c3i_mk*log(tk) - c4i_mk*tk )

END FUNCTION psati_murphykoop

ELEMENTAL FUNCTION psatw_murphykoop_derivative(tk)
  IMPLICIT NONE
  REAL(KIND=wp)             :: psatw_murphykoop_derivative
  REAL(KIND=wp), intent(IN) :: tk

  ! Derivative of Eq. (10) of Murphy and Koop (2005)
  psatw_murphykoop_derivative = &
       & (c2w_mk/tk**2 - c3w_mk/tk + c4w_mk) * exp( c1w_mk - c2w_mk/tk - c3w_mk*log(tk) + c4w_mk*tk)  &
       & + cosh( xi_mk*(tk-t0_mk))**(-2) * (c5w_mk - c6w_mk/tk - c7w_mk*log(tk) + c8w_mk*tk) &
       & + tanh( xi_mk*(tk-t0_mk)) * (c6w_mk/tk**2 - c7w_mk/tk + c8w_mk)

END FUNCTION psatw_murphykoop_derivative

ELEMENTAL FUNCTION psati_murphykoop_derivative(tk)
  IMPLICIT NONE
  REAL(KIND=wp)             :: psati_murphykoop_derivative
  REAL(KIND=wp), intent(IN) :: tk

  ! Derivative of Eq. (7) of Murphy and Koop (2005)
  psati_murphykoop_derivative = (c2i_mk/tk**2 + c3i_mk/tk - c4i_mk) * psati_murphykoop(tk)

END FUNCTION psati_murphykoop_derivative


END MODULE mo_thdyn_functions
