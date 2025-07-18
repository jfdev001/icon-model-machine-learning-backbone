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
!  This module contains the saturation adjustment routines.
!
! Routines (module procedure)
!
!     - satad_V_3D
!       Corrects the temperature, the specific humidity and the cloud water
!      content for condensation/evaporation.

MODULE mo_satad


USE mo_kind,               ONLY: wp
USE mo_thdyn_functions,    ONLY: latent_heat_vaporization, qsat_rho, dqsatdT_rho
USE mo_physical_constants, ONLY: cvd
USE mo_nwp_tuning_config,  ONLY: supsat_limfac => tune_supsat_limfac

!** temporary workaround until USE statements in ART have been adjusted
USE mo_thdyn_functions,    ONLY: sat_pres_water, sat_pres_ice
!**

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: satad_v_3D
  PUBLIC  :: satad_v_3D_gpu

!** temporary workaround until USE statements in ART have been adjusted
  PUBLIC  :: sat_pres_water, sat_pres_ice
!**

CONTAINS

SUBROUTINE satad_v_3D (maxiter, tol, te, qve, qce,    & ! IN, INOUT
                       rhotot, w,                     & ! IN
                       idim, kdim, ilo, iup, klo, kup ) ! IN

  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   This routine corrects the temperature (te), the specific humidity (qve),
  !   the cloud water content (qce) and the pressure (pe) for condensation/evaporation.
  !   Pressure adapts itself in ICON but has to be rediagnosed in COSMO
  !
  ! Method:
  !   Saturation adjustment at constant total density (adjustment of T and p accordingly)
  !   assuming chemical equilibrium of water and vapor. For the heat capacity of
  !   of the total system (dry air, vapor, and hydrometeors) the value of dry air
  !   is taken, which is a common approximation and introduces only a small error.
  !
  ! Inout fields: te, qve, qce
  !
  ! Input only fields: rhotot
  !
  ! Outputs:  - count  :  number of iterations needed, maximum of all iterated gridpoints
  !
  ! Description of input/inout - fields:
  !
  ! te:           abs. temperature (adjusted during satad)           [K]
  ! qve:          specific vapor content (adjusted during satad)     [-]
  ! qce:          specific cloud content (adjusted during satad)     [-]
  ! rhotot:       total density, assumed constant during satad       [kg/m^3]
  ! w:            vertical velocity                                  [m/s]
  !
  !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! Subroutine arguments:
  ! --------------------
  INTEGER,   INTENT (IN)    ::  &
       maxiter,                 & !  Max. number of iterations in the numerical scheme
       idim, kdim,              & !  Dimension of I/O-fields
       ilo, iup, klo, kup         !  start- and end-indices for the computations

  REAL    (KIND=wp),    INTENT (IN) ::  &
       tol                 ! Desired abs. accuracy in K of the adjusted temperature

  REAL    (KIND=wp),    INTENT (INOUT), DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       te      , & ! Temperature on input/ouput
       qve     , & ! Specific humidity on input/output
       qce         ! Specific cloud water content on input/output

  REAL    (KIND=wp),    INTENT (IN),  DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       rhotot    ! density containing dry air and water constituents

  REAL    (KIND=wp),    INTENT (IN),  DIMENSION(:,:) ::  &  !  dim (idim,kdim+1)
       w         ! vertical velocity

  INTEGER :: count              ! number of iterations actually needed;
                                ! maximum over all iterated gridpoints

  ! Local variables:
  ! -------------
  INTEGER :: i, k              !  Loop indices
  INTEGER ::                   &
       nsat,                   & !  Number of saturated gridpoints
       iwrk(idim*kdim),        & !  i-index of saturated gridpoints
       kwrk(idim*kdim),        & !  k-index of saturated gridpoints
       indx                        !, & !  loop index

  REAL    (KIND=wp   ) ::  &
       Ttest(idim,kdim), qtest(idim,kdim), qw(idim,kdim), qwd, qwa, dqwd, fT, dfT, &
       zqwmin              ! Minimum cloud water content for adjustment

  REAL    (KIND=wp),  DIMENSION(idim,kdim) ::  &
       lwdocvd,                &  ! (Temperature-dependent) latent heat of vaporization over cv
       supsatfac                  ! Parameterized supersaturation according to updraft

  REAL (KIND=wp), DIMENSION(idim*kdim) ::  &
       twork, tworkold

  !------------ End of header ----------------------------------------------------
  !-------------------------------------------------------------------------------
  ! Begin Subroutine satad
  !-------------------------------------------------------------------------------

  ! Initialization

  zqwmin = 1.0E-20_wp

  ! Counter for the number of gridpoints which need the Newton iteration:
  nsat = 0

!!!=============================================================================================

    DO k = klo, kup
      DO i = ilo , iup

        ! total content of the species which are changed by the adjustment:
        qw(i,k) = qve(i,k) + qce(i,k)

        ! check, which points will still be subsaturated even
        ! if all the cloud water would have been evaporated.
        ! At such points, the Newton iteration is not necessary and the
        ! adjusted values of T, p, qv and qc can be obtained directly.

        lwdocvd(i,k) = latent_heat_vaporization(te(i,k)) / cvd

        Ttest(i,k) = te(i,k) - lwdocvd(i,k)*qce(i,k)

        qtest(i,k) = qsat_rho(Ttest(i,k), rhotot(i,k)) ! Alberto why do not add supsat here?

        IF (supsat_limfac > 0._wp) THEN
          supsatfac(i,k) = 1._wp + MIN(0.005_wp*MAX(0._wp,0.5_wp*(w(i,k)+w(i,k+1))), &
                                       supsat_limfac*qce(i,k)/MAX(zqwmin,qve(i,k))   )
        ELSE
          supsatfac(i,k) = 1._wp
        ENDIF

      END DO
    END DO

    DO k = klo, kup
      DO i = ilo , iup

        IF (qw(i,k) <= qtest(i,k) ) THEN
          ! In this case, all the cloud water evaporates and there is still (sub)saturation.
          ! The resulting state depends only on the available cloud water and is
          ! not saturated, which enables direct computation of the adjusted variables:
          qve(i,k)  = qw(i,k)
          qce(i,k)  = 0.0_wp
          te(i,k)   = Ttest(i,k)
        ELSE
          ! In this case, the Newton interation is needed
          nsat       = nsat+1
          iwrk(nsat) = i
          kwrk(nsat) = k

          ! Field for the iterated temperature, here set the starting value for the below iteration
          ! to the "old" temperature (as an alternative, the arithmetic mean
          ! between "old" temperature and dew point has been tested, but did
          ! not significantly increase the convergence speed of the iteration):
          twork(nsat)  = te(i,k)

          ! And this is the storage variable for the "old" values in the below iteration:
          ! Add some nonesense increment to the starting, which is sufficient to trigger the
          ! iteration below:
          tworkold(nsat) = twork(nsat) + 10.0_wp
        END IF

      END DO
    END DO

    ! Do the Newton iteration at gridpoints which need it:
    ! ---------------------------------------------------

!!!=============================================================================
!!! ..This is the version for the NEC. On a scalar system, it might be better
!!!   to do the iteration within the "indx"-loop, not outside!
!!!=============================================================================


    IF (nsat > 0) THEN

      count = 0
      DO WHILE (ANY(ABS(twork(1:nsat)-tworkold(1:nsat)) > tol) .AND. count < maxiter)
        DO indx = 1, nsat
          ! The following if-clause is necessary to achieve reproducible results.
          ! If it is not applied, all grid points are iterated until the "worst"
          ! gridpoint has reached convergence, so the result on one gridpoint
          ! depends on some other gridpoint on the same processor grid.
          IF (ABS(twork(indx)-tworkold(indx)) > tol) THEN
            ! Here we still have to iterate ...
            i = iwrk(indx)
            k = kwrk(indx)
            tworkold(indx) = twork(indx)
            qwd  = qsat_rho(twork(indx), rhotot(i,k)) * supsatfac(i,k)
            dqwd = dqsatdT_rho(qwd, twork(indx), rhotot(i,k)) * supsatfac(i,k)
            ! Newton:
            fT = twork(indx) - te(i,k) + lwdocvd(i,k)*(qwd - qve(i,k))
            dfT = 1.0_wp + lwdocvd(i,k)*dqwd
            twork(indx) = twork(indx) - fT / dfT;
          END IF
        END DO
        count = count + 1
      END DO

      ! Distribute the results back to gridpoints:
      ! ------------------------------------------

!$NEC ivdep
!DIR$ IVDEP
      DO indx = 1, nsat
        i = iwrk(indx)
        k = kwrk(indx)
        te (i,k) = twork(indx)

        ! We disregard here the extrapolation of qsat from the second-last iteration
        ! step, which is done in the original routine to exactly preserve the internal energy.
        ! This introduces a small error (see the COSMO-Documentation, Part III).
        qwa = qsat_rho(te(i,k), rhotot(i,k))*supsatfac(i,k)
        qce(i,k) = MAX(qw(i,k) - qwa, zqwmin)
        qve(i,k) = qw(i,k) - qce(i,k)
      ENDDO

    END IF

!!!=============================================================================================

END SUBROUTINE satad_v_3D

SUBROUTINE satad_v_3D_gpu (maxiter, tol, te, qve, qce, & ! IN, INOUT
                       rhotot, w,                      & ! IN
                       idim, kdim, ilo, iup, klo, kup  ) ! IN

  !-------------------------------------------------------------------------------
  !
  ! Description:
  !   This routine corrects the temperature (te), the specific humidity (qve),
  !   the cloud water content (qce) and the pressure (pe) for condensation/evaporation.
  !   Pressure adapts itself in ICON but has to be rediagnosed in COSMO
  !
  ! Method:
  !   Saturation adjustment at constant total density (adjustment of T and p accordingly)
  !   assuming chemical equilibrium of water and vapor. For the heat capacity of
  !   of the total system (dry air, vapor, and hydrometeors) the value of dry air
  !   is taken, which is a common approximation and introduces only a small error.
  !
  ! Inout fields: te, qve, qce
  !
  ! Input only fields: rhotot
  !
  ! Outputs:  - count  :  number of iterations needed, maximum of all iterated gridpoints
  !
  ! Description of input/inout - fields:
  !
  ! te:           abs. temperature (adjusted during satad)           [K]
  ! qve:          specific vapor content (adjusted during satad)     [-]
  ! qce:          specific cloud content (adjusted during satad)     [-]
  ! rhotot:       total density, assumed constant during satad       [kg/m^3]
  ! w:            vertical velocity                                  [m/s]
  !
  !-------------------------------------------------------------------------------

  IMPLICIT NONE

  ! Subroutine arguments:
  ! --------------------
  INTEGER,   INTENT (IN)    ::  &
       maxiter,                 & !  Max. number of iterations in the numerical scheme
       idim, kdim,              & !  Dimension of I/O-fields
       ilo, iup, klo, kup         !  start- and end-indices for the computations

  REAL    (KIND=wp),    INTENT (IN) ::  &
       tol                 ! Desired abs. accuracy in K of the adjusted temperature

  REAL    (KIND=wp),    INTENT (INOUT), DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       te      , & ! Temperature on input/ouput
       qve     , & ! Specific humidity on input/output
       qce         ! Specific cloud water content on input/output

  REAL    (KIND=wp),    INTENT (IN),  DIMENSION(:,:) ::  &  !  dim (idim,kdim)
       rhotot    ! density containing dry air and water constituents

  REAL    (KIND=wp),    INTENT (IN),  DIMENSION(:,:) ::  &  !  dim (idim,kdim+1)
       w         ! vertical velocity

  INTEGER :: count              ! number of iterations actually needed;
                                ! maximum over all iterated gridpoints

  ! Local variables:
  ! -------------
  INTEGER  :: i, k                      !  Loop indices

  REAL    (KIND=wp   ) ::  &
       Ttest, qtest, qw, qwd, qwa, dqwd, fT, dfT, & !, cvvmcl, qd,
       zqwmin              ! Minimum cloud water content for adjustment

  REAL    (KIND=wp) ::  &
       lwdocvd  ! (Temperature-dependent) latent heat of vaporization over cv

  REAL (KIND=wp) ::  &
       twork, tworkold, supsatfac

  LOGICAL :: iter_mask


  !------------ End of header ----------------------------------------------------

  !-------------------------------------------------------------------------------
  ! Begin Subroutine satad
  !-------------------------------------------------------------------------------

  ! Initialization

  zqwmin = 1.0E-20_wp

!!!=============================================================================================

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG VECTOR TILE(128, 1)
  DO k = klo, kup
    DO i = ilo , iup
      ! total content of the species which are changed by the adjustment:
      qw = qve(i,k) + qce(i,k)

      ! check, which points will still be subsaturated even
      ! if all the cloud water would have been evaporated.
      ! At such points, the Newton iteration is not necessary and the
      ! adjusted values of T, p, qv and qc can be obtained directly.

      lwdocvd = latent_heat_vaporization(te(i,k)) / cvd
      Ttest = te(i,k) - lwdocvd*qce(i,k)
      qtest = qsat_rho(Ttest, rhotot(i,k))

      IF (supsat_limfac > 0._wp) THEN
        supsatfac = 1._wp + MIN(0.005_wp*MAX(0._wp,0.5_wp*(w(i,k)+w(i,k+1))), &
                                supsat_limfac*qce(i,k)/MAX(zqwmin,qve(i,k))   )
      ELSE
        supsatfac = 1._wp
      ENDIF

      iter_mask = .FALSE.
      twork = 0.0_wp
      IF (qw <= qtest ) THEN
        ! In this case, all the cloud water evaporates and there is still (sub)saturation.
        ! The resulting state depends only on the available cloud water and is
        ! not saturated, which enables direct computation of the adjusted variables:
        qve(i,k)  = qw
        qce(i,k)  = 0.0_wp
        te(i,k)   = Ttest
      ELSE
        ! Field for the iterated temperature, here set the starting value for the below iteration
        ! to the "old" temperature (as an alternative, the arithmetic mean
        ! between "old" temperature and dew point has been tested, but did
        ! not significantly increase the convergence speed of the iteration):
        twork  = te(i,k)
        ! And this is the storage variable for the "old" values in the below iteration:
        ! Add some nonesense increment to the starting, which is sufficient to trigger the
        ! iteration below:
        tworkold = twork + 10.0_wp
        iter_mask = .TRUE.
      END IF

      ! Do the Newton iteration at gridpoints which need it:
      ! ---------------------------------------------------
      !$ACC LOOP SEQ
      DO count = 1, maxiter
        IF (iter_mask) THEN
          IF (ABS(twork-tworkold) > tol) THEN
            ! Here we still have to iterate ...
            tworkold = twork
            qwd  = qsat_rho(twork, rhotot(i,k))*supsatfac
            dqwd = dqsatdT_rho(qwd, twork, rhotot(i,k) )*supsatfac
            ! Newton:
            fT = twork - te(i,k) + lwdocvd*(qwd - qve(i,k))
            dfT = 1.0_wp + lwdocvd*dqwd
            twork = twork - fT / dfT;
          END IF
        END IF
      END DO !while

      ! Distribute the results back to gridpoints:
      ! ------------------------------------------
      IF(iter_mask) THEN
        te (i,k) = twork
        ! We disregard here the extrapolation of qsat from the second-last iteration
        ! step, which is done in the original routine to exactly preserve the internal energy.
        ! This introduces a small error (see the COSMO-Documentation, Part III).
        qwa = qsat_rho(te(i,k), rhotot(i,k))*supsatfac
        qce(i,k) = MAX(qw - qwa, zqwmin)
        qve(i,k) = qw - qce(i,k)
      ENDIF

    ENDDO !i
  ENDDO !k
  !$ACC END PARALLEL

!!!=============================================================================================
!!!=============================================================================================

END SUBROUTINE satad_v_3D_gpu

END MODULE mo_satad
