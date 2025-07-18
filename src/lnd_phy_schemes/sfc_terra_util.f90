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

!> TERRA utility functions
MODULE sfc_terra_util

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: solve_tridiag
  PUBLIC :: zalfa

  INTERFACE solve_tridiag
    MODULE PROCEDURE solve_tridiag_nomask
    MODULE PROCEDURE solve_tridiag_mask
  END INTERFACE

  REAL(wp), PARAMETER :: zalfa = 1.0_wp !< Implicitness parameter (1: full implicit, 0.5: Cranck-Nicholson).

CONTAINS

SUBROUTINE solve_tridiag_nomask (ivstart, ivend, nvec, nlev, a, b, c, d, out)

  !$ACC ROUTINE GANG

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: nlev

  REAL(wp), INTENT(IN) :: a(nvec, nlev) !< Sub-diagonal (nvec,2:nlev).
  REAL(wp), INTENT(IN) :: b(nvec, nlev) !< Diagonal (nvec,1:nlev).
  REAL(wp), INTENT(INOUT) :: c(nvec, nlev) !< Super-diagonal, modified by elimination (nvec,1:nlev-1).
  REAL(wp), INTENT(INOUT) :: d(nvec, nlev) !< Right-hand side, modified by elimination (nvec,1:nlev).
  REAL(wp), INTENT(INOUT) :: out(nvec, nlev) !< Solution output (nvec,1:nlev).

  REAL(wp) :: denominator

  INTEGER :: i, k

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    c(i,1) = c(i,1)/b(i,1)
    d(i,1) = d(i,1)/b(i,1)
  END DO

  !$ACC LOOP SEQ
  DO k = 2, nlev
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(denominator)
    DO i = ivstart, ivend
      denominator = 1._wp/(b(i,k) - a(i,k)*c(i,k-1))
      c(i,k) = c(i,k) * denominator
      d(i,k) = (d(i,k) - a(i,k)*d(i,k-1)) * denominator
    END DO
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    out(i,nlev) = d(i,nlev)
  END DO

  !$ACC LOOP SEQ
  DO k = nlev-1,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      out(i,k) = d(i,k) - c(i,k)*out(i,k+1)
    END DO
  END DO

END SUBROUTINE solve_tridiag_nomask


SUBROUTINE solve_tridiag_mask (ivstart, ivend, nvec, nlev, a, b, c, d, out, mask)

  !$ACC ROUTINE GANG

  INTEGER, INTENT(IN) :: ivstart
  INTEGER, INTENT(IN) :: ivend
  INTEGER, INTENT(IN) :: nvec
  INTEGER, INTENT(IN) :: nlev

  REAL(wp), INTENT(IN) :: a(nvec, nlev) !< Sub-diagonal (nvec,2:nlev).
  REAL(wp), INTENT(IN) :: b(nvec, nlev) !< Diagonal (nvec,1:nlev).
  REAL(wp), INTENT(INOUT) :: c(nvec, nlev) !< Super-diagonal, modified by elimination (nvec,1:nlev-1).
  REAL(wp), INTENT(INOUT) :: d(nvec, nlev) !< Right-hand side, modified by elimination (nvec,1:nlev).
  REAL(wp), INTENT(INOUT) :: out(nvec, nlev) !< Solution output (nvec,1:nlev).

  !> System is solved where mask is `.TRUE.`.
  LOGICAL, INTENT(IN) :: mask(nvec)

  REAL(wp) :: denominator

  INTEGER :: i, k

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (mask(i)) THEN
      c(i,1) = c(i,1)/b(i,1)
      d(i,1) = d(i,1)/b(i,1)
    END IF
  END DO

  !$ACC LOOP SEQ
  DO k = 2, nlev
    !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(denominator)
    DO i = ivstart, ivend
      IF (mask(i)) THEN
        denominator = 1._wp/(b(i,k) - a(i,k)*c(i,k-1))
        c(i,k) = c(i,k) * denominator
        d(i,k) = (d(i,k) - a(i,k)*d(i,k-1)) * denominator
      END IF
    END DO
  END DO

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  DO i = ivstart, ivend
    IF (mask(i)) THEN
      out(i,nlev) = d(i,nlev)
    END IF
  END DO

  !$ACC LOOP SEQ
  DO k = nlev-1,1,-1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO i = ivstart, ivend
      IF (mask(i)) THEN
        out(i,k) = d(i,k) - c(i,k)*out(i,k+1)
      END IF
    END DO
  END DO

END SUBROUTINE solve_tridiag_mask


END MODULE sfc_terra_util
