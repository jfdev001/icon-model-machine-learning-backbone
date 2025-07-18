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

! @brief calculate indices and weights for a linear interpolation of
!   a zonal climatology to the icon latitudes.
!   Assumption: The climatology is ordered from North to South or
!   South to North, it has equally spaced latitudes but a shift
!   with respect to the poles is allowed that is different from
!   the other spacing, (e.g. Pi/2, Pi/6, 0., -Pi/6, -Pi/2), the shift would be Pi/3.
!   or (-Pi/2, -Pi/6, 0., Pi/6, Pi/2) with also a shift of Pi/3. Latitudes have to
!   be given in radiant. The extrapolation to the poles is done by repeating the value
!   at the next lower latitude.

MODULE mo_latitude_interpolation

  USE mo_kind,                     ONLY: wp
  USE mo_model_domain,             ONLY: p_patch
  USE mo_impl_constants,           ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_loopindices,              ONLY: get_indices_c
  USE mo_fortran_tools,            ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: latitude_weights_li

  CONTAINS

!> SUBROUTINE latitude_weights_li  -- calculate weights and indices for
!             linear latitude interpolation.

  SUBROUTINE latitude_weights_li(jg                   ,jcs                            &
                               & ,kproma              ,kbdim            ,krow         &
                               & ,wgt1_lat            ,wgt2_lat         ,inmw1_lat    &
                               & ,inmw2_lat           ,p_lat_shift      ,p_rdeltalat  &
                               & ,r_lat_clim          ,nlat_clim        ,n_order      &
                               & ,lacc                                                )

    ! n_order=1 if latitudes of climatology are in ascending (S->N), -1 if
    ! latitudes are in descending (N->S) order.
    INTEGER, INTENT(in)               :: jg,        & ! domain index
                                       & jcs,       & ! actual block length (start)
                                       & kproma,    & ! actual block length (end)
                                       & kbdim,     & ! maximal block length
                                       & krow,      & ! block index
                                       & n_order      ! =1 if latitudes in climatology are ordered S->N
                                                      ! =-1 if latitudes in climatology are ordered N->S
    REAL(wp), INTENT(inout)           :: wgt1_lat(kbdim), wgt2_lat(kbdim) ! linear interpolation weights
    INTEGER, INTENT(inout)            :: inmw1_lat(kbdim), inmw2_lat(kbdim) ! linear interpolation indices
    REAL(wp), INTENT(in)              :: p_lat_shift,&! shift of latitudes with respect to pole (see above)
                                       & p_rdeltalat  ! spacing of latitudes in climatology
    INTEGER, INTENT(in)               :: nlat_clim    !number of latitudes minus the values at the poles
    REAL(wp), INTENT(in)              :: r_lat_clim(0:nlat_clim+1)! latitudes of climatology.
                                                      ! ATTENTION: they must contain the poles
                                                      ! r_lat_clim(0)=+-Pi/2, r_lat_clim(nlat_clim+1)=+-Pi/2

    REAL(wp)                          :: zlat
    INTEGER                           :: jc

    LOGICAL, OPTIONAL, INTENT(IN)     :: lacc
    LOGICAL                           :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(zlat)
    DO jc = jcs, kproma
        zlat=p_patch(jg)%cells%center(jc,krow)%lat

        inmw1_lat(jc)=MAX(INT(n_order*(zlat-p_lat_shift)*p_rdeltalat+1),0)
        inmw2_lat(jc)=inmw1_lat(jc)+1
        wgt2_lat(jc)=n_order*(zlat-r_lat_clim(inmw1_lat(jc)))*p_rdeltalat
        wgt1_lat(jc)=1.0_wp-wgt2_lat(jc)
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE latitude_weights_li

END MODULE mo_latitude_interpolation
