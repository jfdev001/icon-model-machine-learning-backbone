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
! @brief Common field definition for the atmosphere-ocean coupling


!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_ocean_coupling_common

  USE mo_coupling_utils,  ONLY: cpl_def_field
  USE mo_grid_config,     ONLY: n_dom

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_ocean_coupling_common' ! Output of module for debug

  PUBLIC :: construct_atmo_ocean_coupling_common
  PUBLIC :: destruct_atmo_ocean_coupling_common
  PUBLIC :: t_out_field_id, t_in_field_id
  PUBLIC :: out_field_ids, in_field_ids

  TYPE t_out_field_id
    INTEGER :: umfl
    INTEGER :: vmfl
    INTEGER :: freshflx
    INTEGER :: heatflx
    INTEGER :: seaice_atm
    INTEGER :: sp10m
    INTEGER :: co2_vmr
    INTEGER :: pres_msl
  END TYPE t_out_field_id

  TYPE t_in_field_id
    INTEGER :: sst
    INTEGER :: oce_u
    INTEGER :: oce_v
    INTEGER :: seaice_oce
    INTEGER :: co2_flx
  END TYPE t_in_field_id

  TYPE(t_out_field_id) :: out_field_ids
  TYPE(t_in_field_id), ALLOCATABLE :: in_field_ids(:)

CONTAINS

  !>
  !! Registers fields required for the coupling between atmo and
  !! ocean
  !!
  !! This subroutine is called from construct_atmo_ocean_coupling.
  !!
  SUBROUTINE construct_atmo_ocean_coupling_common( &
    comp_id, cell_point_id, cell_mask_id, timestepstring, &
    use_ocean_velocity)

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id(0:n_dom)
    INTEGER, INTENT(IN) :: cell_mask_id(0:n_dom)
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    LOGICAL, INTENT(IN) :: use_ocean_velocity

    INTEGER :: jg

    ! define outgoing fields
    ! (in case of nested coupling; cell_point_id(0) contains the combined grid)
    jg = MERGE(0, 1, n_dom > 1)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "surface_downward_eastward_stress", 2, out_field_ids%umfl)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "surface_downward_northward_stress", 2, out_field_ids%vmfl)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "surface_fresh_water_flux", 3, out_field_ids%freshflx)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "total_heat_flux", 4, out_field_ids%heatflx)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "atmosphere_sea_ice_bundle", 2, out_field_ids%seaice_atm)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "10m_wind_speed", 1, out_field_ids%sp10m)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "co2_mixing_ratio", 1, out_field_ids%co2_vmr)

    CALL cpl_def_field( &
      comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
      "sea_level_pressure", 1, out_field_ids%pres_msl)

    ALLOCATE(in_field_ids(n_dom))

    ! define incoming fields
    ! (data is received on all domains)
    DO jg = 1, n_dom
      CALL cpl_def_field( &
        comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
        "sea_surface_temperature", 1, in_field_ids(jg)%sst)

      IF (use_ocean_velocity) THEN
        CALL cpl_def_field( &
          comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
          "eastward_sea_water_velocity", 1, in_field_ids(jg)%oce_u)

        CALL cpl_def_field( &
          comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
          "northward_sea_water_velocity", 1, in_field_ids(jg)%oce_v)
      ELSE
        in_field_ids(jg)%oce_u = -1
        in_field_ids(jg)%oce_v = -1
      END IF

      CALL cpl_def_field( &
        comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
        "ocean_sea_ice_bundle", 3, in_field_ids(jg)%seaice_oce)

      CALL cpl_def_field( &
        comp_id, cell_point_id(jg), cell_mask_id(jg), timestepstring, &
        "co2_flux", 1, in_field_ids(jg)%co2_flx)
    END DO

  END SUBROUTINE construct_atmo_ocean_coupling_common

  SUBROUTINE destruct_atmo_ocean_coupling_common()

    DEALLOCATE(in_field_ids)

  END SUBROUTINE destruct_atmo_ocean_coupling_common

END MODULE mo_atmo_ocean_coupling_common
