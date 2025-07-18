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
! Interface between atmosphere physics and the CLEO microphysics model, through a coupler
!

MODULE mo_atmo_cleo_coupling

  USE mo_kind,               ONLY: wp
  USE mo_model_domain,       ONLY: t_patch
  USE mo_run_config,         ONLY: iqv, iqc
  USE mo_fortran_tools,      ONLY: assert_acc_host_only
  USE mo_coupling_utils,     ONLY: cpl_def_field, cpl_put_field, &
                                   cpl_get_field_collection_size

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_atmo_cleo_coupling_post_sync
  PUBLIC :: couple_atmo_to_cleo

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_atmo_cleo_coupling'

  INTEGER :: field_id_temperature
  INTEGER :: field_id_pressure
  INTEGER :: field_id_qvap
  INTEGER :: field_id_qcond
  INTEGER :: field_id_eastward_wind
  INTEGER :: field_id_northward_wind
  INTEGER :: field_id_vertical_wind
  INTEGER :: horizontal_fields_collection_size
  INTEGER :: vertical_wind_collection_size
  INTEGER :: local_horizontal_cells

CONTAINS

  !>
  !! Registers fields required for the coupling between atmosphere and cleo
  !!
  !! This subroutine is called from construct_atmo_coupling.
  !!
  SUBROUTINE construct_atmo_cleo_coupling_post_sync(comp_id, cell_point_id, &
                                                    horizontal_cells, timestepstring)

    INTEGER, INTENT(IN)          :: comp_id
    INTEGER, INTENT(IN)          :: cell_point_id
    INTEGER, INTENT(IN)          :: horizontal_cells
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//':construct_atmo_cleo_coupling_post_sync'

    local_horizontal_cells = horizontal_cells
    horizontal_fields_collection_size = cpl_get_field_collection_size(routine, "cleo", &
                                                                      "cleo_grid", "temperature")
    vertical_wind_collection_size = cpl_get_field_collection_size(routine, "cleo", &
                                                                  "cleo_grid", "vertical_wind")

    CALL cpl_def_field(comp_id, cell_point_id, timestepstring, &
                       "temperature", horizontal_fields_collection_size, field_id_temperature)

    CALL cpl_def_field(comp_id, cell_point_id, timestepstring, &
                       "pressure", horizontal_fields_collection_size, field_id_pressure)

    CALL cpl_def_field(comp_id, cell_point_id, timestepstring, &
                       "qvap", horizontal_fields_collection_size, field_id_qvap)

    CALL cpl_def_field(comp_id, cell_point_id, timestepstring, &
                       "qcond", horizontal_fields_collection_size, field_id_qcond)

    CALL cpl_def_field(comp_id, cell_point_id, timestepstring, &
                       "eastward_wind", horizontal_fields_collection_size, field_id_eastward_wind)

    CALL cpl_def_field(comp_id, cell_point_id, timestepstring, &
                       "northward_wind", horizontal_fields_collection_size, field_id_northward_wind)

    CALL cpl_def_field(comp_id, cell_point_id, timestepstring, &
                       "vertical_wind", vertical_wind_collection_size, field_id_vertical_wind)

  END SUBROUTINE construct_atmo_cleo_coupling_post_sync

  !>
  !! Exchange fields between atmosphere and cleo model
  !!
  !! This subroutine is called from nwp_nh_interface.
  !!
  SUBROUTINE couple_atmo_to_cleo(temperature, pressure, tracers_data, &
                                 eastward_wind, northward_wind, vertical_wind, lacc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':couple_atmo_to_cleo'

    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: temperature(:, :, :)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: pressure(:, :, :)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: tracers_data(:, :, :, :)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: eastward_wind(:, :, :)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: northward_wind(:, :, :)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN) :: vertical_wind(:, :, :)
    LOGICAL,  OPTIONAL,           INTENT(IN) :: lacc           ! If true, use openacc

    INTEGER :: info
    INTEGER :: ierror

    CALL assert_acc_host_only('couple_atmo_to_cleo', lacc)

    CALL cpl_put_field(routine, field_id_temperature, "temperature", &
         &             temperature(1:local_horizontal_cells,         &
         &                         1:horizontal_fields_collection_size, 1))

    CALL cpl_put_field(routine, field_id_pressure, "pressure", &
         &             pressure(1:local_horizontal_cells,      &
         &                      1:horizontal_fields_collection_size, 1))

    CALL cpl_put_field(routine, field_id_qvap, "qvap",                   &
         &             tracers_data(1:local_horizontal_cells,            &
         &                          1:horizontal_fields_collection_size, &
                                    1, iqv))

    CALL cpl_put_field(routine, field_id_qcond, "qcond",                 &
         &             tracers_data(1:local_horizontal_cells,            &
         &                          1:horizontal_fields_collection_size, &
                                    1, iqc))

    CALL cpl_put_field(routine, field_id_vertical_wind, "vertical_wind", &
         &             vertical_wind(1:local_horizontal_cells,           &
         &                           1:vertical_wind_collection_size, 1))

    CALL cpl_put_field(routine, field_id_eastward_wind, "eastward_wind", &
         &             eastward_wind(1:local_horizontal_cells,           &
         &                           1:horizontal_fields_collection_size, 1))

    CALL cpl_put_field(routine, field_id_northward_wind, "northward_wind", &
         &             northward_wind(1:local_horizontal_cells,            &
         &                            1:horizontal_fields_collection_size, 1))
  END SUBROUTINE couple_atmo_to_cleo

END MODULE mo_atmo_cleo_coupling
