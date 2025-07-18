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

! @brief Interface between atmosphere physics and the ocean, through a coupler

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_ocean_coupling

  USE mo_kind,            ONLY: wp
  USE mo_model_domain,    ONLY: t_patch
  USE mo_ext_data_types,  ONLY: t_external_data
#ifndef __NO_AES__
  USE mo_aes_phy_memory,  ONLY: prm_field
#endif
  USE mo_parallel_config, ONLY: nproma
  USE mo_grid_config,     ONLY: n_dom
  USE mo_impl_constants,  ONLY: inwp, iaes, SUCCESS
  USE mo_mpi,             ONLY: p_pe_work, p_comm_work, p_lor
  USE mo_run_config,      ONLY: iforcing
  USE mo_util_dbg_prnt,   ONLY: dbg_print
  USE mo_exception,       ONLY: finish
  USE mo_coupling_utils,  ONLY: cpl_def_cell_field_mask
  USE mo_atmo_ocean_coupling_common, ONLY: construct_atmo_ocean_coupling_common, &
                                           destruct_atmo_ocean_coupling_common
  USE mo_aes_ocean_coupling, ONLY: construct_aes_ocean_coupling, &
                                   destruct_aes_ocean_coupling

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_ocean_coupling' ! Output of module for debug

  PUBLIC :: construct_atmo_ocean_coupling
  PUBLIC :: destruct_atmo_ocean_coupling

CONTAINS

  !>
  !! Registers fields required for the coupling between atmo and
  !! ocean
  !!
  !! This subroutine is called from construct_atmo_coupling.
  !!
  SUBROUTINE construct_atmo_ocean_coupling( &
    p_patch, ext_data, comp_id, grid_id, cell_point_id, timestepstring, &
    use_ocean_velocity)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    TYPE(t_external_data), INTENT(IN) :: ext_data(:)
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: grid_id(0:)
    INTEGER, INTENT(IN) :: cell_point_id(0:)
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring
    LOGICAL, INTENT(IN) :: use_ocean_velocity

    LOGICAL, ALLOCATABLE :: is_valid(:)
    INTEGER :: cell_mask_id(0:n_dom)
    LOGICAL :: use_mask(1:n_dom)

    INTEGER :: error, jg, total_num_cell_blocks, total_num_cells

    CHARACTER(LEN=*), PARAMETER   :: &
      routine = str_module // ':construct_atmo_ocean_coupling'

    IF (n_dom > 1) THEN

      total_num_cell_blocks = 0
      DO jg = 1, n_dom
        total_num_cell_blocks = &
          total_num_cell_blocks + p_patch(jg)%nblks_c
      END DO

      ALLOCATE(is_valid(nproma * total_num_cell_blocks), STAT = error)
      IF(error /= SUCCESS) &
        CALL finish(routine, "memory allocation failure for is_valid")

      total_num_cells = 0
      DO jg = 1, n_dom
        CALL def_cell_field_mask( &
          p_patch(jg), ext_data(jg), grid_id(jg), &
          cell_mask_id(jg), use_mask(jg), &
          is_valid( &
            total_num_cells + 1:total_num_cells + nproma * p_patch(jg)%nblks_c))
        total_num_cells = total_num_cells + p_patch(jg)%n_patch_cells
      END DO


      CALL cpl_def_cell_field_mask( &
        routine, grid_id(0), is_valid(1:total_num_cells), cell_mask_id(0))

      DEALLOCATE (is_valid, STAT = error)
      IF(error /= SUCCESS) &
        CALL finish(routine, "Deallocation failed for is_valid")

    ELSE

      jg = 1
      cell_mask_id = -1
      use_mask = .FALSE.

      ALLOCATE(is_valid(nproma * p_patch(jg)%nblks_c), STAT = error)
      IF(error /= SUCCESS) &
        CALL finish(routine, "memory allocation failure for is_valid")

      CALL def_cell_field_mask( &
        p_patch(jg), ext_data(jg), grid_id(jg), &
        cell_mask_id(jg), use_mask(jg), is_valid)

      DEALLOCATE (is_valid, STAT = error)
      IF(error /= SUCCESS) &
        CALL finish(routine, "Deallocation failed for is_valid")
    END IF

    CALL construct_atmo_ocean_coupling_common( &
      comp_id, cell_point_id, cell_mask_id, timestepstring, use_ocean_velocity)

    IF (iforcing == iaes) THEN
      CALL construct_aes_ocean_coupling(p_patch, use_mask, use_ocean_velocity)
    END IF

  CONTAINS

    SUBROUTINE def_cell_field_mask( &
      p_patch, ext_data, grid_id, cell_mask_id, use_mask, is_valid)

      TYPE(t_patch), INTENT(IN) :: p_patch
      TYPE(t_external_data), INTENT(IN) :: ext_data
      INTEGER, INTENT(IN) :: grid_id
      INTEGER, INTENT(OUT) :: cell_mask_id
      LOGICAL, INTENT(OUT) :: use_mask
      LOGICAL, INTENT(INOUT) :: is_valid(:)

      INTEGER :: jg, jb, jc, error
      REAL(wp), ALLOCATABLE :: lsmnolake(:,:)

      REAL(wp), PARAMETER :: eps = 1.E-10_wp

      CHARACTER(LEN=*), PARAMETER   :: &
        routine = str_module // ':def_cell_field_mask'

      jg = p_patch%id

      ! The integer land-sea mask:
      !          -2: inner ocean
      !          -1: boundary ocean
      !           1: boundary land
      !           2: inner land
      !
      ! The (fractional) mask which is used in the AES physics is prm_field(1)%lsmask(:,:).
      !
      ! The logical mask for the coupler must be generated from the fractional mask by setting
      !   only those gridpoints to land that have no ocean part at all (lsf<1 is ocean).
      ! The logical mask is then set to .FALSE. for land points to exclude them from mapping by yac.
      ! These points are not touched by yac.
      !
      ALLOCATE(lsmnolake(nproma,p_patch%nblks_c), STAT = error)
      IF(error /= SUCCESS) &
        CALL finish(routine, "memory allocation failure for lsmnolake")

      SELECT CASE( iforcing ) !{{{

        CASE ( inwp )

          use_mask = .TRUE.

          !ICON_OMP_PARALLEL PRIVATE(jb,jc)
            !ICON_OMP_WORKSHARE
            is_valid(:) = .FALSE.
            !ICON_OMP_END_WORKSHARE

            !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
            DO jb = 1, p_patch%nblks_c
              DO jc = 1, ext_data%atm%list_sea%ncount(jb)
                is_valid((jb-1)*nproma + ext_data%atm%list_sea%idx(jc,jb)) = .TRUE.
              END DO
            END DO
            !ICON_OMP_END_DO
          !ICON_OMP_END_PARALLEL

          CALL dbg_print('AtmFrame: fr_land',ext_data%atm%fr_land,str_module,3,in_subset=p_patch%cells%owned)
          CALL dbg_print('AtmFrame: fr_lake',ext_data%atm%fr_lake,str_module,3,in_subset=p_patch%cells%owned)

        CASE ( iaes )
#ifdef __NO_AES__
          CALL finish (routine, &
              & 'coupled model needs aes; remove --disable-aes and reconfigure')
#else
          !ICON_OMP_PARALLEL_DO PRIVATE(jb,jc) ICON_OMP_RUNTIME_SCHEDULE
          DO jb = 1, p_patch%nblks_c
            DO jc = 1, nproma
              !  slo: caution - lsmask includes alake, must be added to refetch pure lsm:
              lsmnolake(jc, jb) = prm_field(jg)%lsmask(jc,jb) + prm_field(jg)%alake(jc,jb)
            ENDDO
          ENDDO
          !ICON_OMP_END_PARALLEL_DO

          use_mask = &
            p_lor(ANY(lsmnolake(1:nproma,1:p_patch%nblks_c) /= 0.0), p_comm_work)

          !
          ! Define cell_mask_ids(1): all ocean and coastal points are valid
          !   This is the standard for the coupling of atmospheric fields listed below
          !
          IF ( use_mask ) THEN
            !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc) ICON_OMP_RUNTIME_SCHEDULE
            DO jb = 1, p_patch%nblks_c
              DO jc = 1, nproma

                IF ( lsmnolake(jc, jb) .LT. (1.0_wp - eps) ) THEN
                  ! ocean point (fraction of ocean is >0., lsmnolake .lt. 1.) is valid
                  is_valid((jb-1)*nproma+jc) = .TRUE.
                ELSE
                  ! land point (fraction of land is one, no sea water, lsmnolake=1.) is undef
                  is_valid((jb-1)*nproma+jc) = .FALSE.
                ENDIF

              ENDDO
            ENDDO
            !ICON_OMP_END_PARALLEL_DO
          ELSE
            !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc) ICON_OMP_RUNTIME_SCHEDULE
            DO jc = 1, p_patch%nblks_c * nproma
              is_valid(jc) = .TRUE.
            ENDDO
            !ICON_OMP_END_PARALLEL_DO
          ENDIF
#endif
        CASE DEFAULT

          CALL finish ( &
            routine,  'Please mask handling for new forcing in. Thank you!')

      END SELECT !}}}

      DEALLOCATE (lsmnolake, STAT = error)
      IF(error /= SUCCESS) &
        CALL finish(routine, "Deallocation failed for lsmnolake")

      CALL cpl_def_cell_field_mask( &
        routine, grid_id, is_valid, cell_mask_id )

    END SUBROUTINE def_cell_field_mask

  END SUBROUTINE construct_atmo_ocean_coupling

  SUBROUTINE destruct_atmo_ocean_coupling ()

    CHARACTER(LEN=*), PARAMETER :: &
      routine = str_module // ':destruct_atmo_ocean_coupling'

    SELECT CASE( iforcing ) !{{{

      CASE ( inwp )
      CASE ( iaes )
#ifdef __NO_AES__
        CALL finish ( &
          routine, 'coupled model needs aes; remove --disable-aes and reconfigure')
#else
        CALL destruct_aes_ocean_coupling()
#endif
      CASE DEFAULT

        CALL finish ( &
          routine, 'Please mask handling for new forcing. Thank you!')

    END SELECT !}}}

    CALL destruct_atmo_ocean_coupling_common()

  END SUBROUTINE destruct_atmo_ocean_coupling

END MODULE mo_atmo_ocean_coupling
