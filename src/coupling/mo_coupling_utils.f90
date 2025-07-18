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

! Set of routines shared by various coupling related modules

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_coupling_utils

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, message_text, warning, finish
  USE mo_model_domain,    ONLY: t_patch
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info
  USE mo_parallel_config, ONLY: nproma
  USE mo_grid_config,     ONLY: n_dom
  USE mo_run_config,      ONLY: ltimer
  USE mo_master_control,  ONLY: get_my_process_name
  USE mo_time_config,     ONLY: time_config
  USE mo_mpi,             ONLY: p_pe_work
  USE mo_loopindices,     ONLY: get_indices_c
  USE mo_fortran_tools,   ONLY: swap
  USE mtime,              ONLY: datetime, newdatetime, deallocateDatetime, &
    &                           datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_timer,           ONLY: timer_start, timer_stop, &
    &                           timer_coupling_nop, &
    &                           timer_coupling_put_reduce, &
    &                           timer_coupling_put, &
    &                           timer_coupling_get, &
    &                           timer_coupling_very_1stget, &
    &                           timer_coupling_1stget, &
    &                           timer_coupling_init, &
    &                           timer_coupling_init_def_comp, &
    &                           timer_coupling_init_enddef
  USE mo_impl_constants,  ONLY: MAX_CHAR_LENGTH
  USE mo_util_table,      ONLY: t_table, initialize_table, add_table_column, &
    &                           set_table_entry, print_table, finalize_table
#ifdef YAC_coupling
  USE mo_mpi,             ONLY: p_comm_yac, p_comm_work, my_process_is_stdio
  USE yac,                ONLY: yac_finit, yac_finit_comm, &
    &                           yac_fread_config_yaml, &
    &                           yac_ffinalize, yac_fdef_comp, &
    &                           yac_fdef_comps, yac_fdef_grid, &
    &                           yac_fget_version, yac_fdef_mask, &
    &                           yac_fdef_points, yac_fdef_field, &
    &                           yac_fdef_datetime, yac_fdef_field_mask, &
    &                           yac_fset_global_index, yac_fset_core_mask, &
    &                           yac_fget_action, yac_fupdate, &
    &                           yac_dble_ptr, yac_fput, yac_fget, &
    &                           yac_fget_field_collection_size, &
    &                           yac_fget_collection_size_from_field_id, &
    &                           yac_fget_field_datetime, &
    &                           yac_fsync_def, yac_fenddef, &
    &                           yac_fget_grid_size, &
    &                           YAC_LOCATION_CELL, &
    &                           YAC_LOCATION_CORNER, &
    &                           YAC_LOCATION_EDGE, &
    &                           YAC_TIME_UNIT_ISO_FORMAT, &
    &                           YAC_ACTION_NONE, &
    &                           YAC_ACTION_REDUCTION, &
    &                           YAC_ACTION_COUPLING, &
    &                           YAC_ACTION_PUT_FOR_RESTART, &
    &                           YAC_ACTION_GET_FOR_RESTART, &
    &                           YAC_ACTION_OUT_OF_BOUND, &
    &                           yac_string, &
    &                           yac_fget_comp_names, yac_fget_grid_names, yac_fget_field_names, &
    &                           yac_fget_field_id, yac_fget_field_timestep
  USE mpi
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cpl_construct
  PUBLIC :: cpl_destruct
  PUBLIC :: cpl_is_initialised
  PUBLIC :: cpl_get_instance_id
  PUBLIC :: cpl_config_file_exists
  PUBLIC :: cpl_def_main
  PUBLIC :: cpl_def_main_dummy
  PUBLIC :: cpl_def_cell_field_mask
  PUBLIC :: cpl_def_field
  PUBLIC :: cpl_get_field
  PUBLIC :: cpl_get_field_collection_size
  PUBLIC :: cpl_get_field_datetime
  PUBLIC :: cpl_put_field
  PUBLIC :: cpl_sync_def
  PUBLIC :: cpl_enddef
  PUBLIC :: cpl_write_config_info

  CHARACTER(LEN=*), PARAMETER :: yaml_filename = "coupling.yaml"

  LOGICAL :: yac_is_initialised = .FALSE.
  INTEGER :: yac_instance_id = -1

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_coupling_utils'

  ! register the main component (and optionally the output component)
  INTERFACE cpl_def_main
    MODULE PROCEDURE cpl_def_main_without_output
    MODULE PROCEDURE cpl_def_main_with_output
  END INTERFACE cpl_def_main

  ! registers a field to the coupler
  INTERFACE cpl_def_field
    MODULE PROCEDURE cpl_def_field_no_mask
    MODULE PROCEDURE cpl_def_field_mask
  END INTERFACE cpl_def_field

  ! receives cell-based field data through the coupler
  INTERFACE cpl_get_field
    MODULE PROCEDURE cpl_get_field_idx_lev_blk
    MODULE PROCEDURE cpl_get_field_idx_blk_collection
    MODULE PROCEDURE cpl_get_field_n_collection
  END INTERFACE cpl_get_field

  ! sends cell-based field data through the coupler
  INTERFACE cpl_put_field
    MODULE PROCEDURE cpl_put_field_idx_blk_collection
    MODULE PROCEDURE cpl_put_field_n_collection
  END INTERFACE cpl_put_field

  INTERFACE cpl_get_field_collection_size
    MODULE PROCEDURE get_field_collection_size_from_name
    MODULE PROCEDURE get_field_collection_size_from_id
  END INTERFACE cpl_get_field_collection_size

  TYPE logical_2d_ptr
    LOGICAL, ALLOCATABLE :: p(:,:)
  END TYPE logical_2d_ptr

CONTAINS

  SUBROUTINE cpl_construct()

#if !defined NOMPI && defined YAC_coupling

    INTEGER :: global_rank, ierror

    IF (ltimer) CALL timer_start (timer_coupling_init)

    yac_is_initialised = .TRUE.

    CALL yac_finit_comm ( p_comm_yac, yac_instance_id )

    IF (ltimer) CALL timer_stop(timer_coupling_init)

#endif

  END SUBROUTINE cpl_construct

  SUBROUTINE cpl_destruct()

#if !defined NOMPI && defined YAC_coupling
    IF (yac_is_initialised) CALL yac_ffinalize( yac_instance_id )
#endif

  END SUBROUTINE cpl_destruct

  FUNCTION cpl_config_file_exists()

    LOGICAL :: cpl_config_file_exists

    LOGICAL, SAVE :: config_files_exist = .FALSE.
    LOGICAL, SAVE :: config_files_have_been_checked = .FALSE.

    LOGICAL :: yaml_exists

    IF (config_files_have_been_checked) THEN

      cpl_config_file_exists = config_files_exist

    ELSE

      INQUIRE(FILE=TRIM(ADJUSTL(yaml_filename)), EXIST=yaml_exists)

      config_files_have_been_checked = .TRUE.
      config_files_exist = yaml_exists
      cpl_config_file_exists = config_files_exist

    END IF

  END FUNCTION cpl_config_file_exists

  FUNCTION cpl_is_initialised()

    LOGICAL :: cpl_is_initialised

    cpl_is_initialised = yac_is_initialised

  END FUNCTION cpl_is_initialised

  FUNCTION cpl_get_instance_id()

    INTEGER :: cpl_get_instance_id

    CHARACTER(*), PARAMETER :: &
      routine = modname // ":cpl_get_instance_id"

    IF (.NOT. yac_is_initialised) &
      CALL finish(routine, "YAC has not been initialised")

    cpl_get_instance_id = yac_instance_id

  END FUNCTION cpl_get_instance_id


  !
  ! Write information about coupler configuration to stdout
  !
  SUBROUTINE cpl_write_config_info(caller, my_comp_name, my_grid_name)
    CHARACTER(LEN=*), INTENT(IN) :: caller             ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: my_comp_name       ! component for which the config info is requested
    CHARACTER(LEN=*), INTENT(IN) :: my_grid_name       ! grid for which the config info is requested

#ifdef YAC_coupling
    ! local
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':cpl_write_config_info'

    INTEGER :: nbr_comp, nbr_grid, nbr_field
    INTEGER :: i, my_comp_id, my_grid_id
    INTEGER :: field_id
    CHARACTER(len=:), ALLOCATABLE:: field_datetime ! coupling startdate
    CHARACTER(len=:), ALLOCATABLE:: field_timestep
    INTEGER, PARAMETER :: UNDEF = -1
    TYPE(yac_string), ALLOCATABLE :: comp_names(:)
    TYPE(yac_string), ALLOCATABLE :: grid_names(:)
    TYPE(yac_string), ALLOCATABLE :: field_names(:)

    ! get names of all registered components
    comp_names = yac_fget_comp_names(yac_instance_id)
    nbr_comp = SIZE(comp_names)
    ! get names of all registered grids
    grid_names = yac_fget_grid_names(yac_instance_id)
    nbr_grid = SIZE(grid_names)

    ! check availability of my component
    my_comp_id = UNDEF
    DO i=1,nbr_comp
      IF (TRIM(comp_names(i)%string) == TRIM(my_comp_name)) THEN
        my_comp_id = i
        EXIT
      ENDIF
    ENDDO
    IF (my_comp_id==UNDEF) CALL finish(routine, "No matching component found")

    ! check availability of my grid
    my_grid_id = UNDEF
    DO i=1,nbr_grid
      IF (TRIM(grid_names(i)%string) == TRIM(my_grid_name)) THEN
        my_grid_id = i
        EXIT
      ENDIF
    ENDDO
    IF (my_grid_id==UNDEF) CALL finish(routine, "No matching grid name found")

    ! get all registered fields for my component and my grid
    field_names = yac_fget_field_names(yac_instance_id, TRIM(my_comp_name), TRIM(my_grid_name))
    nbr_field = SIZE(field_names)

    ! In the following we assume that all registered coupling fields have the same
    ! field datetime. This is ensured by the routine construct_X_Y_coupling_finalize,
    ! stored in the module mo_X_Y_coupling. Hence, it is sufficient at this point
    ! to request the datetime information for a single field.
    field_id = yac_fget_field_id(yac_instance_id, TRIM(my_comp_name), TRIM(my_grid_name), &
      &                          TRIM(field_names(1)%string))

    ! get coupling startdate
    field_datetime = yac_fget_field_datetime(field_id)
    ! get field timestep (not to be confused with the coupling period!)
    field_timestep = yac_fget_field_timestep(field_id)

    ! write configuration details to stdout
    CALL cpl_write_config()

    CONTAINS

    !
    ! Write YAC configuration details to stdout for given component
    !
    SUBROUTINE cpl_write_config

      TYPE(t_table) :: table
      CHARACTER(LEN=*), PARAMETER :: colAtt  = 'Attribute',  &
        &                            colVal  = 'Value'
      INTEGER :: irow

      ! will only be executed by stdio process
      IF(.NOT. my_process_is_stdio()) RETURN

      ! poor man's table header
      WRITE (0,*) " " ! newline
      WRITE(message_text,'(a,a)') 'YAC configuration details for ', TRIM(my_comp_name)
      CALL message('', message_text)

      ! table-based output
      CALL initialize_table(table)
      CALL add_table_column(table, colAtt)
      CALL add_table_column(table, colVal)


      irow = 1
      CALL set_table_entry(table, irow, colAtt, 'Components')
      CALL set_table_entry(table, irow, colVal, TRIM(comp_names(1)%string))
      !
      DO i=2,nbr_comp
        irow = irow+1
        CALL set_table_entry(table, irow, colAtt, ' ')
        CALL set_table_entry(table, irow, colVal, TRIM(comp_names(i)%string))
      ENDDO

      irow = irow+1
      CALL set_table_entry(table, irow, colAtt, 'Grids')
      CALL set_table_entry(table, irow, colVal, TRIM(grid_names(1)%string))
      !
      DO i=2,nbr_grid
        irow = irow+1
        CALL set_table_entry(table, irow, colAtt, ' ')
        CALL set_table_entry(table, irow, colVal, TRIM(grid_names(i)%string))
      ENDDO

      irow = irow+1
      CALL set_table_entry(table, irow, colAtt, 'My component')
      CALL set_table_entry(table, irow, colVal, TRIM(my_comp_name))

      irow = irow+1
      CALL set_table_entry(table, irow, colAtt, 'My grid')
      CALL set_table_entry(table, irow, colVal, TRIM(my_grid_name))

      irow = irow+1
      CALL set_table_entry(table, irow, colAtt, 'Coupling fields')
      CALL set_table_entry(table, irow, colVal, TRIM(field_names(1)%string))
      !
      DO i=2,nbr_field
        irow = irow+1
        CALL set_table_entry(table, irow, colAtt, ' ')
        CALL set_table_entry(table, irow, colVal, TRIM(field_names(i)%string))
      ENDDO

      irow = irow+1
      CALL set_table_entry(table, irow, colAtt, 'Coupling startdate')
      CALL set_table_entry(table, irow, colVal, TRIM(field_datetime))

      irow = irow+1
      CALL set_table_entry(table, irow, colAtt, 'Field timestep')
      CALL set_table_entry(table, irow, colVal, TRIM(field_timestep))

      CALL print_table(table, opt_delimiter=' | ')
      CALL finalize_table(table)

      WRITE (0,*) " " ! newline

    END SUBROUTINE cpl_write_config

! YAC_coupling
#endif
  END SUBROUTINE cpl_write_config_info


  SUBROUTINE def_patch( &
    p_patch, grid_name, grid_id, cell_point_id, vertex_point_id)

    TYPE(t_patch), INTENT(IN) :: p_patch              ! basic patch
    CHARACTER(LEN=*), INTENT(IN) :: grid_name         ! name of the grid
    INTEGER, INTENT(OUT) :: grid_id                   ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id             ! cell coordinate id
    INTEGER, OPTIONAL, INTENT(OUT) :: vertex_point_id ! vertex coordinate id (only with output)

    INTEGER :: jb, jc, jv, nblks, nn

    REAL(wp), ALLOCATABLE :: buffer_lon(:)
    REAL(wp), ALLOCATABLE :: buffer_lat(:)
    INTEGER,  ALLOCATABLE :: buffer_c(:,:)

    LOGICAL,  ALLOCATABLE :: is_valid(:)

#ifdef YAC_coupling

    ! Extract cell information
    !
    ! cartesian coordinates of cell vertices are stored in
    ! p_patch%verts%cartesian(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes in rad.

    nblks = MAX(p_patch%nblks_c,p_patch%nblks_v)

    ALLOCATE(                     &
      buffer_lon(nproma * nblks), &
      buffer_lat(nproma * nblks), &
      buffer_c(3, nproma * nblks))

!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, jv, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, p_patch%nblks_v
      DO jv = 1, nproma
        nn = (jb-1)*nproma+jv
        buffer_lon(nn) = p_patch%verts%vertex(jv,jb)%lon
        buffer_lat(nn) = p_patch%verts%vertex(jv,jb)%lat
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

!ICON_OMP_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, p_patch%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_c(1,nn) = (p_patch%cells%vertex_blk(jc,jb,1)-1)*nproma + &
          &               p_patch%cells%vertex_idx(jc,jb,1)
        buffer_c(2,nn) = (p_patch%cells%vertex_blk(jc,jb,2)-1)*nproma + &
          &               p_patch%cells%vertex_idx(jc,jb,2)
        buffer_c(3,nn) = (p_patch%cells%vertex_blk(jc,jb,3)-1)*nproma + &
                          p_patch%cells%vertex_idx(jc,jb,3)
      ENDDO
    ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! Definition of unstructured horizontal grid
    CALL yac_fdef_grid(                                &
      & grid_name             = TRIM(grid_name),       & !in
      & nbr_vertices          = p_patch%n_patch_verts, & !in
      & nbr_cells             = p_patch%n_patch_cells, & !in
      & nbr_vertices_per_cell = 3,                     & !in
      & x_vertices            = buffer_lon,            & !in
      & y_vertices            = buffer_lat,            & !in
      & cell_to_vertex        = buffer_c,              & !in
      & grid_id               = grid_id)                 !out

    IF (PRESENT(vertex_point_id)) THEN

      ! the output component also defines fields located at the vertices
      CALL yac_fdef_points(    &
        grid_id,               & !in
        p_patch%n_patch_verts, & !in
        YAC_LOCATION_CORNER,   & !in
        buffer_lon,            & !in
        buffer_lat,            & !in
        vertex_point_id)         !out
    END IF

    ! Define cell center points
    !
    ! cartesian coordinates of cell centers are stored in
    ! p_patch%cells%cartesian_center(:,:)%x(1:3)
    ! Here we use the longitudes and latitudes.

    !ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = 1, p_patch%nblks_c
      DO jc = 1, nproma
        nn = (jb-1)*nproma+jc
        buffer_lon(nn) = p_patch%cells%center(jc,jb)%lon
        buffer_lat(nn) = p_patch%cells%center(jc,jb)%lat
      ENDDO
    ENDDO
    !ICON_OMP_END_PARALLEL_DO

    ! Definition of cell center points
    CALL yac_fdef_points (     &
      & grid_id,               & !in
      & p_patch%n_patch_cells, & !in
      & YAC_LOCATION_CELL,     & !in
      & buffer_lon,            & !in
      & buffer_lat,            & !in
      & cell_point_id )          !out

    DEALLOCATE (buffer_lon, buffer_lat, buffer_c)

    nblks = &
      MAX(p_patch%nblks_c, p_patch%nblks_v)
    ALLOCATE(is_valid(nproma*nblks))

    ! set global indices and core masks
    CALL set_basic_info( &
      p_patch%cells%decomp_info, YAC_LOCATION_CELL, .TRUE.)
    CALL set_basic_info( &
      p_patch%verts%decomp_info, YAC_LOCATION_CORNER, .FALSE.)
!     CALL set_basic_info( &
!       p_patch%edges%decomp_info, YAC_LOCATION_EDGE)

    DEALLOCATE (is_valid)

  CONTAINS

    SUBROUTINE set_basic_info(decomp_info, location, set_core_mask)

      TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info
      INTEGER, INTENT(IN) :: location
      LOGICAL, INTENT(IN) :: set_core_mask

      INTEGER :: i, point_count

      point_count = SIZE(decomp_info%glb_index)

      CALL yac_fset_global_index ( &
        decomp_info%glb_index, location, grid_id)

      ! Generate core mask
!ICON_OMP_PARALLEL_DO PRIVATE(i) ICON_OMP_RUNTIME_SCHEDULE
      DO i = 1, point_count
        is_valid(i) = p_pe_work == decomp_info%owner_local(i)
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      IF (set_core_mask) THEN
        ! Define core mask
        CALL yac_fset_core_mask ( is_valid, location, grid_id )
      END IF

    END SUBROUTINE set_basic_info

! YAC_coupling
#endif

  END SUBROUTINE def_patch

  !----------------------------------------------------------------------------
  ! This routine generates a grid which is a combination of the main patch
  ! and all nested patches. It can be used to provide the collective data from
  ! all patches to other components.
  !
  ! Detailed explanation:
  !
  ! The following represents a main patch with:
  !
  !                             15
  !                             / \
  !                            /   \
  !                          29    30
  !                          /  16   \
  !                         /         \
  !                       13----28----14
  !                       / \         / \
  !                      /   \  14   /   \
  !                    24    25    26    27
  !                    /  13   \   /  15   \
  !                   /         \ /         \
  !                 10----22----11----23----12
  !                 / \         / \         / \
  !                /   \   9   /   \   11  /   \
  !              16    17    18    19    20    21
  !              /   8   \   /  10   \   /  12   \
  !             /         \ /         \ /         \
  !            6-----13----7----14-----8-----15----9
  !           / \         / \         / \         / \
  !          /   \   2   /   \   4   /   \   6   /   \
  !         5     6     7     8     9    10    11    12
  !        /   1   \   /   3   \   /   5   \   /   7   \
  !       /         \ /         \ /         \ /         \
  !      1-----1-----2-----2-----3-----3-----4-----4-----5
  !
  ! In addition to the main patch a nested patch is defined, which exacly
  ! covers the cell 10, 11, 12, and 15 of the main patch.
  ! A nested patch consists of a prognostic core and 4 rows of boundary cells
  ! around it. The boundary cells are ignored, as they are covered by the parent patch.
  !
  !                             15
  !                             / \
  !                            /   \
  !                          13----14
  !                          / \   / \
  !                         /   \ /   \
  !                       10----11----12
  !                       / \   / \   / \
  !                      /   \ /   \ /   \
  !                     6-----7-----8-----9
  !                    / \   / \   / \   / \
  !                   /   \ /   \ /   \ /   \
  !                  1-----2-----3-----4-----5
  !
  ! For parallelisation all patches are independently distributed among the
  ! processes of a component. Therefore, if cells from one patch are covered
  ! by a nested patch, the cells of both patches are not necessarly on the same
  ! process.
  !
  ! In the example we assume that there is only a single process.
  !----------------------------------------------------------------------------
#ifdef HAVE_YAXT
  SUBROUTINE def_patch_combined(p_patch, grid_name, grid_id, cell_point_id)

    USE mo_model_domain,       ONLY: p_patch_local_parent
    USE mo_math_types,         ONLY: t_cartesian_coordinates
    USE mo_impl_constants,     ONLY: start_prog_cells, end_prog_cells
    USE mo_impl_constants_grf, ONLY: grf_bdyintp_end_c, grf_bdyintp_end_e
    USE mo_parallel_config,    ONLY: nproma, idx_no, blk_no, idx_1d
    USE yaxt,                  ONLY: xt_idxlist, xt_xmap, xt_redist, &
                                     xt_idxvec_new, xt_idxlist_delete, &
                                     xt_xmap_dist_dir_new, xt_xmap_delete, &
                                     xt_redist_p2p_new, xt_redist_p2p_off_new, &
                                     xt_redist_collection_new, &
                                     xt_redist_s_exchange, &
                                     xt_redist_s_exchange1, xt_redist_delete
    USE iso_c_binding,         ONLY: C_LOC, C_NULL_PTR, C_PTR

    TYPE(t_patch), INTENT(IN) :: p_patch(1:n_dom) ! array of nested patches
    CHARACTER(LEN=*), INTENT(IN) :: grid_name     ! name of the nested grid
    INTEGER, INTENT(OUT) :: grid_id               ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id         ! cell coordinate id

    TYPE(logical_2d_ptr), ALLOCATABLE :: dom_cell_core_mask(:)
    TYPE(logical_2d_ptr), ALLOCATABLE :: dom_vertex_child_mask(:)
    TYPE(logical_2d_ptr), ALLOCATABLE :: dom_edge_child_mask(:)

    INTEGER, ALLOCATABLE :: nbr_vertices_per_cell(:)
    REAL(wp), ALLOCATABLE, TARGET :: vertex_lon(:)
    REAL(wp), ALLOCATABLE, TARGET :: vertex_lat(:)
    REAL(wp), ALLOCATABLE :: cell_lon(:)
    REAL(wp), ALLOCATABLE :: cell_lat(:)
    INTEGER,  ALLOCATABLE :: cell_to_vertex(:)
    INTEGER,  ALLOCATABLE :: cell_global_ids(:)
    INTEGER, ALLOCATABLE, TARGET :: vertex_global_ids(:)

    LOGICAL,  ALLOCATABLE :: is_valid_cell(:)

    INTEGER :: total_num_cells, total_num_verts, total_num_cell_blk
    INTEGER :: cells_offset, verts_offset
    INTEGER :: global_id_offset_c, global_id_offset_v, global_id_offset_e
    INTEGER :: cell_to_vertex_offset
    INTEGER :: nest_vert_offset
    INTEGER, ALLOCATABLE:: parent_global_id_offset_v(:)
    INTEGER, ALLOCATABLE :: parent_global_id_offset_e(:)

    INTEGER :: dom_id, parent_dom_id, i
    INTEGER :: i_c, i_v, i_e, j_c, j_v, j_e
    INTEGER :: jc_c, jb_c, jc_v, jb_v, jc_e, jb_e
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    LOGICAL :: cell_is_valid
    INTEGER :: curr_num_nest_edges
    INTEGER :: start_vertex_idx, start_vertex_blk
    INTEGER :: end_vertex_idx, end_vertex_blk

    INTEGER :: edge_idx(2), edge_blk(2)
    INTEGER :: vertex_idx(4), vertex_blk(4), vertex_1d(4)
    INTEGER :: vertex_order(3)
    TYPE(t_cartesian_coordinates) :: vertex_coord
    INTEGER :: nest_edge_idx(3)
    INTEGER :: num_nest_verts, num_nest_edges

    INTEGER :: local_parent_edge_idx, &
               local_parent_edge_blk, &
               local_parent_edge_vertex_idx(2), &
               local_parent_edge_vertex_blk(2), &
               local_parent_edge_vertex_1d(2), &
               local_parent_edge_vertex_global_ids(2)
    TYPE(t_cartesian_coordinates) :: local_parent_edge_vertex_coord(2)
    INTEGER, ALLOCATABLE :: local_parent_vertex_global_ids(:), &
                            local_parent_edge_global_ids(:), &
                            local_parent_vertex_child_offsets(:), &
                            local_parent_edge_middle_vertex_offsets(:)
    INTEGER :: total_num_local_parent_edges, &
               total_num_local_parent_verts

    INTEGER, ALLOCATABLE :: required_vertices(:), &
                            required_vertices_offsets(:), &
                            required_edges(:), &
                            required_edge_middle_vertex_offsets(:)

    TYPE(xt_idxlist) :: src_idxlist, tgt_idxlist
    TYPE(xt_xmap) :: xmap_edge_vertex_data, xmap_vertex_data
    TYPE(xt_redist) :: redist_edge_vertex_data_int, &
                       redist_edge_vertex_data_dble, &
                       redist_vertex_data_int, &
                       redist_vertex_data_dble, &
                       redist_vertex_data
    TYPE(C_PTR) :: data_ptr(6)
    INTEGER, TARGET :: dummy

    CHARACTER(LEN=*), PARAMETER :: method_name = "def_patch_combined"

#ifdef YAC_coupling

    !--------------------------------------------------------------------------
    ! Count the total number of cells and vertices in all patches on each
    ! process
    !--------------------------------------------------------------------------

    total_num_cells = 0
    total_num_cell_blk = 0
    total_num_verts = 0
    DO dom_id = 1, n_dom
      total_num_cells = &
        total_num_cells + p_patch(dom_id)%n_patch_cells
      total_num_cell_blk = &
        total_num_cell_blk + p_patch(dom_id)%nblks_c
      total_num_verts = &
        total_num_verts + p_patch(dom_id)%n_patch_verts
    END DO

    !--------------------------------------------------------------------------
    ! Generate cell core mask
    ! (deactivate all halo cells and cell that
    !  have children)
    !
    ! In the example the cells 10, 11, 12, and 15
    ! would have to be masked out:
    !                             15
    !                             / \
    !                            /   \
    !                           /     \
    !                          /       \
    !                         /         \
    !                       13----------14
    !                       / \         / \
    !                      /   \       /   \
    !                     /     \     /     \
    !                    /       \   /  xx   \
    !                   /         \ /         \
    !                 10----------11----------12
    !                 / \         / \         / \
    !                /   \       /   \   xx  /   \
    !               /     \     /     \     /     \
    !              /       \   /  xx   \   /  xx   \
    !             /         \ /         \ /         \
    !            6-----------7-----------8-----------9
    !           / \         / \         / \         / \
    !          /   \       /   \       /   \       /   \
    !         /     \     /     \     /     \     /     \
    !        /       \   /       \   /       \   /       \
    !       /         \ /         \ /         \ /         \
    !      1-----------2-----------3-----------4-----------5
    !
    ! These cells will be marked out in the coarse grid, because the data in
    ! these parts will be provided by the nests.
    !
    ! dom_cell_core_mask: contains a core mask for each patch
    ! is_valid_cell: combined core mask for all patches
    !--------------------------------------------------------------------------

    ALLOCATE(dom_cell_core_mask(n_dom))
    ALLOCATE(is_valid_cell(total_num_cells))
    i = 0
    is_valid_cell = .FALSE.

    ! for all domains
    DO dom_id = 1, n_dom
      ALLOCATE( &
        dom_cell_core_mask(dom_id)%p(nproma, p_patch(dom_id)%nblks_c))
      dom_cell_core_mask(dom_id)%p = .FALSE.

      ! prognostic domain
      i_startblk = p_patch(dom_id)%cells%start_block(start_prog_cells)
      i_endblk   = p_patch(dom_id)%cells%end_block(end_prog_cells)

      DO jb_c = i_startblk, i_endblk

        CALL get_indices_c( &
          p_patch(dom_id), jb_c, i_startblk, i_endblk, &
          i_startidx, i_endidx, start_prog_cells, end_prog_cells)

        DO jc_c = i_startidx, i_endidx

          ! determine if current cell is owned by the local process
          ! and has no children,  If it corresponds to a boundary
          ! child cell (refin_ctrl = -1, -2) then keep it.
          cell_is_valid = &
            (p_patch(dom_id)%cells%decomp_info%decomp_domain(jc_c, jb_c) &
             == 0) .AND. &
            ((p_patch(dom_id)%cells%refin_ctrl(jc_c, jb_c) >= grf_bdyintp_end_c))
          dom_cell_core_mask(dom_id)%p(jc_c, jb_c) = cell_is_valid
          is_valid_cell(i + idx_1d(jc_c, jb_c)) = cell_is_valid
        END DO
      END DO

      i = i + p_patch(dom_id)%n_patch_cells
    END DO

    !--------------------------------------------------------------------------
    ! Determine which edges and vertices of a patch are located on the border
    ! of a nest. In addition, they are counted for all patches.
    !
    ! In the example only the main domain has a nested patch. The respective
    ! edge and vertices are marked with "x"
    !
    !                             15
    !                             / \
    !                            /   \
    !                           /     \
    !                          /  16   \
    !                         /         \
    !                       13-----------x
    !                       / \         x \
    !                      /   \  14   x   \
    !                     /     \     x     \
    !                    /  13   \   x  15   \
    !                   /         \ x         \
    !                 10-----------x----------12
    !                 / \         x \         / \
    !                /   \   9   x   \   11  /   \
    !               /     \     x     \     /     \
    !              /   8   \   x  10   \   /  12   \
    !             /         \ x         \ /         \
    !            6-----------xxxxxxxxxxxxxxxxxxxxxxxxx
    !           / \         / \         / \         / \
    !          /   \   2   /   \   4   /   \   6   /   \
    !         /     \     /     \     /     \     /     \
    !        /   1   \   /   3   \   /   5   \   /   7   \
    !       /         \ /         \ /         \ /         \
    !      1-----------2-----------3-----------4-----------5
    !
    ! In the coarse grid, these vertices and edges require special handling
    ! in order to have a consistent grid across these borders from the coarse
    ! to the finer patch. Therefore, the vertices in the coarse grid that are
    ! on the border will be replaced (coordinates and global ids) by the ones
    ! in the nest. In the middle of the border edge an additional vertex will
    ! be introduced.
    !
    ! In the base patch cell 4 consists of the vertices 3, 8, and 7. In the
    ! combined grid it will consists of 3, 3*, 2*, and 1*
    ! (*: vertices of the nest).
    !
    ! dom_vertex_child_mask: contains for the vertices of each patch a mask
    !                        marking the respective vertices
    ! dom_edge_child_mask: contains for the edge of each patch a mask
    !                      marking the respective edges
    !--------------------------------------------------------------------------

    num_nest_edges = 0
    num_nest_verts = 0
    ALLOCATE(dom_vertex_child_mask(n_dom), dom_edge_child_mask(n_dom))

    ! for all domains
    DO dom_id = 1, n_dom

      ALLOCATE( &
        dom_vertex_child_mask(dom_id)%p(nproma, p_patch(dom_id)%nblks_v), &
        dom_edge_child_mask(dom_id)%p(nproma, p_patch(dom_id)%nblks_e))
      dom_vertex_child_mask(dom_id)%p = .FALSE.
      dom_edge_child_mask(dom_id)%p = .FALSE.

      ! if the current domain has nests
      IF (p_patch(dom_id)%n_childdom > 0) THEN

        ! for all cell of the current domain
        DO jb_c = 1, p_patch(dom_id)%nblks_c
          DO jc_c = 1, nproma

            ! if the current cell is active
            IF (dom_cell_core_mask(dom_id)%p(jc_c, jb_c)) THEN

              ! for all edges of the current cell
              DO j_e = 1, 3

                jc_e = p_patch(dom_id)%cells%edge_idx(jc_c, jb_c, j_e)
                jb_e = p_patch(dom_id)%cells%edge_blk(jc_c, jb_c, j_e)

                ! determine whether this edge has children
                dom_edge_child_mask(dom_id)%p(jc_e, jb_e) = &
                  p_patch(dom_id)%edges%refin_ctrl(jc_e, jb_e) == grf_bdyintp_end_e

              END DO
            END IF
          END DO
        END DO

        ! for all cell of the current domain
        DO jb_c = 1, p_patch(dom_id)%nblks_c
          DO jc_c = 1, nproma

            ! if the current cell is active
            IF (dom_cell_core_mask(dom_id)%p(jc_c, jb_c)) THEN

              ! for all vertices of the current cell
              DO j_v = 1, 3

                jc_v = p_patch(dom_id)%cells%vertex_idx(jc_c, jb_c, j_v)
                jb_v = p_patch(dom_id)%cells%vertex_blk(jc_c, jb_c, j_v)

                ! determine whether this vertex has children
                dom_vertex_child_mask(dom_id)%p(jc_v, jb_v) = &
                  p_patch(dom_id)%verts%refin_ctrl(jc_v, jb_v) == &
                  grf_bdyintp_end_c - 1

              END DO
            END IF
          END DO
        END DO

        num_nest_edges = &
          num_nest_edges + COUNT(dom_edge_child_mask(dom_id)%p)
        num_nest_verts = &
          num_nest_verts + COUNT(dom_vertex_child_mask(dom_id)%p)
      END IF
    END DO

    !--------------------------------------------------------------------------
    ! Get global ids of all nest edges and generate an offset which is the
    ! total number of edges in all patches on the local process plus i - 1
    ! (where i it the i'th nested edge; "- 1" because offsets are "0"-based).
    !
    ! For each patch (except for the base patch) there is a local parent patch,
    ! which coveres the same area as the respective patch, but has the
    ! resolution of the parent patch. The number of edges in all local
    ! parent patches is computed here.
    !
    ! These edges will have to split into two edges with an additional vertex
    ! in the middle for which we will require the coordinates and the global id
    ! from the nested patch.
    !
    ! The global ids for cells, vertices, and edges for each patch start of "1"
    ! for each patch. In order to make them unique across patches, the id of
    ! each patch are offset by the global number of cells/vertices/edges of
    ! all patches with lower domain ids.
    !
    !
    ! In the example:
    ! required_edges(1:4) = (/14,15,18,26/)
    ! required_edge_middle_vertex_offsets(1:4) = (/60,61,62,63/)
    ! total_num_local_parent_edges = 9
    !--------------------------------------------------------------------------

    ALLOCATE(required_edges(num_nest_edges), &
             required_edge_middle_vertex_offsets(num_nest_edges))
    i_e = 0
    global_id_offset_e = 0
    total_num_local_parent_edges = 0

    ! for all domains
    DO dom_id = 1, n_dom

      ! for all cells of the current domain
      DO jb_c = 1, p_patch(dom_id)%nblks_c
        DO jc_c = 1, nproma

          ! if the current cell is active
          IF (dom_cell_core_mask(dom_id)%p(jc_c, jb_c)) THEN

            ! for all edges of the current cell
            DO j_e = 1, 3

              jc_e = p_patch(dom_id)%cells%edge_idx(jc_c, jb_c, j_e)
              jb_e = p_patch(dom_id)%cells%edge_blk(jc_c, jb_c, j_e)

              ! determine whether this edge has children
              IF (dom_edge_child_mask(dom_id)%p(jc_e, jb_e)) THEN
                i_e = i_e + 1
                required_edges(i_e) = &
                  p_patch(dom_id)%edges%decomp_info%glb_index( &
                    idx_1d(jc_e, jb_e)) + global_id_offset_e
                required_edge_middle_vertex_offsets(i_e) = &
                  total_num_verts + i_e - 1

              END IF
            END DO
          END IF
        END DO
      END DO

      global_id_offset_e = &
        global_id_offset_e + p_patch(dom_id)%n_patch_edges_g

      ! if the current domain is not the first patch, which has
      ! not local parent patch
      IF (dom_id /= 1) THEN
        total_num_local_parent_edges = &
          total_num_local_parent_edges + &
          p_patch_local_parent(dom_id)%n_patch_edges
      END IF
    END DO

    !--------------------------------------------------------------------------
    ! Get global ids of all nest vertices and generate an offset which is the
    ! offset of the respective vertex in the patch plus the total number
    ! vertices in all patches with lower domain ids.
    !
    ! The number of vertices in all local parent patches is computed here.
    !
    ! These vertices will have to be replaced (coordinates and global id) by
    ! their counterpart in the nested patch.
    !
    ! In the example:
    ! required_vertices(1:5) = (/7,8,9,11,14/)
    ! required_vertices(1:5) = (/6,7,8,10,14/)
    ! total_num_local_parent_verts = 6
    !--------------------------------------------------------------------------

    ALLOCATE(required_vertices(num_nest_verts), &
             required_vertices_offsets(num_nest_verts))
    i = 0
    global_id_offset_v = 0
    total_num_local_parent_verts = 0
    verts_offset = 0

    ! for all domains
    DO dom_id = 1, n_dom

      ! for all vertices of the current domain
      DO jb_v = 1, p_patch(dom_id)%nblks_v
        DO jc_v = 1, nproma

          ! determine whether this vertex has children
          IF (dom_vertex_child_mask(dom_id)%p(jc_v, jb_v)) THEN

            i = i + 1

            required_vertices(i) = &
              p_patch(dom_id)%verts%decomp_info%glb_index( &
                idx_1d(jc_v, jb_v)) + global_id_offset_v
            required_vertices_offsets(i) = &
              verts_offset + idx_1d(jc_v, jb_v) - 1

          END IF
        END DO
      END DO

      global_id_offset_v = &
        global_id_offset_v + p_patch(dom_id)%n_patch_verts_g
      verts_offset = verts_offset + p_patch(dom_id)%n_patch_verts

      ! if the current domain is not the first patch, which has
      ! not local parent patch
      IF (dom_id /= 1) THEN
        total_num_local_parent_verts = &
          total_num_local_parent_verts + &
          p_patch_local_parent(dom_id)%n_patch_verts
      END IF
    END DO

    DO dom_id = 1, n_dom
      DEALLOCATE(dom_vertex_child_mask(dom_id)%p)
    END DO
    DEALLOCATE(dom_vertex_child_mask)

    !--------------------------------------------------------------------------
    ! For all vertices that have a matching vertex in a parent patch, get the
    ! global id of the associated parent vertex and the offset of the
    ! respective vertex in the combined vertices of all patches.
    ! For all vertices that are in the middle of an edge in a matching parent
    ! patch, get the global id of the parent edge and the offset of the
    ! respective vertex in the combined vertices of all patches.
    !
    ! Instead of directly determining the child data that is to be used to
    ! replace border edge vertices and introduce the addition edge middle
    ! points, it is easier to generate a list of all potential candiates on
    ! all processes.
    !
    ! In the example:
    ! local_parent_vertex_global_ids(1:6) = (/7,8,9, 11,12, 15/)
    ! local_parent_vertex_child_offsets(1:6) = (/14,16,18, 23,25, 28/)
    ! local_parent_edge_global_ids(1:9) = (/14,15, 18,19,20,21, 23, 26,27/)
    ! local_parent_edge_middle_vertex_offsets(1:9) =
    !   (/15,17, 19,20,21,22, 24, 26,27/)
    !--------------------------------------------------------------------------

    ALLOCATE(local_parent_vertex_global_ids(total_num_local_parent_verts), &
             local_parent_vertex_child_offsets(total_num_local_parent_verts), &
             local_parent_edge_global_ids(total_num_local_parent_edges), &
             local_parent_edge_middle_vertex_offsets(total_num_local_parent_edges))
    local_parent_vertex_global_ids = -1
    local_parent_vertex_child_offsets = -1
    local_parent_edge_global_ids = -1
    local_parent_edge_middle_vertex_offsets = -1
    global_id_offset_v = p_patch(1)%n_patch_verts_g
    global_id_offset_e = p_patch(1)%n_patch_edges_g
    ALLOCATE(parent_global_id_offset_v(n_dom), &
             parent_global_id_offset_e(n_dom))
    parent_global_id_offset_v(1) = 0
    parent_global_id_offset_e(1) = 0
    i_e = 0
    i_v = 0
    total_num_verts = p_patch(1)%n_patch_verts
    ALLOCATE(vertex_global_ids(4))

    ! for all domains (except for the first one, because it has
    ! no local parent patch)
    DO dom_id = 2, n_dom

      parent_dom_id = p_patch(dom_id)%parent_id

      ! for all edges in the local parent patch
      DO j_e = 1, p_patch_local_parent(dom_id)%n_patch_edges

        local_parent_edge_idx = idx_no(j_e)
        local_parent_edge_blk = blk_no(j_e)

        ! two edge indices matching current local parent edge
        edge_idx(1:2) = &
          p_patch_local_parent(dom_id)%edges%child_idx( &
            local_parent_edge_idx, local_parent_edge_blk, 1:2)
        edge_blk(1:2) = &
          p_patch_local_parent(dom_id)%edges%child_blk( &
            local_parent_edge_idx, local_parent_edge_blk, 1:2)

        ! outer edges have no valid child information
        IF (edge_idx(1) <= 0) CYCLE

        ! edge vertex indices
        ! (four vertices, with the middle one being duplicated)
        vertex_idx(1:2) = &
          p_patch(dom_id)%edges%vertex_idx(edge_idx(1), edge_blk(1), 1:2)
        vertex_blk(1:2) = &
          p_patch(dom_id)%edges%vertex_blk(edge_idx(1), edge_blk(1), 1:2)
        vertex_idx(3:4) = &
          p_patch(dom_id)%edges%vertex_idx(edge_idx(2), edge_blk(2), 1:2)
        vertex_blk(3:4) = &
          p_patch(dom_id)%edges%vertex_blk(edge_idx(2), edge_blk(2), 1:2)
        vertex_1d(1:4) = idx_1d(vertex_idx(1:4), vertex_blk(1:4))

        ! get global ids of edge vertices
        vertex_global_ids(1:4) = &
          p_patch(dom_id)%verts%decomp_info%glb_index( &
            vertex_1d(1:4)) + global_id_offset_v

        ! get global ids of local parent edge vertices
        local_parent_edge_vertex_idx(1:2) = &
          p_patch_local_parent(dom_id)%edges%vertex_idx( &
            local_parent_edge_idx, local_parent_edge_blk, 1:2)
        local_parent_edge_vertex_blk(1:2) = &
          p_patch_local_parent(dom_id)%edges%vertex_blk(&
            local_parent_edge_idx, local_parent_edge_blk, 1:2)
        local_parent_edge_vertex_1d(1:2) = &
          idx_1d( &
            local_parent_edge_vertex_idx(1:2), &
            local_parent_edge_vertex_blk(1:2))
        local_parent_edge_vertex_global_ids(1:2) = &
          p_patch_local_parent(dom_id)%verts%decomp_info%glb_index( &
            local_parent_edge_vertex_1d(1:2)) + &
          parent_global_id_offset_v(parent_dom_id)
        local_parent_edge_vertex_coord(1)%x(:) = &
          p_patch_local_parent(dom_id)%verts%cartesian( &
            local_parent_edge_vertex_idx(1), &
            local_parent_edge_vertex_blk(1))%x(:)
        local_parent_edge_vertex_coord(2)%x(:) = &
          p_patch_local_parent(dom_id)%verts%cartesian( &
            local_parent_edge_vertex_idx(2), &
            local_parent_edge_vertex_blk(2))%x(:)

        ! the dublicated global vertex id is the edge middle vertex
        IF (vertex_global_ids(1) == vertex_global_ids(2)) THEN
          vertex_order(1:3) = (/1,3,4/)
        ELSE IF (vertex_global_ids(1) == vertex_global_ids(3)) THEN
          vertex_order(1:3) = (/1,2,4/)
        ELSE IF (vertex_global_ids(1) == vertex_global_ids(4)) THEN
          vertex_order(1:3) = (/1,2,3/)
        ELSE IF (vertex_global_ids(2) == vertex_global_ids(3)) THEN
          vertex_order(1:3) = (/2,1,4/)
        ELSE IF (vertex_global_ids(2) == vertex_global_ids(4)) THEN
          vertex_order(1:3) = (/2,1,3/)
        ELSE IF (vertex_global_ids(3) == vertex_global_ids(4)) THEN
          vertex_order(1:3) = (/3,1,2/)
         ELSE
          CALL finish( &
            method_name, 'Unable to determine global id of edge middle vertex')
        END IF

        ! get coordinates of the first outer edge vertex
        vertex_coord%x(:) = &
          p_patch(dom_id)%verts%cartesian( &
            vertex_idx(vertex_order(2)), vertex_blk(vertex_order(2)))%x(:)

        ! determine whether the first local parent edge vertex matches
        ! the first outer edge vertex (swap otherwise)
        IF (SUM(ABS(vertex_coord%x(:) - &
                    local_parent_edge_vertex_coord(1)%x(:))) > &
            SUM(ABS(vertex_coord%x(:) - &
                    local_parent_edge_vertex_coord(2)%x(:)))) THEN
          CALL swap(vertex_order(2), vertex_order(3))
        END IF

        ! get global id of current local parent edge
        local_parent_edge_global_ids(i_e + j_e) = &
          p_patch_local_parent(dom_id)%edges%decomp_info%glb_index(j_e) + &
          parent_global_id_offset_e(parent_dom_id)

        ! get offset of edge middle vertex
        local_parent_edge_middle_vertex_offsets(i_e + j_e) = &
          total_num_verts + vertex_1d(vertex_order(1)) - 1

        ! get global id of current local parent edge vertices
        local_parent_vertex_global_ids( &
          i_v + local_parent_edge_vertex_1d(1)) = &
            local_parent_edge_vertex_global_ids(1)
        local_parent_vertex_global_ids( &
          i_v + local_parent_edge_vertex_1d(2)) = &
            local_parent_edge_vertex_global_ids(2)

        ! get global id of matching vertices
        local_parent_vertex_child_offsets( &
          i_v + local_parent_edge_vertex_1d(1)) = &
          total_num_verts + vertex_1d(vertex_order(2)) - 1
        local_parent_vertex_child_offsets( &
          i_v + local_parent_edge_vertex_1d(2)) = &
          total_num_verts + vertex_1d(vertex_order(3)) - 1

      END DO

      parent_global_id_offset_v(dom_id) = global_id_offset_v
      parent_global_id_offset_e(dom_id) = global_id_offset_e
      global_id_offset_v = &
        global_id_offset_v + p_patch(dom_id)%n_patch_verts_g
      global_id_offset_e = &
        global_id_offset_e + p_patch(dom_id)%n_patch_edges_g
      i_v = i_v + p_patch_local_parent(dom_id)%n_patch_verts
      i_e = i_e + p_patch_local_parent(dom_id)%n_patch_edges
      total_num_verts = total_num_verts + p_patch(dom_id)%n_patch_verts
    END DO

    DEALLOCATE(vertex_global_ids)
    DEALLOCATE(parent_global_id_offset_v, &
               parent_global_id_offset_e)

    !--------------------------------------------------------------------------
    ! Generate basic information required for the
    ! definition of the combined grid in YAC
    !
    ! The grid will contain all cells from all patches (some of them will be
    ! masked out, if the are covered by a nested patch (on any process))
    !
    ! nbr_vertices_per_cell: Number of vertices for each cell (cell that border
    !                        to a nest, will have additional vertices)
    !                        (masked out cell have zero vertices)
    ! cell_to_vertex: Offsets for the vertices of all cells (there are a total
    !                 of three times the number of cells plus the number of all
    !                 nest border edges)
    ! vertes_lon/lat: Coordinates for all grid vertices (vertices at the middle
    !                 of nest border edges are added to the end of this array)
    ! vertex_global_ids: global ids of all vertices (vertices at the nest
    !                    border edges will be replaces with global ids of the
    !                    respective vertices in the nested patch at a later
    !                    stage)
    !
    ! In the example:
    !  nbr_vertices_per_cell(1:32) =
    !    (/3,3,3,4,3,4,3, 3,4,0,0,0, 3,4,0, 3,
    !      3,3,3,3,3,3,3, 3,3,3,3,3, 3,3,3, 3/)
    !  cell_to_vertex(1:100) =
    !     (/0,1,5, 1,6,5, 1,2,6, 2,7,30,6, 2,3,7, 3,8,31,7, 3,4,8,
    !       5,6,9, 6,32,10,9,
    !       9,10,12, 10,32,13,12,
    !       12,13,15,
    !       14,16,20, 16,21,20, 16,17,21, 17,22,21, 17,18,22, 18,23,22, 18,19,23,
    !       .../)
    !  vertes_lon/lat(1:34) = ...
    !  vertex_global_ids(1:34) =
    !    (/1,2,3,4,5, 6,7,8,9, 10,11,12, 13,14, 15,
    !      16,17,18,19,20, 21,22,23,24, 25,26,27, 28,29, 30,
    !      *, *, *, */)
    !   *: unset values
    !--------------------------------------------------------------------------

    ALLOCATE(nbr_vertices_per_cell(total_num_cells))
    ALLOCATE(cell_to_vertex(3 * total_num_cells + num_nest_edges))
    ALLOCATE(vertex_lon(total_num_verts + num_nest_edges), &
             vertex_lat(total_num_verts + num_nest_edges))
    ALLOCATE(vertex_global_ids(total_num_verts + num_nest_edges))

    i_v = 0
    i_c = 0
    cell_to_vertex_offset = 0
    nest_vert_offset = total_num_verts
    global_id_offset_v = 0
    verts_offset = 0

    ! for all domains
    DO dom_id = 1, n_dom

      ! for all vertices of the current domain
      DO j_v = 1, p_patch(dom_id)%n_patch_verts

        jc_v = idx_no(j_v)
        jb_v = blk_no(j_v)

        ! get coordinates and global ids of current vertex
        i_v = i_v + 1
        vertex_lon(i_v) = p_patch(dom_id)%verts%vertex(jc_v,jb_v)%lon
        vertex_lat(i_v) = p_patch(dom_id)%verts%vertex(jc_v,jb_v)%lat
        vertex_global_ids(i_v) = &
          p_patch(dom_id)%verts%decomp_info%glb_index(j_v) + &
          global_id_offset_v
      END DO

      ! for all cells of the current domain
      DO j_c = 1, p_patch(dom_id)%n_patch_cells

        jc_c = idx_no(j_c)
        jb_c = blk_no(j_c)
        i_c = i_c + 1

        ! if the current cell is active
        IF (dom_cell_core_mask(dom_id)%p(jc_c, jb_c)) THEN

          curr_num_nest_edges = 0

          ! for all edges of the current cell
          DO j_e = 1, 3

            jc_e = p_patch(dom_id)%cells%edge_idx(jc_c, jb_c, j_e)
            jb_e = p_patch(dom_id)%cells%edge_blk(jc_c, jb_c, j_e)

            ! determine whether this edge has children
            IF (dom_edge_child_mask(dom_id)%p(jc_e, jb_e)) THEN

              curr_num_nest_edges = curr_num_nest_edges + 1
              nest_edge_idx(curr_num_nest_edges) = j_e
            END IF
          END DO

          nbr_vertices_per_cell(i_c) = 3

          ! if the current cell does not have a common edge
          ! with a nest cell
          IF (curr_num_nest_edges == 0) THEN

            ! add all three cell vertices
            cell_to_vertex( &
              1+cell_to_vertex_offset:3+cell_to_vertex_offset) = &
              idx_1d( &
                p_patch(dom_id)%cells%vertex_idx(jc_c, jb_c, 1:3), &
                p_patch(dom_id)%cells%vertex_blk(jc_c, jb_c, 1:3)) + &
              verts_offset
            cell_to_vertex_offset = cell_to_vertex_offset + 3

          ! if the current cell has a common edge with a nest cell
          ELSE

            ! get cell vertex indices
            vertex_idx(1:3) = &
              p_patch(dom_id)%cells%vertex_idx(jc_c, jb_c, 1:3)
            vertex_blk(1:3) = &
              p_patch(dom_id)%cells%vertex_blk(jc_c, jb_c, 1:3)

            ! copy the first cell vertex index to the last position
            ! for convenience
            vertex_idx(4) = vertex_idx(1)
            vertex_blk(4) = vertex_blk(1)

            ! for all three vertices of the current cell
            DO j_v = 1, 3

              start_vertex_idx = vertex_idx(j_v)
              start_vertex_blk = vertex_blk(j_v)
              end_vertex_idx = vertex_idx(j_v + 1)
              end_vertex_blk = vertex_blk(j_v + 1)

              ! for all edges of the current cell
              DO j_e = 1, 3

                jc_e = p_patch(dom_id)%cells%edge_idx(jc_c, jb_c, j_e)
                jb_e = p_patch(dom_id)%cells%edge_blk(jc_c, jb_c, j_e)

                ! if the current edge has the start and end vertex
                IF (((p_patch(dom_id)%edges%vertex_idx(jc_e, jb_e, 1) == &
                      start_vertex_idx) .AND. &
                     (p_patch(dom_id)%edges%vertex_blk(jc_e, jb_e, 1) == &
                      start_vertex_blk) .AND. &
                     (p_patch(dom_id)%edges%vertex_idx(jc_e, jb_e, 2) == &
                      end_vertex_idx) .AND. &
                     (p_patch(dom_id)%edges%vertex_blk(jc_e, jb_e, 2) == &
                      end_vertex_blk)) .OR. &
                    ((p_patch(dom_id)%edges%vertex_idx(jc_e, jb_e, 2) == &
                      start_vertex_idx) .AND. &
                     (p_patch(dom_id)%edges%vertex_blk(jc_e, jb_e, 2) == &
                      start_vertex_blk) .AND. &
                     (p_patch(dom_id)%edges%vertex_idx(jc_e, jb_e, 1) == &
                      end_vertex_idx) .AND. &
                     (p_patch(dom_id)%edges%vertex_blk(jc_e, jb_e, 1) == &
                      end_vertex_blk))) THEN
                  EXIT
                END IF
              END DO

              IF (j_e == 4) THEN
                CALL finish(method_name, 'Inconsistent grid data')
              END IF

              ! add the current start vertex
              cell_to_vertex_offset = cell_to_vertex_offset + 1
              cell_to_vertex(cell_to_vertex_offset) = &
                idx_1d(start_vertex_idx, start_vertex_blk) + &
                verts_offset

              ! if the current edge has a common border with a nest
              IF (dom_edge_child_mask(dom_id)%p(jc_e, jb_e)) THEN

                ! add edge middle vertex
                ! (since additional edge middle vertices for each cell are
                !  stored in the order of the cell edges, we have to use
                !  nest_edge_idx to ensure the correct edge middle vertex
                !  is referenced)
                nbr_vertices_per_cell(i_c) = nbr_vertices_per_cell(i_c) + 1
                cell_to_vertex_offset = cell_to_vertex_offset + 1
                cell_to_vertex(cell_to_vertex_offset) = &
                  nest_vert_offset + &
                  FINDLOC(nest_edge_idx(1:curr_num_nest_edges), j_e, 1)

              END IF

            END DO

            ! set offset to the next set of edge middle vertices
            nest_vert_offset = nest_vert_offset + curr_num_nest_edges

          END IF

        ! if the current cell is not active
        ELSE

          nbr_vertices_per_cell(i_c) = 0

        END IF
      END DO

      global_id_offset_v = &
        global_id_offset_v + p_patch(dom_id)%n_patch_verts_g
      verts_offset = verts_offset + p_patch(dom_id)%n_patch_verts
    END DO

    DO dom_id = 1, n_dom
      DEALLOCATE(dom_edge_child_mask(dom_id)%p)
    END DO
    DEALLOCATE(dom_edge_child_mask)

    !---------------------------------------------------------------------------
    ! Get vertex data (coordinates and global ids) for all vertices that are
    ! located on nest border edges
    !
    ! Above a lists with all required vertices and their respective location in
    ! vertex_global_ids/vertex_lon/vertex_lat have been created. In addition,
    ! lists for all potential vertices in the nested patches have been created
    ! as well. This makes it a perfect task for YAXT, where you provide a list
    ! data you have and the data you want. Using this information YAXT can
    ! generate the matching communication matrix. This can be used to generate
    ! exchange patterns for the actual data, which can be used to exchanges at
    ! once in a single call.
    !
    ! In the example (after the exchange):
    !  vertes_lon/lat(1:34) = ...
    !  vertex_global_ids(1:34) =
    !    (/1,2,3,4,5, 6,16,18,20, 10,25,27, 13,30, 15,
    !      16,17,18,19,20, 21,22,23,24, 25,26,27, 28,29, 30,
    !      17, 19, 21, 28/)
    !---------------------------------------------------------------------------

    ! build exchange map for all nest edges
    tgt_idxlist = xt_idxvec_new(required_edges)
    src_idxlist = xt_idxvec_new(local_parent_edge_global_ids)
    xmap_edge_vertex_data = &
      xt_xmap_dist_dir_new(src_idxlist, tgt_idxlist, p_comm_work)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(tgt_idxlist)
    DEALLOCATE(required_edges)
    DEALLOCATE(local_parent_edge_global_ids)

    ! build exchange map for all nest vertices
    tgt_idxlist = xt_idxvec_new(required_vertices)
    src_idxlist = xt_idxvec_new(local_parent_vertex_global_ids)
    xmap_vertex_data = &
      xt_xmap_dist_dir_new(src_idxlist, tgt_idxlist, p_comm_work)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(tgt_idxlist)
    DEALLOCATE(required_vertices)
    DEALLOCATE(local_parent_vertex_global_ids)

    ! build redistribution for nest edge and vertex data
    redist_edge_vertex_data_int = &
      xt_redist_p2p_off_new( &
        xmap_edge_vertex_data, &
        local_parent_edge_middle_vertex_offsets, &
        required_edge_middle_vertex_offsets, MPI_INTEGER)
    redist_edge_vertex_data_dble = &
      xt_redist_p2p_off_new( &
        xmap_edge_vertex_data, &
        local_parent_edge_middle_vertex_offsets, &
        required_edge_middle_vertex_offsets, MPI_DOUBLE_PRECISION)
    redist_vertex_data_int = &
      xt_redist_p2p_off_new( &
        xmap_vertex_data, &
        local_parent_vertex_child_offsets, &
        required_vertices_offsets, MPI_INTEGER)
    redist_vertex_data_dble = &
      xt_redist_p2p_off_new( &
        xmap_vertex_data, &
        local_parent_vertex_child_offsets, &
        required_vertices_offsets, MPI_DOUBLE_PRECISION)
    CALL xt_xmap_delete(xmap_vertex_data)
    CALL xt_xmap_delete(xmap_edge_vertex_data)
    DEALLOCATE(required_vertices_offsets)
    DEALLOCATE(required_edge_middle_vertex_offsets)
    DEALLOCATE(local_parent_edge_middle_vertex_offsets)
    DEALLOCATE(local_parent_vertex_child_offsets)
    redist_vertex_data = &
      xt_redist_collection_new( &
        (/redist_edge_vertex_data_int, &    ! edge middle vertex global ids
          redist_edge_vertex_data_dble, &   ! edge middle vertex lon
          redist_edge_vertex_data_dble, &   ! edge middle vertex lat
          redist_vertex_data_int, &         ! vertex global ids
          redist_vertex_data_dble, &        ! vertex lon
          redist_vertex_data_dble/), &      ! vertex lat
          p_comm_work)
    CALL xt_redist_delete(redist_edge_vertex_data_int)
    CALL xt_redist_delete(redist_edge_vertex_data_dble)
    CALL xt_redist_delete(redist_vertex_data_int)
    CALL xt_redist_delete(redist_vertex_data_dble)

    ! exchange vertex information
    IF (total_num_verts > 0) THEN
      data_ptr(1) = C_LOC(vertex_global_ids(1))
      data_ptr(2) = C_LOC(vertex_lon(1))
      data_ptr(3) = C_LOC(vertex_lat(1))
      data_ptr(4) = C_LOC(vertex_global_ids(1))
      data_ptr(5) = C_LOC(vertex_lon(1))
      data_ptr(6) = C_LOC(vertex_lat(1))
    ELSE
      data_ptr(:) = C_LOC(dummy)
    END IF
    CALL xt_redist_s_exchange( &
      redist_vertex_data, data_ptr, data_ptr)
    CALL xt_redist_delete(redist_vertex_data)

    !--------------------------------------------------------------------------
    ! The combined grid data is now provided to YAC
    !
    ! In the example the combined grid looks like this:
    !
    !                             15
    !                             / \
    !                            /   \
    !                           /     \
    !                          /       \
    !                         /         \
    !                       13----------30
    !                       / \         / \
    !                      /   \       /   \
    !                     /     \    28----29
    !                    /       \   / \   / \
    !                   /         \ /   \ /   \
    !                 10----------25----26----27
    !                 / \         / \   / \   / \
    !                /   \       /   \ /   \ /   \
    !               /     \    21----22----23----24
    !              /       \   / \   / \   / \   / \
    !             /         \ /   \ /   \ /   \ /   \
    !            6----------16----17----18----19----20
    !           / \         / \         / \         / \
    !          /   \       /   \       /   \       /   \
    !         /     \     /     \     /     \     /     \
    !        /       \   /       \   /       \   /       \
    !       /         \ /         \ /         \ /         \
    !      1-----------2-----------3-----------4-----------5
    !
    !--------------------------------------------------------------------------

    CALL yac_fdef_grid( &
      TRIM(grid_name), total_num_verts + num_nest_edges, &
      total_num_cells, &
      SUM(nbr_vertices_per_cell(1:total_num_cells)), &
      nbr_vertices_per_cell(1:total_num_cells), &
      vertex_lon, vertex_lat, cell_to_vertex, grid_id)
    DO dom_id = 1, n_dom
      DEALLOCATE(dom_cell_core_mask(dom_id)%p)
    END DO
    DEALLOCATE(nbr_vertices_per_cell, &
               vertex_lon, vertex_lat, &
               cell_to_vertex)
    DEALLOCATE(dom_cell_core_mask)

    !----------------------------------------------
    ! define cell center points
    ! (needed e.g. for patch recovery and nearest
    !  neighbour interpolation)
    !----------------------------------------------

    ALLOCATE(cell_lon(total_num_cells), cell_lat(total_num_cells))
    cells_offset = 0
    DO dom_id = 1, n_dom
      DO i = 1, p_patch(dom_id)%n_patch_cells
        jc_c = idx_no(i)
        jb_c = blk_no(i)
        cell_lon(i + cells_offset) = &
          p_patch(dom_id)%cells%center(jc_c, jb_c)%lon
        cell_lat(i + cells_offset) = &
          p_patch(dom_id)%cells%center(jc_c, jb_c)%lat
      END DO
      cells_offset = &
        cells_offset + p_patch(dom_id)%n_patch_cells
    END DO
    CALL yac_fdef_points ( &
      grid_id, total_num_cells, YAC_LOCATION_CELL, &
      cell_lon, cell_lat, cell_point_id )

    DEALLOCATE (cell_lon, cell_lat)

    !----------------------------------------------
    ! define global vertex ids
    !----------------------------------------------

    CALL yac_fset_global_index ( &
      vertex_global_ids(1:(total_num_verts+num_nest_edges)), &
      YAC_LOCATION_CORNER, grid_id)
    DEALLOCATE(vertex_global_ids)

    !----------------------------------------------
    ! define global cell ids
    !----------------------------------------------

    ALLOCATE(cell_global_ids(total_num_cells))
    cells_offset = 0
    global_id_offset_c = 0
    DO dom_id = 1, n_dom
      cell_global_ids( &
        1+cells_offset:p_patch(dom_id)%n_patch_cells+cells_offset) = &
        p_patch(dom_id)%cells%decomp_info%glb_index( &
          1:p_patch(dom_id)%n_patch_cells) + global_id_offset_c
      cells_offset = &
        cells_offset + p_patch(dom_id)%n_patch_cells
      global_id_offset_c = &
        global_id_offset_c + p_patch(dom_id)%n_patch_cells_g
    END DO

    CALL yac_fset_global_index ( &
      cell_global_ids(1:total_num_cells), YAC_LOCATION_CELL, grid_id)
    DEALLOCATE(cell_global_ids)

    !----------------------------------------------
    ! define cell core mask
    !----------------------------------------------

    CALL yac_fset_core_mask( &
      is_valid_cell(1:total_num_cells), YAC_LOCATION_CELL, grid_id)
    DEALLOCATE (is_valid_cell)

! YAC_coupling
#endif

  END SUBROUTINE def_patch_combined
  !-----------------------------------------------------------------------------

#else

  !-----------------------------------------------------------------------------
  ! Does not HAVE_YAXT, call finish
  SUBROUTINE def_patch_combined(p_patch, grid_name, grid_id, cell_point_id)

    TYPE(t_patch), INTENT(IN) :: p_patch(1:n_dom) ! array of nested patches
    CHARACTER(LEN=*), INTENT(IN) :: grid_name     ! name of the nested grid
    INTEGER, INTENT(OUT) :: grid_id               ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id         ! cell coordinate id

    CALL finish("def_patch_combined", "requires the yaxt library")

  END SUBROUTINE def_patch_combined
#endif

  ! registers the main component, grid and points
  SUBROUTINE def_main(                       &
    caller, p_patch, grid_name, with_output, &
    comp_id, output_comp_id, grid_id,        &
    cell_point_id, vertex_point_id)

    CHARACTER(LEN=*), INTENT(IN) :: caller             ! name of the calling routine (for debugging)
    TYPE(t_patch), INTENT(IN) :: p_patch(1:n_dom)      ! basic patch
    CHARACTER(LEN=*), INTENT(IN) :: grid_name          ! name of the grid
    LOGICAL, INTENT(IN) :: with_output                 ! should the output component be registered as well
    INTEGER, INTENT(OUT) :: comp_id                    ! component id
    INTEGER, INTENT(OUT) :: output_comp_id             ! component id of the output
    INTEGER, INTENT(OUT) :: grid_id(0:n_dom)           ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id(0:n_dom)     ! cell coordinate id
    INTEGER, INTENT(OUT) :: vertex_point_id(0:n_dom)   ! vertex coordinate id (only with output)

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: startdatestring
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: stopdatestring
    CHARACTER(LEN=MAX_CHAR_LENGTH)      :: comp_names(2)

    INTEGER :: comp_ids(2), jg
    INTEGER :: comp_comm, comp_rank, ierror

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':def_main', &
      'built without coupling support.')
#else

    ! Print the YAC version
    CALL message( &
      'Running ICON ' // TRIM(get_my_process_name()) // &
      ' in coupled mode with YAC version ', TRIM(yac_fget_version()) )

    ! Overwrite job start and end date with component data
    CALL datetimeToString(time_config%tc_startdate, startdatestring)
    CALL datetimeToString(time_config%tc_stopdate, stopdatestring)
    CALL yac_fdef_datetime (                  &
      yac_instance_id,                        &
      start_datetime = TRIM(startdatestring), & !in
      end_datetime   = TRIM(stopdatestring)   ) !in

    ! Inform the coupler about what we are
    IF (ltimer) CALL timer_start(timer_coupling_init_def_comp)
    IF (with_output) THEN
      comp_names(1)=TRIM(get_my_process_name())
      comp_names(2)=TRIM(get_my_process_name())//"_output"
      CALL yac_fdef_comps (                       &
        yac_instance_id,                          &
         comp_names,                              & !in
         2,                                       & !in
         comp_ids )                                 !out
      comp_id = comp_ids(1)
      output_comp_id = comp_ids(2)
    ELSE
      CALL yac_fdef_comp (           &
        yac_instance_id,             &
        TRIM(get_my_process_name()), & !in
        comp_id )                      !out
      output_comp_id = -1
    END IF
    IF (ltimer) CALL timer_stop(timer_coupling_init_def_comp)

    ! root process of the component reads in the configuration file
    CALL yac_fget_comp_comm(comp_id, comp_comm)
    CALL MPI_COMM_RANK(comp_comm, comp_rank, ierror)
    IF (comp_rank == 0 .AND. cpl_config_file_exists()) &
      CALL yac_fread_config_yaml( yac_instance_id, TRIM(yaml_filename) )
    CALL MPI_Comm_free(comp_comm, ierror)

    grid_id = -1
    cell_point_id = -1
    vertex_point_id = -1

    IF (n_dom > 1) THEN
      CALL def_patch_combined( &
        p_patch, grid_name, grid_id(0), cell_point_id(0))
      DO jg = 1, n_dom
        IF (with_output) THEN
          CALL def_patch( &
            p_patch(jg), TRIM(grid_name) // "_" // TRIM(int2str(jg)), &
            grid_id(jg), cell_point_id(jg), vertex_point_id(jg))
        ELSE
          CALL def_patch( &
            p_patch(jg), TRIM(grid_name) // "_" // TRIM(int2str(jg)), &
            grid_id(jg), cell_point_id(jg))
        END IF
      END DO
    ELSE
      jg = 1
      IF (with_output) THEN
        CALL def_patch( &
          p_patch(jg), grid_name, grid_id(jg), cell_point_id(jg), &
          vertex_point_id(jg))
      ELSE
        CALL def_patch(p_patch(jg), grid_name, grid_id(jg), cell_point_id(jg))
      END IF
    END IF

  CONTAINS

    CHARACTER(LEN=16) FUNCTION int2str(i)
        INTEGER, INTENT(IN) :: i
        WRITE (int2str, *) i
        int2str = ADJUSTL(int2str)
    END FUNCTION int2str

! YAC_coupling
#endif

  END SUBROUTINE def_main

  SUBROUTINE cpl_def_main_without_output( &
    caller, p_patch, grid_name,           &
    comp_id, grid_id, cell_point_id)

    TYPE(t_patch), INTENT(IN) :: p_patch(1:n_dom)  ! basic patch
    CHARACTER(LEN=*), INTENT(IN) :: caller         ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: grid_name      ! name of the grid
    INTEGER, INTENT(OUT) :: comp_id                ! component id
    INTEGER, INTENT(OUT) :: grid_id(0:n_dom)       ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id(0:n_dom) ! cell coordinate id

    INTEGER :: dummy_output_comp_id, dummy_vertex_point_id(0:n_dom)

    CALL def_main(                                        &
      caller // ':cpl_def_main_without_output', p_patch,  &
      grid_name, .FALSE., comp_id, dummy_output_comp_id,  &
      grid_id, cell_point_id, dummy_vertex_point_id)

  END SUBROUTINE cpl_def_main_without_output

  SUBROUTINE cpl_def_main_with_output( &
    caller, p_patch, grid_name,        &
    comp_id, output_comp_id, grid_id,  &
    cell_point_id, vertex_point_id)

    TYPE(t_patch), INTENT(IN) :: p_patch(1:n_dom)    ! basic patch
    CHARACTER(LEN=*), INTENT(IN) :: caller           ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: grid_name        ! name of the grid
    INTEGER, INTENT(OUT) :: comp_id                  ! component id
    INTEGER, INTENT(OUT) :: output_comp_id           ! component id of the output
    INTEGER, INTENT(OUT) :: grid_id(0:n_dom)         ! grid id
    INTEGER, INTENT(OUT) :: cell_point_id(0:n_dom)   ! cell coordinate id
    INTEGER, INTENT(OUT) :: vertex_point_id(0:n_dom) ! vertex coordinate id (only with output)

    CALL def_main(                                    &
      caller // ':cpl_def_main_with_output', p_patch, &
      grid_name, .TRUE., comp_id, output_comp_id,     &
      grid_id, cell_point_id, vertex_point_id)

  END SUBROUTINE cpl_def_main_with_output

  ! registers a dummy main component
  SUBROUTINE cpl_def_main_dummy(caller, comp_name)

    CHARACTER(LEN=*), INTENT(IN) :: caller    ! name of the calling routine (for debugging)
    CHARACTER(LEN=*), INTENT(IN) :: comp_name ! component name

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_def_main_dummy', 'built without coupling support.')
#else

    INTEGER :: comp_id

    IF (ltimer) CALL timer_start(timer_coupling_init_def_comp)
    CALL yac_fdef_comp( &
      yac_instance_id, TRIM(comp_name), comp_id)
    IF (ltimer) CALL timer_stop(timer_coupling_init_def_comp)
#endif

  END SUBROUTINE cpl_def_main_dummy

  ! registers a field mask for a grid
  SUBROUTINE cpl_def_cell_field_mask( &
    caller, grid_id, is_valid, mask_id)

    USE, INTRINSIC :: iso_c_binding, ONLY : c_size_t, c_int


    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN)  :: grid_id        ! grid identifier
    LOGICAL, INTENT(IN)  :: is_valid(:)    ! mask values
    INTEGER, INTENT(OUT) :: mask_id        ! mask identifier

#ifdef YAC_coupling
    CALL yac_fdef_mask (                  &
      & grid_id,                          &
      & INT(yac_fget_grid_size(           &
      &     YAC_LOCATION_CELL, grid_id)), &
      & YAC_LOCATION_CELL,                &
      & is_valid,                         &
      & mask_id )
#endif

  END SUBROUTINE cpl_def_cell_field_mask

  ! registers a field to the coupler without a mask
  SUBROUTINE cpl_def_field_no_mask( &
    comp_id, cell_point_id, timestepstring, &
    field_name, collection_size, field_id)

    INTEGER, INTENT(IN) :: comp_id                 ! component id
    INTEGER, INTENT(IN) :: cell_point_id           ! cell coordinate id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring ! time step of the field
    CHARACTER(LEN=*), INTENT(IN) :: field_name     ! name of the field
    INTEGER, INTENT(IN) :: collection_size         ! number of levels/bundle size
    INTEGER, INTENT(OUT) :: field_id               ! id of the field

#ifdef YAC_coupling
    CALL yac_fdef_field (                           &
      & field_name      = TRIM(field_name),         & !in
      & component_id    = comp_id,                  & !in
      & point_ids       = (/cell_point_id/),        & !in
      & num_pointsets   = 1,                        & !in
      & collection_size = collection_size,          & !in
      & timestep        = timestepstring,           & !in
      & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
      & field_id        = field_id )                  !out
#endif

  END SUBROUTINE cpl_def_field_no_mask

  ! registers a field to the coupler with a mask
  SUBROUTINE cpl_def_field_mask( &
    comp_id, cell_point_id, cell_mask_id, timestepstring, &
    field_name, collection_size, field_id)

    INTEGER, INTENT(IN) :: comp_id                 ! component id
    INTEGER, INTENT(IN) :: cell_point_id           ! cell coordinate id
    INTEGER, INTENT(IN) :: cell_mask_id            ! cell mask id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring ! time step of the field
    CHARACTER(LEN=*), INTENT(IN) :: field_name     ! name of the field
    INTEGER, INTENT(IN) :: collection_size         ! number of levels/bundle size
    INTEGER, INTENT(OUT) :: field_id               ! id of the field

#ifdef YAC_coupling
    CALL yac_fdef_field_mask (                      &
      & field_name      = TRIM(field_name),         & !in
      & component_id    = comp_id,                  & !in
      & point_ids       = (/cell_point_id/),        & !in
      & mask_ids        = (/cell_mask_id/),         & !in
      & num_pointsets   = 1,                        & !in
      & collection_size = collection_size,          & !in
      & timestep        = timestepstring,           & !in
      & time_unit       = YAC_TIME_UNIT_ISO_FORMAT, & !in
      & field_id        = field_id )                  !out
#endif

  END SUBROUTINE cpl_def_field_mask

  ! gets the collection size of a field
  ! (only works after the respective field has been defined and
  !  its information has been distributed among all processes either
  !  by a call to yac_fsync_def or yac_fenddef)
  FUNCTION get_field_collection_size_from_name( &
    caller, comp_name, grid_name, field_name)

    CHARACTER(LEN=*), INTENT(IN) :: comp_name  ! name of the component
    CHARACTER(LEN=*), INTENT(IN) :: grid_name  ! name of the grid
    CHARACTER(LEN=*), INTENT(IN) :: field_name ! name of the field
    CHARACTER(LEN=*), INTENT(IN) :: caller     ! name of the calling routine (for debugging)

    INTEGER :: get_field_collection_size_from_name

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':get_field_collection_size_from_name', &
      'built without coupling support.')
#else

    get_field_collection_size_from_name = &
      yac_fget_field_collection_size(yac_instance_id, &
        comp_name, grid_name, field_name)

! YAC_coupling
#endif

  END FUNCTION get_field_collection_size_from_name

  ! gets the collection size of a field
  FUNCTION get_field_collection_size_from_id(caller, field_id)

    INTEGER, INTENT(IN) :: field_id  ! id of the field
    CHARACTER(LEN=*), INTENT(IN) :: caller     ! name of the calling routine (for debugging)

    INTEGER :: get_field_collection_size_from_id

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':get_field_collection_size_from_id', &
      'built without coupling support.')
#else

    get_field_collection_size_from_id = &
      yac_fget_collection_size_from_field_id(field_id)

! YAC_coupling
#endif

  END FUNCTION get_field_collection_size_from_id


  ! gets the current time of a coupled field
  ! (only works after the respective field and the associated couplings have
  !  been definied and its information has been distributed among all processes
  !  either by a call to yac_fsync_def or yac_fenddef)
  !
  ! it is initialised with:
  !   field_datetime = start_datetime + lag * field_timestep
  ! and incremented with every put/get/exchange/update call by:
  !   field_datetime = field_datetime + field_timestep
  FUNCTION cpl_get_field_datetime(caller, field_id)

    INTEGER, INTENT(IN) :: field_id
    CHARACTER(LEN=*), INTENT(IN) :: caller     ! name of the calling routine (for debugging)

    TYPE(datetime), POINTER :: field_datetime_ptr

    TYPE(datetime) :: cpl_get_field_datetime

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_datetime', &
      'built without coupling support.')
#else

    field_datetime_ptr => &
      newdatetime(TRIM(yac_fget_field_datetime(field_id)))

    cpl_get_field_datetime = field_datetime_ptr

    CALL deallocateDatetime(field_datetime_ptr)

! YAC_coupling
#endif

  END FUNCTION cpl_get_field_datetime

#ifdef YAC_coupling

  ! basic routine for sending a field through the coupler
  SUBROUTINE put( &
    caller, field_id, field_name, field, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller
    INTEGER, INTENT(IN) :: field_id
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    TYPE(yac_dble_ptr), INTENT(IN) :: field(:, :)
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart

    INTEGER :: num_pointsets, collection_size, put_timer
    INTEGER :: info, ierr

    num_pointsets = SIZE(field, 1)
    collection_size = SIZE(field, 2)

    CALL yac_fget_action(field_id, info)

    IF (ltimer) THEN

      SELECT CASE (info)
        CASE (YAC_ACTION_NONE,YAC_ACTION_OUT_OF_BOUND)
          put_timer = timer_coupling_nop
        CASE (YAC_ACTION_REDUCTION)
          put_timer = timer_coupling_put_reduce
        CASE DEFAULT
          put_timer = timer_coupling_put
      END SELECT

      CALL timer_start(put_timer)
    END IF

    IF ((info == YAC_ACTION_NONE) .OR. (info == YAC_ACTION_OUT_OF_BOUND)) THEN

      ! update internal clock without an actual put
      CALL yac_fupdate(field_id)

    ELSE

      CALL yac_fput( &
        field_id, num_pointsets, collection_size, field, info, ierr)

    END IF

    IF (ltimer) CALL timer_stop(put_timer)

    IF ( info == YAC_ACTION_PUT_FOR_RESTART ) THEN
      CALL message( &
        caller // ':put', &
        'YAC says it is put for restart - ' // TRIM(field_name))
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      CALL warning( &
        caller // ':put', &
        'YAC says put called after end of run - ' // TRIM(field_name))
    ENDIF

    IF (PRESENT(write_restart)) &
      write_restart = (info == YAC_ACTION_PUT_FOR_RESTART)

  END SUBROUTINE put

! YAC_coupling
#endif

  ! sends one or more fields through the coupler
  ! remark:
  !   * field data has the dimensions (nidx,nblk)
  !     with nidx * nblk >= num_points
  !   * number of provided fields has to match the collection size
  !     associated with the provided field id
  SUBROUTINE cpl_put_field_idx_blk_collection( &
    caller, field_collection_id, field_collection_name, num_points, &
    field_1, field_2, field_3, field_4, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller                            ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_collection_id                        ! field id of the field collection
    CHARACTER(LEN=*), INTENT(IN) :: field_collection_name             ! name of the field collection (for debugging)
    INTEGER, INTENT(IN) :: num_points                                 ! number of points in the field data (e.g. number of cells)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN):: field_1(:,:)           ! field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(IN):: field_2(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(IN):: field_3(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(IN):: field_4(:,:) ! optional field data
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart                   ! .TRUE. if it was the last valid put

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_put_field_idx_blk_collection', &
      'built without coupling support.')
#else

    INTEGER :: collection_size

    TYPE(yac_dble_ptr) :: send_field_collection(1, 4)

    collection_size = 1
    send_field_collection(1,collection_size)%p(1:num_points) => &
      field_1(:,:)
    IF (PRESENT(field_2)) THEN
      collection_size = collection_size + 1
      send_field_collection(1,collection_size)%p(1:num_points) => &
        field_2(:,:)
    END IF
    IF (PRESENT(field_3)) THEN
      collection_size = collection_size + 1
      send_field_collection(1,collection_size)%p(1:num_points) => &
        field_3(:,:)
    END IF
    IF (PRESENT(field_4)) THEN
      collection_size = collection_size + 1
      send_field_collection(1,collection_size)%p(1:num_points) => &
        field_4(:,:)
    END IF

    CALL put( &
      caller // ':cpl_put_field_idx_blk_collection', field_collection_id, &
      field_collection_name, send_field_collection(:,1:collection_size), &
      write_restart)

! YAC_coupling
#endif

  END SUBROUTINE cpl_put_field_idx_blk_collection

  ! sends collection of fields through the coupler
  ! remark:
  !   * field data has the dimensions (num_points,collection_size)
  !   * number of provided fields has to match the collection size
  !     associated with the provided field id
  SUBROUTINE cpl_put_field_n_collection( &
    caller, field_collection_id, field_collection_name, &
    field_collection, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller                            ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_collection_id                        ! field id of the field collection
    CHARACTER(LEN=*), INTENT(IN) :: field_collection_name             ! name of the field collection (for debugging)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(IN):: field_collection(:,:)  ! field data
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart                   ! .TRUE. if it was the last valid put

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_put_field_n_collection', &
      'built without coupling support.')
#else

    INTEGER :: num_points, collection_size
    INTEGER :: i

    TYPE(yac_dble_ptr), ALLOCATABLE :: send_field_collection(:,:)

    num_points = SIZE(field_collection, 1)
    collection_size = SIZE(field_collection, 2)

    ALLOCATE(send_field_collection(1,collection_size))

    DO i = 1, collection_size
      send_field_collection(1,i)%p(1:num_points) => field_collection(:,i)
    END DO

    CALL put( &
      caller // ':cpl_put_field_n_collection', field_collection_id, &
      field_collection_name, send_field_collection(:,1:collection_size), &
      write_restart)

! YAC_coupling
#endif

  END SUBROUTINE cpl_put_field_n_collection

#ifdef YAC_coupling

  ! basic routine for receiving a field through the coupler
  SUBROUTINE get( &
    caller, field_id, field_name, field, &
    first_get, received_data, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller
    INTEGER, INTENT(IN) :: field_id
    CHARACTER(LEN=*), INTENT(IN) :: field_name
    TYPE(yac_dble_ptr), INTENT(INOUT) :: field(:)
    LOGICAL, INTENT(IN) :: first_get
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart

    INTEGER :: collection_size
    INTEGER :: get_timer, info, ierr

    ! Skip time measurement of the very first yac_fget
    ! as this will measure mainly the wait time caused
    ! by the initialisation of the model components
    ! and does not tell us much about the load balancing
    ! in subsequent calls.
    LOGICAL, SAVE :: lyac_very_1st_get = .TRUE.

    collection_size = SIZE(field, 1)

    CALL yac_fget_action(field_id, info)

    IF (ltimer) THEN

      SELECT CASE (info)
        CASE (YAC_ACTION_NONE,YAC_ACTION_OUT_OF_BOUND)
          get_timer = timer_coupling_nop
        CASE DEFAULT
          IF (first_get) THEN
            get_timer = &
              MERGE( &
                timer_coupling_very_1stget, &
                timer_coupling_1stget, &
                lyac_very_1st_get)
          ELSE
            get_timer = timer_coupling_get
          END IF
          lyac_very_1st_get = .FALSE.
      END SELECT

      CALL timer_start(get_timer)
    END IF

    IF ((info == YAC_ACTION_NONE) .OR. (info == YAC_ACTION_OUT_OF_BOUND)) THEN

      ! update internal clock without an actual get
      CALL yac_fupdate(field_id)

    ELSE

      CALL yac_fget(field_id, collection_size, field, info, ierr)

    END IF

    IF (ltimer) CALL timer_stop(get_timer)

    IF ( info == YAC_ACTION_GET_FOR_RESTART ) THEN
      CALL message( &
        caller // ':get', &
        'YAC says it is get for restart - ' // TRIM(field_name))
    ENDIF
    IF ( info == YAC_ACTION_OUT_OF_BOUND ) THEN
      CALL warning( &
        caller // ':get', &
        'YAC says get called after end of run - ' // TRIM(field_name))
    ENDIF

    IF (PRESENT(received_data)) &
      received_data = &
        (info == YAC_ACTION_COUPLING) .OR. (info == YAC_ACTION_GET_FOR_RESTART)
    IF (PRESENT(write_restart)) &
      write_restart = (info == YAC_ACTION_GET_FOR_RESTART)

  END SUBROUTINE get

! YAC_coupling
#endif

  ! receives one or more fields through the coupler
  ! remark:
  !   * field data has the dimensions (num points, collection size)
  !   * all field data is provided in a single contiguous buffer
  !   * collection size has to match the one associated with the provided
  !     field id
  !   * depending on the field and coupling timestep, no data my actually
  !     be received by this call
  SUBROUTINE cpl_get_field_n_collection( &
    caller, field_collection_id, field_collection_name, &
    field_collection, first_get, received_data, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller                ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_collection_id            ! field id of the field collection
    CHARACTER(LEN=*), INTENT(IN) :: field_collection_name ! name of the field collection (for debugging)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT):: &
      field_collection(:,:)                               ! field data
    LOGICAL, OPTIONAL, INTENT(IN) :: first_get            ! is first get of timestep
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data       ! .TRUE. if data was received by this call
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart       ! .TRUE. if it was the last valid get

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_n_collection', &
      'built without coupling support.')
#else

    INTEGER :: num_points, collection_size
    INTEGER :: i
    LOGICAL :: lfirst_get

    TYPE(yac_dble_ptr) :: recv_field_collection(SIZE(field_collection,2))

    num_points = SIZE(field_collection, 1)
    collection_size = SIZE(field_collection, 2)

    DO i = 1, collection_size
      recv_field_collection(i)%p(1:num_points) => field_collection(:,i)
    END DO

    lfirst_get = .FALSE.
    IF (PRESENT(first_get)) lfirst_get = first_get

    CALL get( &
      caller // ':cpl_get_field_n_collection', field_collection_id, &
      field_collection_name, recv_field_collection, &
      lfirst_get, received_data, write_restart)

! YAC_coupling
#endif

  END SUBROUTINE cpl_get_field_n_collection

  ! receives one or more fields through the coupler
  ! remark:
  !   * field data has the dimensions (nidx,nblk)
  !     with nidx * nblk >= num_points
  !   * number of provided fields has to match the collection size
  !     associated with the provided field id
  !   * depending on the field and coupling timestep, no data my actually
  !     be received by this call
  SUBROUTINE cpl_get_field_idx_blk_collection( &
    caller, field_collection_id, field_collection_name, num_points, &
    field_1, field_2, field_3, field_4, first_get, received_data, &
    write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller                               ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_collection_id                           ! field id of the field collection
    CHARACTER(LEN=*), INTENT(IN) :: field_collection_name                ! name of the field collection (for debugging)
    INTEGER, INTENT(IN) :: num_points                                    ! number of points in the field data (e.g. number of cells)
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT):: field_1(:,:)           ! field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(INOUT):: field_2(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(INOUT):: field_3(:,:) ! optional field data
    REAL(wp), CONTIGUOUS, TARGET, OPTIONAL, INTENT(INOUT):: field_4(:,:) ! optional field data
    LOGICAL, OPTIONAL, INTENT(IN) :: first_get                           ! is first get of timestep
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data                      ! .TRUE. if data was received by this call
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart                      ! .TRUE. if it was the last valid get

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_idx_blk_collection', &
      'built without coupling support.')
#else

    INTEGER :: collection_size
    LOGICAL :: lfirst_get

    TYPE(yac_dble_ptr) :: recv_field_collection(4)

    collection_size = 1
    recv_field_collection(collection_size)%p(1:num_points) => &
      field_1(:,:)
    IF (PRESENT(field_2)) THEN
      collection_size = collection_size + 1
      recv_field_collection(collection_size)%p(1:num_points) => &
        field_2(:,:)
    END IF
    IF (PRESENT(field_3)) THEN
      collection_size = collection_size + 1
      recv_field_collection(collection_size)%p(1:num_points) => &
        field_3(:,:)
    END IF
    IF (PRESENT(field_4)) THEN
      collection_size = collection_size + 1
      recv_field_collection(collection_size)%p(1:num_points) => &
        field_4(:,:)
    END IF

    lfirst_get = .FALSE.
    IF (PRESENT(first_get)) lfirst_get = first_get

    CALL get( &
      caller // ':cpl_get_field_idx_blk_collection', field_collection_id, &
      field_collection_name, recv_field_collection(1:collection_size), &
      lfirst_get, received_data, write_restart)

! YAC_coupling
#endif

  END SUBROUTINE cpl_get_field_idx_blk_collection

  ! receives multiple levels of a single field through the coupler
  ! remark:
  !   * field data has the dimensions (nidx,nlev,nblk)
  !     with nidx * nblk >= num_points
  !   * receive buffer has the dimensions (num points, nlev_)
  !     with nlev_ >= nlev
  !   * number of levels match the collection size associated with
  !     the provided field id
  !   * depending on the field and coupling timestep, no data my actually
  !     be received by this call
  SUBROUTINE cpl_get_field_idx_lev_blk( &
    caller, field_id, field_name, field, recv_buf, scale_factor, &
    first_get, received_data, write_restart)

    CHARACTER(LEN=*), INTENT(IN) :: caller          ! name of the calling routine (for debugging)
    INTEGER, INTENT(IN) :: field_id                 ! field id of the field
    CHARACTER(LEN=*), INTENT(IN) :: field_name      ! name of the field (for debugging)
    REAL(wp), INTENT(INOUT) :: field(:,:,:)         ! field data
    REAL(wp), CONTIGUOUS, TARGET, INTENT(INOUT) :: &
      recv_buf(:,:)                                 ! contiguous temporary buffer used by this routine
    REAL(wp), OPTIONAL, INTENT(IN) :: scale_factor  ! optional: multiply whole field by this factor
                                                    ! (only if data was received)
    LOGICAL, OPTIONAL, INTENT(OUT) :: received_data ! .TRUE. if data was received by this call
    LOGICAL, OPTIONAL, INTENT(IN) :: first_get      ! is first get of timestep
    LOGICAL, OPTIONAL, INTENT(OUT) :: write_restart ! .TRUE. if it was the last valid get

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_get_field_idx_lev_blk', &
      'built without coupling support.')
#else

    INTEGER :: nidx, nlev, nblk, num_points
    LOGICAL :: coupling
    INTEGER :: i, j, k
    LOGICAL :: lfirst_get

    TYPE(yac_dble_ptr) :: recv_field(SIZE(field, 2))

    nidx = SIZE(field, 1)
    nlev = SIZE(field, 2)
    nblk = SIZE(field, 3)
    num_points = SIZE(recv_buf, 1)

    IF (nlev > SIZE(recv_buf, 2)) &
      CALL finish( &
        TRIM(caller) // ':cpl_get_field_idx_lev_blk', &
        'insufficient recv_buf size')

    ! MoHa:
    !   remarks:
    !     * Since independent fields may use different field mask
    !       cells written by yac_fget may be different as well. Therefore
    !       it is easier if the caller provides the receive buffer. Otherwise
    !       we would have to reinitialise it every time.

    DO i = 1, nlev
      recv_field(i)%p => recv_buf(:,i)

    lfirst_get = .FALSE.
    IF (PRESENT(first_get)) lfirst_get = first_get
    END DO

    CALL get( &
      caller // ':cpl_get_field_idx_lev_blk', field_id, field_name, recv_field, &
      lfirst_get, coupling, write_restart)

    ! if data was received
    IF (coupling) THEN

      ! unpack data
      IF (PRESENT(scale_factor)) THEN
        DO i = 1, nblk
          DO j = 1, nlev
            DO k = 1, nidx
              IF ((i-1) * nidx + k > num_points) CYCLE
              field(k,j,i) = scale_factor * recv_buf((i-1)*nidx+k,j)
            END DO
          END DO
        END DO
      ELSE
        DO i = 1, nblk
          DO j = 1, nlev
            DO k = 1, nidx
              IF ((i-1) * nidx + k > num_points) CYCLE
              field(k,j,i) = recv_buf((i-1)*nidx+k,j)
            END DO
          END DO
        END DO
      END IF

    END IF

    IF (PRESENT(received_data)) received_data = coupling

! YAC_coupling
#endif

  END SUBROUTINE cpl_get_field_idx_lev_blk

  SUBROUTINE cpl_sync_def(caller)

    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_sync_def', 'built without coupling support.')
#else
    CALL yac_fsync_def(yac_instance_id)
#endif

  END SUBROUTINE

  SUBROUTINE cpl_enddef(caller)

    CHARACTER(LEN=*), INTENT(IN) :: caller ! name of the calling routine (for debugging)

#ifndef YAC_coupling
    CALL finish( &
      TRIM(caller) // ':cpl_enddef', 'built without coupling support.')
#else
    IF (ltimer) CALL timer_start(timer_coupling_init_enddef)
    CALL yac_fenddef(yac_instance_id)
    IF (ltimer) CALL timer_stop(timer_coupling_init_enddef)
#endif

  END SUBROUTINE

END MODULE mo_coupling_utils
