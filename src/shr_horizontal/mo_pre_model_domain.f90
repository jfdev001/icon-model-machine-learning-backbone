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

MODULE mo_pre_model_domain
  USE mo_impl_constants,          ONLY: max_dom
  USE mo_io_units,                ONLY: filename_max
  USE mo_util_uuid_types,         ONLY: t_uuid
  USE ppm_extents,                ONLY: extent
  USE ppm_distributed_array,      ONLY: dist_mult_array
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pre_grid_cells
  PUBLIC :: t_pre_grid_edges
  !PUBLIC :: t_pre_grid_vertices
  PUBLIC :: t_pre_patch

  ! index of distributed sub-arrays in t_pre_grid_cells%dist
  INTEGER, PUBLIC, PARAMETER :: &
       c_num_edges = 8, &    ! number of edges connected to cell
       c_parent = 2, &    ! index of parent triangle:
       ! parent child index, number of current cell in parent's child_idx/child_blk:
       ! indices of child triangles:
       ! index2=1,4
       c_child = 3, &
       ! physical domain ID of triangles
       ! (may differ from the "normal" domain ID in case of domain merging):
       c_phys_id = 4, &
       ! indices of triangles next to each cell:
       ! index2=1,3
       c_neighbor = 5, &
       ! indices of edges of triangle:
       ! index2=1,3
       c_edge = 6, &
       ! indices of verts of triangle:
       ! index2=1,3
       c_vertex = 7, &
       ! cell geometry
       ! longitude & latitude of centers of triangular cells
       c_center = 1, &
       ! refinement control flag
       c_refin_ctrl = 9


  TYPE t_pre_grid_cells

    ! extents of the local chunk of the distributed arrays
    TYPE(extent) :: local_chunk(1,1)

    INTEGER :: max_connectivity
    TYPE(dist_mult_array) :: dist

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell (defined in mo_impl_constants)
    INTEGER, ALLOCATABLE :: start(:)

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell
    INTEGER, ALLOCATABLE :: end(:)

  END TYPE t_pre_grid_cells

    ! index of distributed sub-arrays in t_pre_grid_edges%dist
  INTEGER, PUBLIC, PARAMETER :: &
       ! index of parent edge:
       e_parent = 1, &
       ! indices of child edges:
       ! index2=1,4
       e_child = 2, &
       ! indices of adjacent cells:
       ! index2=1,2
       e_cell = 3, &
       ! refinement control flag
       e_refin_ctrl = 4

  TYPE t_pre_grid_edges

    ! extents of the local chunk of the distributed arrays
    TYPE(extent) :: local_chunk(1,1)

    TYPE(dist_mult_array) :: dist

    ! list of start indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge (defined in mo_impl_constants)
    INTEGER, ALLOCATABLE :: start(:)

    ! list of end indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge
    INTEGER, ALLOCATABLE :: end(:)

  END TYPE t_pre_grid_edges

  ! index of distributed sub-arrays in t_pre_grid_vertices%dist
  INTEGER, PUBLIC, PARAMETER :: &
       ! line indices of cells around each vertex:
       ! index2=1,6
       v_cell = 3, &
       ! number of edges connected to vertex
       v_num_edges = 2, &
       ! longitude & latitude of vertex:
       ! index2=1,2
       v_vertex = 1, &
       ! refinement control flag
       v_refin_ctrl = 4

  TYPE t_pre_grid_vertices

    ! extents of the local chunk of the distributed arrays
    TYPE(extent) :: local_chunk(1,1)

    INTEGER :: max_connectivity

    TYPE(dist_mult_array) :: dist

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert (defined in mo_impl_constants)
    INTEGER, ALLOCATABLE :: start(:)

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert
    INTEGER, ALLOCATABLE :: end(:)

  END TYPE t_pre_grid_vertices

  TYPE t_pre_patch

    !
    ! !  level in grid hierarchy on which patch lives
    !
    CHARACTER(LEN=filename_max) :: grid_filename, grid_filename_grfinfo
    !
    ! uuid of grid
    TYPE(t_uuid) :: grid_uuid
    !
    ! grid level
    INTEGER :: level
    !
    ! nest level = grid level - start level
    INTEGER :: nest_level
    !
    ! domain ID of current domain
    INTEGER :: id

    !-------------------------------------
    !> The grid domain geometry parameters
    ! cell type =3 or 6
    INTEGER :: cell_type

    !
    ! domain ID of parent domain
    INTEGER :: parent_id
    !
    ! child domain index of current domain as seen from parent domain
    ! In other words: I am the nth child of my parents (n=parent_child_index)
    INTEGER :: parent_child_index
    !
    ! list of child domain ID's
    INTEGER :: child_id(max_dom)
    !
    ! actual number of child domains
    INTEGER :: n_childdom
    !
    ! total number of child domains in the calling tree (over all nest levels)
    INTEGER :: n_chd_total
    !
    ! corresponding list of child domain ID's
    INTEGER :: child_id_list(max_dom)
    !
    ! maximum number of child domains
    INTEGER :: max_childdom
    !
    ! ! number of cells, edges and vertices in the global patch
    INTEGER :: n_patch_cells_g
    INTEGER :: n_patch_edges_g
    INTEGER :: n_patch_verts_g
    !
    ! ! values for the blocking
    !
    INTEGER :: alloc_cell_blocks  ! number of allocated cell blocks
    ! number of blocks and chunk length in last block
    ! ... for the cells
    INTEGER :: nblks_c
    INTEGER :: npromz_c
    ! ... for the edges
    INTEGER :: nblks_e
    INTEGER :: npromz_e
    ! ... for the vertices
    INTEGER :: nblks_v
    INTEGER :: npromz_v
    !
    ! ! vertical full and half levels
    !
    ! number of full and half levels
    INTEGER :: nlev
    INTEGER :: nlevp1
    !
    ! half level of parent domain (jg-1) that coincides
    ! with the upper margin of the current domain jg
    INTEGER :: nshift
    !
    ! total shift of model top with respect to global domain
    INTEGER :: nshift_total
    !
    ! the same information seen from the parent level (duplication needed to simplify flow control)
    INTEGER :: nshift_child
    !
    ! ! grid information on the patch
    !
    TYPE(t_pre_grid_cells) :: cells
    TYPE(t_pre_grid_edges) :: edges
    TYPE(t_pre_grid_vertices) :: verts
    !
    ! Most of the data of the cells, edges and verts in t_pre_patch is stored in
    ! dist_mult_array. For performance reasons the process in p_comm_work are
    ! split into one or multiple contiguous chunks. Together the processes of
    ! each chunk has a complete copy of the cells, edges and verts data. The
    ! number of chunks is set by the name list parameter num_dist_array_replicas
    !
    ! ! communicator containing all process of the process chunk
    !
    INTEGER :: dist_array_comm
    !
    ! ! process ranks (+1) of the chunk the local process is a part of
    !
    TYPE(extent) :: dist_array_pes

  END TYPE t_pre_patch


END MODULE mo_pre_model_domain
