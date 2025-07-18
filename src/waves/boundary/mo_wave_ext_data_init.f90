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

! Initialization/reading reading of external datasets
!
! This module contains read and initialization routines for the external data state.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_wave_ext_data_init

  USE mo_kind,                ONLY: wp, vp
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, SUCCESS
  USE mo_io_units,            ONLY: filename_max
  USE mo_io_config,           ONLY: default_read_method
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: n_dom, nroot
  USE mo_extpar_config,       ONLY: extpar_filename, generate_filename
  USE mo_read_interface,      ONLY: openInputFile, closeFile, t_stream_id, on_cells, read_2D
  USE mo_master_config,       ONLY: getModelBaseDir
  USE mo_intp,                ONLY: cells2edges_scalar
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_fortran_tools,       ONLY: copy, init
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_process_topo,        ONLY: compute_smooth_topo

  USE mo_wave_ext_data_types, ONLY: t_external_wave
  USE mo_wave_config,         ONLY: wave_config
  USE mo_math_gradients,      ONLY: grad_green_gauss_cell

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_data_init'

  PUBLIC :: init_wave_ext_data
  PUBLIC :: cells2edges_bathymetry
  PUBLIC :: compute_depth_gradient
  PUBLIC :: init_coastedge_list

CONTAINS

  !>
  !! read ext data state, and do some postprocessing which includes
  !! - bathymetry smoothing
  !! - computation of the bathymetry gradient
  !! - initialization of water depth
  !!
  SUBROUTINE init_wave_ext_data (p_patch, p_int_state, wave_ext_data)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_wave_ext_data'

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)
    TYPE(t_int_state),     INTENT(IN)    :: p_int_state(:)
    TYPE(t_external_wave), INTENT(INOUT) :: wave_ext_data(:)

    INTEGER :: jg
    INTEGER :: ist

    CHARACTER(len=20) :: depth_str

    REAL(wp), ALLOCATABLE :: topo_smt_c(:,:)


    DO jg = 1, n_dom
      IF (wave_config(jg)%depth > 0.0_wp) THEN
        ! set depth to constant
        write(depth_str,'(f9.2)') wave_config(jg)%depth
        CALL message('###  Run with constant depth of',TRIM(depth_str)//', m')
        wave_ext_data(jg)%bathymetry_c(:,:) = wave_config(jg)%depth
      ELSE

        ! read bathymetry on cells from file.
        ! Includes limiting to the range [depth_min,depth_max].
        !
        CALL read_ext_data_wave(p_patch(jg), wave_ext_data(jg)%bathymetry_c(:,:))


        ! ensure that the bathymetric height stays in the range [depth_min,depth_max]
        !
        CALL set_bathymetry_limits(p_patch      = p_patch(jg),                       & !in
          &                        depth_min    = wave_config(jg)%depth_min,         & !in
          &                        depth_max    = wave_config(jg)%depth_max,         & !in
          &                        bathymetry_c = wave_ext_data(jg)%bathymetry_c(:,:)) !inout

        ! smooth bathymetry
        IF (wave_config(jg)%niter_smooth > 0) THEN
          ALLOCATE(topo_smt_c(&
               SIZE(wave_ext_data(jg)%bathymetry_c,1), &
               SIZE(wave_ext_data(jg)%bathymetry_c,2)),&
               stat=ist)
          IF (ist/=SUCCESS) CALL finish(routine, &
               'allocation of topo_smt_c array failed')

          CALL compute_smooth_topo(p_patch(jg), p_int_state(jg), &
               wave_ext_data(jg)%bathymetry_c, &
               wave_config(jg)%niter_smooth, &
               topo_smt_c)

!$OMP PARALLEL
          CALL copy(src=topo_smt_c,dest=wave_ext_data(jg)%bathymetry_c, lacc=.FALSE.)
!$OMP END PARALLEL

          DEALLOCATE(topo_smt_c,stat=ist)
          IF (ist/=SUCCESS) CALL finish(routine, &
               'deallocation of topo_smt_c array failed')
        END IF
      END IF


      ! cell2edge interpolation of bathymetry
      CALL cells2edges_bathymetry(p_patch      = p_patch(jg),                    &
        &                         p_int_state  = p_int_state(jg),                &
        &                         bathymetry_c = wave_ext_data(jg)%bathymetry_c, &
        &                         bathymetry_e = wave_ext_data(jg)%bathymetry_e)


      ! init water depth
!$OMP PARALLEL
      CALL copy(src=wave_ext_data(jg)%bathymetry_c, dest=wave_ext_data(jg)%depth_c, lacc=.FALSE.)
      CALL copy(src=wave_ext_data(jg)%bathymetry_e, dest=wave_ext_data(jg)%depth_e, lacc=.FALSE.)
!$OMP END PARALLEL


      ! calculate depth gradient
      CALL compute_depth_gradient(p_patch          = p_patch(jg),                     & !in
                                  p_int_state      = p_int_state(jg),                 & !in
                                  depth_c          = wave_ext_data(jg)%depth_c,       & !in
                                  geo_depth_grad_c = wave_ext_data(jg)%geo_depth_grad_c)!out
    END DO  !jg


    CALL message(routine,'finished.')

  END SUBROUTINE init_wave_ext_data


  !>
  !! Ensure that the bathymetric height stays in the user-defined range [depth_min,depth_max]
  !!
  SUBROUTINE set_bathymetry_limits(p_patch, depth_min, depth_max, bathymetry_c)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(IN)    :: depth_min       ! lower limit for bathymetry [m]
    REAL(wp),           INTENT(IN)    :: depth_max       ! upper limit for bathymetry [m]
    REAL(wp),           INTENT(INOUT) :: bathymetry_c(:,:)

    INTEGER :: jb, jc
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx,i_endidx

    rl_start   = 1
    rl_end     = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! set depth limits depth_min and depth_max
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb=i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        bathymetry_c(jc,jb) = MAX(bathymetry_c(jc,jb),depth_min)
        bathymetry_c(jc,jb) = MIN(bathymetry_c(jc,jb),depth_max)
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE set_bathymetry_limits


  !>
  !! cell2edge interpolation for bathymetry
  !!
  SUBROUTINE cells2edges_bathymetry(p_patch, p_int_state, bathymetry_c, bathymetry_e)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    TYPE(t_int_state),  INTENT(IN)    :: p_int_state
    REAL(wp),           INTENT(IN)    :: bathymetry_c(:,:)
    REAL(wp),           INTENT(INOUT) :: bathymetry_e(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':cells2edges_bathymetry'

    REAL(wp):: bath_c_3d(SIZE(bathymetry_c,1),1,SIZE(bathymetry_c,2))
    REAL(wp):: bath_e_3d(SIZE(bathymetry_e,1),1,SIZE(bathymetry_e,2))


!$OMP PARALLEL
    CALL copy(src=bathymetry_c, dest=bath_c_3d(:,1,:), lacc=.FALSE.)
    CALL init(bath_e_3d(:,:,:), lacc=.FALSE.)
!$OMP END PARALLEL

    CALL cells2edges_scalar(p_cell_in      = bath_c_3d(:,:,:),    &
      &                     ptr_patch      = p_patch,             &
      &                     c_int          = p_int_state%c_lin_e, &
      &                     p_edge_out     = bath_e_3d(:,:,:),    &
      &                     lacc           = .FALSE.,             &
      &                     opt_has_latbcs = .TRUE.,              & ! set values at boundary edges
      &                     opt_fill_latbc = .TRUE.               )

!$OMP PARALLEL
    CALL copy(src=bath_e_3d(:,1,:), dest=bathymetry_e, lacc=.FALSE.)
!$OMP END PARALLEL

  END SUBROUTINE cells2edges_bathymetry


  !>
  !! calculate bathymetry gradient
  !!
  !!
  SUBROUTINE compute_depth_gradient(p_patch, p_int_state, depth_c, geo_depth_grad_c)

    TYPE(t_patch),     INTENT(IN)    :: p_patch
    TYPE(t_int_state), INTENT(IN)    :: p_int_state
    REAL(wp),          INTENT(IN)    :: depth_c(:,:)
    REAL(wp),          INTENT(INOUT) :: geo_depth_grad_c(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_depth_gradient'

    REAL(wp) :: depth_c_3d(SIZE(depth_c,1),1,SIZE(depth_c,2))
    REAL(vp) :: depth_grad_c_4d(SIZE(geo_depth_grad_c,1),SIZE(geo_depth_grad_c,2),1,SIZE(geo_depth_grad_c,3))

!$OMP PARALLEL
    CALL copy(src=depth_c, dest=depth_c_3d(:,1,:), lacc=.FALSE.)
    CALL init(depth_grad_c_4d(:,:,:,:), lacc=.FALSE.)
!$OMP END PARALLEL

    CALL grad_green_gauss_cell(depth_c_3d, p_patch, p_int_state, depth_grad_c_4d, &
         &                     lacc=.FALSE., opt_slev=1, opt_elev=1, &
         &                     opt_rlstart=2, opt_rlend=min_rlcell)

!$OMP PARALLEL
    CALL copy(src=depth_grad_c_4d(:,:,1,:), dest=geo_depth_grad_c(:,:,:), lacc=.FALSE.)
!$OMP END PARALLEL

  END SUBROUTINE compute_depth_gradient


  !>
  !! read bathymetric height on cell centers from external parameter file
  !! and ensure that the field falls within the range [depth_min,depth_max]
  !!
  SUBROUTINE read_ext_data_wave(p_patch, bathymetry_c)

    TYPE(t_patch),  INTENT(IN)    :: p_patch
    REAL(wp),       INTENT(INOUT) :: bathymetry_c(:,:)  ! bathymetric height on cells

    TYPE(t_stream_id) :: stream_id

    CHARACTER(filename_max) :: extpar_file

    extpar_file = generate_filename(extpar_filename, getModelBaseDir(), &
      &                             TRIM(p_patch%grid_filename),        &
      &                             nroot,                              &
      &                             p_patch%level, p_patch%id)

    CALL openInputFile(stream_id, extpar_file, p_patch, default_read_method)

    CALL read_2D(stream_id, on_cells, 'z', bathymetry_c)

    CALL closeFile(stream_id)

  END SUBROUTINE read_ext_data_wave


  !>
  !! computes index list for coastal edge points (edges%refin_ctrl = 1) and stores
  !! corresponding edge orientations
  !!
  SUBROUTINE init_coastedge_list (p_patch)

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)

    ! local variables
    INTEGER :: jg, jb, jc, jce
    INTEGER :: eidx, eblk            !< edge index and block
    INTEGER :: i_rlstart_c, i_rlend_c
    INTEGER :: i_startblk_c, i_endblk_c
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: cnt, npts

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_coastedge_list'

    DO jg = 1, n_dom

      ! Allocation
      npts = COUNT(p_patch(jg)%edges%refin_ctrl(:,:) == 1)

      ALLOCATE(wave_config(jg)%idx_coastedges(npts),  &
               wave_config(jg)%blk_coastedges(npts),  &
               wave_config(jg)%orient_coastedges(npts))

      ! set up loop over boundary cells (refine_c_ctrl==1)
      i_rlstart_c  = 1
      i_rlend_c    = 1
      i_startblk_c = p_patch(jg)%cells%start_block(i_rlstart_c)
      i_endblk_c   = p_patch(jg)%cells%end_block(i_rlend_c)

      cnt = 0

      DO jb = i_startblk_c, i_endblk_c

        CALL get_indices_c(p_patch(jg), jb, i_startblk_c, i_endblk_c, &
        &                  i_startidx_c, i_endidx_c, i_rlstart_c, i_rlend_c)

        !
        DO jc = i_startidx_c, i_endidx_c

          ! build list of coastline edges (refin_e_ctrl==1)
          ! and store edge orientation.
          !
          DO jce =1,3
            eidx = p_patch(jg)%cells%edge_idx(jc,jb,jce)
            eblk = p_patch(jg)%cells%edge_blk(jc,jb,jce)

            IF (p_patch(jg)%edges%refin_ctrl(eidx,eblk) == 1) THEN
              ! coastline edge found
              cnt = cnt + 1
              wave_config(jg)%idx_coastedges(cnt) = eidx
              wave_config(jg)%blk_coastedges(cnt) = eblk
              wave_config(jg)%orient_coastedges(cnt) = p_patch(jg)%cells%edge_orientation(jc,jb,jce)
            ENDIF
          ENDDO

        ENDDO  !jc

      ENDDO  !jb

      ! Store size of index list
      wave_config(jg)%n_coastedges = cnt
      IF (cnt /= npts) CALL finish(routine, 'mismatch in number of coastal edge points')

    ENDDO ! jg

  END SUBROUTINE init_coastedge_list

END MODULE mo_wave_ext_data_init
