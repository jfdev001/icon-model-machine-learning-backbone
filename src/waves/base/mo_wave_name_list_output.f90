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
! wave-specific name list output routines

MODULE mo_wave_name_list_output

  USE mo_exception,      ONLY: finish
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_mpi,            ONLY: p_bcast, p_comm_work_2_io, my_process_is_io
  USE mo_wave_config,    ONLY: t_wave_config, wave_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: replicate_wave_data_on_io_procs

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_name_list_output'

CONTAINS

  !
  ! Replicate wave-specific data on IO procs
  !
  SUBROUTINE replicate_wave_data_on_io_procs (n_dom_out, bcast_root)

    INTEGER, INTENT(IN) :: n_dom_out   ! number of output domains
    INTEGER, INTENT(IN) :: bcast_root  ! Broadcast root for intercommunicator broadcasts
    ! local
    CHARACTER(*), PARAMETER :: routine = modname//'::replicate_wave_data_on_io_procs'
    INTEGER :: jg
    INTEGER :: ist
    TYPE(t_wave_config), POINTER :: wc =>NULL()     ! convenience pointer

    DO jg=1,n_dom_out
      ! convenience pointer
      wc => wave_config(jg)
      IF (my_process_is_io()) THEN
        ALLOCATE(wc%freqs(wc%nfreqs), wc%dirs(wc%ndirs), stat=ist)
        IF (ist/=SUCCESS) CALL finish(routine, "allocation for wc%freqs and wc%dirs failed on IO PE")
      ENDIF
      CALL p_bcast(wc%freqs(:), bcast_root, p_comm_work_2_io)
      CALL p_bcast(wc%dirs(:) , bcast_root, p_comm_work_2_io)
    ENDDO

  END SUBROUTINE replicate_wave_data_on_io_procs

END MODULE mo_wave_name_list_output
