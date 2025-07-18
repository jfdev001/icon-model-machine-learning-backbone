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

! Unit test to ensure the following functions from `mo_slice_array` work as expected:
!   get_3d_general
!   get_2d_general

PROGRAM test_slice_array

  USE mo_kind,        ONLY: wp, xwp, dp, sp
  USE mo_exception,   ONLY: finish
  USE mo_var,         ONLY: t_var
  USE mo_slice_array, ONLY: get_2d_general, get_3d_general

  IMPLICIT NONE

  LOGICAL, PARAMETER :: debug = .FALSE. ! Not yet used
  CHARACTER(*), PARAMETER :: modname = 'test_slice_array'

  TYPE(t_var) :: var
  REAL(dp), TARGET, ALLOCATABLE :: r_data(:,:,:,:,:)
  REAL(sp), TARGET, ALLOCATABLE :: s_data(:,:,:,:,:)
  INTEGER,  TARGET, ALLOCATABLE :: i_data(:,:,:,:,:)
  INTEGER, PARAMETER :: data_shape(5) = (/5,4,3,2,1/)
  INTEGER, PARAMETER :: data_size = PRODUCT(data_shape)
  INTEGER :: init_data(data_size)
  INTEGER :: i

  ALLOCATE( r_data(data_shape(1),data_shape(2),data_shape(3), &
    &              data_shape(4),data_shape(5)) )
  DO i = 1, PRODUCT(data_shape)
    init_data(i) = i
  END DO
  r_data = RESHAPE( REAL(init_data,dp), SHAPE=data_shape)
  s_data = RESHAPE( REAL(init_data+data_size,sp), SHAPE=data_shape)
  i_data = RESHAPE( init_data+2*data_size, SHAPE=data_shape)

  var%r_ptr => r_data
  var%s_ptr => s_data
  var%i_ptr => i_data
  CALL var%set_auxiliary_pointers()

  IF (debug) THEN
    WRITE (*,*) "data_shape: ", data_shape
    !WRITE (*,*) "data: ", init_data
    !WRITE (*,*) "r_ptr: ", var%r_ptr
  END IF

  CALL check_get_ptr(var, r_data, s_data, i_data)
  CALL check_get_3d_general(var, r_data, s_data, i_data)
  CALL check_get_2d_general(var, r_data, s_data, i_data)


  CONTAINS
    SUBROUTINE check_get_ptr(var, r_data, s_data, i_data)
      TYPE(t_var), INTENT(IN) :: var
      REAL(dp), INTENT(IN) :: r_data(:,:,:,:,:)
      REAL(sp), INTENT(IN) :: s_data(:,:,:,:,:)
      INTEGER,  INTENT(IN) :: i_data(:,:,:,:,:)

      IF ( ANY(r_data /= var%get_ptr(1._dp)) ) THEN
        CALL finish(modname,"get_ptr(1._dp)")
      END IF

      IF ( ANY(s_data /= var%get_ptr(1._sp)) ) THEN
        CALL finish(modname,"get_ptr(1._sp)")
      END IF

    #ifdef __SINGLE_PRECISION
      IF ( .NOT. ASSOCIATED(var%get_ptr(1._wp),var%get_ptr(1._sp)) ) THEN
        CALL finish(modname,"wp=sp: get_ptr(1._wp)")
      END IF
      IF ( .NOT. ASSOCIATED(var%get_ptr(1._xwp),var%get_ptr(1._dp)) ) THEN
        CALL finish(modname,"wp=sp: get_ptr(1._xwp)")
      END IF
    #else
      IF ( .NOT. ASSOCIATED(var%get_ptr(1._wp),var%get_ptr(1._dp)) ) THEN
        CALL finish(modname,"wp=dp: get_ptr(1._wp)")
      END IF
      IF ( .NOT. ASSOCIATED(var%get_ptr(1._xwp),var%get_ptr(1._sp)) ) THEN
        CALL finish(modname,"wp=dp: get_ptr(1._xwp)")
      END IF
    #endif

      IF ( ANY(i_data /= var%get_ptr(1)) ) THEN
        CALL finish(modname,"get_ptr(1)")
      END IF

    END SUBROUTINE check_get_ptr

    SUBROUTINE check_get_3d_general(var, r_data, s_data, i_data)
      TYPE(t_var), INTENT(IN) :: var
      REAL(dp), INTENT(IN) :: r_data(:,:,:,:,:)
      REAL(sp), INTENT(IN) :: s_data(:,:,:,:,:)
      INTEGER,  INTENT(IN) :: i_data(:,:,:,:,:)

      ! Local vars
      REAL(dp), POINTER :: r_ptr_3d(:,:,:)
      REAL(sp), POINTER :: s_ptr_3d(:,:,:)
      INTEGER,  POINTER :: i_ptr_3d(:,:,:)

      REAL(wp), POINTER :: wp_ptr_3d(:,:,:)
      REAL(vp), POINTER :: vp_ptr_3d(:,:,:)

      r_ptr_3d => get_3d_general(var%r_ptr, (/3,2/), (/.TRUE.,.FALSE.,.TRUE.,.FALSE.,.FALSE./))
      IF ( ANY(r_data(3,:,2,:,:) /= r_ptr_3d) ) THEN
        CALL finish(modname,"get_3d_general(r_ptr)")
      END IF

      s_ptr_3d => get_3d_general(var%s_ptr, (/3,2/), (/.TRUE.,.FALSE.,.TRUE.,.FALSE.,.FALSE./))
      IF ( ANY(s_data(3,:,2,:,:) /= s_ptr_3d) ) THEN
        CALL finish(modname,"get_3d_general(s_ptr)")
      END IF

      i_ptr_3d => get_3d_general(var%i_ptr, (/3,2/), (/.TRUE.,.FALSE.,.TRUE.,.FALSE.,.FALSE./))
      IF ( ANY(i_data(3,:,2,:,:) /= i_ptr_3d) ) THEN
        CALL finish(modname,"get_3d_general(i_ptr)")
      END IF

      wp_ptr_3d => get_3d_general(var%wp_ptr, (/3,2/), (/.TRUE.,.FALSE.,.TRUE.,.FALSE.,.FALSE./))
#ifdef __SINGLE_PRECISION
      IF ( ANY(s_data(3,:,2,:,:) /= wp_ptr_3d) ) THEN
#else
      IF ( ANY(r_data(3,:,2,:,:) /= wp_ptr_3d) ) THEN
#endif
        CALL finish(modname,"get_3d_general(wp_ptr)")
      END IF

      vp_ptr_3d => get_3d_general(var%vp_ptr, (/3,2/), (/.TRUE.,.FALSE.,.TRUE.,.FALSE.,.FALSE./))
#if defined __SINGLE_PRECISION || defined __MIXED_PRECISION
      IF ( ANY(s_data(3,:,2,:,:) /= vp_ptr_3d) ) THEN
#else
      IF ( ANY(r_data(3,:,2,:,:) /= vp_ptr_3d) ) THEN
#endif
        CALL finish(modname,"get_3d_general(vp_ptr)")
      END IF

    END SUBROUTINE check_get_3d_general

    SUBROUTINE check_get_2d_general(var, r_data, s_data, i_data)
      TYPE(t_var), INTENT(IN) :: var
      REAL(dp), INTENT(IN) :: r_data(:,:,:,:,:)
      REAL(sp), INTENT(IN) :: s_data(:,:,:,:,:)
      INTEGER,  INTENT(IN) :: i_data(:,:,:,:,:)

      ! Local vars
      REAL(dp), POINTER :: r_ptr_2d(:,:)
      REAL(sp), POINTER :: s_ptr_2d(:,:)
      INTEGER,  POINTER :: i_ptr_2d(:,:)

      REAL(wp), POINTER :: wp_ptr_2d(:,:)
      REAL(vp), POINTER :: vp_ptr_2d(:,:)

      r_ptr_2d => get_2d_general(var%r_ptr, (/3,2,2/), (/.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE./))
      IF ( ANY(r_data(3,:,2,2,:) /= r_ptr_2d) ) THEN
        CALL finish(modname,"get_2d_general(r_ptr)")
      END IF

      s_ptr_2d => get_2d_general(var%s_ptr, (/3,2,2/), (/.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE./))
      IF ( ANY(s_data(3,:,2,2,:) /= s_ptr_2d) ) THEN
        CALL finish(modname,"get_2d_general(s_ptr)")
      END IF

      i_ptr_2d => get_2d_general(var%i_ptr, (/3,2,2/), (/.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE./))
      IF ( ANY(i_data(3,:,2,2,:) /= i_ptr_2d) ) THEN
        CALL finish(modname,"get_2d_general(i_ptr)")
      END IF

      wp_ptr_2d => get_2d_general(var%wp_ptr, (/3,2,2/), (/.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE./))
#ifdef __SINGLE_PRECISION
      IF ( ANY(s_data(3,:,2,2,:) /= wp_ptr_2d) ) THEN
#else
      IF ( ANY(r_data(3,:,2,2,:) /= wp_ptr_2d) ) THEN
#endif
        CALL finish(modname,"get_2d_general(wp_ptr)")
      END IF

      vp_ptr_2d => get_2d_general(var%vp_ptr, (/3,2,2/), (/.TRUE.,.FALSE.,.TRUE.,.TRUE.,.FALSE./))
#if defined __SINGLE_PRECISION || defined __MIXED_PRECISION
      IF ( ANY(s_data(3,:,2,2,:) /= vp_ptr_2d) ) THEN
#else
      IF ( ANY(r_data(3,:,2,2,:) /= vp_ptr_2d) ) THEN
#endif
        CALL finish(modname,"get_2d_general(vp_ptr)")
      END IF

    END SUBROUTINE check_get_2d_general

END PROGRAM test_slice_array
