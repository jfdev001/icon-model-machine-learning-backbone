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

MODULE mo_slice_array

  USE mo_kind,               ONLY: dp, sp
  USE mo_exception,          ONLY: message, finish, message_text

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_2d_general, get_3d_general

  INTERFACE get_2d_general
    MODULE PROCEDURE get_r2d_general
    MODULE PROCEDURE get_s2d_general
    MODULE PROCEDURE get_i2d_general
  END INTERFACE get_2d_general

  INTERFACE get_3d_general
    MODULE PROCEDURE get_r3d_general
    MODULE PROCEDURE get_s3d_general
    MODULE PROCEDURE get_i3d_general
  END INTERFACE get_3d_general

  CHARACTER(*), PARAMETER :: modname = "mo_slice_array"

CONTAINS

    ! ---------------------------------------------------------------
  !! Usage:
  ! squash_dims := mask for which dimensions to compress
  ! indexes := index corresponding to each compressed dimension
  ! We use the following notation:
  !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
  !    var%s_ptr(:,1,:,:,2) == get_3d_general(var%s_ptr, in_indexes=(1,2), squash_dims=(F,T,F,F,T))
  !
  !    var%r_ptr(:,1,3,2,:) == get_2d_general(var%r_ptr, in_indexes=(1,3,2), squash_dims=(F,T,T,T,F))
  !    var%s_ptr(:,1,3,:,2) == get_2d_general(var%s_ptr, in_indexes=(1,3,2), squash_dims=(F,T,T,F,T))
  FUNCTION get_r3d_general(ptr_5d, in_indexes, squash_dims) RESULT(ptr)

    REAL(dp), POINTER :: ptr_5d (:,:,:,:,:)
    INTEGER, INTENT(IN) :: in_indexes(2)
    LOGICAL, INTENT(IN) :: squash_dims(5)

    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":get_r3d_general"
    INTEGER :: nsquashed_dims
    REAL(dp), POINTER :: ptr    (:,:,:)
    REAL(dp), POINTER :: ptr_4d (:,:,:,:)
    INTEGER :: indexes(SIZE(in_indexes))

    ! Since ptr compressed in reverse, we must also reverse indices
    ! Eg: We want the following notation:
    !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
    indexes(:) = in_indexes(SIZE(in_indexes):1:-1)

    nsquashed_dims=0
    IF (squash_dims(5)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(4)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr => ptr_4d(:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(3)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,indexes(nsquashed_dims+1),:,:)
      ELSE
        ptr => ptr_4d(:,:,indexes(nsquashed_dims+1),:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(2)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_4d(:,indexes(nsquashed_dims+1),:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(1)) THEN
      IF (nsquashed_dims==0) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_4d(indexes(nsquashed_dims+1),:,:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF

    IF (nsquashed_dims/=2) THEN ! Something went wrong :)
      CALL finish(routine, "incorrect state")
    ENDIF

  END FUNCTION get_r3d_general

  FUNCTION get_s3d_general(ptr_5d, in_indexes, squash_dims) RESULT(ptr)

    REAL(sp), POINTER, INTENT(IN) :: ptr_5d (:,:,:,:,:)
    INTEGER, INTENT(IN) :: in_indexes(2)
    LOGICAL, INTENT(IN) :: squash_dims(5)

    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":get_s3d_general"
    INTEGER :: nsquashed_dims
    REAL(sp), POINTER :: ptr    (:,:,:)
    REAL(sp), POINTER :: ptr_4d (:,:,:,:)
    INTEGER :: indexes(SIZE(in_indexes))

    ! Since ptr compressed in reverse, we must also reverse indices
    ! Eg: We want the following notation:
    !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
    indexes(:) = in_indexes(SIZE(in_indexes):1:-1)

    nsquashed_dims=0
    IF (squash_dims(5)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(4)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr => ptr_4d(:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(3)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,indexes(nsquashed_dims+1),:,:)
      ELSE
        ptr => ptr_4d(:,:,indexes(nsquashed_dims+1),:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(2)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_4d(:,indexes(nsquashed_dims+1),:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(1)) THEN
      IF (nsquashed_dims==0) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_4d(indexes(nsquashed_dims+1),:,:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF

    IF (nsquashed_dims/=2) THEN ! Something went wrong :)
      CALL finish(routine, "incorrect state")
    ENDIF

  END FUNCTION get_s3d_general

  FUNCTION get_i3d_general(ptr_5d, in_indexes, squash_dims) RESULT(ptr)

    INTEGER, POINTER, INTENT(IN) :: ptr_5d (:,:,:,:,:)
    INTEGER, INTENT(IN) :: in_indexes(2)
    LOGICAL, INTENT(IN) :: squash_dims(5)

    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":get_i3d_general"
    INTEGER :: nsquashed_dims
    INTEGER, POINTER :: ptr    (:,:,:)
    INTEGER, POINTER :: ptr_4d (:,:,:,:)
    INTEGER :: indexes(SIZE(in_indexes))

    ! Since ptr compressed in reverse, we must also reverse indices
    ! Eg: We want the following notation:
    !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
    indexes(:) = in_indexes(SIZE(in_indexes):1:-1)

    nsquashed_dims=0
    IF (squash_dims(5)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(4)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr => ptr_4d(:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(3)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,indexes(nsquashed_dims+1),:,:)
      ELSE
        ptr => ptr_4d(:,:,indexes(nsquashed_dims+1),:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(2)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_4d(:,indexes(nsquashed_dims+1),:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(1)) THEN
      IF (nsquashed_dims==0) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_4d(indexes(nsquashed_dims+1),:,:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF

    IF (nsquashed_dims/=2) THEN ! Something went wrong :)
      CALL finish(routine, "incorrect state")
    ENDIF

  END FUNCTION get_i3d_general

  !-----------------------------------------------------------------------------

  !! Usage:
  ! squash_dims := mask for which dimensions to compress
  ! indexes := index corresponding to each compressed dimension
  ! We use the following notation:
  !    var%r_ptr(:,1,3,2,:) == get_2d_general(var%r_ptr, in_indexes=(1,3,2), squash_dims=(F,T,T,T,F))
  !    var%s_ptr(:,1,3,:,2) == get_2d_general(var%s_ptr, in_indexes=(1,3,2), squash_dims=(F,T,T,F,T))
  !
  !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
  !    var%s_ptr(:,1,:,:,2) == get_3d_general(var%s_ptr, in_indexes=(1,2), squash_dims=(F,T,F,F,T))
  ! May be able to simplify by calling get_r3d
  FUNCTION get_r2d_general(ptr_5d, in_indexes, squash_dims) RESULT(ptr)

    REAL(dp), POINTER, INTENT(IN) :: ptr_5d (:,:,:,:,:)
    INTEGER, INTENT(IN) :: in_indexes(3)
    LOGICAL, INTENT(IN) :: squash_dims(5)

    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":get_r2d_general"
    INTEGER :: nsquashed_dims
    REAL(dp), POINTER :: ptr    (:,:)
    REAL(dp), POINTER :: ptr_3d (:,:,:)
    REAL(dp), POINTER :: ptr_4d (:,:,:,:)
    INTEGER :: indexes(SIZE(in_indexes))

    ! Since ptr compressed in reverse, we must also reverse indices
    ! Eg: We want the following notation:
    !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
    indexes(:) = in_indexes(SIZE(in_indexes):1:-1)

    nsquashed_dims=0
    IF (squash_dims(5)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(4)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr_3d => ptr_4d(:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(3)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,indexes(nsquashed_dims+1),:,:)
      ELSE IF (nsquashed_dims==1) THEN
        ptr_3d => ptr_4d(:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr => ptr_3d(:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(2)) THEN
      IF (nsquashed_dims==0) THEN
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(:,indexes(nsquashed_dims+1),:,:,:)
      ELSE IF (nsquashed_dims==1) THEN
        ptr_3d => ptr_4d(:,indexes(nsquashed_dims+1),:,:)
      ELSE
        ptr => ptr_3d(:,indexes(nsquashed_dims+1),:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(1)) THEN
      IF (nsquashed_dims==0) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(indexes(nsquashed_dims+1),:,:,:)
      ELSE IF (nsquashed_dims==1) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        ptr_3d => ptr_4d(indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_3d(indexes(nsquashed_dims+1),:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF

    IF (nsquashed_dims/=3) THEN ! Something went wrong :)
      CALL finish(routine, "incorrect state")
    ENDIF

    END FUNCTION get_r2d_general

  ! May be able to simplify by calling get_s3d
  FUNCTION get_s2d_general(ptr_5d, in_indexes, squash_dims) RESULT(ptr)

    REAL(sp), POINTER, INTENT(IN):: ptr_5d (:,:,:,:,:)
    INTEGER, INTENT(IN) :: in_indexes(3)
    LOGICAL, INTENT(IN) :: squash_dims(5)

    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":get_s2d_general"
    INTEGER :: nsquashed_dims
    REAL(sp), POINTER :: ptr    (:,:)
    REAL(sp), POINTER :: ptr_3d (:,:,:)
    REAL(sp), POINTER :: ptr_4d (:,:,:,:)
    INTEGER :: indexes(SIZE(in_indexes))

    ! Since ptr compressed in reverse, we must also reverse indices
    ! Eg: We want the following notation:
    !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
    indexes(:) = in_indexes(SIZE(in_indexes):1:-1)

    nsquashed_dims=0
    IF (squash_dims(5)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(4)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr_3d => ptr_4d(:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(3)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,indexes(nsquashed_dims+1),:,:)
      ELSE IF (nsquashed_dims==1) THEN
        ptr_3d => ptr_4d(:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr => ptr_3d(:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(2)) THEN
      IF (nsquashed_dims==0) THEN
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(:,indexes(nsquashed_dims+1),:,:,:)
      ELSE IF (nsquashed_dims==1) THEN
        ptr_3d => ptr_4d(:,indexes(nsquashed_dims+1),:,:)
      ELSE
        ptr => ptr_3d(:,indexes(nsquashed_dims+1),:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(1)) THEN
      IF (nsquashed_dims==0) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(indexes(nsquashed_dims+1),:,:,:)
      ELSE IF (nsquashed_dims==1) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        ptr_3d => ptr_4d(indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_3d(indexes(nsquashed_dims+1),:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF

    IF (nsquashed_dims/=3) THEN ! Something went wrong :)
      CALL finish(routine, "incorrect state")
    ENDIF

  END FUNCTION get_s2d_general

  FUNCTION get_i2d_general(ptr_5d, in_indexes, squash_dims) RESULT(ptr)

    INTEGER, POINTER, INTENT(IN):: ptr_5d (:,:,:,:,:)
    INTEGER, INTENT(IN) :: in_indexes(3)
    LOGICAL, INTENT(IN) :: squash_dims(5)

    ! Local vars
    CHARACTER(*), PARAMETER :: routine = modname//":get_i2d_general"
    INTEGER :: nsquashed_dims
    INTEGER, POINTER :: ptr    (:,:)
    INTEGER, POINTER :: ptr_3d (:,:,:)
    INTEGER, POINTER :: ptr_4d (:,:,:,:)
    INTEGER :: indexes(SIZE(in_indexes))

    ! Since ptr compressed in reverse, we must also reverse indices
    ! Eg: We want the following notation:
    !    var%r_ptr(:,1,3,:,:) == get_3d_general(var%r_ptr, in_indexes=(1,3), squash_dims=(F,T,T,F,F))
    indexes(:) = in_indexes(SIZE(in_indexes):1:-1)

    nsquashed_dims=0
    IF (squash_dims(5)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(4)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr_3d => ptr_4d(:,:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(3)) THEN
      IF (nsquashed_dims==0) THEN
        ptr_4d => ptr_5d(:,:,indexes(nsquashed_dims+1),:,:)
      ELSE IF (nsquashed_dims==1) THEN
        ptr_3d => ptr_4d(:,:,indexes(nsquashed_dims+1),:)
      ELSE
        ptr => ptr_3d(:,:,indexes(nsquashed_dims+1))
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(2)) THEN
      IF (nsquashed_dims==0) THEN
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(:,indexes(nsquashed_dims+1),:,:,:)
      ELSE IF (nsquashed_dims==1) THEN
        ptr_3d => ptr_4d(:,indexes(nsquashed_dims+1),:,:)
      ELSE
        ptr => ptr_3d(:,indexes(nsquashed_dims+1),:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF
    IF (squash_dims(1)) THEN
      IF (nsquashed_dims==0) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        !ptr_4d => ptr_5d(indexes(nsquashed_dims+1),:,:,:)
      ELSE IF (nsquashed_dims==1) THEN ! Something went wrong :)
        CALL finish(routine, "incorrect state")
        ptr_3d => ptr_4d(indexes(nsquashed_dims+1),:,:,:)
      ELSE
        ptr => ptr_3d(indexes(nsquashed_dims+1),:,:)
      ENDIF
      nsquashed_dims = nsquashed_dims + 1
    ENDIF

    IF (nsquashed_dims/=3) THEN ! Something went wrong :)
      CALL finish(routine, "incorrect state")
    ENDIF

    END FUNCTION get_i2d_general

END MODULE mo_slice_array
