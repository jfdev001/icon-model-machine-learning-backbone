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

program test_index_list

  use mo_index_list
  use mo_kind

  implicit none

  integer,     parameter :: n  = 200000
  integer,     parameter :: nb = 7
  integer(i4), target    :: conditions(n,nb)
  integer,     target    :: indices(n,nb), dev_indices(n,nb)
  integer                :: nvalid(nb), dev_nvalid(nb)
  real                   :: harvest(n,nb)

  integer :: i, b, ic
  logical :: lacc

#ifdef _OPENACC
  lacc = .TRUE.
#else
  lacc = .FALSE.
#endif

  call random_number(harvest)
  conditions = int(harvest * 2)

  !$ACC DATA COPYIN(conditions) CREATE(dev_indices, dev_nvalid)

  ! Test the non-batched version

  nvalid(1) = 0
  do i = 1, n
    if (conditions(i,1) /= 0) then
      nvalid(1) = nvalid(1) + 1
      indices(nvalid(1),1) = i
    end if
  end do


  call generate_index_list(conditions(:,1), dev_indices(:,1), 1, n, dev_nvalid(1), lacc, 1)
  !$ACC UPDATE HOST(dev_indices(:,1)) ASYNC(1)
  !$ACC WAIT(1)

  print *, "CHECK NON-BATCHED: ", nvalid(1) == dev_nvalid(1), all(indices(:nvalid(1),1) == dev_indices(:nvalid(1),1))

  ! Test the non-batched async version

  call generate_index_list(conditions(:,1), dev_indices(:,1), 1, n, dev_nvalid(1), lacc, 1, opt_acc_copy_to_host = .FALSE.)
  !$ACC UPDATE HOST(dev_indices(:,1), dev_nvalid(1)) ASYNC(1)
  !$ACC WAIT(1)

  print *, "CHECK NON-BATCHED ASYNC: ", nvalid(1) == dev_nvalid(1), all(indices(:nvalid(1),1) == dev_indices(:nvalid(1),1))

  ! Test the non-batched version with a shift

  ic = 42
  nvalid(1) = 0
  do i = ic, n
    if (conditions(i,1) /= 0) then
      nvalid(1) = nvalid(1) + 1
      indices(nvalid(1),1) = i
    end if
  end do

  call generate_index_list(conditions(:,1), dev_indices(:,1), ic, n, dev_nvalid(1), lacc, 1)
  !$ACC UPDATE HOST(dev_indices(:,1)) ASYNC(1)
  !$ACC WAIT(1)

  print *, "CHECK NON-BATCHED SHIFTED: ", nvalid(1) == dev_nvalid(1), all(indices(:nvalid(1),1) == dev_indices(:nvalid(1),1))

  ! Test the batched version

  indices(:,:) = 0
  !$ACC KERNELS
  dev_indices(:,:) = 0
  !$ACC END KERNELS

  nvalid = 0
  do b = 1, nb
    do i = 1, n
      if (conditions(i,b) /= 0) then
        nvalid(b) = nvalid(b) + 1
        indices(nvalid(b),b) = i
      end if
    end do
  end do

  call generate_index_list_batched(conditions, dev_indices, 1, n, dev_nvalid, &
    &   lacc=.TRUE., opt_acc_async_queue=1)
  !$ACC UPDATE HOST(dev_indices, dev_nvalid) ASYNC(1)
  !$ACC WAIT(1)

  print *, "CHECK BATCHED: ", all(nvalid == dev_nvalid), all(indices == dev_indices)

  ! Test the batched shifted version

  indices(:,:) = 0
  !$ACC KERNELS
  dev_indices(:,:) = 0
  !$ACC END KERNELS

  ic = 142
  nvalid = 0
  do b = 1, nb
    do i = ic, n
      if (conditions(i,b) /= 0) then
        nvalid(b) = nvalid(b) + 1
        indices(nvalid(b),b) = i
      end if
    end do
  end do

  call generate_index_list_batched(conditions, dev_indices, ic, n, dev_nvalid, lacc, 1)
  !$ACC UPDATE HOST(dev_indices, dev_nvalid) ASYNC(1)
  !$ACC WAIT(1)

  print *, "CHECK BATCHED SHIFTED: ", all(nvalid == dev_nvalid), all(indices == dev_indices)

  !$ACC END DATA


end program test_index_list
