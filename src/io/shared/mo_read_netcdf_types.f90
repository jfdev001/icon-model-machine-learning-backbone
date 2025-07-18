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

! This module provides data structures for reading a NetCDF file in a distributed way.

MODULE mo_read_netcdf_types

  USE mo_kind, ONLY: wp
  USE mo_communication_types, ONLY: t_comm_pattern

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_distrib_read_data
  PUBLIC :: t_alloc_2d
  PUBLIC :: t_alloc_2d_int
  PUBLIC :: t_alloc_3d
  PUBLIC :: t_alloc_3d_int

  TYPE t_distrib_read_data
    INTEGER :: basic_data_index = -1
    CLASS(t_comm_pattern), POINTER :: pat => NULL()
  END TYPE t_distrib_read_data

  TYPE t_alloc_2d
    REAL(wp), ALLOCATABLE :: a(:,:)
  END TYPE t_alloc_2d

  TYPE t_alloc_2d_int
    INTEGER, ALLOCATABLE :: a(:,:)
  END TYPE t_alloc_2d_int

  TYPE t_alloc_3d
    REAL(wp), ALLOCATABLE :: a(:,:,:)
  END TYPE t_alloc_3d

  TYPE t_alloc_3d_int
    INTEGER, ALLOCATABLE :: a(:,:,:)
  END TYPE t_alloc_3d_int

END MODULE mo_read_netcdf_types
