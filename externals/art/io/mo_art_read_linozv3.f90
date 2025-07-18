!---------------------------------------------------------
!---Read LINOZ tables from 'DPMLDT_NOYP.dat' for linearized Ozonechemistry
!---------------------------------------------------------
!
! mo_art_read_linozv3 
! Adapation of mo_art_read_linoz to LinozV3
!
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

MODULE mo_art_read_linozv3
#ifdef __ICON_ART
  USE mo_art_chem_data,                 ONLY: t_art_linozv3

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_linoz_readv3
  PUBLIC  :: art_linoz_deallocatev3

CONTAINS

SUBROUTINE art_linoz_readv3(art_linozv3, nproma, nlev)

!<
! Adaption of SUBROUTINE art_linoz_read to LinozV3
!>

  IMPLICIT NONE 

  INTEGER, INTENT(in) :: &
    &  nproma, nlev              !< dimensions of arrays


  TYPE(t_art_linozv3), INTENT(inout)    :: &
    &  art_linozv3                  !< Pointer to ART chem fields

  CHARACTER(LEN=50) TITEL
  CHARACTER(LEN=120) linozinifile1
  CHARACTER(LEN=120) linozinifile2

  INTEGER n,m,k,l


  IF (.NOT. ALLOCATED (art_linozv3%tparm_max)) ALLOCATE(art_linozv3%tparm_max(25,18,12,9))
  linozinifile1 = 'Linoz2000v3maxnew.txt'

  OPEN(61,file=linozinifile1, STATUS = 'old')
  READ(61,901) TITEL
  DO n=1,9
    READ(61,901) TITEL
    DO m=1,12
      DO k=1,18
        READ(61,902) (art_linozv3%tparm_max(l,k,m,n),l=1,25)
      ENDDO
    ENDDO
  ENDDO
 
  901 FORMAT(A50)
  902 FORMAT(20x,6e10.3/(8e10.3))
  CLOSE(61)


  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab1)) ALLOCATE(art_linozv3%linozv3_tab1(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab2)) ALLOCATE(art_linozv3%linozv3_tab2(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab3)) ALLOCATE(art_linozv3%linozv3_tab3(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab4)) ALLOCATE(art_linozv3%linozv3_tab4(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab5)) ALLOCATE(art_linozv3%linozv3_tab5(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab6)) ALLOCATE(art_linozv3%linozv3_tab6(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab7)) ALLOCATE(art_linozv3%linozv3_tab7(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab8)) ALLOCATE(art_linozv3%linozv3_tab8(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab9)) ALLOCATE(art_linozv3%linozv3_tab9(nproma,nlev))
  art_linozv3%is_init = .TRUE.

  IF (.NOT. ALLOCATED (art_linozv3%tparm_min)) ALLOCATE(art_linozv3%tparm_min(25,18,12,9))
  linozinifile2 = 'Linoz2000v3minnew.txt'

  OPEN(62,file=linozinifile2, STATUS = 'old')
  READ(62,903) TITEL
  DO n=1,9
    READ(62,903) TITEL
    DO m=1,12
      DO k=1,18
        READ(62,904) (art_linozv3%tparm_min(l,k,m,n),l=1,25)
      ENDDO
    ENDDO
  ENDDO

  903 FORMAT(A50)
  904 FORMAT(20x,6e10.3/(8e10.3))
  CLOSE(61)

  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab1)) ALLOCATE(art_linozv3%linozv3_tab1(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab2)) ALLOCATE(art_linozv3%linozv3_tab2(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab3)) ALLOCATE(art_linozv3%linozv3_tab3(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab4)) ALLOCATE(art_linozv3%linozv3_tab4(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab5)) ALLOCATE(art_linozv3%linozv3_tab5(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab6)) ALLOCATE(art_linozv3%linozv3_tab6(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab7)) ALLOCATE(art_linozv3%linozv3_tab7(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab8)) ALLOCATE(art_linozv3%linozv3_tab8(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linozv3%linozv3_tab9)) ALLOCATE(art_linozv3%linozv3_tab9(nproma,nlev))
  art_linozv3%is_init = .TRUE.
 
END SUBROUTINE art_linoz_readv3

SUBROUTINE art_linoz_deallocatev3(art_linozv3)
!<
! Adapation of SUBROUTINE linoz_read to LinozV3
!>
  TYPE(t_art_linozv3), INTENT(inout)    :: &
    &  art_linozv3                    !< Pointer to ART chem fields

  DEALLOCATE(art_linozv3%linozv3_tab1)
  DEALLOCATE(art_linozv3%linozv3_tab2)
  DEALLOCATE(art_linozv3%linozv3_tab3)
  DEALLOCATE(art_linozv3%linozv3_tab4)
  DEALLOCATE(art_linozv3%linozv3_tab5)
  DEALLOCATE(art_linozv3%linozv3_tab6)
  DEALLOCATE(art_linozv3%linozv3_tab7)
  DEALLOCATE(art_linozv3%linozv3_tab8)
  DEALLOCATE(art_linozv3%linozv3_tab9)
  DEALLOCATE(art_linozv3%tparm_max)
  DEALLOCATE(art_linozv3%linozv3_tab1)
  DEALLOCATE(art_linozv3%linozv3_tab2)
  DEALLOCATE(art_linozv3%linozv3_tab3)
  DEALLOCATE(art_linozv3%linozv3_tab4)
  DEALLOCATE(art_linozv3%linozv3_tab5)
  DEALLOCATE(art_linozv3%linozv3_tab6)
  DEALLOCATE(art_linozv3%linozv3_tab7)
  DEALLOCATE(art_linozv3%linozv3_tab8)
  DEALLOCATE(art_linozv3%linozv3_tab9)
  DEALLOCATE(art_linozv3%tparm_min)
  art_linozv3%is_init = .FALSE.

END SUBROUTINE art_linoz_deallocatev3
#endif
END MODULE mo_art_read_linozv3
