!
! mo_art_read_opt_props
! This module provides reading routines of optical properties
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

MODULE mo_art_read_opt_props
  USE mo_kind,                          ONLY: wp
  USE mo_grid_config,                   ONLY: nroot
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mtime,                            ONLY: newDatetime
  USE mo_key_value_store,               ONLY: t_key_value_store
  ! ART
  USE mo_art_impl_constants,            ONLY: IART_XMLTAGLEN,            &
    &                                         IART_VARNAMELEN
  USE mo_art_config,                    ONLY: IART_PATH_LEN,             &
    &                                         art_config
  USE netcdf,                           ONLY: NF90_NOERR,                &
    &                                         NF90_NOWRITE,              &
    &                                         nf90_strerror,             &
    &                                         nf90_open,                 &
    &                                         nf90_close,                &
    &                                         nf90_inq_dimid,            &
    &                                         nf90_inq_varid,            &
    &                                         nf90_inquire_dimension,    &
    &                                         nf90_get_var

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: routine = 'mo_art_read_opt_props'

  public :: art_open_optprops_nc
  public :: art_close_optprops_nc
  
  PUBLIC :: art_read_prognostic_optprops
  PUBLIC :: art_read_aeronet_optprops
  PUBLIC :: art_read_satellite_optprops
  PUBLIC :: art_read_singlewave_optprops
  PUBLIC :: art_read_diagnostic_optprops
  PUBLIC :: art_read_dia_factors
  PUBLIC :: art_read_vardia_optprops

CONTAINS

  !
  ! private helper subroutines
  !
  
  SUBROUTINE art_read_optprops_ext_ssa_asy_bsc(ncid, var_name, wave_dim_name, ext_coeff, ssa_coeff, asy_coeff, bsc_coeff)
  !<
  ! SUBROUTINE art_read_optprops_ext_ssa_asy_bsc
  ! This subroutine reads metadata for optical parameters of aerosols
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-01-10
  ! Modifications:
  !>
    !arguments
    INTEGER, INTENT(IN)               :: ncid
    CHARACTER(LEN=*), INTENT(IN)      :: var_name, wave_dim_name
    REAL(wp), POINTER, INTENT(INOUT)  :: ext_coeff(:),  &
                                      &  ssa_coeff(:),  &
                                      &  asy_coeff(:),  &
                                      &  bsc_coeff(:)
    !local variables
    INTEGER          :: ierror, dimid, varid
    INTEGER          :: i, numwaves
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: temp_coeff
    CHARACTER(LEN=*), PARAMETER :: sroutine=TRIM(routine)//':art_read_optprops'
    CHARACTER(LEN=256) :: ssa_text

    ierror = nf90_inq_dimid(ncid, name=TRIM(wave_dim_name), dimid=dimid)
    IF (ierror .eq. NF90_NOERR) THEN
      ierror = nf90_inquire_dimension(ncid, dimid=dimid, len=numwaves)
    ENDIF
    
    IF (ierror /= NF90_NOERR) THEN
      CALL finish(TRIM(sroutine),'Could not read dimension "'//TRIM(wave_dim_name)//'" from the NetCDF optical properties database\\n'//nf90_strerror(ierror))
    ENDIF
    
    ! Checking if numwaves corresponds to length of first dimension of coefficient-fields
    ierror = 0
    IF (ASSOCIATED(ext_coeff)) THEN
      IF (numwaves /= size(ext_coeff,1)) THEN
        ierror = ierror+1
      END IF
    END IF
    IF (ASSOCIATED(ssa_coeff)) THEN
      IF (numwaves /= size(ssa_coeff,1)) THEN
        ierror = ierror+2
      END IF
    END IF
    IF (ASSOCIATED(asy_coeff)) THEN
      IF (numwaves /= size(asy_coeff,1)) THEN
        ierror = ierror+4
      END IF
    END IF
    IF (ASSOCIATED(bsc_coeff)) THEN
      IF (numwaves /= size(bsc_coeff,1)) THEN
        ierror = ierror+8
      END IF
    END IF
    IF (ierror /= 0) THEN
      WRITE(message_text,'(A,4I1)') 'coeff-field-dim doesnt match number of '// &
        &                           'wavebands in NetCDF file, error-code=',    &
        &                           MOD(ierror/8,2), MOD(ierror/4,2),           &
        &                           MOD(ierror/2,2), MOD(ierror,2)
      CALL finish(TRIM(sroutine),message_text)
    ENDIF

    !create temp container for ext/ssa/asy/bsc coefficients
    ALLOCATE( temp_coeff(1:numwaves, 1:4) )
    ierror = nf90_inq_varid(ncid, name=TRIM(var_name), varid=varid)
    IF (ierror .eq. NF90_NOERR) THEN
      ierror = nf90_get_var(ncid, varid, temp_coeff(:,:), start=(/1,1/), count=(/numwaves,4/))
    END IF
    IF (ierror /= NF90_NOERR) THEN
      DEALLOCATE( temp_coeff )
      CALL finish(TRIM(sroutine),'Could not read coefficients for '//TRIM(var_name)//'\\n'//nf90_strerror(ierror) )
    ENDIF

    !copy read data to extinction coefficients
    IF (ASSOCIATED(ext_coeff)) THEN
      DO i=1,numwaves
        ext_coeff(i) = temp_coeff(i,1)
      END DO
    END IF

    !copy read data to single scattering albedo coefficients
    ierror = 0
    IF (ASSOCIATED(ssa_coeff)) THEN
      DO i=1,numwaves
        ssa_coeff(i) = temp_coeff(i,2)
        IF (ssa_coeff(i)<0.0_wp .or. ssa_coeff(i)>1.0_wp) THEN
           WRITE(ssa_text,'(A36,F10.5,A4,I3)') 'found SSA coefficient out of range: ',ssa_coeff(i),' at ',i
           CALL message(TRIM(sroutine),ssa_text)
          ierror = ierror + 1
        END IF
      END DO
    END IF
    IF (ierror /= 0) THEN
      DEALLOCATE( temp_coeff )
      CALL message(TRIM(sroutine),'found SSA coefficient(s) out of range:'//message_text)
      CALL finish(TRIM(sroutine),'At least one ssa coefficient is out of valid range [0.0 ... 1.0] for mode '//TRIM(var_name))
    END IF

    !copy read data to asymmetry coefficients
    IF (ASSOCIATED(asy_coeff)) THEN
      DO i=1,numwaves
        asy_coeff(i) = temp_coeff(i,3)
      END DO
    END IF

    !copy read data to extinction coefficients
    IF (ASSOCIATED(bsc_coeff)) THEN
      DO i=1,numwaves
        bsc_coeff(i) = temp_coeff(i,4)
      END DO
    END IF

    DEALLOCATE( temp_coeff )
  END SUBROUTINE art_read_optprops_ext_ssa_asy_bsc

  !
  ! public subroutines
  !

  SUBROUTINE art_open_optprops_nc(jg, ncid)
  !<
  ! SUBROUTINE art_open_optprops_nc
  ! Opens the NetCDF file, which contains the prescribed optical properties
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-04-07
  !>
    !arguments
    INTEGER, INTENT(IN)  ::         jg       !< domain ID
    INTEGER, INTENT(OUT) ::         ncid     !< access identifier of the NetCDF file
    !local variable
    INTEGER ::                      ierror
    CHARACTER(LEN=IART_PATH_LEN) :: filepath !< path to the NetCDF file

    filepath = TRIM(art_config(jg)%cart_opt_props_nc)
    ierror = nf90_open(TRIM(filepath), NF90_NOWRITE, ncid)
    IF (ierror /= NF90_NOERR) THEN
      CALL finish('mo_art_read_opt_props:art_open_optprops_nc','Could not open NetCDF file '//TRIM(filepath)//'\\n'//nf90_strerror(ierror))
    ENDIF

  END SUBROUTINE art_open_optprops_nc

  SUBROUTINE art_close_optprops_nc(ncid)
  !<
  ! SUBROUTINE art_close_optprops_nc
  ! Closes the NetCDF file, which contains the prescribed optical properties
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-04-07
  !>
    !arguments
    INTEGER, INTENT(IN) :: ncid    !< access identifier of the NetCDF file
    !local variable
    INTEGER :: ierror

    ierror = nf90_close(ncid)
    IF (ierror /= NF90_NOERR) THEN
      CALL finish('mo_art_read_opt_props:art_open_optprops_nc','Could not close NetCDF file.'//'\\n'//nf90_strerror(ierror))
    END IF

  END SUBROUTINE art_close_optprops_nc

  SUBROUTINE art_read_prognostic_optprops(ncid, mode_name, ext_coeff, ssa_coeff, asy_coeff)
  !<
  ! SUBROUTINE art_read_prognostic_optprops
  ! Reads only extinction, SSA and backscatter coefficients for the prognostic radiation calculation of aerosols.
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-01-10
  !>
    !arguments
    INTEGER, INTENT(IN)               :: ncid
    CHARACTER(LEN=*), INTENT(IN)      :: mode_name
    REAL(wp), POINTER, INTENT(INOUT)  :: ext_coeff(:),  &
                                      &  ssa_coeff(:),  &
                                      &  asy_coeff(:)
    !local arguments
    REAL(wp), POINTER :: unused_coeff(:)

    unused_coeff => NULL() !< mark this coefficient unused by pointing to null (de-association)

    CALL art_read_optprops_ext_ssa_asy_bsc(ncid, &
        &                  TRIM(mode_name)//'_prognostic', 'wave_prognostic', &
        &                  ext_coeff, ssa_coeff, asy_coeff, unused_coeff)
  END SUBROUTINE art_read_prognostic_optprops

  SUBROUTINE art_read_aeronet_optprops(ncid, mode_name, ext_coeff)
  !<
  ! SUBROUTINE art_read_prognostic_optprops
  ! Reads only extinction coefficients for the diagnostic radiation calculation of aerosols at AERONET wavelengths.
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-01-10
  !>
    IMPLICIT NONE
    !arguments
    INTEGER, INTENT(IN)               :: ncid
    CHARACTER(LEN=*), INTENT(IN)      :: mode_name
    REAL(wp), POINTER, INTENT(INOUT)  :: ext_coeff(:)
    !local arguments
    REAL(wp), POINTER :: unused_coeff(:)

    unused_coeff => NULL() !< mark this coefficient unused by pointing to null (de-association)

    CALL art_read_optprops_ext_ssa_asy_bsc(ncid, &
        &                  TRIM(mode_name)//'_aeronet', 'wave_aeronet', &
        &                  ext_coeff, unused_coeff, unused_coeff, unused_coeff)
  END SUBROUTINE art_read_aeronet_optprops

  SUBROUTINE art_read_satellite_optprops(ncid, mode_name, ext_coeff, bsc_coeff)
  !<
  ! SUBROUTINE art_read_prognostic_optprops
  ! Reads only extinction, SSA and backscatter coefficients for the prognostic radiation calculation of aerosols.
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-01-10
  !>
    !arguments
    INTEGER, INTENT(IN)               :: ncid
    CHARACTER(LEN=*), INTENT(IN)      :: mode_name
    REAL(wp), POINTER, INTENT(INOUT)  :: ext_coeff(:),  &
                                      &  bsc_coeff(:)
    !local arguments
    REAL(wp), POINTER :: unused_coeff(:)

    unused_coeff => NULL() !< mark this coefficient unused by pointing to null (de-association)

    CALL art_read_optprops_ext_ssa_asy_bsc(ncid, &
        &                  TRIM(mode_name)//'_satellite', 'wave_satellite', &
        &                  ext_coeff, unused_coeff, unused_coeff, bsc_coeff)
  END SUBROUTINE art_read_satellite_optprops

  SUBROUTINE art_read_singlewave_optprops(ncid, mode_name, ext_coeff)
  !<
  ! SUBROUTINE art_read_prognostic_optprops
  ! Reads only extinction coefficients for the diagnostic radiation calculation of aerosols at a single wavelength.
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-01-10
  !>
    !arguments
    INTEGER, INTENT(IN)               :: ncid
    CHARACTER(LEN=*), INTENT(IN)      :: mode_name
    REAL(wp), POINTER, INTENT(INOUT)  :: ext_coeff(:)
    !local arguments
    REAL(wp), POINTER :: unused_coeff(:)

    unused_coeff => NULL() !< mark this coefficient unused by pointing to null (de-association)

    CALL art_read_optprops_ext_ssa_asy_bsc(ncid, &
        &                  TRIM(mode_name)//'_singlewave', 'wave_single', &
        &                  ext_coeff, unused_coeff, unused_coeff, unused_coeff)
  END SUBROUTINE art_read_singlewave_optprops

  SUBROUTINE art_read_diagnostic_optprops(ncid, mode_name,              &
    &                                     ext_aeronet, ext_snglwave,    &
    &                                     ext_satellite, bsc_satellite)
  !<
  ! SUBROUTINE art_read_diagnostic_optprops
  ! Reads all optical properties needed for the radiation diagnostics
  ! Part of Module: mo_art_read_opt_props
  ! Author: Enrico P. Metzner, KIT
  ! Initial Release: 2025-02-10
  !>
    !arguments
    INTEGER, INTENT(IN)               :: ncid
    CHARACTER(LEN=*), INTENT(IN)      :: mode_name
    REAL(wp), POINTER, INTENT(INOUT)  :: ext_aeronet(:),   &
      &                                  ext_snglwave(:),  &
      &                                  ext_satellite(:), &
      &                                  bsc_satellite(:)
    !local variables
    REAL(wp), POINTER ::                 unused_coeff(:)
    
    unused_coeff => NULL() !< mark this coefficient unused by pointing to null (de-association)

    !read AERONET diagnostic opt. props.
    CALL art_read_optprops_ext_ssa_asy_bsc(ncid, &
        &                  TRIM(mode_name)//'_aeronet', 'wave_aeronet', &
        &                  ext_aeronet, unused_coeff, unused_coeff, unused_coeff)
    
    !read diagnostic opt. props. for a single wavelength (550nm)
    CALL art_read_optprops_ext_ssa_asy_bsc(ncid, &
        &                  TRIM(mode_name)//'_singlewave', 'wave_single', &
        &                  ext_snglwave, unused_coeff, unused_coeff, unused_coeff)
    
    !read satellite diagnostic opt. props.
    CALL art_read_optprops_ext_ssa_asy_bsc(ncid, &
        &                  TRIM(mode_name)//'_satellite', 'wave_satellite', &
        &                  ext_satellite, unused_coeff, unused_coeff, bsc_satellite)

  END SUBROUTINE art_read_diagnostic_optprops

  SUBROUTINE art_read_dia_factors(ncid, dia_min_factor, dia_max_factor)
!<
! SUBROUTINE art_read_dia_factors
! This subroutine reads metadata for optical parameters of aerosols
! Part of Module: mo_art_read_opt_props
! Author: Enrico P. Metzner, KIT
! Initial Release: 2025-02-12
! Modifications:
!>
    INTEGER, INTENT(IN) :: &
      &     ncid                          !< fileaccess-id for the netcdf database of optical properties
    REAL(wp), POINTER, INTENT(INOUT) :: &
      &     dia_min_factor(:),          &
      &     dia_max_factor(:)
    !local variables
    INTEGER                             :: ierror, varid
    REAL(wp), DIMENSION(:), ALLOCATABLE :: temp
    CHARACTER(LEN=*), PARAMETER         :: sroutine=TRIM(routine)//':art_read_dia_factors'

    ALLOCATE( temp(2) )

    ierror = nf90_inq_varid(ncid, name='dia_factors', varid=varid)
    IF (ierror == NF90_NOERR) THEN
      ierror = nf90_get_var(ncid, varid, temp(:), start=(/1/), count=(/2/))
    END IF
    IF (ierror /= NF90_NOERR) THEN
      DEALLOCATE( temp )
      CALL finish(TRIM(sroutine),'Could not read dia factors')
    END IF

    dia_min_factor(1) = temp(1)
    dia_max_factor(1) = temp(2)

    DEALLOCATE( temp )

  END SUBROUTINE art_read_dia_factors

  SUBROUTINE art_read_vardia_optprops(ncid, mode_name,                       &
    &                                 ext_params, ssa_params, asy_params,    &
    &                                 ext_default, ssa_default, asy_default)
!<
! SUBROUTINE art_read_vardia_optprops
! This subroutine reads metadata for optical parameters of aerosols
! Part of Module: mo_art_read_opt_props
! Author: Enrico P. Metzner, KIT
! Initial Release: 2025-02-12
! Modifications:
!>
    INTEGER, INTENT(IN)              :: ncid
    CHARACTER(LEN=*), INTENT(IN)     :: mode_name
    REAL(wp), POINTER, INTENT(INOUT) :: ext_params(:,:),  &
      &                                 ssa_params(:,:),  &
      &                                 asy_params(:,:),  &
      &                                 ext_default(:,:), &
      &                                 ssa_default(:,:), &
      &                                 asy_default(:,:)
    !local variables
    INTEGER                                 :: i, p, ierror, varid, npoly, num_waves
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: params, defaults
    CHARACTER(LEN=*), PARAMETER             :: sroutine=TRIM(routine)//':art_read_vardia_optprops'

    npoly = 4
    num_waves = 30
    
    ALLOCATE( params(npoly,num_waves,4) )
    ierror = nf90_inq_varid(ncid, TRIM(mode_name)//'_prog_params', varid=varid)
    IF (ierror == NF90_NOERR) THEN
      ierror = nf90_get_var(ncid, varid, params(:,:,:), &
        &                   start=(/1,1,1/), count=(/npoly,num_waves,4/))
    END IF
    IF (ierror /= NF90_NOERR) THEN
      DEALLOCATE( params )
      CALL finish(TRIM(sroutine), 'Could not read var dia coefficients for mode '// &
        &         TRIM(mode_name)//'\n'//nf90_strerror(ierror))
    END IF
    DO i=1,num_waves
      DO p=1,npoly
        ext_params(i,p) = params(p,i,1)
        ssa_params(i,p) = params(p,i,2)
        asy_params(i,p) = params(p,i,3)
      END DO
    END DO
    DEALLOCATE( params )

    ALLOCATE( defaults(2,num_waves,4) )
    ierror = nf90_inq_varid(ncid, TRIM(mode_name)//'_prog_defaults', varid=varid)
    IF (ierror == NF90_NOERR) THEN
      ierror = nf90_get_var(ncid, varid, defaults(:,:,:), &
        &                   start=(/1,1,1/), count=(/2,num_waves,4/))
    END IF
    IF (ierror /= NF90_NOERR) THEN
      DEALLOCATE( defaults )
      CALL finish(TRIM(sroutine), 'Could not read var dia default range for mode '// &
        &         TRIM(mode_name)//'\n'//nf90_strerror(ierror))
    END IF
    DO i=1,num_waves
      ext_default(i,1) = defaults(1,i,1)
      ext_default(i,2) = defaults(2,i,1)
      ssa_default(i,1) = defaults(1,i,2)
      ssa_default(i,2) = defaults(2,i,2)
      asy_default(i,1) = defaults(1,i,3)
      asy_default(i,2) = defaults(2,i,3)
    END DO
    DEALLOCATE( defaults )

  END SUBROUTINE art_read_vardia_optprops

END MODULE mo_art_read_opt_props
