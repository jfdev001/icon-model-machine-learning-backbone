!
! mo_art_gas_depo
! This module provides the calculation of the deposition
! velocities of gases
! Based on M. Baer - Parameterization of Trace Gas Dry Deposition 
! for a Regional Mesoscale Diffusion Model
! (1992)   
!
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

MODULE mo_art_gas_depo
#ifdef __ICON_ART
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_fortran_tools,                 ONLY: t_ptr_tracer,init
! ART
  USE mo_art_aerosol_utilities,         ONLY: calc_aerodynamic_resistance,calc_ustar
  USE mo_art_clipping,                  ONLY: art_clip_gt
  
IMPLICIT NONE
  
  PRIVATE
    
  PUBLIC :: art_calc_v_drydep, art_calc_loss_wetdepo

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_v_drydep(temp,temp_sfc,u,v,rho,rh,tcm,tch,qv,qc,qr,theta_v,gz0,dz, &
  &                       dyn_visc,istart,iend,tr_name,vdmol1,vdmol2,reac,heff,pabs,   &
  &                       luse_class,luse_fraction,vdep)
!<
! SUBROUTINE art_calc_v_drydep
! This subroutine calculates dry deposition velocities
! Based on: Baer (1992)
! Part of Module: mo_art_gas_depo
! Author: Sven Werchner, KIT
! Initial Release: 2019-07-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)   :: &
    &  temp(:),             & !< Temperature in lowest layer [K]
    &  temp_sfc(:),         & !< Surface temperature [K]
    &  u(:),v(:),           & !< Horizontal wind in lowest layer [m s-1]
    &  rho(:),              & !< Air density [kg m-3]
    &  rh(:),               & !< Relative humidity [%]
    &  tcm(:), tch(:),      & !< Transfer coefficients for momentum and heat in lowest layer [--]
    &  qv(:),qc(:),qr(:),   & !< Mass mixing ratios for water vapor / cloud water / rain water [kg kg-1]
    &  theta_v(:),          & !< Virtual potential temperature [K]
    &  gz0(:),              & !< roughness length * g of the vertically not resolved canopy [m2 s-2]
    &  dz(:),               & !< Layer height [m]
    &  dyn_visc(:),         & !< Dynamic viscosity
    &  vdmol1,              & !< ratio of mol. diffusion coeff. D_H2O/D_x 
    &  vdmol2,              & !< ratio of mol. diffusion coeff. D_i/D_x with D_i the molecular diffusivity
                              !< of the species for which rsmin was determined (usually H2O or CO2, with current 
                              !< parametrization D_i = D_H2O (-> vdmol1 = vdmol2 in this case)) (see Baer, Eq. 12)
    &  reac, heff,          & !<
    &  pabs(:),             & !< photosynthetic active radiation
    &  luse_fraction(:,:)     !< land use fraction in area
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end of inner loop (nproma)
    &  luse_class(:,:)        !< land use class in tile
  CHARACTER(*),INTENT(in):: &
    &  tr_name
  REAL(wp), INTENT(out)  :: &
    &  vdep(:)                !< Deposition velocities for 0th / 3rd moment
! Local variables
  REAL(wp)               :: &
    &  ustar,               & !< Friction velocity
    &  raerody_rb,rc(7),    & !< aerodynamical+boundary layer, canopy resistance (ncat=7)
    &  rcbase_inv(7),       & !< inverse values of basis for canopy resistance
    &  lstar_inv,           & !< 1./Monin-Obukhov
    &  luse_dep(7),         & !< deposition: land use fraction
    &  rst_inv(7),          & !< inverse values of stomata resistance
    &  rst_inv100(7),       & !< inverse values of stomata resistance for RH>=100%
    &  nu,                  & !< nu = dynamic viscosity / density
    &  rmes,                &
    &  rst,                 &
    &  rfak,                &
    &  a1,a2,               &
    &  ftemp,               &
    &  tdiff
  REAL(wp) ::              &
    & rclso2(7),           &
    & rcgso2(7),           &
    & rcgo3(7),            &
    & xlai(7),             &
    & rclo3(7),            &
    & tmin(7),             &
    & topt(7),             &
    & tmax(7),             &
    & rsmin(7),            &
    & bpar(7) 
  REAL(wp),PARAMETER    ::     &
    & kap_inv = 2.5_wp,        & !< 1./von Karman
    & rst0    = 1.e-15_wp,     &
    & f3      = 1._wp/3._wp,   &
    & f107    = 1.e-7_wp,      &
    & f1000   = 0.001_wp,      &
    & f2000   = 0.0005_wp,     &
    & f3000   = 1._wp/3000._wp,&
    & f5000   = 0.0002_wp
  INTEGER,PARAMETER ::     &
    & ncat = 7
  LOGICAL ::                &
    & l_so2, l_o3, l_h2so4
  INTEGER                :: &
    &  jc,i,n_tiles            !< Loop index

  
  CALL init(vdep, lacc=.FALSE.)
  
  CALL init(luse_dep, lacc=.FALSE.)

  l_so2   = .FALSE.
  l_o3    = .FALSE.
  l_h2so4 = .FALSE.
  
  SELECT CASE(TRIM(tr_name))
    CASE('TRSO2','SO2')
      l_so2   = .TRUE.
    CASE('TRO3', 'O3')
      l_o3    = .TRUE.
    CASE('TRH2SO4', 'H2SO4')
      l_h2so4 = .TRUE.
    CASE DEFAULT
      ! Nothing to do
  END SELECT
  
  raerody_rb = 0._wp
  lstar_inv  = 0._wp

  ! preparing
  rclso2 =       (/9999._wp,7000._wp,4000._wp,4000._wp,       &  
    &              4000._wp,4000._wp,9999._wp/)
  rcgso2 =       (/ 500._wp, 300._wp, 500._wp, 500._wp,       &  
    &               150._wp,  75._wp,  0.1_wp/)
  rcgo3  =       (/ 300._wp, 250._wp, 200._wp, 200._wp,       &  
    &               150._wp,1000._wp,2000._wp/)
  xlai   =       (/   0._wp,  0.3_wp,  3.5_wp,  3.5_wp,       &  
    &                 3._wp,   3._wp,   0._wp/)
  rclo3  =       (/9999._wp,3000._wp, 700._wp, 700._wp,       &  
    &               500._wp,1000._wp,9999._wp/)
  tmin   =       (/   5._wp,   5._wp,   0._wp,   0._wp,       &  
    &                 5._wp,   5._wp,   5._wp/)
  topt   =       (/  25._wp,  25._wp,  15._wp,  15._wp,       &  
    &                25._wp,  25._wp,  25._wp/)
  tmax   =       (/  45._wp,  45._wp,  40._wp,  40._wp,       &  
    &                45._wp,  45._wp,  45._wp/)
  rsmin  =       (/   1._wp, 350._wp, 350._wp, 350._wp,       &  
    &               150._wp, 120._wp,   1._wp/)
  bpar   =       (/   0._wp,  25._wp,  25._wp,  25._wp,       &  
    &                40._wp,  40._wp,   0._wp/)

  ! preparation
  rmes = 1._wp/(heff/3000._wp+100._wp*reac)
  
  ! rc_base^-1 = xlai/rcut + 1/rsoil  [independent of jc]
  ! only used, if not H2SO4
  IF(.NOT.l_h2so4) THEN
    DO i=1,ncat
      rcbase_inv(i) = xlai(i)/rclso2(i)*(heff/4.e+04_wp+reac) + MAX(heff/4.e+04_wp/rcgso2(i)+reac/rcgo3(i),1.e-4_wp)
      ! jc-independent, default stomata (used when RH>95%)
      IF (l_so2) THEN
        rst_inv(i)    = 1.e-2_wp
        rst_inv100(i) = 1._wp/(f5000+f3/rclso2(i))
      ELSE IF (l_o3) THEN
        rst_inv(i)= f3000 + f3 / rclo3(i)
        rst_inv100(i) = 1._wp/(f1000+f3/rclo3(i))
      ELSE
        rst_inv(i) = f3/rclso2(i)+f107*heff           &
          &         + reac/(1._wp/(f2000+f3/rclo3(i)))
      ENDIF
    ENDDO
  ENDIF

  ! jc-loop
  DO jc = istart, iend
    n_tiles = SIZE(luse_fraction(jc,:))
    CALL calc_luse_depo(luse_class(jc,:),luse_fraction(jc,:),n_tiles,luse_dep) 
    ! Eq. 3.71 precalculate several factors 
    CALL calc_ustar(ustar,tcm(jc),u(jc),v(jc))
    CALL calc_aerodynamic_resistance(ustar, theta_v(jc), gz0(jc), qv(jc), qc(jc),   &
      &                              qr(jc), temp(jc), temp_sfc(jc),dz(jc),tch(jc), &
      &                              tcm(jc),raerody_rb,lstar_inv)
    ! calculate (and add) surface boundary layer resistance
    raerody_rb = raerody_rb + 2._wp*kap_inv/ustar * (0.15_wp*4._wp*vdmol1)**0.667_wp

    ! calculate canopy resistance
    IF(l_h2so4) THEN ! special case for H2SO4
      rc = 1._wp/(0.002_wp*ustar)
      IF(lstar_inv < 0.0_wp) THEN
        rc = rc/(1._wp+(-300._wp*lstar_inv)**0.667_wp)
      ENDIF
    ELSE
      ! ... by calculating (inverse of) stomata resistance
      IF (rh(jc) >= 100._wp) THEN
        IF (l_so2.OR.l_o3) THEN
          rst_inv = rst_inv100
        ENDIF
      ELSE IF(rh(jc) <= 95._wp) THEN
        tdiff = temp(jc)-273.15_wp
        DO  i=1,ncat
          ftemp=1.e-5_wp
          IF(tdiff < tmax(i).AND.tdiff > tmin(i)) THEN
            a1 = tmax(i)-topt(i)
            a2 = topt(i)-tmin(i)
            ftemp = ( tdiff-tmin(i))/a2            &
             &    * ((tmax(i)-tdiff)/a1)**(a1/a2)
          END IF
          IF (pabs(jc) > 0.0_wp) THEN
            rfak = rsmin(i)*(1._wp+bpar(i)/pabs(jc))/ftemp
            rst = rfak*vdmol2+rmes
            rst = MIN(rst,10000._wp)
            rst_inv(i) = 1._wp/rst
          ELSE  ! use default value
            rst_inv(i) = rst0
          END IF
        END DO
      END IF
      ! ... and combining with rcbase_inv
      rc = 1._wp / (rst_inv*xlai + rcbase_inv)
    ENDIF
    DO i=1,ncat
      vdep(jc) = vdep(jc) + luse_dep(i) / (raerody_rb + rc(i))

    ENDDO
!   IF (TRIM(tr_name)=="DMSO") THEN    !HEIKE, haben wir DMSO?
!     vdep(jc)=1.e-02_wp*luse_dep(7)
!   ENDIF
      
!     vdep(jc) = 0._wp
  ENDDO !jc
  
  CALL art_clip_gt(vdep(:),0.01_wp) 
  
END SUBROUTINE art_calc_v_drydep

SUBROUTINE calc_luse_depo(luse_class,luse_fraction,n_tiles,luse_dep)
 INTEGER, INTENT(IN)       :: &
  & luse_class(:),            & !< landuse classes in tile (ntiles)
  & n_tiles
 REAL(wp), INTENT(IN)      :: &
  & luse_fraction(:)          !< landuse fractions (ntiles)
 REAL(wp),     INTENT(OUT)     ::                      &
  & luse_dep(:)

! internal variables
 REAL(wp)     :: &
  & luse(23)
 INTEGER      :: &
  & i

!-----------------------------------------------------------------------
! Transform GLOBCOVER2009 landuse classes to deposition landuse classes BWU
! (according to K. Nester (1996))

! landuse GLOBCOVER2009                              | landuse deposition
! luse [1]                                           | BWU [1]
! --------------------------------------------------------------------------
! 01 irrigated croplands                               01 urban area (dense)
! 02 rainfed croplands                                 02 rural area
! 03 mosaic cropland (50-70%) - vegetation (20-50%)    03 forest
! 04 mosaic vegetation (50-70%) - cropland (20-50%)    04 grassland
! 05 closed broadleaved evergreen forest               05 cropland dry
! 06 closed broadleaved deciduous forest               06 cropland wet
! 07 open broadleaved deciduous forest                 07 water bodies
! 08 closed needleleaved evergreen forest              08 water bodies
! 09 open needleleaved decid. or evergr. forest
! 10 mixed broadleaved and needleleaved forest
! 11 mosaic shrubland (50-70%) - grassland (20-50%)
! 12 mosaic grassland (50-70%) - shrubland (20-50%)
! 13 closed to open shrubland
! 14 closed to open herbaceous vegetation
! 15 sparse vegetation
! 16 closed to open forest regulary flooded
! 17 closed forest or shrubland permanently flooded
! 18 closed to open grassland regularly flooded
! 19 artificial surfaces
! 20 bare areas
! 21 water bodies
! 22 permanent snow and ice
! 23 undefined

! CALL init(luse_dep)
! CALL init(luse)
  luse = 0._wp

  DO i=1,n_tiles
    IF(luse_class(i)==-1) CYCLE
    ! care for lsnowtile case - where to put the snow
    ! current version: basically ignoring snow and adding both fractions
    ! together
    ! TODO: search for better treatment
    luse(luse_class(i)) = luse(luse_class(i)) + luse_fraction(i)
!print *, 'lanu',i,luse(luse_class(i)),luse_fraction(i),luse_class(i)
  ENDDO

! landuse class attribution GLOBCOVER2009 -> deposition, KD June 2015
! 01 urban area (dense)
  luse_dep(1) = luse(19) * 0.5_wp

! 02 rural area
  luse_dep(2) = luse(19) * 0.5_wp

! 03 forest
  luse_dep(3) = luse( 5) + luse( 6) + luse( 7) +  &
    &           luse( 8) + luse( 9) + luse(10) +  &
    &           luse(16) + luse(17)
! 04 grassland
  luse_dep(4) = luse(11) + luse(12) + luse(13) +  &
   &            luse(14) + luse(15) * 0.5_wp   +  &
   &            luse(18)

! 05 cropland dry
  luse_dep(5) = luse(1) * 0.5_wp + luse(2) * 0.5_wp + &
   &            luse(3) * 0.5_wp + luse(4) * 0.5_wp + &
   &            luse(20)

! 06 cropland wet
  luse_dep(6) = luse(1) * 0.5_wp + luse(2) * 0.5_wp + &
   &            luse(3) * 0.5_wp + luse(4) * 0.5_wp + &
   &            luse(15) * 0.5_wp
! 07 water bodies
  luse_dep(7) = luse(21) + luse(22)

END SUBROUTINE calc_luse_depo

SUBROUTINE art_calc_loss_wetdepo(tracer, qr, qc, dtime, rho,  &
                                    & istart, iend, nlev, &
                                    & rr_conv_3d,         &
                                    & rr_conv, rr_gsp,    &
                                    & w_in, w_sub,        &
                                    & loss_tracer)
!<
! SUBROUTINE art_calc_loss_wetdepo
! This subroutine calculates loss of tracer due to wet deposition
! Based on: Simpson et al. The EMEP MSC-W chemical transport model
!   - technical description (2012)
! Part of Module: mo_art_gas_depo
!>

! ---- in/out ---------------------------------------------------------

REAL(wp), INTENT(IN) ::                     &
  &                 qr(:,:),                & !< specific rain content [kg/kg]
  &                 qc(:,:),                & !< specific cloud water content [kg/kg]
  &                 dtime,                  & !< time step
  &                 rho(:,:),               & !< density of air 
  &                 rr_conv_3d(:,:),        & !< grid-scale rain rate [kg/(m2 s)]
  &                 rr_gsp(:),              & !< grid-scale rain rate at surface [kg/(m2 s)]
  &                 rr_conv(:),             & !< convective rain rate at surface [kg/(m2 s)]
  &                 w_in,                   & !< in-cloud scavenging ratio
  &                 w_sub                     !< sub-cloud scavenging ratio

REAL(wp),INTENT(INOUT) ::                   &
  &                 tracer(:,:)             !< mmr of species

INTEGER, INTENT(IN) ::                      &
  &                 istart, iend,           & !< start and end of block
  &                 nlev                      !< number of vertical levels

REAL(wp), INTENT(INOUT) ::                        &
  &                 loss_tracer(:,:)

! ---- local ----------------------------------------------------------
INTEGER             :: &
  &  jc,               & !< counter for nproma loop
  &  jk                  !< counter for vertical loop

REAL(wp)            :: &
  &  lwp_kgm3,         & !< liquid water content from precipitation [kg/m3]
  &  precip              !< rain rate (gsp + conv) at surface [kg/(m2 s)]

REAL(wp),PARAMETER    ::      &
  &  hs = 1000.0_wp,          & !< characteristic scavenging depth [m]
  &  rho_w = 1000.0_wp,       & !< density water [kg/m3]
  &  v_rain = 10.0_wp           !< vertical velocity of rain droplets [m/s]

REAL(wp) :: &
  &  tracer_in

! ---- begin -----------------------------------------------------------

DO jk=1,nlev
  DO jc=istart, iend

    tracer_in = tracer(jc,jk)

    ! local liquid water content from precipitation
    lwp_kgm3 = qr(jc,jk) * rho(jc,jk) + (rr_conv_3d(jc,jk) / v_rain)
    
    ! If liquid water is present, do washout
    IF (lwp_kgm3 > 1.0E-15_wp) THEN

      ! rain rate at surface
      precip = rr_gsp(jc) + rr_conv(jc)

      ! If cloud is present, do in-cloud washout, else below cloud
      IF (qc(jc,jk) > 1.0E-15_wp) THEN 

        ! Equation (71) in Simpson et al, ACP, 2012
        tracer(jc,jk) = tracer(jc,jk) * EXP(-dtime*w_in*precip/(hs*rho_w)) 

      ELSE ! below cloud

        ! Equation (72) in Simpson et al, ACP, 2012
        tracer(jc,jk) = tracer(jc,jk) * EXP(-dtime*w_sub*precip/(hs*rho_w))

      ENDIF 

    ENDIF ! precip present
    loss_tracer(jc,jk) =  tracer_in - tracer(jc,jk)
  ENDDO ! jc
ENDDO ! jk



END SUBROUTINE art_calc_loss_wetdepo

#endif
END MODULE mo_art_gas_depo
