!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          SWRAD: MM5 SW (Dudhia) scheme from WRF
!!
!! @author Team SCALE
!!
!! @par History
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_rd_mm5sw
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_atmos_phy_rd_common, only: &
     I_SW, &
     I_LW, &
     I_dn, &
     I_up

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SWRAD
  public :: swinit

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  REAL(RP)  :: CSSCA

  CONTAINS

!------------------------------------------------------------------
   SUBROUTINE SWRAD(dt,RTHRATEN,SDOWN3D,GSW,XLAT,XLONG,ALBEDO,    &
                    rho_phy,T3D,P3D,pi3D,dz8w,                    &
                    solins,cosSZA,                                &
                    QV3D,QC3D,QR3D,QI3D,QS3D,QG3D,                &
                    F_QV,F_QC,F_QR,F_QI,F_QS,F_QG,                &
                    icloud,warm_rain                              )

    use scale_const, only:  &
       GRAV  => CONST_GRAV, & ! 9.81
       Rdry  => CONST_Rdry, & ! 287.
       CPdry => CONST_CPdry   ! 1004.6
    use scale_atmos_solarins, only: &
       SOLARINS => ATMOS_SOLARINS_constant ! 1360.xxx (W/m2)
!------------------------------------------------------------------
   IMPLICIT NONE
!------------------------------------------------------------------

   real(DP), INTENT(IN)      :: dt
   real(RP), INTENT(INOUT)   :: RTHRATEN(IA,KA,JA)
   real(RP), INTENT(INOUT)   :: SDOWN3D(IA,KA,JA)
   real(RP), INTENT(INOUT)   :: GSW(IA,JA)

   real(RP), INTENT(IN)      :: XLAT(IA,JA), XLONG(IA,JA), ALBEDO(IA, JA)
   real(RP), INTENT(IN)      :: rho_phy(IA,KA,JA)
   real(RP), INTENT(IN)      :: P3D(IA,KA,JA),  &
                                T3D(IA,KA,JA),  &
                                pi3D(IA,KA,JA), &
                                dz8w(IA,KA,JA)
   real(RP), INTENT(IN)      :: solins(IA,JA),cosSZA(IA,JA)
   real(RP), OPTIONAL, INTENT(IN) :: QV3D (IA,KA,JA), &
                                     QC3D (IA,KA,JA), &
                                     QR3D (IA,KA,JA), &
                                     QI3D (IA,KA,JA), &
                                     QS3D (IA,KA,JA), &
                                     QG3D (IA,KA,JA)

   LOGICAL, OPTIONAL, INTENT(IN ) :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG
   INTEGER,  INTENT(IN   )        :: icloud
   LOGICAL,  INTENT(IN   )        :: warm_rain

   !real, INTENT(IN    )          :: RADFRQ,DEGRAD,XTIME,DECLIN


   integer   :: its,ite,jts,jte,kts,kte
   real(RP)  ::  R, CP, G, SOLCON

! LOCAL VARS

   real(RP), DIMENSION( KS:KE ) :: TTEN1D, &
                                   RHO01D, &
                                   P1D, &
                                   DZ, &
                                   T1D, &
                                   QV1D, &
                                   QC1D, &
                                   QR1D, &
                                   QI1D, &
                                   QS1D, &
                                   QG1D

   real(RP)                     :: SDOWN1D(KS:KE+1)
!
   real(RP) :: XLAT0,XLONG0,ALB0,GSW0,cosSZA0,solins0
   real(RP) :: aer_dry1(KS:KE),aer_water1(KS:KE)

!
   INTEGER :: i,j,K,NK
   LOGICAL :: predicate , do_topo_shading


!------------------------------------------------------------------

   R  = Rdry
   CP = CPdry
   G  = GRAV
   SOLCON = SOLARINS
   its = IS ; ite = IE
   jts = JS ; jte = JE
   kts = KS ; kte = KE


   j_loop: DO J=jts,jte
   i_loop: DO I=its,ite

! reverse vars
         DO K=KS,KE
            QV1D(K)=0.
            QC1D(K)=0.
            QR1D(K)=0.
            QI1D(K)=0.
            QS1D(K)=0.
            QG1D(K)=0.
         ENDDO

         DO k=KS,KE
            !! NK=kme-1-K+kms
            NK=KE-(k-KS)
            TTEN1D(K)=0.

            T1D(K)=T3D(I,NK,J)
            P1D(K)=P3D(I,NK,J)
            RHO01D(K)=rho_phy(I,NK,J)
            DZ(K)=dz8w(I,NK,J)
         ENDDO

         !IF( PRESENT(pm2_5_dry) .AND. PRESENT(pm2_5_water) )THEN
         !   DO K=kts,kte
         !      NK=kme-1-K+kms
         !      aer_dry1(k)   = pm2_5_dry(i,nk,j)
         !      aer_water1(k) = pm2_5_water(i,nk,j)
         !   ENDDO
         !ELSE
             do k=KS,KE
               aer_dry1(k)   = 0.0
               aer_water1(k) = 0.0
             enddo
         !ENDIF

         IF (PRESENT(F_QV) .AND. PRESENT(QV3D)) THEN
            IF (F_QV) THEN
               do k=KS,KE
                  !NK=kme-1-K+kms
                  NK=KE-(k-KS)
                  QV1D(K)=QV3D(I,NK,J)
                  QV1D(K)=max(0.0_RP,QV1D(K))
               ENDDO
            ENDIF
         ENDIF

         IF (PRESENT(F_QC) .AND. PRESENT(QC3D)) THEN
            IF (F_QC) THEN
               do k=KS,KE
                  !NK=kme-1-K+kms
                  NK=KE-(k-KS)
                  QC1D(K)=QC3D(I,NK,J)
                  QC1D(K)=max(0.0_RP,QC1D(K))
               enddo
            ENDIF
         ENDIF

         IF (PRESENT(F_QR) .AND. PRESENT(QR3D)) THEN
            IF (F_QR) THEN
               do k=kts,kte
                  !NK=kme-1-K+kms
                  NK=KE-(k-KS)
                  QR1D(K)=QR3D(I,NK,J)
                  QR1D(K)=max(0.0_RP,QR1D(K))
               enddo
            ENDIF
         ENDIF

!
         IF ( PRESENT( F_QI ) ) THEN
            predicate = F_QI
         ELSE
            predicate = .FALSE.
         ENDIF

         IF ( predicate .AND. PRESENT( QI3D ) ) THEN
            do k=KS,KE
               !NK=kme-1-K+kms
                NK=KE-(k-KS)
                QI1D(K)=QI3D(I,NK,J)
                QI1D(K)=max(0.0_RP,QI1D(K))
            enddo
         ELSE
            IF (.not. warm_rain) THEN
               do k=KS,KE
               IF(T1D(K) .lt. 273.15) THEN
                  QI1D(K)=QC1D(K)
                  QC1D(K)=0.
                  QS1D(K)=QR1D(K)
                  QR1D(K)=0.
               ENDIF
               enddo
            ENDIF
         ENDIF

         IF (PRESENT(F_QS) .AND. PRESENT(QS3D)) THEN
            IF (F_QS) THEN
               do k=KS,KE
                 !NK=kme-1-K+kms
                  NK=KE-(k-KS)
                  QS1D(K)=QS3D(I,NK,J)
                  QS1D(K)=max(0.0_RP,QS1D(K))
               enddo
            ENDIF
         ENDIF

         IF (PRESENT(F_QG) .AND. PRESENT(QG3D)) THEN
            IF (F_QG) THEN
               do k=KS,KE
                 !NK=kme-1-K+kms
                  NK=KE-(k-KS)
                  QG1D(K)=QG3D(I,NK,J)
                  QG1D(K)=max(0.0_RP,QG1D(K))
               enddo
            ENDIF
         ENDIF

         XLAT0=XLAT(I,J)
         XLONG0=XLONG(I,J)
         solins0=solins(I,J)
         cosSZA0=cosSZA(I,J)
         ALB0=ALBEDO(I,J)
! slope code removed - factor now done in surface driver
           CALL SWPARA(TTEN1D,SDOWN1D,GSW0,ALB0,cosSZA0,                &
                       T1D,QV1D,QC1D,QR1D,QI1D,QS1D,QG1D,P1D,   &
                       RHO01D,DZ,                               &
                       R,CP,G,solins0,                          &
                       XLAT0,XLONG0,                            &
                       icloud,aer_dry1,aer_water1,              &
                       kts,kte      )
         GSW(I,J)=GSW0

!         DO K=kts,kte
!            NK=kme-1-K+kms
!            NK=KA-1-K+1
!            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+TTEN1D(NK)/pi3D(I,K,J)
!            SDOWN3D(I,K,J)=SDOWN1D(NK)
!            if(k==kte)then
!               NK=KA-1-(K+1)+1
!	       SDOWN3D(I,K+1,J)=SDOWN1D(NK)
!            endif
!         ENDDO

         do k=KS,KE
            NK=KE-(k-KS)
            RTHRATEN(I,K,J)=RTHRATEN(I,K,J)+TTEN1D(NK)/pi3D(I,K,J)
            SDOWN3D(I,K,J)=SDOWN1D(NK)
            !print *,KS,KE,K,"<=",NK
         enddo
            k=KS-1 ; NK=KE+1
            SDOWN3D(I,KS-1,J)=SDOWN1D(KE+1)
            !print *,KS,KE,K,"<=",NK
!
   ENDDO i_loop
   ENDDO j_loop

   return
   END SUBROUTINE SWRAD

!------------------------------------------------------------------
   SUBROUTINE SWPARA(TTEN,SDOWN,GSW,ALBEDO,cosSZA,  &
                     T,QV,QC,QR,QI,QS,QG,P,         &
                     RHO0, DZ,             	    &
                     R,CP,G,solins,                 &
                     XXLAT,XXLON,                  &
                     ICLOUD,aer_dry1,aer_water1,    &
                     kts,kte                        )

    use scale_time, only:       &
       NOWDATE => TIME_NOWDATE    !< current time [YYYY MM DD HH MM SS]
    use scale_const, only: &
       D2R    => CONST_D2R         ! degree to radian
    use scale_calendar, only: &
         CALENDAR_getDayOfYear,  &
         CALENDAR_ymd2absday,    &
         CALENDAR_hms2abssec,    &
         I_year, I_month, I_day, &
         I_hour, I_min, I_sec

!------------------------------------------------------------------
!     TO CALCULATE SHORT-WAVE ABSORPTION AND SCATTERING IN CLEAR
!     AIR AND REFLECTION AND ABSORPTION IN CLOUD LAYERS (STEPHENS,
!     1984)
!     CHANGES:
!       REDUCE EFFECTS OF ICE CLOUDS AND PRECIP ON LIQUID WATER PATH
!       ADD EFFECT OF GRAUPEL
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN ) ::                 kts,kte

  real(RP), DIMENSION( kts:kte ), INTENT(IN   )  ::     &
                                              RHO0, &
                                              T,  &
                                              P,  &
                                              DZ, &
                                              QV, &
                                              QC, &
                                              QR, &
                                              QI, &
                                              QS, &
                                              QG

   real(RP), DIMENSION( kts:kte ), INTENT(INOUT) ::  TTEN
   real(RP), DIMENSION( kts:kte+1 ), INTENT(OUT) ::  SDOWN
!
   real(RP), INTENT(IN  )   ::  R,CP,G
   real(RP), INTENT(IN)     ::  ALBEDO,cosSZA,solins
   integer , INTENT(IN)     ::  icloud
   real(RP), INTENT(INOUT)  ::  GSW
   real(RP), INTENT(IN  )   ::  XXLAT,XXLON

   !real, INTENT(IN  )      :: GMT,RADFRQ,XTIME,XLAT,XLONG,DEGRAD
!
! For slope-dependent radiation
!   INTEGER, OPTIONAL, INTENT(IN) :: slope_rad,shadow
!   real, OPTIONAL,    INTENT(IN) :: slp_azi,slope
!

! LOCAL VARS

   real(RP)  :: XLWP( kts:kte )
   real(RP)  :: XATP( kts:kte )
   real(RP)  :: XWVP( kts:kte )
   real(RP)  :: aer_dry1( kts:kte )
   real(RP)  :: aer_water1( kts:kte )
   real(RP)  :: RO( kts:kte )
!
   real(RP)  ::  ALBTAB(4,5), ABSTAB(4,5)
   real(RP)  ::  XMUVAL(4)

   real(RP) :: beta
   real(RP) :: DayOfYear
!------------------------------------------------------------------

      DATA ALBTAB/0.0,0.0,0.0,0.0, &
           69.0,58.0,40.0,15.0,    &
           90.0,80.0,70.0,60.0,    &
           94.0,90.0,82.0,78.0,    &
           96.0,92.0,85.0,80.0/

      DATA ABSTAB/0.0,0.0,0.0,0.0, &
           0.0,2.5,4.0,5.0,       &
           0.0,2.6,7.0,10.0,      &
           0.0,3.3,10.0,14.0,     &
           0.0,3.7,10.0,15.0/

      DATA XMUVAL/0.0,0.2,0.5,1.0/

      real(RP) :: bext340, absc, alba, alw, csza,dabsa,dsca,dabs
      real(RP) :: bexth2o, dscld, ff,oldalb,oldabs,oldabc
      real(RP) :: soltop, totabs, ugcm, uv,xabs,xabsa,wv
      real(RP) :: wgm, xalb, xi, xsca, xmu,xabsc,trans0,yj
      real(RP) :: HRANG,XT24,TLOCTM,tloc,dsec,lon,lat
      real(RP) :: DECLIN
      real(RP) :: ww
      INTEGER :: iil,ii,jjl,ju,k,iu

! For slope-dependent radiation

      real(RP) :: diffuse_frac, corr_fac, csza_slp

      GSW=0.0
      SDOWN=0.0
      bext340=5.E-6_RP
      bexth2o=5.E-6_RP
      SOLTOP=solins

      !XT24=MOD(XTIME+RADFRQ*0.5,1440.)
      !TLOCTM=GMT+XT24/60.+XLONG/15.
      !HRANG=15.*(TLOCTM-12.)*DEGRAD
      !XXLAT=XLAT*DEGRAD
      !CSZA=SIN(XXLAT)*SIN(DECLIN)+COS(XXLAT)*COS(DECLIN)*COS(HRANG)

      ! local time
      LAT = XXLAT / D2R
      LON = XXLON / D2R

      tloc   = mod((NOWDATE(4)    + int(LON/15.0_RP)),24 )
      dsec   = real(NOWDATE(5)*60 + NOWDATE(6)) / 60.0_RP /60.0_RP
      TLOCTM = real(NOWDATE(4))   + LON/15.0_RP + dsec
      HRANG  = 15.*(TLOCTM-12.)*D2R
      !
      !call CALENDAR_getDayOfYear( DayOfYear, now_date(I_year) )
      !IF( DayOfYear.GE.80.)SXLONG=DPD*( DayOfYear-80.)
      !IF( DayOfYear.LT.80.)SXLONG=DPD*( DayOfYear+285.)
      !ARG=SIN(23.5*D2R)*SIN(SXLONG*D2R)
      !DECLIN=ASIN(ARG)

      CSZA=cosSZA

!     RETURN IF NIGHT
      IF(CSZA .LE. 1.0E-9_RP)GOTO 7
!
      DO K=kts, kte

! P in the unit of 10mb
         RO(K)=P(K)/(R*T(K))
         XWVP(K)=RO(K)*QV(K)*DZ(K)*1000.
! KG/M**2
         XATP(K)=RO(K)*DZ(K)
      ENDDO
!
!     G/M**2
!     REDUCE WEIGHT OF LIQUID AND ICE IN SHORT-WAVE SCHEME
!     ADD GRAUPEL EFFECT (ASSUMED SAME AS RAIN)
!
      IF (ICLOUD.EQ.0)THEN
         DO K=kts, kte
            XLWP(K)=0.
         ENDDO
      ELSE
         DO K=kts, kte
            XLWP(K)=RO(K)*1000.*DZ(K)*(QC(K)+0.1*QI(K)+0.05* &
                    QR(K)+0.02*QS(K)+0.05*QG(K))
         ENDDO
      ENDIF
!
      XMU=CSZA
      ! SDOWN(1)=SOLTOP*XMU !adachi
      !SDOWN(kts)=SOLTOP*XMU
      SDOWN(kts)=solins
!     SET WW (G/M**2) LIQUID WATER PATH INTEGRATED DOWN
!     SET UV (G/M**2) WATER VAPOR PATH INTEGRATED DOWN
      WW=0.
      UV=0.
      OLDALB=0.
      OLDABC=0.
      TOTABS=0.
!     CONTRIBUTIONS DUE TO CLEAR AIR AND CLOUD
      DSCA=0.
      DABS=0.
      DSCLD=0.
!
! CONTRIBUTION DUE TO AEROSOLS (FOR CHEMISTRY)
      DABSA=0.
!
      DO 200 K=kts,kte
         WW=WW+XLWP(K)
         UV=UV+XWVP(K)
!     WGM IS WW/COS(THETA) (G/M**2)
!     UGCM IS UV/COS(THETA) (G/CM**2)
         WGM=WW/XMU
         UGCM=UV*0.0001/XMU
!
         OLDABS=TOTABS
!     WATER VAPOR ABSORPTION AS IN LACIS AND HANSEN (1974)
         TOTABS=2.9*UGCM/((1.+141.5*UGCM)**0.635+5.925*UGCM)
!     APPROXIMATE RAYLEIGH + AEROSOL SCATTERING
!        XSCA=1.E-5*XATP(K)/XMU
!          XSCA=(1.E-5*XATP(K)+aer_dry1(K)*bext340+aer_water1(K)*bexth2o)/XMU
         beta=0.4*(1.0-XMU)+0.1
!     CSSCA - CLEAR-SKY SCATTERING SET FROM NAMELIST SWRAD_SCAT
         XSCA=(cssca*XATP(K)+beta*aer_dry1(K)*bext340*DZ(K) &
              +beta*aer_water1(K)*bexth2o*DZ(K))/XMU

!     LAYER VAPOR ABSORPTION DONE FIRST
      !   XABS=(TOTABS-OLDABS)*(SDOWN(1)-DSCLD-DSCA-DABSA)/SDOWN(K)
          XABS=(TOTABS-OLDABS)*(SDOWN(kts)-DSCLD-DSCA-DABSA)/SDOWN(K)
!rs   AEROSOL ABSORB (would be elemental carbon). So far XABSA = 0.
         XABSA=0.
         IF(XABS.LT.0.)XABS=0.
!
         ALW=log10(WGM+1.0_RP)
         IF(ALW.GT.3.999)ALW=3.999_RP
!
         DO II=1,3
            IF(XMU.GT.XMUVAL(II))THEN
              IIL=II
              IU=II+1
              XI=(XMU-XMUVAL(II))/(XMUVAL(II+1)-XMUVAL(II))+FLOAT(IIL)
            ENDIF
         ENDDO
!
         JJL = int(ALW)+1
         JU  = JJL+1
         YJ  = ALW+1.
!     CLOUD ALBEDO
         ALBA=(ALBTAB(IU,JU)*(XI-IIL)*(YJ-JJL)   &
              +ALBTAB(IIL,JU)*(IU-XI)*(YJ-JJL)   &
              +ALBTAB(IU,JJL)*(XI-IIL)*(JU-YJ)   &
              +ALBTAB(IIL,JJL)*(IU-XI)*(JU-YJ))  &
             /((IU-IIL)*(JU-JJL))
!     CLOUD ABSORPTION
         ABSC=(ABSTAB(IU,JU)*(XI-IIL)*(YJ-JJL)   &
              +ABSTAB(IIL,JU)*(IU-XI)*(YJ-JJL)   &
              +ABSTAB(IU,JJL)*(XI-IIL)*(JU-YJ)   &
              +ABSTAB(IIL,JJL)*(IU-XI)*(JU-YJ))  &
             /((IU-IIL)*(JU-JJL))
!     LAYER ALBEDO AND ABSORPTION
         !XALB=(ALBA-OLDALB)*(SDOWN(1)-DSCA-DABS)/SDOWN(K)
         !XABSC=(ABSC-OLDABC)*(SDOWN(1)-DSCA-DABS)/SDOWN(K)
          XALB=(ALBA-OLDALB)*(SDOWN(kts)-DSCA-DABS)/SDOWN(K)
          XABSC=(ABSC-OLDABC)*(SDOWN(kts)-DSCA-DABS)/SDOWN(K)
         IF(XALB.LT.0.)XALB=0.
         IF(XABSC.LT.0.)XABSC=0.
         DSCLD=DSCLD+(XALB+XABSC)*SDOWN(K)*0.01
         DSCA=DSCA+XSCA*SDOWN(K)
         DABS=DABS+XABS*SDOWN(K)
         DABSA=DABSA+XABSA*SDOWN(K)
         OLDALB=ALBA
         OLDABC=ABSC
!     LAYER TRANSMISSIVITY
         TRANS0=100.0_RP-XALB-XABSC-XABS*100.0_RP-XSCA*100.0_RP
         IF(TRANS0.LT.1.)THEN
           FF=99.0_RP/(XALB+XABSC+XABS*100.0_RP+XSCA*100.0_RP)
           XALB=XALB*FF
           XABSC=XABSC*FF
           XABS=XABS*FF
           XSCA=XSCA*FF
           TRANS0=1.0_RP
         ENDIF
         SDOWN(K+1)=max(1.0E-9_RP,SDOWN(K)*TRANS0*0.01_RP)
         TTEN(K)=SDOWN(K)*(XABSC+XABS*100._RP+XABSA*100._RP)*0.01_RP/(RO(K)*CP*DZ(K))
  200   CONTINUE
!
        GSW=(1.-ALBEDO)*SDOWN(kte+1)

!    IF (PRESENT(slope_rad)) THEN
!! Slope-dependent solar radiation part
!
!      if (slope_rad.eq.1) then
!
!!  Parameterize diffuse fraction of global solar radiation as a function of the ratio between TOA radiation and surface global radiation
!
!        diffuse_frac = min(1.,1/(max(0.1,2.1-2.8*log(log(SDOWN(kts)/max(SDOWN(kte+1),1.e-3))))))
!        if ((slope.eq.0).or.(diffuse_frac.eq.1).or.(csza.lt.1.e-2)) then  ! no topographic effects when all radiation is diffuse or the sun is too close to the horizon
!        corr_fac = 1
!        goto 140
!        endif
!
!! cosine of zenith angle over sloping topography
!
!        csza_slp = ((SIN(XXLAT)*COS(HRANG))*                                          &
!                    (-cos(slp_azi)*sin(slope))-SIN(HRANG)*(sin(slp_azi)*sin(slope))+  &
!                    (COS(XXLAT)*COS(HRANG))*cos(slope))*                              &
!                   COS(DECLIN)+(COS(XXLAT)*(cos(slp_azi)*sin(slope))+                 &
!                   SIN(XXLAT)*cos(slope))*SIN(DECLIN)
!
!        csza_slp = 0
!        IF(csza_slp.LE.1.E-4) csza_slp = 0
!
! Topographic shading
!
!        if (shadow.eq.1) csza_slp = 0
!
! Correction factor for sloping topography; the diffuse fraction of solar radiation is assumed to be unaffected by the slope
!        corr_fac = diffuse_frac + (1-diffuse_frac)*csza_slp/csza
!
! 140	continue
!
!        GSW=(1.-ALBEDO)*SDOWN(kte+1)*corr_fac
!
!      endif
!    ENDIF

      !print *,"A1",TTEN(KE),SDOWN(KE),GSW,ALBEDO,cosSZA
      !print *,"A2",QV(KE),QC(KE),QR(KE),QI(KE),QS(KE),QG(KE)
      !print *,"A3",T(KE),P(KE),RHO0(KE),DZ(KE)
      !print *,"A4",R,CP,G,solins
      !print *,"A5",XXLAT,XXLON,XXLAT/D2R,XXLON/D2R
      !print *,"A6",ICLOUD,aer_dry1(KE),aer_water1(KE)
      !print *,"A7",kts,kte

    7 CONTINUE
      return
   END SUBROUTINE SWPARA

!====================================================================
   SUBROUTINE swinit
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
!   LOGICAL , INTENT(IN)           :: allowed_to_read
!   INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,  &
!                                     ims, ime, jms, jme, kms, kme,  &
!                                     its, ite, jts, jte, kts, kte
!   real , INTENT(IN)              :: swrad_scat

   LOGICAL            :: allowed_to_read = .true.
   real               :: swrad_scat = 1


!     CSSCA - CLEAR-SKY SCATTERING SET FROM NAMELIST SWRAD_SCAT
   cssca = swrad_scat * 1.e-5

   END SUBROUTINE swinit

end module scale_atmos_phy_rd_mm5sw
