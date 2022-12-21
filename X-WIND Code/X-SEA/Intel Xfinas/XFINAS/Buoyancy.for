
C     =========================================================================         
      SUBROUTINE BUOYANCYFORCE (NELEMENT,WATERDEPTH,SEABED,RHOW,ARATIO,BUOYANCYF)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      DIMENSION XBAR(NELE),YBAR(NELE),ZBAR(NELE),WEIGHT(NELE)
      DIMENSION ATMIDX(100),ATMIDY(100),ATMIDZ(100),DIAMX1(100),DIAMY1(100),DIAMZ1(100)
      DIMENSION DIAMX2(100),DIAMY2(100),DIAMZ2(100),ATEPPERL(100)
      DIMENSION XTBAR(100),YTBAR(100),ZTBAR(100),WTEIGHT(100),WEIGHTF(100)
      DIMENSION BPG(12),BWG(12),GPL(12),GPW(12)
      DIMENSION AREAT(100),TOTALWEIGHTT(NELE),XTG(NELE),YTG(NELE),ZTG(NELE)
      DIMENSION TOTALWEIGHTF(NELE)
      
      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
      
      ! OUTPUT FOR CALCULATION
     	COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      
      OFFSELECT = 1D0
      IF (KSPEC.NE.1) THEN
      CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,VB,HM,HW,HC,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                  VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,
     1                  WKF,CBF,WFC,TIME)
      ELSEIF (KSPEC.EQ.1)THEN
      CALL SELECTSPECTRUM (OFFSELECT,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
      ENDIF
      
      ! WATER DENSITY 
      !RHOW = 1025D0 > INPUT FROM OFFSHORE PARAMETER
      IF (SEABED.LT.0.0D0)THEN     ! 
      WATERDEPTH = WATERDEPTH + SEABED
      ELSEIF (SEABED.EQ.0.0D0)THEN ! 
      WATERDEPTH = WATERDEPTH
      ELSEIF (SEABED.GT.0.0D0)THEN !
      WATERDEPTH = WATERDEPTH - SEABED
      ENDIF
      
!      ! XG = CENTER OF GRAVITY, X POSITION
!      ! YG = CENTER OF GRAVITY, Y POSITION
!      ! ZG = CENTER OF GRAVITY, Z POSITION
!      CALL COG (XG,YG,ZG,ID,MSF)
!      ! XB = CENTER OF BUOYANCY, X POSITION
!      ! YB = CENTER OF BUOYANCY, Y POSITION
!      ! ZB = CENTER OF BUOYANCY, Z POSITION
!      CALL COB (WATERDEPTH,SEABED,RHOW,VGV,XB,YB,ZB)
      
C     -----------------------------------------------------------------------------------------------
      !  --- INPUT  ---
      !  NELEMENT  = NUMBER OF ELEMENT
      !  NOP       = 0 > ANALYSIS  BUOYANCY FORCE
      !  --- OUTPUT ---
      !  ALENGTH   = ELEMENT LENGTH
      NOP = 0
      CALL CALLENGTH (NELEMENT,VGV,0,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                ,ATMIDX,ATMIDY,ATMIDZ)
C     ------------------------------------------------------------------------------------------------

      
C     ------------------------------------------------------------------------------------------------
      !  --- INPUT  ---
      !  NELEMENT  = NUMBER OF ELEMENT
      !  --- OUTPUT ---
      !  DENSITY   = MASS DENSITY     
      CALL SELECTPROP (NELEMENT,AMODULUS,POSSION,DENSITY,NSTANDARD,FU,FY)
C     ------------------------------------------------------------------------------------------------
      
C     ------------------------------------------------------------------------------------------------  
      !  --- INPUT  ---
      !  NELEMENT  = NUMBER OF ELEMENT
      !  --- OUTPUT ---
      !  DOUT      = OUTER RIDUS A
      !  DIN       = INNER RIDUS A
      !  DOUT1     = OUTER RIDUS B
      !  DIN1      = INNER RIDUS B
      !
      !         ---           <<< ( B )
      !        l   l
      !       l     l              **********************
      !      l       l               CIRCULAR SECTION
      !     l         l            **********************
      !    l           l
      !   l             l 
      !   ---------------     <<< ( A )
      !
           
      CALL MAXSECTION (NELEMENT,DOUT,DIN,DOUT1,DIN1) ! OUTPUT IS RIDUS
      
C     ------------------------------------------------------------------------------------------------
      IF (DOUT.NE.DOUT1)THEN ! TAPPER SECTION ( IN CASE CIRCULAR SECTION )

           DO I = 1,100
           BL          = ALENGTH/100D0
           ATEPPERL(I) = ATEPPERL(I-1)+BL
           A           = (DOUT-DOUT1)/ALENGTH
           B           = (DIN-DIN1)/ALENGTH
           DIAMY1(I)   = DOUT-A*ATEPPERL(I)
           DIAMY2(I)   = DIN-B*ATEPPERL(I)
           
           AREA1       = 3.141592654D0*(((DIAMY1(I)*2D0)**2)-((DIAMY2(I)*2D0)**2))/4D0
           AREAT(I)    = ABS(AREA1)
           
           WTEIGHT(I)  = AREAT(I)*RHOW*BL*ARATIO
           ENDDO
           
           TOTALWEIGHTT(NELEMENT) = SUM(WTEIGHT)
           BUOYANCYF              = TOTALWEIGHTT(NELEMENT)
      
      
      ELSEIF (DOUT.EQ.DOUT1)THEN ! NORMAL SECTION
      ! ARATIO = MEMBER FLOOD RATIO
      DIAM3     = DOUT
      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      !AREA      = 3.141592654D0*((DOUT*2D0)**2D0-(DIN*2D0)**2D0)/4D0
      !IF (NSECTION.EQ.1)
      !IF (NSECTION.EQ.2)
      IF (NSECTION.EQ.21) AREA = 2D0*3.141592654D0*(DPIPE/2D0)
      BUOYANCYF = RHOW*(AREA*ALENGTH*ARATIO)
      ENDIF
      
      
      END
      
C     =========================================================================         
      SUBROUTINE FRAME_AREA (NELEMENT,GROWTH,AREA_OUT,OPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*2 OPT
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      
      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      IF (OPT.EQ."FN") THEN
      DPIPE_MARINE_OUT = DPIPE + GROWTH*2D0
      DPIPE_MARINE_IN  = DPIPE - TPIPE*2D0
      AREA_OUT = 3.14159265359D0*DPIPE_MARINE_OUT*DPIPE_MARINE_OUT/4D0 - 3.14159265359D0*DPIPE_MARINE_IN*DPIPE_MARINE_IN/4D0
      ELSEIF (OPT.EQ."UN") THEN
      DPIPE_MARINE = DPIPE + GROWTH*2D0
      AREA_OUT = 3.14159265359D0*DPIPE_MARINE*DPIPE_MARINE/4D0
      ELSEIF (OPT.EQ."MA") THEN ! MARINE GROWTH AREA
      DPIPE_MARINE_OUT = DPIPE + GROWTH*2D0
      DPIPE_MARINE_IN  = DPIPE 
      AREA_OUT = 3.14159265359D0*DPIPE_MARINE_OUT*DPIPE_MARINE_OUT/4D0 - 3.14159265359D0*DPIPE_MARINE_IN*DPIPE_MARINE_IN/4D0   
      ENDIF
      
      IF (NSECTION.NE.21) AREA_OUT = AREA
      
      END