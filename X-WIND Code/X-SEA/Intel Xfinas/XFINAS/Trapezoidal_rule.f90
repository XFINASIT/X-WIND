      SUBROUTINE TAPPEROFFSHORETAPIZOIDAL (DIAMOUTCAL,MLE,XYZ,ILOOPXX,INTERVALTPZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
      COMMON /OFFAREA/ AREA
     
      DIMENSION XYZ(NCO*NNM,NELE),VR(3),BPG(12),BWG(12),GPL(12),GPW(12),DIAM3(INTERVALTPZ)
      
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREAT,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
	
      CALL CALLENGTH (MLE,VREW,0,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                ,ATMIDX,ATMIDY,ATMIDZ)

      CALL MAXSECTION (MLE,DOUT,DIN,DOUT1,DIN2)  
      
      IF(DOUT.EQ.DOUT1.AND.DIN.EQ.DIN2)THEN
      CALL ARRAYELEMENT (MLE,NSECTION)   
      DIAMOUTCAL = DPIPE/2D0
      !DIAMOUTCAL = DOUT
      
      ELSE
          
          DIFFDIA = (  DOUT1 - DOUT  )

          DIAMOUTCAL = DOUT + DIFFDIA/INTERVALTPZ * ( ILOOPXX - 1 )

      ENDIF
      
      
      RETURN
      END SUBROUTINE
!================================================================================================  
      SUBROUTINE TAPPEROFFSHORETAPIZOIDAL_WIND (DIAMOUTCAL,MLE,XYZ,ILOOPXX,INTERVALTPZ,IELETPZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
      COMMON /OFFAREA/ AREA
     
      DIMENSION XYZ(NCO*NNM,NELE),VR(3),BPG(12),BWG(12),GPL(12),GPW(12),DIAM3(INTERVALTPZ)
      
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREAT,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
	
      CALL CALLENGTH (MLE,VREW,0,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                ,ATMIDX,ATMIDY,ATMIDZ)

      CALL MAXSECTION (MLE,DOUT,DIN,DOUT1,DIN2)  
      
      IF(DOUT.EQ.DOUT1.AND.DIN.EQ.DIN2)THEN
      CALL ARRAYELEMENT (MLE,NSECTION)   
      DIAMOUTCAL = DPIPE/2D0
      !DIAMOUTCAL = DOUT
      
      ELSE
          
          DIFFDIA = (  DOUT1 - DOUT  )

          IF (IELETPZ.EQ.1) DIAMOUTCAL = DOUT + DIFFDIA/INTERVALTPZ * ( ILOOPXX - 1 )
          IF (IELETPZ.EQ.2) DIAMOUTCAL = DOUT + DIFFDIA/INTERVALTPZ * ( ILOOPXX )

      ENDIF
      
      
      RETURN
      END SUBROUTINE
!================================================================================================      
      SUBROUTINE TRAPIZOIDALCENTROID_YYWAVE( FX1,FY1,FZ1,
     1                                FX2,FY2,FZ2,
     1                              MLE,XYZ,INDEXSTRIP,INTERVALTPZ,CENTROID) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
      COMMON /OFFAREA/ AREA
      COMMON /MGRAV/ NGRAV   
     
      DIMENSION XYZ(NCO*NNM,NELE)
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
       ! SELECTCASE(NGRAV)
       ! CASE(1)
       !   DIFFINTERVAL = DiffXELE/INTERVALTPZ
       ! CASE(2)
       !   
       !   PROJECTEDLENGTH = SQRT(  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
       !     
       !   !DIFFINTERVAL = DiffYELE/INTERVALTPZ
       !   
       !   DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
       !   
       !   !ELEVELEMENTLOCAL1 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
       !   !ELEVELEMENTLOCAL2 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
       !   ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
       !   ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
       !   
       ! CASE(3)
       !   
       !   PROJECTEDLENGTH = SQRT(  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
       !     
       !   !DIFFINTERVAL = DiffYELE/INTERVALTPZ
       !   
       !   DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
       !   
       !   !ELEVELEMENTLOCAL1 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
       !   !ELEVELEMENTLOCAL2 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
       !   ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
       !   ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
       ! ENDSELECT
      
      PROJECTEDLENGTH = SQRT(DiffXELE**2D0 +  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
      ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      
      IF(FY1.GT.FY2)THEN
          
          BB = FY1
          AA = FY2
          HH = DIFFINTERVAL
          
          ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
          
          CENTROID = ELEVELEMENTLOCAL1 + ZBAR
                   
      ELSEIF(FY1.EQ.FY2)THEN
          
          CENTROID = ELEVELEMENTLOCAL1 + DIFFINTERVAL / 2.0D0
          
      ELSEIF(FY1.LT.FY2)THEN

          BB = FY2
          AA = FY1
          HH = DIFFINTERVAL
          
          ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
          
          CENTROID = ELEVELEMENTLOCAL1 + HH - ZBAR
            
      ENDIF
      
      RETURN
      END SUBROUTINE
!================================================================================================
      SUBROUTINE TRAPIZOIDALCENTROID_XXWAVE( FX1,FY1,FZ1,
     1                                FX2,FY2,FZ2,
     1                              MLE,XYZ,INDEXSTRIP,INTERVALTPZ,CENTROID) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
      COMMON /OFFAREA/ AREA
      COMMON /MGRAV/ NGRAV   
     
      DIMENSION XYZ(NCO*NNM,NELE)
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      ! SELECTCASE(NGRAV)
      ! CASE(1)
      !   DIFFINTERVAL = DiffXELE/INTERVALTPZ
      !   
      ! CASE(2)
      !     
      !   PROJECTEDLENGTH = SQRT(  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      !     
      !   !DIFFINTERVAL = DiffYELE/INTERVALTPZ
      !   
      !   DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
      !   
      !   !ELEVELEMENTLOCAL1 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !   !ELEVELEMENTLOCAL2 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
      !   ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !   ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      !   
      ! CASE(3)
      !     
      !   PROJECTEDLENGTH = SQRT(  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      !   
      !   !DIFFINTERVAL = DiffZELE/INTERVALTPZ
      !   
      !   DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
      !   
      !   !ELEVELEMENTLOCAL1 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !   !ELEVELEMENTLOCAL2 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
      !   ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !   ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      ! ENDSELECT
      ! 
      ! PROJECTEDLENGTH = SQRT(DiffXELE**2D0 +  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      ! DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
      ! ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      ! ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      
      PROJECTEDLENGTH = SQRT(DiffXELE**2D0 +  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
      ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
          
          
       IF(FX1.GT.FX2)THEN
           
           BB = FX1
           AA = FX2
           HH = DIFFINTERVAL
           
           ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
           
           CENTROID = ELEVELEMENTLOCAL1 + ZBAR
                    
       ELSEIF(FX1.EQ.FX2)THEN
           
           CENTROID = ELEVELEMENTLOCAL1 + DIFFINTERVAL / 2.0D0
           
       ELSEIF(FX1.LT.FX2)THEN
      
           BB = FX2
           AA = FX1
           HH = DIFFINTERVAL
           
           ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
           
           CENTROID = ELEVELEMENTLOCAL1 + HH - ZBAR
             
       ENDIF
      
      RETURN
      END SUBROUTINE
!================================================================================================      

      SUBROUTINE TRAPIZOIDALCENTROID_ZZWAVE( FZ1,
     1                                FZ2,
     1                              MLE,XYZ,INDEXSTRIP,INTERVALTPZ,RCENTROID,CENTROID) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
      COMMON /OFFAREA/ AREA
      COMMON /MGRAV/ NGRAV   
     
      DIMENSION XYZ(NCO*NNM,NELE)
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      !SELECTCASE(NGRAV)
      !CASE(1)
      !  DIFFINTERVAL = DiffXELE/INTERVALTPZ
      !CASE(2)
      !  DIFFINTERVAL = SQRT(DiffXELE**2.0D0 + DiffYELE**2.0D0 )/INTERVALTPZ
      !  !ELEVELEMENTLOCAL1 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !  !ELEVELEMENTLOCAL2 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
      !  ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !  ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      !CASE(3)
      !  DIFFINTERVAL = SQRT(DiffXELE**2.0D0 + DiffYELE**2.0D0 )/INTERVALTPZ
      !  !ELEVELEMENTLOCAL1 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !  !ELEVELEMENTLOCAL2 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
      !  ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      !  ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      !ENDSELECT
      
      PROJECTEDLENGTH = SQRT(DiffXELE**2D0 +  DiffYELE**2.0D0 + DiffZELE**2.0D0 )
      DIFFINTERVAL = PROJECTEDLENGTH/INTERVALTPZ
      ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
      ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      
      IF(FZ1.GT.FZ2)THEN
          
          BB = FZ1
          AA = FZ2
          HH = DIFFINTERVAL
          
          ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
          
          RCENTROID = ELEVELEMENTLOCAL1 + ZBAR
          
          CENTROID = ZBAR
                   
      ELSEIF(FZ1.EQ.FZ2)THEN
          
          RCENTROID = ELEVELEMENTLOCAL1 + DIFFINTERVAL / 2.0D0
          
          CENTROID = DIFFINTERVAL / 2.0D0
          
      ELSEIF(FZ1.LT.FZ2)THEN

          BB = FZ2
          AA = FZ1
          HH = DIFFINTERVAL
          
          ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
          
          RCENTROID = ELEVELEMENTLOCAL1 + HH - ZBAR
          
          CENTROID = HH - ZBAR
            
      ENDIF
      
      RETURN
      END SUBROUTINE
!================================================================================================      
      SUBROUTINE TRAPIZOIDALFORCEandMOMENTtransformation( FX1,FY1,FZ1,
     1                              FX2,FY2,FZ2,
     1                              FX,FY,FZ,
     1                              RREX1,RREX2,RMEX1,RMEX2,
     1                              RREY1,RREY2,RMEY1,RMEY2,
     1                              RREZ1,RREZ2,RMEZ1,RMEZ2,
     1                              MLE,XYZ,INDEXSTRIP,INTERVALTPZ,CENTROID) 
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
      COMMON /OFFAREA/ AREA
      COMMON /MGRAV/ NGRAV   
     
      DIMENSION XYZ(NCO*NNM,NELE)
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      SELECTCASE(NGRAV)
      CASE(1)
        DIFFINTERVAL = DiffXELE/INTERVALTPZ
      CASE(2)
        DIFFINTERVAL = DiffYELE/INTERVALTPZ
        !ELEVELEMENTLOCAL1 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
        !ELEVELEMENTLOCAL2 = XYZ(2,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
        ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
        ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      CASE(3)
        DIFFINTERVAL = DiffZELE/INTERVALTPZ
        !ELEVELEMENTLOCAL1 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP - 1 )
        !ELEVELEMENTLOCAL2 = XYZ(3,MLE) + DIFFINTERVAL*(INDEXSTRIP  )
        ELEVELEMENTLOCAL1 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP - 1 )
        ELEVELEMENTLOCAL2 = 0.0D0 + DIFFINTERVAL*(INDEXSTRIP  )
      ENDSELECT
      

      
      IF(FX1.GT.FX2)THEN
          
          BB = FX1
          AA = FX2
          HH = DIFFINTERVAL
          
          ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
          
          CENTROID = ELEVELEMENTLOCAL1 + ZBAR
                   
      ELSEIF(FX1.EQ.FX2)THEN
          
          CENTROID = ELEVELEMENTLOCAL1 + DIFFINTERVAL / 2.0D0
          
      ELSEIF(FX1.LT.FX2)THEN

          BB = FX2
          AA = FX1
          HH = DIFFINTERVAL
          
          ZBAR = HH/3.0D0*(BB+2.0D0*AA)/(BB+AA)
          
          CENTROID = ELEVELEMENTLOCAL1 + HH - ZBAR
            
      ENDIF
      
      ALENGTH = DiffZELE
      
      FACTOR_RR1 = (CENTROID**2) * ( CENTROID + 3.0D0 * (ALENGTH-CENTROID) ) /ABS(ALENGTH**3.0D0)
      FACTOR_RR2 = ((ALENGTH-CENTROID)**2) * ( 3.0D0 * CENTROID +  (ALENGTH-CENTROID) ) /ABS(ALENGTH**3.0D0)
      
      FACTOR_RM1 = CENTROID * ((ALENGTH-CENTROID)**2)/ABS(ALENGTH**2.0D0)
      FACTOR_RM2 = -(CENTROID**2)  * (ALENGTH-CENTROID)/ABS(ALENGTH**2.0D0)
      
      RREX1 =  FX * FACTOR_RR1
      RREX2 =  FX * FACTOR_RR2
      
      RMEX1 =  FX * FACTOR_RM1
      RMEX2 =  FX * FACTOR_RM2
      
      RREY1 =  FY * FACTOR_RR1
      RREY2 =  FY * FACTOR_RR2
      
      RMEY1 =  FY * FACTOR_RM1
      RMEY2 =  FY * FACTOR_RM2
      
      RREZ1 =  FZ * FACTOR_RR1
      RREZ2 =  FZ * FACTOR_RR2
      
      RMEZ1 =  FZ * FACTOR_RM1
      RMEZ2 =  FZ * FACTOR_RM2
      
      
      !RRX = RREX1 + RREX2
      !RMX = RMEX1 + RMEX2
      

      RETURN
      END SUBROUTINE
!================================================================================================      
      
      SUBROUTINE TRAPIZOIDAL_REACTION_LCS( ILCAS, LCS , NLCS, OPTION )
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPTION
      
      COMMON /TRAPIZOIDALLCS/ LCSANA(100), NLCSDATA(100)
      COMMON /TRAPIZOIDALLCSS/ LCSTOTAL
      

          
      IF(OPTION.EQ.'WRT')THEN
              NLCSDATA(ILCAS) = 1
      ENDIF
      
      IF(OPTION.EQ.'RED')THEN
          
          NLCS = NLCSDATA(ILCAS) 

      ENDIF
      
      RETURN
      END SUBROUTINE
!================================================================================================     
!================================================================================================
      
!      SUBROUTINE TRAPIZOIDAL_TRANSFORMATION(VR,VS,VT,OPTION,WFLOCAL)
!      
!      IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER*4 (I-N)
!      CHARACTER*3 OPTION
!      
!      COMMON /TRAPIZOIDALLCS/ LCSANA(100), NLCSDATA(100)
!      COMMON /TRAPIZOIDALLCSS/ LCSTOTAL
!
!C	==================================================================
!C	TRAPIZOIDAL RULE   
!      COMMON / WAVEREACFIXPARAMETER /  WAVEREACFIX(9999,7,2)
!C	==================================================================      
!      COMMON / ELEMENT_LOCALWAVEREAC /  LOCALELEMENTWAVEVECTOR(9999,9)
!C     SAVE ELEMENT NUMBER      
!      COMMON /NELEM/ IEL  
!      
!      DIMENSION TT(14,14),FIXLRWAVEFRAME(14),TTT(14,14)
!      DIMENSION VR(3),VS(3),VT(3)
!      DIMENSION WFLOCAL(14)
!      
!      WFLOCAL = 0.0D0
!          
!      IF(OPTION.EQ.'RED')THEN
!          
!          VR(1:3)   = LOCALELEMENTWAVEVECTOR(IEL,1:3)
!          VS(1:3)   = LOCALELEMENTWAVEVECTOR(IEL,4:6)
!          VT(1:3)   = LOCALELEMENTWAVEVECTOR(IEL,7:9)
!
!      CALL TT1A (VR,VS,VT,TT)
! 
!      FIXLRWAVEFRAME(1)	= WAVEREACFIX(IEL,1,1)   
!      FIXLRWAVEFRAME(2)	= WAVEREACFIX(IEL,2,1)   
!      FIXLRWAVEFRAME(3)	= WAVEREACFIX(IEL,3,1)   
!      FIXLRWAVEFRAME(4)	= WAVEREACFIX(IEL,4,1)   
!      FIXLRWAVEFRAME(5)	= WAVEREACFIX(IEL,5,1)   
!      FIXLRWAVEFRAME(6)	= WAVEREACFIX(IEL,6,1)   
!      FIXLRWAVEFRAME(7)	= WAVEREACFIX(IEL,7,1)   
!      FIXLRWAVEFRAME(8)	= WAVEREACFIX(IEL,1,2)   
!      FIXLRWAVEFRAME(9)	= WAVEREACFIX(IEL,2,2)   
!      FIXLRWAVEFRAME(10)	= WAVEREACFIX(IEL,3,2)   
!      FIXLRWAVEFRAME(11)	= WAVEREACFIX(IEL,4,2)   
!      FIXLRWAVEFRAME(12)	= WAVEREACFIX(IEL,5,2)   
!      FIXLRWAVEFRAME(13)	= WAVEREACFIX(IEL,6,2)   
!      FIXLRWAVEFRAME(14)	= WAVEREACFIX(IEL,7,2) 
!
!      !TTT = TRANSPOSE(TT)
!      NEF = 14    
!      DO IEF = 1,NEF
!	!FIXLR(IEF) = 0.0
!	!FIXLO(IEF) = 0.0
!      WFLOCAL(IEF) = 0.0 
!	DO JEF = 1,NEF
!	WFLOCAL(IEF) = WFLOCAL(IEF) + TT(JEF,IEF)*FIXLRWAVEFRAME(JEF)  !VARY FIXEND
!      !FIXLR(IEF) = FIXLR(IEF) + TT(IEF,JEF)*FIXEN(JEF)  !VARY FIXEND
!	!FIXLO(IEF) = FIXLO(IEF) + TT(IEF,JEF)*FIXEO(JEF)  !CONT FIXEND
!	ENDDO
!      ENDDO
!
!      RETURN
!      ELSEIF(OPTION.EQ.'WRT')THEN
!          
!      LOCALELEMENTWAVEVECTOR(IEL,1:3) = VR(1:3)
!      LOCALELEMENTWAVEVECTOR(IEL,4:6) = VS(1:3)    
!      LOCALELEMENTWAVEVECTOR(IEL,7:9) = VT(1:3)
!      
!      ENDIF
!      
!      RETURN
!      END SUBROUTINE
!      
!================================================================================================  
!================================================================================================ 
      
      SUBROUTINE ELEMENT_CENTROID_X(INTERVALTPZ,INDEXTAPPER,WAVEF,CENTROIDHH,OPTION,TOTALFORCE,CENTROIDOUT)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      CHARACTER*3 OPTION
      
      DIMENSION SUMFX(999),SUMFY(999),SUMFZ(999)
      
      COMMON /TRAPIZ_CENTROIDI/ SUMFX,SUMFY,SUMFZ
      
      
      TOTALVALUE = 0
      
      IF(OPTION.EQ.'WRT') THEN
      
      SUMFX(INDEXTAPPER) = WAVEF * CENTROIDHH
      
      ELSEIF(OPTION.EQ.'RED') THEN
          
          DO I = 1, INTERVALTPZ
              
              TOTALVALUE = TOTALVALUE + SUMFX(I)
              
          ENDDO
          
           CENTROIDOUT = TOTALVALUE / TOTALFORCE
           
           IF(TOTALFORCE.EQ.0) CENTROIDOUT = 0.0D0
           

      ENDIF
      
      CHANA = 3
      
      RETURN
      END SUBROUTINE
      
!================================================================================================ 
!================================================================================================
      
      SUBROUTINE ELEMENT_CENTROID_Y(INTERVALTPZ,INDEXTAPPER,WAVEF,CENTROIDHH,OPTION,TOTALFORCE,CENTROIDOUT)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      CHARACTER*3 OPTION
      
      DIMENSION SUMFX(999),SUMFY(999),SUMFZ(999)
      
      COMMON /TRAPIZ_CENTROIDI/ SUMFX,SUMFY,SUMFZ
      
      
      TOTALVALUE = 0
      
      IF(OPTION.EQ.'WRT') THEN
      
      SUMFY(INDEXTAPPER) = WAVEF * CENTROIDHH
      
      ELSEIF(OPTION.EQ.'RED') THEN
          
          DO I = 1, INTERVALTPZ
              
              TOTALVALUE = TOTALVALUE + SUMFY(I)
              
          ENDDO
          
           CENTROIDOUT = TOTALVALUE / TOTALFORCE
           
           IF(TOTALFORCE.EQ.0) CENTROIDOUT = 0.0D0

      ENDIF
      
      CHANA = 3
      
      RETURN
      END SUBROUTINE
      
!================================================================================================
!================================================================================================
      
      SUBROUTINE ELEMENT_CENTROID_Z(INTERVALTPZ,INDEXTAPPER,WAVEF,CENTROIDHH,OPTION,TOTALFORCE,CENTROIDOUT)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      CHARACTER*3 OPTION
      
      DIMENSION SUMFX(999),SUMFY(999),SUMFZ(999)
      
      COMMON /TRAPIZ_CENTROIDI/ SUMFX,SUMFY,SUMFZ
      
      
      TOTALVALUE = 0
      
      IF(OPTION.EQ.'WRT') THEN
      
      SUMFZ(INDEXTAPPER) = WAVEF * CENTROIDHH
      
      ELSEIF(OPTION.EQ.'RED') THEN
          
          DO I = 1, INTERVALTPZ
              
              TOTALVALUE = TOTALVALUE + SUMFZ(I)
              
          ENDDO
          
           CENTROIDOUT = TOTALVALUE / TOTALFORCE
           
           IF(TOTALFORCE.EQ.0) CENTROIDOUT = 0.0D0

      ENDIF
      
      CHANA = 3
      
      RETURN
      END SUBROUTINE
!================================================================================================ 
!================================================================================================
      SUBROUTINE ELEMENT_RATIO(XYZ,OPTION,MLE,CENTROID,RATIO)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      CHARACTER*1 OPTION
      
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
      DIMENSION XYZ(NCO*NNM,NELE)
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      IF(OPTION.EQ.'X')THEN
      
      PROJECTLENGTH =  SQRT( DiffYELE**2.0D0 + DiffZELE**2.0D0 )   
      RATIO = CENTROID / PROJECTLENGTH
      IF( PROJECTLENGTH .EQ. 0 ) RATIO = 0.0D0
      
      ELSEIF(OPTION.EQ.'Y')THEN
          
      PROJECTLENGTH =  SQRT( DiffXELE**2.0D0 + DiffZELE**2.0D0 )   
      RATIO = CENTROID / PROJECTLENGTH
      IF( PROJECTLENGTH .EQ. 0 ) RATIO = 0.0D0
      
      ELSEIF(OPTION.EQ.'Z')THEN
      
      PROJECTLENGTH =  SQRT( DiffXELE**2.0D0 + DiffYELE**2.0D0 )
      RATIO = CENTROID / PROJECTLENGTH
      IF( PROJECTLENGTH .EQ. 0 ) RATIO = 0.0D0
          
      ENDIF
      
      RETURN
      END SUBROUTINE
!================================================================================================
!================================================================================================
      
      SUBROUTINE CUTTING_FORCE(XYZ,OPTION,MLE,FORCE)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      CHARACTER*1 OPTION
      
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
      DIMENSION XYZ(NCO*NNM,NELE)
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      IF(OPTION.EQ.'X')THEN
          
          IF( DiffYELE.EQ.0.0D0 .AND. DiffZELE.EQ.0.0D0 ) FORCE = 0.0D0
          

      ELSEIF(OPTION.EQ.'Y')THEN
          
          IF( DiffXELE.EQ.0.0D0 .AND. DiffZELE.EQ.0.0D0 ) FORCE = 0.0D0
          
      ELSEIF(OPTION.EQ.'Z')THEN
          
          IF( DiffXELE.EQ.0.0D0 .AND. DiffYELE.EQ.0.0D0 ) FORCE = 0.0D0
          
      ENDIF
      
      RETURN
      END SUBROUTINE

!================================================================================================
!================================================================================================
      SUBROUTINE  FORCE_EACH_SEGMENT(XYZ,OPTION,MLE,F1,F2,FORCEOUT,INDEXTAPPER )
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
      CHARACTER*1 OPTION
      
      COMMON /MGRAV/ NGRAV 
            
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      
      DIMENSION XYZ(NCO*NNM,NELE)
      
      PI  = 3.141592653589793
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
            
      AINCLINEFACTORX = ATAN(DiffZELE/DiffXELE)
      DEGREE_AINCLINEFACTORX = AINCLINEFACTORX*180.0D0/PI
      
      FACTOR = 1.0D0/SIN(AINCLINEFACTORX)
      
      
      
!      SELECTCASE(NGRAV)
!      CASE(3)    
!          
!      CALL INCLINEFACTOR(XYZ,MLE,OPTION,FACTOR)
!      
!      FACTOR2 = SQRT(DiffZELE**2.0D0 + DiffXELE**2.0D0)/ABS(DiffZELE)
!      !IF(DiffZELE.EQ.0.0D0
!      
!      WRITE(*,*) 'FACTOR', FACTOR,FACTOR2,DEGREE_AINCLINEFACTORX
!      
!      ENDSELECT
      
      
      IF(OPTION.EQ.'X')THEN
              
              !PROJECLENGTH = SQRT(DiffYELE**2.0D0 + DiffZELE**2.0D0 )
              
              PROJECLENGTH = SQRT(DiffXELE**2.0D0 + DiffYELE**2.0D0 + DiffZELE**2.0D0 )
              
              !AINCLINE_SHAPE = 1/SIN( ATAN(DiffZELE/DiffXELE))
              
              !CALL INCLINEFACTOR(XYZ,MLE,'X',FACTOR)
              
              FORCEOUT = 0.5D0 * ( PROJECLENGTH / INDEXTAPPER ) * ( F1 + F2 ) !* FACTOR  !AINCLINE_SHAPE**2

          
      ELSEIF(OPTION.EQ.'Y')THEN
          
              
              !PROJECLENGTH = SQRT(DiffXELE**2.0D0 + DiffZELE**2.0D0 )
              
              PROJECLENGTH = SQRT(DiffXELE**2.0D0 + DiffYELE**2.0D0 + DiffZELE**2.0D0 )
              
              FORCEOUT = 0.5D0 * ( PROJECLENGTH / INDEXTAPPER ) * ( F1 + F2 )
                        
          
      ELSEIF(OPTION.EQ.'Z')THEN
          
               !PROJECLENGTH = SQRT(DiffXELE**2.0D0 + DiffYELE**2.0D0 )
               
               PROJECLENGTH = SQRT(DiffXELE**2.0D0 + DiffYELE**2.0D0 + DiffZELE**2.0D0 )
         
              FORCEOUT = 0.5D0 * ( PROJECLENGTH / INDEXTAPPER ) * ( F1 + F2 )
              
      ENDIF    
      
      CHANA = 3

      RETURN
      END SUBROUTINE

!================================================================================================
!================================================================================================