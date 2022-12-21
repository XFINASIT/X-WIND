C     =========================================================================  
      SUBROUTINE DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM1,PERIOD,WVHIGHT,RK
     1                          ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                          ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
       
      IF (NFUNCTION.LE.3.0D0)THEN ! PREVENT ERROR FUNCTION FOR USER DEFINED
      ! -------- CALCULATE CD AND CM FOR AUTOMATIC PARTS ( BASE ON DNV )----------
      ! --------- DNV-OS-J101 JULY 2011 PAGE 74-75 - SEC.4  ---------
        ! NFUNCTION = 1 >> AUTOMATIC CACULATION ( DNV )
        ! NFUNCTION = 2 >> AUTOMATIC CACULATION ( API - ROUGH )
        ! NFUNCTION = 3 >> AUTOMATIC CACULATION ( API - SMOOTH )
        ! NFUNCTION = 4 >> USER DEFINED
        IF (NFUNCTION.EQ.1.0D0)THEN
           IF ((ROUGH/DIAM1).LT.0.0001)THEN
           CD  =  0.65 ! SMOOTH SURFACE
           ELSEIF ((ROUGH/DIAM1).GT.0.0001.AND.(ROUGH/DIAM1).LT.0.01)THEN
           CD  =  (29D0+4D0*LOG10(ROUGH/DIAM1))/20D0 
           ELSEIF ((ROUGH/DIAM1).GT.0.01)THEN
           CD  =  1.05 ! ROUGH SURFACE
           ENDIF
           IF (IWAVE.EQ.1.0D0)THEN ! AIRY WAVE THEORY
           ! GENERATE VELOCITY AT SURFACE 
           ! SET X POSITION =0.0, Y POSITION = WATER DEPTH, AND TIME =0.0
           UXMAX = 0.50D0*OMEGA*WVHIGHT*COSH(RK*WDEPTH)/SINH(RK*WDEPTH)*COS(RK*(0) - OMEGA*(TIME)) 
           UYMAX = 0.50D0*OMEGA*WVHIGHT*SINH(RK*WDEPTH)/SINH(RK*WDEPTH)*SIN(RK*(0) - OMEGA*(TIME)) 
           ! AKC = KEULEGAN CARPENTER SEE DNV-OS-J101 JULY 2011 PAGE 74-75 
           AKC   = UXMAX*PERIOD/DIAM1 
              IF (AKC.LT.3.0D0)THEN
              CM   =  2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  =  2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  =  1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM=CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM=CM2
                ENDIF
              ENDIF
           ! FOR WRITE STREAM FUNCTION FILE SEE XFINAS.FOR
           NSTREAMFUNCTION=0
           ELSEIF (IWAVE.EQ.2.0D0)THEN ! STOKE'S WAVE THEORY
           ! GENERATE VELOCITY AT SURFACE
           ! SET X POSITION =0.0, Y POSITION = WATER DEPTH, AND TIME =0.0
                 CALL PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55)     
           CALL STOKE_COEFFICIENT(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
C	     -----------------------------------
           G1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
           G2 = 2.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
           G3 = 3.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
           G4 = 4.0D0*(AVAL**4.0)*G44
           G5 = 5.0D0*(AVAL**5.0)*G55
C	     -----------------------------------
C	     DRAG FORCE        
           UXMAX  =  (G1*COSH(1.0D0*RK*WDEPTH)/SINH(1.0D0*RK*WDEPTH)*COS(1.0D0*(RK*(0)
     1     - OMEGA*(TIME)))
     1     + G2*COSH(2.0D0*RK*WDEPTH)/SINH(2.0D0*RK*WDEPTH)*COS(2.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G3*COSH(3.0D0*RK*WDEPTH)/SINH(3.0D0*RK*WDEPTH)*COS(3.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G4*COSH(4.0D0*RK*WDEPTH)/SINH(4.0D0*RK*WDEPTH)*COS(4.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G5*COSH(5.0D0*RK*WDEPTH)/SINH(5.0D0*RK*WDEPTH)*COS(5.0D0*(RK*(0) 
     1     - OMEGA*(TIME))))
     1      *(OMEGA/RK)
           UYMAX  =  (G1*SINH(1.0D0*RK*Y)/SINH(1.0D0*RK*WDEPTH)*SIN(1.0D0*(RK*(0) 
     1     - OMEGA*(TIME))) 
     1     + G2*SINH(2.0D0*RK*WDEPTH)/SINH(2.0D0*RK*WDEPTH)*SIN(2.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G3*SINH(3.0D0*RK*WDEPTH)/SINH(3.0D0*RK*WDEPTH)*SIN(3.0D0*(RK*(0) 
     1     - OMEGA*(TIME))) 
     1     + G4*SINH(4.0D0*RK*WDEPTH)/SINH(4.0D0*RK*WDEPTH)*SIN(4.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G5*SINH(5.0D0*RK*WDEPTH)/SINH(5.0D0*RK*WDEPTH)*SIN(5.0D0*(RK*(0) 
     1     - OMEGA*(TIME))))
     1      *(OMEGA/RK)
           
      !==================================
      !    MODIFY BY CHANA 23 FEB 2016
      !==================================
      
      !WVHIGHT = H
      !WDEPTH  = HW
      !GRAV = G
           
      chana = 3
           
      RK=(2*(AVAL+((AVAL**3)*F33)+((AVAL**5)*(F35+F55))))/WVHIGHT 
      
      CWAVE = SQRT(GRAV/RK*(1.0D0 + C1*(AVAL**2.0) + C2*AVAL**4.0)*TANH(RK*WDEPTH))
      
      chana = (OMEGA/RK)
      
      FG1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
      FG2 = 1.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
      FG3 = 1.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
      FG4 = 1.0D0*(AVAL**4.0)*G44
      FG5 = 1.0D0*(AVAL**5.0)*G55
      
      Y = WDEPTH
      X = 0.0d0
      
      UXMAX = CWAVE*( 1.0D0*FG1*COSH(1.0D0*RK*Y)*COS(1.0D0*(RK*X - OMEGA*TT))
     1            + 2.0D0*FG2*COSH(2.0D0*RK*Y)*COS(2.0D0*(RK*X - OMEGA*TT))
     1            + 3.0D0*FG3*COSH(3.0D0*RK*Y)*COS(3.0D0*(RK*X - OMEGA*TT))
     1            + 4.0D0*FG4*COSH(4.0D0*RK*Y)*COS(4.0D0*(RK*X - OMEGA*TT))
     1            + 5.0D0*FG5*COSH(5.0D0*RK*Y)*COS(5.0D0*(RK*X - OMEGA*TT)))
      
      UYMAX = CWAVE*( 1.0D0*FG1*SINH(1.0D0*RK*Y)*SIN(1.0D0*(RK*X - OMEGA*TT))
     1            + 2.0D0*FG2*SINH(2.0D0*RK*Y)*SIN(2.0D0*(RK*X - OMEGA*TT))
     1            + 3.0D0*FG3*SINH(3.0D0*RK*Y)*SIN(3.0D0*(RK*X - OMEGA*TT))
     1            + 4.0D0*FG4*SINH(4.0D0*RK*Y)*SIN(4.0D0*(RK*X - OMEGA*TT))
     1            + 5.0D0*FG5*SINH(5.0D0*RK*Y)*SIN(5.0D0*(RK*X - OMEGA*TT)))
      
      chana = 3
      
           AKC  =  UXMAX*PERIOD/DIAM1
              IF (AKC.LT.3.0D0)THEN
              CM  =  2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  =  2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  =  1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM =  CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM =  CM2
                ENDIF
              ENDIF 
           ! FOR WRITE STREAM FUNCTION FILE SEE XFINAS.FOR
           NSTREAMFUNCTION=0
           ELSEIF (IWAVE.EQ.3.0D0)THEN
           ! GENERATE VELOCITY AT SURFACE
           ! SET X POSITION =0.0, Y POSITION = WATER DEPTH, AND TIME =0.0
            CALL VELOC (0,WDEPTH,UXMAX,UYMAX,DUMAX,DVMAX,TIME,ARAMDA)
            AKC  =  UXMAX*PERIOD/DIAM1
           ! FOR WRITE STREAM FUNCTION FILE SEE XFINAS.FOR
           NSTREAMFUNCTION=1
              IF (AKC.LT.3.0D0)THEN
              CM   =  2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  =  2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  =  1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM =  CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM =  CM2
                ENDIF
              ENDIF 
           ELSEIF (IWAVE.EQ.4.0D0)THEN
           ZZETAA = wavenumber * X - 2d0 * AKKK/PERIOD * TIME
           CALL Jacobian_elliptic_function (ZZETAA,Ammm,CN,CN2,SSN,ddn)
           CALL Cnoidal_Wave_Kinematic  (WVHIGHT,WDEPTH,GRAV,PERIOD,TIME,WDEPTH,wavenumber,AKKK,ALAMDA,Ammm,AEEE,CN,SSN,ddn,ZZETAA
     1                                 ,UXMAX,UYMAX,DUDT,DVDTM)
              AKC = UXMAX*PERIOD/DIAM1
              IF (AKC.LT.3.0D0)THEN
              CM   = 2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  = 2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  = 1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM = CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM = CM2
                ENDIF
              ENDIF 
     
           ELSEIF (IWAVE.EQ.5.0D0)THEN
           ! SOLITARY WAVE    
           Z = WDEPTH
           Call Solitary_Wave (WVHIGHT,WDEPTH,GRAV,TIME,0d0,Z,UXMAX,UYMAX,au,av,WNU,q)
              
              AKC = UXMAX*PERIOD/DIAM1
              IF (AKC.LT.3.0D0)THEN
              CM   = 2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  = 2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  = 1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM = CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM = CM2
                ENDIF
              ENDIF 
           
           ENDIF
        ELSEIF (NFUNCTION.EQ.2.0D0.OR.NFUNCTION.EQ.3.0D0)THEN ! BASE ON API
           ! API 2A-WSD(RP 2A-WSD) PAGE 15
           IF (NFUNCTION.EQ.2.0D0)THEN ! ROUGH SURFACE
           CD  = 1.05D0
           CM  = 1.20D0
           ELSEIF (NFUNCTION.EQ.3.0D0)THEN ! SMOOTH SURFACE
           CD  = 0.65D0
           CM  = 1.20D0
           ENDIF
        ENDIF ! (ENDIF : END AUTOMATIC GENERATE CD,CM FUNCTION >> 1,DNV 2,API  )
      ENDIF   ! (ENDIF : PREVENT ERROR FUNCTION FOR USER DEFINED )

        END
        
C	=======================================================================
      SUBROUTINE TAPPEROFFSHORE (DIAM3,MLE,XYZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
      COMMON /OFFAREA/ AREA
     
      DIMENSION XYZ(NCO*NNM,NELE),VR(3),BPG(12),BWG(12),GPL(12),GPW(12),DIAM3(12)
      
      DO IGR = 1,12
	IF(IGR.EQ.1  )THEN
	GPL(IGR) = -1.0D0
	ENDIF
	
	IF(IGR.EQ.12 )THEN 
	GPL(IGR) =  1.0D0 
	ENDIF
	
	IF(IGR.NE.1.AND.IGR.NE.12)THEN
	GPL(IGR) =  GAUSP(IGR-1,12-2)
	ENDIF
	
	IF(IGR.EQ.1  )THEN
	GPW(IGR) =  0.0D0
	ENDIF
	
	IF(IGR.EQ.12 )THEN 
	GPW(IGR) =  0.0D0 
	ENDIF
	
	IF(IGR.NE.1.AND.IGR.NE.12)THEN
	GPW(IGR) =  GAUSW(IGR-1,12-2)
	ENDIF
	
	ENDDO
      
      VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)
	
      DO IGR = 1,12
	RI = GPL(IGR)  !GAUSP(IGR,NGR)
	RW = GPW(IGR)  !GAUSW(IGR,NGR)
	BPG(IGR) = 0.5*ELN*(1.0 + RI)
	BWG(IGR) = 0.5*ELN*RW
      ENDDO
      
     
	
      CALL CALLENGTH (MLE,VREW,0,ALENGTH,AMIDX,AMIDY,AMIDZ,BUOCYX,BUOCYY,BUOCYZ,NOP
     1                ,ATMIDX,ATMIDY,ATMIDZ)
      
      
      CALL MAXSECTION (MLE,DOUT,DIN,DOUT1,DIN2)  
      

      
      DO I=1,12
       IF (DOUT.NE.DOUT1)THEN ! TAPPER SECTION
       A        = (DOUT-DOUT1)/ALENGTH
       DIAM3(I) = DOUT-A*BPG(I)
       ELSEIF (DOUT.EQ.DOUT1)THEN ! NORMAL SECTION
       DIAM3(I) = DOUT
       ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE

C	=======================================================================
C	=======================================================================
	SUBROUTINE OFFREFAXIS(N,XXR,YYR,ZZR,V,NAME) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 NAME
      
C     INPUT      
C     N = COEFFICIENT GROUP NUMBER
C     OUTPUT
C     REFERENCE POSITION H1
C     REFERENCE POSITION H2
C     REFERENCE AXIS V(3)
      
      DIMENSION V(3)
      
      CALL RELFILL(NAME,XXR ,1,N,0)
      CALL RELFILL(NAME,YYR ,2,N,0)
      CALL RELFILL(NAME,ZZR ,3,N,0)
      CALL RELFILL(NAME,V(1),4,N,0)
      CALL RELFILL(NAME,V(2),5,N,0)
      CALL RELFILL(NAME,V(3),6,N,0)      
      
	RETURN
	END
C	=======================================================================
C	=======================FRAME THERMAL LOAD==============================
C	=======================================================================
	SUBROUTINE OFFSHSTEP(TIME,ITIME,NTIME,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 NAME,OPER
      
C     READ & CALL TIME DATA FOR OFFSHORE LOAD
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC  
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL
      COMMON /TIME/ DDT,CTIM,NINC
      
      COMMON /ChanaSpectrum/ Start , AInterval , Aend   !Chana 2014 RESPONSE SPECTRUM FOR ISOLOP 1
      COMMON /ChanaTImeSeries/ DT  
      COMMON /Spectrum_STEP/ NSTEPFA !Chana 2014 RESPONSE SPECTRUM FOR ISOLOP 1
      
      ALLOCATABLE TLIST(:)
      
      
      NAME = 'OFTM' !OFFSHORE TIME DATA
      
      
C     ============================================   
      SELECTCASE(OPER)
C     ============================================   
      
C     --------------------------------------------
      CASE('READ') !READ & GENERATE TIME DATA
C     --------------------------------------------

C     INPTM = INPUT METHOD 
C           = 0 FOR MANUALLY SPECIFY NUMBER OF TIME STEP (NTIME)
C           = 1 FOR AUTOMATICALLY CALCULATE STEP BY DATA INPUT FROM DYNAMIC ANALYSIS
      
      READ(ITI,*)
      READ(ITI,*) INPTM,NTIME,DT 
      
      IF(NTIME.EQ.0) THEN
      WRITE(*,*) 'ERROR IN OFFSHORE LOAD ANALYSIS'
          SELECTCASE(INPTM)
          CASE(0) 
            WRITE(*,*) 'PLEASE CHOOSE AT LEASE ONE TIME STEP OF ANALYSIS'
          CASE(1) 
            WRITE(*,*) 'NUMBER OF ANALYSIS TIME STEP SHOULD BE GREATER THAN ZERO'
          ENDSELECT   
      STOP 
      ENDIF
      
      
      ALLOCATE(TLIST(NTIME)) 
      
      
      SELECTCASE(INPTM)
      CASE(0) 
          READ(ITI,*)
          DO I = 1,NTIME
            READ(ITI,*) IT,TLIST(I) !TIME TO CALCULATE OFFSHORE LOAD
          ENDDO
      CASE(1)
          TIM = 0.0
          DO I = 1,NTIME
              TIM = TIM + DT
              TLIST(I) = TIM !TIME TO CALCULATE OFFSHORE LOAD
          ENDDO
      ENDSELECT
     
      
C     --------------------------------
C     ALLOCATE STORE DATA ---- NTIME     
	CALL DEFNREL(NAME,KTIME,1,NTIME)

C     --------------------------------
C     STORE DATA 
      DO I = 1,NTIME
	CALL RELFILL(NAME,TLIST(I),1,I,1)
	ENDDO
C     --------------------------------

      DEALLOCATE(TLIST)
      
      CASE('RSTOR')
      

C     --------------------------------------------    
      CASE('CALL') !CALL TIME DATA
C     --------------------------------------------
      
	CALL RELFILL(NAME,TIME,1,ITIME,0) !CALL TIME


      
C     --------------------------------------------    
      CASE('CALT') !CALL NTIME NUMBER OF TOTAL TIME STEP FOR OFFSHORE LOAD (MUST BE SAME WITH TIME STEP OF DYNAMIC ANALYSIS)
C     --------------------------------------------
    
	CALL LOCATN (NAME,KTIME,NMAX,NTIME,2)  !CALL NTIME
	
      CASE ('REDD')
    
      READ(ITI,*) 
      READ(ITI,*) INPTM,NTIME,DT 
      
      IF (ISOLOP.EQ.5)THEN
      NTIME   =  ((ABS(STEPEND)-ABS(STEPSTAT))/STEPINCR)+1
      DT      =  STEPINCR
C     FATIGUE ANALYSIS    
      Start     = STEPSTART
      AInterval = DT
      Aend      = STEPEND
      NSTEPFA   = NTIME
      IF(NTIME.EQ.0) THEN
      WRITE(*,*) 'ERROR IN OFFSHORE LOAD ANALYSIS'
          SELECTCASE(INPTM)
          CASE(0) 
            WRITE(*,*) 'PLEASE CHOOSE AT LEASE ONE TIME STEP OF ANALYSIS'
          CASE(1) 
            WRITE(*,*) 'NUMBER OF ANALYSIS TIME STEP SHOULD BE GREATER THAN ZERO'
          ENDSELECT   
      STOP 
      ENDIF
      ENDIF
      
      ALLOCATE(TLIST(NTIME))
       
      SELECTCASE(INPTM)
      CASE(0) 
          READ(ITI,*)
          DO I = 1,NTIME
            READ(ITI,*) IT,TLIST(I) !TIME TO CALCULATE OFFSHORE LOAD
          ENDDO
      CASE(1)
          TIM = 0.0
          DO I = 1,NTIME
              TIM = TIM + DT
              TLIST(I) = TIM !TIME TO CALCULATE OFFSHORE LOAD
          ENDDO
      ENDSELECT
      ! ----- OLD DATA 15-01-2018 --------
         ! TIM = STEPSTAT
         ! DO I = 1,NTIME
         !     IF (I.EQ.1) THEN 
         !     TIM = STEPSTAT
         !     ENDIF
         !     IF (I.GT.1) THEN 
         !     TIM = TIM + DT
         !     ENDIF
         !     TLIST(I) = TIM !TIME TO CALCULATE OFFSHORE LOAD
         ! ENDDO
      ! ---------------------------------
      
C     --------------------------------
C     ALLOCATE STORE DATA ---- NTIME     
	CALL DEFNREL(NAME,KTIME,1,NTIME)

C     --------------------------------
C     STORE DATA 
      DO I = 1,NTIME
	CALL RELFILL(NAME,TLIST(I),1,I,1)
	ENDDO
C     --------------------------------
	
	
	DEALLOCATE(TLIST)
C     ============================================   
      ENDSELECT
C     ============================================   
      
      
	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE OFFSHFILE(NEQ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      NFLV = 317
      NBIT = 2*NEQ
	CALL DIROPEN(NFLV,NBIT) 
	
      NFLC = 318
      NBIT = 2*NEQ
	CALL DIROPEN(NFLC,NBIT) 
	
	NFLC = 319
	CALL DIROPEN(NFLC,NBIT) 
      
	NFLC = 320
	CALL DIROPEN(NFLC,NBIT) 
	RETURN
	END
C	=======================================================================