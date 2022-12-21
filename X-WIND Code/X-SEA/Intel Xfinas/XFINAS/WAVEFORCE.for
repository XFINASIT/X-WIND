C	=======================================================================
C	=======================================================================
	SUBROUTINE OFFSHFORC(R,ICASE,ITIME,NEQ,TYP,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP,OPER
C     STORE & CALL LOAD VECTOR FOR EACH TIME STEP FOR OFFSHORE LOAD
	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      DIMENSION R(1)
      ALLOCATABLE RR(:)
      
      IF(KOFFL.EQ.0.AND.KSPEC.EQ.0) RETURN
        IF (OPRES.EQ.0) CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME HERE
        IF (OPRES.EQ.1) CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME HERE
        IF (OPRES.EQ.2) CALL SELECTGRAPH (1,NTIME,0,0)
      ALLOCATE(RR(NEQ))
      
      NFLV = 317
      NFLC = 318
      NFLS = 319  ! SPECTRAL
      NFLC = 320  ! FATIGUE
      
      SELECTCASE(TYP)
      CASE('VARY') !VARY LOAD INDEX
      NFL = NFLV
      CASE('CONT') !CONSTANT LOAD INDEX
      NFL = NFLC
      CASE('SPEC')
      NFL = NFLS
      ENDSELECT
      
      IF (KFATIGUE.NE.1)THEN
      LRCD = ITIME + NTIME*(ICASE-1)
      ELSEIF (KFATIGUE.EQ.1)THEN
      LRCD = ICASE
      ! LRCD = ITIME + NTIME*(ICASE-1)
      ENDIF
      
      
      SELECTCASE(OPER)
      
      CASE('WRIT')
      WRITE(NFL,REC=LRCD) R(1:NEQ)
      CASE('READ','RADD')
       READ(NFL,REC=LRCD,ERR=10) RR(1:NEQ)
      GOTO 20
10    CONTINUE
      RR(1:NEQ) = 0.0D0
      WRITE(NFL,REC=LRCD) RR(1:NEQ) !DEFAULT TO ZERO IF DATA DOES NOT EXIST    
      
      CASE('REDT')
      READ(NFL,REC=LRCD) R(1:NEQ)
      ENDSELECT
      
      
20    CONTINUE


      SELECTCASE(OPER)
      CASE('READ')
        R(1:NEQ) = RR(1:NEQ)
      CASE('RADD')
        R(1:NEQ) = R(1:NEQ) + RR(1:NEQ)
      ENDSELECT



      DEALLOCATE(RR)
      
	RETURN
	END
C	=======================================================================
      SUBROUTINE WAVE_PRESSURE(OMEGA,RATIO,H,HW,RK,RHO,CI,CD,D,X,Y,TM,IWAVE,ORDER,VR,AVAL,G,
     1                         VTIDE0,VWIND0,PTX,PTY,PTZ,VDUM,NC,H0,VCURRENTP,PW,LC,VCL,PERIOD,LWCASE,CS,
     1                         WAVEKINEMATICFACTOR,CURRENTBLOCKAGEFACTOR,WAVEKINEMATICFACTORBREAKING,WAVEBREAKINGKINEMATICFACTOR,
     1                         WAVEFORCECOEFFICIENTS,YLMIN,IORRE,IRWAVE)
 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	---------------------------------------------------
      DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
      DIMENSION AX(3),AY(3),ATOL(3),ANOL(3),TRANN(14,14),TRANS(3,3)
      DIMENSION VWS(3),VWZ(3),VDUM(3),VCL(5)
      DIMENSION LWCASE(7)
      PI = 3.141592653589793d0
C	---------------------------------------------------
C                       VELOCITY FORCE
C	---------------------------------------------------    
C	1) DRAG VELOCITY      
C       - FOR STREAM FUNCTION : OPERATION ALL IN THIS SUBROUTINE  
      IF (LWCASE(1).EQ.1) THEN
      CALL DRAG_VELOCITY(OMEGA,RATIO,H,HW,RK,X,Y,TM,IWAVE,ORDER,AVAL,UDX,UDY,
     1                   AXSTREAM,AYSTREAM,G,PERIOD,IORRE,IRWAVE,AXRAN,AYRAN)
       
      ! ----- ACCELERATION FROM STREM FUNCTION ------   
        AAX = AXSTREAM
        AAY = AYSTREAM  
      ! --------------------------------------------- 
      ENDIF
      
      ! ----- ACCELERATION FROM RANDOM WAVE ------ 
      IF (IORRE.EQ.2) THEN
        AAX=AXRAN
        AAY=AYRAN     
      ENDIF
      
C     ---------------------------------------------------      
C                GENERATE CURRENT PROFILES 
C     ---------------------------------------------------
C	2) CURRENT VELOCITY
      IF (NC.EQ.1D0)THEN ! DNV
       UT    = VTIDE0*((HW+Y)/HW)**(1.0/7.0)
       UW    = VWIND0*((H0+Y)/H0)
       UC    = UT+UW
      ELSEIF (NC.EQ.3D0)THEN !IEC
       UT    = VTIDE0*((HW+Y)/HW)**(1.0/7.0)
       UW    = VWIND0*((1+Y)/20D0)
       UC    = UT+UW  
      ELSEIF (NC.EQ.4)THEN ! POWER LAW
       UC    = VCURRENTP*(((HW+Y)/HW)**(PW))
      ELSEIF (NC.EQ.5.OR.NC.EQ.2)THEN ! API AND LINEAR CURRENT
       AM1   =  (VCL(1)-VCL(2))/(((HW*4D0)/4D0)-((HW*3D0)/4D0))
       AM2   =  (VCL(2)-VCL(3))/(((HW*3D0)/4D0)-((HW*2D0)/4D0))
       AM3   =  (VCL(3)-VCL(4))/(((HW*2D0)/4D0)-((HW*1D0)/4D0))
       AM4   =  (VCL(4)-VCL(5))/((HW*1D0)/4D0)
       A2    =   HW*3.0D0/4.0D0
       A3    =   HW*2.0D0/4.0D0
       A4    =   HW/4.0D0
         IF     (Y.LE.HW.AND.Y.GE.HW*3D0/4D0)         THEN
         UC  =  VCL(2)+(AM1*(Y-(HW*3D0/4D0)))
         ELSEIF (Y.LE.HW*3D0/4D0.AND.Y.GE.HW*2D0/4D0) THEN 
         UC  =  VCL(3)+(AM2*(Y-(HW*2D0/4D0)))
         ELSEIF (Y.LE.HW*2D0/4D0.AND.Y.GE.HW*1D0/4D0) THEN 
         UC  =  VCL(4)+(AM3*(Y-(HW*1D0/4D0)))
         ELSEIF (Y.LE.HW*1D0/4D0.AND.Y.GE.0.0D0)     THEN
         UC  =  VCL(5)+(AM4*(Y-(HW/4D0)))
         ENDIF 
      ENDIF
          
      
C     ---------------------------------------------------         
C	      TOTAL VELOCITY SELECT BY LOAD CASE     
C     ---------------------------------------------------
C     1) COMBINATION WAVE VELOCITY + CURRENT VELOCITY  
      IF (LWCASE(1).EQ.1.AND.LWCASE(2).EQ.1.0)THEN
      UUX = UDX*WAVEKINEMATICFACTOR + UC*CURRENTBLOCKAGEFACTOR
C     2) PURE WAVE VELOCITY 
      ELSEIF (LWCASE(1).EQ.1)THEN
      UUX = UDX*WAVEKINEMATICFACTOR
C     3) PURE CURRENT VELOCITY
      ELSEIF (LWCASE(2).EQ.1.0)THEN
      UUX = UC*CURRENTBLOCKAGEFACTOR
C     4) PURE WAVE BREAKLING ( PREPARE FOR CALCULATE WAVE VELOCITY BELOW MEAN SEA LEVEL )  
      ELSEIF (LWCASE(5).EQ.1.OR.LWCASE(6).EQ.1)THEN ! BREAKING WAVE
      UUX = UDX*WAVEKINEMATICFACTORBREAKING
      ENDIF
               
      
C     MAGNITUDE OF VELOCITY
      VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 
      
      
C	DRAG FORCE ( STD = SUM TOTAL DRAG FORCE ) Pamin transferforce
      
      
      CALL Diffraction_Force_Parameter(ak,Aka,Cm,Ch,Delta)   
      
      STD = 0.50D0*RHO*CD            
      PDX = STD*VV*UUX  
      PDY = STD*VV*UUY
      PDZ = STD*VV*UUZ 
      
      PDX = STD*ABS(UUX)*UUX  
      PDY = STD*ABS(UUY)*UUY
      PDZ = STD*ABS(UUZ)*UUZ 
      
      IF( IWAVE == 6 ) then
      
      CALL SHELL_Diffraction_Pressure(OMEGA,RATIO,H,HW,RK,RHO,CI,CD,D,X,Y,TM,IWAVE,ORDER,VR,AVAL,G,
     1                         VTIDE0,VWIND0,PDX,PDY,PDZ,VDUM,NC,H0,VCURRENTP,PW,LC,VCL,PERIOD,LWCASE,CS,
     1                         WAVEKINEMATICFACTOR,CURRENTBLOCKAGEFACTOR,WAVEKINEMATICFACTORBREAKING,WAVEBREAKINGKINEMATICFACTOR,
     1                         WAVEFORCECOEFFICIENTS,YLMIN)
            
      endif
      
C	NORMAL VELOCITY
      VX(1) = 1.0D0
      VX(2) = 0.0D0
      VX(3) = 0.0D0      
      VY(1) = 0.0D0
      VY(2) = 1.0D0
      VY(3) = 0.0D0     
      DO IV = 1,3      
      VX(IV)   = UUX*VX(IV)
      VY(IV)   = UUY*VY(IV)      
      VTOL(IV) = VX(IV) + VY(IV)
      ENDDO
      CALL VECTORCROSS(VTOL,VR,VZ)
      CALL VECTORCROSS(VR,VZ,VNOL)      
      CALL VECTORUNIT(VNOL,VDUM,VV)             
      UN = VNOL(1)  
      VN = VNOL(2)  
      WN = VNOL(3)    
        
                        
C	DRAG FORCE 
      STD = 0.50D0*RHO*CD            
      PDX = STD*VV*UN  
      PDY = STD*VV*VN
      PDZ = STD*VV*WN 
      
      IF( IWAVE == 6 ) then
      
      PDX = PDX*UN  
      PDY = PDY*VN
      PDZ = PDZ*WN 

            
      endif
      
      
      
      
C     ---------------------------------------------------
C                  PLUNGING AND SURGING       teay
C     --------------------------------------------------- 
      IF (LWCASE(5).EQ.1.OR.LWCASE(6).EQ.(1))THEN
        IF (LWCASE(1).EQ.0)THEN ! DISABLE WAVE ANALYSIS
           IF (Y.GT.YLMIN)THEN  
              IF (LWCASE(5).EQ.1)THEN
C	         FOR PLUNGING
               PDXP = 0.5D0*RHO*CS*((UDX*WAVEBREAKINGKINEMATICFACTOR)**2)
               PDYP = 0.0D0
               PDZP = 0.0D0 
              ELSEIF (LWCASE(6).EQ.1)THEN
C              FOR SURGING
               PDXS = 0.5D0*RHO*CS*((UDX*WAVEBREAKINGKINEMATICFACTOR)**2)
               PDYS = 0.0D0
               PDZS = 0.0D0
              ENDIF
           ENDIF 
C              PDX  = WAVE FORCE + PLUNGING + SURGING           
               PDX  = PDXP + PDXS + PDX
               PDY  = PDYP + PDYS + PDY
               PDZ  = PDZP + PDZS + PDZ
        ENDIF
      ENDIF      
       
C	---------------------------------------------------
C                    ACCELERATION FORCE  
C	---------------------------------------------------  
C	INERTIA FORCE        Parmain
      IF (IORRE.EQ.1) CALL DRAG_ACCELERATION(OMEGA,RATIO,H,HW,RK,X,Y,TM,IWAVE,ORDER,AVAL,G,AAX,AAY)        
      AAZ = 0.0D0
      
      IF (LWCASE(2).EQ.1.0D0.AND.LWCASE(1).EQ.0.0)THEN
      AAX = 0.0D0
      AAY = 0.0D0
      ENDIF       
      
C	INNERTIA FORCE 
      STI = RHO*CI*PI*D/4.0D0  
      PIX = STI*AAX
      PIY = STI*AAY
      PIZ = STI*AAZ  
      
C	NORMAL ACCELERATION
      AX(1) = 1.0D0
      AX(2) = 0.0D0
      AX(3) = 0.0D0
      AY(1) = 0.0D0
      AY(2) = 1.0D0
      AY(3) = 0.0D0
      DO IV = 1,3            
      AX(IV)   = AAX*AX(IV)
      AY(IV)   = AAY*AY(IV)      
      ATOL(IV) = AX(IV) + AY(IV)
      ENDDO
      CALL VECTORCROSS(ATOL,VR,AZ)
      CALL VECTORCROSS(VR,AZ,ANOL)                
      ANX = ANOL(1)  
      ANY = ANOL(2)  
      ANZ = ANOL(3)
      
C	--------------------------------------   
C     GRAVITY PROBLEM   
      !IF ((Y >= 1.085).AND.(Y <= 18.0)) THEN
      !D1 = -0.76850D0*(Y - 1.085) + 19.50D0      
      !ELSE 
      !D1 = 6.50D0
      !ENDIF
      !D2 = D1 - 1.0D0
C	--------------------------------------
C     FOUNDATION PROBLEM
      !IF ((Y >= 0.0).AND.(Y <= 10.0)) THEN
      !D1 = -0.65*Y + 14.0D0 
      !D2 = -0.15*Y + 8.0D0   
      !ELSE 
      !D1 = 6.50D0
      !D2 = D1 - 0.20D0
      !ENDIF      
C	--------------------------------------  
C     TYPE C
      IF ((Y >= 0.0).AND.(Y < 2.0)) THEN
      D1 = 16.0D0
      D2 = D1 - 4.10D0
      ELSEIF ((Y > 2.0).AND.(Y <= 9.0)) THEN
      D1 = -0.75*(Y-2.0) + 13.50D0 
      D2 = D1 - 1.10D0 
      ELSEIF ((Y > 9.0).AND.(Y <= 10.0)) THEN
      D1 = -0.75*(Y-9.0) + 8.250D0 
      D2 = 5.30D0         
      ELSE 
      D1 = 6.50D0
      D2 = D1 - 0.10D0
      ENDIF      
C	-------------------------------------- 
C     TYPE F
      !IF ((Y >= 0.0).AND.(Y <= 10.0)) THEN
      !D = 1.50D0     
      !STI = RHO*CI*PI*D/4.0D0  
      !ELSE 
      !D1 = 6.50D0
      !D2 = D1 - 0.10D0      
      !STI = RHO*CI*PI*(D1*D1-D2*D2)/(4.0D0*D1) 
      !ENDIF      
C	-------------------------------------- 
      
      
C	INNERTIA FORCE      
      STI = RHO*CI
      IF( IWAVE == 6 ) STI = RHO*Cm
      PIX = STI*ANX
      PIY = STI*ANY
      PIZ = STI*ANZ   
      
C	---------------------------------------------------
C                    TOTAL WAVE FORCE
C	---------------------------------------------------  
!      IF( IWAVE == 6 ) GOTO 900
      PTX = ( PDX + PIX ) * WAVEFORCECOEFFICIENTS
      PTY = ( PDY + PIY ) * WAVEFORCECOEFFICIENTS
      PTZ = ( PDZ + PIZ ) * WAVEFORCECOEFFICIENTS

      RETURN
C	--------------------------------------------------- 
C      WRITE (300,1) X,Y,PTX,PTY,PTZ
C1     FORMAT (E12.5,E12.5,E12.5,E12.5,E12.5) 
      
C     ====================================================
c      DIFFRECTION WAVE THOERY 
C     ====================================================
900      CALL SHELL_Diffraction_Pressure(OMEGA,RATIO,H,HW,RK,RHO,CI,CD,D,X,Y,TM,IWAVE,ORDER,VR,AVAL,G,
     1                         VTIDE0,VWIND0,PTX,PTY,PTZ,VDUM,NC,H0,VCURRENTP,PW,LC,VCL,PERIOD,LWCASE,CS,
     1                         WAVEKINEMATICFACTOR,CURRENTBLOCKAGEFACTOR,WAVEKINEMATICFACTORBREAKING,WAVEBREAKINGKINEMATICFACTOR,
     1                         WAVEFORCECOEFFICIENTS,YLMIN)
   
      ACVET = 10D0

      RETURN
	END      
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE DRAG_VELOCITY(OMEGA,RATIO,H,HW,RK,X,Y,TT,IWAVE,ORDER,AVAL,UDX,UDY,
     1                        AXSTREAM,AYSTREAM,G,PERIOD,IORRE,IRWAVE,AUX,AUY)
    
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N) 
      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
C	-----------------------------------------------------------------------       
      
      IF (IORRE.EQ.1)THEN ! REGULAR WAVE
      IF (IWAVE == 1 . OR . IWAVE == 6) THEN
      ! AIRY WAVE THEORY ( PROVE 13-02-2013 )
                  
      UDX = 0.50D0*OMEGA*H*COSH(RK*Y)/SINH(RK*HW)*COS(RK*X - OMEGA*TT) 
      UDY = 0.50D0*OMEGA*H*SINH(RK*Y)/SINH(RK*HW)*SIN(RK*X - OMEGA*TT) 
            
      ELSEIF (IWAVE == 2) THEN  
      ! STOKE WAVE THEORY  ( PROVE 13-02-2013 )
      CALL PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55)    
      
      
      CALL STOKE_COEFFICIENT(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
                
C	-----------------------------------
      G1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
      G2 = 2.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
      G3 = 3.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
      G4 = 4.0D0*(AVAL**4.0)*G44
      G5 = 5.0D0*(AVAL**5.0)*G55
C	-----------------------------------
C	DRAG FORCE        
      UDXold =(G1*COSH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW)*COS(1.0D0*(RK*X
     1   - OMEGA*TT))
     1   + G2*COSH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW)*COS(2.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G3*COSH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW)*COS(3.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G4*COSH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW)*COS(4.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G5*COSH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW)*COS(5.0D0*(RK*X 
     1   - OMEGA*TT)))
     1    *(OMEGA/RK)
      
      UDY_old=(G1*SINH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW)*SIN(1.0D0*(RK*X 
     1   - OMEGA*TT)) 
     1   + G2*SINH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW)*SIN(2.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G3*SINH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW)*SIN(3.0D0*(RK*X 
     1   - OMEGA*TT)) 
     1   + G4*SINH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW)*SIN(4.0D0*(RK*X 
     1   - OMEGA*TT))
     1   + G5*SINH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW)*SIN(5.0D0*(RK*X 
     1   - OMEGA*TT)))
     1    *(OMEGA/RK)  
      
      !==================================
      !    MODIFY BY CHANA 23 FEB 2016
      !==================================
      
      WVHIGHT = H
      WDEPTH  = HW
      GRAV = G
      RK=(2*(AVAL+((AVAL**3)*F33)+((AVAL**5)*(F35+F55))))/WVHIGHT 
      
      CWAVE = SQRT(GRAV/RK*(1.0D0 + C1*(AVAL**2.0) + C2*AVAL**4.0)*TANH(RK*WDEPTH))
      
      chana = (OMEGA/RK)
      
      FG1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
      FG2 = 1.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
      FG3 = 1.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
      FG4 = 1.0D0*(AVAL**4.0)*G44
      FG5 = 1.0D0*(AVAL**5.0)*G55
      
      UDX = CWAVE*( 1.0D0*FG1*COSH(1.0D0*RK*Y)*COS(1.0D0*(RK*X - OMEGA*TT))
     1            + 2.0D0*FG2*COSH(2.0D0*RK*Y)*COS(2.0D0*(RK*X - OMEGA*TT))
     1            + 3.0D0*FG3*COSH(3.0D0*RK*Y)*COS(3.0D0*(RK*X - OMEGA*TT))
     1            + 4.0D0*FG4*COSH(4.0D0*RK*Y)*COS(4.0D0*(RK*X - OMEGA*TT))
     1            + 5.0D0*FG5*COSH(5.0D0*RK*Y)*COS(5.0D0*(RK*X - OMEGA*TT)))
      
      UDY = CWAVE*( 1.0D0*FG1*SINH(1.0D0*RK*Y)*SIN(1.0D0*(RK*X - OMEGA*TT))
     1            + 2.0D0*FG2*SINH(2.0D0*RK*Y)*SIN(2.0D0*(RK*X - OMEGA*TT))
     1            + 3.0D0*FG3*SINH(3.0D0*RK*Y)*SIN(3.0D0*(RK*X - OMEGA*TT))
     1            + 4.0D0*FG4*SINH(4.0D0*RK*Y)*SIN(4.0D0*(RK*X - OMEGA*TT))
     1            + 5.0D0*FG5*SINH(5.0D0*RK*Y)*SIN(5.0D0*(RK*X - OMEGA*TT)))
      
      

C	-----------------------------------
      ELSEIF (IWAVE == 3)THEN
      ! STREAM FUNCTION
      ! ONLY STREAM FUNCTION GENERATE VELOCITY AND ACCELERATION IN ONE SUBROTINE
        CALL VELOC (X,Y,UDX,UDY,DUDT,DVDTM,TT,RAMDA)
        AXSTREAM=DUDT
        AYSTREAM=DVDTM
      ELSEIF (IWAVE == 4)THEN
        CALL Cnoidal_wave_Length (H,HW,G,PERIOD,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
        RK=wavenumber
        OMEGA=2d0*AKKK/PERIOD
        ZZETAA = wavenumber * X - 2d0 * AKKK/PERIOD * TT
        CALL Jacobian_elliptic_function (ZZETAA,Ammm,CN,CN2,SSN,ddn)
        CALL Cnoidal_Wave_Kinematic  (H,HW,G,PERIOD,X,Y,wavenumber,AKKK,ALAMDA,Ammm,AEEE,CN,SSN,ddn,ZZETAA,UDX,
     1UDY,DUDT,DVDTM)
        AXSTREAM=DUDT
        AYSTREAM=DVDTM
        ELSEIF (IWAVE == 5)THEN
        call Solitary_Wave (H,HW,G,TT,X,Y,UDX,UDY,DUDT,DVDTM,Surface_Elevation,q)
        AXSTREAM=DUDT
        AYSTREAM=DVDTM
        
        ELSEIF (IWAVE == 7)THEN    
        CALL PMSpectrum_Time_domain      
        call Generate_Spectrum_Time( TT,RK,X,OMEGA  , SumWaveSpmTCos , SumWaveSpmTSin , SumWaveSjsTCos , SumWaveSjsTSin  )
        
         write(987655,*) TT , SumWaveSpmTCos , SumWaveSpmTSin
         
        UDX = 0.50D0*OMEGA*H*COSH(RK*Y)/SINH(RK*HW)*SumWaveSpmTCos  !COS(RK*X - OMEGA*TT) 
        UDY = 0.50D0*OMEGA*H*SINH(RK*Y)/SINH(RK*HW)*SumWaveSpmTSin  !SIN(RK*X - OMEGA*TT)     
        
       ELSEIF (IWAVE == 8)THEN    
        CALL PMSpectrum_Time_domain     
        call Generate_Spectrum_Time( TT,RK,X,OMEGA  , SumWaveSpmTCos , SumWaveSpmTSin , SumWaveSjsTCos , SumWaveSjsTSin  )
      
        write(987655,*) TT , SumWaveSjsTCos , SumWaveSjsTSin
        
      UDX = 0.50D0*OMEGA*H*COSH(RK*Y)/SINH(RK*HW)*SumWaveSjsTCos    !COS(RK*X - OMEGA*TT) 
      UDY = 0.50D0*OMEGA*H*SINH(RK*Y)/SINH(RK*HW)*SumWaveSjsTSin    !SIN(RK*X - OMEGA*TT) 
        
      ENDIF  
       
      ELSEIF (IORRE.EQ.2)THEN ! IRRGULA WAVE
      ! RANDOM WAVE  12-01-2021
      ! GENERATE WAVE VELOCITY AND ACCELERATION
      IF (IRWAVE.EQ.1)THEN
      CALL Pierson_Moskowitz_RANDOM ("DAI",TT,offselect,WNU,X,Y,UDX,UDY,AUX,AUY)  
      ELSEIF (IRWAVE.EQ.2) THEN
      CALL Jonswap_RANDOM ("DAI",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY) 
      ELSEIF (IRWAVE.EQ.3)THEN
      CALL Ochi_RANDOM ("DAI",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY) 
      ELSEIF (IRWAVE.EQ.4)THEN
      CALL Bretschneider_RANDOM ("DAI",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY) 
      ELSEIF (IRWAVE.EQ.5)THEN
      CALL TMA_SPectrum_RANDOM ("DAI",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY)     
      ENDIF
      
      ENDIF
                
      
C     -----------------------------------------------------------------------   
      RETURN
	END   
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE DRAG_ACCELERATION(OMEGA,RATIO,H,HW,RK,X,Y,TT,IWAVE,ORDER,AVAL,G,AAX,AAY)   
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------  

C	-----------------------------------------------------------------------
      If (IWAVE == 1) THEN    
      ! AIRY WAVE THEORY  
 
      AAX = + 0.50D0*(OMEGA*OMEGA)*H*COSH(RK*Y)/
     1      SINH(RK*HW)*SIN(RK*X - OMEGA*TT) 
      AAY = - 0.50D0*(OMEGA*OMEGA)*H*SINH(RK*Y)/
     1      SINH(RK*HW)*COS(RK*X - OMEGA*TT)   
          
      ELSEIF (IWAVE == 2) THEN
      ! STOKE FIFTH ORDER THEORY
      CALL PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55)
      CALL PARAMETER_C(RATIO,C1,C2,C3,C4)
      
      CALL STOKE_COEFFICIENT(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
C	--------------------------------------------
      G1 = 1.0D0*(AVAL*G11 + (AVAL**3)*G13 + (AVAL**5)*G15)
      G2 = 2.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
      G3 = 3.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
      G4 = 4.0D0*(AVAL**4.0)*G44
      G5 = 5.0D0*(AVAL**5.0)*G55           
C	--------------------------------------------
      U1 = G1*(COSH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW))
      U2 = G2*(COSH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW))
      U3 = G3*(COSH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW))
      U4 = G4*(COSH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW))
      U5 = G5*(COSH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW))
C	--------------------------------------------
      V1 = G1*(SINH(1.0D0*RK*Y)/SINH(1.0D0*RK*HW))
      V2 = G2*(SINH(2.0D0*RK*Y)/SINH(2.0D0*RK*HW))
      V3 = G3*(SINH(3.0D0*RK*Y)/SINH(3.0D0*RK*HW))
      V4 = G4*(SINH(4.0D0*RK*Y)/SINH(4.0D0*RK*HW))
      V5 = G5*(SINH(5.0D0*RK*Y)/SINH(5.0D0*RK*HW))
C	---------------------------------------------
      R1 =  2.0D0*U1 - 1.0D0*U1*U2 - 1.0D0*V1*V2 - 1.0D0*U2*U3 
     1      - 1.0D0*V2*V3
      R2 =  4.0D0*U2 - 1.0D0*U1*U1 + 1.0D0*V1*V1 - 2.0D0*U1*U3 
     1      - 2.0D0*V1*V3
      R3 =  6.0D0*U3 - 3.0D0*U1*U2 + 3.0D0*V1*V2 - 3.0D0*U1*U4 
     1      - 3.0D0*V1*V4
      R4 =  8.0D0*U4 - 2.0D0*U2*U2 + 2.0D0*V2*V2 - 4.0D0*U1*U3 
     1      + 4.0D0*V1*V3
      R5 = 10.0D0*U5 - 5.0D0*U1*U4 - 5.0D0*U2*U3 + 5.0D0*V1*V4 
     1      + 5.0D0*V2*V3
C	-----------------------------------
      S0 = -2.0D0*U1*V1 + 4.0D0*U2*V2
      S1 =  2.0D0*V1    - 3.0D0*U1*V2  - 3.0D0*U2*V1  - 5.0D0*U2*V3 
     1     -5.0D0*U3*V2
      S2 =  4.0D0*V2    - 4.0D0*U1*V3  - 4.0D0*U3*V1
      S3 =  6.0D0*V3    - 1.0D0*U1*V2  + 1.0D0*U2*V1  - 5.0D0*U1*V4 
     1     -5.0D0*U4*V1
      S4 =  8.0D0*V4    - 2.0D0*U1*V3  + 2.0D0*U3*V1  + 4.0D0*U2*V2
      S5 = 10.0D0*V5    - 3.0D0*U1*V4  + 3.0D0*U4*V1  - 1.0D0*U2*V3 
     1     +1.0D0*U3*V2
      C = SQRT(G/RK*(1.0D0 + (AVAL**2.0)*C1 + (AVAL**4.0)*C2)*TANH(RK*HW))
      AAX_OLD =  (0.50D0*RK*C**2.0)*(R1*SIN(1.0D0*(RK*X - OMEGA*TT)) +
     1                           R2*SIN(2.0D0*(RK*X - OMEGA*TT)) +
     1                           R3*SIN(3.0D0*(RK*X - OMEGA*TT)) +
     1                           R4*SIN(4.0D0*(RK*X - OMEGA*TT)) +
     1                           R5*SIN(5.0D0*(RK*X - OMEGA*TT)))       
      AAY_OLD = -(0.50D0*RK*C**2.0)*(S1*COS(1.0D0*(RK*X - OMEGA*TT)) +
     1                           S2*COS(2.0D0*(RK*X - OMEGA*TT)) +
     1                           S3*COS(3.0D0*(RK*X - OMEGA*TT)) +
     1                           S4*COS(4.0D0*(RK*X - OMEGA*TT)) +
     1                           S5*COS(5.0D0*(RK*X - OMEGA*TT)))
      
      
      !==================================
      !    MODIFY BY CHANA 23 FEB 2016
      !==================================
      
      WVHIGHT = H
      WDEPTH  = HW
      GRAV = G
      RK=(2*(AVAL+((AVAL**3)*F33)+((AVAL**5)*(F35+F55))))/WVHIGHT 
      
      CWAVE = SQRT(GRAV/RK*(1.0D0 + C1*(AVAL**2.0) + C2*AVAL**4.0)*TANH(RK*WDEPTH))
      
      chana = (OMEGA/RK)
      
      CALL STOKE_COEFFICIENT(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
      
      FG1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
      FG2 = 1.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
      FG3 = 1.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
      FG4 = 1.0D0*(AVAL**4.0)*G44
      FG5 = 1.0D0*(AVAL**5.0)*G55
      
      AAX = OMEGA*CWAVE*( (1.0D0**2)*FG1*COSH(1.0D0*RK*Y)*SIN(1.0D0*(RK*X - OMEGA*TT))
     1                     + (2.0D0**2)*FG2*COSH(2.0D0*RK*Y)*SIN(2.0D0*(RK*X - OMEGA*TT))
     1                     + (3.0D0**2)*FG3*COSH(3.0D0*RK*Y)*SIN(3.0D0*(RK*X - OMEGA*TT))
     1                     + (4.0D0**2)*FG4*COSH(4.0D0*RK*Y)*SIN(4.0D0*(RK*X - OMEGA*TT))
     1                     + (5.0D0**2)*FG5*COSH(5.0D0*RK*Y)*SIN(5.0D0*(RK*X - OMEGA*TT)))
      
      AAY = -OMEGA*CWAVE*( (1.0D0**2)*FG1*SINH(1.0D0*RK*Y)*COS(1.0D0*(RK*X - OMEGA*TT))
     1                      + (2.0D0**2)*FG2*SINH(2.0D0*RK*Y)*COS(2.0D0*(RK*X - OMEGA*TT))
     1                      + (3.0D0**2)*FG3*SINH(3.0D0*RK*Y)*COS(3.0D0*(RK*X - OMEGA*TT))
     1                      + (4.0D0**2)*FG4*SINH(4.0D0*RK*Y)*COS(4.0D0*(RK*X - OMEGA*TT))
     1                      + (5.0D0**2)*FG5*SINH(5.0D0*RK*Y)*COS(5.0D0*(RK*X - OMEGA*TT)))
      
      chana = 3
      
      ELSEIF(IWAVE == 3)THEN
      ! STREAM FUNCTION
      ! STREAM ACCELERATION GENERATE BY SUBROUTINE DRAG_VELOCITY
      ELSEIF(IWAVE == 7)THEN
      call Generate_Spectrum_Time( TT,RK,X,OMEGA  , SumWaveSpmTCos , SumWaveSpmTSin , SumWaveSjsTCos , SumWaveSjsTSin  )             
      AAX = + 0.50D0*(OMEGA*OMEGA)*H*COSH(RK*Y)/
     1      SINH(RK*HW)*SumWaveSpmTSin  !SIN(RK*X - OMEGA*TT) 
      AAY = - 0.50D0*(OMEGA*OMEGA)*H*SINH(RK*Y)/
     1      SINH(RK*HW)*SumWaveSpmTCos   !COS(RK*X - OMEGA*TT)        
      ELSEIF(IWAVE == 8)THEN      
      call Generate_Spectrum_Time( TT,RK,X,OMEGA  , SumWaveSpmTCos , SumWaveSpmTSin , SumWaveSjsTCos , SumWaveSjsTSin  )                 
      AAX = + 0.50D0*(OMEGA*OMEGA)*H*COSH(RK*Y)/
     1      SINH(RK*HW)*SumWaveSjsTCos  !SIN(RK*X - OMEGA*TT) 
      AAY = - 0.50D0*(OMEGA*OMEGA)*H*SINH(RK*Y)/
     1      SINH(RK*HW)*SumWaveSjsTSin  !COS(RK*X - OMEGA*TT)   
          
      ENDIF      
      
C     ------------------------------------------------------------------   
      RETURN
	END   
C	=======================================================================
C	=======================OFFSHORE LOAD PARAMETER=========================
C	=======================================================================   
      SUBROUTINE PARAMETER_A(RK,H,RATIO,AVAL)   
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------

      CALL PARAMETER_F(RATIO,F22,F24,F33,F35,F44,F55)
C     OMEGA
      R = RK*H/2.0D0

C     THE INITIAL VALUES
      AVAL = 0.10D0
      NITER = 100

C     ITERATION 
      DISP = AVAL
      
      DO II = 2,NITER
C     FUNCTION 
      Y = (F35 + F55)*(AVAL**5.0) + F33*(AVAL**3.0) + AVAL
C     DIRIVERTIVE FUNCTION 
      DY = 5.0D0*(F35 + F55)*(AVAL**4.0) + 3.0*F33*(AVAL**2.0) + 1.0D0
            
C     ---------------------------------------------
      F = R - Y
C     ---------------------------------------------      
      DELTA = F/DY
C     ---------------------------------------------
C     Update the initial value
      DISP = DISP + DELTA
      AVAL = DISP
            
C     Check all values---------------------------    
      YREAL = (F35 + F55)*(AVAL**5.0) + F33*(AVAL**3.0) + AVAL
         
      ERR = ABS((R - YREAL)/R*100.0D0)
      IF(ERR <= 0.0001) THEN
      EXIT 
      ENDIF    

      ENDDO      

      CHANA = 3

C     ------------------------------------------------------------------   
      RETURN
	END   
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE PARAMETER_F(RATIO,F22,F24,F33,F35,F44,F55)  
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------
      DIMENSION F(9,6),Y3(1,6)

      F(1:9,1:6) = 0.0D0

C	F(22)-------------------------------------------------------------------
      F(1,1) = 3.8920D0
      F(2,1) = 1.5390D0
      F(3,1) = 0.9270D0
      F(4,1) = 0.6990D0
      F(5,1) = 0.5990D0
      F(6,1) = 0.5510D0
      F(7,1) = 0.5270D0
      F(8,1) = 0.5070D0
      F(9,1) = 0.5020D0
C	F(24)-------------------------------------------------------------------
      F(1,2) = -28.610D0
      F(2,2) = 1.3440D0
      F(3,2) = 1.3980D0
      F(4,2) = 1.0640D0
      F(5,2) = 0.8930D0
      F(6,2) = 0.8040D0
      F(7,2) = 0.7590D0
      F(8,2) = 0.7220D0
      F(9,2) = 0.7120D0
C	F(33)-------------------------------------------------------------------
      F(1,3) = 13.090D0
      F(2,3) = 2.3810D0
      F(3,3) = 0.9960D0
      F(4,3) = 0.6300D0
      F(5,3) = 0.4950D0
      F(6,3) = 0.4350D0
      F(7,3) = 0.4100D0
      F(8,3) = 0.3840D0
      F(9,3) = 0.3770D0
C	F(35)-------------------------------------------------------------------
      F(1,4) = -138.60D0
      F(2,4) = 6.9350D0
      F(3,4) = 3.6790D0
      F(4,4) = 2.2440D0
      F(5,4) = 1.6850D0
      F(6,4) = 1.4380D0
      F(7,4) = 1.3300D0
      F(8,4) = 1.2300D0
      F(9,4) = 1.2050D0
C	F(44)-------------------------------------------------------------------
      F(1,5) = 44.990D0
      F(2,5) = 4.1470D0
      F(3,5) = 1.2590D0
      F(4,5) = 0.6760D0
      F(5,5) = 0.4840D0
      F(6,5) = 0.4070D0
      F(7,5) = 0.3710D0
      F(8,5) = 0.3440D0
      F(9,5) = 0.3370D0
C	F(55)-------------------------------------------------------------------
      F(1,6) = 163.80D0
      F(2,6) = 7.9350D0
      F(3,6) = 1.7340D0
      F(4,6) = 0.7970D0
      F(5,6) = 0.5250D0
      F(6,6) = 0.4200D0
      F(7,6) = 0.3730D0
      F(8,6) = 0.3390D0
      F(9,6) = 0.3290D0
C	------------------------------------------------------------------------
      Y3(1,1:6) = 0.0D0
      
      IF ((RATIO < 0.10)) THEN
!      WRITE (6,10) 
!   10 FORMAT (2X,'H/L < 0.10 ( NOT INCLUDE IN PROGRAMS )')
      ELSEIF ((RATIO >= 0.10).AND.(RATIO < 0.15)) THEN    
        X1 = 0.100D0
        X2 = 0.150D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(1,II)
        Y2 = F(2,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO       
      ELSEIF ((RATIO >= 0.15).AND.(RATIO < 0.20)) THEN
        X1 = 0.150D0
        X2 = 0.200D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(2,II)
        Y2 = F(3,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.20).AND.(RATIO < 0.25)) THEN
        X1 = 0.200D0
        X2 = 0.250D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(3,II)
        Y2 = F(4,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.25).AND.(RATIO < 0.30)) THEN
        X1 = 0.250D0
        X2 = 0.300D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(4,II)
        Y2 = F(5,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.30).AND.(RATIO < 0.35)) THEN
        X1 = 0.300D0
        X2 = 0.350D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(5,II)
        Y2 = F(6,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.35).AND.(RATIO < 0.40)) THEN
        X1 = 0.350D0
        X2 = 0.400D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(6,II)
        Y2 = F(7,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.40).AND.(RATIO < 0.50)) THEN
        X1 = 0.400D0
        X2 = 0.500D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(7,II)
        Y2 = F(8,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.50).AND.(RATIO <= 0.60)) THEN
        X1 = 0.500D0
        X2 = 0.600D0
        X3 = RATIO
        DO II = 1,6
        Y1 = F(8,II)
        Y2 = F(9,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO       
      ELSEIF (RATIO > 0.60D0)THEN 
!      WRITE (6,20)
!   20 FORMAT (2X,'H/L >0.6 ( NOT INCLUDE IN PROGRAMS )') 
      ENDIF  

C	------------------------------------------------------------------------

      F22 = Y3(1,1) 
      F24 = Y3(1,2) 
      F33 = Y3(1,3) 
      F35 = Y3(1,4) 
      F44 = Y3(1,5) 
      F55 = Y3(1,6) 
      
C	------------------------------------------------------------------------ 
      RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------
      DIMENSION G(9,9),Y3(1,9)
            
      G(1,1:9) = 0.0D0

C	G(11)-------------------------------------------------------------------
      G(1,1) = 1.0D0
      G(2,1) = 1.0D0
      G(3,1) = 1.0D0
      G(4,1) = 1.0D0
      G(5,1) = 1.0D0
      G(6,1) = 1.0D0
      G(7,1) = 1.0D0
      G(8,1) = 1.0D0
      G(9,1) = 1.0D0
C	G(13)-------------------------------------------------------------------
      G(1,2) = -7.3940D0
      G(2,2) = -2.3200D0
      G(3,2) = -1.2630D0
      G(4,2) = -0.9110D0
      G(5,2) = -0.7650D0
      G(6,2) = -0.6960D0
      G(7,2) = -0.6620D0
      G(8,2) = -0.6350D0
      G(9,2) = -0.6280D0  
         
C	G(15)-------------------------------------------------------------------
      G(1,3) = -12.730D0
      G(2,3) = -4.8640D0
      G(3,3) = -2.2660D0
      G(4,3) = -1.4150D0
      G(5,3) = -1.0770D0
      G(6,3) = -0.9250D0
      G(7,3) = -0.8500D0
      G(8,3) = -0.7900D0
      G(9,3) = -0.7770D0   
  
C	G(22)-------------------------------------------------------------------
      G(1,4) = 2.9960D0
      G(2,4) = 0.8600D0
      G(3,4) = 0.3260D0
      G(4,4) = 0.1540D0
      G(5,4) = 0.0760D0
      G(6,4) = 0.0380D0
      G(7,4) = 0.0200D0
      G(8,4) = 0.0060D0
      G(9,4) = 0.0020D0   
  
C	G(24)-------------------------------------------------------------------
      G(1,5) = -48.140D0
      G(2,5) = -0.9070D0
      G(3,5) = 0.6800D0
      G(4,5) = 0.6730D0
      G(5,5) = 0.6010D0
      G(6,5) = 0.5560D0
      G(7,5) = 0.5280D0
      G(8,5) = 0.5030D0
      G(9,5) = 0.5020D0
   
C	G(33)-------------------------------------------------------------------
      G(1,6) = 5.9420D0
      G(2,6) = 0.3100D0
      G(3,6) = -0.0170D0
      G(4,6) = -0.0300D0
      G(5,6) = -0.0200D0
      G(6,6) = -0.0120D0
      G(7,6) = -0.0060D0
      G(8,6) = -0.0020D0
      G(9,6) = -0.0010D0   
  
C	G(35)-------------------------------------------------------------------
      G(1,7) = -121.70D0
      G(2,7) = 2.84300D0
      G(3,7) = 1.09300D0
      G(4,7) = 0.44000D0
      G(5,7) = 0.23100D0
      G(6,7) = 0.15200D0
      G(7,7) = 0.11700D0
      G(8,7) = 0.09200D0
      G(9,7) = 0.08600D0
C	G(44)-------------------------------------------------------------------
      G(1,8) = 7.67100D0
      G(2,8) = -0.16700D0
      G(3,8) = -0.04400D0
      G(4,8) = -0.00500D0
      G(5,8) = 0.00200D0
      G(6,8) = 0.00200D0
      G(7,8) = 0.00100D0
      G(8,8) = 0.00000D0
      G(9,8) = 0.00000D0
C	G(55)-------------------------------------------------------------------
      G(1,9) = 0.89200D0
      G(2,9) = -0.2570D0
      G(3,9) = 0.0060D0
      G(4,9) = 0.0050D0
      G(5,9) = 0.0010D0
      G(6,9) = 0.0000D0
      G(7,9) = 0.0000D0
      G(8,9) = 0.0000D0
      G(9,9) = 0.0000D0
C	------------------------------------------------------------------------

      Y3(1,1:9) = 0.0D0
      
      IF (RATIO < 0.10) THEN
!      WRITE (6,30) 
!   30 FORMAT (2X,'H/L < 0.10 ( NOT INCLUDE IN PROGRAMS )')
      ELSEIF ((RATIO >= 0.10).AND.(RATIO < 0.15)) THEN    
        X1 = 0.100D0
        X2 = 0.150D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(1,II)
        Y2 = G(2,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO       
      ELSEIF ((RATIO >= 0.15).AND.(RATIO < 0.20)) THEN
        X1 = 0.150D0
        X2 = 0.200D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(2,II)
        Y2 = G(3,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO
      ELSEIF ((RATIO >= 0.20).AND.(RATIO < 0.25)) THEN
        X1 = 0.200D0
        X2 = 0.250D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(3,II)
        Y2 = G(4,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO  
      ELSEIF ((RATIO >= 0.25).AND.(RATIO < 0.30)) THEN
        X1 = 0.250D0
        X2 = 0.300D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(4,II)
        Y2 = G(5,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.30).AND.(RATIO < 0.35)) THEN
        X1 = 0.300D0
        X2 = 0.350D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(5,II)
        Y2 = G(6,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.35).AND.(RATIO < 0.40)) THEN
        X1 = 0.350D0
        X2 = 0.400D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(6,II)
        Y2 = G(7,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.40).AND.(RATIO < 0.50)) THEN
        X1 = 0.400D0
        X2 = 0.500D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(7,II)
        Y2 = G(8,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.50).AND.(RATIO <= 0.60)) THEN
        X1 = 0.500D0
        X2 = 0.600D0
        X3 = RATIO
        DO II = 1,9
        Y1 = G(8,II)
        Y2 = G(9,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO       
      ELSEIF (RATIO > 0.60D0)THEN 
!      WRITE (6,40)
!   40 FORMAT (2X,'H/L >0.6 ( NOT INCLUDE IN PROGRAMS )') 
      ENDIF  

C	------------------------------------------------------------------------
      G11 = Y3(1,1) 
      G13 = Y3(1,2) 
      G15 = Y3(1,3) 
      G22 = Y3(1,4)
      G24 = Y3(1,5) 
      G33 = Y3(1,6) 
      G35 = Y3(1,7) 
      G44 = Y3(1,8) 
      G55 = Y3(1,9)   
      
      
C	------------------------------------------------------------------------ 
      RETURN
	END	
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE PARAMETER_C(RATIO,C1,C2,C3,C4)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------
      DIMENSION C(9,4),Y3(1,6)
            
      
      C(9,4) = 0.0D0

C	C(11)-------------------------------------------------------------------
      C(1,1) = 8.7910D0
      C(2,1) = 2.6460D0
      C(3,1) = 1.5490D0
      C(4,1) = 1.2290D0
      C(5,1) = 1.1070D0
      C(6,1) = 1.0550D0
      C(7,1) = 1.0270D0
      C(8,1) = 1.0080D0
      C(9,1) = 1.0020D0
  
C	C(13)-------------------------------------------------------------------
      C(1,2) = 383.70D0
      C(2,2) = 19.820D0
      C(3,2) = 5.0440D0
      C(4,2) = 2.5680D0
      C(5,2) = 1.8330D0
      C(6,2) = 1.5320D0
      C(7,2) = 1.3930D0
      C(8,2) = 1.2830D0
      C(9,2) = 1.2400D0   
           
C	C(15)-------------------------------------------------------------------
      C(1,3) = -0.3100D0
      C(2,3) = -0.1550D0
      C(3,3) = -0.0820D0
      C(4,3) = -0.0430D0
      C(5,3) = -0.0230D0
      C(6,3) = -0.0120D0
      C(7,3) = -0.0070D0
      C(8,3) = -0.0010D0
      C(9,3) = -0.0010D0   
     
C	C(22)-------------------------------------------------------------------
      C(1,4) = -0.0600D0
      C(2,4) = 0.2570D0
      C(3,4) = 0.0770D0
      C(4,4) = 0.0280D0
      C(5,4) = 0.0100D0
      C(6,4) = 0.0040D0
      C(7,4) = 0.0020D0
      C(8,4) = 0.0000D0
      C(9,4) = 0.0000D0 

C	------------------------------------------------------------------------

      Y3(1,1:4) = 0.0D0
      
      IF (RATIO < 0.10) THEN
!      WRITE (6,50) 
!   50 FORMAT (2X,'H/L < 0.10 ( NOT INCLUDE IN PROGRAMS )')
      ELSEIF ((RATIO >= 0.10).AND.(RATIO < 0.15)) THEN    
        X1 = 0.100D0
        X2 = 0.150D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(1,II)
        Y2 = C(2,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO       
      ELSEIF ((RATIO >= 0.15).AND.(RATIO < 0.20)) THEN
        X1 = 0.150D0
        X2 = 0.200D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(2,II)
        Y2 = C(3,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO
      ELSEIF ((RATIO >= 0.20).AND.(RATIO < 0.25)) THEN
        X1 = 0.200D0
        X2 = 0.250D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(3,II)
        Y2 = C(4,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.25).AND.(RATIO < 0.30)) THEN
        X1 = 0.250D0
        X2 = 0.300D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(4,II)
        Y2 = C(5,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.30).AND.(RATIO < 0.35)) THEN
        X1 = 0.300D0
        X2 = 0.350D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(5,II)
        Y2 = C(6,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.35).AND.(RATIO < 0.40)) THEN
        X1 = 0.350D0
        X2 = 0.400D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(6,II)
        Y2 = C(7,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO   
      ELSEIF ((RATIO >= 0.40).AND.(RATIO < 0.50)) THEN
        X1 = 0.400D0
        X2 = 0.500D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(7,II)
        Y2 = C(8,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO  
      ELSEIF ((RATIO >= 0.50).AND.(RATIO <= 0.60)) THEN
        X1 = 0.50D0
        X2 = 0.60D0
        X3 = RATIO
        DO II = 1,4
        Y1 = C(8,II)
        Y2 = C(9,II)
        Y3(1,II) = Y2 - (X2 - X3)*(Y2 - Y1)/(X2 - X1)
        ENDDO       
      ELSEIF (RATIO > 0.60D0)THEN 
!      WRITE (6,60)
!   60 FORMAT (2X,'H/L >0.6 ( NOT INCLUDE IN PROGRAMS )') 
      ENDIF  

C	------------------------------------------------------------------------

      C1 = Y3(1,1) 
      C2 = Y3(1,2) 
      C3 = Y3(1,3) 
      C4 = Y3(1,4) 

C	------------------------------------------------------------------------ 
      RETURN
	END
C	=======================================================================
C	=======================================================================
C	======================================================================= 
      SUBROUTINE GUASSDATA(NGP,GP,GW)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------
      DIMENSION GP(NGP),GW(NGP)


      IF (NGP == 1) THEN      
      GP(1) = 0.00D0
      GW(1) = 2.00D0
    
      ELSEIF (NGP == 2) THEN
      GP(1) = -SQRT(1.0/3.0)
      GP(2) =  SQRT(1.0/3.0)
    
      GW(1) = 1.0D0
      GW(2) = 1.0D0
    
      ELSEIF (NGP == 3) THEN
      GP(1) = -SQRT(0.6)
      GP(2) =  0.0D0
      GP(3) =  SQRT(0.6)
    
      GW(1) = 5.0D0/9.0
      GW(2) = 8.0D0/9.0
      GW(3) = 5.0D0/9.0
    
      ELSEIF (NGP == 4) THEN
      GP(1) = -0.8611363120D0
      GP(2) = -0.3399810440D0
      GP(3) =  0.3399810440D0
      GP(4) =  0.8611363120D0
    
      GW(1) =  0.34785480D0
      GW(2) =  0.65214520D0
      GW(3) =  0.65214520D0
      GW(4) =  0.34785480D0
    
      ELSEIF (NGP == 5) THEN
      GP(1) = -0.9061798460D0
      GP(2) = -0.5384693100D0
      GP(3) =  0.0
      GP(4) =  0.5384693100D0
      GP(5) =  0.9061798460D0
    
      GW(1) =  0.23692690D0
      GW(2) =  0.47862870D0
      GW(3) =  0.56888890D0   
      GW(4) =  0.47862870D0
      GW(5) =  0.23692690D0
    
      ELSEIF (NGP == 6) THEN
      GP(1) = -0.9324695140D0
      GP(2) = -0.6612093860D0
      GP(3) = -0.2386191860D0
      GP(4) =  0.2386191860D0
      GP(5) =  0.6612093860D0
      GP(6) =  0.9324695140D0
    
      GW(1) =  0.17132450D0
      GW(2) =  0.36076160D0
      GW(3) =  0.46791390D0    
      GW(4) =  0.46791390D0
      GW(5) =  0.36076160D0
      GW(6) =  0.17132450D0
    
      ELSEIF (NGP == 7) THEN     
      GP(1) = -0.949107910D0
      GP(2) = -0.741531190D0
      GP(3) = -0.405845150D0
      GP(4) =  0.00D0
      GP(5) =  0.405845150D0
      GP(6) =  0.741531190D0
      GP(7) =  0.949107910D0

      GW(1) =  0.129484970D0
      GW(2) =  0.279705390D0
      GW(3) =  0.381830050D0
      GW(4) =  0.417959180D0
      GW(5) =  0.381830050D0
      GW(6) =  0.279705390D0
      GW(7) =  0.129484970D0 
    
      ELSEIF (NGP == 8) THEN
      GP(1) = -0.960289860D0
      GP(2) = -0.796666480D0
      GP(3) = -0.525532410D0
      GP(4) = -0.183434640D0
      GP(5) =  0.183434640D0
      GP(6) =  0.525532410D0
      GP(7) =  0.796666480D0
      GP(8) =  0.960289860D0
    
      GW(1) =  0.101228540D0
      GW(2) =  0.222381030D0
      GW(3) =  0.313706650D0
      GW(4) =  0.362683780D0
      GW(5) =  0.362683780D0
      GW(6) =  0.313706650D0
      GW(7) =  0.222381030D0
      GW(8) =  0.101228540D0  
    
      ELSEIF (NGP == 9) THEN    
      GP(1) = -0.9681602395076261
      GP(2) = -0.8360311073266358
      GP(3) = -0.6133714327005904
      GP(4) = -0.3242534234038089
      GP(5) =  0.0D0
      GP(6) =  0.324253423403808
      GP(7) =  0.613371432700590
      GP(8) =  0.836031107326635
      GP(9) =  0.968160239507626

      GW(1) =  0.0812743883615744
      GW(2) =  0.1806481606948570
      GW(3) =  0.2606106964029350
      GW(4) =  0.3123470770400020
      GW(5) =  0.3302393550012590
      GW(6) =  0.3123470770400020
      GW(7) =  0.2606106964029350
      GW(8) =  0.1806481606948570
      GW(9) =  0.0812743883615744

      ELSEIF (NGP == 10) THEN       
      GP(1)  = -0.973906530D0
      GP(2)  = -0.865063370D0
      GP(3)  = -0.679409570D0
      GP(4)  = -0.433395390D0
      GP(5)  = -0.148874340D0
      GP(6)  =  0.148874340D0
      GP(7)  =  0.433395390D0
      GP(8)  =  0.679409570D0
      GP(9)  =  0.865063370D0
      GP(10) =  0.973906530D0

      GW(1)  =  0.066671340D0
      GW(2)  =  0.149451350D0
      GW(3)  =  0.219086360D0
      GW(4)  =  0.269266720D0
      GW(5)  =  0.295524220D0
      GW(6)  =  0.295524220D0
      GW(7)  =  0.269266720D0
      GW(8)  =  0.219086360D0
      GW(9)  =  0.149451350D0
      GW(10) =  0.066671340D0    
      ENDIF      


C	------------------------------------------------------------------------ 
      RETURN
	END
C	=======================================================================
C	=======================OFFSHORE LOAD PARAMETER=========================
C	=======================================================================   
      SUBROUTINE NEWTON_RAPHSON(HW,G,OMEGA,RK)   
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------
      
C	OMEGA
      R = OMEGA**2.0

C	HE INITIAL VALUES
      RK = 1.0D0  
      NITER = 100   

C	ITERATION 
      DISP = RK
      ERR = 100.0D0

      DO II = 2,NITER  
C	FUNCTION 
      Y = G*RK*TANH(RK*HW)
C	DIRIVERTIVE FUNCTION 
      DY = G*TANH(RK*HW) + G*RK*(1.0 - TANH(RK*HW)**2.0)*HW
    
C	---------------------------------------------
      F = R-Y
C	---------------------------------------------
      DELTA = F/DY
C	---------------------------------------------
C	Update the initial value
      DISP = DISP + DELTA
      RK = DISP
    
C	Check all values---------------------------    
      YREAL = G*RK*TANH(RK*HW)   
    
      ERR = ABS((R - YREAL)/R*100.0D0)
      IF (ERR <= 0.0001) THEN
        EXIT 
      ENDIF    

      ENDDO


C     ------------------------------------------------------------------   
      RETURN
	END   
C	=======================================================================

C	=======================================================================
C	=======================================================================	
      SUBROUTINE INCL_WATERSURFACE_STOKE(HW,AVAL,RATIO,OMEGA,TT,RK,PXYZ,XL) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------
      DIMENSION PXYZ(2,3) 

      X1 = PXYZ(1,1)
      Y1 = PXYZ(1,2)
      X2 = PXYZ(2,1)      
      Y2 = PXYZ(2,2)

      CALL PARAMETER_F(RATIO,F22,F24,F33,F35,F44,F55)
      
      F1 = AVAL
      F2 = (AVAL**2)*F22 + (AVAL**4)*F24
      F3 = (AVAL**3)*F33 + (AVAL**5)*F35
      F4 = (AVAL**4)*F44
      F5 = (AVAL**5)*F55
C	OMEGA
      R = Y1 - HW

C	THE INITIAL VALUES
      XL = 1.0D0
      NITER = 100

C	ITERATION
      DISP = XL
      ERR = 100.0D0

      DO II = 2,NITER   
C	FUNCTION 
      AA = (Y2 - Y1)/(X2 - X1)
      Y = -AA*(XL - X1) + (1.0D0/RK)*(F1*COS(1.0D0*(RK*XL - OMEGA*TT)) +
     1                                F2*COS(2.0D0*(RK*XL - OMEGA*TT)) +
     1                                F3*COS(3.0D0*(RK*XL - OMEGA*TT)) +
     1                                F4*COS(4.0D0*(RK*XL - OMEGA*TT)) +
     1                                F5*COS(5.0D0*(RK*XL - OMEGA*TT)))
C	DIRIVERTIVE FUNCTION 
      DY = -AA - (1.0D0*F1*SIN(1.0D0*(RK*XL - OMEGA*TT)) + 
     1            2.0D0*F2*SIN(2.0D0*(RK*XL - OMEGA*TT)) + 
     1            3.0D0*F3*SIN(3.0D0*(RK*XL - OMEGA*TT)) + 
     1            4.0D0*F4*SIN(4.0D0*(RK*XL - OMEGA*TT)) + 
     1            5.0D0*F5*SIN(5.0D0*(RK*XL - OMEGA*TT)))
    
C	---------------------------------------------
      F = R-Y
C	---------------------------------------------
      DELTA = F/DY
C	---------------------------------------------
C	Update the initial value
      DISP = DISP + DELTA
      XL = DISP
    
C	Check all values---------------------------    
      YREAL = -AA*(XL - X1) + 
     1        (1.0D0/RK)*(F1*COS(1.0D0*(RK*XL - OMEGA*TT)) +
     1                    F2*COS(2.0D0*(RK*XL - OMEGA*TT)) + 
     1                    F3*COS(3.0D0*(RK*XL - OMEGA*TT)) + 
     1                    F4*COS(4.0D0*(RK*XL - OMEGA*TT)) + 
     1                    F5*COS(5.0D0*(RK*XL - OMEGA*TT)))
    
      ERR = ABS((R - YREAL)/R*100.0D0)
      IF (ERR <= 0.0001) THEN
        EXIT 
      ENDIF    

      ENDDO
      

C     ------------------------------------------------------------------   
      RETURN
	END  


C	===========================================================================
C	===========================================================================
C	===========================================================================	
      SUBROUTINE VECTORDOT(AVAL,B,UA,UB,DOTV)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	---------------------------------------------------------------------------     
      DIMENSION UA(1,3),UB(1,3),AVAL(1,3),B(1,3) 
      
C	MAGNITUDE OF VECTORS
      DA = SQRT(AVAL(1,1)*AVAL(1,1) + AVAL(1,2)*AVAL(1,2) + AVAL(1,3)*AVAL(1,3))
      DB = SQRT(B(1,1)*B(1,1) + B(1,2)*B(1,2) + B(1,3)*B(1,3))

C	UNIT VECTORS
      UA(1,1:3) = 0.0D0
      UB(1,1:3) = 0.0D0

      UA(1,1) = AVAL(1,1)/DA
      UA(1,2) = AVAL(1,2)/DA
      UA(1,3) = AVAL(1,3)/DA

      UB(1,1) = B(1,1)/DB
      UB(1,2) = B(1,2)/DB
      UB(1,3) = B(1,3)/DB

C	DOT PRODUCT
      DOTV = AVAL(1,1)*B(1,1) + AVAL(1,2)*B(1,2) + AVAL(1,3)*B(1,3)

C	ANGLE
      ANG = ACOS(DOTV/(DA*DB))      
      
      
C     -----------------------------------------------------------------------   
      RETURN
	END  
C	=======================================================================
C	=======================================================================
C	=======================================================================	
      SUBROUTINE VECTORUNIT(AVAL,UA,DA)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------  
      DIMENSION UA(1,3),AVAL(1,3)     
        
C	MAGNITUDE OF VECTORS
      DA = SQRT(AVAL(1,1)*AVAL(1,1) + AVAL(1,2)*AVAL(1,2) + AVAL(1,3)*AVAL(1,3))

C	UNIT VECTORS
      UA(1,1) = AVAL(1,1)/DA
      UA(1,2) = AVAL(1,2)/DA
      UA(1,3) = AVAL(1,3)/DA    
      
C     -----------------------------------------------------------------------   
      RETURN
	END   
C	=======================================================================
C	=======================================================================
C	=======================================================================	
      SUBROUTINE VECTORCROSS(AVAL,B,C)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)      
C	-----------------------------------------------------------------------  
      DIMENSION AVAL(1,3),B(1,3),C(1,3)     

      A1 = AVAL(1,1)
      A2 = AVAL(1,2)
      A3 = AVAL(1,3)

      B1 = B(1,1)
      B2 = B(1,2)
      B3 = B(1,3)

      C1 = A2*B3 - B2*A3
      C2 = A3*B1 - B3*A1
      C3 = A1*B2 - B1*A2

      C(1,1) = C1
      C(1,2) = C2
      C(1,3) = C3
      
      
      
      
      
C     -----------------------------------------------------------------------------------   
      RETURN
	END  
C	===================================================================================
C	===================================================================================
C	===================================================================================  

C	=======================================================================      
      SUBROUTINE WAVEBREAKLING
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)

C     ----------- INPUT ----------- 
C     DIAM           = NORMINAL DIAMETER
C     FCURLING       = CURLING FACTOR
C     UDX            = WAVE VELOCITY
C     GRAMMA         = THE ANGLE NETWEEN DIRECTION OF MOTION OF THE MASS OF WATER  
C     T1             = TIME            
     
      
      
      
      
C     -------------------------------------------------------
C                        WAVE BREAKLING 
C     -------------------------------------------------------    
      PI       = 3.141592653589793D0
      RADIUS   = DIAM/2D0
      FCURLING = 0.5D0
      
      A1     = RADIUS/(8D0*UDX*COS(GRAMMA))
      A2     = 12D0*RADIUS/(32D0*UDX*COS(GRAMMA))
      A3     = 3D0*RADIUS/(32D0*UDX*COS(GRAMMA))
c      T2     = T1-(RADIUS/(32D0*UDX*COS(GRAMMA))
      
      IF (T1.GT.0D0.AND.T1.LT.A1)THEN
      FTEST1 = (TANH(SQRT(1D0-(T1*UDX*COS(GRAMMA)/(4D0*RADIUS))))**-1D0)
      FTEST2 = 2*PI-(2*SQRT(T1*UDX*COS(GRAMMA)/RADIUS))
      F1     = FCURLING*HB*ROW*RADIUS*UDX*UDX*(COS(GRAMMA*(FTEST2*FTEST1))**2D0)
      ENDIF
      IF (T1.GE.A3.AND.T2.LE.A2)THEN
      FTEST1 = (TANH(SQRT(1D0-(T1*UDX*COS(GRAMMA)/(4D0*RADIUS))*(SQRT(6D0*UDX*COS(GRAMMA)*T2/RADIUS))))**-1D0)
      FTEST2 = ((8D0/3D0)*UDX*COS(GRAMMA)*T2)**(1D0/4D0)
      FTEST3 = PI*SQRT(1D0/(6D0*UDX*COS(GRAMMA)*T2/(RADIUS)))
      F1     = FCURLING*HB*ROW*RADIUS*UDX*UDX*(COS(FTEST3-FTEST2*FTEST1)**2D0)
      ENDIF
      
      
      END
c     =================================================================================
c     =================================================================================
      SUBROUTINE WAVE_SURFACE_PROFILE_FREQ 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION VERO(10),ACC(10)
      ! THIS SUBOUTINE IS USED FOR CALCULATE WAVE SURFACE PROFILE WHICH RELATED WITH SPECTRUM
      
      ! INPUT DATA 
      NTFREQ = TIME_INCREA/TOTAL_TIME
      
      DO I = 1,NTFREQ
      FREQ = 1/TIME_INCREA 
      ! WAVE SPECTRUM
      IF (IRWAVE.EQ.1)THEN     ! PIERSON MOSKOWITZ
      CALL Pierson_Moskowitz_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,WH1,WGV,SNA,SDG,TAP,UCURRENT)
      ELSEIF (IRWAVE.EQ.2)THEN ! JONSWAP 
      CALL Jonswap_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,WH1,WGV,SNA,SDG,TAP,UCURRENT)
      ELSEIF (IRWAVE.EQ.3)THEN ! OCHI 
      CALL Ochi_SPectrum (WVH1,PER1,SHAPE1,WVH2,PER2,SHAPE2,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,WH1,WGV,SNA,SDG,TAP,UCURRENT)
      ELSEIF (IRWAVE.EQ.4)THEN ! BRSTSCHNEIDER
      CALL  Bretschneider_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,WH1,WGV,SNA,SDG,TAP,UCURRENT)
      ELSEIF (IRWAVE.EQ.5)THEN ! TMA
      CALL TMA_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,WH1,WGV,SNA,SDG,TAP,UCURRENT)
      ELSEIF (IRWAVE.EQ.6)THEN ! USER DEFIND
      CALL User_Define_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,WH1,WGV,SNA,SDG,TAP,UCURRENT)
      ENDIF
      
      ! ZERO UPCORSSING
      OMEGA = SQRT(GRAV*AK*TANH(AK*WVHIGHT))
      AJ    = SQRT(SPM*0.5)
      SURFACE_PROFILE = SQRT(AJ)*COS(OMEGA)
      
      ! CALCULATE THE FORCE BY USING AIRY WAVE
C      VERO(I) = 0.5D0*OMEGA*WVHIGHT*(COSH(AK*Y/SINH(AK*WVHIGHT))*COS(AK*X-OMEGA*TIME)
C      ACC(I)  = 0.5D0*OMEGA*WVHIGHT*(SINH(AK*Y/SINH(AK*WVHIGHT))*SIN(AK*X-OMEGA*TIME)
      
      ! PRINTING THE RESULT OF SURFACE PROFILE
      
      
      
      ENDDO

      
      
      
      END
c     =================================================================================
c     =================================================================================