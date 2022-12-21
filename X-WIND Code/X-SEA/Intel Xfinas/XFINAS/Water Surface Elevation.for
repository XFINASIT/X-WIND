      Subroutine  PRINT_WATER_SURFACE_ELEVATION (XYZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION XYZ(1)
      DIMENSION AMAX_MIN(6)
      CHARACTER*200 NAME_PARAMETER,NAME_SPECTRUM,NAME
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
      COMMON /SEA_PARA/ NSPECTRUM
      COMMON /NAME_OFFSHORE_PARAMETER/ NAME_PARAMETER(1000),NAME_SPECTRUM(1000)
      COMMON /RANDOM_NUMBER_STORE/ RADNUM(1000,100)
      COMMON /MGRAV/ NGRAV  
      
      PI = 3.141592654d0
      
      AMAX_MIN = 0.D0
      ! CALCULATE MID OF STRUCTURE AT WATER SURFACE.
      CALL MID_STRUCTURE (XYZ,XX,YY,ZZ,AMAX_MIN) 
      
      ! LINEAR STATIC AND LINEAR DYNNAMIC CONTROL
      IF (ISOLOP.EQ.1.OR.ISOLOP.EQ.5) THEN
      
      DO I =1,NUM_OF_OFFSHORE_PARAMETER
      offselect = I    
      CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,VB,HM,HW,HC,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                  VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,
     1                  WKF,CBF,WFC,ATIME)
      
      ! SELECT CASE IWAVE ( 1:AIRY WAVE THEORY 2:STOKE FIFTH ORDER THEORY   
      IF (IWAVE == 1 .or. IWAVE == 7 .or. IWAVE == 8 ) THEN
          ! --- AIRY WAVE THEORY ---
          ! --- INPUT DATA ---
          ! WDEPTH  = WATER DEPTH
          ! GRAV    = GAVITY ACCLERATION
          ! OMEGA   = ANGULAR ACCELERATION
          ! --- OUTPUT DATA ---
          ! RK      = WAVE NUMBER      
          OMEGA = 2.0*PI/PERIOD      
          CALL Airy_Wave_Number(GRAV,OMEGA,WDEPTH,RK)               
          Check_Wave_Number = (OMEGA**2.0d0)/(GRAV*tanh(RK*WDEPTH))
          
          RAMDA = 2.0D0*PI/RK    
          AVAL = 1.0D0
          RATIO = 1.0D0   ! SET CONDITON     
      ELSEIF (IWAVE == 2) THEN
          ! --- STOKE FIFTH ORDER THEORY ---
          ! --- INPUT DATA ---
          ! GRAV,WDEPTH,PERIOD,WVHIGHT,RK,RAMDA,
          ! --- OUTPUT DATA ---
          ! ADUM    
             
          !!- CALL STOKES_WAVELENGTH(GRAV,WDEPTH,PERIOD,WVHIGHT,RK,RAMDA)
          CALL STOKES_WAVELENGTH_Time_Modify(AVAL,RAMDA) ! RK = ASSUME EQUAL, AVAL
          RATIO = WDEPTH/RAMDA
          
         
          CALL STOKE_COEFFICIENT(RATIO,A11,A13,A15,A22,A24,A33,A35,A44,A55,B22,B24,B33,B35,B44,B55,
     1                             C1X,C2X,C3X,C4X)

          RKx=2.0d0*22.0/7.0d0/RAMDA
          ! --- PARAMETER_F ---
            ! --- INPUT DATA ---
            ! RATIO     = DIMETER/WAVE LENGTH
            ! --- OUTPUT DATA ---
            ! F22,F24,F33,F35,F44,F55
          CALL PARAMETER_F(RATIO,F22,F24,F33,F35,F44,F55) ! HAVE STOP FUNCTION 
          !	GET ALL PARAMETERS
            ! --- PARAMETER_G ---
            ! --- INPUT DATA ---
            ! RATIO     = DIMETER/WAVE LENGTH
            ! --- OUTPUT DATA ---
            ! G11,G13,G15,G22,G24,G33,G35,G44,G55    
          CALL PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55)
          !	FREE-SURFACE WATER DEFLECTION (NU)
            ! --- INPUT DATA ---
            ! RATIO   = DIMETER/WAVE LENGTH
            ! RK,WV,RATIO
            ! --- OUTPUT DATA ---
            ! AVAL    = WAVE-HEIGHT PARAMETER
          CALL PARAMETER_C(RATIO,C1,C2,C3,C4)      
          
 !         CALL STOKE_COEFFICIENT(RATIO,A11,A13,A15,A22,A24,A33,A35,A44,A55,F22,F24,F33,F35,F44,F55,
 !    1                             C1,C2,C3,C4)
          
          
          CALL STOKE_COEFFICIENT(RATIO,A11,A13,A15,A22,A24,A33,A35,A44,A55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
          
          ! --- MODIFINE ON 22-08-2012 ----
          RK=(2*(AVAL+((AVAL**3)*F33)+((AVAL**5)*(F35+F55))))/WVHIGHT 
          ! -----------------------          
          
          CALL STOKE_COEFFICIENT(RATIO,A11,A13,A15,A22,A24,A33,A35,A44,A55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
                    
          
          !	---------------------------------
          !	FREE-SURFACE WATER DEFLECTION (NU)
          !CALL PARAMETER_A(RK,WVHIGHT,RATIO,AVAL)     
          OMEGA = SQRT(GRAV*RK*(1.0D0 + C1*(AVAL**2.0) + C2*AVAL**4.0)*TANH(RK*WDEPTH))
          
          
          CALL PARAMETER_F(RATIO,F22,F24,F33,F35,F44,F55) ! HAVE STOP FUNCTION 
          
          CALL STOKE_COEFFICIENT(RATIO,A11,A13,A15,A22,A24,A33,A35,A44,A55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
                    
          !	---------------------------------
          F1 = AVAL
          F2 = (AVAL**2.0)*F22 + (AVAL**4.0)*F24
          F3 = (AVAL**3.0)*F33 + (AVAL**5.0)*F35
          F4 = (AVAL**4.0)*F44
          F5 = (AVAL**5.0)*F55
          !	---------------------------------   
          
      ELSEIF (IWAVE == 3) THEN            
      ! FOR STEAMFUNCTION DON'T HAVE WAVE LENTGTH AND WAVE NUMBER
      ! NO NEED
      ENDIF    
      

           ! INITIAL CONDITION SET
           H1M  = 0.D0
           IF (NGRAV.EQ.1)THEN
               H1M = ZZ
           ELSEIF (NGRAV.EQ.2)THEN
               H1M = XX
           ELSEIF (NGRAV.EQ.3)THEN
               H1M = XX
           ENDIF
           
           IF (ISOLOP.EQ.1)THEN ! LINEAR STATIC CONDITION
           TIME  = 0.0D0
           NTIME = 61D0
           ELSEIF (ISOLOP.EQ.5) THEN ! LINEAR DYNAMIC
           CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') 
           ENDIF
           WNU  = 0.0D0
           NAME = TRIM(NAME_PARAMETER(I))
           NLENGTH = LEN_TRIM(NAME_PARAMETER(I))
           WRITE (67,100) NAME(1:NLENGTH)
           WRITE (67,102)
100        FORMAT (A)
           
           DO JJ = 1,NTIME ! TIME VARY
           IF (ISOLOP.EQ.1) THEN
             TIMEINVERVAL = PERIOD/20D0
             IF (JJ.EQ.1) TIME = 0.D0
             IF (JJ.NE.1) TIME = TIME + TIMEINVERVAL
           ELSEIF (ISOLOP.EQ.5) THEN
           CALL OFFSHSTEP(TIME,JJ,NTIME,'CALL')    
           ENDIF
           ! WATER SURFACE ELEVATION
           IF (IWAVE == 1  ) THEN 
           WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
           
           ELSEIF (IWAVE == 2) THEN
           WNU = (1.0D0/RK)*(F1*COS(1.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F2*COS(2.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F3*COS(3.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F4*COS(4.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F5*COS(5.0D0*(RK*H1M - OMEGA*TIME)))
           
           ELSEIF (IWAVE == 3 )THEN
           CALL STREAMWAVEFUNCTION (WVHIGHT,PERIOD,WDEPTH,TIME,WNU,GRAV,ORDER)
           
    	     ! FIND RAMDA AND WAVE NUMBER
    	     CALL VELOC (XC,YC,U,V,DUDT,DVDT,TT,ARAMDA)
    	     RAMDA = ARAMDA
    	     RK=2D0*3.141592654D0/RAMDA
           
    	     ELSEIF (IWAVE.EQ.4)THEN
    	     call Cnoidal_Wave_Crest  (WVHIGHT,WDEPTH,GRAV,PERIOD,TIME,Wave_Crest)
    	     call Cnoidal_wave_Length (WVHIGHT,WDEPTH,GRAV,PERIOD,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
    	     RK=wavenumber
    	     RAMDA=ALAMDA
    	     WNU = Wave_Crest
           
    	     ELSEIF (IWAVE.EQ.5)THEN
    	     Call Solitary_Wave (WVHIGHT,WDEPTH,GRAV,TIME,0d0,z,u,v,au,av,WNU,q)
           
           
           ELSEIF (IWAVE.EQ.6)THEN
           WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
           
           ELSEIF ( IWAVE == 7 .or. IWAVE == 8 ) THEN 
           WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
           
           ENDIF
           
           
           IF (WVHIGHT.EQ.0.OR.PERIOD.EQ.0.0D0) WNU = 0.0D0
           ELEVATION_DEPTH = WNU + WDEPTH
           WRITE (67,101) TIME,WNU,ELEVATION_DEPTH
101        FORMAT (E12.5,2X,E12.5,2X,E12.5)          
           
          ENDDO
          
           WRITE (67,104)
           WRITE (67,103)
           WRITE (67,110)
           
           IF (NGRAV.EQ.1)THEN
              ADATA_LENGTH = AMAX_MIN(5) - AMAX_MIN(6)
              H1M = AMAX_MIN(6)
              IF (ADATA_LENGTH.EQ.0) THEN
                  ADATA_LENGTH = 20D0
                  H1M = -10
              ENDIF
           ELSEIF (NGRAV.EQ.2)THEN
              ADATA_LENGTH = AMAX_MIN(1) - AMAX_MIN(2)
              H1M = AMAX_MIN(2)
              IF (ADATA_LENGTH.EQ.0) THEN
                  ADATA_LENGTH = 20D0
                  H1M = -10
              ENDIF
           ELSEIF (NGRAV.EQ.3)THEN
              ADATA_LENGTH = AMAX_MIN(1) - AMAX_MIN(2)
              H1M = AMAX_MIN(2)
              IF (ADATA_LENGTH.EQ.0) THEN
                  ADATA_LENGTH = 20D0
                  H1M = -10
              ENDIF
           ENDIF
           ADATA_INTERVAL = ADATA_LENGTH/60
           
          DO JJ = 1,61
           IF (JJ.EQ.1) H1M = H1M
           IF (JJ.GT.1) THEN
               H1M = H1M + ADATA_INTERVAL
           ENDIF
           ! WATER SURFACE ELEVATION
           IF (IWAVE == 1  ) THEN 
           WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
           
           ELSEIF (IWAVE == 2) THEN
           WNU = (1.0D0/RK)*(F1*COS(1.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F2*COS(2.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F3*COS(3.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F4*COS(4.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                       F5*COS(5.0D0*(RK*H1M - OMEGA*TIME)))
           
           ELSEIF (IWAVE == 3 )THEN
           CALL STREAMWAVEFUNCTION (WVHIGHT,PERIOD,WDEPTH,TIME,WNU,GRAV,ORDER)
           
    	     ! FIND RAMDA AND WAVE NUMBER
    	     CALL VELOC (XC,YC,U,V,DUDT,DVDT,TT,ARAMDA)
    	     RAMDA = ARAMDA
    	     RK=2D0*3.141592654D0/RAMDA
           
    	     ELSEIF (IWAVE.EQ.4)THEN
    	     call Cnoidal_Wave_Crest  (WVHIGHT,WDEPTH,GRAV,PERIOD,TIME,Wave_Crest)
    	     call Cnoidal_wave_Length (WVHIGHT,WDEPTH,GRAV,PERIOD,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
    	     RK=wavenumber
    	     RAMDA=ALAMDA
    	     WNU = Wave_Crest
           
    	     ELSEIF (IWAVE.EQ.5)THEN
    	     Call Solitary_Wave (WVHIGHT,WDEPTH,GRAV,TIME,0d0,z,u,v,au,av,WNU,q)
           
           
           ELSEIF (IWAVE.EQ.6)THEN
           WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
           
           ELSEIF ( IWAVE == 7 .or. IWAVE == 8 ) THEN 
           WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
           
           ENDIF
           
           
           IF (WVHIGHT.EQ.0.OR.PERIOD.EQ.0.0D0) WNU = 0.0D0
           ELEVATION_DEPTH = WNU + WDEPTH
           WRITE (67,101) H1M,WNU,ELEVATION_DEPTH 
           
          ENDDO
          
          WRITE (67,104)
          WRITE (67,103)
           
           ! FORMAT LOCATION
102        FORMAT ("Time    Water Surface Elevation(WSE)    WSE + Water Depth")
103        FORMAT ("")  
104        FORMAT ("END")
110        FORMAT ("COORDINATE(@W.L)    Water Surface Elevation(WSE)    WSE + Water Depth")
           
      ENDDO
      
      
      DO I =1,NSPECTRUM
      OFFSELECT = I    
      CALL SELECTSPECTRUM (OFFSELECT,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
      
          ! INITIAL CONDITION SET
               H1M  = 0.D0
           IF (NGRAV.EQ.1)THEN
               H1M = ZZ
           ELSEIF (NGRAV.EQ.2)THEN
               H1M = XX
           ELSEIF (NGRAV.EQ.3)THEN
               H1M = XX
           ENDIF
           IF (ISOLOP.EQ.1)THEN ! LINEAR STATIC CONDITION
           TIME  = 0.0D0
           NTIME = 61D0
           ELSEIF (ISOLOP.EQ.5) THEN ! LINEAR DYNAMIC
           CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') 
           ENDIF
           WNU  = 0.0D0
           NAME = TRIM(NAME_SPECTRUM(I))
           NLENGTH = LEN_TRIM(NAME_SPECTRUM(I))
           WRITE (67,100) NAME(1:NLENGTH)
           WRITE (67,102)
           
         DO JJ = 1,NTIME ! TIME VARY
           IF (ISOLOP.EQ.1) THEN
             TIMEINVERVAL = PERIOD/20D0
             IF (JJ.EQ.1) TIME = 0.D0
             IF (JJ.NE.1) TIME = TIME + TIMEINVERVAL
           ELSEIF (ISOLOP.EQ.5) THEN
           CALL OFFSHSTEP(TIME,JJ,NTIME,'CALL')    
           ENDIF
           IF (IRWAVE.EQ.1)THEN
           CALL Pierson_Moskowitz_RANDOM ("SUR",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY)  
           ELSEIF (IRWAVE.EQ.2) THEN
           CALL Jonswap_RANDOM ("SUR",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY) 
           ELSEIF (IRWAVE.EQ.3)THEN
           CALL Ochi_RANDOM ("SUR",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY) 
           ELSEIF (IRWAVE.EQ.4)THEN
           CALL Bretschneider_RANDOM ("SUR",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY) 
           ELSEIF (IRWAVE.EQ.5)THEN
           CALL TMA_SPectrum_RANDOM ("SUR",TIME,offselect,WNU,X,Y,UDX,UDY,AUX,AUY)     
           ENDIF
           
           IF (WVHIGHT.EQ.0.OR.PERIOD.EQ.0.0D0) WNU = 0.0D0
           ELEVATION_DEPTH = WNU + WDEPTH
           WRITE (67,101) TIME,WNU,ELEVATION_DEPTH
         ENDDO
         
           WRITE (67,104)
           WRITE (67,103)
           WRITE (67,110)
           
          IF (NGRAV.EQ.1)THEN
              ADATA_LENGTH = AMAX_MIN(5) - AMAX_MIN(6)
              H1M = AMAX_MIN(6)
           ELSEIF (NGRAV.EQ.2)THEN
              ADATA_LENGTH = AMAX_MIN(1) - AMAX_MIN(2)
              H1M = AMAX_MIN(2)
           ELSEIF (NGRAV.EQ.3)THEN
              ADATA_LENGTH = AMAX_MIN(1) - AMAX_MIN(2)
              H1M = AMAX_MIN(2)
           ENDIF
           ADATA_INTERVAL = ADATA_LENGTH/60
           
          
         DO JJ = 1,61
           IF (JJ.EQ.1) H1M = H1M
           IF (JJ.GT.1) THEN
               H1M = H1M + ADATA_INTERVAL
           ENDIF
           IF (IRWAVE.EQ.1)THEN
           CALL Pierson_Moskowitz_RANDOM ("SUR",TIME,offselect,WNU,H1M,Y,UDX,UDY,AUX,AUY)  
           ELSEIF (IRWAVE.EQ.2) THEN
           CALL Jonswap_RANDOM ("SUR",TIME,offselect,WNU,H1M,Y,UDX,UDY,AUX,AUY) 
           ELSEIF (IRWAVE.EQ.3)THEN
           CALL Ochi_RANDOM ("SUR",TIME,offselect,WNU,H1M,Y,UDX,UDY,AUX,AUY) 
           ELSEIF (IRWAVE.EQ.4)THEN
           CALL Bretschneider_RANDOM ("SUR",TIME,offselect,WNU,H1M,Y,UDX,UDY,AUX,AUY) 
           ELSEIF (IRWAVE.EQ.5)THEN
           CALL TMA_SPectrum_RANDOM ("SUR",TIME,offselect,WNU,H1M,Y,UDX,UDY,AUX,AUY)     
           ENDIF
           
           IF (WVHIGHT.EQ.0.OR.PERIOD.EQ.0.0D0) WNU = 0.0D0
           ELEVATION_DEPTH = WNU + WDEPTH
           WRITE (67,101) H1M,WNU,ELEVATION_DEPTH
         ENDDO
         
           WRITE (67,104)
           WRITE (67,103)
      ENDDO
          
      ENDIF
          END
C	=====================================================================
          
      SUBROUTINE  MID_STRUCTURE (XYZ,XX,YY,ZZ,AMAX_MIN) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION XYZ(1)
      DIMENSION AMAX_MIN(6)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      
      AXX = 0.0D0
      AYY = 0.0D0
      AZZ = 0.0D0
      
      INDEX_START = 1D0
      INDEX_END = NSN
      AXX = SUM(XYZ(INDEX_START:INDEX_END))
      AMAX_MIN(1) = MAXVAL(XYZ(INDEX_START:INDEX_END))*5.0D0
      AMAX_MIN(2) = MINVAL(XYZ(INDEX_START:INDEX_END))*5.0D0 
      
      INDEX_START = NSN+1
      INDEX_END = NSN*2D0
      AYY = SUM(XYZ(INDEX_START:INDEX_END))
      AMAX_MIN(3) = MAXVAL(XYZ(INDEX_START:INDEX_END))*5.0D0
      AMAX_MIN(4) = MINVAL(XYZ(INDEX_START:INDEX_END))*5.0D0
      
      INDEX_START = NSN*2D0 + 1D0
      INDEX_END = NSN*3D0
      AZZ = SUM(XYZ(INDEX_START:INDEX_END))
      AMAX_MIN(5) = MAXVAL(XYZ(INDEX_START:INDEX_END))*5.0D0
      AMAX_MIN(6) = MINVAL(XYZ(INDEX_START:INDEX_END))*5.0D0
      
      XX = AXX/NSN
      YY = AYY/NSN
      ZZ = AZZ/NSN
      IF (XX.LT.0.00001) XX = 0.D0
      IF (XX.LT.-0.00001) XX = 0.D0
      IF (YY.LT.0.00001) YY = 0.D0
      IF (YY.LT.-0.00001) YY = 0.D0
      IF (ZZ.LT.0.00001) ZZ = 0.D0
      IF (ZZ.LT.-0.00001) ZZ = 0.D0
      
      
      
      END