      SUBROUTINE READSPECTRUM (NSPECTRUM)      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      
      CHARACTER*20 NAME_READ
      CHARACTER*200 NAME_PARAMETER,NAME_SPECTRUM
      CHARACTER*250 NAME
      
      
      DIMENSION WAVE(3),WIND(3)
      DIMENSION VH1(3),VH2(3),VGV(3)
      COMMON /INOU/ ITI 
      COMMON /MGRAV/ NGRAV
      COMMON/ RES_OFFSHOREX / SEABEDRX(1000),WVHIGHTRX(1000),WDEPTHRX(1000),H1POSRX(1000),H2POSRX(1000),
     1                        NRESWAVERX(1000),PERIODRX(1000),GRAVRX(1000),RHOWRX(1000),RHOARX(1000),WKFRX(1000),WFCRX(1000),
     1                        FREQENCYX(1000),DAMPINGX(1000),
     1                        WVHIGHTRX1(1000),WVHIGHTRX2(1000),PERIODRX1(1000),PERIODRX2(1000),SHAPEX1(1000),SHAPEX2(1000),
     1                        TAPX(1000),UCURRENTX(1000),
     1                        NSWINDX(1000),UMEANWINDSPECX(1000),POWERLAWSPECX(1000),ROUGHSPECX(1000),REFERSPECX(1000),
     1                        TURFACTORX(1000),SCALEFACTORX(1000),STDFACTORX(1000),
     1                        POINT(1000,1000),SFREQ(1000,1000),AMPLI(1000,1000),
     1                        VGX(1000,3),V1X(1000,3),V2X(1000,3),VWX(1000,3),
     1                        FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000)
C    1                        NOPSPECX(1000),TIMESTARTX(1000),TIMEINX(1000),TIMEENDX(1000),
      COMMON /NAME_OFFSHORE_PARAMETER/ NAME_PARAMETER(1000),NAME_SPECTRUM(1000)
     
      DO I = 1,NSPECTRUM
      READ(ITI,*)! --------------------------------------------------------------------------
      READ(ITI,'(A)') NAME(1:250)!OFFSHORE LOAD PARAMETERS CASE ---
      NAME_SPECTRUM(I) = NAME(38:250)
      READ(ITI,*)! --------------------------------------------------------------------------
      READ(ITI,*)!---------------------------- WAVE PARAMETERS------------------------------
      
      READ(ITI,*) SEABEDRX(I) ,NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) WVHIGHTRX(I) ,NAME_READ,NAME_READ,NAME_READ   
      READ(ITI,*) WDEPTHRX(I) ,NAME_READ,NAME_READ,NAME_READ                    
      READ(ITI,*) H1POSRX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) H2POSRX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) NRESWAVERX(I) ,NAME_READ,NAME_READ,NAME_READ      
      READ(ITI,*) PERIODRX(I) ,NAME_READ,NAME_READ,NAME_READ        
      READ(ITI,*) GRAVRX(I) ,NAME_READ,NAME_READ,NAME_READ          
      READ(ITI,*) RHOWRX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) RHOARX(I) ,NAME_READ,NAME_READ,NAME_READ         
      READ(ITI,*) WKFRX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) WFCRX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      
      READ(ITI,*) WVHIGHTRX1(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) WVHIGHTRX2(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) PERIODRX1(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) PERIODRX2(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) SHAPEX1(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) SHAPEX2(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ 
      
      !READ(ITI,*) FREQENCYX(I) , NAME_READ,NAME_READ,NAME_READ
      !READ(ITI,*) DAMPINGX(I)  , NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) WAVEANGLE    ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      
   
      READ(ITI,*)
      READ(ITI,*) FSpectrumStartX(I),NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) FSpectrumEndX(I),NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) NUMBEROFRANDOMWANVEXInterationX(I),NAME_READ,NAME_READ,NAME_READ,NAME_READ

      
      READ(ITI,*)!"---------------------------- CURRENT PARAMETERS-----------------------------"
      READ(ITI,*) UCURRENTX(I)  , NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) TAPX(I)       , NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*)!"----------------------------- WIND PARAMETERS------------------------------"
      
      
      READ(ITI,*) NSWINDX(I)        ,NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) UMEANWINDSPECX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) POWERLAWSPECX(I)  ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) ROUGHSPECX(I)     ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) REFERSPECX(I)     ,NAME_READ,NAME_READ,NAME_READ
      
      READ(ITI,*) TURFACTORX(I)   ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) SCALEFACTORX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) STDFACTORX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ   
      READ(ITI,*) WINDANGLE   ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      
      READ(ITI,*)
      

      
!      !READ(ITI,*) SEABEDRX(I),WVHIGHTRX(I),WDEPTHRX(I),H1POSRX(I),H2POSRX(I)
!      !READ(ITI,*) NRESWAVERX(I),PERIODRX(I),GRAVRX(I),RHOWRX(I),RHOARX(I),WKFRX(I),WFCRX(I)
!      READ(ITI,*) WVHIGHTRX1(I),WVHIGHTRX2(I),PERIODRX1(I),PERIODRX2(I),SHAPEX1(I),SHAPEX2(I) 
!      READ(ITI,*) FREQENCYX(I),DAMPINGX(I)  
!      
!      READ(ITI,*) UCURRENTX(I),TAPX(I)
!      READ(ITI,*) NSWINDX(I),UMEANWINDSPECX(I),POWERLAWSPECX(I),ROUGHSPECX(I),REFERSPECX(I)
!      READ(ITI,*) TURFACTORX(I),SCALEFACTORX(I),STDFACTORX(I)
!C      READ(ITI,*) NOPSPECX(I),TIMESTARTX(I),TIMEINX(I),TIMEENDX(I)
!      READ(ITI,*) WAVEANGLE
!      READ(ITI,*) WINDANGLE
  
      

C    --------------------------------      
C      IF (NOPSPECX(I).EQ.2)THEN
C      READ(ITI,*)
C      READ(ITI,*)
C      READ(ITI,*) NSTEP
C      READ(ITI,*)
      
C       DO J = 1,NSTEP
C       READ (ITI,*) POINT(I,J),SFREQ(I,J),AMPLI(I,J)
C       ENDDO
C      ENDIF
C    --------------------------------
      
      
      GRAVRX(I) = ABS(GRAVRX(I))
C     ================================================================ 
C        MODIFY BY TOEY 30-08-2012 FOR WAVE ANGLE AND WIND ANGLE
C     ================================================================
      CALL WAVE_VECTOR (WAVEANGLE,WAVE)
      CALL WIND_VECTOR (WINDANGLE,WIND)
      
      WVZETA = 0.0D0
      WVTIME = 0.0D0
      GRAV = ABS(GRAV)
      
	WAVE(NGRAV) = 0.0D0 !SET WAVE VECTOR IN GRAVITY DIRECTION TO BE ZERO
	WIND(NGRAV) = 0.0D0 !SET WIND VECTOR IN GRAVITY DIRECTION TO BE ZERO

	CALL SCALEN(WAVE,WAVE,DUM,3)
	CALL SCALEN(WIND,WIND,DUM,3)
	
	VGV(1:3)   = 0.0D0
	VGV(NGRAV) = 1.0D0
	VH1(1:3)   = WAVE(1:3)
	
      CALL VECPRD (VGV,VH1,VH2)
	CALL SCALEN(VH2,VH2,DUM,3)

      CALL INTFILL ('UNIT',IUL,1,2,0)
	FAL  = UNITFAC('LENG',IUL,0)
	Z0   = Z0*FAL  !Z0 in (m) unit
C     -------------------------

      VGX(I,1:3)    =   VGV(1:3) 
      V1X(I,1:3)    =   VH1(1:3)
      V2X(I,1:3)    =   VH2(1:3)
      VWX(I,1:3)    =   WIND(1:3)
     
      ENDDO
      
      
               
      END
      
C   ==========================================================================================================           
      SUBROUTINE READGRAPH
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /INOU/ ITI 
      COMMON /GRAPHRES/ POINTX(500),FREQX(500),AVALUEX(500)
      COMMON /NGRES/ DETA_FREQ,NBERX
      ! NONLINER GRAPH
      READ (ITI,*)
      READ (ITI,*)
      READ (ITI,*) NLOOP
      READ (ITI,*)

      
      DO I = 1,NLOOP
      READ (ITI,*) POINTX(I),FREQX(I),AVALUEX(I)
      ENDDO
      
      READ (ITI,*)
      READ (ITI,*) NBERX,DETA_FREQ
      
      END
C   ==========================================================================================================           
      SUBROUTINE SELECTGRAPH (NFUNCTION,NUMBER,NPOI,FREQ)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      COMMON /INOU/ ITI 
      COMMON /GRAPHRES/ POINTX(500),FREQX(500),AVALUEX(500)
      COMMON /NGRES/ DETA_FREQ,NBERX
      
      IF (NFUNCTION.EQ.1) NUMBER = NBERX
      
      IF (NFUNCTION.EQ.2) THEN
         IF (NPOI.EQ.1) THEN
           IF (FREQX(1).EQ.0)THEN
           FREQ = 0.0D0
           ELSEIF (FREQX(1).NE.0)THEN
           FREQ = FREQX(1)    
           ENDIF
         ENDIF
         IF (NPOI.GT.1) THEN
         FREQ = FREQ + DETA_FREQ
         ENDIF
      ENDIF
      
      IF (NFUNCTION.EQ.3) THEN
         
      ENDIF
      
      END
C   ==========================================================================================================      
      SUBROUTINE SELECTSPECTRUM (OFFSELECT,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,NRESWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                           FREQENCY,VG,V1,V2,VWIND,DAMPING,NSWIND,UMEANWINDSPEC,POWERLAWSPEC,ROUGHSPEC,REFERSPEC,FWIND,
     1                           TAP,UCURRENT,WVHIGHTR1,WVHIGHTR2,PERIODR1,PERIODR2,SHAPE1,SHAPE2)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION VWAVE(3)
      DIMENSION FWIND(3)
      DIMENSION VG(3),V1(3),V2(3),VWIND(3)
      COMMON/ RES_OFFSHOREX / SEABEDRX(1000),WVHIGHTRX(1000),WDEPTHRX(1000),H1POSRX(1000),H2POSRX(1000),
     1                        NRESWAVERX(1000),PERIODRX(1000),GRAVRX(1000),RHOWRX(1000),RHOARX(1000),WKFRX(1000),WFCRX(1000),
     1                        FREQENCYX(1000),DAMPINGX(1000),
     1                        WVHIGHTRX1(1000),WVHIGHTRX2(1000),PERIODRX1(1000),PERIODRX2(1000),SHAPEX1(1000),SHAPEX2(1000),
     1                        TAPX(1000),UCURRENTX(1000),
     1                        NSWINDX(1000),UMEANWINDSPECX(1000),POWERLAWSPECX(1000),ROUGHSPECX(1000),REFERSPECX(1000),
     1                        TURFACTORX(1000),SCALEFACTORX(1000),STDFACTORX(1000),
     1                        POINT(1000,1000),SFREQ(1000,1000),AMPLI(1000,1000),
     1                        VGX(1000,3),V1X(1000,3),V2X(1000,3),VWX(1000,3),
     1                        FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000)
C    1                        NOPSPECX(1000),TIMESTARTX(1000),TIMEINX(1000),TIMEENDX(1000),
     
      SEABED    = SEABEDRX(OFFSELECT)
      WVHIGHT   = WVHIGHTRX(OFFSELECT)
      WDEPTH    = WDEPTHRX(OFFSELECT)
      H1POS     = H1POSRX(OFFSELECT)
      H2POS     = H2POSRX(OFFSELECT)
      NRESWAVE  = NRESWAVERX(OFFSELECT)
      PERIOD    = PERIODRX(OFFSELECT)
      GRAV      = GRAVRX(OFFSELECT)
      RHOW      = RHOWRX(OFFSELECT)
      RHOA      = RHOARX(OFFSELECT)
      WKF       = WKFRX(OFFSELECT)
      WFC       = WFCRX(OFFSELECT)
      FREQENCY  = FREQENCYX(OFFSELECT)
      DAMPING   = DAMPINGX(OFFSELECT)
      VG(1:3)   = VGX(OFFSELECT,1:3)
      V1(1:3)   = V1X(OFFSELECT,1:3)
      V2(1:3)   = V2X(OFFSELECT,1:3)
      VWIND(1:3)= VWX(OFFSELECT,1:3)
      
      WVHIGHTR1 = WVHIGHTRX1(OFFSELECT)
      WVHIGHTR2 = WVHIGHTRX2(OFFSELECT)
      PERIODR1  = PERIODRX1(OFFSELECT)
      PERIODR2  = PERIODRX2(OFFSELECT)
      SHAPE1    = SHAPEX1(OFFSELECT)
      SHAPE2    = SHAPEX2(OFFSELECT)
      
      TAP       = TAPX(OFFSELECT)
      UCURRENT  = UCURRENTX(OFFSELECT)
      
      NSWIND         = NSWINDX(OFFSELECT)
      UMEANWINDSPEC  = UMEANWINDSPECX(OFFSELECT)
      POWERLAWSPEC  = POWERLAWSPECX(OFFSELECT)
      ROUGHSPEC     = ROUGHSPECX(OFFSELECT)
      REFERSPEC     = REFERSPECX(OFFSELECT)
      
      FWIND(1)   = TURFACTORX(OFFSELECT)
      FWIND(2)   = SCALEFACTORX(OFFSELECT)
      FWIND(3)   = STDFACTORX(OFFSELECT)
      
!      NOPSPEC    = FSpectrumStartX(OFFSELECT)
!      TIMESTART  = FSpectrumEndX(OFFSELECT)
!      TIMEIN     = NUMBEROFRANDOMWANVEXInterationX(OFFSELECT)

C      NOPSPEC    = NOPSPECX(OFFSELECT)
C      TIMESTART  = TIMESTARTX(OFFSELECT)
C      TIMEIN     = TIMEINX(OFFSELECT)
C      TIMEEND    = TIMEENDX(OFFSELECT)
      
      END

C   ==========================================================================================================        
      SUBROUTINE NODALSPEC (NTLAP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)     
      COMMON /INOU/ ITI 
      COMMON /SPECNODE/ STEPFREQ(1000,1000),FORCE(1000,1000),NODE(1000),NLAP(1000)
      
      
      DO I = 1,NTLAP
      READ (ITI,*)
      READ (ITI,*) NLAP(I),NODE(I)
      READ (ITI,*) ! POINT NO.     FREQUENCY     FORCE
      
         DO J = 1,NLAP(I)
         READ (ITI,*) NUM,STEPFREQ(I,J),FORCE(I,J)
         ENDDO
      
      ENDDO
 
      CALL NODALDFORCESPEC
      END     
 
C   ==========================================================================================================               
      SUBROUTINE NODALDFORCESPEC
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)  
      COMMON /SPECNODE/ STEPFREQ(1000,1000),FORCE(1000,1000),NODE(1000),NLAP(1000)
      
      DO I = 1,ITIME
      CALL SELECTGRAPH (2,0,ITIME,TIME)
      ! INTRERPORATION FUNCTION
          DO J = 1,NLAP(I)
          STEP = STEPFREQ(I,J)-STEPFREQ(I,J+1)
          IF (TIME.LE.STEPFREQ(I,J+1).AND.TIME.GE.STEPFREQ(I,J))THEN
             
          ENDIF
          ENDDO
      ENDDO
      
      END   
C   ==========================================================================================================          

      SUBROUTINE SPECTURMPARAMETER (MLE,NFUNCTION,ROUGH,CD,CM,CSA,NDI,DIAM1,NGROWTH,GROWTH,OFFSELECT_X,LWCASE,ILCN,ILCC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)      
      DIMENSION LWCASE(8)
      COMMON /INOU/ ITI  
      READ(ITI,*) MLE,NFUNCTION,ROUGH,CD,CM,CSA,NDI,DIAM1,NGROWTH,GROWTH,OFFSELECT_X,LWCASE(8),
     1              LWCASE(5),LWCASE(6),LWCASE(7),ILCN,ILCC !READ ONLY FIRST TIME STEP
      
      END
C     ======================================================================================================         
      SUBROUTINE STOREMODESHAPE (RV,NR)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)  
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON / STOREMODE / RVEC(100000,100),NMOD
      DIMENSION RV(NEQ,NR)
      RVEC(1:NEQ,1:NR) = RV(1:NEQ,1:NR)
      NMOD = NR
      END
      
C     ======================================================================================================            
      SUBROUTINE MUTIFORCE (ITIME,ICASE,ISPEC,NOFFL)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      CHARACTER*4 OPER
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
     	COMMON /NUMB/ HED(20),MODEX,NRE,NSN1,NEG,NBS,NLS,NLA,
     1              NSC1,NSF1,IDOF1(9),LCS,ISOLOP,LSYMM
     	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 
      COMMON / STOREMODE / RVECT(100000,100),NR
     	
     	DIMENSION R(NEQ),SHAPET(NEQ),C(NEQ),CV(NEQ),RVEC(NEQ),RV( NEQ,NR),RU(NEQ)
     	DIMENSION CVEC(NEQ),AVEC(NEQ)
     	
      RV(1:NEQ,1:NR) = RVECT(1:NEQ,1:NR) 

C     ==========================
C             INITIAL SET
C     ========================== 
      ICOM = 1 ! SRSS METHOD
      
C     ==========================      
C            CLEAR MARTIX
C     ==========================
      DO K = 1,NEQ
      AVEC(K) = 0.0D0 
      ENDDO
      
      SELECTCASE (ICOM)
C     ===========================
C       COMBINATION SRSS METHOD
C     ===========================
      CASE (1)
      DO I = 1,NEQ
         DO J = 1,NR
         AVEC(I) = (ABS(AVEC(I))) + (ABS(RV(I,J)))
         CVEC(I) = SQRT(AVEC(I)**2D0)
         ENDDO
         RVEC(I) = CVEC(I)
      ENDDO 
      
C     ===========================
C       COMBINATION CQC METHOD
C     ===========================      
      CASE (2)
      
      ENDSELECT
      
C     ===========================
C       READING OFFSHORE DATA
C     ===========================      
      OPER  = 'REDT'
      CALL OFFSHFORC(R,ICASE,ITIME,NEQ,'SPEC',OPER)
      
C      WRITE (300,1) R(1:NEQ)

C     --------------------------
C     MODAL MATRIX *FORCE MATRIX 
C     --------------------------  
      DO I=1,NEQ
      SHAPET(I)   = RVEC(I)*R(I) 
      CV(I)       = RVEC(I)*R(I)     
C      WRITE (299,1) SHAPET(I)  
      ENDDO
      
C     ==========================
C        UPDATE OFFSHORE DATA
C     ==========================  
      IF(ISPEC.GT.1)THEN
   	OPER  = 'REDT'
	CALL OFFSHFORC(RU,ICASE,ITIME,NEQ,'VARY',OPER)   
	DO I=1,NEQ
	C(I) = RU(I)+ CV(I)
	ENDDO  
	ENDIF

C     ==========================
C        WRITE OFFSHORE DATA
C     ==========================  
      IF (ISPEC.EQ.1)THEN    
      OPER  = 'WRIT'
      CALL OFFSHFORC(CV,ICASE,ITIME,NEQ,'VARY',OPER) 
      ELSEIF (ISPEC.GT.1)THEN
      OPER  = 'WRIT'
      CALL OFFSHFORC(C,ICASE,ITIME,NEQ,'VARY',OPER) 
      ENDIF
C      WRITE (299,1) C(1:NEQ)
      
1     FORMAT (E12.5)
        
      END

C     ======================================================================================================            
      SUBROUTINE MUTIFORCETOEY (ITIME,NTOEY,ISPEC,NOFFL)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      CHARACTER*4 OPER
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
     	COMMON /NUMB/ HED(20),MODEX,NRE,NSN1,NEG,NBS,NLS,NLA,
     1              NSC1,NSF1,IDOF1(9),LCS,ISOLOP,LSYMM
     	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 
      COMMON / STOREMODE / RVECT(100000,100),NR
     	
     	DIMENSION R(NEQ),SHAPET(NEQ),C(NEQ),CV(NEQ),RVEC(NEQ),RV(NEQ,NR),RU(NEQ)
     	DIMENSION CVEC(NEQ),AVEC(NEQ)
     	
      RV(1:NEQ,1:NR) = RVECT(1:NEQ,1:NR) 

C     ==========================
C             INITIAL SET
C     ========================== 
      ICOM = 1 ! SRSS METHOD
      
C     ==========================      
C            CLEAR MARTIX
C     ==========================
      DO K = 1,NEQ
      AVEC(K) = 0.0D0 
      ENDDO
      
      SELECTCASE (ICOM)
C     ===========================
C       COMBINATION SRSS METHOD
C     ===========================
      CASE (1)
      DO I = 1,NEQ
         DO J = 1,NR
         AVEC(I) = (ABS(AVEC(I))) + (ABS(RV(I,J)))
         CVEC(I) = SQRT(AVEC(I)**2D0)
         ENDDO
         RVEC(I) = CVEC(I)
      ENDDO 
      
C     ===========================
C       COMBINATION CQC METHOD
C     ===========================      
      CASE (2)
      
      ENDSELECT
      
C     ===========================
C       READING OFFSHORE DATA
C     ===========================      
      OPER  = 'REDT'
      CALL OFFSHFORC(R,NTOEY,ITIME,NEQ,'SPEC',OPER)
      

C     --------------------------
C     MODAL MATRIX *FORCE MATRIX 
C     --------------------------  
      DO I=1,NEQ
      CV(I)       = RVEC(I)*R(I)      
      ENDDO
      
C     ============================
C        UPDATE OFFSHORE DATA
C        MODE SHAPE SUPERPOSITION
C     ============================ 
      IF(ISPEC.GT.1)THEN
   	OPER  = 'REDT'
	CALL OFFSHFORC(RU,NTOEY,ITIME,NEQ,'VARY',OPER)   
	DO I=1,NEQ
	! SUPERPOSITION 
	C(I) = RU(I)+ CV(I)
	ENDDO  
	ENDIF

C     ==========================
C        WRITE OFFSHORE DATA
C     ==========================  
      IF (ISPEC.EQ.1)THEN    
      OPER  = 'WRIT'
      CALL OFFSHFORC(CV,NTOEY,ITIME,NEQ,'VARY',OPER) 
      ELSEIF (ISPEC.GT.1)THEN
      OPER  = 'WRIT'
      CALL OFFSHFORC(C,NTOEY,ITIME,NEQ,'VARY',OPER) 
      ENDIF
      
      END
C     ======================================================================================================   
      Subroutine  RANDOM_NUMBER_ROUTINE(FREQ_DATA,ANUMBR_RANDOM)
      USE IFPORT
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION ANUMBR_RANDOM(FREQ_DATA)
      CALL RANDOM_NUMBER(ANUMBR_RANDOM)
      END
        
C     ======================================================================================================      
      Subroutine  Pierson_Moskowitz_SPectrum(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SFUU,SFVV,SNA,SDG,Tapp,Ucurent)
        USE IFPORT
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        COMMON / WAVESPECTRUMPLOT / NWSPECTRUMPLOT

        DIMENSION SNA(10),SDG(10)
        
        DIMENSION SPM_TEST(400)
        
        DIMENSION ANUMBR_RANDOM(50)
        DIMENSION W(50),WT(50),PHI(50)
        DIMENSION SURFACE_ELEVATION(1600)
        
        !INTEGER :: ISEED = 43566
        !REAL*8  :: RI(5)!,IIIN(5)
        !INTEGER :: IIIN(5)
        
        !DIMENSION NRR(10)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        Uapp = Ucurent
        Fp=1/Tp
        !f=0.15
        ! F  = 0.
        !Hs = 5
        !SPM_TEST = 0.
        !N_FREQ = 400
        !DO I = 1,N_FREQ
        !f = 0.1d0    
        

        SPM=0.3125d0*(Hs**2)*(Fp**4D0)*(f**-5d0)*exp(-1.25D0*(Fp/f)**4D0)
        
       ! FREQ_START    = 0.2D0
       ! FREQ_END      = 2.5D0
       ! FREQ_DATA     = 50D0
       ! 
       ! FREQ_INTEVAL  = (FREQ_END-FREQ_START)/(FREQ_DATA-1)
       ! 
       ! DELTA_W       = FREQ_INTEVAL
       ! ! ------ FROM THE STACK OVERFLOW ------ (RANDOM NUMBER)
       ! CALL RANDOM_NUMBER_ROUTINE (FREQ_DATA,ANUMBR_RANDOM)
       ! 
       ! DO JJ = 1,FREQ_DATA
       ! IF (JJ.EQ.1) W(JJ)  = FREQ_START
       ! IF (JJ.NE.1) W(JJ)  = W(JJ-1) + DELTA_W
       ! ENDDO
       ! 
       ! DO JJ =  1,FREQ_DATA
       ! WT(JJ)  = W(JJ) + DELTA_W*ANUMBR_RANDOM(JJ)
       ! PHI(JJ) = 2D0*3.141592654D0*(ANUMBR_RANDOM(JJ)-0.5D0)
       ! ENDDO
       ! 
       ! ! -------------------------------------
        
        !SURFACE_ELEVATION = 0.
        !DO JJ = 1,FREQ_DATA
        !DO KK = 1,1600
        !IF (KK.EQ.1) TIME_VECTOR = 0.0D0
        !IF (KK.NE.1) TIME_VECTOR = TIME_VECTOR + 0.25D0
        !F           = W(JJ)
        !SPM         = 0.3125d0*(Hs**2)*(Fp**4)*(f**-5d0)*exp(-1.25*(Fp/f)**4)
        !AMPLITUDE   = SQRT(2D0*SPM*DELTA_W)
        !SURFACE_ELEVATION(KK) = SURFACE_ELEVATION(KK) + AMPLITUDE*(COS(F*TIME_VECTOR+PHI(JJ)))
        !ENDDO
        !ENDDO
        
        !WRITE (300,'(E12.5)') SURFACE_ELEVATION
        
        !SPM_A=0.0081D0*(G**2)*((3.141592654D0*2)**-4)*(f**-5)*EXP(-1.25*((Fp/f)**4))
        IF (F.EQ.0.0) SPM = 0.0D0
        
        !SPM_TEST(I)= SPM
        !F = F + 0.0025D0
        !ENDDO
        
       ! CALL FAST_FOURIER_FORWARD_TEST (N_FREQ,SPM_TEST)
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
C       Wave Current Inter Action , Hedges et al.(1985)
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SPMU= BB/AA*SPM
        else
        SPMU= SPM
        endif
        
        Pi=3.141592654
        
        !SPMU       = 0.93842310D0
        !TP         = 20D0
        !RAMDA      = 1588.46D0D
        !WAVENUMBER = 2*PI/RAMDA
        
        AL = 2D0* Pi / WAVENUMBER
        
        CONSTANT = ((2D0*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)
        
        SUUX=((2D0*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SVVX=((2D0*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SAUUX=((2D0*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SPMU
        SAVVX=((2D0*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SPMU
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SPMU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SPMU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4D0)**2d0)*(8D0/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4D0)**2d0)*(8D0/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
        
        IF ( NWSPECTRUMPLOT .EQ. 1 ) THEN
        WRITE(*,10) f , SPM , AHH1
!        WRITE(5041,10) f , SPM , AHH1  ! TOEY 10/2021
10      FORMAT(F10.3,2X,F12.5,2X,F12.5,2X,'Pierson Moskowitz')
        NWSPECTRUMPLOT = 2
        ENDIF
            
 
        SFUU = sqrt(SFUU*AHH1)
        SFVV = sqrt(SFVV*AHH1)
    
        
        RETURN
        END SUBROUTINE
        
C     ======================================================================================================           
        !**********************************************************
        !For Solid Element
        !**********************************************************
        
        Subroutine  Pierson_Moskowitz_SPectrum_Solid(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SNA,SDG,VNOL,VR,DIAM,
     1                                               PTX,PTY,PTZ,Tapp,Ucurent)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
         Uapp = Ucurent
        Fp=1/Tp
        !f=0.15
        
        SPM=0.3125d0*(Hs**2)*(Fp**4)*(f**-5d0)*exp(-1.25*(Fp/f)**4)
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
C       Wave Current Inter Action , Hedges et al.(1985)
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SPMU= BB/AA*SPM
        else
        SPMU= SPM
        endif
        
        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SPMU
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SPMU
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SPMU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SPMU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SPMU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        
        Solid_SFUU  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SUU)
        Solid_SFUAU = ((1d0)*SAUU)
        
        Solid_SFVV  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SVV)
        Solid_SFVAV = ((1d0)*SAVV)
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)

C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )

        RETURN
        END SUBROUTINE
        
C     ======================================================================================================            
      Subroutine  Jonswap_SPectrum(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SFUU,SFVV,SNA,SDG,Tapp,Ucurent)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION SNA(10),SDG(10)  
      COMMON / WAVESPECTRUMPLOT / NWSPECTRUMPLOT
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        Uapp = Ucurent
        Fp=1/Tp
        !f=0.15
        
        if (f.LE.Fp) then
        St=0.07d0
        else
        St=0.09d0
        endif
        
        Alpha=exp(-((f-Fp)**2)/(2d0*St**2d0*Fp**2))
        
        if(Tp/sqrt(Hs).LE.3.6d0)then
        Gamma=5
        elseif(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) then
        Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
        elseif(Tp/sqrt(Hs).GT.5.0d0)then
        Gamma=1d0
        endif
        
        !Gamma=3.3d0 
        
        SJS=0.3125d0*(Hs**2)*(1/Fp)*((f/Fp)**-5d0)*exp(-1.25*(Fp/f)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((f/Fp-1)/st)**2)
        
      IF ( NWSPECTRUMPLOT .EQ. 1 ) THEN
        WRITE(*,10) f , SJS
!        WRITE(5041,10) f , SJS   ! TOEY 10/2021
10      FORMAT(F10.3,2X,F12.5,2X,'Jonswap_SPectrum')
        NWSPECTRUMPLOT = 2
      ENDIF
                
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SJSU= BB/AA*SJS
        else
        SJSU= SJS
        endif
        
        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SJSU
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SJSU
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SJSU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SJSU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        SFUU=sqrt(SFUU*AHH1)
        SFVV=sqrt(SFVV*AHH1)
        
        RETURN
        END SUBROUTINE
       
C     ======================================================================================================           
        Subroutine  Jonswap_SPectrum_Solid(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SNA,SDG,VNOL,VR,DIAM,
     1                                               PTX,PTY,PTZ,Tapp,Ucurent)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------

        Uapp = Ucurent
        Fp=1/Tp
        !f=0.15
        
        if (f.LE.Fp) then
        St=0.07d0
        else
        St=0.09d0
        endif
        
        Alpha=exp(-((f-Fp)**2)/(2d0*St**2d0*Fp**2))
        
        if(Tp/sqrt(Hs).LE.3.6d0)then
        Gamma=5
        elseif(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) then
        Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
        elseif(Tp/sqrt(Hs).GT.5.0d0)then
        Gamma=1d0
        endif
        
        
        SJS=0.3125d0*(Hs**2)*(1/Fp)*((f/Fp)**-5d0)*exp(-1.25*(Fp/f)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((f/Fp-1)/st)**2)
        
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SJSU= BB/AA*SJS
        else
        SJSU= SJS
        endif

        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SJSU
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SJSU
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SJSU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SJSU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SJSU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        
        Solid_SFUU  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SUU)
        Solid_SFUAU = ((1d0)*SAUU)
        
        Solid_SFVV  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SVV)
        Solid_SFVAV = ((1d0)*SAVV)
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        END SUBROUTINE
        
C     ======================================================================================================            
        Subroutine  Ochi_SPectrum(Hs1,Tp1,AShape1,Hs2,Tp2,AShape2,G,d,ROH,Cd,CI,Dia,Z,f,SFUU,SFVV,SNA,SDG,Tapp,Ucurent)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(10),SDG(10)
        
        COMMON / WAVESPECTRUMPLOT / NWSPECTRUMPLOT
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        Uapp = Ucurent
        Fp=1/Tp1
 		
        Pi=3.141592653589793d0
        Call GAMMA_FUNCTION(AShape1,GAM1)
        Call GAMMA_FUNCTION(AShape2,GAM2)
    
        OMEGAM1=2d0*Pi/Tp1
        OMEGAM2=2d0*Pi/Tp2
        OMEGA=2d0*Pi*f
        SOH1=(Hs1**2)/(4d0*GAM1)*(((4d0*AShape1+1d0)/4d0*OMEGAM1**4d0)**AShape1)*(1d0/(OMEGA**(4d0*AShape1+1d0)))
     1       *EXP(((4d0*AShape1+1d0)/-4d0)*(OMEGAM1/OMEGA)**4d0)
        SOH2=(Hs2**2)/(4d0*GAM2)*(((4d0*AShape2+1d0)/4d0*OMEGAM2**4d0)**AShape2)*(1d0/(OMEGA**(4d0*AShape2+1d0)))
     1       *EXP(((4d0*AShape2+1d0)/-4d0)*(OMEGAM2/OMEGA)**4d0)
        
        SOH=SOH1+SOH2
        
        IF ( NWSPECTRUMPLOT .EQ. 1 ) THEN
        WRITE(*,10) f , SOH
!       WRITE(5041,10) f , SOH   ! TOEY 10/2021
10      FORMAT(F10.3,2X,F12.5,2X,'Ochi Spectrum')
        NWSPECTRUMPLOT = 2
        ENDIF
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP1
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
C       Wave Current Inter Action , Hedges et al.(1985)
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SOHU= BB/AA*SOH
        else
        SOHU= SOH
        endif

        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SOHU
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SOHU
        
        
        SUU=((G*Tp1/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SVV=((G*Tp1/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SAUU=((G*Tp1/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SOHU
        SAVV=((G*Tp1/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SOHU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs1/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs1/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        SFUU = sqrt(SFUU*AHH1)
        SFVV = sqrt(SFVV*AHH1)
        
        
        RETURN
        END SUBROUTINE
        
C     ======================================================================================================    		
        Subroutine  Ochi_SPectrum_Solid(Hs1,Tp1,AShape1,Hs2,Tp2,AShape2,G,d,ROH,Cd,CI,Dia,Z,f,SNA,SDG,VNOL,VR,DIAM,PTX,PTY,PTZ
     1                                  ,Tapp,Ucurent)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        Uapp = Ucurent
        Fp=1/Tp1
 		
		Pi=3.141592653589793d0
        Call GAMMA_FUNCTION(AShape1,GAM1)
        Call GAMMA_FUNCTION(AShape2,GAM2)
    
        OMEGAM1=2d0*Pi/Tp1
        OMEGAM2=2d0*Pi/Tp2
        OMEGA=2d0*Pi*f
        SOH1=(Hs1**2)/(4d0*GAM1)*(((4d0*AShape1+1d0)/4d0*OMEGAM1**4d0)**AShape1)*(1d0/(OMEGA**(4d0*AShape1+1d0)))
     1   *EXP(((4d0*AShape1+1d0)/-4d0)*(OMEGAM1/OMEGA)**4d0)
        SOH2=(Hs2**2)/(4d0*GAM2)*(((4d0*AShape2+1d0)/4d0*OMEGAM2**4d0)**AShape2)*(1d0/(OMEGA**(4d0*AShape2+1d0)))
     1   *EXP(((4d0*AShape2+1d0)/-4d0)*(OMEGAM2/OMEGA)**4d0)
        SOH=SOH1+SOH2
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP1
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
C       Wave Current Inter Action , Hedges et al.(1985)
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SOHU= BB/AA*SOH
        else
        SOHU= SOH
        endif

        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SOHU
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SOHU
        
        
        SUU=((G*Tp1/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SVV=((G*Tp1/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SOHU
        SAUU=((G*Tp1/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SOHU
        SAVV=((G*Tp1/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SOHU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs1/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs1/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        
        Solid_SFUU  = ((1d0)*((Hs1/4)**2d0)*(8/Pi)*SUU)
        Solid_SFUAU = ((1d0)*SAUU)
        
        Solid_SFVV  = ((1d0)*((Hs1/4)**2d0)*(8/Pi)*SVV)
        Solid_SFVAV = ((1d0)*SAVV)
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)

C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 
              
C	NORMAL VELOCITY
C        VX(1) = 1.0D0
C        VX(2) = 0.0D0
C        VX(3) = 0.0D0      
C        VY(1) = 0.0D0
C        VY(2) = 1.0D0
C        VY(3) = 0.0D0     
C        DO IV = 1,3      
C        VX(IV)   = UUX*VX(IV)
C        VY(IV)   = UUY*VY(IV)      
C        VTOL(IV) = VX(IV) + VY(IV)
C        ENDDO
C        CALL VECTORCROSS(VTOL,VR,VZ)
C        CALL VECTORCROSS(VR,VZ,VNOL)      
C        CALL VECTORUNIT(VNOL,VDUM,VV)
        
        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        END SUBROUTINE
        
C     ======================================================================================================           
        Subroutine  Bretschneider_SPectrum(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SFUU,SFVV,SNA,SDG,Tapp,Ucurent)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (i-n)
      DIMENSION SNA(10),SDG(10)
      
      COMMON / WAVESPECTRUMPLOT / NWSPECTRUMPLOT
      
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
      
       ! f=1.0d0
        Uapp = Ucurent
        pi=3.141592653589793d0
        Wm=2d0*pi/Tp
        W=2*pi*f
        SBS=1.25d0/4d0*(Wm**4d0)/(W**-5d0)*(Hs**2d0)*exp(-1.25d0*(Wm/W)**4d0)
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
        AL = 2d0* Pi / WAVENUMBER
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SBSU= BB/AA*SBS
        else
        SBSU= SBS
        endif
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SBSU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SBSU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SBSU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SBSU
      
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
             
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
        
        IF ( NWSPECTRUMPLOT .EQ. 1 ) THEN
        WRITE(*,10) f , SBS , AHH1
!        WRITE(5041,10) f , SBS , AHH1 ! TOEY 10/2021
10      FORMAT(F10.3,2X,F12.5,2X,F12.5,2X,'Bretschneider SPectrum')
        NWSPECTRUMPLOT = 2
        ENDIF
 
        SFUU=sqrt(SFUU*AHH1)
        SFVV=sqrt(SFVV*AHH1)
        
        RETURN
        END SUBROUTINE
      
C     ======================================================================================================           
        Subroutine  Bretschneider_SPectrum_Solid(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SNA,SDG,VNOL,VR,DIAM,
     1                                               PTX,PTY,PTZ,Tapp,Ucurent)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        Uapp = Ucurent 
        pi=3.141592653589793d0
        Wm=2d0*pi/Tp
        W=2*pi*f
        SBS=1.25d0/4d0*(Wm**4d0)/(W**5d0)*(Hs**2d0)*exp(-1.25d0*(Wm/W)**4d0)
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        SBSU= BB/AA*SBS
        else
        SBSU= SBS
        endif

        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SBSU
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SBSU
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SBSU
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SBSU
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SBSU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SBSU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SBSU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SBSU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        
        Solid_SFUU  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SUU)
        Solid_SFUAU = ((1d0)*SAUU)
        
        Solid_SFVV  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SVV)
        Solid_SFVAV = ((1d0)*SAVV)
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        END SUBROUTINE
        
C     ======================================================================================================           
        Subroutine TMA_SPectrum(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SFUU,SFVV,SNA,SDG,Tapp,Ucurent)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION SNA(10),SDG(10)
              
        COMMON / WAVESPECTRUMPLOT / NWSPECTRUMPLOT
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        
        pi=3.141592654d0
        Uapp = Ucurent
        
        
!        call Airywavenumber( d,Tp,wavenumber,G )
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592653589793d0
        OMEGA = 2.0D0*PI/Tp
        CALL NEWTON_RAPHSON (d,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
               
        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
       
        wavenumber=WAVENUMBER
       
        SJS=0d0
        Fp=1/Tp
        !f=0.075d0
        
        if (f.LE.Fp) then
        St=0.07d0
        else
        St=0.09d0
        endif
        
        Alpha=exp(-((f-Fp)**2)/(2d0*St**2d0*Fp**2))
        
        if(Tp/sqrt(Hs).LE.3.6d0)then
        Gamma=5
        elseif(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) then
        Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
        elseif(Tp/sqrt(Hs).GT.5.0d0)then
        Gamma=1d0
        endif
        
        
        
        SJS=0.3125d0*(Hs**2)*(1/Fp)*((f/Fp)**-5d0)*exp(-1.25*(Fp/f)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((f/Fp-1)/st)**2)
        
        
        
        Phi=((cosh(wavenumber*d))**2d0)/((sinh(wavenumber*d))**2d0+ ((2d0*pi*f)**2)*d/G)
        
        !Phi=1.0D0
        
        TMA=SJS*Phi
C       ------------------------------------------------------------        
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
        Pi    = 3.141592653589793d0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        TMAU= BB/AA*TMA
        else
        TMAU= TMA
        endif
        
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*TMAU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*TMAU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*TMAU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*TMAU
      
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
             
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
        
        IF ( NWSPECTRUMPLOT .EQ. 1 ) THEN
        WRITE(*,10) f , TMA , AHH1 , SJS, Phi
!        WRITE(5041,10) f , TMA , AHH1 , SJS, Phi ! TOEY 2021
10      FORMAT(F10.3,2X,F12.5,2X,F12.5,2X,'TMA SPectrum',2X,F12.5,2X,F12.5)
        NWSPECTRUMPLOT = 2
        ENDIF
                
        SFUU=sqrt(SFUU*AHH1)
        SFVV=sqrt(SFVV*AHH1)
        
        RETURN
        END SUBROUTINE
        
C     ======================================================================================================            
        Subroutine  TMA_SPectrum_Solid(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SNA,SDG,VNOL,VR,DIAM,
     1                                               PTX,PTY,PTZ,Tapp,Ucurent)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        pi=3.141592653589793d0
        Uapp = Ucurent
        
!        call Airywavenumber( d,Tp,wavenumber,G )
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/Tp
        CALL NEWTON_RAPHSON(d,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        
        wavenumber=WAVENUMBER
        
        SJS=0d0
        Fp=1/Tp
        !f=0.075d0
        
        if (f.LE.Fp) then
        St=0.07d0
        else
        St=0.09d0
        endif
        
        Alpha=exp(-((f-Fp)**2)/(2d0*St**2d0*Fp**2))
        
        if(Tp/sqrt(Hs).LE.3.6d0)then
        Gamma=5
        elseif(Tp/sqrt(Hs).GT.3.6d0.and.Tp/sqrt(Hs).LE.5.0d0) then
        Gamma=exp(5.75-1.15*Tp/sqrt(Hs))
        elseif(Tp/sqrt(Hs).GT.5.0d0)then
        Gamma=1d0
        endif
        
        
        SJS=0.3125d0*(Hs**2)*(1/Fp)*((f/Fp)**-5d0)*exp(-1.25*(Fp/f)**4)*(1-0.287d0*log(Gamma))*Gamma**exp(-0.5d0*((f/Fp-1)/st)**2)
        
        Phi=((cosh(wavenumber*d))**2d0)/((sinh(wavenumber*d))**2d0+ ((2d0*pi*f)**2)*d/G)
        
        TMA=SJS*Phi
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        if(Uapp.GT.0) Then
C                           CALCULATE APPRENT WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA_app = 2.0D0*PI/Tapp
        CALL NEWTON_RAPHSON(D,G,OMEGA,aKapp)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        BB=OMEGA*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*aKapp*d)))
        AA=2d0*aKapp*(Ucurent+OMEGA/2d0/WAVENUMBER*(1d0+(2d0*WAVENUMBER*d)/(sinh(2d0*WAVENUMBER*d) )  ) )
        TMAU= BB/AA*TMA
        else
        TMAU= TMA
        endif
        
        
        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*TMAU
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*TMAU
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*TMAU
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*TMAU
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*TMAU
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*TMAU
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*TMAU
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*TMAU
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        
        Solid_SFUU  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SUU)
        Solid_SFUAU = ((1d0)*SAUU)
        
        Solid_SFVV  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SVV)
        Solid_SFVAV = ((1d0)*SAVV)
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        end

C     ======================================================================================================            
        Subroutine  User_Define_SPectrum(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SFUU,SFVV,SNA,SDG,Tapp,Ucurent)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(10),SDG(10)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        Uapp = Ucurent
        Fp=1/Tp
 
        
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592653589793d0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        
        Pi=3.141592654
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SUSER
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SUSER
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SUSER
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SUSER
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        SFUU=sqrt(SFUU*AHH1)
        SFVV=sqrt(SFVV*AHH1)
        
        RETURN
        end

C     ======================================================================================================            
        !*********************************************************
        !For Solid Element
        !**********************************************************
        
        Subroutine  User_Define_SPectrum_SOILD(Hs,Tp,G,d,ROH,Cd,CI,Dia,Z,f,SNA,SDG,VNOL,VR,DIAM,
     1                                               PTX,PTY,PTZ,Tapp,Ucurent)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        
       !-----------------------------------------------------------
       ! Definition
       !-----------------
       ! Input
       !-------
       ! Hs = wave height
       ! Tp = wave period
       ! G  = Gravity Acc.
       ! d  = water depth
       ! ROH = water density
       ! Cd  = Dreg Coefficient
       ! CI  = Innertia Coefficient
       ! Dia = Section Diameter
       ! Z   = Vertical Elevation of wave force
       ! f   = abitary frequency
       ! SNA = NATURAL PERIOD OF STRUCTURE
       ! SDG = STRUCTURE DAMPING
       !----------
       ! Output
       !--------
       ! SFUU = forec in horizontial direction
       ! SFVV = force in vertical direction
       !----------------------------------------------------------
        Uapp = Ucurent
        Fp=1/Tp
 
C       ----------------------------------------------------------
C                           CALCULATE WAVE NUMBER
C       ----------------------------------------------------------
        Pi    = 3.141592654D0
        OMEGA = 2.0D0*PI/TP
        CALL NEWTON_RAPHSON(D,G,OMEGA,WAVENUMBER)
C       Call Airywavenumber( d  , Tp , WAVENUMBER , G )
C       ----------------------------------------------------------
        
        Pi=3.141592653589793d0
        
        AL = 2* Pi / WAVENUMBER
        
        SUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SAUUX=((2*Pi*Fp*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SUSER
        SAVVX=((2*Pi*Fp*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2)*SUSER
        
        
        SUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*SUSER
        SAUU=((G*Tp/AL*cosh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SUSER
        SAVV=((G*Tp/AL*sinh(WAVENUMBER*Z)/cosh(WAVENUMBER*d))**2d0)*((2*Pi*Fp)**2d0)*SUSER
        
        SFUU=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SUU+((ROH*CI*Pi/4d0*Dia**2)**2)*SAUU
        SFVV=((0.5d0*ROH*Cd*Dia)**2)*((Hs/4)**2d0)*(8/Pi)*SVV+((ROH*CI*Pi/4d0*Dia**2)**2)*SAVV
        
        
        Solid_SFUU  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SUU)
        Solid_SFUAU = ((1d0)*SAUU)
        
        Solid_SFVV  = ((1d0)*((Hs/4)**2d0)*(8/Pi)*SVV)
        Solid_SFVAV = ((1d0)*SAVV)
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        end
        
C     ======================================================================================================            
        !**************************************************************************
        ! Wind Spectrum
        !**************************************************************************
        Subroutine Kaimal_Spectrum(zo,Alpha,Umean,RHIGH,G,ROHAIR,Cd,Dia,Z,f,SFUU,SNA,FWIND)

        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION SNA(10),SDG(10)
        DIMENSION FWIND(3)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo         = Foughness parameter(m)
        ! Alpha      = Power-Law Exponent
        ! U10        = Mean Wind Velocity
        ! AIo        = Turbulence Model
        ! ALu        = Integral Length Scale  
        ! SD         = Standaud Deviation = U10 * Lu 
        ! FWIND(1)   = Turbulence Model FACTOR
        ! FWIND(2)   = Integral Length Scale FACTOR
        ! FWIND(3)   = Standaud Deviation FACTOR
        !*********************************************
        FN = SNA(1)
        
        IF (Z.eq.0D0)THEN
        ALu=300D0
        ELSE
        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        ENDIF
        AIo = 1d0/log(z/zo)*FWIND(2)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        IF (U.eq.0D0)THEN
        SKaimal=0D0
        ELSE
        SKaimal=(SD**2d0)*6.868d0*(ALu/U*f)/(1d0+10.32d0*f*ALu/U)**(5d0/3d0) 
        ENDIF
      !  SKaimal=(SD**2d0)*6.868d0*(ALu/U*f)/(1d0+10.32d0*f*ALu/U)**(5d0/3d0) 
        
        SFF=((0.5d0*ROHAIR*(U**2d0)*Dia*Cd)**2d0)*SKaimal
        
        !**********************
        ! Transfer function
        !**********************
 
        AfRatio=f/fn

        AH= 1d0/( (1-AfRatio**2d0)**2d0+(2d0*DampWind*AfRatio)**2d0)

        SFUU=sqrt(AH*SFF)
        
        end

C     ======================================================================================================            
        Subroutine Kaimal_Spectrum_Solid(zo,Alpha,Umean,RHIGH,G,ROH,Cd,Dia,Z,f,SNA,FWIND,SDG,VNOL,VR,DIAM,PTX,PTY,PTZ)
                                     
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        DIMENSION FWIND(3)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo         = Foughness parameter(m)
        ! Alpha      = Power-Law Exponent
        ! U10        = Mean Wind Velocity
        ! AIo        = Turbulence Model
        ! ALu        = Integral Length Scale  
        ! SD         = Standaud Deviation = U10 * Lu 
        ! FWIND(1)   = Turbulence Model FACTOR
        ! FWIND(2)   = Integral Length Scale FACTOR
        ! FWIND(3)   = Standaud Deviation FACTOR
        !*********************************************
        FN = SNA(1)
        
        IF (Z.eq.0D0)THEN
        ALu=300D0
        ELSE
        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        ENDIF
        AIo = 1d0/log(z/zo)*FWIND(2)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        IF (U.eq.0D0)THEN
        SKaimal=0D0
        ELSE
        SKaimal=(SD**2d0)*6.868d0*(ALu/U*f)/(1d0+10.32d0*f*ALu/U)**(5d0/3d0) 
        ENDIF
      !  SKaimal=(SD**2d0)*6.868d0*(ALu/U*f)/(1d0+10.32d0*f*ALu/U)**(5d0/3d0) 
        
      !  SFF=((0.5d0*ROHAIR*(U**2d0)*Dia*Cd)**2d0)*SKaimal
        
        
        Solid_SFUU  = SQRT( (U**2d0)* SKaimal )                     
        Solid_SFUAU = 0D0 
        
        Solid_SFVV  = 0D0 
        Solid_SFVAV = 0D0 
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        end

C     ======================================================================================================            
        Subroutine Kaimal_IEC_Spectrum(zo,Alpha,Umean,RHIGH,G,ROHAIR,Cd,Dia,Z,f,SFUU,SNA,FWIND)
        !IEC 61400-1 page 70
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION FWIND(3)
        DIMENSION SNA(10),SDG(10)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
        FN = SNA(1)
        
        if (z.LE.60d0)then
        ALu = 8.1d0*(0.7d0*z )*FWIND(1)
        else
        ALu = 8.1d0*(42d0)*FWIND(2)
        endif
        
        
        AIo = 1d0/log(z/zo)
        
        ! Calculation Wind Profile
        
        !Z=10D0
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
       ! SKaimal_IEC=(SD**2d0)*6.868d0*(ALu/U*f)/(1d0+10.32d0*f*ALu/U)**(5d0/3d0) 
       IF (U.EQ.0D0) THEN
        SKaimal_IEC=0D0
       ELSE       
        SKaimal_IEC=(SD**2d0)*6.868d0*(ALu/U*f)/(1d0+10.32d0*f*ALu/U)**(5d0/3d0)
        ENDIF
        
        SFF=((0.5d0*ROHAIR*(U**2d0)*Dia*Cd)**2d0)*SKaimal_IEC
        
        !**********************
        ! Transfer function
        !**********************
 
        AfRatio=f/fn

        AH= 1d0/( (1-AfRatio**2d0)**2d0+(2d0*DampWind*AfRatio)**2d0)

        SFUU=sqrt(AH*SFF)

        RETURN
        END SUBROUTINE
        
C     ======================================================================================================            
        Subroutine Kaimal_IEC_Spectrum_Solid(zo,Alpha,Umean,RHIGH,G,ROH,Cd,Dia,Z,f,SNA,FWIND,SDG,VNOL,VR,DIAM,PTX,PTY,PTZ)
        !IEC 61400-1 page 70
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        DIMENSION FWIND(3)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
        FN = SNA(1)
        
        if (z.LE.60d0)then
        ALu = 8.1d0*(0.7d0*z )*FWIND(1)
        else
        ALu = 8.1d0*(42d0)*FWIND(2)
        endif
        
        
        AIo = 1d0/log(z/zo)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        SKaimal_IEC=(SD**2d0)*6.868d0*(ALu/U*f)/(1d0+10.32d0*f*ALu/U)**(5d0/3d0) 
        
       
        Solid_SFUU  = SQRT( (U**2d0)* SKaimal_IEC )                     
        Solid_SFUAU = 0D0 
        
        Solid_SFVV  = 0D0 
        Solid_SFVAV = 0D0 
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN

        END SUBROUTINE
        
        
C     ======================================================================================================           
        Subroutine Von_Karman_Spectrum(zo,Alpha,Umean,RHIGH,G,ROHAIR,Cd,Dia,Z,f,SFUU,SNA,FWIND)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION FWIND(3)
        DIMENSION SNA(10),SDG(10)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
         FN = SNA(1)
        
        IF (Z.eq.0D0)THEN
        ALu=300D0
        ELSE
        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        ENDIF
c        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1) 
        AIo = 1d0/log(z/zo)*FWIND(2)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        
        IF(U.EQ.0D0)THEN
        Svk=0D0
        ELSE
        Svk= (SD**2d0)*4.0d0*(ALu/U*f)/(1d0+70.8d0*(f*ALu/U)**2d0)**(5d0/6d0)
        ENDIF
        
        SFF=((0.5d0*ROHAIR*(U**2d0)*Dia*Cd)**2d0)*Svk
        
        !**********************
        ! Transfer function
        !**********************
 
        AfRatio=f/fn

        AH= 1d0/( (1-AfRatio**2d0)**2d0+(2d0*DampWind*AfRatio)**2d0)

        SFUU=sqrt(AH*SFF)
        
        
        end

C     ======================================================================================================            
        Subroutine Von_Karman_Spectrum_Solid(zo,Alpha,Umean,RHIGH,G,ROH,Cd,Dia,Z,f,SNA,FWIND,SDG,VNOL,VR,DIAM,PTX,PTY,PTZ)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        DIMENSION FWIND(3)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
         FN = SNA(1)
        
        IF (Z.eq.0D0)THEN
        ALu=300D0
        ELSE
        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        ENDIF
c        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1) 
        AIo = 1d0/log(z/zo)*FWIND(2)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        
        IF(U.EQ.0D0)THEN
        Svk=0D0
        ELSE
        Svk= (SD**2d0)*4.0d0*(ALu/U*f)/(1d0+70.8d0*(f*ALu/U)**2d0)**(5d0/6d0)
        ENDIF
        
        Solid_SFUU  = SQRT( (U**2d0)* Svk )                     
        Solid_SFUAU = 0D0 
        
        Solid_SFVV  = 0D0 
        Solid_SFVAV = 0D0 
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        
        end
        

C     ======================================================================================================            
        Subroutine Davenport_Spectrum(zo,Alpha,Umean,RHIGH,G,ROHAIR,Cd,Dia,Z,f,SFUU,SNA,FWIND)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION FWIND(3)
        DIMENSION SNA(10),SDG(10)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
        FN = SNA(1)
                IF (Z.eq.0D0)THEN
        ALu=300D0
        ELSE
        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        ENDIF
c        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        AIo = 1d0/log(z/zo)*FWIND(2)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        IF(U.EQ.0D0)THEN
        SDavenport=0D0
        ELSE
        SDavenport = (SD**2d0)*2d0/3d0*((1200d0/U*f)**2d0)/(1d0+(f*1200d0/U)**2d0)**(4d0/3d0)
        ENDIF
        SFF=((0.5d0*ROHAIR*(U**2d0)*Dia*Cd)**2d0)*SDavenport
        
        !**********************
        ! Transfer function
        !**********************
 
        AfRatio=f/fn

        AH= 1d0/( (1-AfRatio**2d0)**2d0+(2d0*DampWind*AfRatio)**2d0)

        SFUU=sqrt(AH*SFF)
        
        
        end

C     ======================================================================================================            
        Subroutine Davenport_Spectrum_Solid(zo,Alpha,Umean,RHIGH,G,ROH,Cd,Dia,Z,f,SNA,FWIND,SDG,VNOL,VR,DIAM,PTX,PTY,PTZ)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        DIMENSION FWIND(3)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
        FN = SNA(1)
                IF (Z.eq.0D0)THEN
        ALu=300D0
        ELSE
        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        ENDIF
c        ALu = 300d0*(z/300d0)**(0.46d0+0.074*log(zo))*FWIND(1)
        AIo = 1d0/log(z/zo)*FWIND(2)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        IF(U.EQ.0D0)THEN
        SDavenport=0D0
        ELSE
        SDavenport = (SD**2d0)*2d0/3d0*((1200d0/U*f)**2d0)/(1d0+(f*1200d0/U)**2d0)**(4d0/3d0)
        ENDIF
         
        Solid_SFUU  = SQRT( (U**2d0)* SDavenport )                     
        Solid_SFUAU = 0D0 
        
        Solid_SFVV  = 0D0 
        Solid_SFVAV = 0D0 
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        
        end

C     ======================================================================================================            
        Subroutine Eurocode_Spectrum(zo,Alpha,Umean,RHIGH,G,ROHAIR,Cd,Dia,Z,f,SFUU,SNA,FWIND)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION FWIND(3)
        DIMENSION SNA(10),SDG(10)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
        
        Zmin=1d0 ! for Minimum Height Define  Eurocode part 4 page 20
                 ! Table 4.1: Sea of Coastal area exposed to the open sea
                 
        Ak = 1d0 ! is the turbulence factor
        C0 = 1d0 ! is the orgraphy factor 
                 ! Eurocode part 4 page 22
        FN = SNA(1)     
        
        if (z.GT.Zmin)THEN
        ALu = 300d0*(z/300d0)**(0.67d0+0.05*log(zo))*FWIND(2) 
        else
        ALu = 300d0*(Zmin/300d0)**(0.67d0+0.05*log(zo))*FWIND(2)
        endif
        
        AIo = Ak/(C0*log(z/zo))*FWIND(1)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        IF(U.EQ.0D0)THEN
        SEurocode=0D0
        ELSE
        SEurocode=(SD**2d0)*6.8d0*(ALu/U*f)/(1d0+10.2d0*f*ALu/U)**(5d0/3d0) 
        ENDIF
        SFF=((0.5d0*ROHAIR*(U**2d0)*Dia*Cd)**2d0)*SEurocode
        
        !**********************
        ! Transfer function
        !**********************
 
        AfRatio=f/fn

        AH= 1d0/( (1-AfRatio**2d0)**2d0+(2d0*DampWind*AfRatio)**2d0)

        SFUU=sqrt(AH*SFF)
        
        
        end
        
        
        Subroutine Eurocode_Spectrum_Solid(zo,Alpha,Umean,RHIGH,G,ROH,Cd,Dia,Z,f,SNA,FWIND,SDG,VNOL,VR,DIAM,PTX,PTY,PTZ)
        
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION SNA(1),SDG(1)
        DIMENSION VX(3),VY(3),VZ(3),VTOL(3),VNOL(3),VT(3),VR(3)
        DIMENSION AX(3),AY(3),ATOL(3),ANOL(3)
        DIMENSION FWIND(3)
        
        !*********************************************
        ! Roudhness Length (Zo)
        !*********************************************
        ! zo = Foughness parameter(m)
        ! Alpha = Power-Law Exponent
        ! U10 =  Mean Wind Velocity
        ! AIo = Turbulence Model
        ! ALu = Integral Length Scale  
        ! SD = Standaud Deviation = U10 * Lu 
        !*********************************************
        
        Zmin=1d0 ! for Minimum Height Define  Eurocode part 4 page 20
                 ! Table 4.1: Sea of Coastal area exposed to the open sea
                 
        Ak = 1d0 ! is the turbulence factor
        C0 = 1d0 ! is the orgraphy factor 
                 ! Eurocode part 4 page 22
        FN = SNA(1)     
        
        if (z.GT.Zmin)THEN
        ALu = 300d0*(z/300d0)**(0.67d0+0.05*log(zo))*FWIND(2) 
        else
        ALu = 300d0*(Zmin/300d0)**(0.67d0+0.05*log(zo))*FWIND(2)
        endif
        
        AIo = Ak/(C0*log(z/zo))*FWIND(1)
        
        ! Calculation Wind Profile
        
        if (g.LT.12d0)then
        U=Umean*(z/RHIGH)**Alpha
        else
        U=Umean*(z/RHIGH)**Alpha
        endif
        SD = AIo*U*FWIND(3)
        
        IF(U.EQ.0D0)THEN
        SEurocode=0D0
        ELSE
        SEurocode=(SD**2d0)*6.8d0*(ALu/U*f)/(1d0+10.2d0*f*ALu/U)**(5d0/3d0) 
        ENDIF
         
        Solid_SFUU  = SQRT( (U**2d0)* SEurocode )                     
        Solid_SFUAU = 0D0 
        
        Solid_SFVV  = 0D0 
        Solid_SFVAV = 0D0 
        
        
        ! Transfer Function
        
        Fn1=SNA(1) !1D0/0.0638D0 !EIGVPeriod(1)
        !Fn2=SNA(2)
        
        DAM1=SDG(1)
        !DAM2=SDG(2)
        
        FreRATIO1=f/Fn1
       ! FreRATIO2=f/Fn2
        
        AHH1 = 1d0/((1-FreRATIO1**2d0)**2d0+ (2d0*DAM1*FreRATIO1)**2d0)
       ! AHH2 = 1d0/((1-FreRATIO2**2d0)**2d0+ (2d0*DAM2*FreRATIO2)**2d0)
        
       ! AHH = AHH1 + AHH2
 
        UUX=sqrt(Solid_SFUU*AHH1)
        UUY=sqrt(Solid_SFUAU*AHH1)
        
        AAX=sqrt(Solid_SFVV*AHH1)
        AAY=sqrt(Solid_SFVAV*AHH1)
        
C     MAGNITUDE OF VELOCITY
        VV = SQRT(UUX*UUX + UUY*UUY + UUZ*UUZ) 

        ! FORCE DATA BY TOEY 
        VNOL(1) = UUX
        VNOL(2) = UUY
        VNOL(3) = 0.0D0
              
        UN = VNOL(1)
        VN = VNOL(2)  
        WN = VNOL(3)  
        
C	DRAG FORCE 
        STD = 0.50D0*ROH*CD            
        PDX = STD*VV*UN  
        PDY = STD*VV*VN
        PDZ = STD*VV*WN 
        
C	INNERTIA FORCE 
        STI = ROH*CI*PI*DIAM/4.0D0  
        PIX = STI*AAX
        PIY = STI*AAY
        PIZ = STI*AAZ  
        
        PTX = ( PDX + PIX )
        PTY = ( PDY + PIY )
        PTZ = ( PDZ + PIZ )
        
        RETURN
        
        
        end
        
        
C     ======================================================================================================    	  
	  SUBROUTINE GAMMA_FUNCTION(X,GA)
C       Purpose: Compute the gamma function ?(x)
C       Input :  x  --- Argument of ?(x)
C                       ( x is not equal to 0,-1,-2,??? )
C       Output:  GA --- ?(x)
C       ==================================================
C
        !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (i-n)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     1      -0.6558780715202538D0, -0.420026350340952D-1,
     1          0.1665386113822915D0,-.421977345555443D-1,
     1          -.96219715278770D-2, .72189432466630D-2,
     1          -.11651675918591D-2, -.2152416741149D-3,
     1          .1280502823882D-3, -.201348547807D-4,
     1          -.12504934821D-5, .11330272320D-5,
     1          -.2056338417D-6, .61160950D-8,
     1          .50020075D-8, -.11812746D-8,
     1          .1043427D-9, .77823D-11,
     1          -.36968D-11, .51D-12,
     1          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
C     ======================================================================================================             
        
      Subroutine  PMSpectrum_Time_domain
  
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
   
      
      COMMON /ChanaTImeSeries/ DT 
      
      COMMON/ Spectrum_Time_domain / AmSpm(100000),AmSjs(100000),fcal(100000),Phi(100000)
      COMMON/ Spectrum_Time_domain2 / Loop_Spectrum
       
      COMMON/chana_OFFSHORE/ Irandom 
       
     	COMMON /OFFSHORR/ SEABED,WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,
	1                  RWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
	1                  VWIND0,VTIDE,H0,
	1                  AP,SP,CS,
	1                  RHIGH,UH,ALPHA,Z0,
	1                  VG(3),V1(3),V2(3),VW(3),WVZETA,WVTIME,PEAKWLEV,
     1                  POWERLAW,VCURRENTL(5),VCURRENTAPI,UHAPI,NCURRENT,NWIND,
     1                  FACTOR,VCURRENTP,AVERAGE,UHD,
     1                  WKF,CBF,WFC,
     1                  FSpectrumStart,FSpectrumEnd
      COMMON /OFFSHORI/ NUMBEROFRANDOMWANVEXInteration
       
       if ( Irandom == 1) Then 
       
      
      OPEN(UNIT=987654  ,FILE='Spectrum.txt'     ,STATUS='UNKNOWN')
      OPEN(UNIT=987655  ,FILE='SpectrumTime.txt' ,STATUS='UNKNOWN')
      
      
      Hw = WVHIGHT
      Tw = PERIOD
      Pi = 3.141592654d0
      
      time = 0.5
      
      call OFFSHSTEP(TIME,ITIME,NTIME,'CALL')
      
      
      call OFFSHSTEP(TIME,ITIME,NTIME,'CALT')
      
       fstart = FSpectrumStart
       fend   = FSpectrumEnd
       interation = NUMBEROFRANDOMWANVEXInteration
       
       Loop_Spectrum =  interation
        
        
            fint = ( fend - fstart ) / interation
              fp = 1.0d0/Tw
              Tp = Tw
               
        chana = 234d0
        
            f = 0.05d0
        
        Do i=1,interation
        
        CALL RANDOM_NUMBER(random) 
        
        fi = fi + fint 
        f =  fi + fint * random
        
        fcal(i) = f
        
        !***************************
        ! Pierson-Moskowitz
        !***************************
        
        Spm =  0.3125d0 *(Hw**2)* (fp**4) * (f**-5) *exp(-1.25d0*((fp/f)**4))
        
        !***************************
        ! Jonswap Parameter
        !***************************   
        
        if (f <= fp )then
            Sigma = 0.07d0
        elseif ( f > fp ) then
            Sigma = 0.09d0
        endif
        
        Alpha = exp( - ( ( f - fp)**2 ) / (2 * ( Sigma * fp ) **2 ) )
        
        Gammafactor = Tp / sqrt( Hw )
        
        if ( Gammafactor <= 3.6d0 ) Then
        
        Gramma = 5d0
            
        elseif ( Gammafactor > 3.6d0 .AND. Gammafactor < 5.0d0 ) Then
            
        Gramma = exp(5.75d0 - 1.15d0 * Gammafactor )
            
        elseif ( Gammafactor > 5.0d0 ) Then
          
        Gramma = 1d0    
            
        endif
               
        Cjs = 1.0d0 - 0.287d0 * LOG(Gramma)
       
        Sjs =  Cjs * Spm * ( Gramma**Alpha )     
        
        CALL RANDOM_NUMBER(random) 
    
        Phi(i)  = 2d0 * Pi *( random - 0.5d0 )       
        AmSpm(i) = sqrt( 2d0 * Spm * fint )
        AmSjs(i) = sqrt( 2d0 * Sjs * fint )

        !write(*,*)  f ,Spm , Sjs
        write(987654,*)  f ,Spm , Sjs

        ENDDO    
        
       
       Tstart = 0 
       Tinterval = 1d0/SFR  
       Ttotal =  ( Tend - Tstart )* SFR 
       Time = Tstart 
        
   ! call OFFSHSTEP(TIME,ITIME,NTIME,'CALT')  
   !   
   ! Time = DT
   ! 
   ! Do j= 1, NTIME
   !     
   !     chana = 20d0
   !     
   !       call Generate_Spectrum_Time( Time,0,0,0, SumWaveSpmTCos , SumWaveSpmTSin , SumWaveSjsTCos , SumWaveSjsTSin  )
   !     
   !               
   !     
   !      write(987655,*) Time , SumWaveSpmTCos , SumWaveSpmTSin 
   !      
   !     Time = Time + DT
   !     
   ! Enddo   
   ! 
   !  stop
       
        Irandom = 2d0
       
       Endif
      
      return

      end subroutine 
      
!================================================================================================================      
!================================================================================================================      
!================================================================================================================      
                  
            
      Subroutine  Generate_Spectrum_Time( Time,RK,X,OMEGA, SumWaveSpmTCos , SumWaveSpmTSin , SumWaveSjsTCos , SumWaveSjsTSin  )
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /ChanaTImeSeries/ DT 
      
      COMMON/ Spectrum_Time_domain /  AmSpm(100000),AmSjs(100000),fcal(100000),Phi(100000) 
      
      COMMON/ Spectrum_Time_domain2 / Loop_Spectrum
      
!       SKL  = 80d0    ! use scale model to reduce time consume in calculation..!!
!       Tend = 1560d0  ! example : about 3 hours for model scale 1:50
!       SFR  = 25d0    ! Sample Frequency (Hz)   
!       
!
!       
!       call OFFSHSTEP(TIME,ITIME,NTIME,'CALT')
      
       Pi = 3.141592654d0
       
!       Do j= 1, Ttotal
!           
!           WaveSpmT = 0d0 
!           WaveSjsT = 0d0 
!           
!           Do k = 1,interation 
!               
!               WaveSpm = AmSpm(k) * cos(  ( 2d0*Pi*fcal(k) ) * Time - Phi(k) )     
!               WaveSjs = AmSjs(k) * cos(  ( 2d0*Pi*fcal(k) ) * Time - Phi(k) )    
!               
!               WaveSpmT = WaveSpmT + WaveSpm 
!               
!               WaveSjsT = WaveSjsT + WaveSjs
!               
!           Enddo     
!           
!            write(987655,*) Time , WaveSpmT , WaveSjsT 
!            
!           Time = Time + Tinterval
!           
!       Enddo   
       
       
           
           SumWaveSpmTCos = 0d0
           SumWaveSpmTSin = 0d0
           
           SumWaveSjsTCos = 0d0
           SumWaveSjsTSin = 0d0
           
           Do k = 1,Loop_Spectrum
               
               WaveSpmCos = AmSpm(k) * cos( RK*X -   ( 2d0*Pi*fcal(k) ) * Time + Phi(k) )
               WaveSpmSin = AmSpm(k) * sin( RK*X -   ( 2d0*Pi*fcal(k) ) * Time + Phi(k) )
               
               WaveSjsCos = AmSjs(k) * cos( RK*X -  ( 2d0*Pi*fcal(k) ) * Time + Phi(k) )
               WaveSjsSin = AmSjs(k) * sin( RK*X -  ( 2d0*Pi*fcal(k) ) * Time + Phi(k) )
               
               SumWaveSpmTCos = SumWaveSpmTCos + WaveSpmCos
               SumWaveSpmTSin = SumWaveSpmTSin + WaveSpmSin
               
               SumWaveSjsTCos = SumWaveSjsTCos + WaveSjsCos
               SumWaveSjsTSin = SumWaveSjsTSin + WaveSjsSin
               
               
           Enddo     
           


      
      return
            
      end  subroutine 

      
!================================================================================================================      
!================================================================================================================      
!================================================================================================================      
                        
      
      
      