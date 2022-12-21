C	=======================================================================
C	=======================OFFSHORE LOAD PARAMETER=========================
C	=======================================================================
      SUBROUTINE OFFSPARA 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

	COMMON /MGRAV/ NGRAV
	
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      
    	COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
	
	COMMON /OFFSHORE/ SEABED,WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,
	1                  RWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
	1                  VWIND0,VTIDE,H0,
	1                  AP,SP,CS,VB(3),HM,HW,HC,
	1                  RHIGH,UH,ALPHA,Z0,
	1                  VG(3),V1(3),V2(3),VW(3),WVZETA,WVTIME,PEAKWLEV,
     1                  POWERLAW,VCURRENTL,VCURRENTAPI,UHAPI,NCURRENT,NWIND
     1                  CHARNC,AKAMAN,AVERAGE,UHD,NROUGHNESS,								   
     1                  FSpectrumStart,FSpectrumEnd,NUMBEROFRANDOMWANVEXInteration,
     1                  TIME
	      
C	=======================================================================
      
      CALL  OFFSPARA_READ
      
      CALL  OFFSPARA_READ_COLLECT_DATA

      IF (ICONTROLSPEC.NE.1) CALL OFFSHSTEP(TIME,ITIME,NTIME,'READ') 
      IF (ICONTROLSPEC.EQ.1) CALL OFFSHSTEP(TIME,ITIME,NTIME,'REDD') 
      
      CALL STRTALIGN 
      
      IF (OPRES.EQ.2) CALL READGRAPH
       
      RETURN
      END SUBROUTINE
C	=======================================================================
      
C	=======================================================================
C	================  OFFSHORE PARAMETER READ DATA   ======================
C	=======================================================================
C     --- READ DATA FUNCTION FROM DATA INPUT ---
C     INPUT  -- READ OFFSHORE PARAMETER FROM DATA INPUT 
C     OUTPUT -- COMMON /OFFSHOREx ( ALL )
      SUBROUTINE OFFSPARA_READ
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*20 NAME_SCATTER,NAME_HEIGHT,NAME_PERIOD,NAME_MAGNI
      CHARACTER*20 NAME_READ
      CHARACTER*200 NAME_PARAMETER,NAME_SPECTRUM
      CHARACTER*250 NAME
C	=======================================================================
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /MGRAV/ NGRAV
	COMMON /OFFSHOREx/ SEABEDx(1000),WVHIGHTx(1000),WDEPTHx(1000),THIGHTx(1000),H1POSx(1000),H2POSx(1000),
	1                  RWAVEx(1000),ORDERx(1000),PERIODx(1000),GRAVx(1000),RHOWx(1000),RHOAx(1000),
	1                  VWIND0x(1000),VTIDEx(1000),H0x(1000),
	1                  APx(1000),SPx(1000),CSx(1000),HMX(1000),HWX(1000),HCX(1000),VBREAKINGX1(1000),VBREAKINGX2(1000),
     1                  VBREAKINGX3(1000),IBREAKING(1000),
	1                  RHIGHx(1000),UHx(1000),ALPHAx(1000),Z0x(1000),
     1                  VWAVEx1(1000),VWAVEx2(1000),VWAVEx3(1000), VWINDx1(1000),VWINDx2(1000),VWINDx3(1000),
     1                  NCURRENTX(1000),POWERLAWX(1000),VCURRENTLX(1000,5),VCURRENTAPIX(1000),UHAPIX(1000),NWINDX(1000),
     1                  FACTORX(1000),VCURRENTPX(1000),AVERAGEX(1000),UHDX(1000),
     1                  WKFX(1000),CBFX(1000),WFCX(1000),
     1                  Offshoreparameterx(1000),
     1                  FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000),
     1                  TIMEX(1000)

      

      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      COMMON /FATIGUE_OFFSHORE/ NWAVESCAT, NWINDSCAT
      COMMON /WAVE_SCATTER/ NUMSCAT(1000),HSCAT(1000),PERIODSCAT(1000),AMENITUDESCAT(1000)
      
      COMMON /NAME_OFFSHORE_PARAMETER/ NAME_PARAMETER(1000),NAME_SPECTRUM(1000)
      
C	=======================================================================
      COMMON /XVWAVE/ VWAVE(3),VWIND(3)
      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
      COMMON /SEA_PARA/ NSPECTRUM
      COMMON /stroke_k_and_Ramda/RKx(1000),RAMDAx(1000),ACOEFFICIENT(1000)
      DIMENSION WAVE(3),WIND(3)
      DIMENSION BREAKING_VE(3)
C	=======================================================================  
      ! MODIFY 5-01-2013
      ! OCCURRENTX(I)     =  TYPES OF CURRENT PROFILE
      ! POWERLAWX(I)      =  EXPONENT POWER-LAW MODEL FOR WIND SPEED PROFILE
      ! VCURRENTLX(I)     =  VELOCITY LINEAR CURRENT
      ! VCURRENTAPIX(I)   =  CURRENT VELOCITY (API)
      ! VWIND0x(i)        =  WIND-GENERATE CURRENT VELOCITY (DNV)
      ! VWIND0X(I)        =  WIND GENERATE CURRENT
      ! VTIDEX(I)         =  TIDE CURRENT AT STILL WATER LEVEL
      ! UHX(I)            =  MEAN WIND SPEED ( DNV ) >> 10 MIN
      ! UHAPIX(I)         =  MEAN WIND SPEED ( API ) >> 1  HOUR
      ! NWIND(I)          =  TYPE OF WIND LOAD
      ! ROUGHNESSX(I)     =  ROUGHNESS PARAMETER FOR USER DEFINED 
      ! AVERAGEX(I)       =  AVERAGING TIME (API)
      ! UHDX(I)           =  WIND VELOCITY ( USER DEFINED )
      ! FACTORX(I)        =  WIND-GENERATE CURRENT FACTOR (K) SEE DNV OS-J101 JULY 2011 PAGE 56. AND 61400-3 IEC:2009
      ! VCURRENTPX(I)     =  CURRENT VELOCITY ON POWER-LAW PROFILE
      ! Z0X(I)            =  REFERENCE HIGHT FOR WIND
      ! H0X(I)            =  REFERENCE DEPTH FOR WIND-GENERATE CURRENT ( DNV = 50 m , IEC = 20 m )
      ! ALPHAX(I)         =  WIND PROFILE PARAMETER
      ! RWAVEX(I)         =  WAVE OPTION 
      !                      1 >> AIRY WAVE THEORY
      !                      2 >> STROKE'S FIFTH ORTHEORY
      !                      3 >> STREAM FUNCTION
      ! GRAVX(I)          =  GRAVITY ACCLERATION 
      ! RHOWx(i)          =  WATER DENSITY
      ! RHOAx(i)          =  AIR DENSITY
      ! ORDERX(I)         =  STREMA FUNCTION ORDER 
      ! WAVEANGLE         =  WAVE ANGLE PROPAGATE X-DIRECTION
      ! WINDANGLE         =  WIND ANGLE PROPAGATE X-DIRECTION
     
      IF(KFTTD.EQ.1) THEN
      WRITE (*,*) "FATIGUE ANALYSYS NOT AVAILABLE IN THIS VERSION"
      STOP
      IF (KFATL.EQ.2)THEN
      WRITE (*,*) "FATIGUE ANALYSYS NOT AVAILABLE IN THIS VERSION"
      STOP
      ENDIF
          
          
      READ(ITI,*) 
      READ(ITI,*) NUM_OF_OFFSHORE_PARAMETER,NSPECTRUM
      

      DO i=1,NUM_OF_OFFSHORE_PARAMETER 
      Offshoreparameterx(i) = i
      
      ! PREVENT ERROR EFFECT ON THIGHTX(I)
      THIGHTX(I)=0.0D0
      
      READ(ITI,*) !---------------------------------------
      READ(ITI,'(A)') NAME(1:250)!OFFSHORE LOAD PARAMETERS CASE ---
      NAME_PARAMETER(I) = NAME(34:250)
      READ(ITI,*) !---------------------------------------
      READ(ITI,*) !---------------------------- WAVE PARAMETERS------------------------------
      
      READ(ITI,*) SEABEDx(i) ,NAME_READ,NAME_READ,NAME_READ       
      READ(ITI,*) WVHIGHTx(i),NAME_READ,NAME_READ,NAME_READ       
      READ(ITI,*) WDEPTHx(i) ,NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) TIMEX(I)   !
      READ(ITI,*) H1POSx(i)  ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ  
      READ(ITI,*) H2POSx(i)  ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ  
      READ(ITI,*) RWAVEx(i)  ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) ORDERx(i)  ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) PERIODx(i) ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) GRAVx(i)   ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) RHOWx(i)   ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) RHOAx(i)   ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) WKFX(i)    ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) WFCX(i)    ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*)!---------------------------- CURRENT PARAMETERS-----------------------------
      READ(ITI,*) NCURRENTX(I) ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) POWERLAWX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VWIND0x(i)   ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VTIDEX(I)    ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VCURRENTPX(I),NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) H0X(I)       ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) CBFX(I)      ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VCURRENTLX(I,1),NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VCURRENTLX(I,2),NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VCURRENTLX(I,3),NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VCURRENTLX(I,4),NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) VCURRENTLX(I,5),NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      
C     ================================================================================      
C        MUTIPLY FACTOR ( SEE DNV OS-J101 JULY 2011 PAGE 56. AND 61400-3 IEC:2009 )
C     ================================================================================
      READ(ITI,*)!------------------------------- BREAKING PARAMETERS-----------------------------
      !READ(ITI,*) IBREAKING(I)
      READ(ITI,*) APx(i)
      READ(ITI,*) SPx(i)
      READ(ITI,*) CSx(i)
      !READ(ITI,*) HMX(I)
      !READ(ITI,*) HWX(I)
      !READ(ITI,*) HCX(I)
      !READ(ITI,*) BREAKING
      IBREAKING(I) = 0.0d0
      HMX(I)    = 0.0d0
      HWX(I)    = 0.0d0
      HCX(I)    = 0.0d0
      BREAKING  = 0.0d0
      READ(ITI,*)!------------------------------- WIND PARAMETERS-----------------------------
      READ(ITI,*) NWINDX(I) ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) Z0X(I)    ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) UHX(I)    ,NAME_READ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) RHIGHx(I) ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) ALPHAX(I) ,NAME_READ,NAME_READ,NAME_READ,NAME_READ
      READ(ITI,*) UHAPIX(I) ,NAME_READ ,NAME_READ ,NAME_READ ,NAME_READ
      READ(ITI,*) AVERAGEX(I),NAME_READ ,NAME_READ ,NAME_READ
      READ(ITI,*) UHDX(I)   ,NAME_READ,NAME_READ,NAME_READ,NAME_READ      
      
      chana = 3
      
      
!      READ(ITI,*) 
!      READ(ITI,*) SEABEDx(i),WVHIGHTx(i),WDEPTHx(i),H1POSx(i),H2POSx(i)
!      READ(ITI,*) RWAVEx(i),ORDERx(i),PERIODx(i),GRAVx(i),RHOWx(i),RHOAx(i),WKFX(I),WFCX(I)
!      READ(ITI,*) NCURRENTX(I),POWERLAWX(I),VWIND0x(i),VTIDEX(I),VCURRENTPX(I),H0X(I),CBFX(I)
!      READ(ITI,*) VCURRENTLX(I,1:5)
      
C     ================================================================================      
C        MUTIPLY FACTOR ( SEE DNV OS-J101 JULY 2011 PAGE 56. AND 61400-3 IEC:2009 )
C     ================================================================================
      
!      write(*,*) 'Hello Chana'
 !     stop

C      VWIND0X(I)  = VWIND0X(I)*FACTORX(I)
 !     READ(ITI,*) APx(i),SPx(i),CSx(i)
!      READ(ITI,*) NWINDX(I),Z0X(I),UHX(I),RHIGHx(I),ALPHAX(I),UHAPIX(I),AVERAGEX(I),UHDX(I)

C     ================================================================================      
C        NUMBER OF RANDOM WAVE SIMULATION BY CHANA
C     ================================================================================
!     
!     READ(ITI,*)
!     READ(ITI,*) FSpectrumStartX(I),NAME_READ,NAME_READ,NAME_READ 
!     READ(ITI,*) FSpectrumEndX(I),NAME_READ,NAME_READ,NAME_READ
!     READ(ITI,*) NUMBEROFRANDOMWANVEXInterationX(I),NAME_READ,NAME_READ,NAME_READ,NAME_READ

C     ================================================   
C         MODIFY BY TOEY 30-08-2012 FOR WAVE ANGLE 
C     ================================================  
      
      READ(ITI,*)
      READ(ITI,*) WAVEANGLE ,NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) WINDANGLE ,NAME_READ,NAME_READ,NAME_READ 
      READ(ITI,*) 
      
      chana = 3
       
        CALL WAVE_VECTOR (WAVEANGLE,WAVE)
        VWAVEX1(I) = WAVE(1)
        VWAVEx2(I) = WAVE(2)
        VWAVEx3(I) = WAVE(3)
     
        CALL WIND_VECTOR (WINDANGLE,WIND)
        VWINDx1(I) = WIND(1)
        VWINDx2(I) = WIND(2)
        VWINDx3(I) = WIND(3)
        
        CALL WIND_VECTOR (BREAKING,BREAKING_VE)
        VBREAKINGx1(I) = BREAKING_VE(1)
        VBREAKINGx2(I) = BREAKING_VE(2)
        VBREAKINGx3(I) = BREAKING_VE(3)
        
        
      IF( KSPEC .NE. 1 ) THEN
          
      !==============================  
      IF (RWAVEx(i) == 1) THEN
      !==============================
          
      WRITE(ITO,200)
      WRITE(10,200)
      WRITE(ITO,201) I  
      WRITE(10,201) I  
      WRITE(ITO,102)
      WRITE(10,102)
102   FORMAT("WAVE THEORY ************         AIRY")   
      WRITE(ITO,209) WVHIGHTx(i)
      WRITE(10,209) WVHIGHTx(i)
      WRITE(ITO,203) WDEPTHx(i)
      WRITE(10,203) WDEPTHx(i)
      WRITE(ITO,204) PERIODx(i)
      WRITE(10,204) PERIODx(i)
      
      PI = 3.14159D0
      OMEGA = 2.0*PI/PERIODx(i)     
      
      CALL NEWTON_RAPHSON(WDEPTHx(i),GRAVx(i),OMEGA,RK)
      
      CALL Airy_Wave_Number(GRAVx(i),OMEGA,WDEPTHx(i),RK)
      
      RAMDAAIRY = 2.0*PI/RK
      WRITE(ITO,205) RAMDAAIRY
      WRITE(10,205) RAMDAAIRY
      WRITE(ITO,103) RK
      WRITE(10,103) RK
103   FORMAT("WAVE NUMBER ************ ",F12.4)   
      WRITE(ITO,207) RAMDAAIRY/PERIODx(i)  
      WRITE(10,207) RAMDAAIRY/PERIODx(i) 
      WRITE(ITO,208) SEABEDx(i)
      WRITE(10,208) SEABEDx(i)
      WRITE(ITO,104) (WDEPTHx(i)+0.5d0*WVHIGHTx(i))
      WRITE(10,104) (WDEPTHx(i)+0.5d0*WVHIGHTx(i))
104   FORMAT("CREST WATER DEPTH ****** ",F12.4)   
      WRITE(ITO,105) (WDEPTHx(i)-0.5d0*WVHIGHTx(i))
      WRITE(10,105) (WDEPTHx(i)-0.5d0*WVHIGHTx(i))
105   FORMAT("TROUGH WATER DEPTH ***** ",F12.4) 
      WRITE(ITO,200)
      WRITE(10,200) 
      
      WRITE(ITO,*)
      WRITE(10,*)
      
      !==============================
      ELSEIF (RWAVEx(i) == 2) THEN
      !==============================
          
      CALL STOKES_WAVELENGTH(abs(GRAVx(i)),WDEPTHx(i),PERIODx(i),WVHIGHTx(i),aLamda2,RAMDA)
      PI        = 3.14D0
      ACOEFFICIENT(i) = aLamda2 ! STOKE_COEFFICIENT
      DDLL =  WDEPTHx(i)/RAMDA
      RKx(i)    = 2*pi/RAMDA
      RAMDAx(i) = RAMDA
      
      WRITE(ITO,200)
      WRITE(10,200)
200   FORMAT("=========================================================")   
      WRITE(ITO,201)I
      WRITE(10,201)I
201   FORMAT("OFFSHORE PARAMETER ",I4)
      WRITE(ITO,202)
      WRITE(10,202)
202   FORMAT("WAVE THEORY ************   STOKES 5TH") 
      WRITE(ITO,209)WVHIGHTx(i)
      WRITE(10,209)WVHIGHTx(i)
209   FORMAT("WAVE HEIGHT ************ ",F12.4)   
      WRITE(ITO,203)WDEPTHx(i)
      WRITE(10,203)WDEPTHx(i)
203   FORMAT("WATER DEPTH ************ ",F12.4)       
      WRITE(ITO,204)PERIODx(i)
      WRITE(10,204)PERIODx(i)
204   FORMAT("WAVE PERIOD ************ ",F12.4)
      WRITE(ITO,205)RAMDAx(i)
      WRITE(10,205)RAMDAx(i)
205   FORMAT("WAVE LENGTH ************ ",F12.4)
      WRITE(ITO,206)ACOEFFICIENT(i)
      WRITE(10,206)ACOEFFICIENT(i)
206   FORMAT("WAVE COEFFICIENT ******* ",F12.5)  
      WRITE(ITO,210)DDLL
      WRITE(10,210)DDLL
210   FORMAT("D/L RATIO0 ************* ",F12.5)       
      WRITE(ITO,207)RAMDAx(i)/PERIODx(i)
      WRITE(10,207)RAMDAx(i)/PERIODx(i)
207   FORMAT("WAVE CELERITY ********** ",F12.5) 
      WRITE(ITO,208)SEABEDx(i)
      WRITE(10,208)SEABEDx(i)
208   FORMAT("MUDLINE ELEVATION ****** ",F12.5) 
      WRITE(ITO,200)    
      WRITE(10,200) 
      
      WRITE(ITO,*)
      WRITE(10,*)
      
      !==============================
      ELSEIF (RWAVEx(i) == 3) THEN
      !==============================
          
      WRITE(ITO,200)
      WRITE(10,200)
      WRITE(ITO,201) I  
      WRITE(10,201)  I  
      WRITE(ITO,302)
      WRITE(10,302)
302   FORMAT("WAVE THEORY ************     STREAM FUNCTIN ")   
      WRITE(ITO,209)WVHIGHTx(i)
      WRITE(10,209)WVHIGHTx(i)
      WRITE(ITO,203)WDEPTHx(i)
      WRITE(10,203)WDEPTHx(i)
      WRITE(ITO,204)PERIODx(i)   
      WRITE(10,204)PERIODx(i) 
      WRITE(ITO,303)ORDERx(i)
      WRITE(10,303)ORDERx(i)
303   FORMAT("STREAM ORDER ***********     ",F5.0)  
      ! WNU = WAVE CREST
      TIMEXSS = 0.0D0
      CALL STREAMWAVEFUNCTION (WVHIGHTx(i),PERIODx(i),WDEPTHx(i),TIMEXSS,WNU,abs(GRAVx(i)),ORDERx(i))
      
    	! FIND RAMDA AND WAVE NUMBER
    	CALL VELOC (0,0,U,V,DUDT,DVDT,TT,ARAMDA)
      WRITE(ITO,205)ARAMDA
      WRITE(10,205)ARAMDA
      WRITE(ITO,207)ARAMDA/PERIODx(i)
      WRITE(10,207)ARAMDA/PERIODx(i)
      WRITE(ITO,208)SEABEDx(i)
      WRITE(10,208)SEABEDx(i)
      WRITE(ITO,407) (WNU+WDEPTHx(i))
      WRITE(10,407) (WNU+WDEPTHx(i))
      WRITE(ITO,408) (WNU+WDEPTHx(i)-WVHIGHTx(i))
      WRITE(10,408) (WNU+WDEPTHx(i)-WVHIGHTx(i))
      WRITE(ITO,200)  
      WRITE(10,200)  
      
      WRITE(ITO,*)
      WRITE(10,*) 
          
          

      
      !==============================
      ELSEIF (RWAVEx(i) == 4) THEN
      !==============================
      
      WRITE(ITO,200)
      WRITE(10,200)
      WRITE(ITO,201) I  
      WRITE(10,201) I  
      WRITE(ITO,402)
      WRITE(10,402)
402   FORMAT("WAVE THEORY ************     CNOIDAL WAVE")   
      WRITE(ITO,209)WVHIGHTx(i)
      WRITE(10,209)WVHIGHTx(i)
      WRITE(ITO,203)WDEPTHx(i)
      WRITE(10,203)WDEPTHx(i)
      WRITE(ITO,204)PERIODx(i)
      WRITE(10,204)PERIODx(i)
      TIMEXCA = 0.0D0
      call Cnoidal_Wave_Crest  (WVHIGHTx(i),WDEPTHx(i),abs(GRAVx(i)),PERIODx(i),TIMEXCA,Wave_Crest)
    	CALL Cnoidal_wave_Length (WVHIGHTx(i),WDEPTHx(i),abs(GRAVx(i)),PERIODx(i),wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
      WRITE(ITO,205)ALAMDA
      WRITE(10,205)ALAMDA
      WRITE(ITO,207)ALAMDA/PERIODx(i)
      WRITE(10,207)ALAMDA/PERIODx(i)
      WRITE(ITO,405)WAVEANGLE
      WRITE(10,405)WAVEANGLE
405   FORMAT("WAVE ANGLE ************* ",f12.5)  
      WRITE(ITO,406) SEABEDx(i)
      WRITE(10,406) SEABEDx(i)
406   FORMAT("MUDLINE ELEVATION ****** ",f12.5) 
      WRITE(ITO,407) (Wave_Crest+WDEPTHx(i))
      WRITE(10,407) (Wave_Crest+WDEPTHx(i))
407   FORMAT("MAXIMUM WAVE CREST ***** ",f12.5) 
      WRITE(ITO,408) (Wave_Crest+WDEPTHx(i)-WVHIGHTx(i))
      WRITE(10,408) (Wave_Crest+WDEPTHx(i)-WVHIGHTx(i))
408   FORMAT("MINIMUM WAVE CREST ***** ",f12.5) 
      WRITE(ITO,409) Ammm 
      WRITE(10,409) Ammm 
409   FORMAT("ELLIPTIC FUNCTIN MODULOUS     ***** ",f12.5)   
      WRITE(ITO,410) AKKK 
      WRITE(10,410) AKKK 
410   FORMAT("FIRST KIND ELLIPTIC INTEGRAL  ***** ",f12.5)  
      WRITE(ITO,411) AEEE 
      WRITE(10,411) AEEE 
411   FORMAT("SECOND KIND ELLIPTIC INTEGRAL ***** ",f12.5)        
      
      CALL CNOIDALSURFACEELEVPLOT(WVHIGHTx(i),WDEPTHx(i),PERIODx(i),Wave_Crest,Ammm)
      
      
      WRITE(ITO,200)
      WRITE(10,200) 
      
      WRITE(ITO,*)   
      WRITE(10,*) 
      
      chana = 3
      
      !CALL CNOIEAL_CALL_MODULOUS(WVHIGHTx(i),WDEPTHx(i),abs(GRAVx(i)),AMODOUS)

      !CALL JACOBIAN_ELIPFUNCTION

      ENDIF
      
      ENDIF
      
      

      ENDDO
      
      
      CALL STOKES_WAVELENGTH_SURFACE_PLOT
           

C     ================================================   
C                 READ SPECTRUM PARAMETER
C     ================================================   
      CALL READSPECTRUM (NSPECTRUM)
      
      WRITE(ITO,*)'READ OFFSHORE LOAD PARAMETERS'
      WRITE(10,*)'READ OFFSHORE LOAD PARAMETERS'
      
      
      
      ENDIF
      
      
      RETURN
      END SUBROUTINE
C	=======================================================================
C     =======================================================================
C     ================== GENERATE WAVE AND WIND VECTOR ======================
C     =======================================================================
      SUBROUTINE WAVE_VECTOR (WAVEANGLE,WAVE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION ANGLE(100),WAVE(3)
      COMMON /MGRAV/ NGRAV  
      PI=3.141592684
      ! INPUT - WAVEANGLE THIS FROM TCL AND TK
      ! OUTPUT - VWAVEX1,VWAVEx2,VWAVEx3 > WAVE VECTOR
      IF (WAVEANGLE.LE.360D0)THEN
      GOTO 500
      ELSEIF (WAVEANGLE.GT.360D0)THEN
        ICALCULATE=WAVEANGLE/360D0
        DO I=1,ICALCULATE
        WAVEANGLE=WAVEANGLE-360
        IF (WAVEANGLE.LT.360D0)THEN
        GOTO 500
        ENDIF
        ENDDO
      ENDIF
      
500   IF (WAVEANGLE.EQ.0.0D0.OR.WAVEANGLE.EQ.360D0)THEN
      XWAVE=1.0D0
      YWAVE=0.0D0
      ZWAVE=0.0D0
      ELSEIF (WAVEANGLE.GT.0D0.AND.WAVEANGLE.LT.90D0)THEN
      RAD=PI/180D0
      XWAVE=COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=SIN(WAVEANGLE*RAD)
      ELSEIF (WAVEANGLE.EQ.90D0)THEN
      XWAVE=0.0D0
      YWAVE=0.0D0
      ZWAVE=1.0D0
      ELSEIF (WAVEANGLE.GT.90D0.AND.WAVEANGLE.LT.180D0)THEN
      RAD=PI/180D0
      XWAVE=1D0*COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
      ELSEIF (WAVEANGLE.EQ.180)THEN
      XWAVE=-1.0D0
      YWAVE=0.0D0
      ZWAVE=0.0D0
      ELSEIF (WAVEANGLE.GT.180D0.AND.WAVEANGLE.LT.270D0)THEN
      RAD=PI/180D0
      XWAVE=1D0*COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
      ELSEIF (WAVEANGLE.EQ.270D0)THEN
      XWAVE=0.0D0
      YWAVE=0.0D0
      ZWAVE=-1.0D0
      ELSEIF (WAVEANGLE.GT.270D0.AND.WAVEANGLE.LT.360)THEN
      RAD=PI/180D0
      XWAVE=1D0*COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
      ENDIF
      
      SELECTCASE(NGRAV)
      CASE(1)
      WAVE(1)=YWAVE
      WAVE(2)=XWAVE
      WAVE(3)=ZWAVE    
      CASE(2)   
      WAVE(1)=XWAVE
      WAVE(2)=YWAVE
      WAVE(3)=ZWAVE
      CASE(3)
      WAVE(1)=XWAVE
      WAVE(3)=YWAVE
      WAVE(2)=ZWAVE
      ENDSELECT
      
      RETURN
      END SUBROUTINE
C     ------------------------------------------
      SUBROUTINE WIND_VECTOR (WINDANGLE,WIND)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION WIND(3)
       !INPUT - WINDANGLE THIS FROM TCL AND TK
       !OUTPUT - VWINDx1,VWINDx2,VWINDx3 > WIND VECTOR
      PI=3.141592684
      WAVEANGLE=WINDANGLE
      IF (WAVEANGLE.LE.360D0)THEN
      GOTO 500
      ELSEIF (WAVEANGLE.GT.360D0)THEN
        ICALCULATE=WAVEANGLE/360D0
        DO I=1,ICALCULATE
        WAVEANGLE=WAVEANGLE-360
        IF (WAVEANGLE.LT.360D0)THEN
        GOTO 500
        ENDIF
        ENDDO
      ENDIF
      
500   IF (WAVEANGLE.EQ.0.0D0.OR.WAVEANGLE.EQ.360D0)THEN
      XWAVE=1.0D0
      YWAVE=0.0D0
      ZWAVE=0.0D0
      ELSEIF (WAVEANGLE.GT.0D0.AND.WAVEANGLE.LT.90D0)THEN
      RAD=PI/180D0
      XWAVE=COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=SIN(WAVEANGLE*RAD)
      ELSEIF (WAVEANGLE.EQ.90D0)THEN
      XWAVE=0.0D0
      YWAVE=0.0D0
      ZWAVE=1.0D0
      ELSEIF (WAVEANGLE.GT.90D0.AND.WAVEANGLE.LT.180D0)THEN
      RAD=PI/180D0
      XWAVE=1D0*COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
      ELSEIF (WAVEANGLE.EQ.180)THEN
      XWAVE=-1.0D0
      YWAVE=0.0D0
      ZWAVE=0.0D0
      ELSEIF (WAVEANGLE.GT.180D0.AND.WAVEANGLE.LT.270D0)THEN
      RAD=PI/180D0
      XWAVE=1D0*COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
      ELSEIF (WAVEANGLE.EQ.270D0)THEN
      XWAVE=0.0D0
      YWAVE=0.0D0
      ZWAVE=-1.0D0
      ELSEIF (WAVEANGLE.GT.270D0.AND.WAVEANGLE.LT.360)THEN
      RAD=PI/180D0
      XWAVE=1D0*COS(WAVEANGLE*RAD)
      YWAVE=0.0D0
      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
      ENDIF
      WIND(1)=XWAVE
      WIND(2)=YWAVE
      WIND(3)=ZWAVE
      END
C	=======================================================================
C     =============== Modify stroke ramda for reducing time  ================
C     =======================================================================
      SUBROUTINE STOKES_WAVELENGTH_Time_Modify(ACOEFFICIENTSTOKE,RAMDA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
	COMMON /OFFSHOREy_Data_Selection/ SEABEDy(1000),WVHIGHTy(1000),WDEPTHy(1000),THIGHTy(1000),H1POSy(1000),H2POSy(1000),
	1                  RWAVEy(1000),ORDERx(1000),PERIODy(1000),GRAVy(1000),RHOWy(1000),RHOAy(1000),
	1                  VWIND0y(1000),VTIDEy(1000),H0y(1000),
	1                  APy(1000),SPy(1000),CSy(1000),HMY(1000),HWY(1000),HCY(1000),VBy(1000,3),
	1                  RHIGHy(1000),UHy(1000),ALPHAy(1000),Z0y(1000),
	1                  WVZETAy(1000),WVTIMEy(1000),PEAKWLEVy(1000),VGy(1000,3),V1y(1000,3),V2y(1000,3),VWy(1000,3),
	1                  NCURRENTY(1000),POWERLAWY(1000),VCURRENTLY(1000,5),VCURRENTAPIY(1000),UHAPIX(1000),NWINDX(1000),
     1                  VCURRENTPY(1000),AVERAGEY(1000),UHDY(1000),
     1                  WKFY(1000),CBFY(1000),WFCY(1000),
     1                  FSpectrumStartY(1000), FSpectrumEndY(1000), NUMBEROFRANDOMWANVEXInterationY(1000),
     1                  TIMEy(1000)
      
	COMMON /stroke_k_and_Ramda/RKx(1000),RAMDAx(1000),ACOEFFICIENT(1000)
	
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	
	Do i=1,NUM_OF_OFFSHORE_PARAMETER
		
	IF (i == offselect) THEN
	ACOEFFICIENTSTOKE   =   ACOEFFICIENT(i)
	RAMDA     =   RAMDAx(i)
      
      CHANA = 3
	  
	endif
	enddo
	
      END
      
C	=======================================================================
C	===================  OFFSPARA_READ_SELECTION 2  =======================
C	=======================================================================
C     --- SELECT OFFSHORE PARAMETER CASE ---
C     INPUT --- COMMON/OFFSHOREx ( ALL) 
      
      SUBROUTINE OFFSPARA_READ_COLLECT_DATA ! (OFF_SELECT)  
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /MGRAV/ NGRAV
C	=======================================================================
      
	COMMON /OFFSHOREx/ SEABEDx(1000),WVHIGHTx(1000),WDEPTHx(1000),THIGHTx(1000),H1POSx(1000),H2POSx(1000),
	1                  RWAVEx(1000),ORDERx(1000),PERIODx(1000),GRAVx(1000),RHOWx(1000),RHOAx(1000),
	1                  VWIND0x(1000),VTIDEx(1000),H0x(1000),
	1                  APx(1000),SPx(1000),CSx(1000),HMX(1000),HWX(1000),HCX(1000),VBREAKINGX1(1000),VBREAKINGX2(1000),
     1                  VBREAKINGX3(1000),IBREAKING(1000),
	1                  RHIGHx(1000),UHx(1000),ALPHAx(1000),Z0x(1000),
     1                  VWAVEx1(1000),VWAVEx2(1000),VWAVEx3(1000),VWINDx1(1000),VWINDx2(1000),VWINDx3(1000),
     1                  NCURRENTX(1000),POWERLAWX(1000),VCURRENTLX(1000,5),VCURRENTAPIX(1000),UHAPIX(1000),NWINDX(1000),
     1                  FACTORX(1000),VCURRENTPX(1000),AVERAGEX(1000),UHDX(1000),
     1                  WKFX(1000),CBFX(1000),WFCX(1000),
     1                  Offshoreparameterx(1000),
     1                  FSpectrumStartX(1000), FSpectrumEndX(1000), NUMBEROFRANDOMWANVEXInterationX(1000),
     1                  TIMEX(1000)
     
     	COMMON /OFFSHORE/ SEABED,WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,
	1                  RWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
	1                  VWIND0,VTIDE,H0,
	1                  AP,SP,CS,VB(3),HM,HW,HC,
	1                  RHIGH,UH,ALPHA,Z0,
	1                  VG(3),V1(3),V2(3),VW(3),WVZETA,WVTIME,PEAKWLEV,
     1                  POWERLAW,VCURRENTL(5),VCURRENTAPI,UHAPI,NCURRENT,NWIND,
     1                  FACTOR,VCURRENTP,AVERAGE,UHD,
     1                  WKF,CBF,WFC,							   
     1                  FSpectrumStart, FSpectrumEnd, NUMBEROFRANDOMWANVEXInteration,
     1                  TIME
	
	COMMON /OFFSHOREy_Data_Selection/ SEABEDy(1000),WVHIGHTy(1000),WDEPTHy(1000),THIGHTy(1000),H1POSy(1000),H2POSy(1000),
	1                  RWAVEy(1000),ORDERy(1000),PERIODy(1000),GRAVy(1000),RHOWy(1000),RHOAy(1000),
	1                  VWIND0y(1000),VTIDEy(1000),H0y(1000),
	1                  APy(1000),SPy(1000),CSy(1000),HMY(1000),HWY(1000),HCY(1000),VBy(1000,3),
	1                  RHIGHy(1000),UHy(1000),ALPHAy(1000),Z0y(1000),
	1                  WVZETAy(1000),WVTIMEy(1000),PEAKWLEVy(1000),VGy(1000,3),V1y(1000,3),V2y(1000,3),VWy(1000,3),
	1                  NCURRENTY(1000),POWERLAWY(1000),VCURRENTLY(1000,5),VCURRENTAPIY(1000),UHAPIY(1000),NWINDY(1000),
     1                  FACTORY(1000),VCURRENTPY(1000),AVERAGEY(1000),UHDY(1000),
     1                  WKFY(1000),CBFY(1000),WFCY(1000),
     1                  FSpectrumStartY(1000), FSpectrumEndY(1000), NUMBEROFRANDOMWANVEXInterationY(1000),
     1                  TIMEy(1000)
      
	COMMON /XVWAVE/ VWAVE(3),VWIND(3),VBREAKING(3)
	
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
C     =======================================================================
	
	DIMENSION VH1(3),VH2(3),VGV(3)!,VG(3),V1(3),V2(3),VW(3)
	
	DO i=1,NUM_OF_OFFSHORE_PARAMETER
	
      Offshoreparameter     =   Offshoreparameterx(i)
      TIME                  =   TIMEX(i)  
      SEABED                =   SEABEDx(i)
      WVHIGHT               =   WVHIGHTx(i)
      WDEPTH                =   WDEPTHx(i)
      THIGHT                =   THIGHTx(i)
      H1POS                 =   H1POSx(i)
      H2POS                 =   H2POSx(i)
      RWAVE                 =   RWAVEx(i)
      ORDER                 =   ORDERx(i)
      PERIOD                =   PERIODx(i)
      GRAV                  =   GRAVx(i)
      RHOW                  =   RHOWx(i)
      RHOA                  =   RHOAx(i)
      WKF                   =   WKFX(I)
      WFC                   =   WFCX(I)
      ! ---------- Generate Time domail form RANDOM ----------------
      
      FSpectrumStart        = FSpectrumStartX(I) 
      FSpectrumEnd          = FSpectrumEndX(I) 
      NUMBEROFRANDOMWANVEXInteration  = NUMBEROFRANDOMWANVEXInterationX(I) 
      
      
      ! --------------- CURRENT ----------------
      VWIND0                =   VWIND0x(i)
      NCURRENT              =   NCURRENTX(I)
      POWERLAW              =   POWERLAWX(I)
      VCURRENTL(1:5)        =   VCURRENTLX(I,1:5)   
      VCURRENTAPI           =   VCURRENTAPIX(I)
      VCURRENTP             =   VCURRENTPX(I)
      VTIDE                 =   VTIDEX(i)
      H0                    =   H0x(i)
      FACTOR                =   FACTORX(I)
      CBF                   =   CBFX(I)
      ! ----------------------------------------
      AP                    =   APx(i)
      SP                    =   SPx(i)
      CS                    =   CSx(i)
      RHIGH                 =   RHIGHx(i)
      HM                    =   HMX(I)
      HW                    =   HWX(I)
      HC                    =   HCX(I)
      VBREAKING(1)         =   VBREAKINGX1(I)
      VBREAKING(2)         =   VBREAKINGX2(I)
      VBREAKING(3)         =   VBREAKINGX3(I)
      ! ----------------- WIND ------------------
      AVERAGE               =   AVERAGEX(I)
      UHD                   =   UHDX(I)
      NWIND                 =   NWINDX(I)
      UH                    =   UHx(i)
      UHAPI                 =   UHAPIX(I)
      ALPHA                 =   ALPHAx(i)
      Z0                    =   Z0x(i)
      ! -----------------------------------------
      
     
      VWAVE(1)              =   VWAVEx1(i)
      VWAVE(2)              =   VWAVEx2(i)
      VWAVE(3)              =   VWAVEx3(i)
      
      VWIND(1)              =   VWINDx1(i)
      VWIND(2)              =   VWINDx2(i)
      VWIND(3)              =   VWINDx3(i)
      
      WVZETA = 0.0D0
      WVTIME = 0.0D0
      GRAV = ABS(GRAV)
      
	VWAVE(NGRAV)     = 0.0D0 
	VWIND(NGRAV)     = 0.0D0 
      VBREAKING(NGRAV) = 0.0D0 

	CALL SCALEN(VWAVE,VWAVE,DUM,3)
	CALL SCALEN(VWIND,VWIND,DUM,3)
      CALL SCALEN(VBREAKING,VBREAKING,DUM,3)
	
	VGV(1:3) = 0.0D0
	VGV(NGRAV) = 1.0D0
	VH1(1:3) = VWAVE(1:3)
      CALL VECPRD (VGV,VH1,VH2)
	CALL SCALEN(VH2,VH2,DUM,3)

      
	CALL INTFILL ('UNIT',IUL,1,2,0)
	FAL  = UNITFAC('LENG',IUL,0)
	Z0   = Z0*FAL  !Z0 in (m) unit
C     -------------------------
      VG(1:3) = VGV(1:3)
      V1(1:3) = VH1(1:3)
      V2(1:3) = VH2(1:3)
      VW(1:3) = VWIND(1:3)
      VB(1:3) = VBREAKING(1:3)
      
      SEABEDY(i)        =   SEABED
      TIMEy(i)        =   TIME
      WVHIGHTy(i)       =   WVHIGHT
      WDEPTHy(i)        =   WDEPTH
      ORDERY(I)         =   ORDER  
      THIGHTy(i)        =   THIGHT
      H1POSy(i)         =   H1POS
      H2POSy(i)         =   H2POS
      RWAVEy(i)         =   RWAVE
      ORDERY(I)         =   ORDER
      PERIODy(i)        =   PERIOD
      GRAVy(i)          =   GRAV
      RHOWy(i)          =   RHOW
      RHOAy(i)          =   RHOA
      WKFY(I)           =   WKF
      WFCY(I)           =   WFC
      ! ---------- Generate Time domail form RANDOM ----------------
      
      FSpectrumStartY(I)        = FSpectrumStart 
      FSpectrumEndY(I)          = FSpectrumEnd
      NUMBEROFRANDOMWANVEXInterationY(I)  = NUMBEROFRANDOMWANVEXInteration
      
      ! ------------ CURRENT ------------
      NCURRENTY(I)      =   NCURRENT
      POWERLAWY(I)      =   POWERLAW
      VCURRENTLY(I,1:5) =   VCURRENTL(1:5) 
      VCURRENTAPIY(I)   =   VCURRENTAPI
      VCURRENTPY(I)     =   VCURRENTP
      VWIND0y(i)        =   VWIND0
      VTIDEy(i)         =   VTIDE
      H0y(i)            =   H0
      FACTORY(I)        =   FACTOR
      CBFY(I)           =   CBF
      ! ---------------------------------
      APy(i)            =   AP
      SPy(i)            =   SP
      CSy(i)            =   CS
      HMY(I)            =   HM
      HWY(I)            =   HW
      HCY(I)            =   HC
      ! ------------- WIND --------------
      AVERAGEY(I)       =   AVERAGE
      UHDY(I)           =   UHD
      RHIGHy(i)         =   RHIGH
      NWINDY(I)         =   NWIND
      UHy(i)            =   UH
      UHAPIY(I)         =   UHAPI
      Z0y(i)            =   Z0
      ! --------------------------------
      
      ALPHAy(i) = ALPHA
      WVZETAy(i) = WVZETA
      WVTIMEy(i) = WVTIME
      PEAKWLEVy(i) = PEAKWLEV
      
      VGy(i,1:3)    =   VG(1:3) 
      V1y(i,1:3)    =   V1(1:3)
      V2y(i,1:3)    =   V2(1:3)
      VWy(i,1:3)    =   VW(1:3)
      VBy(i,1:3)    =   VB(1:3)
      
      ENDDO 
      
      END
C	=======================================================================
C	=======================================================================
      SUBROUTINE OFFSPARA_PARA_SELECT_DATA ! (OFF_SELECT)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
	COMMON /OFFSHOREy_Data_Selection/ SEABEDy(1000),WVHIGHTy(1000),WDEPTHy(1000),THIGHTy(1000),H1POSy(1000),H2POSy(1000),
	1                  RWAVEy(1000),ORDERY(1000),PERIODy(1000),GRAVy(1000),RHOWy(1000),RHOAy(1000),
	1                  VWIND0y(1000),VTIDEy(1000),H0y(1000),
	1                  APy(1000),SPy(1000),CSy(1000),HMY(1000),HWY(1000),HCY(1000),VBy(1000,3),
	1                  RHIGHy(1000),UHy(1000),ALPHAy(1000),Z0y(1000),
	1                  WVZETAy(1000),WVTIMEy(1000),PEAKWLEVy(1000),VGy(1000,3),V1y(1000,3),V2y(1000,3),VWy(1000,3),
     1                  NCURRENTY(1000),POWERLAWY(1000),VCURRENTLY(1000,5),VCURRENTAPIY(1000),UHAPIY(1000),NWINDY(1000),
     1                  FACTORY(1000),VCURRENTPY(1000),AVERAGEY(1000),UHDY(1000),
     1                  WKFY(1000),CBFY(1000),WFCY(1000),    
     1                  FSpectrumStartY(1000), FSpectrumEndY(1000), NUMBEROFRANDOMWANVEXInterationY(1000),
     1                  TIMEy(1000)
      
     	COMMON /OFFSHORE/ SEABED,WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,
	1                  RWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
	1                  VWIND0,VTIDE,H0,
	1                  AP,SP,CS,VB(3),HM,HW,HC,
	1                  RHIGH,UH,ALPHA,Z0,
	1                  VG(3),V1(3),V2(3),VW(3),WVZETA,WVTIME,PEAKWLEV,
     1                  POWERLAW,VCURRENTL(5),VCURRENTAPI,UHAPI,NCURRENT,NWIND,
     1                  FACTOR,VCURRENTP,AVERAGE,UHD,
     1                  WKF,CBF,WFC,					   
     1                  FSpectrumStart, FSpectrumEnd, NUMBEROFRANDOMWANVEXInteration,
     1                  TIME
	      
	
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	
	
	i=offselect
	
	SEABED       =   SEABEDy(i)
      TIME         =   TIMEY(i)
      WVHIGHT      =   WVHIGHTy(i)
      WDEPTH       =   WDEPTHy(i)
      THIGHT       =   THIGHTy(i)
      H1POS        =   H1POSy(i)
      H2POS        =   H2POSy(i)
      RWAVE        =   RWAVEy(i)
      ORDER        =   ORDERY(i)
      PERIOD       =   PERIODy(i)
      GRAV         =   GRAVy(i)
      RHOW         =   RHOWy(i)
      RHOA         =   RHOAy(i)
      WKF          =   WKFY(I)
      WFC          =   WFCY(I)
      
      ! ---------- Generate Time domail form RANDOM ----------------
      
      FSpectrumStart        = FSpectrumStartY(I) 
      FSpectrumEnd          = FSpectrumEndY(I) 
      NUMBEROFRANDOMWANVEXInteration  = NUMBEROFRANDOMWANVEXInterationY(I)    
      
      ! ----------- CURRENT ----------
      VWIND0            =   VWIND0y(i)
      VTIDE             =   VTIDEy(i)
      NCURRENT          =   NCURRENTY(I)
      POWERLAW          =   POWERLAWY(I)
      VCURRENTL(1:5)    =   VCURRENTLY(I,1:5)   
      VCURRENTAPI       =   VCURRENTAPIY(I)
      VCURRENTP         =   VCURRENTPY(I)
      H0                =   H0y(i)
      FACTORY           =   FACTORY(I)
      CBF               =   CBFY(I)
      
      ! ------------------------------
      
      AP           =   APy(i)
      SP           =   SPy(i)
      CS           =   CSy(i)
      HM           =   HMY(I)
      HW           =   HWY(I)
      VB(1)        = VBy(i,1)
      VB(2)        = VBy(i,2)
      VB(3)        = VBy(i,3)
      
      ! ------------ WIND ------------
      AVERAGE        =   AVERAGEY(I)
      UHD            =   UHDY(I)
      RHIGH          =   RHIGHy(i)
      NWIND          =   NWINDY(I)
      UH             =   UHy(i)
      UHAPI          =   UHAPIY(I)
      ALPHA          =   ALPHAy(i)
      Z0             =   Z0y(i)
      ! ------------------------------
      
      
      WVZETA       =   WVZETAy(i)
      WVTIME       =   WVTIMEy(i)
      PEAKWLEV     =   PEAKWLEVy(i)
      
      VG(1)   = VGy(i,1)
      VG(2)   = VGy(i,2)
      VG(3)   = VGy(i,3)
      
      V1(1)   = V1y(i,1)
      V1(2)   = V1y(i,2)
      V1(3)   = V1y(i,3)
      
      V2(1)   = V2y(i,1)
      V2(2)   = V2y(i,2)
      V2(3)   = V2y(i,3)
            
      VW(1)   = VWy(i,1)
      VW(2)   = VWy(i,2)
      VW(3)   = VWy(i,3)

      END
C	=======================================================================
C	=========================== OFFSPARA_CALL =============================
C	=======================================================================
      SUBROUTINE OFFSPARA_CALL (WVH,WDE,THI,H1P,H2P,IWV,OR,PED,GRA,RHO,RHA,
     1                  WVZ,VTI,VWI,H0V,APV,SPV,CSV,VBD,HMV,HWV,HCV,RHI,UHV,ALP,Z0V,WVT,
     1                  VGV,VH1,VH2,VWD,PEAK,SEAC,NCURRENTV,POWERLAWC,VCURRENTLV,
     1                  VCURRENTAPIV,UHVAPI,NWINDO,AVE,UD,VCURRENTPV,FACTORV,
     1                  WKFV,CBFV,WFCV,TIME_OUT)
      
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /OFFSHORE/ SEABED,WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,
	1                  RWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
	1                  VWIND0,VTIDE,H0,
	1                  AP,SP,CS,VB(3),HM,HW,HC,
	1                  RHIGH,UH,ALPHA,Z0,
	1                  VG(3),V1(3),V2(3),VW(3),WVZETA,WVTIME,PEAKWLEV,
	1                  POWERLAW,VCURRENTL(5),VCURRENTAPI,UHAPI,NCURRENT,NWIND,
     1                  FACTOR,VCURRENTP,AVERAGE,UHD,
     1                  WKF,CBF,WFC,					   
     1                  FSpectrumStart, FSpectrumEnd, NUMBEROFRANDOMWANVEXInteration,
     1                  TIME
	
      COMMON/INDEXSCATDIAGRAM/ INDEXSCAT
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER 
      
      DIMENSION VH1(3),VH2(3),VGV(3),VWD(3),VBD(3),VCURRENTLV(5)
      
      DIMENSION SPECTRALNATURALV(10),SPECTRALDAMPINGV(10)
      
      CALL OFFSPARA_PARA_SELECT_DATA
      
      
      
      VGV(1:3) = VG(1:3)
      VH1(1:3) = V1(1:3)
      VH2(1:3) = V2(1:3)  
      VWD(1:3) = VW(1:3) 
      VBD(1:3) = VB(1:3)
      
      SEAC = 0.
      TIME_OUT  = 0.0D0
      WVH  = 0.
      WDE  = 0.
      THI  = 0.
      H1P  = 0.
      H2P  = 0.
      IWV  = 0.
      OR   = 0.
      PED  = 0.
      GRA  = 0.
      RHO  = 0.
      RHA  = 0.
      WKFV = 0.
      WFCV = 0.
      
      SEAC              =  SEABED
      TIME_OUT          =  TIME      
      WVH               =  WVHIGHT
      WDE               =  WDEPTH
      THI               =  THIGHT
      H1P               =  H1POS
      H2P               =  H2POS
      IWV               =  INT(RWAVE)
      OR                =  ORDER
      PED               =  PERIOD
      GRA               =  GRAV
      RHO               =  RHOW
      RHA               =  RHOA
      WKFV              =  WKF
      WFCV              =  WFC
      ! ------------ CURRENT ------------
      NCURRENTV         =  NCURRENT
      POWERLAWC         =  POWERLAW
      VCURRENTLV(1:5)   =  VCURRENTL(1:5) 
      VCURRENTAPIV      =  VCURRENTAPI
      VCURRENTPV        =  VCURRENTP
      VWI               =  VWIND0
      VTI               =  VTIDE
      H0V               =  H0
      FACTORV           =  FACTOR
      CBFV              =  CBF
      ! ---------------------------------
      APV               =  AP
      SPV               =  SP
      CSV               =  CS
      HMV               =  HM
      HWV               =  HW
      HCV               =  HC
      ! ------------- WIND --------------
      AVE               =  AVERAGE
      UD                =  UHD
      RHI               =  RHIGH
      NWINDO            =  NWIND
      UHV               =  UH
      UHVAPI            =  UHAPI
      ALP               =  ALPHA
      Z0V               =  Z0
      ! ---------------------------------
      WVZ               =  WVZETA
      WVT               =  WVTIME
      PEAK              =  PEAKWLEV
      
      SEAC              =  SEABED
      WVH               =  WVHIGHT
      WDE               =  WDEPTH
      THI               =  THIGHT
      H1P               =  H1POS
      H2P               =  H2POS
      IWV               =  INT(RWAVE)
      OR                =  ORDER
      PED               =  PERIOD
      GRA               =  GRAV
      RHO               =  RHOW
      RHA               =  RHOA
      WKFV              =  WKF
      WFCV              =  WFC
          
    
      
	RETURN

      END