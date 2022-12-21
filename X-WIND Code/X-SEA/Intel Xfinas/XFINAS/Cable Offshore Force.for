C	=======================================================================
      SUBROUTINE CABLEOFFSHORE (LM,XYZ,IGSET,PROPG,MTSET,PROPM,IGIDM,
	1				         MLOAD,LRPIN,IHSET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	=======================================================================
      
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /MGRAV/ NGRAV
	
	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
	
	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
	
      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN

	COMMON /LCSS/ ILCN,ILCC
	
	COMMON /WARNING/ WARNING,RAMDA(100),RK
	
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION,NOFFSHORE
	
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	
	COMMON / EIGVPED / EIGVFREQ(1000),EIGVPERI(1000)
      COMMON / STOREMODE / RVECT(100000,100),NMOD
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      !===============================================================================
      !  FATIGUE 
      !===============================================================================
      COMMON/WAVE_SCATTER/ NUMSCAT(1000),HSCAT(1000),PERIODSCAT(1000),AMENITUDESCAT(1000)
      COMMON/FATIGUE_OFFSHORE/ NWAVESCAT, NWINDSCAT
      COMMON /WAVEFRAMEELEMENTVECT/ ELEMENTFXYZ(6), NXWAVEFRAME
      COMMON/chana_OFFSHORE/ Irandom 
C	==================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH
C	=======================================================================
	COMMON A(7000000),IA(6000000)

	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION XYZ(NCO*NNM,NELE),LM(NEF,NELE),R(NEQ)
	DIMENSION FIXF(7,2),IGIDM(NELE)
	DIMENSION FRMFOC(7,2),LRPIN(14,1),IHSET(1),IPIN(14)
	DIMENSION VR(3),WW(6),MLOAD(20)

C     FOR OFFSHORE LOAD
      DIMENSION VWIND(3),VH1(3),VH2(3),VGV(3),LWCASE(10)
      DIMENSION FWIND(3)
C     FOR OFFSHORE LOAD -- PRAMIN 
      DIMENSION PXYZ(2,3),FLIST(5),VREW(3)
      DIMENSION Z(100),VCURRENTL(5),DIAM3(12),SNA(10),SDG(10)
      
      DIMENSION RMA(NEQ),VFOMA(3),FRMFOCM(7,2),FIXFM(7,2)
      DIMENSION SVELO(NEQ),SACCE(NEQ),FORCE(NEQ)
      DIMENSION VELOINTER(3),ACCEINTER(3)
      DIMENSION VELO(3),ACCE(3)
      DIMENSION IEQ1(3),IEQ2(3)
      DIMENSION RL(NNF,NNM)
      
	ALLOCATABLE GPL(:),GPW(:),BPG(:),BWG(:)
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:)
      
      NOFFL = MLOAD(1)
      PI    = 3.141592654
      
C	======================================================================

C	======================================================================
C	OFFSHORE LOAD
C	======================================================================
	IF(NOFFL.EQ.0D0) GOTO 117

	ALLOCATE(RVAL(25,NOFFL),IVAL(20,NOFFL))

C     CLEAR MATRIX
      RVAL(1:25,1:NOFFL) = 0
      IVAL(1:20,1:NOFFL) = 0
      
      GROWTH =0.0D0
C     ----------------------------------------

C     ========================= READ INPUT DATA AND STORE ON TEMPORARY VARIABLE =========================
      READ (ITI,*) 
	DO IOFFL = 1,NOFFL
!        READ(ITI,*) MLE,NFUNCTION,ROUGH,CA,CD,CM,CSA,NDI,DIAM1,OFFSELECT_X,NGROWTH,GROWTH,LWCASE(1),LWCASE(2),
!     1              LWCASE(5),LWCASE(6),LWCASE(7),ILCN,ILCC !READ ONLY FIRST TIME STEP
!        
       READ(ITI,*) MLE,NFUNCTION,ROUGH,CA,CD,CM,CSA,NWAVESEGOPT,NWAVESEG,NSURFWSEGOPT,NSURFWSEGMT,NDI,DIAM1,NGROWTH,
     1     GROWTH,NIRRE,OFFSELECT_X,LWCASE(1),LWCASE(2),LWCASE(5),LWCASE(6),LWCASE(7),ILCN,ILCC !READ ONLY FIRST TIME STEP
        
        IVAL(1 ,IOFFL)  = MLE
        IVAL(2 ,IOFFL)  = LWCASE(1)
        IVAL(3 ,IOFFL)  = LWCASE(2)
        IVAL(4 ,IOFFL)  = 0.0D0     
        IVAL(5 ,IOFFL)  = 0.0D0    
        IVAL(6 ,IOFFL)  = LWCASE(5)
        IVAL(7 ,IOFFL)  = LWCASE(6)
        IVAL(8 ,IOFFL)  = LWCASE(7)
        IVAL(9 ,IOFFL)  = ILCN
        IVAL(10,IOFFL)  = ILCC
        IVAL(11,IOFFL)  = NIRRE
             
        RVAL(1 ,IOFFL) = CD
        RVAL(2 ,IOFFL) = CM
        ! ADD 06-05-2015
        RVAL(4,IOFFL)  = CA 
              
        
C       ---- CALCULATE MARINE GROWTH ----	         
        IF (NDI.EQ.2.0)THEN ! USER DEFINED
        DIAM2          =  DIAM1 ! DIAM1+(GROWTH*2)
        ENDIF
C        ---------------------------------
        
        RVAL(3 ,IOFFL) = DIAM2
        RVAL(19,IOFFL) = NDI
        RVAL(20,IOFFL) = GROWTH
        RVAL(21,IOFFL) = OFFSELECT_X
        RVAL(22,IOFFL) = NFUNCTION
        RVAL(23,IOFFL) = ROUGH
        RVAL(24,IOFFL) = DIAM1 
        RVAL(25,IOFFL) = CSA 
        

        
    ! -------------------------------  FOR ERROR MESSAGE MARINE GROWTH ------------------------------------ 
!      CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
!     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
!     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
!     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,3,STOPFUNCTION,
!     1                        NGROWTH,GROWTH)
!      IF (STOPFUNCTION.EQ.1) GOTO 101
	ENDDO

C     =====================================================================================================

   ! -------------------------------  FOR ERROR MESSAGE LOAD CASE ------------------------------- 
!101    CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
!     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
!     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
!     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,1,STOPFUNCTION,
!     1                        NGROWTH,GROWTH)
        
	IF(KOFFL.EQ.0) THEN
	  DEALLOCATE(RVAL,IVAL)
	  GOTO 117
      ENDIF
	! -----------------------
      
      !IF (ISOLOP.EQ.5.OR.ISOLOP.EQ.6) THEN
      !CALL STOREOFFSHOREELEMENT (NOFFL,IVAL,RVAL,LRPIN,"STOR")   
      !GOTO 117
      !ENDIF
  

	ALLOCATE(GPL(12),GPW(12),BPG(12),BWG(12),RRV(NEQ),RRC(NEQ))

	DO IGR = 1,12
	IF(IGR.EQ.1  ) GPL(IGR) = -1.0D0
	IF(IGR.EQ.12 ) GPL(IGR) =  1.0D0 
	IF(IGR.NE.1.AND.IGR.NE.12) GPL(IGR) =  GAUSP(IGR-1,12-2)
	IF(IGR.EQ.1  ) GPW(IGR) =  0.0D0
	IF(IGR.EQ.12 ) GPW(IGR) =  0.0D0 
	IF(IGR.NE.1.AND.IGR.NE.12) GPW(IGR) =  GAUSW(IGR-1,12-2)
	ENDDO

      IF(KOFFL.NE.0) THEN 
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') 
      ENDIF

      DO 7499 ILCAS = 1,LCS 
          
          
      Irandom = 1.0d0  
	      
      DO 7500 ITIME = 1,NTIME

      RRV(1:NEQ) = 0.0D0 
      RRC(1:NEQ) = 0.0D0 
      
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') 
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')
      CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','READ')


	DO 127 IOFFL = 1,NOFFL

        WFX         = 0.
        WFY         = 0.
        WFZ         = 0.
        MLE         = IVAL(1 ,IOFFL)
        LWCASE(1)   = IVAL(2 ,IOFFL)
        LWCASE(2)   = IVAL(3 ,IOFFL)
        LWCASE(3)   = IVAL(4 ,IOFFL)
        LWCASE(4)   = IVAL(5 ,IOFFL)
        LWCASE(5)   = IVAL(6 ,IOFFL)
        LWCASE(6)   = IVAL(7 ,IOFFL)
        LWCASE(7)   = IVAL(8 ,IOFFL)
        !LWCASE(8)   = IVAL(11 ,IOFFL)
        ILCN        = IVAL(9 ,IOFFL)
        ILCC        = IVAL(10,IOFFL)
        IORRE       = IVAL(11,IOFFL)
        GROWTH      = RVAL(20,IOFFL)
        NDI         = RVAL(19,IOFFL)
        COA         = RVAL(4,IOFFL) 
        
        NFUNCTION   = RVAL(22,IOFFL)
        ROUGH       = RVAL(23,IOFFL)
        
        ISET        = IGSET(MLE)
        DPG         = PROPG(2,ISET)
        DIAM        = DPG + (GROWTH*2D0) ! ADD MARINE GROWTH
        DIAM1       = RVAL(24,IOFFL)     ! NORMAL DIAMETER
        
        CSA         = RVAL(25,IOFFL) 
        IF (NFUNCTION.EQ.4.0D0)THEN
        CD          = RVAL(1 ,IOFFL)
        CM          = RVAL(2 ,IOFFL)
        ELSEIF (NFUNCTION.EQ.1.0.OR.NFUNCTION.EQ.2.0.OR.NFUNCTION.EQ.3.0D0)THEN
        CD          = 0 
        CM          = 0
        ENDIF
        
        ! --------------------------------------------
        
        
C   ==================================================================   
C   ======== Chana Modifile offshore loadcase selecttion
C   ======== 5 / July / 2012 
C   ==================================================================   
      offselect = RVAL(21,IOFFL) 
      
      IF (IORRE.EQ.1)THEN ! REGULAR WAVE 12-01-2021
      CALL OFFSPARA_READ_COLLECT_DATA
        
      CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,VB,HM,HW,HC,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                  VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,
     1                  WKF,CBF,WFC,ATIME)
      WARNING=GRAV ! WARINING : SENT GRAVITY TO XFINAS.FOR. USED FOR OFFSHORE WARNING MESSAGE
      ! LINEAR STATIC TIME AT ANALYSIS
      IF (ISOLOP.EQ.1) TIME = ATIME
      
      
C   ================================================================== 
 
   ! -------------------------------  FOR ERROR MESSAGE LOAD PARAMETER ------------------------------- 
!      CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
!     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
!     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
!     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,2,STOPFUNCTION,
!     1                        NGROWTH,GROWTH)
8     FORMAT ('') ! BLANK

      ! --- STOP FUNCTION ---
!      IF (STOPFUNCTION.EQ.1D0)THEN
!      STOP
!      ENDIF
      ! --- END STOP FUNCTION ---
      
      IF (LWCASE(2).EQ.0.0D0)THEN
      VTIDE        = 0.0D0
      VWIND0       = 0.0D0
      FACTOR       = 0.0D0
      VCURRENTAPI  = 0.0D0
      POWERLAW     = 0.0D0
      VCURRENTP    = 0.0D0
      DO I=1,5
      VCURRENTL(I) = 0.0D0
      ENDDO
      ENDIF
      
      ELSEIF (IORRE.EQ.2)THEN ! RANDOOM WAVE 12-01-2021
      CALL SELECTSPECTRUM (OFFSELECT,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
      ENDIF
C   ========================================================================================   

      IF(ILCN.NE.ILCAS.AND.ILCC.NE.ILCAS) GOTO 127
        
      ! READING DIAMETER
      CALL CABLE_DIAMETER (IA(LMS+MLE),A(LGP),DIAM)
      
C	LOC = LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR BODY LOAD)
	
	FRMFOC = 0.0D0

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 207
	ENDIF
	ENDDO
	EXIT
207	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)

C     ===================================================	
C     PRAMIN WAVE LOAD MODIFICATION
C     ===================================================   
      
C	WAVE REFERNCE COORDINATE
      H1REF = H1POS
      H2REF = H2POS
      HGREF = SEABED	     

      PXYZ(1,1) = XYZ(1,MLE)
      PXYZ(1,2) = XYZ(2,MLE)
      PXYZ(1,3) = XYZ(3,MLE)
      PXYZ(2,1) = XYZ(4,MLE)
      PXYZ(2,2) = XYZ(5,MLE)
      PXYZ(2,3) = XYZ(6,MLE)
      ELEMENTFXYZ(1:6) = XYZ(1:6,MLE)

	XMID = 0.5*(PXYZ(1,1) + PXYZ(2,1))
      YMID = 0.5*(PXYZ(1,2) + PXYZ(2,2))
      ZMID = 0.5*(PXYZ(1,3) + PXYZ(2,3))

      HM1 = VH1(1)*XMID + VH1(2)*YMID + VH1(3)*ZMID 
      HMG = VGV(1)*XMID + VGV(2)*YMID + VGV(3)*ZMID 
      HM2 = VH2(1)*XMID + VH2(2)*YMID + VH2(3)*ZMID
C     ----------------------------------------------------

C     ----------------------------------------------------
C     POSITION OF COORDINATES REFER TO WAVE REFERNCE COORDINATE (IN WAVE DIRECTION SYSTEM)
      H1M =  HM1 - H1REF !H1 DISTANCE REFER TO WAVE REFERNCE COORDINATE
      HGM =  HMG - HGREF !GRAVITY DISTANCE REFER TO WAVE REFERNCE COORDINATE
      H2M =  HM2 - H2REF !H2 DISTANCE REFER TO WAVE REFERNCE COORDINATE 
C     ----------------------------------------------------

C     ----------------------------------------------------
C     STRUCTURAL ALIGNMENT REFERENCE VECTOR IN WAVE DIRECTION SYSTEM 
      VREW(1) = VH1(1)*VR(1) + VH1(2)*VR(2) + VH1(3)*VR(3)
      VREW(2) = VGV(1)*VR(1) + VGV(2)*VR(2) + VGV(3)*VR(3) 
      VREW(3) = VH2(1)*VR(1) + VH2(2)*VR(2) + VH2(3)*VR(3) 
C     ----------------------------------------------------
      IF (IORRE.EQ.1) THEN
C     FIND (OMEGA,WAVENUMBER(RK),RAMDA,A,RATIO)  
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
          CALL NEWTON_RAPHSON(WDEPTH,GRAV,OMEGA,RK)
          IF (RK.LE.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,30)
          WRITE (*,31)
30        FORMAT ('************** OFFSHORE WARNING MESSAGE ( AIRY WAVE THEORY FOR FRAME ) **************')
31        FORMAT ('- PLEASE CHECK WAVE PERIOD AND THEORY CONDITION ( WAVE NUMBER < 0 )')
          ENDIF
          RAMDA(ILCAS) = 2.0*PI/RK    
          AVAL = 1.0D0
          RATIO = 1.0D0   ! SET CONDITON     
      ELSEIF (IWAVE == 2) THEN
          ! --- STOKE FIFTH ORDER THEORY ---
          ! --- INPUT DATA ---
          ! GRAV,WDEPTH,PERIOD,WVHIGHT,RK,RAMDA,
          ! --- OUTPUT DATA ---
          ! ADUM    
             
          !!- CALL STOKES_WAVELENGTH(GRAV,WDEPTH,PERIOD,WVHIGHT,RK,RAMDA)
          CALL STOKES_WAVELENGTH_Time_Modify(AVAL,RAMDA(ILCAS)) ! RK = ASSUME EQUAL, AVAL
          RATIO = WDEPTH/RAMDA(ILCAS)
          RK=AVAL
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
          
          ! --- MODIFINE ON 22-08-2012 ----
          RK=(2*(RK+((RK**3)*F33)+((RK**5)*(F35+F55))))/WVHIGHT 
          IF (RK.LE.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,32)
          WRITE (*,33)
32        FORMAT ('************** OFFSHORE WARNING MESSAGE ( STOKE FIFTH ORDER THEORY FOR FRAME) **************')
33        FORMAT ('- PLEASE CHECK WAVE PERIOD AND THEORY CONDITION ( WAVE NUMBER < 0 )')
          ENDIF
          ! -----------------------          
          
          !	---------------------------------
          !	FREE-SURFACE WATER DEFLECTION (NU)
          CALL PARAMETER_A(RK,WVHIGHT,RATIO,AVAL)     
          OMEGA = SQRT(GRAV*RK*(1.0D0 + C1*(AVAL**2.0) + C2*AVAL**4.0)*TANH(RK*WDEPTH))
          !	---------------------------------
          F1 = AVAL
          F2 = (AVAL**2.0)*F22 + (AVAL**4.0)*F24
          F3 = (AVAL**3.0)*F33 + (AVAL**5.0)*F35
          F4 = (AVAL**4.0)*F44
          F5 = (AVAL**5.0)*F55
          !	---------------------------------   
          FLIST(1) = F1
          FLIST(2) = F2
          FLIST(3) = F3
          FLIST(4) = F4
          FLIST(5) = F5  
      ELSEIF (IWAVE == 3) THEN            
      ! FOR STEAMFUNCTION DON'T HAVE WAVE LENTGTH AND WAVE NUMBER
      ! NO NEED
      ENDIF   
            
      
           
C	---------------------------------------------------
C	WATER SURFACE LEVEL
C     MAXIMUM WATER SURFACE, SO APPLY ON THE STRUCTURE    
      ! IWAVE  = 1 >> AIRY WAVE THEORY
      ! IWAVE  = 2 >> STOKE'S FIFTH ORDER THEORY
      ! IWAVE  = 3 >> STREAM FUNCTION WAVE THEORY
      ! IWAVE  = 4 >> CNOIDAL WAVE THEORY
      ! IWAVE  = 5 >> SOLITATY WAVE THEORY  
      ! IWAVE  = 6 >> DIFFRACTION THEORY
      
        
           
      IF (IWAVE == 1  ) THEN 
      WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
      
      ELSEIF (IWAVE == 2) THEN
      WNU = (1.0D0/RK)*(F1*COS(1.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F2*COS(2.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F3*COS(3.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F4*COS(4.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F5*COS(5.0D0*(RK*H1M - OMEGA*TIME)))
     
      ELSEIF (IWAVE == 3 )THEN
      CALL STREAMWAVEFUNCTION (WVHIGHT,PERIOD,WDEPTH,TIME,
    	1                         WNU,GRAV,ORDER)
    	
    	! FIND RAMDA AND WAVE NUMBER
    	CALL VELOC (XC,YC,U,V,DUDT,DVDT,TT,ARAMDA)
    	RAMDA(ILCAS) = ARAMDA
    	RK=2D0*3.141592654D0/RAMDA(ILCAS)
    	
    	ELSEIF (IWAVE.EQ.4)THEN
    	call Cnoidal_Wave_Crest  (WVHIGHT,WDEPTH,GRAV,PERIOD,TIME,Wave_Crest)
    	call Cnoidal_wave_Length (WVHIGHT,WDEPTH,GRAV,PERIOD,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
    	RK=wavenumber
    	RAMDA(ILCAS)=ALAMDA
    	WNU = Wave_Crest
    	
    	ELSEIF (IWAVE.EQ.5)THEN
    	Call Solitary_Wave (WVHIGHT,WDEPTH,GRAV,TIME,0d0,z,u,v,au,av,WNU,q)
 
      
      ELSEIF (IWAVE.EQ.6)THEN
          
      WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
      
      ELSEIF ( IWAVE == 7 .or. IWAVE == 8 ) THEN 
          
      WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
      
      ! WAVE ELEVATION 
      ENDIF
      
      ELSEIF (IORRE.EQ.2)THEN
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
      ENDIF  
      
C	---------------------------------------------------
C	WATER SURFACE PROFILE MAX
      YL    = WDEPTH + WNU
      IF (WNU.GE.0.0D0) YLMIN = WDEPTH - WNU
      IF (WNU.LT.0.0D0) YLMIN = YL
      
C     WATER SURFACE PROFILE MIN
      IF (IWAVE.EQ.1 .or. IWAVE.EQ.7 .or. IWAVE.EQ.8 )THEN
      YLMIN = WDEPTH - WNU
      ELSEIF (IWAVE.EQ.2.OR.IWAVE.EQ.3.OR.IWAVE.EQ.4.OR.IWAVE.EQ.5)THEN
      YLMIN = (WDEPTH + WNU) - WVHIGHT
      ENDIF
      
      
      ! ----------------- GENERATE DRAG AND INTIA COEFFICIENT -----------------
      ! INPUT  >> NDI = 1 : AUTOMATIC NORMINAL DIAMETER
      !        >> NDI = 2 : USER DEFINED NORMINAL DIAMETER
      ! OUTPUT >> CD  = DRAG COEFFICIENR
      !        >> CM  = INTIA COEFFICIENT
      IF (NDI.EQ.2.0D0)THEN ! USER DEFINED NORMINAL DIAMETER
          CALL DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM1,PERIOD,WVHIGHT,RK
     1                        ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                        ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
      ENDIF ! (ENDIF : USER DEFINED NORMINAL DIAMETER ) 
      ! -----------------------------------------------------------------------  
     
C     ===================================================	
C     ===================================================	

C     BPG = GAUSS POINT IN EACH ELEMENT
	DO IGR = 1,12
	RI = GPL(IGR)  !GAUSP(IGR,NGR)
	RW = GPW(IGR)  !GAUSW(IGR,NGR)
	BPG(IGR) = 0.5*ELN*(1.0 + RI)
	BWG(IGR) = 0.5*ELN*RW
	ENDDO
	
C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ----------------------------------------------------------
      DO 307 IGR = 1,12

C     GAUSS LOCATION ALONG ELEMENT AXIS
	BXD = BPG(IGR)
	BXW = BWG(IGR)

      XR = XYZ(1,MLE) + VR(1)*BXD
      YR = XYZ(2,MLE) + VR(2)*BXD
      ZR = XYZ(3,MLE) + VR(3)*BXD
      SELECTCASE(NGRAV)
      CASE(1)
        RLEV = XR
      CASE(2)
        RLEV = YR
      CASE(3)
        RLEV = ZR
      ENDSELECT
      

C     GAUSS COORDINATE IN WAVE DIRECTION SYSTEM	
      H1R = VH1(1)*XR + VH1(2)*YR + VH1(3)*ZR
      HGR = VGV(1)*XR + VGV(2)*YR + VGV(3)*ZR
      H2R = VH2(1)*XR + VH2(2)*YR + VH2(3)*ZR

C     GAUSS COORDINATE RELATIVE TO WAVE REFERENCE POINT (IN WAVE DIRECTION SYSTEM)
      H1R =  H1R - H1REF
      HGM =  HGR - HGREF
      H2M =  H2R - H2REF
      
      ! ----------------- GENERATE DRAG AND INTIA COEFFICIENT -----------------
      ! INPUT  >> NDI = 1 : AUTOMATIC NORMINAL DIAMETER
      !        >> NDI = 2 : USER DEFINED NORMINAL DIAMETER
      ! OUTPUT >> CD  = DRAG COEFFICIENR
      !        >> CM  = INTIA COEFFICIENT
      ! ****** AUTOMATIC GENERATE BASE ON GAUSS POINT *****
      IF (NDI.EQ.1.0D0)THEN ! AUTOMATIC NORMINAL DIAMETER
      CALL DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM,PERIOD,WVHIGHT,RK
     1                          ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                          ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
      ENDIF   ! (ENDIF : AUTOMATIC NORMINAL DIAMETER ) 
      ! -----------------------------------------------------------------------  
      

      IF (LWCASE(2).EQ.1.0D0.AND.LWCASE(1).EQ.0.0)THEN
      YL = WDEPTH
      ENDIF 
      
      CALL MARINE_GROWTH_CAL (RVAL(21,IOFFL),RLEV,GROWTH,DENSITY,IORRE) 
      DIAM = DIAM + (GROWTH*2D0)
      
      KFLACAL = 0
      
      WW(1:6) = 0.0D0
      
      DO 407 LCASE = 1,8
      ICASE = LWCASE(LCASE)
      IF(ICASE.EQ.0) GOTO 407
     
      ! SELECT CASE WITH LOAD FOR ANALYSIS ( LWCASE )
      SELECTCASE(LCASE)
      CASE(1,2,3) !WAVE LOAD - CURRENT LOAD - WATER TIDE LOAD    
      
      IF(KFLACAL.EQ.1) GOTO 407  !IF CURRENT&TIDAL LOAD ALREADY CALCULATE TOGETHER WITH WAVE LOAD, THEN JUMP
      
      IF(LWCASE(2).EQ.0) VWIND0 = 0.0D0
      IF(LWCASE(3).EQ.0) VTIDE  = 0.0D0
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407
          IF(WLEV.GT.YL) GOTO 407
      
      IF ( IWAVE . NE . 6D0 ) THEN
      CALL WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGMWAVEINP,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,IORRE,IRWAVE)
      
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2
      KFLACAL = 1
      
      ELSE
          ! CHANA
      CALL DIFFRECTION_WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,H2M,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA)
      
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2
      KFLACAL = 1
      
      ELEVATION_X = H1R
      ELEVATION_Y = HGM
      
!        WRITE(791,791)  ELEVATION_Y ,',',WFX
791     FORMAT(F12.5,A,E12.5)
      

      ENDIF

      CASE(5,6) !WAVE BREAKING LOAD  PLUNGING&SURGING
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407
          IF(WLEV.GT.YL) GOTO 407
      CALL WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGMWAVEINP,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,IORRE,IRWAVE)
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2
      
      CASE(7) !WIND LOAD 
          WLEV = RLEV - SEABED - WDEPTH
          IF(WLEV.LT.0.0D0) GOTO 407 
          ! ------------------------------------------------------------------------------
          ! GENERATE WIND LOAD ( DNV,API,IEC64100-1 )
          ! NWINDO = 1 ; LOGARITHMIC PROFILE   (DNV) 
          ! NWINDO = 2 ; POWER LAW PROFILE     (DNV)
          ! NWINDO = 3 ; WIND PROFILE AND GUST (API)
          ! NWINDO = 4 ; IEC 614000-1
          ! NWINDO = 5 ; USER DEFINED
            IF (NWINDO.EQ.1)THEN ! DNV > LOGARITHMIC PROFILE
            UZT  = UH*(LOG(WLEV/Z0)/LOG(RHIGH/Z0))
            
             !  ---- MATERMATIC PROTECT ----
             IF (WLEV.EQ.0.0D0)THEN 
             UZT = 0.0D0
             ENDIF
             !  ---------------------------
             
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
            
            ELSEIF (NWINDO.EQ.2.OR.NWINDO.EQ.5)THEN ! DNV > POWER LAW PROFILE NWINDO >> 2, USER DEFIND NWINDO >> 5 
                 IF (NWINDO.EQ.5)THEN ! SAME EQUATION FROM USER DEFINED AND POWER LAW IN DNV
                 UH = UHD             ! UH  = MEAN WIND SPEED ON DNV 
                 ENDIF                ! UHD = MEAN WIND SPEED ON USER DEFINED
            CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA(ILCAS),GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     2                        CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
     
            ELSEIF (NWINDO.EQ.3)THEN ! API > WIND PROFILE AND GUST
            C   = 0.0573D0*SQRT(1D0+0.0457D0*UHAPI)
            UZ  = UHAPI*(1+(C*LOG(WLEV/RHIGH)))
            ZI  = 0.06D0*(1+0.0131D0*UHAPI)*((WLEV/RHIGH)**(-0.22))
            UZT = UZ*(1-0.41D0*ZI*LOG(AVERAGE/3600D0))
              
              !  ---- MATERMATIC PROTECT ----
              IF (WLEV.EQ.0.0D0)THEN 
              UZT = 0.0D0
              ENDIF
              !  ---------------------------
              
            WFF = 0.5D0*RHOA*DIAM*UZT*UZT
      
            ELSEIF (NWINDO.EQ.4)THEN ! IEC 614000-1
            UZT = UHD*((WLEV/RHIGH)**ALPHA)
            WFF = 0.5D0*RHOA*DIAM*UZT*UZT
            ENDIF
          ! --------------------------------------------------------------------------------
      WFX = VWIND(1)*WFF*CSA ; WFY = VWIND(2)*WFF*CSA ; WFZ = VWIND(3)*WFF*CSA ;
      
      ENDSELECT
      
      
      ! --- TOTAL WAVE FORCE * GAUSS ---
      ! *** THIS DATA WILL SENT TO CALCULATION DISPLACMENT ***
	WW(1) = WW(1) + WFX*BXW
	WW(2) = WW(2) + WFX*BXW
	WW(3) = WW(3) + WFY*BXW
	WW(4) = WW(4) + WFY*BXW
	WW(5) = WW(5) + WFZ*BXW
	WW(6) = WW(6) + WFZ*BXW
      
          WFX = 0.
          WFY = 0.
          WFZ = 0.
	
407   CONTINUE
C     -------------------------- 	


C	-----------------------------------------------------
	DO I = 1,3
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	DISA = BXD
	DISB = BXD
	IDIR = I
	IND = 2 !PROJECTED LOAD
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
      

	
	CALL FRAMFIX(RP,RPB,DISA,DISB,ECR,ECS,ECT,IDIR,0,VR,
	1			   ELN,FIXF,RANG,LOC,LOE,IMOM,IPIN)
	
	DO J = 1,7
	DO K = 1,2
	FRMFOC(J,K) = FRMFOC(J,K) + FIXF(J,K)
	ENDDO
	ENDDO
	
      ENDDO
      
      RL(1:3,1) = FRMFOC(1:3,1)
      RL(1:3,2) = FRMFOC(1:3,2)
      
      CALL FIXSOL_OFF (RL,KEG,MLE,NEF,NNM) ! FOR FIXED FORCE OFFSHORE
C	-----------------------------------------------------
307   CONTINUE
C	-----------------------------------------------------
C	-----------------------------------------------------

C     MUST BE ADDED LATER  FOR FRAME FIXEND FORCED
C	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

	IK = 0
	DO I = 1,2
	DO J = 1,3
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + FRMFOC(J,I) !VARY LOAD
      IF (IEQ.NE.0.AND.ILCC.GT.0) RRC(IEQ) = RRC(IEQ) + FRMFOC(J,I) !CONSTANT LOAD
	ENDDO
	ENDDO

127	CONTINUE

      IF(KOFFL.NE.0) THEN
        
          
        CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
        CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','WRIT') !STORE CONSTANT LOAD
        IREAD=1D0
        
!        DO I= 1,20 
!            
!        WRITE(790,790)  I ,',',IREAD,',',RRV(IREAD)
!790     FORMAT(I5,A,I5,A,E12.5)    
!        IREAD=IREAD + 6D0    
!        ENDDO
        
      ENDIF      

7500  CONTINUE
      
7499  CONTINUE
      

	DEALLOCATE(GPL,GPW,BPG,BWG,RRV,RRC,RVAL,IVAL)
      
	
117   CONTINUE
C     ------------------------------
      END
      
      
      
C	=======================================================================
      SUBROUTINE CABLE_OFFSHORE_BOUN (LM,XYZ,IGSET,PROPG,MTSET,PROPM,IGIDM,
	1				         MLOAD,LRPIN,IHSET,NBOFL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,                      
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,                  
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT  
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
      
      COMMON /MGRAV/ NGRAV
      
      COMMON A(9000000),IA(9000000)
      
      ALLOCATABLE RRV(:)
      ALLOCATABLE NMLE(:),AOFFM(:),ARAT(:),NILCN(:),NILCC(:)
      DIMENSION LRPIN(14,1),IPIN(14),WW(6),VR(3)
      DIMENSION FIXF(7,2),FRMFOC(7,2)
      DIMENSION XYZ(NCO*NNM,NELE)
      DIMENSION IGIDM(NELE),IHSET(1),LM(NEF,NELE)
      DIMENSION RL(NNF,NNM)
      
      
      ALLOCATE(RRV(NEQ))
      ALLOCATE(NMLE(NBOFL),AOFFM(NBOFL),ARAT(NBOFL),NILCN(NBOFL),NILCC(NBOFL))
C      CALL CLEARA (RRV,NEQ)
      
      IF(KOFFL.NE.0) THEN !CALL ALL OFFSHORE LOAD PARAMETERS
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      ENDIF
      
      READ(ITI,*) ! READ HEAD
      
      DO KK = 1,NBOFL
      READ(ITI,*) MLE,offselect,ARATIO,ILCN,ILCC    
      NMLE(KK)  = MLE
      AOFFM(KK) = offselect
      ARAT(KK)  = ARATIO
      NILCN(KK) = ILCN
      NILCC(KK) = ILCC
      ENDDO
      
      DO 130 ILCAS = 1,LCS !LOOP OVER LOAD CASE NUBER
      
      DO 130 ITIME = 1,NTIME !LOOP OVER OFFSHORE TIME STEP

      RRV(1:NEQ) = 0.0D0 !INITIALIZE VARY OFFSHORE LOAD
      
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')

      
      DO 132 II = 1,NBOFL ! BUOYANCY FORCE NUMBER
          
       MLE       = NMLE(II) 
       offselect = AOFFM(II)
       ARATIO    = ARAT(II) 
       ILCN      = NILCN(II)
       ILCC      = NILCC(II)
      
      IF (ILCAS.NE.ILCN) GOTO 132
     
      CALL OFFSPARA_READ_COLLECT_DATA
        
      CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,VB,HM,HW,HC,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                  VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,
     1                  WKF,CBF,WFC,TIME)
      
      ! READING DIAMETER
      CALL CABLE_DIAMETER (IA(LMS+MLE),A(LGP),DIAM)

        
	FRMFOC = 0.0D0

	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR SELFWEIGHT LOAD)

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 210
	ENDIF
	ENDDO
	EXIT
210	CONTINUE

      VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
      ALENGTH = SQRT(VR(1)**2D0 + VR(2)**2D0 + VR(3)**2D0)
	CALL SCALEN(VR,VR,ELN,3)
	
      !XR = XYZ(1,MLE) + VR(1)*BXD
      !YR = XYZ(2,MLE) + VR(2)*BXD
      !ZR = XYZ(3,MLE) + VR(3)*BXD
      
      XR = (XYZ(4,MLE)+XYZ(1,MLE))*0.5D0
	YR = (XYZ(5,MLE)+XYZ(2,MLE))*0.5D0
	ZR = (XYZ(6,MLE)+XYZ(3,MLE))*0.5D0
      
      SELECTCASE(NGRAV)
      CASE(1)
        RLEV = XR
      CASE(2)
        RLEV = YR
      CASE(3)
        RLEV = ZR
      ENDSELECT
      
      
      ! CHECKING POSITION CABLE
      IF (SEABED.LT.0.0D0)THEN     
      WATER_ELEVATION = WDEPTH + SEABED
      ELSEIF (SEABED.EQ.0.0D0)THEN 
      WATER_ELEVATION = WDEPTH
      ELSEIF (SEABED.GT.0.0D0)THEN 
      WATER_ELEVATION = WDEPTH - SEABED
      ENDIF
      
      PI    = 3.141592654
      AREA  = PI*(DIAM**2D0)/4D0
      
      ! PGV
      ! VOLUMN OF ELEMENT IS CONVERT TO NODAL LOAD. THE ELEMENT LENGTH IS APPILED IN DISTRIBUTE FUNCTION
      
      IF (RLEV.GT.WATER_ELEVATION) BUOYANCYF = 0.0D0
      IF (RLEV.LE.WATER_ELEVATION) THEN
      !DEPTH_ELE = WATER_ELEVATION - RLEV     
      BUOYANCYF = RHOW*GRAV*AREA    !DEPTH_ELE*AREA !*ALENGTH
      ENDIF
      

      WW = 0.
      SELECTCASE (NGRAV)
      CASE (1)
      WW(1) = BUOYANCYF
      CASE (2)
      WW(2) = BUOYANCYF  
      CASE (3)
      WW(3) = BUOYANCYF
      ENDSELECT


C	-----------------------------------------------------
	DO I = 1,3
	RP   = WW(I)
	RPB  = RP
	DISA = 0.0D0
	DISB = ELN
	IDIR = I
	IND  = 1
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
	CALL FRAMFIX(RP,RPB,DISA,DISB,ECR,ECS,ECT,IDIR,0,VR,
	1			   ELN,FIXF,RANG,LOC,LOE,IMOM,IPIN)
	DO J = 1,7
	DO K = 1,2
	FRMFOC(J,K) = FRMFOC(J,K) + FIXF(J,K)
	ENDDO
	ENDDO
	
      ENDDO
      RL(1:3,1) = FRMFOC(1:3,1)
      RL(1:3,2) = FRMFOC(1:3,2)
      CALL FIXSOL_OFF (RL,KEG,MLE,NEF,NNM) ! FOR FIXED FORCE OFFSHORE
C	-----------------------------------------------------

    	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

      IK = 0
	DO I = 1,2
	DO J = 1,3
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + FRMFOC(J,I)
	ENDDO
	ENDDO
      
132	CONTINUE
      
      IF(KOFFL.NE.0) THEN
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
      ENDIF
      
130   CONTINUE

      DEALLOCATE(RRV)
      DEALLOCATE (NMLE,AOFFM,ARAT,NILCN,NILCC)
120   CONTINUE
      
      END
C	=======================================================================
      SUBROUTINE CABLE_DIAMETER (MSET,PROPG,DIAM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      DIMENSION PROPG(NGP)
      PI    = 3.141592654
      DIAM = SQRT(4D0*PROPG(2)/PI)
      END