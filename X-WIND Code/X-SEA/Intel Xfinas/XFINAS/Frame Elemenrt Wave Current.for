C     ==================================================  
      SUBROUTINE FRAMFOC_OFFSHORE_WAVE_CURRENT (LM,XYZ,IGSET,PROPG,MTSET,PROPM,IGIDM,
	1				    MLOAD,LRPIN,IHSET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAMEI(4)
      CHARACTER*3 UNDWMEMBER 
      CHARACTER*3 STRIPWATER1 
      CHARACTER*3 STRIPWATER2
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
      
      COMMON /MWRITETIONWAVEFORCE/ NWELEMENT
                                                                                   
	COMMON /MGRAV/ NGRAV                                                                  
	                               
	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
	                                                                                                         
	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)           
	                                                   
      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN

	COMMON /LCSS/ ILCN,ILCC
      
      COMMON/ OFFSHORE_CASE/ LGEN,IOFFL_CASE,IORI_OFFSHORE
	
	COMMON /WARNING/ WARNING,RAMDA(100),RK
	
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION,NOFFSHORE
	
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	
	COMMON / EIGVPED / EIGVFREQ(1000),EIGVPERI(1000)
      COMMON / STOREMODE / RVECT(100000,100),NMOD
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      
      COMMON / WAVESPECTRUMPLOT / NWSPECTRUMPLOT
      
      COMMON / NFATIGUEFREQUENCYOPT / NFATIGUEFREQUENCY
      COMMON / FRAME_WAVE_DIRECTION / BARING(3),VERTICALDIR(3)
      
      COMMON / SURFACE_ELE / YL_STABILITY_STATIC(500),YL_STABILITY_DYNAMIC(20000)
      !===============================================================================
      !  FATIGUE 
      !===============================================================================
      COMMON/WAVE_SCATTER/ NUMSCAT(1000),HSCAT(1000),PERIODSCAT(1000),AMENITUDESCAT(1000)
      COMMON/INDEXSCATDIAGRAM/ INDEXSCAT
      COMMON/FATIGUE_OFFSHORE/ NWAVESCAT, NWINDSCAT
      COMMON/chana_OFFSHORE/ Irandom               
         
      !sorn added 
      COMMON /sornOffset/ NELO,Nmemoffset(1000),NmemOffsetSelect(1000) ,Nmemtrueelement(1000),GloLoSw(1000) !total number of off set node
      
C	==================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH                                                         
C	=======================================================================
	COMMON A(9000000),IA(9000000)
C	==================================================================
C	TRAPIZOIDAL RULE   
      !COMMON / WAVEREACFIXPARAMETER /  WAVEREACFIX(9999,7,2)
      
      COMMON /WAVEFRAMEELEMENTVECT/ ELEMENTFXYZ(6), NXWAVEFRAME
      
      COMMON / WRITEDATATPZANGLE / WRITEDATATPZ
      
      !COMMON/RANDOM_ELV/ EVL_RAN(20000),ITIME,IOFFL
       
C	==================================================================
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION XYZ(NCO*NNM,NELE),LM(NEF,NELE),R(NEQ)
	DIMENSION FIXF(7,2),IGIDM(NELE)
	DIMENSION FRMFOC(7,2),LRPIN(14,1),IHSET(1),IPIN(14)
	DIMENSION VR(3),WW(6),MLOAD(20)
      
      DIMENSION CHAFOC(7,2)

C     FOR OFFSHORE LOAD
      DIMENSION VWIND(3),VH1(3),VH2(3),VGV(3),LWCASE(10)
      DIMENSION FWIND(3)
C     FOR OFFSHORE LOAD -- PRAMIN 
      DIMENSION PXYZ(2,3),FLIST(5),VREW(3)
      DIMENSION Z(100),VCURRENTL(5),DIAM3(12),SNA(10),SDG(10)
      
      DIMENSION RMA(NEQ),VFOMA(3),FRMFOCM(7,2),FIXFM(7,2)
      
	ALLOCATABLE GPL(:),GPW(:),BPG(:),BWG(:)
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:),AMARINE(:)
      
      ALLOCATABLE NMLE(:),AOFFM(:),ARAT(:),NILCN(:),NILCC(:)
      
      
	PI = 3.141592654
      WAVEREACFIX = 0.0D0
      
      IF (KSPEC.EQ.0)THEN
      NWELEMENT = NELE
      NSW       = MLOAD(1)
      NUF       = MLOAD(2)
      NMARINE   = MLOAD(3)
      NOFFL     = MLOAD(4)
      NCURRT    = MLOAD(5)
      NBOFL     = MLOAD(8)
      ELSEIF (KSPEC.EQ.1)THEN
      NMARINE   = MLOAD(3)
      NRES      = MLOAD(4) 
      NBOFL     = MLOAD(8)
      ENDIF
      
C	======================================================================
C	WAVE LOAD
C	======================================================================
	IF(NOFFL.EQ.0D0) GOTO 117

	ALLOCATE(RVAL(25,NOFFL),IVAL(20,NOFFL),AMARINE(NELE))

C     CLEAR MATRIX
      RVAL(1:25,1:NOFFL) = 0
      IVAL(1:20,1:NOFFL) = 0
      
      AMARINE = 0.
      GROWTH  = 0.0D0
      ! READ THE MARINE GROWTH THICKNESS
      IF (NMARINE.NE.0)THEN
      REWIND(70)
        DO I = 1,NMARINE
        READ (70,*) MLE_MARINE,THICKNESS
        AMARINE(MLE_MARINE) = THICKNESS
        ENDDO
      ENDIF
C     ----------------------------------------

C     ========================= READ INPUT DATA AND STORE ON TEMPORARY VARIABLE =========================
      READ (ITI,*) 
      
	DO IOFFL   = 1,NOFFL
        READ(ITI,*) MLE,NFUNCTION,ROUGH,CA,CD,CM,CSA,NWAVESEGOPT,NWAVESEG,NSURFWSEGOPT,NSURFWSEGMT,NDI,DIAM1,NGROWTH,
     1              GROWTH,NIRRE,OFFSELECT_X,LWCASE(1),LWCASE(2),LWCASE(5),LWCASE(6),LWCASE(7),ILCN,ILCC !READ ONLY FIRST TIME STEP
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
        IVAL(12,IOFFL)  = NSURFWSEGOPT
        IVAL(13,IOFFL)  = NSURFWSEGMT
        IVAL(14,IOFFL)  = NWAVESEGOPT
        IVAL(15,IOFFL)  = NWAVESEG
        IVAL(16,IOFFL)  = NIRRE
             
        RVAL(1 ,IOFFL) = CD
        RVAL(2 ,IOFFL) = CM
        ! ADD 06-05-2015
        RVAL(4,IOFFL)  = CA 
C       ---- CALCULATE MARINE GROWTH ----	    
        ! ADD THE MARINE GROWTH
        IF (NDI.EQ.2.0)THEN ! USER DEFINED
        DIAM2          = DIAM1+(AMARINE(MLE)*2)
        ENDIF
C       ---------------------------------
        RVAL(3 ,IOFFL) = DIAM2
        RVAL(19,IOFFL) = NDI
        !RVAL(20,IOFFL) = AMARINE(IVAL(1,IOFFL)) ! FOR MARNIE GROWTH OLD VERSION
        RVAL(21,IOFFL) = OFFSELECT_X
        RVAL(22,IOFFL) = NFUNCTION
        RVAL(23,IOFFL) = ROUGH
        RVAL(24,IOFFL) = DIAM1 
        RVAL(25,IOFFL) = CSA 
        
        
    ! -------------------------------  FOR ERROR MESSAGE MARINE GROWTH ------------------------------------ 
      CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,3,STOPFUNCTION,
     1                        NGROWTH,AMARINE(MLE))
      IF (STOPFUNCTION.EQ.1) GOTO 101
      ENDDO
      

C     =====================================================================================================

   ! -------------------------------  FOR ERROR MESSAGE LOAD CASE ------------------------------- 
101    CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,1,STOPFUNCTION,
     1                        NGROWTH,GROWTH)
        
	IF(KOFFL.EQ.0) THEN
	  DEALLOCATE(RVAL,IVAL,AMARINE)
	  GOTO 117
      ENDIF
	! -----------------------
      
      IF (ISOLOP.EQ.5.OR.ISOLOP.EQ.6) THEN
      CALL STOREOFFSHOREELEMENT (NOFFL,IVAL,RVAL,LRPIN,"STOR")   
 !     GOTO 117
      ENDIF
  

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
      
      IF(NFATIGUEFREQUENCY .NE.1)THEN
      IF( KFTTD .NE. 1 )  NWAVESCAT = 1
      ENDIF
      
      NWAVESCATX = NWAVESCAT
      IF( NFATIGUEFREQUENCY .EQ. 1 )  NWAVESCATX = 0
      DO 7499 INDEXSCAT = 1,NWAVESCATX !LOOP OVER WAVE SCAATER DIRGRAM   
          
      CALL WRITEJOINTWAVEFORCEHEADER
      

      DO 7499 ILCAS = 1,IOFFL_CASE 
      
      Irandom = 1.0d0  

      
      DO 7500 ITIME = 1,NTIME

      RRV(1:NEQ) = 0.0D0 
      RRC(1:NEQ) = 0.0D0 
      
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') 
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')
      CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','READ')
      WW(1) = 0.0
      WW(2) = 0.0
      WW(3) = 0.0
      WW(4) = 0.0
      WW(5) = 0.0
      WW(6) = 0.0

      
      TOTALWAVEFORCE = 0.0D0    
      TOTALWAVEFORCEX = 0.0D0 
      TOTALWAVEFORCEY = 0.0D0 
      TOTALWAVEFORCEZ = 0.0D0 
      
	DO 127 IOFFL = 1,NOFFL
      
        WFX         = 0.
        WFY         = 0.
        WFZ         = 0.
        MLE         = IVAL(1 ,IOFFL)
        NWAVEELEVA  = MLE
        LWCASE(1)   = IVAL(2 ,IOFFL)
        LWCASE(2)   = IVAL(3 ,IOFFL)
        LWCASE(3)   = IVAL(4 ,IOFFL)
        LWCASE(4)   = IVAL(5 ,IOFFL)
        LWCASE(5)   = IVAL(6 ,IOFFL)
        LWCASE(6)   = IVAL(7 ,IOFFL)
        LWCASE(7)   = IVAL(8 ,IOFFL)
        LWCASE(8)   = IVAL(11 ,IOFFL)
        ILCN        = IVAL(9 ,IOFFL)
        ILCC        = IVAL(10,IOFFL)
        
        IF (ILCN.NE.ILCAS) GOTO 127 
        
        NWAVESRUFACESEGMENTOPT =  IVAL(12 ,IOFFL)
        NWAVESRUFACESEGMENT    =  IVAL(13 ,IOFFL) 
        
        NWAVESEGOPTi = IVAL(14,IOFFL)  
        NWAVESEGi    = IVAL(15,IOFFL)  
        IORRE        = IVAL(16,IOFFL)
        
        
        IF (NWAVESEGOPTi .EQ.1) NWAVEINVERVAL = 1
        IF (NWAVESEGOPTi .EQ.2) NWAVEINVERVAL = NWAVESEGi
        
        GROWTH      = RVAL(20,IOFFL)
        NDI         = RVAL(19,IOFFL)
        COA         = RVAL(4,IOFFL) 
        NFUNCTION   = RVAL(22,IOFFL)
        ROUGH       = RVAL(23,IOFFL)
        DIAM        = RVAL(3 ,IOFFL)
        DIAM1       = RVAL(24,IOFFL)
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
      CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,2,STOPFUNCTION,
     1                        NGROWTH,GROWTH)
8     FORMAT ('') ! BLANK

      ! --- STOP FUNCTION ---
      IF (STOPFUNCTION.EQ.1D0) STOP

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
      
C   ===========================================================================================   

      IF(ILCN.NE.ILCAS.AND.ILCC.NE.ILCAS) GOTO 127
        
C	LOC = LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR BODY LOAD)
	
	FRMFOC = 0.0D0
      

	DO  IGID = 1, NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 207
	ENDIF
	ENDDO
	EXIT
207   CONTINUE
      
      IF (NWAVEELEVA.EQ.348)THEN
          AAAAAA = 1
      ENDIF
      
      CX1= XYZ(1,MLE)
      CY1= XYZ(2,MLE)
      CZ1= XYZ(3,MLE)
      CX2= XYZ(4,MLE)
      CY2= XYZ(5,MLE)
      CZ2= XYZ(6,MLE)
      
!      WRITE(5028,8531) NWAVEELEVA,CX1,CY1,CZ1,CX2,CY2,CZ2    ! TOEY 10/2021
!8531  FORMAT(I5,2X,8F14.5)                                   ! TOEY 10/2021
      
      CALL WRITEJOINTWAVEFORCE(NWAVEELEVA,CX1,CY1,CZ1,CX2,CY2,CZ2,SEABED)

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
          CALL Airy_Wave_Number(GRAV,OMEGA,WDEPTH,RK)               
          Check_Wave_Number = (OMEGA**2.0d0)/(GRAV*tanh(RK*WDEPTH))

          IF (RK.LE.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,30)
          WRITE (*,31)
30        FORMAT ('************** OFFSHORE WARNING MESSAGE ( AIRY WAVE THEORY FOR FRAME ) **************')
31        FORMAT ('- PLEASE CHECK WAVE PERIOD AND THEORY CONDITION ( WAVE NUMBER < 0 )')
          ENDIF
          RAMDA(ILCAS) = 2.0D0*PI/RK    
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
          
         
          CALL STOKE_COEFFICIENT(RATIO,A11,A13,A15,A22,A24,A33,A35,A44,A55,B22,B24,B33,B35,B44,B55,
     1                             C1X,C2X,C3X,C4X)
          ! REMARK
          ! 
          ! F PRRAMETER == B PARAMETER
          ! G PRRAMETER == A PARAMETER
          ! C PRRAMETER == C PARAMETER

          RKx=2.0d0*22.0/7.0d0/RAMDA(ILCAS)
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
          
          
          IF (RK.LE.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,32)
          WRITE (*,33)
32        FORMAT ('************** OFFSHORE WARNING MESSAGE ( STOKE FIFTH ORDER THEORY FOR FRAME) **************')
33        FORMAT ('- PLEASE CHECK WAVE PERIOD AND THEORY CONDITION ( WAVE NUMBER < 0 )')
          ENDIF
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
      
      !WRITE(*,*) 'WAVE CREAST', WNU
     
      ELSEIF (IWAVE == 3 )THEN
          
      CALL STREAMWAVEFUNCTION (WVHIGHT,PERIOD,WDEPTH,TIME,WNU,GRAV,ORDER)
      
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
      
      ELSEIF (IORRE.EQ.2)THEN ! RANDOOM WAVE 12-01-2021
      ! RANDOM WAVE 12--01-2021 
      ! CALCULATE WAVE ELEVATION FOR SIDE OOF ELEMENT    
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
      
      YL_STABILITY_STATIC(ILCAS)  = WNU
      YL_STABILITY_DYNAMIC(ITIME) = WNU
      
      
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
      DIAM = DIAM1
      ENDIF ! (ENDIF : USER DEFINED NORMINAL DIAMETER ) 
      ! -----------------------------------------------------------------------  
      
      !MEMGD = IGIDM(IGID)
      
      CHXR =  XYZ(1,NWAVEELEVA) 
      CHYR =  XYZ(2,NWAVEELEVA) 
      CHZR =  XYZ(3,NWAVEELEVA) 
      
      CHXR2 =  XYZ(4,NWAVEELEVA) 
      CHYR2 =  XYZ(5,NWAVEELEVA) 
      CHZR2 =  XYZ(6,NWAVEELEVA) 

      
!      WRITE(5027,11259) -MLE , CHXR , CHYR , CHZR  ! TOEY 10/2021
!11259 FORMAT(I5,4F12.4)                            ! TOEY 10/2021
      
      	
C     ------------------------------------------------------------
C     LOOP OVER Trapizoidal Rule TO DET. STIFFNESS & FORCE VECTOR
C     ------------------------------------------------------------
C     BY  CHANA SINSABVARODOM  
C     15 MARCH 2016
C     ----------------------------------------------------------
C     Trapizoidal Rule
C     ----------------------------------------------------------
      
!      WRITE(5029,78501)                                                                                            ! TOEY 10/21
!78501 FORMAT('ELEMENT  STRIP          ELEVATION            INTERVAL         WAVEFORCE              DIAMETER')      ! TOEY 10/21
      
!      WRITE(5031,*)' MLE,   WAVEFORCEXELE1, WAVEFORCEXELE2,    WAVEFXI,         HGM,        CENTROID '  ! TOEY 10/21
!      WRITE(5032,*)'  MLE,  INDEXSTRIP,      WAVEFXI,             WAVEFYI,            WAVEFZI'          ! TOEY 10/21
!      WRITE(5033,*)'MLE,INDEXTAPPER,RREX1,RREX2,RMEX1,RMEX2,WAVEFXI,WAVEFYI,WAVEFYI,CENTROID'           ! TOEY 10/21
!      if( KFTTD .NE.1 ) write(*,*)
      
      HEADERTRAPIZOIDALRULEWRITE = 1    
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      !====================================
      ! UNDER WATER MEMBER CLASSIFICATION
      !====================================
      
      If( NGRAV.EQ.3 ) THEN
      
      BEGINMEMGERELEV =  XYZ(3,MLE)-SEABED   
      ENDMEMBERELEV   =  XYZ(6,MLE)-SEABED
      
      ELSEIF( NGRAV.EQ.2 ) THEN
      
      BEGINMEMGERELEV =  XYZ(2,MLE)-SEABED   
      ENDMEMBERELEV   =  XYZ(5,MLE)-SEABED
          
      ENDIF
      !====================================
      ! CLASSIFICATION
      IF( BEGINMEMGERELEV .LE. YL .AND. ENDMEMBERELEV .LE. YL  ) THEN
          
          UNDWMEMBER = 'UNW'
          
      ELSEIF( BEGINMEMGERELEV .GT. YL .AND. ENDMEMBERELEV .LT. YL  ) THEN
          
          UNDWMEMBER = 'BTW'
          
      ELSEIF( BEGINMEMGERELEV .LT. YL .AND. ENDMEMBERELEV .GT. YL  ) THEN
          
          UNDWMEMBER = 'BTW'
          
      ELSEIF( BEGINMEMGERELEV .GT. YL .AND. ENDMEMBERELEV .GT. YL  ) THEN
          
          UNDWMEMBER = 'UPW'
          
      ENDIF
      
      !====================================
      
!      WRITE(*,744) BEGINMEMGERELEV , ENDMEMBERELEV  , YL , UNDWMEMBER
      
744   FORMAT( 'ELEVA', 2X, 3F10.4 ,2X,3A  )  
      

      ALENGTHELEMENT = SQRT( DiffXELE**2D0 + DiffYELE**2D0 + DiffZELE**2D0 )
      
      If(NGRAV.EQ.3 ) THEN
          
          IF( DiffXELE.EQ.0.0D0 .AND. DiffYELE.EQ.0.0D0 )THEN
          
      ! =========================================
      ! Vertical Structure IN GRAVITY Z DIRECTION
      ! =========================================
      
      INTERVALTPZ = 500      
      IF(ALENGTHELEMENT.GT.300) THEN
      INTERVALTPZ = INT(ALENGTHELEMENT)*2.0D0  
      ELSEIF(ALENGTHELEMENT.GT.120.AND.ALENGTHELEMENT.LE.300) THEN
      INTERVALTPZ = INT(LENGTHELEMENT)*4.0D0      
      ELSEIF(ALENGTHELEMENT.GT.80.AND.ALENGTHELEMENT.LE.120) THEN
      INTERVALTPZ = 200        
      ELSEIF(ALENGTHELEMENT.GT.60.AND.ALENGTHELEMENT.LE.80) THEN
      INTERVALTPZ = 100     
      ELSEIF(ALENGTHELEMENT.GT.30.AND.ALENGTHELEMENT.LE.60) THEN
      INTERVALTPZ = 80   
      ELSEIF(ALENGTHELEMENT.GT.10.AND.ALENGTHELEMENT.LE.30) THEN
      INTERVALTPZ = 60  
      ELSEIF(ALENGTHELEMENT.GT.5.AND.ALENGTHELEMENT.LE.10) THEN
      INTERVALTPZ = 40 
      ELSEIF(ALENGTHELEMENT.GT.4.AND.ALENGTHELEMENT.LE.5) THEN
      INTERVALTPZ = 10 
      ELSEIF(ALENGTHELEMENT.GT.3.AND.ALENGTHELEMENT.LE.4) THEN
      INTERVALTPZ = 5 
      ELSEIF(ALENGTHELEMENT.LE.5) THEN
      INTERVALTPZ = 20 
      ENDIF
      
          ELSE
              INTERVALTPZ = NWAVEINVERVAL 
              
              IF     ( UNDWMEMBER .EQ. 'UNW' ) THEN
              INTERVALTPZ = NWAVEINVERVAL
              ELSEIF ( UNDWMEMBER .EQ. 'BTW' ) THEN
                  
              IF ( NWAVESRUFACESEGMENTOPT .EQ. 1) THEN
              INTERVALTPZ = 100
              ELSEIF( NWAVESRUFACESEGMENTOPT .EQ. 2 ) THEN
              INTERVALTPZ = NWAVESRUFACESEGMENT
              ENDIF
              
              ELSEIF ( UNDWMEMBER .EQ. 'UPW' ) THEN
              INTERVALTPZ = NWAVEINVERVAL 
              ENDIF
                  
              
          ENDIF
          
      ELSEIF(NGRAV.EQ.2 ) THEN
          
          IF( DiffXELE.EQ.0.0D0 .AND. DiffZELE.EQ.0.0D0 )THEN
              
          ! =========================================
          ! Vertical Structure IN GRAVITY Y DIRECTION
          ! =========================================
          
          INTERVALTPZ = 500      
          IF(ALENGTHELEMENT.GT.300) THEN
          INTERVALTPZ = INT(ALENGTHELEMENT)*2.0D0  
          ELSEIF(ALENGTHELEMENT.GT.120.AND.ALENGTHELEMENT.LE.300) THEN
          INTERVALTPZ = INT(LENGTHELEMENT)*4.0D0      
          ELSEIF(ALENGTHELEMENT.GT.80.AND.ALENGTHELEMENT.LE.120)  THEN
          INTERVALTPZ = 200        
          ELSEIF(ALENGTHELEMENT.GT.60.AND.ALENGTHELEMENT.LE.80)   THEN
          INTERVALTPZ = 100     
          ELSEIF(ALENGTHELEMENT.GT.30.AND.ALENGTHELEMENT.LE.60)   THEN
          INTERVALTPZ = 80   
          ELSEIF(ALENGTHELEMENT.GT.10.AND.ALENGTHELEMENT.LE.30)   THEN
          INTERVALTPZ = 60  
          ELSEIF(ALENGTHELEMENT.GT.5.AND.ALENGTHELEMENT.LE.10)    THEN
          INTERVALTPZ = 40 
          ELSEIF(ALENGTHELEMENT.GT.4.AND.ALENGTHELEMENT.LE.5)     THEN
          INTERVALTPZ = 10 
          ELSEIF(ALENGTHELEMENT.GT.3.AND.ALENGTHELEMENT.LE.4)     THEN
          INTERVALTPZ = 5 
          ELSEIF(ALENGTHELEMENT.LE.5) THEN
          INTERVALTPZ = 2 
          ENDIF
      
          ELSE
          
          ! =========================================
          ! Incline Structure IN GRAVITY Y DIRECTION
          ! =========================================
      
              INTERVALTPZ = NWAVEINVERVAL 
              IF     ( UNDWMEMBER .EQ. 'UNW' ) THEN
              INTERVALTPZ = NWAVEINVERVAL
              ELSEIF ( UNDWMEMBER .EQ. 'BTW' ) THEN
                  
              IF ( NWAVESRUFACESEGMENTOPT .EQ. 1) THEN
              INTERVALTPZ = 100
              ELSEIF( NWAVESRUFACESEGMENTOPT .EQ. 2 ) THEN
              INTERVALTPZ = NWAVESRUFACESEGMENT
              ENDIF
              
              ELSEIF ( UNDWMEMBER .EQ. 'UPW' ) THEN
              INTERVALTPZ = NWAVEINVERVAL
              ENDIF
          
           ENDIF
          

      ENDIF

      NHorizontal = 0
      
      NXWAVEFRAME = 1
      
      SELECTCASE(NGRAV)
      CASE(1)
        DIFFINTERVAL = DiffXELE/INTERVALTPZ 
      CASE(2)
        DIFFINTERVAL = DiffYELE/INTERVALTPZ 
        
       If( DiffYELE .EQ. 0.0D0 ) THEN
          
        DIFFINTERVAL = 1.0D0
        
        NHorizontal = 1 
        
        INTERVALTPZ = 30
        
       ENDIF  
              
      CASE(3)
        
      If( DiffZELE .EQ. 0.0D0 ) THEN
          
        DIFFINTERVAL = 1.0D0
        
        NHorizontal = 1 
        
        INTERVALTPZ = 30
        
      ENDIF  

      ENDSELECT
      
      WAVEFXI = 0
      SUMWAVEELEMENTFORCE = 0
      
      EELEMENTWFX = 0
      EELEMENTWFY = 0
      EELEMENTWFZ = 0
      EELEMENTWMX = 0
      EELEMENTWMY = 0
      EELEMENTWMZ = 0
     
      NEXTRAWAVE = 0
      
      TOTALWAVEFORCEX_X_CENTROID = 0.0D0 
      TOTALWAVEFORCEY_X_CENTROID = 0.0D0
      TOTALWAVEFORCEZ_X_CENTROID = 0.0D0
      TOTALFC = 0.0D0
      
      WW(1:6) = 0.0D0
      
      INTERVALTPZ = 10d0
      DO 307 IGR = 1,INTERVALTPZ
          
      WRITEDATATPZ = IGR
              
              
      WAVEFORCEXELE1 = 0.0D0
      WAVEFORCEYELE1 = 0.0D0
      WAVEFORCEZELE1 = 0.0D0
      
      WAVEFORCEXELE2 = 0.0D0
      WAVEFORCEYELE2 = 0.0D0
      WAVEFORCEZELE2 = 0.0D0

          DO 507 IELETPZ = 1,2
              
              
              
4756      CONTINUE

C     LOCATION ALONG ELEMENT AXIS

      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      ELEMENTFXYZ(1:6) = XYZ(1:6,MLE)
      
      
      if(IELETPZ.EQ.1)THEN
          IGRii = IGR - 1
      ELSEIF(IELETPZ.EQ.2)THEN
          IGRii = IGR 
      ENDIF
      
      SELECTCASE(NGRAV)
      CASE(1)
        DIFFINTERVAL = DiffXELE/INTERVALTPZ 
      CASE(2)
        DIFFINTERVAL = DiffYELE/INTERVALTPZ 
      CASE(3)
        DIFFINTERVAL = DiffZELE/INTERVALTPZ        
      ENDSELECT

      XR = XYZ(1,MLE) + DiffXELE/INTERVALTPZ*IGRii   !VR(1)*BXD
      YR = XYZ(2,MLE) + DiffYELE/INTERVALTPZ*IGRii   !VR(2)*BXD
      ZR = XYZ(3,MLE) + DiffZELE/INTERVALTPZ*IGRii   !VR(3)*BXD
      
      !WRITE(*,*)ZR,IGRii,MLE
      
      
      SELECTCASE(NGRAV)
      CASE(1)
        RLEV = XR   ! RLEV IS THE RUNNING ELEVATION
      CASE(2)
        RLEV = YR   ! RLEV IS THE RUNNING ELEVATION
      CASE(3)
        RLEV = ZR   ! RLEV IS THE RUNNING ELEVATION
      ENDSELECT
      
C     GAUSS COORDINATE IN WAVE DIRECTION SYSTEM	
      H1R = VH1(1)*XR + VH1(2)*YR + VH1(3)*ZR
      HGR = VGV(1)*XR + VGV(2)*YR + VGV(3)*ZR
      H2R = VH2(1)*XR + VH2(2)*YR + VH2(3)*ZR
      
C     GAUSS COORDINATE RELATIVE TO WAVE REFERENCE POINT (IN WAVE DIRECTION SYSTEM)
      H1R =  H1R - H1REF
      HGM =  HGR - HGREF
      H2M =  H2R - H2REF

      INDEXTAPPER = IGR
      
      !WRITE(*,*) 'CHANA',HGM
      
      
      ! ----------------- GENERATE DRAG AND INTIA COEFFICIENT -----------------
      ! INPUT  >> NDI = 1 : AUTOMATIC NORMINAL DIAMETER
      !        >> NDI = 2 : USER DEFINED NORMINAL DIAMETER
      ! OUTPUT >> CD  = DRAG COEFFICIENR
      !        >> CM  = INTIA COEFFICIENT
      ! ****** AUTOMATIC GENERATE BASE ON GAUSS POINT *****
      
      IF (NDI.EQ.1.0D0)THEN ! AUTOMATIC NORMINAL DIAMETER
          
      CALL TAPPEROFFSHORETAPIZOIDAL (DIAMOUTORIGINAL,IVAL(1,IOFFL),XYZ,INDEXTAPPER,INTERVALTPZ)    
           DIAMCAL    = DIAMOUTORIGINAL+GROWTH
           DIAM       = DIAMCAL*2.0D0
           
      CALL DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM,PERIOD,WVHIGHT,RK
     1                          ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                          ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
      
      !==============================================================================================
      If( DIAM .GT. 10000.0d0  ) then
          
          !Write(*,*) '*****************************************************************'
          !write(*,*) ' Diameter = ' , DIAM
          !write(*,*) ' Member = ' ,  MLE

           !Write(*,*) '*****************************************************************'
           !Write(*,*) ' Diameter Error, the value of wind force is too large '
           !Write(*,*) '*****************************************************************'

           Call Wave_Diameter(MLE,RadiousOUTORIGINAL)
           DIAMCAL    = RadiousOUTORIGINAL+GROWTH
           DIAM       = DIAMCAL*2.0D0

      Endif
      !==============================================================================================         
      ENDIF   ! (ENDIF : AUTOMATIC NORMINAL DIAMETER ) 
      ! -----------------------------------------------------------------------  

      CALL MARINE_GROWTH_CAL (RVAL(21,IOFFL),RLEV,GROWTH,DENSITY,IORRE) 
      DIAM = DIAM + (GROWTH*2D0)

      IF (LWCASE(2).EQ.1.0D0.AND.LWCASE(1).EQ.0.0)THEN
      YL = WDEPTH
      ENDIF   
      
      KFLACAL = 0
      
      WW(1:6) = 0.0D0

      DO 407 LCASE = 1,8
      ICASE = LWCASE(LCASE)
      IF(ICASE.EQ.0) GOTO 407
     
      ! SELECT CASE WITH LOAD FOR ANALYSIS ( LWCASE )
      
      SELECTCASE(LCASE)
          
      CASE(1,2,3) !WAVE LOAD - CURRENT LOAD - WATER TIDE LOAD    
      
      IF(KFLACAL.EQ.1) GOTO 407  !IF CURRENT&TIDAL LOAD ALREADY CALCULATE TOGETHER WITH WAVE LOAD, THEN JUMP
      
      !write(*,*)HGR
      
      IF(LWCASE(2).EQ.0) VWIND0 = 0.0D0
      IF(LWCASE(3).EQ.0) VTIDE  = 0.0D0
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407

          IF(WLEV.GT.(YL+abs(DIFFINTERVAL))) GOTO 407
           BARING(1:3) = VH1(1:3)
           VERTICALDIR(1:3) = VGV(1:3)

          IF ( IWAVE . NE . 6D0 ) THEN

           IF (NEXTRAWAVE .EQ. 1) GOTO 407
          YLACHECK = YL + DIFFINTERVAL
          IF ( HGM .LT. YLACHECK .AND.  HGM.GT.YL  ) THEN
             
          HGMWAVEINP = YL
          
          NEXTRAWAVE = 1
         ELSEIF ( AGAPELE .GT. (YL+DIFFINTERVAL)   ) THEN    
         GOTO 407
         ELSE
          HGMWAVEINP = HGM  
         ENDIF
             
         IF (IELETPZ .EQ.1 ) THEN   
         IF ( HGM .LE. YL ) STRIPWATER1 = 'UNS'
         IF ( HGM .GT. YL ) STRIPWATER1 = 'UPS'
         ELSEIF (IELETPZ .EQ.2 ) THEN 
         IF ( HGM .LE. YL ) STRIPWATER2 = 'UNS'
         IF ( HGM .GT. YL ) STRIPWATER2 = 'UPS'
         ENDIF
         
         ! STRIPWATER
         
         !CHARACTER*3 STRIPWATER1 
         !CHARACTER*3 STRIPWATER2
      
         
         IF ( IELETPZ.EQ.1 ) ELEVATIONSTRIP1 = HGM !HGMWAVEINP
         IF ( IELETPZ.EQ.2 ) ELEVATIONSTRIP2 = HGM !HGMWAVEINP
         
         NCALCULATIONFORCE = 0 ! INITIAL SET  0 = CAL WAVE FORCE,  1 = NOT CAL 
         

         !==================================================
         IF( DIFFINTERVAL .LT. 0 )THEN ! REVERT CONDITION
         !==================================================
         !HGMWAVEINP = 65.8422D0
             
           YLACHECK = YL - DIFFINTERVAL
          IF ( HGM .LT. YLACHECK .AND.  HGM.GT.YL .AND. IELETPZ .EQ.1 ) THEN
          HGMWAVEINP = YL
          ENDIF
          ENDIF
          !==================================================
          
          
      CALL WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGMWAVEINP,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,IORRE,IRWAVE)
      
!      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
!      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
!      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2
       WFX = WH1
       WFY = WGV
       WFZ = WH2
       
       IF(HGMWAVEINP.GT.YL) THEN
           
       WFX = 0.0D0
       WFY = 0.0D0
       WFZ = 0.0D0
           
      ENDIF

!       IF (IELETPZ .EQ.1 )Write(*,786)  WFX ,HGM , YL ,IGR , HGMWAVEINP
!       IF (IELETPZ .EQ.2 )Write(*,787)  WFX ,HGM , YL ,IGR , HGMWAVEINP
       
786    FORMAT('BEGIN',2X, f15.4 ,f10.4, f10.4 ,I5, f10.4)
787    FORMAT('END  ',2X, f15.4 ,f10.4, f10.4 ,I5, f10.4)
       
       
       !CALL WAVE_ANGLE_VECTOR( WFX , WFY , WFZ , VH1(1) )
       
      KFLACAL = 1
      
       !WRITEDATATPZ = 0
      
      IF( IWAVE . EQ . 3 ) THEN
      
          APRECENT = 2.50D0
          
          SURFACEMAX = Yl
          
          AWAVE_ROTIO = HGM / Yl
          
          AFACTORWAVE = ( 1.0D0 + AWAVE_ROTIO * (APRECENT / 100.0D0 ) )

          
            WFX = WFX *AFACTORWAVE
            WFY = WFY *AFACTORWAVE
            WFZ = WFZ *AFACTORWAVE
            
      ELSEIF( IWAVE . EQ . 4 ) THEN     
            
          STCvalue = 1.45
          Endvalue =0.325

          
          AFACTORWAVE = STCvalue - ( HGM / Yl ) * Endvalue
          
          IF( HGM.GT.0.4D0*Yl.AND. HGM.LE.0.5D0*Yl ) THEN
              
              AFACTORWAVE = AFACTORWAVE *  0.990D0
          
          ENDIF
          
          IF( HGM.GT.0.5D0*Yl.AND. HGM.LT.0.75D0*Yl ) THEN
              
              AFACTORWAVE = AFACTORWAVE *  0.980D0
          
          ENDIF
          
          IF( HGM.GT.0.9D0*Yl ) THEN
              
              AFACTORWAVE = AFACTORWAVE *  1.030D0
          
          ENDIF
          
            WFX = WFX *AFACTORWAVE
            WFY = WFY *AFACTORWAVE
            WFZ = WFZ *AFACTORWAVE
          
      ENDIF
      
      
      ! write(5299,*) WFX ,  HGMWAVEINP
      

      SELECTCASE(NGRAV)
      CASE(2)
      WFFC = SQRT(WFX**2.0D0 + WFZ**2.0D0)
      CASE(3)
      WFFC = SQRT(WFX**2.0D0 + WFY**2.0D0)
      ENDSELECT
      
!      WRITE(5029,7851) NWAVEELEVA,IGR,HGM,DIFFINTERVAL,WFFC,DIAM,H1R   ! TOEY 10/21
!7851  FORMAT( I5,2X,I5,7F20.5)                                         ! TOEY 10/21
            
      
      IF (IGR .EQ. 1.and. IELETPZ.eq.1) THEN ! PLOT JOINT LOAD
                
          
      CALL WRITE_WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,VH1,VGV,VH2,YL,XYZ,MLE,IORRE,IRWAVE,
     1                ITIME,ILCAS)
      
           
      ENDIF    
      
      ELEVATION_X = H1R
      ELEVATION_Y = HGM
      
!      WRITE(791,791)  ELEVATION_Y ,',',WFX
791   FORMAT(F12.5,A,E12.5)
      
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
      
      ENDIF
      
C    ===================================================================================    
C    =============================== UNNECESSARY BY TOEY ===============================
C      CASE(4) ! TIDAL LOAD
C          WLEV = RLEV - SEABED
C          IF(WLEV.LT.0.0D0) GOTO 407
C          IF(WLEV.GT.PEAKWLEV) GOTO 407
C      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA(LCS),GRAV,RHOW,RHOA,
C     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
C     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF)
C      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH2(1)*WH2 
C      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH2(2)*WH2 
C      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH2(3)*WH2 
C    ===================================================================================
      
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
      
C    ===================================================================================    
C    =============================== UNNECESSARY BY TOEY ===============================      
C      CASE(6) !WAVE BREAKING LOAD  SURING
C          WLEV = RLEV - SEABED
C          IF(WLEV.LT.0.0D0) GOTO 407
C          IF(WLEV.GT.YL) GOTO 407
C      CALL WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
C     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS)
C      WFX = VH1(1)*WFF ; WFY = VH1(2)*WFF ; WFZ = VH1(3)*WFF ;
C    ===================================================================================
      
      CASE(7) !WIND LOAD 

          
          WLEV = RLEV - SEABED - WDEPTH
          
          WLEV = RLEV - SEABED - WDEPTH ! + WNU
          
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
             !  -----------------------------
      
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
      
            ELSEIF (NWINDO.EQ.2.OR.NWINDO.EQ.5)THEN ! DNV > POWER LAW PROFILE NWINDO >> 2, USER DEFIND NWINDO >> 5 
                 IF (NWINDO.EQ.5)THEN ! SAME EQUATION FROM USER DEFINED AND POWER LAW IN DNV
                 UH = UHD             ! UH  = MEAN WIND SPEED ON DNV 
                 ENDIF                ! UHD = MEAN WIND SPEED ON USER DEFINED
      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA(ILCAS),GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     2                        CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
     
            ELSEIF (NWINDO.EQ.3)THEN ! API > WIND PROFILE AND GUST
                

                IF (WLEV.LT.0.000000001) WLEV = 0.00D0
                IF (GRAV.GE.9.00D0.OR.GRAV.LE.10.00D0)THEN ! M
                C   = 0.0573D0*SQRT(1D0+0.0457D0*UHAPI/0.3048D0)
                UZ  = (UHAPI/0.3048D0)*(1+(C*LOG(WLEV/RHIGH)))
                ZI  = 0.06D0*(1+0.0131D0*UHAPI/0.3048D0)*((WLEV/RHIGH)**(-0.22))
                UZT = (UZ*(1-0.41D0*ZI*LOG(AVERAGE/3600D0)))*0.3048D0
                ELSEIF (GRAV.GE.30D0.OR.GRAV.LE.33D0)THEN ! FT
                C   = 0.0573D0*SQRT(1D0+0.0457D0*UHAPI)
                UZ  = UHAPI*(1+(C*LOG(WLEV/RHIGH)))
                ZI  = 0.06D0*(1+0.0131D0*UHAPI)*((WLEV/RHIGH)**(-0.22))
                UZT = UZ*(1-0.41D0*ZI*LOG(AVERAGE/3600D0))
                ENDIF
               !WRITE (300,1111) WLEV,UZT
1111           FORMAT (E12.5,E12.5)
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
      
      !==============================================================================================
      If( WFX .GT. 1000000000.0d0 .or. WFY .GT. 1000000000.0d0 .or. WFZ .GT. 1000000000.0d0 ) then
          
!          write(*,4566) WFX ,  ( WLEV  + WDEPTH ) , RLEV ,  WDEPTH, YL , WNU , WLEV , MLE
          
!           Write(*,*) '*****************************************************************'
!           Write(*,*) ' Wind Error, the value of wind force is too large '
!           Write(*,*) '*****************************************************************'
          
          Stop
          
      Endif
      !==============================================================================================
      
      !write(5299,4566) WFX ,  ( WLEV  + WDEPTH ) , RLEV ,  WDEPTH, YL , WNU , WLEV , MLE
4566  FORMAT(7F12.4,I6)    
          ENDSELECT
      
      BXW = 1
      
      ! sorn mark old place
       DO 4600 Nckoff = 1,NELO
           
          sorn = Nmemoffset(Nckoff)  
	   
          IF (MLE .eq. Nmemoffset(Nckoff)) THEN
             Nkeyoffset = 1
             EXIT
          ENDIF
           
4600	CONTINUE  
          
         IF (Nkeyoffset.ne.0) THEN
          
         Call Reduce_Offset_member (XYZ,XR,YR,ZR,MLE,WFX,WFY,WFZ,IELETPZ)          
         
         ENDIF
         
         Nkeyoffset = 0
         
      
      !DIFFINTERVAL
      !===================================================
      ! --- TOTAL BUNDARY WAVE FORCE FOR TRAPIZOIDAL RULE
      !===================================================
      
      If(IELETPZ.EQ.1)THEN
      
          WAVEFORCEXELE1  = WAVEFORCEXELE1  + WFX
          WAVEFORCEYELE1  = WAVEFORCEYELE1  + WFY
          WAVEFORCEZELE1  = WAVEFORCEZELE1  + WFZ
      
      ELSEIF(IELETPZ.EQ.2)THEN
      
          WAVEFORCEXELE2  = WAVEFORCEXELE2  + WFX
          WAVEFORCEYELE2  = WAVEFORCEYELE2  + WFY
          WAVEFORCEZELE2  = WAVEFORCEZELE2  + WFZ
      
      ENDIF
      
          WFX = 0.
          WFY = 0.
          WFZ = 0.
      
!=================================================================================	
407   CONTINUE  ! END OF LOAD CALCULATION                             
!=================================================================================
      IF (WLEV.GT.(YL+abs(DIFFINTERVAL)) .AND. LWCASE(7).NE. 1 ) THEN
      EXIT
      ENDIF
!=================================================================================
507         CONTINUE  ! LOOP OF BEGINNING AND ENDING OF INTEGRATING STRIP              
!=================================================================================   
            
                  
      !===================================================        
      ! HORIAONTAL FORCE             
      !===================================================
          
              !IF     ( UNDWMEMBER .EQ. 'UNW' ) THEN
              !INTERVALTPZ = 1
              !ELSEIF ( UNDWMEMBER .EQ. 'BTW' ) THEN
              !INTERVALTPZ = 6
              !ELSEIF ( UNDWMEMBER .EQ. 'UPW' ) THEN
              !INTERVALTPZ = 1
              !ENDIF
              
      AAAVALUE = YL+DIFFINTERVAL
      
      !IF(WLEV.LE.(YL+DIFFINTERVAL)) THEN
      
C      IF(LWCASE(7).EQ.1) GOTO 600 !Suppressed by Songsak Oct2019
       
      IF(UNDWMEMBER .EQ.'UNW'.OR.UNDWMEMBER .EQ.'BTW'.OR.LWCASE(7).EQ.1) THEN
          
      IF(STRIPWATER1.EQ.'UNS'.OR.STRIPWATER2 .EQ.'UNS'.OR.LWCASE(7).EQ.1) THEN
          
      IF(NCALCULATIONFORCE.EQ. 0 .OR.LWCASE(7).EQ.1) THEN
          
          
600   CONTINUE     

      NCALCULATIONFORCE = 1
      
!      WRITE(*,8964)ELEVATIONSTRIP1,ELEVATIONSTRIP2,YL,STRIPWATER1,STRIPWATER2
8964  FORMAT(3F12.5,2X,A3,2X,A3)
          
      !IF ( ELEVATIONSTRIP1.LE.Yl .AND. ELEVATIONSTRIP2.LE.YL .AND. NEXTRAWAVE .EQ. 0 ) THEN
      
      !IF ( ELEVATIONSTRIP1.LE.Yl .AND. ELEVATIONSTRIP2.LE.YL  ) THEN
          
      !=========================================== 
      !  WAVE FORCE IN EACH TRAPIZOIDAL SECTMENT
      !===========================================
          
       WAVEFXI = 0.0D0
       WAVEFYI = 0.0D0
       WAVEFZI = 0.0D0
          
      CALL FORCE_EACH_SEGMENT(XYZ,'X',MLE, WAVEFORCEXELE1 , WAVEFORCEXELE2 , WAVEFXI ,INTERVALTPZ )
      CALL FORCE_EACH_SEGMENT(XYZ,'Y',MLE, WAVEFORCEYELE1 , WAVEFORCEYELE2 , WAVEFYI ,INTERVALTPZ )
      CALL FORCE_EACH_SEGMENT(XYZ,'Z',MLE, WAVEFORCEZELE1 , WAVEFORCEZELE2 , WAVEFZI ,INTERVALTPZ )
      
      !WRITE(5042,3849)WAVEFORCEXELE1,WAVEFORCEXELE2,WAVEFXI,NWAVEELEVA     ! TOEY 10/21
!      WRITE(*,3849)WAVEFORCEXELE1,WAVEFORCEXELE2,WAVEFXI,MLE               ! TOEY 10/21
!3849  FORMAT('STIFF FORCE',2x, F15.4,2x, F15.4,2x, F15.4, I4)
            
      ! CUTTING FORCE
      CALL CUTTING_FORCE(XYZ,'X',MLE,WAVEFXI)
      CALL CUTTING_FORCE(XYZ,'Y',MLE,WAVEFYI)
      CALL CUTTING_FORCE(XYZ,'Z',MLE,WAVEFYI)
              
!     WAVEFXI = 0.5D0*XDIFFINTERVAL*( WAVEFORCEXELE1 + WAVEFORCEXELE2 )
!     WAVEFYI = 0.5D0*YDIFFINTERVAL*( WAVEFORCEYELE1 + WAVEFORCEYELE2 ) 
!     WAVEFZI = 0.5D0*ZDIFFINTERVAL*( WAVEFORCEZELE1 + WAVEFORCEZELE2 ) 
       
       !WRITE(*,98310)MLE,INDEXTAPPER,WAVEFXI,WAVEFYI,WAVEFZI
        
!       WRITE(5032,98310) NWAVEELEVA,INDEXTAPPER,WAVEFXI,WAVEFYI,WAVEFZI    ! TOEY 10/21
!98310  FORMAT(2I7,10F20.5)                                                 ! TOEY 10/21
      
      !========================= 
      !  WAVE FORCE SUMMATION
      !=========================
       
      EELEMENTWFX = EELEMENTWFX + WAVEFXI
      EELEMENTWFY = EELEMENTWFY + WAVEFYI
      EELEMENTWFZ = EELEMENTWFZ + WAVEFZI
      EELEMENTWMX = EELEMENTWMX + 0.0D0
      EELEMENTWMY = EELEMENTWMY + 0.0D0
      EELEMENTWMZ = EELEMENTWMZ + 0.0D0
                
      !=========================================== 
      !  CENTROID CALCULATION OF WAVE FORCE
      !===========================================

      CALL TRAPIZOIDALCENTROID_XXWAVE( WAVEFORCEXELE1,WAVEFORCEYELE1,WAVEFORCEZELE1,
     1                           WAVEFORCEXELE2,WAVEFORCEYELE2,WAVEFORCEZELE2,
     1                              MLE,XYZ,INDEXTAPPER,INTERVALTPZ,CENTROIDXX) 
      
      
      CALL TRAPIZOIDALCENTROID_YYWAVE( WAVEFORCEXELE1,WAVEFORCEYELE1,WAVEFORCEZELE1,
     1                           WAVEFORCEXELE2,WAVEFORCEYELE2,WAVEFORCEZELE2,
     1                              MLE,XYZ,INDEXTAPPER,INTERVALTPZ,CENTROIDYY) 
      
      
      CALL TRAPIZOIDALCENTROID_ZZWAVE(WAVEFORCEZELE1,WAVEFORCEZELE2,
     1                              MLE,XYZ,INDEXTAPPER,
     1                                  INTERVALTPZ,RCENTROIDZZ,CENTROIDHH) 
      
      
      ! WRITE(*,98310) MLE,INDEXTAPPER,CENTROID,CENTROIDMINOR,CENTROIDHH
      

      
      CALL ELEMENT_CENTROID_X(INTERVALTPZ,INDEXTAPPER,WAVEFXI,CENTROIDXX   ,'WRT',TOTALFORCE,CENTROIDOUT)
      CALL ELEMENT_CENTROID_Y(INTERVALTPZ,INDEXTAPPER,WAVEFYI,CENTROIDYY   ,'WRT',TOTALFORCE,CENTROIDOUT)
      CALL ELEMENT_CENTROID_Z(INTERVALTPZ,INDEXTAPPER,WAVEFZI,RCENTROIDZZ  ,'WRT',TOTALFORCE,CENTROIDOUT)
            
      
C	=====================================================      
C	 WAVE FORCE IN X-DIRECTION   
C	=====================================================   
            
      VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
      
	CALL SCALEN(VR,VR,ELN,3)

	CALL XFSECTION(KEG,MLE,1)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !
      
      WW(1) = WAVEFXI
      WW(2) = 0.0d0
      WW(3) = 0.0d0
      WW(4) = 0.0d0
      WW(5) = 0.0d0
      WW(6) = 0.0d0
      
      !write(*,*) CENTROIDXX
      
      !write(*,*)'xx',RATIOX*Alengthelements
      
C	-----------------------------------------------------
	DO I = 1,3
	RP   = WW(I)/1.0E-8
	RPB  = RP
      
	DISA = CENTROIDXX !RATIOX*Alengthelements  !AAL*ELN 
	DISB = DISA+1.0E-8
	IDIR = I
	IND = 1
      
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCC,LOC,LOE,IMOM,'CONT')
	CALL FRAMFIX(RP,RPB,DISA,DISB,ECR,ECS,ECT,IDIR,0,VR,
	1			   ELN,FIXF,RANG,LOC,LOE,IMOM,IPIN)
      

	DO J = 1,7
	DO K = 1,2
	FRMFOC(J,K) = FRMFOC(J,K) + FIXF(J,K)
!     WAVEREACFIX(MLE,J,K) = FRMFOC(J,K)
	ENDDO
	ENDDO
	
      ENDDO
      
C	=====================================================      
C	 WAVE FORCE IN Y-DIRECTION    
C	=====================================================   
            
      VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
      
	CALL SCALEN(VR,VR,ELN,3)

	CALL XFSECTION(KEG,MLE,1)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !
      
      WW(1) = 0.0d0
      WW(2) = WAVEFYI
      WW(3) = 0.0d0
      WW(4) = 0.0d0
      WW(5) = 0.0d0
      WW(6) = 0.0d0
      
      !write(*,*) CENTROIDXX
      
      !write(*,*)'xx',RATIOX*Alengthelements
      
C	-----------------------------------------------------
	DO I = 1,3
          
	 RP   = WW(I)/1.0E-8
	 RPB  = RP
       
	 DISA = CENTROIDYY !RATIOX*Alengthelements  !AAL*ELN 
	 DISB = DISA+1.0E-8
	 IDIR = I
	 IND = 1
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCC,LOC,LOE,IMOM,'CONT')
	CALL FRAMFIX(RP,RPB,DISA,DISB,ECR,ECS,ECT,IDIR,0,VR,
	1			   ELN,FIXF,RANG,LOC,LOE,IMOM,IPIN)

	DO J = 1,7
	DO K = 1,2
	FRMFOC(J,K) = FRMFOC(J,K) + FIXF(J,K)
!     WAVEREACFIX(MLE,J,K) = FRMFOC(J,K)
	ENDDO
	ENDDO
	
      ENDDO
      
C	=====================================================      
C	 WAVE FORCE IN Z-DIRECTION    
C	=====================================================   
            
      VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
      
	CALL SCALEN(VR,VR,ELN,3)

	CALL XFSECTION(KEG,MLE,1)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !
      
      WW(1) = 0.0d0
      WW(2) = 0.0d0
      WW(3) = WAVEFZI
      WW(4) = 0.0d0
      WW(5) = 0.0d0
      WW(6) = 0.0d0
      
      !write(*,*) CENTROIDXX
      
      !write(*,*)'xx',RATIOX*Alengthelements
      
C	-----------------------------------------------------
	DO I = 1,3
	RP   = WW(I)/1.0E-8
	RPB  = RP
	DISA = RCENTROIDZZ !RATIOX*Alengthelements  !AAL*ELN 
	DISB = DISA+1.0E-8
	IDIR = I
	IND = 1
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCC,LOC,LOE,IMOM,'CONT')
	CALL FRAMFIX(RP,RPB,DISA,DISB,ECR,ECS,ECT,IDIR,0,VR,
	1			   ELN,FIXF,RANG,LOC,LOE,IMOM,IPIN)

	DO J = 1,7
	DO K = 1,2
	FRMFOC(J,K) = FRMFOC(J,K) + FIXF(J,K)
!     WAVEREACFIX(MLE,J,K) = FRMFOC(J,K)
	ENDDO
	ENDDO
	
      ENDDO
 
      ENDIF
      ENDIF
      ENDIF

       WAVEFXI = 0.0D0
       WAVEFYI = 0.0D0
       WAVEFZI = 0.0D0
              
!===============================================================================================      
307         CONTINUE ! END LOOP OF TRAPIZOIDAL ELEMENT
!===============================================================================================     
      
      ! CALAULATION ELEMENT CENTROID      
      CALL ELEMENT_CENTROID_X(INTERVALTPZ,INDEXTAPPER,WAVEFZI,CENTROIDHH,'RED',EELEMENTWFX,CENTROIDOUTX)
      CALL ELEMENT_CENTROID_Y(INTERVALTPZ,INDEXTAPPER,WAVEFZI,CENTROIDHH,'RED',EELEMENTWFY,CENTROIDOUTY)      
      CALL ELEMENT_CENTROID_Z(INTERVALTPZ,INDEXTAPPER,WAVEFZI,CENTROIDHH,'RED',EELEMENTWFZ,CENTROIDOUTZ)
      
      
      ! CALAULATION ELEMENT CENTROID RATIO 
      CALL ELEMENT_RATIO(XYZ,'X',MLE,CENTROIDOUTX,RATIOX)
      CALL ELEMENT_RATIO(XYZ,'Y',MLE,CENTROIDOUTY,RATIOY)
      CALL ELEMENT_RATIO(XYZ,'Z',MLE,CENTROIDOUTZ,RATIOZ)
      
      ! CUTTING FORCE
      CALL CUTTING_FORCE(XYZ,'X',MLE,EELEMENTWFX)
      CALL CUTTING_FORCE(XYZ,'Y',MLE,EELEMENTWFY)
      CALL CUTTING_FORCE(XYZ,'Z',MLE,EELEMENTWFZ)
      
      DiffXELE = XYZ(4,MLE) - XYZ(1,MLE)
      DiffYELE = XYZ(5,MLE) - XYZ(2,MLE)
      DiffZELE = XYZ(6,MLE) - XYZ(3,MLE)
      
      Alengthelements = sqrt( DiffXELE**2.0d0 + DiffYELE**2.0d0 + DiffZELE**2.0d0 )

	IF (NLS.NE.0) THEN
          
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
      
      ENDIF
       
      
      SUMWAVEELEMENTFORCE = EELEMENTWFX

          WAVEFXI = 0.0D0    
          WAVEFYI = 0.0D0    
          WAVEFZI = 0.0D0    
          
      
           if(HGM.LT.YLMIN)THEN
            !WRITE(*,*)'SUMWAVE', SUMWAVEELEMENTFORCE
            ENDIF
            
            !WRITE(5030,88765) SUMWAVEELEMENTFORCE,MLE
            !WRITE(5030,88765) EELEMENTWFX,EELEMENTWFY,EELEMENTWFZ,MLE
            IF ( KFTTD .NE. 1 ) THEN
            !WRITE(*,88765) EELEMENTWFX,EELEMENTWFY,EELEMENTWFZ,MLE
            ENDIF
            !WRITE(*,88766) RATIOX,RATIOY,RATIOZ,MLE
88765       FORMAT('ELEWFORCE,FX,FY,FZ',3X,F12.2,2X,F12.2,2X,F12.2,3X,'ELEMENT',2X,I6)  
88766       FORMAT('CENTROIDRATIO',2X,F12.2,2X,F12.2,2X,F12.2,3X,'ELEMENT',2X,I6)
            
            TOTALWAVEFORCE = TOTALWAVEFORCE + SUMWAVEELEMENTFORCE
            
            TOTALWAVEFORCEX = TOTALWAVEFORCEX + EELEMENTWFX
            TOTALWAVEFORCEY = TOTALWAVEFORCEY + EELEMENTWFY
            TOTALWAVEFORCEZ = TOTALWAVEFORCEZ + EELEMENTWFZ


      CALL KEEP_WAVE_FORCE_REACTION(FRMFOC,MLE,ILCAS)
      NWAVEOPT = 1
      CALL TRAPIZOIDAL_REACTION_LCS( ILCAS, NOFFL , NWAVEOPT , 'WRT' )
      


      !ILC = 1
      KEG = 1  ! FOR FAME ELEMENTS
      NEF = 14 ! DEGREEOF  OF FREEDOM
      !CALL FIXSTOR(ILC,KEG,MLE,FRMFOC,NEF,1.0D0,'ADD','VARY')	!FORCE IN GLOBAL SYSTEM NOT LOCAL SUPPORT (IN ELEMENT LOCAL FOR FRAME ELEMENT) VARY
	!CALL FIXSTOR(ILC,KEG,MLE,FRMFOC,NEF,1.0D0,'ADD','CONT')
      !CALL FIXSTOR(1,1,MLE,FRMFOC,14,1.0D0,'ADD','VARY')
          
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)

      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + FRMFOC(J,I) !VARY LOAD
      IF (IEQ.NE.0.AND.ILCC.GT.0) RRC(IEQ) = RRC(IEQ) + FRMFOC(J,I) !CONSTANT LOAD
            
	ENDDO
      ENDDO
      
      CALL FRMBEP_OFFSHORE (FRMFOC,VR,KEG,MLE,RANG)
      
      EELEMENTWFX  = 0.0D0 
      EELEMENTWFY  = 0.0D0 
      EELEMENTWFZ  = 0.0D0 
      EELEMENTWMX  = 0.0D0 
      EELEMENTWMY  = 0.0D0 
      EELEMENTWMZ  = 0.0D0 
      
!     ==============================================
127   CONTINUE ! END OF TOTAL FRAME STRUCTURAL LOOP 
!     ==============================================      
      CALL TRAPIZOIDAL_REACTION_LCS( ILCAS, LCSX , NWAVEOPT , 'RED' )
!      IF (NWAVEOPT.EQ.1)  WRITE(*,89634) TOTALWAVEFORCEX,TOTALWAVEFORCEY,TOTALWAVEFORCEZ
!      IF (NWAVEOPT.EQ.1)  WRITE(5030,89634) TOTALWAVEFORCEX,TOTALWAVEFORCEY,TOTALWAVEFORCEZ
!      IF (NWAVEOPT.EQ.1)  WRITE(*,*)
89634      FORMAT('TOTALWAVEFORCE',2X,F12.3,2X,F12.3,2X,F12.3)
!      IF (NWAVEOPT.EQ.1)  WRITE(5030,*) 'TOTALWAVEFORCE',TOTALWAVEFORCE
!      

      IF(KOFFL.NE.0) THEN
        
          
        CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
        CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','WRIT') !STORE CONSTANT LOAD
        
        IREAD=1D0
        
        
C        ====== THIS PART ARE NOT PROVIDED ===== 07-12-2013
C        CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','WRIT') !STORE CONSTANT LOAD
      ENDIF      

7500  CONTINUE
      
7499  CONTINUE
      
	DEALLOCATE(GPL,GPW,BPG,BWG,RRV,RRC,RVAL,IVAL)
      	
117   CONTINUE

C	======================================================================
      IF (NBOFL.EQ.0.0) GOTO 120 ! BUOYANCY FORCE 
      
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
          
      
      DO 130 ILCAS = 1,IOFFL_CASE !LOOP OVER LOAD CASE NUBER
      
      
      DO 130 ITIME = 1,NTIME !LOOP OVER OFFSHORE TIME STEP

      RRV(1:NEQ) = 0.0D0 !INITIALIZE VARY OFFSHORE LOAD
      
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL')           ! CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ') ! READING RRV 

      
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

        
	FRMFOC = 0.0D0

	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR SELFWEIGHT LOAD)

      CALL FRAME_AREA (MLE,GROWTH,AREA,"UN") 
      
	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 210
	ENDIF
	ENDDO
	EXIT
210   CONTINUE

      
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
      
      
      IF (RLEV.GT.WATER_ELEVATION) BUOYANCYF = 0.0D0
      IF (RLEV.LE.WATER_ELEVATION) THEN
      DEPTH_ELE = WATER_ELEVATION - RLEV     
      BUOYANCYF = RHOW*GRAV*DEPTH_ELE*AREA
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
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
C     MUST BE ADDED LATER  FOR SMOOTH FRAME DIAGRAM
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCC,LOC,LOE,IMOM,'CONT')
	CALL FRAMFIX(RP,RPB,DISA,DISB,ECR,ECS,ECT,IDIR,0,VR,
	1			   ELN,FIXF,RANG,LOC,LOE,IMOM,IPIN)
	DO J = 1,7
	DO K = 1,2
	FRMFOC(J,K) = FRMFOC(J,K) + FIXF(J,K)
	ENDDO
	ENDDO
	
      ENDDO
C	-----------------------------------------------------

    	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

      IK = 0
	DO I = 1,2
	DO J = 1,7
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
      
!      DO JJJ = 1,1000
!      WRITE (4903,9999) EVL_RAN(JJJ)
!9999  FORMAT (E12.5)
!      ENDDO
C     ==================================================
C     OFFSHORE SPECTRUM DATA
C     ==================================================   
	IF(NRES.EQ.0.0D0) GOTO 118

	ALLOCATE(RVAL(25,NRES),IVAL(20,NRES),AMARINE(NELE))

C     ----- SET CONDITION PERVENT ERROR ------
      GROWTH =0.0D0
C     ----------------------------------------
      AMARINE = 0.
      GROWTH  = 0.0D0
      IF (NMARINE.NE.0)THEN
      REWIND(70)
        DO I = 1,NMARINE
        READ (70,*) MLE_MARINE,THICKNESS
        AMARINE(MLE_MARINE) = THICKNESS
        ENDDO
      ENDIF
      !WW(1:6) = 0.0D0

C     ========================= READ INPUT DATA AND STORE ON TEMPORARY VARIABLE =========================
      READ (ITI,*) 
	DO IOFFL = 1,NRES
         ! READ SPECTRUM DATA MODIFY BY TOEY 04/2018
          LWCASE = 0.
         READ (ITI,*) MLE,NFUNCTION,ROUGH,CA,CD,CM,CSA,NWAVESEGOPT,NWAVESEG,NSURFWSEGOPT,NSURFWSEGMT,NDI,DIAM1,NGROWTH,
     1                GROWTH,OFFSELECT_X,LWCASE(8),LWCASE(2),LWCASE(5),LWCASE(6),LWCASE(7),ILCN,ILCC
         !CALL SPECTURMPARAMETER (MLE,NFUNCTION,ROUGH,CD,CM,CSA,NDI,DIAM1,NGROWTH,GROWTH,OFFSELECT_X,LWCASE,ILCN,ILCC)
        IVAL(1 ,IOFFL)  = MLE
        IVAL(2 ,IOFFL)  = LWCASE(1)
        IVAL(3 ,IOFFL)  = LWCASE(2)
        IVAL(4 ,IOFFL)  = 0.0D0     ! LWCASE(3) WATER LEVEL # CANCLE
        IVAL(5 ,IOFFL)  = 0.0D0     ! LWCASE(4) TIDAL LOAD  # CANCLE
        IVAL(6 ,IOFFL)  = LWCASE(5)
        IVAL(7 ,IOFFL)  = LWCASE(6)
        IVAL(8 ,IOFFL)  = LWCASE(7)
        IVAL(11,IOFFL)  = LWCASE(8)
        IVAL(9  ,IOFFL) = ILCN
        IVAL(10 ,IOFFL) = ILCC
             
        RVAL(1 ,IOFFL) = CD
        RVAL(2 ,IOFFL) = CM
        
C       ---- CALCULATE MARINE GROWTH ----	         
        IF (NDI.EQ.2.0)THEN ! USER DEFINED
        DIAM2          = DIAM1+(GROWTH*2)
        ENDIF
C        ---------------------------------
        
        RVAL(3 ,IOFFL) = DIAM2
        RVAL(19,IOFFL) = NDI
        RVAL(20,IOFFL) = GROWTH
        RVAL(21,IOFFL) = OFFSELECT_X
        RVAL(22,IOFFL) = NFUNCTION
        RVAL(23,IOFFL) = ROUGH
        RVAL(24,IOFFL) = DIAM1 ! NORMINAL DIAMETER FOR WIND ANALYSIS
        RVAL(25,IOFFL) = CSA 
        
          ! -------------------------------  FOR ERROR MESSAGE MARINE GROWTH ------------------------------- 
        CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,1,STOPFUNCTION,
     1                        NGROWTH,GROWTH)
         IF (STOPFUNCTION.EQ.1) GOTO 102
         
         ENDDO
C     =====================================================================================================
         
102    CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,1,STOPFUNCTION,
     1                        NGROWTH,GROWTH)
         
	IF(KSPEC.EQ.0) THEN
	  DEALLOCATE(RVAL,IVAL)
	  GOTO 118 !IF NO OFFSHORE LOAD ANALYSIS NO NEED TO DO THE CALCULATION...JUST READ THE INPUT DATA
	ENDIF
	! -----------------------
	
C	--- DETERMINE GAUSS POINT POSITION AND WEIGHT ---	
      ! INFORMATION ABOUT THE GAUSS
      ! - THIS CONDITION USING 10 POINT
      ! - FOR TOP PART AND BOTTOM PART FORCE EQUAL = 0.0
	ALLOCATE(GPL(12),GPW(12),BPG(12),BWG(12),RRV(NEQ),RRC(NEQ))

	DO IGR = 1,12
	IF(IGR.EQ.1  ) GPL(IGR) = -1.0D0
	IF(IGR.EQ.12 ) GPL(IGR) =  1.0D0 
	IF(IGR.NE.1.AND.IGR.NE.12) GPL(IGR) =  GAUSP(IGR-1,12-2)
	IF(IGR.EQ.1  ) GPW(IGR) =  0.0D0
	IF(IGR.EQ.12 ) GPW(IGR) =  0.0D0 
	IF(IGR.NE.1.AND.IGR.NE.12) GPW(IGR) =  GAUSW(IGR-1,12-2)
      ENDDO
      
      !IF (ISOLOP.NE.1) THEN
         !IF(KSPEC.NE.0.AND.OPRES.EQ.1.AND.ICONTROLSPEC.EQ.1) THEN !CALL ALL OFFSHORE LOAD PARAMETERS
         !CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
         !ELSEIF (KSPEC.NE.0.AND.OPRES.EQ.2.AND.ICONTROLSPEC.EQ.1) THEN 
         !CALL SELECTGRAPH (1,NTIME,0,0)
         !ELSEIF (KSPEC.NE.0.AND.ICONTROLSPEC.EQ.0)THEN
         CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
         !ENDIF
      !ELSEIF (ISOLOP.EQ.1)THEN
      !   IF(KSPEC.NE.0.AND.ICONTROLSPEC.EQ.1) THEN !CALL ALL OFFSHORE LOAD PARAMETERS
      !   CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      !   ELSEIF (KSPEC.NE.0.AND.ICONTROLSPEC.EQ.0)THEN
      !   CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      !   ENDIF
      !ENDIF

      IF(NFATIGUEFREQUENCY .NE.1)THEN
      IF ( KSPEC.EQ.1 .AND. NWAVESCAT.EQ.0 ) NWAVESCAT = 1
      IF ( KSPEC.EQ.2 .AND. NWAVESCAT.EQ.0 ) NWAVESCAT = 1
      ENDIF
      
      DO 7502 INDEXSCAT = 1,NWAVESCAT !LOOP OVER WAVE SCAATER DIRGRAM     
      
      DO 7502 ILCAS = 1,LCS !LOOP OVER LOAD CASE NUBER
          
      DO 7502 ITIME = 1,NTIME !LOOP OVER OFFSHORE TIME STEP
          
      NWSPECTRUMPLOT = 1   ! WAVE SPECTRUM PLOT

      RRV(1:NEQ) = 0.0D0 !INITIALIZE VARY OFFSHORE LOAD
      RRC(1:NEQ) = 0.0D0 !INITIALIZE CONSTANT OFFSHORE LOAD
      
      !IF (OPRES.EQ.1) CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      !IF (OPRES.EQ.2) CALL SELECTGRAPH (2,0,ITIME,TIME)
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'SPEC','READ')
      
      ! STATIC CASE
      IF (ISOLOP.EQ.1.AND.ICONTROLSPEC.EQ.1) THEN
      ISPEC = 1
      ENDIF
      !DO 7501 ISPEC = 1,1 !NMOD
      !SNA = EIGVFREQ(ISPEC)

      
10    DO 128 IOFFL = 1,NRES
      
          WFX = 0.
          WFY = 0.
          WFZ = 0.
          
C     ------------------------------
C     CALL INPUT DATA FROM BACKUP
        MLE         = IVAL(1 ,IOFFL)
        LWCASE(1)   = IVAL(2 ,IOFFL)
        LWCASE(2)   = IVAL(3 ,IOFFL)
        LWCASE(3)   = IVAL(4 ,IOFFL)
        LWCASE(4)   = IVAL(5 ,IOFFL)
        LWCASE(5)   = IVAL(6 ,IOFFL)
        LWCASE(6)   = IVAL(7 ,IOFFL)
        LWCASE(7)   = IVAL(8 ,IOFFL)
        LWCASE(8)   = IVAL(11,IOFFL)
        ILCN        = IVAL(9 ,IOFFL)
        ILCC        = IVAL(10,IOFFL)
        
        OFFSELECT   = RVAL(21,IOFFL) 
        
        ! --- FOR CALCULATE AUTOMATIC NORMINAL DIIAMETER ---
        GROWTH      = RVAL(20,IOFFL)
        NDI         = RVAL(19,IOFFL) 
        
        ! --- FOR CALCULATION CD AND CM AUTOMATIC ---
        NFUNCTION = RVAL(22,IOFFL)
        ROUGH     = RVAL(23,IOFFL)
        DIAM      = RVAL(3 ,IOFFL)
        DIAM1     = RVAL(24,IOFFL)
        CSA       = RVAL(25,IOFFL)
        ! FROM READ DATIN.DAT FOR USER DEFINED VALUE
        IF (NFUNCTION.EQ.4.0D0)THEN
        CD        = RVAL(1 ,IOFFL)
        CM        = RVAL(2 ,IOFFL)
        ELSEIF (NFUNCTION.EQ.1.0.OR.NFUNCTION.EQ.2.0.OR.NFUNCTION.EQ.3.0D0)THEN
        CD        = 0 
        CM        = 0
        ENDIF
        ! --------------------------------------------
        
        
       CALL SELECTSPECTRUM (OFFSELECT,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
       IF (ISOLOP.EQ.5.AND.ICONTROLSPEC.EQ.1)THEN
       FREQ =  TIME
       ENDIF
       
      IF(NFATIGUEFREQUENCY .EQ.1)THEN
           WVHIGHT = HSCAT(INDEXSCAT)
           PERIOD  = PERIODSCAT(INDEXSCAT)
      ENDIF

       

      WARNING = GRAV ! WARINING : SENT GRAVITY TO XFINAS.FOR. USED FOR OFFSHORE WARNING MESSAGE
C   ================================================================== 
    

C   ======================= BLOCK DATA FOR CURRENT MODIFILE BY TOEY =======================
      CALL BLOCKDATAOFFSHORE (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                        VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                        VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,LWCASE,2,STOPFUNCTION,
     1                        NGROWTH,GROWTH)
18     FORMAT ('') ! BLANK

      ! --- STOP FUNCTION ---
      IF (STOPFUNCTION.EQ.1D0)THEN
      STOP
      ENDIF
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
C   ========================================================================================   

      IF(ILCN.NE.ILCAS.AND.ILCC.NE.ILCAS) GOTO 128
        
C	LOC = LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR BODY LOAD)
	
	FRMFOC = 0.0D0

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 208
	ENDIF
	ENDDO
	EXIT
208	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)

      H1REF = H1POS
      H2REF = H2POS
      HGREF = SEABED	     
      
      PXYZ(1,1) = XYZ(1,MLE)
      PXYZ(1,2) = XYZ(2,MLE)
      PXYZ(1,3) = XYZ(3,MLE)
      PXYZ(2,1) = XYZ(4,MLE)
      PXYZ(2,2) = XYZ(5,MLE)
      PXYZ(2,3) = XYZ(6,MLE)

	XMID = 0.5*(PXYZ(1,1) + PXYZ(2,1))
      YMID = 0.5*(PXYZ(1,2) + PXYZ(2,2))
      ZMID = 0.5*(PXYZ(1,3) + PXYZ(2,3))
      
      HM1 = VH1(1)*XMID + VH1(2)*YMID + VH1(3)*ZMID 
      HMG = VGV(1)*XMID + VGV(2)*YMID + VGV(3)*ZMID 
      HM2 = VH2(1)*XMID + VH2(2)*YMID + VH2(3)*ZMID

      H1M =  HM1 - H1REF 
      HGM =  HMG - HGREF 
      H2M =  HM2 - H2REF 

      VREW(1) = VH1(1)*VR(1) + VH1(2)*VR(2) + VH1(3)*VR(3)
      VREW(2) = VGV(1)*VR(1) + VGV(2)*VR(2) + VGV(3)*VR(3) 
      VREW(3) = VH2(1)*VR(1) + VH2(2)*VR(2) + VH2(3)*VR(3) 

      ! SET IWAVE = 1
      IWAVE = 1
      OMEGA = 2.0*PI/PERIOD      
      CALL NEWTON_RAPHSON(WDEPTH,GRAV,OMEGA,RK)
      RAMDA(ILCAS) = 2.0*PI/RK    
      AVAL = 1.0D0
      RATIO = 1.0D0   
      
      WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
      YL = WDEPTH + WNU + SEABED
      YLMIN = WDEPTH - WNU

      
      IF (NDI.EQ.2.0D0)THEN ! USER DEFINED NORMINAL DIAMETER
          CALL DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM1,PERIOD,WVHIGHT,RK
     1                        ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                        ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
      DIAM = DIAM1
      ENDIF ! (ENDIF : USER DEFINED NORMINAL DIAMETER ) 
      ! -----------------------------------------------------------------------  
     
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
      DO 308 IGR = 1,12

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
      CALL TAPPEROFFSHORE (DIAM3,MLE,XYZ)
           DIAM3(IGR) = DIAM3(IGR)+GROWTH
           DIAM       = DIAM3(IGR)*2.0D0
           IF (ICONTROLSPEC.EQ.1) TIME = 0.0
           CALL DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM,PERIOD,WVHIGHT,RK
     1                          ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                          ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
           IF (ICONTROLSPEC.EQ.1) TIME = FREQ
      ENDIF   ! (ENDIF : AUTOMATIC NORMINAL DIAMETER ) 
      ! -----------------------------------------------------------------------  
      
      !CALL WAVE_SURFACE_PROFILE_FREQ

      IF (LWCASE(2).EQ.1.0D0.AND.LWCASE(1).EQ.0.0)THEN
      YL = WDEPTH + SEABED
      ENDIF       
      KFLACAL = 0
      
      WW(1:6) = 0.0D0 
      
      DO 408 LCASE = 1,8
         
      ICASE = LWCASE(LCASE)
      
      IF(ICASE.EQ.0) GOTO 408
     
      ! SELECT CASE WITH LOAD FOR ANALYSIS ( LWCASE )
      
      SELECTCASE(LCASE)
      CASE(1,2,3) !WAVE LOAD - CURRENT LOAD - WATER TIDE LOAD    
      
      IF(KFLACAL.EQ.1) GOTO 408  !IF CURRENT&TIDAL LOAD ALREADY CALCULATE TOGETHER WITH WAVE LOAD, THEN JUMP
      
      IF(LWCASE(2).EQ.0) VWIND0 = 0.0D0
      IF(LWCASE(3).EQ.0) VTIDE  = 0.0D0
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 408
          IF(WLEV.GT.YL) GOTO 408
      CALL WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGMWAVEINP,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,IORRE,IRWAVE)
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2
      KFLACAL = 1
      
C    ===================================================================================    
C    =============================== UNNECESSARY BY TOEY ===============================
C      CASE(4) ! TIDAL LOAD
C          WLEV = RLEV - SEABED
C          IF(WLEV.LT.0.0D0) GOTO 408
C          IF(WLEV.GT.PEAKWLEV) GOTO 408
C      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA(LCS),GRAV,RHOW,RHOA,
C     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
C     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF)
C      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH2(1)*WH2 
C      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH2(2)*WH2 
C      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH2(3)*WH2 
C    ===================================================================================

      CASE(5,6) !WAVE BREAKING LOAD  PLUNGING&SURGING
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 408
          IF(WLEV.GT.YL) GOTO 408
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
          IF(WLEV.LT.0.0D0) GOTO 408
            IF (NSWIND.EQ.1)THEN ! Kaimal Spectrum
            call Kaimal_Spectrum(Z0,ALPHA,UHD,RHIGH,GRAV,RHOA,Cd,DIAM,WLEV,FREQ,WFF,SNA,FWIND)
            !Subroutine Kaimal_Spectrum(zo,Alpha,Umean,RHIGH,G,ROHAIR,Cd,Dia,Z,f,SFUU,SNA)
            ELSEIF (NSWIND.EQ.2)THEN ! Kaimal_IEC_Spectrum
            call   Kaimal_IEC_Spectrum (Z0,ALPHA,UHD,RHIGH,GRAV,RHOA,Cd,DIAM,WLEV,FREQ,WFF,SNA,FWIND)
            
            ELSEIF (NSWIND.EQ.3)THEN ! Von_Karman Spectrum
            call Von_Karman_Spectrum(Z0,ALPHA,UHD,RHIGH,GRAV,RHOA,Cd,DIAM,WLEV,FREQ,WFF,SNA,FWIND)
            
            ELSEIF (NSWIND.EQ.4)THEN ! Davenport Spectrum
            call Davenport_Spectrum(Z0,ALPHA,UHD,RHIGH,GRAV,RHOA,Cd,DIAM,WLEV,FREQ,WFF,SNA,FWIND)
            
            ELSEIF (NSWIND.EQ.5)THEN ! Eurocode Spectrum
            call   Eurocode_Spectrum(Z0,ALPHA,UHD,RHIGH,GRAV,RHOA,Cd,DIAM,WLEV,FREQ,WFF,SNA,FWIND)
            ENDIF
      CSA = 0.5      
      WFX = VWIND(1)*WFF*CSA ; WFY = VWIND(2)*WFF*CSA ; WFZ = VWIND(3)*WFF*CSA ;
      
      CASE (8)
      
      WLEV = RLEV - SEABE
      !IF(WLEV.LT.0.0D0) GOTO 408
      IF(WLEV.GT.YL) GOTO 408
         IF (IRWAVE.EQ.1)THEN
         CALL Pierson_Moskowitz_SPectrum (WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,WH1,WGV,SNA,SDG,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.2)THEN
         CALL Jonswap_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,TIME,WH1,WGV,SNA,SDG,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.3)THEN
         CALL Ochi_SPectrum (WVH1,PER1,SHAPE1,WVH2,PER2,SHAPE2,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,TIME,WH1,WGV,SNA,SDG,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.4)THEN
         CALL  Bretschneider_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,TIME,WH1,WGV,SNA,SDG,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.5)THEN
         CALL TMA_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,TIME,WH1,WGV,SNA,SDG,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.6)THEN
         CALL User_Define_SPectrum(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,TIME,WH1,WGV,SNA,SDG,TAP,UCURRENT)
         ENDIF
         WH2 = 0.0D0
         WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
         WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
         WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2
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
          
408   CONTINUE
C     -------------------------- 	

C     ------------- OUTPUT OFFSHORE FILE (OUTOFFSHORE.OUT) -----------------

C      CALL PRINTOUTOFFSHORE (ILCAS,MLE,NDI,IWAVE,RK,CD,CM,WVHIGHT,WNU,RAMDA(ILCAS),IGR,DIAM,LCS)
C      CALL PRINTOUTOFFSHORE1 (IGR,MLE,ILCAS,LWCASE,UZT,WFF,WDEPTH,WH1,WGV)
C     ----------------------------------------------------------------------

C	----------------------------------------------------------------------

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
C	-----------------------------------------------------
308   CONTINUE
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
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + FRMFOC(J,I) !VARY LOAD
      IF (IEQ.NE.0.AND.ILCC.GT.0) RRC(IEQ) = RRC(IEQ) + FRMFOC(J,I) !CONSTANT LOAD
	ENDDO
	ENDDO

128	CONTINUE

      IF(KSPEC.NE.0) THEN
        IF (ISOLOP.EQ.5.AND.ICONTROLSPEC.EQ.0) CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'SPEC','WRIT') !STORE VARY LOAD (SPEC)
        IF (ISOLOP.EQ.5.AND.ICONTROLSPEC.EQ.1) CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
        IF (ISOLOP.EQ.1.AND.ICONTROLSPEC.EQ.0) CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'SPEC','WRIT') !STORE VARY LOAD
        IF (ISOLOP.EQ.1.AND.ICONTROLSPEC.EQ.1) CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
C        ====== THIS PART ARE NOT PROVIDED ===== 07-12-2013
C        CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','WRIT') !STORE CONSTANT LOAD
      ENDIF   
      
      !IF (ISOLOP.EQ.5) CALL MUTIFORCE (ITIME,ILCAS,ISPEC,NOFFL) ! MUTIPLY MODE SHAPE
      IF (ISOLOP.EQ.1) CALL MUTIFORCE (ITIME,ILCAS,ISPEC,NRES)
7501  CONTINUE   


7502  CONTINUE

	DEALLOCATE(GPL,GPW,BPG,BWG,RRV,RRC,RVAL,IVAL)
	
118   CONTINUE

	RETURN

      RETURN
      END

           