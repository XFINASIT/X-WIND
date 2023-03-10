      ! DEAD AND SELF WEIGHT LOAD
       SUBROUTINE FRAMFOC_OFFSHORE (LM,XYZ,IGSET,PROPG,MTSET,PROPM,IGIDM,
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
      
      DIMENSION WW_MARINE(6)
      
	ALLOCATABLE GPL(:),GPW(:),BPG(:),BWG(:)
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:),AMARINE(:)
      
	PI = 3.141592654
      WAVEREACFIX = 0.0D0
      
      NWELEMENT = NELE
      NSW      = MLOAD(1)
      NUF      = MLOAD(2)
      NMARINE  = MLOAD(3)
      NOFFL    = MLOAD(4)
      NBOFL    = MLOAD(8)
      
C	======================================================================
C	SELF WEIGHT LOADFRAMFOC
	IF(NSW.NE.0) READ (ITI,*)
	DO 110 ISW = 1,NSW

	CALL CLEARA (R,NEQ)
	FRMFOC = 0.0D0

	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR SELFWEIGHT LOAD)
	READ(ITI,*) MLE,GX,GY,GZ,NMARINE,NMARINE_GROUP,ILCN,ILCC
      
      IF (NMARINE.NE.0)THEN
      CALL MARINE_GROWTH_WEIGHT (MLE,IGIDM,XYZ,NMARINE_GROUP,GX,GY,GZ,WW_MARINE,ARATIO,BRATIO)
      ENDIF
      
	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 200
	ENDIF
	ENDDO
	EXIT
200	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)


	INAME(1:4) = [5,0,1,KEG] !XSEF
	CALL ICONC(INAME,NAMEI)
	CALL MINTFIL(NAMEI,NGAS,1,3,0)

	CALL XFSECTION(KEG,MLE,1)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG  ,1,5 ,0) !
	CALL MRELFIL(NAMEI,AREAD1,1,61,0) !Equivalent AREA-Density I

	CALL XFSECTION(KEG,MLE,NGAS)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,AREAD2,1,61,0) !Equivalent AREA-Density J

	UWT1 = AREAD1
	UWT2 = AREAD2 
      
	WW(1) = GX*UWT1
	WW(2) = GX*UWT2
	WW(3) = GY*UWT1
	WW(4) = GY*UWT2
	WW(5) = GZ*UWT1
	WW(6) = GZ*UWT2
      
      
C	-----------------------------------------------------
      IF (NMARINE.EQ.0) THEN
	DO I = 1,3
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	DISA = 0.0
	DISB = ELN
	IDIR = I
	IND  = 1

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
	ENDDO
	ENDDO
      
      ENDDO
C	-----------------------------------------------------
      ELSEIF (NMARINE.EQ.1) THEN
      DO I = 1,3
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	DISA = 0.0
	DISB = ELN
	IDIR = I
	IND  = 1

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
	ENDDO
	ENDDO
      
      ENDDO
      
      ! MARINE WEIGHT
      DO I = 1,3
	RP   = WW_MARINE(2*I-1)
	RPB  = WW_MARINE(2*I-0)
	DISA = ELN*ARATIO
	DISB = ELN*BRATIO
	IDIR = I
	IND  = 1

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
	ENDDO
	ENDDO
      
      ENDDO
C	-----------------------------------------------------
      ELSEIF (NMARINE.EQ.2) THEN

      DO I = 1,3
	RP   = WW_MARINE(2*I-1)
	RPB  = WW_MARINE(2*I-0)
	DISA = ELN*ARATIO
	DISB = ELN*BRATIO
	IDIR = I
	IND  = 1

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
	ENDDO
	ENDDO
      
      ENDDO
      
      
      ENDIF
C	-----------------------------------------------------
	CALL FRMBEP_OFFSHORE (FRMFOC,VR,KEG,MLE,RANG)
      
      
	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF
      
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) R(IEQ) = R(IEQ) + FRMFOC(J,I)
	ENDDO
	ENDDO

	CALL LDASEM_OFFSHORE (R)

110	CONTINUE

C	======================================================================
C	DEAD LOAD
	IF(NUF.NE.0) READ (ITI,*)
	DO 111 IUF = 1,NUF

	CALL CLEARA (R,NEQ)
      
	FRMFOC = 0.0D0

C	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
C	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
C	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
C	ECR , ECS, ECT !ECCENTRICITY
	READ(ITI,*) MLE,WW(1),WW(2),WW(3),ECR,ECS,ECT,NFOMA,VFOMA(1:3),GRAV,LOC,LOE,IMOM,ILCN,ILCC 


	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 201
	ENDIF
	ENDDO
	EXIT
201	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)


	CALL XFSECTION(KEG,MLE,1)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !
	

C	-----------------------------------------------------
	DO I = 1,3
	RP   = WW(I)
	RPB  = RP
	DISA = 0.0D0
	DISB = ELN
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
	ENDDO
	ENDDO
	
      ENDDO
C	-----------------------------------------------------
	CALL FRMBEP_OFFSHORE (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) R(IEQ) = R(IEQ) + FRMFOC(J,I)
	ENDDO
      ENDDO
      
	CALL LDASEM_OFFSHORE (R)
      
C     ==================================================      
C     ----- FORCE TO MASS POWER BY TOEY 29-12-2014 -----
C     ==================================================
      IF (NFOMA.EQ.1) THEN
      ! CLEAR FORCE MATRIX
      CALL CLEARA (R,NEQ)
      ! CALL FORCE MASS MATRIX
      CALL FORCEMASSVECTOR (RMA,'READD')
      FRMFOCM = 0.0D0
      DO IMA  = 1,3
      DO IMN  = 1,3
      
      RPM        = WW(IMA)*VFOMA(IMN)
      RPBM       = RPM
      DISA       = 0.0D0
      DISB       = ELN
      IDIR       = IMN
	IND        = 1
      IPIN(1:14) = 0
      IHET       = IHSET(MLE)
          
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
                                           
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF_OFFSHORE(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCC,LOC,LOE,IMOM,'CONT')
	CALL FRAMFIX(RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IDIR,0,VR,
	1			   ELN,FIXFM,RANG,LOC,LOE,IMOM,IPIN)
	DO J = 1,7
	DO K = 1,2
	FRMFOCM(J,K) = FRMFOCM(J,K) + FIXFM(J,K)
	ENDDO
	ENDDO
      
      ENDDO
      ENDDO
C	-----------------------------------------------------
	CALL FRMBEP_OFFSHORE (FRMFOCM,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOCM,NSF,NNF,5)
	ENDIF

	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) R(IEQ) = R(IEQ) + FRMFOCM(J,I)
	ENDDO
      ENDDO
      
      ! UPDATE FORCE MASS MATRIX
      NT_MASS_INDEX = 1.
      DO NT_MASS = 1,NEQ/6D0
      RMA(NT_MASS_INDEX) = RMA(NT_MASS_INDEX) + (RMA(NT_MASS_INDEX)/ABS(GRAV))
      NT_MASS_INDEX = NT_MASS_INDEX + 6D0
      ENDDO
      NT_MASS_INDEX = 2.
      DO NT_MASS = 1,NEQ/6D0
      RMA(NT_MASS_INDEX) = RMA(NT_MASS_INDEX) + (RMA(NT_MASS_INDEX)/ABS(GRAV))
      NT_MASS_INDEX = NT_MASS_INDEX + 6D0
      ENDDO
      NT_MASS_INDEX = 3.
      DO NT_MASS = 1,NEQ/6D0
      RMA(NT_MASS_INDEX) = RMA(NT_MASS_INDEX) + (RMA(NT_MASS_INDEX)/ABS(GRAV))
      NT_MASS_INDEX = NT_MASS_INDEX + 6D0
      ENDDO
      CALL FORCEMASSVECTOR (RMA,'WRITE')
      
      ENDIF
      
C     ==================================================  
      
111   CONTINUE
      RETURN

      END