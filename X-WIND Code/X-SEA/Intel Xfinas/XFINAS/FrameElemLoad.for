C	=======================================================================
C	=======================FRAME THERMAL LOAD==============================
C	=======================================================================
      SUBROUTINE FRAMFOC (LM,XYZ,IGSET,PROPG,MTSET,PROPM,IGIDM,
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
      
C     =================================================================== sorn vlue for offset
      COMMON /sornOffset/ NELO,Nmemoffset(1000),NmemOffsetSelect(1000) ,Nmemtrueelement(1000),GloLoSw(1000) !total number of off set node
C     ================================================================== 
       
C	==================================================================
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019 (REDUCE THE SIZE OF R VECTOR)
      ALLOCATABLE MAL(:) !SONGSAK TWEAK SPEED OCT2019 (REDUCE THE SIZE OF R VECTOR)
      
      DIMENSION R(NEQ)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION XYZ(NCO*NNM,NELE),LM(NEF,NELE)
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
      
C     SORN OFFSET
      DIMENSION SXYZ(6), SVR(3), AAVR(3), BBVR(3) !sorn for calculate new XYZ
    
      
	ALLOCATABLE GPL(:),GPW(:),BPG(:),BWG(:)
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:)
      
	PI = 3.141592654
      WAVEREACFIX = 0.0D0
      
      NWELEMENT = NELE
      
	NSW     = MLOAD(1)	
	NUF     = MLOAD(2)
	NDP     = MLOAD(3)
	NPJ     = MLOAD(4)
	NLD     = MLOAD(5)
	NPN     = MLOAD(6)
	NBD     = MLOAD(7)
	NOFFL   = MLOAD(9)
	NRES    = MLOAD(10)
      NFAT    = MLOAD(11)
	NBOFL   = MLOAD(12)
      NHIROI  = MLOAD(13)

      
C	======================================================================
C	SELF WEIGHT LOADFRAMFOC
	IF(NSW.NE.0) READ (ITI,*)
      
      ALLOCATE(RAL(NEF),MAL(NEF))    
      
	DO 110 ISW = 1,NSW
          
      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0
      
	FRMFOC = 0.0D0

	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR SELFWEIGHT LOAD)
	READ(ITI,*) MLE,GX,GY,GZ,ILCN,ILCC

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
	DO I = 1,3
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	DISA = 0.0D0
	DISB = ELN
	IDIR = I
	IND  = 1
      
      
      DO isorn = 1,NELO
           IF (MLE == Nmemoffset(isorn) ) THEN        
C          calculate length 
C          ===========================select type of offset   
               
           OFtype = NmemOffsetSelect(MLE)
                                       
           SXYZ(1) = XYZ(1,MLE) + A(LOP+(isorn-1)*6)                   
           SXYZ(2) = XYZ(2,MLE) + A(LOP+(isorn-1)*6+1)                  
           SXYZ(3) = XYZ(3,MLE) + A(LOP+(isorn-1)*6+2)                  
           SXYZ(4) = XYZ(4,MLE) + A(LOP+(isorn-1)*6+3)
           SXYZ(5) = XYZ(5,MLE) + A(LOP+(isorn-1)*6+4)
           SXYZ(6) = XYZ(6,MLE) + A(LOP+(isorn-1)*6+5)
 
C         ==========================Calculate AAL BBL
           
           AAVR(1) = SXYZ(1) - XYZ(1,MLE)
           AAVR(2) = SXYZ(2) - XYZ(2,MLE)
           AAVR(3) = SXYZ(3) - XYZ(3,MLE)
           CALL SCALEN(AAVR,AAVR,AAELN,3)
           
           BBVR(1) = XYZ(4,MLE) - SXYZ(4)
           BBVR(2) = XYZ(5,MLE) - SXYZ(5)
           BBVR(3) = XYZ(6,MLE) - SXYZ(6)
           CALL SCALEN(BBVR,BBVR,BBELN,3)
           
           sorn = 3
           
           AAL = (AAELN)/(ELN)
           BBL = (ELN-BBELN)/(ELN)          
          
C          change propertie of uniform load 
      
              RP   = WW(2*I-1)
              RPB  = WW(2*I-0)
              DISA = AAL*ELN
              DISB = BBL*ELN
              IDIR = I
              IND  = 1  
         ENDIF    
      ENDDO

      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)

C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)
      
      
	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
      ENDIF
      
      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + FRMFOC(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)

110	CONTINUE

      DEALLOCATE(RAL,MAL)
      
C	======================================================================
C	UNIFORM LOAD
	IF(NUF.NE.0) READ (ITI,*)
      
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 111 IUF = 1,NUF

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	 
      
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
      
      DO isorn = 1,NELO

           IF (MLE == Nmemoffset(isorn) ) THEN
                             
C          calculate length 
C          ===========================select type of offset        
           OFtype = NmemOffsetSelect(MLE)
                                       
           SXYZ(1) = XYZ(1,MLE) + A(LOP+(isorn-1)*6)                   
           SXYZ(2) = XYZ(2,MLE) + A(LOP+(isorn-1)*6+1)                  
           SXYZ(3) = XYZ(3,MLE) + A(LOP+(isorn-1)*6+2)                  
           SXYZ(4) = XYZ(4,MLE) + A(LOP+(isorn-1)*6+3)
           SXYZ(5) = XYZ(5,MLE) + A(LOP+(isorn-1)*6+4)
           SXYZ(6) = XYZ(6,MLE) + A(LOP+(isorn-1)*6+5)

C         ==========================Calculate AAL BBL
           
           AAVR(1) = SXYZ(1) - XYZ(1,MLE)
           AAVR(2) = SXYZ(2) - XYZ(2,MLE)
           AAVR(3) = SXYZ(3) - XYZ(3,MLE)
           CALL SCALEN(AAVR,AAVR,AAELN,3)
           
           BBVR(1) = XYZ(4,MLE) - SXYZ(4)
           BBVR(2) = XYZ(5,MLE) - SXYZ(5)
           BBVR(3) = XYZ(6,MLE) - SXYZ(6)
           CALL SCALEN(BBVR,BBVR,BBELN,3)
           
           AAL = (AAELN)/(ELN)
           BBL = (ELN-BBELN)/(ELN)          

C          change propertie of uniform load 

              RP   = WW(I)
              RPB  = WW(I)
              DISA = AAL*ELN
              DISB = BBL*ELN
              IDIR = I
              IND  = 1
         ENDIF
      ENDDO
       
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + FRMFOC(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
        
      
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
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOCM,VR,KEG,MLE,RANG)

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
      
      DEALLOCATE(RAL,MAL)   
C	======================================================================
C	======================================================================
C	DISTRIBUTED PARTIAL GLOBAL LOAD
	IF(NDP.NE.0) READ (ITI,*)
      
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 112 IDP = 1,NDP

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	 
      
	FRMFOC = 0.0D0


C	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
C	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
C	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
C	ECR , ECS, ECT !ECCENTRICITY
	READ(ITI,*) MLE,WW(1),WW(2),WW(3),WW(4),WW(5),WW(6),
	1			AAL,BBL,ECR,ECS,ECT,NFOMA,VFOMA(1:3),GRAV,LOC,LOE,IMOM,ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 202
	ENDIF
	ENDDO
	EXIT
202	CONTINUE

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
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	DISA = AAL*ELN
	DISB = BBL*ELN
	IDIR = I
	IND  = 1
      
      DO isorn = 1,NELO
          
           IF (MLE == Nmemoffset(isorn) ) THEN
 
C          calculate length 
C          ===========================select type of offset        
           OFtype = NmemOffsetSelect(MLE)
                                       
           SXYZ(1) = XYZ(1,MLE) + A(LOP+(isorn-1)*6)                   
           SXYZ(2) = XYZ(2,MLE) + A(LOP+(isorn-1)*6+1)                  
           SXYZ(3) = XYZ(3,MLE) + A(LOP+(isorn-1)*6+2)                  
           SXYZ(4) = XYZ(4,MLE) + A(LOP+(isorn-1)*6+3)
           SXYZ(5) = XYZ(5,MLE) + A(LOP+(isorn-1)*6+4)
           SXYZ(6) = XYZ(6,MLE) + A(LOP+(isorn-1)*6+5)
 
C         ==========================Calculate AAL BBL
           
           AAVR(1) = SXYZ(1) - XYZ(1,MLE)
           AAVR(2) = SXYZ(2) - XYZ(2,MLE)
           AAVR(3) = SXYZ(3) - XYZ(3,MLE)
           CALL SCALEN(AAVR,AAVR,AAELN,3)
           
           BBVR(1) = XYZ(4,MLE) - SXYZ(4)
           BBVR(2) = XYZ(5,MLE) - SXYZ(5)
           BBVR(3) = XYZ(6,MLE) - SXYZ(6)
           CALL SCALEN(BBVR,BBVR,BBELN,3)
           
           AALof = (AAELN)/(ELN)
           BBLof = (ELN-BBELN)/(ELN)          
           
C         compare 

C           reset
             
           WWA = 0
           WWB = 0 
           SMaxAAL = 0
           SMinBBL = 0
           
C         find max AAL
           if (AALof >= AAL) then
               SMaxAAL = AALof
           else if (ALLof < AAL) then
               SMaxAAL = AAL     
           endif
           
C           find min BBL
           if (BBLof >= BBL) then
               SMinBBL = BBL
           else if (ALLof < AAL) then
               SMinBBL = BBLof     
           endif

C         interpulated force

           SDIFHig = WW(2*I-0) - WW(2*I-1)
           SDIFDis = (BBL-AAL)*ELN 
           SFacHD   = SDIFHig / SDIFDis
           
           WWA = WW(2*I-1) + (SMaxAAL-AAL)*ELN*SFacHD 
           WWB = WW(2*I-0) - (BBL-SMinBBL)*ELN*SFacHD 
           
           
C          change propertie of uniform load 
           
           	
              RP   = WWA
              RPB  = WWB
              DISA = SMaxAAL*ELN
              DISB = SMinBBL*ELN
              IDIR = I
              IND  = 1
         ENDIF
      ENDDO
      
	

      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + FRMFOC(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
	   
C     ==================================================      
C     ----- FORCE TO MASS POWER BY TOEY 29-12-2014 -----
C     ==================================================
      IF (NFOMA.EQ.1) THEN
      CALL CLEARA (R,NEQ)
      CALL FORCEMASSVECTOR (RMA,'READD')
      FRMFOCM = 0.0D0
      
	DO I  = 1,3
      DO JK = 1,3
	RPM   = WW(2*I-1)*VFOMA(JK)
	RPBM  = WW(2*I-0)*VFOMA(JK)
	DISA  = AAL*ELN
	DISB  = BBL*ELN
	IDIR  = JK
	IND   = 1
	

      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOCM,VR,KEG,MLE,RANG)

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
C     ======================================================================  
112   CONTINUE
      
      DEALLOCATE(RAL,MAL)   
C	======================================================================
      
C	DISTRIBUTED PARTIAL GLOBAL PROJECT LOAD
	IF(NPJ.NE.0) READ (ITI,*) 
      
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 113 IPJ = 1,NPJ 

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	  
      
	FRMFOC = 0.0D0 


C	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
C	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
C	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
C	ECR , ECS, ECT !ECCENTRICITY
	READ(ITI,*) MLE,WW(1),WW(2),WW(3),WW(4),WW(5),WW(6), 
	1			AAL,BBL,ECR,ECS,ECT,NFOMA,VFOMA(1:3),GRAV,LOC,LOE,IMOM,ILCN,ILCC 

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 203
	ENDIF
	ENDDO
	EXIT
203	CONTINUE

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
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	DISA = AAL*ELN
	DISB = BBL*ELN
	IDIR = I
	IND = 2
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCC,LOC,LOE,IMOM,'CONT')
	CALL FRAMFIX(RP,RPB,DISA,DISB,ECR,ECS,ECT,IDIR,1,VR,
	1			   ELN,FIXF,RANG,LOC,LOE,IMOM,IPIN)
	DO J = 1,7
	DO K = 1,2
	FRMFOC(J,K) = FRMFOC(J,K) + FIXF(J,K)
	ENDDO
	ENDDO
	
	ENDDO
C	-----------------------------------------------------
	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF
	  
      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + FRMFOC(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      
C     ==================================================      
C     ----- FORCE TO MASS POWER BY TOEY 29-12-2014 -----
C     ==================================================
      IF (NFOMA.EQ.1) THEN
      CALL CLEARA (R,NEQ)
      CALL FORCEMASSVECTOR (RMA,'READD')
	DO I  = 1,3
      DO JK = 1,3
	RPM   = WW(2*I-1)*VFOMA(JK)
	RPBM  = WW(2*I-0)*VFOMA(JK)
	DISA  = AAL*ELN
	DISB  = BBL*ELN
	IDIR  = JK
	IND   = 1
	

      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOCM,VR,KEG,MLE,RANG)

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

C     ======================================================================  
113   CONTINUE
      
      DEALLOCATE(RAL,MAL)
C	======================================================================
      
C	LINEAR VARYING GLOBAL LOAD
	IF(NLD.NE.0) READ (ITI,*)
      
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 114 ILD = 1,NLD

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	  
      
	FRMFOC = 0.0D0


C	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
C	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
C	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
C	ECR , ECS, ECT !ECCENTRICITY
C	IDIR  coord  TX  TY  TZ  ECC  ILCN
	READ(ITI,*) MLE,LDIR,CC,TX,TY,TZ,ECR,ECS,ECT,NFOMA,VFOMA(1:3),GRAV,LOC,LOE,ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 204
	ENDIF
	ENDDO
	EXIT
204	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)

	CALL XFSECTION(KEG,MLE,1)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !
	

C	-----------------------------------------------------
	TEST1 = CC - XYZ(LDIR+0,MLE)
	TEST2 = CC - XYZ(LDIR+3,MLE)
	IF(TEST1.GE.0.0.AND.TEST2.GE.0.0) THEN
	DISA = 0.0D0
	DISB = ELN
	WW(1) = TX*TEST1
	WW(2) = TX*TEST2
	WW(3) = TY*TEST1
	WW(4) = TY*TEST2
	WW(5) = TZ*TEST1
	WW(6) = TZ*TEST2
	ENDIF
	IF(TEST1.LT.0.0.AND.TEST2.GE.0.0) THEN
	X1 = XYZ(LDIR+0,MLE)
	X2 = XYZ(LDIR+3,MLE)
	DISA = ELN*(1.0-(CC-X2)/(X1-X2))
	DISB = ELN
	WW(1) = 0.0D0
	WW(2) = TX*TEST2
	WW(3) = 0.0D0
	WW(4) = TY*TEST2
	WW(5) = 0.0D0
	WW(6) = TZ*TEST2
	ENDIF
	IF(TEST1.GE.0.0.AND.TEST2.LT.0.0) THEN
	X1 = XYZ(LDIR+0,MLE)
	X2 = XYZ(LDIR+3,MLE)
	DISA = 0.0D0
	DISB = ELN*((CC-X1)/(X2-X1))
	WW(1) = TX*TEST1
	WW(2) = 0.0D0
	WW(3) = TY*TEST1
	WW(4) = 0.0D0
	WW(5) = TZ*TEST1
	WW(6) = 0.0D0
	ENDIF
	IF(TEST1.LT.0.0.AND.TEST2.LT.0.0) THEN
	DISA = 0.0D0
	DISB = ELN
	WW(1) = 0.0D0
	WW(2) = 0.0D0
	WW(3) = 0.0D0
	WW(4) = 0.0D0
	WW(5) = 0.0D0
	WW(6) = 0.0D0
	ENDIF
C	-----------------------------------------------------
	DO I = 1,3
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	IDIR = I
	IND = 1
      
      
      DO isorn = 1,NELO
          
           IF (MLE == Nmemoffset(isorn) ) THEN
 
C          calculate length 
C          ===========================select type of offset        
           OFtype = NmemOffsetSelect(MLE)
                                       
           SXYZ(1) = XYZ(1,MLE) + A(LOP+(isorn-1)*6)                   
           SXYZ(2) = XYZ(2,MLE) + A(LOP+(isorn-1)*6+1)                  
           SXYZ(3) = XYZ(3,MLE) + A(LOP+(isorn-1)*6+2)                  
           SXYZ(4) = XYZ(4,MLE) + A(LOP+(isorn-1)*6+3)
           SXYZ(5) = XYZ(5,MLE) + A(LOP+(isorn-1)*6+4)
           SXYZ(6) = XYZ(6,MLE) + A(LOP+(isorn-1)*6+5)
 
C         ==========================Calculate AAL BBL
           
           AAVR(1) = SXYZ(1) - XYZ(1,MLE)
           AAVR(2) = SXYZ(2) - XYZ(2,MLE)
           AAVR(3) = SXYZ(3) - XYZ(3,MLE)
           CALL SCALEN(AAVR,AAVR,AAELN,3)
           
           BBVR(1) = XYZ(4,MLE) - SXYZ(4)
           BBVR(2) = XYZ(5,MLE) - SXYZ(5)
           BBVR(3) = XYZ(6,MLE) - SXYZ(6)
           CALL SCALEN(BBVR,BBVR,BBELN,3)
           
           AALof = (AAELN)/(ELN)
           BBLof = (ELN-BBELN)/(ELN)          
           
C         compare 

C           reset
             
           WWA = 0
           WWB = 0 
           SMaxAAL = 0
           SMinBBL = 0
           
C         find max AAL
           if (AALof >= AAL) then
               SMaxAAL = AALof
           else if (ALLof < AAL) then
               SMaxAAL = AAL     
           endif
           
C           find min BBL
           if (BBLof >= BBL) then
               SMinBBL = BBL
           else if (ALLof < AAL) then
               SMinBBL = BBLof     
           endif

C         interpulated force

           SDIFHig = WW(2*I-0) - WW(2*I-1)
           SDIFDis = (BBL-AAL)*ELN 
           SFacHD   = SDIFHig / SDIFDis
           
           WWA = WW(2*I-1) + (SMaxAAL-AAL)*ELN*SFacHD 
           WWB = WW(2*I-0) - (BBL-SMinBBL)*ELN*SFacHD 
           
           
C          change propertie of uniform load 
           
           	
              RP   = WWA
              RPB  = WWB
              DISA = SMaxAAL*ELN
              DISB = SMinBBL*ELN
              IDIR = I
              IND  = 1
         ENDIF
      ENDDO
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + FRMFOC(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
	  
C     ==================================================      
C     ----- FORCE TO MASS POWER BY TOEY 29-12-2014 -----
C     ==================================================
      IF (NFOMA.EQ.1) THEN
      CALL CLEARA (R,NEQ)
      CALL FORCEMASSVECTOR (RMA,'READD')
	DO I  = 1,3
      DO JK = 1,3
	RPM   = WW(2*I-1)*VFOMA(JK)
	RPBM  = WW(2*I-0)*VFOMA(JK)
	IDIR  = JK
	IND   = 1
	

      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOCM,VR,KEG,MLE,RANG)

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


114   CONTINUE
      
      DEALLOCATE(RAL,MAL)
C	======================================================================
      
C	GLOBAL POINT LOAD
	IF(NPN.NE.0) READ (ITI,*)
      
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO 115 IPN = 1,NPN

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	  
	  
	FRMFOC = 0.0D0

C	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
C	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
C	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
C	ECR , ECS, ECT !ECCENTRICITY
	READ(ITI,*) MLE,WW(1),WW(2),WW(3),AAL,ECR,ECS,ECT,NFOMA,VFOMA(1:3),GRAV,LOC,LOE,IMOM,ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 205
	ENDIF
	ENDDO
	EXIT
205	CONTINUE

	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)

	CALL XFSECTION(KEG,MLE,1)
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0) !

      !WW(1) = 0.0d0
      
      
C	-----------------------------------------------------

	DO I = 1,3
	RP   = WW(I)/1.0E-8
	RPB  = RP
	DISA = AAL*ELN
	DISB = DISA+1.0E-8
	IDIR = I
	IND = 1
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
      
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + FRMFOC(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
	  
C     ==================================================      
C     ----- FORCE TO MASS POWER BY TOEY 29-12-2014 -----
C     ==================================================
      IF (NFOMA.EQ.1) THEN
      CALL CLEARA (R,NEQ)
      CALL FORCEMASSVECTOR (RMA,'READD')
	DO I  = 1,3
      DO JK = 1,3
	RPM   = (WW(I)/1.0E-8)*VFOMA(JK)
	RPBM  = RPM
	DISA  = AAL*ELN
	DISB  = DISA+1.0E-8
	IDIR  = JK
	IND   = 1

      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RPM,RPBM,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOCM,VR,KEG,MLE,RANG)

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

C     ======================================================================  
115   CONTINUE
      
      DEALLOCATE(RAL,MAL)
C	======================================================================
      
C	BODY LOAD
	IF(NBD.NE.0) READ (ITI,*)
      
      ALLOCATE(RAL(NEF),MAL(NEF))

	DO 116 IBD = 1,NBD

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	  
	  
	FRMFOC = 0.0D0


C	LOC = LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
	ECR = 0.0D0 ; ECS = 0.0D0 ; ECT = 0.0D0 !ECCENTRICITY (NOTHING FOR BODY LOAD)
	READ(ITI,*) MLE,GX,GY,GZ,LOC,ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MLE) THEN
	MLE = IGID 
	GOTO 206
	ENDIF
	ENDDO
	EXIT
206	CONTINUE

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
	DO I = 1,3
	RP   = WW(2*I-1)
	RPB  = WW(2*I-0)
	DISA = 0.0D0
	DISB = ELN
	IDIR = I
	IND = 1
	
      IPIN(1:14) = 0
      IHET = IHSET(MLE)
      IF(IHET.NE.0) IPIN(1:14) = LRPIN(1:14,IHET)
      
C	DOING SMOOTH LINE DIAGRAM CALCULATION
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
	1			 XYZ,RANG,NP_SMH,ILCN,LOC,LOE,IMOM,'VARY')
	CALL SMHUNIF(MLE,IDIR,RP,RPB,DISA,DISB,ECR,ECS,ECT,IND,
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
	CALL FRMBEP (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
	ENDIF

      IEFL = 0
	IK = 0
	DO I = 1,2
	DO J = 1,7
	IK  = IK + 1
      IEQ = LM(IK,MLE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + FRMFOC(J,I)
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	ENDDO

	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
	  
116   CONTINUE
      
      DEALLOCATE(RAL,MAL)
C	======================================================================
	RETURN

      END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE FRAMFIX(W1,W2,AL,BL,ECR,ECS,ECT,IDR,IPRO,VR,
	1				     ELN,RG,ANG,LOC,LOE,IMOM,LREAS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	==============================================================
	DIMENSION VR(3),VS(3),VT(3),VL(3),RG(7,2)
	DIMENSION COEF(6),TRANS(14,14)
	DIMENSION FIXG(14),FIXD(14)
	DIMENSION ECCMT(3,3)
	DIMENSION VPRO(3)
	DIMENSION TRANH(14,14),LREAS(14)
	DIMENSION TRANE(3,3)

C	FOR LOCAL LOAD
	VL(1:3) = 0.0
	IF(LOC.NE.0) VL(IDR) = 1.0

	RANG = ANG

C	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
C	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
C	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
C	ECR , ECS, ECT !ECCENTRICITY

	ECCMT = 0.0
	ECCMT(1,1) =  0.0D0
	ECCMT(1,2) = -ECT
	ECCMT(1,3) =  ECS
	ECCMT(2,1) =  ECT
	ECCMT(2,2) =  0.0D0
	ECCMT(2,3) = -ECR
	ECCMT(3,1) = -ECS
	ECCMT(3,2) =  ECR
	ECCMT(3,3) =  0.0D0

	DO I = 1,3
	VPRO(I) = 0.0
	ENDDO

	VPRO(IDR) = 1.0

	COST = VR(1)*VPRO(1) +  VR(2)*VPRO(2) + VR(3)*VPRO(3)
	COST = SQRT(1.0 - COST*COST)

	IF(IPRO.EQ.1) THEN
	W1 = W1*COST
	W2 = W2*COST
	ENDIF

	IF(AL.EQ.BL) THEN
C	DX = ELN*1.0E-8
C	BL = BL+DX
C	W2 = W1/DX
C	W1 = W1/DX
	W2 = W1	
	ENDIF

	CALL FMVEVR(VR,VS,VT)

	CALL ROMBAC(VR,VS,VT,RANG)

	CALL TRANLG(VR,VS,VT,TRANS)

C     IF ECC IS IN GLOBAL THEN DO THE TRANSFORMATION
	IF(LOE.EQ.0) THEN !TRANSFORM THE ECC FROM GLOBAL TO LOCAL (E-XYZ TO E-RST)
	    TRANE(1:3,1) = VR(1:3)
	    TRANE(1:3,2) = VS(1:3)
	    TRANE(1:3,3) = VT(1:3)
	    ECCMT = MATMUL(TRANSPOSE(TRANE),MATMUL(ECCMT,TRANE))
	ENDIF
	
	IF(LOC.EQ.0) THEN !GLOBAL
	W1R = VR(IDR)*W1
	W1S = VS(IDR)*W1
	W1T = VT(IDR)*W1
	ELSE !LOCAL
	W1R = VL(1)*W1
	W1S = VL(2)*W1
	W1T = VL(3)*W1
	ENDIF

	W1MR = ECCMT(1,1)*W1R + ECCMT(1,2)*W1S + ECCMT(1,3)*W1T
	W1MS = ECCMT(2,1)*W1R + ECCMT(2,2)*W1S + ECCMT(2,3)*W1T 
	W1MT = ECCMT(3,1)*W1R + ECCMT(3,2)*W1S + ECCMT(3,3)*W1T 

	IF(LOC.EQ.0) THEN !GLOBAL
	W2R = VR(IDR)*W2
	W2S = VS(IDR)*W2
	W2T = VT(IDR)*W2
	ELSE !LOCAL
	W2R = VL(1)*W2
	W2S = VL(2)*W2
	W2T = VL(3)*W2
	ENDIF

	W2MR = ECCMT(1,1)*W2R + ECCMT(1,2)*W2S + ECCMT(1,3)*W2T
	W2MS = ECCMT(2,1)*W2R + ECCMT(2,2)*W2S + ECCMT(2,3)*W2T 
	W2MT = ECCMT(3,1)*W2R + ECCMT(3,2)*W2S + ECCMT(3,3)*W2T 

	FIXD = 0.0

	IF(IMOM.EQ.1) THEN
	W1MR = W1R
	W1MS = W1S
	W1MT = W1T

	W2MR = W2R
	W2MS = W2S
	W2MT = W2T
	GOTO 100
	ENDIF


C	FROM CONCENTRIC LOAD

C	LOCAL AXIAL FORCE
	CALL FXCONT(W1R,W2R,AL,BL,ELN,COEF)
	FIXD(1) = COEF(5)
	FIXD(8) = COEF(6)


C	SHEAR IN S-AXIS AND MOMENT IN T-AXIS
	CALL FXCONT(W1S,W2S,AL,BL,ELN,COEF)
	FIXD(2)  = COEF(1)
	FIXD(6)  = COEF(2)
	FIXD(9)  = COEF(3)
	FIXD(13) = COEF(4)	

C	SHEAR IN T-AXIS AND MOMENT IN S-AXIS
	CALL FXCONT(W1T,W2T,AL,BL,ELN,COEF)
	FIXD(3)  = COEF(1)
	FIXD(5)  =-COEF(2)
	FIXD(10) = COEF(3)
	FIXD(12) =-COEF(4)	


100	CONTINUE
C	FROM ECCENTRICITY MOMENT

C	LOCAL AXIAL MOMENT
	CALL FMCONT(W1MR,W2MR,AL,BL,ELN,COEF)
	FIXD(4)  = FIXD(4)  + COEF(5)
	FIXD(11) = FIXD(11) + COEF(6)


C	SHEAR IN T-AXIS AND MOMENT IN S-AXIS
	CALL FMCONT(W1MS,W2MS,AL,BL,ELN,COEF)
	FIXD(3)  = FIXD(3)  - COEF(1)
	FIXD(5)  = FIXD(5)  + COEF(2)
	FIXD(10) = FIXD(10) - COEF(3)
	FIXD(12) = FIXD(12) + COEF(4)


C	SHEAR IN S-AXIS AND MOMENT IN T-AXIS
	CALL FMCONT(W1MT,W2MT,AL,BL,ELN,COEF)
	FIXD(2)  = FIXD(2)  + COEF(1)
	FIXD(6)  = FIXD(6)  + COEF(2)
	FIXD(9)  = FIXD(9)  + COEF(3)
	FIXD(13) = FIXD(13) + COEF(4)

C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
	CALL TRNHIG(TRANH,ELN,LREAS)
	CALL TRNMUL(TRANH,FIXD,2)


	FIXG = MATMUL(TRANS,FIXD)

	DO I = 1,7
	RG(I,1) = 0.0
	RG(I,2) = 0.0
	RG(I,1) = FIXG(I)
	RG(I,2) = FIXG(I+7)
	ENDDO


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
C	=======================================================================
C	START BEAM&FRAME SMOOTH LINE DIAGRAM BY INCLUDE INTERNAL POINT OF 
C	PRINTING STRESS&RESULTANT FORCE SONGSAK FEB2007
C	=======================================================================
	SUBROUTINE SMHZERO
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	=======================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH
      COMMON A(9000000),IA(9000000)
C	=======================================================================
C	-----------------------------------------------------------------------
	DIMENSION RESUL(7)
C	-----------------------------------------------------------------------

	IF(ITYPE.NE.5) RETURN

	CALL INTFILL('OGRF',NFL,13,KEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM

	RESUL(1:7) = 0.0D0
C	--------------------------------------------
C	WRITE TO THE FILE
C	--------------------------------------------
	DO ILC = 1,LCS
	DO IEL = 1,NELE
	DO IP  = 1,NP_SMH
	IRO = IP + NP_SMH*(IEL-1) + NP_SMH*NELE*(ILC-1)
	IRC = 2*IRO-1
	WRITE(NFL,REC=IRC) RESUL(1:7)
	IRC = 2*IRO-0
	WRITE(NFL,REC=IRC) RESUL(1:7)
	ENDDO
	ENDDO
	ENDDO
C	--------------------------------------------
 

	RETURN

      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
C	=======================================================================
C	START BEAM&FRAME SMOOTH LINE DIAGRAM BY INCLUDE INTERNAL POINT OF 
C	PRINTING STRESS&RESULTANT FORCE SONGSAK FEB2007
C	=======================================================================
	SUBROUTINE SMHZERO_OFFSHORE
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	=======================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH
C	=======================================================================
C	-----------------------------------------------------------------------
	DIMENSION RESUL(7)
C	-----------------------------------------------------------------------

	IF(ITYPE.NE.5) RETURN

	CALL INTFILL('OGRF',NFL,20,KEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM

	RESUL(1:7) = 0.0D0
C	--------------------------------------------
C	WRITE TO THE FILE
C	--------------------------------------------
	DO ILC = 1,LCS
	DO IEL = 1,NELE
	DO IP  = 1,NP_SMH
	IRO = IP + NP_SMH*(IEL-1) + NP_SMH*NELE*(ILC-1)
	IRC = 2*IRO-1
	WRITE(NFL,REC=IRC) RESUL(1:7)
	IRC = 2*IRO-0
	WRITE(NFL,REC=IRC) RESUL(1:7)
	ENDDO
	ENDDO
	ENDDO
C	--------------------------------------------
 

	RETURN

	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SMHUNIF(MLE,IDIR,W1,W2,AL,BL,ECR,ECS,ECT,IND,
	1				     XYZ,RANG,NP,ILC,LOC,LOE,IMOM,TYP) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C	=======================================================================
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	=======================================================================
	DIMENSION VR(3),VS(3),VT(3),XYZ(NCO*NNM,NELE)
	DIMENSION RESUL(7),ECCMT(3,3),VPRO(3),VL(3)
	DIMENSION TRANE(3,3)

	IF(ITYPE.NE.5) RETURN
	IF(ILC.LE.0)   RETURN

	CALL INTFILL('OGRF',NFL,13,KEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM


C	FOR LOCAL LOAD
	VL(1:3) = 0.0
	IF(LOC.NE.0) VL(IDIR) = 1.0

C	-----------------------------------------------------------
C	DETERMINE LOCAL BASE VECTOR
C	-----------------------------------------------------------
	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)

	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)
C	-----------------------------------------------------------
C	ECCENTRICITY EFFECT
C	-----------------------------------------------------------
	ECCMT = 0.0
	ECCMT(1,1) =  0.0D0
	ECCMT(1,2) = -ECT
	ECCMT(1,3) =  ECS
	ECCMT(2,1) =  ECT
	ECCMT(2,2) =  0.0D0
	ECCMT(2,3) = -ECR
	ECCMT(3,1) = -ECS
	ECCMT(3,2) =  ECR
	ECCMT(3,3) =  0.0D0
C     IF ECC IS IN GLOBAL THEN DO THE TRANSFORMATION
	IF(LOE.EQ.0) THEN !TRANSFORM THE ECC FROM GLOBAL TO LOCAL (E-XYZ TO E-RST)
	    TRANE(1:3,1) = VR(1:3)
	    TRANE(1:3,2) = VS(1:3)
	    TRANE(1:3,3) = VT(1:3)
	    ECCMT = MATMUL(TRANSPOSE(TRANE),MATMUL(ECCMT,TRANE))
      ENDIF
C	-----------------------------------------------------------
C	MODIFY PROJECTED LOAD
C	-----------------------------------------------------------
	VPRO(1:3)  = 0.0
	VPRO(IDIR) = 1.0D0

	COST = VR(1)*VPRO(1) +  VR(2)*VPRO(2) + VR(3)*VPRO(3)
	COST = SQRT(1.0 - COST*COST)
	IF(IND.EQ.2) THEN
	W1 = W1*COST
	W2 = W2*COST
	ENDIF

C	-----------------------------------------------------------
C	MODIFY FULL UNIFORM LOAD
C	-----------------------------------------------------------
	IF(IND.EQ.4) THEN
	W2 = W1
	AL = 0.0
	BL = ELN
	ENDIF
C	-----------------------------------------------------------
C	MODIFY SELFWEIGHT LOAD
C	-----------------------------------------------------------
	IF(IND.EQ.3) THEN
	AL = 0.0
	BL = ELN
	ENDIF
C	-----------------------------------------------------------
C	LOAD VALUE IN LOCAL SYSTEM
C	-----------------------------------------------------------
	IF(LOC.EQ.0) THEN
	W1R = VR(IDIR)*W1 !FORCE
	W1S = VS(IDIR)*W1
	W1T = VT(IDIR)*W1
	ELSE
	W1R = VL(1)*W1 !FORCE
	W1S = VL(2)*W1
	W1T = VL(3)*W1
	ENDIF

	W1MR = ECCMT(1,1)*W1R + ECCMT(1,2)*W1S + ECCMT(1,3)*W1T !MOMENT
	W1MS = ECCMT(2,1)*W1R + ECCMT(2,2)*W1S + ECCMT(2,3)*W1T 
	W1MT = ECCMT(3,1)*W1R + ECCMT(3,2)*W1S + ECCMT(3,3)*W1T 

	IF(LOC.EQ.0) THEN
	W2R = VR(IDIR)*W2 !FORCE
	W2S = VS(IDIR)*W2
	W2T = VT(IDIR)*W2
	ELSE
	W2R = VL(1)*W2 !FORCE
	W2S = VL(2)*W2
	W2T = VL(3)*W2
	ENDIF

	W2MR = ECCMT(1,1)*W2R + ECCMT(1,2)*W2S + ECCMT(1,3)*W2T !MOMENT
	W2MS = ECCMT(2,1)*W2R + ECCMT(2,2)*W2S + ECCMT(2,3)*W2T 
	W2MT = ECCMT(3,1)*W2R + ECCMT(3,2)*W2S + ECCMT(3,3)*W2T 
C	-----------------------------------------------------------

	IF(IMOM.EQ.1) THEN
	W1MR = W1R
	W1MS = W1S
	W1MT = W1T
	W1R = 0.0D0
	W1S = 0.0D0
	W1T = 0.0D0

	W2MR = W2R
	W2MS = W2S
	W2MT = W2T
	W2R = 0.0D0
	W2S = 0.0D0
	W2T = 0.0D0
	ENDIF
	
	
	DL = BL - AL

	X   = 0.0
	DS  = ELN/(NP+1)
	DO 1000 IP = 1,NP

	IRC = 0
	IRO = IP + NP*(MLE-1) + NP*NELE*(ILC-1)
	IF(TYP.EQ.'VARY') IRC = 2*IRO-1	      !VARY LOAD LINE NUMBER 2*N-1
	IF(TYP.EQ.'CONT') IRC = 2*IRO-0     !CONT LOAD LINE NUMBER 2*N
	IF(IRC.EQ.0) RETURN
	READ(NFL,REC=IRC) RESUL(1:7)

	X = X + DS

	XM1 = VMCAULY(X,AL)
	XM2 = VMCAULY(X,BL)

	IF(AL.EQ.BL) GOTO 100 !FOR POINT LOAD


C	FROM WMR
	CALL TRAPIZ(XM1,DL,W1MR,W2MR,W1MR,DMR1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MR,W2MR,W2MR,DMR2,DUM1)
	RESUL(4) = RESUL(4) + DMR1-DMR2  !TORSION

C	FROM WMS
	CALL TRAPIZ(DL ,DL,W1MS,W2MS,W1MS,D1,DUM1)
	RT = -D1/ELN
	CALL TRAPIZ(XM1,DL,W1MS,W2MS,W1MS,D1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MS,W2MS,W2MS,D2,DUM1)
	RMS = RT*X + D1 - D2 
	RESUL(3) = RESUL(3)       !SHEAR T
	RESUL(5) = RESUL(5) + RMS !BENDING S

C	FROM WMT
	CALL TRAPIZ(DL ,DL,W1MT,W2MT,W1MT,D1,DUM1)
	RS = -D1/ELN
	CALL TRAPIZ(XM1,DL,W1MT,W2MT,W1MT,D1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MT,W2MT,W2MT,D2,DUM1)
	RMT = RS*X + D1 - D2
	RESUL(2) = RESUL(2)       !SHEAR S
	RESUL(6) = RESUL(6) - RMT !BENDING T



C	FROM WR
	CALL TRAPIZ(XM1,DL,W1R,W2R,W1R,RF1,DUM1)
	CALL TRAPIZ(XM2,DL,W1R,W2R,W2R,RF2,DUM1)
	RESUL(1) = RESUL(1) + RF1 - RF2  !AXIAL

C	FROM WS
	CALL TRAPIZ(DL ,DL,W1S,W2S,W1S,D1,DIS)
	D2 = D1*(ELN-AL-DIS)
	RS =-D2/ELN
	CALL TRAPIZ(XM1,DL,W1S,W2S,W1S,D1,DIS1)
	CALL TRAPIZ(XM2,DL,W1S,W2S,W2S,D2,DIS2)
	D3 = D1*(X  -AL-DIS1)
	D4 = D2*(X  -BL-DIS2)
	RMT = RS*X + D3 - D4
	RESUL(2) = RESUL(2) + D1 - D2		 !SHEAR S
	RESUL(6) = RESUL(6) + RMT          !BENDING T


C	FROM WT
	CALL TRAPIZ(DL ,DL,W1T,W2T,W1T,D1,DIS)
	D2 = D1*(ELN-AL-DIS)
	RT =-D2/ELN
	CALL TRAPIZ(XM1,DL,W1T,W2T,W1T,D1,DIS1)
	CALL TRAPIZ(XM2,DL,W1T,W2T,W2T,D2,DIS2)
	D3 = D1*(X  -AL-DIS1)
	D4 = D2*(X  -BL-DIS2)
	RMS = RT*X + D3 - D4
	RESUL(3) = RESUL(3) + D1 - D2		 !SHEAR T
	RESUL(5) = RESUL(5) - RMS          !BENDING S

	GOTO 200

100	CONTINUE

	DUM = 0.0
	IF(XM1.GT.0.0) DUM = 1.0 

C	FROM P MR
	RESUL(4) = RESUL(4) + W1MR*DUM  !TORSION

C	FROM P MS
	RT = -W1MS/ELN
	RMS = RT*X + W1MS*DUM
	RESUL(3) = RESUL(3)       !SHEAR T
	RESUL(5) = RESUL(5) + RMS !BENDING S

C	FROM P MT
	RS = -W1MT/ELN
	RMT = RS*X + W1MT*DUM
	RESUL(2) = RESUL(2)       !SHEAR S
	RESUL(6) = RESUL(6) - RMT !BENDING T


C	FROM P R
	RESUL(1) = RESUL(1) + W1R*DUM  !AXIAL

C	FROM P S
	D2 = W1S*(ELN-AL)
	RS =-D2/ELN
	RMT = RS*X   + W1S*XM1
	RESUL(2) = RESUL(2) + W1S*DUM       !SHEAR S
	RESUL(6) = RESUL(6) + RMT           !BENDING T

C	FROM P T
	D2 = W1T*(ELN-AL)
	RT =-D2/ELN
	RMS = RT*X   + W1T*XM1
	RESUL(3) = RESUL(3) + W1T*DUM       !SHEAR T
	RESUL(5) = RESUL(5) - RMS           !BENDING S


200	CONTINUE


	IRC = 0
	IRO = IP + NP*(MLE-1) + NP*NELE*(ILC-1)
	IF(TYP.EQ.'VARY') IRC = 2*IRO-1	      !VARY LOAD LINE NUMBER 2*N-1
	IF(TYP.EQ.'CONT') IRC = 2*IRO-0     !CONT LOAD LINE NUMBER 2*N
	IF(IRC.EQ.0) RETURN
	WRITE(NFL,REC=IRC) RESUL(1:7)
	
1000	CONTINUE
C	-----------------------------------------------------------	

	RETURN
      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SMHUNIF_OFFSHORE(MLE,IDIR,W1,W2,AL,BL,ECR,ECS,ECT,IND,
	1				     XYZ,RANG,NP,ILC,LOC,LOE,IMOM,TYP) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C	=======================================================================
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	=======================================================================
	DIMENSION VR(3),VS(3),VT(3),XYZ(NCO*NNM,NELE)
	DIMENSION RESUL(7),ECCMT(3,3),VPRO(3),VL(3)
	DIMENSION TRANE(3,3)

	IF(ITYPE.NE.5) RETURN
	IF(ILC.LE.0)   RETURN

	CALL INTFILL('OGRF',NFL,20,KEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM


C	FOR LOCAL LOAD
	VL(1:3) = 0.0
	IF(LOC.NE.0) VL(IDIR) = 1.0

C	-----------------------------------------------------------
C	DETERMINE LOCAL BASE VECTOR
C	-----------------------------------------------------------
	VR(1) = XYZ(4,MLE)-XYZ(1,MLE)
	VR(2) = XYZ(5,MLE)-XYZ(2,MLE)
	VR(3) = XYZ(6,MLE)-XYZ(3,MLE)
	CALL SCALEN(VR,VR,ELN,3)

	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)
C	-----------------------------------------------------------
C	ECCENTRICITY EFFECT
C	-----------------------------------------------------------
	ECCMT = 0.0
	ECCMT(1,1) =  0.0D0
	ECCMT(1,2) = -ECT
	ECCMT(1,3) =  ECS
	ECCMT(2,1) =  ECT
	ECCMT(2,2) =  0.0D0
	ECCMT(2,3) = -ECR
	ECCMT(3,1) = -ECS
	ECCMT(3,2) =  ECR
	ECCMT(3,3) =  0.0D0
C     IF ECC IS IN GLOBAL THEN DO THE TRANSFORMATION
	IF(LOE.EQ.0) THEN !TRANSFORM THE ECC FROM GLOBAL TO LOCAL (E-XYZ TO E-RST)
	    TRANE(1:3,1) = VR(1:3)
	    TRANE(1:3,2) = VS(1:3)
	    TRANE(1:3,3) = VT(1:3)
	    ECCMT = MATMUL(TRANSPOSE(TRANE),MATMUL(ECCMT,TRANE))
      ENDIF
C	-----------------------------------------------------------
C	MODIFY PROJECTED LOAD
C	-----------------------------------------------------------
	VPRO(1:3)  = 0.0
	VPRO(IDIR) = 1.0D0

	COST = VR(1)*VPRO(1) +  VR(2)*VPRO(2) + VR(3)*VPRO(3)
	COST = SQRT(1.0 - COST*COST)
	IF(IND.EQ.2) THEN
	W1 = W1*COST
	W2 = W2*COST
	ENDIF

C	-----------------------------------------------------------
C	MODIFY FULL UNIFORM LOAD
C	-----------------------------------------------------------
	IF(IND.EQ.4) THEN
	W2 = W1
	AL = 0.0
	BL = ELN
	ENDIF
C	-----------------------------------------------------------
C	MODIFY SELFWEIGHT LOAD
C	-----------------------------------------------------------
	IF(IND.EQ.3) THEN
	AL = 0.0
	BL = ELN
	ENDIF
C	-----------------------------------------------------------
C	LOAD VALUE IN LOCAL SYSTEM
C	-----------------------------------------------------------
	IF(LOC.EQ.0) THEN
	W1R = VR(IDIR)*W1 !FORCE
	W1S = VS(IDIR)*W1
	W1T = VT(IDIR)*W1
	ELSE
	W1R = VL(1)*W1 !FORCE
	W1S = VL(2)*W1
	W1T = VL(3)*W1
	ENDIF

	W1MR = ECCMT(1,1)*W1R + ECCMT(1,2)*W1S + ECCMT(1,3)*W1T !MOMENT
	W1MS = ECCMT(2,1)*W1R + ECCMT(2,2)*W1S + ECCMT(2,3)*W1T 
	W1MT = ECCMT(3,1)*W1R + ECCMT(3,2)*W1S + ECCMT(3,3)*W1T 

	IF(LOC.EQ.0) THEN
	W2R = VR(IDIR)*W2 !FORCE
	W2S = VS(IDIR)*W2
	W2T = VT(IDIR)*W2
	ELSE
	W2R = VL(1)*W2 !FORCE
	W2S = VL(2)*W2
	W2T = VL(3)*W2
	ENDIF

	W2MR = ECCMT(1,1)*W2R + ECCMT(1,2)*W2S + ECCMT(1,3)*W2T !MOMENT
	W2MS = ECCMT(2,1)*W2R + ECCMT(2,2)*W2S + ECCMT(2,3)*W2T 
	W2MT = ECCMT(3,1)*W2R + ECCMT(3,2)*W2S + ECCMT(3,3)*W2T 
C	-----------------------------------------------------------

	IF(IMOM.EQ.1) THEN
	W1MR = W1R
	W1MS = W1S
	W1MT = W1T
	W1R = 0.0D0
	W1S = 0.0D0
	W1T = 0.0D0

	W2MR = W2R
	W2MS = W2S
	W2MT = W2T
	W2R = 0.0D0
	W2S = 0.0D0
	W2T = 0.0D0
	ENDIF
	
	
	DL = BL - AL

	X   = 0.0
	DS  = ELN/(NP+1)
	DO 1000 IP = 1,NP

	IRC = 0
	IRO = IP + NP*(MLE-1) + NP*NELE*(ILC-1)
	IF(TYP.EQ.'VARY') IRC = 2*IRO-1	      !VARY LOAD LINE NUMBER 2*N-1
	IF(TYP.EQ.'CONT') IRC = 2*IRO-0     !CONT LOAD LINE NUMBER 2*N
	IF(IRC.EQ.0) RETURN
	READ(NFL,REC=IRC) RESUL(1:7)

	X = X + DS

	XM1 = VMCAULY(X,AL)
	XM2 = VMCAULY(X,BL)

	IF(AL.EQ.BL) GOTO 100 !FOR POINT LOAD


C	FROM WMR
	CALL TRAPIZ(XM1,DL,W1MR,W2MR,W1MR,DMR1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MR,W2MR,W2MR,DMR2,DUM1)
	RESUL(4) = RESUL(4) + DMR1-DMR2  !TORSION

C	FROM WMS
	CALL TRAPIZ(DL ,DL,W1MS,W2MS,W1MS,D1,DUM1)
	RT = -D1/ELN
	CALL TRAPIZ(XM1,DL,W1MS,W2MS,W1MS,D1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MS,W2MS,W2MS,D2,DUM1)
	RMS = RT*X + D1 - D2 
	RESUL(3) = RESUL(3)       !SHEAR T
	RESUL(5) = RESUL(5) + RMS !BENDING S

C	FROM WMT
	CALL TRAPIZ(DL ,DL,W1MT,W2MT,W1MT,D1,DUM1)
	RS = -D1/ELN
	CALL TRAPIZ(XM1,DL,W1MT,W2MT,W1MT,D1,DUM1)
	CALL TRAPIZ(XM2,DL,W1MT,W2MT,W2MT,D2,DUM1)
	RMT = RS*X + D1 - D2
	RESUL(2) = RESUL(2)       !SHEAR S
	RESUL(6) = RESUL(6) - RMT !BENDING T



C	FROM WR
	CALL TRAPIZ(XM1,DL,W1R,W2R,W1R,RF1,DUM1)
	CALL TRAPIZ(XM2,DL,W1R,W2R,W2R,RF2,DUM1)
	RESUL(1) = RESUL(1) + RF1 - RF2  !AXIAL

C	FROM WS
	CALL TRAPIZ(DL ,DL,W1S,W2S,W1S,D1,DIS)
	D2 = D1*(ELN-AL-DIS)
	RS =-D2/ELN
	CALL TRAPIZ(XM1,DL,W1S,W2S,W1S,D1,DIS1)
	CALL TRAPIZ(XM2,DL,W1S,W2S,W2S,D2,DIS2)
	D3 = D1*(X  -AL-DIS1)
	D4 = D2*(X  -BL-DIS2)
	RMT = RS*X + D3 - D4
	RESUL(2) = RESUL(2) + D1 - D2		 !SHEAR S
	RESUL(6) = RESUL(6) + RMT          !BENDING T


C	FROM WT
	CALL TRAPIZ(DL ,DL,W1T,W2T,W1T,D1,DIS)
	D2 = D1*(ELN-AL-DIS)
	RT =-D2/ELN
	CALL TRAPIZ(XM1,DL,W1T,W2T,W1T,D1,DIS1)
	CALL TRAPIZ(XM2,DL,W1T,W2T,W2T,D2,DIS2)
	D3 = D1*(X  -AL-DIS1)
	D4 = D2*(X  -BL-DIS2)
	RMS = RT*X + D3 - D4
	RESUL(3) = RESUL(3) + D1 - D2		 !SHEAR T
	RESUL(5) = RESUL(5) - RMS          !BENDING S

	GOTO 200

100	CONTINUE

	DUM = 0.0
	IF(XM1.GT.0.0) DUM = 1.0 

C	FROM P MR
	RESUL(4) = RESUL(4) + W1MR*DUM  !TORSION

C	FROM P MS
	RT = -W1MS/ELN
	RMS = RT*X + W1MS*DUM
	RESUL(3) = RESUL(3)       !SHEAR T
	RESUL(5) = RESUL(5) + RMS !BENDING S

C	FROM P MT
	RS = -W1MT/ELN
	RMT = RS*X + W1MT*DUM
	RESUL(2) = RESUL(2)       !SHEAR S
	RESUL(6) = RESUL(6) - RMT !BENDING T


C	FROM P R
	RESUL(1) = RESUL(1) + W1R*DUM  !AXIAL

C	FROM P S
	D2 = W1S*(ELN-AL)
	RS =-D2/ELN
	RMT = RS*X   + W1S*XM1
	RESUL(2) = RESUL(2) + W1S*DUM       !SHEAR S
	RESUL(6) = RESUL(6) + RMT           !BENDING T

C	FROM P T
	D2 = W1T*(ELN-AL)
	RT =-D2/ELN
	RMS = RT*X   + W1T*XM1
	RESUL(3) = RESUL(3) + W1T*DUM       !SHEAR T
	RESUL(5) = RESUL(5) - RMS           !BENDING S


200	CONTINUE


	IRC = 0
	IRO = IP + NP*(MLE-1) + NP*NELE*(ILC-1)
	IF(TYP.EQ.'VARY') IRC = 2*IRO-1	      !VARY LOAD LINE NUMBER 2*N-1
	IF(TYP.EQ.'CONT') IRC = 2*IRO-0     !CONT LOAD LINE NUMBER 2*N
	IF(IRC.EQ.0) RETURN
	WRITE(NFL,REC=IRC) RESUL(1:7)
	
1000	CONTINUE
C	-----------------------------------------------------------	

	RETURN
      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	FUNCTION VMCAULY(X,A)          !MACAULAY FUNCTION
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	VMCAULY = X - A
	IF(VMCAULY.LT.0.0) VMCAULY = 0.0

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE TRAPIZ(X,D,A,B,C,AREA,DISA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	B > A

	W1 = C
	W2 = C + (B-A)*X/D

	AREA = 0.5*X*(W1+W2)

	A1 = X*W1
	A2 = 0.5*(W2-W1)*X
	Q1 = A1*0.5*X
	Q2 = A2*2.0*X/3.0

	DISA = 0.0
	IF(AREA.NE.0.0) DISA = (Q1+Q2)/AREA

	RETURN
	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SMHFAC(RHO,FIN,IEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	=======================================================================
C	LOADCASE
	COMMON /STCAS/ ILC
C	=======================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH
C	-----------------------------------------------------------------------
	DIMENSION FIN(1),FN1(7),FN2(7),FF(7),FO(7),FN(7)
C	-----------------------------------------------------------------------

	NP  = NP_SMH  
	NP2 = NP_SMH + 2  !2 IS FOR BOTH END

	CALL INTFILL('OGRF',NFL,13,KEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM


	DS = 1.0/(NP+1)

	FN1(1:7) = FIN(1:7 )
	FN2(1:7) = FIN(8:14)

	X = 0.0

	IP2 = 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (-FN1(I),I=1,7)
	
	
	DO 100 IP  = 1,NP
	X  = X + DS
	H1 = 1.0 - X
	H2 = X

	FF(1:7) = 0.0D0
	FO(1:7) = 0.0D0
	IRO = IP + NP*(IEL-1) + NP*NELE*(ILC-1)
	IF(ILC.NE.0) IRC = 2*IRO-1	      !VARY LOAD LINE NUMBER 2*N-1
	IF(ILC.NE.0) READ(NFL,REC=IRC) (FF(I),I=1,7)
	IF(ILC.NE.0) IRC = 2*IRO-0	      !CONT LOAD LINE NUMBER 2*N-1
	IF(ILC.NE.0) READ(NFL,REC=IRC) (FO(I),I=1,7)

	FF(1:7) = RHO*FF(1:7) + FO(1:7)

	FF(1) =- FF(1) - FN1(1)

	FF(2) =- FF(2) - FN1(2)

	FF(3) =- FF(3) - FN1(3)

	FF(4) =- FF(4) - FN1(4)

	FN(5) =- H1*FN1(5) + H2*FN2(5)
	FF(5) =  FF(5) + FN(5)

	FN(6) =- H1*FN1(6) + H2*FN2(6)
	FF(6) =  FF(6) + FN(6)

	FN(7) =- H1*FN1(7) + H2*FN2(7)
	FF(7) =  FF(7) + FN(7)

	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FF(I),I=1,7)

100	CONTINUE

	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FN2(I),I=1,7)


	RETURN

      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SMHFAC_OFFSHORE(RHO,FIN,IEL,ILOFF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	=======================================================================
C	LOADCASE
	COMMON /STCAS/ ILC
C	=======================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH
C	-----------------------------------------------------------------------
	DIMENSION FIN(1),FN1(7),FN2(7),FF(7),FO(7),FN(7)
C	-----------------------------------------------------------------------

	NP  = NP_SMH  
	NP2 = NP_SMH + 2  !2 IS FOR BOTH END

	CALL INTFILL('OGRF',NFL,20,KEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM


	DS = 1.0/(NP+1)

	FN1(1:7) = FIN(1:7 )
	FN2(1:7) = FIN(8:14)

	X = 0.0

	IP2 = 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (-FN1(I),I=1,7)
	
	
	DO 100 IP  = 1,NP
	X  = X + DS
	H1 = 1.0 - X
	H2 = X

	FF(1:7) = 0.0D0
	FO(1:7) = 0.0D0
	IRO = IP + NP*(IEL-1) + NP*NELE*(ILOFF-1)
	IF(ILOFF.NE.0) IRC = 2*IRO-1	      !VARY LOAD LINE NUMBER 2*N-1
	IF(ILOFF.NE.0) READ(NFL,REC=IRC) (FF(I),I=1,7)
	IF(ILOFF.NE.0) IRC = 2*IRO-0	      !CONT LOAD LINE NUMBER 2*N-1
	IF(ILOFF.NE.0) READ(NFL,REC=IRC) (FO(I),I=1,7)

	FF(1:7) = RHO*FF(1:7) + FO(1:7)

	FF(1) =- FF(1) - FN1(1)

	FF(2) =- FF(2) - FN1(2)

	FF(3) =- FF(3) - FN1(3)

	FF(4) =- FF(4) - FN1(4)

	FN(5) =- H1*FN1(5) + H2*FN2(5)
	FF(5) =  FF(5) + FN(5)

	FN(6) =- H1*FN1(6) + H2*FN2(6)
	FF(6) =  FF(6) + FN(6)

	FN(7) =- H1*FN1(7) + H2*FN2(7)
	FF(7) =  FF(7) + FN(7)

	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FF(I),I=1,7)

100	CONTINUE

	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FN2(I),I=1,7)


	RETURN

	END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE LNHFAC_OFFSHORE(FIN,IEG,RHO,IEL,NELE,LCS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /LANPRIN/  LAN_PRN,ILPL,MLANE,NPL,LAN_OPT
	COMMON /LANEFIX/ ALF1(50000),ALF2(50000),ALF3(50000),ILF(10000)
C	=======================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH
C	-----------------------------------------------------------------------
	DIMENSION FIN(1),FN1(7),FN2(7),FF(7),FN(7)
C	-----------------------------------------------------------------------

	NP  = NP_SMH  
	NP2 = NP_SMH + 2  !2 IS FOR BOTH END
	CALL INTFILL('OGRF',NFL,20,IEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM
	
	NMI = 0
	NMS = 0

	DS = 1.0/(NP+1)
	DO 2000 IPL = 1,NPL

	NEFX = ILF(NMI+4)

	IF(IPL.NE.ILPL) GOTO 1500

	FN1(1:7) = FIN(1:7 )
	FN2(1:7) = FIN(8:14)
	X = 0.0D0

	IP2 = 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (-FN1(I),I=1,7)
	

	DO 100 IP  = 1,NP
	X  = X + DS
	H1 = 1.0 - X
	H2 = X
	
C	-------------------------
	FF(1:7) = 0.0D0
	IF(LAN_OPT.EQ.0) THEN
C	-----------------
	DO IEFX = 1,NEFX
	NUMI = 1 +  4*(IEFX-1) + NMI
	IGM = ILF(NUMI+0)
	KEG = ILF(NUMI+1)
	MLE = ILF(NUMI+2)
	NN  = 7*(IP-1) + 7*NP*(IEFX-1) + NMS
	IF(IEL.EQ.MLE.AND.IEG.EQ.KEG) THEN
	DO J = 1,6
	FF(J) = FF(J) + ALF3(NN+J)
	ENDDO
	ENDIF
	ENDDO
C	-----------------
	ENDIF
C	-------------------------


	FF(1) =- RHO*FF(1) - FN1(1)

	FF(2) =- RHO*FF(2) - FN1(2)

	FF(3) =- RHO*FF(3) - FN1(3)

	FF(4) =- RHO*FF(4) - FN1(4)

	FN(5) =- H1*FN1(5) + H2*FN2(5)
	FF(5) =  RHO*FF(5) + FN(5)

	FN(6) =- H1*FN1(6) + H2*FN2(6)
	FF(6) =  RHO*FF(6) + FN(6)

	FN(7) =- H1*FN1(7) + H2*FN2(7)
	FF(7) =  RHO*FF(7) + FN(7)

	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FF(I),I=1,7)

100	CONTINUE


	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FN2(I),I=1,7)

	
1500	CONTINUE

	NMI = NMI +  4*NEFX
	NMS = NMS +  7*NP*NEFX
2000	CONTINUE


	RETURN

      END
C	=======================================================================
	SUBROUTINE LNHFAC(FIN,IEG,RHO,IEL,NELE,LCS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /LANPRIN/  LAN_PRN,ILPL,MLANE,NPL,LAN_OPT
	COMMON /LANEFIX/ ALF1(50000),ALF2(50000),ALF3(50000),ILF(10000)
C	=======================================================================
C	SMOOTH LINE DIAGRAM FOR 2 NODE BEAM AND FRAME (NUMBER OF STATION POINT)
	COMMON /BF_SMOTH/ NP_SMH
C	-----------------------------------------------------------------------
	DIMENSION FIN(1),FN1(7),FN2(7),FF(7),FN(7)
C	-----------------------------------------------------------------------

	NP  = NP_SMH  
	NP2 = NP_SMH + 2  !2 IS FOR BOTH END
	CALL INTFILL('OGRF',NFL,13,IEG,0) !ONLY FOR FRAME ELEMENT FOR SMOOOTH LINE DIAGRAM
	
	NMI = 0
	NMS = 0

	DS = 1.0/(NP+1)
	DO 2000 IPL = 1,NPL

	NEFX = ILF(NMI+4)

	IF(IPL.NE.ILPL) GOTO 1500

	FN1(1:7) = FIN(1:7 )
	FN2(1:7) = FIN(8:14)
	X = 0.0D0

	IP2 = 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (-FN1(I),I=1,7)
	

	DO 100 IP  = 1,NP
	X  = X + DS
	H1 = 1.0 - X
	H2 = X
	
C	-------------------------
	FF(1:7) = 0.0D0
	IF(LAN_OPT.EQ.0) THEN
C	-----------------
	DO IEFX = 1,NEFX
	NUMI = 1 +  4*(IEFX-1) + NMI
	IGM = ILF(NUMI+0)
	KEG = ILF(NUMI+1)
	MLE = ILF(NUMI+2)
	NN  = 7*(IP-1) + 7*NP*(IEFX-1) + NMS
	IF(IEL.EQ.MLE.AND.IEG.EQ.KEG) THEN
	DO J = 1,6
	FF(J) = FF(J) + ALF3(NN+J)
	ENDDO
	ENDIF
	ENDDO
C	-----------------
	ENDIF
C	-------------------------


	FF(1) =- RHO*FF(1) - FN1(1)

	FF(2) =- RHO*FF(2) - FN1(2)

	FF(3) =- RHO*FF(3) - FN1(3)

	FF(4) =- RHO*FF(4) - FN1(4)

	FN(5) =- H1*FN1(5) + H2*FN2(5)
	FF(5) =  RHO*FF(5) + FN(5)

	FN(6) =- H1*FN1(6) + H2*FN2(6)
	FF(6) =  RHO*FF(6) + FN(6)

	FN(7) =- H1*FN1(7) + H2*FN2(7)
	FF(7) =  RHO*FF(7) + FN(7)

	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FF(I),I=1,7)

100	CONTINUE


	IP2 = IP2 + 1
	IRC = IP2 + NP2*(IEL-1) + 2*NP*NELE*LCS
	WRITE(NFL,REC=IRC) (FN2(I),I=1,7)

	
1500	CONTINUE

	NMI = NMI +  4*NEFX
	NMS = NMS +  7*NP*NEFX
2000	CONTINUE


	RETURN

	END


C	=======================================================================
	SUBROUTINE STRTALIGN
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 NAME
C     READ STRUCTURAL ALIGNMENT DATA  

      COMMON /MGRAV/ NGRAV     
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
     
      DIMENSION V(3)

      READ(ITI,*)
      READ(ITI,*) NDATA
      

C     ALLOCATE STORE DATA ---- 6,NDATA   H1,H2,HG,V(1:3)     
	CALL DEFNREL('OREF',KDATA,6,NDATA)
	

      
      NMAX = 0
      DO I = 1,NDATA
      
      READ(ITI,*)
      READ(ITI,*) XXR,YYR,ZZR,V(1:3)  !REFERENCE POSITION & ALIGNMENT VECTOR
      ELN = SQRT( V(1)*V(1) +  V(2)*V(2) +  V(3)*V(3) )
      IF(ELN.EQ.0.0D0) THEN
      V(1:3) = 0.0D0
      V(NGRAV) = 1.0D0
      ELSE
      V(1:3) = V(1:3)/ELN 
      ENDIF
      ! BLOCK 12-01-2021
      !CALL RELFILL('OREF',XXR ,1,I,1)
      !CALL RELFILL('OREF',YYR ,2,I,1)
      !CALL RELFILL('OREF',ZZR ,3,I,1)
      !CALL RELFILL('OREF',V(1),4,I,1)
      !CALL RELFILL('OREF',V(2),5,I,1)
      !CALL RELFILL('OREF',V(3),6,I,1)
          
      ENDDO
      
      
	RETURN
      END
C	=======================================================================
C	======================================================================= 
      SUBROUTINE STOREOFFSHOREELEMENT (NNOFFL,IIVAL,RRVAL,LLRPIN,OPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPT
      DIMENSION IIVAL(20,NNOFFL),RRVAL(25,NNOFFL),LLRPIN(14,1)
C     STOREAGE OFFSHORE FRAME ELEMENT FOR LINEAR DYNAMIC ANALYSIS
      COMMON /OFFELE/ RVAL(25,200000),NOFFL,IVAL(20,200000),LRPIN(14,1)
      
      IF (OPT.EQ."STOR")THEN
      NOFFL               = NNOFFL
      IVAL(1:20,1:NOFFL)  = IIVAL(1:20,1:NOFFL)
      RVAL(1:25,1:NOFFL)  = RRVAL(1:25,1:NOFFL) 
      LRPIN(1:14,1)       = LLRPIN(1:14,1)
      ELSEIF (OPT.EQ."CALN")THEN
      NNOFFL              = NOFFL
      ELSEIF (OPT.EQ."CALT")THEN
      IIVAL(1:20,1:NOFFL) = IVAL(1:20,1:NOFFL)
      RRVAL(1:25,1:NOFFL) = RVAL(1:25,1:NOFFL) 
      LLRPIN(1:14,1)      = LRPIN(1:14,1)
      ENDIF
      
      
      END
C	=======================================================================
C	=======================================================================