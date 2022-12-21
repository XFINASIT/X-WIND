      SUBROUTINE HIROI (LM,XYZ,IGSET,PROPG,MTSET,PROPM,IGIDM,
	1				    MLOAD,LRPIN,IHSET)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON A(9000000),IA(9000000)
      
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
      
      COMMON /BF_SMOTH/ NP_SMH 
      
      COMMON /MGRAV/ NGRAV 
      
      COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
      
      COMMON /LCSS/ ILCN,ILCC
      
      DIMENSION   INAME(4)
      
      DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION XYZ(NCO*NNM,NELE),LM(NEF,NELE),R(NEQ)
	DIMENSION FIXF(7,2),IGIDM(NELE)
	DIMENSION FRMFOC(7,2),LRPIN(14,1),IHSET(1),IPIN(14)
	DIMENSION VR(3),WW(6),MLOAD(20)
      
      DIMENSION RMA(NEQ),VFOMA(3),FRMFOCM(7,2),FIXFM(7,2)
      
      DIMENSION LWCASE(7),FORCE_B(3),FORCE_A(3),WAVE(3)
      
      DIMENSION VB(3)
      
      DIMENSION NODE(2)
      
      
      NHIROI   = MLOAD(9)
      
C	======================================================================
C	HIROI FORMULA
	IF(NHIROI.EQ.0D0) GOTO 112
      READ (ITI,*)
	DO 111 IUF = 1,NHIROI

	CALL CLEARA (R,NEQ)
      
	FRMFOC = 0.0D0
      ! ===== INPUT DATA =====
      ! MLE               = NUMBER OF ELEMENT
      ! NOPTION           = SURFACE OPTION CONTROL (0,1)
      ! SURFACE_INPUT     = USER-DEFINED SURFACE VALUE
      ! OFFSELECT         = OFFSHORE PARAMETER 
      ! ILCN              = LOAD CASE NUMBER
      ! ILCC              = CONSTANT LOAD NUMBER

	!READ(ITI,*) MLE,WW(1),WW(2),WW(3),ECR,ECS,ECT,NFOMA,VFOMA(1:3),GRAV,LOC,LOE,IMOM,ILCN,ILCC 
      READ(ITI,*) MLE,offselect,NOPTION,SURFACE,ILCN,ILCC
      
      ! CALL GEOMETRY
      ! =======================================================================
      ! AUTOMATIC SURFACE PROPERTIES IS PROVIDED RECTAGULA AND CIRCULAR SECTION
      ! =======================================================================
      IF (NOPTION.EQ.1) CALL GEOMETRY_SURFACE (MLE,SURFACE)
      IF (NOPTION.EQ.0) SURFACE = SURFACE_INPUT
      
      CALL CALLNUMNODE_F (MLE,NODE(1),NODE(2)) 

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
	CALL MRELFIL(NAMEI,RANG ,1,5 ,0,'THIR') !
      
      ! CALLING OFFSHORE PARAMETER
      CALL OFFSPARA_READ_COLLECT_DATA
        
      CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,VB,HM,HW,HC,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                  VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,
     1                  WKF,CBF,WFC,TIME)
      
      ! MAXIMUM ELEVATION FOR CALCULATE FORCE
      IF (SEABED.LE.0.0D0) ELEVATION = (WDEPTH+1.25D0*WVHIGHT) + SEABED
      IF (SEABED.GT.0.0D0) ELEVATION = (WDEPTH+1.25D0*WVHIGHT) - SEABED
      
      COOR_X_1 = XYZ(1,MLE)
      COOR_X_2 = XYZ(4,MLE)
      COOR_Y_1 = XYZ(2,MLE)
      COOR_Y_2 = XYZ(5,MLE)
      COOR_Z_1 = XYZ(3,MLE)
      COOR_Z_2 = XYZ(6,MLE)
      
      IF     (NGRAV.EQ.1) THEN ! X-DIRECTION
             IF (COOR_X_1.GT.COOR_X_2) THEN
                 COOR_MAX = COOR_X_1
                 COOR_MIN = COOR_X_2
                 NMETHOD  = 2
             ENDIF
             IF (COOR_X_1.LE.COOR_X_2) THEN
                 COOR_MAX = COOR_X_2
                 COOR_MIN = COOR_X_1
                 NMETHOD  = 1
             ENDIF
      ELSEIF (NGRAV.EQ.2) THEN ! Y-DIRECTION
             IF (COOR_Y_1.GT.COOR_Y_2) THEN
                 COOR_MAX = COOR_Y_1
                 COOR_MIN = COOR_Y_2
                 NMETHOD  = 2
             ENDIF
             IF (COOR_Y_1.LE.COOR_Y_2) THEN
                 COOR_MAX = COOR_Y_2
                 COOR_MIN = COOR_Y_1
                 NMETHOD  = 1
             ENDIF
      ELSEIF (NGRAV.EQ.3) THEN ! Z-DIRECTION
             IF (COOR_Z_1.GT.COOR_Z_2) THEN
                 COOR_MAX = COOR_Z_1
                 COOR_MIN = COOR_Z_2
                 NMETHOD  = 2
             ENDIF
             IF (COOR_Z_1.LE.COOR_Z_2) THEN
                 COOR_MAX = COOR_Z_2
                 COOR_MIN = COOR_Z_1
                 NMETHOD  = 1
             ENDIF!
      ENDIF
      
      AAL = 1.0D0
      BBL = 1.0D0
      IF (NMETHOD.EQ.1) THEN
      IF    (COOR_MAX.GT.ELEVATION)THEN
            ! DOUBLE CHECK DATA
            IF (COOR_MIN.LE.ELEVATION) THEN
            ! CUT THE ELEMENT LENGTH
            ELN_TRUE = ELEVATION - COOR_MIN
            BBL = ELN_TRUE/ELN
            ELSEIF (COOR_MIN.GT.ELEVATION)THEN
            ! END OF IF LOOP
            ENDIF
      ELSEIF (COOR_MAX.LE.ELEVATION)THEN
      ! END OF IF LOOP    
      ENDIF
      
      ELSEIF (NMETHOD.EQ.2)THEN
      IF    (COOR_MAX.GT.ELEVATION)THEN
            ! DOUBLE CHECK DATA
            IF (COOR_MIN.LE.ELEVATION) THEN
            ! CUT THE ELEMENT LENGTH
            ELN_TRUE = ELEVATION - COOR_MIN
            AAL = ELN_TRUE/ELN
            ELSEIF (COOR_MIN.GT.ELEVATION)THEN
            ! END OF IF LOOP
            ENDIF
      ELSEIF (COOR_MAX.LE.ELEVATION)THEN
      ! END OF IF LOOP    
      ENDIF     
      ENDIF
          
      CALL BREAKING_ANGLE (WAVEANGLE,VB)
      CALL WAVE_VECTOR (WAVEANGLE,WAVE)
      
      DO I = 1,3
      FORCE_A(I) = 1.5D0*RHOW*WVHIGHT*SURFACE*WAVE(I)
      FORCE_B(I) = 1.5D0*RHOW*WVHIGHT*SURFACE*WAVE(I)
      ENDDO
      
      IF (ELEVATION.LT.0.0D0) THEN
      FORCE_A(1) = 0.0D0
      FORCE_B(1) = 0.0D0
      ENDIF
       

C	-----------------------------------------------------
	DO I = 1,3
	RP   = FORCE_A(I)
	RPB  = FORCE_B(I)
      IF (NMETHOD.EQ.2)THEN
      DISA = ELN*(1D0-AAL)
	DISB = ELN*BBL
      ELSEIF (NMETHOD.EQ.1)THEN
      DISA = 0.0D0
	DISB = ELN*BBL
      ENDIF 
	IDIR = I
	IND = 2
       
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
	CALL FRMBEP_OFFSHORE (FRMFOC,VR,KEG,MLE,RANG)

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
112   END
C     ==================================================
C     =================================================
      SUBROUTINE GEOMETRY_SURFACE (NELEMENT,SURFACE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      
      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      
      IF (NSECTION.EQ.1)  SURFACE = BREC
      IF (NSECTION.EQ.20) SURFACE = BBOX
      IF (NSECTION.EQ.15) SURFACE = DROUND/2.0D0
      IF (NSECTION.EQ.21) SURFACE = DPIPE/2.0D0
      
      END
C     ==================================================
C     =================================================
      SUBROUTINE GEOMETRY_SURFACE_B (NELEMENT,SURFACE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON / SECTIONDETAIL / BFISHAPE,TFISHAPE,DISHAPE,BWISHAPE,TWISHAPE,ROOTISHAPE
     1                        ,BFHSHAPE,TFHSHAPE,DHSHAPE,BWHSHAPE,TWHSHAPE,ROOTHSHAPE
     1                        ,BANGLE,TANGLE,HANGLE,THANGLE,AXBAR
     1                        ,BFCHANNEL,TFCHANNEL,DCHANNEL,TWCHANNEL,HCHANNEL,ROOTCSHAPE,CXBAR
     1                        ,BFTSHAPE,TFTSHAPE,DTSHAPE,BWTSHAPE,TWTSHAPE,TYBAR
     1                        ,BBOX,TFBOX,DBOX,HBOX,TWBOX,ROOTBOX
     1                        ,DPIPE,TPIPE
     1                        ,DROUND
     1                        ,BREC,HREC
     1                        ,SECTIONT,SECTIONS,AJ,GSECTION,CW,AREA,AIS,AIT,ARGS,ARGT
     1                        ,PLASTICT,PLASTICS,AMODULUS
      
      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      
      IF (NSECTION.EQ.1)  SURFACE = HREC
      IF (NSECTION.EQ.20) SURFACE = DBOX
      IF (NSECTION.EQ.15) SURFACE = DROUND/2.0D0
      IF (NSECTION.EQ.21) SURFACE = DPIPE/2.0D0
      
      END