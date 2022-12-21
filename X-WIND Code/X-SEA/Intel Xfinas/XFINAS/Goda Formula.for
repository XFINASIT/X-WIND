      SUBROUTINE GODA (LM,XYZ,IGSET,PROPG,MTSET,PROPM,IGIDM,
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
      
      DIMENSION LWCASE(7)
      
      DIMENSION ALPHA2_1(2)
      
      DIMENSION VB(3),WAVE(3)
      
      DIMENSION NODE(2)
      
      DIMENSION FORCE_A(3),FORCE_B(3)
      
      NGODA   = MLOAD(10)
      
C	======================================================================
C	GODA FORMULA
	IF(NGODA.EQ.0D0) GOTO 112
      READ (ITI,*)
	DO 111 IUF = 1,NGODA

	CALL CLEARA (R,NEQ)
      
	FRMFOC = 0.0D0
      ! ===== INPUT DATA =====
      ! MLE               = NUMBER OF ELEMENT
      ! NOPTION           = SURFACE OPTION CONTROL (0,1)
      ! SURFACE_INPUT     = USER-DEFINED SURFACE VALUE
      ! OFFSELECT         = OFFSHORE PARAMETER 
      ! ILCN              = LOAD CASE NUMBER
      ! ILCC              = CONSTANT LOAD NUMBER

C	LOC = 0 !LOCAL LOAD FLAG 0=GLOBAL 1=LOCAL
C	LOE = 0 !LOCAL ECC  FLAG 0=GLOBAL 1=LOCAL
C	IMOM = 0 !MOMENT LOAD FLAG 0=FORCE 1=MOMENT
C	ECR , ECS, ECT !ECCENTRICITY
	!READ(ITI,*) MLE,WW(1),WW(2),WW(3),ECR,ECS,ECT,NFOMA,VFOMA(1:3),GRAV,LOC,LOE,IMOM,ILCN,ILCC 
      READ(ITI,*) MLE,offselect,NOPTION,SURFACE,ILCN,ILCC
      
      CALL ARRAYELEMENT (NELEMENT,NSECTION)
      ! CALL GEOMETRY
      ! =======================================================================
      ! AUTOMATIC SURFACE PROPERTIES IS PROVIDED RECTAGULA AND CIRCULAR SECTION
      ! =======================================================================
      IF (NOPTION.EQ.1) CALL GEOMETRY_SURFACE (MLE,SURFACE)
      IF (NOPTION.EQ.1) CALL GEOMETRY_SURFACE_B (MLE,SURFACE_B)
      IF (NOPTION.EQ.0) SURFACE = SURFACE_INPUT
      
      CALL CALLNUMNODE_F (MLE,NODE(1),NODE(2)) 
      !CALL CALLNUMNODE (MLE,NODE)
      
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
      
      ! =========== INPUT DATA =============
      ! WAVE_ANGLE = WAVE ANGLE
      ! WVHIGHT    = WAVE HEIGHT
      ! GRAV       = GRAVITATIONAL
      ! WDEPTH     = WATER DEPTH
      ! HS         = SHALLOW MOUND
      ! HW         = 
      ! HC         = STRUCTURE ABOVE MSL
      ! RHOW       = WATER DENSITY
      HC     = 2.5D0
      ALPHA2 = 0.
      PI     = 3.141592654
      CALL BREAKING_ANGLE (WAVE_ANGLE,VB)
      CALL WAVE_VECTOR (WAVE_ANGLE,WAVE)
      
      HW                = HW + HC
      WVHIGHT_MAX       = WVHIGHT*1.8D0
      SURFACE_ELEVATION = 0.75D0*(1.0+COS(WAVE_ANGLE))*WVHIGHT_MAX
      OMEGA = 2.0*PI/PERIOD  
      CALL Airy_Wave_Number(GRAV,OMEGA,WDEPTH,RK)
      RK = 0.272D0
      ALPHA1      = 0.6D0+0.5D0*(((2D0*RK*WDEPTH)/(SINH(2D0*RK*WDEPTH)))**2D0)
      ALPHA2_1(1) = ((HM)/(3D0*WDEPTH))*((WVHIGHT_MAX/(WDEPTH-HM))**2D0)
      ALPHA2_1(2) = 2D0*(WDEPTH-HM)/WVHIGHT_MAX
      ALPHA2      = MINVAL(ALPHA2_1)
      ALPHA3      = 1D0-(((HW-HC)/WDEPTH)*(1D0-(1D0/COSH(RK*(WDEPTH)))))

      P1 = 0.5D0*(1+COS(WAVE_ANGLE))*(ALPHA1+ALPHA2*(COS(WAVE_ANGLE))**2D0)*RHOW*WVHIGHT_MAX*SURFACE
      P1_P2 = 0.5D0*(1+COS(WAVE_ANGLE))*(ALPHA1+ALPHA2*(COS(WAVE_ANGLE))**2D0)*RHOW*WVHIGHT_MAX
      P3 = ALPHA3*P1_P2*SURFACE
      PU = 0.5D0*(1+COS(WAVE_ANGLE))*ALPHA1*ALPHA3*RHOW*WVHIGHT_MAX*SURFACE_B
      
      ! ! CALAULATE LENGTH OF THE STRUCTURE ABOVE MSL
      ! IF (NOPTION.EQ.1)     THEN ! AUTOMATIC
      ! !  **** LIMIT OF THE PPROGRAM ****
      ! ! FOR A SINGLE ELEMENT 
      !     IF (NGRAV.EQ.1) THEN
      !        HC = WDEPTH - ABS(VR(1))
      !     ELSEIF (NGRAV.EQ.2)THEN
      !        HC = WDEPTH - ABS(VR(2))
      !     ELSEIF (NGRAV.EQ.3)THEN
      !        HC = WDEPTH - ABS(VR(3))
      !     ENDIF
      ! ELSEIF (NOPTION.EQ.2) THEN ! USER DEFINED
      !       HC  = HC_INPUT
      ! ENDIF
      
      IF (SURFACE_ELEVATION.GT.HC)THEN
      P2 = (1D0-HC/SURFACE_ELEVATION)*P1_P2*SURFACE
      ELSEIF (SURFACE_ELEVATION.LE.0)THEN
      P2 = 0.0D0    
      ENDIF
      
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
      FORCE_D = 0.0D0
      CALL BREAKING_ANGLE (WAVEANGLE,VB)
      CALL WAVE_VECTOR (WAVEANGLE,WAVE)
      TOTAL_LENGTH = HM + HW + HC
      
      IF (COOR_MIN.LE.WDEPTH.AND.COOR_MAX.LE.WDEPTH)THEN
         IF (COOR_MIN.LT.HM.AND.COOR_MAX.LT.HM)THEN
         FORCE_A = 0.
         FORCE_B = 0.    
         ELSEIF (COOR_MIN.GE.HM.AND.COOR_MAX.GE.HM)THEN
         DO I = 1,3
         FORCE_A(I) = P3*WAVE(I)
         FORCE_B(I) = P1*WAVE(I)
         ENDDO    
          IF (NMETHOD.EQ.1)THEN
          AAL = 0.0D0
          ELSEIF (NMETHOD.EQ.2)THEN
          BBL = 0.0D0
          ENDIF
         ENDIF
         FORCE_D = PU*0.5D0
      ELSEIF (COOR_MIN.GE.WDEPTH.AND.COOR_MAX.GE.WDEPTH)THEN
          IF (COOR_MIN.GT.TOTAL_LENGTH.AND.COOR_MAX.GT.TOTAL_LENGTH)THEN
          FORCE_A = 0.
          FORCE_B = 0.
          ELSEIF (COOR_MIN.LE.TOTAL_LENGTH.AND.COOR_MAX.LE.TOTAL_LENGTH)THEN
          DO I = 1,3
          FORCE_A(I) = P1*WAVE(I)
          FORCE_B(I) = P2*WAVE(I)
          ENDDO
          IF (NMETHOD.EQ.1)THEN
          AAL = 0.0D0
          ELSEIF (NMETHOD.EQ.2)THEN
          BBL = 0.0D0
          ENDIF
          ELSE
          FORCE_A = 0.
          FORCE_B = 0.    
          ENDIF
      ENDIF
      
!      ! INTERPORATION
!      HMC   = WDEPTH - HM
!      HWC   = HC + WDEPTH
!      HW_HM = HW + HM
!      TOTAL_LENGTH = HM + HW + HC
!      
!      IF (COOR_MIN.LT.HM) THEN
!          ! CHACK COOR_MAX
!          IF (COOR_MAX.GT.HM)THEN
!             IF (COOR_MAX.LE.HW_HM) THEN
!                  IF (NMETHOD.EQ.1)THEN
!                  ELN_TRUE = COOR_MAX - HM
!                  BBL = ELN_TRUE/ELN
!                  ELSEIF (NMETHOD.EQ.2)THEN
!                  ELN_TRUE = COOR_MAX - HM
!                  AAL = ELN_TRUE/ELN 
!                  ENDIF
!             ELSEIF (COOR_MAX.GT.HW_HM.AND.COOR_MAX.LE.TOTAL_LENGTH) THEN
!             ! FULL FORCE CAPACITY   
!             DO I = 1,3
!             FORCE_A(I) = P3*WAVE(I)*POINT
!             FORCE_B(I) = P1*WAVE(I)*POINT
!             ENDDO    
!
!             IF (NMETHOD.EQ.1)THEN
!             BBL = 0.0D0
!             ELSEIF (NMETHOD.EQ.2)THEN
!             ELN_TRUE = COOR_MAX - HM
!             AAL = 0.0D0
!             ENDIF
!             
!             ! REMAIN
!             POINT_LENGTH = HC
!             POINT_FORCE  = P1-P2
!             POINT_RATIO  = POINT_FORCE/POINT_LENGTH
!             POINT        = P1 - (COOR_MAX - HW_HM )*POINT_RATIO
!             DO I = 1,3
!             FORCE_B(I) = P1*WAVE(I)*POINT
!             FORCE_C(I) = P2*WAVE(I)*POINT
!             ENDDO
!             
!             IF (NMETHOD.EQ.1)THEN
!             ELN_TRUE = COOR_MAX - HM
!             BBL = ELN_TRUE/ELN
!             ELSEIF (NMETHOD.EQ.2)THEN
!             ELN_TRUE = COOR_MAX - HM
!             AAL = ELN_TRUE/ELN 
!             ENDIF
!             
!             ENDIF
!          ELSE IF (COOR_MAX.LE.HM)THEN
!          FORCE_A = 0.
!          FORCE_B = 0.
!          AAL     = 0D0
!          BBL     = 1.0D0
!          ENDIF
!          
!      ELSEIF (COOR_MIN.GE.HM.AND.COOR_MIN.LE.HW_HM) THEN
!
!      ELSEIF (COOR_MIN.GE.HW_HM.AND.COOR_MIN.LE.TOTAL_LENGTH) THEN
!          
!      ENDIF

      
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
      
      IF (IUF.EQ.NGODA) THEN
      ENDIF
	CALL FRMBEP_OFFSHORE (FRMFOC,VR,KEG,MLE,RANG)

	IF (NLS.NE.0) THEN
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MLE),A(LES),A(LED),
     1             A(LEI),FRMFOC,NSF,NNF,5)
      ENDIF
      
      IF (FORCE_D.NE.0.0)THEN
        IF (NGRAV.EQ.1)THEN
        IEQ    = LM(1,MLE)  
        R(IEQ) = R(IEQ) + FORCE_D
        ELSEIF (NGRAV.EQ.2)THEN
        IEQ    = LM(2,MLE)  
        R(IEQ) = R(IEQ) + FORCE_D
        ELSEIF (NGRAV.EQ.3)THEN
        IEQ    = LM(3,MLE)  
        R(IEQ) = R(IEQ) + FORCE_D 
        ENDIF
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

C     ------------------------------------------
      SUBROUTINE BREAKING_ANGLE (BREAKING_VA,VB)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /MGRAV/ NGRAV 
      DIMENSION VB(3)

      IF (NGRAV.EQ.1)     THEN  
      XWAVE  = VB(2)
      YWAVE  = VB(1)
      ZWAVE  = VB(3)
      ELSEIF (NGRAV.EQ.2) THEN
      XWAVE  = VB(1)
      YWAVE  = VB(2)
      ZWAVE  = VB(3)
      ELSEIF (NGRAV.EQ.3)THEN
      XWAVE  = VB(1)
      YWAVE  = VB(3)
      ZWAVE  = VB(2)   
      ENDIF
      
      PI     = 3.141592654
      
      BREAKING_VA_1=ACOS(XWAVE)*180D0/PI
      BREAKING_VA_2=ASIN(ZWAVE) *180D0/PI
      
      ! DOUBLE CHECK VALUE
      IF (BREAKING_VA_1.NE.BREAKING_VA_2) THEN
      BREAKING_VA = (BREAKING_VA_1 + BREAKING_VA_2)/2D0 
      ELSEIF (BREAKING_VA_1.EQ.BREAKING_VA_2) THEN
      BREAKING_VA = BREAKING_VA_1 
      ENDIF
      
      
!      IF (WAVEANGLE.EQ.0.0D0.OR.WAVEANGLE.EQ.360D0)THEN
!      XWAVE=1.0D0
!      YWAVE=0.0D0
!      ZWAVE=0.0D0
!      ELSEIF (WAVEANGLE.GT.0D0.AND.WAVEANGLE.LT.90D0)THEN
!      RAD=PI/180D0
!      XWAVE=COS(WAVEANGLE*RAD)
!      YWAVE=0.0D0
!      ZWAVE=SIN(WAVEANGLE*RAD)
!      ELSEIF (WAVEANGLE.EQ.90D0)THEN
!      XWAVE=0.0D0
!      YWAVE=0.0D0
!      ZWAVE=1.0D0
!      ELSEIF (WAVEANGLE.GT.90D0.AND.WAVEANGLE.LT.180D0)THEN
!      RAD=PI/180D0
!      XWAVE=1D0*COS(WAVEANGLE*RAD)
!      YWAVE=0.0D0
!      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
!      ELSEIF (WAVEANGLE.EQ.180)THEN
!      XWAVE=-1.0D0
!      YWAVE=0.0D0
!      ZWAVE=0.0D0
!      ELSEIF (WAVEANGLE.GT.180D0.AND.WAVEANGLE.LT.270D0)THEN
!      RAD=PI/180D0
!      XWAVE=1D0*COS(WAVEANGLE*RAD)
!      YWAVE=0.0D0
!      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
!      ELSEIF (WAVEANGLE.EQ.270D0)THEN
!      XWAVE=0.0D0
!      YWAVE=0.0D0
!      ZWAVE=-1.0D0
!      ELSEIF (WAVEANGLE.GT.270D0.AND.WAVEANGLE.LT.360)THEN
!      RAD=PI/180D0
!      XWAVE=1D0*COS(WAVEANGLE*RAD)
!      YWAVE=0.0D0
!      ZWAVE=1D0*SIN(WAVEANGLE*RAD)
!      ENDIF
!      WIND(1)=XWAVE
!      WIND(2)=YWAVE
!      WIND(3)=ZWAVE
      END