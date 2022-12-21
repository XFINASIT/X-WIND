C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SHEPSRE (PROPG,IGSET,XYZ,NODEX,LM,
	1					MGP,MXY,MEX,MEF,NLOAD,IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CONVERTES UNIFORMLY DISTRIBUTED GLOBAL TRACTIONS OR NORMAL
C     PRESSURES TO EQUIVALENT NODAL LOADS
C	-----------------------------------
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     IGSET(NELE)      = GEOMETRIC SET NUMBER
C     XYZ(MXY,NELE)    = NODAL COORDINATES FOR ELEMENTS
C     NODEX(NEX,NELE)  = LOCATIONS OF EXCESSIVE NODES
C     LM(NEF,NELE)     = EQUATION NUMBERS FOR ELEMENT D.O.F.
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================
C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)
      
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION H(9),P(2,9),XJI(4),FA(4),G(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,9),IGPOS(9),IGIDM(NELE)

	DIMENSION LHOPT(NELE),WTABZ(NELE)

C
      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)

C     -----------------------------------
      ALLOCATE(RAL(MEF),MAL(MEF))
C     -----------------------------------
      DO 900  ILOAD=1,NLOAD
      
      RAL(1:MEF) = 0.0D0
      MAL(1:MEF) = 0

	READ(ITI,*) IELE,LOPT,PREST,ZREF,ILCN,ILCC

	CALL ELEREODER(IGIDM,NELE,IELE)

      ISET = IGSET(IELE)
	THICK = PROPG(2,ISET)

      CALL CLEARA (RL,3*NNO)

	MMGR = 3
	MMGS = 3
	IF(NNO.EQ.3.OR.NNO.EQ.6) MMGS = 1   !SHELL 3 NODE

C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	SELECTCASE(NNO)
      CASE(3)!SHELL 3 NODE
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNO)
      CASE(6)!SHELL 3 NODE ONATE
      H = 0.0D0
      P = 0.0D0
      NNOM = 3
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNOM)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNOM)
	CASE DEFAULT
      IF(NNO.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX(1,IELE),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,IELE),NNO)            !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,IELE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	ENDSELECT
C     --------------------------------------------------------

	IF(LOPT.EQ.1) GOTO 501

C	UNIFORM NORMAL PRESSURE LOAD
	DO 500  INO=1,NNO
      FAC = PREST*H(INO)*WT*DET
      DO 480  I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 500
 480  RL(I,INO) = RL(I,INO) + VT(ICO)*FAC
 500  CONTINUE
	GOTO 800

 501	CONTINUE
C	HYDROSTATIC LOAD
	WALEV = ZREF
C	SHELL HYDROSTATIC
	GALEV = 0.0D0
	DO INO = 1,NNO
	JJ = 3*(INO-1) + NGRAV
	GALEV = GALEV + H(INO)*XYZ(JJ,IELE)
	ENDDO

	IF(GALEV.GT.WALEV) GOTO 800
	DTH = WALEV-GALEV

	DO 502 INO=1,NNO
	FAC =-H(INO)*WT*DET*PREST*DTH
      DO 481 I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 502
 481  RL(I,INO) = RL(I,INO) + VT(ICO)*FAC
 502  CONTINUE

 800  CONTINUE

	CALL FIXSHE (RL,KEG,IELE,NNO,3)
C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,4)
C
 700  IEFL = 0
      DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 600  CONTINUE

	CALL LDASEM_NEW (RAL,MAL,IEFL)

 900  CONTINUE

      DEALLOCATE(RAL,MAL)
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SHEPERE (PROPG,IGSET,XYZ,NODEX,LM,
	1					MGP,MXY,MEX,MEF,NLOAD,IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CONVERTES UNIFORMLY DISTRIBUTED GLOBAL TRACTIONS OR NORMAL
C     PRESSURES TO EQUIVALENT NODAL LOADS
C	-----------------------------------
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     IGSET(NELE)      = GEOMETRIC SET NUMBER
C     XYZ(MXY,NELE)    = NODAL COORDINATES FOR ELEMENTS
C     NODEX(NEX,NELE)  = LOCATIONS OF EXCESSIVE NODES
C     LM(NEF,NELE)     = EQUATION NUMBERS FOR ELEMENT D.O.F.
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================
C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)
      
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION H(9),P(2,9),XJI(4),FA(4),G(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,9),IGPOS(9)

	DIMENSION LHOPT(NELE),WTABZ(NELE)


      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)

C     -----------------------------------
      ALLOCATE(RAL(MEF),MAL(MEF))
C     -----------------------------------
      DO 900  ILOAD=1,NLOAD

      RAL(1:MEF) = 0.0D0
      MAL(1:MEF) = 0

	READ(ITI,*) IELE,IFACE,LOPT,PREST,ZREF,ILCN,ILCC

	CALL ELEREODER(IGIDM,NELE,IELE)

      ISET = IGSET(IELE)
	THICK = PROPG(2,ISET)

      CALL CLEARA (RL,3*NNO)

	MMGR = 3
C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 805  IGR=1,MMGR

	IF(NNO.EQ.3.OR.NNO.EQ.6) GOTO 601

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE2D(IFACE,IGR,RI,SI,WT,MMGR)
      IF(NNO.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX(1,IELE),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,IELE),NNO)             !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,IELE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R= 1.0
C	FACE 2; R=-1.0
C	FACE 3: S= 1.0
C	FACE 4: S=-1.0
C	---------------------
	VR = 0.0D0
	VS = 0.0D0
C	---------------------------------------------
	SELECTCASE (IFACE)
	CASE(4)
	DO II = 1,NNO
	VR(1) = VR(1) + P(2,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + P(2,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + P(2,II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(2)
	DO II = 1,NNO
	VR(1) = VR(1) - P(2,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) - P(2,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) - P(2,II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(1)
	DO II = 1,NNO
	VR(1) = VR(1) - P(1,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) - P(1,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) - P(1,II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(3)
	DO II = 1,NNO
	VR(1) = VR(1) + P(1,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + P(1,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + P(1,II)*XYZ(3*II-0,IELE)
	ENDDO
	ENDSELECT

	GOTO 602
C	---------------------------------------------
601	CONTINUE
C	SHELL 3 NODE EDGE LOADING
	RI = GLOC(IGR,MMGR)
	WT =  GWT(IGR,MMGR)
	H = 0.0D0  !SHAPE FUNTION FOR LINE INTEGRATION
	G = 0.0D0  !SHAPE FUNTION DERIVATIVE FOR LINE INTEGRATION
	CALL SHAP1D3(RI,H,G,2)
	CALL JACO1D3(IFACE,XYZ(1,IELE),G,DET,IELE,2)
	
	SELECTCASE(NNO)
      CASE(3)!SHELL 3 NODE
	  CALL RESHPT3(H,G,IFACE,NNO)
      CASE(6)!SHELL 3 NODE ONATE
        NNOM = 3
	  CALL RESHPT3(H,G,IFACE,NNOM)
      ENDSELECT
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R + S = 1.0
C	FACE 2; R= 0.0
C	FACE 3: S= 0.0
C	---------------------
	II = 2  !V.2-1
	JJ = 1
	VR(1) = XYZ(3*II-2,IELE) - XYZ(3*JJ-2,IELE)
	VR(2) = XYZ(3*II-1,IELE) - XYZ(3*JJ-1,IELE)
	VR(3) = XYZ(3*II-0,IELE) - XYZ(3*JJ-0,IELE)
	CALL SCALEN(VR,VR,DUM,3)
	II = 3  !V.3-1
	JJ = 1
	VS(1) = XYZ(3*II-2,IELE) - XYZ(3*JJ-2,IELE)
	VS(2) = XYZ(3*II-1,IELE) - XYZ(3*JJ-1,IELE)
	VS(3) = XYZ(3*II-0,IELE) - XYZ(3*JJ-0,IELE)
	CALL SCALEN(VS,VS,DUM,3)
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
	VR = 0.0D0
	VS = 0.0D0

	SELECTCASE (IFACE)
	CASE(2)
	DO II = 1,NNO
	VR(1) = VR(1) + G(II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + G(II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + G(II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(3)
	DO II = 1,NNO
	VR(1) = VR(1) - G(II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) - G(II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) - G(II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(1)
	DO II = 1,NNO
	VR(1) = VR(1) + G(II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + G(II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + G(II)*XYZ(3*II-0,IELE)
	ENDDO
	ENDSELECT

C	---------------------------------------------
602	CONTINUE
C	---------------------------------------------
	CALL SCALEN(VR,VR,DET,3)
	CALL VECPRD(VR,VT,VS)
	CALL SCALEN(VS,VS,DUM,3)
C     ---------------------------------------------

	IF(LOPT.EQ.1) GOTO 701

C	UNIFORM NORMAL PRESSURE LOAD
	DO 680  INO=1,NNO
      FAC = PREST*H(INO)*WT*DET*THICK
      DO 620  I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 680
 620  RL(I,INO) = RL(I,INO) + VS(ICO)*FAC
 680  CONTINUE

	GOTO 805

 701	CONTINUE
C	HYDROSTATIC LOAD
	WALEV = ZREF
C	SHELL HYDROSTATIC
	GALEV = 0.0D0
	DO INO = 1,NNO
	JJ = 3*(INO-1) + NGRAV
	GALEV = GALEV + H(INO)*XYZ(JJ,IELE)
	ENDDO

	IF(GALEV.GT.WALEV) GOTO 805
	DTH = WALEV-GALEV

	DO 705 INO=1,NNO
	FAC = -H(INO)*WT*DET*PREST*DTH*THICK
      DO 702 I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 705
 702  RL(I,INO) = RL(I,INO) + VS(ICO)*FAC
 705  CONTINUE


 805  CONTINUE
	
	CALL FIXSHE (RL,KEG,IELE,NNO,3)

C	-------------------------------------------------------
C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,4)
C
 700  IEFL = 0
      DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 600  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)

 900  CONTINUE

      DEALLOCATE(RAL,MAL)

C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SHETRAC (PROPG,IGSET,XYZ,NODEX,LM,
	1					MGP,MXY,MEX,MEF,NLOAD,IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CONVERTES UNIFORMLY DISTRIBUTED GLOBAL TRACTIONS OR NORMAL
C     PRESSURES TO EQUIVALENT NODAL LOADS
C	-----------------------------------
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     IGSET(NELE)      = GEOMETRIC SET NUMBER
C     XYZ(MXY,NELE)    = NODAL COORDINATES FOR ELEMENTS
C     NODEX(NEX,NELE)  = LOCATIONS OF EXCESSIVE NODES
C     LM(NEF,NELE)     = EQUATION NUMBERS FOR ELEMENT D.O.F.
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================
C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)      
      
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION H(9),P(2,9),XJI(4),FA(4),G(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,9),IGPOS(9),GG(3)

	DIMENSION LHOPT(NELE),WTABZ(NELE)


      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)

C     -----------------------------------
      ALLOCATE(RAL(MEF),MAL(MEF))
C     -----------------------------------
      DO 900  ILOAD=1,NLOAD

      RAL(1:MEF) = 0.0D0
      MAL(1:MEF) = 0

	READ(ITI,*) IELE,IFACE,GG(1),GG(2),GG(3),ILCN,ILCC

	CALL ELEREODER(IGIDM,NELE,IELE)

      ISET = IGSET(IELE)
	THICK = PROPG(2,ISET)

      CALL CLEARA (RL,3*NNO)

	MMGR = 3
C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 805  IGR=1,MMGR

	IF(NNO.EQ.3.OR.NNO.EQ.6) GOTO 601

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE2D(IFACE,IGR,RI,SI,WT,MMGR)
      IF(NNO.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX(1,IELE),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,IELE),NNO)             !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,IELE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R= 1.0
C	FACE 2; R=-1.0
C	FACE 3: S= 1.0
C	FACE 4: S=-1.0
C	---------------------
	VR = 0.0D0
	VS = 0.0D0
C	---------------------------------------------
	SELECTCASE (IFACE)
	CASE(4)
	DO II = 1,NNO
	VR(1) = VR(1) + P(2,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + P(2,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + P(2,II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(2)
	DO II = 1,NNO
	VR(1) = VR(1) - P(2,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) - P(2,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) - P(2,II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(1)
	DO II = 1,NNO
	VR(1) = VR(1) - P(1,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) - P(1,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) - P(1,II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(3)
	DO II = 1,NNO
	VR(1) = VR(1) + P(1,II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + P(1,II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + P(1,II)*XYZ(3*II-0,IELE)
	ENDDO
	ENDSELECT

	GOTO 602
C	---------------------------------------------
601	CONTINUE
C	SHELL 3 NODE EDGE LOADING
	RI = GLOC(IGR,MMGR)
	WT =  GWT(IGR,MMGR)
	H = 0.0D0  !SHAPE FUNTION FOR LINE INTEGRATION
	G = 0.0D0  !SHAPE FUNTION DERIVATIVE FOR LINE INTEGRATION
	CALL SHAP1D3(RI,H,G,2)
	CALL JACO1D3(IFACE,XYZ(1,IELE),G,DET,IELE,2)
	
	SELECTCASE(NNO)
      CASE(3)!SHELL 3 NODE
	  CALL RESHPT3(H,G,IFACE,NNO)
      CASE(6)!SHELL 3 NODE ONATE
        NNOM = 3
	  CALL RESHPT3(H,G,IFACE,NNOM)
      ENDSELECT
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R + S = 1.0
C	FACE 2; R= 0.0
C	FACE 3: S= 0.0
C	---------------------
	II = 2  !V.2-1
	JJ = 1
	VR(1) = XYZ(3*II-2,IELE) - XYZ(3*JJ-2,IELE)
	VR(2) = XYZ(3*II-1,IELE) - XYZ(3*JJ-1,IELE)
	VR(3) = XYZ(3*II-0,IELE) - XYZ(3*JJ-0,IELE)
	CALL SCALEN(VR,VR,DUM,3)
	II = 3  !V.3-1
	JJ = 1
	VS(1) = XYZ(3*II-2,IELE) - XYZ(3*JJ-2,IELE)
	VS(2) = XYZ(3*II-1,IELE) - XYZ(3*JJ-1,IELE)
	VS(3) = XYZ(3*II-0,IELE) - XYZ(3*JJ-0,IELE)
	CALL SCALEN(VS,VS,DUM,3)
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
	VR = 0.0D0
	VS = 0.0D0

	SELECTCASE (IFACE)
	CASE(2)
	DO II = 1,NNO
	VR(1) = VR(1) + G(II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + G(II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + G(II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(3)
	DO II = 1,NNO
	VR(1) = VR(1) - G(II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) - G(II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) - G(II)*XYZ(3*II-0,IELE)
	ENDDO
	CASE(1)
	DO II = 1,NNO
	VR(1) = VR(1) + G(II)*XYZ(3*II-2,IELE)
	VR(2) = VR(2) + G(II)*XYZ(3*II-1,IELE)
	VR(3) = VR(3) + G(II)*XYZ(3*II-0,IELE)
	ENDDO
	ENDSELECT

C	---------------------------------------------
602	CONTINUE
C	---------------------------------------------
	CALL SCALEN(VR,VR,DET,3)
	CALL VECPRD(VR,VT,VS)
	CALL SCALEN(VS,VS,DUM,3)
C     ---------------------------------------------


C	UNIFORM NORMAL PRESSURE LOAD
	DO 680  INO=1,NNO
      FAC = H(INO)*WT*DET*THICK
      DO 620  I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 680
 620  RL(I,INO) = RL(I,INO) + GG(ICO)*FAC
 680  CONTINUE


 805  CONTINUE


	CALL FIXSHE (RL,KEG,IELE,NNO,3)

C	-------------------------------------------------------
C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,4)
C
 700  IEFL = 0
      DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 600  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)

 900  CONTINUE

      DEALLOCATE(RAL,MAL)
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SHEPURE (PROPG,IGSET,XYZ,NODEX,LM,
	1					MGP,MXY,MEX,MEF,NLOAD,IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CONVERTES UNIFORMLY DISTRIBUTED GLOBAL TRACTIONS OR NORMAL
C     PRESSURES TO EQUIVALENT NODAL LOADS
C	-----------------------------------
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     IGSET(NELE)      = GEOMETRIC SET NUMBER
C     XYZ(MXY,NELE)    = NODAL COORDINATES FOR ELEMENTS
C     NODEX(NEX,NELE)  = LOCATIONS OF EXCESSIVE NODES
C     LM(NEF,NELE)     = EQUATION NUMBERS FOR ELEMENT D.O.F.
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC


C	==================================================================
C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)
      
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION H(9),P(2,9),XJI(4),FA(4),G(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,9),IGPOS(9),TT(3)

	DIMENSION LHOPT(NELE),WTABZ(NELE)

	DIMENSION HC1(3,9)  !C1 SHAPE FUNCTION


      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)

C     -----------------------------------
      ALLOCATE(RAL(MEF),MAL(MEF))
C     -----------------------------------
      DO 900  ILOAD=1,NLOAD

      RAL(1:MEF) = 0.0D0
      MAL(1:MEF) = 0

	READ(ITI,*) IELE,TT(1:3),ILCN,ILCC

	CALL ELEREODER(IGIDM,NELE,IELE)

      ISET = IGSET(IELE)
	THICK = PROPG(2,ISET)

      CALL CLEARA (RL,3*NNO)

	MMGR = 3
	MMGS = 3
	IF(NNO.EQ.3.OR.NNO.EQ.6) MMGS = 1   !SHELL 3 NODE

C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
      SELECTCASE(NNO)
      CASE(3)!SHELL 3 NODE
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNO)
      CASE(6)!SHELL 3 NODE ONATE
      H = 0.0D0
      P = 0.0D0
      NNOM = 3
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNOM)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNOM)
	CASE DEFAULT
      IF(NNO.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX(1,IELE),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,IELE),NNO)            !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,IELE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	ENDSELECT
C     --------------------------------------------------------


      DO 390  INO=1,NNO
      FAC = H(INO)*WT*DET
      DO 380  I=1,3
	ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 390
 380  RL(I,INO) = RL(I,INO) + TT(ICO)*FAC
 390  CONTINUE

 800  CONTINUE

	CALL FIXSHE (RL,KEG,IELE,NNO,3)
		
C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,4)
C
 700  IEFL = 0
      DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 600  CONTINUE
      
	CALL LDASEM_NEW (RAL,MAL,IEFL)

 900  CONTINUE
      
      DEALLOCATE(RAL,MAL)
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SHEOFFL (PROPG,IGSET,XYZ,NODEX,LM,
	1					MGP,MXY,MEX,MEF,NLOAD,IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CONVERTES OFFSHORE LOAD TO EQUIVALENT NODAL LOADS
C	-----------------------------------
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     IGSET(NELE)      = GEOMETRIC SET NUMBER
C     XYZ(MXY,NELE)    = NODAL COORDINATES FOR ELEMENTS
C     NODEX(NEX,NELE)  = LOCATIONS OF EXCESSIVE NODES
C     LM(NEF,NELE)     = EQUATION NUMBERS FOR ELEMENT D.O.F.
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================
C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

C     ============================================================================
C     OFFSHORE PART COMMON 
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION
	COMMON /WARNING/ WARNING,RAMDA(100),RK
      COMMON /DIFFRACTIONCOORDINATE / DIFFRACT_X , DIFFRACT_Y , DIFFRACT_GRAVITY 
C     ============================================================================

      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION R(NEQ),H(9),P(2,9),XJI(4),FA(4),G(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,9),IGPOS(9),IGIDM(NELE)

C     FOR OFFSHORE LOAD
      DIMENSION VWIND(3),VH1(3),VH2(3),VGV(3),LWCASE(10),TT(3),XYZCP(3)
      
C     FOR OFFSHORE LOAD -- PRAMIN 
      DIMENSION FLIST(5),VREF(3),VREW(3),VNOL(3),VNOLG(3),VCURRENTL(5)
      
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:)
	      

      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)


	IF(NLOAD.EQ.0) GOTO 1000	
      	
	PI = 3.141592654

	ALLOCATE(RVAL(25,NLOAD),IVAL(20,NLOAD))
	
C     READ INPUT DATA AND STORE ON TEMPORARY VARIABLE
	DO ILOAD=1,NLOAD
        READ(ITI,*) IELE,NUMCP,NFUNCTION,ROUGH,CD,CM,CSA,DIAM1,SHFAC,GUST_FACTOR,NGROWTH,GROWTH,NIRRE,OFFSELECT_X,LWCASE(1)
     1             ,LWCASE(2),LWCASE(5),LWCASE(6),LWCASE(7),ILCN,ILCC 
        IVAL(1 ,ILOAD) = IELE
        IVAL(2 ,ILOAD) = NUMCP
        IVAL(3 ,ILOAD) = LWCASE(1) ! WAVE LOAD
        IVAL(4 ,ILOAD) = LWCASE(2) ! CURRENT LOAD
        IVAL(5 ,ILOAD) = 0.0 !LWCASE(3) ! TIDAL LOAD (CANCEL)
        IVAL(6 ,ILOAD) = 0.0 !LWCASE(4) ! LEVEL LOAD (CANCEL)
        IVAL(7 ,ILOAD) = LWCASE(5) ! PLUNGING LOAD (BREAKLING WAVE)
        IVAL(8 ,ILOAD) = LWCASE(6) ! SURING LOAD (BREAKLING WAVE)
        IVAL(9 ,ILOAD) = LWCASE(7) ! WIND LOAD
        IVAL(10,ILOAD) = ILCN
        IVAL(11,ILOAD) = ILCC
        IVAL(12,ILOAD) = NIRRE
        
        RVAL(1 ,ILOAD) = CD
        RVAL(2 ,ILOAD) = CM
        RVAL(3 ,ILOAD) = GUST_FACTOR
        IF (LWCASE(7).NE.1) GUST_FACTOR = 1.0D0 ! PROTECT OTHER CASE =  0.0D0
        DIAM2=DIAM1+(GROWTH*2D0)
        RVAL(12 ,ILOAD) = DIAM2   !NORMINAL DIAMETER
        RVAL(13 ,ILOAD) = SHFAC  !SHAPE FACTOR = 0.25*PI for circular section
        RVAL(21,ILOAD)  = OFFSELECT_X  ! CASE SELECT
        RVAL(22,ILOAD)  = NFUNCTION
        RVAL(23,ILOAD)  = ROUGH
        RVAL(24,ILOAD)  = DIAM1
        RVAL(25,ILOAD)  = CSA 
        ! -------------------------------  FOR ERROR MESSAGE MARINE GROWTH ------------------------------- 
!             IF (NGROWTH.EQ.0.0D0)THEN
!                IF (GROWTH.NE.0.0D0)THEN
!                   WRITE (*,59)
!                   WRITE (*,60)
!                   WRITE (*,61)
!                   WRITE (*,59)
!                   WRITE (*,62)
!                   WRITE (*,63)
!                   WRITE (*,64)
!59                 FORMAT ('')                   
!60                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( FRAME ELEMENT )----------------')
!61                 FORMAT ('- INSIDE MARINE GROWTH FUNCTION HAVE SOME VALUE. II WILL EFFECT ON THE RESULT')
!62                 FORMAT ('**** COMMENT *****')
!63                 FORMAT (' CALCULATE     = SELECT MARINE GROWTH BLOCK')
!64                 FORMAT (' NOT CALAULATE = UNSELECT MARINE GROWTH BLOCK')
!                  STOPFUNCTION=1D0 
!                  GOTO 101
!                ELSEIF (GROWTH.EQ.0.0D0)THEN
!                ! CONTINUE 
!                ENDIF
!             ELSEIF (NGROWTH.NE.0.0D0)THEN
!                IF (GROWTH.GT.0.0D0)THEN
!                ! CONTINUE
!                ELSEIF (GROWTH.LE.0.0D0)THEN
!                   WRITE (*,65)
!                   WRITE (*,66)
!                   WRITE (*,67)
!                   WRITE (*,65)
!                   WRITE (*,68)
!                   WRITE (*,69)
!                   WRITE (*,70)
!65                 FORMAT ('')
!66                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( FRAME ELEMENT )----------------')
!67                 FORMAT ('- MARINE GROWTH MUST BE GRATER THAN ZERO ')
!68                 FORMAT ('**** COMMENT *****')
!69                 FORMAT (' CALCULATE     = SELECT MARINE GROWTH BLOCK')
!70                 FORMAT (' NOT CALAULATE = UNSELECT MARINE GROWTH BLOCK')
!                   ! STOPFUNCTION USE FOR STOP PROGRAM WHEN ERROR MESSAGE APPEAR
!                   STOPFUNCTION=1D0 
!                   ! GOTO 101 : OUT READ INPUT FOR STOP THE PROGRAM
!                   GOTO 101
!                ENDIF
!             ENDIF
         ! -------------------------------------------------------------------------------------------------- 
	ENDDO
		
	STOPFUNCTION=0.0D0 ! SET CONDITION PROTECT ERROR EFFECT
!	! ---- WARNING MESSAGE ---
!101   IF (LWCASE(1).EQ.0.AND.LWCASE(2).EQ.0.AND.LWCASE(3).EQ.0.AND.LWCASE(4).EQ.0.AND.LWCASE(5).EQ.0.AND.LWCASE(6).EQ.0
!     1    .AND.LWCASE(7).EQ.0.0D0)THEN
!      WRITE (*,8)
!      WRITE (*,34)
!      WRITE (*,35)
!34    FORMAT ('---------------- OFFSHORE WARNING MESSAGE ----------------')
!35    FORMAT ('- PLEASE CHECK OFFSHORE LOAD TYPES')
!      STOP
!      ENDIF
!      IF (NFUNCTION.EQ.3)THEN
!      IF (CD.LE.0.0D0.OR.CM.LE.0.0D0.OR.DIAM2.LE.0.0D0)THEN
!      WRITE (*,8)
!      WRITE (*,36)
!36    FORMAT ('---------------- OFFSHORE ERROR MESSAGE ----------------')
!          IF (CD.LE.0.0D0)THEN
!          WRITE (*,37)
!37        FORMAT ('- DRAG COEFFICIENT MUST BE GRATER THAN ZERO')
!          ENDIF
!          IF (CM.LE.0.0D0)THEN
!          WRITE (*,38)
!38        FORMAT ('- INERTIA COEFFICIENT MUST BE GRATER THAN ZERO')
!          ENDIF
!          IF (DIAM2.LE.0.0D0)THEN
!          WRITE (*,39)
!39        FORMAT ('- NORMINAL DIMAETER MUST BE GRATER THAN ZERO')
!          ENDIF
!      STOPFUNCTION=1D0    
!      ENDIF
!      ENDIF
	
	
	IF(KOFFL.EQ.0) THEN
	  DEALLOCATE(RVAL,IVAL)
	  GOTO 1000 !IF NO OFFSHORE LOAd ANALYSIS NO NEED TO DO THE CALCULATION...JUST READ THE INPUT DATA
	ENDIF
	
      
	ALLOCATE(RRV(NEQ),RRC(NEQ))

      IF(KOFFL.NE.0) THEN !CALL ALL OFFSHORE LOAD PARAMETERS
    !  CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,PERIOD,GRAV,RHOW,RHOA,
    ! 1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
    ! 1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED)
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      ENDIF  
      	
C     -----------------------------------
C     -----------------------------------

      DO 7500 ILCAS = 1,LCS !LOOP OVER LOAD CASE NUBER
      
      DO 7500 ITIME = 1,NTIME !LOOP OVER OFFSHORE TIME STEP

      RRV(1:NEQ) = 0.0D0 !INITIALIZE VARY OFFSHORE LOAD
      RRC(1:NEQ) = 0.0D0 !INITIALIZE CONSTANT OFFSHORE LOAD
      
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')
      CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','READ')
            
      DO 900  ILOAD=1,NLOAD

C     ------------------------------
C     CALL INPUT DATA FROM BACKUP
        IELE        = IVAL(1 ,ILOAD)
        NUMCP       = IVAL(2 ,ILOAD)
        LWCASE(1)   = IVAL(3 ,ILOAD)
        LWCASE(2)   = IVAL(4 ,ILOAD)
        LWCASE(3)   = IVAL(5 ,ILOAD)
        LWCASE(4)   = IVAL(6 ,ILOAD)
        LWCASE(5)   = IVAL(7 ,ILOAD)
        LWCASE(6)   = IVAL(8 ,ILOAD)
        LWCASE(7)   = IVAL(9 ,ILOAD)
        ILCN        = IVAL(10,ILOAD)
        ILCC        = IVAL(11,ILOAD)
        IORRE       = IVAL(12,ILOAD)
        
        GUSTFACTOR  = RVAL(3 ,ILOAD) 
        NFUNCTION   = RVAL(22,ILOAD) 
        ROUGH       = RVAL(23,ILOAD)
        DIAM        = RVAL(12,ILOAD)
        DIAM1       = RVAL(24,ILOAD)
        SHFAC       = RVAL(13,ILOAD)
        CSA         = RVAL(25,ILOAD)
        
        IF (NFUNCTION.EQ.2.0D0)THEN
        CD          = RVAL(1 ,ILOAD)
        CM          = RVAL(2 ,ILOAD)
        ELSEIF (NFUNCTION.EQ.1.0)THEN
        CD          = 0D0 
        CM          = 0D0
        ENDIF
        


C    ==================================================================   
C    ======== Chana Modifile offshore loadcase selecttion
C    ======== 5 / July / 2012 
C    ==================================================================   
      offselect = RVAL(21,ILOAD) 

      IF (IORRE.EQ.1) THEN
      CALL OFFSPARA_READ_COLLECT_DATA
       
      CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,VB,HM,HW,HC,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                  VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,
     1                  WKF,CBF,WFC,TIME)
      WARNING=GRAV    
     
C    ================================================================== 
       ! --- ERROR MESSAGE  --- 20 / JULY / 2012
!       IF (LWCASE(1).EQ.1.0D0)THEN ! WAVE LOAD
!          IF (WVHIGHT.EQ.0.0D0.OR.WDEPTH.EQ.0.0D0.OR.PERIOD.EQ.0.0D0)THEN
!          WRITE (*,8)
!          WRITE (*,4)
!4         FORMAT ('----- OFFSHORE ERROR MESSAGE ( WAVE LOAD FOR SHELL )-- LOAD PARAMETER',F4.0)
!             IF (WVHIGHT.EQ.0.0D0)THEN
!             WRITE (*,1)
!1            FORMAT ('- WAVE HEIGHT MUST BE GREATER THAN ZERO')
!             ENDIF
!             IF (WDEPTH.EQ.0.0D0)THEN
!             WRITE (*,2)
!2            FORMAT ('- WATER DEPTH MUST BE GREATER THAN ZERO')
!             ENDIF
!             IF (PERIOD.EQ.0.0D0)THEN
!             WRITE (*,3)
!3            FORMAT ('- PERIOD MUST BE GREATER THAN ZERO')
!             ENDIF
!          STOPFUNCTION=1D0
!          ENDIF
!      ENDIF
      
8     FORMAT ('') ! BLANK

!      IF (LWCASE(2).EQ.1.0D0)THEN ! CURRENT LOAD
!      ! NCURRENT = 1 ; DNV
!      ! NCURRENT = 2 ; API
!      ! NCURRENT = 3 ; IEC
!      ! NCURRENT = 4 ; POWER LAW
!      ! NCURRENT = 5 ; LINEAR CURRENT 5 VALUE
!      
!         IF (NCURRENT.EQ.1.OR.NCURRENT.EQ.3)THEN
!            IF (VTIDE.EQ.0.0D0.OR.VWIND0.EQ.0.0D0.OR.FACTOR.EQ.0.0D0)THEN
!               WRITE (*,8)
!               WRITE (*,7) OFFSELECT
!7              FORMAT ('----- OFFSHORE ERROR MESSAGE FOR SHELL ( CURRENT PAMETERS FOR FRAME )-- LOAD PARAMETER',F4.0)
!               IF (VTIDE.EQ.0.0D0)THEN
!               WRITE (*,5)
!5              FORMAT ('- WIND-GENERATE CURRENT VELOCITY ( DNV,IEC ) MUST BE GREATER THAN ZERO')
!               ENDIF
!               IF (VWIND0.EQ.0.0D0)THEN
!               WRITE (*,6)
!6              FORMAT ('- WIND-GENERATE TIDAL VELOCITY ( DNV,IEC ) MUST BE GREATER THAN ZERO')
!               ENDIF
!               IF (FACTOR.EQ.0.0D0)THEN
!               WRITE (*,24)             
!24             FORMAT ('- WIND-GENERATE CURRENT ( DNV,IEC ) MUST BE GREATER THAN ZERO')
!               ENDIF
!               STOPFUNCTION=1D0
!            ENDIF
!            
!          ELSEIF (NCURRENT.EQ.2)THEN
!             IF (VCURRENTAPI.EQ.0.0D0)THEN
!                 WRITE (*,8)
!                 WRITE (*,25) OFFSELECT
!25               FORMAT ('----- OFFSHORE ERROR MESSAGE FOR SHELL ( CURRENT PAMETERS )-- LOAD PARAMETER',F4.0)
!                 WRITE (*,26)
!26               FORMAT ('- CURRENT VELOCITY (API) MUST BE GREATER THAN ZERO')
!             STOPFUNCTION=1D0
!             ENDIF
!             
!          ELSEIF (NCURRENT.EQ.4)THEN
!             IF (POWERLAW.EQ.0.0D0.OR.VCURRENTP.EQ.0.0D0)THEN
!                 WRITE (*,8)
!                 WRITE (*,27) OFFSELECT
!27               FORMAT ('----- OFFSHORE ERROR MESSAGE FOR SHELL ( CURRENT PAMETERS )-- LOAD PARAMETER',F4.0)
!                 WRITE (*,28)
!28               FORMAT ('- CURRENT VELOCITY ( POWER LAW ) MUST BE GREATER THAN ZERO')
!             STOPFUNCTION=1D0
!             ENDIF
!          ELSEIF (NCURRENT.EQ.5)THEN
!             IF (VCURRENTL(1).EQ.0.0D0)THEN
!                WRITE (*,8)
!                WRITE (*,29) OFFSELECT
!29              FORMAT ('----- OFFSHORE ERROR MESSAGE FOR SHELL ( CURRENT PAMETERS )-- LOAD PARAMETER',F4.0)
!                WRITE (*,40)
!40              FORMAT ('- CURRENT VELOCITY ( LINEAR CURRENT ) MUST BE GREATER THAN ZERO')
!             STOPFUNCTION=1D0
!             ENDIF
!          ENDIF
!      ENDIF
      
!      IF (LWCASE(3).EQ.1.0D0)THEN ! TIDAL LOAD
!        IF (THIGHT.EQ.0.0D0)THEN
!            WRITE (*,8)
!            WRITE (*,13) OFFSELECT
!            WRITE (*,14)
!13          FORMAT ('----- OFFSHORE ERROR MESSAGE ( TIDAL LOAD FOR FRAME )-- LOAD PARAMETER',F4.0)
!14          FORMAT ('- TIDAL HEIGHT FROM NORMAL WATER SURFACE MUST BE GREATER THAN ZERO')
!            STOPFUNCTION=1D0
!        ENDIF
!      ENDIF
      
!      IF (LWCASE(4).EQ.1.0D01)THEN ! PLUNGING LOAD
      
!      ENDIF
      
!      IF (LWCASE(4).EQ.1.0D01)THEN ! PLUNGING LOAD
!      
!      ENDIF
!      
!      IF (LWCASE(7).EQ.1.0D0)THEN ! WIND LOAD
!          ! NWINDO = 1 ; LOGARITHMIC PROFILE   (DNV) 
!          ! NWINDO = 2 ; POWER LAW PROFILE     (DNV)
!          ! NWINDO = 3 ; WIND PROFILE AND GUST (API)
!          ! NWINDO = 4 ; IEC 614000-1
!          ! NWINDO = 5 ; USER DEFINED
!       IF(NWINDO.EQ.1.OR.NWINDO.EQ.2)THEN
!         IF (ALPHA.EQ.0.0D0.OR.UH.EQ.0.0D0.OR.RHIGH.EQ.0.0D0)THEN
!         WRITE (*,8)
!         WRITE (*,9) OFFSELECT
!9        FORMAT ('----- OFFSHORE ERROR MESSAGE ( WIND LOAD FOR SHELL )-- LOAD PARAMETER',F4.0)
!            IF(UH.EQ.0.0D0)THEN !
!            WRITE (*,10)
!10          FORMAT ('- MEAN WIND SPEED (DNV) MUST BE GREATER THAN ZERO')
!            ENDIF
!            IF (ALPHA.EQ.0.0)THEN
!            WRITE (*,11)
!11          FORMAT ('- EXPONENT IN POWER-LAW MODEL FOR WIND SPEED PROFILE MUST BE GREATER THAN ZERO')
!            ENDIF
!            IF (RHIGH.EQ.0.0D0)THEN
!            WRITE (*,15)
!15          FORMAT ('- REFERENCE HIGHT MUST BE GREATER THAN ZERO')  
!            ENDIF
!         STOPFUNCTION=1D0
!         ENDIF
!        ELSEIF (NWINDO.EQ.3)THEN
!           IF (UHAPI.EQ.0.0D0.OR.AVERAGE.EQ.0.0D0.OR.RHIGH.EQ.0.0D0)THEN
!           WRITE (*,8)
!           WRITE (*,16) OFFSELECT
!16         FORMAT ('----- OFFSHORE ERROR MESSAGE ( WIND LOAD FOR SHELL )-- LOAD PARAMETER',F4.0)
!           IF (UHAPI.EQ.0.0D0)THEN
!           WRITE (*,17)
!17         FORMAT ('- MEAN WIND SPEED (API) MUST BE GREATER THAN ZERO')             
!           ENDIF
!           IF (AVERAGE.EQ.0.0D0)THEN
!           WRITE (*,18)
!18         FORMAT ('- AVERAGING TIME MUST BE GREATER THAN ZERO AND LESS THAN 3600 SEC')             
!           ENDIF
!           IF (RHIGH.EQ.0.0D0)THEN
!           WRITE (*,19)
!19         FORMAT ('- REFERENCE HIGHT MUST BE GREATER THAN ZERO')  
!           ENDIF
!           STOPFUNCTION=1D0
!           ENDIF
!        ELSEIF (NWINDO.EQ.4.OR.NWINDO.EQ.5)THEN
!         IF (ALPHA.EQ.0.0D0.OR.UHD.EQ.0.0D0.OR.RHIGH.EQ.0.0D0)THEN
!         WRITE (*,8)
!         WRITE (*,20) OFFSELECT 
!20        FORMAT ('----- OFFSHORE ERROR MESSAGE ( WIND LOAD FOR SHELL )-- LOAD PARAMETER',F4.0)
!            IF(UHD.EQ.0.0D0)THEN !
!            WRITE (*,21)
!21          FORMAT ('- MEAN WIND SPEED MUST BE GREATER THAN ZERO')
!            ENDIF
!            IF (ALPHA.EQ.0.0)THEN
!            WRITE (*,21)
!22          FORMAT ('- EXPONENT IN POWER-LAW MODEL FOR WIND SPEED PROFILE MUST BE GREATER THAN ZERO')
!            ENDIF
!            IF (RHIGH.EQ.0.0D0)THEN
!            WRITE (*,23)
!23          FORMAT ('- REFERENCE HIGHT MUST BE GREATER THAN ZERO')  
!            ENDIF   
!            STOPFUNCTION=1D0      
!            ENDIF
!        ENDIF
!      ENDIF
      ! ----------------------------------------------------------------------------------------------
      
      ! --- STOP FUNCTION ---
      ! IF (STOPFUNCTION.EQ.1D0)THEN
      ! STOP
      ! ENDIF
      ! --- END STOP FUNCTION -- 
C     ------------------------------
      ELSEIF (IORRE.EQ.2)THEN
              
      CALL SELECTSPECTRUM (offselect,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
      ENDIF
      
      IF(ILCN.NE.ILCAS.AND.ILCC.NE.ILCAS) GOTO 900
      
      
	
	CALL ELEREODER(IGIDM,NELE,IELE)

      ISET = IGSET(IELE)
	THICK = PROPG(2,ISET)


C     ===================================================	
C     PRAMIN WAVE LOAD MODIFICATION
C     ===================================================

C	WAVE REFERNCE COORDINATE
      H1REF = H1POS
      H2REF = H2POS
      HGREF = SEABED		

C     ----------------------------------------------------
C     ELEMENT MID COORDINATES
	XMID = 0.0D0
      YMID = 0.0D0
      ZMID = 0.0D0
	DO INO = 1,NNO
	  !DO I = 1,3
	  JX = 3*(INO-1) + 1
	  JY = 3*(INO-1) + 2
	  JZ = 3*(INO-1) + 3
	  XMID = XMID + XYZ(JX,IELE)
	  YMID = YMID + XYZ(JY,IELE)
	  ZMID = ZMID + XYZ(JZ,IELE)
	  !ENDDO
	ENDDO
      XMID = XMID/NNO
      YMID = YMID/NNO
      ZMID = ZMID/NNO

C     ELEMENT MID COORDINATES IN WAVE DIRECTION SYSTEM
      HM1 = VH1(1)*XMID + VH1(2)*YMID + VH1(3)*ZMID 
      HMG = VGV(1)*XMID + VGV(2)*YMID + VGV(3)*ZMID 
      HM2 = VH2(1)*XMID + VH2(2)*YMID + VH2(3)*ZMID 
C     ----------------------------------------------------
      
      
! C     ----------------------------------------------------
! C     STRUCTURAL ALIGNMENT REFERENCE PROPERTIES
!       CALL OFFREFAXIS(NUMCP,XXR,YYR,ZZR,VREF,'OREF')
!       
! C     STRUCTURAL ALIGNMENT REFERENCE POSITION IN WAVE DIRECTION SYSTEM 
!       HH1 = VH1(1)*XXR + VH1(2)*YYR + VH1(3)*ZZR 
!       HHG = VGV(1)*XXR + VGV(2)*YYR + VGV(3)*ZZR 
!       HH2 = VH2(1)*XXR + VH2(2)*YYR + VH2(3)*ZZR 
! 
! C     STRUCTURAL ALIGNMENT REFERENCE VECTOR IN WAVE DIRECTION SYSTEM 
!       VREW(1) = VH1(1)*VREF(1) + VH1(2)*VREF(2) + VH1(3)*VREF(3)
!       VREW(2) = VGV(1)*VREF(1) + VGV(2)*VREF(2) + VGV(3)*VREF(3) 
!       VREW(3) = VH2(1)*VREF(1) + VH2(2)*VREF(2) + VH2(3)*VREF(3) 
! C     ----------------------------------------------------
!       
!       
! C     ----------------------------------------------------
! C     POSITION OF COORDINATES REFER TO WAVE REFERNCE COORDINATE (IN WAVE DIRECTION SYSTEM) CORESPONDING TO STRUCTURAL ALIGNMENT
!       HIG = HMG - HHG   !VERTICAL DISTANCE FROM STRUCTURAL REFERENCE COORDINATE
!       RLN = HIG/VREW(2) !LENGTH ALONG ALIGNMENT VECTOR
!       H1M =  HH1 + RLN*VREW(1) - H1REF !H1 DISTANCE REFER TO WAVE REFERNCE COORDINATE
!       HGM =  HHG + RLN*VREW(2) - HGREF !GRAVITY DISTANCE REFER TO WAVE REFERNCE COORDINATE
!       H2M =  HH2 + RLN*VREW(3) - H2REF !H2 DISTANCE REFER TO WAVE REFERNCE COORDINATE
! C     ----------------------------------------------------
      
      H1M = VH1(1)*XMID + VH1(2)*YMID + VH1(3)*ZMID 
      HGM = VGV(1)*XMID + VGV(2)*YMID + VGV(3)*ZMID 
      H2M = VH2(1)*XMID + VH2(2)*YMID + VH2(3)*ZMID

      	
      
C     FIND (OMEGA,WAVENUMBER(RK),RAMDA,A,RATIO) 
      !========================
      ! Call Wave length
      !========================
      IF (IWAVE == 1 . OR . IWAVE == 6 . OR . IWAVE == 7 ) THEN
          ORDER=0.0D0 ! PROTECT ERROR FROM STREAM FUNCTION
          OMEGA = 2.0*PI/PERIOD      
          CALL NEWTON_RAPHSON(WDEPTH,GRAV,OMEGA,RK)
          RAMDA = 2.0*PI/RK    
          IF (RK.LE.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,30)
          WRITE (*,31)
30        FORMAT ('************** OFFSHORE WARNING MESSAGE ( AIRY WAVE THEORY FOR SHELL ) **************')
31        FORMAT ('- PLEASE CHECK WAVE PERIOD AND THEORY CONDITION ( WAVE NUMBER < 0 )')
          STOP
          ENDIF
          AVAL = 1.0D0
          RATIO = 1.0D0                
      ELSEIF (IWAVE == 2) THEN
          ORDER=0.0D0 ! PROTECT ERROR FROM STREAM FUNCTION      
          !!!-CALL STOKES_WAVELENGTH(GRAV,WDEPTH,PERIOD,WVHIGHT,RK,RAMDA)
          CALL STOKES_WAVELENGTH_Time_Modify(RK,RAMDA)
          
          RATIO = WDEPTH/RAMDA(LCS)
          
          
          
          !	GET ALL PARAMETERS    
          CALL PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55)
          CALL PARAMETER_F(RATIO,F22,F24,F33,F35,F44,F55)
          CALL PARAMETER_C(RATIO,C1,C2,C3,C4)
          
          CALL STOKE_COEFFICIENT(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
          
          ! --- MODIFINE ON 22-08-2012 ----
          AVAL=RK
          RK=(2*(RK+((RK**3)*F33)+((RK**5)*(F35+F55))))/WVHIGHT 
          IF (RK.LE.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,32)
          WRITE (*,33)
32        FORMAT ('************** OFFSHORE WARNING MESSAGE ( STOKE FIFTH ORDER THEORY FOR SHELL ) **************')
33        FORMAT ('- PLEASE CHECK WAVE PERIOD AND THEORY CONDITION ( WAVE NUMBER < 0 )')
          STOP
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
      ! ****** NO NEED ******                 
      ENDIF    
   
              
C	---------------------------------------------------
C	WATER SURFACE LEVEL
C     MAXIMUM WATER SURFACE, SO APPLY ON THE STRUCTURE    
      ! IWAVE  = 1 >> AIRY WAVE THEORY
      ! IWAVE  = 2 >> STOKE'S FIFTH ORDER THEORY
      ! IWAVE  = 3 >> STREAM FUNCTION WAVE THEORY 
      ! IWAVE  = 4 >> CNOIDAL WAVE THEORY
      ! IWAVE  = 5 >> SOLITATY WAVE THEORY!  
              
      !==============================
      ! Call WATER SURFACE LEVELTION
      !==============================
      IF (IORRE.EQ.1) THEN 
      IF (IWAVE == 1) THEN 
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
    	
    	ELSEIF (IWAVE.EQ.4)THEN   
    	call Cnoidal_Wave_Crest  (WVHIGHT,WDEPTH,GRAV,PERIOD,TIME,Wave_Crest)
    	call Cnoidal_wave_Length (WVHIGHT,WDEPTH,GRAV,PERIOD,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
    	RK=wavenumber
    	RAMDA=ALAMDA
    	WNU = Wave_Crest    
    		
    	ELSEIF (IWAVE.EQ.5)THEN
    	Call Solitary_Wave (WVHIGHT,WDEPTH,GRAV,TIME,0d0,z,u,v,au,av,Surface_Elevation,q)
    	WNU = Surface_Elevation
      
      ELSEIF (IWAVE.EQ.6)THEN ! Diffrection
          
      WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
 
!    	WRITE (*,75)
!    	WRITE (*,76)
!    	WRITE (*,77)
!    	WRITE (*,76)
!75    FORMAT ('')
!76    FORMAT ('**************************************')
!77    FORMAT ('SOLITATY WAVE THEORY ON DEVELOPMENT') 
!      STOP     
      ENDIF     
C	---------------------------------------------------
      ELSEIF (IORRE.EQ.2) THEN
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
C	WATER SURFACE PROFILE MAX
      YL = WDEPTH + WNU 
C     WATER SURFACE PROFILE MIN
      YLMIN = WDEPTH - WNU     
      
      IF (WNU.GE.0.0D0) YLMIN = WDEPTH - WNU
      IF (WNU.LT.0.0D0) YLMIN = YL
      
      IF (IWAVE.EQ.1 .OR.IWAVE.EQ.6 )THEN
      YLMIN = WDEPTH - WNU
      ELSEIF (IWAVE.EQ.2.OR.IWAVE.EQ.3.OR.IWAVE.EQ.4.OR.IWAVE.EQ.5)THEN
      YLMIN = (WDEPTH + WNU) - WVHIGHT
      ENDIF
      
      
      
      IF (NFUNCTION.LE.3.0D0)THEN ! PREVENT ERROR FUNCTION FOR USER DEFINED
      ! -------- CALCULATE CD AND CM FOR AUTOMATIC PARTS ( BASE ON DNV )----------
      ! --------- DNV-OS-J101 JULY 2011 PAGE 74-75 - SEC.4  ---------
        ! NFUNCTION = 1 >> AUTOMATIC CACULATION ( DNV )
        ! NFUNCTION = 2 >> AUTOMATIC CACULATION ( API - ROUGH )
        ! NFUNCTION = 3 >> AUTOMATIC CACULATION ( API - SMOOTH )
        ! NFUNCTION = 4 >> USER DEFINED
        IF (NFUNCTION.EQ.1.0D0)THEN
           IF ((ROUGH/DIAM1).LT.0.0001)THEN
           CD  =  0.65 ! SMOOTH SURFACE
           ELSEIF ((ROUGH/DIAM1).GT.0.0001.AND.(ROUGH/DIAM1).LT.0.01)THEN
           CD  =  (29D0+4D0*LOG10(ROUGH/DIAM1))/20D0 
           ELSEIF ((ROUGH/DIAM1).GT.0.01)THEN
           CD  =  1.05 ! ROUGH SURFACE
           ENDIF
           
           !==============================
           ! FOR CD AND CM CALCULATION
           !==============================
           IF (IWAVE.EQ.1.0D0 .OR. IWAVE.EQ.6.0D0 )THEN ! AIRY WAVE THEORY
           ! GENERATE VELOCITY AT SURFACE 
           ! SET X POSITION =0.0, Y POSITION = WATER DEPTH, AND TIME =0.0
           UXMAX = 0.50D0*OMEGA*WVHIGHT*COSH(RK*WDEPTH)/SINH(RK*WDEPTH)*COS(RK*(0) - OMEGA*(TIME)) 
           UYMAX = 0.50D0*OMEGA*WVHIGHT*SINH(RK*WDEPTH)/SINH(RK*WDEPTH)*SIN(RK*(0) - OMEGA*(TIME)) 
           ! AKC = KEULEGAN CARPENTER SEE DNV-OS-J101 JULY 2011 PAGE 74-75 
           AKC   = UXMAX*PERIOD/DIAM1 
              IF (AKC.LT.3.0D0)THEN
              CM   =  2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  =  2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  =  1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM=CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM=CM2
                ENDIF
              ENDIF
           ! FOR WRITE STREAM FUNCTION FILE SEE XFINAS.FOR
           NSTREAMFUNCTION=0
           ELSEIF (IWAVE.EQ.2.0D0)THEN ! STOKE'S WAVE THEORY
           ! GENERATE VELOCITY AT SURFACE
           ! SET X POSITION =0.0, Y POSITION = WATER DEPTH, AND TIME =0.0
           CALL PARAMETER_G(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55)     
           CALL STOKE_COEFFICIENT(RATIO,G11,G13,G15,G22,G24,G33,G35,G44,G55,F22,F24,F33,F35,F44,F55,
     1                             C1,C2,C3,C4)
C	     -----------------------------------
           G1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
           G2 = 2.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
           G3 = 3.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
           G4 = 4.0D0*(AVAL**4.0)*G44
           G5 = 5.0D0*(AVAL**5.0)*G55
C	     -----------------------------------
C	     DRAG FORCE        
           UXMAX  =  (G1*COSH(1.0D0*RK*WDEPTH)/SINH(1.0D0*RK*WDEPTH)*COS(1.0D0*(RK*(0)
     1     - OMEGA*(TIME)))
     1     + G2*COSH(2.0D0*RK*WDEPTH)/SINH(2.0D0*RK*WDEPTH)*COS(2.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G3*COSH(3.0D0*RK*WDEPTH)/SINH(3.0D0*RK*WDEPTH)*COS(3.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G4*COSH(4.0D0*RK*WDEPTH)/SINH(4.0D0*RK*WDEPTH)*COS(4.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G5*COSH(5.0D0*RK*WDEPTH)/SINH(5.0D0*RK*WDEPTH)*COS(5.0D0*(RK*(0) 
     1     - OMEGA*(TIME))))
     1      *(OMEGA/RK)
           UYMAX  =  (G1*SINH(1.0D0*RK*Y)/SINH(1.0D0*RK*WDEPTH)*SIN(1.0D0*(RK*(0) 
     1     - OMEGA*(TIME))) 
     1     + G2*SINH(2.0D0*RK*WDEPTH)/SINH(2.0D0*RK*WDEPTH)*SIN(2.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G3*SINH(3.0D0*RK*WDEPTH)/SINH(3.0D0*RK*WDEPTH)*SIN(3.0D0*(RK*(0) 
     1     - OMEGA*(TIME))) 
     1     + G4*SINH(4.0D0*RK*WDEPTH)/SINH(4.0D0*RK*WDEPTH)*SIN(4.0D0*(RK*(0) 
     1     - OMEGA*(TIME)))
     1     + G5*SINH(5.0D0*RK*WDEPTH)/SINH(5.0D0*RK*WDEPTH)*SIN(5.0D0*(RK*(0) 
     1     - OMEGA*(TIME))))
     1      *(OMEGA/RK)
           
      !==================================
      !    MODIFY BY CHANA 23 FEB 2016
      !==================================
      
      !WVHIGHT = H
      !WDEPTH  = HW
      !GRAV = G
           
       Y = WDEPTH
       X =  0.0d0
       TT = TIME
       
      RK=(2*(AVAL+((AVAL**3)*F33)+((AVAL**5)*(F35+F55))))/WVHIGHT 
      
      CWAVE = SQRT(GRAV/RK*(1.0D0 + C1*(AVAL**2.0) + C2*AVAL**4.0)*TANH(RK*WDEPTH))
      
      chana = (OMEGA/RK)
      
      FG1 = 1.0D0*(AVAL*G11 + (AVAL**3.0)*G13 + (AVAL**5.0)*G15)
      FG2 = 1.0D0*((AVAL**2.0)*G22 + (AVAL**4.0)*G24)
      FG3 = 1.0D0*((AVAL**3.0)*G33 + (AVAL**5.0)*G35)
      FG4 = 1.0D0*(AVAL**4.0)*G44
      FG5 = 1.0D0*(AVAL**5.0)*G55
      
      TT = TIME 
      
      UXMAX = CWAVE *( 1.0D0*FG1*COSH(1.0D0*RK*Y)*COS(1.0D0*(RK*X - OMEGA*TIME))
     1               + 2.0D0*FG2*COSH(2.0D0*RK*Y)*COS(2.0D0*(RK*X - OMEGA*TIME))
     1               + 3.0D0*FG3*COSH(3.0D0*RK*Y)*COS(3.0D0*(RK*X - OMEGA*TIME))
     1               + 4.0D0*FG4*COSH(4.0D0*RK*Y)*COS(4.0D0*(RK*X - OMEGA*TIME))
     1               + 5.0D0*FG5*COSH(5.0D0*RK*Y)*COS(5.0D0*(RK*X - OMEGA*TIME)))
      
      UYMAX = CWAVE*( 1.0D0*FG1*SINH(1.0D0*RK*Y)*SIN(1.0D0*(RK*X - OMEGA*TIME))
     1              + 2.0D0*FG2*SINH(2.0D0*RK*Y)*SIN(2.0D0*(RK*X - OMEGA*TIME))
     1              + 3.0D0*FG3*SINH(3.0D0*RK*Y)*SIN(3.0D0*(RK*X - OMEGA*TIME))
     1              + 4.0D0*FG4*SINH(4.0D0*RK*Y)*SIN(4.0D0*(RK*X - OMEGA*TIME))
     1              + 5.0D0*FG5*SINH(5.0D0*RK*Y)*SIN(5.0D0*(RK*X - OMEGA*TIME)))
      
      CHANA = 3
      
           AKC  =  UXMAX*PERIOD/DIAM1
              IF (AKC.LT.3.0D0)THEN
              CM  =  2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  =  2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  =  1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM =  CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM =  CM2
                ENDIF
              ENDIF 
           ! FOR WRITE STREAM FUNCTION FILE SEE XFINAS.FOR
           NSTREAMFUNCTION=0
           ELSEIF (IWAVE.EQ.3.0D0)THEN
           ! GENERATE VELOCITY AT SURFACE
           ! SET X POSITION =0.0, Y POSITION = WATER DEPTH, AND TIME =0.0
            CALL VELOC (0,WDEPTH,UXMAX,UYMAX,DUMAX,DVMAX,TIME,ARAMDA)
           AKC=UXMAX*PERIOD/DIAM1
           ! FOR WRITE STREAM FUNCTION FILE SEE XFINAS.FOR
           NSTREAMFUNCTION=1
              IF (AKC.LT.3.0D0)THEN
              CM   =  2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  =  2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  =  1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM =  CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM =  CM2
                ENDIF
              ENDIF 
           ELSEIF (IWAVE.EQ.4.0D0)THEN ! CANOIDAL WAVE 
           ZZETAA = wavenumber * X - 2d0 * AKKK/PERIOD * TIME
           CALL Jacobian_elliptic_function (ZZETAA,Ammm,CN,CN2,SSN,ddn)
           CALL Cnoidal_Wave_Kinematic  (WVHIGHT,WDEPTH,GRAV,PERIOD,TIME,WDEPTH,wavenumber,AKKK,ALAMDA,Ammm,AEEE,CN,SSN,ddn,ZZETAA
     1                                 ,UXMAX,UYMAX,DUDT,DVDTM)
              AKC = UXMAX*PERIOD/DIAM1
              IF (AKC.LT.3.0D0)THEN
              CM   = 2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  = 2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  = 1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM = CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM = CM2
                ENDIF
              ENDIF 
     
           ELSEIF (IWAVE.EQ.5.0D0)THEN
           ! SOLITARY WAVE    
           Z = WDEPTH
           Call Solitary_Wave (WVHIGHT,WDEPTH,GRAV,TIME,0d0,Z,UXMAX,UYMAX,au,av,WNU,q)
              
              AKC = UXMAX*PERIOD/DIAM1
              IF (AKC.LT.3.0D0)THEN
              CM   = 2.0D0
              ELSEIF (AKC.GT.3.0D0)THEN
              CM1  = 2.0D0-(0.044D0*(AKC-3.0D0))
              CM2  = 1.6D0-(CD-0.5D0)
                IF (CM1.GT.CM2)THEN
                CM = CM1
                ELSEIF (CM2.GT.CM1)THEN
                CM = CM2
                ENDIF
              ENDIF 
           
           ENDIF
        ELSEIF (NFUNCTION.EQ.2.0D0.OR.NFUNCTION.EQ.3.0D0)THEN ! BASE ON API
           ! API 2A-WSD(RP 2A-WSD) PAGE 15
           IF (NFUNCTION.EQ.2.0D0)THEN ! ROUGH SURFACE
           CD  = 1.05D0
           CM  = 1.20D0
           ELSEIF (NFUNCTION.EQ.3.0D0)THEN ! SMOOTH SURFACE
           CD  = 0.65D0
           CM  = 1.20D0
           ENDIF
        ENDIF ! (ENDIF : END AUTOMATIC GENERATE CD,CM FUNCTION >> 1,DNV 2,API  )
      ENDIF   ! (ENDIF : PREVENT ERROR FUNCTION FOR USER DEFINED )
C      ACM(LCS)=CM
C      ACD(LCS)=CD
    ! ----------------------------------------------------------------------------------------------------------------------
    
    
C     MODIFY CM WITH NORMINAL DIAMETER AND SHAPE FACTOR        
        CM = CM * DIAM
      

C     ===================================================	
C     ===================================================	

      CALL CLEARA (RL,3*NNO)

	MMGR = 3
	MMGS = 3
	IF(NNO.EQ.3.OR.NNO.EQ.6) THEN
	MMGR = 3
	MMGS = 1   !SHELL 3 NODE
	ENDIF

C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
      SELECTCASE(NNO)
      CASE(3)!SHELL 3 NODE
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNO)
      CASE(6)!SHELL 3 NODE ONATE
      H = 0.0D0
      P = 0.0D0
      NNOM = 3
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNOM)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNOM)
	CASE DEFAULT
      IF(NNO.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX(1,IELE),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,IELE),NNO)            !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,IELE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	ENDSELECT
C     --------------------------------------------------------
      XYZCP(1:3) = 0.0D0
	RLEV = 0.0D0
	DO INO = 1,NNO
	  JJ = 3*(INO-1) + NGRAV
	  RLEV = RLEV + H(INO)*XYZ(JJ,IELE)
	  DO I = 1,3
	      JJ = 3*(INO-1) + I
	      XYZCP(I) = XYZCP(I) + H(INO)*XYZ(JJ,IELE)
	  ENDDO
      ENDDO
      
C     GAUSS COORDINATE IN WAVE DIRECTION SYSTEM	
      H1R = VH1(1)*XYZCP(1) + VH1(2)*XYZCP(2) + VH1(3)*XYZCP(3)
      HGR = VGV(1)*XYZCP(1) + VGV(2)*XYZCP(2) + VGV(3)*XYZCP(3)
      H2R = VH2(1)*XYZCP(1) + VH2(2)*XYZCP(2) + VH2(3)*XYZCP(3)

C     GAUSS COORDINATE RELATIVE TO WAVE REFERENCE POINT (IN WAVE DIRECTION SYSTEM)
      H1R =  H1R - H1REF
      HGM =  HGR - HGREF
      H2M =  H2R - H2REF

      DIFFRACT_X = H1R
      DIFFRACT_Y = H2M
      DIFFRACT_GRAVITY = HGM
      
      
C     --------------------------    
      KFLACAL = 0
      
      TT(1:3) = 0.0D0
      
      VNOLG(1:3) = VH1(1:3)
      
      IF (LWCASE(2).EQ.1.0D0.AND.LWCASE(1).EQ.0.0)THEN
      YL = WDEPTH
      ENDIF       
              
      DO 407 LCASE = 1,7
      ICASE = LWCASE(LCASE)
      IF(ICASE.EQ.0) GOTO 407
      
      SELECTCASE(LCASE)
      CASE(1,2,3) !WAVE LOAD - CURRENT LOAD - WATER TIDE LOAD    
      
      IF(KFLACAL.EQ.1) GOTO 407  !IF CURRENT&TIDAL LOAD ALREADY CALCULATE TOGETHER WITH WAVE LOAD, THEN JUMP
        
      IF(LWCASE(2).EQ.0) VWIND0 = 0.0D0
      IF(LWCASE(3).EQ.0) VTIDE  = 0.0D0
      
          WLEV = RLEV - (SEABED + WDEPTH)
          IF(WLEV.GT.0.0D0) GOTO 407
          !IF(WLEV.GT.YL) GOTO 407
          
      CALL WAVE_PRESSURE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                   VTIDE,VWIND0,WH1,WGV,WH2,VNOL,NCURRENT,H0,VCURRENTP,POWERLAW,LC,VCURRENTL,PERIOD,LWCASE,CS,
     1                   WKF,CBF,AP,SP,
     1                   WFC,YLMIN,IORRE,IRWAVE)
      
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 
      
C     TRANSFORM WAVE FORCE VECTOR TO GLOBAL   VNOL  TO VNOLG      
      VNOLG(1) = VH1(1)*VNOL(1) + VGV(1)*VNOL(2) + VH1(1)*VNOL(3)
      VNOLG(2) = VH1(2)*VNOL(1) + VGV(2)*VNOL(2) + VH1(2)*VNOL(3)
      VNOLG(3) = VH1(3)*VNOL(1) + VGV(3)*VNOL(2) + VH1(3)*VNOL(3) 

      FACA = ABS(VT(1)*VNOLG(1)+VT(2)*VNOLG(2)+VT(3)*VNOLG(3))  !PROJECTION AREA IN FORCE DIRECTION
      WFX = 0.5*WFX*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA
      WFY = 0.5*WFY*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
      WFZ = 0.5*WFZ*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
            
      KFLACAL = 1
!    ---------------------- REMOVE BY TOEY 2019 ----------------------  
!      CASE(4) !WATER LEVEL LOAD   
!          WLEV = RLEV - SEABED
!          IF(WLEV.LT.0.0D0) GOTO 407
!          IF(WLEV.GT.PEAKWLEV) GOTO 407
!      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
!     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
!     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF)
!      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
!      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
!      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 
!    ---------------------- REMOVE BY TOEY 2019 ----------------------  
      
      CASE(5,6) !WAVE BREAKING LOAD  PLUNGING&SURGING
          WLEV = RLEV - (SEABED + WDEPTH)
          IF(WLEV.GT.0.0D0) GOTO 407
      CALL WAVE_FORCE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGMWAVEINP,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                VTIDE,VWIND0,WH1,WGV,WH2,NCURRENT,H0,VCURRENTP,POWERLAW,VCURRENTL,PERIOD,LWCASE,CS,
     1                WKF,CBF,AP,SP,
     1                WFC,YLMIN,
     1                VELO,ACCE,COA,IORRE,IRWAVE)
      
!    ---------------------- REMOVE BY TOEY 2019 ----------------------  
!      CASE(5) !WAVE BREAKING LOAD  PLUNGING 
!          WLEV = RLEV - SEABED
!          IF(WLEV.LT.0.0D0) GOTO 407
!          IF(WLEV.GT.PEAKWLEV) GOTO 407
!      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
!     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
!     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
!      WFX = VGV(1)*WFF ; WFY = VGV(2)*WFF ; WFZ = VGV(3)*WFF ;
!      CASE(6) !WAVE BREAKING LOAD  SURING
!          WLEV = RLEV - SEABED
!          IF(WLEV.LT.0.0D0) GOTO 407
!          IF(WLEV.GT.PEAKWLEV) GOTO 407
!      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
!     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
!     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
!    ---------------------- REMOVE BY TOEY 2019 ----------------------  
      
      WFX = VH1(1)*WFF ; WFY = VH1(2)*WFF ; WFZ = VH1(3)*WFF ;
      CASE(7) !WIND LOAD 
          WLEV = RLEV - SEABED - WDEPTH
          IF(WLEV.LT.0.0D0) GOTO 407 
          ! -------------------------------------------------------------------------------------------------------
          ! GENERATE WIND LOAD ( DNV,API,IEC64100-1 )
          ! NWINDO = 1 ; LOGARITHMIC PROFILE   (DNV) 
          ! NWINDO = 2 ; POWER LAW PROFILE     (DNV)
          ! NWINDO = 3 ; WIND PROFILE AND GUST (API)
          ! NWINDO = 4 ; IEC 614000-1
          ! NWINDO = 5 ; USER DEFINED
            IF (NWINDO.EQ.1)THEN ! DNV > LOGARITHMIC PROFILE
            UZT=UH*(LOG(WLEV/Z0)/LOG(RHIGH/Z0))
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
            
            ELSEIF (NWINDO.EQ.2.OR.NWINDO.EQ.5)THEN ! DNV > POWER LAW PROFILE NWINDO >> 2, USER DEFIND NWINDO >> 5 
                 IF (NWINDO.EQ.5)THEN ! SAME EQUATION FROM USER DEFINED AND POWER LAW IN DNV
                 UH=UHD               ! UH  = MEAN WIND SPEED ON DNV 
                 ENDIF                ! UHD = MEAN WIND SPEED ON USER DEFINED
            CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     2                        CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
     
            ELSEIF (NWINDO.EQ.3)THEN ! API > WIND PROFILE AND GUST
            C=0.0573D0*SQRT(1D0+0.0457D0*UHAPI)
            UZ=UHAPI*(1+(C*LOG(WLEV/32.8D0)))
            ZI=0.06D0*(1+0.0131D0*UHAPI)*((WLEV/32.8D0)**(-0.22))
            UZT=UZ*(1-0.41D0*ZI*LOG(AVERAGE/3600D0))
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
            
            ELSEIF (NWINDO.EQ.4)THEN ! IEC 614000-1
            UZT=UHD*((WLEV/RHIGH)**ALPHA)
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
            ENDIF
         ! ----------------------------------------------------------------------------------------------------------
         ! CSA = SHAPE COEFFICIENT FOR WIND LOAD
      
      WFX = VWIND(1)*WFF * SHFAC * GUSTFACTOR
      WFY = VWIND(2)*WFF * SHFAC * GUSTFACTOR
      WFZ = VWIND(3)*WFF * SHFAC * GUSTFACTOR ;
      ENDSELECT
      
	TT(1) = TT(1) + WFX
	TT(2) = TT(2) + WFY
	TT(3) = TT(3) + WFZ
	
407   CONTINUE
C     -------------------------- 

      

	DO 502 INO=1,NNO
	FAC = H(INO)*WT*DET
      DO 481 I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 502
 481  RL(I,INO) = RL(I,INO) + TT(ICO)*FAC
 502  CONTINUE

 800  CONTINUE


C     MUST BE ADDED LATER  FOR SHELL FIXEND FORCED
C	CALL FIXSHE (RL,KEG,IELE,NNO,3)

C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,4)
C
 700  DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + RL(INF,INO) !VARY LOAD
      IF (IEQ.NE.0.AND.ILCC.GT.0) RRC(IEQ) = RRC(IEQ) + RL(INF,INO) !CONSTANT LOAD
 600  CONTINUE

      

 900  CONTINUE

      
      IF(KOFFL.NE.0) THEN
        CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
        CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','WRIT') !STORE CONSTANT LOAD
      ENDIF      

7500  CONTINUE

	DEALLOCATE(RRV,RRC,RVAL,IVAL)
C
1000  RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SHEOFFLSPEC (PROPG,IGSET,XYZ,NODEX,LM,
	1					MGP,MXY,MEX,MEF,NLOAD,IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CONVERTES OFFSHORE LOAD TO EQUIVALENT NODAL LOADS
C	-----------------------------------
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     IGSET(NELE)      = GEOMETRIC SET NUMBER
C     XYZ(MXY,NELE)    = NODAL COORDINATES FOR ELEMENTS
C     NODEX(NEX,NELE)  = LOCATIONS OF EXCESSIVE NODES
C     LM(NEF,NELE)     = EQUATION NUMBERS FOR ELEMENT D.O.F.
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================
C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

C     ============================================================================
C     OFFSHORE PART COMMON 
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION
	COMMON /WARNING/ WARNING,RAMDA(100),RK
      COMMON / EIGVPED / EIGVFREQ(1000),EIGVPERI(1000)
      COMMON / STOREMODE / RVECT(100000,100),NMOD
C     ============================================================================

      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION R(NEQ),H(9),P(2,9),XJI(4),FA(4),G(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,9),IGPOS(9),IGIDM(NELE)

C     FOR OFFSHORE LOAD
      DIMENSION VWIND(3),VH1(3),VH2(3),VGV(3),LWCASE(10),TT(3),XYZCP(3)
      
C     FOR OFFSHORE LOAD -- PRAMIN 
      DIMENSION FLIST(5),VREF(3),VREW(3),VNOL(3),VNOLG(3),VCURRENTL(5)
      
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:)
	      

      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)


	IF(NLOAD.EQ.0) GOTO 1000	
      	
	PI = 3.141592654

	ALLOCATE(RVAL(25,NLOAD),IVAL(20,NLOAD))
	
C     READ INPUT DATA AND STORE ON TEMPORARY VARIABLE
	DO ILOAD=1,NLOAD
        READ(ITI,*) IELE,NUMCP,NFUNCTION,ROUGH,CD,CM,CSA,DIAM1,SHFAC,NGROWTH,GROWTH,OFFSELECT_X,LWCASE(1:4),ILCN,ILCC 
        ! ----------------------
        LWCASE(5)    = LWCASE(2)
        LWCASE(6)    = LWCASE(3)
        LWCASE(7)    = LWCASE(4)
        ! ---------------------
        
        IVAL(1 ,ILOAD) = IELE
        IVAL(2 ,ILOAD) = NUMCP
        IVAL(3 ,ILOAD) = LWCASE(1)
        IVAL(4 ,ILOAD) = LWCASE(2)
        IVAL(5 ,ILOAD) = LWCASE(3)
        IVAL(6 ,ILOAD) = LWCASE(4)
        IVAL(7 ,ILOAD) = LWCASE(5)
        IVAL(8 ,ILOAD) = LWCASE(6)
        IVAL(9 ,ILOAD) = LWCASE(7)
        IVAL(10,ILOAD) = ILCN
        IVAL(11,ILOAD) = ILCC
        
        RVAL(1 ,ILOAD) = CD
        RVAL(2 ,ILOAD) = CM
        DIAM2=DIAM1+(GROWTH*2D0)
        RVAL(12 ,ILOAD) = DIAM2   !NORMINAL DIAMETER
        RVAL(13 ,ILOAD) = SHFAC  !SHAPE FACTOR = 0.25*PI for circular section
        RVAL(21,ILOAD)  = OFFSELECT_X  ! CASE SELECT
        RVAL(22,ILOAD)  = NFUNCTION
        RVAL(23,ILOAD)  = ROUGH
        RVAL(24,ILOAD)  = DIAM1
        RVAL(25,ILOAD)  = CSA 
        ! -------------------------------  FOR ERROR MESSAGE MARINE GROWTH ------------------------------- 
             IF (NGROWTH.EQ.0.0D0)THEN
                IF (GROWTH.NE.0.0D0)THEN
                   WRITE (*,59)
                   WRITE (*,60)
                   WRITE (*,61)
                   WRITE (*,59)
                   WRITE (*,62)
                   WRITE (*,63)
                   WRITE (*,64)
59                 FORMAT ('')                   
60                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( FRAME ELEMENT )----------------')
61                 FORMAT ('- INSIDE MARINE GROWTH FUNCTION HAVE SOME VALUE. II WILL EFFECT ON THE RESULT')
62                 FORMAT ('**** COMMENT *****')
63                 FORMAT (' CALCULATE     = SELECT MARINE GROWTH BLOCK')
64                 FORMAT (' NOT CALAULATE = UNSELECT MARINE GROWTH BLOCK')
                  STOPFUNCTION=1D0 
                  GOTO 101
                ELSEIF (GROWTH.EQ.0.0D0)THEN
                ! CONTINUE 
                ENDIF
             ELSEIF (NGROWTH.NE.0.0D0)THEN
                IF (GROWTH.GT.0.0D0)THEN
                ! CONTINUE
                ELSEIF (GROWTH.LE.0.0D0)THEN
                   WRITE (*,65)
                   WRITE (*,66)
                   WRITE (*,67)
                   WRITE (*,65)
                   WRITE (*,68)
                   WRITE (*,69)
                   WRITE (*,70)
65                 FORMAT ('')
66                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( FRAME ELEMENT )----------------')
67                 FORMAT ('- MARINE GROWTH MUST BE GRATER THAN ZERO ')
68                 FORMAT ('**** COMMENT *****')
69                 FORMAT (' CALCULATE     = SELECT MARINE GROWTH BLOCK')
70                 FORMAT (' NOT CALAULATE = UNSELECT MARINE GROWTH BLOCK')
                   ! STOPFUNCTION USE FOR STOP PROGRAM WHEN ERROR MESSAGE APPEAR
                   STOPFUNCTION=1D0 
                   ! GOTO 101 : OUT READ INPUT FOR STOP THE PROGRAM
                   GOTO 101
                ENDIF
             ENDIF
         ! -------------------------------------------------------------------------------------------------- 
	ENDDO
		
	STOPFUNCTION=0.0D0 ! SET CONDITION PROTECT ERROR EFFECT
	! ---- WARNING MESSAGE ---
101   IF (LWCASE(1).EQ.0.AND.LWCASE(2).EQ.0.AND.LWCASE(3).EQ.0.AND.LWCASE(4).EQ.0.AND.LWCASE(5).EQ.0.AND.LWCASE(6).EQ.0
     1    .AND.LWCASE(7).EQ.0.0D0)THEN
      WRITE (*,8)
      WRITE (*,34)
      WRITE (*,35)
34    FORMAT ('---------------- OFFSHORE WARNING MESSAGE ----------------')
35    FORMAT ('- PLEASE CHECK OFFSHORE LOAD TYPES')
      STOP
      ENDIF
      IF (NFUNCTION.EQ.3)THEN
      IF (CD.LE.0.0D0.OR.CM.LE.0.0D0.OR.DIAM2.LE.0.0D0)THEN
      WRITE (*,8)
      WRITE (*,36)
36    FORMAT ('---------------- OFFSHORE ERROR MESSAGE ----------------')
          IF (CD.LE.0.0D0)THEN
          WRITE (*,37)
37        FORMAT ('- DRAG COEFFICIENT MUST BE GRATER THAN ZERO')
          ENDIF
          IF (CM.LE.0.0D0)THEN
          WRITE (*,38)
38        FORMAT ('- INERTIA COEFFICIENT MUST BE GRATER THAN ZERO')
          ENDIF
          IF (DIAM2.LE.0.0D0)THEN
          WRITE (*,39)
39        FORMAT ('- NORMINAL DIMAETER MUST BE GRATER THAN ZERO')
          ENDIF
      STOPFUNCTION=1D0    
      ENDIF
      ENDIF
	
	
	IF(KOFFL.EQ.0) THEN
	  DEALLOCATE(RVAL,IVAL)
	  GOTO 1000 !IF NO OFFSHORE LOAd ANALYSIS NO NEED TO DO THE CALCULATION...JUST READ THE INPUT DATA
	ENDIF
	
      
	ALLOCATE(RRV(NEQ),RRC(NEQ))

      IF(KOFFL.NE.0) THEN !CALL ALL OFFSHORE LOAD PARAMETERS
    !  CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,PERIOD,GRAV,RHOW,RHOA,
    ! 1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
    ! 1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED)
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      ENDIF  
      	
C     -----------------------------------
C     -----------------------------------

      DO 7500 ILCAS = 1,LCS !LOOP OVER LOAD CASE NUBER
      
      DO 7500 ITIME = 1,NTIME !LOOP OVER OFFSHORE TIME STEP

      RRV(1:NEQ) = 0.0D0 !INITIALIZE VARY OFFSHORE LOAD
      RRC(1:NEQ) = 0.0D0 !INITIALIZE CONSTANT OFFSHORE LOAD
      
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')
      CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','READ')
            
      DO 7501 ISPEC = 1,NMOD
      SNA = EIGVFREQ(ISPEC)
            
      DO 900  ILOAD=1,NLOAD

C     ------------------------------
C     CALL INPUT DATA FROM BACKUP
        IELE        = IVAL(1 ,ILOAD)
        NUMCP       = IVAL(2 ,ILOAD)
        LWCASE(1)   = IVAL(3 ,ILOAD)
        LWCASE(2)   = IVAL(4 ,ILOAD)
        LWCASE(3)   = IVAL(5 ,ILOAD)
        LWCASE(4)   = IVAL(6 ,ILOAD)
        LWCASE(5)   = IVAL(7 ,ILOAD)
        LWCASE(6)   = IVAL(8 ,ILOAD)
        LWCASE(7)   = IVAL(9 ,ILOAD)
        ILCN        = IVAL(10,ILOAD)
        ILCC        = IVAL(11,ILOAD)
        
        NFUNCTION   = RVAL(22,ILOAD) 
        ROUGH       = RVAL(23,ILOAD)
        DIAM        = RVAL(12,ILOAD)
        DIAM1       = RVAL(24,ILOAD)
        SHFAC       = RVAL(13,ILOAD)
        CSA         = RVAL(25,ILOAD)
        offselect   = RVAL(21,ILOAD) 
        
        IF (NFUNCTION.EQ.2.0D0)THEN
        CD          = RVAL(1 ,ILOAD)
        CM          = RVAL(2 ,ILOAD)
        ELSEIF (NFUNCTION.EQ.1.0)THEN
        CD          = 0D0 
        CM          = 0D0
        ENDIF
        
       CALL SELECTSPECTRUM (OFFSELECT,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
       IF (ISOLOP.EQ.1.AND.ICONTROLSPEC.EQ.1)THEN
       FREQ =  TIME
       ENDIF
       
      WARNING=GRAV    
C    ================================================================== 
8     FORMAT("")
      ! --- STOP FUNCTION ---
      IF (STOPFUNCTION.EQ.1D0)THEN
      STOP
      ENDIF
      ! --- END STOP FUNCTION -- 
C     ------------------------------
      IF(ILCN.NE.ILCAS.AND.ILCC.NE.ILCAS) GOTO 900
      
      
	
	CALL ELEREODER(IGIDM,NELE,IELE)

      ISET = IGSET(IELE)
	THICK = PROPG(2,ISET)


C     ===================================================	
C     PRAMIN WAVE LOAD MODIFICATION
C     ===================================================

C	WAVE REFERNCE COORDINATE
      H1REF = H1POS
      H2REF = H2POS
      HGREF = SEABED		

C     ----------------------------------------------------
C     ELEMENT MID COORDINATES
	XMID = 0.0D0
      YMID = 0.0D0
      ZMID = 0.0D0
	DO INO = 1,NNO
	  DO I = 1,3
	      JX = 3*(INO-1) + 1
	      JY = 3*(INO-1) + 2
	      JZ = 3*(INO-1) + 3
	      XMID = XMID + XYZ(JX,IELE)
	      YMID = YMID + XYZ(JY,IELE)
	      ZMID = ZMID + XYZ(JZ,IELE)
	  ENDDO
	ENDDO
      XMID = XMID/NNO
      YMID = YMID/NNO
      ZMID = ZMID/NNO

C     ELEMENT MID COORDINATES IN WAVE DIRECTION SYSTEM
      HM1 = VH1(1)*XMID + VH1(2)*YMID + VH1(3)*ZMID 
      HMG = VGV(1)*XMID + VGV(2)*YMID + VGV(3)*ZMID 
      HM2 = VH2(1)*XMID + VH2(2)*YMID + VH2(3)*ZMID 
C     ----------------------------------------------------
      
      
C     ----------------------------------------------------
C     STRUCTURAL ALIGNMENT REFERENCE PROPERTIES
      CALL OFFREFAXIS(NUMCP,XXR,YYR,ZZR,VREF,'OREF')
      
C     STRUCTURAL ALIGNMENT REFERENCE POSITION IN WAVE DIRECTION SYSTEM 
      HH1 = VH1(1)*XXR + VH1(2)*YYR + VH1(3)*ZZR 
      HHG = VGV(1)*XXR + VGV(2)*YYR + VGV(3)*ZZR 
      HH2 = VH2(1)*XXR + VH2(2)*YYR + VH2(3)*ZZR 

C     STRUCTURAL ALIGNMENT REFERENCE VECTOR IN WAVE DIRECTION SYSTEM 
      VREW(1) = VH1(1)*VREF(1) + VH1(2)*VREF(2) + VH1(3)*VREF(3)
      VREW(2) = VGV(1)*VREF(1) + VGV(2)*VREF(2) + VGV(3)*VREF(3) 
      VREW(3) = VH2(1)*VREF(1) + VH2(2)*VREF(2) + VH2(3)*VREF(3) 
C     ----------------------------------------------------
      
      
C     ----------------------------------------------------
C     POSITION OF COORDINATES REFER TO WAVE REFERNCE COORDINATE (IN WAVE DIRECTION SYSTEM) CORESPONDING TO STRUCTURAL ALIGNMENT
      HIG = HMG - HHG   !VERTICAL DISTANCE FROM STRUCTURAL REFERENCE COORDINATE
      RLN = HIG/VREW(2) !LENGTH ALONG ALIGNMENT VECTOR
      H1M =  HH1 + RLN*VREW(1) - H1REF !H1 DISTANCE REFER TO WAVE REFERNCE COORDINATE
      HGM =  HHG + RLN*VREW(2) - HGREF !GRAVITY DISTANCE REFER TO WAVE REFERNCE COORDINATE
      H2M =  HH2 + RLN*VREW(3) - H2REF !H2 DISTANCE REFER TO WAVE REFERNCE COORDINATE
C     ----------------------------------------------------

      	
   ! SET IWAVE = 1
      IWAVE = 1
      OMEGA = 2.0*PI/PERIOD      
      CALL NEWTON_RAPHSON(WDEPTH,GRAV,OMEGA,RK)
      RAMDA(ILCAS) = 2.0*PI/RK    
      AVAL = 1.0D0
      RATIO = 1.0D0   
      
      WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*0D0)
      YL = WDEPTH + WNU
      YLMIN = WDEPTH - WNU
            
         CALL DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM1,PERIOD,WVHIGHT,RK
     1                          ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                          ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
    

C     MODIFY CM WITH NORMINAL DIAMETER AND SHAPE FACTOR        
        CM = CM*DIAM*SHFAC
      

C     ===================================================	
C     ===================================================	

      CALL CLEARA (RL,3*NNO)

	MMGR = 3
	MMGS = 3
	IF(NNO.EQ.3.OR.NNO.EQ.6) THEN
	MMGR = 3
	MMGS = 1   !SHELL 3 NODE
	ENDIF

C	-------------------------------------------------------
C     GAUSS POINT LOOP
C	-------------------------------------------------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
      SELECTCASE(NNO)
      CASE(3)!SHELL 3 NODE
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNO)
      CASE(6)!SHELL 3 NODE ONATE
      H = 0.0D0
      P = 0.0D0
      NNOM = 3
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNOM)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNOM)
	CASE DEFAULT
      IF(NNO.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX(1,IELE),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,IELE),NNO)            !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,IELE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	ENDSELECT
C     --------------------------------------------------------
      XYZCP(1:3) = 0.0D0
	RLEV = 0.0D0
	DO INO = 1,NNO
	  JJ = 3*(INO-1) + NGRAV
	  RLEV = RLEV + H(INO)*XYZ(JJ,IELE)
	  DO I = 1,3
	      JJ = 3*(INO-1) + I
	      XYZCP(I) = XYZCP(I) + H(INO)*XYZ(JJ,IELE)
	  ENDDO
	ENDDO

C     GAUSS COORDINATE IN WAVE DIRECTION SYSTEM	
      H1R = VH1(1)*XYZCP(1) + VH1(2)*XYZCP(2) + VH1(3)*XYZCP(3)
      HGR = VGV(1)*XYZCP(1) + VGV(2)*XYZCP(2) + VGV(3)*XYZCP(3)
      H2R = VH2(1)*XYZCP(1) + VH2(2)*XYZCP(2) + VH2(3)*XYZCP(3)

C     GAUSS COORDINATE RELATIVE TO WAVE REFERENCE POINT (IN WAVE DIRECTION SYSTEM)
      H1R =  H1R - H1REF
      HGM =  HGR - HGREF
      H2M =  H2R - H2REF

      
C     --------------------------    
      KFLACAL = 0
      
      TT(1:3) = 0.0D0
      
      VNOLG(1:3) = VH1(1:3)
      
      IF (LWCASE(2).EQ.1.0D0.AND.LWCASE(1).EQ.0.0)THEN
      YL = WDEPTH
      ENDIF       
              
      DO 407 LCASE = 1,7
      ICASE = LWCASE(LCASE)
      IF(ICASE.EQ.0) GOTO 407
      
      SELECTCASE(LCASE)
      CASE(1,2,3) !WAVE LOAD - CURRENT LOAD - WATER TIDE LOAD    
      
      IF(KFLACAL.EQ.1) GOTO 407  !IF CURRENT&TIDAL LOAD ALREADY CALCULATE TOGETHER WITH WAVE LOAD, THEN JUMP
        
      IF(LWCASE(2).EQ.0) VWIND0 = 0.0D0
      IF(LWCASE(3).EQ.0) VTIDE  = 0.0D0
      
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407
          IF(WLEV.GT.YL) GOTO 407
          
C      CALL WAVE_PRESSURE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
C     1                   VTIDE,VWIND0,WH1,WGV,WH2,VNOL,NCURRENT,H0,VCURRENTP,POWERLAW,LC,VCURRENTL,PERIOD,LWCASE,CS,
C     1                   WKF,CBF,AP,SP,
C     1                   WFC,YLMIN,IORRE,IRWAVE)

         IF (IRWAVE.EQ.1)THEN
         CALL  Pierson_Moskowitz_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                         WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.2)THEN
         CALL Jonswap_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                         WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.3)THEN
         CALL Ochi_SPectrum_Solid(WVH1,PER1,SHAPE1,WVH2,PER2,SHAPE2,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                            WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.4)THEN
         CALL  Bretschneider_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                         WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.5)THEN
         CALL TMA_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                         WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.6)THEN
         CALL User_Define_SPectrum_SOILD(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                         WH1,WGV,WH2,TAP,UCURRENT)
         ENDIF
         
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 
      
C     TRANSFORM WAVE FORCE VECTOR TO GLOBAL   VNOL  TO VNOLG      
      VNOLG(1) = VH1(1)*VNOL(1) + VGV(1)*VNOL(2) + VH1(1)*VNOL(3)
      VNOLG(2) = VH1(2)*VNOL(1) + VGV(2)*VNOL(2) + VH1(2)*VNOL(3)
      VNOLG(3) = VH1(3)*VNOL(1) + VGV(3)*VNOL(2) + VH1(3)*VNOL(3) 

      FACA = ABS(VT(1)*VNOLG(1)+VT(2)*VNOLG(2)+VT(3)*VNOLG(3))  !PROJECTION AREA IN FORCE DIRECTION
      WFX = 0.5*WFX*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA
      WFY = 0.5*WFY*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
      WFZ = 0.5*WFZ*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
            
      KFLACAL = 1
      
      CASE(4) !WATER LEVEL LOAD   
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407
          IF(WLEV.GT.PEAKWLEV) GOTO 407
      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF)
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 
      CASE(5) !WAVE BREAKING LOAD  PLUNGING 
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407
          IF(WLEV.GT.PEAKWLEV) GOTO 407
      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
      WFX = VGV(1)*WFF ; WFY = VGV(2)*WFF ; WFZ = VGV(3)*WFF ;
      CASE(6) !WAVE BREAKING LOAD  SURING
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407
          IF(WLEV.GT.PEAKWLEV) GOTO 407
      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
      WFX = VH1(1)*WFF ; WFY = VH1(2)*WFF ; WFZ = VH1(3)*WFF ;
      CASE(7) !WIND LOAD 
          WLEV = RLEV - SEABED - WDEPTH
          IF(WLEV.LT.0.0D0) GOTO 407 
          ! -------------------------------------------------------------------------------------------------------
          ! GENERATE WIND LOAD ( DNV,API,IEC64100-1 )
          ! NWINDO = 1 ; LOGARITHMIC PROFILE   (DNV) 
          ! NWINDO = 2 ; POWER LAW PROFILE     (DNV)
          ! NWINDO = 3 ; WIND PROFILE AND GUST (API)
          ! NWINDO = 4 ; IEC 614000-1
          ! NWINDO = 5 ; USER DEFINED
            IF (NWINDO.EQ.1)THEN ! DNV > LOGARITHMIC PROFILE
            UZT=UH*(LOG(WLEV/Z0)/LOG(RHIGH/Z0))
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
            
            ELSEIF (NWINDO.EQ.2.OR.NWINDO.EQ.5)THEN ! DNV > POWER LAW PROFILE NWINDO >> 2, USER DEFIND NWINDO >> 5 
                 IF (NWINDO.EQ.5)THEN ! SAME EQUATION FROM USER DEFINED AND POWER LAW IN DNV
                 UH=UHD               ! UH  = MEAN WIND SPEED ON DNV 
                 ENDIF                ! UHD = MEAN WIND SPEED ON USER DEFINED
            CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
     1                        WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
     2                        CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
     
            ELSEIF (NWINDO.EQ.3)THEN ! API > WIND PROFILE AND GUST
            C=0.0573D0*SQRT(1D0+0.0457D0*UHAPI)
            UZ=UHAPI*(1+(C*LOG(WLEV/32.8D0)))
            ZI=0.06D0*(1+0.0131D0*UHAPI)*((WLEV/32.8D0)**(-0.22))
            UZT=UZ*(1-0.41D0*ZI*LOG(AVERAGE/3600D0))
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
            
            ELSEIF (NWINDO.EQ.4)THEN ! IEC 614000-1
            UZT=UHD*((WLEV/RHIGH)**ALPHA)
            WFF=0.5D0*RHOA*DIAM*UZT*UZT
            ENDIF
         ! ----------------------------------------------------------------------------------------------------------
         ! CSA = SHAPE COEFFICIENT FOR WIND LOAD
      WFX = VWIND(1)*WFF*CSA ; WFY = VWIND(2)*WFF*CSA ; WFZ = VWIND(3)*WFF*CSA ;
      ENDSELECT
      
	TT(1) = TT(1) + WFX
	TT(2) = TT(2) + WFY
	TT(3) = TT(3) + WFZ
	
407   CONTINUE
C     -------------------------- 

      

	DO 502 INO=1,NNO
	FAC = H(INO)*WT*DET
      DO 481 I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 502
 481  RL(I,INO) = RL(I,INO) + TT(ICO)*FAC
 502  CONTINUE

 800  CONTINUE


C     MUST BE ADDED LATER  FOR SHELL FIXEND FORCED
C	CALL FIXSHE (RL,KEG,IELE,NNO,3)

C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,4)
C
 700  DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + RL(INF,INO) !VARY LOAD
      IF (IEQ.NE.0.AND.ILCC.GT.0) RRC(IEQ) = RRC(IEQ) + RL(INF,INO) !CONSTANT LOAD
 600  CONTINUE

      

 900  CONTINUE

      
      IF(KOFFL.NE.0) THEN
        IF (ISOLOP.EQ.1.AND.ICONTROL.EQ.0) CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'SPEC','WRIT') !STORE VARY LOAD
        IF (ISOLOP.EQ.1.AND.ICONTROL.EQ.1) CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
C        ====== THIS PART ARE NOT PROVIDED ===== 07-12-2013
C        CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','WRIT') !STORE CONSTANT LOAD
      ENDIF      
      
      IF (ISOLOP.EQ.1) CALL MUTIFORCE (ITIME,ILCAS,ISPEC,NOFFL)

7501  CONTINUE
7500  CONTINUE

	DEALLOCATE(RRV,RRC,RVAL,IVAL)
C
1000  RETURN
      END
C
C	=====================================================================
C	=====================================================================

	SUBROUTINE SHAP4C1(H,RI,SI)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	C1 SHAPE FUNCTION FOR 4 NODE PLATE ELEMENT
	DIMENSION H(3,1),D(12)
	
C	W, Theta R, Theta S
      D(1)  =  0.25D0+0.375D0*RI+0.375D0*SI+0.5D0*RI*SI-0.125D0*RI**3
     #	    -0.125D0*SI**3-0.125D0*RI**3*SI-0.125D0*RI*SI**3
      D(2)  = -0.125D0-0.125D0*RI-0.125D0*SI-0.125D0*RI*SI+0.125D0*SI**2
     #        +0.125D0*SI**2*RI+0.125D0*SI**3+0.125D0*RI*SI**3
      D(3)  = -0.125D0-0.125D0*RI-0.125D0*SI+0.125D0*RI**2-0.125D0*RI*SI
     #        +0.125D0*RI**3+0.125D0*RI**2*SI+0.125D0*RI**3*SI
      D(4)  =  0.25D0-0.375D0*RI+0.375D0*SI-0.5D0*RI*SI+0.125D0*RI**3
     #        -0.125D0*SI**3+0.125D0*RI**3*SI+0.125D0*RI*SI**3
      D(5)  = -0.125D0+0.125D0*RI-0.125D0*SI+0.125D0*RI*SI+0.125D0*SI**2
     #        -0.125D0*SI**2*RI+0.125D0*SI**3-0.125D0*RI*SI**3
      D(6)  =  0.125D0-0.125D0*RI+0.125D0*SI-0.125D0*RI**2-0.125D0*RI*SI
     #        +0.125D0*RI**3-0.125D0*RI**2*SI+0.125D0*RI**3*SI
      D(7)  =  0.25D0-0.375D0*RI-0.375D0*SI+0.5D0*RI*SI+0.125D0*RI**3
     #        +0.125D0*SI**3-0.125D0*RI**3*SI-0.125D0*RI*SI**3
      D(8)  =  0.125D0-0.125D0*RI-0.125D0*SI+0.125D0*RI*SI-0.125D0*SI**2
     #        +0.125D0*SI**2*RI+0.125D0*SI**3-0.125D0*RI*SI**3
      D(9)  =  0.125D0-0.125D0*RI-0.125D0*SI-0.125D0*RI**2+0.125D0*RI*SI
     #        +0.125D0*RI**3+0.125D0*RI**2*SI-0.125D0*RI**3*SI
      D(10) =  0.25D0+0.375D0*RI-0.375D0*SI-0.5D0*RI*SI-0.125D0*RI**3
     #        +0.125D0*SI**3+0.125D0*RI**3*SI+0.125D0*RI*SI**3
      D(11) =  0.125D0+0.125D0*RI-0.125D0*SI-0.125D0*RI*SI-0.125D0*SI**2
     #        -0.125D0*SI**2*RI+0.125D0*SI**3+0.125D0*RI*SI**3
      D(12) = -0.125D0-0.125D0*RI+0.125D0*SI+0.125D0*RI**2+0.125D0*RI*SI
     #        +0.125D0*RI**3-0.125D0*RI**2*SI-0.125D0*RI**3*SI


C	REVERSE ROTATION ABOUT SI
	D(3)  = -1.0D0*D(3 )
	D(6)  = -1.0D0*D(6 )
	D(9)  = -1.0D0*D(9 )
	D(12) = -1.0D0*D(12)

	K = 0
	DO I = 1,4
	DO J = 1,3
	K = K + 1
	H(J,I) = D(K)
	ENDDO
	ENDDO


      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
