C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE ELGLOL (LEST,NEGLO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C	ELEMENT GLOBAL LOAD
C	---	ELEMENT POINT LOAD
C	---	ELEMENT LINE  LOAD
C	---	ELEMENT PATCH LOAD
C     -----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /LCSS/ ILCN,ILCC

	COMMON /GiDEle/ LGID 

      COMMON A(9000000),IA(9000000)

	ALLOCATABLE XYP(:,:),PVL(:,:),LOD(:) 
	ALLOCATABLE XLUMP(:,:),PLUMP(:,:) 
	ALLOCATABLE IFDN(:) 


      DIMENSION LEST(1),LTYPE(20),COORP(3),GG(3),POS(3)


	IF(NEGLO.LE.0) RETURN
	
	READ(ITI,*)

	KLDN = 100
C     -----------------------------------------------------
	DO 5000	IEGLO = 1,NEGLO
	

C	KLN = NUMBER OF LOAD POINT   KLM = NUMBER OF ELEMENT IN PATCH
	READ(ITI,*) II,LTYP,KLN,KLM,ILCN,ILCC


C	============================
	SELECTCASE(LTYP)

C	IOPT  --->  0 = POINT LOAD   1 = LINE LOAD   2 = PATCH LOAD
C	----------------
	CASE(1) !PLANE POINT LOAD
	NTYPE = 2
	LTYPE(1:2) = [6,8]
	IOPT  = 0 
	NCOL  = 2   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 2   !NUMBER OF NATURAL COORDINATE
C	----------------
	CASE(2) !PLANE LINE LOAD
	NTYPE = 2
	LTYPE(1:2) = [6,8]
	IOPT  = 1 
	NCOL  = 2   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 2   !NUMBER OF NATURAL COORDINATE
C	----------------
	CASE(3) !PLANE PATCH LOAD
	NTYPE = 2
	LTYPE(1:2) = [6,8]
	IOPT  = 2
	NCOL  = 2   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 2   !NUMBER OF NATURAL COORDINATE
C	----------------

C	----------------
	CASE(4) !SHELL POINT LOAD
	NTYPE = 1
	LTYPE(1) = 9
	IOPT  = 0 
	NCOL  = 3   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 2   !NUMBER OF NATURAL COORDINATE
C	----------------
	CASE(5) !SHELL LINE LOAD
	NTYPE = 1
	LTYPE(1) = 9
	IOPT  = 1 
	NCOL  = 3   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 2   !NUMBER OF NATURAL COORDINATE
C	----------------
	CASE(6) !SHELL PATCH LOAD
	NTYPE = 1
	LTYPE(1) = 9
	IOPT  = 2 
	NCOL  = 3   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 2   !NUMBER OF NATURAL COORDINATE
C	----------------

C	----------------
	CASE(7) !SOLID POINT LOAD
	NTYPE = 2
	LTYPE(1:2) = [10,11]
	IOPT  = 0 
	NCOL  = 3   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 3   !NUMBER OF NATURAL COORDINATE
C	----------------
	CASE(8) !SOLID LINE LOAD
	NTYPE = 2
	LTYPE(1:2) = [10,11]
	IOPT  = 1 
	NCOL  = 3   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 3   !NUMBER OF NATURAL COORDINATE
C	----------------
	CASE(9) !SOLID PATCH LOAD
	NTYPE = 2
	LTYPE(1:2) = [10,11]
	IOPT  = 2 
	NCOL  = 3   !NUMBER OF GLOBAL  COORDINATE
	NPR   = 3   !NUMBER OF NATURAL COORDINATE
C	----------------

	CASE DEFAULT
	RETURN

	ENDSELECT
C	============================

	ALLOCATE( XYP(NCOL,KLN),PVL(NCOL,KLN),LOD(KLM) )
	ALLOCATE( XLUMP(NCOL,KLDN),PLUMP(NCOL,KLDN) )	
	ALLOCATE( IFDN(KLDN) )


C	READ LOAD POSITIONS AND VALUES
	DO J = 1,KLN
	READ(ITI,*)     PVL(1:NCOL,J),XYP(1:NCOL,J)
	ENDDO

C	READ ELEMENT IN PATCH
	DO ILM = 1,KLM
	READ(ITI,*)     LOD(ILM)
	ENDDO	


C     ----------------------------------------------

C	GENERATE LUMP LOAD USING GAUSS POINT METHOD	
	CALL PATCHL(XLUMP,PLUMP,XYP,PVL,NCOL,KLN,MLDN,IOPT)


	IFDN(1:MLDN) = 0
C     ----------------------------------------------
	DO 1000 IEG = 1,NEG
	KEG = IEG
	NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     
      CALL MOVLEV (2)

C	SCAN FOR ELEMENT TYPE
	IFIND = 0
	DO I = 1,NTYPE
	IF(ITYPE.EQ.LTYPE(I)) IFIND = 1
	ENDDO
	IF(IFIND.EQ.0) GOTO 1000
	

C     ------------------------
	DO 500 ILDN = 1,MLDN   !MLDN
	IF(IFDN(ILDN).EQ.1) GOTO 500

	COORP(1:NCOL) = XLUMP(1:NCOL,ILDN)


	TOL   = 1.001
	IFIND = 0
	
	DO 400 ILM = 1,KLM

	IELE = LOD(ILM)
	
	CALL ELEREODER(IA(LGID),NELE,IELE)

	CALL PSEARCH(POS,A(LCO),IA(LEX),COORP,NPR,NCO,NNM,NEX,IELE,
	1			 SLEN)

	ITEST = 0
	DO IPR = 1,NPR
	IF(ABS(POS(IPR)).GT.TOL) ITEST = 1
	IF(SLEN.GT.0.01) ITEST = 1
	ENDDO

C	+++++++++++++++++++++
	IF(ITEST.EQ.0) THEN
	IFDN(ILDN) = 1
	

	GG(1:3) = 0.0D0
	DO ICOL = 1,NCOL
	GG(ICOL) = PLUMP(ICOL,ILDN)
	ENDDO

	SELECTCASE(ITYPE)
	CASE(6,8)
	CALL MEPOIN (IELE,GG,POS,A(LCO),IA(LGS),A(LGP),
	1			 IA(LLM),A(LMP),IA(LMS)) 

	CASE(9)
	CALL SHPOIN (IELE,GG,POS,A(LCO),IA(LGS),A(LGP),
	1			 IA(LLM),A(LMP),IA(LMS),IA(LEX)) 


	CASE(10,11)
	CALL SOPOIN (IELE,GG,POS,A(LCO),IA(LGS),A(LGP),
	1			 IA(LLM),A(LMP),IA(LMS)) 

	ENDSELECT

	EXIT
	ENDIF
C	+++++++++++++++++++++


400	CONTINUE



500	CONTINUE
C     ------------------------

	ITEST = 0
	DO ILDN = 1,MLDN
	ITEST = ITEST + IFDN(ILDN)
	ENDDO
	IF(ITEST.EQ.MLDN) EXIT
	

1000	CONTINUE
C     ----------------------------------------------

	DEALLOCATE( XYP,PVL,LOD )
	DEALLOCATE( XLUMP,PLUMP )
	DEALLOCATE( IFDN )	

5000	CONTINUE
C     -----------------------------------------------------


      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SOPOIN (IELE,GG,POS,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     READS AND GENERATES NODAL LOADS, ADDS CONTRIBUTIONS TO LOAD V
C     -------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
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

	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C	==================================================================

	COMMON A(9000000),IA(9000000)

      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)
      
	DIMENSION XJ(3,3),XJI(3,3)
	DIMENSION RL(NNF,NNM),XYZ(NCO*NNM,NELE),LM(NEF,NELE),
	1		  H(NNM),P(3,NNM),XY(3,NNM),AJ(9)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION GG(3),POS(3)

C     -----------------------------------
      ALLOCATE(RAL(NEF),MAL(NEF))
C     -----------------------------------
      
C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
	RI = POS(1)
	SI = POS(2)
	TI = POS(3)
      
      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0

	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,IELE)
	ENDDO
	ENDDO

      CALL CLEARA (RL,NNF*NNM)

C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	IF(NNM.EQ.4.OR.NNM.EQ.10) GOTO 250

      CALL SHAP3D8 (RI,SI,TI,H,P) 

	XJ = MATMUL(P,TRANSPOSE(XY))
      K=1
      DO 110 J=1,3
      DO 110 I=1,3
	AJ(K)=XJ(I,J)
	K=K+1
110   CONTINUE

      DET = AJ(1)*AJ(5)*AJ(9) + AJ(4)*AJ(8)*AJ(3) + AJ(7)*AJ(2)*AJ(6)
     1    - AJ(7)*AJ(5)*AJ(3) - AJ(4)*AJ(2)*AJ(9) - AJ(1)*AJ(8)*AJ(6)
	
	GOTO 251
250	CONTINUE

	CALL GAUSST (RI,SI,TI,WT,IPT,MGR,1)
	CALL SHAP3DT(RI,SI,TI,H,P,NNM)
	CALL JACO3D(XY,P,XJ,XJI,DET,IELE,NNM)
	DET = DET/6.0

251	CONTINUE

	
      DO 390  INO=1,NNM
      FAC = H(INO)
      DO 380 I = 1,3
	RL(I,INO) = RL(I,INO) + GG(I)*FAC
 380	CONTINUE
 390  CONTINUE



	CALL FIXSOL (RL,KEG,IELE,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
 711  CONTINUE

      
      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)
      
      DEALLOCATE(RAL,MAL)
      

      RETURN

      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SHPOIN (IELE,GG,POS,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET,NODEX) 
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

      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	==================================================================
C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)
      
      DIMENSION PROPG(NGP,1),IGSET(1),XYZ(NCO*NNM,1),NODEX(NEX,1)
	DIMENSION LM(NEF,1)
      DIMENSION H(9),P(2,9),XJI(4),FA(4),G(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,9),IGPOS(9),GG(3),POS(2)


      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)

C     -----------------------------------
      ALLOCATE(RAL(NEF),MAL(NEF))
C     -----------------------------------
      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0

	RI = POS(1)
	SI = POS(2)

      ISET = IGSET(IELE)
	THICK = PROPG(2,ISET)

      CALL CLEARA (RL,3*NNO)

C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	IF(NNO.NE.3) THEN
      IF(NNO.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX(1,IELE),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,IELE),NNO)            !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,IELE),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	ELSEIF(NNO.EQ.3) THEN                                             !SHELL 3 NODE
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,IELE),P,VR,VS,VT,FA,XJI,DET,IELE,NNO)
	ENDIF
C     --------------------------------------------------------

      DO 390  INO=1,NNO
      FAC = H(INO)
      DO 380  I=1,3
	ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 390
 380  RL(I,INO) = RL(I,INO) + GG(ICO)*FAC
 390  CONTINUE

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
	
      DEALLOCATE(RAL,MAL)

C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE MEPOIN (IELE,GG,POS,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     READS AND GENERATES NODAL LOADS, ADDS CONTRIBUTIONS TO LOAD V
C     -------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
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

	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C	==================================================================

	COMMON A(9000000),IA(9000000)

      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)
      
	DIMENSION RL(NNF,NNM),IGSET(NELE),XYZ(NCO*NNM,NELE),LM(NEF,NELE),
	1		  H(NNM),P(2,NNM),XY(2,NNM),XJ(2,2)

	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1)
	DIMENSION GG(2),NODEX(4),POS(2)
C     ---------------
C     INITIALISATION
C     --------------
	NODEX(1:4) = [5,6,7,8]

C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
	RI = POS(1)
	SI = POS(2)

C     -----------------------------------
      ALLOCATE(RAL(NEF),MAL(NEF))
C     -----------------------------------
      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0
      

 202  KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,2
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,IELE)
	ENDDO
	ENDDO

      CALL CLEARA (RL,NNF*NNM)

C     ----------------
C     GAUSS POINT LOOP
C     ----------------

	CALL SHAP2D (RI,SI,H,P,NODEX,NNM)

	XJ = MATMUL(P,TRANSPOSE(XY))
     

	XBAR = 0.0D0
	DO INO = 1,NNM
	XBAR = XBAR + H(INO)*XY(1,INO)
	ENDDO
	IF(ISTYP.NE.0.AND.ISTYP.NE.3) XBAR = 1.0D0
	
C     ----------------------------------------
C     DETERMINANT OF THE JACOBIAN MATRIX (DET)
C     ----------------------------------------
      DET = XJ(1,1)*XJ(2,2) - XJ(1,2)*XJ(2,1) 

      DO 390  INO=1,NNM
      FAC = H(INO)
      DO 380 I = 1,2
	RL(I,INO) = RL(I,INO) + GG(I)*FAC*XBAR
 380	CONTINUE
 390  CONTINUE

	
	CALL FIXMRE (RL,KEG,IELE,NEF,NNM)


	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
 711	CONTINUE


      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,2
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)
      
      DEALLOCATE(RAL,MAL)


      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================

	SUBROUTINE PATCHL(XLUMP,PLUMP,COORD,PVALU,NCO,NNM,MLDN,IOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT

	DIMENSION COORD(NCO,1),PVALU(NCO,1)

	DIMENSION XLUMP(NCO,1),PLUMP(NCO,1)

	DIMENSION X12(NCO,2),P12(NCO,2)
	DIMENSION X14(NCO,4),P14(NCO,4)
	DIMENSION H4(4),P4(2,4)
	DIMENSION H1(2),P1(2)

	DIMENSION DXYZ(NCO),GXYZ(NCO),PXYZ(NCO)

	DIMENSION VR(3),VS(3),VT(3),XJI(2,2),FA(2,2)


	PLUMP(1:NCO,1:NNM) = 0.0D0

	SELECTCASE(IOPT)

	CASE(0)
C	-----------------------------------------
	MLDN = NNM
	XLUMP(1:NCO,1:MLDN) = COORD(1:NCO,1:MLDN)
	PLUMP(1:NCO,1:MLDN) = PVALU(1:NCO,1:MLDN)
C	-----------------------------------------

	
	CASE(1)
C	-----------------------------------------
	MLDN = 0
	MGR  = 7
	DO ISP = 1,NNM-1

	X12(1:NCO,1) = COORD(1:NCO,ISP+0)
	X12(1:NCO,2) = COORD(1:NCO,ISP+1)
	P12(1:NCO,1) = PVALU(1:NCO,ISP+0)
	P12(1:NCO,2) = PVALU(1:NCO,ISP+1)


	DO IGR = 1,MGR
	RI = GLOC(IGR,MGR)
	WR =  GWT(IGR,MGR)
	MLDN = MLDN + 1    !!!

	CALl POIN1D(RI,H1,P1)

	DVOL = 0.0D0
	DO ICO = 1,NCO
	DXYZ(ICO) = P1(1)*X12(ICO,1) + P1(2)*X12(ICO,2)
	GXYZ(ICO) = H1(1)*X12(ICO,1) + H1(2)*X12(ICO,2)
	PXYZ(ICO) = H1(1)*P12(ICO,1) + H1(2)*P12(ICO,2)
	DVOL = DVOL + DXYZ(ICO)*DXYZ(ICO)
	ENDDO
	DVOL = SQRT(DVOL)
	DVOL = DVOL*WR

	DO ICO = 1,NCO
	XLUMP(ICO,MLDN) = GXYZ(ICO)
	PLUMP(ICO,MLDN) = PXYZ(ICO)*DVOL
	ENDDO

	ENDDO
	ENDDO
C	-----------------------------------------

	
	CASE(2)
C	-----------------------------------------
	DO INM = 1,4
	X14(1:NCO,1:INM) = COORD(1:NCO,1:INM)
	P14(1:NCO,1:INM) = PVALU(1:NCO,1:INM)
	ENDDO

	MLDN = 0
	MGR  = 7
	MGS  = 7
	DO IGR = 1,MGR
	RI = GLOC(IGR,MGR)
	WR =  GWT(IGR,MGR)
	DO IGS = 1,MGS
	SI = GLOC(IGS,MGS)
	WS =  GWT(IGS,MGS)
	WW = WR*WS
	MLDN = MLDN + 1    !!!

	CALl POIN2D(RI,SI,H4,P4)
	CALL POIN2J(NNM,X14,P4,VR,VS,VT,XJI,DET,RI,SI,NCO)
	DVOL = DET*WW

	DO ICO = 1,NCO
	XVAL = 0.0D0
	PVAL = 0.0D0
	DO INM = 1,4
	XVAL = XVAL + H4(INM)*COORD(ICO,INM)
	PVAL = PVAL + H4(INM)*PVALU(ICO,INM)
	ENDDO
	PLUMP(ICO,MLDN) = PVAL*DVOL
	XLUMP(ICO,MLDN) = XVAL
	ENDDO

	ENDDO
	ENDDO
C	-----------------------------------------

	ENDSELECT


	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================

      SUBROUTINE POIN1D (R,H,P)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      DIMENSION  H(2),P(2)

	H(1) = 0.5*(1.0-R)
	H(2) = 0.5*(1.0+R)

	P(1) =-0.5
	P(2) = 0.5

      RETURN
      END
C
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE POIN2D (R,S,H,P)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
      DIMENSION  H(1),P(2,1),PP(2,2,4)
C
      RP  = 1.0+R
      SP  = 1.0+S
      RM  = 1.0-R
      SM  = 1.0-S
      R2  = 1.0-R*R
      S2  = 1.0-S*S
C     ---------------------------------------------
C     INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     FOR A FOUR NODE ELEMENT
C     ---------------------------------------------
      H(1)   = 0.25*RP*SP
      H(2)   = 0.25*RM*SP
      H(3)   = 0.25*RM*SM
      H(4)   = 0.25*RP*SM
      P(1,1) = 0.25*SP
      P(1,2) = -P(1,1)
      P(1,3) = -0.25*SM
      P(1,4) = -P(1,3)
      P(2,1) = 0.25*RP
      P(2,2) = 0.25*RM
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)

	PP(1,1,1) = 0.0
	PP(1,1,2) = 0.0
	PP(1,1,3) = 0.0
	PP(1,1,4) = 0.0

	PP(1,2,1) = 0.25
	PP(1,2,2) =-0.25
	PP(1,2,3) = 0.25
	PP(1,2,4) =-0.25

	PP(2,1,1) = 0.25
	PP(2,1,2) =-0.25
	PP(2,1,3) = 0.25
	PP(2,1,4) =-0.25

	PP(2,2,1) = 0.0
	PP(2,2,2) = 0.0
	PP(2,2,3) = 0.0
	PP(2,2,4) = 0.0


      RETURN
      END
C
C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE POIN2J(NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,NCO)
	IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------
      DIMENSION COORD(NCO,1),HD(2,1),VR(3),VS(3),VT(3),XJI(4)
      DIMENSION COVR(3),COVS(3),RV(3),SV(3)
C
      CALL CLEARA (COVR,3)
      CALL CLEARA (COVS,3)
      DO 20 I=1,NNO
      DO 20 J=1,NCO
      COVR(J)=COVR(J)+HD(1,I)*COORD(J,I)
   20 COVS(J)=COVS(J)+HD(2,I)*COORD(J,I)
      CALL VECPRD (COVR,COVS,VT)
      CALL SCALEN (VT,VT,DET,3)
      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)

      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)

      RETURN

      END


C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE PSEARCH(POS,COORD,NODEX,COORP,NPR,NCO,NNM,NEX,IELE,
	1				   SLEN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION COORD(NCO*NNM,1),COORP(NCO),NODEX(NEX,1)
	DIMENSION F(NPR),HH(NPR,NPR),HHI(NPR,NPR),POS(NPR)

	NTER = 10

C	ASSUMED THE INITIAL VALUE OF RI,SI,TI
	POS(1:NPR) = 0.0D0
	 
	DO 100 ITER = 1,NTER

	CALL PDERIV(POS,COORD(1,IELE),COORP,F,HH,SLEN,NPR,NCO,NNM,
	1			NODEX(1,IELE))


	CALL INVMATRIX(HH,HHI,NPR)

	DO I = 1,NPR
	DO J = 1,NPR
	POS(I) = POS(I) - HHI(I,J)*F(J)
	ENDDO
	ENDDO


C	WRITE(*,*) IELE,COORP
C	PAUSE

100	CONTINUE


	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================

	SUBROUTINE PDERIV(POS,COORD,COORP,F,HH,SLEN,NPR,NCO,NNM,NODEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION H(NNM),P(NPR,NNM),COORD(NCO,NNM),COORP(NCO),POS(NPR)
	DIMENSION PP(NPR,NPR,NNM),X(NCO,NNM),XP(NCO,NNM+1)
	DIMENSION F(NPR),HH(NPR,NPR),HP(NNM+1),NODEX(1)

C	INITIALIZE THE SECONDERIVATIVE OF SHAPE FUNCTION
	PP = 0.0D0

C	+++++++++++++++++++++++++++
	SELECTCASE(NPR)

	CASE(1)
	RI = POS(1)
	CALL POIN1D (RI,H,P)

	CASE(2)
	RI = POS(1)
	SI = POS(2)
	IF(NNM.NE.3) THEN
		IF(NNM.NE.9) CALL SHAP2D (RI,SI,H,P,NODEX,NNM)
		IF(NNM.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX,NNM)				!9 NODE ELEMENT
	ELSEIF(NNM.EQ.3) THEN											!3 NODE ELEMENT
		CALL SHAP2D3(RI,SI,H,P,NNM)
	ENDIF

	CASE(3)
	RI = POS(1)
	SI = POS(2)
	TI = POS(3)
	IF(NNM.EQ.4.OR.NNM.EQ.10) GOTO 250
		CALL SHAP3D8 (RI,SI,TI,H,P) 
		GOTO 251
250	CONTINUE
		CALL SHAP3DT (RI,SI,TI,H,P,NNM)
251	CONTINUE

	ENDSELECT
C	+++++++++++++++++++++++++++
	
	HP(1+NNM) = -1.0D0
	HP(1:NNM) = H(1:NNM)


	 X(1:NCO,1:NNM) = COORD(1:NCO,1:NNM)
	XP(1:NCO,1:NNM) = COORD(1:NCO,1:NNM)

	XP(1:NCO,1+NNM) = COORP(1:NCO)

	DO I = 1,NPR
	F(I) = 0.0D0
	DO M = 1,NCO
	DO K = 1,NNM+1
	DO N = 1,NNM
	F(I) = F(I) + 2.0*HP(K)*XP(M,K)*P(I,N)*X(M,N)
	ENDDO
	ENDDO
	ENDDO
	ENDDO


	DO I = 1,NPR
	DO J = 1,NPR
	HH(I,J) = 0.0D0
	DO M = 1,NCO

	DO K = 1,NNM
	DO N = 1,NNM
	HH(I,J) = HH(I,J) + 2.0*P(J,K)*X(M,K)*P(I,N)*X(M,N)
	ENDDO
	ENDDO
	DO K = 1,NNM+1
	DO N = 1,NNM
	HH(I,J) = HH(I,J) + 2.0*HP(K)*XP(M,K)*PP(I,J,N)*X(M,N)
	ENDDO
	ENDDO

	ENDDO
	ENDDO
	ENDDO


	SLEN = 0.0D0
	DO M = 1,NCO
	DO K = 1,NNM+1
	DO N = 1,NNM+1
	SLEN = SLEN + HP(K)*XP(M,K)*HP(N)*XP(M,N)
	ENDDO
	ENDDO
	ENDDO
	
	IF(SLEN.GT.0.0D0) SLEN = SQRT(SLEN)



	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================
