C	=======================================================================
C	=================== START OF SOLID SELFWEIGHT LOAD ====================
C	=======================================================================
      SUBROUTINE SOLWTH (NGV,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET,IGIDM) 
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

C	NEXT LINE ADDED BY GILSON - MARCH2004 (LOAD INPUT)
	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C	==================================================================

	COMMON A(9000000),IA(9000000)

      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
	DIMENSION XJ(3,3),XJI(3,3)
	DIMENSION RL(NNF,NNM),XYZ(NCO*NNM,NELE),LM(NEF,NELE),
	1		  H(NNM),P(3,NNM),XY(3,NNM),GLOC(2),GWT(2),AJ(9)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION IGIDM(NELE),GG(3)

      
      ALLOCATE(RAL(NEF),MAL(NEF))
      
C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
      IGV = 0

 200	READ (ITI,*) MEMBA,GG(1),GG(2),GG(3),ILCN,ILCC
	IGV = IGV + 1

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

	MSET = MTSET(MEMBA)
	DENS = PROPM(5,MSET)

	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO

      RL = 0.0

	MGR = 2
	MGS = 2
	MGT = 2

	IF(NNM.EQ.4.OR.NNM.EQ.10) THEN
	MGR = 4
	IF(NNM.EQ.10) MGR = 4
	MGS = 1
	MGT = 1
	ENDIF

	GLOC(1) = -1.0/SQRT(3.0)
	GLOC(2) =  1.0/SQRT(3.0)
	GWT(1) = 1.0
	GWT(2) = 1.0

	
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	IPT = 0
 121  DO 800  IGR=1,MGR
      RI = GLOC(IGR)
      DO 800  IGS=1,MGS
      SI = GLOC(IGS)
	DO 800  IGT=1,MGT
	TI = GLOC(IGT)
      WT = GWT(IGR)*GWT(IGS)*GWT(IGT)
	IPT = IPT + 1
	
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
	CALL JACO3D(XY,P,XJ,XJI,DET,MEMBA,NNM)
	DET = DET/6.0

251	CONTINUE

	
      DO 390  INO=1,NNM
      FAC = H(INO)*WT*DET
      DO 380 I = 1,3
	RL(I,INO) = RL(I,INO) + DENS*GG(I)*FAC
 380	CONTINUE
	
 390  CONTINUE

	
 800  CONTINUE

	CALL FIXSOL (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE

      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)

 290  IF (IGV.LT.NGV) GOTO 200


      DEALLOCATE(RAL,MAL)	  

      RETURN

      END
C
C	=======================================================================
C	=======================================================================
      SUBROUTINE SOLWTH_OFFSHORE (NGV,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET,IGIDM) 
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

C	NEXT LINE ADDED BY GILSON - MARCH2004 (LOAD INPUT)
	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C	==================================================================

	COMMON A(9000000),IA(9000000)

C
      DIMENSION R(NEQ)
	DIMENSION XJ(3,3),XJI(3,3)
	DIMENSION RL(NNF,NNM),XYZ(NCO*NNM,NELE),LM(NEF,NELE),
	1		  H(NNM),P(3,NNM),XY(3,NNM),GLOC(2),GWT(2),AJ(9)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(NELE)
	DIMENSION IGIDM(NELE),GG(3)

C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
      IGV = 0

 200	READ (ITI,*) MEMBA,GG(1),GG(2),GG(3),ILCN,ILCC
	IGV = IGV + 1

      CALL CLEARA (R,NEQ)

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

	MSET = MTSET(MEMBA)
	DENS = PROPM(5,MSET)

	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO

      RL = 0.0

	MGR = 2
	MGS = 2
	MGT = 2

	IF(NNM.EQ.4.OR.NNM.EQ.10) THEN
	MGR = 4
	IF(NNM.EQ.10) MGR = 4
	MGS = 1
	MGT = 1
	ENDIF

	GLOC(1) = -1.0/SQRT(3.0)
	GLOC(2) =  1.0/SQRT(3.0)
	GWT(1) = 1.0
	GWT(2) = 1.0

	
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	IPT = 0
 121  DO 800  IGR=1,MGR
      RI = GLOC(IGR)
      DO 800  IGS=1,MGS
      SI = GLOC(IGS)
	DO 800  IGT=1,MGT
	TI = GLOC(IGT)
      WT = GWT(IGR)*GWT(IGS)*GWT(IGT)
	IPT = IPT + 1
	
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
	CALL JACO3D(XY,P,XJ,XJI,DET,MEMBA,NNM)
	DET = DET/6.0

251	CONTINUE

	
      DO 390  INO=1,NNM
      FAC = H(INO)*WT*DET
      DO 380 I = 1,3
	RL(I,INO) = RL(I,INO) + DENS*GG(I)*FAC
 380	CONTINUE
	
 390  CONTINUE

	
 800  CONTINUE
      
      CALL FIXSOL_OFF (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE

      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) R(IEQ) = R(IEQ) + RL(INF,INO)
 700  CONTINUE

	CALL LDASEM_OFFSHORE (R)

 290  IF (IGV.LT.NGV) GOTO 200



      RETURN

      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE SHAP3D8 (R,S,T,H,P)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     AT THE NODAL POINTS OF A CURVILINEAR ISOPARAMETRIC HEXAHEDRON
C     (8 NODES)
C	---------------
C                    NODE NUMBERING CONVENTION
C				   -------------------------
C                                  2
C                                 0
C                               . . .                   T.     S.
C                             .   .   .                  .    .
C                           .     .     .                .  .
C                         .       0       .              .
C                       .       .  6 .      .              .
C                     .       .        .      .              .
C                3  .       .            .      .              .R
C                 0       .                .      .
C                 . .   .                    .      .  1
C                 .   .                        .      0
C                 .     .                        .  . .
C                 0       .                       .   . 
C                7  .       .                   .     .
C                     .       .               .       0
C                       .       .           .       .  5
C                         .       .       .       .
C                           .       .4  .       .
C                             .       0       .
C                               .     .     . 
C                                 .   .   .
C                                   . . .
C                                     0
C                                      8
C	--------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     R,S,T      = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     H(NNO)     = INTERPOLATION (SHAPE) FUNCTIONS
C     P(3,NNO)   = FUNCTION DERIVATIVES WITH RESPECT TO R,S,T RESP.
C     NODEX(NEX) = POSITIONS OF MIDSIDE NODES (EXCESS NODES)
C     NNO        = NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ----------------------------------------------------------------
      DIMENSION H(8),P(3,8)


      RP  = 1.+R
      SP  = 1.+S
      TP  = 1.+T
      RM  = 1.-R
      SM  = 1.-S
      TM  = 1.-T
      RR  = 1.-R*R
      SS  = 1.-S*S
      TT  = 1.-T*T
C     ----------------------------------------------------------
C     INTERPOLATION FUNCTIONS AND DERIVATIVES FOR A 8 NODE BRICK
C     ----------------------------------------------------------
	H(1)   = .125*RP*SP*TP
      H(2)   = .125*RM*SP*TP
      H(3)   = .125*RM*SM*TP
      H(4)   = .125*RP*SM*TP
      H(5)   = .125*RP*SP*TM
      H(6)   = .125*RM*SP*TM
      H(7)   = .125*RM*SM*TM
      H(8)   = .125*RP*SM*TM
C
      P(1,1) = .125*SP*TP
      P(1,2) = -P(1,1)
      P(1,3) = -.125*SM*TP
      P(1,4) = -P(1,3)
      P(1,5) = .125*SP*TM
      P(1,6) = -P(1,5)
      P(1,7) = -.125*SM*TM
      P(1,8) = -P(1,7)
C
      P(2,1) = .125*RP*TP
      P(2,2) = .125*RM*TP
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)
      P(2,5) = .125*RP*TM
      P(2,6) = .125*RM*TM
      P(2,7) = -P(2,6)
      P(2,8) = -P(2,5)
C
      P(3,1) = .125*RP*SP
      P(3,2) = .125*RM*SP
      P(3,3) = .125*RM*SM
      P(3,4) = .125*RP*SM
      P(3,5) = -P(3,1)
      P(3,6) = -P(3,2)
      P(3,7) = -P(3,3)
      P(3,8) = -P(3,4)


      RETURN

      END
C	=======================================================================
C	==================== END OF SOLID SELFWEIGHT LOAD =====================
C	=======================================================================


C	=======================================================================
C	================= START OF MEMBRANE SELFWEIGHT LOAD ===================
C	=======================================================================
      SUBROUTINE MEMWTH (NGV,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET,IGIDM) 
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

      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
	DIMENSION RL(NNF,NNM),IGSET(NELE),XYZ(NCO*NNM,NELE),LM(NEF,NELE),
	1		  H(NNM),P(2,NNM),XY(2,NNM),GLOC(2),GWT(2),XJ(2,2)

	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1)
	DIMENSION IGIDM(NELE),GLD(2),NODEX(4)
      
      ALLOCATE(RAL(NEF),MAL(NEF))    
      
C     ---------------
C     INITIALISATION
C     --------------
	NODEX(1:4) = [5,6,7,8]

C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
      IGV = 0

 200	READ (ITI,*,END=202) MEMBA,GLD(1),GLD(2),ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	IGV = IGV + 1

	MSET = MTSET(MEMBA)
	DENS = PROPM(5,MSET)

	MSET = IGSET(MEMBA)
	THICK = PROPG(2,MSET)

	IF(DENS.EQ.0.0D0) THEN
	WRITE(*,*) 'MATERIAL DENSITY IS EQUAL TO ZERO AT ELEMENT',MEMBA
	ENDIF

 202  KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,2
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO


      RL = 0.0

	MGR = 2
	MGS = 2

	GLOC(1) = -1.0/SQRT(3.0)
	GLOC(2) =  1.0/SQRT(3.0)
	GWT(1) = 1.0
	GWT(2) = 1.0

	
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	IPT = 0
 121  DO 800  IGR=1,2
      RI = GLOC(IGR)
      DO 800  IGS=1,2
      SI = GLOC(IGS)
      WT = GWT(IGR)*GWT(IGS)
	IPT = IPT + 1


C	CALL SHAP2D4 (RI,SI,H,P) 
	CALL SHAP2D  (RI,SI,H,P,NODEX,NNM)

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
      FAC = H(INO)*WT*DET*THICK
      DO 380 I = 1,2
	GLOAD = GLD(I)
	RL(I,INO) = RL(I,INO) + DENS*GLOAD*FAC*XBAR
 380	CONTINUE
 390  CONTINUE

	
 800  CONTINUE

	CALL FIXMRE (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE

      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,2
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)


 290  IF (IGV.LT.NGV) GOTO 200


      DEALLOCATE(RAL,MAL)	  

      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE SHAP2D4 (R,S,H,P)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION H(4),P(2,4)


      RP  = 1.+R
      SP  = 1.+S
      RM  = 1.-R
      SM  = 1.-S
C     ----------------------------------------------------------
C     INTERPOLATION FUNCTIONS AND DERIVATIVES FOR A 4 NODE
C     ----------------------------------------------------------
	H(1)   = .25*RP*SP
      H(2)   = .25*RM*SP
      H(3)   = .25*RM*SM
      H(4)   = .25*RP*SM

C
      P(1,1) =  0.25*SP
      P(1,2) = -P(1,1)
      P(1,3) = -0.25*SM
      P(1,4) = -P(1,3)

C
      P(2,1) =  0.25*RP
      P(2,2) =  0.25*RM
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)


      RETURN

      END
C	=======================================================================
C	================== END OF MEMBRANE SELFWEIGHT LOAD ====================
C	=======================================================================
	SUBROUTINE MEMPRS (NGV,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET,IGIDM) 
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
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

	COMMON A(9000000),IA(9000000)
      
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
	DIMENSION RL(NNF,NNM),IGSET(NELE),XYZ(NCO*NNM,NELE),
	1		  H(NNM),P(2,NNM),XY(2,NNM),GLOC(2),GWT(2),XJ(2,2)
	DIMENSION VR(3),VS(3),VT(3),LM(NEF,NELE)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1)
	DIMENSION IGIDM(NELE),NODEX(4) 
	  
      ALLOCATE(RAL(NEF),MAL(NEF))
      
C     ---------------
C     INITIALISATION
C     ---------------
	NODEX(1:4) = [5,6,7,8]

C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
      IGV = 0

200	READ (ITI,*) MEMBA,NFACE,LOPT,PREST,ZREF,ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	CALL CLEARA (RL,NEF)

	MSET  = IGSET(MEMBA)
	THICK = PROPG(2,MSET)

	IGV = IGV + 1

202	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,2
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO



	MGR = 2
	MGS = 2
	GLOC(1) = -1.0/SQRT(3.0)
	GLOC(2) =  1.0/SQRT(3.0)
	GWT(1) = 1.0
	GWT(2) = 1.0

	
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	IPT = 0
 121  DO 800  IGR=1,2
      RI = GLOC(IGR)
      WT = GWT(IGR)
	IPT = IPT + 1
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R= 1.0
C	FACE 2; R=-1.0
C	FACE 3: S= 1.0
C	FACE 4: S=-1.0
C	---------------------
C	---------------------------------------------
	VR = 0.0
	VS = 0.0
	VT = 0.0
	VS(1) = 0.0
	VS(2) = 0.0
	VS(3) = 1.0
C	---------------------------------------------
	SELECTCASE (NFACE)
	CASE(4)
	GR = 1.0D0
	GS = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ  = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(2,1)**2+XJ(2,2)**2)
	VR(1) = XJ(2,1)
	VR(2) = XJ(2,2)
	CASE(1)
	GS = 1.0D0
	GR = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(1,1)**2+XJ(1,2)**2)
	VR(1) = -XJ(1,1)
	VR(2) = -XJ(1,2)
	CASE(2)
	GR = -1.0D0
	GS = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(2,1)**2+XJ(2,2)**2)
	VR(1) = -XJ(2,1)
	VR(2) = -XJ(2,2)
	CASE(3)
	GS = -1.0D0
	GR = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(1,1)**2+XJ(1,2)**2)
	VR(1) = XJ(1,1)
	VR(2) = XJ(1,2)
	ENDSELECT
C	---------------------------------------------
	CALL SCALEN(VR,VR,DUM,3)
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
C	---------------------------------------------
	XBAR = 0.0D0
	DO INO = 1,NNM
	XBAR = XBAR + H(INO)*XY(1,INO)
	ENDDO
	IF(ISTYP.NE.0.AND.ISTYP.NE.3) XBAR = 1.0D0
C	---------------------------------------------
		
	IF(LOPT.EQ.0) THEN     !--------------------------
      DO INO=1,NNM
      FAC = H(INO)*WT*DET*PREST*THICK*XBAR
      DO I = 1,2
	RL(I,INO) = RL(I,INO) + VT(I)*FAC
	ENDDO
	ENDDO
	ELSEIF(LOPT.EQ.1) THEN !--------------------------
	GALEV = 0.0
	DO INO=1,NNM
	GALEV = GALEV + H(INO)*XY(NGRAV,INO)
	ENDDO
	IF(GALEV.GT.ZREF) GOTO 789
	DTH = ZREF-GALEV
	DO INO=1,NNM
      FAC = -H(INO)*WT*DET*PREST*DTH*THICK*XBAR
      DO I = 1,2
	RL(I,INO) = RL(I,INO) + VT(I)*FAC
	ENDDO
	ENDDO
789	CONTINUE
	ENDIF                  !--------------------------
C	---------------------------------------------
C	---------------------------------------------	
 800  CONTINUE
C	---------------------------------------------

	CALL FIXMRE (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE

      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,2
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)

 290  IF (IGV.LT.NGV) GOTO 200

      DEALLOCATE(RAL,MAL)	  

      RETURN

      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================	
	SUBROUTINE MEMTRAC (NGV,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET,IGIDM) 
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
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

	COMMON A(9000000),IA(9000000)
	  
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
	DIMENSION RL(NNF,NNM),IGSET(NELE),XYZ(NCO*NNM,NELE),
	1		  H(NNM),P(2,NNM),XY(2,NNM),GLOC(2),GWT(2),XJ(2,2)
	DIMENSION VR(3),VS(3),VT(3),LM(NEF,NELE)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1)
	DIMENSION IGIDM(NELE),NODEX(4),GG(3)
	  
      ALLOCATE(RAL(NEF),MAL(NEF))
      
C     ---------------
C     INITIALISATION
C     ---------------
	NODEX(1:4) = [5,6,7,8]

C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
      IGV = 0

200	READ (ITI,*) MEMBA,NFACE,GG(1),GG(2),ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	CALL CLEARA (RL,NEF)

	MSET  = IGSET(MEMBA)
	THICK = PROPG(2,MSET)

	IGV = IGV + 1

202	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,2
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO



	MGR = 2
	MGS = 2
	GLOC(1) = -1.0/SQRT(3.0)
	GLOC(2) =  1.0/SQRT(3.0)
	GWT(1) = 1.0
	GWT(2) = 1.0

	
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	IPT = 0
 121  DO 800  IGR=1,2
      RI = GLOC(IGR)
      WT = GWT(IGR)
	IPT = IPT + 1
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R= 1.0
C	FACE 2; R=-1.0
C	FACE 3: S= 1.0
C	FACE 4: S=-1.0
C	---------------------
C	---------------------------------------------
	VR = 0.0
	VS = 0.0
	VT = 0.0
	VS(1) = 0.0
	VS(2) = 0.0
	VS(3) = 1.0
C	---------------------------------------------
	SELECTCASE (NFACE)
	CASE(4)
	GR = 1.0D0
	GS = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ  = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(2,1)**2+XJ(2,2)**2)
	VR(1) = XJ(2,1)
	VR(2) = XJ(2,2)
	CASE(1)
	GS = 1.0D0
	GR = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(1,1)**2+XJ(1,2)**2)
	VR(1) = -XJ(1,1)
	VR(2) = -XJ(1,2)
	CASE(2)
	GR = -1.0D0
	GS = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(2,1)**2+XJ(2,2)**2)
	VR(1) = -XJ(2,1)
	VR(2) = -XJ(2,2)
	CASE(3)
	GS = -1.0D0
	GR = RI 
C	CALL SHAP2D4 (GR,GS,H,P) 
	CALL SHAP2D  (GR,GS,H,P,NODEX,NNM) 
	XJ = MATMUL(P,TRANSPOSE(XY))
	DET = SQRT(XJ(1,1)**2+XJ(1,2)**2)
	VR(1) = XJ(1,1)
	VR(2) = XJ(1,2)
	ENDSELECT
C	---------------------------------------------
	CALL SCALEN(VR,VR,DUM,3)
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
C	---------------------------------------------
	XBAR = 0.0D0
	DO INO = 1,NNM
	XBAR = XBAR + H(INO)*XY(1,INO)
	ENDDO
	IF(ISTYP.NE.0.AND.ISTYP.NE.3) XBAR = 1.0D0
C	---------------------------------------------

      DO INO=1,NNM
      FAC = H(INO)*WT*DET*THICK*XBAR
      DO I = 1,2
	RL(I,INO) = RL(I,INO) + GG(I)*FAC
	ENDDO
	ENDDO

C	---------------------------------------------
C	---------------------------------------------	
 800  CONTINUE
C	---------------------------------------------

	CALL FIXMRE (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711   CONTINUE
      
      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,2
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)
      
 290  IF (IGV.LT.NGV) GOTO 200


      DEALLOCATE(RAL,MAL)	

      RETURN

      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================	
	SUBROUTINE SEFTRS (NGV,XYZ,IGSET,PROPG,LM,
	1				   PROPM,MTSET,IGIDM,OPT) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPT
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

      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
      DIMENSION LM(NEF,NELE)
	DIMENSION RL(3,2),IGSET(NELE),XYZ(NCO*NNM,NELE),
	1		  XY(3,2)

	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1)
	DIMENSION IGIDM(NELE),GLD(3)
	  
      ALLOCATE(RAL(NEF),MAL(NEF))
      
C     --------------------
C     INPUT OF NODAL LOADS
C     --------------------
      IGV = 0


 200	READ (ITI,*) MEMBA,GLD(1),GLD(2),GLD(3),ILCN,ILCC

	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	IGV = IGV + 1

	MSET = MTSET(MEMBA)
	DENS = PROPM(5,MSET)

	MSET = IGSET(MEMBA)
	ARET = PROPG(2,MSET)

	IF(DENS.EQ.0.0D0) THEN
	WRITE(*,*) 'MATERIAL DENSITY IS EQUAL TO ZERO AT ELEMENT',MEMBA
	ENDIF

	KKAK = 0
	DO IKAK = 1,2
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO

	DX1 = XY(1,2)-XY(1,1)
	DX2 = XY(2,2)-XY(2,1)
	DX3 = XY(3,2)-XY(3,1)
	ELN = SQRT(DX1**2 + DX2**2 + DX3**2)
	WT = 0.5*ELN

      RL = 0.0
	
      DO 390 INO=1,2
      DO 380 I = 1,3
	GLOAD = GLD(I)
	RL(I,INO) = RL(I,INO) + WT*DENS*GLOAD*ARET
 380	CONTINUE
 390  CONTINUE

	CALL FIXTRS (RL,KEG,MEMBA,NEF,NNO)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE

      IEFL = 0
      DO 700  INO=1,2
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	IF (OPT.EQ."GEN") CALL LDASEM_NEW (RAL,MAL,IEFL) ! GENERAL LOAD STORAGE
      IF (OPT.EQ."OFF") CALL LDASEM_OFFSHORE (R)       ! OFFSHORE LOAD STORAGE

 290  IF (IGV.LT.NGV) GOTO 200

      DEALLOCATE(RAL,MAL)	  
C
      RETURN

      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE SEFSHS (PROPG,IGSET,XYZ,NODEX,LM,PROPM,MTSET,
     1                   MGP,MXY,MEX,MEF,IGIDM,NGV)
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
C     RT(3,NELE)       = GLOBAL TRACTIONS IN X,Y,Z
C     RP(NELE)         = NORMAL PRESSURES
C     R(NEQ)           = LOAD VECTOR (EQUIVALENT NODAL LOADS)
C     IND              = FLAG FOR TRACTIONS (IND=2),PRESSURES (IND=3)
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
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
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	==================================================================
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION RL(3,9),H(9),P(2,9),XJI(4),FA(4)
      DIMENSION VR(3),VS(3),VT(3)
	DIMENSION PROPM(NMP,1),MTSET(1),GLD(3),IGIDM(1)

C     -----------------------------------
      ALLOCATE(RAL(NEF),MAL(NEF))
C     -----------------------------------
      IGV = 0

200	READ (ITI,*) MEMBA,GLD(1),GLD(2),GLD(3),ILCN,ILCC
	IGV = IGV + 1


	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

      ISET = IGSET(MEMBA)
	THCK = PROPG(2,ISET)

	MSET = MTSET(MEMBA)
	DENS = PROPM(5,MSET)
	
C
      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
      CALL CLEARA (RL,27)


	MMGR = 3
	MMGS = 3
	IF(NNO.EQ.3) MMGS = 1   !SHELL 3 NODE

C     ----------------
C     GAUSS POINT LOOP
C     ----------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	IF(NNO.NE.3) THEN
      CALL SHAP2D (RI,SI,H,P,NODEX(1,MEMBA),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,MEMBA),NNO)            !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,MEMBA),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	ELSEIF(NNO.EQ.3) THEN                                             !SHELL 3 NODE
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,MEMBA),P,VR,VS,VT,FA,XJI,DET,MEMBA,NNO)
	ENDIF
C     --------------------------------------------------------
	
C     --------------------------------------------------------
C     UNIFORMLY DISTRIBUTED GLOBAL TRACTION (RT); CONSERVATIVE
C     --------------------------------------------------------
	DO 390  INO=1,NNO
      FAC = H(INO)*WT*DET*DENS*THCK
      DO 380  I=1,3
 380  RL(I,INO) = RL(I,INO) + FAC*GLD(I)
 390  CONTINUE


 800  CONTINUE

	CALL FIXSHE (RL,KEG,MEMBA,NNO,3)

C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1             A(LEI),RL,NSF,NNF,4)
C
 700  IEFL = 0
      DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 600  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)

 290  IF (IGV.LT.NGV) GOTO 200

      DEALLOCATE(RAL,MAL)	 
C
      RETURN
      END

C	=======================================================================
C	=======================================================================
      SUBROUTINE SEFSHS_OFF (PROPG,IGSET,XYZ,NODEX,LM,PROPM,MTSET,
     1                   MGP,MXY,MEX,MEF,IGIDM,NGV)
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
C     RT(3,NELE)       = GLOBAL TRACTIONS IN X,Y,Z
C     RP(NELE)         = NORMAL PRESSURES
C     R(NEQ)           = LOAD VECTOR (EQUIVALENT NODAL LOADS)
C     IND              = FLAG FOR TRACTIONS (IND=2),PRESSURES (IND=3)
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
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
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON A(9000000),IA(9000000)

	COMMON /LCSS/ ILCN,ILCC

C	==================================================================
C
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),NODEX(MEX,1),LM(MEF,1)
      DIMENSION RL(3,9),H(9),P(2,9),XJI(4),FA(4)
      DIMENSION VR(3),VS(3),VT(3)
	DIMENSION PROPM(NMP,1),MTSET(1),GLD(3),IGIDM(1)
	DIMENSION R(NEQ)

C     -----------------------------------
C     -----------------------------------
      IGV = 0

200	READ (ITI,*) MEMBA,GLD(1),GLD(2),GLD(3),ILCN,ILCC
	IGV = IGV + 1


	DO  IGID = 1,NELE
	MEMGD = IGIDM(IGID) 
	IF(MEMGD.EQ.MEMBA) THEN
	MEMBA = IGID 
	GOTO 201
	ENDIF
	ENDDO
201	CONTINUE

      ISET = IGSET(MEMBA)
	THCK = PROPG(2,ISET)

	MSET = MTSET(MEMBA)
	DENS = PROPM(5,MSET)
	
C
	CALL CLEARA (R,NEQ)
      CALL CLEARA (RL,27)


	MMGR = 3
	MMGS = 3
	IF(NNO.EQ.3) MMGS = 1   !SHELL 3 NODE

C     ----------------
C     GAUSS POINT LOOP
C     ----------------
      DO 800  IGR=1,MMGR
      RI = GLOC(IGR,MMGR)
      DO 800  IGS=1,MMGS
      SI = GLOC(IGS,MMGS)
      WT = GWT(IGR,MMGR)*GWT(IGS,MMGS)
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
	IF(NNO.NE.3) THEN
      CALL SHAP2D (RI,SI,H,P,NODEX(1,MEMBA),NNO)
	IF(NNO.EQ.9) CALL SHAP2D9(RI,SI,H,P,NODEX(1,MEMBA),NNO)            !9 NODE ELEMENT
	CALL SHJACO (NNO,XYZ(1,MEMBA),P,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FA)
	ELSEIF(NNO.EQ.3) THEN                                             !SHELL 3 NODE
	CALL GAUSST(RI,SI,TI,WT,IGR,MMGR,0)
	CALL SHAP2D3(RI,SI,H,P,NNO)
	CALL JACO2D3(XYZ(1,MEMBA),P,VR,VS,VT,FA,XJI,DET,MEMBA,NNO)
	ENDIF
C     --------------------------------------------------------
	
C     --------------------------------------------------------
C     UNIFORMLY DISTRIBUTED GLOBAL TRACTION (RT); CONSERVATIVE
C     --------------------------------------------------------
	DO 390  INO=1,NNO
      FAC = H(INO)*WT*DET*DENS*THCK
      DO 380  I=1,3
 380  RL(I,INO) = RL(I,INO) + FAC*GLD(I)
 390  CONTINUE


 800  CONTINUE

	CALL FIXSHE_OFFSHORE (RL,KEG,MEMBA,NNO,3)

C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1             A(LEI),RL,NSF,NNF,4)
C
 700  DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,MEMBA)
 600  IF (IEQ.NE.0) R(IEQ) = R(IEQ) + RL(INF,INO)


	CALL LDASEM_OFFSHORE (R)



 290  IF (IGV.LT.NGV) GOTO 200

 
C
      RETURN
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================

C	=======================================================================
C	================= START OF TETRAHEDRA SURFACE LOAD ====================
C	=======================================================================
	SUBROUTINE SOLTRAG(NGV,XYZ,LM,IGIDM) 
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

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
C	NEXT LINE ADDED BY GILSON - MARCH2004 (LOAD INPUT)
	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================

	COMMON A(9000000),IA(9000000)

      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
	DIMENSION XJ(3,3),XJI(3,3)
	DIMENSION RL(3,NNM),XYZ(NCO*NNM,NELE),LM(NEF,NELE)
	DIMENSION H(NNM),P(3,NNM),XY(3,NNM)
	DIMENSION IGIDM(NELE),TT(3)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)
	  
      ALLOCATE(RAL(NEF),MAL(NEF))
      
C     ---------------
C     INITIALISATION
C     --------------

C	FACE FOR TETRAHEDRA ELEMENT
	NN(4,1:6) = [2,3,4,8,9,10]
	NN(1,1:6) = [1,3,4,6,9,7]
	NN(2,1:6) = [1,2,4,5,10,7]
	NN(3,1:6) = [1,2,3,5,8,6]


      IGV = 0

 200	READ (ITI,*) MEMBA,IFACE,TT(1),TT(2),TT(3),ILCN,ILCC
	IGV = IGV + 1

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	CALL ELEREODER(IGIDM,NELE,MEMBA)

	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO

	CALL CLEARA(RL,3*NNM)

	MG1 = 3
	MG2 = 3
	IF(NNM.EQ.4.OR.NNM.EQ.10) THEN
	MG1 = 3
	IF(NNM.EQ.10) MG1 = 4
	MG2 = 1
	MMM = 3
	IF(NNM.EQ.10) MMM = 6
	ENDIF


	IPT = 0
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	DO 800  II=1,MG1
      DO 800  JJ=1,MG2
	IPT = IPT + 1
	
	IF(NNM.EQ.4.OR.NNM.EQ.10) GOTO 250

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE3D(IFACE,II,JJ,RI,SI,TI,WT,MG2,MG1)
C	SHAPE FUNCTION
	CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO3D(XY,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL SOLIDVECT(IFACE,XJ,DET,VT)



	GOTO 251
250	CONTINUE

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GAUSST(RI,SI,TI,WT,II,MG1,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H2,G,MMM)
C	JACOBIAN
	DO I = 1,3
	COVR(I) = 0.0D0
	COVS(I) = 0.0D0
	DO J = 1,MMM
	MM = NN(IFACE,J)
	COVR(I) = COVR(I) + G(1,J)*XY(I,MM)
	COVS(I) = COVS(I) + G(2,J)*XY(I,MM)
	ENDDO
	ENDDO
	CALL VECPRD(COVR,COVS,VT)
	CALL SCALEN(VT,VT,DET,3)
	IF(IFACE.EQ.2.OR.IFACE.EQ.4) VT = -1.0*VT
C	TRIANGULAR AREA
	DET = 0.5*DET
C	REARRANGE SHAPE FUNCTION
	H(1:NNM) = 0.0D0
	DO I = 1,MMM
	MM = NN(IFACE,I)
	H(MM) = H2(I)
	ENDDO

251	CONTINUE

	
      DO 390  INO=1,NNM
      FAC = H(INO)*WT*DET
      DO 380 I = 1,3
	RL(I,INO) = RL(I,INO) + TT(I)*FAC
 380	CONTINUE
 390  CONTINUE


 800  CONTINUE

	CALL FIXSOL (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE

      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)     
      
 290  IF (IGV.LT.NGV) GOTO 200

      DEALLOCATE(RAL,MAL)	  
      
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE SOLTRAT(NGV,XYZ,LM,IGIDM) 
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

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
C	NEXT LINE ADDED BY GILSON - MARCH2004 (LOAD INPUT)
	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================

	COMMON A(9000000),IA(9000000)

      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
C
	DIMENSION XJ(3,3),XJI(3,3)
	DIMENSION RL(3,NNM),XYZ(NCO*NNM,NELE),LM(NEF,NELE)
	DIMENSION H(NNM),P(3,NNM),XY(3,NNM)
	DIMENSION IGIDM(NELE)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)
	  
      ALLOCATE(RAL(NEF),MAL(NEF))
      
C     ---------------
C     INITIALISATION
C     --------------

C	FACE FOR TETRAHEDRA ELEMENT
	NN(4,1:6) = [2,3,4,8,9,10]
	NN(1,1:6) = [1,3,4,6,9,7]
	NN(2,1:6) = [1,2,4,5,10,7]
	NN(3,1:6) = [1,2,3,5,8,6]


      IGV = 0

 200	READ (ITI,*) MEMBA,IFACE,FLOD,LOHT,WALEV,ILCN,ILCC
	IGV = IGV + 1

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	CALL ELEREODER(IGIDM,NELE,MEMBA)

	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO

	CALL CLEARA(RL,3*NNM)

	MG1 = 3
	MG2 = 3
	IF(NNM.EQ.4.OR.NNM.EQ.10) THEN
	MG1 = 3
	IF(NNM.EQ.10) MG1 = 4
	MG2 = 1
	MMM = 3
	IF(NNM.EQ.10) MMM = 6
	ENDIF

	IPT = 0
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	DO 800  II=1,MG1
      DO 800  JJ=1,MG2
	IPT = IPT + 1
	
	IF(NNM.EQ.4.OR.NNM.EQ.10) GOTO 250

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE3D(IFACE,II,JJ,RI,SI,TI,WT,MG2,MG1)
C	SHAPE FUNCTION
	CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO3D(XY,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL SOLIDVECT(IFACE,XJ,DET,VT)


	GOTO 251
250	CONTINUE

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GAUSST(RI,SI,TI,WT,II,MG1,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H2,G,MMM)
C	JACOBIAN
	DO I = 1,3
	COVR(I) = 0.0D0
	COVS(I) = 0.0D0
	DO J = 1,MMM
	MM = NN(IFACE,J)
	COVR(I) = COVR(I) + G(1,J)*XY(I,MM)
	COVS(I) = COVS(I) + G(2,J)*XY(I,MM)
	ENDDO
	ENDDO
	CALL VECPRD(COVR,COVS,VT)
	CALL SCALEN(VT,VT,DET,3)
	IF(IFACE.EQ.2.OR.IFACE.EQ.4) VT = -1.0*VT
C	TRIANGULAR AREA
	DET = 0.5*DET
C	REARRANGE SHAPE FUNCTION
	H(1:NNM) = 0.0D0
	DO I = 1,MMM
	MM = NN(IFACE,I)
	H(MM) = H2(I)
	ENDDO

251	CONTINUE

	FAC1 = FLOD*WT*DET

	IF(LOHT.EQ.1)  THEN
	GALEV = 0.0D0
	DO INO = 1,NNM
	GALEV = GALEV + H(INO)*XY(NGRAV,INO)
	ENDDO
	IF(GALEV.GT.WALEV) GOTO 800
	DTH = WALEV-GALEV
	FAC1 =-FLOD*DTH*WT*DET
	ENDIF


      DO 390  INO=1,NNM
      FAC = H(INO)*FAC1
      DO 380 I = 1,3
	RL(I,INO) = RL(I,INO) + VT(I)*FAC
 380	CONTINUE
 390  CONTINUE

	
 800  CONTINUE

	CALL FIXSOL (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE

      IEFL = 0
      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RAL(IEFL) + RL(INF,INO)
          MAL(IEFL) = IEQ
      ENDIF
 700  CONTINUE
 
	CALL LDASEM_NEW (RAL,MAL,IEFL)      

 290  IF (IGV.LT.NGV) GOTO 200

      DEALLOCATE(RAL,MAL)	

      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE SOLOFFL(NGV,XYZ,LM,IGIDM) 
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

	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
C	NEXT LINE ADDED BY GILSON - MARCH2004 (LOAD INPUT)
	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV
	
	
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION
	
	COMMON /WARNING/ WARNING,RAMDA(100),RK
      
       COMMON /DIFFRACTIONCOORDINATE / DIFFRACT_X , DIFFRACT_Y , DIFFRACT_GRAVITY 

C	==================================================================

	COMMON A(9000000),IA(9000000)

C
      DIMENSION R(NEQ)
	DIMENSION XJ(3,3),XJI(3,3)
	DIMENSION RL(3,NNM),XYZ(NCO*NNM,NELE),LM(NEF,NELE)
	DIMENSION H(NNM),P(3,NNM),XY(3,NNM)
	DIMENSION IGIDM(NELE)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)
	
C     FOR OFFSHORE LOAD
      DIMENSION VWIND(3),VH1(3),VH2(3),VGV(3),LWCASE(10),TT(3),XYZCP(3),VCURRENTL(5)

C     FOR OFFSHORE LOAD -- PRAMIN 
      DIMENSION FLIST(5),VREF(3),VREW(3),VNOL(3),VNOLG(3)
      
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:)
	
	      	
	IF(NGV.EQ.0) GOTO 1000
		
	PI = 3.141592654	
	
      ! ****************************************************
      ! *                                                  *
      ! *     THE DETAIL SEE IN FRAMEELEMLOAD.FOR          *
      ! *                                                  *
      ! ****************************************************
      
	ALLOCATE(RVAL(21,NGV),IVAL(20,NGV))
	
C     READ INPUT DATA AND STORE ON TEMPORARY VARIABLE
	DO IGV = 1,NGV
        READ(ITI,*) MEMBA,IFACE,NUMCP,NFUNCTION,ROUGH,CD,CM,CSA,DIAM1,SHFAC,GUST_FACTOR,NGROWTH,GROWTH,NIRRE,OFFSELECT_X,LWCASE(1)
     1             ,LWCASE(2),LWCASE(5),LWCASE(6),LWCASE(7),ILCN,ILCC !READ ONLY FIRST TIME STEP
        IVAL(1 ,IGV) = MEMBA
        IVAL(2 ,IGV) = IFACE
        IVAL(3 ,IGV) = NUMCP
        IVAL(4 ,IGV) = LWCASE(1)
        IVAL(5 ,IGV) = LWCASE(2)
        IVAL(6 ,IGV) = LWCASE(3)
        IVAL(7 ,IGV) = LWCASE(4)
        IVAL(8 ,IGV) = LWCASE(5)
        IVAL(9 ,IGV) = LWCASE(6)
        IVAL(10,IGV) = LWCASE(7)
        IVAL(11,IGV) = ILCN
        IVAL(12,IGV) = ILCC
        IVAL(12,IGV) = NIRRE
        
        RVAL(1 ,IGV) = CD
        RVAL(2 ,IGV) = CM
        ! ---- CALCULATE MARINE GROWTH ----
        DIAM2=DIAM1+(GROWTH*2)
        ! ---------------------------------
        RVAL(3 ,IGV) = DIAM2   !NORMINAL DIAMETER
        RVAL(4 ,IGV) = SHFAC  !SHAPE FACTOR = 0.25*PI for circular section
        RVAL(5 ,IGV) = GUST_FACTOR 
        IF (LWCASE(7).NE.1) GUST_FACTOR = 1.0D0 ! PROTECT OTHER CASE =  0.0D0
        RVAL(21,IGV) = OFFSELECT_X  ! CASE SELECT
        RVAL(13,IGV) = NFUNCTION
        RVAL(14,IGV) = ROUGH
        RVAL(15,IGV) = DIAM1
        RVAL(16,IGV) = CSA 
        
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
!60                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( SOLID ELEMENT )----------------')
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
!66                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( SOLID ELEMENT )----------------')
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
         ! --------------------------------------------------------------------------------------------------------------
        
	ENDDO
C     =====================================================================================================
	
101	STOPFUNCTION=0.0D0 ! SET CONDITION PROTECT ERROR EFFECT
	
	IF(KOFFL.EQ.0) THEN
	  DEALLOCATE(RVAL,IVAL)
	  GOTO 1000 !IF NO OFFSHORE LOAd ANALYSIS NO NEED TO DO THE CALCULATION...JUST READ THE INPUT DATA
	ENDIF
  	
  	
	ALLOCATE(RRV(NEQ),RRC(NEQ))
	

      IF(KOFFL.NE.0) THEN !CALL ALL OFFSHORE LOAD PARAMETERS
      !  --- OLD : READING DATA ------------------------------------------------------------
      !  CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,PERIOD,GRAV,RHOW,RHOA,
      ! 1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
      ! 1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED)
      ! -----------------------------------------------------------------------------------
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      ENDIF   
         
      
C	FACE FOR TETRAHEDRA ELEMENT
	NN(4,1:6) = [2,3,4,8,9,10]
	NN(1,1:6) = [1,3,4,6,9,7]
	NN(2,1:6) = [1,2,4,5,10,7]
	NN(3,1:6) = [1,2,3,5,8,6]

      DO 7500 ILCAS = 1,LCS !LOOP OVER LOAD CASE NUBER
      
      DO 7500 ITIME = 1,NTIME !LOOP OVER OFFSHORE TIME STEP
        
      RRV(1:NEQ) = 0.0D0 !INITIALIZE VARY OFFSHORE LOAD
      RRC(1:NEQ) = 0.0D0 !INITIALIZE CONSTANT OFFSHORE LOAD
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')
      CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','READ')
      
      DO 290 IGV = 1,NGV
 
C     ------------------------------
C     CALL INPUT DATA FROM BACKUP
        MEMBA       = IVAL(1 ,IGV)
        IFACE       = IVAL(2 ,IGV)
        NUMCP       = IVAL(3 ,IGV)
        LWCASE(1)   = IVAL(4 ,IGV)
        LWCASE(2)   = IVAL(5 ,IGV)
        LWCASE(3)   = IVAL(6 ,IGV)
        LWCASE(4)   = IVAL(7 ,IGV)
        LWCASE(5)   = IVAL(8 ,IGV)
        LWCASE(6)   = IVAL(9 ,IGV)
        LWCASE(7)   = IVAL(10,IGV)
        ILCN        = IVAL(11,IGV)
        ILCC        = IVAL(12,IGV)
        IORRE       = IVAL(12,IGV)
        
        CD         = RVAL(1 ,IGV)
        CM         = RVAL(2 ,IGV)
        DIAM       = RVAL(3 ,IGV)
        SHFAC      = RVAL(4 ,IGV)
        GUSTFACTOR = RVAL(5 ,IGV)
        NFUNCTION  = RVAL(13,IGV)
        ROUGH      = RVAL(14,IGV)
        DIAM1      = RVAL(15,IGV)
        CSA        = RVAL(16,IGV)  
        IF (NFUNCTION.EQ.2.0D0)THEN
        CD         = RVAL(1 ,IGV)
        CM         = RVAL(2 ,IGV)
        ELSEIF (NFUNCTION.EQ.1.0)THEN
        CD         = 0 
        CM         = 0
        ENDIF
        
C    ==================================================================   
C    ======== Chana Modifile offshore loadcase selecttion
C    ======== 5 / July / 2012 
C    ==================================================================   
      offselect = RVAL(21,IGV) 
      
      IF (IORRE.EQ.1)THEN
      CALL OFFSPARA_READ_COLLECT_DATA
       
      CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,ORDER,PERIOD,GRAV,RHOW,RHOA,
     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,VB,HM,HW,HC,RHIGH,UH,ALPHA,Z0,WVTIME,
     1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED,NCURRENT,POWERLAW,VCURRENTL,
     1                  VCURRENTAPI,UHAPI,NWINDO,AVERAGE,UHD,VCURRENTP,FACTOR,
     1                  WKF,CBF,WFC,TIME)
      WARNING=GRAV   
C    ==================================================================

      ! --- STOP FUNCTION ---
      !IF (STOPFUNCTION.EQ.1D0)THEN
      !STOP
      !ENDIF
      ! --- END STOP FUNCTION ---
      ELSEIF (IORRE.EQ.2)THEN
      CALL SELECTSPECTRUM (offselect,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
      ENDIF
      
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
        
      IF(ILCN.NE.ILCAS.AND.ILCC.NE.ILCAS) GOTO 290


	CALL ELEREODER(IGIDM,NELE,MEMBA)

	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO


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
	!DO INO = 1,NNO
      XMID = XMID + XY(1,INO)
      YMID = YMID + XY(2,INO)
      ZMID = ZMID + XY(3,INO)
	!ENDDO
      XMID = XMID/NNO
      YMID = YMID/NNO
      ZMID = ZMID/NNO
	
!C     ELEMENT MID COORDINATES IN WAVE DIRECTION SYSTEM
!      HM1 = VH1(1)*XMID + VH1(2)*YMID + VH1(3)*ZMID 
!      HMG = VGV(1)*XMID + VGV(2)*YMID + VGV(3)*ZMID 
!      HM2 = VH2(1)*XMID + VH2(2)*YMID + VH2(3)*ZMID 
!C     ----------------------------------------------------	
!    
!C     ----------------------------------------------------
!C     STRUCTURAL ALIGNMENT REFERENCE PROPERTIES
!      CALL OFFREFAXIS(NUMCP,XXR,YYR,ZZR,VREF,'OREF')
!      
!C     STRUCTURAL ALIGNMENT REFERENCE POSITION IN WAVE DIRECTION SYSTEM 
!      HH1 = VH1(1)*XXR + VH1(2)*YYR + VH1(3)*ZZR 
!      HHG = VGV(1)*XXR + VGV(2)*YYR + VGV(3)*ZZR 
!      HH2 = VH2(1)*XXR + VH2(2)*YYR + VH2(3)*ZZR 
!
!C     STRUCTURAL ALIGNMENT REFERENCE VECTOR IN WAVE DIRECTION SYSTEM 
!      VREW(1) = VH1(1)*VREF(1) + VH1(2)*VREF(2) + VH1(3)*VREF(3)
!      VREW(2) = VGV(1)*VREF(1) + VGV(2)*VREF(2) + VGV(3)*VREF(3) 
!      VREW(3) = VH2(1)*VREF(1) + VH2(2)*VREF(2) + VH2(3)*VREF(3) 
!C     ----------------------------------------------------
!
!C     ----------------------------------------------------
!C     POSITION OF COORDINATES REFER TO WAVE REFERNCE COORDINATE (IN WAVE DIRECTION SYSTEM) CORESPONDING TO STRUCTURAL ALIGNMENT
!      HIG = HMG - HHG   !VERTICAL DISTANCE FROM STRUCTURAL REFERENCE COORDINATE
!      RLN = HIG/VREW(2) !LENGTH ALONG ALIGNMENT VECTOR
!      H1M =  HH1 + RLN*VREW(1) - H1REF !H1 DISTANCE REFER TO WAVE REFERNCE COORDINATE
!      HGM =  HHG + RLN*VREW(2) - HGREF !GRAVITY DISTANCE REFER TO WAVE REFERNCE COORDINATE
!      H2M =  HH2 + RLN*VREW(3) - H2REF !H2 DISTANCE REFER TO WAVE REFERNCE COORDINATE

C     ----------------------------------------------------


      H1M = VH1(1)*XMID + VH1(2)*YMID + VH1(3)*ZMID 
      HGM = VGV(1)*XMID + VGV(2)*YMID + VGV(3)*ZMID 
      H2M = VH2(1)*XMID + VH2(2)*YMID + VH2(3)*ZMID
      
C     FIND (OMEGA,WAVENUMBER(RK),RAMDA,A,RATIO) 
      IF (IWAVE == 1 . OR . IWAVE == 6 . OR . IWAVE == 7 ) THEN
          ORDER=0.0D0 ! PROTECT ERROR FROM STREAM FUNCTION
          OMEGA = 2.0*PI/PERIOD      
          CALL NEWTON_RAPHSON(WDEPTH,GRAV,OMEGA,RK)
          RAMDA = 2.0*PI/RK    
          IF (RK.LE.0.0D0)THEN
          WRITE (*,8)
          WRITE (*,30)
          WRITE (*,31)
30        FORMAT ('************** OFFSHORE WARNING MESSAGE ( AIRY WAVE THEORY FOR SOLID ) **************')
31        FORMAT ('- PLEASE CHECK WAVE PERIOD AND THEORY CONDITION ( WAVE NUMBER < 0 )')
          STOP
          ENDIF
          AVAL = 1.0D0
          RATIO = 1.0D0                
      ELSEIF (IWAVE == 2) THEN
          ORDER=0.0D0 ! PROTECT ERROR FROM STREAM FUNCTION      
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
32        FORMAT ('************** OFFSHORE WARNING MESSAGE ( STOKE FIFTH ORDER THEORY FOR SOLID ) **************')
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
      ENDIF    
8     FORMAT ('')              
C	---------------------------------------------------
C	WATER SURFACE LEVEL   
      IF (IORRE.EQ.1) THEN
      IF (IWAVE == 1) THEN
      WNU = 0.50D0*WVHIGHT*COS(RK*H1M - OMEGA*TIME)
      ELSEIF (IWAVE == 2) THEN
      WNU = (1.0D0/RK)*(F1*COS(1.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F2*COS(2.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F3*COS(3.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F4*COS(4.0D0*(RK*H1M - OMEGA*TIME)) + 
     1                  F5*COS(5.0D0*(RK*H1M - OMEGA*TIME))) 
      ELSEIF (IWAVE == 3) THEN
      CALL STREAMWAVEFUNCTION (WVHIGHT,PERIOD,WDEPTH,TIME,
	1                         WNU,GRAV,ORDER)     
          
      ELSEIF (IWAVE.EQ.4)THEN   
    	call Cnoidal_Wave_Crest  (WVHIGHT,WDEPTH,GRAV,PERIOD,TIME,Wave_Crest)
    	call Cnoidal_wave_Length (WVHIGHT,WDEPTH,GRAV,PERIOD,wavenumber,ALAMDA,Ammm ,AKKK ,AEEE)
    	RK=wavenumber
    	RAMDA=ALAMDA
    	WNU = Wave_Crest  
    	ELSEIF (IWAVE.EQ.4)THEN   
      Call Solitary_Wave (WVHIGHT,WDEPTH,GRAV,TIME,0d0,z,u,v,au,av,Surface_Elevation,q)
    	WNU = Surface_Elevation
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
      YL = WDEPTH + WNU 
      
      IF (WNU.GE.0.0D0) YLMIN = WDEPTH - WNU
      IF (WNU.LT.0.0D0) YLMIN = YL
      
C     WATER SURFACE PROFILE MIN
      IF (IWAVE.EQ.1 .OR.IWAVE.EQ.6 )THEN
      YLMIN = WDEPTH - WNU
      ELSEIF (IWAVE.EQ.2.OR.IWAVE.EQ.3.OR.IWAVE.EQ.4.OR.IWAVE.EQ.5)THEN
      YLMIN = (WDEPTH + WNU) - WVHIGHT
      ENDIF
            
      CALL DAMCOEFFICIENT (NFUNCTION,ROUGH,DIAM,PERIOD,WVHIGHT,RK
     1                          ,WDEPTH,TIME,RATIO,OMEGA,GRAV,IWAVE,AVAL,WAVENUMBER
     1                          ,AKKK,X,AMMM,ALAMDA,AEEE,CD,CM)
      
C     MODIFY CM WITH NORMINAL DIAMETER AND SHAPE FACTOR        
      CM = CM*D
      
C     ===================================================	
C     ===================================================	


	CALL CLEARA(RL,3*NNM)

	MG1 = 3 !FOR ACCURACY OF INTEGRATION
	MG2 = 3 !FOR ACCURACY OF INTEGRATION
	IF(NNM.EQ.4.OR.NNM.EQ.10) THEN
	MG1 = 3
	IF(NNM.EQ.10) MG1 = 4
	MG2 = 1
	MMM = 3
	IF(NNM.EQ.10) MMM = 6
	ENDIF

	IPT = 0
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	DO 800  II=1,MG1
      DO 800  JJ=1,MG2
	IPT = IPT + 1
	
	IF(NNM.EQ.4.OR.NNM.EQ.10) GOTO 250 ! 4-NODE OR 10-NODE SOLID ELEMENT

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE3D(IFACE,II,JJ,RI,SI,TI,WT,MG2,MG1)
C	SHAPE FUNCTION
	CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO3D(XY,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL SOLIDVECT(IFACE,XJ,DET,VT)


	GOTO 251
250	CONTINUE

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GAUSST(RI,SI,TI,WT,II,MG1,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H2,G,MMM)
C	JACOBIAN
	DO I = 1,3
	COVR(I) = 0.0D0
	COVS(I) = 0.0D0
	DO J = 1,MMM
	MM = NN(IFACE,J)
	COVR(I) = COVR(I) + G(1,J)*XY(I,MM)
	COVS(I) = COVS(I) + G(2,J)*XY(I,MM)
	ENDDO
	ENDDO
	CALL VECPRD(COVR,COVS,VT)
	CALL SCALEN(VT,VT,DET,3)
	IF(IFACE.EQ.2.OR.IFACE.EQ.4) VT = -1.0*VT
C	TRIANGULAR AREA
	DET = 0.5*DET
C	REARRANGE SHAPE FUNCTION
	H(1:NNM) = 0.0D0
	DO I = 1,MMM
	MM = NN(IFACE,I)
	H(MM) = H2(I)
	ENDDO

251	CONTINUE

      XYZCP(1:3) = 0.0D0
	RLEV = 0.0D0
	DO INO = 1,NNM
	  RLEV = RLEV + H(INO)*XY(NGRAV,INO)
	  DO I = 1,3
	      XYZCP(I) = XYZCP(I) + H(INO)*XY(I,INO)
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

C     Diffrection Coordinate 
      
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
     
      CALL WAVE_PRESSURE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                   VTIDE,VWIND0,WH1,WGV,WH2,VNOL,NCURRENT,H0,VCURRENTP,POWERLAW,LC,VCURRENTL,PERIOD,LWCASE,CS,
     1                   WKF,CBF,AP,SP,
     1                   WFC,YLMIN,IORRE,IRWAVE)
      
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 
      
!      WRITE(57,98765) WFX , HGM , IWAVE , YL      ! TOEY 10/2021
!98765 format(F12.3,3X,F12.3,4x,I3,4X,F12.3)       ! TOEY 10/2021
      
!      if ( IWAVE .EQ. 6 ) WRITE(58,98765) WFX , HGM , IWAVE , YL 
!      if ( IWAVE .NE. 6 ) WRITE(59,98765) WFX , HGM , IWAVE , YL 
      
      

C      WIFX = VH1(1)*WIH1 + VGV(1)*WIGV + VH1(1)*WIH2 
C      WIFY = VH1(2)*WIH1 + VGV(2)*WIGV + VH1(2)*WIH2 
C      WIFZ = VH1(3)*WIH1 + VGV(3)*WIGV + VH1(3)*WIH2 
      
C     TRANSFORM WAVE FORCE VECTOR TO GLOBAL   VNOL  TO VNOLG      
      VNOLG(1) = VH1(1)*VNOL(1) + VGV(1)*VNOL(2) + VH1(1)*VNOL(3)
      VNOLG(2) = VH1(2)*VNOL(1) + VGV(2)*VNOL(2) + VH1(2)*VNOL(3)
      VNOLG(3) = VH1(3)*VNOL(1) + VGV(3)*VNOL(2) + VH1(3)*VNOL(3) 

      FACA = ABS(VT(1)*VNOLG(1)+VT(2)*VNOLG(2)+VT(3)*VNOLG(3))  !PROJECTION AREA IN FORCE DIRECTION
      WFX = 0.5*WFX*FACA !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA
      WFY = 0.5*WFY*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
      WFZ = 0.5*WFZ*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
C     -------------------------- 
      KFLACAL = 1
      
C    ===================================================================================    
C    =============================== UNNECESSARY BY TOEY ===============================  
C      CASE(4) !WATER LEVEL LOAD   
C          WLEV = RLEV - SEABED
C          IF(WLEV.LT.0.0D0) GOTO 407
C          IF(WLEV.GT.PEAKWLEV) GOTO 407
C      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
C     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
C     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF)
C      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH2(1)*WH2 
C      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH2(2)*WH2 
C      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH2(3)*WH2 
C    ===================================================================================

      CASE(5,6) !WAVE BREAKING LOAD  PLUNGING 
          WLEV = RLEV - (SEABED + WDEPTH)
          IF(WLEV.GT.0.0D0) GOTO 407
          
      CALL WAVE_PRESSURE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                   VTIDE,VWIND0,WH1,WGV,WH2,VNOL,NCURRENT,H0,VCURRENTP,POWERLAW,LC,VCURRENTL,PERIOD,LWCASE,CS,
     1                   WKF,CBF,AP,SP,
     1                   WFC,YLMIN,IORRE,IRWAVE)
      
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 

     
      
      WIFX = VH1(1)*WIH1 + VGV(1)*WIGV + VH1(1)*WIH2 
      WIFY = VH1(2)*WIH1 + VGV(2)*WIGV + VH1(2)*WIH2 
      WIFZ = VH1(3)*WIH1 + VGV(3)*WIGV + VH1(3)*WIH2 
      
C     TRANSFORM WAVE FORCE VECTOR TO GLOBAL   VNOL  TO VNOLG      
      VNOLG(1) = VH1(1)*VNOL(1) + VGV(1)*VNOL(2) + VH1(1)*VNOL(3)
      VNOLG(2) = VH1(2)*VNOL(1) + VGV(2)*VNOL(2) + VH1(2)*VNOL(3)
      VNOLG(3) = VH1(3)*VNOL(1) + VGV(3)*VNOL(2) + VH1(3)*VNOL(3) 

      FACA = ABS(VT(1)*VNOLG(1)+VT(2)*VNOLG(2)+VT(3)*VNOLG(3))  !PROJECTION AREA IN FORCE DIRECTION
      WFX = 0.5*WFX*FACA !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA
      WFY = 0.5*WFY*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
      WFZ = 0.5*WFZ*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 

C    ===================================================================================    
C    =============================== UNNECESSARY BY TOEY ===============================       
C      CASE(6) !WAVE BREAKING LOAD  SURING
C          WLEV = RLEV - SEABED
C          IF(WLEV.LT.0.0D0) GOTO 407
C          IF(WLEV.GT.PEAKWLEV) GOTO 407
C      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
C     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
C     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
C      WFX = VH1(1)*WFF ; WFY = VH1(2)*WFF ; WFZ = VH1(3)*WFF ;
C    ===================================================================================
     
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
      WFX = VWIND(1)*WFF*CSA*SHFAC*GUSTFACTOR ; WFY = VWIND(2)*WFF*CSA*SHFAC*GUSTFACTOR ; WFZ = VWIND(3)*WFF*CSA*SHFAC*GUSTFACTOR ;
      ENDSELECT
      
	TT(1) = TT(1) + WFX
	TT(2) = TT(2) + WFY
	TT(3) = TT(3) + WFZ
	
407   CONTINUE
C     -------------------------- 


      DO 390  INO=1,NNM
      FAC = H(INO)*WT*DET
      DO 380 I = 1,3
	RL(I,INO) = RL(I,INO) + TT(I)*FAC
 380	CONTINUE
 390  CONTINUE

	
 800  CONTINUE

C     MUST BE ADDED LATER  FOR SOLID FIXEND FORCED
C	CALL FIXSOL (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE


      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + RL(INF,INO) !VARY LOAD
      IF (IEQ.NE.0.AND.ILCC.GT.0) RRC(IEQ) = RRC(IEQ) + RL(INF,INO) !CONSTANT LOAD
 700  CONTINUE


 290  CONTINUE


      IF(KOFFL.NE.0) THEN
        CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','WRIT') !STORE VARY LOAD
        CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','WRIT') !STORE CONSTANT LOAD
      ENDIF      

7500  CONTINUE

	DEALLOCATE(RRV,RRC,RVAL,IVAL)

1000  RETURN
C
      END
      
C	=======================================================================
C	=======================================================================
      SUBROUTINE SOLOFFLRES(NGV,XYZ,LM,IGIDM) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     READS AND GENERATES NODAL LOADS, ADDS CONTRIBUTIONS TO LOAD V
C     -------------------------------------------------------------
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

	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
C	NEXT LINE ADDED BY GILSON - MARCH2004 (LOAD INPUT)
	COMMON /LCSS/ ILCN,ILCC

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV
	
	
	COMMON /offshoreselectx_data_correction/ offselect,NUM_OF_OFFSHORE_PARAMETER
	
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION
	
	COMMON / EIGVPED / EIGVFREQ(1000),EIGVPERI(1000)
      COMMON / STOREMODE / RVECT(100000,100),NMOD
	
	COMMON /WARNING/ WARNING,RAMDA(100),RK

C	==================================================================

	COMMON A(9000000),IA(9000000)

C
      DIMENSION R(NEQ)
	DIMENSION XJ(3,3),XJI(3,3)
	DIMENSION RL(3,NNM),XYZ(NCO*NNM,NELE),LM(NEF,NELE)
	DIMENSION H(NNM),P(3,NNM),XY(3,NNM)
	DIMENSION IGIDM(NELE)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)
	
C     FOR OFFSHORE LOAD
      DIMENSION VWIND(3),VH1(3),VH2(3),VGV(3),LWCASE(10),TT(3),XYZCP(3),VCURRENTL(5)

C     FOR OFFSHORE LOAD -- PRAMIN 
      DIMENSION FLIST(5),VREF(3),VREW(3),VNOL(3),VNOLG(3)
      
	ALLOCATABLE RRV(:),RRC(:),RVAL(:,:),IVAL(:,:)
	
	      	
	IF(NGV.EQ.0) GOTO 1000
		
	PI = 3.141592654	
	
      ! ****************************************************
      ! *                                                  *
      ! *     THE DETAIL SEE IN FRAMEELEMLOAD.FOR          *
      ! *                                                  *
      ! ****************************************************
      
	ALLOCATE(RVAL(21,NGV),IVAL(20,NGV))
	
C     READ INPUT DATA AND STORE ON TEMPORARY VARIABLE
	DO IGV = 1,NGV
        READ(ITI,*) MEMBA,IFACE,NUMCP,NFUNCTION,ROUGH,CD,CM,CSA,DIAM1,SHFAC,NGROWTH,GROWTH,OFFSELECT_X,LWCASE(1:4),ILCN,ILCC !READ ONLY FIRST TIME STEP
        LWCASE(5)    = LWCASE(2)
        LWCASE(6)    = LWCASE(3)
        LWCASE(7)    = LWCASE(4)
        IVAL(1 ,IGV) = MEMBA
        IVAL(2 ,IGV) = IFACE
        IVAL(3 ,IGV) = NUMCP
        IVAL(4 ,IGV) = LWCASE(1)
        IVAL(5 ,IGV) = LWCASE(2)
        IVAL(6 ,IGV) = 0.0
        IVAL(7 ,IGV) = 0.0
        IVAL(8 ,IGV) = LWCASE(5)
        IVAL(9 ,IGV) = LWCASE(6)
        IVAL(10,IGV) = LWCASE(7)
        IVAL(11,IGV) = ILCN
        IVAL(12,IGV) = ILCC
        
        RVAL(1 ,IGV) = CD
        RVAL(2 ,IGV) = CM
        ! ---- CALCULATE MARINE GROWTH ----
        DIAM2=DIAM1+(GROWTH*2)
        ! ---------------------------------
        RVAL(3 ,IGV) = DIAM2   !NORMINAL DIAMETER
        RVAL(4 ,IGV) = SHFAC  !SHAPE FACTOR = 0.25*PI for circular section
        RVAL(21,IGV) = OFFSELECT_X  ! CASE SELECT
        RVAL(13,IGV) = NFUNCTION
        RVAL(14,IGV) = ROUGH
        RVAL(15,IGV) = DIAM1
        RVAL(16,IGV) = CSA 
        
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
60                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( SOLID ELEMENT )----------------')
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
66                 FORMAT ('---------------- OFFSHORE WARNING MESSAGE ( SOLID ELEMENT )----------------')
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
         ! --------------------------------------------------------------------------------------------------------------
        
	ENDDO
C     =====================================================================================================
	
101	STOPFUNCTION=0.0D0 ! SET CONDITION PROTECT ERROR EFFECT
	
	IF(KOFFL.EQ.0) THEN
	  DEALLOCATE(RVAL,IVAL)
	  GOTO 1000 !IF NO OFFSHORE LOAd ANALYSIS NO NEED TO DO THE CALCULATION...JUST READ THE INPUT DATA
	ENDIF
  	
  	
	ALLOCATE(RRV(NEQ),RRC(NEQ))
	

      IF(KOFFL.NE.0) THEN !CALL ALL OFFSHORE LOAD PARAMETERS
      !  --- OLD : READING DATA ------------------------------------------------------------
      !  CALL OFFSPARA_CALL (WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,IWAVE,PERIOD,GRAV,RHOW,RHOA,
      ! 1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
      ! 1                  VGV,VH1,VH2,VWIND,PEAKWLEV,SEABED)
      ! -----------------------------------------------------------------------------------
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALT') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      ENDIF   
         
      
C	FACE FOR TETRAHEDRA ELEMENT
	NN(4,1:6) = [2,3,4,8,9,10]
	NN(1,1:6) = [1,3,4,6,9,7]
	NN(2,1:6) = [1,2,4,5,10,7]
	NN(3,1:6) = [1,2,3,5,8,6]

      DO 7500 ILCAS = 1,LCS !LOOP OVER LOAD CASE NUBER
      
      DO 7500 ITIME = 1,NTIME !LOOP OVER OFFSHORE TIME STEP
        
      RRV(1:NEQ) = 0.0D0 !INITIALIZE VARY OFFSHORE LOAD
      RRC(1:NEQ) = 0.0D0 !INITIALIZE CONSTANT OFFSHORE LOAD
      CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL') !CALL NTIME (NUMBER OF TIME STEP FOR OFFSHORE LOAD GENERATION)
      CALL OFFSHFORC(RRV,ILCAS,ITIME,NEQ,'VARY','READ')
      CALL OFFSHFORC(RRC,ILCAS,ITIME,NEQ,'CONT','READ')
     
      IF (ISOLOP.EQ.1.AND.ICONTROLSPEC.EQ.1) ISPEC = 1
      
      DO 7501 ISPEC = 1,NMOD
      SNA = EIGVFREQ(ISPEC)
      
      DO 290 IGV = 1,NGV
 
C     ------------------------------
C     CALL INPUT DATA FROM BACKUP
        MEMBA       = IVAL(1 ,IGV)
        IFACE       = IVAL(2 ,IGV)
        NUMCP       = IVAL(3 ,IGV)
        LWCASE(1)   = IVAL(4 ,IGV)
        LWCASE(2)   = IVAL(5 ,IGV)
        LWCASE(3)   = IVAL(6 ,IGV)
        LWCASE(4)   = IVAL(7 ,IGV)
        LWCASE(5)   = IVAL(8 ,IGV)
        LWCASE(6)   = IVAL(9 ,IGV)
        LWCASE(7)   = IVAL(10,IGV)
        ILCN        = IVAL(11,IGV)
        ILCC        = IVAL(12,IGV)
        
        CD         = RVAL(1 ,IGV)
        CM         = RVAL(2 ,IGV)
        DIAM       = RVAL(3 ,IGV)
        SHFAC      = RVAL(4 ,IGV)
        NFUNCTION  = RVAL(13,IGV)
        ROUGH      = RVAL(14,IGV)
        DIAM1      = RVAL(15,IGV)
        CSA        = RVAL(16,IGV)  
        offselect  = RVAL(21,IGV) 
        IF (NFUNCTION.EQ.2.0D0)THEN
        CD         = RVAL(1 ,IGV)
        CM         = RVAL(2 ,IGV)
        ELSEIF (NFUNCTION.EQ.1.0)THEN
        CD         = 0 
        CM         = 0
        ENDIF    
       
       
       CALL SELECTSPECTRUM (OFFSELECT,SEABED,WVHIGHT,WDEPTH,H1POS,H2POS,IRWAVE,PERIOD,GRAV,RHOW,RHOA,WKF,WFC,
     1                      FREQ,VGV,VH1,VH2,VWIND,SDG,NSWIND,UHD,ALPHA,Z0,RHIGH,FWIND,TAP,UCURRENT,WVH1,WVH2,
     1                      PER1,PER2,SHAPE1,SHAPE2)
       IF (ISOLOP.EQ.1.AND.ICONTROLSPEC.EQ.1)THEN
       FREQ =  TIME
       ENDIF
      WARNING=GRAV   
C    ==================================================================

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
        
      IF(ILCN.NE.ILCAS.AND.ILCC.NE.ILCAS) GOTO 290


	CALL ELEREODER(IGIDM,NELE,MEMBA)

	KKAK = 0
	DO IKAK = 1,NNM
	DO JKAK = 1,3
	KKAK = KKAK + 1
	XY(JKAK,IKAK) = XYZ(KKAK,MEMBA)
	ENDDO
	ENDDO


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
        XMID = XMID + XY(1,INO)
        YMID = YMID + XY(2,INO)
        ZMID = ZMID + XY(3,INO)
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
      CM = CM*SHFAC*DIAM
      
C     ===================================================	
C     ===================================================	


	CALL CLEARA(RL,3*NNM)

	MG1 = 3 !FOR ACCURACY OF INTEGRATION
	MG2 = 3 !FOR ACCURACY OF INTEGRATION
	IF(NNM.EQ.4.OR.NNM.EQ.10) THEN
	MG1 = 3
	IF(NNM.EQ.10) MG1 = 4
	MG2 = 1
	MMM = 3
	IF(NNM.EQ.10) MMM = 6
	ENDIF

	IPT = 0
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
	DO 800  II=1,MG1
      DO 800  JJ=1,MG2
	IPT = IPT + 1
	
	IF(NNM.EQ.4.OR.NNM.EQ.10) GOTO 250 ! 4-NODE OR 10-NODE SOLID ELEMENT

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE3D(IFACE,II,JJ,RI,SI,TI,WT,MG2,MG1)
C	SHAPE FUNCTION
	CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO3D(XY,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL SOLIDVECT(IFACE,XJ,DET,VT)


	GOTO 251
250	CONTINUE

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GAUSST(RI,SI,TI,WT,II,MG1,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H2,G,MMM)
C	JACOBIAN
	DO I = 1,3
	COVR(I) = 0.0D0
	COVS(I) = 0.0D0
	DO J = 1,MMM
	MM = NN(IFACE,J)
	COVR(I) = COVR(I) + G(1,J)*XY(I,MM)
	COVS(I) = COVS(I) + G(2,J)*XY(I,MM)
	ENDDO
	ENDDO
	CALL VECPRD(COVR,COVS,VT)
	CALL SCALEN(VT,VT,DET,3)
	IF(IFACE.EQ.2.OR.IFACE.EQ.4) VT = -1.0*VT
C	TRIANGULAR AREA
	DET = 0.5*DET
C	REARRANGE SHAPE FUNCTION
	H(1:NNM) = 0.0D0
	DO I = 1,MMM
	MM = NN(IFACE,I)
	H(MM) = H2(I)
	ENDDO

251	CONTINUE

      XYZCP(1:3) = 0.0D0
	RLEV = 0.0D0
	DO INO = 1,NNM
	  RLEV = RLEV + H(INO)*XY(NGRAV,INO)
	  DO I = 1,3
	      XYZCP(I) = XYZCP(I) + H(INO)*XY(I,INO)
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
          
C          CALL WAVE_PRESSURE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
C     1                   VTIDE,VWIND0,WH1,WGV,WH2,VNOL,NCURRENT,H0,VCURRENTP,POWERLAW,LC,VCURRENTL,PERIOD,LWCASE,CS,
C     1                   WKF,CBF,AP,SP,
C     1                   WFC,YLMIN,IORRE,IRWAVE)
     
         IF (IRWAVE.EQ.1)THEN
         CALL  Pierson_Moskowitz_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                         WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.2)THEN
         CALL Jonswap_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                               WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.3)THEN
         CALL Ochi_SPectrum_Solid(WVH1,PER1,SHAPE1,WVH2,PER2,SHAPE2,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                            WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.4)THEN
         CALL  Bretschneider_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                      WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.5)THEN
         CALL TMA_SPectrum_Solid(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                           WH1,WGV,WH2,TAP,UCURRENT)
         ELSEIF (IRWAVE.EQ.6)THEN
         CALL User_Define_SPectrum_SOILD(WVHIGHT,PERIOD,GRAV,WDEPTH,RHOW,CD,CM,DIAM,HGM,FREQ,SNA,SDG,VNOL,VREW,DIAM,
     1                                    WH1,WGV,WH2,TAP,UCURRENT)
         ENDIF
         
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 

      WIFX = VH1(1)*WIH1 + VGV(1)*WIGV + VH1(1)*WIH2 
      WIFY = VH1(2)*WIH1 + VGV(2)*WIGV + VH1(2)*WIH2 
      WIFZ = VH1(3)*WIH1 + VGV(3)*WIGV + VH1(3)*WIH2 
      
C     TRANSFORM WAVE FORCE VECTOR TO GLOBAL   VNOL  TO VNOLG      
      VNOLG(1) = VH1(1)*VNOL(1) + VGV(1)*VNOL(2) + VH1(1)*VNOL(3)
      VNOLG(2) = VH1(2)*VNOL(1) + VGV(2)*VNOL(2) + VH1(2)*VNOL(3)
      VNOLG(3) = VH1(3)*VNOL(1) + VGV(3)*VNOL(2) + VH1(3)*VNOL(3) 

      FACA = ABS(VT(1)*VNOLG(1)+VT(2)*VNOLG(2)+VT(3)*VNOLG(3))  !PROJECTION AREA IN FORCE DIRECTION
      WFX = 0.5*WFX*FACA !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA
      WFY = 0.5*WFY*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
      WFZ = 0.5*WFZ*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
C     -------------------------- 
      KFLACAL = 1
      
C    ===================================================================================    
C    =============================== UNNECESSARY BY TOEY ===============================  
C      CASE(4) !WATER LEVEL LOAD   
C          WLEV = RLEV - SEABED
C          IF(WLEV.LT.0.0D0) GOTO 407
C          IF(WLEV.GT.PEAKWLEV) GOTO 407
C      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
C     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
C     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF)
C      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH2(1)*WH2 
C      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH2(2)*WH2 
C      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH2(3)*WH2 
C    ===================================================================================

      CASE(5,6) !WAVE BREAKING LOAD  PLUNGING 
          WLEV = RLEV - SEABED
          IF(WLEV.LT.0.0D0) GOTO 407
          IF(WLEV.GT.YL) GOTO 407
      CALL WAVE_PRESSURE(OMEGA,RATIO,WVHIGHT,WDEPTH,RK,RHOW,CM,CD,DIAM,H1R,HGM,TIME,IWAVE,ORDER,VREW,AVAL,GRAV,
     1                   VTIDE,VWIND0,WH1,WGV,WH2,VNOL,NCURRENT,H0,VCURRENTP,POWERLAW,LC,VCURRENTL,PERIOD,LWCASE,CS,
     1                   WKF,CBF,AP,SP,
     1                   WFC,YLMIN,IORRE,IRWAVE)
      
      WFX = VH1(1)*WH1 + VGV(1)*WGV + VH1(1)*WH2 
      WFY = VH1(2)*WH1 + VGV(2)*WGV + VH1(2)*WH2 
      WFZ = VH1(3)*WH1 + VGV(3)*WGV + VH1(3)*WH2 

      WIFX = VH1(1)*WIH1 + VGV(1)*WIGV + VH1(1)*WIH2 
      WIFY = VH1(2)*WIH1 + VGV(2)*WIGV + VH1(2)*WIH2 
      WIFZ = VH1(3)*WIH1 + VGV(3)*WIGV + VH1(3)*WIH2 
      
C     TRANSFORM WAVE FORCE VECTOR TO GLOBAL   VNOL  TO VNOLG      
      VNOLG(1) = VH1(1)*VNOL(1) + VGV(1)*VNOL(2) + VH1(1)*VNOL(3)
      VNOLG(2) = VH1(2)*VNOL(1) + VGV(2)*VNOL(2) + VH1(2)*VNOL(3)
      VNOLG(3) = VH1(3)*VNOL(1) + VGV(3)*VNOL(2) + VH1(3)*VNOL(3) 

      FACA = ABS(VT(1)*VNOLG(1)+VT(2)*VNOLG(2)+VT(3)*VNOLG(3))  !PROJECTION AREA IN FORCE DIRECTION
      WFX = 0.5*WFX*FACA !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA
      WFY = 0.5*WFY*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 
      WFZ = 0.5*WFZ*FACA  !0.5 FOR 2 SIDE OF FORCE ATTACKING AREA 

C    ===================================================================================    
C    =============================== UNNECESSARY BY TOEY ===============================       
C      CASE(6) !WAVE BREAKING LOAD  SURING
C          WLEV = RLEV - SEABED
C          IF(WLEV.LT.0.0D0) GOTO 407
C          IF(WLEV.GT.PEAKWLEV) GOTO 407
C      CALL WAVE_LOADING(WVHIGHT,WDEPTH,THIGHT,H1POS,H2POS,RAMDA,GRAV,RHOW,RHOA,
C     1                  WVZETA,VTIDE,VWIND0,H0,AP,SP,CS,RHIGH,UH,ALPHA,Z0,WVTIME,
C     2                  CD,CM,DIAM,LCASE,WLEV,WH1,WGV,WH2,WFF) 
C      WFX = VH1(1)*WFF ; WFY = VH1(2)*WFF ; WFZ = VH1(3)*WFF ;
C    ===================================================================================
     
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
      WFX = VWIND(1)*WFF*CSA ; WFY = VWIND(2)*WFF*CSA ; WFZ = VWIND(3)*WFF*CSA ;
      ENDSELECT
      
	TT(1) = TT(1) + WFX
	TT(2) = TT(2) + WFY
	TT(3) = TT(3) + WFZ
	
407   CONTINUE
C     -------------------------- 


      DO 390  INO=1,NNM
      FAC = H(INO)*WT*DET
      DO 380 I = 1,3
	RL(I,INO) = RL(I,INO) + TT(I)*FAC
 380	CONTINUE
 390  CONTINUE

	
 800  CONTINUE

C     MUST BE ADDED LATER  FOR SOLID FIXEND FORCED
C	CALL FIXSOL (RL,KEG,MEMBA,NEF,NNM)

	IF (NLS.EQ.0) GOTO 711
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEMBA),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,5)
711	CONTINUE


      DO 700  INO=1,NNM
      II = (INO-1)*NNF
      DO 700  INF=1,3
      IEQ = LM(II+INF,MEMBA)
      IF (IEQ.NE.0.AND.ILCN.GT.0) RRV(IEQ) = RRV(IEQ) + RL(INF,INO) !VARY LOAD
      IF (IEQ.NE.0.AND.ILCC.GT.0) RRC(IEQ) = RRC(IEQ) + RL(INF,INO) !CONSTANT LOAD
 700  CONTINUE


 290  CONTINUE


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

1000  RETURN
C
      END
C     =======================================================================
C     =======================================================================
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SOLIDVECT(NFACE,XJ,DET,VC)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION XJ(3,3),VA(3),VB(3),VC(3)

	VA(1:3) = 0.0D0
	VB(1:3) = 0.0D0
	VC(1:3) = 0.0D0

	SELECT CASE(NFACE)

	CASE(2,4)  !5,3

C	PARALLEL VECTOR COMPONENT
C	XJ(2,1),XJ(2,2),XJ(2,3) - ETA DIRECTION
C	XJ(3,1),XJ(3,2),XJ(3,3) - ZETA DIRECTION

	VA(1) = XJ(2,1)
	VA(2) = XJ(2,2)
	VA(3) = XJ(2,3)

	VB(1) = XJ(3,1)
	VB(2) = XJ(3,2)
	VB(3) = XJ(3,3)


	CASE(3,5)  !2,4

C	PARALLEL VECTOR COMPONENT
C	XJ(1,1),XJ(1,2),XJ(1,3) - XI DIRECTION
C	XJ(3,1),XJ(3,2),XJ(3,3) - ZETA DIRECTION

	VA(1) = XJ(1,1)
	VA(2) = XJ(1,2)
	VA(3) = XJ(1,3)

	VB(1) = XJ(3,1)
	VB(2) = XJ(3,2)
	VB(3) = XJ(3,3)

	CASE(6,1)  !6,1

C	PARALLEL VECTOR COMPONENT
C	XJ(1,1),XJ(1,2),XJ(1,3) - XI DIRECTION
C	XJ(3,1),XJ(3,2),XJ(3,3) - ETA DIRECTION

	VA(1) = XJ(1,1)
	VA(2) = XJ(1,2)
	VA(3) = XJ(1,3)

	VB(1) = XJ(2,1)
	VB(2) = XJ(2,2)
	VB(3) = XJ(2,3)
	
	END SELECT


	SELECT CASE(NFACE)

	CASE(5,4,6)
	CALL VECPRD (VA,VB,VC)	

	CASE(3,2,1)
	CALL VECPRD (VB,VA,VC)

	END SELECT

	CALL SCALEN (VC,VC,DET,3)



	RETURN

	END

C	=====================================================================
C	=====================================================================	
C	=====================================================================
