C	=======================================================================
C	======================START OF 8 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE HEATLOD3D(ID,XYZ,LM,PROPM,MTSET,IGIDM,TAMBT,
	1					 HSOC,HCON,TLAMT,NLL,IND) 
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
      ALLOCATABLE RAL(:)
      ALLOCATABLE MAL(:)
      
	PARAMETER(NFACE=6)
	DIMENSION HSOC(NELE*NPT,LCS),HCON(NELE*NFACE,LCS)
	DIMENSION TLAMT(NELE*NFACE,LCS)
	DIMENSION HSOW(NPT),HCOW(NFACE),TLAMW(NFACE),TAMBT(2,LCS)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE),PROPM(NMP,1),MTSET(1)

	CALL CLEARA(HSOC,NELE*NPT*LCS)
	CALL CLEARA(HCON,NELE*NFACE*LCS)
	CALL CLEARI(TLAMT,NELE*NFACE*LCS)


	NHSOC = 0
	NHCON = 0
	NHAMT = 0
	SELECTCASE (IND)
	CASE(1) 
	NHSOC = NLL
	CASE(2) 
	NHCON = NLL
	CASE(3) 
	NHAMT = NLL
      ENDSELECT
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHSOC
	READ(ITI,*) MEMBA,HSO,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	HSOW(1:NPT) = HSO
	CALL ELEREODER(IGIDMEM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL HEATSTOR(HSOC(1,ILCN),HSOW,NELE,NPT,MEMBA)
	CALL HEATL3D1(XYZ(1,MEMBA),PROPM(1,MSET),HSOW,RL)
	CALL FIXSOL(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHCON
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	HCOW = 0.0D0
	HCOW(IFACE) = DUM1
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL HEATSTOR(HCON(1,ILCN),HCOW,NELE,NFACE,MEMBA)
	CALL HEATL3D2(XYZ(1,MEMBA),PROPM(1,MSET),HCOW,RL)
	CALL FIXSOL(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)	  
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHAMT
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	TLAMW = 0.0D0
	TLAMW(IFACE) = DUM1
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	TMP  = TAMBT(1,ILCN)
	LFG  = TAMBT(2,ILCN)
	IF(LFG.EQ.0) GOTO 100
	CALL HEATSTOR(TLAMT(1,ILCN),TLAMW,NELE,NFACE,MEMBA)
	CALL HEATL3D3(XYZ(1,MEMBA),PROPM(1,MSET),TLAMW,TMP,RL)
	CALL FIXSOL(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
100	CONTINUE	
	ENDDO
      
      DEALLOCATE(RAL,MAL)	 
C	==================================================================
	

C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL3D1(COORD,PROPM,HSOC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),HSOC(NPT)
	DIMENSION NODEX(13)
	DIMENSION H(NNM),P(NCO,NNM),XJ(NCO,NCO),XJI(NCO,NCO)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

	NODEX(1:13) = [9,10,11,12,13,14,15,16,17,18,19,20,21]
C	=====================================================

C	==========================================
C	SETTING ELEMENT GAUSSIAN INTEGRATION POINT
C	STANDARD GAUSS INTEGRATION FOR 8 NODES SOLID ELEMENT
C	(MGR,MGS,MGT) = (2,2,2)
C	=======================
      MGR = 2
      MGS = 2
	MGT = 2

	IPT = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS  FOR HEAT SOURCE (CAN MOVE TO ELEMENT LOOP IF CHANGE WITH THE TIME)
C	=========================================
      DO 900  IGR=1,MGR
      RI = GLOC(IGR,MGR)
	DO 900  IGS=1,MGS
	SI = GLOC(IGS,MGS)
	DO 900  IGT=1,MGT
	TI = GLOC(IGT,MGT)
	IPT = IPT + 1
	WT = GWT(IGR,MGR)*GWT(IGS,MGS)*GWT(IGT,MGT)
C     ====================================
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ===================================================
	CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNM)
	CALL JACO3D(COORD,P,XJ,XJI,DET,MEL,NNM)

	DVOL = WT*DET
C	========================================
C	COMPUTE HEAT SOURCE LOAD VECTOR
C	========================================

	DO I  = 1,NNM
	RH(I) = RH(I) + HSOC(IPT)*H(I)*DVOL
	ENDDO

 900  CONTINUE

C	==========================================


C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL3D2(COORD,PROPM,QVC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT



	PARAMETER (NFACE=6)

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),QVC(NFACE)
	DIMENSION NODEX(13)
	DIMENSION H(NNM),P(NCO,NNM),XJ(NCO,NCO),XJI(NCO,NCO)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF


	RH = 0.0D0
	REF= 0.0D0

	NODEX(1:13) = [9,10,11,12,13,14,15,16,17,18,19,20,21]
C	=====================================================

	MG1   = 2
	MG2   = 2

	DO 1000 IFACE = 1,NFACE

	QVCF = QVC(IFACE)

	DO 1000 II = 1,MG1
	DO 1000 JJ = 1,MG2

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE3D(IFACE,II,JJ,RI,SI,TI,WT,MG1,MG2)
C	SHAPE FUNCTION
	CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO3D(COORD,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL SOLIDVECTOR(IFACE,XJ,DET)

	DFACE = WT*DET
	
	DO I  = 1,NNM
	RH(I) = RH(I) + QVCF*H(I)*DFACE
	ENDDO

1000	CONTINUE

C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL3D3(COORD,PROPM,TCD,TMP,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	PARAMETER (NFACE=6)

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),TCD(NFACE)
	DIMENSION NODEX(13)
	DIMENSION H(NNM),P(NCO,NNM),XJ(NCO,NCO),XJI(NCO,NCO)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

	NODEX(1:13) = [9,10,11,12,13,14,15,16,17,18,19,20,21]
C	=====================================================
C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

	MG1   = 2
	MG2   = 2

	DO 1000 IFACE = 1,NFACE
	FAC = TCD(IFACE)
	TCD(IFACE) = 1.0D0

	DO 1000 II = 1,MG1
	DO 1000 JJ = 1,MG2

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE3D(IFACE,II,JJ,RI,SI,TI,WT,MG1,MG2)

C	SHAPE FUNCTION
	CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO3D(COORD,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL SOLIDVECTOR(IFACE,XJ,DET)

	DFACE = WT*DET

	DO I  = 1,NNM
	RH(I) = RH(I) + FAC*HCV*TMP*H(I)*DFACE
	ENDDO
	

1000	CONTINUE

C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================







C	=======================================================================
C	======================START OF 4 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE HEATLOD2D(ID,XYZ,LM,PROPM,MTSET,PROPG,IGSET,
	1					 IGIDM,TAMBT,HSOC,HCON,TLAMT,NLL,IND) 
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
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
      
	PARAMETER(NFACE=6)
	DIMENSION HSOC(NELE*NPT,LCS),HCON(NELE*NFACE,LCS)
	DIMENSION TLAMT(NELE*NFACE,LCS)
	DIMENSION HSOW(NPT),HCOW(NFACE),TLAMW(NFACE),TAMBT(2,LCS)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(1)

	CALL CLEARA(HSOC,NELE*NPT*LCS)
	CALL CLEARA(HCON,NELE*NFACE*LCS)
	CALL CLEARI(TLAMT,NELE*NFACE*LCS)


	NHSOC = 0
	NHCON = 0
	NHAMT = 0


	SELECTCASE (IND)
	CASE(1) 
	NHSOC = NLL
	CASE(2) 
	NHCON = NLL
	CASE(3) 
	NHCON = NLL
	CASE(4) 
	NHAMT = NLL
	CASE(5) 
	NHAMT = NLL
	ENDSELECT
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHSOC
	READ(ITI,*) MEMBA,HSO,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	HSOW(1:NPT) = HSO
	CALL ELEREODER(IGIDMEM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL HEATSTOR(HSOC(1,ILCN),HSOW,NELE,NPT,MEMBA)
	CALL HEATL2D1(XYZ(1,MEMBA),PROPM(1,MSET),HSOW,RL)
	CALL FIXMRE(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)	 
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHCON

	SELECTCASE (IND)
	CASE(2) 
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC
	HCOW = 0.0D0
	HCOW(IFACE) = DUM1
	CASE(3) 
	READ(ITI,*) MEMBA,DUM1,ILCN,ILCC
	HCOW = 0.0D0
	HCOW(5) = DUM1
	ENDSELECT

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	NSET = IGSET(MEMBA)
	CALL HEATSTOR(HCON(1,ILCN),HCOW,NELE,NFACE,MEMBA)
	CALL HEATL2D2(XYZ(1,MEMBA),PROPM(1,MSET),PROPG(1,NSET),
	1			  HCOW,RL)
	CALL FIXMRE(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)	
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHAMT

	SELECTCASE (IND)
	CASE(4) 
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC
	TLAMW = 0.0D0
	TLAMW(IFACE) = DUM1
	CASE(5) 
	READ(ITI,*) MEMBA,DUM1,ILCN,ILCC
	TLAMW = 0.0D0
	TLAMW(5) = DUM1
	ENDSELECT

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	NSET = IGSET(MEMBA)
	TMP  = TAMBT(1,ILCN)
	LFG  = TAMBT(2,ILCN)
	IF(LFG.EQ.0) GOTO 100
	CALL HEATSTOR(TLAMT(1,ILCN),TLAMW,NELE,NFACE,MEMBA)
	CALL HEATL2D3(XYZ(1,MEMBA),PROPM(1,MSET),PROPG(1,NSET),
	1				  TLAMW,TMP,RL)
	CALL FIXMRE(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	 
100	CONTINUE	
	ENDDO
      
      DEALLOCATE(RAL,MAL)	  
C	==================================================================
	

C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL2D1(COORD,PROPM,HSOC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),HSOC(NPT)
	DIMENSION NODEX(4)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

	NODEX(1:4) = [5,6,7,8]
C	=====================================================

C	==========================================
C	SETTING ELEMENT GAUSSIAN INTEGRATION POINT
C	STANDARD GAUSS INTEGRATION FOR 8 NODES SOLID ELEMENT
C	(MGR,MGS,MGT) = (2,2,2)
C	=======================
      MGR = 2
      MGS = 2
	MGT = 2

	IPT = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS  FOR HEAT SOURCE (CAN MOVE TO ELEMENT LOOP IF CHANGE WITH THE TIME)
C	=========================================
      DO 900  IGR=1,MGR
      RI = GLOC(IGR,MGR)
	DO 900  IGS=1,MGS
	SI = GLOC(IGS,MGS)
	IPT = IPT + 1
	WT = GWT(IGR,MGR)*GWT(IGS,MGS)
C     ====================================
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ===================================================
	CALL SHAP2D(RI,SI,H,P,NODEX,NNM)
	CALL JACO2H(COORD,P,XJ,XJI,DET,MEL,NNM)

	DVOL = WT*DET
C	========================================
C	COMPUTE HEAT SOURCE LOAD VECTOR
C	========================================
	DO I  = 1,NNM
	RH(I) = RH(I) + HSOC(IPT)*H(I)*DVOL
	ENDDO

 900  CONTINUE

C	==========================================


C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL2D2(COORD,PROPM,PROPG,QVC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT



	PARAMETER (NFACE=6)

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),PROPG(1),QVC(NFACE)
	DIMENSION NODEX(4)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

C	THICKNESS
	TH = PROPG(2)

	RH = 0.0D0
	REF= 0.0D0

	NODEX(1:4) = [5,6,7,8]
C	=====================================================

	MG1   = 2
	MGR   = 2
	MGS   = 2

	DO 1200 IFACE = 1,NFACE

	QVCF = QVC(IFACE)

	SELECT CASE(IFACE)
	CASE(1,2,3,4)
C	-------------------------------------------------------
	DO 1000 II = 1,MG1
C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE2D(IFACE,II,RI,SI,WT,MG1)
C	SHAPE FUNCTION
	CALL SHAP2D(RI,SI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO2H(COORD,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL MEMBVECTOR(IFACE,XJ,DET)
	DFACE = WT*TH*DET
	DO I  = 1,NNM
	RH(I) = RH(I) + QVCF*H(I)*DFACE
	ENDDO
1000	CONTINUE
C	-------------------------------------------------------
	CASE(5,6)
C	-------------------------------------------------------
	DO 1100  IGR=1,MGR
      RI = GLOC(IGR,MGR)
	DO 1100  IGS=1,MGS
	SI = GLOC(IGS,MGS)
	WT = GWT(IGR,MGR)*GWT(IGS,MGS)
C	SHAPE FUNCTION
	CALL SHAP2D(RI,SI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO2H(COORD,P,XJ,XJI,DET,MEL,NNM)
	DFACE = WT*DET
	DO I  = 1,NNM
	RH(I) = RH(I) + QVCF*H(I)*DFACE
	ENDDO
1100	CONTINUE
C	-------------------------------------------------------
	END SELECT


1200	CONTINUE

C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL2D3(COORD,PROPM,PROPG,TCD,TMP,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	PARAMETER (NFACE=6)

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),PROPG(1),TCD(NFACE)
	DIMENSION NODEX(4)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

	NODEX(1:4) = [5,6,7,8]
C	=====================================================
C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

C	THICKNESS
	TH = PROPG(2)

	MG1   = 2
	MGR   = 2
	MGS   = 2

	DO 1200 IFACE = 1,NFACE
	FAC = TCD(IFACE)
	TCD(IFACE) = 1.0D0
	IF(FAC.EQ.0.0) GOTO 1200

	SELECT CASE(IFACE)
	CASE(1,2,3,4)
C	-------------------------------------------------------
	DO 1000 II = 1,MG1
C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GFACE2D(IFACE,II,RI,SI,WT,MG1)
C	SHAPE FUNCTION
	CALL SHAP2D(RI,SI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO2H(COORD,P,XJ,XJI,DET,MEL,NNM)
C	DETERMINE THE AREA FACTOR ACCORDING TO FACE
	CALL MEMBVECTOR(IFACE,XJ,DET)
	DFACE = WT*TH*DET
	DO I  = 1,NNM
	RH(I) = RH(I) + FAC*HCV*TMP*H(I)*DFACE
	ENDDO
1000	CONTINUE
C	-------------------------------------------------------
	CASE(5,6)
C	-------------------------------------------------------
	DO 1100  IGR=1,MGR
      RI = GLOC(IGR,MGR)
	DO 1100  IGS=1,MGS
	SI = GLOC(IGS,MGS)
	WT = GWT(IGR,MGR)*GWT(IGS,MGS)
C	SHAPE FUNCTION
	CALL SHAP2D(RI,SI,H,P,NODEX,NNM)
C	JACOBIAN
	CALL JACO2H(COORD,P,XJ,XJI,DET,MEL,NNM)
	DFACE = WT*DET
	DO I  = 1,NNM
	RH(I) = RH(I) + FAC*HCV*TMP*H(I)*DFACE
	ENDDO
1100	CONTINUE
C	-------------------------------------------------------
	END SELECT

1200	CONTINUE

C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATSTOR(H1,H2,NELE,N,IEL)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION H1(N,NELE),H2(N)

	H1(1:N,IEL) = H1(1:N,IEL) + H2(1:N)


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATAMBT(TAMBT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

	DIMENSION TAMBT(2,LCS)

	TAMBT = 0.0

	READ(ITI,*)
	READ(ITI,*) NAMBT

	DO I = 1,NAMBT
	READ(ITI,*) ILC,TMP
	TAMBT(1,ILC) = TMP
	TAMBT(2,ILC) = 1.0
	ENDDO


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================



C	=======================================================================
C	======================START OF 3 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE HEATLOD2DT(ID,XYZ,LM,PROPM,MTSET,PROPG,IGSET,
	1					 IGIDM,TAMBT,HSOC,HCON,TLAMT,NLL,IND) 
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
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
      
	PARAMETER(NFACE=5)
	DIMENSION HSOC(NELE*NPT,LCS),HCON(NELE*NFACE,LCS)
	DIMENSION TLAMT(NELE*NFACE,LCS)
	DIMENSION HSOW(NPT),HCOW(NFACE),TLAMW(NFACE),TAMBT(2,LCS)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(1)

	CALL CLEARA(HSOC,NELE*NPT*LCS)
	CALL CLEARA(HCON,NELE*NFACE*LCS)
	CALL CLEARI(TLAMT,NELE*NFACE*LCS)

	NHSOC = 0
	NHCON = 0
	NHAMT = 0

	SELECTCASE (IND)
	CASE(1) 
	NHSOC = NLL
	CASE(2) 
	NHCON = NLL
	CASE(3) 
	NHCON = NLL
	CASE(4) 
	NHAMT = NLL
	CASE(5) 
	NHAMT = NLL
	ENDSELECT
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHSOC
	READ(ITI,*) MEMBA,HSO,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	HSOW(1:NPT) = HSO
	CALL ELEREODER(IGIDMEM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL HEATSTOR(HSOC(1,ILCN),HSOW,NELE,NPT,MEMBA)
	CALL HEATL2D1T(XYZ(1,MEMBA),PROPM(1,MSET),HSOW,RL)
	CALL FIXMRE(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)	  
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHCON

	SELECTCASE (IND)
	CASE(2) 
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC
	HCOW = 0.0D0
	HCOW(IFACE) = DUM1
	CASE(3) 
	READ(ITI,*) MEMBA,DUM1,ILCN,ILCC
	HCOW = 0.0D0
	HCOW(5) = DUM1
	ENDSELECT

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	

	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	NSET = IGSET(MEMBA)
	CALL HEATSTOR(HCON(1,ILCN),HCOW,NELE,NFACE,MEMBA)
	CALL HEATL2D2T(XYZ(1,MEMBA),PROPM(1,MSET),PROPG(1,NSET),
	1			  HCOW,RL)
	CALL FIXMRE(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)	  
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHAMT

	SELECTCASE (IND)
	CASE(4) 
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC
	TLAMW = 0.0D0
	TLAMW(IFACE) = DUM1
	CASE(5) 
	READ(ITI,*) MEMBA,DUM1,ILCN,ILCC
	TLAMW = 0.0D0
	TLAMW(5) = DUM1
	ENDSELECT

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	NSET = IGSET(MEMBA)
	TMP  = TAMBT(1,ILCN)
	LFG  = TAMBT(2,ILCN)
	IF(LFG.EQ.0) GOTO 100
	CALL HEATSTOR(TLAMT(1,ILCN),TLAMW,NELE,NFACE,MEMBA)
	CALL HEATL2D3T(XYZ(1,MEMBA),PROPM(1,MSET),PROPG(1,NSET),
	1				  TLAMW,TMP,RL)
	CALL FIXMRE(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
100	CONTINUE	
	ENDDO
      
      DEALLOCATE(RAL,MAL)	
C	==================================================================
	

C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL2D1T(COORD,PROPM,HSOC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),HSOC(NPT)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0
C	==========================================
C	SETTING ELEMENT GAUSSIAN INTEGRATION POINT
C	STANDARD GAUSS INTEGRATION FOR 8 NODES SOLID ELEMENT
C	(MGR,MGS,MGT) = (2,2,2)
C	=======================
      MGR = 3

	IPT = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS  FOR HEAT SOURCE (CAN MOVE TO ELEMENT LOOP IF CHANGE WITH THE TIME)
C	=========================================
      DO 900  IGR=1,MGR

	CALL GAUSST(RI,SI,TI,WT,IGR,MGR,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H,P,NNM)
C	JACOBIAN
	CALL JACO2HT(COORD,P,XJ,XJI,DET,MEL,NNM)

	DVOL = WT*DET
C	========================================
C	COMPUTE HEAT SOURCE LOAD VECTOR
C	========================================
	DO I  = 1,NNM
	RH(I) = RH(I) + HSOC(IPT)*H(I)*DVOL
	ENDDO

 900  CONTINUE

C	==========================================

C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL2D2T(COORD,PROPM,PROPG,QVC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	PARAMETER (NFACE=5)

	DIMENSION COORD(NCO,NNM)
	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),PROPG(1),QVC(NFACE)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

C	THICKNESS
	TH = PROPG(2)

	RH = 0.0D0
	REF= 0.0D0

C	=====================================================
	MG1   = 2
	MGR   = 3

	DO 1200 IFACE = 1,NFACE

	QVCF = QVC(IFACE)

	SELECT CASE(IFACE)
	CASE(1,2,3)
C	-------------------------------------------------------
	DO 1000 II = 1,MG1

C	DEFINE GAUSS LOACATION AND WEIGTH
	RI = GLOC(II,MG1)
	WT =  GWT(II,MG1)
C	SHAPE FUNCTION
	CALL SHAP1D(RI,H,P,2)
C	JACOBIAN
	CALL JACO1H(IFACE,COORD,P,DET,MEL,2)
C	REARRANGE POSITION OF SHAPE FUNCTION
	CALL RESHPT(H,IFACE,NNM)

	DFACE = WT*TH*DET
	DO I  = 1,NNM
	RH(I) = RH(I) + QVCF*H(I)*DFACE
	ENDDO
1000	CONTINUE
C	-------------------------------------------------------
	CASE(4,5)
C	-------------------------------------------------------
	DO 1100  IGR=1,MGR

	CALL GAUSST(RI,SI,TI,WT,IGR,MGR,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H,P,NNM)
C	JACOBIAN
	CALL JACO2HT(COORD,P,XJ,XJI,DET,MEL,NNM)

	DFACE = WT*DET
	DO I  = 1,NNM
	RH(I) = RH(I) + QVCF*H(I)*DFACE
	ENDDO
1100	CONTINUE
C	-------------------------------------------------------
	END SELECT


1200	CONTINUE

	
C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL2D3T(COORD,PROPM,PROPG,TCD,TMP,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	PARAMETER (NFACE=5)

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),PROPG(1),TCD(NFACE)
	DIMENSION H(NNM),P(NCO,NNM),XJ(NCO,NCO),XJI(NCO,NCO)

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

C	=====================================================
C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

C	THICKNESS
	TH = PROPG(2)

	MG1   = 2
	MGR   = 3

	DO 1200 IFACE = 1,NFACE
	FAC = TCD(IFACE)
	TCD(IFACE) = 1.0D0
	IF(FAC.EQ.0.0) GOTO 1200

	SELECT CASE(IFACE)
	CASE(1,2,3)
C	-------------------------------------------------------
	DO 1000 II = 1,MG1

C	DEFINE GAUSS LOACATION AND WEIGTH
	RI = GLOC(II,MG1)
	WT =  GWT(II,MG1)
C	SHAPE FUNCTION
	CALL SHAP1D(RI,H,P,2)
C	JACOBIAN
	CALL JACO1H(IFACE,COORD,P,DET,MEL,2)
C	REARRANGE POSITION OF SHAPE FUNCTION
	CALL RESHPT(H,IFACE,NNM)

	DFACE = WT*TH*DET

	DO I  = 1,NNM
	RH(I) = RH(I) + FAC*HCV*TMP*H(I)*DFACE
	ENDDO
1000	CONTINUE
C	-------------------------------------------------------
	CASE(4,5)
C	-------------------------------------------------------
	DO 1100  IGR=1,MGR

	CALL GAUSST(RI,SI,TI,WT,IGR,MGR,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H,P,NNM)
C	JACOBIAN
	CALL JACO2HT(COORD,P,XJ,XJI,DET,MEL,NNM)

	DFACE = WT*DET
	DO I  = 1,NNM
	RH(I) = RH(I) + FAC*HCV*TMP*H(I)*DFACE
	ENDDO
1100	CONTINUE
C	-------------------------------------------------------
	END SELECT

1200	CONTINUE

C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================




C	=======================================================================
C	===================START OF TETRAHEDRA HEAT ELEMENT====================
C	=======================================================================
	SUBROUTINE HEATLOD3DT(ID,XYZ,LM,PROPM,MTSET,IGIDM,TAMBT,
	1					  HSOC,HCON,TLAMT,NLL,IND) 
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
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)

	PARAMETER(NFACE=4)
	DIMENSION HSOC(NELE*NPT,LCS),HCON(NELE*NFACE,LCS)
	DIMENSION TLAMT(NELE*NFACE,LCS)
	DIMENSION HSOW(NPT),HCOW(NFACE),TLAMW(NFACE),TAMBT(2,LCS)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE),PROPM(NMP,1),MTSET(1)

	CALL CLEARA(HSOC,NELE*NPT*LCS)
	CALL CLEARA(HCON,NELE*NFACE*LCS)
	CALL CLEARI(TLAMT,NELE*NFACE*LCS)


	NHSOC = 0
	NHCON = 0
	NHAMT = 0
	SELECTCASE (IND)
	CASE(1) 
	NHSOC = NLL
	CASE(2) 
	NHCON = NLL
	CASE(3) 
	NHAMT = NLL
	ENDSELECT
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHSOC
	READ(ITI,*) MEMBA,HSO,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	HSOW(1:NPT) = HSO
	CALL ELEREODER(IGIDMEM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL HEATSTOR(HSOC(1,ILCN),HSOW,NELE,NPT,MEMBA)
	CALL HEATL3D1T(XYZ(1,MEMBA),PROPM(1,MSET),HSOW,RL)
	CALL FIXSOL(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)	 
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHCON
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	HCOW = 0.0D0
	HCOW(IFACE) = DUM1
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL HEATSTOR(HCON(1,ILCN),HCOW,NELE,NFACE,MEMBA)
	CALL HEATL3D2T(XYZ(1,MEMBA),PROPM(1,MSET),HCOW,RL)
	CALL FIXSOL(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	  
      ENDDO
      
      DEALLOCATE(RAL,MAL)	 
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NHAMT
	READ(ITI,*) MEMBA,IFACE,DUM1,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	TLAMW = 0.0D0
	TLAMW(IFACE) = DUM1
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	TMP  = TAMBT(1,ILCN)
	LFG  = TAMBT(2,ILCN)
	IF(LFG.EQ.0) GOTO 100
	CALL HEATSTOR(TLAMT(1,ILCN),TLAMW,NELE,NFACE,MEMBA)
	CALL HEATL3D3T(XYZ(1,MEMBA),PROPM(1,MSET),TLAMW,TMP,RL)
	CALL FIXSOL(RL,KEG,MEMBA,NEF,NNM)
C	MECHANICAL DOF IS NOT YET INCLUDE SO NO NEED TO TRANSFORM FOR SKEW SUPPORT
	IEFL = 0
	DO J = 1,NEF
	IEQ = LM(J,MEMBA)
      IF (IEQ.NE.0) THEN
          IEFL = IEFL + 1
          RAL(IEFL) = RL(J) 
          MAL(IEFL) = IEQ
      ENDIF
	ENDDO
	CALL LDASEM_NEW (RAL,MAL,IEFL)	
100	CONTINUE	
	ENDDO
      
      DEALLOCATE(RAL,MAL)	
C	==================================================================
	
C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL3D1T(COORD,PROPM,HSOC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),HSOC(NPT)
	DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

C	=====================================================
      MGR = 4

	IPT = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS  FOR HEAT SOURCE (CAN MOVE TO ELEMENT LOOP IF CHANGE WITH THE TIME)
C	=========================================
      DO 900  IGR=1,MGR
	IPT = IPT + 1

	CALL GAUSST (RI,SI,TI,WT,IPT,MGR,1)
C     ====================================
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ===================================================
	CALL SHAP3DT(RI,SI,TI,H,P,NNM)
	CALL JACO3D(COORD,P,XJ,XJI,DET,MEL,NNM)
C	MODIFY DET DUE TO TETRAHEDRA
	DET = DET/6.0

	DVOL = WT*DET
C	========================================
C	COMPUTE HEAT SOURCE LOAD VECTOR
C	========================================

	DO I  = 1,NNM
	RH(I) = RH(I) + HSOC(IPT)*H(I)*DVOL
	ENDDO

 900  CONTINUE

C	==========================================


C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL3D2T(COORD,PROPM,QVC,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT



	PARAMETER (NFACE=4)

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),QVC(NFACE)
	DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

C	=====================================================

	MG1   = 3

	NN(4,1:6) = [2,3,4,8,9,10]
	NN(1,1:6) = [1,3,4,6,9,7]
	NN(2,1:6) = [1,2,4,5,10,7]
	NN(3,1:6) = [1,2,3,5,8,6]

	DO 1000 IFACE = 1,NFACE

	QVCF = QVC(IFACE)

	DO 1000 II = 1,MG1

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GAUSST(RI,SI,TI,WT,II,MG1,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H2,G,3)
C	JACOBIAN
	DO I = 1,3
	COVR(I) = 0.0D0
	COVS(I) = 0.0D0
	DO J = 1,3
	MM = NN(IFACE,J)
	COVR(I) = COVR(I) + G(1,J)*COORD(I,MM)
	COVS(I) = COVS(I) + G(2,J)*COORD(I,MM)
	ENDDO
	ENDDO
	CALL VECPRD(COVR,COVS,VT)
	CALL SCALEN(VT,VT,DET,3)
C	TRIANGULAR AREA
	DET = 0.5*DET
C	REARRANGE SHAPE FUNCTION
	H(1:NNM) = 0.0D0
	DO I = 1,3
	MM = NN(IFACE,I)
	H(MM) = H2(I)
	ENDDO

	DFACE = WT*DET
	
	DO I  = 1,NNM
	RH(I) = RH(I) + QVCF*H(I)*DFACE
	ENDDO

1000	CONTINUE

C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATL3D3T(COORD,PROPM,TCD,TMP,REF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT


	PARAMETER (NFACE=4)

	DIMENSION COORD(NCO,NNM)

	DIMENSION RH(NNM),REF(NEF)
	DIMENSION PROPM(1),TCD(NFACE)
	DIMENSION NODEX(13)
	DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	RH = 0.0D0
	REF= 0.0D0

C	=====================================================
C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

	MG1   = 3

	NN(4,1:6) = [2,3,4,8,9,10]
	NN(1,1:6) = [1,3,4,6,9,7]
	NN(2,1:6) = [1,2,4,5,10,7]
	NN(3,1:6) = [1,2,3,5,8,6]


	DO 1000 IFACE = 1,NFACE
	FAC = TCD(IFACE)
	TCD(IFACE) = 1.0D0

	DO 1000 II = 1,MG1

C	DEFINE GAUSS LOACATION AND WEIGTH
	CALL GAUSST(RI,SI,TI,WT,II,MG1,0)
C	SHAPE FUNCTION
	CALL SHAP2DT(RI,SI,H2,G,3)
C	JACOBIAN
	DO I = 1,3
	COVR(I) = 0.0D0
	COVS(I) = 0.0D0
	DO J = 1,3
	MM = NN(IFACE,J)
	COVR(I) = COVR(I) + G(1,J)*COORD(I,MM)
	COVS(I) = COVS(I) + G(2,J)*COORD(I,MM)
	ENDDO
	ENDDO
	CALL VECPRD(COVR,COVS,VT)
	CALL SCALEN(VT,VT,DET,3)
C	TRIANGULAR AREA
	DET = 0.5*DET
C	REARRANGE SHAPE FUNCTION
	H(1:NNM) = 0.0D0
	DO I = 1,3
	MM = NN(IFACE,I)
	H(MM) = H2(I)
	ENDDO

	DFACE = WT*DET

	DO I  = 1,NNM
	RH(I) = RH(I) + FAC*HCV*TMP*H(I)*DFACE
	ENDDO


1000	CONTINUE

	
C	==============================
C	ASSEMBLE THE TEMPERATURE
C	==============================
	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE RHSPRN(ID,FOV,IDOF,ISO,NSN,NSF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION ID(NSF,NSN),FOV(1),IDOF(9),D(9)

	M = 0
	DO I = 1,NSF
	IF(IDOF(I).GT.M) M = IDOF(I)
	ENDDO

	DO I = 1,NSN
	D = 0.0
	DO J = 1,NSF
	IEQ = ID(J,I)
	NN = IDOF(J)
	IF(IEQ.NE.0.0) D(NN) = FOV(IEQ)
	ENDDO
	WRITE(ISO,100) I,D(1:M)
	ENDDO


100	FORMAT(I5,X,9E15.6)

	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
