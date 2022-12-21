C	=======================================================================
C	======================START OF 8 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE SEEPLOD3D(ID,XYZ,LM,PROPM,MTSET,IGIDM,NLL,IND) 
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
      ALLOCATABLE RAL(:) !SONGSAK TWEAK SPEED OCT2019
      ALLOCATABLE MAL(:)
	  
	PARAMETER(NFACE=6)
	DIMENSION FACEF(NFACE)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE),PROPM(NMP,1),MTSET(1)

	NFLOW = 0
	SELECTCASE (IND)
	CASE(1) 
	NFLOW = NLL
	ENDSELECT
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NFLOW
	READ(ITI,*) MEMBA,IFACE,FLOW,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	FACEF = 0.0D0
	FACEF(IFACE) = FLOW
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL SEEPL3D1(XYZ(1,MEMBA),PROPM(1,MSET),FACEF,RL)
	CALL   FIXSOL(RL,KEG,MEMBA,NEF,NNM)
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


C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SEEPL3D1(COORD,PROPM,QVC,REF)
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
	NMTEM = 1                !NUMBER OF FLOW DOF


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
	DO I = 1,NNM             !FOR FLOW DOF
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
	SUBROUTINE SEEPLOD2D(ID,XYZ,LM,PROPM,MTSET,PROPG,IGSET,
	1					 IGIDM,NLL,IND) 
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
	DIMENSION FACEF(NFACE)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(1)


	NFLOW = 0
	SELECTCASE (IND)
	CASE(1) 
	NFLOW = NLL
	ENDSELECT

C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NFLOW

	READ(ITI,*) MEMBA,IFACE,FLOW,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	FACEF = 0.0D0
	FACEF(IFACE) = FLOW
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	NSET = IGSET(MEMBA)
	CALL SEEPL2D1(XYZ(1,MEMBA),PROPM(1,MSET),PROPG(1,NSET),
	1			  FACEF,RL)
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

C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SEEPL2D1(COORD,PROPM,PROPG,QVC,REF)
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
	NMTEM = 1                !NUMBER OF FLOW DOF

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
	DO I = 1,NNM             !FOR FLOW DOF
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
C	======================START OF 3 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE SEEPLOD2DT(ID,XYZ,LM,PROPM,MTSET,PROPG,IGSET,
	1					  IGIDM,NLL,IND) 
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

	DIMENSION FACEF(NFACE)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE)
	DIMENSION PROPM(NMP,1),MTSET(1),PROPG(NGP,1),IGSET(1)


	NFLOW = 0
	SELECTCASE (IND)
	CASE(1) 
	NFLOW = NLL
	ENDSELECT
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NFLOW

	READ(ITI,*) MEMBA,IFACE,FLOW,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	FACEF = 0.0D0
	FACEF(IFACE) = FLOW

	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	NSET = IGSET(MEMBA)
	CALL SEEPL2D1T(XYZ(1,MEMBA),PROPM(1,MSET),PROPG(1,NSET),
	1			   FACEF,RL)
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

	

C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SEEPL2D1T(COORD,PROPM,PROPG,QVC,REF)
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
	NMTEM = 1                !NUMBER OF FLOW DOF

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
	DO I = 1,NNM             !FOR FLOW DOF
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
	SUBROUTINE SEEPLOD3DT(ID,XYZ,LM,PROPM,MTSET,IGIDM,NLL,IND) 
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

	DIMENSION FACEF(NFACE)

	DIMENSION IGIDM(1),XYZ(NCO*NNM,NELE),RL(NEF),RR(NEQ)
      DIMENSION ID(NSF,NSN),LM(NEF,NELE),PROPM(NMP,1),MTSET(1)


	NFLOW = 0
	SELECTCASE (IND)
	CASE(1) 
	NFLOW = NLL
	ENDSELECT
C	==================================================================
      ALLOCATE(RAL(NEF),MAL(NEF))
      
	DO I = 1,NFLOW
	READ(ITI,*) MEMBA,IFACE,FLOW,ILCN,ILCC

      RAL(1:NEF) = 0.0D0
      MAL(1:NEF) = 0	
      
	FACEF = 0.0D0
	FACEF(IFACE) = FLOW
	CALL ELEREODER(IGIDM,NELE,MEMBA)
	MSET = MTSET(MEMBA)
	CALL SEEPL3D1T(XYZ(1,MEMBA),PROPM(1,MSET),FACEF,RL)
	CALL    FIXSOL(RL,KEG,MEMBA,NEF,NNM)
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

	
C
      RETURN
C
      END
C
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SEEPL3D1T(COORD,PROPM,QVC,REF)
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
	NMTEM = 1                !NUMBER OF FLOW DOF

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
	DO I = 1,NNM             !FOR FLOW DOF
	II = (NMCHA+NMTEM)*I
	REF(II) = REF(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
