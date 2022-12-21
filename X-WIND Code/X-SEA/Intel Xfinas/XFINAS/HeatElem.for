C	SOLID HEAT TRANSFER LIBRARY
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATELE3D(COORD,PROPM,EDIS,S,RE,TCD,TMP,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION COORD(1),PROPM(1),EDIS(1),RE(1),TCD(1),TMP(1) 
	DIMENSION WA(1)

	SELECT CASE(ISTYP)

	CASE(3)
	CALL HEATELE3(COORD,PROPM,EDIS,S,RE,TCD,TMP,WA)

	CASE(4)
	CALL HEATELE3T(COORD,PROPM,EDIS,S,RE,TCD,TMP,WA)

	END SELECT


	RETURN
	END


C	=======================================================================
C	======================START OF 8 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE HEATELE3(COORD,PROPM,EDIS,S,RE,TCD,TMP,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

C	PARAMETER FOR TRANSIENT HEAT tRANSFER SONGSAK MAR2007
	COMMON /SHEAT1/ PIMPLT,KSTEAD

	PARAMETER (NFACE=6)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM),TMP(2)
	DIMENSION HM(3)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1),TCD(NFACE)
	DIMENSION BH(3,NNM),NODEX(13),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)


      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	FAC1 = PIMPLT
	FAC2 = 1.0 - FAC1
C	TRANSFER NODAL VALUE INTO TEMPERATURE PART
	DO I = 1,NNM
	II = (NMCHA+NMTEM)*I
	ETEMP(I) = EDIS(II)
	ENDDO

	S  = 0.0D0
	RE = 0.0D0

	SH = 0.0D0
	RH = 0.0D0

	
	NODEX(1:13) = [9,10,11,12,13,14,15,16,17,18,19,20,21]
C	==========================================
C	SETTING ELEMENT GAUSSIAN INTEGRATION POINT
C	STANDARD GAUSS INTEGRATION FOR 8 NODES SOLID ELEMENT
C	(MGR,MGS,MGT) = (2,2,2)
C	=======================

      MGR = 2
      MGS = 2
	MGT = 2
C	===============================================
C	CALL HEAT CONDUCTION MATRIX (CONSTITUTIVE HEAT)
C	===============================================
	DO I = 1,3
	HM(I)  = PROPM(5+I)
	ENDDO

C	HEAT CAPACITY COEFICIENT
	HCP = PROPM(11)

C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

	IPT = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS
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
C	COMPUTE HEAT CONDUCTION STIFFNESS MATRIX
C	========================================
	BH = MATMUL(XJI,P)

	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,3
	SH(I,J) = SH(I,J) + FAC1*HM(K)*BH(K,I)*BH(K,J)*DVOL

	IF(KSTEAD.EQ.1) THEN
	RH(I)   = RH(I)   + FAC2*HM(K)*BH(K,I)*BH(K,J)*DVOL*ETEMP(J)
	ENDIF

	ENDDO
	ENDDO
	ENDDO

	DO I = 1,NNM
	DO J = 1,NNM

	IF(KSTEAD.EQ.1) THEN
	SH(I,J) = SH(I,J) + HCP*H(I)*H(J)*DVOL/DTINC
	RH(I)   = RH(I)   - HCP*H(I)*H(J)*DVOL*ETEMP(J)/DTINC
	ENDIF

	ENDDO
	ENDDO

 900  CONTINUE
 	
C	WRITE(*,*) 'KAE',ETEMP
C	==========================================
	LFG = TMP(2)
	IF(LFG.EQ.0) GOTO 1200

	MG1   = 2
	MG2   = 2

	DO 1100 IFACE = 1,NFACE
	
	IF(TCD(IFACE).EQ.0.0) GOTO 1100

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
	
	DO I = 1,NNM
	DO J = 1,NNM
	SH(I,J)  = SH(I,J) + FAC1*HCV*H(I)*H(J)*DFACE*TCD(IFACE)
	IF(KSTEAD.EQ.1) THEN
	RH(I)    = RH(I)   + FAC2*HCV*H(I)*H(J)*DFACE*ETEMP(J)*TCD(IFACE)
	ENDIF

	ENDDO
	ENDDO

1000	CONTINUE

1100	CONTINUE


C	==============================
C	TRANSFER LOCAL HEAT CONDUCTION 
C	TO GLOBAL ELEMENT STIFFNESS
C	==============================
1200	CONTINUE

	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	RE(II) = RE(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE GFACE3D(IFACE,II,JJ,RI,SI,TI,WT,MG1,MG2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------------------------------------------------
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT

C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R= 1.0
C	FACE 2; R=-1.0
C	FACE 3: S= 1.0
C	FACE 4: S=-1.0
C	FACE 5: T= 1.0
C	FACE 6: T=-1.0
C	---------------------

	SELECT CASE (IFACE)
	
	CASE(2)  !5
	RI = 1.0D0
	SI = GLOC(II,MG1)
	TI = GLOC(JJ,MG2)
	WT = GWT(II,MG1)*GWT(JJ,MG2)

	CASE(4)  !3
	RI =-1.0D0
	SI = GLOC(II,MG1)
	TI = GLOC(JJ,MG2)
	WT = GWT(II,MG1)*GWT(JJ,MG2)

	CASE(3)  !2
	SI = 1.0D0
	RI = GLOC(II,MG1)
	TI = GLOC(JJ,MG2)
	WT = GWT(II,MG1)*GWT(JJ,MG2)

	CASE(5)  !4
	SI =-1.0D0
	RI = GLOC(II,MG1)
	TI = GLOC(JJ,MG2)
	WT = GWT(II,MG1)*GWT(JJ,MG2)

	CASE(6)  !6
	TI = 1.0D0
	RI = GLOC(II,MG1)
	SI = GLOC(JJ,MG2)
	WT = GWT(II,MG1)*GWT(JJ,MG2)

	CASE(1)  !1
	TI =-1.0D0
	RI = GLOC(II,MG1)
	SI = GLOC(JJ,MG2)
	WT = GWT(II,MG1)*GWT(JJ,MG2)

	END SELECT
	

	RETURN

	END


C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE SOLIDVECTOR(NFACE,XJ,DET)
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


	CALL VECPRD (VA,VB,VC)
	CALL SCALEN (VC,VC,DET,3)

	RETURN

	END

C	=====================================================================
C	=====================================================================	
C	=====================================================================

C	MEMBRANE HEAT TRANSFER LIBRARY
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE HEATELE2D(COORD,PROPM,PROPG,EDIS,S,RE,TCD,TMP,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION COORD(1),PROPM(1),PROPG(1),EDIS(1),RE(1),TCD(1),TMP(1) 
	DIMENSION WA(1)

	SELECT CASE(ISTYP)

	CASE(3)
	CALL HEATELE2(COORD,PROPM,PROPG,EDIS,S,RE,TCD,TMP,WA)

	CASE(4)
	CALL HEATELE2T(COORD,PROPM,PROPG,EDIS,S,RE,TCD,TMP,WA)

	END SELECT


	RETURN
	END


C	=======================================================================
C	======================START OF 4 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE HEATELE2(COORD,PROPM,PROPG,EDIS,S,RE,TCD,TMP,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

C	PARAMETER FOR TRANSIENT HEAT tRANSFER SONGSAK MAR2007
	COMMON /SHEAT1/ PIMPLT,KSTEAD

	PARAMETER (NFACE=6)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM),TMP(2)
	DIMENSION HM(2)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1),PROPG(1),TCD(NFACE)
	DIMENSION BH(2,NNM),NODEX(4),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)
	DIMENSION VXY(2),TMAT(2,2) !FOR TRANSFORMATION OF FLOW VELOCITY OF SEEPAGE

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	FAC1 = PIMPLT
	FAC2 = 1.0 - FAC1
C	TRANSFER NODAL VALUE INTO TEMPERATURE PART
	DO I = 1,NNM
	II = (NMCHA+NMTEM)*I
	ETEMP(I) = EDIS(II)
	ENDDO

	S  = 0.0D0
	RE = 0.0D0

	SH = 0.0D0
	RH = 0.0D0

	
	NODEX(1:4) = [5,6,7,8]
C	==========================================
C	SETTING ELEMENT GAUSSIAN INTEGRATION POINT
C	STANDARD GAUSS INTEGRATION FOR 8 NODES SOLID ELEMENT
C	(MGR,MGS,MGT) = (2,2,2)
C	=======================

      MGR = 2
      MGS = 2
	MGT = 2
	IF(NNM.EQ.8) THEN
	MGR = 3
      MGS = 3
	MGT = 3
	ENDIF
C	===============================================
C	CALL HEAT CONDUCTION MATRIX (CONSTITUTIVE HEAT)
C	===============================================
	DO I = 1,2
	HM(I)  = PROPM(5+I)
	ENDDO

C	HEAT CAPACITY COEFICIENT
	HCP = PROPM(11)

C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

C	THICKNESS
	TH  = PROPG(2)
	
	IPT = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS
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
	DVOL = WT*TH*DET

C	========================================
C	COMPUTE HEAT CONDUCTION STIFFNESS MATRIX
C	========================================
	BH = MATMUL(XJI,P)

	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,2
	SH(I,J) = SH(I,J) + FAC1*HM(K)*BH(K,I)*BH(K,J)*DVOL

	IF(KSTEAD.EQ.1) THEN
	RH(I)   = RH(I)   + FAC2*HM(K)*BH(K,I)*BH(K,J)*DVOL*ETEMP(J)
	ENDIF

	ENDDO
	ENDDO
	ENDDO

	DO I = 1,NNM
	DO J = 1,NNM

	IF(KSTEAD.EQ.1) THEN
	SH(I,J) = SH(I,J) + HCP*H(I)*H(J)*DVOL/DTINC
	RH(I)   = RH(I)   - HCP*H(I)*H(J)*DVOL*ETEMP(J)/DTINC
	ENDIF

	ENDDO
	ENDDO

 900  CONTINUE
 	
C	==========================================
	LFG = TMP(2)
	IF(LFG.EQ.0) GOTO 1300

	MG1   = 2
	DO 1200 IFACE = 1,NFACE
	IF(TCD(IFACE).EQ.0.0) GOTO 1200

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
	DO I = 1,NNM
	DO J = 1,NNM
	SH(I,J)  = SH(I,J) + FAC1*HCV*H(I)*H(J)*DFACE*TCD(IFACE)
	IF(KSTEAD.EQ.1) THEN
	RH(I)    = RH(I)   + FAC2*HCV*H(I)*H(J)*DFACE*ETEMP(J)*TCD(IFACE)
	ENDIF
	ENDDO
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
	DO I = 1,NNM
	DO J = 1,NNM
	SH(I,J)  = SH(I,J) + FAC1*HCV*H(I)*H(J)*DFACE*TCD(IFACE)
	IF(KSTEAD.EQ.1) THEN
	RH(I)    = RH(I)   + FAC2*HCV*H(I)*H(J)*DFACE*ETEMP(J)*TCD(IFACE)
	ENDIF
	ENDDO
	ENDDO
1100	CONTINUE
C	-------------------------------------------------------
	END SELECT

1200	CONTINUE
C	==============================
C	TRANSFER LOCAL HEAT CONDUCTION 
C	TO GLOBAL ELEMENT STIFFNESS
C	==============================
1300	CONTINUE

	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	RE(II) = RE(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE GFACE2D(IFACE,II,RI,SI,WT,MG1)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------------------------------------------------
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT

C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R= 1.0
C	FACE 2; R=-1.0
C	FACE 3: S= 1.0
C	FACE 4: S=-1.0
C	---------------------

	SELECT CASE (IFACE)
	
	CASE(4)
	RI = 1.0D0
	SI = GLOC(II,MG1)
	WT = GWT(II,MG1)

	CASE(2)
	RI =-1.0D0
	SI = GLOC(II,MG1)
	WT = GWT(II,MG1)

	CASE(1)
	SI = 1.0D0
	RI = GLOC(II,MG1)
	WT = GWT(II,MG1)

	CASE(3)
	SI =-1.0D0
	RI = GLOC(II,MG1)
	WT = GWT(II,MG1)

	END SELECT
	

	RETURN

	END


C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE MEMBVECTOR(NFACE,XJ,DET)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION XJ(2,2),VA(2)

	VA(1:2) = 0.0D0

	SELECT CASE(NFACE)

	CASE(4,2)

	VA(1) = XJ(2,1)
	VA(2) = XJ(2,2)

	CASE(1,3)

	VA(1) = XJ(1,1)
	VA(2) = XJ(1,2)

	END SELECT

	DET = SQRT(VA(1)*VA(1) + VA(2)*VA(2))

	RETURN

	END

C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE JACO2H1(XY,P,XJ,XJI,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(2,1),P(2,NNO),XJ(2,2),XJI(4)
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  I=1,2
      DO 100  J=1,2
      DUM = 0.0
      DO 90   K=1,NNO
 90   DUM = DUM + P(I,K)*XY(J,K)
 100  XJ(I,J) = DUM
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = XJ(1,1)*XJ(2,2) - XJ(2,1)*XJ(1,2)
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')
C     -----------------------------
C     INVERSE (XJI) OF THE JACOBIAN
C     -----------------------------
      DUM = 1.0/DET
      XJI(1) =  XJ(2,2)*DUM
      XJI(2) = -XJ(2,1)*DUM
      XJI(3) = -XJ(1,2)*DUM
      XJI(4) =  XJ(1,1)*DUM
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE JACO2H(XY,P,XJ,XJI,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------
      DIMENSION XY(3,NNO),P(2,NNO),VR(3),VS(3),VT(3),XJI(2,2)
      DIMENSION COVR(3),COVS(3),RV(3),SV(3),XJ(2,2)
C
      CALL CLEARA (COVR,3)
      CALL CLEARA (COVS,3)
      DO 20 I=1,NNO
      DO 20 J=1,3
      COVR(J)=COVR(J)+P(1,I)*XY(J,I)
   20 COVS(J)=COVS(J)+P(2,I)*XY(J,I)

      CALL VECPRD (COVR,COVS,VT)
      CALL SCALEN (VT,VT,DET,3)

      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')

      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)
      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)


      CALL SCAPRD (COVR,VR,F1,3)
      CALL SCAPRD (COVS,VR,F2,3)
      CALL SCAPRD (COVR,VS,F3,3)
      CALL SCAPRD (COVS,VS,F4,3)

      XJ(1,1)=F1
	XJ(2,1)=F2
	XJ(1,2)=F3
	XJ(2,2)=F4
	XJI(1,1)= F4/DET
      XJI(2,1)=-F2/DET
      XJI(1,2)=-F3/DET
      XJI(2,2)= F1/DET

      RETURN

      END

C	=====================================================================
C	=====================================================================	
C	=====================================================================




C	=====================================================================
C	=====================START OF 3 NODE HEAT ELEMENT====================	
C	=====================================================================
	SUBROUTINE HEATELE2T(COORD,PROPM,PROPG,EDIS,S,RE,TCD,TMP,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

C	PARAMETER FOR TRANSIENT HEAT tRANSFER SONGSAK MAR2007
	COMMON /SHEAT1/ PIMPLT,KSTEAD

	PARAMETER (NFACE=5)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM),TMP(2)
	DIMENSION HM(2)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1),PROPG(1),TCD(NFACE)
	DIMENSION BH(2,NNM),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)
	DIMENSION VXY(2),TMAT(2,2) !FOR TRANSFORMATION OF FLOW VELOCITY OF SEEPAGE

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	FAC1 = PIMPLT
	FAC2 = 1.0 - FAC1
C	TRANSFER NODAL VALUE INTO TEMPERATURE PART
	DO I = 1,NNM
	II = (NMCHA+NMTEM)*I
	ETEMP(I) = EDIS(II)
	ENDDO

	S  = 0.0D0
	RE = 0.0D0

	SH = 0.0D0
	RH = 0.0D0

C	==========================================
C	SETTING ELEMENT GAUSSIAN INTEGRATION POINT
C	STANDARD GAUSS INTEGRATION FOR 8 NODES SOLID ELEMENT
C	(MGR,MGS,MGT) = (2,2,2)
C	=======================

      MGR = 3
C	===============================================
C	CALL HEAT CONDUCTION MATRIX (CONSTITUTIVE HEAT)
C	===============================================
	DO I = 1,2
	HM(I)  = PROPM(5+I)
	ENDDO

C	HEAT CAPACITY COEFICIENT
	HCP = PROPM(11)

C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

C	THICKNESS
	TH  = PROPG(2)
	
	AREA = 0.0D0
	IPT  = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS
C	=========================================
      DO 900  IGR=1,MGR
	IPT = IPT + 1

      CALL GAUSST(RI,SI,TI,WT,IGR,MGR,0)
C     ====================================
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     =================================================== 
	CALL SHAP2DT(RI,SI,H,P,NNM) 
	CALL JACO2HT(COORD,P,XJ,XJI,DET,MEL,NNM)
	DVOL = WT*TH*DET

C	========================================
C	COMPUTE HEAT CONDUCTION STIFFNESS MATRIX
C	========================================
	BH = MATMUL(XJI,P)


	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,2
	SH(I,J) = SH(I,J) + FAC1*HM(K)*BH(K,I)*BH(K,J)*DVOL

	IF(KSTEAD.EQ.1) THEN
	RH(I)   = RH(I)   + FAC2*HM(K)*BH(K,I)*BH(K,J)*DVOL*ETEMP(J)
	ENDIF

	ENDDO
	ENDDO
	ENDDO

	DO I = 1,NNM
	DO J = 1,NNM

	IF(KSTEAD.EQ.1) THEN
	SH(I,J) = SH(I,J) + HCP*H(I)*H(J)*DVOL/DTINC
	RH(I)   = RH(I)   - HCP*H(I)*H(J)*DVOL*ETEMP(J)/DTINC
	ENDIF

	ENDDO
	ENDDO

 900  CONTINUE
 	
C	==========================================
	LFG = TMP(2)
	IF(LFG.EQ.0) GOTO 1300


	MG1   = 2
	DO 1200 IFACE = 1,NFACE
	IF(TCD(IFACE).EQ.0.0) GOTO 1200

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


	DO I = 1,NNM
	DO J = 1,NNM
	SH(I,J)  = SH(I,J) + FAC1*HCV*H(I)*H(J)*DFACE*TCD(IFACE)
	IF(KSTEAD.EQ.1) THEN
	RH(I)    = RH(I)   + FAC2*HCV*H(I)*H(J)*DFACE*ETEMP(J)*TCD(IFACE)
	ENDIF
	ENDDO
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
	DO I = 1,NNM
	DO J = 1,NNM
	SH(I,J)  = SH(I,J) + FAC1*HCV*H(I)*H(J)*DFACE*TCD(IFACE)
	IF(KSTEAD.EQ.1) THEN
	RH(I)    = RH(I)   + FAC2*HCV*H(I)*H(J)*DFACE*ETEMP(J)*TCD(IFACE)
	ENDIF
	ENDDO
	ENDDO
1100	CONTINUE
C	-------------------------------------------------------
	END SELECT

1200	CONTINUE


C	==============================
C	TRANSFER LOCAL HEAT CONDUCTION 
C	TO GLOBAL ELEMENT STIFFNESS
C	==============================
1300	CONTINUE

	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	RE(II) = RE(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE SHAP1D(R,H,P,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------------------------------------------------
      DIMENSION H(NNO),P(NNO)

	H = 0.0
	P = 0.0

	SELECT CASE (NNO)
	
	CASE(2)
	H(1) = 0.5*(1.0 + R)
	H(2) = 0.5*(1.0 - R)
	P(1) = 0.5
	P(2) =-0.5

	CASE(3)
	H(1) =  0.5*R + 0.5*R*R
	H(2) =  1.0   -     R*R
	H(3) = -0.5*R + 0.5*R*R
	P(1) =  0.5   + 1.0*R
	P(2) =        - 2.0*R
	P(3) = -0.5   + 1.0*R


	END SELECT
	

	RETURN

	END


C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE JACO1H(IFACE,XY,P,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(3,NNO),P(NNO),XJ(3),NN(3)
C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R + S = 1.0
C	FACE 2; R= 0.0
C	FACE 3: S= 0.0
C	---------------------

	SELECT CASE(IFACE)
	CASE(2)
	IF(NNO.EQ.2) NN(1:2) = [3,2]
	CASE(3)
	IF(NNO.EQ.2) NN(1:2) = [3,1]
	CASE(1)
	IF(NNO.EQ.2) NN(1:2) = [2,1] 
	END SELECT
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  J=1,3
      DUM = 0.0
      DO 90   K=1,NNO
	KK  = NN(K)
 90   DUM = DUM + P(K)*XY(J,KK)
 100  XJ(J) = DUM
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = SQRT(XJ(1)*XJ(1) + XJ(2)*XJ(2) + XJ(3)*XJ(3)) 
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')   
C
      RETURN
      END


C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE SHAP2DT(R,S,H,P,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION  H(NNO),P(2,NNO)
C
      H = 0.0D0
	P = 0.0D0
	
	X1 = 1.0-R-S
	X2 = R
	X3 = S
	P1R = -1.0D0
	P1S = -1.0D0
	P2R =  1.0D0
	P2S =  0.0D0
	P3R =  0.0D0
	P3S =  1.0D0

	SELECT CASE (NNO)

	CASE(3)
	H(1) = X1
	H(2) = X2
	H(3) = X3
	P(1,1) =  P1R
	P(2,1) =  P1S
	P(1,2) =  P2R
	P(2,2) =  P2S
	P(1,3) =  P3R
	P(2,3) =  P3S

	CASE(6)
	H(1) = X1*(2.0*X1-1.0)
	H(2) = X2*(2.0*X2-1.0)
	H(3) = X3*(2.0*X3-1.0)
	H(4) = 4.0*X1*X2
	H(5) = 4.0*X2*X3
	H(6) = 4.0*X1*X3

	P(1,1) = (4.0*X1-1.0)*P1R
	P(1,2) = (4.0*X2-1.0)*P2R
	P(1,3) = (4.0*X3-1.0)*P3R
	P(1,4) = 4.0*X1*P2R + 4.0*X2*P1R
	P(1,5) = 4.0*X2*P3R + 4.0*X3*P2R
	P(1,6) = 4.0*X1*P3R + 4.0*X3*P1R

	P(2,1) = (4.0*X1-1.0)*P1S
	P(2,2) = (4.0*X2-1.0)*P2S
	P(2,3) = (4.0*X3-1.0)*P3S
	P(2,4) = 4.0*X1*P2S + 4.0*X2*P1S
	P(2,5) = 4.0*X2*P3S + 4.0*X3*P2S
	P(2,6) = 4.0*X1*P3S + 4.0*X3*P1S

	END SELECT

      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE JACO2HT(XY,P,XJ,XJI,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(3,NNO),P(2,NNO),XJ(2,2),XJI(4)
	DIMENSION RS(2,NNO),VR(3),VS(3),VT(3)

	VR(1) = XY(1,2) - XY(1,1)
	VR(2) = XY(2,2) - XY(2,1)
	VR(3) = XY(3,2) - XY(3,1)
	CALL SCALEN(VR,VR,DUM,3)
	VS(1) = XY(1,3) - XY(1,1)
	VS(2) = XY(2,3) - XY(2,1)
	VS(3) = XY(3,3) - XY(3,1)
	CALL SCALEN(VS,VS,DUM,3)
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
	CALL VECPRD(VT,VR,VS)

	DO I = 1,NNO
	RS(1,I) = XY(1,I)*VR(1) + XY(2,I)*VR(2) + XY(3,I)*VR(3) 
	RS(2,I) = XY(1,I)*VS(1) + XY(2,I)*VS(2) + XY(3,I)*VS(3) 
	ENDDO
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  I=1,2
      DO 100  J=1,2
      DUM = 0.0
      DO 90   K=1,NNO
 90   DUM = DUM + P(I,K)*RS(J,K)
 100  XJ(I,J) = DUM
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = XJ(1,1)*XJ(2,2) - XJ(2,1)*XJ(1,2)
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')
C     -----------------------------
C     INVERSE (XJI) OF THE JACOBIAN
C     -----------------------------
      DUM = 1.0/DET
      XJI(1) =  XJ(2,2)*DUM
      XJI(2) = -XJ(2,1)*DUM
      XJI(3) = -XJ(1,2)*DUM
      XJI(4) =  XJ(1,1)*DUM

C	MODIFIED FOR TRIANGULAR SHAPE
	DET = 0.5*DET      
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE GAUSST (R,S,T,W,IPT,NPT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION  RI(4),SI(4),TI(4),WI(4)

	SELECT CASe(IND)

	CASE(0)

	IF(NPT.EQ.1) THEN
	RI(1) = 1.0/3.0
	SI(1) = 1.0/3.0	
	WI(1)  = 1.0
	ELSEIF(NPT.EQ.3) THEN
	RI(1) = 2.0/3.0
	SI(1) = 1.0/6.0
	RI(2) = 1.0/6.0
	SI(2) = 1.0/6.0
	RI(3) = 1.0/6.0
	SI(3) = 2.0/3.0
	WI(1:3)= 1.0/3.0
	ELSEIF(NPT.EQ.4) THEN
	RI(1) = 3.0/5.0
	SI(1) = 1.0/5.0
	RI(2) = 1.0/5.0
	SI(2) = 1.0/5.0
	RI(3) = 1.0/5.0
	SI(3) = 3.0/5.0
	RI(4) = 1.0/3.0
	SI(4) = 1.0/3.0
	WI(1:3)= 25.0/48.0
	WI(4)  =-27.0/48.0
	ENDIF

	CASE(1)

	IF(NPT.EQ.1) THEN
	RI(1) = 1.0/4.0
	SI(1) = 1.0/4.0	
	TI(1) = 1.0/4.0
	WI(1) = 1.0D0
	ELSEIF(NPT.EQ.4) THEN
	A = 5.0 + 3.0*SQRT(5.0)
	A = A/20.0
	B = 5.0 - 1.0*SQRT(5.0)
	B = B/20.0
	RI(1) = A
	SI(1) = B
	TI(1) = B
	RI(2) = B
	SI(2) = B
	TI(2) = B
	RI(3) = B
	SI(3) = B
	TI(3) = A
	RI(4) = B
	SI(4) = A
	TI(4) = B
	WI(1:4)= 1.0/4.0
	ELSEIF(NPT.EQ.5) THEN
	A = 0.25
	B = 0.50
	C = 1.0/6.0
	RI(1) = A
	SI(1) = A
	TI(1) = A

	RI(2) = B
	SI(2) = C
	TI(2) = C

	RI(3) = C
	SI(3) = C
	TI(3) = C

	RI(4) = C
	SI(4) = C
	TI(4) = B

	RI(5) = C
	SI(5) = B
	TI(5) = C
	WI(1)  = -4.0/5.0
	WI(2:5) = 9.0/20.0
	ENDIF

	END SELECT

	R = RI(IPT)
	S = SI(IPT)
	T = TI(IPT)
	W = WI(IPT)

      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE RESHPT(H,IFACE,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------------------------------------------------
      DIMENSION H(NNO),HH(NNO)

C	---------------------
C	DEFINE EACH SURFACE
C	FACE 1; R + S = 1.0
C	FACE 2; R= 0.0
C	FACE 3: S= 0.0
C	---------------------

	HH = 0.0

	SELECT CASE(IFACE)

	CASE(2)
	IF(NNO.EQ.3) THEN
	HH(3) = H(1)
	HH(2) = H(2)
	ENDIF

	CASE(3)
	IF(NNO.EQ.3) THEN
	HH(3) = H(1)
	HH(1) = H(2)
	ENDIF

	CASE(1)
	IF(NNO.EQ.3) THEN
	HH(2) = H(1)
	HH(1) = H(2)
	ENDIF

	END SELECT

	H(1:NNO) = HH(1:NNO)
	

	RETURN

	END



C	=====================================================================
C	=====================================================================	
C	=====================================================================

C	=======================================================================
C	==============START OF 4 NODE TETRAHEDRA HEAT ELEMENT==================
C	=======================================================================
	SUBROUTINE HEATELE3T(COORD,PROPM,EDIS,S,RE,TCD,TMP,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

C	PARAMETER FOR TRANSIENT HEAT tRANSFER SONGSAK MAR2007
	COMMON /SHEAT1/ PIMPLT,KSTEAD

	PARAMETER (NFACE=4)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM),TMP(2)
	DIMENSION HM(3)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1),TCD(NFACE)
	DIMENSION BH(3,NNM),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF TEMPERATURE DOF

	FAC1 = PIMPLT
	FAC2 = 1.0 - FAC1
C	TRANSFER NODAL VALUE INTO TEMPERATURE PART
	DO I = 1,NNM
	II = (NMCHA+NMTEM)*I
	ETEMP(I) = EDIS(II)
	ENDDO

	S  = 0.0D0
	RE = 0.0D0

	SH = 0.0D0
	RH = 0.0D0

C	==========================================
C	SETTING ELEMENT GAUSSIAN INTEGRATION POINT
C	STANDARD GAUSS INTEGRATION FOR 8 NODES SOLID ELEMENT
C	(MGR,MGS,MGT) = (2,2,2)
C	=======================

      MGR = 4
C	===============================================
C	CALL HEAT CONDUCTION MATRIX (CONSTITUTIVE HEAT)
C	===============================================
	DO I = 1,3
	HM(I)  = PROPM(5+I)
	ENDDO

C	HEAT CAPACITY COEFICIENT
	HCP = PROPM(11)

C	HEAT CONDUCTION COEFICIENT
	HCV = PROPM(12)

	IPT = 0
C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS
C	=========================================
      DO 900  IGR=1,MGR
	IPT = IPT +1 

	
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
C	COMPUTE HEAT CONDUCTION STIFFNESS MATRIX
C	========================================
	BH = MATMUL(XJI,P)


	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,3
	SH(I,J) = SH(I,J) + FAC1*HM(K)*BH(K,I)*BH(K,J)*DVOL

	IF(KSTEAD.EQ.1) THEN
	RH(I)   = RH(I)   + FAC2*HM(K)*BH(K,I)*BH(K,J)*DVOL*ETEMP(J)
	ENDIF

	ENDDO
	ENDDO
	ENDDO

	DO I = 1,NNM
	DO J = 1,NNM

	IF(KSTEAD.EQ.1) THEN
	SH(I,J) = SH(I,J) + HCP*H(I)*H(J)*DVOL/DTINC
	RH(I)   = RH(I)   - HCP*H(I)*H(J)*DVOL*ETEMP(J)/DTINC
	ENDIF

	ENDDO
	ENDDO

 900  CONTINUE
 	
C	WRITE(*,*) 'KAE',ETEMP
	
C	==========================================
	LFG = TMP(2)
	IF(LFG.EQ.0) GOTO 1200

	MG1   = 3

	NN(4,1:6) = [2,3,4,8,9,10]
	NN(1,1:6) = [1,3,4,6,9,7]
	NN(2,1:6) = [1,2,4,5,10,7]
	NN(3,1:6) = [1,2,3,5,8,6]


	DO 1100 IFACE = 1,NFACE
	
	IF(TCD(IFACE).EQ.0.0) GOTO 1100

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
	
	DO I = 1,NNM
	DO J = 1,NNM
	SH(I,J)  = SH(I,J) + FAC1*HCV*H(I)*H(J)*DFACE*TCD(IFACE)
	IF(KSTEAD.EQ.1) THEN
	RH(I)    = RH(I)   + FAC2*HCV*H(I)*H(J)*DFACE*ETEMP(J)*TCD(IFACE)
	ENDIF

	ENDDO
	ENDDO

1000	CONTINUE

1100	CONTINUE


C	==============================
C	TRANSFER LOCAL HEAT CONDUCTION 
C	TO GLOBAL ELEMENT STIFFNESS
C	==============================
1200	CONTINUE

	DO I = 1,NNM             !FOR TEMPERATURE DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	RE(II) = RE(II) + RH(I)
	ENDDO
C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SHAP3DT(R,S,T,H,P,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION  H(NNO),P(3,NNO)
C
      H = 0.0D0
	P = 0.0D0
	
	X1 = 1.0-R-S-T
	X2 = R
	X3 = S
	X4 = T
	P1R = -1.0D0
	P1S = -1.0D0
	P1T = -1.0D0
	P2R =  1.0D0
	P2S =  0.0D0
	P2T =  0.0D0
	P3R =  0.0D0
	P3S =  1.0D0
	P3T =  0.0D0
	P4R =  0.0D0
	P4S =  0.0D0
	P4T =  1.0D0

	SELECT CASE (NNO)

	CASE(4)
	H(1) = X1
	H(2) = X2
	H(3) = X3
	H(4) = X4
	P(1,1) =  P1R
	P(2,1) =  P1S
	P(3,1) =  P1T
	P(1,2) =  P2R
	P(2,2) =  P2S
	P(3,2) =  P2T
	P(1,3) =  P3R
	P(2,3) =  P3S
	P(3,3) =  P3T
	P(1,4) =  P4R
	P(2,4) =  P4S
	P(3,4) =  P4T

	CASE(10)
	H(1) = X1*(2.0*X1-1.0)
	H(2) = X2*(2.0*X2-1.0)
	H(3) = X3*(2.0*X3-1.0)
	H(4) = X4*(2.0*X4-1.0)

	H(5) = 4.0*X1*X2
	H(6) = 4.0*X1*X3
	H(7) = 4.0*X1*X4

	H(8) = 4.0*X2*X3
	H(9) = 4.0*X3*X4
	H(10)= 4.0*X2*X4

	P(1,1) = (4.0*X1-1.0)*P1R
	P(1,2) = (4.0*X2-1.0)*P2R
	P(1,3) = (4.0*X3-1.0)*P3R
	P(1,4) = (4.0*X4-1.0)*P4R

	P(1,5) = 4.0*X1*P2R + 4.0*X2*P1R
	P(1,6) = 4.0*X1*P3R + 4.0*X3*P1R
	P(1,7) = 4.0*X1*P4R + 4.0*X4*P1R

	P(1,8) = 4.0*X2*P3R + 4.0*X3*P2R
	P(1,9) = 4.0*X3*P4R + 4.0*X4*P3R
	P(1,10)= 4.0*X2*P4R + 4.0*X4*P2R

	P(2,1) = (4.0*X1-1.0)*P1S
	P(2,2) = (4.0*X2-1.0)*P2S
	P(2,3) = (4.0*X3-1.0)*P3S
	P(2,4) = (4.0*X4-1.0)*P4S

	P(2,5) = 4.0*X1*P2S + 4.0*X2*P1S
	P(2,6) = 4.0*X1*P3S + 4.0*X3*P1S
	P(2,7) = 4.0*X1*P4S + 4.0*X4*P1S

	P(2,8) = 4.0*X2*P3S + 4.0*X3*P2S
	P(2,9) = 4.0*X3*P4S + 4.0*X4*P3S
	P(2,10)= 4.0*X2*P4S + 4.0*X4*P2S

	P(3,1) = (4.0*X1-1.0)*P1T
	P(3,2) = (4.0*X2-1.0)*P2T
	P(3,3) = (4.0*X3-1.0)*P3T
	P(3,4) = (4.0*X4-1.0)*P4T

	P(3,5) = 4.0*X1*P2T + 4.0*X2*P1T
	P(3,6) = 4.0*X1*P3T + 4.0*X3*P1T
	P(3,7) = 4.0*X1*P4T + 4.0*X4*P1T

	P(3,8) = 4.0*X2*P3T + 4.0*X3*P2T
	P(3,9) = 4.0*X3*P4T + 4.0*X4*P3T
	P(3,10)= 4.0*X2*P4T + 4.0*X4*P2T

	END SELECT

      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================




C	=======================================================================
C	==============START OF 4 NODE TETRAHEDRA SOLID ELEMENT=================
C	=======================================================================
	SUBROUTINE SOLIDEL3T(PROPM,PROPG,WA,S,COORD,EDIS,EDISI,RE,MWG,
     +			  	     FIN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK


      DIMENSION PROPM(1),PROPG(1),WA(MWG,1),S((NEF*NEF+NEF)/2)
      DIMENSION EDIS(NEF),EDISI(NEF),RE(NEF),SG((NEF*NEF+NEF)/2)
 
      DIMENSION H(NNM),P(3,NNM),XJI(3,3),B(3*NNM),DISD(9)
	DIMENSION BM(6,NEF),DE(6,6),DEP(6,6)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6)
	DIMENSION SK(NEF,NEF)

	DIMENSION FIN(NEF),COORD(3,NNM)
	DIMENSION DFVP(NEF)
	DIMENSION XJ(3,3)


	CALL HOKLAW_S (PROPM,PROPG,1)

	SG   = 0.0D0
	SK   = 0.0D0
	RE   = 0.0D0
	DFVP = 0.0D0

      MGR = 4
	IF(NNM.EQ.10) MGR = 4
	IPT = 0
	

C     =========================================
C     LOOP OVER VOLUME INTEGRATION GAUSS POINTS
C	=========================================
      DO 900  IGR=1,MGR
	IPT = IPT +1 

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


	IF (ITASK.NE.5)  GOTO 50
      CALL SOMASS_S (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
      GOTO 900

50    CONTINUE
C	========================================
C	COMPUTE STRAIN-DISPLACEMENT MATRIX
C	========================================
	CALL BMATSLT (P,XJI,B,BM,NNO)
	CALL DMATSLD_S (DE)

	
	IF (NLOPT+ITASK.EQ.1)  GOTO 700
C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------
      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,6
100   TAU(I) = WA(I,IPT)*DVOL
      GOTO 800

	
C     -------------------------------------
C      COMPUTE DISPLACEMENT GRADIENT (DISD)
C     -------------------------------------

 200  CONTINUE

	DISD = 0.0D0

      DO IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDIS(IEF)   ! H,x *u
      DISD(2) = DISD(2) + B(JEF)*EDIS(JEF)   ! H,y *v
      DISD(3) = DISD(3) + B(KEF)*EDIS(KEF)   ! H,z *w
      DISD(4) = DISD(4) + B(JEF)*EDIS(IEF)   ! H,y *u
      DISD(5) = DISD(5) + B(KEF)*EDIS(IEF)   ! H,z *u
      DISD(6) = DISD(6) + B(IEF)*EDIS(JEF)   ! H,x *v 
      DISD(7) = DISD(7) + B(KEF)*EDIS(JEF)   ! H,z *v
      DISD(8) = DISD(8) + B(IEF)*EDIS(KEF)   ! H,x *w
      DISD(9) = DISD(9) + B(JEF)*EDIS(KEF)   ! H,y *w
	END DO

C     ------------------------------
C     LINEAR COMPATIBLE STRAIN TERMS
C     ------------------------------
      STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)
      STRAIN(4) = DISD(4) + DISD(6)
      STRAIN(5) = DISD(5) + DISD(8)
      STRAIN(6) = DISD(7) + DISD(9)

C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT QUADRATIC STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400
      QSTRAI(1) = .5*(DISD(1)*DISD(1)+DISD(6)*DISD(6)+DISD(8)*DISD(8))
      QSTRAI(2) = .5*(DISD(4)*DISD(4)+DISD(2)*DISD(2)+DISD(9)*DISD(9))
      QSTRAI(3) = .5*(DISD(5)*DISD(5)+DISD(7)*DISD(7)+DISD(3)*DISD(3))
      QSTRAI(4) =     DISD(1)*DISD(4)+DISD(6)*DISD(2)+DISD(8)*DISD(9)
      QSTRAI(5) =     DISD(1)*DISD(5)+DISD(6)*DISD(7)+DISD(8)*DISD(3)
      QSTRAI(6) =     DISD(4)*DISD(5)+DISD(2)*DISD(7)+DISD(9)*DISD(3)

C	----------------------
C	SUBTRACT ALMASI STRAIN
C	----------------------
      DO 390  I=1,6
 390  STRAIN(I) = STRAIN(I) - QSTRAI(I)

C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------

400   GOTO (410,420,430,440,450,460,470,470),MTMOD

410   CONTINUE

C	----------------------------------------------------------
C	LINEAR ELASTIC HOOKE'S CONSTITUTUVE STRESS-STRAIN RELATION
C	----------------------------------------------------------

	CALL SOLSIG_S (STRAIN,STRESS)
      DO 415  I=1,6
415   WA(I,IPT) = STRESS(I)
      GOTO 500

420	CONTINUE
	
	WRITE(*,*) 'MATERIAL NUMBER IS NOT AVAILIBLE NOW'
	STOP

	GOTO 500


430	CONTINUE

	CALL DMMISE (WA(1,IPT),WA(7,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DEP)

      GOTO 500


440   CONTINUE

C	CALL MOHRCRI3D(WA(1,IPT),WA(7,IPT),STRAIN,STRESS,DEP,
C	1			   PROPM,WA(13,IPT))
	CALL MOHR3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			 STRAIN,DEP,PROPM)

	GO TO 500

450   CONTINUE

	CALL DRUG3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			 STRAIN,DEP,PROPM,6,0)

      GOTO 500


460   CONTINUE

	CALL DAMAGE(STRAIN,STRESS,DEP,PROPM,WA(1,IPT))

      GOTO 500


470	CONTINUE 

	CALL PLVIS(WA(1,IPT),WA(7,IPT),WA(14,IPT),
	1		   DFVP,STRESS,STRAIN,PROPM,
     2		   BM,DVOL,WA(13,IPT),DEP)

	GO TO 500
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------

 500  CONTINUE

	DO I=1,6
      TAU(I) = STRESS(I)*DVOL
	END DO

	IF (ITASK.LE.3) GOTO 520
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900
 520	CONTINUE  

C     ---------------------------------------
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     ---------------------------------------
	RE = RE + MATMUL(TRANSPOSE(BM),STRESS)*DVOL

	IF (IFEIG.EQ.0.AND.ISOLOP.EQ.4) GOTO 800
C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------
 151  IF (IFREF) 900,700,900
 700  IF (NLOPT.NE.0) GOTO 750

	
C	****************
C	LINEAR ANALYSIS
C	****************
	DEP = 0.
	SK  = 0.

	
	CALL DMATSLD_S (DEP)



	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL

	  KSK = 0
	    DO ISK =1,NEF
	      DO JSK =ISK,NEF
	         KSK = KSK+1 
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO



	GOTO 790

C	*************************
C	...GEOMETRIC NONLINEARITY
C	*************************

 750  IF (MTMOD.LE.2)THEN

	DEP =0.
	SK = 0.

	CALL DMATSLD_S (DEP)

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,NEF
	      DO JSK =ISK,NEF
	         KSK = KSK+1 
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	END IF	

C	************************
C	...MATERIAL NONLINEARITY
C	************************
	IF (MTMOD.GT.2)THEN

	IF(NLOPT+ITASK.EQ.1) THEN


	IF (MTMOD.EQ.3) THEN

	DEP =0.

	CALL DMMISE (WA(1,IPT),WA(7,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DEP)

	ENDIF


	IF (MTMOD.EQ.4) THEN

	DEP =0.0

	CALL MOHRCRI3D(WA(1,IPT),WA(7,IPT),STRAIN,STRESS,DEP,
	1			   PROPM,WA(13,IPT))

	ENDIF

	IF (MTMOD.EQ.5) THEN

	DEP =0.

	CALL DRUG3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			STRAIN,DEP,PROPM,6,0)

	ENDIF


	IF (MTMOD.EQ.6) THEN

	DEP =0.

	CALL DAMAGE(STRAIN,STRESS,DEP,PROPM,WA(1,IPT))

	ENDIF


	IF (MTMOD.EQ.7) THEN

	DEP =0.

	CALL PLVIS(WA(1,IPT),WA(7,IPT),WA(14,IPT),
	1		   DFVP,STRESS,STRAIN,PROPM,
     2		   BM,DVOL,WA(13,IPT),DEP)

	ENDIF


	ENDIF ! <--- ENDIF NLOPT+ITASK.EQ.1


	SK = 0.

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,NEF
	      DO JSK =ISK,NEF
	         KSK = KSK+1 
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	END IF


C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------

790   IF (NLOPT.LE.1) GOTO 810
800	CONTINUE

C	------------------------------------------------
C	ADDITIONAL COMPUTE NONLINEAR GEOMATRIC STIFFNESS
C	------------------------------------------------
	CALL KNLSTIFT (SG,TAU,B,NNM,NEF)

 810  CONTINUE 


 900  CONTINUE


	IF(ITASK.EQ.5) GOTO 150   !MASS MATRIX
	IF(IFEIG.EQ.0.AND.ISOLOP.EQ.4) GOTO 950

950	KGG = 0
	DO I = 1,NEF
	DO J = I,NEF
	KGG = KGG + 1
	S(KGG) = S(KGG) + SG(KGG)
	ENDDO
	ENDDO
	 
C     ------------------------------------
C     FOR NONLINEAR ANALYSIS OF EAS METHOD
C     ------------------------------------
	IF (ITASK+NLOPT.EQ.1) GOTO 150

	DO I=1,NEF
	 RE(I) = RE(I) + DFVP(I)
	END DO

150	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF

c	KGG = 0
c	DO I = 1,NEF
c	DO J = I,NEF
c	KGG = KGG + 1
c	IF(I.EQ.J) WRITE(*,*) I,NEF,S(KGG)
c	ENDDO
c	ENDDO
c	PAUSE


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE BMATSLT (P,XJI,B,BM,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION P(3,NNO),XJI(3,3),B(3*NNO),BM(6,3*NNO)


      CALL CLEARA (B ,3*NNO)
	CALL CLEARA (BM,6*3*NNO)

C	------------------------------------------------
C	 COMPUTE COMPACTED STRESS-STRAIN GRADIENT VECTOR FOR GEOMETRIC STIFFNESS CALCULATION
C	------------------------------------------------
      K = 1
	KK= 0
      DO 100  INO=1,NNO
      DO 50   I=1,3
      B(K)   = B(K)   + XJI(1,I)*P(I,INO)
      B(K+1) = B(K+1) + XJI(2,I)*P(I,INO)
      B(K+2) = B(K+2) + XJI(3,I)*P(I,INO)
 50   CONTINUE
      K = K+3
	KK=KK+3
 100  CONTINUE

C	----------------------------------------
C	 COMPUTE STRESS-STRAIN MATRIX
C	----------------------------------------
	DO INO=1,NNO
	HX =0.0D0
	HY =0.0D0
	HZ =0.0D0
	DO I=1,3
        HX  = HX  + XJI(1,I)*P(I,INO)
        HY  = HY  + XJI(2,I)*P(I,INO)
        HZ  = HZ  + XJI(3,I)*P(I,INO)
	END DO

C	--------------------------
C  	 COMPUTE MATRIX [B](6,24)
C	--------------------------

	BM(1,3*INO-2) = HX
	BM(2,3*INO-1) = HY
	BM(3,3*INO-0) = HZ

	BM(4,3*INO-2) = HY
	BM(4,3*INO-1) = HX

	BM(5,3*INO-2) = HZ
	BM(5,3*INO-0) = HX

	BM(6,3*INO-1) = HZ
	BM(6,3*INO-0) = HY

	END DO

      RETURN
      END
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE KNLSTIFT (S,SIG,B,NNM,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
		
C     ----------------------------------------------------------------
C     ADDS INITIAL STRESS CONTRIBUTION TO ELEMENT STIFFNESS MATRIX
C	------------------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     S(NWS)    = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     SIG(6)    = GAUSS POINT STRESSES (MULTIPLIED BY DVOL=WT*DET)
C     B(3*NNO)  = COMPACTED STRAIN DISPLACEMENT MATRIX
C     NEF       = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ----------------------------------------------------------------

      DIMENSION S(1),SIG(6),B(3*NNM)
	DIMENSION Q(9,9),G(9,NEF),SNL(NEF,NEF)


C	...SETTING MATRIX [G]

	G   = 0.
	Q   = 0.
	SNL = 0.

	DO I =1,NNM

	G(1,3*I-2) = B(3*I-2) !Hx[I]
	G(2,3*I-1) = B(3*I-2) !Hx[I]
	G(3,3*I-0) = B(3*I-2) !Hx[I]

	G(4,3*I-2) = B(3*I-1) !Hy[I]
	G(5,3*I-1) = B(3*I-1) !Hy[I]
	G(6,3*I-0) = B(3*I-1) !Hy[I]

	G(7,3*I-2) = B(3*I-0) !Hz[I]
	G(8,3*I-1) = B(3*I-0) !Hz[I]
	G(9,3*I-0) = B(3*I-0) !Hz[I]

	END DO

C	...SETTING MATRIX [Q]

      Q(1,1)=SIG(1)
	Q(2,2)=SIG(1)
	Q(3,3)=SIG(1)

      Q(4,4)=SIG(2)
	Q(5,5)=SIG(2)
	Q(6,6)=SIG(2)

	Q(7,7)=SIG(3)
	Q(8,8)=SIG(3)
	Q(9,9)=SIG(3)

	Q(4,1)=SIG(4)
	Q(5,2)=SIG(4)
	Q(6,3)=SIG(4)

	Q(7,1)=SIG(5)
	Q(8,2)=SIG(5)
	Q(9,3)=SIG(5)

	Q(7,4)=SIG(6)
	Q(8,5)=SIG(6)
	Q(9,6)=SIG(6)

C	FILLED FOR SYMMETRY

	DO I=1,9
	  DO J =I,9
	    Q(I,J)=Q(J,I)
	  END DO
	END DO


C	COMPUTE NON LINEAR GEOMETRIC STIFFNESS

	SNL = MATMUL(TRANSPOSE(G),MATMUL(Q,G))

	K = 0
	DO IPT =1,NEF
	 DO JPT =IPT,NEF 
	     K = K+1
	  S(K) = S(K)+SNL(IPT,JPT)
	END DO
	END DO


	RETURN 
	END


C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE TRNMAT2Q(XY,P,NNO,TMAT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------
      DIMENSION XY(3,NNO),P(2,NNO),VR(3),VS(3),VT(3),TMAT(2,2)
      DIMENSION COVR(3),COVS(3),RV(3),SV(3)
C
      CALL CLEARA (COVR,3)
      CALL CLEARA (COVS,3)
      DO 20 I=1,NNO
      DO 20 J=1,3
      COVR(J)=COVR(J)+P(1,I)*XY(J,I)
   20 COVS(J)=COVS(J)+P(2,I)*XY(J,I)

      CALL VECPRD (COVR,COVS,VT)
      CALL SCALEN (VT,VT,DET,3)

      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)
      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)


      TMAT(1:2,1) = VR(1:2)
      TMAT(1:2,2) = VS(1:2)

      RETURN

      END

C	=====================================================================
C	=====================================================================	
C	=====================================================================
	SUBROUTINE TRNMAT2T(XY,P,NNO,TMAT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(3,NNO),P(2,NNO),TMAT(2,2)
	DIMENSION VR(3),VS(3),VT(3)

	VR(1) = XY(1,2) - XY(1,1)
	VR(2) = XY(2,2) - XY(2,1)
	VR(3) = XY(3,2) - XY(3,1)
	CALL SCALEN(VR,VR,DUM,3)
	VS(1) = XY(1,3) - XY(1,1)
	VS(2) = XY(2,3) - XY(2,1)
	VS(3) = XY(3,3) - XY(3,1)
	CALL SCALEN(VS,VS,DUM,3)
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
	CALL VECPRD(VT,VR,VS)

      TMAT(1:2,1) = VR(1:2)
      TMAT(1:2,2) = VS(1:2)
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================	
C	=====================================================================

