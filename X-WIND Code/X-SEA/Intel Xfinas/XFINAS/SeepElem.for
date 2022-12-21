C	SOLID SEEPAGE LIBRARY
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SEEPELE3D(COORD,PROPM,EDIS,S,RE,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION COORD(1),PROPM(1),EDIS(1),RE(1)
	DIMENSION WA(1)

	SELECT CASE(ISTYP)

	CASE(3)
	CALL  SEEPELE3(COORD,PROPM,EDIS,S,RE,WA)

	CASE(4)
	CALL SEEPELE3T(COORD,PROPM,EDIS,S,RE,WA)

	END SELECT


	RETURN
	END


C	=======================================================================
C	======================START OF 8 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE SEEPELE3(COORD,PROPM,EDIS,S,RE,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	COMMON /MGRAV/ NGRAV

	PARAMETER (NFACE=6)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM)
	DIMENSION HM(3),HMP(3)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1)
	DIMENSION BH(3,NNM),NODEX(13),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)


      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF FLOW DOF


C	TRANSFER NODAL VALUE INTO FLOW PART
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
C	CALL CONDUCTIVITY MATRIX (CONSTITUTIVE POTENTIAL)
C	===============================================
	DO I = 1,3
	HM(I)  = PROPM(5+I)
	HMP(I) = HM(I)
	ENDDO

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

C	---------------------------------------
	IF(KFRES.EQ.1.AND.KSTEP.GT.1) THEN  !KFRES = FREE SURFACE ITERATION OPTION FOR 1
	GLEV = 0.0D0
	PTEN = 0.0D0
	DO INM = 1,NNM
	GLEV = GLEV + H(INM)*COORD(NGRAV,INM)
	PTEN = PTEN + H(INM)*EDIS(4*INM)
	ENDDO
	IF(GLEV.GT.PTEN) HM(1:3) = HMP(1:3)*1.0E-6
	ENDIF

	DO ICOM = 1,3
	WA(ICOM,IPT) = 0.0D0
	DO INM = 1,NNM
	WA(ICOM,IPT) = WA(ICOM,IPT) - HM(ICOM)*BH(ICOM,INM)*EDIS(4*INM)
	ENDDO
	ENDDO
C	---------------------------------------

	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,3
	SH(I,J) = SH(I,J) + HM(K)*BH(K,I)*BH(K,J)*DVOL



	ENDDO
	ENDDO
	ENDDO


 900  CONTINUE
C	====================================================
C	TRANSFER LOCAL STIFFNESS TO GLOBAL ELEMENT STIFFNESS
C	====================================================

	DO I = 1,NNM             !FOR FLOW DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	ENDDO

C	====================================================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================


C	MEMBRANE SEEPAGE LIBRARY
C	=======================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE SEEPELE2D(COORD,PROPM,PROPG,EDIS,S,RE,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION COORD(1),PROPM(1),PROPG(1),EDIS(1),RE(1),TCD(1),TMP(1) 
	DIMENSION WA(1)

	SELECT CASE(ISTYP)

	CASE(3)
	CALL  SEEPELE2(COORD,PROPM,PROPG,EDIS,S,RE,WA)

	CASE(4)
	CALL SEEPELE2T(COORD,PROPM,PROPG,EDIS,S,RE,WA)

	END SELECT


	RETURN
	END


C	=======================================================================
C	======================START OF 4 NODE HEAT ELEMENT=====================
C	=======================================================================
	SUBROUTINE SEEPELE2(COORD,PROPM,PROPG,EDIS,S,RE,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	COMMON /MGRAV/ NGRAV

	PARAMETER (NFACE=6)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM)
	DIMENSION HM(2),HMP(2)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1),PROPG(1)
	DIMENSION BH(2,NNM),NODEX(4),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)
	DIMENSION VXY(2),TMAT(2,2) !FOR TRANSFORMATION OF FLOW VELOCITY OF SEEPAGE

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF FLOW DOF

C	TRANSFER NODAL VALUE INTO FLOW PART
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
	HMP(I) = HM(I)
	ENDDO

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

C	---------------------------------------
	IF(KFRES.EQ.1.AND.KSTEP.GT.1) THEN  !KFRES = FREE SURFACE ITERATION OPTION FOR 1
	GLEV = 0.0D0
	PTEN = 0.0D0
	DO INM = 1,NNM
	GLEV = GLEV + H(INM)*COORD(NGRAV,INM)
	PTEN = PTEN + H(INM)*EDIS(3*INM)
	ENDDO
	IF(GLEV.GT.PTEN) HM(1:2) = HMP(1:2)*1.0E-6
	ENDIF

	CALL TRNMAT2Q(COORD,P,NNM,TMAT)
	DO ICOM = 1,2
	VXY(ICOM) = 0.0D0
	DO INM = 1,NNM
	VXY(ICOM) = VXY(ICOM) - HM(ICOM)*BH(ICOM,INM)*EDIS(3*INM)
	ENDDO
	ENDDO
	DO IRS = 1,2
	WA(IRS,IPT) = 0.0D0
	DO JRS = 1,2
	WA(IRS,IPT) = WA(IRS,IPT) + TMAT(IRS,JRS)*VXY(JRS)
	ENDDO
	ENDDO
C	---------------------------------------

	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,2
	SH(I,J) = SH(I,J) + HM(K)*BH(K,I)*BH(K,J)*DVOL
	ENDDO
	ENDDO
	ENDDO


 900  CONTINUE

C	============================================================
C	TRANSFER LOCAL STIFFNESS TO GLOBAL ELEMENT STIFFNESS
C	============================================================
	DO I = 1,NNM             !FOR FLOW DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	ENDDO
C	==============================


	RETURN
	END


C	=====================================================================
C	==================START OF 3 NODE SEEPAGE ELEMENT====================	
C	=====================================================================
	SUBROUTINE SEEPELE2T(COORD,PROPM,PROPG,EDIS,S,RE,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES


	COMMON /MGRAV/ NGRAV

	PARAMETER (NFACE=5)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM)
	DIMENSION HM(2),HMP(2)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1),PROPG(1)
	DIMENSION BH(2,NNM),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(2,NNM),XJ(2,2),XJI(2,2)
	DIMENSION VXY(2),TMAT(2,2) !FOR TRANSFORMATION OF FLOW VELOCITY OF SEEPAGE

      NMCHA = 2                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF FLOW DOF

C	TRANSFER NODAL VALUE INTO FLOW PART
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
	HMP(I) = HM(I)
	ENDDO

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


C	========================================
	IF(KFRES.EQ.1.AND.KSTEP.GT.1) THEN  !KFRES = FREE SURFACE ITERATION OPTION FOR 1
	GLEV = 0.0D0
	PTEN = 0.0D0
	DO INM = 1,NNM
	GLEV = GLEV + H(INM)*COORD(NGRAV,INM)
	PTEN = PTEN + H(INM)*EDIS(3*INM)
	ENDDO
	IF(GLEV.GT.PTEN) HM(1:2) = HMP(1:2)*1.0E-6
	ENDIF

	CALL TRNMAT2T(COORD,P,NNM,TMAT)
	DO ICOM = 1,2
	VXY(ICOM) = 0.0D0
	DO INM = 1,NNM
	VXY(ICOM) = VXY(ICOM) - HM(ICOM)*BH(ICOM,INM)*EDIS(3*INM)
	ENDDO
	ENDDO
	DO IRS = 1,2
	WA(IRS,IPT) = 0.0D0
	DO JRS = 1,2
	WA(IRS,IPT) = WA(IRS,IPT) + TMAT(IRS,JRS)*VXY(JRS)
	ENDDO
	ENDDO
C	========================================

	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,2
	SH(I,J) = SH(I,J) + HM(K)*BH(K,I)*BH(K,J)*DVOL
	ENDDO
	ENDDO
	ENDDO


 900  CONTINUE
 	
C	============================================================
C	TRANSFER LOCAL STIFFNESS TO GLOBAL ELEMENT STIFFNESS
C	============================================================

	DO I = 1,NNM             !FOR FLOW DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	ENDDO

C	==============================


	RETURN
	END

C	=======================================================================
C	==============START OF 4 NODE TETRAHEDRA SEEP ELEMENT==================
C	=======================================================================
	SUBROUTINE SEEPELE3T(COORD,PROPM,EDIS,S,RE,WA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET

      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
 	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	COMMON /MGRAV/ NGRAV

	PARAMETER (NFACE=4)

	DIMENSION COORD(3,NNM),S((NEF*NEF+NEF)/2),WA(NWG,1)
	DIMENSION SH(NNM,NNM)
	DIMENSION HM(3),HMP(3)
	DIMENSION RE(NEF),RH(NNM)
	DIMENSION PROPM(1)
	DIMENSION BH(3,NNM),EDIS(NEF),ETEMP(NNM)
	DIMENSION H(NNM),P(3,NNM),XJ(3,3),XJI(3,3)
	DIMENSION H2(NNM),G(2,NNM),NN(4,6),COVR(3),COVS(3),VT(3)

      NMCHA = 3                !NUMBER OF MECHANICAL DOF   
	NMTEM = 1                !NUMBER OF FLOW DOF

C	TRANSFER NODAL VALUE INTO FLOW PART
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
	HMP(I) = HM(I)
	ENDDO

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

C	---------------------------------------
	IF(KFRES.EQ.1.AND.KSTEP.GT.1) THEN  !KFRES = FREE SURFACE ITERATION OPTION FOR 1
	GLEV = 0.0D0
	PTEN = 0.0D0
	DO INM = 1,NNM
	GLEV = GLEV + H(INM)*COORD(NGRAV,INM)
	PTEN = PTEN + H(INM)*EDIS(4*INM)
	ENDDO
	IF(GLEV.GT.PTEN) HM(1:3) = HMP(1:3)*1.0E-6
	ENDIF

	DO ICOM = 1,3
	WA(ICOM,IPT) = 0.0D0
	DO INM = 1,NNM
	WA(ICOM,IPT) = WA(ICOM,IPT) - HM(ICOM)*BH(ICOM,INM)*EDIS(4*INM)
	ENDDO
	ENDDO
C	---------------------------------------

	DO I = 1,NNM
	DO J = 1,NNM
	DO K = 1,3
	SH(I,J) = SH(I,J) + HM(K)*BH(K,I)*BH(K,J)*DVOL
	ENDDO
	ENDDO
	ENDDO


 900  CONTINUE
 
C	======================================================
C	TRANSFER LOCAL STIFFNESS TO GLOBAL ELEMENT STIFFNESS
C	======================================================

	DO I = 1,NNM             !FOR FLOW DOF
	II = (NMCHA+NMTEM)*I
	DO J = I,NNM
	JJ = (NMCHA+NMTEM)*J
	LL = II - 1
	KK = (LL*LL+LL)/2 + (NEF-LL)*LL + (JJ-II+1)
	S(KK)  =  S(KK) + SH(I,J)
	ENDDO
	ENDDO

C	==============================


	RETURN
	END

C	=======================================================================
C	=======================================================================
C	=======================================================================

