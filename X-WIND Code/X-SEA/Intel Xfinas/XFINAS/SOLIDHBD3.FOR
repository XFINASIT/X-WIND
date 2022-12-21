      SUBROUTINE CONSOLEAS3D(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,
	1                       MWG,ALPHA,SEL,SEDI,FIN,SEQ,RH,RW)
C
C	----------------------------------------------------------------
C	THIS SUBROUTINE USED TO COMPUTE STIFFNSS MATRIX OF CONSOLIDATION
C	BASED ON ENGINEERING STRAIN NLOPT.EQ.1 AND ENHANCED ASSUMED STRAIN
C	APPLICATION TO MATERIAL MODEL
C	-----------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
C
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
C
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
C#####ADDED 31-08-08
	COMMON /MMENH/ MM,MM1,MM2,NDIMC

C
C#####ADDED 09-09-08 ADDED PARAMETER TIMEF FOR CONSOLIDATION
	COMMON /TFAC/ DFTIM,TIME0,TIMEF
C
      DIMENSION PROPM(NMP),PROPG(NGP),NODEX(NEX),WA(MWG,1),SS(300)
	1         ,COORD(3,NNO)
C
      DIMENSION EDIS(32),EDISI(32),RE(32),FIN(32)
C
      DIMENSION S(528),ST(32,32)
C
      DIMENSION DP(6,6),H(21),P(3,21),XJI(9),B(63),DISD(9)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6)
C
      DIMENSION RES(MM1)
C
C     --------------------------------
C     ENHANCED ASSUMED STRAIN VARIABLE
C	--------------------------------
C	SED (24,24) =  ENHANCED STIFFNESS MATRIX
C	SEDI(24,24) =  INVERSE OF ENHANCED STIFFNESS MATRIX
C	SEL (24,24) =  COUPLE COMPATIBLE-ENHANCED STIFFNESS MATRIX
C     ----------------------------------------------------------

	DIMENSION BM(6,24),CE(6,6),SK(24,24)
C
	DIMENSION GE(6,24)
	DIMENSION T0(6,6),TT(6,6),TT0(6,6)   
	DIMENSION XJ(3,3),XJ0(3,3)
	DIMENSION EAS(6),ALPHA(24)
      DIMENSION SED(24,24),SEL(24,24),SEDI(24,24)
	DIMENSION SE(24,24),SG(300)
C
	DIMENSION RH(MM),REH(MM)
C
C	-----------------------------------------------------------
C	ADDITIONAL MATRIX FOR ENHANCED ASSUMED STRAIN CONSOLIDATION
C	-----------------------------------------------------------

	DIMENSION SEQ(MM,MM2) ! COUPLE ENHANCED STRAIN-PRESSURE
	DIMENSION SKQ(MM1,MM2) ! COUPLE SOLID-PRESSURE
	DIMENSION SP(MM2,MM2)  ! PRESSURE-PRESSURE
C
	DIMENSION PM(3,3),HP(1,MM2) !HP PRESSURE SHAPE FUNCTION    
	DIMENSION EM(6,1),QM(3,MM2) 
C
	DIMENSION EDISU(MM1),EDISP(MM2)
C
	DIMENSION PRESS(6)
C
	DIMENSION RW(MM),RP(MM1),RS(MM1)
	DIMENSION FAU(MM1),FAP(MM2)
	DIMENSION RPR(MM2)
C
C	STARTING SUBROUTINE
C	
	EDISU(1:MM1) = 0.0D0
	EDISP(1:MM2) = 0.0D0
C
	CALL EXTRACTDISP38(EDIS,EDISU)  !ELEMENT DISPLCEMENT
	CALL EXTRACTPRES38(EDIS,EDISP)  !PORE PRESSURE
C
C	-------------------
C	INITIALIZED VAIABLE
C	------------------- 

	HP(1,1:MM2) = 0.0D0

	RH(1:MM)= 0.0D0
	RW(1:MM)= 0.0D0

C	-------------------
C	SUB. ELEMENT GROUPS
C	-------------------
C	SK = SOLID STIFFNESS MATRIX
C	SP = PRESSURE STIFFNESS MATRIX
C	SS = UPPER-TRIANGULAR SOLID STIFFNESS
C	SG = UPPER-TRIANGULAR SOLID GEOMETRIC STIFFNESS
C	-----------------------------------------------
C
	SK (1:MM1,1:MM1) = 0.0D0
	SKQ(1:MM1,1:MM2) = 0.0D0
	SP (1:MM2,1:MM2) = 0.0D0
	SS(1:300)= 0.0D0
	SG(1:300)= 0.0D0
C
C	------------------------------
C	SUB. ENHANCED STIFFNESS MATRIX
C	------------------------------
C	SE   = ENHANCED SOLID STIFFNESS MATRIX  
C	SED  = ENHANCED STIFFNESS
C	SEDI = INVERSE OF [SED]
C	SEL  = COUPLE ENHANCED-SOLID STIFFNESS MATX.
C	SEQ  = COUPLE ENHANCED-PRES. STIFFNESS MATX.
C
	SE  (1:MM1,1:MM1)  = 0.0D0
	SED (1:MM,1:MM)  = 0.0D0
	SEDI(1:MM,1:MM)  = 0.0D0
	SEL (1:MM,1:MM1)  = 0.0D0
	SEQ (1:MM,1:MM2)  = 0.0D0
C
C	-----------------
C	ELEMENT STIFFNESS
C	-----------------
C	S   = ELEMENT STIFFNESS STORE IN UPPER-TRIANGULAR
C	ST  = ELEMENT STIFFNESS STORE IN MATRIX FORM
C	--------------------------------------------
C
	S(1:528)	   = 0.0D0
	ST(1:32,1:32)  = 0.0D0
C
	RES(1:24) = 0.0D0
	RP (1:24) = 0.0D0
	RS(1:MM1)  = 0.0D0
C
	RE(1:32)  = 0.0D0
C
	TAU(1:6)	 = 0.0D0
	EAS(1:6)	 = 0.0D0
	GE(1:6,1:MM)	 = 0.0D0
	B(1:24) = 0.0D0
C
C*****ADDITIONAL FOR BIOT
C
	PM(1:3,1:3)	= 0.0D0
	PRESS(1:6)  = 0.0D0
	FAU(1:MM1)  = 0.0D0 
	FAP(1:MM2)	= 0.0D0
	RES(1:MM1)	= 0.0D0
	RPR(1:MM2)  = 0.0D0
C
C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------
C
      CALL HOKLAW (PROPM,PROPG,1)
      CALL PERMEA (PROPM,PM)  !PM PERMEABILITY MATRIX (CONSTANT PROPERTIES)
C
C	--------------------------------------
C	DEFINES NORMAL STRESS COMPONENT MATRIX
C	--------------------------------------
C
	EM(1:6,1) = 0.0D0
	EM(1:3,1) = 1.0D0 !EM NORMAL COMPONENTS
C
C	------------------
C	INTERGATION POINTS
C	------------------
C
      MGR = NGR
      MGS = NGS
      MGT = NGT
C
	SED = 0.0D0
C
C     ------------------------
C     COMPUTE TRANSFORM MATRIX
C     ------------------------
C
      CALL SHAP3D (0.0D0,0.0D0,0.0D0,H,P,0,8)   ! RETURN H,P
      CALL JACO3D (COORD,P,XJ0,XJI,DET0,MEL,8)
C
      CALL MATRIXT(XJ0,T0)
      CALL INVMATRIX(T0,TT,6)

      TT0 = TRANSPOSE(TT)
C
C	------------------------------
C	COMPLETE TRANSFORMATION MATRIX
C	------------------------------
C
      IF (ITASK.NE.5) GOTO 10
C
      MGR = 3
      MGS = 3
      MGT = 3
C
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------

 10   IPT = 0
      
      DO 900  IGR=1,MGR
      RI = GLOC(IGR,MGR)
      DO 900  IGS=1,MGS
      SI = GLOC(IGS,MGS)
      DO 900  IGT=1,MGT
      TI = GLOC(IGT,MGT)
      WT = GWT(IGR,MGR)*GWT(IGS,MGS)*GWT(IGT,MGT)
      IPT = IPT+1

C     ---------------------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ---------------------------------------------------
      
      CALL SHAP3D (RI,SI,TI,H,P,NODEX,NNO)
      CALL JACO3D (COORD,P,XJ,XJI,DET,MEL,NNO)
      DVOL = WT*DET
	
C
C	--------------------------------------
C	TRANSFER SHAPE FUNCTION TO MATRIX FORM
C	--------------------------------------
C
	HP(1,1:MM2) = H(1:8)
C
C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
      IF (ITASK.NE.5)  GOTO 50
      CALL SOMASS (SS,H,PROPM(5),DVOL,NNO,MM1,IPT)
      GOTO 900
C     ----------------------------------------
C     COMPACTED STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------
 50   CONTINUE
C
	CALL SOBMAT (P,XJI,B,NNO)
 	CALL BMSLD  (P,XJI,BM,8)
	CALL CSOLID (CE)
	CALL GRADP  (BM,QM)
C
C     --------------------------     
C     FOR ENHANCED STRAIN METHOD 
C     COMPUTE MATRIX [GE]
C     --------------------------
C
	CALL SOLIDGE(RI,SI,TI,DET0,DET,TT0,24,GE)
C
C     ------------------------
C	ENHANCED STIFFNESS SOLID
C     ------------------------
C
	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(CE,GE))*DVOL
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(CE,BM))*DVOL
C
C	-----------------------
C	COUPLING SOLID-PRESSURE
C	-----------------------
C
	SEQ = SEQ +  MATMUL(TRANSPOSE(GE),MATMUL(EM,HP))*DVOL
	SKQ = SKQ +  MATMUL(TRANSPOSE(BM),MATMUL(EM,HP))*DVOL
C
C	--------
C	PRESSURE
C	--------
C
	SP  = SP  +  MATMUL(TRANSPOSE(QM),MATMUL(PM,QM))*DVOL
C
C     NEXT LINE SUPPRESSED FOR NONLINEAR BASED ALGORITHM
C     IF (NLOPT+ITASK.EQ.1)  GOTO 700
C
C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------
C
      IF (ITASK.NE.4)  GOTO 200
C
      DO 100  I=1,6
 100  TAU(I) = WA(I,IPT)*DVOL
      GOTO 800
C
C     -------------------------------
C     DISPLACEMENT DERIVATIVES (DISD)
C     -------------------------------
C
200   CONTINUE
C
	CALL CLEARA (DISD,9)
C
      DO 290  IEF=1,MM1,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDISU(IEF)
      DISD(2) = DISD(2) + B(JEF)*EDISU(JEF)
      DISD(3) = DISD(3) + B(KEF)*EDISU(KEF)
      DISD(4) = DISD(4) + B(JEF)*EDISU(IEF)
      DISD(5) = DISD(5) + B(KEF)*EDISU(IEF)
      DISD(6) = DISD(6) + B(IEF)*EDISU(JEF)
      DISD(7) = DISD(7) + B(KEF)*EDISU(JEF)
      DISD(8) = DISD(8) + B(IEF)*EDISU(KEF)
290   DISD(9) = DISD(9) + B(JEF)*EDISU(KEF)
C
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
C      
      STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)
      STRAIN(4) = DISD(4) + DISD(6)
      STRAIN(5) = DISD(5) + DISD(8)
      STRAIN(6) = DISD(7) + DISD(9)
C
C     -------------------------
C     FOR ENHANCED STRAIN TEMRS
C     -------------------------
C
	EAS = MATMUL(GE,ALPHA)
C
C	--------------------------------------
C	COMPUTE TOTAL LINEAR COMPATIBLE STRAIN
C	-------------------------------------- 
C
	DO I=1,6
       STRAIN(I) = STRAIN(I) - EAS(I)
	END DO

C
C     ----------------------------------------------------------------
C     FOR NLOPT.GT.2 SUBTRACT QUADRATIC STRAIN TERMS (ALMANSI STRAINS)
C	UPDATED LAGRAGIAN CONFIGURATION
C     -------------------------------
C
C    IF (NLOPT.LE.1) GOTO 400
C
      QSTRAI(1) = .5*(DISD(1)*DISD(1)+DISD(6)*DISD(6)+DISD(8)*DISD(8))
      QSTRAI(2) = .5*(DISD(4)*DISD(4)+DISD(2)*DISD(2)+DISD(9)*DISD(9))
      QSTRAI(3) = .5*(DISD(5)*DISD(5)+DISD(7)*DISD(7)+DISD(3)*DISD(3))
      QSTRAI(4) =     DISD(1)*DISD(4)+DISD(6)*DISD(2)+DISD(8)*DISD(9)
      QSTRAI(5) =     DISD(1)*DISD(5)+DISD(6)*DISD(7)+DISD(8)*DISD(3)
      QSTRAI(6) =     DISD(4)*DISD(5)+DISD(2)*DISD(7)+DISD(9)*DISD(3)
C
      DO 390  I=1,6
C 390  STRAIN(I) = STRAIN(I)-QSTRAI(I)
C#####ABOVE LINE SUPPRESSED BECAUSE OF ENGNEERING STRAIN ASSUMPTION
390   STRAIN(I)  = STRAIN(I)
C
C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------
400   IF (MTMOD-3)  405,450,460
405   CONTINUE
C
C	------------------------------
C	ELASTIC ISOTROPIC CONSTITUTIVE
C	------------------------------
C
	PRESS(1:6)  = 0.0D0
	STRESS(1:6) = 0.0D0
C
	PRESS  = MATMUL(EM,MATMUL(HP,EDISP))  ! CURRENT CAUCHY PRESSURE
	STRESS = MATMUL(CE,STRAIN)			  ! CURRENT CAUCHY STRESS

      DO 410  I=1,6
 410  WA(I,IPT) = STRESS(I)  !EFFECTIVE STRESS

      GOTO 500
C
 450	CONTINUE
C
C	-------------------------------------
C	ELASTO-PLASTIC VON-MISES CONSTITUTIVE
C	-------------------------------------
C
	CALL MISE3D (WA(1,IPT),WA(7,IPT),WA(13,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DP)
      GOTO 500
C
 460  CONTINUE
	CALL MOHR3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			 STRAIN,DP,PROPM)

	GOTO 500

 500  CONTINUE
C
C	------------------------------------------------
C	GEOMETRIC STIFFNESS REQUIRES TOTAL CAUCHY STRESS
C	------------------------------------------------
C
	DO I=1,6
	TAU(I) = (STRESS(I)+PRESS(I))*DVOL
	END DO
C
      IF (ITASK.LE.2) GOTO 520
C
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900
C
 520  CONTINUE
C
C	-------------------------------------------
C	-COMPATIBLE EQUILIBRIUM FORCE
C	-ENHANCED ASSUMED STRAIN  EQUILIBRIUM FORCE
C	-ENHANCED PRESSURE EQUILIBRIUM
C	-------------------------------------------
C
	RS  = RS + MATMUL(TRANSPOSE(BM),STRESS)*DVOL
	RP  = RP + MATMUL(TRANSPOSE(BM),PRESS )*DVOL
C
	RH = RH + MATMUL(TRANSPOSE(GE) ,STRESS)*DVOL
	RW = RW + MATMUL(TRANSPOSE(GE) ,PRESS )*DVOL
C
C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------
C
      IF (IFREF) 900,700,900
C
 700  IF (NLOPT.NE.0) GOTO 750
C
C	------------------------------------------------
C	LINEAR ELASTICS NLOPT = 0, REFORMATION STIFFNESS
C	------------------------------------------------
C
	SK = (MATMUL(TRANSPOSE(BM),MATMUL(CE,BM)))*DVOL
	KSK = 0
	DO ISK =1,24
	DO JSK =ISK,24
	KSK = KSK+1 
	SS(KSK)=SS(KSK)+SK(ISK,JSK)
	END DO
	END DO
C
      GOTO 790
C
 750  CONTINUE

	IF (MTMOD.LE.2)THEN
C
C	----------------------------------------------------------------------
C	NONLINEAR ELASTICS NLOPT = 3 OR NLOPT = 1 ELASTICS MATERTIAL, REFORMATION STIFFNESS
C	----------------------------------------------------------------------
C
	SK = (MATMUL(TRANSPOSE(BM),MATMUL(CE,BM)))*DVOL
	KSK = 0
	DO ISK =1,24
	DO JSK =ISK,24
	KSK = KSK+1 
	SS(KSK)=SS(KSK)+SK(ISK,JSK)
	END DO
	END DO
C
      END IF
C
	IF (MTMOD.GT.2) THEN
C
C	-----------------------------------------------------
C	ELASTO-PLASTIC CONSTITUTIVE MATERIALS SOLID STIFFNESS 
C	-----------------------------------------------------
C
	SK = (MATMUL(TRANSPOSE(BM),MATMUL(DP,BM)))*DVOL
	KSK = 0
	DO ISK =1,24
	DO JSK =ISK,24
	KSK = KSK+1 
	SS(KSK)=SS(KSK)+SK(ISK,JSK)
	END DO
	END DO
C
	END IF
C
C     -------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------
C
 790  IF (NLOPT.LE.1) GOTO 810
C
 800  CONTINUE
C
C     -------------------------------------------------
C     COMPUTE GEOMETRICAL NONLINEAR STIFFNESS NLOPT = 3
C     -------------------------------------------------
C
	CALL KNLSTIFF(SG,TAU,B)
C
 810  CONTINUE
 900  CONTINUE
C
      IF (ITASK.EQ.5) GOTO 950

	CALL INVMATRIX(SED,SEDI,24)
C
C     -----------------------------------------------------
C	COMPLETE COMPATIBLE-ENHANCED ELEMENT STIFFNESS MATRIX 
C	-----------------------------------------------------
C
	SE =  MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,SEL))
C
C---------------------------------------------------------
C	-EQUILIBRIUM FORCE     ENHANCED ASSUMED STRAIN SOLID
C	-EQUILIBRIUM PRESSURE  ENHANCED ASSUMED STRAIN PRESURE
C-----------------------------------------------------------

	FAU = MATMUL(MATMUL(TRANSPOSE(SEL),SEDI),RH+RW)
	FAP = MATMUL(MATMUL(TRANSPOSE(SEQ),SEDI),RH+RW)
C
C	----------------------------------------
C	COMPLETE CORRECT EQUILIBRIUM FORCE SOLID
C	----------------------------------------
C
	DO ISD=1,MM1
	RES(ISD) = RS(ISD) + RP(ISD)-FAU(ISD)
	END DO
C
C	----------------------------------------------
C	PRESSURE EQUILIBRIUM FOR NEXT INCREMENTAL STEP
C	CURRENT CAUCHY PRESSURE 
C	-----------------------
C
	RPR = MATMUL(SP,EDISP)
C
	DO IPR=1,MM2
	RPR(IPR) = -1.0D0*DFTIM*RPR(IPR)-FAP(IPR)
	END DO
C
C	---------------------------------------------------------
C	COMPUTE LINEAR ELEMENT STIFFNESS MATRIX UPPER-TRIAGNGULAR
C	---------------------------------------------------------
C


	KSK = 0
	DO ISK =1,24
	DO JSK = ISK,24
	KSK = KSK+1 
	SS(KSK)=SS(KSK)-SE(ISK,JSK)
	END DO
	END DO
C
C	-----------------------------------------
C	COUPLING SOLID-PRESSURE STINFFNESS MATRIX
C	-----------------------------------------
C
	SKQ = SKQ - MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,SEQ))
C
C	-------------------------
C	PRESSURE STIFFNESS MATRIX
C	------------------------- 
C
	SP  = -1.0D0*SP*DFTIM - MATMUL(TRANSPOSE(SEQ),MATMUL(SEDI,SEQ))
C
C	-------------------------------------------------------
C	ADDITIONAL GEOMETRIC STIFFNESS MATRIX UPPER-TRIAGNGULAR
C	-------------------------------------------------------
C
	KSK = 0
	DO ISK =1,24
	DO JSK = ISK,24
	KSK = KSK+1 
C	SS(KSK)=SS(KSK) + SG(KSK)
C#####ABOVE LINE SUPPRESSED BECAUSE OF ENGNEERING STRAIN ASSUMPTION
      SS(KSK)=SS(KSK)
	END DO
	END DO
C
C	----------------------------
C	FILL MATRIX UPPER TRIANGULAR
C	----------------------------
C
950	SK(1:MM1,1:MM1) = 0.0D0
	KSK = 0
	DO ISK =1,MM1
	DO JSK = ISK,MM1
	KSK = KSK+1 
	SK(ISK,JSK)=SS(KSK)
	END DO
	END DO
C
	ST(1:32,1:32) = 0.0D0
C
C	-------------------------
C	FILL COMPLETE MATRIX FORM 
C	-------------------------
C
	DO ISK=1,MM1
	DO JSK=1,MM1
	SK(JSK,ISK) = SK(ISK,JSK)
	END DO
	END DO
C
C	-------------------------------------------------------------
C	ASSEMBLE SK,SKQ,TRANS(SKQ),SP TO TOTAL ELEMENT STIFFNESS [ST]
C	-------------------------------------------------------------
C
	DO I =1,24
      DO J =1,24
	ST(I,J) = SK(I,J)
	END DO
	END DO
C
C	-------------------------------------------------------------------
C	FILL SUB. MATRIX COUPLE-SOLID STIFFNESS TO ELEMENT STIFFNESS MATRIX
C	-------------------------------------------------------------------
C
	DO I =1,24
	DO J =1,8
	ST(I,J+24) = SKQ(I,J)  !SKQ
	ST(J+24,I) = SKQ(I,J)  !TRANSPOSE(SKQ) 
	END DO
	END DO
C
C	---------------------------------------------------------------
C	FILL SUB. MATRIX PRESSURE STIFFNESS TO ELEMENT STIFFNESS MATRIX
C	---------------------------------------------------------------
C
	DO I=1,8
	DO J=1,8
	ST(I+24,J+24) = SP(I,J)
	END DO
	END DO
C
C	------------------------------------------------------------------
C	ARRANGE TO UPPER-ROW WISE MATRIX S(528) COMPLETE ELEMENT STIFFNESS
C	------------------------------------------------------------------
C
	CALL UPPERTRIA(S,ST)			!ELEMENT STIFFNESS MATRIX
C
C	-------------------------------
C	REARRANGE THE EQUILIBRIUM FORCE
C	-------------------------------
C
	CALL EQUIFORCE38(RES,RPR,RE)    !ELEMENT EQUILIBRIUM FORCE

	FIN(1:32) = RE(1:32)


c	KSK = 0
c	DO ISK =1,32
c	DO JSK =ISK,32
c	KSK = KSK+1 
c	IF(ISK.EQ.JSK) WRITE(*,*) ISK,S(KSK)
c	WRITE(110,*) ISK,JSK,S(KSK)
c	END DO
c	END DO
c
c	PAUSE

C
      RETURN
      END
