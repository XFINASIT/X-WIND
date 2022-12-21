C	----------------------------------------------------------------
C	----------------------------------------------------------------
C	----------------------------------------------------------------
      SUBROUTINE CONSOLEAS2D(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,
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

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
C
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
C
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

C	EAS
	COMMON /MMENH/ MM,MM1,MM2,NDIMC

C	FOR CONSOLIDATION
	COMMON /TFAC/ DFTIM,TIME0,TIMEF
C
      DIMENSION PROPM(1),PROPG(1),NODEX(NEX),WA(MWG,1),
	1          COORD(NCO,NNO)
C
      DIMENSION EDIS(1),EDISI(1),RE(1),FIN(1)
C
      DIMENSION S(1),ST(NEF,NEF)
C
      DIMENSION DMM(6,6),DE4(6,6),CE(4,4),H(8),P(2,8)
	DIMENSION B(4,16),DISD(5)
      DIMENSION STRAIN(4),QSTRAI(4),STRESS(4),TAU(4)
	DIMENSION XJ(2,2),XJ0(2,2),XJI(2,2),KK(4)

C     --------------------------------
C     ENHANCED ASSUMED STRAIN VARIABLE
C	--------------------------------
C	SED (24,24) =  ENHANCED STIFFNESS MATRIX
C	SEDI(24,24) =  INVERSE OF ENHANCED STIFFNESS MATRIX
C	SEL (24,24) =  COUPLE COMPATIBLE-ENHANCED STIFFNESS MATRIX
C     ----------------------------------------------------------

	DIMENSION BM(4,MM1),GE(4,MM)
	DIMENSION T0(3,3),TT(3,3),TT0(3,3)   
	DIMENSION EAS(4),ALPHA(MM)
      DIMENSION SED(MM,MM),SEL(MM,MM1),SEDI(MM,MM)
	DIMENSION SE(MM1,MM1),SK(MM1,MM1),SS((MM1*MM1+MM1)/2)
C
	DIMENSION RH(MM),REH(MM)
C
C	-----------------------------------------------------------
C	ADDITIONAL MATRIX FOR ENHANCED ASSUMED STRAIN CONSOLIDATION
C	-----------------------------------------------------------

	DIMENSION SEQ(MM ,MM2)  ! COUPLE ENHANCED STRAIN-PRESSURE
	DIMENSION SKQ(MM1,MM2)  ! COUPLE SOLID-PRESSURE
	DIMENSION  SP(MM2,MM2)  ! PRESSURE-PRESSURE
C
	DIMENSION PM(2,2),HP(1,MM2) !HP PRESSURE SHAPE FUNCTION    
	DIMENSION EM(4,1),QM(2,MM2) 
C
	DIMENSION EDISU(MM1),EDISP(MM2)
C
	DIMENSION PRESS(4)
C
	DIMENSION RW(MM),RP(MM1),RS(MM1)
	DIMENSION FAU(MM1),FAP(MM2)
	DIMENSION RES(MM1),RPR(MM2)
C
C	STARTING SUBROUTINE
C
	TH = PROPG(2)

	ISTYP1 = ISTYP
C	ISTYP FOR PLANE STRESS - PLANE STRAIN FLAG
	ISTYP = ISTYP - 3 !SONGSAK SEP2005
	
	
	EDISU(1:MM1) = 0.0D0
	EDISP(1:MM2) = 0.0D0

	CALL EXTRACTDISP (EDIS,EDISU)  !ELEMENT DISPLCEMENT
	CALL EXTRACTPRES (EDIS,EDISP)  !PORE PRESSURE

C	-------------------
C	INITIALIZED VAIABLE
C	------------------- 
	HP(1,1:MM2) = 0.0D0

	RH(1:MM) = 0.0D0
	RW(1:MM) = 0.0D0
C	-------------------
C	SUB. ELEMENT GROUPS
C	-------------------
C	SK = SOLID STIFFNESS MATRIX
C	SP = PRESSURE STIFFNESS MATRIX
C	SS = UPPER-TRIANGULAR SOLID STIFFNESS
C	-----------------------------------------------
	SK (1:MM1,1:MM1) = 0.0D0
	SKQ(1:MM1,1:MM2) = 0.0D0
	SP (1:MM2,1:MM2) = 0.0D0
	SS(1:(MM1*MM1+MM1)/2)= 0.0D0
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
	SE  (1:MM1,1:MM1) = 0.0D0
	SED (1:MM,1:MM)   = 0.0D0
	SEDI(1:MM,1:MM)   = 0.0D0
	SEL (1:MM,1:MM1)  = 0.0D0
	SEQ (1:MM,1:MM2)  = 0.0D0

C	-----------------
C	ELEMENT STIFFNESS
C	-----------------
C	S   = ELEMENT STIFFNESS STORE IN UPPER-TRIANGULAR
C	ST  = ELEMENT STIFFNESS STORE IN MATRIX FORM
C	--------------------------------------------
	 S(1:(NEF*NEF+NEF)/2) = 0.0D0
	ST(1:NEF,1:NEF)       = 0.0D0

	RP(1:MM1) = 0.0D0
	RS(1:MM1) = 0.0D0

	RE(1:NEF)  = 0.0D0

	TAU(1:4)	 = 0.0D0
	EAS(1:4)	 = 0.0D0
	GE(1:4,1:MM) = 0.0D0
	 B(1:4,1:16) = 0.0D0
	BM(1:4,1:MM1)= 0.0D0

C	ADDITIONAL FOR BIOT
	PM(1:2,1:2)	= 0.0D0
	PRESS(1:4)  = 0.0D0
	FAU(1:MM1)  = 0.0D0 
	FAP(1:MM2)	= 0.0D0
	RES(1:MM1)	= 0.0D0
	RPR(1:MM2)  = 0.0D0

C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------
      CALL HOKLAW   (PROPM,PROPG,1)
      CALL PERMEA2D (PROPM,PM)  !PM PERMEABILITY MATRIX (CONSTANT PROPERTIES)
	CALL DELA3D   (DMM,ISTYP)

	KK(1:4) = [1,2,3,6]
	DO I = 1,4
	II = KK(I)
	DO J = 1,4
	JJ = KK(J)
	CE(I,J) = DMM(II,JJ)
	ENDDO
	ENDDO


C	FOR AXISSYMETRIC
	DD1 = COORD(1,1)
	DD2 = COORD(1,2)
	DD3 = COORD(1,3)
	DD4 = COORD(1,4)
	RTA0 = 0.25*( DD1 + DD2 + DD3 + DD4 )  ! A0 = 1  1  1  1
	RTA1 = 0.25*( DD1 - DD2 - DD3 + DD4 )  ! A1 = 1 -1 -1  1
	RTA2 = 0.25*( DD1 + DD2 - DD3 - DD4 )  ! A2 = 1  1 -1 -1
	ETA1 = RTA1/RTA0/3.0
	ETA2 = RTA2/RTA0/3.0

C	--------------------------------------
C	DEFINES NORMAL STRESS COMPONENT MATRIX
C	--------------------------------------
	EM(1:4,1) = 0.0D0
	EM(1:2,1) = 1.0D0 !EM NORMAL COMPONENTS

C	------------------
C	INTERGATION POINTS
C	------------------
      MGR = NGR
      MGS = NGS

C     ------------------------
C     COMPUTE TRANSFORM MATRIX
C     ------------------------
C
      CALL SHAP2D (0.0D0,0.0D0,H,P,NODEX,NNO)			! RETURN H,P
	CALL JACOB2D(MEL,NNO,COORD,P,XJ0,XJI,DET0,NCO)	
C
	CALL TTOMAT(T0,XJ0)
      CALL INVMATRIX(T0,TT,3)

      TT0 = TRANSPOSE(TT)

C	------------------------------
C	COMPLETE TRANSFORMATION MATRIX
C	------------------------------
      IF (ITASK.NE.5) GOTO 10

      MGR = 3
      MGS = 3

C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------

 10   IPT = 0
      
      DO 900  IGR=1,MGR
      RI = GLOC(IGR,MGR)
      DO 900  IGS=1,MGS
      SI = GLOC(IGS,MGS)
      WT = GWT(IGR,MGR)*GWT(IGS,MGS)
      IPT = IPT+1

C     ---------------------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ---------------------------------------------------
      
      CALL SHAP2D (RI,SI,H,P,NODEX,NNO)
	CALL JACOB2D(MEL,NNO,COORD,P,XJ,XJI,DET,NCO)
	CALL MEBMAT (COORD,H,P,XJI,B,ISTYP,XBAR,NNO)
	
	DO I = 1,4
	BM(I,1:MM1) = B(I,1:MM1)
	ENDDO
		
      IF (ISTYP.NE.0) XBAR = TH
      DVOL = WT*DET*XBAR

C	--------------------------------------
C	TRANSFER SHAPE FUNCTION TO MATRIX FORM
C	--------------------------------------
	HP(1,1:MM2) = H(1:MM2)

C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
      IF (ITASK.NE.5)  GOTO 50
      CALL MASS2D (SS,H,PROPM(5),DVOL,NNO,MM1,IPT)
      GOTO 900
 50   CONTINUE

C     ----------------------------------------
C     COMPACTED PRESSURE - DISPLACEMENT MATRIX (QM)
C     ----------------------------------------
	CALL GRADP2D (BM,QM,NNO)

C     ---------------------------------------------     
C     FOR ENHANCED STRAIN METHOD 
C     COMPUTE MATRIX [GE] AND TRANSPOSE [GE] MATRIX
C     ---------------------------------------------
	R12 = H(1)*DD1 + H(2)*DD2 + H(3)*DD3 + H(4)*DD4
C	PLANE
	IF(ISTYP.NE.0) CALL GEMATRIX1(GE,RI,SI,TT0,MM,DET0,DET,MEL)	
C	AXISYMMETRIC
	IF(ISTYP.EQ.0) CALL GEMATRIX2(GE,RI,SI,TT0,MM,DET0,DET,MEL,
	1				              R12,ETA1,ETA2)

C     ------------------------
C	ENHANCED STIFFNESS SOLID
C     ------------------------
C
	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(CE,GE))*DVOL
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(CE,BM))*DVOL

C	-----------------------
C	COUPLING SOLID-PRESSURE
C	-----------------------
	SEQ = SEQ +  MATMUL(TRANSPOSE(GE),MATMUL(EM,HP))*DVOL
	SKQ = SKQ +  MATMUL(TRANSPOSE(BM),MATMUL(EM,HP))*DVOL

C	--------
C	PRESSURE
C	--------
	SP  = SP  +  MATMUL(TRANSPOSE(QM),MATMUL(PM,QM))*DVOL


	IF (NLOPT+ITASK.EQ.1)  GOTO 700
C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------
      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,4
 100  TAU(I) = WA(I,IPT)*DVOL
      GOTO 800

C     -------------------------------
C     DISPLACEMENT DERIVATIVES (DISD)
C     -------------------------------
 200  CALL CLEARA (DISD,5)
      DO 250  J=1,MM1,2
      J1 = J+1
      DISD(1) = DISD(1) + B(1,J) *EDISU(J)
      DISD(2) = DISD(2) + B(2,J1)*EDISU(J1)
      DISD(3) = DISD(3) + B(3,J) *EDISU(J)
 250  DISD(4) = DISD(4) + B(3,J1)*EDISU(J1)
      IF (ISTYP.GE.1)  GOTO 300
      DO 290  IEF=1,MM1,2
 290  DISD(5) = DISD(5) + B(4,IEF)*EDISU(IEF)
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
 300  STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)+DISD(4)
      STRAIN(4) = 0.0D0
      IF (ISTYP.EQ.0)  STRAIN(4) = DISD(5)

C     -------------------------
C     FOR ENHANCED STRAIN TEMRS
C     -------------------------
	EAS = MATMUL(GE,ALPHA)

C	--------------------------------------
C	COMPUTE TOTAL LINEAR COMPATIBLE STRAIN
C	-------------------------------------- 
	DO I=1,4
	STRAIN(I) = STRAIN(I) - EAS(I)
	END DO

C     --------------------------------------------------------------
C     FOR NLOPT>1 SUBSTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     --------------------------------------------------------------
      IF (NLOPT.LE.1)  GOTO 400
      STRAIN(1) = STRAIN(1) - 0.5*(DISD(1)*DISD(1) + DISD(4)*DISD(4))
      STRAIN(2) = STRAIN(2) - 0.5*(DISD(2)*DISD(2) + DISD(3)*DISD(3))
      STRAIN(3) = STRAIN(3) -     (DISD(3)*DISD(1) + DISD(4)*DISD(2))
      IF (ISTYP.EQ.0)  STRAIN(4) = STRAIN(4) - 0.5*(DISD(5)*DISD(5))

C	------------------------------
C	PRESSURE AT GUASS POINT
C	------------------------------
400	PRESS(1:4)  = 0.0D0
	PRESS = MATMUL(EM,MATMUL(HP,EDISP))  ! CURRENT CAUCHY PRESSURE

	STRESS(1:4) = 0.0D0
C     ------------------------------------
C     COMPUTE AND STORE GAUSS POINT STRESSES
C     ------------------------------------
	GOTO (410,410,430,440,450,440,450),MTMOD

410	CALL MEMSIG(STRAIN,STRESS,ISTYP) 
	DO I=1,4
	WA(I,IPT) = STRESS(I)
	END DO
	GOTO 480 

430	CALL VMISES(WA(1,IPT),WA(5,IPT),WA(9,IPT),STRAIN,STRESS,DE4
	1			,ISTYP)

	GOTO 460

440	CONTINUE

	WA(1:2,IPT) = WA(1:2,IPT) - 100.0 
	WA(  4,IPT) = WA(  4,IPT) - 100.0 

	CALL MOHR2D(WA(1,IPT),WA(5,IPT),WA(9,IPT),STRESS,
	1			STRAIN,DE4,PROPM,WA(10,IPT),WA(11,IPT),ISTYP)

	WA(1:2,IPT) = WA(1:2,IPT) + 100.0  
	WA(  4,IPT) = WA(  4,IPT) + 100.0
	STRESS(1:2) = STRESS(1:2) + 100.0  
	STRESS(  4) = STRESS(  4) + 100.0

	GOTO 460

450   CONTINUE

	CALL DRUG2D(WA(1,IPT),WA(5,IPT),WA(9,IPT),STRESS,
	1			STRAIN,DE4,PROPM,WA(10,IPT),WA(11,IPT),ISTYP)

	GOTO 460


460   CONTINUE


	CALL REARDEP(DMM,DE4)
C     ----------------------------------------
C     ADJUST THICKNESS (NLOPT>1 .AND. ISTYP=2)
C     ----------------------------------------
480	IF (ITASK.EQ.3)  GOTO 500

C	FOR PLANE GOTO 500 DIRECTLY
      IF (NLOPT.EQ.1.OR.ISTYP.LT.2)  GOTO 500

C	FOR PLANE STRESS CASE
      THICK = TH
      STRAIN(4) = -(DMM(6,1)*STRAIN(1)+DMM(6,2)*STRAIN(2)
     1            + DMM(6,3)*STRAIN(3)) /DMM(6,6)

      EXT  = 1. - 2.*STRAIN(4)
      EXT  = 1. - SQRT(EXT)

      THICK = THICK*(1.-EXT)
      DVOL = WT*DET*THICK

C	------------------------------------------------
C	GEOMETRIC STIFFNESS REQUIRES TOTAL CAUCHY STRESS
C	------------------------------------------------
500	DO I=1,4
	TAU(I) = (STRESS(I)+PRESS(I))*DVOL
	END DO


      IF (ITASK.LE.3.AND.ISOLOP.NE.4)  GOTO 520
      IF (IFEIG.EQ.0)  GOTO 800
      GOTO 900
C	-------------------------------------------
C	-COMPATIBLE EQUILIBRIUM FORCE
C	-ENHANCED ASSUMED STRAIN  EQUILIBRIUM FORCE
C	-ENHANCED PRESSURE EQUILIBRIUM
C	-------------------------------------------
 520  RS = RS + MATMUL(TRANSPOSE(BM),STRESS)*DVOL
	RP = RP + MATMUL(TRANSPOSE(BM),PRESS )*DVOL

	RH = RH + MATMUL(TRANSPOSE(GE),STRESS)*DVOL
	RW = RW + MATMUL(TRANSPOSE(GE),PRESS )*DVOL

C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------
      IF (IFREF) 900,700,900

 700  CONTINUE
	CALL MEK0PL (SS,DMM,B,DVOL,MM1,ISTYP)


C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
	IF (NLOPT.LE.1) GOTO 810
 800  CALL MEKSIG (SS,TAU,B,MM1,ISTYP)
 810  CONTINUE 


 900  CONTINUE
C     --------------------------------------------------------
C     --------------------------------------------------------

      IF (ITASK.EQ.5) GOTO 950

C     --------------------------------------------------------
C	INVERSE OF EAS MATRIX
C     --------------------------------------------------------
	CALL INVMATRIX(SED,SEDI,MM)

C     -----------------------------------------------------
C	COMPLETE COMPATIBLE-ENHANCED ELEMENT STIFFNESS MATRIX 
C	-----------------------------------------------------
	SE =  MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,SEL))

C	---------------------------------------------------------
C	-EQUILIBRIUM FORCE     ENHANCED ASSUMED STRAIN SOLID
C	-EQUILIBRIUM PRESSURE  ENHANCED ASSUMED STRAIN PRESURE
C	-----------------------------------------------------------
	FAU = MATMUL(MATMUL(TRANSPOSE(SEL),SEDI),RH+RW)
	FAP = MATMUL(MATMUL(TRANSPOSE(SEQ),SEDI),RH+RW)

C	----------------------------------------
C	COMPLETE CORRECT EQUILIBRIUM FORCE SOLID
C	----------------------------------------
	DO ISD=1,MM1
	RES(ISD) = RS(ISD) + RP(ISD) - FAU(ISD)
	END DO

C	----------------------------------------------
C	PRESSURE EQUILIBRIUM FOR NEXT INCREMENTAL STEP
C	CURRENT CAUCHY PRESSURE 
C	-----------------------
	RPR = MATMUL(SP,EDISP)

	DO IPR=1,MM2
	RPR(IPR) = -1.0D0*DFTIM*RPR(IPR)-FAP(IPR)
	END DO

C	---------------------------------------------------------
C	COMPUTE LINEAR ELEMENT STIFFNESS MATRIX UPPER-TRIAGNGULAR
C	---------------------------------------------------------
	KSK = 0
	DO ISK = 1  ,MM1
	DO JSK = ISK,MM1
	KSK = KSK+1 
	SS(KSK) = SS(KSK) - SE(ISK,JSK)
	END DO
	END DO

C	-----------------------------------------
C	COUPLING SOLID-PRESSURE STINFFNESS MATRIX
C	-----------------------------------------
	SKQ = SKQ - MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,SEQ))

C	-------------------------
C	PRESSURE STIFFNESS MATRIX
C	------------------------- 
	SP  = -1.0D0*SP*DFTIM - MATMUL(TRANSPOSE(SEQ),MATMUL(SEDI,SEQ))

C	----------------------------
C	FILL MATRIX UPPER TRIANGULAR
C	----------------------------
950	SK(1:MM1,1:MM1) = 0.0D0
	KSK = 0
	DO ISK = 1  ,MM1
	DO JSK = ISK,MM1
	KSK = KSK+1 
	SK(ISK,JSK)=SS(KSK)
	END DO
	END DO

C	-------------------------
C	FILL COMPLETE MATRIX FORM 
C	-------------------------
	DO ISK=1,MM1
	DO JSK=1,MM1
	SK(JSK,ISK) = SK(ISK,JSK)
	END DO
	END DO

C	-------------------------------------------------------------
C	ASSEMBLE SK,SKQ,TRANS(SKQ),SP TO TOTAL ELEMENT STIFFNESS [ST]
C	-------------------------------------------------------------
	ST(1:NEF,1:NEF) = 0.0D0
	DO I =1,MM1
      DO J =1,MM1
	ST(I,J) = SK(I,J)
	END DO
	END DO

C	-------------------------------------------------------------------
C	FILL SUB. MATRIX COUPLE-SOLID STIFFNESS TO ELEMENT STIFFNESS MATRIX
C	-------------------------------------------------------------------
	DO I =1,MM1
	DO J =1,MM2
	ST(I,J+MM1) = SKQ(I,J)  !SKQ
	ST(J+MM1,I) = SKQ(I,J)  !TRANSPOSE(SKQ) 
	END DO
	END DO

C	---------------------------------------------------------------
C	FILL SUB. MATRIX PRESSURE STIFFNESS TO ELEMENT STIFFNESS MATRIX
C	---------------------------------------------------------------
	DO I=1,MM2
	DO J=1,MM2
	ST(I+MM1,J+MM1) = SP(I,J)
	END DO
	END DO

C	------------------------------------------------------------------
C	ARRANGE TO UPPER-ROW WISE MATRIX S(528) COMPLETE ELEMENT STIFFNESS
C	------------------------------------------------------------------
	CALL UPPERTRIA2D(S,ST,NNO,NEF)			!ELEMENT STIFFNESS MATRIX

C	KEF = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	ENDDO
C	ENDDO
C	KEE = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	KEE = KEE + 1
C	S(KEF) = S(KEE)
C	ENDDO
C	ENDDO


C	-------------------------------
C	REARRANGE THE EQUILIBRIUM FORCE
C	-------------------------------
	CALL EQUIFORCE2D(RES,RPR,RE,NNO)    !ELEMENT EQUILIBRIUM FORCE


	IF (ITASK.EQ.3) FIN(1:NEF) = RE(1:NEF)


C	RECALL ISTYP
	ISTYP = ISTYP1 


      RETURN
      END


C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE EQUIFORCE2D(RES,RPR,RE,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------
C	EQUILIBRIUM LOAD FROM SOLID-PRESSURE TO ELEMENT LOAD
C	----------------------------------------------------
	DIMENSION RES(1),RPR(1),RE(1)
	
	DO I=1,NNO
	RE(3*I-2) = RES(2*I-1)
	RE(3*I-1) = RES(2*I-0)
	RE(3*I-0) = RPR(I)
	END DO

	RETURN
	END


C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE GRADP2D(BM,QM,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION BM(4,1),QM(2,1)

	DO I=1,NNO
	QM(1,I) = BM(1,2*I-1) 
	QM(2,I) = BM(2,2*I-0) 
	END DO

	RETURN
	END


C	=================================================================
C	=================================================================
C	=================================================================

	SUBROUTINE PERMEA2D(PROPM,PM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	------------------------------------------------------------
C	SETTING BASIC PROPERTIES OF PERMEABILITY MATRIX PER UNIT WT.
C	------------------------------------------------------------
	DIMENSION PROPM(1),PM(2,2)


	PM(1:2,1:2) = 0.0D0

C	INSERTING VALUE	
	PM(1,1) = PROPM(6)	
	PM(2,2) = PROPM(7)		

	RETURN
	END


C	=================================================================
C	=================================================================
C	=================================================================