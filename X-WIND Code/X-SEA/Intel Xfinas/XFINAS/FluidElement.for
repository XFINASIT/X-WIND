C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
	SUBROUTINE FLUELE2D (PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,
     +		             ELOD,MWG,ALPHA,SEL,SEDI,FIN,HINFC)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     ------------------------------------
C     CALLS THE APPROPRIATE ELEMENT MODULE
C     ------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /MMENH/ MM,MM1,MM2,NDIMC

      DIMENSION PROPM(1),PROPG(1),NODEX(1),WA(1),S(1),COORD(1),EDIS(1)
      DIMENSION EDISI(1),ELOD(1),FIN(1)

	DIMENSION ALPHA(MM,1),SEL(MM,MM1),SEDI(MM,MM),HINFC(MM)


	NEFU = NEF - NNO


	CALL FLUEAS2D(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG
     1			 ,ALPHA,SEL,SEDI,FIN,MEL,HINFC,NEFU)


	RETURN

      END


C	=======================================================================
C	=======================================================================
C	=======================================================================
      SUBROUTINE FLUEAS2D(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,MWG
     1			       ,ALPHA,SEL,SEDI,FIN,IEL,RH,NEFU)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     ------------------------------------------------------------------
C     MAIN PROGRAM FOR THE 2-D PLANE PROBLEM
C     EVALUATES THE TANGENTIAL STIFFNESS MATRIX,STRAINS AND STRESSES
C	FOR 4-NODES EAS(ENHANCED ASSUMED STRAIN METHODS)
C	---------------
C     INPUT VARIABLES
C	---------------
C     PROPM(NMP)    = MATERIAL PROPERTIES (YM,PR,YLD,HP,DEN)
C     PROPG(NGP)    = GEOMETRIC PROPERTIES (NNO)
C     NODEX(NEX)    = LOCATIONS OF EXCESS NODES (MIDSIDE NODES)
C     WA(MWG,NPT)   = WORKING ARRAY (6 STRESSES + (6 STRAINS,YLD,IPEL))
C     COORD(2,NNO)  = CURRENT NODAL COORDINATES X,Y,Z
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     EDISI(NEF)    = CURRENT NODAL DISPLACEMENT INCREMENTS
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     S(NWS)        = ELEMENT STIFFNESS MATRIX (UPPER TRIANG.ROW-WISE)
C     RE(NEF)       = EQUILIBRIUM LOADS AT ELEMENT NODES
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/
C	--------------------------------
C     NAME(2)       = NAME OF ELEMENT MODULE
C     ITYPE         = CODE NUMBER FOR ELEMENT MODULE
C     ISTYP         = ELEMENT SUBTYPE

C     NLOPT         = CODE FOR NONLINEAR OPTION
C     NLOPT=0         LINEAR ANALYSIS
C     NLOPT=1         MATERIALLY NONLINEAR ONLY
C     NLOPT=2,3       TOTAL LAGRANGIAN,UPDATED LAGRANGIAN

C     MTMOD         = CODE FOR MATERIAL MODULE
C     MTMOD=1         LINEAR ELASTIC,ISOTROPIC
C     MTMOD=2         LINEAR ELASTIC,ORTHOTROPIC
C     MTMOD=3         ELASTO-PLASTIC,VON-MISES
C     MTMOD=5         CONCRETE WITH CRACKING

C     NSINC         = FACTOR CONTROLLING NUMBER OF SUBINCREMENTS
C     ITOLEY        = TOLERANCE ON YIELD FUNCTION
C     NELE          = NUMBER OF ELEMENTS IN THIS GROUP
C     NMPS          = NUMBER OF MATERIAL PROPERTY SETS
C     NGPS          = NUMBER OF GEOMETRIC PROPERTY SETS
C     NMP           = NUMBER OF MATERIAL PROPERTIES PER SET
C     NGP           = NUMBER OF GEOMETRIC PROPERTIES PER SET
C     NNM           = MAXIMUM NUMBER OF NODES FOR ANY ONE ELEMENT
C     NEX           = MAXIMUM NUMBER OF EXCESS NODES
C     NCO           = NUMBER OF NODAL COORDINATES
C     NNF           = NUMBER OF NODAL DEGREES OF FREEDOM
C     NEF           = MAXIMUM NUMBER OF ELEMENT DEGREES OF FREEDOM
C     NWG           = NUMBER OF STORAGE LOCATIONS AT EACH GAUSS POINT
C     NPT           = NUMBER OF GAUSS POINTS
C     NWA           = SIZE OF WORKING ARRAY
C     NWS           = SIZE OF ELEMENT STIFFNESS MATRIX
C     MEL           = CURRENT ELEMENT NUMBER
C     NNO           = NUMBER OF NODES FOR THIS ELEMENT (4)
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR THIS ELEMENT
C     NELTOT        = TOTAL NUMBER OF ELEMENTS (ALL GROUPS)
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /GAUS/
C	--------------------------------
C     GLOC(4,4)     = NATURAL GAUSS POINT COORDINATES (1*1 TO 4*4)
C     GWT (4,4)     = GAUSS POINT WEIGHTS
C     NGR,NGS,NGT   = NUMBER OF GAUSS POINTS IN RN,SN,TN DIRECTION
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /FLAG/
C	--------------------------------
C     IFPRI,ISPRI   = FLAG FOR PRINTING DISPL.OR STRESSES (ISPRI=0)
C     IFPLO         = FLAG FOR PLOT OUTPUT
C     IFREF         = FLAG FOR REFORMATION OF STIFFNESS (IFREF=0)
C     IFEIG         = FLAG FOR EIGENVALUE SOLUTION (IFEIG=0)
C     ITASK = 1       FIRST ENTRY INTO ELEMENT MODULE
C     ITASK = 2       ENTRY DURING EQUILIBRIUM ITERATIONS
C     ITASK = 3       ENTRY TO WORK OUT STRESSES (LAST STEP ONLY)
C     ITASK = 4       ENTRY TO DETERMINE GEOMETRIC STIFF.MATRIX ONLY
C     KSTEP           CURRENT STEP NUMBER
C     KITE            CURRENT ITERATION NUMBER
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DP(36)        = ELASTIC OR ELASTO-PLASTIC MATERIAL MATRIX
C     H(21)         = INTERPOLATION FUNCTIONS
C     HD(3,21)      = SHAPE FUNCTION DERIVATIVES WITH RESP.TO R,S,T
C     XJ(3,3)       = JACOBIAN MATRIX
C     XJI(9)        = INVERSE OF THE JACOBIAN MATRIX
C     B(2*NNO)      = COMPRESSED STRAIN-DISPLACEMENT MATRIX
C     DISD(9)       = DISPLACEMENT DERIVATIVES
C     EPS(6)        = GAUSS POINT STRAINS
C     EPSQ(6)       = QUADRATIC PART OF GAUSS POINT STRAINS
C     SIG(6)        = GAUSS POINT STRESSES
C     IPEL          = SECTION PLASTICITY INDICATER (1=ELASTIC,2=PL)
C     DET           = DETERMINANT OF THE INVERSE JACOBIAN
C     DVOL          = INTEGRATION CONSTANT
C     ----------------
C     VARIABLES OF EAS METHOD
C     ----------------
C     MM             =EAS TERMS OF ALPHA(9, 15,21 OR 30)
C     XJO(3,3)       =JACOBIAN MATRIX AT THE ORIGIN (R,S,T=0.0)
C     XJ(3,3)        =JACOBIAN MATRIX
C     TTO(6,6)       =COEFFICIENT MATRIX OF ENHANCED STRAIN
C     TM(6,MM)       =COEFFICIENT MATRIX OF ENHANCED STRAIN
C     SED(MM,MM)     =MATRIX OF EAS, TRANSPOSE(M)*E*M 
C     SEL(MM,24)     =MATRIX OF EAS, TRANSPOSE(M)*E*B
C     SEDI(MM,MM)    =INVERSE OF MATRIX SED
C     EAS(6,1)       =ENHANCED STRAIN VECTOR
C     ALPHA(MM,1)    =ENHANCED TERMS ALPHA
C     TMT(MM,6)      =TRANSPOSE OF MATRIX TM(6,MM)
C     RE1(24)        =EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     RE2(24,1)      =EQUILIBRIUM FORCE OF EAS METHOD
C     ----------------------------------------------------------------


C	BASIC COMMON FOR SUBROUTINE CONTROL VARIABLE INPUT 
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST

	COMMON /MMENH/ MM,MM1,MM2,NDIMC

C	DIMENSION FOR BASIC PARAMETER
      DIMENSION PROPM(1),PROPG(1),NODEX(1),WA(MWG,1),S(1),COORD(NCO,1)
      DIMENSION EDIS(NEF),EDISI(NEF),RE(NEF),EDISR(NEFU),EDISP(NNO)
	DIMENSION H(8),P(2,8),B(4,16),BP(2,16),HP(2,8)
      DIMENSION DISD(5),STRAIN(4),STRESS(4),TAU(4)
      DIMENSION D(6,6),DE4(6,6),DMM(6,6)
	DIMENSION RTRAN(NEF,NEF),SE(NEF,NEF)
	DIMENSION SD(NEFU,NEFU),SS((NEFU*NEFU+NEFU)/2)
C     ------------
C      EAS METHOD
C     ------------    
	DIMENSION XJ(2,2),XJI(2,2),XJO(2,2),TTO(3,3),TT(3,3)
	DIMENSION EAS(4),ALPHA(MM),KK(4)
      DIMENSION SED(MM,MM),SEL(MM,NEFU),SEDI(MM,MM),RH(MM),SEP(NEFU,NNO)
	DIMENSION GE(4,MM),EM1(4,MM),EM2(4,NEFU),EM3(MM)


C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)

	DIMENSION FIN(NEF)


C	REARRANGE MATRIX
	RTRAN(1:NEF,1:NEF) = 0.0D0
	RTRAN(1 ,1 ) = 1.0D0
	RTRAN(2 ,2 ) = 1.0D0
	RTRAN(3 ,4 ) = 1.0D0
	RTRAN(4 ,5 ) = 1.0D0
	RTRAN(5 ,7 ) = 1.0D0
	RTRAN(6 ,8 ) = 1.0D0
	RTRAN(7 ,10) = 1.0D0
	RTRAN(8 ,11) = 1.0D0
	RTRAN(9 ,3 ) = 1.0D0
	RTRAN(10,6 ) = 1.0D0
	RTRAN(11,9 ) = 1.0D0
	RTRAN(12,12) = 1.0D0


	TH = PROPG(2)

	ISTYP1 = ISTYP

C	ISTYP FOR PLANE STRESS - PLANE STRAIN FLAG
	ISTYP = ISTYP - 3 !SONGSAK SEP2005


	EDISR(1:NEFU) = 0.0D0
	EDISP(1:NNO ) = 0.0D0

	CALL EXTRACTDISP (EDIS,EDISR)  !ELEMENT DISPLCEMENT
	CALL EXTRACTPRES (EDIS,EDISP)  !PORE PRESSURE

C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------
      CALL HOKLAW (PROPM,PROPG,ISTYP)
	CALL DELA3D (D,ISTYP)
	RHO = PROPM(5)

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
C     -------------------------------------------------------------
C     COMPUTE JACOBIAN MATRIX AT THE ELEMENT CENTER (R,S)=(0.0,0.0)
C     -------------------------------------------------------------

      CALL SHAP2D(0.0D0,0.0D0,H,P,NODEX,NNO)
	CALL JACOB2D(MEL,NNO,COORD,P,XJO,XJI,DETO,NCO)	

C     -------------------------------------------------
C     COMPUTE [TTO] MATRIX (FOR EAS METHOD) USING [XJO]
C     -------------------------------------------------

C	COMPUTE [TTO]**(-1)

	CALL TTOMAT(TTO,XJO)
      CALL INVMAT1(TTO,TT,3)

C	---------------------
C	COMPUTE [TTO]**(-T)
C	---------------------

	TTO = TRANSPOSE(TT)

C     -------------------------------
C     OBTAIN LINEAR STRESS-STRAIN LAW
C     -------------------------------
      ISR = 3
      IST = 3
      IF (ISTYP.LE.1) IST = 4
      IF (ISTYP.EQ.0) ISR = 4
      MGR = NGR
      MGS = NGS

      IF (ITASK.EQ.5) THEN
      MGR = 3
      MGS = 3
	ENDIF

	IPT = 0
C     -----------------------------------------
C      INITIALIZATION OF EAS VARIABLE MATRICES
C     -----------------------------------------
	RH  = 0.0D0
	SED = 0.0D0
	SEL = 0.0D0
	SEP = 0.0D0
	SS  = 0.0D0
	SD  = 0.0D0
	SE  = 0.0D0
C     ------------------------
C      LOOP OVER GAUSS POINTS
C     ------------------------

      DO 900  IGR=1,MGR
	RI = GLOC(IGR,MGR)
	DO 900  IGS=1,MGS
	SI = GLOC(IGS,MGS)
      
	WT  = GWT(IGR,MGR)*GWT(IGS,MGS)
	IPT = IPT+1

C     ---------------------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ---------------------------------------------------

      CALL SHAP2D (RI,SI,H,P,NODEX,NNO)
	CALL JACOB2D(MEL,NNO,COORD,P,XJ,XJI,DET,NCO)
	CALL MEBMAT (COORD,H,P,XJI,B,ISTYP,XBAR,NNO)

	HP(1,1:8) = H(1:8)
	HP(2,1:8) = H(1:8)
	BP(1,1:16)  = B(1,1:16)
	BP(2,1:16)  = B(2,1:16)
	B(1:2,1:16) = 0.0D0
	B(  4,1:16) = 0.0D0


      IF (ISTYP.NE.0) XBAR = TH
      DVOL = WT*DET*XBAR
C     -------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX
C     -------------------------------
      IF (ITASK.NE.5)GOTO 50
      CALL MASS2D (SS,H,RHO,DVOL,NNO,NEFU,IPT)
      GOTO 900
50	CONTINUE

C     -------------------------------
C     ADD CONTRIBUTION TO DAMP MATRIX
C     -------------------------------
	CALL DAMP2D (SD,H,P,RHO,DVOL,NNO,NEFU,IPT,EDISR)

C     ---------------------------------------------     
C     FOR ENHANCED STRAIN METHOD 
C     COMPUTE MATRIX [GE] AND TRANSPOSE [GE] MATRIX
C     ---------------------------------------------
	R12 = H(1)*DD1 + H(2)*DD2 + H(3)*DD3 + H(4)*DD4
C	PLANE
	IF(ISTYP.NE.0) CALL GEMATFLU1(GE,RI,SI,TTO,MM,DETO,DET,MEL)	
C	AXISYMMETRIC
	IF(ISTYP.EQ.0) CALL GEMATFLU2(GE,RI,SI,TTO,MM,DETO,DET,MEL,
	1				              R12,ETA1,ETA2)

C     ------------------------------------
C     D=TRANSPOSE(GE)*E*GE = [D](MP,MP)
C     L=TRANSPOSE(GE)*E*BB = [L](MP,2*NNO)
C     ------------------------------------
	KK(1:4) = [1,2,3,6]
	DO I = 1,ISR
	N = KK(I)
	DO J = 1,MM
	EM1(I,J) = 0.0D0
	DO K = 1,ISR
	IF(K.LE.3) EM1(I,J) = EM1(I,J) + D(I,K)*GE(K,J)
	IF(K.EQ.4) EM1(I,J) = EM1(I,J) + D(N,6)*GE(4,J)
	ENDDO
	ENDDO
	DO J = 1,NEFU
	EM2(I,J) = 0.0D0
	DO K = 1,ISR
	IF(K.LE.3) EM2(I,J) = EM2(I,J) + D(I,K)*B(K,J)
	IF(K.EQ.4) EM2(I,J) = EM2(I,J) + D(N,6)*B(4,J)
	ENDDO
	ENDDO
	ENDDO


	DO I = 1,MM
	DO J = 1,MM
	DO K = 1,ISR
	SED(I,J) = SED(I,J) + GE(K,I)*EM1(K,J)*DVOL
	ENDDO
	ENDDO
	DO J = 1,NEFU
	DO K = 1,ISR
	SEL(I,J) = SEL(I,J) + GE(K,I)*EM2(K,J)*DVOL
	ENDDO
	ENDDO
	ENDDO

	DO I = 1,NEFU
	DO J = 1,NNO
	DO K = 1,2
	SEP(I,J) = SEP(I,J) + BP(K,I)*HP(K,J)*DVOL
	ENDDO
	ENDDO
	ENDDO
C     ----------------------------------------
C     ----------------------------------------
	IF (ITASK.EQ.1)  GOTO 700

C     -------------------------------
C     DISPLACEMENT DERIVATIVES (DISD)
C     -------------------------------

200	CALL CLEARA (DISD,5)
      DO 250  J=1,NEFU,2
      J1 = J+1
      DISD(1) = DISD(1) + B(1,J) *EDISR(J)
      DISD(2) = DISD(2) + B(2,J1)*EDISR(J1)
      DISD(3) = DISD(3) + B(3,J) *EDISR(J)
250	DISD(4) = DISD(4) + B(3,J1)*EDISR(J1)
      IF (ISTYP.GE.1)  GOTO 300
      DO 290  IEF=1,NEFU,2
290	DISD(5) = DISD(5) + B(4,IEF)*EDISR(IEF)
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
300	STRAIN(1) = 0.0D0
      STRAIN(2) = 0.0D0
      STRAIN(3) = DISD(3)+DISD(4)
      STRAIN(4) = 0.0D0

C     ---------------------------------
C      FOR EAS TERM (EPS(EAS)=GE*ALPHA)
C     ---------------------------------
	DO I = 1,ISR
	EAS(I) = 0.0
	DO J = 1,MM
	EAS(I) = EAS(I) + GE(I,J)*ALPHA(J)
	ENDDO
	ENDDO

	DO I=1,ISR
	STRAIN(I)=STRAIN(I) - EAS(I)  
	ENDDO 

C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------
400	STRESS(1:4) = 0.0D0
	KK(1:4) = [1,2,3,6]
	DO I = 1,4
	II = KK(I)
	DO J = 1,4
	JJ = KK(J)
	STRESS(I) = STRESS(I) + D(II,JJ)*STRAIN(J)
	ENDDO
	ENDDO

	DO I=1,4
	WA(I,IPT) = STRESS(I)
	END DO


C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
	DO 510  I=1,IST
510	TAU(I) = STRESS(I)*DVOL
C
520	DO 550  J=1,NEFU,2
      J1 = J+1
      RE(J)  = RE(J)  + B(1,J) *TAU(1) + B(3,J) *TAU(3)
550	RE(J1) = RE(J1) + B(2,J1)*TAU(2) + B(3,J1)*TAU(3)
      IF (ISTYP.NE.0)  GOTO 595
      DO 590  IEF=1,NEFU,2
590	RE(IEF) = RE(IEF) + B(4,IEF)*TAU(4)

C     ------------------------------------------------------
C      COMPUTE ENHANCED EQUILIBRIUM FORCE {h} FOR EAS METHOD
C     ------------------------------------------------------
595	CONTINUE
	DO I = 1,MM
	DO J = 1,ISR
	RH(I) = RH(I) + GE(J,I)*TAU(J)
	ENDDO
	ENDDO

C     ------------------------------------------------------
C     FIND LINEAR CONTRIBUTION TO STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------------------
600	IF (IFREF) 900,700,900

700	CALL MEK0PL (SS,DMM,B,DVOL,NEFU,ISTYP)
	
900	CONTINUE
C     --------------------------------------------------------
C     --------------------------------------------------------


	IF(ITASK.EQ.5) GOTO 920   !MASS & DAMP MATRIX

C	SOLVE FOR EAS TERM
	CALL INVMAT1(SED,SEDI,MM)
C     ------------------------------------
C     FOR ENHANCED STRAIN METHOD
C     K=TRANSPOSE(L)*INVERSE(D)*L
C     K(FINAL)=K.Linear(PURE DISP.)-K(EAS)
C     ------------------------------------
	CALL MATEAS(SEL,SEDI,SS,MM,NNO,NEFU)

C     ------------------------------------
C     FOR NONLINEAR ANALYSIS OF EAS METHOD
C     ------------------------------------
	IF(ITASK.EQ.1) GO TO 920

	DO I = 1,MM
	EM3(I) = 0.0D0
	DO J = 1,MM
	EM3(I) = EM3(I) + SEDI(I,J)*RH(J)
	ENDDO
	ENDDO
	
	DO I = 1,NEFU
	DUM = 0.0D0
	DO J = 1,MM
	DUM = DUM + SEL(J,I)*EM3(J)
	ENDDO
	RE(I) = RE(I) + DUM
	ENDDO
	RE(NEFU+1:NEFU+NNO) = EDISP(1:NNO)

	RE = MATMUL(TRANSPOSE(RTRAN),RE)


920	CONTINUE
	KEF = 0
	DO IEF = 1  ,NEFU
	DO JEF = IEF,NEFU
	KEF = KEF + 1
	SE(IEF,JEF) = SS(KEF)
	ENDDO
	ENDDO

	DO IEF = 1  ,NEFU
	DO INO = 1  ,NNO
	SE(IEF,INO+NEFU) = SEP(IEF,INO)
	ENDDO
	ENDDO
	DO IEF = 1  ,NEF
	DO JEF = IEF,NEF
	SE(JEF,IEF) = SE(IEF,JEF)
	ENDDO
	ENDDO
	DO IEF = 1  ,NEFU
	DO JEF = 1  ,NEFU
	SE(IEF,JEF) = SE(IEF,JEF) + SD(IEF,JEF)
	ENDDO
	ENDDO

	SE = MATMUL(TRANSPOSE(RTRAN),MATMUL(SE,RTRAN))
	KEF = 0
	DO IEF = 1  ,NEF
	DO JEF = IEF,NEF
	KEF = KEF + 1
	S(KEF) = SE(IEF,JEF) 
	ENDDO
	ENDDO


	MEF = KEF

	DO IEF = 1  ,NEF
	DO JEF = IEF,NEF
	KEF = KEF + 1
	S(KEF) = SE(JEF,IEF) 
	ENDDO
	ENDDO

	KEF = 0
	DO IEF = 1  ,NEF
	DO JEF = IEF,NEF
	KEF = KEF + 1
C	WRITE(*,*) 'K1',S(KEF)-S(MEF+KEF)
	ENDDO
	ENDDO

C	WRITE(*,*) SE(1,[3,6,9,12])
C	PAUSE


	IF (ITASK.EQ.3) THEN
	DO I = 1,NEF
	FIN(I) = RE(I)
	ENDDO
	ENDIF

C	RECALL ISTYP
	ISTYP = ISTYP1 !SONGSAK SEP2005
	

	RETURN
	END

C	==============================================================
C	==============================================================
C	==============================================================
      SUBROUTINE DAMP2D (DP,H,P,DE,DVOL,NNO,NEF,IPT,EDIS)
      IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------------------
C	--------------------------------------------------------------
C
      DIMENSION DP(NEF,NEF),H(1),P(2,1),EDIS(1),ED(2,2)
	DIMENSION HH(2,NEF),PP(2,NEF)
C
	FAC=DVOL*DE

	DO I = 1,NNO
	HH(1,2*I-1) = H(I)
	HH(2,2*I-0) = H(I)
	PP(1,2*I-1) = P(1,I)
	PP(2,2*I-0) = P(2,I)
	ENDDO

	ED(1:2,1:2) = 0.0D0
	DO I = 1,2
	DO J = 1,NEF
	ED(I,I) = ED(I,I) + HH(I,J)*EDIS(J)
	ENDDO
	ENDDO

	DP = MATMUL(TRANSPOSE(HH),MATMUL(ED,PP))*FAC

C
      END

C
C	==============================================================
C	==============================================================
C	==============================================================
      SUBROUTINE GEMATFLU1(GE,RI,SI,TT,MP,DETJO,DETJ,MEL)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	-------------------------------------------------------------------------
C		FOR ENHANCED STRAIN TERM OF ENHANCED STRAIN METHOD

C		COMPUTE COEFFICIENT MATRIX EM(3,MP)=[TT]*[M]

C		COMPUTE ENHANCED MATRIX [GE] = FAC*[EM(I,J)] 

C		IN MASTER OR LOCAL COORDINATE SYSTEM

C		OUTPUT MATRIX [GE]

C	-------------------------------------------------------------------------

	DIMENSION EM(4,MP),TT(3,3),GE(4,MP),GM(4,MP)
      
	DO 10 I =1,4
		DO 10 J =1,MP
			EM(I,J) = 0.0
			GE(I,J) = 0.0
			GM(I,J) = 0.0 
10	CONTINUE

C	...ENHANCED 2 TERM

	IF (MP.EQ.2) THEN
	  GM(3,1) = RI
	  GM(3,2) = SI  
	  GOTO 100    
	END IF

C	...ENHANCED 3 TERM

	IF (MP.EQ.3) THEN
        GM(3,1) = RI
        GM(3,2) = SI
        GM(3,3) = RI*SI
        GOTO 100    
	END IF

100	CONTINUE

	DO 110 I = 1,3
	DO 110 J = 1,MP
	EM(I,J) = 0.0D0
	DO 110 K = 1,3
110	EM(I,J) = EM(I,J) + TT(I,K)*GM(K,J)

	CONST = (DETJO/DETJ) 

	DO 200 IEM = 1,4
		DO 200 JEM = 1,MP 
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
200	CONTINUE


	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE GEMATFLU2 (GE,RI,SI,TT,MP,DETJO,DETJ,MEL,R12,ETA1,ETA2)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	-------------------------------------------------------------------------
C		FOR ENHANCED STRAIN TERM OF ENHANCED STRAIN METHOD

C		COMPUTE COEFFICIENT MATRIX EM(3,MP)=[TT]*[M]

C		COMPUTE ENHANCED MATRIX [GE] = FAC*[EM(I,J)] 

C		IN MASTER OR LOCAL COORDINATE SYSTEM

C		OUTPUT MATRIX [GE]

C	-------------------------------------------------------------------------

	DIMENSION EM(4,MP),TT(3,3),GE(4,MP),GM(4,MP)
      
	DO 10 I =1,4
		DO 10 J =1,MP
			EM(I,J) = 0.0
			GE(I,J) = 0.0
			GM(I,J) = 0.0 
10	CONTINUE



C	...ENHANCED 2 TERM

	IF (MP.EQ.2) THEN
        GM(3,1) = RI - ETA1
        GM(3,2) = SI - ETA2
	END IF

	DO 100 I = 1,3
	DO 100 J = 1,MP
	EM(I,J) = 0.0D0
	DO 100 K = 1,3
100	EM(I,J) = EM(I,J) + TT(I,K)*GM(K,J)
	EM(4,5) = GM(4,5)

	CONST = (DETJO/DETJ) 


	DO 200 IEM = 1,4
		DO 200 JEM = 1,MP 
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
200	CONTINUE


	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
