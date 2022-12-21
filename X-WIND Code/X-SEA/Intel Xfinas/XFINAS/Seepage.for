C	=================================================================
      SUBROUTINE SEEPAG2D(PROPM,PROPG,NODEX,WA,S,COORD,
	1                  EDIS,EDISI,RE,MWG,FIN)
C	FIN - ADDED TO PREVIOUS LINE BY GILSON - JUL2003 (INT FORCE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     COMPUTES ELEMENT FORCE VECTOR BALANCED IN CURRENT DISPLACEMENT
C     FIELD. IN CASE IFREF=0 THE STIFFNESS MATRIX IS UPDATED FOR
C	--------------------------------------
C     MATERIAL NONLINEARITIES    (NOLIN=1,3)
C     GEOMETRIC NONLINEARITIES   (NOLIN=2,3)
C     FINDING LINEAR STIFFNESS   (NOLIN=0)
C	---------------
C     INPUT VARIABLES
C	---------------
C     COORD(3,8)     = NODAL COORDINATES X,Y Z
C     PROP(5)        = MEMBER PROPERTIES (SEE ROUTINE CONSTI)
C     NODEX(NNO-4)   = POSITIONS OF MIDSIDE NODES
C     EDIS(24)       = TOTAL DISPLACEMENTS AT ELEMENT NODES
C     EDISI(24)      = NODAL ELEMENT DISPLACEMENT INCREMENTS
C     IE             = ELEMENT NUMBER
C     NI             = NUMBER OF NODES FOR THIS ELEMENT
C     MW             = DIMENSION OF WORKING ARRAY STORING 4 STRESSES
C                      (PLUS 4 STRAINS AND A PLASTICITY FLAG NOLIN=1,3)
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     RE(16)         = ELEMENT FORCE VECTOR (TWO EQUIL. LOADS PER NODE)
C     S(300)         = NEW TANGENTIAL STIFFNESS MATRIX (UPPER TRIANG)
C     WA(MW,NGR*NGS) = WORKING ARRAY STORING STRESS,(STRAIN) COMPONENTS
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DISD(5)        = DISPLACEMENT DERIVATIVES
C     STRAIN(4)      = LINEAR OR NONLINEAR STRAIN VECTOR
C     STRESS(4)      = LINEAR OR NONLINEAR STRESS VECTOR
C     D(6,6)         = STRESS - STRAIN MATRIX
C	NLOPT(INPUT DATA)	= 0-LINEAR ANALYSIS
C						= 1-MATERIAL NONLINEAR ONLY
C						= 2-TOTAL LAGRANGIAN
C						= 3-UPDATED LAGRANGIAN
C	ISOLOP(INPUT DATA)	= 1-LINEAR ANALYSIS
C						= 2-NONLINEAR STATIC
C						= 3-EIGN VALUES(FREEQUANCY ANALYSIS)
C						= 4-EIGN VALUES(LINEAR BUCKLING ANALYSIS)
C	IFREF				= 0-FLAG FOR REFORMATION OF STIFFNESS
C	DPER(2,2)			= PERMIABILITY MATRIX
C	ITASK				= 1-FIRST ENTRY INTO ELEMENT MODULUS
C						= 2-ENTRY DURING EQULIBRIUM ITERATIONS
C						= 3-ENTRY TO WORK OUT STRESS
C						= 4-ENTRY TO DETERMINE GEOMETRIC STIFFNESS
C						= 5-CONTRIBUTION TO MASS MATRIX
C						= 6-CONTRIBUTION TO DAMPING MATRIX
C	DTINC				= Dt TIME INCREMENT
C	NTSTEP				= TOTAL TIME STEPS
C	KTSTEP				= CURRENT TIME STEP
C	CTIME				= CURRENT TIME	
C     ----------------------------------------------------------------

      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES

C	ADD THIS LINE FOR ITASK SONGSAK MAY2006
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
      DIMENSION PROPM(5),PROPG(2),WA(MWG,1),S(300),COORD(16),EDIS(24)
      DIMENSION EDISI(24),RE(24),NODEX(1),D(3,3),H(8),P(2,8),XJI(2,2)
      DIMENSION B(3,16),DISD(4),STRAIN(4),STRESS(4),TAU(4)
	DIMENSION SM(16,16),CC(16,4),SP(4,4),SS(20,20)
	DIMENSION HF(4),PF(2,4),PERM(3,3)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
C     ----------------------------
C     INITIALIZE INTEGRATION RULE
C     ---------------------------
	CALL CLEAR (16,16,SM)
	CALL CLEAR (16,4,CC)
	CALL CLEAR (4,4,SP)
	CALL CLEAR (20,20,SS)
C	CONSIDER 8 NODE MEMBRAIN SOLID AND 4 NODE PORE PRESSURE
C	LET CONSIDER MGR = 2,MGS=2
      MGR = NGR
      MGS = NGS
c	write(100,*)coord
C	LOAD COMMON BLOCK /HOOK/
C	-------------------------------------
      CALL SHOKLAW (PROPM,PROPG,PERM,WADEN)
	

C     -------------------------------
C     OBTAIN LINEAR STRESS-STRAIN LAW
C     -------------------------------
      CALL SDELA3D (D)
c	print*,th
c	write(100,*)d

	IPT = 0
C     -----------------------------------------------------------------
C     LOOP OVER GAUSS POINTS,	START GAUSS LOOP
C	-----------------------------------------------------------------
c	WRITE(100,*) coord
      DO 900  IGR=1,MGR	
      RI = GLOC(IGR,NGR)	
      DO 900  IGS=1,MGS	
      SI = GLOC(IGS,NGS)	
      WT = GWT(IGR,NGR) * GWT(IGS,NGS)
	IPT = IPT + 1	
c	print*,wt,ri,si
C     ----------------------------------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P), JACOBIAN (XJ), ITS INVERSE
C     (XJI) AND DETERMINANT (DET), DISPLACEMENT DERIVATIVES-MATRIX (B)
C     ----------------------------------------------------------------
C	CALCULATE H,P,XJ,XJI,DET AND B FOR SOLID PART 8 NODE ELEMENT
C	------------------------------------------------------------
      CALL SSHAPD (RI,SI,H,P,8)
      CALL SJACOD (MEL,8,COORD,P,XJI,DET)
      CALL SMEBMAT (P,XJI,B,8)
      DVOL = WT*DET*TH

C	CALCULATE THE STRESS AT GAUSS POINTS AND STORE THEM
	IF(KTSTEP.LE.1) GOTO 100
C     DISPLACEMENT DERIVATIVES
C     ------------------------
      CALL CLEARA (DISD,4)
	K = 0
      DO 250  J=1,16,2
      J1 = J+1
	K = K+1
      DISD(1) = DISD(1) + B(1,J) *EDIS(K)
      DISD(2) = DISD(2) + B(2,J1)*EDIS(K+1)
      DISD(3) = DISD(3) + B(3,J) *EDIS(K)
      DISD(4) = DISD(4) + B(3,J1)*EDIS(K+1)
	K = K+2
 250	CONTINUE

C     LINEAR STRAIN TERMS
C     -------------------
      STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)+DISD(4)
      STRAIN(4) = 0.
C     COMPUTE AND STORE STRESSES
C     ------------------------------------
      CALL SMEMSIG (STRAIN,STRESS)
      DO 415  I=1,4
 415  WA(I,IPT) = STRESS(I)
C	WRITE(100,*) IPT,MEL,MWG
C	WRITE(100,*)STRAIN
C	PRINT LAST STEP STESS ONLY
C	--------------------------
	IF(KTSTEP.GT.NTSTEP) GOTO 900

C	CALCULATE HERE SOLID CONTRIBUTION (KM)
C	--------------------------------------
  100 CALL SMEKPL (SM,D,B,DVOL)

C	CALCULATE HF,PF,XJF,XJIF,DETF FOR FLUID PART 4 NODE ELEMENT
C	-----------------------------------------------------------
      CALL SSHAPD (RI,SI,HF,PF,4)

C	ADD HERE CONTRIBUTION FROM COPLING TERM (C)
C	-------------------------------------------
	CALL SCOUPL (CC,HF,B,DVOL)

C	ADD HERE CONTRIBUTION FROM SEEPAGE (KP)
C	---------------------------------------
	CALL SEEPCON (SP,XJI,PF,PERM,DVOL)

C	ADD HERE NON LINEAR CONTRIBUTION, IF NLOPT>1
 900  CONTINUE  
C	END OF GAUSS POINT LOOP
C	-----------------------
C	ASSEMBLE (KM,C,KP) IN GLOBAL MATRIX
C	-----------------------------------
	CALL SASSEMB (SS,SM,CC,SP)

C	ARRANGE SS AND MAKE UPPER TRANGULAR MATRIX (S)
C	-------------------------------------------
	CALL SARRENG (S,SS)

C	-------------------------------
C	FIND THE RESULTANT FORCE VECTOR 
C	-------------------------------
C	CONTRIBUTION FROM THE Nth TIME STEP AND CONTRIBUTION
C	FROM THE FIXED PORE PRESSURE
C	----------------------------------------------------
	CALL SNTHRE (SS,SM,CC,SP,RE,EDIS)
c	WRITE(100,*) ss
C	WRITE(7,*)
C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN !CALCULATE STERSS AT THE LAST STEP
	DO 2000 I = 1,NEF
	FIN(I) = RE(I)
2000	CONTINUE
	ENDIF

      RETURN
      END
C
      SUBROUTINE SHOKLAW (PROPM,PROPG,PERM,WADEN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     DETERMINES ALL THE MATERIAL CONSTANTS IN COMMON BLOCK /HOOK/
C	------------------------------------------------------------
C     PROPM(NMP)    = MATERIAL PROPERTIES
C     PROPG(NGP)    = GEOMETRIC PROPERTIES
C     IND           = FLAG FOR PLANE STRAIN(1),OR PLANE STRESS(2)
C     ------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION PROPM(7),PROPG(*),PERM(3,3)
C
	CALL CLEAR (3,3,PERM)

      YM  = PROPM(1)
      PR  = PROPM(2)
      YLD = PROPM(3)
      HP  = PROPM(4)
      DEN = PROPM(5)

	WADEN = PROPM(9)

	XK  = PROPM(6)/WADEN
	YK  = PROPM(7)/WADEN
	ZK  = PROPM(8)/WADEN
C
      D2  = PR/(PR-1.)
      C2  = YM/(1.+PR)
      B2  = PR/(1.-2.*PR)
      A2  = (1.-PR)/(1.-2.*PR)
      BM  = YM/(3.-6.*PR)
      C1  = .5*C2
      D1  = C1/1.2
C
C     ----------------------
C     PLANE STRAIN CONDITION
C     ----------------------
      B1 = B2*C2
      A1 = C2+B1

C     -----------------------
C     NODAL ELEMENT THICKNESS
C     -----------------------
      TH  = PROPG(2)

C	------------------------
C	PERMIABILITY MATRIX(has called here Kx/gamaW
C	------------------------
	PERM(1,1) = XK
	PERM(2,2) = YK
	PERM(3,3) = ZK


      RETURN
      END
C
C	===============================================================
      SUBROUTINE SDELA3D (D)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------
C     GENERATES THE THREE-DIMENSIONAL STRESS-STRAIN MATRIX FOR A
C     LINEAR ELASTIC,ISOTROPIC MATERIAL
C	---------------------------------
C     D(6,6)  = ELASTIC,ISOTROPIC STRESS-STRAIN MATRIX
C     ----------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION D(3,3)
C
      CALL CLEAR (3,3,D)
C     ----------------------
C     PLANE STRAIN CONDITION
C     ----------------------
      D(1,1)  = A1
      D(1,2)  = B1
      D(2,1)  = B1
      D(2,2)  = A1
      D(3,3)  = C1

      RETURN
      END
C	================================================================
      SUBROUTINE SSHAPD (R,S,H,P,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THERE DERIVATIVES
C     AT THE NODAL POINTS OF A 4 TO 8 NODE,ISOPARAMETRIC QUADRILATERAL
C	----------------------------------------------------------------
C                         NODE NUMBERING CONVENTION
C                         -------------------------
C
C                   2             5             1
C
C                     0 . . . . . 0  . . . . . 0
C                     .                        .
C                     .           S            .
C                     .           **            .
C                   6 0           .  > R       0 8
C                     .                        .
C                     .                        .
C                     .                        .
C                     0 . . . . . 0  . . ... . 0
C
C                   3             7              4
C
C     R,S        = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     H(8)       = INTERPOLATION (SHAPE) FUNCTIONS
C     P(2,8)     = FUNCTION DERIVATIVES WITH RESPECT TO R,S RESP.
C     NODEX(NEX) = POSITION OF MIDSIDE (EXCESSIVE) NODES
C     NNO        = NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ----------------------------------------------------------------
      DIMENSION  H(NNO),P(2,NNO)

	CALL CLEAR (2,NNO,P)
	CALL CLEARA(H,NNO)

	IF (NNO.GT.4) GOTO 100
      RP  = 1.0+R
      SP  = 1.0+S
      RM  = 1.0-R
      SM  = 1.0-S
C     ---------------------------------------------
C     INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     FOR A FOUR NODE ELEMENT
C	H1 = 0.25(1+r)(1+s)
C	H2 = 0.25(1-r)(1+s)
C	H3 = 0.25(1-r)(1-s)
C	H4 = 0.25(1+r)(1-s)
C	P(1,*) DERIVATIVE wrt r DIRECTION
C	P(2,*) DERIVATIVE wrt s DIRECTION
C     ---------------------------------------------
      H(1)   = 0.25*RP*SP
      H(2)   = 0.25*RM*SP
      H(3)   = 0.25*RM*SM
      H(4)   = 0.25*RP*SM
C	SHAPEFUNCTION DERIVATIVES wrt r
      P(1,1) = 0.25*SP
      P(1,2) = -P(1,1)
      P(1,3) = -0.25*SM
      P(1,4) = -P(1,3)
C	SHAPE FUNCTION DERIVATIVES wrt s
      P(2,1) = 0.25*RP
      P(2,2) = 0.25*RM
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)
	RETURN
C     ----------------------------------------------
C	IF 8-NODE
C	H1=0.25(1+r)(1+s)(r+s-1)
C	H2=0.25(1-r)(1+s)(-r+s-1)
C	H3=0.25(1-r)(1-s)(-r-s-1)
C	H4=0.25(1+r)(1-s)(r-s-1)
C	H5=0.5(1-r**2)(1+s)
C	H6=0.5(1-s**2)(1-r)
C	H7=0.5(1-r**2)(1-s)
C	H8=0.5(1-s**2)(1+r)
C     ----------------------------------------------
 100  RP  = 1.0+R
      SP  = 1.0+S
      RM  = 1.0-R
      SM  = 1.0-S
      R2  = 1.0-R*R
      S2  = 1.0-S*S
      H(1)   = 0.25*RP*SP*(R+S-1.0)
      H(2)   = 0.25*RM*SP*(-R+S-1.0)
      H(3)   = 0.25*RM*SM*(-R-S-1.0)
      H(4)   = 0.25*RP*SM*(R-S-1.0)
      H(5)   = 0.5*R2*SP
      H(6)   = 0.5*S2*RM
      H(7)   = 0.5*R2*SM
      H(8)   = 0.5*S2*RP
C	SHAPE FUNCTION DERIVATIVES wrt r
      P(1,1) = 0.25*SP*(2.0*R+S)
      P(1,2) = 0.25*SP*(2.0*R-S)
      P(1,3) = 0.25*SM*(2.0*R+S)
      P(1,4) = 0.25*SM*(2.0*R-S)
	P(1,5) = -R*SP
      P(1,6) = -0.5*S2
      P(1,7) = -R*SM
      P(1,8) = 0.5*S2
C	SHAPE FUNCTION DERIVATIVES wrt s
      P(2,1) = 0.25*RP*(2.0*S+R)
      P(2,2) = 0.25*RM*(2.0*S-R)
      P(2,3) = 0.25*RM*(2.0*S+R)
      P(2,4) = 0.25*RP*(2.0*S-R)
      P(2,5) = 0.5*R2
      P(2,6) = -S*RM
      P(2,7) = -.5*R2
      P(2,8) = -S*RP

	RETURN
      END
C
C	=================================================================
C     ------------------------------------------
C	XY - COORDINETES
C	P  - SHAPE FUNCTION DERIVATIVES
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      SUBROUTINE SJACOD (MEL,NNO,XY,P,XJI,DET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION XY(2,NNO),P(2,NNO),XJ(2,2),XJI(2,2)
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
	xj = matmul(p,transpose(xy))
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = XJ(1,1)*XJ(2,2) - XJ(2,1)*XJ(1,2)
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')
C     -----------------------------
C     INVERSE (XJI) OF THE JACOBIAN (COLUMN WISE ARRANGEMENT)
C     -----------------------------
      DUM = 1.0/DET
      XJI(1,1) =  XJ(2,2)*DUM
      XJI(2,1) = -XJ(2,1)*DUM
      XJI(1,2) = -XJ(1,2)*DUM
      XJI(2,2) =  XJ(1,1)*DUM

      RETURN
      END
C
C	================================================================
      SUBROUTINE SMEBMAT (P,XJI,B,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------
C     EVALUATES THE GLOBAL LINEAR STRAIN-DISPLACEMENT MATRIX B
C     FOR A QUADRILATERAL MEMBRANE ELEMENT
C	------------------------------------
C     INPUT VARIABLES
C	---------------
C     H(NNO)     =  SHAPE FUNCTIONS
C     P(2,NNO)   =  SHAPE FUNCTION DERIVATIVES WITH RESP.TO R,S
C						P(1,NNO) W.R.S.TO r
C						P(2,NNO) W.R.S.TO s
C     XJI(4  )   =  INVERSE OF THE JACOBIAN MATRIX(COLUMN WISE)
C     B(4,NNO)   =  STRAIN-DISPLACEMENT MATRIX
C     ISTYP      =  ELEMENT SUBTYPE
C     XBAR       =  RADIUS AT GAUSS POINT
C     NNO        =  NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ---------------------------------------------------------
      DIMENSION P(2,NNO),XJI(2,2),B(3,2*NNO)
	DIMENSION XJIA(2,NNO)
      CALL CLEAR (3,2*NNO,B)

C	CALCULATE DERIVATIVES wrt X AND Y
	XJIA = MATMUL(XJI,P)

C	CALCULATE B MATRIX

	DO I = 1,NNO
	K = 2*I
	L = K-1
	X = XJIA(1,I)
	Y = XJIA(2,I)
	B(1,L) = X
	B(2,K) = Y
	B(3,K) = X
	B(3,L) = Y
	END DO

      RETURN
      END
C
C	================================================================
      SUBROUTINE SMEMSIG (STRAIN,STRESS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CONVERTS GLOBAL MEMBRANE STRAINS TO GLOBAL MEMBRANE STRESSES
C	------------------------------------------------------------
C     STRAIN(4)      = GIVEN STRAIN VECTOR
C     STRESS(4)      = RETURNED STRESS VECTOR
C     ISTYP          = FLAG (AXISYMMETRIC=0,PL.STRAIN=1,PL.STRESS=2)
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION STRAIN(4),STRESS(4)
C
      STRESS(1) = A1*STRAIN(1) + B1*STRAIN(2)
      STRESS(2) = B1*STRAIN(1) + A1*STRAIN(2)
      STRESS(3) = C1*STRAIN(3)
      STRESS(4) = B1*(STRAIN(1)+STRAIN(2))
      STRAIN(4) = 0.
      RETURN
      END
C
C	===============================================================
C     ----------------------------------------------------------------
C     ADDS LINEAR CONTRIBUTION TO ELEMENT STIFFNESS FROM SOLID PART
C	-----------------------------------------------------------
C     S(136)     = ELEMENT STIFFNESS MATRIX IN UPPER TRIANGULAR FORM(8NODE)
C			   = 36 FOR 4 NODE ELEMENT
C     D(6,6)     = STRESS - STRAIN MATRIX (3X3 REPRESENT THE 2D CONSTITUTIVE)
C     B(3,2*NNO) = LINEAR STRAIN-DISPLACEMENT MATRIX
C     FAC        = DET*WT*TH
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ISTYP      = ELEMENT SUPTYPE (0=AXISYM.,1/2=PLANE STRAIN/STRESS)
C     ----------------------------------------------------------------
      SUBROUTINE SMEKPL (SM,D,B,FAC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION B(3,16),DB(4),SM(16,16),D(3,3)

	SM = SM + MATMUL(MATMUL(TRANSPOSE(B),D),B)*FAC	

	RETURN
      END
C
C	================================================================
C     ----------------------------------------------------------------
C     CALCULATE COUPLING PART OT THE LINEAR STIFFNESS MATRIX
C	-----------------------------------------------------------
C     CC(16,8)   = COUPLING TERM
C     B(3,2*NNO) = LINEAR STRAIN-DISPLACEMENT MATRIX
C     FAC        = DET*WT*TH
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C	H(8)	   = SHAPE FUNCTION
C     ----------------------------------------------------------------
      SUBROUTINE SCOUPL (CC,H,B,FAC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION B(3,16),CC(16,8),H(4)

	DO J = 1,4
	DO I = 1,16
	CC(I,J) = CC(I,J) + (B(1,I)+B(2,I))*H(J)*FAC
	END DO
	END DO

	RETURN
      END
C
C	=================================================================
C     ----------------------------------------------------------------
C     CALCULATE THE SEEPAGE CONTRY BUTION TO LINEAR STIFFNESS MATRIX
C	-----------------------------------------------------------
C     SP(4,4)    = CONTRIBUTION FROM WATER
C     DPER(2,2)  = REMIABILITY MATRIX
C     P(2,8)	   = SHAPE FUNCTION DERIVATIVES
C	XJI(2,2)   = INVERSE JACOBIAN
C     FAC        = DET*WT*TH
C     NNO        = NUMBER OF NODE PER ELEMENT
C     ----------------------------------------------------------------
      SUBROUTINE SEEPCON (SP,XJI,P,DPER1,FAC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION SP(4,4),XJI(2,2),P(2,4),DPER(2,2),PJI(2,4)
	DIMENSION DPER1(3,3)

	DPER = 0.0
	DPER(1,1) = DPER1(1,1)
	DPER(2,2) = DPER1(2,2)
C	CALCULATE HERE J~*P
C	-------------------
	PJI = MATMUL(XJI,P)

	SP = SP + MATMUL(MATMUL(TRANSPOSE(PJI),DPER),PJI)*FAC

	RETURN
      END
C
C	================================================================
C	ASSEMBLY OF THE ALL THE TERMS CALCULATED ABOVE IN LEFT HAND SIDE
C	USE LINER INTERPOLATION IN TIME USING FINITE DIFFERENCE AND
C	USING CRANK-NICOLSON TYPE APPROXIMATION THETA=1/2 IN BOTH EQUATIONS
C	AND TAKE CONSTANT TIME INCREMENT Dt
C	===============================================================
      SUBROUTINE SASSEMB (SS,SM,CC,SP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /SEEP/ NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	DIMENSION SS(20,20),SM(16,16),CC(16,4),SP(4,4)

C	ADD THE CONTRIBUTION FROM SOLID PART TO THE STIFNESS MATRIX
C	-----------------------------------------------------------
	DO I = 1,16
	DO J = 1,16
	SS(I,J) = SM(I,J)
	END DO
	END DO

C	ADD THE COUPLING TERM TO THE STIFNESS MATRIX
C	--------------------------------------------
	DO I = 1,16
	DO J = 1,4
	SS(I,J+16) = CC(I,J)
	SS(J+16,I) = CC(I,J)
	END DO
	END DO

C	ADD THE CONTRIBUTION FROM POREWATER PART TO THE STIFNESS MATRIX
C	----------------------------------------------------------------
	DO I = 1,4
	DO J = 1,4
	SS(I+16,J+16) = -DTINC*SP(I,J)
	END DO
	END DO

	RETURN
	END

C	============================================================
      SUBROUTINE SARRENG (S,SS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATE COUPLING PART OT THE LINEAR STIFFNESS MATRIX
C	-----------------------------------------------------------
C     CC(16,8)   = COUPLING TERM
C     B(3,2*NNO) = LINEAR STRAIN-DISPLACEMENT MATRIX
C     FAC        = DET*WT*TH
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C	H(8)	   = SHAPE FUNCTION
C     ----------------------------------------------------------------
      DIMENSION S(300),SS(20,20),IAR8(20,24),SS1(24,24)

	CALL CLEAR (20,24,IAR8)
	CALL CLEAR (24,24,SS1)

C	8 NODE ELEMENT
C	----------------------------------------
	IAR8(1,1)   = 1
	IAR8(2,2)   = 1
	IAR8(3,4)   = 1
	IAR8(4,5)   = 1
	IAR8(5,7)   = 1
	IAR8(6,8)   = 1
	IAR8(7,10)  = 1
	IAR8(8,11)  = 1
	IAR8(9,13)  = 1
	IAR8(10,14) = 1
	IAR8(11,16) = 1
	IAR8(12,17) = 1
	IAR8(13,19) = 1
	IAR8(14,20) = 1
	IAR8(15,22) = 1
	IAR8(16,23) = 1

	IAR8(17,3)  = 1
	IAR8(18,6)  = 1
	IAR8(19,9)  = 1
	IAR8(20,12) = 1
	
	SS1 = SS1 + MATMUL(MATMUL(TRANSPOSE(IAR8),SS),IAR8)

C	MAKE AND ARRANGE A UPPER TRANGULAR MATRIX ROW WISE
C	--------------------------------------------------
	II = 0
	DO I = 1,24
	DO J = I,24
	II = II + 1
	S(II) = S(II) + SS1(I,J)
	END DO
	END DO

	RETURN
      END
C
C	================================================================
C	FIND THE RESULTANT FORCE VECTOR AND ADD THE
C	CONTRIBUTION FROM THE Nth TIME STEP
C	ASSEMBLY OF THE ALL THE TERMS CALCULATED ABOVE IN RIGHT HAND SIDE
C	USE LINER INTERPOLATION IN TIME USING FINITE DIFFERENCE AND
C	USING CRANK-NICOLSON TYPE APPROXIMATION THETA=1/2 IN BOTH EQUATIONS
C	AND TAKE CONSTANT TIME INCREMENT Dt
C	----------------------------------------------------------------
      SUBROUTINE SNTHRE (SS,SM,CC,SP,RE,EDIS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /SEEP/ NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	DIMENSION SS1(20,20),SM(16,16),CC(16,4),SP(4,4),SS2(24,24)
	DIMENSION RE(24),EDIS(24),SS(20,20),SS3(24,24)

C	INTIALIZE INTERNAL VARIABLES
C	----------------------------
	CALL CLEAR (20,20,SS1)

C	ASSEMBLE COUPLING TERM
C	----------------------
	DO I = 1,16
	DO J = 1,4
	SS1(J+16,I) = CC(I,J)
	END DO
	END DO

C	ARRANGE THE SS MATRIX ACCORDING TO THE ORDER OF 
C	DISPLACEMENT VECTOR(U1,V1,Uw1,......Un,Vn,Uwn...) n=1,NNO
C	---------------------------------------------------------
	CALL SARRE1 (SS2,SS1)

C	CALCULATE THE RESULTANT FORCE VECTORE
C	-------------------------------------	
	re = -matmul(ss2,edis)
c	print*,re

C	CALCULATE THE RESULTANT FORCE VECTORE FROM SPECIFIED PORE PRESSURE
C	------------------------------------------------------------------
C	CALCULATE THE ADITIONAL FORCE DUE TO PRE SPECIFIED PORE PRESUURE
C	------------------------------------------------------------------
	CALL SARRE1 (SS3,SS)
		

	RETURN
	END

C	================================================================
      SUBROUTINE SARRE1 (SS1,SS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATE COUPLING PART OT THE LINEAR STIFFNESS MATRIX
C	-----------------------------------------------------------
C     CC(16,8)   = COUPLING TERM
C     B(3,2*NNO) = LINEAR STRAIN-DISPLACEMENT MATRIX
C     FAC        = DET*WT*TH
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C	H(8)	   = SHAPE FUNCTION
C     ----------------------------------------------------------------
      DIMENSION SS(20,20),IAR8(20,24),SS1(24,24)

	CALL CLEAR (20,24,IAR8)
	CALL CLEAR (24,24,SS1)

C	8 NODE ELEMENT
C	--------------------------------------
	IAR8(1,1)   = 1
	IAR8(2,2)   = 1
	IAR8(3,4)   = 1
	IAR8(4,5)   = 1
	IAR8(5,7)   = 1
	IAR8(6,8)   = 1
	IAR8(7,10)  = 1
	IAR8(8,11)  = 1
	IAR8(9,13)  = 1
	IAR8(10,14) = 1
	IAR8(11,16) = 1
	IAR8(12,17) = 1
	IAR8(13,19) = 1
	IAR8(14,20) = 1
	IAR8(15,22) = 1
	IAR8(16,23) = 1

	IAR8(17,3)  = 1
	IAR8(18,6)  = 1
	IAR8(19,9) = 1
	IAR8(20,12) = 1
	
	SS1 = SS1 + MATMUL(MATMUL(TRANSPOSE(IAR8),SS),IAR8)

	RETURN
      END
C




C	===============================================================
C	2 D - MEMBRAIN ELEMENTS PRESUURE LOADING
C	===============================================================

      SUBROUTINE LOADVE2DM (PROPG,IGSET,XYZ,NODEX,LM,RT,RP,R,
     1                   MGP,MXY,MEX,MEF,NSID,WTABZ)
C	----------------------------------------------------------------
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
C	DC(3)			 = DIRECTION COSINES OF LOAD VECTOR
C	NFCE             = FACE NUMBER OF THE ELEMENT
C	WTABZ            = WATER TABLE IN THE Z DIRECTION
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
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
	COMMON /DENSI/ DENT
C
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),
     +		  NODEX(MEX,1),LM(MEF,1),XJ(9),DC(3,1)
      DIMENSION RT(3,1),RP(NELE),R(1),H(8),P(3,8),XJI(9)
      DIMENSION RL(3,8),NNOD(3)
	DIMENSION NSID(6,NELE),WTABZ(NELE)

C     -----------------------------------
C     ELEMENT LOOP,SKIP UNLOADED ELEMENTS
C     -----------------------------------
      DO 900  IELE=1,NELE
	DO 11 I=1,3
	DC(I,1)=0.
 11	CONTINUE
	INDEL = 0
 110  IF (RP(IELE).NE.0.) INDEL = 1
 120  IF (INDEL.EQ.0) GOTO 900
C

	DO 2000 ISIDE = 1,4
	NSIDE = NSID(ISIDE,IELE)
	IF (NSIDE.EQ.0) GOTO 2000
      CALL CLEARA (RL,24)
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
C	---------------------
C	DEFINE EACH SIDE
C	NSIDE 1: 1 2 OR 1 5 2
C	NSIDE 2: 2 3 OR 2 6 3
C	NSIDE 3: 3 4 OR 3 7 4
C	NSIDE 4: 4 1 OR 4 8 1
C	---------------------
C	FOR 4 NODE
C	FOR 8 NODE
      GOTO (41,42,43,44), NSIDE
 41	NNOD(1) = 1
	NNOD(2) = 2
	NNOD(3) = 5
	GOTO 49
 42	NNOD(1) = 2
	NNOD(2) = 3
	NNOD(3) = 6
	GOTO 49
 43	NNOD(1) = 3
	NNOD(2) = 4
	NNOD(3) = 7
	GOTO 49
 44	NNOD(1) = 4
	NNOD(2) = 1
	NNOD(3) = 8
  49	CONTINUE
C	CALCULATE THE CORDINATE OF THE NODES HAVING PRESSURE LOADS
C	-----------------
C	FOR HYDROSTATIC PRESSURE LOADING NORMAL TO THE SIDE
C	ZREF:	Y-COR ORDINATE OF THE FREE SURFACE
C	AREA:   AREA OF THE ELEMENT INTERACTING WITH WATER
C	DC:     DIRECTION COSINES OF THE NORMAL TO THE SURFACE
C	NSIDE:  SIDE NUMBER INTERACTING WITH WATER
C	DENT:	DENSITY OF WATER
	YREF = WTABZ(IELE)

C	CALL DIR2D(XYZ(1,IEL),NNO,NNOD,YREF,RL)  !RELEASE ERROR SONGSAK MAY2006  IEL --> IELE
	CALL DIR2D(XYZ(1,IELE),NNO,NNOD,YREF,RL)


C	END OF HYDROSTATIC PRESSUR LOADING
	II = 0
 	DO 600  INO=1,NNO
      DO 600  INF=1,NNF
	II = II+1
      IEQ = LM(II,IELE)
 600  IF (IEQ.NE.0) R(IEQ) = R(IEQ) + RL(INF,INO)
 2000 CONTINUE
C	END LOOP OVER LOADED SIDES
 900  CONTINUE
C	END LOOP OVER ELEMENTS
C	NEXT LINE ADDED BY GILSON - JAN2004 (LOAD INPUT)
	CALL LDASEM (R)
C
      RETURN
      END
C
C	==========================================================
C	CALCULATE HERE LENGHT AND DIRECTION COSINES OY OUT WORD
C	NORMAL TO THE SIDE HAVING PRESSURE LOADINGS
C	==========================================================
	SUBROUTINE DIR2D(XYZ,NNO,NNOD,YREF,RL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION XYZ(2,NNO),NNOD(3),AB(3),BC(3),RL(3,8)
	COMMON /DENSI/ DENT

C	FOR 4 NODE ELEMENT
	IF (NNO.EQ.4) THEN
C	CALCULATE THE CO ORDINATES OF THE NODES AT CORESPONDING SIDE
	XCA = XYZ(1,NNOD(1))
	YCA = XYZ(2,NNOD(1))	
	XCB = XYZ(1,NNOD(2))
	YCB = XYZ(2,NNOD(2))
C	CALCULATE LENGTH
	SLENGTH = SQRT((XCB-XCA)**2+(YCB-YCA)**2)
C	CALCULATE DIRECTION COSINES
C	CALCULATE THE VECTOR AB
	AB(1) = XCB-XCA
	AB(2) = YCB-YCA
	AB(3) = 0.0
	DCX = AB(2)/SLENGTH
	DCY = -AB(1)/SLENGTH
C	CALCULATE THE PRESSURE AT THE ENDS
C	CONSIDER Y-IS IN NEGATIVE GRAVITY DIRECTION
	PA = (YCA-YREF)*9.81*DENT
	IF((YCA-YREF).GT.0.0) PA = 0
	PB = (YCB-YREF)*9.81*DENT
	IF((YCB-YREF).GT.0.0) PB = 0
C	CALCULATE THE RESULTANT FORCE AT A,B ALONG NORAMAL TO THE LINE
	IF(ABS(PB).GE.ABS(PA)) THEN
	RA = SLENGTH*(2*PA+PB)/6.0
	RB = SLENGTH*(PA+2*PB)/6.0
	ELSE
	RA = SLENGTH*(2*PB+PA)/6.0
	RB = SLENGTH*(PB+2*PA)/6.0
	END IF
C	CALCULATE THE COMPONENTS OF LOADS AT POINT A ALONG X AND Y DIR
	RAX = RA*DCX
	RAY = RA*DCY
C	CALCULATE THE COMPONENTS OF LOADS AT POINT B ALONG X AND Y DIR
	RBX = RB*DCX
	RBY = RB*DCY
C	ASSEMBLE THIS INTO RL
	RL(1,NNOD(1)) = RAX +RL(1,NNOD(1))
	RL(2,NNOD(1)) = RAY +RL(2,NNOD(1))
	RL(1,NNOD(2)) = RBX +RL(1,NNOD(2))
	RL(2,NNOD(2)) = RBY +RL(2,NNOD(2))

	END IF

C	FOR 8 NODE ELEMENT
	IF (NNO.EQ.8) THEN
C	CALCULATE THE CO ORDINATES OF THE NODES AT CORESPONDING SIDE
	XCA = XYZ(1,NNOD(1))
	YCA = XYZ(2,NNOD(1))	
	XCB = XYZ(1,NNOD(3))
	YCB = XYZ(2,NNOD(3))
	XCC = XYZ(1,NNOD(2))
	YCC = XYZ(2,NNOD(2))
C	CALCULATE THE LENGTH
	SLENGTH1 = SQRT((XCB-XCA)**2+(YCB-YCA)**2)
	SLENGTH2 = SQRT((XCC-XCB)**2+(YCC-YCB)**2)
	SLENGTH	= SLENGTH1 + SLENGTH2
C	CALCULATE THE DIRECTION COSIGNS
C	CALCULATE THE VECTOR AB
	AB(1) = XCB-XCA
	AB(2) = YCB-YCA
	AB(3) = 0.0
	DCX1 = AB(2)/SLENGTH1
	DCY1 = -AB(1)/SLENGTH1
C	CALCULATE THE VECTOR BC
	BC(1) = XCC-XCB
	BC(2) = YCC-YCB
	BC(3) = 0.0
	DCX2 = BC(2)/SLENGTH2
	DCY2 = -BC(1)/SLENGTH2
C	CALCULATE THE PRESSURE AT THE ENDS
C	CONSIDER Y-IS IN NEGATIVE GRAVITY DIRECTION
	PA = (YCA-YREF)*9.81*DENT
	IF((YCA-YREF).GT.0.0) PA = 0.0
	PB = (YCB-YREF)*9.81*DENT
	IF((YCB-YREF).GT.0.0) PB = 0.0
	PC = (YCC-YREF)*9.81*DENT
	IF((YCC-YREF).GT.0.0) PC = 0.0
C	CALCULATE THE RESULTANT FORCE AT A,B ALONG NORAMAL TO THE LINE
	IF(ABS(PB).GE.ABS(PA)) THEN
	RA = SLENGTH1*(2*PA+PB)/6.0
	RB1 = SLENGTH1*(PA+2*PB)/6.0
	ELSE
	RA = SLENGTH1*(2*PB+PA)/6.0
	RB1 = SLENGTH1*(PB+2*PA)/6.0
	END IF


C	CALCULATE THE RESULTANT FORCE AT B,C ALONG NORAMAL TO THE LINE
	IF(ABS(PC).GE.ABS(PB)) THEN
	RB2 = SLENGTH2*(2*PB+PC)/6.0 
	RC = SLENGTH2*(PB+2*PC)/6.0
	ELSE
	RB2 = SLENGTH2*(2*PC+PB)/6.0 
	RC = SLENGTH2*(PC+2*PB)/6.0
	END IF

C	CALCULATE THE COMPONENTS OF LOADS AT POINT A ALONG X AND Y DIR
	RAX = RA*DCX1
	RAY = RA*DCY1
C	CALCULATE THE COMPONENTS OF LOADS AT POINT B1 ALONG X AND Y DIR
	RBX = RB1*DCX1
	RBY = RB1*DCY1

C	CALCULATE THE COMPONENTS OF LOADS AT POINT B2 ALONG X AND Y DIR
	RBX = RB2*DCX2 + RBX
	RBY = RB2*DCY2 + RBY
C	CALCULATE THE COMPONENTS OF LOADS AT POINT C ALONG X AND Y DIR
	RCX = RC*DCX2
	RCY = RC*DCY2
C	ASSEMBLE THIS INTO RL
	RL(1,NNOD(1)) = RAX +RL(1,NNOD(1))
	RL(2,NNOD(1)) = RAY +RL(2,NNOD(1))
	RL(1,NNOD(2)) = RCX +RL(1,NNOD(2))
	RL(2,NNOD(2)) = RCY +RL(2,NNOD(2))
	RL(1,NNOD(3)) = RBX +RL(1,NNOD(3))
	RL(2,NNOD(3)) = RBY +RL(2,NNOD(3))

	END IF

	RETURN
	END
C
C	===============================================================
C	===============================================================
C	===============================================================

