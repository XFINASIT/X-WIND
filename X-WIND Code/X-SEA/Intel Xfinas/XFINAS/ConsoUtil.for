
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE EXTRACTDISP38(EDIS,EDISU)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C
C	-----------------------------------
C	SUBROUTINE USING FOR RE-ORDER FOR 
C	DISPACEMENT VECTOR IN SOLID ELEMENT
C	-----------------------------------
C
	DIMENSION EDIS(1),EDISU(1)

	EDISU(1:24) = 0.0D0

C	LOCAL NODE-1
	EDISU(1)  = EDIS(1)
	EDISU(2)  = EDIS(2)
	EDISU(3)  = EDIS(3)
C	LOCAL NODE-2
	EDISU(4)  = EDIS(5)
	EDISU(5)  = EDIS(6)
	EDISU(6)  = EDIS(7)
C	LOCAL NODE-3
	EDISU(7)  = EDIS(9)
	EDISU(8)  = EDIS(10)
	EDISU(9)  = EDIS(11)
C	LOCAL NODE-4
	EDISU(10) = EDIS(13)
	EDISU(11) = EDIS(14)
	EDISU(12) = EDIS(15)
C	LOCAL NODE-5
	EDISU(13) = EDIS(17)
	EDISU(14) = EDIS(18)
	EDISU(15) = EDIS(19)
C	LOCAL NODE-6
	EDISU(16) = EDIS(21)
	EDISU(17) = EDIS(22)
	EDISU(18) = EDIS(23)
C	LOCAL NODE-7
	EDISU(19) = EDIS(25)
	EDISU(20) = EDIS(26)
	EDISU(21) = EDIS(27)
C	LOCAL NODE-8
	EDISU(22) = EDIS(29)
	EDISU(23) = EDIS(30)
	EDISU(24) = EDIS(31)

	RETURN
	END

C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE EXTRACTPRES38(EDIS,EDISP)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	---------------------------------
C	SUBROUTINE USING FOR RE-ORDER FOR 
C	PRESSURE VECTOR IN SOLID ELEMENT
C	--------------------------------

	DIMENSION EDIS(1),EDISP(1)

	EDISP(1:8) = 0.0D0

C	LOCAL NODE-1
	EDISP(1)  = EDIS(4)
C	LOCAL NODE-2
	EDISP(2)  = EDIS(8)
C	LOCAL NODE-3
	EDISP(3)  = EDIS(12)
C	LOCAL NODE-4
	EDISP(4)  = EDIS(16)
C	LOCAL NODE-5
	EDISP(5)  = EDIS(20)
C	LOCAL NODE-6
	EDISP(6)  = EDIS(24)
C	LOCAL NODE-7
	EDISP(7)  = EDIS(28)
C	LOCAL NODE-8
	EDISP(8)  = EDIS(32)

	RETURN
	END
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE PERMEA(PROPM,PM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C
C	------------------------------------------------------------
C	SETTING BASIC PROPERTIES OF PERMEABILITY MATRIX PER UNIT WT.
C	------------------------------------------------------------
C
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
C
	DIMENSION PROPM(NMP),PM(3,3)
C
	PM(1:3,1:3) = 0.0D0
C
C	INSERTING VALUE	
C
	PM(1,1) = PROPM(6)	
	PM(2,2) = PROPM(7)	
	PM(3,3) = PROPM(8)	
C
	RETURN
	END

C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE BMSLD(P,XJI,BM,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     ----------------------------------------------------------------
C     EVALUATES THE GLOBAL LINEAR STRAIN-DISPLACEMENT MATRIX FOR A
C     CURVILINEAR ISOPARAMETRIC HEXAHEDRON (8 NODES)
C	----------------------------------------------------
C     INPUT OUTPUT VARIABLES
C	----------------------
C     P(3,NNO)    = SHAPE FUNCTION DERIVATIVES WITH RESPECT TO R,S,T
C     XJI(3,3)    = INVERSE OF THE JACOBIAN MATRIX
C     BM(6,24)   = COMPRESSED LINEAR STRAIN DISPLACEMENT MATRIX
C     NNO         = NUMBER OF NODES USED TO DESCRIBE THIS ELEMENT
C     ----------------------------------------------------------------

      DIMENSION P(3,8),XJI(3,3),BM(6,24),B(24)
C
      CALL CLEARA (B,24)

C	------------------------------------------------
C	 COMPUTE COMPACTED STRESS-STRAIN GRADIENT VECTOR
C	------------------------------------------------
C
      K = 1
C
      DO INO=1,NNO
      DO I=1,3
      B(K)   = B(K)   + XJI(1,I)*P(I,INO)
      B(K+1) = B(K+1) + XJI(2,I)*P(I,INO)
      B(K+2) = B(K+2) + XJI(3,I)*P(I,INO)
	END DO
      K = K+3
	END DO
C
C	----------------------------------------
C	 COMPUTE COMPACTED STRESS-STRAIN MATRIX
C	----------------------------------------

	BM(1:6,1:24) = 0.0D0
	DO INO = 1,8
	HX =0.0D0
	HY =0.0D0
	HZ =0.0D0
	DO I=1,3
      HX  = HX  + XJI(1,I)*P(I,INO)
      HY  = HY  + XJI(2,I)*P(I,INO)
      HZ  = HZ  + XJI(3,I)*P(I,INO)
	END DO
C
C	--------------------------
C  	 COMPUTE MATRIX [B](6,24)
C	--------------------------
C
C	NORMAL STRAIN
C
	BM(1,3*INO-2) = HX
	BM(2,3*INO-1) = HY
	BM(3,3*INO-0) = HZ
C
C	SHEAR STRAIN 12
C
	BM(4,3*INO-2) = HY
	BM(4,3*INO-1) = HX
C
C	SHEAR STRAIN 13
C
	BM(5,3*INO-2) = HZ
	BM(5,3*INO-0) = HX
C
C	SHEAR STRAIN 23
C
	BM(6,3*INO-1) = HZ
	BM(6,3*INO-0) = HY
C
	END DO

      RETURN
      END
C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE CSOLID(C)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C	--------------------------------------------
C	COMPUTE LINEAR ISOTROPIC CONSTITUTIVE MATRIX
C	--------------------------------------------
C   
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
	DIMENSION C(6,6)
C      
	DO I =1,6
	DO J =1,6
      	C(I,J)=0.0D0
     	END DO
	END DO
C
C	LINEAR CONSITITUTIVE 
C
      C(1,1)  = A1
      C(1,2)  = B1
      C(1,3)  = B1
      C(2,1)  = B1
      C(2,2)  = A1
      C(2,3)  = B1
      C(3,1)  = B1
      C(3,2)  = B1
      C(3,3)  = A1

      C(4,4)  = C1
      C(5,5)  = C1
      C(6,6)  = C1

      RETURN
      END
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE GRADP(BM,QM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C
	DIMENSION BM(6,24),QM(3,8)
C
	DO I=1,8
	QM(1,I) = BM(1,3*I-2) !(1, 1...4...7...10...13...16...19...22)
	QM(2,I) = BM(2,3*I-1) !(2, 2...5...8...11...14...17...20...23)
	QM(3,I) = BM(3,3*I-0) !(3, 3...6...9...12...15...18...21...24)
	END DO
C
	RETURN
	END
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE SOLIDGE(R,S,T,DETJO,DETJ,TT,MM,GE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	-----------------------------------------------
C	SUBROUTINE COMPUTE THE COMPATIBLE ENHANCED TERM
C	-----------------------------------------------

	DIMENSION GM(6,MM),GE(6,MM),TT(6,6),EM(6,MM)

	DO I =1,6
	 DO J=1,MM
		GM(I,J) = 0.
	    GE(I,J) = 0.
		EM(I,J) = 0.
	 END DO
	END DO

	GM(1,1) =  -2.0*R
	GM(2,2) =  -2.0*S
	GM(3,3) =  -2.0*T
	GM(4,4) =  -2.0*R
	GM(4,5) =  -2.0*S

	IF(MM.EQ.5)GOTO 100

	GM(5,6) =  -2.0*R
	GM(5,7) =  -2.0*T
	GM(6,8) =  -2.0*S
	GM(6,9) =  -2.0*T


	IF(MM.EQ.9)GOTO 100

	GM(4,10) = -3.0*R*S
	GM(4,11) = -3.0*S*T
	GM(4,12) = -3.0*R*T

	GM(5,13) = -3.0*R*S
	GM(5,14) = -3.0*S*T
	GM(5,15) = -3.0*R*T

	GM(6,16) = -3.0*R*S
	GM(6,17) = -3.0*S*T
	GM(6,18) = -3.0*R*T

	GM(1,19) = -3.0*R*T
	GM(1,20) = -3.0*R*S

	GM(2,21) = -3.0*R*S
	GM(2,22) = -3.0*S*T

	GM(3,23) = -3.0*S*T
	GM(3,24) = -3.0*R*T

	IF(MM.EQ.24)GOTO 100

100	CONTINUE

	EM = MATMUL(TT,GM)
      CONST = (DETJO/DETJ)

C	COMPUTE [GE] MATRIX

	DO IEM = 1,6
		DO  JEM = 1,MM
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
		END DO
	END DO

	RETURN
	END
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE EQUIFORCE38(RES,RPR,RE)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	
C	---------------------------------------------------------------------------
C	SUBROUTINE USED FOR RE-EQUILIBRIUM LOAD FROM SOLID-PRESSURE TO ELEMENT LOAD
C	---------------------------------------------------------------------------
C
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
C	
      COMMON  /EASP/  MB,MG,MQ
C
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
	DIMENSION RES(MB),RPR(MQ),RE(NEF)
C
	RE(1:32) = 0.0D0
C	
	DO I=1,8
	RE(4*I-3) = RES(3*I-2)
	RE(4*I-2) = RES(3*I-1)
	RE(4*I-1) = RES(3*I-0)
	RE(4*I-0) = RPR(I)
	END DO
C
	RETURN
	END
C	=================================================================
C	=================================================================
C	=================================================================

	SUBROUTINE DISVECTOR(ID,ISN,MSF,INEQ,A,DIS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	=============================================================
C	SUBROUTINE TAKES THE PORE PRESSURE OUT AND STORE TO ARRAY:DIS
C	PURPOSE
C	COMPUTE : EUCLIDIAN NORM VECTOR
C	ALLOW ONLY DISPLACEMENT DEGREE OF FREEDOM
C	=========================================
C	DESCRIPTION
C	===========

C     IEQ   = ID(LIMEQ(2),LIMEQ(1))
C	ID(1,2) : 1 MEANS LIMEQ(2) - DIRECTION 
C					  LIMEQ(2) = 1 ; DISPLACEMENT IN X-DIRECTION
C					  LIMEQ(2) = 2 ; DISPLACEMENT IN Y-DIRECTION
C					  LIMEQ(2) = 3 ; DISPLACEMENT IN Z-DIRECTION
C			: 2 MEANS LIMEQ(1) NODE NUMBER

C	======================================

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT

      DIMENSION ID(MSF,1),DIS(INEQ),A(INEQ)

	DIS = 0.0

C	========================================
C	DISPLACEMENT IN X-DIRECTION LIMEQ(2) = 1
C	========================================

	IF(IDOF(1).EQ.1)THEN
	DO I=1,ISN
	IEQ1 = ID(1,I)
	IF(ID(1,I).NE.0)THEN
	DIS(IEQ1) = A(IEQ1)
	END IF
	END DO
	END IF

C	========================================
C	DISPLACEMENT IN Y-DIRECTION LIMEQ(2) = 2
C	========================================

	IF(IDOF(2).EQ.2)THEN
	DO I=1,ISN
	IEQ2 = ID(2,I)
	IF(ID(2,I).NE.0)THEN
	DIS(IEQ2) = A(IEQ2)
	END IF
	END DO
	END IF

	IF(IDOF(3).EQ.9)RETURN

C	========================================
C	DISPLACEMENT IN Z-DIRECTION LIMEQ(2) = 3
C	THIS BLOCK CONDITION FOR 3-D 
C	============================

	IF(IDOF(3).EQ.3.AND.IDOF(4).EQ.9)THEN
	DO I=1,ISN
	IEQ3 = ID(3,I)
	IF(ID(3,I).NE.0)THEN
	DIS(IEQ3) = A(IEQ3)
	END IF
	END DO
	END IF

	RETURN
	END

C	======================================================================
C	======================================================================
C	======================================================================

	SUBROUTINE EXTRACTDISP3D(EDIST,EDIS1)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4(I,N)

	DIMENSION EDIST(32),EDIS1(24)

	EDIS1 = 0.0

	EDIS1(1)  = EDIST(1)
	EDIS1(2)  = EDIST(2)
	EDIS1(3)  = EDIST(3)

	EDIS1(4)  = EDIST(5)
	EDIS1(5)  = EDIST(6)
	EDIS1(6)  = EDIST(7)

	EDIS1(7)  = EDIST(9)
	EDIS1(8)  = EDIST(10)
	EDIS1(9)  = EDIST(11)

	EDIS1(10) = EDIST(13)
	EDIS1(11) = EDIST(14)
	EDIS1(12) = EDIST(15)

	EDIS1(13) = EDIST(17)
	EDIS1(14) = EDIST(18)
	EDIS1(15) = EDIST(19)

	EDIS1(16) = EDIST(21)
	EDIS1(17) = EDIST(22)
	EDIS1(18) = EDIST(23)

	EDIS1(19) = EDIST(25)
	EDIS1(20) = EDIST(26)
	EDIS1(21) = EDIST(27)

	EDIS1(22) = EDIST(29)
	EDIS1(23) = EDIST(30)
	EDIS1(24) = EDIST(31)

	RETURN
	END

C	======================================================================
C	======================================================================
C	======================================================================

	SUBROUTINE EXTRACTPRES3D(EDIST,EDIS2)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4(I,N)

	DIMENSION EDIST(32),EDIS2(8)

	EDIS2 = 0.0

	EDIS2(1) = EDIST(4)
	EDIS2(2) = EDIST(8)
	EDIS2(3) = EDIST(12)
	EDIS2(4) = EDIST(16)
	EDIS2(5) = EDIST(20)
	EDIS2(6) = EDIST(24)
	EDIS2(7) = EDIST(28)
	EDIS2(8) = EDIST(32)


	RETURN
	END

C

C	======================================================================
C	======================================================================
C	======================================================================
	SUBROUTINE EXTRACTDISP(EDIST,EDIS1)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4(I,N)
C	======================================================
C	THIS SUBROUTINE USING TO EXTRACT DISPLACEMENT DOF FROM
C	TOTAL ELEMENT DISPLACEMENT VECTOR APPEAR IN SUBROUTINE ELLOOP
C	=============================================================

C	   VARIABLE			DESCRIPTION
C	   ========			===========
C	   EDIS(12)			FOR FOUR NODES MEMBRANE ELEMENT (3 DOF/NODE)
C						8 DISPLACEMENT VECTORS
C						4 PRESSURE VECTORS

C	   EDIS1(8)			8 DISPLACEMENT VECTOR COMPONENTS

C	==================================================== 

	DIMENSION EDIST(12),EDIS1(8)

	EDIS1 = 0.0

	EDIS1(1) = EDIST(1)
	EDIS1(2) = EDIST(2)

	EDIS1(3) = EDIST(4)
	EDIS1(4) = EDIST(5)
	
	EDIS1(5) = EDIST(7)
	EDIS1(6) = EDIST(8)
	
	EDIS1(7) = EDIST(10)
	EDIS1(8) = EDIST(11)

	RETURN
	END

C	======================================================================
C	======================================================================
C	======================================================================
	SUBROUTINE EXTRACTPRES(EDIST,EDIS2)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4(I,N)
C	==================================================
C	THIS SUBROUTINE USING TO EXTRACT PRESSURE DOF FROM
C	TOTAL ELEMENT DISPLACEMENT VECTOR APPEAR IN SUBROUTINE ELLOOP
C	=============================================================

C	  VARIABLE			DESCRIPTION
C	  ========			===========
C	   EDIS(12)			FOR FOUR NODES MEMBRANE ELEMENT (3 DOF/NODE)
C						8 DISPLACEMENT VECTORS
C						4 PRESSURE VECTORS

C	   EDIS2(4)			4 PRESSURE VECTOR COMPONENTS

C	================================================ 

	DIMENSION EDIST(12),EDIS2(4)

	EDIS2 = 0.0

	EDIS2(1) = EDIST(3)
	EDIS2(2) = EDIST(6)

	EDIS2(3) = EDIST(9)
	EDIS2(4) = EDIST(12)

	RETURN
	END

C	======================================================================
C	======================================================================
C	======================================================================
      SUBROUTINE KNLSTIFF (S,SIG,B)
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

      DIMENSION S(1),SIG(6),B(63)
	DIMENSION Q(9,9),G(9,24),SNL(24,24)


C	...SETTING MATRIX [G]

	G   =0.
	Q   = 0.
	SNL = 0.

	DO I =1,8

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
	DO IPT =1,24
	 DO JPT =IPT,24 
	     K = K+1
	  S(K) = S(K)+SNL(IPT,JPT)
	END DO
	END DO


	RETURN 
	END

C	======================================================================
C	======================================================================
C	======================================================================

	SUBROUTINE UPPERTRIA(S,SS)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SS(32,32),SS1(32,32),S(528),IAR4(32,32)

c	----------------------------------
C	Initialized Upper Triangular Array
c	----------------------------------

	S = 0.0
	SS1 = 0.0
	IAR4 = 0.0

	IAR4(1,1)    = 1
	IAR4(2,2)    = 1
	IAR4(3,3)    = 1
	IAR4(4,5)    = 1
	IAR4(5,6)    = 1
	IAR4(6,7)    = 1
	IAR4(7,9)    = 1
	IAR4(8,10)   = 1
	IAR4(9,11)   = 1
	IAR4(10,13)  = 1
	IAR4(11,14)  = 1
	IAR4(12,15)  = 1
	IAR4(13,17)  = 1
	IAR4(14,18)  = 1
	IAR4(15,19)  = 1
	IAR4(16,21)  = 1
	IAR4(17,22)  = 1
	IAR4(18,23)  = 1
	IAR4(19,25)  = 1
	IAR4(20,26)  = 1
	IAR4(21,27)  = 1
	IAR4(22,29)  = 1
	IAR4(23,30)  = 1
	IAR4(24,31)  = 1

	IAR4(25,4)   = 1
	IAR4(26,8)   = 1
	IAR4(27,12)  = 1
	IAR4(28,16)  = 1
	IAR4(29,20)  = 1
	IAR4(30,24)  = 1
	IAR4(31,28)  = 1
	IAR4(32,32)  = 1

	SS1 = SS1 + MATMUL(MATMUL(TRANSPOSE(IAR4),SS),IAR4)

C	MAKE AND ARRANGE A UPPER TRANGULAR MATRIX ROW WISE
C	--------------------------------------------------


C	DO I=1,32
C	DO J=I,32
C	WRITE(*,*)I,J,SS1(I,J)
C	ENDDO
C	PAUSE
C	END DO



	II = 0
	DO I = 1,32
	DO J = I,32
	II = II + 1
	S(II) = S(II) + SS1(I,J)
	END DO
	END DO


	RETURN
	END
C	======================================================================
C	======================================================================
C	======================================================================

	SUBROUTINE UPPERTRIA2D(S,SS,NNO,NEF)
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SS(NEF,NEF),SS1(NEF,NEF),S(1),AR(NEF,NEF)

c	----------------------------------
C	Initialized Upper Triangular Array
c	----------------------------------
	AR(1:NEF,1:NEF) = 0.0D0

	DO I = 1,NNO
	AR(2*I-1,3*I-2) = 1
	AR(2*I-0,3*I-1) = 1

	AR(I+2*NNO,3*I-0) = 1
	ENDDO

	SS1 = MATMUL(MATMUL(TRANSPOSE(AR),SS),AR)

	II = 0
	DO I = 1,NEF
	DO J = I,NEF
	II = II + 1
	S(II) = SS1(I,J)
	END DO
	END DO


	RETURN
	END
C	======================================================================
C	======================================================================
C	======================================================================