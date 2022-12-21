C=====================================================================
      SUBROUTINE DELA3D (D,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------
C     GENERATES THE THREE-DIMENSIONAL STRESS-STRAIN MATRIX FOR A
C     LINEAR ELASTIC,ISOTROPIC MATERIAL
C	---------------------------------
C     D(6,6)  = ELASTIC,ISOTROPIC STRESS-STRAIN MATRIX
C     IND     = FLAG (AXISYMMETRIC=0,PL.STRAIN=1,PL.STRESS=2)
C     ----------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION D(36)
C
      CALL CLEARA (D,36)
      IF (IND-1)  100,100,200
C     ----------------------
C     PLANE STRAIN CONDITION
C     ----------------------
 100  D(1)  = A1
      D(8)  = A1
      D(36) = A1
      D(2)  = B1
      D(6)  = B1
      D(7)  = B1
      D(12) = B1
      D(31) = B1
      D(32) = B1
      D(15) = C1
      D(22) = C1
      D(29) = C1
      RETURN
C     ----------------------
C     PLANE STRESS CONDITION
C     ----------------------
 200  D(1)  = A1
      D(8)  = A1
      D(2)  = B1
      D(7)  = B1
      D(15) = C1
      D(22) = D1
      D(29) = D1
      D(6)  = B2*C2
      D(12) = B2*C2
      D(31) = D(6)
      D(32) = D(6)
      D(36) = A2*C2
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MEBMAT (XY,H,P,XJI,B,ISTYP,XBAR,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------
C     EVALUATES THE GLOBAL LINEAR STRAIN-DISPLACEMENT MATRIX B
C     FOR A QUADRILATERAL MEMBRANE ELEMENT
C	------------------------------------
C     INPUT VARIABLES
C	---------------
C     XY(2,NNO)  =  COORDINATES FOR ELEMENT NODES
C     H(NNO)     =  SHAPE FUNCTIONS
C     P(2,NNO)   =  SHAPE FUNCTION DERIVATIVES WITH RESP.TO R,S
C     XJI(2,2)   =  INVERSE OF THE JACOBIAN MATRIX
C     B(4,NNO)   =  STRAIN-DISPLACEMENT MATRIX
C     ISTYP      =  ELEMENT SUBTYPE
C     XBAR       =  RADIUS AT GAUSS POINT
C     NNO        =  NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ---------------------------------------------------------
      DIMENSION XY(2,1),H(8),P(2,8),XJI(4),B(64)
      CALL CLEARA (B,64)
C
      K = 1
      DO 100  INO=1,NNO
      HIX = XJI(1)*P(1,INO) + XJI(3)*P(2,INO)
      HIY = XJI(2)*P(1,INO) + XJI(4)*P(2,INO)
C
      B(K)   = HIX
      B(K+2) = HIY
      B(K+5) = HIY
      B(K+6) = HIX
 100  K = K+8
      IF (ISTYP.NE.0)  RETURN
C     ---------------------------------------------
C     CALCULATE THE RADIUS AT THE GAUSS POINT
C     FOR ZERO RADIUS EQUATE HOOP TO RADIAL STRAINS
C     ---------------------------------------------
      XBAR = 0.0
      DO 200  INO=1,NNO
 200  XBAR = XBAR + H(INO)*XY(1,INO)
      IF (ABS(XBAR).GT.0.0000001)  GOTO 300
      KB = 4
      DO 250  INO=1,NNO
      B(KB) = B(KB-3)
 250  KB = KB+8
      RETURN
C     ---------------
C     NON-ZERO RADIUS
C     ---------------
 300  DUM = 1./XBAR
      KB  = 4
      DO 350  INO=1,NNO
      B(KB) = H(INO)*DUM
 350  KB = KB+8
C
C
      RETURN
      END
C
C=====================================================================

      SUBROUTINE MEK0PL (S,D,B,FAC,NEF,ISTYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ADDS LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX S(NWS)
C	-----------------------------------------------------------
C     S(136)     = ELEMENT STIFFNESS MATRIX IN UPPER TRIANGULAR FORM
C     D(6,6)     = STRESS - STRAIN MATRIX
C     B(3,2*NNO) = LINEAR STRAIN-DISPLACEMENT MATRIX
C     FAC        = DET*WT*TH
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ISTYP      = ELEMENT SUPTYPE (0=AXISYM.,1/2=PLANE STRAIN/STRESS)
C     ----------------------------------------------------------------
      DIMENSION S(136),D(6,6),B(4,16),DB(4)
C     -----------------------------
C     PLANE STRAIN AND PLANE STRESS
C     -----------------------------
      IF (ISTYP.EQ.0)  GOTO 300
      KL = 1
      DO 100  J=1,NEF,2
      DO 50  K=1,3
 50   DB(K) = (D(K,1)*B(1,J) + D(K,3)*B(3,J)) * FAC
      DO 80   I=J,NEF,2
      S(KL) = S(KL) + B(1,I)*DB(1) + B(3,I)*DB(3)
      KL = KL+1
      S(KL) = S(KL) + B(2,I+1)*DB(2) + B(3,I+1)*DB(3)
 80   KL = KL+1
 100  KL = KL+NEF-J
C
      KS = NEF+1
      DO 200  J=2,NEF,2
      DO 110  K=1,3
 110  DB(K) = (D(K,2)*B(2,J) + D(K,3)*B(3,J)) * FAC
      S(KS) = S(KS) + B(2,J)*DB(2) + B(3,J)*DB(3)
      K  = J+1
      KS = KS+1
      IF (NEF-K)  210,210,150
 150  DO 180  I=K,NEF,2
      S(KS) = S(KS) + B(1,I)*DB(1) + B(3,I)*DB(3)
      KS = KS+1
      S(KS) = S(KS) + B(2,I+1)*DB(2) + B(3,I+1)*DB(3)
 180  KS = KS+1
 200  KS = KS+NEF-J
 210  RETURN
C     ------------
C     AXISYMMETRIC
C     ------------
  300 KL = 1
      DO 101  J=1,NEF,2
      DO 51  K=1,4
      M = K
      IF (K.GT.3) M=6
 51   DB(K) = (D(M,1)*B(1,J) + D(M,3)*B(3,J) + D(M,6)*B(4,J)) * FAC
      DO 81   I=J,NEF,2
      S(KL) = S(KL) + B(1,I)*DB(1) + B(3,I)*DB(3) + B(4,I)*DB(4)
      KL = KL+1
      S(KL) = S(KL) + B(2,I+1)*DB(2) + B(3,I+1)*DB(3)
 81   KL = KL+1
 101  KL = KL+NEF-J
C
      KS = NEF+1
      DO 201  J=2,NEF,2
      DO 111  K=1,4
      M = K
      IF (K.GT.3) M=6
 111  DB(K) = (D(M,2)*B(2,J) + D(M,3)*B(3,J)) * FAC
      S(KS) = S(KS) + B(2,J)*DB(2) + B(3,J)*DB(3)
      K  = J+1
      KS = KS+1
      IF (NEF-K) 410,410,151
 151  DO 181  I=K,NEF,2
      S(KS) = S(KS) + B(1,I)*DB(1) + B(3,I)*DB(3) +B(4,I)*DB(4)
      KS = KS+1
      S(KS) = S(KS) + B(2,I+1)*DB(2) + B(3,I+1)*DB(3)
 181  KS = KS+1
 201  KS = KS+NEF-J
C
 410  RETURN
      END
C
C=====================================================================
      SUBROUTINE JACO2D (MEL,NNO,XY,P,XJI,DET,MCO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(MCO,1),P(2,8),XJ(2,2),XJI(4)
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
C=====================================================================
      SUBROUTINE MEMSIG (STRAIN,STRESS,ISTYP)
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

      IF (ISTYP-1)  100,200,300
C     ------------
C     AXISYMMETRIC
C     ------------
 100  STRESS(1) = STRESS(1) + B1*STRAIN(4)
      STRESS(2) = STRESS(2) + B1*STRAIN(4)
      STRESS(4) = B1*(STRAIN(1)+STRAIN(2)) + A1*STRAIN(4)
      RETURN
C     ------------
C     PLANE STRAIN
C     ------------
 200  STRESS(4) = B1*(STRAIN(1)+STRAIN(2))
      STRAIN(4) = 0.
      RETURN
C     ------------
C     PLANE STRESS
C     ------------
 300  STRESS(4) = 0.
      STRAIN(4) = D2*(STRAIN(1)+STRAIN(2))
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MEKSIG (S,TAU,B,NEF,ISTYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ADDS NONLINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX S(NWS)
C	--------------------------------------------------------------
C     S(136)     = ELEMENT STIFFNESS MATRIX (UPPER TRIANG.ROW-WISE)
C     TAU(4)     = STRESSES TIMES WT*DET*THICK
C     B(3,2*NNO) = STRAIN-DISPLACEMENT MATRIX
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ISTYP      = ELEMENT SUBTYPE (0=AXISYM.,1/2=PLANE STRAIN/STRESS)
C     ----------------------------------------------------------------
      DIMENSION S(136),TAU(4),B(4,16),TB(3)
C
      KL = 1
      DO 200  J=1,NEF,2
      TB(1) = TAU(1)*B(1,J) + TAU(3)*B(3,J)
      TB(2) = TAU(3)*B(1,J) + TAU(2)*B(3,J)
C
      KS = KL
      DO 100  I=J,NEF,2
      KSS = KS+NEF-J+1
      DUM = B(1,I)*TB(1) + B(3,I)*TB(2)
      S(KS)  = S(KS) + DUM
      S(KSS) = S(KSS)+ DUM
 100  KS = KS+2
 200  KL = KL + 2*NEF - 2*J + 1
C     --------------------------------
C     CONTRIBUTION FROM HOOP DIRECTION
C     --------------------------------
      IF (ISTYP.NE.0)  RETURN
      KL = 1
      DO 400  J=1,NEF,2
      TB(3) = TAU(4)*B(4,J)
      DO 300  I=J,NEF,2
      S(KL) = S(KL) + TB(3)*B(4,I)
 300  KL = KL+2
 400  KL = KL + NEF - J
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE VMISES (SIG,EPS,IPEL,STRAIN,STRESS,DP,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     VON MISES YIELD CRITARION
C	CONVERTS ELASTO-PLASTIC STRAIN INCREMENTS INTO STRESS INCREMENTS
C	----------------------------------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     SIG(4)    = STRESSES AT THE END OF THE PREVIOUS UPDATE
C     EPS(4)    = STRAINS AT THE END OF THE PREVIOUS UPDATE
C     YIELD     = YIELD STRESS AT END OF PREVIOUS UPDATE
C     IPEL      = PLASTICITY FLAG (1=ELASTIC,2=PLASTIC)
C     STRAIN(4) = CURRENT TOTAL STRAIN
C     STRESS(4) = RETURNED TOTAL STRESS
C     DP(6,6)   = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     IND       = FLAG (AXISYM.=0,PLANE STRAIN=1,PLANE STRESS=2,3)
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DELEPS(4) = INCREMENT IN STRAIN
C     DELSIG(4) = INCREMENT IN STRESSES, ASSUMING ELASTIC BEHAVIOR
C     DEPS(4)   = SUBINCREMENT IN STRAINS (EQUIVALENCED WITH DELEPS)
C     RATIO     = PART OF STRESS TAKEN ELASTICLY
C     PLAMDA    = PLASTIC STRAIN RATE MULTIPLIER
C     SX,SY,SS, = DEVIATORIC STRESSES
C     FT        = VON MISES YIELD FUNCTION
C	------------------------------
C     VARIABLES IN /ELEM/ AND /HOOK/
C	------------------------------
C     NSINC     = SPECIFIED FIXED NUMBER OF SUBINCREMENTS
C     ITOLEY    = TOLERANCE ON EXCESSIVE YIELD FUNCTION
C     YLD       = CURRENT UPDATED YIELD STRESS
C     ISR       = NUMBER OF STRAIN COMPONENTS
C     IST       = NUMBER OF STRESS COMPONENTS
C     ----------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION SIG(4),EPS(4),STRAIN(4),STRESS(4),DP(6,6)
      DIMENSION DELEPS(4),DELSIG(4),DEPS(4)
      EQUIVALENCE (DELEPS(1),DEPS(1))
C     --------------------------------------------------------------
C     1. CALCULATE INCREMENTAL STRAINS
C     2. CALCULATE INCREMENTAL STRESSES,ASSUMING AN ELASTIC BEHAVIOR
C     3. CALCULATE TOTAL STRESSES,ASSUMING AN ELASTIC BEHAVIOR
C     --------------------------------------------------------------
C	NEXT LINE CHANGED BY GILSON - JUNE2003 (PLASTICITY)
C	DFLOATI(ITOLEY)=1.E-10
	DFLOATI(ITOLEY)=1.E10
	DFF = DFLOATI(ITOLEY)
      DO 110  I=1,ISR
 110  DELEPS(I) = STRAIN(I)-EPS(I)
      CALL MEMSIG (DELEPS,DELSIG,IND)
      STRESS(4) = 0.0D0
      DO 300  I=1,IST
 300  STRESS(I) = SIG(I) + DELSIG(I)
C     ----------------------------------------------------------------
C     4. CHECK WHETHER STATE OF STRESS FALLS OUTSIDE THE YIELD SURFACE
C     ----------------------------------------------------------------
      SM = (STRESS(1) + STRESS(2) + STRESS(4)) / 3.
      SX = STRESS(1)-SM
      SY = STRESS(2)-SM
      SS = STRESS(3)
      SZ = STRESS(4)-SM
      FT = 0.5*(SX*SX + SY*SY + SZ*SZ) + SS*SS - YLD*YLD/3.
      IF (FT)  410,410,450
C     -------------------------------------------------------
C     STATE OF STRESS WITHIN YIELD SURFACE - ELASTIC BEHAVIOR
C     -------------------------------------------------------
 410  IPEL = 1
      IF (IND.EQ.2)  STRAIN(4) = EPS(4) + DELEPS(4)
      GOTO 700
C     ------------------------------------------------------------
C     STATE OF STRESS OUTSIDE YIELD SURFACE - PLASTIC BEHAVIOR
C     DETERMINE PART OF STRAIN TAKEN ELASTICLY AND ADD STRESSES
C     DUE TO THE ELASTIC STRAIN INCREMENT TO THE PREVIOUS STRESSES
C     ------------------------------------------------------------
 450  IPEL = 2
      DO 460  I=1,IST
 460  STRESS(I) = SIG(I)
C
      SM = (STRESS(1)+ STRESS(2) + STRESS(4)) / 3.
      SX = STRESS(1)-SM
      SY = STRESS(2)-SM
      SS = STRESS(3)
      SZ = STRESS(4)-SM
      DM = (DELSIG(1) + DELSIG(2) + DELSIG(4)) / 3.
      DX = DELSIG(1)-DM
      DY = DELSIG(2)-DM
      DS = DELSIG(3)
      DZ = DELSIG(4)-DM
C
      A = DX*DX + DY*DY + 2.*DS*DS + DZ*DZ
      B = SX*DX + SY*DY + 2.*SS*DS + SZ*DZ
      E = SX*SX + SY*SY + 2.*SS*SS + SZ*SZ - 2.*YLD*YLD/3.
	SQQ = B*B-A*E
	SQA = ABS(SQQ)
      RATIO = (-B + SQRT(SQA)) / A
C
      DO 480  I=1,IST
 480  STRESS(I) = SIG(I) + RATIO*DELSIG(I)
      IF (IND.EQ.2)  STRAIN(4) = EPS(4) + RATIO*DELEPS(4)
C     --------------------------------------------------
C     5. DETERMINE NUMBER AND MAGNITUDE OF SUBINCREMENTS
C     --------------------------------------------------
 500	FFTT = SQRT(FT)/YLD
	NFF = INT(FFTT)
	NSINC = 20*NFF + 1
      NSINC = 20.*SQRT(FT)/YLD + 1
      IF (NSINC.GT.30)  NSINC = 30
      FACT = (1.0-RATIO) / NSINC
      DO 510  I=1,ISR
 510  DEPS(I) = FACT*DELEPS(I)
C     ----------------------------------------
C     6. CALCULATION OF ELASTOPLASTIC STRESSES
C     ----------------------------------------
      DO 690  ISINC=1,NSINC
      CALL MEDELP (STRESS,STRAIN,DEPS,DP,IND)
      DO 610  I=1,3
      DO 610  J=1,3
 610  STRESS(I) = STRESS(I) + DP(I,J)*DEPS(J)
      IF (IND.NE.0)  GOTO 630
      DO 620  I=1,3
 620  STRESS(I) = STRESS(I) + DP(I,6)*DEPS(4)
      STRESS(4) = STRESS(4) + DP(6,1)*DEPS(1) + DP(6,2)*DEPS(2) +
     1                        DP(6,3)*DEPS(3) + DP(6,6)*DEPS(4)
C	----------------------------------------
C     FORCE STRESS STATE BACK ON YIELD SURFACE
C	----------------------------------------
 630  SM = (STRESS(1) + STRESS(2) + STRESS(4)) / 3.
      SX = STRESS(1)-SM
      SY = STRESS(2)-SM
      SS = STRESS(3)
      SZ = STRESS(4)-SM
      FTA= 0.5*(SX*SX + SY*SY + SZ*SZ) + SS*SS
      FTB= YLD*YLD/3.
	ABCD=(1./DFF)
      IF (FTA-FTB.LE.ABCD)  GOTO 690
C
      IF (IND.GE.2)  GOTO 650
      COEF = -1. + SQRT(FTB/FTA)
      STRESS(1) = STRESS(1) + COEF*SX
      STRESS(2) = STRESS(2) + COEF*SY
      STRESS(3) = STRESS(3) + COEF*SS
      STRESS(4) = STRESS(4) + COEF*SZ
      GOTO 690
C
C	CHANGES MADE TO NEXT FOUR LINES BY GILSON - JULY2003 (PLASTICITY)
c 650  COEF = SQRT(FTB/FTA)
c      STRESS(1) = STRESS(1)*COEF
c      STRESS(2) = STRESS(2)*COEF
c      STRESS(3) = STRESS(3)*COEF
650   COEF = -1. + SQRT(FTB/FTA)
      STRESS(1) = STRESS(1) + COEF*SX
      STRESS(2) = STRESS(2) + COEF*SY
      STRESS(3) = STRESS(3) + COEF*SS
      STRAIN(4) = STRAIN(4) + (COEF-1.)*SM/BM
C
 690  CONTINUE
C     --------------------------------------
C     7. UPDATING STRESSES,STRAINS AND YIELD
C     --------------------------------------
 700  DO 710  I=1,IST
 710  SIG(I) = STRESS(I)
      DO 720  I=1,ISR
 720  EPS(I) = STRAIN(I)
      IF (IND.EQ.2)  EPS(4) = STRAIN(4)
C     ------------------------------------
C     8. FORM THE MATERIAL LAW (MATRIX DP)
C     ------------------------------------
      IF (IPEL.EQ.1)  CALL DELA3D (DP,IND)
C	SUNIL 09/01/01 REMEOVE 2 FROM FOLLOWING LINE
C      IF (IPEL.EQ.2)  CALL MEDELP (STRESS,STRAIN,DEPS,DP,IND,2)
      IF (IPEL.EQ.2)  CALL MEDELP (STRESS,STRAIN,DEPS,DP,IND)
C
C

      RETURN
      END
C
C=====================================================================
      SUBROUTINE MEDELP (STRESS,STRAIN,DEPS,DP,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     FORMS THE ELASTO-PLASTC STRESS-STRAIN MATRIX
C	--------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     STRESS(4)  = CURRENT TOTAL STRESS
C     STRAIN(4)  = CURRENT TOTAL STRAIN
C     DEPS(4)    = SUBINCREMENT IN STRAIN
C     DP(6,6)    = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     IND        = FLAG (AXISYM.=0,PLANE STRAIN=1,PLANE STRESS=2)
C	-------------------------------------------------------
C     FOR VARIABLES IN COMMON BLOCK /HOOK/ SEE ROUTINE HOKLAW
C     -----------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION STRESS(4),STRAIN(4),DEPS(4),DP(36)
C     --------------------
C     ELASTIC CONTRIBUTION
C     --------------------
      DP(1)  = A2*C2
      DP(8)  = A2*C2
      DP(36) = A2*C2
      DP(2)  = B2*C2
      DP(6)  = B2*C2
      DP(7)  = B2*C2
      DP(12) = DP(2)
      DP(31) = DP(2)
      DP(32) = DP(2)
      DP(3)  = 0.0
      DP(9)  = 0.0
      DP(13) = 0.0
      DP(14) = 0.0
      DP(18) = 0.0
      DP(33) = 0.0
      DP(15) = C1
C
      SM = (STRESS(1)+STRESS(2)+STRESS(4))/3.
      SX = STRESS(1)-SM
      SY = STRESS(2)-SM
      SS = STRESS(3)
      SZ = STRESS(4)-SM
      BETA = 1.5/YLD/YLD
C     -----------------------------
C     CHECK FOR UNLOADING (LAMDA<0)
C     -----------------------------
      IF (IND-1)  50,20,30
 20   WP = SX*DEPS(1) + SY*DEPS(2) + SS*DEPS(3)
      GOTO 100
C
 30   CBETSZ = C2*BETA*SZ
      DP6  = DP(6) - CBETSZ*SX
      DP12 = DP(12)- CBETSZ*SY
      DP18 =       - CBETSZ*SS
      DP36 = DP(36)- CBETSZ*SZ
      DEPS(4) = (-DP6*DEPS(1) - DP12*DEPS(2) - DP18*DEPS(3)) / DP36
C
 50   WP = SX*DEPS(1) + SY*DEPS(2) + SS*DEPS(3) + SZ*DEPS(4)
 100  IF (WP.LT.0.)  GOTO 200
C     -----------------------------
C     SUBTRACT PLASTIC CONTRIBUTION
C     -----------------------------
      C2BETA = C2*BETA
      DP(1)  = DP(1)  - C2BETA*SX*SX
      DP(2)  = DP(2)  - C2BETA*SX*SY
      DP(3)  =        - C2BETA*SX*SS
      DP(6)  = DP(6)  - C2BETA*SX*SZ
      DP(7)  = DP(2)
      DP(8)  = DP(8)  - C2BETA*SY*SY
      DP(9)  =        - C2BETA*SY*SS
      DP(12) = DP(12) - C2BETA*SY*SZ
      DP(13) = DP(3)
      DP(14) = DP(9)
      DP(15) = DP(15) - C2BETA*SS*SS
      DP(18) =        - C2BETA*SS*SZ
      DP(31) = DP(6)
      DP(32) = DP(12)
      DP(33) = DP(18)
      DP(36) = DP(36) - C2BETA*SZ*SZ
 200  IF (IND.LE.1)  RETURN
C     ---------------------------------
C     MODIFY DP MATRIX FOR PLANE STRESS
C     ---------------------------------
      DO 250  I=1,3
      A = DP(30+I)/DP(36)
      DO 250  J=I,3
      KR = (J-1)*6 + I
      KC = (I-1)*6 + J
      DP(KR) = DP(KR) - DP(30+J)*A
 250  DP(KC) = DP(KR)
C     -------------------------------------
C     SET STRAIN INCREMENT AND TOTAL STRAIN
C     -------------------------------------
      IF (WP.LT.0.)  DEPS(4) = D2 * (DEPS(1)+DEPS(2))
      STRAIN(4) = STRAIN(4)+DEPS(4)
C
      RETURN
      END
C
C=====================================================================
C=====================================================================
C	END MEMBRANE ELEMENT (PLANE STRESS, STRAIN & AXISYMMETRIC)
C=====================================================================




C=====================================================================
C	START COMMON ROUTINES
C=====================================================================
C
C=====================================================================
      SUBROUTINE HOKLAW (PROPM,PROPG,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     DETERMINES ALL THE MATERIAL CONSTANTS IN COMMON BLOCK /HOOK/
C	------------------------------------------------------------
C     PROPM(NMP)    = MATERIAL PROPERTIES
C     PROPG(NGP)    = GEOMETRIC PROPERTIES
C     IND           = FLAG FOR PLANE STRAIN(1),OR PLANE STRESS(2)
C     ------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C	INTRODUCE BY DE SILVA
	COMMON /HARD/  HP,DEN
C
      DIMENSION PROPM(1),PROPG(1)
C
      YM  = PROPM(1)
      PR  = PROPM(2)
      YLD = PROPM(3)
      HP  = PROPM(4)
      DEN = PROPM(5)
C
      D2  = PR/(PR-1.)
      C2  = YM/(1.+PR)
      B2  = PR/(1.-2.*PR)
      A2  = (1.-PR)/(1.-2.*PR)
      BM  = YM/(3.-6.*PR)
      C1  = .5*C2
      D1  = C1/1.2
C
      IF (IND.EQ.2)  GOTO 200
C     ----------------------
C     PLANE STRAIN CONDITION
C     ----------------------
      B1 = B2*C2
      A1 = C2+B1
      GOTO 500
C     ----------------------
C     PLANE STRESS CONDITION
C     ----------------------
 200  A1 = YM/(1.-PR*PR)
      B1 = A1*PR
C     -----------------------
C     NODAL ELEMENT THICKNESS
C     -----------------------
 500  TH  = PROPG(2)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHAP2D (R,S,H,P,NODEX,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THERE DERIVATIVES
C     AT THE NODAL POINTS OF A 4 TO 8 NODE,ISOPARAMETRIC QUADRILATERAL
C	----------------------------------------------------------------
C                         NODE NUMBERING CONVENTION
C                         -------------------------
C
C                   2                 5                 1
C
C                     0 . . . . . . . 0 . . . . . . . 0
C                     .                               .
C                     .                               .
C                     .               S               .
C                     .               .               .
C                     .               .               .
C                   6 0               . . . R         0 8
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     0 . . . . . . . 0 . . . . . . . 0
C
C                   3                 7                 4
C
C     R,S        = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     H(8)       = INTERPOLATION (SHAPE) FUNCTIONS
C     P(2,8)     = FUNCTION DERIVATIVES WITH RESPECT TO R,S RESP.
C     NODEX(NEX) = POSITION OF MIDSIDE (EXCESSIVE) NODES
C     NNO        = NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ----------------------------------------------------------------
      DIMENSION  H(8),P(2,8),NODEX(1),IPERM(4)
      DATA (IPERM(I), I=1,4)  /2,3,4,1/
C
      NMI = NNO-4
      RP  = 1.0+R
      SP  = 1.0+S
      RM  = 1.0-R
      SM  = 1.0-S
      R2  = 1.0-R*R
      S2  = 1.0-S*S
C     ---------------------------------------------
C     INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     FOR A FOUR NODE ELEMENT
C     ---------------------------------------------
      H(1)   = 0.25*RP*SP
      H(2)   = 0.25*RM*SP
      H(3)   = 0.25*RM*SM
      H(4)   = 0.25*RP*SM
      P(1,1) = 0.25*SP
      P(1,2) = -P(1,1)
      P(1,3) = -0.25*SM
      P(1,4) = -P(1,3)
      P(2,1) = 0.25*RP
      P(2,2) = 0.25*RM
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)
      IF (NNO.EQ.4)  RETURN
C     -------------------------------------
C     ADD DEGREES OF FREEDOM IN EXCESS OF 4
C     -------------------------------------
      DO 100  IMI=1,NMI
      NN = NODEX(IMI)-4
      GOTO (50,60,70,80), NN
C      CALL GOTOER!NON-EXIST SUBROUTINE
 50   H(5)   = 0.5*R2*SP
      P(1,5) = -R*SP
      P(2,5) = 0.5*R2
      GOTO 100
 60   H(6)   = 0.5*RM*S2
      P(1,6) = -0.5*S2
      P(2,6) = -RM*S
      GOTO 100
 70   H(7)   = 0.5*R2*SM
      P(1,7) = -R*SM
      P(2,7) = -0.5*R2
      GOTO 100
 80   H(8)   = 0.5*RP*S2
      P(1,8) = 0.5*S2
      P(2,8) = -RP*S
 100  CONTINUE
C     ----------------------------------------------
C     CORRECT FUNCTIONS AND DERIVATIVES IF 5 OR MORE
C     NODES ARE USED TO DESCRIBE THE ELENENT
C     ----------------------------------------------
      DO 200  IMI=1,NMI
      IN = NODEX(IMI)
      I1 = IN-4
      I2 = IPERM(I1)
      H(I1)      = H(I1)-0.5*H(IN)
      H(I2)      = H(I2)-0.5*H(IN)
      H(IMI+4)   = H(IN)
      DO 200  J=1,2
      P(J,I1)    = P(J,I1)-0.5*P(J,IN)
      P(J,I2)    = P(J,I2)-0.5*P(J,IN)
 200  P(J,IMI+4) = P(J,IN)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SCAPRD (AV,BV,SC,NV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------
C     FORMS SCALAR PRODUCT SC=AV.BV
C     -----------------------------
      DIMENSION AV(3),BV(3)
C
      SC=0.
      DO 20 I=1,NV
   20 SC=SC+AV(I)*BV(I)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SCALEN (AV,BV,SC,NV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------
C     FORMS SCALAR PRODUCT SC=AV.AV AND STORES
C     NORMALISED COMPONENTS OF AV IN BV
C     ----------------------------------------
C      DIMENSION AV(3),BV(3)
      DIMENSION AV(NV),BV(NV) !CHANGED BY GILSON - OCT2004
C
      SC=0.0D0
      DO 30 I=1,NV
   30 SC=SC+AV(I)*AV(I)
      IF (SC.EQ.0.0D0) RETURN
      SC=DSQRT(SC)
      DO 40 I=1,NV
   40 BV(I)=AV(I)/SC
      RETURN
      END
C
C=====================================================================
      SUBROUTINE VECPRD (AV,BV,CV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------
C     FORMS VECTOR PRODUCT CV=AV*BV
C     -----------------------------
      DIMENSION AV(3),BV(3),CV(3)
C
      CV(1)=AV(2)*BV(3)-AV(3)*BV(2)
      CV(2)=AV(3)*BV(1)-AV(1)*BV(3)
      CV(3)=AV(1)*BV(2)-AV(2)*BV(1)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE ROTVEC (VA,RV,RR,VB)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     ROTATES VA THROUGH AN ANGLE RR (RADS) ABOUT A FIXED AXIS RV
C     INTO NEW ORIENTATION VB (VA,VB AND RV ARE UNIT VECTORS)
C     -----------------------------------------------------------
      DIMENSION VA(3),RV(3),VB(3),VV(3),VVV(3)
C
      SN=SIN(RR)
      CS=1.-COS(RR)
      CALL VECPRD (RV,VA,VV)
      CALL VECPRD (RV,VV,VVV)
      DO 10 I=1,3
   10 VB(I)=VA(I)+SN*VV(I)+CS*VVV(I)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SCAVEC (AV,BV,SC,NV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------
C     VECTOR SCALING - BV(I)=SC*AV(I)
C     -------------------------------
      DIMENSION AV(3),BV(3)
C
      DO 10 I=1,NV
   10 BV(I)=SC*AV(I)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE ADDVEC(AV,BV,CV)
	IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C     ---------------------------------
C     VECTOR ADDITION CV(I)=AV(I)+BV(I)
C     ---------------------------------
      DIMENSION AV(3),BV(3),CV(3)
C
      DO 10 I=1,3
   10 CV(I)=AV(I)+BV(I)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SUBVEC(AV,BV,CV)
	IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C     -----------------------------------
C     VECTOR SUBRACTION CV(I)=AV(I)-BV(I)
C     -----------------------------------
      DIMENSION AV(3),BV(3),CV(3)
C
      DO 20 I=1,3
   20 CV(I)=AV(I)-BV(I)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,ID,f)
	IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
C
C     --------------------------------------------------------------
C     EVALUATES LOCAL DIRECTION COSINE VECTORS, JACOBIAN
C     DETERMINANT, JACOBIAN INVERSE AND BASE VECTOR PARAMETERS
C
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     COORD(3,NNO)      = CURRENT NODAL COORDINATES X,Y,Z
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     XJI(4)            = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     DET               = JACOBIAN DETERMINANT
C     RR,SS             = SQUARED BASE VECTOR LENGTHS
C     SNA               = SIN OF ANGLE SUBTENDED BY BASE VECTORS
C     ID                = FLAG
C	f				  = jacobian coefficient stored in coulmn	
C     COVR(3),COVS(3)   = COVARIENT BASE VECTORS ALONG RN,SN
C     --------------------------------------------------------------
C
      DIMENSION COORD(3,*),HD(2,8),VR(3),VS(3),VT(3),XJI(4)
      DIMENSION COVR(3),COVS(3),RV(3),SV(3),f(4)
C
      CALL CLEARA (COVR,3)
      CALL CLEARA (COVS,3)
      DO 20 I=1,NNO
      DO 20 J=1,3
      COVR(J)=COVR(J)+HD(1,I)*COORD(J,I)
   20 COVS(J)=COVS(J)+HD(2,I)*COORD(J,I)
      CALL VECPRD (COVR,COVS,VT)
      CALL SCALEN (VT,VT,DET,3)
      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)
C	NEXT LINE - CHANGED BY GILSON - JULY2002
C      CALL ADDVEC (SV,VS,VS,3)
      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)
      IF (ID.EQ.2) RETURN
      RR=RL*RL
      SS=SL*SL
      SNA=DET/(RL*SL)
      CALL SCAPRD (COVR,VR,F1,3)
      CALL SCAPRD (COVS,VR,F2,3)
      CALL SCAPRD (COVR,VS,F3,3)
      CALL SCAPRD (COVS,VS,F4,3)
      f(1)=f1
	f(2)=f2
	f(3)=f3
	f(4)=f4
	XJI(1)= F4/DET
      XJI(2)=-F2/DET
      XJI(3)=-F3/DET
      XJI(4)= F1/DET
      SNA=F4/SL
      RETURN
C  f() added by Hari for jacobian------
      END
c
      SUBROUTINE SHBMATrr (NNO,H,HD,VR,VS,VT,XJI,HR,HS,brr,brt,f,k)
	IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)

c      -----------------------------------------------------------------
c      EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
c
c      NNO               = NUMBER OF NODES FOR ELEMENT
C     H(8)              = SHAPE FUNCTIONS
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     XJI(4)            = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)            = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     -----------------------------------------------------------------
C
      DIMENSION H(8),HD(2,8),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION brr(6,48), brt(6,48),f(4)
      l=1
c      M=1
c      N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      F1=A2*VR(J)
      F2=A1*VS(J)
      BM=A1*VR(J)
      BM1=A2*VS(J)
      BM2=F1+F2
      BM3=A1*VT(J)
      BM4=A2*VT(J)
c      B(N)=F2
c      B(N+1)=-F1
c      B(N+2)=B(M+1)-B(M)
      BN3=A3*VS(J)
      BN4=-A3*VR(J)
	brr(k,l)=f(1)*f(1)*BM+f(3)*f(3)*BM1+2.0*f(1)*f(3)*BM2
	brt(k,l)=f(1)*BM3 + f(3)*BM4
	brt(k,l+3)=f(1)*BN3 + f(3)*BN4
c      M=M+5
c      N=N+5
	l=l+1
   10 CONTINUE
c      M=M+15
c      N=N+15
	l=l+3
   20 CONTINUE
      RETURN
C**  -------assumed strain for err, ert----------- **
      END
c
      SUBROUTINE SHBMATss (NNO,H,HD,VR,VS,VT,XJI,HR,HS,bss,bst,f,k)
	IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)

c      -----------------------------------------------------------------
c      EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
c
c      NNO               = NUMBER OF NODES FOR ELEMENT
C     H(8)              = SHAPE FUNCTIONS
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     XJI(4)            = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)            = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     -----------------------------------------------------------------
C
      DIMENSION H(8),HD(2,8),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION bss(6,48),bst(6,48),f(4)
      l=1
c      M=1
c      N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      F1=A2*VR(J)
      F2=A1*VS(J)
      BM=A1*VR(J)
      BM1=A2*VS(J)
      BM2=F1+F2
      BM3=A1*VT(J)
      BM4=A2*VT(J)
c      B(N)=F2
c      B(N+1)=-F1
c      B(N+2)=B(M+1)-B(M)
      BN3=A3*VS(J)
      BN4=-A3*VR(J)
	bss(k,l)=f(2)*f(2)*BM+f(4)*f(4)*BM1+2.0*f(2)*f(4)*BM2
	bst(k,l)=f(2)*BM3 + f(4)*BM4
	bst(k,l+3)=f(2)*BN3+f(4)*BN4
c     M=M+5
c     N=N+5
	l=l+1
   10 CONTINUE
c     M=M+15
c     N=N+15
	l=l+3
   20 CONTINUE
      RETURN
C**  ---assumed strain for ess and est---**
      END
c
      SUBROUTINE SHBMATrs (NNO,H,HD,VR,VS,VT,XJI,HR,HS,brs,f,k)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

c      -----------------------------------------------------------------
c      EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
c
c      NNO               = NUMBER OF NODES FOR ELEMENT
C     H(8)              = SHAPE FUNCTIONS
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     XJI(4)            = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)            = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     -----------------------------------------------------------------
C
      DIMENSION H(8),HD(2,8),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION brs(4,48),f(4)
      l=1
c     M=1
c     N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      F1=A2*VR(J)
      F2=A1*VS(J)
      BM=A1*VR(J)
      BM1=A2*VS(J)
      BM2=F1+F2
      BM3=A1*VT(J)
      BM4=A2*VT(J)
c      B(N)=F2
c      B(N+1)=-F1
c      B(N+2)=B(M+1)-B(M)
c      B(N+3)=A3*VS(J)
c      B(N+4)=-A3*VR(J)
      brs(k,l)=f(1)*f(2)*BM+f(3)*f(4)*BM1+(f(1)*f(4)+f(2)*f(3))*BM2
c     M=M+5
c      N=N+5
      l=l+1
   10 CONTINUE
c     M=M+15
c      N=N+15
	l=l+3
   20 CONTINUE
      RETURN
C**  -------assumed strain for ers ---- **
      END	  	   
c
      SUBROUTINE SHBMAT8(NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,brr,bss,brs,brt
     +                	,bst,rn,sn)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
c      -----------------------------------------------------------------
c      EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
c
c      NNO               = NUMBER OF NODES FOR ELEMENT
C     H(8)              = SHAPE FUNCTIONS
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     XJI(4)            = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)            = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     -----------------------------------------------------------------
C
      DIMENSION H(8),HD(2,8),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION B(240),brr(6,48),bss(6,48),brs(4,48),brt(6,48),bst(6,48)
C      M=1
      N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
C      A1=HR(I)
C      A2=HS(I)
C      A3=H(I)
      DO 10 J=1,3
C      F1=A2*VR(J)
C      F2=A1*VS(J)
      B(N)=VS(J)*HR(I)
      B(N+1)=-VR(J)*HS(I)
      B(N+2)=VS(J)*HS(I)-VR(J)*HR(I)
C      M=M+5
      N=N+5
   10 CONTINUE
C      M=M+15
      N=N+15
   20 CONTINUE
      a=sqrt(3.)
      r1 =0.25*(1+a*rn)* sn*(sn+1)
	r2 =0.25*(1-a*rn)* sn*(sn+1)
	r3 = 0.5*(1+a*rn)* (1-sn*sn)
	r4 = 0.5*(1-a*rn)* (1-sn*sn)
	r5 =0.25*(1+a*rn)* sn*(sn-1)
	r6 = 0.25*(1-a*rn)*sn*(sn-1)
c
	s1 =0.25*(1+a*sn)* rn*(rn+1)
	s2 =0.25*(1-a*sn)* rn*(rn+1)
	s3 = 0.5*(1+a*sn)* (1-rn*rn)
	s4 = 0.5*(1-a*sn)* (1-rn*rn)
	s5 =0.25*(1+a*sn)* rn*(rn-1)
	s6 = 0.25*(1-a*sn)*rn*(rn-1)
c	
	rs1 = 0.25*(1+a*rn)*(1+a*sn)
	rs2 = 0.25*(1-a*rn)*(1+a*sn)
	rs3 = 0.25*(1+a*rn)*(1-a*sn)
	rs4 = 0.25*(1-a*rn)*(1-a*sn)
c
	M=1
	L=1
	DO 31 N=1,8
	DO 30 J=1,3
	brrn=r1*brr(1,l) +r2*brr(2,l) + r3*brr(3,l)
     +    +r4*brr(4,l) +r5*brr(5,l) + r6*brr(6,l)           	 
c
     	bssn=s1*bss(1,l) +s2*bss(2,l) + s3*bss(3,l)
     +    +s4*bss(4,l) +s5*bss(5,l) + s6*bss(6,l)           	 
c     
	brsn=rs1*brs(1,l) +rs2*brs(2,l) +rs3*brs(3,l) + rs4*brs(4,l)
c
	B(M) =brrn*xji(1)*xji(1)+bssn*xji(3)*xji(3)+brsn*2.0*xji(1)*xji(3)
	B(M+1)=brrn*xji(2)*xji(2)+bssn*xji(4)*xji(4)+brsn*2.*xji(2)*xji(4)
	B(M+2)=brrn*xji(1)*xji(2)+bssn*xji(3)*xji(4)+brsn*(xji(1)*xji(4)
     +                                                 +xji(2)*xji(3))
      M=M+5 
      l=l+1
  30	continue
      M=M+15
      l=l+3
  31	continue
      M=1
	do 32 L=1,48
	brtn=r1*brt(1,l) +r2*brt(2,l) + r3*brt(3,l)
     +    +r4*brt(4,l) +r5*brt(5,l) + r6*brt(6,l)           	 
c
     	bstn=s1*bst(1,l) +s2*bst(2,l) + s3*bst(3,l)
     +    +s4*bst(4,l) +s5*bst(5,l) + s6*bst(6,l) 
      B(M+3)=xji(1)*brtn+xji(3)*bstn
	B(M+4)=xji(2)*brtn+xji(4)*bstn
      M=M+5
  32  continue
      RETURN
C** ---------------------  **
      END
C=====================================================================
      SUBROUTINE SHGEOA (S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	BY GILSON - MAY2002
C     ---------------------------------------------------------------
C     EVALUATES GEOMETRIC CONTRIBUTION TO TANGENTIAL STIFFNESS MATRIX
C     THIS VERSION IS BASED ON FULL GREEN'S STRAIN EXPANSION
C	------------------------------------------------------
C     S(1176)           = GEOMETRIC STIFFNESS STORED UPPER TRIANGULAR
C                         ROW-WISE
C     SIGR(8)           = WEIGHTED STRESS-RESULTANT VECTOR
C     H(8)              = SHAPE FUNCTIONS
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     TH                = GAUSS POINT THICKNESS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     ---------------------------------------------------------------
      DIMENSION S(*),SIGR(8),H(8),HR(8),HS(8),VR(3),VS(3),VT(3),F(6)
	DIMENSION RTR(9),TRT(9),STS(9),TST(9)

C	--------------------
C	CONTRIBUTION OF KG11
C	--------------------
	IP=1
	DO 10 I=1,NNO
	 DO 20 J=1,3
	  DO 30 K=I,NNO
	   FF=HR(I)*HR(K)
	   GG=HS(I)*HS(K)
	   FG=HR(I)*HS(K)
	   GF=HS(I)*HR(K)
	   S(IP)=S(IP)+SIGR(1)*FF+SIGR(2)*GG+SIGR(3)*(FG+GF)
	   IP=IP+6
30	  CONTINUE
	  IP=IP+1-J
20	 CONTINUE
	 IP=IP+(NNO-I)*(18)+6
10	CONTINUE

C	--------------------
C	CONTRIBUTION OF KG12
C	--------------------
	IP=4
	DO 40 I=1,NNO
	 DO 50 J=1,3
	  DO 60 K=I,NNO
	   FF=HR(I)*HR(K)
	   GG=HS(I)*HS(K)
	   FG=HR(I)*HS(K)
	   GF=HS(I)*HR(K)
	   FH=HR(I)*H(K)
	   GH=HS(I)*H(K)
	   SELECT CASE (J)
	    CASE(1)
	     FVF1 =  0.0
	     FVF2 =  FF*VT(3)
	     FVF3 = -FF*VT(2)
	     GVG1 =  0.0
	     GVG2 =  GG*VT(3)
	     GVG3 = -GG*VT(2)
	     FVG1 =  0.0
	     FVG2 =  FG*VT(3)
	     FVG3 = -FG*VT(2)
	     GVF1 =  0.0
	     GVF2 =  GF*VT(3)
	     GVF3 = -GF*VT(2)
	     FVH1 =  0.0
	     FVH2 =  FH*VT(3)
	     FVH3 = -FH*VT(2)
	     GVH1 =  0.0
	     GVH2 =  GH*VT(3)
	     GVH3 = -GH*VT(2)
	    CASE(2)
	     FVF1 = -FF*VT(3)
	     FVF2 =  0.0
	     FVF3 =  FF*VT(1)
	     GVG1 = -GG*VT(3)
	     GVG2 =  0.0
	     GVG3 =  GG*VT(1)
	     FVG1 = -FG*VT(3)
	     FVG2 =  0.0
	     FVG3 =  FG*VT(1)
	     GVF1 = -GF*VT(3)
	     GVF2 =  0.0
	     GVF3 =  GF*VT(1)
	     FVH1 = -FH*VT(3)
	     FVH2 =  0.0
	     FVH3 =  FH*VT(1)
	     GVH1 = -GH*VT(3)
	     GVH2 =  0.0
	     GVH3 =  GH*VT(1)
	    CASE(3)
	     FVF1 =  FF*VT(2)
	     FVF2 = -FF*VT(1)
	     FVF3 =  0.0
	     GVG1 =  GG*VT(2)
	     GVG2 = -GG*VT(1)
	     GVG3 =  0.0
	     FVG1 =  FG*VT(2)
	     FVG2 = -FG*VT(1)
	     FVG3 =  0.0
	     GVF1 =  GF*VT(2)
	     GVF2 = -GF*VT(1)
	     GVF3 =  0.0
	     FVH1 =  FH*VT(2)
	     FVH2 = -FH*VT(1)
	     FVH3 =  0.0
	     GVH1 =  GH*VT(2)
	     GVH2 = -GH*VT(1)
	     GVH3 =  0.0
	   END SELECT
	   S(IP)=S(IP)+SIGR(4)*FVF1+SIGR(5)*GVG1+SIGR(6)*(FVG1+GVF1)+
	1               SIGR(7)*FVH1+SIGR(8)*GVH1
	   S(IP+1)=S(IP+1)+SIGR(4)*FVF2+SIGR(5)*GVG2+SIGR(6)*(FVG2+GVF2)+
	1                   SIGR(7)*FVH2+SIGR(8)*GVH2
	   S(IP+2)=S(IP+2)+SIGR(4)*FVF3+SIGR(5)*GVG3+SIGR(6)*(FVG3+GVF3)+
	1                   SIGR(7)*FVH3+SIGR(8)*GVH3
	   IP=IP+6
60	  CONTINUE
	  IP=IP-J
50	 CONTINUE
	 IP=IP+(NNO-I)*(18)+9
40	CONTINUE

C	--------------------
C	CONTRIBUTION OF KG22
C	--------------------
C	INITIAL POINTER
	IP=(NNO-1)*18+16
C	PRODUCT OF COORDINATE VECTORS - ROW-WISE(LEFT-RIGHT,TOP-BOTTOM)
	IN = 0
	DO 100 II=1,3
	  DO 100 JJ=1,3
	    IN=IN+1
	    RTR(IN) = VR(II)*VT(JJ)+VR(JJ)*VT(II)
	    TRT(IN) = VT(II)*VR(JJ)+VT(JJ)*VR(II)
	    STS(IN) = VS(II)*VT(JJ)+VS(JJ)*VT(II)
	    TST(IN) = VT(II)*VS(JJ)+VT(JJ)*VS(II)
100	CONTINUE
C	ASSEMBLE KG22 IN S(*)
	DO 70 I=1,NNO
C	  6X6 SUBMATRIXES ALONG THE DIAGONAL
C	  POINTERS
	  IP1 = IP
	  IP2 = IP+3+6*(NNO-I)
	  IP3 = IP2+2+6*(NNO-I)
C	  PRODUCT OF SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES
	  FH=HR(I)*H(I)
	  GH=HS(I)*H(I)
	  HH=H(I)*H(I)
C	  GEOMETRIC STIFFNESS
	  DO 110 II=1,3
	    S(IP1+II-1) = S(IP1+II-1)+0.5*SIGR(4)*FH*(RTR(II)+TRT(II))+
	1                              0.5*SIGR(5)*GH*(STS(II)+TST(II))+
	2                              0.5*SIGR(6)*FH*(STS(II)+TST(II))+
	3                              0.5*SIGR(6)*GH*(RTR(II)+TRT(II))+
	4                          0.5*HH*(SIGR(7)*RTR(II)+SIGR(8)*STS(II))
110	  CONTINUE
	  DO 120 II=1,2
	   S(IP2+II-1) = S(IP2+II-1)+0.5*SIGR(4)*FH*(RTR(II+4)+TRT(II+4))+
	1                             0.5*SIGR(5)*GH*(STS(II+4)+TST(II+4))+
	2                             0.5*SIGR(6)*FH*(STS(II+4)+TST(II+4))+
	3                             0.5*SIGR(6)*GH*(RTR(II+4)+TRT(II+4))+
	4                      0.5*HH*(SIGR(7)*RTR(II+4)+SIGR(8)*STS(II+4))
120	  CONTINUE
	  S(IP3) = S(IP3)+0.5*SIGR(4)*FH*(RTR(9)+TRT(9))+
	1                  0.5*SIGR(5)*GH*(STS(9)+TST(9))+
	2                  0.5*SIGR(6)*FH*(STS(9)+TST(9))+
	3                  0.5*SIGR(6)*GH*(RTR(9)+TRT(9))+
	4                  0.5*HH*(SIGR(7)*RTR(9)+SIGR(8)*STS(9))
C	  OFF DIAGONAL 6X6 SUBMATRICES
	  K=I+1
	  IF (K.NE.NNO) THEN
	  DO 90 J=K,NNO
C	  POINTERS
	    IP1 = IP1+6
	    IP2 = IP2+6
	    IP3 = IP3+6
C	    PRODUCT OF SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES
	    HF=H(I)*HR(J)
	    FH=HR(I)*H(J)
	    HG=H(I)*HS(J)
	    GH=HS(I)*H(J)
	    HH=H(I)*H(J)
C	    GEOMETRIC STIFFNESS
	    DO 130 II=1,3
	      S(IP1+II-1) = S(IP1+II-1)+
	1                    0.5*SIGR(4)*(HF*RTR(II)+FH*TRT(II))+
	2	                0.5*SIGR(5)*(HG*STS(II)+GH*TST(II))+
	3                    0.5*SIGR(6)*(HF*STS(II)+FH*TST(II))+
	4                    0.5*SIGR(6)*(HG*RTR(II)+GH*TRT(II))+
	5                    0.5*HH*(SIGR(7)*RTR(II)+SIGR(8)*STS(II))
130	    CONTINUE
	    DO 140 II=1,3
	      S(IP2+II-2) = S(IP2+II-2)+
	1                    0.5*SIGR(4)*(HF*RTR(II+3)+FH*TRT(II+3))+
	2	                0.5*SIGR(5)*(HG*STS(II+3)+GH*TST(II+3))+
	3                    0.5*SIGR(6)*(HF*STS(II+3)+FH*TST(II+3))+
	4                    0.5*SIGR(6)*(HG*RTR(II+3)+GH*TRT(II+3))+
	5                    0.5*HH*(SIGR(7)*RTR(II+3)+SIGR(8)*STS(II+3))
140	    CONTINUE
	    DO 150 II=1,3
	      S(IP3+II-3) = S(IP3+II-3)+
	1                    0.5*SIGR(4)*(HF*RTR(II+6)+FH*TRT(II+6))+
	2	                0.5*SIGR(5)*(HG*STS(II+6)+GH*TST(II+6))+
	3                    0.5*SIGR(6)*(HF*STS(II+6)+FH*TST(II+6))+
	4                    0.5*SIGR(6)*(HG*RTR(II+6)+GH*TRT(II+6))+
	5                    0.5*HH*(SIGR(7)*RTR(II+6)+SIGR(8)*STS(II+6))
150	    CONTINUE

90	  CONTINUE
	  ENDIF

	  IP=IP+18*(2*NNO-2*I-1)+21
70	CONTINUE

      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHGEO (S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     EVALUATES GEOMETRIC CONTRIBUTION TO TANGENTIAL STIFFNESS MATRIX
C     THIS VERSION IS BASED ON FULL GREEN'S STRAIN EXPANSION
C	------------------------------------------------------
C     S(1176)           = GEOMETRIC STIFFNESS STORED UPPER TRIANGULAR
C                         ROW-WISE
C     SIGR(8)           = WEIGHTED STRESS-RESULTANT VECTOR
C     H(8)              = SHAPE FUNCTIONS
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     TH                = GAUSS POINT THICKNESS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     ---------------------------------------------------------------
      DIMENSION S(*),SIGR(8),H(8),HR(8),HS(8),VR(3),VS(3),VT(3),F(6)
	DIMENSION RTR(9),TRT(9),STS(9),TST(9)
C
      F1=VR(1)*VR(1)+VS(1)*VS(1)
      F2=VR(1)*VR(2)+VS(1)*VS(2)
      F3=VR(1)*VR(3)+VS(1)*VS(3)
      F4=VR(2)*VR(2)+VS(2)*VS(2)
      F5=VR(2)*VR(3)+VS(2)*VS(3)
      F6=VR(3)*VR(3)+VS(3)*VS(3)
      FB=TH*TH/12.
      NS=1
      LROW=6*NNO
      DO 30 I=1,NNO
      A1=SIGR(1)*HR(I)+SIGR(3)*HS(I)
      A2=SIGR(3)*HR(I)+SIGR(2)*HS(I)
      A3=SIGR(4)*HR(I)+SIGR(6)*HS(I)
      A4=SIGR(6)*HR(I)+SIGR(5)*HS(I)
      A5=SIGR(7)*HR(I)+SIGR(8)*HS(I)
      HH=.5*H(I)
      HHR=.5*HR(I)
      HHS=.5*HS(I)
      B1=SIGR(7)*HH+SIGR(4)*HHR+SIGR(6)*HHS
      B2=SIGR(8)*HH+SIGR(6)*HHR+SIGR(5)*HHS
      B3=SIGR(4)*HH
      B4=SIGR(5)*HH
      B5=SIGR(6)*HH
      DO 20 J=I,NNO
      N=NS
      A=HR(J)*A1+HS(J)*A2
      E=HR(J)*A3+HS(J)*A4
      B=H(I)*(SIGR(7)*HR(J)+SIGR(8)*HS(J))+E
      C=H(J)*A5+E
      D=A*FB
      S1=H(J)*B1+HR(J)*B3+HS(J)*B5
      S2=H(J)*B2+HR(J)*B5+HS(J)*B4
      M=1
      DO 10 II=1,3
      DO 10 JJ=II,3
      F(M)=S1*(VR(II)*VT(JJ)+VT(II)*VR(JJ))
     1    +S2*(VS(II)*VT(JJ)+VT(II)*VS(JJ))
   10 M=M+1
C     ----------------------------------------------
C     UPPER TRIANGULAR PORTION OF EACH 6*6 PARTITION
C     ----------------------------------------------
      S(N)=S(N)+A
      S(N+4)=S(N+4)+C*VT(3)
      S(N+5)=S(N+5)-C*VT(2)
      N=N+LROW
      S(N)=S(N)+A
      S(N+2)=S(N+2)-C*VT(3)
      S(N+4)=S(N+4)+C*VT(1)
      N=N+LROW-1
      S(N)=S(N)+A
      S(N+1)=S(N+1)+C*VT(2)
      S(N+2)=S(N+2)-C*VT(1)
      N=N+LROW-2
C      S(N)=S(N)+D*F1+F(1)
C      S(N+1)=S(N+1)+D*F2+F(2)
C      S(N+2)=S(N+2)+D*F3+F(3)
      N=N+LROW-3
C      S(N)=S(N)+D*F4+F(4)
C      S(N+1)=S(N+1)+D*F5+F(5)
      N=N+LROW-4
C      S(N)=S(N)+D*F6+F(6)
      IF (I.EQ.J) GOTO 20
C     --------------------------------------------------------------
C     FILL IN REMAINING COEFFICIENTS FOR OFF-DIAGONAL 6*6 PARTITIONS
C     --------------------------------------------------------------
      N=NS+3*(LROW-2)
      S(N+1)=S(N+1)-B*VT(3)
      S(N+2)=S(N+2)+B*VT(2)
      N=N+LROW-4
      S(N)=S(N)+B*VT(3)
      S(N+2)=S(N+2)-B*VT(1)
C      S(N+3)=S(N+3)+D*F2+F(2)
      N=N+LROW-5
      S(N)=S(N)-B*VT(2)
      S(N+1)=S(N+1)+B*VT(1)
C      S(N+3)=S(N+3)+D*F3+F(3)
C      S(N+4)=S(N+4)+D*F5+F(5)
   20 NS=NS+6
      NS=N+6
   30 LROW=LROW-6

C	--------------------
C	CONTRIBUTION OF KG22
C	--------------------
C	INITIAL POINTER
	IP=(NNO-1)*18+16
C	PRODUCT OF COORDINATE VECTORS - ROW-WISE(LEFT-RIGHT,TOP-BOTTOM)
	IN = 0
	DO 100 II=1,3
	  DO 100 JJ=1,3
	    IN=IN+1
	    RTR(IN) = VR(II)*VT(JJ)+VR(JJ)*VT(II)
	    TRT(IN) = VT(II)*VR(JJ)+VT(JJ)*VR(II)
	    STS(IN) = VS(II)*VT(JJ)+VS(JJ)*VT(II)
	    TST(IN) = VT(II)*VS(JJ)+VT(JJ)*VS(II)
100	CONTINUE
C	ASSEMBLE KG22 IN S(*)
	DO 70 I=1,NNO
C	  6X6 SUBMATRIXES ALONG THE DIAGONAL
C	  POINTERS
	  IP1 = IP
	  IP2 = IP+3+6*(NNO-I)
	  IP3 = IP2+2+6*(NNO-I)
C	  PRODUCT OF SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES
	  FH=HR(I)*H(I)
	  GH=HS(I)*H(I)
	  HH=H(I)*H(I)
C	  GEOMETRIC STIFFNESS
	  DO 110 II=1,3
	    S(IP1+II-1) = S(IP1+II-1)+0.5*SIGR(4)*FH*(RTR(II)+TRT(II))+
	1                              0.5*SIGR(5)*GH*(STS(II)+TST(II))+
	2                              0.5*SIGR(6)*FH*(STS(II)+TST(II))+
	3                              0.5*SIGR(6)*GH*(RTR(II)+TRT(II))+
	4                          0.5*HH*(SIGR(7)*RTR(II)+SIGR(8)*STS(II))
110	  CONTINUE
	  DO 120 II=1,2
	   S(IP2+II-1) = S(IP2+II-1)+0.5*SIGR(4)*FH*(RTR(II+4)+TRT(II+4))+
	1                             0.5*SIGR(5)*GH*(STS(II+4)+TST(II+4))+
	2                             0.5*SIGR(6)*FH*(STS(II+4)+TST(II+4))+
	3                             0.5*SIGR(6)*GH*(RTR(II+4)+TRT(II+4))+
	4                      0.5*HH*(SIGR(7)*RTR(II+4)+SIGR(8)*STS(II+4))
120	  CONTINUE
	  S(IP3) = S(IP3)+0.5*SIGR(4)*FH*(RTR(9)+TRT(9))+
	1                  0.5*SIGR(5)*GH*(STS(9)+TST(9))+
	2                  0.5*SIGR(6)*FH*(STS(9)+TST(9))+
	3                  0.5*SIGR(6)*GH*(RTR(9)+TRT(9))+
	4                  0.5*HH*(SIGR(7)*RTR(9)+SIGR(8)*STS(9))
C	  OFF DIAGONAL 6X6 SUBMATRICES
	  K=I+1
	  IF (K.NE.NNO) THEN
	  DO 90 J=K,NNO
C	  POINTERS
	    IP1 = IP1+6
	    IP2 = IP2+6
	    IP3 = IP3+6
C	    PRODUCT OF SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES
	    HF=H(I)*HR(J)
	    FH=HR(I)*H(J)
	    HG=H(I)*HS(J)
	    GH=HS(I)*H(J)
	    HH=H(I)*H(J)
C	    GEOMETRIC STIFFNESS
	    DO 130 II=1,3
	      S(IP1+II-1) = S(IP1+II-1)+
	1                    0.5*SIGR(4)*(HF*RTR(II)+FH*TRT(II))+
	2	                0.5*SIGR(5)*(HG*STS(II)+GH*TST(II))+
	3                    0.5*SIGR(6)*(HF*STS(II)+FH*TST(II))+
	4                    0.5*SIGR(6)*(HG*RTR(II)+GH*TRT(II))+
	5                    0.5*HH*(SIGR(7)*RTR(II)+SIGR(8)*STS(II))
130	    CONTINUE
	    DO 140 II=1,3
	      S(IP2+II-2) = S(IP2+II-2)+
	1                    0.5*SIGR(4)*(HF*RTR(II+3)+FH*TRT(II+3))+
	2	                0.5*SIGR(5)*(HG*STS(II+3)+GH*TST(II+3))+
	3                    0.5*SIGR(6)*(HF*STS(II+3)+FH*TST(II+3))+
	4                    0.5*SIGR(6)*(HG*RTR(II+3)+GH*TRT(II+3))+
	5                    0.5*HH*(SIGR(7)*RTR(II+3)+SIGR(8)*STS(II+3))
140	    CONTINUE
	    DO 150 II=1,3
	      S(IP3+II-3) = S(IP3+II-3)+
	1                    0.5*SIGR(4)*(HF*RTR(II+6)+FH*TRT(II+6))+
	2	                0.5*SIGR(5)*(HG*STS(II+6)+GH*TST(II+6))+
	3                    0.5*SIGR(6)*(HF*STS(II+6)+FH*TST(II+6))+
	4                    0.5*SIGR(6)*(HG*RTR(II+6)+GH*TRT(II+6))+
	5                    0.5*HH*(SIGR(7)*RTR(II+6)+SIGR(8)*STS(II+6))
150	    CONTINUE

90	  CONTINUE
	  ENDIF

	  IP=IP+18*(2*NNO-2*I-1)+21
70	CONTINUE



      RETURN
      END
C
C=====================================================================
      SUBROUTINE MASS2D (CM,H,DE,DVOL,NNO,NEF,IPT)
C	PROGRAMED BY GILSON - DEC2001
      IMPLICIT REAL*8(A-H,O-Z)
C	------------------------------------------------------------------
C	PURPOSE: CALCULATE THE MASS MATRIX OF A QUADRILATERAL ELEMENT
C
C	DESCRIPTION OF VARIABLES
C		CM		= ELEMENT MASS MATRIX
C		DE		= DENSITY
C		DVOL	= VOLUME FACTOR
C		FAC		= MASS FACTOR
C		H,D		= SHAPE FUNCTION
C		IPT		= INTEGRATION POINT NUMBER
C		NNO		= NUMBER OF ELEMENT NODES
C		IMASS	= FLAG FOR TYPE OF MASS MATRIX
C					1,2 - CONSISTENTLY LUMPED MASS MATRIX
C					3   - CONSISTENT MASS MATRIX
C		ISTYP	= TYPE OF MEMBRANE ELEMENT
C					0 - AXISYMMETRIC
C					1 - PLANE STRAIN
C					2 - PLANE STRESS
C		KL,KS	= POINTERS
C		NEF		= NUMBER OF DEGREES OF FREEDOM (2 * NUMBER OF NODES)
C	------------------------------------------------------------------
      COMMON /DYNA/ CDEN,IMASS
C
      DIMENSION CM(1),D(16),H(8)
C
	FAC=DVOL*DE
C
C	---------------
C     CONSISTENT MASS
C	---------------
      IF (IMASS.LT.3) GO TO 320
      DO 200 I = 1,NNO
      D(2*I-1) = H(I)
  200 D(2*I) = H(I)
      KL=1
      DO 300 I=1,NEF,2
      DO 301 J=I,NEF,2
      CM(KL)=CM(KL) + D(I)*D(J)*FAC
  301 KL=KL + 2
  300 KL=KL + NEF - I
C
      KL=1
      DO 401 I=1,NEF,2
      KS=KL + NEF - I + 1
      DO 400 J=I,NEF,2
      CM(KS)=CM(KL)
      KS=KS + 2
  400 KL=KL + 2
  401 KL=KL + NEF - I
C
      RETURN
C
C	-----------
C     LUMPED MASS
C	-----------
  320 DO 325 I=1,NEF,2
      FACM=FAC/NNO
  325 CM(I)=CM(I) + FACM
C
C
  335 DO 340 I=1,NEF,2
  340 CM(I+1)=CM(I)
C
      RETURN
C
C
      END

C
C=====================================================================
