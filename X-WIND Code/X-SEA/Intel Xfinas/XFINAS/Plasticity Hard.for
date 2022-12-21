C	=============================================================
C	START ISOTROPIC LINEAR HARDENNING SUB ROUTINES 
C			( IVANOVES AND VON MISSES )
C
C	=============================================================
C	VON MISSES HARDENNING ROUTINES
C	=============================================================
      SUBROUTINE MLAYERH (WA,STRAIN,STRESS,DP,DVOL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     SET UP LOOP OVER NUMBER OF LAYERS
C     CALL VMISES TO COMPUTE ELASTO-PLASTIC STRESSES AND STRESS RES.
C     INTEGRATES ELESTO-PLASTIC MATERIAL MATRIX DP OVER THE DEPTH
C	-----------------------------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     WA(8,NGT+1) = WORKING ARRAY STORING STRESSES AND STRAINS
C     STRAIN(8)   = CURRENT TOTAL STRAINS
C     STRESS(8)   = RETURNED STRESS RESULTANTS
C     DP(8,8)     = MATRIX CONTAINING ELASTO-PLASTIC RIGIDITIES
C     DVOL        = WT*DET
C	----------------
C     LOCAL VARIABLES
C	----------------
C     EPSZ(4)     = STRAINS AT LEVEL Z
C     SIGZ(4)     = STRESSES AT LEVEL Z
C     CP(6,6)     = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     DEPTH(26)   = Z-COORDINATES OF LAYER BOUNDARIES
C     --------------------------------------------------------------
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /FLTH/  THICK
	
C
      DIMENSION WA(8,1),STRAIN(8),STRESS(8),DP(8,8)
      DIMENSION EPSZ(4),SIGZ(4),CP(6,6),DEPTH(26)
      EQUIVALENCE (IPE,APEL)
C
      THICK = TH
      IND = 3
      ISR = 3
      IST = 3
      IPEL = 1
      CALL CLEARA (STRESS,8)
      IF (IFREF.EQ.0) CALL CLEARA (DP,64)
      IF (NGT.GT.1)  GOTO 400
C     --------------------------
C     PURE MEMBRANE STRESS STATE
C     --------------------------
      DO 110  I=1,3
 110  EPSZ(I) = STRAIN(I)
      CALL VMISES (WA(1,1),WA(4,1),WA(8,1),EPSZ,SIGZ,CP,IND)

      DO 120  I=1,3
      STRESS(I) = SIGZ(I)*THICK
 120  WA(8+I,1) = STRESS(I)
      IF (IFREF.NE.0)  RETURN
C     -------------------------
C     ELASTO-PLASTIC RIGIDITIES
C     -------------------------
      TH3 = THICK*THICK*THICK/12.
      FACM = DVOL*THICK
      FACB = DVOL*TH3
      DO 150  I=1,3
      DO 150  J=I,3
      K = I+3
      L = J+3
      DP(I,J) = CP(I,J)*FACM
      DP(K,L) = CP(I,J)*FACB
      DP(J,I) = DP(I,J)
 150  DP(L,K) = DP(K,L)
      DP(7,7) = D1*FACM
      DP(8,8) = D1*FACM
      RETURN
C     ------------------------------------
C	MULTI-LAYER SOLUTION, INITIALISATION
C     ------------------------------------
c	large strain contribution(update thickness)
c      EPSZ(4) = -(D(6)*STRAIN(1)+D(12)*STRAIN(2)+D(18)*STRAIN(3))
c     1            /D(36)
c      EXT  = 1. - 2.*STRAIN(4)
c      EXT  = 1. - SQRT(EXT)
C
c      THICK = THICK*(1.-EXT)

 400  DZ = THICK/(NGT-1)
      Z  = -THICK/2-DZ
C     -------------------------------------------------------
C     LOOP OVER NUMBER OF LAYERS,CALCULATE STRAINS AT LEVEL Z
C     -------------------------------------------------------
C	NGT-NUMBER OF GOUSS POINT IN Z DIRECTION
      DO 800  IGT=1,NGT
      Z = Z+DZ
      DO 550  I=1,3
 550  EPSZ(I) = STRAIN(I) + Z*STRAIN(I+3)
C     ----------------------
C     COMPUTE LAYER STRESSES
C     ----------------------
      CALL VMISESH (WA(1,IGT),WA(4,IGT),EPSZ,SIGZ,
	1			  CP,IND,WA(8,IGT))

C	FOR HARDENNING IT CAN BE USED WA(7,IGT) TO UPDATE YIELD STRESS
C	CHECK FOR ITS RELIABILITY.....................................2003/15
      APEL = WA(8,IGT)
      IF (IPE.EQ.2)  IPEL = 2
C     -------------------------------------------
C     ADD LAYER CONTRIBUTION TO STRESS RESULTANTS
C     -------------------------------------------
      WTMEM=DZ
      ZZ = Z
      IF (IGT.NE.1)  GOTO 620
      WTMEM = DZ/2.
      ZZ = Z+DZ/3.
 620  IF (IGT.NE.NGT)  GOTO 650
      WTMEM = DZ/2.
      ZZ = Z-DZ/3.
 650  WTBEN = ZZ*WTMEM
      DO 690  I=1,3
      J = I+3
      STRESS(I) = STRESS(I) + WTMEM*SIGZ(I)
 690  STRESS(J) = STRESS(J) + WTBEN*SIGZ(I)
      IF (IFREF.NE.0)  GOTO 800
C     -------------------------
C     INTEGRATE MATERIAL MATRIX
C     -------------------------
      FACM = WTMEM*DVOL
      FACC = Z*FACM
      FACB = Z*FACC
      DO 750  I=1,3
      DO 750  J=I,3
      K = I+3
      L = J+3
      DP(I,J) = DP(I,J) + CP(I,J)*FACM
 750  DP(K,L) = DP(K,L) + CP(I,J)*FACB
      DO 760  I=1,3
      DO 760  J=1,3
      L = J+3
 760  DP(I,L) = DP(I,L) + CP(I,J)*FACC
C
 800  CONTINUE


c	large strain contribution(update thickness)
c      EPSZ(4) = -(D(6)*EPSZ(1)+D(12)*EPSZ(2)+D(18)*EPSZ(3))
c     1            /D(36)
c      EXT  = 1. - 2.*EPSZ(4)
c      EXT  = 1. - SQRT(EXT)
C
c      THICK = THICK*(1.-EXT)

C     ---------------------------------------------------
C     FILL IN LOWER TRIANGLE OF RIGIDITY MATRIX (IFREF=0)
C     ---------------------------------------------------
      IF (IFREF.NE.0)  GOTO 900
      IF (IPEL.EQ.1)   GOTO 850
      DO 830  I=1,6
      DO 830  J=I,6
 830  DP(J,I) = DP(I,J)
      DP(7,7) = D1*THICK*DVOL
      DP(8,8) = D1*THICK*DVOL
      GOTO 900
C
 850  	SBB=1.
	SAA=1.
	CALL SHDELA (DP,DVOL,SBB,SAA)
C     ------------------------------------
C     STORE STRESS RESULTANTS AT END OF WA
C     ------------------------------------
C	THIS IS FOR PRINTING STERSS RESULTANTS( BENDING MOMENT AND NORMAL FORCE)
 900  STRESS(7) = D1*STRAIN(7)*THICK
      STRESS(8) = D1*STRAIN(8)*THICK
      DO 990  I=1,8
 990  WA(I,NGT+1) = STRESS(I)
C
C
      RETURN
      END
C	=====================================================================
      SUBROUTINE VMISESH (SIG,EPS,STRAIN,STRESS,DP,IND,YLD1)
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
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
	COMMON /HARD/  HP,DEN
      DIMENSION SIG(4),EPS(4),STRAIN(4),STRESS(4),DP(6,6)
      DIMENSION DELEPS(4),DELSIG(4),DEPS(4)
      EQUIVALENCE (DELEPS(1),DEPS(1))
C     --------------------------------------------------------------
C     1. CALCULATE INCREMENTAL STRAINS
C     2. CALCULATE INCREMENTAL STRESSES,ASSUMING AN ELASTIC BEHAVIOR
C     3. CALCULATE TOTAL STRESSES,ASSUMING AN ELASTIC BEHAVIOR
C     --------------------------------------------------------------
C	FOR THIS MOMENT LET CONSIDER HP=HP1
	HP1=HP
	DFF = 1.E-6

	IF(YLD1.EQ.0.0D0) YLD1 = YLD

C	CALCULATE HERE INCREMENTAL STRAIN
      DO 110  I=1,ISR
 110  DELEPS(I) = STRAIN(I)-EPS(I)
C	CALCULATE HERE INCREMENTAL STRESS USING INCREMENTAL ELASTIC STRESS STRAIN RELATION
      CALL MEMSIG (DELEPS,DELSIG,IND)
      STRESS(4) = 0.
C	CALCULATE TRIAL STRESS USING STRESS AT PRIVIOUS UPDATE AND STRESS INCREMENT.
      DO 300  I=1,IST
 300  STRESS(I) = SIG(I) + DELSIG(I)
C     ----------------------------------------------------------------
C     4. CHECK WHETHER STATE OF STRESS FALLS OUTSIDE THE YIELD SURFACE
C     ----------------------------------------------------------------
C	CALCULATE DEVIATORIC STERSS OF TRIAL STERSS
      SM = (STRESS(1) + STRESS(2) + STRESS(4)) / 3.
      SX = STRESS(1)-SM
      SY = STRESS(2)-SM
      SS = STRESS(3)
      SZ = STRESS(4)-SM
      FT = 0.5*(SX*SX + SY*SY + SZ*SZ) + SS*SS - YLD1*YLD1/3.
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
C	SINCE STRESS FALLS OUTSIDE ADMISIBLE RANGE, FIND THE INTERSECTION POINT WHERE STRESS
C	ON THE YEILD SURFACE 
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
      E = SX*SX + SY*SY + 2.*SS*SS + SZ*SZ - 2.*YLD1*YLD1/3.
	SQQ = B*B-A*E
	SQA = ABS(SQQ)
      RATIO = (-B + SQRT(SQA)) / A

C
C	CALCULATE STRESS INCREMENT PURELY DUE TO ELASTIC BEHAVIOUR AND THEN 
C	CALCULATE TITAL STERSS PUELY DUE TO ELASTIC BEHAVIOUR
      DO 480  I=1,IST
 480  STRESS(I) = SIG(I) + RATIO*DELSIG(I)
      IF (IND.EQ.2)  STRAIN(4) = EPS(4) + RATIO*DELEPS(4)
C     --------------------------------------------------
C     5. DETERMINE NUMBER AND MAGNITUDE OF SUBINCREMENTS
C     --------------------------------------------------
 500	FFTT = SQRT(FT)/YLD1
c	NFF = INT(FFTT)
c	NSINC = 20*NFF + 1
      NSINC = 20.*SQRT(FT)/YLD1 + 1
      IF (NSINC.GT.30)  NSINC = 30
	IF (NSINC.LT.30)  NSINC = 5
      FACT = (1.0-RATIO) / NSINC
      DO 510  I=1,ISR
 510  DEPS(I) = FACT*DELEPS(I)
C     ----------------------------------------
C     6. CALCULATION OF ELASTOPLASTIC STRESSES
C     ----------------------------------------
C	START SUB INCREMENT 
      DO 695  ISINC=1,NSINC

      CALL MEDELPH (STRESS,STRAIN,DEPS,DP,IND,YLD1)
C	HERE STERSS IS TOTAL STRESS
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
      FTA= 1.5*(SX*SX + SY*SY + SZ*SZ + 2*SS*SS)
      FTB= YLD1*YLD1

      IF (ABS(FTA-FTB).LE.DFF)  GOTO 690
C
      IF (IND.GE.2)  GOTO 650
C	PLAIN STRAIN
      COEF =-(1.5*C2/(1.5*C2+HP1))*(1.-SQRT(FTB/FTA))
      STRESS(1) = STRESS(1) + COEF*SX
      STRESS(2) = STRESS(2) + COEF*SY
      STRESS(3) = STRESS(3) + COEF*SS
      STRESS(4) = STRESS(4) + COEF*SZ
C	CALCULATE LAMBDA
	SM = (STRESS(1) + STRESS(2) + STRESS(4)) / 3.
      SX = STRESS(1)-SM
      SY = STRESS(2)-SM
      SS = STRESS(3)
      SZ = STRESS(4)-SM

      GOTO 690
C
C	RETURN MAPPING ALGORITHUM FOR PLAIN STRESS

 650  COEF =-(1.5*C2/(1.5*C2+HP1))*(1.-SQRT(FTB/FTA))
      STRESS(1) = STRESS(1) + COEF*SX
      STRESS(2) = STRESS(2) + COEF*SY
      STRESS(3) = STRESS(3) + COEF*SS
	STRAIN(4) = STRAIN(4) + (COEF-1.)*SM/BM

C
C	CALCULATE LAMBDA
	SM = (STRESS(1) + STRESS(2) + STRESS(4)) / 3.
      SX = STRESS(1)-SM
      SY = STRESS(2)-SM
      SS = STRESS(3)
      SZ = STRESS(4)-SM

 690  CONTINUE
	
C	YLD1 HAS TO PASS AND STORE INSIDE A ARRAY
C	AND INITIALIZE IT TO YLD
      FTA1= 1.5*(SX*SX + SY*SY + SZ*SZ + 2*SS*SS)
	SIDEP=SX*DEPS(1)+SY*DEPS(2)+SS*DEPS(3)+SZ*DEPS(4)

	PLAM=(1.5*C2/(1.5*C2+HP1))*SIDEP/SQRT(FTA1)
	YLD1=YLD1+HP1*PLAM
 

695	CONTINUE
C	END OF SUB INCREMENTS

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
      IF (IPEL.EQ.2)  CALL MEDELPH (STRESS,STRAIN,DEPS,DP,IND,YLD1)
C
C
      RETURN
      END
C
C	=====================================================================
      SUBROUTINE MEDELPH (STRESS,STRAIN,DEPS,DP,IND,YLD1)
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
	COMMON /HARD/  HP,DEN
C
      DIMENSION STRESS(4),STRAIN(4),DEPS(4),DP(36)
C	FOR THIS MOMENT LET CONSIDER HP=HP1
	HP1=HP
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
C	SIGE=1.5*(SX*SX+SY*SY+SZ*SZ+2*SS*SS) 
c	AND APPLYING YIELD CRITERIA
C	YLD HAS TO UPDATE
      SIGE = YLD1*YLD1
	BETA1=9*C2*C2/4/(1.5*C2+HP1)
C     -----------------------------
C     CHECK FOR UNLOADING (LAMDA<0)
C     -----------------------------
C	CHECK FOR STABILITY
      IF (IND-1)  50,20,30
 20   WP = SX*DEPS(1) + SY*DEPS(2) + SS*DEPS(3)
      GOTO 100
C
C	CALCULATE INCREMENTAL STRAIN IN Z DIRECTION FOR PLAIN STRESS
 30   CBETSZ = BETA1*SZ/SIGE
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
C	CALCULATE ELASTO PLASTIC MODULUS MATRIX
      C2BETA = BETA1/SIGE
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
C	==========================================================
C	END OF VON MISS HARDENNING ROUTINES
C	==========================================================
C
C	STARTING OF IVANOVES HARDENNING ROUTINES
C	PLASTIC CURVATURE UPDATE IS INEFFECTIVE AT THE MOMENT
C
C	==========================================================
      SUBROUTINE IVANOVH (SIG,EPS,IPEL,STRAIN,STRESS,DP,DVOL,YLD2,EPC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CONVERTS ELASTO-PLASTIC STRAINS INTO STRESSES USING FLOW RULE
C     AND NORMALITY RULE (BASED ON IVANOVS YIELD CRITERION)
C	-----------------------------------------------------
C     INPUT VARIABLES
C	---------------
C     SIG(8)    = STRESS RESULTANTS AT THE END OF THE PREVIOUS UPDATE
C     EPS(8)    = STRAINS AT THE END OF THE PREVIOUS UPDATE
C     IPEL      = PLASTICITY FLAG (1=ELASTIC,2=PLASTIC)
C     STRAIN(8) = CURRENT TOTAL STRAINS
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     STRESS(8) = CURRENT TOTAL STRESS RESULTANTS
C     DP(8,8) = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DELEPS(8) = INCREMENT IN STRAINS
C     DELSIG(8) = INCREMENT IN STRESS RESULTANTS
C     DEPS(8)   = SUBINCREMENT IN STRAINS (EQUIVALENCED WITH DELEPS)
C     DEPSE(8)  = ELASTIC PART OF STRAIN INCREMENTS
C     DF(6)     = YIELD FUNCTION DERIVATIVES DF/DNX,DF/DNY,DF/DNXY,..
C     RATIO     = PART OF STRESS TAKEN ELASTICLY
C     NSINC     = NUMBER OF SUBINCREMENTS
C     PLAMDA    = PLASTIC STRAIN RATE MULTIPLIER
C     /HOOK/    = FOR VARIABLES IN COMMON BLOCK HOOK SEE CONSTI
C     ----------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION SIG(8),EPS(8),STRAIN(8),STRESS(8),DP(8,8)
      DIMENSION DELEPS(8),DELSIG(8),DEPS(8),DEPSE(8),DF(6)
      EQUIVALENCE (DELEPS(1),DEPS(1))
C     --------------------------------------------------------------
C     1. CALCULATE INCREMENTAL STRAINS
C     2. CALCULATE INCREMENTAL STRESS RES. ASSUMING ELASTIC BEHAVIOR
C     3. CALCULATE TOTAL STRESS RESULTANTS,ASSUMING ELASTIC BEHAVIOR
C     --------------------------------------------------------------

!	HAS TO PASS EPC AND YLD2

	IF(YLD2.EQ.0.0D0) YLD2 = YLD

	EPC = 0.
      DO 110  I=1,8
 110  DELEPS(I) = STRAIN(I)-EPS(I)
      CALL MPSIGA (DELEPS,DELSIG)
      DO 300  I=1,8
 300  STRESS(I) = SIG(I) + DELSIG(I)
C     ----------------------------------------------------------------
C     4. CHECK WHETHER STATE OF STRESS FALLS OUTSIDE THE YIELD SURFACE
C     ----------------------------------------------------------------
      CALL YIELDFH (STRESS(1),STRESS(2),STRESS(3),STRESS(4),STRESS(5),
     1 STRESS(6),UYN,UYM,UYN2,UYM2,QT,QM,QTM,Q,R,S,FT,2,YLD2,EPC,GAMA)
      IF (FT)  410,410,450
C     -------------------------------------------------------
C     STATE OF STRESS WITHIN YIELD SURFACE - ELASTIC BEHAVIOR
C     -------------------------------------------------------
 410  IPEL = 1
      GOTO 700
C     ------------------------------------------------------------
C     STATE OF STRESS OUTSIDE YIELD SURFACE - PLASTIC BEHAVIOR
C     DETERMINE PART OF STRAIN TAKEN ELASTICLY AND ADD STRESSES
C     DUE TO THE ELASTIC STRAIN INCREMENT TO THE PREVIOUS STRESSES
C     ------------------------------------------------------------
 450  IPEL = 2
      CALL YIELDFH (SIG(1),SIG(2),SIG(3),SIG(4),SIG(5),SIG(6),
     1     UYN,UYM,UYN2,UYM2,QT,QM,QTM,Q,R,S,YFT,2,YLD2,EPC,GAMA)
      IF (YFT.LT.0)  GOTO 470
      RATIO = 0.
      DO 460  I=1,8
 460  STRESS(I) = SIG(I)
      GOTO 500
 470  CALL TRANSIH (SIG(1),SIG(2),SIG(3),SIG(4),SIG(5),SIG(6),
     1             DELSIG(1),DELSIG(2),DELSIG(3),DELSIG(4),DELSIG(5),
     2             DELSIG(6),UYN,UYM,UYN2,UYM2,RATIO,YLD2)
      DO 480  I=1,8
 480  STRESS(I) = SIG(I) + RATIO*DELSIG(I)
C     --------------------------------------------------
C     5. DETERMINE NUMBER AND MAGNITUDE OF SUBINCREMENTS
C     --------------------------------------------------
 500  NSINC = 20.*SQRT(FT) + 1
      IF (NSINC.GT.30)  NSINC = 30
	IF (NSINC.LT.5)   NSINC = 5
      FACT = (1.0-RATIO) / NSINC
      DO 510  I=1,8
 510  DEPS(I) = FACT*DELEPS(I)
C     ----------------------------------------
C     6. CALCULATION OF ELASTOPLASTIC STRESSES 
C     ----------------------------------------

C	STARTING SUB ITERATION........

      DO 690  ISINC=1,NSINC
      CALL LAMDAPH (STRESS(1),STRESS(2),STRESS(3),STRESS(4),STRESS(5),
     1          STRESS(6),DEPS,PEQSTR,DF,DENOM,YLD2,EPC,GAMA,PLAMDA)
      IF (PEQSTR.LT.0)  PEQSTR = 0.
	IF (PLAMDA.LT.0)  PLAMDA = 0.
C
      DO 610  I=1,6
 610  DEPSE(I) = DEPS(I) - PLAMDA*DF(I)
      DEPSE(7) = DEPS(7)
      DEPSE(8) = DEPS(8)
      CALL MPSIGA (DEPSE,DELSIG)
      DO 620  I=1,8
 620  STRESS(I) = STRESS(I) + DELSIG(I)

C     FORCE STRESS STATE BACK ON YIELD SURFACE

 690  CONTINUE

C     -----------------------------------------
C     7. UPDATING STRESS RESULTANTS AND STRAINS
C     -----------------------------------------
 700  DO 710  I=1,8
      SIG(I) = STRESS(I)
 710  EPS(I) = STRAIN(I)
C     -------------------------------------------
C     8. FORM THE MATERIAL LAW (MATRIX DP) IF THE
C        STIFNESS IS TO BE REFORMED (IFREF=0)
C     -------------------------------------------
      IF (IFREF.NE.0) RETURN
      IF (IPEL.EQ.1) RETURN
      CALL LAMDAPH (STRESS(1),STRESS(2),STRESS(3),STRESS(4),STRESS(5),
     1             STRESS(6),DEPS,PLAMDA,DF,DENOM,YLD2,EPC,GAMA,PLAMDA)
      IF (PLAMDA.GT.0.) CALL MPDELPH (DF,DVOL/DENOM,DP,DP)
C
      RETURN
      END
C
C	=====================================================================
      SUBROUTINE YIELDFH (NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     1                   QT,QM,QTM,Q,R,S,FT,IND,YLD2,EPC,GAMA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     RETURNS CURRENT VALUE OF IVANOVS YIELD FUNCTION
C	-----------------------------------------------
C     NX,NY,NXY,MX,MY,MXY= STRESS RESULTANTS
C     UYN,UYM,UYN2,UYM2  = UNIAXIAL YIELD FORCE/MOMENT PER UNIT WIDTH
C     QT,QM,QTM          = QUADRATIC NON-DIMENSIONAL STRESS INTENSITIE
C     Q,R,S              = COMPONENTS IN IVANOVS YIELD FUNCTION
C     FT                 = VALUE OF IVANOV'S YIELD FUNCTION
C     IND                = INDICATES WHETHER FT IS REQUIRED (IND=2)
C	GAMA			   = CORRECTION FACTOR
C	EPC				   = EQUIVALENT PLASTIC CURVATURE
C	YLD2			   = UPDATED YIELD STRESS
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      REAL*8 NX,NY,NXY,MX,MY,MXY
C
!	HAS TO PASS EPC

C	CALCULATE GAMA HERE USING EQUIVALANT PLASTIC CURVATURE
	GAMA = 1.-0.4*EXP(-2.6*SQRT(EPC))
	IF(EPC.LE.1.E-6) GAMA = 1.

      UYN  = TH
      UYM  = GAMA*TH*TH/4.
      UYN2 = UYN*UYN
      UYM2 = UYM*UYM
C
      QT =(NX*NX + NY*NY - NX*NY + 3.*NXY*NXY) / UYN2
      QM =(MX*MX + MY*MY - MX*MY + 3.*MXY*MXY) / UYM2
      QTM=(MX*NX + MY*NY - MX*NY/2. - MY*NX/2. + 3.*MXY*NXY)/(UYN*UYM)
C

!	HAS TO PASS YLD2
      ENTRY FTVALUH(NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     +             QT,QM,QTM,Q,R,S,FT,IND,YLD2)
      Q   = QT + .48*QM
      R   = SQRT(QM*QM/4. + QTM*QTM)
      S   = QT*QM - QTM*QTM
C
      IF (IND.EQ.1)  RETURN
      FT = -YLD2*YLD2
      IF (Q.EQ.0.)  RETURN
      FT = QT + QM/2. + R - .25*S/Q - YLD2*YLD2
C
      RETURN
      END
C
C	=====================================================================
      SUBROUTINE TRANSIH (NX,NY,NXY,MX,MY,MXY,DNX,DNY,DNXY,
     1                   DMX,DMY,DMXY,UYN,UYM,UYN2,UYM2,RATIO,YLD2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     RETURNS RATIO OF STRAIN INCREMENTS TAKEN ELASTICLY
C	--------------------------------------------------
C     NX,NY,NXY,MX,MY,MXY      = STRESS RESULTANTS AT PREVIOUS UPDATE
C     DNX,DNY,DNXY,DMX,DMY,DMXY= INCREMENTS IN STRESS RESULTANTS
C     UYN,UYM,UYN2,UYM2        = UNIAXIAL YIELD FORCE/UNIT WIDTH
C     RATIO                    = PART OF STRAIN INCR. TAKEN ELASTICLY
C     ----------------------------------------------------------------
	
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C
      REAL*8 NX,NY,NXY,MX,MY,MXY
C     -------------------------------------
C     QUADRATIC IN PLANE STRESS INTENSITIES
C     N = AM*RATIO**2 + BM*RATIO + CM
C     -------------------------------------
	IF (ITOLEY .EQ.0) ITOLEY = 100000
      AM = DNX*DNX + DNY*DNY - DNX*DNY + 3.*DNXY*DNXY
      BM = 2.*(NX*DNX + NY*DNY) - NX*DNY - DNX*NY + 6.*NXY*DNXY
      CM = NX*NX + NY*NY - NX*NY + 3.*NXY*NXY
C     ------------------------------------
C     QUADRATIC BENDING STRESS INTENSITIES
C     M = AB*RATIO**2 + BB*RATIO + CB
C     ------------------------------------
      AB = DMX*DMX + DMY*DMY - DMX*DMY + 3.*DMXY*DMXY
      BB = 2.*(MX*DMX + MY*DMY) - MX*DMY - DMX*MY + 6.*MXY*DMXY
      CB = MX*MX + MY*MY - MX*MY + 3.*MXY*MXY
C     --------------------------------------------------------
C     COMPUTE RATIO IN CASE OF PURE STRETCHING OR PURE BENDING
C     --------------------------------------------------------
      IF (AM+CM.NE.0. .AND. AB+CB.NE.0.)  GOTO 200
	UYN2 = UYN2*YLD2*YLD2
	UYM2 = UYM2*YLD2*YLD2

      RATIO = -BM-BB+SQRT(BM*BM-4.*AM*(CM-UYN2)+BB*BB-4.*AB*(CB-UYM2))
      RATIO = RATIO/2./(AM+AB)
      RETURN
C     -----------------------------------
C     MIXED QUADRATIC STRESS INTENSITIES
C     MN = AMB*RATIO**2 + BMB*RATIO + CMB
C     -----------------------------------
 200  AMB = DMX*DNX + DMY*DNY - .5*(DMX*DNY+DMY*DNX) + 3.*DMXY*DNXY
      BMB = MX*DNX + DMX*NX + MY*DNY + DMY*NY - .5*(MX*DNY+DMX*NY+
     1      MY*DNX+DMY*NX) + 3.*(MXY*DNXY+DMXY*NXY)
      CMB = MX*NX + MY*NY - .5*(MX*NY+MY*NX) + 3.*MXY*NXY
C     ---------------------------------------
C     NUMERICAL LOOP TO FIND RATIO FOR FT = 0
C     ---------------------------------------
      GAP = 1.0
      FAC = 1.0
      RATIO = 0.0
C
 300  GAP    = GAP/2.
      RATIO  = RATIO + FAC*GAP
      RATIO2 = RATIO*RATIO
      FAC    = 1.
C
      QT  = (AM *RATIO2 + BM *RATIO + CM) / UYN2
      QM  = (AB *RATIO2 + BB *RATIO + CB) / UYM2
      QTM = (AMB*RATIO2 + BMB*RATIO + CMB) / (UYN*UYM)

!	HAS TO PASS YLD2

      CALL FTVALUH (NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     1             QT,QM,QTM,Q,R,S,FT,2,YLD2)
C
      IF (FT.GT.0.)  FAC = -1.
	ABCD=1./1000000.
      IF (ABS(FT)-ABCD)  400,400,300
C
  400 RETURN
      END
C
C	=====================================================================
      SUBROUTINE LAMDAPH (NX,NY,NXY,MX,MY,MXY,DEPS,PEQSTR,DF,B,
	1					YLD2,EPC,GAMA,PLAMDA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     RETURNS THE PLASTIC STRAIN RATE MULTIPLIER AND THE YIELD
C     FUNCTION DERIVATIVES
C	--------------------
C     NX,NY,NXY,MX,MY,MXY  = STRESS RESULTANTS
C     DEPS(8)              = SUBINCREMENT OF STRAIN
C     DF(6)                = YIELD FUNCTION DERIVATIVES
C     PLAMDA               = PLASTIC STRAIN RATE MULTIPLIER
C     B                    = DENOMINATOR OF PLAMDA
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
	COMMON /HARD/  HP,DEN
      DIMENSION DEPS(8),DQ(6),DF(6),DFD(6)
C
      REAL*8 NX,NY,NXY,MX,MY,MXY
C     -----------------------------------------------------
C     UNIAXIAL YIELD FORCE AND MOMENT, QUADRATIC STRESS
C     INTENSITIES N,M,NM, CONSTANTS Q,R,S OF YIELD FUNCTION
C     -----------------------------------------------------
      CALL YIELDFH (NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     1             QT,QM,QTM,Q,R,S,FT,1,YLD2,EPC,GAMA)
C     -----------------------------------------------------
C     CHECK FOR DISCONTINUITY OF YIELD FUNCTION DERIVATIVES
C     -----------------------------------------------------
      HR = 0.
      IF (R.GT.1.E-04)  HR = 1./R
C
      Q2 = Q*Q
      Q4 = 4.*Q
C	D = FIRST COIFICEANT OF PARTIAL DERIVATIVE OF YIELD FUNCTION W.R.T. NORMAL FORCE
      C = 1. - QM/Q4 + S/(4.*Q2)

C	D = FIRST COIFICEANT OF PARTIAL DERIVATIVE OF YIELD FUNCTION W.R.T. MOMENT
      D = .5 + HR*QM/4. - QT/Q4 + .12*S/Q2

      UYNM = UYM/UYN
      E = QTM * (HR+.5/Q)
	E1 = E
      F = E/UYNM
      E = E*UYNM
C     -------------------------------------------
C     VECTOR OF DERIVATIVES OF STRESS INTENSITIES
C     -------------------------------------------
      DQ(1) = (2.*NX-NY) / UYN2
      DQ(2) = (2.*NY-NX) / UYN2
      DQ(3) = 6.*NXY/UYN2
      DQ(4) = (MX-0.5*MY) / UYM2
      DQ(5) = (MY-0.5*MX) / UYM2
      DQ(6) = 3.*MXY/UYM2
C     ------------------------------------
C     VECTOR OF YIELD FUNCTION DERIVATIVES
C     ------------------------------------
      DO 100  I=1,3
      DF(I)   = C*DQ(I) + E*DQ(I+3)
 100  DF(I+3) = F*DQ(I) + D*DQ(I+3)

!	GAMA HAS TO PASS
C	--------------------------------------------
C	DERIVATIVES OF STRESS INTENSITIES W.R.T GAMA
C	--------------------------------------------
	DQM = -32.*(MX*MX + MY*MY - MX*MY + 3.*MXY*MXY)/GAMA**3/TH**4	
	DQTM = 4.*(MX*NX + MY*NY - MX*NY/2. - MY*NX/2. + 3.*MXY*NXY)
	1					/GAMA**2/TH**3
C	--------------------------------------------
C	DERIVATIVE OF YIELD FUNCTION W.R.T GAMA	
C	---------------------------------------------
	DFG = D*DQM+E1*DQTM

!	HAS TO PASS EPC
C	-----------------------------------------------
C	DERIVATIVE OF COIFICIENT W.R.T EPC
C	-------------------------------------------------
	GAEP = 0.52*TH*YM*EXP(-2.6*SQRT(EPC))/3./YLD/SQRT(EPC)
	IF(EPC.LE.1.E-6) GAEP = 0.
C	---------------------------------------------
C	CALCULATE B
C	------------------------------------------
	BXL = 2./SQRT(3.)*SQRT(DF(4)*DF(4)+DF(5)*DF(5)+
	1						DF(4)*DF(5)+DF(6)*DF(6)/4.)		

!	HAS TO ADD HARDDENING PARAMETER HP (COMMON BLOCK FOR HARDENNING)
!	HAS TO PASS YLD2
C	-------------------------------------------------------
C	CALCULATE HERE ADDITINAL PART DUE TO HARDENNING
C	----------------------------------------------------------
	ECOIF = NX*DF(1)+NY*DF(2)+NXY*DF(3)+
	1				MX*DF(4)+MY*DF(5)+MXY*DF(6)

	AL = 2.*HP*ECOIF
	BAL = AL-BXL*DFG*GAEP

C     ------------------------------
C     PLASTIC STRAIN RATE MULTIPLIER
C     ------------------------------
      TH3 = TH*TH*TH/12.
      DFD(1) = (A1*DF(1) + B1*DF(2)) * TH
      DFD(2) = (B1*DF(1) + A1*DF(2)) * TH
      DFD(3) = C1*DF(3)*TH
      DFD(4) = (A1*DF(4) + B1*DF(5)) * TH3
      DFD(5) = (B1*DF(4) + A1*DF(5)) * TH3
      DFD(6) = C1*DF(6)*TH3
C


      A = 0.
      B = 0.
      DO 200  I=1,6
      A = A + DFD(I)*DEPS(I)
 200  B = B + DFD(I)*DF(I)

	B = B + BAL
	PLAMDA = A/B
      PEQSTR = ECOIF*PLAMDA/YLD2
C
C	-----------------------------------------------------------
C	UPDATE HERE YIELD FUNCTION AND EQUIVALANT PLASTIC CURVATURE
C	-----------------------------------------------------------
	YLD2 = YLD2+HP*PEQSTR
C	EPC  = EPC + BXL*PLAMDA

      RETURN
      END
C
C	=====================================================================
      SUBROUTINE MPDELPH (DF,FAC,DP,CP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     FORMS THE ELASTO-PLASTIC STRESS - STRAIN MATRIX
C	-----------------------------------------------
C     DF(6)    = YIELD FUNCTION DERIVATIVES
C     FAC      = DVOL/DENOMINATOR OF STRAIN RATE MULTIPLIER
C     DP(64)   = STRESS - STRAIN MATRIX
C     CP(8,8)  = DP(64)
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION DF(6),DP(64),CP(8,8),FFT(6,6),FD(6)
C     -----------------------------------------------
C     1. PRODUCT MATRIX OF YIELD FUNCTION DERIVATIVES
C        FFT(I,J) = DF(I)*DF(J)
C     -----------------------------------------------
      DO 50  I=1,6
      DO 50  J=I,6
      FFT(I,J) = DF(I)*DF(J)
 50   FFT(J,I) = FFT(I,J)
C     ---------------------------------------------------------
C     2. MATRIX MULTIPLICATION [DP] = [DP] - [D]*[FFT]*[D]*FACT
C     ---------------------------------------------------------
      TH3   = TH*TH*TH/12.
      FACM  = TH*TH*FAC
      FACB  = TH3*TH3*FAC
      FACMB = TH*TH3*FAC
C
      DO 100  I=1,6
 100  FD(I)  = FFT(I,1)*A1 + FFT(I,2)*B1
      DP(1)  = DP(1)  - (A1*FD(1) + B1*FD(2)) * FACM
      DP(2)  = DP(2)  - (B1*FD(1) + A1*FD(2)) * FACM
      DP(3)  =        -  C1*FD(3)*FACM
      DP(4)  =        - (A1*FD(4) + B1*FD(5)) * FACMB
      DP(5)  =        - (B1*FD(4) + A1*FD(5)) * FACMB
      DP(6)  =        -  C1*FD(6)*FACMB
C
      DO 200  I=1,6
 200  FD(I)  = FFT(I,1)*B1 + FFT(I,2)*A1
      DP(10) = DP(10) - (B1*FD(1) + A1*FD(2)) * FACM
      DP(11) =        -  C1*FD(3)*FACM
      DP(12) =        - (A1*FD(4) + B1*FD(5)) * FACMB
      DP(13) =        - (B1*FD(4) + A1*FD(5)) * FACMB
      DP(14) =        -  C1*FD(6)*FACMB
C
      DO 300  I=3,6
 300  FD(I)  = FFT(I,3)*C1
      DP(19) = DP(19) -  C1*FD(3)*FACM
      DP(20) =        - (A1*FD(4) + B1*FD(5)) * FACMB
      DP(21) =        - (B1*FD(4) + A1*FD(5)) * FACMB
      DP(22) =        -  C1*FD(6)*FACMB
C
      DO 400  I=4,6
 400  FD(I)  = FFT(I,4)*A1 + FFT(I,5)*B1
      DP(28) = DP(28) - (A1*FD(4) + B1*FD(5)) * FACB
      DP(29) = DP(29) - (B1*FD(4) + A1*FD(5)) * FACB
      DP(30) =        -  C1*FD(6)*FACB
C
      DO 500  I=4,6
 500  FD(I)  = FFT(I,4)*B1 + FFT(I,5)*A1
      DP(37) = DP(37) - (B1*FD(4) + A1*FD(5)) * FACB
      DP(38) =        -  C1*FD(6)*FACB
C
      FD(6)  = FFT(6,6)*C1
      DP(46) = DP(46) -  C1*FD(6)*FACB
C     -----------------------------------------------
C     3. FILL IN SYMMETRIC ELEMENTS OF UPPER TRIANGLE
C     -----------------------------------------------
 900  DO 950  I=1,6
      DO 950  J=I,6
 950  CP(I,J) = CP(J,I)
C
      RETURN
      END
C
C	==============================================================
C	END OF IVANOVES HARDENNING ROUTINES
C	==============================================================
