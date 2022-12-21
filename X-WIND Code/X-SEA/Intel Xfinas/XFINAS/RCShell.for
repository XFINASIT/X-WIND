C	==================================================================
C	========= EPF MODEL REINFORCED CONCRETE SHELL ELEMENT ============
C	==================================================================
	SUBROUTINE RCLAYRE(WA,STNCR,STSCR,IELE,TH,DP,DVOL,MSET,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	PRODUCED BY SONGSAK - FEB2005
C	==================================================================
C	COMPUTE THE GAUSS POINT STRESSES AND MATERIAL RIGIDITY FOR
C	REINFORCED CONCRETE SHELL ELEMENT USING LAYERED MODEL
C	==================================================================
C	INPUT VARIABLE
C	==============
C	STNCR(8)     = CURRENT GAUSS POINT STRAIN (MEMBRANE-BENDING-SHEAR)
C	TH			 = GAUSS POINT THICKNESS
C	DVOL		 = GAUSS POINT INTEGRATION VOLUME
C
C	==================================================================
C	MSTAT    = MATERIAL NUMBER
C				 1 = ELASTIC BEHAVIOR CONCRETE
C				 2 = CRACKED CONCRETE (CRACK IN ONE DIRECTION)
C				 3 = CRACKED CONCRETE (CRACK IN TWO DIRECTION)
C				 4 = ELASTIC BEHAVIOR OF STEEL
C				 5 = ELASTO-PLASTIC BEHAVIOR STEEL
C				 6 = CRUSHING OF CONCRETE
C	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (C)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C	
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
C	=============================================================
	DIMENSION WA(43,1),STSPS(5),STNPS(5),STNCR(8),STSCR(8),
	1		  DMATX(5,5),STNZ(5),DSTNZ(5),DSTSZ(5),SIGMA(5),
     2		  DIREC(2),EPSTN(5),GRTSN(2),DP(8,8),CP(5,5),
     3		  WOKPR(3),EQIST(2),EQISN(3),ELSNP(5),EQCRK(6)


C	WRITE(*,*) 'KAK'

	VOLM  = DVOL
	LLAYR = MSET
	NLAYR = MATLR(MSET,MNLYR+1)  !MNLYR
	ZETAD  = -1.0

	DO 5 I = 1,8
	STSCR(I) = 0.0
	DO 5 J = 1,8
5	DP(I,J)  = 0.0
	
C	------------------------------------------------------------------
C	LOOP OVER CONCRETE & STEEL LAYER
C	------------------------------------------------------------------
      DO 200 ILAYR = 1,NLAYR
	LPROP = MATLR(LLAYR,ILAYR)
      IDTCS = PROPS(LPROP,10)

	DO 8 I = 1,5
	J = I+5
	STNPS(I) = WA(I,ILAYR)
8	STSPS(I) = WA(J,ILAYR)

	DIREC(1) = WA(11,ILAYR)
	DIREC(2) = WA(12,ILAYR)

	EPSTN(1) = WA(13,ILAYR)
	EPSTN(2) = WA(14,ILAYR)
	EPSTN(3) = WA(15,ILAYR)
	EPSTN(4) = WA(16,ILAYR)
	EPSTN(5) = WA(17,ILAYR)

	ELSNP(1) = WA(18,ILAYR)
	ELSNP(2) = WA(19,ILAYR)
	ELSNP(3) = WA(20,ILAYR)
	ELSNP(4) = WA(21,ILAYR)
	ELSNP(5) = WA(22,ILAYR)

	GRTSN(1) = WA(23,ILAYR)
	GRTSN(2) = WA(24,ILAYR)

	EQIST(1) = WA(25,ILAYR)
	EQIST(2) = WA(26,ILAYR)

	EQISN(1) = WA(27,ILAYR)
	EQISN(2) = WA(28,ILAYR)
	EQISN(3) = WA(29,ILAYR)

	MSTAT    = WA(30,ILAYR)

	DEO2     = WA(31,ILAYR)
	DEMX2    = WA(32,ILAYR)

	WOKPR(1) = WA(33,ILAYR)
	WOKPR(2) = WA(34,ILAYR)
	WOKPR(3) = WA(35,ILAYR)

	EQCRK(1) = WA(36,ILAYR)
	EQCRK(2) = WA(37,ILAYR)
	EQCRK(3) = WA(38,ILAYR)
	EQCRK(4) = WA(39,ILAYR)
	EQCRK(5) = WA(40,ILAYR)
	EQCRK(6) = WA(41,ILAYR)

	FLAG_ANG = WA(42,ILAYR)
	IF(FLAG_ANG.NE.0.0) REFANG = WA(43,ILAYR)

	IF(IDTCS.EQ.0) THEN
	DZETA = PROPS(LPROP,3)
	ZETA  = ZETAD + DZETA/2.0
	ZETAG = 0.5*TH*ZETA
	DZETG = 0.5*TH*DZETA
	ELSE
	DZETA = PROPS(LPROP,3)
	ZETA  = PROPS(LPROP,6)
	ZETAG = 0.5*TH*ZETA
	DZETG = 0.5*TH*DZETA
	ENDIF

C	------------------------------------------------------------------
C	COMPUTE STRAINS AT THE MIDDLE OF LAYER
C	------------------------------------------------------------------
	DO 10 I = 1,3
	STNZ(I)  = STNCR(I) + ZETAG*STNCR(I+3)
10	DSTNZ(I) = STNZ(I)  - STNPS(I)
	STNZ(4)  = STNCR(7)
	STNZ(5)  = STNCR(8)
	DSTNZ(4) = STNCR(7) - STNPS(4)
	DSTNZ(5) = STNCR(8) - STNPS(5)
	
	IF(IDTCS.EQ.1) DIREC(1) = PROPS(LPROP,7)
	IF(IDTCS.EQ.1) GO TO 100

C	------------------------------------------------------------------
C	COMPUTE STRESSES AT THE MIDDLE OF CONCRETE LAYER
C	------------------------------------------------------------------

	TENST = PROPS(LPROP,5)
	UNIAX = PROPS(LPROP,6)
	POISN = PROPS(LPROP,2)
	YOUNG = PROPS(LPROP,1)


	NSTAT = MSTAT
	ANGLE = DIREC(1)
	STRA1 = WOKPR(1)
	STRA2 = WOKPR(2)
	IF(NSTAT.EQ.1) THEN
	CALL	MODUL1(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,YOUNG,  0.0,
     2			  POISN)
	ELSE
	CALL	MODUL1(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,	0.0,  0.0,
     2			  0.0)
	ENDIF

	DO 15 I = 1,5
	DSTSZ(I) = 0.0
	DO 15 J = 1,5
15	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)

	DO 20 I = 1,5
	SIGMA(I) = 0.0
20	SIGMA(I) = STSPS(I) + DSTSZ(I)

C	------------------------------------------------------------------
C	COMPUTE MAX & MIN STRESS IN THE PRINCIPAL DIRECTION
C	------------------------------------------------------------------
	CALL	PRISTE(SIGMA,SGMAX,SGMIN,1)

C	------------------------------------------------------------------
C	CHECK FOR CRAKING OF CONCRETE
C	------------------------------------------------------------------
	FKO   = WOKPR(3)
	RF    = FKO*FKO*FKO
	CRIT  = 0.0
	CKDUM = SGMAX*SGMIN
	IF(CKDUM.LT.0.0) THEN
	CRIT = SGMAX/RF/TENST
	ENDIF
	IF(SGMAX.GT.0.0.AND.SGMIN.GT.0.0) THEN
	CRIT = (SGMAX/RF/TENST) + 0.26*(SGMIN/SGMAX)
	ENDIF
	IF(CRIT.GT.1.0) GO TO 30
C	IF(SGMAX.GT.TENST) GO TO 30
	IF(NSTAT.EQ.1) GO TO 40

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
30	CALL	CCRACKE(STSPS,NSTAT,SIGMA,DIREC,
	1			   EPSTN,GRTSN,SGMAX,LPROP,
	2			   STNZ,WOKPR)

	IF(NSTAT.EQ.3) GO TO 95
C	------------------------------------------------------------------
C	COMPRESION BEHAVIOR OF CONCRETE
C	------------------------------------------------------------------
40	CALL	CYIELDE(STSPS,NSTAT,SIGMA,DIREC,
	1			   EPSTN,LPROP,EQIST,EQISN,
	2			   STNZ,WOKPR,DEMX2,ELSNP,
	3			   EQCRK,KUNLO)



	GO TO 95

C	------------------------------------------------------------------
C	COMPUTE STRESS AT THE MIDDLE OF STEEL LAYER
C	------------------------------------------------------------------
100	NSTAT = MSTAT
	ANGLE = PROPS(LPROP,7)+REFANG

	CALL	MODUL1(DMATX,LPROP,	4,ANGLE,
	1			  0.0, 0.0, 0.0, 0.0, 0.0)

	DO 60 I = 1,5
	DSTSZ(I) = 0.0
	DO 60 J = 1,5
60	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)

	CALL	STELRE(STSPS,NSTAT,EPSTN,DSTSZ,LPROP,REFANG)

95	CONTINUE

C	------------------------------------------------------------------
C	UPDATE THE LAYER STRAIN
C	------------------------------------------------------------------
	DO 96 I = 1,5 
96	STNPS(I) = STNZ(I)
	
	MSTAT = NSTAT

C	------------------------------------------------------------------
C	STORE THE CURRENT STRESSES IN WA
C	------------------------------------------------------------------
	DO 97 I = 1,5
	J = I+5
	WA(I,ILAYR)  = STNPS(I)
97	WA(J,ILAYR)  = STSPS(I)

	WA(11,ILAYR) = DIREC(1)
	WA(12,ILAYR) = DIREC(2)

	WA(13,ILAYR) = EPSTN(1)
	WA(14,ILAYR) = EPSTN(2)
	WA(15,ILAYR) = EPSTN(3)
	WA(16,ILAYR) = EPSTN(4)
	WA(17,ILAYR) = EPSTN(5)

	WA(18,ILAYR) = ELSNP(1)
	WA(19,ILAYR) = ELSNP(2)
	WA(20,ILAYR) = ELSNP(3)
	WA(21,ILAYR) = ELSNP(4)
	WA(22,ILAYR) = ELSNP(5)

	WA(23,ILAYR) = GRTSN(1)
	WA(24,ILAYR) = GRTSN(2)

	WA(25,ILAYR) = EQIST(1)
	WA(26,ILAYR) = EQIST(2)

	WA(27,ILAYR) = EQISN(1)
	WA(28,ILAYR) = EQISN(2)
	WA(29,ILAYR) = EQISN(3)

	WA(30,ILAYR) = MSTAT

	WA(31,ILAYR) = DEO2
	WA(32,ILAYR) = DEMX2

	WA(33,ILAYR) = WOKPR(1)
	WA(34,ILAYR) = WOKPR(2)
	WA(35,ILAYR) = WOKPR(3)

	WA(36,ILAYR) = EQCRK(1)
	WA(37,ILAYR) = EQCRK(2)
	WA(38,ILAYR) = EQCRK(3)
	WA(39,ILAYR) = EQCRK(4)
	WA(40,ILAYR) = EQCRK(5)
	WA(41,ILAYR) = EQCRK(6)

	IF(FLAG_ANG.EQ.0.0) THEN
	WA(42,ILAYR) = 1.0
	WA(43,ILAYR) = REFANG 
	ENDIF
	
C	------------------------------------------------------------------
C	CONTRIBUTION OF STRESSES OVER THICKNESS
C	MEMBRANE - BENDING - SHEAR
C	------------------------------------------------------------------
	DO 150 I = 1,3
	J = I+3
	STSCR(I) = STSCR(I) + STSPS(I)*DZETG
150	STSCR(J) = STSCR(J) + STSPS(I)*DZETG*ZETAG
	STSCR(7) = STSCR(7) + STSPS(4)*DZETG
	STSCR(8) = STSCR(8) + STSPS(5)*DZETG

	CALL	CONSTVE(STSPS,GRTSN,NSTAT,EPSTN,
	1			   DIREC,LPROP,CP,KUNLO,WOKPR,REFANG)

C	------------------------------------------------------------------
C	INITIALIZE ELASTO-PLASTIC RIGIDITY MATRIX	
C	------------------------------------------------------------------
	FA = DZETG*VOLM
	FC = FA*ZETAG
	FB = FC*ZETAG
	FE = DZETG*DZETG*DZETG/12.0

C	------------------------------------------------------------------
C	MEMBRANE & BENDING RIGIDITY
C	------------------------------------------------------------------
	DO 300 I = 1,3
	DO 300 J = 1,3
	K = I+3
	L = J+3
	DP(I,J) = DP(I,J) + CP(I,J)*FA
300	DP(K,L) = DP(K,L) + CP(I,J)*(FB+FE)

C	------------------------------------------------------------------
C	SHEAR RIGIDITY
C	------------------------------------------------------------------
	DO 350 I = 4,5
	DO 350 J = 4,5
	K = I+3
	L = J+3
350	DP(K,L) = DP(K,L) + CP(I,J)*FA

C	------------------------------------------------------------------
C	COUPLE MEMBRANE AND BENDING RIGIDITY
C	(UPPER PART OF TRIANGLE)
C	------------------------------------------------------------------
	DO 380 I = 1,3
	DO 380 J = 1,3
	L = J+3
380	DP(I,L) = DP(I,L) + CP(I,J)*FC

C	------------------------------------------------------------------
C	FILL IN LOWER PART OF RIGIDITY MATRIX
C	------------------------------------------------------------------
	DO 400 I = 1,6
	DO 400 J = 1,6
400	DP(J,I) = DP(I,J)

	IF(IDTCS.EQ.0) THEN
	ZETAD  = ZETA + DZETA/2.0
	ENDIF

200	CONTINUE


	RETURN
	END


C	==================================================================

	SUBROUTINE CYIELDE(STSPS,NSTAT,SIGMA,DIREC,
	1				  EPSTN,LPROP,EQIST,EQISN,
	2			      STNZ,WOKPR,DEMX2,ELSNP,
	3			      EQCRK,KUNLO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	COMPRESSIVE BEHAVIOUR OF CONCRETE (EPF MODEL)
C	=============================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
	DIMENSION STSPS(5),SIGMA(5),DMATX(5,5),STNZ(5),DSTNZ(5),
	1		  DSTSZ(5),DIREC(2),EPSTN(5),STRSP(5),STRNP(5),
     2		  SGTOT(5),TRAMX(5,5),WOKPR(3),ELSTN(5),EQIST(2),
	3		  EQISN(3),ELSNP(5),DELSN(5),EQCRK(6),DPSTN(5)
C	==================================================================
	
	YOUNG = PROPS(LPROP,1)
	POISN = PROPS(LPROP,2)
	UNIAX = PROPS(LPROP,6)
	UNSTN = PROPS(LPROP,7)
	UNSTN = -UNSTN ! SONGSAK MAR2005
	ANGL1 = DIREC(1)
	GAMMA = 1.570796327
	IF(DIREC(1).GE.0.0) GAMMA = -GAMMA
	ANGL2 = ANGL1 + GAMMA

C	==================================================================
C	INITIALIZE THE CONSTANTS
C	==================================================================
	EOPIM = 2.0*UNIAX/YOUNG
	EO    = 2.0
	PN    = 2.0

	SMAX  = EQIST(1)
	SO    = EQIST(2)
	EMAX  = EQISN(1)
	EEMAX = EQISN(2)
	EEO   = EQISN(3)

C	==================================================================
C	CHECK FOR CRUSHING OF CONCRETE
C	==================================================================	
	IF(NSTAT.EQ.1) THEN
	CALL INVA2E(STNZ,ESN)
	IF(ESN.GT.ABS(UNSTN)) NSTAT = 6
	ELSE
	CALL	PRISTE(STNZ,ESN1,ESN2,2)
	IF(ESN1.LE.UNSTN.OR.ESN2.LE.UNSTN) NSTAT = 6
	ENDIF
	IF(NSTAT.EQ.6) GO TO 900
	IF(NSTAT.NE.1) GO TO 400

C	==================================================================
C	COMPUTE THE TOTAL ELASTIC STRAIN
C	==================================================================
	DO 5 I = 1,5
	ELSTN(I) = 0.0
5	ELSTN(I) = STNZ(I) - EPSTN(I)

C	==================================================================
C	COMPUTE THE DAMAGED POISSON RATIO
C	==================================================================	
	POISM = POISN
	IF(EMAX.GT.0.5) POISM = POISN*(1.8*(EMAX-0.5) + 1.0)
	IF(POISM.GE.0.4999) POISM = 0.4999

C	==================================================================
C	COMPUTE THE EQUIVALENT STRAIN
C	==================================================================	
	CALL	EFECTE(ELSTN,EED,EOPIM,POISM,SE)


C	==================================================================
C	COMPUTE THE EQUIVALENT PLASTIC STRAIN
C	==================================================================
	EP = EMAX - EEMAX

C	==================================================================
C	COMPUTE THE EQUIVALENT STRESS CORESPONDING TO LOADING STAGE
C	==================================================================
	IF(EED.GE.EEMAX) GO TO 200
	IF(EED.LT.EEMAX.AND.EED.GE.EEO) GO TO 150

C	==================================================================
C	UNLOADING STAGE
C	==================================================================
	EMAX = EEMAX + EP
	CALL	FRCPRE(EMAX,FKO)
	SLOP = FKO*FKO
	ADUM = (EED/EEO)**PN
	BDUM = SO/(FKO*EO*EEO) - SLOP
	ALPA = SLOP + BDUM*ADUM
	S    = FKO*EO*EED*ALPA
	EEO  = EED
	KUNLO = 3

	GO TO 300

C	==================================================================
C	RELOADING STAGE
C	==================================================================
150	CONTINUE

	S = SO + (SMAX - SO)*(EED - EEO)/(EEMAX - EEO)
	EEO  = EED
	KUNLO = 2

	GO TO 300

C	==================================================================
C	LOADING STAGE
C	==================================================================
200	CONTINUE

	EMAX  = EED + EP	
	CALL	FRCPRE(EMAX,FKO)
	ADUM  = 1.0 - EXP(-0.35*EMAX)
	EEO   = (20.0/7.0)*ADUM
	S     = FKO*EO*EEO
	KUNLO = 1

300	CONTINUE

C	==================================================================
C	UPDATE THE MAXIMUM AND CURRENT STRESS & STRAIN
C	==================================================================
	SO = S
	IF(SO.GE.SMAX) SMAX = SO
	IF(EEO.GE.EEMAX) EEMAX = EEO

C	==================================================================
C	COMPUTE THE FRACTURED YOUNG MODULUS AND STRESS-STRAIN MATRIX
C	==================================================================
	IF(SE.EQ.0.0) THEN
	YOUNM = YOUNG
	ELSE
	YOUNM = (UNIAX/EOPIM)*(S/EED)*(EED/SE)
	ENDIF
	
	CALL	MODUL1(DMATX,LPROP,NSTAT,  0.0,
	1			  0.0, 0.0,   YOUNM, 0.0,
     2			  POISM)

C	==================================================================
C	COMPUTE THE INCREMENTAL ELASTIC STRAIN
C	==================================================================
	DO 10 I = 1,5
	DELSN(I) = 0.0
10	DELSN(I) = ELSTN(I) - ELSNP(I)

C	==================================================================
C	COMPUTE THE INCREMENTAL PLASTIC STRAIN
C	==================================================================	
	CALL	FLOWSE(ELSTN,EPSTN,DPSTN,DELSN,EMAX,
	1		      DEMX2,LPROP,POISM)

C	==================================================================
C	COMPUTE THE TOTAL PLASTIC STRAIN
C	==================================================================	
	DO 20 I  = 1,5
20	EPSTN(I) = EPSTN(I) + DPSTN(I)

C	==================================================================
C	COMPUTE THE TOTAL ELASTIC STRAIN
C	==================================================================
	DO 30 I  = 1,5
	ELSTN(I) = 0.0
30	ELSTN(I) = STNZ(I)  - EPSTN(I)

C	==================================================================
C	COMPUTE THE TOTAL STRESSES
C	==================================================================	
	DO 350 I = 1,5
	STSPS(I) = 0.0
	DO 350 J = 1,5
350	STSPS(I) = STSPS(I) + DMATX(I,J)*ELSTN(J)

C	==================================================================
C	UPDATE THE VARIABLES
C	==================================================================
	DO 385 I = 1,5
385	ELSNP(I) = ELSTN(I)

	EQIST(1) = SMAX
	EQIST(2) = SO
	EQISN(1) = EMAX
	EQISN(2) = EEMAX
	EQISN(3) = EEO

	WOKPR(3) = FKO

	GO TO 800

C	==================================================================
C	START THE MODULE FOR CRACKED CONCRETE IN COMPRESSION
C	==================================================================
400	CONTINUE

C	==================================================================
C	COMPUTE STRAIN IN CRACK DIRECTION
C	==================================================================
	CALL	TRANSE(TRAMX,ANGL2)
	SNP = 0.0
	SNT = 0.0
	DO 450 J = 1,5
	SNT = SNT + TRAMX(2,J)*STNZ(J)	
450	SNP = SNP + TRAMX(1,J)*STNZ(J)	

C	==================================================================
C	CHECK FOR TENSION
C	==================================================================
	IF(SNP.GE.0.0) GO TO 710

C	==================================================================
C	CALL THE VARIABLES FROM PREVIOUS LOOP
C	==================================================================
	ECMAX = EQCRK(1)
	EOCUR = EQCRK(2)
	SIGMX = EQCRK(3)
	SIGCR = EQCRK(4)
	EEP   = EQCRK(5)
	ETMAX = EQCRK(6)

C	==================================================================
C	STORE THE MAXIMUM STRAIN NORMAL TO THE CRACK
C	==================================================================
	IF(SNT.GE.ETMAX) ETMAX = SNT

	CONFK = ABS(ECMAX)/EOPIM

	BETA  = 1.75

C	==================================================================
C	COMPUTE THE STRESS CORESPONDING TO LOADING STAGE
C	==================================================================
	IF(SNP.LE.ECMAX) GO TO 600                             !NEGATIVE SIGN
	IF(SNP.GT.ECMAX.AND.SNP.LE.EOCUR) GO TO 500            !NEGATIVE SIGN

C	==================================================================
C	UNLOADING STAGE
C	==================================================================
	CALL	FRCPRE(CONFK,FKO)
	
	SLOP = FKO*FKO
	ADUM = ((SNP-EEP)/(EOCUR-EEP))**PN
	BDUM = SIGCR/(FKO*YOUNG*(EOCUR-EEP)) - SLOP
	ALPA = SLOP + BDUM*ADUM

C	SIGG = FKO*YOUNG*(SNP-EPP)*ALPA    !RELEASE ERROR SONGSAK MAY2006 EPP --> EEP
	SIGG = FKO*YOUNG*(SNP-EEP)*ALPA

	GO TO 700

C	==================================================================
C	RELOADING STAGE
C	==================================================================
500	CONTINUE

	SIGG = SIGMX - (SIGMX - SIGCR)*(ECMAX - SNP)/(ECMAX - EOCUR)

	GO TO 700

C	==================================================================
C	LOADING STAGE
C	==================================================================
600	CONTINUE

	ECMAX = SNP	
	CONFK = ABS(SNP)/EOPIM	
	CALL	FRCPRE(CONFK,FKO)
	EVALU = 2.718281828
	ADUM  = 1.0 - EVALU**(-0.35*CONFK)
	EEP   = CONFK - (20.0/7.0)*ADUM
	EEP   = -EEP*BETA*EOPIM
	SIGG  = FKO*YOUNG*(SNP-EEP)	

700	CONTINUE

C	==================================================================
C	STORE THE MAXIMUM AND CURRENT STRESS & STRAIN
C	==================================================================
	EQCRK(1) = ECMAX
	EQCRK(2) = SNP
	IF(ABS(SIGG).GE.ABS(SIGMX)) EQCRK(3) = SIGG
	EQCRK(4) = SIGG
	EQCRK(5) = EEP
	EQCRK(6) = ETMAX

C	==================================================================
C	REDUCTION OF COMPRESSIVE STRESS DUE TO CRACK
C	==================================================================
	OMEGA = 1.0
	IF(ETMAX.GE.0.001) OMEGA = -100.0*ETMAX + 1.1
	IF(ETMAX.GE.0.005) OMEGA = 0.6
	SIGMA(2) = SIGG*OMEGA

710	CONTINUE
	
C	==================================================================
C	TRANSFER STRESSES TO GLOBAL DIRECTION
C	==================================================================
	CALL	TRANSE(TRAMX,ANGL1)
	DO 750 I = 1,5
	STSPS(I) = 0.0
	DO 750 J = 1,5
750	STSPS(I) = STSPS(I) + TRAMX(J,I)*SIGMA(J)


800	CONTINUE

	GO TO 1000

C	==================================================================
C	CRUSHING OF CONCRETE
C	==================================================================
900	CONTINUE
	DO 950 I = 1,5
950	STSPS(I) = 0.0
	
1000	CONTINUE


	RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE CCRACKE(STSPS,NSTAT,SIGMA,DIREC,
	1				  EPSTN,GRTSN,SGMAX,LPROP,
	2			      STNZ,WOKPR)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CRACKING BEHAVIOUR OF CONCRETE
C	=============================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
	DIMENSION STSPS(5),SIGMA(5),DMATX(5,5),STNZ(5),DSTNZ(5),
	1		  DSTSZ(5),DIREC(2),EPSTN(5),GRTSN(2),STRSP(5),
     2		  STRNP(5),SGTOT(5),TRAMX(5,5),WOKPR(3)

	
	YOUNG = PROPS(LPROP,1)
	TENST = PROPS(LPROP,5)
	UNIAX = PROPS(LPROP,6)
	TENSN = PROPS(LPROP,8)
	CONSE = PROPS(LPROP,9)
	ECRAK = TENST/YOUNG


	IF(NSTAT.EQ.2.OR.NSTAT.EQ.3) GO TO 30
	IF(DIREC(1).NE.0.0) GO TO 25

C	------------------------------------------------------------------
C	FOR THE FIRST CRACKING COMPUTE DIRECTION OF THE CRACK
C	------------------------------------------------------------------
	DUMMY = SGMAX - SIGMA(2)
	IF(DUMMY.EQ.0.0) DUMMY = 0.1E-20
	DIREC(1) = ATAN(SIGMA(3)/DUMMY)

C	------------------------------------------------------------------
C	COMPUTE TRANSFORMATION MATRIX
C	------------------------------------------------------------------
25	NSTAT = 2
30    ANGLE = DIREC(1)
	CALL	TRANSE(TRAMX,ANGLE)
	C = COS(ANGLE)
	S = SIN(ANGLE)

C	------------------------------------------------------------------
C	COMPUTE TWO PRINCIPAL STRESSES
C	------------------------------------------------------------------
	STRSP(1) = C*C*SIGMA(1) + S*S*SIGMA(2) + 2.0*C*S*SIGMA(3)
	STRSP(2) = SIGMA(1) + SIGMA(2) - STRSP(1)

C	------------------------------------------------------------------
C	TEST FOR OPENING OF CRACKS IN BOTH DIRECTIONS
C	------------------------------------------------------------------
	IF(NSTAT.EQ.2.AND.STRSP(2).GE.TENST) NSTAT = 3

C	------------------------------------------------------------------
C	TRANSFORM STRAINS INTO THE CRACK DIRECTION
C	------------------------------------------------------------------
	DO 40 I  = 1,5
	STRNP(I) = 0.0
	DO 40 J  = 1,5
40	STRNP(I) = STRNP(I) + TRAMX(I,J)*STNZ(J)

	STRA1 = STRNP(1)
	STRA2 = STRNP(2)
	THETA = 0.0
	CALL	MODUL1(DMATX,LPROP,NSTAT,THETA,
	1			  STRA1,STRA2, 0.0, 0.0,
     2			  0.0)

	DO 50 I  = 1,5
	SGTOT(I) = 0.0
	DO 50 J  = 1,5
50	SGTOT(I) = SGTOT(I) + DMATX(I,J)*STRNP(J)


C	------------------------------------------------------------------
C	MODIFY THE STRESSES CORRESPONDING TO TENSION STIFFENING	
C	------------------------------------------------------------------

	IF(STRA1.GE.TENSN) GO TO 52
	IF(ABS(STRA1).GE.GRTSN(1)) GO TO 51
	DUMM1    = TENST*(ECRAK/ABS(GRTSN(1)))**CONSE
	IF(DUMM1.LT.0.0) DUMM1 = 0.0
	SGTOT(1) = DUMM1*STRA1/GRTSN(1)
	GO TO 52
51	SGTOT(1) = TENST*(ECRAK/ABS(STRA1))**CONSE

52	IF(NSTAT.EQ.2) GO TO 55
	IF(STRA2.GE.TENSN) GO TO 55
	IF(ABS(STRA2).GE.GRTSN(2)) GO TO 54
	DUMM2    = TENST*(ECRAK/ABS(GRTSN(2)))**CONSE
	IF(DUMM2.LT.0.0) DUMM2 = 0.0
	SGTOT(2) = DUMM2*STRA2/GRTSN(2)
	GO TO 55
54	SGTOT(2) = TENST*(ECRAK/ABS(STRA2))**CONSE
55    CONTINUE

	FKO   = WOKPR(3)
	RF    = FKO*FKO*FKO
	CRIT  = TENST*RF
	IF(SGTOT(1).GT.CRIT) SGTOT(1) = CRIT
	IF(NSTAT.EQ.2) GO TO 60
	IF(SGTOT(2).GT.CRIT) SGTOT(2) = CRIT
60	CONTINUE

C	------------------------------------------------------------------
C	UPDATE THE MAXIMUM TENSILE STRAIN
C	------------------------------------------------------------------
	IF(STRA1.GT.GRTSN(1)) GRTSN(1) = STRA1
	IF(NSTAT.EQ.3.AND.STRA2.GT.GRTSN(2)) GRTSN(2) = STRA2

C	------------------------------------------------------------------
C	STORE EFFECTIVE STRESS & STRAIN
C	------------------------------------------------------------------
	WOKPR(1) = STRA1
	IF(NSTAT.EQ.3) WOKPR(2) = STRA2

	IF(NSTAT.EQ.2) GO TO 72

C	------------------------------------------------------------------
C	CHECK FOR CRACK CLOSE IN 2 DIRECTION
C	------------------------------------------------------------------
	IF(STRA2.GT.0.0) GO TO 72

C	------------------------------------------------------------------
C	IF THE CRACK IS CLOSE CHANGE THE MATERIAL NUMBER
C	------------------------------------------------------------------
	WOKPR(2) = 0.0
	NSTAT = 2

C	------------------------------------------------------------------
C	CHECK FOR CRACK CLOSE IN 1 DIRECTION
C	------------------------------------------------------------------
72	IF(STRA1.GT.0.0) GO TO 76
	IF(NSTAT.EQ.2) GO TO 74

C	------------------------------------------------------------------
C	IF THE CRACK IS CLOSE IN 1 DIRECTION 
C	BUT STILL OPEN IN 2 DIRECTION 
C	CHANGE THE DIRECTION OF CRACK AND MATERIAL NUMBER
C	------------------------------------------------------------------
	WOKPR(1) = WOKPR(2)
	WOKPR(2) = 0.0
	GAMMA = 1.570796327
	IF(DIREC(1).GE.0.0) GAMMA = -GAMMA
	DIREC(1) = DIREC(1) + GAMMA
	NSTAT = 2
	GO TO 76

C	------------------------------------------------------------------
C	CRACK ARE CLOSED IN ALL DIRECTION
C	------------------------------------------------------------------
74	WOKPR(1) = 0.0
	NSTAT = 1


76	CONTINUE

C	==================================================================
C	IF TWO DIRECTIONS ARE CRACKED --> TRANSFER STRESSES TO GLOBAL SYSTEM
C	==================================================================
	DO 70 I = 1,5
	SIGMA(I) = 0.0
	DO 70 J = 1,5
	IF(NSTAT.EQ.3) THEN
	SIGMA(I) = SIGMA(I) + TRAMX(J,I)*SGTOT(J)
	ENDIF
70	CONTINUE

C	==================================================================
C	IF TWO DIRECTIONS ARE CRACKED --> COMPUTE THE REAL STRESSES
C	IF ONE DIRECTIONS ARE CRACKED --> STORE STRESSES INTO TEMPORARY ARRAY
C	==================================================================
	DO 80 I = 1,5
	IF(NSTAT.EQ.3) THEN
	STSPS(I) = SIGMA(I)
	ELSE
	SIGMA(I) = SGTOT(I)
	ENDIF
80	CONTINUE

C	==================================================================
C	COMPUTE THE ANGLE NORMAL TO THE CRACK
C	==================================================================	
	GAMMA = 1.570796327
	IF(DIREC(1).GE.0.0) GAMMA = -GAMMA
	DIREC(2) = DIREC(1) + GAMMA


	RETURN
	END
C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE	STELRE(STSPS,NSTAT,EPSTN,DSTSZ,LPROP,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	ELASTIC & PLASTIC BEHAVIOR OF REINFORCEMENT STEEL
C	=============================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
	DIMENSION	STSPS(5),EPSTN(5),DSTSZ(5)


	ANGLE = PROPS(LPROP,7)+REFANG
	YOUNG = PROPS(LPROP,1)
	HARDS = PROPS(LPROP,2)
	YOUN2 = YOUNG*HARDS/(YOUNG+HARDS)
	YIELD = PROPS(LPROP,5)
	C  = COS(ANGLE)
	S  = SIN(ANGLE)
	C1 = C*C
	C2 = S*S
	C3 = 2*C*S
	STREP = C1*DSTSZ(1) + C2*DSTSZ(2) + C3*DSTSZ(3)
	STRPP = C1*STSPS(1) + C2*STSPS(2) + C3*STSPS(3)
	STOTP = STRPP + STREP

	PREYS = YIELD + HARDS*EPSTN(1)
	IF(ABS(STRPP).GE.PREYS) GO TO 20
	ESCUR = ABS(STOTP) - PREYS
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(STREP)
	GO TO 30

20	IF(STOTP.GT.0.0.AND.STREP.LT.0.0) GO TO 35
	IF(STOTP.LT.0.0.AND.STREP.GT.0.0) GO TO 35
	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	STRPP = STRPP + REDUC*STREP + RFACT*STREP*(1-YOUNG/(YOUNG+HARDS))
	EPSTN(1) = EPSTN(1) + RFACT*STREP/(YOUNG+HARDS)
	NSTAT = 5
	EPSTN(2) = 0.0
	GO TO 50

35	NSTAT = 4
40	STRPP = STRPP + STREP
50	STSPS(1) = C1*STRPP
	STSPS(2) = C2*STRPP
	STSPS(3) = 0.5*C3*STRPP
	STSPS(4) = 0.0
	STSPS(5) = 0.0


	RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE	CONSTVE(STSPS,GRTSN,NSTAT,EPSTN,
	1				   DIREC,LPROP,DMATX,KUNLO,
     2				   WOKPR,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CONSTITUTIVE LAW FOR EACH LAYER
C	=============================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
	DIMENSION	STSPS(5),GRTSN(2),DIREC(2),EPSTN(5),
	1			DMATX(5,5),WOKPR(3)


	ANGLE = DIREC(1)
	IDTCS = PROPS(LPROP,10)

	IF(IDTCS.EQ.1) ANGLE = PROPS(LPROP,7)+REFANG
	IF(IDTCS.EQ.1) GO TO 100

C	------------------------------------------------------------------
C	CONTRIBUTION OF CONCRETE RIGIDITY
C	------------------------------------------------------------------

	CONSE = PROPS(LPROP,9)
	TENSN = PROPS(LPROP,8)
	TENST = PROPS(LPROP,5)
	YOUNG = PROPS(LPROP,1)
	POISN = PROPS(LPROP,2)
	ECRAK = TENST/YOUNG
	STRA1 = WOKPR(1)
	STRA2 = WOKPR(2)

C	------------------------------------------------------------------
C	CHECK FOR YIELDING OF CONCRETE	
C	------------------------------------------------------------------
	IF(NSTAT.NE.1) GO TO 80

C	------------------------------------------------------------------
C	CHECK FOR UNLOADING
C	------------------------------------------------------------------
	IF(KUNLO.NE.1) GO TO 80

C	------------------------------------------------------------------
C	LOADING OF CONCRETE (COMPRESSION)
C	------------------------------------------------------------------
	CALL	MODUL1(DMATX,LPROP,NSTAT,    0.0,
	1			  0.0,	0.0,  YOUNG,	0.0,
     2			  POISN)

	
	GO TO 90

80	CONTINUE

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
	YOUNA = 0.0
	YOUNB = 0.0
	IF(NSTAT.EQ.1) THEN
	YOUNA = YOUNG
	GO TO 83
	ENDIF

C	YOUNA = 0.0
C	YOUNB = 0.0
	IF(NSTAT.EQ.2) GO TO 81
	IF(STRA2.EQ.GRTSN(2)) GO TO 81
	IF(GRTSN(2).GT.TENSN) GO TO 81
	IF(GRTSN(2).EQ.0.0)   GO TO 81
	YOUNB = (TENST*(ECRAK/STRA2)**CONSE)/GRTSN(2)

81	IF(STRA1.EQ.GRTSN(1)) GO TO 83
	IF(GRTSN(1).GT.TENSN) GO TO 83
	IF(GRTSN(1).EQ.0.0)   GO TO 83
	YOUNA = (TENST*(ECRAK/STRA1)**CONSE)/GRTSN(1)

	

83	CALL	MODUL1(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,YOUNA,YOUNB,
     2			  POISN)


	GO TO 90

C	------------------------------------------------------------------
C	CONTRIBUTION OF STEEL RIGIDITY
C	------------------------------------------------------------------
100	CONTINUE

	CALL	MODUL1(DMATX,LPROP,NSTAT,ANGLE,
	1			  0.0,	0.0,	0.0,  0.0,
     2			  0.0)

90	CONTINUE


	RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE EFECTE(STN,EFT,EOP,POISM,SE)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	EVALUATE THE EQUIVALENT STRAIN
C	==================================================================
	DIMENSION	STN(5)


	CALL	PRISTE(STN,E1,E2,2)
	

	EM  = SQRT(2.0)*(E1+E2)/2.0
	YD  = SQRT(2.0)*(E1-E2)/2.0
	
	EFT = SQRT( (0.62*EM/EOP)**2.0 + (0.98*YD/EOP)**2.0 )
		
	SE  = SQRT( (0.60*EM/EOP/(1.0-POISM))**2.0 + 
	1		    (1.30*YD/EOP/(1.0+POISM))**2.0 )	
  

	RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE FRCPRE(EMAX,FKO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	COMPUTE FRACTURE PARAMETER
C	==================================================================

	AVALU = 1.0 - EXP(-1.25*EMAX)
	BVALU = -0.73*EMAX
	FKO   = EXP(BVALU*AVALU)

	RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE FLOWSE(ELSTN,EPSTN,DPSTN,DELSN,EMAX,
	1		         DEMX2,LPROP,POISM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	DETERMINE FLOW MATRIX
C	=============================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
	DIMENSION ELSTN(5),EPSTN(5),ELSND(5),DELTA(5),
	1		  ELSTD(5),DPSTN(5),DELSN(5),DDLSN(5)

	YOUNG = PROPS(LPROP,1)
	POISN = PROPS(LPROP,2)
	UNIAX = PROPS(LPROP,6) 

	EO    = 1.6*(1.0+POISN)*UNIAX/YOUNG

	CALL DEVIAE(ELSTN,ELSND,SE1,DE2,2,POISM)
	CALL DEVIAE(DELSN,DDLSN,DSE1,DDE2,2,POISM)

	EZZ  = (-ELSND(1)-ELSND(2))*POISM
	DEZZ = (-DDLSN(1)-DDLSN(2))*POISM
	DELJ2 = 0.0

	DO 5 I = 1,5
5	DELJ2 = DELJ2 + ELSTN(I)*DELSN(I)
	DELJ2 = DELJ2 + EZZ*DEZZ
	IF(DE2.NE.0.0) DELJ2 = 0.5*DELJ2/DE2


	IF(DELJ2.GE.0.0.AND.DE2.GE.DEMX2) THEN
	UP = 1.0
	ELSE
	UP = 0.0
	ENDIF

	IF(DE2.GT.DEMX2) DEMX2 = DE2

	CONS1 = DE2*DE2/EO/EO
	CONS1 = (27.0/10.0)*CONS1

	DCON0 = 1.0 - 2.0*POISN
	DCON0 = -DCON0/SQRT(3.0)/(1.0+POISN)

	DCON1 = SQRT(2.0)*SE1 + 0.38*EO
	DCON1 = DCON1 / 0.28 / EO

	XVALU = SQRT(3.0)*SE1/DE2

	IF(XVALU.GT.-1.0.AND.XVALU.LT.1.0) THEN
	PVALU = -0.5*SIN(1.570796*XVALU) + 0.5
	ENDIF
	IF(XVALU.LE.-1.0) PVALU = 1.0
	IF(XVALU.GE.1.0)  PVALU = 0.0

	CALL	FRCPRE(EMAX,FKO)

	DVALU = DCON0*(2.0*FKO)**2.0 +
	1		DCON1*(1.0-4.0*FKO*FKO)
	IF(FKO.GE.0.5) DVALU = DCON0

	DVALU = DVALU*PVALU
	
	DELTA(1) = 1.0
	DELTA(2) = 1.0
	DELTA(3) = 0.0
	DELTA(4) = 0.0
	DELTA(5) = 0.0

	DO 10 I = 1,5
	CUM1 = 0.0
	IF(DE2.NE.0.0) CUM1 = ELSND(I)/DE2
	CUM2 = DELTA(I)*DVALU
	CUM3 = CONS1*UP
	DPSTN(I) = 0.0
	DPSTN(I) = (CUM1+CUM2)*CUM3*DELJ2
10	CONTINUE
	

	RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE DEVIAE(STN,STND,ST1,DT2,LGF,POISM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	COMPUTE THE STRESS INVARIANT
C	=============================
C	LGF = FLAG(1=STRESS,2=STRAIN)
C	==================================================================
	DIMENSION STN(5),STND(5)
	
	IF(LGF.EQ.1) THEN
	ST1 = STN(1) + STN(2) + 0.0
	ELSE
	STNZZ = -STN(1)*POISM - STN(2)*POISM
	ST1   = STN(1) + STN(2) + STNZZ
	ENDIF

	ST1 = ST1/3.0

	STND(1) = STN(1) - ST1
	STND(2) = STN(2) - ST1
	STND(3) = STN(3) 
	STND(4) = STN(4)
	STND(5) = STN(5)

	IF(LGF.EQ.1) THEN
	DT2 = 0.5*( STND(1)*STND(1) + STND(2)*STND(2) + 
	1	  STND(3)*STND(3) + STND(4)*STND(4) + STND(5)*STND(5) )
	ELSE
	STNZZ = STNZZ - ST1
	DT2 = 0.5*( STND(1)*STND(1) + STND(2)*STND(2) + STNZZ*STNZZ +
	1	0.25*(STND(3)*STND(3) + STND(4)*STND(4) + STND(5)*STND(5)))
	ENDIF

	DT2 = SQRT(DT2) 


	RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE	MODUL1(DMATX,LPROP,NSTAT,ANGLE,
	1				  STRA1,STRA2,YOUNA,YOUNB,
     2				  POISM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CONSTRUCT THE CONSTITUTIVE MATRIX ACCORDING TO MATERIAL NUMBER
C	=============================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
	DIMENSION DMATX(5,5),GASHM(5,5),TRAMX(5,5)


	DO 10 I = 1,5
	DO 10 J = 1,5
10	DMATX(I,J) = 0.0
	
	YOUNG = PROPS(LPROP,1)
	GO TO (1,2,2,4,5,6),NSTAT

C	------------------------------------------------------------------
C	ELASTIC & PLASTIC CONCRETE BEHAVIOR
C	------------------------------------------------------------------
1	COEF1 = 5.0/6.0
	CONS1 = YOUNA/(1.0-POISM*POISM)
	CONS2 = YOUNA/(2.0*(1.0+POISM))
	DMATX(1,1) = CONS1
	DMATX(2,2) = CONS1
	DMATX(1,2) = POISM*CONS1
	DMATX(2,1) = POISM*CONS1
	DMATX(3,3) = CONS2
	DMATX(4,4) = COEF1*CONS2
	DMATX(5,5) = COEF1*CONS2
	GO TO 100

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
2	POISN = PROPS(LPROP,2)
	GMODU = YOUNG/(2.0*(1.0+POISN))
	CONSG = 0.25
	GTENS = 0.004
	YOUN2 = YOUNG
	GMOD2 = GMODU*5.0/6.0
	YOUN1 = YOUNA
	GMOD1 = CONSG*GMODU*(1.0-ABS(STRA1)/GTENS)
	IF(GMOD1.LE.0.0) GMOD1=0.0
	IF(NSTAT.EQ.2) GO TO 35
	YOUN2 = YOUNB
	GMOD2 = CONSG*GMODU*(1.0-ABS(STRA2)/GTENS)
	IF(GMOD2.LE.0.0) GMOD1=0.0
	GMOD3 = GMOD1
	IF(GMOD2.LE.GMOD1) GMOD3 = GMOD2
35    DMATX(1,1) = YOUN1
	DMATX(2,2) = YOUN2
	DMATX(3,3) = GMOD1
	IF(NSTAT.EQ.3) DMATX(3,3) = 0.5*GMOD3
	DMATX(4,4) = GMOD1
	DMATX(5,5) = GMOD2
	GO TO 50

C	------------------------------------------------------------------
C	ELASTIC STEEL BEHAVIOR
C	------------------------------------------------------------------
4	DMATX(1,1) = YOUNG
	GO TO 50
C	------------------------------------------------------------------
C	ELASTIC-PLASTIC STEEL BEHAVIOR
C	------------------------------------------------------------------
5     DMATX(1,1) = PROPS(LPROP,2)

C	------------------------------------------------------------------
C	FORM THE TRANSFORMATION MATRIX
C	------------------------------------------------------------------
50    IF(ABS(ANGLE).LT.0.001) GO TO 100
	CALL TRANSE(TRAMX,ANGLE)

	DO 60 I = 1,5
	DO 60 J = 1,5
	GASHM(I,J) = 0.0
	DO 60 K = 1,5
60	GASHM(I,J) = GASHM(I,J) + DMATX(I,K)*TRAMX(K,J)

	DO 70 I = 1,5
	DO 70 J = 1,5
	DMATX(I,J) = 0.0
	DO 70 K = 1,5
70	DMATX(I,J) = DMATX(I,J) + GASHM(K,J)*TRAMX(K,I)


C	------------------------------------------------------------------
C	CRUSHING OF CONCRETE
C	------------------------------------------------------------------
6	CONTINUE

100   RETURN
	END


C	==================================================================
C	==================================================================
	SUBROUTINE PRISTE(STN,SGMAX,SGMIN,LGF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	DETERMINE THE PRINCIPAL STRESS
C	==============================
C	LGF = FLAG(1=STRESS,2=STRAIN)
C	==================================================================
	DIMENSION ST(5),STN(5)

	ST(1) = STN(1)
	ST(2) = STN(2)
	ST(3) = STN(3)
	ST(4) = STN(4)
	ST(5) = STN(5)
	IF(LGF.EQ.2) THEN
	ST(3) = 0.5*ST(3)
	ST(4) = 0.5*ST(4)
	ST(5) = 0.5*ST(5)
	ENDIF

	GASH1 = (ST(1) + ST(2))*0.5
	GASH2 = (ST(1) - ST(2))*0.5
	GASH3 = SQRT(GASH2*GASH2 + ST(3)*ST(3))
	SGMAX = GASH1 + GASH3
	SGMIN = ST(1) + ST(2) - SGMAX


	RETURN
	END


C	==================================================================
	SUBROUTINE TRANSE(TRAMX,ANGLE)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CONSTRUCT THE TRANSFORMATION MATRIX
C	==================================================================
	DIMENSION TRAMX(5,5)

	DO 10 I = 1,5
	DO 10 J = 1,5
10    TRAMX(I,J) = 0.0
	C = COS(ANGLE)
	S = SIN(ANGLE)
	TRAMX(1,1) = C*C
	TRAMX(2,2) = C*C
	TRAMX(1,2) = S*S
	TRAMX(2,1) = S*S
	TRAMX(1,3) = C*S
	TRAMX(3,1) = -2.0*C*S
	TRAMX(2,3) = -1.0*C*S
	TRAMX(3,2) = 2.0*C*S
	TRAMX(3,3) = C*C - S*S
	TRAMX(4,4) = C
	TRAMX(4,5) = S
	TRAMX(5,4) = -S
	TRAMX(5,5) = C

	
	RETURN 
	END

C	==================================================================
	SUBROUTINE INVA2E(SN,CRUSH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	EVALUATE THE CRUSHING FUNCTION FOR THE GIVEN STRAIN
C	==================================================================
	DIMENSION SN(5)

	BETA  = 1.355
	ALPHA = 0.355
	C1 = BETA*( SN(1)*SN(1) + SN(2)*SN(2) - SN(1)*SN(2)+
	1	 0.75*(SN(3)*SN(3) + SN(4)*SN(4) + SN(5)*SN(5)) )
	C2 = ALPHA*(SN(1) + SN(2))*0.5
	C3 = SQRT(C2*C2 + C1)
	CRUSH = C2 + C3

	RETURN
	END


C	==================================================================
C	========================= END OF EPF MODEL =======================
C	==================================================================





C	==================================================================
C	==================================================================
C	===== PLASTICITY MODEL OF REINFORCED CONCRETE SHELL ELEMENT ======
C	==================================================================
	SUBROUTINE RCLAYR(WA,STNCR,STSCR,IELE,TH,DP,DVOL,MSET,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	PRODUCED BY SONGSAK - DECEMBER 2004
C	==================================================================
C	COMPUTE THE GAUSS POINT STRESSES AND MATERIAL RIGIDITY FOR
C	REINFORCED CONCRETE SHELL ELEMENT USING LAYERED MODEL
C	==================================================================
C	INPUT VARIABLE
C	==============
C	STNCR(8)     = CURRENT GAUSS POINT STRAIN (MEMBRANE-BENDING-SHEAR)
C	TH			 = GAUSS POINT THICKNESS
C	DVOL		 = GAUSS POINT INTEGRATION VOLUME
C
C	==================================================================
C	INITAILIZE VARIABLE
C	===================
C	STNPS(5) = GAUSS POINT STRAIN AT THE END OF PREVIOUS UPDATE
C	STSPS(5) = GAUSS POINT STRESS AT THE END OF PREVIOUS UPDATE
C	DIREC(2) = DIRECTION OF CRACK (1 & 2-DIRECTION)
C			   (MEASURE ANTICLOCKWISE IN RADIANS)
C	EFFST(1) = GAUSS POINT EFFECTIVE STRESS
C	EPSTN(2) = GAUSS POINT EFFECTIVE STRAIN (1 & 2-DIRECTION)
C	GRTST(2) = GAUSS POINT GREATEST STRAIN (1 & 2-DIRECTION)
C			   (USED FOR TENSION STIFFENING MODEL)
C	MSTAT    = MATERIAL NUMBER
C				 1 = ELASTIC BEHAVIOR CONCRETE
C				 2 = CRACKED CONCRETE (CRACK IN ONE DIRECTION)
C				 3 = CRACKED CONCRETE (CRACK IN TWO DIRECTION)
C				 4 = ELASTO-PLASTIC BEHAVIOR CONCRETE
C				 5 = ELASTIC BEHAVIOR OF STEEL
C				 6 = ELASTO-PLASTIC BEHAVIOR STEEL
C				 7 = CONCRETE CRUSHING
C
C	==================================================================
C	OUTPUT VARIABLE
C	===============
C	STSCR(8)     = CURRENT GAUSS POINT STRESS (MEMBRANE-BENDING-SHEAR)
C	DP(8,8)		 = MATERIAL RIGIDITY MATRIX
C
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C	
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION WA(20,1),STSPS(5),STNPS(5),STNCR(8),STSCR(8),
	1		  DMATX(5,5),STNZ(5),DSTNZ(5),GRTST(2),EPSTN(2),
	2		  DSTSZ(5),SIGMA(5),DIREC(2),EFFST(1),DP(8,8),
	3		  CP(5,5),HDMAT(100,2)



	VOLM  = DVOL
	LLAYR = MSET
	NLAYR = MATLR(MSET,MNLYR+1) !MNLYR
	ZETAD  = -1.0

	DO 250 I = 1,8
	STSCR(I) = 0.0
	DO 250 J = 1,8
250	DP(I,J)  = 0.0



      
C	------------------------------------------------------------------
C	LOOP OVER CONCRETE & STEEL LAYER
C	------------------------------------------------------------------
      DO 200 ILAYR = 1,NLAYR
	LPROP = MATLR(LLAYR,ILAYR)
      IDTCS = PROPS(LPROP,10)

	DO 8 I = 1,5
	J = I+5
	STNPS(I) = WA(I,ILAYR)
8	STSPS(I) = WA(J,ILAYR)
	DIREC(1) = WA(11,ILAYR)
	DIREC(2) = WA(12,ILAYR)
	EFFST(1) = WA(13,ILAYR)
	EPSTN(1) = WA(14,ILAYR)
	EPSTN(2) = WA(15,ILAYR)
	GRTST(1) = WA(16,ILAYR)
	GRTST(2) = WA(17,ILAYR)
	MSTAT    = WA(18,ILAYR)

	FLAG_ANG = WA(19,ILAYR)
	IF(FLAG_ANG.NE.0.0) REFANG = WA(20,ILAYR)
		
	IF(IDTCS.EQ.0) THEN
	DZETA = PROPS(LPROP,3)
	ZETA  = ZETAD + DZETA/2.0
	ZETAG = 0.5*TH*ZETA
	DZETG = 0.5*TH*DZETA
	ELSE
	DZETA = PROPS(LPROP,3)
	ZETA  = PROPS(LPROP,6)
	ZETAG = 0.5*TH*ZETA
	DZETG = 0.5*TH*DZETA
	ENDIF

C	------------------------------------------------------------------
C	COMPUTE STRAINS AT THE MIDDLE OF LAYER
C	------------------------------------------------------------------
	DO 10 I = 1,3
	STNZ(I)  = STNCR(I) + ZETAG*STNCR(I+3)
10	DSTNZ(I) = STNZ(I)  - STNPS(I)
	STNZ(4)  = STNCR(7)
	STNZ(5)  = STNCR(8)
	DSTNZ(4) = STNCR(7) - STNPS(4)
	DSTNZ(5) = STNCR(8) - STNPS(5)
	
	IF(IDTCS.EQ.1) DIREC(1) = PROPS(LPROP,7)
	IF(IDTCS.EQ.1) GO TO 100

      IF(MOPTN.EQ.2)	CALL HDATM(HDMAT,LPROP)

C	------------------------------------------------------------------
C	COMPUTE STRESSES AT THE MIDDLE OF CONCRETE LAYER
C	------------------------------------------------------------------
	NSTAT = MSTAT
	ANGLE = DIREC(1)
	STRA1 = EPSTN(1)
	STRA2 = EPSTN(2)
	CALL	MODUL(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,	0.0,	0.0)
	DO 15 I = 1,5
	DSTSZ(I) = 0.0
	DO 15 J = 1,5
15	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)

	DO 20 I = 1,5
20	SIGMA(I) = STSPS(I) + DSTSZ(I)

C	------------------------------------------------------------------
C	COMPUTE MAX & MIN STRESS IN THE PRINCIPAL DIRECTION
C	------------------------------------------------------------------
	CALL	PRIST(SIGMA,SGMAX,SGMIN)

	TENST = PROPS(LPROP,5)
	UNIAX = PROPS(LPROP,6)

C	------------------------------------------------------------------
C	CHECK FOR CRAKING OF CONCRETE
C	------------------------------------------------------------------
	IF(SGMAX.GT.TENST) GO TO 30
	IF(NSTAT.EQ.1.OR.NSTAT.EQ.4) GO TO 40
	IF(NSTAT.EQ.7) GO TO 40

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
30	CALL	CCRACK(STSPS,STNPS,STNCR,EFFST,
	1			   SIGMA,NSTAT,DIREC,EPSTN,
	2			   GRTST,SGMAX,LPROP,YIELD,
	3			   STNZ,HDMAT)

C	------------------------------------------------------------------
C	CHECK FOR YIELDING OF CONCRETE	
C	------------------------------------------------------------------
	IF(MOPTN.EQ.2) THEN
	CALL	HARDN(EPSTN,LPROP,HARDS,YIELS,HDMAT)
	ELSE
	YIELS = UNIAX
	HARDS = 0.0
	ENDIF

	IF(YIELD.LE.YIELS) GO TO 50
	DO 35 I = 1,5
35	DSTSZ(I) = SIGMA(I) - STSPS(I)

C	------------------------------------------------------------------
C	COMPRESION BEHAVIOR OF CONCRETE
C	------------------------------------------------------------------
40	CALL	CYIELD(STSPS,STNPS,STNCR,EFFST,
	1			   NSTAT,DIREC,EPSTN,STNZ,
	2			   DSTSZ,SIGMA,LPROP,HDMAT)
	
50	CONTINUE

	GO TO 95

C	------------------------------------------------------------------
C	COMPUTE STRESS AT THE MIDDLE OF STEEL LAYER
C	------------------------------------------------------------------
100	NSTAT = MSTAT
	ANGLE = PROPS(LPROP,7)+REFANG

	CALL	MODUL(DMATX,LPROP,	5,ANGLE,
	1			  0.0,	0.0,	0.0,	0.0)
	DO 60 I = 1,5
	DSTSZ(I) = 0.0
	DO 60 J = 1,5
60	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)
	CALL	STELR(STSPS,NSTAT,DIREC,EFFST,EPSTN,
	1			  DSTSZ,LPROP,REFANG)

95	CONTINUE

C	UPDATE THE MATERIAL NUMBER
	MSTAT = NSTAT

C	UPDATE THE LAYER STRAIN
	DO 96 I = 1,5 
96	STNPS(I) = STNZ(I)
	

C	STORE THE CURRENT STRESSES IN WA
	DO 97 I = 1,5
	J = I+5
	WA(I,ILAYR)  = STNPS(I)
97	WA(J,ILAYR)  = STSPS(I)
	WA(11,ILAYR) = DIREC(1)
	WA(12,ILAYR) = DIREC(2)
	WA(13,ILAYR) = EFFST(1)
	WA(14,ILAYR) = EPSTN(1)
	WA(15,ILAYR) = EPSTN(2)
	WA(16,ILAYR) = GRTST(1)
	WA(17,ILAYR) = GRTST(2)
	WA(18,ILAYR) = MSTAT

	IF(FLAG_ANG.EQ.0.0) THEN
	WA(19,ILAYR) = 1.0
	WA(20,ILAYR) = REFANG 
	ENDIF
	
C	CONTRIBUTION OF STRESSES OVER THICKNESS
C	MEMBRANE - BENDING - SHEAR
	DO 150 I = 1,3
	J = I+3
	STSCR(I) = STSCR(I) + STSPS(I)*DZETG
150	STSCR(J) = STSCR(J) + STSPS(I)*DZETG*ZETAG
	STSCR(7) = STSCR(7) + STSPS(4)*DZETG
	STSCR(8) = STSCR(8) + STSPS(5)*DZETG

	CALL	CONSTV(STSPS,GRTST,NSTAT,EPSTN,DIREC,
	1			   LPROP,CP,HDMAT,REFANG)

C	INITIALIZE ELASTO-PLASTIC RIGIDITY MATRIX	
	FA = DZETG*VOLM
	FC = FA*ZETAG
	FB = FC*ZETAG
	FE = DZETG*DZETG*DZETG/12.0

C	MEMBRANE & BENDING RIGIDITY
	DO 300 I = 1,3
	DO 300 J = 1,3
	K = I+3
	L = J+3
	DP(I,J) = DP(I,J) + CP(I,J)*FA
300	DP(K,L) = DP(K,L) + CP(I,J)*(FB+FE)

C	SHEAR RIGIDITY
	DO 350 I = 4,5
	DO 350 J = 4,5
	K = I+3
	L = J+3
350	DP(K,L) = DP(K,L) + CP(I,J)*FA

C	COUPLE MEMBRANE AND BENDING RIGIDITY
C	(UPPER PART OF TRIANGLE)
	DO 380 I = 1,3
	DO 380 J = 1,3
	L = J+3
380	DP(I,L) = DP(I,L) + CP(I,J)*FC

C	FILL IN LOWER PART OF RIGIDITY MATRIX
	DO 400 I = 1,6
	DO 400 J = 1,6
400	DP(J,I) = DP(I,J)


	IF(IDTCS.EQ.0) THEN
	ZETAD  = ZETA + DZETA/2.0
	ENDIF
      
200	CONTINUE


	RETURN
	END


C	==================================================================
	SUBROUTINE CYIELD(STSPS,STNPS,STNCR,EFFST,
	1				  NSTAT,DIREC,EPSTN,STNZ,
	2				  DSTSZ,SIGMA,LPROP,HDMAT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	COMPRESSION BEHAVIOR OF CONCRETE INCLUDED 
C	ELASTIC,ELASTO-PLASTIC AND CRUSHING 
C	==================================================================
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION STSPS(5),STNPS(5),STNCR(8),STNZ(5),EFFST(1),
	1		  DSTSZ(5),SIGMA(5),DIREC(2),EPSTN(2),AVECT(5),
	2		  DVECT(5),HDMAT(100,2)


	IF(NSTAT.NE.2) GO TO 10
C	------------------------------------------------------------------
C	BACK UP DATA FOR CRACKED CONCRETE
C	------------------------------------------------------------------
	DUMM1    = EPSTN(1)
	EPSTN(1) = EPSTN(2)
	DUMM2    = EFFST(1)
	EFFST(1) = DIREC(2)

10    CONTINUE
	IF(NSTAT.EQ.7) GO TO 160
	UNIAX = PROPS(LPROP,6)
	UNSTN = PROPS(LPROP,7)
	YOUNG = PROPS(LPROP,1)

C	------------------------------------------------------------------
C	TEST FOR CRUSHING
C	------------------------------------------------------------------
	EPSTN(1) = ABS(EPSTN(1))
	CALL	INVA2(STNZ,CRUSH)
	IF(NSTAT.EQ.2) CRUSH = EFFST(1)/YOUNG + EPSTN(1)
	IF(CRUSH.GE.UNSTN) GO TO 160

	IF(MOPTN.EQ.2) THEN
	CALL	HARDN(EPSTN,LPROP,HARDS,YIELS,HDMAT)
	ELSE
	YIELS = UNIAX
	HARDS = 0.0
	ENDIF
C	------------------------------------------------------------------
C	REDUCE STRESSES TO THE YIELD SURFACE
C	------------------------------------------------------------------
	PREYS = YIELS
	CALL	INVAR(SIGMA,YIELD)
	ESPRE = EFFST(1) - PREYS


C	CHECK FOR PREVIOUS YIELDING OF CONCRETE
	IF(ESPRE.GE.0.0) GO TO 55
	ESCUR = YIELD - PREYS
C	CHECK FOR EFFECTIVE STRESS SATISFY YIELD CRITERIA
	IF(ESCUR.LE.0.0) GO TO 60
	CALL	CONFAC(STSPS,DSTSZ,PREYS,RFACT)
	GO TO 70
55	ESCUR = YIELD - EFFST(1)
C	CHECK FOR UNLOADED
	IF(ESCUR.LE.0.0) GO TO 60
	RFACT = 1.0

C	COMPUTE NUMBER OF INCREMENTAL
70	MSTEP = ESCUR*8.0/UNIAX + 1.0
	IF(MSTEP.LT.30) MSTEP = 30
	ASTEP = MSTEP
	REDUC = 1.0 - RFACT

	DO 80 I = 1,5
	STSPS(I) = STSPS(I) + REDUC*DSTSZ(I)
80    DSTSZ(I) = RFACT*DSTSZ(I)/ASTEP
	
C	------------------------------------------------------------------
C	LOOP OVER EACH INCREMENTAL
C	------------------------------------------------------------------
	DO 90 ISTEP = 1,MSTEP
	CALL	INVAR(STSPS,YIELD)
	CALL    FLOWS(ABETA,AVECT,DVECT,LPROP,STSPS,EPSTN,HDMAT)
	ADUMM = 0.0
	DO 100 I = 1,5
100	ADUMM = ADUMM + AVECT(I)*DSTSZ(I)
	DLAMD = ADUMM*ABETA
	IF(DLAMD.LT.0.0) DLAMD = 0.0
	BDUMM = 0.0
	DO 110 I = 1,5
	BDUMM = BDUMM + AVECT(I)*STSPS(I)
110	STSPS(I) = STSPS(I) + DSTSZ(I)
	1					 - DLAMD*DVECT(I)
	EPSTN(1) = EPSTN(1) + DLAMD*BDUMM/YIELD
90	CONTINUE
	
	CALL	INVAR(STSPS,YIELD)

	IF(MOPTN.EQ.2) THEN
	CALL	HARDN(EPSTN,LPROP,HARDS,YIELS,HDMAT)
	ELSE
	YIELS = UNIAX
	HARDS = 0.0
	ENDIF

	CURYS = YIELS
	BRING = CURYS/YIELD
	DO 130 I = 1,5
130	STSPS(I) = BRING*STSPS(I)
	MSTAT = 4

	GO TO 190

160	MSTAT = 7
	DO 165 I = 1,5
165	STSPS(I) = 0.0
	EFFST(1) = 0.0

	GO TO 190

C	------------------------------------------------------------------
C	ELASTIC BEHAVIOR (WITH UNLOADING)
C	------------------------------------------------------------------
60	MSTAT = 1
	DO 180 I = 1,5
180	STSPS(I) = STSPS(I) + DSTSZ(I)
	EFFST(1) = YIELD
	IF(EPSTN(1).EQ.0.0.OR.ESCUR.EQ.0.0) GO TO 190
	EPSTN(1) = -EPSTN(1)

190   CONTINUE

	IF(NSTAT.NE.2) GO TO 200

C	------------------------------------------------------------------
C	RETURN BACK UP DATA FOR CRACKED CONCRETE
C	------------------------------------------------------------------
	EPSTN(2) = EPSTN(1)
	EPSTN(1) = DUMM1
	DIREC(2) = EFFST(1)
	EFFST(1) = DUMM2

	GO TO 220

200	IF(MSTAT.EQ.1) THEN
	NSTAT = NSTAT
	ELSE
	NSTAT = MSTAT
	ENDIF
		    	
220	CONTINUE

	RETURN
	END


C	=================================================================
	SUBROUTINE CCRACK(STSPS,STNPS,STNCR,EFFST,
	1				  SIGMA,NSTAT,DIREC,EPSTN,
	2				  GRTST,SGMAX,LPROP,YIELD,
	3				  STNZ,HDMAT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	TENSION BEHAVIOR OF CONCRETE INCLUDED 
C	ELASTIC CRACKING AND TENSION STIFFENING
C	==================================================================
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C
C	==================================================================
C	OUTPUT VARIABLE
C	==============
C	YIELD = CURRENT YIELD FUNCTION
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION STSPS(5),STNPS(5),STNCR(8),EFFST(1),SIGMA(5),
	1		  DIREC(2),EPSTN(2),GRTST(2),STNZ(5),TRAMX(5,5),
	2		  STRSP(5),STRNP(5),SGTOT(5),DMATX(5,5),HDMAT(100,2)



	IF(NSTAT.EQ.2.OR.NSTAT.EQ.3) GO TO 30
	IF(NSTAT.EQ.1) GO TO 5

C	------------------------------------------------------------------
C	BACK UP DATA FOR YIELDED CONCRETE
C	------------------------------------------------------------------
	EPSTN(2) = EPSTN(1)
	DIREC(2) = EFFST(1)

5	EPSTN(1) = 0.0
	IF(DIREC(1).NE.0.0) GO TO 25

C	------------------------------------------------------------------
C	FOR THE FIRST CRACKING COMPUTE DIRECTION OF THE CRACK
C	------------------------------------------------------------------
	DUMMY = SGMAX - SIGMA(2)
	IF(DUMMY.EQ.0.0) DUMMY = 0.1E-20
	DIREC(1) = ATAN(SIGMA(3)/DUMMY)

25	NSTAT = 2
30    ANGLE = DIREC(1)
	CALL	TRANS(TRAMX,ANGLE)
	C = COS(ANGLE)
	S = SIN(ANGLE)
C	COMPUTE TWO PRINCIPAL STRESSES
	STRSP(1) = C*C*SIGMA(1) + S*S*SIGMA(2) + 2.0*C*S*SIGMA(3)
	STRSP(2) = SIGMA(1) + SIGMA(2) - STRSP(1)

	TENST = PROPS(LPROP,5)
	UNIAX = PROPS(LPROP,6)

C	------------------------------------------------------------------
C	TEST FOR OPENING OF CRACKS IN TWO DIRECTION
C	------------------------------------------------------------------
	IF(NSTAT.EQ.2.AND.STRSP(2).GE.TENST) NSTAT = 3

	DO 40 I = 1,5
	STRNP(I) = 0.0
	DO 40 J = 1,5
40	STRNP(I) = STRNP(I) + TRAMX(I,J)*STNZ(J)

	CONSE = PROPS(LPROP,9)
	TENSN = PROPS(LPROP,8)
	STRA1 = STRNP(1)
	STRA2 = STRNP(2)
	THETA = 0.0
	CALL	MODUL(DMATX,LPROP,NSTAT,THETA,
	1			  STRA1,STRA2,	0.0,	0.0)

	DO 50 I = 1,5
	SGTOT(I) = 0.0
	DO 50 J = 1,5
50	SGTOT(I) = SGTOT(I) + DMATX(I,J)*STRNP(J)

C	------------------------------------------------------------------
C	MODIFY THE STRESSES CORRESPONDING TO TENSION STIFFENING	
C	------------------------------------------------------------------
	IF(STRA1.GE.TENSN) GO TO 52
	IF(ABS(STRA1).GE.GRTST(1)) GO TO 51
	DUMM1 = CONSE*TENST*(1.0-GRTST(1)/TENSN)
	IF(DUMM1.LT.0.0) DUMM1=0.0
	SGTOT(1) = DUMM1*STRA1/GRTST(1)
	GO TO 52
51	SGTOT(1) = CONSE*TENST*(1.0-STRA1/TENSN)

52	IF(NSTAT.EQ.2) GO TO 55
	IF(STRA2.GE.TENSN) GO TO 55
	IF(ABS(STRA2).GE.GRTST(2)) GO TO 54
	DUMM2 = CONSE*TENST*(1.0-GRTST(2)/TENSN)
	IF(DUMM2.LT.0.0) DUMM2=0.0
	SGTOT(2) = DUMM2*STRA2/GRTST(2)
	GO TO 55
54	SGTOT(2) = CONSE*TENST*(1.0-STRA2/TENSN)
	
55    CONTINUE

C	UPDATE THE MAXIMUM TENSILE STRAIN
	IF(STRA1.GT.GRTST(1)) GRTST(1) = STRA1
	IF(NSTAT.EQ.3.AND.STRA2.GT.GRTST(2)) GRTST(2) = STRA2

C	STORE EFFECTIVE STRESS & STRAIN
	EPSTN(1) = STRA1
	IF(NSTAT.EQ.3) EPSTN(2) = STRA2
	EFFST(1) = SGTOT(1)

	IF(NSTAT.EQ.2) GO TO 72

C	CHECK FOR CRACK CLOSE IN 2 DIRECTION
	IF(STRA2.GT.0.0) GO TO 72

C	IF THE CRACK IS CLOSE CHANGE THE MATERIAL NUMBER
	EPSTN(2) = 0.0
	NSTAT = 2

C	CHECK FOR CRACK CLOSE IN 1 DIRECTION
72	IF(STRA1.GT.0.0) GO TO 76
	IF(NSTAT.EQ.2) GO TO 74

C	IF THE CRACK IS CLOSE IN 1 DIRECTION 
C	BUT STILL OPEN IN 2 DIRECTION 
C	CHANGE THE DIRECTION OF CRACK AND MATERIAL NUMBER
	EPSTN(1) = EPSTN(2)
	EPSTN(2) = 0.0
	GAMMA = 1.570796327
	IF(DIREC(1).GE.0.0) GAMMA = -GAMMA
	DIREC(1) = DIREC(1) + GAMMA
	NSTAT = 2
	GO TO 76

74	EPSTN(1) = 0.0
	NSTAT = 1
	IF(EPSTN(2).EQ.0.0) GO TO 76
	EPSTN(1) = EPSTN(2)
	EFFST(1) = DIREC(2)

76	CONTINUE

	DO 70 I = 1,5
	SIGMA(I) = 0.0
	DO 70 J = 1,5
70	SIGMA(I) = SIGMA(I) + TRAMX(J,I)*SGTOT(J)

	CALL	INVAR(SIGMA,YIELD)

	IF(MOPTN.EQ.2) THEN
	CALL	HARDN(EPSTN,LPROP,HARDS,YIELS,HDMAT)
	ELSE
	YIELS = UNIAX
	HARDS = 0.0
	ENDIF

	IF(YIELD.GT.YIELS) GO TO 90

	DO 80 I = 1,5
80	STSPS(I) = SIGMA(I)

90	CONTINUE
	
	IF(NSTAT.EQ.3) THEN
	GAMMA = 1.570796327
	IF(DIREC(1).GE.0.0) GAMMA = -GAMMA
	DIREC(2) = DIREC(1) + GAMMA
	ENDIF

	RETURN
	END


C	==================================================================
	SUBROUTINE	STELR(STSPS,NSTAT,DIREC,EFFST,EPSTN,
	1				  DSTSZ,LPROP,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	ELASTIC & PLASTIC BEHAVIOR OF REINFORCEMENT STEEL
C	==================================================================
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION	STSPS(5),DIREC(2),DSTSZ(5),
	1			EFFST(1),EPSTN(2)


	ANGLE = PROPS(LPROP,7)+REFANG
	YOUNG = PROPS(LPROP,1)
	HARDS = PROPS(LPROP,2)
	YOUN2 = YOUNG*HARDS/(YOUNG+HARDS)
	YIELD = PROPS(LPROP,5)
	
      
	C  = COS(ANGLE)
	S  = SIN(ANGLE)
	C1 = C*C
	C2 = S*S
	C3 = 2*C*S
	STREP = C1*DSTSZ(1) + C2*DSTSZ(2) + C3*DSTSZ(3)
	XGASH = STREP
	STRPP = C1*STSPS(1) + C2*STSPS(2) + C3*STSPS(3)
	STOTP = STRPP + STREP

	PREYS = YIELD + HARDS*ABS(EPSTN(1))
	IF(ABS(STRPP).GE.PREYS) GO TO 20
	ESCUR = ABS(STOTP) - PREYS
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(STREP)
	GO TO 30

20	IF(STOTP.GT.0.0.AND.STREP.LT.0.0) GO TO 35
	IF(STOTP.LT.0.0.AND.STREP.GT.0.0) GO TO 35
	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	STRPP = STRPP + REDUC*STREP + RFACT*XGASH*(1-YOUNG/(YOUNG+HARDS))
	EPSTN(1) = EPSTN(1) + RFACT*XGASH/(YOUNG+HARDS)
	NSTAT  = 6
	EPSTN(2) = 0.0
	GO TO 50

35	NSTAT  = 5
40	STRPP = STRPP + STREP
50	STSPS(1) = C1*STRPP
	STSPS(2) = C2*STRPP
	STSPS(3) = 0.5*C3*STRPP
	STSPS(4) = 0.0
	STSPS(5) = 0.0
	EFFST(1) = STRPP

	RETURN
	END


C	==================================================================
	SUBROUTINE	MODUL(DMATX,LPROP,NSTAT,ANGLE,
	1				  STRA1,STRA2,YOUNA,YOUNB)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CONSTRUCT THE STRESS-STRAIN MATRIX ACCORDING TO 
C	MATERIAL IDENTIFICATION NUMBER 
C	==================================================================
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION DMATX(5,5),GASHM(5,5),TRAMX(5,5)


	DO 10 I = 1,5
	DO 10 J = 1,5
10	DMATX(I,J) = 0.0
	
	YOUNG = PROPS(LPROP,1)
	GO TO (1,2,2,1,5,6,7),NSTAT

C	------------------------------------------------------------------
C	ELASTIC & PLASTIC CONCRETE BEHAVIOR
C	------------------------------------------------------------------
1	POISN = PROPS(LPROP,2)
	COEF1 = 5.0/6.0
	CONS1 = YOUNG/(1.0-POISN*POISN)
	CONS2 = YOUNG/(2.0*(1.0+POISN))
	DMATX(1,1) = CONS1
	DMATX(2,2) = CONS1
	DMATX(1,2) = POISN*CONS1
	DMATX(2,1) = POISN*CONS1
	DMATX(3,3) = CONS2
	DMATX(4,4) = COEF1*CONS2
	DMATX(5,5) = COEF1*CONS2
	GO TO 100

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
2	POISN = PROPS(LPROP,2)
	GMODU = YOUNG/(2.0*(1.0+POISN))
	CONSG = 0.25
	GTENS = 0.004
	YOUN2 = YOUNG
	GMOD2 = GMODU*5.0/6.0
	YOUN1 = YOUNA
	GMOD1 = CONSG*GMODU*(1.0-ABS(STRA1)/GTENS)
	IF(GMOD1.LE.0.0) GMOD1=0.0
	IF(NSTAT.EQ.2) GO TO 35
	YOUN2 = YOUNB
	GMOD2 = CONSG*GMODU*(1.0-ABS(STRA2)/GTENS)
	IF(GMOD2.LE.0.0) GMOD1=0.0
	GMOD3 = GMOD1
	IF(GMOD2.LE.GMOD1) GMOD3 = GMOD2
35    DMATX(1,1) = YOUN1
	DMATX(2,2) = YOUN2
	DMATX(3,3) = GMOD1
	IF(NSTAT.EQ.3) DMATX(3,3) = 0.5*GMOD3
	DMATX(4,4) = GMOD1
	DMATX(5,5) = GMOD2
	GO TO 50

C	------------------------------------------------------------------
C	ELASTIC STEEL BEHAVIOR
C	------------------------------------------------------------------
5	DMATX(1,1) = YOUNG
	GO TO 50
C	------------------------------------------------------------------
C	ELASTIC-PLASTIC STEEL BEHAVIOR
C	------------------------------------------------------------------
6     DMATX(1,1) = PROPS(LPROP,2)

C	------------------------------------------------------------------
C	FORM THE TRANSFORMATION MATRIX
C	------------------------------------------------------------------
50    IF(ABS(ANGLE).LT.0.001) GO TO 100
	CALL TRANS(TRAMX,ANGLE)

	DO 60 I = 1,5
	DO 60 J = 1,5
	GASHM(I,J) = 0.0
	DO 60 K = 1,5
60	GASHM(I,J) = GASHM(I,J) + DMATX(I,K)*TRAMX(K,J)

	DO 70 I = 1,5
	DO 70 J = 1,5
	DMATX(I,J) = 0.0
	DO 70 K = 1,5
70	DMATX(I,J) = DMATX(I,J) + GASHM(K,J)*TRAMX(K,I)


C	------------------------------------------------------------------
C	CRUSHING OF CONCRETE
C	------------------------------------------------------------------
7	CONTINUE

100   RETURN
	END


C	==================================================================
	SUBROUTINE	CONSTV(STSPS,GRTST,NSTAT,EPSTN,DIREC,
	1				   LPROP,DMATX,HDMAT,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CONSTRUCT CONSTITUTIVE MATRIX FOR CONCRETE AND REINFORCING STEEL
C	==================================================================
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION	STSPS(5),GRTST(2),DIREC(2),ST(5),
	1			AVECT(5),DVECT(5),DMATX(5,5),EPSTN(2),
	2			HDMAT(100,2)

	ANGLE = DIREC(1)
	IDTCS = PROPS(LPROP,10)

	IF(IDTCS.EQ.1) ANGLE = PROPS(LPROP,7)+REFANG
	IF(IDTCS.EQ.1) GO TO 100

C	------------------------------------------------------------------
C	CONTRIBUTION OF CONCRETE RIGIDITY
C	------------------------------------------------------------------

	CONSE = PROPS(LPROP,9)
	TENSN = PROPS(LPROP,8)
	TENST = PROPS(LPROP,5)
	STRA1 = EPSTN(1)
	STRA2 = EPSTN(2)

C	CHECK FOR YIELDING OF CONCRETE	
	IF(NSTAT.NE.4) GO TO 80

C	CHECK FOR UNLOADING
	IF(EPSTN(1).LE.0.0) GO TO 80

C	------------------------------------------------------------------
C	YIELDING OF CONCRETE
C	------------------------------------------------------------------
	CALL	MODUL(DMATX,LPROP,NSTAT,ANGLE,
	1			  0.0,	0.0,	0.0,	0.0)

	DO 50 I = 1,5
50	ST(I) = STSPS(I)
	CALL	FLOWS(ABETA,AVECT,DVECT,LPROP,ST,EPSTN,HDMAT)

	DO 70 I = 1,5
	DO 70 J = 1,5
70	DMATX(I,J) = DMATX(I,J) - ABETA*DVECT(I)*DVECT(J)
	
	GO TO 90

80	CONTINUE

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
	YOUNA = 0.0
	YOUNB = 0.0
	IF(NSTAT.LT.2.OR.NSTAT.GT.3) GO TO 83
C	YOUNA = 0.0
C	YOUNB = 0.0
	IF(NSTAT.EQ.2) GO TO 81
	IF(STRA2.EQ.GRTST(2)) GO TO 81
	IF(GRTST(2).GT.TENSN) GO TO 81
	IF(GRTST(2).EQ.0.0) GO TO 81
	YOUNB = CONSE*TENST*(1-GRTST(2)/TENSN)/GRTST(2)

81	IF(STRA1.EQ.GRTST(1)) GO TO 83
	IF(GRTST(1).GT.TENSN) GO TO 83
	IF(GRTST(1).EQ.0.0) GO TO 83
	YOUNA = CONSE*TENST*(1-GRTST(1)/TENSN)/GRTST(1)

83	CALL	MODUL(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,YOUNA,YOUNB)

	GO TO 90

C	------------------------------------------------------------------
C	CONTRIBUTION OF STEEL RIGIDITY
C	------------------------------------------------------------------
100	CONTINUE

	CALL	MODUL(DMATX,LPROP,NSTAT,ANGLE,
	1			  0.0,	0.0,	0.0,	0.0)

90	CONTINUE


	RETURN
	END


C	==================================================================
	SUBROUTINE FLOWS(ABETA,AVECT,DVECT,LPROP,ST,EPSTN,HDMAT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	COMPUTE FLOW VECTOR 
C	==================================================================
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION AVECT(5),DMATX(5,5),DVECT(5),ST(5),EPSTN(2),HDMAT(100,2)

	IF(MOPTN.EQ.2) THEN
	CALL	HARDN(EPSTN,LPROP,HARDS,YIELS,HDMAT)
	ELSE
C	YIELS = UNIAX  !NOT USE
	HARDS = 0.0
	ENDIF
	
	BETA = 1.355
	ALPHA = 0.355
	C1 = ALPHA*0.5
	C2 = C1*C1 + BETA
	C3 = 2.0*C1*C1 - BETA
	C4 = BETA*3.0
	AFUNC = SQRT(C2*(ST(1)*ST(1) + ST(2)*ST(2)) + C3*ST(1)*ST(2) + 
	1		C4*(ST(3)*ST(3) + ST(4)*ST(4) + ST(5)*ST(5)))
	AFUNC = 2.0*AFUNC
	AVECT(1) = C1+ (2.0*C2*ST(1) + C3*ST(2))/AFUNC
	AVECT(2) = C1+ (2.0*C2*ST(2) + C3*ST(1))/AFUNC
	AVECT(3) = 2.0*C4*ST(3)/AFUNC
	AVECT(4) = 2.0*C4*ST(4)/AFUNC
	AVECT(5) = 2.0*C4*ST(5)/AFUNC

	CALL	MODUL(DMATX,LPROP,	1,	0.0,
	1			  0.0,	0.0,	0.0,	0.0)

	DO 10 I = 1,5
	DVECT(I) = 0.0
	DO 10 J = 1,5
10	DVECT(I) = DVECT(I) + DMATX(I,J)*AVECT(J)

	DENOM = HARDS
	DO 20 I = 1,5
20	DENOM = DENOM + AVECT(I)*DVECT(I)
	ABETA = 1.0/DENOM

	RETURN
	END


C	==================================================================
	SUBROUTINE INVAR(ST,YIELD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	EVALUATE THE YIELD FUNCTION FOR THE GIVEN STRESSES
C	==================================================================
	DIMENSION ST(5)

	BETA  = 1.355
	ALPHA = 0.355
	C1 = BETA*( ST(1)*ST(1) + ST(2)*ST(2) - ST(1)*ST(2)+
	1	 3.0*(ST(3)*ST(3) + ST(4)*ST(4) + ST(5)*ST(5)) )
	C2 = ALPHA*(ST(1) + ST(2))*0.5
	C3 = SQRT(C2*C2 + C1)
	YIELD = C2 + C3

	RETURN
	END


C	==================================================================
	SUBROUTINE INVA2(SN,CRUSH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	EVALUATE THE CRUSHING FUNCTION FOR THE GIVEN STRAIN
C	==================================================================
	DIMENSION SN(5)

	BETA  = 1.355
	ALPHA = 0.355
	C1 = BETA*( SN(1)*SN(1) + SN(2)*SN(2) - SN(1)*SN(2)+
	1	 0.75*(SN(3)*SN(3) + SN(4)*SN(4) + SN(5)*SN(5)) )
	C2 = ALPHA*(SN(1) + SN(2))*0.5
	C3 = SQRT(C2*C2 + C1)
	CRUSH = C2 + C3

	RETURN
	END


C	==================================================================
	SUBROUTINE PRIST(ST,SGMAX,SGMIN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	EVALUATE THE PRINCIPAL STRESS OF 2D PROBLEM
C	==================================================================
	DIMENSION ST(5)
	GASH1 = (ST(1) + ST(2))*0.5
	GASH2 = (ST(1) - ST(2))*0.5
	GASH3 = SQRT(GASH2*GASH2 + ST(3)*ST(3))
	SGMAX = GASH1 + GASH3
	SGMIN = ST(1) + ST(2) - SGMAX

	RETURN
	END


C	==================================================================
	SUBROUTINE TRANS(TRAMX,ANGLE)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	COMPUTE TRANSFORMATION MATRIX OF 2D PROBLEM
C	==================================================================
	DIMENSION TRAMX(5,5)

	DO 10 I = 1,3
	DO 10 J = 1,3
10    TRAMX(I,J) = 0.0
	C = COS(ANGLE)
	S = SIN(ANGLE)
	TRAMX(1,1) = C*C
	TRAMX(2,2) = C*C
	TRAMX(1,2) = S*S
	TRAMX(2,1) = S*S
	TRAMX(1,3) = C*S
	TRAMX(3,1) = -2.0*C*S
	TRAMX(2,3) = -1.0*C*S
	TRAMX(3,2) = 2.0*C*S
	TRAMX(3,3) = C*C - S*S
	TRAMX(4,4) = C
	TRAMX(4,5) = S
	TRAMX(5,4) = -S
	TRAMX(5,5) = C
	
	RETURN 
	END


C	==================================================================
	SUBROUTINE	HARDN(EPSTN,LPROP,HARDS,YIELS,HDMAT)
	IMPLICIT REAL*8	(A-H,O-Z)
	IMPLICIT INTEGER*4	(I-N)
C	==================================================================
C	COMPUTE HARDENING PARAMETER AND CURRENT YIELD STRENGTH
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION HDMAT(100,2),EPSTN(2)

	UNIAX = PROPS(LPROP,6)
	YOUNG = PROPS(LPROP,1)
	UNISN = PROPS(LPROP,7)

	EMAXS = 2.0*UNIAX/YOUNG
	EPMIN = 0.02668*UNIAX/YOUNG
	EPMAX = EMAXS - (UNIAX/YOUNG)
	DEP   = (EPMAX-EPMIN)/(NPONT-1)

	IF(EPSTN(1).LE.0.0) THEN
	HARDS = 0.0
	YIELS = 0.3*UNIAX
	RETURN
	ENDIF

	DO 300 I   = 1,NPONT
	IF(EPSTN(1).GT.HDMAT(I,1)) GO TO 300
	HARDS = (HDMAT(I,2) - HDMAT(I-1,2)) 
	1		/(HDMAT(I,1) - HDMAT(I-1,1))
	YIELS = HDMAT(I-1,2) + HARDS*(EPSTN(1)-HDMAT(I-1,1))
	IF(HARDS.LT.0.0)   HARDS = 0.0
	IF(YIELS.GT.UNIAX) YIELS = UNIAX
	RETURN
300   CONTINUE

	HARDS = 0.0
	YIELS = UNIAX
C	CONST = (0.8*UNIAX - UNIAX)/( (UNISN-2.0*UNIAX/YOUNG)-
C	1		EPMAX)
C	YIELS = UNIAX + CONST*(EPSTN(1)-EPMAX)

	RETURN
	END


C	==================================================================
	SUBROUTINE	HDATM(HDMAT,LPROP)
	IMPLICIT REAL*8	(A-H,O-Z)
	IMPLICIT INTEGER*4	(I-N)
C	==================================================================
C	GENERATE THE PLASTIC STRAIN-STRESS CURVE
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION HDMAT(100,2)

	UNIAX = PROPS(LPROP,6)
	YOUNG = PROPS(LPROP,1)

	EMAXS = 2*UNIAX/YOUNG
	EPMIN = 0.02668*UNIAX/YOUNG
	EPMAX = EMAXS - (UNIAX/YOUNG)
	DEP   = (EPMAX-EPMIN)/(NPONT-1)
	
	DO 100 I   = 1,NPONT
	HDMAT(I,1) = 0.0
100	HDMAT(I,2) = 0.0

	DO 200 I   = 1,NPONT
	EP = EPMIN + DEP*(I-1)
	HDMAT(I,1) = EP - EPMIN
200	HDMAT(I,2) = -YOUNG*EP + 2*SQRT(YOUNG*UNIAX*EP)


	RETURN
	END


C	==================================================================
	SUBROUTINE	CONFAC(STSPS,DSTSZ,PREYS,RFACT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	COMPUTE THE FACTOR INDICATE PORTION OF ELASTIC-PLASTIC (RATIO) 
C	==================================================================
	DIMENSION STSPS(5),DSTSZ(5)


	ACONS = 3.0*DSTSZ(3)*DSTSZ(3) + 3.0*DSTSZ(4)*DSTSZ(4) +
	1		3.0*DSTSZ(5)*DSTSZ(5) + 1.355*DSTSZ(1)*DSTSZ(1) +
	2		1.355*DSTSZ(2)*DSTSZ(2) - 1.355*DSTSZ(1)*DSTSZ(2)

	BCONS = 2.71*STSPS(1)*DSTSZ(1) + 2.71*STSPS(2)*DSTSZ(2) +
	1		6.0*STSPS(3)*DSTSZ(3) + 0.355*PREYS*(DSTSZ(1)+DSTSZ(2))+
	2		6.0*STSPS(5)*DSTSZ(5) + 6.0*STSPS(4)*DSTSZ(4) +
	3		-1.355*STSPS(1)*DSTSZ(2) -1.355*STSPS(2)*DSTSZ(1)

	CCONS = 1.355*STSPS(1)*STSPS(1) + 0.355*PREYS*(STSPS(1)+STSPS(2))
	1		+ 3.0*STSPS(5)*STSPS(5) + 1.355*STSPS(2)*STSPS(2) +
	2		3.0*STSPS(3)*STSPS(3) + 3.0*STSPS(4)*STSPS(4) -
	3		1.355*STSPS(1)*STSPS(2) - 1.0*PREYS*PREYS

	DCONS = BCONS*BCONS - 4*ACONS*CCONS
	IF(DCONS.LT.0.0) DCONS = 0.0

	RFACT = 1 - (-BCONS + SQRT(DCONS))/2.0/ACONS


	RETURN
	END


C	==================================================================
C	=============== END OF LAYERED MODEL OF RC-SHELL =================
C	==================================================================



C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INSTAT (MWA,MTSET,IEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------------
C     READS, GENERATES THE MATERIAL NUMBER OF GAUSS POINT (MSTAT)
C	PRODUCED BY SONGSAK FEB2004
C	------------------------------------------------------------------
C     WA(NWA,NELE)  = WORKING ARRAY STORING STRESSES AND STRAINS
C     ------------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT

C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
C
      DIMENSION WA(MWA),MTSET(1)


	DO 150 IELE = 1,NELE

	CALL ADREWT(IEG,IELE,WA,'RED')
      
      IPT = 0
	LPROP = MTSET(IELE)
	DO 100 IGR = 1,2
	DO 100 IGS = 1,2
	IPT = IPT + 1
	DO 100 ILYR = 1,MNLYR
	NMATS = MATLR(LPROP,ILYR)

	IF(NMATS.GT.0) THEN
	IDCST = PROPS(NMATS,10)
	IF(IDCST.EQ.0) WA(20*(ILYR-1)+18+NWG*(IPT-1)) = 1.0  !MSTAT
	IF(IDCST.EQ.1) WA(20*(ILYR-1)+18+NWG*(IPT-1)) = 5.0  !MSTAT
	ENDIF

100	CONTINUE

	CALL ADREWT(IEG,IELE,WA,'WRT')
150	CONTINUE
	
C
      RETURN
      END
C
C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE INSTATE (MWA,MTSET,IEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------------
C     READS, GENERATES THE MATERIAL NUMBER OF GAUSS POINT (MSTAT)
C	PRODUCED BY SONGSAK FEB2004
C     ------------------------------------------------------------------
C     WA(NWA,NELE)  = WORKING ARRAY STORING STRESSES AND STRAINS
C     ------------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT

C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
C
      DIMENSION WA(MWA),MTSET(1)

	DO 150 IELE = 1,NELE
	CALL ADREWT(IEG,IELE,WA,'RED')

	IPT = 0
	LPROP = MTSET(IELE)
	DO 100 IGR = 1,2
	DO 100 IGS = 1,2
	IPT = IPT + 1
	DO 100 ILYR = 1,MNLYR
	NMATS = MATLR(LPROP,ILYR)

	IF(NMATS.GT.0) THEN
	IDCST = PROPS(NMATS,10)
	IF(IDCST.EQ.0) WA(43*(ILYR-1)+30+NWG*(IPT-1)) = 1.0  !MSTAT
	IF(IDCST.EQ.1) WA(43*(ILYR-1)+30+NWG*(IPT-1)) = 4.0  !MSTAT
	ENDIF

100	CONTINUE	

	CALL ADREWT(IEG,IELE,WA,'WRT')
150	CONTINUE


      RETURN
      END
C	==================================================================
C	==================================================================
C	==================================================================
C	==================================================================
C	==================================================================
C	===== DENSITY OF 4 NODE CONCRETE SHELL ELEMENT =====
C	==================================================================
	SUBROUTINE CONMSS(TH,MSET,RHO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================


	LLAYR = MSET
	NLAYR = MATLR(MSET,MNLYR+1) !MNLYR


	RHO = 0.0

C	------------------------------------------------------------------
C	LOOP OVER CONCRETE & STEEL LAYER
C	------------------------------------------------------------------
      DO 200 ILAYR = 1,NLAYR
	LPROP = MATLR(LLAYR,ILAYR)
      IDTCS = PROPS(LPROP,10)
	DENS = PROPS(LPROP,4)



	DZETA = PROPS(LPROP,3)
	DZETG = 0.5*TH*DZETA

C	------------------------------------------------------------------
C	COMPUTE STRAINS AT THE MIDDLE OF LAYER
C	------------------------------------------------------------------

	IF(IDTCS.EQ.1) GO TO 100


C	------------------------------------------------------------------
C	COMPUTE MASS PER UNIT AREA AT THE MIDDLE OF CONCRETE LAYER
C	------------------------------------------------------------------

	RHO = RHO + DENS*DZETG


	GO TO 95

C	------------------------------------------------------------------
C	COMPUTE MASS PER UNIT AREA AT THE MIDDLE OF STEEL LAYER
C	------------------------------------------------------------------
100	CONTINUE

	RHO = RHO + DENS*DZETG

95	CONTINUE


200	CONTINUE

	RHO = RHO/TH


	RETURN
	END

C	==================================================================
C	==================================================================
C	==================================================================
C	==================================================================










C	=====================================================================
      SUBROUTINE MTPROPE (PROPM,MTSET,MMP,IGIDMEM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------
C     READ MATERIAL PROPERTIES AND ASSIGN PROPERTY SET NUMBERS 
C	TO EACH ELEMENTS OF THE CURRENT ELEMENT GROUP
C	---------------------------------------------
C     PROPM(NMP,NMPS)  = MATERIAL PROPERTIES
C     MTSET(NELE)      = MATERIAL PROPERTY SET NUMBERS
C     ------------------------------------------------
      LOGICAL PROMPT,ERROR
C
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /LOGO/ PROMPT,ERROR,ITEST

C	INTRODUCED BY DE SILVA
	COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT

C	RC SHELL COMMON BLOCK BY SONGSAK OCT2005
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT,MSTEP
C	==================================================================

C      DIMENSION A(3,3),BC(3,3),D(3,3),GG(2,2)
      DIMENSION PROPM(MMP,1),MTSET(1),PROP(50),IGIDMEM(nele)

	DIMENSION LSTEEL(200)  !SONGSAK RC-SHELL DEC2006
C     -----------------------------------------------
C     READ MATERIAL PROPERTIES AND SET DEFAUTE VALUES
C     -----------------------------------------------
c	IGIDMEM = 0
      IF (PROMPT)  WRITE (ITO,1000)  NMP,NMPS
      IF (PROMPT)  WRITE (10,1000)  NMP,NMPS
      IMPS = 0

	LSTEEL = 0  !SONGSAK RC-SHELL DEC2006

C	NEXT IF ADDED BY SONGSAK RC SHELL
	IF (ITYPE.EQ.9.AND.MTMOD.EQ.6) GO TO 500

C	=======================================================
C	CONCRETE WITH CRACKING
C	PRODUCED BY SONGSAK DEC2004
C	=======================================================
500	CONTINUE

	READ(ITI,*)
	READ(ITI,*)
	READ(ITI,*)	MOPTN,MLAYR,MMATS

	READ(ITI,*)
	DO 520 I = 1,MLAYR
	MATLR(I,1:MNLYR) = 0
	READ(ITI,*) NLAYER,MATLR(I,1:NLAYER)
520	MATLR(I,MNLYR+1) = NLAYER

	READ(ITI,*)
	DO 530 I = 1,MMATS
530	READ(ITI,*) (PROPS(I,J),J=1,10)


	READ(ITI,*)

	DO IELET = 1,NELE
      READ (ITI,*)  MEMGID,IPTRN
	IGIDMEM(IELET) = MEMGID
	MTSET(IELET) = IPTRN
	ENDDO

	WRITE(ISO,5000) MLAYR,MMATS
	WRITE(ISO,5050)

	DO 561 I = 1,MLAYR
	NLAYER = MATLR(I,MNLYR+1)
	WRITE(ISO,5100)	I,NLAYER
	WRITE(ISO,5200) (J,J=1,NLAYER)
	WRITE(ISO,5201) (MATLR(I,J),J=1,NLAYER)
	NSTEEL = 0
	ZETAZ = 0.0
	DO 570 J = 1,NLAYER
	NMAT  = MATLR(I,J)
	IDTCS = PROPS(NMAT,10)
	IF(IDTCS.EQ.0) THEN
	ZETAZ = PROPS(NMAT,3) + ZETAZ
	ELSE
	NSTEEL = NSTEEL + 1
	PII = 3.141592654

	IF(LSTEEL(NMAT).EQ.0) THEN
	PROPS(NMAT,7) = PROPS(NMAT,7)*PII/180.0
	ENDIF
	LSTEEL(NMAT) = 1

	ENDIF
570	CONTINUE

	WRITE(ISO,5202) ZETAZ,NSTEEL

	IF(ZETAZ.LT.1.98.OR.ZETAZ.GT.2.02) THEN
	WRITE(*,5600) I
	STOP
	ENDIF

561	CONTINUE

	WRITE(ISO,5300) 
	DO 580 I = 1,MMATS
	IF(PROPS(I,10).EQ.1) GO TO 575
	WRITE(ISO,5400) I,(PROPS(I,J),J=1,9)
	GO TO 580
575	WRITE(ISO,5500) I,(PROPS(I,J),J=1,9)
580	CONTINUE

C	=======================================================
C	OUTPUT FORMAT FOR REINFORCED CONCRETE
C	=======================================================

5000	FORMAT(//28X,24(1H*)/28X,1H*,22X,1H*/
     1        28X,24H* REINFORCED CONCRETE  */
     2        28X,24H* MATERIAL INFORMATION */
	3		28X,1H*,22X,1H*/
     3        28X,24(1H*)/,
	3		/23X,'* ELASTO-PLASTIC FRACTURE MODEL *'/
     4	   //14X,'NUMBER OF LAYERED PATTERN. . .MLAYR =',I5/,
	5	   14X,'NUMBER OF MATERIAL TYPE. . . .MMATS =',I5/)
5050	FORMAT(14X,'MATERIAL IDENTIFICATION OF EACH LAYER PATERN
	1. . .MATLR(MLAYR,MNLYR)'/)
C5100	FORMAT(/14X,'PATTERN NO.     =',I5,4X/,
C     1	   14X,'NUMBER OF LAYER =',I5/,
C     2	   14X,'PATERN LAYOUT')
C5200	FORMAT(10X,100I5)
5100	FORMAT(//14X,'PATTERN NO. . . .=',I5,4X/,
     1	   14X,'NUMBER OF LAYER .=',I5/)
5200	FORMAT(14X,'LAYER  NUMBER',2X,100I4/)
5201	FORMAT(14X,'PATERN LAYOUT',2X,100I4)
5202	FORMAT(14X,'TOTAL ZETA(2)',2X,F9.4,5X,
	1	   'NUMBER OF STEEL LAYER =',I5)
5300	FORMAT(//14X,'MATERIAL PROPERTY. . .PROPS(MMATS,10)'/)
5400	FORMAT(14X,'[LAYER TYPE = CONCRETE]'/,
     1	   14X,'MATERIAL NUMBER. . . . . . . . . . . .=',I4,/
	1	   14X,'YOUNG MODULUS. . . . . . . . . . . . .=',E15.6,/
	2	   14X,'POISSON RATIO. . . . . . . . . . . . .=',F10.5,/
	3	   14X,'LAYER THICKNESS. . . . . . . . . . . .=',F10.5,/
	4	   14X,'DENSITY. . . . . . . . . . . . . . . .=',F10.5,/
	5	   14X,'ULTIMATE TENSILE STRENGTH. . . . . . .=',E15.6,/
	6	   14X,'ULTIMATE COMPRESSION STRENGTH. . . . .=',E15.6,/
	7	   14X,'ULTIMATE COMPRESSIVE STRAIN. . . . . .=',F10.5,/
	8	   14X,'TENSION STIFFENING PARAMETER(Em) . . .=',F10.5,/
	9	   14X,'TENSION STIFFENING PARAMETER(C). . . .=',F10.5/)
 
5500	FORMAT(14X,'[LAYER TYPE = STEEL]'/,
	1	   14X,'MATERIAL NUMBER. . . . . . . . . . . .=',I4,/
	1	   14X,'YOUNG MODULUS. . . . . . . . . . . . .=',E15.6,/
	2	   14X,'ELASTO PLASTIC YOUNG MODULUS . . . . .=',E15.6,/
	3	   14X,'LAYER THICKNESS. . . . . . . . . . . .=',F10.5,/
	4	   14X,'DENSITY. . . . . . . . . . . . . . . .=',F10.5,/
	5	   14X,'YIELD STRENGTH OF STEEL. . . . . . . .=',E15.6,/
	6	   14X,'DISTANCE OF CENTROID TO MID SURFACE  .=',F10.5,/
	7	   14X,'ANGLE (STEEL - LOCAL X AXIS) (RAD) . .=',F10.5,/
	8	   14X,'CURRENT NOT USED . . . . . . . . . . .=',F10.5,/
	9	   14X,'CURRENT NOT USED . . . . . . . . . . .=',F10.5/)

5600	FORMAT(//12X,'*SUM OF LAYER',I3,' IS NOT EQUAL TO THICKNESS*'//)
C	=======================================================
C	END OF OUTPUT FORMAT FOR REINFORCED CONCRETE
C	=======================================================



 1000 FORMAT (//18X,24HINPUT OF MATERIAL PROP.(,I2,15H) FOR ANY SET (,I2
     2,          1H)/18X,44(1H-)/20H SET-NUM  PROPERTIES/1X,19(1H-))
 1100 FORMAT (I4,(6E12.5))
 2000 FORMAT (//21X,35HMATERIAL PROPERTIES FOR SET NUMBER ,I2/
     1         21X,37(1H-)//
     2         6H PROP(,I2,32H)  YOUNG'S MODULUS . . . .  YM =,E14.6/
     3           6H PROP(,I2,32H)  POISSON'S RATIO . . . .  PR =,E14.6/
     4             6H PROP(,I2,32H)  YIELD STRESS  . . . . . YLD =,E14.6
     5/              6H PROP(,I2,32H)  STRAIN HARDENING MODULE  HP =,E14
     6.6/             (6H PROP(,I2,32H)  . . . . . . . . . . . . . . =,E
     614.6/))
C***+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 2005	FORMAT (///21X,35HMATERIAL PROPERTIES FOR SET NUMBER ,I2//
     1         21X,37(1H-)/21X,33HMATERIAL PROPERTIES FOR CONCRETE /
     1         21X,37(1H-)//
     2  6H PROP(,I2,36H) POISSON'S RATIO . . . . . . . PR =,E14.6/
     3  6H PROP(,I2,36H) YIELD STRENGTH OF HOOP. . . . FHY=,E14.6/ 
     4  6H PROP(,I2,36H) COMPRESSIVE STRENGTH . . . . .FC =,E14.6/          
     5  6H PROP(,I2,36H) DIAMETER OF HOOP. . . . . . . DS =,E14.6/    
     6  6H PROP(,I2,36H) LONGITUDINAL SPACE OF HOOPS . SH =,E14.6/  
     7  6H PROP(,I2,36H) CONCRETE MODULUS . . . . . . .YMC=,E14.6/
	8//21X,38HMATERIAL PROPERTIES FOR REINFORCEMENT /
     9         21X,37(1H-)//
     2  6H PROP(,I2,36H)  YOUNG'S MODULUS . . . . . . YMS =,E14.6/
     3  6H PROP(,I2,36H)  POISSON'S RATIO . . . . . . PRS =,E14.6/
     4  6H PROP(,I2,36H)  YIELD STRESS . . . . . . . YLDS =,E14.6/
     5  6H PROP(,I2,36H)  STRAIN-HARDENING STRAIN . . SHS =,E14.6/ 
     6  6H PROP(,I2,36H)  ULTIMATE STRAIN . . . . . . UTE =,E14.6/
     7  6H PROP(,I2,36H)  ULTIMATE STRESS . . . . .. .UTS =,E14.6/)
 2010 FORMAT (//21X,35HMATERIAL PROPERTIES FOR SET NUMBER ,I2/
     1         21X,37(1H-)//
     1  ,10X,'LAMINATE LONGITUDUNAL MODULUS     (0)       E11  =',E15.7/
     1  ,10X,'LAMINATE TRANSVERSE MODULUS       (90)      E22  =',E15.7/
     1  ,10x,'LAMINATE TRANSVERSE MODULUS                 E33  =',E15.7/
     1  ,10X,'POISSON RATIO  12                          MU12  =',F10.5/
     1  ,10x,'POISSON RATIO  23                          MU23  =',F10.5/
     1  ,10X,'POISSON RATIO  31                          MU21  =',F10.5/
     1  ,10x,'SHEAR MODULUS                               G12  =',E15.7/
     1  ,10X,'SHEAR MODULUS                               G13  =',E15.7/
     1 ,10X,'SHEAR MODULUS                               G23  =',E15.7/)

 3000 FORMAT (//16X,22HASSIGN MATERIAL SETS (,I2,15H) TO ELEMENTS (
     1         ,I3,1H)/16X,43(1H-)//21H ELEMA ELEMB INC ISET/1X,20(1H-))
 3100 FORMAT (4I5)
C
15000	CONTINUE

      RETURN
      END
C
C=====================================================================









C=====================================================================
      SUBROUTINE MPPNRCH_OLD(PROPM,PROPG,MTSET,IGSET,WA,FIN,MMP,MGP,MWA)
C	'FIN' ADDED TO PREVIOUS LINE BY GILSON - JUL2003 (INT FORCE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PRINTS GLOBAL STRESS RESULTANTS AT THE GAUSS POINTS AS STORED
C     IN WORKING ARRAY WA, EVALUATES AND PRINTS YIELD FUNCTION
C	--------------------------------------------------------
C     PROPM(NMP,NMPS)  = MATERIAL PROPERTIES
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     MTSET(NELE)      = MATERIAL SET NUMBERS
C     IGSET(NELE)      = GEOMETRIC SET NUMBERS
C     WA(MWA,NELE)     = WORKING ARRAY STORING STRESSES (STRAINS,FLAG)
C
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /CONC/
C	--------------------------------
C	MOPTN		   = MATERIAL MODEL OPTION
C				     (1 = PERFECT-PLASTICITY,2 = STRAIN-HARDENING)
C	MLAYR		   = NUMBER OF LAYERED PATTERN
C	MMATS		   = NUMBER OF MATERIAL PROPERTY
C	MNLYR          = NUMBER OF LAYER IN ONE SECTION
C	NPONT          = NUMBER OF POINT TO GENERATE 
C					 STRESS - PLASTIC STRAIN CURVE	
C	MATLR(MLAYR,50)= MATERIAL IDENTIFICATION NUMBER FOR EACH
C					 LAYER FROM BOTTOM TO TOP
C	PROPS(MMATS,10)= MATERIAL PROPERTY
C		:FOR CONCRETE MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = POISSON RATIO
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = CONCRETE ULTIMATE TENSILE STRENGTH
C			PROPS(MMATS,6) = CONCRETE ULTIMATE COMPRESSION STRENGTH
C			PROPS(MMATS,7) = CONCRETE ULTIMATE COMPRESSIVE STRAIN
C			PROPS(MMATS,8) = TENSION STIFFENING PARAMETER (Em)
C			PROPS(MMATS,9) = TENSION STIFFENING PARAMETER (alpha)
C			PROPS(MMATS,10)= FLAG FOR CONCRETE MATERIAL = 0
C
C		:FOR STEEL MATEIAL
C			PROPS(MMATS,1) = INITIAL YOUNG MODULUS
C			PROPS(MMATS,2) = ELASTO-PLASTIC YOUNG MODULUS
C			PROPS(MMATS,3) = LAYER THICKNESS EXPRESSED 
C							 IN THE NORMALIZED COORDINATE
C			PROPS(MMATS,4) = MATERIAL DENSITY
C			PROPS(MMATS,5) = YIELD STRENGTH OF STEEL
C			PROPS(MMATS,6) = CURRENT NOT USED = 0
C			PROPS(MMATS,7) = ANGLE BETWEEN THE REINFORCEMENT AND THE 
C							 X'-AXIS (MEASURE ANTICLOCKWISE IN RADIANS)
C			PROPS(MMATS,8) = CURRENT NOT USED = 0
C			PROPS(MMATS,9) = CURRENT NOT USED = 0
C			PROPS(MMATS,10)= FLAG FOR STEEL MATERIAL = 1
C
C     ----------------------------------------------------------------
	CHARACTER*1 YFLAG(20),STATE(2)
C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
	CHARACTER*10 NSNAME(20)
C
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/ A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C	NEXT COMMON ADDED BY GILSON - JUL2003 (INT FORCE)
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

C	NEXT COMMON ADDED BY SONGSAK - DEC2004 (CONCRETE)
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT,MSTEP
C	==================================================================
	COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT
C
      DIMENSION PROPM(MMP,1),PROPG(MGP,1),MTSET(1),IGSET(1),WA(MWA,1)

C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(1),NPOGS(27)

      REAL*8 N,M,MN
C	SUNIL 18/01/01 NEXT LINE REMOVED
C      EQUIVALENCE (APEL,IPEL)
      DATA STATE /'.','P'/
C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
      DATA NSNAME /'   (PLANE ','TRUSS)    ','   (SPACE ','TRUSS)    ',
	1             '       (GR','ID)       ','    (PLANE',' FRAME)   ',
     2             '    (SPACE',' FRAME)   ','  (PLANE M','EMBRANE)  ',
     3             '  (MEMBRAN','E-PLATE)  ','   (PLATE-','BENDING)  ',
     4             '    (SHELL',',DOF=6)   ','    (3-D S','OLID)     '/
C
C	NEXT LINE CHANGED BY GILSON - JUL2003 (INT FORCE)
C      WRITE (ISO,1000) KEG,NAME
	WRITE (ISO,1000) KEG,NSNAME(NAME(1)),NSNAME(NAME(2))
C     ---------------------------------------
C     ELASTIC CASE,IN PLANE STRESS RESULTANTS
C     ---------------------------------------

C     ---------------------------------------
C	NEXT LINE CHANGED BY SONGSAK - DEC2004
	IF(MTMOD.EQ.5) GO TO 800
C     ---------------------------------------	

      IF (MTMOD.GT.2) GOTO 200
      WRITE (ISO,1100)
      DO 100  IEL=1,NELE
      WRITE (ISO,1200) IEL
      I1 = 1
      DO 100  IPT=1,NPT
      I2 = I1+2
      WRITE (ISO,1300) IPT,(WA(I,IEL), I=I1,I2)
 100  I1 = I1+NWG
      GOTO 300
C     ----------------------------------------------
C     ELASTO-PLASTIC CASE,IN PLANE STRESS RESULTANTS
C     ----------------------------------------------
 200  IF (MTMOD.GT.3)  GOTO 500
      WRITE (ISO,2100)
      DO 290  IEL=1,NELE
      WRITE (ISO,2200)  IEL
      MSET = MTSET(IEL)
      ISET = IGSET(IEL)
      YLD  = PROPM(3,MSET)
      TH   = PROPG(2,ISET)
      I1 = 1
      DO 290  IPT=1,NPT
      I2 = I1+2
C	SUNIL 18/01/01 NEXT LINE REMOVED AND 2ND LINE MODIFIED
C      APEL = WA(I1+16,IEL)
      IPEL = WA(I1+16,IEL)
      CALL YIELDF (WA(I1,IEL),WA(I1+1,IEL),WA(I1+2,IEL),WA(I1+3,IEL),
     1             WA(I1+4,IEL),WA(I1+5,IEL),UYN,UYM,UYN2,UYM2,
     2             N,M,MN,Q,R,S,FT,2)
      WRITE (ISO,2300) IPT,STATE(IPEL),(WA(I,IEL),I=I1,I2),FT
 290  I1 = I1+NWG
C     --------------------------------------------------
C     ELASTIC AND PLASTIC CASE,BENDING STRESS RESULTANTS
C     --------------------------------------------------
 300  WRITE (ISO,3100)
      DO 390  IEL=1,NELE
      WRITE (ISO,3200) IEL
      I1 = 4
      DO 390  IPT=1,NPT
      I2 = I1+4
      WRITE (ISO,3300) IPT,(WA(I,IEL),I=I1,I2)
 390  I1 = I1+NWG
C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
	GOTO 700
      RETURN
C     --------------------------
C     ELASTO-PLASTIC MULTI LAYER
C     --------------------------
 500  IFPRI = IFPR(8)
 505  IF (IFPRI.LT.2) WRITE (ISO,5000)
      IF (IFPRI.GE.2) WRITE (ISO,6000)
      DO 510  I=1,20
 510  YFLAG(I) = ' '
      DO 590  IEL=1,NELE
      WRITE (ISO,2200) IEL
      I1 = 1
      DO 590  IPT=1,NPT
      DO 520  IGT=1,NGT
C	SUNIL 18/01/01 NEXT LINE REMOVED AND 2ND LINE MODIFIED
C      APEL = WA(I1+8*IGT-1,IEL)
      IPEL = WA(I1+8*IGT-1,IEL)
 520  YFLAG(IGT) = STATE(IPEL)
      IF (IFPRI.GE.2) GOTO 600
      K = NGT-1
      IF (IFPRI.EQ.1) K=1
      IF (NGT.EQ.1) K=1
      DO 550  IGT=1,NGT,K
      I3 = I1 + 8*(IGT-1)
      I4 = I3+2
      EQSIG = DSQRT(WA(I3,IEL)*WA(I3,IEL) + WA(I3+1,IEL)*WA(I3+1,IEL)
     1       - WA(I3,IEL)*WA(I3+1,IEL) + 3.*WA(I4,IEL)*WA(I4,IEL))
      IF (IGT.GT.1) GOTO 530
      WRITE (ISO,5200) IPT,YFLAG,(WA(I,IEL),I=I3,I4),EQSIG
      GOTO 550
 530  WRITE (ISO,5300) (WA(I,IEL),I=I3,I4),EQSIG
 550  CONTINUE
      GOTO 590
C
 600  IA = I1+NWG-9
      IE = IA+7
      WRITE (ISO,6200) IPT,YFLAG,(WA(I,IEL),I=IA,IE)
 590  I1=I1+NWG
      IFPRI = IFPRI-3
      IF (IFPRI.GE.0) WRITE (ISO,6100)
      IF (IFPRI.GE.0) GOTO 505

	GO TO 880

C	=============================================
C	REINFORCED CONCRETE ADDED BY SONGSAK FEB2005
C	=============================================
C	STRESS IN REINFORCED CONCRETE LAYER (MTMOD=5)
C	---------------------------------------------
800	IF(ITASK.EQ.3) GO TO 700 !CONTINUE

	IF(NPT.EQ.4) THEN
	NPOGS(1) = 4
	NPOGS(2) = 2
	NPOGS(3) = 1
	NPOGS(4) = 3
	ELSEIF(NPT.EQ.9) THEN
	NPOGS(1) = 9
	NPOGS(2) = 3
	NPOGS(3) = 1
	NPOGS(4) = 7
	NPOGS(5) = 6
	NPOGS(6) = 2
	NPOGS(7) = 4
	NPOGS(8) = 8
	NPOGS(9) = 5
	ENDIF

	MSTEP = (MSTEP + 1)
	JSTEP = MSTEP*NPRIN


	NPRINT = 0

805	CONTINUE

	NPRINT = NPRINT + 1

	IF(NPRINT.EQ.1) THEN
	WRITE (101,8000) NPT
	WRITE (101,8100) JSTEP
	ELSEIF(NPRINT.EQ.2) THEN
	WRITE (101,9100) JSTEP
	ENDIF


	DO 870 IELE  = 1,NELE
C	IGASH = 0
	IEOUT = MTSET(IELE + NELE) !s
	LPROP = MTSET(IELE)
	NLAYER = MATLR(LPROP,MNLYR+1)

	DO 870 IPT = 1,NPT

	IGASH = ((MNLYR)*(NPOGS(IPT)-1))

	IF(IELE.NE.IEOUT) GO TO 810 !s
	WRITE(ISO,7700)
	WRITE(ISO,7500) IELE,NPOGS(IPT),JSTEP
	WRITE(ISO,7800)
810	CONTINUE !s
	DO 870 INLYR = 1,MNLYR
	IF(IELE.NE.IEOUT) GO TO 870 !s

c	WA(I,ILAYR)  = STNPS(I) I = 1-5
c	WA(J,ILAYR)  = STSPS(I) J = 6-10
c	WA(11,ILAYR) = DIREC(1)
c	WA(12,ILAYR) = DIREC(2)
c	WA(13,ILAYR) = EFFST(1)
c	WA(14,ILAYR) = EPSTN(1)
c	WA(15,ILAYR) = EPSTN(2)
c	WA(16,ILAYR) = GRTST(1)
c	WA(17,ILAYR) = GRTST(2)
c	WA(18,ILAYR) = MSTAT

	IGASH = IGASH+1	
	DUMM1 = 20*IGASH - 12
	DUMM2 = 20*IGASH - 11
	DUMM3 = 20*IGASH - 10
	DUMM4 = 20*IGASH - 9
	DUMM5 = 10*IGASH - 8

	DUMM6 = 20*IGASH - 7
	DUMM7 = 20*IGASH - 6
	DUMM8 = 20*IGASH - 0

	DUMM9  = 20*IGASH - 2
	DUMM10 = 20*IGASH - 1

	MSTAT = WA(DUMM8,IELE)
	PII = 3.141569
	ANGLE1 = WA(DUMM6,IELE)*180./PII
	ANGLE2 = WA(DUMM7,IELE)*180./PII

	IF(MSTAT.EQ.5.OR.MSTAT.EQ.6) THEN
	SSTTC = 'STEEL'
	ELSE
	SSTTC = 'CONCRETE'
	ENDIF

	IF(INLYR.LE.NLAYER) THEN
	WRITE(ISO,7600)	INLYR,SSTTC,MSTAT,ANGLE1,ANGLE2,
	1				WA(DUMM1,IELE),WA(DUMM2,IELE),WA(DUMM3,IELE),
     2				WA(DUMM4,IELE),WA(DUMM5,IELE)

	ENDIF


	LAYNO = 1
	IF(INLYR.EQ.LAYNO)THEN  

	IF(NPRINT.EQ.1) THEN

	IF(INLYR.EQ.LAYNO)THEN
	SIGMAZ = 0.0
	IF(IPT.EQ.1) THEN
	WRITE (101,8200) IELE,WA(DUMM1,IELE),WA(DUMM2,IELE),SIGMAZ,
     2				WA(DUMM3,IELE),WA(DUMM4,IELE),WA(DUMM5,IELE)
	ELSE
	WRITE (101,8300) WA(DUMM1,IELE),WA(DUMM2,IELE),SIGMAZ,
     2				WA(DUMM3,IELE),WA(DUMM4,IELE),WA(DUMM5,IELE)
	ENDIF
	ENDIF

	ELSEIF(NPRINT.EQ.2) THEN

	IF(INLYR.EQ.LAYNO)THEN
	SIGMAZ = 0.0
	IF(IPT.EQ.1) THEN
	WRITE (101,9200) IELE,WA(DUMM9,IELE),WA(DUMM10,IELE),SIGMAZ
	ELSE
	WRITE (101,9300)      WA(DUMM9,IELE),WA(DUMM10,IELE),SIGMAZ
	ENDIF
	ENDIF

	ENDIF

	ENDIF



870	CONTINUE

	
	IF(NPRINT.EQ.1) THEN
	WRITE (101,8400)
	ELSEIF(NPRINT.EQ.2) THEN
	WRITE (101,9400)
	ENDIF

	IF(NPRINT.LT.2) GO TO 805

880	CONTINUE
C	=============================================

C	------------------------------------------
C	PRINT NODAL FORCES
C	PROGRAMMED BY GILSON - JUL2003 (INT FORCE)
C	------------------------------------------
700	IF (ITASK.NE.3) RETURN
	WRITE (ISO,7000) KEG,NSNAME(NAME(1)),NSNAME(NAME(2))
      WRITE (ISO,7100)

C     ------------------
C     LOOP OVER ELEMENTS
C     ------------------
	IC = 1
      DO 710 IE=1,NELE
        WRITE (ISO,7200) IE
C	  LOOP OVER NODES
C	  ---------------
	  DO 720 IN = 1,NNO
          WRITE (ISO,7300) IN,FIN(IC),FIN(IC+1),FIN(IC+2)
	    WRITE (ISO,7400) FIN(IC+3),FIN(IC+4),FIN(IC+5)
	    IC = IC+6
720	  CONTINUE
710	CONTINUE

C
 1000 FORMAT (1H1///28X,24(1H*)/28X,1H*,22X,1H*/
     1        28X,24H* GAUSS POINT STRESSES */
     2        28X,19H* FOR ELEMENT GROUP,I3,2H */
     3         28X,2H* ,2A10,2H */28X,1H*,22X,1H*/28X,24(1H*)//)
 1100 FORMAT (13X,43HELEMENT  GAUSS PT.   NORMAL-X    NORMAL-Y
     1        ,12H  SHEAR-XY  /13X,53(1H-)//)
 1200 FORMAT (12X,I6)
 1300 FORMAT (22X,I4,5X,3E12.4)
 2100 FORMAT (48H ELEM. GPT.  STR.STATE    NORMAL-X    NORMAL-Y
     1        ,26H  SHEAR-XY    YIELD FUNCT./1X,79(1H-)//)
 2200 FORMAT (I4)
 2300 FORMAT (4X,I5,5X,A1,7X,3E12.4,3X,E12.4)
 3100 FORMAT (1H1///44H ELEMENT  GAUSS-PT.   MOMENT-X    MOMENT-Y
     1        ,36H  TWIST-XY    SHEAR-XZ    SHEAR-YZ  /1X,79(1H-)//)
 3200 FORMAT (I6)
 3300 FORMAT (11X,I4,4X,5E12.4)
 5000 FORMAT (22H ELEM. GPT.  STR.STATE,14X,20HSIGMA-X1    SIGMA-X2,4X,
     1        20HTAU-X12    EQ.STRESS/1X,79(1H-)//)
 5200 FORMAT (4X,I5,3X,20A1,4E12.4)
 5300 FORMAT (32X,4E12.4)
 6000 FORMAT (22H ELEM. GPT.  STR.STATE,14X,12HNORMAL-X
     1        ,48HNORMAL-Y    SHEAR-XY    MOMENT-X    MOMENT-Y
     2         ,32HTWIST-XY    SHEAR-XZ    SHEAR-YZ/1X,127(1H-)//)
 6100 FORMAT (1H1///)
 6200 FORMAT (4X,I5,3X,20A1,8E12.4)
C
C	NEXT 'FORMAT' LINES ADDED BY GILSON - JUL2003 (INT FORCE)
7000	FORMAT (1H1///28X,27(1H*)/28X,1H*,25X,1H*/
     +        28X,27H*     INTERNAL FORCES     */
     +        28X,27H* LOCAL STRESS RESULTANTS */
     +        28X,20H*  FOR ELEMENT GROUP,I3,4H   */
     +        28X,3H*  ,2A10,4H   */
     +        28X,1H*,25X,1H*/28X,27(1H*)//)
7100	FORMAT (15X,50HELEM    NODE     FORCE-X      FORCE-Y      FORCE-Z/
     +32X,34HMOMENT-X     MOMENT-Y     MOMENT-Z/14X,
     +54(1H-)//)
7200	FORMAT (14X,I5)
7300	FORMAT (21X,I5,3X,3E13.4)
7400	FORMAT (29X,4E13.4)

C	NEXT 'FORMAT' LINES ADDED BY SONGSAK - DEC2004 (CONCRETE)
7500	FORMAT (2X,I5,2X,I5,4X,I5/)
7600	FORMAT (2X,I5,2X,A9,X,I5,4X,F9.3,4X,F9.3,4X,5E12.4)
7700	FORMAT(//2X,23H ELEM     IPT     MSTEP)
7800	FORMAT(2X,24HLAYER   MAT.TYPE   MSTAT,6X,'ANGLE 1',5X,'ANGLE 2'
	1	   ,6X,'SIGMA X',5X,'SIGMA Y',5X,'SIGMA XY',5X,'SIGMA XZ',
	2	   3X,'SIGMA YZ')

8000	FORMAT('GaussPoints "Shell" Elemtype Quadrilateral'/,
	1	   'Number of Gauss Points:',I5/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

8100  FORMAT('Result "Stress"',2x,'"XFinas"',2x,i3,5x,'Matrix',
     1         2x,'OnGaussPoints',' "Shell" '/,
     2 'ComponentNames  "Sigma-r", "Sigma-s", "Sigma-t", "Tau-rs"',
     3 ', "Tau-rt", "Tau-st"'/'Values') 

8200  FORMAT (I5,X,6E15.4)
8300  FORMAT (6X,6E15.4)
8400	FORMAT ('End Values'/)


9100  FORMAT('Result "Crack Strain"',2x,'"XFinas"',2x,i3,5x,'Vector',
     1         2x,'OnGaussPoints',' "Shell" '/,
     2'ComponentNames  "1st-Crack", "2nd-Crack", "------"
     3 ',/'Values') 

9200  FORMAT (I5,X,3E15.4)
9300  FORMAT (6X,3E15.4)
9400	FORMAT ('End Values'/)

      RETURN
      END
C
C=====================================================================
      SUBROUTINE MPPNRCE_OLD(PROPM,PROPG,MTSET,IGSET,WA,FIN,MMP,MGP,MWA)
C	'FIN' ADDED TO PREVIOUS LINE BY GILSON - JUL2003 (INT FORCE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PRINTS GLOBAL STRESS RESULTANTS AT THE GAUSS POINTS AS STORED
C     IN WORKING ARRAY WA, EVALUATES AND PRINTS YIELD FUNCTION
C	--------------------------------------------------------
C     PROPM(NMP,NMPS)  = MATERIAL PROPERTIES
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     MTSET(NELE)      = MATERIAL SET NUMBERS
C     IGSET(NELE)      = GEOMETRIC SET NUMBERS
C     WA(MWA,NELE)     = WORKING ARRAY STORING STRESSES (STRAINS,FLAG)
C     ----------------------------------------------------------------
	CHARACTER*1 YFLAG(20),STATE(2)
C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
	CHARACTER*10 NSNAME(20)
C
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/ A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C	NEXT COMMON ADDED BY GILSON - JUL2003 (INT FORCE)
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

C	NEXT COMMON ADDED BY SONGSAK - DEC2004 (CONCRETE)
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT,MSTEP
C	==================================================================

	COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT
C
      DIMENSION PROPM(MMP,1),PROPG(MGP,1),MTSET(1),IGSET(1),WA(MWA,1)

C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(1),NPOGS(27)

      REAL*8 N,M,MN
C	SUNIL 18/01/01 NEXT LINE REMOVED
C      EQUIVALENCE (APEL,IPEL)
      DATA STATE /'.','P'/
C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
      DATA NSNAME /'   (PLANE ','TRUSS)    ','   (SPACE ','TRUSS)    ',
	1             '       (GR','ID)       ','    (PLANE',' FRAME)   ',
     2             '    (SPACE',' FRAME)   ','  (PLANE M','EMBRANE)  ',
     3             '  (MEMBRAN','E-PLATE)  ','   (PLATE-','BENDING)  ',
     4             '    (SHELL',',DOF=6)   ','    (3-D S','OLID)     '/
C
C	NEXT LINE CHANGED BY GILSON - JUL2003 (INT FORCE)
C      WRITE (ISO,1000) KEG,NAME
	WRITE (ISO,1000) KEG,NSNAME(NAME(1)),NSNAME(NAME(2))
C     ---------------------------------------
C     ELASTIC CASE,IN PLANE STRESS RESULTANTS
C     ---------------------------------------
C	NEXT LINE ADDED BY SONGSAK APR2005
	IF(MTMOD.EQ.6) GO TO 800

      IF (MTMOD.GT.2) GOTO 200
      WRITE (ISO,1100)
      DO 100  IEL=1,NELE
      WRITE (ISO,1200) IEL
      I1 = 1
      DO 100  IPT=1,NPT
      I2 = I1+2
      WRITE (ISO,1300) IPT,(WA(I,IEL), I=I1,I2)
 100  I1 = I1+NWG
      GOTO 300
C     ----------------------------------------------
C     ELASTO-PLASTIC CASE,IN PLANE STRESS RESULTANTS
C     ----------------------------------------------
 200  IF (MTMOD.GT.3)  GOTO 500
      WRITE (ISO,2100)
      DO 290  IEL=1,NELE
      WRITE (ISO,2200)  IEL
      MSET = MTSET(IEL)
      ISET = IGSET(IEL)
      YLD  = PROPM(3,MSET)
      TH   = PROPG(2,ISET)
      I1 = 1
      DO 290  IPT=1,NPT
      I2 = I1+2
C	SUNIL 18/01/01 NEXT LINE REMOVED AND 2ND LINE MODIFIED
C      APEL = WA(I1+16,IEL)
      IPEL = WA(I1+16,IEL)
      CALL YIELDF (WA(I1,IEL),WA(I1+1,IEL),WA(I1+2,IEL),WA(I1+3,IEL),
     1             WA(I1+4,IEL),WA(I1+5,IEL),UYN,UYM,UYN2,UYM2,
     2             N,M,MN,Q,R,S,FT,2)
      WRITE (ISO,2300) IPT,STATE(IPEL),(WA(I,IEL),I=I1,I2),FT
 290  I1 = I1+NWG
C     --------------------------------------------------
C     ELASTIC AND PLASTIC CASE,BENDING STRESS RESULTANTS
C     --------------------------------------------------
 300  WRITE (ISO,3100)
      DO 390  IEL=1,NELE
      WRITE (ISO,3200) IEL
      I1 = 4
      DO 390  IPT=1,NPT
      I2 = I1+4
      WRITE (ISO,3300) IPT,(WA(I,IEL),I=I1,I2)
 390  I1 = I1+NWG
C	NEXT LINE ADDED BY GILSON - JUL2003 (INT FORCE)
	GOTO 700
      RETURN
C     --------------------------
C     ELASTO-PLASTIC MULTI LAYER
C     --------------------------
 500  IFPRI = IFPR(8)
 505  IF (IFPRI.LT.2) WRITE (ISO,5000)
      IF (IFPRI.GE.2) WRITE (ISO,6000)
      DO 510  I=1,20
 510  YFLAG(I) = ' '
      DO 590  IEL=1,NELE
      WRITE (ISO,2200) IEL
      I1 = 1
      DO 590  IPT=1,NPT
      DO 520  IGT=1,NGT
C	SUNIL 18/01/01 NEXT LINE REMOVED AND 2ND LINE MODIFIED
C      APEL = WA(I1+8*IGT-1,IEL)
      IPEL = WA(I1+8*IGT-1,IEL)
 520  YFLAG(IGT) = STATE(IPEL)
      IF (IFPRI.GE.2) GOTO 600
      K = NGT-1
      IF (IFPRI.EQ.1) K=1
      IF (NGT.EQ.1) K=1
      DO 550  IGT=1,NGT,K
      I3 = I1 + 8*(IGT-1)
      I4 = I3+2
      EQSIG = DSQRT(WA(I3,IEL)*WA(I3,IEL) + WA(I3+1,IEL)*WA(I3+1,IEL)
     1       - WA(I3,IEL)*WA(I3+1,IEL) + 3.*WA(I4,IEL)*WA(I4,IEL))
      IF (IGT.GT.1) GOTO 530
      WRITE (ISO,5200) IPT,YFLAG,(WA(I,IEL),I=I3,I4),EQSIG
      GOTO 550
 530  WRITE (ISO,5300) (WA(I,IEL),I=I3,I4),EQSIG
 550  CONTINUE
      GOTO 590
C
 600  IA = I1+NWG-9
      IE = IA+7
      WRITE (ISO,6200) IPT,YFLAG,(WA(I,IEL),I=IA,IE)
 590  I1=I1+NWG
      IFPRI = IFPRI-3
      IF (IFPRI.GE.0) WRITE (ISO,6100)
      IF (IFPRI.GE.0) GOTO 505

	GO TO 880

C	=============================================
C	REINFORCED CONCRETE ADDED BY SONGSAK
C	=============================================
C	STRESS IN REINFORCED CONCRETE LAYER (MTMOD=5)
C	---------------------------------------------
800	IF(ITASK.EQ.3) GO TO 700 !CONTINUE

	IF(NPT.EQ.4) THEN
	NPOGS(1) = 4
	NPOGS(2) = 2
	NPOGS(3) = 1
	NPOGS(4) = 3
	ELSEIF(NPT.EQ.9) THEN
	NPOGS(1) = 9
	NPOGS(2) = 3
	NPOGS(3) = 1
	NPOGS(4) = 7
	NPOGS(5) = 6
	NPOGS(6) = 2
	NPOGS(7) = 4
	NPOGS(8) = 8
	NPOGS(9) = 5
	ENDIF

	MSTEP = (MSTEP + 1)
	JSTEP = MSTEP*NPRIN


	NPRINT = 0

805	CONTINUE

	NPRINT = NPRINT + 1

	IF(NPRINT.EQ.1) THEN
	WRITE (101,8000) NPT
	WRITE (101,8100) JSTEP
	ELSEIF(NPRINT.EQ.2) THEN
	WRITE (101,9100) JSTEP
	ENDIF


	DO 870 IELE  = 1,NELE
C	IGASH = 0
	IEOUT = MTSET(IELE + NELE) !s

	LPROP = MTSET(IELE) ! p
	NLAYER = MATLR(LPROP,MNLYR+1)

	DO 870 IPT = 1,NPT

	IGASH = ((MNLYR)*(NPOGS(IPT)-1))

	IF(IELE.NE.IEOUT) GO TO 810 !s
	WRITE(ISO,7700) 
	WRITE(ISO,7500) IELE,NPOGS(IPT),JSTEP
	WRITE(ISO,7800)
810	CONTINUE !s
	DO 870 INLYR = 1,MNLYR
	IF(IELE.NE.IEOUT) GO TO 870 !s

c	STNPS(I) = WA(I,ILAYR)  I = 1-5
c	STSPS(I) = WA(J,ILAYR)  J = 6-10
c
c	DIREC(1) = WA(11,ILAYR)
c	DIREC(2) = WA(12,ILAYR)
c
c	EPSTN(1) = WA(13,ILAYR)
c	EPSTN(2) = WA(14,ILAYR)
c	EPSTN(3) = WA(15,ILAYR)
c	EPSTN(4) = WA(16,ILAYR)
c	EPSTN(5) = WA(17,ILAYR)
c
c	ELSNP(1) = WA(18,ILAYR)
c	ELSNP(2) = WA(19,ILAYR)
c	ELSNP(3) = WA(20,ILAYR)
c	ELSNP(4) = WA(21,ILAYR)
c	ELSNP(5) = WA(22,ILAYR)
c
c	GRTSN(1) = WA(23,ILAYR)
c	GRTSN(2) = WA(24,ILAYR)
c
c	EQIST(1) = WA(25,ILAYR)
c	EQIST(2) = WA(26,ILAYR)
c
c	EQISN(1) = WA(27,ILAYR)
c	EQISN(2) = WA(28,ILAYR)
c	EQISN(3) = WA(29,ILAYR)
c
c	MSTAT    = WA(30,ILAYR)
c
c	DEO2     = WA(31,ILAYR)
c	DEMX2    = WA(32,ILAYR)
c
c	WOKPR(1) = WA(33,ILAYR)
c	WOKPR(2) = WA(34,ILAYR)
c	WOKPR(3) = WA(35,ILAYR)
c
c	EQCRK(1) = WA(36,ILAYR)
c	EQCRK(2) = WA(37,ILAYR)
c	EQCRK(3) = WA(38,ILAYR)
c	EQCRK(4) = WA(39,ILAYR)
c	EQCRK(5) = WA(40,ILAYR)
c	EQCRK(6) = WA(41,ILAYR)


	IGASH = IGASH+1 
	DUMM1 = 43*IGASH - 35
	DUMM2 = 43*IGASH - 34
	DUMM3 = 43*IGASH - 33
	DUMM4 = 43*IGASH - 32
	DUMM5 = 43*IGASH - 31
c	DUMM1 = 43*IGASH - 28
c	DUMM2 = 43*IGASH - 27
c	DUMM3 = 43*IGASH - 26
c	DUMM4 = 43*IGASH - 25
c	DUMM5 = 43*IGASH - 24 

	DUMM6 = 43*IGASH - 30
	DUMM7 = 43*IGASH - 29
	DUMM8 = 43*IGASH - 11

c	DUMM10 = 43*IGASH - 28 !p
c	DUMM9  = 43*IGASH - 8  !c

	DUMM9  = 43*IGASH - 18
	DUMM10 = 43*IGASH - 17 

	PII = 3.141569
	ANGLE1 = WA(DUMM6,IELE)*180./PII
	ANGLE2 = WA(DUMM7,IELE)*180./PII
	MSTAT  = WA(DUMM8,IELE)

	IF(MSTAT.EQ.4.OR.MSTAT.EQ.5) THEN
	SSTTC = 'STEEL'
	ELSE
	SSTTC = 'CONCRETE'
	ENDIF

	IF(INLYR.LE.NLAYER) THEN
	WRITE(ISO,7600)	INLYR,SSTTC,MSTAT,ANGLE1,ANGLE2,
	1				WA(DUMM1,IELE),WA(DUMM2,IELE),WA(DUMM3,IELE),
     2				WA(DUMM4,IELE),WA(DUMM5,IELE)
	ENDIF


	LAYNO = 1
	IF(INLYR.EQ.LAYNO)THEN  

	IF(NPRINT.EQ.1) THEN

	IF(INLYR.EQ.LAYNO)THEN
	SIGMAZ = 0.0
	IF(IPT.EQ.1) THEN
	WRITE (101,8200) IELE,WA(DUMM1,IELE),WA(DUMM2,IELE),SIGMAZ,
     2				WA(DUMM3,IELE),WA(DUMM4,IELE),WA(DUMM5,IELE)
	ELSE
	WRITE (101,8300) WA(DUMM1,IELE),WA(DUMM2,IELE),SIGMAZ,
     2				WA(DUMM3,IELE),WA(DUMM4,IELE),WA(DUMM5,IELE)
	ENDIF
	ENDIF

	ELSEIF(NPRINT.EQ.2) THEN

	IF(INLYR.EQ.LAYNO)THEN
	SIGMAZ = 0.0
	IF(IPT.EQ.1) THEN
	WRITE (101,9200) IELE,WA(DUMM9,IELE),WA(DUMM10,IELE),SIGMAZ
	ELSE
	WRITE (101,9300)      WA(DUMM9,IELE),WA(DUMM10,IELE),SIGMAZ
	ENDIF
	ENDIF

	ENDIF

	ENDIF



870	CONTINUE
C	=============================================


	IF(NPRINT.EQ.1) THEN
	WRITE (101,8400)
	ELSEIF(NPRINT.EQ.2) THEN
	WRITE (101,9400)
	ENDIF

	IF(NPRINT.LT.2) GO TO 805


880	CONTINUE
C	=============================================

C	------------------------------------------
C	PRINT NODAL FORCES
C	PROGRAMMED BY GILSON - JUL2003 (INT FORCE)
C	------------------------------------------
700	IF (ITASK.NE.3) RETURN
	WRITE (ISO,7000) KEG,NSNAME(NAME(1)),NSNAME(NAME(2))
      WRITE (ISO,7100)

C     ------------------
C     LOOP OVER ELEMENTS
C     ------------------
	IC = 1
      DO 710 IE=1,NELE
        WRITE (ISO,7200) IE
C	  LOOP OVER NODES
C	  ---------------
	  DO 720 IN = 1,NNO
          WRITE (ISO,7300) IN,FIN(IC),FIN(IC+1),FIN(IC+2)
	    WRITE (ISO,7400) FIN(IC+3),FIN(IC+4),FIN(IC+5)
	    IC = IC+6
720	  CONTINUE
710	CONTINUE

C
 1000 FORMAT (1H1///28X,24(1H*)/28X,1H*,22X,1H*/
     1        28X,24H* GAUSS POINT STRESSES */
     2        28X,19H* FOR ELEMENT GROUP,I3,2H */
     3         28X,2H* ,2A10,2H */28X,1H*,22X,1H*/28X,24(1H*)//)
 1100 FORMAT (13X,43HELEMENT  GAUSS PT.   NORMAL-X    NORMAL-Y
     1        ,12H  SHEAR-XY  /13X,53(1H-)//)
 1200 FORMAT (12X,I6)
 1300 FORMAT (22X,I4,5X,3E12.4)
 2100 FORMAT (48H ELEM. GPT.  STR.STATE    NORMAL-X    NORMAL-Y
     1        ,26H  SHEAR-XY    YIELD FUNCT./1X,79(1H-)//)
 2200 FORMAT (I4)
 2300 FORMAT (4X,I5,5X,A1,7X,3E12.4,3X,E12.4)
 3100 FORMAT (1H1///44H ELEMENT  GAUSS-PT.   MOMENT-X    MOMENT-Y
     1        ,36H  TWIST-XY    SHEAR-XZ    SHEAR-YZ  /1X,79(1H-)//)
 3200 FORMAT (I6)
 3300 FORMAT (11X,I4,4X,5E12.4)
 5000 FORMAT (22H ELEM. GPT.  STR.STATE,14X,20HSIGMA-X1    SIGMA-X2,4X,
     1        20HTAU-X12    EQ.STRESS/1X,79(1H-)//)
 5200 FORMAT (4X,I5,3X,20A1,4E12.4)
 5300 FORMAT (32X,4E12.4)
 6000 FORMAT (22H ELEM. GPT.  STR.STATE,14X,12HNORMAL-X
     1        ,48HNORMAL-Y    SHEAR-XY    MOMENT-X    MOMENT-Y
     2         ,32HTWIST-XY    SHEAR-XZ    SHEAR-YZ/1X,127(1H-)//)
 6100 FORMAT (1H1///)
 6200 FORMAT (4X,I5,3X,20A1,8E12.4)
C
C	NEXT 'FORMAT' LINES ADDED BY GILSON - JUL2003 (INT FORCE)
7000	FORMAT (1H1///28X,27(1H*)/28X,1H*,25X,1H*/
     +        28X,27H*     INTERNAL FORCES     */
     +        28X,27H* LOCAL STRESS RESULTANTS */
     +        28X,20H*  FOR ELEMENT GROUP,I3,4H   */
     +        28X,3H*  ,2A10,4H   */
     +        28X,1H*,25X,1H*/28X,27(1H*)//)
7100	FORMAT (15X,50HELEM    NODE     FORCE-X      FORCE-Y      FORCE-Z/
     +32X,34HMOMENT-X     MOMENT-Y     MOMENT-Z/14X,
     +54(1H-)//)
7200	FORMAT (14X,I5)
7300	FORMAT (21X,I5,3X,3E13.4)
7400	FORMAT (29X,4E13.4)

C	NEXT 'FORMAT' LINES ADDED BY SONGSAK - DEC2004 (CONCRETE)
7500	FORMAT (2X,I5,2X,I5,4X,I5/)
7600	FORMAT (2X,I5,2X,A9,X,I5,4X,F9.3,4X,F9.3,4X,5E12.4)
7700	FORMAT(//2X,23H ELEM     IPT     MSTEP)
7800	FORMAT(2X,24HLAYER   MAT.TYPE   MSTAT,6X,'ANGLE 1',5X,'ANGLE 2'
	1	   ,6X,'SIGMA X',5X,'SIGMA Y',5X,'SIGMA XY',5X,'SIGMA XZ',
	2	   3X,'SIGMA YZ')
7900	FORMAT(10X,E14.6) ! p
7910	FORMAT(4X,I5) ! p

8000	FORMAT('GaussPoints "Shell" Elemtype Quadrilateral'/,
	1	   'Number of Gauss Points:',I5/,
	2	   'Natural Coordinates: internal'/,
	3	   'end gausspoints')

8100  FORMAT('Result "Stress"',2x,'"XFinas"',2x,i3,5x,'Matrix',
     1         2x,'OnGaussPoints',' "Shell" '/,
     2 'ComponentNames  "Sigma-r", "Sigma-s", "Sigma-t", "Tau-rs"',
     3 ', "Tau-rt", "Tau-st"'/'Values') 

8200  FORMAT (I5,X,6E15.4)
8300  FORMAT (6X,6E15.4)
8400	FORMAT ('End Values'/)


9100  FORMAT('Result "Crack Strain"',2x,'"XFinas"',2x,i3,5x,'Vector',
     1         2x,'OnGaussPoints',' "Shell" '/,
     2'ComponentNames  "1st-Crack", "2nd-Crack", "------"
     3 ',/'Values') 

9200  FORMAT (I5,X,3E15.4)
9300  FORMAT (6X,3E15.4)
9400	FORMAT ('End Values'/)

      RETURN
      END


C
C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE SHSTRS(VR,VS,VT,STNCR,STSFB,PROPM,TH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	DIMENSION STSFB(6,1)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION ZETA(3),STNCR(8),PROPM(1)
	DIMENSION DMATX(6,6),STNZ(6),STSCR(6)

	DMATX(1:6,1:6) = 0.0D0

	ZETA(1:3) = [1.0D0,0.0D0,-1.0D0]

	DO 200 ILAYR = 1,3

	ZETAG = 0.5*TH*ZETA(ILAYR)
C	------------------------------------------------------------------
C	COMPUTE STRAINS AT THE MIDDLE OF LAYER
C	------------------------------------------------------------------
	STNZ(1)  = STNCR(1) + ZETAG*STNCR(4)
	STNZ(2)  = STNCR(2) + ZETAG*STNCR(5)
	STNZ(4)  = STNCR(3) + ZETAG*STNCR(6)
	STNZ(3)  = 0.0D0
	STNZ(5)  = STNCR(7)
	STNZ(6)  = STNCR(8)

	YOUNG = PROPM(1)
	POISN = PROPM(2)

	COEF1 = 5.0/6.0
	CONS1 = YOUNG/(1.0-POISN*POISN)
	CONS2 = YOUNG/(2.0*(1.0+POISN))
	DMATX(1,1) = CONS1
	DMATX(2,2) = CONS1
	DMATX(1,2) = POISN*CONS1
	DMATX(2,1) = POISN*CONS1
	DMATX(3,3) = YOUNG
	DMATX(4,4) = CONS2
	DMATX(5,5) = COEF1*CONS2
	DMATX(6,6) = COEF1*CONS2

	DO 15 I = 1,6
	STSCR(I) = 0.0D0
	DO 15 J = 1,6
15	STSCR(I) = STSCR(I) + DMATX(I,J)*STNZ(J)

C	CALL STSTRAN(STSCR,VR,VS,VT,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS

	STSFB(1:6,ILAYR) = STSCR(1:6)

200	CONTINUE


	RETURN
	END


C	=================================================================
C	=================================================================
C	=================================================================
	SUBROUTINE STSTRAN(STSCR,VR,VS,VT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=============================================================

	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION STSCR(6),STS(3,3),TRANS(3,3)

	STS(1:3,1:3) = 0.0D0

	STS(1,1) = STSCR(1)
	STS(2,2) = STSCR(2)
	STS(3,3) = STSCR(3)

	STS(1,2) = STSCR(4)
	STS(2,1) = STSCR(4)

	STS(1,3) = STSCR(5)
	STS(3,1) = STSCR(5)

	STS(2,3) = STSCR(6)
	STS(3,2) = STSCR(6)

	TRANS(1:3,1) =  VR(1:3)
	TRANS(1:3,2) =  VS(1:3)
	TRANS(1:3,3) =  VT(1:3)

	SELECTCASE(IND)
	CASE(0)   !FROM LOCAL  TO GLOBAL
	STS = MATMUL(TRANS,MATMUL(STS,TRANSPOSE(TRANS)))
	CASE(1)   !FROM GLOBAL TO LOCAL
	STS = MATMUL(TRANSPOSE(TRANS),MATMUL(STS,TRANS))
	ENDSELECT

	STSCR(1) =  STS(1,1)
	STSCR(2) =  STS(2,2)
	STSCR(3) =  STS(3,3)
	STSCR(4) =  STS(1,2)
	STSCR(5) =  STS(1,3)
	STSCR(6) =  STS(2,3)
	

	RETURN
	END


C	=================================================================
C	=================================================================
C	=================================================================