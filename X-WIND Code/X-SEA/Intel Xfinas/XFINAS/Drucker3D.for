C	===================================================================
C	===================================================================
	SUBROUTINE DRUG3D(SIGP,EPSP,EPSTN,SIG,EPS,DEP,PROPM,NSIG,
	1				  IFLAG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	===================================================================
	COMMON /DRUGP/ ALPI1
	COMMON /DRSOL/ COH,PHI,IIP

	COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION SIGP(6),EPSP(6),DEP(6,6)
	DIMENSION SIG(6),EPS(6)
	DIMENSION AVECT(6),DVECT(6)
	DIMENSION DSIG(6),DEPS(6)
	DIMENSION PROPM(*),DEPM(6,6)
	
	PI  = 3.141592654 

	YOUNG = PROPM(1)
	POISN = PROPM(2)

	UNIAX = PROPM(3)
	HARDP = PROPM(4)

	ALPI1 = PROPM(6) !ALPHA FOR DRUCKER-PRAGER

	COH  = PROPM(10)
	ANG  = PROPM(11)
	IIP  = PROPM(12) !0 = INNER, 1 = OUTTER
	PHI = ANG*PI/180.
	

	CALL HARPM(EPSTN,UNIAX,HARDP,YIELD) !!!

	DO 10 I = 1,6
10	DEPS(I) = EPS(I) - EPSP(I)

	CALL DMATX(YOUNG,POISN,DEPS,DSIG,6)

	CALL MODUE(YOUNG,POISN,DEP,6)

	AFIRT   = 0.0
	DO 20 I = 1,6
	SIG(I)  = SIGP(I) + DSIG(I)
20	AFIRT   = AFIRT + SIG(I)

	IF(AFIRT.EQ.0.0) GO TO 300 !FIRST ENTRY STEP

	ASCND   = 0.0
	DO I = 1,6
	ASCND   = ASCND + SIGP(I)
	ENDDO

	IF(ASCND.EQ.0.0) THEN !SECOND ENTRY STEP
	PREY = 0.0
	ELSE
	CALL YSURF(SIGP,PREY,AVECT,6)
	ENDIF

	CALL YSURF(SIG ,CURY,AVECT,6)


	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP

	


C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,6
30	SIG(I) = SIGP(I)

	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RFACT(SIG,DSIG,YIELD,RFAC,6)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL STNUM(CURY,YIELD,NSTEP,ASTEP)

	DO 40 I = 1,6
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,6
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP


	DO 200 ISTEP = 1,NSTEP
	
	CALL YSURF(SIG,CURY,AVECT,6)


	CALL DMATX(YOUNG,POISN,AVECT,DVECT,6)

	ABETA = 0.0
	ADUMM = 0.0
	DO I = 1,6
	ADUMM = ADUMM + AVECT(I)*DSIG(I)
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DLAMD = ADUMM*ABETA
	IF(DLAMD.LT.0.0) DLAMD = 0.0

	BDUMM = 0.0
	DO 100 I = 1,6
	BDUMM = BDUMM + AVECT(I)*SIG(I)
	SIG(I) = SIG(I) + DSIG(I)
	1					 - DLAMD*DVECT(I)

100	CONTINUE

	EPSTN  = EPSTN  + DLAMD*BDUMM/CURY

	CALL HARPM(EPSTN,UNIAX,HARDP,YIELD) !!!


200	CONTINUE

	CALL YSURF(SIG,CURY,AVECT,6)

	BRING = 1.0
	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 250 I = 1,6
250	SIG(I) = BRING*SIG(I)


	CALL YSURF(SIG,CURY,AVECT,6)

	CALL DMATX(YOUNG,POISN,AVECT,DVECT,6)

	ABETA = 0.0
	DO I = 1,6
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DO 280 I = 1,6
	DO 280 J = 1,6
280	DEP(I,J) = DEP(I,J) - DVECT(I)*DVECT(J)*ABETA


C	======================================================================
C	END OF ELASTIC-PLASTIC BEHAVIOR
C	======================================================================	
	
300	CONTINUE
	
	DO 350 I = 1,6
	EPSP(I) = EPS(I)
350	SIGP(I) = SIG(I)


	RETURN

	END

C	======================================================================	
C	======================================================================	
C	======================================================================	
	SUBROUTINE DMATX(YOUNG,POISN,EP,SG,NSIG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SG(6),EP(6),STIFG(6,6)

C	=======================================================================


	DO 110 I = 1,6
	DO 110 J = 1,6
110	STIFG(I,J) = 0.0

	STIFG(1,1) = 1.0 - POISN
	STIFG(2,2) = 1.0 - POISN
	STIFG(3,3) = 1.0 - POISN
	STIFG(1,2) = POISN
	STIFG(1,3) = POISN
	STIFG(2,3) = POISN
	STIFG(2,1) = POISN
	STIFG(3,1) = POISN
	STIFG(3,2) = POISN
	STIFG(4,4) = (1.0 - 2.0*POISN)/2.0
	STIFG(5,5) = (1.0 - 2.0*POISN)/2.0
	STIFG(6,6) = (1.0 - 2.0*POISN)/2.0
		 
	DO 150 I = 1,6
	DO 150 J = 1,6
150	STIFG(I,J) = STIFG(I,J)*YOUNG/(1.0+POISN)/(1.0-2.0*POISN)


	DO 200 I = 1,6
	SG(I) = 0.0
	DO 200 J = 1,6
200	SG(I) = SG(I) + STIFG(I,J)*EP(J)

	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================	
	SUBROUTINE MODUE(YOUNG,POISN,STIFG,NSIG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION STIFG(6,6)

C	=======================================================================

	DO 110 I = 1,6
	DO 110 J = 1,6
110	STIFG(I,J) = 0.0

	STIFG(1,1) = 1.0 - POISN
	STIFG(2,2) = 1.0 - POISN
	STIFG(3,3) = 1.0 - POISN
	STIFG(1,2) = POISN
	STIFG(1,3) = POISN
	STIFG(2,3) = POISN
	STIFG(2,1) = POISN
	STIFG(3,1) = POISN
	STIFG(3,2) = POISN
	STIFG(4,4) = (1.0 - 2.0*POISN)/2.0
	STIFG(5,5) = (1.0 - 2.0*POISN)/2.0
	STIFG(6,6) = (1.0 - 2.0*POISN)/2.0
		 
	DO 150 I = 1,6
	DO 150 J = 1,6
150	STIFG(I,J) = STIFG(I,J)*YOUNG/(1.0+POISN)/(1.0-2.0*POISN)


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE YSURF(SG,FN,AVECT,NSIG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)


	COMMON /DRUGP/ ALPI1


	DIMENSION SG(6),AVECT(6),PST(3)
	DIMENSION AV1(6,1),AV2(6,1),AV3(6,1)
	DIMENSION AM1(6,6),AM2(6,6),AM3(6,6)
	


C	======================================================================	
C	PLASTIC POTENTIAL DERIVATIVE MATRIX ( ASSOCIATE FLOW RULE )
C	======================================================================	
	SX  = SG(1)
	SY  = SG(2)
	SZ  = SG(3)
	TXY = SG(4)
	TXZ = SG(5)
	TYZ = SG(6)

	AV1 = 0.0
	AV2 = 0.0
	AV3 = 0.0
	AVECT  = 0.0

	CALL INVAT(SG,PST,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)

C	FLOW VECTOR COMPONENT {AV1}

	AV1(1,1) = 1.     !SX
	AV1(2,1) = 1.     !SY
	AV1(3,1) = 1.     !SZ
	AV1(4,1) = 0.     !TXY
	AV1(5,1) = 0.     !TXZ
	AV1(6,1) = 0.     !TYZ

C	FLOW VECTOR COMPONENT {AV2}

	PM = (1./3.)*(SX+SY+SZ)

	AV2(1,1) = (SX-PM)/(2.*SQRT(VARJ2))    !SX
	AV2(2,1) = (SY-PM)/(2.*SQRT(VARJ2))    !SY
	AV2(3,1) = (SZ-PM)/(2.*SQRT(VARJ2))    !SZ
	AV2(4,1) = (2.*TXY)/(2.*SQRT(VARJ2))   !TXY
	AV2(5,1) = (2.*TXZ)/(2.*SQRT(VARJ2))   !TXZ
	AV2(6,1) = (2.*TYZ)/(2.*SQRT(VARJ2))   !TYZ


C	FLOW VECTOR COMPONENT {AV3}

	AV3(1,1) = (SY-PM)*(SZ-PM)-TYZ*TYZ+VARJ2/3.
	AV3(2,1) = (SX-PM)*(SZ-PM)-TXZ*TXZ+VARJ2/3.
	AV3(3,1) = (SX-PM)*(SY-PM)-TXY*TXY+VARJ2/3.
	AV3(4,1) = 2.*(TYZ*TXZ-(SZ-PM)*TXY)
	AV3(5,1) = 2.*(TXY*TYZ-(SY-PM)*TXZ)
	AV3(6,1) = 2.*(TXZ*TXY-(SX-PM)*TYZ)

C	======================================================================	
C	END OF PLASTIC POTENTIAL DERIVATIVE MATRIX
C	======================================================================


C	=====================================	
C	YIELD FUNCTION FLOW VECTOR COEFICIENT
C	=====================================

C	=========================
C	VON-MISE
C	=========================
C	FN  = SQRT(VARJ2)
C	CF1 = 0.0
C	CF2 = 1.0
C	CF3 = 0.0
C	=========================
C	END VON-MISE
C	=========================

C	=========================
C	DRUCKER-PRACKER
C	=========================
	FN  = (ALPI1*VARI1) + SQRT(VARJ2)
	CF1 = ALPI1 
	CF2 = 1.0 
	CF3 = 0.0
C	=========================
C	END DRUCKER-PRACKER
C	=========================


C	============================================	
C	END OF YIELD FUNCTION FLOW VECTOR COEFICIENT
C	============================================


C	=========================	
C	FLOW VECTOR
C	=========================
	DO I = 1,6
	AVECT(I) = 0.0
	AVECT(I) = CF1*AV1(I,1) + CF2*AV2(I,1) + CF3*AV3(I,1)
	ENDDO



	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE RFACT(SG,DSG,YIELD,RFAC,NSIG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /DRUGP/ ALPI1

	DIMENSION SG(6),DSG(6)
	DIMENSION SGR(6),AVECT(6)


	NITER = 20
	RFAC1 = 0.0
	RFAC2 = 1.0
	TOL = 0.00001
	RFACO = RFAC2
	DO ITER = 1,NITER
	
	RFAC = 0.5*(RFAC1 + RFAC2)
	DO I = 1,6
	SGR(I) = SG(I) + RFAC*DSG(I)
	ENDDO
	CALL YSURF(SGR,FR,AVECT,6)
	TEST = FR - YIELD
	IF(TEST.LE.0.0) RFAC1 = RFAC
	IF(TEST.GT.0.0) RFAC2 = RFAC
	TEST = (RFAC-RFACO)/RFACO
	IF(ABS(TEST).LE.TOL) GO TO 200
	RFACO = RFAC
	ENDDO

200	CONTINUE

	
	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE STNUM(CURY,YIELD,NSTEP,ASTEP)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	=======================================================================

	NSTEP = 30

	IF(YIELD.NE.0.0D0) THEN
	ESCUR = CURY - YIELD
	NSTEP = ESCUR*8.0/YIELD + 1.0
	ENDIF

	IF(NSTEP.GT.30) NSTEP = 30
	ASTEP = NSTEP


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE HARPM(EPSTN,UNIAX,HARDP,YIELD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /DRUGP/ ALPI1
	COMMON /DRSOL/ COH,PHI,IIP
C	=======================================================================


	YIELD = UNIAX + HARDP*EPSTN

	YIELD = YIELD/SQRT(3.0) ! FOR VON-MISE & DRUCKER-PRAGER

	IF(COH.GT.0.0.OR.PHI.GT.0.0) THEN
	FACT = (3.0+SIN(PHI))				       ! MOHR COULOMB INNER
	IF(IIP.EQ.1) FACT = (3.0-SIN(PHI))         ! MOHR COULOMB OUTTER
	
	FACT1 = (3.0+SIN(PHI))
	FACT2 = (3.0-SIN(PHI))
	FACT  = 0.5*(FACT1 + FACT2)

	YIELD = 6.0*COH*COS(PHI)/ SQRT(3.0) /FACT
	ALPI1 =     2.0*SIN(PHI)/ SQRT(3.0) /FACT  
	ENDIF


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE INVAT(SIGMA,PSTRES,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	======================================================================
C	SUBROUTINE TO COMPUTE THE PRINCIPLE STRESS AND INVARIANT
C	======================================================================

	DIMENSION SIGMA(6),PSTRES(3),STMAT(3,3)



	PI = 3.141592654

	SX  = SIGMA(1)
	SY  = SIGMA(2)
	SZ  = SIGMA(3)
	TXY = SIGMA(4)
	TXZ = SIGMA(5)
	TYZ = SIGMA(6)

	STMAT(1,1) = SX
	STMAT(1,2) = TXY
	STMAT(1,3) = TXZ

	STMAT(2,1) = TXY
	STMAT(2,2) = SY
	STMAT(2,3) = TYZ

	STMAT(3,1) = TXZ
	STMAT(3,2) = TYZ
	STMAT(3,3) = SZ

C	==========================
C	S1 > S3 > S2 = MAX/MIN/MID
C	==========================
	
	VARI1 = 0.0
	VARI2 = 0.0
	VARI3 = 0.0

	DO I = 1,3
	VARI1 = VARI1 + STMAT(I,I)
	DO J = 1,3
	VARI2 = VARI2 + (0.5*STMAT(I,J)*STMAT(I,J))
	DO K = 1,3
	VARI3 = VARI3 + (STMAT(I,J)*STMAT(J,K)*STMAT(K,I)/3.0)
	ENDDO
	ENDDO
	ENDDO

	STMAT(1,1) = STMAT(1,1) - (VARI1/3.0)
	STMAT(2,2) = STMAT(2,2) - (VARI1/3.0)
	STMAT(3,3) = STMAT(3,3) - (VARI1/3.0)

	VARJ1 = 0.0
	VARJ2 = 0.0
	VARJ3 = 0.0

	DO I = 1,3
	VARJ1 = VARJ1 + STMAT(I,I)
	DO J = 1,3
	VARJ2 = VARJ2 + (0.5*STMAT(I,J)*STMAT(I,J))
	DO K = 1,3
	VARJ3 = VARJ3 + (STMAT(I,J)*STMAT(J,K)*STMAT(K,I)/3.0)
	ENDDO
	ENDDO
	ENDDO

C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE

	STEFF = SQRT(VARJ2)

	IF(STEFF.EQ.0.0) GOTO 10

	SINT3 = -3.0*SQRT(3.0)*VARJ3/(2.0*VARJ2*STEFF)

	IF(SINT3.GT.1.0D0) SINT3 = 1.0D0

	GOTO 20

10	SINT3 = 0.0D0

20	CONTINUE
	IF(SINT3.LT.-1.0D0) SINT3 = -1.0D0	
	IF(SINT3.GT. 1.0D0) SINT3 =  1.0D0

	THETA = (1.0/3.0)*ASIN(SINT3)
	
C	...COMPUTE THE PRICIPLE STRESSES

	PSTRES(1)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+2.*PI/3.)+VARI1/3. !MAX		 
	PSTRES(2)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+4.*PI/3.)+VARI1/3. !MIN	 
	PSTRES(3)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA         )+VARI1/3. !MID

C	==========================
C	S1 > S3 > S2 = MAX/MIN/MID
C	==========================

	RETURN

	END

C	======================================================================
C	======================================================================
C	======================================================================
	SUBROUTINE MODUE2D(DEPM,DEP,IFLAG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION DEP(6,6),DEPM(6,6),LDA(4)
	DIMENSION DK2(3),DK3(3)

	LDA(1) = 1
	LDA(2) = 2
	LDA(3) = 4
	LDA(4) = 3

	IF(IFLAG.EQ.1) THEN !PLANE STRAIN

	DEP = 0.0

	DEP(1,1) = DEPM(1,1)
	DEP(1,2) = DEPM(1,2)
	DEP(2,1) = DEPM(2,1)
	DEP(2,2) = DEPM(2,2)
	DEP(3,3) = DEPM(4,4)
	
	ENDIF


	IF(IFLAG.EQ.2) THEN !PLANE STRESS

	DO I = 1,4
	NUMI = LDA(I)
	DO J = 1,4
	NUMJ = LDA(J)
	DEP(I,J) = DEPM(NUMI,NUMJ)
	ENDDO
	ENDDO

	DK2(1) = DEP(1,4)
	DK2(2) = DEP(2,4)
	DK2(3) = DEP(3,4)
	DK3(1) = DEP(4,1)
	DK3(2) = DEP(4,2)
	DK3(3) = DEP(4,3)
	DK4    = DEP(4,4)
	DO I = 1,3
	DO J = 1,3
	DEP(I,J) = DEP(I,J) - DK2(I)*DK3(J)*(1.0/DK4)
	ENDDO
	ENDDO

	
	ENDIF	
	


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================


	SUBROUTINE DRUG2D(SIGPM,EPSPM,EPSTN,SIGM,EPSM,DEP,PROPM,
	1				  STSZ,STNZ,IFLAG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /DRUGP/ ALPI1
	COMMON /DRSOL/ COH,PHI,IIP

	COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

	DIMENSION SIGPM(4),EPSPM(4)
	DIMENSION SIGM(4),EPSM(4)
	DIMENSION SIGP(6),EPSP(6),DEP(6,6)
	DIMENSION SIG(6),EPS(6)
	DIMENSION AVECT(6),DVECT(6)
	DIMENSION DSIG(6),DEPS(6)
	DIMENSION PROPM(*),DEPM(6,6)
	
	PI  = 3.141592654 

	YOUNG = PROPM(1)
	POISN = PROPM(2)

	UNIAX = PROPM(3)
	HARDP = PROPM(4)


	ALPI1 = PROPM(6) !ALPHA FOR DRUCKER-PRAGER


	COH  = PROPM(10)
	ANG  = PROPM(11)
	IIP  = PROPM(12) !0 = INNER, 1 = OUTTER
	PHI  = ANG*PI/180.


	CALL HARPM(EPSTN,UNIAX,HARDP,YIELD) !!!

	EPSTNO = EPSTN


	SIGP(1) = SIGPM(1)
	SIGP(2) = SIGPM(2)
	SIGP(3) = STSZ
	SIGP(4) = SIGPM(3)
	SIGP(5) = 0.0
	SIGP(6) = 0.0

	EPSP(1) = EPSPM(1)
	EPSP(2) = EPSPM(2)
	EPSP(3) = STNZ       !0.0
	EPSP(4) = EPSPM(3)
	EPSP(5) = 0.0
	EPSP(6) = 0.0

	EPS(1)  = EPSM(1)
	EPS(2)  = EPSM(2)
	EPS(3)  = EPSM(4)    !0.0
	EPS(4)  = EPSM(3)
	EPS(5)  = 0.0
	EPS(6)  = 0.0

C	PLANE STRESS
	IF(IFLAG.EQ.2) THEN
	SIGP(3) = 0.0D0 
	EPSP(3) = -POISN*(EPSPM(1) + EPSPM(2))/(1.0-POISN) 
	EPS(3)  = -POISN*( EPSM(1) +  EPSM(2))/(1.0-POISN) 
	ENDIF


	DO 10 I = 1,6
10	DEPS(I) = EPS(I) - EPSP(I)


	CALL DMATX(YOUNG,POISN,DEPS,DSIG,6)

	CALL MODUE(YOUNG,POISN,DEP,6)

	AFIRT   = 0.0
	DO 20 I = 1,6
	SIG(I)  = SIGP(I) + DSIG(I)
20	AFIRT   = AFIRT + SIG(I)

	IF(AFIRT.EQ.0.0) GO TO 300 !FIRST ENTRY STEP

	ASCND   = 0.0
	DO I = 1,6
	ASCND   = ASCND + SIGP(I)
	ENDDO

	IF(ASCND.EQ.0.0) THEN !SECOND ENTRY STEP
	PREY = 0.0
	ELSE
	CALL YSURF(SIGP,PREY,AVECT,6)
	ENDIF

	CALL YSURF(SIG ,CURY,AVECT,6)


	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP


C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,6
30	SIG(I) = SIGP(I)

	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RFACT(SIG,DSIG,YIELD,RFAC,6)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL STNUM(CURY,YIELD,NSTEP,ASTEP)

	DO 40 I = 1,6
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,6
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP


	DO 200 ISTEP = 1,NSTEP
	
	CALL YSURF(SIG,CURY,AVECT,6)

C	PLANE STRESS
	IF(IFLAG.EQ.2) THEN
	AVECT(3)  = -POISN*(AVECT(1) + AVECT(2))/(1.0-POISN) 
	ENDIF

	CALL DMATX(YOUNG,POISN,AVECT,DVECT,6)

	ABETA = 0.0
	ADUMM = 0.0
	DO I = 1,6
	ADUMM = ADUMM + AVECT(I)*DSIG(I)
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DLAMD = ADUMM*ABETA
	IF(DLAMD.LT.0.0) DLAMD = 0.0

	BDUMM = 0.0
	DO 100 I = 1,6
	BDUMM = BDUMM + AVECT(I)*SIG(I)
	SIG(I) = SIG(I) + DSIG(I)
	1					 - DLAMD*DVECT(I)

100	CONTINUE


	EPSTN  = EPSTN  + DLAMD*BDUMM/CURY

	CALL HARPM(EPSTN,UNIAX,HARDP,YIELD) !!!
	

200	CONTINUE

	
	CALL YSURF(SIG,CURY,AVECT,6)



	BRING = 1.0
	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 250 I = 1,6
250	SIG(I) = BRING*SIG(I)


	CALL YSURF(SIG,CURY,AVECT,6)


	CALL DMATX(YOUNG,POISN,AVECT,DVECT,6)

	ABETA = 0.0
	DO I = 1,6
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DO 280 I = 1,6
	DO 280 J = 1,6
280	DEP(I,J) = DEP(I,J) !- DVECT(I)*DVECT(J)*ABETA


C	======================================================================
C	END OF ELASTIC-PLASTIC BEHAVIOR
C	======================================================================	
	
300	CONTINUE
	
	DO 350 I = 1,6
	EPSP(I) = EPS(I)
350	SIGP(I) = SIG(I)


	SIGM(1)  = SIG(1)
	SIGM(2)  = SIG(2)
	SIGM(3)  = SIG(4)
	SIGM(4)  = SIG(3)
	STSZ     = SIG(3)

	SIGPM(1) = SIGP(1)
	SIGPM(2) = SIGP(2)
	SIGPM(3) = SIGP(4)
	SIGPM(4) = SIGP(3)

	EPSPM(1) = EPSP(1)
	EPSPM(2) = EPSP(2)
	EPSPM(3) = EPSP(4)
	EPSPM(4) = EPSP(3) !0.0


	EPSM(1)  = EPS(1)
	EPSM(2)  = EPS(2)
	EPSM(3)  = EPS(4)
	EPSM(4)  = EPS(3)  !0.0
	STNZ     = EPS(3)  !0.0

	DO I = 1,6
	DO J = 1,6
	DEPM(I,J) = DEP(I,J)
	ENDDO
	ENDDO

	CALL MODUE2D(DEPM,DEP,IFLAG)



500	CONTINUE

	RETURN

	END

C	======================================================================
C	======================================================================	
C	======================================================================
	SUBROUTINE POSVMS(SG,FN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SG(6),AVECT(6),PST(3)

	
C	======================================================================	
C	PLASTIC POTENTIAL DERIVATIVE MATRIX ( ASSOCIATE FLOW RULE )
C	======================================================================	
	AV1 = 0.0
	AV2 = 0.0
	AV3 = 0.0
	AVECT  = 0.0

	SIGTO = 0.0D0
	DO I = 1,6
	SIGTO = SIGTO + SG(I)
	ENDDO
	
	VARJ2 = 0.0D0
	IF(SIGTO.NE.0.0D0) THEN
	CALL INVAT(SG,PST,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)
	ENDIF


C	=====================================	
C	YIELD FUNCTION FLOW VECTOR COEFICIENT
C	=====================================

C	=========================
C	VON-MISE
C	=========================
	FN  = SQRT(VARJ2)
C	=========================
C	END VON-MISE
C	=========================


	RETURN

	END	

C	======================================================================	
C	======================================================================	