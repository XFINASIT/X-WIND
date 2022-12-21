	SUBROUTINE	PLVIS(SIGP,EPSP,VVEP,DFVP,
	1			      SIG,EPS,PROPM,BMAT,DVOL,EPSTN,
     2				  DEP)

c	SUBROUTINE	PLVIS(SIGP,EPSP,VVEP,EVP,FVP,DFVP,
c	1			      SIG,EPS,PROPM,BMAT,DVOL,EPSTN,
c    2				  DEP)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
C     3D VISCO-PLASTIC ANALYSIS  (DRUKER-PRAGER CRITERIA)
C	=======================================================================
	COMMON /SEEP/  NTSTEP,KTSTEP,CTIME,DTINC,KFRES
C	=======================================================================
	COMMON /DRUGP/ ALPI1
C	=======================================================================
	DIMENSION SIGP(6),EPSP(6),PST(3),SIG(6),EPS(6)

	DIMENSION DEPS(6),DELPN(6),DSIG(6),DEP(6,6),DEVP(6)

	DIMENSION FVP(24),DFVP(24),VVEP(6),EVP(6),VVEP1(6)

	DIMENSION PROPM(*),BMAT(6,24),DFVPG(24)

	DIMENSION AVECT(6)

	YOUNG = PROPM(1)
	POISN = PROPM(2)
	UNIAX = PROPM(3)
	HARDP = PROPM(4)
	ALPI1 = PROPM(14) !ALPHA FOR DRUCKER-PRAGER
	AMM   = PROPM(15) !FLUIDITY PARAMETER GAMMA
	PMAN  = PROPM(16) !POWER OF FLOW FUNCTION (EQ 8.8-8.9 IN OWEN-HINTON)

C	COH   = PROPM(10)
C	ANG   = PROPM(11)
C	BETA  = PROPM(12)
	
	DO I = 1,6
	DEPS(I)  = EPS(I)  - EPSP(I)
	DELPN(I) = DEPS(I) - (VVEP(I)*DTINC)
	ENDDO


C	=============================================
C	=============================================

	CALL HARPM(EPSTN,UNIAX,HARDP,YIELD) !!!


	CALL DMATX(YOUNG,POISN,DELPN,DSIG,6)

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

	CURY1 = CURY

	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

C	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP (NO UNLOAD IN VISCO PLASTICITY)

	CALL RETPL(SIGP,SIG,DSIG,YIELD,PREY,
	1		   CURY,YOUNG,POISN,EPSTN,DEP,AVECT,UNIAX)

	CALL FNVIS(CURY1,YIELD,AMM,PMAN,3,QF)


	DO I = 1,6
	VVEP1(I) = QF*AVECT(I)
	ENDDO

	GO TO 400

C	=============================================
C	=============================================

300	CONTINUE

	QF = 0.0
	DO I = 1,6
	VVEP1(I) = 0.0
	ENDDO

400	CONTINUE

	
	TFACT = 1.0  !CAN BE CHANGE
	DO I = 1,6
	DEVP(I) = DTINC*( (1.0-TFACT)*VVEP(I) + TFACT*VVEP1(I) )

	EVP(I) = EVP(I) + DEVP(I)

	VVEP(I) = VVEP1(I)
	ENDDO

	
	DFVPG = MATMUL(TRANSPOSE(BMAT),MATMUL(DEP,DEVP))
	DO I = 1,6
	DFVP(I) =  DFVP(I) - DFVPG(I)*DVOL
	FVP(I)  =  FVP(I)  - DFVPG(I)
	ENDDO



C	STORE PREVIOUS STRESS & STRAIN
	DO I = 1,6
	SIGP(I) = SIG(I)
	EPSP(I) = EPS(I)
	ENDDO


	RETURN


	END


C	UPDATE TIME INCREMENTAL ( OWEN HINTON ) !MUST BE PLACE OUTSIDE THE ELEMENT LOOP
C	FAK   = 1.5
C	DTINC = FAK*DTINC

C	LIMITATION OF DELTA T
C	DTINI = (1.0+POISN)*Fy/GAMM/YOUNG                      !FOR TRESCA
C	DTINI = 4.0*(1.0+POISN)*Fy/3.0/GAMM/YOUNG              !FOR VON-MISES
C	DTINI = 4.0*(1.0+POISN)*(1.0-2.0*POISN)*(C*COS(PHI))
C	1        /(1.0-2.0*POISN+SIN(PHI)*SIN(PHI))/GAMM/YOUNG  !FOR MOHR-COULOMB


C	======================================================================
C	======================================================================
	SUBROUTINE RETPL(SIGP,SIG,DSIG,YIELD,PREY,
	1				 CURY,YOUNG,POISN,EPSTN,DEP,AVECT,UNIAX)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	======================================================================
	DIMENSION SIGP(6),SIG(6),DSIG(6),AVECT(6),DVECT(6),DEP(6,6)
C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================
	NPLOP = 0   ! 0 = NO RETURN STRESS , 1 = RETURN STRESS
	IF(NPLOP.EQ.0) GO TO 255  !!

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,6
30	SIG(I) = SIGP(I)

	

	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RFACT(SIG,DSIG,YIELD,RFAC,6)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL STNUM(CURY,YIELD,NSTEP,ASTEP)
C	NSTEP = 100
C	ASTEP = 1.0*NSTEP

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


255	CONTINUE

	CALL YSURF(SIG,CURY,AVECT,6)

	IF(NPLOP.EQ.0) GO TO 300  !!


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

	RETURN

	END

C	======================================================================	
C	======================================================================
C	======================================================================
	SUBROUTINE FNVIS(F,Fo,AMM,PMAN,IOPT,QF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	======================================================================
C	EVALUATE THE TIME DERIVATIVE FUNCTION ( OWEN & HINTON )
C	======================================================================

	RATIO = (F-Fo)/Fo

C	write(*,*) ratio

	IF(IOPT.EQ.1) QF = (EXP(PMAN*RATIO) - 1.0)

	IF(IOPT.EQ.2) QF = (RATIO**PMAN)

	IF(IOPT.EQ.3) QF = F

	QF = AMM*QF



	RETURN

	END
C	======================================================================
C	======================================================================


