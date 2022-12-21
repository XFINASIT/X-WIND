C	===================================================================
C	MATERIAL MODEL USED IN ZERO THICKNESS CONTACT INTERFACE ELEMENT
C	===================================================================
	SUBROUTINE THIN3D(SIGP,EPSP,SIG,EPS,DEP,PROPM,WORKH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	===================================================================
	COMMON /THIN3/ GIIF,COH,PHI,ALPAB

	DIMENSION SIGP(3),EPSP(3),DEP(3,3)
	DIMENSION SIG(3),EPS(3)
	DIMENSION AVECT(3),DVECT(3)
	DIMENSION DSIG(3),DEPS(3)
	DIMENSION PROPM(*)


	YOUNG = PROPM(1)
	POISN = PROPM(2)


	GIIF = PROPM(7)
	HARDP = 0.0

	COH  = PROPM(10)
	ANG  = PROPM(11)

	PI   = 3.141592654 
	PHI  = ANG*PI/180.0
	

	CALL HTHIN(UNIAX,YIELD,WORKH) !!!
	

	DO 10 I = 1,3
10	DEPS(I) = EPS(I) - EPSP(I)


	CALL DTHIN(YOUNG,POISN,DEPS,DSIG)

	CALL MTHIN(YOUNG,POISN,DEP)


	AFIRT   = 0.0
	DO 20 I = 1,3
	SIG(I)  = SIGP(I) + DSIG(I)
20	AFIRT   = AFIRT + SIG(I)

	IF(AFIRT.EQ.0.0) GO TO 300 !FIRST ENTRY STEP

	ASCND   = 0.0
	DO I = 1,3
	ASCND   = ASCND + SIGP(I)
	ENDDO

	IF(ASCND.EQ.0.0) THEN !SECOND ENTRY STEP
	PREY = 0.0
	ELSE
	CALL YTHIN(SIGP,PREY,AVECT)
	ENDIF

	CALL YTHIN(SIG ,CURY,AVECT)


	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP

	


C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,3
30	SIG(I) = SIGP(I)


	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RTHIN(SIG,DSIG,YIELD,RFAC)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL STHIN(CURY,YIELD,NSTEP,ASTEP)


	DO 40 I = 1,3
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,3
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP



	DO 200 ISTEP = 1,NSTEP

	
	CALL YTHIN(SIG ,CURY,AVECT)

	CALL DTHIN(YOUNG,POISN,AVECT,DVECT)

	ABETA = 0.0
	ADUMM = 0.0
	DO I = 1,3
	ADUMM = ADUMM + AVECT(I)*DSIG(I)
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DLAMD = ADUMM*ABETA
	IF(DLAMD.LT.0.0) DLAMD = 0.0


	DO 100 I = 1,3  !BACK TO SURFACE
	SIG(I) = SIG(I) + DSIG(I)
	1					 - DLAMD*DVECT(I)


100	CONTINUE


	CALL HTHIN(UNIAX,YIELD,WORKH) !!!

	
	CALL YTHIN(SIG ,CURY,AVECT)

	
	BRING = 1.0  !BRING IN THE EVERY LOOP
 	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 150 I = 1,3
150	SIG(I) = BRING*SIG(I)
	
	DO I = 1,3
	WORKH = WORKH + DLAMD*AVECT(I)*SIG(I)
	ENDDO

200	CONTINUE



	CALL YTHIN(SIG ,CURY,AVECT)


	CALL DTHIN(YOUNG,POISN,AVECT,DVECT)

	ABETA = 0.0
	DO I = 1,3
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DO 280 I = 1,3
	DO 280 J = 1,3
280	DEP(I,J) = DEP(I,J) !- DVECT(I)*DVECT(J)*ABETA

C	======================================================================
C	END OF ELASTIC-PLASTIC BEHAVIOR
C	======================================================================	
	
300	CONTINUE

	
	DO 350 I = 1,3
	EPSP(I) = EPS(I)
350	SIGP(I) = SIG(I)


	RETURN

	END

C	======================================================================	
C	======================================================================	
C	======================================================================	
	SUBROUTINE DTHIN(YOUNG,POISN,EP,SG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SG(3),EP(3),STIFG(3,3)
C	=======================================================================



	DO 110 I = 1,3
	DO 110 J = 1,3
110	STIFG(I,J) = 0.0
		 

	STIFG(1,1) = YOUNG/(2.0+2.0*POISN)
	STIFG(2,2) = YOUNG/(2.0+2.0*POISN)
	STIFG(3,3) = YOUNG


	DO 200 I = 1,3
	SG(I) = 0.0
	DO 200 J = 1,3
200	SG(I) = SG(I) + STIFG(I,J)*EP(J)

	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================	
	SUBROUTINE MTHIN(YOUNG,POISN,STIFG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION STIFG(3,3)
C	=======================================================================

	DO 110 I = 1,3
	DO 110 J = 1,3
110	STIFG(I,J) = 0.0
	 
	STIFG(1,1) = YOUNG/(2.0+2.0*POISN)
	STIFG(2,2) = YOUNG/(2.0+2.0*POISN)
	STIFG(3,3) = YOUNG

	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE YTHIN(SG,FN,AVECT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /THIN3/ GIIF,COH,PHI,ALPAB

	DIMENSION SG(3),AVECT(3)

	SL  = SG(1)
	SM  = SG(2)
	SN  = SG(3)


	VARI1 = SN
	VARJ2 = SL*SL + SM*SM + 
	1	    0.5*(SN - SN/3.0)*(SN - SN/3.0) +
	2		0.5*(-SN/3.0)*(-SN/3.0)         +
	3		0.5*(-SN/3.0)*(-SN/3.0)         



	AVECT  = 0.0


	FN  = (ALPAB*VARI1) + SQRT(VARJ2)
	CF1 = ALPAB 
	CF2 = 0.5/SQRT(VARJ2) 



C	=========================	
C	FLOW VECTOR
C	=========================
	AVECT(1) = CF2*2.0*SL
	AVECT(2) = CF2*2.0*SM
	AVECT(3) = (CF2*2.0*SN/3.0 + ALPAB)



	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE RTHIN(SG,DSG,YIELD,RFAC)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SG(3),DSG(3)
	DIMENSION SGR(3),AVECT(3)


	NITER = 30
	RFAC1 = 0.0
	RFAC2 = 1.0
	TOL = 0.00001
	RFACO = RFAC2
	DO ITER = 1,NITER
	
	RFAC = 0.5*(RFAC1 + RFAC2)
	DO I = 1,3
	SGR(I) = SG(I) + RFAC*DSG(I)
	ENDDO
	CALL YTHIN(SGR,FR,AVECT)


	TEST = FR - YIELD
	IF(TEST.LE.0.0) RFAC1 = RFAC
	IF(TEST.GT.0.0) RFAC2 = RFAC
	TEST = (RFAC-RFACO)/RFACO
	IF(ABS(TEST).LE.TOL) GO TO 200
	RFACO = RFAC
	ENDDO

100	CONTINUE

200	CONTINUE


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE STHIN(CURY,YIELD,NSTEP,ASTEP)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	=======================================================================

	ESCUR = CURY - YIELD
	NSTEP = ESCUR*8.0/YIELD + 1.0
	IF(NSTEP.GT.30) NSTEP = 30
	ASTEP = NSTEP

	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE HTHIN(UNIAX,YIELD,WORKH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /THIN3/ GIIF,COH,PHI,ALPAB
C	=======================================================================



	COHE =   COH*(1.0-WORKH/GIIF)
	IF(COHE.LT.0.0) COHE = 0.0


	YIELD1 = 6.0*COHE*COS(PHI)/ SQRT(3.0)/(3.0-SIN(PHI)) ! MOHR COULOMB EXTER
	ALPAB1 = 2.0*SIN(PHI)/ SQRT(3.0) /(3.0-SIN(PHI))     ! MOHR COULOMB EXTER

	YIELD2 = 6.0*COHE*COS(PHI)/ SQRT(3.0)/(3.0+SIN(PHI)) ! MOHR COULOMB INTER
	ALPAB2 = 2.0*SIN(PHI)/ SQRT(3.0) /(3.0+SIN(PHI))     ! MOHR COULOMB INTER

	YIELD = YIELD1 !0.5*(YIELD1 + YIELD2) !
	ALPAB = ALPAB1 !0.5*(ALPAB1 + ALPAB2) !



	RETURN

	END	

C	======================================================================	

