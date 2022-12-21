C	===================================================================
C	===================================================================
C	======================================================================

	SUBROUTINE MOHR2D(SIGPM,EPSPM,EPSTN,SIGM,EPSM,DEP,PROPM,
	1				  STSZ,STNZ,IFLAG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /SOILP/ COH,PHI

	DIMENSION SIGPM(4),EPSPM(4)
	DIMENSION SIGM(4),EPSM(4)
	DIMENSION SIGP(6),EPSP(6),DEP(6,6)
	DIMENSION SIG(6),EPS(6)
	DIMENSION AVECT(6),DVECT(6)
	DIMENSION DSIG(6),DEPS(6)
	DIMENSION PROPM(1),DEPM(6,6)
	
	PI  = 3.141592654 

	YOUNG = PROPM(1)
	POISN = PROPM(2)

	COH  = PROPM(10)
	ANG  = PROPM(11)
	PHI  = ANG*PI/180.0


	HARDP = 0.0D0

	CALL HARPMOHR(YIELD) !!!


	SIGP(1) = SIGPM(1)
	SIGP(2) = SIGPM(2)
	SIGP(3) = SIGPM(4) !STSZ
	SIGP(4) = SIGPM(3)
	SIGP(5) = 0.0
	SIGP(6) = 0.0

	EPSP(1) = EPSPM(1)
	EPSP(2) = EPSPM(2)
	EPSP(3) = EPSPM(4) !STNZ       !0.0
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
	CALL YSURMOR(SIGP,PREY,AVECT)
	ENDIF

	CALL YSURMOR(SIG ,CURY,AVECT)


	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP


C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,6
30	SIG(I) = SIGP(I)

	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RFACTMOR(SIG,DSIG,YIELD,RFAC)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL STNUM(CURY,YIELD,NSTEP,ASTEP)

	DO 40 I = 1,6
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,6
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP


	DO 200 ISTEP = 1,NSTEP
	
	CALL YSURMOR(SIG,CURY,AVECT)

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

	CALL HARPMOHR(YIELD) !!!
	

200	CONTINUE

	
	CALL YSURMOR(SIG,CURY,AVECT)



	BRING = 1.0
	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 250 I = 1,6
250	SIG(I) = BRING*SIG(I)


	CALL YSURMOR(SIG,CURY,AVECT)


	CALL DMATX(YOUNG,POISN,AVECT,DVECT,6)

	ABETA = 0.0
	DO I = 1,6
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO


	IF(ABETA.NE.0.0D0) ABETA = 1.0 / (ABETA + HARDP)

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
	SUBROUTINE MOHR3D(SIGP,EPSP,EPSTN,SIG,EPS,DEP,PROPM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	===================================================================
	COMMON /SOILP/ COH,PHI


	DIMENSION SIGP(6),EPSP(6),DEP(6,6)
	DIMENSION SIG(6),EPS(6)
	DIMENSION AVECT(6),DVECT(6)
	DIMENSION DSIG(6),DEPS(6)
	DIMENSION PROPM(1),DEPM(6,6)
	
	PI  = 3.141592654 

	YOUNG = PROPM(1)
	POISN = PROPM(2)

	COH  = PROPM(10)
	ANG  = PROPM(11)
	PHI  = ANG*PI/180.0
	
	HARDP = 0.0D0

	CALL HARPMOHR(YIELD) !!!

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
	CALL YSURMOR(SIGP,PREY,AVECT)
	ENDIF

	CALL YSURMOR(SIG ,CURY,AVECT)


	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP

	


C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,6
30	SIG(I) = SIGP(I)

	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RFACTMOR(SIG,DSIG,YIELD,RFAC)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL STNUM(CURY,YIELD,NSTEP,ASTEP)

	DO 40 I = 1,6
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,6
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP


	DO 200 ISTEP = 1,NSTEP
	
	CALL YSURMOR(SIG,CURY,AVECT)


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

	CALL HARPMOHR(YIELD) !!!


200	CONTINUE

	CALL YSURMOR(SIG,CURY,AVECT)

	BRING = 1.0
	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 250 I = 1,6
250	SIG(I) = BRING*SIG(I)


	CALL YSURMOR(SIG,CURY,AVECT)

	CALL DMATX(YOUNG,POISN,AVECT,DVECT,6)

	ABETA = 0.0
	DO I = 1,6
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	IF(ABETA.NE.0.0D0) ABETA = 1.0 / (ABETA + HARDP)

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
	SUBROUTINE YSURMOR(SIG,FUNC,AVECT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)


	COMMON /SOILP/ COH,PHI


	DIMENSION SIG(6),AVECT(6),PST(3)
	DIMENSION AV1(6),AV2(6),AV3(6)
	

	CALL INVAT(SIG,PST,VI1,VI2,VI3,VJ2,VJ3,THETA)


	SX  = SIG(1)
	SY  = SIG(2)
	SZ  = SIG(3)
	TXY = SIG(4)
	TXZ = SIG(5)
	TYZ = SIG(6)

	SM  = (SX+SY+SZ)/3.0D0

	SXP = SX - SM
	SYP = SY - SM
	SZP = SZ - SM
	TXYP= TXY
	TXZP= TXZ
	TYZP= TYZ

C	======================================================================	
C	PLASTIC POTENTIAL DERIVATIVE MATRIX
C	======================================================================	
	AV1(1:6) = [1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0]

	IF(VI1.EQ.0.0.AND.VJ2.EQ.0.0) THEN
	FUNC = 0.0D0
	AVECT(1:6) = [0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0]
	RETURN
	ENDIF

	AV2(1:3) = [SXP,SYP,SZP]
	AV2(4) = 2.0*TXYP
	AV2(5) = 2.0*TXZP
	AV2(6) = 2.0*TYZP
	AV2 = 0.5*AV2/SQRT(VJ2)

	D1 = SYP*SZP - TYZ*TYZ + VJ2/3.0
	D2 = SXP*SZP - TXZ*TXZ + VJ2/3.0
	D3 = SXP*SYP - TXY*TXY + VJ2/3.0
	D4 = 2.0D0*(TXZ*TXY - SXP*TYZ)
	D5 = 2.0D0*(TXZ*TYZ - SYP*TXZ)
	D6 = 2.0D0*(TYZ*TXZ - SZP*TXY)

	AV3(1:6) = [D1,D2,D3,D4,D5,D6]
C	======================================================================	
C	END OF PLASTIC POTENTIAL DERIVATIVE MATRIX
C	======================================================================


C	=====================================	
C	YIELD FUNCTION FLOW VECTOR COEFICIENT
C	=====================================

C	=========================
C	MOHR-COULOMB
C	=========================
	FUNC  = VI1*SIN(PHI)/3.0D0 + 
	1		SQRT(VJ2)*(COS(THETA)-SIN(THETA)*SIN(PHI)/SQRT(3.0))

C	WRITE(*,*) FUNC,VI1*SIN(PHI)/3.0D0,
C	1	SQRT(VJ2)*(COS(THETA)-SIN(THETA)*SIN(PHI)/SQRT(3.0))
C	PAUSE

	CF1 = SIN(PHI)/3.0
	ABTHE = ABS(THETA*57.29577951308)
	IF(ABTHE.LT.29.0) GOTO 30
	CF3 = 0.0D0
	PLUMI = 1.0D0
	IF(THETA.GT.0.0D0) PLUMI = -1.0D0
	CF2 = 0.5*( SQRT(3.0) + PLUMI*CF1*SQRT(3.0) )
	GOTO 40 
30	CF2 = COS(THETA)*( 1.0+TAN(THETA)*TAN(3.0*THETA) + 
	1	               SIN(PHI)*(TAN(3.0*THETA)-TAN(THETA))/
     2				   SQRT(3.0) )
	CF3 = (SQRT(3.0)*SIN(THETA) + COS(THETA)*SIN(PHI)) /
	1	  (2.0*VJ2*COS(3.0*THETA))
	GOTO 40
C	=========================
C	END MOHR-COULOMB
C	=========================


C	============================================	
C	END OF YIELD FUNCTION FLOW VECTOR COEFICIENT
C	============================================


C	=========================	
C	FLOW VECTOR
C	=========================
40	DO I = 1,6
	AVECT(I) = CF1*AV1(I) + CF2*AV2(I) + CF3*AV3(I)
	ENDDO


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE RFACTMOR(SG,DSG,YIELD,RFAC)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

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
	CALL YSURMOR(SGR,FR,AVECT)
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
	SUBROUTINE HARPMOHR(YIELD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	COMMON /SOILP/ COH,PHI

	YIELD = COH*COS(PHI) ! Tresca, Mohr Coulomb

	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
