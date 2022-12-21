C	===================================================================
C	MATERIAL MODEL USED IN ZERO THICKNESS CONTACT INTERFACE ELEMENT
C	===================================================================
	SUBROUTINE BOND3D(SIGP,EPSP,WORKH,SIG,EPS,DEP,PROPM,MTMOD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	===================================================================
	COMMON /BOND3/ GFII,GFI,RCOH,RPHI,ALPAB,MTMD

	DIMENSION SIGP(3),EPSP(3),DEP(3,3)
	DIMENSION SIG(3),EPS(3)
	DIMENSION AVECT(3),DVECT(3)
	DIMENSION DSIG(3),DEPS(3),DEPCR(3)
	DIMENSION PROPM(*)

	MTMD = MTMOD

	YOUNG = PROPM(7)
	YOUNS = PROPM(8)
	UNIAX = PROPM(9)
	RCOH  = PROPM(10)
	ANG   = PROPM(11)

	PI    = 3.141592654 
	RPHI  = ANG*PI/180.0
	

	CALL HBOND(WORKH,UNIAX,YIELD) !!!


	DO 10 I = 1,3
10	DEPS(I) = EPS(I) - EPSP(I)


	CALL DBOND(YOUNG,YOUNS,DEPS,DSIG,WORKH)

	CALL MOBON(YOUNG,YOUNS,DEP,WORKH)

	AFIRT   = 0.0
	DO 20 I = 1,3
	SIG(I)  = SIGP(I) + DSIG(I)
20	AFIRT   = AFIRT + SIG(I)

	IF(AFIRT.EQ.0.0.OR.MTMD.EQ.1) GO TO 300 !FIRST ENTRY STEP OR ELASTIC MATERIAL

	ASCND   = 0.0
	DO I = 1,3
	ASCND   = ASCND + SIGP(I)
	ENDDO

	IF(ASCND.EQ.0.0) THEN !SECOND ENTRY STEP
	PREY = 0.0
	ELSE
	CALL YBOND(SIGP,PREY,AVECT,WORKH)
	ENDIF

	CALL YBOND(SIG ,CURY,AVECT,WORKH)



	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP

C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,3
30	SIG(I) = SIGP(I)


	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RBOND(SIG,DSIG,YIELD,RFAC,WORKH)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL SBOND(CURY,YIELD,NSTEP,ASTEP)


	DO 40 I = 1,3
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,3
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP

	DEPCR = 0.0 
	DO 200 ISTEP = 1,NSTEP

	
	CALL YBOND(SIG ,CURY,AVECT,WORKH)

	CALL DBOND(YOUNG,YOUNS,AVECT,DVECT,WORKH)

	ABETA = 0.0
	ADUMM = 0.0
	DO I = 1,3
	ADUMM = ADUMM + AVECT(I)*DSIG(I)
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA)

	DLAMD = ADUMM*ABETA
	IF(DLAMD.LT.0.0) DLAMD = 0.0


	DO 100 I = 1,3  !BACK TO SURFACE
	SIG(I) = SIG(I) + DSIG(I)
	1					 - DLAMD*DVECT(I)

100	CONTINUE

	SL  = SIG(1)
	SM  = SIG(2)
	SN  = SIG(3)

	DEPCR(1) = DEPCR(1) + DLAMD*AVECT(1)
	DEPCR(2) = DEPCR(2) + DLAMD*AVECT(2)
	DEPCR(3) = DEPCR(3) + DLAMD*AVECT(3)
	



	CALL HBOND(WORKH,UNIAX,YIELD) !!!


C200	CONTINUE



	CALL YBOND(SIG ,CURY,AVECT,WORKH)

	

	BRING = 1.0  !BRING IN THE EVERY LOOP
 	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 150 I = 1,3
150	SIG(I) = BRING*SIG(I)



200	CONTINUE


	DO I = 1,3  ! SLIP ENERGY

	TEST = SIG(I)*DEPCR(I)
	IF(TEST.GE.0.0)THEN

	IF(SN.GT.0.0) THEN
	WORKH = WORKH + SIG(I)*DEPCR(I)                   
	ELSE
	IF(I.LE.2) THEN
	WORKH = WORKH + SIG(I)*DEPCR(I)*
	1		        (1.0+SN*TAN(RPHI)/SQRT(SL*SL + SM*SM))
	ENDIF
	ENDIF

	ENDIF

	ENDDO
C	WORKH = 0.0	


	CALL YBOND(SIG ,CURY,AVECT,WORKH)


	CALL DBOND(YOUNG,YOUNS,AVECT,DVECT,WORKH)

	ABETA = 0.0
	DO I = 1,3
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA)

	DO 280 I = 1,3
	DO 280 J = 1,3
280	DEP(I,J) = DEP(I,J) - DVECT(I)*DVECT(J)*ABETA

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
	SUBROUTINE DBOND(YOUNG,YOUNS,EP,SG,WORKH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SG(3),EP(3),STIFG(3,3)

	COMMON /BOND3/ GFII,GFI,RCOH,RPHI,ALPAB,MTMD
C	=======================================================================



	DO 110 I = 1,3
	DO 110 J = 1,3
110	STIFG(I,J) = 0.0


	STIFG(1,1) = YOUNS
	STIFG(2,2) = YOUNS
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
	SUBROUTINE MOBON(YOUNG,YOUNS,STIFG,WORKH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION STIFG(3,3)

	COMMON /BOND3/ GFII,GFI,RCOH,RPHI,ALPAB,MTMD
C	=======================================================================


	DO 110 I = 1,3
	DO 110 J = 1,3
110	STIFG(I,J) = 0.0


		 
	STIFG(1,1) = YOUNS
	STIFG(2,2) = YOUNS
	STIFG(3,3) = YOUNG

	

	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE YBOND(SG,FN,AVECT,WORKH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /BOND3/ GFII,GFI,RCOH,RPHI,ALPAB,MTMD

	DIMENSION SG(3),AVECT(3)

	SL  = SG(1)
	SM  = SG(2)
	SN  = SG(3)


	VARI1 = SN
	VARJ2 = SL*SL + SM*SM + 
	1	    0.5*(SN - SN/3.0)*(SN - SN/3.0) +
	2		0.5*(-SN/3.0)*(-SN/3.0)         +
	3		0.5*(-SN/3.0)*(-SN/3.0)    
	
C     MODIFY DRAGER-PRAGER -- NO CONTRIBUTION OF NORMAL STRESS IN THE SHEAR TERMS	
	VARJ2 = SL*SL + SM*SM      

	IF(WORKH.GE.GFII) WORKH = GFII
	COHE =   RCOH !*(1.0-(WORKH/GFII))

	AVECT  = 0.0


	IF(SN.GT.0.0) THEN
	CFD = 1.0
	ELSE
	CFD = 1.0 !COHE/RCOH
	ENDIF


	IF(MTMD.EQ.3) THEN
	FN  = (ALPAB*VARI1) + SQRT(VARJ2)
	CF1 = ALPAB 
	CF2 = 0.5/SQRT(VARJ2) 
	ELSEIF(MTMD.EQ.4) THEN
	TAN2 = TAN(RPHI)*TAN(RPHI)
	FN = SL*SL + SM*SM  + 2.0*COHE*TAN(RPHI)*SN - TAN2*SN*SN
	ENDIF

C	=========================	
C	FLOW VECTOR
C	=========================
	IF(MTMD.EQ.3) THEN
	AVECT(1) = CF2*2.0*SL
	AVECT(2) = CF2*2.0*SM
	AVECT(3) = (CF2*2.0*SN/3.0 + ALPAB) * CFD
	
	AVECT(3) = (ALPAB) * CFD !MODIFIED DRAGER PRAGER   NO CONTRIBUTION OF NORMAL STRESS ON FLOW VECTOR
	ELSEIF(MTMD.EQ.4) THEN
	AVECT(1) = 2.0*SL
	AVECT(2) = 2.0*SM
	AVECT(3) = 2.0*TAN(RPHI)*(COHE-SN*TAN(RPHI))
	
	AVECT(3) = 0.0 !MODIFIED HYPERBOLIC   NO CONTRIBUTION OF NORMAL STRESS ON FLOW VECTOR
	ENDIF


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE RBOND(SG,DSG,YIELD,RFAC,WORKH)
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
	CALL YBOND(SGR,FR,AVECT,WORKH)


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
	SUBROUTINE SBOND(CURY,YIELD,NSTEP,ASTEP)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	=======================================================================

	ESCUR = CURY - YIELD
	ASTEP = ESCUR*8.0/YIELD + 1.0
	NSTEP = ASTEP
	IF(NSTEP.GT.30) NSTEP = 30
	ASTEP = NSTEP


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE HBOND(WORKH,UNIAX,YIELD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	COMMON /BOND3/ GFII,GFI,RCOH,RPHI,ALPAB,MTMD
C	=======================================================================

C	IF(WORKH.GE.GFII) WORKH = GFII
C	IF(WORKH.GE.GFI)  WORKH = GFI

C	COHE =   RCOH *(1.0-(WORKH/GFII))
C	UNIE =  UNIAX *(1.0-(WORKH/GFI ))
C	IF(COHE.LE.0.0D0) COHE = 0.0D0
	COHE =   RCOH
	UNIE =  UNIAX

	IF(MTMD.EQ.3) THEN

	YIELD1 = 6.0*COHE*COS(RPHI)/ SQRT(3.0)/(3.0-SIN(RPHI)) ! MOHR COULOMB EXTER
	ALPAB1 = 2.0*SIN(RPHI)/ SQRT(3.0) /(3.0-SIN(RPHI))     ! MOHR COULOMB EXTER

	YIELD2 = 6.0*COHE*COS(RPHI)/ SQRT(3.0)/(3.0+SIN(RPHI)) ! MOHR COULOMB INTER
	ALPAB2 = 2.0*SIN(RPHI)/ SQRT(3.0) /(3.0+SIN(RPHI))     ! MOHR COULOMB INTER

	YIELD = 0.5*(YIELD1 + YIELD2) !
	ALPAB = 0.5*(ALPAB1 + ALPAB2) !

C     MODIFIIED DRAGER PRAGER	
	YIELD = COHE
	ALPAB = TAN(RPHI)

	ELSEIF(MTMD.EQ.4) THEN


	YIELD = COHE*COHE - (COHE-UNIE*TAN(RPHI))*(COHE-UNIE*TAN(RPHI))
	ALPAB = 0.0D0


	ENDIF



	RETURN

	END	

C	======================================================================	
C	===================================================================
C	MATERIAL MODEL USED IN ZERO THICKNESS CONTACT INTERFACE ELEMENT 2D
C	===================================================================
	SUBROUTINE BOND2D(SGP,ESP,WORKH,SGG,EPP,DE,PROPM,MTMOD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	===================================================================
	COMMON /BOND3/ GFII,GFI,RCOH,RPHI,ALPAB,MTMD

	DIMENSION SGP(2),ESP(2),DE(2,2),SGG(2),EPP(2)
	DIMENSION SIGP(3),EPSP(3),DEP(3,3)
	DIMENSION SIG(3),EPS(3)
	DIMENSION AVECT(3),DVECT(3)
	DIMENSION DSIG(3),DEPS(3),DEPCR(3)
	DIMENSION PROPM(*)

C	WRITE(*,*) WORKH,GFII
C	PAUSE
	MTMD = MTMOD

	YOUNG = PROPM(7)
	YOUNS = PROPM(8)
	UNIAX = PROPM(9)
	RCOH  = PROPM(10)
	ANG   = PROPM(11)

	PI    = 3.141592654 
	RPHI  = ANG*PI/180.0


	SIGP(1) = SGP(1)
	SIGP(2) = 0.0D0
	SIGP(3) = SGP(2)
	EPSP(1) = ESP(1)
	EPSP(2) = 0.0D0
	EPSP(3) = ESP(2)
	EPS(1)  = EPP(1)
	EPS(2)  = 0.0D0
	EPS(3)  = EPP(2)


	CALL HBOND(WORKH,UNIAX,YIELD) !!!

C	WRITE(*,*) WORKH,YIELD


	DO 10 I = 1,3
10	DEPS(I) = EPS(I) - EPSP(I)



	CALL DBOND(YOUNG,YOUNS,DEPS,DSIG,WORKH)

	CALL MOBON(YOUNG,YOUNS,DEP,WORKH)

	AFIRT   = 0.0
	DO 20 I = 1,3
	SIG(I)  = SIGP(I) + DSIG(I)
20	AFIRT   = AFIRT + SIG(I)

	IF(AFIRT.EQ.0.0.OR.MTMD.EQ.1) GO TO 300 !FIRST ENTRY STEP OR ELASTIC MATERIAL

	ASCND   = 0.0
	DO I = 1,3
	ASCND   = ASCND + SIGP(I)
	ENDDO

	IF(ASCND.EQ.0.0) THEN !SECOND ENTRY STEP
	PREY = 0.0
	ELSE
	CALL YBOND(SIGP,PREY,AVECT,WORKH)
	ENDIF

	CALL YBOND(SIG ,CURY,AVECT,WORKH)



	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP



C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,3
30	SIG(I) = SIGP(I)



	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL RBOND(SIG,DSIG,YIELD,RFAC,WORKH)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF


	CALL SBOND(CURY,YIELD,NSTEP,ASTEP)

				


	DO 40 I = 1,3
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,3
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP

	DEPCR = 0.0 
	DO 200 ISTEP = 1,NSTEP

	
	CALL YBOND(SIG ,CURY,AVECT,WORKH)

	CALL DBOND(YOUNG,YOUNS,AVECT,DVECT,WORKH)

	ABETA = 0.0
	ADUMM = 0.0
	DO I = 1,3
	ADUMM = ADUMM + AVECT(I)*DSIG(I)
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA)

	DLAMD = ADUMM*ABETA
	IF(DLAMD.LT.0.0) DLAMD = 0.0


	DO 100 I = 1,3  !BACK TO SURFACE
	SIG(I) = SIG(I) + DSIG(I)
	1					 - DLAMD*DVECT(I)

100	CONTINUE

	SL  = SIG(1)
	SM  = SIG(2)
	SN  = SIG(3)

	DEPCR(1) = DEPCR(1) + DLAMD*AVECT(1)
	DEPCR(2) = DEPCR(2) + DLAMD*AVECT(2)
	DEPCR(3) = DEPCR(3) + DLAMD*AVECT(3)
	



	CALL HBOND(WORKH,UNIAX,YIELD) !!!


C200	CONTINUE



	CALL YBOND(SIG ,CURY,AVECT,WORKH)

	

	BRING = 1.0  !BRING IN THE EVERY LOOP
 	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 150 I = 1,3
150	SIG(I) = BRING*SIG(I)



200	CONTINUE


	DO I = 1,3  ! SLIP ENERGY

	TEST = SIG(I)*DEPCR(I)
	IF(TEST.GE.0.0)THEN

	IF(SN.GT.0.0) THEN
	WORKH = WORKH + SIG(I)*DEPCR(I)                   
	ELSE
	IF(I.LE.2) THEN
	WORKH = WORKH + SIG(I)*DEPCR(I)*
	1		        (1.0+SN*TAN(RPHI)/SQRT(SL*SL + SM*SM))
	ENDIF
	ENDIF

	ENDIF

	ENDDO
C	WORKH = 0.0	



	CALL YBOND(SIG ,CURY,AVECT,WORKH)


	CALL DBOND(YOUNG,YOUNS,AVECT,DVECT,WORKH)

	ABETA = 0.0
	DO I = 1,3
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA)

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


	SGP(1) = SIGP(1)
	SGP(3) = SIGP(3)
	ESP(1) = EPSP(1)
	ESP(2) = EPSP(3)

	SGG(1) = SIG(1)
	SGG(2) = SIG(3)

	DE(1,1) = DEP(1,1)
	DE(1,2) = DEP(1,3)
	DE(2,1) = DEP(3,1)
	DE(2,2) = DEP(3,3)

	

	RETURN

	END

C	======================================================================	
C	======================================================================	
C	======================================================================
