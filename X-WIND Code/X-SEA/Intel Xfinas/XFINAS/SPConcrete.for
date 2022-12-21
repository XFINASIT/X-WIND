C	=============================================================================
	SUBROUTINE SPCONT(EPS,IPT,DR,SIGR,NMP,PROPM,SPCEPS,FSTRN,
	1				  NLOPT,NFIB,ISEC,IEG,NPT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	=============================================================================
C	=============================================================================
	DIMENSION EPS(7),DR(49),SIGR(7)
	DIMENSION DEPS(7),WAG(4),PROPM(NMP,1)
	DIMENSION SPCEPS(NPT,8),FSTRN(3*NPT,1)
	DIMENSION SDATA(20)
	DIMENSION AA(10)
C	=============================================================================

	III = 1
	INAME(1:4) = [5,1,ISEC,IEG] !XSEH
	CALL ICONC(INAME,NAMEI)
	DO II = 1,7
	CALL MRELFIL(NAMEI,AA(II),III,II,0) !Store Section Props 
	ENDDO
	NSTS  = INT(AA(2))
	NTOUR = INT(AA(6))
	NSEGT = INT(AA(7))

C	================================================
C	COMPUTE THE GENERALIZED STRAIN
C	================================================
	DO I = 1,7
	EPSP = SPCEPS(IPT,I)
	DEPS(I) = EPS(I) - EPSP
	SPCEPS(IPT,I) = EPS(I)
	ENDDO
C	================================================


C	INITIALIZED THE ARRAY
	DR = 0.0
	SIGR = 0.0


C	AREA = 0.0
C	================================================
	DO 1000 IFIB = 1,NFIB  !LOOP OVER ALL FIBER
C	================================================

	NSTAT = 0

	III = 1 + NSEGT + (1+NTOUR+NFIB+NSTS)*(IPT-1) + 1 + NTOUR + IFIB
	DO II = 1,7
	CALL MRELFIL(NAMEI,AA(II),III,II,0) !Store Section Props 
	ENDDO
	AFIBR = AA(1) ; SCOR = AA(2) ; TCOR = AA(3) ; FMT = AA(4)
	FICS  = AA(5) ; WP   = AA(6) ; TR = AA(7)

	MTSET = INT(FMT)
	IDTCS = INT(FICS)

C	===========================
C	CALL THE MATERIAL PARAMETER
C	===========================
	SELECTCASE (IDTCS)
	CASE(0)
	FC   = PROPM(14,MTSET)
	FCR  = PROPM(15,MTSET)
	PRC  = PROPM(2,MTSET)
	YMC  = PROPM(1,MTSET)	
	CASE(1)
	YMS  = PROPM(1,MTSET)
	YMSH = PROPM(4,MTSET)
	PRS  = PROPM(2,MTSET)
	YLDS = PROPM(3,MTSET)
	ENDSELECT

	EUC = 0.0038
	EO  = 2.0*FC/YMC


C	================================================
C	COMPUTE FIBER STRAIN
C	================================================
	FEPS  = EPS(1)  + TCOR*EPS(5)  - SCOR*EPS(6)  + WP*EPS(7)
	DFEPS = DEPS(1) + TCOR*DEPS(5) - SCOR*DEPS(6) + WP*DEPS(7)
C	================================================

C	================================================
C	NEXT BLOCK FOR STEEL BAR (IDTCS = 1)
C	================================================
	IF(IDTCS.EQ.1) GO TO 500 ! STEEL
C	================================================

	
C	AREA = AREA + AFIBR
C	=========================
C	START CONCRETE MODULE
C	=========================
	STSPS = FSTRN(1+3*(IPT-1),IFIB)
	YIELD = FSTRN(2+3*(IPT-1),IFIB)

	
	EM = FEPS
	IF(EM.LE.0.0) MSTAT = 1
	IF(EM.GT.0.0) MSTAT = 2


	IF(MSTAT.EQ.1) THEN !COMPRESSION
	
	DSIG  = YMC*DFEPS
	STSCR = STSPS + DSIG


	IF(ABS(STSCR).GE.YIELD) THEN  !LOADING <---MATERIAL STAGE

	NSTAT = FSTRN(3+3*(IPT-1),IFIB)


	IF(NLOPT.EQ.0) THEN
	YMCC = YMC
	YIELD =  ABS(STSCR)
	STSPS = STSCR
	STSCR = STSCR
	GOTO 111
	ENDIF

		
	EM = ABS(EM)
	IF(EM.LE.EO) THEN
	CON1  = EM/EO
	STSCR = FC*CON1*(2.0-CON1)
	YMCC  = YMC*(1.0-CON1)
	ENDIF
		
	IF(EM.GT.EO.AND.EM.LT.EUC) THEN
	STSCR = FC  - (0.15*FC*(EM - EO)/(EUC - EO))
	YMCC  = 0.0
	ENDIF

	IF(EM.GE.EUC.OR.NSTAT.EQ.7) THEN
	STSCR = 0.3*FC
	EUCM = 10.0*EUC
	STSCR = (0.85*FC) - (0.85*FC*(EM - EUC)/(EUCM - EUC))
	IF(EM.GE.EUCM) STSCR = 0.0
	YMCC  = 0.0
	FSTRN(3+3*(IPT-1),IFIB) = 7
	ENDIF

C	CONTROL STABILITY OF THE SECTION 
C	IF(YMCC.LE.0.5*YMC) YMCC = 0.5*YMC

	YIELD =  ABS(STSCR)

	STSPS = -STSCR
	STSCR = -STSCR

111	CONTINUE

	ELSEIF(ABS(STSCR).LT.YIELD) THEN ! ELASTIC (RELOAD & UNLOAD)

	STSCR = STSCR
	STSPS = STSCR
	YMCC  = YMC


	IF(STSCR.GT.0.0) THEN
	NSTAT = FSTRN(3+3*(IPT-1),IFIB)

	IF(STSCR.GT.FCR) THEN !CRACK
	STSCR = 0.0
	YMCC  = 0.0 !0.5*YMC
	FSTRN(3+3*(IPT-1),IFIB) = 3
	ENDIF

	IF(NSTAT.EQ.3) THEN
	STSCR = 0.0
	YMCC  = 0.0 !0.5*YMC
	ENDIF

	IF(NSTAT.EQ.7) THEN
	STSCR = 0.0
	YMCC  = 0.0
	FSTRN(3+3*(IPT-1),IFIB) = 7
	ENDIF

	ENDIF


	ENDIF !END MATERIAL STAGE
	

	ELSEIF(MSTAT.EQ.2) THEN !TENSION

	NSTAT = FSTRN(3+3*(IPT-1),IFIB)

	ETU = FCR/YMC
	YMCC = YMC
	DSIG = YMC*DFEPS
	STSCR = STSPS + DSIG

	IF(NLOPT.EQ.0) THEN
	YMCC = YMC
	GOTO 222
	ENDIF


	IF(STSCR.GT.FCR) THEN !CRACK
	STSCR = 0.0
	YMCC  = 0.0 !0.5*YMC
	FSTRN(3+3*(IPT-1),IFIB) = 3
	ENDIF

	IF(NSTAT.EQ.3) THEN
	STSCR = 0.0
	YMCC  = 0.0 !0.5*YMC !0.0
	ENDIF

	IF(FEPS.GT.ETU) THEN !TENSION STIFFENING
	STSCR = FCR*((ETU/FEPS)**0.9)
C	IF(FEPS.GT.5.0*ETU) STSCR = 0.0
	YMCC  = 0.0 !0.5*YMC
	FSTRN(3+3*(IPT-1),IFIB) = 3

	ENDIF
	

	IF(NSTAT.EQ.7) THEN
	STSCR = 0.0
	YMCC  = 0.0
	FSTRN(3+3*(IPT-1),IFIB) = 7
	ENDIF

222	CONTINUE

	STSPS = STSCR

	ENDIF !END TENSION AND COMPRESSION

	YMG = 0.5*YMC/(1.0 + PRC)

C	CONTROL STABILITY OF THE SECTION 
CC	IF(YMCC.LE.0.5*YMC) YMCC = 0.5*YMC
C	YMCC = YMC


	GO TO 700

C	=========================
C	END CONCRETE MODULE
C	=========================

500	CONTINUE

C	=========================
C	START STEEL MODULE
C	=========================	

	
	STSPS = FSTRN(1+3*(IPT-1),IFIB)
	EPSTN = FSTRN(2+3*(IPT-1),IFIB)

	CALL SPCSTL(STSPS,DFEPS,STSCR,YMS,YMSH,EPSTN,YLDS,
	1			YMCC,NLOPT) 
	

	STSPS = STSCR
	YIELD = EPSTN

	YMG = 0.5*YMS/(1.0 + PRS)


C	CONTROL STABILITY OF THE SECTION
CC	IF(YMCC.LE.0.5*YMS) YMCC = 0.5*YMS
C	YMCC = YMS

C	=========================
C	END STEEL MODULE
C	=========================
	
	
700	CONTINUE

C	============================
C	STORE THE WORKING ARRAY
C	============================
	FSTRN(1+3*(IPT-1),IFIB) = STSPS
	FSTRN(2+3*(IPT-1),IFIB) = YIELD

C	============================
C	INTEGRATING OVER THE SECTION
C	============================
C	RIGIDITY MATRIX 
	DR(1)  = DR(1)  + YMCC*AFIBR										!EA
	DR(25) = DR(25) + YMG*(SCOR*SCOR + TCOR*TCOR)*AFIBR - YMG*TR*AFIBR	!GJr

	DR(29) = DR(29) + TCOR*YMCC*AFIBR									!EQs
	DR(33) = DR(33) + TCOR*TCOR*YMCC*AFIBR								!EIs
	
	DR(36) = DR(36) - SCOR*YMCC*AFIBR									!-EQt
	DR(40) = DR(40) - SCOR*TCOR*YMCC*AFIBR								!-EIst
	DR(41) = DR(41) + SCOR*SCOR*YMCC*AFIBR								!EIt

	DR(43) = DR(43) + WP*YMCC*AFIBR										!EQw
	DR(47) = DR(47) + WP*TCOR*YMCC*AFIBR								!EItw
	DR(48) = DR(48) - WP*SCOR*YMCC*AFIBR								!-EIsw
	DR(49) = DR(49) + WP*WP*YMCC*AFIBR									!EIw


C	STRESS CONTRIBUTION
	SIGR(1) = SIGR(1) + STSCR*AFIBR				!P
	SIGR(5) = SIGR(5) + TCOR*STSCR*AFIBR		!Ms
	SIGR(6) = SIGR(6) - SCOR*STSCR*AFIBR		!Mt
	SIGR(7) = SIGR(7) + WP*STSCR*AFIBR			!Mw
	SIGR(4) = SIGR(4) + YMG*(SCOR*SCOR + TCOR*TCOR - TR)*EPS(4)*AFIBR !TORSIONAL MOMENT

C	================================================
1000	CONTINUE  !END OF LOOP OVER ALL FIBER
C	================================================


	RETURN

	END

C	=============================================================================
C	=============================================================================
C	=============================================================================
C	=============================================================================
	SUBROUTINE SPCSTL(SIGP,DEPS,SIG,YOUNG,YOUN2,EPSTN,YIELD,
	1				  YMODU,NLOPT) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================

	HARDS = YOUNG*YOUN2/(YOUNG-YOUN2)
	DSIG  = DEPS*YOUNG

	SIGO = SIGP

	SIG  = SIGO + DSIG

	PREY  = YIELD + HARDS*EPSTN

	IF(NLOPT.EQ.0) GOTO 40

	IF(ABS(SIGO).GE.PREY) GO TO 20
	ESCUR = ABS(SIG) - PREY
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(DSIG)
	GO TO 30

20	IF(SIG.GT.0.0.AND.DSIG.LT.0.0) GO TO 40
	IF(SIG.LT.0.0.AND.DSIG.GT.0.0) GO TO 40
	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	SIGO  = SIGO + REDUC*DSIG + 
	1	    RFACT*DSIG*(1.0-YOUNG/(YOUNG+HARDS))
	EPSTN = EPSTN + RFACT*DSIG/(YOUNG+HARDS)
	EPSTN = ABS(EPSTN)
	YMODU = YOUN2

	GO TO 50	

40	YMODU = YOUNG
	SIGO  = SIGO + DSIG
50	SIGP  = SIGO

	SIG  = SIGP


	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================

	