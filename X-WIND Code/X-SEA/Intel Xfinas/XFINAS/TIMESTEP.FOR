C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================
	SUBROUTINE TIMESTEP

C	==================================================================
C	PURPOSE:	TO PERFORM NONLINEAR TIMESTEP ANALYSIS
C	==================================================================

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON A(9000000),IA(9000000)

	DIMENSION LMEQ1(10),LMEQ2(10)
C	-------------------------------------------------------------------
C	-------------------------------------------------------------------
	DIMENSION AA(NWK),DREF(NEQ)
C	-------------------------------------------------------------------

	INDPD = KPOSD

	ITOPT = 3
	NTIME = 400
	NITER = 20
	DTIM  = 0.001
	TOLD  = 0.001

	LMEQ1(1:2) = [2,1]
	LMEQ2(1:2) = [5,5]
	LMQ1 = LMEQ1(1)
	LMQ2 = LMEQ2(1)

C	TEST WHETHER ROTATION DOF HAS BEEN INVOLVE OR NOT
	ITESR = 0
	DO I = 1,9
	IF(IDOF(I).GE.4.AND.IDOF(I).LE.6) ITESR = 1
	ENDDO

C     --------------
C     INITIALISATION
C     --------------
	RHO   = 0.0D0
	CTIM  = 0.0D0
	ISTEP = 0
	CALL CLEARA (A(LDT),NEQ)

C	--------------------------------
	DO 5000 ITIME = 1,NTIME
C	--------------------------------
	CTIM = CTIM + DTIM

	IF(ITOPT.EQ.1) THEN    !FOR LOAD CONTROL
      DRHO = A(LLF+ITIME-1)  !RFAC(ITIME)
      RHO  = RHO + DRHO
	ENDIF

	CALL CLEARA (A(LDI),NEQ)
C	--------------------------------
	DO 3000	IITER = 1,NITER
C	--------------------------------

	CALL VECNORM(A(LDI),A(LDI),DNORM1,NEQ)

      ITASK = 2  !CALL STIFF AND INTERNAL FORCES	(YES)
	IF(IITER.EQ.1) ITASK = 1

	IFREF = 0  !STRUCTURE STIFFNESS				(YES)
      IFEIG = 1  !EIGEN STIFFNESS					(NO)

      IFPRI = 0  !STEP   PRINT					(NO)
      ISPRI = 1  !STRESS PRINT					(NO)

	CALL CLEARA (A(LRE),NEQ)
      CALL GRLOOP (IA(LEL),KSC)

	CALL MOVE(A(LLV),A(LDL),NEQ)
C	ADDITIVE MOMENT (MAKE CONSERVATIVE)
	IF(ITESR.EQ.1) CALL MOTRAN (IA(LID),A(LDT),A(LDI),A(LDL),NSF)

C     ------------------------------------------------
C     PERFORM --> EXTERNAL LOAD = RHO * REFERENCE LOAD
C     ------------------------------------------------
	CALL LOADFC(IA(LID),A(LLF),DREF,A(LDL),A(LDK),
     1            LMQ1,LMQ2,RHO,NEQ,ITIME,ITOPT,0,KONEQ,DRHO) 
C     ------------------------------------------------

	FAC = 1.0D0
	CALL ADDROW(A(LDL),A(LRE),A(LDL),FAC,NEQ)

      CALL COLSOL (IA(LMA),AA,A(LDK),A(LDL),1,INDPD,'STIF','TEMP')  !FACTORIZE
	CALL COLSOL (IA(LMA),AA,A(LDK),A(LDL),2,INDPD,'TEMP','TEMP')  !BACK SUBSTITUTION

C	---------------------------------------------------------------------
	IF(ITOPT.EQ.1) GOTO 2500
	CALL MOVE(A(LLV),DREF,NEQ)
	CALL COLSOL(IA(LMA),AA,A(LDK),DREF,2,INDPD,'TEMP','TEMP')		!BACK SUBSTITUTION	

	IF(IITER.EQ.1) THEN  !FIRST STEP
	CALL LOADFC(IA(LID),A(LLF),DREF,A(LDL),A(LDK),
	1            LMQ1,LMQ2,RHO,NEQ,ITIME,ITOPT,1,KONEQ,DRHO) 
	ENDIF

	IF(IITER.GT.1) THEN  !ITERATION
	CALL LOADFC(IA(LID),A(LLF),DREF,A(LDL),A(LDK),
	1            LMQ1,LMQ2,RHO,NEQ,ITIME,ITOPT,2,KONEQ,DRHO) 
	ENDIF
2500	CONTINUE
C	---------------------------------------------------------------------

	FAC = 1.0D0
	CALL ADDROW(A(LDI),A(LDL),A(LDI),FAC,NEQ)
C	ADD ROTATION DUE TO LARGE ROTATION
	IF(ITESR.EQ.1) CALL ADDROT(IA(LID),A(LDT),A(LDI),A(LDL),NSF,ISO)

	CALL VECNORM(A(LDI),A(LDI),DNORM2,NEQ)
	IF((DNORM2-DNORM1)/DNORM2.LT.TOLD) GOTO 3100

C	--------------------------------
3000	CONTINUE
C	--------------------------------

3100	FAC = 1.0D0
	CALL ADDROW(A(LDT),A(LDI),A(LDT),FAC,NEQ)

	CALL TIMPRN(ITIME,IITER,IA(LID),A(LDT),RHO,LMEQ1,LMEQ2,LMQ1,LMQ2,
	1			DRHO)

C	--------------------------------
5000	CONTINUE
C	--------------------------------


      RETURN
	END

C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================
      SUBROUTINE TIMPRN (KSTEP,ITER,ID,DISP,RHO,LMEQ1,LMEQ2,LMQ1,LMQ2,
	1				   DRHO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC      
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION ID(NSF,1),DISP(1)
	DIMENSION LMEQ1(1),LMEQ2(1),DELTA(10)

	II = 0
	DO I = 1,10
	DELTA(I) = 0.0D0
	IF(LMEQ1(I).NE.0.AND.LMEQ2(I).NE.0) THEN
      IEQ   = ID(LMEQ1(I),LMEQ2(I))
      DELTA(I) = DISP(IEQ)
	II = II + 1
	ENDIF
	ENDDO

	WRITE (*  ,1300) KSTEP,ITER,LMQ1,LMQ2,DELTA(1:II),RHO,DRHO
	WRITE (110,1300) KSTEP,ITER,LMQ1,LMQ2,DELTA(1:II),RHO,DRHO


 1300 FORMAT (2X,I5,X,I3,2X,2I5,3X,10E13.4)

	RETURN
      END
C
C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================
      SUBROUTINE LOADFC (ID,RFAC,RE,DISPL,DD,LMQ1,LMQ2,RHO,NEQ,KSTEP,
	1				   ITOPT,IENTRY,KONEQ,DRHO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     SELECTS INITIAL LOAD INCREMENT (IENTRY<2) OR LOAD INCREMENT
C     FOR NEXT ITERATIVE CYCLE (IENTRY=3)
C	-----------------------------------
C     RFAC(NSTEP)  = PRESCRIBED INCR. LOAD FACTORS / DISP. INCREMENTS
C                    (RFAC(KSTEP)=0 INDICATES AUTOMATIC SELECTION)
C     R(NEQ)       = REFERENCE LOAD VECTOR
C     RE(NEQ)      = DISPLACEMENT FIELD CORRESPONDING TO A REFERENCE
C                    LOAD OF ARBITRARY MAGNITUDE (DIRECTION OF MOTION)
C     DISPL(NEQ)   = EXTERNAL APPLIED REFERENCE LOAD VECTOR (IENTRY=0)
C                    LAST DISP. INCREMENT DUE TO RES. LOAD (IENTRY>0)
C     DD(NEQ)      = DIAGONAL TERMS OF INVERTED STIFFNESS MATRIX
C     KSTEP        = STEP COUNTER
C     IENTRY       = ENTRY POINT
C     IENTRY = 0     SET INITIAL LOAD INCREMENT BEFORE BACK-SUBST.
C     IENTRY = 1     SELECT INITIAL LOAD INCREMENT AFTER BACK-SUBST.
C     IENTRY = 2     SELECT INCREMENT FOR NEXT ITERATIVE CYCLE
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
C
      DIMENSION ID(NSF,1),RFAC(1),RE(1),DISPL(1),REF(NEQ)
	DIMENSION DD(1)
C
      DU = 0.0D0
      IF (IENTRY-1) 100,200,500
C     ------------------------------------------------
C     PERFORM --> EXTERNAL LOAD = RHO * REFERENCE LOAD
C     ------------------------------------------------
100	CALL VECMUL (DISPL,DISPL,DISPL,RHO,RNORM,NEQ)  

      RETURN
C     -----------------------------------------------------------
C     CHECK VALIDITY OF PRESCIBED DOF AND INITIALIZE EQUATION NO.
C     -----------------------------------------------------------
 200  NODE = LMQ2
      IPOS = LMQ1
      IF (NODE.LT.1 .OR. NODE.GT.NSN) GOTO 210
      IF (IPOS.LT.1 .OR. IPOS.GT.NSF) GOTO 210
      KONEQ = ID(IPOS,NODE)
      IF (KONEQ.NE.0) GOTO 250
C
 210  IF (ITOPT.EQ.2) ITOPT = 3
      CALL VECMAX (RE,NEQ,REMAX,KONEQ)
      DO 220  ISN=1,NSN
      DO 220  ISF=1,NSF
      IF (KONEQ.NE.ID(ISF,ISN)) GOTO 220
      LMQ2 = ISN
      LMQ1 = ISF
      GOTO 250
 220  CONTINUE

C     -------------------------------------------------------
C     CHOOSE APPROPRIATE SIGN FOR INITIAL TANGENTIAL SOLUTION
C     -------------------------------------------------------
 250  IF (ITOPT.LE.2) GOTO 300
      S = 1.0D0
      DO 270  IEQ=1,NEQ
      IF (DD(IEQ)) 280,270,270
 270  CONTINUE
      GOTO 290
 280  S = -1.0D0
C     ---------------------------------
C     SELECT APPROPRIATE NODE (ITOPT=3)
C     ---------------------------------
 290  REMAX = 0.0D0
	DO 277 ISN=1,NSN
	DO 276 ISF=1,NSF
	IEQ = ID(ISF,ISN)
	IPO = IDOF(ISF)
	IF (IEQ .LE. 0 .OR. IPO .LE. 0 .OR. IPO .GT. 3) GOTO 276
	IF (DABS(RE(IEQ)) .LT. DABS(REMAX)) GOTO 276
	REMAX = RE(IEQ)
	KONEQ = IEQ 
 276	CONTINUE
 277	CONTINUE
C     -----------------------------------------------------------
C     DETERMINE INCREMENTAL DISP./NORM TO BE APPLIED AT THIS STEP
C     -----------------------------------------------------------
 300  DU = RFAC(KSTEP)
      IF (ITOPT.EQ.3) DU = SIGN(DU,REMAX*S)
C     ---------------------------------------------------------
C     SELECT LOAD INCREMENT FOR NEXT ITERATIVE CYCLE
C     SET FINAL DISPLACEMENT INCREMENTS DISPL = DISPL + DRHO*RE
C     ---------------------------------------------------------
 500  DISUNB = DISPL(KONEQ)
      DISTAN = RE(KONEQ)
      DRHO = (DU-DISUNB)/DISTAN

C
 590  RHO = RHO + DRHO
      CALL VECMUL (RE,RE,REF,DRHO,DNORM,NEQ)
      CALL VECADD (REF,DISPL,DISPL,DRHO,DNORM,NEQ)

C	IF(IENTRY.EQ.2) THEN
C	WRITE(*,111) 'KAK',DRHO,DISTAN,DISUNB,RFAC(KSTEP)
C	WRITE(*,*) 
C111	FORMAT(A3,4E13.5)
C	PAUSE
C	ENDIF

C
      RETURN
      END
C
C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================
