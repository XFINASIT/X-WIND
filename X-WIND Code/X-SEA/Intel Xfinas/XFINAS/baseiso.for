C=====================================================================
	SUBROUTINE MPWISE(W,DRATIO,COFF,F,G,YT,YDOT,XI,XDOT,
	1			      NMODE,DELT,NSTEP,NOUT,TOL,NITER,FMODG,VECTIM)
  
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	(Programmed25-10-2002), Revised Nov19-02, and Jan22-2003
C	-------------------------------------------------------------------
C	PROGRAM TO SOLVE THE MODAL EQUATION OF MOTION FOR "NON-PROPORTIONAL
C	DAMPING CASES" BY "PSEUDO-FORCE METHOD" AND PIECE-WISE EXACT METHOD
C	(using Formulas in pp.112"Dynamics of Structures",Clough & Penzien)
C	-------------------------------------------------------------------
C	INPUT ARGUMENTS:
C	----------------
C         W      = MODAL FREQUENCIES ARRAY
C         DRATIO = MODAL DAMPING RATIOS ARRAY
C         COFF   = OFF-DIAGONAL MODAL DAMPING MATRIX
C         F      = MODAL LOADING AMPLITUDES
C         G      = DISCRETIZED EXCITATION RECORD
C         XI     = INITIAL CONDITIONS FOR MODAL DISPLACEMENTS
C		XDOT   = INITIAL CONDITIONS FOR MODAL VELOCITY
C		NMODE  = NUMBER OF MODAL EQUATIONS
C		DELT   = CONSTANT TIME STEP FOR COMPUTATION,= SUB-TIME INTEVAL
C         NSTEP  = TOTAL NUMBER OF TIME STEPS FOR OUTPUT REPORT
C		NOUT   = NUMBER OF SUBSTEPS BETWEEN OUTPUT REPORT
C		TOL    = SPECIFIED TOLERANCE FOR ITERATIVE PROCESS
C		NITER  = MAXIMUM NUMBER OF ITERATIONS
C		ET     = WORKING ARRAY OF LENGTH 2*NMODE
C	----------------
C	OUPUT ARGUMENTS:
C	----------------
C		Y      = MODAL DISPLACMENTS
C		YDOT   = MODAL VELOCITY IN THE LAST STEP
C	-------------------------------------------------------------	
	DIMENSION W(NMODE),DRATIO(NMODE),COFF(NMODE,NMODE),F(NMODE)
      DIMENSION G(NSTEP*NOUT+1),Y(NMODE,NSTEP+1),YDOT(NMODE)
      DIMENSION XI(NMODE),XDOT(NMODE),ET(NMODE,2)
      DIMENSION YT(NSTEP+1,NMODE)
	DIMENSION FMODG(NMODE,3),VECTIM(3,1),AAG(3),BBG(3)
C	------------------------
C	(A)-INITIAL COMPUTATIONS
C	------------------------
	ZERO = 0.0
	CALL CLEARMAT(ET,NMODE,2)

	CALL CLEARMAT(Y,NSTEP+1,NMODE)
	CALL CLEARA(YDOT,NMODE)
c
	DO 100 NM = 1,NMODE
		WDT = DELT*W(NM)*SQRT(1.-DRATIO(NM)**2)
		ENM = EXP(-W(NM)*DRATIO(NM)*DELT)
		ET(NM,1) = ENM*COS(WDT)
		ET(NM,2) = ENM*SIN(WDT)
  100 CONTINUE
C	------------------
C	(B)- FOR EACH STEP
C	------------------
	DO 700 NS = 1,NSTEP
	DO 700 NO = 1,NOUT
C
C	(B1)- Set initial conditions for iteration process
	DO 200 NM = 1,NMODE
		YDOT(NM) = XDOT(NM)
  200 CONTINUE
C	(B2)-Compute loading
C	next line changed Jan.13,2003 due to violating default NC
C	NC = (NS-1)*NOUT + NO
	NCC = (NS-1)*NOUT + NO

	AA = G(NCC)
	BB = (G(NCC+1) - G(NCC))/DELT

C	SONGSAK NEW EARTHQUAKE SEP2007
	AAG(1) = VECTIM(1,NCC)
	AAG(2) = VECTIM(2,NCC)
	AAG(3) = VECTIM(3,NCC)
	BBG(1) = (VECTIM(1,NCC+1)-VECTIM(1,NCC))/DELT
	BBG(2) = (VECTIM(2,NCC+1)-VECTIM(2,NCC))/DELT
	BBG(3) = (VECTIM(3,NCC+1)-VECTIM(2,NCC))/DELT

C	------------------------------------
C	(C)- ITERATE THROUGH NUMBER OF MODES
C	------------------------------------
	DO 600 NI = 1,NITER
	UNBAL1 = ZERO
C	
	DO 500 NM = 1,NMODE
		A = F(NM)*AA 
		B = F(NM)*BB

C	SONGSAK NEW EARTHQUAKE SEP2007
	A = A+FMODG(NM,1)*AAG(1)+ FMODG(NM,2)*AAG(2)+FMODG(NM,3)*AAG(3)
	B = B+FMODG(NM,1)*BBG(1)+ FMODG(NM,2)*BBG(2)+FMODG(NM,3)*BBG(3)
C
C		(C1)- compute off-diagonal damping force
		C = DOT(COFF(1,NM),XDOT,NMODE)
		D = DOT(COFF(1,NM),YDOT,NMODE)
		D = (D-C)/DELT
C
C		(C2)- Update loading
		ASTAR  = A-C
		BSTAR  = B-D
C
C		(C3)- Compute coefficients
		WUD = W(NM)
		DR = DRATIO(NM)
		WD = WUD * DSQRT(1.-DR*DR)
C	
		A0 = ASTAR/(WUD*WUD)-(2.* DR*BSTAR)/(WUD**3)
		A1 = BSTAR/(WUD * WUD)
		A2 = XI(NM) - A0
		A3 = (XDOT(NM) + DR*WUD*A2-A1)/WD
C
C		(C4)- Compute new values for displacement and velocity
		VEL = A1 + (WD*A3 - DR*WUD*A2)*ET(NM,1)
	1         - (WD*A2 + DR*WUD*A3)*ET(NM,2)
		Y(NM,NS) = A0 + A1*DELT + A2*ET(NM,1) + A3*ET(NM,2)
C		
C		(C5)- Check for convergence
		UNBAL2 = DABS((VEL - YDOT(NM))/VEL)
		YDOT(NM) = VEL
		UNBAL1 = DMAX1(UNBAL1,UNBAL2)
C
  500	CONTINUE
C	
	IF (UNBAL1.LE.TOL) GOTO 400
C
  600 CONTINUE
	WRITE (*,2000) NS,NO,NITER,UNBAL1
 2000	FORMAT ('WARNING AT STEP =',I5,' SUBSTEP ='I5,/
     1        'AFTER ITERARION =',I5,' UNBALNCE ='E12.5)
C     --------------------------------------------------     	
C	(D)- SET INITIAL CONDITIONS FOR THE NEXT TIME STEP
C	--------------------------------------------------
  400 DO 300 NM = 1,NMODE
         XI(NM)   = Y(NM,NS)
	   XDOT(NM) = YDOT(NM)
  300 CONTINUE
C
  700 CONTINUE
C
	CALL CLEARMAT(YT,NSTEP+1,NMODE)
	YT = TRANSPOSE(Y)
C	
C	WRITE(100,*) '-----------------YT-----------------'
C	WRITE(100,*) YT
C     -------------------
	RETURN
	END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE ADDDAMP(ID,MAXA,SD,DAPS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

	DIMENSION ID(NSF,NSN),SD(6,NSN),MAXA(1),DAPS(NSN)


	DO II = 1,NSN
	DO J  = 1,NSF
	JJ = IDOF(J)
	IF (JJ.LE.6) THEN
	IEQ   = ID(JJ,II)
	SPFAC = SD(JJ,II)
	IF(IEQ.GT.0) THEN
	MI = MAXA(IEQ)

C	A(MI) = A(MI) + SPFAC*DAPS(II)
	VALV = SPFAC*DAPS(II)
	CALL MESTIF(IEQ,VALV,1,1,'WRT')

	ENDIF
	ENDIF
	ENDDO
	ENDDO



	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE INPSD2(SD2,NODE,IOPT2,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

	DIMENSION SD2(7,6)

	IF(IND.EQ.1) THEN
	DO II = 1,6
	SD2(1,II) = NODE
	SD2(6,II) = IOPT2

	ENDDO
	ENDIF

	

	IF(IND.EQ.2) THEN
	DO II = 1,6
	READ(ITI,*) (SD2(JJ,II),JJ=2,3)   !READ ONLY 2 COMPONENT STIFF2 AND FSY  (4,5 FOR DMAX AND QFY SEE ALSO IN FOCSPIG)
	ENDDO
	ENDIF


	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SPGPAS(EPSP,SIGP,EPS,SIG,YOUNG,YOUN2,EPSTN,YIELD,
	1				  KMODU,EPTTN) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
C	THIS MODULE DEVELOPED BY SONGSAK FOR BILINEAR MODEL OF BASE ISOLATION
C	=====================================================================

      
      HARDS = 0.0D0
	IF(YOUNG.NE.YOUN2) HARDS = YOUNG*YOUN2/(YOUNG-YOUN2)
	
	DEPS  = EPS - EPSP
	DSIG  = DEPS*YOUNG

	SIGO = SIGP

	SIG  = SIGO + DSIG 

	PREY  = YIELD !+ HARDS*EPSTN
	
      IF(YOUNG.EQ.0.0D0) GOTO 40

	IF(ABS(SIGO - HARDS*EPTTN).GE.PREY) GO TO 20
	ESCUR = ABS(SIG - HARDS*EPTTN) - PREY
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(DSIG)
	GO TO 30


20	IF(SIG - HARDS*EPTTN.GT.0.0.AND.DSIG.LT.0.0) GO TO 40
	IF(SIG - HARDS*EPTTN.LT.0.0.AND.DSIG.GT.0.0) GO TO 40


	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	SIGO  = SIGO + REDUC*DSIG + 
	1	    RFACT*DSIG*(1.0-YOUNG/(YOUNG+HARDS))
	EPSTN = EPSTN + RFACT*DSIG/(YOUNG+HARDS)
	EPTTN = EPTTN + RFACT*DSIG/(YOUNG+HARDS)
	EPSTN = ABS(EPSTN)
	KMODU = 2 !YOUN2

	GO TO 50	

40	KMODU = 1 !YOUNG
	SIGO  = SIGO + DSIG
50	SIGP  = SIGO

	EPSP = EPS

	SIG  = SIGP

C	IF(KMODU.EQ.1) THEN
C	WRITE(*,*) KMODU,SIG
C	PAUSE
C	ENDIF


	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================