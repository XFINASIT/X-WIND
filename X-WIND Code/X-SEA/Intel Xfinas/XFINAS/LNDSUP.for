C
C=====================================================================
      SUBROUTINE LNDSUP(W,AA,BB,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------------
C	PROGRAM FOR LINEAR DYNAMIC ANALYSIS BY VECTOR SUPERPOSITION METHODS
C		DAMPING + PROPORTIONAL
C			    + NONPROPORTIONAL
C		VECTOR BASES
C				+ SUBSPACE EIGENVECTORS
C				+ LANCZOS EIGENVECTORS
C				+ LD RITZ VECTORS
C	(Changed to Subroutine April 5,2003)
C     -------------------------------------------------------------------
	CHARACTER*80 TITLE
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
c     COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
c     +              NRED,KPOSD,DETK,DET1,DAVR
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /DYNA/ IMASS
      COMMON /INCO/ A0,A1,A2,A3,A4,A5,A6,A7,A8,ALFA,BETA,
	1              A11,A12,ALF,IOPT
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C     COMMON /ITER/ NSTEP,DLMAX,NPRIN,NDRAW,KONEQ,NIREF,ITEMAX,RTOL,
C     1              ITOPT,NUMREF,NUMITE,ITETOT,RHO,RHOP,ICONV,NOLIN,
C     2              LIMEQ(2),KSTEP,ETOL,RHOPREV
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C     COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
	COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON A(9000000),IA(9000000)

C	NEW EARTHQUAKE ANALYSIS SONGSAK SEP2007
	COMMON /LNDEQK/ LDEK,LPR2,LDESTP
	COMMON /LNDEQR/ EQAMP(3),EQGAP(3)
C	---------------------------------------------------------------
	DIMENSION  W(1)
      DIMENSION  RVEC(NEQ,NROOT),RVECT(NROOT,NEQ),REIGV(NROOT)
      DIMENSION  RFREQ(NROOT),RSTIF(NROOT,NROOT),RMAS(NROOT,NROOT)
      DIMENSION  RPDAM(NROOT,NROOT)
	DIMENSION  RDAM(NROOT,NROOT),COFF(NROOT,NROOT),DRATIO(NROOT)
      DIMENSION  RFAC(1),RLV(NEQ),FF(NROOT),F(NROOT),FMOD(NROOT)	  
	DIMENSION  TIMLOAD(NINC+1),XI(NROOT),XDOT(NROOT),CTIMA(NINC+1)
	DIMENSION  Y(NINC+1,NROOT),YD(NINC+1,NROOT),YDD(NINC+1,NROOT)     
      DIMENSION  DISP(NINC+1,NEQ),VELO(NINC+1,NEQ),ACCE(NINC+1,NEQ)
	DIMENSION  DLOAD(NINC*10+1),DNDIAG(NEQ)
      DIMENSION  RVER(NEQ,NROOT+2),REIR(NROOT+2),STDISM(NINC+1,NEQ)
	DIMENSION  YDOT(NROOT)
	DIMENSION  FMODG(NROOT,3),VECTIM(3,NINC*10+1)
	DIMENSION  SEF(NEQ,3),PER(10000),RES(30000)
	DIMENSION  RMATD(NROOT,NROOT)
	DIMENSION  RLO(NEQ)  !CONSTANT LOAd VECTOR
C     =================================================
	DIMENSION AA(1),BB(1),CC(1)


C     21-10-2002
C     ----------------------------------------------------------------------
C     LINEAR DYNAMICS ANALYSIS USING SUPERPOSITION OF EIGEN OR LD-RITZ BASIS
C     ----------------------------------------------------------------------
C	ARGUMENTS:
C	----------
C	NEQ = SIZE OF GLOBAL STIFFNESS/MASS MATRIX = NO.OF SYSTEM D.O.F
C	NROOT  = IVPRT    = NUMBER OF MODES TAKEN INTO ACCOUNT TO SUPERPOSE 
C	RVEC(NEQ,NROOT)   = REDUCED MATRIX OF EIGENVECTOR TO SUPERPOSE
C	RVECT(NROOT,NEQ)  = TRANSPOSED MATRIX OF RVEC
C	REIGV(NROOT)      = REDUCED EIGENVALUE ARRAY, CORRESPONDING TO [RVEC]
C	RFREQ(NROOT)      = REDUCED FREQUENCY ARRAY
C	RMAS(Nroot,Nroot) = MODAL MASS MATRIX []
C	RDAM(Nroot,Nroot) = MODAL DAMPING MATRIX []
C	DRATIO(NROOT)     = MODAL DAMPING RATIOS ARRAY (CXI)
C
C	XI(NROOT)         = INITIAL CONDITIONS FOR MODAL DISPLACEMENTS
C	XDOT(NROOT)       = INITIAL CONDITIONS FOR MODAL VELOCITY 
C	------------------------------------------------------------------
	KSTEP = 1
C	-----------------
C     SET CONTROL FLAGS
C     -----------------
      IFPRI = 0
      IFPLO = 0
C     ----------------
C     FORM MASS MATRIX
C     ----------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 5
      CALL GRLOOP (IA(LEL),KSC)

C	--------------------------------------------------------
C	SONGSAK NEW EARTHQUAKE SEP2007
	IF(LDEK.NE.0) THEN
	CALL SEIFVEC(IA(LMA),IA(LID),IDOF,NSN,NSF,NEQ,NWM,SEF,
	1			 PER,RES,NUMT,NEAQ,BB,'MASS')
	ENDIF

C     -------------------------
C     FORM DAMPING MATRIX: 
C     -------------------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 6
      CALL GRLOOP (IA(LEL),KSC)
 
C     ------------------------------------------
C     FORM TANGENTIAL STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------
	IFEIG = 1
	IFREF = 0
      ISPRI = 1
      ITASK = 1
      CALL GRLOOP (IA(LEL),KSC)

      NUMREF = NUMREF + 1
C
C	------------------------------
C	Clear full matrices and arrays
C	------------------------------
	CALL CLEARMAT(RVEC,NEQ,NROOT)
	CALL CLEARMAT(RVER,NEQ,NROOT+2)

	CALL CLEARMAT(RMAS,NROOT,NROOT)
	CALL CLEARMAT(RPDAM,NROOT,NROOT)
	CALL CLEARMAT(RDAM,NROOT,NROOT)

	CALL CLEARA(REIGV,NROOT)
	CALL CLEARA(REIR,NROOT+2)
	CALL CLEARA(RFREQ,NROOT)
C	--------------------------------------
C	CALL REFERRENCE LOAD VECTOR {RLV(NEQ)}
C	--------------------------------------
	CALL CLEARA (RLV,NEQ)
	CALL VECADD (RLV,A(LLV),RLV,1D0,1D0,NEQ)
C	--------------------------------------
C	CALL CONSTANT   LOAD VECTOR {RLO(NEQ)}
C	--------------------------------------    
	CALL CLEARA (RLO,NEQ)
	CALL VECADD (RLO,A(LLO),RLO,1D0,1D0,NEQ)    
	        
C	---------------------------------------------------------
C	OPTIONS FOR EIGEN SOLVER (Subspace ISOLV=1, Ritz ISOLV=4)
C	---------------------------------------------------------
	IF(ISOLV.EQ.1) THEN
	CALL STABIL(AA,BB,'STIF','MASS')
	  CALL MPFCAL(IA(LMA),IA(LID),A(LER),A(LEV),A(LDIM),
	1			    A(LFRE),NROOT,NITEM,BB,'MASS',KSC)
	ENDIF
C	-------------------	
	IF(ISOLV.EQ.2) THEN
	CALL LANC(W,IA(LID),IA(LMA),N11,N10,AA,BB,'STIF','MASS')
	  CALL MPFCAL(IA(LMA),IA(LID),W(N10),W(N11),A(LDIM),
	1			   A(LFRE),NROOT,NITEM,BB,'MASS',KSC)
	ENDIF
C	-------------------
	IF(ISOLV.EQ.4) THEN
	NRIZV = NROOT
	NCR   = NRIZV*(NRIZV+1)/2			
	CALL RITZ(IA(LID),IA(LMA),A(LDK),RLV,
	1		  RVER,REIR,NRIZV,NCR,AA,BB,'STIF','MASS')

	DO 605 J=1,NROOT
	  DO 604 I=1,NEQ
		RVEC(I,J) = RVER(I,J)
  604	  CONTINUE
	    REIGV(J)  = REIR(J)	
  605	CONTINUE

	GOTO 614	
	ENDIF

C	========================================
C	DETERMINE COEFFICIENTS OF MODAL EQUATION
C	----------------------------------------
C	FORM MATRIX OF COLUMN-WISE EIGENVECTOR	
C	----------------------------------------
	II=0
	DO 610 J=1,NROOT
	  DO 607 I=1,NEQ
		IF (ISOLV.EQ.1) RVEC(I,J) = A(LER+II)
		IF (ISOLV.EQ.2) RVEC(I,J) = W(N10+II)
	    II=II+1
  607	  CONTINUE		
  610	CONTINUE
C
C	---------------------------------------------
C	FORM EIGENVALUES VECTOR CORRESPONDING TO RVEC 
C	---------------------------------------------
	JJ=0
	DO 612 J=1,NROOT
	  IF (ISOLV.EQ.1) REIGV(J) = A(LEV+JJ)
	  IF (ISOLV.EQ.2) REIGV(J) = W(N11+JJ)  	
	  JJ=JJ+1	
  612	CONTINUE

C	------------------------------
C	FORM REDUCED FREQUENCY ARRAYS
C	-----------------------------
  614	DO 615 IROW=1,NROOT
	    RFREQ(IROW) = DSQRT(REIGV(IROW))
  615	CONTINUE

C	------------------------------
C	PRINT EIGENVALUE AND MODE SHAPE
C	-----------------------------
	CALL MODPRN  

C	--------------------------------------------------------------------
C	FORM MODAL MASS MATRIX [RMAS(Nroot,Nroot)] = [RVECT]T *[MASS]*[RVEC]
C	--------------------------------------------------------------------
	CALL VFVMUL (IA(LMA),BB,RVEC,RMAS,NEQ,NROOT,1,'MASS')
C	---------------------------------------------
C	FORM MODAL PROPORTIONAL DAMPING MATRIX 
C	[RPDAM(Nroot,Nroot)]=ALFA*[RMAS]+BETA*[RSTIF]
C     ---------------------------------------------
C	In Case BETA.NE.0.0, form Modal stiffness matrix RSTIF, if BETA= 0.0, 
C				not necessary to form RSTIF, and proceed to Mass Proport.
C
	IF (BETA.NE.0.0) GOTO 617
	GOTO 625

  617 CALL CLEARMAT(RSTIF,NROOT,NROOT)
	CALL VFVMUL (IA(LMA),AA,RVEC,RSTIF,NEQ,NROOT,1,'STIF')
C
C	Form Stiff proportional damping to Modal Damping
C	------------------------	
	DO 622 I=1,NROOT
		DO 620 J=1,NROOT
			RPDAM(I,J) = BETA*RSTIF(I,J)			 
  620		CONTINUE
  622 CONTINUE
C
C	Form Mass proportional damping to Modal Damping
C	------------------------

  625	IF (ALFA.NE.0.0) THEN 
	DO 630 I=1,NROOT
		DO 627 J=1,NROOT
			RPDAM(I,J) = RPDAM(I,J) + ALFA*RMAS(I,J)				 
  627		CONTINUE
  630 CONTINUE
  	ENDIF

C	--------------------------------------------------------------------
C	FORM MODAL MATERIAL DAMPING/ISOLATOR MATRIX 
C		[RMATD(Nroot,Nroot)] = [RVECT]T *[MATDAM]*[RVEC]
C	--------------------------------------------------------------------
	CALL CLEARMAT(RMATD,NROOT,NROOT)

	CALL VFVMUL (IA(LMA),CC,RVEC,RMATD,NEQ,NROOT,1,'DAMP')	
	
C	--------------------------------------------------------
C	FORM TOTAL MODAL DAMPING MATRIX [RDAM] = [RPDAM]
C	--------------------------------------------------------
C	NOTE: This damping matrix sum help account for off-diagonal 
C	   	  terms in the modal proportional damping matrix. That will 
C	 	  increase the accuracy of the solution using iterative	
C	Jan30/03
	DO 637 I=1,NROOT
		DO 635 J= 1,NROOT
			RDAM(I,J) = RPDAM(I,J) + RMATD(I,J)
  635		CONTINUE
  637 CONTINUE
C	----------------------------------------
C	FORM TOTAL OFF-DIAGONAL MODAL DAMPING MATRIX
C	----------------------------------------
  	CALL CLEARMAT(COFF,NROOT,NROOT)

	DO 639 I=1,NROOT
		DO 638 J= 1,NROOT
			COFF(I,J) = RDAM(I,J)
  638		CONTINUE
		COFF(I,I) = 0.0 
  639 CONTINUE	

C	--------------------------------------------------
C	COMPUTE MODAL DAMPING RATIOS ARRAY {DRATIO(NROOT)}
C	--------------------------------------------------
  640	DO I=1,NROOT
	DRATIO(I) = RDAM(I,I)/(2*RMAS(I,I)*RFREQ(I))
	END DO

	CALL SDPROP (RFREQ,DRATIO,RMAS,NROOT)

C	=================================================
C	DETERMINE LOADING AMPLITUDE VECTOR AND LOAD CURVE
C	-------------------------------------------------
C	--------------------------
C	CALL LOAD COEFFICIENT RFAC
C	--------------------------
	JJ=0
	DO 645 J=1,NSTEP
	  RFAC(J) = A(LLF+JJ)
	  JJ=JJ+1	
  645	CONTINUE
	RHO = RFAC(KSTEP)
C	----------------------------------------
C	FORM LOADING AMPLITUDE VECTOR {F(NROOT)}
C	----------------------------------------
C	Transpose of eigenvector matrix = RVECT
C     ----------------
	RVECT = TRANSPOSE(RVEC)

	FF(1:NROOT) = 0.0D0 
	CALL MATMULT(RVECT,RLV,FF,NROOT,NEQ,1,CONST,100)
	CALL  VECMUL(FF,F,F,RHO,VNORM,NROOT)

C	ADD CONSTANT LOAD PART
	FF(1:NROOT) = 0.0D0 
	CALL MATMULT(RVECT,RLO,FF,NROOT,NEQ,1,CONST,100)
	CALL  VECADD(F,FF,F,1D0,1D0,NROOT)    
	     
C	-------------------------------------------------
C	FORM MODAL LOADING AMPLITUDE VECTOR {FMOD(NROOT)}
C	-------------------------------------------------
	DO I=1,NROOT
	FMOD(I) = F(I)/RMAS(I,I)
	END DO

C	-------------------------------------------------
C	-------------------------------------------------
C	SONGSAK NEW EARTHQUAKE SEP2007
  	CALL CLEARMAT(FMODG,NROOT,3)
	CALL CLEARMAT(VECTIM,3,(NINC*10+1))	
	IF(LDEK.NE.0) THEN
	CALL MATMULT(RVECT,SEF,FMODG,NROOT,NEQ,3,CONST,100)
	DO I = 1,NROOT
	DO J = 1,NEAQ
	FMODG(I,J) = FMODG(I,J)/RMAS(I,I)
	ENDDO
	ENDDO
	ENDIF
C	-------------------------------------------------
C	-------------------------------------------------

C	======================================================
C	BEGIN THE LOOP TO SOLVE NMODE OF SINGLE MODAL EQUATION
C	------------------------------------------------------
C	Initialization of Response Matrices
c	-----------------------------------
	CALL CLEARMAT (DISP,(NINC+1),NEQ)
	CALL CLEARMAT (VELO,(NINC+1),NEQ)
	CALL CLEARMAT (ACCE,(NINC+1),NEQ)
C	-------------------------------------------------
C	Options for method to solve single modal equation
c	-------------------------------------------------
C	When ISOLOOP = 6 (Mode superposition method) then:
c		+ IOPT = 1: Using piece wise exact method to solve for 
c					non-proportional or proportional damping
c		+ IOPT = 2: Using NEWMARK for the proportional damping eq.
C
C		+ IOPT = 3: Using WILSON-theta for the proportional damping eq.
C					(not available now)
c

	IF (IOPT.EQ.1) GOTO 650
	IF (IOPT.EQ.2) GOTO 665
C
C	=================================================================
C	(16Dec,2002,)THIS PART IS TO SOLVE SINGLE MODAL EQUATION
C	 USING MPWISE- PIECE WISE EXACT METHOD
C	Revised Jan 22,2003
C	------------------------------------------
C	FORM LOAD DISCRETIZATION ARRAY PROPLD(CTIM)
C	------------------------------------------
C	Maximum time of calculation
C	---------------------------	
  650	TMAX = NINC*DDT
C
C	Number of Sub-time interval (No. of Sub-interval for one DDT,=NOUT)
C	---------------------------
	NOUT = 10
	NSSIN = NINC*NOUT
C
C	Sub-time inteval
C	----------------
C	SDDT = DDT/NOUT
C	Previous line changed to the next lines 3Nov2003 by NguyenDV
	FOUT = DFLOAT(NOUT)
	DELT = DDT*(1.0/FOUT)
C
C	-----------
	CTIM = 0.0
	PROPLD = 0.0D0
	DO 655 I=1,NSSIN+1

	IF(LDEK.EQ.0) CALL LOADCURV(CTIM,2,NSTEP,A(LPR),PROPLD)
		DLOAD(I) = PROPLD

	IF(LDEK.NE.0) CALL SEIFVEP(CTIM,PER,RES,NUMT,NEAQ,VECTIM(1,I))  !SONGSAK NEW EARTHQUAKE SEP2007

		CTIM = CTIM + DELT

  655	CONTINUE
C	-----------------
C	Initial conditions of Modal Displacements and Velocities
C	-----------------
C	DO I=1,NROOT
C	  XI(I)  = 0.
C	  XDOT(I) = 0.
C	END DO
	CALL CLEARA(XI,NROOT)
	CALL CLEARA(XDOT,NROOT)

	TOL = 0.1E-6
C	----------------
	CALL MPWISE(RFREQ,DRATIO,COFF,FMOD,DLOAD,Y,YDOT,XI,XDOT,
	1			NROOT,DELT,NINC,NOUT,TOL,2*NC,FMODG,VECTIM)

C	--------------------------------------------------------
C	SUPERPOSE THE MODAL RESPONSE [DISP]=[RVECT]*[MODAL DISP]
C	--------------------------------------------------------
C	Initialize the TIMLOAD(NINC) for static correction
C	----------------------------
  	CTIM = 0.0
	DO 660 I=1,NINC+1
		CTIMA(I) = CTIM	
	    CALL LOADCURV(CTIM,2,NSTEP,A(LPR),PROPLD)
		TIMLOAD(I) = PROPLD
		CTIM = CTIM + DDT
  660 CONTINUE
C	--------------------------------------------
C	Form displacement matrix [DISP(NINC+1,NEQ)]
C	--------------------------------------------
	CALL MATMULT(Y,RVECT,DISP,(NINC+1),NROOT,NEQ,CONST,100)
C	------------------
C	Acount for Static Effects of (NEQ-NROOT) higher modes
C	Added Feb 18,2003	
C	------------------
	CALL CLEARMAT(STDISM,(NINC+1),NEQ)

	IF(ISOLV.EQ.1)CALL STACOR(IA(LMA),A(LDK),
     1					RLV,RVEC,TIMLOAD,STDISM,AA,'STIF')
C
	CALL MSPOUT (IA(LID),CTIMA,DISP,VELO,ACCE,STDISM)
C	-------------------------------------------
C	End of the piece wise exact method to solve for 
C	non-proportional or proportional damping eq.
C	-------------------------------------------	
C	RETURN
	GOTO 1000
C	===============================================================
C	THE FOLLOWING PART USING "NEWMARK TIME-INTEGRATION" METHOD
C	TO SOLVE SINGLE MODAL EQUATION
C	----------------------------------------------
C	FORM LOAD DISCRETIZATION ARRAY TIMLOAD(NINC+1)
C	----------------------------------------------	
  665	CTIM = 0.0
	PROPLD = 0.0D0
	DO 670 I=1,NINC+1

		CTIMA(I) = CTIM	

	IF(LDEK.EQ.0) CALL LOADCURV(CTIM,2,NSTEP,A(LPR),PROPLD)
		TIMLOAD(I) = PROPLD

	IF(LDEK.NE.0) CALL SEIFVEP(CTIM,PER,RES,NUMT,NEAQ,VECTIM(1,I))  !SONGSAK NEW EARTHQUAKE SEP2007

		CTIM = CTIM + DDT

  670	CONTINUE

C
  	CALL MNMARK(REIGV,FMOD,TIMLOAD,DRATIO,NROOT,NINC,DDT,Y,YD,YDD,
	1			FMODG,VECTIM)	        	  
C	==============================================================
C
C	--------------------------------------------------------
C	SUPERPOSE THE MODAL RESPONSE [DISP]=[RVECT]*[MODAL DISP]
C	--------------------------------------------------------
C	Form displacement matrix [DISP(NINC+1,NEQ)]
C	--------------------------------------------
	CALL MATMULT(Y,RVECT,DISP,(NINC+1),NROOT,NEQ,CONST,100)

C	Acount for Static Effects of (NEQ-NROOT) higher modes
C	Added Feb 18,2003	
C	------------------
	CALL CLEARMAT(STDISM,NINC+1,NEQ)

	IF(ISOLV.EQ.1)CALL STACOR (IA(LMA),A(LDK),
     1					RLV,RVEC,TIMLOAD,STDISM,AA,'STIF')
	
C	--------------------
C	Form Velocity matrix [VELO(NSSIN+1,NEQ)]
C	--------------------------------------------
	CALL MATMULT(YD,RVECT,VELO,(NINC+1),NROOT,NEQ,CONST,100)
C
C	Form Acceleration matrix [ACCE(NSSIN+1,NEQ)]
C	--------------------------------------------
	CALL MATMULT(YDD,RVECT,ACCE,(NINC+1),NROOT,NEQ,CONST,100)
C
C	-----------------------------------------------------------
C	EXTRACT DISP,VELO & ACCE AT SELECTED D.O.F OF SELECTED NODE
C	-----------------------------------------------------------
	CALL MSPOUT(IA(LID),CTIMA,DISP,VELO,ACCE,STDISM)

C	===============================	

1000	CONTINUE

C	--------------------------------------------------------------
C	ADD THIS BLOCK FOR PRINTING STRESS OF ALL STEP SONGSAK JUL2007

      ITASK = 3
      IFREF = 1
	IFEIG = 1
      ISPRI = 0
	DO INC=1,NINC+1

	CALL CLEROUT

	A(LDT:LDT+NEQ-1) = DISP(INC,1:NEQ)
C	NEW OUTPUT SONGSAK JUL2007
	CALL DISOUT(IA(LID),A(LDT))

      CALL GRLOOP (IA(LEL),KSC)

C	NEW OUTPUT SONGSAK JUL2007
	CALL PRNFLAG('ELEM','LINK','GSUP','LSUP','DISP','GSPG','LSPG')
	CALL  PRNOUT('STND','PONE','NONE',INC)

	ENDDO
C	--------------------------------------------------------------


      RETURN

	END
C     --------------------------------------------------------------------
C     END OF LINEAR DYNAMICS USING SUPERPOSITION OF EIGEN OR LD-RITZ BASIS
C
C=====================================================================
	SUBROUTINE MSPOUT (ID,CTIMA,DISP,VELO,ACCE,STDISM)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	Programmed Dec 20,2002, Revised & Changed to routine Jan 11,2003
C	--------------------------------------------------------------
C	PROGRAM TO EXTRACT DISP,VELO & ACCE FROM MODE SUPERPOSITION  
C	METHOD, AT A SELECTED D.O.F OF A SELECTED NODE, AND PRINT OUT
C	IN AN OUTJOB FILE 
C	--------------------------------------------------------------
C	INPUT ARGUMENTS:
C	----------------
C		ID    = 
C		CTIMA = DISCRETIZATION ARRAY OF TIME [CTIMA(NINC+1)]
C		DISP  = DISPLACMENTS MATRIX          [DISP(NINC+1,NEQ)]
C		VELO  = VELOCITY MATRIX              [VELO(NINC+1,NMOD)]
C		ACCE  = ACCELERATION MATRIX          [ACCE(NINC+1,NMOD)]
C
C	--------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

C	ADD THIS LINE FOR ITO SONGSAK MAY2006
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)


C      COMMON /ITER/ NSTEP,DLMAX,NPRIN,NDRAW,KONEQ,NIREF,ITEMAX,RTOL,
C     1              ITOPT,NUMREF,NUMITE,ITETOT,RHO,RHOP,ICONV,NOLIN,
C     2              LIMEQ(2),KSTEP,ETOL,RHOPREV
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
C	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /TIME/ DDT,CTIM,NINC
	COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C	-------------------------
      DIMENSION ID(NSF,1),CTIMA(NINC+1),DISP(NINC+1,NEQ),
     *		  VELO(NINC+1,NEQ),ACCE(NINC+1,NEQ),
     *		  DISNOD(NINC+1),VELNOD(NINC+1),ACCNOD(NINC+1),
     *		  STDNOD(NINC+1),STDISM(NINC+1,NEQ)	      
C	-----------------------------------------------------------
C	EXTRACT DISP,VELO & ACCE AT SELECTED D.O.F OF SELECTED NODE
C	-----------------------------------------------------------
	CALL CLEARA (DISNOD,(NINC+1))
	CALL CLEARA (VELNOD,(NINC+1))
	CALL CLEARA (ACCNOD,(NINC+1))

	CALL CLEARA (STDNOD,(NINC+1))
C	--------------------------------------
C	Select the D.O.F of a selected node to print out the response
C	--------------------------------------	
C	KDOF = 10
	KDOF = ID(LIMEQ(2),LIMEQ(1))
C	---------------------------------
C	Print out the response w.r.t time
C	---------------------------------
	IF(ISOLV.EQ.1) WRITE(100,691)
	IF(ISOLV.EQ.2) WRITE(100,692)

	IF(ISOLV.EQ.4) WRITE(100,694)
		
	DO IROW = 1,NINC+1

	DISNOD(IROW) = DISP(IROW,KDOF) + STDISM(IROW,KDOF)
	DISNOD(IROW) = DISP(IROW,KDOF)
	VELNOD(IROW) = VELO(IROW,KDOF)
	ACCNOD(IROW) = ACCE(IROW,KDOF)

	STDNOD(IROW) = STDISM(IROW,KDOF)

	WRITE(100,695)CTIMA(IROW),DISNOD(IROW),VELNOD(IROW),ACCNOD(IROW),
	1			  STDNOD(IROW)            
	WRITE(ITO,695)CTIMA(IROW),DISNOD(IROW),VELNOD(IROW),ACCNOD(IROW),
	1			  STDNOD(IROW)
      WRITE(10,695)CTIMA(IROW),DISNOD(IROW),VELNOD(IROW),ACCNOD(IROW),
	1			  STDNOD(IROW)
				 
	ENDDO

C
  691	FORMAT (//1H#,7X,45(1H*)/1H#,7X,1H*,43X,1H*/1H#,
     1  7X,45H* JOB PROGRESS BY MODE SUPERPOSITION METHOD */
     2  1H#,7X,1H*,4X,34H USING SUBSPACE EIGENVECTORS BASIS,5X,1H*/
     3  1H#,7X,45(1H*)/1H#,7X,'TIME',8X,'DISP',8X,'VELO',8X,'ACCE',
     4  6X,'STDISP')
C  
  692	FORMAT (//1H#,7X,45(1H*)/1H#,7X,1H*,43X,1H*/1H#,
     1  7X,45H* JOB PROGRESS BY MODE SUPERPOSITION METHOD */
     2  1H#,7X,1H*,5X,33H USING LANCZOS EIGENVECTORS BASIS,5X,1H*/
     3  1H#,7X,45(1H*)/1H#,7X,'TIME',8X,'DISP',8X,'VELO',8X,'ACCE',
     4  6X,'STDISP')
C
  694	FORMAT (//1H#,7X,45(1H*)/1H#,7X,1H*,43X,1H*/1H#,
     1  7X,45H* JOB PROGRESS BY MODE SUPERPOSITION METHOD */
     2  1H#,7X,1H*,5X,32H USING LOAD-DEPENDENT RITZ BASIS,6X,1H*/
     3  1H#,7X,45(1H*)/1H#,7X,'TIME',8X,'DISP',8X,'VELO',8X,'ACCE',
     4  6X,'STDISP')

  695	FORMAT(5E12.4)
                                                                       
	RETURN
      END
C
C=====================================================================
	SUBROUTINE MNMARK(REIGV,FMOD,TIMLOAD,DRATIO,NMOD,NINC,DDT,
     1	        	  Y,YD,YDD,FMODG,VECTIM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	
C
C     (17-Dec-2002, Revised and Changed to subroutine Jan 12,2003)
C	--------------------------------------------------------------
C	SOLVE THE NUMBER OF MODAL EQUATIONS OF MOTION BY NEWMARK TIME- 
C	INTEGRATION METHOD (Formulas in pp.177"Dynamics of Struc",A.K.Chopra)
C	--------------------------------------------------------------
C	INPUT ARGUMENTS:
C	----------------
C		REIGV  = EIGENALUES ARRAY
C         DRATIO = MODAL DAMPING RATIOS ARRAY
C		FMOD   = MODAL LOADING AMPLITUDE VECTOR {FMOD(NROOT)}
C		TIMLOAD= LOAD DISCRETIZATION ARRAY PROPLD(CTIM)
C		NMOD   = NUMBER OF MODES TAKEN INTO SUPERPOSITION
C		NINC   = NUMBER OF TIME INTERVAL DDT
C		DDT    = CONSTANT TIME STEP FOR COMPUTATION
C	---------------
C	LOCAL VARIABLES
C	---------------
C		DLMAT  = MODAL DISCRETIZED LOAD MATRIX
C		DELP   = LOAD INCREMENT
C		DELPBAR= EQUIVALENT LOAD INCREMENT
C		DELY   = DISPLACEMENT INCREMENT
C		DELYD  = VELOCITY INCREMENT
C		DELYDD = ACCELERATION INCREMENT
C	----------------
C	OUPUT ARGUMENTS:
C	----------------
C		Y    = MODAL DISPLACMENTS MATRIX [Y(NINC+1,NMOD)]
C		YD   = MODAL VELOCITY MATRIX     [YD(NINC+1,NMOD)]
C		YDD  = MODAL ACCELERATION MATRIX [YDD(NINC+1,NMOD)]
C	-------------------------------------------------------------
	DIMENSION REIGV(NMOD),DRATIO(NMOD),FMOD(NMOD),TIMLOAD(NINC+1),
	1		  FMODT(1,NMOD),DLMAT(NINC+1,NMOD),CM(NMOD),
	2		  ACOF(NMOD),BCOF(NMOD),CTIMA(NINC+1),PO(NMOD),	  	  
     3          DELY(NINC,NMOD),DELYD(NINC,NMOD),DELYDD(NINC,NMOD),	  
     4		  STIFBAR(NMOD),DELP(NINC,NMOD),DPBAR(NINC,NMOD),
     5		  Y(NINC+1,NMOD),YD(NINC+1,NMOD),YDD(NINC+1,NMOD)	
	DIMENSION FMODG(NMOD,3),VECTIM(3,1),DLMAG(NMOD,NINC+1)     
C	-------------
C	INIALIZATION
C	-------------
C	Coefficents of Newmark Liear Acceleration Method
C	-----------------------------------------------
	A11 = 0.5
	A12 = 0.25
c	A12 = 0.16666666667
C	-----------------
C	Initial conditions of Modal Displacements and Velocities
C	-----------------
	DO I=1,NMOD
	  Y(1,I)  = 0
	  YD(1,I) = 0
	END DO
C	------------------------------------------------
C	Form DLMAT((NINC+1),NROOT) matrix of total modal 
C	force discretized with time
C	------------------------------------------------
	CALL MATRAN(FMOD,FMODT,NMOD,1,1)
	CALL MATMULT(TIMLOAD,FMODT,DLMAT,(NINC+1),1,NMOD,CONST,100)

C	===============================================
C	SONGSAK NEW EARTHQUAKE SEP2007
	CALL CLEARMAT(DLMAG,NMOD,(NINC+1))
	CALL MATMULT(FMODG,VECTIM,DLMAG,NMOD,3,(NINC+1),CONST,100)
	DO I = 1,NINC+1
	DO J = 1,NMOD
	DLMAT(I,J) = DLMAT(I,J) + DLMAG(J,I)
	ENDDO
	ENDDO
C	===============================================

C	========================	
C	LOOP OVER NUMBER OF MODE
C	------------------------
	DO 685 J=1,NMOD
C	------------------------
C	CALCULATION OF PARAMETER
C	------------------------
C	Equivalent modal damping array CM(NMOD):
C
	CM(J) = 2*DRATIO(J)*SQRT(REIGV(J))
C
C	Modal Initial force (at CTIM = 0)
c
	PO(J) = DLMAT(1,J)
C
C	Modal Initial acceleration (at CTIM = 0 or I=1)
C
	YDD(1,J) = PO(J) - CM(J)*YD(1,J) - REIGV(J)*Y(1,J)
C
C	Equivalennt modal stiffness STIFBAR(NROOT):
C
	STIFBAR(J) = REIGV(J) + (1/A12)*(1/(DDT**2)) +
	1			 (A11/A12)*CM(J)*(1/DDT)
C
	ACOF(J) = (1/A12)*(1/DDT) + (A11/A12)*CM(J)
	BCOF(J) = (1/(2*A12)) + DDT*((A11/(2*A12))-1)*CM(J)
C	------------------------------
C	CALCUALTION FOR EACH TIME STEP
C	------------------------------
	DO 680 I = 1,NINC
C	I=1
	DELP(I,J) = DLMAT(I+1,J)-DLMAT(I,J)
C
C	Equivalennt modal stiffness:
C
	DPBAR(I,J) = DELP(I,J)+ACOF(J)*YD(I,J)+BCOF(J)*YDD(I,J)
C	
C	Increments of dis, velo. and acceleration:

	DELY(I,J)  = DPBAR(I,J)/STIFBAR(J)

	DELYD(I,J) = (A11/A12)*((DELY(I,J)/DDT) - YD(I,J))+
	1			  DDT*(1 - (A11/(2*A12)))*YDD(I,J)
C
	DELYDD(I,J)= (1/(A12*(DDT**2)))*DELY(I,J) - 
	1			 (1/(A12*DDT))*YD(I,J)- 
     2			 (1/(2*A12))*YDD(I,J)
C
C	Updation of new dis, velo. and acceleration:

	Y(I+1,J)  = Y(I,J) + DELY(I,J)	 
	YD(I+1,J) = YD(I,J) + DELYD(I,J)
	YDD(I+1,J)= YDD(I,J) + DELYDD(I,J)
  
  680 CONTINUE
  685	CONTINUE
C

C	==========================
	RETURN
      END
C
C=====================================================================
	SUBROUTINE CLEARMAT (AA,IA,JA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	Dec,18,2002
C     -----------------------------------------------------
C     PROGRAM TO CLEAR FLOATING FULL FORM MATRIX [A(IA,JA)]
C     -----------------------------------------------------
      DIMENSION AA(IA,JA)
      DO 10 I=1,IA
	DO 10 J=1,JA
   10    AA(I,J) = 0.0
      
	RETURN
      END
C
C=====================================================================
	SUBROUTINE CLEARMATI (IAA,IA,JA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	Dec,18,2002
C     -----------------------------------------------------
C     PROGRAM TO CLEAR FLOATING FULL FORM MATRIX [A(IA,JA)]
C     -----------------------------------------------------
      DIMENSION IAA(IA,JA)
      DO 10 I=1,IA
	DO 10 J=1,JA
   10    IAA(I,J) = 0
      
	RETURN
      END
C=====================================================================
C=====================================================================
	SUBROUTINE STACOR(MAXA,DD,SL,RVEC,TIMLOAD,STDISM,AA,TYP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C	Programmed Feb 18,2003 
C	--------------------------------------------------------------
C	PROGRAM TO COMPUTE THE EFFECTS OF HIGHER MODES THAT ARE NOT 
C			TAKEN INTO SUPERPOSITION BY STATIC CORRECTION METHOD
C	--------------------------------------------------------------
C	INPUT ARGUMENTS:
C	----------------
C		AA     =  STIFFNESS MATRIX [AA(NEQ,NEQ)] STORED IN ARRAY
C		BB     =  MASS MATRIX [BB(NEQ,NEQ)] STORED IN ARRAY
C		DD     =  DIAGONAL TERMS OF GLOBAL STIFFNESS DD[NEQ]
C		SL	   =  STATIC LOAD VECTOR = RFAC*REFFERNCE LOAD	 
C		RVEC   =  MATRIX OF COLUMN-WISE EIGEMVECTORS																																																						
C	----------------
C	OUTPUT ARGUMENTS:
C	----------------
C		STDISM = MATRIX CONTAINS STATIC EFFECT HISTORY
C				 OF (NEQ - NROOT)HIGHER MODES [NINC+1,NEQ]
C	--------------------------------------------------------------
C	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /TIME/ DDT,CTIM,NINC
	COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C	------------------------------------------------------------
	DIMENSION MAXA(NEQ),DD(NEQ),SL(NEQ),U(NEQ),
	1		  RVEC(NEQ,NROOT),RSTIF(NROOT,NROOT),PHI(NEQ),SCAPHI(NEQ),
     2		  DISLOW(NEQ),STADIS(NEQ),TIMLOAD(NINC+1),STDISA(NEQ),
	3		  STDISM(NINC+1,NEQ)
	DIMENSION AA(1)
C     --------------
C     INITIALISATION
C     --------------
      INDPD  = KPOSD
C	------------------------------------------------------------
C	EVALUATE TOTAL STATIC DISPLACEMENT RESPONSE OF ALL NEQ-MODES
C	------------------------------------------------------------	
	U(1:NEQ) = SL(1:NEQ)
	CALL COLSOL (MAXA,AA,DD,U,1,INDPD, TYP  ,'TEMP')
	CALL COLSOL (MAXA,AA,DD,U,2,INDPD,'TEMP','TEMP')        !GETTING DISPLACEMENT

C	----------------------------
C	FORM MODAL STIFFNESS MATRIX
C	----------------------------	
	CALL CLEARMAT(RSTIF,NROOT,NROOT)

	CALL VFVMUL (MAXA,AA,RVEC,RSTIF,NEQ,NROOT,0,TYP)
C	------------------------------------------	
C	COMPUTE STATIC EFFECT OF NROOT LOWER MODES
C	------------------------------------------	
	CALL CLEARA(DISLOW,NEQ)

	DO 50 IMODE = 1,NROOT

	RKN = RSTIF(IMODE,IMODE)
	CALL CLEARA(PHI,NEQ)
C	CALL CLEARA(PHIT,NEQ)

	DO 10 I = 1,NEQ
		PHI(I) = RVEC(I,IMODE)
C	    PHIT = TRANSPOSE(PHI)
C	    PHIT(1,I) = PHIT(1,I)/RKN
   10 CONTINUE

	CALL SCAPRD(PHI,SL,PRN,NEQ)
	PRKN = PRN/RKN

	CALL CLEARA(SCAPHI,NEQ)
	CALL SCAVEC(PHI,SCAPHI,PRKN,NEQ)

	DO 20 J=1,NEQ
   20	DISLOW(J) = DISLOW(J) + SCAPHI(J)
		
C	CALL MATMULT(PHI,PHIT,FN,NEQ,1,NEQ,CONST,100)   
   50 CONTINUE
C	-------------------------------------------------	
C	COMPUTE STATIC EFFECT OF (NEQ-NROOT) HIGHER MODES
C	-------------------------------------------------
	CALL CLEARA(STADIS,NEQ)

	DO 60 J=1,NEQ
   60	STADIS(J) = U(J) - DISLOW(J)
C	----------------------------------------------------------------
C	FORM MATRIX OF STATIC EFFECT HISTORY OF (NEQ-NROOT) HIGHER MODES
C	----------------------------------------------------------------
	CALL CLEARMAT(STDISM,(NINC+1),NEQ)

	DO 80 I=1,NINC+1

		TLOAD =TIMLOAD(I) 

		CALL SCAVEC(STADIS,STDISA,TLOAD,NEQ)
	    DO 70 J=1,NEQ
   70			STDISM(I,J) = STDISA(J)
   80 CONTINUE
		                                                                   
	RETURN

      END
C
C=======================================================================
	SUBROUTINE RITZ (ID,MAXA,DD,SL,EVECRZ,EVALRZ,NRV,NCR,AA,BB,
	1				 TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP1,TYP2
C	Nov 19,2002, Revised Dec,28,2002
C     -----------------------------------------------------------------
C     THIS PROGRAM IS USED TO
C	 + COMPUTE A LOAD DEPENDENT RITZ VECTOR BASIS (A BASIS CONCERNING
C	   THE SPATIAL DISTRIBUTION OF LOADS) FOR RAYLEIGH-RITZ METHOD
C	   OR RITZ METHOD 
C	 + FORM NUMBER OF APPROXIMATING EIGENVECTORS AND EIGENVALUES BASED
C	   ON RITZ VECTORS, WHICH ARE COMPATIBLE TO THE SPATIAL 
C	   DISTRIBUTION OF LOAD(S) OVER THE STRUCTURE
C     -----------------------------------------------------------------
C	INPUT ARGUMENTS
C	---------------
C		AA    =  STIFFNESS MATRIX [AA(NEQ,NEQ)] STORED IN ARRAY
C		BB    =  MASS MATRIX [BB(NEQ,NEQ)] STORED IN ARRAY
C		DD    =  DIAGONAL TERMS OF GLOBAL STIFFNESS DD[NEQ]
C		SL	  =  STATIC LOAD VECTOR = RFAC*REFFERNCE LOAD	 
C		T	  =  TEMPORARY STORAGE
C		NEQ   =  NUMBER OF EQUATION																																																												
C		NRV	  =  NUMBER OF RITZ VECTOR TO BE GENERATED,(TO BE = NROOT)
C		NSA   =  OPTION FOR STATIC RESIDUAL
C					NSA = 0; LAST VECTOR IS STATIC RESIDUAL
C					NSA = 1; NO STATIC RESIDUAL IN THE BASIS
C							
C		NWM   =  NUMBER OF COEFFICIENTS IN AA
C		NCR   =  NUMBER OF COEFFICIENTS IN ST,RM (RITZ-REDUCED 
C					STIFFNESS,MASS) IN UPPER TRIANGULAR FORMS
C	---------------
C	LOCAL VARIABLES
C	---------------
C		U	  =  STATIC RESPONSE- UPDATED VECTOR		
C		RIZ   =  RITZ COLUMN-WISE MATRIX (RITZ BASIS)[(NEQ,NRV)]
C	---------------
C	OUPUT ARGUMENTS
C	---------------	
C		EVALRZ(NRV)     = APPROXIMATING EIGENVALUES ARRAY
C		EVECRZ(NEQ,NRV) = APPROXIMATING EIGENVECTORS MATRIX
C	-----------------------------------------------------------------
C     COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
C
C	-----------------------------------------------------------
	DIMENSION MAXA(NEQ1),DD(NEQ),SL(NEQ),U(NEQ),RR(NRV,NRV),ID(1),
	1		  RIZ(NEQ,NRV),REVEC(NRV,NRV),EVECRZ(NEQ,NRV),
	2		  EVALRZ(NRV),UM(NEQ),T(NEQ),RM(NCR),RST(NCR),D(NRV)
	DIMENSION AA(1),BB(1)
C
C     --------------
C     INITIALISATION
C     --------------
      INDPD  = KPOSD
      TOLJ   = 1.0E-5
      NSMAX  = 12
	NSA    = 0
C	-------------------------------------
C	EVALUATE STATIC DISPLACEMENT RESPONSE
C	-------------------------------------
	U(1:NEQ) = SL(1:NEQ)
	CALL COLSOL (MAXA,AA,DD,U,1,INDPD, TYP1 ,'TEMP')
	CALL COLSOL (MAXA,AA,DD,U,2,INDPD,'TEMP','TEMP')        !GETTING DISPLACEMENT

C	=====================================
C	LOOP TO GENERATE "NRV" RITZ VECTORS
C	-------------------------------------
	DO 50 I=1,NRV

C	---------------------
C	Use static vector as last vector
C	---------------------
	IF ((I.EQ.NRV).AND.(NSA.EQ.0)) THEN

	DO 10 J=1,NEQ
   10	RIZ(J,I) = U(J)

	GO TO 30
	ENDIF

C	--------------------
C	SOLVE FOR NEW VECTOR
C	--------------------
C	Use previous disp.vector U(i-1) to form static load 
C	              {UM} = [M]{U}(i-1) for next step
C	------------------------------------------------
	CALL MAMULT(MAXA,BB,U,UM,TYP2,'STD') 

C	----------------
C	Solve static inertia eqs [K]{U}(i) = {UM}
C	----------------
	U(1:NEQ) = UM(1:NEQ)
	CALL COLSOL (MAXA,AA,DD,U,2,INDPD,'TEMP','TEMP')        !GETTING DISPLACEMENT

C	The output disp.vector U = {U}(i) is still not ortho. and not normal
C	
C	Pre-assign {U} as the ith Ritz vector
C	----------------	
	DO 25 J=1,NEQ
   25	RIZ(J,I) = U(J)
C	----------------
C	Orthogonalize the ith Ritz vector against all (i-1) 
C	previous Ritz vectors then normalize it w.r.t mass matrix
C	-------------
   30	CALL ORTRIZ(MAXA,RIZ,I,NEQ,NRV,NWM,BB,TYP2)

C	--------------------
C	Update static vector
C	--------------------
	CALL MAMULT(MAXA,BB,U,T,TYP2,'STD')
C
	CJ = DOT(RIZ(1,I),T,NEQ)

	DO 40 J=1,NEQ
   40 U(J) = U(J) - CJ*RIZ(J,I)
C
   50 CONTINUE
C
C	========================================================
C	PROJECT GLOBAL MASS AND STIFFNESS MATRIX ONTO RITZ BASIS	
C	--------------------------------------------------------		
C	Form Ritz-reduced stiffness matrix (stored in an array)
C	{RST(NWA)} = [RIZ]T *[STIFF]*[RIZ]
C	----------------------------
	CALL VFVMUL  (MAXA,AA,RIZ,RR,NEQ,NRV,0,TYP1)
	K = 0
	DO I = 1,NRV
	DO J = I,NRV
	K = K + 1
	RST(K) = RR(J,I)
	ENDDO
	ENDDO
C
C	-----------------------------------
C	Form Ritz-reduced mass matrix (stored in an array)
C	{RM(NWA)} = [RIZ]T *[MASS]*[RIZ]
C	------------------------------------
	CALL VFVMUL  (MAXA,BB,RIZ,RR,NEQ,NRV,0,TYP2)
	K = 0
	DO I = 1,NRV
	DO J = I,NRV
	K = K + 1
	RM(K) = RR(J,I)
	ENDDO
	ENDDO

C	=====================================================
C	SOLVE THE GENERAL EIGENPROBLEM OF RITZ-REDUCED SYSTEM
C	-----------------------------------------------------	
	CALL JACOBI(RST,RM,REVEC,EVALRZ,D,NRV,NCR,TOLJ,NSMAX)

C	--------------------------------------------------
C	FORM FINAL APPROXIMATING EIGENVECTORS MATRIX
C	 [EVECRZ(NEQ,NRV)]=[RIZ(NEQ,NRV)]*[REVEC(NRV,NRV)]
C	--------------------------------------------------
	CALL MATMULT(RIZ,REVEC,EVECRZ,NEQ,NRV,NRV,CONST,100)

C		EVALRZ(NRV)     = APPROXIMATING EIGENVALUES ARRAY
C		EVECRZ(NEQ,NRV) = APPROXIMATING EIGENVECTORS MATRIX	
	DO IRV = 1,NRV
	CALL RELFILL('EIVV',EVALRZ(IRV),1,IRV,1) !STORE EIGEN VALUE
	CALL MODOUT(ID,EVECRZ(1,IRV),IRV)        !STORE MODE SHAPE TO FILE
	ENDDO


	RETURN
      END
C	--------------------------------------
C	END OF SUBROUTINE RITZ,January,10,2003
C
C=======================================================================
C=====================================================================
	SUBROUTINE ORTRIZ (MAXA,RIZ,I,NEQ,NRV,NWM,BB,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C	Dec,28,2002
C     ---------------------------------------------------------------
C     PROGRAM TO 
C		+ ORTHOGONALIZE VECTOR "I" AGAINST "I-1" VECTORS W.R.T MASS
C		+ AND NORMALIZE VECTOR W.R.T MASS
C     ---------------------------------------------------------------
C	INPUT ARGUMENTS
C	---------------
C		BB    =  GLOBAL MASS MATRIX STORED IN ARRAY
C		NWM   =  NUMBER OF COEFFICIENTS IN BB
CC	---------------
C	OUPUT ARGUMENTS
C	---------------
C		RIZ(NEQ,I) = Ith RITZ VECTOR ORTHONORMALIZED TO PREVIOUS (I-1) VECTORS
C	---------------
	DIMENSION MAXA(*),RIZ(NEQ,NRV),T(NEQ),TT(NEQ),BB(1),SUMI(1)
C	---------------------------------------------------------
C	ORTHOGONALIZE VECTOR "I" AGAINST "I-1" VECTORS W.R.T MASS
C	---------------------------------------------------------
	IF (I.EQ.1) GO TO 55
	JJ = I-1
C
	DO 20 J=1,NEQ
   20	TT(J) = RIZ(J,I) 

   	CALL MAMULT(MAXA,BB,TT,T,TYP,'STD') 
C
	DO 50 J=1,JJ
		CJ = DOT(RIZ(1,J),T,NEQ)
		DO 45 L=1,NEQ
   45		RIZ(L,I) = RIZ(L,I) - CJ*RIZ(L,J)
   50 CONTINUE
C
C	----------------
C	NORMALIZE VECTOR
C	----------------
   55 SUMI(1) = 0.0

	DO 60 J=1,NEQ
   60	TT(J) = RIZ(J,I)

	CALL VFVMUL (MAXA,BB,TT,SUMI,NEQ,1,0,TYP)

	SUM = DSQRT(SUMI(1))

	DO 70 J=1,NEQ
   70 RIZ(J,I) = RIZ(J,I)/SUM					

	RETURN
      END
C
C=====================================================================

	SUBROUTINE SUMDAM (MAXA,AA,BB,CC,ALFA,BETA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C	Modified from RAYDAM subroutine by NguyenDV, Mar27-03
C     ----------------------------------------------------------------
C     FORM TOTAL SYSTEM DAMPING MATRIX [C] = [Cp] + [Cn] FROM 
C		+ MATERIAL DAMPING MATRIX [Cp] (RAYLEIGH - PROPORTIONAL)
C		+ EXTERNALLY NON-PROPORTIONAL DAMPING MATRIX [Cn]
C	OUTPUT OF DAMPING MATRIX IN [CC]
C	----------------------------------------------------------------
C     MAXA(NEQ1)= DIAGONAL ADDRESSES
C     AA(NWM)	  = STIFFNESS MATRIX STRORED IN COMPACTED FORM
C     BB(NWM)	  = MASS STRORED IN COMPACTED FORM
C	ALFA	  = MASS PROPORTIONAL CONSTANT
C	BETA	  = STIFFNESS PROPORTINAL CONSTANT
C	(NOTE: DAMPING MATRIX IS LUMPED FOR LUMPED MASS MATRIX IS USED)
C	-------
C	OUTPUT
C	-------
C	CC(NWK)	  = DAMPING STRORED IN COMPACTED FORM
C     ----------------------------------------------------------------
C	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C    +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL

C	----------------------------------------------------------------
	DIMENSION MAXA(NEQ1),AA(NWK),BB(NWM),CC(NWK),
     +		   DNDIAG(NEQ),DNARR(NWM)
C	
	CALL CLEARA(CC,NWK)
C	------------------------------
C	FORM MASS PROPORTIONAL DAMPING
C	------------------------------
	IF (ALFA.NE.0.0) GOTO 100
	GOTO 200

100	DO I = 1,NWM
		CC(I) = ALFA*BB(I)
	ENDDO

200	IF (BETA.NE.0.0) GOTO 300
	GOTO 500
C	-----------------------------------
C	FORM STIFFNESS PROPORTIONAL DAMPING
C	-----------------------------------
300	IF (NWM.EQ.NEQ) GOTO 400
C	------------------------------
C	FOR CONSTISTENT DAMPING MATRIX
C	------------------------------
	DO I = 1,NWK
		CC(I) = CC(I) + BETA*AA(I)
	ENDDO
	GOTO 500
C	-------------------------
C	FOR LUMPED DAMPING MATRIX
C	-------------------------
400	DO I = 1,NEQ
		JDIA = MAXA(I)
		CC(I) = CC(I) + BETA*AA(JDIA)
	ENDDO

C	--------------------
500	RETURN

	END
C
C====================================================================