C ====================================================================
      SUBROUTINE LNDTIM_SOIL_PILE(AA,BB,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     PROGRAM FOR LINEAR DYNAMICS BY DIRECT TIME INTEGRATION METHODS
C		+ ARBITRARY DYNAMIC LOAD (NORMAL LOADING)
C		+ SEISMIC LOADING BY DISCRETIZED GROUND ACCELERATION RECORD
C	Modified & Changed to subroutine,July25,2003 by NguyenDV
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      
      COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV,KFTTD,NSUPER 

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /DYNA/ CDEN,IMASS
      COMMON /INCO/ A0,A1,A2,A3,A4,A5,A6,A7,A8,ALFA,BETA,
	1              A11,A12,ALF,IOPT
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

C----	Next Common Block added 12Dec03 by NguyenDV 
	COMMON /SPBC/ NSS,NLSS

	COMMON /SEIS/ ISEIDA,ISEIOP,IADIR,IAFM,ICOMB,IRES,NPICK,
	1			  NSTIM,NEAQ,NAV,NSAV
	COMMON /LOCW/ N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,NSTI,NAX,NAY,
     +			  NAZ,NSTIN,NAXI,NAYI,NAZI
      
      COMMON / LOCALFORCE / DESIGNAXIAL(100000),DESIGNSHEART(100000),DESIGNSHEARS(100000)
     1                     ,DESIGNTORSION(100000),DESIGNMOMENTT(100000),DESIGNMOMENTS(100000)
C	Next commnon block added 15Apr07 by NguyenDV
	COMMON /RMAX/ DISMAX(2),VELMAX(2),ACCMAX(2)	   

C	ADDED BY SONGSAK APR2006 REACTION	
	COMMON /REACT/ LRC,LRCT,MFQ,LRID

C	NEW EARTHQUAKE ANALYSIS SONGSAK SEP2007
	COMMON /LNDEQK/ LDEK,LPR2,LDESTP
	COMMON /LNDEQR/ EQAMP(3),EQGAP(3)
      
      COMMON /COUANA/ COUPLEANALYSIS

      COMMON A(9000000),IA(9000000)
      
C     RELATIVE MOTION TOEY MAY 2015
      COMMON /GiDEle/ LGID 
      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN
      
      COMMON /SOILDYNA/ NSOIL
      COMMON / SOILRE/ NSOILRETURN
      COMMON/COUNT/ICOUNT,NCOUNT

C	Changed by Mr Van Nov06
C	COMMON /MEMW/ W(5000000)
	COMMON /MEMW/ W(7000000),IW(7000000)
      
      COMMON /DYNA_STEP/ INC
      
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      
      ! OFFSHORE LOADCASE BY TOEY 04/2018
      COMMON/OFFSHORE_CASE/ LGEN,IOFFL,IORI_OFFSHORE
C	----------------------------------------------------
	DIMENSION SEF(NEQ,3),PER(10000),RES(30000)

	DIMENSION AA(1),BB(1),CC(1)

	DIMENSION RFV(NEQ),RFC(NEQ)
      
C     --------------
C     INITIALISATION
C     --------------
      ! COMBINATION OF OFFSHORE LOADCASE AND GENREAL LOADCASE
      AA = 0.0D0
      BB = 0.0D0
      CC = 0.0D0
      IF (IOFFL.NE.0) CALL COMBINE_OFF_GEN (A(LLV),A(LLOF))
      
      REWIND(552)
      REWIND (550)
      READ (550,*) IFAST
      IF (KFAST.EQ.2) CALL DIRECTORY_OF_FILE_FAST
      
      NUMITE = 0
      NUMREF = 0
      INDPD  = KPOSD
C
C 500	KSTEP = 1
	KSTEP = 1
C     -----------------
C     SET CONTROL FLAGS
C     -----------------
      IFPRI = 0
      IFPLO = 0
      
C     ----------------
C      FORM MASS MATRIX
C     ----------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 5
      CALL GRLOOP (IA(LEL),KSC)
C	
C     -------------------------
C     FORM DAMPING MATRIX: 
C-----Added 16Nov03 by NguyenDV
C     -------------------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 6
      CALL GRLOOP (IA(LEL),KSC)

C     ------------------------------------------
C     FORM TANGENTIAL STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------
  505	IFEIG = 1
	IFREF = 0
      ISPRI = 1
      ITASK = 1
      CALL GRLOOP (IA(LEL),KSC)
      
      !CALL MASS_AND_STIFFNESS_BLOCK
      
      !CALL CPU_TIME (TIM1)

      NUMREF = NUMREF + 1
C	--------------------------------------------------------
C	SONGSAK NEW EARTHQUAKE SEP2007
	IF(LDEK.NE.0) THEN
	CALL SEIFVEC(IA(LMA),IA(LID),IDOF,NSN,NSF,NEQ,NWM,SEF,
	1			 PER,RES,NUMT,NEAQ,BB,'MASS')
	ENDIF

C     --------------------------------------------------------
C     FORM RAYLEIGH DAMPING MATRIX (ALFA OR/AND BETA .NE. 0.0)
C     --------------------------------------------------------
	CALL RAYDAM (IA(LMA),ALFA,BETA,NWK,NWM,NEQ,AA,BB,CC,
	1			 'STIF','MASS','DAMP')
C     ----------------------------------------------------------
C     FORM EFFECTIVE LINEAR STIFFNESS MATRIX [K]=[K]+A0[M]+A1[C]
C     ----------------------------------------------------------
	CALL EFSTIF (IA(LMA),NWK,NWM,NEQ,AA,BB,CC,
	1			 'STIF','MASS','DAMP','EFTF')
C     ----------------------------------------------------------
C     FACTORIZE EFFECTIVE LINEAR STIFFNESS MATRIX [K]=[L][D][L]T
C     ----------------------------------------------------------
      CALL COLSOL (IA(LMA),AA,A(LDK),A(LEF),1,INDPD,'EFTF','TEMP')
C     --------------------------------------------------
C     UPDATE INITIAL DISP., VELO. AND ACCE. (TIME T=0.0)
C     --------------------------------------------------
	IF (ISOLOP.EQ.9) THEN
		WRITE (100,1100)
		WRITE (ITO,1100)
          WRITE (10,1100)
	ELSE
		WRITE (100,1000)
		WRITE (ITO,1000)
          WRITE (10,1000)
	ENDIF

!      IF (COUPLEANALYSIS.EQ.3) CALL FORCEANDMOMENT (AM,BKG,0,0,0,0,IEL,ITASK,"STOR")
      
C	-------------------------
C	START TIME INCREMENT LOOP
C	-------------------------

      IF (KFAST.EQ.2) CALL FASTTOPFORCE (IA(LID),NSF,NODEA,NUMCASE,NFASTPARA,NTFASTPARA,RFC,NEQ,1,'REDT')
	CTIM = 0.0
      
C     -----------------------------------------------
      ! THIS PHASE FOR GENERATED TOTAL TIME STEP
      ! OPRES = 1 ; FREQ. LENGTH
      ! OPRES = 2 ; FREQ. HISTORY
	IF (KSPEC.EQ.1.AND.ICONTROLSPEC.EQ.1)THEN ! CASE; FREQ. ANALYSIS OF WAVE FORCE
	  IF (OPRES.EQ.1) CALL OFFSHSTEP(TIME,ITIME,NINC,'CALT') !READ TIME DATA FOR OFFSHORE LOAD 
        IF (OPRES.EQ.2) CALL SELECTGRAPH (1,NINC,0,0)
      ELSEIF (KSPEC.EQ.0.AND.ICONTROLSPEC.EQ.1) THEN ! CASE; NORMAL FREQUENCY 
        IF (OPRES.EQ.1)THEN
        NINC = ((STEPEND-STEPSTAT)/STEPINCR) + 1     
        ELSEIF (OPERS.EQ.2)THEN
        IF (OPRES.EQ.2) CALL SELECTGRAPH (1,NINC,0,0)    
        ENDIF
      ENDIF
C     -----------------------------------------------
      
	DO 5000 INC = 1,NINC
      
      NSOILRETURN = 0
      ICOUNT      = 0
      IF (ICONTROLSPEC.EQ.0) CTIM = CTIM + DDT
      
      IF (ICONTROLSPEC.EQ.1.AND.OPRES.EQ.1) CALL FREQ_STEP_INCR (INC,CTIM)
      IF (ICONTROLSPEC.EQ.1.AND.OPRES.EQ.2) CALL SELECTGRAPH (2,0,INC,CTIM)
      
100   IFEIG = 1
	IFREF = 0
      ISPRI = 1
      ITASK = 1
      CALL GRLOOP (IA(LEL),KSC)

C     ----------------
C     FORM MASS MATRIX
C     ----------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 5
      CALL GRLOOP (IA(LEL),KSC)
      
      TIMEDYNA =  CTIM ! CHANA
      
C     ----------------------------------------------
C     FORM CURRENT APPLIED FORCE VECTOR AT TIME T+DT
C     ----------------------------------------------
	CALL CPU_TIME (TIM1)
      
      REVATV = 2
      IF (COUPLEANALYSIS.EQ.2.OR.COUPLEANALYSIS.EQ.3) REVATV = 1
      
2100	IF(LDEK.EQ.0) THEN
	RFV(1:NEQ) = 0.0D0
	RFC(1:NEQ) = 0.0D0
      
      ! GENERAL LOADCASE FOR OFFSHORE (OLD VERSION BEFORCE 04/2018)
      CALL OFFSHFORC(RFV,1,INC,NEQ,'VARY','RADD') !VARY OFFSHORE LOAD     -- LOADCASE 1
      CALL OFFSHFORC(RFC,1,INC,NEQ,'CONT','RADD') !CONSTANT OFFSHORE LOAD -- LOADCASE 1
      CALL FASTTOPFORCE_DUMMY (IFAST,KFAST,NODEFAST,NTFASTPARA,NUMCASE,NFASTPARA)
      IF (KFAST.EQ.2) CALL FASTTOPFORCE (IA(LID),NSF,NODEA,NUMCASE,NFASTPARA,NTFASTPARA,RFC,NEQ,INC,'REDD')
      
      IF(KOREV.EQ.1) THEN  ! REALATIVE MOTION OF OFFSHORE 
      CALL RELATIVEMOTION (IA(LLM),A(LCO),IA(LGS),A(LGP),IA(LMS),A(LMP),IA(LGID),MFRMLD,IA(LHG),IA(LHGN),INC,A(LVL),A(LAL),RFV) 
      ENDIF
      
      REXTERNAL = 0.
      CALL RHSDYNA (A(LLF),A(LLV),A(LLO),RFV,RFC,REXTERNAL,A(LEF),A(LPR),A(LRT1),NEQ,2)
	ENDIF
	
	
C	--------------------------------------------------------
C	SONGSAK NEW EARTHQUAKE SEP2007
	IF(LDEK.NE.0) THEN
	CALL SEIFNEW(CTIM,PER,RES,NUMT,NEQ,NEAQ,SEF,A(LEF))
	ENDIF

C     -----------------------------------------------
C     FORM CURRENT EFFECTIVE LOAD VECTOR AT TIME T+DT
C     -----------------------------------------------
C      REPLACE FAST DISPLACEMENT,ROTATION,ACCELECTION,VELOCITY
      IF (REVATV.EQ.1) CALL WAITRESULT (IA(LID),NSF,A(LVL),A(LAL),A(LDT),"REAF")
  525 CALL EFLOAD (IA(LMA),A(LEF),A(LDT),A(LVL),A(LAL),A(LDU),
	1			 NEQ,NWK,NWM,AA,BB,CC,'STIF','MASS','DAMP')  
      IF (REVATV.EQ.1)  CALL WAITRESULT (IA(LID),NSF,A(LVL),A(LAL),A(LDU),"REAT")
C     CALCURATE REACTION FORCE AT INTERFACE POINT
      IF (REVATV.EQ.1) CALL FORCEANDMOMENT (AM,BKG,IA(LID),NSF,A(LVL),A(LAL),NELEMENT,ITASK,"CALU")

      CALL CPU_TIME (TIM2)
      TIM(10) = TIM(10) + (TIM2-TIM1)
C     ------------------------------------
C     SOLVE FOR DISPLACEMENTS AT TIME T+DT
C     ------------------------------------
	CALL COLSOL (IA(LMA),AA,A(LDK),A(LEF),2,INDPD,'TEMP','TEMP')
      
      IF (NSOILRETURN.NE.1) THEN
      IF (NSOIL.EQ.1) GOTO 100
      CALL CLEROUT
      ENDIF
      
C     ------------------------------------------
C     UPDATE DISP., VELO. AND ACCE. AT TIME T+DT
C     ------------------------------------------
	CALL MOVE (A(LDT),A(LDU),NEQ)
	CALL MOVE (A(LEF),A(LDT),NEQ)

      CALL CURDVA (IA(LID),A(LDT),A(LVL),A(LAL),A(LDU),NEQ,NSF,REVATV)
!      IF (REVATV.EQ.2) CALL WAITRESULT (IA(LID),NSF,A(LVL),A(LAL),A(LDT),A(LVL),"REAF")
      !IF (REVATV.EQ.2) GOTO 2100
      
C     -----------------------------------------------
C     CALCULATE AND PRINT NEW DISPLACEMENTS (IFPRI=0)
C     -----------------------------------------------
	CALL CLEROUT
590   CALL NEWDIS (IA(LID),A(LDI),A(LDT),NSF)		  


C	--------------------------------------------------------------
C	ADD THIS BLOCK FOR PRINTING STRESS OF ALL STEP SONGSAK JUL2007
      ITASK = 3
      IFREF = 1
	IFEIG = 1
      ISPRI = 0
      CALL GRLOOP (IA(LEL),KSC)

C	NEW OUTPUT SONGSAK JUL2007
      ! PRINT OUTPUT DATA OF DYNAMIC ANALYSIS
	CALL PRNFLAG('ELEM','LINK','GSUP','LSUP','DISP','GSPG','LSPG')
	CALL PRNOUT('STND','PONE','NONE',INC)
      
C	--------------------------------------------------------------

5000	CONTINUE  !MOVED HERE BY SONGSAK TO PRINT THE DISPLACEMENT AT EVERY TIME STEP
C	-------------------------
C	END TIME INCREMENT LOOP
C	-------------------------

C	-----------------------------------------------
C	OUTPUT OF EXTREMA RESPONSE FOR OVERALL ANALYSIS
C	Added 15Apr07 by NguyenDV
C	-----------------------------------------------
	WRITE(100,1500)	
      DO 800 I=1,2
        WRITE(100,1700) DISMAX(I),VELMAX(I),ACCMAX(I)
  800 CONTINUE

 1500 FORMAT (//5X,'EXTREMA VALUES (MAX & MIN) OF RESPONSES OVER WHOLE
     +TIME-HISTORY AT CURRENT DOFs'/
     +/26X,'DISP.EXTR.',2X,'VELO.EXTR.',2X,'ACCE.EXTR.')
 1700 FORMAT(24X,3(E12.4)) 
C     ------------------------------------------
C     CALCULATE AND PRINT STRESSES FOR LAST STEP
C     ------------------------------------------
      ITASK = 3
      IFREF = 1
      ISPRI = 0

      CALL GRLOOP (IA(LEL),KSC)
      
      IF (ICONTROLSPEC.EQ.1) CALL TRANFORM (NELE) 
      
      
      ! THE PRESENT VERSION OF X-SEA, WHICH IS USING GH_BLADED.EXE FOR ANALIZE FATIGUE 
C      CALL  FATIGUE_TIME_DOMAIN_CALAULATION ! OLD VERSION OF FATIGUE

      IF (KFAST.EQ.2) CALL CLEAR_DIRECTORY_OF_FILE_FAST
      
      RETURN
C
 1000	FORMAT (1H#,32X,16(1H*)/1H#,32X,1H*,14X,1H*/1H#,
     1  32X,16H* JOB PROGRESS */1H#,32X,1H*,14X,1H*/1H#,32X,16(1H*)/
     2  1H#,8X,'TIME',9X,'RHO',8X,'DISP',8X,'VELO',8X,'ACCE')

 1100	FORMAT (//1H#,5X,53(1H*)/
     + 1H#,5X,1H*,51X,1H*/
     + 1H#,5X,52H* LINEAR SEISMIC ANALYSIS BY DIRECT TIME-INTEGRATION,
     + 1H*/1H#,5X,1H*,51X,1H*/
     + 1H#,5X,1H*,19X,13H JOB PROGRESS,19X,1H*/
     + 1H#,5X,53(1H*)//
     + 1H#,8X,'TIME',9X,'RHO',8X,'DISP',8X,'VELO',8X,'ACCE')

      END
C
C=====================================================================s