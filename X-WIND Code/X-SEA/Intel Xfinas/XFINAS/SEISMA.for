C=====================================================================
      SUBROUTINE LNDTIM(AA,BB,CC)
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

	DIMENSION RFV(NEQ),RFC(NEQ),RFVO(NEQ),RFCO(NEQ),REXTERNAL(NEQ)
      
      DIMENSION ALOAD_GEN(NEQ)
C     ----------------------------------------------------
      
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
      
      CALL READ_EXTERNALFORCE_TIMESERIES (IA(LID))
C     -----------------------------------------------
	DO 5000 INC = 1,NINC

      IF (ICONTROLSPEC.EQ.0) CTIM = CTIM + DDT
      
      IF (ICONTROLSPEC.EQ.1.AND.OPRES.EQ.1) CALL FREQ_STEP_INCR (INC,CTIM)
      IF (ICONTROLSPEC.EQ.1.AND.OPRES.EQ.2) CALL SELECTGRAPH (2,0,INC,CTIM)
      
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

      DO IOFFCASE =1,IOFFL
      RFVO = 0.
      RFCO = 0.
      CALL OFFSHFORC(RFVO,IOFFCASE,INC,NEQ,'VARY','RADD') !VARY OFFSHORE LOAD     -- LOADCASE 1
      CALL OFFSHFORC(RFCO,IOFFCASE,INC,NEQ,'CONT','RADD') !CONSTANT OFFSHORE LOAD -- LOADCASE 1
      RFV(1:NEQ) = RFV(1:NEQ) + RFVO(1:NEQ)
      RFC(1:NEQ) = RFC(1:NEQ) + RFCO(1:NEQ)
      ENDDO
      CALL FASTTOPFORCE_DUMMY (IFAST,KFAST,NODEFAST,NTFASTPARA,NUMCASE,NFASTPARA)
      IF (KFAST.EQ.2) CALL FASTTOPFORCE (IA(LID),NSF,NODEA,NUMCASE,NFASTPARA,NTFASTPARA,RFC,NEQ,INC,'REDD')
      
      IF(KOREV.EQ.1) THEN  ! REALATIVE MOTION OF OFFSHORE 
      CALL RELATIVEMOTION (IA(LLM),A(LCO),IA(LGS),A(LGP),IA(LMS),A(LMP),IA(LGID),MFRMLD,IA(LHG),IA(LHGN),INC,A(LVL),A(LAL),RFV) 
      ENDIF
      
      REXTERNAL = 0.
      CALL EXTERNAL_FORCE (INC,REXTERNAL)
      
      !ALOAD_GEN = 0.
      !INDEX_IGEN_BEG = 0.
      !INDEX_IGEN_END = 0.
      !DO IGEN_CASE = 1,LGEN
      !INDEX_IGEN_BEG   = INDEX_IGEN_BEG + NEQ
      !INDEX_IGEN_END   = INDEX_IGEN_END + NEQ
      !IF (IGEN_CASE.EQ.1) INDEX_IGEN_BEG = LLV
      !IF (IGEN_CASE.EQ.1) INDEX_IGEN_END = LLV + NEQ - 1
      !ALOAD_GEN(1:NEQ)  = ALOAD_GEN(1:NEQ) + A(INDEX_IGEN_BEG:INDEX_IGEN_END)
      !ENDDO
      !! UPDATE 2020 FOR AUTOCAD
      !CALL RHSDYNA (A(LLF),ALOAD_GEN,A(LLO),RFV,RFC,A(LEF),A(LPR),A(LRT1),NEQ,2)
      
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
      !IF (REVATV.EQ.1) CALL FORCEANDMOMENT (AM,BKG,IA(LID),NSF,A(LVL),A(LAL),NELEMENT,ITASK,"CALU")

      CALL CPU_TIME (TIM2)
      TIM(10) = TIM(10) + (TIM2-TIM1)
C     ------------------------------------
C     SOLVE FOR DISPLACEMENTS AT TIME T+DT
C     ------------------------------------
	CALL COLSOL (IA(LMA),AA,A(LDK),A(LEF),2,INDPD,'TEMP','TEMP')
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
      IF (REVATV.EQ.1) CALL FORCEANDMOMENT (AM,BKG,IA(LID),NSF,A(LVL),A(LAL),NELEMENT,ITASK,"CALU")
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
C=====================================================================

C=====================================================================
	SUBROUTINE SDPROP (RFREQ,DRATIO,RMAS,NROOT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	Programmed Oct3,2003 by NguyenDV
C	--------------------------------------------------------------
C	PROGRAM TO PRINT STRUCTURAL DYNAMIC PROPERTIES
C		+ NATURAL CIRCULAR FREQUENCIES (Rad)
C		+ MODAL PERIOD (Sec.)
C		+ MODAL DAMPING RATIO
C		+ MODAL MASS
C	--------------------------------------------------------------
C	INPUT ARGUMENTS:
C	----------------
C		RFREQ
C		DRATIO
C		RMAS
C		NROOT
C	-------------------------------------------------
	DIMENSION RFREQ(NROOT),RPERI(NROOT),DRATIO(NROOT),
     +		  RMAS(NROOT,NROOT)
C
	WRITE (100,1000)

	DO IMOD = 1,NROOT
		RPERI(IMOD) = (2.*3.1416)/RFREQ(IMOD)
		WRITE (100,2000)IMOD,RFREQ(IMOD),RPERI(IMOD),DRATIO(IMOD),
     +					RMAS(IMOD,IMOD)
	ENDDO

 1000	FORMAT (//1H#,17X,33(1H*)/1H#,17X,1H*,31X,1H*/1H#,
     117X,33H* STRUCTURAL DYNAMIC PROPERTIES */
     21H#,17X,33(1H*)/
     31H#,1X,'MODE NO.',2X,'CIR.FREQ(Rad)',3X,'PERIOD(Sec.)',2X,
     4'MODAL DAMP.RATIO',3X,'MODAL MASS.')

 2000 FORMAT (3X,I3,6X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6,4X,E13.6)
                                                                       
	RETURN
      END
C
C=====================================================================
C=====================================================================
      SUBROUTINE INTCOS
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO CALCULATE CONSTANTS FOR TIME INTEGRATION METHODS 
C		(FOR ISOLOP = 9 :SEISMIC ANALYSIS ONLY) 
C	Modified from INTCON, by NguyenDV,Sep02/03 for Seismic Analysis
C	----------------------------------------------------------------
C     IOPT		TIME INTEGRATION OPTIONS
C					1 = WILSONS THETA METHOD
C					2 = NEWMARK METHOD
C					3 = CENTRAL DIFFERENCE METHOD
C					4 = HILBER-HUGHES-TAYLOR METHOD
C     DDT			TIME INCREMENT
C	CTIM		CURRENT TIME
C	NINC		NUMBER OF TIME INCREMENTS
C	IDAMP		DAMPING MATRIX OPTION
C	IMASSN		NODAL MASS OPTION
C	IDAMPN		NODAL DAMPING OPTION
C	A0-A8		TIME INTEGRATION CONSTANTS
C	ALFA		RAYLEIGH DAMPING PARAMETER (ALPHA)
C	BETA		RAYLEIGH DAMPING PARAMETER (BETA)
C	A11			DELTA FOR NEWMARK (OR THETA FOR WILSON THETA)
C	A12			ALPHA  FOR NEWPARK (OR THETA*DDT FOR WILSON THETA)
C	--------------------------------------------------------------
      COMMON /DYNA/ IMASS
      COMMON /TIME/ DDT,CTIM,NINC
      COMMON /INCO/ A0,A1,A2,A3,A4,A5,A6,A7,A8,ALFA,BETA,
	1              A11,A12,ALF,IOPT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
c
	COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
C
      DATA ZERO/1.E-12/
      DATA THMIN/1.39/,THMAX/2.01/
C     --------------------------------------------------------------
      GOTO (100,200,300,231,400),IOPT
C	-------------------
C     WILSON THETA METHOD
C	-------------------                           
  100	IF (A11.LE.0.0) A11=1.4
      TH = A11
      IF ((TH.GT.THMIN).AND.(TH.LT.THMAX)) GOTO 110
      WRITE (ISO,1002)
      STOP
  110	CONTINUE
      DTA=DDT*TH
      A0 = 6.0/(DTA*DTA)
      A1 = 3.0/DTA
      A2 = 2.0*A1
      A3 = 0.5*DTA
      A4 = A0/TH
      A5 =-A2/TH
      A6 = 1.0-3.0/TH
      A7 = 0.5*DDT
      A8 = DDT*DDT/6.0
	A12= DTA
	WRITE (ISO,1500) NINC,DDT,ALFA,BETA,A11
	RETURN
C	--------------
C     NEWMARK METHOD
C	--------------
  200	IF (A11.LE.0.0) A11 = 0.5
      IF (A12.LE.0.0) A12 = 0.25*((A11 + 0.50)**2)
      DELT = A11
      ALPA = A12
      IF (DELT.LT.0.50) GOTO 210
      IF (DELT.GT.0.55) GOTO 220                                       
      D1=0.5*(DELT+0.5)
      D2=D1*D1
      IF (ALPA.LE.ZERO) GOTO 210
      IF (ALPA.LT.D2) GOTO 230
      IF (ALPA.GT.D1) GOTO 220
      GOTO 240
  210	CONTINUE
      WRITE (ISO,1002)
      STOP
  220	CONTINUE
      WRITE (ISO,1004)                                                 
      GOTO 240
  230	CONTINUE
      WRITE (ISO,1005)
C	FOLLOWING LINES WERE ADDED BY GILSON - OCT2001
	GOTO 240
C	---------------------------
C	HILBER-HUGHES-TAYLOR METHOD
C	---------------------------
  231	IF (ALF.GE.-0.333.OR.ALF.LE.0.0) GOTO 232
	WRITE (ISO,1002)
	STOP
  232	IF (A11.EQ.0.0) THEN
	  DELT = (1.0-2.0*ALF)/2.0
	ELSE
	  DELT = A11
	ENDIF
	IF (A12.EQ.0.0) THEN
	  ALPA = 0.25*(1.0-ALF)**2
	ELSE
	  ALPA = A12
	ENDIF
C	---------------------------------
C	COMMON TO NEWMARK AND HHT METHODS
C	---------------------------------
240	CONTINUE
      DEAL=DELT/ALPA
      A0 = 1.0/(ALPA*DDT*DDT)
      A1 = DEAL/DDT
      A2 = 1.0/(ALPA*DDT)
      A3 = 0.5/ALPA-1.0
      A4 = DEAL-1.0
      A5 = DDT*(0.5*DEAL-1.0)
      A6 = DDT*(1.0 - DELT)
      A7 = DELT*DDT
	WRITE (ISO,2000) NINC,DDT,ALFA,BETA,A11,A12
	RETURN
C	--------------------------
C     CENTRAL DIFFENRENCE METHOD
C	--------------------------
  300	A0 = 1/DDT/DDT
	A1 = 1/2.0/DDT
	A2 = 2.0*A0
	A3 = 1/A2
	WRITE (ISO,3000) NINC,DDT,ALFA,BETA

  400	RETURN
C                                                                       
 1002 FORMAT(//,10X,'ERROR; TIME INTEGRATION PARAMETER NOT ADMISSIBLE')
 1003 FORMAT(//,10X,'ERROR; TIME STEP TOO SMALL')
 1004 FORMAT(//,10X,'WARNNING; TIME INTEGRATION PARAMETER SUSPICIOUS')
 1005 FORMAT(//,10X,'CONDITIONALLY STABLE NEWMARK ALGORITHM')

 1500	FORMAT(/19X,'TIME INTEGRATION PARAMETER FOR WILSON-THETA'//,
     1 11X,' INCR',2X,'TIME-STEP',2X,'RAY DAM-ALPHA',2X,'RAY DAM-BETA',
     2 6X,'THETA',/12X,55(1H-)/10X,I6,E11.4,4X,E11.4,3X,E11.4,E11.4)
 2000	FORMAT(/21X,'TIME INTEGRATION PARAMETER FOR NEWMARK'//,
     1 1X,' INCR',2X,'TIME-STEP',2X,'RAY DAM-ALPHA',2X,'RAY DAM-BETA',
     2 2X,'NEWMARK-DELTA',2X,'NEWMARK-ALPHA'/2X,74(1H-)/,
     3 I6,E11.4,4X,E11.4,3X,E11.4,4X,E11.4,4X,E11.4)
 3000	FORMAT(/15X,'TIME INTEGRATION PARAMETER FOR CENTRAL DIFFERENCE'//,
     1 16X,' INCR',2X,'TIME-STEP',2X,'RAY DAM-ALPHA',2X,'RAY DAM-BETA',
     2 /17X,44(1H-)/15X,I6,E11.4,4X,E11.4,3X,E11.4)
C                                                                       
	END
C
C=====================================================================
C=====================================================================
      SUBROUTINE REDYPST
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------
C     TO READ ADDITIONAL DYNAMIC PROPERTIES OF STRUCTURES:
C	- DRAT1,DRAT2 : Modal damping ratios of 1st and 2nd modes
C	This program added 08Apr07 by Nguyen
C	-----------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
c	COMMON /DYPR/ DRAT1,DRAT2 
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW	
C	-----------------------------------------------------
C	READ MODAL DAMPING RATIOS
	READ(ITI,*) 
c	READ(ITI,*) DRAT1,DRAT2
	READ(ITI,*) ZET1,ZET2
	WRITE(ITO,*)'READ MODAL DAMPING RATIOS'
      WRITE(10,*)'READ MODAL DAMPING RATIOS'

c      WRITE(ISO,1020) DRAT1,DRAT2
      WRITE(ISO,1020) ZET1,ZET2

 1020 FORMAT(//,1X,'INPUTTED MODAL DAMPING RATIOS OF STRUCTURE',/
     +1X,'Modal damping ratio of 1st mode = ',F10.5/
     +1X,'Modal damping ratio of 2nd mode = ',F10.5)

      RETURN	
	END

C=====================================================================
      SUBROUTINE VECNORM(A,B,VNORM,N)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------------------
C	PROGRAM TO FORM NORM PRODUCT OF TWO VECTORS
C	 VNORM = SUM(A(I)*B(I),I=1,N) = (A(1)*B(1)+...+(A(N)*B(N)
C	--------------------------------------------------------------
	DIMENSION A(1),B(1)

      VNORM = 0.
      DO 100  I=1,N
  100		VNORM = VNORM + A(I)*B(I)
      RETURN
	END
C
C=====================================================================
      SUBROUTINE VECFACT(A,FACT,N)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------------------
C	PROGRAM TO MULTIPLY A VECTOR WITH A FACTOR
C	 A(I) = A(I)*FACT
C	--------------------------------------------------------------
	DIMENSION A(1)

      DO 100  I=1,N
  100		A(I)= A(I)*FACT
      RETURN
	END
C
C	=====================================================================
C	======================================================SONGSAK AUG2007
C	=====================================================================
      SUBROUTINE RESSPE1(KND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------------

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
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
      COMMON /LOCW/ N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,NSTI,NAX,NAY,
     +			  NAZ,NSTIN,NAXI,NAYI,NAZI

	COMMON /MEMW/ W(7000000),IW(7000000)
      COMMON A(9000000),IA(9000000)

C	MOVING LOAD PRINTING FLAG SONGSAK JUL2007
	COMMON /LANPRIN/  LAN_PRN,ILPL,MLANE,NPL,LAN_OPT

	COMMON /STCAS/ ILC

C	----------------------------------------------------
	DIMENSION PER(10000),RES(30000),ANGAP(3),RMAS(NROOT,NROOT)
	DIMENSION CVEC(NEQ,NROOT),SFREQ(NROOT)
C	----------------------------------------------------
	DIMENSION AA(NWK),BB(NWM)

	SELECTCASE(KND)

	CASE(0)
C	READ THE EARTHQUEAK RESPONSE SPECTRUM DATA  
	CALL RESPDATA(PER,RES,ANGAP,IEQK,NUMT,0)
		 		
	READ(ITI,*) 
	READ(ITI,*) NUMR

	CALL  DEMREL('SEQK',KSEQK,NUMR,5) 
	DO IUMR = 1,NUMR
	READ(ITI,*) MLC,IEQK,ICOMB,ZET1,ZET2
	FMLC   = FLOAT(MLC)
	FIEQK  = FLOAT(IEQK)
	FICOMB = FLOAT(ICOMB)
	CALL MRELFIL('SEQK',FMLC  ,IUMR,1,1)
	CALL MRELFIL('SEQK',FIEQK ,IUMR,2,1)
	CALL MRELFIL('SEQK',FICOMB,IUMR,3,1)
	CALL MRELFIL('SEQK',ZET1  ,IUMR,4,1)
	CALL MRELFIL('SEQK',ZET2  ,IUMR,5,1)
	ENDDO

	CALL DEFNINT('RESP',KRESP,1,20)   !RESPONSE SPECTRUM ANALYSIS DATA   
	CALL INTZERO('RESP')
	NFIL = 851						  !RESPONSE SPECTRUM OUTPUT FILE
	CALL INTFILL('RESP',NFIL,1,1,1)     
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS
	NFIL = 852						  !DUMMY FOR MODe COMB FILE
	CALL INTFILL('RESP',NFIL,1,2,1)     
	CALL INIOPER(NFIL,0,0.0D0,'NONE','NEWF')  !DIRECT ACCESS

	RETURN



	CASE(1)

C	SONGSAK NEW OUTPUT JUL2007
	LAN_PRN = 0  ! 2 is no fixend force for response spectrum analysis

	ILC   = 0   !NO LOAD CASE FOR RESPONSE SPECTRUM  (THIS RELATED TO FIXEND FORCES AND TENDON)
C     -----------------
C     SET CONTROL FLAGS
C     -----------------
      IFPRI = 0
      IFPLO = 0
C	==========================================

C     ----------------
C     FORM MASS MATRIX
C     ----------------
	IFEIG = 0
	IFREF = 1
      ISPRI = 1
      ITASK = 5
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
C     -------------------------------------
C     SOLVE EIGEN VALUES FOR FREE VIBRATION
C     -------------------------------------
	IF (ISOLV.EQ.1) THEN
	CALL STABIL(AA,BB,'STIF','MASS')  !SUBSPACE
	  CALL MPFCAL(IA(LMA),IA(LID),A(LER),A(LEV),A(LDIM),
	1			    A(LFRE),NROOT,NITEM,BB,'MASS',KSC)
	ENDIF
	IF (ISOLV.EQ.2) THEN
	CALL LANC(W,IA(LID),IA(LMA),N11,N10,AA,BB,'STIF','MASS')
	  CALL MPFCAL(IA(LMA),IA(LID),W(N10),W(N11),A(LDIM),
	1			   A(LFRE),NROOT,NITEM,BB,'MASS',KSC)
	ENDIF
C	==========================================

	CALL MODPRN  !PRINT EIGENVALUE AND MODE SHAPE

	SELECTCASE(ISOLV)
	CASE(1)              !SUBSPACE ITERATION
	NEV = LEV 
	NER = LER 
	CASE(2)              !LANCZOS
	NEV = N11
	NER = N10 
	ENDSELECT

C	MODAL MASS = VT*M*V
	CALL VFVMUL (IA(LMA),BB,A(NER),RMAS,NEQ,NROOT,0,'MASS')
			 		
	CALL INTFILL('RESP',NFIL1,1,1,0) 

C	CALLING NUMR
	CALL  LOCATM('SEQK',KSEQK,NUMR,NCMR,2)

	DO 1000 I = 1,NUMR
	CALL CLEARA(A(LDL),NEQ)

C	READ(ITI,*) MLC,IEQK,ICOMB,ZET1,ZET2
	CALL MRELFIL('SEQK',FMLC  ,I,1,0)
	CALL MRELFIL('SEQK',FIEQK ,I,2,0)
	CALL MRELFIL('SEQK',FICOMB,I,3,0)
	CALL MRELFIL('SEQK',ZET1  ,I,4,0)
	CALL MRELFIL('SEQK',ZET2  ,I,5,0)
	MLC   = INT(FMLC)
	IEQK  = INT(FIEQK)
	ICOMB = INT(FICOMB)

	CALL RESPDATA(PER,RES,ANGAP,IEQK,NUMT,1)
	CALL  RESSPE2(PER,RES,ANGAP,NUMT,IA(LMA),IA(LID),A(NEV),A(NER),
	1			  RMAS,IDOF,NSN,NSF,NEQ,NROOT,NWM,ICOMB,CVEC,SFREQ,
     2			  BB,'MASS')

	DO IROOT = 1,NROOT
	DO II = 1,NEQ
	A(LDL+II-1) = CVEC(II,IROOT)
	ENDDO


	CALL CLEROUT  !CLEAR OUTPUT DATA ARRAY (XFINAS)

	CALL DISOUT(IA(LID),A(LDL))
C     ------------------------------------------
	CALL MOVE (A(LDL),A(LDT),NEQ)

      ITASK = 3
      IFREF = 1
      ISPRI = 0
      CALL GRLOOP (IA(LEL),KSC)

	CALL INTFILL('RESP',NFIL1,1,1,0)  
	LP1 = IROOT
	FACTOR = 1.0D0
	CALL INIOPER(NFIL1,LP1,FACTOR,'REPC','OLDF')  
	ENDDO

	CALL RESCOME(NROOT,SFREQ,ICOMB)
	CALL BACKLCS(MLC,'COMB')


1000	CONTINUE


	ENDSELECT

	RETURN
	END
C
C	=====================================================================
C	======================================================SONGSAK AUG2007
C	=====================================================================
      SUBROUTINE RESSPE2(PER,RES,ANGAP,NUMT,MAXA,ID,SPER,RVEC,RMAS,IDOF,
	1				   NSN,NSF,NEQ,NR,NWM,ICOMB,CVEC,SFREQ,BB,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP,OPER
C     -------------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN1,NEG,NBS,NLS,NLA,
     1              NSC1,NSF1,IDOF1(9),LCS,ISOLOP,LSYMM
	COMMON /MGRAV/ NGRAV
C     -------------------------------------------------------------------
	DIMENSION PER(NUMT),RES(NUMT,3),SPER(NR),SFREQ(NR),RESTR(3,NR)
	DIMENSION VH1(3),VH2(3),VG(3),ANGAP(3),SCOF(NEQ,3),SEF(NEQ,3)
	DIMENSION PARF(3,NR),SD(NR),SF(NEQ),RVEC(NEQ,NR),RMAS(NR,NR)
	DIMENSION ID(NSF,1),CVEC(NEQ,NR),IDOF(9),MAXA(1)
	DIMENSION BB(1)
	DIMENSION R(NEQ),SHAPET(NEQ,NR),C(NEQ)
C     -------------------------------------------------------------------
	NEAQ = 3   !THREE COMPONENT  H1,H2,VERTICAL

	 VG(1:3)  = 0.0D0
	VH1(1:3)  = ANGAP(1:3)
	 VG(NGRAV) = 1.0D0
	VH1(NGRAV) = 0.0D0
      CALL SCALEN (VH1,VH2,DUM,3)	
	VH1(1:3) = VH2(1:3)
      CALL VECPRD (VH1,VG,VH2)

	CALL CLEARMAT(SCOF,NEQ,NEAQ)
	DO J = 1,NSN
	DO II = 1,NEAQ
	IF(IDOF(II).EQ.II) THEN
	IEQ = ID(II,J)
	IF(IEQ.NE.0) THEN
	VH1(II) = 1D0
	VH2(II) = 1D0
	VG(II)  = 1D0
	SCOF(IEQ,1) = VH1(II)
	SCOF(IEQ,2) = VH2(II)
	SCOF(IEQ,3) =  VG(II)
	ENDIF
	ENDIF
	ENDDO
	ENDDO
	

	DO IR = 1,NR
	SP = DSQRT(SPER(IR))
	SFREQ(IR) = SP
	SP = (2.*3.14159)/SP
	DO IE = 1,NEAQ
	CALL RINTERPL(SP,RESTR(IE,IR),PER,RES(1,IE),NUMT)
      	RESTR(IE,IR) = RESTR(IE,IR)/SPER(IR)               !FOR ACCELERATION SPECTRUM
	ENDDO
	ENDDO


C	M*Ixyz MASS INFLUENCE MATRIX  SEF(NEQ,NEAQ)
	CALL CLEARMAT(SEF,NEQ,NEAQ)
	DO I=1,NEAQ
	  	CALL MAMULT(MAXA,BB,SCOF(1,I),SEF(1,I),TYP,'STD')
	ENDDO

C	
	CALL CLEARMAT(PARF,NEAQ,NR)
	CALL CLEARA(SD,NR)

	DO J=1,NR
C     MATRIX MULTIPLICATION C(IA,JB) = C + A(IA,JA)*B(JA,JB)
	CALL MATMULT (RVEC(1,J),SEF,PARF(1,J),1,NEQ,NEAQ,0.0,2)
	DO K=1,NEAQ
	PARF(K,J) = PARF(K,J)/RMAS(J,J)
	ENDDO
	CALL VECNORM(PARF(1,J),RESTR(1,J),SD(J),NEAQ)
	DO I=1,NEQ
  	CVEC(I,J) = RVEC(I,J)*SD(J)
	ENDDO
	ENDDO
	
C     VECTOR FORCE CONTROL
      IF (VH1(1).EQ.1.AND.VH1(2).EQ.0.AND.VH1(3).EQ.0)THEN ! VX DIRECTION
      NDO  = NEQ/6.0
      DO K = 1,3
        IF (K.EQ.1) NDIR = 3
        IF (K.EQ.2) NDIR = 4
        IF (K.EQ.3) NDIR = 5
        DO J = 1,NR
          DO M = 1,NDO
          CVEC(NDIR,J) = 0.0D0
          NDIR = NDIR + 6.0D0
          ENDDO
          M    = 0
          IF (K.EQ.1) NDIR = 3
          IF (K.EQ.2) NDIR = 4
          IF (K.EQ.3) NDIR = 5
        ENDDO
        J = 0
      ENDDO
      ELSEIF (VH1(1).EQ.0.AND.VH1(2).EQ.1.AND.VH1(3).EQ.0)THEN ! VY DIRECTION
      NDO  = NEQ/6.0
      DO K = 1,3
        IF (K.EQ.1) NDIR = 1
        IF (K.EQ.2) NDIR = 3
        IF (K.EQ.3) NDIR = 4
        DO J = 1,NR
          DO M = 1,NDO
          CVEC(NDIR,J) = 0.0D0
          NDIR = NDIR + 6.0D0
          ENDDO
          M    = 0
          IF (K.EQ.1) NDIR = 1
          IF (K.EQ.2) NDIR = 3
          IF (K.EQ.3) NDIR = 4
        ENDDO
        J = 0
      ENDDO
      ELSEIF (VH1(1).EQ.0.AND.VH1(2).EQ.1.AND.VH1(3).EQ.0)THEN ! VZ DIRECTION
            NDO  = NEQ/6.0
      DO K = 1,3
        IF (K.EQ.1) NDIR = 1
        IF (K.EQ.2) NDIR = 4
        IF (K.EQ.3) NDIR = 6
        DO J = 1,NR
          DO M = 1,NDO
          CVEC(NDIR,J) = 0.0D0
          NDIR = NDIR + 6.0D0
          ENDDO
          M    = 0
          IF (K.EQ.1) NDIR = 1
          IF (K.EQ.2) NDIR = 4
          IF (K.EQ.3) NDIR = 6
        ENDDO
        J = 0
      ENDDO
      ENDIF


	RETURN
	END
C
C	=====================================================================
C	======================================================SONGSAK AUG2007
C	=====================================================================

      SUBROUTINE RESPDATA(PER,RES,ANGAP,IEQK,NUMT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /EREQK/ IERQK(100,2),AERQK(10000),NUMEK,NELOC1,NELOC2
	DIMENSION PER(1),RES(1),ANGAP(3),EARTH(3)

	SELECTCASE(IND)

	CASE(0)
	IERQK(1,2) = 0	
	READ(ITI,*) 
	READ(ITI,*) NUMEK


	NELOC1 = 1
	NELOC2 = NELOC1 + 3*NUMEK

	DO I = 1,NUMEK
	NN = NELOC1 + 3*(I-1)
	READ(ITI,*) NUM,NUMT,NOPTION,ADEGREE,AERQK(NN+0),AERQK(NN+1),AERQK(NN+2)
	   
	   IF (NOPTION.EQ.2)THEN
	   CALL WIND_VECTOR (ADEGREE,EARTH)
	   AERQK(NN+0) = EARTH(1)
	   AERQK(NN+1) = EARTH(2)
	   AERQK(NN+2) = EARTH(3)
	   ENDIF

	IERQK(I,1)   = NUMT
	IERQK(I+1,2) = NUMT + IERQK(I,1)
	DO J = 1,NUMT
	MM1 = NELOC2 + 4*IERQK(I,2) + 0*NUMT + J
	MM2 = NELOC2 + 4*IERQK(I,2) + 1*NUMT + J
	MM3 = NELOC2 + 4*IERQK(I,2) + 2*NUMT + J
	MM4 = NELOC2 + 4*IERQK(I,2) + 3*NUMT + J
	READ(ITI,*) NUM,AERQK(MM1),AERQK(MM2),AERQK(MM3),AERQK(MM4)
	ENDDO
	ENDDO

	CASE(1)
	NUMT = IERQK(IEQK,1)
	NN   = NELOC1 + 3*(IEQK-1)
	ANGAP(1) = AERQK(NN+0)
	ANGAP(2) = AERQK(NN+1)
	ANGAP(3) = AERQK(NN+2)
	MM   = NELOC2 + 4*IERQK(IEQK,2)
	DO I = 1,NUMT
	MM1  = MM + 0*NUMT + I
	MM2  = MM + 1*NUMT + I
	MM3  = MM + 2*NUMT + I
	MM4  = MM + 3*NUMT + I
	PER(I)        = AERQK(MM1)
	RES(I+0*NUMT) = AERQK(MM2)
	RES(I+1*NUMT) = AERQK(MM3)
	RES(I+2*NUMT) = AERQK(MM4)
	ENDDO

	ENDSELECT


	RETURN
	END
C	=====================================================================
C	======================================================SONGSAK AUG2007
C	=====================================================================
	SUBROUTINE RINTERPL(XX,VV,VEC,DAT,NUMT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	PRODUCED BY SONGSAK JUN2007
C	LINEAR INTERPOLATION OF DATA
C	=============================================================
	DIMENSION VEC(1),DAT(1)

	NP = NUMT/4       !NUM 1ST PICK
	IF(NP.GT.NUMT) NP = NUMT
	I1 = 1
	I2 = NP

	DO I = 1,NUMT
	
	IF(I2.GT.NUMT) THEN
	I2 = NUMT
	I1 = I2 - NP + 1
	GOTO 5
	ENDIF

	TT1 = VEC(I1)
	TT2 = VEC(I2)

	IF(XX.GE.TT1.AND.XX.LE.TT2) GOTO 5

	I1 = I2
	I2 = I2 + NP

	ENDDO

	I1 = 1
	I2 = NUMT

5	CONTINUE


	DO I = I1,I2
	TT1 = VEC(I  )
	TT2 = VEC(I+1)

	IF(XX.GE.TT1.AND.XX.LE.TT2) THEN
	IT = I
	GOTO 10
	ENDIF 
	ENDDO

	VV = 0.0D0
	GOTO 20

10	CONTINUE

	DT = XX -TT1
	TT = TT2-TT1
	H1 = 1.0 - DT/TT
	H2 = DT/TT

	VV = DAT(IT)*H1 + DAT(IT+1)*H2

20	CONTINUE


	RETURN
	END

C	=====================================================================
C	======================================================SONGSAK AUG2007
C	=====================================================================
	SUBROUTINE NODALMS(RMASS,ID)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	PRODUCED BY SONGSAK JUN2007
C	ADDED NODAL MASS
	COMMON /NMAS/ NCM,NCC  !THIS COMMON IS ORIGINALLY USED BY BENNY
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
      COMMON /ADD_MASS/ ND_MASS(500),VALV_MASS(6,500)
C	=============================================================
	DIMENSION RMASS(1),ID(NSF,1),VALV(6)


	READ(ITI,*) 
	READ(ITI,*) NCM

	IF(NCM.LE.0) RETURN
	IF(NCM.GT.0) READ(ITI,*)

	CALL CLEARA(RMASS,NEQ)

	DO ICM = 1,NCM
	READ(ITI,*) ND,VALV(1:6)
      ND_MASS(ICM)       = ND         ! STORAGE MASS VALUE FOR COB AND COG
      VALV_MASS(1:6,ICM) =  VALV(1:6) ! STORAGE MASS VALUE FOR COB AND COG

	DO J = 1,6
	KK = 0
	DO K = 1,NSF
	IF(J.EQ.IDOF(K)) THEN
	KK = K
	EXIT
	ENDIF
	ENDDO
	IF(KK.NE.0) THEN
	IEQ = ID(KK,ND)
	IF(IEQ.NE.0) RMASS(IEQ) = RMASS(IEQ) + VALV(J)
	ENDIF
	ENDDO

	ENDDO
	

	RETURN
	END

C	=====================================================================
C	======================================================SONGSAK AUG2007
C	=====================================================================
	SUBROUTINE NODALDP(RDAMP,ID)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	PRODUCED BY SONGSAK JUN2007
C	ADDED NODAL DAMPER
	COMMON /NMAS/ NCM,NCC  !THIS COMMON IS ORIGINALLY USED BY BENNY
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
C	=============================================================
	DIMENSION RDAMP(1),ID(NSF,1),VALV(6)


	READ(ITI,*) 
	READ(ITI,*) NCC

	IF(NCC.LE.0) RETURN
	IF(NCC.GT.0) READ(ITI,*)

	CALL CLEARA(RDAMP,NEQ)

	DO ICM = 1,NCC
	READ(ITI,*) ND,VALV(1:6)

	DO J = 1,6
	KK = 0
	DO K = 1,NSF
	IF(J.EQ.IDOF(K)) THEN
	KK = K
	EXIT
	ENDIF
	ENDDO
	IF(KK.NE.0) THEN
	IEQ = ID(KK,ND)
	IF(IEQ.NE.0) RDAMP(IEQ) = RDAMP(IEQ) + VALV(J)
	ENDIF
	ENDDO

	ENDDO
	

	RETURN
	END

C	=====================================================================
C	======================================================SONGSAK AUG2007
C	=====================================================================
	SUBROUTINE RESCOME(NMOD,RFREQ,ICOMB)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /BRI3/ H4,ECC,ZET1,ZET2,RDM,RDK,NELW
C	------------------------------------------------------------------	
	DIMENSION RFREQ(NMOD),DRAT(NMOD),CMC(NMOD,NMOD)
C	-----------------------------------------------------------------
	CALL INTFILL('RESP',NFIL1,1,1,0) 
	CALL INTFILL('RESP',NFIL2,1,2,0) 

	FACTOR = 0.0D0
	CALL INIOPER(NFIL2,1,FACTOR,'CLEF','OLDF') 

C	Modal damping ratio array
	DRAT(1) = ZET1
	DO I=1,NMOD
		DRAT(I)=ZET2
	ENDDO
C
	GOTO (100,200),ICOMB
C	-----------------------------------
C-----1. MODAL COMBINATION BY CQC-METHOD
C	-----------------------------------		
  100 CONTINUE !COMBINE MODAL RESPONSES BY CQC METHOD'	
C	CALCULATE THE CROSS MODAL CORRELATION COEFICIENTS MATRIX
C	Form upper triangular matrix
  	CALL CLEARMAT(CMC,NMOD,NMOD)



	DO 150 I=1,NMOD
		DO 150 J=1,NMOD

			IF(I.EQ.J) THEN
				CMC(I,J) = 1.0
				GOTO 150
			ENDIF

			QIJ   = RFREQ(I)/RFREQ(J)
			DAMI = DRAT(I)
			DAMJ = DRAT(J)
			DIJ  = DAMI*DAMJ
			Q2 = QIJ*QIJ
			
			AA = 8.0*(SQRT(DIJ))*(DAMI + QIJ*DAMJ)*Q2*SQRT(QIJ)
			BB = ((1.0 - Q2)**2) + 4.0*DIJ*QIJ*(1.0 + Q2) +
     1			 4.0*(DAMI*DAMI + DAMJ*DAMJ)*Q2

			CMC(I,J) = AA/BB
			
  150	CONTINUE
C	----------------
C	Form full matrix
C	----------------
C	COMBINE FOR THE STRUCTURE RESPONSE OF ALL DOFs
	
			
	DO 175 I=1,NMOD 
		DO 175 J=1,NMOD
		FACTOR = ABS(CMC(I,J))
		CALL INIOPER(NFIL1,I,FACTOR,'CALL','OLDF') 
		FACTOR = 1.0D0
		CALL INIOPER(NFIL1,J,FACTOR,'MULT','OLDF') 
		FACTOR = 1.0D0
		CALL INIOPER(NFIL2,1,FACTOR,'COMB','OLDF') 
  175	CONTINUE 

	FACTOR = 1.0D0
	CALL INIOPER(NFIL2,1,FACTOR,'SQRT','OLDF') 

	GOTO 500
C	-----------------------------------
C-----2. MODAL COMBINATION BY SRSS-METHOD
C	-----------------------------------
  200	CONTINUE !'COMBINE MODAL RESPONSES BY SRSS METHOD'
		
	DO 210 I=1,NMOD
	FACTOR = 1.0D0
	CALL INIOPER(NFIL1,I,FACTOR,'CALL','OLDF') 
	FACTOR = 1.0D0
	CALL INIOPER(NFIL1,I,FACTOR,'MULT','OLDF') 
	FACTOR = 1.0D0
	CALL INIOPER(NFIL2,1,FACTOR,'COMB','OLDF') 
  210	CONTINUE			

	FACTOR = 1.0D0
	CALL INIOPER(NFIL2,1,FACTOR,'SQRT','OLDF') 


500	CONTINUE

	FACTOR = 1.0D0
	CALL INIOPER(NFIL2,1,FACTOR,'CALL','OLDF') 
	

      RETURN
      END


C	=====================================================================
C	======================================================SONGSAK SEP2007
C	=====================================================================
	SUBROUTINE SEIFNEW(TIM,PER,RES,NUMT,NEQ,NEAQ,SEF,SEFOR)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------------
	COMMON /LNDEQK/ LDEK,LPR2,LDESTP
	COMMON /LNDEQR/ EQAMP(3),EQGAP(3)
	COMMON /EREQK/ IERQK(100,2),AERQK(10000),NUMEK,NELOC1,NELOC2
C     -------------------------------------------------------------------
	DIMENSION PER(NUMT),RES(NUMT,3),SEF(NEQ,3),VEQ(3),SEFOR(1)
C     -------------------------------------------------------------------
	IF(LDEK.EQ.0) RETURN

C	==============================================
	VEQ(1:3) = 0.0D0
	IF(LDEK.EQ.1) THEN
	DO IE = 1,NEAQ
	CALL RINTERPL(TIM,VEQ(IE),PER,RES(1,IE),NUMT)
	ENDDO
	ENDIF
	IF(LDEK.EQ.2) THEN
	CALL LOADQUAK(TIM,2,LDESTP,AERQK(1),VV)
	VEQ(1:3) = VV*EQAMP(1:3)
	ENDIF
C	==============================================

	DO I = 1,NEQ
	SEFOR(I) = 0.0D0
	DO J = 1,NEAQ
	SEFOR(I) = SEFOR(I) + SEF(I,J)*VEQ(J)
	ENDDO
	ENDDO

	RETURN
	END
C
C	=====================================================================
C	======================================================SONGSAK SEP2007
C	=====================================================================
	SUBROUTINE SEIFVEP(TIM,PER,RES,NUMT,NEAQ,VEQ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------------
	COMMON /LNDEQK/ LDEK,LPR2,LDESTP
	COMMON /LNDEQR/ EQAMP(3),EQGAP(3)
	COMMON /EREQK/ IERQK(100,2),AERQK(10000),NUMEK,NELOC1,NELOC2
C     -------------------------------------------------------------------
	DIMENSION PER(NUMT),RES(NUMT,3),VEQ(3)
C     -------------------------------------------------------------------
	IF(LDEK.EQ.0) RETURN

C	==============================================
	VEQ(1:3) = 0.0D0
	IF(LDEK.EQ.1) THEN
	DO IE = 1,NEAQ
	CALL RINTERPL(TIM,VEQ(IE),PER,RES(1,IE),NUMT)
	ENDDO
	ENDIF
	IF(LDEK.EQ.2) THEN
	CALL LOADQUAK(TIM,2,LDESTP,AERQK(1),VV)
	VEQ(1:3) = VV*EQAMP(1:3)
	ENDIF
C	==============================================



	RETURN
	END
C
C	=====================================================================
C	======================================================SONGSAK SEP2007
C	=====================================================================
	SUBROUTINE SEIFVEC(MAXA,ID,IDOF,NSN,NSF,NEQ,NWM,SEF,
	1				   PER,RES,NUMT,NEAQ,BB,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C     -------------------------------------------------------------------
	COMMON /MGRAV/ NGRAV
	COMMON /LNDEQK/ LDEK,LPR2,LDESTP
	COMMON /LNDEQR/ EQAMP(3),EQGAP(3)
C     -------------------------------------------------------------------
	DIMENSION VH1(3),VH2(3),VG(3),ANGAP(3),SCOF(NEQ,3),SEF(NEQ,3)
	DIMENSION ID(NSF,1),IDOF(9),MAXA(1),PER(1),RES(1)
	DIMENSION BB(1)  !MASS MATRIC
C     -------------------------------------------------------------------

	IF(LDEK.EQ.1) CALL RESPDATA(PER,RES,ANGAP,1,NUMT,1)
	IF(LDEK.EQ.2) ANGAP(1:3) = EQGAP(1:3)

	NEAQ = 3   !THREE COMPONENT  H1,H2,VERTICAL

	 VG(1:3)  = 0.0D0
	VH1(1:3)  = ANGAP(1:3)
	 VG(NGRAV) = 1.0D0
	VH1(NGRAV) = 0.0D0
      CALL SCALEN (VH1,VH2,DUM,3)	
	VH1(1:3) = VH2(1:3)
      CALL VECPRD (VH1,VG,VH2)

	CALL CLEARMAT(SCOF,NEQ,NEAQ)
	DO J = 1,NSN
	DO II = 1,NEAQ
	IF(IDOF(II).EQ.II) THEN
	IEQ = ID(II,J)
	IF(IEQ.NE.0) THEN
	SCOF(IEQ,1) = VH1(II)
	SCOF(IEQ,2) = VH2(II)
	SCOF(IEQ,3) =  VG(II)
	ENDIF
	ENDIF
	ENDDO
	ENDDO

C	M*Ixyz MASS INFLUENCE MATRIX  SEF(NEQ,NEAQ)
	CALL CLEARMAT(SEF,NEQ,NEAQ)
	DO I=1,NEAQ
	CALL MAMULT(MAXA,BB,SCOF(1,I),SEF(1,I),TYP,'STD')
	ENDDO


	RETURN
	END
C
C	=====================================================================
C	======================================================SONGSAK SEP2007
C	=====================================================================

      SUBROUTINE READQUAK
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	NEW EARTHQUAKE ANALYSIS SONGSAK SEP2007
	COMMON /LNDEQK/ LDEK,LPR2,LDESTP
	COMMON /LNDEQR/ EQAMP(3),EQGAP(3)
	COMMON /EREQK/ IERQK(100,2),AERQK(10000),NUMEK,NELOC1,NELOC2

	DIMENSION PER(1),RES(1),ANGAP(3)

C	MOVE TO READ IN NLMODE
C	READ(ITI,*)
C	READ(ITI,*) LDEK   !0 = NOT INCLUDE GROUND MOTION, 1=TIME HISTORY ,2=DYNAMIC CURVE GROUND MOTION

	IF(LDEK.EQ.0) RETURN

	IF(LDEK.EQ.1) CALL RESPDATA(PER,RES,ANGAP,1,NUMT,0)

	IF(LDEK.EQ.2) THEN
	READ(ITI,*)
	READ(ITI,*) EQAMP(1:3),EQGAP(1:3)
	READ(ITI,*)
	CALL LOADQUAK (T,1,LDESTP,AERQK(1),PROPLD)
	ENDIF

	RETURN
	END

C	=====================================================================
C	======================================================SONGSAK SEP2007
C	=====================================================================

      SUBROUTINE LOADQUAK (T,J1,J2,A,PROPLD)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SUNIL 01/05/11 PROGRAMMED BY SUNIL
C	----------------------------------------------------------------
C	READ PROPORTIONAL LOAD CURVES FOR DYNAMIC ANALYSIS AND
C	COMPUTE PROPLD AT EACH TIME T; IT IS POSSIBLE TO GENERATE
C	ANY ARBITRAY LOAD USING THE FOURIER COEFFICIENTS
C	----------------------------------------------------------------
C	T			= CURRENT TIME
C	J1			= FLAG (1-READ LOAD CUREVES, 2-COMPUTE PROPLD)
C	J2			= NUMBER OF LOAD CURVES 
C	A(1,J2)		= STARTING TIMES OF EACH LOAD CURVE
C	A(2,J2)		= ENDING TIMES OF EACH LOAD CURVE
C	A(6,J2)		= COEFFICEIENTS OF FORMULAR
C				  A(3,J2)-A(5,J3) ARE IN FORCE UNITS
C				  A(6,J2) IS IN Rad/Sec. AND A(7,J2) IS IN Deg.
C	PROPLD		= PROPORTIONAL LOAD {A3+A4*T+A5(SIN(A6*T+A7))}
C	----------------------------------------------------------------
	CHARACTER*80 TITLE

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      DIMENSION A(7,1)

      GOTO (100,200), J1
C	-------------------------------
C	INPUT PROPORTIONAL LOADS CURVES
C	-------------------------------
100   CONTINUE
      DO 110 I=1,J2
		READ (ITI,*,END=102) (A(JK,I),JK=1,7)
110	CONTINUE
	RETURN
102	CALL ERRORS (9,I,J2,'LOAD CURVE')
C	-----------------------
C	COMPUTE VALUE AT TIME T
C	-----------------------
200	PROPLD = 0.0
      IF (T.LT.A(1,1) .OR. T.GT.A(2,J2)) GOTO 500
      I = 1
220	IF (T.GE.A(1,I) .AND. T.LT.A(2,I)) GOTO 230
      IF (I.GE.J2) GOTO 250
      I = I+1
      GOTO 220

230	CON = 3.141592654/180.0
	PROPLD = PROPLD + 
	1         A(3,I) + A(4,I)*T + A(5,I)*(DSIN(A(6,I)*T+CON*A(7,I)))
	I = I+1
	GOTO 220

250	RETURN 

500   RETURN

      END
C
C	=====================================================================
C	======================================================SONGSAK SEP2007
C	=====================================================================
	SUBROUTINE VFVMUL (MAXA,FM,V,R,NF,NR,MTYPE,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C     ----------------------------------------------------------------
C     PROGRAM TO TRIMULT MATRIX [R] = [V]T*[FM]*[V],[M] STORED IN ARRAY 
C	[V(NF,NR)] IS FULL FORM MATRIX: PRODUCT R STORED IN MATRIX
C     ----------------------------------------------------------------
C	INPUT ARGUMENTS
C	---------------
C	[FM(NF,NF)] = STIFF/MASS/DAMP MATRIX STORED IN ARRAY
C	[V(NF,NR)]  = FULL FORM MATRIX (EX. EIGENVECTOR MATRIX)
C	MAXA	    = DIAGONAL ADDRESS ARRAY OF FM
C	MTYPE		= R OPTION, = 1: FULL OUPUT MATRIX R, 
C						    = OTHER THAN 1: UPPER TRIAGULAR MATRIX	
C	---------------
C	OUPUT ARGUMENTS
C	---------------
C	R(NR,NR) = REDUCED MATRIX (PROJECTION OF [FM] ONTO [V])
C				   + IN FULL  SYMMETRIC  FORM IF MTYPE = 1	
C				   + IN UPPER TRIANGULAR FORM IF MTYPE = OTHER THAN 1
C	---------------
	DIMENSION FM(1),V(NF,1),R(NR,NR),RMJ(NF)
C	---------------------------------
C	FORM UPPER TRIANGULAR MATRIX OF R 
C	---------------------------------
	

	IJ = 0
	DO 25  J=1,NR

		CALL MAMULT(MAXA,FM,V(1,J),RMJ,TYP,'STD')
		DO 20  I=J,NR
		     BRT = 0.0
	         DO 10  K=1,NF
   10				BRT = BRT + V(K,I)*RMJ(K)
			 IJ  = IJ+1
   20		     R(I,J) = BRT

   25 CONTINUE

C	-----------------------------------------------------------	
C	FORM THE FULL SYMMETRIC MATRIX FROM UPPER TRIANGULAR MATRIX
C	-----------------------------------------------------------	
	IF(MTYPE.EQ.1) GOTO 30
	GO TO 45
   30	DO 40 I=1,NR
		DO 35 J=1,NR
		   IF (I.LT.J) GO TO 40
	       R(J,I) = R(I,J)
   35		CONTINUE
   40 CONTINUE

   45	RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE VFVMULD (MAXA,FM,V,DIAR,NF,NR,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C     ----------------------------------------------------------------
C     PROGRAM TO TRIMULT MATRIX [R] = [V]T*[FM]*[V],[M] STORED IN ARRAY 
C	[V(NF,NR)] IS FULL FORM MATRIX: PRODUCT R STORED IN MATRIX
C     ----------------------------------------------------------------
C	INPUT ARGUMENTS
C	---------------
C	[FM(NF,NF)] = STIFF/MASS/DAMP MATRIX STORED IN ARRAY
C	[V(NF,NR)]  = FULL FORM MATRIX (EX. EIGENVECTOR MATRIX)
C	MAXA	    = DIAGONAL ADDRESS ARRAY OF FM
C	---------------
C	OUPUT ARGUMENTS
C	---------------
C	DIAR(I) = R(I,I)   VECTOR WITH DIAGONAL ELEMENTS OF R 
C	---------------
	DIMENSION FM(1),V(NF,1),R(NR,NR),DIAR(NR),MAXA(1) 


C	------------------------------
C	FORM UPPER TRIANGULAR R MATRIX
C	------------------------------
	CALL VFVMUL (MAXA,FM,V,R,NF,NR,0,TYP)

C	---------------------------------------
C	EXTRACT DIAGONAL ELEMENTS FORM R MATRIX
C	---------------------------------------
   	DO 50 I=1,NR
	DIAR(I) = R(I,I)
   50 CONTINUE


	RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE FREQ_STEP_INCR (INC,CTIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND
      
      IF (ICONTROLSPEC.NE.1) RETURN
      
      IF (OPRES.EQ.1)THEN
          IF (INC.EQ.1)THEN
              IF (STEPSTAT.EQ.0.0D0) CTIM = 0D0
              IF (STEPSTAT.NE.0.0D0) CTIM = STEPSTAT
          ELSEIF (INC.NE.1)THEN
              CTIM = CTIM + STEPINCR
          ENDIF
      ENDIF    

 
      
      END
C	=====================================================================
C	=====================================================================
      SUBROUTINE COMBINE_OFF_GEN (R,R_OFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION R(NEQ),R_OFF(NEQ)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      ! GENERAL LOADCASE + OFFSHORE LOADCASE
      DO I = 1,NEQ
      R(I) = R(I)+R_OFF(I)
      ENDDO
      END