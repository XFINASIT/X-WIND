C	=====================================================================
C	=====================================================================
C	=====================================================================
      PROGRAM  FINAS       
	USE IFPORT  
C      USE MKL_PARDISO_PRIVATE
C      USE MKL_DSS
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      INCLUDE 'mkl.fi'
      !OGICAL(4) MAKEDIR
      CHARACTER*200   PATH
      CHARACTER*8 DATE_CLOCK


C     ----------------------------------------------------------------
C
C                    * * * *   X F I N A S   * * * *
C
C	AN INTERACTIVE FINITE ELEMENT STRUCTURAL ANALYSIS PROGRAM WITH 
C	            LINEAR AND NONLINEAR DYNAMIC OPTIONS
C            By  Imperial College,London & Dr. K.D. Kim in AIT 
C     ----------------------------------------------------------------
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI
      COMMON A(9000000),IA(9000000)
      
	COMMON /WLAS/ MERW,MEIW,NLASI,NLASW
	! ------------ COMMON FOR OFFSHORE ------------------------------
	COMMON /WARNING/ WARNING,RAMDA(100),RK
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION
	! ---------------------------------------------------------------
	
	COMMON /MEMW/ W(7000000),IW(7000000)
      
      COMMON /COUANA/ COUPLEANALYSIS
      

	MEMA = 9000000
	MEMI = 9000000
	MERW = 7000000
	MEIW = 7000000
      
      Call OPENXSEADATA (1)
      
      
      ! THE COUPLED ANALYSIS WAS WORK COMPATIBILITY WITH X-OFFSHORE + FAST BY ACCESS THE FILE NAME "NCOUPLE.DAT"
      ! COUPLEANALYSIS = 1 ( UNCOUPLED ANALYSIS )
      ! COUPLEANALYSIS = 2 ( SEMI-COUPLED ANALYSIS )
      ! COUPLEANALYSIS = 3 ( COUPLED ANALYSIS )
      READ (551,*,IOSTAT=IREASON) COUPLEANALYSIS
      ! CLEAR DIRECTORY OF COUPLED ANALYSIS. THIS IS USED FOR PREVENT READING OLD DATA
      IF (COUPLEANALYSIS.EQ.2D0.OR.COUPLEANALYSIS.EQ.3D0) CALL DIRECTORY_OF_FILE_FAST_CLEAR
      
C	CONSTRUCTION STAGE ANALYSIS SONGSAK NOV2007
      CALL INTCONS

C	MEMORY MANAGER INITIALIZE
      CALL INTMEMM

C	BLOCK SOLVER SONGSAK JAN2007
	CALL ADMEMY

C	SECURITY CHECK SONGSAK MAR2006
	CALL SECURE
	
C     INITIALIZE MEMORY INDEX FOR W,IW ARRAY
      CALL LOCAW(5,0) 
      
C     --------------------------------
C     INPUT AND OUTPUT OF SYSTEM DATA
C     --------------------------------
      CALL INPUT1
      
      CALL SYSOUT_MEM_USAGE
      
CC      CALL MWCOUNT
      
 !    CALL McF_DIFFRECTION
 !    CALL WAVE_DIFFRECTION_XFINAS_FORCE
              
      IF(MODEX.EQ.1)GOTO 900
C     -----------------------------------
C     LINEAR OR NONLINEAR STATIC ANALYSIS
C     -----------------------------------
      CALL ITERAT1(W,0)

C     CLEAR ALL MEMORY IN PARDISO
      CALL CLRPARD
 
C     OFFSHORE WARNING SIGN
C      CALL OFFSHORE_WARN
 
 1200 CONTINUE
      
CC      CALL MWCOUNT
      
	WRITE(*,*) '-----------END OF XFINAS CALCULATION-----------'
	CALL CPU_TIME (TOEY1)
C     ---------------------------------------------------------------
C     PRINTS OUT SYSTEM DATA, SOLUTION TIMES AND STORAGE REQUIREMENTS
C     ---------------------------------------------------------------
 900  CALL SYSOUT

      
      ! DELETE FILE
      CALL CLOSE_FILE
C	---------------------------------------------------------
C     CLOSE OUTPUT FILE FOR EXCEL REPORT
!      CLOSE(50)
!      CLOSE(51)
!      CLOSE(52)
!      CLOSE(53)
!      CLOSE(54)
!      CLOSE(55)
!      CLOSE(106)
!      CLOSE(107)
!      CLOSE(108)
!      CLOSE(109)
C	---------------------------------------------------------
!      CLOSE(1000)
C
!      STOP
100   END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SYSOUT_MEM_USAGE
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     PRINTS OUT SYSTEM DATA, SOLUTION TIMES AND STORAGE REQUIREMENTS
C     ---------------------------------------------------------------

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /BOPT/ IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      
C	Next common /WLAS/ created Jul04 to store the last pointers in IW and W-arrays
	COMMON /WLAS/ MERW,MEIW,NLASI,NLASW      
      
C     COMMON FOR IN-CORE ARRAY STORAGE IN MemManager.FOR      
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(100000)
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1000000)
      
C     COMMON FOR IN-CORE ARRAY STORAGE IN INOUTUTIL.FOR  
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(2000000)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(6000000)
      
      MTMMI = MEMMI(1)
	NTMMI = MEMMI(4)
      MTMMR = MEMMR(1)
	NTMMR = MEMMR(4)
      
      MTTMI = MTOTI(1)
	NTTMI = MTOTI(4)
      MTTMR = MTOTR(1)
	NTTMR = MTOTR(4)
      
C     -----------
C     SYSTEM DATA
C     -----------
      WRITE (ITO,200)
      WRITE (ITO,210)  NEQ,NWK,NWM,MK,BM,
     1                 NELEA,NELEI,MEMA,LASTA,MEMI,LASTI,MERW,NLASW,MEIW,NLASI,
     2                 MTMMI,NTMMI,MTMMR,NTMMR,MTTMI,NTTMI,MTTMR,NTTMR
      
      WRITE (10,200)
      WRITE (10,210)   NEQ,NWK,NWM,MK,BM,
     1                 NELEA,NELEI,MEMA,LASTA,MEMI,LASTI,MERW,NLASW,MEIW,NLASI,
     2                 MTMMI,NTMMI,MTMMR,NTMMR,MTTMI,NTTMI,MTTMR,NTTMR

 200  FORMAT  (1H1///,29X,21(1H*)/29X,1H*,19X,1H*/
     1         29X,21H* TOTAL SYSTEM DATA */
     2         29X,1H*,19X,1H*/29X,21(1H*)//)
 210  FORMAT (14X,'NUMBER OF EQUATIONS . . . . . . . . NDOFS =',I8/  
     +        14X,'FINAL NO. OF COEF.IN MATRIX A .   NWK =',I12/    
     +        14X,'FINAL NO. OF COEF.IN MATRIX B .   NWM =',I12/    
     +        14X,'FINAL MAXIMUM SEMI-BANDWIDTH  . . .    MK =',I8/    
     +        14X,'FINAL MEAN SEMI-BANDWIDTH . . . . .    BM =',F8.2//
     +        14X,'MAX. LENGTH OF ELEMENT INFO IN A. . NELEA =',I8/    
     +        14X,'MAX. LENGTH OF ELEMENT INFO IN IA . NELEI =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK A .  MEMA =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK A . LASTA =',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK IA.  MEMI =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK IA. LASTI =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK W .  MERW =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK W . NLASW =',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK IW.  MEIW =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK IW. NLASI =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /MEMTINT/ =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /MEMTINT/ =',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /MEMTREL/ =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /MEMTREL/ =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /CONSTINT/=',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /CONSTINT/=',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /CONSTREL/=',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /CONSTREL/=',I8//)

      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SYSOUT
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     PRINTS OUT SYSTEM DATA, SOLUTION TIMES AND STORAGE REQUIREMENTS
C     ---------------------------------------------------------------

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /BOPT/ IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      
C	Next common /WLAS/ created Jul04 to store the last pointers in IW and W-arrays
	COMMON /WLAS/ MERW,MEIW,NLASI,NLASW   
      
C     COMMON FOR IN-CORE ARRAY STORAGE IN MemManager.FOR      
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(100000)
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1000000)
      
C     COMMON FOR IN-CORE ARRAY STORAGE IN INOUTUTIL.FOR  
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(2000000)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(6000000)
      
      MTMMI = MEMMI(1)
	NTMMI = MEMMI(4)
      MTMMR = MEMMR(1)
	NTMMR = MEMMR(4)
      
      MTTMI = MTOTI(1)
	NTTMI = MTOTI(4)
      MTTMR = MTOTR(1)
	NTTMR = MTOTR(4)    
      
C     -----------
C     SYSTEM DATA
C     -----------
C      CALL MEMORY (MEMO) !NON-EXIST SUBROUTINE
      WRITE (ISO,200)
      WRITE (ISO,210)  NEQ,IPF1,IBW1,NWK,NWM,MK,BM,
     1                 NELEA,NELEI,MEMA,LASTA,MEMI,LASTI,MERW,NLASW,MEIW,NLASI,
     2                 MTMMI,NTMMI,MTMMR,NTMMR,MTTMI,NTTMI,MTTMR,NTTMR

C     --------------
C     SOLUTION TIMES
C     --------------
C	SUNIL 12/01/01 CHANGE TIME ROUTINE
C      TIM(11) = TIM(11)-TIM(13)-TIM(12)

	CALL CPU_TIME (TIM(20))

      WRITE (ISO,300)
      WRITE (ISO,310) (TIM(I), I=1,6)
      WRITE (ISO,320) NSTEP,LCS,TIM(11),TIM(10),TIM(13),TIM(14),
     +                TIM(15),TIM(16),TIM(17)
      WRITE (ISO,330) NEIG,TIM(18),TIM(20)

      
      WRITE (ITO,300)
      WRITE (ITO,310) (TIM(I), I=1,6)
      WRITE (ITO,320) NSTEP,LCS,TIM(11),TIM(10),TIM(13),TIM(14),
     +                TIM(15),TIM(16),TIM(17)
      WRITE (ITO,330) NEIG,TIM(18),TIM(20)
      
      WRITE (10,300)
      WRITE (10,310) (TIM(I), I=1,6)
      WRITE (10,320) NSTEP,LCS,TIM(11),TIM(10),TIM(13),TIM(14),
     +                TIM(15),TIM(16),TIM(17)
      WRITE (10,330) NEIG,TIM(18),TIM(20)
C     ----------------------------
C     FLUSH BUFFERS
C     RELEASE DYNAMIC STORAGE AREA
C     ----------------------------
      REWIND ISO
      REWIND NDATI

 900	CLOSE(UNIT=100)

      DO 910  IEG=1,NEG
      NELEMI = 10 + IEG
      NELEMA = 30 + IEG
	CLOSE(UNIT=NELEMI,STATUS='DELETE')
	CLOSE(UNIT=NELEMA,STATUS='DELETE')
 910  CONTINUE

	IF (ISOLOP.EQ.11) THEN
	CLOSE(UNIT=500)
	ELSE
	CLOSE(UNIT=500,STATUS='DELETE')
	ENDIF
C
 200  FORMAT  (1H1///,29X,21(1H*)/29X,1H*,19X,1H*/
     1         29X,21H* TOTAL SYSTEM DATA */
     2         29X,1H*,19X,1H*/29X,21(1H*)//)
 210  FORMAT (14X,'NUMBER OF EQUATIONS . . . . . . . . NDOFS =',I8/
     +        14X,'ORIGINAL NO.OF COEF.IN MATRIX A . .  IPF1 =',I8/    
     +        14X,'ORIGINAL MAXIMUM SEMI-BANDWIDTH . .  IBW1 =',I8/    
     +        14X,'FINAL NO. OF COEF.IN MATRIX A .   NWK =',I12/    
     +        14X,'FINAL NO. OF COEF.IN MATRIX B .   NWM =',I12/    
     +        14X,'FINAL MAXIMUM SEMI-BANDWIDTH  . . .    MK =',I8/    
     +        14X,'FINAL MEAN SEMI-BANDWIDTH . . . . .    BM =',F8.2//
     +        14X,'MAX. LENGTH OF ELEMENT INFO IN A. . NELEA =',I8/    
     +        14X,'MAX. LENGTH OF ELEMENT INFO IN IA . NELEI =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK A .  MEMA =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK A . LASTA =',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK IA.  MEMI =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK IA. LASTI =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK W .  MERW =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK W . NLASW =',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK IW.  MEIW =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK IW. NLASI =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /MEMTINT/ =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /MEMTINT/ =',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /MEMTREL/ =',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /MEMTREL/ =',I8//    
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /CONSTINT/=',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /CONSTINT/=',I8/
     +        14X,'AVAILABLE SPACE IN COMMON BLOCK /CONSTREL/=',I8/    
     +        14X,'SPACED USED     IN COMMON BLOCK /CONSTREL/=',I8)
 300  FORMAT  (///,31X,18(1H*)/31X,1H*,16X,1H*/
     1         31X,18H* SOLUTION TIMES */
     2         31X,1H*,16X,1H*,/31X,18(1H*)//)
 310  FORMAT (26X,28HINPUT PHASE (SEGMENT INPUT1)/26X,28(1H-)/
     +        14X,'CONTROL INFORMATION AND MESH GENERATION   =',F9.2/
     +        14X,'GENERATION OF ELEMENT INFORMATION. . . .  =',F9.2/
     +        14X,'BANDWIDTH OPTIMISATION . . . . . . . . .  =',F9.2/
     +        14X,'STIFFNESS ADDRESSES AND BLOCK SIZES  . .  =',F9.2/
     +        14X,'INPUT AND GENERATION OF LOADING DATA . .  =',F9.2/
     +        14X,'TOTAL SOLUTION TIME (INPUT PHASE)  . . .  =',F9.2)
 320  FORMAT (//21X,23HSTEP-BY-STEP SOLUTION (,I2,12H LOAD STEPS)/
     +     21X,37(1H-)/
     +     14X,'NUMBER OF LOAD CASES . . . . . . . . . .  =',I9//
     +     14X,'ELEMENT CALCULATION OF STRES.&STIF . . .  =',F9.2/
     +     14X,'ASSEMBLE OF EFFECTIVE LOAD VECTOR  . . .  =',F9.2/
     +     14X,'ASSEMBLY OF COMPACTED GLOBAL STIFFNESS .  =',F9.2/
     +     14X,'    FACTORIZATIONS OF STIFFNESS BLOCKS .  =',F9.2/
     +     14X,'    REDUCTIONS,BACK-SUBSTITUTIONS OF LOADS=',F9.2/
     +     14X,'    EQUILIBRIUM ITERATIONS. . . . . . . . =',F9.2/
     +     14X,'PRINT OUT OF DISPLACEMENTS AND STRESSES   =',F9.2/)
 330  FORMAT(/22X,36HEIGENVALUE ANALYSIS (SEGMENT STABIL)/22X,36(1H-)/
     +     14X,I3,' SOLUTIONS OF EIGENVALUE PROBLEM  . .  =',F9.2//
     +     14X,'TOTAL SOLUTION TIME. . . . . . . . . . .  =',F9.2)

      !STOP
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE LOCATI (IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     COMPUTES THE FIRST ADDRESSES OF ALL THE FLEXIBLE ARRAYS IN A
C	------------------------------------------------------------
C     STORES THESE ADDRESSES IN COMMON BLOCK /LOCA/
C	---------------------------------------------
C     ARRAY A     ADDRESS  CONTENTS
C	-------     -------  --------
C	DIRCOS(9,NLA)   LDC  DIRECTIONS COSINES OF LOCAL AXES
C     XY(NSN,NSC)     LXY  NODAL CO-ORDINATES X,Y,Z
C     PROPM(NMP,NMPS) LMP  MAT. PROPERTIES FOR CURRENT ELEMENT GROUP
C     PROPG(NGP,NGPS) LGP  GEO. PROPERTIES FOR CURRENT ELEMENT GROUP
C     PROPO(  6,NOPS) LOP  OFF. PROPERTIES FOR CURRENT ELEMENT GROUP
C	COORD(NCO*NNM,NELE)
C					LCO  NODAL COORDINATE FOR CURRENT ELEMENT GROUP
C     S(NWS)          LES  ELEMENT STIFFNESS MATRIX
C	COOR(NCO,NNM)	LEC  COORINATES FOR SINGLE ELEMENT
C     EDIS(NEF)       LED  NODAL DISPLACEMENTS OF ELEMENT
C     EDISI(NEF)      LEI  NODAL DISPLACEMENT INCREMENTS OF ELEMENT
C     ELOD(NEF)       LEE  ELEMENT EQUILIBRIUM LOADS
C     RFAC(NSTEP)     LLF  LOAD FACTORS TO DEFINE THE CURRENT LOAD V.
C     R(NEQ)          LLV  1. COLUMN HEIGHTS IN COMP.STIFFNESS MATRIX
C                          2. CURRENT LOAD VECTOR
C     RE(NEQ)         LRE  1. CURRENT LOAD VECTOR
C                          2. RESIDUAL LOAD VECTOR
C                          3. ADDITIONAL DISPLACEMENT INCREMENT
C     DISPI(NEQ)      LDI  SUM OF DISPL. INCREMENTS FOR CURRENT STEP
C     DISLI(NEQ)      LDL  DISPLACEMENT INCREMENT FOR LAST ITERATION
C     DISP(NEQ)       LDT  VECTOR OF TOTAL DISPLACEMENTS
C     VELO(NEQ)       LVL  VECTOR OF TOTAL VELOCITIES (DYNAMIC)
C     ACCE(NEQ)       LAL  VECTOR OF TOTAL ACCELETATIONS (DYNAMIC)
C     EFLO(NEQ)       LEF  EFFECTIVE LOAND VECTOR
C     DUMY(NEQ)       LDU  DUMMY VETOR
C     PRLO(7,NSTEP)	LPR  PROPROTIONAL LOAD TABLE (DYNAMIC)
C     D(NEQ)          LDK  DIAGONAL TERMS OF GLOBAL STIFFNESS MATRIX
C	AMV(NMV,3)      LMV  NORMALIZED MATERIAL REFERENCE AXIS

C     ARRAY IA    ADDRESS  CONTENTS
C	--------    -------  --------
C     ID(NSF,NSN)     LID  BOUNDARY CODES AND EQUATION NUMBERS
C	IDSET(NSN)      LDS  LOCAL AXES SET NUMBERS (IF NLA.GT.0)
C     LEST(NEG*2)     LEL  STORAGE LENGTHS OF ARRAYS IA & A 
C						 OF ELEMENT GROUP DATA INFORMATION
C	NCON(NSN,MCON+7)LCH  NODAL CONNECTIVITY TABLE (MCON.LE.MAXCON.LE.51)
C	IA(40)          LNU  STORAGE FOR CURRENT ELEMENT GROUP INFORMATION
C     MTSET(NELE)     LMS  MAT. PROPERTY SET NUMBERS FOR CURRENT ELENENT GROUP
C	IGSET(NELE)     LGS  GEO. PROPERTY SET NUMBERS FOR CURRENT ELENENT GROUP
C     IOSET(NELE)     LOS  OFF. PROPERTY SET NUMBERS FOR CURRENT ELENENT GROUP
C	NODEX(NEX,NELE) LEX  EXCESS NODE NUMBERS FOR CURRENT ELEMENT GROUP
C     LM(NEF,NELE)    LLM  ELEMENT CONNECTIVITY FOR CURRENT ELEMENT GROUP
C						 & GLOBAL EQUATION NUMBER FOR ELEMENT DOF
C     LC(NEF,NELE)    LCN  ELEMENT CONNECTIVITY 
C     MAXA(NEQ1)      LMA  ADDRESSES OF DIAGONAL ELEMENTS IN STIFFNESS
C     IAX(NELE)       LXI  ELEMENT MATERIAL AXIS NUMBER

C     NELEA                MAX. LENGTH OF ELEMENT GROUP DATA IN ARRAY A
C     NELEI                MAX. LENGTH OF ELEMENT GROUP DATA IN ARRAY IA
C	MEMA				 AVA. MEMEROY IN ARRAY A
C	MEMI				 AVA. MEMEROR IN ARRAY IA
C	LASTA				 LAST ADDRESS USED IN ARRAY A
C	LASTI                LAST ADDRESS USED IN ARRAY IA
C	ISTOR				 MAX. LENGTH OF STIFFNESS BLOCK

C     RT(NK,NE)       LTL  CONSERVATIVE GLOBAL TRACTIONS    (NTL>0)
C     RP(NE)          LPL  NON CONSERVATIVE NORMAL PRESSURE (NPL>0)
C
C     MAXFL                MAX. FIELD LENGTH AVAILABLE
C     --------------------------------------------------------------
	CHARACTER*10 MESSAGE
      LOGICAL      PROMPT,ERROR
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5              LCN,LDIM,LFRE,LSFC,LLOF
C     ----------------------------------------------
C	FOLLOWING COMMON BLOCKS ADDED BY DESILVA- 2003OCT
      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN
      COMMON /OFFS/ NOPS,NHIGE
C	ADDED BY DE SILVA 1/05/2004 SEEPAGE ANALYSIS
	COMMON /SPBC/ NSS,NLSS

C	NSS		= NUMBER OF SPRING BOUNDARY COND
C	NLSS	= NUMBER OF LOCAL SPRING BOUNDARY COND
	COMMON /GiDEle/ LGID

C----	Next Common Block added 13Nov03 by NguyenDV to store the spring 
C	damping factors(LSC)
      COMMON /LOCD/ LSC

	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /LOGO/ PROMPT,ERROR,ITEST
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /BOPT/ IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU

      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET

      COMMON /DYNA/ CDEN,IMASS
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT

C	NEXT COMMON ADDED BY GILSON - MARCH2004 (GRID ANALYSIS)
	COMMON /ILAN/ ILANE,NLN(12),IELN(12),IELL(750,12),LSPANS(12,15)
     +             ,NLSPAN(12),IGRD
      COMMON A(9000000),IA(9000000)

C     -----------------------------------------------------
C	REACTION BY SONGSAK APR2006
	COMMON /REACT/ LRC,LRCT,MFQ,LRID
C     -----------------------------------------------------
C	==================================================================
C	CABLE PRETENSION OPTIMIZATION
	COMMON /CB556/ LCPZ,NCABZ,KBOPZ,NCOBZ
C	==================================================================
C	CABLE PRETENSION OPTMIZATION SONGSAK NOV2006 STORE GID ELEMENT NUMBER
	COMMON /CB558/ LGOPZ
C	==================================================================
C	JOINT GLOBAL CONSTRAINT SONGSAK JAN2007 
	COMMON /LC007/ LC107,LC207,LC307,LC407
C	==================================================================
C	COMMON BLOCK FOR HEAT SONGSAK MAR2007
	COMMON /SHEAT/ NHAEF,LHEAT1,LHEAT2,LHEAT3,LHEAT4
C	==================================================================
	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE
C	==================================================================
C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN
C	==================================================================

      GOTO (100,200,300,400,500,600,700,800,900),IND

C     -----------------------------------------------------
C     FIND ALLOCATED CORE SIZE (AS GIVEN ON JOB CARD)
C     COMPUTE LENGTH OF DYNAMIC AREA
C     ASSIGNE APPROPRIATE SIZE TO FAST DYNAMIC LOADER
C     REQUEST FIXED POSITION BLOCK
C     FIND AVAILABLE SPACE MEMA IN BLANK COMMON BLOCK A(1)
C     -----------------------------------------------------
C     ---------------------------------------------------------------
C     EQUATION NUMBERS, STRUCTURE COORDINATES AND CONNECTIVITY MATRIX
C     ---------------------------------------------------------------
100	MESSAGE = 'COORDINATE'

C	IA-ARRAY POINTER
	LID   = 1

	LRID  = LID  + NSF*NSN       !REACTION SONGSAK APR2006

	LC307 = LRID  + NSF*NSN         !JOINT CONSTRAINT SONGSAK JAN2007
	LC407 = LC307 + 11*LC107        !JOINT CONSTRAINT SONGSAK JAN2007
	LDS   = LC407 + LC107*LC207     !JOINT CONSTRAINT SONGSAK JAN2007

      LEL   = LDS + NSN
	LCH   = LEL + NEG*2

	LASTI = LCH + NSN*NSF
	IF (MAXCON.GT.NSF) LASTI = LCH + NSN*(MAXCON+7) + 1

C	A-ARRAY POINTER
      LDC = 1
      LXY = LDC + 9*NLA
      IF (NLA.LE.0) LXY = LDC

	LHEAT1= LXY   + NSN*NSC     !FOR HEAT AMBIENT TEMP SONGSAK MAR2007
	LSS   = LHEAT1+ 2*LCS 
C	-----------------------------------------------------
C	Changed by Songsak Nov2006 Bilinear spring support Add LSS2,LSS3
C	LSS2 = INPUT DATA FOR BILINEAR SPRING SUPPORT
C	LSS3 = WORKING ARRAY FOR BILINEAR SPRING SUPPORT
	LSS2 = LSS  +   6*NSN
	LSS3 = LSS2 + 7*6*(NSS+NLSS)   
	LDIM = LSS3 + 6*5*(NSS+NLSS)
	A(LSS2:LDIM-1) = 0.0D0
C	-----------------------------------------------------
      LFRE = LDIM + NROOT   !LDIM STORE DIAGONAL OF MODAL MASS AND LFRE STORE MODAL FREQUENCY
      LSC  = LFRE + NROOT
      
	LASTA= LSC + NSN

      MESSAGE = 'CONNECTION'

      RETURN


C     ----------------------------------------------
C     ELEMENT GROUP DATA
C     ----------------------------------------------
 200  MESSAGE = 'ELEM. DATA'

C	IA-ARRAY POINTER
C	----------------
	LNU = LASTI
	LMS = LNU + 100

C	-----------------------------------------------------------------
	LGID = LMS + NELE
C	-----------------------------------------------------------------

	LGS  = LGID + NELE
	IF(ITYPE.EQ.5) THEN
		LOS = LGS + NELE
		LHG = LOS + NELE
		LHGN= LHG + 14*NHIGE
		LEX = LHGN+ NELE
	ELSEIF(ITYPE.EQ.17) THEN	!TENDON APR2009
		LOS  = LGS + NELE
		LHG  = LOS
		LHGN = LHG
		LTDN = LHGN
		LEX  = LTDN + NTEND*10   !TENDON DATA
	ELSE
		LOS = LGS + NELE
		LHG = LOS
		LHGN= LHG
		LEX = LHGN
	ENDIF

      LSFC= LEX + NEX*NELE
      IF(ITYPE.EQ.9.AND.ISTYP.EQ.10) THEN !SHELL 3 NODE WITH NO ROTATION BY ONATE (SUPPORT ASSIGNED IN ELEMENT FACE)
        LLM = LSFC+ 2*NEFC*NELE !NEFC = NUMBER OF ELEMENT FACE = NEFC  (SEE ELEMIN)
      ELSE
        LLM = LSFC 
      ENDIF
      
C     LLM = ACTIVE EQUATION      
	LRC = LLM + NEF*NELE   !REACTION EQUATION
	LCN = LRC + NEF*NELE   !ELEMENT CONNECTIVITIES

	LXI = LCN + NEF*NELE

	LMA = LXI + NELE

C	A-ARRAY POINTER
C	---------------
990	LMP = LASTA

	LGP = LMP+NMP*NMPS

	IF(ITYPE.EQ.5) THEN
		LOP = LGP + NGP*NGPS     
		LCO = LOP + 6*NOPS
	ELSE
		LOP = LGP + NGP*NGPS
		LCO = LOP
	END IF


	LHEAT2 = LCO    + NCO*NNM*NELE         !ADDED FOR HEAT SONGSAK MAR2007
	LHEAT3 = LHEAT2 + NELE*NPT*LCS         !ADDED FOR HEAT SONGSAK MAR2007
	LHEAT4 = LHEAT3 + NELE*NHAEF*LCS       !ADDED FOR HEAT SONGSAK MAR2007
	LES    = LHEAT4 + NELE*NHAEF*LCS       
C     ----------------------------------------------
      LEC = LES + NWS
      
C     TOEY WORK
C     LTO = LEC + NWS
C     ----------      

	LMV = LEC + NCO*NNM
	LED = LMV + 3*NMV
      LEI = LED + NEF
      LEE = LEI + NEF

C	IF (NLOPT.EQ.0) LEI = LED
	LLF = LEE + NEF

      IF (IFPR(4).GT.0) GOTO 910

      RETURN
C     ----------------------
C     BANDWIDTH OPTIMISATION
C     ----------------------
 300  MESSAGE = 'BAND OPTIM'
C	LCON FOR BAND OPTIMITIATION SUBROUTINE
      LCON = NSN

C     -----------------------------------
C     SET LAST POSITION OF ARRAY IA AND A
C     -----------------------------------
	LASTA = LASTA + NELEA
	LASTI = LASTI + NELEI

      RETURN
C     -------------------------
C     DIAGONAL ADDRESSES IN [K]
C     -------------------------
 400  MESSAGE = 'ADDRESSES '
C	IA-ARRAY POINTER
	LMA   = LASTI
      LASTI = LMA + NEQ1

      RETURN
C     ----------------------------------------------
C     NUMBER OF COLUMN, COUPLING OF STIFFNESS BLOCKS
C     ----------------------------------------------
 500  MESSAGE = 'BLOCK NUM.'

      MESSAGE = 'STIF.BLOCK'

      RETURN
 600  MESSAGE = 'BLOCK INF.'

      RETURN
C	----------------------------
C     LOAD FACTORS AND LOAD VECTOR
C     ----------------------------
 700  NBLOCK = 1
	ISTOR  = NWK
C	-------------------------------------------
C	IA-ARRAY POINTER
	LGOPZ = LASTI
        
	LASTI = LGOPZ+NCABZ		!LGOPZ ADDED BY SONGSAK NOV2006 FO OPTIMIZATION
C	-------------------------------------------

C	A-ARRAY POINTER
      LLF = LASTA
C	-----------------------------------------------------------------
	LCM  = LLF  + NSTEP   ! ADDED NODAL MASS SONGSAK AUG2007 SEE ALSO Subroutine NODALMS
	LCC  = LCM  + NEQ     ! ADDED NODAL DAMP SONGSAK AUG2007 SEE ALSO Subroutine NODALDP
	LRCT = LCC  + NEQ     ! FOR REACTION ADDED BY SONGSAK APR2006

	LCPZ  = LRCT + MFQ                      !ADD LCPZ FOR CABLE OPTIMIZATION SONGSAK NOV02006 
	LSTL  = LCPZ + 2*NEQ*NCABZ + NCABZ		!ADD LCPZ FOR CABLE OPTIMIZATION SONGSAK NOV02006 
	A(LCPZ:LSTL-1) = 0.0D0
C	-----------------------------------------------------------------
	LLV   = LSTL + NEQ*LCS                  !SETTLEMENT SONGSAK JUL2008 
	A(LSTL:LLV-1 ) = 0.0D0
C	-----------------------------------------------------------------

C     ---------------------------------
C	DISPLACEMENT AND LOAD VECTOR
C     ---------------------------------
	LLO   = LLV + NEQ*LCS                       !LLV = EXTERNAL LOAD VECTOR
      LLOF  = LLO + NEQ*LCS                       !LLOF = EXTERNAL LOAD VECTOR FOR OFFSHORE
      LRE   = LLOF + NEQ*LCS						!LLO = EXTERNAL LOAD VECTOR CONSTANT
      LDI   = LRE + NEQ
      LDL   = LDI + NEQ
      LDT   = LDL + NEQ
      LASTA = LDT + NEQ

C     ----------------------------------------------
C     ADDITIONAL STORAGE FOR LINEAR DYNAMIC ANALYSIS
C	----------------------------------------------
725	SELECT CASE (ISOLOP)

C	----------------------------
	CASE(5,6,7,8,9,10,11,15,16)
C	----------------------------
C	ADDITIONAL STORRAGE FOR DYNAMIC ANALYSIS

	LDT = LASTA
	LVL = LDT + NEQ
	LAL = LVL + NEQ
	LEF = LAL + NEQ
	LDU = LEF + NEQ
	LPR = LDU + NEQ

C	NEXT LINE ADDED BY GILSON - JULY2002
	LRV   = LPR   + 7*NSTEP
	LRT1  = LRV   + NEQ
	LRET  = LRT1  + NEQ
	LRET1 = LRET  + NEQ

	LDPT  = LRET1 + NEQ

	LVL1  = LDPT  + NEQ

	LASTA = LVL1  + NEQ

C	----------------------------
	CASE(1,2,3,4,12,13,14,18,19)  !	ADDED ISOLOP.EQ.18 CONSOLIDATION SACHARUCK DEC2006
C	----------------------------

750	LVL = LLV
	LAL = LLV
	LEF = LLV
	LDU = LLV
	LPR = LLV

	ENDSELECT

C     ----------------------------------------
C     STIFFNESS DIAGONAL (AFTER FACTERIZATION)
C     ----------------------------------------
775	LDK    = LASTA
	LASTA  = LDK + NEQ
C	-------------------------------------------------------

	RETURN
C     -------------------
C     EIGENVALUE ANALYSIS
C     -------------------
 800  MESSAGE = 'EIGENVALUE'
	LER = LASTA

      LTT = LER + NEQ*NC
      LWV = LTT + NEQ
      LAR = LWV + NEQ
      LBR = LAR + NNC
      LVE = LBR + NNC
      LEV = LVE + NC*NC
      LDD = LEV + NC
      LRT = LDD + NC
      LBU = LRT + NC
      LBC = LBU + NC
      MAXA= LBC + NC
      LASTA = MAX0(MAXA,LASTA)

      IF (MODEX.LE.1) RETURN

      RETURN

C     ----------------------------------
C     PRINT LOCATIONS OF FLEXIBLE ARRAYS
C     ----------------------------------
 900  IF (IFPR(4).EQ.0) RETURN
 910  IF (IFPR(4).EQ.0) GOTO 950
      NEG1   = NEG+1
      INDDIM = NBLOCK+NELTOT+1
      JXYDIM = MAX0(NSF,NSC)
      JCHDIM = MAX0(NSF,MAXCON)
      JBBDIM = ISTOR
      IF (NBLOCK.EQ.1 .AND.IEIG.EQ.0) JBBDIM = 0
      IF (NBLOCK.EQ.1 .AND.IEIG.GT.0) JBBDIM = NWM
      IF (MODEX.EQ.1)                 JBBDIM = 0
      NXY    = NCO*NNM
      MAXFL  = -1

      IF (LASTA.GT.MEMA) GOTO 920
      IF (IND.LE.3)  RETURN
 920  WRITE (ISO,3000)  LASTA,MEMA,MAXFL
      IF (LASTA.LE.MEMA)  RETURN
C     -----------------------------------
C     INSUFFICIENT STORAGE, PRINT MESSAGE
C     -----------------------------------
 950  MAXFL = -1


 3000 FORMAT (/48H LAST ADDRESS USED IN COMMON BLOCK A    LASTA = I6/   
     5 48H LAST ADMISSIBLE ADDRESS IN A           MEMA  = I6/           
     8 48H CENTRAL MEMORY ALLOCATION              MAXFL = I6)           
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE CORDIN (XY,LM,COORD,MSN,MEF,MAXN,MXY)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     SETS THE COORDINATES FOR ANY ELEMENT IN THIS GROUP
C	--------------------------------------------------
C     XY(NSN,NSC)         = COORDINATES FOR ENTIRE STRUCTURE
C     LM(NEF,NELE)        = ELEMENT COONECTIONS (NODAL INCIDENCES)
C     COORD(NCO*NNM,NELE) = ELEMENT GROUP COORDINATES
C     ------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C
      DIMENSION XY(MSN,1),LM(MEF,1),COORD(MXY,1)
C
      DO 200  IEL=1,NELE
      LOC = 0
      DO 200  INM=1,MAXN
      NODE = LM(INM,IEL)
      
      DO 100  ICO=1,NCO
      LOC = LOC+1
 100  IF (NODE.NE.0) COORD(LOC,IEL) = XY(NODE,ICO)

 200  CONTINUE


C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ELEOUT (LM,MTSET,IGSET,IAX,MEF,MAXN)
C	IAX - ADDED TO NEXT LINE BY GILSON - SEPT2002
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     WRITES ELEMENT CONNECTION TABLES TO OUTPUT
C	------------------------------------------
C     LM(NEF,NELE)  = ELEMENT CONNECTIVITY ARRAY
C     MTSET(NELE)   = MATERIAL SET NUMBERS
C     IGSET(NELE)   = GEOMETRIC SET NUMBERS
C	IAX(NELE)	  = MATERIAL AXIS FOR COMPOSIT
C     ------------------------------------------
      CHARACTER*4 ANODE
C
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	IAX - ADDED BY GILSON - SEPT2002
      DIMENSION LM(MEF,1),MTSET(1),IGSET(1),IAX(1)
C
      
      RETURN !SONGSAK TURNOFF PRINTING TO OUTPUT.OUT HERE...NO ONE WANT TO READ LARGE DATA OF ELEMENT CONNECTIVITY
      
      
      ANODE = 'NODE'
C	NEXT LINE ADDED BY GILSON - JUL2003 (MATERIAL AXIS)
C      WRITE (ISO,2000) (ANODE,I,I=1,10)


C	==================================================================
C	NEXT LINE ADDED BY SONGSAK--OCT2005 RC SHELL
	IF(ITYPE.EQ.9) THEN
	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6.AND.ISOLOP.EQ.2) GO TO 50
	ENDIF
	IF(ITYPE.EQ.10.AND.ISTYP.EQ.4) THEN
	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6) GO TO 50
	ENDIF
C	==================================================================

	IF (MTMOD.EQ.2) THEN
        WRITE (ISO,2000) (ANODE,I,I=1,10)
	ELSE
        WRITE (ISO,1000) (ANODE,I,I=1,10)
	ENDIF

	GO TO 90

C	==================================================================
C	NEXT LINE ADDED BY SONGSAK--OCT2005 RC SHELL
 50	WRITE (ISO,5000) (ANODE,I,I=1,10)
C	==================================================================

90	CONTINUE


C	CONDITION TO NEXT, ADDED BY GILSON SEPT2002
	IF (MTMOD.NE.2) THEN
        DO 100  IELE=1,NELE
 100    WRITE (ISO,1100) IELE,MTSET(IELE),IGSET(IELE),
     1                   (LM(INM,IELE),INM=1,NNM)
	ENDIF
C	NEXT IF BLOCK ADDED BY GILSON - SEPT2002
	IF (MTMOD.EQ.2) THEN
        DO 200  IELE=1,NELE
 200    WRITE (ISO,2200) IELE,MTSET(IELE),IAX(IELE),IGSET(IELE),
     1                   (LM(INM,IELE),INM=1,MAXN)
	ENDIF
C
 1000 FORMAT (//23X,35(1H*)/23X,1H*,33X,1H*/
     1        23X,35H*    NODAL ELEMENT CONNECTIONS    */
     2        23X,35H* MATERIAL, GEOMETRIC SET NUMBERS */
     3        23X,1H*,33X,1H*/23X,35(1H*)//
     4        20H ELEM M-SET G-SET   ,9(A4,I1,1X),A4,I2/1X,79(1H-)/)
 1100 FORMAT (3(I4,2X),10I6/(18X,10I6))
C
 2000 FORMAT (//23X,35(1H*)/23X,1H*,33X,1H*/
     1        23X,35H*    NODAL ELEMENT CONNECTIONS    */
     2        23X,35H* MATERIAL, GEOMETRIC SET NUMBERS */
     3        23X,1H*,33X,1H*/23X,35(1H*)//
     4        26H ELEM M-SET M-AXI G-SET   ,9(A4,I1,1X),
     5        A4,I2/1X,79(1H-)/)
 2200 FORMAT (4(I4,2X),10I6/(18X,10I6))

 5000 FORMAT (//23X,35(1H*)/23X,1H*,33X,1H*/
     1        23X,35H*    NODAL ELEMENT CONNECTIONS    */
     2        23X,35H* MATERIAL, GEOMETRIC SET NUMBERS */
     3        23X,1H*,33X,1H*/23X,35(1H*)//
     4        22H ELEM PATTERN G-SET   ,9(A4,I1,1X),A4,I2/1X,79(1H-)/)
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE CONECT (NCON,LM,MSN,MEF,MAXN,NELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     --------------------------------------------------------------
C     SETS UP NODAL CONNECTIVITY TABLE, WHICH IS LATER NEEDED FOR
C     BANDWIDTH OPTIMISATION. NOTE ZERO ENTRIES IN INCIDENCES ARRAY
C     (LM) INDICATE THAT THIS NODE WAS NOT USED TO DEFINE AN ELEMENT
C	--------------------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     NCON(NSN,50)  = NODAL CONNECTIVITY TABLE (MAX.CONNECTIONS = 50)
C                     NCON(I,J) = NODE NUMBER OF JTH CONNECTION TO
C                     NODE NUMBER I. A CONNECTION OF A NODE TO ITSELF
C                     IS NOT LISTED.
C     LM(NEF,NELE)  = ELEMENT CONNECTIVITY ARRAY
C     MSN = NSN     = TOTAL NUMBER OF NODES IN STRUCTURE
C     MEF = NEF     = NUMBER OF ELEMENT DEGREES OF FREEDOM
C     MAXN          = MAX. PERMITTED NUM OF NODES FOR ANY ONE ELEMENT
C     NELE          = NUMBER OF ELEMENTS IN THIS GROUP
C     MAXCON        = MAXIMUM PERMITTED NUMBER OF COUPLING NODES
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /BOPT/ IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C
      DIMENSION NCON(MSN,1),LM(MEF,1)
C      DIMENSION NCON(MSN,1),LM(MEF,NELE)
C
      IF (IREDU.GT.0) GOTO 200
C     -------------------------------------------
C     NO BANDWIDTH OPTIMISATION REQUESTED
C     FLAG ACTIVE NODES BY SETTING NCON(IACT) = 1
C     -------------------------------------------
      DO 100  IEL=1,NELE
      DO 100  INO=1,MAXN
      NODE = LM(INO,IEL)
      IF (NODE.EQ.0) GOTO 100
      NCON(NODE,1) = 1
 100  CONTINUE
      RETURN
C     ---------------------------------
C     LOOP OVER ALL ELEMENTS IN GROUP
C     LOOP OVER NUMBER OF ELEMENT NODES
C     ---------------------------------
 200  DO 900  IEL=1,NELE
      DO 900  INO=1,MAXN
      NOMAS = LM(INO,IEL)
      IF (NOMAS.EQ.0) GOTO 900
C     ----------------------------------------------------------
C     LOOP OVER NODES WHICH ARE COUPLED WITH MASTER NODE (NOMAS)
C     (CONNECTIONS OF A NODE TO ITSELF ARE NOT LISTED)
C     ----------------------------------------------------------
      DO 800  JNO=1,MAXN
      NOCOP = LM(JNO,IEL)
      IF (NOCOP.EQ.NOMAS) GOTO 800
      IF (NOCOP.EQ.0) GOTO 800
C     -----------------------------------------
C     FIND POSITION OF COUPLING NODE (NOCOP) IN
C     CONNECTIVITY TABLE (ROW NOMAS)
C     -----------------------------------------
      DO 500  JPOS=1,MAXCON
      IF (NCON(NOMAS,JPOS).EQ.0) NCON(NOMAS,JPOS) = NOCOP
      IF (NOCOP-NCON(NOMAS,JPOS)) 410,800,500
C     ------------------------------------------------------------
C     INSERT NEW COUPLING NODE AND SHIFT HIGHER COUPLING POSITIONS
C     ------------------------------------------------------------
 410  IF (JPOS.EQ.MAXCON) GOTO 510
      DO 490  KPOS=JPOS,MAXCON
      NOTEMP = NCON(NOMAS,KPOS)
      NCON(NOMAS,KPOS) = NOCOP
      NOCOP = NOTEMP
      IF (KPOS.GT.MCON) MCON = KPOS
 490  IF (NOTEMP.EQ.0) GOTO 800
C
 500  CONTINUE
 510  CALL ERRORS (20,IEL,MAXCON,'CONECT    ')
C     -------------------------------------------
C     END OF LOOP OVER COUPLING NODES
C     END OF LOOP OVER ELEMENT NODES AND ELEMENTS
C     -------------------------------------------
 800  CONTINUE
 900  CONTINUE
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE NEQUAT (ID,NCON,IOLD,MSF,MSN,LCON1,LCON2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     REMOVES INACTIVE NODES FROM SYSTEM OF SIMULTANEOUS EQUATIONS
C     COMPUTES NUMBER OF UNKNOWNS AND DEFINES IDENTIFICATION ARRAY
C	------------------------------------------------------------
C     ID(NSF,NSN)  = BOUNDARY CODES AT INPUT,RETURNS EQUATION NUMBERS
C     NCON(NSN,50) = NODAL CONNECTIVITY TABLE
C     IOLD(NSN)    = OLD NODE NUMBER FOR RENUMBERED NODE I
C     MSF = NSF    = NUMBER OF NODAL DEGREES OF FREEDOM FOR STRUCTURE
C     MSN = NSN    = NUMBER OF STRUCTURE NODES
C     NEQ          = NUMBER OF EQUATIONS (NEQ1 = NEQ+1)
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
C	==================================================================
C	JOINT GLOBAL CONSTRAINT SONGSAK JAN2007 
	COMMON /LC007/ LC107,LC207,LC307,LC407
C	==================================================================

C	CHANGED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
C	DIMENSION ID(MSF,1),NCON(MSN,1),IOLD(1),NUMEQ(7)
      DIMENSION ID(MSF,1),NCON(MSN,1),IOLD(1),NUMEQ(9)


	DIMENSION LCON1(1),LCON2(1)
C
C     -------------------------------------
C     SUPPRESS EQUATIONS FOR INACTIVE NODES
C     -------------------------------------
      DO 190  ISN=1,NSN
      IF (NCON(ISN,1).GT.0) GOTO 190
      DO 150  ISF=1,NSF
 150  ID(ISF,ISN) = 1
 190  CONTINUE
C     ----------------------------------------------------------
C     ASSIGN EQUATION NUMBERS
C     THE SEQUENCE IS GIVEN BY THE ORDER OF THE RENUMBERED NODES
C     ----------------------------------------------------------
      NEQ = 0
      DO 290  ISN=1,NSN
      NODE = IOLD(ISN)
      DO 290  ISF=1,NSF
      IF (ID(ISF,NODE)) 250,220,250
 220  NEQ = NEQ+1
      ID(ISF,NODE) = NEQ
      GOTO 290
 250  ID(ISF,NODE) = 0
 290  CONTINUE

C	JOINT GLOBAL CONSTRAINT SONGSAK JAN2007
	IF(LC107.NE.0) THEN
	CALL JONSTRN(ID,IDOF,MSF,MSN,NEQ,NCEQ,LCON1,LCON2)
	NEQ = NCEQ
	ENDIF

      NEQ1 = NEQ+1
C     --------------------------------------
C     PRINT EQUATION NUMBERS (IFPR(2)=10/11)
C     --------------------------------------
      IFPR1 = IFPR(2)/10
      
      RETURN  !SONGSAK SKIP PRINTING IN OUTPUT.OUT OCT2019... NO ONE WANT TO SEE LARGE OUTPUT ON THAT FILE
      
C	IF (IFPR1.EQ.0) RETURN    !SUPPRESSED BY SACHARUCK MAR2007  (SONGSAK IMPLEMENTOR)
      WRITE (ISO,1000)

C	CHANGED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
C	CALL CLEARI (NUMEQ,7)
	CALL CLEARI (NUMEQ,9)

      DO 390  ISN=1,NSN
      DO 350  ISF=1,NSF
      IPO = IDOF(ISF)
 350  NUMEQ(IPO) = ID(ISF,ISN)
 390  WRITE (ISO,1100) ISN,NUMEQ
C

C	CHANGED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
C1000 FORMAT (1H1///,30X,20(1H*)/30X,1H*,18X,1H*/
C	1        30X,20H* EQUATION NUMBERS */
C	2        30X,1H*,18X,1H*/30X,20(1H*)////
C	3        22X,9HJOINT NR.,10X,17HDEGREE OF FREEDOM//
C	4        22X,2HIJ,9X,2HDX,3X,2HDY,3X,2HDZ,3X,2HRX,3X,2HRY,3X,2HRZ/
C	5        22X,38(1H-)//)
C1100 FORMAT (19X,I5,6X,7I5)
 1000 FORMAT (1H1///,30X,20(1H*)/30X,1H*,18X,1H*/
	1        30X,20H* EQUATION NUMBERS */
     2        30X,1H*,18X,1H*/30X,20(1H*)////
     3        22X,9HJOINT NR.,10X,17HDEGREE OF FREEDOM//
     4        22X,2HIJ,9X,2HDX,3X,2HDY,3X,2HDZ,3X,2HRX,3X,2HRY,3X,2HRZ,
     5		3X,2HWp,3X,2HT.,3X,2HP./22X,52(1H-)//)
 1100 FORMAT (19X,I5,6X,9I9)


C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE READEL (LEST)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     LOOPS OVER NUMBER OF ELEMENT GROUPS AND READS ELEMENT INFORMATION
C     FROM RANDOM ACCESS TAPE NELEM.
C     CALLS COLHEI TO ESTABLISH GLOBAL EQUATION NUMBERS FOR ELEMENT
C     DEGREES OF FREEDOM AND TO DETERMINE COLUMN HEIGHTS OR BANDWIDTH
C     OF GLOBAL COMPACTED STIFFNESS MATRIX
C	------------------------------------
C     LEST(2*NEG) = LENGTHS OF ELEMENT BLOCK INFORMATION (ARRAY A & IA)
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C	==================================================================
C	ADDED BY SONGSAK APR2006 REACTION	
	COMMON /REACT/ LRC,LRCT,MFQ,LRID
C	==================================================================

	COMMON A(9000000),IA(9000000)
C
      DIMENSION LEST(1)
C
      DO 900  IEG=1,NEG
      NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     

      CALL MOVLEV (2)
	
C
      CALL COLHEI (IA(LID),IA(LCH),IA(LLM),NSF,NEF,IA(LRC),IA(LRID))    ! ADDED BY SONGSAK APR2006 FOR REACTION
      REWIND NELEMI
      REWIND NELEMA
C      WRITE (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      WRITE (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      WRITE (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      WRITE (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     
 900  CONTINUE
C
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE COLHEI (ID,MHT,LM,MSF,MEF,LMRCT,IDRCT) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     DETERMINES GLOBAL EQUATION NUMBERS FOR ELEMENT DEGREES OF FREEDOM
C     UPDATES COLUMN HEIGHTS AND BANDWIDTH
C	------------------------------------
C     ID(NSF,NSN)   = EQUATION NUMBERS
C     MHT(NEQ)      = COLUMN HEIGHTS AND BANDWIDTH OF GLOBAL COMP.[K]
C     LM(NEF,NELE)  = NODAL ELEMENT CONNECTIONS AT ENTRY, RETURNS
C     MSF = NSF     = NUMBER OF NODAL DEGREES OF FREEDOM FOR STRUCTURE
C     MEF = NEF     = NUMBER OF ELEMENT DEGREES OF FREEDOM
C     NNF           = NUMBER OF NODAL DEGREES OF FREEDOM FOR THIS GROUP
C     IDOF(NSF)     = GLOBAL POSITIONS OF NODAL STRUCTURE D.O.F.
C     NELE          = NUMBER OF ELEMENTS IN THIS GROUP
C                     GLOBAL EQUATION NUMBER FOR ANY ELEMENT D.O.F.
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON A(9000000),IA(9000000)

C	ADDED BY SONGSAK APR2006 REACTION	
	COMMON /REACT/ LRC,LRCT,MFQ,LRID

C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)

C	ADDED BY SONGSAK APR2006 REACTION
	DIMENSION IDRCT(MSF,NSN),LMRCT(NEF,NELE)

      DIMENSION ID(MSF,1),MHT(1),LM(MEF,1),IGPOS(9),INCI(200)
	DIMENSION IEPOS(9),IMPOS(9)

C     ---------------------------------------------------------------------
C	GLOBAL BOUNDARY EQUATION NUMBER FOR REACTION ADDED BY SONGSAK APR2006
C     ---------------------------------------------------------------------
	MRC = 0
	DO ISN = 1,NSN
	DO ISF = 1,NSF
	IDRCT(ISF,ISN) = 0
	IF(ID(ISF,ISN).EQ.0) THEN
	MRC = MRC + 1
	IDRCT(ISF,ISN) = MRC
	ENDIF
	ENDDO
	ENDDO
	MFQ = MRC


C     ---------------------------------------------------
C     UNPACK VARIABLE LINKF
C     SET GLOBAL POSITIONS IGPOS FOR NODAL ELEMENT D.O.F.
C     ---------------------------------------------------
      LINK = NFLINK(ITYPE,ISTYP)

	IEPOS(1:NNF) = 0
      DO INF=1,NNF
      IEPOS(INF) = LINK/10**(NNF-INF)
	LINK = LINK - 10**(NNF-INF)*IEPOS(INF)
	ENDDO

	IMPOS(1:9) = 0
      DO INF=1,NNF
	II = 0
	DO I = 1,9
	IF(IEPOS(INF).EQ.IDOF(I)) II =  I
	ENDDO
	IF(II.NE.0) IMPOS(II) = INF
	ENDDO

      LINK = NFLINK(ITYPE,ISTYP)
	IGPOS(1:NNF) = 0
      IGF  = 0
      DO 190  INF=1,NNF
      IGPOS(INF) = LINK/10**(NNF-INF)
 150  IGF = IGF+1
	IF (IGF.GT.9) THEN
	IGPOS(INF) = 0
	GOTO 195
	ENDIF
      IF (IDOF(IGF)-IGPOS(INF)) 150,170,170
 170  LINK = LINK - 10**(NNF-INF)*IGPOS(INF)
 190  IGPOS(INF) = IGF
 195	CONTINUE


C     --------------------------------------------
C     RELATE GLOBAL EQUATION NUM TO ELEMENT D.O.F.
C     --------------------------------------------
      DO 500  IEL=1,NELE

      DO 210  IDF=1,NEF
      INCI(IDF)      = LM(IDF,IEL)
	LMRCT(IDF,IEL) = 0							   !REACTION SONGSAK APR2006
 210  LM(IDF,IEL)    = 0

      KF = -NNF
      DO 290  IDF=1,NEF
      NODE = INCI(IDF)
      KF = KF+NNF
      IF (NODE.EQ.0)  GOTO 290
      
      DO 250  INF=1,9
	IMF = IMPOS(INF)
	IF (IMF.LE.0)  GOTO 250
	IGF = IGPOS(IMF)
	LMRCT(KF+IMF,IEL) = IDRCT(IGF,NODE)            !REACTION SONGSAK APR2006
	   LM(KF+IMF,IEL) =    ID(IGF,NODE)
 250	CONTINUE
 
 290  CONTINUE

C     ---------------------------------------------
C     UPDATE COLUMN HEIGHTS OF GLOBAL COMPACTED [K]
C     ---------------------------------------------
      MEQ = 0
      DO 390  IDF=1,NEF
      IEQ = LM(IDF,IEL)
      IF (IEQ)     310,390,310
 310  IF (MEQ.EQ.0) MEQ = IEQ
      IF (MEQ-IEQ) 320,390,390
 320  MEQ = IEQ
 390  CONTINUE
C
      DO 400  IDF=1,NEF
      IEQ = LM(IDF,IEL)
      IF (IEQ.EQ.0) GOTO 400
      KHT = MEQ-IEQ
      IF (KHT.GT.MHT(IEQ)) MHT(IEQ) = KHT
 400	CONTINUE

 500  CONTINUE
C
      IFPR1 = IFPR(2)/10
      IFPRI = IFPR(2) - 10*IFPR1
      IF (IFPRI.EQ.0) RETURN
      WRITE (ISO,1200)
      CALL MATOUT (iA(LLM),NEF,NELE,10,5,'I',0,12,1,'COLUMN LM ')
C
 1200 FORMAT (////,23X,31(1H*)/23X,1H*,29X,1H*/
     1        23X,31H* COLUMN AND ROW NUMBERS (LM) */
     2        23X,1H*,29X,1H*/23X,31(1H*))

C
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ADRESS (MAXA,MHT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SUNIL 01/05/11 MODIFIED
C     ----------------------------------------------------------------
C     CALCULATES ADDRESSES OF DIAGONAL ELEMENTS IN GLOBAL COMPACTED
C     STIFFNESS MATRIX
C	----------------
C     MAXA(NEQ1)  = ADDRESSES OF DIAGONAL ELEMENTS IN STIFFNESS
C     MHT(NEQ)    = COLUMN HEIGHTS
C     NEQ         = NUMBER OF EQUATIONS (NEQ1=NEQ+1)
C     MK          = MAXIMUM COLUMN HEIGHT
C     BM          = MEAN BANDWIDTH
C     NWK         = NUMBER OF COEFF. BELOW SKYLINE OF STIFFNESS MATRIX
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /DYNA/ CDEN,IMASS
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C	NEXT COMMON ADDED BY GILSON - MAY2003 (ARC-LENGTH)
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
	3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
C
C	-------------------------------------------------------------------
      DIMENSION MAXA(1),MHT(1)
C
	NWM = 0
      MAXA(1) = 1
      MK      = 0
      DO 250  IEQ=1,NEQ
      IF (MHT(IEQ).GT.MK) MK = MHT(IEQ)
 250  MAXA(IEQ+1) = MAXA(IEQ) + MHT(IEQ) + 1
      MK  = MK+1
      NWK = MAXA(NEQ1) - MAXA(1)
      BM  = FLOAT(NWK)/FLOAT(NEQ)
C
	NWM = NWK
C
	IFPR2 = IFPR(3)/10
      IF (IFPR2.EQ.0) RETURN
      WRITE (ISO,1000)
      WRITE (ISO,1100) (IEQ,MHT(IEQ),MAXA(IEQ), IEQ=1,NEQ)
      WRITE (ISO,1200) MAXA(NEQ+1)
C
 1000 FORMAT (////18X,41(1H*)/18X,1H*,39X,1H*/
     1       18X,41H* COLUMN  HEIGHTS  AND  BANDWIDTH (MHT) */
     2       18X,41H* ADDRESSES OF DIAGONAL ELEMENTS (MAXA) */
     3       18X,1H*,39X,1H*/18X,41(1H*)////
     4       22X,11HCOLUMN NR.I,5X,6HMHT(I),5X,7HMAXA(I)/22X,34(1H-)//)
 1100 FORMAT ((25X,3(I5,7X)))
 1200 FORMAT (49X,I5)
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE REDUCE(NDSTK, IOLD, RENUM, NDEG, LVL, LVLS1,
     * LVLS2, CCSTOR, NR)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C  SUBROUTINE REDUCE DETERMINES A ROW AND COLUMN PERMUTATION WHICH,
C  WHEN APPLIED TO A GIVEN SPARSE MATRIX, PRODUCES A PERMUTED
C  MATRIX WITH A SMALLER BANDWIDTH AND PROFILE.
C  THE INPUT ARRAY IS A CONNECTION TABLE WHICH REPRESENTS THE
C  INDICES OF THE NONZERO ELEMENTS OF THE MATRIX, A.  THE ALGO-
C  RITHM IS DESCRIBED IN TERMS OF THE ADJACENCY GRAPH WHICH
C  HAS THE CHARACTERISTIC THAT THERE IS AN EDGE (CONNECTION)
C  BETWEEN NODES I AND J IF A(I,J) .NE. 0 AND I .NE. J.
C  DIMENSIONING INFORMATION--THE FOLLOWING INTEGER ARRAYS MUST BE
C  DIMENSIONED IN THE CALLING ROUTINE.
C    NDSTK(NR,D1)        D1 IS .GE. MAXIMUM DEGREE OF ALL NODES.
C    IOLD(D2)            D2 AND NR ARE .GE. THE TOTAL NUMBER OF
C    RENUM(D2+1)         NODES IN THE GRAPH.
C    NDEG(D2)            STORAGE REQUIREMENTS CAN BE SIGNIFICANTLY
C    LVL(D2)             DECREASED FOR IBM 360 AND 370 COMPUTERS
C    LVLS1(D2)           BY REPLACING INTEGER NDSTK BY
C    LVLS2(D2)           INTEGER*2 NDSTK IN SUBROUTINES REDUCE,
C    CCSTOR(D2)          DGREE, FNDIAM, TREE AND NUMBER.
C  COMMON INFORMATION--THE FOLLOWING COMMON BLOCK MUST BE IN THE
C  CALLING ROUTINE.
C     COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,IREDU
C  EXPLANATION OF INPUT VARIABLES--
C    NDSTK-     CONNECTION TABLE REPRESENTING GRAPH.
C               NDSTK(I,J)=NODE NUMBER OF JTH CONNECTION TO NODE
C               NUMBER I.  A CONNECTION OF A NODE TO ITSELF IS NOT
C               LISTED.  EXTRA POSITIONS MUST HAVE ZERO FILL.
C    NR-        ROW DIMENSION ASSIGNED NDSTK IN CALLING PROGRAM.
C    IOLD(I)-   NUMBERING OF ITH NODE UPON INPUT.
C               IF NO NUMBERING EXISTS THEN IOLD(I)=I.
C    LCON-      NUMBER OF NODES IN GRAPH (EQUAL TO ORDER OF MATRIX).
C    MCON-      MAXIMUM DEGREE OF ANY NODE IN THE GRAPH.
C  EXPLANATION OF OUTPUT VARIABLES--
C    RENUM(I)-  THE NEW NUMBER FOR THE ITH NODE.
C    NDEG(I)-   THE DEGREE OF THE ITH NODE.
C    IBW2-      THE BANDWIDTH AFTER RENUMBERING.
C    IPF2-      THE PROFILE AFTER RENUMBERING.
C    IDPTH-     NUMBER OF LEVELS IN REDUCE LEVEL STRUCTURE.
C  THE FOLLOWING ONLY HAVE MEANING IF THE GRAPH WAS CONNECTED--
C    LVL(I)-    INDEX INTO LVLS1 TO THE FIRST NODE IN LEVEL I.
C               LVL(I+1)-LVL(I)= NUMBER OF NODES IN ITH LEVEL
C    LVLS1-     NODE NUMBERS LISTED BY LEVEL.
C    LVLS2(I)-  THE LEVEL ASSIGNED TO NODE I BY REDUCE.
C  WORKING STORAGE VARIABLE--
C    CCSTOR
C  LOCAL STORAGE--
C    COMMON/CC/-SUBROUTINES REDUCE, SORT2 AND PIKLVL ASSUME THAT
C               THE GRAPH HAS AT MOST 50 CONNECTED COMPONENTS.
C               SUBROUTINE FNDIAM ASSUMES THAT THERE ARE AT MOST
C               100 NODES IN THE LAST LEVEL.
C    COMMON/LVLW/-SUBROUTINES SETUP AND PIKLVL ASSUME THAT THERE
C               ARE AT MOST 100 LEVELS.
C USE INTEGER*2 NDSTK  WITH AN IBM 360 OR 370.
      INTEGER*4 NDSTK
      INTEGER*4 STNODE, RVNODE, RENUM, XC, SORT2, STNUM, CCSTOR,
     * SIZE1, STPT, SBNUM
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /BOPT/ IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C IT IS ASSUMED THAT THE GRAPH HAS AT MOST 50 CONNECTED COMPONENTS.
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level)
      COMMON /CC/ XC, SIZE1(100), STPT(100)
      COMMON /LVLW/ NHIGH(200), NLOW(200), NACUM(200)
      DIMENSION CCSTOR(1), IOLD(1)
      DIMENSION NDSTK(NR,1), LVL(1), LVLS1(1), LVLS2(1), RENUM(1),
     * NDEG(1)
C
C     -----------------------------------------
C     INITIALISATION, PRINT CONNECTIVITY MATRIX (IFPR(7)>0)
C     RENUM(I)=0 INDICATES NODE I IS UNNUMBERED
C     SBNUM = LOW END OF AVAILABLE NUMBERS FOR RENUMBERING
C     STNUM = HIGH END OF AVAILABLE NUMBERS FOR RENUMBERING
C     -----------------------------------------------------
C	SUNIL 20/01/01 NEXT LINE REMOVED (NEVER USED)
C      NAN   = NSN
      IBW2  = 0
      IPF2  = 0
      SBNUM = 1
      STNUM = LCON
      DO 10  I=1,LCON
      IOLD(I) = I
 10   RENUM(I) = 0
      IF (IREDU.EQ.0)  RETURN
      IF (IFPR(1).LE.0)  GOTO 110
      WRITE (ISO,1000)
      CALL MATOUT (NDSTK,LCON,MCON,20,6,'I',0,5,1,'CONNECTION')
C
C     -----------------------------------------------------------
C     COMPUTE DEGREE OF EACH NODE, ORIGINAL BANDWIDTH AND PROFILE
C     NUMBER THE NODES OF DEGREE ZERO
C     -----------------------------------------------------------
 110  CALL DEGREE (NDSTK,IOLD,NDEG,LCON)
      DO 200  ICON=1,LCON
      IF (NDEG(ICON).GT.0)  GOTO 200
      RENUM(ICON) = STNUM
      STNUM = STNUM-1
 200  CONTINUE
C	SUNIL 20/01/01 NEXT LINE REMOVED (NEVER USED)
C      NAN = STNUM
C FIND AN UNNUMBERED NODE OF MIN DEGREE TO START ON
   30 LOWDG = MCON + 1
      NFLG = 1
      ISDIR = 1
      DO 40 I=1,LCON
        IF (NDEG(I).GE.LOWDG) GOTO 40
        IF (RENUM(I).GT.0) GOTO 40
        LOWDG = NDEG(I)
        STNODE = I
   40 CONTINUE
C FIND PSEUDO-DIAMETER AND ASSOCIATED LEVEL STRUCTURES.
C STNODE AND RVNODE ARE THE ENDS OF THE DIAM AND LVLS1 AND LVLS2
C ARE THE RESPECTIVE LEVEL STRUCTURES.
      CALL FNDIAM(STNODE, RVNODE, NDSTK, NR, NDEG, LVL, LVLS1,
     * LVLS2, CCSTOR, IDFLT)
      IF (IREDU.EQ.2)  GOTO 50
      IF (NDEG(STNODE).LE.NDEG(RVNODE)) GOTO 50
C NFLG INDICATES THE END TO BEGIN NUMBERING ON
      NFLG = -1
      STNODE = RVNODE
   50 CALL SETUP(LVL, LVLS1, LVLS2)
C FIND ALL THE CONNECTED COMPONENTS  (XC COUNTS THEM)
      XC = 0
      LROOT = 1
      LVLN = 1
      DO 60 I=1,LCON
        IF (LVL(I).NE.0) GOTO 60
        XC = XC + 1
        STPT(XC) = LROOT
        CALL TREE(I, NDSTK, NR, LVL, CCSTOR, NDEG, LVLWTH, LVLBOT,
     *   LVLN, MAXLW, LCON)
        SIZE1(XC) = LVLBOT + LVLWTH - LROOT
        LROOT = LVLBOT + LVLWTH
        LVLN = LROOT
   60 CONTINUE
      IF (SORT2(DMY).EQ.0) GOTO 70
      CALL PIKLVL(LVLS1, LVLS2, CCSTOR, IDFLT, ISDIR)
C ON RETURN FROM PIKLVL, ISDIR INDICATES THE DIRECTION THE LARGEST
C COMPONENT FELL.  ISDIR IS MODIFIED NOW TO INDICATE THE NUMBERING
C DIRECTION.  NUM IS SET TO THE PROPER VALUE FOR THIS DIRECTION.
   70 ISDIR = ISDIR*NFLG
      NUM = SBNUM
C
C     -------------------------------------------------------
C     BANDWIDTH OPTIMISATION MORE IMPORTANT THAN PROFILE RED.
C     -------------------------------------------------------
      IF (IREDU.EQ.2)  GOTO 90
      IF (ISDIR.LT.0) NUM = STNUM
      CALL NUMBER(STNODE, NUM, NDSTK, LVLS2, NDEG, RENUM, LVLS1,
     * LVL, NR, NFLG, CCSTOR, ISDIR)
C UPDATE STNUM OR SBNUM AFTER NUMBERING
      IF (ISDIR.LT.0) STNUM = NUM
      IF (ISDIR.GT.0) SBNUM = NUM
      IF (SBNUM.LE.STNUM) GOTO 30
      IF (IBW2.LE.IBW1) GOTO 900
C     -----------------------------------------
C     ORIGINAL NUMBERING IS BETTER THAN NEW ONE
C     -----------------------------------------
 750  DO 800  ICON=1,LCON
 800  RENUM(ICON) = IOLD(ICON)
C	SUNIL 20/01/01 NEXT LINE REMOVED (NEVER USED)
C      NAN  = NSN
      IBW2 = IBW1
      IPF2 = IPF1
      GOTO 900
C
C     ----------------------------------------------------
C     PROFILE REDUCTION MORE IMPORTANT THAN BANDWIDTH OPT.
C     ----------------------------------------------------
 90   CALL PROFIT (NR,NDSTK,RENUM,NDEG,LVLS2,LVLS1,LVL,NUM)
      SBNUM = NUM
      IF (SBNUM.LE.STNUM)  GOTO 30
C
C     ----------------------------------------------------
C     CHECK WHETHER PROFILE IS IMPROVED BY REVERSING NUMBERING
C     CHECK WHETHER ORIGINAL NUMBERING IS BETTER THAN NEW ONE
C     --------------------------------------------------------
      CALL CHECK (IBW2,IPF2,RENUM,NDSTK,NR,NDEG,LVL)
      IF (IPF2.LE.IPF1)  GOTO 900
      GOTO 750
C     -----------------------------------------------
C     FIND INVERSE RELATIONSHIP NEW NODE TO OLD NODE
C     PRINT NEW NODES AFTER RENUMBERING (IFPR(7)>0)
C     ----------------------------------------------
 900  DO 950  ICON=1,LCON
      INODE = RENUM(ICON)
 950  IOLD(INODE) = ICON
      IF (IFPR(1).LE.0)  GOTO 960
      WRITE (ISO,1100)  (I,IOLD(I),RENUM(I),NDEG(I),I=1,LCON)
      WRITE (ISO,1200)  IBW1,IPF1,IBW2,IPF2,IDPTH
 960  IBW1 = IBW1*NSF
      IPF1 = IPF1*NSF*NSF - LCON*(NSF*NSF-NSF)/2
C
 1000 FORMAT (1H1//20X,40(1H*)/20X,1H*,38X,1H*/
     1        20X,40H* BANDWIDTH OPTIMISATION (GIBBS ET AL) */
     2        20X,40H* CONNECTIVITY MATRIX,RENUMBERED NODES */
     3        20X,1H*,38X,1H*/20X,40(1H*)//)
 1100 FORMAT (///35H NODAL MAPPING MATRICES IOLD,RENUM /
     1        29H DEGREES FOR EACH NODE  NDEG //
     2        3X,1HI,3X,4HIOLD,2X,5HRENUM,2X,4HNDEG/1X,23(1H-)//
     3        (4I6))
 1200 FORMAT (//15X,45HMAX.SEMI-BANDWIDTH BEFORE RENUMBERING IBW1 = ,I5/
     1          15X,45HCOEF.BELOW SKYLINE BEFORE RENUMBERING IPF1 = ,I5/
     2           15X,45HMAX.SEMI-BANDWIDTH AFTER  RENUMBERING IBW2 = ,I5
     3/           15X,45HCOEF.BELOW SKYLINE AFTER  RENUMBERING IPF2 = ,I
     45/           15X,45HNO. OF LEVELS IN REDUCE LEVEL STRUCT.IDPTH = ,
     4I5)
C
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE DEGREE (NCON,IOLD,NDEG,MR)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C     --------------------------------------------------------------
C     COMPUTES THE DEGREE OF EACH NODE AND BANDWIDTH AND PROFILE OF
C     ORIGINAL OR INPUT RENUMBERING
C
C     NCON(LCON,MCON)  = NODAL CONNECTIVITY TABLE
C     IOLD(LCON)       = ORIGINAL NUMBERING
C     NDEG(LCON)       = DEGREE OF EACH NODE
C     IPF1             = ORIGINAL PROFILE
C     IBW1             = ORIGINAL BANDWIDTH
C     --------------------------------------------------------------
C
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
      DIMENSION NCON(MR,1),IOLD(1),NDEG(1)
C
      IBW1 = 0
      IPF1 = 0
      DO 300  ICON=1,LCON
      NDEG(ICON) = 0
      IRW = 0
      DO 200  JCON=1,MCON
      NODE = NCON(ICON,JCON)
      IF (NODE)  210,210,110
 110  NDEG(ICON) = NDEG(ICON) + 1
      IDIF = IOLD(ICON) - IOLD(NODE)
 200  IF (IRW.LT.IDIF)  IRW = IDIF
C
 210  IPF1 = IPF1+IRW
 300  IF (IRW.GT.IBW1)  IBW1 = IRW
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE FNDIAM(SND1, SND2, NDSTK, NR, NDEG, LVL, LVLS1,
     * LVLS2, IWK, IDFLT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C  FNDIAM IS THE CONTROL PROCEDURE FOR FINDING THE PSEUDO-DIAMETER OF
C  NDSTK AS WELL AS THE LEVEL STRUCTURE FROM EACH END
C  SND1-        ON INPUT THIS IS THE NODE NUMBER OF THE FIRST
C               ATTEMPT AT FINDING A DIAMETER.  ON OUTPUT IT
C               CONTAINS THE ACTUAL NUMBER USED.
C  SND2-        ON OUTPUT CONTAINS OTHER END OF DIAMETER
C  LVLS1-       ARRAY CONTAINING LEVEL STRUCTURE WITH SND1 AS ROOT
C  LVLS2-       ARRAY CONTAINING LEVEL STRUCTURE WITH SND2 AS ROOT
C  IDFLT-       FLAG USED IN PICKING FINAL LEVEL STRUCTURE, SET
C               =1 IF WIDTH OF LVLS1 .LE. WIDTH OF LVLS2, OTHERWISE =2
C  LVL,IWK-     WORKING STORAGE
C USE INTEGER*2 NDSTK  WITH AN IBM 360 OR 370.
      INTEGER*4 NDSTK
      INTEGER*4 FLAG, SND, SND1, SND2
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C IT IS ASSUMED THAT THE LAST LEVEL HAS AT MOST 100 NODES.
C	SUNIL 18/10/10 CHANGE COMMON BLOCK CC TO CC1 SINCE CC HAVING
C	A DIFFERENT NUMBER OF ARGUMENTS
C      COMMON /CC/ NDLST(100)
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level)
      COMMON /CC1/ NDLST(200)
      DIMENSION NDSTK(NR,1), NDEG(1), LVL(1), LVLS1(1), LVLS2(1),
     * IWK(1)
      FLAG = 0
      MTW2 = LCON
      SND = SND1
C ZERO LVL TO INDICATE ALL NODES ARE AVAILABLE TO TREE
   10 DO 20 I=1,LCON
        LVL(I) = 0
   20 CONTINUE
      LVLN = 1
C DROP A TREE FROM SND
      CALL TREE(SND, NDSTK, NR, LVL, IWK, NDEG, LVLWTH, LVLBOT,
     * LVLN, MAXLW, MTW2)
      IF (FLAG.GE.1) GOTO 50
      FLAG = 1
   30 IDPTH = LVLN - 1
      MTW1 = MAXLW
C COPY LEVEL STRUCTURE INTO LVLS1
      DO 40 I=1,LCON
        LVLS1(I) = LVL(I)
   40 CONTINUE
      NDXN = 1
      NDXL = 0
      MTW2 = LCON
C SORT LAST LEVEL BY DEGREE  AND STORE IN NDLST
      CALL SORTDG(NDLST, IWK(LVLBOT), NDXL, LVLWTH, NDEG)
      SND = NDLST(1)
      GOTO 10
   50 IF (IDPTH.GE.LVLN-1) GOTO 60
C START AGAIN WITH NEW STARTING NODE
      SND1 = SND
      GOTO 30
   60 IF (MAXLW.GE.MTW2) GOTO 80
      MTW2 = MAXLW
      SND2 = SND
C STORE NARROWEST REVERSE LEVEL STRUCTURE IN LVLS2
      DO 70 I=1,LCON
        LVLS2(I) = LVL(I)
   70 CONTINUE
   80 IF (NDXN.EQ.NDXL) GOTO 90
C TRY NEXT NODE IN NDLST
      NDXN = NDXN + 1
      SND = NDLST(NDXN)
      GOTO 10
   90 IDFLT = 1
      IF (MTW2.LE.MTW1) IDFLT = 2
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE TREE(IROOT, NDSTK, NR, LVL, IWK, NDEG, LVLWTH,
     * LVLBOT, LVLN, MAXLW, IBORT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C  TREE DROPS A TREE IN NDSTK FROM IROOT
C  LVL-         ARRAY INDICATING AVAILABLE NODES IN NDSTK WITH ZERO
C               ENTRIES. TREE ENTERS LEVEL NUMBERS ASSIGNED
C               DURING EXECUTION OF THIS PROCEDURE
C  IWK-         ON OUTPUT CONTAINS NODE NUMBERS USED IN TREE
C               ARRANGED BY LEVELS (IWK(LVLN) CONTAINS IROOT
C               AND IWK(LVLBOT+LVLWTH-1) CONTAINS LAST NODE ENTERED)
C  LVLWTH-      ON OUTPUT CONTAINS WIDTH OF LAST LEVEL
C  LVLBOT-      ON OUTPUT CONTAINS INDEX INTO IWK OF FIRST
C               NODE IN LAST LEVEL
C  MAXLW-       ON OUTPUT CONTAINS THE MAXIMUM LEVEL WIDTH
C  LVLN-        ON INPUT THE FIRST AVAILABLE LOCATION IN IWK
C               USUALLY ONE BUT IF IWK IS USED TO STORE PREVIOUS
C               CONNECTED COMPONENTS, LVLN IS NEXT AVAILABLE LOCATION.
C               ON OUTPUT THE TOTAL NUMBER OF LEVELS + 1
C  IBORT-       INPUT PARAM WHICH TRIGGERS EARLY RETURN IF
C               MAXLW BECOMES .GE. IBORT
C USE INTEGER*2 NDSTK  WITH AN IBM 360 OR 370.
      INTEGER*4 NDSTK
      DIMENSION NDSTK(NR,1), LVL(*), IWK(*), NDEG(*)
      MAXLW = 0
      ITOP = LVLN
      INOW = LVLN
      LVLBOT = LVLN
      LVLTOP = LVLN + 1
      LVLN = 1
      LVL(IROOT) = 1
      IWK(ITOP) = IROOT
   10 LVLN = LVLN + 1
   20 IWKNOW = IWK(INOW)
      NDROW = NDEG(IWKNOW)
      DO 30 J=1,NDROW
        ITEST = NDSTK(IWKNOW,J)
        IF (LVL(ITEST).NE.0) GOTO 30
        LVL(ITEST) = LVLN
        ITOP = ITOP + 1
        IWK(ITOP) = ITEST
   30 CONTINUE
      INOW = INOW + 1
      IF (INOW.LT.LVLTOP) GOTO 20
      LVLWTH = LVLTOP - LVLBOT
      IF (MAXLW.LT.LVLWTH) MAXLW = LVLWTH
      IF (MAXLW.GE.IBORT) RETURN
      IF (ITOP.LT.LVLTOP) RETURN
      LVLBOT = INOW
      LVLTOP = ITOP + 1
      GOTO 10
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SORTDG(STK1, STK2, X1, X2, NDEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C SORTDG SORTS STK2 BY DEGREE OF THE NODE AND ADDS IT TO THE END
C OF STK1 IN ORDER OF LOWEST TO HIGHEST DEGREE.  X1 AND X2 ARE THE
C NUMBER OF NODES IN STK1 AND STK2 RESPECTIVELY.
      INTEGER*4 X1, X2, STK1, STK2, TEMP
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
      DIMENSION NDEG(1), STK1(1), STK2(1)
      IND = X2
   10 ITEST = 0
      IND = IND - 1
      IF (IND.LT.1) GOTO 30
      DO 20 I=1,IND
        J = I + 1
        ISTK2 = STK2(I)
        JSTK2 = STK2(J)
        IF (NDEG(ISTK2).LE.NDEG(JSTK2)) GOTO 20
        ITEST = 1
        TEMP = STK2(I)
        STK2(I) = STK2(J)
        STK2(J) = TEMP
   20 CONTINUE
      IF (ITEST.EQ.1) GOTO 10
   30 DO 40 I=1,X2
        X1 = X1 + 1
        STK1(X1) = STK2(I)
   40 CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SETUP(LVL, LVLS1, LVLS2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C SETUP COMPUTES THE REVERSE LEVELING INFO FROM LVLS2 AND STORES
C IT INTO LVLS2.  NACUM(I) IS INITIALIZED TO NODES/ITH LEVEL FOR NODES
C ON THE PSEUDO-DIAMETER OF THE GRAPH.  LVL IS INITIALIZED TO NON-
C ZERO FOR NODES ON THE PSEUDO-DIAM AND NODES IN A DIFFERENT
C COMPONENT OF THE GRAPH.
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C IT IS ASSUMED THAT THERE ARE AT MOST 100 LEVELS.
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level)
      COMMON /LVLW/ NHIGH(200), NLOW(200), NACUM(200)
      DIMENSION LVL(1), LVLS1(1), LVLS2(1)
      DO 10 I=1,IDPTH
        NACUM(I) = 0
   10 CONTINUE
      DO 30 I=1,LCON
        LVL(I) = 1
        LVLS2(I) = IDPTH + 1 - LVLS2(I)
        ITEMP = LVLS2(I)
        IF (ITEMP.GT.IDPTH) GOTO 30
        IF (ITEMP.NE.LVLS1(I)) GOTO 20
        NACUM(ITEMP) = NACUM(ITEMP) + 1
        GOTO 30
   20   LVL(I) = 0
   30 CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      INTEGER FUNCTION SORT2(DMY)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C SORT2 SORTS SIZE1 AND STPT INTO DESCENDING ORDER ACCORDING TO
C VALUES OF SIZE1. XC=NUMBER OF ENTRIES IN EACH ARRAY
      INTEGER*4 TEMP, SIZE1, STPT, XC
C IT IS ASSUMED THAT THE GRAPH HAS AT MOST 50 CONNECTED COMPONENTS.
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level)
      COMMON /CC/ XC, SIZE1(100), STPT(100)
      SORT2 = 0
      IF (XC.EQ.0) RETURN
      SORT2 = 1
      IND = XC
   10 ITEST = 0
      IND = IND - 1
      IF (IND.LT.1) RETURN
      DO 20 I=1,IND
        J = I + 1
        IF (SIZE1(I).GE.SIZE1(J)) GOTO 20
        ITEST = 1
        TEMP = SIZE1(I)
        SIZE1(I) = SIZE1(J)
        SIZE1(J) = TEMP
        TEMP = STPT(I)
        STPT(I) = STPT(J)
        STPT(J) = TEMP
   20 CONTINUE
      IF (ITEST.EQ.1) GOTO 10
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE PIKLVL(LVLS1, LVLS2, CCSTOR, IDFLT, ISDIR)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C PIKLVL CHOOSES THE LEVEL STRUCTURE  USED IN NUMBERING GRAPH
C LVLS1-       ON INPUT CONTAINS FORWARD LEVELING INFO
C LVLS2-       ON INPUT CONTAINS REVERSE LEVELING INFO
C              ON OUTPUT THE FINAL LEVEL STRUCTURE CHOSEN
C CCSTOR-      ON INPUT CONTAINS CONNECTED COMPONENT INFO
C IDFLT-       ON INPUT =1 IF WDTH LVLS1.LE.WDTH LVLS2, =2 OTHERWISE
C NHIGH        KEEPS TRACK OF LEVEL WIDTHS FOR HIGH NUMBERING
C NLOW-        KEEPS TRACK OF LEVEL WIDTHS FOR LOW NUMBERING
C NACUM-       KEEPS TRACK OF LEVEL WIDTHS FOR CHOSEN LEVEL STRUCTURE
C XC-          NUMBER OF CONNECTED COMPONENTS
C SIZE1(I)-    SIZE OF I'TH CONNECTED COMPONENT
C STPT(I)-     INDEX INTO CCSTORE OF 1ST NODE IN I'TH CON COMPT
C ISDIR-       FLAG WHICH INDICATES WHICH WAY THE LARGEST CONNECTED
C              COMPONENT FELL.  =+1 IF LOW AND -1 IF HIGH
      INTEGER*4 CCSTOR, SIZE1, STPT, XC, END1
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C IT IS ASSUMED THAT THE GRAPH HAS AT MOST 50 COMPONENTS AND
C THAT THERE ARE AT MOST 100 LEVELS.
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level)
      COMMON /LVLW/ NHIGH(200), NLOW(200), NACUM(200)
      COMMON /CC/ XC, SIZE1(100), STPT(100)
      DIMENSION LVLS1(1), LVLS2(1), CCSTOR(1)
C FOR EACH CONNECTED COMPONENT DO
      DO 80 I=1,XC
        J = STPT(I)
        END1 = SIZE1(I) + J - 1
C SET NHIGH AND NLOW EQUAL TO NACUM
        DO 10 K=1,IDPTH
          NHIGH(K) = NACUM(K)
          NLOW(K) = NACUM(K)
   10   CONTINUE
C UPDATE NHIGH AND NLOW FOR EACH NODE IN CONNECTED COMPONENT
        DO 20 K=J,END1
          INODE = CCSTOR(K)
          LVLNH = LVLS1(INODE)
          NHIGH(LVLNH) = NHIGH(LVLNH) + 1
          LVLNL = LVLS2(INODE)
          NLOW(LVLNL) = NLOW(LVLNL) + 1
   20   CONTINUE
        MAX1 = 0
        MAX2 = 0
C SET MAX1=LARGEST NEW NUMBER IN NHIGH
C SET MAX2=LARGEST NEW NUMBER IN NLOW
        DO 30 K=1,IDPTH
          IF (2*NACUM(K).EQ.NLOW(K)+NHIGH(K)) GOTO 30
          IF (NHIGH(K).GT.MAX1) MAX1 = NHIGH(K)
          IF (NLOW(K).GT.MAX2) MAX2 = NLOW(K)
   30   CONTINUE
C SET IT= NUMBER OF LEVEL STRUCTURE TO BE USED
        IT = 1
        IF (MAX1.GT.MAX2) IT = 2
        IF (MAX1.EQ.MAX2) IT = IDFLT
        IF (IT.EQ.2) GOTO 60
        IF (I.EQ.1) ISDIR = -1
C COPY LVLS1 INTO LVLS2 FOR EACH NODE IN CONNECTED COMPONENT
        DO 40 K=J,END1
          INODE = CCSTOR(K)
          LVLS2(INODE) = LVLS1(INODE)
   40   CONTINUE
C UPDATE NACUM TO BE THE SAME AS NHIGH
        DO 50 K=1,IDPTH
          NACUM(K) = NHIGH(K)
   50   CONTINUE
        GOTO 80
C UPDATE NACUM TO BE THE SAME AS NLOW
   60   DO 70 K=1,IDPTH
          NACUM(K) = NLOW(K)
   70   CONTINUE
   80 CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE NUMBER(SND, NUM, NDSTK, LVLS2, NDEG, RENUM, LVLST,
     * LSTPT, NR, NFLG, IPFA, ISDIR)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C  NUMBER PRODUCES THE NUMBERING OF THE GRAPH FOR MIN BANDWIDTH
C  SND-         ON INPUT THE NODE TO BEGIN NUMBERING ON
C  NUM-         ON INPUT AND OUTPUT, THE NEXT AVAILABLE NUMBER
C  LVLS2-       THE LEVEL STRUCTURE TO BE USED IN NUMBERING
C  RENUM-       THE ARRAY USED TO STORE THE NEW NUMBERING
C  LVLST-       ON OUTPUT CONTAINS LEVEL STRUCTURE
C  LSTPT(I)-    ON OUTPUT, INDEX INTO LVLST TO FIRST NODE IN ITH LVL
C               LSTPT(I+1) - LSTPT(I) = NUMBER OF NODES IN ITH LVL
C  NFLG-        =+1 IF SND IS FORWARD END OF PSEUDO-DIAM
C               =-1 IF SND IS REVERSE END OF PSEUDO-DIAM
C  IBW2-        BANDWIDTH OF NEW NUMBERING COMPUTED BY NUMBER
C  IPF2-        PROFILE OF NEW NUMBERING COMPUTED BY NUMBER
C  IPFA-        WORKING STORAGE USED TO COMPUTE PROFILE AND BANDWIDTH
C  ISDIR-       INDICATES STEP DIRECTION USED IN NUMBERING(+1 OR -1)
C USE INTEGER*2 NDSTK  WITH AN IBM 360 OR 370.
      INTEGER*4 NDSTK
      INTEGER*4 SND, STKA, STKB, STKC, STKD, XA, XB, XC, XD, CX, END1,
     * RENUM, TEST
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C THE STORAGE IN COMMON BLOCKS CC AND LVLW IS NOW FREE AND CAN
C BE USED FOR STACKS.
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level)
      COMMON /LVLW/ STKA(200), STKB(200), STKC(200)
C	SUNIL 18/10/10 CHANGE COMMON BLOCK CC TO CC1 SINCE CC HAVING
C	A DIFFERENT NUMBER OF ARGUMENTS
C      COMMON /CC/ STKD(100)
      COMMON /CC1/ STKD(200)
      DIMENSION IPFA(1)
      DIMENSION NDSTK(NR,1), LVLS2(1), NDEG(1), RENUM(1), LVLST(1),
     * LSTPT(1)
C SET UP LVLST AND LSTPT FROM LVLS2
      DO 10 I=1,LCON
        IPFA(I) = 0
   10 CONTINUE
      NSTPT = 1
      DO 30 I=1,IDPTH
        LSTPT(I) = NSTPT
        DO 20 J=1,LCON
          IF (LVLS2(J).NE.I) GOTO 20
          LVLST(NSTPT) = J
          NSTPT = NSTPT + 1
   20   CONTINUE
   30 CONTINUE
      LSTPT(IDPTH+1) = NSTPT
C STKA, STKB, STKC AND STKD ARE STACKS WITH POINTERS
C XA,XB,XC, AND XD.  CX IS A SPECIAL POINTER INTO STKC WHICH
C INDICATES THE PARTICULAR NODE BEING PROCESSED.
C LVLN KEEPS TRACK OF THE LEVEL WE ARE WORKING AT.
C INITIALLY STKC CONTAINS ONLY THE INITIAL NODE, SND.
      LVLN = 0
      IF (NFLG.LT.0) LVLN = IDPTH + 1
      XC = 1
      STKC(XC) = SND
   40 CX = 1
      XD = 0
      LVLN = LVLN + NFLG
      LST = LSTPT(LVLN)
      LND = LSTPT(LVLN+1) - 1
C BEGIN PROCESSING NODE STKC(CX)
   50 IPRO = STKC(CX)
      RENUM(IPRO) = NUM
      NUM = NUM + ISDIR
      END1 = NDEG(IPRO)
      XA = 0
      XB = 0
C CHECK ALL ADJACENT NODES
      DO 80 I=1,END1
        TEST = NDSTK(IPRO,I)
        INX = RENUM(TEST)
C ONLY NODES NOT NUMBERED OR ALREADY ON A STACK ARE ADDED
        IF (INX.EQ.0) GOTO 60
        IF (INX.LT.0) GOTO 80
C DO PRELIMINARY BANDWIDTH AND PROFILE CALCULATIONS
        NBW = (RENUM(IPRO)-INX)*ISDIR
        IF (ISDIR.GT.0) INX = RENUM(IPRO)
        IF (IPFA(INX).LT.NBW) IPFA(INX) = NBW
        GOTO 80
   60   RENUM(TEST) = -1
C PUT NODES ON SAME LEVEL ON STKA, ALL OTHERS ON STKB
        IF (LVLS2(TEST).EQ.LVLS2(IPRO)) GOTO 70
        XB = XB + 1
        STKB(XB) = TEST
        GOTO 80
   70   XA = XA + 1
        STKA(XA) = TEST
   80 CONTINUE
C SORT STKA AND STKB INTO INCREASING DEGREE AND ADD STKA TO STKC
C AND STKB TO STKD
      IF (XA.EQ.0) GOTO 100
      IF (XA.EQ.1) GOTO 90
      CALL SORTDG(STKC, STKA, XC, XA, NDEG)
      GOTO 100
   90 XC = XC + 1
      STKC(XC) = STKA(XA)
  100 IF (XB.EQ.0) GOTO 120
      IF (XB.EQ.1) GOTO 110
      CALL SORTDG(STKD, STKB, XD, XB, NDEG)
      GOTO 120
  110 XD = XD + 1
      STKD(XD) = STKB(XB)
C BE SURE TO PROCESS ALL NODES IN STKC
  120 CX = CX + 1
      IF (XC.GE.CX) GOTO 50
C WHEN STKC IS EXHAUSTED LOOK FOR MIN DEGREE NODE IN SAME LEVEL
C WHICH HAS NOT BEEN PROCESSED
      MAX = MCON + 1
      SND = LCON + 1
      DO 130 I=LST,LND
        TEST = LVLST(I)
        IF (RENUM(TEST).NE.0) GOTO 130
        IF (NDEG(TEST).GE.MAX) GOTO 130
        RENUM(SND) = 0
        RENUM(TEST) = -1
        MAX = NDEG(TEST)
        SND = TEST
  130 CONTINUE
      IF (SND.EQ.LCON+1) GOTO 140
      XC = XC + 1
      STKC(XC) = SND
      GOTO 50
C IF STKD IS EMPTY WE ARE DONE, OTHERWISE COPY STKD ONTO STKC
C AND BEGIN PROCESSING NEW STKC
  140 IF (XD.EQ.0) GOTO 160
      DO 150 I=1,XD
        STKC(I) = STKD(I)
  150 CONTINUE
      XC = XD
      GOTO 40
C DO FINAL BANDWIDTH AND PROFILE CALCULATIONS
  160 DO 170 I=1,LCON
        IF (IPFA(I).GT.IBW2) IBW2 = IPFA(I)
        IPF2 = IPF2 + IPFA(I)
  170 CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE PROFIT(NR, NDSTK, NEW, NDEG, LVLS2, LVLST, LSTPT,
     * NXTNUM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C SUBROUTINE PROFIT NUMBERS LEVEL BY LEVEL WITH CONSECUTIVE INTEGERS
C USING A MODIFIED VERSION OF KING*S ALGORITHM.
C   NR-        ROW DIMENSION OF CONNECTION TABLE.
C   NDSTK-     THE CONNECTION TABLE.
C   NEW-       VECTOR TO STORE THE NEW NUMBERING.
C   NDEG(I)-   THE DEGREE OF NODE I.
C   LVLS2-     THE LEVEL STRUCTURE PRODUCED BY PIKLVL.
C              LVLS2(I) = J MEANS VERTEX I HAS BEEN
C              PLACED IN LEVEL J.
C   LVLST-     ON OUTPUT,  CONTAINS THE LEVEL STRUCTURE USED.
C              LVLST(LSTPT(I)),...,LVLST(LSTPT(I+1)-1) ARE
C              THE VERTICES IN LEVEL I.
C   LSTPT(I)-  ON OUTPUT, INDEX INTO LVLST TO FIRST NODE IN LEVEL I.
C              LSTPT(I+1) - LSTPT(I) = NUMBER OF NODES IN I*TH LEVEL.
C   NXTNUM-    ON INPUT AND OUTPUT,  THE NEXT AVAILABLE NUMBER.
C ON IBM 360 OR 370 USE INTEGER * 2 NDSTK.
      INTEGER*4 NDSTK
      INTEGER*4 NEW(1), NDEG(1), LVLS2(1), LVLST(1), LSTPT(1)
      DIMENSION NDSTK(NR,1)
C COMMON AREA GRA HOLDS VITAL INFORMATION ABOUT THE GRAPH
C N-         THE NUMBER OF NODES
C IDPTH-     THE NUMBER OF LEVELS FOUND BY PIKLVL.
C MCON-      MAXIMUM DEGREE OF GRAPH -- COLUMN DIMENSION OF NDSTK.
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
C IT IS ASSUMED THAT NO LEVEL HAS MORE THAN 100 NODES.
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level)
      COMMON /LVLW/ S2(200), S3(200), Q(200)
C	SUNIL 18/10/10 CHANGE COMMON BLOCK CC TO CC1 SINCE CC HAVING
C	A DIFFERENT NUMBER OF ARGUMENTS
C      COMMON /CC/ CONECT(100)
      COMMON /CC1/ CONECT(200)
      INTEGER*4 S2, S3, Q, CONECT, S2SZE, S3SZE, QPTR, CONSZE
C SET UP LVLST AND LSTPT FROM LVLS2.
      NSTPT = 1
      DO 20 I=1,IDPTH
        LSTPT(I) = NSTPT
        DO 10 J=1,LCON
          IF (LVLS2(J).NE.I) GOTO 10
          LVLST(NSTPT) = J
          NSTPT = NSTPT + 1
   10   CONTINUE
   20 CONTINUE
      LSTPT(IDPTH+1) = NSTPT
C ******************  STEP P0   ****************************************
C  S2 IS THE FIRST LEVEL.
      LEVEL = 1
      CALL FORMLV(S2, S2SZE, LSTPT, LVLST, LEVEL)
C ******************  STEP P1   ****************************************
C  S3 IS THE LEVEL ADJACENT TO THE LEVEL S2.
C  Q IS A QUEUE USED TO RETAIN THE ORDER IN WHICH THE ELEMENTS OF S3
C  ARE REMOVED.  Q EVENTUALLY BECOMES THE NEW S2 AND IS ORDERED
C  ACCORDING TO KING*S CRITERA.
   30 CALL FORMLV(S3, S3SZE, LSTPT, LVLST, LEVEL+1)
      QPTR = 0
C ******************  STEP P2   ****************************************
C  FIND THE NODE M IN S2 WHICH IS ADJACENT TO THE FEWEST NODES IN S3.
   40 M = MINCON(S2,S2SZE,S3,S3SZE,CONECT,CONSZE,NDSTK,NR,NDEG)
C ******************  STEP P3   ****************************************
C  NUMBER M AND REMOVE IT FROM S2.
      NEW(M) = NXTNUM
      NXTNUM = NXTNUM + 1
      CALL DELETE(S2, S2SZE, M)
      IF (CONSZE.LE.0) GOTO 60
C THE ELEMENTS OF CONLST ARE TO BE REMOVED FROM S3 AND PLACED INTO
C Q.
      DO 50 I=1,CONSZE
        QPTR = QPTR + 1
        Q(QPTR) = CONECT(I)
        CALL DELETE(S3, S3SZE, CONECT(I))
   50 CONTINUE
C ******************  STEP P4   ****************************************
   60 IF (S2SZE.LE.0) GOTO 80
C ******************  STEP P5   ****************************************
      IF (S3SZE.GT.0) GOTO 40
C ******************  STEP P6   ****************************************
C  S3 IS EMPTY, BUT S2 IS NOT.  RENUMBER THE NODES WHICH REMAIN IN S2.
      DO 70 I=1,S2SZE
        NS2 = S2(I)
        NEW(NS2) = NXTNUM
        NXTNUM = NXTNUM + 1
   70 CONTINUE
      GOTO 100
C ******************  STEP P7   ****************************************
   80 IF (S3SZE.LE.0) GOTO 100
C S2 IS EMPTY, BUT S3 IS NOT.  MOVE S3*S REMAINING NODES INTO Q.
      DO 90 I=1,S3SZE
        QPTR = QPTR + 1
        Q(QPTR) = S3(I)
   90 CONTINUE
C ******************  STEP P8   ****************************************
  100 LEVEL = LEVEL + 1
      IF (LEVEL.GE.IDPTH) GOTO 120
C S2 BECOMES THE OLD Q SINCE BOTH S2 AND S3 ARE EMPTY.
      DO 110 I=1,QPTR
        S2(I) = Q(I)
  110 CONTINUE
      S2SZE = QPTR
      GOTO 30
C ******************  STEP P9   ****************************************
C  LAST LEVEL IS ORDERED IN Q,  SO NUMBER IT BEFORE RETURNING.
  120 DO 130 I=1,QPTR
        IQ = Q(I)
        NEW(IQ) = NXTNUM
        NXTNUM = NXTNUM + 1
  130 CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      FUNCTION MINCON(X, XSZE, Y, YSZE, CONLST, CONSZE, NDSTK, NR,
     * NDEG)
	IMPLICIT REAL*8 (A-H,O-Z)
C FUNCTION MINCON RETURNS AS ITS FUNCTIONAL VALUE A VERTEX X(I) SUCH
C THAT THE NUMBER OF CONNECTIONS FROM X(I) TO THE SET Y IS A MINIMUM.
C THE VERTICES OF Y WHICH ARE ADJACENT TO X(I) ARE PLACED IN
C CONLST(1), CONLST(2),...,CONLST(CONSZE).
C USE INTEGER * 2 NDSTK ON IBM 360 OR 370.
      INTEGER*4 NDSTK
      DIMENSION NDSTK(NR,1)
      INTEGER*4 X(1), XSZE, Y(1), YSZE, CONLST(1), CONSZE, NDEG(1)
C IT IS ASSUMED THAT NO LEVEL HAS MORE THAN 100 VERTICES.
c	sunil 01/04/23 change sizes of /cc/ to 100 and /lvlw/ to 200
c	(making the graph has at most 100 connected components and
c	 at most 200 nodes at the last level), due to change of above
c	 two common blocks smlst(100) change to smlst(200)
      INTEGER*4 SMLST(200)
      CONSZE = YSZE + 1
      DO 50 I=1,XSZE
        LSTSZE = 0
        IX = X(I)
        IROWDG = NDEG(IX)
        DO 20 J=1,YSZE
          DO 10 K=1,IROWDG
            IX = X(I)
            IF (NDSTK(IX,K).NE.Y(J)) GOTO 10
            SMLST(LSTSZE+1) = Y(J)
            LSTSZE = LSTSZE + 1
            IF (LSTSZE.GE.CONSZE) GOTO 50
            GOTO 20
   10     CONTINUE
   20   CONTINUE
        IF (LSTSZE.GT.0) GOTO 30
C WE HAVE FOUND A VERTEX IN X WHICH IS NOT CONNECTED TO ANY VERTEX
C IN Y
        MINCON = X(I)
        CONSZE = 0
        RETURN
C WE HAVE FOUND A VERTEX X(I) WITH FEWEST CONNECTIONS (NONZERO) TO Y
C SO FAR.  SAVE THE ELEMENTS OF Y WHICH CONNECT TO X(I) IN CONLST AND
C SAVE X(I) AS THE FUNCTIONAL VALUE.
   30   CONSZE = LSTSZE
        DO 40 J=1,LSTSZE
          CONLST(J) = SMLST(J)
   40   CONTINUE
        MINCON = X(I)
   50 CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE DELETE(SET, SETSZE, ELEMNT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SUBROUTINE DELETE REMOVES ELEMNT FROM THE SET SET IF ELEMNT
C	IS IN SET.  OTHERWISE,  IT ISSUES A DIAGNOSTIC.
      INTEGER*4 SET(1), SETSZE, ELEMNT
      IF (SETSZE.GT.1) GOTO 10
      IF (SETSZE.EQ.1 .AND. SET(1).NE.ELEMNT) GOTO 30
      SETSZE = 0
      RETURN
   10 DO 20 I=1,SETSZE
        IF (SET(I).EQ.ELEMNT) GOTO 40
   20 CONTINUE
   30 WRITE (6,99999) ELEMNT, (SET(I),I=1,SETSZE)
      RETURN
   40 SETSZE = SETSZE - 1
      DO 50 J=I,SETSZE
        SET(J) = SET(J+1)
   50 CONTINUE
      RETURN
99999 FORMAT (10H0ERROR -- , I6, 8H NOT IN , (20I5))
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE FORMLV(SET, SETSZE, LSTPT, LVLST, LEVEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	FORMLVL COPIES LEVEL(LEVEL) INTO SET.
      INTEGER*4 SET(1), SETSZE, LSTPT(1), LVLST(1), UPPER
      LOWER = LSTPT(LEVEL)
      UPPER = LSTPT(LEVEL+1) - 1
      SETSZE = 1
      DO 10 I=LOWER,UPPER
        SET(SETSZE) = LVLST(I)
        SETSZE = SETSZE + 1
   10 CONTINUE
      SETSZE = SETSZE - 1
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE CHECK(BESTBW, BESTPF, RENUM, NDSTK, NR, NDEG, IWK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SUBROUTINE CHECK TESTS TO SEE IF REVERSED NUMBERING GIVES BETTER
C	PROFILE THAN PROFIT.  IF IT DOES, THEN RENUM IS REVERSED AND BESTPF
C	IS SET TO THE SMALLEST OF RENUM AND REVERSED RENUM.
C	USE INTEGER * 2 NDSTK ON IBM 360 OR 370
      INTEGER*4 NDSTK
      DIMENSION NDSTK(NR,1)
      INTEGER*4 BESTBW, BESTPF, RENUM(1), NDEG(1), IWK(1)
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
      DO 10 I=1,LCON
        IWK(I) = LCON - RENUM(I) + 1
   10 CONTINUE
      CALL BAND(BESTBW, BESTPF, RENUM, NDSTK, NR, NDEG)
      CALL BAND(IBW, IPF, IWK, NDSTK, NR, NDEG)
      IF (IPF.GE.BESTPF) RETURN
      DO 20 I=1,LCON
        RENUM(I) = IWK(I)
   20 CONTINUE
      BESTPF = IPF
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE BAND(IBW, IPF, NEW, NDSTK, NR, NDEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SUBROUTINE BAND COMPUTES THE BANDWIDTH IBW AND THE PROFILE
C	IPF OF THE GRAPH REPRESENTED BY NDSTK USING THE NUMBERING NEW.
C	ON IBM 360 OR 370 USE INTEGER * 2 NDSTK.
      INTEGER*4 NDSTK
      DIMENSION NDSTK(NR,1), NEW(1), NDEG(1)
      COMMON /BOPT/  IBW1,IPF1,IBW2,IPF2,LCON,MCON,IDPTH,MAXCON,IREDU
      IPF = 0
      IBW = 0
      DO 20 K=1,LCON
        IEND = NDEG(K)
        IF (IEND.EQ.0) GOTO 20
        NBW = 0
        DO 10 J=1,IEND
          IDUMMY = NDSTK(K,J)
          NTEST = NEW(K) - NEW(IDUMMY)
          IF (NTEST.LE.NBW) GOTO 10
          NBW = NTEST
   10   CONTINUE
        IPF = IPF + NBW
        IF (NBW.GT.IBW) IBW = NBW
   20 CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE LOCRES (ID,IDSET,DIRCOS,LM,ESTI,EDIS,EDISI,ELOD,
     1                   MSF,MNF,MODE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIt INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     SCANS NODES TO DETERMINE WHICH HAVE LOCAL (SKEWED) COORDINATES
C     AND CALLS ROUTINES TO PERFORM APPROPRIATE TRANSFORMATIONS
C	---------------------------------------------------------
C     MODE = 1  TRANSFORMS ELEMENT DISPLACEMENT VECTOR INTO GLOBAL
C     MODE = 2  TRANSFORMS ELEMENT LOAD VECTOR INTO LOCAL  AND
C        IF IFREF = 0  TRANSFORMS ELEMENT STIFFNESS MATRIX INTO GLOBAL
C     MODE = 3  TRANSFORMS NODAL DISPLACEMENT VECTOR INTO GLOBAL
C               PRIOR TO OUTPUT
C     MODE = 4  TRANSFORMS APPLIED TRACTION AND PRESSURE LOADS
C               INTO LOCAL AT SKEW NODES
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C      COMMON /ITER/ NSTEP,DLMAX,NPRIN,NDRAW,KONEQ,NIREF,ITEMAX,RTOL,
C     1              ITOPT,NUMREF,NUMITE,ITETOT,RHO,RHOP,ICONV,NOLIN,
C     2              LIMEQ(2),KSTEP,ETOL,RHOPREV
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON A(9000000),IA(9000000)

C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)
C
      DIMENSION ID(MSF,1),IDSET(1),DIRCOS(9,1),LM(MNF,1),ESTI(1)
      DIMENSION EDIS(1),EDISI(1),ELOD(1)

C	DIMENSION COSN(6,6),TCOS(6,6),S(6,6),ST(6,6)
	DIMENSION COSN(81),TCOS(81),S(81),ST(81)

      DIMENSION IGPOS(9),IDPOS(9)
C
C
      IF (ITASK.EQ.5) RETURN
      IF ((NOLIN.EQ.0).AND.(ITASK.NE.3).AND.(MODE.EQ.1)) RETURN
C
      LNF = MNF

      IF (MODE.EQ.4) LNF = 3
C     ----------------------------------------------------------------------------
C     TRANSFORM NODAL DISPLACEMENT VECTOR FOR A NODE D(1-9) INTO GLOBAL FOR OUTPUT
C     ----------------------------------------------------------------------------
      IF (MODE.NE.3) GOTO 20
      
      ISN = MNF
      ISET = IDSET(ISN)
      IF (ISET.EQ.0) RETURN
      DO 40  I=1,6
  40  IGPOS(I) = I
      CALL COSMAT (DIRCOS(1,ISET),COSN,IGPOS,6)
      CALL MATRAN (COSN,TCOS,6,6,1)
      CALL TREVEC (EDIS(1),TCOS,6)
      RETURN
      
C     ----------------------------------------------------
C     IGPOS = GLOBAL D.O.F. NOS. POSSIBLE FOR THIS ELEMENT
C     IDPOS = POSITIONS OF THESE D.O.F.'S IN ID
C     ----------------------------------------------------
  20  LINK = NFLINK(ITYPE,ISTYP)
      IGF = 0
      DO 300  INF = 1,NNF
      IGPOS(INF) = LINK/10**(NNF-INF)
 320  IGF = IGF+1
	IF(IGF.GT.9) GOTO 301
      IF (IDOF(IGF)-IGPOS(INF)) 320,340,340
 340  LINK = LINK-10**(NNF-INF)*IGPOS(INF)
 300  IDPOS(INF) = IGF
 301	CONTINUE

      LINK = NFLINK(ITYPE,ISTYP)
      DO INF=1,NNF
      IGPOS(INF) = LINK/10**(NNF-INF)
      IF ((IGPOS(INF).GT.6)) LNF = LNF-1    !MODIFY LNF HERE IF DOF NOT SUBJECT TO TRANSFORMATION
	LINK = LINK-10**(NNF-INF)*IGPOS(INF)
	ENDDO
C     ---------------------------------------------------
C     SCAN LM FOR FIRST EQN. NO. FOR EACH NODE OF ELEMENT
C     ---------------------------------------------------
      DO 100  INO = 1,NNO
      DO 120  INF = 1,NNF
      IF (LM(INF,INO).EQ.0) GOTO 120
      IEQ   = LM(INF,INO)
	IGF   = IGPOS(INF)
	IDROW = IDOFCALL(IDOF,IGF)
      GOTO 140
 120  CONTINUE
	GOTO 100
	
C     --------------------------------
C     SCAN ID FOR AXES SET NO. OF NODE
C     --------------------------------
 140  DO 200  ISN = 1,NSN
      IF (ID(IDROW,ISN).EQ.IEQ) GOTO 220
 200  CONTINUE
 220  ISET = IDSET(ISN)
      IF (ISET.EQ.0) GOTO 100
C     -------------------------------------------
C     EXTRACT APPROPRIATE DIRECTION COSINE MATRIX
C     -------------------------------------------
      CALL COSMAT (DIRCOS(1,ISET),COSN,IGPOS,LNF)
C
      NADD = (INO-1)*MNF+1
      IF (MODE-2) 400,500,600
C     ----------------------------------------------------------
C     TRANSFORM ELEMENT DISPLACEMENT VECTOR INTO GLOBAL (MODE=1)
C     ----------------------------------------------------------
 400  CALL MATRAN (COSN,TCOS,LNF,LNF,1)
      CALL TREVEC (EDIS(NADD),TCOS,LNF)
      GOTO 100
      
C     ------------------------------------------------------
C     TRANSFORM ELEMENT LOAD VECTOR INTO GLOBAL AND
C     ELEMENT STIFFNESS OR INITIAL STRESS MATRIX IF REQUIRED
C     ------------------------------------------------------
 500  CONTINUE
	CALL TREVEC (ELOD(NADD),COSN,LNF)
      IF (ITASK.GE.3) GOTO 520
      IF (IFREF.NE.0) GOTO 100
      GOTO 540
 520  IF (ITASK.EQ.3 .AND. IFEIG.NE.0) GOTO 100
 540  CALL MATRAN (COSN,TCOS,LNF,LNF,1)
      CALL TRESTI (ESTI,COSN,TCOS,MNF,LNF,NNO,INO,S,ST)

	IF(LSYMM.EQ.1) THEN
	NEF2 = (NEF*NEF+NEF)/2 + 1 
	CALL TRESTI (ESTI(NEF2),COSN,TCOS,MNF,LNF,NNO,INO,S,ST)
	ENDIF

      GOTO 100
C     ------------------------------------
C     TRANSFORM EXTERNAL DISTRIBUTED LOADS
C     INTO LOCAL AT SKEW NODES (MODE=4)
C     ------------------------------------
 600  IF(MODE.EQ.5) GOTO 700
	NADD = (INO-1)*LNF + 1
      IF (NOLIN.NE.0) CALL TREVEC (ELOD(NADD),COSN,LNF)
	GOTO 100

 700  NADD = (INO-1)*MNF + 1
      CALL TREVEC (ELOD(NADD),COSN,LNF)
C
 100  CONTINUE
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE TRESTI_OLD (ESTI,COSN,TCOS,MNF,LNF,NNO,INO,S,ST)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------
C     TRANSFORM ELEMENT STIFFNESS MATRIX TO LOCAL AXES
C     ------------------------------------------------
      DIMENSION ESTI(1),COSN(LNF,1),TCOS(LNF,1)
      DIMENSION S(LNF,1),ST(LNF,1)
      NEF = MNF*NNO
      LNF2 = LNF*LNF
C     ------------------------------------------------
C     EXTRACT AND TRANSFORM LEADING DIAGONAL SUBMATRIX
C     ------------------------------------------------
      DO 100  J = 1,LNF
      INF = (INO-1)*MNF+J
      IJ = LNF+1-J
      DO 100  I = 1,IJ
      IPOS = NEF*(INF-1)-(INF-1)*(INF-2)/2+I
      II = I+J-1
      S(J,II) = ESTI(IPOS)
      IF (II.GT.J)  S(II,J) = S(J,II)
 100  CONTINUE
      CALL CLEARA (ST,LNF2)
      CALL MATMULT (S,TCOS,ST,LNF,LNF,LNF,0.,1)
      CALL CLEARA (S,LNF2)
      CALL MATMULT (COSN,ST,S,LNF,LNF,LNF,0.,1)
      DO 120  J = 1,LNF
      INF = (INO-1)*MNF+J
      IJ = LNF+1-J
      DO 120  I = 1,IJ
      IPOS = NEF*(INF-1)-(INF-1)*(INF-2)/2+I
      II = I+J-1
 120  ESTI(IPOS) = S(J,II)
      IF (INO.EQ.1) GOTO 300
C     -------------------------------------------
C     EXTRACT AND TRANSFORM SUBMATRICES UP COLUMN
C     -------------------------------------------
      INO1 = INO-1
      DO 200  JNO = 1,INO1
      DO 220  J = 1,LNF
      JJ = (JNO-1)*MNF+J
      MPOS = 1+NEF*(JJ-1)-(JJ-1)*(JJ-2)/2+(INO-JNO)*MNF-J+1
      DO 220  I = 1,LNF
      IPOS = MPOS-1+I
 220  S(J,I) = ESTI(IPOS)
      CALL CLEARA (ST,LNF2)
      CALL MATMULT (S,TCOS,ST,LNF,LNF,LNF,0.,1)
      DO 240  J = 1,LNF
      JJ = (JNO-1)*MNF+J
      MPOS = 1+NEF*(JJ-1)-(JJ-1)*(JJ-2)/2+(INO-JNO)*MNF-J+1
      DO 240  I = 1,LNF
      IPOS = MPOS-1+I
 240  ESTI(IPOS) = ST(J,I)
 200  CONTINUE
 300  IF (INO.EQ.NNO) RETURN
C     -------------------------------------------
C     EXTRACT AND TRANSFORM SUBMATRICES ALONG ROW
C     -------------------------------------------
      NNO1 = NNO-INO
      DO 400  JNO = 1,NNO1
      DO 420  J = 1,LNF
      INF = (INO-1)*MNF+J
      MPOS = 1+NEF*(INF-1)-(INF-1)*(INF-2)/2+MNF*JNO+1-J
      DO 420  I = 1,LNF
      IPOS = MPOS-1+I
 420  S(J,I) = ESTI (IPOS)
      CALL CLEARA (ST,LNF2)
      CALL MATMULT (COSN,S,ST,LNF,LNF,LNF,0.,1)
      DO 440  J = 1,LNF
      INF = (INO-1)*MNF+J
      MPOS = 1+NEF*(INF-1)-(INF-1)*(INF-2)/2+MNF*JNO+1-J
      DO 440  I = 1,LNF
      IPOS = MPOS-1+I
 440  ESTI(IPOS) = ST(J,I)
 400  CONTINUE
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================

      SUBROUTINE TRESTI (ESTI,COSN,TCOS,MNF,LNF,NNO,INO,S,ST)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------
C     TRANSFORM ELEMENT STIFFNESS MATRIX TO LOCAL AXES
C     ------------------------------------------------
      DIMENSION ESTI(1),COSN(LNF,1),TCOS(LNF,1)
      DIMENSION S(LNF,1),ST(LNF,1)
      
      ALLOCATABLE SS(:,:),TMAT(:,:)
      
      NEF = MNF*NNO
      
      ALLOCATE(SS(NEF,NEF),TMAT(NEF,NEF))
      
      TMAT(1:NEF,1:NEF) = 0.0D0
      
      KEF = 0
      DO IEF = 1,NEF
      TMAT(IEF,IEF) = 1.0D0
      DO JEF = IEF,NEF
      KEF = KEF + 1
      SS(IEF,JEF) = ESTI(KEF)
      IF(JEF.GT.IEF) SS(JEF,IEF) = SS(IEF,JEF)
      ENDDO
      ENDDO
      
      NMF = MNF*(INO-1)
      DO INF = 1,LNF
      DO JNF = 1,LNF
      TMAT(INF+NMF,JNF+NMF) = TCOS(INF,JNF)
      ENDDO
      ENDDO
      
      SS = MATMUL(TRANSPOSE(TMAT),MATMUL(SS,TMAT))
      
      KEF = 0
      DO IEF = 1,NEF
      DO JEF = IEF,NEF
      KEF = KEF + 1
      ESTI(KEF) = SS(IEF,JEF)
      ENDDO
      ENDDO

      DEALLOCATE(SS,TMAT)
      
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================

      SUBROUTINE TREVEC (EVEC,COSN,LNF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------
C     TRANSFORM NODAL SUBVECTOR
C     -------------------------
      DIMENSION EVEC(1),COSN(LNF,1),EL(6)
      CALL CLEARA (EL,LNF)
      CALL MATMULT (COSN,EVEC,EL,LNF,LNF,1,0.,1)
      DO 100  I = 1,LNF
 100  EVEC(I) = EL(I)
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE COSMAT (DIRCOS,COSN,IGPOS,LNF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------
C     EXTRACT APPROPRIATE TERMS FROM DIRECTION COSINE MATRIX
C     ------------------------------------------------------
      DIMENSION DIRCOS(1),IGPOS(1),COSN(LNF,1)
      DO 110  I = 1,LNF
      N = IGPOS(I)
      DO 110  J = 1,LNF
      M = IGPOS(J)
      IF ((N.LT.4).AND.(M.LT.4)) GOTO 130
      IF ((N.GE.4).AND.(M.GE.4)) GOTO 140
      COSN(N,M) = 0.
      GOTO 110
 130  NM = (N-1)*3+M
      COSN(N,M) = DIRCOS(NM)
      GOTO 110
 140  NM = (N-4)*3+M-3
      COSN(N,M) = DIRCOS(NM)
 110  CONTINUE
      RETURN
      END
C
C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE MODULE (PROPM,PROPG,PROPO,NODEX,WA,WA2,AMV,S,
     +		COORD,EDIS,EDISI,ELOD,ALPHA,SEL,SEDI,FIN,
     +		IPIN,ISET,MSET,HINFC,FIXEN,FIXLR,FIXEO,FIXLO,
     +		CABFF,SEQ,HINFP,TCHET,TAMBT,PMATRL,ISFAC,DISLI,
     +        FIXEN_OFF,FIXLR_OFF,FIXEO_OFF,FIXLO_OFF)


	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
C     ------------------------------------
C     CALLS THE APPROPRIATE ELEMENT MODULE
C     ------------------------------------
      
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


C	INTRODUCE BY DE SILVA
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT

C	COMMON FOR EAS ADDED BY SONGSAK FEB2006
	COMMON /MMENH/ MM,MM1,MM2,NDIMC

C	COMMON BLOCK FOR HEAT SONGSAK MAR2007
	COMMON /SHEAT/ NHAEF,LHEAT1,LHEAT2,LHEAT3,LHEAT4


C
      DIMENSION PROPM(1),PROPG(1),NODEX(1),WA(1),S(1),COORD(1),EDIS(1) 
C	ANG,AMV - ADDED TO NEXT LINE BY GILSON - SEPT2002
C	PROPO,FEF - ADDED TO NEXT LINE BY DE SILVA - OCT2003
      DIMENSION EDISI(1),ELOD(1),AMV(3),PROPO(1)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(1)


C	CHANGED BY SONGSAK FEB2006
	DIMENSION ALPHA(MM,1),SEL(MM,MM1),SEDI(MM,MM),RH(MM,1),
	1		  HINFC(MM)

	DIMENSION SEQ(MM,MM2),HINFP(MM) !CONSO

      DIMENSION WA2(1) !---------NEW WORKING ARRAY BY BJ

C	HINFP : EQUIVALENT EAS LOADING DUE TO PRESSURE
C	SEQ   : ENHANCED ASSUMED STRAIN COUPLING FLUID MATRIX
	

      DIMENSION IPIN(14)
	DIMENSION FIXEN(NEF),FIXLR(NEF),FIXEO(NEF),FIXLO(NEF),CABFF(1)
	DIMENSION TCHET(NHAEF,NELE),TAMBT(1)        !HEAT CONDUCTION FLAC FOR ELEMENT (NFACE)
	DIMENSION PMATRL(1)

      DIMENSION ISFAC(1)
      
      DIMENSION DISLI(1)
C	-------------------------------------------------------------------

      GOTO (110,120,130,140,150,160,170,180,190,200,
	1	  210,220,230,240,250,260,270,280,290,300),ITYPE

C	 2D TRUSS
C	--------------------------------------------------
110	CONTINUE
	RETURN 

C	 3D TRUSS
C	--------------------------------------------------------------
120	CALL TRUS3DM (PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,FIN
	1			 ,MSET,CABFF,DISLI)

	RETURN

 130  RETURN
 140  RETURN


c	--------------------------------------------------------------
 150  CALL BEAM21 (PROPM,PROPG,PROPO,WA,S,COORD,EDIS,
	1			EDISI,ELOD,NWG,FIN,IPIN,NCF,ISET,
     2			FIXEN,FIXLR,FIXEO,FIXLO,PMATRL,
     3            FIXEN_OFF,FIXLR_OFF,FIXEO_OFF,FIXLO_OFF)			
     	
	RETURN


c	--------------------------------------------------------------
 160  CONTINUE
	IF (MTYP.EQ.0) THEN !STANDARD MEMBRANE AND EAS

	WRITE (*,*) "MEMBRANE NOT AVAILABLE IN THIS VERSION"
      STOP
	ELSE !MEMBRANE WITH DRILLING DEGREE OF FREEDOM

      WRITE (*,*) "MEMBRANE NOT AVAILABLE IN THIS VERSION"
      STOP
	ENDIF

      RETURN


C	PLATE ELEMENT
C	--------------------------------------------------------------
 170  CONTINUE
	CONTINUE
      RETURN

C	2D SEEPAGE ANALYSIS COUPLING ELEMENT
C	--------------------------------------------------------------
 180  CONTINUE
	WRITE (*,*) "2D SEEPAGE NOT AVAILABLE IN THIS VERSION"
      STOP
	RETURN

C	SHELL ELEMENT
c	--------------------------------------------
C	ANG,AMV - ADDED TO NEXT LINE BY GILSON - SEPT2002
C	FIN - ADDED TO NEXT LINE BY GILSON - JUL2003 (INT FORCE)
 190  CONTINUE
      WRITE (*,*) "SHELL NOT AVAILABLE IN THIS VERSION"
      STOP
      RETURN

C	SOLID ELEMENT
C	--------------------------------------------------
C	FIN - ADDED TO NEXT LINE BY GILSON - JUL2003 (INT FORCE)
 200  CONTINUE
      WRITE (*,*) "SOLID ELEMENT NOT AVAILABLE IN THIS VERSION"
      STOP
C	--------------------------------------------------
	RETURN

C	==================================================
C	ADDITIONAL ELEMENT LIBRARY
C	ADDED BY SACHARUCK MAR2007 (SONGSAK IMPLEMENTOR)
C	==================================================
C	--------------------------------------------------
210   CONTINUE
C	SOLID CONSOLIDATION 
      WRITE (*,*) "SOLID CONSOLIDATION ELEMENT NOT AVAILABLE IN THIS VERSION"
      STOP
	RETURN
C	--------------------------------------------------
220	CONTINUE
	RETURN
C	--------------------------------------------------
230   CONTINUE
      WRITE (*,*) "NOT AVAILABLE IN THIS VERSION"
	RETURN
C	--------------------------------------------------
240   CONTINUE
      WRITE (*,*) "NOT AVAILABLE IN THIS VERSION"
		
	RETURN
C	--------------------------------------------------
250	CONTINUE
      WRITE (*,*) "NOT AVAILABLE IN THIS VERSION"
	RETURN
C	--------------------------------------------------
260   CONTINUE
      WRITE (*,*) "NOT AVAILABLE IN THIS VERSION"
	RETURN
C	--------------------------------------------------
270   CONTINUE  !TENDON
      WRITE (*,*) "NOT AVAILABLE IN THIS VERSION"
	RETURN
C	--------------------------------------------------
280	CONTINUE
	RETURN
C	--------------------------------------------------
290	CONTINUE
	RETURN
C	--------------------------------------------------
300   CONTINUE
      WRITE (*,*) "NOT AVAILABLE IN THIS VERSION"
	RETURN
C	--------------------------------------------------
	
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ELCORD (XYZ,LM,DISPI,DISP,COORD,EDIS,EDISI,MEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     DETERMINES DISPLACEMENTS, DISP.INCREMENTS AND COORDINATES AT
C     FINITE ELEMENT NODES
C	--------------------
C     XYZ(NXY,NELE)    = COORDINATES FOR ELEMENT GROUP
C     LM(NEF,NELE)     = GLOBAL EQUATION NUMBERS FOR ELEMENT D.O.F.
C     DISPI(NEQ)       = ACCUMULATED DISPL.INCREMENTS FOR CURRENT STEP
C     DISP(NEQ)        = TOTAL DISPLACEMENT VECTOR
C     COORD(NCO,NNM)   = (CURRENT) COORDINATES FOR ELEMENT NODES
C     EDIS(NEF)        = TOTAL DISPLACEMENTS AT ELEMENT NODES
C     EDISI(NEF)       = INCREMENTAL DISPLACEMENTS AT ELEMENT NODES
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C      COMMON /ITER/ NSTEP,DLMAX,NPRIN,NDRAW,KONEQ,NIREF,ITEMAX,RTOL,
C     1              ITOPT,NUMREF,NUMITE,ITETOT,RHO,RHOP,ICONV,NOLIN,
C     2              LIMEQ(2),KSTEP,ETOL,RHOPREV
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC, !temp, can be suppressed
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON A(9000000),IA(9000000)
C
c      DIMENSION XYZ(1),LM(MEF,1),DISPI(1),DISP(1)
      DIMENSION XYZ(1),LM(MEF,1),DISPI(neq),DISP(neq)
C      DIMENSION COORD(1),EDIS(1),EDISI(1)
      DIMENSION COORD(12),EDIS(NEF),EDISI(NEF)
C     --------------
C     INITIALISATION
C     --------------
      NXY = NCO*NNM
      CALL CLEARA (COORD,NXY)
      CALL CLEARA (EDIS,NEF)
      CALL CLEARA (EDISI,NEF)

C     -------------------
C     NODAL DISPLACEMENTS
C     -------------------
      IF (ITASK.EQ.5) GOTO 300

	IF (NOLIN.EQ.0 .AND. ITASK.EQ.1) GOTO 300

	DO 200  IDF=1,NEF
      IEQ = LM(IDF,MEL)
      IF (IEQ.EQ.0) GOTO 200
      EDISI(IDF) = DISPI(IEQ)
      EDIS(IDF)  = DISP(IEQ)
      IF (ITASK.NE.2) GOTO 200
	EDIS(IDF)  = EDIS(IDF) + DISPI(IEQ)
 200  CONTINUE

	
      IF (NLS.LE.0) GOTO 300
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,MEL),A(LES),EDIS,
     1             EDISI,A(LEE),NSF,NNF,1)
C     ------------------------------------------------
C     NODAL COORDINATES (TOTAL LAGRANGIAN FORMULATION)
C     ------------------------------------------------
 300  K = (MEL-1)*NXY
      COORD(1:NXY) = XYZ(K+1:K+NXY)
CC      DO 390  I=1,NXY
CC 390  COORD(I) = XYZ(K+I)
C     ---------------------------------------------------
C     UPDATE COORDINATES (UPDATED LAGRANGIAN FORMULATION)
C     ---------------------------------------------------
      IF (NLOPT.LT.3) RETURN
      DO 400  I=1,3
      IPO = IDOF(I)
      IF (IPO.LE.0 .OR. IPO.GT.3) RETURN
      KDI = I-NNF
      DO 400  KCO=IPO,NXY,NCO
      KDI = KDI+NNF
 400  COORD(KCO) = COORD(KCO) + EDIS(KDI)

C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE EQUITE(R,RO,RFAC,RE,DISPI,DISLI,DISP,ID,KSTRA,KITE,AA,
	1				  LTYP,ILC)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 LTYP
C     ----------------------------------------------------------------
C     ITERATES FOR EQUILIBRIUM
C	------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     R(NEQ)      = REFERENCE LOAD VECTOR
C     RFAC(NSTEP) = LOAD FACTORS TO DEFINE LOADING LEVEL
C     RE(NEQ)     = EQUILIBRIUM LOADS BALANCED IN CURRENT DISP.FIELD
C     DISPI(NEQ)  = SUM OF DISPL.INCREMENTS AT CURRENT LOAD STEP
C     DISLI(NEQ)  = 1.  NON-CONSERVATIVE AND TOTAL VECTOR
C                   2.  RESIDUAL LOAD VECTOR
C                   3.  ADDITIONAL LAST INCREMENT IN DISPLACEMENT
C     DISP(NEQ)   = TOTAL DISPLACEMENT VECTOR
C     KSTEP       = STEP NUMBER
C	-------------------------------------------
C     VARIABLES IN COMMON BLOCK /ITER/ AND /FLAG/
C	-------------------------------------------
C     NIREF       = NUMBER OF ITERATIONS BETWEEN REFORMING STIFFNESS
C     ITEMAX      = MAXIMUM NUMBER OF ITERATIONS PERMITTED
C     RTOL        = RELATIVE TOLERANCE TO MEASURE CONVERGENCE
C     ITOPT       = ITERATION OPTION (1=LOAD CONTROLED,LOAD FIXED,
C                                     2=LOAD INITIATED,LOAD FREE,
C                                     3=DISPLACEMENT CONTROLED)
C     NUMREF      = NUMBER OF STIFFNESS REFORMATIONS
C     NUMITE      = NUMBER OF PERFORMED EQUILIBRIUM ITERATIONS
C     RHO         = LOAD FACTOR WHICH DEFINES LOADING LEVEL WHICH
C                   IS BEST ADAPTED TO CURRENT DISP.FIELD (DISP)
C     ICONV       = FLAG FOR CONVERGENCE (1) OR DIVERGENCE (-1)
C     NOLIN       = FLAG FOR TYPE OF NONLINEAR ANALYSIS
C                   1 = MATERIALLY NONLINEAR ONLY
C                   2 = GEOMETRICALLY NONLINEAR ONLY
C                   3 = MATERIALLY AND GEOMETRICALLY NONLINEAR
C     ISPRI       = FLAG TO SUPPRESS STRESS OUTPUT (ISPRI=1)
C     IFREF       = FLAG FOR STIFFNESS REFORMATION (IFREF=0)
C	---------------
C     LOCAL VARIABLES
C	---------------
C     TOL         = RTOL*DNORM MAXIMUM ADMISSIBLE NORM OF DISPI
C     RRNORM      = NORM OF REFERENCE LOAD VECTOR
C     RENORM      = NORM OF RESIDUAL- TIMES REFERENCE LOAD VECTOR
C     DNORM       = NORM OF TOTAL DISPLACEMENT VECTOR DISP(ITE)
C     DINORM      = NORM OF LAST INCREMENT IN DISPLACEMENT (DISLI)
C     DISNRM      = NORM OF INITIAL FIRST DISPLACEMENT INCREMENT
C     -----------------------------------------------------------------
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

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
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON A(9000000),IA(9000000)
	

      DIMENSION R(NEQ),RO(NEQ),RFAC(NSTEP),RE(NEQ),DISPI(NEQ),
	1          DISLI(NEQ),DISP(NEQ)

	DIMENSION ID(NSF,NSN)
	DIMENSION AA(1)

	IF(LTYP.EQ.'CONT') ILPT = 0  !FOR CONSTANT LOAD ANALYSIS WE DONT NEED TO FIND THE LOAD FACTOR
	IF(LTYP.EQ.'VARY') ILPT = 1  !FOR VARY LOAD
C     --------------
C     INITIALISATION
C     --------------
      CALL CPU_TIME (TIM1)
      KITE  = 0
      ITASK = 2
      ISPRI = 1
      INDPD = KPOSD
      IFREF = 1
      IF (KSTRA.GE.2)  IFREF = 0
      
	CALL VENORM (DISPI,DISPI,DISPI,RHO,DISNRM,NEQ)

      GOTO 110
C     -----------------------
C     START OF ITERATION LOOP
C     -----------------------
	
 100  IFREF = 1

      IF (KSTRA.GE.3)  IFREF = 0
 110  KITE = KITE+1

      IF (IFPR(9).GE.1)  WRITE (ISO,1000)  KSTEP,KITE
C     -------------------------------------------------------
C     FIND EQUIL.LOADS BALANCED IN CURRENT DISPL.FIELD (DISP)
C     FIND NEW TANGENTIAL STIFFNESS MATRIX (IFREF=0)
C     -------------------------------------------------------
      CALL CLEARA (RE,NEQ)

	IF (ITYPE .NE. 5 .OR. ITYPE .NE. 9 )GOTO 120
	IF (NLOPT .EQ. 1) GOTO 120
	CALL ADDROT (IA(LID),DISP,DISPI,DISLI,NSF,ISO)
 120  CALL GRLOOP (IA(LEL),KSC)

      IF (IFPR(6)+IFPR(9).EQ.2) CALL MATOUT (RE,NEQ,1,1,22,'E',
	1									   15,10,2,'EQUIL LOAD')
C     ----------------------------------------------------------
C     FIND CURRENT REFERENCE LOAD VECTOR AND LOADING LEVEL RHO
C     WHICH IS BEST ADAPTED TO CURRENT DISPLACEMENT FIELD (DISP)
C     ----------------------------------------------------------
	CALL CLEARA (DISLI,NEQ)
	
C     ----------------------------------------------------------
C              READ OFFSHORE FORCE WITH LINK ELEMENT 
C     ----------------------------------------------------------	
	IF (KOFFL.EQ.1.OR.KSPEC.EQ.1)THEN
      CALL OFFSHFORC(R,ILC,1,NEQ,'VARY','RADD')
      ENDIF

      CALL VECADD (R,DISLI,DISLI,RHO,RRNORM,NEQ)		  !DRIVEN LOAD

	IF (ITYPE .NE. 5) GOTO 200
	IF (NLOPT .EQ. 1) GOTO 200
C	CALL MOTRAN (IA(LID),DISP,DISPI,DISLI,NSF)
C     ------------------------------------------------------------
C     CALCULATE RESIDUAL LOAD VECTOR DISLI(I) = RHO*DISLI(I)-RE(I)
C     ------------------------------------------------------------
 200  CALL VECMUL (DISLI,DISLI,DISLI,RHO,RRNORM,NEQ)

	CALL VECADD (DISLI,RE,DISLI,RHO,RRNORM,NEQ)       !UNBALANCE LOAD = DRIVEN LOAD - RESISTING LOAD

	CALL VECADD (DISLI,RO,DISLI,RHO,RRNORM,NEQ)		  !ADD CONTSTANT LOAD VECTOR TO UNBALANCE LOAD

      IF (IFPR(6)+IFPR(9).EQ.2) CALL MATOUT (DISLI,NEQ,1,1,22,'E',
	1									   15,10,2,'RESID LOAD')
      CALL CPU_TIME (TIM2)
      TIM(16) = TIM(16) + (TIM2-TIM1)
C     ------------------------------------------------------
C     ASSEMBLE GLOBAL COMPACTED STIFFNESS BLOCKS (IFREF=0)
C     TRIANGULARIZE EFFECTIVE STIFFNESS MATRIX (IFREF=0)
C     CALCULATE ADDITIONAL INCREMENT IN DISPLACEMENT (DISLI)
C     ------------------------------------------------------
      IF (IFREF.NE.0)  GOTO 310
      NUMREF = NUMREF+1

	CALL COLSOL (IA(LMA),AA,A(LDK),DISLI,1,INDPD,'STIF','TEMP')
310	CALL COLSOL (IA(LMA),AA,A(LDK),DISLI,2,INDPD,'TEMP','TEMP')

      IF (ITOPT.EQ.1)  GOTO 320
	IF(ILPT.EQ.0) GOTO 320												!FOR CONSTANT LOAD ANALYSIS WE DONT NEED TO FIND THE LOAD FACTOR
      CALL CLEARA (RE,NEQ)
      CALL VECADD (R ,RE,RE,RHO,DNORM,NEQ)								!REFERENCE LOAD
C	CALL VECADD (RO,RE,RE,RHO,DNORM,NEQ)								!ADD CONSTANT LAOD TO REFERENCE LOAD
      CALL COLSOL (IA(LMA),AA,A(LDK),RE,2,INDPD,'TEMP','TEMP')
      CALL LOADUP (IA(LID),RFAC,RE,DISLI,A(LDK),2,NSF)
 320  CALL CPU_TIME (TIM1)
	IF (IFPR(7)+IFPR(9).EQ.2) CALL MATOUT (DISLI,NEQ,1,1,22,'E',
	1									   15,10,2,'DISPL INCR')
C     -----------------------------------------------
C     UPDATE DISP.INCR.(DISPI) AND TOTAL DISPL.(DISP)
C     CHECK FOR CONVERGENCE OR DIVERGENCE
C     -----------------------------------------------
	CALL VECADD (DISLI,DISPI,DISPI,RHO,DNORM,NEQ)
	CALL ADDNRM (DISP ,DISPI,DISPI,RHO,DNORM,NEQ)
	CALL VENORM (DISLI,DISLI,DISLI,RHO,DINORM,NEQ)

      IF (IFPR(9).LT.1)  GOTO 350
      CALL VENORM (DISPI,DISPI,DISPI,RHO,DINRM,NEQ)
      DTNRM = DSQRT(DNORM)
      DANRM = DSQRT(DISNRM)
      DLNRM = DSQRT(DINORM)
      DINRM = DSQRT(DINRM)
      WRITE (ISO,2000) RHO,DTNRM,DANRM,DLNRM,DINRM
C
 350  TOL = DNORM*RTOL

      IF (DINORM.LT.TOL)			GOTO 500

	IF (DINORM.GT.500.*DISNRM)  GOTO 360  !Changed from 100.*DISNRM to 500.*DISNRM by Songsak Nov2006

      IF (KITE.LT.4)				GOTO 400
      IF (DINORM.LT.DISNRM)		GOTO 400

 360  ICONV = -2
      GOTO 900

 400  IF (KITE.LT.ITEMAX)			GOTO 100
      ICONV = -1
      GOTO 900
C
 500  ICONV = 1
C
 900  NUMITE = NUMITE+KITE
      ITETOT = ITETOT+KITE
      CALL CPU_TIME (TIM2)
      TIM(16) = TIM(16) + (TIM2-TIM1)
C
 1000 FORMAT (/////,22X,31HL O A D   S T E P   N U M B E R,I4//
     1        22X,31HI T E R A T I O N   N U M B E R,I4//22X,35(1H*)//)
 2000 FORMAT (//25X,30HNORMS TO CHECK FOR CONVERGENCE/25X,30(1H-)//
     1        5X,40HFACTOR TO DEFINE LOADING LEVEL    RHO = ,E20.12/
     2         5X,40HTOTAL DISPLACEMENT VECTOR . . . DTNRM = ,E20.12/
     3          5X,40HAPPLIED DISPLACEMENT VECTOR . . DANRM = ,E20.12/
     4           5X,40HLAST DISPLACEMENT INCREMENT . . DLNRM = ,E20.12/
     5            5X,40HSUM OF DISP.INC. FOR STEP . . . DINRM = ,E20.12)
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE RETAKE (ID,RFAC,DISPI,DISLI,DISP,DD,
     +                   KRECO,KSTRA,KITE,MSF,AA,ITIME,NFATIGUELOOPOUT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CONTROLS CHANGE OF NONLINEAR SOLUTION METHOD
C     INITIALISES CONTROL FLAGS AND ARRAYS TO RETAKE A STEP IN CASE
C     CONVERGENCE FAILURE OCCURED
C     ---------------------------
C     INPUT, OUTPUT VARIABLES
C	-----------------------
C     RFAC(NSTEP) = LOAD OR DISPLACEMENT INCREMENTS
C     KRECO       = RECOVERY COUNTER
C     KSTRA       = RETURNS STRATEGY FOR NEXT STEP (ICONV > 0) OR
C                   FOR CURRENT STEP (ICONV < 0)
C     KITE        = NUMBER OF ITERATIONS USED TO ESTABLISH EQUILIBRIUM
C	IFFLAG      = FLAG FOR COMPOSITE PLY FAILURE CHECKING
C			  2 = FOR NO PLY DAMAGE CHECKING
C			  1 =    PLY FAILURE OCCURRED
C			  0 = NO PLY FAILURE OCCURRED
C     ----------------------------------------------------------------
      LOGICAL PROMPT,ERROR
C
      CHARACTER*2 MARK
	CHARACTER*6 IFREE(9)
	CHARACTER*6 IPO
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /LOGO/ PROMPT,ERROR,ITEST
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      
      COMMON /NSCON/ NOUTS
      COMMON /RESO/ OPRES,STEPSTAT,STEPINCR,STEPEND

C	THIS COMMON LINE ADDED BY DE SILVA FOR SEEPAGE 30/04/2004
	COMMON /SEEP/ NTSTEP,KTSTEP,CTIME,DTINC,KFRES
      
      COMMON /STCAS/ ILC

      COMMON A(9000000),IA(9000000)
C
      DIMENSION ID(MSF,1),RFAC(1),DISPI(1),DISLI(1),DISP(NEQ)
      DIMENSION DD(1),AA(1)
      DATA IFREE /'X-DISP','Y-DISP','Z-DISP','X-ROTA','Y-ROTA',
     +            'Z-ROTA','WARP. ','TEMP. ','PRES. '/
C     --------------------------------------------
C     INITIALISATION, PRINT TITLE FOR JOB PROGRESS
C     --------------------------------------------
      LSTRA = KSTRA
      KREF  = KSTRA
      IF (KSTRA.EQ.0 .AND. KSTEP.EQ.1) KREF = 1
      IF (KSTRA.EQ.3) KREF = KITE +1
      IF (KSTRA.EQ.4) KREF = 1
C
      IF (KSTEP.GT.1 .OR. KRECO.NE.0) GOTO 10
      
      !IF (ICONTROLSPEC.EQ.0)THEN
      WRITE (100,1000)
      WRITE (ITO,1000)
      WRITE (10,1000)
      !ELSEIF (ICONTROLSPEC.EQ.1)THEN
      !    IF (NOUTS.EQ.0)THEN
      !    WRITE (100,1001)
      !    WRITE (100,1002)
      !    WRITE (ITO,1001)
      !    WRITE (ITO,1002)
      !    NOUTS = 1
      !    ELSEIF (NOUTS.EQ.1)THEN
      !    ENDIF
      !ENDIF

 50   IF (ICONV.LT.0) KRECO = 4

 10	IF (IFFLAG.EQ.1) GOTO 200  !IF PLY FAILURE OCCURRED COMPOSITE

	IF (ICONV) 200,100,100     ! 200= CONVERGENCE FAILURE 100=CONVERGENCE ACHIEVED
C     -----------------------------------------------------
C     CONVERGENCE ACHIEVED, SELECT LOWER STRATEGY (KRECO=0)
C     ADJUST SIZE OF NEXT INCREMENT ACCORDING TO KITE
C     -----------------------------------------------------
 100  RHOP = RHO
      MARK = '  '
      IF (KRECO.EQ.4) MARK = '**'
      IF (KRECO.EQ.0 .OR. KITE.GE.ITEMAX/2) GOTO 110
      KRECO = KRECO-1
      IF (KRECO.GT.0 .OR. KSTRA.EQ.NIREF  ) GOTO 110
      KSTRA = KSTRA-1
      KRECO = 4
C
 110  RSCALE = 1.0
      IF (RFAC(KSTEP+1).NE.0.0.OR.FLOAT(KSTEP).EQ.FLOAT(NSTEP))GOTO 400
      IF (ICONV.EQ.0)  GOTO 120
      RSCALE = SQRT(0.4*FLOAT(ITEMAX)/FLOAT(KITE))
 120  DRHO = RFAC(KSTEP)*RSCALE
      DRHOMX = RFAC(1)*DLMAX
      IF (DABS(DRHO).GT.DABS(DRHOMX)) DRHO = DRHOMX
      RFAC(KSTEP+1) = DRHO
      GOTO 400

C     ------------------------------------------------------
C     INITIALISE TOTAL LOAD LEVEL RHO
C     REDUCE LOADING (DISPLACEMENT) INCREMENT AND CLEAR LOAD
C     FACTORS > KSTEP
C     ------------------------------------------------------
 200  FAC = 4.0
      IF (KRECO.EQ.4 .AND. NUMITE.EQ.KITE) FAC = 1.0

	IF (IFFLAG.EQ.1) FAC     = 1.0						!FOR COMPOSITE PLY DAMAGE---NO MODIFICATION OF LOAD STEP SIZE
	IF (IFFLAG.NE.1) KSTRA   = KSTRA  + 1				!IF NO COMPOSITE CHECK (STANDARD)
	IF (IFFLAG.NE.1) RTOL    = RTOL   + 10.0*RTOL		!IF NO COMPOSITE CHECK (STANDARD) ADD HERE BY SONGSAK TO INCREASE THE TOLERENCE IF CONVERGENCE FAILURE OCCURRED
	IF (RTOL.GT.0.1) RTOL = 0.1
	IF (IFFLAG.NE.1) ITEMAX  = ITEMAX + 2*ITEMAX		!IF NO COMPOSITE CHECK (STANDARD) ADD HERE BY SONGSAK TO INCREASE THE NUMBER OF ITERATION IF CONVERGENCE FAILURE OCCURRED

      KRECO = 4
      IFEIG = 1
      ISPRI = 1
      RHO   = RHOP
      MARK  = '  '
      RFAC(KSTEP) = RFAC(KSTEP)/FAC
      IF (NIREF.EQ.0) NIREF = 1
      IF (KSTEP.EQ.NSTEP) GOTO 400
      KSTEP1 = KSTEP+1
      DO 350  ISTEP=KSTEP1,NSTEP
 350  RFAC(ISTEP) = 0.0
C     --------------------------------------------------------
C     PRINT JOB SUMMERY ONTO TAPE DATIN (TAPE 1)
C     --------------------------------------------------------
 400  NODE = LIMEQ(1)
      IPOS = LIMEQ(2)
      IF (ITOPT.NE.3) GOTO 450
C
      DO 410  ISN=1,NSN
      DO 410  ISF=1,NSF
      IF (KONEQ.NE.ID(ISF,ISN)) GOTO 410
      IPOS = IDOF(ISF)
	IF (IPOS.LE.0.OR.IPOS.GT.3) GOTO 410
      NODE = ISN
      GOTO 450
 410  CONTINUE
C
 450  DELTA = 0.
      IPO   = IFREE(IPOS)
      IEQ   = ID(LIMEQ(2),LIMEQ(1))
      DELTA = DISP(IEQ)
      !IF(DELTA.LT.1.0E-5) DELTA(3) = 0.0D0
      IF (NOLIN.NE.0.AND.ISOLOP.NE.1) DELTA = DELTA + DISPI(IEQ)

	IF (IFFLAG.EQ.1) GOTO 505    !FOR COMPOSITE PLY DAMAGE

 500  IF (ICONV+1) 510,520,530

C	FOR COMPOSITE
 505	WRITE (100,1250) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,DETK
	WRITE (ITO,1250) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,DETK
      WRITE (10,1250) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,DETK
	IF (IFFLAG.EQ.1.AND.KSTEP.GE.NSTEP) GOTO 530
	GOTO 590

C	DIVERGENCE OCCURRED  ICONV = -2
 510  WRITE (ITO,1100) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,DETK
      WRITE (10,1100) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,DETK
      GOTO 590

C	CONVERGENCE FAILURE  ICONV = -1
 520  WRITE (ITO,1200) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,DETK
      WRITE (10,1200) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,DETK
      GOTO 590

C	CONVERGENCE ARCHIEVE ICONV >= 0
      !IF (ICONTROLSPEC.EQ.0)THEN
530   IF(LSTRA.NE.4) THEN ! SONGSAK
	WRITE (100,1300) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,
     1                 DETK,DELTA,RHO
	ENDIF
	WRITE (ITO,1300) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,
     1                 DETK,DELTA,RHO
      WRITE (10,1300) MARK,KSTEP,LSTRA,KREF,KITE,NODE,IPO,
     1                 DETK,DELTA,RHO
      !ELSEIF (ICONTROLSPEC.EQ.1)THEN
      !IF (OPRES.EQ.1) CALL OFFSHSTEP(TIME,ITIME,NTIME,'CALL:') !READ TIME DATA FOR OFFSHORE LOAD 
      !IF (OPRES.EQ.2) CALL SELECTGRAPH (2,0,ITIME,TIME)
      IF(LSTRA.NE.4) THEN ! SONGSAK
	!IF (ICONTROLSPEC.NE.1) 
      WRITE (100,1301) ITIME,ILC,NODE,TIME,DELTA    
	ENDIF
  	!WRITE (ITO,1301) ITIME,ILC,NODE,TIME,DELTA     
      !ENDIF 

 590  RETURN

C
 1000 FORMAT (//32X,16(1H*)/32X,1H*,14X,1H*/
     1        32X,16H* JOB PROGRESS */32X,1H*,14X,1H*/32X,16(1H*)//
     2        42H KSTEP  KSTRA   KREF   KITE   NODE    IPOS
     3        ,38H      DETK        DELTA           RHO /1X,79(1H-)/)
1001  FORMAT (//32X,16(1H*)/32X,1H*,14X,1H*/
     1        32X,16H* JOB PROGRESS */32X,1H*,14X,1H*/32X,16(1H*)//)
1002  FORMAT  (" STEP     LOADCASE   NODE    FREQ.    DISPLACEMENT ")
 1100 FORMAT (1X,A1,I3,4I7,4X,A6,E11.3,2X,24H * DIVERGENCE OCCURRED *)
 1200 FORMAT (1X,A1,I3,4I7,4X,A6,E11.3,2X,24H * CONVERGENCE FAILURE *)
C	NEXT FORMAT LINE ADDED BY GILSON - JULY2002
 1250 FORMAT (1X,A1,I3,4I7,4X,A6,E11.3,2X,24H * PLY FAILURE OCCURRED*)
 1300 FORMAT (1X,A1,I4,4I7,3X,A6,E11.3,2E13.4)
 1301 FORMAT (1X,I4,5X,I4,6X,I4,1X,E11.3,2X,E11.3)
c	next 2 lines added by gilson - sept2004 (graph)
1400	FORMAT ('#X:"DISPLACEMENT" Y:"LOAD FACTOR"')
1500	FORMAT (E15.6,2x,E15.6)
C
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE LOADUP (ID,RFAC,RE,DISLI,DD,IENTRY,MSF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     SELECTS INITIAL LOAD INCREMENT (IENTRY<2) OR LOAD INCREMENT
C     FOR NEXT ITERATIVE CYCLE (IENTRY=3)
C	-----------------------------------
C     RFAC(NSTEP)  = PRESCRIBED INCR. LOAD FACTORS / DISP. INCREMENTS
C                    (RFAC(KSTEP)=0 INDICATES AUTOMATIC SELECTION)
C     RE(NEQ)      = DISPLACEMENT FIELD CORRESPONDING TO A REFERENCE
C                    LOAD OF ARBITRARY MAGNITUDE (DIRECTION OF MOTION)
C     DISLI(NEQ)   = EXTERNAL APPLIED REFERENCE LOAD VECTOR (IENTRY=0)
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
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
      DIMENSION ID(MSF,NSN),RFAC(1),RE(1),DISLI(NEQ)
	DIMENSION DD(1)
C
      DU = 0.
      IF (IENTRY-1) 100,200,500
C     ------------------------------------
C     SET INITIAL LOAD INCREMENT (ITOPT=1)
C     ------------------------------------
100	IF (ITOPT.GT.1) GOTO 150  !FOR DISPLACEMENT CONTROL AND AUTOMATIC CONTROL

      DRHO = RFAC(KSTEP)
      RHO  = RHO + DRHO

150	CALL VECMUL (DISLI,DISLI,DISLI,RHO,RNORM,NEQ)  

      RETURN
C     -----------------------------------------------------------
C     CHECK VALIDITY OF PRESCIBED DOF AND INITIALIZE EQUATION NO.
C     -----------------------------------------------------------
 200  IF (KSTEP.GT.1 .OR. ICONV.LT.0) GOTO 230
      NODE = LIMEQ(1)
      IPOS = LIMEQ(2)
      IF (NODE.LT.1 .OR. NODE.GT.NSN) GOTO 210
      IF (IPOS.LT.1 .OR. IPOS.GT.NSF) GOTO 210
      KONEQ = ID(IPOS,NODE)
      IF (KONEQ.NE.0) GOTO 230
C
 210  IF (ITOPT.EQ.2) ITOPT = 3
      CALL VECMAX (RE,NEQ,REMAX,KONEQ)
      DO 220  ISN=1,NSN
      DO 220  ISF=1,NSF
      IF (KONEQ.NE.ID(ISF,ISN)) GOTO 220
      LIMEQ(1) = ISN
      LIMEQ(2) = ISF
      GOTO 230
 220  CONTINUE

C     -------------------------------------------------------
C     CHOOSE APPROPRIATE SIGN FOR INITIAL TANGENTIAL SOLUTION
C     -------------------------------------------------------
 230  IF (ITOPT.LE.2) GOTO 300
      S = 1.0D0	
      
C     PIVOT FLAG  0=POSITIVE DEF  1=NEGATIVE DEF  ONLY APPEARED IN PARDISO,NEWSOLVER,LOADUP
      CALL MINTFIL('BLOK',IPIVT ,1,6 ,0)  
      IF(IPIVT.EQ.1) GOTO 280
      
      DO 270  IEQ=1,NEQ
      IF (DD(IEQ)) 280,270,270
 270  CONTINUE
      GOTO 290
 280  S = -1.0D0
C     ---------------------------------
C     SELECT APPROPRIATE NODE (ITOPT=3)
C     ---------------------------------
 290  IF (ITOPT.NE.3) GOTO 300
	REMAX = 0.
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
C     SET FINAL DISPLACEMENT INCREMENTS DISLI = DISLI + DRHO*RE
C     ---------------------------------------------------------
 500  DISUNB = DISLI(KONEQ)
      DISTAN = RE(KONEQ)
	DRHO   = (DU-DISUNB)/DISTAN
      
 590  RHO = RHO + DRHO
      CALL VECMUL (RE,RE,RE,DRHO,DNORM,NEQ)
      CALL VECADD (RE,DISLI,DISLI,DRHO,DNORM,NEQ)
C
      IF (IFPR(9).EQ.3) WRITE (ISO,1000) KONEQ,DU,DISUNB,DISTAN,DRHO,RHO

C
 1000 FORMAT (//21X,37HAUTOMATIC SELECTION OF LOAD INCREMENT/
     1        21X,37(1H-)//
     2        5X,'PRESCRIBED OR SELECTED EQUATION    KEQ =',I20/
     3        5X,'PRESCRIBED DISPLACEMENT PARAMETER   DU =',E20.12/
     4        5X,'LAST DISP. INCREMENT / NORM      DISLI =',E20.12/
     5        5X,'TANGENTIAL DISP. FIELD / NORM    RE    =',E20.12/
     6        5X,'INCREMENTAL CHANGE OF LOAD LEVEL DRHO  =',E20.12/ 
     7        5X,'CURRENT TOTAL LOADING LEVEL      RHO   =',E20.12)
C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE LDFACR (ID,RFAC,RE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
C
      DIMENSION ID(NSF,1),RFAC(1),RE(1)

	IF (ITOPT.NE.3) RETURN
C     -----------------------------------------------------------
C     CHECK VALIDITY OF PRESCIBED DOF AND INITIALIZE EQUATION NO.
C     -----------------------------------------------------------
 200  NODE = LIMEQ(1)
      IPOS = LIMEQ(2)
      IF (NODE.LT.1 .OR. NODE.GT.NSN) GOTO 210
      IF (IPOS.LT.1 .OR. IPOS.GT.NSF) GOTO 210
      KONEQ = ID(IPOS,NODE)
      IF (KONEQ.NE.0) GOTO 230
C
 210  CALL VECMAX (RE,NEQ,REMAX,KONEQ)
      DO 220  ISN=1,NSN
      DO 220  ISF=1,NSF
      IF (KONEQ.NE.ID(ISF,ISN)) GOTO 220
      LIMEQ(1) = ISN
      LIMEQ(2) = ISF
      GOTO 230
 220  CONTINUE

C     ---------------------------------
C     SELECT APPROPRIATE NODE (ITOPT=3)
C     ---------------------------------
 230  IF (ITOPT.NE.3) GOTO 300
	REMAX = 0.
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
 300  RFAC(1) = RFAC(1) + ABS(RE(KONEQ))
	

      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE NEWDIS (ID,DISPI,DISP,MSF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------
C     PRINTS NEW DISPLACEMENTS FOR STEP KSTEP
C     ---------------------------------------
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON /PLOT/ LODE(2,10)


      DIMENSION ID(MSF,1),DISPI(1),DISP(1)

     
C     -----------------------------------------------------
C     UPDATE TOTAL DISPLACEMENTS DISP(I) = DISP(I)+DISPI(I)
C     -----------------------------------------------------
	CALL VECADD (DISP,DISPI,DISP,1.0,DNORM,NEQ)

C     ------------------------------
C     PRINT HEADER AND DISPLACEMENTS
C     ------------------------------
 100  IF (IFPRI.GT.0)  GOTO 200

C	NEW OUTPUT SONGSAK JUL2007
	CALL DISOUT(ID,DISP)

C     --------------------------------------------------------
C     WRITE DATA REQUIRED FOR GRAPHICAL OUTPUT ONTO TAPE NPLOT
C     --------------------------------------------------------
 200  RETURN

C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE COLSOL_ORIGINAL (MAXA,A,D,V,IND,INDPD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SUNIL 06/02/01 COMPLETELY MODIFIED
C     ----------------------------------------------------------------
C     SOLVE FINITE ELEMENT STATIC EQUILIBRIUM EQUATIONS USING IN CORE
C	SOLVER BASED ON COMPACTED STORAGE AND COLUMN REDUCTION SCHEME
C	-------------------------------------------------------------
C     INPUT VARIABLES IN ARGUMENT LIST
C	--------------------------------
C     MAXA(NEQ1)    = ADDRESSES OF DIAGONAL ELEMENTS IN A
C     A(ISTOR)      = STIFFNESS BLOCK TO BE REDUCED
C     B(ISTOR)      = COUPLING BLOCK
C     D(NEQ)        = DIAGONAL ELEMENTS OF GLOBAL STIFFNESS
C     V(NEQ)        = RIGHT-HAND-SIDE LOAD VECTOR
C     IND=1         = FLAG FOR FACTORISATION OF STIFFNESS MATRIX OR
C     IND=2         = REDUCTION AND BACK-SUBSTITUTION OF LOAD VECTOR
C     INDPD         = FLAG FOR FACTORISATION OF NON-POSITIVE DEFINITE A
C	---------------------------------
C     OUTPUT VARIABLES IN ARGUMENT LIST
C	---------------------------------
C     A(ISTOR)      = L AND FACTORS OF STIFFNESS BLOCK
C     D(NEQ)        = D FACTORS OF GLOBAL STIFFNESS MATRIX
C     V(NEQ)        = VECTOR OF GENERALIZED DISPLACEMENTS
C	-------------------------------------------
C     VARIABLES IN COMMON BLOCK /SOLU/ AND /INOU/
C	-------------------------------------------
C     NEQ           = NUMBER OF EQUATIONS (NEQ1=NEQ+1)
C     NBLOCK        = NUMBER OF STIFFNESS BLOCKS
C     ISTOR         = MAXIMUM LENGTH OF ANY ONE BLOCK
C     DETK          = DETERMINANT OF GLOBAL STIFFNESS MATRIX
C     NKFAC         = UNIT NUMBER OF RANDOM ACCESS FILE STORING
C                     FACTORIZED STIFFNESS BLOCKS
C     LSIGN         = STABILITY INDICATER (+1=STABEL, -1=INSTABEL)
C	---------------
C     LOCAL VARIABLES
C	---------------
C     NCOLUM        = NUMBER OF COLUMNS FOR THIS BLOCK
C     IBLO,JBLO     = LOOP COUNTER FOR REDUCTION/COUPLING BLOCK
C     ICOL,JCOL     = LOOP COUNTER FOR REDUCTION/COUPLING COLUMN
C     KBLO,KCOL     = BLOCK/COLUMN COUNTER
C     JBF,JBL       = FIRST,LAST COUPLED BLOCK
C     JCF,JCL       = FIRST,LAST COUPLED COLUMN
C     KDIA,KLOW,KUPP= ADDRESSES OF DIAGONAL,FIRST AND LAST ELEMENT IN
C                     COLUMN TO BE REDUCED
C     KHEI,KMOD,LCOP= TOTAL,EFFECTIVE HEIGHT AND NUM.OF COUPLED ELEM.
C     KTOP          = TOP ELEMENT TO BE REDUCED
C     LDIA,LHEI     = DIAGONAL ELEMENT AND HEIGHT OF COUPLING COLUMN
C     IEQ           = EQUATION NUMBER
C     -----------------------------------------------------------------
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
C      COMMON /ITER/ NSTEP,DLMAX,NPRIN,NDRAW,KONEQ,NIREF,ITEMAX,RTOL,
C     1              ITOPT,NUMREF,NUMITE,ITETOT,RHO,RHOP,ICONV,NOLIN,
C     2              LIMEQ(2),KSTEP,ETOL,RHOPREV
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
C	DIMENSION MAXA(1),A(1),B(1),D(1),V(1)
	DIMENSION MAXA(NEQ1),A(NWK),D(NEQ),V(NEQ)
C
	CALL CPU_TIME (TIM1)
	NCOLUM = NEQ
      IF (IND-2) 400,800,800
C     --------------------------------------------------------
C     FACTORIZE STIFFNESS MATRIX A (A = L*D*L'T DECOMPOSISION)
C     --------------------------------------------------------
C
C     -----------------------------------------
C     REDUCE BLOCK, LOOP OVER NUMBER OF COLUMNS
C     -----------------------------------------
 400  DO 690  ICOL=1,NCOLUM
      KDIA = MAXA(ICOL)
      KLOW = KDIA+1
      KUPP = MAXA(ICOL+1)-1
      KCOL = ICOL
      KHEI = KUPP-KLOW
      KMOD = MIN0(KHEI,ICOL-1)
      IF (KMOD) 600,500,410
C     ---------------------------------------------
C     LOOP OVER NUMBER OF COUPLING COLUMNS JCF, JCL
C     ---------------------------------------------
 410  LCOP = 0
      JCF  = ICOL-KMOD
      JCL  = ICOL-1
      KTOP = KLOW+KMOD
      IF (ICOL-1.LT.KHEI) LCOP = KHEI-ICOL+1
      DO 490  JCOL=JCF,JCL
      LCOP = LCOP+1
      KTOP = KTOP-1
      LDIA = MAXA(JCOL)
      LHEI = MAXA(JCOL+1)-1-LDIA
      IF (LHEI) 490,490,430
 430  NL = MIN0(LCOP,LHEI)
      C = 0.0
      DO 450  IL=1,NL
 450  C = C + A(LDIA+IL)*A(KTOP+IL)
      A(KTOP) = A(KTOP) - C
 490  CONTINUE
C     -------------------------------------------
C     FINAL COLUMN TERMS LIJ = GIJ/DII  AND
C     DIAGONAL TERMS     DJJ = KJJ - SUM(IRJ*GRJ)
C     -------------------------------------------
 500  IEQ = KCOL
      SUM = 0.0
      DO 590  KA=KLOW,KUPP
      IEQ = IEQ-1
      C = A(KA)/D(IEQ)
      SUM = SUM + C*A(KA)
 590  A(KA) = C
      A(KDIA) = A(KDIA) - SUM
C     ----------------------------------------------------------
C     SET DIAGONAL TERMS IN D AND TEST WHETHER POSITIVE DEFINITE
C     ----------------------------------------------------------
 600  D(KCOL) = A(KDIA)
      PIV = A(KDIA)
      IF (PIV-STOL)       610,610,690
 610  IF (DABS(PIV)-STOL) 650,650,620
 620  IF (INDPD-1)        630,690,690
 630  CALL ERRORS (12,KCOL,PIV,' SOLUTION ')
 650  IF (INDPD-1)       630,630,660
 660  PIV = STOL
      IF (PIV.EQ.0.0) PIV = -1.0D-16
      D(KCOL) = PIV
      A(KDIA) = PIV
 690  CONTINUE
C     ------------------------------------------------
C     END OF LOOP
C     CALCULATE DETERMINANT OF STIFFNESS MATRIX (DETK)
C     ------------------------------------------------
 695  CONTINUE !IF (ITASK.EQ.1) WRITE (NPLOT) (A(I),I=1,ISTOR)
 700	IF (ITASK-2) 705,790,795
C
 705  IF (KSTEP.GT.1) GOTO 750
      DETK=1.0
      GOTO 790
 750  DETK=1.0
      POWER = 1.0/DFLOAT(NEQ)
      DO 720 IEQ=1,NEQ
C
      DETK=DETK*DABS(D(IEQ))**power
 720  IF (D(IEQ).LT. 0.0) DETK = -DETK
 790  NFAC = NFAC+1
C
	CALL CPU_TIME (TIM2)
      TIM(14) = TIM(14) + (TIM2-TIM1)
 795  RETURN
C     --------------------------------
C     FORWARD REDUCTION OF LOAD VECTOR
C     --------------------------------
C
C     ------------------------------------------------
C     LOOP OVER NUMBER OF COLUMNS, VI = RI-SUM(IRI*VR)
C     ------------------------------------------------
 800	DO 880  ICOL=1,NCOLUM
      KLOW = MAXA(ICOL)+1
      KUPP = MAXA(ICOL+1)-1
      IF (KUPP-KLOW) 880,830,830
 830  IEQ = ICOL
      KV = IEQ
      C = 0.0
      DO 850  KA=KLOW,KUPP
      KV = KV-1
 850  C = C + A(KA)*V(KV)
      V(IEQ) = V(IEQ) - C
 880  CONTINUE
C     -----------------------------
C     BACK SUBSTITUTION , V=[D]-1*V
C     -----------------------------
 890	DO 900  IEQ=1,NEQ
 900  V(IEQ) = V(IEQ)/D(IEQ)
C     ---------------------------
C     LOOP OVER NUMBER OF COLUMNS
C     ---------------------------
 910	KCOL = NCOLUM
      DO 980  ICOL=1,NCOLUM
      KLOW = MAXA(KCOL)+1
      KUPP = MAXA(KCOL+1)-1
      IF (KUPP-KLOW) 980,930,930
 930  IEQ = KCOL
      KV = IEQ
      DO 950  KA=KLOW,KUPP
      KV = KV-1
 950  V(KV) = V(KV) - A(KA)*V(IEQ)
 980  KCOL = KCOL-1
C
 990  IF (ITASK.GT.2) RETURN
      NRED = NRED+1
C
	CALL CPU_TIME (TIM2)
      TIM(15) = TIM(15) + (TIM2-TIM1)
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ADDBAN (LM,MAXA,S,NEQF,NEQL,NPRE,NEF,NWS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ASSEMBLES UPPER TRIANGULAR ELEMENT STIFFNESS INTO COMPACTED
C     GLOBAL STIFFNESS BLOCK (CONTRIBUTIONS BETWEEN NEQF AND NEQL)
C	------------------------------------------------------------
C     LM(NEF)       = EQUATION NUMBERS FOR ELEMENT DEGREE OF FREEDOMS
C     MAXA(NEQ1)    = ADDRESSES OF DIAGONAL ELEMENTS IN A
C     A(ISTOR)      = GLOBAL COMPACTED STIFFNESS BLOCK
C     S(NWS)        = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     NEQF,NEQL     = FIRST AND LAST EQUATION CONTAINED IN BLOCK
C     NPRE          = NUMBER OF PREVIOUS ELEMENTS IN A
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ----------------------------------------------------------------
C
      DIMENSION LM(1),MAXA(1),S(1),KK(3)

	CALL MESTIF(LM,S,NWS,NEF,'WRT')  ! STORAGE STIFFNESS
	RETURN

C
      NDI = 0
      DO 200 IEF=1,NEF
      II = LM(IEF)
      IF (II.LT.NEQF .OR. II.GT.NEQL) GOTO 200
      MI = MAXA(II) - NPRE
      KS = IEF
      DO 220  JEF=1,NEF
      JJ = LM(JEF)
      IF (JJ) 220,220,110
 110  IJ = II-JJ
      IF (IJ) 220,210,210
 210  KK = MI+IJ
      KSS = KS
      IF (JEF.GE.IEF) KSS = JEF+NDI
C	A(KK) = A(KK)+S(KSS)
 220  KS = KS+NEF-JEF
 200  NDI = NDI+NEF-IEF

C
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ADDDIA (LM,DIA,S,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ASSEMBLES DIAGONAL TERM OF STIFFNESS 
C     ----------------------------------------------------------------
C
      DIMENSION LM(1),DIA(1),S(1)

	K = 0
	DO I = 1,NEF
	IEQ = LM(I)
	DO J = I,NEF
	K = K + 1
	IF(I.EQ.J.AND.IEQ.NE.0) DIA(IEQ) = DIA(IEQ) + S(K)
	ENDDO
	ENDDO


      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE MOVLEV(NT)
   	IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN
	COMMON /GiDEle/ LGID
      COMMON /OFFS/ NOPS,NHIGE

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT

C	ADDED BY SONGSAK APR2006 REACTION	
	COMMON /REACT/ LRC,LRCT,MFQ,LRID

C	==================================================================
C	COMMON BLOCK FOR HEAT SONGSAK MAR2007
	COMMON /SHEAT/ NHAEF,LHEAT1,LHEAT2,LHEAT3,LHEAT4
C	==================================================================
C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN
C	==================================================================

      COMMON A(9000000),IA(9000000)
C
      IF (NT .EQ. 2) GOTO 70

      IA(LNU)    = NAME(1)
      IA(LNU+1)  = NAME(2)
      IA(LNU+2)  = ITYPE
      IA(LNU+3)  = ISTYP
      IA(LNU+4)  = NLOPT
      IA(LNU+5)  = MTMOD
      IA(LNU+6)  = NSINC
      IA(LNU+7)  = ITOLEY
      IA(LNU+8)  = NELE
      IA(LNU+9)  = NMPS
      IA(LNU+10) = NGPS
      IA(LNU+11) = NMP
      IA(LNU+12) = NGP
      IA(LNU+13) = NNM
      IA(LNU+14) = NEX
      IA(LNU+15) = NCO
      IA(LNU+16) = NNF
      IA(LNU+17) = NEF
      IA(LNU+18) = NWG
      IA(LNU+19) = NPT
      IA(LNU+20) = NWA
      IA(LNU+21) = NWS

      IA(LNU+22) = NGR
      IA(LNU+23) = NGS
      IA(LNU+24) = NGT

      IA(LNU+25) = LNU
      IA(LNU+26) = LMP
      IA(LNU+27) = LGP
      IA(LNU+28) = LMS
      IA(LNU+29) = LGS
      IA(LNU+30) = LCO
      IA(LNU+31) = LEX
      IA(LNU+32) = LLM
      IA(LNU+33) = LCN

      IA(LNU+34) = LES
      IA(LNU+35) = LEC
      IA(LNU+36) = LED
      IA(LNU+37) = LEI
C	NEXT LINE CHANGED TO NEXT TWO BY GILSON - AUG2002
      IA(LNU+38) = LXI
      IA(LNU+39) = LEE

C	NEXT 5 LINES, ADDED BY GILSON - OCT2003 (FRAME ELEM)
	IA(LNU+40) = LOP
	IA(LNU+41) = LOS
	IA(LNU+42) = LSS

C	NEXT LINE ADDED BY GILSON - MARCH2004 (BEAM POST)
	IA(LNU+43) = LHG
	IA(LNU+44) = LHGN
	IA(LNU+45) = NHIGE

	IA(LNU+46) = LGID
	IA(LNU+47) = NOPS

	IA(LNU+48) = ISECT

	IA(LNU+49) = LRC      !ADDED BY SONGSAK APR2006 REACTION

	IA(LNU+50) = NHAEF    !ADDED BY SONGSAK MAR2007 
	IA(LNU+51) = LHEAT2   !ADDED BY SONGSAK MAR2007 
	IA(LNU+52) = LHEAT3   !ADDED BY SONGSAK MAR2007 
	IA(LNU+53) = LHEAT4   !ADDED BY SONGSAK MAR2007 

	IA(LNU+54) = LMV		!ADDED BY GILSON  MAY2007 KOREA  COMPOSITE
	IA(LNU+55) = LXI		!ADDED BY GILSON  MAY2007 KOREA  COMPOSITE


	IA(LNU+56) = NTEND	!TENDON DATA SONGSAK APR2009
	IA(LNU+57) = LTDN		!TENDON DATA SONGSAK APR2009

      IA(LNU+58) = NNO		!NUMBER OF ELEMENT NODE = NNM   (SEE ELEMIN)

      IA(LNU+59) = NEFC		!NUMBER OF ELEMENT FACE = NEFC  (SEE ELEMIN)
      IA(LNU+60) = LSFC		!STORE SUPPORT FLAG FOR EACH ELEMENT FACE (SEE SHELL 3 NODE ONATE)
      
      GOTO 80     

70    NAME(1)= IA(LNU)
      NAME(2)= IA(LNU+1)
      ITYPE  = IA(LNU+2)
      ISTYP  = IA(LNU+3)
      NLOPT  = IA(LNU+4)
      MTMOD  = IA(LNU+5)
      NSINC  = IA(LNU+6)
      ITOLEY = IA(LNU+7)
      NELE   = IA(LNU+8)
      NMPS   = IA(LNU+9)
      NGPS   = IA(LNU+10)
      NMP    = IA(LNU+11)
      NGP    = IA(LNU+12)
      NNM    = IA(LNU+13)
      NEX    = IA(LNU+14)
      NCO    = IA(LNU+15)
      NNF    = IA(LNU+16)
      NEF    = IA(LNU+17)
      NWG    = IA(LNU+18)
      NPT    = IA(LNU+19)
      NWA    = IA(LNU+20)
      NWS    = IA(LNU+21) 

      NGR    = IA(LNU+22)
      NGS    = IA(LNU+23)
      NGT    = IA(LNU+24)

      LNU    = IA(LNU+25)
      LMP    = IA(LNU+26)
      LGP    = IA(LNU+27)
      LMS    = IA(LNU+28)
      LGS    = IA(LNU+29)
      LCO    = IA(LNU+30)
      LEX    = IA(LNU+31)
      LLM    = IA(LNU+32)
      LCN    = IA(LNU+33)

      LES    = IA(LNU+34)
      LEC    = IA(LNU+35)
      LED    = IA(LNU+36)
      LEI    = IA(LNU+37)
C	NEXT LINE CHANGED TO NEXT TWO BY GILSON - AUG2002
      LXI    = IA(LNU+38)
      LEE    = IA(LNU+39)

C	NEXT 5 LINES, ADDED BY GILSON - OCT2003 (FRAME ELEM)
	LOP    = IA(LNU+40)
	LOS    = IA(LNU+41)
	LSS    = IA(LNU+42)

C	NEXT LINE ADDED BY GILSON - MARCH2004 (BEAM POST)
	LHG	 = IA(LNU+43) 
	LHGN	 = IA(LNU+44)
	NHIGE  = IA(LNU+45)
  
	LGID   = IA(LNU+46)
	NOPS   = IA(LNU+47)

	ISECT  = IA(LNU+48)

	LRC    = IA(LNU+49)      !ADDED BY SONGSAK APR2006 REACTION

	NHAEF  = IA(LNU+50)     !ADDED BY SONGSAK MAR2007 
	LHEAT2 = IA(LNU+51)     !ADDED BY SONGSAK MAR2007 
	LHEAT3 = IA(LNU+52)     !ADDED BY SONGSAK MAR2007 
	LHEAT4 = IA(LNU+53)     !ADDED BY SONGSAK MAR2007 

	LMV	 = IA(LNU+54)      !ADDED BY GILSON MAY2007 KOREA   COMPOSITE
	LXI    = IA(LNU+55)      !ADDED BY GILSON MAY2007 KOREA   COMPOSITE

	NTEND  = IA(LNU+56) 	!TENDON DATA SONGSAK APR2009
	LTDN   = IA(LNU+57) 	!TENDON DATA SONGSAK APR2009

      NNO	 = IA(LNU+58)	!NUMBER OF ELEMENT NODE = NNM   (SEE ELEMIN)

      NEFC   = IA(LNU+59)   !NUMBER OF ELEMENT FACE = NEFC  (SEE ELEMIN)
      LSFC   = IA(LNU+60)   !STORE SUPPORT FLAG FOR EACH ELEMENT FACE (SEE SHELL 3 NODE ONATE)
      
80    CONTINUE


	RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ROTCOM(VA,VB,VV)
	IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C     --------------------------------------------------------------
C     COMBINES TWO EULER ROTATION VECTORS VA AND VB TO GIVE A SINGLE
C     EQUIVALENT ROTATION VECTOR VV ( ALL COMPONENTS ARE IN RADS)
C     --------------------------------------------------------------
      DIMENSION VA(3),VB(3),VV(3),AV(3),BV(3),CV(3)
C
      ARC=6.2831853071796
      CALL ADDVEC (VA,VB,VV)
      CALL SCALEN (VA,AV,RA,3)
      CALL SCALEN (VB,BV,RB,3)
      IF (RA*RB.EQ.0.) RETURN
      CALL SCAPRD (AV,BV,AB,3)
      IF (DABS(AB).GE.1.) RETURN
      RA=.5*RA
      RB=.5*RB
      SNA=SIN(RA)
      SNB=SIN(RB)
      CSA=COS(RA)
      CSB=COS(RB)
      IF (CSA.EQ.0.) CSA=1.E-20
      IF (CSB.EQ.0.) CSB=1.E-20
      TNA=SNA/CSA
      TNB=SNB/CSB
      CALL SCALEN (VV,CV,RO,3)
      CALL SCAVEC (AV,AV,TNA,3)
      CALL SCAVEC (BV,BV,TNB,3)
      CALL SCAPRD (AV,BV,F,3)
      IF (F.EQ.1.) F=.999999999999
      F=1.0D0/(1.0D0-F)
      CALL VECPRD (AV,BV,VV)
      CALL ADDVEC (AV,BV,BV)
      CALL SUBVEC (BV,VV,VV)
      CALL SCAVEC (VV,VV,F,3)
      CALL SCALEN (VV,VV,F,3)
      RR=2.0D0*ATAN(F)
      CALL SCAPRD (CV,VV,CS,3)
      IF (CS.GE.0.) GOTO 20
      RR=-RR
      CALL SCAVEC (VV,VV,-1.,3)
   20 RD=RO-RR
      N=(RD+SIGN(ARC/2.01,RD))/ARC
      RR=RR+N*ARC
      CALL SCAVEC (VV,VV,RR,3)
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE ADDROT (ID,DISP,DISPI,DISLI,MSF,ISO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------
C     UPDATES NODAL ROTATIONS (FINITE ROTATION THEORY)
C	------------------------------------------------
C     ID(NSF,NSN) = EQUATION NUMBERS
C     DISP(NEQ)   = DISPLACEMENT VECTOR AT START OF STEP
C     DISPI(NEQ)  = ACCUMULATED DISPLACEMENT INCREMENT
C     DISLI(NEQ)  = CURRENT DISPLACEMENT INCREMENT
C     MSF         = MAXIMUM NUMBER OF NODAL D.O.F.
C     --------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
C
      DIMENSION ID(MSF,*),DISP(*),DISPI(*),DISLI(*),VA(3),VB(3),VV(3)
C
      DO 100 ISN=1,NSN
      CALL CLEARA (VA,3)
      CALL CLEARA (VB,3)
      DO 20 ISF=4,6
C	------------------------------------------------------
	IFF = 0
	DO II=1,9
	IF(IDOF(II).EQ.ISF) IFF = II
	ENDDO
	IEQ = 0
      IF(IFF.NE.0) IEQ = ID(IFF,ISN)
C	IEQ=ID(ISF,ISN)   
C	------------------------------------------------------
      IF (IEQ.EQ.0) GOTO 20
      L=ISF-3
      VA(L)=DISP(IEQ)+DISPI(IEQ)-DISLI(IEQ)
      VB(L)=DISLI(IEQ)
   20 CONTINUE
      CALL ROTCOM (VA,VB,VV)
      DO 60 ISF=4,6
C	------------------------------------------------------
	IFF = 0
	DO II=1,9
	IF(IDOF(II).EQ.ISF) IFF = II
	ENDDO
	IEQ = 0
      IF(IFF.NE.0) IEQ = ID(IFF,ISN)
C	IEQ=ID(ISF,ISN)   
C	------------------------------------------------------

      IF (IEQ.EQ.0) GOTO 60
      DISPI(IEQ)=VV(ISF-3)-DISP(IEQ)
   60 CONTINUE
  100 CONTINUE
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MOTRAN (ID,DISP,DISPI,DISLI,MSF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     EXACT TRANSFORMATION OF CONSERVATIVE NODAL ROTATION VECTORS
C     (SEMITANGENTIAL MODEL)
C	----------------------
C     ID(NSF,NSN)  = EQUATION NUMBERS
C     DISP(NEQ)    = DISPLACEMENT VECTOR AT START OF STEP
C     DISPI(NEQ)   = ACCUMULATED DISPLACEMENT INCREMENT
C     DISLI(NEQ)   = CURRENT DISPLACEMENT INCREMENT
C     MSF          = MAXIMUM NUMBER OF NODAL D.O.F.
C     --------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
C
      DIMENSION ID(MSF,*),DISP(*),DISPI(*),DISLI(*),ROT(3),ROV(3)
      DIMENSION VM(3),VMT(3),FC(6),FS(3),F(3)
	

      DO 100 ISN=1,NSN  !LOOP OVER NODE

      CALL CLEARA (ROT,3)
      CALL CLEARA (VM,3)

      DO 20 ISF=4,6     !LOOP OVER ROTATION D.O.F    
	
C	------------------------------------------------------
	IFF = 0
	DO II=1,9
	IF(IDOF(II).EQ.ISF) IFF = II
	ENDDO

	IEQ = 0
      IF(IFF.NE.0) IEQ = ID(IFF,ISN)
C	IEQ=ID(ISF,ISN)   !CALL EQUATION NUMBER
C	------------------------------------------------------

      IF (IEQ.EQ.0) GOTO 20 !AT THE BOUNDARY CONDITION

      L      = ISF-3                !BACK TO NUMBER 1-3 BUT STIL ROTATION D.O.F
      ROT(L) = DISP(IEQ)+DISPI(IEQ) !CURRENT TOTAL ROTATION
      VM(L)  = DISLI(IEQ)           !CURRENT INCREMENTAL ROTATION

   20 CONTINUE

      CALL SCAPRD (VM,VM,DUM,3)     !INCREMENTAL ROTATION (SCALAR)
      IF (DUM.EQ.0.) GOTO 100       !IF NO INCREMENTAL ROTATION
      CALL SCALEN (ROT,ROV,RO,3)    !TOTAL ROTATION PSUEDO VECTOR 

      CS = COS(RO)                  !COSINE
      SS = 0.5*SIN(RO)			  !SINE
      CC = 0.5*(1.0 - CS)           
 
      DO 30 I=1,3
   30 FS(I) = ROV(I)*SS

      N=1
      DO 40 I = 1,3
      F(I)  = ROV(I)*CC
      DO 40 J = I,3
      FC(N) = F(I)*ROV(J)
   40 N=N+1

      VMT(1)= (FC(4)+FC(6)+CS)*VM(1)-(FC(2)+FS(3))*VM(2)
     1       -(FC(3)-FS(2))*VM(3)
      VMT(2)=-(FC(2)-FS(3))*VM(1)+(FC(1)+FC(6)+CS)*VM(2)
     1       -(FC(5)+FS(1))*VM(3)
      VMT(3)=-(FC(3)+FS(2))*VM(1)-(FC(5)-FS(1))*VM(2)
     1       +(FC(1)+FC(4)+CS)*VM(3)

      DO 60 ISF=4,6
C	------------------------------------------------------
	IFF = 0
	DO II=1,9
	IF(IDOF(II).EQ.ISF) IFF = II
	ENDDO

	IEQ = 0
      IF(IFF.NE.0) IEQ = ID(IFF,ISN)
C	IEQ=ID(ISF,ISN)   
C	------------------------------------------------------
      IF (IEQ.EQ.0) GOTO 60
      DISLI(IEQ)=VMT(ISF-3)  !NEW INCREMENTAL ROTATION
   60 CONTINUE

  100 CONTINUE

      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FMVECB (X1,Y1,Z1,X2,Y2,Z2,VR,VS,VT,ELN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	=====================================================================
	DIMENSION VR(3),VS(3),VT(3)


	DX = X2-X1
	DY = Y2-Y1
	DZ = Z2-Z1

	ELN = SQRT(DX*DX + DY*DY + DZ*DZ)

	
	VR = 0.0D0 
	VS = 0.0D0
	VT = 0.0D0

	VR(1) = DX/ELN
	VR(2) = DY/ELN
	VR(3) = DZ/ELN

	CALL FMVEVR (VR,VS,VT)

	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FMVEVR (VR,VS,VT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV
C	=====================================================================
	DIMENSION VR(3),VS(3),VT(3),VRP(3),VG(3)


	VG = 0.0D0
	VG(NGRAV) = 1.0D0
	DOT = VG(1)*VR(1) + VG(2)*VR(2) + VG(3)*VR(3)
	DOT = ABS(1.0D0-ABS(DOT))
	TOL = 1.0E-8

	IF(DOT.GE.TOL) THEN
	VS = VG
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
	CALL VECPRD(VT,VR,VS)
	CALL SCALEN(VS,VS,DUM,3)
	ENDIF

	IF(DOT.LT.TOL) THEN
	SELECTCASE (NGRAV)
	CASE(1) 
	VS(1:3) = [0.0,-1.0,0.0]
	CASE(2)
	VS(1:3) = [-1.0,0.0,0.0]
	CASE(3)
	VS(1:3) = [-1.0,0.0,0.0]
	ENDSELECT
	CALL VECPRD(VR,VS,VT)
	CALL SCALEN(VT,VT,DUM,3)
	CALL VECPRD(VT,VR,VS)
	CALL SCALEN(VS,VS,DUM,3)
	ENDIF

	DO I = 1,3
	IF(ABS(VR(I)).LT.1.0E-6) VR(I) = 0.0D0
	IF(ABS(VS(I)).LT.1.0E-6) VS(I) = 0.0D0
	IF(ABS(VT(I)).LT.1.0E-6) VT(I) = 0.0D0
	ENDDO



	RETURN

	END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE ROMBAC (VR,VS,VT,RANG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION ROMAT(3,3),ROVCT(3,3),VECTM(3,3)

	RANG = RANG*0.01745329252  !ROTATION ANGLE

C	ROTATION MATRIX 
	ROMAT(1,1) = 1.0
	ROMAT(2,1) = 0.0
	ROMAT(3,1) = 0.0
	ROMAT(1,2) = 0.0
	ROMAT(2,2) = COS(RANG)
	ROMAT(3,2) = SIN(RANG)
	ROMAT(1,3) = 0.0
	ROMAT(2,3) =-SIN(RANG)
	ROMAT(3,3) = COS(RANG)

C	LOCAL VECTOR MATRIX
	ROVCT(1,1) = VR(1)
	ROVCT(2,1) = VR(2)
	ROVCT(3,1) = VR(3)
	ROVCT(1,2) = VS(1)
	ROVCT(2,2) = VS(2)
	ROVCT(3,2) = VS(3)
	ROVCT(1,3) = VT(1)
	ROVCT(2,3) = VT(2)
	ROVCT(3,3) = VT(3)

	VECTM = MATMUL(TRANSPOSE(ROMAT),TRANSPOSE(ROVCT))

C	NEW LOCAL BASE VECTOR
	VR(1) = VECTM(1,1)
	VR(2) = VECTM(1,2)
	VR(3) = VECTM(1,3)
	VS(1) = VECTM(2,1)
	VS(2) = VECTM(2,2)
	VS(3) = VECTM(2,3)
	VT(1) = VECTM(3,1)
	VT(2) = VECTM(3,2)
	VT(3) = VECTM(3,3)

	DO I = 1,3
	IF(ABS(VR(I)).LT.1.0E-6) VR(I) = 0.0D0
	IF(ABS(VS(I)).LT.1.0E-6) VS(I) = 0.0D0
	IF(ABS(VT(I)).LT.1.0E-6) VT(I) = 0.0D0
	ENDDO

	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE ADDSPIG(ID,MAXA,SD,SD2,SDWOK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

C	NEW SPRING SUPPORT SONGSAK JUL2006
	COMMON /SPBC/ NSS,NLSS

	DIMENSION ID(NSF,NSN),SD(6,NSN),MAXA(1)
	DIMENSION SD2(42,NSS),SDWOK(30,NSS)
	DIMENSION ST2(7,6),SDW(5,6)


C	DO II = 1,NSN
C	DO JJ = 1,NSF
C	IF (JJ.LE.6) THEN
C	IEQ   = ID(JJ,II)
C	SPFAC = SD(JJ,II)
C	IF(IEQ.GT.0) THEN
C	MI = MAXA(IEQ)
C	A(MI) = A(MI) + SPFAC
C	ENDIF
C	ENDIF
C	ENDDO
C	ENDDO


	DO ISS =1,NSS+NLSS
	CALL CALFOS(SD2(1,ISS),SDWOK(1,ISS),ST2,SDW)
	DO J  = 1,NSF
	JJ = IDOF(J)
	IF (JJ.LE.6) THEN
		II = ST2(1,JJ)
		IEQ   = ID(JJ,II)
		YOUNG = SD(JJ,II)

	IF(IEQ.GT.0) THEN

			YOUN2 = ST2(2,JJ)
			KFLAG = SDW(5,JJ)
			SPFAC = YOUNG

			IF(KFLAG.EQ.2) THEN
			SPFAC = YOUN2
			ENDIF

	MI = MAXA(IEQ)

C	A(MI) = A(MI) + SPFAC
	CALL MESTIF(IEQ,SPFAC,1,1,'WRT')

	ENDIF

	ENDIF
	ENDDO
	ENDDO



	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FOCSPIG(ID,SD,DISP,DISPI,RE,ITASK,NEQ,SD2,SDWOK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

C	NEW SPRING SUPPORT SONGSAK JUL2006
	COMMON /SPBC/ NSS,NLSS

	DIMENSION ID(NSF,NSN),SD(6,NSN),DISP(NEQ),DISPI(NEQ),DIP(NEQ)
	DIMENSION RE(NEQ)
	DIMENSION SD2(42,NSS),SDWOK(30,NSS)
	DIMENSION ST2(7,6),SDW(5,6)

	DO IEQ = 1,NEQ
	DIP(IEQ) = DISP(IEQ)
	IF(ITASK.EQ.2) DIP(IEQ) = DISP(IEQ)+DISPI(IEQ) 
	ENDDO


C	DO II = 1,NSN
C	DO JJ = 1,NSF
C	IF (JJ.LE.6) THEN
C	IEQ   = ID(JJ,II)
C	SPFAC = SD(JJ,II)
C	IF(IEQ.GT.0) THEN
C	DD = DIP(IEQ)
C	RE(IEQ) = RE(IEQ) - DD*SPFAC
C	ENDIF
C	ENDIF
C	ENDDO
C	ENDDO



	DO ISS =1,NSS+NLSS
	CALL CALFOS(SD2(1,ISS),SDWOK(1,ISS),ST2,SDW)

		DO J  = 1,NSF
		JJ = IDOF(J)
		IF (JJ.LE.6) THEN
		II = ST2(1,JJ)
		IOPT2 = ST2(6,JJ)
		IEQ   = ID(JJ,II)
		YOUNG = SD(JJ,II)

	


			IF(IEQ.GT.0) THEN
			DD = DIP(IEQ)

			YOUN2 = ST2(2,JJ)
			FY    = ST2(3,JJ)
			DMAX  = ST2(4,JJ) 
			QFY   = ST2(5,JJ) 

			EPSP  = SDW(1,JJ)
			SIGP  = SDW(2,JJ) 
			EPSTN = SDW(3,JJ) 
			EPTTN = SDW(4,JJ) 


			IF(IOPT2.EQ.0) THEN
			RE(IEQ) = RE(IEQ) - DD*YOUNG
			KFLAG = 1
			ELSE
			CALL SPGPAS(SDW(1,JJ),SDW(2,JJ),DD,RE(IEQ),YOUNG,YOUN2,
 	1				    SDW(3,JJ),ST2(3,JJ),KFLAG,SDW(4,JJ)) 
			ENDIF

			SDW(5,JJ) = KFLAG 	
				
			ENDIF

		ENDIF
		ENDDO
	CALL CALFOS2(SD2(1,ISS),SDWOK(1,ISS),ST2,SDW)
	ENDDO


	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE CALFOS(SD2,SDWOK,ST2,SDW)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION SD2(7,6),SDWOK(5,6)
	DIMENSION ST2(7,6),SDW(5,6)


	ST2 = SD2
	SDW = SDWOK


	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE CALFOS2(SD2,SDWOK,ST2,SDW)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION SD2(7,6),SDWOK(5,6)
	DIMENSION ST2(7,6),SDW(5,6)


	SD2   = ST2
	SDWOK = SDW


	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE TRANLG(VR,VS,VT,TRANS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION TRANS(14,14),VR(3),VS(3),VT(3)
C	=====================================================================
	DO I = 1,14
	DO J = 1,14
	TRANS(I,J) = 0.0
	ENDDO
	ENDDO

	DO I = 1,3
	TRANS(I  ,1) = VR(I)
	TRANS(I  ,2) = VS(I)
	TRANS(I  ,3) = VT(I)
	TRANS(I+3,4) = VR(I)
	TRANS(I+3,5) = VS(I)
	TRANS(I+3,6) = VT(I)

	TRANS(I+7 ,8)  = VR(I)
	TRANS(I+7 ,9)  = VS(I)
	TRANS(I+7 ,10) = VT(I)
	TRANS(I+10,11) = VR(I)
	TRANS(I+10,12) = VS(I)
	TRANS(I+10,13) = VT(I)
	ENDDO

	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FXCONT(W1,W2,A,B,ELN,COEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION CA(6),CB(6),DA(6),DB(6),C(6),D(6),COEF(6)
C	=====================================================================
	
	ELN2 = ELN*ELN
	ELN3 = ELN*ELN2

      IF(A.EQ.B) THEN
	ACON = W1
	BCON = 0.0D0
      ELSE
	ACON = (A*W2-B*W1) / (A-B)
	BCON =   (W1-W2)   / (A-B)
	ENDIF

C	---------------------------------------------------------
	DM = A

	DM2 = DM*DM
	DM3 = DM*DM2
	DM4 = DM*DM3

	C(1) = DM*(2.0D0*ELN3-2.0D0*DM2*ELN+DM3)/2.0D0/ELN3

	C(2) = DM2*(6.0D0*ELN2-8.0D0*DM*ELN+3.0D0*DM2)/12.0D0/ELN2

	C(3) = DM3*(-2.0D0*ELN+DM)/2.0D0/ELN3

	C(4) = DM3*(-4.0D0*ELN+3.0D0*DM)/12.0D0/ELN2

	C(5) = DM*(-2.0D0*ELN+DM)/2.0D0/ELN

	C(6) = DM2/2.0D0/ELN

	D(1) = DM2*(8.0D0*DM3-15.0D0*DM2*ELN+10.0D0*ELN3)/20.0D0/ELN3

	D(2) = DM3*(6.0D0*DM2-15.0D0*DM*ELN+10.0D0*ELN2)/30.0D0/ELN2

	D(3) = DM4*(8.0D0*DM-15.0D0*ELN)/20.0D0/ELN3

	D(4) = DM4*(4.0D0*DM-5.0D0*ELN)/20.0D0/ELN2

	D(5) = DM2*(2.0D0*DM-3.0D0*ELN)/6.0D0/ELN

	D(6) = DM3/3.0D0/ELN



      IF(A.EQ.B) THEN
	C(1) = DM*(-4.0D0*DM*ELN+3.0D0*DM2)/2.0D0/ELN3 + 
	1(2.0D0*ELN3-2.0D0*DM2*ELN+DM3)/2.0D0/ELN3

	C(2) = DM2*(-8.0D0*ELN+6.0D0*DM)/12.0D0/ELN2 + 
	12.0D0*DM*(6.0D0*ELN2-8.0D0*DM*ELN+3.0D0*DM2)/12.0D0/ELN2

	C(3) = DM3*(1.0D0)/2.0D0/ELN3 + 3.0D0*DM2*(-2.0D0*ELN+DM)/2.0D0/ELN3

	C(4) = DM3*(3.0D0)/12.0D0/ELN2 + 
	13.0D0*DM2*(-4.0D0*ELN+3.0D0*DM)/12.0D0/ELN2

	C(5) = DM*(1.0D0)/2.0D0/ELN + (-2.0D0*ELN+DM)/2.0D0/ELN

	C(6) = 2.0D0*DM/2.0D0/ELN

	D(1) = DM2*(24.0D0*DM2-30.0D0*DM*ELN)/20.0D0/ELN3 + 
	12.0D0*DM*(8.0D0*DM3-15.0D0*DM2*ELN+10.0D0*ELN3)/20.0D0/ELN3

	D(2) = DM3*(12.0D0*DM-15.0D0*ELN)/30.0D0/ELN2 + 
	13.0D0*DM2*(6.0D0*DM2-15.0D0*DM*ELN+10.0D0*ELN2)/30.0D0/ELN2

	D(3) = DM4*(8.0D0)/20.0D0/ELN3 + 
	14.0D0*DM3*(8.0D0*DM-15.0D0*ELN)/20.0D0/ELN3

	D(4) = DM4*(4.0D0)/20.0D0/ELN2 + 
	14.0D0*DM3*(4.0D0*DM-5.0D0*ELN)/20.0D0/ELN2

	D(5) = DM2*(2.0D0)/6.0D0/ELN + 2.0D0*DM*(2.0D0*DM-3.0D0*ELN)/6.0D0/ELN

	D(6) = 3.0D0*DM2/3.0D0/ELN
	ENDIF
	
	
	
	DO I = 1,6
	CA(I) = C(I)
	DA(I) = D(I)
	ENDDO
C	---------------------------------------------------------
C	---------------------------------------------------------
	DM = B

	DM2 = DM*DM
	DM3 = DM*DM2
	DM4 = DM*DM3

	C(1) = DM*(2.0D0*ELN3-2.0D0*DM2*ELN+DM3)/2.0D0/ELN3

	C(2) = DM2*(6.0D0*ELN2-8.0D0*DM*ELN+3.0D0*DM2)/12.0D0/ELN2

	C(3) = DM3*(-2.0D0*ELN+DM)/2.0D0/ELN3

	C(4) = DM3*(-4.0D0*ELN+3.0D0*DM)/12.0D0/ELN2

	C(5) = DM*(-2.0D0*ELN+DM)/2.0D0/ELN

	C(6) = DM2/2.0D0/ELN

	D(1) = DM2*(8.0D0*DM3-15.0D0*DM2*ELN+10.0D0*ELN3)/20.0D0/ELN3

	D(2) = DM3*(6.0D0*DM2-15.0D0*DM*ELN+10.0D0*ELN2)/30.0D0/ELN2

	D(3) = DM4*(8.0D0*DM-15.0D0*ELN)/20.0D0/ELN3

	D(4) = DM4*(4.0D0*DM-5.0D0*ELN)/20.0D0/ELN2

	D(5) = DM2*(2.0D0*DM-3.0D0*ELN)/6.0D0/ELN

	D(6) = DM3/3.0D0/ELN

	DO I = 1,6
	CB(I) = C(I)
	DB(I) = D(I)
	ENDDO
C	---------------------------------------------------------

	COEF(1) = ( CB(1) - CA(1) )*ACON + ( DB(1) - DA(1) )*BCON
	COEF(2) = ( CB(2) - CA(2) )*ACON + ( DB(2) - DA(2) )*BCON
	COEF(3) = (-CB(3) + CA(3) )*ACON + (-DB(3) + DA(3) )*BCON
	COEF(4) = ( CB(4) - CA(4) )*ACON + ( DB(4) - DA(4) )*BCON

	COEF(5) = (-CB(5) + CA(5) )*ACON + (-DB(5) + DA(5) )*BCON
	COEF(6) = ( CB(6) - CA(6) )*ACON + ( DB(6) - DA(6) )*BCON


      IF(A.EQ.B) THEN
	COEF(1) = ( CA(1) )*ACON 
	COEF(2) = ( CA(2) )*ACON 
	COEF(3) = (-CA(3) )*ACON 
	COEF(4) = ( CA(4) )*ACON 

	COEF(5) = (-CA(5) )*ACON 
	COEF(6) = ( CA(6) )*ACON 
	ENDIF


	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FMCONT(W1,W2,A,B,ELN,COEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION CA(6),CB(6),DA(6),DB(6),C(6),D(6),COEF(6)
C	=====================================================================
	
	ELN2 = ELN*ELN
	ELN3 = ELN*ELN2


      IF(A.EQ.B) THEN
	ACON = W1
	BCON = 0.0D0
      ELSE
	ACON = (A*W2-B*W1) / (A-B)
	BCON =   (W1-W2)   / (A-B)
	ENDIF
	
C	---------------------------------------------------------
	DM = A

	DM2 = DM*DM
	DM3 = DM*DM2
	DM4 = DM*DM3

	C(1) = DM2*(2.0D0*DM-3.0D0*ELN)/ELN3

	C(2) = DM*(ELN2-2.0D0*DM*ELN+DM2)/ELN2

	C(3) = DM2*(2.0D0*DM-3.0D0*ELN)/ELN3

	C(4) = DM2*(-ELN+DM)/ELN2

	C(5) = DM*(-2.0D0*ELN+DM)/2.0D0/ELN

	C(6) = DM2/2.0D0/ELN

	D(1) = DM3*(3.0D0*DM-4.0D0*ELN)/2.0D0/ELN3

	D(2) = DM2*(9.0D0*DM2-16.0D0*DM*ELN+6.0D0*ELN2)/12.0D0/ELN2

	D(3) = DM3*(3.0D0*DM-4.0D0*ELN)/2.0D0/ELN3

	D(4) = DM3*(9.0D0*DM-8.0D0*ELN)/12.0D0/ELN2

	D(5) = DM2*(2.0D0*DM-3.0D0*ELN)/6.0D0/ELN

	D(6) = DM3/3.0D0/ELN



      IF(A.EQ.B) THEN
 	C(1) = DM2*(2.0D0)/ELN3 + 2.0D0*DM*(2.0D0*DM-3.0D0*ELN)/ELN3

	C(2) = DM*(-2.0D0*ELN+2.0D0*DM)/ELN2 + (ELN2-2.0D0*DM*ELN+DM2)/ELN2

	C(3) = DM2*(2.0D0)/ELN3 + 2.0D0*DM*(2.0D0*DM-3.0D0*ELN)/ELN3

	C(4) = DM2*(1.0D0)/ELN2 + 2.0D0*DM*(-ELN+DM)/ELN2

	C(5) = DM*(1.0D0)/2.0D0/ELN + (-2.0D0*ELN+DM)/2.0D0/ELN

	C(6) = 2.0D0*DM/2.0D0/ELN

	D(1) = DM3*(3.0D0)/2.0D0/ELN3 + 
	13.0D0*DM2*(3.0D0*DM-4.0D0*ELN)/2.0D0/ELN3

	D(2) = DM2*(18.0D0*DM-16.0D0*ELN)/12.0D0/ELN2 + 
	12.0D0*DM*(9.0D0*DM2-16.0D0*DM*ELN+6.0D0*ELN2)/12.0D0/ELN2

	D(3) = DM3*(3.0D0)/2.0D0/ELN3 + 
	13.0D0*DM2*(3.0D0*DM-4.0D0*ELN)/2.0D0/ELN3

	D(4) = DM3*(9.0D0)/12.0D0/ELN2 + 
	13.0D0*DM2*(9.0D0*DM-8.0D0*ELN)/12.0D0/ELN2

	D(5) = DM2*(2.0D0)/6.0D0/ELN + 2.0D0*DM*(2.0D0*DM-3.0D0*ELN)/6.0D0/ELN

	D(6) = 3.0D0*DM2/3.0D0/ELN     
      ENDIF
      
	DO I = 1,6
	CA(I) = C(I)
	DA(I) = D(I)
	ENDDO
C	---------------------------------------------------------
C	---------------------------------------------------------
	DM = B

	DM2 = DM*DM
	DM3 = DM*DM2
	DM4 = DM*DM3

	C(1) = DM2*(2.0D0*DM-3.0D0*ELN)/ELN3

	C(2) = DM*(ELN2-2.0D0*DM*ELN+DM2)/ELN2

	C(3) = DM2*(2.0D0*DM-3.0D0*ELN)/ELN3

	C(4) = DM2*(-ELN+DM)/ELN2

	C(5) = DM*(-2.0D0*ELN+DM)/2.0D0/ELN

	C(6) = DM2/2.0D0/ELN

	D(1) = DM3*(3.0D0*DM-4.0D0*ELN)/2.0D0/ELN3

	D(2) = DM2*(9.0D0*DM2-16.0D0*DM*ELN+6.0D0*ELN2)/12.0D0/ELN2

	D(3) = DM3*(3.0D0*DM-4.0D0*ELN)/2.0D0/ELN3

	D(4) = DM3*(9.0D0*DM-8.0D0*ELN)/12.0D0/ELN2

	D(5) = DM2*(2.0D0*DM-3.0D0*ELN)/6.0D0/ELN

	D(6) = DM3/3.0D0/ELN

	DO I = 1,6
	CB(I) = C(I)
	DB(I) = D(I)
	ENDDO
C	---------------------------------------------------------

	COEF(1) = ( CB(1) - CA(1) )*ACON + ( DB(1) - DA(1) )*BCON
	COEF(2) = ( CB(2) - CA(2) )*ACON + ( DB(2) - DA(2) )*BCON
	COEF(3) = (-CB(3) + CA(3) )*ACON + (-DB(3) + DA(3) )*BCON
	COEF(4) = ( CB(4) - CA(4) )*ACON + ( DB(4) - DA(4) )*BCON

	COEF(5) = (-CB(5) + CA(5) )*ACON + (-DB(5) + DA(5) )*BCON
	COEF(6) = ( CB(6) - CA(6) )*ACON + ( DB(6) - DA(6) )*BCON

      IF(A.EQ.B) THEN
	COEF(1) = ( CA(1) )*ACON 
	COEF(2) = ( CA(2) )*ACON 
	COEF(3) = (-CA(3) )*ACON 
	COEF(4) = ( CA(4) )*ACON 

	COEF(5) = (-CA(5) )*ACON 
	COEF(6) = ( CA(6) )*ACON 
	ENDIF
	
	
	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FIXHNG (RG,VR,ANG,ELN,LREAS,IHSET,IELE,NELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION VR(3),VS(3),VT(3),TRANS(7,7),TRANH(14,14)
	DIMENSION RM(7,3),RL(7,3),RG(7,3),RN(14),RNN(14)
	DIMENSION LREAS(14,1),IHSET(1),IPIN(14)

	
	RANG = ANG

	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

	TRANS = 0.0
	DO I = 1,3
	TRANS(I,1) = VR(I)
	TRANS(I,2) = VS(I)
	TRANS(I,3) = VT(I)
	TRANS(I+3,4) = VR(I)
	TRANS(I+3,5) = VS(I)
	TRANS(I+3,6) = VT(I)
	ENDDO

	RL = 0.0
	DO INM = 1,2

	DO IFF = 1,6
	DO JFF = 1,6
	RL(IFF,INM) = RL(IFF,INM) + TRANS(JFF,IFF)*RG(JFF,INM) !GLOBAL TO LOCAL
	ENDDO
	ENDDO

	ENDDO

	DO I = 1,7
	RN(I)   = RL(I,1)
	RN(I+7) = RL(I,2)
	ENDDO

C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
      IPIN(1:14) = 0
      IHET = IHSET(IELE)
      IF(IHET.NE.0) IPIN(1:14) = LREAS(1:14,IHET)
	CALL TRNHIG(TRANH,ELN,IPIN)


	RNN = MATMUL(TRANSPOSE(TRANH),RN)
	DO I = 1,7
	RM(I,1) = RNN(I)
	RM(I,2) = RNN(I+7)
	ENDDO
C	------------------------------------------------------------

	RG = 0.0
	DO INM = 1,2

	DO IFF = 1,6
	DO JFF = 1,6
	RG(IFF,INM) = RG(IFF,INM) + TRANS(IFF,JFF)*RM(JFF,INM) !BACK TO GLOBAL
	ENDDO
	ENDDO

	ENDDO


	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LAGROT (ROT,XYZ,DIS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION ROT(3),XYZ(3),DIS(3),SMAT1(3,3),SMAT2(3,3),RMAT(3,3)

	PI  = 3.1415926
	PI2 = 6.2831853
	

	RX  = ROT(1)*PI/180.0
	RY  = ROT(2)*PI/180.0
	RZ  = ROT(3)*PI/180.0

	ANG  = SQRT ( RX*RX + RY*RY + RZ*RZ )
	ANGO = ANG

	DIS(1:3) = 0.0D0
	

	IF(ANG.LE.1.0E-12) RETURN

	DO WHILE (ANG > PI2) 
	ANG = ANG - PI2
	ENDDO

	RX  = ANG*RX/ANGO
	RY  = ANG*RY/ANGO
	RZ  = ANG*RZ/ANGO

	ANG1 = SIN(ANG)/ANG
	ANG2 = (1.0-COS(ANG))/ANG**2.0

	ZR = 0.0D0 
	SMAT1(1,1:3) = [ ZR,-RZ, RY]
	SMAT1(2,1:3) = [ RZ, ZR,-RX]
	SMAT1(3,1:3) = [-RY, RX, ZR]

	SMAT2 = MATMUL(SMAT1,SMAT1)
	RMAT = ANG1*SMAT1 + ANG2*SMAT2

	DIS = MATMUL(RMAT,XYZ)



	RETURN
	END	
C	=====================================================================
C     SUBROUTINE OFFSHORE_WARN ( WARNING )
C	=====================================================================
	SUBROUTINE OFFSHORE_WARN	
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /WARNING/ WARNING,RAMDA(100),RK
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /OFFSHOREOUT/ UXMAX,UYMAX,NSTREAMFUNCTION,NFUNCTION
1	FORMAT ('')
2	FORMAT ('                     ***************************************')
3     FORMAT ('                     *                                     *')
4	FORMAT ('                     *       WARING OFFSHORE ANALYSIS      *')
5     FORMAT ('                     *                                     *')
6     FORMAT ('                     *        PLEASE CHECK BASE UNIT       *')
7     FORMAT ('                     *                                     *')
8	FORMAT ('                     ***************************************')
9     FORMAT ('                       LENGTH : m  MASS  : kg => FORCE : N  ')	
10    FORMAT ('                       LENGTH : cm  MASS : kg => FORCE : N ')	
11    FORMAT ('                       LENGTH : mm  MASS : kg => FORCE : N ')
12    FORMAT ('                       LENGTH : ft  MASS : slug => FORCE : lb')
13	FORMAT (' -- THE RESULT OFFSHORE ANALYSIS ---')
14    FORMAT ('    WAVE LENGTH         = ',F12.5)
15    FORMAT ('    WAVE NUMBER         = ',F12.5)
16    FORMAT ('    VELOCITY AT SURFACE = ',F12.5)
17    FORMAT ('    OFFSHORE PARAMETER',I3,'')
	   IF (WARNING.EQ.9.806D0)THEN
	WRITE (*,1)
	WRITE (*,2)
	WRITE (*,3)
	WRITE (*,4)
	WRITE (*,5)
	WRITE (*,6)
	WRITE (*,7)
	WRITE (*,8)
      WRITE (*,1)
      WRITE (*,9)
      WRITE (*,1)
         IF (NSTREAMFUNCTION.EQ.0.0D0)THEN
         DO I=1,LCS
         WRITE (*,17) I
         WRITE (*,14) RAMDA(I)
         ENDDO
         WRITE (*,15) RK
             IF (NFUNCTION.EQ.3)THEN
             WRITE (*,16) UXMAX
             ENDIF
         ELSEIF(NSTREAMFUNCTION.EQ.1.0D0)THEN
             IF (NFUNCTION.EQ.3)THEN
             WRITE (*,16) UXMAX
             ELSEIF (NFUNCTION.NE.3)THEN
             DO I=1,LCS
             WRITE (*,17) I
             WRITE (*,14) RAMDA(I)
             ENDDO
             WRITE (*,15) RK
             ENDIF
         ENDIF
      WRITE (*,1)
         ELSEIF (WARNING.EQ.9806D0)THEN
      WRITE (*,1)
	WRITE (*,2)
	WRITE (*,3)
	WRITE (*,4)
	WRITE (*,5)
	WRITE (*,6)
	WRITE (*,7)
	WRITE (*,8)
	WRITE (*,1)
      WRITE (*,10)
      WRITE (*,1)
         IF (NSTREAMFUNCTION.EQ.0.0D0)THEN
         DO I=1,LCS
         WRITE (*,17) I
         WRITE (*,14) RAMDA(I)
         ENDDO
         WRITE (*,15) RK
             IF (NFUNCTION.EQ.3)THEN
             WRITE (*,16) UXMAX
             ENDIF
         ELSEIF(NSTREAMFUNCTION.EQ.1.0D0)THEN
             IF (NFUNCTION.EQ.3)THEN
             WRITE (*,16) UXMAX
             ELSEIF (NFUNCTION.NE.3)THEN
             DO I=1,LCS
             WRITE (*,17) I
             WRITE (*,14) RAMDA(I)
             ENDDO
             WRITE (*,15) RK
             ENDIF
         ENDIF
      WRITE (*,1) 
         ELSEIF (WARNING.EQ.9806D0)THEN
      WRITE (*,1)
	WRITE (*,2)
	WRITE (*,3)
	WRITE (*,4)
	WRITE (*,5)
	WRITE (*,6)
	WRITE (*,7)
	WRITE (*,8)
	WRITE (*,1)
      WRITE (*,11)
      WRITE (*,1)
         IF (NSTREAMFUNCTION.EQ.0.0D0)THEN
         DO I=1,LCS
         WRITE (*,17) I
         WRITE (*,14) RAMDA
         ENDDO
         WRITE (*,15) RK
             IF (NFUNCTION.LT.3)THEN
             WRITE (*,16) UXMAX
             ENDIF
         ELSEIF(NSTREAMFUNCTION.EQ.1.0D0)THEN
             IF (NFUNCTION.LT.3)THEN
             WRITE (*,16) UXMAX
             ELSEIF (NFUNCTION.NE.3)THEN
             DO I=1,LCS
             WRITE (*,17) I
             WRITE (*,14) RAMDA(I)
             ENDDO
             WRITE (*,15) RK
             ENDIF
         ENDIF
      WRITE (*,1) 
	   ELSEIF (WARNING.EQ.32.172D0)THEN
      WRITE (*,1)
	WRITE (*,2)
	WRITE (*,3)
	WRITE (*,4)
	WRITE (*,5)
	WRITE (*,6)
	WRITE (*,7)
	WRITE (*,8)
	WRITE (*,1)
      WRITE (*,12)
      WRITE (*,1)
         IF (NSTREAMFUNCTION.EQ.0.0D0)THEN
         DO I=1,LCS
         WRITE (*,14) RAMDA(I)
         ENDDO
         WRITE (*,15) RK
             IF (NFUNCTION.EQ.3)THEN
             WRITE (*,16) UXMAX
             ENDIF
         ELSEIF(NSTREAMFUNCTION.EQ.1.0D0)THEN
             IF (NFUNCTION.EQ.3)THEN
             WRITE (*,16) UXMAX
             ELSEIF (NFUNCTION.NE.3)THEN
             DO I=1,LCS
             WRITE (*,17) I
             WRITE (*,14) RAMDA(I)
             ENDDO
             WRITE (*,15) RK
             ENDIF
         ENDIF
      WRITE (*,1)
	  ENDIF		
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================