C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CONCRETE
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

      DIMENSION NI(3)

C	---------------------------------------------------------------------
C     CONCRETE PARAMETERS MODULE  (BUILD CREEP COMPLIANCE TABLES)    -
C	---------------------------------------------------------------------                                         
      NCPT    = 1
	NI(1:3) = [32,48,4] !NA,NT,NR
      IPRN    = 0
	CALL FREBUF
	CALL FREECH
      CALL FREINT('N',NCPT,1)
      CALL FREINT('T',NI  ,3)
      CALL FREINT('P',IPRN,1)

      CALL DEFNINT('%XNI',MXNI,1,4)  !MAX--NA NT NR NCPC	

	NA = NI(1)
	NT = NI(2)
	NR = NI(3)

	CALL INTFILL ('%XNI',NA,1,1,1)
	CALL INTFILL ('%XNI',NT,1,2,1)
	CALL INTFILL ('%XNI',NR,1,3,1)

C	FILE FOR PRINTING  NOT=ISO
	CALL INTFILL('%IOL',NOT,1,3,0)

      WRITE(NOT,2001) NCPT,NA,NT,NR

      CALL DEFNREL('%FNC',KFNC,20,NCPT) !MECHANICAL DATA DATA INPUT FOR EACH MAT No.
	CALL RELZERO('%FNC')
      CALL DEFNINT('%INC',KINC,5,NCPT)  !MECHANICAL DATA DATA INPUT FOR EACH MAT No.
	CALL INTZERO('%INC')

      CALL DEFNINT('%CPT',KCPT,7,NCPT)  !DATA INPUT FOR EACH MAT No. (CREEP)
	CALL INTZERO('%CPT')
      CALL DEFNREL('%GAM',KGAM,1,NCPT)  !FIRST RETADATION TIME
	CALL RELZERO('%GAM')

      NCPC = 5 + NR + NR     
	CALL INTFILL('%XNI',NCPC,1,4,1)                                              
      CALL DEFNREL('%PRM',KPRM,NA*NCPC,NCPT)  !Parameter for Time Dependent Analysis
	CALL RELZERO('%PRM')


      DO 60 NPTP=1,NCPT
      IF(IPRN.GE.2) WRITE(NOT,2002) NPTP


      CALL CNCTYP (RCONDT(KFNC),ICONDT(KINC),
	1			 ICONDT(KCPT),RCONDT(KGAM),
     2			 NPTP,NPIN,IPRN,NI,NOT)

	CALL INTFILL('%CPT',KM  ,1,NPIN,0)
	CALL INTFILL('%CPT',NAT ,2,NPIN,0)
	CALL INTFILL('%CPT',NTT ,3,NPIN,0)
	CALL INTFILL('%CPT',NRT ,4,NPIN,0)
	CALL INTFILL('%CPT',KL  ,5,NPIN,0)
	CALL INTFILL('%CPT',IGEN,6,NPIN,0)

	IF(KM.EQ.4) GOTO 60  !LINEAR -NO TIME DEPENDENT DATA

C-----FIND CREEP COMPLIANCE COEFFICIENTS-----
      CALL DEFNREL('&LSC',KLSC,NTT,NRT)     !TEMP  Time matrix  G
      CALL DEFNREL('&CCF',KCCF,NRT,KL )     !TEMP  Creep Coef   A 
      CALL DEFNREL('&TIM',KTIM,1  ,NTT)     !TEMP  TIME AFTER LOADING (OBSERVATION TIME)
      CALL DEFNREL('&CRP',KCRP,NTT,2  )     !TEMP  Function J


	NUM = KPRM + NA*NCPC*(NPIN-1)
      CALL CRPCFT (RCONDT(KFNC),ICONDT(KINC),
	1			 ICONDT(KCPT),RCONDT(NUM),RCONDT(KTIM),RCONDT(KCRP),
     2             RCONDT(KLSC),RCONDT(KCCF),RCONDT(KGAM),
     3             IPRN,NA,NT,NR,NPIN,KM,NAT,NTT,NRT,KL,IGEN,NOT)



      CALL DELTREL('&LSC')
      CALL DELTREL('&CCF')
      CALL DELTREL('&TIM')
      CALL DELTREL('&CRP')

   60 CONTINUE

C	------------------------------------
C	FILE AND PRINT SUMMARY
C	------------------------------------
	IF(KM.EQ.4) RETURN  !LINEAR -NO TIME DEPENDENT DATA
      CALL CRPPRN (ICONDT(KCPT),RCONDT(KPRM),NA,NCPC,NCPT,IPRN,NOT)

      RETURN

 2000 FORMAT (1X,79('='),/,' =  C O N C R E T E   T I M E   P A R A M',
     1 ' E T E R   S P E C I F I C A T I O N  =',/,1X,79('='))
 2001 FORMAT (/,1X,'NUMBER OF CONCRETE TYPES  . .',I5,
     1        /,1X,'MAXIMUM   LOADING AGES  . . .',I5,
     2        /,1X,'MAXIMUM   CREEP TIMES . . . .',I5,
     3        /,1X,'MAXIMUM   RETARDATION TIMES .',I5,/,/)
 2002 FORMAT (/,' CONCRETE PARAMETER INPUT SET',I5,/,1X,33('-'),/)
C
      END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE CRPCFT (FCNC,ICNC,ICPT,PRAM,CTIM,CREP,ALSC,CCOF,GAM0,     
     1                   IPRN,NA,NT,NR,NPIN,KM,NAT,NTT,NRT,KL,IGEN,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	------------------------------------
C     BUILD CREEP COEFFICIENT TABLES    -
C	------------------------------------
C
C     ICPT = MODEL SIZE : TYPE,NAT,NRT,KL FOR EACH MODEL
C     PRAM = PARAMETERS ARRAY : TAU, EMOD, SHRN, A(I,J) . . .
C     CTIM = OBSERVATION TIMES
C     CREP = CREEP VALUES AT OBSERVATION TIMES
C     ALSC = LEAST SQUARES COEFFICIENT MATRIX
C     CCOF = CREEP COEFFICIENTS
C     GAM0 = SMALLEST RETARDATION TIMES
C     NA,NAT = NUMBER OF LOADING AGES (MAX, THIS SET)
C     NT,NTT = NUMBER OF OBSRVN TIMES (MAX, THIS SET)
C     NR,NRT = NUMBER OF RETRDN TIMES (MAX, THIS SET)
C     KM = MODEL TYPE (1=LAB,2=ACI,3=CEB, 4=LINEAR, 5=KOREAN, 6=CHINA JTG D62-2004)
C     IGEN = AGE/TIME INPUT FLAG (0=GENERATED ,1=INPUT)
C     KL = NUMBER OF LOADING TERMS
C
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)


      DIMENSION ICPT(7,1),PRAM(NA,1),CTIM(1),CREP(NTT,1),
     1          ALSC(NTT,NRT),CCOF(NRT,KL),GAM0(1)
     	DIMENSION ATA(NRT,NRT)                       !NRxNR LEAST SQUARE MATRIX
	DIMENSION FCNC(20,1),ICNC(5,1)               !CONC PROPERITES
	DIMENSION ACI(12),CEB(6),COR(6),CHN(6)       !DATA INPUT FOR ACI & CEB & KOREAN & CHINA

	

	CALL INTFILL ('UNIT',IUF,1,1,0)
	CALL INTFILL ('UNIT',IUL,1,2,0)
	CALL INTFILL ('UNIT',IUM,1,4,0)

	MM   = ICNC(1,NPIN) !0 = Linear, 1 = Concrete, 2 = Concrete with nl stress-strain
	FC28 = FCNC(1,NPIN)  !28 Day Strength
	R28  = FCNC(15,NPIN) !28 Day tens Strength
	CU   = FCNC(2,NPIN) !Ult. Creep Coef
	SU   = FCNC(3,NPIN) !Ult. Shrnk Strain
	WDN  = FCNC(4,NPIN) !WEIGHT DENSITY

	IF(MM.EQ.0) EE = FCNC(14,NPIN)
C	------------------------------------
C	READ IN THE MODEL PARAMETERS 
C	------------------------------------
      GO TO (11,12,13,100,14,15),KM  !LAB ACI CEB LIN KOR CHN

   11 GO TO 20

   12 CALL ACIINP (ACI,IPRN,NOT)                                             
      GO TO 16
                       
   13 CALL CEBINP (CEB,IPRN,NOT)                                     
      GO TO 16 
                       
   14 CALL KORINP (COR,IPRN,NOT)                                     
      GO TO 16                         
                       
   15 CALL CHNINP (CHN,IPRN,NOT)                                     
      GO TO 16   
C	------------------------------------
C	GENERATE LOADING AGES
C	------------------------------------
   16 IF(IGEN.EQ.0) THEN !AUTOMATIC GENERATE LOADING TIME
	TMIN = 1.0D0  !DAYS
	TMAX = 2.0E4  !DAYS
      CALL GENTIM (PRAM(1,1),TMIN,TMAX,NAT)
      IF(NAT.EQ.1) PRAM(1,1) = 28.0D+0
	ENDIF
C	------------------------------------
20	CONTINUE
C	------------------------------------
C	GENERATE PARAMETER TABLES FOR EACH LOADING AGE 
C	------------------------------------
	DO 95 NATN=1,NAT

      IF(IGEN.NE.0) THEN
C	MANNUAL INPUT OF TIME
      CALL FREBUF
      CALL FREREL('A',TAU ,1  )  !LOADING TIME (ONE VALUE)
      PRAM(NATN,1) = TAU
      CALL FREREL('T',CTIM,NTT)  !TIMES AFTER LOADING (NTT VALUES)
      GO TO 35
	ENDIF

C	------------------------------------
C	GENERATE TIMES AFTER LOADING 
C	------------------------------------
	TAU = PRAM(NATN,1)
	TMIN = TAU +0.0D0    !DAYS
	TMAX = TMIN+3.0E4    !DAYS
      CALL GENTIM (CTIM,TMIN,TMAX,NTT)
C	------------------------------------

C	------------------------------------
C	FIND CREEP CURVE EMOD SHRINKAGE FACTORS FOR THIS LOADING AGE
C	------------------------------------
   35 GO TO (41,42,43,100,44,45),KM

   41 CALL LABCRP (TAU,NTT,CREP,EMOD,SHRN,FPC,FPT)                   
      GO TO 46

   42 CALL ACICRP (TAU,CTIM,NTT,ACI,FC28,IUF,IUL,
	1			 CREP,EMOD,SHRN,FPC,FPT,CU,SU,
     2			 IUM,WDN,R28,MM)  
      GO TO 46

   43 CALL CEBCRP (TAU,CTIM,NTT,CEB,FC28,IUF,IUL,
	1			 CREP,EMOD,SHRN,FPC,FPT,CU,SU,MM)                     
      GO TO 46

   44 CALL KORCRP (TAU,CTIM,NTT,COR,FC28,IUF,IUL,
	1			 CREP,EMOD,SHRN,FPC,FPT,CU,SU,MM)                     
      GO TO 46

   45 CALL CHNCRP (TAU,CTIM,NTT,CHN,FC28,IUF,IUL,
	1			 CREP,EMOD,SHRN,FPC,FPT,CU,SU,MM)                     
      GO TO 46
	    	                       
C	------------------------------------
C	FORM THE TIME CONSTANT MATRIX
C	------------------------------------
   46 RTIM = GAM0(NPIN) / 10.0D+0
      DO 50 J=1,NRT
      RTIM = RTIM * 10.0D+0                           !RETARDATION TIME
      DO 50 I=1,NTT
   50 ALSC(I,J) = 1.0D+0 - DEXP((TAU-CTIM(I))/RTIM)   !TIME MATRIX

	IF(MM.EQ.0) EMOD = EE  !LINEAR WITH SPECIFIC ELASTIC MODULUS
C	------------------------------------
C	PERFORM LEAST SQUARES SOLUTION 
C	------------------------------------
      CALL LSTSQR (ALSC,CREP,CCOF,ATA,NRT,NTT,KL,NOT)

C	------------------------------------
C	ESTIMATE ERROR, PRINT IT ALL OUT
C	------------------------------------
      IF(IPRN.GE.2)
     1CALL CRPERR (TAU,ALSC,CCOF,CREP,CTIM,
     2             EMOD,SHRN,FPC,FPT,
     3			 NRT,NTT,KL,NATN,NOT)

C	------------------------------------
C	STORE INTO BIG ARRAY
C	------------------------------------
60	PRAM(NATN,2) = EMOD
      PRAM(NATN,3) =-SHRN  !MINUS FOR CONTRACTION
      PRAM(NATN,4) = FPC
      PRAM(NATN,5) = FPT

	
      II = 5                                                               
      DO 80 K=1,KL
      DO 80 J=1,NRT
      II = II + 1
	PRAM(NATN,II) = CCOF(J,K)
80	CONTINUE
C	------------------------------------

95	CONTINUE


100   RETURN
C
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE GENTIM (TIM,TS,TF,NTIM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------
C     GENERATE LOG-LINEAR TIMES BETWEEN TS AND TF IN TIM  -
C	----------------------------------------------------------
      DIMENSION TIM(NTIM)
C
      IF(NTIM.GT.1) GO TO 20
      TIM(1) = TF
      RETURN
C
   20 TSL = DLOG10(TS)
      TFL = DLOG10(TF)
      TLI = (TFL - TSL) / (NTIM - 1)
C
      DO 40 NN=1,NTIM
      TML = TSL + TLI * (NN - 1)
   40 TIM(NN) = 10.0D+0**TML
C
      RETURN
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE LABCRP (TAU,NTT,CREP,EMOD,SHRN,FPC,FPT)     
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)         
C	----------------------------------------------------------
C     INPUT LABORATORY TEST DATA FOR CREEP / SHRINKAGE    -
C	----------------------------------------------------------
      DIMENSION CREP(1)

	CALL FREBUF
	CALL FREECH
      CALL FREREL('E',EMOD,1)
      CALL FREREL('S',SHRN,1)
      CALL FREREL('F',FPC ,1)                                              
      CALL FREREL('R',FPT ,1)     
	                                         
      CALL FREBUF
      CALL FREREL('C',CREP,NTT)  !CREEP COEF.
C
      RETURN
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ACIINP (ACI,IPRN,NOT)     
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)                   
C	----------------------------------------------------------
C     INPUT COEFFICIENTS FOR ACI CREEP / SHRINKAGE MODEL  -
C	----------------------------------------------------------
	DIMENSION ACI(1)

      A	=  4.00D+0
      B	=  8.50D-1
      C	=  1.25D+0
      D	=  1.18D-1
      E	=  1.00D+0
      F	=  3.50D+1
      T	=  7.00D+0    !days
      W	=  150.0D0    !pcf
      R	=  0.0D0      !TENSILE STRENGTH
	H	=  40.0D0     !% 
	TH	=  6.0D0      !     
	SL  =  2.7D0      !
	ACI(1:12) = [A,B,C,D,E,F,T,W,R,H,TH,SL]
C
      CALL FREREL ('A',ACI(1),1)
      CALL FREREL ('B',ACI(2),1)
      CALL FREREL ('C',ACI(3),1)
      CALL FREREL ('D',ACI(4),1)
      CALL FREREL ('E',ACI(5),1)
      CALL FREREL ('F',ACI(6),1)
      CALL FREREL ('J',ACI(7),1)
	CALL FREBUF
	CALL FREECH
C	CALL FREREL ('W',ACI(8),1)  USE WDN INSTEAD
C	CALL FREREL ('X',ACI(9),1) 
      CALL FREREL ('H',ACI(10),1) 
      CALL FREREL ('S',ACI(11),1) 
      CALL FREREL ('L',ACI(12),1) 

      IF(IPRN.GT.0) WRITE(NOT,1000) ACI(1:7),ACI(10:12)                     

      RETURN

 1000 FORMAT(/10X,'AGE FACT. A =',E15.6/,
	1        10X,'AGE FACT. B =',E15.6/,
	1        10X,'CRP FACT. C =',E15.6/,
	1        10X,'CRP FACT. D =',E15.6/,
	1        10X,'SRK FACT. E =',E15.6/,
	1        10X,'SRK FACT. F =',E15.6/,
	1        10X,'SRK TIME  J =',E15.6,' Days'/,
C	1        10X,'RHO       W =',E15.6,' pcf'/,
C	1        10X,'TEN STRNG.X =',E15.6,' sqrt(lb)/in'/,
	1        10X,'HUMID     H =',E15.6,' %'/,
	1        10X,'NORM SIZE S =',E15.6/,
	1        10X,'SLUMP     L =',E15.6/)    

      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ACICRP (TAU,CTIM,NTT,ACI,FC28,IUF,IUL,          
     1                   CREP,EMOD,SHRN,FPC,FPT,CU,SU,
     2				   IUM,WDN,R28,MM)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------
C     GENERATE CREEP AND SHRINKAGE VALUES FOR ACI MODEL    -
C	-----------------------------------------------------------
	DIMENSION ACI(1)
      DIMENSION CTIM(1),CREP(1)

      A	=  ACI(1)
      B	=  ACI(2)
      C	=  ACI(3)
      D	=  ACI(4)
      E	=  ACI(5)
      F	=  ACI(6)
      T	=  ACI(7)
	W	=  ACI(8)       !pcf    
C	R28	=  ACI(9)       !TENSILE STRENGTH 
	H	=  ACI(10)      !% 
	TH	=  ACI(11)          
	SL  =  ACI(12)       


	FAL  = UNITFAC('LENG',IUL,4)
	TH   = TH*FAL
	SL   = SL*FAL

	FAM  = UNITFAC('MASS',IUM,4)  !LB
	FAL  = UNITFAC('LENG',IUL,5)  !FT
	IF(WDN.LE.0.0D0) WDN = W  !USE DEFAULT FROM ACI(8)
	W    = WDN*FAM/(FAL*FAL*FAL)

C
      FPC  = TAU*FC28 / (A + B * TAU)


	FAF = UNITFAC('FORC',IUF,4)
	FAL = UNITFAC('LENG',IUL,4)
	FC_psi = FPC*FAF/(FAL*FAL)


      FPT  = TAU*R28 / (A + B * TAU)                                             
      EMOD = 33.0D0 * W**1.5D+0 * DSQRT(FC_psi)    !W=pcf

      IF(MM.EQ.0) THEN
      FPCC  = FC28
	FCC_psi = FPCC*FAF/(FAL*FAL)                     
      EMOD = 33.0D0 * W**1.5D+0 * DSQRT(FCC_psi)    !W=pcf
      ENDIF

	FAF = UNITFAC('FORC',4,IUF)
	FAL = UNITFAC('LENG',4,IUL)
	EMOD= EMOD*FAF/(FAL*FAL)

      SHRN = 0.0D+0
      IF(TAU.GT.T) SHRN = SU*(TAU - T)**E / (F + (TAU - T)**E)
C	----------------------------------------------
C	CREEP CORRECTION FACTORS
C	----------------------------------------------
	RHK = 1.27-0.0067*H                              !HUMID H=%

	TK  = 1.0D0
	IF(TH.GT.6.0) TK  = 1.0D0-0.0167*(TH-6.0D0)      !THICKNESS   TH=in

	SLM = 0.81 + 0.07*(SL)                           !SLUMP  SLUMP=in

	ALK = C/TAU**D                                   !LOADING AGE  TAU=days
C	----------------------------------------------
	CRF = RHK*TK*ALK*SLM/EMOD   
C	----------------------------------------------

      DO 50 I=1,NTT
	TT = (CTIM(I)-TAU)**0.6D+0
	CREP(I) = CU*CRF*TT / (10.0D+0 + TT)   !C(Ti,TAU)
50	CONTINUE


      RETURN
      END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================               
      SUBROUTINE CEBINP (CEB,IPRN,NOT)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------
C     INPUT COEFFICIENTS FOR CEB/FIP CREEP / SHRINKAGE MODEL    -
C	----------------------------------------------------------------
	DIMENSION CEB(1)


C-----INTERPRET INPUT DATA-----
      HUMD = 70.0D+0 !%
      AREA = 30.0D+0 !mm**2
      PRIM =  1.0D+0 !mm
      TCEB =  3.0D+0 !Days

C-----STRENGTH COEFFICIENT-----
	S   = 0.25

C-----SHRINKAGE COEFFICIENTS-----
	BSC = 5.0

	CEB(1:6) = [HUMD,AREA,PRIM,TCEB,S,BSC]
C
      CALL FREREL('H',CEB(1),1)
      CALL FREREL('A',CEB(2),1)
      CALL FREREL('P',CEB(3),1)
      CALL FREREL('J',CEB(4),1)
      CALL FREREL('S',CEB(5),1)    !CEMENT TYPE FACTOR
	CALL FREREL('C',CEB(6),1)    !CEMENT TYPE FACTOR

      IF(IPRN.GT.0) WRITE(NOT,1000) CEB(1:6)

      RETURN

C
 1000 FORMAT(/10X,'HUMID      H =',E15.6,' %'/,
	1        10X,'AREA       A =',E15.6/,
	1        10X,'PERIM      P =',E15.6/,
	1        10X,'SRK TIME   J =',E15.6,' Days'/,
	1        10X,'TYPE COEF. S =',E15.6/,
	1        10X,'SRNK COEF. C =',E15.6/)  

      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
              
      SUBROUTINE KORINP (COR,IPRN,NOT)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------
C     INPUT COEFFICIENTS FOR KOREAN CREEP / SHRINKAGE MODEL    -
C	----------------------------------------------------------------
	DIMENSION COR(1)


C-----INTERPRET INPUT DATA-----
      HUMD = 70.0D+0 !%
      AREA = 30.0D+0 !mm**2
      PRIM =  1.0D+0 !mm
      TCEB =  3.0D+0 !Days

C-----STRENGTH COEFFICIENT-----
	S   = 0.25

C-----SHRINKAGE COEFFICIENTS-----
	BSC = 5.0

	COR(1:6) = [HUMD,AREA,PRIM,TCEB,S,BSC]
C
      CALL FREREL('H',COR(1),1)
      CALL FREREL('A',COR(2),1)
      CALL FREREL('P',COR(3),1)
      CALL FREREL('J',COR(4),1)
      CALL FREREL('S',COR(5),1)    !CEMENT TYPE FACTOR
	CALL FREREL('C',COR(6),1)    !CEMENT TYPE FACTOR

      IF(IPRN.GT.0) WRITE(NOT,1000) COR(1:6)

      RETURN

C
 1000 FORMAT(/10X,'HUMID      H =',E15.6,' %'/,
	1        10X,'AREA       A =',E15.6/,
	1        10X,'PERIM      P =',E15.6/,
	1        10X,'SRK TIME   J =',E15.6,' Days'/,
	1        10X,'TYPE COEF. S =',E15.6/,
	1        10X,'SRNK COEF. C =',E15.6/)  

      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	======================================================================= 
      SUBROUTINE CHNINP (CHN,IPRN,NOT)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------
C     INPUT COEFFICIENTS FOR CHINA JTG D62-2004 CREEP / SHRINKAGE MODEL    -
C	----------------------------------------------------------------
	DIMENSION CHN(1)


C-----INTERPRET INPUT DATA-----
      HUMD = 70.0D+0 !%
      AREA = 30.0D+0 !mm**2
      PRIM =  1.0D+0 !mm
      TCEB =  3.0D+0 !Days

C-----STRENGTH COEFFICIENT-----
	S   = 0.25

C-----SHRINKAGE COEFFICIENTS-----
	BSC = 5.0

	CHN(1:6) = [HUMD,AREA,PRIM,TCEB,S,BSC]
C
      CALL FREREL('H',CHN(1),1)
      CALL FREREL('A',CHN(2),1)
      CALL FREREL('P',CHN(3),1)
      CALL FREREL('J',CHN(4),1)
      CALL FREREL('S',CHN(5),1)    !CEMENT TYPE FACTOR
	CALL FREREL('C',CHN(6),1)    !CEMENT TYPE FACTOR

      IF(IPRN.GT.0) WRITE(NOT,1000) CHN(1:6)

      RETURN

C
 1000 FORMAT(/10X,'HUMID      H =',E15.6,' %'/,
	1        10X,'AREA       A =',E15.6/,
	1        10X,'PERIM      P =',E15.6/,
	1        10X,'SRK TIME   J =',E15.6,' Days'/,
	1        10X,'TYPE COEF. S =',E15.6/,
	1        10X,'SRNK COEF. C =',E15.6/)  

      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
                         
      SUBROUTINE CEBCRP (TAU,CTIM,NTT,CEB,FC28,IUF,IUL,
	1                   CREP,EMOD,SHRN,FPC,FPT,CU,SU,MM)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------
C     GENERATE CREEP VALUES ACCORDING TO CEB-FIP MODEL    -
C	----------------------------------------------------------
	DIMENSION CEB(1)
      DIMENSION CTIM(NTT),CREP(NTT,2)

	HUMD = CEB(1)
	AREA = CEB(2)
	PRIM = CEB(3)
	TCEB = CEB(4)
	S    = CEB(5)
	BSC  = CEB(6)

C	NOTATIONAL SIZE
	FAL  = UNITFAC('LENG',IUL,2)
	AREA = AREA*(FAL*FAL)
	PRIM = PRIM*FAL
	THIK = (AREA+AREA)/PRIM !mm

C-----COMPRESSIVE STRENGTH-----
	IF (TAU.GT.0.0D0) THEN
	TFNC = S*(1.0-SQRT(28.0D0/TAU))
	ELSE
	TFNC = 0.0D0
	ENDIF

	DFPC= 8.0  !Mpa
	FAF = UNITFAC('FORC',6,IUF)
	FAL = UNITFAC('LENG',0,IUL)
	DF_Mpa = DFPC*FAF/(FAL*FAL)
	FCM = FC28 + DF_Mpa

	BCC = DEXP(TFNC)
	FPC = BCC*FCM

C-----TENSILE STRENGTH-----
	Fct  = 1.40 !Mpa mean
	Fcko = 10.0 !Mpa
	FAF  = UNITFAC('FORC',6,IUF)
	FAL  = UNITFAC('LENG',0,IUL)
	Fct  = Fct*FAF/(FAL*FAL)
	Fcko = Fcko*FAF/(FAL*FAL)
	FPT  = Fct*(FPC/Fcko)**(2.0/3.0)

C-----ELASTIC MODULUS-----
	Fcmo = 10.0   !Mpa
	Eco  = 2.15E4 !Mpa
	FAF  = UNITFAC('FORC',6,IUF)
	FAL  = UNITFAC('LENG',0,IUL)
	Fcmo = Fcmo*FAF/(FAL*FAL)
	Eco  = Eco*FAF/(FAL*FAL)
	EMOD = SQRT(BCC)*Eco*(FCM/Fcmo)**(1.0/3.0)

      IF(MM.EQ.0) EMOD = Eco*(FCM/Fcmo)**(1.0/3.0)

C-----SHRINKAGE FACTOR-----
	SHRN = 0.0D+0
      IF(TAU.GT.TCEB)   THEN
	BRH = 0.25
	IF (HUMD.LT.99.0) BRH = -1.55*(1.0-(HUMD/100.0)**3.0)
	EPSH = 160.0+10.0*BSC*(9.0-FCM/Fcmo)
	EPSH = EPSH*1.0E-6
	BS = ( (TAU-TCEB)/(350.0*(THIK/100.0)**2.0+(TAU-TCEB)) )**0.5  !THIK = mm
	SHRN = -BRH*BS*EPSH*SU                                            !MUST BE MINUS SIGN DRYING
	ENDIF

C
C-----CREEP COEFFICIENTS-----

C	HUMIDITY AND THICKNESS FACTOR/EXPOSURE FACTOR
	BH = 150.0*( 1.0+ (1.20*HUMD/100.0)**18.0 )*THIK/100.0+250.0
	IF (BH.GT.1500.0) BH = 1500.0

C	Humidity and Thickness Factor/Exposure Factor
	PRH = 1.0+ (1.0-HUMD/100.0) /( 0.46* (THIK/100.0)**(1.0/3.0) )

	BFCM = 5.30/SQRT(FCM/Fcmo)

      DO 10 I=1,NTT
		TTAU = CTIM(I)-TAU
		BCTT0 = (TTAU/(BH+TTAU))**0.3
		BTOT = 1.0/(0.10+TAU**0.20)
		CREP(I,1) =  BCTT0*PRH*BFCM*BTOT*CU/EMOD  !Total
		CREP(I,2) =  BCTT0*PRH*BFCM*BTOT*CU/EMOD  !Recovery  Same with total (new code)
10	CONTINUE
C
      RETURN
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      FUNCTION CEBBFT (TAU,THIC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------
C     EVALUATE FLOW FUNCTION FOR CEB-FIP MODEL    -
C	--------------------------------------------------
      THIK = THIC / 10.0D+0
      F1 = 1.0 - 0.12 * THIK**0.3
      F2 = 1.0D+11 / (TAU**3 + 1.0D+11)
      TTF1 = TAU**F1
C
      BB = 1.0D+0 + THIK**1.7D+0 / 26600.0D+0
      BB = BB * TTF1 / (TTF1 + 21.5D+0)
      BB = BB + 4.0D+0 * TAU**3 / (TAU**5 + 2430.0D+0)
      BB = BB**F2
C
      CEBBFT = BB
C
      RETURN
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
                  
      SUBROUTINE KORCRP (TAU,CTIM,NTT,COR,FC28,IUF,IUL,
	1                   CREP,EMOD,SHRN,FPC,FPT,CU,SU,MM)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------
C     GENERATE CREEP VALUES ACCORDING TO KOREAN MODEL    -
C	----------------------------------------------------------
	DIMENSION COR(1)
      DIMENSION CTIM(NTT),CREP(NTT,2)

	HUMD = COR(1)
	AREA = COR(2)
	PRIM = COR(3)
	TCOR = COR(4)
	S    = COR(5)
	BSC  = COR(6)

C	NOTATIONAL SIZE
	FAL  = UNITFAC('LENG',IUL,2)
	AREA = AREA*(FAL*FAL)
	PRIM = PRIM*FAL
	THIK = (AREA+AREA)/PRIM !mm

C-----COMPRESSIVE STRENGTH-----
	IF (TAU.GT.0.0D0) THEN
	TFNC = S*(1.0-SQRT(28.0D0/TAU))
	ELSE
	TFNC = 0.0D0
	ENDIF

	DFPC= 8.0  !Mpa
	FAF = UNITFAC('FORC',6,IUF)
	FAL = UNITFAC('LENG',0,IUL)
	DF_Mpa = DFPC*FAF/(FAL*FAL)
	FCM = FC28 + DF_Mpa

	BCC = DEXP(TFNC)
	FPC = BCC*FCM

C-----TENSILE STRENGTH-----
	Fct  = 1.40 !Mpa mean
	Fcko = 10.0 !Mpa
	FAF  = UNITFAC('FORC',6,IUF)
	FAL  = UNITFAC('LENG',0,IUL)
	Fct  = Fct*FAF/(FAL*FAL)
	Fcko = Fcko*FAF/(FAL*FAL)
	FPT  = Fct*(FPC/Fcko)**(2.0/3.0)

C-----ELASTIC MODULUS-----
	Fcmo = 10.0   !Mpa
	Eco  = 2.1545E4 !Mpa
	FAF  = UNITFAC('FORC',6,IUF)
	FAL  = UNITFAC('LENG',0,IUL)
	Fcmo = Fcmo*FAF/(FAL*FAL)
	Eco  = Eco*FAF/(FAL*FAL)
	EMOD = SQRT(BCC)*Eco*(FCM/Fcmo)**(1.0/3.0)

      IF(MM.EQ.0) EMOD = Eco*(FCM/Fcmo)**(1.0/3.0)

C-----SHRINKAGE FACTOR-----
	SHRN = 0.0D+0
      IF(TAU.GT.TCOR)   THEN
	BRH = 0.25
	IF (HUMD.LT.99.0) BRH = -1.55*(1.0-(HUMD/100.0)**3.0)
	EPSH = 160.0+10.0*BSC*(9.0-FCM/Fcmo)
	EPSH = EPSH*1.0E-6
	BS = ( (TAU-TCOR)/(350.0*(THIK/100.0)**2.0+(TAU-TCOR)) )**0.5  !THIK = mm
	SHRN = -BRH*BS*EPSH*SU                                            !MUST BE MINUS SIGN DRYING
	ENDIF

C
C-----CREEP COEFFICIENTS-----

C	HUMIDITY AND THICKNESS FACTOR/EXPOSURE FACTOR
	BH = 150.0*( 1.0+ (1.20*HUMD/100.0)**18.0 )*THIK/100.0+250.0
	IF (BH.GT.1500.0) BH = 1500.0

C	Humidity and Thickness Factor/Exposure Factor
	PRH = 1.0+ (1.0-HUMD/100.0) /( 0.46* (THIK/100.0)**(1.0/3.0) )

	BFCM = 5.3125/SQRT(FCM/Fcmo)

      DO 10 I=1,NTT
		TTAU = CTIM(I)-TAU
		BCTT0 = (TTAU/(BH+TTAU))**0.3
		BTOT = 1.0/(0.10+TAU**0.20)
		CREP(I,1) =  BCTT0*PRH*BFCM*BTOT*CU/EMOD  !Total
		CREP(I,2) =  BCTT0*PRH*BFCM*BTOT*CU/EMOD  !Recovery  Same with total (new code)
10	CONTINUE
C
      RETURN
      END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CHNCRP (TAU,CTIM,NTT,CHN,FC28,IUF,IUL,
	1                   CREP,EMOD,SHRN,FPC,FPT,CU,SU,MM)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------
C     GENERATE CREEP VALUES ACCORDING TO CHINA JTG D62-2004 MODEL    -
C	----------------------------------------------------------
	DIMENSION CHN(1)
      DIMENSION CTIM(NTT),CREP(NTT,2)

	HUMD = CHN(1)
	AREA = CHN(2)
	PRIM = CHN(3)
	TCOR = CHN(4)
	S    = CHN(5)
	BSC  = CHN(6)

C	NOTATIONAL SIZE
	FAL  = UNITFAC('LENG',IUL,2)
	AREA = AREA*(FAL*FAL)
	PRIM = PRIM*FAL
	THIK = (AREA+AREA)/PRIM !mm

C-----COMPRESSIVE STRENGTH-----
	IF (TAU.GT.0.0D0) THEN
	TFNC = S*(1.0-SQRT(28.0D0/TAU))
	ELSE
	TFNC = 0.0D0
	ENDIF

	DFPC= 8.0  !Mpa
	FAF = UNITFAC('FORC',6,IUF)
	FAL = UNITFAC('LENG',0,IUL)
	DF_Mpa = DFPC*FAF/(FAL*FAL)
	FCM = FC28 + DF_Mpa

	BCC = DEXP(TFNC)
	FPC = BCC*FCM

C-----TENSILE STRENGTH-----
	Fct  = 1.40 !Mpa mean
	Fcko = 10.0 !Mpa
	FAF  = UNITFAC('FORC',6,IUF)
	FAL  = UNITFAC('LENG',0,IUL)
	Fct  = Fct*FAF/(FAL*FAL)
	Fcko = Fcko*FAF/(FAL*FAL)
	FPT  = Fct*(FPC/Fcko)**(2.0/3.0)

C-----ELASTIC MODULUS-----
	Fcmo = 10.0   !Mpa
	Eco  = 2.1545E4 !Mpa
	FAF  = UNITFAC('FORC',6,IUF)
	FAL  = UNITFAC('LENG',0,IUL)
	Fcmo = Fcmo*FAF/(FAL*FAL)
	Eco  = Eco*FAF/(FAL*FAL)
	EMOD = SQRT(BCC)*Eco*(FCM/Fcmo)**(1.0/3.0)

      IF(MM.EQ.0) EMOD = Eco*(FCM/Fcmo)**(1.0/3.0)

C-----SHRINKAGE FACTOR-----
	SHRN = 0.0D+0
      IF(TAU.GT.TCOR)   THEN
	BRH = 0.25
	IF (HUMD.LT.99.0) BRH = -1.55*(1.0-(HUMD/100.0)**3.0)
	EPSH = 160.0+10.0*BSC*(9.0-FCM/Fcmo)
	EPSH = EPSH*1.0E-6
	BS = ( (TAU-TCOR)/(350.0*(THIK/100.0)**2.0+(TAU-TCOR)) )**0.5  !THIK = mm
	SHRN = -BRH*BS*EPSH*SU                                     !MUST BE MINUS SIGN DRYING  
	ENDIF

C
C-----CREEP COEFFICIENTS-----

C	HUMIDITY AND THICKNESS FACTOR/EXPOSURE FACTOR
	BH = 150.0*( 1.0+ (1.20*HUMD/100.0)**18.0 )*THIK/100.0+250.0
	IF (BH.GT.1500.0) BH = 1500.0

C	Humidity and Thickness Factor/Exposure Factor
	PRH = 1.0+ (1.0-HUMD/100.0) /( 0.46* (THIK/100.0)**(1.0/3.0) )

	BFCM = 5.3125/SQRT(FCM/Fcmo)

      DO 10 I=1,NTT
		TTAU = CTIM(I)-TAU
		BCTT0 = (TTAU/(BH+TTAU))**0.3
		BTOT = 1.0/(0.10+TAU**0.20)
		CREP(I,1) =  0.9*BCTT0*PRH*BFCM*BTOT*CU/EMOD  !Total  --- 0.9 FROM CHINA CODE 
		CREP(I,2) =  0.9*BCTT0*PRH*BFCM*BTOT*CU/EMOD  !Recovery  Same with total (new code)  --- 0.9 FROM CHINA CODE 
10	CONTINUE
C
      RETURN
      END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CRPERR (TAU,ALSC,CCOF,CREP,CTIM,
     1                   EMOD,SHRN,FPC,FPT,
     1				   NRT,NTT,KL,NATN,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------
C     PRINT CREEP MODEL DIAGNOSIS 
C	----------------------------------------------
      DIMENSION ALSC(NTT,NRT),CCOF(NRT,KL),CREP(NTT,KL),CTIM(NTT),
     1          CRAT(NTT,KL),CERR(NTT,KL),CEST(NTT,KL)

C	----------------------------------------------
C	PRINT CREEP COEFFICIENTS
C	----------------------------------------------
      WRITE(NOT,1000) NATN,TAU,EMOD,SHRN,FPC,FPT
      DO 10 I=1,NRT
   10 WRITE(NOT,1001) I,(CCOF(I,L),L=1,KL)

C	----------------------------------------------
C	COMPARE CREEP DATA WITH LEAST SQUARES FIT
C	----------------------------------------------
      DO 20 L=1,KL
      DO 20 I=1,NTT

      EST = 0.0D+0
      DO 15 J=1,NRT
   15 EST = EST + ALSC(I,J) * CCOF(J,L)   !G*A

      CERR(I,L) = CREP(I,L) - EST
      CRAT(I,L) = 0.0D+0
      IF(CREP(I,L).EQ.0.0D+0) GO TO 20
      CRAT(I,L) = 100.0 * CERR(I,L) / CREP(I,L)
   20 CEST(I,L) = EST

C	----------------------------------------------
C	PRINT COMPARISON FOR EACH OBSERVATION TIME
C	----------------------------------------------
      WRITE(NOT,1200) NATN,TAU
      N2 = 0
      DO 50 I=1,NTT
      IF(N2.EQ.NTT) GO TO 60
      N1 = N2 + 1
      N2 = N2 + 6
      N2 = MIN0(N2,NTT)

                  WRITE(NOT,1201) (M,M=N1,N2)
                  WRITE(NOT,1202) (CTIM(M),M=N1,N2)
                  WRITE(NOT,1204) (CREP(M,1),M=N1,N2)
      IF(KL.EQ.2) WRITE(NOT,1205) (CREP(M,2),M=N1,N2)
                  WRITE(NOT,1206) (CEST(M,1),M=N1,N2)
      IF(KL.EQ.2) WRITE(NOT,1207) (CEST(M,2),M=N1,N2)
                  WRITE(NOT,1208) (CERR (M,1),M=N1,N2)
      IF(KL.EQ.2) WRITE(NOT,1209) (CERR (M,2),M=N1,N2)
                  WRITE(NOT,1210) (CRAT (M,1),M=N1,N2)
      IF(KL.EQ.2) WRITE(NOT,1211) (CRAT (M,2),M=N1,N2)

   50 CONTINUE
C	----------------------------------------------

   60 RETURN
C
 1000 FORMAT(/,' CONCRETE PARAMETERS FOR LOADING AGE',I5,', DAY =',
     1 F11.2,//,
     1 4X,'ELASTIC MODULUS  =',1PE11.3,/,
     1 4X,'SHRINKAGE STRAIN =',1PE11.3,/,
     1 4X,'COMP. STRENGTH   =',1PE11.3,/,
     1 4X,'TENS. STRENGTH   =',1PE11.3,//,
     1 4X,'CREEP COMPLIANCE COEFFICIENTS    I    A(I)')

 1001 FORMAT(33X,I5,1P2E13.6)
 1200 FORMAT(/,' CREEP MODEL DIAGNOSTICS FOR LOADING AGE',I5,', DAY =',
     1 F11.2)
 1201 FORMAT(/,'  TIME NUMBER',I6,5I11)
 1202 FORMAT(  '  TIME (DAYS)',6F11.4)
 1204 FORMAT(  '  CREEP VALUE',1P6E11.3)
 1205 FORMAT(  '   (RECOVERY)',1P6E11.3)
 1206 FORMAT(  '  CREEP ESTIM',1P6E11.3)
 1207 FORMAT(  '   (RECOVERY)',1P6E11.3)
 1208 FORMAT(  '  CREEP ERROR',1P6E11.3)
 1209 FORMAT(  '   (RECOVERY)',1P6E11.3)
 1210 FORMAT(  '  PERCENT ERR',6F11.4)
 1211 FORMAT(  '   (RECOVERY)',6F11.4)
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CRPPRN (ICPT,PRAM,NA,NCPC,NCPT,IPRN,NOT)       
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------
C     LIST ALL COEFFICIENTS FOR ALL PARAMETER TYPES    -
C	-------------------------------------------------------
      CHARACTER*1 TITLE(19,6)
      DIMENSION ICPT(7,NCPT),PRAM(NA,NCPC,NCPT)
      DATA TITLE
     1           /'E','X','P','E','R','I','M','E','N','T','A','L',' ',
     2    'D','A','T','A',' ',' ',
     3            'A','C','I',' ',
     4    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S',
     5            'C','E','B',' ',
     6    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S',
     7            'L','I','E','A','R',' ','A','N','A','L','Y','S','I',
     8    'S',' ',' ',' ',' ',' ',
     9            'K','O','R',' ',
     1    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S',
     2            'J','T','G',' ',
     3    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S'/
C	----------------------------------------------
C	FOR EACH PARAMETER TYPE 
C	----------------------------------------------
      DO 60 NPIN=1,NCPT
C
      KM  = ICPT(1,NPIN)
      NAT = ICPT(2,NPIN)
      NTT = ICPT(3,NPIN)
      NRT = ICPT(4,NPIN)
      KL  = ICPT(5,NPIN)
C
      IF(IPRN.LT.1) GO TO 60
C
      WRITE(NOT,1000) NPIN,(TITLE(K,KM),K=1,19),NAT,NRT
C	----------------------------------------------
C	PRINT PARAMETERS AT EACH LOADING AGE
C	----------------------------------------------
      N2 = 0
      DO 50 I=1,NAT
      IF(N2.EQ.NAT) GO TO 60
      N1 = N2 + 1
      N2 = N2 + 6
      N2 = MIN0(N2,NAT)
C
      WRITE(NOT,1201) (M,M=N1,N2)
      WRITE(NOT,1202) (PRAM(M,1,NPIN),M=N1,N2)
      WRITE(NOT,1204) (PRAM(M,2,NPIN),M=N1,N2)
      WRITE(NOT,1205) (PRAM(M,3,NPIN),M=N1,N2)
      WRITE(NOT,1207) (PRAM(M,4,NPIN),M=N1,N2)
      WRITE(NOT,1208) (PRAM(M,5,NPIN),M=N1,N2)                            
      II = 5                                                               
      DO 40 K=1,KL
      DO 40 J=1,NRT
      II = II + 1
   40 WRITE(NOT,1206) J,K,(PRAM(M,II,NPIN),M=N1,N2)

   50 CONTINUE

C	----------------------------------------------
   60 CONTINUE
C	----------------------------------------------


      RETURN
C
 1000 FORMAT(//,' CONCRETE MODEL',I5,' BASED ON ',19A1,/,1X,48('-'),//,
     1       5X,' NUMBER OF LOADING AGES      =',I5,/,
     2       5X,' NUMBER OF RETARDATION TIMES =',I5)
 1201 FORMAT(/,'  AGE NUMBER ',I6,5I11)
 1202 FORMAT(  '  AGE (DAYS) ',6F11.4)
 1204 FORMAT(  '  EMOD       ',1P6E11.3)
 1205 FORMAT(  '  SHRN       ',1P6E11.3)
 1206 FORMAT(  '  COEF A(',I1,',',I1,')',1P6E11.3)
 1207 FORMAT(  '  FPC        ',1p6e11.3)                                   
 1208 FORMAT(  '  FPT        ',1p6e11.3)                                   
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE LSTSQR (ALSC,BRHS,XLSS,ATA,NRT,NTT,KL,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------
C     SIMPLE LEAST SQUARES SOLUTION   -
C	--------------------------------------
      DIMENSION ALSC(NTT,NRT),BRHS(NTT,KL),XLSS(NRT,KL),ATA(NRT,NRT)

C
C-----FORM THE LEAST SQUARES COEFFICIENT MATRIX-----
C
      DO 50 I=1,NRT
      DO 50 J=I,NRT
      BB = 0.0D+0
      DO 40 II=1,NTT
   40 BB = BB + ALSC(II,I) * ALSC(II,J)
      ATA(I,J) = BB
   50 ATA(J,I) = BB
C
C-----FORM THE RIGHT HAND SIDE VECTOR-----
C
      DO 70 L=1,KL
      DO 70 I=1,NRT
      BB = 0.0D+0
      DO 60 J=1,NTT
   60 BB = BB + ALSC(J,I) * BRHS(J,L)
   70 XLSS(I,L) = BB
C
C-----SOLVE SYMMETRICAL NORMAL EQUATIONS-----
C
      CALL SYMSOL (ATA,XLSS,NRT,0,KL,NOT)
C
      RETURN
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SYMSOL (A,B,NN,KEY,NL,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	------------------------------------------------
C     SOLVE SYMMETRIC SIMULTANEOUS EQUATIONS    -
C	------------------------------------------------
      DIMENSION A(NN,NN),B(NN,NL)

C
C-----TRIANGULAR FACTORIZATION-----
C
      NM = NN - 1
      IF(KEY.EQ.2) GO TO 50
C
      DO 40 N=1,NM
      D = A(N,N)
      IF(D.EQ.0.0) WRITE(NOT,2000) N
      N1 = N + 1
      DO 30 J=N1,NN
      IF(A(N,J).EQ.0.0D+0) GO TO 30
      AA = A(N,J) / D
      A(N,J) = AA
      DO 20 I=J,NN
      A(I,J) = A(I,J) - A(I,N) * AA
   20 A(J,I) = A(I,J)
   30 CONTINUE
   40 CONTINUE
C
C-----FORWARD REDUCTION-----
C
   50 IF(KEY.EQ.1) RETURN
      DO 90 L=1,NL
      DO 70 N=1,NM
      B(N,L) = B(N,L) / A(N,N)
      N1 = N + 1
      DO 70 I=N1,NN
   70 B( I,L) = B( I,L) - A(I,N) * B(N,L)
      B(NN,L) = B(NN,L) / A(NN,NN)
C
C-----BACK SUBSTITUTION-----
C
      N  = NN
      DO 85 I=1,NM
      N1 = N
      N  = N - 1
      DO 85 J=N1,NN
   85 B(N,L) = B(N,L) - A(N,J) * B(J,L)
C
   90 CONTINUE
      RETURN
C
 2000 FORMAT (/,' WARNING: ZERO DIAGONAL IN EQUATION NUMBER',I5)
      END
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CNCTYP (FCNC,ICNC,ICPT,GAM0,NPTP,NPIN,IPRN,NI,NOT)
C	----------------------------------------
C     PARSE PARAMETER TYPE DATA LINE    -
C	----------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 TITLE(19,6)
      LOGICAL TEST
      DIMENSION FCNC(20,1),ICNC(5,1),ICPT(7,1),GAM0(1),NI(3),NIT(3)
      DATA TITLE 
     1           /'E','X','P','E','R','I','M','E','N','T','A','L',' ',
     2    'D','A','T','A',' ',' ',
     3            'A','C','I',' ',
     4    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S',
     5            'C','E','B',' ',
     6    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S',
     7            'L','I','E','A','R',' ','A','N','A','L','Y','S','I',
     8    'S',' ',' ',' ',' ',' ',
     9            'K','O','R',' ',
     1    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S',
     2            'J','T','G',' ',
     3    'R','E','C','O','M','M','E','N','D','A','T','I','O','N','S'/

C	-------------------------------------------------------------
C	START CONCRETE PROPERTIES INPUT
C	-------------------------------------------------------------
      CALL FREBUF
	CALL FREECH
      CALL FREINT(' ',NPIN,1)  !READ CONCRETE NO.


	M = 0
      CALL FREINT('M',M,1)     !0 = Linear, 1 = Concrete, 2 = Concrete with nl stress-strain
	ICNC(1,NPIN) = M

	IF(IPRN.GT.0) THEN
	IF(M.EQ.0) WRITE(NOT,501) NPIN,M 
	IF(M.EQ.1) WRITE(NOT,502) NPIN,M 
	IF(M.EQ.2) WRITE(NOT,503) NPIN,M 
	ENDIF

      FCNC( 2,NPIN) = 1.0D0       !Cu
      FCNC( 3,NPIN) = 3.1E-3		!Su
      FCNC( 7,NPIN) = 1.0D0       !B
      FCNC( 8,NPIN) = 3.1E-3		!Eu
      FCNC( 9,NPIN) = 3.0E-3      !Ee
      FCNC(10,NPIN) = 0.85        !Be
      FCNC(11,NPIN) = 1.0D0       !r1      
      FCNC(12,NPIN) = 1.0D0       !r2

      CALL FREREL('W',FCNC(4,NPIN),1)  !WEIGHT DENSITY
      CALL FREREL('A',FCNC(5,NPIN),1)  !ALFA THERMAL

	CALL FREREL('F',FCNC(1,NPIN),1)  !Fcp

	FCNC(14:15,NPIN) = -1.0D0
	CALL FREREL('E',FCNC(14,NPIN),1) !YOUNG  *must read if KM=4* or M = 0
	CALL FREREL('X',FCNC(15,NPIN),1) !Fct    *must read if KM=4*
	CALL FREREL('V',FCNC(16,NPIN),1) !POISSON RATIO

      CALL FREREL('C',FCNC(2,NPIN),1)   !Cu
      CALL FREREL('S',FCNC(3,NPIN),1)   !Su
      CALL FREREL('O',FCNC(11,NPIN),1)  !r1
      CALL FREREL('P',FCNC(12,NPIN),1)  !r2


      IF(FCNC(11,NPIN).EQ.1.0D0) THEN
      C1 = 1.0D0
      C2 = 0.0D0
	GOTO 45
      ENDIF
      C1 = (FCNC(12,NPIN) - FCNC(11,NPIN)) / (1.0D0 - FCNC(11,NPIN))
      C2 =  FCNC(11,NPIN) * (1.0D0 - C1)
   45 FCNC(13,NPIN) = FCNC(11,NPIN)   !r1
      FCNC(11,NPIN) = C1
      FCNC(12,NPIN) = C2
	
	IF(IPRN.GT.0) THEN 
	WRITE(NOT,1001) 
	1	FCNC(4:5,NPIN),FCNC(1,NPIN),FCNC(2:3,NPIN),FCNC(11:12,NPIN)
	IF(FCNC(14,NPIN).NE.-1.0D0) WRITE(NOT,1002) FCNC(14,NPIN)
	IF(FCNC(15,NPIN).NE.-1.0D0) WRITE(NOT,1004) FCNC(15,NPIN)
	ENDIF

	IF(M.EQ.2) THEN
      CALL FREBUF
      CALL FREINT('C',FCNC( 7,NPIN),1)
      CALL FREINT('U',FCNC( 8,NPIN),1)
      CALL FREINT('E',FCNC( 9,NPIN),1)
      CALL FREINT('R',FCNC(10,NPIN),1)
      CALL FREINT('K',ICNC( 2,NPIN),1)   !TENSION STIFFENING FLAG
	IF(IPRN.GT.0) WRITE(NOT,1003) FCNC(7:10,NPIN),ICNC(2,NPIN)
	ENDIF


C	-------------------------------------------------------------
C	START CONCRETE PARAMETER INPUT (CREEP&SHRINKAGE)
C	-------------------------------------------------------------
      CALL FREBUF
	CALL FREECH

C-----FIND MODEL TYPE (KM), NUMBER OF LOADING TERMS (KL)-----
C	KM = 1  !LAB
C	KM = 2  !ACI
C	KM = 3  !CEB
C	KM = 4  !LINEAR
C	KM = 5  !KOREAN CODE
C	KM = 6  !CHINA JTG CODE
	KM = 2
      CALL FREINT('M',KM,1)

C	KL = NUMBER OF CREEP COEF SET.  2 FOR CEB (CREEP AND CREEP ROCOVERY)	
	KL = 1
	IF(KM.EQ.3) KL = 2 !FOR CEB WE NEED TWO TIME OF LEAST SQUARE SOLVE
	IF(KM.EQ.5) KL = 2 !FOR KOR WE NEED TWO TIME OF LEAST SQUARE SOLVE
	IF(KM.EQ.6) KL = 2 !FOR CHN WE NEED TWO TIME OF LEAST SQUARE SOLVE

C-----FIND NUMBER OF LDG AGES, CREEP TIMES, RETARDATION TIMES-----

	NA = NI(1)
	NT = NI(2)
	NR = NI(3)

	NIT(1:3) = [NA,NT,NR] !NA,NT,NR
      CALL FREINT('T',NIT,3)
	NAT = NIT(1)
	NTT = NIT(2)
	NRT = NIT(3)
	NAT = MIN0 (NAT,NA)
	NTT = MIN0 (NTT,NT)
	NRT = MIN0 (NRT,NR)
	
	ICRP = 1  !DEFAULT C-C
      CALL FREINT('L',ICRP,1)     !1=C-C 2=C-L 3=L-L  CREEP INTEGRATION FLAG


C-----SEE IF GENERATED OR INPUT AGES AND TIMES (IGEN)-----
      IGEN = 0
      CALL FREINT('G',IGEN,1)
      IF(KM.EQ.1) IGEN = 1  !LAB

C-----READ FIRST RETARDATION TIME-----
      GAM0(NPIN) = 5.0D+0
      CALL FREREL('R',GAM0(NPIN),1)

	FCNC(6,NPIN) = GAM0(NPIN)  !STORE THE FIRST RETARDATION TIME

C-----STORE MODEL SIZE AND TYPE
      ICPT(1,NPIN) = KM
      ICPT(2,NPIN) = NAT
      ICPT(3,NPIN) = NTT
      ICPT(4,NPIN) = NRT
      ICPT(5,NPIN) = KL
      ICPT(6,NPIN) = IGEN
      ICPT(7,NPIN) = ICRP

C-----PRINT IT OUT-----
      IF(IPRN.LT.2) RETURN
      IF(KM.NE.4) WRITE(NOT,2001) NPIN,NAT,NTT,NRT,GAM0(NPIN)
      WRITE(NOT,2002) (TITLE(K,KM),K=1,19)
      IF(KL.EQ.1) WRITE(NOT,2003)
      IF(KL.EQ.2) WRITE(NOT,2004)

C	-------------------------------------------------------------


      RETURN


501	FORMAT (/,10X,'    CONCRETE TYPE NO . . . . . .',I5,
	1		/,10X,'(M) CONCRETE STRESS-STRAIN TYPE.',I5,
     1        /,10X,'    LINEAR WITH SPECIFIC ELASTIC MODULUS'/)
502	FORMAT (/,10X,'    CONCRETE TYPE NO . . . . . .',I5,
	1		/,10X,'(M) CONCRETE STRESS-STRAIN TYPE.',I5,
     1        /,10X,'    LINEAR STRESS-STRAIN RELATION'/)
503	FORMAT (/,10X,'    CONCRETE TYPE NO . . . . . .',I5,
	1		/,10X,'(M) CONCRETE STRESS-STRAIN TYPE.',I5,
     1        /,10X,'    NONLINEAR STRESS-STRAIN RELATION'/)

1001	FORMAT (/,10X,'(W) MASS DENSITY . . . . . . . .',E12.5,
     1        /,10X,'(A) THERMAL EXPANSION COEF . . .',E12.5,
     1        /,10X,'(F) 28 DAYS COMP. STRENGTH . . .',E12.5,
     1        /,10X,'(C) ULTIMATE CREEP COEF. . . . .',E12.5,
     1        /,10X,'(S) ULTIMATE SHRINKAGE COEF. . .',E12.5,
     1        /,10X,'(C1)NONLINEAR CREEP FACTOR . . .',E12.5,
     1        /,10X,'(C2)NONLINEAR CREEP FACTOR . . .',E12.5/)

1002	FORMAT (/,10X,'(E) 28 DAYS YOUNG MODULUS. . . .',E12.5/)

1003	FORMAT (/,10X,'(C) Beta . . . . . . . . . . . .',E12.5,
     1        /,10X,'(U) Eu . . . . . . . . . . . . .',E12.5,
     1        /,10X,'(E) Ee . . . . . . . . . . . . .',E12.5,
     1        /,10X,'(R) Beta-e . . . . . . . . . . .',E12.5,
     1        /,10X,'(K) TENSION STIFFENING FLAG. . .',I5/)

1004	FORMAT (/,10X,'(X) TENSILE STRENGTH . . . . . .',E12.5/)

 2001 FORMAT (/,10X,'CONCRETE TYPE NUMBER  . . . .',I5,
     1        /,10X,'NUMBER OF LOADING AGES  . . .',I5,
     2        /,10X,'NUMBER OF CREEP TIMES . . . .',I5,
     3        /,10X,'NUMBER OF RETARDATION TIMES .',I5,
     4        /,10X,'SMALLEST  RETARDATION TIME  .',F9.3,' DAYS')
 2002 FORMAT (/,10X,'MODEL BASED ON ',19A1)
 2003 FORMAT (16X,'WITH TOTAL CREEP USED FOR BOTH CREEP ',
     1            'AND CREEP RECOVERY')
 2004 FORMAT (16X,'WITH DIFFERENT FUNCTIONS USED FOR CREEP ',
     1            'AND CREEP RECOVERY')


      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE TMPCNC(AGEC,MATN,EMOD,SHRN,FPC,FPT,RHO,POISN,GG,CCOF) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)	

      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

	DIMENSION CCOF(1)  !(NRT,KL)



C	MAXIMUM NUMBER OF LOADING TIME
	CALL INTFILL('%XNI',NA  ,1,1,0)
	CALL INTFILL('%XNI',MCPC,1,4,0)
	CALL  LOCATN('%PRM',KPRM,MPRM,NCPT,2)
	CALL  LOCATN('%FNC',KFNC,MFNC,NCPT,2)

	CALL INTFILL('%CPT',KM ,1,MATN,0)
	CALL INTFILL('%CPT',NAT,2,MATN,0)
	CALL INTFILL('%CPT',NRT,4,MATN,0)
	CALL INTFILL('%CPT',KL ,5,MATN,0)
	CALL INCCNC(AGEC,RCONDT(KPRM),RCONDT(KFNC),MCPC,MATN,
     1            KL,NA,NAT,NRT,KM,EMOD,SHRN,FPC,FPT,CCOF) 
	CALL RELFILL('%FNC',RHO  ,4 ,MATN,0)
	CALL RELFILL('%FNC',POISN,16,MATN,0)
	GG = EMOD/(2.0 + 2.0*POISN)

	RETURN
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE INCCNC(TAU,PRAM,FCNC,NCPC,NPIN,KL,NA,NAT,NRT,KM,
     1                  EMOD,SHRN,FPC,FPT,CCOF) 
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	------------------------------------------------------
C     GET CONCRETE PARAMETERS FROM TABLE FOR EACH FIBER    -
C	------------------------------------------------------

      DIMENSION PRAM(NA,NCPC,1),FCNC(20,1),CCOF(NRT,KL)
	ALLOCATABLE T(:)

	NC = 5 + NRT * KL
	ALLOCATE (T(NC))
	T(1:NC) = 0.0D0

	IF(KM.EQ.4) GOTO 55  !NO TIME DEPENDENT VARIATION
C
C-----FIND INTERPOLATION FACTOR-----
C
      RAT = 0.0D+0
      NAC = NAT - 1
      DO 20 N=1,NAC
      IF(TAU.LT.PRAM(N,1,NPIN).OR.TAU.GT.PRAM(N+1,1,NPIN)) GO TO 20
      N1 = N
      N2 = N + 1
      RAT=(TAU - PRAM(N1,1,NPIN)) / (PRAM(N2,1,NPIN) - PRAM(N1,1,NPIN))  !RAT=(T-T1)/(T2-T1)
      GO TO 40
20	CONTINUE
C
      IF(TAU.GT.PRAM(NAT,1,NPIN)) GO TO 30
      N1 = 1
      N2 = 1
      GO TO 40
30	N1 = NAT
      N2 = NAT
C
C-----INTERPOLATE AND FACTOR FOR CURRENT PARAMETERS-----
C
40	CONTINUE
	DO 50 K=2,NC
50	T(K) = PRAM(N1,K,NPIN) + RAT * (PRAM(N2,K,NPIN) - PRAM(N1,K,NPIN))  !T = P1 + RAT*(P2-P1)

55	CONTINUE

	SELECT CASE (KM) 
	CASE (1,2,3,5,6) !LAB ACI CEB KOR CHN
	EMOD = T(2) 
	SHRN = T(3)
	FPC  = T(4) 
	FPT  = T(5) 
	CASE (4)     !LIN
	EMOD = FCNC(14,NPIN) 
	SHRN = 0.0D0 
	FPC  = FCNC(1 ,NPIN) 
	FPT  = FCNC(15,NPIN)
	END SELECT

      II = 5
      DO 60 L=1,KL
      DO 60 N=1,NRT
	II = II + 1
	CCOF(N,L) = T(II)
60	CONTINUE


	DEALLOCATE (T)

C
      RETURN
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE INPMSL(FMSL,NMSL,NS,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	--------------------------------------------------
C     INPUT MATERIAL PROPERTIES FOR MILD STEEL    -
C	--------------------------------------------------
      DIMENSION FMSL(NS,1)

C
C     FMSL(1,M) = EMOD1
C     FMSL(2,M) = EMOD2
C     FMSL(3,M) = FSY
C     FMSL(4,M) = ESU
C     FMSL(5,M) = THERM EXPN COEFF
C	FMSL(6,M) = UNIT WEIGHT
C	FMSL(7,M) = TYPE OF MATERIAL BEHAVIOUR 1-LINEAR 2-BILINEAR
C	FMSL(8,M) = POISN
C	FMSL(9,M) = MASS DENSITY
	
      WRITE(NOT,1000)

C	READ MATERIAL PROPERTIES-----

      DO 50 N=1,NMSL
C
      CALL FREBUF

      CALL FREINT('N',M,1)
      CALL FREREL('E',FMSL(1,M),1)
      CALL FREREL('T',FMSL(2,M),1)
      CALL FREREL('Y',FMSL(3,M),1)
      CALL FREREL('U',FMSL(4,M),1)
      CALL FREREL('A',FMSL(5,M),1)
	CALL FREREL('W',FMSL(6,M),1) 
	CALL FREREL('M',FMSL(7,M),1) 
	CALL FREREL('D',FMSL(9,M),1)  !MASS DENSITY 
	CALL FREREL('V',FMSL(8,M),1)  !POISN

	IF (FMSL(7,M).EQ.1.) THEN 
		FMSL(2,M) = FMSL(1,M)
		FMSL(3,M) = 0.0D0
		FMSL(4,M) = 0.0D0
	ENDIF

50	CONTINUE


C-----PRINT MATERIAL PROPERTIES-----
      WRITE (NOT,1100)

      WRITE (NOT,2000) (M,(FMSL(I,M),I=1,6),M=1,NMSL)

      RETURN

 1000 FORMAT (//,' MILD STEEL PROPERTIES INPUT',/,1X,27('='))
 1100 FORMAT (/,'  STEEL NO    EMOD1      EMOD2       FSY        ESU',
     1 '       ALPHA','      WEIGHT',/)
 2000 FORMAT (I10,1P6E11.3)
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE INPSCB(FSPP,FSSC,NMSL,NPSL,NS,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	----------------------------------------------------------
C     INPUT MATERIAL PROPERTIES FOR CABLE
C	----------------------------------------------------------


      DIMENSION FSPP(NS,1),FSSC(NPSL,3,NMSL)

C
C     FSPP(1,M) = .1 PERCENT OFFSET YEILD STRESS.
C     FSPP(2,M) = UNIT WEIGHT.
C     FSPP(3,M) = ALFA, COEFF. OF THERMAL EXPANSION.
C     FSPP(4,M) = RELAXATION COEFFICIENT.
C     FSPP(5,M) = NUM OF STRESS-STRAIN POINTS
C	FSPP(6,M) = POISN
C	FSPP(7,M) = MASS DENSITY
C

C     FSSC(I,M) = STRESS AT POINT I.
C     FSSC(I,M) = STRAIN AT POINT I.
C     FSSC(I,M) = MODULUS FROM POINT (I-1), I.


      WRITE(NOT,1000)


C	-----READ MATERIAL PROPERTIES-----

      DO 50 N=1,NMSL

      CALL FREBUF
      CALL FREINT('N',M,1)
      CALL FREREL('F',FSPP(1,M),1)
      CALL FREREL('W',FSPP(2,M),1)
      CALL FREREL('A',FSPP(3,M),1)
      CALL FREREL('R',FSPP(4,M),1)
      CALL FREREL('D',FSPP(7,M),1) !MASS DENSITY 
      CALL FREREL('V',FSPP(6,M),1) !POISN

	FSPP(5,M) = 1.0
      CALL FREREL('P',FSPP(5,M),1) 
	NP = INT(FSPP(5,M))

      CALL FREBUF
      CALL FREREL('S',FSSC(1,1,M),NP)
      CALL FREBUF
      CALL FREREL('E',FSSC(1,2,M),NP)

      FSSC(1,3,M) = FSSC(1,1,M)/FSSC(1,2,M)
      DO 40 I=2,NP
      FSSC(I,3,M) = (FSSC(I,1,M)-FSSC(I-1,1,M))/
	1			  (FSSC(I,2,M)-FSSC(I-1,2,M))
40    CONTINUE

50    CONTINUE

C	-----PRINT MATERIAL PROPERTIES-----

      WRITE (NOT,1100)
      WRITE (NOT,2000) (M,(FSPP(I,M),I=1,4),M=1,NMSL)

      DO 60 M=1,NMSL
      WRITE (NOT,3000) M
	NP = INT(FSPP(5,M))
   60 WRITE (NOT,3100) (I,FSSC(I,1,M),FSSC(I,2,M),FSSC(I,3,M),I=1,NP)


      RETURN

 1000 FORMAT (//,' STAY PROPERTIES INPUT        ',/,1X,22('='),/)
 1100 FORMAT (   '  TYPE      YEILD      WEIGHT       ALFA       RLXS ',
     1         /,
     2 '    NO.    STRESS     DENSITY ',/)


 2000 FORMAT (I7,1P4E11.3)
 3000 FORMAT ( /,' STRESS-STRAIN CURVE POINTS FOR CABLE STAY ',
     1         'STEEL NO.',I3/1X,56('-'),
     2         /,4X,'POINT',12X,'STRESS',6X,'STRAIN',5X,'MODULUS'/)
 3100 FORMAT (I9,6X,1P3E12.3)
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE INPSPS(FSPP,FSSC,NMSL,NPSL,NS,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------
C     INPUT MATERIAL PROPERTIES FOR PRESTRESS
C	----------------------------------------------------------
      DIMENSION FSPP(NS,1),FSSC(NPSL,3,NMSL)

C     FSPP(1,M) = .1 PERCENT OFFSET YEILD STRESS.
C     FSPP(2,M) = UNIT WEIGHT.
C     FSPP(3,M) = ALFA, COEFF. OF THERMAL EXPANSION.
C     FSPP(4,M) = RELAXATION COEFFICIENT.
C     FSPP(5,M) = FRICTION COEF.
C     FSPP(6,M) = WOBBLE   COEF.
C     FSPP(7,M) = STRESS-STRAIN POINTS
C     FSPP(8,M) = POISN
C     FSPP(9,M) = MASS DENSITY
C
C     FSPP(10,M) = CHINA CODE TENSIONING FACTOR = 1.0 FOR NORMAL-TENSIONED = 0.9 FOR OVER-TENSIONED
C     FSPP(11,M) = CHINA CODE RELAXATION FACTOR = 1.0 FOR NORMAL-RELAXATION STRAND = 0.3 FOR LOW-RELAXATION STRAND

C     FSSC(I,M) = STRESS AT POINT I.
C     FSSC(I,M) = STRAIN AT POINT I.
C     FSSC(I,M) = MODULUS FROM POINT (I-1), I.


      WRITE(NOT,1000)


C	-----READ MATERIAL PROPERTIES-----

      DO 50 N=1,NMSL

      CALL FREBUF
      CALL FREINT('N',M,1)
      CALL FREREL('F',FSPP(1,M),1)
      CALL FREREL('D',FSPP(2,M),1)  
      CALL FREREL('A',FSPP(3,M),1)
      CALL FREREL('R',FSPP(4,M),1)
      CALL FREREL('U',FSPP(5,M),1)
      CALL FREREL('W',FSPP(6,M),1)
      CALL FREREL('V',FSPP(8,M),1)  !POISN
      CALL FREREL('H',FSPP(9,M),1)  !MASS DENSITY  
      
      FSPP(10,M) = 0.0D0
      FSPP(11,M) = 0.0D0
      CALL FREREL('M',FSPP(10,M),1) !CHINA CODE TENSIONING FACTOR
      CALL FREREL('L',FSPP(11,M),1) !CHINA CODE RELAXATION FACTOR

	FSPP(7,M) = 1.0
      CALL FREREL('P',FSPP(7,M),1) 
	NP = INT(FSPP(7,M))

      CALL FREBUF
      CALL FREREL('S',FSSC(1,1,M),NP)
      CALL FREBUF
      CALL FREREL('E',FSSC(1,2,M),NP)

      FSSC(1,3,M) = FSSC(1,1,M)/FSSC(1,2,M)
      DO 40 I=2,NP
      FSSC(I,3,M) = (FSSC(I,1,M)-FSSC(I-1,1,M))/
	1			  (FSSC(I,2,M)-FSSC(I-1,2,M))
40    CONTINUE

50    CONTINUE

C	-----PRINT MATERIAL PROPERTIES-----

      WRITE (NOT,1100)
      WRITE (NOT,2000) (M,(FSPP(I,M),I=1,6),M=1,NMSL)

      DO 60 M=1,NMSL
      WRITE (NOT,3000) M
	NP = INT(FSPP(7,M))
   60 WRITE (NOT,3100) (I,FSSC(I,1,M),FSSC(I,2,M),FSSC(I,3,M),I=1,NP)


      RETURN

 1000 FORMAT (//,' PRESTRESSING STEEL PROPERTIES INPUT',/,1X,22('='),/)
 1100 FORMAT (   '  TYPE      YIELD      WEIGHT       ALFA       RLXS',
     1           '      FRIC      WOBBLE '/,
     2           '    NO.    STRESS     DENSITY                      ',
	3		   '      COEF.     COEF.'/)


 2000 FORMAT (I7,1P6E11.3)
 3000 FORMAT ( /,' STRESS-STRAIN CURVE POINTS FOR PRESTESSING ',
     1         'STEEL NO.',I3/1X,56('-'),
     2         /,4X,'POINT',12X,'STRESS',6X,'STRAIN',5X,'MODULUS'/)
 3100 FORMAT (I9,6X,1P3E12.3)
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE INPSRK(FSRK,FSKP,NSKP,NPS,NS,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------
C     INPUT SHRINKAGE PATTERN
C	----------------------------------------------------------
      DIMENSION FSRK(NS,NSKP),FSKP(NPS,2,NSKP)

C	FSRK(1  ,M) = NUMBER OF POINT

C     FSKP(I,1,M) = NORMALIZE LEVEL ALONG THINKNESS
C     FSKP(I,2,M) = SHRINKAGE FACTOR


      WRITE(NOT,1000)


C	-----READ MATERIAL PROPERTIES-----

      DO 50 N=1,NSKP

      CALL FREBUF
      CALL FREINT('N',M,1)
      CALL FREREL('P',FSRK(1,M),1)

	NP = INT(FSRK(1,M))

      CALL FREBUF
      CALL FREREL('L',FSKP(1,1,M),NP)
      CALL FREBUF
      CALL FREREL('F',FSKP(1,2,M),NP)

50    CONTINUE

C	-----PRINT MATERIAL PROPERTIES-----

      DO 60 M=1,NSKP
      WRITE (NOT,3000) M
	NP = INT(FSRK(1,M))
   60 WRITE (NOT,3100) (I,FSKP(I,1,M),FSKP(I,2,M),I=1,NP)


      RETURN

 1000 FORMAT (//,' SHRINKAGE PATTERN INPUT'/)

 3000 FORMAT (   4X,'SHRINKAGE PATTERN NO. =',I5/,
	1           4X,'POINT',12X,'LEVEL',6X,'FACTOR'/)
 3100 FORMAT (I9,6X,1P2E12.3)

      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE INPTMP(FTEM,FTMP,NTMP,NPS,NS,NOT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------
C     INPUT TEMPERATURE PATTERN
C	----------------------------------------------------------
      DIMENSION FTEM(NS,NTMP),FTMP(NPS,2,NTMP)

C	FTEM(1  ,M) = NUMBER OF POINT

C     FTMP(I,1,M) = NORMALIZE LEVEL ALONG THINKNESS
C     FTMP(I,2,M) = TEMPERATURE FACTOR


      WRITE(NOT,1000)


C	-----READ MATERIAL PROPERTIES-----

      DO 50 N=1,NTMP

      CALL FREBUF
      CALL FREINT('N',M,1)
      CALL FREREL('P',FTEM(1,M),1)

	NP = INT(FTEM(1,M))

      CALL FREBUF
      CALL FREREL('L',FTMP(1,1,M),NP)
      CALL FREBUF
      CALL FREREL('F',FTMP(1,2,M),NP)

50    CONTINUE

C	-----PRINT MATERIAL PROPERTIES-----

      DO 60 M=1,NTMP
      WRITE (NOT,3000) M
	NP = INT(FTEM(1,M))
   60 WRITE (NOT,3100) (I,FTMP(I,1,M),FTMP(I,2,M),I=1,NP)


      RETURN

 1000 FORMAT (//,' TEMPERATURE PATTERN INPUT'/)

 3000 FORMAT (   4X,'TEMPERATURE PATTERN NO. =',I5/,
	1           4X,'POINT',12X,'LEVEL',6X,'FACTOR'/)
 3100 FORMAT (I9,6X,1P2E12.3)

      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================







