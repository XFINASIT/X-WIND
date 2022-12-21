C	FILD	Non Mechanical Load
C	FDLI	Selfweight Load
C	FXLD	External Joint Load
C	FUBI	Unbalance Joint Load
C	FGRV	Gravity value in each direction


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE MESHINP
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	---------------------------
C	INPUT MODULE    -
C	---------------------------

      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

      COMMON A(9000000),IA(9000000)


C	HEADER
	CALL TITLE 


	CALL INTFILL('%IOL',ITO,1,2,0)  !SCREEN FILE

C	CALL FILE FOR PRINTING  NOT=ISO
	CALL INTFILL('%IOL',NOT,1,3,0)

	CALL INTFILL('%NUB',LEL  ,1,10,0)  !LEL  ELEMENT GROUP POINTER


C	CONCRETE MATERIAL INPUT
	CALL CONCRETE


C	STEEL MATERIAL INPUT	
	NMSL = 0
	CALL FREBUF
      CALL FREINT  ('N',NMSL,1)
	IF(NMSL.GT.0) THEN
	NS = 10
      CALL DEFNREL('%MSL',KMSL,NS,NMSL)
	CALL RELZERO('%MSL')
      CALL  INPMSL(RCONDT(KMSL),NMSL,NS,NOT)
	ENDIF

C	CABLE STEEL MATERIAL INPUT	
	NMSL = 0
	NPSL = 5
	NS   = 10
	CALL FREBUF
      CALL FREINT  ('N',NMSL,1)
	IF(NMSL.GT.0) THEN
      CALL FREINT  ('P',NPSL,1)
      CALL DEFNREL('%CBS',KCBS,NS,NMSL)
      CALL DEFNREL('%CBP',KCBP,NPSL*3,NMSL)  !STRESS-STRAIN-YOUNG
	CALL RELZERO('%CBS')
	CALL RELZERO('%CBP')
	CALL  INPSCB(RCONDT(KCBS),RCONDT(KCBP),NMSL,NPSL,NS,NOT)
	ENDIF

C	PRESTRESSING STEEL MATERIAL INPUT	
	NMSL = 0
	NPSL = 5
	NS   = 12
	CALL FREBUF
      CALL FREINT  ('N',NMSL,1)
	IF(NMSL.GT.0) THEN
      CALL FREINT  ('P',NPSL,1)
      CALL DEFNREL('%CPS',KCPS,NS,NMSL)
      CALL DEFNREL('%CPP',KCPP,NPSL*3,NMSL)  !STRESS-STRAIN-YOUNG
	CALL RELZERO('%CPS')
	CALL RELZERO('%CPP')
	CALL  INPSPS(RCONDT(KCPS),RCONDT(KCPP),NMSL,NPSL,NS,NOT)
	ENDIF


C	FRAME SECTION INPUT
	CALL FSECTION

C	SHELL SECTION INPUT
	CALL SSECTION

C	SOLID SECTION INPUT
	CALL SOSECTION		
	
C	TRUSS SECTION INPUT
	CALL TSECTION

C	PCABLE SECTION INPUT
	CALL PSECTION

C	CCABLE SECTION INPUT
	CALL CSECTION

C	LINK PROPERTIES INPUT
	CALL LINKCSA

C	TRAVELER PROPERTIES INPUT
	CALL TRAVCSA

C	SPRING SUPPORT
	CALL SPIGINIT

	WRITE(ITO,*) 'INITIALIZE ELEMENT DATA'
      WRITE(10,*) 'INITIALIZE ELEMENT DATA'
C	ELEMENT DATA INPUT
	CALL INPGRP(IA(LEL))
	WRITE(ITO,*) 'FINISH INITIALIZE ELEMENT DATA'
      WRITE(10,*) 'FINISH INITIALIZE ELEMENT DATA'

	WRITE(ITO,*) 'INITIALIZE TENDON DATA'
      WRITE(10,*) 'INITIALIZE TENDON DATA'
C	FRAME TENDON DATA INPUT
	CALL INPTENF
	WRITE(ITO,*) 'FINISH INITIALIZE TENDON DATA'
      WRITE(10,*) 'FINISH INITIALIZE TENDON DATA'

C	--------------------------------------


      RETURN
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE TITLE
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	--------------------------------------------------
C     PRINT PROGRAM HEADER AND TITLE IN OUTPUT    -
C	--------------------------------------------------

C	FILE FOR PRINTING  NOT=ISO
	CALL INTFILL('%IOL',NOT,1,3,0)

C	HEADER
      WRITE(NOT,3000)
 
	RETURN
 3000 FORMAT(//1x,79('=')/1x,'=',77x,'=',
     1        /1x,'=',31x,'  XFINAS-FEM   ',31x,'='/1x,'=',77x,'=',
     2        /1x,'=',12x, 'NONLINEAR ANALYSIS OF SEGMENTALLY ',
     3         'ERECTED STRUCTURES',13x,'=',
     4        /1x,'=',77x,'=',/1x,'=',25x,'VERSION 4.00,      DEC 2007',
     5        25x,'=',/1x,'=',77x,'='/1x,79('=')//)

      END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE INPXYZ (XYZ,NSC,NSN) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	-----------------------------------
C     TRANSFER NODE COORDINATES    -
C	-----------------------------------

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON A(9000000),IA(9000000)
      DIMENSION XYZ(NSC,1)

C-----TRANSFER COORDINATE DATA-----
	IXY = LXY
	DO 20 I = 1,NSC
	DO 10 J = 1,NSN
	XYZ(I,J) = A(IXY)
	IXY = IXY+1
10	CONTINUE
20	CONTINUE

	RETURN
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE INPLAX (VLAX,NLX,NLA) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	-----------------------------------
C     LOCAL AXIS
C	-----------------------------------

      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON A(9000000),IA(9000000)
      DIMENSION VLAX(NLX,1)

C-----TRANSFER COORDINATE DATA-----
	IDC = LDC
	DO 20 I = 1,NLA
	DO 10 J = 1,NLX
	VLAX(J,I) = A(IDC)
	IDC = IDC+1
10	CONTINUE
20	CONTINUE

	RETURN
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE INPLID (ID,NSF,NSN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	----------------------------------------------------
C     FULL EQUATION NUMBERING  (NO SUPPORT CONDITION HERE)
C	----------------------------------------------------
      DIMENSION ID(NSF,NSN)

	IEQ = 0
      DO ISN=1,NSN
      DO ISF=1,NSF
	IEQ = IEQ + 1
	ID(ISF,ISN) = IEQ
	ENDDO
	ENDDO 

      RETURN
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE INPSEQ (NSEQ,NSN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------
C     INPUT NODE SEQUENCE 
C	(CAN USE IF GID NODE ARE NOT THE SAME WITH STRUC. NODE)
C	-------------------------------------------------------
      DIMENSION NSEQ(1)

      DO 10 J=1,NSN
   10 NSEQ(J) = J

      RETURN
      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE INPGRP(LEST)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	COMMON /GiDEle/ LGID

      COMMON A(9000000),IA(9000000)

      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

      COMMON /LOCO/ LOP,LOS,LSS,LSS2,LSS3,LHG,LHGN
      COMMON /OFFS/ NOPS,NHIGE

      DIMENSION LEST(1)
C	----------------------------------------------------------------------------------------
C	----------------------------------------------------------------------------------------
C	GROUP LEVEL
C
C	ELEMENT LEVEL
C	CONNECTIVITY	- (N1,NELE)	 				--> I NNM NODES
C	CONNECTIVITY	- (N2,NELE)	 				--> I ORI EQU NO. + UPD EQU NO. + UPD RES EQU NO. -------UPD = UPDATED N2 = 3*NEF
C	SECTION			- (N3,NELE)	 				--> I IGIDM + SECTION FLAG + CREEP INT FLAG + TOTAL EXIST FLAG + SECTYPE + NUMBER OF PART + FIRST INSTALL DAY + (RELEASE FLAG FOR FRAME 14)
C	PART			- (N4,NPRT,NELE)	        --> R EXISTING FLAG + NUM FIBER + CAST DATE + INSTALL DATE 
C	PART			- (N5,NFIB,NGAS,NPRT,NELE)	--> R WORKING ARRAY FOR EACH GAUSS AND FIBER
C	NODAL			- (N6,NELE)				    --> R INITIAL DISP + RESISTING + RESISTING NON-MCHA + FIXEND + (6 OFFSET FOR FRAME) 
C	GAUSS           - (N7,NGAS,NELE)            --> R WORKING ARRAY FOR A GAUSS POINT
C	FIXEND			- (N8,NELE)					--> R DEAD, TRAVELER, SETTLEMENT, TEMPERATURE, CONSTRUCTION, TENDON PR, TENDON SC, CREEP+SHRN PR, CREEP+SHRN SC, TOTAL
C	----------------------------------------------------------------------------------------
C	----------------------------------------------------------------------------------------
	

C	CALLING NEG
	CALL INTFILL('%NUB',NEG,1,5,0)

C	ADDRESS FOR GROUP
	CALL DEFNINT('%GRP',NGRP,30,NEG+3+1)  !GROUP DATA   !+3 FOR FRAME+SHELL+SOLID TENDON + 1 FOR LINK
	CALL DEFNINT('%GRF',NGRF,40,NEG+3+1)  !GROUP FILE

C	GID ELEMENT MAPPING ARRAY   
C	NGIDM = 100  
	CALL INTFILL('%IOL',ITI,1,1,0)   !ITI FILE	
	CALL FREBUF
	CALL FREINT('N',NGIDM,1)   !GID ELEMENT NUMBER FOR CONSTRUCTION ANALYSIS

C	----------------------------------
C	FIND THE MAXIMUM NUMBER OF GID ELEMENT FOR STORAGE
	NGIDX = 0
	DO IEG = 1,NEG
	NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG    )-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)      
      
      CALL MOVLEV (2)
	DO IELE = 1,NELE
	IGM = IA(LGID+IELE-1)
	IF(IGM.GT.NGIDX) NGIDX = IGM
	ENDDO
	ENDDO
C	----------------------------------

	CALL DEFNINT('%GDM',KGDM,2,NGIDX) 
	CALL INTZERO('%GDM')

C	----------------------------------------
C	LOOP OVER ELEMENT GROUP 
C	----------------------------------------	
	DO 2500 IEG = 1,NEG
	NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG    )-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     
      CALL MOVLEV (2)


	SELECTCASE(ITYPE)

	CASE(2)
	NEL  = NELE
	MPRT = 1
	MGAS = 1
	MFIB = 1 
	N1   = NNM
	N2   = NEF*3
	N3   = 10
	N4   = 10
	N5   = 20
	N6   = NEF*4
	N7   = 100        !WORKING ARRAY FOR A GAUSS POINT
	N8	 = NEF*10*3   !3 For Current Old and Incremental
	N9	 = NEF*10*3   !3 For Current Old and Incremental
	N10  = NEF		  !FOR INITIAL DISP OF EACH PART JAN09
	NGR  = MGAS
	NGS  = 1
	NGT  = 1

	CASE(5) ! FRAME
	NEL  = NELE
	CALL INTFILL('$FSC',MPRT,1,4,0) !MAX PART
	CALL INTFILL('$FSC',MGAS,1,6,0)	!NUM GAUS
	CALL INTFILL('$FSC',MFIB,1,7,0) !MAX FIBER IN EACH PART
	N1   = NNM
	N2   = NEF*3
	N3   = 7 + 14          !14 FOR RELEASE FLAG
	N4   = 10 
	N5   = 200
	N6   = NEF*4 + 6       !6 FOR OFFSET  
	N7   = 20 + 8 + 8      !20 FOR RIGIDITY  8 FOR MEC STRESSES  8 FOR NON-MEC STRESSES
	N8	 = NEF*10*3   !3 For Current Old and Incremental
	N9	 = NEF*10*3   !3 For Current Old and Incremental
	N10  = NEF		  !FOR INITIAL DISP OF EACH PART JAN09
	NGR  = MGAS
	NGS  = 1
	NGT  = 1

C	STRESS CALCULATION OF SECTION PART OF FRAME ELEMENT
	MXSTRP = 12  !MAXIMUM STRESS POINT PER PART   NUM-S-T
	CALL DEFNREL('+FCS',KFCS,MXSTRP*2*MPRT*10,NELE)   !SECTION ROTATION DATA  2 FOR NODE I-J
	CALL RELZERO('+FCS')

	CASE(6)
	NEL  = NELE
	MPRT = 1
	MGAS = NPT
	MFIB = 1 
	N1   = NNM
	N2   = NEF*3
	N3   = 10
	N4   = 10
	N5   = 2
	N6   = NEF*4
	N7   = 3        !WORKING ARRAY FOR A GAUSS POINT
	N8	 = NEF*10*3   !3 For Current Old and Incremental
	N9	 = NEF*10*3   !3 For Current Old and Incremental
	N10  = NEF		  !FOR INITIAL DISP OF EACH PART JAN09

	CASE(9) ! SHELL
	NGR  = 2
	NGS  = 2
	NGT  = 1
	NEL  = NELE
	MPRT = 1
	MGAS = NGR*NGS
	MFIB = 14         !MAXIMUM LAYER LIMITED TO 22
	N1   = NNM
	N2   = NEF*3
	N3   = 10
	N4   = 10
	N5   = 230
	N6   = NEF*4
	N7   = 40 + 8 + 8 + 18       !WORKING ARRAY FOR A GAUSS POINT (18 FOR 2 SET OF VR VS VT)  1 FOR DIRLLING STRAIN
	N8	 = NEF*10*3   !3 For Current Old and Incremental
	N9	 = NEF*10*3   !3 For Current Old and Incremental
	N10  = NEF		  !FOR INITIAL DISP OF EACH PART JAN09

	N7 = N7 + 1         !1 FOR DRILLING OLD STRAIN
	N7 = N7 + 1         !1 FOR DRILLING OLD STRESS
	N7 = N7 + 1         !1 FOR DRILLING STRESS
	
	N5 = N5 + 5*11		!FOR LOAD COMBINATION OF FIBER STRESSES		10 LOADCASE  JAN09 LOADCOMB  PLUS 1 case for reference case
	N7 = N7 + 8*11		!FOR LOAD COMBINATION OF STRESS RESULTANTS	10 LOADCASE  JAN09 LOADCOMB  PLUS 1 case for reference case 


	CASE(10)   ! SOLID ELEMENT BY PRAMIN DEC 2010
	
	NGR  = 2           !NUMBER OF GAUSS FOR R-DIRECTION
	NGS  = 2           !NUMBER OF GAUSS FOR S-DIRECTION 
	NGT  = 2           !NUMBER OF GAUSS FOR T-DIRECTION
	NEL  = NELE        !NUMBER OF ELEMENT
	MPRT = 1           !MAXIMUM PART (MPRT = 1,THERE IS NO ANY PART) 
	MGAS = NGR*NGS*NGT !MAXIMUM GAUSS POINT (8 GAUSS POINTS) 
	MFIB = 1           !MAXIMUM LAYER (MFIB = 1, THERE IS NO ANY LAYER)
	N1   = NNM         !NUMBER OF NODE / ELEMENT (8 NODES)
	N2   = NEF*3       !3TIMES OF NUMBER OF DOF PER ELEMENT (NEF=24 DOFS)
	N3   = 10          !JUST WANT TO KEEP 10 DATA 
	N4   = 10          !JUST WANT TO KEEP 10 DATA 
	N5   = 300         !JUST WANT TO KEEP 300 DATA 
	N6   = NEF*4       !4TIMES OF NUMBER OF DOF PER ELEMENT (NEF=24 DOFS)
	N7   = 12 + 6 + 6 + 18  !WORKING ARRAY FOR A GAUSS POINT AND 18 FOR 2 SET (VR VS VT)
	N8	 = NEF*10*3   !3 For Current Old and Incremental
	N9	 = NEF*10*3   !3 For Current Old and Incremental
	N10  = NEF		  !FOR INITIAL DISP OF EACH PART JAN09

	N7 = N7 + 6*11		!FOR LOAD COMBINATION OF GAUSS POINT STRESSES		10 LOADCASE  JAN09 LOADCOMB  PLUS 1 case for reference case
	

	CASE(17) ! TENDON ELEMENT
	NEL  = 0 !NELE -- NO TENDON ELEMENT HERE
	MPRT = 1
	MGAS = NPT
	MFIB = 1 
	N1   = NNM
	N2   = NEF*3
	N3   = 10
	N4   = 10
	N5   = 6
	N6   = NEF*4
	N7   = 3          !WORKING ARRAY FOR A GAUSS POINT
	N8	 = NEF*10*3   !3 For Current Old and Incremental
	N9	 = NEF*10*3   !3 For Current Old and Incremental
	N10  = NEF		  !FOR INITIAL DISP OF EACH PART JAN09

	ENDSELECT


C     FILE FOR INTEGER DATA PER ONE ELEMENT
      NDX1 = 1
      NDX2 = NDX1 + N1
      NDX3 = NDX2 + N2
      LASI = NDX3 + N3
      
C     FILE FOR REAL DATA PER ONE ELEMENT
      NDX4 = 1
      NDX5 = NDX4 + N4*MPRT
      NDX6 = NDX5 + N5*MFIB*MGAS*MPRT
      NDX7 = NDX6 + N6
      NDX8 = NDX7 + N7*MGAS
      NDX9 = NDX8 + N8
      NDX10= NDX9 + N9
      LASR = NDX10+ N10*MPRT

C	OPEN FILE FOR STORE GROUP DATA
	NFIL = 3000
	NFLL = 2
	NFLI = 1 + NFLL*(IEG-1) + NFIL !INTEGER FILE
	NFLR = 2 + NFLL*(IEG-1) + NFIL !REAL FILE

C	GROUP DATA
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,1) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,1) !
	CALL INTFILL('%GRP',NEL  ,3 ,IEG,1) !NEL -- INSTEAD OF NELE  --- FOR TENDON LET NEL = 0
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,1) !
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,1) !
	CALL INTFILL('%GRP',MFIB ,6 ,IEG,1) !
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,1) !
	CALL INTFILL('%GRP',NNM  ,8 ,IEG,1) !
	CALL INTFILL('%GRP',NNF  ,9 ,IEG,1) !
	CALL INTFILL('%GRP',NGR  ,10,IEG,1) !
	CALL INTFILL('%GRP',NGS  ,11,IEG,1) !
	CALL INTFILL('%GRP',NGT  ,12,IEG,1) !
C	GROUP FILE
	CALL INTFILL('%GRF',N1   ,1 ,IEG,1) !
	CALL INTFILL('%GRF',N2   ,2 ,IEG,1) !
	CALL INTFILL('%GRF',N3   ,3 ,IEG,1) !
	CALL INTFILL('%GRF',N4   ,4 ,IEG,1) !
	CALL INTFILL('%GRF',N5   ,5 ,IEG,1) !
	CALL INTFILL('%GRF',N6   ,6 ,IEG,1) !
	CALL INTFILL('%GRF',N7   ,7 ,IEG,1) !
	CALL INTFILL('%GRF',N8   ,8 ,IEG,1) !
	CALL INTFILL('%GRF',N9   ,9 ,IEG,1) !
	CALL INTFILL('%GRF',N10  ,10,IEG,1) !

	CALL INTFILL('%GRF',NDX1 ,11,IEG,1) !
	CALL INTFILL('%GRF',NDX2 ,12,IEG,1) !
	CALL INTFILL('%GRF',NDX3 ,13,IEG,1) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,1) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,1) !
	CALL INTFILL('%GRF',NDX6 ,16,IEG,1) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,1) !
	CALL INTFILL('%GRF',NDX8 ,18,IEG,1) !
	CALL INTFILL('%GRF',NDX9 ,19,IEG,1) !
	CALL INTFILL('%GRF',NDX10,20,IEG,1) !
	
	CALL INTFILL('%GRF',NFLI ,21,IEG,1) !
	CALL INTFILL('%GRF',NFLR ,22,IEG,1) !
	
	CALL INTFILL('%GRF',LASI ,23,IEG,1) !
	CALL INTFILL('%GRF',LASR ,24,IEG,1) !


	CALL DIROPEN(NFLI,2*LASI)     !INTEGER FILE
	CALL DIROPEN(NFLR,2*LASR)     !REAL FILE


	DO 2000 IEL = 1,NEL !NEL -- INSTEAD OF NELE  --- FOR TENDON LET NEL = 0

	IGM = IA(LGID+IEL-1)

C	STORE GID ELEMENT MAPPING DATA
	CALL INTFILL('%GDM',IEG,1,IGM,1)
	CALL INTFILL('%GDM',IEL,2,IGM,1)

2000	CONTINUE

2500	CONTINUE

C	----------------------------------------
C	LOOP OVER ELEMENT GROUP 
C	----------------------------------------	
	DO 5000 IEM = 1,NGIDM

	ISEC = 0
	ICRP = 1
	CALL FREBUF
	CALL FREINT('E',IGM ,1)   !GID ELEMENT NO
	CALL FREINT('S',ISEC,1)   !SECTION FLAG

	CALL INTFILL('%GDM',IEG,1,IGM,0)
	CALL INTFILL('%GDM',IEL,2,IGM,0)

	NELEMI = 10 + IEG
      NELEMA = 30 + IEG
      REWIND NELEMI
      REWIND NELEMA
C      READ (NELEMI) (IA(NLNU),NLNU=LNU,LNU + LEST(IEG)-1)
C      READ (NELEMA) ( A(NLNU),NLNU=LMP,LMP + LEST(IEG+NEG)-1)
      READ (NELEMI) IA(LNU:LNU + LEST(IEG    )-1)
      READ (NELEMA)  A(LMP:LMP + LEST(IEG+NEG)-1)     
      CALL MOVLEV (2)

	CALL LOCATN('$LID',KLID,NSF,NSN,1)

	CALL INPELM(IA(LCN),ICONDT(KLID),NSF,NSN,
	1			IA(LHG),IA(LHGN),IA(LOS),A(LOP),NEF,
     1			IEG,IEL,IGM,ISEC)                !LHG LHGN LOS LOP FOR HINGE AND OFFSET 


5000	CONTINUE
C	----------------------------------------

C	----------------------------------------
C	LOOP OVER ELEMENT GROUP 
C	----------------------------------------	
	DO IEG = 1,NEG

C	ELEMENT WORKING ARRAY INITIALIZATION
	CALL  ELLWOK(IEG)


	ENDDO
C	----------------------------------------


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE INPELM(LM,ID,NSF,NSN,
	1				  LREAS,IHSET,IOSET,PROPO,NEF,
     1				  IEG,IEL,IGM,ISEC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	ELEMENT LINK POSITION AND DOF 
	COMMON /EFLINK/ NFLINK(30,30)

      DIMENSION LM(NEF,1),ID(NSF,NSN),IDOF(9),IGPOS(9)
	DIMENSION NPM(10),NPI(10)                     !TEMPORARY ARRAY

	DIMENSION LREAS(14,1),IHSET(1),IOSET(1),PROPO(6,1)   
	
	ALLOCATABLE IF1(:),IF2(:),IF3(:),AF4(:),AF5(:),AF6(:),AF7(:)       !TEMPORARY ARRAY
	ALLOCATABLE AF8(:),AF9(:),AF10(:)

	ALLOCATABLE IAX(:),AAX(:) !TEMPORARY ARRAY


	IF(IEL.LE.0) GOTO 5000

C	----------------------------------------
C	CALLING IDOF
C	----------------------------------------
	DO I = 1,9
	CALL INTFILL('%DOF',IDOF(I),1,I,0)
	ENDDO

C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',MPRT ,4 ,IEG,0) !
	CALL INTFILL('%GRP',MGAS ,5 ,IEG,0) !
	CALL INTFILL('%GRP',MFIB ,6 ,IEG,0) !
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	CALL INTFILL('%GRP',NNM  ,8 ,IEG,0) !
	CALL INTFILL('%GRP',NNF  ,9 ,IEG,0) !
	CALL INTFILL('%GRP',NGR  ,10,IEG,0) !
	CALL INTFILL('%GRP',NGS  ,11,IEG,0) !
	CALL INTFILL('%GRP',NGT  ,12,IEG,0) !
C	GROUP FILE
	CALL INTFILL('%GRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('%GRF',N2   ,2 ,IEG,0) !
	CALL INTFILL('%GRF',N3   ,3 ,IEG,0) !
	CALL INTFILL('%GRF',N4   ,4 ,IEG,0) !
	CALL INTFILL('%GRF',N5   ,5 ,IEG,0) !
	CALL INTFILL('%GRF',N6   ,6 ,IEG,0) !
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !
	CALL INTFILL('%GRF',N8   ,8 ,IEG,0) !
	CALL INTFILL('%GRF',N9   ,9 ,IEG,0) !
	CALL INTFILL('%GRF',N10  ,10,IEG,0) !

	CALL INTFILL('%GRF',NDX1 ,11,IEG,0) !
	CALL INTFILL('%GRF',NDX2 ,12,IEG,0) !
	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX6 ,16,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
	CALL INTFILL('%GRF',NDX8 ,18,IEG,0) !
	CALL INTFILL('%GRF',NDX9 ,19,IEG,0) !
	CALL INTFILL('%GRF',NDX10,20,IEG,0) !
	
	CALL INTFILL('%GRF',NFLI ,21,IEG,0) !
	CALL INTFILL('%GRF',NFLR ,22,IEG,0) !
	
	CALL INTFILL('%GRF',LASI ,23,IEG,0) !
	CALL INTFILL('%GRF',LASR ,24,IEG,0) !
	
	ALLOCATE (IF1(N1),IF2(N2),IF3(N3),AF4(N4),AF5(N5),AF6(N6),AF7(N7))
	ALLOCATE (AF8(N8),AF9(N9),AF10(N10))

	ALLOCATE (IAX(LASI),AAX(LASR))
      IAX(1:LASI) = 0
      AAX(1:LASR) = 0.0D0
      
C	----------------------------------------
C     UNPACK VARIABLE LINKF
C	----------------------------------------
      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLPOS(IDOF,LINK,NNF,IGPOS)
C	----------------------------------------


C	----------------------------------------
C	FIRST FILE INTG
	IF1(1:N1) = LM(1:N1,IEL)                        !CONNECTIVITY  N1 = NNM 
	NPM(1)  = N1
	NPI(1)  = N1
	CALL CALLINTG(IF1,IAX,NDX1,NPM,NPI,1,'WRIT')              !INTEGER DATA

C	----------------------------------------
C	SECOND FILE	INTG
	IF2(1:N2) = 0
	DO INN = 1,NNM
	NIJ = LM(INN,IEL)
	DO INF = 1,NNF
	IGF = IGPOS(INF)
	IF(IGF.NE.0) IF2(INF+NNF*(INN-1)) = ID(IGF,NIJ) !ORI EQU NO. + UPD EQU NO. + UPD RES EQU NO. -------UPD = UPDATED N2 = 3*NEF
	ENDDO
	ENDDO
	NPM(1) = N2
	NPI(1) = N2
	CALL CALLINTG(IF2,IAX,NDX2,NPM,NPI,1,'WRIT')              !INTEGER DATA

C	----------------------------------------
C	EIGHTH FILE  REAL	
	AF8(1:N8) = 0.0D0                               !FIXEND FORCE
	NPM(1)  = N8
	NPI(1)  = N8
	CALL CALLREAL(AF8,AAX,NDX8,NPM,NPI,1,'WRIT')  
C	----------------------------------------
C	NINTH FILE  REAL	
	AF9(1:N9) = 0.0D0                               !FIXEND FORCE
	NPM(1)  = N9
	NPI(1)  = N9
	CALL CALLREAL(AF9,AAX,NDX9,NPM,NPI,1,'WRIT') 
C	----------------------------------------


C	========================================
	SELECT CASE (ITYPE)

C	========================
	CASE(2)  !TRUSS CABLE 
C	----------------------------------------
C	THIRD FILE INTG
	IF3(1) = IGM                             !GID ELEMENT NUMBER     
	IF3(2) = ISEC                                   !SECTION NO.
C	IF3(3) = ICRP                                   !CREEP INTEGRATION FLAG  (ALREADY IN MATERIAL)
	IF3(4) = 0	                                    !EXISTING FLAG 
C	IF3(5) = IMAT                                   !NOT USe FOr TRUSS AND CABLE (ALREADY IN SECTION)
	IF3(6) = MPRT                                   !NUMBER OF PART NOT USE FOR NON-FRAME ELEMENT
	IF3(7) = -1                                     !1st INSTALLATION DAY
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	FOURTH FILE  REAL  FOR EACH PART
	DO IPRT = 1,MPRT 
	AF4(1) = 0.0 !PART EXISTING FLAG
	AF4(2) = FLOAT(MFIB)
	AF4(3) = 0.0 !CASTING DAY
	AF4(4) = 0.0 !INSTALLATION DAY
C	AF4(5) = 0.0 !MASS DENSITY (NOT USE)
	AF4(6) = 0.0 !FORCE FOR CABLE PRE-TENSION OR INITIAL LENGTH OF CATENARY
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')
	ENDDO
C	----------------------------------------
C	FIFTH FILE  REAL
	AF5(1:N5) = 0.0D0                              !WORKING ARRAY FOR EACH FIBER  NOT USE FOR NON-FRAME ELEMENT
	DO IPRT = 1,MPRT 
	DO IGAS = 1,MGAS
	DO IFIB = 1,MFIB
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,IFIB,IGAS,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')
	ENDDO
	ENDDO
	ENDDO
C	----------------------------------------
C	SIXTH FILE  REAL	
	AF6(1:N6) = 0.0D0                              !INITIAL DISP. + RESISTING FORCE + RESISTING FORCE NON-MCHA + FIXEND FORCE
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	SEVENTH FILE  REAL  WORKING ARRAY FOR EACH GAUSS	
	AF7(1:N7) = 0.0D0                             !WORKING ARRAY 
	DO IGAS = 1,MGAS
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IGAS]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT') 
	ENDDO

C	========================
	CASE(5)  !FRAME

	CALL INTFILL('FPT1',M   ,1,ISEC,0)    !SEC TYPE FLAG -- READY OR FIBER          
	CALL INTFILL('FPT1',NPRT,2,ISEC,0)    !NUMBER OF PART

C	----------------------------------------
C	THIRD FILE INTG
	IF3(1) = IGM                             !GID ELEMENT NUMBER             
	IF3(2) = ISEC                                   !SECTION NO.
C	IF3(3) = ICRP                                   !CREEP INTEGRATION FLAG  (ALREADY IN MATERIAL)
	IF3(4) = 0	                                    !EXISTING FLAG
	IF3(5) = M                                      !SEC TYPE FLAG -- READY OR FIBER
	IF3(6) = NPRT								    !NUMBER OF PART
	IF3(7) = -1                                     !1st INSTALLATION DAY
	IHET = IHSET(IEL)                               !OFFSET NO.
	IF3(8:21) = 0
	IF(IHET.NE.0) IF3(8:21) = LREAS(1:14,IHET)                     !7-21 14 RELEASE FLAG
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	FOURTH FILE  REAL
	DO IPRT = 1,NPRT 
	CALL INTFILL('FPT2',NPIB,IPRT,ISEC,0)           !NUMBER OF FIBER FOR THIS PART
	AF4(1) = 0.0D0			!PART EXISTING FLAG
	AF4(2) = FLOAT(NPIB)    !NUMBER OF FIBER			FOR THIS PART
	AF4(3) = 0.0D0			!CASTING DAY				FOR THIS PART
	AF4(4) = 0.0D0			!INSTALLATION DAY			FOR THIS PART
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')
	ENDDO
C	----------------------------------------
C	FIFTH FILE  REAL
	AF5(1:N5) = 0.0D0                               !WORKING ARRAY FOR EACH FIBER
	DO IPRT = 1,NPRT 
	CALL INTFILL('FPT2',NPIB,IPRT,ISEC,0)           !NUMBER OF FIBER FOR THIS PART
	DO IGAS = 1,MGAS
	DO IFIB = 1,NPIB
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,IFIB,IGAS,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')
	ENDDO
	ENDDO
	ENDDO
C	----------------------------------------
C	SIXTH FILE  REAL	
	AF6(1:N6) = 0.0D0                               !INITIAL DISP. + RESISTING FORCE + RESISTING FORCE - NON MCHA + FIXEND FORCE + 6 OFFSET COMPONENT 
	NUM1 = 1 + 4*NEF
	NUM2 = 6 + 4*NEF
	IOET = IOSET(IEL)                               !OFFSET NO.
	IF(IOET.NE.0) AF6(NUM1:NUM2) = PROPO(1:6,IOET)  !OFFSET VALUE.
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	SEVENTH FILE  REAL STRESS AND RIGIDITY FOR EACH GAUSS	
	AF7(1:N7) = 0.0D0                               !8 FOR RIGIDITY  8 FOR MEC STRESSES  8 FOR NON-MEC STRESSES
	DO IGAS = 1,MGAS
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IGAS]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT') 
	ENDDO
C	----------------------------------------
C	TENTH FILE  REAL
	AF10(1:N10) = 0.0D0                             !INITIAL DISP FOR EACH PART
	DO IPRT = 1,NPRT 
	NPM(1:2)  = [N10,MPRT]
	NPI(1:2)  = [N10,IPRT]
	CALL CALLREAL(AF10,AAX,NDX10,NPM,NPI,2,'WRIT') 
	ENDDO
C	========================
	CASE(9)  !SHELL

C	----------------------------------------
C	THIRD FILE INTG
	IF3(1) = IGM									!GID ELEMENT NUMBER             
	IF3(2) = ISEC                                   !SECTION NO.
C	IF3(3) = ICRP                                   !CREEP INTEGRATION FLAG  (ALREADY IN MATERIAL)
	IF3(4) = 0	                                    !EXISTING FLAG
C	IF3(5) = IMAT                                   !MATERIAL No. - NOT USE FOR SHELL (ALREADY IN SECTION)
	IF3(6) = MPRT								    !NUMBER OF PART (ONLY 1 PART FOR SHELL)
	IF3(7) = -1                                     !1st INSTALLATION DAY
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'WRIT')

C	----------------------------------------
C	FOURTH FILE  REAL
	DO IPRT = 1,MPRT 
	CALL RELFILL('$SSC',FLAY,4,ISEC,0)           !NUMBER OF LAYER FOR THIS PART (ONLY 1 PART FOR SHELL)
	AF4(1) = 0.0D0			!PART EXISTING FLAG
	AF4(2) = FLAY			!NUMBER OF LAYER			FOR THIS PART
	AF4(3) = 0.0D0			!CASTING DAY				FOR THIS PART
	AF4(4) = 0.0D0			!INSTALLATION DAY			FOR THIS PART
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')
	ENDDO

C	----------------------------------------
C	FIFTH FILE  REAL
	AF5(1:N5) = 0.0D0                               !WORKING ARRAY FOR EACH FIBER
	DO IPRT = 1,MPRT 
	CALL RELFILL('$SSC',FLAY,4,ISEC,0)				!NUMBER OF LAYER FOR THIS PART (ONLY 1 PART FOR SHELL)
	NLAY = INT(FLAY)
	DO IGAS = 1,MGAS
	DO ILAY = 1,NLAY
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,ILAY,IGAS,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')
	ENDDO
	ENDDO
	ENDDO
C	----------------------------------------
C	SIXTH FILE  REAL	
	AF6(1:N6) = 0.0D0                               !INITIAL DISP. + RESISTING FORCE + RESISTING FORCE - NON MCHA + FIXEND FORCE 
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	SEVENTH FILE  REAL  WORKING ARRAY FOR EACH GAUSS
	AF7(1:N7) = 0.0D0                               !40 FOR RIGIDITY  8 FOR MEC STRESSES  8 FOR NON-MEC STRESSES
	DO IGAS = 1,MGAS
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IGAS]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT') 
	ENDDO
	


C	========================
	CASE(10)  !SOLID BY PRAMIN DEC 2010
C	========================
C	THIRD FILE INTG : GENERAL SECTION DATA (ELEMENT,SECTION,CREEP,TYPE)
	IF3(1) = IGM                  !GID ELEMENT NUMBER     
	IF3(2) = ISEC                 !SECTION NO.
C	IF3(3) = ICRP                 !CREEP INTEGRATION FLAG  (ALREADY IN MATERIAL)
	IF3(4) = 0	                !EXISTING FLAG 
C	IF3(5) = IMAT                 !NOT USe FOr TRUSS AND CABLE (ALREADY IN SECTION)
	IF3(6) = MPRT                 !NUMBER OF PART NOT USE FOR NON-FRAME ELEMENT
	IF3(7) = -1                   !1st INSTALLATION DAY
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	FOURTH FILE  REAL  FOR EACH PART : CASTING DAY, INSTALLATION DAY
	DO IPRT = 1,MPRT 
	AF4(1) = 0.0D0 !PART EXISTING FLAG
	AF4(2) = FLOAT(MFIB) 
	AF4(3) = 0.0D0 !CASTING DAY
	AF4(4) = 0.0D0 !INSTALLATION DAY
	!AF4(5) = 0.0D0 !MASS DENSITY (NOT USE)
	!AF4(6) = 0.0D0 !FORCE FOR CABLE PRE-TENSION OR INITIAL LENGTH OF CATENARY
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')
	ENDDO
C	----------------------------------------
C	FIFTH FILE  REAL : WORKING ARRAY FOR LAYER/FIBER (STRAIN,CREEP,SHRINKAGE)
	AF5(1:N5) = 0.0D0         !WORKING ARRAY FOR EACH FIBER  NOT USE FOR NON-FRAME ELEMENT
	DO IPRT = 1,MPRT 
	DO IGAS = 1,MGAS
	DO IFIB = 1,MFIB 
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,IFIB,IGAS,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')
	ENDDO
	ENDDO
	ENDDO
C	----------------------------------------
C	SIXTH FILE  REAL : INITIAL DISP, RESISTING FORCE
	AF6(1:N6) = 0.0D0         !INITIAL DISP. + RESISTING FORCE + RESISTING FORCE NON-MCHA + FIXEND FORCE
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	SEVENTH FILE  REAL  WORKING ARRAY FOR EACH GAUSS(STRESS,STRIAN,RIGIDITY)
	AF7(1:N7) = 0.0D0         !WORKING ARRAY 
	DO IGAS = 1,MGAS
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IGAS]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT') 
	ENDDO	
	
	

C	========================
	CASE DEFAULT  !PLANE SHELL SOLID  !SOON SHELL WILL MODIFIED
C	----------------------------------------
C	THIRD FILE INTG
	IF3(1) = IGM                             !GID ELEMENT NUMBER     
	IF3(2) = ISEC                                   !KEEP HERE BUT NOT USE FOR NON-FRAME ELEMENT
C	IF3(3) = ICRP                                   !CREEP INTEGRATION FLAG
	IF3(4) = 0	                                    !EXISTING FLAG 
C	IF3(5) = IMAT                                   !MATERIAL FLAG FOR NON-FRAME ELEMENT
	IF3(6) = MPRT                                   !NUMBER OF PART NOT USE FOR NON-FRAME ELEMENT
	IF3(7) = -1                                     !1st INSTALLATION DAY
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	FOURTH FILE  REAL  FOR EACH PART
	DO IPRT = 1,MPRT 
	AF4(1) = 0.0 !PART EXISTING FLAG
	AF4(2) = FLOAT(MFIB)
	AF4(3) = 0.0 !CASTING DAY
	AF4(4) = 0.0 !INSTALLATION DAY
	AF4(5) = 0.0 !MASS DENSITY
	NPM(1:2) = [N4,MPRT]
	NPI(1:2) = [N4,IPRT]
	CALL CALLREAL(AF4,AAX,NDX4,NPM,NPI,2,'WRIT')
	ENDDO
C	----------------------------------------
C	FIFTH FILE  REAL
	AF5(1:N5) = 0.0D0                              !WORKING ARRAY FOR EACH FIBER  NOT USE FOR NON-FRAME ELEMENT
	DO IPRT = 1,MPRT 
	DO IGAS = 1,MGAS
	DO IFIB = 1,MFIB
	NPM(1:4) = [N5,MFIB,MGAS,MPRT]
	NPI(1:4) = [N5,IFIB,IGAS,IPRT]
	CALL CALLREAL(AF5,AAX,NDX5,NPM,NPI,4,'WRIT')
	ENDDO
	ENDDO
	ENDDO
C	----------------------------------------
C	SIXTH FILE  REAL	
	AF6(1:N6) = 0.0D0                              !INITIAL DISP. + RESISTING FORCE + RESISTING FORCE NON-MCHA + FIXEND FORCE
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')
C	----------------------------------------
C	SEVENTH FILE  REAL  WORKING ARRAY FOR EACH GAUSS	
	AF7(1:N7) = 0.0D0                             !WORKING ARRAY 
	DO IGAS = 1,MGAS
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IGAS]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'WRIT') 
	ENDDO

	ENDSELECT
C	========================================


      WRITE(NFLI,REC=IEL) IAX
      WRITE(NFLR,REC=IEL) AAX

	DEALLOCATE (IAX,AAX)
	DEALLOCATE (IF1,IF2,IF3,AF4,AF5,AF6,AF7)
	DEALLOCATE (AF8,AF9,AF10)
C	----------------------------------------
C	----------------------------------------
5000	CONTINUE




	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================



