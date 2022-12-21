C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================
	SUBROUTINE TDGSPAN(XYZ,SPAN,NS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

	DIMENSION SPAN(10,1),VR(3),VX(3),OFF(3),XYZ(NSN,1)

	ALLOCATABLE TPT(:,:),NST(:),TSTP(:,:),SLT(:)
	


C	READ(ITI,*) NP,NS
	CALL FREBUF
C	CALL FREECH
      CALL FREINT('P',NP,1)
      CALL FREINT('S',NS,1) !NUM. SPAN NODES


	ALLOCATE(TPT(3,NP),NST(NS),TSTP(2,NP),SLT(NS))

	TSTP(1:2,1:NP) = 0.0D0
	 TPT(1:3,1:NP) = 0.0D0
	DO IP = 1,NP

C	READ(ITI,*) IPT,TPT(1:3,IPT)
	CALL FREBUF
C	CALL FREECH
      CALL FREINT('N',IPT,1)	
      CALL FREREL('P',TPT(1,IPT),3)
		
	ENDDO


	TELN = 0.0D0
	DO IP = 1,NP-1	
	DO I = 1,3
      VR(I) = TPT(I,IP+1) - TPT(I,IP)
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)            !GET LENGTH OF ELEMENT HERE
	TELN = TELN + ELN
	TSTP(1,IP+1) = TELN
	ENDDO

	READ(ITI,*) NST(1:NS)  !READ SPAN NODES


	SLT(1:NS) = 0.0D0
	SELN = 0.0D0
	DO IS = 1,NS-1
	N1 = NST(IS  )
	N2 = NST(IS+1)
	VR(1:3) = 0.0D0
	DO ISC = 1,NSC
	C1 = XYZ(N1,ISC)  !GETTING HERE NODAL COORDINATE
	C2 = XYZ(N2,ISC)  !GETTING HERE NODAL COORDINATE
	VR(ISC) = C2 - C1
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)            !GET LENGTH OF ELEMENT HERE
	SELN = SELN + ELN
	SLT(IS+1) = SELN
	ENDDO	

C	----------------

	FACTOR = TELN/SELN
	DO 100 IS = 1,NS
	NN = NST(IS)
	SELN = SLT(IS)*FACTOR
	TSTP(2,1:NP) = TPT(1,1:NP)   !X
	CALL INTERPOL(SELN,VALX,TSTP,NP,IV)
	TSTP(2,1:NP) = TPT(2,1:NP)   !Y
	CALL INTERPOL(SELN,VALY,TSTP,NP,IV)
	TSTP(2,1:NP) = TPT(3,1:NP)   !Z
	CALL INTERPOL(SELN,VALZ,TSTP,NP,IV)
	FNN = FLOAT(NN)

	SPAN(1:4,IS) = [FNN,VALX,VALY,VALZ]  !NODE & COORDINATE

	VR(1:3) = 0.0D0
	DO ISC = 1,NSC
	VR(ISC) = XYZ(NN,ISC)  !GETTING HERE NODAL COORDINATE
	ENDDO
	OFF(1) = VALX - VR(1)
	OFF(2) = VALY - VR(2)
	OFF(3) = VALZ - VR(3)

	SPAN(5:7,IS) = OFF(1:3)     !OFFSET


100	CONTINUE	
	

	DEALLOCATE(TPT,NST,TSTP,SLT)



	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE TDJACK(NJACK,LTPROP,PROPM,PROPG) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN
C	==================================================================
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION LTPROP(10,NTEND)
	DIMENSION FJAK(2),FDIL(2)
	DIMENSION PROPM(1),PROPG(1)


	IF (NJACK.GT.0) READ(ITI,*)
	DO 1000 I = 1,NJACK
	READ(ITI,*) IT,IP,IJ,FA,FB,DA,DB,VLOS,ILC 

C	IP   0 = POST  1 = PRE
C	IJ   0 = STRESS RATIO  1 = STRESS VALUE   2 = FORCE VALUE
C	FA = FORCE A  FB = FORCE B
C	DA = SLIP A   DB = SLIP B
C	VLOS = LUMP LOSS IN PERCENT
C	ILC = LOAD CASE
	
	FJAK(1:2) = [FA,FB]
	FDIL(1:2) = [DA,DB]
      CALL TENLOS (LTPROP,PROPM,PROPG,NMP,NGP,LCS,
	1			 IT,ILC,IP,IJ,FJAK,FDIL,VLOS)

1000	CONTINUE



	RETURN
	END


C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================

	SUBROUTINE TDINIT(LTPROP,WA,IEL,ILC,KSTEP,IMAT,IGEO,LOPT) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 LOPT
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN
C	==================================================================
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION LTPROP(10,NTEND)
	DIMENSION WA(1)

	IF(LOPT.EQ.'CALL') IOPT = 0
	IF(LOPT.EQ.'RECD') IOPT = 1


      MTSEG = 0
	DO ITEND = 1,NTEND
	  NTSEG = LTPROP(2,ITEND)	!Number of tendon segment
	  MTSEG = MTSEG + NTSEG
	ENDDO
	
C	------------------------------------------
	SELECTCASE(IOPT)

	CASE(0)

      CALL DEFNINT('#TED',KTED,1,5)
      
	IF(ILC.LE.0) THEN
	WA(1) = 0.0D0 !FORCE AFTER JACKING
	WA(2) = 0.0D0 !JACKING FORCE
	RETURN
	ENDIF

	KEL = 0
	DO ITEND = 1,NTEND
	NSEG = LTPROP(2,ITEND)	!Number of tendon segment
	DO ISEG = 1,NSEG
	KEL = KEL + 1
	IF(KEL.EQ.IEL) THEN
	NT = ITEND
	NS = ISEG
	EXIT
	ENDIF
	ENDDO
	IF(KEL.EQ.IEL) EXIT
	ENDDO

      
	NNOD = LTPROP(1,NT)	!Number of tendon node
	NSEG = LTPROP(2,NT)	!Number of tendon segment
	IFIL = LTPROP(3,NT)	!file number
	MTREC= LTPROP(4,NT)	!last record number (just before this tendon)
	ITDYP= LTPROP(5,NT)	!ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)
	IMAT = LTPROP(6,NT)	!MAT No.
	IGEO = LTPROP(7,NT)	!GEO No.
	IFILJ= LTPROP(8,NT)	!file number (store force after jacking)


      CALL INTFILL('#TED',NT,1,1,1) !CURRENT TENDON NUMBER
      CALL INTFILL('#TED',NS,1,2,1) !CURRENT TENDON SEGMENT
	ISEGR = NNOD + NS
	IREC = MTREC+ISEGR
      CALL INTFILL('#TED',IREC ,1,3,1) !CURRENT TENDON SEGMENT RECORD
      CALL INTFILL('#TED',ITDYP,1,4,1) !ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)
      CALL INTFILL('#TED',IFIL ,1,5,1) !file number
      
	
	MPSEG = 0
	DO ITEND = 1,NTEND
	  NTSEG = LTPROP(2,ITEND)	!Number of tendon segment
	  IF(NT.EQ.ITEND) EXIT
	  MPSEG = MPSEG + NTSEG
	ENDDO
	IREC = MPSEG*LCS + NSEG*(ILC-1)
	IRECO= MTSEG*LCS + MTSEG*LCS*(KSTEP-1) + MPSEG*LCS + NSEG*(ILC-1)
		
      FJAK = 0.0D0
      FINI = 0.0D0
      READ(IFILJ,REC=IREC+NS,ERR=10) FJAK,FINI   !FORCE IMMEDIATE AFTER JACKING , JACKING FORCE WITHOUT LOSS
      
10	WA(1) = FJAK !FORCE AFTER JACKING
	WA(2) = FINI !JACKING FORCE

C	------------------------------------------
	CASE(1)

	IF(KSTEP.GT.0) THEN

	IF(ILC.LE.0) THEN
	RETURN
	ENDIF

      CALL INTFILL('#TED',NT,1,1,0) !CURRENT TENDON NUMBER
      CALL INTFILL('#TED',NS,1,2,0) !CURRENT TENDON SEGMENT

	NNOD = LTPROP(1,NT)	!Number of tendon node
	NSEG = LTPROP(2,NT)	!Number of tendon segment
	IFIL = LTPROP(3,NT)	!file number
	MTREC= LTPROP(4,NT)	!last record number (just before this tendon)
	ITDYP= LTPROP(5,NT)	!ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)
	IMAT = LTPROP(6,NT)	!MAT No.
	IGEO = LTPROP(7,NT)	!GEO No.
	IFILJ= LTPROP(8,NT)	!file number (store force after jacking)

	MPSEG = 0
	DO ITEND = 1,NTEND
	  NTSEG = LTPROP(2,ITEND)	!Number of tendon segment
	  IF(NT.EQ.ITEND) EXIT
	  MPSEG = MPSEG + NTSEG
	ENDDO
	IREC = MPSEG*LCS + NSEG*(ILC-1)
	IRECO= MTSEG*LCS + MTSEG*LCS*(KSTEP-1) + MPSEG*LCS + NSEG*(ILC-1)	

	FORC = WA(3) !FORCE DURING CALCULATION
	SIGA = WA(4) !STRESS DURING CALCULATION
      WRITE(IFILJ,REC=IRECO+NS) FORC,SIGA   !FORCE DURING CALCULATION , STRES DURING CALCULATION


	ENDIF

      CALL DELTINT('#TED')
	ENDSELECT
C	------------------------------------------


	RETURN
	END


C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================
      SUBROUTINE TEDELEM (PROPM,PROPG,S,COORD,EDIS,RE,FIN,WA)
C	-------------------------------------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION PROPM(1),PROPG(1),WA(1),RE(1),S(1)
	DIMENSION COORD(3,1),COORP(3,2),EDIS(1),FIN(1)
	DIMENSION ROT(3,3),DROT(3,2),VR(3),OFFS(3),OFFE(3)
	DIMENSION BMAT(3,NEF),BM(1,NEF),AKG(NEF,NEF)
	
	DIMENSION TOFF(3,2),COORI(3,2),EDIST(NEF)
      DIMENSION COORN1(3),COORN2(3),TOFFS1(3),TOFFS2(3),NODE1(27),NODE2(27),HF1(27),HF2(27)
      DIMENSION SH1(8),SH2(8),P(3,8)

	ALLOCATABLE FFPS(:,:)


      CALL INTFILL('#TED',NT   ,1,1,0) !CURRENT TENDON NUMBER
      CALL INTFILL('#TED',NS   ,1,2,0) !CURRENT TENDON SEGMENT
      CALL INTFILL('#TED',IREC ,1,3,0) !CURRENT TENDON SEGMENT RECORD
      CALL INTFILL('#TED',ITDYP,1,4,0) !ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)
      CALL INTFILL('#TED',IFIL ,1,5,0) !file number
      

      SELECTCASE(ITDYP)
      CASE(0)
      READ(IFIL,REC=IREC) SELN,ANGF,ANGB,
     1                      NNS1,XT1,YT1,ZT1,OFFX1,OFFY1,OFFZ1,
     2                      NNS2,XT2,YT2,ZT2,OFFX2,OFFY2,OFFZ2
      TOFF(1:3,1) = [OFFX1,OFFY1,OFFZ1]
      TOFF(1:3,2) = [OFFX2,OFFY2,OFFZ2]
      CASE(5,9,10) !SOLID
      READ(IFIL,REC=IREC) SELN,ANGF,ANGB,
     1                    IEG1,IEL1,ITYPE1,ISTYP1,NNM1,
	1                    COORN1(1:3),XT1,YT1,ZT1,NODE1(1:NNM1),TOFFS1(1:3),HF1(1:NNM1),
     2                    IEG2,IEL2,ITYPE2,ISTYP2,NNM2,
	2                    COORN2(1:3),XT2,YT2,ZT2,NODE2(1:NNM2),TOFFS2(1:3),HF2(1:NNM2)
      TOFF(1:3,1) = TOFFS1(1:3)
      TOFF(1:3,2) = TOFFS2(1:3)
      ENDSELECT
      COORI(1:3,1) = [XT1,YT1,ZT1]
      COORI(1:3,2) = [XT2,YT2,ZT2]
                 
                 
	AREA = PROPG(2)
	MP   = INT(PROPM(14))
	FORJ = WA(1)

C	NO USE
	FINI = WA(2)
	FOLD = WA(3)
	SOLD = WA(4)

	SIGO = FORJ/AREA

C     B-MATRIX      
      EDIST(1:NEF) = 0.0D0
	IF(NLOPT.GT.1) EDIST(1:NEF) = EDIS(1:NEF)
      SELECTCASE(ITDYP)
      CASE(0) !TENDON ON NODE
        CALL TDNBMAT(COORI,TOFF,EDIST,BM,TELN,'GNL')
      CASE(5,9,10) !TENDON ON ELEMENT
        CALL TDNBMATE(COORI,TOFF,HF1,HF2,NNM,NNF,NEF,EDIST,BM,TELN,'GNL')
      ENDSELECT
      
      ALLOCATE(FFPS(MP,3))

	CALL TEDYUG(PROPM(15),FFPS,SIGO,EMCO,MP)


	RE(1:NEF) = 0.0D0
	 S(1:NWS) = 0.0D0

	IF(EMCO.LE.0.0D0) GOTO 900
C     ------------------------------------------------------------
C     STRAIN TERMS 
C     ------------------------------------------------------------
	EMC = EMCO                   !JACKING STRAIN INITIAL
	DO I = 1,NEF
	EMC = EMC + BM(1,I)*EDIS(I)  !TOTAL CURRENT MECHANICAL STRAIN
	ENDDO
	
	IF(EMC.LE.0.0D0) THEN
	EMOD = 0.0D0
	SIGA = 0.0D0
	GOTO 66
	ENDIF

      DO 65 IP = 1,MP
		IF(EMC.GT.FFPS(IP,2)) GOTO 65
		EMOD   = FFPS(IP,3)               !CURRENT YOUNG
        SIGA   = EMOD * (EMC - FFPS(IP,2)) + FFPS(IP,1) !TOTAL STRESS
		GOTO 66
65	CONTINUE

66	FORC   = SIGA*AREA
	IF(FORC.LT.0.0D0) FORC = 0.0D0

      
	WA(3) = FORC  !STORE THE FORCE DURING CALCULATION HERE FOR PRINTING
	WA(4) = SIGA  !STORE THE STRES DURING CALCULATION HERE FOR PRINTING

	AKG  = AREA*TELN*(EMOD+SIGA)*MATMUL(TRANSPOSE(BM),BM)

	K = 0
	DO I = 1,NEF
	    RE(I) = RE(I) + TELN*FORC*BM(1,I)
	    DO J = I,NEF
	        K = K + 1
	        S(K) = S(K) + AKG(I,J)
	    ENDDO
	ENDDO

c      KK = 0
c      DO II = 1,NNM
c      WRITE(*,*) RE(1+KK:3+KK)
c      KK = KK + 3
c      ENDDO
c      PAUSE
      
c      IF(ITDYP.EQ.9.AND.FORJ.NE.0.0) THEN
C      WRITE(110,444) mel,RE(6 ),RE(12),RE(18),RE(24),FORC
C      WRITE(110,445)     RE(30),RE(36),RE(42),RE(48),FORC
c      WRITE(110,444) mel,RE(5 ),RE(11),RE(17),RE(23),FORC
c      WRITE(110,445)     RE(29),RE(35),RE(41),RE(47),FORC
c      ENDIF
c 444  FORMAT(I5,2X,5E16.4)
c 445  FORMAT(5X,2X,5E16.4)
      
	
900	CONTINUE

C	IF(MEL.EQ.1.OR.MEL.EQ.38.OR.MEL.EQ.75.
C	1	OR.MEL.EQ.112.OR.MEL.EQ.149.OR.MEL.EQ.186) THEN
C	WRITE(*,*) MEL,RE(1),RE(10),OFFS(3),OFFE(3)
C	PAUSE
C	ENDIF


	DEALLOCATE(FFPS)


	RETURN
	END


C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================
      SUBROUTINE TEDYUG(PROPM,FFPS,SIGO,EMCO,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=================================================================
C	GENERATE THE STRESS-STRAIN CURVE FOR CABLE MATERIAL
C	GENERATE THE STRESS-APPARENT STRAIN CURVE FOR CABLE MATERIAL
C	=================================================================

	DIMENSION FFPS(MP,3),PROPM(MP,2)


	DO I = 1,MP
	FFPS(I,1) = PROPM(I,1)     !STRESS
	FFPS(I,2) = PROPM(I,2)     !STRAIN
	ENDDO

C	STRESS-STRAIN CURVE
	S1 = 0.0
	E1 = 0.0
	DO I = 1,MP
	S2 = FFPS(I,1)
	E2 = FFPS(I,2)
	EY = (S2-S1)/(E2-E1)
	FFPS(I,3) = EY
	S1 = S2
	E1 = E2
	ENDDO


C	SIGO = JACKING STRESS
C
C-----DETERMINE MATERIAL STATE AFTER INSTALLATION-----
C
      DO 65 IP = 1,MP
		IF(SIGO.GT.FFPS(IP,1)) GOTO 65
		GOTO 66
65	CONTINUE
66	EMOD  = FFPS(IP,3)               !CURRENT YOUNG
      EMCO  = (SIGO - FFPS(IP,1)) / EMOD + FFPS(IP,2) !JACKING STRAIN



	RETURN

	END

C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================



	SUBROUTINE TENMAT (A,B,R1,R2,R3,ST1,ST2,ST3,SPANL,LOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------
C     -----------------------------
      DIMENSION A(9,9),B(9)
C
	A = 0.0
	B = 0.0
	DO I = 1,9
	A(I,I) = 1.0
	ENDDO
	
	IF(LOPT.EQ.0) THEN	
      A(1,1) = 1.0
	A(2,2) = 1.0
	A(3,2) = 1.0
	A(3,3) = 2.0*R1
	A(3,5) = -1.0

	A(4,1) = 1.0
	A(4,2) = R1
	A(4,3) = R1*R1
	A(4,4) = -1.0

	A(5,4) = 1.0
	A(5,5) = R2-R1
	A(5,6) = (R2-R1)*(R2-R1)

	A(6,4) = 1.0
	A(6,5) = SPANL-R3-R1
	A(6,6) = (SPANL-R3-R1)*(SPANL-R3-R1)
	A(6,7) = -1.0

	A(7,7) = 1.0
	A(7,8) = R3
	A(7,9) = R3*R3

	A(8,5) = 1.0
	A(8,6) = 2.0*(SPANL-R3-R1)
	A(8,8) = -1.0

	A(9,8) = 1.0
	A(9,9) = 2.0*R3

	B(1) = ST1
	B(5) = ST2
	B(7) = ST3
	ELSEIF(LOPT.EQ.1) THEN
	A(5,4) = 1.0
	A(5,5) = R2-R1
	A(5,6) = (R2-R1)*(R2-R1)

	A(6,4) = 1.0
	A(6,5) = SPANL-R3-R1
	A(6,6) = (SPANL-R3-R1)*(SPANL-R3-R1)
	A(6,7) = -1.0

	A(7,7) = 1.0
	A(7,8) = R3
	A(7,9) = R3*R3

	A(8,5) = 1.0
	A(8,6) = 2.0*(SPANL-R3-R1)
	A(8,8) = -1.0

	A(9,8) = 1.0
	A(9,9) = 2.0*R3

	B(4) = ST1
	B(5) = ST2
	B(7) = ST3
	ELSEIF(LOPT.EQ.2) THEN
	A(1,1) = 1.0
	A(2,2) = 1.0
	A(3,2) = 1.0
	A(3,3) = 2.0*R1
	A(3,5) = -1.0

	A(4,1) = 1.0
	A(4,2) = R1
	A(4,3) = R1*R1
	A(4,4) = -1.0

	A(5,4) = 1.0
	A(5,5) = R2-R1
	A(5,6) = (R2-R1)*(R2-R1)

	A(6,4) = 1.0
	A(6,5) = SPANL-R3-R1
	A(6,6) = (SPANL-R3-R1)*(SPANL-R3-R1)

	B(1) = ST1
	B(5) = ST2
	B(6) = ST3
	ELSEIF(LOPT.EQ.3) THEN
	A(5,4) = 1.0
	A(5,5) = R2-R1
	A(5,6) = (R2-R1)*(R2-R1)

	A(6,4) = 1.0
	A(6,5) = SPANL-R3-R1
	A(6,6) = (SPANL-R3-R1)*(SPANL-R3-R1)

	B(4) = ST1
	B(5) = ST2
	B(6) = ST3
	ELSEIF(LOPT.EQ.4) THEN
	A(3,1) = 1.0
	A(3,2) = R1
	A(3,3) = R1*R1

	B(1) = ST1
	IF(ST3.EQ.0.0) THEN
	B(3) = ST2
	B(4) = ST2
	B(7) = ST2
	ENDIF
	IF(ST2.EQ.0.0) THEN
	B(3) = ST3
	B(4) = ST3
	B(7) = ST3
	ENDIF

	ELSEIF(LOPT.EQ.5) THEN
	A(8,8) = 1.0
	A(8,9) = 2.0*R3
	A(9,7) = 1.0
	A(9,8) = R3
	A(9,9) = R3*R3

	B(9) = ST3
	IF(ST2.EQ.0.0) THEN
	B(1) = ST1
	B(4) = ST1
	B(7) = ST1
	ENDIF
	IF(ST1.EQ.0.0) THEN
	B(1) = ST2
	B(4) = ST2
	B(7) = ST2
	ENDIF

	ENDIF

      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE LOPREB (IDSET,DIRCOS,ELOD,NOD1,NOD2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIt INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION IDSET(1),DIRCOS(9,1),ELOD(1)
      DIMENSION COSN(6,6),IGPOS(7)

      DO 300  INF = 1,6
      IGPOS(INF) = INF
 300  CONTINUE
C     ---------------------------------------------------
C     SCAN LM FOR FIRST EQN. NO. FOR EACH NODE OF ELEMENT
C     ---------------------------------------------------
      DO 100  INX = 1,2
	IF(INX.EQ.1) INO = NOD1
	IF(INX.EQ.2) INO = NOD2

 220  ISET = IDSET(INO)
      IF (ISET.EQ.0) GOTO 100
C     -------------------------------------------
C     EXTRACT APPROPRIATE DIRECTION COSINE MATRIX
C     -------------------------------------------
      CALL COSMAT (DIRCOS(1,ISET),COSN,IGPOS,6)
      
 600  NADD = 7*(INX-1) + 1
      CALL TREVEC (ELOD(NADD),COSN,6)

 100  CONTINUE



      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE TENLOS (LTPROP,PROPM,PROPG,NMP,NGP,LCS,
	1				   KTEND,ILC,KPOS,IJ,FJAK,FDIL,VLOS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	IJ  JACKING INPUT METHOD   0 = STRESS RATIO  1 = STRESS VALUE   2 = FORCE VALUE

C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN

	DIMENSION LTPROP(10,NTEND),PROPM(NMP,1),PROPG(NGP,1)
C	-----------------------------------------------
C     COMPUTE INITIAL FORCES IN TENDON SEGMENTS     -
C	-----------------------------------------------

      LOGICAL SDON
      DIMENSION FJAK(2),FDIL(2)
	ALLOCATABLE FFPS(:,:),STFC(:),STF(:),STFR(:),FOCJ(:)
	ALLOCATABLE TEN(:),ANG(:,:)

	NNOD = LTPROP(1,KTEND)	!Number of tendon node
	NSEG = LTPROP(2,KTEND)	!Number of tendon segment
	IFIL = LTPROP(3,KTEND)	!file number
	MTREC= LTPROP(4,KTEND)	!last record number (just before this tendon)
	ITDYP= LTPROP(5,KTEND)	!ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)
	IMAT = LTPROP(6,KTEND)	!MAT No.
	IGEO = LTPROP(7,KTEND)	!GEO No.
	IFILJ= LTPROP(8,KTEND)	!file number (store force after jacking)

      IF(NSEG.LE.0) RETURN
      
	AREA = PROPG(2,IGEO)  !TENDON AREA
		
	FPYD = PROPM(6 ,IMAT)	!YIELD STRESS
	CFTR = PROPM(8 ,IMAT)	!FRICTION
	WFTR = PROPM(9 ,IMAT)	!WOBBLE
	SNPT = PROPM(14,IMAT)	!STRESS-STRAIN POINTS


	SELECTCASE(IJ)
	CASE(0)	!STRESS RATIO
	STSI = FJAK(1)*FPYD
	STSJ = FJAK(2)*FPYD
	FJAK(1) = STSI*AREA
	FJAK(2) = STSJ*AREA
	CASE(1)	!STRESS VALUE
	STSI = FJAK(1)
	STSJ = FJAK(2)
	FJAK(1) = STSI*AREA
	FJAK(2) = STSJ*AREA
	ENDSELECT


C	GET YOUNG MODULUS
	NSPT = INT(SNPT)
	ALLOCATE(FFPS(NSPT,3))
	SIGO = 0.0D0  !NO NEED TO KNOW THE JACKING STRAIN HERE
	CALL TEDYUG(PROPM(15,IMAT),FFPS,SIGO,EMCO,NSPT)
	YUNG = FFPS(3,1)

	ALLOCATE(STFC(NSEG),STF(NSEG),TEN(NSEG),ANG(2,NSEG),STFR(NSEG))
	ALLOCATE(FOCJ(NSEG))

	IDUM = NSEG
	DO ISEG = 1,NSEG
	
	ISEGR = NNOD + ISEG
      SELECTCASE(ITDYP)
      CASE(0)
      READ(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB
      CASE(5,9,10) !SOLID
      READ(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB
      ENDSELECT
      
	TEN(ISEG  ) = SELN	!SEGMENT LENGTH
	ANG(1,ISEG) = ANGF	!ANGLE FORWARD 
	ANG(2,IDUM) = ANGB	!ANGLE BACKWARD
	IDUM = IDUM - 1
	ENDDO

C	FOR PRETENSION -- NO FRICTION AND ANCHORAGE SLIP LOSS
	IF(KPOS.EQ.1) THEN
	FJAK(2) = FJAK(1)
	CFTR = 0.0D0
	WFTR = 0.0D0
	FDIL(1:2) = 0.0D0
	ENDIF



C
C-----STRESS EACH END OF THE TENDON AS REQUIRED-----
C
25	CALL RZERO (STFC,NSEG)    !INITIALIZE SEGMENT FORCE
	CALL RZERO (STF ,NSEG)    !INITIALIZE SEGMENT FORCE
C

      IDIR = 0
      DO IEND=1,2
	MSG1 = IABS(IDIR-1)
	FLST = FJAK(IEND)
	TANG = 0.0D+0
	TLEN = 0.0D+0
	TLST = 0.0D+0
	DO MS=1,NSEG
	MSEG = IABS(IDIR-MS)
	TANG = TANG + ANG(IEND,MSEG)
	TLEN = TLEN + TEN(MSEG)
	FNXT = FJAK(IEND) * DEXP(-CFTR*TANG - WFTR*TLEN)
	FAVG = (FLST + FNXT) / 2.0D+0                        !AVERAGE SEGMENT FORCE

	IF(IEND.EQ.1) STFR(MSEG) = FAVG
	IF(IEND.EQ.2) THEN
	IF(FAVG.GT.STFR(MSEG))  STFR(MSEG) = FAVG
	ENDIF

	IF(IEND.EQ.1) FOCJ(MSEG) = FJAK(IEND)
	IF(IEND.EQ.2) THEN
	IF(FJAK(IEND).GT.FOCJ(MSEG))  FOCJ(MSEG) = FJAK(IEND)
	ENDIF

	TLST = TLEN
	FLST = FNXT
	ENDDO
	IDIR = NSEG + 1
	ENDDO

C
	NEND = 0

      IDIR = 0
      DO 60 IEND=1,2

	CALL RZERO (STFC,NSEG) !!!

		MSG1 = IABS(IDIR-1)
		IF(FJAK(IEND).LE.STFC(MSG1)) GOTO 60
C
C		-----INITIALIZE FOR FRICTION AND ANCHOR SLIP LOSSES-----
C
		AREQ = FDIL(IEND) * YUNG * AREA / 2.0D+0          !AREA REQUIRED (HALF)
		APRV = 0.0D+0
		TANG = 0.0D+0
		TLEN = 0.0D+0
		TLST = 0.0D+0
		SDON = .FALSE.
		IF(AREQ.LE.0.0D+0) SDON = .TRUE.

C
C		-----FIND FORCES IN EACH TENDON SEGMENT-----
C
		FLST = FJAK(IEND)
		DO 45 MS=1,NSEG

			MSEG = IABS(IDIR-MS)
C
C			-----FRICTION LOSSES-----
C
			TANG = TANG + ANG(IEND,MSEG)
			TLEN = TLEN + TEN(MSEG)
			FNXT = FJAK(IEND) * DEXP(-CFTR*TANG - WFTR*TLEN)
			FAVG = (FLST + FNXT) / 2.0D+0                        !AVERAGE SEGMENT FORCE

C
			IF (FAVG.LT.STFC(MSEG)) GOTO 55
			STFC(MSEG) = FAVG
C
C			-----DRAW-IN LOSSES-----
C
			IF (SDON) GOTO 35
			ATRP = (FLST-FNXT) * (TLST+TLEN) / 2.0D+0            !TRAPIZOIDAL AREA (HALF)
			APRV = APRV + ATRP                                   !ACCUMULATE TRAPIZOIDAL AREA (HALF)
			IF (APRV.LT.AREQ) GOTO 35
C
C			-----STATIONARY POINT FOUND-----
C
			SDON = .TRUE.
			AREQ = AREQ - APRV + ATRP                                  !FRACMENT OF REQUIRED AREA (HALF) DA
			SLOP = TEN(MSEG) / (FLST - FNXT)                           !DL/DF
			OFST = (TLST - DSQRT(TLST*TLST + 2.0D+0*SLOP*AREQ)) / SLOP !SOLVE EQUATION 0.5*DF*(2*TLST+SLOP*DF) = DA (SOLVE FOR NEGATIVE DF)
			FSTA = FLST + OFST                                         !AJUST F = F-DF
C
C			-----ADJUST PRIOR SEGMENTS FOR DRAW-IN LOSSES-----
C
			FST2 = FSTA + FSTA
			DO 30 ME=1,MS
				MSAG = IABS(IDIR-ME)
				IF(STFC(MSAG).LE.FSTA) GOTO 30
				SKAK = DMAX1((FST2-STFC(MSAG)),0.0D0)  !ADJUST DUE TO TRIANGULAR SHAPE  ! S - 2(S-F)
				IF(SKAK.LT.STFR(MSAG)) STFR(MSAG)= SKAK !!!
30			CONTINUE

C
C			-----UPDATE TOTALS-----
C
35			TLST = TLEN
			FLST = FNXT
C

45		CONTINUE
C
C		-----STATIONARY POINT AT FAR END OF TENDON-----
C
		IF(SDON) GOTO 55

		NEND = NEND + 1
		SDON = .TRUE.
		AREQ = AREQ - APRV
		OFST = 2.0D0*AREQ / TLEN
		FSTA = STFC(MSEG) - OFST
C
		FST2 = FSTA + FSTA
		DO 50 MSAG=1,NSEG
			SKAK = DMAX1((FST2-STFC(MSAG)),0.0D0) 
			IF(SKAK.LT.STFR(MSAG)) STFR(MSAG)= SKAK !!!
50		CONTINUE

		IF(NEND.EQ.2) THEN
			DO MSAG=1,NSEG
				SKAK = DMIN1(STF(MSAG),STFC(MSAG))
				IF(SKAK.LT.STFR(MSAG)) STFR(MSAG)= SKAK !!!
			ENDDO			
		ENDIF
		STF(1:NSEG) = STFC(1:NSEG)

C
C-----THIS END STRESSING COMPLETE - CHECK FOR DRAW-IN COMPLETE-----
C
55		IF(SDON) GOTO 60
		FDIL(1) = 0.0D+0
		FDIL(2) = 0.0D+0
		WRITE(*,2000)
		GOTO 25
C
60	IDIR = NSEG + 1


C
C-----BUILD TENDON SEGMENT RECORDS-----
C
	MPSEG = 0
	DO ITEND = 1,NTEND
	  NTSEG = LTPROP(2,ITEND)	!Number of tendon segment
	  IF(ITEND.EQ.KTEND) EXIT
	  MPSEG = MPSEG + NTSEG
	ENDDO
	IREC = MPSEG*LCS + NSEG*(ILC-1)
	
      DO 70 ISEG=1,NSEG

	PLOS = VLOS/100.0                     !LUMP LOSS  DF/FI
      FINI = FOCJ(ISEG)				        !INITIAL FORCE
      FORC = STFR(ISEG) - PLOS*FINI	        !CURRENT MECHANICAL FORCE
	SIGA = FORC / AREA                    !CURRENT MECHANICAL STRESS

      WRITE(IFILJ,REC=IREC+ISEG) FORC,FINI
70	CONTINUE

C
	DEALLOCATE(FFPS)
	DEALLOCATE(STFC,STF,TEN,ANG,STFR,FOCJ)


      RETURN


2000	FORMAT(44X,'** WARNING: ERROR IN DRAW-IN LOSSES',/,
     1       44X,'   UNABLE TO COMPUTE DRAW-IN LOSSES',/,
     2       44X,'     FOR OVERLAPPING DRAW-IN ZONES.',/,
     3       44X,'   STRESS TENDON FROM ONE END ONLY',/,
     4       44X,'     OR ADJUST DRAW-IN DISTANCE.',/,
     5       44X,'   FOR THIS ANALYSIS THE EFFECTS OF',/,
     6       44X,'     ANCHOR DRAW-IN ARE *NEGLECTED*')

      END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

	SUBROUTINE TDOUTP(LTPROP,ISTEP,INM,IND) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*200 NAME

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

C	COMMON BLOCK FOR TEND SONGSAK APR2009 
	COMMON /TENDON/ NTEND,LTDN
C	==================================================================

      DIMENSION LTPROP(10,NTEND)


	NOT = 103  !TEDOUT.FLAVIA.RES


	ILC = ISTEP
	IF(IND.EQ.1) ILC = 1

	JSTEP = ISTEP
	CALL PRNLCHD(NAME,NAML,JSTEP,INM,IND)


	WRITE(NOT,2000) NAME(1:NAML),JSTEP


      MTSEG = 0
	DO ITEND = 1,NTEND
	  NTSEG = LTPROP(2,ITEND)	!Number of tendon segment
	  MTSEG = MTSEG + NTSEG
	ENDDO
	

	DO 80 ITEND = 1,NTEND
	
	NNOD = LTPROP(1,ITEND)	!Number of tendon node
	NSEG = LTPROP(2,ITEND)	!Number of tendon segment
	IFIL = LTPROP(3,ITEND)	!file number
	MTREC= LTPROP(4,ITEND)	!last record number (just before this tendon)
	ITDYP= LTPROP(5,ITEND)	!ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)
	IMAT = LTPROP(6,ITEND)	!MAT No.
	IGEO = LTPROP(7,ITEND)	!GEO No.
	IFILJ= LTPROP(8,ITEND)	!file number (store force after jacking)
	IEGT = LTPROP(9,ITEND)	!tendon group number

	MPSEG = 0
	DO JTEND = 1,NTEND
	  NTSEG = LTPROP(2,JTEND)	!Number of tendon segment
	  IF(JTEND.EQ.ITEND) EXIT
	  MPSEG = MPSEG + NTSEG
	ENDDO
	IREC = MPSEG*LCS + NSEG*(ILC-1)
	IRECO= MTSEG*LCS + MTSEG*LCS*(JSTEP-1) + MPSEG*LCS + NSEG*(ILC-1)
	
	ITEST = 0	
	DO ISEG = 1,NSEG
      READ(IFILJ,REC=IREC+ISEG) FJAK,FINI
	IF(FINI.NE.0.0D0) THEN		!IF NO JACK FORCE (MEANS NO INSTALL)
	ITEST = 1 
	EXIT
	ENDIF
	ENDDO
	IF(ITEST.EQ.0) GOTO 80


	
	TOLEN = 0.0D0
      DO 70 ISEG = 1,NSEG

	ISEGR = NNOD + ISEG
      SELECTCASE(ITDYP)
      CASE(0)
      READ(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB     !GET SEGMENT LENGTH HERE
      CASE(5,9,10) !SOLID
      READ(IFIL,REC=MTREC+ISEGR) SELN,ANGF,ANGB
      ENDSELECT

      READ(IFILJ,REC=IREC +ISEG) FJAK,FINI   !FORCE IMMEDIATE AFTER JACKING , JACKING FORCE WITHOUT LOSS
      READ(IFILJ,REC=IRECO+ISEG) FORC,SIGA   !FORCE DURING CALCULATION , STRES DURING CALCULATION

	FLOS = 0.0D0
	IF(FINI.NE.0.0D0) THEN
	FLOS=(FINI-FORC)/FINI*100.0D0
	ENDIF

	TOLEN = TOLEN + SELN
	TLN   = TOLEN - 0.5*SELN

	CALL PRNTHD(NAME,NAML,ITEND,IEGT)

	IF(ISEG.EQ.1) THEN
	WRITE(NOT,2001) NAME(1:NAML),TLN,FORC,SIGA,FLOS
	ELSE
	WRITE(NOT,2002)              TLN,FORC,SIGA,FLOS
	ENDIF

70	CONTINUE


80	CONTINUE
C
	WRITE(NOT,2003)


2000	FORMAT(/'Result "Tendon Output"',2X,
	1	    '"',A,'"',
	1		2X,I6,2X/,
	2		'ComponentNames  "R Distance" "Force" "Stress" "Loss %"'/,
	3		'Values')
2001	FORMAT(1X,'{',A,'}',2X,10E16.5)
2002	FORMAT(1X,5X,2X,10E16.5)
2003	FORMAT('End Values'/)

	RETURN
	END


C	=====================================================================
C	==========NONLINEAR TIME STEP ANALYSIS SONGSAK JUN2008===============
C	=====================================================================

	SUBROUTINE TDNBMAT(COORI,OFFI,EDIS,BM,ELN,OPTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPTN
      DIMENSION COORI(3,2),EDIS(12),OFFI(3,2)
      DIMENSION COORD(3,2),DISP(3,2),ROT(3,2),DROT(3,2)
	DIMENSION VR(3),OFFS(3,2),BMAT(3,12),BM(12)
	


C     OPTION TO UPDATE TENDON COORDINATES      
      SELECTCASE(OPTN)
      CASE('LIN') !LINEAR SMALL DEFORMATION
      COORD(1:3,1) = COORI(1:3,1)  !CALCULATE NEW TENDON COORDINATE
      COORD(1:3,2) = COORI(1:3,2)  
	 OFFS(1:3,1) =  OFFI(1:3,1)  !CALCULATE NEW TENDON OFFSET VALUE
	 OFFS(1:3,2) =  OFFI(1:3,2)  
	CASE('GNL') !GEOMETRIC NONLINEAR LARGE DEFORMATION 
C     
      DISP(1:3,1) = EDIS(1 :3 ) ! i-Node Translation
       ROT(1:3,1) = EDIS(4 :6 ) ! i-Node Rotation
      DISP(1:3,2) = EDIS(7 :9 ) ! j-Node Translation
       ROT(1:3,2) = EDIS(10:12) ! j-Node Rotation

C     DISPLACEMENT DUE TO ROTATION      
      DO I = 1,2
        CALL ROTFUNC(DROT(1,I),OFFI(1,I),ROT(1,I))
      ENDDO
      	
      COORD(1:3,1) = COORI(1:3,1) + DISP(1:3,1) + DROT(1:3,1)   !CALCULATE NEW TENDON COORDINATE
      COORD(1:3,2) = COORI(1:3,2) + DISP(1:3,2) + DROT(1:3,2)
	 OFFS(1:3,1) =  OFFI(1:3,1) + DROT(1:3,1)                 !CALCULATE NEW TENDON OFFSET VALUE
	 OFFS(1:3,2) =  OFFI(1:3,2) + DROT(1:3,2)
	ENDSELECT 

C     TENDON DIRECTIONAL VECTOR      
      VR(1:3) = COORD(1:3,2)-COORD(1:3,1)
	CALL SCALEN(VR,VR,ELN,3)  !GET LENGTH OF SEGMENT HERE
C     	
      FACL = 1.0/ELN !  1/L


C     STRAIN-DISP. MATRIX--B
	BMAT = 0.0D0
	
C     FOR i-Node
	BMAT(1,1) =  1.0D0
	BMAT(2,2) =  1.0D0
	BMAT(3,3) =  1.0D0
	BMAT(2,4) = -1.0D0*OFFS(3,1)
	BMAT(3,4) =  1.0D0*OFFS(2,1)
	BMAT(1,5) =  1.0D0*OFFS(3,1)
	BMAT(3,5) = -1.0D0*OFFS(1,1)
	BMAT(1,6) = -1.0D0*OFFS(2,1)
	BMAT(2,6) =  1.0D0*OFFS(1,1)
	BMAT(1:3,1:6) = -1.0D0*BMAT(1:3,1:6)  !Uj-Ui

C     FOR j-Node
	BMAT(1,7) =  1.0D0
	BMAT(2,8) =  1.0D0
	BMAT(3,9) =  1.0D0
	BMAT(2,10)= -1.0D0*OFFS(3,2)
	BMAT(3,10)=  1.0D0*OFFS(2,2)
	BMAT(1,11)=  1.0D0*OFFS(3,2)
	BMAT(3,11)= -1.0D0*OFFS(1,2)
	BMAT(1,12)= -1.0D0*OFFS(2,2)
	BMAT(2,12)=  1.0D0*OFFS(1,2)

C     
	DO I = 1,12
	  BM(I) = 0.0D0
	  DO J = 1,3
	      BM(I) = BM(I) + FACL*VR(J)*BMAT(J,I) !e = (Urj-Uri) / L
	  ENDDO
	ENDDO
	

	RETURN
	END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

	SUBROUTINE TDNBMATE(COORI,OFFI,H1,H2,NNM,NNF,NEFT,EDIS,BB,ELN,OPTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*3 OPTN
      DIMENSION COORI(3,2),EDIS(1),OFFI(3,2),H1(1),H2(1)
      DIMENSION COORD(3,2),DISP(3,2),ROT(3,2),DROT(3,2)
	DIMENSION VR(3),OFFS(3,2),BMAT(3,12),BM(12)
	DIMENSION BE(12,NEFT),BB(NEFT)

      NNM2 = NNM/2
      
      BE = 0.0D0
      DO I = 1,NNM2
          DO J = 1,3 !TRANSLATION DOF
            NUM = NNF-J
            BE(J,NNF*I-NUM) = H1(I)
          ENDDO
          IF(NNF.GT.3) THEN !INCLUDE ROTATION DOF
              DO J = 1,3
                NUM = NNF-J-3
                BE(J+3,NNF*I-NUM) = H1(I)
              ENDDO
          ENDIF
      ENDDO
      NUM1 = NNF*NNM2
      DO I = 1,NNM2
          DO J = 1,3 !TRANSLATION DOF
            NUM  = NNF-J
            BE(J+6,NNF*I-NUM+NUM1) = H2(I)
          ENDDO
          IF(NNF.GT.3) THEN !INCLUDE ROTATION DOF
              DO J = 1,3
                NUM = NNF-J-3
                BE(J+9,NNF*I-NUM+NUM1) = H2(I)
              ENDDO
          ENDIF
      ENDDO
      

C     OPTION TO UPDATE TENDON COORDINATES      
      SELECTCASE(OPTN)
      CASE('LIN') !LINEAR SMALL DEFORMATION
      COORD(1:3,1) = COORI(1:3,1)  !CALCULATE NEW TENDON COORDINATE
      COORD(1:3,2) = COORI(1:3,2)  
	 OFFS(1:3,1) =  OFFI(1:3,1)  !CALCULATE NEW TENDON OFFSET VALUE
	 OFFS(1:3,2) =  OFFI(1:3,2)  
	CASE('GNL') !GEOMETRIC NONLINEAR LARGE DEFORMATION   	
C     
      DISP(1:3,1) = 0.0 ! i-Node Translation
       ROT(1:3,1) = 0.0 ! i-Node Rotation
      DISP(1:3,2) = 0.0 ! j-Node Translation
       ROT(1:3,2) = 0.0 ! j-Node Rotation
      DO I = 1,NEFT
          DO J = 1,3
            DISP(J,1) = DISP(J,1) + BE(J+0,I)*EDIS(I)
            DISP(J,2) = DISP(J,2) + BE(J+6,I)*EDIS(I)
             ROT(J,1) =  ROT(J,1) + BE(J+3,I)*EDIS(I)
             ROT(J,2) =  ROT(J,2) + BE(J+9,I)*EDIS(I)
          ENDDO  
      ENDDO 

C     DISPLACEMENT DUE TO ROTATION      
      DO I = 1,2
        CALL ROTFUNC(DROT(1,I),OFFI(1,I),ROT(1,I))
      ENDDO
      	
      COORD(1:3,1) = COORI(1:3,1) + DISP(1:3,1) + DROT(1:3,1)   !CALCULATE NEW TENDON COORDINATE
      COORD(1:3,2) = COORI(1:3,2) + DISP(1:3,2) + DROT(1:3,2)
	 OFFS(1:3,1) =  OFFI(1:3,1) + DROT(1:3,1)                 !CALCULATE NEW TENDON OFFSET VALUE
	 OFFS(1:3,2) =  OFFI(1:3,2) + DROT(1:3,2)
	ENDSELECT 
	
C     TENDON DIRECTIONAL VECTOR      
      VR(1:3) = COORD(1:3,2)-COORD(1:3,1)
	CALL SCALEN(VR,VR,ELN,3)  !GET LENGTH OF SEGMENT HERE
C     	
      FACL = 1.0/ELN !  1/L


C     STRAIN-DISP. MATRIX--B
	BMAT = 0.0D0
	
C     FOR i-Node
	BMAT(1,1) =  1.0D0
	BMAT(2,2) =  1.0D0
	BMAT(3,3) =  1.0D0
	BMAT(2,4) = -1.0D0*OFFS(3,1)
	BMAT(3,4) =  1.0D0*OFFS(2,1)
	BMAT(1,5) =  1.0D0*OFFS(3,1)
	BMAT(3,5) = -1.0D0*OFFS(1,1)
	BMAT(1,6) = -1.0D0*OFFS(2,1)
	BMAT(2,6) =  1.0D0*OFFS(1,1)
	BMAT(1:3,1:6) = -1.0D0*BMAT(1:3,1:6)  !Uj-Ui

C     FOR j-Node
	BMAT(1,7) =  1.0D0
	BMAT(2,8) =  1.0D0
	BMAT(3,9) =  1.0D0
	BMAT(2,10)= -1.0D0*OFFS(3,2)
	BMAT(3,10)=  1.0D0*OFFS(2,2)
	BMAT(1,11)=  1.0D0*OFFS(3,2)
	BMAT(3,11)= -1.0D0*OFFS(1,2)
	BMAT(1,12)= -1.0D0*OFFS(2,2)
	BMAT(2,12)=  1.0D0*OFFS(1,2)

C     
	DO I = 1,12
	  BM(I) = 0.0D0
	  DO J = 1,3
	      BM(I) = BM(I) + FACL*VR(J)*BMAT(J,I) !e = (Urj-Uri) / L
	  ENDDO
	ENDDO
	
	DO I = 1,NEFT
	    BB(I) = 0.0D0
	    DO J = 1,12
	      BB(I) = BB(I) + BM(J)*BE(J,I)
	    ENDDO
	ENDDO

	RETURN
	END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

	SUBROUTINE ROTFUNC(DROT,TOFF,ROT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION DROT(3),TOFF(3),ROT(3),ROTM(3,3)
    
	ROTM = 0.0D0
	ROTM(2,1) = -1.0D0*TOFF(3)
	ROTM(3,1) =  1.0D0*TOFF(2)
	ROTM(1,2) =  1.0D0*TOFF(3)
	ROTM(3,2) = -1.0D0*TOFF(1)
	ROTM(1,3) = -1.0D0*TOFF(2)
	ROTM(2,3) =  1.0D0*TOFF(1)

	DO I = 1,3
	    DROT(I) = 0.0D0
	    DO J = 1,3
	      DROT(I) = DROT(I) + ROTM(I,J)*ROT(J)
	    ENDDO
	ENDDO

	 
	RETURN
	END 

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================


C    ********************************************************************************
C    *															                    *	
C    *    SUBROUTINE : SHAPE FUNCTION						                        *
C    *													    		                *	
C    *    ======================================================================    *	
C    *	  Purpose:	To create SHAPE FUNCTION for each element                       *	
C    *													                            *	
C    *															                    *
C    ********************************************************************************        
      
      SUBROUTINE SHAPE8(xi,eta,zeta,SH,P)        
     
C	-----------------------------------------------------------------------     
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
      DIMENSION SH(1,8),P(3,8)    
C	----------------------------------------------------------------------- 

      ! Shape Function
      SH(1,1) = 0.125*(1.0 + xi)*(1.0 + eta)*(1.0 + zeta)
      SH(1,2) = 0.125*(1.0 - xi)*(1.0 + eta)*(1.0 + zeta)
      SH(1,3) = 0.125*(1.0 - xi)*(1.0 - eta)*(1.0 + zeta)
      SH(1,4) = 0.125*(1.0 + xi)*(1.0 - eta)*(1.0 + zeta)
      SH(1,5) = 0.125*(1.0 + xi)*(1.0 + eta)*(1.0 - zeta)
      SH(1,6) = 0.125*(1.0 - xi)*(1.0 + eta)*(1.0 - zeta)
      SH(1,7) = 0.125*(1.0 - xi)*(1.0 - eta)*(1.0 - zeta)
      SH(1,8) = 0.125*(1.0 + xi)*(1.0 - eta)*(1.0 - zeta)
            
      ! Diff with xi                     
      P(1,1) = +0.125*(1.0 + eta)*(1.0 + zeta)
      P(1,2) = -0.125*(1.0 + eta)*(1.0 + zeta)
      P(1,3) = -0.125*(1.0 - eta)*(1.0 + zeta)
      P(1,4) = +0.125*(1.0 - eta)*(1.0 + zeta)
      P(1,5) = +0.125*(1.0 + eta)*(1.0 - zeta)
      P(1,6) = -0.125*(1.0 + eta)*(1.0 - zeta)
      P(1,7) = -0.125*(1.0 - eta)*(1.0 - zeta)
      P(1,8) = +0.125*(1.0 - eta)*(1.0 - zeta)
            
      ! Diff with eta
      P(2,1) = +0.125*(1.0 + xi)*(1.0 + zeta)
      P(2,2) = +0.125*(1.0 - xi)*(1.0 + zeta)
      P(2,3) = -0.125*(1.0 - xi)*(1.0 + zeta)
      P(2,4) = -0.125*(1.0 + xi)*(1.0 + zeta)
      P(2,5) = +0.125*(1.0 + xi)*(1.0 - zeta)
      P(2,6) = +0.125*(1.0 - xi)*(1.0 - zeta)
      P(2,7) = -0.125*(1.0 - xi)*(1.0 - zeta)
      P(2,8) = -0.125*(1.0 + xi)*(1.0 - zeta)
            
      ! Diff with zeta
      P(3,1) = +0.125*(1.0 + xi)*(1.0 + eta)
      P(3,2) = +0.125*(1.0 - xi)*(1.0 + eta)
      P(3,3) = +0.125*(1.0 - xi)*(1.0 - eta)
      P(3,4) = +0.125*(1.0 + xi)*(1.0 - eta)
      P(3,5) = -0.125*(1.0 + xi)*(1.0 + eta)
      P(3,6) = -0.125*(1.0 - xi)*(1.0 + eta)
      P(3,7) = -0.125*(1.0 - xi)*(1.0 - eta)
      P(3,8) = -0.125*(1.0 + xi)*(1.0 - eta)      
      
      
C     ------------------------------------------------------------------   
      RETURN
	END       
      
C    ********************************************************************************
C    *					  									                        *	
C    *                      THE END OF THE SUBROUTINE PROGRAM                       *
C    *															                    *
C    ********************************************************************************

