C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE TRVLOOP
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	FORCE AND STIFFNESS
	ALLOCATABLE RE(:),SE(:),STLD(:),EXDI(:)
	ALLOCATABLE KEQ(:),KRQ(:)



      CALL LOCATN('-GSF',KGSF,NUM,NWK ,2) 
	CALL LOCATN('-LMA',KLMA,NUM,NEQ1,1)  


	CALL INTFILL('%NUB',NSN,1,1,0)
	CALL INTFILL('%NUB',NSF,1,2,0)
	CALL INTFILL('%NUB',NSC,1,4,0)


C	SOLUTION CONTROL ARRAY
      CALL INTFILL('%SOL',LOPT1,1,1,0)
      CALL INTFILL('%SOL',LOPT2,1,2,0)
      CALL INTFILL('%SOL',LOPT3,1,3,0)
      CALL INTFILL('%SOL',LOPT4,1,4,0)
      CALL INTFILL('%SOL',LOPT5,1,5,0)

C	------------------------------------
C	LOOP OVER LINK ELEMENT
      CALL LOCATN ('+TRV',KTRV,LWRK,NELE,2) !CALLING NELE
C	------------------------------------
	DO 5000 IEL = 1,NELE


	CALL RELFILL('+TRV',FNM ,1,IEL,0)  !NUM NODE
	CALL RELFILL('+TRV',FNF ,4,IEL,0)
	NNM   = INT(FNM) !NUM NODE
	NNF   = INT(FNF)
	NEF   = NNF*NNM
	NEK = (NEF*NEF+NEF)/2

	IF(NNM.LE.0) GOTO 5000

	ALLOCATE(RE(NEF),SE(NEK),STLD(NEF),EXDI(NEF))
	ALLOCATE(KEQ(NEF),KRQ(NEF))

	SE(1:NEK)   = 0.0D0
	RE(1:NEF)   = 0.0D0
	EXDI(1:NEF) = 0.0D0
	STLD(1:NEF) = 0.0D0
	KEQ(1:NEF)  = 0
	KRQ(1:NEF)  = 0


C	------------------------------------
C	IF(LOPT2.GT.0) 
C		NOTHING FOR TRAVELER
C	ENDIF
C	------------------------------------

C	------------------------------------
C	IF(LOPT1.GT.0) THEN !UPDAT MECHANICAL
C		NOTHING FOR TRAVELER
C	ENDIF
C	------------------------------------


C	LSTF  0 = NO COMPUTE OF STIFFNESS   1 = COMPUTE STIFFNESS
C	LFOC  0 = NO COMPUTE OF ELM FORCE   1 = COMPUTE ELM FORCE
	LSTF = 0
	IF(LOPT4.GT.0.OR.LOPT5.GT.0) LSTF = 1
C	------------------------------------
	IF(LOPT3.EQ.1) THEN
C		NOTHING FOR TRAVELER NON MECHANICAL
	ENDIF
C	------------------------------------
	IF(LOPT3.EQ.2) THEN
C	--------------------------
	LFOC = 2
	CALL TRVCAL(IEL,SE,RE,KRQ,KEQ,LSTF,LFOC)
	DO IEF = 1,NEF
	IEQ  = KRQ(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-EIS',RE(IEF),IEQ,1,2)    !STORE MECHANICAL FORCE
	ENDDO
C	------------------------------------
	ENDIF
C	------------------------------------
	IF(LOPT3.EQ.3) THEN
C	--------------------------
C		NOTHING FOR TRAVELER NON MECHANICAL
C	--------------------------
	LFOC = 2
	CALL TRVCAL(IEL,SE,RE,KRQ,KEQ,LSTF,LFOC)
	DO IEF = 1,NEF
	IEQ  = KRQ(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-EIS',RE(IEF),IEQ,1,2)    !STORE MECHANICAL FORCE
	ENDDO
C	------------------------------------
	ENDIF
C	------------------------------------
	IF(LOPT3.EQ.4) THEN
C	--------------------------
	LFOC = 2
	CALL TRVCAL(IEL,SE,RE,KRQ,KEQ,LSTF,LFOC)
	DO IEF = 1,NEF
	IEQ  = KRQ(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-EIS',RE(IEF),IEQ,1,2)    !STORE MECHANICAL FORCE
	ENDDO
C	------------------------------------
	ENDIF
C	------------------------------------


	IF(LOPT3.EQ.0.AND.LSTF.GT.0) THEN
	LFOC = 0
	CALL TRVCAL(IEL,SE,RE,KRQ,KEQ,LSTF,LFOC) !STIFFNESS COMPUTATION HERE IF LOPT3 EQ 0
	ENDIF


C	ASSEMBLE OF STIFFNESS
	IF(LOPT4.GT.0) THEN
	CALL STIFLAX(KRQ,SE,NEF,NNM)  !TRANSFORMATION DUE TO LOCAL SUPPORT
	CALL ASBSTF(SE,KEQ,NEF) 
	ENDIF

C	------------------------------------
C	FORCE FROM INCREMENTAL PRESCRIBED DISP.
	IF(LOPT5.GT.0) THEN
	DO IEF = 1,NEF
	IEQ  = KRQ(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('$XDI',EXDI(IEF),IEQ,1,0)  !GETTING HERE INCREMENTAL NODAL PRESCRIBED DISPLACEMENT
	ENDDO
	CALL MULTRI(SE,EXDI,STLD,NEF)
	DO IEF = 1,NEF
	IEQ  = KRQ(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-STL',STLD(IEF),IEQ,1,2)  !STORE FORCE FROM INCREMENTAL NODAL PRESCRIBED DISPLACEMENT
	ENDDO
	ENDIF
C	------------------------------------

	DEALLOCATE(RE,SE,STLD,EXDI)
	DEALLOCATE(KEQ,KRQ)

C	------------------------------------
5000	CONTINUE
C	------------------------------------



	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE TRAVCSA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	FILE FOR PRINTING  NOT=ISO
	CALL INTFILL('%IOL',NOT,1,3,0)

	CALL FREBUF
	CALL FREECH
      CALL FREINT('S',NTRS,1)
      CALL FREINT('N',NTRV,1)
      CALL FREINT('M',MAXM,1) !MAX NODE

	IF(NTRV.LE.0) RETURN
	

C	INPUT OF LINK ELEMENT
	LLEN = 20
      CALL DEFNREL('+TRS',KTRS,LLEN,NTRS) 
	CALL RELZERO('+TRS')

	DO 1000 ITRS = 1,NTRS

	CALL FREBUF
	CALL FREECH
      CALL FREINT('N',N,1)  !NO.
      CALL FREREL('A',AREA,1)
      CALL FREREL('B',QS  ,1)
      CALL FREREL('C',QT  ,1)
      CALL FREREL('D',SIS ,1)
      CALL FREREL('E',SIT ,1)
      CALL FREREL('F',SIST,1)
      CALL FREREL('G',TOR ,1)
      CALL FREREL('H',WEIG,1)
      CALL FREREL('I',YUNG,1)
      CALL FREREL('J',GUNG,1)
      CALL FREREL('R',RANG,1)

	CALL RELFILL('+TRS',AREA,1 ,N,1)
	CALL RELFILL('+TRS',QS  ,2 ,N,1)
	CALL RELFILL('+TRS',QT  ,3 ,N,1)
	CALL RELFILL('+TRS',SIS ,4 ,N,1)
	CALL RELFILL('+TRS',SIT ,5 ,N,1)
	CALL RELFILL('+TRS',SIST,6 ,N,1)
	CALL RELFILL('+TRS',TOR ,7 ,N,1)
	CALL RELFILL('+TRS',WEIG,8 ,N,1)
	CALL RELFILL('+TRS',YUNG,9 ,N,1)
	CALL RELFILL('+TRS',GUNG,10,N,1)
	CALL RELFILL('+TRS',RANG,11,N,1)

1000	CONTINUE	



	NNF  = 6
	LINK = 123456
	MAXF = NNF*MAXM


C	DATA + NOD NO + ORI EQUA + UPD EQUA + REA EQUA + ERE DISP + FIXEND
	LWRK = 10+MAXF+MAXF+MAXF+MAXF+MAXF+MAXF+MAXF
      CALL DEFNREL('+TRV',KTRV,LWRK,NTRV) 
	CALL RELZERO('+TRV')

	DO 2000 ITRV = 1,NTRV

C	READ WHEN INSTALL
C	CALL RELFILL('+TRV',  ,1,N,1)  !NUM NOD  
C	CALL RELFILL('+TRV',  ,2,N,1)  !MAT SET  
C	CALL RELFILL('+TRV',  ,3,N,1)  !EXISTING FLAG

	FNNF = FLOAT(NNF)
	FLNK = FLOAT(LINK)
	CALL RELFILL('+TRV',FNNF,4,ITRV,1)  !NNF
	CALL RELFILL('+TRV',FLNK,5,ITRV,1)  !LINK


2000	CONTINUE	



	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE TRVEQNR(MHT,IEMQ,IEMR)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      DIMENSION MHT(1)
      DIMENSION IEMQ(1),IEMR(1)


	ALLOCATABLE LEQN(:)


      CALL LOCATN ('+TRV',KTRV,LWRK,NELE,2) 
C	-----------------------------------------
C	LOOP OVER TENDON ELEMENT
C	-----------------------------------------
	DO 1600 IEL = 1,NELE

	CALL RELFILL('+TRV',FNM ,1,IEL,0)  !NUM NODE
C	CALL RELFILL('+TRV',FMT ,2,IEL,0)  !MATERIAL SET
	CALL RELFILL('+TRV',XIST,3,IEL,0)
	CALL RELFILL('+TRV',FNF ,4,IEL,0)
	CALL RELFILL('+TRV',FNK ,5,IEL,0)
	NNM   = INT(FNM) !NUM NODE
	NNF   = INT(FNF)
	NEF   = NNF*NNM
	LINK  = INT(FNK)
	IXIST = INT(XIST)


	IF (IXIST.EQ. 0) GOTO 1600 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED

	ALLOCATE(LEQN(NEF))


C	-----------------------------------------
C	LOOP OVER ELEMENT DOF
C	-----------------------------------------
	LEQN(1:NEF) = 0
	DO 1000 IEF = 1,NEF
	NFI = 10 + NNM + IEF 
	NUQ = 10 + NNM + NEF + IEF                  !UPDATE EQUATION NUMBER
	NUR = 10 + NNM + NEF + NEF + IEF            !UPDATE REACTION NUMBER

	CALL RELFILL('+TRV',FNFI,NFI,IEL,0)
	J = INT(FNFI)  !ORIGINAL EQUA

	FNUQ  = 0.0D0
	FNUR  = 0.0D0
	IF(J.GT.0) THEN
	FNUQ  = FLOAT(IEMQ(J))
	FNUR  = FLOAT(IEMR(J))
	CALL RELFILL('+TRV',FNUQ,NUQ,IEL,1)
	CALL RELFILL('+TRV',FNUR,NUR,IEL,1)
	LEQN(IEF) = IEMQ(J)   !TO CALCULATE COLUMN HEIGHT
	ENDIF

C	-----------------------------------------
1000	CONTINUE
C	-----------------------------------------


C     ------------------------------------
C     UPDATE COLUMN HEIGHTS [MHT]
C     ------------------------------------
      MEQ = 0
      DO 390  IEF=1,NEF
      IEQ = LEQN(IEF)
      IF (IEQ)     310,390,310
310   IF (MEQ.EQ.0) MEQ = IEQ
      IF (MEQ-IEQ) 320,390,390
320	MEQ = IEQ
390	CONTINUE

      DO 400  IEF=1,NEF
      IEQ = LEQN(IEF)
      IF (IEQ.EQ.0) GOTO 400
      KHT = MEQ-IEQ
      IF (KHT.GT.MHT(IEQ)) MHT(IEQ) = KHT
400	CONTINUE




	DEALLOCATE(LEQN)
C	-----------------------------------------
1600	CONTINUE
C	-----------------------------------------


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE TRVCAL(IEL,S,RE,KRQ,KEQ,LSF,LFC)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C	KEQ = CURRENT EQUATION NUMBER
C	KRQ = ORIGINAL EQUATION NUMBER
C     ------------------------------------------------------------- 

	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)

      DIMENSION S(1),RE(1),KEQ(1),KRQ(1)



	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION TT(12,12)
	DIMENSION BPG(10),BWG(10),BBX(12)
	DIMENSION NPM(10),NPI(10),EDISL(12)


	DIMENSION EPS(8),DR(6,6),SIGR(8),RL(12)
	DIMENSION AKG(12,12),AKL(12,12),BMSG(6,12)


	ALLOCATABLE IE(:),EIJ(:),LEQN(:),LEQ(:)
	ALLOCATABLE COORD(:,:),EDIS(:)



	CALL RELFILL('+TRV',FNM ,1,IEL,0)  !NUM NODE
	CALL RELFILL('+TRV',FMT ,2,IEL,0)  !MATERIAL SET
	CALL RELFILL('+TRV',XIST,3,IEL,0)
	CALL RELFILL('+TRV',FNF ,4,IEL,0)
	CALL RELFILL('+TRV',FNK ,5,IEL,0)
	NNM   = INT(FNM) !NUM NODE
	MATN  = INT(FMT)
	NNF   = INT(FNF)
	NEF   = NNF*NNM
	LINK  = INT(FNK)
	IXIST = INT(XIST)


	KEQ(1:NEF) = 0
	KRQ(1:NEF) = 0


	IF (IXIST.EQ. 0) GOTO 1600 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED


	ALLOCATE(IE(NEF),EIJ(NNM),LEQN(NEF),LEQ(NEF))
	ALLOCATE(COORD(3,NNM),EDIS(NEF))


C	CALLING NODE NUM
	DO INM = 1,NNM
	NUM = 10+INM
	CALL RELFILL('+TRV',EIJ(INM),NUM,IEL,0)
	ENDDO



C	CONSTRUCT ELEMENT LOCAL EQUATION NUMBER
	DO IEF = 1,NEF

	LEQN(IEF) = IEF

	NUM = 10 + NNM + IEF
	CALL RELFILL('+TRV',FQ,NUM,IEL,0) !ORIGINAL EQUATION NUMBER
	KRQ(IEF) = INT(FQ)

	NUM = 10 + NNM + NEF + IEF
	CALL RELFILL('+TRV',FQ,NUM,IEL,0) !CURRENT EQUATION NUMBER
	KEQ(IEF) = INT(FQ)

	ENDDO

C	--------------------------
C	EXTRACT NODAL COORDINATE
	CALL INTFILL('%NUB',NSC,1,4,0)
	DO INM = 1,NNM
	NOD = INT(EIJ(INM))
	COORD(1:3,INM) = 0.0D0
	DO ISC = 1,NSC
	CALL RELFILL('$XYZ',COORD(ISC,INM),ISC,NOD,0)  !GETTING HERE NODAL COORDINATE
	ENDDO
	ENDDO

C	EXTRACT NODAL DISPLACEMENT
	EDIS(1:NEF) = 0.0D0
	DO IEF = 1,NEF
	IEQ = KRQ(IEF)
	CALL RELFILL('$XDS',EDIS(IEF),IEQ,1,0)    !GETTING HERE NODAL DISPLACEMENT

C	DO INM = 1,NNM
C	NN = 6*(INM-1)
C	DO ISC = 1,NSC
C	COORD(ISC,INM) = COORD(ISC,INM) + EDIS(ISC+NN)
C	ENDDO
C	ENDDO

C	ERECTED DISP.
	NUM = 10 + NNM + NEF + NEF + NEF + IEF
	CALL RELFILL('+TRV',EDISO,NUM,IEL,0)      !GETTING HERE NODAL DISPLACEMENT
	EDIS(IEF) = EDIS(IEF) - EDISO

	ENDDO
C	--------------------------

	CALL RELFILL('+TRS',AREA,1 ,MATN,0)
	CALL RELFILL('+TRS',QS  ,2 ,MATN,0)
	CALL RELFILL('+TRS',QT  ,3 ,MATN,0)
	CALL RELFILL('+TRS',SIS ,4 ,MATN,0)
	CALL RELFILL('+TRS',SIT ,5 ,MATN,0)
	CALL RELFILL('+TRS',SIST,6 ,MATN,0)
	CALL RELFILL('+TRS',TOR ,7 ,MATN,0)
	CALL RELFILL('+TRS',WEIG,8 ,MATN,0)
	CALL RELFILL('+TRS',YUNG,9 ,MATN,0)
	CALL RELFILL('+TRS',GUNG,10,MATN,0)
	CALL RELFILL('+TRS',RANG,11,MATN,0)

	AE   = AREA*YUNG
	QSE  = QS  *YUNG
	QTE  = QT  *YUNG
	SIE  = SIS *YUNG
	TIE  = SIT *YUNG
	SITE = SIST*YUNG
	PJL  = TOR *GUNG 

C	-------------------------------------
C	LOOP OVER PART
	DO 500 IPRT = 1,NNM-1

	RL(1:12) = 0.0D0
	AKG(1:12,1:12) = 0.0D0

	N1 = INT(EIJ(IPRT  ))
	N2 = INT(EIJ(IPRT+1))
	
	LEQ(1:NEF) = 0
	M1 = 1 +NNF*(IPRT-1)
	M2 = 12+NNF*(IPRT-1)
	LEQ(M1:M2) = LEQN(M1:M2)

	VR = 0.0D0
C	LOCAL VECTOR
	DO I = 1,3
	VR(I) = COORD(I,IPRT+1) - COORD(I,IPRT) 
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)

	IF(ELN.EQ.0) VR(1:3) = [1.0D0,0.0D0,0.0D0]
	CALL FMVEVR (VR,VS,VT)
	CALL ROMBAC (VR,VS,VT,RANG)

C	DETERMINE GAUSS POINT POSITION AND WEIGHT
	NGR = 2
	DO IGR = 1,NGR
	RI = GAUSP(IGR,NGR)
	RW = GAUSW(IGR,NGR)
	BPG(IGR) = 0.5*ELN*(1.0 + RI)
	BWG(IGR) = 0.5*ELN*RW
	ENDDO


	TT(1:12,1:12) = 0.0D0
      TT(1,1)   = VR(1)
      TT(1,2)   = VS(1)
      TT(1,3)   = VT(1)
      TT(2,1)   = VR(2)
      TT(2,2)   = VS(2)
      TT(2,3)   = VT(2)
      TT(3,1)   = VR(3)
      TT(3,2)   = VS(3)
      TT(3,3)   = VT(3)
      TT(4,4)   = VR(1)
      TT(4,5)   = VS(1)
      TT(4,6)   = VT(1)
      TT(5,4)   = VR(2)
      TT(5,5)   = VS(2)
      TT(5,6)   = VT(2)
      TT(6,4)   = VR(3)
      TT(6,5)   = VS(3)
      TT(6,6)   = VT(3)

      TT(7,7)   = VR(1)
      TT(7,8)   = VS(1)
      TT(7,9)   = VT(1)
      TT(8,7)   = VR(2)
      TT(8,8)   = VS(2)
      TT(8,9)   = VT(2)
      TT(9,7)   = VR(3)
      TT(9,8)   = VS(3)
      TT(9,9)   = VT(3)
      TT(10,10) = VR(1)
      TT(10,11) = VS(1)
      TT(10,12) = VT(1)
      TT(11,10) = VR(2)
      TT(11,11) = VS(2)
      TT(11,12) = VT(2)
      TT(12,10) = VR(3)
      TT(12,11) = VS(3)
      TT(12,12) = VT(3)

	
	EDISL(1:12) = 0.0D0
	II = 0
	DO IEF = 1,NEF
	IF(LEQ(IEF).NE.0) II = II + 1
	JJ = 0
	DO JEF = 1,NEF
	IF(LEQ(JEF).NE.0) JJ = JJ + 1
	IF(LEQ(IEF).NE.0.AND.LEQ(JEF).NE.0) THEN
	EDISL(II) = EDISL(II) + TT(JJ,II) * EDIS(JEF)
	ENDIF
	ENDDO
	ENDDO


C	-------------------------------------
C	LOOP OVER GAUSS
	DO 400 IGR = 1,NGR


C	INITILIZE ARRAY
	 EPS(1:8) = 0.0D0				!STRAIN
	SIGR(1:8) = 0.0D0               !RESULTANT FORCE

C     GAUSS LOCATION ALONG ELEMENT AXIS
	BXD = BPG(IGR)



C	BMATRIX
	FBM1  =  1.0/ELN                             
      FBM2  = -6.0*(1.0-2.0*BXD/ELN)/ELN**2.0    
	FBM3  = -6.0*(1.0-2.0*BXD/ELN)/ELN**2.0
	FBM4  =  1.0/ELN 
	FBM5  = -2.0*(2.0-3.0*BXD/ELN)/ELN
	FBM6  = -2.0*(2.0-3.0*BXD/ELN)/ELN
	FBM7  =  1.0/ELN
	FBM8  =  6.0*(1.0-2.0*BXD/ELN)/ELN**2.0 
	FBM9  =  6.0*(1.0-2.0*BXD/ELN)/ELN**2.0 
	FBM10 =  1.0/ELN
	FBM11 = -2.0*(1.0-3.0*BXD/ELN)/ELN
	FBM12 = -2.0*(1.0-3.0*BXD/ELN)/ELN
		            
	BBX(1) = -FBM1                               
      BBX(2) = +FBM2                
      BBX(3) = -FBM3        
      BBX(6) = -FBM4      !TORSION                                                          
      BBX(4) = +FBM5                   
      BBX(5) = +FBM6               
      BBX(7) = +FBM7                                
      BBX(8) = +FBM8               
      BBX(9 )= -FBM9              
      BBX(12)= +FBM10     !TORSION                                                                   
      BBX(10)= +FBM11
      BBX(11)= +FBM12

C     ------------------------------------------------------------
C     STRAIN TERMS MEMBRANE - BENDING - TORSION
C     ------------------------------------------------------------
	EPS(1)= BBX(1)*EDISL(1)+BBX(7)*EDISL(7)
      EPS(4)= BBX(3)*EDISL(3)+BBX(4)*EDISL(5)+BBX(9)*EDISL(9)+
     +		BBX(10)*EDISL(11)
      EPS(5)= BBX(2)*EDISL(2)+BBX(5)*EDISL(6)+BBX(8)*EDISL(8)+
     +		BBX(11)*EDISL(12)
	EPS(6)= BBX(6)*EDISL(4)+BBX(12)*EDISL(10)


	SIGR(1) = AE*EPS(1) + QSE*EPS(4) - QTE*EPS(5)   !AXIAL FORCE
	SIGR(4) = QSE*EPS(1)+ SIE*EPS(4) - SITE*EPS(5)  !MOMENT ABOUT S-AXIS (MY FOR XY PLANE)
	SIGR(5) =-QTE*EPS(1)- SITE*EPS(4)+ TIE*EPS(5)   !MOMENT ABOUT T-AXIS (MZ FOR XY PLANE)
	SIGR(6) = PJL*EPS(6)                            !TORSIONAL MOMENT


	DR(1:6,1:6) = 0.0D0
	DR(1,1) = +AE     !AE
	DR(4,1) = +QSE    !QsE
	DR(5,1) = -QTE    !QtE
	DR(4,4) = +SIE    !IsE
	DR(5,4) = -SITE   !IstE
	DR(5,5) = +TIE    !ItE
	DR(1,4) = +QSE    !QsE
	DR(1,5) = -QTE    !QtE
	DR(4,5) = -SITE   !IstE
	DR(6,6) = +PJL    !GJ


	BMSG(1:6,1:12) = 0.0D0
	BMSG(1,1)  =  BBX(1)
	BMSG(1,7)  =  BBX(7)
	BMSG(4,3)  =  BBX(3)
	BMSG(4,5)  =  BBX(4)
	BMSG(4,9 ) =  BBX(9 )
	BMSG(4,11) =  BBX(10)
	BMSG(5,2)  =  BBX(2)
	BMSG(5,6)  =  BBX(5)
	BMSG(5,8)  =  BBX(8)
	BMSG(5,12) =  BBX(11)
	BMSG(6,4)  =  BBX(6)
	BMSG(6,10) =  BBX(12)



	AKL = MATMUL(TRANSPOSE(BMSG),MATMUL(DR,BMSG))

	DO I = 1,12
	DO J = 1,12
	AKG(I,J) = AKG(I,J) + BWG(IGR)*AKL(I,J)
	ENDDO
	ENDDO

C	MEMBRANE
	RL(1)  = RL(1)  + BWG(IGR)*SIGR(1)*BBX(1) 
	RL(7)  = RL(7)  + BWG(IGR)*SIGR(1)*BBX(7) 
C	BENDING
      RL(2)  = RL(2)  + BWG(IGR)*SIGR(5)*BBX(2) 
      RL(3)  = RL(3)  + BWG(IGR)*SIGR(4)*BBX(3) 
      RL(5)  = RL(5)  + BWG(IGR)*SIGR(4)*BBX(4) 
      RL(6)  = RL(6)  + BWG(IGR)*SIGR(5)*BBX(5) 
      RL(8)  = RL(8)  + BWG(IGR)*SIGR(5)*BBX(8) 
      RL(9)  = RL(9) + BWG(IGR)*SIGR(4)*BBX(9) 
      RL(11) = RL(11) + BWG(IGR)*SIGR(4)*BBX(10)
      RL(12) = RL(12) + BWG(IGR)*SIGR(5)*BBX(11)
C	TORSION
	RL(4)  = RL(4)  + BWG(IGR)*SIGR(6)*BBX(6)
	RL(10) = RL(10) + BWG(IGR)*SIGR(6)*BBX(12) 


C	----------------------------------------
400	CONTINUE  !END GAUSS
C	----------------------------------------



C	TRANSFORMATION TO GLOBAL SYSTEM
	RL   = MATMUL(TT,RL)
	AKG  = MATMUL(MATMUL(TT,AKG),TRANSPOSE(TT))


	IF(LFC.GT.0) THEN
	K = 0
	II = 0
	DO I = 1,NEF
	K = K + 1
	IF(LEQ(I).NE.0) THEN
	II = II + 1
	RE(K) = RE(K) + RL(II)
	ENDIF
	ENDDO
	ENDIF
	

	IF(LSF.GT.0) THEN
	K = 0
	II = 0
	DO I = 1,NEF
	IF(LEQ(I).NE.0) II = II + 1
	JJ = II-1
	DO J = I,NEF
	IF(LEQ(J).NE.0) JJ = JJ + 1
	K = K + 1
	IF(LEQ(I).NE.0.AND.LEQ(J).NE.0) THEN
	S(K) = S(K) + AKG(II,JJ)
	ENDIF
	ENDDO
	ENDDO
	ENDIF
	
	
	


C	----------------------------------------
500	CONTINUE  !END PART
C	----------------------------------------


	DEALLOCATE(IE,EIJ,LEQN,LEQ)
	DEALLOCATE(COORD,EDIS)

1600	CONTINUE




	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================     
	SUBROUTINE TRAVDIS(NNM,ND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	UPDATE NODAL DISP AT THE TIME OF INSTALLATION OF TRAVELLER
	DIMENSION ND(1),NNX(2,NNM),NNO(2,NNM)
	DIMENSION IDOF(9)

	DIMENSION XYZO(3),XYZ(3),OFFS(3),BMAT(6,6)
	DIMENSION DISP(6),DISU(6),DISX(2,6)
C	----------------------------------------
C	CALLING IDOF
C	----------------------------------------
	DO I = 1,9
	CALL INTFILL('%DOF',IDOF(I),1,I,0)
	ENDDO


	CALL INTFILL('%NUB',NSN,1,1,0)
	CALL INTFILL('%NUB',NSF,1,2,0)	
	CALL INTFILL('%NUB',NSC,1,4,0)

	NX = 0  
	NN = 0
	NNX(1:2,1:NNM) = 0
	NNO(1:2,1:NNM) = 0
	DO 2000 INM = 1,NNM
	ISN = ND(INM)
	CALL INTFILL('$NXF',IX,1,ISN,0) !IX = EXISTING FLAG
	CALL INTFILL('$NXF',1 ,1,ISN,2) !UPDATE EXISTING NODE FLAG HERE
	IF(IX.GT.0) THEN
	NX = NX + 1
	NNX(1:2,NX) = [INM,ISN]
	ELSE
	NN = NN + 1
	NNO(1:2,NN) = [INM,ISN]
	ENDIF

2000	CONTINUE

C	RETURN HERE AFTER UPDATE NODAL EXISTING FLAG  
C	NO DISPLACEMENT UPDATE FOR TRAVELER
	RETURN
C	----------------------
	IF(NX.EQ.1) THEN  !FOR ONE EXISTING NODE
C	----------------------

	INX = NNX(2,1)
	XYZO(1:3) = 0.0D0 
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZO(ISC),ISC,INX,0) 
	ENDDO

	DO 2500 IN = 1,NN
	ISN = NNO(2,IN)

	XYZ(1:3) = 0.0D0 
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZ(ISC) ,ISC,ISN,0) 
	ENDDO
	OFFS(1:3) = 0.0D0
	DO ISC = 1,3
	OFFS(ISC) = XYZ(ISC) - XYZO(ISC)
	ENDDO


	BMAT(1:6,1:6) = 0.0D0
	BMAT(1,1) =  1.0D0
	BMAT(2,2) =  1.0D0
	BMAT(3,3) =  1.0D0
	BMAT(2,4) = -1.0D0*OFFS(3)
	BMAT(3,4) =  1.0D0*OFFS(2)
	BMAT(1,5) =  1.0D0*OFFS(3)
	BMAT(3,5) = -1.0D0*OFFS(1)
	BMAT(1,6) = -1.0D0*OFFS(2)
	BMAT(2,6) =  1.0D0*OFFS(1)
	BMAT(4,4) =  1.0D0
	BMAT(5,5) =  1.0D0
	BMAT(6,6) =  1.0D0


	DISP(1:6) = 0.0D0
	DO ISF = 1,NSF
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	CALL RELFILL('$XDS',EDISO,ISF,INX,0)  !CALLING TOTAL DISP.
	DISP(IDOF(ISF)) = EDISO
	ENDIF
	ENDDO

	DISU = MATMUL(BMAT,DISP)
	DO ISF = 1,NSF
      CALL RELFILL('$XDS',0.0D0,ISF,ISN,1) !INITIALIZE
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	CALL RELFILL('$XDS',DISU(IDOF(ISF)),ISF,ISN,1) !UPDATE NODAL DISP
	ENDIF
	ENDDO

2500	CONTINUE

C	----------------------
	ENDIF
C	----------------------


C	----------------------
	IF(NX.EQ.2) THEN  !FOR TWO EXISTING NODE (USING LINEAR INTERPOLATION+EXTRAPOLATION)
C	----------------------
	DO 3000 IX = 1,NX-1

	INX1 = NNX(2,IX)
	INX2 = NNX(2,IX+1)

	XYZO(1:3) = 0.0D0 
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZ2,ISC,INX2,0) 
      CALL RELFILL('$XYZ',XYZ1,ISC,INX1,0) 
	XYZO(ISC) = XYZ2 - XYZ1
	ENDDO	
	ELN = SQRT(XYZO(1)*XYZO(1) + XYZO(2)*XYZO(2) + XYZO(3)*XYZO(3))


	DISX(1:2,1:6) = 0.0D0
	DO ISF = 1,NSF
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	CALL RELFILL('$XDS',EDISO,ISF,INX1,0)  !CALLING TOTAL DISP.
	DISX(1,IDOF(ISF)) = EDISO
	CALL RELFILL('$XDS',EDISO,ISF,INX2,0)  !CALLING TOTAL DISP.
	DISX(2,IDOF(ISF)) = EDISO
	ENDIF
	ENDDO


C	-------------------
	DO 2800 IN = 1,NN
	INO = NNO(1,IN)
	ISN = NNO(2,IN)
	
C	-------

	XYZ(1:3) = 0.0D0
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZ1,ISC,INX1,0) 
      CALL RELFILL('$XYZ',XYZ2,ISC,ISN ,0) 
	XYZ(ISC) = XYZ2 - XYZ1
	ENDDO

	SIGN = XYZ(1)*XYZO(1) + XYZ(2)*XYZO(2) + XYZ(3)*XYZO(3)
	IF(SIGN.NE.0.0D0) THEN
	SIGN = SIGN/ABS(SIGN)
	ELSE
	SIGN = 1.0D0
	ENDIF

	SEGL = SQRT(XYZ(1)*XYZ(1) + XYZ(2)*XYZ(2) + XYZ(3)*XYZ(3))
	FAC = SIGN*SEGL/ELN


	DO ISF = 1,NSF
      CALL RELFILL('$XDS',0.0D0,ISF,ISN,1) !INITIALIZE
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN

C	if(IDOF(ISF).EQ.1.OR.IDOF(ISF).eq.2.OR.IDOF(ISF).EQ.6) then
	EDIS = DISX(1,IDOF(ISF)) + 
	1       FAC*(DISX(2,IDOF(ISF)) - DISX(1,IDOF(ISF)))
      CALL RELFILL('$XDS',EDIS,ISF,ISN,1)  !UPDATE NODAL DISP
C	endif

	ENDIF
	ENDDO

C	-------

2800	CONTINUE
C	-------------------

3000	CONTINUE

C	----------------------
	ENDIF
C	----------------------

C	----------------------
	IF(NX.GT.2) THEN  !FOR MORE THAN TWO EXISTING NODE (USING LINEAR INTERPOLATION)
C	----------------------
	DO 3500 IX = 1,NX-1

	INX1 = NNX(2,IX)
	INX2 = NNX(2,IX+1)

	XYZO(1:3) = 0.0D0 
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZ2,ISC,INX2,0) 
      CALL RELFILL('$XYZ',XYZ1,ISC,INX1,0) 
	XYZO(ISC) = XYZ2 - XYZ1
	ENDDO	
	ELN = SQRT(XYZO(1)*XYZO(1) + XYZO(2)*XYZO(2) + XYZO(3)*XYZO(3))


	DISX(1:2,1:6) = 0.0D0
	DO ISF = 1,NSF
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	CALL RELFILL('$XDS',EDISO,ISF,INX1,0)  !CALLING TOTAL DISP.
	DISX(1,IDOF(ISF)) = EDISO
	CALL RELFILL('$XDS',EDISO,ISF,INX2,0)  !CALLING TOTAL DISP.
	DISX(2,IDOF(ISF)) = EDISO
	ENDIF
	ENDDO

C	-------------------
	DO 3200 IN = 1,NN
	INO = NNO(1,IN)
	ISN = NNO(2,IN)
	
C	-------
	IF(INO.GE.NNX(1,IX).AND.INO.LE.NNX(1,IX+1)) THEN  !(ONLY NODE IN THE INTERVAL WILL BE UPDATE)

	XYZ(1:3) = 0.0D0
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZ1,ISC,INX1,0) 
      CALL RELFILL('$XYZ',XYZ2,ISC,ISN ,0) 
	XYZ(ISC) = XYZ2 - XYZ1
	ENDDO

	SIGN = XYZ(1)*XYZO(1) + XYZ(2)*XYZO(2) + XYZ(3)*XYZO(3)
	IF(SIGN.NE.0.0D0) THEN
	SIGN = SIGN/ABS(SIGN)
	ELSE
	SIGN = 1.0D0
	ENDIF

	SEGL = SQRT(XYZ(1)*XYZ(1) + XYZ(2)*XYZ(2) + XYZ(3)*XYZ(3))
	FAC = SIGN*SEGL/ELN

	DO ISF = 1,NSF
      CALL RELFILL('$XDS',0.0D0,ISF,ISN,1) !INITIALIZE
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	EDIS = DISX(1,IDOF(ISF)) + 
	1       FAC*(DISX(2,IDOF(ISF)) - DISX(1,IDOF(ISF)))
      CALL RELFILL('$XDS',EDIS,ISF,ISN,1)  !UPDATE NODAL DISP
	ENDIF
	ENDDO

	ENDIF
C	-------

3200	CONTINUE
C	-------------------


3500	CONTINUE

C	----------------------
	ENDIF
C	----------------------



	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE UPDTDIS
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	UPDATE NODAL DISP USING EXISTING NODE
	DIMENSION IDOF(9)

	DIMENSION XYZO(3),XYZ(3),OFFS(3),BMAT(6,6)
	DIMENSION DISP(6),DISU(6),DISX(2,6),VV(3)

      DIMENSION EDISX(6),DISJ(6)
      
	ALLOCATABLE ND(:),NNX(:,:),NNO(:,:)

C	READ COMMAND
C	CALL FREBUF
C	CALL FREECH

      NNM = 0
      CALL FREINT('N',NNM,1)	  !NUMBER OF GROUP NODES

	IF(NNM.LE.0) RETURN

	ALLOCATE(ND(NNM),NNX(2,NNM),NNO(2,NNM))
      ND(1:NNM) = 0
      CALL FREINT('D',ND,NNM)	  !NODES NUMBER IN GROUP


C	----------------------------------------
C	CALLING IDOF
C	----------------------------------------
	DO I = 1,9
	CALL INTFILL('%DOF',IDOF(I),1,I,0)
	ENDDO


	CALL INTFILL('%NUB',NSN,1,1,0)
	CALL INTFILL('%NUB',NSF,1,2,0)	
	CALL INTFILL('%NUB',NSC,1,4,0)

	NX = 0  
	NN = 0
	NNX(1:2,1:NNM) = 0
	NNO(1:2,1:NNM) = 0
	DO 2000 INM = 1,NNM
	ISN = ND(INM)
	CALL INTFILL('$NXF',IX,1,ISN,0) !IX = EXISTING FLAG
	IF(IX.GT.0) THEN
	NX = NX + 1
	NNX(1:2,NX) = [INM,ISN]
	ELSE
	NN = NN + 1
	NNO(1:2,NN) = [INM,ISN]
	ENDIF

2000	CONTINUE


C	----------------------
	IF(NX.EQ.1) THEN  !FOR ONE EXISTING NODE
C	----------------------

	INX = NNX(2,1)
	XYZO(1:3) = 0.0D0 
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZO(ISC),ISC,INX,0) 
	ENDDO

	DO 2500 IN = 1,NN
	ISN = NNO(2,IN)

	XYZ(1:3) = 0.0D0 
	DO ISC = 1,NSC
      CALL RELFILL('$XYZ',XYZ(ISC) ,ISC,ISN,0) 
	ENDDO
	OFFS(1:3) = 0.0D0
	DO ISC = 1,3
	OFFS(ISC) = XYZ(ISC) - XYZO(ISC)
	ENDDO


	BMAT(1:6,1:6) = 0.0D0
	BMAT(1,1) =  1.0D0
	BMAT(2,2) =  1.0D0
	BMAT(3,3) =  1.0D0
	BMAT(2,4) = -1.0D0*OFFS(3)
	BMAT(3,4) =  1.0D0*OFFS(2)
	BMAT(1,5) =  1.0D0*OFFS(3)
	BMAT(3,5) = -1.0D0*OFFS(1)
	BMAT(1,6) = -1.0D0*OFFS(2)
	BMAT(2,6) =  1.0D0*OFFS(1)
	BMAT(4,4) =  1.0D0
	BMAT(5,5) =  1.0D0
	BMAT(6,6) =  1.0D0


	DISP(1:6) = 0.0D0
	DO ISF = 1,NSF
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	CALL RELFILL('$XDS',EDISO,ISF,INX,0)  !CALLING TOTAL DISP.
	DISP(IDOF(ISF)) = EDISO
	ENDIF
	ENDDO

	DISU = MATMUL(BMAT,DISP)
	DO ISF = 1,NSF
      CALL RELFILL('$XDS',0.0D0,ISF,ISN,1) !INITIALIZE
	IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	CALL RELFILL('$XDS',DISU(IDOF(ISF)),ISF,ISN,1) !UPDATE NODAL DISP
	ENDIF
	ENDDO

2500	CONTINUE

C	----------------------
	ENDIF
C	----------------------



	IF(NX.GE.2) THEN  !FOR TWO EXISTING NODE (USING LAGRANGE INTERPOLATION+EXTRAPOLATION)
C	-------------------------
	DO 6000 IN = 1,NN
C	-------------------------
	INO = NNO(1,IN)
	ISN = NNO(2,IN)

	XYZ(1:3) = 0.0D0
	DO ISC = 1,NSC
        CALL RELFILL('$XYZ',XYZ(ISC),ISC,ISN,0) 
	ENDDO
	
	IF(NX.EQ.2) THEN
	  INX = NNX(2,1)
	  JNX = NNX(2,2)
	  VV(1:3) = 0.0D0
	  DO ISC = 1,NSC
	      CALL RELFILL('$XYZ',XYZI,ISC,INX,0) 
	      CALL RELFILL('$XYZ',XYZJ,ISC,JNX,0)
	      VV(ISC) = XYZJ - XYZI
	  ENDDO
	  VELN = SQRT(VV(1)*VV(1)+VV(2)*VV(2)+VV(3)*VV(3))
	  IF(VELN.NE.0.0D0) CALL SCALEN(VV,VV,VELN,3)
	ENDIF
	
C	-------------------
      EDISX(1:6) = 0.0D0
	DO JX = 1,NX
	    JNX = NNX(2,JX)
	    DISJ(1:6) = 0.0D0
	    DO ISF = 1,NSF
	        IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	            CALL RELFILL('$XDS',EDISO,ISF,JNX,0)  !CALLING TOTAL DISP.
	            DISJ(IDOF(ISF)) = EDISO
	        ENDIF
	    ENDDO
    	
	    HHH = 0.0D0 !LAGRANGE POLYNOMIAL
	    DO IX = 1,NX
	      IF(IX.NE.JX) THEN
	        INX = NNX(2,IX)
              DEVID = 0.0D0 ; DENOM = 0.0D0 ; PROJT = 0.0D0 
	        DO ISC = 1,NSC
                CALL RELFILL('$XYZ',XYZI,ISC,INX,0) 
                CALL RELFILL('$XYZ',XYZJ,ISC,JNX,0) 
                VD1 = (XYZ(ISC) - XYZI) ; VD2 = (XYZJ - XYZI)  
                IF(NX.EQ.2) THEN
                    IF(VELN.NE.0.0D0) THEN
                        DEVID = DEVID + VD1*VV(ISC)
                        DENOM = DENOM + VD2*VV(ISC)
                    ELSE
                        DEVID = DEVID + VD1*VD1
                        DENOM = DENOM + VD2*VD2
                    ENDIF    
                ELSE
                    DEVID = DEVID + VD1*VD1
                    DENOM = DENOM + VD2*VD2
                ENDIF
                PROJT = PROJT + VD1*VD2
	        ENDDO	
	        IF(NX.EQ.2) THEN
                    IF(VELN.NE.0.0D0) THEN
	                  DENOM = ABS(DENOM) 
	                  DEVID = ABS(DEVID) 
	              ELSE
	                  DENOM = SQRT(DENOM) 
	                  DEVID = SQRT(DEVID) 
	              ENDIF
	        ELSE
	          DENOM = SQRT(DENOM) 
	          DEVID = SQRT(DEVID) 
	        ENDIF  
	        VAL = DSIGN(DEVID/DENOM,PROJT)
	        IF(HHH.EQ.0.0D0) THEN
	          HHH = VAL
	        ELSE  
                HHH = HHH*VAL
              ENDIF  
            ENDIF  
          ENDDO

          EDISX(1:6) = EDISX(1:6) + DISJ(1:6)*HHH
      ENDDO
C	-------------------

	DO ISF = 1,NSF
          CALL RELFILL('$XDS',0.0D0,ISF,ISN,1) !INITIALIZE
	    IF(IDOF(ISF).GE.1.AND.IDOF(ISF).LE.6) THEN
	        EDIS = EDISX(IDOF(ISF))
              CALL RELFILL('$XDS',EDIS,ISF,ISN,1)  !UPDATE NODAL DISP
	    ENDIF
	ENDDO

C	-------------------------
6000	CONTINUE
C	-------------------------
      ENDIF



	DEALLOCATE(ND,NNX,NNO)


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================