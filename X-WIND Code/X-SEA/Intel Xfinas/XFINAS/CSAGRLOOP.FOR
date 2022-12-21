C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE GRLCON
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	CALLING NEG
	CALL INTFILL('%NUB',NEG,1,5,0)


      CALL INTFILL('%SOL',LOPT1,1,1,0)
      CALL INTFILL('%SOL',LOPT2,1,2,0)
      CALL INTFILL('%SOL',LOPT3,1,3,0)
      CALL INTFILL('%SOL',LOPT4,1,4,0)
      CALL INTFILL('%SOL',LOPT5,1,5,0)
      CALL INTFILL('%SOL',LOPT6,1,6,0)

C	----------------------------------------
C	LOOP OVER ELEMENT GROUP 
C	----------------------------------------	
	DO 5000 IEG = 1,NEG


C	STRESS AND RIGIDITY & MATERIAL STATE
	IF(LOPT1.NE.0.OR.LOPT2.NE.0) THEN
	CALL ELLUPD(IEG)
	ENDIF

C	ELEMENT FORCE AND STIFFNESS
	IF(LOPT3.NE.0.OR.LOPT4.NE.0.OR.LOPT5.NE.0) THEN
	CALL ELLEFF(IEG)
	ENDIF


C	----------------------------------------
5000	CONTINUE
C	----------------------------------------



	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ELLEFF(IEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	ELEMENT FORCE AND STIFFNESS
	DIMENSION NPM(10),NPI(10)
	ALLOCATABLE IF1(:),IF2(:),AF6(:)
	ALLOCATABLE COORD(:,:)
	ALLOCATABLE LMEQ(:),LREQ(:),RE(:),RL(:),SE(:),STLD(:),EXDI(:)

	ALLOCATABLE IAX(:),AAX(:)



      CALL LOCATN('-GSF',KGSF,NUM,NWK ,2) 
	CALL LOCATN('-LMA',KLMA,NUM,NEQ1,1)  


	CALL INTFILL('%NUB',NSN,1,1,0)
	CALL INTFILL('%NUB',NSF,1,2,0)
	CALL INTFILL('%NUB',NSC,1,4,0)


C	SOLUTION CONTROL ARRAY
      CALL INTFILL('%SOL',LOPT3,1,3,0)
      CALL INTFILL('%SOL',LOPT4,1,4,0)
      CALL INTFILL('%SOL',LOPT5,1,5,0)

	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) !
	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) !
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	CALL INTFILL('%GRP',NNM  ,8 ,IEG,0) !
	
	CALL INTFILL('%GRF',N1   ,1 ,IEG,0) !
	CALL INTFILL('%GRF',N2   ,2 ,IEG,0) !
	CALL INTFILL('%GRF',N6   ,6 ,IEG,0) !

	CALL INTFILL('%GRF',NDX1 ,11,IEG,0) !
	CALL INTFILL('%GRF',NDX2 ,12,IEG,0) !
	CALL INTFILL('%GRF',NDX6 ,16,IEG,0) !
	
	CALL INTFILL('%GRF',NFLI ,21,IEG,0) !
	CALL INTFILL('%GRF',NFLR ,22,IEG,0) !
	
	CALL INTFILL('%GRF',LASI ,23,IEG,0) !
	CALL INTFILL('%GRF',LASR ,24,IEG,0) !
	

	NEK = (NEF*NEF+NEF)/2

	ALLOCATE(LMEQ(NEF),LREQ(NEF),RE(NEF),RL(NEF),SE(NEK))
	ALLOCATE(COORD(3,NNM))
	ALLOCATE(IF1(N1),IF2(N2),AF6(N6))
	ALLOCATE(STLD(NEF),EXDI(NEF))
	ALLOCATE(IAX(LASI),AAX(LASR))

C	------------------------------------
C	LOOP OVER ELEMENT
C	------------------------------------
	DO 5000 IEL = 1,NELE

      READ(NFLI,REC=NELE+IEL,ERR=5) IXIST !READ EXISTING FLAG TO LAST RECORD OF ELEMENT INTEGER DATA (INCREASE THE SPEED PURPOSE)
      GOTO 6
5     IXIST = 0
6     IF(IXIST.LE.0) GOTO 5000
      
      READ(NFLI,REC=IEL) IAX
      READ(NFLR,REC=IEL) AAX

	SE(1:NEK) = 0.0D0
	RE(1:NEF) = 0.0D0
	RL(1:NEF) = 0.0D0
	EXDI(1:NEF) = 0.0D0
	STLD(1:NEF) = 0.0D0

C	--------------------------
C	EXTRACT NODAL COORDINATE
	NPM(1) = N1
	NPI(1) = N1
	CALL CALLINTG(IF1,IAX,NDX1,NPM,NPI,1,'READ')  
	DO INM = 1,NNM
	NOD = IF1(INM)
	COORD(1:3,INM) = 0.0D0
	DO ISC = 1,NSC
	CALL RELFILL('$XYZ',COORD(ISC,INM),ISC,NOD,0)  !GETTING HERE NODAL COORDINATE
	ENDDO
	ENDDO
C	--------------------------

	NPM(1) = N2
	NPI(1) = N2
	CALL CALLINTG(IF2,IAX,NDX2,NPM,NPI,1,'READ')  
	DO IEF = 1,NEF
	LREQ(IEF) = IF2(IEF)       !ORIGINAL EQUA
	LMEQ(IEF) = IF2(IEF+NEF)   !IF2(NEF+1:NEF+NEF) = CURRENT EQAUATION
	ENDDO


C	LSTF  0 = NO COMPUTE OF STIFFNESS   1 = COMPUTE STIFFNESS
C	LFOC  0 = NO COMPUTE OF ELM FORCE   1 = COMPUTE ELM FORCE
	LSTF = 0
	IF(LOPT4.GT.0.OR.LOPT5.GT.0) LSTF = 1
C	------------------------------------
	IF(LOPT3.EQ.1) THEN
C	--------------------------
	LFOC = 1
	CALL MDUOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)
	DO IEF = 1,NEF
	IEQ  = IF2(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-NML',RE(IEF),IEQ,1,2)    !STORE NON-MECHANICAL FORCE
	ENDDO
C	--------------------------
C	BACK UP NODAL FORCE FROM NON-MECHA
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')       
	NUM1 = 1   + NEF + NEF
	NUM2 = NEF + NEF + NEF
	AF6(NUM1:NUM2) = AF6(NUM1:NUM2) + RL(1:NEF)    
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')  
C	--------------------------
C	------------------------------------
	ENDIF
C	------------------------------------
	IF(LOPT3.EQ.2) THEN
C	--------------------------
	LFOC = 2
	CALL MDUOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)
	DO IEF = 1,NEF
	IEQ  = IF2(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-EIS',RE(IEF),IEQ,1,2)    !STORE MECHANICAL FORCE
	ENDDO

C	--------------------------
C	BACK UP NODAL FORCE FROM MECHA
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')        
	NUM1 = 1   + NEF 
	NUM2 = NEF + NEF 
	AF6(NUM1:NUM2) = RL(1:NEF)  
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')         
C	--------------------------
C	------------------------------------
	ENDIF
C	------------------------------------
	IF(LOPT3.EQ.3) THEN
C	--------------------------
	LFOC = 1
	CALL MDUOPT3(IEG,IEL,COORD,SE,RE,RL,0   ,LFOC,IAX,AAX)  !NO STIFFNESS COMPUTATION HERE
	DO IEF = 1,NEF
	IEQ  = IF2(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-NML',RE(IEF),IEQ,1,2)    !STORE NON-MECHANICAL FORCE 
	ENDDO
C	--------------------------
C	BACK UP NODAL FORCE FROM NON-MECHA
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')        
	NUM1 = 1   + NEF + NEF
	NUM2 = NEF + NEF + NEF
	AF6(NUM1:NUM2) = AF6(NUM1:NUM2) + RL(1:NEF)
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT') 
C	--------------------------
C	--------------------------
	LFOC = 2
	RE(1:NEF)   = 0.0D0
	CALL MDUOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)
	DO IEF = 1,NEF
	IEQ  = IF2(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-EIS',RE(IEF),IEQ,1,2)    !STORE MECHANICAL FORCE
	ENDDO
C	--------------------------
C	BACK UP NODAL FORCE FROM MECHA
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')       
	NUM1 = 1   + NEF 
	NUM2 = NEF + NEF 
	AF6(NUM1:NUM2) = RL(1:NEF)
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'WRIT')      
C	--------------------------
C	------------------------------------
	ENDIF
C	------------------------------------
	IF(LOPT3.EQ.4) THEN     !ASSEMBLY ONLY MECHANICAL FORCE NO BACK UP NO UPDATE
C	--------------------------
	LFOC = 2
	CALL MDUOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)
	DO IEF = 1,NEF
	IEQ  = IF2(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-EIS',RE(IEF),IEQ,1,2)    !STORE MECHANICAL FORCE
	ENDDO
C	------------------------------------
	ENDIF
C	------------------------------------


	IF(LOPT3.EQ.0.AND.LSTF.GT.0) THEN
	LFOC = 0
	CALL MDUOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)  !STIFFNESS COMPUTATION HERE IF LOPT3 EQ 0
	ENDIF


C	ASSEMBLE OF STIFFNESS
	IF(LOPT4.GT.0) THEN
	CALL STIFLAX(LREQ,SE,NEF,NNM)  !TRANSFORMATION DUE TO LOCAL SUPPORT
	CALL ASBSTF(SE,LMEQ,NEF) 
	ENDIF

C	------------------------------------
C	FORCE FROM INCREMENTAL PRESCRIBED DISP.
	IF(LOPT5.GT.0) THEN
	DO IEF = 1,NEF
	IEQ  = IF2(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('$XDI',EXDI(IEF),IEQ,1,0)  !GETTING HERE INCREMENTAL NODAL PRESCRIBED DISPLACEMENT
	ENDDO
	CALL MULTRI(SE,EXDI,STLD,NEF)
	DO IEF = 1,NEF
	IEQ  = IF2(IEF)   !IF2(1:NEF) = ORIGINAL EQAUATION
	CALL RELFILL('-STL',STLD(IEF),IEQ,1,2)  !STORE FORCE FROM INCREMENTAL NODAL PRESCRIBED DISPLACEMENT
	ENDDO
	ENDIF
C	------------------------------------

      WRITE(NFLI,REC=IEL) IAX
      WRITE(NFLR,REC=IEL) AAX

C	------------------------------------
5000	CONTINUE
C	------------------------------------

      
	DEALLOCATE(IAX,AAX)
	
	DEALLOCATE(LMEQ,LREQ,RE,RL,SE)
	DEALLOCATE(COORD)
	DEALLOCATE(IF1,IF2,AF6)
	DEALLOCATE(STLD,EXDI)


	RETURN
	END




C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ASBSTF(S,LM,NEF)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION S(1),LM(1)
C	STIFFNESS ASSEMBLY

	NWS = (NEF*NEF+NEF)/2
	CALL MESTIF(LM,S,NWS,NEF,'WRT')
	RETURN

C	NEQF = 1
C	NEQL = NEQ
	NPRE = 0

      NDI = 0
      DO 200 IEF=1,NEF
      II = LM(IEF)
C	IF (II.LT.NEQF .OR. II.GT.NEQL) GOTO 200

C	MI = MAXA(II) - NPRE
      CALL INTFILL('-LMA',MAXA,1,II,0) 
      MI = MAXA - NPRE

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
      CALL RELFILL('-GSF',S(KSS),1,KK,2) 

 220  KS = KS+NEF-JEF
 200  NDI = NDI+NEF-IEF


	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE MULTRI(A,B,R,NEQ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	MATRIX MULTIPLY FOR UPPER TRAIANGULAR SYM. MATRIX  
C	=============================================================
	DIMENSION A(1),B(1),R(1),C(NEQ,NEQ)
	DIMENSION MAXA(NEQ+1)


	K = 0
	DO I = 1,NEQ
	DO J = I,NEQ
	K = K + 1
	C(I,J) = A(K)
	ENDDO
	ENDDO
	DO I = 1,NEQ
	DO J = 1,NEQ
	C(J,I) = C(I,J)
	ENDDO
	ENDDO
	
	DO I = 1,NEQ
	R(I) = 0.0D0
	DO J = 1,NEQ
	R(I) = R(I) + C(I,J)*B(J)
	ENDDO
	ENDDO
	
	RETURN


C	NEK+NEK+NEK OPERATION
C	ROW
	K = 1
	DO I = 1,NEQ
	MAXA(I) = K
	DO J = I,NEQ
	K = K + 1
	ENDDO
	ENDDO
	MAXA(NEQ+1) = K

C	=====================================================
	DO IEQ = 1,NEQ
	R(IEQ) = 0.0D0
	LHIG = MAXA(IEQ+1)-MAXA(IEQ)
	MMM  = IEQ + LHIG - 1 
C	-----------------------------------------------------
	DO JEQ = IEQ,MMM                           !UPPER TRIANGULAR
	LEN  = JEQ - IEQ
	NUM  = MAXA(IEQ) + LEN
	R(IEQ) = R(IEQ) + A(NUM)*B(JEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDDO
C	-----------------------------------------------------
	DO JEQ = 1,IEQ-1                           !LOWER TRIANGULAR
	LHIG = MAXA(JEQ+1)-MAXA(JEQ)
	LHIG = LHIG - 1
	LEN  = JEQ - IEQ
	IF(LEN.LE.LHIG) THEN
	NUM  = MAXA(JEQ) + LEN
	R(IEQ) = R(IEQ) + A(NUM)*B(JEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDDO
C	-----------------------------------------------------
	ENDDO
C	=====================================================



C	NEK+NEK+NEK OPERATION
C	COL
	K = 1
	DO I = 1,NEQ
	MAXA(I) = K
	DO J = 1,I
	K = K + 1
	ENDDO
	ENDDO
	MAXA(NEQ+1) = K

C	=====================================================
	DO IEQ = 1,NEQ
	R(IEQ) = 0.0D0
C	-----------------------------------------------------
	DO JEQ = IEQ,NEQ	!UPPER TRIANGULAR
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1

	IF(LHEI.LE.IEQ) THEN
	R(IEQ) = R(IEQ) + A(NUM)*B(JEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF

	ENDDO
C	-----------------------------------------------------
	NUM = MAXA(IEQ+1) - MAXA(IEQ)
	LCL = MAXA(IEQ)-1
	DO J = 2,NUM		!LOWER TRIANGULAR
	JEQ = IEQ - J + 1

	IF(JEQ.GT.0) THEN
	NUM2 = LCL+J
	R(IEQ) = R(IEQ) + A(NUM2)*B(JEQ)           !R(IEQ) + A(LCL+J)*B(JEQ)
	ENDIF

	ENDDO
C	-----------------------------------------------------

C	-----------------------------------------------------
	ENDDO
C	=====================================================


	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================

      SUBROUTINE MDUOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION SE(1),RE(1),RL(1),COORD(3,1),IAX(1),AAX(1)
C	LSTF  0 = NO COMPUTE OF STIFFNESS   1 = COMPUTE STIFFNESS
C	LFOC  0 = NO COMPUTE OF ELM FORCE   1 = COMPUTE ELM FORCE

	CALL INTFILL('%GRP',ITYPE,1 ,IEG,0) 

C	------------------------------------
	SELECTCASE(ITYPE)
C	------------------------------------

	CASE(2)  !TRUSS ELEMENT
	CALL LINOPT3(IEG,IEL,COORD,SE,RE,LSTF,LFOC,IAX,AAX)
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	RL(1:NEF) = RE(1:NEF)
	
	CASE(5)  !FRAME ELEMENT
	CALL FMEOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)

	CASE(9)  !SHELL ELEMENT
	CALL SHEOPT3(IEG,IEL,COORD,SE,RE,LSTF,LFOC,IAX,AAX)
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	RL(1:NEF) = RE(1:NEF)
	
	CASE(10) !SOLID ELEMENT
	CALL SOEOPT3(IEG,IEL,COORD,SE,RE,LSTF,LFOC,IAX,AAX)
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !
	RL(1:NEF) = RE(1:NEF)
	
C	------------------------------------
	ENDSELECT
C	------------------------------------


	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FMEOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION SE(1),RE(1),RL(1),COORD(3,1),IAX(1),AAX(1)
C	LSTF  0 = NO COMPUTE OF STIFFNESS   1 = COMPUTE STIFFNESS
C	LFOC  0 = NO COMPUTE OF ELM FORCE   1 = COMPUTE ELM FORCE

	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) 

	SELECTCASE(ISTYP)

	CASE(5)
	CALL FRMOPT3(IEG,IEL,COORD,SE,RE,RL,LSTF,LFOC,IAX,AAX)

c	IF(IEL.EQ.10.AND.LFOC.EQ.2) THEN
c	WRITE(*,*) IEL,RE(1:14)
c	PAUSE
c	ENDIF


	ENDSELECT



	RETURN
	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE SOEOPT3(IEG,IEL,COORD,SE,RE,LSTF,LFOC,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION SE(1),RE(1),COORD(3,1),IAX(1),AAX(1)
C	LSTF  0 = NO COMPUTE OF STIFFNESS   1 = COMPUTE STIFFNESS
C	LFOC  0 = NO COMPUTE OF ELM FORCE   1 = COMPUTE ELM FORCE

	CALL INTFILL('%GRP',ISTYP,2 ,IEG,0) 

	SELECTCASE(ISTYP)

	CASE(1) ! SOLICO 
	CALL SOCOPT3(IEG,IEL,COORD,SE,RE,LSTF,LFOC,IAX,AAX)

	CASE(6) ! SOLIDEAS
	CALL SODOPT3(IEG,IEL,COORD,SE,RE,LSTF,LFOC,IAX,AAX)

	ENDSELECT



	RETURN
	END
	


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FRMOPT3(IEG,IEL,COORD,S,RE,RL,LSTF,LFOC,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	LTMOPT(3)
C	 0 = NO-ASSEMBLE OF RESISTING LOAD
C	 1 = ASSEMBLE OF NON-MECHANICAL INTERNAL LOAD
C	 2 = ASSEMBLE OF MECHANICAL INTERNAL LOAD
C	 3 = ASSEMBLE OF NON-MECHANICAL & MECHANICAL INTERNAL LOAD
C	 4 = ELEMENT FORCE + ASSEMBLE REACTION(FROM MECHANICAL) 

	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)


	DIMENSION DSIGM(5,5),QMAT(5,5),FMAT(5,5),TRANM(5,5),TRNFLX(14,5)
	DIMENSION IJ(4)

	DIMENSION PROPA(20)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION TT(14,14),TRANH(14,14),LREAS(14),TRANO(14,14),PROPO(6)
	DIMENSION BPG(10),BWG(10),BBX(14)
	DIMENSION NPM(10),NPI(10)

	DIMENSION BMATX(7,14)


	DIMENSION SFMAT(5,5),FSMAT(5,5),TRNSFX(14,5),EDIF(14),RLS(14)
	

	DIMENSION S(1),DR(7,7),SIGR(7),RE(14),RL(14),COORD(3,1)
	DIMENSION AKG(14,14),AKL(14,14),BMSG(7,14)
	
	DIMENSION IAX(1),AAX(1)

	ALLOCATABLE IF3(:),AF6(:),AF7(:)
	ALLOCATABLE COORO(:,:)
	ALLOCATABLE GPL(:),GPW(:)
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

	CALL INTFILL('%GRF',NDX1 ,11,IEG,0) !
	CALL INTFILL('%GRF',NDX2 ,12,IEG,0) !
	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX4 ,14,IEG,0) !
	CALL INTFILL('%GRF',NDX5 ,15,IEG,0) !
	CALL INTFILL('%GRF',NDX6 ,16,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
C	----------------------------------------
	ALLOCATE(IF3(N3),AF6(N6),AF7(N7))
	ALLOCATE(COORO(3,NNM))

C	---------------------------------------------------------
	ALLOCATE(GPL(NGR),GPW(NGR))

	DO IGR = 1,NGR
	IF(IGR.EQ.1  ) GPL(IGR) = -1.0D0
	IF(IGR.EQ.NGR) GPL(IGR) =  1.0D0 
	IF(IGR.NE.1.AND.IGR.NE.NGR) GPL(IGR) =  GAUSP(IGR-1,NGR-2)
	IF(IGR.EQ.1  ) GPW(IGR) =  0.0D0
	IF(IGR.EQ.NGR) GPW(IGR) =  0.0D0 
	IF(IGR.NE.1.AND.IGR.NE.NGR) GPW(IGR) =  GAUSW(IGR-1,NGR-2)
	ENDDO

C	---------------------------------------------------------	
	NPM(1) = N6
	NPI(1) = N6
	CALL CALLREAL(AF6,AAX,NDX6,NPM,NPI,1,'READ')
C	GET OFFSET VALUE
	DO I = 1,6
	PROPO(I) = AF6(I+NEF*4)
	ENDDO

C	---------------------------------------------------------
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')

	ISEC = IF3(2)     !SECTION SET
	
	IXIST = IF3(4)  !TOTAL EXISTING FLAG 
	IF(IXIST.EQ.0) GOTO 8000 !NO EXISTED ELEMENT 
			
C	GET RELEASE FLAG
	DO I = 1,14
	LREAS(I) = IF3(7+I)
	ENDDO

C	----------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	CALL TRNOOF(TRANO,PROPO)
	CALL TRNOFT1(COORD,COORO,EDIS,PROPO,TRANO)


C	CALCULATE ORIGINAL S BASE VECTOR
	DO I = 1,3
	VR(I) = COORO(I,2) - COORO(I,1)
	ENDDO
	CALL SCALEN(VR,VR,ELN,3)            !GET LENGTH OF ELEMENT HERE
	CALL FMVEVR(VR,VS,VT)
	CALL RELFILL('$RNG',RANG,1,ISEC,0)  !SECTION ROTATION
	CALL ROMBAC(VR,VS,VT,RANG)


C	OBTAIN END RELEASE TRANSFORMATION MATRIX
	CALL TRNHIG(TRANH,ELN,LREAS)

C	DETERMINE GAUSS POINT POSITION AND WEIGHT
	DO IGR = 1,NGR
	RI = GPL(IGR)  !GAUSP(IGR,NGR)
	RW = GPW(IGR)  !GAUSW(IGR,NGR)
	BPG(IGR) = 0.5*ELN*(1.0 + RI)
	BWG(IGR) = 0.5*ELN*RW
	ENDDO

C	INITILIZE ARRAY
	S(1:105)		= 0.0D0       !STIFF
	DR(1:7,1:7)		= 0.0D0       !RIGIDITY
	SIGR(1:7)		= 0.0D0       !RESULTANT STRESS
	RE(1:14)		= 0.0D0       !INTERNAL FORCE
	RL(1:14)		= 0.0D0       !INTERNAL FORCE LOCAL
	AKG(1:14,1:14)	= 0.0D0       !FULL STIFF


	TRANM = 0.0D0
	QMAT  = 0.0D0
	FMAT  = 0.0D0
	DISX  = 0.0D0
C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ----------------------------------------------------------
      DO 400 IGR = 1,NGR

	NPM(1:2) = [N7,NGR]
	NPI(1:2) = [N7,IGR]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  !CALLING RIGIDITY AND STRESS (MCHA & NON-MCHA)
	N7_R = 1           !POINTER OF RIGIDITY
	N7_S = 1 + 20       !POINTER OF MECHANICAL STRESS
	N7_N = 1 + 20 + 8   !POINTER OF NON-MECHANICAL STRESS

C	---------------------------
C	IF(LSTF.GT.0) THEN
	DO I = 1,20
	PROPA(I) = AF7(I+N7_R-1) !RIGIDITY
	ENDDO
C	ENDIF
	
	IF(LFOC.GT.0) THEN
	SELECTCASE(LFOC)

	CASE(1)
	DO I = 1,7
	SIGR(I) = AF7(I+N7_N-1)  !NON-MECHA
	ENDDO

	CASE(2)
	DO I = 1,7
	SIGR(I) = AF7(I+N7_S-1)  !MECHA
	ENDDO

	ENDSELECT

	ENDIF
C	---------------------------

C     GAUSS LOCATION ALONG ELEMENT AXIS
	BXD = BPG(IGR)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
C	CALL BBXFRM(BBX,ELN,BXD)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
	CALL BBXFRMW(BMATX,ELN,BXD)

C	!!!
C	IF(LSTF.GT.0) THEN
	EA   = PROPA(1)
	EQs  = PROPA(2)
	EQt  = PROPA(3)
	EIs  = PROPA(4)
	EIt  = PROPA(5)
	EIst = PROPA(6)
	GJr  = PROPA(7)
	RHOA = PROPA(8)
	EIw  = PROPA(9)
	EIsw = PROPA(10)
	EItw = PROPA(11)
	EQw  = PROPA(12)

	ASS  = PROPA(13)
	ATT  = PROPA(14)

	DR(1,1) = EA
	DR(4,4) = GJr

	DR(1,5) = EQs
	DR(5,5) = EIs

	DR(1,6) =-EQt
	DR(5,6) =-EIst
	DR(6,6) = EIt

	DR(1,7) = EQw
	DR(5,7) = EItw
	DR(6,7) =-EIsw
	DR(7,7) = EIw


	DO I = 1,7
	DO J = I,7
	DR(J,I) = DR(I,J)
	ENDDO
	ENDDO


C	---------------------------
C	D11 = EA 
C	D12 = EQt
C	D13 = EQs
C	D21 = EQt 
C	D22 = EIt
C	D23 = EIst
C	D31 = EQs
C	D32 = EIst
C	D33 = EIs

	D11 = DR(1,1) 
	D12 = DR(1,6)  
	D13 = DR(1,5)
	D21 = D12 
	D22 = DR(6,6)
	D23 =-DR(5,6)
	D31 = D13
	D32 = D23
	D33 = DR(5,5)

	DSIGM = 0.0D0

	DETS = D11*D22*D33+D21*D13*D32+D31*D12*D23
	1	  -D11*D23*D32-D21*D12*D33-D31*D13*D22
	DSIGM(1,1) = (D22*D33-D23*D32)/DETS
	DSIGM(1,2) = (D13*D32-D12*D33)/DETS
	DSIGM(1,3) = (D12*D23-D13*D22)/DETS

	DSIGM(2,1) = (D23*D31-D21*D33)/DETS
	DSIGM(2,2) = (D11*D33-D13*D31)/DETS
	DSIGM(2,3) = (D13*D21-D11*D23)/DETS

	DSIGM(3,1) = (D21*D32-D22*D31)/DETS
	DSIGM(3,2) = (D12*D31-D11*D32)/DETS
	DSIGM(3,3) = (D11*D22-D12*D21)/DETS

	IF(ASS.GT.0.0D0) DSIGM(4,4) = 1.0/ASS
	IF(ATT.GT.0.0D0) DSIGM(5,5) = 1.0/ATT

	DSIGM = DSIGM*BWG(IGR)
 
	TRANM(1,1) = 1.0D0
	TRANM(2,2) = ELN-BXD
	TRANM(3,3) =-1.0D0*(ELN-BXD)
	TRANM(3,4) = 1.0D0
	TRANM(2,5) = 1.0D0
	TRANM(4,2) = 1.0D0
	TRANM(5,3) = 1.0D0

	FMAT = FMAT + MATMUL(TRANSPOSE(TRANM),MATMUL(DSIGM,TRANM))

C	---------------------------



	AKL = MATMUL(TRANSPOSE(BMATX),MATMUL(DR,BMATX))

	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = AKG(I,J) + BWG(IGR)*AKL(I,J)
	ENDDO
	ENDDO
	
C	ENDIF
C	!!!

	IF(LFOC.GT.0) THEN
C	INTERNAL FORCE
	RL = RL + BWG(IGR)*MATMUL(TRANSPOSE(BMATX),SIGR)
	ENDIF


400	CONTINUE


c	IF(FMAT(1,1).NE.0.0) THEN
c	DO I =1,5
c	WRITE(*,*) FMAT(1:5,I)
c	ENDDO
c	PAUSE
c	ENDIF
C	------------------------------------------------------------
C	!!!
C	IF(LSTF.GT.0) THEN

C	CALCULATE DISP. CORRESPONDING TO INTERNAL FORCE RL 
C	FIRST GET THE STIFFNESS MATRIX AT THE END OF FRAME
	II = 0
	DO I = 8,13
	IF(I.NE.11) II = II + 1
	JJ = 0
	DO J = 8,13
	IF(J.NE.11) JJ = JJ + 1
	IF(I.NE.11.AND.J.NE.11) THEN
	SFMAT(II,JJ) = AKG(I,J)
	ENDIF
	ENDDO
	ENDDO
C	INVERSE TO GET FLEXIBILITY MATRIX
	FSMAT = 0.0
	CALL INVMATF(SFMAT,FSMAT,5,IB)
C	TRANSFORMATION MATRIX (TO GET ONLY FORCE AT THE END OF FRAME ELEM.)
	TRNSFX = 0.0D0
	TRNSFX(8 ,1) = 1.0D0
	TRNSFX(9 ,2) = 1.0D0
	TRNSFX(10,3) = 1.0D0
	TRNSFX(12,4) = 1.0D0
	TRNSFX(13,5) = 1.0D0
C	GET THE DISP. CORRESPONDING TO INTERNAL FORCE RL 
	EDIF = MATMUL(TRNSFX,MATMUL(MATMUL(FSMAT,TRANSPOSE(TRNSFX)),RL))


C	INCLUDE SHEAR EFFECT INTO STIFFNESS

	QMAT = 0.0
C	CALL INVMATRIX(FMAT,QMAT,5)
	CALL INVMATF(FMAT,QMAT,5,IB)
	IF(IB.EQ.1) GOTO 305

	
	TRNFLX = 0.0D0
	TRNFLX(1,1) = -1.0D0
	TRNFLX(2,2) = -1.0D0
	TRNFLX(3,3) = -1.0D0
	TRNFLX(5,4) = -1.0D0
	TRNFLX(6,5) = -1.0D0

	TRNFLX(5,3) =  ELN
	TRNFLX(6,2) = -ELN

	TRNFLX(8 ,1) =  1.0D0
	TRNFLX(9 ,2) =  1.0D0
	TRNFLX(10,3) =  1.0D0
	TRNFLX(12,4) =  1.0D0
	TRNFLX(13,5) =  1.0D0

C	STIFFNESS MATRIX (ALREADY INCLUDE SHEAR)
	AKL = MATMUL(TRNFLX,MATMUL(QMAT,TRANSPOSE(TRNFLX)))

	IJ(1:4) = [4,7,11,14]

C	CORRECTION OF INTERNAL FORCE DUE TO SHEAR DEFORMATION
	IF(LFOC.GT.0) THEN
	RLS = MATMUL(AKL,EDIF)
	DO II = 1,4
	I = IJ(II)
	RLS(I) = RL(I)   !TORSION AND WARPING USE THE ORIGINAL ONE
	ENDDO
	RL = RLS
	ENDIF


	DO II = 1,4
	I = IJ(II)
	DO J = 1,14
	AKL(I,J) = AKG(I,J)   !TORSION AND WARPING USE THE ORIGINAL ONE
	AKL(J,I) = AKG(J,I)   !TORSION AND WARPING USE THE ORIGINAL ONE
	ENDDO
	ENDDO
	AKG = AKL


305	CONTINUE

C	ENDIF
C	------------------------------------------------------------


C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
	IF(LFOC.GT.0) CALL TRNMUL(TRANH,RL,2)
	IF(LSTF.GT.0) CALL TRNMUM(TRANH,AKG)
C	------------------------------------------------------------

C	TRANSFORMATION TO GLOBAL SYSTEM
	CALL TT1A (VR,VS,VT,TT)
	IF(LFOC.GT.0) RE   = MATMUL(TT,RL)
	IF(LSTF.GT.0) AKG  = MATMUL(MATMUL(TT,AKG),TRANSPOSE(TT))

C	----------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	IF(LFOC.GT.0) CALL TRNMUL(TRANO,RE,2)
	IF(LSTF.GT.0) CALL TRNMUM(TRANO,AKG)
C	------------------------------------------------------------

	IF(LSTF.GT.0) THEN
	K = 0
	DO I = 1,14
	DO J = I,14
	K = K + 1
	S(K) = S(K) + AKG(I,J)
	ENDDO
	ENDDO
	ENDIF

      
8000	DEALLOCATE(IF3,AF6,AF7)
	DEALLOCATE(COORO)
	DEALLOCATE(GPL,GPW)


	RETURN

	END



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SODOPT3(IEG,IEL,COORD,S,RE,LSTF,LFOC,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      	

C	LTMOPT(3)
C	 0 = NO-ASSEMBLE OF RESISTING LOAD
C	 1 = ASSEMBLE OF NON-MECHANICAL INTERNAL LOAD
C	 2 = ASSEMBLE OF MECHANICAL INTERNAL LOAD
C	 3 = ASSEMBLE OF NON-MECHANICAL & MECHANICAL INTERNAL LOAD
C	 4 = ELEMENT FORCE + ASSEMBLE REACTION(FROM MECHANICAL) 

	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
C     SOLID ELEMENT -------------------------------------------------		 
C     COMMON /MMENH/ MM,MM1,MM2,NDIMC  
      COMMON /CSAMM/ MM,MM1,MM2,NDIMC     
C	BASIC MATRIX FOR COMPUTATION
      DIMENSION DP(6,6),H(21),HD(3,21),XJ(3,3),XJI(3,3),B(63),DISD(9)
	DIMENSION BM(6,24),DE(6,6),DEP(6,6),STRESS(6)
	DIMENSION COORD(3,1),EDIS(1),IAX(1),AAX(1)
	DIMENSION EPS(6),NPM(10),NPI(10)   
	DIMENSION SK(24,24),SKe(24,24)	
C     EAS METHOD -----------------------
	DIMENSION GE(6,MM)
	DIMENSION T0(6,6),TT(6,6),TT0(6,6)   
	DIMENSION XJO(3,3),TM(6,MM)
	DIMENSION EAS(6),ALPHA(MM),RE1(24),RE2(24),HINFC(MM)
      DIMENSION SED(MM,MM),SEL(MM,24),SEDI(MM,MM),RH(MM)
C     ---------------------------------
	DIMENSION SE(1),DR(6,6),SIGR(6),RE(24),S(300)
C     ---------------------------------

	ALLOCATABLE IF3(:),AF7(:)
	ALLOCATABLE NODEX(:)
C     ---------------------------------	
	     
      

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
	CALL INTFILL('%GRF',N3   ,3 ,IEG,0) !
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
C	---------------------------------------------------------

	ALLOCATE(IF3(N3),AF7(N7))
	NDEX = NNM - 8 + 1	
	ALLOCATE(NODEX(NDEX))
	
C	CALLING SECTION SET NUMBER
C	IF3-R(SECTION)-------------------------------------------
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')
	ISEC  = IF3(2)     !SECTION SET
	IXIST = IF3(4)     !TOTAL EXISTING FLAG 	
	IF(IXIST.EQ.0) GOTO 8000 !NO EXISTED ELEMENT 
C	---------------------------------------------------------
C     INITIALIZE ADTA FOR he
      DO I = 1,MM
	HINFC(I) = 0.0
	ENDDO

C	---------------------------------------------------------
C     EAS METHOD
C	---------------------------------------------------------	
C     COMPUTE JACOBIAN MATRIX AT THE NATURAL CENTER (R=0,S=0,T=0)
      CALL SHAP3D_S (0.0D0,0.0D0,0.0D0,H,HD,NODEX,NNM)
	CALL JACO3D_S (COORD,HD,XJO,XJI,DETO,IEL,NNM)	
	
C     COMPUTE T MATRIX (FOR EAS METHOD)	
      CALL MATRIXT_S (XJO,T0)      
      CALL INVMATRIX (T0,TT,6)      
      TT0 = TRANSPOSE(TT)
      
C	---------------------------------------------------------
      
C	INITILIZE ARRAY 
      SK   = 0.0D0
	SEL  = 0.0D0
	SED  = 0.0D0
	SEDI  = 0.0D0
	RH   = 0.0D0
	RE1  = 0.0D0
	RE2  = 0.0D0	
	
		
	S(1:300)		= 0.0D0       !STIFF
	DR(1:6,1:6)	= 0.0D0       !RIGIDITY
	SIGR(1:6)		= 0.0D0       !RESULTANT STRESS
	RE(1:24)		= 0.0D0       !INTERNAL FORCE
	

C     ==========================================================
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ==========================================================
C     SETTING INDEX FOR GAUSS POINT
	IPT = 0	
      DO 400 IGR = 1,NGR
      DO 400 IGS = 1,NGS
      DO 400 IGT = 1,NGT
	RI = GAUSP(IGR,NGR)
	SI = GAUSP(IGS,NGS)
	TI = GAUSP(IGT,NGT)
	WT = GAUSW(IGR,NGR)*GAUSW(IGS,NGS)*GAUSW(IGT,NGT)
	IPT = IPT + 1
	
C     POINTER FOR RIGIDITY AND MECHANICAL STRESS
C	AF7 -----------------------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  
	N7_R = 1            !POINTER OF RIGIDITY
	N7_S = 1 + 12       !POINTER OF MECHANICAL STRESS (6 STRESSES)
	N7_N = 1 + 12 + 6   !POINTER OF NON-MECHANICAL STRESS (6 STRESSES)
C	---------------------------------------------------------
		
C	---------------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C	--------------------------------------------------------- 
      CALL SHAP3D_S (RI,SI,TI,H,HD,NODEX,NNM)
      CALL JACO3D_S (COORD,HD,XJ,XJI,DET,IEL,NNM)          
      DVOL = WT*DET
C     STRAIN-DISPLACEMENT MATRIX      
      CALL BMATSLD_S (HD,XJI,B,BM,NNM)
C     FOR ENHANCED STRAIN METHOD AND COMPUTE MATRIX [GE]
	!CALL SOLIDGE_S (RI,SI,TI,DETO,DET,TT0,MM,GE)	 
	CALL MTEAS7(RI,SI,TI,TT0,DETO,DET,GE,MM)        
C     MATERIAL PROPERTIES (LINEAR)       
	CALL DMATSLD_S (DE)	      

C	---------------------------------------------------------
C     CALL RIGIDITY
C	---------------------------------------------------------
C     LATF 0 = NO COMPUTE STIFFNESS,  LATF 1 = COMPUTE STIFFNESS      
	IF(LSTF.GT.0.OR.LFOC.GT.0) THEN
	K = 0
	DO I = 1,3
	DO J = 1,3
	K = K + 1
	DR(I,J) = AF7(K+N7_R-1) !RIGIDITY AND INTEGRATION TERM
	ENDDO
	ENDDO
	DO I = 4,6	
	K = K + 1
	DR(I,I) = AF7(K+N7_R-1) !RIGIDITY AND INTEGRATION TERM	
	ENDDO
	ENDIF

C	---------------------------------------------------------
C     CALL MECHANICAL AND NON-MECHANICALSTRESS
C	---------------------------------------------------------
	IF(LFOC.GT.0) THEN
	
	SELECTCASE(LFOC)

	CASE(1) !NON-MECHANICAL STRESS
	DO I = 1,6
	SIGR(I) = AF7(I+N7_N-1)
	ENDDO

	CASE(2) !MECHANICAL STRESS
	DO I = 1,6
	SIGR(I) = AF7(I+N7_S-1) 
	ENDDO
	
	ENDSELECT
	ENDIF
	
	
	
C	---------------------------------------------------------	
C     CALCULATE THE ELEMENT STIFFNESS     
C	---------------------------------------------------------
      IF(LSTF.GT.0.OR.LFOC.GT.0) THEN		
C     -------------------------------------
C	COMPUTE ENHANCED STIFFNESS MATRIX SED
C	COMPUTE ENHANCED COUPLING  MATRIX SEL
C     -------------------------------------
C     THE STIFFNESS - Kaa 
	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(DR,GE))*DVOL
C	THE STIFFNESS - Kau
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(DR,BM))*DVOL

C     ----------------------------------------      
C     THE STIFFNESS - Kuu
      SK = (MATMUL(TRANSPOSE(BM),MATMUL(DR,BM)))*DVOL      
      K = 0
      DO ISK =1,24
      DO JSK =ISK,24
      K = K+1 
      S(K)=S(K)+SK(ISK,JSK)      
      ENDDO
      ENDDO	    
      
      
      ENDIF
C	---------------------------------------------------------

      
	
C	---------------------------------------------------------
C	CALCULATE INTERNAL FORCE 
C	---------------------------------------------------------
      IF(LFOC.GT.0) THEN
      
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT - Fe
	RE1 = RE1 + MATMUL(TRANSPOSE(BM),SIGR)*DVOL

C     COMPUTE EQUILIBRIUM FORCE FOR EAS METHOD - he
	RH = RH + MATMUL(TRANSPOSE(GE),SIGR)*DVOL
		
	ENDIF
		
C	---------------------------------------------------------
400	CONTINUE  ! GAUSS LOOP 
  
      
    
C     THE STATIC CONDENSATION 
C	---------------------------------------------------------
C     STIFFNESS FROM EAS
C	---------------------------------------------------------
      !IF(LSTF.GT.0) THEN
C     Kaa -> Kaa**-1 
      CALL INVMATRIX (SED,SEDI,MM)
C     Kuu = Kuu - Kua[Kaa**-1]Kau      
      CALL MATRIXES_S (SEL,SEDI,S,MM)
      !ENDIF 

C	---------------------------------------------------------      
C	INTERNAL FORCE FROM EAS
C	---------------------------------------------------------
      IF(LFOC.GT.0) THEN
C     Kua[Kaa**-1]he
      !CALL INVMATRIX (SED,SEDI,MM)
      !CALL MATRIXES_S (SEL,SEDI,S,MM)
	RE2 = MATMUL(MATMUL(TRANSPOSE(SEL),SEDI),RH)
C     Fe = Fe - Kua[Kaa**-1]he	
	DO I=1,24
	RE1(I) = RE1(I) - RE2(I)	
	RE(I)  = RE(I) + RE1(I) 
	END DO
	ENDIF
C     ----------------------------------------------------	
	
C     UPDATE - he
	DO IHH = 1,MM
	HINFC(IHH) = RH(IHH)
	ENDDO	
C	---------------------------------------------------------

8000	DEALLOCATE(IF3,AF7)
	DEALLOCATE(NODEX)


	RETURN
	END
	
	
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SOCOPT3(IEG,IEL,COORDI,S,RE,LSTF,LFOC,IAX,AAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      

C	LTMOPT(3)
C	 0 = NO-ASSEMBLE OF RESISTING LOAD
C	 1 = ASSEMBLE OF NON-MECHANICAL INTERNAL LOAD
C	 2 = ASSEMBLE OF MECHANICAL INTERNAL LOAD
C	 3 = ASSEMBLE OF NON-MECHANICAL & MECHANICAL INTERNAL LOAD
C	 4 = ELEMENT FORCE + ASSEMBLE REACTION(FROM MECHANICAL) 

	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
C     SOLID ELEMENT -------------------------------------------------		 
C      COMMON /MMENH/ MM,MM1,MM2,NDIMC 
      COMMON /CSAMM/ MM,MM1,MM2,NDIMC   
C	BASIC MATRIX FOR COMPUTATION
      DIMENSION DP(6,6),H(21),HD(3,21),XJ(3,3),XJI(3,3),B(63),DISD(9)
	DIMENSION BM(6,24),DE(6,6),DEP(6,6),STRESS(6)
	DIMENSION COORD(3,8),EDIS(1),IAX(1),AAX(1)
	DIMENSION EPS(6),NPM(10),NPI(10)   
	DIMENSION SK(24,24),SKe(24,24)	
C     EAS METHOD -----------------------
	DIMENSION GE(6,MM)
	DIMENSION T0(6,6),TT(6,6),TT0(6,6)   
	DIMENSION XJO(3,3),TM(6,MM)
	DIMENSION EAS(6),ALPHA(MM),RE1(24),RE2(24),HINFC(MM)
      DIMENSION SED(MM,MM),SEL(MM,24),SEDI(MM,MM),RH(MM),D(6,6)
      DIMENSION SEL1(MM,24)	
C     ---------------------------------
	DIMENSION SE(1),DR(6,6),SIGR(6),RE(24),S(300)
C     ---------------------------------

      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION COORDI(3,8)	
	DIMENSION XY(3,8),BB(6,24)   
	DIMENSION EDISO(24),EDISI(24)
	DIMENSION EDISLO(24)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
C	DIMENSION FIN(NEF)
	DIMENSION REP(24),SEP1(24,24),SEP2(24,24),EDISP(24),MPNOD(8),MPDIS(24),MPGUS(8)
C     REMOVE RIGID BODY MOTION
      DIMENSION COORDO(3,8),REDIS(24)
	
C     ---------------------------------
	ALLOCATABLE IF3(:),AF7(:)
	ALLOCATABLE NODEX(:)
C     ---------------------------------		     
      

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
	CALL INTFILL('%GRF',N3   ,3 ,IEG,0) !
	CALL INTFILL('%GRF',N7   ,7 ,IEG,0) !

	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	CALL INTFILL('%GRF',NDX7 ,17,IEG,0) !
C	---------------------------------------------------------

	ALLOCATE(IF3(N3),AF7(N7))
	NDEX = NNM - 8 + 1	
	ALLOCATE(NODEX(NDEX))
	
C	CALLING SECTION SET NUMBER
C	IF3-R(SECTION)-------------------------------------------
	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')
	ISEC  = IF3(2)     !SECTION SET
	IXIST = IF3(4)     !TOTAL EXISTING FLAG 	
	IF(IXIST.EQ.0) GOTO 8000 !NO EXISTED ELEMENT 
C	---------------------------------------------------------
C     INITIALIZE ADTA FOR he
      DO I = 1,MM
	HINFC(I) = 0.0
	ENDDO
	
C	----------------------------------------------------------
C	DETECTED THE THIN DIRECTION
C	----------------------------------------------------------
	XYZ1 = SQRT( (COORDI(1,2)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,2)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,2)-COORDI(3,1))**2.0 )
	XYZ2 = SQRT( (COORDI(1,4)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,4)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,4)-COORDI(3,1))**2.0 )
	XYZ3 = SQRT( (COORDI(1,5)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,5)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,5)-COORDI(3,1))**2.0 )

	NMX = 0
	IF(XYZ1.LT.XYZ2.AND.XYZ1.LT.XYZ3) NMX = 1
	IF(XYZ2.LT.XYZ1.AND.XYZ2.LT.XYZ3) NMX = 2
	IF(XYZ3.LT.XYZ1.AND.XYZ3.LT.XYZ2) NMX = 3

	IF(NMX.EQ.0) THEN
	IF(XYZ1.EQ.XYZ2) NMX = 1
	IF(XYZ1.EQ.XYZ3) NMX = 1
	IF(XYZ2.EQ.XYZ3) NMX = 2
	ENDIF

	CALL SLCOFACE(NMX,MPNOD,MPDIS,MPGUS)
	
C     FLIPPING COORDINATE
      DO I = 1,NNM
      NOD = MPNOD(I)
      COORD(1:3,I) = COORDI(1:3,NOD)
      ENDDO
      
C     FLIPPING DISPLACEMENT
      EDISP(1:24) = EDIS(1:24)
      DO I = 1,24
      IEF = MPDIS(I)
      EDISO(IEF) = EDISP(I)
      ENDDO
C     ------------------------------------------------------------
C	FIND THE LOCAL VECTORS,LOCAL COORDINATES AND DISPLACEMENTS
C     -----------------------------------------------------------   
C     COMPUTE JACOBIAN MATRIX AT THE NATURAL CENTER (R=0,S=0,T=0)   
      CALL SHAP3D (0.0D0,0.0D0,0.0D0,H,HD,NODEX,NNM)
C     CALL SHAP3D_S (0.0D0,0.0D0,0.0D0,H,HD,NODEX,NNM)      
  	CALL SOLIRST(COORD,HD,VR,VS,VT)   ! LOCAL VECTORS
  	
C     REMOVE RIGID BODY MOTION
  	REDIS(1:24) = EDISO(1:24)
      IF(NLOPT.EQ.3) THEN
          K = 0
          DO I = 1,NNM
              DO J = 1,3
                  K = K + 1
                  COORDO(J,I) = COORD(J,I) - EDISO(K)
              ENDDO
          ENDDO
          CALL SOMDSP(COORD,COORDO,EDISO,REDIS,H,HD,VR,VS,VT,NNM)
      ENDIF
C     REDIS = DISPLACEMENT REMOVING RIGID BOBY      
      CALL SCOOR(COORD,REDIS,XY,EDISLO,VR,VS,VT) !REDIS INSTEAD OF EDIS HERE 
C     ----------------------------------------------------------     	

C	---------------------------------------------------------
C     EAS METHOD
C	---------------------------------------------------------	
C     COMPUTE JACOBIAN MATRIX AT THE NATURAL CENTER (R=0,S=0,T=0)
      CALL JACO3D (XY,HD,XJO,XJI,DETO,IEL,NNM)	
	!CALL JACO3D_S (COORD,HD,XJO,XJI,DETO,IEL,NNM)	
	
C     COMPUTE T MATRIX (FOR EAS METHOD)	
      CALL TNEAS6(XJO,TT0)   
C	---------------------------------------------------------
      
C	INITILIZE ARRAY 
      SK   = 0.0D0
	!SEL  = 0.0D0
	!SED  = 0.0D0
	!RH   = 0.0D0
	RE1  = 0.0D0
	RE2  = 0.0D0	
	
C     SET ZEROS STIFFNESS FOR EAS	
	SED  = 0.0D0                                
	SEDI = 0.0D0                                 
	SEL  = 0.0D0                                 
	RH   = 0.0D0                               
	SG = 0.0D0 	
	
		
	S(1:300)		= 0.0D0       !STIFF
	DR(1:6,1:6)	= 0.0D0       !RIGIDITY
	SIGR(1:6)		= 0.0D0       !RESULTANT STRESS
	RE(1:24)		= 0.0D0       !INTERNAL FORCE
	

C     ==========================================================
C     LOOP OVER GAUSS TO DET. STIFFNESS & FORCE VECTOR
C     ==========================================================
C     SETTING INDEX FOR GAUSS POINT
	IPT = 0	
      DO 400 IGR = 1,NGR
      DO 400 IGS = 1,NGS
      DO 400 IGT = 1,NGT
	RI = GAUSP(IGR,NGR)
	SI = GAUSP(IGS,NGS)
	TI = GAUSP(IGT,NGT)
	WT = GAUSW(IGR,NGR)*GAUSW(IGS,NGS)*GAUSW(IGT,NGT)
	IPT = IPT + 1
	
C     POINTER FOR RIGIDITY AND MECHANICAL STRESS
C	AF7 -----------------------------------------------------
	NPM(1:2) = [N7,MGAS]
	NPI(1:2) = [N7,IPT ]
	CALL CALLREAL(AF7,AAX,NDX7,NPM,NPI,2,'READ')  
	N7_R = 1            !POINTER OF RIGIDITY
	N7_S = 1 + 12       !POINTER OF MECHANICAL STRESS (6 STRESSES)
	N7_N = 1 + 12 + 6   !POINTER OF NON-MECHANICAL STRESS (6 STRESSES)
C	---------------------------------------------------------
		
C	---------------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C	--------------------------------------------------------- 
      !CALL SHAP3D_S (RI,SI,TI,H,HD,NODEX,NNM)
      CALL SHAP3D (RI,SI,TI,H,HD,NODEX,NNM)
      !CALL JACO3D_S (COORD,HD,XJ,XJI,DET,IEL,NNM) 
      CALL JACO3D (XY,HD,XJ,XJI,DET,IEL,NNM)         
      DVOL = WT*DET
C     STRAIN-DISPLACEMENT MATRIX      
      !CALL BMATSLD_S (HD,XJI,B,BM,NNM)
50    CALL SOBMAT (HD,XJI,B,NNM)      
C     FOR ENHANCED STRAIN METHOD AND COMPUTE MATRIX [GE]
	!CALL SOLIDGE_S (RI,SI,TI,DETO,DET,TT0,MM,GE)
	CALL MTEAS7(RI,SI,TI,TT0,DETO,DET,GE,MM)	         
C     MATERIAL PROPERTIES (LINEAR)       
	!CALL DMATSLD_S (DE)	
C     ASSUMED STRAIN-DISPLACEMENT MATRIX (BB)
      CALL SOBMATANS(XY,BM,B,RI,SI,TI,NNM,IEL,NODEX)   ! NEW ANS METHOD			      

C	---------------------------------------------------------
C     CALL RIGIDITY
C	---------------------------------------------------------
C     LATF 0 = NO COMPUTE STIFFNESS,  LATF 1 = COMPUTE STIFFNESS      
	IF(LSTF.GT.0.OR.LFOC.GT.0) THEN
	K = 0
	DO I = 1,3
	DO J = 1,3
	K = K + 1
	DR(I,J) = AF7(K+N7_R-1) !RIGIDITY AND INTEGRATION TERM
	ENDDO
	ENDDO
	DO I = 4,6	
	K = K + 1
	DR(I,I) = AF7(K+N7_R-1) !RIGIDITY AND INTEGRATION TERM	
	ENDDO
	ENDIF

C	---------------------------------------------------------
C     CALL MECHANICAL AND NON-MECHANICALSTRESS
C	---------------------------------------------------------
	IF(LFOC.GT.0) THEN
	
	SELECTCASE(LFOC)

	CASE(1) !NON-MECHANICAL STRESS
	DO I = 1,6
	SIGR(I) = AF7(I+N7_N-1)
	ENDDO

	CASE(2) !MECHANICAL STRESS
	DO I = 1,6
	SIGR(I) = AF7(I+N7_S-1) 
	ENDDO
	
	ENDSELECT
	ENDIF
	
	
	
C	---------------------------------------------------------	
C     CALCULATE THE ELEMENT STIFFNESS     
C	---------------------------------------------------------
      IF(LSTF.GT.0.OR.LFOC.GT.0) THEN		
C     -------------------------------------
C	COMPUTE ENHANCED STIFFNESS MATRIX SED
C	COMPUTE ENHANCED COUPLING  MATRIX SEL
C     -------------------------------------
C     THE STIFFNESS - Kaa 
	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(DR,GE))*DVOL
C	THE STIFFNESS - Kau
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(DR,BM))*DVOL

C     ----------------------------------------      
C     THE STIFFNESS - Kuu
      SK = (MATMUL(TRANSPOSE(BM),MATMUL(DR,BM)))*DVOL      
      K = 0
      DO ISK =1,24
      DO JSK =ISK,24
      K = K+1 
      S(K)=S(K)+SK(ISK,JSK)      
      ENDDO
      ENDDO	    
      
      
      ENDIF
C	---------------------------------------------------------

      
	
C	---------------------------------------------------------
C	CALCULATE INTERNAL FORCE 
C	---------------------------------------------------------
      IF(LFOC.GT.0) THEN
      
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT - Fe
	RE1 = RE1 + MATMUL(TRANSPOSE(BM),SIGR)*DVOL

C     COMPUTE EQUILIBRIUM FORCE FOR EAS METHOD - he
	RH = RH + MATMUL(TRANSPOSE(GE),SIGR)*DVOL
		
	ENDIF
		
C	---------------------------------------------------------
400	CONTINUE  ! GAUSS LOOP 
  
      
    
C     THE STATIC CONDENSATION 
C	---------------------------------------------------------
C     STIFFNESS FROM EAS
C	---------------------------------------------------------
      !IF(LSTF.GT.0) THEN
C     Kaa -> Kaa**-1 
      CALL INVMATRIX (SED,SEDI,MM)
C     Kuu = Kuu - Kua[Kaa**-1]Kau      
      CALL MATRIXES_S (SEL,SEDI,S,MM)
      !ENDIF 

C	---------------------------------------------------------      
C	INTERNAL FORCE FROM EAS
C	---------------------------------------------------------
      IF(LFOC.GT.0) THEN
C     Kua[Kaa**-1]he
      !CALL INVMATRIX (SED,SEDI,MM)
      !CALL MATRIXES_S (SEL,SEDI,S,MM)
	RE2 = MATMUL(MATMUL(TRANSPOSE(SEL),SEDI),RH)
C     Fe = Fe - Kua[Kaa**-1]he	
	DO I=1,24
	RE1(I) = RE1(I) - RE2(I)	
	RE(I)  = RE(I) + RE1(I) 
	END DO
	ENDIF
C     ----------------------------------------------------	
	
C     UPDATE - he
	DO IHH = 1,MM
	HINFC(IHH) = RH(IHH)
	ENDDO	
C	---------------------------------------------------------

      CALL STIFFTRANS(VR,VS,VT,S,RE)	
	CALL STIFFTRANE(VR,VS,VT,MM,SEL) !TRANSFORM EAS MATRIX TO GLOBAL
C     ----------------------------------------------------------------  

C     ---------------------------------------
C     FLIPPING DUE TO THICKNESS DETECTION
C     ---------------------------------------
      REP(1:24)  = RE(1:24)
      DO IEF = 1,24
          II = MPDIS(IEF)
          RE(IEF) = REP(II)
      ENDDO
      
      SEP1(1:24,1:24) =  0.0D0
      SEP2(1:24,1:24) =  0.0D0
      KEF = 0
      DO IEF = 1,24
          DO JEF = IEF,24
              KEF = KEF+1
              SEP1(IEF,JEF) = S(KEF)
              SEP1(JEF,IEF) = S(KEF)
              SEP2(IEF,JEF) = S(KEF)
              SEP2(JEF,IEF) = S(KEF)
          ENDDO
      ENDDO

      DO IEF = 1,24
          II = MPDIS(IEF)
          DO JEF = 1,24
              JJ = MPDIS(JEF)
              SEP1(IEF,JEF) = SEP2(II,JJ)
          ENDDO
      ENDDO
      
C     REARRANGE EAS MATRIX      
      SEL1 = SEL
      DO JEF = 1,24
          JJ = MPDIS(JEF)
          SEL(1:MM,JEF) = SEL1(1:MM,JJ)
      ENDDO
      
      KEF = 0
      DO IEF = 1,24
          DO JEF = IEF,24
              KEF = KEF+1
              S(KEF) = SEP1(IEF,JEF)
          ENDDO
      ENDDO        
  
          
C     ---------------------------------------     

8000	DEALLOCATE(IF3,AF7)
	DEALLOCATE(NODEX)


	RETURN
	END	


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================





