C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE STIFLAX(LNQ,S,NEF,NNM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION LNQ(1),S(1)
	DIMENSION NF(NEF),NL(NEF),LCL(NEF)
	DIMENSION MAXA(NEF+1),MMN(NNM),MLA(NNM)

	DIMENSION AA(NEF*NEF),A(NEF,NEF)
C	STIFFNESS TRANSFORMATION DUE TO LOCAL SUPPORT TT*S*T
	ALLOCATABLE MNF(:,:),B(:,:)


	CALL INTFILL('%NUB',NSF,1,2,0)


	K = 0
	DO I = 1,NEF
	DO J = I,NEF
	K = K + 1
	A(I,J) = S(K)
	A(J,I) = S(K)
	ENDDO
	ENDDO

	K = 0
	DO I = 1,NEF
	IJ = I
	DO J = 1,I
	K = K + 1
	AA(K) = A(I,IJ)
	IJ = IJ-1
	ENDDO
	ENDDO

C	ROW WISE
	K = 1
	DO I = 1,NEF
	MAXA(I) = K
	DO J = I,NEF
	K = K + 1
	ENDDO
	ENDDO
	MAXA(NEF+1) = K

	K = 1
	DO I = 1,NEF
	MAXA(I) = K
	DO J = 1,I
	K = K + 1
	ENDDO
	ENDDO
	MAXA(NEF+1) = K


	NNF = NEF/NNM
	ALLOCATE(MNF(NNF,NNM))
C	----------------------------
C	MAPPING MATRIX
C	----------------------------
	IEF  = 0
	MCL  = 0
	DO INM = 1,NNM
	MMN(INM) = 0
	MLA(INM) = 0
	MNF(1:NNF,INM) = 0
	DO INF = 1,NNF
	IEF = IEF + 1
	IEQ = LNQ(IEF)                    !ORIGINAL EQUATION
	IF(IEQ.LE.0) GOTO 10
	ISN = (IEQ-1)/NSF
	ISF = IEQ - NSF*ISN
	ISN = ISN+1
	CALL INTFILL('$NBC',IAX,1,ISN,0)  !CALLING LAX SET
	CALL INTFILL('$MBC',IBC,IEQ,1,0)  !CALLING BOUNDARY CONDITION FLAG
	IF(ISF.LE.6.AND.IBC.GT.0.AND.IAX.NE.0) THEN    !FOR LOCAL SUPPORT
	MCL  = 1
	MMN(INM)  = 1
	MLA(INM)  = IAX
	ENDIF
	IF(ISF.LE.6) THEN
	MNF(INF,INM) = ISF
	ENDIF
10	CONTINUE
	ENDDO
	ENDDO

C	----------------------------
	IF(MCL.LE.0) THEN
	DEALLOCATE(MNF)
	RETURN
	ENDIF
C	----------------------------


	NEFO = 0  
	DO INM = 1,NNM
	IFG = MMN(INM)
	IF(IFG.LE.0) GOTO 100
	DO INF = 1,NNF
	ISF = MNF(INF,INM)
	IF(ISF.GT.0) NEFO = NEFO + 1
	ENDDO
100	CONTINUE
	ENDDO


	ALLOCATE(B(NEFO,NEFO))
	B(1:NEFO,1:NEFO) = 0.0D0
	LCL(1:NEFO) = 0


C	----------------------------
C	TRANSFORMATION MATRIX
C	----------------------------
	NEFO = 0  
	DO INM = 1,NNM
	NEFJ = NEFO
	IFG = MMN(INM)
	IF(IFG.LE.0) GOTO 200
	IAX = MLA(INM)
	DO INF = 1,NNF
	ISF = MNF(INF,INM)

	IF(ISF.EQ.0) GOTO 150
	NEFO = NEFO + 1
	LCL(NEFO) = INF+NNF*(INM-1)
	NN   = NEFJ
	NN3  = 0
	DO JNF = 1,NNF
	JSF = MNF(JNF,INM)
	IF(JSF.EQ.0) GOTO 120

	IF(JSF.LE.3) NN3 = NN3+1
	
	IF((ISF.GE.1.AND.ISF.LE.3).AND.(JSF.GE.1.AND.JSF.LE.3)) THEN
	NN = NN + 1
	ISF3 = ISF
	JSF3 = JSF
	IF(ISF3.GT.3) ISF3 = ISF3-3
	IF(JSF3.GT.3) JSF3 = JSF3-3
	IDM = ISF3 + 3*(JSF3-1)
	CALL RELFILL('$LAX',VAX,IDM,IAX,0)
	B(NEFO,NN) = VAX
	ENDIF

	IF((ISF.GE.4.AND.ISF.LE.6).AND.(JSF.GE.4.AND.JSF.LE.6)) THEN
	NN = NN + 1
	ISF3 = ISF
	JSF3 = JSF
	IF(ISF3.GT.3) ISF3 = ISF3-3
	IF(JSF3.GT.3) JSF3 = JSF3-3
	IDM = ISF3 + 3*(JSF3-1)
	CALL RELFILL('$LAX',VAX,IDM,IAX,0)
	B(NEFO,NN+NN3) = VAX
	ENDIF
		
120	CONTINUE
	ENDDO

150	CONTINUE
	ENDDO
200	CONTINUE
	ENDDO
C	----------------------------




C	PERFORM TRANSFORMATION
	CALL TRANSYM1(MAXA,AA,B,NEF,NEFO,LCL)


	K = 0
	DO I = 1,NEF
	IJ = I
	DO J = 1,I
	K = K + 1
	A(I ,IJ) = AA(K)
	A(IJ,I ) = AA(K)
	IJ = IJ-1
	ENDDO
	ENDDO


	K = 0
	DO I = 1,NEF
	DO J = I,NEF
	K = K + 1
	S(K) = A(I,J)
	ENDDO
	ENDDO


	DEALLOCATE(MNF)
	DEALLOCATE(B)



	RETURN

	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE FORCLAX(ISN,C,NSF,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION C(1)
	DIMENSION NL(NSF),LCL(NSF)
C	FORCE TRANSFORMATION DUE TO LOCAL SUPPORT TT*S*T
	ALLOCATABLE MNF(:),B(:,:)


	CALL INTFILL('%NUB',NSF,1,2,0)


	ALLOCATE(MNF(NSF))
C	----------------------------
C	MAPPING MATRIX
C	----------------------------
	MCL  = 0
	MMN = 0
	MLA = 0
	MNF(1:NSF) = 0
	DO ISF = 1,NSF
	CALL INTFILL('$NBC',IAX,1,ISN,0)    !CALLING LAX SET
	CALL INTFILL('$MBC',IBC,ISF,ISN,0)  !CALLING BOUNDARY CONDITION FLAG
	IF(ISF.LE.6.AND.IBC.GT.0.AND.IAX.NE.0) THEN    !FOR LOCAL SUPPORT
	MCL  = 1
	MMN  = 1
	MLA  = IAX
	ENDIF
	IF(ISF.LE.6) THEN
	MNF(ISF) = ISF
	ENDIF
	ENDDO


C	----------------------------
	IF(MCL.LE.0) THEN
	DEALLOCATE(MNF)
	RETURN
	ENDIF
C	----------------------------


	NEFO = 0  
	DO IS = 1,NSF
	ISF = MNF(IS)
	IF(ISF.GT.0) NEFO = NEFO + 1
	ENDDO


	ALLOCATE(B(NEFO,NEFO))
	B(1:NEFO,1:NEFO) = 0.0D0
	LCL(1:NEFO) = 0


C	----------------------------
C	TRANSFORMATION MATRIX
C	----------------------------
	NEFO = 0  
	NEFJ = NEFO
	IAX  = MLA
	DO IS = 1,NSF
	ISF = MNF(IS)

	IF(ISF.EQ.0) GOTO 150
	NEFO = NEFO + 1
	LCL(NEFO) = ISF
	NN   = NEFJ
	NN3  = 0
	DO JS = 1,NSF
	JSF = MNF(JS)
	IF(JSF.EQ.0) GOTO 120

	IF(JSF.LE.3) NN3 = NN3+1
	
	IF((ISF.GE.1.AND.ISF.LE.3).AND.(JSF.GE.1.AND.JSF.LE.3)) THEN
	NN = NN + 1
	ISF3 = ISF
	JSF3 = JSF
	IF(ISF3.GT.3) ISF3 = ISF3-3
	IF(JSF3.GT.3) JSF3 = JSF3-3
	IDM = ISF3 + 3*(JSF3-1)
	CALL RELFILL('$LAX',VAX,IDM,IAX,0)
	B(NEFO,NN) = VAX
	ENDIF

	IF((ISF.GE.4.AND.ISF.LE.6).AND.(JSF.GE.4.AND.JSF.LE.6)) THEN
	NN = NN + 1
	ISF3 = ISF
	JSF3 = JSF
	IF(ISF3.GT.3) ISF3 = ISF3-3
	IF(JSF3.GT.3) JSF3 = JSF3-3
	IDM = ISF3 + 3*(JSF3-1)
	CALL RELFILL('$LAX',VAX,IDM,IAX,0)
	B(NEFO,NN+NN3) = VAX
	ENDIF
		
120	CONTINUE
	ENDDO

150	CONTINUE
	ENDDO


C	----------------------------



C	WRITE(*,*) 'K1',IND,C(1:NSF)

C	PERFORM TRANSFORMATION
	CALL TRANCOL(B,C,NSF,NEFO,LCL,IND)  


	DEALLOCATE(MNF,B)
C	WRITE(*,*) 'K2',IND,C(1:NSF)
C	PAUSE

	RETURN

	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE TRANSYM1(MAXA,A,B,NEQ,MCL,LCL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	MATRIX MULTIPLY FOR UPPER TRAIANGULAR SYM. MATRIX  Akl = Bik*Aij*Bjl
C	=============================================================
	DIMENSION A(1),MAXA(1)
	DIMENSION B(MCL,MCL),LCL(MCL)

	DIMENSION LRO(MCL)

	ALLOCATABLE FR(:)
	MRO = MCL
	NUM = 2*(NEQ*(MCL+MRO) + MRO*MCL)
	ALLOCATE(FR(NUM))
	LRO(1:MCL) = LCL(1:MCL)

	JJ = 0
C	PERFORM   Aik = Aij*Bjk
	DO KEQ = 1,MCL
C	=====================================================
	DO IEQ = 1,NEQ
	R = 0.0D0
C	-----------------------------------------------------
	DO JE  = 1,MRO	!UPPER TRIANGULAR + DIAGONAL
	JEQ = LRO(JE)
	IF(JEQ.GE.IEQ) THEN
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1
	IF(LHEI.LE.IEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + A(NUM)*B(JEB,KEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	IF(IEQ.GT.JEQ) THEN
	NUM = MAXA(IEQ) + (IEQ - JEQ)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.JEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + B(JEB,KEQ)*A(NUM)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	ENDDO
C	-----------------------------------------------------
	KEB  = LCL(KEQ)
	IF(KEB.GE.IEQ) THEN
	NUM  = MAXA(KEB) + (KEB - IEQ)
	KHIGH = MAXA(KEB+1) - MAXA(KEB)
	LHEI = KEB - KHIGH + 1
	IF(LHEI.LE.IEQ) THEN
	JEB = MAPPING(LRO,IEQ,MRO)
C	IF(JEB.EQ.0) THEN
	JJ = JJ + 1
	FR(2*JJ-1)  = FLOAT(NUM)
	FR(2*JJ-0)  = R
C	ENDIF
	ENDIF
	ENDIF
C	-----------------------------------------------------
	ENDDO
C	=====================================================
	ENDDO



C	PERFORM   Aki = Bjk*Aji
	DO KEQ = 1,MCL
C	=====================================================
	DO IEQ =  1,NEQ
	R = 0.0D0
C	-----------------------------------------------------
	DO JE  =  1,MRO	!UPPER TRIANGULAR + DIAGONAL
	JEQ = LRO(JE)
	IF(IEQ.GE.JEQ) THEN
	NUM = MAXA(IEQ) + (IEQ - JEQ)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.JEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + B(JEB,KEQ)*A(NUM)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	IF(JEQ.GT.IEQ) THEN
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1
	IF(LHEI.LE.IEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + A(NUM)*B(JEB,KEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	ENDDO
C	-----------------------------------------------------
	KEB  = LCL(KEQ)
	IF(IEQ.GE.KEB) THEN
	NUM  = MAXA(IEQ) + (IEQ - KEB)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.KEB) THEN
	IEB = MAPPING(LCL,IEQ,MCL)
C	IF(IEB.EQ.0) THEN
	JJ = JJ + 1
	FR(2*JJ-1)  = FLOAT(NUM)
	FR(2*JJ-0)  = R
C	ENDIF
	ENDIF
	ENDIF
C	-----------------------------------------------------
	ENDDO
C	=====================================================
	ENDDO
	



C	PERFORM THE SYMETRICAL PART  Akl = Bik*Aij*Bjl           
	DO KEQ = 1,MCL
	DO LEQ = 1,MCL
	R = 0.0D0
C	=====================================================
	DO IE = 1,MRO
	IEQ = LRO(IE)
C	-----------------------------------------------------
	DO JE = 1,MCL	!UPPER TRIANGULAR + DIAGONAL
	JEQ = LCL(JE)
	IF(JEQ.GE.IEQ) THEN
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1
	IF(LHEI.LE.IEQ) R = R + B(IE,KEQ)*A(NUM)*B(JE,LEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	IF(JEQ.LT.IEQ) THEN
	NUM = MAXA(IEQ) + (IEQ - JEQ)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.JEQ) R = R + B(IE,KEQ)*A(NUM)*B(JE,LEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDDO
C	-----------------------------------------------------
	ENDDO
C	-----------------------------------------------------
	KEB  = LCL(KEQ)
	LEB  = LCL(LEQ)
	IF(LEB.GE.KEB) THEN
	NUM  = MAXA(LEB) + (LEB - KEB)
	KHIGH = MAXA(LEB+1) - MAXA(LEB)
	LHEI = LEB - KHIGH + 1
	IF(LHEI.LE.KEB) THEN
	JJ = JJ + 1
	FR(2*JJ-1)  = FLOAT(NUM)
	FR(2*JJ-0)  = R
	ENDIF
	ENDIF
C	=====================================================
	ENDDO
	ENDDO


C	FILL IN A
	DO II = 1,JJ
	NUM = INT(FR(2*II-1))
	R   = FR(2*II-0)
	A(NUM) = R
	ENDDO


	DEALLOCATE(FR)

	RETURN

	END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
	SUBROUTINE TRANSYM2(MAXA,A,B,C,NEQ,MCL,LCL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	MATRIX MULTIPLY FOR UPPER TRAIANGULAR SYM. MATRIX    
C	Akl = Bik*Aij*Bjl
C	Ck  = Bik*Ci
C	=============================================================
	DIMENSION A(1),C(1),MAXA(1)
	DIMENSION B(MCL,MCL),LCL(MCL)

	DIMENSION LRO(MCL),FC(MCL)

	ALLOCATABLE FR(:)
	MRO = MCL
	NUM = 2*(NEQ*(MCL+MRO) + MRO*MCL)
	ALLOCATE(FR(NUM))
	LRO(1:MCL) = LCL(1:MCL)

	JJ = 0
C	PERFORM   Aik = Aij*Bjk
	DO KEQ = 1,MCL
C	=====================================================
	DO IEQ = 1,NEQ
	R = 0.0D0
C	-----------------------------------------------------
	DO JE  = 1,MRO	!UPPER TRIANGULAR + DIAGONAL
	JEQ = LRO(JE)
	IF(JEQ.GE.IEQ) THEN
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1
	IF(LHEI.LE.IEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + A(NUM)*B(JEB,KEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	IF(IEQ.GT.JEQ) THEN
	NUM = MAXA(IEQ) + (IEQ - JEQ)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.JEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + B(JEB,KEQ)*A(NUM)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	ENDDO
C	-----------------------------------------------------
	KEB  = LCL(KEQ)
	IF(KEB.GE.IEQ) THEN
	NUM  = MAXA(KEB) + (KEB - IEQ)
	KHIGH = MAXA(KEB+1) - MAXA(KEB)
	LHEI = KEB - KHIGH + 1
	IF(LHEI.LE.IEQ) THEN
	JEB = MAPPING(LRO,IEQ,MRO)
C	IF(JEB.EQ.0) THEN
	JJ = JJ + 1
	FR(2*JJ-1)  = FLOAT(NUM)
	FR(2*JJ-0)  = R
C	ENDIF
	ENDIF
	ENDIF
C	-----------------------------------------------------
	ENDDO
C	=====================================================
	ENDDO



C	PERFORM   Aki = Bjk*Aji  Ck=Bik*Ci
	DO KEQ = 1,MCL
	FC(KEQ) = 0.0D0
C	=====================================================
	DO IEQ =  1,NEQ
	R = 0.0D0
C	-----------------------------------------------------
	DO JE  =  1,MRO	!UPPER TRIANGULAR + DIAGONAL
	JEQ = LRO(JE)
	IF(IEQ.GE.JEQ) THEN
	NUM = MAXA(IEQ) + (IEQ - JEQ)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.JEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + B(JEB,KEQ)*A(NUM)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	IF(JEQ.GT.IEQ) THEN
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1
	IF(LHEI.LE.IEQ) THEN
	JEB  = MAPPING(LRO,JEQ,MRO)
	IF(JEB.NE.0) R = R + A(NUM)*B(JEB,KEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDIF
	ENDDO
C	-----------------------------------------------------
	KEB  = LCL(KEQ)
	IF(IEQ.GE.KEB) THEN
	NUM  = MAXA(IEQ) + (IEQ - KEB)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.KEB) THEN
	IEB = MAPPING(LCL,IEQ,MCL)
C	IF(IEB.EQ.0) THEN
	JJ = JJ + 1
	FR(2*JJ-1)  = FLOAT(NUM)
	FR(2*JJ-0)  = R
C	ENDIF
	ENDIF
	ENDIF
C	-----------------------------------------------------
	IEB  = MAPPING(LRO,IEQ,MRO)
	IF(IEB.NE.0) FC(KEQ) = FC(KEQ) + B(IEB,KEQ)*C(IEQ)
C	-----------------------------------------------------
	ENDDO
C	=====================================================
	ENDDO
	
C	FILL IN C
	DO KEQ = 1,MCL
	KEB = LCL(KEQ)
	C(KEB) = FC(KEQ)
	ENDDO


C	PERFORM THE SYMETRICAL PART  Akl = Bik*Aij*Bjl           
	DO KEQ = 1,MCL
	DO LEQ = 1,MCL
	R = 0.0D0
C	=====================================================
	DO IE = 1,MRO
	IEQ = LRO(IE)
C	-----------------------------------------------------
	DO JE = 1,MCL	!UPPER TRIANGULAR + DIAGONAL
	JEQ = LCL(JE)
	IF(JEQ.GE.IEQ) THEN
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1
	IF(LHEI.LE.IEQ) R = R + B(IE,KEQ)*A(NUM)*B(JE,LEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	IF(JEQ.LT.IEQ) THEN
	NUM = MAXA(IEQ) + (IEQ - JEQ)
	KHIGH = MAXA(IEQ+1) - MAXA(IEQ)
	LHEI = IEQ - KHIGH + 1
	IF(LHEI.LE.JEQ) R = R + B(IE,KEQ)*A(NUM)*B(JE,LEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF
	ENDDO
C	-----------------------------------------------------
	ENDDO
C	-----------------------------------------------------
	KEB  = LCL(KEQ)
	LEB  = LCL(LEQ)
	IF(LEB.GE.KEB) THEN
	NUM  = MAXA(LEB) + (LEB - KEB)
	KHIGH = MAXA(LEB+1) - MAXA(LEB)
	LHEI = LEB - KHIGH + 1
	IF(LHEI.LE.KEB) THEN
	JJ = JJ + 1
	FR(2*JJ-1)  = FLOAT(NUM)
	FR(2*JJ-0)  = R
	ENDIF
	ENDIF
C	=====================================================
	ENDDO
	ENDDO


C	FILL IN A
	DO II = 1,JJ
	NUM = INT(FR(2*II-1))
	R   = FR(2*II-0)
	A(NUM) = R
	ENDDO

	DEALLOCATE(FR)


	RETURN

	END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

	SUBROUTINE TRANCOL(B,C,NEQ,MCL,LCL,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	MATRIX MULTIPLY FOR UPPER TRAIANGULAR SYM. MATRIX    
C	=============================================================
	DIMENSION C(1)
	DIMENSION B(MCL,MCL),LCL(MCL)
	DIMENSION LRO(MCL),FC(MCL)

	MRO = MCL
	LRO(1:MCL) = LCL(1:MCL)


	SELECTCASE(IND)
	CASE(0) 

C	PERFORM   Ck=Bik*Ci
	DO KEQ = 1,MCL
	FC(KEQ) = 0.0D0
C	=====================================================
	DO IEQ =  1,NEQ
C	-----------------------------------------------------
	IEB  = MAPPING(LRO,IEQ,MRO)
	IF(IEB.NE.0) FC(KEQ) = FC(KEQ) + B(IEB,KEQ)*C(IEQ)
C	-----------------------------------------------------
	ENDDO
C	=====================================================
	ENDDO

C	FILL IN C
	DO KEQ = 1,MCL
	KEB = LCL(KEQ)
	C(KEB) = FC(KEQ)
	ENDDO


	CASE(1)
C	PERFORM   Ck=Bki*Ci
	DO KEQ = 1,MRO
	FC(KEQ) = 0.0D0
C	=====================================================
	DO IEQ =  1,NEQ
C	-----------------------------------------------------
	IEB  = MAPPING(LCL,IEQ,MCL)
	IF(IEB.NE.0) FC(KEQ) = FC(KEQ) + B(KEQ,IEB)*C(IEQ)
C	-----------------------------------------------------
	ENDDO
C	=====================================================
	ENDDO

C	FILL IN C
	DO KEQ = 1,MRO
	KEB = LRO(KEQ)
	C(KEB) = FC(KEQ)
	ENDDO	


	ENDSELECT


	RETURN
	END

C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
