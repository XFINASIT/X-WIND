C	=======================================================================
C	=======================================================================
	SUBROUTINE SETLOD (NEQ,NDUM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
	COMMON /REACT/ LRC,LRCT,MFQ,LRID


	NDUM = MFQ*NEQ


C	-----------------------------------------------------------------------

	RETURN
	END

C	=======================================================================
C	=======================================================================
	SUBROUTINE INPSET(NSET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE
C	-----------------------------------------------------------------------
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	-----------------------------------------------------------------------
	NSETT = NSET

	IF(NSETT.NE.0) READ(ITI,*)

	CALL  DEMREL('SETL',KSETL,NSETT,10) 

	DO ISETT = 1,NSETT
	READ(ITI,*) NODS,IDIR,VALU,ILN,ILO,IFG 
	
	FNODS = FLOAT(NODS)
	FIDIR = FLOAT(IDIR)
	FILCN = FLOAT(ILN )
	FILCC = FLOAT(ILO )
	FIFG  = FLOAT(IFG )

	CALL MRELFIL('SETL',FNODS,ISETT,1,1)
	CALL MRELFIL('SETL',FIDIR,ISETT,2,1)
	CALL MRELFIL('SETL',VALU ,ISETT,3,1)
	CALL MRELFIL('SETL',FILCN,ISETT,4,1)
	CALL MRELFIL('SETL',FILCC,ISETT,5,1)
	CALL MRELFIL('SETL',FIFG ,ISETT,6,1)

	ENDDO
	
C	-----------------------------------------------------------------------

	RETURN
	END

C	=======================================================================
C	=======================================================================
	SUBROUTINE INPSETOFF(NSET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE FOR OFFSHORE
C	-----------------------------------------------------------------------
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C	-----------------------------------------------------------------------
	NSETT = NSET

	IF(NSETT.NE.0) READ(ITI,*)

	CALL  DEMREL('SETF',KSETL,NSETT,10) 

	DO ISETT = 1,NSETT
	READ(ITI,*) NODS,IDIR,VALU,ILN,ILO,IFG 
	
	FNODS = FLOAT(NODS)
	FIDIR = FLOAT(IDIR)
	FILCN = FLOAT(ILN )
	FILCC = FLOAT(ILO )
	FIFG  = FLOAT(IFG )

	CALL MRELFIL('SETF',FNODS,ISETT,1,1)
	CALL MRELFIL('SETF',FIDIR,ISETT,2,1)
	CALL MRELFIL('SETF',VALU ,ISETT,3,1)
	CALL MRELFIL('SETF',FILCN,ISETT,4,1)
	CALL MRELFIL('SETF',FILCC,ISETT,5,1)
	CALL MRELFIL('SETF',FIFG ,ISETT,6,1)

	ENDDO
	
C	-----------------------------------------------------------------------

	RETURN
	END

C	=======================================================================
C	=======================================================================
	SUBROUTINE SETLAD(SETIF,VEC,ILC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL
C	-----------------------------------------------------------------------
	DIMENSION SETIF(NEQ,1)
C	-----------------------------------------------------------------------
	DIMENSION VEC(NEQ)

	IF(ILC.LE.0) RETURN

	DO I = 1,NEQ
	VEC(I) = VEC(I) + SETIF(I,ILC)
	ENDDO

C	-----------------------------------------------------------------------

	RETURN
	END

C	=======================================================================
C	=======================================================================
	SUBROUTINE SETLODR(SETIF,ID,IDRCT,LM,LMRCT,S,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	-----------------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL
C	-----------------------------------------------------------------------
	COMMON /REACT/ LRC,LRCT,MFQ,LRID
	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE
	COMMON /STCAS/ ILC
C	-----------------------------------------------------------------------
	DIMENSION SETIF(NEQ,LCS),IDRCT(NSF,NSN),ESTF(NEF,NEF),S(1)
	DIMENSION ID(NSF,NSN)
C	-----------------------------------------------------------------------
	DIMENSION LM(NEF),LMRCT(NEF)

      IF(NSETT.LE.0) RETURN
      
C	CALLING THE ELEMENT STIFFNESS
	ESTF = 0.0
	K = 0
	DO I = 1,NEF
	DO J = I,NEF
	K = K+1
	ESTF(I,J) = S(K)
	ESTF(J,I) = S(K)
	ENDDO
	ENDDO

	IF(LSYMM.EQ.1) THEN
	DO I = 1,NEF
	DO J = I,NEF
	K = K+1
	ESTF(J,I) = S(K)
	ENDDO
	ENDDO
	ENDIF

C	----------------------------------------------------------------------
	DO 100 II = 1,NSETT               !LOOP OVER THE NODAL SETTLEMENT DATA
C	----------------------------------------------------------------------
	CALL MRELFIL('SETL',FNODS,II,1,0) !CALLING THE SETTLEMENT DATA
	CALL MRELFIL('SETL',FIDIR,II,2,0)
	CALL MRELFIL('SETL',VALU ,II,3,0)
	CALL MRELFIL('SETL',FILCN,II,4,0)
	CALL MRELFIL('SETL',FILCC,II,5,0)
	CALL MRELFIL('SETL',FIFG ,II,6,0)
	NODS = INT(FNODS)
	IDIR = INT(FIDIR)
	ILN  = INT(FILCN)
	ILO  = INT(FILCC)
	IFG  = INT(FIFG )

	IDIR = IDOFCALL(IDOF,IDIR)
	IFQ  = IDRCT(IDIR,NODS)
	IKQ  =    ID(IDIR,NODS)

	IF(IFG.NE.0  ) GOTO 100

	DO IEF = 1,NEF                !STORE THE SETTLEMENT LOAD INTO ARRAY
	IEQ = LM(IEF)
	IF(IEQ.EQ.0  ) GOTO 50
	DO JEF = 1,NEF
	JFQ = LMRCT(JEF)
	JEQ =    LM(JEF)

	IF(JFQ.EQ.IFQ.AND.JFQ.NE.0) THEN
	IF(ILN.NE.0) SETIF(IEQ,ILN) = SETIF(IEQ,ILN) - ESTF(IEF,JEF)*VALU
	IF(ILO.NE.0) SETIF(IEQ,ILO) = SETIF(IEQ,ILO) - ESTF(IEF,JEF)*VALU
	ENDIF

	IF(JEQ.EQ.IKQ.AND.JEQ.NE.0) THEN
	IF(ILN.NE.0) SETIF(IEQ,ILN) = SETIF(IEQ,ILN) - ESTF(IEF,JEF)*VALU
	IF(ILO.NE.0) SETIF(IEQ,ILO) = SETIF(IEQ,ILO) - ESTF(IEF,JEF)*VALU
	ENDIF

	ENDDO
50	CONTINUE
	ENDDO
C	----------------------------------------------------------------------
100	CONTINUE                         !END LOOP OVER THE NODAL SETTLEMENT DATA
C	----------------------------------------------------------------------


C	----------------------------------------------------------------------
	DO 200 II = 1,NSETT               !LOOP OVER THE NODAL SETTLEMENT DATA
C	----------------------------------------------------------------------
	CALL MRELFIL('SETL',FNODS,II,1,0) !CALLING THE SETTLEMENT DATA
	CALL MRELFIL('SETL',FIDIR,II,2,0)
	CALL MRELFIL('SETL',VALU ,II,3,0)
	CALL MRELFIL('SETL',FILCN,II,4,0)
	CALL MRELFIL('SETL',FILCC,II,5,0)
	CALL MRELFIL('SETL',FIFG ,II,6,0)
	NODS = INT(FNODS)
	IDIR = INT(FIDIR)
	ILN  = INT(FILCN)
	ILO  = INT(FILCC)
	IFG  = INT(FIFG )

	IDIR = IDOFCALL(IDOF,IDIR)
	IFQ  = IDRCT(IDIR,NODS)
	IKQ  =    ID(IDIR,NODS)

	IF(ILN.NE.ILC.AND.ILO.NE.ILC) GOTO 200
	IF(IFG.NE.0  ) GOTO 200
	IF(IKQ.EQ.0  ) GOTO 200

	DO IEF = 1,NEF                
	IEQ = LM(IEF)

	IF(IEQ.NE.0.AND.IEQ.EQ.IKQ) THEN
	ESTF(1:IEF,IEF)   = 0.0D0               !FILL IN ZERO
	ESTF(IEF,IEF:NEF) = 0.0D0
	ENDIF

	ENDDO
C	----------------------------------------------------------------------
200	CONTINUE                        !END LOOP OVER THE NODAL SETTLEMENT DATA
C	----------------------------------------------------------------------
	K = 0
	DO I = 1,NEF
	DO J = I,NEF
	K = K+1
	S(K) = ESTF(I,J)
	ENDDO
	ENDDO

	IF(LSYMM.EQ.1) THEN
	DO I = 1,NEF
	DO J = I,NEF
	K = K+1
	S(K) = ESTF(J,I)
	ENDDO
	ENDDO
	ENDIF


	RETURN
	END

C	=======================================================================
C	=======================================================================
	SUBROUTINE SETELM (EDIS,LMRCT,IDRCT,NEF,LM,ID)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL				!GETTING NEQ
	COMMON /REACT/ LRC,LRCT,MFQ,LRID							!GETTING MFQ

	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE
	COMMON /STCAS/ ILC
      
      ALLOCATABLE SETFIX(:),SETFRE(:)
C	-----------------------------------------------------------------------
	DIMENSION LMRCT(NEF),IDRCT(NSF,NSN),EDIS(NEF)
	DIMENSION LM(NEF),ID(NSF,NSN) ! FOR HEAT ANALYSIS MAR2007
C	-----------------------------------------------------------------------

	IF(ILC.EQ.0) RETURN

      
      IF(NSETT.GT.0) THEN
          
	ALLOCATE(SETFIX(MFQ),SETFRE(NEQ))
	SETFIX(1:MFQ) = 0.0D0
	SETFRE(1:NEQ) = 0.0D0

	DO II = 1,NSETT

	CALL MRELFIL('SETL',FNODS,II,1,0) !CALLING THE SETTLEMENT DATA
	CALL MRELFIL('SETL',FIDIR,II,2,0)
	CALL MRELFIL('SETL',VALU ,II,3,0)
	CALL MRELFIL('SETL',FILCN,II,4,0)
	CALL MRELFIL('SETL',FILCC,II,5,0)
	CALL MRELFIL('SETL',FIFG ,II,6,0)
	NODS = INT(FNODS)
	IDIR = INT(FIDIR)
	ILN  = INT(FILCN)
	ILO  = INT(FILCC)
	IFG  = INT(FIFG )

	IF(ILN.NE.ILC.AND.ILO.NE.ILC) GOTO 10
	IF(IFG.NE.0  ) GOTO 10

	IDIR = IDOFCALL(IDOF,IDIR)
	JFQ  = IDRCT(IDIR,NODS)
	JEQ  =    ID(IDIR,NODS)
	
	IF(JFQ.NE.0) SETFIX(JFQ) = SETFIX(JFQ) + VALU
	IF(JEQ.NE.0) SETFRE(JEQ) = SETFRE(JEQ) + VALU

10	CONTINUE

	ENDDO

	DO IDF = 1,NEF
	    IFQ = LMRCT(IDF)
	    IEQ = LM(IDF)
	
	    IF(IFQ.NE.0.AND.SETFIX(IFQ).NE.0.0D0) EDIS(IDF) = SETFIX(IFQ)
	    IF(IEQ.NE.0.AND.SETFRE(IEQ).NE.0.0D0) EDIS(IDF) = SETFRE(IEQ)
      ENDDO
          
      DEALLOCATE(SETFIX,SETFRE)
          
      ENDIF !ENDIF NSETT > 0



	RETURN
	END

C	=======================================================================
C	=======================================================================
	SUBROUTINE SETPRN (D,IDRCT,ID,ISN,ISF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE
	COMMON /STCAS/ ILC
C	-----------------------------------------------------------------------
	DIMENSION ID(NSF,NSN),IDRCT(NSF,NSN),D(9)
C	-----------------------------------------------------------------------

	IF(ILC.EQ.0) RETURN

	IFQ = IDRCT(ISF,ISN)
	IEQ =    ID(ISF,ISN)

	DO 100 II = 1,NSETT

	CALL MRELFIL('SETL',FNODS,II,1,0) !CALLING THE SETTLEMENT DATA
	CALL MRELFIL('SETL',FIDIR,II,2,0)
	CALL MRELFIL('SETL',VALU ,II,3,0)
	CALL MRELFIL('SETL',FILCN,II,4,0)
	CALL MRELFIL('SETL',FILCC,II,5,0)
	CALL MRELFIL('SETL',FIFG ,II,6,0)
	NODS = INT(FNODS)
	IDIR = INT(FIDIR)
	ILN  = INT(FILCN)
	ILO  = INT(FILCC)
	IFG  = INT(FIFG )

	IDIR = IDOFCALL(IDOF,IDIR)
	JFQ  = IDRCT(IDIR,NODS)
	JEQ  =    ID(IDIR,NODS)

	VAV = 0.0D0
	IF(ILN.EQ.ILC) VAV = VAV + VALU
	IF(ILO.EQ.ILC) VAV = VAV + VALU
	IF(ILN.NE.ILC.AND.ILO.NE.ILC) GOTO 100

	IF(JFQ.EQ.IFQ.AND.JFQ.NE.0.AND.IFG.EQ.0) THEN
	LL = IDOF(ISF)
	D(LL) = VAV
	ENDIF
	IF(JEQ.EQ.IEQ.AND.JEQ.NE.0.AND.IFG.EQ.0) THEN
	LL = IDOF(ISF)
	D(LL) = VAV
	ENDIF

100	CONTINUE




	RETURN
	END

C	=======================================================================
C	=======================================================================
	SUBROUTINE INTLAD(ID,DIS,ILC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-----------------------------------------------------------------------
C	-----------------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL
C	-----------------------------------------------------------------------
	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE
C	-----------------------------------------------------------------------
	DIMENSION ID(NSF,NSN),DIS(NEQ)
C	-----------------------------------------------------------------------

	IF(ILC.EQ.0) RETURN

	DO 100 II = 1,NSETT
	
	CALL MRELFIL('SETL',FNODS,II,1,0) !CALLING THE SETTLEMENT DATA
	CALL MRELFIL('SETL',FIDIR,II,2,0)
	CALL MRELFIL('SETL',VALU ,II,3,0)
	CALL MRELFIL('SETL',FILCN,II,4,0)
	CALL MRELFIL('SETL',FILCC,II,5,0)
	CALL MRELFIL('SETL',FIFG ,II,6,0)
	NODS = INT(FNODS)
	IDIR = INT(FIDIR)
	ILN  = INT(FILCN)
	ILO  = INT(FILCC)
	IFG  = INT(FIFG )

	IDIR = IDOFCALL(IDOF,IDIR)
	IEQ  = ID(IDIR,NODS)

	VAV = 0.0D0
	IF(ILN.EQ.ILC) VAV = VAV + VALU
	IF(ILO.EQ.ILC) VAV = VAV + VALU
	IF(ILN.NE.ILC.AND.ILO.NE.ILC) GOTO 100

	IF(IFG.EQ.1.AND.IEQ.NE.0) THEN
	DIS(IEQ) = VALU
	ENDIF

100	CONTINUE
C	-----------------------------------------------------------------------

	RETURN
	END

C	=======================================================================
C	=======================================================================
	FUNCTION IDOFCALL(IDOF,NUM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION IDOF(9)

	DO II = 1,9
	IF(NUM.EQ.IDOF(II)) THEN
	IDOFCALL = II
	RETURN
	ENDIF
	ENDDO

	RETURN
	END

C	=======================================================================
C	=======================================================================
