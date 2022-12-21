C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE ELMEATV(IATV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      DIMENSION IATV(1)
      DIMENSION NPM(10),NPI(10)
	ALLOCATABLE IF2(:),IF3(:)
	ALLOCATABLE IAX(:),AAX(:)


C	CALLING NEG
	CALL INTFILL('%NUB',NEG,1,5,0)

C	-----------------------------------------
C	LOOP OVER GROUP
C	-----------------------------------------
	DO 1600 IEG = 1,NEG

C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('%GRP',NELE ,3 ,IEG,0) !
	CALL INTFILL('%GRP',NEF  ,7 ,IEG,0) !

C	GROUP FILE
	CALL INTFILL('%GRF',N2   ,2 ,IEG,0) !
	CALL INTFILL('%GRF',N3   ,3 ,IEG,0) !

	CALL INTFILL('%GRF',NDX2 ,12,IEG,0) !
	CALL INTFILL('%GRF',NDX3 ,13,IEG,0) !
	
	CALL INTFILL('%GRF',NFLI ,21,IEG,0) !
	CALL INTFILL('%GRF',NFLR ,22,IEG,0) !
	
	CALL INTFILL('%GRF',LASI ,23,IEG,0) !
	CALL INTFILL('%GRF',LASR ,24,IEG,0) !


	ALLOCATE(IF2(N2),IF3(N3))
	ALLOCATE(IAX(LASI),AAX(LASR))
	
C	-----------------------------------------
C	LOOP OVER ELEMENT
C	-----------------------------------------
	DO 1500 IEL = 1,NELE

      READ(NFLI,REC=NELE+IEL,ERR=5) IXIST !READ EXISTING FLAG TO LAST RECORD OF ELEMENT INTEGER DATA (INCREASE THE SPEED PURPOSE)
      GOTO 6
5     IXIST = 0
6     IF(IXIST.LE.0) GOTO 1500

      READ(NFLI,REC=IEL) IAX
      READ(NFLR,REC=IEL) AAX

	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(IF3,IAX,NDX3,NPM,NPI,1,'READ')
	IXIST = IF3(4)  !TOTAL EXISTING FLAG 
	IF (IXIST.EQ. 0) GOTO 1500 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED

	NPM(1) = N2
	NPI(1) = N2
	CALL CALLINTG(IF2,IAX,NDX2,NPM,NPI,1,'READ')              !INTEGER DATA


C	-----------------------------------------
C	LOOP OVER ELEMENT DOF
C	-----------------------------------------
	DO 1000 IEF = 1,NEF
	    JEQ = IF2(IEF) !ORIGINAL EQUATION NUMBER
	    IF(JEQ.GT.0) THEN
            IATV(JEQ) = 0    !FREE THAT DOF IF IT OCCUPIED BY ACTIVE ELEMENT
	    ENDIF
1000	CONTINUE
C	-----------------------------------------


C	-----------------------------------------
1500	CONTINUE
C	-----------------------------------------

	DEALLOCATE(IAX,AAX)
	DEALLOCATE(IF2,IF3)
	
C	-----------------------------------------
1600	CONTINUE
C	-----------------------------------------



	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE TDNEATV(IATV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /TENDONC/ NEGT  !NUMBER OF TENDON GROUP
	
      DIMENSION IATV(1)
      DIMENSION NPM(10),NPI(10)
	ALLOCATABLE IF2(:),AF3(:)
	ALLOCATABLE IAX(:),AAX(:)


C	-----------------------------------------
C	LOOP OVER GROUP
C	-----------------------------------------

      DO 8000 IEGT = 1,NEGT
      
C	----------------------------------------
C	GROUP DATA
C	----------------------------------------
	CALL INTFILL('TGRP',NTEND,1 ,IEGT,0) !
	CALL INTFILL('TGRP',MSEG ,2 ,IEGT,0) !MAX NUM OF TENDON SEGMENT
	CALL INTFILL('TGRP',NEFT ,3 ,IEGT,0) !
	CALL INTFILL('TGRP',NNMT ,4 ,IEGT,0) !
	CALL INTFILL('TGRP',NNFT ,5 ,IEGT,0) !
	CALL INTFILL('TGRP',LINK ,6 ,IEGT,0) !
	CALL INTFILL('TGRP',ITDYP,7 ,IEGT,0) !ITDYP = TENDON OPTION (0=ON NODE, OTHERWISE=ON ELEM SEPECIFY BY ITDYP)	
C	GROUP FILE	
	CALL INTFILL('TGRF',N2   ,2 ,IEGT,0) !
	CALL INTFILL('TGRF',N3   ,3 ,IEGT,0) !

	CALL INTFILL('TGRF',NDX2 ,12,IEGT,0) !
	CALL INTFILL('TGRF',NDX3 ,13,IEGT,0) !
	
	CALL INTFILL('TGRF',NFLI ,21,IEGT,0) !
	CALL INTFILL('TGRF',NFLR ,22,IEGT,0) !
	
	CALL INTFILL('TGRF',LASI ,23,IEGT,0) !
	CALL INTFILL('TGRF',LASR ,24,IEGT,0) ! 
	
C	-----------------------------------------

	ALLOCATE(IF2(N2),AF3(N3))
	ALLOCATE(IAX(LASI),AAX(LASR))
	
C	-----------------------------------------
C	LOOP OVER TENDON ELEMENT
C	-----------------------------------------
	DO 1600 ITEND = 1,NTEND

      READ(NFLI,REC=NTEND+ITEND,ERR=5) IXIST !READ EXISTING FLAG TO LAST RECORD OF ELEMENT INTEGER DATA (INCREASE THE SPEED PURPOSE)
      GOTO 6
5     IXIST = 0
6     IF(IXIST.LE.0) GOTO 1600

      READ(NFLI,REC=ITEND) IAX
      READ(NFLR,REC=ITEND) AAX

	NPM(1) = N3
	NPI(1) = N3
	CALL CALLINTG(AF3,AAX,NDX3,NPM,NPI,1,'READ')
	IXIST = INT(AF3(5))  !EXISTING FLAG 
	NSEG  = INT(AF3(4))  !Number of tendon segment
	IF (IXIST.EQ. 0) GOTO 1600 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED

C	-----------------------------------------
C	LOOP OVER TENDON SEGMENT
C	-----------------------------------------
	DO 1500 ISEG = 1,NSEG

	NPM(1:2) = [N2,MSEG]
	NPI(1:2) = [N2,ISEG]
	CALL CALLINTG(IF2,IAX,NDX2,NPM,NPI,2,'READ')              !INTEGER DATA

C	-----------------------------------------
C	LOOP OVER ELEMENT DOF
C	-----------------------------------------

	DO 1000 IEFT = 1,NEFT

	    JEQ = IF2(IEFT) !ORIGINAL EQUATION NUMBER
	    IF(JEQ.GT.0) THEN
            IATV(JEQ) = 0    !FREE THAT DOF IF IT OCCUPIED BY ACTIVE ELEMENT
	    ENDIF

C	-----------------------------------------
1000	CONTINUE
C	-----------------------------------------

C	-----------------------------------------
1500	CONTINUE
C	-----------------------------------------

C	-----------------------------------------
1600	CONTINUE
C	-----------------------------------------

	DEALLOCATE(IAX,AAX)
	DEALLOCATE(IF2,AF3)
	
8000  CONTINUE


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE LNKEATV(IATV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      DIMENSION IATV(1)


      CALL LOCATN ('+LNK',KLNK,LNKWRK,NELE,2) !CALLING NELE
C	-----------------------------------------
C	LOOP OVER TENDON ELEMENT
C	-----------------------------------------
	DO 1600 IEL = 1,NELE

	CALL RELFILL('+LNK',XIST,4,IEL,0)
	IXIST = INT(XIST)

	IF (IXIST.EQ. 0) GOTO 1600 !ELEMENT HAS NEVER EXISTED & CURRENTLY REMOVED

	CALL RELFILL('+LNK',FNEF,6,IEL,0)  !NEF
	NEF = INT(FNEF)

C	-----------------------------------------
C	LOOP OVER ELEMENT DOF
C	-----------------------------------------
	DO 1000 IEF = 1,NEF
	    NFI = 10 + IEF
	    CALL RELFILL('+LNK',FNFI,NFI,IEL,0)
	    JEQ = INT(FNFI) !ORIGINAL EQUATION NUMBER
	    IF(JEQ.GT.0) THEN
            IATV(JEQ) = 0    !FREE THAT DOF IF IT OCCUPIED BY ACTIVE ELEMENT
	    ENDIF
C	-----------------------------------------
1000	CONTINUE
C	-----------------------------------------


C	-----------------------------------------
1600	CONTINUE
C	-----------------------------------------


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE TRVEATV(IATV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      DIMENSION IATV(1)



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

C	-----------------------------------------
C	LOOP OVER ELEMENT DOF
C	-----------------------------------------
	DO 1000 IEF = 1,NEF
	    NFI = 10 + NNM + IEF 
	    CALL RELFILL('+TRV',FNFI,NFI,IEL,0)
	    JEQ = INT(FNFI) !ORIGINAL EQUATION NUMBER
	    IF(JEQ.GT.0) THEN
            IATV(JEQ) = 0    !FREE THAT DOF IF IT OCCUPIED BY ACTIVE ELEMENT
	    ENDIF
C	-----------------------------------------
1000	CONTINUE
C	-----------------------------------------

C	-----------------------------------------
1600	CONTINUE
C	-----------------------------------------


	RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== SONGSAK NOV2007 =============
C	=======================================================================

