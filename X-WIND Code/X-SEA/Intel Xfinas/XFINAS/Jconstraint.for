C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE JONSTRN(ID,IGOF,NSF,NSN,NEQ,KEQ,LCON1,LCON2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	==================================================================
C	JOINT GLOBAL CONSTRAINT SONGSAK JAN2007 
	COMMON /LC007/ LC107,LC207,LC307,LC407
C	==================================================================
C	LC107 = NUMBER CONSTRAINT GROUP
C	LC207 = MAXIMUM NUMBER OF NODE IN ONE GROUP
C	LC307 = ARRAY POINTOR FOR LCON1
C	LC407 = ARRAY POINTOR FOR LCON2
C	LCON1 = CONSTRAINT DATA (JGP,JND,ICOF)
C	LCON1 = CONSTRAINT DATA (NODE NUMBER)
C	==================================================================
	DIMENSION JCOUT(NEQ),ICOF(9),NODE(LC207),JCGRP(NEQ,LC107)
	DIMENSION ID(NSF,NSN),IGOF(9),IDQ(NSF,NSN),JCONT(NEQ),KAK(NEQ)
	DIMENSION LACTV(NEQ)
	DIMENSION LCON1(11,LC107),LCON2(LC207,LC107)

	NCX   = LC107
	MAXCN = LC207


	JCGRP = 0
	IDQ   = 0
C	----------------------------------------------------------------
	NCT = 0
	DO I = 1,NCX

	JGP = LCON1(1,I)                 !GROUP NO.
	JND = LCON1(2,I)                 !NUMBER OF NODE IN THIS GROUP
	ICOF(1:9)   = LCON1(3:11 ,I)     !CONSTRAINT DOF (1 = CONSTRAINT)
	NODE(1:JND) = LCON2(1:JND,I)     !NODE NUMBER IN THIS GROUP

	ISF = 0
	DO K = 1,9 
	IG  = IGOF(K)
	IF(IG.EQ.0) GOTO 10
	ISF = ISF + 1 
	KCP = ICOF(IG)
	IF(KCP.EQ.1) NCT = NCT + 1
	DO J = 1,JND
	ND  = NODE(J)
	IEQ = ID(ISF,ND)
	IF(KCP.EQ.1.AND.IEQ.NE.0) JCGRP(IEQ,I) = NCT
	ENDDO
10	CONTINUE
	ENDDO
	
	ENDDO

C	WRITE(*,*) NCX,MAXCN
C	DO I = 1,NEQ
C	WRITE(110,*) JCGRP(I,1:NCX)
C	ENDDO
C	PAUSE
C	WRITE(*,*) JCGRP(1:NEQ,1)
	JCOUT = 0
C	----------------------------------------------------------------
	DO I = 1,NCX
	DO J = 1,NEQ
	IF(JCGRP(J,I).NE.0.AND.JCOUT(J).EQ.0) JCOUT(J) = JCGRP(J,I)
	ENDDO
	ENDDO	

C	DO I = 1,NEQ
C	WRITE(110,*) JCOUT(I)
C	ENDDO
C	PAUSE
C	----------------------------------------------------------------
	DO I = 1,NCX

	JGP = LCON1(1,I)
	JND = LCON1(2,I)
	ICOF(1:9)   = LCON1(3:11 ,I)
	NODE(1:JND) = LCON2(1:JND,I)

	ISF = 0
	DO K = 1,9 !NSF
	IG  = IGOF(K)
	IF(IG.EQ.0) GOTO 20

	ISF = ISF + 1 
	KCP = ICOF(IG)

	MIN = NCT
	DO J = 1,JND
	ND  = NODE(J)
	IEQ = ID(ISF,ND)
	IF(KCP.EQ.1.AND.IEQ.NE.0) THEN
	IF(MIN.GT.JCOUT(IEQ)) MIN = JCOUT(IEQ)
	ENDIF
	ENDDO

	DO J = 1,JND
	ND  = NODE(J)
	IEQ = ID(ISF,ND)
	IF(KCP.EQ.1.AND.IEQ.NE.0) THEN
	JCOUT(IEQ) = MIN
	ENDIF
	ENDDO

20	CONTINUE

	ENDDO
	
	ENDDO

C	----------------------------------------------------------------
C	THIS BLOCK ADD TO FIND THE CONSTRAINT THAT FIX TO INACTIVE DOF SUCH AS SUPPORT
	LACTV(1:NEQ) = 0
	DO I = 1,NCX

	JGP = LCON1(1,I)
	JND = LCON1(2,I)
	ICOF(1:9)   = LCON1(3:11 ,I)
	NODE(1:JND) = LCON2(1:JND,I)

	ISF = 0
	DO K = 1,9 !NSF
	IG  = IGOF(K)
	IF(IG.EQ.0) GOTO 30

	ISF = ISF + 1 
	KCP = ICOF(IG)

C	SCAN FOR ACTIVE DOF 
	IEQA = 0
	DO J = 1,JND
	ND  = NODE(J)
	IEQ = ID(ISF,ND)
	IF(KCP.EQ.1.AND.IEQ.NE.0) THEN
	IEQA = IEQ  !PICK ACTIVE DOF
	EXIT
	ENDIF
	ENDDO

C	SCAN WHETHER THERE IS INACTIVE DOF (SUPPORT) INSIDE THIS GROUP OR NOT
	IACTV = 0 
	DO J = 1,JND
	ND  = NODE(J)
	IEQ = ID(ISF,ND)
	IF(KCP.EQ.1.AND.IEQ.EQ.0) THEN
	IACTV = 1
	EXIT
	ENDIF
	ENDDO


	IF(IEQA.NE.0.AND.IACTV.EQ.1) THEN
	NUM = JCOUT(IEQA)
	DO IQ = 1,NEQ
	IF(JCOUT(IQ).EQ.NUM) LACTV(IQ) = 1 !SET DOF IN CURRENT CONSTRAINT GROUP TO INACTIVE DUE TO CONSTRAINT TO SUPPORT
	ENDDO
	ENDIF


30	CONTINUE

	ENDDO
	
	ENDDO	

C	----------------------------------------------------------------

	MAK = 0
	KAK = 0
	JCONT = 0
	DO I = 1,NEQ
	KI = JCOUT(I)
	IACTV = LACTV(I)
	IF(IACTV.EQ.1) GOTO 40  !IF THIS DOF IS INACTIVE GOTO 40

	CALL CSEARH(KI,KAK,NEQ,IND)
	IF(IND.EQ.1) THEN

	MAK = MAK + 1
	IF(KI.NE.0) THEN
	DO J = 1,NEQ
	KJ = JCOUT(J)
	IF(KJ.EQ.KI) JCONT(J) = MAK
	ENDDO
	ENDIF

	ENDIF

40	CONTINUE
	ENDDO
	JCOUT(1:NEQ) = JCONT(1:NEQ)
C	----------------------------------------------------------------

	KEQ = 0
	DO I = 1,NSN
	DO J = 1,NSF

	IEQ = ID(J,I)
	IF(IEQ.NE.0) THEN
	KCP = JCOUT(IEQ)
	IACTV = LACTV(IEQ)
	IF(KCP.EQ.0.AND.IACTV.EQ.0) KEQ = KEQ + 1    !FOR ACTIVE DOF
	IF(KCP.EQ.0.AND.IACTV.EQ.0) IDQ(J,I) = KEQ   !FOR ACTIVE DOF
	ENDIF
	ENDDO
	ENDDO

C	----------------------------------------------------------------
	MAQ = 0
	DO I = 1,NSN
	DO J = 1,NSF
	IEQ = ID(J,I)
	IF(IEQ.NE.0) THEN
	KCP = JCOUT(IEQ)
	IF(KCP.GT.MAQ) THEN
	MAQ = KCP
	ENDIF
	IF(KCP.NE.0) IDQ(J,I) = KEQ + KCP
	ENDIF
	ID(J,I) = IDQ(J,I)
	ENDDO
	ENDDO
	KEQ = KEQ + MAQ
C	----------------------------------------------------------------
	

	RETURN

	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE CSEARH(KI,KAK,NEQ,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION KAK(NEQ)

	IND = 1
	DO I = 1,NEQ
	IF(KI.EQ.KAK(I)) IND = 0
	ENDDO
	
	DO I = 1,NEQ
	IF(KAK(I).EQ.0) THEN
	KAK(I) = KI
	RETURN
	ENDIF
	ENDDO


	RETURN

	END
C	=============================================================
C	=============================================================
C	=============================================================

