C	=================================================================
C	=================================================================
C	=================================================================

	SUBROUTINE PARLWRIT(LP,KREC,AOUT,NUMDAT,LERR,OPER)   
      USE OMP_LIB
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPER
      REAL*8 AOUT(NUMDAT)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      LERR = 0
      
      SELECTCASE(OPER)
      CASE('READ')
      

!$OMP PARALLEL SHARED(LP,KREC,AOUT,NUMDAT,LERR)
	IF(LP.GT.0)  READ(LP,REC=KREC,ERR=10) AOUT
	GOTO 20
10    LERR = 1
20    CONTINUE
!$OMP END PARALLEL
      GOTO 1000

      CASE('WRIT')

!$OMP PARALLEL SHARED(LP,KREC,AOUT,NUMDAT,LERR)
	IF(LP.GT.0) WRITE(LP,REC=KREC,ERR=30) AOUT
	GOTO 40
30    LERR = 2
40    CONTINUE
!$OMP END PARALLEL
      GOTO 1000

      ENDSELECT
      
      
1000  CONTINUE

      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
