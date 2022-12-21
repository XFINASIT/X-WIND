C     =========================================================
C     =========================================================
C     =========================================================
        SUBROUTINE PROGBAR(I,IMAX)
	  IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        COMMON /PRGBAR/ TIMES
        CHARACTER(LEN=1) :: BAR, BACK
        
        RETURN
        
        BACK = CHAR(8) !BACK SPACE
        BAR  = '='
        BARX = '>'
        BARD = '.'

        ITO = 6
        

        IF(I.EQ.1) THEN
            TIMES = SECOND() 
            CALL PRINBAR()
        ENDIF
        
        TIMER = 0.0D0
        IF(TIMES.NE.0.0D0) THEN
            TIMEN = SECOND() 
            TIMER = (TIMEN-TIMES)*(FLOAT(IMAX-I)/FLOAT(I))/60.0 !REMAINING TIME IN MINUTE
        ENDIF
        
        CALL DELTBAR()
        CALL PRINBAR()
C       CALL SLEEP(1) !PAUSE FOR A GIVE SECOND
        IF(I.EQ.IMAX) THEN
            TIMEN = SECOND() 
            TIMEE = (TIMEN - TIMES)/60.0
            WRITE(ITO,'(A)') '| complete'
            WRITE(ITO,'(2X,1A12,1F7.2,1A3)') 'TIME ELAPSE=',TIMEE,'MIN'
            TIMES = 0.0D0
        ENDIF

        
        CONTAINS
        
        SUBROUTINE PRINBAR
C     print the percentage and the bar
C        WRITE(ITO,'(2X,1I3,1A1,2X,1A1,256A1)',ADVANCE='NO')
C     1              100*I/IMAX,'%','|', (BAR, K=1,50*I/IMAX)
        WRITE(ITO,'(2X,1A10,1F7.2,1A3,2X,1I3,1A1,2X,1A1,38A1)',
     1              ADVANCE='NO')
     1              'TIME LEFT=',TIMER,'MIN',100*I/IMAX,'%','|', 
     2     (BAR, K=1,38*(I-1)/IMAX), BARX, (BARD, K=1,38*(IMAX-I)/IMAX)
        CLOSE(ITO)
        OPEN(ITO)
        END SUBROUTINE PRINBAR
        
        SUBROUTINE DELTBAR
C     delete the bar and the percentage
        WRITE(ITO,'(69A1)',ADVANCE='NO') (BACK, K=1,38+31) !9 = THERE ARE 9 CHARACTER BEFORE REACH PROGRESS BAR (FROM LEFT)
        END SUBROUTINE DELTBAR
        
        END SUBROUTINE PROGBAR
C     =========================================================
C     =========================================================
C     =========================================================
        SUBROUTINE PROGBAR_O(I,IMAX)
	  IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        COMMON /PRGBAR/ TIMES
        CHARACTER(LEN=1) :: BAR, BACK
        BACK = CHAR(8) !BACK SPACE
        BAR  = '='

        ITO = 6
        

        IF(I.EQ.1) THEN
            TIMES = SECOND() 
            CALL PRINBAR()
        ENDIF
        
        TIMER = 0.0D0
        IF(TIMES.NE.0.0D0) THEN
            TIMEN = SECOND() 
            TIMER = (TIMEN-TIMES)*(FLOAT(IMAX-I)/FLOAT(I))/60.0 !REMAINING TIME IN MINUTE
        ENDIF
        
        CALL DELTBAR()
        CALL PRINBAR()
C       CALL SLEEP(1) !PAUSE FOR A GIVE SECOND
        IF(I.EQ.IMAX) THEN
            WRITE(ITO,'(A)') '| complete.'
            TIMES = 0.0D0
        ENDIF

        
        CONTAINS
        
        SUBROUTINE PRINBAR
C     print the percentage and the bar
C        WRITE(ITO,'(2X,1I3,1A1,2X,1A1,256A1)',ADVANCE='NO')
C     1              100*I/IMAX,'%','|', (BAR, K=1,50*I/IMAX)
        WRITE(ITO,'(2X,1A10,1F7.2,1A3,2X,1I3,1A1,2X,1A1,256A1)',
     1              ADVANCE='NO')
     1 'TIME LEFT=',TIMER,'MIN',100*I/IMAX,'%','|', (BAR, K=1,38*I/IMAX)
        CLOSE(ITO)
        OPEN(ITO)
        END SUBROUTINE PRINBAR
        
        SUBROUTINE DELTBAR
C     delete the bar and the percentage
        WRITE(ITO,'(256A1)',ADVANCE='NO') (BACK, K=1,(38*I/IMAX)+32) !9 = THERE ARE 9 CHARACTER BEFORE REACH PROGRESS BAR (FROM LEFT)
        END SUBROUTINE DELTBAR
        
        END SUBROUTINE PROGBAR_O
C     =========================================================
C     =========================================================
C     =========================================================
