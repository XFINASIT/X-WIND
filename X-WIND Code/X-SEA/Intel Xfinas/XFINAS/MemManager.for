C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE INTMEMM
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(100000)
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1000000)

C	FOR INTEGER & CHARACTER
	MTOTM = 100000
	NTOTN = 1
	NUMA  = 0
      NEXT  = 1  !FIRST LOCATION OF ARRAY IN IMEMM
	MEMMI(1) = MTOTM
	MEMMI(2) = NTOTN
	MEMMI(3) = NUMA
	MEMMI(4) = NEXT

C	FOR REAL
	MTOTM = 1000000
	NTOTN = 1
	NUMA  = 0
      NEXT  = 1  !FIRST LOCATION OF ARRAY IN RMEMM
	MEMMR(1) = MTOTM
	MEMMR(2) = NTOTN
	MEMMR(3) = NUMA
	MEMMR(4) = NEXT

      RETURN
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE DEMREL(NAME,NA,NR,NC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----ALLOCATE A REAL ARRAY-----
      CHARACTER*1 NAME(4)
      NP = 2 !2 MEANS REAL
      CALL DEFIMM(NAME,NA,NR,NC,NP)
      RETURN
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE DEMINT(NAME,NA,NR,NC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----ALLOCATE AN INTEGER ARRAY-----
      CHARACTER*1 NAME(4)
      NP = 1 !1 MEANS INTEGER
      CALL DEFIMM(NAME,NA,NR,NC,NP)
      RETURN
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE DEMCHA(NAME,NA,NR,NC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----ALLOCATE A CHARACTER ARRAY-----
      CHARACTER*1 NAME(4)
      NP = 3 !3 MEANS CHARACTER
      CALL DEFIMM(NAME,NA,NR,NC,NP)
      RETURN
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE DEFIMM(NAME,NA,NR,NC,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----DEFINE AND RESERVE STORAGE FOR ARRAY -----
      CHARACTER*1  NAME(4)
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(1)
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)
C----------------------------------------------------
C WHERE NAME = NAME OF ARRAY - 4 LOGICALS MAXIMUM
C       NA = LOCATION OF ARRAY IF IN BLANK COMMON
C       NR = NUMBER OF ROWS
C       NC = NUMBER OF COLUMNS
C       MTOT = END OF DIRECTORY
C       NUMA = NUMBER OF ARRAYS IN DATA BASE
C       NEXT = NEXT AVAILABLE STORAGE LOCATION
C       IDIR = START OF DIRECTORY IN BLANK COMMON
C       IP   = NUMBER OF LOGICALS CONTAINED IN DATA TYPE
C       LENR = NUMBER OF LOGICALS IN PHYSICAL RECORD
C       NP = TYPE OF DATA
C          = 1 INTEGER DATA
C          = 2 REAL DATA
C          = 3 LOGICAL (CHARACTER) DATA
C------DIRECTORY DEFINITION FOR CORE OR SEQUENTIAL FILES
C      IDIR(1,N) = NAME OF ARRAY  - INAME (4 CHAR.)
C      IDIR(5,N) = NUMBER OF ROWS    - NR
C      IDIR(6,N) = NUMBER OF COLUMNS - NC
C      IDIR(7,N) = TYPE OF DATA      - NP
C      IDIR(8,N) = INCORE ADDRESS    - NA
C      IDIR(9,N) = SIZE OF ARRAY
C      IDIR(10,N) = 0 IF IN CORE STORAGE
C----------------------------------------------------
C-----EVALUATE STORAGE REQUIREMENTS -----
	DIMENSION IP(3)

C-----SET MACHINE PRECISION-----
C
      IP(1) = 1 !4
      IP(2) = 1 !8
      IP(3) = 1

      IF((NP.NE.1).AND.(NP.NE.3)) NP = 2
      NSIZE = (NR*NC*IP(NP) -1)/(IP(1)*2)
      NSIZE = NSIZE*2 + 2

C----------------------------------------------------

	SELECTCASE(NP)
	    CASE(1,3) !INT CHA
	    MTOTM = MEMMI(1) ; NTOTN = MEMMI(2) 
	    NUMA  = MEMMI(3) ; NEXT  = MEMMI(4)
	    CASE(2) !REL
	    MTOTM = MEMMR(1) ; NTOTN = MEMMR(2)
	    NUMA  = MEMMR(3) ; NEXT  = MEMMR(4)
	ENDSELECT

C----------------------------------------------------

      NA   = NEXT
      NEXT = NEXT + NSIZE
C-----SET UP NEW DIRECTORY -----
      NUMA  = NUMA + 1
      I     = NTOTN
      NTOTN = NTOTN + 10
C-----CHECK STORAGE LIMITS -----
CC      WRITE (*,*) 'TYPE:', NP
CC      WRITE (*,2000) NAME,NEXT,MTOTM
      
      IF(MTOTM.GE.NEXT) GO TO 100
      WRITE (*,2000) NAME,NEXT,MTOTM
      STOP
      
  100 CONTINUE


C----------------------------------------------------

	SELECTCASE(NP)
	    CASE(1,3) !INT CHA
	    CALL ICONI(NAME,NMEMI(I))
          NMEMI(I+4:I+9) = [NR,NC,NP,NA,NSIZE,0]
	    MEMMI(2:4)     = [NTOTN,NUMA,NEXT]
          
	    CASE(2) !REL
	    CALL ICONI(NAME,NMEMR(I))
          NMEMR(I+4:I+9) = [NR,NC,NP,NA,NSIZE,0]
	    MEMMR(2:4)     = [NTOTN,NUMA,NEXT]      
	ENDSELECT

C----------------------------------------------------


 900  RETURN
C
 2000 FORMAT(/,' ALLOCATION ERROR : ARRAY=',4A1,/,
     1' STORAGE REQUIRED  =',I7,/,' STORAGE AVAILABLE =',I7)
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE DEMREN(NAME1,NAME2,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME1,NAME2
      DIMENSION NAME1(4),NAME2(4),INAME(4)
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(1)
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)
C----------------------------------------------------
C	RENAME BLOCK DATA

C----------------------------------------------------

	SELECTCASE(NP)
	    CASE(1,3) !INT CHA
	        NUMA  = MEMMI(3)
              CALL ICONI(NAME1,INAME)
              I = IFINDL(INAME,0,NUMA,NMEMI)
              IF(I.EQ.0) GO TO 900
	        CALL ICONI(NAME2,NMEMI(I))
	    CASE(2) !REL
	        NUMA  = MEMMR(3)
              CALL ICONI(NAME1,INAME)
              I = IFINDL(INAME,0,NUMA,NMEMR)
              IF(I.EQ.0) GO TO 900
	        CALL ICONI(NAME2,NMEMR(I))
	ENDSELECT

C----------------------------------------------------


900	RETURN

      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE LOCATM(NAME,NA,NR,NC,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(1)
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NR = 0
	NC = 0

C	-----------------------------------------------

	SELECTCASE(NP)
	    CASE(1,3) !INT CHA
	        NUMA  = MEMMI(3)
              CALL ICONI(NAME,INAME)
              I = IFINDL(INAME,0,NUMA,NMEMI)
              IF(I.EQ.0) GO TO 900
              NA = NMEMI(I+7)
              NR = NMEMI(I+4)
              NC = NMEMI(I+5)
              NT = NMEMI(I+6)
	    CASE(2) !REL
	        NUMA  = MEMMR(3)
              CALL ICONI(NAME,INAME)
              I = IFINDL(INAME,0,NUMA,NMEMR)
              IF(I.EQ.0) GO TO 900
              NA = NMEMR(I+7)
              NR = NMEMR(I+4)
              NC = NMEMR(I+5)
              NT = NMEMR(I+6)
	ENDSELECT

C	-----------------------------------------------


  900 RETURN
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE DELMINT(NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(1)
      DIMENSION NAME(4),INAME(4)

C	-----------------------------------------------

	MTOTM = MEMMI(1)
	NTOTN = MEMMI(2)
	NUMA  = MEMMI(3)
	NEXT  = MEMMI(4)

C	-----------------------------------------------

C-----DELETE ARRAY FROM STORAGE -----
  100 CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NMEMI)
      IF(I.EQ.0) GO TO 900
C-----CHECK ON STORAGE LOCATION -----
  200 CONTINUE
	NSIZE = NMEMI(I+8)
C-----SET SIZE OF ARRAY -----
      NEXT = NEXT - NSIZE
      NUMA = NUMA - 1
      NA   = NMEMI(I+7)
C-----CHECK IF === IN-CORE OR DIRECT ACCESS -----
      IF(NA.GT.0) GO TO 500
      WRITE(*,1000) NAME
      GO TO 800
  500 IF(NA.EQ.NEXT) GO TO 800
C-----COMPACT STORAGE -----
      NNXT = NEXT - 1
      II   = NA + NSIZE
      JJ   = NNXT-NSIZE
      IMEMM(NA:JJ) = IMEMM(II:NNXT)

C-----COMPACT AND UPDATE DIRECTORY -----
  800 NTOTN = NTOTN - 10
	NA    = NTOTN - I
      IF(NA.EQ.0) GO TO 900
      NA = NA/10
      DO 860 K=1,NA
      II = I
      DO 850 J=1,10
      NMEMI(II) = NMEMI(II+10)
  850 II = II + 1
	
      NMEB = NMEMI(I+7)
      IF(NMEB.LE.0) GO TO 860
      NMEB = NMEMI(I+9)
      IF(NMEB.EQ.0) THEN
	NMEMI(I+7) = NMEMI(I+7) - NSIZE
	ENDIF

  860 I = I + 10
C
  900 CONTINUE

C	-----------------------------------------------

	MEMMI(2)   = NTOTN
	MEMMI(3)   = NUMA
	MEMMI(4)   = NEXT
	
C	-----------------------------------------------

	RETURN
      
 1000 FORMAT(/,' NAME ',4A1,' IS BEING USED FOR AN === IN-CORE FILE',/)
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE DELMREL(NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)
      DIMENSION NAME(4),INAME(4)

C	-----------------------------------------------

	MTOTM = MEMMR(1)
	NTOTN = MEMMR(2)
	NUMA  = MEMMR(3)
	NEXT  = MEMMR(4)

C	-----------------------------------------------


C-----DELETE ARRAY FROM STORAGE -----
  100 CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NMEMR)
      IF(I.EQ.0) GO TO 900
C-----CHECK ON STORAGE LOCATION -----
  200 CONTINUE
	NSIZE = NMEMR(I+8)
C-----SET SIZE OF ARRAY -----
      NEXT = NEXT - NSIZE
      NUMA = NUMA - 1
      NA   = NMEMR(I+7)
C-----CHECK IF === IN-CORE OR DIRECT ACCESS -----
      IF(NA.GT.0) GO TO 500
      WRITE(*,1000) NAME
      GO TO 800
  500 IF(NA.EQ.NEXT) GO TO 800
C-----COMPACT STORAGE -----
      NNXT = NEXT - 1
      II   = NA + NSIZE
      JJ   = NNXT-NSIZE
      RMEMM(NA:JJ) = RMEMM(II:NNXT)

C-----COMPACT AND UPDATE DIRECTORY -----
  800 NTOTN = NTOTN - 10
	NA    = NTOTN - I
      IF(NA.EQ.0) GO TO 900
      NA = NA/10
      DO 860 K=1,NA
      II = I
      DO 850 J=1,10
      NMEMR(II) = NMEMR(II+10)
  850 II = II + 1
  
      NMEB = NMEMR(I+7)
      IF(NMEB.LE.0) GO TO 860
      NMEB = NMEMR(I+9)
      IF(NMEB.EQ.0) THEN
	NMEMR(I+7) = NMEMR(I+7) - NSIZE
	ENDIF
  860 I = I + 10
C
  900 CONTINUE

C	-----------------------------------------------

	MEMMR(2)   = NTOTN
	MEMMR(3)   = NUMA
	MEMMR(4)   = NEXT

C	-----------------------------------------------

	RETURN

 1000 FORMAT(/,' NAME ',4A1,' IS BEING USED FOR AN === IN-CORE FILE',/)
      END

C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
      SUBROUTINE MRELZER(NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)
      DIMENSION NAME(4),INAME(4)

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MEMMR(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NMEMR)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NMEMR(I+7)
      NR = NMEMR(I+4)
      NC = NMEMR(I+5)
      NP = NMEMR(I+6)

	KK  = NA
      NSIZE = NR*NC
      RMEMM(KK:KK+NSIZE-1) = 0.0D0

      
  900 RETURN
      END
C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================

      SUBROUTINE MINTZER(NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(1)
      DIMENSION NAME(4),INAME(4)
      
C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MEMMI(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NMEMI)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NMEMI(I+7)
      NR = NMEMI(I+4)
      NC = NMEMI(I+5)
      NP = NMEMI(I+6)

	KK  = NA
      NSIZE = NR*NC
      IMEMM(KK:KK+NSIZE-1) = 0

      
  900 RETURN
      END
C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================

      SUBROUTINE MRELFIL (NAME,VALV,MR,MC,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)
      DIMENSION NAME(4),INAME(4)

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MEMMR(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NMEMR)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NMEMR(I+7)
      NR = NMEMR(I+4)
      NC = NMEMR(I+5)
      NP = NMEMR(I+6)

	NN = MR + NR*(MC-1) - 1
	IF(NN.LT.0) GO TO 900
	NA = NA + NN
      
	CALL RFILL(RMEMM(NA),VALV,IND)

      MEMMR(10) = MEMMR(10) + 1 !COUNTING HOW MANY TIME THIS SUB. HAS BEEN CALLED ... REAL OPERATION

  900 RETURN
      END


C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================

      SUBROUTINE MRELFILA (NAME,VALV,MR,MCi,MCj,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)
      
      ALLOCATABLE MAPP(:)
      
      DIMENSION NAME(4),INAME(4)
      DIMENSION VALV(1)

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MEMMR(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NMEMR)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NMEMR(I+7)
      NR = NMEMR(I+4)
      NC = NMEMR(I+5)
      NP = NMEMR(I+6)

	NNi = MR + NR*(MCi-1) - 1
	NNj = MR + NR*(MCj-1) - 1
	IF(NNi.LT.0.OR.NNj.LT.0) GO TO 900
	NAi = NA + NNi
	NAj = NA + NNj
      
      NLEN = (MCj-MCi)+1
      ALLOCATE(MAPP(NLEN))
      MAPP = 0
      MAPP(1:NLEN) = MR + NR*([MCi:MCj]-1)
      
	CALL RFILLA(RMEMM(NA),VALV(1),MAPP,NLEN,IND)

      DEALLOCATE(MAPP)
      
      MEMMR(10) = MEMMR(10) + 1 !COUNTING HOW MANY TIME THIS SUB HAS BEEN CALLED ... REAL OPERATION
      
      
  900 RETURN
      END


C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================

      SUBROUTINE MINTFIL (NAME,IVAL,MR,MC,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(1)
      DIMENSION NAME(4),INAME(4)

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MEMMI(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NMEMI)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NMEMI(I+7)
      NR = NMEMI(I+4)
      NC = NMEMI(I+5)
      NP = NMEMI(I+6)

	NN = MR + NR*(MC-1) - 1
	IF(NN.LT.0) GO TO 900
	NA = NA + NN
      
	CALL IFILL(IMEMM(NA),IVAL,IND)

      MEMMI(10) = MEMMI(10) + 1 !COUNTING HOW MANY TIME THIS SUB HAS BEEN CALLED ... INTEGER OPERATION

      
  900 RETURN
      END


C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================

      SUBROUTINE ICONC(INAME,NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME(4)
      DIMENSION INAME(4)
C-----CONVERT INTEGER TO CHARACTER DATA -----
      DO I = 1,4
50        IF(INAME(I).GT.255) THEN
            INAME(I-1) = INAME(I-1) + 20
            INAME(I  ) = INAME(I  ) - 255
          ENDIF
          IF(INAME(I).GT.255) GOTO 50
      ENDDO
      
      DO 100 I = 1,4
  100 NAME(I) = CHAR( INAME(I) )
      RETURN
      END
C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================

      SUBROUTINE ICON4(NAME,N1,N2,N3,N4)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME(4)
      DIMENSION INAME(4)
	INAME(1:4) = [N1,N2,N3,N4]
C-----CONVERT INTEGER TO CHARACTER DATA -----
      DO 100 I = 1,4
  100 NAME(I) = CHAR( INAME(I) )
      RETURN
      END
C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================

      SUBROUTINE MWCOUNT
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      COMMON /MEMFILE/ LFILI,LFILR
      COMMON /MEMTINT/ MEMMI(10),NMEMI(10000),IMEMM(1)
      COMMON /MEMTREL/ MEMMR(10),NMEMR(10000),RMEMM(1)
      DIMENSION NAME(4),INAME(4)


      WRITE(*,*) 'INTG OPER. IN MEM.MANAGER ',MEMMI(10)

      WRITE(*,*) 'REAL OPER. IN MEM.MANAGER ',MEMMR(10)
      
  900 RETURN
      END


C	=======================================================================
C	========== IN-CORE ARRAY UTILITY === SONGSAK NOV2007 / OCT2019 ========
C	=======================================================================
