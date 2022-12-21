C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE INTCONS
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(2000000)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(6000000)

C	FOR INTEGER & CHARACTER
	MTOTM = 2000000
	NTOTN = 1
	NUMA  = 0
      NEXT  = 1  !FIRST LOCATION OF ARRAY IN ICONDT
	MTOTI(1) = MTOTM
	MTOTI(2) = NTOTN
	MTOTI(3) = NUMA
	MTOTI(4) = NEXT

C	FOR REAL
	MTOTM = 6000000
	NTOTN = 1
	NUMA  = 0
      NEXT  = 1  !FIRST LOCATION OF ARRAY IN RCONDT
	MTOTR(1) = MTOTM
	MTOTR(2) = NTOTN
	MTOTR(3) = NUMA
	MTOTR(4) = NEXT



      RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FREBUF
C-----PUT NEXT DATA LINE IN BUFFER - MULTILINE INPUT VERSION-----
      CHARACTER*1 LINE,SP,BS,SC,C,MLIN,ML,MLINS

      COMMON /MLINE/ MLIN(161),MLINS(161)
      COMMON /CLINE/ LINE(160)
      COMMON /ILINE/ II

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DATA ML/'!'/
      DATA SP/' '/,SC/';'/,C/'C'/,BS/'\'/

	NIN = ITI
	NTM = ITO

C-----CLEAR LINE BUFFER-----
   10 DO 15 K=1,160
   15 LINE(K) = SP
C-----SEE IF MULTI-LINE BUFFER IS EMPTY-----
      IF(MLIN(1).EQ.ML) GO TO 50
C-----READ LINE OF FREE FIELD DATA ----
   20 J = 2
      JJ= 81
   30 READ (NIN,1001,ERR=100,END=100) (MLIN(K),K=J,JJ)
C-----CHECK FOR ADDITIONAL LINE -----
   40 DO 45 K=J,JJ
      IF(MLIN(K).NE.BS) GO TO 45
      J = K
      JJ= K+79
      IF(JJ.GT.161) JJ = 161
      GO TO 30
   45 CONTINUE
C-----CHECK FOR COMMENT LINE (;; MEANS NEGLECT PHYSICAL LINE)-----
      IF(MLIN(2).EQ.SC.AND.MLIN(3).EQ.SC) GO TO 20
C-----TRANSFER LINE TO SINGLE-LINE BUFFER-----
   50 DO 60 K=1,159
      IF(MLIN(K+1).NE.ML) GO TO 60
      II = K
      GO TO 65
   60 LINE(K) = MLIN(K+1)
      II = 159
      KK = 1
      GO TO 75
C-----SHIFT MULTI-LINE BUFFER-----
   65 JJ = II + 1
      KK = 1
      DO 70 K=JJ,161
      MLIN(KK) = MLIN(K)
   70 KK = KK + 1
   75 DO 80 K=KK,161
   80 MLIN(K) = SP
C-----CHECK FOR COMMENT LINE (; MEANS NEGLECT LOGICAL LINE)-----
      CALL FREUPP
      IF(LINE(1).NE.SC) GO TO 900
      WRITE (NTM,2000) (LINE(I),I=1,II)
      GO TO 10
C-----ERROR IN READ -----
  100 WRITE(NTM,2001)
C
  900 RETURN
C
 1000 FORMAT (80A1)
 1001 FORMAT (80A1)
 2000 FORMAT (1X,80A1)
 2001 FORMAT(/,' * END OF FILE OR ERROR IN READING INPUT LINE *',/)
 3000 FORMAT (' :>',\) 

      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FREUPP
      CHARACTER*1 LINE,SP,SC

      COMMON /CLINE/ LINE(160)
      COMMON /ILINE/ II
      DATA SP/' '/,SC/';'/
C-----CHECK END OF LINE AND CONVERT TO UPPER CASE ---
      ISP = ICHAR(SP)
      DO 80 I=1,II
      IF(LINE(I).NE.SC) GO TO 75
      IF(I.EQ.1) GO TO 75
      II = I
      GO TO 85
   75 NN = ICHAR( LINE(I) )
      IF(NN.GT.96) LINE(I)= CHAR( NN - ISP )

   80 CONTINUE
C-----REMOVE BLANKS AT END OF LINE -----
   85 IF(LINE(II).NE.SP) GO TO 95
      IF(II.EQ.1) GO TO 95
      II=II-1
      GO TO 85
C
   95 RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FREECH
C-----ECHO OF FREE FIELD INFORMATION TO TERMINAL----
      CHARACTER*1 LINE

      COMMON /ILINE/ II
      COMMON /CLINE/ LINE(160)

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

	NTM = ITO

	WRITE(NTM,1000) (LINE(I),I=1,II)

      RETURN
 1000 FORMAT (1X,78A1)
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE FREINT(IC,IDATA,NUM)
      CHARACTER*1 LINE,IC,BLK,COMMA,I0,EQS,NEG,COLIN,LNE
      DIMENSION IDATA(72)
      COMMON /CLINE/ LINE(160)
      COMMON /ILINE/ II

      DATA BLK/' '/,COMMA/','/,I0/'0'/,EQS/'='/,NEG/'-'/,COLIN/':'/
C-----FIND INTEGER STRING ----
   90 I=0
      IF(IC.EQ.BLK) GO TO 200
      DO 100 I=1,II
      IF((LINE(I).EQ.IC).AND.(LINE(I+1).EQ.EQS)) GO TO 200
  100 CONTINUE
      RETURN
C-----ZERO INTEGER STRING -----
  200 DO 210 J=1,NUM
  210 IDATA(J)=0
      IF(LINE(I+1).EQ.EQS) I=I+1
      DO 250 J=1,NUM
      ISIGN = 1
C-----SKIP BLANKS BETWEEN INTEGERS -----
  215 IF(LINE(I+1).NE.BLK) GO TO 220
      I=I+1
      IF(I.GT.II) GO TO 900
      GO TO 215
  220 I=I+1
      IF(I.GT.II) GO TO 230
C-----CHECK FOR SIGN -----
      LNE = LINE(I)
      IF(LNE.NE.NEG) GO TO 225
      ISIGN = -1
      GO TO 220
C-----EXTRACT INTEGER -----
  225 IF(LNE.EQ.BLK) GO TO 230
      IF(LNE.EQ.COMMA) GO TO 230
      IF(LNE.EQ.COLIN) GO TO 230
      NN = ICHAR( LNE ) - ICHAR( I0 )
      IF((NN.LT.0).OR.(NN.GT.9)) GO TO 900
      IDATA(J)=10*IDATA(J)+NN
      GO TO 220
C-----SET SIGN -----
  230 IDATA(J) = IDATA(J)*ISIGN
  250 CONTINUE
  900 RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FREREL(IC,DATA,NUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 LINE,BLK,IC,MUL,DIV,ADD,SUB,EQS,E
      DIMENSION DATA(30)                                                   
      COMMON /CLINE/ LINE(160)
      COMMON /ILINE/ II
      DATA BLK/' '/,MUL/'*'/,DIV/'/'/,ADD/'+'/,SUB/'-'/,E/'E'/,EQS/'='/
C-----FIND REAL STRING -----
   90 I=0
      IF(IC.EQ.BLK) GO TO 250
      DO 100 I=1,II
      IF((LINE(I).EQ.IC).AND.(LINE(I+1).EQ.EQS)) GO TO 250
  100 CONTINUE
      RETURN
C---- EXTRACT REAL DATA -----
  250 DO 260 J=1,NUM
  260 DATA(J)=0.0D+0
      DO 300 J=1,NUM
      JJ=0
  270 IF(I.GT.II) GO TO 300
      CALL FRERDT(I,XX,NN)
      IF(JJ.NE.0) GO TO 275
      DATA(J) = XX
      GO TO 290
C-----ARITHMETRIC STATEMENT -----
  275 IF(JJ.EQ.1) DATA(J)=DATA(J)*XX
      IF(JJ.EQ.2) DATA(J)=DATA(J)/XX
      IF(JJ.EQ.3) DATA(J)=DATA(J)+XX
      IF(JJ.EQ.4) DATA(J)=DATA(J)-XX
      IF(JJ.NE.5) GO TO 290
C-----EXPONENTIAL DATA -----
      JJ = DABS(XX)
      IF(JJ.EQ.0) GO TO 290
      DO 280 K=1,JJ
      IF(XX.LT.0.0) DATA(J) = DATA(J)/10.
      IF(XX.GT.0.0) DATA(J) = DATA(J)*10.
  280 CONTINUE
C-----SET TYPE OF STATEMENT -----
  290 JJ=0
      IF(LINE(I).EQ.MUL) JJ=1
      IF(LINE(I).EQ.DIV) JJ=2
      IF(LINE(I).EQ.ADD) JJ=3
      IF(LINE(I).EQ.SUB) JJ=4
      IF(LINE(I).EQ.E)   JJ=5
      IF(LINE(I+1).EQ.EQS) JJ=0
      IF(JJ.NE.0) GO TO 270
      IF(NN.GT.9) RETURN
  300 CONTINUE
C
      RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE  FRERDT(I,XX,NN)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 LINE,BLK,COMMA,I0,DOT,NEG,EQS,ADD,E
      COMMON /CLINE/ LINE(160)
      COMMON /ILINE/ II
      DATA BLK/' '/,COMMA/','/,I0/'0'/,DOT/'.'/,NEG/'-'/,EQS/'='/
	DATA ADD/'+'/,E/'E'/ !ADDED DEC2005
C-----CONVERT STRING TO REAL FLOATING POINT NUMBER --
      IF(LINE(I+1).EQ.EQS) I=I+1
	IF((LINE(I).EQ.E).AND.(LINE(I+1).EQ.ADD)) I=I+1
      Y=0
      IS=1
      XX=0.0
      IF(LINE(I+1).NE.NEG) GO TO 267
      IS=-1
      I=I+1
  267 IF(LINE(I+1).NE.BLK) GO TO 270
      I=I+1
      IF(I.GT.II) GO TO 300
      GO TO 267
  270 I=I+1
      IF(I.GT.II) GO TO 300
      IF((LINE(I).EQ.BLK).AND.(LINE(I+1).EQ.BLK)) GO TO 270
      NN = ICHAR( LINE(I) ) - ICHAR( I0 )
      XN=ISIGN(NN,IS)
      IF(LINE(I).NE.DOT) GO TO 275
      Y=1.0D+0
      GO TO 270
  275 IF(LINE(I).EQ.BLK) GO TO 300
      IF(LINE(I).EQ.COMMA) GO TO 300
      IF((NN.LT.0).OR.(NN.GT.9)) GO TO 300
      IF(Y.EQ.0.0D+0) GO TO 280
      Y=Y/10.D+0
      XN=XN*Y
      XX=XX+XN
      GO TO 270
  280 XX=10.0D+0*XX+XN
      GO TO 270
  300 RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE FRECHA(IC,IDATA,NC,NUM)
      CHARACTER*1 LINE,IC,BLK,COMMA,I0,EQS,COLIN,IDATA
      DIMENSION IDATA(NC,NUM)
      COMMON /CLINE/ LINE(160)
      COMMON /ILINE/ II
      DATA BLK/' '/,COMMA/','/,I0/'0'/,EQS/'='/,COLIN/':'/
C-----FIND HOLLERITH STRING -----
   90 I=0
      IF(IC.EQ.BLK) GO TO 200
      DO 100 I=2,II
      IF((LINE(I-1).EQ.IC).AND.(LINE(I).EQ.EQS)) GO TO 200
  100 CONTINUE
      RETURN
C-----EXTRACT HOLLERITH DATA -----
  200 DO 210 J=1,NUM
      DO 210 N=1,NC
  210 IDATA(N,J)=BLK
C
      DO 300 J=1,NUM
  260 I = I + 1
      IF(I.GT.II) GO TO 400
      IF(LINE(I).EQ.COMMA) GO TO 260
      IF(LINE(I).EQ.BLK) GO TO 260
      DO 290 N=1,NC
      IF(LINE(I).EQ.COLIN) GO TO 300
      IF(LINE(I).EQ.BLK)   GO TO 300
      IF(LINE(I).EQ.COMMA) GO TO 300
      IDATA(N,J) = LINE(I)
      IF(N.EQ.NC) GO TO 290
      I = I + 1
  290 CONTINUE
  300 CONTINUE
C
  400 CONTINUE

	RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      FUNCTION FRETES (OP,NOP,NC)
C     RETURN TRUE IF STRINGS MATCH    -
      LOGICAL FRETES
      CHARACTER*1 OP(NC),NOP(NC)
      FRETES = .FALSE.
      DO 10 I=1,NC
      IF(NOP(I).NE.OP(I)) RETURN
   10 CONTINUE
      FRETES = .TRUE.
      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE DEFNREL(NAME,NA,NR,NC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----ALLOCATE A REAL ARRAY-----
      CHARACTER*1 NAME(4)
      NP = 2 !2 MEANS REAL
      CALL DEFINN(NAME,NA,NR,NC,NP)
      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE DEFNINT(NAME,NA,NR,NC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----ALLOCATE AN INTEGER ARRAY-----
      CHARACTER*1 NAME(4)
      NP = 1 !1 MEANS INTEGER
      CALL DEFINN(NAME,NA,NR,NC,NP)
      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE DEFNCHA(NAME,NA,NR,NC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----ALLOCATE A CHARACTER ARRAY-----
      CHARACTER*1 NAME(4)
      NP = 3 !3 MEANS CHARACTER
      CALL DEFINN(NAME,NA,NR,NC,NP)
      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE DEFINN(NAME,NA,NR,NC,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----DEFINE AND RESERVE STORAGE FOR ARRAY -----
      CHARACTER*1  NAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)
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
	DIMENSION IP(3),INAME(4)
      
      ALLOCATABLE IA(:),AA(:)
C
C-----SET MACHINE PRECISION-----
C
      IP(1) = 1 !4
      IP(2) = 1 !8
      IP(3) = 1

      IF((NP.NE.1).AND.(NP.NE.3)) NP = 2
      NSIZE = (NR*NC*IP(NP) -1)/(IP(1)*2)
      NSIZE = NSIZE*2 + 2
      
	SELECTCASE(NP)
	CASE(1,3) !INT CHA
	MTOTM = MTOTI(1)
	NTOTN = MTOTI(2)
	NUMA  = MTOTI(3)
	NEXT  = MTOTI(4)
	
      CALL ICONI(NAME,INAME)
      ITEST = IFINDL(INAME,0,NUMA,NCONI)
      IF(ITEST.NE.0) GOTO 900 !BLOCK ALREADY EXIST
      
	CASE(2) !REL
	MTOTM = MTOTR(1)
	NTOTN = MTOTR(2)
	NUMA  = MTOTR(3)
	NEXT  = MTOTR(4)
	
      CALL ICONI(NAME,INAME)
      ITEST = IFINDL(INAME,0,NUMA,NCONR)
      IF(ITEST.NE.0) GOTO 900 !BLOCK ALREADY EXIST
      ENDSELECT

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

	SELECTCASE(NP)
	CASE(1,3) !INT CHA
	CALL ICONI(NAME,NCONI(I))
      NCONI(I+4) = NR
      NCONI(I+5) = NC
      NCONI(I+6) = NP
      NCONI(I+7) = NA
      NCONI(I+8) = NSIZE
      NCONI(I+9) = 0
	MTOTI(2)   = NTOTN
	MTOTI(3)   = NUMA
	MTOTI(4)   = NEXT
      
	CASE(2) !REL
	CALL ICONI(NAME,NCONR(I))
      NCONR(I+4) = NR
      NCONR(I+5) = NC
      NCONR(I+6) = NP
      NCONR(I+7) = NA
      NCONR(I+8) = NSIZE
      NCONR(I+9) = 0
	MTOTR(2)   = NTOTN
	MTOTR(3)   = NUMA
	MTOTR(4)   = NEXT

	ENDSELECT

 900  RETURN
C
 2000 FORMAT(/,' ALLOCATION ERROR: ARRAY=',4A1,/,' STORAGE REQUIRED  =',
     1 I17,/,' STORAGE AVAILABLE =',I17)
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE DEFREN(NAME1,NAME2,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----DEFINE AND RESERVE STORAGE FOR ARRAY -----
      CHARACTER*1 NAME1,NAME2
      DIMENSION NAME1(4),NAME2(4),INAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)
C----------------------------------------------------
C	RENAME BLOCK DATA

	SELECTCASE(NP)

	CASE(1,3) !INT CHA
	NUMA  = MTOTI(3)
	ITEST = 0

      CALL ICONI(NAME1,INAME)
      I1 = IFINDL(INAME,0,NUMA,NCONI)
      IF(I1.EQ.0) ITEST = 1

      CALL ICONI(NAME2,INAME)
      I2 = IFINDL(INAME,0,NUMA,NCONI)
      IF(I2.EQ.0) ITEST = 1

	IF(ITEST.EQ.1) GO TO 900

	CALL ICONI(NAME2,NCONI(I1))
	CALL ICONI(NAME1,NCONI(I2))


	CASE(2) !REL
	NUMA  = MTOTR(3)
	ITEST = 0

      CALL ICONI(NAME1,INAME)
      I1 = IFINDL(INAME,0,NUMA,NCONR)
      IF(I1.EQ.0) ITEST = 1

      CALL ICONI(NAME2,INAME)
      I2 = IFINDL(INAME,0,NUMA,NCONR)
      IF(I2.EQ.0) ITEST = 1

	IF(ITEST.EQ.1) GO TO 900

	CALL ICONI(NAME2,NCONR(I1))
	CALL ICONI(NAME1,NCONR(I2))

	ENDSELECT



900	RETURN

      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE LOCATN(NAME,NA,NR,NC,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)


C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NR = 0
	NC = 0

	SELECTCASE(NP)

	CASE(1,3) !INT CHA
	NUMA  = MTOTI(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
      IF(I.EQ.0) GO TO 900
      NA = NCONI(I+7)
      NR = NCONI(I+4)
      NC = NCONI(I+5)
      NT = NCONI(I+6)

	CASE(2) !REL
	NUMA  = MTOTR(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
      IF(I.EQ.0) GO TO 900
      NA = NCONR(I+7)
      NR = NCONR(I+4)
      NC = NCONR(I+5)
      NT = NCONR(I+6)

	ENDSELECT



  900 RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE LOCADD(NAME1,NAME2,FAC,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME1,NAME2
      DIMENSION NAME1(4),NAME2(4),INAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)


C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NR = 0
	NC = 0

	SELECTCASE(NP)

	CASE(1,3) !INT CHA
	NUMA  = MTOTI(3)
      CALL ICONI(NAME1,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
      IF(I.EQ.0) GO TO 900
      NA1 = NCONI(I+7)
      NR1 = NCONI(I+4)
      NC1 = NCONI(I+5)
      NT  = NCONI(I+6)
      CALL ICONI(NAME2,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
      IF(I.EQ.0) GO TO 900
      NA2 = NCONI(I+7)
      NR2 = NCONI(I+4)
      NC2 = NCONI(I+5)
      NT  = NCONI(I+6)
	IF(NR1.NE.NR2.OR.NC1.NE.NC2) GOTO 900
	IFAC = INT(FAC)
	DO J = 1,NR2*NC2
	J1 = NA1 + J - 1
	J2 = NA2 + J - 1
	ICONDT(J1) = ICONDT(J1) + ICONDT(J2)*IFAC
	ENDDO

	CASE(2) !REL
	NUMA  = MTOTR(3)
      CALL ICONI(NAME1,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
      IF(I.EQ.0) GO TO 900
      NA1 = NCONR(I+7)
      NR1 = NCONR(I+4)
      NC1 = NCONR(I+5)
      NT  = NCONR(I+6)
      CALL ICONI(NAME2,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
      IF(I.EQ.0) GO TO 900
      NA2 = NCONR(I+7)
      NR2 = NCONR(I+4)
      NC2 = NCONR(I+5)
      NT  = NCONR(I+6)
	IF(NR1.NE.NR2.OR.NC1.NE.NC2) GOTO 900
	DO J = 1,NR2*NC2
	J1 = NA1 + J - 1
	J2 = NA2 + J - 1
	RCONDT(J1) = RCONDT(J1) + RCONDT(J2)*FAC
	ENDDO


	ENDSELECT



  900 RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SERHINT(NAME,I)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)

	MTOTM = MTOTI(1)
	NTOTN = MTOTI(2)
	NUMA  = MTOTI(3)
	NEXT  = MTOTI(4)
	
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
C     IF I = 0 MEANS NOT FOUND

      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE SERHREL(NAME,I)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

	MTOTM = MTOTR(1)
	NTOTN = MTOTR(2)
	NUMA  = MTOTR(3)
	NEXT  = MTOTR(4)
	
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
C     IF I = 0 MEANS NOT FOUND

      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE DELTINT(NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)

	MTOTM = MTOTI(1)
	NTOTN = MTOTI(2)
	NUMA  = MTOTI(3)
	NEXT  = MTOTI(4)

C-----DELETE ARRAY FROM STORAGE -----
  100 CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
      IF(I.EQ.0) GO TO 900
C-----CHECK ON STORAGE LOCATION -----
  200 NSIZE = NCONI(I+8)
C-----SET SIZE OF ARRAY -----
      NEXT = NEXT - NSIZE
      NUMA = NUMA - 1
      NA   = NCONI(I+7)
C-----CHECK IF OUT OF CORE OR DIRECT ACCESS -----
      IF(NA.GT.0) GO TO 500
      WRITE(*,1000) NAME
      GO TO 800
  500 IF(NA.EQ.NEXT) GO TO 800
C-----COMPACT STORAGE -----
      NNXT = NEXT - 1
      II   = NA + NSIZE
      JJ   = NNXT-NSIZE
      ICONDT(NA:JJ) = ICONDT(II:NNXT)
CC      DO 700 J=NA,NNXT
CC      ICONDT(J) = ICONDT(II)
CC  700 II = II + 1
      
C-----COMPACT AND UPDATE DIRECTORY -----
  800 NTOTN = NTOTN - 10
	NA    = NTOTN - I
      IF(NA.EQ.0) GO TO 900
      NA = NA/10
      DO 860 K=1,NA
      II = I
      DO 850 J=1,10
      NCONI(II) = NCONI(II+10)
  850 II = II + 1
      IF(NCONI(I+7).LE.0) GO TO 860
      IF(NCONI(I+9).EQ.0) NCONI(I+7) = NCONI(I+7) - NSIZE
  860 I = I + 10
C
  900 CONTINUE

	MTOTI(2)   = NTOTN
	MTOTI(3)   = NUMA
	MTOTI(4)   = NEXT
	
	RETURN
 1000 FORMAT(/,' NAME ',4A1,' IS BEING USED FOR AN OUT OF CORE FILE',/)
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE DELTREL(NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

	MTOTM = MTOTR(1)
	NTOTN = MTOTR(2)
	NUMA  = MTOTR(3)
	NEXT  = MTOTR(4)

C-----DELETE ARRAY FROM STORAGE -----
  100 CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
      IF(I.EQ.0) GO TO 900
C-----CHECK ON STORAGE LOCATION -----
  200 NSIZE = NCONR(I+8)
C-----SET SIZE OF ARRAY -----
      NEXT = NEXT - NSIZE
      NUMA = NUMA - 1
      NA   = NCONR(I+7)
C-----CHECK IF OUT OF CORE OR DIRECT ACCESS -----
      IF(NA.GT.0) GO TO 500
      WRITE(*,1000) NAME
      GO TO 800
  500 IF(NA.EQ.NEXT) GO TO 800
C-----COMPACT STORAGE -----
      NNXT = NEXT - 1
      II   = NA + NSIZE
      JJ   = NNXT-NSIZE
      RCONDT(NA:JJ) = RCONDT(II:NNXT)
CC      DO 700 J=NA,NNXT
CC      RCONDT(J) = RCONDT(II)
CC  700 II = II + 1
      
C-----COMPACT AND UPDATE DIRECTORY -----
  800 NTOTN = NTOTN - 10
	NA    = NTOTN - I
      IF(NA.EQ.0) GO TO 900
      NA = NA/10
      DO 860 K=1,NA
      II = I
      DO 850 J=1,10
      NCONR(II) = NCONR(II+10)
  850 II = II + 1
      IF(NCONR(I+7).LE.0) GO TO 860
      IF(NCONR(I+9).EQ.0) NCONR(I+7) = NCONR(I+7) - NSIZE
  860 I = I + 10
C
  900 CONTINUE

	MTOTR(2)   = NTOTN
	MTOTR(3)   = NUMA
	MTOTR(4)   = NEXT
	
	RETURN

 1000 FORMAT(/,' NAME ',4A1,' IS BEING USED FOR AN OUT OF CORE FILE',/)
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE ICONI(NAME,INAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME(4)
      DIMENSION INAME(4)
C-----CONVERT LOGICALS TO INTEGER DATA -----
      DO 100 I = 1,4
  100 INAME(I) = ICHAR( NAME(I) )
      RETURN
      END
C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      FUNCTION IFINDL(INAME,LUN,NUMA,IA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION INAME(4),IA(1)
C-----FIND ARRAY LOCATION -----
      I = 1
      DO 100 N=1,NUMA
      IF(LUN.NE.IA(I+9)) GO TO 100
      IF (INAME(1).NE.IA(I  )) GO TO 100
      IF (INAME(2).NE.IA(I+1)) GO TO 100
      IF (INAME(3).NE.IA(I+2)) GO TO 100
      IF (INAME(4).EQ.IA(I+3)) GO TO 200
  100 I = I + 10
      I = 0
  200 IFINDL = I
      RETURN
      END
C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE RELZERO (NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MTOTR(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NCONR(I+7)
      NR = NCONR(I+4)
      NC = NCONR(I+5)
      NP = NCONR(I+6)
	KK  = NA
      
      NSIZE = NR*NC
      RCONDT(KK:KK+NSIZE-1) = 0.0D0

C	DO IR = 1,NR
C	DO IC = 1,NC
C	RCONDT(KK)  = 0.0D0              
C	KK = KK + 1
C	ENDDO
C	ENDDO
      
  900 RETURN
      END
C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE INTZERO (NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MTOTI(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NCONI(I+7)
      NR = NCONI(I+4)
      NC = NCONI(I+5)
      NP = NCONI(I+6)
	KK  = NA
      
      NSIZE = NR*NC
      ICONDT(KK:KK+NSIZE-1) = 0
      
C	DO IR = 1,NR
C	DO IC = 1,NC
C	ICONDT(KK)  = 0
C	KK = KK + 1
C	ENDDO
C	ENDDO
      
  900 RETURN
      END
C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE RELFILL (NAME,VALV,MR,MC,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)
C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MTOTR(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NCONR(I+7)
      NR = NCONR(I+4)
      NC = NCONR(I+5)
      NP = NCONR(I+6)
	NN = MR + NR*(MC-1) - 1
	IF(NN.LT.0) THEN
	IF(IND.EQ.0) VALV = 0.0D0
	GO TO 900
	ENDIF
	NA = NA + NN
	CALL RFILL(RCONDT(NA),VALV,IND)
  900 RETURN
      END

      
C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE INTFILL (NAME,IVAL,MR,MC,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)
C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MTOTI(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NCONI(I+7)
      NR = NCONI(I+4)
      NC = NCONI(I+5)
      NP = NCONI(I+6)
	NN = MR + NR*(MC-1) - 1
	IF(NN.LT.0) THEN
	IF(IND.EQ.0) IVAL = 0
	GO TO 900
	ENDIF
	NA = NA + NN
	CALL IFILL(ICONDT(NA),IVAL,IND)
  900 RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE RFILL (FA,VALV,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	SELECTCASE(IND)
	
	CASE(0)
      VALV = FA

	CASE(1)
      FA = VALV

	CASE(2)
      FA = FA + VALV

	CASE(3)
      FA = FA - VALV

	CASE(4)
      FA = FA * VALV

	CASE(5)
      FA = FA / VALV

	CASE(6)
      FA = VALV/FA

	CASE(7)
      FA = ABS(VALV)

	ENDSELECT


      RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE RFILLA (FA,VALV,MAPP,N,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION FA(N),VALV(N)
      DIMENSION MAPP(N)
      
	SELECTCASE(IND)
	
	CASE(0)
      VALV(1:N) = FA(MAPP(1:N))

	CASE(1)
      FA(MAPP(1:N)) = VALV(1:N)

	CASE(2)
      FA(MAPP(1:N)) = FA(MAPP(1:N)) + VALV(1:N)

	CASE(3)
      FA(MAPP(1:N)) = FA(MAPP(1:N)) - VALV(1:N)

	CASE(4)
      FA(MAPP(1:N)) = FA(MAPP(1:N)) * VALV(1:N)

	CASE(5)
      FA(MAPP(1:N)) = FA(MAPP(1:N)) / VALV(1:N)

	CASE(6)
      FA(MAPP(1:N)) = VALV(1:N)/FA(MAPP(1:N))

	CASE(7)
      FA(MAPP(1:N)) = ABS(VALV(1:N))

	ENDSELECT


      RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE IFILL (IA,IVAL,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


	SELECTCASE(IND)
	
	CASE(0)
      IVAL = IA

	CASE(1)
      IA = IVAL

	CASE(2)
      IA = IA + IVAL

	CASE(3)
      IA = IA - IVAL

	CASE(4)
      IA = IA * IVAL

	ENDSELECT


      RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

	FUNCTION UNITFAC(NAME,IO,IN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=============================================================================
      CHARACTER*4 NAME
C	=============================================================================

	UNITFAC = 1.0
	IF(NAME.EQ.'FORC') INF = 0 !FORC
	IF(NAME.EQ.'LENG') INF = 1 !LENG
	IF(NAME.EQ.'ROTA') INF = 2 !ROTA
	IF(NAME.EQ.'MASS') INF = 3 !MASS
	IF(NAME.EQ.'TIME') INF = 4 !TIME
	IF(NAME.EQ.'TEMP') INF = 5 !TEMP
	GOTO 10
	RETURN
10	CONTINUE

	DO I = 1,2
	
	INU = IO
	IF(I.EQ.2) INU = IN

	SELECT CASE (INF)

C	----------------------------------------------
	CASE(0) !FORCE

	SELECTCASE (INU)
	CASE(0) !N
	FAC = 1.0
	CASE(1) !KN
	FAC = 0.001
	CASE(2) !TON
	FAC = 1.0E-4 !1.019368E-4
	CASE(3) !KG
	FAC = 0.1019368
	CASE(4) !LB
	FAC = 0.224809
	CASE(5) !KIPS
	FAC = 0.224809E-3
	CASE(6) !MN
	FAC = 1.0E-6
	CASE(7) !GN
	FAC = 1.0E-9
	CASE(8) !KP
	FAC = 9.806
	ENDSELECT
C	----------------------------------------------
C	----------------------------------------------
	CASE(1) !LENGTH

	SELECTCASE (INU)
	CASE(0) !M
	FAC = 1.0
	CASE(1) !CM
	FAC = 100.0
	CASE(2) !MM
	FAC = 1000.0
	CASE(3) !KM
	FAC = 1.0E-3
	CASE(4) !INCH
	FAC = 39.370079
	CASE(5) !FT
	FAC = 3.28084
	ENDSELECT
C	----------------------------------------------
C	----------------------------------------------
	CASE(2) !ROTATION

	SELECTCASE (INU)
	CASE(0) !RAD
	FAC = 1.0
	CASE(1) !DEG
	FAC = 57.29577951
	ENDSELECT
C	----------------------------------------------
C	----------------------------------------------
	CASE(3) !MASS

	SELECTCASE (INU)
	CASE(0) !KG
	FAC = 1.0
	CASE(1) !G
	FAC = 1000.0
	CASE(2) !TON
	FAC = 0.001
	CASE(3) !SLUG
	FAC = 0.068522
	CASE(4) !LB
	FAC = 2.204623
	ENDSELECT
C	----------------------------------------------
C	----------------------------------------------
	CASE(4) !TIME

	SELECTCASE (INU)
	CASE(0) !S
	FAC = 1.0
	CASE(1) !MS
	FAC = 1000.0
	CASE(2) !MIN
	FAC = 0.01666667
	CASE(3) !HOURS
	FAC = 2.777777777E-4
	CASE(4) !DAY
	FAC = 2.3148148E-5
	ENDSELECT
C	----------------------------------------------


	END SELECT

	IF(I.EQ.1) FACO = FAC
	IF(I.EQ.2) FACN = FAC

	ENDDO

	UNITFAC = FACN/FACO

	RETURN

	END
C	=============================================================================
C	=============================================================================
C	=============================================================================
      SUBROUTINE DIROPEN(LUN,LRECL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----OPEN THE UNFORMATTED DIRECT ACCESS FILE "FIN.FXT"-----
C        LRECL = NUMBER OF BYTES IN A RECORD
      
      IREC = 0
100   IREC = IREC + 4
      IF(IREC.LT.LRECL) GOTO 100
      IREC = IREC + INT(0.1*IREC)

      OPEN (LUN,STATUS='SCRATCH',ACCESS='DIRECT',
     1      FORM='UNFORMATTED',RECL=IREC)

C
      RETURN
      END
C	=============================================================================
C	=============================================================================
C	=============================================================================
      SUBROUTINE PUTREAL(NFL,AA,NM,NI,NN,IND)  !REAL DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION AA(1),NM(NN),NI(NN)
	ALLOCATABLE AA1(:)

	
	N   = NM(1)
	NV  = NI(1)
      NRC = NDEX(NM,NI,NN,'INDXR')

	ALLOCATE(AA1(N))
	
	SELECT CASE(IND)
	CASE(0)
	AA(1:N) = 0
      READ(NFL,REC=NRC) AA1
	AA(1:NV) = AA1(1:NV)

	CASE(1)
      READ(NFL,REC=NRC) AA1
	AA1(1:NV) = AA(1:NV)
      WRITE(NFL,REC=NRC) AA1

	CASE(2)
	AA1(1:N ) = 0.0D0
	AA1(1:NV) = AA(1:NV)
      WRITE(NFL,REC=NRC) AA1

	ENDSELECT

	DEALLOCATE(AA1)
	
	RETURN
	END

C	=============================================================================
C	=============================================================================
C	=============================================================================
      SUBROUTINE PUTINTG(NFL,IA,NM,NI,NN,IND)  !INTEGER DATA
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION IA(1),NM(NN),NI(NN)
	ALLOCATABLE IA1(:)

	N   = NM(1)
	NV  = NI(1)
      NRC = NDEX(NM,NI,NN,'INDXR')

	ALLOCATE(IA1(N))
	
	SELECT CASE(IND)
	CASE(0)
	IA(1:N) = 0
      READ(NFL,REC=NRC) IA1
	IA(1:NV) = IA1(1:NV)

	CASE(1)
      READ(NFL,REC=NRC) IA1
	IA1(1:NV) = IA(1:NV)
      WRITE(NFL,REC=NRC) IA1

	CASE(2)
	IA1(1:N ) = 0.0D0
	IA1(1:NV) = IA(1:NV)
      WRITE(NFL,REC=NRC) IA1

	ENDSELECT

	DEALLOCATE(IA1)
	
	RETURN
	END

C	=============================================================================
C	=============================================================================
C	=============================================================================
      SUBROUTINE RELCALL (NAME,VALV,NM,NI,NN,NV,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4),VALV(1),NM(NN),NI(NN)
      COMMON /CONSTREL/ MTOTR(10),NCONR(5000),RCONDT(1)

	N   = NM(1)
	NRC = NI(1) - 1
	DO I = 2,NN
	M = 1
	DO J = 1,I-1
	M = M*NM(J)
	ENDDO
	NRC = NRC + M*(NI(I)-1)
	ENDDO

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MTOTR(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONR)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NCONR(I+7)+NRC
      NR = NCONR(I+4)
      NC = NCONR(I+5)
      NP = NCONR(I+6)
	CALL RCALL(RCONDT(NA),VALV,NV,IND)
  900 RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================


      SUBROUTINE INTCALL (NAME,IALV,NM,NI,NN,NV,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 NAME
      DIMENSION NAME(4),INAME(4),IALV(1),NM(NN),NI(NN)
      COMMON /CONSTINT/ MTOTI(10),NCONI(5000),ICONDT(1)

	N   = NM(1)
	NRC = NI(1) - 1
	DO I = 2,NN
	M = 1
	DO J = 1,I-1
	M = M*NM(J)
	ENDDO
	NRC = NRC + M*(NI(I)-1)
	ENDDO

C-----LOCATE AND RETURN PROPERTIES OF ARRAY -----
      NA = 0
	NUMA  = MTOTI(3)
      CALL ICONI(NAME,INAME)
      I = IFINDL(INAME,0,NUMA,NCONI)
      IF(I.EQ.0) GO TO 900
C-----RETURN ARRAY PROPERTIES -----
      NA = NCONI(I+7)+NRC
      NR = NCONI(I+4)
      NC = NCONI(I+5)
      NP = NCONI(I+6)
	CALL ICALL(ICONDT(NA),IALV,NV,IND)
  900 RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE RCALL (FA,VALV,NN,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION FA(1),VALV(1)

	DO I = 1,NN
	
	SELECTCASE(IND)
	
	CASE(0)
      VALV(I) = FA(I)

	CASE(1)
      FA(I) = VALV(I)

	CASE(2)
      FA(I) = FA(I) + VALV(I)

	CASE(3)
      FA(I) = FA(I) - VALV(I)

	CASE(4)
      FA(I) = FA(I) * VALV(I)

	CASE(5)
      FA(I) = FA(I) / VALV(I)

	CASE(6)
      FA(I) = VALV(I)/FA(I)

	CASE(7)
      FA(I) = ABS(VALV(I))

	CASE(10)
      FA(I) = 0.0D0

	ENDSELECT

	ENDDO


      RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================

      SUBROUTINE ICALL (IA,IALV,NN,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION IA(1),IALV(1)

	DO I = 1,NN
	
	SELECTCASE(IND)
	
	CASE(0)
      IALV(I) = IA(I)

	CASE(1)
      IA(I) = IALV(I)

	CASE(2)
      IA(I) = IA(I) + IALV(I)

	CASE(3)
      IA(I) = IA(I) - IALV(I)

	CASE(4)
      IA(I) = IA(I) * IALV(I)

	CASE(7)
      IA(I) = ABS(IALV(I))

	CASE(10)
      IA(I) = 0

	ENDSELECT

	ENDDO


      RETURN
      END


C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      FUNCTION MAPPING(MP,L,N)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION MP(1)

	MAPPING = 0
	DO I = 1,N
	IF(MP(I).EQ.L) THEN
	MAPPING = I
	RETURN
	ENDIF 
	ENDDO


      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      FUNCTION TESTC (OP,NOP,NP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------
C     RETURN TRUE IF STRINGS MATCH    -
C	--------------------------------------

      LOGICAL TESTC
      CHARACTER*1 OP(8),NOP(8)
      TESTC = .FALSE.
      DO 10 I=1,NP
      IF(NOP(I).NE.OP(I)) RETURN
   10 CONTINUE
      TESTC = .TRUE.


      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CALLINTG(IA,IAX,NDX,NM,NI,NN,OPTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPTN
C	--------------------------------------
C     DATA OPERATION    -
C	--------------------------------------
      DIMENSION NM(1),NI(1),IAX(1),IA(1)

      NDR = (NDX - 1) + NDEX(NM,NI,NN,'INDXN')

      N = NM(1)
      SELECTCASE(OPTN)
      CASE('READ')
        IA(1:N) = IAX(NDR:NDR+N-1)
      CASE('WRIT')
        IAX(NDR:NDR+N-1) = IA(1:N)
      ENDSELECT


      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      SUBROUTINE CALLREAL(AA,AAX,NDX,NM,NI,NN,OPTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPTN
C	--------------------------------------
C     DATA OPERATION    -
C	--------------------------------------
      DIMENSION NM(1),NI(1),AAX(1),AA(1)

      NDR = (NDX - 1) + NDEX(NM,NI,NN,'INDXN')

      N = NM(1)
      SELECTCASE(OPTN)
      CASE('READ')
        AA(1:N) = AAX(NDR:NDR+N-1)
      CASE('WRIT')
        AAX(NDR:NDR+N-1) = AA(1:N)
      ENDSELECT


      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================
      FUNCTION NDEX(NM,NI,NN,OPTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*5 OPTN
C	--------------------------------------
C     RETURN INDEX OF DATA    -
C	--------------------------------------
      DIMENSION NM(1),NI(1)

	N   = NM(1)
	NV  = NI(1)
	NRC = 1
	DO I = 2,NN
	M = 1
	DO J = 1,I-1
	M = M*NM(J)
	ENDDO
	NRC = NRC + M*(NI(I)-1)/N
	ENDDO
	
	NDEX = -1
	SELECTCASE(OPTN)
	CASE('INDXN') !RETURN INDEX NUMBER
	  NDEX = N*(NRC-1) + 1
	CASE('INDXR') !RETURN RECORD NUMBER
	  NDEX = NRC
	ENDSELECT


      RETURN
      END

C	=======================================================================
C	=== FREE FIELD INPUT & OUTPUT UTILITY === SONGSAK NOV2007 =============
C	=======================================================================




