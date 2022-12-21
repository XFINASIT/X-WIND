C
C=====================================================================
      BLOCK DATA GAUSSP
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     INITIALIZES COMMON BLOCKS INOU,SOLU,ITER,EIGN
C     SETS SAMPLING POINT POSITIONS AND WEIGHTING FACTORS FOR
C     GAUSS NUMERICAL INTEGRATION ORDERS 1,2,3,4
C	------------------------------------------
C     VARIABLES IN COMMON BLOCK /GAUS/
C     GLOC(4,4)    = SAMPLING POINT POSITIONS
C     GWT (4,4)    = WEIGHTING FACTORS
C     NGR          = ORDER OF NUMERICAL INTEGRATION IN R DIRECTION
C     NGS          = ORDER OF NUMERICAL INTEGRATION IN S DIRECTION
C     NGT          = ORDER OF NUMERICAL INTEGRATION IN T DIRECTION
C     ------------------------------------------------------------
      COMMON /INOU/  ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
C      COMMON /SOLU/  NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +               NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /EIGN/  NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
C
      DATA ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,IFPR,IFPL
C	SUNIL 01/05/03 CHANGED ITI TO 1
C     1     /5,6,7,1,4,8,9,20*0/
     1     /1,6,7,1,4,8,9,20*0/
C
      DATA NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,NRED
     1     /4*0,0.0,3*0,0.0,2*0/
      DATA NSTEP,DLMAX,NPRIN,NDRAW,KONEQ,NIREF,ITEMAX,RTOL,ITOPT,
     1     NUMREF,NUMITE,ITETOT,RHO,RHOP,ICONV,NOLIN
     2     /1,1.0,1,1,2,15,15,1.0E-6,1,3*0,2*0.0,2*0/
C
      DATA NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG
     1     /6*0,2*0.0,2*0/
C
      DATA ((GLOC(I,J),I=1,10),J=1,8) /
     1            10*0.,
     2           -.57735026918963,  .57735026918963,  8*0.,
     3           -.77459666924148,  .00000000000000,  .77459666924148,
     3           7*0.,
     4           -.86113631159405, -.33998104358486,  .33998104358486,
     4            .86113631159405,  6*0.,
     5           -.90617984593866, -.53846931010568,  .00000000000000,
     5            .53846931010568,  .90617984593866,  5*0.,
     6           -.93246951420315, -.66120938646627, -.23861918608320,
     6            .23861918608320,  .66120938646627,  .93246951420315,
     6            4*0.,
     7           -.94910791234276, -.74153118559939, -.40584515137740,
     7            .00000000000000,  .40584515137740,  .74153118559939,
     7            .94910791234276,  3*0.,
     8           -.96028985649754, -.79666647741363, -.52553240991633,
     8           -.18343464249565,  .18343464249565,  .52553240991633,
     8            .79666647741363,  .96028985649754,  2*0./
      DATA ((GLOC(I,J),I=1,10),J=9,10) /
     9           -.96816023950763, -.83603110732664, -.61337143270059,
     9           -.32425342340381,  .00000000000000,  .32425342340381,
     9            .61337143270059,  .83603110732664,  .96816023950763,
     9            .00000000000000,
     1           -.97390652851717, -.86506336668899, -.67940956829902,
     1           -.43339539412925, -.14887433898163,  .14887433898163,
     1            .43339539412925,  .67940956829902,  .86506336668899,
     1            .97390652851717/
      DATA ((GWT(I,J),I=1,10),J=1,8) /
     1           2.00000000000000,  9*0.,
     2           1.00000000000000, 1.00000000000000, 8*0.,
     3            .55555555555556,  .88888888888889,  .55555555555556,
     3            7*0.,
     4            .34785484513745,  .65214515486255,  .65214515486255,
     4            .34785484513745,  6*0.,
     5            .23692688505619,  .47862867049937,  .56888888888889,
     5            .47862867049937,  .23692688505619,  5*0.,
     6            .17132449237917,  .36076157304814,  .46791393457269,
     6            .46791393457269,  .36076157304814,  .17132449237917,
     6            4*0.,
     7            .12948496616887,  .27970539148928,  .38183005050512,
     7            .41795918367347,  .38183005050512,  .27970539148928,
     7            .12948496616887,  3*0.,
     8            .10122853629038,  .22238103445337,  .31370664587789,
     8            .36268378337836,  .36268378337836,  .31370664587789,
     8            .22238103445337,  .10122853629038,  2*0./
      DATA ((GWT(I,J),I=1,10),J=9,10) /
     9            .08127438836157,  .18064816069486,  .26061069640294,
     9            .31234707704000,  .33023935500126,  .31234707704000,
     9            .26061069640294,  .18064816069486,  .08127438836157,
     9            .00000000000000,
     1            .06667134430869,  .14945134915058,  .21908636251598,
     1            .26926671931000,  .29552422471475,  .29552422471475,
     1            .26926671931000,  .21908636251598,  .14945134915058,
     1            .06667134430869/
      DATA NGR,NGS,NGT  /2,2,1/
C
      END
C
C=====================================================================
      SUBROUTINE CLEARA (A,N)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------
C     PROGRAM TO CLEAR FLOATING POINT ARRAY A
C     ---------------------------------------
      DIMENSION A(*)
      
      A(1:N) = 0.0D0

      RETURN
      END
C
C=====================================================================
      SUBROUTINE CLEARI (IA,N)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------
C     PROGRAM TO CLEAR INTEGER ARRAY IA(N)
C     ------------------------------------
      DIMENSION IA(*)
      
      IA(1:N) = 0

      RETURN
      END
C
C=====================================================================
      SUBROUTINE ERRORS (K,NN,NUM,LPLACE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------
C     ERROR MESSAGES
C     --------------
C	SUNIL 07/01/01 ADD FOLLOWING LINE
	CHARACTER*10 LPLACE

      COMMON /INOU/  ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
C
C	The next line changed to the second by NguyenDV, Jan27,2003
C     IF (K.NE.7) WRITE (ITO,3000)
	IF (K.NE.7.OR.K.EQ.26) WRITE (ITO,3000)
      IF (K.NE.7.OR.K.EQ.26) WRITE (10,3000)
C
      IF ((K.LE.9 ).OR.(K.GE.20)) WRITE (ITO,3100) LPLACE,NUM
      IF ((K.LE.9 ).OR.(K.GE.20)) WRITE (10,3100) LPLACE,NUM
      IF ((K.GE.10).AND.(K.LE.19)) WRITE (ITO,3100) LPLACE
      IF ((K.GE.10).AND.(K.LE.19)) WRITE (10,3100) LPLACE
C
C     POINTER (26,27) ADDED BY NguyenDV, January 27,2003 
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
     1      ,21,22,23,24,26,27), K
 1    WRITE (ITO,100)
      WRITE (10,100)
      GOTO 50
 2    WRITE (ITO,200)
      WRITE (10,200)
      GOTO 50
 3    WRITE (ITO,300)
      WRITE (10,300)
      GOTO 50
 4    WRITE (ITO,400)
      WRITE (10,400)
      GOTO 50
 5    WRITE (ITO,500)
      WRITE (10,500)
      GOTO 50
 6    WRITE (ITO,600)
      WRITE (10,600)
      GOTO 50
 7    WRITE (ITO,700)
      WRITE (10,700)
      GOTO 75
 8    WRITE (ITO,800)
      WRITE (10,800)
      GOTO 50
 9    WRITE (ITO,900)
      WRITE (10,900)
      GOTO 25
 10   WRITE (ITO,1000) LPLACE,NN,NUM
      WRITE (10,1000) LPLACE,NN,NUM
      GOTO 25
 11   WRITE (ITO,1100) NN,NUM
      WRITE (10,1100) NN,NUM
      GOTO 25
 12   CALL ERRORK (NN) !!! WRITE (ITO,1200) NN,NUM
      GOTO 25
 13   WRITE (ITO,1300) NUM,NN
      WRITE (10,1300) NUM,NN
      GOTO 25
 14   WRITE (ITO,1400)
      WRITE (10,1400)
      GOTO 25
 15   WRITE (ITO,1500) NUM
      WRITE (10,1500) NUM
      GOTO 25
 16   WRITE (ITO,1600) NN,NUM
      WRITE (10,1600) NN,NUM
      GOTO 25
 17   WRITE (ITO,1700) NN,NUM
      WRITE (10,1700) NN,NUM
      GOTO 25
 18   WRITE (ITO,1800) NUM
      WRITE (10,1800) NUM
      GOTO 25
 19   WRITE (ITO,1900) NUM
      WRITE (10,1900) NUM
      GOTO 75
 20   WRITE (ITO,2000) NUM
      WRITE (10,2000) NUM
 21   WRITE (ITO,2100)
      WRITE (10,2100)
      GOTO 50
 22   WRITE (ITO,2200)
      WRITE (10,2200)
      GOTO 50
 23   WRITE (ITO,2300)
      WRITE (10,2300)
      GOTO 50
 24   WRITE (ITO,2400) NN,NUM
      WRITE (10,2400) NN,NUM
      GOTO 25
C
 25   CALL SYSOUT

C	THE NEXT FOUR LINES ADDED BY NguyenDV, January 27,2003
 26   WRITE (ITO,2600)
      WRITE (10,2600)
      GOTO 75
 27   WRITE (ITO,2700) NN,NUM
      WRITE (10,2700) NN,NUM
      GOTO 25
C
 50   NN   = NN+1
C	SUNIL 09/02/01 ADDED NEXT TWO LINE
 30	STOP
 75   RETURN
C
 3000 FORMAT (/,27X,23(1H*)/27X,1H*,21X,1H*/
     1        27X,23H*  ERROR  MESSAGE(S)  */
     2        27X,1H*,21X,1H*,/27X,23(1H*)/)
 3100 FORMAT (23X,9HINPUT OF ,A10,10H / NUMBER ,I3/23X,32(1H-))
 3200 FORMAT (26X,15HCOMPUTATION OF ,A10/26X,25(1H-))
C
 100  FORMAT (15X,'***         IMPOSSIBLE CONTROL NUMBER         ***'/)
 200  FORMAT (15X,'***  JOINT NO. BIGGER THAN MAXIMUM SPECIFIED  ***'/)
 300  FORMAT (15X,'*** ELEMENT NO. BIGGER THAN MAXIMUM SPECIFIED ***'/)
 400  FORMAT (15X,'*** SECTION NO. BIGGER THAN MAXIMUM SPECIFIED ***'/)
 500  FORMAT (15X,'***     BOUNDARY CONDITION UNEQUAL 1 OR O     ***'/)
 600  FORMAT (15X,'***  DEGREE OF FREEDOM FIXED OR OUT OF BOUND  ***'/)
 700  FORMAT (15X,'*** WARNING, LOADED SUPPORT, LOAD IGNORED     ***'/)
 800  FORMAT (15X,'***   IMPROPER TYPE OF DATA (INTEGER/REAL)    ***'/)
 900  FORMAT (15X,'***       END OF FILE DATIN ENCOUNTERD        ***'/)
 1000 FORMAT (22X,23HNOT ENOUGH STORAGE FOR ,A10,/22X,33(1H-)/
     1         16X,44HALLOCATED FIELD LENGTH AS GIVEN ON JOB CARD ,I6/
     2          16X,44HREQUIRED FIELD LENGTH (DECIMAL)). . . . . . ,I6//
     2)
 1100 FORMAT (22X,23HNOT ENOUGH STORAGE FOR ,A10/22X,33(1H-)/
     1         19X,37HMAXIMUM COLUMN LENGTH . . . . MK    =,I5/
     2          19X,37HMAXIMUM BLOCK SIZE PERMITTED  ISTOR =,I5//)
 1200 FORMAT (15X,48H*** STIFFNESS MATRIX NOT POSITIVE DEFINITE. .***///
     1         21X,30HNONPOSITIVE PIVOT FOR EQUATION,I4//
     2         21X,14HPIVOT = . . . ,E20.12)
 1300 FORMAT (15X,48H*** STURM SEQUENCE CHECK FAILED BECAUSE OF   ***//
     1         15X,40H*** MULTIPLIER GROWTH FOR COLUMN NUMBER ,I4,4H ***
     2//         15X,24H*** MULTIPLIER . . . . =,E20.8,4H ***)
 1400 FORMAT (15X,48H*** NC LARGER THAN NO OF DEGREES OF FREEDOM  ***///
     +)
 1500 FORMAT (15X,39H*** ZERO JACOBIAN DETERMINANT FOR ELEM.,I3,3X,3H***
     +)
 1600 FORMAT (15X,48H** RESIDUAL LOAD LARGER THAN INCREMENTAL LOAD **//
     1         15X,43H** NUMBER OF CURRENT LOAD STEP . . KSTEP = ,I2,3H
     2**/        15X,43H** NUMBER OF CURRENT ITERATION . . . ITE = ,I2,3
     2H **)
 1700 FORMAT (15X,21H** NO CONVERGENCE FOR,I3,24H PERMITTED ITERATIONS *
     1*/        15X,43H** NUMBER OF CURRENT LOAD STEP . . KSTEP = ,I2,3H
     1 **)
 1800 FORMAT (15X,28H** INADMISSIBLE UNIT NUMBER(,I3,17H) FOR I/O FILE *
     +*)
 1900 FORMAT (15X,'** WARNING,COORD.OF NODES >',I4,' NOT SPECIFIED **'/)
 2000 FORMAT (15X,'***  MAXIMUM NUMBER OF CONNECTED NODES >',I3,' ***'/)
 2100 FORMAT (15X,'***  AXES SET N.BIGGER THAN MAXIMUM SPECIFIED  ***'/)
 2200 FORMAT (15X,'***   LOCAL X AND Y AXES ARE NOT ORTHOGONAL    ***'/)
 2300 FORMAT (15X,'***  NODE RESTRAINED IN A DIFFERENT AXES SET   ***'/)
 2400 FORMAT (10X,'*** GENERATED N0. OF JOINT LOADS ',I5,' > SPECIFIED',
	1        I5,' ***'/)
C
C	THE NEXT LINES (2600,2700) ADDED BY NguyenDV, January 27,2003
 2600 FORMAT (15X,'*** WARNING, DAMPER PLACE AT SUPPORT, DAMPER IGNORED 
     1         ***'/)
 2700 FORMAT (10X,'*** GENERATED N0.OF JOINT DAMPERS ',I5,' > SPECIFIED',
	1        I5,' ***'/) 
C
      END
C
C=====================================================================
      SUBROUTINE MATOUT (A,IA,JA,K1,K2,K3,K4,K5,MODE,TITLE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------
C     TO PRINT THE MATRIX A(IA,JA) WITH THE K SPECIFIED
C               K1 = NUMBER OF COLUMNS
C               K2 = COLUMN WIDTH
C               K3 = TYPE OF FORMAT  (I,O,A,R / F,G,E)
C               K4 = NUMBER OF DIGITS AFTER DECIMAL POINT
C               K5 = NUMBER OF PRECEEDING BLANKS
C     ---------------------------------------------------
      CHARACTER*60 FORM(4)
C	SUNIL 16/01/01 ADD FOLLOWING TWO LINES
	CHARACTER*1 K3
	CHARACTER*10 TITLE
C
      COMMON /INOU/  ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
      DIMENSION A(IA,1)
C	SUNIL 09/01/01 REMOVE FOLLOWING LINE
C      INTEGER*4 TITLE

      FORM(1)='         '
      FORM(2)='         '
      FORM(3)='         '
      FORM(4)='         '
      WRITE (ISO,1000) TITLE,IA,JA
      I5 = K5+K2/2+3
      I6 = K2-3
      IF (I6.LT.0) I6=0
      ENCODE (20,2000,FORM(3)) I5,K1,I6
c      WRITE (FORM(3),2000) I5,K1,I6
      GOTO (100,200,500), MODE
C      CALL GOTOER!NON-EXIST SUBROUTINE
C
C     -----------------------------------------
C     FORMAT (I4,*X,***)     FOR INTEGER OR BCD
C     -----------------------------------------
c 100  WRITE (FORM(1),3000) K5,K1,K3,K2
 100  ENCODE (20,3000,FORM(1)) K5,K1,K3,K2
      GOTO 300
C     ----------------------------------------------
C     FORMAT (1X,I3,*X,***.*/)    FOR FLOATING POINT
C     ----------------------------------------------
c 200  WRITE (FORM(1),5000) K5,K1,K3,K2,K4
 200  ENCODE (20,5000,FORM(1)) K5,K1,K3,K2,K4
C
 300  NB = FLOAT(JA)/FLOAT(K1) + 0.999
      J1 = 1
      J2 = K1
      IF (J2.GT.JA) J2 = JA
C
      DO 310  IB=1,NB
      WRITE (ISO,FORM(3)) (J,J=J1,J2)
      WRITE (ISO,2100)
      DO 320  I=1,IA
 320  WRITE (ISO,FORM(1)) I,(A(I,J),J=J1,J2)
      J1 = J1+K1
      J2 = J2+K1
 310  IF (J2.GT.JA) J2 = JA
      RETURN
C     ---------------------------------------------------
C     MATRIX A IN UPPER TRIANGULAR FORM (SIZE IA(IA+1)/2)
C     ---------------------------------------------------
 500  IF (JA.GT.K1) GOTO 600
      J1  = 1
      J2  = JA
      WRITE (ISO,FORM(3)) (J,J=J1,J2)
      WRITE (ISO,2100)
      DO 550  I=1,JA
      I5 = (I-1)*K2+K5
      I1 = JA+1-I
      ENCODE (20,5000,FORM(1)) I5,I1,K3,K2,K4
c      WRITE (FORM(1),5000) I5,I1,K3,K2,K4
      WRITE (ISO,FORM(1)) I, (A(K,1),K=J1,J2)
      J1 = J2+1
 550  J2 = J2+I1-1
      RETURN
C
c 600  WRITE (FORM(1),5000) K5,K1,K3,K2,K4
 600  ENCODE (20,5000,FORM(1)) K5,K1,K3,K2,K4
C
      KS = 1
      KE = JA
      WRITE (ISO,5500)
      DO 650  I=1,IA
      DO 630  KL=KS,KE,K1
      KLP = KL+K1-1
      IF (KLP.GT.KE) KLP = KE
 630  WRITE (ISO,FORM(1)) KL,(A(J,1),J=KL,KLP)
      WRITE (ISO,6000)
      KS = KE+1
 650  KE = KE+JA-I
C
C
 1000 FORMAT (////1X,A10,' MATRIX',I4,' ROWS',I4,' COLUMNS'/1X,38(1H-))
 2000 FORMAT ('(//',I2,2HX,,I2,4H(I3,,I2,3HX)))
 2100 FORMAT (/)
 3000 FORMAT (4H(I4,,I2,2HX,,I2,A1,I2,1H))
 5000 FORMAT (7H(1X,I4,,I2,1HX,I2,A1,I2,1H.,I2,1H))
 5500 FORMAT (///)
 6000 FORMAT (1H0)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE VECTOR (A,B,C,FACT,VNORM,IDIM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PERFORMES VECTOR OPERATIONS
C     ----------------------------------------------------------------
      DIMENSION A(100),B(100),C(100)
C     --------------------------------
C     VECTOR ADDITION C(I) = A(I)+B(I)
C     --------------------------------
      ENTRY VECADD(A,B,C,FACT,VNORM,IDIM)
      C(1:IDIM) =  A(1:IDIM)+B(1:IDIM)!*FACT
      RETURN
      DO 100  I=1,IDIM
 100  C(I) = A(I)+B(I)!*FACT
      RETURN
C     -----------------------------------
C     VECTOR SUBTRACTION C(I) = A(I)-B(I)
C     -----------------------------------
      ENTRY VECSUB(A,B,C,FACT,VNORM,IDIM)
      C(1:IDIM) =  A(1:IDIM)-B(1:IDIM)
      RETURN
      DO 200  I=1,IDIM
 200  C(I) = A(I)-B(I)
      RETURN
C     --------------------------------------------
C     MULTIPLIE VECTOR BY SCALAR C(I) = FACT*A(I)
C     --------------------------------------------
      ENTRY VECMUL(A,B,C,FACT,VNORM,IDIM)
      C(1:IDIM) =  FACT*A(1:IDIM)
      RETURN
      DO 300  I=1,IDIM
 300  C(I) = FACT*A(I)
      RETURN
C     ------------------------------------
C     VECTOR NORM VNORM = SUM OF A(I)*B(I)
C     ------------------------------------
      ENTRY VENORM(A,B,C,FACT,VNORM,IDIM)
      VNORM = 0.
      DO 400  I=1,IDIM
 400  VNORM = VNORM + A(I)*B(I)
      RETURN
C     ------------------------------------
C     ADDITIVE NORM VNORM = [A(I)+B(I)]**2
C     ------------------------------------
      ENTRY ADDNRM(A,B,C,FACT,VNORM,IDIM)
      VNORM = 0.
      DO 500  I=1,IDIM
      SUM = A(I)+B(I)
 500  VNORM = VNORM + SUM*SUM
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE VECMAX (V,NEQ,VMAX,IPOS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     FINDS THE POSITION AND VALUE OF THE ABSOLUTE LARGEST COEFF.IN V
C	---------------------------------------------------------------
C     V(NEQ)    = VECTOR TO BE SCANED FOR MAXIMUM
C     NEQ       = SIZE OF VECTOR
C     VMAX      = ABSOLUTE LARGEST COEFFICIENT
C     IPOS      = POSITION OF LARGEST COEFFICIENT
C     ---------------------------------------------------------------
      DIMENSION V(1)
C
      AMAX = 0.
      DO 100  IEQ=1,NEQ
      VABS = ABS(V(IEQ))
      IF (VABS.LT.AMAX)  GOTO 100
      AMAX = VABS
      VMAX = V(IEQ)
      IPOS = IEQ
 100  CONTINUE
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MATMULT (A,B,C,IA,JA,JB,CONST,MODE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------
C     MATRIX MULTIPLICATION C(IA,JB) = C + A(IA,JA)*B(JA,JB)
C     ------------------------------------------------------
      DIMENSION A(1),B(1),C(1)
      GOTO (100,200,300,400),MODE
C     --------------------------
C     FULL A AND FULL B MATRICES
C     --------------------------
 100  KCC = 1
      KBB = -JA
      DO 190  J=1,JB
      KBB = KBB+JA
      DO 190  I=1,IA
      KA = I-IA
      KB = KBB
      V = 0.
C
      DO 150  K=1,JA
      KB = KB+1
      KA = KA+IA
 150  V = V + A(KA)*B(KB)
      C(KCC) = C(KCC) + V
 190  KCC = KCC + 1
      RETURN
C     ----------------------------
C     SPARSE A AND FULL B MATRICES
C     ----------------------------
 200  KA = 1
      DO 290  J=1,JA
      DO 290  I=1,IA
      V = A(KA)
      IF (V.EQ.0.0)  GOTO 290
      KC = I-IA
      KB = J-JA
C
      DO 250  K=1,JB
      KC = KC+IA
      KB = KB+JA
 250  C(KC) = C(KC) + V*B(KB)
 290  KA = KA+1
      RETURN
C     ----------------------------
C     FULL A AND SPARSE B MATRICES
C     ----------------------------
 300  KBB = 1
      KCC = -IA
      DO 390  J=1,JB
      KCC = KCC+IA
      KAA = -IA
      DO 390  I=1,JA
      V = B(KBB)
      KC = KCC
      KAA = KAA+IA
      KA = KAA
      IF (V.EQ.0.0)  GOTO 390
C
      DO 350  K=1,IA
      KC = KC+1
      KA = KA+1
 350  C(KC) = C(KC) + V*A(KA)
 390  KBB = KBB+1
      RETURN
C     -----------------------------------------------
C     SPARSE A AND FULL B,C UPPER TRIANGULAR ROW-WISE
C     -----------------------------------------------
 400  KA = 1
      KC = -JA-1
      DO 490  J=1,JA
      KC = KC+JA-J+2
      DO 490  I=1,IA
      IF (A(KA).EQ.0.0)  GOTO 490
      V = A(KA)*CONST
      KCC = KC
      KB = KA-IA
C
      DO 450  K=J,JA
      KCC = KCC+1
      KB = KB+IA
 450  C(KCC) = C(KCC) + V*B(KB)
 490  KA = KA+1
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MATRAN (A,B,IS,JS,MODE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------
C     TRANSPOSES MATRIX A
C     -------------------
      DIMENSION A(IS,1),B(JS,1)
C     -------------------------------------------------
C     MATRIX A(IS,JS) TRANSPOSED AND STORED AS B(JS,IS)
C     -------------------------------------------------
      IF (MODE.NE.1)  GOTO 10
      DO 20  J = 1,JS
      DO 20  I = 1,IS
  20  B(J,I) = A(I,J)
      RETURN
C     --------------------------------------------
C     SQUARE MATRIX A TRANSPOSED AND REPLACED IN A
C     --------------------------------------------
  10  DO 30  J = 1,JS
      DO 30  I = 1,IS
      AA = A(I,J)
      A(I,J) = A(J,I)
  30  A(J,I) = AA
      RETURN
      END
C
C=====================================================================
C
C=====================================================================
C
C=====================================================================
	SUBROUTINE GIDORDER(IGIDMEM,NELE,MEMBA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION IGIDMEM(NELE)

	DO  IGIDE = 1,NELE
		MEMGID = IGIDMEM(IGIDE) 
		IF(MEMGID.EQ.MEMBA) THEN
			MEMBA = IGIDE
			RETURN
		ENDIF
	ENDDO
	MEMBA = 0

      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE CLEAR (M,N,A)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------
C     PROGRAM TO CLEAR FLOATING POINT ARRAY A(I,J)=0.0
C     ------------------------------------------------
      DIMENSION A(M,*)
C
	DO 10 I=1,M
	DO 10 J=1,N
10	A(I,J)=0.0
	RETURN
	END
C
C	=====================================================================

      SUBROUTINE ADDROW(AV,BV,CV,FAC,NN)
	IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
C     ---------------------------------
C     VECTOR ADDITION CV(I)=AV(I)+BV(I)
C     ---------------------------------
      DIMENSION AV(1),BV(1),CV(1)
C     
      DO 10 I=1,NN
C      AV(I)=0.0
   10 CV(I)=AV(I)+FAC*BV(I)
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================

      SUBROUTINE EDOFLPOS(IDOF,NEFLNK,NNF,IGPOS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION IDOF(9),IGPOS(NNF),LEFLNK(NNF)

C	----------------------------------------
C     UNPACK VARIABLE LINKF
C	----------------------------------------
      LINK = NEFLNK
      
      CALL EDOFLINK(LINK,NNF,LEFLNK)
      
	IGPOS(1:NNF) = 0
      DO INF = 1,NNF
        ITEST = LEFLNK(INF)
        DO IGF = 1,9
            IF(IDOF(IGF).EQ.ITEST) IGPOS(INF) = IGF
        ENDDO
      ENDDO
      
	
      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================

      SUBROUTINE EDOFLINK(NEFLNK,NNF,LEFLNK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION LEFLNK(NNF)

C     ----------------------------------------------
C     UNPACK ELEMENT DOF VARIABLE -- NEFLNK
C     ----------------------------------------------
      LINK = NEFLNK

      DO INF=1,NNF
          LEFLNK(INF) = LINK/10**(NNF-INF)
	    LINK = LINK - 10**(NNF-INF)*LEFLNK(INF)
	ENDDO
	
      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================

      SUBROUTINE NESUPPRS(MCONT,IDOF,NSF,NSN,LINKF,NNF,NEF,NNM,NELE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      DIMENSION MCONT(NEF,1),IDOF(1)
      DIMENSION LEFLNK(NNF),IDSP(NSF,NSN)

      REWIND(8)
      READ(8) IDSP
      
      CALL EDOFLINK(LINKF,NNF,LEFLNK)
      
      DO IELE = 1,NELE
          DO  INM = 1,NNM
              NOD = MCONT(INM,IELE)
              IF(NOD.LE.0) GOTO 100
              DO INF = 1,NNF
                  ISF = 0
                  DO JSF = 1,NSF
                      IF(LEFLNK(INF).EQ.IDOF(JSF)) THEN
                          ISF = JSF
                          EXIT
                      ENDIF
                  ENDDO
                  IF(ISF.NE.0) IDSP(ISF,NOD) = 1
              ENDDO
100           CONTINUE
          ENDDO
      ENDDO
  
      REWIND(8)
      WRITE(8) IDSP

      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================

      SUBROUTINE NNSUPPRS(ID,NSF,NSN,OPTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPTN
      
      DIMENSION ID(NSF,1),IDSP(NSF,NSN)

      SELECTCASE(OPTN)
      
      CASE('INIT')
          OPEN (UNIT=8,STATUS='SCRATCH',FORM='UNFORMATTED')
          IDSP = 0
          WRITE(8) IDSP

      CASE('SUPP')
          REWIND(8)
          READ(8) IDSP
          DO ISN = 1,NSN
              DO ISF = 1,NSF
                  IF(IDSP(ISF,ISN).EQ.0) ID(ISF,ISN) = 1
              ENDDO
          ENDDO
          CLOSE(UNIT=8,STATUS='DELETE')
          
      ENDSELECT

      RETURN
      END
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE INVERSE3 (AI,AI2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION AI(3,3),AI2(3,3)

	DEL=AI(1,1)*(AI(2,2)*AI(3,3)-AI(2,3)*AI(3,2))-
     1	AI(1,2)*(AI(2,1)*AI(3,3)-AI(3,1)*AI(2,3))+
     2	AI(1,3)*(AI(2,1)*AI(3,2)-AI(2,2)*AI(3,1)) 

      AI2(1,1)=(AI(2,2)*AI(3,3)-AI(2,3)*AI(3,2))/DEL
	AI2(1,2)=(-1.0)*(AI(1,2)*AI(3,3)-AI(1,3)*AI(3,2))/DEL
      AI2(1,3)=(AI(1,2)*AI(2,3)-AI(1,3)*AI(2,2))/DEL
	AI2(2,1)=(-1.0)*(AI(2,1)*AI(3,3)-AI(2,3)*AI(3,1))/DEL
      AI2(2,2)=(AI(1,1)*AI(3,3)-AI(1,3)*AI(3,1))/DEL
	AI2(2,3)=(-1.0)*(AI(1,1)*AI(2,3)-AI(1,3)*AI(2,1))/DEL
      AI2(3,1)=(AI(2,1)*AI(3,2)-AI(2,2)*AI(3,1))/DEL
	AI2(3,2)=(-1.0)*(AI(1,1)*AI(3,2)-AI(1,2)*AI(3,1))/DEL
      AI2(3,3)=(AI(1,1)*AI(2,2)-AI(1,2)*AI(2,1))/DEL

	RETURN
	END

C
C     ======================================================================
C     ======================================================================
C     ======================================================================
      SUBROUTINE INVERSE2(AI,AI2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION AI(2,2),AI2(2,2)

	DEL=AI(1,1)*AI(2,2)-AI(1,2)*AI(2,1)

	AI2(1,1)=AI(2,2)/DEL
	AI2(1,2)=(-1.0)*AI(1,2)/DEL
	AI2(2,1)=(-1.0)*AI(2,1)/DEL
	AI2(2,2)=AI(1,1)/DEL

	RETURN
	END
      
C
C     ======================================================================
C     ======================================================================
C     ======================================================================
	SUBROUTINE CANGLE(VR,VS,VM,ANG)
C	-------------------------------------
C	TO COMPUTE FOR INITIAL MATERIAL ANGLE
C	  PROGRAMMED BY GILSON
C	-------------------------------------
C	ANG  - INITIAL MATERIAL ANGLE
C	DVEC - DUMMY VECTOR
C	DUM  - DUMMY VARIABLE
C	R    - PROJECTION OF VM ON VR
C	S    - PROJECTION OF VM ON VS
C	VM   - MATERIAL BASE VECTOR
C	VR   - LOCAL BASE VECTOR IN R DIR
C	VS   - LOCAL BASE VECTOR IN S DIR
C	-------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION VR(3),VS(3),VM(3)
	DIMENSION DVEC(3)

C	------------------------------------
C	COMPUTE FOR PRINCIPAL MATERIAL ANGLE
C	------------------------------------
	CALL VENORM(VM,VS,DVEC,DUM,S,3)
	CALL VENORM(VM,VR,DVEC,DUM,R,3)
	ANG = ATAN(S/R)

	RETURN
	END
C
C     ======================================================================
C     ======================================================================
C     ======================================================================

      SUBROUTINE ERRORK (MEQ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PRINT THE ERROR MESSAGE IN CASE OF NON-POSITIVE DEFINITION
C     ----------------------------------------------------------------
	CHARACTER*6 IFREE(9)
	CHARACTER*6 IPO
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /INOU/  ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)

      COMMON A(9000000),IA(9000000)
C     ----------------------------------------------------------------
      DATA IFREE /'X-DISP','Y-DISP','Z-DISP','X-ROTA','Y-ROTA',
     +            'Z-ROTA','WARP. ','TEMP. ','POTEN.'/
C     ----------------------------------------------------------------

	K = 0
	DO ISN = 1,NSN
	DO ISF = 1,NSF
	K = K + 1
	IEQ = IA(LID+K-1)
	IF(IEQ.EQ.MEQ) THEN

	IDF = IDOF(ISF)
	IPO = IFREE(IDF)

	WRITE(ITO,1200) ISN,IPO
      WRITE(10,1200) ISN,IPO

	RETURN
	ENDIF
	ENDDO
	ENDDO


 1200 FORMAT (15X,48H*** STIFFNESS MATRIX NOT POSITIVE DEFINITE. .***///
     1         21X,30HNONPOSITIVE PIVOT FOR NODE No.,I4//
     2         21X,30HAT DEGREE OF FREEDOM . . . . .,A6)


	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
