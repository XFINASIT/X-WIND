
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE ADDMAPS(S,SS,MAPS,NWK,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      CHARACTER*3 OPER
      
      DIMENSION S(1),SS(1),MAPS(1)

	IND = -1
	IF(OPER.EQ.'ADD') IND = 0 
	IF(OPER.EQ.'SUB') IND = 1
	IF(IND.EQ.-1) RETURN
      
	SELECT CASE(IND)

C	---------------------------------------
	CASE(0)
      S(1:NWK) = S(1:NWK) + SS(MAPS(1:NWK))
C	---------------------------------------
	CASE(1)
      S(1:NWK) = S(1:NWK) - SS(MAPS(1:NWK))
C	---------------------------------------

	ENDSELECT	 
      
      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SHKLINEAS (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD,FACT)
C
C     ------------------------------------------------------------
C     EVALUATES LINEAR CONTRIBUTION TO TANGENTIAL STIFFNESS MATRIX
C
C     S(1176)  = STIFFNESS MATRIX STORED UPPER-TRIANGULAR ROW-WISE
C     DR(64)   = ELASTO-PLASTIC RIGIDITY MATRIX STORED COLUMN-WISE
C     B(240)   = STRAIN-DISPLACEMENT MATRIX
C     BDRL(48) = DERIVATIVE OPERATORS FOR DRILLING STRAIN
C     NNO      = NUMBER OF NODES FOR ELEMENT
C     NEF      = NUMBER OF D.O.F FOR ELEMENT
C     IPEL     = SECTION PLASTICITY INDICATOR (1=EL,2=EL-PL)
C     ------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      DIMENSION S(*),DR(64),B(240),BDRL(48)
C
      N=1
      I=4
      FDRL=FACT*DR(19)
C     -----------------------------
C     TRANSVERSE SHEAR CONTRIBUTION
C     -----------------------------
      DO 30 NROW=1,NEF
      A1=DR(55)*B(I) + DR(56)*B(I+1)
      A2=DR(63)*B(I) + DR(64)*B(I+1)
      A3=FDRL*BDRL(NROW)
      J=I
      DO 20 NCOL=NROW,NEF
      S(N)=S(N)+A1*B(J)+A2*B(J+1)+A3*BDRL(NCOL)
      N=N+1
   20 J=J+5
   30 I=I+5
      LROW=NEF
      I=1
      N1=1
      N2=3*NEF-2
      IF (IPEL.EQ.2.OR.MTMOD.EQ.2) GO TO 120
C     -----------------------
C     ELASTIC RIGIDITY MATRIX
C     -----------------------
      DO 100 IR=1,NNO
      DO 90 IRR=1,3
      A1=DR(1)*B(I)+DR(2)*B(I+1)
      A2=DR(2)*B(I)+DR(10)*B(I+1)
      A3=DR(19)*B(I+2)
      A4=DR(28)*B(I+15)+DR(29)*B(I+16)
      A5=DR(29)*B(I+15)+DR(37)*B(I+16)
      A6=DR(46)*B(I+17)
      J=I
      IF=4-IRR
C     ------------------------------------------------
C     UPPER TRIANGULAR PART OF DIAGONAL 3*3 PARTITIONS
C     ------------------------------------------------
      DO 40 JR=1,IF
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
   40 J=J+5
C     ----------------------------------------
C     INTERMEDIATE OFF-DIAGONAL 3*3 PARTITIONS
C     ----------------------------------------
      N1=N1+3
      IF (IR.EQ.NNO) GO TO 70
      NB=NNO-IR
      N2=N2+3
      J=J+15
      DO 60 JB=1,NB
      DO 50 JR=1,3
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
   50 J=J+5
      N1=N1+3
      N2=N2+3
   60 J=J+15
      N2=N2-3
   70 I=I+5
   90 LROW=LROW-1
      I=I+15
      N1=N2
      LROW=LROW-3
  100 N2=N2+3*LROW-3
      RETURN
C     --------------------------------------
C     MULTILAYERED ANISOTROPIC COMPOSITE AND
C     ELASTO-PLASTIC RIGIDITY MATRIX
C     --------------------------------------
  120 DO 200 IR=1,NNO
      DO 190 IRR=1,3
      A1=DR(1)*B(I)+DR(2)*B(I+1)+DR(3)*B(I+2)
      A2=DR(2)*B(I)+DR(10)*B(I+1)+DR(11)*B(I+2)
      A3=DR(3)*B(I)+DR(11)*B(I+1)+DR(19)*B(I+2)
      A4=DR(28)*B(I+15)+DR(29)*B(I+16)+DR(30)*B(I+17)
      A5=DR(29)*B(I+15)+DR(37)*B(I+16)+DR(38)*B(I+17)
      A6=DR(30)*B(I+15)+DR(38)*B(I+16)+DR(46)*B(I+17)
      A7=DR(4)*B(I)+DR(12)*B(I+1)+DR(20)*B(I+2)
      A8=DR(5)*B(I)+DR(13)*B(I+1)+DR(21)*B(I+2)
      A9=DR(6)*B(I)+DR(14)*B(I+1)+DR(22)*B(I+2)
      J=I
      IF=4-IRR
C     ------------------------------------------------
C     UPPER TRIANGULAR PART OF DIAGONAL 6*6 PARTITIONS
C     ------------------------------------------------
      DO 140 JR=1,IF
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
  140 J=J+5
      S(N1)=S(N1)+A7*B(J)+A8*B(J+1)+A9*B(J+2)
      S(N1+1)=S(N1+1)+A7*B(J+5)+A8*B(J+6)+A9*B(J+7)
      S(N1+2)=S(N1+2)+A7*B(J+10)+A8*B(J+11)+A9*B(J+12)
C     ---------------------------
C     OFF-DIAGONAL 6*6 PARTITIONS
C     ---------------------------
      N1=N1+3
      IF (IR.EQ.NNO) GO TO 170
      NB=NNO-IR
      N2=N2+3
      J=J+15
      A10=DR( 4)*B(I+15)+DR( 5)*B(I+16)+DR( 6)*B(I+17)
      A11=DR(12)*B(I+15)+DR(13)*B(I+16)+DR(14)*B(I+17)
      A12=DR(20)*B(I+15)+DR(21)*B(I+16)+DR(22)*B(I+17)
      DO 160 JB=1,NB
      DO 150 JR=1,3
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      S(N1+3)=S(N1+3)+A7*B(J+15)+A8*B(J+16)+A9*B(J+17)
      S(N2-3)=S(N2-3)+A10*B(J)+A11*B(J+1)+A12*B(J+2)
      N1=N1+1
      N2=N2+1
  150 J=J+5
      N1=N1+3
      N2=N2+3
  160 J=J+15
      N2=N2-3
  170 I=I+5
  190 LROW=LROW-1
      I=I+15
      N1=N2
      LROW=LROW-3
  200 N2=N2+3*LROW-3
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================

