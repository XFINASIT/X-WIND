C	=============================================================
C	=============================================================
      SUBROUTINE SHJACO2 (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,ID,f,JDIRTH)
	IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
C
C     --------------------------------------------------------------
C     EVALUATES LOCAL DIRECTION COSINE VECTORS, JACOBIAN
C     DETERMINANT, JACOBIAN INVERSE AND BASE VECTOR PARAMETERS
C
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     COORD(3,NNO)      = CURRENT NODAL COORDINATES X,Y,Z
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     XJI(4)            = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     DET               = JACOBIAN DETERMINANT
C     RR,SS             = SQUARED BASE VECTOR LENGTHS
C     SNA               = SIN OF ANGLE SUBTENDED BY BASE VECTORS
C     ID                = FLAG
C	f				  = jacobian coefficient stored in coulmn	
C     COVR(3),COVS(3)   = COVARIENT BASE VECTORS ALONG RN,SN
C     --------------------------------------------------------------
C
      DIMENSION COORD(3,*),HD(2,8),VR(3),VS(3),VT(3),XJI(4)
      DIMENSION COVR(3),COVS(3),RV(3),SV(3),f(4)
C
      CALL CLEARA (COVR,3)
      CALL CLEARA (COVS,3)
      DO 20 I=1,NNO
      DO 20 J=1,3
      COVR(J)=COVR(J)+HD(1,I)*COORD(J,I)
   20 COVS(J)=COVS(J)+HD(2,I)*COORD(J,I)
   
      DO II = 1,3
        CR = COVR(II) 
        CS = COVS(II)
        IF(CR.EQ.0.0.AND.CS.EQ.0.0) JDIRTH = II
      ENDDO
      
      CALL VECPRD (COVR,COVS,VT)
      CALL SCALEN (VT,VT,DET,3)
      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)
C	NEXT LINE - CHANGED BY GILSON - JULY2002
C      CALL ADDVEC (SV,VS,VS,3)
      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)
      IF (ID.EQ.2) RETURN
      RR=RL*RL
      SS=SL*SL
      SNA=DET/(RL*SL)
      CALL SCAPRD (COVR,VR,F1,3)
      CALL SCAPRD (COVS,VR,F2,3)
      CALL SCAPRD (COVR,VS,F3,3)
      CALL SCAPRD (COVS,VS,F4,3)
      f(1)=f1
	f(2)=f2
	f(3)=f3
	f(4)=f4
	XJI(1)= F4/DET
      XJI(2)=-F2/DET
      XJI(3)=-F3/DET
      XJI(4)= F1/DET
      SNA=F4/SL
      RETURN
C  f() added by Hari for jacobian------
      END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SHKLINEAS2 (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD,FACT,JSTIFTRS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION S(*),DR(64),B(240),BDRL(48)
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
      
      DIMENSION JSTIFTRS(3) ! FOR SHEAR WALL BY BJ
      DIMENSION JSTIFTRS2(12) ! FOR SHEAR WALL BY BJ
C
      JSTIFTRS2(1:3) = JSTIFTRS(1:3)
      
      NN = 1
      DO NNODE = 1,4
        JSTIFTRS2(NN) = JSTIFTRS(1) + 6*(NNODE-1)
        JSTIFTRS2(NN+1) = JSTIFTRS(2) + 6*(NNODE-1)
        JSTIFTRS2(NN+2) = JSTIFTRS(3) + 6*(NNODE-1)
        NN = NN + 3
      ENDDO
      
      N=1
      I=4
      FDRL=FACT*DR(19)
C     -----------------------------
C     TRANSVERSE SHEAR CONTRIBUTION
C     -----------------------------
      DO 30 NROW=1,NEF
      JROW = 0
      A1=DR(55)*B(I) + DR(56)*B(I+1)
      A2=DR(63)*B(I) + DR(64)*B(I+1)
      A3=FDRL*BDRL(NROW)
      J=I
        DO JJ = 1,12
          JST = JSTIFTRS2(JJ)
          IF(NROW.EQ.JST) JROW = 1
        ENDDO
        
      DO 20 NCOL=NROW,NEF
      JCOL = 0
      
        DO KK = 1,12
          KST = JSTIFTRS2(KK)
          IF(NCOL.EQ.KST) JCOL = 1
        ENDDO
      
        IF(JROW.EQ.1.OR.JCOL.EQ.1) THEN
            S(N) = S(N) + 0.0
        ELSE
            S(N)=S(N)+A1*B(J)+A2*B(J+1)+A3*BDRL(NCOL)
        ENDIF
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
      SUBROUTINE WALLSTIFF (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD,FACT,JSTIFTRS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
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

C     SAVE IEL
      COMMON /INELE/ IEL     
      
      DIMENSION S(*),DR(64),B(240),BDRL(48)      
      DIMENSION JSTIFTRS(3) ! FOR SHEAR WALL BY BJ
      DIMENSION JSTIFTRS2(12) ! FOR SHEAR WALL BY BJ      
      
      
c      JSTIFTRS2(1:3) = JSTIFTRS(1:3)
      
C      JSTIFTRS(1) = 3  !w
C      JSTIFTRS(2) = 4  !theta x
C      JSTIFTRS(3) = 5  !theta y

      NN = 1
      DO NNODE = 1,4
        JSTIFTRS2(NN) = JSTIFTRS(1) + 6*(NNODE-1)   !w for each node
        JSTIFTRS2(NN+1) = JSTIFTRS(2) + 6*(NNODE-1) !theta x for each node
        JSTIFTRS2(NN+2) = JSTIFTRS(3) + 6*(NNODE-1) !theta y for each node
        NN = NN + 3
      ENDDO
      
      N=1
      I=4
C     -----------------------------
C     TRANSVERSE SHEAR CONTRIBUTION
C     -----------------------------
      DO 30 NROW=1,NEF
        JROW = 0
        J=I
        DO JJ = 1,12
          JST = JSTIFTRS2(JJ)
          IF(NROW.EQ.JST) JROW = 1
        ENDDO
        
        DO 20 NCOL=NROW,NEF
            JCOL = 0
            IF(JROW.EQ.1) GO TO 10
            DO KK = 1,12
                KST = JSTIFTRS2(KK)
                IF(NCOL.EQ.KST) JCOL = 1
            ENDDO
      
10        IF(JROW.EQ.1.OR.JCOL.EQ.1) THEN
            IF(NCOL.EQ.NROW) THEN
                S(N) = 1.0
            ELSE
                S(N) = 0.0
            ENDIF
        ENDIF
      N=N+1
   20 J=J+5
   30 I=I+5      
      
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SHDELA2 (DR,DVOL,SLA,SLB)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATES ELASTIC RIGIDITIES FOR MEMBRANE-PLATE OR SHELL ELEM.
C	---------------------------------------------------------------
C     DR(8,8) = ELASTIC MEMBRANE,BENDING AND SHEAR RIGIDITIES
C     DVOL    = WT*DET
C     SLA,SLB = SHEAR LOCKING CONSTRAINT FACTORS
C	----------------------------------------------------------------
C	PLEASE NOT SLA,SLB ARE 1.0 MODIFIED ASSUMED STRAIN METHOD, 
C	NO FACTOR REQUIRED
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION DR(64)
C
      CALL CLEARA (DR,64)

C	WRITE(*,*) PR
C	PAUSE

      FACT = 1.0


      TH3  = TH*TH*TH/12.
      FACM = DVOL*TH
CBJ      FACB = DVOL*TH3
      FACB = DVOL*TH3
C
      DR(1)  = FACM*A1 !MEMBRANE
      DR(10) = FACM*A1 
      DR(2)  = FACM*B1 
      DR(9)  = FACM*B1 
      DR(19) = FACM*C1*FACT
      DR(28) = FACB*A1*FACT !BENDING
      DR(37) = FACB*A1*FACT
      DR(29) = FACB*B1*FACT
      DR(36) = FACB*B1*FACT
      DR(46) = FACB*C1
      DR(55) = FACM*SLA*D1 !SHEAR
      DR(64) = FACM*SLB*D1

C
      RETURN
      END
C     =============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SHBMAT11 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,RN,SN,BA,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------
C	B matrix containing transverse shear terms SHBMAT1, Hari
c     ----------------------------------------------------------------
c      EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
c	-------------------------------------------
c     NNO              = NUMBER OF NODES FOR ELEMENT
C     H(4)             = SHAPE FUNCTIONS
C     HD(2,4)          = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(4),HS(4)      = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3)= LOCAL DIRECTION COSINE VECTORS
C     XJI(4)           = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(120)           = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     BA(4,120)        = B VALUES AT SAMPLI NG POINTS, ASSUME STRAIN
C     SN1,SN2,SN3,SN4  = strain interpolation function in natural 
C						coordinate system
C    -----------------------------------------------------------------
      DIMENSION H(4),HD(2,4),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION B(120),BA(4,120)

C	---------------------------------------------------------------
C	COMPUTE INTERPOLATION VALUES AT GAUSS POINTS, 4 SAMPLING POINTS
C	---------------------------------------------------------------
	SN1 = 0.5*(1+RN)
	SN2 = 0.5*(1-RN)
	SN3 = 0.5*(1+SN)
	SN4 = 0.5*(1-SN)
	M = 4
	DO 20 J=1,4
	DO 20 I=1,6
	B(M)=0.0
	B(M+1)=0.0
	B(M)  = BA(3,M)*SN3  +BA(4,M)*SN4
 	B(M+1)= BA(1,M+1)*SN1+BA(2,M+1)*SN2
	M=M+5
  20	CONTINUE
	
C	--------------------------------------------------------------
C	NOW TRANSFER BS FROM NATURAL COORDINATE SYSTEM TO LOCAL SYSTEM
C	NOTE ALL BA(120) VALUES ARE NOT USED
C	------------------------------------
      M=4
	DO 30 I=1,NEF
	B(M)   = XJI(1)*B(M) +XJI(3)*B(M+1)
CBJ	B(M+1) = XJI(2)*B(M) +XJI(4)*B(M+1) ! ROTATION TERM OF SHEAR
  	M=M+5
  30	CONTINUE
C	----------------------------------
C	OTHER TERMS SEE ANOTHER SUBROUTINE
C	----------------------------------
      RETURN
      END 
C     =============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SHBMAT22 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
C	------------------------------------------
C     NNO              = NUMBER OF NODES FOR ELEMENT
C     H(8)             = SHAPE FUNCTIONS
C     HD(2,8)          = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)      = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3)= LOCAL DIRECTION COSINE VECTORS
C     XJI(4)           = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)           = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     ----------------------------------------------------------------
      DIMENSION H(8),HD(2,8),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION B(240)

      M=1
      N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      F1=A2*VR(J)
      F2=A1*VS(J)
      B(M)=A1*VR(J)
      B(M+1)=A2*VS(J)
CBJ      B(M+2)=F1+F2

CBJ      B(N)=F2
CBJ      B(N+1)=-F1
      B(N+2)=B(M+1)-B(M)

      M=M+5
      N=N+5
   10 CONTINUE
      M=M+15
      N=N+15
   20 CONTINUE
      RETURN
      END
C
C     =============================================================
C	=============================================================
C	=============================================================