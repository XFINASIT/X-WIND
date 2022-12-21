
      SUBROUTINE CHAGEDIS(EDIS,EDIS2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION EDIS(16),EDIS2(12)

	K=0
	DO I=1,16,4
	K=K+1
	EDIS2(K)   = EDIS(I)
	EDIS2(K+1) = EDIS(I+1)
	EDIS2(K+2) = EDIS(I+2)
	K=K+2
	END DO

      RETURN
      END
C	============================================================
C
C	PURPOSE:	TO SOLVE FOR THE GLOBAL STIFFNESS MATRIX T*KL*TT
C
C	============================================================
      SUBROUTINE SOLDSK(ST,TT,SK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	------------------------------------------------------------
C	INPUT VARIABLES
C	ST(300) = ELEMENT STIFFNESS MATRIX IN LOCAL COORDINATES
C	TT(24,24) = TRANSPOSE OF THE TRANSFORMATION MATRIX
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIAVBLES
C	S(*)		= GLOBAL STIFFNESS MATRIX IN ROW FORMAT
C	------------------------------------------------------------
	DIMENSION ST(78),TT(12,12)
	DIMENSION S(78),SK(12,12)

      t3 = TT(1,1)*ST(1)+TT(2,1)*ST(2)
      t7 = TT(1,1)*ST(2)+TT(2,1)*ST(13)
      t18 = TT(1,1)*ST(4)+TT(2,1)*ST(15)
      t22 = TT(1,1)*ST(5)+TT(2,1)*ST(16)
      t33 = TT(1,1)*ST(7)+TT(2,1)*ST(18)
      t37 = TT(1,1)*ST(8)+TT(2,1)*ST(19)
      t48 = TT(1,1)*ST(10)+TT(2,1)*ST(21)
      t52 = TT(1,1)*ST(11)+TT(2,1)*ST(22)
      t75 = TT(1,2)*ST(4)+TT(2,2)*ST(15)
      t79 = TT(1,2)*ST(5)+TT(2,2)*ST(16)
      t90 = TT(1,2)*ST(7)+TT(2,2)*ST(18)
      t94 = TT(1,2)*ST(8)+TT(2,2)*ST(19)
      t105 = TT(1,2)*ST(10)+TT(2,2)*ST(21)
      t109 = TT(1,2)*ST(11)+TT(2,2)*ST(22)
      t138 = TT(1,1)*ST(34)+TT(2,1)*ST(35)
      t142 = TT(1,1)*ST(35)+TT(2,1)*ST(43)
      t153 = TT(1,1)*ST(37)+TT(2,1)*ST(45)
      t157 = TT(1,1)*ST(38)+TT(2,1)*ST(46)
      t168 = TT(1,1)*ST(40)+TT(2,1)*ST(48)
      t172 = TT(1,1)*ST(41)+TT(2,1)*ST(49)
      t195 = TT(1,2)*ST(37)+TT(2,2)*ST(45)
      t199 = TT(1,2)*ST(38)+TT(2,2)*ST(46)
      t210 = TT(1,2)*ST(40)+TT(2,2)*ST(48)
      t214 = TT(1,2)*ST(41)+TT(2,2)*ST(49)
      t237 = TT(1,1)*ST(58)+TT(2,1)*ST(59)
      t241 = TT(1,1)*ST(59)+TT(2,1)*ST(64)
      t252 = TT(1,1)*ST(61)+TT(2,1)*ST(66)
      t256 = TT(1,1)*ST(62)+TT(2,1)*ST(67)
      t279 = TT(1,2)*ST(61)+TT(2,2)*ST(66)
      t283 = TT(1,2)*ST(62)+TT(2,2)*ST(67)
      t300 = TT(1,1)*ST(73)+TT(2,1)*ST(74)
      t304 = TT(1,1)*ST(74)+TT(2,1)*ST(76)
      S(1) = t3*TT(1,1)+t7*TT(2,1)
      S(2) = t3*TT(1,2)+t7*TT(2,2)
      S(3) = TT(1,1)*ST(3)+TT(2,1)*ST(14)
      S(4) = t18*TT(1,1)+t22*TT(2,1)
      S(5) = t18*TT(1,2)+t22*TT(2,2)
      S(6) = TT(1,1)*ST(6)+TT(2,1)*ST(17)
      S(7) = t33*TT(1,1)+t37*TT(2,1)
      S(8) = t33*TT(1,2)+t37*TT(2,2)
      S(9) = TT(1,1)*ST(9)+TT(2,1)*ST(20)
      S(10) = t48*TT(1,1)+t52*TT(2,1)
      S(11) = t48*TT(1,2)+t52*TT(2,2)
      S(12) = TT(1,1)*ST(12)+TT(2,1)*ST(23)
      S(13) = (TT(1,2)*ST(1)+TT(2,2)*ST(2))*TT(1,2)+(TT(1,2)*ST(2)+TT(2,
     #2)*ST(13))*TT(2,2)
      S(14) = TT(1,2)*ST(3)+TT(2,2)*ST(14)
      S(15) = t75*TT(1,1)+t79*TT(2,1)
      S(16) = t75*TT(1,2)+t79*TT(2,2)
      S(17) = TT(1,2)*ST(6)+TT(2,2)*ST(17)
      S(18) = t90*TT(1,1)+t94*TT(2,1)
      S(19) = t90*TT(1,2)+t94*TT(2,2)
      S(20) = TT(1,2)*ST(9)+TT(2,2)*ST(20)
      S(21) = t105*TT(1,1)+t109*TT(2,1)
      S(22) = t105*TT(1,2)+t109*TT(2,2)
      S(23) = TT(1,2)*ST(12)+TT(2,2)*ST(23)
      S(24) = ST(24)
      S(25) = TT(1,1)*ST(25)+TT(2,1)*ST(26)
      S(26) = TT(1,2)*ST(25)+TT(2,2)*ST(26)
      S(27) = ST(27)
      S(28) = TT(1,1)*ST(28)+TT(2,1)*ST(29)
      S(29) = TT(1,2)*ST(28)+TT(2,2)*ST(29)
      S(30) = ST(30)
      S(31) = TT(1,1)*ST(31)+TT(2,1)*ST(32)
      S(32) = TT(1,2)*ST(31)+TT(2,2)*ST(32)
      S(33) = ST(33)
      S(34) = t138*TT(1,1)+t142*TT(2,1)
      S(35) = t138*TT(1,2)+t142*TT(2,2)
      S(36) = TT(1,1)*ST(36)+TT(2,1)*ST(44)
      S(37) = t153*TT(1,1)+t157*TT(2,1)
      S(38) = t153*TT(1,2)+t157*TT(2,2)
      S(39) = TT(1,1)*ST(39)+TT(2,1)*ST(47)
      S(40) = t168*TT(1,1)+t172*TT(2,1)
      S(41) = t168*TT(1,2)+t172*TT(2,2)
      S(42) = TT(1,1)*ST(42)+TT(2,1)*ST(50)
      S(43) = (TT(1,2)*ST(34)+TT(2,2)*ST(35))*TT(1,2)+(TT(1,2)*ST(35)+TT
     #(2,2)*ST(43))*TT(2,2)
      S(44) = TT(1,2)*ST(36)+TT(2,2)*ST(44)
      S(45) = t195*TT(1,1)+t199*TT(2,1)
      S(46) = t195*TT(1,2)+t199*TT(2,2)
      S(47) = TT(1,2)*ST(39)+TT(2,2)*ST(47)
      S(48) = t210*TT(1,1)+t214*TT(2,1)
      S(49) = t210*TT(1,2)+t214*TT(2,2)
      S(50) = TT(1,2)*ST(42)+TT(2,2)*ST(50)
      S(51) = ST(51)
      S(52) = TT(1,1)*ST(52)+TT(2,1)*ST(53)
      S(53) = TT(1,2)*ST(52)+TT(2,2)*ST(53)
      S(54) = ST(54)
      S(55) = TT(1,1)*ST(55)+TT(2,1)*ST(56)
      S(56) = TT(1,2)*ST(55)+TT(2,2)*ST(56)
      S(57) = ST(57)
      S(58) = t237*TT(1,1)+t241*TT(2,1)
      S(59) = t237*TT(1,2)+t241*TT(2,2)
      S(60) = TT(1,1)*ST(60)+TT(2,1)*ST(65)
      S(61) = t252*TT(1,1)+t256*TT(2,1)
      S(62) = t252*TT(1,2)+t256*TT(2,2)
      S(63) = TT(1,1)*ST(63)+TT(2,1)*ST(68)
      S(64) = (TT(1,2)*ST(58)+TT(2,2)*ST(59))*TT(1,2)+(TT(1,2)*ST(59)+TT
     #(2,2)*ST(64))*TT(2,2)
      S(65) = TT(1,2)*ST(60)+TT(2,2)*ST(65)
      S(66) = t279*TT(1,1)+t283*TT(2,1)
      S(67) = t279*TT(1,2)+t283*TT(2,2)
      S(68) = TT(1,2)*ST(63)+TT(2,2)*ST(68)
      S(69) = ST(69)
      S(70) = TT(1,1)*ST(70)+TT(2,1)*ST(71)
      S(71) = TT(1,2)*ST(70)+TT(2,2)*ST(71)
      S(72) = ST(72)
      S(73) = t300*TT(1,1)+t304*TT(2,1)
      S(74) = t300*TT(1,2)+t304*TT(2,2)
      S(75) = TT(1,1)*ST(75)+TT(2,1)*ST(77)
      S(76) = (TT(1,2)*ST(73)+TT(2,2)*ST(74))*TT(1,2)+(TT(1,2)*ST(74)+TT
     #(2,2)*ST(76))*TT(2,2)
      S(77) = TT(1,2)*ST(75)+TT(2,2)*ST(77)
      S(78) = ST(78)

C	FORMULATE SK
C	----------------------
	KK=0
	DO I=1,12
	DO J=I,12
	KK=KK+1
	SK(I,J)	= S(KK)
	SK(J,I) = SK(I,J)
	END DO
	END DO

      RETURN
      END
C	================================================================
C     ----------------------------------------------------------------
C     CALCULATE COUPLING PART OT THE LINEAR STIFFNESS MATRIX
C	-----------------------------------------------------------
C     CC(12,4)   = COUPLING TERM
C     B(2,3*NNO) = LINEAR STRAIN-DISPLACEMENT MATRIX
C     FAC        = DET*WT*TH
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C	H(4)	   = SHAPE FUNCTION
C     ----------------------------------------------------------------
      SUBROUTINE COUPL4 (CC,H,B,FAC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION B(2,12),CC(12,4),H(4)

	DO J = 1,4
	DO I = 1,12
	CC(I,J) = CC(I,J) + (B(1,I)+B(2,I))*H(J)*FAC
	END DO
	END DO

	RETURN
      END
C
C	================================================================
C	ASSEMBLY OF THE ALL THE TERMS CALCULATED ABOVE IN LEFT HAND SIDE
C	USE LINER INTERPOLATION IN TIME USING FINITE DIFFERENCE AND
C	USING CRANK-NICOLSON TYPE APPROXIMATION THETA=1/2 IN BOTH EQUATIONS
C	AND TAKE CONSTANT TIME INCREMENT Dt
C	===============================================================
      SUBROUTINE ASSEMB4 (SS,SM,CC,SP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /SEEP/ NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	DIMENSION SS(16,16),SM(12,12),CC(12,4),SP(4,4)

C	ADD THE CONTRIBUTION FROM SOLID PART TO THE STIFNESS MATRIX
C	-----------------------------------------------------------
	DO I = 1,12
	DO J = 1,12
	SS(I,J) = SM(I,J)
	END DO
	END DO

C	ADD THE COUPLING TERM TO THE STIFNESS MATRIX
C	--------------------------------------------
	DO I = 1,12
	DO J = 1,4
	SS(I,J+12) = CC(I,J)
	SS(J+12,I) = CC(I,J)
	END DO
	END DO

C	ADD THE CONTRIBUTION FROM POREWATER PART TO THE STIFNESS MATRIX
C	----------------------------------------------------------------
	DO I = 1,4
	DO J = 1,4
	SS(I+12,J+12) = -DTINC*SP(I,J)
	END DO
	END DO

	RETURN
	END
C	============================================================
      SUBROUTINE ARRENG4 (S,SS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATE COUPLING PART OT THE LINEAR STIFFNESS MATRIX
C	-----------------------------------------------------------
C     CC(16,8)   = COUPLING TERM
C     B(3,2*NNO) = LINEAR STRAIN-DISPLACEMENT MATRIX
C     FAC        = DET*WT*TH
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C	H(8)	   = SHAPE FUNCTION
C     ----------------------------------------------------------------
      DIMENSION S(136),SS(16,16),IAR4(16,16),SS1(16,16)

	CALL CLEAR (16,16,IAR4)
	CALL CLEAR (16,16,SS1)

C	4	 NODE ELEMENT
C	----------------------------------------
	IAR4(1,1)   = 1
	IAR4(2,2)   = 1
	IAR4(3,3)   = 1
	IAR4(4,5)   = 1
	IAR4(5,6)   = 1
	IAR4(6,7)   = 1
	IAR4(7,9)   = 1
	IAR4(8,10)  = 1
	IAR4(9,11)  = 1
	IAR4(10,13) = 1
	IAR4(11,14) = 1
	IAR4(12,15) = 1
	IAR4(13,4)  = 1
	IAR4(14,8)  = 1
	IAR4(15,12) = 1
	IAR4(16,16) = 1
	
	SS1 = SS1 + MATMUL(MATMUL(TRANSPOSE(IAR4),SS),IAR4)

C	MAKE AND ARRANGE A UPPER TRANGULAR MATRIX ROW WISE
C	--------------------------------------------------
	II = 0
	DO I = 1,16
	DO J = I,16
	II = II + 1
	S(II) = S(II) + SS1(I,J)
	END DO
	END DO

	RETURN
      END
C
C	================================================================
C	FIND THE RESULTANT FORCE VECTOR AND ADD THE
C	CONTRIBUTION FROM THE Nth TIME STEP
C	ASSEMBLY OF THE ALL THE TERMS CALCULATED ABOVE IN RIGHT HAND SIDE
C	USE LINER INTERPOLATION IN TIME USING FINITE DIFFERENCE AND
C	USING CRANK-NICOLSON TYPE APPROXIMATION THETA=1/2 IN BOTH EQUATIONS
C	AND TAKE CONSTANT TIME INCREMENT Dt
C	----------------------------------------------------------------
      SUBROUTINE RESNTH (SS,SM,CC,SP,RE,EDIS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /SEEP/ NTSTEP,KTSTEP,CTIME,DTINC,KFRES

	DIMENSION SS1(16,16),SM(12,12),CC(12,4),SP(4,4),SS2(16,16)
	DIMENSION RE(16),EDIS(16),SS(16,16),SS3(16,16)

C	INTIALIZE INTERNAL VARIABLES
C	----------------------------
	CALL CLEAR (16,16,SS1)

C	ASSEMBLE COUPLING TERM
C	----------------------
	DO I = 1,12
	DO J = 1,4
	SS1(J+12,I) = CC(I,J)
	END DO
	END DO

C	ARRANGE THE SS MATRIX ACCORDING TO THE ORDER OF 
C	DISPLACEMENT VECTOR(U1,V1,Uw1,......Un,Vn,Uwn...) n=1,NNO
C	---------------------------------------------------------
	CALL ARRE2 (SS2,SS1)

C	CALCULATE THE RESULTANT FORCE VECTORE
C	-------------------------------------	
	re = -matmul(ss2,edis)

C	CALCULATE THE RESULTANT FORCE VECTORE FROM SPECIFIED PORE PRESSURE
C	------------------------------------------------------------------
C	CALCULATE THE ADITIONAL FORCE DUE TO PRE SPECIFIED PORE PRESUURE
C	------------------------------------------------------------------
	CALL ARRE2 (SS3,SS)
		

	RETURN
	END
C
C	================================================================
      SUBROUTINE ARRE2 (SS1,SS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATE COUPLING PART OT THE LINEAR STIFFNESS MATRIX
C	-----------------------------------------------------------
C     CC(12,4)   = COUPLING TERM
C     ----------------------------------------------------------------
      DIMENSION SS(16,16),IAR4(16,16),SS1(16,16)

	CALL CLEAR (16,16,IAR4)
	CALL CLEAR (16,16,SS1)

C	4	 NODE ELEMENT
C	--------------------------------------
	IAR4(1,1)   = 1
	IAR4(2,2)   = 1
	IAR4(3,3)   = 1
	IAR4(4,5)   = 1
	IAR4(5,6)   = 1
	IAR4(6,7)   = 1
	IAR4(7,9)   = 1
	IAR4(8,10)  = 1
	IAR4(9,11)  = 1
	IAR4(10,13) = 1
	IAR4(11,14) = 1
	IAR4(12,15) = 1
	IAR4(13,4)  = 1
	IAR4(14,8)  = 1
	IAR4(15,12) = 1
	IAR4(16,16) = 1
	
	SS1 = SS1 + MATMUL(MATMUL(TRANSPOSE(IAR4),SS),IAR4)
	
	RETURN
      END
C
