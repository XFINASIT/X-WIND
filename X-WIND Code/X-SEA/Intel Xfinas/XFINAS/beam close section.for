C	=====================================================================
C	=====================================================================
C	=====================================================================

	SUBROUTINE INVERT(NCONS,IERR,D)
C	-------------------------------
C	SUBROUTINE TO INVERT MATRIX
C	-------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
	IMPLICIT INTEGER*4(I-N)

	DIMENSION D(NCONS,NCONS),D1(NCONS,2*NCONS)
C
   	RMAX=0.D0
	DO I=1,NCONS
		D1(I,I+NCONS)=1.D0
		IF (RMAX.LT.DABS(D(I,I)))  RMAX=DABS(D(I,I))
		 DO J=1,NCONS
			D1(I,J)=D(I,J)
			IF (I.NE.J) D1(I,J+NCONS)=0.D0
		 END DO
	END DO
c
	RLIM=1.D-20*RMAX
	DO I=1,NCONS
		IF(DABS(D1(I,I)).LT.RLIM) THEN
			IERR=1
			RETURN
		END IF

		DO J=I+1,NCONS
			RM=D1(I,J)/D1(I,I)
			DO K=J,2*NCONS
				D1(J,K)=D1(J,K)-RM*D1(I,K)
			END DO
		END DO

	END DO
c
	DO I=NCONS,1,-1
		DO K=1,NCONS
			DO J=I+1,NCONS
				D1(I,K+NCONS)=D1(I,K+NCONS)-(D1(I,J)*D1(J,K+NCONS))
			END DO
			D1(I,K+NCONS)=D1(I,K+NCONS)/D1(I,I)
		END DO
	END DO
c
	IERR=0

	DO I=1,NCONS
		DO J=1,NCONS
			D(I,J)=D1(I,J+NCONS)
		END DO
	END DO
C
	RETURN
	END


C=====================================================================
