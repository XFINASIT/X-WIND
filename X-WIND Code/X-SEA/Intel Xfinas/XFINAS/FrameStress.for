C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE FSTPOIN(PROST,ISSET)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

	DIMENSION PROST(NBP,NBPS),ISSET(1)


C     --------------------------
C     READ STRESS POINT DATA
C     --------------------------
	READ(ITI,*)
	DO IBPS = 1,NBPS
	READ(ITI,*)
	READ
	ENDDO



C     --------------------------
C     ASSIGN TO ELEMENT 
C     --------------------------
	READ(ITI,*)
	DO IE = 1,NELE
	READ (ITI,*) MLE,ISSET(IE)
	ENDDO



	RETURN

	END
C	=================================================================
C	=================================================================
C	=================================================================
