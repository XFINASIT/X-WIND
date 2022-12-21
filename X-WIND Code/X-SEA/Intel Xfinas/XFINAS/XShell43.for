
C-----------------------------------------------------------------------------
C
	SUBROUTINE COUPLEMEM(PROPM,GT1,GT2,GT3,GT01,GT02,GT03,DM1,DB1) 
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NDF,NWG,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


	DIMENSION PROPM(NMP),DMN(3,3)
	DIMENSION DM1(3,3),DB1(3,3)

	DIMENSION GT1(3),GT2(3),GT3(3)
	DIMENSION GT01(3),GT02(3),GT03(3)

	DMN(1:3,1:3) = 0.0D0

C	MATERIAL PROPERTIES 

      EL  = PROPM(1) !(Young's Modulus)
      PR  = PROPM(2) !(Poission's Ratio)
C
C	COMPUTE FACTOR
C
	A = EL/(1.0-PR*PR)
	B = EL*PR/(1.0-PR*PR)
	C = EL/(2.0*(1.0+PR))
C
	GB11 = GT1(1)*GT1(1)+GT1(2)*GT1(2)+GT1(3)*GT1(3)
	GB12 = GT1(1)*GT2(1)+GT1(2)*GT2(2)+GT1(3)*GT2(3)
	GB22 = GT2(1)*GT2(1)+GT2(2)*GT2(2)+GT2(3)*GT2(3)
	GB33 = GT3(1)*GT3(1)+GT3(2)*GT3(2)+GT3(3)*GT3(3)
C
	DMN(1,1)=B*(GT01(1)*GT01(1)+GT01(2)*GT01(2)+GT01(3)*GT01(3))*GB11+
	1      2.*C*(GT01(1)*GT1(1)+GT01(2)*GT1(2)+GT01(3)*GT1(3))*
	2           (GT01(1)*GT1(1)+GT01(2)*GT1(2)+GT01(3)*GT1(3))
C
	DMN(1,2)=B*(GT01(1)*GT01(1)+GT01(2)*GT01(2)+GT01(3)*GT01(3))*GB22+
	1      2.*C*(GT01(1)*GT2(1)+GT01(2)*GT2(2)+GT01(3)*GT2(3))*
	2           (GT01(1)*GT2(1)+GT01(2)*GT2(2)+GT01(3)*GT2(3))
C
	DMN(1,3)=B*(GT01(1)*GT01(1)+GT01(2)*GT01(2)+GT01(3)*GT01(3))*GB12+
	1      2.*C*(GT01(1)*GT2(1)+GT01(2)*GT2(2)+GT01(3)*GT2(3))*
	2           (GT01(1)*GT1(1)+GT01(2)*GT1(2)+GT01(3)*GT1(3))
C
	DMN(2,1)=B*(GT02(1)*GT02(1)+GT02(2)*GT02(2)+GT02(3)*GT02(3))*GB11+
	1      2.*C*(GT02(1)*GT1(1)+GT02(2)*GT1(2)+GT02(3)*GT1(3))*
	2           (GT02(1)*GT1(1)+GT02(2)*GT1(2)+GT02(3)*GT1(3))
C
	DMN(2,2)=B*(GT02(1)*GT02(1)+GT02(2)*GT02(2)+GT02(3)*GT02(3))*GB22+
	1      2.*C*(GT02(1)*GT2(1)+GT02(2)*GT2(2)+GT02(3)*GT2(3))*
	2           (GT02(1)*GT2(1)+GT02(2)*GT2(2)+GT02(3)*GT2(3))
C
	DMN(2,3)=B*(GT02(1)*GT02(1)+GT02(2)*GT02(2)+GT02(3)*GT02(3))*GB12+
	1      2.*C*(GT02(1)*GT1(1)+GT02(2)*GT1(2)+GT02(3)*GT1(3))*
	2           (GT02(1)*GT2(1)+GT02(2)*GT2(2)+GT02(3)*GT2(3))
C
	DMN(3,1)=B*(GT01(1)*GT02(1)+GT01(2)*GT02(2)+GT01(3)*GT02(3))*GB11+
	1      2.*C*(GT01(1)*GT1(1)+GT01(2)*GT1(2)+GT01(3)*GT1(3))*
	2           (GT02(1)*GT1(1)+GT02(2)*GT1(2)+GT02(3)*GT1(3))
C
	DMN(3,2)=B*(GT01(1)*GT02(1)+GT01(2)*GT02(2)+GT01(3)*GT02(3))*GB22+
	1      2.*C*(GT01(1)*GT2(1)+GT01(2)*GT2(2)+GT01(3)*GT2(3))*
	2           (GT02(1)*GT2(1)+GT02(2)*GT2(2)+GT02(3)*GT2(3))
C
	DMN(3,3)=B*(GT01(1)*GT02(1)+GT01(2)*GT02(2)+GT01(3)*GT02(3))*GB12+
	1         C*(GT01(1)*GT1(1)+GT01(2)*GT1(2)+GT01(3)*GT1(3))*
	2           (GT02(1)*GT2(1)+GT02(2)*GT2(2)+GT02(3)*GT2(3))+
	3         C*(GT01(1)*GT2(1)+GT01(2)*GT2(2)+GT01(3)*GT2(3))*
	4           (GT02(1)*GT1(1)+GT02(2)*GT1(2)+GT02(3)*GT1(3))
C
C	PRE-INTEGRATION COUPLING MEMBRANE
C
	DM1 = 0.0D0
	DM1 =  2.0D0*DMN
	DB1 = (2.0D0/3.0D0)*DMN
C
	RETURN
	END


	SUBROUTINE ENHANCEDCON(PROPM,GT1,GT2,GT3,DM0,DB0) 
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NDF,NWG,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


	DIMENSION PROPM(NMP),DMN(3,3)
	DIMENSION DM0(3,3),DB0(3,3)

	DIMENSION GT1(3),GT2(3),GT3(3)

	DMN(1:3,1:3) = 0.0D0

C	MATERIAL PROPERTIES 

      EL  = PROPM(1) !(Young's Modulus)
      PR  = PROPM(2) !(Poission's Ratio)
C
C	COMPUTE FACTOR
C
	A = EL/(1.0-PR*PR)
	B = EL*PR/(1.0-PR*PR)
	C = EL/(2.0*(1.0+PR))
C
	GB11 = GT1(1)*GT1(1)+GT1(2)*GT1(2)+GT1(3)*GT1(3)
	GB12 = GT1(1)*GT2(1)+GT1(2)*GT2(2)+GT1(3)*GT2(3)
	GB22 = GT2(1)*GT2(1)+GT2(2)*GT2(2)+GT2(3)*GT2(3)
	GB33 = GT3(1)*GT3(1)+GT3(2)*GT3(2)+GT3(3)*GT3(3)
C
C	NATURAL MEMBRANE
C
	DMN(1,1) = A*GB11*GB11
	DMN(1,2) = B*GB11*GB22+2.0*C*GB12*GB12
	DMN(1,3) = A*GB11*GB12
	DMN(2,1) = B*GB11*GB22+2.0*C*GB12*GB12
	DMN(2,2) = A*GB22*GB22
	DMN(2,3) = A*GB12*GB22
	DMN(3,1) = A*GB11*GB12
	DMN(3,2) = A*GB12*GB22 
	DMN(3,3) = B*GB12*GB12+C*(GB11*GB22+GB12*GB12)
C
C	NATURAL PRE-INTEGRATION MEMBRANE
	DM0 = 0.0D0
	DM0 =  2.0D0*DMN
	DB0 = (2.0D0/3.0D0)*DMN
C
	RETURN
	END



C	---------------------------
C	DEFINE THE SURFACE DIRECTOR
C	---------------------------

	SUBROUTINE SURFACEDIRECTOR(PROPG,COORD,DI,TH)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	--------	-----------
C	VARIABLE	DESCRIPTION
C	--------	-----------
C	COORD		ELEMENT COORDINATES
C	P			SHAPE INTERPOLATION FUNCTION DERIVATIVE
C				WITH RESPECT TO xi,eta
C	RP,SP		NODAL POINTS IN NATURAL COORDINATES
C	COVR,COVS	COVARIANT BASES IN xi,eta DIRECTION
C	VECT		VECTOR PRODUCTS 
C	DI			SURFACE DIRECTOR
C	PROPG	    GEOMETRIC PROPERTIES
C				PROPG(1) - NODE
C				PROPG(2) - ELEMENT THICKNESS
C---------------------------------------------
C
	DIMENSION COORD(3,4),P(2,4),H(4)
	DIMENSION DI(3,4),PROPG(2)
C
      DIMENSION COVR(3),COVS(3),VECT(3)
	DIMENSION RP(4),SP(4)
C
	RP(1) =  1.0D0
	RP(2) = -1.0D0
	RP(3) = -1.0D0
	RP(4) =  1.0D0
C
	SP(1) =  1.0D0
	SP(2) =  1.0D0
	SP(3) = -1.0D0
	SP(4) = -1.0D0
C
	TH = PROPG(2)
C
	DO K= 1,4
	H(1:4)     = 0.0D0
	P(1:2,1:4) = 0.0D0
	CALL SHAPE4EAS4N(RP(K),SP(K),H,P,4)
	COVR(1:3) =0.0D0
	COVS(1:3) =0.0D0
      DO I=1,4
      DO J=1,3
      COVR(J)=COVR(J)+P(1,I)*COORD(J,I)
      COVS(J)=COVS(J)+P(2,I)*COORD(J,I)
	END DO
	END DO
C
	VECT(1:3) = 0.0D0
      CALL VECPRD (COVR,COVS,VECT)	
      CALL SCALEN (VECT,VECT,DET,3) ! NORMALIZED VECTOR
C
	DI(1,K) = (TH/2.)*VECT(1)
	DI(2,K) = (TH/2.)*VECT(2)
	DI(3,K) = (TH/2.)*VECT(3)
	END DO
C
	RETURN
	END



	SUBROUTINE COVBASE(COORD,DI,P,H,G1,G2,G3,D1,D2,D)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C
	DIMENSION COORD(3,4),DI(3,4)
	DIMENSION H(4),P(2,4)
	DIMENSION G1(3),G2(3),G3(3),D1(3),D2(3),D(3)
C
C	INITIALIZED VARIABLE
C
	G1(1:3) = 0.0D0
	G2(1:3) = 0.0D0
	G3(1:3) = 0.0D0
C
	D1(1:3) = 0.0D0
	D2(1:3) = 0.0D0
	 D(1:3) = 0.0D0
C
      DO I=1,4
      DO J=1,3
C	COVARIANT BASIS DUE TO MID-SURFACE
      G1(J)=G1(J)+P(1,I)*COORD(J,I)
      G2(J)=G2(J)+P(2,I)*COORD(J,I)   
C	GRADIENTS OF DIRECTORS
	D1(J)=D1(J)+P(1,I)*DI(J,I) 
      D2(J)=D2(J)+P(2,I)*DI(J,I) 
	END DO
	END DO
C	DIRECTOR
      DO I=1,4
      D(1)  = D(1)+H(I)*DI(1,I)
      D(2)  = D(2)+H(I)*DI(2,I)
      D(3)  = D(3)+H(I)*DI(3,I)
	END DO

	G3(1:3) = D(1:3)

	RETURN
	END


      SUBROUTINE SHAPE4EAS4N (XI,ETA,H,P,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     =======================================================
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND DERIVATIVES
C     AT THE NODAL POINTS OF  4 ,ISOPARAMETRIC QUADRILATERAL
C	======================================================

C                         NODE NUMBERING CONVENTION
C                         =========================
C
C                   2                                  1
C
C                     0 . . . . . . . . . . . . . . . 0
C                     .                               .
C                     .                               .
C                     .               S               .
C                     .               .               .
C                     .               .               .
C                     .               . . . R         .
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     0 . . . . . . . . . . . . . . . 0
C
C                   3                                   4
C
C     =====================================================

C	    VARIABLE				DESCRIPTION
C	    ========				===========

C		H(1,4)		4 NODES SHAPE FUNCTION		
C	
C		P(2,NNO)	DERIVATIVE OF SHAPE FUNCTION WITH RESPECT TO xi,eta

C						WITH RESPECT TO xi
C							P(1,1)	= d(H[1])/d_xi
C							P(1,2)	= d(H[2])/d_xi
C							P(1,3)	= d(H[3])/d_xi
C							P(1,4)	= d(H[4])/d_xi

C						WITH RESPECT TO eta
C							P(2,1)	= d(H[1])/d_eta
C							P(2,2)	= d(H[2])/d_eta
C							P(2,3)	= d(H[3])/d_eta
C							P(2,4)	= d(H[4])/d_eta

C     NNO        = NUMBER OF NODES USED TO DESCRIBE ELEMENT

C	=====================================================

      DIMENSION  H(1,4),P(2,4)

C	INITIALIZED VARIABLE

	H(1,1:4)   =  0.0D0
	P(1:2,1:4) =  0.0D0

	IF(NNO.EQ.4)THEN

      RP  = 1.0+XI
      SP  = 1.0+ETA
      RM  = 1.0-XI
      SM  = 1.0-ETA

C     =============================================
C     INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     FOR A FOUR NODE ELEMENT
C     =======================

      H(1,1)   = 0.25*RP*SP
      H(1,2)   = 0.25*RM*SP
      H(1,3)   = 0.25*RM*SM
      H(1,4)   = 0.25*RP*SM

C	DERIVATIVE WITH RESPECT TO - xi

      P(1,1) = 0.25*SP
      P(1,2) = -P(1,1)
      P(1,3) = -0.25*SM
      P(1,4) = -P(1,3)

C	DERIVATIVE WITH RESPECT TO - eta

      P(2,1) = 0.25*RP
      P(2,2) = 0.25*RM
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)

	ELSE
	WRITE(*,*)'!!!THIS ELEMENT NNO IS NOT EQUAL TO 4 CHECK SUB.SHAPE4EAS4N'
	END IF

	RETURN
	END


	SUBROUTINE TRSHEAR13(G1,G2,DI,D,H,P,B)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION G1(3),G2(3)
	DIMENSION P(2,4),H(4),B(1,24)
	DIMENSION DI(3,4),D(3)
C	TRANSLATION DEGREES OF FREDOM
      DO I = 1,4
	DO J = 1,3
	B(1,6*I-(6-J)) = D(J)*P(1,I)
	END DO
	END DO
C	ROTATIONAL DEGREES OF FREEDOM 
	DO I=1,4
	B(1,6*I-(3-1)) = H(I)*(G1(3)*DI(2,I)-G1(2)*DI(3,I))
	B(1,6*I-(3-2)) = H(I)*(G1(1)*DI(3,I)-G1(3)*DI(1,I))
	B(1,6*I-(3-3)) = H(I)*(G1(2)*DI(1,I)-G1(1)*DI(2,I))
	END DO
	RETURN
	END


	SUBROUTINE TRSHEAR23(G1,G2,DI,D,H,P,B)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION G1(3),G2(3)
	DIMENSION P(2,4),H(4),B(1,24)
	DIMENSION DI(3,4),D(3)
C	TRANSLATION DEGREES OF FREDOM
      DO I = 1,4
	DO J = 1,3
	B(1,6*I-(6-J)) = D(J)*P(2,I)
	END DO
	END DO
C	ROTATIONAL DEGREES OF FREEDOM 
	DO I=1,4
	B(1,6*I-(3-1)) = H(I)*(G2(3)*DI(2,I)-G2(2)*DI(3,I))
	B(1,6*I-(3-2)) = H(I)*(G2(1)*DI(3,I)-G2(3)*DI(1,I))
	B(1,6*I-(3-3)) = H(I)*(G2(2)*DI(1,I)-G2(1)*DI(2,I))
	END DO
	RETURN
	END



	SUBROUTINE CTVBASE(G1,G2,G3,D,GT1,GT2,GT3,DET) 
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION G1(3),G2(3),G3(3),D(3)
	DIMENSION GT1(3),GT2(3),GT3(3)

	DIMENSION GB(3,3),GBI(3,3)

	GB(1:3,1:3) = 0.0D0
	GBI(1:3,1:3) = 0.0D0

C	LOAD COVARIANT BASE TO COVARIANT METRIC TENSOR

	GB(1,1:3) = G1(1:3)
	GB(2,1:3) = G2(1:3)
	GB(3,1:3) = G3(1:3)

	CALL INVMAT1(GB,GBI,3)

	GT1 = 0.0D0
	GT2 = 0.0D0
	GT3 = 0.0D0

	GT1(1:3) = GBI(1:3,1)
	GT2(1:3) = GBI(1:3,2)
	GT3(1:3) = GBI(1:3,3)

	DET = GB(1,1)*GB(2,2)*GB(3,3)+GB(1,2)*GB(2,3)*GB(3,1)
	1     +GB(1,3)*GB(2,1)*GB(3,2)-GB(3,1)*GB(2,2)*GB(1,3)
	2     -GB(3,2)*GB(2,3)*GB(1,1)-GB(3,3)*GB(2,1)*GB(1,2)



1500	FORMAT(10F12.5)

	RETURN
	END


	SUBROUTINE STRAIN11(G1,G2,P,B)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION G1(3),G2(3)
	DIMENSION P(2,4),B(1,24)
	B(1,1:24) = 0.0D0
      DO I = 1,4
	DO J = 1,3
C	E11 NORMAL STRAIN
	B(1,6*I-(6-J)) = P(1,I)*G1(J)
	END DO
	END DO
	RETURN
	END



	SUBROUTINE STRAIN22(G1,G2,P,B)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	
	DIMENSION G1(3),G2(3)
	DIMENSION P(2,4),B(1,24)
	B(1,1:24) = 0.0D0
      DO I = 1,4
	DO J = 1,3
C	E22 NORMAL STRAIN
	B(1,6*I-(6-J)) = P(2,I)*G2(J)
	END DO
	END DO
	RETURN
	END


	SUBROUTINE GAMMA12(G1,G2,P,B)
	
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION G1(3),G2(3)
	DIMENSION P(2,4),B(1,24)
	B(1,1:24) = 0.0D0
      DO I = 1,4
	DO J = 1,3
C	E22 NORMAL STRAIN
	B(1,6*I-(6-J)) = P(1,I)*G2(J)+ G1(J)*P(2,I)
	END DO
	END DO
	RETURN
	END



	SUBROUTINE NATURALCONS(PROPM,GT1,GT2,GT3,DM,DB,DS) 
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NDF,NWG,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C
	DIMENSION PROPM(NMP),DMN(3,3),DSN(2,2)
	DIMENSION DM(3,3),DB(3,3),DS(2,2) 
	DIMENSION GT1(3),GT2(3),GT3(3)
C
	DMN(1:3,1:3) = 0.0D0
	DSN(1:2,1:2) = 0.0D0
C
C	MATERIAL PROPERTIES 
C
      EL  = PROPM(1) !(Young's Modulus)
      PR  = PROPM(2) !(Poission's Ratio)
C
C	COMPUTE FACTOR
C
	A = EL/(1.0-PR*PR)
	B = EL*PR/(1.0-PR*PR)
	C = EL/(2.0*(1.0+PR))
C
	GB11 = GT1(1)*GT1(1)+GT1(2)*GT1(2)+GT1(3)*GT1(3)
	GB12 = GT1(1)*GT2(1)+GT1(2)*GT2(2)+GT1(3)*GT2(3)
	GB22 = GT2(1)*GT2(1)+GT2(2)*GT2(2)+GT2(3)*GT2(3)
	GB33 = GT3(1)*GT3(1)+GT3(2)*GT3(2)+GT3(3)*GT3(3)

C	NATURAL MEMBRANE

	DMN(1,1) = A*GB11*GB11
	DMN(1,2) = B*GB11*GB22+2.0*C*GB12*GB12
	DMN(1,3) = A*GB11*GB12
	DMN(2,1) = B*GB11*GB22+2.0*C*GB12*GB12
	DMN(2,2) = A*GB22*GB22
	DMN(2,3) = A*GB12*GB22
	DMN(3,1) = A*GB11*GB12
	DMN(3,2) = A*GB12*GB22 
	DMN(3,3) = B*GB12*GB12+C*(GB11*GB22+GB12*GB12)
C
C	NATURAL TRANSVERSE SHEAR
C
	DSN(1,1) = (5.0/6.0)*C*GB11*GB33
	DSN(1,2) = (5.0/6.0)*C*GB12*GB33
	DSN(2,1) = (5.0/6.0)*C*GB12*GB33
	DSN(2,2) = (5.0/6.0)*C*GB22*GB33
C
C	NATURAL PRE-INTEGRATION MEMBRANE
C
	DM =  2.0D0*DMN
	DB = (2.0D0/3.0D0)*DMN
	DS =  2.0D0*DSN
C
	RETURN
	END




	SUBROUTINE CURVATURE(DI,D1,D2,D,G1,G2,P,BD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION DI(3,4),D1(3),D2(3),D(3)
	DIMENSION G1(3),G2(3)
	DIMENSION P(2,4),BD(3,24)

	BD(1:3,1:24) = 0.0D0

C	CONTRIBUTION OF TRANSLATION DEGREES OF FREDOM
      DO I = 1,4
	DO J = 1,3
	BD(1,6*I-(6-J)) = P(1,I)*D1(J)
	BD(2,6*I-(6-J)) = P(2,I)*D2(J)
	BD(3,6*I-(6-J)) = P(1,I)*D2(J)+P(2,I)*D1(J)
	END DO
	END DO


C	BENDING DUE TO ROTATIONAL DEGREES OF FREEDOM 
	DO I=1,4

C	FIRST ROW

	BD(1,6*I-(3-1)) = P(1,I)*(DI(2,I)*G1(3)-DI(3,I)*G1(2))
	BD(1,6*I-(3-2)) = P(1,I)*(DI(3,I)*G1(1)-DI(1,I)*G1(3))
	BD(1,6*I-(3-3)) = P(1,I)*(DI(1,I)*G1(2)-DI(2,I)*G1(1))

C	SECOND ROW

	BD(2,6*I-(3-1)) = P(2,I)*(DI(2,I)*G2(3)-DI(3,I)*G2(2))
	BD(2,6*I-(3-2)) = P(2,I)*(DI(3,I)*G2(1)-DI(1,I)*G2(3))
	BD(2,6*I-(3-3)) = P(2,I)*(DI(1,I)*G2(2)-DI(2,I)*G2(1))

C	THIRD ROW

	BD(3,6*I-(3-1)) = DI(2,I)*(G2(3)*P(1,I)+G1(3)*P(2,I))
	1           -1.0D0*DI(3,I)*(G2(2)*P(1,I)+G1(2)*P(2,I))

	BD(3,6*I-(3-2)) = DI(3,I)*(G2(1)*P(1,I)+G1(1)*P(2,I))
	1           -1.0D0*DI(1,I)*(G2(3)*P(1,I)+G1(3)*P(2,I))
      
	BD(3,6*I-(3-3)) = DI(1,I)*(G2(2)*P(1,I)+G1(2)*P(2,I))
	1           -1.0D0*DI(2,I)*(G2(1)*P(1,I)+G1(1)*P(2,I))

	END DO

	RETURN
	END



	SUBROUTINE DRILING(XJI,P,H,VR,VS,VT,BD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION XJI(2,2)
	DIMENSION P(2,4),H(4)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION BD(1,24)

      DO 10 I = 1,4

C	LOOP OVER ELEMENT NODAL

      fi = XJI(1,1)*P(1,I)+XJI(1,2)*P(2,I)
      gi = XJI(2,1)*P(1,I)+XJI(2,2)*P(2,I)
     	
C	LOOP OVER STRAIN-DISPLACMENT (VECTOR COMPONENTS)

	DO 20 J = 1,3

C	DRILLING DOF

	BD(1,6*I-(6-J)) = (1.0D0/2.0D0)*(gi*VR(J)-fi*VS(J)) ! u,v,w (COLUMN 1 2 3)
	BD(1,6*I-(3-J)) = H(I)*VT(J)						! theta_x,theta_Y,theta_Z (COLUMN 4 5 6)


20	CONTINUE	
10	CONTINUE

	RETURN
	END 



	SUBROUTINE JACOBSH4(COORD,P,XJ,XJI,DET,VR,VS,VT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION COORD(3,4),P(2,4)
	DIMENSION XJ(2,2),XJI(2,2)
	DIMENSION VR(3),VS(3),VT(3)


      DIMENSION COVR(3),COVS(3),RV(3),SV(3)
C
      CALL CLEARA (COVR,3)
      CALL CLEARA (COVS,3)


      DO 20 I=1,4
      DO 20 J=1,3
      COVR(J)=COVR(J)+P(1,I)*COORD(J,I)
      COVS(J)=COVS(J)+P(2,I)*COORD(J,I)
20	CONTINUE

C*******************************

      CALL VECPRD (COVR,COVS,VT)
      CALL SCALEN (VT,VT,DET,3)

c	complete VT

      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)
      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)


      CALL SCAPRD (COVR,VR,F1,3)
      CALL SCAPRD (COVS,VR,F2,3)
      CALL SCAPRD (COVR,VS,F3,3)
      CALL SCAPRD (COVS,VS,F4,3)


C*********************************

      XJ(1,1) = F1
	XJ(2,1) = F2
	XJ(1,2) = F3
	XJ(2,2) = F4

	XJI(1,1)= F4/DET
      XJI(2,1)=-F2/DET
      XJI(1,2)=-F3/DET
      XJI(2,2)= F1/DET


      RETURN
      END



