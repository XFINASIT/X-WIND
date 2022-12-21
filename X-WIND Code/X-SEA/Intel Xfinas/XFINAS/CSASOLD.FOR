
C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
	SUBROUTINE SOYSURF(SG6,FN,AVEC6,ALPI1)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	YIELD SURFACE OF DRUCKER-PRAGER FOR SHELL LAYER

	DIMENSION SG6(6),AVEC6(6)
	DIMENSION SG(6),AVECT(6),PST(3)
	DIMENSION AV1(6,1),AV2(6,1),AV3(6,1)

C	======================================================================	
C	PLASTIC POTENTIAL DERIVATIVE MATRIX ( ASSOCIATE FLOW RULE )
C	======================================================================	
	SG(1) = SG6(1)
	SG(2) = SG6(2)
	SG(3) = SG6(3)
	SG(4) = SG6(4)
	SG(5) = SG6(5)
	SG(6) = SG6(6)

	TEST = 0.0D0
	DO I = 1,6
	TEST = TEST + SG(I)*SG(I)
	ENDDO
	TEST = SQRT(TEST)

	IF(TEST.EQ.0.0D0) THEN
	AVEC6(1:6) = 0.0D0
	FN = 0.0D0
	RETURN
	ENDIF


	SX  = SG(1)
	SY  = SG(2)
	SZ  = SG(3)
	TXY = SG(4)
	TXZ = SG(5)
	TYZ = SG(6)

	AV1 = 0.0
	AV2 = 0.0
	AV3 = 0.0
	AVECT  = 0.0

	CALL INVAT(SG,PST,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)

C	FLOW VECTOR COMPONENT {AV1}

	AV1(1,1) = 1.     !SX
	AV1(2,1) = 1.     !SY
	AV1(3,1) = 1.     !SZ
	AV1(4,1) = 0.     !TXY
	AV1(5,1) = 0.     !TXZ
	AV1(6,1) = 0.     !TYZ

C	FLOW VECTOR COMPONENT {AV2}

	PM = (1./3.)*(SX+SY+SZ)

	AV2(1,1) = (SX-PM)/(2.*SQRT(VARJ2))    !SX
	AV2(2,1) = (SY-PM)/(2.*SQRT(VARJ2))    !SY
	AV2(3,1) = (SZ-PM)/(2.*SQRT(VARJ2))    !SZ
	AV2(4,1) = (2.*TXY)/(2.*SQRT(VARJ2))   !TXY
	AV2(5,1) = (2.*TXZ)/(2.*SQRT(VARJ2))   !TXZ
	AV2(6,1) = (2.*TYZ)/(2.*SQRT(VARJ2))   !TYZ


C	FLOW VECTOR COMPONENT {AV3}

	AV3(1,1) = (SY-PM)*(SZ-PM)-TYZ*TYZ+VARJ2/3.
	AV3(2,1) = (SX-PM)*(SZ-PM)-TXZ*TXZ+VARJ2/3.
	AV3(3,1) = (SX-PM)*(SY-PM)-TXY*TXY+VARJ2/3.
	AV3(4,1) = 2.*(TYZ*TXZ-(SZ-PM)*TXY)
	AV3(5,1) = 2.*(TXY*TYZ-(SY-PM)*TXZ)
	AV3(6,1) = 2.*(TXZ*TXY-(SX-PM)*TYZ)

C	======================================================================	
C	END OF PLASTIC POTENTIAL DERIVATIVE MATRIX
C	======================================================================


C	=====================================	
C	YIELD FUNCTION FLOW VECTOR COEFICIENT
C	=====================================

C	=========================
C	VON-MISE
C	=========================
C	FN  = SQRT(VARJ2)
C	CF1 = 0.0
C	CF2 = 1.0
C	CF3 = 0.0
C	=========================
C	END VON-MISE
C	=========================

C	=========================
C	DRUCKER-PRACKER
C	=========================
	FN  = SQRT(3.0)*((ALPI1*VARI1) + SQRT(VARJ2))
	CF1 = ALPI1 
	CF2 = 1.0 
	CF3 = 0.0
C	=========================
C	END DRUCKER-PRACKER
C	=========================


C	============================================	
C	END OF YIELD FUNCTION FLOW VECTOR COEFICIENT
C	============================================


C	=========================	
C	FLOW VECTOR
C	=========================
	DO I = 1,6
	AVECT(I) = 0.0
	AVECT(I) = CF1*AV1(I,1) + CF2*AV2(I,1) + CF3*AV3(I,1)
	ENDDO

	AVEC6(1) = AVECT(1)
	AVEC6(2) = AVECT(2)
	AVEC6(3) = AVECT(3)
	AVEC6(4) = AVECT(4)
	AVEC6(5) = AVECT(5)
	AVEC6(6) = AVECT(6)


	RETURN
	END	



C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
	SUBROUTINE SOMODUC(DMATX,YOUNA,POISM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CONSTRUCT THE CONSTITUTIVE MATRIX ACCORDING TO MATERIAL NUMBER
C	=============================================================
	DIMENSION DMATX(6,6),GASHM(6,6),TRAMX(6,6)

            
	DO 10 I = 1,6
	DO 10 J = 1,6
10	DMATX(I,J) = 0.0
	
C	------------------------------------------------------------------
C	ELASTIC & PLASTIC CONCRETE BEHAVIOR
C	------------------------------------------------------------------	
	FAC = YOUNA / ((1.0 + POISM)*(1.0 - 2*POISM))
	
      DMATX(1,1) = (1.0 - POISM)*FAC
      DMATX(1,2) = POISM*FAC
      DMATX(1,3) = POISM*FAC

      DMATX(2,1) = POISM*FAC
      DMATX(2,2) = (1.0 - POISM)*FAC
      DMATX(2,3) = POISM*FAC

      DMATX(3,1) = POISM*FAC
      DMATX(3,2) = POISM*FAC
      DMATX(3,3) = (1.0 - POISM)*FAC

      DMATX(4,4) = 0.5*FAC*(1.0 - 2*POISM)
      DMATX(5,5) = 0.5*FAC*(1.0 - 2*POISM)
      DMATX(6,6) = 0.5*FAC*(1.0 - 2*POISM)   
    
	GO TO 100

C	------------------------------------------------------------------

100   RETURN
	END


C	=======================================================================
C	=== CONSTRUCTION ANALYSIS =============== PRAMIN OCT 2010 =============
C	=======================================================================
	SUBROUTINE SCMODUC(DMATX,YOUNA,POISM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	CONSTRUCT THE CONSTITUTIVE MATRIX ACCORDING TO MATERIAL NUMBER
C	=============================================================
	DIMENSION DMATX(6,6),GASHM(6,6),TRAMX(6,6)

            
	DO 10 I = 1,6
	DO 10 J = 1,6
10	DMATX(I,J) = 0.0
	
C	------------------------------------------------------------------
C	ELASTIC & PLASTIC CONCRETE BEHAVIOR (PLANE STRESS)
C	------------------------------------------------------------------    
      SMU=YOUNA/(2*(1+POISM))
      SLA=YOUNA*POISM/(1-POISM*POISM)

	DMATX(1,1)=SLA+2*SMU
	DMATX(1,2)=SLA	
      DMATX(2,2)=SLA+2*SMU	
	DMATX(2,1)=SLA	
	
	DMATX(3,3)=YOUNA
	DMATX(4,4)=SMU
	DMATX(5,5)=SMU
	DMATX(6,6)=SMU
    
	GO TO 100

C	------------------------------------------------------------------

100   RETURN
	END


C	======================================================================	
C	======================================================================	
C	======================================================================