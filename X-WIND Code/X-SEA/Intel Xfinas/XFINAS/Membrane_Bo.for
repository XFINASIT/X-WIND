
C	==============================================================
C	==============================================================
C	==============================================================
C	==============================================================


	SUBROUTINE MOHRCRI(SIG,EPS,STRAIN,STRESS,DEP,PROPM,WK)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	****************************************************

C	SUBROUTINE COMPUTE SOIL PLASTIC CONSTITUTIVE MATRIX
C		COMPUTE CORRECTION STRESS AND STRAIN

C	****************************************************	

C	SIGMA(6)  = STRESS AT THE END OF THE PREVIOUS UPDATE
C	EPS(6)    = STRAIN AT THE END OF THE PREVIOUS UPDATE
C	EPSTN	  = EPSILON PLASTIC STRAIN
C	STRAIN(6) = TOTAL STRAIN AT CURRENT STAGE
C	STRESS(6) = TOTAL STRESS AT CURRENT STAGE

C	LOCAL VARIABLE

C	SXX = STRESS AT X-DIRECTION
C	SYY = STRESS AT Y-DIRECTION
C	SZZ = STRESS AT Z-DIRECTION
C	TXY = STRESS AT XY-DIRECTION
C	TXZ = STRESS AT XZ-DIRECTION
C	TYZ = STRESS AT YZ-DIRECTION

C	PRES = PREVIOUS STRESS SURFACE F:=(SIG)
C	PREY = PREVIOUS YIELD SURFACE
C	CURS = CURRENT STRESS SURFACE F:=(STRESS)

C	AV   = {a} VECTOR
C	DA   = {d} VECTOR = [D]{a}

C	DELSIG = INCREMENTAL STRESS = {DS} = [DES]{DE}
C	DELEPS = INCREMENTAL STRAIN = {DE} = {E}(i)-{E}(i-1)

C	PSI  = HARDENING AT CALCULATION STAGE
C	EV   =	e-VOID RATIO
C	BLS  = BULK DENSITY OF SOIL
C	GS   = SHEAR MODULUS OF SOIL

C	P	 = MEAN NORMAL EFFECTIVE STRESS

C	Y    = YIELD SURFACE OF MOHR COULOMB
C	ANG  = PRICTION ANGLE INPUT(DEGREE)
C	PHI  = PRICTION ANGLE(RAD)
	
C	RFACT = RATIO FACTOR (ADJUST INSIDE STRESS TO YIELD SURFACE)
C	REDUC = REDUCTION FACTOR OUT SIDE YIELD SURFACE 

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT


      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

	COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST

c	COMMON /MOHR/  COH,ANG,BETA

c	COMMON /CRIP/  OCR,VO,DEN,VGL,SWL,SM,OK

	DIMENSION SIG(6),EPS(6),STRAIN(6),STRESS(6)

	DIMENSION DELSIG(6),DELEPS(6)

	DIMENSION DD(3,3),DV(9),PROPM(*)

	DIMENSION DES(4,4),DPS(4,4),DEP(6,6)

	DIMENSION PRESP(3) ! PRICIPLE STRESS OF PRES
	DIMENSION CURSP(3) ! PRICIPLE STRESS OF CURS


	DIMENSION AV1(4,1),AV2(4,1),AV3(4,1),AV(4,1)
	DIMENSION DA(4,1),CONST(1,1),WK(1)

C	...COMPUTE SOIL PARAMETER

C	WRITE(*,*) 'SONG',PR
	COH  = PROPM(10)
	ANG  = PROPM(11)
	BETA = PROPM(12)

	PI  = 3.141592654 
	PHI = ANG*PI/180.

	HK = YM/(3.*(1.-2.*PR))
	HG = YM/(2.*(1.+PR))

C	**************************************
C	INITIALIZE PLASTIC CONSTITUTIVE MATRIX
C	************************************** 

	DPS   = 0.
	RFACT = 0.
	
	PRESP = 0.
	CURSP = 0.
	EPSTN = 0.

C	CHANGE SIGN FOR POSITIVE COMPRESSION DEFINITION CRITERION

	SIG(1) =  1.0*SIG(1)
	SIG(2) =  1.0*SIG(2)
	SIG(3) =  1.0*SIG(3)
	SIG(4) =  1.0*SIG(4)


C	**************************************************
C	  COMPUTE THE STRESS PREVIOUS YIELD SURFACE (PRES)
C	**************************************************

	P = (1.0/3.0)*(SIG(1)+SIG(2)+SIG(4))

C	...FIRST ENTRY MODULE

	
	IF(P.EQ.0.)THEN
	WK(2) = COH*COS(PHI)
	END IF


	IF(P.EQ.0.0) GOTO 100 !FIRST ENTRY PREVIOUS STRESS = 0.0


C	***************************************************
C	COMPUTE STRESSES INVARIANT AND DEVIATORIC INVARIANT
C	FOR PREVIOUS STRESSES
C	***************************************************

	SX  = SIG(1)
	SY  = SIG(2)
	TXY = SIG(3)
	SZ  = SIG(4)

	TXZ = 0.
	TYZ = 0.

C	...COMPUTE THE STRESS INVARIANT

	PART1 = SX*SY*SZ
	PART2 = 2*TXY*TXZ*TYZ
	PART3 = SX*TYZ*TYZ
	PART4 = SZ*TXY*TXY
	PART5 = SY*TXZ*TXZ
	
	VARI1 = SX+SY+SZ
	VARI2 = SX*SY+SX*SZ+SY*SZ-TXY*TXY-TXZ*TXZ-TYZ*TYZ
	VARI3 = PART1+PART2-PART3-PART4-PART5

C	...COMPUTE DEVIATORIC STRESS INVARIANT

	VARJ2 = (VARI1*VARI1)/3.0-VARI2
	VARJ3 = 2.0*(VARI1*VARI1*VARI1)/27.0-(VARI1*VARI2)/3.0+VARI3

C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE

	SUBV = SQRT(VARJ2*VARJ2*VARJ2)

	IF(SUBV.EQ.0.0)THEN
	WRITE(*,*)' !!! PRINCIPLE STRESS DEVIDED BY ZERO'
	STOP
	END IF

	VALUE = -1.5*SQRT(3.0)*VARJ3/SUBV

	IF(VALUE.GT.1.0) VALUE = 1.0
	IF(VALUE.LT.-1.0) VALUE = -1.0
	
	THETA = (1.0/3.0)*ASIN(VALUE)

C	...COMPUTE PREVIOUS VECTOR YIELD SURFACE

	TERM1 = (1./3.)*VARI1*SIN(PHI)
	TERM2 = SQRT(VARJ2)*COS(THETA)
	TERM3 = SQRT(VARJ2)*(1./SQRT(3.))*SIN(THETA)*SIN(PHI)

	PRES = TERM1+TERM2-TERM3

C	************************************************
C		COMPUTE CURRENT YIELD SURFACE
C	************************************************

C	SETTING HARDENING PARAMETER 

C	WRITE(*,*)SEPSTN


C	...READ PLASTIC STRAIN FROM PREVIOUS SESSION

	PYLD  = WK(2)
	PREY  = PYLD
	YIELD = COH*COS(PHI) + PYLD

	SEPSTN = WK(1)

C	*******************************************
C	COMPUTE CURRENT STRESS YIELD SURFACE VECTOR
C	*******************************************

C	...UPDATING INCREMENTAL STRESSES

	DO I = 1,3
	DELEPS(I) = STRAIN(I) - EPS(I)
	END DO

	DELSIG(1) = A1*DELEPS(1)+B1*DELEPS(2) !DSXX
	DELSIG(2) = B1*DELEPS(1)+A1*DELEPS(2) !DSYY
	DELSIG(3) = C1*DELEPS(3)			  !DTXY
	DELSIG(4) = B1*DELEPS(1)+B1*DELEPS(2) !DSZZ

C	...CHANGE SIGN FOR POSITIVE CRITERION DEFINITION
	
c	DELSIG(1) = -1.0*DELSIG(1)
c	DELSIG(2) = -1.0*DELSIG(2)
c	DELSIG(3) =  1.0*DELSIG(3)
c	DELSIG(4) = -1.0*DELSIG(4)

	DELSIG(1) =  1.0*DELSIG(1)
	DELSIG(2) =  1.0*DELSIG(2)
	DELSIG(3) =  1.0*DELSIG(3)
	DELSIG(4) =  1.0*DELSIG(4)
C	...COMPUTE CURRENT STRESS (TOTAL STRESS)

	DO I =1,4
	STRESS(I) = SIG(I) + DELSIG(I)
	END DO

	SX  = STRESS(1)
	SY  = STRESS(2)
	TXY = STRESS(3)
	SZ  = STRESS(4)

	TXZ = 0.
	TYZ = 0.

C	...COMPUTE THE STRESS INVARIANT

	PART1 = SX*SY*SZ
	PART2 = 2*TXY*TXZ*TYZ
	PART3 = SX*TYZ*TYZ
	PART4 = SZ*TXY*TXY
	PART5 = SY*TXZ*TXZ
	
	VARI1 = SX+SY+SZ
	VARI2 = SX*SY+SX*SZ+SY*SZ-TXY*TXY-TXZ*TXZ-TYZ*TYZ
	VARI3 = PART1+PART2-PART3-PART4-PART5

C	...COMPUTE DEVIATORIC STRESS INVARIANT

	VARJ2 = (VARI1*VARI1)/3.0-VARI2
	VARJ3 = 2.0*(VARI1*VARI1*VARI1)/27.0-(VARI1*VARI2)/3.0+VARI3

C	****************************************
C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE
C	****************************************

	SUBV = SQRT(VARJ2*VARJ2*VARJ2)

	IF(SUBV.EQ.0.0)THEN
	WRITE(*,*)' !!! PRINCIPLE STRESS DEVIDED BY ZERO'
	STOP
	END IF

	VALUE = -1.5*SQRT(3.0)*VARJ3/SUBV

	IF(VALUE.GT.1.0) VALUE = 1.0
	IF(VALUE.LT.-1.0) VALUE = -1.0
	
	THETA = (1.0/3.0)*ASIN(VALUE)

C	...COMPUTE CURRENT VECTOR YIELD SURFACE

	TERM1 = (1./3.)*VARI1*SIN(PHI)
	TERM2 = SQRT(VARJ2)*COS(THETA)
	TERM3 = SQRT(VARJ2)*(1./SQRT(3.))*SIN(THETA)*SIN(PHI)

	CURS = TERM1+TERM2-TERM3



C	**********************************************
C	CHECK THE CURRENT STRESS OUTSIDE YIELD SURFACE
C	**********************************************

C	...PREVIOUS YIELD

	IF(PRES.GE.PREY)THEN
	   IF(CURS.LE.PREY) GOTO 100 ! ELASTIC CASE	
	    RFACT = 0.0
	END IF

	RFACT = 0.

C	...CROSSING FROM YIELD SURFACE

C	IF(PRES.LT.PREY)THEN
C	   IF(CURS.LE.PREY)GOTO 100 ! ELASTIC CASE	
C		CALL RFMOHR(RFACT,SIG,DELSIG,PHI,YIELD)   	   
C	END IF


	IF(PRES.LT.PREY)THEN
		IF(CURS.GE.PREY)THEN 
C		 CALL RFMOHR(RFACT,SIG,DELSIG,PHI,YIELD) !(EXACT INCORRECT)
		 CALL RFMOHR_NEW(RFACT,SIG,DELSIG,PHI,PREY) !SONGSAK NOV2005
		ELSE
	     GOTO 100 
		END IF 	   
	END IF

C	*********************************
C	REDUCTION STRESS TO YIELD SURFACE
C	*********************************

	MSTEP = 30*(CURS/PREY)+1
	ASTEP = 1.0*MSTEP

	REDUC = 1.0-RFACT

C	...COMPUTE INCREMENTAL STRESS CLOSED TO THE YIELD SURFACE

	DO I =1,4
	STRESS(I) = SIG(I)+RFACT*DELSIG(I)
	DELSIG(I) = REDUC*DELSIG(I)/ASTEP
	END DO

C	...COMPUTE SUB INCREMENTAL STRESSES

	DO ISTEP = 1,MSTEP

C	****************************************************
C	...COMPUTE CURRENT STRESS SURFACE AT YIELD CRITERION
C	****************************************************


	DO IDM = 1,4
	  DO JDM = 1,4 
	    DES(IDM,JDM) = 0.0	  
	  END DO
	END DO

C	-----------------
C	...COMPUTE [DES]
C	-----------------

	CALL DMATRIX(DV,DD,1)

	DO I =1,3
	 DO J =1,3
        DES(I,J)=DD(I,J)
	 END DO
	END DO


	DES(1,4)  = B1
	DES(2,4)  = B1

	DES(4,1)  = B1
	DES(4,2)  = B1
      DES(4,4)  = A1

C	***************************
C	    COMPUTE FLOW VECTOR
C	***************************

	SX  = STRESS(1)
	SY  = STRESS(2)
	TXY = STRESS(3)
	SZ  = STRESS(4)

	TXZ = 0.
	TYZ = 0.

C	...COMPUTE THE STRESS INVARIANT

	PART1 = SX*SY*SZ
	PART2 = 2*TXY*TXZ*TYZ
	PART3 = SX*TYZ*TYZ
	PART4 = SZ*TXY*TXY
	PART5 = SY*TXZ*TXZ
	
	VARI1 = SX+SY+SZ
	VARI2 = SX*SY+SX*SZ+SY*SZ-TXY*TXY-TXZ*TXZ-TYZ*TYZ
	VARI3 = PART1+PART2-PART3-PART4-PART5

C	...COMPUTE DEVIATORIC STRESS INVARIANT

	VARJ2 = (VARI1*VARI1)/3.0-VARI2
	VARJ3 = 2.0*(VARI1*VARI1*VARI1)/27.0-(VARI1*VARI2)/3.0+VARI3

C	****************************************
C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE
C	****************************************

	SUBV = SQRT(VARJ2*VARJ2*VARJ2)

	IF(SUBV.EQ.0.0)THEN
	WRITE(*,*)' !!! PRINCIPLE STRESS DEVIDED BY ZERO'
	STOP
	END IF

	VALUE = -1.5*SQRT(3.0)*VARJ3/SUBV

	IF(VALUE.GT.1.0) VALUE = 1.0
	IF(VALUE.LT.-1.0) VALUE = -1.0
	
	THETA = (1.0/3.0)*ASIN(VALUE)

C	FLOW VECTOR COMPONENT {AV1}

	AV1(1,1) = 1.
	AV1(2,1) = 1.
	AV1(3,1) = 0.
	AV1(4,1) = 1.

C	FLOW VECTOR COMPONENT {AV2}

	PM = (1./3.)*(SX+SY+SZ)

	AV2(1,1) = (SX-PM)/(2.*SQRT(VARJ2))
	AV2(2,1) = (SY-PM)/(2.*SQRT(VARJ2))
	AV2(3,1) = (2.*TXY)/(2.*SQRT(VARJ2))
	AV2(4,1) = (SZ-PM)/(2.*SQRT(VARJ2))

C	FLOW VECTOR COMPONENT {AV3}

	AV3(1,1) = (SY-PM)*(SZ-PM)-TYZ*TYZ+VARJ2/3.
	AV3(2,1) = (SX-PM)*(SZ-PM)-TXZ*TXZ+VARJ2/3.
	AV3(3,1) = 2.*(TYZ*TXZ-(SZ-PM)*TXY)
	AV3(4,1) = (SX-PM)*(SY-PM)-TXY*TXY+VARJ2/3.

C	------------------------
C	COMPUTE {a} FLOW VECTOR
C	------------------------

	CF1 = (1./3.)*SIN(PHI)

C	...CONVERT PHI TO DEGREE PHI*180/Pi,REMOVE SINGULAR POINT
C	THETA ANGLE FROM PRICIPLE STRESS

	ABSTHE = ABS(THETA*57.2957795130824)

	IF(ABSTHE.LT.29.0) GOTO 5

C	WRITE(*,*) 'KAKKKK'
		
	CF3 = 0.0
	PLUMI = 3.0
	IF(THETA.GT.0.0) PLUMI = -3.0
	CF2 = 0.5*(SQRT(3.)+PLUMI*CF1*SQRT(3.))

	GOTO 10

5	CONTINUE

	FACF12 = 1.+TAN(THETA)*TAN(3.*THETA)
	FACF22 = CF1*(TAN(3.*THETA)-TAN(THETA))*SQRT(3.)

	CF2 = COS(THETA)*(FACF12+FACF22)

	FACF3= SQRT(3.)*SIN(THETA)+3.*CF1*COS(THETA)
	CF3 =  FACF3/(2.*VARJ2*COS(3.*THETA))

10	CONTINUE

	AV = CF1*AV1+CF2*AV2+CF3*AV3

C	-----------------------------------------------------------
C	... COMPUTE CONSTITUTIVE VECTOR {d} VECTOR DA = D VECTOR A
C	-----------------------------------------------------------

	CONST = MATMUL(TRANSPOSE(AV),MATMUL(DES,AV))

	CALL HARDMOHR(PSI,SEPSTN)

C	PSI = 0.

	SUB = PSI + CONST(1,1)
	DA =  MATMUL(DES,AV)

C	********************************************************
C	...COMPUTE DLAMDA = (1/SUB)*MATMUL(TRANSPOSE(AV),DELSIG)
C		{a}**T{dS}/({a}**T*[DES]*{a})
C	********************************************************

	ATS = 0.0
	DO I=1,4
	ATS = ATS+AV(I,1)*DELSIG(I)
	END DO

	DLAMD = ATS*(1.0/SUB)

	IF(DLAMD.LT.0.0) DLAMD = 0.0

C	...UPDATING TOTAL STRESS RETURN TO CLOSED YIELD SURFACE

	ATSIG = 0.0
	DO I=1,4
	ATSIG = ATSIG+AV(I,1)*STRESS(I)
	STRESS(I) = STRESS(I)+DELSIG(I)-DLAMD*DA(I,1)
	END DO

C	*************************************
C	...BRING IT BACK TO THE YIELD SURFACE
C	*************************************

	SX  = STRESS(1)
	SY  = STRESS(2)
	TXY = STRESS(3)
	SZ  = STRESS(4)

	TXZ = 0.
	TYZ = 0.

C	...COMPUTE THE STRESS INVARIANT

	PART1 = SX*SY*SZ
	PART2 = 2*TXY*TXZ*TYZ
	PART3 = SX*TYZ*TYZ
	PART4 = SZ*TXY*TXY
	PART5 = SY*TXZ*TXZ
	
	VARI1 = SX+SY+SZ
	VARI2 = SX*SY+SX*SZ+SY*SZ-TXY*TXY-TXZ*TXZ-TYZ*TYZ
	VARI3 = PART1+PART2-PART3-PART4-PART5

C	...COMPUTE DEVIATORIC STRESS INVARIANT

	VARJ2 = (VARI1*VARI1)/3.0-VARI2
	VARJ3 = 2.0*(VARI1*VARI1*VARI1)/27.0-(VARI1*VARI2)/3.0+VARI3

C	****************************************
C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE
C	****************************************

	SUBV = SQRT(VARJ2*VARJ2*VARJ2)

	IF(SUBV.EQ.0.0)THEN
	WRITE(*,*)' !!! PRINCIPLE STRESS DEVIDED BY ZERO'
	STOP
	END IF

	VALUE = -1.5*SQRT(3.0)*VARJ3/SUBV

	IF(VALUE.GT.1.0) VALUE = 1.0
	IF(VALUE.LT.-1.0) VALUE = -1.0
	
	THETA = (1.0/3.0)*ASIN(VALUE)

C	...COMPUTE CURRENT VECTOR YIELD SURFACE

	TERM1 = (1./3.)*VARI1*SIN(PHI)
	TERM2 = SQRT(VARJ2)*COS(THETA)
	TERM3 = SQRT(VARJ2)*(1./SQRT(3.))*SIN(THETA)*SIN(PHI)

	CURS = TERM1+TERM2-TERM3

C	...SETTING RETURN STRESS BACK COEFFICIENT FOR EACH INCRE. STEP

	CRETRN = 1.0

	IF(CURS.GT.PREY)THEN
	CRETRN = PREY/CURS
	END IF

	DO I =1,4
	STRESS(I) = CRETRN*STRESS(I)
	END DO

C	************************
C	 COMPUTE PLASTIC STRAIN
C	************************

C	WRITE(720,1000)(STRESS(IS),IS=1,4),PRES,PREY,CURS,RFACT

	EPSTN = EPSTN + DLAMD*ATSIG/CURS

C	************************
	END DO  ! END OF MSTEP
C	************************ 

C	SET THE UPDATING CURRENT YIELD SURFACE AND FLOW VECTOR

	SX  = STRESS(1)
	SY  = STRESS(2)
	TXY = STRESS(3)
	SZ  = STRESS(4)

	TXZ = 0.
	TYZ = 0.

C	...COMPUTE THE STRESS INVARIANT

	PART1 = SX*SY*SZ
	PART2 = 2.0*TXY*TXZ*TYZ
	PART3 = SX*TYZ*TYZ
	PART4 = SZ*TXY*TXY
	PART5 = SY*TXZ*TXZ
	
	VARI1 = SX+SY+SZ
	VARI2 = SX*SY+SX*SZ+SY*SZ-TXY*TXY-TXZ*TXZ-TYZ*TYZ
	VARI3 = PART1+PART2-PART3-PART4-PART5

C	...COMPUTE DEVIATORIC STRESS INVARIANT

	VARJ2 = (VARI1*VARI1)/3.0-VARI2
	VARJ3 = 2.0*(VARI1*VARI1*VARI1)/27.0-(VARI1*VARI2)/3.0+VARI3

C	****************************************
C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE
C	****************************************

	SUBV = SQRT(VARJ2*VARJ2*VARJ2)

	IF(SUBV.EQ.0.0)THEN
	WRITE(*,*)' !!! PRINCIPLE STRESS DEVIDED BY ZERO'
	STOP
	END IF

	VALUE = -1.5*SQRT(3.0)*VARJ3/SUBV

	IF(VALUE.GT.1.0) VALUE = 1.0
	IF(VALUE.LT.-1.0) VALUE = -1.0
	
	THETA = (1.0/3.0)*ASIN(VALUE)

C	FLOW VECTOR COMPONENT {AV1}

	AV1(1,1) = 1.
	AV1(2,1) = 1.
	AV1(3,1) = 0.
	AV1(4,1) = 1.

C	FLOW VECTOR COMPONENT {AV2}

	PM = (1./3.)*(SX+SY+SZ)

	AV2(1,1) = (SX-PM)/(2.*SQRT(VARJ2))
	AV2(2,1) = (SY-PM)/(2.*SQRT(VARJ2))
	AV2(3,1) = (2.*TXY)/(2.*SQRT(VARJ2))
	AV2(4,1) = (SZ-PM)/(2.*SQRT(VARJ2))

C	FLOW VECTOR COMPONENT {AV3}

	AV3(1,1) = (SY-PM)*(SZ-PM)-TYZ*TYZ+VARJ2/3.
	AV3(2,1) = (SX-PM)*(SZ-PM)-TXZ*TXZ+VARJ2/3.
	AV3(3,1) = 2.*(TYZ*TXZ-(SZ-PM)*TXY)
	AV3(4,1) = (SX-PM)*(SY-PM)-TXY*TXY+VARJ2/3.

C	------------------------
C	COMPUTE {a} FLOW VECTOR
C	------------------------

	CF1 = (1./3.)*SIN(PHI)

C	CONVERT PHI TO DEGREE PHI*180/Pi,REMOVE SINGULAR POINT

C	THETA ANGLE FROM PRICIPLE STRESS

	ABSTHE = ABS(THETA*57.2957795130824)

	IF(ABSTHE.LT.29.)GOTO 50
		
	CF3 = 0.0
	PLUMI = 3.0
	IF(THETA.GT.0.0)PLUMI = -3.0
	CF2 = 0.5*(SQRT(3.)+PLUMI*CF1*SQRT(3.))

	GOTO 60

50	CONTINUE

	FACF12 = 1.+TAN(THETA)*TAN(3.*THETA)
	FACF22 = CF1*(TAN(3.*THETA)-TAN(THETA))*SQRT(3.)

	CF2 = COS(THETA)*(FACF12+FACF22)

	FACF3= SQRT(3.)*SIN(THETA)+3.*CF1*COS(THETA)
	CF3 =  FACF3/(2.*VARJ2*COS(3.*THETA))

60	CONTINUE

	AV = CF1*AV1+CF2*AV2+CF3*AV3

C	-----------------------------------------------------------
C	... COMPUTE CONSTITUTIVE VECTOR {d} VECTOR DA = D VECTOR A
C	-----------------------------------------------------------

	CONST = MATMUL(TRANSPOSE(AV),MATMUL(DES,AV))

	CALL HARDMOHR(PSI,SEPSTN)

	SUB = PSI + CONST(1,1)

	DA =  MATMUL(DES,AV)

C	-----------------------------------------
C	... COMPUTE PLASTIC CONSTITUTIVE MATRIX
C	-----------------------------------------

	DPS =  (1./SUB)*(MATMUL(DA,TRANSPOSE(DA)))

C	------------------------------------------
C	COMPUTE ELASTO-PLASTIC CONSTITUTIVE MATRIX
C	------------------------------------------

	DO I = 1,6
	 DO J = 1,6
		DEP(I,J)=DES(I,J) !-DPS(I,J)
	 END DO
	END DO

C	...UPDATING NEW SURFACE

C	WRITE(*,*)EPSTN
C	PAUSE

	PVEPSTN = WK(1)

	SEPSTN = PVEPSTN + EPSTN

	PYLD = COH + PSI*(SEPSTN-PVEPSTN)

C	WRITE(*,*)PYLD

C	IF(ITASK.GT.1) THEN
	WK(1) = SEPSTN
	WK(2) = PYLD
C	ENDIF


C	*******************************
C	STORE STRESS TO PREVIOUS UPDATE
C	*******************************

C	CHANGED DIRECTION BACK

c	SIG(1) = -1.0*STRESS(1)
c	SIG(2) = -1.0*STRESS(2)
c	SIG(3) =  1.0*STRESS(3)
c	SIG(4) = -1.0*STRESS(4)

c	STRESS(1) = -1.0*STRESS(1)
c	STRESS(2) = -1.0*STRESS(2)
c	STRESS(3) =  1.0*STRESS(3)
c	STRESS(4) = -1.0*STRESS(4)

C	IF(ITASK.GT.1) THEN	
	SIG(1) =  1.0*STRESS(1)
	SIG(2) =  1.0*STRESS(2)
	SIG(3) =  1.0*STRESS(3)
	SIG(4) =  1.0*STRESS(4)
C	ENDIF

	STRESS(1) =  1.0*STRESS(1)
	STRESS(2) =  1.0*STRESS(2)
	STRESS(3) =  1.0*STRESS(3)
	STRESS(4) =  1.0*STRESS(4)

C	UPDATING STRAIN

C	IF(ITASK.GT.1) THEN	
	DO I =1,4	
	EPS(I) = STRAIN(I)
	END DO
C	ENDIF

C	WRITE(720,1000)(STRESS(IS),IS=1,4),PRES,PREY,CURS,RFACT
C	WRITE(720,*)''

	GOTO 200

100	CONTINUE

C	*************************
C	ELASTIC SOIL CONSTITUTIVE
C	*************************

	CALL DMATRIX(DV,DD,1)

	DO I =1,3
	 DO J =1,3
        DES(I,J)=DD(I,J)
	 END DO
	END DO

	DES(1,4)  = B1
	DES(2,4)  = B1

	DES(4,1)  = B1
	DES(4,2)  = B1
      DES(4,4)  = A1

	DO I = 1,3
	DELEPS(I) = STRAIN(I) - EPS(I)
	END DO

	DELSIG(1) = A1*DELEPS(1)+B1*DELEPS(2) !DSXX
	DELSIG(2) = B1*DELEPS(1)+A1*DELEPS(2) !DSYY
	DELSIG(3) = C1*DELEPS(3)			  !DTXY
	DELSIG(4) = B1*DELEPS(1)+B1*DELEPS(2) !DSZZ

C	----------------------------------
C	CHANGE SIGN FOR POSITIVE CRITERION
C	----------------------------------
	
c	DELSIG(1) = -1.0*DELSIG(1)
c	DELSIG(2) = -1.0*DELSIG(2)
c	DELSIG(3) =  1.0*DELSIG(3)
c	DELSIG(4) = -1.0*DELSIG(4)
	
	DELSIG(1) =  1.0*DELSIG(1)
	DELSIG(2) =  1.0*DELSIG(2)
	DELSIG(3) =  1.0*DELSIG(3)
	DELSIG(4) =  1.0*DELSIG(4)
C	*************************************
C	COMPUTE CURRENT STRESS (TOTAL STRESS)
C	*************************************

	DO I =1,4
	STRESS(I) = SIG(I) + DELSIG(I)
	END DO

	DO I = 1,6
	 DO J = 1,6
		DEP(I,J)=DES(I,J)-DPS(I,J)
	 END DO
	END DO

C	*************
C	ELASTIC STATE
C	*************

C	CHANGED DIRECTION BACK

45	CONTINUE

c	SIG(1) = -1.0*STRESS(1)
c	SIG(2) = -1.0*STRESS(2)
c	SIG(3) =  1.0*STRESS(3)
c	SIG(4) = -1.0*STRESS(4)

c	STRESS(1) = -1.0*STRESS(1)
c	STRESS(2) = -1.0*STRESS(2)
c	STRESS(3) =  1.0*STRESS(3)
c	STRESS(4) = -1.0*STRESS(4)


	SIG(1) =  1.0*STRESS(1)
	SIG(2) =  1.0*STRESS(2)
	SIG(3) =  1.0*STRESS(3)
	SIG(4) =  1.0*STRESS(4)


	STRESS(1) =  1.0*STRESS(1)
	STRESS(2) =  1.0*STRESS(2)
	STRESS(3) =  1.0*STRESS(3)
	STRESS(4) =  1.0*STRESS(4)

C	UPDATING STRAIN

	DO I =1,4	
	EPS(I) = STRAIN(I)
	END DO


	GOTO 200

200	CONTINUE ! RETURN VALUE


1000	FORMAT(10E14.6)
1001	FORMAT('FIRST STAGE')

	RETURN
	END


C	************************************************************************
C	************************************************************************
C	************************************************************************


	SUBROUTINE RFMOHR(RFACT,SIGI,DELSIGI,PHI,Y)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	*****************************************

C	SUBROUTINE COMPUTE MOHR COULOMB R-FACTOR

C	SIGI    = PREVIOUS STRESS INPUT
C	DELSIGI = DELTA STRESS INPUT
C	PHI	    = PRICTION ANGLE(RAD)
C	Y		= YIELD SURFACE
 
C	*****************************************

	DIMENSION SIGI(4),DELSIGI(4)

	SX  = SIGI(1)
	SY  = SIGI(2)
	TXY = SIGI(3)
	SZ  = SIGI(4)

	DSX  = DELSIGI(1)
	DSY  = DELSIGI(2)
	DTXY = DELSIGI(3)
	DSZ  = DELSIGI(4)

C	--------------------------
C	COMPUTE CONSTANT COMPONENT
C	--------------------------

C	...COMPUTE AF1 COMPONENT

	FAF1 = (DSX*DSX+DSY*DSY)*(SIN(PHI)*SIN(PHI)-1.0)
	FAF2 = 2.0*DSX*DSY*(SIN(PHI)*SIN(PHI)+1.0)
	FAF3 = 4.0*DTXY*DTXY

	AF1=FAF1+FAF2-FAF3

C	...COMPUTE BF1 COMPONENT

	FBF1 = 2.*(SX*DSX+SY*DSY+2*TXY*DTXY)*(SIN(PHI)*SIN(PHI)-1.)
	FBF2 = 2.*(SX*DSY+SY*DSX-2.*TXY*DTXY)*(SIN(PHI)*SIN(PHI)+1.)
	FBF3 = 2*Y*(DSX+DSY)*SIN(PHI)

	BF1 = FBF1+FBF2-FBF3

C	...COMPUTE CF1 COMPONENT

	FCF1 = (SX*SX+SY*SY)*(SIN(PHI)*SIN(PHI)-1.)
	FCF2 = 2.*(SX*SY)*(SIN(PHI)*SIN(PHI)+1.)
	FCF3 = -4.*TXY*TXY-2.*Y*SIN(PHI)*(SX+SY)+Y*Y

	CF1 = FCF1+FCF2+FCF3

	R1 = (-BF1-SQRT(BF1*BF1-4.*AF1*CF1))/(2.*AF1)
	R2 = (-BF1+SQRT(BF1*BF1-4.*AF1*CF1))/(2.*AF1)

C	IF((0.0.LE.R).AND.(R.LE.1.0))THEN

	RFACT = R1

C	ELSE
C
C	WRITE(*,*)'!!! PROBLEM ROOT'
C	WRITE(720,100)R1,R2
C	STOP


C	END IF

C	WRITE(720,*)R
C	WRITE(*,*)AF1,BF1,CF1,R
C	WRITE(720,100)R1,R2
100	FORMAT(4E14.5)
	RETURN
	END

C	COH = 4.
C	ANG = 30.

C	PI  = 3.141592654 
C	PHI = ANG*PI/180.
C	Y   = 2*COH*COS(PHI)

C	SX  =  1.
C	SY  =  -2.
C	TXY =  0.
C	SZ  =  0.

C	DSX  = 6.
C	DSY  = 4.
C	DTXY = 0.
C	DSZ  = 1.5


C	=============================================================
	SUBROUTINE HARDMOHR(HD,TEPSTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)


c	COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
c	COMMON /CRIP/  OCR,VO,DEN,VGL,SWL,SM,OK
c	COMMON /MOHR/  COH,ANG,BETA

	

C	ALPHA = 1000.

C	TERM1 = 1.0-EXP(-ALPHA*TEPSTN)
C	TERM2 = -BETA*TEPSTN

C	HD = YM*(EXP(TERM1*TERM2)) !*(1.-COS(PHI))

c	WRITE(*,*) 'SONG',HD
c100	FORMAT(5E14.5)

	HD = 0.0


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================

C	=============================================================
C=====================================================================
	
	SUBROUTINE HOKLAW12(PROPM,PROPG,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     DETERMINES ALL THE MATERIAL CONSTANTS IN COMMON BLOCK /HOOK/
C	------------------------------------------------------------
C     PROPM(NMP)    = MATERIAL PROPERTIES
C     PROPG(NGP)    = GEOMETRIC PROPERTIES
C     IND           = FLAG FOR PLANE STRAIN(1),OR PLANE STRESS(2)
C     ------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION PROPM(5),PROPG(*)
C
      YM  = PROPM(1)
      PR  = PROPM(2)
      YLD = PROPM(3)
      HP  = PROPM(4)
      DEN = PROPM(5)
C
      D2  = PR/(PR-1.)
      C2  = YM/(1.+PR)
      B2  = PR/(1.-2.*PR)
      A2  = (1.-PR)/(1.-2.*PR)
      BM  = YM/(3.-6.*PR)
      C1  = .5*C2
      D1  = C1/1.2
C
      IF (IND.EQ.2)  GOTO 200
C     ----------------------
C     PLANE STRAIN CONDITION
C     ----------------------
      B1 = B2*C2
      A1 = C2+B1
      GOTO 500
C     ----------------------
C     PLANE STRESS CONDITION
C     ----------------------
 200  A1 = YM/(1.-PR*PR)
      B1 = A1*PR
C     -----------------------
C     NODAL ELEMENT THICKNESS
C     -----------------------
 500  TH  = PROPG(2)
C
      RETURN

C	RETURN A1,B1,C1
C  	
C	RE ARRANGE TO |  A1 B1 0  | AT SUBROUTINE DMATRIX
C		          |  B1 A1 0  |	[DM]
C		          |  0  0  C1 |  

      END

C	=============================================================
C	==============================================================

      SUBROUTINE JACOB2D(MEL,NNO,XY,P,XJ,XJI,DET,MCO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------

      DIMENSION XY(MCO,NNO),P(2,8),XJ(2,2),XJI(4)

C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  I=1,2
        DO 100  J=1,2
         DUM = 0.0
          DO 90   K=1,NNO
 		  DUM = DUM + P(I,K)*XY(J,K)
 90	    CONTINUE
      XJ(I,J) = DUM
 100	CONTINUE

C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = XJ(1,1)*XJ(2,2) - XJ(2,1)*XJ(1,2)
      IF (ABS(DET).LT.1.0E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')
C     -----------------------------
C     INVERSE (XJI) OF THE JACOBIAN
C     -----------------------------

      DUM = 1.0/DET
      XJI(1) =  XJ(2,2)*DUM
      XJI(2) = -XJ(2,1)*DUM
      XJI(3) = -XJ(1,2)*DUM
      XJI(4) =  XJ(1,1)*DUM
C
      RETURN

C	RETURN
C	XJI(1) = TAU11 , XJI(3) = TAU12
C	XJI(2) = TAU21 , XJI(4) = TAU22
 
C	XJ = | J11 J12 |
C		 | J21 J22 |
      END

C=====================================================================

      SUBROUTINE TTOMAT(TTO,XJO)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------------------
C		FOR ENHANCED STRAIN METHOD
C     COMPUTE MATRIX T0 (HERE SET TTO) THE TRANSFORM MATRIX
C     XJO  = JACOBIAN MATRIX AT THE ORIGIN (R,S) = (0,0)
C     TTO  = TRANSFORMATION MATRIX [TTO](0,0) AT CENTER OF ELEMENT.  
C	-------------------------------------------------------------------------

	DIMENSION TTO(3,3),XJO(2,2)

C	...INITIALIZATION [TTO]

      DO 5 I=1,3
	  DO 5 J=1,3
         TTO(I,J)=0.0
5	CONTINUE

C	... FOR TTO(1,1)  TTO(1,2) TTO(2,1) TTO(2,2)

      DO 10 I=1,2
	 DO 10 J=1,2
	      TTO(I,J)=(XJO(J,I))**2
10	CONTINUE

C	...FOR TTO(3,1)  TTO(3,2) 

      TTO(3,1)=XJO(1,1)*XJO(1,2)
	TTO(3,2)=XJO(2,1)*XJO(2,2)
	
C	...FOR TTO(1,3)  TTO(2,3) 

	TTO(1,3)=2*XJO(1,1)*XJO(2,1)
	TTO(2,3)=2*XJO(1,2)*XJO(2,2)
	
C	...FOR TTO(3,3)
      	
      TTO(3,3)=XJO(1,1)*XJO(2,2)+XJO(2,1)*XJO(1,2)

	RETURN

C			        | j11**2    j21**2    2*j11*j21       |
C	[TTO](0,0)	=	| j12**2    j22**2    2*j12*j22       | 
C			        | j11*j12  j21*j22  j11*j22+j21*j12 |(0,0)
C
	END

C=====================================================================
      SUBROUTINE DMATRIX(D,DM,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------
C     GENERATES THE TWO-DIMENSIONAL STRESS-STRAIN MATRIX FOR A
C     LINEAR ELASTIC,ELASTO-PLASTIC ISOTROPIC MATERIAL
C	----------------------------------------------------------
C     DM(3,3)  = ELASTIC,ELASTO-PLASTIC,ISOTROPIC STRESS-STRAIN MATRIX

C		DM	=	| X1 X2 X3 |
C				| Y1 Y2 Y3 | MATRIX POSITION
C				| Z1 Z2 Z3 | 

C	{D} = [DM] REARRAGE TO VECTOR {D} 
C	TO BE INSERTED TO SUBROUTINE KPLNMAT 

C	RETURN D(1..9) = { X1 X2 X3 Y1 Y2 Y3 Z1 Z2 Z3 }
C	RETURN DM(3,3)

C     IND     = FLAG (PL.STRAIN=1,PL.STRESS=2)
C     ----------------------------------------------------------

      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION D(9),DE(3,3),DP(3,3),DM(3,3)
C
      CALL CLEARA (D,9)
C	---------------
C	INITIALIZE [D] 
C	---------------

	DO 10 IDM =1,9
	D(IDM) = 0.0
10	CONTINUE

	DO 15 IDM = 1,3
	  DO 15 JDM = 1,3 
		DM(IDM,JDM) = 0.0
	    DE(IDM,JDM) = 0.0
	    DP(IDM,JDM) = 0.0
15	CONTINUE

C	--------------------------------
C		ELASTICS BEHAVIOUR
C		PLANE STRAIN CONDITION
C     --------------------------------
      DM(1,1)  = A1
      DM(1,2)  = B1
      DM(2,1)  = B1
      DM(2,2)  = A1
      DM(3,3)  = C1

C		DM	=	| A1 B1  0 |
C				| B1 A1  0 | MATRIX PARAMETER
C				|  0  0 C1 | 
      GOTO 300

C	--------------------------------
C		ELASTO-PLASTIC BEHAVIOUR
C		PLANE STRAIN CONDITION
C     --------------------------------

C	ELASTICS BEHAVIOUR
	
      DE(1,1)  = A1
      DE(1,2)  = B1
      DE(2,1)  = B1
      DE(2,2)  = A1
      DE(3,3)  = C1

C		DE	=	| A1 B1  0 |
C				| B1 A1  0 | ELASTIC MATRIX PARAMETER
C				|  0  0 C1 | 

C	COMPUTE CONSTITUTIVE PLASTICS MATRIX [DP]	


C	I WILL WRITE LATER
C	CALL VECTA RETURN VALUE{VA}
C	CALL VECTD RETURN VALUE{VD}
C	CALL DPMATRIX RETURN VALUE[DP]



C	----------------------------------------------
C	SUBTRACT BY CONSTITUTIVE PLASTIC MATRIX [DP]
C	COMPUTE [DM] = [DE]-[DP]
C	----------------------------------------------
	DO 400 IDEP = 1,3
	 DO 400 JDEP = 1,3
		DM(IDEP,JDEP) = DE(IDEP,JDEP)-DP(IDEP,JDEP)
400	CONTINUE

      GOTO 300

C	--------------------------------
C		CONVERT [DM] TO {D}
C	--------------------------------

300	CONTINUE

	KDX = 0
	DO 500 IDX = 1,3
	 DO 500 JDX = 1,3
		KDX = KDX + 1
		D(KDX) = DM(IDX,JDX)
500	CONTINUE


C	WRITE(*,*) 'SONG',DM(1,2),DM(2,1)
      RETURN

C	RETURN
C	[DM](3,3) 
C	D = { X1 X2 X3 Y1 Y2 Y3 Z1 Z2 Z3 } GENERAL(ELASTO PLASTIC OR ELASTIC)
C	EXAMPLE
C	D = { A1 B1 0  B1 A1 0  0  0  C1 } VALUE FOR AN ELASTICS

      END


C=====================================================================
      SUBROUTINE BMATRIX (XY,H,P,XJI,B,BB,ISTYP,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------
C     EVALUATES THE GLOBAL LINEAR STRAIN-DISPLACEMENT MATRIX B
C     FOR A QUADRILATERAL MEMBRANE ELEMENT
C	------------------------------------
C     INPUT VARIABLES
C	---------------
C     XY(2,NNO)  =  COORDINATES FOR ELEMENT NODES
C     H(NNO)     =  SHAPE FUNCTIONS
C     P(2,NNO)   =  SHAPE FUNCTION DERIVATIVES WITH RESP.TO R,S
C     XJI(2,2)   =  INVERSE OF THE JACOBIAN MATRIX
C     B(1,2*NNO) =  STRAIN-DISPLACEMENT MATRIX(ARRAY)
C	BB(3,2*NNO)=  STRAIN-DISPLACEMENT MATRIX
C     ISTYP      =  ELEMENT SUBTYPE
C     NNO        =  NUMBER OF NODES USED TO DESCRIBE ELEMENT
C	SETTING DEFAULT TO B(8) FOR 4 NODES AND 2 DOF PER NODE
C     ---------------------------------------------------------

      DIMENSION XY(2,NNO),H(8),P(2,8),XJI(4),B(2*NNO),BB(3,2*NNO)

C	------------------------------------
C		INITIALIZE [BM](3,8) FOR 4 NODE
C	------------------------------------

	DO 10 IBM =1,3
	  DO 10 JBM =1,2*NNO
		BB(IBM,JBM) = 0.0
10	CONTINUE

C	---------------------------------
C		   VARIABLE DEFINITION
C	---------------------------------
C	FOR 4 NODES

C	P = | H1,R  H2,R  H3,R  H4,R |
C		| H1,S  H2,S  H3,S  H4,S |
C
C	XJI(1) = TAU11 , XJI(3) = TAU12
C	XJI(2) = TAU21 , XJI(4) = TAU22

C	---------------------------------
      
      DO 100  INO=1,NNO
C
C	HIX = TAU11*HI,R + TAU12*HI,S
C	HIY=  TAU21*HI,R + TAU22*HI,S
C	
		HIX = XJI(1)*P(1,INO) + XJI(3)*P(2,INO)
		HIY = XJI(2)*P(1,INO) + XJI(4)*P(2,INO)

C	COMPUTE {B} MATRIX

		B(2*INO-1) = HIX
		B(2*INO-0) = HIY

 100  CONTINUE

C	---------------------------
C	CONVERT {B} TO MATRIX [BM]
C	---------------------------

	DO 200 IBM =1,NNO
	
	BB(1,2*IBM-1) = B(2*IBM-1) !HIX
	BB(2,2*IBM-0) = B(2*IBM-0) !HIY
	BB(3,2*IBM-1) = B(2*IBM-0) !HIY
	BB(3,2*IBM-0) = B(2*IBM-1) !HIX

200	CONTINUE
	
      RETURN

C	RETURN VALUE
C	B = { H1X H1Y H2X H2Y H3X H3Y H4X H4Y }

C	[BB](3,8) = | H1X   0   H2X   0   H3X   0   H4X   0    |
C		        |  0   H1Y   0   H1Y   0   H3Y   0   H4Y   |
C				| H1Y  H1X  H2Y  H2X  H3Y  H3X  H4Y  H4X   |

	END

C	=========================================================================

      SUBROUTINE GEMATRIX (GE,RI,SI,TT,MP,DETJO,DETJ,MEL)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	-------------------------------------------------------------------------
C		FOR ENHANCED STRAIN TERM OF ENHANCED STRAIN METHOD

C		COMPUTE COEFFICIENT MATRIX EM(3,MP)=[TT]*[M]

C		COMPUTE ENHANCED MATRIX [GE] = FAC*[EM(I,J)] 

C		IN MASTER OR LOCAL COORDINATE SYSTEM

C		OUTPUT MATRIX [GE]

C	-------------------------------------------------------------------------

	DIMENSION EM(3,MP),TT(3,3),GE(3,MP),GM(3,MP)
      
	DO 10 I =1,3
		DO 10 J =1,MP
			EM(I,J) = 0.0
			GE(I,J) = 0.0
			GM(I,J) = 0.0 
10	CONTINUE

C	...ENHANCED 2 TERM

	IF (MP.EQ.2) THEN
	  GM(1,1) = RI
	  GM(2,2) = SI   
	  GOTO 100  
	END IF

C	...ENHANCED 4 TERM

	IF (MP.EQ.4) THEN
	  GM(1,1) = RI
	  GM(2,2) = SI
	  GM(3,3) = RI
	  GM(3,4) = SI  
	  GOTO 100    
	END IF

C	...ENHANCED 7 TERM

	IF (MP.EQ.7) THEN
        GM(1,1) = RI
        GM(2,2) = SI
        GM(3,3) = RI
        GM(3,4) = SI
        GM(1,5) = RI*SI
        GM(2,6) = RI*SI
        GM(3,7) = RI*SI
        GOTO 100    
	END IF


100		EM = MATMUL(TT,GM)
          CONST = (DETJO/DETJ) 

	DO 200 IEM = 1,3
		DO 200 JEM = 1,MP 
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
200	CONTINUE


	RETURN
	END

C	=========================================================================

      SUBROUTINE GEMATRIX1(GE,RI,SI,TT,MP,DETJO,DETJ,MEL)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	-------------------------------------------------------------------------
C		FOR ENHANCED STRAIN TERM OF ENHANCED STRAIN METHOD

C		COMPUTE COEFFICIENT MATRIX EM(3,MP)=[TT]*[M]

C		COMPUTE ENHANCED MATRIX [GE] = FAC*[EM(I,J)] 

C		IN MASTER OR LOCAL COORDINATE SYSTEM

C		OUTPUT MATRIX [GE]

C	-------------------------------------------------------------------------

	DIMENSION EM(4,MP),TT(3,3),GE(4,MP),GM(4,MP)
      
	DO 10 I =1,4
		DO 10 J =1,MP
			EM(I,J) = 0.0
			GE(I,J) = 0.0
			GM(I,J) = 0.0 
10	CONTINUE

C	...ENHANCED 2 TERM

	IF (MP.EQ.2) THEN
	  GM(1,1) = RI
	  GM(2,2) = SI   
	  GOTO 100  
	END IF

C	...ENHANCED 4 TERM

	IF (MP.EQ.4) THEN
	  GM(1,1) = RI
	  GM(2,2) = SI
	  GM(3,3) = RI
	  GM(3,4) = SI  
	  GOTO 100    
	END IF

C	...ENHANCED 7 TERM

	IF (MP.EQ.7) THEN
        GM(1,1) = RI
        GM(2,2) = SI
        GM(3,3) = RI
        GM(3,4) = SI
        GM(1,5) = RI*SI
        GM(2,6) = RI*SI
        GM(3,7) = RI*SI
        GOTO 100    
	END IF

100	CONTINUE

	DO 110 I = 1,3
	DO 110 J = 1,MP
	EM(I,J) = 0.0D0
	DO 110 K = 1,3
110	EM(I,J) = EM(I,J) + TT(I,K)*GM(K,J)

	CONST = (DETJO/DETJ) 

	DO 200 IEM = 1,4
		DO 200 JEM = 1,MP 
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
200	CONTINUE


	RETURN
	END

C=====================================================================
C=====================================================================
	SUBROUTINE GEMATRIX2 (GE,RI,SI,TT,MP,DETJO,DETJ,MEL,R12,ETA1,ETA2)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

C	-------------------------------------------------------------------------
C		FOR ENHANCED STRAIN TERM OF ENHANCED STRAIN METHOD

C		COMPUTE COEFFICIENT MATRIX EM(3,MP)=[TT]*[M]

C		COMPUTE ENHANCED MATRIX [GE] = FAC*[EM(I,J)] 

C		IN MASTER OR LOCAL COORDINATE SYSTEM

C		OUTPUT MATRIX [GE]

C	-------------------------------------------------------------------------

	DIMENSION EM(4,MP),TT(3,3),GE(4,MP),GM(4,MP)
      
	DO 10 I =1,4
		DO 10 J =1,MP
			EM(I,J) = 0.0
			GE(I,J) = 0.0
			GM(I,J) = 0.0 
10	CONTINUE



C	...ENHANCED 5 TERM

	IF (MP.EQ.5) THEN
        GM(1,1) = RI - ETA1
        GM(2,2) = SI - ETA2
        GM(3,3) = RI - ETA1
        GM(3,4) = SI - ETA2
        GM(4,5) = RI*SI*DETJ/R12/DETJO 
	END IF

	DO 100 I = 1,3
	DO 100 J = 1,MP
	EM(I,J) = 0.0D0
	DO 100 K = 1,3
100	EM(I,J) = EM(I,J) + TT(I,K)*GM(K,J)
	EM(4,5) = GM(4,5)

	CONST = (DETJO/DETJ) 


	DO 200 IEM = 1,4
		DO 200 JEM = 1,MP 
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
200	CONTINUE


	RETURN
	END

C=====================================================================
C	==============================================================
      
	SUBROUTINE MATEAS(A,B,S,MM,NNO,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     ------------------------------------
C       TO COMPUTE TOTAL STIFFNESS MATRIX
C           K=K(COMPATIBLE)-K(EAS)
C     ------------------------------------

      DIMENSION A(MM,NEF),B(MM,MM),S(1)
	DIMENSION SE(NEF,NEF)

	
	SE = MATMUL(TRANSPOSE(A),MATMUL(B,A))
	K = 0
	DO I=1,NEF
	 DO J=I,NEF
	   K = K+1
	  S(K) =S(K)-SE(I,J)
	 END DO
	END DO	   

      RETURN
      END

C	==============================================================



C	*******************************************************************

	SUBROUTINE RFMOHR_NEW(RFAC,SG,DSG,PHI,YIELD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	*****************************************

C	SUBROUTINE COMPUTE MOHR COULOMB R-FACTOR

C	SIGI    = PREVIOUS STRESS INPUT
C	DELSIGI = DELTA STRESS INPUT
C	PHI	    = PRICTION ANGLE(RAD)
C	Y		= YIELD SURFACE
 
C	*****************************************


	DIMENSION SG(4),DSG(4)
	DIMENSION SGR(4)


	NITER = 20
	RFAC1 = 0.0
	RFAC2 = 1.0
	TOL = 0.00001
	RFACO = RFAC2
	DO ITER = 1,NITER
	RFAC = 0.5*(RFAC1 + RFAC2)
	DO I = 1,4
	SGR(I) = SG(I) + RFAC*DSG(I)
	ENDDO

	SX  = SGR(1)
	SY  = SGR(2)
	TXY = SGR(3)
	SZ  = SGR(4)

	TXZ = 0.
	TYZ = 0.

C	...COMPUTE THE STRESS INVARIANT

	PART1 = SX*SY*SZ
	PART2 = 2*TXY*TXZ*TYZ
	PART3 = SX*TYZ*TYZ
	PART4 = SZ*TXY*TXY
	PART5 = SY*TXZ*TXZ
	
	VARI1 = SX+SY+SZ
	VARI2 = SX*SY+SX*SZ+SY*SZ-TXY*TXY-TXZ*TXZ-TYZ*TYZ
	VARI3 = PART1+PART2-PART3-PART4-PART5

C	...COMPUTE DEVIATORIC STRESS INVARIANT

	VARJ2 = (VARI1*VARI1)/3.0-VARI2
	VARJ3 = 2.0*(VARI1*VARI1*VARI1)/27.0-(VARI1*VARI2)/3.0+VARI3

C	****************************************
C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE
C	****************************************

	SUBV = SQRT(VARJ2*VARJ2*VARJ2)

	IF(SUBV.EQ.0.0)THEN
	WRITE(*,*)' !!! PRINCIPLE STRESS DEVIDED BY ZERO'
	STOP
	END IF

	VALUE = -1.5*SQRT(3.0)*VARJ3/SUBV

	IF(VALUE.GT.1.0) VALUE = 1.0
	IF(VALUE.LT.-1.0) VALUE = -1.0
	
	THETA = (1.0/3.0)*ASIN(VALUE)

C	...COMPUTE CURRENT VECTOR YIELD SURFACE

	TERM1 = (1./3.)*VARI1*SIN(PHI)
	TERM2 = SQRT(VARJ2)*COS(THETA)
	TERM3 = SQRT(VARJ2)*(1./SQRT(3.))*SIN(THETA)*SIN(PHI)

	FR = TERM1+TERM2-TERM3

	TEST = FR - YIELD
	IF(TEST.LE.0.0) RFAC1 = RFAC
	IF(TEST.GT.0.0) RFAC2 = RFAC
	TEST = (RFAC-RFACO)/RFACO
	IF(ABS(TEST).LE.TOL) GO TO 100
	RFACO = RFAC
	ENDDO

100	CONTINUE

C	WRITE(*,*) RFAC

	IF(RFAC.LT.0.00001) RFAC = 0.0


	RETURN
	END

C	*********************************************************
	SUBROUTINE REARDEP(D,DE4)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	REARRANGE RIGIDITY MATRIX
	DIMENSION D(6,6),DE4(6,6),NN(4)

	NN(1:4) = [1,2,3,6]

	D = 0.0D0
	DO I = 1,4
	NI = NN(I)
	DO J = 1,4
	NJ = NN(J)
	D(NI,NJ) = DE4(I,J)
	ENDDO
	ENDDO


	RETURN
	END

C	*********************************************************
