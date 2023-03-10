C     ======================================================================
C     ======================================================================
C     ======================================================================
	SUBROUTINE COMPINP(PROPM,ISET,IN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C***********************************************************************
C	PURPOSE
C		TO READ IN ALL THE DATA
C	USAGE
C		CALLED IN PROGRAM A8913A
C	DESCRIPTION OF PARAMETERS
C		SEE INPUT DATA SECTION OF DATA ITEM NO. 89013
C		LSMAT - NUMBER OF LAYER TYPES
C		PNL	  - ARRAY CONTAINING LAYER PROPERTIES: LIMITED TO 50
C				STORES THE FOLLOWING, IN ORDER. FOR EACH LAYER K:
C				PNL(K,1)=C,PNL(K,2) =EA ,PNL(K,3) =EB ,PNL(K,4) =EZ
C						   PNL(K,5) =VAB,PNL(K,6) =VAZ,PNL(K,7) =VBZ
C						   PNL(K,8) =GBZ,PNL(K,9) =GAZ,PNL(K,10)=GAB
C						   PNL(K,11)=TL, PNL(K,12)=PHI
C				NOTE THAT ORDER OF INPUT IS
C						 C			 EA				EB			 EZ
C									 VAB			VBZ		     VAZ
C									 GAB			GBZ			 GAZ
C				INPUT ORDER IS SET INTO ARRAY PNL TO ACCOUNT FOR THE
C				DIFFERENCE IN ORDER OF PROPERTIES BETWEEN ABOVE SETS
C		LSSK  - NUMBER OF LAYERS IN STACKING SEQUENCE
C		LAYSS - ARRAY STORING LAY-UP CODE FOR LAMINATES
C				SIGN IS SIGN OF PHI
C		NOTE10 =1 IF MORE THAN 50 MATERIALS IN FILE
C		NOTE11 =1 IF MORE THAN 500 LAYERS SPECIFIED
C		NOTE13 =1 IF STACK CALLS UNDECLARED LAYER
C	SUBROUTINES AND FUNCTIONS CALLED
C	NONE
C	METHOD
C		SEE DATA INPUT SECTION OF DATA ITEM N0.89013
C
C.......................................................................

      DIMENSION PROPM(1)
      
	ALLOCATABLE PNL(:,:),LAYS(:),TLIST(:),ALIST(:),ZT(:),DTMP(:)


C	READ IN NUMBER OF SORTS OF LAYER AND CHECK SIZE LIMIT
	READ(IN,*)
	READ(IN,*) ISET,LSMAT,IDV,IDE,IFC


C	READ INTO ARRAY PNL LAYER TYPE NUMBER AND PROPERTIES
      ALLOCATE(PNL(20,LSMAT))
      ISAND=0
	DO 30 J=1,LSMAT
	
		READ(IN,*) PNL(1 ,J),PNL(2 ,J),PNL(3 ,J),PNL(4 ,J)
		READ(IN,*) PNL(5 ,J),PNL(7 ,J),PNL(6 ,J)
		READ(IN,*) PNL(10,J),PNL(8 ,J),PNL(9 ,J)
        READ(IN,*) PNL(13,J),PNL(14,J),PNL(16,J),PNL(17,J),PNL(18,J),PNL(19,J),PNL(20,J)
        READ(IN,*) PNL(11,J) !MASS DENSITY
      
C	CHECK IF THE SANDWICH ANALYSIS SHOULD CONSIDER CONSTANT SHEAR STRESSES
C	ISAND= 1 FOR CONSTANT SHEAR STRESS ANALYSIS, OTHERWISE 0
	  SAND=PNL(2,J)+PNL(3,J)+PNL(4,J)+PNL(10,J)
	  IF (SAND.EQ.0.0) ISAND=1
		    
30	CONTINUE


C	READ IN STACK LENGTH, THAT IS NUMBER OF LAYERS, AND CHECK SIZE
C	NEXT LINE ADDED BY GILSON - JUN2004
	READ(IN,*)
	READ(IN,*) LSSK

C	READ IN STACK
      ALLOCATE(LAYS(LSSK),TLIST(LSSK),ALIST(LSSK))
	READ(IN,*)
	DO ISSK = 1,LSSK
	  READ(IN,*) LAYS(ISSK),TLIST(ISSK),ALIST(ISSK) !Layer MAT NO, THICKNESS, ANGLE (deg)
      ENDDO

      ALIST(1:LSSK) = ALIST(1:LSSK)*0.017453292519943 !CHANGE ANGLE TO RADIAN

C	0 - NOT HELD FLAT, 1 - HELD FLAT
	READ(IN,*)
	READ(IN,*) IBEND


C	READ AND COMPUTE FOR THERMAL LOADING
      ALLOCATE(DTMP(LSSK))

C	CURING TEMPERATURE
	READ(IN,*)
	READ(IN,*) TCUR
	
C	THERMAL DISTRIBUTION INDICATOR
	READ(IN,*)
56	READ(IN,*) ITEMP
C
	IF(ITEMP.EQ.1) GOTO 66
C	LINEAR VARIATION

	IF(ITEMP.EQ.2) GOTO 58
C	VARIES IN EACH LAYER

C	NO THERMAL LOADING (ITEMP=0) - CHANGE IN TEMPERATURE DTMP=0.0
	DO 93 J1=1,LSSK
	  DTMP(J1) = 0.0D0
93	CONTINUE
	GOTO 74


C	READ THE SERVICE TEMPERATURE IN EACH LAYER
58	READ(IN,*)
	DO 62 J1=1,LSSK
	  READ(IN,*) DTMP(J1)
62	CONTINUE
	GOTO 74


C	READ SERVICE TEMPERATURE FOR TOP, MIDDLE AND BOTTOM FOR LINEAR VARIATION
66	READ(IN,*) ST1,ST2,ST3

C	COMPUTE FOR TOTAL THICKNESS ( S )
74	S=0.0
	DO 150 I=1,LSSK
		DK=TLIST(I)
		S=S+DK
150	CONTINUE

C	COMPUTE FOR MID HEIGHT (DMID)
	DMID=S/2.0

C	COMPUTE FOR LAMINATE DISTANCE FROM MIDPLANE - UPWARD POSITIVE  (BOTTOM TO TOP)
      ALLOCATE(ZT(LSSK))
	S1=0.0
	DO 160 I=1,LSSK
	  J=LSSK-I+1
	  DK=TLIST(J)
	  ZT(J)=S1+DK/2.0-DMID
	  S1=S1+DK
160	CONTINUE


	
	IF (ITEMP.EQ.2) GOTO 84
C	VARIES IN EACH LAYER

	IF (ITEMP.NE.1) GOTO 70
C	NO THERMAL LOADING

C	COMPUTE FOR LAYER TEMPERATURE IN LINEAR VARIATION DISTRIBUTION
	DMID1=-DMID
	ZM1=(ST1-ST2)/DMID
	ZM2=(ST3-ST2)/DMID1
	DO 170 J1=1,LSSK
		IF(ZT(J1).GE.0.0) DTMP(J1)=ZM1*ZT(J1)+ST2
		IF(ZT(J1).LT.0.0) DTMP(J1)=ZM2*ZT(J1)+ST2
170	CONTINUE


C	COMPUTE FOR THE CHANGE IN TEMPERATURE
84	DO 97 ICUR=1,LSSK
		DTMP(ICUR)=TCUR-DTMP(ICUR)
97 	CONTINUE


C	END OF THERMAL LOADING
70	CONTINUE



C	COMPUTE EQUIVALENT MASS DENSITY 
	CDEN = 0.0
	DO 300 K = 1,LSSK
	  KS  = LAYS(K)
	  CDEN = CDEN + PNL(11,KS)*(TLIST(K)/S)
300	CONTINUE

C     --------------------------------------------------------
      NUM0 = 0
      
      NUM1 = NUM0+1 ; NUM2 = NUM+10         !FOR GENERAL DATA
            
      NUM3 = NUM2 +1 ; NUM4 = NUM2+20*LSMAT  !FOR PNL

      NUM5 = NUM4+1 ; NUM6 = NUM4+LSSK      !FOR LAYSS
      
      NUM7 = NUM6 +1 ; NUM8 = NUM6+LSSK      !FOR TLIST
      
      NUM9 = NUM8 +1 ; NUM10= NUM8+LSSK      !FOR ALIST

      NUM11= NUM10+1 ; NUM12= NUM10+LSSK     !FOR ZT
      
      NUM13= NUM12+1 ; NUM14= NUM12+LSSK     !FOR DTMP
      
      PROPM(NUM1:NUM10) = 0.0D0
      
      PROPM(1) = FLOAT(LSSK)    !NUMBER OF LAYER IN STACKS
      PROPM(2) = FLOAT(IBEND)   !HELD FLAT
      PROPM(3) = FLOAT(ITEMP)   !TEMPERATURE CHECK
      PROPM(4) = FLOAT(IFC)     !FAILURE CRITERIA
      PROPM(5) = CDEN           !EQUIVALENT MASS DENSITY
      PROPM(6) = FLOAT(LSMAT)   !NUMBER OF LAYER MATERIALS
      PROPM(7) = FLOAT(ISAND)   !ISAND= 1 FOR CONSTANT SHEAR STRESS ANALYSIS, OTHERWISE 0
      
      NN = NUM4-NUM3
      CALL MOVE(PNL,PROPM(NUM3),NN)
      
      PROPM(NUM5:NUM6)   = FLOAT(LAYS(1:LSSK))  !LAYER MATERIAL LIST
      
      PROPM(NUM7:NUM8)   = TLIST(1:LSSK)        !THICKNESS LIST
      
      PROPM(NUM9:NUM10)  = ALIST(1:LSSK)        !ANGLE LIST
      
      PROPM(NUM11:NUM12) = ZT(1:LSSK)           !LAYER LEVEL LIST
      
      PROPM(NUM13:NUM14) = DTMP(1:LSSK)         !LATER DELTA TEMP. LIST
C     --------------------------------------------------------
      
      
      DEALLOCATE(PNL,LAYS,TLIST,ALIST,ZT,DTMP)



	RETURN
	END
C
C     ======================================================================
C     ======================================================================
C     ======================================================================
	SUBROUTINE COMRGD(EPS,STS,DR,PROPM,TH,PANG,DVOL,ETYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*5 ETYP

C
C	TITLE
C		TRANSVERSE (THROUGH-THE-THICKNESS) SHEAR STIFFNESSES OF FIBRE
C		REINFORCED COMPOSITE LAMINATED PLATES
C
C	PURPOSE
C		TO ESTIMATE THE THROUGH-THE-THICKNESS SHEAR STIFFNESSES OF A
C		LAMINATED COMPOSITE FROM THE PROPERTIES OF THE LAMINATES
C		AND THEIR ORIENTATIONS.
C
C	DESCRIPTION OF MAJOR PARAMETERS IN CONTROLLING PROGRAM
C		THROUGHOUT ALL FLAGS ARE OF THE FORM NOTEXX WHERE XX IS AN
C		IDENTIFYING NUMBER WHICH IS INITILISED AT ZERO AND RAISED
C		TO ONE TO TRIP THE NOTE
C			KOUNT	DUMMY VALUE SET AT UNITY, USED IN CALL OF C8713A
C			LAYS	ARRAY THE FIRST ROW OF WHICH IS SET EQUAL TO ARRAY
C				LAYSS SO THAT C8713A CAN BE USED
C		EQUATION NUMBERS GIVEN RELATE ESDU DOCUMENTATION OF THE
C		SOLUTION
C
      DIMENSION EPS(1),STS(1),PROPM(1),DR(1)
      DIMENSION AA(3,3),BB(3,3),DD(3,3),GG(2,2),DSHELL(8,8),DSOLID(9,9)
      
      ALLOCATABLE PNL(:,:),LAYS(:),TLIST(:),ALIST(:),ZT(:),DTMP(:)
      

C     --------------------------------------------------------
      LSSK  = INT(PROPM(1)) !NUMBER OF LAYER IN STACKS
      IBEND = INT(PROPM(2)) !HELD FLAT
      ITEMP = INT(PROPM(3)) !TEMPERATURE CHECK
      IFC   = INT(PROPM(4)) !FAILURE CRITERIA
      CDEN  = PROPM(5)      !EQUIVALENT MASS DENSITY
      LSMAT = INT(PROPM(6)) !NUMBER OF LAYER MATERIALS
      ISAND = INT(PROPM(7)) !ISAND= 1 FOR CONSTANT SHEAR STRESS ANALYSIS, OTHERWISE 0
C     --------------------------------------------------------
      NUM0 = 0
      
      NUM1 = NUM0+1 ; NUM2 = NUM+10         !FOR GENERAL DATA
            
      NUM3 = NUM2 +1 ; NUM4 = NUM2+20*LSMAT  !FOR PNL

      NUM5 = NUM4+1 ; NUM6 = NUM4+LSSK      !FOR LAYSS
      
      NUM7 = NUM6 +1 ; NUM8 = NUM6+LSSK      !FOR TLIST
      
      NUM9 = NUM8 +1 ; NUM10= NUM8+LSSK      !FOR ALIST

      NUM11= NUM10+1 ; NUM12= NUM10+LSSK     !FOR ZT
      
      NUM13= NUM12+1 ; NUM14= NUM12+LSSK     !FOR DTMP

      ALLOCATE(PNL(20,LSMAT),LAYS(LSSK),TLIST(LSSK),ALIST(LSSK),ZT(LSSK),DTMP(LSSK))
      
      NN = NUM4-NUM3
      CALL MOVE(PROPM(NUM3),PNL,NN)
      
      LAYS(1:LSSK) = INT(PROPM(NUM5:NUM6))  !LAYER MATERIAL LIST
 
      TLIST(1:LSSK)= PROPM(NUM7:NUM8)       !THICKNESS LIST
      
      ALIST(1:LSSK)= PROPM(NUM9:NUM10)      !ANGLE LIST
      
      ZT(1:LSSK)   = PROPM(NUM11:NUM12)     !LAYER LEVEL LIST
      
      DTMP(1:LSSK) = PROPM(NUM13:NUM14)     !LATER DELTA TEMP. LIST
C     --------------------------------------------------------
      ALIST(1:LSSK)= ALIST(1:LSSK) + PANG   !ADD THE ANGLE OF REFERENCE VECTOR
 
	CALL CALCA(EPS,IFC,PNL,LSSK,LAYS,TLIST,ALIST,ZT,TH,ISAND,AA,BB,DD,GG,EE)

		
	SELECTCASE(ETYP)
	
	CASE('SHELL')

      DSHELL(1:8,1:8) = 0.0D0
      DO I = 1,3
          DO J = 1,3
              DSHELL(I+0,J+0) = AA(I,J)
              DSHELL(I+3,J+3) = DD(I,J)
              
              DSHELL(I+0,J+3) = BB(J,I)
              DSHELL(I+3,J+0) = BB(J,I)
          ENDDO
      ENDDO
      DO I = 1,2
          DO J = 1,2
              DSHELL(I+6,J+6) = GG(I,J)
          ENDDO
      ENDDO
      
      STS(1) = AA(1,1)*EPS(1)+AA(1,2)*EPS(2)+AA(1,3)*EPS(3)
     1       + BB(1,1)*EPS(4)+BB(1,2)*EPS(5)+BB(1,3)*EPS(6)
      STS(2) = AA(2,1)*EPS(1)+AA(2,2)*EPS(2)+AA(2,3)*EPS(3)
     1       + BB(2,1)*EPS(4)+BB(2,2)*EPS(5)+BB(2,3)*EPS(6)
      STS(3) = AA(3,1)*EPS(1)+AA(3,2)*EPS(2)+AA(3,3)*EPS(3)
     1       + BB(3,1)*EPS(4)+BB(3,2)*EPS(5)+BB(3,3)*EPS(6)
      STS(4) = BB(1,1)*EPS(1)+BB(1,2)*EPS(2)+BB(1,3)*EPS(3)
     1       + DD(1,1)*EPS(4)+DD(1,2)*EPS(5)+DD(1,3)*EPS(6)
      STS(5) = BB(2,1)*EPS(1)+BB(2,2)*EPS(2)+BB(2,3)*EPS(3)
     1       + DD(2,1)*EPS(4)+DD(2,2)*EPS(5)+DD(2,3)*EPS(6)
      STS(6) = BB(3,1)*EPS(1)+BB(3,2)*EPS(2)+BB(3,3)*EPS(3)
     1       + DD(3,1)*EPS(4)+DD(3,2)*EPS(5)+DD(3,3)*EPS(6)        
      STS(7) = GG(1,1)*EPS(7)+GG(1,2)*EPS(8)
      STS(8) = GG(1,1)*EPS(7)+GG(2,2)*EPS(8)
      
	DR(1:64)=0.0D0
	K = 0
      DO I = 1,8
          DO J = 1,8
              K = K + 1
              DR(K) = DSHELL(J,I)*DVOL
          ENDDO
      ENDDO
	
	CASE('SOLID')

      STS(1) = AA(1,1)*EPS(1)+AA(1,2)*EPS(2)+AA(1,3)*EPS(3)
     1       + BB(1,1)*EPS(4)+BB(1,2)*EPS(5)+BB(1,3)*EPS(6)
      STS(2) = AA(2,1)*EPS(1)+AA(2,2)*EPS(2)+AA(2,3)*EPS(3)
     1       + BB(2,1)*EPS(4)+BB(2,2)*EPS(5)+BB(2,3)*EPS(6)
      STS(3) = AA(3,1)*EPS(1)+AA(3,2)*EPS(2)+AA(3,3)*EPS(3)
     1       + BB(3,1)*EPS(4)+BB(3,2)*EPS(5)+BB(3,3)*EPS(6)
      STS(4) = BB(1,1)*EPS(1)+BB(1,2)*EPS(2)+BB(1,3)*EPS(3)
     1       + DD(1,1)*EPS(4)+DD(1,2)*EPS(5)+DD(1,3)*EPS(6)
      STS(5) = BB(2,1)*EPS(1)+BB(2,2)*EPS(2)+BB(2,3)*EPS(3)
     1       + DD(2,1)*EPS(4)+DD(2,2)*EPS(5)+DD(2,3)*EPS(6)
      STS(6) = BB(3,1)*EPS(1)+BB(3,2)*EPS(2)+BB(3,3)*EPS(3)
     1       + DD(3,1)*EPS(4)+DD(3,2)*EPS(5)+DD(3,3)*EPS(6) 
      STS(7) = EE*EPS(7)     
      STS(8) = GG(1,1)*EPS(8)+GG(1,2)*EPS(9)
      STS(9) = GG(1,1)*EPS(8)+GG(2,2)*EPS(9)

      DSOLID(1:9,1:9) = 0.0D0
      DO I = 1,3
          DO J = 1,3
              DSOLID(I+0,J+0) = AA(I,J)
              DSOLID(I+3,J+3) = DD(I,J)
              
              DSOLID(I+0,J+3) = BB(J,I)
              DSOLID(I+3,J+0) = BB(J,I)
          ENDDO
      ENDDO
      DO I = 1,2
          DO J = 1,2
              DSOLID(I+7,J+7) = GG(I,J)
          ENDDO
      ENDDO
      DSOLID(7,7) = EE
            	
	DR(1:81)=0.0D0
	K = 0
      DO I = 1,9
          DO J = 1,9
              K = K + 1
              DR(K) = DSOLID(J,I)
          ENDDO
      ENDDO
      	
	ENDSELECT
	
	
      DEALLOCATE(PNL,LAYS,TLIST,ALIST,ZT,DTMP)

	RETURN
	END

C     ======================================================================
C     ======================================================================
C     ======================================================================
	SUBROUTINE CALCA(EPS,IFC,PNL,LSSK,LAYS,TLIST,ALIST,ZT,TH,ISAND,AA,BB,DD,GG,EE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
C SUBROUTINE CALCA
C
C PURPOSE
C	TO CALCULATE THE STIFFNESS MATRICES AND BUILT IN STRESSES
C	OF THE LAMINATE
C
C USAGE
C	CALL CALCA(REM1,REM2,REM6,LSSK,IOUT)
C
C DESCRIPTION OF PARAMETERS
C	AA,BB,DD - IN-PLANE, BENDING-STRECHING COUPLING AND BENDING
C					SIFFNESS MATRICES
C	BK, BBK - TRANSFORMED STIFFNESS MATRICES
C	CK		- LAYER STIFFNESS MATRIX
C	CON     - COMPLIANCE MATRIX ABD**-1
C	DK		- LAYER THICKNESS
C	EK#		- ELASTIC MODULI
C	G12		- IN-PLANE SHEAR MODULI
C	LAND#	- COEFFICIENT OF THERMAL EXPANSION
C	LSSK	- N0. OF LAYERS IN THE LAMINATE
C	MTP,NTP - TEMPERATURE PLANE STRESS RESULTANTS
C	QBAR	- LAYER TRANSFORMED REDUCED STIFFNESS MATRIX
C	REM#	- MATRIX THAT RECORDS LAYER FAILURE
C	THETAK	- LAYER FIRBER ANGLE
C	TK,TTK	- TRANSFORMATION MATRICES
C	TOR#	- POISSON'S RATIO
C	ZTSAVE	- LAYER DISTANCE TO THE MID-PLANE
C
C SUBROUTINES AND FUNCTION PROGRAMS REQUIRED
C	UNSMDT
C	UNSMIN
C	MULT1
C

	DIMENSION EPS(6),REM(3,LSSK),LAYS(1),TLIST(1),ALIST(1),PNL(20,1),ZT(1)
	DIMENSION CK(3,3),QK(2,2),BK(3,3)
	DIMENSION CT(3,3),QT(2,2),STN(3),STS(3)
	
	DIMENSION AA(3,3),BB(3,3),DD(3,3),GG(2,2)
	
      AA(1:3,1:3)=0.0
      BB(1:3,1:3)=0.0
      DD(1:3,1:3)=0.0
      EE = 0.0D0
      
C     ----------------------
	DO 1000 II=1,LSSK
C     ----------------------
C
		KS=LAYS(II)
C		
		PHI = ALIST(II)
		ETH = PNL(4 ,KS)
		CALL QIJ(CK,QK,CT,QT,PNL(1,KS),PHI,'LOCL')
C
        XT = PNL(16,KS)
        XC = PNL(17,KS)
        YT = PNL(18,KS)
        YC = PNL(19,KS)
         Q = PNL(20,KS)
C		
		ZK = ZT(II)*TH
	  DO I = 1,3
	      STN(I) = EPS(I) + ZK*EPS(I+3)
	  ENDDO
C
        STS = MATMUL(MATMUL(CK,CT),STN)
        CALL COMPFAIL(STS,REM(1,II),XT,XC,YT,YC,Q,IFC)	
        
C		CHECK FOR COMPLETE PLY FAILURE
		PF = REM(1,II)+REM(2,II)+REM(3,II)
		IF(PF.EQ.3.0) GOTO 1000
C
		DT  = TLIST(II)*TH
C
		CALL QIJ(BK,QK,CT,QT,PNL(1,KS),PHI,'GOBL')
C
		TB1 = DT
		TB2 = DT*ZK
		TB3 = DT*ZK*ZK + DT*DT*DT/12.0
C
		DO 30 I=1,3
		DO 30 J=1,3
			AA(I,J) = AA(I,J) + BK(I,J)*TB1
			BB(I,J) = BB(I,J) - BK(I,J)*TB2
			DD(I,J) = DD(I,J) + BK(I,J)*TB3
30		CONTINUE

C     THICK/EE = sum of --> Dt/Et
        EE = EE + DT/ETH
        
C     ----------------------
1000	CONTINUE
C     ----------------------

C     EE = THICK/ (sum of --> Dt/Et)
      EE = TH/EE !THICKNESS STIFFNESS

      CALL COMPSHEAR(REM,PNL,LSSK,LAYS,TLIST,ALIST,ZT,TH,ISAND,GG)     
C
	RETURN
	END

C     ======================================================================
C     ======================================================================
C     ======================================================================
	SUBROUTINE COMPFAIL(STS,REM,XT,XC,YT,YC,Q,IFC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C SUBROUTINE TOPIT
C
C PURPOSE
C	TO CALCULATE FAILURE 
C
C USAGE
C	CALL TOPIT(FAC,R1,R2)
C
C DESCRIPTION OF PARAMETERS
C	F	- FAILURE CRITERION
C	FAC	- MATRIX OF CALCULATED LAYER STRESSES
C	H#, IHMAX - FAILURE MODE CRITERION
C	R1	- FIBRE FAILURE FLAG
C	R2	- LAYER MATRIX FAILURE FLAG (Y DIRECTION)
C	R6  - IN-PLANE SHEAR FAILURE
C	X0,Y0,Z0 - LAYER STRESSES
C	XC,XT,YC,YT,Q,XTER,YTER - MATERIAL STRENGTHS
C
C SUBROUTINES AND FUNCTION PROGRAMS REQUIRED
C	NONE
C
C METHOD
C	THE METHOD EMPLOYED IS MODIFIED PUCK CRITERION
C

	DIMENSION STS(3),REM(3)
C
C	CHANGE STRESS VARIABLES
	X0=STS(1)
	Y0=STS(2)
	Z0=STS(3)
C
C	DETERMINE APPRPRIATE ALLOWABLE STRESS, i.e. COMPRESSION OR TENSION
	IF(X0.GE.0.0) XTER= XT
	IF(X0.LT.0.0) XTER=-XC
	IF(Y0.GE.0.0) YTER= YT
	IF(Y0.LT.0.0) YTER=-YC
C
C	INTIALIZE FLAGS
	REM(1:3)=0.0
	F=0.0
C
C	CHECK FOR FAILURE USING CHOSEN CRITERION
	SELECT CASE (IFC)
C
C	  TSAI-WU CRITERION
	  CASE(1)
		F=(1.0/XT-1.0/XC)*X0+(1.0/YT-1.0/YC)*Y0
	    F=F+((X0)**2)/(XT*XC)+((Y0)**2)/(YT*YC)
	    F=F+(Z0/Q)**2
		IF (F.GT.(1.0)) THEN
C
	    H1=(1.0/XT-1.0/XC)*X0+((X0)**2)/(XT*XC)
		  H2=(1.0/YT-1.0/YC)*Y0+((Y0)**2)/(YT*YC)
		  H6=(Z0/Q)**2
C
C		  IDENTIFY THE MAXIMUM Hi TO DETERMINE FAILURE MODE
		  IHMAX=1
		  IF (H2.GT.H1) IHMAX=2
		  IF (IHMAX.EQ.1) THEN
		    IF (H6.GT.H1) IHMAX=6
		  END IF
		  IF (IHMAX.EQ.2) THEN
			IF (H6.GT.H2) IHMAX=6
		  END IF
C
C		  RAISE APPROPRIATE FLAG BASE ON MODE OF FAILURE
		  IF (IHMAX.EQ.1) REM(1)=1.0
		  IF (IHMAX.EQ.2) REM(2)=1.0
		  IF (IHMAX.EQ.6) REM(3)=1.0
C
		END IF
C
C	  MAXIMUM STRESS(2)/STRAIN(3) CRITERION
	  CASE(2)
		REM(1)=X0/XTER
		REM(2)=Y0/YTER
		REM(3)=ABS(Z0)/Q
		IF(REM(1).GT.1.0D0) REM(1)=1.0D0
		IF(REM(2).GT.1.0D0) REM(2)=1.0D0
		IF(REM(3).GT.1.0D0) REM(3)=1.0D0
C
C	  TSAI-HILL CRITERION
	  CASE(3)
		F=(X0/XTER)**2-X0*Y0/(XTER)**2+(Y0/YTER)**2+(Z0/Q)**2
		IF (F.GT.(1.0)) THEN
C
          H1=(X0/XTER)**2
		  H2=(Y0/YTER)**2
		  H6=(Z0/Q)**2
C
C		  IDENTIFY THE MAXIMUM Hi TO DETERMINE FAILURE MODE
		  IHMAX=1
		  IF (H2.GT.H1) IHMAX=2
		  IF (IHMAX.EQ.1) THEN
		    IF (H6.GT.H1) IHMAX=6
		  END IF
		  IF (IHMAX.EQ.2) THEN
			IF (H6.GT.H2) IHMAX=6
		  END IF
C
C		  RAISE APPROPRIATE FLAG BASE ON MODE OF FAILURE
		  IF (IHMAX.EQ.1) REM(1)=1.0
		  IF (IHMAX.EQ.2) REM(2)=1.0
		  IF (IHMAX.EQ.6) REM(3)=1.0
C
		END IF
C
C	  MODIFIED PUCK CRITERION
	  CASE(4)
		A66=1.0/(Q*Q)
		B2=1.0/YT-1.0/YC
		A22=1.0/(YT*YC)
		R1=X0/XTER
		R2=A22*Y0**2+A66*Z0**2+B2*Y0
		IF (R2.GE.(1.0)) R6=1.0
C
CC	  TSAI-WU CRITERION
	  CASE DEFAULT
		F=(1.0/XT-1.0/XC)*X0+(1.0/YT-1.0/YC)*Y0
	    F=F+((X0)**2)/(XT*XC)+((Y0)**2)/(YT*YC)
	    F=F+(Z0/Q)**2
		IF (F.GT.(1.0)) THEN
C
	      H1=(1.0/XT-1.0/XC)*X0+((X0)**2)/(XT*XC)
		  H2=(1.0/YT-1.0/YC)*Y0+((Y0)**2)/(YT*YC)
		  H6=(Z0/Q)**2
C
C		  IDENTIFY THE MAXIMUM Hi TO DETERMINE FAILURE MODE
		  IHMAX=1
		  IF (H2.GT.H1) IHMAX=2
		  IF (IHMAX.EQ.1) THEN
		    IF (H6.GT.H1) IHMAX=6
		  END IF
		  IF (IHMAX.EQ.2) THEN
			IF (H6.GT.H2) IHMAX=6
		  END IF
C
C		  RAISE APPROPRIATE FLAG BASE ON MODE OF FAILURE
		  IF (IHMAX.EQ.1) REM(1)=1.0
		  IF (IHMAX.EQ.2) REM(2)=1.0
		  IF (IHMAX.EQ.6) REM(3)=1.0
C
		END IF
	END SELECT
C
	RETURN
	END

C     ======================================================================
C     ======================================================================
C     ======================================================================

	SUBROUTINE COMPSHEAR(REM,PNL,LSSK,LAYS,TLIST,ALIST,ZT,TH,ISAND,G)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	PURPOSE
C	   TO CALCULATE STIFFNESSES INCLUDING THE IN?PLANE,
C	   FLEXURAL AND COUPLING STIFFNESSES FROM LAYER PROPERTIES.
C	   NOTE: THIS SUBROUTINE IS IDENTICAL TO SUBROUTINE S05STF IN
C	   DATA ITEM NO. 86020, ES DU PAC 1051. THE SUBROUTINE ALSO
C	   APPEARS IN DATA ITEM N0.88015.
C	USAGE
C	   CALL C8713A(LS,PNL,KOUNT,A,BC,D,LAYS,ZZ,TT)
C	   CALLED IN PROGRAM A8913A, DATA ITEM N0.89013
C	DESCRIPTION OF PARAMETERS
C	   KOUNT -NUMBER OF PLATES) IN DATA SEQUENCE?LIMITED TO 100
C	   LS	 -ARRAY OF NUMBER OF LAYERS IN LAMINATE FOR EACH PLATE
C			  IN SEQUENCE -LIMITED TO 50 PLATES
C	   LAYS	 -ARRAY OF STACKING SEQUENCES (ORDERED WITH LS) 50 OFF
C			  LIMITED TO 500 LAYERS IN EACH
C	   PNL	 -ARRAY STORING LAYER PROPERTIES -LIMITED TO 50. DATA
C			  MUST BE IN THE FOLLOWING ORDER:-	LAYER IDENTIFICATION
C			  NUMBER, EAA,EBB,EZZ,VAB,VAZ,VBZ,GBZ,GAZ,GAB,THICKNESS
C			  ANGLE (DEGREES, ALWAYS POSITIVE)
C	   A	 -IN-PLANE STTFFNFSS MATRTX
C	   BC	 -COUPLED IN?PLANE AND FLEXURAL STIFFNESS MATRIX
C	   D	 -FLEXURAL STIFFNESS MATRIX
C	   ST	 -CUMULATIVE LAYER THICKNESS FROM ONE SIDE
C			  INITIALLY SET TO ZERO
C	   K	 -CURRENT LAYER NUMBER IN CALCULATION
C	   KS	 -CURRENT LAYER TYPE IDENTIFICATION NUMBER
C	   PHI	 -ANGLE OF LAYER AXIS RADIANS
C	   ZK	 -DISTANCE OF LAYER MID-PLANE FROM MID?PLANE OF PLATE
C	   T	 -LAYER THICKNESS
C	   TT	 -TOTAL THICKNESS OF LAMINATE
C	   AT,DT -LARGEST TERMS TN A AND D RESPECTIVELY
C	   AZ	 -TERMS OF A WHICH ARE ZERO
C	   BCZ	 -TERMS OF B WHICH ARE ZERO
C	   DZ	 -TERMS OF D WHICH ARE ZERO
C	   ZZ	 -ARRAY OF LAYER Z VALUES WITH RESPECT TO PLATE MID PLANE
C	   DEN   -DENSITY OF COMPOSITE
C	REMARKS
C	   NONE
C	SUBROUTINES AND FUNCTIONS REQUIRED
C	   NONE
C	METHOD
C	   METHOD FOLLOWS DATA ITEM N0.75002. THE COMMENTS
C	   INDICATE THE EQUATION NUMBERS IN THAT DATA ITEM.
C
      
	DIMENSION TLIST(1),ALIST(1),PNL(20,1),LAYS(1),ZT(1)
	DIMENSION REM(3,1)
	
	DIMENSION BK(3,3),QK(2,2),BTRN(3,3),QTRN(2,2)
	DIMENSION A(3,3),BC(3,3),D(3,3),G(2,2)

C	FIRST SET CUMULATIVE THICKNESS FROM ONE SIDE TO ZERO AND
C	SET MATRICES TO ZERO BEFORE CALCULATION OF STIFFNESSES
       A(1:3,1:3) = 0.0D0
      BC(1:3,1:3) = 0.0D0
       D(1:3,1:3) = 0.0D0
       G(1:2,1:2) = 0.0D0    !SHEAR RIGIDITY IF ISAND=1 (CONSTANT SHEAR)
C
C	LOOP THROUGH LAYERS PICKING APPROPRIATE LAYER ACCORDING TO CODED
C	-------------------
	DO 90 K=1,LSSK
C	-------------------

C	  CHECK TRANVERSE RESISTANCE FAILURE
	  IF ((REM(1,K)+REM(2,K)).EQ.2.0) GOTO 90

C	  BUILD SINGLE LAYER STRESS STRAIN RELATIONSHIP MATRIX CK (EQ(1))
	  KS=LAYS(K)
		      
C     MATERIAL ANGLE W.R.P TO R AXIS (RAD) 
	  PHI= ALIST(K)
        
C	CALCULATE LAYER MATRIX BK (EQ(3))
        CALL QIJ(BK,QK,BTRN,QTRN,PNL(1,KS),PHI,'GOBL') 
        
C	  LAYER DISTANCE FROM NEUTRAL AXIS OF PLATE ZK (BOTTOM TO TOP)
        ZK = ZT(K)*TH

C	  CALCULATE LAYER A, B, D VALUES
	  DT=TLIST(K)*TH
		TB1 = DT
		TB2 = DT*ZK
		TB3 = DT*ZK*ZK + DT*DT*DT/12.0
C
C	MEMBRANE RIGIDITY
		DO 30 I=1,3
		DO 30 J=1,3
			 A(I,J) =  A(I,J) + BK(I,J)*TB1
			BC(I,J) = BC(I,J) - BK(I,J)*TB2
			 D(I,J) =  D(I,J) + BK(I,J)*TB3
30		CONTINUE
C
C	SHEAR RIGIDITY OF SANDWICH PLATE - CONVENTIONAL SOLUTION
        TQ1 = 5.0/6.0*DT
		DO 40 I=1,2
		DO 40 J=1,2
			 G(I,J) = G(I,J) + QK(I,J)*TQ1
40		CONTINUE

C	-------------------
90	CONTINUE
C	-------------------


C	FIND LARGEST TERM IN STIFFNESS SUBMATRICES (AT OF A, DT OF D, AND
C	BCT THE GEOMETRIC MEAN OF AT AND DT)
	AT=0.0
	DT=0.0
	DO 110 I=1,3
	  DO 110 J=1,3
		IF(ABS(A(I,J)).GT.AT) AT=ABS(A(I,J))
		IF(ABS(D(I,J)).GT.DT) DT=ABS(D(I,J))
110	CONTINUE
	BCT=SQRT(AT*DT)

C	ZEROISE NEGLIGIBLE TERMS OF STIFFNESS SUB?MATRICES (MADE ZERO IF
C     6 ORDERS LESS THAN AT, BCT OR DT FOR A, BC AND D RESPECTIVELY)
C	ALSO SUM TERMS WHICH SHOULD BE ZERO FOR SPECIAL ORTHOTROPY
	DO 120 I=1,3
	  DO 120 J=1,3
		IF(((ABS(A(I,J)))/AT).LT.10E-6) A(I,J)=0.0
		IF(((ABS(BC(I,J)))/BCT).LT.10E-6) BC(I,J)=0.0
		IF(((ABS(D(I,J)))/DT).LT.10E-6) D(I,J)=0.0
120	CONTINUE

      IF(ISAND.EQ.0) CALL ROHWER(REM,PNL,LSSK,LAYS,TLIST,ALIST,TH,A,BC,D,G)
      
	RETURN

	END
C
C     ======================================================================
C     ======================================================================
C     ======================================================================
      SUBROUTINE ROHWER(REM,PNL,LSSK,LAYS,TLIST,ALIST,TH,A,BC,D,G)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     A,BC,D  :  PLATE RIGIDITY FROM C8913A
C     PNL     :  INPUT DATA PNL(50,12)
C     G       :  SHEAR STIFFNESS G(2,2) :: A44,A45,A55
C     LAYSS   :  LAYSS - ARRAY STORING LAY-UP CODE FOR LAMINATES
C     LSSK    :  NUMBER OF LAYERS IN STACKING SEQUENCE
	DIMENSION REM(3,1),TLIST(1),ALIST(1),PNL(20,1),LAYS(1)
      DIMENSION A(3,3),BC(3,3),D(3,3),G(2,2)
      DIMENSION A1(3,3),B1(3,3),D1(3,3),A1I(3,3)

      DIMENSION AZ(3,3),BZ(3,3),DZ(3,3),FZ(2,2),GZ(2,2)
      DIMENSION DZI(3,3),GZI(2,2),FZ1(3,3),FZ2(3,3)
      DIMENSION HH(2,2),HHI(2,2)
      DIMENSION TEM(3,3),TEMP(3,3)   
        
	DIMENSION QMAT(3,3),EMAT(2,2),QTRN(3,3),ETRN(2,2)
	DIMENSION THZ(LSSK+1)
	

C     NUMBER OF INTEGRATION POINTS THROUGH THE LAYER THICKNESS 
      NN=100     
	
C
      A1=A
      B1=BC             
      D1=D
C     
      FZ=0.0D0
      DZ=0.0D0
 
	HH=0.0D0
	AZ=0.0D0
	BZ=0.0D0
	
	
C     INVERSE A1 FOR EQ(16)
      CALL INVERSE3 (A1,A1I)

C     EQ. (16)  ::  D* == DZ
      TEM=MATMUL(TRANSPOSE(B1),A1I)
	TEMP=MATMUL(TEM,B1)
	DZ = D1 - TEMP

C     INVERSE DZ FOR EQ.(18) 
      CALL INVERSE3(DZ,DZI)


C     CALCULATE THE DISTANCE OF EACE LAYER FROM THE MID-SURFACE FOR INTEGRATION
      THZ = 0.0D0
      THZ(1) = -TH/2.0D0
      DO I=1,LSSK
        THZ(I+1) = THZ(I) + TLIST(I)*TH
      END DO
	
	
C     DO LOOP FOR NUMBER OF LAYERS
C     EQ. (19) AND EQ. (20) ::  A(Z),B(Z)	
C     -------------------------  
	DO 1000 N=1,LSSK
C     -------------------------

C	  CHECK TRANVERSE RESISTANCE FAILURE
	  IF ((REM(1,K)+REM(2,K)).EQ.2.0) GOTO 1000
	
	  AZ=0.0D0
	  BZ=0.0D0
	  
	  IF(N.NE.1) THEN
            DO I=1,3
               DO J=1,3
		          DO K=1,N-1
		        
                    KS=LAYS(K)
                    PHI= ALIST(K)
		            CALL QIJ(QMAT,EMAT,QTRN,ETRN,PNL(1,KS),PHI,'GOBL') !Q MATRIX IN EQ.(19) AND EQ.(20)
                    AZ(I,J)=AZ(I,J)+QMAT(I,J)*(THZ(K+1)-THZ(K))
                    BZ(I,J)=BZ(I,J)+QMAT(I,J)*(THZ(K+1)*THZ(K+1) - THZ(K)*THZ(K))/2.0D0
                    
                  ENDDO
               ENDDO
            ENDDO 
        ENDIF	  

        KS=LAYS(N)
        PHI= ALIST(N)
		CALL QIJ(QMAT,EMAT,QTRN,ETRN,PNL(1,KS),PHI,'GOBL')

C     INTEGRATION POINT :: NN=100 
C     -------------------
        DO 500 M=1,NN
C     -------------------
  
            T1=THZ(N)+DBLE(M  )*TLIST(N)*TH/DBLE(NN)
            T2=THZ(N)+DBLE(M-1)*TLIST(N)*TH/DBLE(NN)

            DO I=1,3
            DO J=1,3
                AZ(I,J)=AZ(I,J)+QMAT(I,J)*(T1-T2)
                BZ(I,J)=BZ(I,J)+QMAT(I,J)*(T1*T1-T2*T2)/2.0D0
            END DO
            END DO

C     EQ. (18) :: F(Z)	   
            FZ1=0.0D0 ; FZ1=MATMUL(MATMUL(AZ,A1I),B1)-BZ
            FZ2=0.0D0 ; FZ2=MATMUL(FZ1,DZI)

	 
C     EQ. (26) :: f(z)   	   
            FZ(1,1)=FZ2(1,1)
            FZ(1,2)=FZ2(3,2)
            FZ(2,1)=FZ2(3,1)
            FZ(2,2)=FZ2(2,2)
	
C     G MATRIX OF EQ. (33)	::  ROHWER USED TAU XZ : Q55
C                                           TAU YZ : Q44
C     INPUT TAU YZ
C           TAU XZ , THERFORE REARRANGE FOR ROHWER THEORY
            GZ(1,1)=EMAT(2,2) ; GZ(1,2)=EMAT(2,1)
            GZ(2,1)=EMAT(1,2) ; GZ(2,2)=EMAT(1,1)

C     INVERSE G
            CALL INVERSE2(GZ,GZI)
	
C     EQ. (33) :: FT*G-1*F	 
            HH = HH + MATMUL(TRANSPOSE(FZ),MATMUL(GZI,FZ))*(T1-T2)
 
C     -------------------
500     CONTINUE
C     -------------------

C     -------------------------
1000  CONTINUE
C     -------------------------


C     INVERSE H FOR EQ. (33) :: H-1
      CALL INVERSE2(HH,HHI)
	
	
C     SHEAR STIFFNESS OF EQ. (4) 
C     ROHWER USED     TAU XZ : Q55
C                     TAU YZ : Q44	

C     GENERAL FORMAT  TAU YZ : Q44
C                     TAU XZ : Q55 , THEREFORE REARRANGE THE G MATRIX 

C	G(1,1)=HHI(2,2) ; G(1,2)=HHI(2,1)
C	G(2,1)=HHI(1,2) ; G(2,2)=HHI(1,1)


C     NOTE :: XFINAS USE THE ROHWER FORMAT  !!!!!!!!!!!!!!!!!!!!!!!!!        
      G(1,1)=HHI(1,1) ; G(1,2)=HHI(1,2)
      G(2,1)=HHI(2,1) ; G(2,2)=HHI(2,2)


      RETURN
      END
C

C
C     ======================================================================
C     ======================================================================
C     ======================================================================
	SUBROUTINE QIJ(QMAT,EMAT,TT1,TT2,PNL,PHI,OPTN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 OPTN
	DIMENSION PNL(1)
	DIMENSION QMAT(3,3),EMAT(2,2)
	DIMENSION QM(3,3),EM(2,2),TT1(3,3),TT2(2,2)
     

      CC = COS(PHI) 
      SS = SIN(PHI)	
      	  		  			   	  		 						   	     	   
      TT1(1,1)=CC*CC ; TT1(2,1)= SS*SS ; TT1(3,1)=-2.0D0*SS*CC
      TT1(1,2)=SS*SS ; TT1(2,2)= CC*CC ; TT1(3,2)= 2.0D0*SS*CC
      TT1(1,3)=SS*CC ; TT1(2,3)=-SS*CC ; TT1(3,3)= CC*CC-SS*SS
      
      TT2(1,1)= CC ; TT2(2,1)= SS
      TT2(1,2)=-SS ; TT2(2,2)= CC


	Q16=0.0D0
	Q26=0.0D0
	Q45=0.0D0

	ANU2 = PNL(5)*PNL(3)/PNL(2)

	Q11=PNL(2)/(1.0-PNL(5)*ANU2)
	Q22=PNL(3)/(1.0-PNL(5)*ANU2)
	Q12=PNL(3)/(1.0-PNL(5)*ANU2)*PNL(5)
	Q66=PNL(10)
	Q44=PNL(8)
	Q55=PNL(9)

	QM(1,1)=Q11 ; QM(1,2)=Q12 ; QM(1,3)=Q16
	QM(2,1)=Q12 ; QM(2,2)=Q22 ; QM(2,3)=Q26
	QM(3,1)=Q16 ; QM(3,2)=Q26 ; QM(3,3)=Q66
	EM(1,1)=Q44 ; EM(1,2)=Q45
	EM(2,1)=Q45 ; EM(2,2)=Q55

	QMAT=0.0D0
	EMAT=0.0D0

      SELECTCASE(OPTN)
          CASE('GOBL')
	      QMAT=MATMUL(TRANSPOSE(TT1),MATMUL(QM,TT1))
	      EMAT=MATMUL(TRANSPOSE(TT2),MATMUL(EM,TT2))
          CASE('LOCL')
	      QMAT=QM
	      EMAT=EM
      ENDSELECT

	RETURN
	END


C
C     ======================================================================
C     ======================================================================
C     ======================================================================
