C=====================================================================
      SUBROUTINE SPCABLE(PROPM,PROPG,NODEX,SPCFR,S,COORD,EDIS,EDISI,RE,
     #                   MWG,FIN,CABFF)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------------
C     3D PARABOLIC CABLE  (REFERENCE FROM SPC-FRAME)
C     ---------------------------------------------------------------------
C
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV

C	==================================================================
C	CABLE PRETENSION OPTIMIZATION
	COMMON /CB556/ LCPZ,NCABZ,KBOPZ,NCOBZ
C	==================================================================


      DIMENSION PROPM(1),PROPG(5),EDIS(6),EDISI(6)
      DIMENSION COORD(6)
C
	DIMENSION HN(2,6),DISP(2,1),RE(6)
	DIMENSION S(21)

	DIMENSION VR(3),VS(3),VT(3),VRP(3)
	DIMENSION STIF(6,6)
	DIMENSION FIN(NEF)
	DIMENSION YMCUV(20,3),YCCUV(20,2)
	DIMENSION SPCFR(1)

C	CABLE FORCE FOR OPTIMIZATION
	DIMENSION CABFF(NELE)

      DIMENSION COORO(6),REO(6),VRO(3),VSO(3),VTO(3)



C	-------------------------
C	ORIGINAL LENGTH AND ANGLE
C	-------------------------

	COORO(1:6) = COORD(1:6)
	IF(NLOPT.GT.1) COORO(1:6) = COORD(1:6) - EDIS(1:6)

C	--------------------------------
C	FOR LINEAR AND MATERIAL NONLINEAR ONLY - UPDATE THE ELEMENT NODAL COORDINATE BY DISP. TERMS
C	--------------------------------
	IF(NLOPT.LE.1) THEN
	DO I = 1,6
	COORD(I) = COORD(I)+EDIS(I)
	ENDDO
	ENDIF

C	--------------------------------
C	INITIALIZATION OF SOME VARIABLES
C	--------------------------------

	TDEN = PROPM(5)
	WDEN = PROPM(6)

C	-----------------------------------------
C	TO COMPUTE THE CHORD-LENGTH OF EACH CABLE
C	-----------------------------------------
	A1=COORD(4)-COORD(1)
	A2=COORD(5)-COORD(2)
	A3=COORD(6)-COORD(3)
	WL=SQRT(A1**2+A2**2+A3**2)

C	-------------------------
C	ANGLE WITH THE HORIZONTAL
C	-------------------------
	IF(NGRAV.EQ.1) THEN
	ANGL = ATAN(A1/SQRT(A2**2+A3**2))  !X IS VERTICAL AXIS
	ELSEIF(NGRAV.EQ.2) THEN
	ANGL = ATAN(A2/SQRT(A1**2+A3**2))  !Y IS VERTICAL AXIS
	ELSEIF(NGRAV.EQ.3) THEN
	ANGL = ATAN(A3/SQRT(A1**2+A2**2))  !Z IS VERTICAL AXIS
	ENDIF

	WDEN = WDEN*COS(ANGL)

C	---------------------------------------
C	SET VALUES FOR LINEAR STRESS-STRAIN LAW 
C	INITIALISATION OF INTEGRATION RULE
C	---------------------------------------
C	CALL HOKLAW (PROPM,PROPG,1)


	AREA = PROPG(2)

	IFOPT = INT(PROPG(4))  !IFOPT    1=INITIAL FORCE WILL ADD TO LOAD VECTOR    0=INITIAL FORCE WILL NOT ADD TO LOAD VECTOR FEB09

	FORCO  = SPCFR(3) ! INITIAL TENSION FORCE

C	FOR CABLE OPTIMIZATION USING THE FORCE FROM ADS PROGRAM
	IF(KBOPZ.EQ.1.AND.KSTEP.NE.1) FORCO = CABFF(MEL)

	STRAO  = SPCFR(4) ! PREVIOUS STRAIN


	SIGMO = FORCO/AREA


	MP = PROPM(14)
C	----------------------------
C	GENERATE STRESS-STRAIN CURVE
C	----------------------------
C	CALL STGEN(PROPM,YMCUV,YCCUV,SIGMO,WDEN,WL,MP)
	CALL STGEN_MO(PROPM(15),YMCUV,YCCUV,SIGMO,WDEN,WL,MP,STRAO)

C	------------
C	BASE VECTORS
C	------------
	CALL TRVRST(COORD,VR,VS,VT)

C	------------------------
C	FORM MASS MATRIX
C	------------------------
	IF (ITASK.NE.5) GOTO 50
	SSAM = TDEN*AREA*WL
	IMASS = 1
	IF (IMASS.EQ.1) THEN
C	  LUMPED MASS MATRIX
C	  ------------------
	  FACT = SSAM/2.0D0
	  S( 1) = FACT
	  S( 7) = FACT
	  S(12) = FACT
	  S(16) = FACT
	  S(19) = FACT
	  S(21) = FACT
	ELSE
C	  CONSISTENT MASS MATRIX
C	  ----------------------
	  FACT  = SSAM/3.0D0
	  S( 1) = FACT
	  S( 7) = FACT
	  S(12) = FACT
	  S(16) = FACT
	  S(19) = FACT
	  S(21) = FACT
	  FACT  = SSAM/6.0D0
	  S( 4) = FACT
	  S(10) = FACT
	  S(15) = FACT
	ENDIF
	RETURN

50	CONTINUE  !IF (NLOPT+ITASK.EQ.1) GOTO 1000



C	------------
C	ORIGINAL BASE VECTORS
C	------------
	CALL TRVRST(COORO,VRO,VSO,VTO)

	A10 = COORO(4) - COORO(1)
	A20 = COORO(5) - COORO(2)
	A30 = COORO(6) - COORO(3)
	WL0=SQRT(A10**2+A20**2+A30**2)

C	--------------
C	ELEMENT STRAIN
C	--------------
	STRAN  = (WL-WL0)/WL0 !current strain

	IF(NLOPT.LE.1) THEN
	DIS1 = VR(1)*EDIS(1) + VR(2)*EDIS(2) + VR(3)*EDIS(3)
	DIS2 = VR(1)*EDIS(4) + VR(2)*EDIS(5) + VR(3)*EDIS(6)
	STRAN = (DIS2-DIS1)/WL
	ENDIF


	DSTRN  = STRAN 

	CALL CABST(YMCUV,YCCUV,STRAO,DSTRN,WDEN,WL,STRESS,MP,MCODS,
	1		   EMODS)

	SPCFR(4) = STRAO

C	----------------------
C	UPDATE TANGENT MODULUS
C	----------------------
	DO I = 1,6
	RE(I) = 0.0
	DO J = 1,6
	STIF(I,J) = 0.0
	ENDDO
	ENDDO


      IF (MCODS.EQ.1) THEN
        ETAN  = 1.D-20
	  P  = 0.0
	  P1 = 0.0
        GOTO 500
      ELSEIF (MCODS.EQ.2) THEN		
        ETAN  =  EMODS/ (1.0+ EMODS * (WDEN * WL0)**2.0
     1        / (12.0 * STRESS**3.0))
      END IF


C	FOR CABLE OPTIMIZATION USING TRUSS FORCE FOR FIRST ESTIMATE
C	IF(KBOPZ.EQ.1) THEN
C	ETAN   = EMODS
C	STRESS = SIGMO+STRAN*EMODS
C	IF(STRESS.LT.0.0) STRESS = 0.0
C	IF(STRESS.LT.0.0) ETAN   = 0.0
C	ENDIF


C	WRITE(*,*) 'LL',EDIS(4),SIGMO,STRESS
C	--------------------------------
C	COMPUTE INTERNAL RESISTING FORCE
C	--------------------------------
	P  = STRESS*AREA

	IF(STRAN.EQ.0.0D0) P = FORCO 

	P1 = P
	P0 = 0.0D0
	IF(IFOPT.EQ.0) THEN !IFOPT    0=INITIAL FORCE WILL NOT ADD TO LOAD VECTOR    1=INITIAL FORCE WILL DD TO LOAD VECTOR
	P0 = FORCO  
      ENDIF
      
C	TAKE OUT INITIAL FORCE FROM ELEMENT LOAD VECTOR INCASE OF CABLE OPTIMIZATION
	IF(KBOPZ.EQ.1.AND.KSTEP.NE.1) P0 = FORCO

	IF(P1.LT.0.0) P1 = 0.0D0


500	CONTINUE

	DP = P - FORCO

C	WRITE(*,*) MEL,P,DP,FORCO

	SPCFR(1) = P
	SPCFR(2) = DP


	REO(1) = -P0*VRO(1)
	REO(2) = -P0*VRO(2)
	REO(3) = -P0*VRO(3)
	REO(4) =  P0*VRO(1)
	REO(5) =  P0*VRO(2)
	REO(6) =  P0*VRO(3)

	RE(1) = -P1*VR(1) - REO(1)
	RE(2) = -P1*VR(2) - REO(2)
	RE(3) = -P1*VR(3) - REO(3)
	RE(4) =  P1*VR(1) - REO(4)
	RE(5) =  P1*VR(2) - REO(5)
	RE(6) =  P1*VR(3) - REO(6)

	DO 20 I = 1,NEF
	    FIN(I) = RE(I)
20	CONTINUE

C	----------------------------
C	COMPUTE FOR STIFFNESS MATRIX
C	----------------------------
1000	IF (IFREF) 900,750,900

750	CONTINUE

      AEL  = ETAN * AREA / WL
      FL  =  P/WL


      STIF(1,1) = VR(1)**2*AEL+VS(1)**2*FL+VT(1)**2*FL
      STIF(1,2) = VR(1)*AEL*VR(2)+VS(1)*FL*VS(2)+VT(1)*FL*VT(2)
      STIF(1,3) = VR(1)*AEL*VR(3)+VS(1)*FL*VS(3)+VT(1)*FL*VT(3)
      STIF(1,4) = -VR(1)**2*AEL-VS(1)**2*FL-VT(1)**2*FL
      STIF(1,5) = -VR(1)*AEL*VR(2)-VS(1)*FL*VS(2)-VT(1)*FL*VT(2)
      STIF(1,6) = -VR(1)*AEL*VR(3)-VS(1)*FL*VS(3)-VT(1)*FL*VT(3)
      STIF(2,1) = VR(1)*AEL*VR(2)+VS(1)*FL*VS(2)+VT(1)*FL*VT(2)
      STIF(2,2) = VR(2)**2*AEL+VS(2)**2*FL+VT(2)**2*FL
      STIF(2,3) = VR(2)*AEL*VR(3)+VS(2)*FL*VS(3)+VT(2)*FL*VT(3)
      STIF(2,4) = -VR(1)*AEL*VR(2)-VS(1)*FL*VS(2)-VT(1)*FL*VT(2)
      STIF(2,5) = -VR(2)**2*AEL-VS(2)**2*FL-VT(2)**2*FL
      STIF(2,6) = -VR(2)*AEL*VR(3)-VS(2)*FL*VS(3)-VT(2)*FL*VT(3)
      STIF(3,1) = VR(1)*AEL*VR(3)+VS(1)*FL*VS(3)+VT(1)*FL*VT(3)
      STIF(3,2) = VR(2)*AEL*VR(3)+VS(2)*FL*VS(3)+VT(2)*FL*VT(3)
      STIF(3,3) = VR(3)**2*AEL+VS(3)**2*FL+VT(3)**2*FL
      STIF(3,4) = -VR(1)*AEL*VR(3)-VS(1)*FL*VS(3)-VT(1)*FL*VT(3)
      STIF(3,5) = -VR(2)*AEL*VR(3)-VS(2)*FL*VS(3)-VT(2)*FL*VT(3)
      STIF(3,6) = -VR(3)**2*AEL-VS(3)**2*FL-VT(3)**2*FL
      STIF(4,1) = -VR(1)**2*AEL-VS(1)**2*FL-VT(1)**2*FL
      STIF(4,2) = -VR(1)*AEL*VR(2)-VS(1)*FL*VS(2)-VT(1)*FL*VT(2)
      STIF(4,3) = -VR(1)*AEL*VR(3)-VS(1)*FL*VS(3)-VT(1)*FL*VT(3)
      STIF(4,4) = VR(1)**2*AEL+VS(1)**2*FL+VT(1)**2*FL
      STIF(4,5) = VR(1)*AEL*VR(2)+VS(1)*FL*VS(2)+VT(1)*FL*VT(2)
      STIF(4,6) = VR(1)*AEL*VR(3)+VS(1)*FL*VS(3)+VT(1)*FL*VT(3)
      STIF(5,1) = -VR(1)*AEL*VR(2)-VS(1)*FL*VS(2)-VT(1)*FL*VT(2)
      STIF(5,2) = -VR(2)**2*AEL-VS(2)**2*FL-VT(2)**2*FL
      STIF(5,3) = -VR(2)*AEL*VR(3)-VS(2)*FL*VS(3)-VT(2)*FL*VT(3)
      STIF(5,4) = VR(1)*AEL*VR(2)+VS(1)*FL*VS(2)+VT(1)*FL*VT(2)
      STIF(5,5) = VR(2)**2*AEL+VS(2)**2*FL+VT(2)**2*FL
      STIF(5,6) = VR(2)*AEL*VR(3)+VS(2)*FL*VS(3)+VT(2)*FL*VT(3)
      STIF(6,1) = -VR(1)*AEL*VR(3)-VS(1)*FL*VS(3)-VT(1)*FL*VT(3)
      STIF(6,2) = -VR(2)*AEL*VR(3)-VS(2)*FL*VS(3)-VT(2)*FL*VT(3)
      STIF(6,3) = -VR(3)**2*AEL-VS(3)**2*FL-VT(3)**2*FL
      STIF(6,4) = VR(1)*AEL*VR(3)+VS(1)*FL*VS(3)+VT(1)*FL*VT(3)
      STIF(6,5) = VR(2)*AEL*VR(3)+VS(2)*FL*VS(3)+VT(2)*FL*VT(3)
      STIF(6,6) = VR(3)**2*AEL+VS(3)**2*FL+VT(3)**2*FL


C	--------------------------------------------------
C	COMPUTE FULL STIFFNESS MATRIX S(21) UPPER-TRIANGLE
C	--------------------------------------------------
1100	KK=1
	DO 700 I=1,6
	DO 700 J=I,6
	S(KK)=STIF(I,J)
	KK=KK+1
700	CONTINUE
	

C	TIM(12) = TIM(12)+TIM2

900	CONTINUE

C	STORE THE CABLE TOTAL FORCE FOR EACH STEP (ADS OPTIMIZATION)
	IF(KBOPZ.EQ.1) CABFF(MEL) = P


	RETURN


	END
C
C	=================================================================
C	=================================================================
C	=================================================================
C	=================================================================
      SUBROUTINE STGEN(PROPM,YMCUV,YCCUV,SIGMO,WDEN,WL,MP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=================================================================
C	GENERATE THE STRESS-STRAIN CURVE FOR CABLE MATERIAL
C	GENERATE THE STRESS-APPARENT STRAIN CURVE FOR CABLE MATERIAL
C	=================================================================
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV

	DIMENSION PROPM(NMP)
	DIMENSION YMCUV(20,3),YCCUV(20,2)


	DO I = 1,MP
	YMCUV(I,2) = PROPM(5+I)     !STRESS
	YMCUV(I,1) = PROPM(5+I+MP)  !STRAIN
	ENDDO

C	STRESS-STRAIN CURVE
	S1 = 0.0
	E1 = 0.0
	DO I = 1,MP
	E2 = YMCUV(I,1)
	S2 = YMCUV(I,2)
	EY = (S2-S1)/(E2-E1)
	YMCUV(I,3) = EY
	S1 = S2
	E1 = E2
	ENDDO


C	STRESS-APPARENT STRAIN CURVE

	CONST = WDEN*WDEN*WL*WL/24.0

	YCCUV(1,1) = 0.0
	YCCUV(1,2) = SIGMO



	S1 = SIGMO
	DO I = 1,MP
	S2 = YMCUV(I,2)
	E  = YMCUV(I,3)
	
	CN1 = (S2-S1)/E
	CN2 = (1.0/(S1*S1)) - (1.0/(S2*S2))

	YCCUV(I+1,1) = CN1 + (CONST*CN2) + YCCUV(I,1)
	YCCUV(I+1,2) = S2

	S1 = S2

	ENDDO

	
	RETURN

	END

C	======================================================================
C	======================================================================
      SUBROUTINE CABST(YMCUV,YCCUV,STRAO,DSTRN,WDEN,WL,STRESS,MP,MCODS,
	1				 EMODS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV

	DIMENSION PROPM(NMP)
	DIMENSION YMCUV(20,3),YCCUV(20,2)

	CONST = WDEN*WDEN*WL*WL/24.0

	STRNC = STRAO + DSTRN


C	FOR CABLE IN COMPRESSION
	IF(STRNC.LT.0.0) THEN
	STRESS = 0.0
	MCODS = 1
	EMODS = 0.0
	GO TO 1000
	ENDIF

C	START OF CABLE IN TENSION
	MCODS = 2

C	THE FIRST ENTRY STEP
	IF(STRNC.EQ.0.0) THEN
	STRESS = YCCUV(1,2)
	EMODS  = YMCUV(1,3)
	GO TO 1000
	ENDIF


C	PICK THE POINT WHICH THE CURRENT STRAIN LAYING ON
	IE = 0
	DO 10 I = 1,MP+1
	IF(STRNC.GT.YCCUV(I,1)) GO TO 10
	IE = I-1
	E  = YMCUV(IE,3)
	E1 = YCCUV(IE,1)
	S1 = YCCUV(IE,2)
	SE = YCCUV(I ,2)
	EMODS = E
	GO TO 15
10	CONTINUE

15	CONTINUE
		
	
	IF(IE.EQ.0) THEN ! ASSUME PERFECT PLASTICITY FOR RANGE OVER THE DATA
	E  = 1.0E-20
	E1 = YCCUV(MP,1)
	S1 = YCCUV(MP,2)
	SE = S1
	EMODS = 1.0E-20
	ENDIF

	S12 = S1*S1

	DELTL = (STRNC - E1)

	

C	SOLVE THE CURRENT STRESS BY "BI-SECTION ITERATIVE METHOD"
	ITERM = 200
	TOL = 0.0000001
	S1I = S1
	S2I = SE
	STSO = 0.0
	DO I = 1,ITERM

	STS = 0.5*(S1I + S2I)

	CN1 = (STS-S1)/E
	CN2 = (1.0/(S1*S1)) - (1.0/(STS*STS)) 

	FUNC = CN1 + CONST*CN2 - DELTL

	TEST = ABS((STS-STSO)/STS)
	IF(TEST.LE.TOL) THEN
	STRESS = STS
	GO TO 1000
	ENDIF

	STSO = STS
	IF(FUNC.GT.0.0) S2I = STS
	IF(FUNC.LE.0.0) S1I = STS

	ENDDO
	
C	IF SOLUTION IS NOT CONVERGE USE THE LINEAR INTERPOLATION

	STRESS = (STRNC-E1)*E + S1


1000	CONTINUE
	

	RETURN

	END

C	======================================================================
C	======================================================================
C	======================================================================
C	=================================================================
      SUBROUTINE STGEN_MO(PROPM,YMCUV,YCCUV,SIGMO,WDEN,WL,MP,STRAO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=================================================================
C	GENERATE THE STRESS-STRAIN CURVE FOR CABLE MATERIAL
C	GENERATE THE STRESS-APPARENT STRAIN CURVE FOR CABLE MATERIAL
C	=================================================================
	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV

	DIMENSION YMCUV(20,3),YCCUV(20,2),PROPM(MP,2)


	DO I = 1,MP
	YMCUV(I,2) = PROPM(I,1)     !STRESS
	YMCUV(I,1) = PROPM(I,2)     !STRAIN
	ENDDO

C	STRESS-STRAIN CURVE
	S1 = 0.0
	E1 = 0.0
	DO I = 1,MP
	E2 = YMCUV(I,1)
	S2 = YMCUV(I,2)
	EY = (S2-S1)/(E2-E1)
	YMCUV(I,3) = EY
	S1 = S2
	E1 = E2
	ENDDO


C	STRESS-APPARENT STRAIN CURVE

	CONST = WDEN*WDEN*WL*WL/24.0

	YCCUV(1,1) = 0.0
	YCCUV(1,2) = 0.3*SIGMO


	S1 = 0.3*SIGMO
	DO I = 1,MP
	S2 = YMCUV(I,2)
	E  = YMCUV(I,3)
	
	CN1 = (S2-S1)/E
	CN2 = (1.0/(S1*S1)) - (1.0/(S2*S2))

	YCCUV(I+1,1) = CN1 + (CONST*CN2) + YCCUV(I,1)
	YCCUV(I+1,2) = S2

	S1 = S2

	ENDDO


	CN1 = (SIGMO-0.3*SIGMO)/YMCUV(1,3)
	CN2 = (1.0/(0.3*SIGMO*0.3*SIGMO)) - (1.0/(SIGMO*SIGMO))
	STRAO = CN1 + (CONST*CN2)

	
	RETURN

	END

C	======================================================================
C	OPTIMIZATION OF CABLE PRETENSION
C	=============================================================
C	=============================================================
	SUBROUTINE CABLEXF(RHF,CAB,FOC,NCABZ,NEQ)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	DIMENSION RHF(NEQ),CAB(NEQ,NCABZ),FOC(NCABZ)

C	WRITE(*,*) FOC

	RHF = MATMUL(CAB,FOC)

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE OBJECZ(B,ESTN,NEQ,MAXA,AA,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP
C	=============================================================
C	=============================================================
	DIMENSION B(NEQ),R(NEQ),AA(1)
	DIMENSION MAXA(NEQ+1)

	CALL MAMULT(MAXA,AA,B,R,TYP,'STD')  

	ESTN = 0.0D0
	DO IEQ = 1,NEQ
	ESTN = ESTN + 0.5*R(IEQ)*B(IEQ)
	ENDDO

	RETURN

C	DO I = 1,NEQ
C	KHIGH = MAXA(I+1) - MAXA(I)
C	LHEI(I) = I - KHIGH + 1
C	ENDDO
	ESTN = 0.0D0

C	=====================================================
	DO IEQ = 1,NEQ
	R(IEQ) = 0.0D0
C	-----------------------------------------------------
	DO JEQ = IEQ,NEQ	!UPPER TRIANGULAR
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1

	IF(LHEI.LE.IEQ) THEN
	R(IEQ) = R(IEQ) + AA(NUM)*B(JEQ)            !R(IEQ) + A(NUM)*B(JEQ)
	ENDIF

	ENDDO
C	-----------------------------------------------------
	NUM = MAXA(IEQ+1) - MAXA(IEQ)
	LCL = MAXA(IEQ)-1
	DO J = 2,NUM		!LOWER TRIANGULAR
	JEQ = IEQ - J + 1

	IF(JEQ.GT.0) THEN
	NUM2 = LCL+J
	R(IEQ) = R(IEQ) + AA(NUM2)*B(JEQ)           !R(IEQ) + A(LCL+J)*B(JEQ)
	ENDIF

	ENDDO
C	-----------------------------------------------------
	ESTN = ESTN + 0.5*R(IEQ)*B(IEQ)
C	-----------------------------------------------------
	ENDDO
C	=====================================================

	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================

	SUBROUTINE OBJECZ_OLD(DIS,AA,ESTN,NEQ,MAXA)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	DIMENSION DIS(NEQ,1),AA(1),STG(NEQ,NEQ),ETN(1,1),MAXA(NEQ+1)



	K = 0
	STG = 0.0
	DO I = 1,NEQ
	
	NHEI = MAXA(I+1)-MAXA(I)
	NHEI = I-NHEI+1 

	DO J = I,NHEI,-1

	K = K+1

	STG(J,I) = AA(K)

	
	ENDDO
	ENDDO


	DO I = 1,NEQ
	DO J = I,NEQ
	STG(J,I) = STG(I,J)
	ENDDO
	ENDDO

	

	ETN = 0.5*MATMUL(TRANSPOSE(DIS),MATMUL(STG,DIS))

	ESTN = ETN(1,1)



111	FORMAT(100E15.6)	
	
	RETURN
	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE CABCON(NCABZ,NCOBZ)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

C	-------------------------------------------------------------
	COMMON /CB557I/ MCSUM,MCFOC,MCDEF,MCMOM,LS557,LF557,LD557,LM557,
	1			    LR557,NM556
	COMMON /CB557R/ TPOPZ(10000)
C	-------------------------------------------------------------

	READ(ITI,*) 
	READ(ITI,*) NM556,MCSUM,MCFOC,MCDEF,MCMOM

	NCOBZ = MCSUM+2*MCDEF+2*MCMOM             !NUMBER OF CONSTRAIN

	

	LT557 = 1
C	--------------------------------------
	LS557  = LT557
	LF557  = LS557 + MCSUM*(NCABZ+3+20)
	LD557  = LF557 + 3*MCFOC
	LM557  = LD557 + 3*MCDEF
	LR557  = LM557 + 4*MCMOM + 2*MCMOM
	LST557 = LR557 + NM556*NCABZ
C	--------------------------------------
	DO I = LT557,LST557
	TPOPZ(I) = 0.0D0
	ENDDO
C	--------------------------------------
	READ(ITI,*)
	CALL READREC(NM556,TPOPZ(LR557),ITI)	
C	--------------------------------------
	IF(MCSUM.NE.0) THEN
	READ(ITI,*) 
	CALL READSUM(MCSUM,TPOPZ(LS557),NCABZ,ITI)

	ENDIF
C	--------------------------------------	
	IF(MCFOC.NE.0) THEN
	READ(ITI,*) 
	CALL READFOC(MCFOC,TPOPZ(LF557),ITI)

	ENDIF
C	--------------------------------------
	IF(MCDEF.NE.0) THEN
	READ(ITI,*) 
	CALL READDEF(MCDEF,TPOPZ(LD557),ITI)	

	ENDIF
C	--------------------------------------
	IF(MCMOM.NE.0) THEN
	READ(ITI,*) 
	CALL READMOM(MCMOM,TPOPZ(LM557),ITI)	

	ENDIF
C	--------------------------------------

	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE READREC(NM556,MGRUP,ITI)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION MGRUP(NM556,1)

	DO I = 1,NM556
	READ(ITI,*) N,NN
	READ(ITI,*) MGRUP(I,1:NN)
	ENDDO

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE READSUM(MCSUM,PCSUM,NCABZ,ITI)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION PCSUM(MCSUM,NCABZ+3+20)

	DO IGR = 1,MCSUM
	READ(ITI,*) NELM,VALUE,NODE
	PCSUM(IGR,1) = NELM
	PCSUM(IGR,2) = VALUE
	PCSUM(IGR,3) = NODE
	READ(ITI,*) (PCSUM(IGR,II+3),II=1,NODE)           !NODE IN GROUP
	READ(ITI,*) (PCSUM(IGR,II+3+NODE),II=1,NELM)      !CABLE ELEMENT IN GROUP
	ENDDO


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE READFOC(MCFOC,PCFOC,ITI)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION PCFOC(MCFOC,3)


	DO II = 1,MCFOC
	READ(ITI,*) PCFOC(II,1),PCFOC(II,2),PCFOC(II,3)  !ELEM  MIN  MAX
	ENDDO

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE READDEF(MCDEF,PCDEF,ITI)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION PCDEF(MCDEF,3)

	DO II = 1,MCDEF
	READ(ITI,*) PCDEF(II,1),PCDEF(II,2),PCDEF(II,3)  !NODE  DOF  VALUE
	ENDDO

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE READMOM(MCMOM,PCMOM,ITI)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION PCMOM(MCMOM,4)

	DO II = 1,MCMOM
	READ(ITI,*) PCMOM(II,2),PCMOM(II,3),PCMOM(II,4)  !GROUP  ELEM  DOF  VALUE
	ENDDO

C	PCMOM(II,1)

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE CUTCON(IDG,NCOBZ)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	COMMON /CB557I/ MCSUM,MCFOC,MCDEF,MCMOM,LS557,LF557,LD557,LM557,
	1			    LR557,NM556

	DIMENSION IDG(NCOBZ)

	K = 0
	DO I = 1,MCSUM
	K = K+1
	IDG(K) = -1
	ENDDO

C	DO I = MCFOC  !ASSIGN THIS CONSTRAIN TO BE THE UPPER AND LOWER BOUND
C	K = K+1
C	IDG(K) = 0
C	ENDDO

	DO I = 1,2*MCDEF
	K = K+1
	IDG(K) = 0
	ENDDO

	DO I = 1,2*MCMOM
	K = K+1
	IDG(K) = 0
	ENDDO

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE BNDCON(VL556,VU556,PCFOC,IGIDMEM,MCFOC,NCABZ)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION PCFOC(MCFOC,3),IGIDMEM(NCABZ)
	DIMENSION VL556(NCABZ),VU556(NCABZ)


	DO II = 1,NCABZ
	VL556(II) = 0.0D0
	VU556(II) = 1.0E15
	ENDDO

	DO II = 1,MCFOC

	IGD = PCFOC(II,1)
	CALL EGIDNUM(IGIDMEM,IGD,NCABZ,IEL)
	IF(II.EQ.IEL) VL556(II) = PCFOC(II,2)
	IF(II.EQ.IEL) VU556(II) = PCFOC(II,3)

	ENDDO

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE DEFCON(G556,DISP,ID,MCSUM,MCDEF,
	1				  PCDEF,NCABZ,NCOBZ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
C	-------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL


	DIMENSION ID(NSF,NSN),PCDEF(MCDEF,3),G556(NCOBZ),DISP(NEQ)

	
C	NEGATIVE DISPLACEMENT CONSTRAINT
	DO 100 ICDEF = 1,MCDEF

	NODE = PCDEF(ICDEF,1)
	IDF  = PCDEF(ICDEF,2)
	VALV = PCDEF(ICDEF,3)

	G556(ICDEF+MCSUM) = VALV

	IQ = ID(IDF,NODE)    
	DIP = 0.0D0
	IF(IQ.NE.0) DIP = DISP(IQ)        
	G556(ICDEF+MCSUM) = -1.0*(G556(ICDEF+MCSUM) + DIP)

100	CONTINUE

C	POSITIVE DISPLACEMENT CONSTRAINT
	DO 200 ICDEF = 1,MCDEF

	NODE = PCDEF(ICDEF,1)
	IDF  = PCDEF(ICDEF,2)
	VALV = PCDEF(ICDEF,3)

	G556(ICDEF+MCSUM+MCDEF) = -1.0*VALV

	IQ = ID(IDF,NODE)     
	DIP = 0.0D0
	IF(IQ.NE.0) DIP = DISP(IQ)            
	G556(ICDEF+MCSUM+MCDEF) = (G556(ICDEF+MCSUM+MCDEF) + DIP)

200	CONTINUE	

	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE SUMCON(G556,CABFX,CBFOC,ID,MCSUM,PCSUM,IGIDMEM,
	1				  NCABZ,NCOBZ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
C	-------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL


	DIMENSION CABFX(NEQ,NCABZ),ID(NSF,NSN),PCSUM(MCSUM,1)
	DIMENSION CBFOC(NCABZ),RHX(NEQ),CML(NEQ),G556(NCOBZ)
	DIMENSION IGIDMEM(NCABZ)

	
	
	

	DO 100 ICSUM = 1,MCSUM

	NELE = PCSUM(ICSUM,1)
	VALV = PCSUM(ICSUM,2)
	NODE = PCSUM(ICSUM,3)

	RHX = 0.0D0
	DO IEL = 1,NELE
	NUME = PCSUM(ICSUM,IEL+3+NODE)
	CALL EGIDNUM(IGIDMEM,NUME,NCABZ,MEM)
	CALL CABMUL(CABFX,CBFOC,CML,NEQ,NCABZ,MEM)
	RHX = RHX + CML
	ENDDO

	G556(ICSUM) = -1.0*VALV
	DO IOD = 1,NODE
	ND = PCSUM(ICSUM,IOD+3)
	IQ = ID(1,ND)                !X-DIRECTION FORCE
	G556(ICSUM) = G556(ICSUM) + RHX(IQ)
	ENDDO


100	CONTINUE

C	WRITE(*,*) G556

	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE EGIDNUM(IGIDMEM,IGD,NCABZ,IEL)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION IGIDMEM(NCABZ)

	IEL = 0
	DO I = 1,NCABZ
	IGID = IGIDMEM(I)
	IF(IGID.EQ.IGD) IEL = I
	ENDDO

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE CABMUL(CABFX,CBFOC,CML,NEQ,NCABZ,NUME)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION CABFX(NEQ,NCABZ),CML(NEQ),CBFOC(NCABZ)

	DO I = 1,NEQ
	CML(I) = CABFX(I,NUME)*CBFOC(NUME)
	ENDDO
	

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE PRNOBJ(X556,OBJ,IGIDMEM,NCABZ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
C	-------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)


	DIMENSION X556(NCABZ)
	DIMENSION IGIDMEM(NCABZ)

	WRITE(6,140)
	WRITE(6,101) OBJ
	WRITE(6,110) 

	WRITE(ISO,140)
	WRITE(ISO,101) OBJ
	WRITE(ISO,110) 

	WRITE(100,140)
	WRITE(100,101) OBJ
	WRITE(100,110) 

	DO 100 ICZ = 1,NCABZ

	MEM = IGIDMEM(ICZ)

	WRITE(6,120) MEM,X556(ICZ)
	WRITE(ISO,120) MEM,X556(ICZ)
	WRITE(100,120) MEM,X556(ICZ)

100	CONTINUE
	WRITE(6,130) 
	WRITE(ISO,130) 
	WRITE(100,130) 



101	FORMAT(15X,'OPTIMIZED STRAIN ENERGY=',E15.6//)

110	FORMAT(15X,'OPTIMIZED CABLE PRE-TENSION FORCES'//,
	1       15X,'CABLE NO.            PRE-TENSION FORCE')
120	FORMAT(15X,I5,15X,E15.6)
130	FORMAT(' ')

140	FORMAT ('1',//,23X,22('-'),/,24X,'OPTIMIZATION RESULTS',/,23X,
     1		22('-')//)     

	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MOMCON1(IGIDMEM,FIN,NEF,NELE,PCMOM,CNMOM,MEL,KEG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	COMMON /CB557I/ MCSUM,MCFOC,MCDEF,MCMOM,LS557,LF557,LD557,LM557,
	1			    LR557,NM556

	DIMENSION PCMOM(MCMOM,4),IGIDMEM(NELE),FIN(NEF)
	DIMENSION CNMOM(2*MCMOM)


	DO 100 ICMOM = 1,MCMOM

	IGP = PCMOM(ICMOM,1)
	IGD = PCMOM(ICMOM,2)
	IDF = PCMOM(ICMOM,3)
	VAL = PCMOM(ICMOM,4)
	
	IF(KEG.EQ.IGP) THEN
	CALL EGIDNUM(IGIDMEM,IGD,NELE,MEM)
	IF(MEM.EQ.MEL) THEN
	CNMOM(ICMOM) = VAL
	CNMOM(ICMOM) = -1.0*(CNMOM(ICMOM) + FIN(IDF))
	ENDIF
	ENDIF

	

100	CONTINUE

	DO 200 ICMOM = 1,MCMOM

	IGP = PCMOM(ICMOM,1)
	IGD = PCMOM(ICMOM,2)
	IDF = PCMOM(ICMOM,3)
	VAL = PCMOM(ICMOM,4)
	
	IF(KEG.EQ.IGP) THEN
	CALL EGIDNUM(IGIDMEM,IGD,NELE,MEM)
	IF(MEM.EQ.MEL) THEN
	CNMOM(ICMOM+MCMOM) = -1.0*VAL
	CNMOM(ICMOM+MCMOM) = (CNMOM(ICMOM+MCMOM) + FIN(IDF))
	ENDIF
	ENDIF

200	CONTINUE
	

	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MOMCON2(G556,CNMOM,NCOBZ,MCMOM,MCSUM,MCDEF)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------

	DIMENSION G556(NCOBZ),CNMOM(2*MCMOM)

C	WRITE(*,*) CNMOM(1),CNMOM(2)

	KAK = MCSUM + 2*MCDEF 

	DO 100 ICMOM = 1,2*MCMOM

	G556(ICMOM + KAK) = CNMOM(ICMOM)

100	CONTINUE


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE PRNCOZ(G556,PCSUM,PCDEF,PCMOM,NCABZ,
	1				  NCOBZ,MCSUM,MCDEF,MCMOM)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

	DIMENSION G556(NCOBZ)
	DIMENSION PCSUM(MCSUM,NCABZ+3+20),NOM(20)
	DIMENSION PCDEF(MCDEF,3)
	DIMENSION PCMOM(MCMOM,4)

	ICN = 0
	DO IGR = 1,MCSUM
	ICN = ICN + 1
	NODE = PCSUM(IGR,3)
	DO KK = 1,NODE
	NOM(KK) = PCSUM(IGR,KK+3)
	ENDDO
	WRITE(6,110) IGR
	WRITE(6,111) (NOM(II),II=1,NODE)
	WRITE(6,112)  G556(ICN)
	WRITE(ISO,110) IGR
	WRITE(ISO,111) (NOM(II),II=1,NODE)
	WRITE(ISO,112)  G556(ICN)
	WRITE(100,110) IGR
	WRITE(100,111) (NOM(II),II=1,NODE)
	WRITE(100,112)  G556(ICN)
	ENDDO


	DO II = 1,MCDEF
	ICN = ICN + 1
	NODE = PCDEF(II,1)
	MDOF = PCDEF(II,2)
	VALV = PCDEF(II,3)
	IF(VALV.EQ.0.0) VALV = 1.0
	NEGA = ABS(G556(ICN)/VALV*100.0)
	POSI = ABS(G556(ICN+MCDEF)/VALV*100.0)
	IF(G556(ICN).LE.0.0) NEGA = 0.0D0
	IF(G556(ICN+MCDEF).LE.0.0) POSI = 0.0D0
	WRITE(6,116) NODE,MDOF
	WRITE(6,117) POSI,NEGA
	WRITE(ISO,116) NODE,MDOF
	WRITE(ISO,117) POSI,NEGA
	WRITE(100,116) NODE,MDOF
	WRITE(100,117) POSI,NEGA
	ENDDO


	ICN = ICN + MCDEF

	DO II = 1,MCMOM
	ICN = ICN + 1
	IGR  = PCMOM(II,1)
	IEL  = PCMOM(II,2)
	MDOF = PCMOM(II,3)
	VALV = PCMOM(II,4) 
	IF(VALV.EQ.0.0) VALV = 1.0
	NEGA = ABS(G556(ICN)/VALV*100.0)
	POSI = ABS(G556(ICN+MCMOM)/VALV*100.0)
	IF(G556(ICN).LE.0.0) NEGA = 0.0D0
	IF(G556(ICN+MCMOM).LE.0.0) POSI = 0.0D0
	WRITE(6,118) IGR,IEL,MDOF
	WRITE(6,119) POSI,NEGA
	WRITE(ISO,118) IGR,IEL,MDOF
	WRITE(ISO,119) POSI,NEGA
	WRITE(100,118) IGR,IEL,MDOF
	WRITE(100,119) POSI,NEGA
	ENDDO



110	FORMAT(//'HORIZONTAL FORCES CONSTRAINT GROUP',I5)
111	FORMAT(  'NODE NUMBER IN THIS GROUP. . . . .',20I6)
112	FORMAT(  'CONSTRAINT VIOLATED VALUE. . . . .',E15.6)

116	FORMAT(//'DISPLACEMENT CONSTRAINT AT NODE',I5,'  DOF',I5)
117	FORMAT(  'PERCENT CONSTRAINT VIOLATED. . .POSITIVE',E15.6,
	1		/'                                NEGATIVE',E15.6)

118	FORMAT(//'MEMBER FORCE CONSTRAINT OF ELEMENT GROUP',I5,
	1		 '  ELEM NO.',I5'  DOF',I5)
119	FORMAT(  'PERCENT CONSTRAINT VIOLATED. . .POSITIVE',E15.6,
	1		/'                                NEGATIVE',E15.6)


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MOMGRP(IGIDMEM,PCMOM,MEL,KEG,NELE)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	COMMON /CB557I/ MCSUM,MCFOC,MCDEF,MCMOM,LS557,LF557,LD557,LM557,
	1			    LR557,NM556

	DIMENSION PCMOM(MCMOM,4),IGIDMEM(NELE)

	IGD = IGIDMEM(MEL)

	DO ICMOM = 1,MCMOM
	IEM = PCMOM(ICMOM,2)
	IF(IEM.EQ.IGD) PCMOM(ICMOM,1) = KEG
	ENDDO


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE RECVAR(X,CAFOC,MGRUP,VUB,VLB,IGIDM,NMVAR,NCABZ,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------
	DIMENSION CAFOC(1),IGIDM(1),VUB(1),VLB(1),VUBD(NCABZ),VLBD(NCABZ)
	DIMENSION X(1),MGRUP(NMVAR,1)

	SELECT CASE(IND)

	CASE(0)
	VUBD(1:NCABZ) = VUB(1:NCABZ)
	VLBD(1:NCABZ) = VLB(1:NCABZ)
	VUB(1:NCABZ) = 0.0D0
	VLB(1:NCABZ) = 0.0D0
	DO I = 1,NMVAR
	N = MGRUP(I,1)   !ONLY 1st Element in group
	CALL EGIDNUM(IGIDM,N,NCABZ,IEL)
	X(I) = CAFOC(IEL)
	VUB(I) = VUBD(IEL)
	VLB(I) = VLBD(IEL)
	ENDDO

	CASE(1)
	DO I = 1,NMVAR
	DO J = 1,NCABZ
	N = MGRUP(I,J)
	IF(N.NE.0) CALL EGIDNUM(IGIDM,N,NCABZ,IEL)
	IF(N.NE.0) CAFOC(IEL) = X(I)
	ENDDO
	ENDDO


	ENDSELECT


C	WRITE(*,*) MGRUP(1,1:NCABZ)
CC	WRITE(*,*) MGRUP(2,1:NCABZ)
C	WRITE(*,*) MGRUP(3,1:NCABZ)
C	WRITE(*,*) MGRUP(4,1:NCABZ)
C	PAUSE


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================

	SUBROUTINE OBJECZ_OOLD(B,A,ESTN,NEQ,MAXA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	=============================================================
	DIMENSION A(1),B(NEQ),R(NEQ)
	DIMENSION MAXA(NEQ+1)

	
C	DO I = 1,NEQ
C	KHIGH = MAXA(I+1) - MAXA(I)
C	LHEI(I) = I - KHIGH + 1
C	ENDDO
	ESTN = 0.0D0

C	=====================================================
	DO IEQ = 1,NEQ
	R(IEQ) = 0.0D0
C	-----------------------------------------------------
	DO JEQ = IEQ,NEQ	!UPPER TRIANGULAR
	NUM = MAXA(JEQ) + (JEQ - IEQ)
	KHIGH = MAXA(JEQ+1) - MAXA(JEQ)
	LHEI = JEQ - KHIGH + 1
	IF(LHEI.LE.IEQ) R(IEQ) = R(IEQ) + A(NUM)*B(JEQ)
	ENDDO
C	-----------------------------------------------------
	NUM = MAXA(IEQ+1) - MAXA(IEQ)
	LCL = MAXA(IEQ)-1
	DO J = 2,NUM		!LOWER TRIANGULAR
	JEQ = IEQ - J + 1
	IF(JEQ.GT.0) R(IEQ) = R(IEQ) + A(LCL+J)*B(JEQ)
	ENDDO
C	-----------------------------------------------------
	ESTN = ESTN + 0.5*R(IEQ)*B(IEQ)
C	-----------------------------------------------------
	ENDDO
C	=====================================================


	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================