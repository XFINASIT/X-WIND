C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE CNCABLE(PROPM,PROPG,NODEX,SPCFR,S,COORD,EDIS,EDISI,RE,
     #                   MWG,FIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------------
C     3D PARABOLIC CABLE  (REFERENCE FROM SPC-FRAME)
C     ---------------------------------------------------------------------
C
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV


      DIMENSION PROPM(1),PROPG(5),EDIS(6),EDISI(6)
      DIMENSION COORD(6),COD(6)
C
	DIMENSION HN(2,6),DISP(2,1),RE(6)
	DIMENSION S(21)


	DIMENSION STIF(6,6)
	DIMENSION FIN(1)
	DIMENSION SPCFR(1)

	DIMENSION SKG(6,6),SKR(6,6)

	DIMENSION FOC(4),FO1(4),FO2(4)
	DIMENSION VL(3),VH(3),VG(3),VT(3)
	DIMENSION TRANS(6,6),FCO(6)
	DIMENSION CSTIF(6,6)


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
	YOUNG = PROPM(1)

C	---------------------------------------
C	SET VALUES FOR LINEAR STRESS-STRAIN LAW 
C	INITIALISATION OF INTEGRATION RULE
C	---------------------------------------
C	
	AREA  = PROPG(2)

	IFOPT = INT(PROPG(4))  !IFOPT    1=INITIAL FORCE WILL ADD TO LOAD VECTOR    0=INITIAL FORCE WILL NOT ADD TO LOAD VECTOR FEB09
	IFLAG = INT(SPCFR(8))  !FLAG TO CHECK WHETHER THE INITIAL FORCES IS ALREADY STORED OR NOT

	TDEN  = PROPM(5)  !DENSITY
	WDEN  = PROPM(6)  !UNIT WEIGHT

	PDEN  = WDEN*AREA !WDEN*AREA  !WEIGHT PER UNIT LENGTH

	DX = COORD(4) - COORD(1)
	DY = COORD(5) - COORD(2)
	DZ = COORD(6) - COORD(3)

	TOL = 1.0E-8
	IF(ABS(DX).EQ.0.0D0) DX = TOL
	IF(ABS(DY).EQ.0.0D0) DY = TOL
	IF(ABS(DZ).EQ.0.0D0) DZ = TOL
	IF(ABS(DX).LT.TOL) DX = TOL*DX/ABS(DX)
	IF(ABS(DY).LT.TOL) DY = TOL*DY/ABS(DY)
	IF(ABS(DZ).LT.TOL) DZ = TOL*DZ/ABS(DZ)

	ELN		= SQRT(DX*DX + DY*DY + DZ*DZ)


C	------------------------------------------------
C	FORM MASS MATRIX
C	------------------------------------------------
	IF (ITASK.NE.5) GOTO 50
	SSAM = TDEN*AREA*ELN
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
C	------------------------------------------------
50	CONTINUE 
C	------------------------------------------------

	VL(1:3) = [DX,DY,DZ]

	VH = VL
	VH(NGRAV) = 0.0D0
	CALL SCALEN (VH,VH,DUM,3)
	VG(1:3)   = 0.0D0
	VG(NGRAV) = 1.0D0

	CALL VECPRD (VH,VG,VT)

	CALL SCAPRD (VL,VH,Z,3)
	CALL SCAPRD (VL,VG,T,3)

	T = COORD(NGRAV+3) - COORD(NGRAV)

	TRANS = 0.0D0
	TRANS(1,1:3) = VH(1:3)
	TRANS(2,1:3) = VG(1:3)
	TRANS(3,1:3) = VT(1:3)
	TRANS(4,4:6) = VH(1:3)
	TRANS(5,4:6) = VG(1:3)
	TRANS(6,4:6) = VT(1:3)

C	CLNO  = PROPG(3)
	FHOR  = PROPG(3) !INITIAL HORIZONTAL FORCE
	CLNO  = SPCFR(9) !THIS ARRAY STORE THE INITIAL LENGTH
	IF(CLNO.EQ.0.0) THEN
	CALL CNINLEN(ELN,Z,T,AREA,YOUNG,PDEN,FHOR,CLNO)
	SPCFR(9) = CLNO
	ENDIF


	CALL PCAFX2 (Z,T,AREA,YOUNG,CLNO,PDEN,FOC)


	ALPHA = 1.0E-6
	H = Z+ALPHA
	V = T
	CALL PCAFX2 (H,V,AREA,YOUNG,CLNO,PDEN,FO1)

	V = T+ALPHA
	H = Z
	CALL PCAFX2 (H,V,AREA,YOUNG,CLNO,PDEN,FO2)

	A1 = (FO1(3)-FOC(3))/(ALPHA)
	A2 = (FO1(4)-FOC(4))/(ALPHA)
	A3 = (FO2(3)-FOC(3))/(ALPHA)
	A4 = (FO2(4)-FOC(4))/(ALPHA)

	FCO = 0.0D0
	FCO(1:2) = FOC(1:2)
	FCO(4:5) = FOC(3:4)
	RE = MATMUL(TRANSPOSE(TRANS),FCO)
	IF(IFLAG.EQ.0) THEN
	SPCFR(2:7) = RE(1:6)
	SPCFR(8) = 1.0       !UPDATE FLAG TO CHECK WHETHER THE INITIAL FORCES IS ALREADY STORED OR NOT  1=ALREADY STORED
	ENDIF
	IF(IFOPT.EQ.0) RE(1:6) = RE(1:6) - SPCFR(2:7)


	CSTIF = 0.0D0

	CSTIF(1,1) = A1
	CSTIF(2,1) = A2
	CSTIF(1,2) = A3
	CSTIF(2,2) = A4

	CSTIF(4,4) = A1
	CSTIF(5,4) = A2
	CSTIF(4,5) = A3
	CSTIF(5,5) = A4

	CSTIF(1,4) = -A1
	CSTIF(2,4) = -A2
	CSTIF(1,5) = -A3
	CSTIF(2,5) = -A4

	CSTIF(4,1) = -A1
	CSTIF(5,1) = -A2
	CSTIF(4,2) = -A3
	CSTIF(5,2) = -A4


	TOL = 1.0e-6
	SMAX = 0.0D0
	DO I = 1,6
	IF(ABS(CSTIF(I,I)).GT.SMAX) SMAX = ABS(CSTIF(I,I))
	ENDDO
	DO I = 1,6
	IF(ABS(CSTIF(I,I)).LT.TOL) CSTIF(I,I) = CSTIF(I,I) + SMAX
	ENDDO

	SKR = MATMUL(TRANSPOSE(TRANS),MATMUL(CSTIF,TRANS))

	KK = 1
	DO I = 1,6
	FIN(I) = RE(I)
	DO J = I,6
	S(KK) = SKR(I,J)
	KK = KK+1
	ENDDO
	ENDDO

C	STORE THE HORIZONTAL FORCE OF CABLE
	SPCFR(1) = ABS(FOC(1))


	RETURN


	END
C
C	=================================================================
C	=================================================================
C	=================================================================	
	SUBROUTINE PCAFX2 (Z,T,A,E,XLO,WO,FOC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION FOC(4)
C	=====================================================================


	EPS=.00001

      PI=3.14159265
	H=Z
	V=T
	ICODE=3

	TEMP = 0.0D0
	ET   = 0.0D0	
C      
C     ADJUST UNSTRETCHED LENGTH AND LOAD FOR TEMPERATURE
C
	XL=XLO*(1.+ET*TEMP)
	XLIN=XL
	W=WO*XLO/XL

C	
C	INTERCHANGE ORIGIN AND END OF CABLE IF V IS POSITIVE
C
	KK = 0
	IF(V.LE.0.0D0) GOTO 2
	KK = 1
	V = -V
	H = -H
C
C	INITIALIZE LAMDA EQ 17 OR 1000000 OR .20
C
2	CORD = SQRT(H*H+V*V)
	AMBDA = 1000000.0
	IF(H.EQ.0.0D0) GOTO 4
	IF(XL.LE.CORD) GOTO 3
	AMBDA = SQRT(((XL*XL-V*V)/(H*H)-1.0D0)*3.0D0)
	GOTO 4
3	AMBDA = 0.20
C
C	STARTING VALUES OF F1 AND F2 EQ 16 AND 18
C	
4	FO1 = -(W*H)/(2.*AMBDA)
	COT = 1./TANH(AMBDA)
	FO2 = (W/2.)*(-V*COT+XL)
	DF1 = 0.0D0
	DF2 = 0.0D0
C	
C	APPLY CORRECTIONS TO F1 AND F2 EQ 9
C	
5	FO1 = FO1 + DF1
	FO2 = FO2 + DF2
6	FO4 = W*XL - FO2
	FO3 = -FO1
	TI = SQRT(FO1*FO1+FO2*FO2)
	TJ = SQRT(FO3*FO3+FO4*FO4)
	F = FO4+TJ
	FF = TI - FO2
	IF(FF.LT.0.0001) FF = 0.0001
	G = F/FF
	IF(G.LT.0.0001) G = 0.0001
C
C	COMPUTE VALUES OF H AND V EQ 4 AND 5
C
	AAH = (1./W)*LOG(G)+XL/(A*E)
	AH  = -FO1*AAH
	BV  = (TJ*TJ-TI*TI)/(2.*E*A*W)+(TJ-TI)/W

C
C	CALCULATE MISCLOSURE VECTOR AND CHECK CONVERGENCE
C
	CA = H - AH
	CB = V - BV
	ACA = ABS(CA)
	ACB = ABS(CB)


	IF(ACA.LE.EPS.AND.ACB.LE.EPS) GOTO 10
C
C	DETERMINE CORRECTION TERMS EQ 9 AND 11
C
	IF(TJ.LT.0.0001) TJ = 0.0001
	B2 = -(FO2/TI+FO4/TJ)/W - XL/(E*A)
	A1 = -AAH-B2-XL/(E*A)
	A2 = (FO1/W)*(1./TJ-1./TI)
	DET = A1*B2 - A2*A2
	DF1 = (CA*B2-CB*A2)/DET
	DF2 = (A1*CB-A2*CA)/DET
	GOTO 5
C
C	ONCE CONVERGENCE ACHIEVED DETERMINE STRETCHED LENGTH AND
C	ACTUAL END FORCES
C	
10	XLAFST = XL + (FO4*TJ+FO2*TI+FO1*FO1*LOG(G))/(2.*E*A*W)
	FOC(1) = FO1*(1.-2.*FLOAT(KK))
	FOC(2) = FO2 + FLOAT(KK)*(FO4-FO2)
	FOC(3) = FO3*(1.-2.*FLOAT(KK))
	FOC(4) = FO4 + FLOAT(KK)*(FO2-FO4)

18	RETURN
	END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE CNINLEN(ELN,Z,T,AREA,YOUNG,PDEN,FHOR,CLNO)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------------
C     CALCULATE THE INITIAL LENGTH FROM INITIAL HORIZONTAL PRETENSION OF CATENARY CABLE
C     ---------------------------------------------------------------------

      DIMENSION FOC(4)

	NTER = 500
	TOL  = 1.0E-5
C	TRIAL INITIAL LENGTH CL
	CLL = 0.5*ELN
	CLU = 5.0*ELN
	CLO = CLL
	DO ITER = 1,NTER
	CL = 0.5*(CLL + CLU)
	CALL PCAFX2 (Z,T,AREA,YOUNG,CL,PDEN,FOC)

C	GET HORIZONTAL FORCE
	FI = ABS(FOC(1)) !SQRT(FOC(1)*FOC(1) + FOC(2)*FOC(2))
	FJ = ABS(FOC(3)) !SQRT(FOC(3)*FOC(3) + FOC(4)*FOC(4))

C	GET THE FORCE FROM THE MAXIMUM JACKING END
	FF = FI
	IF(FJ.GT.FI) FF = FJ

C	ADJUST THE CABLE LENGTH
	IF(FF.GT.FHOR) THEN		!TOO TIGHT
C	INCREASE CABLE LENGTH
		CLL = CL
	ELSE					!TOO LOOSE
C	DECREASE CABLE LENGTH
		CLU = CL
	ENDIF
	
	TEST = ABS((CL-CLO)*100.0/CLO)
	IF(TEST.LT.TOL) GOTO 100
	CLO = CL

	ENDDO

100	CONTINUE
	CLNO = CL
C	WRITE(*,*) CL,CLNO,FF
C	PAUSE



	RETURN


	END
C
C	=================================================================
C	=================================================================
C	=================================================================
