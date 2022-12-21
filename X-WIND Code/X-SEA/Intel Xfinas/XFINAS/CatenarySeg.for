C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE CNCABLE_NEW_MOORING(PROPM,PROPG,NODEX,SPCFR,S,COORD,EDIS,EDISI,RE,
     #                   MWG,FIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV

C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV
      
      COMMON /OUT_MOOR/ OUTPUT_DATA(2)
      
C     SAVE ELEMENT NUMBER      
      COMMON /NELEM/ IEL   


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
      
      DIMENSION COOR_REF(3)


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

	IFOPT = INT(PROPG(6))  !IFOPT    1=INITIAL FORCE WILL ADD TO LOAD VECTOR    0=INITIAL FORCE WILL NOT ADD TO LOAD VECTOR FEB09
	IFLAG = INT(SPCFR(8))  !FLAG TO CHECK WHETHER THE INITIAL FORCES IS ALREADY STORED OR NOT
      IOUNS = INT(PROPG(4))
      
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
	!FHOR  = PROPG(3) !INITIAL HORIZONTAL FORCE
	CLNO  = SPCFR(9) !THIS ARRAY STORE THE INITIAL LENGTH
      
      !AREA  = 0.00636173D0
      CB        = PROPG(3) 
      IF (IOUNS.EQ.1D0)THEN
      UNSTR     = PROPG(5)
      ELSEIF (IOUNS.EQ.0D0)THEN
      CALL CABLE_LENGTH_ALL_TYPE (UN,UN,UN,UNSTR_OUT,IEL,"FIND",COOR_REF)
      UNSTR     = UNSTR_OUT
      ENDIF
      GRAV      = 9.806D0
      ROLLW     = 1025D0
      RHO       = PROPM(5)  !DENSITY
      
      OMEGA     = GRAV*AREA*(RHO-ROLLW)
      CALL SEABEDCONTACT (FAIRLEADH,FAIRLEADV,ROLLW,RHO,OMEGA,GRAV,CB,AREA,YOUNG,UNSTR,Z,T,XS,YS,STIFFNESS,FOC_OUT,COOR_REF,"INTER")  
      FHOR = SQRT(OUTPUT_DATA(1)**2D0 + OUTPUT_DATA(2)**2D0)
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
      
      
      ! SEABED CONTRACT CHECK
      CALL SEABEDCONTACT (FAIRLEADH,FAIRLEADV,ROLLW,RHO,OMEGA,GRAV,CB,AREA,YOUNG,UNSTR,Z,T,XS,YS,STIFFNESS,IOPSTFF,COOR_REF,"CALC") 
      IF (NGRAV.EQ.1)THEN
        IF (COORD(1).LE.COOR_REF(1)) THEN
           A1 = A1*1D0*10**10
           A2 = A2*1D0*10**10
          ENDIF
          IF (COORD(4).LE.COOR_REF(3)) THEN
           A3 = A1*1D0*10**10
           A4 = A2*1D0*10**10
          ENDIF
      ELSEIF (NGRAV.EQ.2)THEN
       IF (COORD(2).LE.COOR_REF(2)) THEN
           A1 = A1*1D0*10**10
           A2 = A2*1D0*10**10
          ENDIF
          IF (COORD(4).LE.COOR_REF(2)) THEN
           A3 = A1*1D0*10**10
           A4 = A2*1D0*10**10
          ENDIF
      ELSEIF (NGRAV.EQ.3)THEN
          IF (COORD(3).LE.COOR_REF(3)) THEN
           A1 = A1*1D0*10**10
           !A2 = A2*1D0*10**10
          ENDIF
          IF (COORD(6).LE.COOR_REF(3)) THEN
           !A3 = A1*1D0*10**10
           A4 = A2*1D0*10**10
          ENDIF
      ENDIF


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
C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE CNCABLE_MOORING (COORD,EDIS,PROPM,PROPG,RHO,GAMMA,WX,WY,WZ,FHOR,NSEG,NPTS,S,RE,FIN,DISLI,SPCFR)
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


      DIMENSION PROPM(1),PROPG(6),EDIS(6),EDISI(6)
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
      DIMENSION STIFFNESS(4)


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
      CB    = PROPG(3) 
      IF (PROPG(4).EQ.1) CLNO = PROPG(5)

	IFOPT = INT(PROPG(6))  !IFOPT    1=INITIAL FORCE WILL ADD TO LOAD VECTOR    0=INITIAL FORCE WILL NOT ADD TO LOAD VECTOR FEB09
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

      RHO       = TDEN   ! CABLE DENSITY    
      
      !OMEGA     = WEIGHT PER UNIT LENGTH IN THE SUBMERGRD FLUID   
      GRAV      = 9.806D0
      ROLLW     = 1025D0
      OMEGA     = GRAV*AREA*(RHO-ROLLW)
      
      ! CALCULATED FORCE ON FAIR LEAD AT H AND V
      CALL SEABEDCONTACT (FAIRLEADH,FAIRLEADV,ROLLW,RHO,OMEGA,GRAV,CB,AREA,YOUNG,CLNO,Z,T,XS,YS,STIFFNESS,FOC,"INTER")  
      !STIIFFNESS = 100000000000000000.
	A1 = STIFFNESS(1)
	A2 = STIFFNESS(2)
	A3 = STIFFNESS(3)
	A4 = STIFFNESS(4)
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
	SUBROUTINE PCAFX2_MOORING (Z,T,A,E,XLO,WO,W,FOC)
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
6     FO4 = W*XL - FO2
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
	SUBROUTINE CNINLEN_MOORING (ELN,Z,T,AREA,YOUNG,WC,FHOR,CLNO)
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
	CALL PCAFX2 (Z,T,AREA,YOUNG,CL,WC,FOC)

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


	RETURN


	END
C
C	=================================================================
C	=================================================================
C	=================================================================

C	=================================================================
C	=================================================================
C	=================================================================	
	SUBROUTINE CATSEG (NSEG,NPTS,Z,T,A,E,XLO,W,FOC,XCOOR)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION FOC(4),XCOOR(2,NPTS)!XCOOR(2,NSEG)

C	=====================================================================
C     Calculate the coordinate along the cable w.r.t calble 2D plane

          
      PI=3.14159265
	H=Z
	V=T
      
	
	KK = 0
	IF(V.LE.0.0D0) GOTO 2
	KK = 1

2     SUBXL = XLO/(FLOAT(NPTS)-1.0)
      XL = -SUBXL
      
      FO1 = FOC(1) - (FOC(1)-FOC(3))*FLOAT(KK)
      FO2 = FOC(2) - (FOC(2)-FOC(4))*FLOAT(KK)
      
      DO MM = 1,NPTS
      !DO MM = 1,NSEG
      XL = XL+SUBXL
      FO4 = W*XL - FO2
	FO3 = -FO1
	TI = SQRT(FO1*FO1+FO2*FO2)
	TJ = SQRT(FO3*FO3+FO4*FO4)
	F = FO4+TJ
	FF = TI - FO2
	IF(FF.LT.0.0001) FF = 0.0001
	G = F/FF
	IF(G.LT.0.0001) G = 0.0001
C
C	COMPUTE VALUES OF H AND V (Horizontal & Vertical length)
C
	AAH = (1./W)*LOG(G)+XL/(A*E)
	AH  = -FO1*AAH
	BV  = (TJ*TJ-TI*TI)/(2.*E*A*W)+(TJ-TI)/W
	!!NN = MM + (NPTS - 2*MM + 1)*KK
	!!XCOOR(1,NN) = AH + Z*FLOAT(KK)
	!!XCOOR(2,NN) = BV + T*FLOAT(KK)
      
	XCOOR(1,MM) = AH !+ Z*FLOAT(KK)
	XCOOR(2,MM) = BV !+ T*FLOAT(KK)      
            	
      ENDDO
      
      !OPEN(UNIT=10005  ,FILE='CATENARY.dat'       ,STATUS='UNKNOWN'    )
     
      !ITI = 10005
      !WRITE(ITI,*) XCOOR(1:2,1:NSEG)    
      

      RETURN
	END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE CATITER (NSEG,NPTS,FCO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION FCO(6)
C	=====================================================================
C     STATIC EQUILIBRIUM BY BJ 
C       KNOWNS : POSITION OF CONNECT NODE
C       UNKNOWNS 
C        : FX1, FY1, FAX1, FAY1, AND Q1 <= Local unknown in line 1
C        : FX2, FY2, FAX2, FAY2, AND Q2 <= Local unknown in line 2
C        : FX3, FY3, FAX3, FAY3, AND Q3 <= Local unknown in line 3
        
C         Qi = (Ri)T(Ri+1 - Ri)
      
         !1) Solve FAXi & FAYi from FXi & FYi
         !2) Substitute FAXi & FAYi into Equilibrium Equation       
       
       DO KK = 1,NPTS
       
       !FAXi + FAXi
           
       !FH1X + FHA2X + FHA3X = 0.0D0
       !FH1Y + FHA2Y + FHA3Y = 0.0D0    
           
       !FEXT = FV1 + FVA2 + FVA3    
       
       ENDDO

      RETURN
	END


C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SEABEDCONTACT (FAIRLEADH,FAIRLEADV,ROLLW,RHO,OMEGA,GRAV,CB,AREA,YOUNG,CLNO,Z,T,XS,YS,STIFFNESS,FOC,COOR_REF,LTYPE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)   
      CHARACTER*4 LTYPE
C     -----------------------------------------------------------------------------------------------      
C     --- INPUT DATA ---
C     FAIRLEADH = HORIZONTAL FAIRLEAD FORCE 
C     FAIRLEADV = VERTICAL   FAIRLEAD FORCE 
C     ROLLW     = FLUID DENSITY
C     RHO       = Mass density of cable material      
C     GRAV      = GRAVITY ACCLECTION 
C     CB        = Mooring Line Static friction Coeff
C     AREA      = CROSS-SECTIONAL AREA
C     YOUNG     = YOUNG MODULUS
C     CLNO      = UNSTRCTCHED LENGTH
C     T         = VECTOR
      
      
C     --- OUTPUT DATA ---  
C     XS        = CABLE PROFILE IN THE HORIZONTAL DIRECTION
C     YS        = CABLE PROFILE IN THE VERTICAL DIRECTION
       
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C	GRAVTITY DIRECTION ADDED BY SONGSAK MAR2006  
	COMMON /MGRAV/ NGRAV
C     SAVE ELEMENT NUMBER      
      COMMON /NELEM/ IEL           
C     MOORING LINE BY BJ
      COMMON /MOOR/ NMOOR,MOORTASK,MOORNODE(10000) ! SAVE ANCHOR NODE IN EACH ELEMENT          
      COMMON /OUT_MOOR/ OUTPUT_DATA(2)
      
      DIMENSION FOC(4),IOPSTFF(4)
      DIMENSION STIFFNESS(4)
      DIMENSION COOR_REF(3)
C     -------------------------------------------------------------
      ! STARING CATENARY EQUATION
      PI=3.14159265 
          
      W = OMEGA
      E = YOUNG
      A = AREA      
      
      !FAIRLEADH = 30000.0D0
      !FAIRLEADV = 30000.0D0
      !ROLLW     = 1000D0
      !GRAV      = 9.806D0
      !CB        = 1.0D0
C     CLNO      = 10D0
    

C     ALENGTH   = IN TERM OF LENGTH (ASSUME EQ. CLNO)
      !ALENGTH   = CLNO
      EA        = AREA*YOUNG  
      ALC       = CLNO
      IOPSTFF   = 0.

      
      ! FROM EQUATION 3.9B SOLVE FOR H (BY USING MATCAD = 0 SOLVE,H)
      ! ---- UNKNOW ----
      ! FAIRLEADV = ?
      ! FAIRLEADH = ?
      ! --- EXAMPLE INPUT DATA ---
      !CB     = 0.001D0 ! SEA BED FRICTION COEFF.
      
      
      ! STATING CALCULATION OF TRUE LENGTH
      ! SEND      IEL = NUMBER OF ELEMENT
      ! RECEIVE   COOR_H1,COOR_H2,COOR_V
      !CALL CABLE_LENGTH (ALC,IEL,ALENGTH,COOR_H,COOR_V,"CALC") ! FOR FRAME ELEMENT
      
      ! SEND      IEL = NUMBER OF ELEMENT IN GROUP
      ! RECEIVE   ALC = UNSTRECTCHED LINE LENGTH
      !           ALENGTH
      !CALL CABLE_LENGTH (ALC,IEL,ALENGTH,COOR_H,COOR_V,"READ")
      ALENGTH = 0.0D0
      ALB = 0.0D0
      COOR_REF = 0.0D0
      
      CALL CABLE_LENGTH_ALL_TYPE (ALB,COOR_H,COOR_V,ALENGTH,IEL,"FIND",COOR_REF)
      
      ! COOR_H AND COOR_V IS LOCAL LENGTH
      ! 
      IF (LTYPE.EQ."INTE")THEN
          
      ! REQUESTED POSTION VECTOR
      !COOR_V = 186D0
      !COOR_H = 796.914D0
      !ALC    = 835.35D0
      !AREA   = PI*(0.0766D0**2D0)/4D0
      !YOUNG  = (7.536*10**8)/AREA
      !W      = 1065.6251D0
      !ALC     = CLNO
          
      !EA = 384243000.00D0 
      !ALC = 902.2D0
      !COOR_V = 250D0
      !COOR_H = 848.67D0
          
      ! --------------------------
      FAIRLEADV = W*ALC*0.5D0
      !CON_FAIR  = W*ALC
      !FAIRLEADV = 584448.744D0
      !FAIRLEADH = 894066.806D0
      DO 100 I = 1, 50000
      !DO 200 K = 1, 10000
      
      !BB1 = FAIRLEADH/W
      !BB2 = SQRT(1+(FAIRLEADV/FAIRLEADH)**2)-1
      !BB3 = (FAIRLEADV**2)/(2D0*YOUNG*AREA*W)
      !COOR_V_NEW = (BB1*BB2)+BB3
      !FAIRLEADV = 631704.765D0
          
      BB1 = 4D0*(AREA**2)*(YOUNG**2)*(COOR_V**2)*(W**2)
      BB2 = 4D0*(AREA**2)*(YOUNG**2)*(FAIRLEADV**2)
      BB3 = 4D0*AREA*YOUNG*(FAIRLEADV**2)
      BB4 = (FAIRLEADV**4)
      BB5 = 4D0*AREA*YOUNG*(FAIRLEADV**2)
      BB6 = 8D0*(AREA**2)*(YOUNG**2)*COOR_V*W
      FAIRLEADH  = (BB1-BB2-BB3+BB4)/(BB5-BB6) 
      
      ALB   = ALC - FAIRLEADV/W
      
      IF (ALB.LE.0)THEN
           FAIRLEADV = -(ALB-ALC)*W
      ENDIF
      
      ! FROM EQUATION 3.9A SOLVE FOR V 
      AA1 = ALB
      AA2 = (FAIRLEADH/W)*(LOG((FAIRLEADV/FAIRLEADH)+SQRT(1D0+((FAIRLEADV/FAIRLEADH)**2))))
      AA3 = FAIRLEADH*ALC/(YOUNG*AREA)
      AA4 = CB*W/(2D0*(YOUNG*AREA))
      AA5 = ALC-(FAIRLEADV/W)-(FAIRLEADH/(CB*W))
      IF (AA5.GT.0)THEN
      AA6 = ALC-(FAIRLEADV/W)-(FAIRLEADH/(CB*W))
      ELSE
      AA6 = 0D0
      ENDIF
      AA7 = (ALC-(FAIRLEADV/W))**2D0
      
      AA8 = COOR_H
      AA  = AA1+AA2+AA3+AA4*(AA6*AA5-AA7)
      RATIO = abs(((AA8-AA)/AA8)*100D0)
      
      IF (RATIO.GT.0.01D0)THEN ! 
      ! FIND NEW FAIRLEADH AND FAIRLEADV
          IF (RATIO.GT.0.01D0) FAIRLEADV = FAIRLEADV*(RATIO/200) + FAIRLEADV
          IF (RATIO.LE.0.01D0) FAIRLEADV = FAIRLEADV*(RATIO/2000) + FAIRLEADV
      ELSEIF (RATIO.LE.0.01D0)THEN
      ! FINISHED OF CALCULATION
      EXIT
      ENDIF
100   CONTINUE
      
      OUTPUT_DATA(1) = FAIRLEADH
      OUTPUT_DATA(2) = FAIRLEADV

      
      ! POSITION OF ALB
      

      !FOC(1) = FAIRLEADV
      !FOC(2) = FAIRLEADH
C     ALB   = UNSTRETCHED LENGTH OF CABLE LYING ON THE SEABED   
      ALB   = ALC - FAIRLEADV/W
C     IF VERTICAL FORCE (FAIRLEADV) IS LESS THAN THE TOTOAL WEIGHT OF THE CABLE
C     , THEN A PORTION OF THE MOORING LINE WILL REST ON THE SEABED
      CON00 = W/ALC

C     STARTINF CALUCULATE MOORING PROFILE WITH SEABEAD CONTRACT
      CON0  = ALB - FAIRLEADH/(CB*W)   
         
         IF (CON0.GT.0.0D0)THEN
            RAMDA = ALB - FAIRLEADH/(CB*OMEGA) 
         ELSEIF (CON0.NE.0.0D0)THEN
            RAMDA = 0.0D0
         ENDIF
         
      CALL MOORING_FORCE (FAIRLEADH,FAIRLEADV,IEL,W,CB,ALB,COOR_H,COOR_V,HA,VA,ALENGTH,"BEG")

!     ****************** BEGIN POSITION ******************
C     HORIZONTAL CABLE PROFILE (TRUE PROFILE) >> BEGIN POSITION
         CON1 = 0.0D0
         CON2 = ALB - HA/(CB*W*(ALB-ALENGTH))
         CON3 = ALB
         CON4 = COOR_H
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                XS1 = ALENGTH
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN   
                CO1 = (ALB-HA/OMEGA)*RAMDA     
                CO2 = 2D0*ALENGTH*(ALB-HA/OMEGA)
                XS1 = ALENGTH+(CB*OMEGA/(2D0*EA))*(ALENGTH**2D0 - CO2 + CO1)
         ELSEIF (ALENGTH.GT.CON3.AND.ALENGTH.LE.CON4)THEN
                CO1 = (CB*OMEGA/(2*EA))*(RAMDA*(ALB-HA/(CB*OMEGA))-ALB**2D0)
                CO2 = HA*ALENGTH/EA
                CO3 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO4 = (OMEGA*(ALENGTH-ALB)/HA) + CO3
                XS1  = ALB + HA/OMEGA*LOG(CO4) + CO2 + CO1
         ENDIF

C     VERTICAL CABLE PROFILE (TRUE PROFILE) >> BEGIN POSITION
         CON1 = 0.0D0
         CON2 = ALB
         CON3 = COOR_H
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                YS1  = 0.0D0
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN
                CO1 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO2 = OMEGA*((ALENGTH-ALB)**2D0)/(2D0*EA)
                YS1  = HA/OMEGA*(CO1 - 1D0) + CO2
         ENDIF

C     CALCULATE TENSION FOR OUTPUT
         
         
       HA = HA +0.0001D0
       VA = VA +0.0001D0
         
         
C     HORIZONTAL CABLE PROFILE (TRUE PROFILE + DT) >> BEGIN POSITION
         CON1 = 0.0D0
         !CON2 = ALB - HA/(CB*W)
         CON2 = ALB - HA/(CB*W*(ALB-ALENGTH))
         CON3 = ALB
         CON4 = ALC
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                XS2 = ALENGTH
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN   
                CO1 = (ALB-HA/OMEGA)*RAMDA     
                CO2 = 2D0*ALENGTH*(ALB-HA/OMEGA)
                XS2 = ALENGTH+(CB*OMEGA/(2D0*EA))*(ALENGTH**2D0 - CO2 + CO1)
         ELSEIF (ALENGTH.GT.CON3.AND.ALENGTH.LE.CON4)THEN
                CO1 = (CB*OMEGA/(2*EA))*(RAMDA*(ALB-HA/(CB*OMEGA))-ALB**2D0)
                CO2 = HA*ALENGTH/EA
                CO3 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO4 = (OMEGA*(ALENGTH-ALB)/HA) + CO3
                XS2  = ALB + HA/OMEGA*LOG(CO4) + CO2 + CO1
         ENDIF

C     VERTICAL CABLE PROFILE (TRUE PROFILE + DT) >> BEGIN POSITION
         CON1 = 0.0D0
         CON2 = ALB
         CON3 = ALC
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                YS2  = 0.0D0
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN
                CO1 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO2 = OMEGA*((ALENGTH-ALB)**2D0)/(2D0*EA)
                YS2  = HA/OMEGA*(CO1 - 1D0) + CO2
         ENDIF
      
      ! BEGIN   
      STIFFNESS(1) = 0.0001D0/(XS2-XS1)
      STIFFNESS(2) = 0.0001D0/(YS2-YS1)
      IF (XS2-XS1.EQ.0.0D0) STIFFNESS(1) = 10**12
      IF (YS2-YS1.EQ.0.0D0) STIFFNESS(2) = 10**12
      FOC(1)       = HA - 0.0001D0
      FOC(2)       = VA - 0.0001D0
      
      
      
!     ****************** END POSITION ******************
      ! SEND      IEL = NUMBER OF ELEMENT IN GROUP
      ! RECEIVE   ALENGTH = LINE LENGTH OF SECOND POSITION
      !CALL CABLE_LENGTH (ALC,IEL,ALENGTH,COOR_H,COOR_V,"REAR")
      
      CALL MOORING_FORCE (FAIRLEADH,FAIRLEADV,IEL,W,CB,ALB,COOR_H,COOR_V,HA,VA,ALENGTH,"END")
 
C     HORIZONTAL CABLE PROFILE (TRUE PROFILE) >> END POSITION
         CON1 = 0.0D0
         !CON2 = ALB - HA/(CB*W)
         CON2 = ALB - HA/(CB*W*(ALB-ALENGTH))
         CON3 = ALB
         CON4 = ALC
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                XS1 = ALENGTH
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN   
                CO1 = (ALB-HA/OMEGA)*RAMDA     
                CO2 = 2D0*ALENGTH*(ALB-HA/OMEGA)
                XS1 = ALENGTH+(CB*OMEGA/(2D0*EA))*(ALENGTH**2D0 - CO2 + CO1)
         ELSEIF (ALENGTH.GT.CON3.AND.ALENGTH.LE.CON4)THEN
                CO1 = (CB*OMEGA/(2*EA))*(RAMDA*(ALB-HA/(CB*OMEGA))-ALB**2D0)
                CO2 = HA*ALENGTH/EA
                CO3 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO4 = (OMEGA*(ALENGTH-ALB)/HA) + CO3
                XS1  = ALB + HA/OMEGA*LOG(CO4) + CO2 + CO1
         ENDIF

C     VERTICAL CABLE PROFILE (TRUE PROFILE) >> END POSITION
         CON1 = 0.0D0
         CON2 = ALB
         CON3 = ALC
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                YS1  = 0.0D0
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN
                CO1 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO2 = OMEGA*((ALENGTH-ALB)**2D0)/(2D0*EA)
                YS1  = HA/OMEGA*(CO1 - 1D0) + CO2
         ENDIF
         
       HA = HA + 0.0001D0
       VA = VA + 0.0001D0
C     HORIZONTAL CABLE PROFILE (TRUE PROFILE + DT) >> END POSITION
         CON1 = 0.0D0
         !CON2 = ALB - HA/(CB*W)
         CON2 = ALB - HA/(CB*W*(ALB-ALENGTH))
         CON3 = ALB
         CON4 = ALC
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                XS2 = ALENGTH
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN   
                CO1 = (ALB-HA/OMEGA)*RAMDA     
                CO2 = 2D0*ALENGTH*(ALB-HA/OMEGA)
                XS2 = ALENGTH+(CB*OMEGA/(2D0*EA))*(ALENGTH**2D0 - CO2 + CO1)
         ELSEIF (ALENGTH.GT.CON3.AND.ALENGTH.LE.CON4)THEN
                CO1 = (CB*OMEGA/(2*EA))*(RAMDA*(ALB-HA/(CB*OMEGA))-ALB**2D0)
                CO2 = HA*ALENGTH/EA
                CO3 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO4 = (OMEGA*(ALENGTH-ALB)/HA) + CO3
                XS2  = ALB + HA/OMEGA*LOG(CO4) + CO2 + CO1
         ENDIF

C     VERTICAL CABLE PROFILE (TRUE PROFILE + DT) >> END POSITION
         CON1 = 0.0D0
         CON2 = ALB
         CON3 = ALC
         IF     (ALENGTH.GE.CON1.AND.ALENGTH.LE.CON2)THEN
                YS2  = 0.0D0
         ELSEIF (ALENGTH.GT.CON2.AND.ALENGTH.LE.CON3)THEN
                CO1 = SQRT(1+(OMEGA*(ALENGTH-ALB)/HA)**2D0) 
                CO2 = OMEGA*((ALENGTH-ALB)**2D0)/(2D0*EA)
                YS2  = HA/OMEGA*(CO1 - 1D0) + CO2
         ENDIF
      ! END
      STIFFNESS(3) = 0.0001D0/(XS2-XS1)
      STIFFNESS(4) = 0.0001D0/(YS2-YS1)
      IF (XS2-XS1.EQ.0.0D0) STIFFNESS(3) = 10**12
      IF (YS2-YS1.EQ.0.0D0) STIFFNESS(4) = 10**12
      FOC(3)       = HA - 0.0001D0
      FOC(4)       = VA - 0.0001D0

      ELSEIF (LTYPE.EQ."CALC")THEN
                  
      ! CHECK SEABED CONDITION    
      CALL CABLE_LENGTH_ALL_TYPE (ALB,COOR_H,COOR_V,ALENGTH,IEL,"CALC",COOR_REF)
      
      ENDIF

18	RETURN
	END
C	=====================================================================
C	=====================================================================
C	=====================================================================      
	SUBROUTINE SEABED_MOORING (Z,T,A,E,XLO,WO,W,FOC)
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
      
      ALB = XL - F01/W      
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
!	SUBROUTINE CABLE_LENGTH (ALC,NELEMENT_RE,ALENGTH,COOR_H,COOR_V,ITYPE_MOOR)
!	IMPLICIT REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER*4 (I-N)
!      CHARACTER*4 ITYPE_MOOR
!      
!      COMMON A(9000000),IA(9000000)
!      
!      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
!     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
!     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
!      
!      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
!     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
!     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
!     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
!	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
!	5			    LCN,LDIM,LFRE,LSFC,LLOF
!      
!      COMMON /STORE_MOOR_CONNET/ I_MOOR(10000),IN1_MOOR(10000),IN2_MOOR(10000)
!      
!      COMMON /MGRAV/ NGRAV
!      
!      !COMMON /STORE_MOOR/ NSTORE_ELEMENT(2000,50),ALENGTH_OF_MOORING(100),INDEX_GROUP
!      
!      DIMENSION COOR_N1(10000),COOR_N2(10000)
!      DIMENSION NPO_N1(10000),NPO_N2(10000)
!      DIMENSION NSIMI_CONNECT(100),NSIMI_CONNECT_CHECK(100),NSIMI_ELEMENT(100)
!      DIMENSION N_PEAK_POSITION(100),N_PEAK_ELEMENT(100)
!      DIMENSION ASTORE_ELEMENT_N1(100),ASTORE_ELEMENT_N2(100)
!      DIMENSION N_NODE_LAST(2000)
!      DIMENSION AN1_COOR_LOW(100),AN2_COOR_LOW(100),N1_NODE_LOW(100),N2_NODE_LOW(100)
!      
!      
!      IF (ITYPE_MOOR.EQ."CALC") THEN
!      ! INTIAL SETTING
!      COOR_N1 = 0.
!      COOR_N2 = 0.
!      NPO_N1  = 0.
!      NPO_N2  = 0.
!      NSIMI_CONNECT       = 0.
!      NSIMI_CONNECT_CHECK = 0.
!      
!      ! NELE  = NUMBER ELEMENT IN GROUP
!      
!      ! FIND HIGHTEST POSITION OF MOORING LINE
!      ! ***** LIMITED OF THE PROGRAM IS DUE TO THE HIGHTEST POSITION *****
!      DO 100 I = 1,NELE
!          
!          ! FIND CONNECTIVITY POSITION
!          DO 200 J = 1,NELE
!          IF (I_MOOR(J).EQ.I) THEN
!          N_POSITION = J
!          EXIT
!          ENDIF
!200       CONTINUE 
!          
!          ! FIND CONNECTIVITY OF I
!          N1 = IN1_MOOR(N_POSITION)
!          N2 = IN2_MOOR(N_POSITION)
!          
!          IF (NGRAV.EQ.1) THEN     ! X-DIRECTION
!              
!CC          COOR_N1(I) = AX(N1)    
!CC          COOR_N2(I) = AX(N2) 
!	    CALL RELFILL('@XYZ',COOR_N1(I),1,N1,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_N2(I),1,N2,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!          
!          NPO_N1(I)  = N1
!          NPO_N2(I)  = N2
!          
!          ELSEIF (NGRAV.EQ.2)THEN  ! Y-DIRECTION
!              
!CC          COOR_N1(I) = AY(N1)    
!CC          COOR_N2(I) = AY(N2)   
!	    CALL RELFILL('@XYZ',COOR_N1(I),2,N1,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_N2(I),2,N2,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!          
!          NPO_N1(I)  = N1
!          NPO_N2(I)  = N2
!          
!          ELSEIF (NGRAV.EQ.3)THEN  ! Z-DIRECTION
!              
!CC          COOR_N1(I) = AZ(N1)    
!CC          COOR_N2(I) = AZ(N2)
!	    CALL RELFILL('@XYZ',COOR_N1(I),3,N1,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_N2(I),3,N2,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!          
!          NPO_N1(I)  = N1
!          NPO_N2(I)  = N2
!          
!          ENDIF
!      
!100   CONTINUE      
!      ! HIGHTEST POSTION OF MOORING ELEMENT
!      NMAX_POS_N1 = MAXLOC(COOR_N1(1:NELE),DIM=1)
!      NMAX_POS_N2 = MAXLOC(COOR_N2(1:NELE),DIM=1)
!      
!      ! HIGHTEST COORDINATE OF MOORING ELEMENT
!      COOR_MAX_N1 = COOR_N1(NMAX_POS_N1)
!      COOR_MAX_N2 = COOR_N2(NMAX_POS_N2)
!      
!      IF (COOR_MAX_N1.LE.COOR_MAX_N2)THEN
!      COOR_MAX   = COOR_MAX_N2
!      ELSEIF (COOR_MAX_N1.GT.COOR_MAX_N2)THEN
!      COOR_MAX   = COOR_MAX_N1    
!      ENDIF
!      
!      INDEX = 1
!      DO K = 1,NELE
!      ! FIND NUMBER OF PEAK WHICH DUE TO THE POSITION
!      ! THE NUMBER OF CONNECTIVITY HAS BEEN SIMILAR SHOULD BE CHECK IN THIS ROUTINE
!         IF (COOR_MAX.EQ.COOR_N1(K)) THEN
!         NSIMI_CONNECT(INDEX) = NPO_N1(K)  
!         NSIMI_ELEMENT(INDEX) = I_MOOR(K)
!         INDEX = INDEX + 1
!         ENDIF
!         IF (COOR_MAX.EQ.COOR_N2(K)) THEN
!         NSIMI_CONNECT(INDEX) = NPO_N2(K)
!         NSIMI_ELEMENT(INDEX) = I_MOOR(K)
!         INDEX = INDEX + 1
!         ENDIF
!      ENDDO
!      
!      ! INDEX         = TOTAL PEAK NUMBER
!      ! NSIMI_CONNECT = THE PEAK POSITION (NODAL POSITION) WHICH CAN BE SIMILAR VALUE 
!      
!      ! CHECK VALUE OF NSIMI_CONNECT
!      NSIMI_CONNECT_CHECK(1:INDEX) = NSIMI_CONNECT(1:INDEX)
!      DO I =1,INDEX
!          NCONNET = NSIMI_CONNECT_CHECK(I)
!          DO 300 J = 1,INDEX
!             IF (I.EQ.J) THEN 
!             GOTO 300 ! THE VALUE HAS USED TO CHECK (SKIP THIS VALUE)
!             ELSEIF (I.NE.J) THEN
!                 IF (NCONNET.NE.NSIMI_CONNECT(J))THEN
!                 ! DO NOT THING WITH THE VALUE
!                 ELSEIF (NCONNET.EQ.NSIMI_CONNECT(J))THEN
!                 ! DELETE THE VALUE
!                 NSIMI_CONNECT(J) = 0.
!                 NSIMI_ELEMENT(J) = 0.
!                 ENDIF
!             ENDIF
!300      CONTINUE
!      ENDDO
!      
!      ! REARRANGE THE PEAK POSITION IN SINGLE ARRAY
!      INDEX1  = 1 
!      DO I = 1,INDEX
!      IF (NSIMI_CONNECT(I).NE.0) THEN
!      N_PEAK_POSITION(INDEX1) = NSIMI_CONNECT(I)
!      N_PEAK_ELEMENT(INDEX1)  = NSIMI_ELEMENT(I)
!      ENDIF
!      ENDDO
!      
!      ! TRUE LENGTH OF MOORING LINE IN TERM OF GROUP ELEMENT
!      
!      NTYPE   = 0
!      DO I = 1,INDEX1
!      INDEX4  = 1
!       ! FIND CONNECT FROM HIGHTEST POSITION TO LOWEST POSITION
!1000     DO 400 J = 1,NELE
!           IF (N_PEAK_ELEMENT(I).EQ.I_MOOR(J))THEN
!           N1_CONNECT_PEAK = IN1_MOOR(J) ! FIRST POSITION
!           N2_CONNECT_PEAK = IN2_MOOR(J) ! SECOND POSITION
!           EXIT
!           ENDIF
!400      CONTINUE
!      
!      IF (NTYPE.EQ.0)THEN ! FIRST LOOP OF PEAK POSITION 
!      ! FIND NEXT CONNECTION FROM PEAK POSITION
!         NSTORE_ELEMENT(INDEX4,I) = N_PEAK_ELEMENT(I) ! STORAGE THE ELEMENT IN TERM OF GROUP
!         INDEX4 = INDEX4 + 1
!         IF (N_PEAK_POSITION(INDEX1).EQ.N1_CONNECT_PEAK)THEN
!             N_MAIN = N2_CONNECT_PEAK
!             N_PEAK = N1_CONNECT_PEAK
!         ELSEIF (N_PEAK_POSITION(INDEX1).EQ.N2_CONNECT_PEAK)THEN
!             N_MAIN = N1_CONNECT_PEAK
!             N_PEAK = N2_CONNECT_PEAK
!         ENDIF
!      ELSEIF (NTYPE.EQ.1)THEN ! SECOND LOOP OF LAST POSTION
!          
!         DO 500 LLL = 1,NELE
!         IF (N_PEAK_ELEMENT(LLL).EQ.0) THEN
!         EXIT
!         ELSEIF (N_PEAK_ELEMENT(LLL).NE.0) THEN
!         IF (NGRAV.EQ.1) THEN     ! X-DIRECTION
!             
!         ELSEIF (NGRAV.EQ.2) THEN ! Y-DIRECTION
!             
!         ELSEIF (NGRAV.EQ.3) THEN ! Z-DIRECTION
!             
!CC          COOR_NODE1 = AZ(N1_CONNECT_PEAK) 
!CC          COOR_NODE2 = AZ(N2_CONNECT_PEAK) 
!	    CALL RELFILL('@XYZ',COOR_NODE1,3,N1_CONNECT_PEAK,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_NODE2,3,N2_CONNECT_PEAK,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!          
!          IF (COOR_NODE1.GT.COOR_NODE2)THEN
!          N_MAIN = N2_CONNECT_PEAK
!          N_PEAK = N1_CONNECT_PEAK
!          N_NODE_LAST(I) = N_MAIN
!          ELSEIF (COOR_NODE1.LE.COOR_NODE2)THEN
!          N_MAIN = N1_CONNECT_PEAK
!          N_PEAK = N2_CONNECT_PEAK
!          N_NODE_LAST(I) = N_MAIN
!          ENDIF
!          
!         ENDIF
!         ENDIF  
!500     CONTINUE
!          
!      ENDIF
!      
!      IF (NGRAV.EQ.1) THEN     ! X-DIRECTION
!CC         COOR_MAIN_NODE = AX(N_MAIN)
!CC         COOR_PEAK_NODE = AX(N_PEAK)
!	    CALL RELFILL('@XYZ',COOR_MAIN_NODE,1,N_MAIN,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_PEAK_NODE,1,N_PEAK,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.2) THEN ! Y-DIRECTION
!CC         COOR_MAIN_NODE = AY(N_MAIN) 
!CC         COOR_PEAK_NODE = AY(N_PEAK)
!	    CALL RELFILL('@XYZ',COOR_MAIN_NODE,2,N_MAIN,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_PEAK_NODE,2,N_PEAK,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.3) THEN ! Z-DIRECTION
!CC         COOR_MAIN_NODE = AZ(N_MAIN) 
!CC         COOR_PEAK_NODE = AZ(N_PEAK) 
!	    CALL RELFILL('@XYZ',COOR_MAIN_NODE,3,N_MAIN,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_PEAK_NODE,3,N_PEAK,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ENDIF 
!      
!      ! CHECK ERROR OF THE POSITION OF NODE
!      IF (COOR_MAIN_NODE.GT.COOR_PEAK_NODE)THEN
!      ! FAIL > 
!      ELSEIF (COOR_MAIN_NODE.LE.COOR_PEAK_NODE)THEN
!      ! PASS > CONTINUE
!      ENDIF
!      
!      ! FIND NEXT ELEMENT
!      INDEX2 = 1D0
!      INDEX3 = 1D0
!      DO L = 1, NELE
!          IF (N_MAIN.EQ.IN1_MOOR(L))THEN
!              ASTORE_ELEMENT_N1(INDEX2) = I_MOOR(L)
!              INDEX2 = INDEX2 + 1
!          ELSEIF (N_MAIN.EQ.IN2_MOOR(L))THEN    
!              ASTORE_ELEMENT_N2(INDEX3) = I_MOOR(L)
!              INDEX3 = INDEX3 + 1
!          ENDIF
!      ENDDO
!      
!      ! CUTTING REULT
!      DO L =1,INDEX2
!          IF (ASTORE_ELEMENT_N1(L).EQ.N_PEAK_ELEMENT(I)) THEN
!          ASTORE_ELEMENT_N1(L) = 0.0D0
!          ENDIF
!      ENDDO
!      
!      DO L =1,INDEX3
!          IF (ASTORE_ELEMENT_N2(L).EQ.N_PEAK_ELEMENT(I)) THEN
!          ASTORE_ELEMENT_N2(L) = 0.0D0
!          ENDIF
!      ENDDO
!      
!      ! NUMBER OF ELEMENT HAS CONNECTED
!      NELEMENT_NEXT1 = MAXVAL(ASTORE_ELEMENT_N1(1:INDEX2))
!      NELEMENT_NEXT2 = MAXVAL(ASTORE_ELEMENT_N2(1:INDEX3))
!      
!      IF (NELEMENT_NEXT1.EQ.0.AND.NELEMENT_NEXT2.EQ.0) GOTO 1100 ! CONTINUE ON THE NEXT LOOP
!      IF (NELEMENT_NEXT1.GE.NELEMENT_NEXT2) THEN
!      N_PEAK_ELEMENT(I)        = NELEMENT_NEXT1
!      NSTORE_ELEMENT(INDEX4,I) = NELEMENT_NEXT1 ! STORAGE THE ELEMENT IN TERM OF GROUP
!      INDEX4 = INDEX4 + 1
!      NTYPE = 1
!      GOTO 1000
!      ELSEIF (NELEMENT_NEXT1.LE.NELEMENT_NEXT2) THEN
!      N_PEAK_ELEMENT(I)        = NELEMENT_NEXT2
!      ! THIS STORAGE WILL BE USED FOR CALCULATE POSITION OF ELEMENT
!      NSTORE_ELEMENT(INDEX4,I) = NELEMENT_NEXT2 ! STORAGE THE ELEMENT IN TERM OF GROUP
!      INDEX4 = INDEX4 + 1
!      NTYPE = 1
!      GOTO 1000
!      ENDIF
!      
!1100  ENDDO
!      
!      ! N_PEAK_ELEMENT IS BECOME A LAST ELEMENT OF THE MOORING
!      ! N_PEAK_POSITION IS BECOME A PEAK POSITION 
!      ! N_NODE_LAST IS BECOME A LAST POSITION
!      
!      ! MAKE A GROUP OF ELEMENT AND CALCULATE UNSTRENGTH-LENGTH  
!      DO I = 1,INDEX1 
!          ! PEAK POSITION
!          N_NODE_PEAK_POINT = N_PEAK_POSITION(I)
!          ! LAST POSITION
!          N_NODE_LAST_POINT = N_NODE_LAST(I)
!          
!          ! CALCULATE UNSTRENGTH-LENGTH 3D
!CC          COOR_X_PEAK  = AX(N_NODE_PEAK_POINT)   
!CC          COOR_X_LAST  = AX(N_NODE_LAST_POINT)  
!CC          COOR_Y_PEAK  = AY(N_NODE_PEAK_POINT)   
!CC          COOR_Y_LAST  = AY(N_NODE_LAST_POINT)
!CC          COOR_Z_PEAK  = AZ(N_NODE_PEAK_POINT)
!CC          COOR_Z_LAST  = AZ(N_NODE_LAST_POINT)
!	    CALL RELFILL('@XYZ',COOR_X_PEAK,1,N_NODE_PEAK_POINT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_X_LAST,1,N_NODE_LAST_POINT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_Y_PEAK,2,N_NODE_PEAK_POINT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_Y_LAST,2,N_NODE_LAST_POINT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_Z_PEAK,3,N_NODE_PEAK_POINT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!	    CALL RELFILL('@XYZ',COOR_Z_LAST,3,N_NODE_LAST_POINT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!          
!          ALENGTH_OF_MOORING(I) = SQRT( (COOR_X_LAST-COOR_X_PEAK)**2 + (COOR_Y_LAST-COOR_Y_PEAK)**2 + (COOR_Z_LAST-COOR_Z_PEAK)**2)
!      ENDDO
!      
!      ! CALCULATE INDEX OF STORAGE
!      INDEX_GROUP = 0
!      DO I = 1,50
!          IF (NSTORE_ELEMENT(1,I).NE.0) INDEX_GROUP = INDEX_GROUP + 1
!      ENDDO
!      
!      ! CALCULATE THE COORDINATE POSITION SYSTEM
!      ! REQUIRED VECTOR POSITION
!      IF (NGRAV.EQ.1)    THEN ! X-DIRECTION
!      COOR_H1 = COOR_Y_PEAK - COOR_Y_LAST
!      COOR_H2 = COOR_Z_PEAK - COOR_Z_LAST
!      COOR_H  = SQRT(COOR_H1**2+COOR_H2**2)
!      COOR_V  = COOR_X_PEAK - COOR_X_LAST
!      ELSEIF (NGRAV.EQ.2)THEN ! Y-DIRECTION
!      COOR_H1 = COOR_X_PEAK - COOR_X_LAST
!      COOR_H2 = COOR_Z_PEAK - COOR_Z_LAST
!      COOR_H  = SQRT(COOR_H1**2+COOR_H2**2)
!      COOR_V  = COOR_Y_PEAK - COOR_Y_LAST
!      ELSEIF (NGRAV.EQ.3)THEN ! Z-DIRECTION
!      COOR_H1 = COOR_X_LAST - COOR_X_PEAK
!      COOR_H2 = COOR_Y_LAST - COOR_Y_PEAK
!      COOR_H  = SQRT(COOR_H1**2+COOR_H2**2)
!      COOR_V  = COOR_Z_PEAK - COOR_Z_LAST  
!      ENDIF
!          
!      ELSEIF (ITYPE_MOOR.EQ."READ")THEN
!      ! RECICVED NUMBER OF ELEMENT
!      
!      DO 650 J = 1, INDEX_GROUP
!          DO 600 I = 1,NELE
!              IF (NSTORE_ELEMENT(I,J).EQ.0) EXIT
!              IF (NSTORE_ELEMENT(I,J).EQ.NELEMENT_RE)THEN
!              NOPTION = 1
!              EXIT
!              ENDIF
!600       CONTINUE
!      IF (NOPTION.EQ.1) EXIT
!650   CONTINUE
!      
!      ! J IS NUMBER OF GROUP
!      ! SEND UNSTRENGTH-LENGTH OF LINE
!      ALC  = ALENGTH_OF_MOORING(J)
!      
!      ! FIND THE POSITION OF ELEMENT IN MOORING LINE (ALENGTH)
!      N1_COOR_LOW = 0.
!      N2_COOR_LOW = 0.
!      N1_NODE_LOW = 0.
!      N2_NODE_LOW = 0.
!      INDEX10     = 0.
!      DO 700 K = 1,NELE
!          IF (NSTORE_ELEMENT(K,J).EQ.0) EXIT 
!          INDEX10 = INDEX10 + 1
!          ! FIND LOWEST POSITION OF MOORING
!          DO KK = 1,NELE
!           IF (I_MOOR(KK).EQ.NSTORE_ELEMENT(K,J))THEN
!               IF (NGRAV.EQ.1)THEN
!                   
!CC               AN1_COOR_LOW(K) = AX(IN1_MOOR(KK))
!CC               AN2_COOR_LOW(K) = AX(IN2_MOOR(KK))    
!                CALL RELFILL('@XYZ',AN1_COOR_LOW(K),1,IN1_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                CALL RELFILL('@XYZ',AN2_COOR_LOW(K),1,IN2_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!          
!               N1_NODE_LOW(K) = IN1_MOOR(KK)
!               N2_NODE_LOW(K) = IN2_MOOR(KK)
!               ELSEIF (NGRAV.EQ.2)THEN
!                   
!CC               AN1_COOR_LOW(K) = AY(IN1_MOOR(KK))
!CC               AN2_COOR_LOW(K) = AY(IN2_MOOR(KK)) 
!                CALL RELFILL('@XYZ',AN1_COOR_LOW(K),2,IN1_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                CALL RELFILL('@XYZ',AN2_COOR_LOW(K),2,IN2_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!              
!               N1_NODE_LOW(K) = IN1_MOOR(KK)
!               N2_NODE_LOW(K) = IN2_MOOR(KK)
!               ELSEIF (NGRAV.EQ.3)THEN
!                   
!CC               AN1_COOR_LOW(K) = AZ(IN1_MOOR(KK))
!CC               AN2_COOR_LOW(K) = AZ(IN2_MOOR(KK))
!                CALL RELFILL('@XYZ',AN1_COOR_LOW(K),3,IN1_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                CALL RELFILL('@XYZ',AN2_COOR_LOW(K),3,IN2_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                
!               N1_NODE_LOW(K) = IN1_MOOR(KK)
!               N2_NODE_LOW(K) = IN2_MOOR(KK)
!               ENDIF
!           ENDIF
!          ENDDO
!700   CONTINUE
!      
!      AN1_MIN_LOW_VALUE = MINVAL(AN1_COOR_LOW(1:INDEX10))
!      AN2_MIN_LOW_VALUE = MINVAL(AN2_COOR_LOW(1:INDEX10))
!      N1_MIN_LOW_POSITION = MINLOC(AN1_COOR_LOW(1:INDEX10),DIM=1)
!      N2_MIN_LOW_POSITION = MINLOC(AN1_COOR_LOW(1:INDEX10),DIM=1)
!      
!      IF (AN1_MIN_LOW_VALUE.GT.AN2_MIN_LOW_VALUE)THEN
!      AN_LOWEST_VALUE = AN2_MIN_LOW_VALUE
!      N_LOWEST_POSITION = N2_NODE_LOW(N2_MIN_LOW_POSITION)
!      ELSEIF (N1_MIN_LOW_VALUE.LE.N2_MIN_LOW_VALUE)THEN
!      AN_LOWEST_VALUE = AN1_MIN_LOW_VALUE    
!      N_LOWEST_POSITION = N1_NODE_LOW(N1_MIN_LOW_POSITION)
!      ENDIF
!      
!CC      DO 790 I = 1,NNOT
!CC          IF (NND(I).EQ.N_LOWEST_POSITION)THEN
!CC             AX_LOW = AX(I)
!CC             AY_LOW = AY(I)
!CC             AZ_LOW = AZ(I)
!CC          EXIT
!CC          ENDIF
!CC790   CONTINUE
!      
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!	CALL RELFILL('@XYZ',AX_LOW,1,N_LOWEST_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AY_LOW,2,N_LOWEST_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AZ_LOW,3,N_LOWEST_POSITION,0)  !GETTING NODAL COORDINATE	  
!      
!      ! FIND A COMPONENT OF THE ELEMENT
!      ! FOR CALCULATE LENGTH FROM THE ORIGINAL POSITION
!      DO 800 I = 1,NELE
!         IF (I_MOOR(I).EQ.NELEMENT_RE)THEN
!         N1_NODE_ELEMENT = IN1_MOOR(NELEMENT_RE)
!         N2_NODE_ELEMENT = IN2_MOOR(NELEMENT_RE)
!         EXIT
!         ENDIF
!800   CONTINUE
!      
!CC      DO 810 I = 1,NNOT
!CC      IF (NND(I).EQ.N1_NODE_ELEMENT)THEN
!CC         IF (NGRAV.EQ.1)THEN
!CC         AC_N1 = AX(N1_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.2)THEN
!CC         AC_N1 = AY(N1_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.3)THEN
!CC         AC_N1 = AZ(N1_NODE_ELEMENT)
!CC         ENDIF
!CC      EXIT
!CC      ENDIF
!CC810   CONTINUE
!
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!      IF (NGRAV.EQ.1)THEN
!          CALL RELFILL('@XYZ',AC_N1,1,N1_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.2)THEN
!          CALL RELFILL('@XYZ',AC_N1,2,N1_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.3)THEN
!          CALL RELFILL('@XYZ',AC_N1,3,N1_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ENDIF 
!        
!      
!CC      DO 820 I = 1,NNOT
!CC      IF (NND(I).EQ.N1_NODE_ELEMENT)THEN
!CC         IF (NGRAV.EQ.1)THEN
!CC         AC_N2 = AX(N2_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.2)THEN
!CC         AC_N2 = AY(N2_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.3)THEN
!CC         AC_N2 = AZ(N2_NODE_ELEMENT)
!CC         ENDIF
!CC      EXIT
!CC      ENDIF
!CC820   CONTINUE
!      
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!      IF (NGRAV.EQ.1)THEN
!          CALL RELFILL('@XYZ',AC_N2,1,N2_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.2)THEN
!          CALL RELFILL('@XYZ',AC_N2,2,N2_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.3)THEN
!          CALL RELFILL('@XYZ',AC_N2,3,N2_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ENDIF 
!      
!      
!      IF (AC_N1.GT.AC_N2)THEN
!      NNODE_POSITION = N1_NODE_ELEMENT  
!      ELSEIF (AC_N1.LE.AC_N2)THEN
!      NNODE_POSITION = N2_NODE_ELEMENT      
!      ENDIF
!      
!CC      DO 830 I =1,NNOT
!CC         IF (NND(I).EQ.NNODE_POSITION)THEN
!CC         AX_POSITION = AX(I)
!CC         AY_POSITION = AY(I)
!CC         AZ_POSITION = AZ(I)
!CC         EXIT    
!CC         ENDIF
!CC830   CONTINUE
!      
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!	CALL RELFILL('@XYZ',AX_POSITION,1,NNODE_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AY_POSITION,2,NNODE_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AZ_POSITION,3,NNODE_POSITION,0)  !GETTING NODAL COORDINATE	
!      
!      ! CALCULATE ALENGTH 1D PLAN
!      ALENGTH = SQRT((AX_LOW-AX_POSITION)**2+(AY_LOW-AY_POSITION)**2+(AZ_LOW-AZ_POSITION)**2)
!      
!      ELSEIF (ITYPE_MOOR.EQ."REAR") THEN
!          
!      ! FIND THE POSITION OF ELEMENT IN MOORING LINE (ALENGTH)
!      N1_COOR_LOW = 0.
!      N2_COOR_LOW = 0.
!      N1_NODE_LOW = 0.
!      N2_NODE_LOW = 0.
!      INDEX10     = 0.
!      DO 900 K = 1,NELE
!          IF (NSTORE_ELEMENT(K,J).EQ.0) EXIT 
!          INDEX10 = INDEX10 + 1
!          ! FIND LOWEST POSITION OF MOORING
!          DO KK = 1,NELE
!           IF (I_MOOR(KK).EQ.NSTORE_ELEMENT(K,J))THEN
!               IF (NGRAV.EQ.1)THEN
!                   
!CC               AN1_COOR_LOW(K) = AX(IN1_MOOR(KK))
!CC               AN2_COOR_LOW(K) = AX(IN2_MOOR(KK))    
!                CALL RELFILL('@XYZ',AN1_COOR_LOW(K),1,IN1_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                CALL RELFILL('@XYZ',AN2_COOR_LOW(K),1,IN2_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                
!               N1_NODE_LOW(K) = IN1_MOOR(KK)
!               N2_NODE_LOW(K) = IN2_MOOR(KK)
!               ELSEIF (NGRAV.EQ.2)THEN
!                   
!CC               AN1_COOR_LOW(K) = AY(IN1_MOOR(KK))
!CC               AN2_COOR_LOW(K) = AY(IN2_MOOR(KK))   
!                CALL RELFILL('@XYZ',AN1_COOR_LOW(K),2,IN1_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                CALL RELFILL('@XYZ',AN2_COOR_LOW(K),2,IN2_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                
!               N1_NODE_LOW(K) = IN1_MOOR(KK)
!               N2_NODE_LOW(K) = IN2_MOOR(KK)
!               ELSEIF (NGRAV.EQ.3)THEN
!                   
!CC               AN1_COOR_LOW(K) = AZ(IN1_MOOR(KK))
!CC               AN2_COOR_LOW(K) = AZ(IN2_MOOR(KK))
!                CALL RELFILL('@XYZ',AN1_COOR_LOW(K),3,IN1_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                CALL RELFILL('@XYZ',AN2_COOR_LOW(K),3,IN2_MOOR(KK),0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!                
!               N1_NODE_LOW(K) = IN1_MOOR(KK)
!               N2_NODE_LOW(K) = IN2_MOOR(KK)
!               ENDIF
!           ENDIF
!          ENDDO
!900   CONTINUE
!      
!      AN1_MIN_LOW_VALUE = MINVAL(AN1_COOR_LOW(1:INDEX10))
!      AN2_MIN_LOW_VALUE = MINVAL(AN2_COOR_LOW(1:INDEX10))
!      N1_MIN_LOW_POSITION = MINLOC(AN1_COOR_LOW(1:INDEX10),DIM=1)
!      N2_MIN_LOW_POSITION = MINLOC(AN1_COOR_LOW(1:INDEX10),DIM=1)
!      
!      IF (AN1_MIN_LOW_VALUE.GT.AN2_MIN_LOW_VALUE)THEN
!      AN_LOWEST_VALUE = AN2_MIN_LOW_VALUE
!      N_LOWEST_POSITION = N2_NODE_LOW(N2_MIN_LOW_POSITION)
!      ELSEIF (N1_MIN_LOW_VALUE.LE.N2_MIN_LOW_VALUE)THEN
!      AN_LOWEST_VALUE = AN1_MIN_LOW_VALUE    
!      N_LOWEST_POSITION = N1_NODE_LOW(N1_MIN_LOW_POSITION)
!      ENDIF
!      
!CC      DO 990 I = 1,NNOT
!CC          IF (NND(I).EQ.N_LOWEST_POSITION)THEN
!CC             AX_LOW = AX(I)
!CC             AY_LOW = AY(I)
!CC             AZ_LOW = AZ(I)
!CC          EXIT
!CC          ENDIF
!CC990   CONTINUE
!      
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!	CALL RELFILL('@XYZ',AX_LOW,1,N_LOWEST_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AY_LOW,2,N_LOWEST_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AZ_LOW,3,N_LOWEST_POSITION,0)  !GETTING NODAL COORDINATE
!      
!      
!      ! FIND A COMPONENT OF THE ELEMENT
!      ! FOR CALCULATE LENGTH FROM THE ORIGINAL POSITION
!      DO 1005 I = 1,NELE
!         IF (I_MOOR(I).EQ.NELEMENT_RE)THEN
!         N1_NODE_ELEMENT = IN1_MOOR(NELEMENT_RE)
!         N2_NODE_ELEMENT = IN2_MOOR(NELEMENT_RE)
!         EXIT
!         ENDIF
!1005  CONTINUE
!      
!CC      DO 1010 I = 1,NNOT
!CC      IF (NND(I).EQ.N1_NODE_ELEMENT)THEN
!CC         IF (NGRAV.EQ.1)THEN
!CC         AC_N1 = AX(N1_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.2)THEN
!CC         AC_N1 = AY(N1_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.3)THEN
!CC         AC_N1 = AZ(N1_NODE_ELEMENT)
!CC         ENDIF
!CC      EXIT
!CC      ENDIF
!CC1010  CONTINUE
!      
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!      IF (NGRAV.EQ.1)THEN
!          CALL RELFILL('@XYZ',AC_N1,1,N1_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.2)THEN
!          CALL RELFILL('@XYZ',AC_N1,2,N1_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.3)THEN
!          CALL RELFILL('@XYZ',AC_N1,3,N1_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ENDIF 
!      
!      
!CC      DO 1020 I = 1,NNOT
!CC      IF (NND(I).EQ.N2_NODE_ELEMENT)THEN
!CC         IF (NGRAV.EQ.1)THEN
!CC         AC_N2 = AX(N2_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.2)THEN
!CC         AC_N2 = AY(N2_NODE_ELEMENT)
!CC         ELSEIF (NGRAV.EQ.3)THEN
!CC         AC_N2 = AZ(N2_NODE_ELEMENT)
!CC         ENDIF
!CC      EXIT
!CC      ENDIF
!CC1020   CONTINUE
!      
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!      IF (NGRAV.EQ.1)THEN
!          CALL RELFILL('@XYZ',AC_N2,1,N2_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.2)THEN
!          CALL RELFILL('@XYZ',AC_N2,2,N2_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ELSEIF (NGRAV.EQ.3)THEN
!          CALL RELFILL('@XYZ',AC_N2,3,N2_NODE_ELEMENT,0)  !GETTING NODAL COORDINATE ... SONGSAK OCT2019
!      ENDIF 
!      
!      
!      IF (AC_N1.LE.AC_N2)THEN
!      NNODE_POSITION = N1_NODE_ELEMENT  
!      ELSEIF (AC_N1.GT.AC_N2)THEN
!      NNODE_POSITION = N2_NODE_ELEMENT      
!      ENDIF
!      
!CC      DO 1030 I =1,NNOT
!CC         IF (NND(I).EQ.NNODE_POSITION)THEN
!CC         AX_POSITION = AX(I)
!CC         AY_POSITION = AY(I)
!CC         AZ_POSITION = AZ(I)
!CC         EXIT    
!CC         ENDIF
!CC1030   CONTINUE
!      
!C     ----------------------------------------------------------------  
!C     SONGSAK SUPPRESS ABOVE LOOP ... REPLACE BY NEXT BLOCK ... OCT2019    
!	CALL RELFILL('@XYZ',AX_POSITION,1,NNODE_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AY_POSITION,2,NNODE_POSITION,0)  !GETTING NODAL COORDINATE
!	CALL RELFILL('@XYZ',AZ_POSITION,3,NNODE_POSITION,0)  !GETTING NODAL COORDINATE	
!      
!      ! CALCULATE ALENGTH 1D PLAN
!      ALENGTH = SQRT((AX_LOW-AX_POSITION)**2+(AY_LOW-AY_POSITION)**2+(AZ_LOW-AZ_POSITION)**2)
!          
!  
!      ENDIF
!          
!      END
C	=====================================================================
C	=====================================================================      
	SUBROUTINE OUTPUT_MOORING
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /OUT_MOOR/ OUTPUT_DATA(2)
      OPEN(UNIT=1000,FILE='Mooring Summation.dat'             ,STATUS='UNKNOWN')
      
      WRITE (1000,1) OUTPUT_DATA(1)
1     FORMAT ("FAIRLEAD FORCE H",2X,E12.5)
      WRITE (1000,2) OUTPUT_DATA(2)
2     FORMAT ("FAIRLEAD FORCE V",2X,E12.5)
      
      FAIRLEADH  = OUTPUT_DATA(1)
      FAIRLEADV  = OUTPUT_DATA(2)
      
      ! CALCULATE TENSION
      
      END

