C     --------------------------------------------------------------
C     --------------------------------------------------------------
      SUBROUTINE SOLINEW(PROPM,PROPG,NODEX,WA,WA2,S,COORD,EDIS,EDISI,RE,MWG,
	1				   FIN,MSET,SEDI,SEL,ALPHA,HINFC,AMV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     MAIN PROGRAM FOR THE 3-D SOLID---preintegration
C     EVALUATES THE TANGENTIAL STIFFNESS MATRIX,STRAINS AND STRESSES
C     AND EQUILIBRIUM LOADS FOR THE CURVILINEAR ISOPARAMETRIC
C     HEXAHEDRON (8 TO 21 NODES)
C	--------------------------
C     INPUT VARIABLES
C	---------------
C     PROPM(NMP)    = MATERIAL PROPERTIES (YM,PR,YLD,HP,DEN)
C     PROPG(NGP)    = GEOMETRIC PROPERTIES (NNO)
C     NODEX(NEX)    = LOCATIONS OF EXCESS NODES (MIDSIDE NODES)
C     WA(MWG,NPT)   = WORKING ARRAY (6 STRESSES + (6 STRAINS,YLD,IPEL))
C     COORD(3,NNO)  = CURRENT NODAL COORDINATES X,Y,Z
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     EDISI(NEF)    = CURRENT NODAL DISPLACEMENT INCREMENTS
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     S(NWS)        = ELEMENT STIFFNESS MATRIX (UPPER TRIANG.ROW-WISE)
C     RE(NEF)       = EQUILIBRIUM LOADS AT ELEMENT NODES
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/
C	--------------------------------
C     NAME(2)       = NAME OF ELEMENT MODULE
C     ITYPE         = CODE NUMBER FOR ELEMENT MODULE
C     ISTYP         = ELEMENT SUBTYPE
C     NLOPT         = CODE FOR NONLINEAR OPTION
C     NLOPT=0         LINEAR ANALYSIS
C     NLOPT=1         MATERIALLY NONLINEAR ONLY
C     NLOPT=2,3       TOTAL LAGRANGIAN,UPDATED LAGRANGIAN
C     MTMOD         = CODE FOR MATERIAL MODULE
C     MTMOD=1         LINEAR ELASTIC,ISOTROPIC
C     MTMOD=2         LINEAR ELASTIC,ORTHOTROPIC
C     MTMOD=3         ELASTO-PLASTIC VON-MISES
C     MTMOD=5         CONCRETE WITH CRACKING
C     NSINC         = FACTOR CONTROLLING NUMBER OF SUBINCREMENTS
C     ITOLEY        = TOLERANCE ON YIELD FUNCTION
C     NELE          = NUMBER OF ELEMENTS IN THIS GROUP
C     NMPS          = NUMBER OF MATERIAL PROPERTY SETS
C     NGPS          = NUMBER OF GEOMETRIC PROPERTY SETS
C     NMP           = NUMBER OF MATERIAL PROPERTIES PER SET
C     NGP           = NUMBER OF GEOMETRIC PROPERTIES PER SET
C     NNM           = MAXIMUM NUMBER OF NODES FOR ANY ONE ELEMENT
C     NEX           = MAXIMUM NUMBER OF EXCESS NODES
C     NCO           = NUMBER OF NODAL COORDINATES
C     NNF           = NUMBER OF NODAL DEGREES OF FREEDOM
C     NEF           = MAXIMUM NUMBER OF ELEMENT DEGREES OF FREEDOM
C     NWG           = NUMBER OF STORAGE LOCATIONS AT EACH GAUSS POINT
C     NPT           = NUMBER OF GAUSS POINTS
C     NWA           = SIZE OF WORKING ARRAY
C     NWS           = SIZE OF ELEMENT STIFFNESS MATRIX
C     MEL           = CURRENT ELEMENT NUMBER
C     NNO           = NUMBER OF NODES FOR THIS ELEMENT
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR THIS ELEMENT
C     NELTOT        = TOTAL NUMBER OF ELEMENTS (ALL GROUPS)
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /GAUS/
C	--------------------------------
C     GLOC(4,4)     = NATURAL GAUSS POINT COORDINATES (1*1 TO 4*4)
C     GWT (4,4)     = GAUSS POINT WEIGHTS
C     NGR,NGS,NGT   = NUMBER OF GAUSS POINTS IN RN,SN,TN DIRECTION
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /FLAG/
C	--------------------------------
C     IFPRI,ISPRI   = FLAG FOR PRINTING DISPL.OR STRESSES (ISPRI=0)
C     IFPLO         = FLAG FOR PLOT OUTPUT
C     IFREF         = FLAG FOR REFORMATION OF STIFFNESS (IFREF=0)
C     IFEIG         = FLAG FOR EIGENVALUE SOLUTION (IFEIG=0)
C     ITASK = 1       FIRST ENTRY INTO ELEMENT MODULE
C     ITASK = 2       ENTRY DURING EQUILIBRIUM ITERATIONS
C     ITASK = 3       ENTRY TO WORK OUT STRESSES (LAST STEP ONLY)
C     ITASK = 4       ENTRY TO DETERMINE GEOMETRIC STIFF.MATRIX ONLY
C     KSTEP           CURRENT STEP NUMBER
C     KITE            CURRENT ITERATION NUMBER
C     ----------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

	COMMON /DISPL/ IDISP(4,6)

C     COMPOSITE
      COMMON /DYNA/ CDEN,IMASS

C	EAS
	COMMON /MMENH/ MM,MM1,MM2,NDIMC
	
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
     
CB    SAVE LOCAL DISPLACEMENT    
      COMMON /NEWDATA/ SDISL(1000000,24) 
CB    SAVE NUMBER(NMX) FOR LVNUM ACCORDING TO ELEMENT  
      COMMON /SLVNUM/ INMX(1000000,1)     
      
C     SAVE ELEMENT NUMBER      
      COMMON /NELEM/ IEL 
      
CB    THE THICKNESS
      COMMON /THICK/THK
      COMMON /THICKM/THKM(1000000)
      
      
C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(300),COORD(3,8)
      DIMENSION EDIS(24),EDISI(24),RE(24),XJ(3,3),XJI(3,3),COORDI(3,8)
      DIMENSION STRAIN(9),QSTRAI(9),STRESS(9),TAU(9)
C
      DIMENSION VR(3),VS(3),VT(3),SH(4),P(2,4),SHR(4),SHS(4),SHT(4)
      DIMENSION SHR1(4),SHS1(4),SHT1(4)
	DIMENSION XY(3,4),PXY(3,4),DE(9,9)
	DIMENSION STRESSO(9),DEO(9,9),DEO2(9,9),STRAINO(3),STRAINW(3),STRAINB(3),STRESSL(9),STRESS_FIBER(9) ! FOR NEW STRESS BY BJ
	
C
      DIMENSION BB(9,24),U(3,4),PU(3,4),BBL(9,24),DISL(24)
      DIMENSION XJA(3,3),XJB(3,3),XJC(3,3),XJD(3,3)
	DIMENSION PA(4),PB(4),PC(4),PD(4)
	DIMENSION VRA(3),VSA(3),VTA(3)
      DIMENSION VRB(3),VSB(3),VTB(3)
	DIMENSION VRC(3),VSC(3),VTC(3)
	DIMENSION VRD(3),VSD(3),VTD(3)
	DIMENSION SHA(4),SHB(4),SHC(4),SHD(4)
	DIMENSION BAT(4,24),RA(4),SA(4)
	DIMENSION DUD(9),DPUD(9),STRAI(9),DUD1(9),DPUD1(9)
	DIMENSION TAUA(3),BB1(3,24)
	DIMENSION FRS(4),GRS(4)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
	DIMENSION LVNUM(8)

C	EAS
	DIMENSION TTI(9,9),GE(9,MM)
	DIMENSION REAS(24),SEAS(24,24)
	DIMENSION EAS(9),ALPHA(MM)
      DIMENSION SED(MM,MM),SEL(MM,24),SEDI(MM,MM),RH(MM)
	DIMENSION HINFC(MM)

C	MASS
	DIMENSION HMAS(21),PMAS(3,21),XJIM(3,3),XJM(3,3)

C     LARGE ROTATION
      DIMENSION REDIS(24),XYI(3,4),PXYI(3,4),H8(8),P8(3,8)
      
C     DISPLACEMENT U + PU
      DIMENSION REDISS(24),REDISS1(24)

C     COMPOSITE      
      DIMENSION AMV(3)
      
C     TRANSFORMATION MATRIX BY BJ      
      DIMENSION T1(6,6),TLG(3,3),TLGT(6,6),TG1(6,6),TG2(24,24)      
      
      DIMENSION WA2(MWG,NPT) !---------NEW WORKING ARRAY BY BJ            
      
	K = 0
	DO J=1,8
	DO I=1,3
	K = K+1
	COORDI(I,J) = COORD(I,J)
	IF(NLOPT.EQ.3) COORDI(I,J) = COORD(I,J) - EDIS(K)
	ENDDO
	ENDDO

C     ------------------------------------------------------
C	DETECTED THE THIN DIRECTION
C     ------------------------------------------------------
	XYZ1 = SQRT( (COORDI(1,2)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,2)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,2)-COORDI(3,1))**2.0 )
	XYZ2 = SQRT( (COORDI(1,4)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,4)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,4)-COORDI(3,1))**2.0 )
	XYZ3 = SQRT( (COORDI(1,5)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,5)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,5)-COORDI(3,1))**2.0 )

	NMX = 0
	IF(XYZ1.LT.XYZ2.AND.XYZ1.LT.XYZ3) NMX = 1
	IF(XYZ2.LT.XYZ1.AND.XYZ2.LT.XYZ3) NMX = 2
	IF(XYZ3.LT.XYZ1.AND.XYZ3.LT.XYZ2) NMX = 3

	IF(NMX.EQ.0) THEN
	IF(XYZ1.EQ.XYZ2) NMX = 1
	IF(XYZ1.EQ.XYZ3) NMX = 1
	IF(XYZ2.EQ.XYZ3) NMX = 2
	ENDIF

CB    DETECT THICKNESS
      IF(NMX.EQ.1) THK = XYZ1      
      IF(NMX.EQ.2) THK = XYZ2
      IF(NMX.EQ.3) THK = XYZ3
      
      THKM(IEL) = THK !SAVE THINKNESS INTO MATRIX
      
CB    SAVE INMX ACCORDING TO ELEMENT      
      INMX(IEL,1) = NMX
      
	IF(NMX.EQ.1) THEN
	LVNUM(1) = 1
	LVNUM(2) = 4
	LVNUM(3) = 8
	LVNUM(4) = 5
	LVNUM(5) = 2
	LVNUM(6) = 3
	LVNUM(7) = 7
	LVNUM(8) = 6
	ELSEIF(NMX.EQ.2) THEN
	LVNUM(1) = 1
	LVNUM(2) = 5
	LVNUM(3) = 6
	LVNUM(4) = 2
	LVNUM(5) = 4
	LVNUM(6) = 8
	LVNUM(7) = 7
	LVNUM(8) = 3
	ELSEIF(NMX.EQ.3) THEN
	LVNUM(1) = 1
	LVNUM(2) = 2
	LVNUM(3) = 3
	LVNUM(4) = 4
	LVNUM(5) = 5
	LVNUM(6) = 6
	LVNUM(7) = 7
	LVNUM(8) = 8
	ENDIF
 

C     ------------------------------------------------------
C     GET INTEGRATED MATERIAL MATRIX FOR ISOTROPIC MATERIALS
C     ------------------------------------------------------
      CALL HOKLAW (PROPM,PROPG,1)
      CALL DEMAT_SHL(DE,DEO,DEO2)
 

	CALL SLFACE(NMX)

	K = 0
	DO I = 1,24
	RE(I) = 0.0D0
	DO J = I,24
	K = K + 1
	S(K) = 0.0D0
	ENDDO
	ENDDO
	SED = 0.0
	SEDI = 0.0
	SEL = 0.0
	RH  = 0.0

C     ------------------------------------------------------

C     --------------------------------------------
C     MIDSURFACE COORDINATES AND DIRECTION VECTORS
C     --------------------------------------------
      DO I=1,3
	DO J=1,4
	JJ  = LVNUM(J)
	JJJ = LVNUM(J+4)
       XY(I,J)=(COORD(I,JJ)+COORD(I,JJJ))/2.
      PXY(I,J)=(COORD(I,JJ)-COORD(I,JJJ))/2.
	
       XYI(I,J)=(COORDI(I,JJ)+COORDI(I,JJJ))/2.
      PXYI(I,J)=(COORDI(I,JJ)-COORDI(I,JJJ))/2.
      
C	XY(I,J) =(COORD(I,J)+COORD(I,J+4))/2.
C	PXY(I,J)=(COORD(I,J)-COORD(I,J+4))/2.
	ENDDO
	ENDDO

C     --------------------------------------------
C     EAS TRANSFORMATION MATRIX
C     --------------------------------------------
	RN = 0.0
	SN = 0.0
	CALL SHAP8NEW(RN,SN,SH,P)
	CALL COVAB(XY,PXY,SH,P,VR,VS,VT,XJ,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DETO,MEL,0)
	CALL TNEAS8(XJ,TTI)

	BAT(1:4,1:24) = 0.0D0
C     --------------------------------------------
C     CALCULATE JACOBEAN MATRIX AT SAMPLING POINTS
C     --------------------------------------------
C     ----FOR SHEAR PART----
      RI=-1.0
	SI= 0.0
      CALL SHAP8NEW(RI,SI,SHA,P)
	DO I=1,4
	PA(I)=P(2,I)
	ENDDO
      CALL COVAB(XY,PXY,SHA,P,VRA,VSA,VTA,XJA,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	RI= 0.0
	SI=-1.0
	CALL SHAP8NEW(RI,SI,SHB,P)
	DO I=1,4
	PB(I)=P(1,I)
	ENDDO
      CALL COVAB(XY,PXY,SHB,P,VRB,VSB,VTB,XJB,XJI,
     @    	   SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	RI= 1.0
	SI= 0.0
	CALL SHAP8NEW(RI,SI,SHC,P)
	DO I=1,4
	PC(I)=P(2,I)
	ENDDO
      CALL COVAB(XY,PXY,SHC,P,VRC,VSC,VTC,XJC,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	RI= 0.0
	SI= 1.0
	CALL SHAP8NEW(RI,SI,SHD,P)
	DO I=1,4
	PD(I)=P(1,I)
	ENDDO
      CALL COVAB(XY,PXY,SHD,P,VRD,VSD,VTD,XJD,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
C     -----FOR THICKNESS PART-------
      RA(1) = 1.0D0
	SA(1) = 1.0D0
	RA(2) =-1.0D0
	SA(2) = 1.0D0
	RA(3) =-1.0D0
	SA(3) =-1.0D0
	RA(4) = 1.0D0
	SA(4) =-1.0D0
	DO I=1,4
	CALL SHAP8NEW(RA(I),SA(I),SH,P)
	CALL COVAB(XY,PXY,SH,P,VR,VS,VT,XJ,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	CALL SOBATH(BAT,VR,VS,VT,XJ,SH,I)
	ENDDO


      IF (ITASK.NE.5)  GOTO 50
      
	IF (MTMOD.EQ.2) PROPM(5) = CDEN
C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
	MGR = 3
	MGS = 3
	MGT = 3
	IPT = 0
      DO 90  IGR=1,MGR
	RI = GLOC(IGR,MGR)
	DO 90  IGS=1,MGS
	SI = GLOC(IGS,MGS)
	DO 90  IGT=1,MGT
	TI = GLOC(IGT,MGT)
            
	WT = GWT(IGR,MGR)*GWT(IGS,MGS)*GWT(IGT,MGT)
	IPT = IPT+1

C     ---------------------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ---------------------------------------------------

      CALL SHAP3D_S (RI,SI,TI,HMAS,PMAS,NODEX,NNO)
      CALL JACO3D_S (COORD,P,XJM,XJIM,DET,MEL,NNO)
      DVOL = WT*DET

C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
      CALL SOMASS_S (S,HMAS,PROPM(5),DVOL,NNO,NEF,IPT)

90	CONTINUE
C	----------------------------------------
	GOTO 1900
C	----------------------------------------
50    CONTINUE

CB    SAVE LOCAL DISPLACEMENT
   
      SDISL(IEL,1:24) = 0.0

C     ----------------------   
C     LOOP OVER GAUSS POINTS
C     ----------------------
      MGR=2
	MGS=2
      IPT = 0
      DO 900  IGR=1,MGR
      RI = GLOC(IGR,MGR)
      DO 900  IGS=1,MGS
      SI = GLOC(IGS,MGS)
      WT = GWT(IGR,MGR)*GWT(IGS,MGS)
      IPT = IPT+1
C     ------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P)
C     ------------------------------------
      CALL SHAP8NEW(RI,SI,SH,P)	
      
      SHR1 = SHR
      SHS1 = SHS
      SHT1 = SHT
      
C     ------------------------------------------------------------------
C     SETUP THE COROTATIONAL COORDINATE (VR,VS,VT)
C	JACOBIAN (XJ), INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ------------------------------------------------------------------	    
 	CALL COVAB(XY,PXY,SH,P,VR,VS,VT,XJ,XJI,
     @	SHR,SHS,SHT,FRS,GRS,DET,MEL,0)	     
      DVOL = WT*DET	
      
      CALL SHSTORL1_SHL (XY,PROPG,WA(1,IPT),NWG,VR,VS,VT)	!STORE SOLID-SHELL LAX FOR STRESS CALCULATION JAN09      
      CALL SHSTORL1_SHL (XY,PROPG,WA2(1,IPT),NWG,VR,VS,VT)	!STORE SOLID-SHELL LAX FOR STRESS CALCULATION JAN09       

	CALL MTEAS8(RI,SI,TTI,DETO,DET,GE,MM)

C     COMPOSITE MATIRIAL REF. VECTOR ANGLE RESPECTED TO R AXIS
      SELECTCASE(MTMOD)
      CASE(2,5,6)
	CALL CANGLE(VR,VS,AMV,ANG)
	ENDSELECT
	
C     ------------------------------------------------------------------
C     REMOVE RIGID BODY MOTION
C     ------------------------------------------------------------------
  	REDIS(1:24) = EDIS(1:24) 	
      IF(NLOPT.EQ.3) THEN	
        CALL SHAP3D8(RI,SI,0.0D0,H8,P8)
        CALL SOMDSH(COORD,COORDI,XY,XYI,REDIS,H8,P,VR,VS,VT)
      ENDIF


C     --------------------------------------------------------------
C     ASSUMED GLOBAL SHEAR STRAIN-DISPLACEMENT MATRIX (STORED IN BB)
C     --------------------------------------------------------------   
      BB=0.0D0  
      BBL=0.0D0 !BJ LOCAL 
      CALL SOBQ(RI,SI,XJI,XJA,XJB,XJC,XJD,PA,PB,PC,PD,BB,VRA,VSA,VTA,
     +		  VRB,VSB,VTB,VRC,VSC,VTC,VRD,VSD,VTD,SHA,SHB,SHC,SHD,BBL) 
	CALL SOBT(RI,SI,RA,SA,XJI,BB,BAT,BBL)
C     ----------------------------------------
C     GET GLOBAL STRAIN-DISPLACEMENT MATRIX BB
C     ----------------------------------------
      CALL BBMAT(VR,VS,VT,SH,SHR,SHS,SHT,XJI,BB,BBL)
	CALL WARP(VR,VS,VT,FRS,GRS,BB,BBL)

C     -------------------------------------
C     NEXT PART FOR NONLINEAR ANALYSIS ONLY
C     -------------------------------------

C     ----------------------------------------
C     TRANSFORMATIOIN MATRIX FOR DISPLACMEMENT
C     ----------------------------------------
C     DIMENSION T1(6,6),TLG(3,3),TLGT(6,6),TG1(6,6),TG2(24,24)  

C     T1 MATRIX            
      T1(1,1) = 0.5*1
      T1(2,2) = 0.5*1
      T1(3,3) = 0.5*1
      T1(1,4) = 0.5*1
      T1(2,5) = 0.5*1
      T1(3,6) = 0.5*1
      
      T1(4,1) = 0.5*1
      T1(5,2) = 0.5*1
      T1(6,3) = 0.5*1
      T1(4,4) = -0.5*1
      T1(5,5) = -0.5*1
      T1(6,6) = -0.5*1

C     TLG MATRIX

      TLG(1:3,1) = VR
      TLG(1:3,2) = VS
      TLG(1:3,3) = VT

C     TLGT MATRIX

      TLGT(1:3,1:3) = TRANSPOSE(TLG)
      TLGT(4:6,4:6) = TRANSPOSE(TLG)

C     TG1 MATRIX
      TG1 = MATMUL(T1,TLGT)
      
C     TG2 MATRIX
      
      DO I=1,4
        K = (I-1)*6
        KK = K+1
        TG2(KK:KK+5,KK:KK+5) = TG1 
      ENDDO      


C     ------------------
C     LOCAL DISPLACEMENT
C     ------------------
      DO I=1,8
      K=3*I-2
	DISL(K)  =VR(1)*REDIS(K)+VR(2)*REDIS(K+1)+VR(3)*REDIS(K+2)
	DISL(K+1)=VS(1)*REDIS(K)+VS(2)*REDIS(K+1)+VS(3)*REDIS(K+2)
	DISL(K+2)=VT(1)*REDIS(K)+VT(2)*REDIS(K+1)+VT(3)*REDIS(K+2)
	ENDDO

CBJ      DISL = MATMUL(TG2,REDIS)	
      
      DO I=1,3
	DO J=1,4
C	K=3*(J-1)+I
C	L=3*(J+3)+I
	K = IDISP(J,2*I-1)
	L = IDISP(J,2*I-0)

      U(I,J) =(DISL(K)+DISL(L))/2.
      PU(I,J)=(DISL(K)-DISL(L))/2.
CBJ       U(I,J) =(REDIS(K)+REDIS(L))/2.
CBJ       PU(I,J)=(REDIS(K)-REDIS(L))/2.
       
      M = 6*(J-1)
      N = I-1
 
      REDISS1(1+M+N) = REDIS(K)
      REDISS1(4+M+N) = REDIS(L)
      
	ENDDO
	ENDDO
	
CBJ      REDISS = MATMUL(TG2,REDISS1)
	
	DO I=1,4
	 K = 0
	  DO J=1,3
	      K = K + 1 + (6*(I-1))
	      REDISS(K) = U(J,I)
	      REDISS(K+3) = PU(J,I)
	  ENDDO
	ENDDO
      
	
	
	SDISL(IEL,1:24) = DISL(1:24)
	
C     ------------------------------
C     LOCAL DISPLACEMENT DERIVATIVES
C     ------------------------------
      DUD=0.
	DPUD=0.
	DO J=1,3
	DO I=1,4
       DUD(J)   = SHR(I)*U(J,I) + DUD(J)
	 DUD(J+3) = SHS(I)*U(J,I) + DUD(J+3)     
	 DUD(J+6) = SHT(I)*U(J,I) + DUD(J+6)     
	DPUD(J)   = SHR(I)*PU(J,I)+ DPUD(J)
	DPUD(J+3) = SHS(I)*PU(J,I)+ DPUD(J+3)
	DPUD(J+6) =  SH(I)*PU(J,I)+ DPUD(J+6)
	
CBJ	DPUD(J+6) = SHT(I)*PU(J,I)+ DPUD(J+6)
      ENDDO
	ENDDO
	
	DUD1=0.
	DPUD1=0.
	
	DO I=1,4
	  K = 1+(6*(I-1))
	  DUD1 = DUD1 + MATMUL(BB(1:9,K:K+2),REDISS(K:K+2))
	  DPUD1 = DUD1 + MATMUL(BB(1:9,K+3:K+5),REDISS(K+3:K+5)) 
	ENDDO

C     ---------------------
C     LINEAR ASSUMED STRAIN
C     ---------------------
      
      STRAIN=MATMUL(BB,REDIS)
      
C     ---------------
C     ALMANSI STRAINS
C     ---------------
!	ZT=XJI(3,3)
!	STRAIN(1)=DUD(1)
!	STRAIN(2)=DUD(5)
!	STRAIN(3)=DUD(4)+DUD(2)
!	STRAIN(4)=DPUD(1)
!	STRAIN(5)=DPUD(5)
!	STRAIN(6)=DPUD(4)+DPUD(2)
!	STRAIN(7)=ZT*DPUD(9)+DUD(9)
!	STRAIN(8)=ZT*DPUD(7)+DUD(3)+DUD(7)
!	STRAIN(9)=ZT*DPUD(8)+DUD(6)+DUD(8)
      
CBJ      ZT=XJ(3,3)
CBJ      STRAIN(1)=DUD(1)
CBJ      STRAIN(2)=DUD(5)
CBJ      STRAIN(3)=DUD(4)+DUD(2)      
CBJ      STRAIN(4)=ZT*DPUD(1)
CBJ      STRAIN(5)=ZT*DPUD(5)
CBJ      STRAIN(6)=ZT*(DPUD(4)+DPUD(2))
CBJ      STRAIN(7)=ZT*DPUD(9)+DUD(9)	
CBJ      STRAIN(8)=ZT*DPUD(7)+ZT*DPUD(3)+DUD(7)+DUD(3)	
CBJ      STRAIN(9)=ZT*DPUD(8)+ZT*DPUD(6)+DUD(8)+DUD(6)

CBJ      DO I=1,4
CBJ          K = (I-1)*6
CBJ          KK1= K+1
CBJ          KK2 = K+4
CBJ          STRAIN(1) =  STRAIN(1) + (BBL(1,KK1)*REDISS(KK1) + BBL(1,KK1+1)*REDISS(KK1+1))
CBJ          STRAIN(2) =  STRAIN(2) + (BBL(2,KK1)*REDISS(KK1) + BBL(2,KK1+1)*REDISS(KK1+1))
CBJ          STRAIN(3) =  STRAIN(3) + (BBL(3,KK1)*REDISS(KK1) + BBL(3,KK1+1)*REDISS(KK1+1))
        
CBJ          STRAIN(4) =  STRAIN(4) + (BBL(4,KK2)*REDISS(KK2) + BBL(4,KK2+1)*REDISS(KK2+1))
CBJ          STRAIN(5) =  STRAIN(5) + (BBL(5,KK2)*REDISS(KK2) + BBL(5,KK2+1)*REDISS(KK2+1))
CBJ          STRAIN(6) =  STRAIN(6) + (BBL(6,KK2)*REDISS(KK2) + BBL(6,KK2+1)*REDISS(KK2+1))
CBJ      ENDDO

CB      STRAINO = MATMUL(BBL(1:3,:),DISL)	     
CB      STRAINW = MATMUL((ZT*BBL(1:3,:)),DISL)	       
CB      STRAINB = MATMUL(BBL(4:6,:)+(ZT*BBL(4:6,:)),DISL)	   
CB      STRAIN(1:3) = STRAINO + STRAINW
CB      STRAIN(4:6) = STRAINB

CBJ      DO I=1,8
CBJ      K=3*I-2
CBJ	STRAIN(K)  =VR(1)*STRAIN(K)+VR(2)*STRAIN(K+1)+VR(3)*STRAIN(K+2)
CBJ	STRAIN(K+1)=VS(1)*STRAIN(K)+VS(2)*STRAIN(K+1)+VS(3)*STRAIN(K+2)
CBJ	STRAIN(K+2)=VT(1)*STRAIN(K)+VT(2)*STRAIN(K+1)+VT(3)*STRAIN(K+2)
CBJ	ENDDO 


CBJ   STRAIN = MATMUL(BB,REDISS)
C     ---------------------
C     SUBTRACT EAS STRAIN
C     ---------------------
	EAS = MATMUL(GE,ALPHA)
	DO I=1,3
      STRAIN(I)= STRAIN(I) - EAS(I)
	ENDDO

      
C     ----------------------
C     QUADRATIC STRAIN TERMS
C     ----------------------
C     ----MEMBRANE PART-----
      QSTRAI(1)=.5*(DUD(1)*DUD(1)+DUD(2)*DUD(2)+DUD(3)*DUD(3))
	QSTRAI(2)=.5*(DUD(4)*DUD(4)+DUD(5)*DUD(5)+DUD(6)*DUD(6))
	QSTRAI(3)=    DUD(1)*DUD(4)+DUD(2)*DUD(5)+DUD(3)*DUD(6)
C     ----BENDING PART------
      QSTRAI(4)= DUD(1)*DPUD(1)+DUD(2)*DPUD(2)+DUD(3)*DPUD(3)
      QSTRAI(5)= DUD(4)*DPUD(4)+DUD(5)*DPUD(5)+DUD(6)*DPUD(6)
	QSTRAI(6)= DPUD(1)*DUD(4)+DUD(1)*DPUD(4)+DPUD(2)*DUD(5)
     @          +DUD(2)*DPUD(5)+DPUD(3)*DUD(6)+DUD(3)*DPUD(6)
C     ----THICKNESS PART----
      ZT=XJI(3,3)
      QSTRAI(7)=.5*(DUD(7)*DUD(7)+DPUD(7)*DPUD(7)*ZT*ZT+DUD(8)*DUD(8)
     @     +DPUD(8)*DPUD(8)*ZT*ZT+DUD(9)*DUD(9)+DPUD(9)*DPUD(9)*ZT*ZT)
C     ----SHEAR PART--------
      QSTRAI(8)= DUD(1)*(DUD(7)+DPUD(7)*ZT)+DUD(2)*(DUD(8)+DPUD(8)*ZT)
     @           +DUD(3)*(DUD(9)+DPUD(9)*ZT)
	QSTRAI(9)= DUD(4)*(DUD(7)+DPUD(7)*ZT)+DUD(5)*(DUD(8)+DPUD(8)*ZT)
     @           +DUD(6)*(DUD(9)+DPUD(9)*ZT)

C
	IF(NLOPT.EQ.3) THEN
	DO II = 1,9
	STRAIN(II)=STRAIN(II) - QSTRAI(II)
	ENDDO
	ENDIF
C     -----------------------
C     STRESS IN DEFORMED BODY & GET LINEAR STIFFNESS MATRIX
C     -----------------------
      SELECTCASE(MTMOD)
      CASE(1)
C     LINEAR ELASTIC  


      !STRESS=MATMUL(DE,STRAIN)      
      STRESSO = MATMUL(DEO2,STRAIN) !RESULTANT STRESS
      STRESS = MATMUL(DEO,STRAIN) !STRESS
                  	
C	------------------------------------------------------------------
C	COMPUTE FIBER STRESS 
C	------------------------------------------------------------------
	!ZETAG = 0.5*THK 
      !STRESS_FIBER = STRESS
      !STRESS_FIBER(1)  = STRESS(1) + ZETAG*STRESS(4)
	!STRESS_FIBER(2)  = STRESS(2) + ZETAG*STRESS(5)
	!STRESS_FIBER(4)  = STRESS(3) + ZETAG*STRESS(6)
	!STRESS_FIBER(3)  = 0.0D0
	!STRESS_FIBER(5)  = STRESS(8)
	!STRESS_FIBER(6)  = STRESS(9)
      
      !STRESS_FIBER(7)  = STRESS(7)!0.0D0
      !STRESS_FIBER(8)  = STRESS(8)!0.0D0
      !STRESS_FIBER(9)  = STRESS(9)!0.0D0
C     ---------------------------------------------------------------------      
      
      !THK3 = THK*THK*THK/12.
      !STRESSO(1) = STRESS(1)*THK
      !STRESSO(2) = STRESS(2)*THK
      !STRESSO(3) = STRESS(3)*THK
      !STRESSO(4) = STRESS(4)*THK3
      !STRESSO(5) = STRESS(5)*THK3
      !STRESSO(6) = STRESS(6)*THK3
      !STRESSO(7) = STRESS(7)*THK
      !STRESSO(8) = STRESS(8)*THK      
      !STRESSO(9) = STRESS(9)*THK      

      !DO I=1,8
      !    K=3*I-2      
      !    STRESSL(K)  =VR(1)*STRESSO(K)+VS(1)*STRESSO(K+1)+VT(1)*STRESSO(K+2)
	!    STRESSL(K+1)=VR(2)*STRESSO(K)+VS(2)*STRESSO(K+1)+VT(2)*STRESSO(K+2)
	!    STRESSL(K+2)=VR(3)*STRESSO(K)+VS(3)*STRESSO(K+1)+VT(3)*STRESSO(K+2)
      !ENDDO      
CB      STRESSO(1:3) = MATMUL(DEO(1:3,1:3),STRAINO)
CB      STRESSO(4:6) = MATMUL(DEO(4:6,4:6),(STRAINB+STRAINW))

      !IF(IEL.EQ.16) THEN
      !WRITE(*,*) "IEL", IEL
      !WRITE(*,*) "IPT", IPT  
      !WRITE(*,*) "STRAIN"
      !WRITE(*,*) STRAIN(4:6) 
      !WRITE(*,*) "STRESS0"
      !WRITE(*,*) STRESSO(4:6) 
      !ENDIF      
 
      
            DO I=1,9
          WA(I,IPT) = STRESSO(I) !FOR RESULTANT STRESS
          !WA2(I,IPT) = STRESS_FIBER(I) !FOR STRESS         
          WA2(I,IPT) = STRESS(I) !FOR STRESS         
      ENDDO
      
      CASE(2)
C     LAMINATE COMPOSITE

      THK = 2.0D0 !NORMALIZE THICKNESS FOR COMPOSITE

      CALL COMRGD(STRAIN,STRESS,DE,PROPM,THK,ANG,DVOL,'SOLID')
      DO I=1,9
	    WA(I  ,IPT) = STRAIN(I)
          WA(I+9,IPT) = STRESS(I)
      ENDDO
      
      CASE(4)  	
C	ELASTO PLASTIC DRACKER
	NLAYR = NGT
	CALL SOHLET(PROPM,WA(1,IPT),STRAIN,STRESS,NLAYR,DE)
	
      CASE(5,6)
C	RC EPF MODEL
	CALL SRCLER(STRAIN,STRESS,WA(1,IPT),DE,MSET,MTMOD,ANG)
	
	ENDSELECT
      
      IF(IFREF.EQ.0) THEN
	  CALL SOKNEW(S,BB,DE,DVOL)
	ENDIF
	
C     ---------------------
C     INTERNAL FORCE VECTOR
C     ---------------------
      TAU=STRESS*DVOL


	RE = RE + MATMUL(TRANSPOSE(BB),TAU)	
	
CCC	RE = RE + MATMUL(TRANSPOSE(BB(4:6,:)),TAU(4:6))	

	RH = RH + MATMUL(TRANSPOSE(GE),TAU)


C     ----------------------------------------------------------------
C     TAKE INTO ACCOUNT EXTRA STRAIN TERMS IN SHEAR AND THICKNESS
C     ----------------------------------------------------------------
C	CALL BSMAT(VR,VS,VT,SH,SHR,SHS,SHT,XJI,S,DVOL,DE,RE,TAU)

      IF(IFEIG.EQ.0) GOTO 750
      
      IF(NLOPT.NE.3) GOTO 800 
      IF(IFREF.EQ.1) GOTO 800
      
750   CALL SOKG(S,SH,SHR,SHS,SHT,VR,VS,VT,TAU,ZT)

800	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(DE,GE))*DVOL
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(DE,BB))*DVOL


900   CONTINUE 


CDIS      DEALLOCATE(SDISL)
C     ----------------------------------------------------------------
	DO I = 1,MM
	IF(NLOPT.GT.0) HINFC(I) = RH(I)
	ENDDO
	 
	CALL INVMATRIX(SED,SEDI,MM)

	IF(NLOPT.GT.0) REAS = MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,RH))
	SEAS = MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,SEL))

	K = 0
	DO I = 1,24
	IF(NLOPT.GT.0) RE(I) = RE(I) - REAS(I) 
	DO J = I,24
	K = K + 1
	S(K) = S(K) - SEAS(I,J) 
	ENDDO
	ENDDO 
C     ----------------------------------------------------------------
1900	CONTINUE

C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF


      RETURN
      END
C==================================================================
      SUBROUTINE SLFACE(NMX)                 !BLOCK DATA DISPLACEMETS
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ============================================
C     DEFINE THE FACES FOR THE 8 NODE SOLID ELEMET
C     ============================================
      COMMON /DISPL/ IDISP(4,6)

	IDISP = 0

	IF(NMX.EQ.1) THEN
	IDISP(1,1) = 1
	IDISP(1,2) = 4
	IDISP(1,3) = 2
	IDISP(1,4) = 5
	IDISP(1,5) = 3
	IDISP(1,6) = 6
	IDISP(2,1) = 10
	IDISP(2,2) = 7
	IDISP(2,3) = 11
	IDISP(2,4) = 8
	IDISP(2,5) = 12
	IDISP(2,6) = 9
	IDISP(3,1) = 22
	IDISP(3,2) = 19
	IDISP(3,3) = 23
	IDISP(3,4) = 20
	IDISP(3,5) = 24
	IDISP(3,6) = 21
	IDISP(4,1) = 13
	IDISP(4,2) = 16
	IDISP(4,3) = 14
	IDISP(4,4) = 17
	IDISP(4,5) = 15
	IDISP(4,6) = 18
	ELSEIF(NMX.EQ.2) THEN
	IDISP(1,1) = 1
	IDISP(1,2) = 10
	IDISP(1,3) = 2
	IDISP(1,4) = 11
	IDISP(1,5) = 3
	IDISP(1,6) = 12
	IDISP(2,1) = 13
	IDISP(2,2) = 22
	IDISP(2,3) = 14
	IDISP(2,4) = 23
	IDISP(2,5) = 15
	IDISP(2,6) = 24
	IDISP(3,1) = 16
	IDISP(3,2) = 19
	IDISP(3,3) = 17
	IDISP(3,4) = 20
	IDISP(3,5) = 18
	IDISP(3,6) = 21
	IDISP(4,1) = 4
	IDISP(4,2) = 7
	IDISP(4,3) = 5
	IDISP(4,4) = 8
	IDISP(4,5) = 6
	IDISP(4,6) = 9
	ELSEIF(NMX.EQ.3) THEN
	IDISP(1,1) = 1
	IDISP(1,2) = 13
	IDISP(1,3) = 2
	IDISP(1,4) = 14
	IDISP(1,5) = 3
	IDISP(1,6) = 15
	IDISP(2,1) = 4
	IDISP(2,2) = 16
	IDISP(2,3) = 5
	IDISP(2,4) = 17
	IDISP(2,5) = 6
	IDISP(2,6) = 18
	IDISP(3,1) = 7
	IDISP(3,2) = 19
	IDISP(3,3) = 8
	IDISP(3,4) = 20
	IDISP(3,5) = 9
	IDISP(3,6) = 21
	IDISP(4,1) = 10
	IDISP(4,2) = 22
	IDISP(4,3) = 11
	IDISP(4,4) = 23
	IDISP(4,5) = 12
	IDISP(4,6) = 24
	ENDIF

	RETURN
	END
C*****************************************************************************
      SUBROUTINE COVAB(XY,PXY,SH,P,VR,VS,VT,
     +        	     XJ,XJI,SHR,SHS,SHT,FRS,GRS,DET,MEL,ISAMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------
C     GET BASE VECTORS, JACOBEAN MATRIX AND ITS INVERSE
C     -------------------------------------------------
      DIMENSION VR(3),VS(3),VT(3),COVR(3),COVS(3)
	DIMENSION COORD(3,8),SH(4),P(2,4),XJ(3,3),XJI(3,3)
	DIMENSION XY(3,4),PXY(3,4),PP(3)
	DIMENSION SHR(4),SHS(4),SHT(4),DPR(3),DPS(3)
	DIMENSION FRS(4),GRS(4)
	DIMENSION PR(4),PS(4),PHI(2,2)
	
c	DIMENSION RV(3),SV(3)
C
	PP=0.
	DO I=1,3
	DO J=1,4
      PP(I)=PXY(I,J)*SH(J)+PP(I)
	ENDDO
	ENDDO
C     ---------------------------------------------
C     DERIVATIVES OF VECTOR P (WITH RESPECT TO R,S)
C     ---------------------------------------------
      DPR=0.
	DPS=0.
	DO I=1,3
	DO J=1,4
	DPR(I)=DPR(I)+PXY(I,J)*P(1,J)
	DPS(I)=DPS(I)+PXY(I,J)*P(2,J)
	ENDDO
	ENDDO
C     -----------------------------------
C     TANGENTIAL DIRECTIONS TO MIDSURFACE
C     -----------------------------------
      COVR=0.
	COVS=0.
      DO I=1,3
	DO J=1,4
	COVR(I)=COVR(I)+XY(I,J)*P(1,J)
	COVS(I)=COVS(I)+XY(I,J)*P(2,J)
	ENDDO
	ENDDO	
C     -----------------------------
C     LOCAL COORDINATE BASE VECTORS
C     -----------------------------
      CALL VECPRD (COVR,COVS,VT)
	CALL SCALEN (VT,VT,DD,3)
	CALL VECPRD (VT,COVR,VS)
	CALL SCALEN (VS,VS,DD,3)
	CALL SCALEN (COVR,VR,DD,3)	
	
C      CALL VECPRD (COVR,COVS,VT)
C      CALL SCALEN (VT,VT,DD,3)
C      CALL SCALEN (COVR,RV,RL,3)
C      CALL SCALEN (COVS,SV,SL,3)
C      CALL VECPRD (VT,RV,VS)
C      CALL ADDVEC (SV,VS,VS)
C      CALL SCALEN (VS,VS,DD,3)
C      CALL VECPRD (VS,VT,VR)
C      CALL SCALEN (VR,VR,DD,3)
C      CALL SCALEN (VS,VS,DD,3)
C      CALL SCALEN (VT,VT,DD,3)
      	
C     --------------------------------------
C     JACOBIAN MATRIX COMPUTED IN MIDSURFACE
C     --------------------------------------
      XJ(1,1)=DOT_PRODUCT(COVR,VR)
	XJ(1,2)=DOT_PRODUCT(COVR,VS)
	XJ(2,1)=DOT_PRODUCT(COVS,VR)
	XJ(2,2)=DOT_PRODUCT(COVS,VS)
	XJ(1,3)=0.0D0
	XJ(2,3)=0.0D0
	XJ(3,1)=DOT_PRODUCT(VR,PP)
	XJ(3,2)=DOT_PRODUCT(VS,PP)
	XJ(3,3)=DOT_PRODUCT(VT,PP)
C     -------------------------------------------------------
C     FLAG TO COMPUTE JACOBEAN MATRIX ONLY AT SAMPLING POINTS
C     -------------------------------------------------------
      IF(ISAMP.EQ.1)RETURN
C     ------------------------------	
C     DETERMINANT OF JACOBEAN MATRIX
C     ------------------------------
      DET0 = XJ(1,2)*XJ(2,1)-XJ(1,1)*XJ(2,2)
	DET  =-DET0*XJ(3,3)

	IF(DET.LT.1.E-10) CALL ERRORS (15,H,MEL,'JACOB.DET.')
C     --------------------------
C     INVERSE OF JACOBEAN MATRIX
C     -------------------------- 
      XJI(1,3)=0.
	XJI(2,3)=0.
	XJI(3,1)=(-XJ(2,2)*XJ(3,1) + XJ(2,1)*XJ(3,2))/DET
	XJI(3,2)=( XJ(1,2)*XJ(3,1) - XJ(1,1)*XJ(3,2))/DET
	XJI(1,1)=-XJ(2,2)/DET0
	XJI(1,2)= XJ(1,2)/DET0
	XJI(2,1)= XJ(2,1)/DET0
	XJI(2,2)=-XJ(1,1)/DET0
	XJI(3,3)= 1.0D0/XJ(3,3)

C     ---------------------------------------------------------------
C     SHAPE FUNCTION DERIVATIVES WITH RESPECT TO LOCAL COORDINATE R,S
C     ---------------------------------------------------------------
      DO I=1,4
	SHR(I)=XJI(1,1)*P(1,I)+XJI(1,2)*P(2,I)
      SHS(I)=XJI(2,1)*P(1,I)+XJI(2,2)*P(2,I)
	SHT(I)=XJI(3,1)*P(1,I)+XJI(3,2)*P(2,I)
	ENDDO
C     ----------------------------------------------------
C     NEXT COMPUTE THE AFFECT DUE TO WARPED ELEMENT SHAPES
C     ----------------------------------------------------
      XLAM1=-DET0
	XLAM2=DOT_PRODUCT(VR,COVR)*DOT_PRODUCT(VS,DPS)
	XLAM2=XLAM2+DOT_PRODUCT(VR,DPR)*DOT_PRODUCT(VS,COVS)
	XLAM2=XLAM2-DOT_PRODUCT(VR,COVS)*DOT_PRODUCT(VS,DPR)
	XLAM2=XLAM2-DOT_PRODUCT(VR,DPS)*DOT_PRODUCT(VS,COVR)
	BETA1=1./XLAM1
	BETA2=-XLAM2/(XLAM1*XLAM1)
C
      PHI(1,1)= BETA1*DOT_PRODUCT(VS,DPS)+BETA2*DOT_PRODUCT(VS,COVS)
	PHI(1,2)=-BETA1*DOT_PRODUCT(VS,DPR)-BETA2*DOT_PRODUCT(VS,COVR)
	PHI(2,1)=-BETA1*DOT_PRODUCT(VR,DPS)-BETA2*DOT_PRODUCT(VR,COVS)
	PHI(2,2)= BETA1*DOT_PRODUCT(VR,DPR)+BETA2*DOT_PRODUCT(VR,COVR)
C 
      DO I=1,4
      FRS(I)=PHI(1,1)*P(1,I)+PHI(1,2)*P(2,I)
      GRS(I)=PHI(2,1)*P(1,I)+PHI(2,2)*P(2,I)
	ENDDO
	
C
	RETURN
	END
C**************************************************************************
      SUBROUTINE SHAP8NEWP1(R,S,SH1)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     AT THE NODAL POINTS OF A CURVILINEAR ISOPARAMETRIC HEXAHEDRON
C	--------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     R,S       = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     SH(4)     = INTERPOLATION (SHAPE) FUNCTIONS
C     P(2,4)    = FUNCTION DERIVATIVES WITH RESPECT TO R,S,T RESP.
C     ----------------------------------------------------------------
      DIMENSION SH1(4)!,SH2(8)
C
      RP  = 1.+R
      SP  = 1.+S
C      TP  = (1.-T)
      
      RM  = 1.-R
      SM  = 1.-S
C      TM  = (1.-T)
      
C     ----------------------------------------------------------
C     INTERPOLATION FUNCTIONS AND DERIVATIVES FOR A 8 NODE BRICK
C     ----------------------------------------------------------
      SH1(1)   = .25*RP*SP
      SH1(2)   = .25*RM*SP
      SH1(3)   = .25*RM*SM
      SH1(4)   = .25*RP*SM
      

C      SH2(1)   = TP
C      SH2(2)   = -TP
C      SH2(3)   = TP
C      SH2(4)   = -TP      
C      SH2(5)   = TP
C      SH2(6)   = -TP
C      SH2(7)   = TP
C      SH2(8)   = -TP

      RETURN
	END
C**************************************************************************	
      SUBROUTINE SHAP8NEWP(R,S,T,SH)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     AT THE NODAL POINTS OF A CURVILINEAR ISOPARAMETRIC HEXAHEDRON
C	--------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     R,S       = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     SH(4)     = INTERPOLATION (SHAPE) FUNCTIONS
C     P(2,4)    = FUNCTION DERIVATIVES WITH RESPECT TO R,S,T RESP.
C     ----------------------------------------------------------------
      DIMENSION SH(8)
C
      RP  = 1.+R
      SP  = 1.+S
      TP  = (1.-T)
      
      RM  = 1.-R
      SM  = 1.-S
      TM  = (1.-T)
      
C     ----------------------------------------------------------
C     INTERPOLATION FUNCTIONS AND DERIVATIVES FOR A 8 NODE BRICK
C     ----------------------------------------------------------
      SH(1)   = .25*RP*SP*TP
      SH(2)   = .25*RM*SP*TP
      SH(3)   = .25*RM*SM*TP
      SH(4)   = .25*RP*SM*TP
      SH(5)   = .25*RP*SP*TM
      SH(6)   = .25*RM*SP*TM
      SH(7)   = .25*RM*SM*TM
      SH(8)   = .25*RP*SM*TM

      RETURN
	END
C**************************************************************************
      SUBROUTINE SHAP8NEW(R,S,SH,P)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     AT THE NODAL POINTS OF A CURVILINEAR ISOPARAMETRIC HEXAHEDRON
C	--------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     R,S       = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     SH(4)     = INTERPOLATION (SHAPE) FUNCTIONS
C     P(2,4)    = FUNCTION DERIVATIVES WITH RESPECT TO R,S,T RESP.
C     ----------------------------------------------------------------
      DIMENSION P(2,4),SH(4)
C
      RP  = 1.+R
      SP  = 1.+S
      RM  = 1.-R
      SM  = 1.-S
C     ----------------------------------------------------------
C     INTERPOLATION FUNCTIONS AND DERIVATIVES FOR A 8 NODE BRICK
C     ----------------------------------------------------------
	SH(1)   = .25*RP*SP
      SH(2)   = .25*RM*SP
      SH(3)   = .25*RM*SM
      SH(4)   = .25*RP*SM 
C
      P(1,1) = .25*SP  
      P(1,2) =-.25*SP  
      P(1,3) =-.25*SM  
      P(1,4) = .25*SM   
C
      P(2,1) = .25*RP 
      P(2,2) = .25*RM 
      P(2,3) =-.25*RM  
      P(2,4) =-.25*RP  
C
      RETURN
	END
	
	
	
C*****************************************************************************
      SUBROUTINE SOBQ(RI,SI,XJI,XJA,XJB,XJC,XJD,PA,PB,PC,PD,BB,
     +VRA,VSA,VTA,VRB,VSB,VTB,VRC,VSC,VTC,VRD,VSD,VTD,SHA,SHB,SHC,SHD,BBL) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     CALCULATE THE ASSUMED SHEAR STRAIN-DISPLACEMENT MATRIX (BQ)
C     THE CORRESPONDING DISPLACEMENTS TO BQ IS IN GLOBAL COORDINATE
C     -------------------------------------------------------------
      COMMON /DISPL/ IDISP(4,6)
      DIMENSION XJI(3,3),XJA(3,3),XJB(3,3),XJC(3,3),XJD(3,3)
	DIMENSION BB(9,24),PA(4),PB(4),PC(4),PD(4),BBL(9,24)
	DIMENSION VRA(3),VSA(3),VTA(3)
      DIMENSION VRB(3),VSB(3),VTB(3)
	DIMENSION VRC(3),VSC(3),VTC(3)
	DIMENSION VRD(3),VSD(3),VTD(3)
	DIMENSION SHA(4),SHB(4),SHC(4),SHD(4)
C
      DIMENSION T(2,2)
	DIMENSION WORKA(1,3),WORKB(1,3),WORKC(1,3),WORKD(1,3)
	DIMENSION W(2,3),W1(3,3),W2(2,6),WL(2,3),WL2(2,6)
C     -------------------------------------------
C     TRANSFORMATION MATRIX FOR SHEAR STRAINS (T)
C     -------------------------------------------
      T(1,1)=XJI(1,1)*XJI(3,3)
	T(2,2)=XJI(2,2)*XJI(3,3)
	T(1,2)=XJI(1,2)*XJI(3,3)
	T(2,1)=XJI(2,1)*XJI(3,3) 


	T(1,1)=XJI(1,1)*XJI(3,3)
	T(2,2)=XJI(2,2)*XJI(3,3)
	T(1,2)=XJI(2,1)*XJI(3,3)
	T(2,1)=XJI(1,2)*XJI(3,3) 
	
C     ------PART 1------
      DO I=1,4
      WORKB(1,1)=0.
	WORKB(1,2)=0.
	WORKB(1,3)=PB(I)
      WORKD(1,1)=0.
	WORKD(1,2)=0.
	WORKD(1,3)=PD(I)
C
	DO J=1,3
	W1(1,J)=VRB(J)
	W1(2,J)=VSB(J)
	W1(3,J)=VTB(J)
	ENDDO
	WORKB=MATMUL(WORKB,XJB)	
	WORKB=MATMUL(WORKB,W1)
C
	DO J=1,3
	W1(1,J)=VRD(J)
	W1(2,J)=VSD(J)
	W1(3,J)=VTD(J)
	ENDDO
	WORKD=MATMUL(WORKD,XJD)
	WORKD=MATMUL(WORKD,W1)
C
      WORKA(1,1)=0.
	WORKA(1,2)=0.
	WORKA(1,3)=PA(I)
	WORKC(1,1)=0.
	WORKC(1,2)=0.
	WORKC(1,3)=PC(I)
C
      DO J=1,3
	W1(1,J)=VRA(J)
	W1(2,J)=VSA(J)
	W1(3,J)=VTA(J)
	ENDDO
	WORKA=MATMUL(WORKA,XJA)
	WORKA=MATMUL(WORKA,W1)
C
      DO J=1,3
	W1(1,J)=VRC(J)
	W1(2,J)=VSC(J)
	W1(3,J)=VTC(J)
	ENDDO
	WORKC=MATMUL(WORKC,XJC)
	WORKC=MATMUL(WORKC,W1)
C
      DO J=1,3
      W(1,J)=((1.-SI)*WORKB(1,J)+(1.+SI)*WORKD(1,J))*.5
	W(2,J)=((1.+RI)*WORKA(1,J)+(1.-RI)*WORKC(1,J))*.5
	ENDDO	
	
	WL = W ! BJ LOCAL
	W = MATMUL(T,W)	
C
	DO J=1,3
	W2(1,2*J-1)=W(1,J)*.5
	W2(2,2*J-1)=W(2,J)*.5
	W2(1,2*J)=W(1,J)*.5
	W2(2,2*J)=W(2,J)*.5
	ENDDO
	
	DO J=1,3 ! BJ LOCAL
	WL2(1,2*J-1)=WL(1,J)*.5
	WL2(2,2*J-1)=WL(2,J)*.5
	WL2(1,2*J)=WL(1,J)*.5
	WL2(2,2*J)=WL(2,J)*.5
	ENDDO
	
C	
	DO M=1,6
      N=IDISP(I,M)
	DO J=1,2
	BB(J+7,N)=W2(J,M)
      ENDDO
	ENDDO

	DO M=1,6 ! BJ LOCAL
      N=IDISP(I,M)
	DO J=1,2
	BBL(J+7,N)=WL2(J,M)
      ENDDO
	ENDDO
	
C
 	ENDDO
 	
C     -------PART 2--------
      DO I=1,4
	WORKB(1,2)=0.
	WORKB(1,3)=0.
	WORKD(1,2)=0.
	WORKD(1,3)=0.
	WORKB(1,1)=SHB(I)
	WORKD(1,1)=SHD(I)
	DO J=1,3
	W1(1,J)=VRB(J)
	W1(2,J)=VSB(J)
	W1(3,J)=VTB(J)
	ENDDO
	WORKB=MATMUL(WORKB,XJB)	
	WORKB=MATMUL(WORKB,W1)
C
	DO J=1,3
	W1(1,J)=VRD(J)
	W1(2,J)=VSD(J)
	W1(3,J)=VTD(J)
	ENDDO
	WORKD=MATMUL(WORKD,XJD)
	WORKD=MATMUL(WORKD,W1)
C
      WORKA(1,1)=0.
	WORKC(1,1)=0.
	WORKA(1,3)=0.
	WORKC(1,3)=0.
	WORKA(1,2)=SHA(I)
	WORKC(1,2)=SHC(I)
      DO J=1,3
	W1(1,J)=VRA(J)
	W1(2,J)=VSA(J)
	W1(3,J)=VTA(J)
	ENDDO
	WORKA=MATMUL(WORKA,XJA)
	WORKA=MATMUL(WORKA,W1)
C
      DO J=1,3
	W1(1,J)=VRC(J)
	W1(2,J)=VSC(J)
	W1(3,J)=VTC(J)
	ENDDO
	WORKC=MATMUL(WORKC,XJC)
	WORKC=MATMUL(WORKC,W1)
C
      DO J=1,3
      W(1,J)=((1.-SI)*WORKB(1,J)+(1.+SI)*WORKD(1,J))*.5
	W(2,J)=((1.+RI)*WORKA(1,J)+(1.-RI)*WORKC(1,J))*.5
	ENDDO	
	
	WL = W ! BJ LOCAL
	W = MATMUL(T,W)	
		 
	DO J=1,3
	W2(1,2*J-1)=W(1,J)*.5 
	W2(2,2*J-1)=W(2,J)*.5 
	W2(1,2*J)=-W(1,J)*.5 
	W2(2,2*J)=-W(2,J)*.5 
	ENDDO

	DO J=1,3 ! BJ LOCAL
	WL2(1,2*J-1)=WL(1,J)*.5 
	WL2(2,2*J-1)=WL(2,J)*.5 
	WL2(1,2*J)=-WL(1,J)*.5 
	WL2(2,2*J)=-WL(2,J)*.5 
	ENDDO

C
	DO M=1,6
      N=IDISP(I,M)
	DO J=1,2
	BB(J+7,N)=W2(J,M)+BB(J+7,N)
      ENDDO
	ENDDO

	DO M=1,6 ! BJ LOCAL
      N=IDISP(I,M)
	DO J=1,2
	BBL(J+7,N)=WL2(J,M)+BBL(J+7,N)
      ENDDO
	ENDDO


      ENDDO


      RETURN
	END
C*****************************************************************************
      SUBROUTINE BBMAT(VR,VS,VT,SH,SHR,SHS,SHT,XJI,BB,BBL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------
C     CALCULATE GLOBAL STRAIN-DISPLACEMENT MATRIX BB
C     ----------------------------------------------
      COMMON /DISPL/ IDISP(4,6)
      DIMENSION VR(3),VS(3),VT(3)
	DIMENSION SH(4),SHR(4),SHS(4),SHT(4)
	DIMENSION XJI(3,3),BB(9,24)
C
      DIMENSION W(3,3),W1(3,3),WORK(3,6),WW(1,3)
      DIMENSION BBL(9,24),WORKL(3,6),WL(3,3)
C
	DO J=1,3
	W1(1,J)=VR(J)
	W1(2,J)=VS(J)
	W1(3,J)=VT(J)
	ENDDO
	
      DO I=1,4
	W=0.
      W(1,1)=SHR(I)
	W(2,2)=SHS(I)
	W(3,1)=SHS(I)
	W(3,2)=SHR(I)
	
	
      WL = W ! BJ LOCAL *********
	W=MATMUL(W,W1)	
		
	
C     ----MEMBRANE PART-------
      DO J=1,3
	DO K=1,3
      WORK(K,2*J)=W(K,J)*.5
	WORK(K,2*J-1)=W(K,J)*.5
	ENDDO
	ENDDO
	
      DO M=1,6
	N=IDISP(I,M)
	DO J=1,3
      BB(J,N)=WORK(J,M)
	ENDDO
	ENDDO

      DO J=1,3 ! BJ LOCAL
	DO K=1,3
      WORKL(K,2*J)=WL(K,J)*.5
	WORKL(K,2*J-1)=WL(K,J)*.5
	ENDDO
	ENDDO
	
      DO M=1,6
	N=IDISP(I,M)
	DO J=1,3
      BBL(J,N)=WORKL(J,M)
	ENDDO
	ENDDO

	
C     ----BENDING PART-------
      DO J=1,3
	DO K=1,3
      WORK(K,2*J)=-W(K,J)*.5
	WORK(K,2*J-1)=W(K,J)*.5
	ENDDO
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	DO J=1,3
      BB(J+3,N)=WORK(J,M) 
	ENDDO
	ENDDO
	
	DO J=1,3 ! BJ LOCAL
	DO K=1,3
      WORKL(K,2*J)=-WL(K,J)*.5
	WORKL(K,2*J-1)=WL(K,J)*.5
	ENDDO
	
	ENDDO ! BJ LOCAL
      DO M=1,6
	N=IDISP(I,M)
	DO J=1,3
      BBL(J+3,N)=WORKL(J,M) 
	ENDDO
	ENDDO
	
C     ------------------------------------------------------
C     ----SHEAR PART ALREADY DONE IN SUBROUTINE SOBQ--------     
C     ----THICKNESS PART ALREADY DONE IN SUBROUTINE SOBT----
C     ------------------------------------------------------
	ENDDO
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE DEMAT(DE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------
C     GET INTEGRATED MATERIAL MATRIX FOR ISOTROPIC MATERIALS
C     ------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      
CB    THE THICKNESS
      COMMON /THICK/THK
            
	DIMENSION DE(9,9),DMB(3,3)
C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     ------------------------------------------------------------	

      
      SMU=YM/(2.*(1.+PR))
      SLA=YM*PR/(1.-PR*PR)
C
      DMB(1,1)=SLA+2.*SMU
	DMB(1,2)=SLA
	DMB(2,1)=SLA
	DMB(2,2)=SLA+2.*SMU
	DMB(3,3)=SMU
	DMB(1,3)=0.
	DMB(2,3)=0.
	DMB(3,1)=0.
	DMB(3,2)=0.
C      
		
	RETURN
	END
C***************************************************************************************      
      SUBROUTINE DEMAT_SHL(DE,DEO,DEO2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------
C     GET INTEGRATED MATERIAL MATRIX FOR ISOTROPIC MATERIALS
C     ------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      
CB    THE THICKNESS
      COMMON /THICK/THK
            
	DIMENSION DE(9,9),DMB(3,3),DMBO(9,9),DEO(9,9),DMBO2(9,9),DEO2(9,9)
C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     ------------------------------------------------------------	

      
      SMU=YM/(2.*(1.+PR))
      SLA=YM*PR/(1.-PR*PR)
C
      DMB(1,1)=SLA+2.*SMU
	DMB(1,2)=SLA
	DMB(2,1)=SLA
	DMB(2,2)=SLA+2.*SMU
	DMB(3,3)=SMU
	DMB(1,3)=0.
	DMB(2,3)=0.
	DMB(3,1)=0.
	DMB(3,2)=0.
C
      DE=0.
	DO I=1,3
	DO J=1,3
	DE(I,J)=2.*DMB(I,J)
	DE(I+3,J+3)=2./3.*DMB(I,J)
	ENDDO
	ENDDO
	DE(7,7)=2.*YM
	DE(8,8)=2.*SMU*5./6.
	DE(9,9)=2.*SMU*5./6.

CC    REGIDITY MATRIX BEFORE INTEGRATION (FOR STRESS (LIKE FIBER STRESS OF SHELL))	
      DMBO(1,1)=SLA+2.*SMU
	DMBO(1,2)=SLA
	DMBO(2,1)=SLA
	DMBO(2,2)=SLA+2.*SMU
	DMBO(3,3)=SMU
	DMBO(1,3)=0.
	DMBO(2,3)=0.
	DMBO(3,1)=0.
	DMBO(3,2)=0.
	
      DEO = 0.0
	DO I=1,3
	DO J=1,3
      DEO(I,J)=DMBO(I,J)            
      DEO(I+3,J+3)=DMBO(I,J)
	ENDDO
	ENDDO
	
	DEO(7,7)=YM
	DEO(8,8)=SMU*(5.0/6.0) !SHEAR COEFFICIENT (5/6)
	DEO(9,9)=SMU*(5.0/6.0) !SHEAR COEFFICIENT (5/6)  
      
CC    REGIDITY MATRIX BEFORE INTEGRATION (FOR RESULTANT STRESS)	
      DMBO2(1,1)=SLA+2.*SMU
	DMBO2(1,2)=SLA
	DMBO2(2,1)=SLA
	DMBO2(2,2)=SLA+2.*SMU
	DMBO2(3,3)=SMU
	DMBO2(1,3)=0.
	DMBO2(2,3)=0.
	DMBO2(3,1)=0.
	DMBO2(3,2)=0.
	
CB      DEO(1,1)=SLA/PR
CB	DEO(1,2)=SLA
CB	DEO(2,1)=SLA
CB	DEO(2,2)=SLA/PR
CB	DEO(3,3)=SLA*(1-PR)/2
CB	DEO(1,3)=0.
CB	DEO(2,3)=0.
CB	DEO(3,1)=0.
CB	DEO(3,2)=0.	
      DEO2 = 0.0
	DO I=1,3
	DO J=1,3
      DEO2(I,J)=DMBO2(I,J)*THK
      DEO2(I+3,J+3)=THK*THK*DMBO2(I,J)/6.        
      !DEO2(I+3,J+3)=THK*THK*THK*DMBO2(I,J)/12.
      !DEO(I+3,J+3)=DEO(I,J)
	ENDDO
	ENDDO
	
	DEO2(7,7)=YM*THK
	DEO2(8,8)=SMU*THK*(5.0/6.0) !SHEAR COEFFICIENT (5/6)
	DEO2(9,9)=SMU*THK*(5.0/6.0) !SHEAR COEFFICIENT (5/6) 
      
		
	RETURN
	END
C***************************************************************************************
      SUBROUTINE SOKNEW(S,BB,DE,DVOL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ==========================================================
C     CALCULATE LINEAR STIFFNESS MATRIX FOR 8 NODE SOLID ELEMENT
C     ==========================================================
      DIMENSION S(300),STIFF(24,24),BB(9,24),DE(9,9)

	STIFF=MATMUL(MATMUL(TRANSPOSE(BB),DE),BB)*DVOL
	K=0
	DO I=1,24
	DO J=I,24
	K=K+1
	S(K)=STIFF(I,J)+S(K)
	ENDDO
	ENDDO  
	RETURN
	END
C*****************************************************************************************
      SUBROUTINE SOBATH(BA,VR,VS,VT,XJ,SH,II)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------------
C     GET THE NATURAL THICKNESS STRAIN-DISP MATRIX AT THE SAMPLING POINTS
C     -------------------------------------------------------------------
      COMMON /DISPL/ IDISP(4,6)
      DIMENSION BA(4,24),VR(3),VS(3),VT(3),XJ(3,3),SH(4)
	DIMENSION W1(3,3),WW(1,3),WORK(1,6)
C
	DO I=1,3
	W1(1,I)=VR(I)
	W1(2,I)=VS(I)
	W1(3,I)=VT(I)
	ENDDO
	DO I=1,4
      WW(1,1)=0.
	WW(1,2)=0.
	WW(1,3)=SH(I)
	WW=MATMUL(WW,XJ)
      WW=MATMUL(WW,W1)
      DO J=1,3
	WORK(1,2*J)=-WW(1,J)*.5
	WORK(1,2*J-1)=WW(1,J)*.5
      ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BA(II,N)=WORK(1,M)
      ENDDO
C    
	ENDDO
C
      RETURN
	END
C*********************************************************************************
      SUBROUTINE SOBT(RI,SI,RA,SA,XJI,BB,BAT,BBL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C     CALCULATE THE ASSUMED THICKNESS NORMAL STRAIN-DISPLACEMENT MATRIX 
C     THE CORRESPONDING DISPLACEMENTS TO ARE IN GLOBAL COORDINATE
C     -----------------------------------------------------------------
      DIMENSION RA(4),SA(4),XJI(3,3),BB(9,24),BAT(4,24),BBL(9,24)
	DIMENSION WORK(4,24)	
	DO IL=1,24
	DO I=1,4
	TH=.25*(1.+RA(I)*RI)*(1.+SA(I)*SI)
	WORK(I,IL)=TH*BAT(I,IL)
	ENDDO
	ENDDO
	DO IL=1,24
	DO I=1,4
	BB(7,IL)=BB(7,IL)+WORK(I,IL)*XJI(3,3)*XJI(3,3)
	BBL(7,IL)=BBL(7,IL)+WORK(I,IL) ! BJ LOCAL
	ENDDO
	ENDDO
C
      RETURN
	END
C******************************************************************************
C******************************************************************************
      SUBROUTINE SOKG(S,SH,SHR,SHS,SHT,VR,VS,VT,TAU,ZT)     
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------------------
C     CALCULATE GEOMETRIC STIFFNESS MATRIX AND ADD TO TANGENT STIFFNESS MATRIX
C     ------------------------------------------------------------------------
      COMMON /DISPL/ IDISP(4,6)
      DIMENSION S(300),SH(4),SHR(4),SHS(4),SHT(4)
	DIMENSION VR(3),VS(3),VT(3),TAU(9)
	DIMENSION STIFF(24,24)
	DIMENSION BRU(3,24),BSU(3,24),BRP(3,24),BSP(3,24),P(3,24)
	DIMENSION W(3),WORK(6),BTU(3,24)
C      
	DO II=1,3
      DO I=1,4
	IF(II.EQ.1)W=SHT(I)*VR
	IF(II.EQ.2)W=SHT(I)*VS
	IF(II.EQ.3)W=SHT(I)*VT	
      DO J=1,3
      WORK(2*J)=W(J)*.5
	WORK(2*J-1)=W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BTU(II,N)=WORK(M)
      ENDDO
      ENDDO
	ENDDO
C
      DO II=1,3
      DO I=1,4
	IF(II.EQ.1)W=SHR(I)*VR
	IF(II.EQ.2)W=SHR(I)*VS
	IF(II.EQ.3)W=SHR(I)*VT	
      DO J=1,3
      WORK(2*J)=W(J)*.5
	WORK(2*J-1)=W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BRU(II,N)=WORK(M)
      ENDDO
      ENDDO
	ENDDO
C
      DO II=1,3
      DO I=1,4
	IF(II.EQ.1)W=SHR(I)*VR
	IF(II.EQ.2)W=SHR(I)*VS
	IF(II.EQ.3)W=SHR(I)*VT	
      DO J=1,3
      WORK(2*J)=-W(J)*.5
	WORK(2*J-1)=W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BRP(II,N)=WORK(M)
      ENDDO
      ENDDO
	ENDDO
C
      DO II=1,3
      DO I=1,4
	IF(II.EQ.1)W=SHS(I)*VR
	IF(II.EQ.2)W=SHS(I)*VS
	IF(II.EQ.3)W=SHS(I)*VT	
      DO J=1,3
      WORK(2*J)=W(J)*.5
	WORK(2*J-1)=W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BSU(II,N)=WORK(M)
      ENDDO
      ENDDO
	ENDDO
C
      DO II=1,3
      DO I=1,4
	IF(II.EQ.1)W=SHS(I)*VR
	IF(II.EQ.2)W=SHS(I)*VS
	IF(II.EQ.3)W=SHS(I)*VT	
      DO J=1,3
      WORK(2*J)=-W(J)*.5
	WORK(2*J-1)=W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BSP(II,N)=WORK(M)
      ENDDO
      ENDDO
	ENDDO
C
      DO II=1,3
      DO I=1,4
	IF(II.EQ.1)W=SH(I)*VR
	IF(II.EQ.2)W=SH(I)*VS
	IF(II.EQ.3)W=SH(I)*VT	
      DO J=1,3
      WORK(2*J)=-W(J)*.5
	WORK(2*J-1)=W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	P(II,N)=WORK(M)
      ENDDO
      ENDDO
	ENDDO
C
C	-----------------------------------------------------------
C	MEMBRANE
	STIFF=MATMUL(TRANSPOSE(BRU),BRU)*TAU(1) 
	STIFF=STIFF+MATMUL(TRANSPOSE(BSU),BSU)*TAU(2) 
	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),BSU)*TAU(3)+
     1MATMUL(TRANSPOSE(BSU),BRU)*TAU(3)

C	BENDING
	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),BRP)*TAU(4)+
	1MATMUL(TRANSPOSE(BRP),BRU)*TAU(4)
	STIFF=STIFF+MATMUL(TRANSPOSE(BSU),BSP)*TAU(5)+
	1MATMUL(TRANSPOSE(BSP),BSU)*TAU(5)
	STIFF=STIFF+MATMUL(TRANSPOSE(BSU),BRP)*TAU(6)+
	1MATMUL(TRANSPOSE(BRP),BSU)*TAU(6)
	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),BSP)*TAU(6)+
	1MATMUL(TRANSPOSE(BSP),BRU)*TAU(6)

C	THICKNESS
	STIFF=STIFF+MATMUL(TRANSPOSE(BTU),BTU)*TAU(7)
CC	STIFF=STIFF+MATMUL(TRANSPOSE(P),P)*TAU(7)*.5*ZT*ZT

C	SHEAR
	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),BTU)*TAU(8)+
	1MATMUL(TRANSPOSE(BTU),BRU)*TAU(8)
CC	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),P)*TAU(8)*ZT
	STIFF=STIFF+MATMUL(TRANSPOSE(BSU),BTU)*TAU(9)+
	1MATMUL(TRANSPOSE(BTU),BSU)*TAU(9)
CC	STIFF=STIFF+MATMUL(TRANSPOSE(BSU),P)*TAU(9)*ZT
C	-----------------------------------------------------------


C	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),BSP)*TAU(6)
C	STIFF=STIFF+MATMUL(TRANSPOSE(BRP),BSU)*TAU(6)
C	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),BRP)*TAU(4)
C	STIFF=STIFF+MATMUL(TRANSPOSE(BSU),BSP)*TAU(5)
C	STIFF=STIFF+MATMUL(TRANSPOSE(BTU),BTU)*TAU(7)*.5
C	STIFF=STIFF+MATMUL(TRANSPOSE(P),P)*TAU(7)*.5*ZT
C	STIFF=STIFF+MATMUL(TRANSPOSE(BRU),BTU+ZT*P)*TAU(8)
C	STIFF=STIFF+MATMUL(TRANSPOSE(BSU),BTU+ZT*P)*TAU(9)

    	K=0
	DO I=1,24
	DO J=I,24
	K=K+1
	S(K)=STIFF(I,J)+S(K)
	ENDDO
	ENDDO  
C
	RETURN
	END
C********************************************************************
      SUBROUTINE BB1MAT(BB1,BB,SH,SHR,SHS,SHT,VR,VS,VT,XJI)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     GET STRAIN-DISP MATRIX FOR THICK AND SHEAR
C     ------------------------------------------
      COMMON /DISPL/ IDISP(4,6)
      DIMENSION BB1(3,24),SH(4),SHR(4),SHS(4),SHT(4)
	DIMENSION XJI(3,3),W(2,3),T(3,3),WORK(2,6)
	DIMENSION W1(1,3),BB(9,24),VR(3),VS(3),VT(3)
C 
      BB1=0.
      DO I=1,3
      T(1,I)=VR(I)
	T(2,I)=VS(I)
	T(3,I)=VT(I)
	ENDDO

C
      DO I=1,4
	W=0.
	W(1,3)=SHR(I)
	W(2,3)=SHS(I)
	W=MATMUL(W,T) 
	DO J=1,3
	DO K=1,2
	WORK(K,2*J)=W(K,J)*.5
	WORK(K,2*J-1)=W(K,J)*.5
      ENDDO	
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	DO K=1,2
	BB1(K+1,N)=WORK(K,M)
      ENDDO
	ENDDO
	ENDDO
C
      DO I=1,4
	W=0.
	W(1,1)=SH(I)
	W(2,2)=SH(I)
	W=MATMUL(W,T) 
	DO J=1,3
	DO K=1,2
	WORK(K,2*J)=-W(K,J)*.5
	WORK(K,2*J-1)=W(K,J)*.5
      ENDDO	
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	DO K=1,2
	BB1(K+1,N)=WORK(K,M)+BB1(K+1,N)
      ENDDO
	ENDDO
	ENDDO	
C
      DO I=1,4
	W1=0.
	W1(1,3)=SHT(I)
	W1=MATMUL(W1,T) 
	DO J=1,3
	WORK(1,2*J)=W1(1,J)*.5
	WORK(1,2*J-1)=W1(1,J)*.5	
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BB1(1,N)=WORK(1,M) 
	ENDDO
	ENDDO
C
      DO I=1,4
	W1=0.
	W1(1,3)=SH(I)*XJI(3,3)
	W1=MATMUL(W1,T) 
	DO J=1,3
	WORK(1,2*J)=-W1(1,J)*.5
	WORK(1,2*J-1)=W1(1,J)*.5	
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
	BB1(1,N)=WORK(1,M)+BB1(1,N)
	ENDDO
	ENDDO
C
      DO I=1,24
	DO J=1,3
      BB1(J,I)=BB(J+6,I)-BB1(J,I)
	ENDDO
	ENDDO
      RETURN
	END
C****************************************************************************
      SUBROUTINE WARP(VR,VS,VT,FRS,GRS,BB,BBL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------
C     MODIFY BB FOR WARPED ELEMENT SHAPE
C     ----------------------------------
      COMMON /DISPL/ IDISP(4,6)
      DIMENSION VR(3),VS(3),VT(3),FRS(4),GRS(4)
	DIMENSION BB(9,24)
	DIMENSION T(3,3),W(3,3),WORK(3,6)
	DIMENSION BBL(9,24),WORKL(3,6),WL(3,3)
C
      DO I=1,3
      T(1,I)=VR(I)
	T(2,I)=VS(I)
	T(3,I)=VT(I)
	ENDDO
C
      DO I=1,4
	W=0.
	W(1,1)=FRS(I)
	W(2,2)=GRS(I)
	W(3,1)=GRS(I)
	W(3,2)=FRS(I)
	
      WL = W ! BJ LOCAL ************
	W=MATMUL(W,T)
	
	DO J=1,3
	DO K=1,3
	WORK(K,2*J)   = W(K,J)*.5
	WORK(K,2*J-1) = W(K,J)*.5
	ENDDO
      ENDDO
      DO M=1,6
	N=IDISP(I,M)
	DO K=1,3
	BB(K+3,N)= WORK(K,M)+BB(K+3,N)
	ENDDO
      ENDDO  
	
	DO J=1,3 ! BJ LOCAL
	DO K=1,3
	WORKL(K,2*J)   = WL(K,J)*.5
	WORKL(K,2*J-1) = WL(K,J)*.5
	ENDDO
      ENDDO
      DO M=1,6
	N=IDISP(I,M)
	DO K=1,3
	BBL(K+3,N)= WORKL(K,M)+BBL(K+3,N)
	ENDDO
      ENDDO
	ENDDO  
	
C
      RETURN
	END


C	==================================================================
C	==================================================================
C	==================================================================
      SUBROUTINE BSMAT(VR,VS,VT,SH,SHR,SHS,SHT,XJI,S,DVOL,DE,RE,TAU)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------
C     CALCULATE GLOBAL STRAIN-DISPLACEMENT MATRIX BB
C     ----------------------------------------------
      COMMON /DISPL/ IDISP(4,6)

      DIMENSION VR(3),VS(3),VT(3),S(300),DE(9,9),TAU(9)
	DIMENSION SH(4),SHR(4),SHS(4),SHT(4)
	DIMENSION XJI(3,3),BBS(2,24),BBT(1,24),SKG(24,24)
C
      DIMENSION W(3),WORK(6),DS(2,2),DT(1,1)
	DIMENSION RE(24,1),EDIS(24)
	DIMENSION QS(2,1),QT(1,1)
	
	
	QS(1:2,1) = TAU(8:9)
	QT(1,1)   = TAU(7)
	
C     ------------------------------------------------------
C     ------------------------------------------------------
C     ----EXTRA SHEAR PART-------
	DO II = 1,2
      DO I  = 1,4
	IF(II.EQ.1) W = SH(I)*XJI(3,3)*VR
	IF(II.EQ.2) W = SH(I)*XJI(3,3)*VS
      DO J=1,3
	WORK(2*J-1) = W(J)*.5
	WORK(2*J-0) =-W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
      BBS(II,N)=WORK(M) 
	ENDDO
	ENDDO
      ENDDO
      
C     ----EXTRA THICKNESS PART-------
	DO II = 1,1
      DO I  = 1,4
	IF(II.EQ.1) W = SH(I)*XJI(3,3)*VT
      DO J=1,3
	WORK(2*J-1) = W(J)*.5
	WORK(2*J-0) =-W(J)*.5
	ENDDO
      DO M=1,6
	N=IDISP(I,M)
      BBT(II,N)=WORK(M) 
	ENDDO
	ENDDO
      ENDDO      
C     ------------------------------------------------------
C     ------------------------------------------------------

	RE = RE + MATMUL(TRANSPOSE(BBS),QS) + MATMUL(TRANSPOSE(BBT),QT)
	
      RETURN
C     ------------------------------------------------------  
C     ------------------------------------------------------  
      
      
      	
C
C     ------------------------------------------------------      	
      DS = 0.
	DT = 0.

	DT(1,1)=2./3.*DE(7,7)/2.0
	DS(1,1)=2./3.*DE(8,8)/2.0
	DS(2,2)=2./3.*DE(9,9)/2.0
C     ------------------------------------------------------

	SKG = MATMUL(TRANSPOSE(BBS),MATMUL(DS,BBS))*DVOL
	K=0
	DO I=1,24
	DO J=I,24
	K=K+1
	S(K)=SKG(I,J)+S(K)
	ENDDO
	ENDDO  
	
	SKG = MATMUL(TRANSPOSE(BBT),MATMUL(DT,BBT))*DVOL
	K=0
	DO I=1,24
	DO J=I,24
	K=K+1
	S(K)=SKG(I,J)+S(K)
	ENDDO
	ENDDO  



      RETURN
	END

C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE SOHLET(PROPM,WA,EPS,SIG,NLAYR,DE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	==================================================================
	DIMENSION SIG(9),EPS(9),PROPM(1)
	DIMENSION DE(9,9),WA(1)


	CALL SOHLER(PROPM,WA(1),EPS,SIG,NLAYR,WA(10),DE)


      RETURN
	END

C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE SOHLER(PROPM,EPSP,EPS,SIG,NLAYR,WA,DE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	==================================================================
	DIMENSION EPSP(9),EPS(9),PROPM(1)
	DIMENSION SIG(9),DEPS(9),SIGL(6),DEPL(6),DMATX(6,6)
	DIMENSION DE(9,9),WA(7,NLAYR)

	YOUNG = PROPM(1)
	POISN = PROPM(2)
	YIELD = PROPM(3)
	HARDS = PROPM(4)
	ALPHA = 0.0D0 !PROPM(6) ASSUMED VON-MISES
	

	SIG = 0.0
	DE  = 0.0

	DO I = 1,9
	DEPS(I) = EPS(I)-EPSP(I)
	ENDDO


	THLR = 2.0/NLAYR

	DO 900 ILAYR = 1,NLAYR
	XETA = -1.0 + 0.5*THLR + (ILAYR-1)*THLR

	DO I = 1,3
	DEPL(I) = DEPS(I) + XETA*DEPS(I+3)
	ENDDO
	DO I = 4,6
	DEPL(I) = DEPS(I+3)
	ENDDO

	CALL DRUKLR(YOUNG,POISN,YIELD,HARDS,ALPHA,
	1		    WA(1,ILAYR),DEPL,SIGL,DMATX)


	DO 150 I = 1,3
	J = I+3
	SIG(I) = SIG(I) + SIGL(I)*THLR
150	SIG(J) = SIG(J) + SIGL(I)*THLR*XETA
	SIG(8) = SIG(8) + SIGL(4)*THLR
	SIG(9) = SIG(9) + SIGL(5)*THLR
	SIG(7) = SIG(7) + SIGL(6)*THLR

	CALL INGLYR(DMATX,THLR,XETA,DE)


900	CONTINUE

	DO I = 1,9
	EPSP(I) = EPS(I)
	ENDDO

C
      RETURN
	END

C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE INGLYR(CP,DZETG,ZETAG,DP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	==================================================================
	DIMENSION CP(6,6),DP(9,9)
C	-----------------------------------------
C	INITIALIZE ELASTO-PLASTIC RIGIDITY MATRIX	
C	-----------------------------------------
	FA = DZETG
	FC = FA*ZETAG
	FB = FC*ZETAG
	FE = DZETG*DZETG*DZETG/12.0

C	---------------------------
C	MEMBRANE & BENDING RIGIDITY
C	---------------------------
	DO 300 I = 1,3
	DO 300 J = 1,3
	K = I+3
	L = J+3
	DP(I,J) = DP(I,J) + CP(I,J)*FA
300	DP(K,L) = DP(K,L) + CP(I,J)*(FB+FE)

C	--------------
C	SHEAR RIGIDITY
C	--------------
	DO 350 I = 4,5
	DO 350 J = 4,5
	K = I+4
	L = J+4
350	DP(K,L) = DP(K,L) + CP(I,J)*FA

C	------------------------------------
C	COUPLE MEMBRANE AND BENDING RIGIDITY
C	(UPPER PART OF TRIANGLE)
C	------------------------------------
	DO 380 I = 1,3
	DO 380 J = 1,3
	L = J+3
380	DP(I,L) = DP(I,L) + CP(I,J)*FC

C	-------------------------------------
C	FILL IN LOWER PART OF RIGIDITY MATRIX
C	-------------------------------------
	DO 400 I = 1,6
	DO 400 J = 1,6
400	DP(J,I) = DP(I,J)

C	------------------
C	THICKNESS RIGIDITY
C	------------------
	DP(7,7) = DP(7,7) + CP(6,6)*FA

C
      RETURN

	END

C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE DRUKLR(YOUNG,POISN,UNIAX,HARDP,ALPI1,
	1		          WA,DEPS,SIG,DEP)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
	DIMENSION SIGP(6),EPSP(6),DEP(6,6)
	DIMENSION SIG(6),EPS(6)
	DIMENSION AVECT(6),DVECT(6)
	DIMENSION DSIG(6),DEPS(6)
	DIMENSION WA(7),WDEPS(6),WDEP(6,6)
	DIMENSION LVECN(6)

	DO I = 1,6
	SIGP(I) = WA(I)
	ENDDO
	EPSTN   = WA(7)

C	=================================================================
C	FLIP THE ROW AND COLUMN OF STRESSES AND RIGIDITY MATRIX COMPONENT
C	     Sx,Sy,Sxy,Szz,Sxz,Syz --> Sx,Sy,Szz,Sxy,Sxz,Syz
C	=================================================================
	DO I = 1,6
	WDEPS(I) = 0.0
	WDEPS(I) = DEPS(I)
	ENDDO
	DEPS = 0.0
	DEPS(1) = WDEPS(1)
	DEPS(2) = WDEPS(2)
	DEPS(3) = WDEPS(4)
	DEPS(4) = WDEPS(3)
	DEPS(5) = WDEPS(5)
	DEPS(6) = WDEPS(6)
C	=================================================================
	

	
	CALL SHARD(EPSTN,UNIAX,HARDP,YIELD)
	
	CALL SMATD(YOUNG,POISN,DEPS,DSIG,DEP)


	AFIRT   = 0.0
	DO 20 I = 1,6
	SIG(I)  = SIGP(I) + DSIG(I)
20	AFIRT   = AFIRT + SIG(I)

	IF(AFIRT.EQ.0.0) GO TO 300 !FIRST ENTRY STEP

	ASCND   = 0.0
	DO I = 1,6
	ASCND   = ASCND + SIGP(I)
	ENDDO

	IF(ASCND.EQ.0.0) THEN !SECOND ENTRY STEP
	PREY = 0.0
	ELSE
	CALL SSULR(SIGP,PREY,AVECT,ALPI1)
	ENDIF

	CALL SSULR(SIG ,CURY,AVECT,ALPI1)


	IF(CURY.LE.YIELD) GO TO 300  !ELASTIC IN THIS STEP

	IF(CURY.LT.PREY)  GO TO 300  !UNLOADING IN THIS STEP


C	DO III = 1,9
C	WRITE(110,113) (DE(III,JJJ),JJJ=1,9)
C	ENDDO
C113	FORMAT(9E16.5)
	


C	======================================================================
C	START ELASTIC-PLASTIC BEHAVIOR
C	======================================================================

C	RECALL THE PREVIOUS STRESS IN THE CASE OF ELASTIC-PLASTIC BEHAVIOR
	DO 30 I = 1,6
30	SIG(I) = SIGP(I)

	IF(PREY.LT.YIELD) THEN     !ELASTIC IN PREVIOUS STEP
	CALL SFALR(SIG,DSIG,YIELD,RFAC,ALPI1)
	ELSEIF(PREY.GE.YIELD) THEN !ALREADY YIELD IN PREVIOUS STEP
	RFAC = 0.0
	ENDIF

C	WRITE(*,*) RFAC

	CALL STNUM(CURY,YIELD,NSTEP,ASTEP)

	DO 40 I = 1,6
40	SIG(I) = SIGP(I) + RFAC*DSIG(I)

	DO 50 I = 1,6
50	DSIG(I) = (1.0-RFAC)*DSIG(I)/ASTEP


	DO 200 ISTEP = 1,NSTEP
	
	CALL SSULR(SIG,CURY,AVECT,ALPI1)


	CALL SMATD(YOUNG,POISN,AVECT,DVECT,DEP)

	ABETA = 0.0
	ADUMM = 0.0
	DO I = 1,6
	ADUMM = ADUMM + AVECT(I)*DSIG(I)
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DLAMD = ADUMM*ABETA
	IF(DLAMD.LT.0.0) DLAMD = 0.0

	BDUMM = 0.0
	DO 100 I = 1,6
	BDUMM = BDUMM + AVECT(I)*SIG(I)
	SIG(I) = SIG(I) + DSIG(I)
	1					 - DLAMD*DVECT(I)

100	CONTINUE

	EPSTN  = EPSTN  + DLAMD*BDUMM/CURY

	CALL SHARD(EPSTN,UNIAX,HARDP,YIELD) !!!


200	CONTINUE

	CALL SSULR(SIG,CURY,AVECT,ALPI1)

	BRING = 1.0
	IF(CURY.GT.YIELD) BRING = YIELD/CURY
	DO 250 I = 1,6
250	SIG(I) = BRING*SIG(I)


	CALL SSULR(SIG,CURY,AVECT,ALPI1)


	CALL SMATD(YOUNG,POISN,AVECT,DVECT,DEP)

	ABETA = 0.0
	DO I = 1,6
	ABETA = ABETA + AVECT(I)*DVECT(I)
	ENDDO
	ABETA = 1.0 / (ABETA + HARDP)

	DO 280 I = 1,6
	DO 280 J = 1,6
280	DEP(I,J) = DEP(I,J) !- DVECT(I)*DVECT(J)*ABETA


C	======================================================================
C	END OF ELASTIC-PLASTIC BEHAVIOR
C	======================================================================	

300	CONTINUE

	SIG(3) = SIGP(3) + DEPS(3)*YOUNG
	DEP(3,3) = YOUNG
	DEP(3,1) = 0.0D0
	DEP(3,2) = 0.0D0
	DEP(1,3) = 0.0D0
	DEP(2,3) = 0.0D0


	DO 350 I = 1,6
350	WA(I) = SIG(I)
	WA(7) = EPSTN

C	=================================================================
C	FLIP THE ROW AND COLUMN OF STRESSES AND RIGIDITY MATRIX COMPONENT
C	      Sxx,Syy,Szz,Sxy,Sxz,Syz --> Sxx,Syy,Sxy,Sxz,Syz,Szz
C	=================================================================
	DO I = 1,6
	WDEPS(I) = 0.0
	WDEPS(I) = SIG(I)
	ENDDO
	SIG = 0.0
	SIG(1) = WDEPS(1)
	SIG(2) = WDEPS(2)
	SIG(3) = WDEPS(4)
	SIG(4) = WDEPS(5)
	SIG(5) = WDEPS(6)
	SIG(6) = WDEPS(3)

	DO I = 1,6
	DO J = 1,6
	WDEP(I,J) = 0.0
	WDEP(I,J) = DEP(I,J)
	ENDDO
	ENDDO
	LVECN(1) = 1
	LVECN(2) = 2
	LVECN(3) = 4
	LVECN(4) = 5
	LVECN(5) = 6
	LVECN(6) = 3
	DEP = 0.0
	DO I = 1,6
	LNI = LVECN(I)
	DO J = 1,6
	LNJ = LVECN(J)
	DEP(I,J) = WDEP(LNI,LNJ)
	ENDDO
	ENDDO

C	=================================================================
C	=================================================================


	RETURN

	END

C	======================================================================	
C	======================================================================	
C	======================================================================	
	SUBROUTINE SHARD(EPSTN,UNIAX,HARDP,YIELD)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================


	YIELD = UNIAX + HARDP*EPSTN

	YIELD = YIELD/SQRT(3.0)

C	WRITE(*,*) YIELD,UNIAX


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE SMATD(YOUNG,POISN,EP,SG,STIFG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SG(6),EP(6),STIFG(6,6)

C	=======================================================================

	DO 110 I = 1,6
	DO 110 J = 1,6
110	STIFG(I,J) = 0.0

	STIFG(1,1) = 1.0 - POISN
	STIFG(2,2) = 1.0 - POISN
	STIFG(3,3) = 1.0 - POISN
	STIFG(1,2) = POISN
	STIFG(1,3) = POISN
	STIFG(2,3) = POISN
	STIFG(2,1) = POISN
	STIFG(3,1) = POISN
	STIFG(3,2) = POISN
	STIFG(4,4) = (1.0 - 2.0*POISN)/2.0
	STIFG(5,5) = (1.0 - 2.0*POISN)/2.0
	STIFG(6,6) = (1.0 - 2.0*POISN)/2.0
		 
	DO 150 I = 1,6
	DO 150 J = 1,6
150	STIFG(I,J) = STIFG(I,J)*YOUNG/(1.0+POISN)/(1.0-2.0*POISN)	



	DO 170 I = 1,6
	DO 170 J = 1,6
170	STIFG(I,J) = 0.0

	SMU = YOUNG/(2.0*(1.0+POISN))
	SLA = YOUNG*POISN/(1.0-POISN*POISN)

	STIFG(1,1)=SLA+2.*SMU
	STIFG(1,2)=SLA
	STIFG(2,1)=SLA
	STIFG(2,2)=SLA+2.*SMU
	STIFG(3,3)=YOUNG
	STIFG(4,4)=SMU
	STIFG(5,5)=SMU
	STIFG(6,6)=SMU



	DO 200 I = 1,6
	SG(I) = 0.0
	DO 200 J = 1,6
200	SG(I) = SG(I) + STIFG(I,J)*EP(J)

	SG(3) = 0.0D0  !SONGSAK SEP2006 CHECK

	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================	
	SUBROUTINE SSULR(SG,FN,AVECT,ALPI1)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)


	DIMENSION SG(6),AVECT(6),PST(3)
	DIMENSION AV1(6,1),AV2(6,1),AV3(6,1)
	DIMENSION AM1(6,6),AM2(6,6),AM3(6,6)
	

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

	CALL INVSR(SG,PST,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)

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
C	DRUCKER-PRACKER
C	=========================
	FN  = (ALPI1*VARI1) + SQRT(VARJ2)
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



	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE INVSR(SIGMA,PSTRES,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	======================================================================
C	SUBROUTINE TO COMPUTE THE PRINCIPLE STRESS AND INVARIANT
C	======================================================================

	DIMENSION SIGMA(6),PSTRES(3),STMAT(3,3)

	PI = 3.141592654

	SX  = SIGMA(1)
	SY  = SIGMA(2)
	SZ  = SIGMA(3)
	TXY = SIGMA(4)
	TXZ = SIGMA(5)
	TYZ = SIGMA(6)

	STMAT(1,1) = SX
	STMAT(1,2) = TXY
	STMAT(1,3) = TXZ

	STMAT(2,1) = TXY
	STMAT(2,2) = SY
	STMAT(2,3) = TYZ

	STMAT(3,1) = TXZ
	STMAT(3,2) = TYZ
	STMAT(3,3) = SZ

	
	VARI1 = 0.0
	VARI2 = 0.0
	VARI3 = 0.0

	DO I = 1,3
	VARI1 = VARI1 + STMAT(I,I)
	DO J = 1,3
	VARI2 = VARI2 + (0.5*STMAT(I,J)*STMAT(I,J))
	DO K = 1,3
	VARI3 = VARI3 + (STMAT(I,J)*STMAT(J,K)*STMAT(K,I)/3.0)
	ENDDO
	ENDDO
	ENDDO

	STMAT(1,1) = STMAT(1,1) - (VARI1/3.0)
	STMAT(2,2) = STMAT(2,2) - (VARI1/3.0)
	STMAT(3,3) = STMAT(3,3) - (VARI1/3.0)

	VARJ1 = 0.0
	VARJ2 = 0.0
	VARJ3 = 0.0

	DO I = 1,3
	VARJ1 = VARJ1 + STMAT(I,I)
	DO J = 1,3
	VARJ2 = VARJ2 + (0.5*STMAT(I,J)*STMAT(I,J))
	DO K = 1,3
	VARJ3 = VARJ3 + (STMAT(I,J)*STMAT(J,K)*STMAT(K,I)/3.0)
	ENDDO
	ENDDO
	ENDDO

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
	
C	...COMPUTE THE PRICIPLE STRESSES

	PSTRES(1)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+2.*PI/3.)+VARI1/3. !MAX		 
	PSTRES(2)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+4.*PI/3.)+VARI1/3. !MIN	 
	PSTRES(3)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA)+VARI1/3.          !MID

C	==========================
C	S1 > S3 > S2 = MAX/MIN/MID
C	==========================

	
	RETURN

	END

C	======================================================================
C	======================================================================
C	======================================================================
	SUBROUTINE SFALR(SG,DSG,YIELD,RFAC,ALPI1)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION SG(6),DSG(6)
	DIMENSION SGR(6),AVECT(6)


	NITER = 20
	RFAC1 = 0.0
	RFAC2 = 1.0
	TOL = 0.00001
	RFACO = RFAC2
	DO ITER = 1,NITER
	
	RFAC = 0.5*(RFAC1 + RFAC2)
	DO I = 1,6
	SGR(I) = SG(I) + RFAC*DSG(I)
	ENDDO
	CALL SSULR(SGR,FR,AVECT,ALPI1)
	TEST = FR - YIELD
	IF(TEST.LE.0.0) RFAC1 = RFAC
	IF(TEST.GT.0.0) RFAC2 = RFAC
	TEST = (RFAC-RFACO)/RFACO
	IF(ABS(TEST).LE.TOL) GO TO 200
	RFACO = RFAC
	ENDDO

100	CONTINUE

200	CONTINUE


	RETURN

	END	

C	======================================================================	
C	======================================================================	
C	======================================================================
	SUBROUTINE SRCLER(EPS,SIG,WA,DE,MSET,MTMOD,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION STNCR(8),STSCR(8),DP(8,8)
	DIMENSION EPS(9),SIG(9),DE(9,9),WA(1)
	DIMENSION LDP(8)


	LDP(1) = 1
	LDP(2) = 2
	LDP(3) = 3
	LDP(4) = 4
	LDP(5) = 5
	LDP(6) = 6
	LDP(7) = 8
	LDP(8) = 9

C	========================================================
	DO I = 1,9
	DO J = 1,9
	DE(I,J) = 0.0D0
	ENDDO
	ENDDO
C	========================================================
	DO I = 1,6
	STNCR(I) = EPS(I)
	ENDDO
	STNCR(7) = EPS(8)
	STNCR(8) = EPS(9)
C	========================================================

	IF(MTMOD.EQ.5) CALL RCSOLP(WA,STNCR,STSCR,DP,MSET,THYOG,REFANG)
	IF(MTMOD.EQ.6) CALL RCSOLE(WA,STNCR,STSCR,DP,MSET,THYOG,REFANG)

C	========================================================
	DO I = 1,6
	SIG(I) = STSCR(I)
	ENDDO
	SIG(7) = THYOG*EPS(7)
	SIG(8) = STSCR(7)
	SIG(9) = STSCR(8)
C	========================================================
	DO I = 1,8
	II = LDP(I)
	DO J = 1,8
	JJ = LDP(J)
	DE(II,JJ) = DP(I,J)
	ENDDO
	ENDDO
	DE(7,7) = THYOG
C	========================================================

C	WRITE(*,*) THYOG,MSET
C	PAUSE


	RETURN
	END


C	==================================================================
C	========= EPF MODEL REINFORCED CONCRETE SHELL ELEMENT ============
C	==================================================================
	SUBROUTINE RCSOLE(WA,STNCR,STSCR,DP,MSET,THYOG,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	=============================================================
C	=============================================================
	DIMENSION WA(41,1),STSPS(5),STNPS(5),STNCR(8),STSCR(8),
	1		  DMATX(5,5),STNZ(5),DSTNZ(5),DSTSZ(5),SIGMA(5),
     2		  DIREC(2),EPSTN(5),GRTSN(2),DP(8,8),CP(5,5),
     3		  WOKPR(3),EQIST(2),EQISN(3),ELSNP(5),EQCRK(6)

	THYOG = 0.0D0

	LLAYR = MSET
	NLAYR = MATLR(MSET,MNLYR+1)  !MNLYR
	ZETAD  = -1.0

	DO 5 I = 1,8
	STSCR(I) = 0.0
	DO 5 J = 1,8
5	DP(I,J)  = 0.0
	
C	------------------------------------------------------------------
C	LOOP OVER CONCRETE & STEEL LAYER
C	------------------------------------------------------------------
      DO 200 ILAYR = 1,NLAYR
	LPROP = MATLR(LLAYR,ILAYR)
      IDTCS = PROPS(LPROP,10)

	DO 8 I = 1,5
	J = I+5
	STNPS(I) = WA(I,ILAYR)
8	STSPS(I) = WA(J,ILAYR)

	DIREC(1) = WA(11,ILAYR)
	DIREC(2) = WA(12,ILAYR)

	EPSTN(1) = WA(13,ILAYR)
	EPSTN(2) = WA(14,ILAYR)
	EPSTN(3) = WA(15,ILAYR)
	EPSTN(4) = WA(16,ILAYR)
	EPSTN(5) = WA(17,ILAYR)

	ELSNP(1) = WA(18,ILAYR)
	ELSNP(2) = WA(19,ILAYR)
	ELSNP(3) = WA(20,ILAYR)
	ELSNP(4) = WA(21,ILAYR)
	ELSNP(5) = WA(22,ILAYR)

	GRTSN(1) = WA(23,ILAYR)
	GRTSN(2) = WA(24,ILAYR)

	EQIST(1) = WA(25,ILAYR)
	EQIST(2) = WA(26,ILAYR)

	EQISN(1) = WA(27,ILAYR)
	EQISN(2) = WA(28,ILAYR)
	EQISN(3) = WA(29,ILAYR)

	MSTAT    = WA(30,ILAYR)

	DEO2     = WA(31,ILAYR)
	DEMX2    = WA(32,ILAYR)

	WOKPR(1) = WA(33,ILAYR)
	WOKPR(2) = WA(34,ILAYR)
	WOKPR(3) = WA(35,ILAYR)

	EQCRK(1) = WA(36,ILAYR)
	EQCRK(2) = WA(37,ILAYR)
	EQCRK(3) = WA(38,ILAYR)
	EQCRK(4) = WA(39,ILAYR)
	EQCRK(5) = WA(40,ILAYR)
	EQCRK(6) = WA(41,ILAYR)


	IF(IDTCS.EQ.0) THEN
	DZETA = PROPS(LPROP,3)
	ZETA  = ZETAD + DZETA/2.0
	ZETAG = ZETA
	DZETG = DZETA
	ELSE
	DZETA = PROPS(LPROP,3)
	ZETA  = PROPS(LPROP,6)
	IF(ZETAD.LE.0.0) ZETA = -ZETA
	ZETAG = ZETA
	DZETG = DZETA
	ENDIF

C	------------------------------------------------------------------
C	COMPUTE STRAINS AT THE MIDDLE OF LAYER
C	------------------------------------------------------------------
	DO 10 I = 1,3
	STNZ(I)  = STNCR(I) + ZETAG*STNCR(I+3)
10	DSTNZ(I) = STNZ(I)  - STNPS(I)
	STNZ(4)  = STNCR(7)
	STNZ(5)  = STNCR(8)
	DSTNZ(4) = STNCR(7) - STNPS(4)
	DSTNZ(5) = STNCR(8) - STNPS(5)
	
	IF(IDTCS.EQ.1) DIREC(1) = PROPS(LPROP,7)
	IF(IDTCS.EQ.1) GO TO 100

C	------------------------------------------------------------------
C	COMPUTE STRESSES AT THE MIDDLE OF CONCRETE LAYER
C	------------------------------------------------------------------

	TENST = PROPS(LPROP,5)
	UNIAX = PROPS(LPROP,6)
	POISN = PROPS(LPROP,2)
	YOUNG = PROPS(LPROP,1)

	TYOG = PROPS(LPROP,1) !THICKNESS STIFFNESS SEP2006

	NSTAT = MSTAT
	ANGLE = DIREC(1)
	STRA1 = WOKPR(1)
	STRA2 = WOKPR(2)
	IF(NSTAT.EQ.1) THEN
	CALL	MODUL1(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,YOUNG,  0.0,
     2			  POISN)
	ELSE
	CALL	MODUL1(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,	0.0,  0.0,
     2			  0.0)
	ENDIF

	DO 15 I = 1,5
	DSTSZ(I) = 0.0
	DO 15 J = 1,5
15	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)

	DO 20 I = 1,5
	SIGMA(I) = 0.0
20	SIGMA(I) = STSPS(I) + DSTSZ(I)

C	------------------------------------------------------------------
C	COMPUTE MAX & MIN STRESS IN THE PRINCIPAL DIRECTION
C	------------------------------------------------------------------
	CALL	PRISTE(SIGMA,SGMAX,SGMIN,1)

C	------------------------------------------------------------------
C	CHECK FOR CRAKING OF CONCRETE
C	------------------------------------------------------------------
	FKO   = WOKPR(3)
	RF    = FKO*FKO*FKO
	CRIT  = 0.0
	CKDUM = SGMAX*SGMIN
	IF(CKDUM.LT.0.0) THEN
	CRIT = SGMAX/RF/TENST
	ENDIF
	IF(SGMAX.GT.0.0.AND.SGMIN.GT.0.0) THEN
	CRIT = (SGMAX/RF/TENST) + 0.26*(SGMIN/SGMAX)
	ENDIF
	IF(CRIT.GT.1.0) GO TO 30
C	IF(SGMAX.GT.TENST) GO TO 30
	IF(NSTAT.EQ.1) GO TO 40

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
30	CALL	CCRACKE(STSPS,NSTAT,SIGMA,DIREC,
	1			   EPSTN,GRTSN,SGMAX,LPROP,
	2			   STNZ,WOKPR)

	IF(NSTAT.EQ.3) GO TO 95
C	------------------------------------------------------------------
C	COMPRESION BEHAVIOR OF CONCRETE
C	------------------------------------------------------------------
40	CALL	CYIELDE(STSPS,NSTAT,SIGMA,DIREC,
	1			   EPSTN,LPROP,EQIST,EQISN,
	2			   STNZ,WOKPR,DEMX2,ELSNP,
	3			   EQCRK,KUNLO)



	GO TO 95

C	------------------------------------------------------------------
C	COMPUTE STRESS AT THE MIDDLE OF STEEL LAYER
C	------------------------------------------------------------------
100	NSTAT = MSTAT
	ANGLE = PROPS(LPROP,7)+REFANG

	TYOG = PROPS(LPROP,1) !THICKNESS STIFFNESS SEP2006

	CALL	MODUL1(DMATX,LPROP,	4,ANGLE,
	1			  0.0, 0.0, 0.0, 0.0, 0.0)

	DO 60 I = 1,5
	DSTSZ(I) = 0.0
	DO 60 J = 1,5
60	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)

	CALL	STELRE(STSPS,NSTAT,EPSTN,DSTSZ,LPROP,REFANG)

95	CONTINUE

C	------------------------------------------------------------------
C	UPDATE THE LAYER STRAIN
C	------------------------------------------------------------------
	DO 96 I = 1,5 
96	STNPS(I) = STNZ(I)
	
	MSTAT = NSTAT

C	------------------------------------------------------------------
C	STORE THE CURRENT STRESSES IN WA
C	------------------------------------------------------------------
	DO 97 I = 1,5
	J = I+5
	WA(I,ILAYR)  = STNPS(I)
97	WA(J,ILAYR)  = STSPS(I)

	WA(11,ILAYR) = DIREC(1)
	WA(12,ILAYR) = DIREC(2)

	WA(13,ILAYR) = EPSTN(1)
	WA(14,ILAYR) = EPSTN(2)
	WA(15,ILAYR) = EPSTN(3)
	WA(16,ILAYR) = EPSTN(4)
	WA(17,ILAYR) = EPSTN(5)

	WA(18,ILAYR) = ELSNP(1)
	WA(19,ILAYR) = ELSNP(2)
	WA(20,ILAYR) = ELSNP(3)
	WA(21,ILAYR) = ELSNP(4)
	WA(22,ILAYR) = ELSNP(5)

	WA(23,ILAYR) = GRTSN(1)
	WA(24,ILAYR) = GRTSN(2)

	WA(25,ILAYR) = EQIST(1)
	WA(26,ILAYR) = EQIST(2)

	WA(27,ILAYR) = EQISN(1)
	WA(28,ILAYR) = EQISN(2)
	WA(29,ILAYR) = EQISN(3)

	WA(30,ILAYR) = MSTAT

	WA(31,ILAYR) = DEO2
	WA(32,ILAYR) = DEMX2

	WA(33,ILAYR) = WOKPR(1)
	WA(34,ILAYR) = WOKPR(2)
	WA(35,ILAYR) = WOKPR(3)

	WA(36,ILAYR) = EQCRK(1)
	WA(37,ILAYR) = EQCRK(2)
	WA(38,ILAYR) = EQCRK(3)
	WA(39,ILAYR) = EQCRK(4)
	WA(40,ILAYR) = EQCRK(5)
	WA(41,ILAYR) = EQCRK(6)

C	------------------------------------------------------------------
C	CONTRIBUTION OF STRESSES OVER THICKNESS
C	MEMBRANE - BENDING - SHEAR
C	------------------------------------------------------------------
	DO 150 I = 1,3
	J = I+3
	STSCR(I) = STSCR(I) + STSPS(I)*DZETG
150	STSCR(J) = STSCR(J) + STSPS(I)*DZETG*ZETAG
	STSCR(7) = STSCR(7) + STSPS(4)*DZETG
	STSCR(8) = STSCR(8) + STSPS(5)*DZETG

	CALL	CONSTVE(STSPS,GRTSN,NSTAT,EPSTN,
	1			   DIREC,LPROP,CP,KUNLO,WOKPR,REFANG)

C	------------------------------------------------------------------
C	INITIALIZE ELASTO-PLASTIC RIGIDITY MATRIX	
C	------------------------------------------------------------------
	FA = DZETG
	FC = FA*ZETAG
	FB = FC*ZETAG
	FE = DZETG*DZETG*DZETG/12.0

C	------------------------------------------------------------------
C	MEMBRANE & BENDING RIGIDITY
C	------------------------------------------------------------------
	DO 300 I = 1,3
	DO 300 J = 1,3
	K = I+3
	L = J+3
	DP(I,J) = DP(I,J) + CP(I,J)*FA
300	DP(K,L) = DP(K,L) + CP(I,J)*(FB+FE)

C	------------------------------------------------------------------
C	SHEAR RIGIDITY
C	------------------------------------------------------------------
	DO 350 I = 4,5
	DO 350 J = 4,5
	K = I+3
	L = J+3
350	DP(K,L) = DP(K,L) + CP(I,J)*FA

C	------------------------------------------------------------------
C	COUPLE MEMBRANE AND BENDING RIGIDITY
C	(UPPER PART OF TRIANGLE)
C	------------------------------------------------------------------
	DO 380 I = 1,3
	DO 380 J = 1,3
	L = J+3
380	DP(I,L) = DP(I,L) + CP(I,J)*FC

C	------------------------------------------------------------------
C	FILL IN LOWER PART OF RIGIDITY MATRIX
C	------------------------------------------------------------------
	DO 400 I = 1,6
	DO 400 J = 1,6
400	DP(J,I) = DP(I,J)

C	------------------------------------------------------------------
C	THICKNESS RIGIDITY
C	------------------------------------------------------------------
	THYOG = THYOG + FA/TYOG


	IF(IDTCS.EQ.0) THEN
	ZETAD  = ZETA + DZETA/2.0
	ENDIF

200	CONTINUE

      THYOG = 2.0D0/THYOG
      

	RETURN
	END

C	==================================================================
C	==================================================================
C	===== PLASTICITY MODEL OF REINFORCED CONCRETE SHELL ELEMENT ======
C	==================================================================
	SUBROUTINE RCSOLP(WA,STNCR,STSCR,DP,MSET,THYOG,REFANG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	==================================================================
	COMMON /CONC/ MATLR(50,50),PROPS(50,10),MOPTN,MLAYR,MMATS,
	1			  MNLYR,NPONT
C	==================================================================
	DIMENSION WA(18,1),STSPS(5),STNPS(5),STNCR(8),STSCR(8),
	1		  DMATX(5,5),STNZ(5),DSTNZ(5),GRTST(2),EPSTN(2),
	2		  DSTSZ(5),SIGMA(5),DIREC(2),EFFST(1),DP(8,8),
	3		  CP(5,5),HDMAT(100,2)


	THYOG = 0.0D0

	LLAYR = MSET
	NLAYR = MATLR(MSET,MNLYR+1) !MNLYR
	ZETAD  = -1.0

	DO 250 I = 1,8
	STSCR(I) = 0.0
	DO 250 J = 1,8
250	DP(I,J)  = 0.0

C	------------------------------------------------------------------
C	LOOP OVER CONCRETE & STEEL LAYER
C	------------------------------------------------------------------
      DO 200 ILAYR = 1,NLAYR
	LPROP = MATLR(LLAYR,ILAYR)
      IDTCS = PROPS(LPROP,10)

	DO 8 I = 1,5
	J = I+5
	STNPS(I) = WA(I,ILAYR)
8	STSPS(I) = WA(J,ILAYR)
	DIREC(1) = WA(11,ILAYR)
	DIREC(2) = WA(12,ILAYR)
	EFFST(1) = WA(13,ILAYR)
	EPSTN(1) = WA(14,ILAYR)
	EPSTN(2) = WA(15,ILAYR)
	GRTST(1) = WA(16,ILAYR)
	GRTST(2) = WA(17,ILAYR)
	MSTAT    = WA(18,ILAYR)
	
	IF(IDTCS.EQ.0) THEN
	DZETA = PROPS(LPROP,3)
	ZETA  = ZETAD + DZETA/2.0
	ZETAG = ZETA
	DZETG = DZETA
	ELSE
	DZETA = PROPS(LPROP,3)
	ZETA  = PROPS(LPROP,6)
	IF(ZETAD.LE.0.0) ZETA = -ZETA
	ZETAG = ZETA
	DZETG = DZETA
	ENDIF

C	------------------------------------------------------------------
C	COMPUTE STRAINS AT THE MIDDLE OF LAYER
C	------------------------------------------------------------------
	DO 10 I = 1,3
	STNZ(I)  = STNCR(I) + ZETAG*STNCR(I+3)
10	DSTNZ(I) = STNZ(I)  - STNPS(I)
	STNZ(4)  = STNCR(7)
	STNZ(5)  = STNCR(8)
	DSTNZ(4) = STNCR(7) - STNPS(4)
	DSTNZ(5) = STNCR(8) - STNPS(5)
	
	IF(IDTCS.EQ.1) DIREC(1) = PROPS(LPROP,7)
	IF(IDTCS.EQ.1) GO TO 100

      IF(MOPTN.EQ.2)	CALL HDATM(HDMAT,LPROP)

C	------------------------------------------------------------------
C	COMPUTE STRESSES AT THE MIDDLE OF CONCRETE LAYER
C	------------------------------------------------------------------
	TYOG = PROPS(LPROP,1) !THICKNESS STIFFNESS SEP2006

	NSTAT = MSTAT
	ANGLE = DIREC(1)
	STRA1 = EPSTN(1)
	STRA2 = EPSTN(2)
	CALL	MODUL(DMATX,LPROP,NSTAT,ANGLE,
	1			  STRA1,STRA2,	0.0,	0.0)
	DO 15 I = 1,5
	DSTSZ(I) = 0.0
	DO 15 J = 1,5
15	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)

	DO 20 I = 1,5
20	SIGMA(I) = STSPS(I) + DSTSZ(I)

C	------------------------------------------------------------------
C	COMPUTE MAX & MIN STRESS IN THE PRINCIPAL DIRECTION
C	------------------------------------------------------------------
	CALL	PRIST(SIGMA,SGMAX,SGMIN)

	TENST = PROPS(LPROP,5)
	UNIAX = PROPS(LPROP,6)

C	------------------------------------------------------------------
C	CHECK FOR CRAKING OF CONCRETE
C	------------------------------------------------------------------
	IF(SGMAX.GT.TENST) GO TO 30
	IF(NSTAT.EQ.1.OR.NSTAT.EQ.4) GO TO 40
	IF(NSTAT.EQ.7) GO TO 40

C	------------------------------------------------------------------
C	CRACKING OF CONCRETE
C	------------------------------------------------------------------
30	CALL	CCRACK(STSPS,STNPS,STNCR,EFFST,
	1			   SIGMA,NSTAT,DIREC,EPSTN,
	2			   GRTST,SGMAX,LPROP,YIELD,
	3			   STNZ,HDMAT)

C	------------------------------------------------------------------
C	CHECK FOR YIELDING OF CONCRETE	
C	------------------------------------------------------------------
	IF(MOPTN.EQ.2) THEN
	CALL	HARDN(EPSTN,LPROP,HARDS,YIELS,HDMAT)
	ELSE
	YIELS = UNIAX
	HARDS = 0.0
	ENDIF

	IF(YIELD.LE.YIELS) GO TO 50
	DO 35 I = 1,5
35	DSTSZ(I) = SIGMA(I) - STSPS(I)

C	------------------------------------------------------------------
C	COMPRESION BEHAVIOR OF CONCRETE
C	------------------------------------------------------------------
40	CALL	CYIELD(STSPS,STNPS,STNCR,EFFST,
	1			   NSTAT,DIREC,EPSTN,STNZ,
	2			   DSTSZ,SIGMA,LPROP,HDMAT)
	
50	CONTINUE

	GO TO 95

C	------------------------------------------------------------------
C	COMPUTE STRESS AT THE MIDDLE OF STEEL LAYER
C	------------------------------------------------------------------
100	NSTAT = MSTAT
	ANGLE = PROPS(LPROP,7)+REFANG

	TYOG = PROPS(LPROP,1) !THICKNESS STIFFNESS SEP2006

	CALL	MODUL(DMATX,LPROP,	5,ANGLE,
	1			  0.0,	0.0,	0.0,	0.0)
	DO 60 I = 1,5
	DSTSZ(I) = 0.0
	DO 60 J = 1,5
60	DSTSZ(I) = DSTSZ(I) + DMATX(I,J)*DSTNZ(J)
	CALL	STELR(STSPS,NSTAT,DIREC,EFFST,EPSTN,
	1			  DSTSZ,LPROP,REFANG)

95	CONTINUE

C	UPDATE THE MATERIAL NUMBER
	MSTAT = NSTAT

C	UPDATE THE LAYER STRAIN
	DO 96 I = 1,5 
96	STNPS(I) = STNZ(I)
	

C	STORE THE CURRENT STRESSES IN WA
	DO 97 I = 1,5
	J = I+5
	WA(I,ILAYR)  = STNPS(I)
97	WA(J,ILAYR)  = STSPS(I)
	WA(11,ILAYR) = DIREC(1)
	WA(12,ILAYR) = DIREC(2)
	WA(13,ILAYR) = EFFST(1)
	WA(14,ILAYR) = EPSTN(1)
	WA(15,ILAYR) = EPSTN(2)
	WA(16,ILAYR) = GRTST(1)
	WA(17,ILAYR) = GRTST(2)
	WA(18,ILAYR) = MSTAT

C	CONTRIBUTION OF STRESSES OVER THICKNESS
C	MEMBRANE - BENDING - SHEAR
	DO 150 I = 1,3
	J = I+3
	STSCR(I) = STSCR(I) + STSPS(I)*DZETG
150	STSCR(J) = STSCR(J) + STSPS(I)*DZETG*ZETAG
	STSCR(7) = STSCR(7) + STSPS(4)*DZETG
	STSCR(8) = STSCR(8) + STSPS(5)*DZETG

	CALL	CONSTV(STSPS,GRTST,NSTAT,EPSTN,DIREC,
	1			   LPROP,CP,HDMAT,REFANG)

C	INITIALIZE ELASTO-PLASTIC RIGIDITY MATRIX	
	FA = DZETG
	FC = FA*ZETAG
	FB = FC*ZETAG
	FE = DZETG*DZETG*DZETG/12.0

C	MEMBRANE & BENDING RIGIDITY
	DO 300 I = 1,3
	DO 300 J = 1,3
	K = I+3
	L = J+3
	DP(I,J) = DP(I,J) + CP(I,J)*FA
300	DP(K,L) = DP(K,L) + CP(I,J)*(FB+FE)

C	SHEAR RIGIDITY
	DO 350 I = 4,5
	DO 350 J = 4,5
	K = I+3
	L = J+3
350	DP(K,L) = DP(K,L) + CP(I,J)*FA

C	COUPLE MEMBRANE AND BENDING RIGIDITY
C	(UPPER PART OF TRIANGLE)
	DO 380 I = 1,3
	DO 380 J = 1,3
	L = J+3
380	DP(I,L) = DP(I,L) + CP(I,J)*FC

C	FILL IN LOWER PART OF RIGIDITY MATRIX
	DO 400 I = 1,6
	DO 400 J = 1,6
400	DP(J,I) = DP(I,J)

C	------------------------------------------------------------------
C	THICKNESS RIGIDITY
C	------------------------------------------------------------------
	THYOG = THYOG + FA/TYOG


	IF(IDTCS.EQ.0) THEN
	ZETAD  = ZETA + DZETA/2.0
	ENDIF

200	CONTINUE

      THYOG = 2.0D0/THYOG

	RETURN
	END


C	======================================================================	
C	======================================================================	
C	======================================================================
C     --------------------------------------------------------------
C     --------------------------------------------------------------
      SUBROUTINE SOLINEW_OLD(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,
	1				       RE,MWG,FIN,MSET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     MAIN PROGRAM FOR THE 3-D SOLID---preintegration
C     EVALUATES THE TANGENTIAL STIFFNESS MATRIX,STRAINS AND STRESSES
C     AND EQUILIBRIUM LOADS FOR THE CURVILINEAR ISOPARAMETRIC
C     HEXAHEDRON (8 TO 21 NODES)
C	--------------------------
C     INPUT VARIABLES
C	---------------
C     PROPM(NMP)    = MATERIAL PROPERTIES (YM,PR,YLD,HP,DEN)
C     PROPG(NGP)    = GEOMETRIC PROPERTIES (NNO)
C     NODEX(NEX)    = LOCATIONS OF EXCESS NODES (MIDSIDE NODES)
C     WA(MWG,NPT)   = WORKING ARRAY (6 STRESSES + (6 STRAINS,YLD,IPEL))
C     COORD(3,NNO)  = CURRENT NODAL COORDINATES X,Y,Z
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     EDISI(NEF)    = CURRENT NODAL DISPLACEMENT INCREMENTS
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     S(NWS)        = ELEMENT STIFFNESS MATRIX (UPPER TRIANG.ROW-WISE)
C     RE(NEF)       = EQUILIBRIUM LOADS AT ELEMENT NODES
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/
C	--------------------------------
C     NAME(2)       = NAME OF ELEMENT MODULE
C     ITYPE         = CODE NUMBER FOR ELEMENT MODULE
C     ISTYP         = ELEMENT SUBTYPE
C     NLOPT         = CODE FOR NONLINEAR OPTION
C     NLOPT=0         LINEAR ANALYSIS
C     NLOPT=1         MATERIALLY NONLINEAR ONLY
C     NLOPT=2,3       TOTAL LAGRANGIAN,UPDATED LAGRANGIAN
C     MTMOD         = CODE FOR MATERIAL MODULE
C     MTMOD=1         LINEAR ELASTIC,ISOTROPIC
C     MTMOD=2         LINEAR ELASTIC,ORTHOTROPIC
C     MTMOD=3         ELASTO-PLASTIC VON-MISES
C     MTMOD=5         CONCRETE WITH CRACKING
C     NSINC         = FACTOR CONTROLLING NUMBER OF SUBINCREMENTS
C     ITOLEY        = TOLERANCE ON YIELD FUNCTION
C     NELE          = NUMBER OF ELEMENTS IN THIS GROUP
C     NMPS          = NUMBER OF MATERIAL PROPERTY SETS
C     NGPS          = NUMBER OF GEOMETRIC PROPERTY SETS
C     NMP           = NUMBER OF MATERIAL PROPERTIES PER SET
C     NGP           = NUMBER OF GEOMETRIC PROPERTIES PER SET
C     NNM           = MAXIMUM NUMBER OF NODES FOR ANY ONE ELEMENT
C     NEX           = MAXIMUM NUMBER OF EXCESS NODES
C     NCO           = NUMBER OF NODAL COORDINATES
C     NNF           = NUMBER OF NODAL DEGREES OF FREEDOM
C     NEF           = MAXIMUM NUMBER OF ELEMENT DEGREES OF FREEDOM
C     NWG           = NUMBER OF STORAGE LOCATIONS AT EACH GAUSS POINT
C     NPT           = NUMBER OF GAUSS POINTS
C     NWA           = SIZE OF WORKING ARRAY
C     NWS           = SIZE OF ELEMENT STIFFNESS MATRIX
C     MEL           = CURRENT ELEMENT NUMBER
C     NNO           = NUMBER OF NODES FOR THIS ELEMENT
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR THIS ELEMENT
C     NELTOT        = TOTAL NUMBER OF ELEMENTS (ALL GROUPS)
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /GAUS/
C	--------------------------------
C     GLOC(4,4)     = NATURAL GAUSS POINT COORDINATES (1*1 TO 4*4)
C     GWT (4,4)     = GAUSS POINT WEIGHTS
C     NGR,NGS,NGT   = NUMBER OF GAUSS POINTS IN RN,SN,TN DIRECTION
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /FLAG/
C	--------------------------------
C     IFPRI,ISPRI   = FLAG FOR PRINTING DISPL.OR STRESSES (ISPRI=0)
C     IFPLO         = FLAG FOR PLOT OUTPUT
C     IFREF         = FLAG FOR REFORMATION OF STIFFNESS (IFREF=0)
C     IFEIG         = FLAG FOR EIGENVALUE SOLUTION (IFEIG=0)
C     ITASK = 1       FIRST ENTRY INTO ELEMENT MODULE
C     ITASK = 2       ENTRY DURING EQUILIBRIUM ITERATIONS
C     ITASK = 3       ENTRY TO WORK OUT STRESSES (LAST STEP ONLY)
C     ITASK = 4       ENTRY TO DETERMINE GEOMETRIC STIFF.MATRIX ONLY
C     KSTEP           CURRENT STEP NUMBER
C     KITE            CURRENT ITERATION NUMBER
C     ----------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

	COMMON /DISPL/ IDISP(4,6)

C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(300),COORD(3,8)
      DIMENSION EDIS(24),EDISI(24),RE(24),XJ(3,3),XJI(3,3),COORDI(3,8)
      DIMENSION STRAIN(9),QSTRAI(9),STRESS(9),TAU(9)
C
      DIMENSION VR(3),VS(3),VT(3),SH(4),P(2,4),SHR(4),SHS(4),SHT(4)
	DIMENSION XY(3,4),PXY(3,4),DE(9,9)
C
      DIMENSION BB(9,24),DISL(24),U(3,4),PU(3,4)
      DIMENSION XJA(3,3),XJB(3,3),XJC(3,3),XJD(3,3)
	DIMENSION PA(4),PB(4),PC(4),PD(4)
	DIMENSION VRA(3),VSA(3),VTA(3)
      DIMENSION VRB(3),VSB(3),VTB(3)
	DIMENSION VRC(3),VSC(3),VTC(3)
	DIMENSION VRD(3),VSD(3),VTD(3)
	DIMENSION SHA(4),SHB(4),SHC(4),SHD(4)
	DIMENSION BAT(4,24),RA(4),SA(4)
	DIMENSION DUD(9),DPUD(9),STRAI(9)
	DIMENSION TAUA(3),BB1(3,24)
	DIMENSION FRS(4),GRS(4)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
	DIMENSION LVNUM(8)
C     ------------------------------------------------------
C     GET INTEGRATED MATERIAL MATRIX FOR ISOTROPIC MATERIALS
C     ------------------------------------------------------
      CALL HOKLAW (PROPM,PROPG,1)
      CALL DEMAT(DE)


	K = 0
	DO J=1,8
	DO I=1,3
	K = K+1
	COORDI(I,J) = COORD(I,J)
	IF(NLOPT.EQ.3) COORDI(I,J) = COORD(I,J) - EDIS(K)
	ENDDO
	ENDDO

	

C	WRITE(*,*) EDIS

C     ------------------------------------------------------
C	DETECTED THE THIN DIRECTION
C     ------------------------------------------------------
	XYZ1 = SQRT( (COORDI(1,2)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,2)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,2)-COORDI(3,1))**2.0 )
	XYZ2 = SQRT( (COORDI(1,4)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,4)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,4)-COORDI(3,1))**2.0 )
	XYZ3 = SQRT( (COORDI(1,5)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,5)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,5)-COORDI(3,1))**2.0 )

	NMX = 0
	IF(XYZ1.LT.XYZ2.AND.XYZ1.LT.XYZ3) NMX = 1
	IF(XYZ2.LT.XYZ1.AND.XYZ2.LT.XYZ3) NMX = 2
	IF(XYZ3.LT.XYZ1.AND.XYZ3.LT.XYZ2) NMX = 3

	IF(NMX.EQ.0) THEN
	IF(XYZ1.EQ.XYZ2) NMX = 1
	IF(XYZ1.EQ.XYZ3) NMX = 1
	IF(XYZ2.EQ.XYZ3) NMX = 2
	ENDIF

	IF(NMX.EQ.1) THEN
	LVNUM(1) = 1
	LVNUM(2) = 4
	LVNUM(3) = 8
	LVNUM(4) = 5
	LVNUM(5) = 2
	LVNUM(6) = 3
	LVNUM(7) = 7
	LVNUM(8) = 6
	ELSEIF(NMX.EQ.2) THEN
	LVNUM(1) = 1
	LVNUM(2) = 5
	LVNUM(3) = 6
	LVNUM(4) = 2
	LVNUM(5) = 4
	LVNUM(6) = 8
	LVNUM(7) = 7
	LVNUM(8) = 3
	ELSEIF(NMX.EQ.3) THEN
	LVNUM(1) = 1
	LVNUM(2) = 2
	LVNUM(3) = 3
	LVNUM(4) = 4
	LVNUM(5) = 5
	LVNUM(6) = 6
	LVNUM(7) = 7
	LVNUM(8) = 8
	ENDIF


	CALL SLFACE(NMX)

	K = 0
	DO I = 1,24
	RE(I) = 0.0D0
	DO J = I,24
	K = K + 1
	S(K) = 0.0D0
	ENDDO
	ENDDO
	
C     ------------------------------------------------------

C     --------------------------------------------
C     MIDSURFACE COORDINATES AND DIRECTION VECTORS
C     --------------------------------------------
      DO I=1,3
	DO J=1,4
	JJ  = LVNUM(J)
	JJJ = LVNUM(J+4)
      XY(I,J) =(COORD(I,JJ)+COORD(I,JJJ))/2.
      PXY(I,J)=(COORD(I,JJ)-COORD(I,JJJ))/2.
	
C	XY(I,J) =(COORD(I,J)+COORD(I,J+4))/2.
C	PXY(I,J)=(COORD(I,J)-COORD(I,J+4))/2.
	ENDDO
	ENDDO
C     --------------------------------------------
C     CALCULATE JACOBEAN MATRIX AT SAMPLING POINTS
C     --------------------------------------------
C     ----FOR SHEAR PART----
      RI=-1.0
	SI= 0.0
      CALL SHAP8NEW(RI,SI,SHA,P)
	DO I=1,4
	PA(I)=P(2,I)
	ENDDO
      CALL COVAB(XY,PXY,SHA,P,VRA,VSA,VTA,XJA,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	RI= 0.0
	SI=-1.0
	CALL SHAP8NEW(RI,SI,SHB,P)
	DO I=1,4
	PB(I)=P(1,I)
	ENDDO
      CALL COVAB(XY,PXY,SHB,P,VRB,VSB,VTB,XJB,XJI,
     @    	   SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	RI= 1.0
	SI= 0.0
	CALL SHAP8NEW(RI,SI,SHC,P)
	DO I=1,4
	PC(I)=P(2,I)
	ENDDO
      CALL COVAB(XY,PXY,SHC,P,VRC,VSC,VTC,XJC,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	RI= 0.0
	SI= 1.0
	CALL SHAP8NEW(RI,SI,SHD,P)
	DO I=1,4
	PD(I)=P(1,I)
	ENDDO
      CALL COVAB(XY,PXY,SHD,P,VRD,VSD,VTD,XJD,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
C     -----FOR THICKNESS PART-------
      RA(1) = 1.0D0
	SA(1) = 1.0D0
	RA(2) =-1.0D0
	SA(2) = 1.0D0
	RA(3) =-1.0D0
	SA(3) =-1.0D0
	RA(4) = 1.0D0
	SA(4) =-1.0D0
	DO I=1,4
	CALL SHAP8NEW(RA(I),SA(I),SH,P)
	CALL COVAB(XY,PXY,SH,P,VR,VS,VT,XJ,XJI,
     @	       SHR,SHS,SHT,FRS,GRS,DET,MEL,1)
	CALL SOBATH(BAT,VR,VS,VT,XJ,SH,I)
	ENDDO
C     ----------------------   
C     LOOP OVER GAUSS POINTS
C     ----------------------
      MGR=2
	MGS=2
      IPT = 0
      DO 900  IGR=1,MGR
      RI = GLOC(IGR,MGR)
      DO 900  IGS=1,MGS
      SI = GLOC(IGS,MGS)
      WT = GWT(IGR,MGR)*GWT(IGS,MGS)
      IPT = IPT+1
C     ------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P)
C     ------------------------------------
      CALL SHAP8NEW(RI,SI,SH,P)	
C     ------------------------------------------------------------------
C     SETUP THE COROTATIONAL COORDINATE (VR,VS,VT)
C	JACOBIAN (XJ), INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ------------------------------------------------------------------	    
 	CALL COVAB(XY,PXY,SH,P,VR,VS,VT,XJ,XJI,
     @	SHR,SHS,SHT,FRS,GRS,DET,MEL,0)	     
      DVOL = WT*DET	

C     --------------------------------------------------------------
C     ASSUMED GLOBAL SHEAR STRAIN-DISPLACEMENT MATRIX (STORED IN BB)
C     --------------------------------------------------------------   
      BB=0.0D0    
      CALL SOBQ(RI,SI,XJI,XJA,XJB,XJC,XJD,PA,PB,PC,PD,BB,VRA,VSA,VTA,
     +		  VRB,VSB,VTB,VRC,VSC,VTC,VRD,VSD,VTD,SHA,SHB,SHC,SHD) 
	CALL SOBT(RI,SI,RA,SA,XJI,BB,BAT)
C     ----------------------------------------
C     GET GLOBAL STRAIN-DISPLACEMENT MATRIX BB
C     ----------------------------------------
      CALL BBMAT(VR,VS,VT,SH,SHR,SHS,SHT,XJI,BB)
	CALL WARP(VR,VS,VT,FRS,GRS,BB)
C     ---------------------------
C     GET LINEAR STIFFNESS MATRIX
C     ---------------------------  
      IF(MTMOD.EQ.1) CALL SOKNEW(S,BB,DE,DVOL)	
C     -------------------------------------
C     NEXT PART FOR NONLINEAR ANALYSIS ONLY
C     -------------------------------------
      IF(NLOPT.LE.0) GOTO 900
C     ---------------------
C     LINEAR ASSUMED STRAIN
C     ---------------------
      STRAIN=MATMUL(BB,EDIS)
C     ------------------
C     LOCAL DISPLACEMENT
C     ------------------
      DO I=1,8
      K=3*I-2
	DISL(K)  =VR(1)*EDIS(K)+VR(2)*EDIS(K+1)+VR(3)*EDIS(K+2)
	DISL(K+1)=VS(1)*EDIS(K)+VS(2)*EDIS(K+1)+VS(3)*EDIS(K+2)
	DISL(K+2)=VT(1)*EDIS(K)+VT(2)*EDIS(K+1)+VT(3)*EDIS(K+2)
	ENDDO
      DO I=1,3
	DO J=1,4
C	K=3*(J-1)+I
C	L=3*(J+3)+I
	K = IDISP(J,2*I-1)
	L = IDISP(J,2*I-0)

      U(I,J) =(DISL(K)+DISL(L))/2.
      PU(I,J)=(DISL(K)-DISL(L))/2.
	ENDDO
	ENDDO
C     ------------------------------
C     LOCAL DISPLACEMENT DERIVATIVES
C     ------------------------------
      DUD=0.
	DPUD=0.
	DO J=1,3
	DO I=1,4
      DUD(J)   = SHR(I)*U(J,I) + DUD(J)
	DUD(J+3) = SHS(I)*U(J,I) + DUD(J+3)     
	DUD(J+6) = SHT(I)*U(J,I) + DUD(J+6)     
	DPUD(J)  = SHR(I)*PU(J,I)+ DPUD(J)
	DPUD(J+3)= SHS(I)*PU(J,I)+ DPUD(J+3)
	DPUD(J+6)=  SH(I)*PU(J,I)+ DPUD(J+6)
      ENDDO
	ENDDO
C     ----------------------
C     QUADRATIC STRAIN TERMS
C     ----------------------
C     ----MEMBRANE PART-----
      QSTRAI(1)=.5*(DUD(1)*DUD(1)+DUD(2)*DUD(2)+DUD(3)*DUD(3))
	QSTRAI(2)=.5*(DUD(4)*DUD(4)+DUD(5)*DUD(5)+DUD(6)*DUD(6))
	QSTRAI(3)=    DUD(1)*DUD(4)+DUD(2)*DUD(5)+DUD(3)*DUD(6)
C     ----BENDING PART------
      QSTRAI(4)= DUD(1)*DPUD(1)+DUD(2)*DPUD(2)+DUD(3)*DPUD(3)
      QSTRAI(5)= DUD(4)*DPUD(4)+DUD(5)*DPUD(5)+DUD(6)*DPUD(6)
	QSTRAI(6)= DPUD(1)*DUD(4)+DUD(1)*DPUD(4)+DPUD(2)*DUD(5)
     @          +DUD(2)*DPUD(5)+DPUD(3)*DUD(6)+DUD(3)*DPUD(6)
C     ----THICKNESS PART----
      ZT=XJI(3,3)
      QSTRAI(7)=.5*(DUD(7)*DUD(7)+DPUD(7)*DPUD(7)*ZT*ZT+DUD(8)*DUD(8)
     @     +DPUD(8)*DPUD(8)*ZT*ZT+DUD(9)*DUD(9)+DPUD(9)*DPUD(9)*ZT*ZT)
C     ----SHEAR PART--------
      QSTRAI(8)= DUD(1)*(DUD(7)+DPUD(7)*ZT)+DUD(2)*(DUD(8)+DPUD(8)*ZT)
     @           +DUD(3)*(DUD(9)+DPUD(9)*ZT)
	QSTRAI(9)= DUD(4)*(DUD(7)+DPUD(7)*ZT)+DUD(5)*(DUD(8)+DPUD(8)*ZT)
     @           +DUD(6)*(DUD(9)+DPUD(9)*ZT)
C     ---------------
C     ALMANSI STRAINS
C     ---------------
!	STRAI(1)=DUD(1)
!	STRAI(2)=DUD(5)
!	STRAI(3)=DUD(4)+DUD(2)
!	STRAI(4)=DPUD(1)
!	STRAI(5)=DPUD(4)
!	STRAI(6)=DPUD(4)+DPUD(2)
!	STRAI(7)=DPUD(9)
!	STRAI(8)=DPUD(7)+DUD(3)
!	STRAI(9)=DPUD(8)+DUD(6)
C
	IF(NLOPT.EQ.3) THEN
	DO II = 1,9
	STRAIN(II)=STRAIN(II) - QSTRAI(II)
	ENDDO
	ENDIF
C     -----------------------
C     STRESS IN DEFORMED BODY
C     -----------------------
      IF(MTMOD.EQ.1) STRESS=MATMUL(DE,STRAIN)

	
C	ELASTO PLASTIC DRACKER
	IF(MTMOD.EQ.4) THEN
	NLAYR = NPT/4
	CALL SOHLET(PROPM,WA(1,IPT),STRAIN,STRESS,NLAYR,DE)
	CALL SOKNEW(S,BB,DE,DVOL)
	ENDIF
	
C	RC EPF MODEL
	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6) THEN
	CALL SRCLER(STRAIN,STRESS,WA(1,IPT),DE,MSET,MTMOD)
	CALL SOKNEW(S,BB,DE,DVOL)
	ENDIF	

C     ---------------------
C     INTERNAL FORCE VECTOR
C     ---------------------
      TAU=STRESS*DVOL


	RE = RE + MATMUL(TRANSPOSE(BB),TAU)	


C     ----------------------------------------------------------------
C     TAKE INTO ACCOUNT EXTRA STRAIN TERMS IN SHEAR AND THICKNESS
C     ----------------------------------------------------------------
C	CALL BSMAT(VR,VS,VT,SH,SHR,SHS,SHT,XJI,S,DVOL,DE,RE,TAU)


      IF(NLOPT.EQ.3) CALL SOKG(S,SH,SHR,SHS,SHT,VR,VS,VT,TAU,ZT)
900   CONTINUE  
C

C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF



      RETURN
      END
C	==================================================================
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE TNEAS8(FJO,TTI)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION FJO(3,3),TTO(9,9),TTI(9,9)

C	FROM Z-SOIL BOOK
	TTO = 0.0

	DO I = 1,9
	TTO(I,I) = 1.0
	ENDDO

	TTO(1,1) = FJO(1,1)*FJO(1,1)
	TTO(1,2) = FJO(2,1)*FJO(2,1)
	TTO(1,3) = FJO(1,1)*FJO(2,1)

	TTO(2,1) = FJO(1,2)*FJO(1,2)
	TTO(2,2) = FJO(2,2)*FJO(2,2)
	TTO(2,3) = FJO(1,2)*FJO(2,2)

	TTO(3,1) = 2.0*FJO(1,1)*FJO(1,2)
	TTO(3,2) = 2.0*FJO(2,1)*FJO(2,2)
	TTO(3,3) = FJO(1,1)*FJO(2,2) + FJO(1,2)*FJO(2,1)


C	FROM CRISFIELD BOOK
	TTO = 0.0

	DO I = 1,9
	TTO(I,I) = 1.0
	ENDDO


	TTO(1,1) = FJO(1,1)*FJO(1,1)
	TTO(1,2) = FJO(1,2)*FJO(1,2)
	TTO(1,3) = FJO(1,1)*FJO(1,2)

	TTO(2,1) = FJO(2,1)*FJO(2,1)
	TTO(2,2) = FJO(2,2)*FJO(2,2)
	TTO(2,3) = FJO(2,2)*FJO(2,1)

	TTO(3,1) = 2.0*FJO(1,1)*FJO(2,1)
	TTO(3,2) = 2.0*FJO(1,2)*FJO(2,2)
	TTO(3,3) = FJO(1,1)*FJO(2,2) + FJO(1,2)*FJO(2,1)


	CALL INVMATRIX(TTO,TTI,9)


	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE MTEAS8(RN,SN,TTI,DETO,DET,GE,MM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION TTO(9,9),GE(9,MM),TTI(9,9),GEO(9,MM)

	CALL EASGE8(RN,SN,GEO,MM)

	DO I = 1,3
	DO J = 1,3
	TTO(I,J) = DETO*TTI(I,J)/DET
	ENDDO
	ENDDO

	GE = MATMUL(TTO,GEO)

	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE EASGE8(RN,SN,GE,MM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION GE(9,MM)

	GE = 0.0

	IF(MM.EQ.5) THEN
	GE(1,1) = RN
	GE(2,2) = SN
	GE(3,3) = RN
	GE(3,4) = SN
	GE(1,5) =  RN*SN
	GE(2,5) = -RN*SN
	GE(3,5) = RN*RN - SN*SN
	ENDIF

	IF(MM.EQ.7) THEN
	GE(1,1) = RN
	GE(2,2) = SN
	GE(3,3) = RN
	GE(3,4) = SN
	GE(1,5) = RN*SN
	GE(2,6) = RN*SN
	GE(3,7) = RN*SN
	ENDIF


	RETURN

	END

C	=====================================================================
C	=====================================================================
      SUBROUTINE SOMDSH(COORD,COORDI,XY,XYI,REDIS,H8,P,VR,VS,VT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------
C     MODIFIES TOTAL DISPLACEMENT VECTOR BY DEDUCTING
C     RIGID BODY TRANSLATIONS AND ROTATIONS
C	-------------------------------------
C     COORD(3,NNO)      = CURRENT NODAL COORDINATES
C     COORDI(3,NNO)     = INITIAL NODAL COORDINATES
C     EDIS(NEF)         = CURRENT NODAL DISPLACEMENTS
C     REDIS(48)         = COROTATIONAL FORM OF EDIS
C     H(8)              = SHAPE FUNCTIONS
C     HD(2,8)           = SHAPE FUNCTION DERIVATIVES
C     VR(3),VS(3),VT(3) = CURRENT DIRECTION COSINE VECTORS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     ----------------------------------------------------
      DIMENSION COORD(3,1),COORDI(3,1),XY(3,1),XYI(3,1),REDIS(1)
      DIMENSION XYZ(3),XYZO(3),CD(3,8),CDO(3,8)
      DIMENSION H8(1),P(2,1),VR(3),VS(3),VT(3)
      DIMENSION VRO(3),VSO(3),VTO(3),TM(3,3)
C
C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
      CALL CLEARA (XYZO,3)
	CALL CLEARA (XYZ ,3)
	
      DO 140 I=1,8
      DO 140 J=1,3
	XYZO(J)=XYZO(J)+H8(I)*COORDI(J,I)
  140 XYZ (J)=XYZ (J)+H8(I)*COORD (J,I)
  
      DO 150 I=1,8
      DO 150 J=1,3
	CDO(J,I)=COORDI(J,I)-XYZO(J)
  150 CD (J,I)=COORD (J,I)-XYZ (J)
  
  	CALL SOSHCOV(XYI,P,VRO,VSO,VTO)
      
      DO 160 I=1,3
      DO 160 J=1,3
 160  TM(I,J)=VR(I)*VRO(J)+VS(I)*VSO(J)+VT(I)*VTO(J)
 
      K=1
      DO 210 I=1,8
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)
      REDIS(K)=CD(J,I)-TCDO
  210 K=K+1
  
      RETURN
      END
C
C	=====================================================================
C	=====================================================================
      SUBROUTINE SOSHCOV(XY,P,VR,VS,VT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------
C     GET BASE VECTORS, JACOBEAN MATRIX AND ITS INVERSE
C     -------------------------------------------------
      DIMENSION VR(3),VS(3),VT(3),COVR(3),COVS(3)
	DIMENSION P(2,1)
	DIMENSION XY(3,1)
	
C	DIMENSION RV(3),SV(3)

C     -----------------------------------
C     TANGENTIAL DIRECTIONS TO MIDSURFACE
C     -----------------------------------
      COVR=0.
	COVS=0.
      DO I=1,3
	DO J=1,4
	COVR(I)=COVR(I)+XY(I,J)*P(1,J)
	COVS(I)=COVS(I)+XY(I,J)*P(2,J)
	ENDDO
	ENDDO	
C     -----------------------------
C     LOCAL COORDINATE BASE VECTORS
C     -----------------------------
      CALL VECPRD (COVR,COVS,VT)
	CALL SCALEN (VT,VT,DD,3)
	CALL VECPRD (VT,COVR,VS)
	CALL SCALEN (VS,VS,DD,3)
	CALL SCALEN (COVR,VR,DD,3)	

C      CALL VECPRD (COVR,COVS,VT)
C      CALL SCALEN (VT,VT,DD,3)
C      CALL SCALEN (COVR,RV,RL,3)
C      CALL SCALEN (COVS,SV,SL,3)
C      CALL VECPRD (VT,RV,VS)
C      CALL ADDVEC (SV,VS,VS)
C      CALL SCALEN (VS,VS,DD,3)
C      CALL VECPRD (VS,VT,VR)
C      CALL SCALEN (VR,VR,DD,3)
C      CALL SCALEN (VS,VS,DD,3)
C      CALL SCALEN (VT,VT,DD,3)
      
	RETURN
	END
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHSRLAX_SHL (WA,NWG,SIGR,SIGT,NAME)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*6 NAME
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      DIMENSION COORD(3,1),PROPG(1),VRF(3),VSF(3),VTF(3)
      DIMENSION VRO(3),VSO(3),VTO(3),SIGR(1),SIGT(1),SIG(6)

	CALL SHSTORL2_SHL (WA,NWG,VRO,VSO,VTO,VRF,VSF,VTF)

	IOPT = -1
	IF(NAME.EQ.'RESULT') IOPT = 0
	IF(NAME.EQ.'FIBER6') IOPT = 1
	IF(NAME.EQ.'FIBER3') IOPT = 2
	IF(IOPT.EQ.-1) RETURN

	SELECTCASE(IOPT)
	CASE(0)
	CALL STNRLAX_SHL (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	CASE(1)
	CALL STNFLAX_SHL (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)    
	ENDSELECT


      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE STNRLAX_SHL (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      DIMENSION VRF(3),VSF(3),VTF(3)
      DIMENSION VRO(3),VSO(3),VTO(3)
	DIMENSION SIGR(1),SIGT(1),SIG(9)


C	RESULTANT MEMBRANE
	SIG(1) = SIGR(1)
	SIG(2) = SIGR(2)
	SIG(3) = 0.0D0
	SIG(4) = SIGR(3)
	SIG(5) = SIGR(8)
	SIG(6) = SIGR(9)
      
      

	CALL STSTRAN(SIG,VRO,VSO,VTO,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	CALL STSTRAN(SIG,VRF,VSF,VTF,1)  !TRANSFORM GLOBAL STRESS TO LOCAL STRESS                  

	SIGT(1) = SIG(1)
	SIGT(2) = SIG(2)
	SIGT(6) = SIG(5)
	SIGT(8) = SIG(6)
	SIGT(9) = SIG(4)

C	RESULTANT MOMENT
	SIG(1) = SIGR(4)
	SIG(2) = SIGR(5)
	SIG(3) = 0.0D0
	SIG(4) = SIGR(6)
	SIG(5) = 0.0D0
	SIG(6) = 0.0D0
	CALL STSTRAN(SIG,VRO,VSO,VTO,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	CALL STSTRAN(SIG,VRF,VSF,VTF,1)  !TRANSFORM GLOBAL STRESS TO LOCAL STRESS
      
	
	SIGT(3) = SIG(1)
	SIGT(4) = SIG(2)
	SIGT(5) = SIG(4)



      RETURN
      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE STNFLAX_SHL (VRO,VSO,VTO,VRF,VSF,VTF,SIGR,SIGT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     --------------------------------------------------------------

      DIMENSION VRF(3),VSF(3),VTF(3)
      DIMENSION VRO(3),VSO(3),VTO(3)
	DIMENSION SIGR(1),SIGT(1),SIG(6)


	SIG(1:6) = SIGR(1:6)

	CALL STSTRAN(SIG,VRO,VSO,VTO,0)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS
	CALL STSTRAN(SIG,VRF,VSF,VTF,1)  !TRANSFORM LOCAL STRESS TO GLOBAL STRESS

	SIGT(1) = SIG(1)
	SIGT(2) = SIG(2)
	SIGT(3) = SIG(3)
	SIGT(4) = SIG(4)
	SIGT(5) = SIG(5)
	SIGT(6) = SIG(6)


      RETURN
      END



C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHSTORL1_SHL (COORD,PROPG,WA,NWG,VR,VS,VT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION COORD(3,1),PROPG(1)
	DIMENSION WA(1),VR(3),VS(3),VT(3)
	DIMENSION VRN(3),VSN(3),VTN(3)


	THETA = PROPG(11)
	LXOPT = INT(PROPG(12))
	CALL SHNEWBV_SHL(COORD,VRN,VSN,VTN,THETA,LXOPT)


	N1 = NWG-17
	N2 = N1+2
	WA(N1:N2) = VR(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VS(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VT(1:3)

	N1 = NWG-8
	N2 = N1+2
	WA(N1:N2) = VRN(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VSN(1:3)
	N1 = N2+1
	N2 = N1+2
	WA(N1:N2) = VTN(1:3)


	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHSTORL2_SHL (WA,NWG,VR,VS,VT,VRN,VSN,VTN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION WA(1),VR(3),VS(3),VT(3)
	DIMENSION VRN(3),VSN(3),VTN(3)


	N1 = NWG-17
	N2 = N1+2
	VR(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VS(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VT(1:3) = WA(N1:N2)

	N1 = NWG-8
	N2 = N1+2
	VRN(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VSN(1:3) = WA(N1:N2)
	N1 = N2+1
	N2 = N1+2
	VTN(1:3) = WA(N1:N2)



	RETURN
	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHNEWBV_SHL (COORD,VR,VS,VT,THETA,LXOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
C
C     --------------------------------------------------------------
C     CALLING BASE VECTOR FROM PROPG
C     --------------------------------------------------------------
C
      DIMENSION COORD(3,1),VR(3),VS(3),VT(3),V1(3),V2(3),RV(3)

	VR(1:3) = COORD(1:3,2) - COORD(1:3,1)
	VS(1:3) = COORD(1:3,3) - COORD(1:3,1)

      CALL VECPRD (VR,VS,VT)
      CALL SCALEN (VT,VT,DET,3)

	RV(1:3) = [1.0D0,0.0D0,0.0D0]
      CALL VECPRD (VT,RV,V1)
	ELN = SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
	IF(ELN.LT.0.0001) RV(1:3) = [0.0D0,1.0D0,0.0D0]


	IF(LXOPT.EQ.1) THEN
      CALL SCALEN (VR,VR,DET,3)
	VS(1:3) = COORD(1:3,2) - COORD(1:3,3)
      CALL SCALEN (VS,VS,DET,3)
	CALL SCAPRD (VR,RV,DR,3)
	CALL SCAPRD (VS,RV,DS,3)
	RV(1:3) = VR(1:3)
	IF(ABS(DS).GT.ABS(DR)) RV(1:3) = VS(1:3)
	ENDIF


      CALL VECPRD (VT,RV,V1)
      CALL VECPRD (V1,VT,VR)
      CALL VECPRD (VT,VR,VS)
      CALL VECPRD (VR,VS,VT)
      CALL SCALEN (VR,VR,DET,3)
      CALL SCALEN (VS,VS,DET,3)
      CALL SCALEN (VT,VT,DET,3)

	THETR = THETA*0.01745329252
	IF(ABS(THETR).GT.0.0001) THEN
	V1(1:3) = VR(1:3)
	V2(1:3) = VS(1:3)
	V1(1:3) = V1(1:3)*COS(THETR) + V2(1:3)*SIN(THETR)
      CALL VECPRD (VT,V1,VS)
      CALL VECPRD (VS,VT,VR)
      CALL VECPRD (VR,VS,VT)
      CALL SCALEN (VR,VR,DET,3)
      CALL SCALEN (VS,VS,DET,3)
      CALL SCALEN (VT,VT,DET,3)
	ENDIF


      RETURN

      END

C	=====================================================================
C	=====================================================================
C	=====================================================================