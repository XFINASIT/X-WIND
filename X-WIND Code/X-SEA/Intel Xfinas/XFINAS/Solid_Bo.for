C	***************************************************************

      SUBROUTINE XSOLIDHBD(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE
     +				,MWG,ALPHA,SEL,SEDI,FIN,HINFC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     --------------------------------------------------------------
C     MAIN PROGRAM FOR THE 3-D SOLID
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
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DP(64)        = ELASTIC OR ELASTO-PLASTIC MATERIAL MATRIX
C     H(21)         = INTERPOLATION FUNCTIONS
C     HD(3,21)      = SHAPE FUNCTION DERIVATIVES WITH RESP.TO R,S,T
C     XJ(3,3)       = JACOBIAN MATRIX
C     XJI(9)        = INVERSE OF THE JACOBIAN MATRIX
C     B(3*NNO)      = COMPRESSED STRAIN-DISPLACEMENT MATRIX
C     DISD(9)       = DISPLACEMENT DERIVATIVES
C     EPS(6)        = GAUSS POINT STRAINS
C     EPSQ(6)       = QUADRATIC PART OF GAUSS POINT STRAINS
C     SIG(6)        = GAUSS POINT STRESSES
C     IPEL          = SECTION PLASTICITY INDICATER (1=ELASTIC,2=PL)
C     RN,SN,TN      = NON-DIMENSIONAL COORDINATES
C     DET           = DETERMINANT OF THE INVERSE JACOBIAN
C     DVOL          = INTEGRATION CONSTANT
C     ----------------
C     VARIABLES OF EAS METHOD
C     ----------------
C     MM             =EAS TERMS OF ALPHA(9, 15,21 OR 30)
C     XJO(3,3)       =JACOBIAN MATRIX AT THE ORIGIN (R,S,T=0.0)
C     XJ(3,3)        =JACOBIAN MATRIX
C     TTO(6,6)       =COEFFICIENT MATRIX OF ENHANCED STRAIN
C     TM(6,MM)       =COEFFICIENT MATRIX OF ENHANCED STRAIN
C     SED(MM,MM)     =MATRIX OF EAS, TRANSPOSE(M)*E*M
C     SEL(MM,24)     =MATRIX OF EAS, TRANSPOSE(M)*E*B
C     SEDI(MM,MM)    =INVERSE OF MATRIX SED
C     EAS(6,1)       =ENHANCED STRAIN VECTOR
C     ALPHA(MM,1)    =ENHANCED TERMS ALPHA
C     TMT(MM,6)      =TRANSPOSE OF MATRIX TM(6,MM)
C     RE1(24)        =EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     RE2(24,1)      =EQUILIBRIUM FORCE OF EAS METHOD
C     ----------------------------------------------------------------

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK

	COMMON /MMENH/ MM,MM1,MM2,NDIMC
C
C	STREAM INPUT AND OUTPUT FOR WORKING ARRAY

	DIMENSION S(300)

      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,NPT),COORD(3,NNO)
      DIMENSION EDIS(*),EDISI(*),RE(*)

 
C	BASIC MATRIX FOR COMPUTATION

      DIMENSION DP(6,6),H(21),P(3,21),XJI(3,3),B(63),DISD(9)
	DIMENSION BM(6,24),DE(6,6)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6)

	DIMENSION SK(24,24),DEP(6,6),SG(300)

	DIMENSION BC(6,24)
	DIMENSION EDISC(24,1),ESTRAIN(6,1)

	DIMENSION LVNUM(8)
	

C     --------------
C       EAS METHOD
C     --------------

	DIMENSION GE(6,MM)
	DIMENSION T0(6,6),TT(6,6),TT0(6,6)   
	DIMENSION XJ(3,3),XJO(3,3),TM(6,MM)
	DIMENSION EAS(6),ALPHA(MM),RE1(24),RE2(24)
      DIMENSION SED(MM,MM),SEL(MM,24),SEDI(MM,MM),RH(MM)
	DIMENSION HINFC(MM)

	DIMENSION FIN(NEF)

C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------
      CALL HOKLAW_S (PROPM,PROPG,1)

C     ------------------------------------------------------
C	DETECTED THE THIN DIRECTION
C     ------------------------------------------------------
	XYZ1 = SQRT( (COORD(1,2)-COORD(1,1))**2.0 +
	1	         (COORD(2,2)-COORD(2,1))**2.0 +
	2	         (COORD(3,2)-COORD(3,1))**2.0 )
	XYZ2 = SQRT( (COORD(1,4)-COORD(1,1))**2.0 +
	1	         (COORD(2,4)-COORD(2,1))**2.0 +
	2	         (COORD(3,4)-COORD(3,1))**2.0 )
	XYZ3 = SQRT( (COORD(1,5)-COORD(1,1))**2.0 +
	1	         (COORD(2,5)-COORD(2,1))**2.0 +
	2	         (COORD(3,5)-COORD(3,1))**2.0 )

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

C     --------------------------------------------------------------
C     COMPUTE JACOBIAN MATRIX AT THE NATURAL CENTER (R=0.,S=0.,T=0.)
C     --------------------------------------------------------------

      CALL SHAP3D_S (0.0D0,0.0D0,0.0D0,H,P,NODEX,NNO)

	CALL JACO3D_S (COORD,P,XJO,XJI,DETO,MEL,NNO)

C     ----------------
C     COMPUTE T MATRIX (FOR EAS METHOD)
C     ---------------

	CALL MATRIXT_S (XJO,T0)
      CALL INVMATRI_S(T0,TT,6)

	TT0 = TRANSPOSE(TT)

      MGR = NGR
      MGS = NGS
      MGT = NGT
      IF (ITASK.NE.5) GOTO 10
      MGR = 3
      MGS = 3
      MGT = 3

C     -----------------------------
C     SETTING INDEX FOR GAUSS POINT
C     -----------------------------
10    IPT = 0

C     ----------------------------------------
C     INITIALIZATION OF SOME VARIABLE MATRICES
C     ----------------------------------------

	SG  = 0.0
	SK  = 0.
	SEL = 0.
	SED = 0.
	RH  = 0.
	RE1 = 0.
	RE2 = 0.
	BC  = 0.

C	------------------------
C	 LOOP OVER GAUSS POINT
C	------------------------

      DO 900  IGR=1,MGR
        RI = GLOC(IGR,MGR)
         DO 900  IGS=1,MGS
          SI = GLOC(IGS,MGS)
           DO 900  IGT=1,MGT
            TI = GLOC(IGT,MGT)

		  WT = GWT(IGR,MGR)*GWT(IGS,MGS)*GWT(IGT,MGT)
            IPT = IPT+1

C     ---------------------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ---------------------------------------------------

      CALL SHAP3D_S (RI,SI,TI,H,P,NODEX,NNO)
      CALL JACO3D_S (COORD,P,XJ,XJI,DET,MEL,NNO)
      DVOL = WT*DET

C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------

      IF (ITASK.NE.5)  GOTO 50
      CALL SOMASS_S (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
      GOTO 900

C     ----------------------------------------
C     COMPACTED STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------

50    CONTINUE

	CALL BMATSLD1 (P,XJI,B,BC,8)
	CALL DMATSLD1 (DE)
	CALL BMATANSANS (RI,SI,TI,DETO,DET,XJ,TT0,BM,P,NMX)

C     -----------------------------------------
C     FOR ENHANCED STRAIN METHOD
C     COMPUTE MATRIX [GE]
C     -----------------------------------------

	CALL SOLIDGE1(RI,SI,TI,DETO,DET,TT0,MM,GE)

C     -------------------------------------
C	COMPUTE ENHANCED STIFFNESS MATRIX SED
C	COMPUTE ENHANCED COUPLING  MATRIX SEL
C     -------------------------------------
	

	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(DE,GE))*DVOL
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(DE,BM))*DVOL


      IF (NLOPT+ITASK.EQ.1)  GOTO 700

C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------

      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,6
100   TAU(I) = WA(I,IPT)*DVOL
      GOTO 800

C     -------------------------------------
C      COMPUTE DISPLACEMENT GRADIENT (DISD)
C     -------------------------------------

200   CONTINUE

	DISD = 0.

      DO IEF=1,NEF,3

      JEF = IEF+1
      KEF = IEF+2

      DISD(1) = DISD(1) + B(IEF)*EDIS(IEF)   ! H,x *u
      DISD(2) = DISD(2) + B(JEF)*EDIS(JEF)   ! H,y *v
      DISD(3) = DISD(3) + B(KEF)*EDIS(KEF)   ! H,z *w
      DISD(4) = DISD(4) + B(JEF)*EDIS(IEF)   ! H,y *u
      DISD(5) = DISD(5) + B(KEF)*EDIS(IEF)   ! H,z *u
      DISD(6) = DISD(6) + B(IEF)*EDIS(JEF)   ! H,x *v
      DISD(7) = DISD(7) + B(KEF)*EDIS(JEF)   ! H,z *v
      DISD(8) = DISD(8) + B(IEF)*EDIS(KEF)   ! H,x *w
      DISD(9) = DISD(9) + B(JEF)*EDIS(KEF)   ! H,y *w

	END DO

C     ------------------------------
C     LINEAR COMPATIBLE STRAIN TERMS
C     ------------------------------

	STRAIN = 0.

	STRAIN(1)=0.
	STRAIN(2)=0.
	STRAIN(3)=0.
	STRAIN(4)=0.
	STRAIN(5)=0.
	STRAIN(6)=0.

	DO I=1,24
	STRAIN(1)=STRAIN(1)+BM(1,I)*EDIS(I)
      STRAIN(2)=STRAIN(2)+BM(2,I)*EDIS(I)
	STRAIN(3)=STRAIN(3)+BM(3,I)*EDIS(I)
	STRAIN(4)=STRAIN(4)+BM(4,I)*EDIS(I)
      STRAIN(5)=STRAIN(5)+BM(5,I)*EDIS(I)
	STRAIN(6)=STRAIN(6)+BM(6,I)*EDIS(I)
      ENDDO

C     -------------------------------
C     FOR EAS TERM  EAS=[GE]*{ALPHA}
C     -------------------------------

	EAS = MATMUL(GE,ALPHA)

C	-------------------------------
C	COMPUTE TOTAL COMPATIBLE STRAIN
C	-------------------------------

	DO I=1,6
      STRAIN(I)=STRAIN(I)-EAS(I)
	END DO

C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT QUADRATIC STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400

      QSTRAI(1) = .5*(DISD(1)*DISD(1)+DISD(6)*DISD(6)+DISD(8)*DISD(8))
      QSTRAI(2) = .5*(DISD(4)*DISD(4)+DISD(2)*DISD(2)+DISD(9)*DISD(9))
      QSTRAI(3) = .5*(DISD(5)*DISD(5)+DISD(7)*DISD(7)+DISD(3)*DISD(3))
      QSTRAI(4) =     DISD(1)*DISD(4)+DISD(6)*DISD(2)+DISD(8)*DISD(9)
      QSTRAI(5) =     DISD(1)*DISD(5)+DISD(6)*DISD(7)+DISD(8)*DISD(3)
      QSTRAI(6) =     DISD(4)*DISD(5)+DISD(2)*DISD(7)+DISD(9)*DISD(3)

C	----------------------
C	SUBTRACT ALMASI STRAIN
C	----------------------

      DO 390  I=1,6
390   STRAIN(I) = STRAIN(I) - QSTRAI(I)

C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------

400   GOTO (410,420,430,440,450,460,460,460),MTMOD

410   CONTINUE

C	----------------------------------------------------------
C	LINEAR ELASTIC HOOKE'S CONSTITUTUVE STRESS-STRAIN RELATION
C	----------------------------------------------------------

	CALL SOLSIG_S (STRAIN,STRESS)
      DO 415  I=1,6
415   WA(I,IPT) = STRESS(I)
      GOTO 500
C

420	CONTINUE
	
	WRITE(*,*) 'MATERIAL NUMBER IS NOT AVAILIBLE NOW'
	STOP

	GOTO 500


460   CONTINUE

c	WRITE(*,*) 'MATERIAL NUMBER IS NOT AVAILIBLE NOW'

	CALL DAMAGE(STRAIN,STRESS,DEP,PROPM,WA(1,IPT))

      GOTO 500


430	CONTINUE



	CALL DMMISE (WA(1,IPT),WA(7,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DEP)

      GOTO 500


440   CONTINUE

	CALL MOHRCRI3D(WA(1,IPT),WA(7,IPT),STRAIN,STRESS,DEP,
	1			   PROPM,WA(13,IPT))

	GO TO 500

450   CONTINUE

	CALL DRUG3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			 STRAIN,DEP,PROPM,6,0)

      GOTO 500

C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------

500   CONTINUE


	DO I=1,6
      TAU(I) = STRESS(I)*DVOL
	END DO


      IF (ITASK.LE.3) GOTO 520
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900

520	CONTINUE

C     ---------------------------------------
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     ---------------------------------------

	RE1 = RE1 + MATMUL(TRANSPOSE(BM),STRESS)*DVOL

C     ----------------------------------------
C     COMPUTE EQUILIBRIUM FORCE FOR EAS METHOD
C     ----------------------------------------

	RH = RH + MATMUL(TRANSPOSE(GE),STRESS)*DVOL

C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------

151   IF (IFREF) 900,700,900
700   IF (NLOPT.NE.0) GOTO 750

C	****************
C	LINEAR ANALYSIS
C	****************

	DEP = 0.
	SK  = 0.

	CALL DMATSLD1(DEP)

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,24
	      DO JSK =ISK,24
	         KSK = KSK+1
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	GOTO 790

C	*************************
C	...GEOMETRIC NONLINEARITY
C	*************************

750   IF (MTMOD.LE.2)THEN

	DEP =0.
	SK = 0.

	CALL DMATSLD1(DEP)

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,24
	      DO JSK =ISK,24
	         KSK = KSK+1
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	END IF

C	************************
C	...MATERIAL NONLINEARITY
C	************************

	IF (MTMOD.GT.2)THEN

	IF(NLOPT+ITASK.EQ.1) THEN


	IF (MTMOD.EQ.3) THEN

	DEP =0.

	CALL DMMISE (WA(1,IPT),WA(7,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DEP)

	ENDIF


	IF (MTMOD.EQ.4) THEN

	DEP =0.0

	CALL MOHRCRI3D(WA(1,IPT),WA(7,IPT),STRAIN,STRESS,DEP,
	1			   PROPM,WA(13,IPT))

	ENDIF


	IF (MTMOD.EQ.5) THEN

	DEP =0.0

	CALL DRUG3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			STRAIN,DEP,PROPM,6,0)

	ENDIF


	IF (MTMOD.EQ.6) THEN

	DEP =0.


	CALL DAMAGE(STRAIN,STRESS,DEP,PROPM,WA(1,IPT))

	ENDIF



	ENDIF ! <--- ENDIF NLOPT+ITASK.EQ.1



	SK = 0.

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,24
	      DO JSK =ISK,24
	         KSK = KSK+1
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	END IF

	


C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------

790   IF (NLOPT.LE.1) GOTO 810
800	CONTINUE

C	------------------------------------------------
C	ADDITIONAL COMPUTE NONLINEAR GEOMATRIC STIFFNESS
C	------------------------------------------------

	CALL KNLSTIFF(S,TAU,B)

810   CONTINUE !TIM(12) = TIM(12) + TIM2-TIM1
 
900   CONTINUE

	CALL INVMATRI_S(SED,SEDI,MM)

	CALL MATRIXES_S(SEL,SEDI,S,MM)

C	KGG = 0
C	DO I = 1,24
C	DO J = I,24
C	KGG = KGG + 1
C	S(KGG) = S(KGG) + SG(KGG)
C	ENDDO
C	ENDDO



C     ------------------------------------
C     FOR NONLINEAR ANALYSIS OF EAS METHOD
C     ------------------------------------

c  	CALL INVMATRI_S(SED,SEDI,MM)

	RE2 =  MATMUL(MATMUL(TRANSPOSE(SEL),SEDI),RH)

C     ------------------------------
C     FOR ENHANCED STRAIN METHOD
C     K=TRANSPOSE(L)*INVERSE(D)*L
C     K(FINAL)=K(PURE DISP.)-K(EAS)
C     ------------------------------

c      CALL MATRIXES_S(SEL,SEDI,S,MM)

	IF (ITASK+NLOPT.EQ.1) GOTO 150

	DO I=1,24
	 RE(I)=RE1(I)-RE2(I)
	END DO

C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)

150	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF

	

      RETURN
      END

C	*******************************************************

	SUBROUTINE BMATANSANS(RI,SI,TI,DETO,DET,XJ,TT0,BM,P,NMX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	-----------------------------------------

C		ASSUMED NATURAL STRAIN SUBROUTINE
C		STRAIN-DISPLACEMENT RELATIONSHIP
C		ASSUMED IN e11,e22,e33,2e12,2e13,2e23

C	-----------------------------------------

	DIMENSION XJ(3,3)
	DIMENSION BM(6,24),HXI(8),HET(8),HZE(8),XJJ(24,24)
	DIMENSION TT0(6,6)
	DIMENSION P(3,8)
	DIMENSION B(63),TT0B(3,3),BB(3,24),XJJB(24,24)

	DIMENSION BA(6,24)

C	----------------
C	STRAT SUBROUTINE
C	----------------

	BM  = 0.
	HXI = 0.
	HET = 0.
	HZE = 0.
	XJJ = 0.
	B   = 0.

C	------------------------------
C	COMPATIBLE FROM NATURAL STRAIN
C	------------------------------

	DO I =1,8
	 BM(1,3*I-2) = P(1,I)
	 BM(2,3*I-1) = P(2,I)
	 BM(3,3*I-0) = P(3,I)
	END DO

	DO I =1,8
	 BM(4,3*I-2) = 0.
	 BM(4,3*I-1) = 0.
	 BM(4,3*I-0) = 0
	END DO

	DO I =1,8
	 BM(4,3*I-2) = P(2,I)
	 BM(4,3*I-1) = P(1,I)
	END DO

C	***********************************************
C		ASSUMED NATURAL TRANSVERSE SHEAR STRAIN
C	***********************************************

C	---------------------------------
C		TRANSVERSE SHEAR STRAIN 13
C	---------------------------------

C	------------------------------
C		SAMPLING POINT (0,-1,0)
C	------------------------------

	  XI   =  0.
	  ETA  = -1.
	  ZETA =  0.


	HZE(1) =  (1./8.)*(1.+XI)*(1.+ETA)
	HZE(2) =  (1./8.)*(1.-XI)*(1.+ETA)
	HZE(3) =  (1./8.)*(1.-XI)*(1.-ETA)
	HZE(4) =  (1./8.)*(1.+XI)*(1.-ETA)
	HZE(5) = -(1./8.)*(1.+XI)*(1.+ETA)
	HZE(6) = -(1./8.)*(1.-XI)*(1.+ETA)
	HZE(7) = -(1./8.)*(1.-XI)*(1.-ETA)
	HZE(8) = -(1./8.)*(1.+XI)*(1.-ETA)

	DO I=1,8
	 BM(5,3*I-2) = BM(5,3*I-2)+(1./2.)*(1.-SI)*HZE(I)
	END DO

C	------------------------------
C		SAMPLING POINT (0,-1,0)
C	------------------------------

	  XI   =  0.
	  ETA  =  1.
	  ZETA =  0.


	HZE(1) =  (1./8.)*(1.+XI)*(1.+ETA)
	HZE(2) =  (1./8.)*(1.-XI)*(1.+ETA)
	HZE(3) =  (1./8.)*(1.-XI)*(1.-ETA)
	HZE(4) =  (1./8.)*(1.+XI)*(1.-ETA)
	HZE(5) = -(1./8.)*(1.+XI)*(1.+ETA)
	HZE(6) = -(1./8.)*(1.-XI)*(1.+ETA)
	HZE(7) = -(1./8.)*(1.-XI)*(1.-ETA)
	HZE(8) = -(1./8.)*(1.+XI)*(1.-ETA)

	DO I=1,8
	 BM(5,3*I-2) = BM(5,3*I-2)+(1./2.)*(1.+SI)*HZE(I)
	END DO

C	------------------------------
C		SAMPLING POINT (0,-1,0)
C	------------------------------

	  XI   =  0.
	  ETA  = -1.
	  ZETA =  0.

	HXI(1) =  (1./8.)*(1.+ETA)*(1.+ZETA)
	HXI(2) = -(1./8.)*(1.+ETA)*(1.+ZETA)
	HXI(3) = -(1./8.)*(1.-ETA)*(1.+ZETA)
	HXI(4) =  (1./8.)*(1.-ETA)*(1.+ZETA)
	HXI(5) =  (1./8.)*(1.+ETA)*(1.-ZETA)
	HXI(6) = -(1./8.)*(1.+ETA)*(1.-ZETA)
	HXI(7) = -(1./8.)*(1.-ETA)*(1.-ZETA)
	HXI(8) =  (1./8.)*(1.-ETA)*(1.-ZETA)

	DO I=1,8
	 BM(5,3*I-0) = BM(5,3*I-0)+(1./2.)*(1.-SI)*HXI(I)
	END DO

C	------------------------------
C		SAMPLING POINT (0,1,0)
C	------------------------------

	  XI   =  0.
	  ETA  =  1.
	  ZETA =  0.

	HXI(1) =  (1./8.)*(1.+ETA)*(1.+ZETA)
	HXI(2) = -(1./8.)*(1.+ETA)*(1.+ZETA)
	HXI(3) = -(1./8.)*(1.-ETA)*(1.+ZETA)
	HXI(4) =  (1./8.)*(1.-ETA)*(1.+ZETA)
	HXI(5) =  (1./8.)*(1.+ETA)*(1.-ZETA)
	HXI(6) = -(1./8.)*(1.+ETA)*(1.-ZETA)
	HXI(7) = -(1./8.)*(1.-ETA)*(1.-ZETA)
	HXI(8) =  (1./8.)*(1.-ETA)*(1.-ZETA)

	DO I=1,8
	 BM(5,3*I-0) = BM(5,3*I-0)+(1./2.)*(1.+SI)*HXI(I)
	END DO

C	---------------------------------
C		TRANSVERSE SHEAR STRAIN 23
C	---------------------------------

C	------------------------------
C		SAMPLING POINT (-1,0,0)
C	------------------------------

	  XI   = -1.
	  ETA  =  0.
	  ZETA =  0.


	HZE(1) =  (1./8.)*(1.+XI)*(1.+ETA)
	HZE(2) =  (1./8.)*(1.-XI)*(1.+ETA)
	HZE(3) =  (1./8.)*(1.-XI)*(1.-ETA)
	HZE(4) =  (1./8.)*(1.+XI)*(1.-ETA)
	HZE(5) = -(1./8.)*(1.+XI)*(1.+ETA)
	HZE(6) = -(1./8.)*(1.-XI)*(1.+ETA)
	HZE(7) = -(1./8.)*(1.-XI)*(1.-ETA)
	HZE(8) = -(1./8.)*(1.+XI)*(1.-ETA)

	DO I=1,8
	 BM(6,3*I-1) = BM(6,3*I-1)+(1./2.)*(1.-RI)*HZE(I)
	END DO

C	------------------------------
C		SAMPLING POINT (1,0,0)
C	------------------------------

	  XI   =  1.
	  ETA  =  0.
	  ZETA =  0.


	HZE(1) =  (1./8.)*(1.+XI)*(1.+ETA)
	HZE(2) =  (1./8.)*(1.-XI)*(1.+ETA)
	HZE(3) =  (1./8.)*(1.-XI)*(1.-ETA)
	HZE(4) =  (1./8.)*(1.+XI)*(1.-ETA)
	HZE(5) = -(1./8.)*(1.+XI)*(1.+ETA)
	HZE(6) = -(1./8.)*(1.-XI)*(1.+ETA)
	HZE(7) = -(1./8.)*(1.-XI)*(1.-ETA)
	HZE(8) = -(1./8.)*(1.+XI)*(1.-ETA)

	DO I=1,8
	 BM(6,3*I-1) = BM(6,3*I-1)+(1./2.)*(1.+RI)*HZE(I)
	END DO


C	------------------------------
C		SAMPLING POINT (-1,0,0)
C	------------------------------

	  XI   = -1.
	  ETA  =  0.
	  ZETA =  0.

	HET(1) =   (1./8.)*(1.+XI)*(1.+ZETA)
	HET(2) =   (1./8.)*(1.-XI)*(1.+ZETA)
	HET(3) =  -(1./8.)*(1.-XI)*(1.+ZETA)
	HET(4) =  -(1./8.)*(1.+XI)*(1.+ZETA)
	HET(5) =   (1./8.)*(1.+XI)*(1.-ZETA)
	HET(6) =   (1./8.)*(1.-XI)*(1.-ZETA)
	HET(7) =  -(1./8.)*(1.-XI)*(1.-ZETA)
	HET(8) =  -(1./8.)*(1.+XI)*(1.-ZETA)

	DO I=1,8
	 BM(6,3*I-0) = BM(6,3*I-0)+(1./2.)*(1.-RI)*HET(I)
	END DO

C	------------------------------
C		SAMPLING POINT (1,0,0)
C	------------------------------

	  XI   =  1.
	  ETA  =  0.
	  ZETA =  0.

	HET(1) =   (1./8.)*(1.+XI)*(1.+ZETA)
	HET(2) =   (1./8.)*(1.-XI)*(1.+ZETA)
	HET(3) =  -(1./8.)*(1.-XI)*(1.+ZETA)
	HET(4) =  -(1./8.)*(1.+XI)*(1.+ZETA)
	HET(5) =   (1./8.)*(1.+XI)*(1.-ZETA)
	HET(6) =   (1./8.)*(1.-XI)*(1.-ZETA)
	HET(7) =  -(1./8.)*(1.-XI)*(1.-ZETA)
	HET(8) =  -(1./8.)*(1.+XI)*(1.-ZETA)

	DO I=1,8
	 BM(6,3*I-0) = BM(6,3*I-0)+(1./2.)*(1.+RI)*HET(I)
	END DO

C	-----------------------------------
C		COMPUTE JACOBIAN DETERMINANT
C	-----------------------------------

      DO I=1,8
      XJJ(3*i-2,3*i-2) =XJ(1,1)
      XJJ(3*i-2,3*i-1) =XJ(1,2)
      XJJ(3*i-2,3*i)   =XJ(1,3)
      XJJ(3*i-1,3*i-2) =XJ(2,1)
      XJJ(3*i-1,3*i-1) =XJ(2,2)
      XJJ(3*i-1,3*i)   =XJ(2,3)
      XJJ(3*i,3*i-2)   =XJ(3,1)
      XJJ(3*i,3*i-1)   =XJ(3,2)
      XJJ(3*i,3*i)     =XJ(3,3)
	ENDDO


	BM = MATMUL(MATMUL(TT0,BM),XJJ)*(DETO/DET)


	RETURN
	END

C	**************************************************

	SUBROUTINE SOLIDGE1(R,S,T,DETJO,DETJ,TT,MM,GE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	-----------------------------------------------
C	SUBROUTINE COMPUTE THE COMPATIBLE ENHANCED TERM
C	-----------------------------------------------

	DIMENSION GM(6,MM),GE(6,MM),TT(6,6),EM(6,MM)

	DO I =1,6
	 DO J=1,MM
		GM(I,J) = 0.
	    GE(I,J) = 0.
		EM(I,J) = 0.
	 END DO
	END DO

	GM(1,1) =  -2.0*R
	GM(2,2) =  -2.0*S
	GM(3,3) =  -2.0*T
	GM(4,4) =  -2.0*R
	GM(4,5) =  -2.0*S

	IF(MM.EQ.5)GOTO 100

	GM(5,6) =  -2.0*R
	GM(5,7) =  -2.0*T
	GM(6,8) =  -2.0*S
	GM(6,9) =  -2.0*T


	IF(MM.EQ.9)GOTO 100

	GM(4,10) = -3.0*R*S
	GM(4,11) = -3.0*S*T
	GM(4,12) = -3.0*R*T

	GM(5,13) = -3.0*R*S
	GM(5,14) = -3.0*S*T
	GM(5,15) = -3.0*R*T

	GM(6,16) = -3.0*R*S
	GM(6,17) = -3.0*S*T
	GM(6,18) = -3.0*R*T

	GM(1,19) = -3.0*R*T
	GM(1,20) = -3.0*R*S

	GM(2,21) = -3.0*R*S
	GM(2,22) = -3.0*S*T

	GM(3,23) = -3.0*S*T
	GM(3,24) = -3.0*R*T

	IF(MM.EQ.24)GOTO 100

100	CONTINUE

	EM = MATMUL(TT,GM)
      CONST = (DETJO/DETJ)

C	COMPUTE [GE] MATRIX

	DO IEM = 1,6
		DO  JEM = 1,MM
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
		END DO
	END DO

	RETURN
	END

C	*****************************************************************

	SUBROUTINE INVMATRI_S(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------
C     MATRIX INVERSION BY ELIMMINATION WITH PARTIAL PIVOTING
C     AND DETERMINANT OF MATRIX A
C     ROGINAL MATRIX =A, INVERSE MATRIX=B
C
C     A: INPUT SQUARE MATRIX
C     B: INVERSE MATRIX OF A
C     N: SIZE OF A AND B
C     EPS: CONTROL VARIABLE
C     DET: DETERMINANT OF A
C     --------

      DIMENSION A(N,N),B(N,N),C(N,N)

C     CONSTRUCT IDENTITY MATRIX B(I,J)=I
C	MODIFIED SINGULARITY

	EPS=0.0000000000001


	DO 88 I=1,N
	DO 88 J=1,N
	C(I,J)=A(I,J)
  88  CONTINUE
	DO 6 I=1,N
	DO 5 J=1,N
	IF(I-J) 4,3,4
   3  B(I,J)=1.0
      GOTO 5
   4  B(I,J)=0.0
   5  CONTINUE
   6  CONTINUE

C     LOCATE MAXIMUM MAGNITUDE A(I,K) ON OR BELOW MAIN DIAGONAL
      DET=1.0
	DO 45 K=1,N
	IF (K-N) 12,30,30
   12 IMAX=K
      AMAX=ABS(C(K,K))
	KP1=K+1
	DO 20 I=KP1,N
	IF (AMAX-ABS(C(I,K))) 15,20,20
   15 IMAX=I
      AMAX=ABS(C(I,K))
   20 CONTINUE

C     INTERCHANGE ROWS IMAX AND K IF IMAX NOT EQUAL TO K
      IF (IMAX-K) 25,30,25
   25 DO 29 J=1,N
      ATMP=C(IMAX,J)
	C(IMAX,J)=C(K,J)
	C(K,J)=ATMP
	BTMP=B(IMAX,J)
	B(IMAX,J)=B(K,J)
   29 B(K,J)=BTMP
      DET=-DET
   30 CONTINUE

C     TEST FOR SINGULAR MATRIX
      IF (ABS(C(K,K))-EPS) 33,33,35
   33 WRITE(*,*) 'SINGULAR MATRIX EPS- INSIDE SUB INVMATRIX'
	STOP
   35 DET=C(K,K)*DET

C     DIVIDE PIVOT ROW BY ITS MAIN DIAGONAL ELEMENT
      DIV=C(K,K)
	DO 38 J=1,N
	C(K,J)=C(K,J)/DIV
   38 B(K,J)=B(K,J)/DIV

C     REPLACE EACH ROW BY LINEAR COMBINATION WITH PIVOT ROW
      DO 43 I=1,N
	AMULT=C(I,K)
	IF (I-K) 39,43,39
   39 DO 42 J=1,N
      C(I,J)=C(I,J)-AMULT*C(K,J)
   42 B(I,J)=B(I,J)-AMULT*B(K,J)
   43 CONTINUE
   45 CONTINUE
C
      RETURN
	END


C=====================================================================

      SUBROUTINE BMATSLD1 (P,XJI,B,BC,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     EVALUATES THE GLOBAL LINEAR STRAIN-DISPLACEMENT MATRIX FOR A
C     CURVILINEAR ISOPARAMETRIC HEXAHEDRON (8 TO 21 NODES)
C	----------------------------------------------------
C     INPUT OUTPUT VARIABLES
C	----------------------
C     P(3,NNO)    = SHAPE FUNCTION DERIVATIVES WITH RESPECT TO R,S,T
C     XJI(3,3)    = INVERSE OF THE JACOBIAN MATRIX
C     B(3*NNO)    = COMPRESSED LINEAR STRAIN DISPLACEMENT MATRIX
C     NNO         = NUMBER OF NODES USED TO DESCRIBE THIS ELEMENT
C     ----------------------------------------------------------------

      DIMENSION P(3,21),XJI(3,3),B(63),BC(6,24)
C
      CALL CLEARA (B,63)

C	------------------------------------------------
C	 COMPUTE COMPACTED STRESS-STRAIN GRADIENT VECTOR
C	------------------------------------------------

	BC = 0.

      K = 1
	KK= 0

      DO 100  INO=1,NNO
      DO 50   I=1,3
      B(K)   = B(K)   + XJI(1,I)*P(I,INO)
      B(K+1) = B(K+1) + XJI(2,I)*P(I,INO)
      B(K+2) = B(K+2) + XJI(3,I)*P(I,INO)
50    CONTINUE
      K = K+3
	KK=KK+3
100   CONTINUE
C

C	----------------------------------------
C	 COMPUTE COMPACTED STRESS-STRAIN MATRIX
C	----------------------------------------

	DO I=1,6
	 DO J =1,24
	  BC(I,J)=0.
	 END DO
	END DO



	DO INO=1,8

	HX =0.
	HY =0.
	HZ =0.

	DO I=1,3
        HX  = HX  + XJI(1,I)*P(I,INO)
        HY  = HY  + XJI(2,I)*P(I,INO)
        HZ  = HZ  + XJI(3,I)*P(I,INO)
	END DO

C	--------------------------
C  	 COMPUTE MATRIX [B](6,24)
C	--------------------------

	BC(1,3*INO-2) = HX
	BC(2,3*INO-1) = HY
	BC(3,3*INO-0) = HZ

	BC(4,3*INO-2) = HY
	BC(4,3*INO-1) = HX

	BC(5,3*INO-2) = HZ
	BC(5,3*INO-0) = HX

	BC(6,3*INO-1) = HZ
	BC(6,3*INO-0) = HY

	END DO




      RETURN
      END

C	***************************************************************

	SUBROUTINE DMATSLD1(D)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	**********************************
C	COMPUTE LINEAR CONSTITUTIVE MATRIX
C	**********************************

      COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST

	DIMENSION D(6,6)

	DO I =1,6
	 DO J =1,6
	  D(I,J) = 0.
	 END DO
	END DO

C	LINEAR CONSITITUTIVE

      D(1,1)  = A1
      D(1,2)  = B1
      D(1,3)  = B1
      D(2,1)  = B1
      D(2,2)  = A1
      D(2,3)  = B1
      D(3,1)  = B1
      D(3,2)  = B1
      D(3,3)  = A1

      D(4,4) = C1
      D(5,5) = C1
      D(6,6) = C1

      RETURN
      END

C	**************************************************


C=====================================================================


C=====================================================================
	SUBROUTINE HOKLAW_S(PROPM,PROPG,IND)
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

      COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST

c	COMMON /CRIP_S/  OCR,VO,DEN,VGL,SWL,SM,OK

c	COMMON /MOHR_S/  COH,ANG,BETA
C
      DIMENSION PROPM(5),PROPG(*)
C
      YM  = PROPM(1)
      PR  = PROPM(2)
      YLD = PROPM(3)
      HP  = PROPM(4)
      DEN = PROPM(5)

C	CRITICAL STATE MATERIAL

	OCR = PROPM(3)
	VO  = PROPM(4)
	VGL = PROPM(6)
	SWL = PROPM(7)
	SM  = PROPM(8)
	OK  = PROPM(9)

C	MOHR COULOMB

	COH  = PROPM(10)
	ANG  = PROPM(11)
	BETA = PROPM(12)


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
200   A1 = YM/(1.-PR*PR)
      B1 = A1*PR
C     -----------------------
C     NODAL ELEMENT THICKNESS
C     -----------------------
500   TH  = PROPG(2)
C
      RETURN

C	RETURN A1,B1,C1
C
C	RE ARRANGE TO |  A1 B1 0  | AT SUBROUTINE DMATRIX
C		          |  B1 A1 0  |	[DM]
C		          |  0  0  C1 |

      END

C
C=====================================================================
C	==============================================================
      SUBROUTINE MATRIXT_S(XJO,TTO)
C	--------
C     FOR ENHANCED STRAIN METHOD
C     COMPUTE MATRIX T0 (HERE SET TTO) THE TRANSFORM MATRIX
C     XJO  : JACOBIAN MATRIX AT THE ORIGIN (0,0,0)
C	--------
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION TTO(6,6),XJO(3,3)
	DO 5 I=1,6
	DO 5 J=1,6
      TTO(I,J)=0.0
5     CONTINUE
	DO 10 I=1,3
	DO 10 J=1,3
	TTO(I,J)=XJO(J,I)**2
  10  CONTINUE

      TTO(4,1)=XJO(1,1)*XJO(1,2)
	TTO(5,1)=XJO(1,1)*XJO(1,3)
	TTO(6,1)=XJO(1,2)*XJO(1,3)

	TTO(4,2)=XJO(2,1)*XJO(2,2)
	TTO(5,2)=XJO(2,1)*XJO(2,3)
	TTO(6,2)=XJO(2,2)*XJO(2,3)

	TTO(4,3)=XJO(3,1)*XJO(3,2)
	TTO(5,3)=XJO(3,1)*XJO(3,3)
	TTO(6,3)=XJO(3,2)*XJO(3,3)

	TTO(1,4)=2*XJO(1,1)*XJO(2,1)
	TTO(1,5)=2*XJO(1,1)*XJO(3,1)
	TTO(1,6)=2*XJO(2,1)*XJO(3,1)

	TTO(2,4)=2*XJO(1,2)*XJO(2,2)
	TTO(2,5)=2*XJO(1,2)*XJO(3,2)
	TTO(2,6)=2*XJO(2,2)*XJO(3,2)

	TTO(3,4)=2*XJO(1,3)*XJO(2,3)
	TTO(3,5)=2*XJO(1,3)*XJO(3,3)
	TTO(3,6)=2*XJO(2,3)*XJO(3,3)

      TTO(4,4)=XJO(1,1)*XJO(2,2)+XJO(2,1)*XJO(1,2)
	TTO(4,5)=XJO(1,1)*XJO(3,2)+XJO(3,1)*XJO(1,2)
	TTO(4,6)=XJO(2,1)*XJO(3,2)+XJO(3,1)*XJO(2,2)

	TTO(5,4)=XJO(1,1)*XJO(2,3)+XJO(2,1)*XJO(1,3)
	TTO(5,5)=XJO(1,1)*XJO(3,3)+XJO(3,1)*XJO(1,3)
      TTO(5,6)=XJO(2,1)*XJO(3,3)+XJO(3,1)*XJO(2,3)

	TTO(6,4)=XJO(1,2)*XJO(2,3)+XJO(2,2)*XJO(1,3)
	TTO(6,5)=XJO(1,2)*XJO(3,3)+XJO(3,2)*XJO(1,3)
	TTO(6,6)=XJO(2,2)*XJO(3,3)+XJO(3,2)*XJO(2,3)
	RETURN
      END

C	======================================================
C=====================================================================
      SUBROUTINE SHAP3D_S (R,S,T,H,P,NODEX,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     AT THE NODAL POINTS OF A CURVILINEAR ISOPARAMETRIC HEXAHEDRON
C     (8 TO 21 NODES)
C	---------------
C                    NODE NUMBERING CONVENTION
C				   -------------------------
C                                  2
C                                 0
C                               . . .                   T.     S.
C                             .18 0   .                  .    .
C                       10  .     .     .                .  .
C                         0       0       .  9           .
C                       .       .  6 .      0              .
C                     .       .        .      .              .
C                3  .       .            .      .              .R
C                 0       0                0      .
C                 . .   .  14            13  .      .  1
C              19 0   .                        .      0
C                 .     .           0            .  . .
C                 0       .  11   21              .   0 17
C                7  .       0               12  .     .
C                     .       .               0       0
C                       .       .           .       .  5
C                         .       .       .       .
C                           0       .4  .       .
C                         15  .       0       0
C                               .     .     . 16
C                                 .   0 20.
C                                   . . .
C                                     0
C                                      8
C	--------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     R,S,T      = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     H(NNO)     = INTERPOLATION (SHAPE) FUNCTIONS
C     P(3,NNO)   = FUNCTION DERIVATIVES WITH RESPECT TO R,S,T RESP.
C     NODEX(NEX) = POSITIONS OF MIDSIDE NODES (EXCESS NODES)
C     NNO        = NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ----------------------------------------------------------------
      DIMENSION H(21),P(3,21),NODEX(1),IPERM(8)
C
      DATA IPERM  /2,3,4,1,6,7,8,5/
C
      NMI = NNO-8
      RP  = 1.+R
      SP  = 1.+S
      TP  = 1.+T
      RM  = 1.-R
      SM  = 1.-S
      TM  = 1.-T
      RR  = 1.-R*R
      SS  = 1.-S*S
      TT  = 1.-T*T
C     ----------------------------------------------------------
C     INTERPOLATION FUNCTIONS AND DERIVATIVES FOR A 8 NODE BRICK
C     ----------------------------------------------------------
	H(1)   = .125*RP*SP*TP
      H(2)   = .125*RM*SP*TP
      H(3)   = .125*RM*SM*TP
      H(4)   = .125*RP*SM*TP
      H(5)   = .125*RP*SP*TM
      H(6)   = .125*RM*SP*TM
      H(7)   = .125*RM*SM*TM
      H(8)   = .125*RP*SM*TM
C
      P(1,1) = .125*SP*TP
      P(1,2) = -P(1,1)
      P(1,3) = -.125*SM*TP
      P(1,4) = -P(1,3)
      P(1,5) = .125*SP*TM
      P(1,6) = -P(1,5)
      P(1,7) = -.125*SM*TM
      P(1,8) = -P(1,7)
C
      P(2,1) = .125*RP*TP
      P(2,2) = .125*RM*TP
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)
      P(2,5) = .125*RP*TM
      P(2,6) = .125*RM*TM
      P(2,7) = -P(2,6)
      P(2,8) = -P(2,5)
C
      P(3,1) = .125*RP*SP
      P(3,2) = .125*RM*SP
      P(3,3) = .125*RM*SM
      P(3,4) = .125*RP*SM
      P(3,5) = -P(3,1)
      P(3,6) = -P(3,2)
      P(3,7) = -P(3,3)
      P(3,8) = -P(3,4)

      IF (NNO.EQ.8)  RETURN
C     -------------------------------------
C     ADD DEGREES OF FREEDOM IN EXCESS OF 8
C     -------------------------------------
      DO 100  IMI=1,NMI
      NN = NODEX(IMI)-8
      GOTO (9,10,11,12,13,14,15,16,17,18,19,20,21),NN
C      CALL GOTOER !NON-EXIST SUBROUTINE
C
9     H(9)    = .25*RR*SP*TP
      P(1,9)  =-.50*R*SP*TP
      P(2,9)  = .25*RR*TP
      P(3,9)  = .25*RR*SP
      GOTO 100
C
10    H(10)   = .25*RM*SS*TP
      P(1,10) =-.25*SS*TP
      P(2,10) =-.50*RM*S*TP
      P(3,10) = .25*RM*SS
      GOTO 100
C
11    H(11)   = .25*RR*SM*TP
      P(1,11) =-.50*R*SM*TP
      P(2,11) =-.25*RR*TP
      P(3,11) = .25*RR*SM
      GOTO 100
C
12    H(12)   = .25*RP*SS*TP
      P(1,12) = .25*SS*TP
      P(2,12) =-.50*RP*S*TP
      P(3,12) = .25*RP*SS
      GOTO 100
C
13    H(13)   = .25*RR*SP*TM
      P(1,13) =-.50*R*SP*TM
      P(2,13) = .25*RR*TM
      P(3,13) =-.25*RR*SP
      GOTO 100
C
14    H(14)   = .25*RM*SS*TM
      P(1,14) =-.25*SS*TM
      P(2,14) =-.50*RM*S*TM
      P(3,14) =-.25*RM*SS
      GOTO 100
C
15    H(15)   = .25*RR*SM*TM
      P(1,15) =-.50*R*SM*TM
      P(2,15) =-.25*RR*TM
      P(3,15) =-.25*RR*SM
      GOTO 100
C
16    H(16)   = .25*RP*SS*TM
      P(1,16) = .25*SS*TM
      P(2,16) =-.50*RP*S*TM
      P(3,16) =-.25*RP*SS
      GOTO 100
C 
17    H(17)   = .25*RP*SP*TT
      P(1,17) = .25*SP*TT
      P(2,17) = .25*RP*TT
      P(3,17) =-.50*RP*SP*T
      GOTO 100
C
18    H(18)   = .25*RM*SP*TT
      P(1,18) =-.25*SP*TT
      P(2,18) = .25*RM*TT
      P(3,18) =-.50*RM*SP*T
      GOTO 100
C
19    H(19)   = .25*RM*SM*TT
      P(1,19) =-.25*SM*TT
      P(2,19) =-.25*RM*TT
      P(3,19) =-.50*RM*SM*T
      GOTO 100
C
20    H(20)   = .25*RP*SM*TT
      P(1,20) = .25*SM*TT
      P(2,20) =-.25*RP*TT
      P(3,20) =-.50*RP*SM*T
      GOTO 100
C
21    H(21)   = RR*SS*TT
      P(1,21) =-2.*R*SS*TT
      P(2,21) =-2.*S*RR*TT
      P(3,21) =-2.*T*RR*SS
100   CONTINUE
C     ------------------------------------------------------
C     MODIFY FIRS 8 FUNCTIONS IF 9 OR MORE NODES ARE PRESENT
C     ------------------------------------------------------
      DO 290  IMI=1,NMI
      IN = NODEX(IMI)
      IF (IN.GT.16)  GOTO 210
      I1 = IN-8
      I2 = IPERM(I1)
      GOTO 250
210   IF (IN.EQ.21)  GOTO 300
      I1 = IN-16
      I2 = I1+4
C
250   H(I1) = H(I1) - .5*H(IN)
      H(I2) = H(I2) - .5*H(IN)
      H(IMI+8) = H(IN)
      DO 290  J=1,3
      P(J,I1) = P(J,I1) - .5*P(J,IN)
      P(J,I2) = P(J,I2) - .5*P(J,IN)
290   P(J,IMI+8) = P(J,IN)
      RETURN
C     -----------------------------------------------
C     MODIFY FIRST 20 FUNCTIONS IF NODE 21 IS PRESENT
C     -----------------------------------------------
300   DO 390  IMI=1,NMI
      IN = NODEX(IMI)
      IF (IN.GT.16)  GOTO 310
      I1 = IN-8
      I2 = IPERM(I1)
      GOTO 350
310   IF (IN.EQ.21)  GOTO 400
      I1 = IN-16
      I2 = I1+4
C
350   H(I1) = H(I1) + .125*H(21)
      H(I2) = H(I2) + .125*H(21)
      DO 390  J=1,3
      P(J,I1) = P(J,I1) + .125*P(J,21)
390   P(J,I2) = P(J,I2) + .125*P(J,21)
C
400   DO 410  I=1,8
      H(I) = H(I) - .125*H(21)
      DO 410  J=1,3
410   P(J,I) = P(J,I) - .125*P(J,21)
      NN = NMI+7
      IF (NN.EQ.8)  RETURN
      DO 420  I=9,NN
      H(I) = H(I) - .25*H(21)
      DO 420  J=1,3
420   P(J,I) = P(J,I) - .25*P(J,21)
      H(NMI+8) = H(21)
      DO 430  J=1,3
430   P(J,NMI+8) = P(J,21)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE JACO3D_S (XY,P,XJ,XJI,DET,MEL,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     EVALUATES THE JACOBIAN MATRIX AT POINT (R,S,T),ITS DETERMINANT
C     AND INVERSE FOR A CURVILINEAR ISOPARAMETRIC HEXAHEDRON (8-21 NO)
C	----------------------------------------------------------------
C     INPUT OUTPUT VARIABLES
C	----------------------
C     XY(3,NNO)    = COORDINATES OF ELEMENT NODES
C     P(3,NNO)     = SHAPE FUNCTION DERIVATIVES WITH RESPECT TO R,S,T
C     XJI(3,3)     = INVERSE OF THE JACOBIAN MATRIX
C     DET          = DETERMINANT OF THE JACOBIAN MATRIX
C     MEL          = ELEMENT NUMBER
C     NNO          = NUMBER OF NODES USED TO DESCRIBE THIS ELEMENT
C     ----------------------------------------------------------------
      DIMENSION XY(3,NNO),P(3,21),XJ(3,3),AJ(9),XJI(9),INSONG(8)
C
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------


	INSONG(1) = 1
	INSONG(2) = 2
	INSONG(3) = 3
	INSONG(4) = 4
	INSONG(5) = 5
	INSONG(6) = 6
	INSONG(7) = 7
	INSONG(8) = 8

      DO 100  I=1,3
      DO 100  J=1,3
      DUM = 0.
      DO 90  INO=1,NNO
	ISNG = INSONG(INO)
90    DUM = DUM + P(I,INO)*XY(J,ISNG)
100   XJ(I,J) = DUM
      K=1
      DO 110 J=1,3
      DO 110 I=1,3
	AJ(K)=XJ(I,J)
	K=K+1
110   CONTINUE

	
C     ----------------------------------------
C     DETERMINANT OF THE JACOBIAN MATRIX (DET)
C     ----------------------------------------
      DET = AJ(1)*AJ(5)*AJ(9) + AJ(4)*AJ(8)*AJ(3) + AJ(7)*AJ(2)*AJ(6)
     1    - AJ(7)*AJ(5)*AJ(3) - AJ(4)*AJ(2)*AJ(9) - AJ(1)*AJ(8)*AJ(6)

      IF (DET.LT.1.E-8) CALL ERRORS (15,H,MEL,'JACOB.DET.')
C     ------------------------------------
C     INVERSE OF THE JACOBIAN MATRIX (XJI)
C     ------------------------------------
      DUM = 1./DET
      XJI(1) = DUM*( AJ(5)*AJ(9) - AJ(8)*AJ(6))
      XJI(2) = DUM*(-AJ(2)*AJ(9) + AJ(8)*AJ(3))
      XJI(3) = DUM*( AJ(2)*AJ(6) - AJ(5)*AJ(3))
      XJI(4) = DUM*(-AJ(4)*AJ(9) + AJ(7)*AJ(6))
      XJI(5) = DUM*( AJ(1)*AJ(9) - AJ(7)*AJ(3))
      XJI(6) = DUM*(-AJ(1)*AJ(6) + AJ(4)*AJ(3))
      XJI(7) = DUM*( AJ(4)*AJ(8) - AJ(7)*AJ(5))
      XJI(8) = DUM*(-AJ(1)*AJ(8) + AJ(7)*AJ(2))
      XJI(9) = DUM*( AJ(1)*AJ(5) - AJ(4)*AJ(2))
C
      RETURN
      END
C
C=====================================================================
C=====================================================================
      SUBROUTINE SOMASS_S (CM,H,DEN,DVOL,NNO,NEF,IPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------
C     CALCULATES ELEMENT MASS MATRIX FOR THE 3-D CURVILINEAR,
C     ISOPARAMETRIC HEXAHEDRON (8 TO 21 NODES)
C	----------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     CM(NWS)    = ELEMENT MASS MATRIX (UPPER TRIANG. ROW-WISE)
C     H(NNO)     = SHAPE FUNCTIONS
C     DEN        = DENSITY
C     DVOL       = INTEGRATION CONSTANT (WT*DET)
C     NNO        = NUMBER OF NODES USED TO DESCRIBE THIS ELEMENT
C     NEF        = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     IPT        = GAUSS POINT COUNTER
C     ----------------------------------------------------------
      COMMON /DYNA/ CDEN,IMASS

      DIMENSION CM(1),H(21)
C     --------------
C     INITIALISATION
C     --------------
      IF (IPT.GT.1) GOTO 50
      VOL = 0.
50    DVOL = DVOL*DEN
      IF (IMASS-2) 100,100,300
C     ----------------------------------------
C     INTEGRATE LUMPED OR DIAGONAL MASS MATRIX
C     ----------------------------------------
100   VOL = VOL + DVOL
      IF (IMASS.EQ.1) GOTO 500
C
      DO 250  INO=1,NNO
      KM = 3*(INO-1) + 1
250   CM(KM) = CM(KM) + H(INO)*H(INO)
      GOTO 500
C     --------------------------------
C     INTEGRATE CONSISTENT MASS MATRIX
C     --------------------------------
300   KM = 1
      DO 390  INO=1,NNO
      DO 350  JNO=INO,NNO
      CM(KM) = CM(KM) + H(INO)*H(JNO)*DVOL
350   KM = KM+3
390   KM = KM + 6*(NNO-INO) + 3
C     --------------------------------------------
C     FILL IN REMAINING TERMS IN Y AND Z DIRECTION
C     --------------------------------------------
500   IF (IPT.LT.27) RETURN
      IF (IMASS-2) 600,700,800
C     ------------------------------
C     LUMPED OR DIAGONAL MASS MATRIX
C     ------------------------------
600   FACM = VOL/NNO
      DO 650  IEF=1,NEF
650   CM(IEF) = FACM
      RETURN
C
700   SUM = 0.
      DO 710  IEF=1,NEF,3
710   SUM = SUM + CM(IEF)
      FACM = VOL/SUM
      DO 750  IEF=1,NEF,3
      CM(IEF) = CM(IEF)*FACM
      CM(IEF+1) = CM(IEF)
750   CM(IEF+2) = CM(IEF)
      RETURN
C     ----------------------
C     CONSISTENT MASS MATRIX
C     ----------------------
800   KM1 = 1
      DO 890  IEF=1,NEF,3
      KM2 = KM1 + NEF-IEF+1
      KM3 = KM2 + NEF-IEF
      DO 850  JEF=IEF,NEF,3
      CM(KM2) = CM(KM1)
      CM(KM3) = CM(KM1)
      KM1 = KM1+3
      KM2 = KM2+3
850   KM3 = KM3+3
890   KM1 = KM1 + 2*(NEF-IEF) - 1
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MISE3D_S (SIGP,EPSP,IPEL,EPS,SIG,DP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     MAIN ROUTINE FOR ELASTO-PLASTIC MATERIAL BEHAVIOUR OF
C     THE 3-D SOLID ELEMENTS. THE THEORY IS BASED ON THE VON MISES
C     YIELD CONDITION AND USES NORMALITY RULE AND PRANDTL-REUSS FLOW
C     RULE TO CONVERT ELASTO-PLASTIC STRAIN INCREMENTS INTO STRESS INC
C	----------------------------------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     SIGP(6)    = STRESSES AT THE END OF THE PREVIOUS UPDATE
C     EPSP(6)    = STRAINS AT THE END OF THE PREVIOUS UPDATE
C     YIELD      = YIELD STRESS AT THE END OF THE PREVIOUS UPDATE
C     IPEL       = PLASTICITY FLAG (1=ELASTIC,2=PLASTIC)
C     EPS(6)     = CURRENT TOTAL STRAINS
C     SIG(6)     = RETURNED TOTAL STRESSES
C     DP(6,6)    = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DELEPS (6) = STRAIN INCREMENTS
C     DELSIG(6)  = STRESS INCREMENTS
C     DEPS(6)    = STRAIN SUBINCREMENTS (EQUIVALENCED WITH DELEPS)
C     DEPSE(6)   = ELASTIC PART OF STRAIN SUBINCREMENTS
C     DF(6)      = TOTAL DEVIATORIC STRESSES
C     DS(6)      = INCREMENTAL DEVIATORIC STRESSES
C     SX,SY...   = EQUIVALENCED WITH DF(6)
C     DX,DY...   = EQUIVALENCED WITH DS(6)
C     RATIO      = PART OF STRESS TAKEN ELASTICLY
C     PLAMDA     = PLASTIC STRAIN RATE MULTIPLIER
C     FT         = VON MISES YIELD FUNCTION
C     NINC,IINC  = NUMBER (COUNTER) OF SUBINCREMENTS
C	---------------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/,/HOOK/
C	---------------------------------------
C     NSINC      = VARIABLE CONTROLLING NUMBER OF SUBINCREMENTS
C     ITOLEY     = TOLERANCE ON EXCESS YIELD FUNCTION
C     YLD        = CURRENT UPDATED YIELD STRESS
C     ----------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION SIGP(6),EPSP(6),EPS(6),SIG(6),DP(6,6)
      DIMENSION DELEPS(6),DELSIG(6),DEPS(6),DEPSE(6),DF(6),DS(6)
C
      EQUIVALENCE (DELEPS,DEPS)
      EQUIVALENCE (SX,DF(1)),(SY,DF(2)),(SZ,DF(3)),(SXY,DF(4))
      EQUIVALENCE (SXZ,DF(5)),(SYZ,DF(6))
      EQUIVALENCE (DX,DS(1)),(DY,DS(2)),(DZ,DS(3)),(DXY,DS(4))
      EQUIVALENCE (DXZ,DS(5)),(DYZ,DS(6))
C     ------------------------------------------------------------
C     1. CALCULATE INCREMENTAL STRAINS
C     2. CALCULATE INCREMENTAL STRESSES ASSUMING ELASTIC BEHAVIOUR
C     3. CALCULATE TOTAL STRESSES ASSUMING ELASTIC BEHAVIOUR
C     ------------------------------------------------------------
      YLD2 = YLD*YLD
      DO 110  I=1,6
110   DELEPS(I) = EPS(I) - EPSP(I)
      CALL SOLSIG_S (DELEPS,DELSIG)
      DO 130  I=1,6
130   SIG(I) = SIGP(I) + DELSIG(I)
C     --------------------------------------------------
C     4. CALCULATE DEVIATORIC STRESSES AND CHECK WHETHER
C        STATE OF STRESS FALLS OUTSIDE THE YIELD SURFACE
C     --------------------------------------------------
      CALL DEVI3D_S (SIG,DF)
      F = 1.5*(SX*SX + SY*SY + SZ*SZ)+3.*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
      CONTINUE

c	IF(MEL.EQ.1) WRITE(*,*) F,YLD2,NSINC
	
	IF (F-YLD2)  410,410,450
C     --------------------------------------------------------
C     STATE OF STRESS WITHIN YIELD SURFACE - ELASTIC BEHAVIOUR
C     --------------------------------------------------------
410   IPEL = 1
      GOTO 900
C     ---------------------------------------------------------
C     STATE OF STRESS OUTSIDE YIELD SURFACE - PLASTIC BEHAVIOUR
C     ---------------------------------------------------------
450   IPEL = 2
      DO 460  I=1,6
460   SIG(I) = SIGP(I)
C     -------------------------------------------------------
C     5. DETERMINE PART OF STRAIN TAKEN ELASTICLY AND ADD
C        STRESSES DUE TO THE ELASTIC STRAINS TO PREVIOUS STRESSES
C     -----------------------------------------------------------
      CALL DEVI3D_S (SIG,DF)
      CALL DEVI3D_S (DELSIG,DS)
      A = DX*DX + DY*DY + DZ*DZ + 2.*(DXY*DXY + DXZ*DXZ + DYZ*DYZ)
      B = SX*DX + SY*DY + SZ*DZ + 2.*(SXY*DXY + SXZ*DXZ + SYZ*DYZ)
      E = SX*SX + SY*SY + SZ*SZ + 2.*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
     1  - 2.*YLD2/3.

	FY1 = 1.5*(SX*SX + SY*SY + SZ*SZ)+3.*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)

c	RTIO = (SQRT(F)-YLD)/(SQRT(F)-SQRT(FY1))
c	RATIO = 1 - RTIO


      RATIO = (-B + SQRT(B*B-A*E))/A

      DO 480  I=1,6
480   SIG(I) = SIG(I) + RATIO*DELSIG(I)

c	WRITE(*,*) MEL,A
C     --------------------------------------------------
C     6. DETERMINE NUMBER AND MAGNITUDE OF SUBINCREMENTS
C     --------------------------------------------------
600   NINC = NSINC*(SQRT(F)-YLD)/YLD + 1

      IF (NINC.GT.2*NSINC)  NINC = 2*NSINC
      FACT = (1.-RATIO) / NINC
      DO 610  I=1,6
610   DEPS(I) = FACT*DELEPS(I)
C     ----------------------------------------------------
C     7. CALCULATE PLASTIC STRAIN RATE MULTIPLIER AND
C        CONVERT ELASTIC PART OF STRAIN SUBINC. TO STRESSES
C     -----------------------------------------------------
      DO 890  IINC=1,NINC
      CALL DEVI3D_S (SIG,DF)
      PLAMDA = (SX*DEPS(1) + SY*DEPS(2) + SZ*DEPS(3) + SXY*DEPS(4) +
     1          SXZ*DEPS(5) + SYZ*DEPS(6)) * 0.5/YLD2
      IF (PLAMDA.LT.0.)  PLAMDA = 0.
      DO 710  I=1,6
710   DEPSE(I) = DEPS(I) - 3.*PLAMDA*DF(I)
      CALL SOLSIG_S (DEPSE,DELSIG)
      DO 720  I=1,6
720   SIG(I) = SIG(I) + DELSIG(I)
C     ----------------------------------------------
C     8. STRESS CORRECTION BACK TO THE YIELD SURFACE
C     ----------------------------------------------
      CALL DEVI3D_S (SIG,DF)
      F = 1.5*(SX*SX + SY*SY + SZ*SZ)+3.*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
C	ABCD=1./(DFLOAT(ITOLEY))
	ABCD=1./1000000.
      IF (F-YLD2.LE.ABCD)  GOTO 890

c      COEF = SQRT(YLD2)/SQRT(F)
c      DO 810  I=1,6
c810   SIG(I) = COEF*SIG(I)

      COEF = -1. + SQRT(YLD2/F)
      DO 810  I=1,6
810   SIG(I) = SIG(I) + COEF*DF(I)

C
890   CONTINUE
C     ----------------------------------------------
C     9. UPDATE PREVIOUS STRESSES, STRAINS AND YIELD
C     ----------------------------------------------
900   DO 910  I=1,6
      SIGP(I) = SIG(I)
910   EPSP(I) = EPS(I)
C     ------------------------------------
C     FORM THE ELASTO-PLASTIC MATERIAL LAW
C     ------------------------------------
      CALL DELP3D_S (SIG,DEPS,DP,DP,IPEL)  ! CHECK IT (SONGSAK)

	
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DEVI3D_S (SIG,S)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CALCULATES DEVIATORIC STRESSES  S(IJ) = SIG(IJ)-1/3SIG(KK)D(IJ)
C	---------------------------------------------------------------
C     SIG(6)    = STRESSES
C     S(6)      = DEVIATORIC STRESSES
C     ---------------------------------------------------------------
      DIMENSION SIG(6),S(6)
C
      SM = (SIG(1) + SIG(2) + SIG(3)) / 3.
      DO 100  I=1,3
      S(I)   = SIG(I) - SM
100   S(I+3) = SIG(I+3)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DELP3D_S (SIG,DEPS,DP,CP,IPEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     FORMS THE ELASTO-PLASTIC STRESS-STRAIN MATRIX FOR SOLID ELEMENTS
C	----------------------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     SIG(6)     = CURRENT TOTAL STRESSES
C     DEPS(6)    = SUBINCREMENTS IN STRAIN
C     DP(6,6)    = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     CP(36)     = DP(6,6)
C     IPEL       = PLASTICITY FLAG (1=ELASTIC,2=PLASTIC)
C	---------------
C     LOCAL VARIABLES
C	---------------
C     S(6)       = DEVIATORIC STRESSES
C     WP         = INCREMENTAL ELASTO-PLASTIC WORK
C
C     FOR VARIABLES IN COMMON BLOCK /HOOK/ SEE ROUTINE HOKLAW
C     --------------------------------------------------------------
      COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION SIG(6),DEPS(6),DP(6,6),CP(36),S(6)
C     --------------------
C     ELASTIC CONTRIBUTION
C     --------------------
      CALL CLEARA (CP,36)
      CP(1)  = A1
      CP(8)  = A1
      CP(15) = A1
      CP(2)  = B1
      CP(3)  = B1
      CP(9)  = B1
      CP(22) = C1
      CP(29) = C1
      CP(36) = C1
      IF (IPEL.EQ.1)  GOTO 300
C     -----------------------------------------------------
C     CALCULATE DEVIATORIC STRESSES AND CHECK FOR UNLOADING
C     -----------------------------------------------------
      CALL DEVI3D_S (SIG,S)
      WP = 0.
      DO 100  I=1,6
100   WP = WP + S(I)*DEPS(I)
      IF (WP.LT.0.)  GOTO 300
C     ---------------------------------------------
C     DEDUCT PLASTIC CONTRIBUTION IN LOWER TRIANGLE
C     ---------------------------------------------
      C2BETA = C2*1.5/YLD/YLD
      DO 200  I=1,6
      DO 200  J=1,I
200   DP(I,J) = DP(I,J) - C2BETA*S(J)*S(I)
C     --------------------------------------------
C     FILL IN SYMMETRIC ELEMENTS IN UPPER TRIANGLE
C     --------------------------------------------
300   DO 390  I=1,6
      DO 390  J=I,6
390   DP(I,J) = DP(J,I)
C
      RETURN
      END
C
C=====================================================================
C=====================================================================
      SUBROUTINE SOLSIG_S (STRAIN,STRESS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CONVERTS GLOBAL STRAINS TO GLOBAL STRESSES FOR THE 3-D SOLID
C	------------------------------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     STRAIN(6)  = GLOBAL STRAINS  (EPSX,EPSY,EPSZ,EPSXY,EPSXZ,EPSYZ)
C     STRESS(6)  = GLOBAL STRESSES (SIGX,SIGY,SIGZ,SIGXY,SIGXZ,SIGYZ)
C     FOR VARIABLES IN /HOOK/ REFERE TO ROUTINE HOKLAW
C     ----------------------------------------------------------------
      COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION STRAIN(6),STRESS(6)
C
      STRESS(1) = A1*STRAIN(1) + B1*(STRAIN(2)+STRAIN(3))
      STRESS(2) = A1*STRAIN(2) + B1*(STRAIN(1)+STRAIN(3))
      STRESS(3) = A1*STRAIN(3) + B1*(STRAIN(1)+STRAIN(2))
      DO 100  I=4,6
100   STRESS(I) = C1*STRAIN(I)
C
      RETURN
      END
C
C=====================================================================
C	==============================================================
      SUBROUTINE MATRIXES_S(A,B,S,MM)
C     ---------
C     TO COMPUTE TOTAL STIFFNESS MATRIX
C     K=K(LINEAR)-K(EAS)
C     ---------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(MM,24),B(MM,MM),AB(24,MM),S(1)
	DO 10 I=1,24
	DO 10 J=1,MM
	AB(I,J)=0.0
10    CONTINUE
      DO 30 I=1,24
	DO 30 J=1,MM
	DO 30 K=1,MM
	AB(I,J)=AB(I,J)+A(K,I)*B(K,J)
  30  CONTINUE
      KK=1
	DO 40 I=1,24
	DO 40 J=I,24
	DO 50 K=1,MM
	S(KK)=S(KK)-AB(I,K)*A(K,J)
50    CONTINUE
      KK=KK+1
40    CONTINUE
      RETURN
      END


C	==============================================================






C=====================================================================
      SUBROUTINE REARGE_S (SK,SK1)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     EVALUATES THE JACOBIAN MATRIX AT POINT (R,S,T),ITS DETERMINANT
C     AND INVERSE FOR A CURVILINEAR ISOPARAMETRIC HEXAHEDRON (8-21 NO)
C	----------------------------------------------------------------
C     INPUT OUTPUT VARIABLES
C	----------------------
C     XY(3,NNO)    = COORDINATES OF ELEMENT NODES
C     P(3,NNO)     = SHAPE FUNCTION DERIVATIVES WITH RESPECT TO R,S,T
C     XJI(3,3)     = INVERSE OF THE JACOBIAN MATRIX
C     DET          = DETERMINANT OF THE JACOBIAN MATRIX
C     MEL          = ELEMENT NUMBER
C     NNO          = NUMBER OF NODES USED TO DESCRIBE THIS ELEMENT
C     ----------------------------------------------------------------
      DIMENSION SK(24,24),SK1(24,24),INSONG(8)
C

	INSONG(1) = 6
	INSONG(2) = 2
	INSONG(3) = 1
	INSONG(4) = 5
	INSONG(5) = 7
	INSONG(6) = 3
	INSONG(7) = 4
	INSONG(8) = 8

      DO 100  I=1,8
	ING = INSONG(I)
      DO 100  J=1,8
	JNG = INSONG(J)
	SK1(3*I-2,3*J-2) = SK(3*ING-2,3*JNG-2)
	SK1(3*I-1,3*J-1) = SK(3*ING-1,3*JNG-1)
100	SK1(3*I-0,3*J-0) = SK(3*ING-0,3*JNG-0)


C
      RETURN
      END
C
C=====================================================================



C=====================================================================

      SUBROUTINE SOLIDEAS(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,MWG
     +				,ALPHA,SEL,SEDI,FIN,HINFC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     --------------------------------------------------------------
C     MAIN PROGRAM FOR THE 3-D SOLID
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
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DP(64)        = ELASTIC OR ELASTO-PLASTIC MATERIAL MATRIX
C     H(21)         = INTERPOLATION FUNCTIONS
C     HD(3,21)      = SHAPE FUNCTION DERIVATIVES WITH RESP.TO R,S,T
C     XJ(3,3)       = JACOBIAN MATRIX
C     XJI(9)        = INVERSE OF THE JACOBIAN MATRIX
C     B(3*NNO)      = COMPRESSED STRAIN-DISPLACEMENT MATRIX
C     DISD(9)       = DISPLACEMENT DERIVATIVES
C     EPS(6)        = GAUSS POINT STRAINS
C     EPSQ(6)       = QUADRATIC PART OF GAUSS POINT STRAINS
C     SIG(6)        = GAUSS POINT STRESSES
C     IPEL          = SECTION PLASTICITY INDICATER (1=ELASTIC,2=PL)
C     RN,SN,TN      = NON-DIMENSIONAL COORDINATES
C     DET           = DETERMINANT OF THE INVERSE JACOBIAN
C     DVOL          = INTEGRATION CONSTANT
C     ----------------
C     VARIABLES OF EAS METHOD
C     ----------------
C     MM             =EAS TERMS OF ALPHA(9, 15,21 OR 30)
C     XJO(3,3)       =JACOBIAN MATRIX AT THE ORIGIN (R,S,T=0.0)
C     XJ(3,3)        =JACOBIAN MATRIX
C     TTO(6,6)       =COEFFICIENT MATRIX OF ENHANCED STRAIN
C     TM(6,MM)       =COEFFICIENT MATRIX OF ENHANCED STRAIN
C     SED(MM,MM)     =MATRIX OF EAS, TRANSPOSE(M)*E*M 
C     SEL(MM,24)     =MATRIX OF EAS, TRANSPOSE(M)*E*B
C     SEDI(MM,MM)    =INVERSE OF MATRIX SED
C     EAS(6,1)       =ENHANCED STRAIN VECTOR
C     ALPHA(MM,1)    =ENHANCED TERMS ALPHA
C     TMT(MM,6)      =TRANSPOSE OF MATRIX TM(6,MM)
C     RE1(24)        =EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     RE2(24,1)      =EQUILIBRIUM FORCE OF EAS METHOD
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
C

C	...SELECT ENAHNACED INTERPOLATION STRAIN

C     PARAMETER (MM= 9)
C     PARAMETER (MM=12) ! FOR OTHER APPLICATION
C     PARAMETER (MM=15) 
C     PARAMETER (MM=18)
C      PARAMETER (MM=21) ! WANG
C      PARAMETER (MM=24) ! SACHARUCK
	COMMON /MMENH/ MM,MM1,MM2,NDIMC
c     PARAMETER (MM=30) 
C
C	STREAM INPUT AND OUTPUT FOR WORKING ARRAY
C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(*),COORD(*)
      DIMENSION EDIS(*),EDISI(*),RE(*),SG(300)
 
C	BASIC MATRIX FOR COMPUTATION

      DIMENSION DP(6,6),H(21),P(3,21),XJI(3,3),B(63),DISD(9)
	DIMENSION BM(6,24),DE(6,6)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6)
      DIMENSION STRESSLO(6) !------- LOCAL STRESS BY BJ

	DIMENSION SK(24,24),DEP(6,6)

C     --------------
C       EAS METHOD
C     --------------

	DIMENSION GE(6,MM)
	DIMENSION T0(6,6),TT(6,6),TT0(6,6)   
	DIMENSION XJ(3,3),XJO(3,3),TM(6,MM)
	DIMENSION EAS(6),ALPHA(MM),RE1(24),RE2(24),HINFC(MM)
      DIMENSION SED(MM,MM),SEL(MM,24),SEDI(MM,MM),RH(MM)

C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)

	DIMENSION FIN(NEF)

	DIMENSION DFVP(24)


C	WRITE(*,*) DFVP
C		WRITE(*,*) 'KAK',MEL

	DO I = 1,MM
	HINFC(I) = 0.0
	ENDDO

C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------

      CALL HOKLAW_S (PROPM,PROPG,1)

C     --------------------------------------------------------------
C     COMPUTE JACOBIAN MATRIX AT THE NATURAL CENTER (R=0.,S=0.,T=0.)
C     --------------------------------------------------------------

      CALL SHAP3D_S (0.0D0,0.0D0,0.0D0,H,P,NODEX,NNO)
	CALL JACO3D_S (COORD,P,XJO,XJI,DETO,MEL,NNO)

C     ------------------------------------------------------------
C	FIND THE LOCAL VECTORS,LOCAL COORDINATES AND DISPLACEMENTS
C     ------------------------------------------------------------      
  	!CALL SOLIRST(COORD,P,VR,VS,VT)   ! LOCAL VECTORS	-----BY BJ
      
C     ----------------
C     COMPUTE T MATRIX (FOR EAS METHOD)
C     ---------------

	CALL MATRIXT_S (XJO,T0)
      CALL INVMATRIX (T0,TT,6)

	TT0 = TRANSPOSE(TT)

      MGR = NGR
      MGS = NGS
      MGT = NGT
      IF (ITASK.NE.5) GOTO 10
      MGR = 3
      MGS = 3
      MGT = 3

C     -----------------------------
C     SETTING INDEX FOR GAUSS POINT
C     -----------------------------

 10   IPT = 0

C     ----------------------------------------
C     INITIALIZATION OF SOME VARIABLE MATRICES
C     ----------------------------------------
	
	SG   = 0.0
	SK   = 0.0
	SEL  = 0.0
	SED  = 0.0
	RH   = 0.0
	RE1  = 0.0
	RE2  = 0.0
	DFVP = 0.0

C	WRITE(*,*) MEL

C	------------------------
C	 LOOP OVER GAUSS POINT
C	------------------------

      DO 900  IGR=1,MGR
        RI = GLOC(IGR,MGR)
         DO 900  IGS=1,MGS
          SI = GLOC(IGS,MGS)
           DO 900  IGT=1,MGT
            TI = GLOC(IGT,MGT)
            
		  WT = GWT(IGR,MGR)*GWT(IGS,MGS)*GWT(IGT,MGT)
            IPT = IPT+1

C     ---------------------------------------------------
C     SHAPE FUNCTIONS (H), DERIVATIVES (P),
C	INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     ---------------------------------------------------

      CALL SHAP3D_S (RI,SI,TI,H,P,NODEX,NNO)
      CALL JACO3D_S (COORD,P,XJ,XJI,DET,MEL,NNO)
      DVOL = WT*DET

C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------

      IF (ITASK.NE.5)  GOTO 50
      CALL SOMASS_S (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
      GOTO 900

C     ----------------------------------------
C     COMPACTED STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------

50    CONTINUE

	CALL BMATSLD_S (P,XJI,B,BM,NNO)
	CALL DMATSLD_S (DE)

C     -----------------------------------------     
C     FOR ENHANCED STRAIN METHOD 
C     COMPUTE MATRIX [GE]
C     -----------------------------------------

	CALL SOLIDGE_S (RI,SI,TI,DETO,DET,TT0,MM,GE)

C     -------------------------------------
C	COMPUTE ENHANCED STIFFNESS MATRIX SED
C	COMPUTE ENHANCED COUPLING  MATRIX SEL
C     -------------------------------------

	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(DE,GE))*DVOL
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(DE,BM))*DVOL

      IF (NLOPT+ITASK.EQ.1)  GOTO 700

C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------

      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,6
100   TAU(I) = WA(I,IPT)*DVOL
      GOTO 800

C     -------------------------------------
C      COMPUTE DISPLACEMENT GRADIENT (DISD)
C     -------------------------------------

 200  CONTINUE

	DISD = 0.

      DO IEF=1,NEF,3

      JEF = IEF+1
      KEF = IEF+2

      DISD(1) = DISD(1) + B(IEF)*EDIS(IEF)   ! H,x *u
      DISD(2) = DISD(2) + B(JEF)*EDIS(JEF)   ! H,y *v
      DISD(3) = DISD(3) + B(KEF)*EDIS(KEF)   ! H,z *w
      DISD(4) = DISD(4) + B(JEF)*EDIS(IEF)   ! H,y *u
      DISD(5) = DISD(5) + B(KEF)*EDIS(IEF)   ! H,z *u
      DISD(6) = DISD(6) + B(IEF)*EDIS(JEF)   ! H,x *v 
      DISD(7) = DISD(7) + B(KEF)*EDIS(JEF)   ! H,z *v
      DISD(8) = DISD(8) + B(IEF)*EDIS(KEF)   ! H,x *w
      DISD(9) = DISD(9) + B(JEF)*EDIS(KEF)   ! H,y *w

	END DO

C     ------------------------------
C     LINEAR COMPATIBLE STRAIN TERMS
C     ------------------------------

      STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)
      STRAIN(4) = DISD(4) + DISD(6)
      STRAIN(5) = DISD(5) + DISD(8)
      STRAIN(6) = DISD(7) + DISD(9)

C     -------------------------------
C     FOR EAS TERM  EAS=[GE]*{ALPHA}
C     -------------------------------

	EAS = MATMUL(GE,ALPHA)

C	-------------------------------
C	COMPUTE TOTAL COMPATIBLE STRAIN
C	------------------------------- 

	DO I=1,6
      STRAIN(I)=STRAIN(I) - EAS(I)
	END DO

C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT QUADRATIC STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400

      QSTRAI(1) = .5*(DISD(1)*DISD(1)+DISD(6)*DISD(6)+DISD(8)*DISD(8))
      QSTRAI(2) = .5*(DISD(4)*DISD(4)+DISD(2)*DISD(2)+DISD(9)*DISD(9))
      QSTRAI(3) = .5*(DISD(5)*DISD(5)+DISD(7)*DISD(7)+DISD(3)*DISD(3))
      QSTRAI(4) =     DISD(1)*DISD(4)+DISD(6)*DISD(2)+DISD(8)*DISD(9)
      QSTRAI(5) =     DISD(1)*DISD(5)+DISD(6)*DISD(7)+DISD(8)*DISD(3)
      QSTRAI(6) =     DISD(4)*DISD(5)+DISD(2)*DISD(7)+DISD(9)*DISD(3)

C	----------------------
C	SUBTRACT ALMASI STRAIN
C	----------------------

      DO 390  I=1,6
 390  STRAIN(I) = STRAIN(I) - QSTRAI(I)

C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------

400   GOTO (410,420,430,440,450,460,470,470),MTMOD

410   CONTINUE

C	----------------------------------------------------------
C	LINEAR ELASTIC HOOKE'S CONSTITUTUVE STRESS-STRAIN RELATION
C	----------------------------------------------------------

	CALL SOLSIG_S (STRAIN,STRESS)
      
! ----TRANSFORM GLOBAL STRESS TO LOCAL STRESS -------- BY BJ
      !CALL STRESSLOCALTRANS(VR,VS,VT,STRESS,STRESSLO) 
      
      DO 415  I=1,6
415   WA(I,IPT) = STRESS(I)
      GOTO 500

420	CONTINUE
	
	WRITE(*,*) 'MATERIAL NUMBER IS NOT AVAILIBLE NOW'
	STOP

	GOTO 500


430	CONTINUE

	CALL DMMISE (WA(1,IPT),WA(7,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DEP)

      GOTO 500


440   CONTINUE

C	CALL MOHRCRI3D(WA(1,IPT),WA(7,IPT),STRAIN,STRESS,DEP,
C	1			   PROPM,WA(13,IPT))
	CALL MOHR3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			 STRAIN,DEP,PROPM)

	GO TO 500

450   CONTINUE

	CALL DRUG3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			 STRAIN,DEP,PROPM,6,0)

      GOTO 500


460   CONTINUE

	CALL DAMAGE(STRAIN,STRESS,DEP,PROPM,WA(1,IPT))

      GOTO 500


	


470	CONTINUE !IF(WA(14,IPT).NE.0.0) WRITE(*,*) MEL

	CALL PLVIS(WA(1,IPT),WA(7,IPT),WA(14,IPT),
	1		   DFVP,STRESS,STRAIN,PROPM,
     2		   BM,DVOL,WA(13,IPT),DEP)

	GO TO 500
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------

 500  CONTINUE

	DO I=1,6
      TAU(I) = STRESS(I)*DVOL
	END DO

	IF (ITASK.LE.3) GOTO 520
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900
 520	CONTINUE  

C     ---------------------------------------
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     ---------------------------------------

	RE1 = RE1 + MATMUL(TRANSPOSE(BM),STRESS)*DVOL

C     ----------------------------------------
C     COMPUTE EQUILIBRIUM FORCE FOR EAS METHOD
C     ----------------------------------------

	RH = RH + MATMUL(TRANSPOSE(GE),STRESS)*DVOL

	IF (IFEIG.EQ.0.AND.ISOLOP.EQ.4) GOTO 800
C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------

 151  IF (IFREF) 900,700,900
 700  IF (NLOPT.NE.0) GOTO 750

C	****************
C	LINEAR ANALYSIS
C	****************

C	CALL SOK0NL (S,B,NEF,DVOL)

	DEP = 0.
	SK  = 0.

	CALL DMATSLD_S (DEP)

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,24
	      DO JSK =ISK,24
	         KSK = KSK+1 
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	GOTO 790

C	*************************
C	...GEOMETRIC NONLINEARITY
C	*************************

 750  IF (MTMOD.LE.2)THEN

C	CALL SOK0NL (S,B,NEF,DVOL)   ! LINEAR  COMPATIBLE STIFF.

	DEP =0.
	SK = 0.

	CALL DMATSLD_S (DEP)

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,24
	      DO JSK =ISK,24
	         KSK = KSK+1 
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	END IF	

C	************************
C	...MATERIAL NONLINEARITY
C	************************
	IF (MTMOD.GT.2)THEN

	IF(NLOPT+ITASK.EQ.1) THEN


	IF (MTMOD.EQ.3) THEN

	DEP =0.

	CALL DMMISE (WA(1,IPT),WA(7,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DEP)

c	CALL MISE3D_S (WA(1,IPT),WA(7,IPT),WA(13,IPT),
c     1             STRAIN,STRESS,DEP)

	ENDIF


	IF (MTMOD.EQ.4) THEN

	DEP =0.0

	CALL MOHRCRI3D(WA(1,IPT),WA(7,IPT),STRAIN,STRESS,DEP,
	1			   PROPM,WA(13,IPT))

	ENDIF

	IF (MTMOD.EQ.5) THEN

	DEP =0.

	CALL DRUG3D(WA(1,IPT),WA(7,IPT),WA(13,IPT),STRESS,
	1			STRAIN,DEP,PROPM,6,0)

	ENDIF


	IF (MTMOD.EQ.6) THEN

	DEP =0.

	CALL DAMAGE(STRAIN,STRESS,DEP,PROPM,WA(1,IPT))

	ENDIF


	IF (MTMOD.EQ.7) THEN

	DEP =0.

	CALL PLVIS(WA(1,IPT),WA(7,IPT),WA(14,IPT),
	1		   DFVP,STRESS,STRAIN,PROPM,
     2		   BM,DVOL,WA(13,IPT),DEP)

	ENDIF




	ENDIF ! <--- ENDIF NLOPT+ITASK.EQ.1


C	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(DEP,GE))*DVOL
C	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(DEP,BM))*DVOL


	SK = 0.

	  SK = (MATMUL(TRANSPOSE(BM),MATMUL(DEP,BM)))*DVOL
	  KSK = 0
	    DO ISK =1,24
	      DO JSK =ISK,24
	         KSK = KSK+1 
	       S(KSK)=S(KSK)+SK(ISK,JSK)
	      END DO
	    END DO

	END IF


C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------

790   IF (NLOPT.LE.1) GOTO 810
800	CONTINUE

C	------------------------------------------------
C	ADDITIONAL COMPUTE NONLINEAR GEOMATRIC STIFFNESS
C	------------------------------------------------
	CALL KNLSTIFF (SG,TAU,B)


 810  CONTINUE !TIM(12) = TIM(12) + TIM2-TIM1


 900  CONTINUE

	IF(ITASK.EQ.5) GOTO 150   !MASS MATRIX
	IF(IFEIG.EQ.0.AND.ISOLOP.EQ.4) GOTO 950

	CALL INVMATRIX (SED,SEDI,MM)

C     ------------------------------------
C     FOR ENHANCED STRAIN METHOD
C     K=TRANSPOSE(L)*INVERSE(D)*L
C     K(FINAL)=K.Linear(PURE DISP.)-K(EAS)
C     ------------------------------------
	CALL MATRIXES_S (SEL,SEDI,S,MM)


950	KGG = 0
	DO I = 1,24
	DO J = I,24
	KGG = KGG + 1
	S(KGG) = S(KGG) + SG(KGG)
	ENDDO
	ENDDO
	 

C     ------------------------------------
C     FOR NONLINEAR ANALYSIS OF EAS METHOD
C     ------------------------------------
	DO IHH = 1,MM
	HINFC(IHH) = RH(IHH)
	ENDDO

C	INTERNAL FORCE FROM EAS
	RE2 = MATMUL(MATMUL(TRANSPOSE(SEL),SEDI),RH)



	IF (ITASK+NLOPT.EQ.1) GOTO 150

	DO I=1,24
	 RE(I) = RE1(I) - RE2(I) + DFVP(I)
	END DO

C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)

150	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF




      RETURN
      END

C	=============================================================
C	=============================================================



C=====================================================================
C	***************************************************************

	SUBROUTINE DMATSLD_S(D)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
 
C	**********************************
C	COMPUTE LINEAR CONSTITUTIVE MATRIX
C	**********************************
   
      COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST

	DIMENSION D(6,6)

      
	DO I =1,6
	 DO J =1,6
	  D(I,J) = 0.0
	 END DO
	END DO

C	LINEAR CONSITITUTIVE 

      D(1,1)  = A1
      D(1,2)  = B1
      D(1,3)  = B1
      D(2,1)  = B1
      D(2,2)  = A1
      D(2,3)  = B1
      D(3,1)  = B1
      D(3,2)  = B1
      D(3,3)  = A1

      D(4,4) = C1
      D(5,5) = C1
      D(6,6) = C1

      RETURN
      END

C	**************************************************

	SUBROUTINE SOLIDGE_S(R,S,T,DETJO,DETJ,TT,MM,GE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	-----------------------------------------------
C	SUBROUTINE COMPUTE THE COMPATIBLE ENHANCED TERM
C	-----------------------------------------------

	DIMENSION GM(6,MM),GE(6,MM),TT(6,6),EM(6,MM),GM1(6,18)

	

	DO I =1,6
	 DO J=1,MM
		GM(I,J) = 0.
	    GE(I,J) = 0.
		EM(I,J) = 0.
	 END DO
	END DO

	IF(MM.EQ.12) THEN
	GM(1,1) =  -2.0*R
	GM(2,2) =  -2.0*S
	GM(3,3) =  -2.0*T

	GM(4,4) = -3.0*R*S 
	GM(4,5) = -3.0*S*T
	GM(4,6) = -3.0*R*T 
	
	GM(5,7) = -3.0*R*S 
	GM(5,8) = -3.0*S*T
	GM(5,9) = -3.0*R*T

	GM(6,10) = -3.0*R*S 
	GM(6,11) = -3.0*S*T
	GM(6,12) = -3.0*R*T

	GO TO 100

	ENDIF

	GM(1,1) =  -2.0*R
	GM(2,2) =  -2.0*S
	GM(3,3) =  -2.0*T
	GM(4,4) =  -2.0*R
	GM(4,5) =  -2.0*S

	IF(MM.EQ.5)GOTO 100

	GM(5,6) =  -2.0*R
	GM(5,7) =  -2.0*T	
	GM(6,8) =  -2.0*S
	GM(6,9) =  -2.0*T
	
	IF(MM.EQ.9)GOTO 100

	GM(4,10) = -3.0*R*S 
	GM(4,11) = -3.0*S*T
	GM(4,12) = -3.0*R*T 
	
	GM(5,13) = -3.0*R*S 
	GM(5,14) = -3.0*S*T
	GM(5,15) = -3.0*R*T

	GM(6,16) = -3.0*R*S 
	GM(6,17) = -3.0*S*T
	GM(6,18) = -3.0*R*T

	GM(1,19) = -3.0*R*T 
	GM(1,20) = -3.0*R*S
	 
	GM(2,21) = -3.0*R*S 
	GM(2,22) = -3.0*S*T 

	GM(3,23) = -3.0*S*T 
	GM(3,24) = -3.0*R*T

	IF(MM.EQ.24)GOTO 100

100	CONTINUE


	IF(MM.EQ.18) THEN

	DO I =1,6
	DO J=1,18
	GM1(I,J) = 0.0
	END DO
	END DO

c	write(*,*) MM
c	stop

	PM1 = -2.0*R*(1.0-S*S)*(1.0-T*T)
	PM2 = -2.0*S*(1.0-R*R)*(1.0-T*T)
	PM3 = -2.0*T*(1.0-R*R)*(1.0-S*S)


	GM1(1,1) = PM1
	GM1(2,2) = PM2
	GM1(3,3) = PM3
	GM1(4,4) = PM2
	GM1(4,5) = PM1
	GM1(5,6) = PM3
	GM1(5,7) = PM1
	GM1(6,8) = PM3
	GM1(6,9) = PM2

	GM1(1,10) = PM2*PM3
	GM1(2,11) = PM1*PM3
	GM1(3,12) = PM1*PM2
	GM1(4,13) = PM3*PM2
	GM1(4,14) = PM3*PM1
	GM1(5,15) = PM2*PM3
	GM1(5,16) = PM2*PM1
	GM1(6,17) = PM1*PM3
	GM1(6,18) = PM1*PM2


	

	DO I = 1,6
	DO J = 1,18

	GM(I,J) = GM1(I,J)
	ENDDO
	ENDDO

	ENDIF

 

	EM = MATMUL(TT,GM)
      CONST = (DETJO/DETJ)

C	COMPUTE [GE] MATRIX

	DO IEM = 1,6
		DO  JEM = 1,MM 
		GE(IEM,JEM) = CONST*EM(IEM,JEM)
		END DO
	END DO

	RETURN 
	END

C	*****************************************************************
C=====================================================================

      SUBROUTINE BMATSLD_S (P,XJI,B,BM,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     EVALUATES THE GLOBAL LINEAR STRAIN-DISPLACEMENT MATRIX FOR A
C     CURVILINEAR ISOPARAMETRIC HEXAHEDRON (8 TO 21 NODES)
C	----------------------------------------------------
C     INPUT OUTPUT VARIABLES
C	----------------------
C     P(3,NNO)    = SHAPE FUNCTION DERIVATIVES WITH RESPECT TO R,S,T
C     XJI(3,3)    = INVERSE OF THE JACOBIAN MATRIX
C     B(3*NNO)    = COMPRESSED LINEAR STRAIN DISPLACEMENT MATRIX
C     NNO         = NUMBER OF NODES USED TO DESCRIBE THIS ELEMENT
C     ----------------------------------------------------------------

      DIMENSION P(3,21),XJI(3,3),B(63),BM(6,24)
C
      CALL CLEARA (B,63)

C	------------------------------------------------
C	 COMPUTE COMPACTED STRESS-STRAIN GRADIENT VECTOR
C	------------------------------------------------

      K = 1
	KK= 0

      DO 100  INO=1,NNO
      DO 50   I=1,3
      B(K)   = B(K)   + XJI(1,I)*P(I,INO)
      B(K+1) = B(K+1) + XJI(2,I)*P(I,INO)
      B(K+2) = B(K+2) + XJI(3,I)*P(I,INO)
 50   CONTINUE
      K = K+3
	KK=KK+3
 100  CONTINUE
C
 
C	----------------------------------------
C	 COMPUTE COMPACTED STRESS-STRAIN MATRIX
C	----------------------------------------

	DO I=1,6
	 DO J =1,24
	  BM(I,J)=0.
	 END DO
	END DO

	DO INO=1,8

	HX =0.
	HY =0.
	HZ =0.

	DO I=1,3
        HX  = HX  + XJI(1,I)*P(I,INO)
        HY  = HY  + XJI(2,I)*P(I,INO)
        HZ  = HZ  + XJI(3,I)*P(I,INO)
	END DO

C	--------------------------
C  	 COMPUTE MATRIX [B](6,24)
C	--------------------------

	BM(1,3*INO-2) = HX
	BM(2,3*INO-1) = HY
	BM(3,3*INO-0) = HZ

	BM(4,3*INO-2) = HY
	BM(4,3*INO-1) = HX

	BM(5,3*INO-2) = HZ
	BM(5,3*INO-0) = HX

	BM(6,3*INO-1) = HZ
	BM(6,3*INO-0) = HY

	END DO

      RETURN
      END

C	***************************************************************
C	*****************************************************************

	SUBROUTINE INVMATRIX(A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------
C     MATRIX INVERSION BY ELIMMINATION WITH PARTIAL PIVOTING
C     AND DETERMINANT OF MATRIX A
C     ROGINAL MATRIX =A, INVERSE MATRIX=B
C
C     A: INPUT SQUARE MATRIX
C     B: INVERSE MATRIX OF A
C     N: SIZE OF A AND B
C     EPS: CONTROL VARIABLE
C     DET: DETERMINANT OF A
C     --------

      DIMENSION A(N,N),B(N,N),C(N,N)

C     CONSTRUCT IDENTITY MATRIX B(I,J)=I
      EPS = 1.0E-12
	DO 88 I=1,N
	DO 88 J=1,N
	C(I,J)=A(I,J)
  88  CONTINUE	
 	DO 6 I=1,N
	DO 5 J=1,N
	IF(I-J) 4,3,4
   3  B(I,J)=1.0
      GOTO 5 
   4  B(I,J)=0.0
   5  CONTINUE
   6  CONTINUE

C     LOCATE MAXIMUM MAGNITUDE A(I,K) ON OR BELOW MAIN DIAGONAL
      DET=1.0
	DO 45 K=1,N
	IF (K-N) 12,30,30
   12 IMAX=K
      AMAX=ABS(C(K,K))
	KP1=K+1
	DO 20 I=KP1,N
	IF (AMAX-ABS(C(I,K))) 15,20,20
   15 IMAX=I
      AMAX=ABS(C(I,K))
   20 CONTINUE

C     INTERCHANGE ROWS IMAX AND K IF IMAX NOT EQUAL TO K
      IF (IMAX-K) 25,30,25
   25 DO 29 J=1,N
      ATMP=C(IMAX,J)
	C(IMAX,J)=C(K,J)
	C(K,J)=ATMP
	BTMP=B(IMAX,J)
	B(IMAX,J)=B(K,J)
   29 B(K,J)=BTMP
      DET=-DET
   30 CONTINUE

C     TEST FOR SINGULAR MATRIX
C      IF (ABS(C(K,K))-EPS) 33,33,35
C   33 WRITE(*,*) 'SINGULAR MATRIX EPS- INSIDE SUB INVMATRIX'
C	STOP
   35 DET=C(K,K)*DET

C     DIVIDE PIVOT ROW BY ITS MAIN DIAGONAL ELEMENT
      DIV=C(K,K)
	
	IF(DIV.EQ.0.0) DIV = EPS
	IF(ABS(DIV).LT.EPS) DIV = EPS*DIV/ABS(DIV)

	DO 38 J=1,N
	C(K,J)=C(K,J)/DIV
   38 B(K,J)=B(K,J)/DIV

C     REPLACE EACH ROW BY LINEAR COMBINATION WITH PIVOT ROW
      DO 43 I=1,N
	AMULT=C(I,K)
	IF (I-K) 39,43,39
   39 DO 42 J=1,N
      C(I,J)=C(I,J)-AMULT*C(K,J)
   42 B(I,J)=B(I,J)-AMULT*B(K,J)
   43 CONTINUE
   45 CONTINUE
C
      RETURN
	END
C	=========================================================


	SUBROUTINE MOHRCRI3D(SIG,EPS,STRAIN,STRESS,DEP,PROPM,WK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	****************************************************

C	SUBROUTINE COMPUTE SOIL PLASTIC CONSTITUTIVE MATRIX
C		COMPUTE CORRECTION STRESS AND STRAIN

C	****************************************************	


      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT

	COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
c	COMMON /CRIP/  OCR,VO,DEN,VGL,SWL,SM,OK
c	COMMON /MOHR/  COH,ANG,BETA

C	...STRESS AND STRAIN STORE DATA

	DIMENSION SIG(6),EPS(6),STRAIN(6),STRESS(6)
	DIMENSION DELSIG(6),DELEPS(6)

	DIMENSION DD(6,6),DV(36),PROPM(*)

C	...CONSTITUTIVE MATRIX

	DIMENSION DES(6,6),DPS(6,6),DEP(6,6)
	DIMENSION PRESP(3) ! PRICIPLE STRESS OF PRES
	DIMENSION CURSP(3) ! PRICIPLE STRESS OF CURS

	DIMENSION AV1(6,1),AV2(6,1),AV3(6,1),AV(6,1)
	DIMENSION DA(6,1),SONGST(6),WK(2)

C	...COMPUTE SOIL PARAMETER

	IEL = MEL

	COH  = PROPM(10)
	ANG  = PROPM(11)
	BETA = PROPM(12)

	PI  = 3.141592654 
	PHI = ANG*PI/180.

C	WRITE(*,*)COH

c	WRITE(*,*)'RUNNING',IEL

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

C	CHANGE SIGN FOR POSITIVE COMPRESSIVE CRITERION DEFINITION

	CALL	CONVSIG3D(SIG,SIG)

C	**************************************************
C	  COMPUTE THE STRESS PREVIOUS YIELD SURFACE (PRES)
C	**************************************************

	P = (1.0/3.0)*(SIG(1)+SIG(2)+SIG(3))

C	WRITE(*,*)P

C	...FIRST ENTRY MODULE

	IF(P.EQ.0.0)THEN
	WK(2) = COH*COS(PHI)
	END IF


	IF(P.EQ.0.0) GOTO 100 !FIRST ENTRY PREVIOUS STRESS = 0.0

	CALL	SINVAR3D(SIG,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)  
	CALL	YIEDSF3D(VARI1,VARJ2,PHI,THETA,PRES)

C	WRITE(*,*)VARI1

C	************************************************
C		COMPUTE CURRENT YIELD SURFACE
C	************************************************

C	...READ PLASTIC STRAIN FROM PREVIOUS SESSION

	PYLD  = WK(2)
	PREY  = PYLD
	YIELD = COH*COS(PHI) + PYLD

	EPSTN = WK(1)

C	*******************************************
C	COMPUTE CURRENT STRESS YIELD SURFACE VECTOR
C	*******************************************

C	...UPDATING INCREMENTAL STRESSES

	DES = 0.0

	CALL DMATSLD_S(DES)

	DO I = 1,6
	DELEPS(I) = STRAIN(I) - EPS(I)
	END DO

	DELSIG = MATMUL(DES,DELEPS)

C	...CHANGE SIGN FOR POSITIVE CRITERION DEFINITION

	CALL	CONVSIG3D(DELSIG,DELSIG)


C	...COMPUTE CURRENT STRESS (TOTAL STRESS)

	DO I =1,6
	STRESS(I) = SIG(I) + DELSIG(I)
	END DO

	CALL	SINVAR3D(STRESS,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)  
	CALL	YIEDSF3D(VARI1,VARJ2,PHI,THETA,CURS)

C	WRITE(*,*)PREY,CURS

C	**********************************************
C	CHECK THE CURRENT STRESS OUTSIDE YIELD SURFACE
C	**********************************************

C	...PREVIOUS YIELD

	IF(PRES.GE.PREY)THEN
	   IF(CURS.LE.PREY)GOTO 100 ! ELASTIC CASE--UNLOADING	
	    RFACT = 0.0
	END IF

	RFACT = 0.

C	...CROSSING FROM YIELD SURFACE

	IF(PRES.LT.PREY)THEN
	   IF(CURS.GE.PREY)THEN 
C		CALL RFMOHR3D(RFACT,SIG,DELSIG,PHI,YIELD) ! R-FACTOR (EXACT CORRECT)
		CALL RFMOHR3D_NEW(RFACT,SIG,DELSIG,PHI,PREY) ! SONGSAK NOV2005
			IF (RFACT.EQ.1.0) THEN
				RFACT = (CURS - PREY) / (CURS - PRES)
				RFACT = 1.0 - RFACT
			ENDIF
	   ELSE
	    GOTO 100 ! ELASTIC -- STILL INSIDE THE YIELD SURFACE
	   END IF 	   
	END IF

C	WRITE(*,*)SIG
C	PAUSE

C	*********************************
C	REDUCTION STRESS TO YIELD SURFACE
C	*********************************

	MSTEP = 30*(CURS/PREY)+1
	ASTEP = 1.0*MSTEP
	REDUC = 1.0-RFACT

C	...COMPUTE INCREMENTAL STRESS CLOSED TO THE YIELD SURFACE

	DO I =1,6
	STRESS(I) = SIG(I)+RFACT*DELSIG(I)
	DELSIG(I) = REDUC*DELSIG(I)/ASTEP
	END DO

C	...COMPUTE SUB INCREMENTAL STRESSES

	DO ISTEP = 1,MSTEP

C	****************************************************
C	...COMPUTE CURRENT STRESS SURFACE AT YIELD CRITERION
C	****************************************************

	DES = 0.
	
C	...COMPUTE [DES]

	CALL DMATSLD_S(DES)

C	***************************
C	    COMPUTE FLOW VECTOR
C	***************************

	CALL	SINVAR3D(STRESS,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)  
	CALL	YIEDSF3D(VARI1,VARJ2,PHI,THETA,CURS)
	CALL	FLOWS3D(STRESS,AV,DA,DES,EPSTN,CONST,PHI)

C	********************************************************
C	...COMPUTE DLAMDA = (1/SUB)*MATMUL(TRANSPOSE(AV),DELSIG)
C		{a}**T{dS}/({a}**T*[DES]*{a})
C	********************************************************

	ATS = 0.0
	DO I = 1,6
	ATS = ATS + AV(I,1)*DELSIG(I)
	END DO

	DLAMD = ATS*(1.0/CONST)

	IF(DLAMD.LT.0.0) DLAMD = 0.0

C	...UPDATING TOTAL STRESS RETURN TO CLOSED YIELD SURFACE

	ATSIG = 0.0
	DO I=1,6
	ATSIG = ATSIG+AV(I,1)*STRESS(I)
	STRESS(I) = STRESS(I)+DELSIG(I)-DLAMD*DA(I,1)
	END DO

	EPSTN = EPSTN + DLAMD*ATSIG/CURS

C	*************************************
C	...BRING IT BACK TO THE YIELD SURFACE
C	*************************************

	CALL	SINVAR3D(STRESS,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA) 
	CALL	YIEDSF3D(VARI1,VARJ2,PHI,THETA,CURS)

C	...SETTING RETURN STRESS BACK COEFFICIENT FOR EACH INCRE. STEP

	CRETRN = 1.0
	
	CALL HARDMOHR(PSI,EPSTN)

C	SIGP = DLAMD*ATSIG/CURS
C	PREY = PREY + PSI*SIGP

	PREY = PREY + PSI*DLAMD*ATSIG/CURS
	 

	IF(CURS.GT.PREY)THEN
	CRETRN = PREY/CURS
	END IF

	DO I =1,6
	STRESS(I) = CRETRN*STRESS(I)
	END DO

C	***********************************
	END DO ! END OF SUBINCREMENTAL STEP
C	***********************************

C	...SET THE UPDATING CURRENT YIELD SURFACE AND FLOW VECTOR

	CALL	SINVAR3D(STRESS,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)  
	CALL	FLOWS3D(STRESS,AV,DA,DES,EPSTN,CONST,PHI)

C	... COMPUTE PLASTIC CONSTITUTIVE MATRIX

	DPS =  (1./CONST)*(MATMUL(DA,TRANSPOSE(DA)))

C	...COMPUTE ELASTO-PLASTIC CONSTITUTIVE MATRIX

	DO I = 1,6
	 DO J = 1,6
		DEP(I,J)=DES(I,J)-DPS(I,J)
	 END DO
	END DO

C	...UPDATING NEW SURFACE

	WK(1) = EPSTN
	WK(2) = PREY

C	*******************************
C	STORE STRESS TO PREVIOUS UPDATE
C	*******************************

C	...CHANGED DIRECTION BACK TO STORE DATA

	CALL	CONVSIG3D(STRESS,SIG)
	CALL	CONVSIG3D(STRESS,STRESS)

C	...UPDATING STRAIN

	DO I =1,6	
	EPS(I) = STRAIN(I)
	END DO

	GOTO 200

100	CONTINUE

C	*************************
C	ELASTIC SOIL CONSTITUTIVE
C	*************************

	DPS = 0.
	DES = 0.

	CALL DMATSLD_S(DES)

	DO I = 1,6
	DELEPS(I) = STRAIN(I) - EPS(I)
	END DO


	DELSIG = MATMUL(DES,DELEPS)

C	...CHANGE SIGN FOR POSITIVE CRITERION
	
	CALL  CONVSIG3D(DELSIG,DELSIG)
	
C	*************************************
C	COMPUTE CURRENT STRESS (TOTAL STRESS)
C	*************************************

	DO I =1,6
	STRESS(I) = SIG(I) + DELSIG(I)
	END DO

	DO I = 1,6
	 DO J = 1,6
		DEP(I,J)=DES(I,J)-DPS(I,J)
	 END DO
	END DO

C	...CHANGED DIRECTION BACK

	CALL	CONVSIG3D(STRESS,SIG)
	CALL	CONVSIG3D(STRESS,SONGST)
	DO I = 1,6
	STRESS(I) = SONGST(I)
	ENDDO

C	...UPDATING STRAIN

	DO I =1,6	
	EPS(I) = STRAIN(I)
	END DO


	GOTO 200
200	CONTINUE 


C1000	FORMAT(10E14.6)
C1001	FORMAT('FIRST STAGE')

	RETURN
	END

C	****************************************************

	SUBROUTINE CONVSIG3D(ST,STT)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	****************************************************

C	SUBROUTINE CONVERT DIRECTION FOR POSITIVE CROTERION

C	****************************************************

	DIMENSION ST(6),STT(6)

	STT(1) =  1.0*ST(1)   !SX
	STT(2) =  1.0*ST(2)   !SY
	STT(3) =  1.0*ST(3)   !SZ
	STT(4) =  1.0*ST(4)   !TXY
	STT(5) =  1.0*ST(5)   !TXZ
	STT(6) =  1.0*ST(6)   !TYZ

	RETURN
	END

C	*********************************************************

	SUBROUTINE FLOWS3D(ST,AV,DA,DES,EPSTN,CONST,PHI)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	****************************************************

C	SUBROUTINE COMPUTE FLOW VECTOR FOR CONTINUUM ELEMENT

C	****************************************************

	DIMENSION ST(6),CONS(1,1)
	DIMENSION AV(6,1),AV1(6,1),AV2(6,1),AV3(6,1),DA(6,1)
	DIMENSION DES(6,6)

	SX  = ST(1)
	SY  = ST(2)
	SZ  = ST(3)
	TXY = ST(4)
	TXZ = ST(5)
	TYZ = ST(6)

	AV1 = 0.
	AV2 = 0.
	AV3 = 0.
	AV  = 0.

	CALL SINVAR3D(ST,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)

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

C	...COMPUTE {a} FLOW VECTOR

	CF1 = (1./3.)*SIN(PHI)

C	...CONVERT PHI TO DEGREE PHI*180/Pi,REMOVE SINGULAR POINT
C	THETA ANGLE FROM PRICIPLE STRESS

	ABSTHE = ABS(THETA*57.2957795130824)

	IF(ABSTHE.LT.29.) GOTO 5
		
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

	CONS = MATMUL(TRANSPOSE(AV),MATMUL(DES,AV))
	CONST = CONS(1,1)

	CALL HARDMOHR(PSI,EPSTN)

	CONST = PSI + CONST
	DA =  MATMUL(DES,AV)

	RETURN
	END

C	************************************************

	SUBROUTINE SINVAR3D(ST,VARI1,VARI2,VARI3,VARJ2,VARJ3,ANGLE)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	***************************************************
C	COMPUTE STRESSES INVARIANT AND DEVIATORIC INVARIANT
C	FOR PREVIOUS STRESSES
C	***************************************************

	DIMENSION ST(6)


C	...DEFINE STRESSES

	SX  = ST(1)
	SY  = ST(2)
	SZ  = ST(3)
	TXY = ST(4)
	TXZ = ST(5)
	TYZ = ST(6)

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

C	...COMPUTE THE PRICIPLE IMAGINARY ANGLE

	SUBV = SQRT(VARJ2*VARJ2*VARJ2)

	IF(SUBV.EQ.0.0)THEN
	WRITE(*,*)' !!! PRINCIPLE STRESS DEVIDED BY ZERO'
	STOP
	END IF

	VALUE = -1.5*SQRT(3.0)*VARJ3/SUBV

	IF(VALUE.GT.1.0) VALUE = 1.0
	IF(VALUE.LT.-1.0) VALUE = -1.0
	
	ANGLE = (1.0/3.0)*ASIN(VALUE)

	RETURN
	END


C	***********************************************

	SUBROUTINE PRINSTR3D(SIGMA,PSTRES)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	******************************************

C	SUBROUTINE TO COMPUTE THE PRINCIPLE STRESS

C	****************************************** 

	DIMENSION SIGMA(6),PSTRES(3)

C	...SETTING CONSTANT AND DEFINE STRESSES

	PI = 3.141592654

	SX  = SIGMA(1)
	SY  = SIGMA(2)
	SZ  = SIGMA(3)
	TXY = SIGMA(4)
	TXZ = SIGMA(5)
	TYZ = SIGMA(6)

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
	
C	...COMPUTE THE PRICIPLE STRESSES

	PSTRES(1)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+2.*PI/3.)+VARI1/3. !MAX		 
	PSTRES(2)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+4.*PI/3.)+VARI1/3. !MIN	 
	PSTRES(3)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA)+VARI1/3.          !MID

C	**************************
C	S1 > S3 > S2 = MAX/MIN/MID
C	**************************

	RETURN
	END

C	*******************************************************************

	SUBROUTINE RFMOHR3D(RFACT,SIGMA,DSIGMAI,PHI,Y)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	*****************************************

C	SUBROUTINE COMPUTE MOHR COULOMB R-FACTOR

C	SIGI    = PREVIOUS STRESS INPUT
C	DELSIGI = DELTA STRESS INPUT
C	PHI	    = PRICTION ANGLE(RAD)
C	Y		= YIELD SURFACE
 
C	*****************************************

	DIMENSION SIGMAT(6),SIGMA(6),DSIGMAI(6)
	DIMENSION PSTRES(3)

C	...INITIALIZE

	SIGMAT = 0.
	PSTRES = 0.

C	...COMPUTE TOTAL STRESSES

	SIGMAT = SIGMA + DSIGMAI

	SX  = SIGMAT(1)
	SY  = SIGMAT(2)
	SZ  = SIGMAT(3)
	TXY = SIGMAT(4)
	TXZ = SIGMAT(5)
	TYZ = SIGMAT(6)

C	...SETTING CONTROL VARIBLE

	TOL = 1.0E-05
	K   = 0
	ITERATION = 50

C	*************************
C	   ITERATION SCHEME
C	*************************

100	CONTINUE

	K = K + 1

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
	
C	...COMPUTE THE PRICIPLE STRESSES

	PSTRES(1)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+2.*PI/3.)+VARI1/3. !MAX		 
	PSTRES(2)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA+4.*PI/3.)+VARI1/3. !MIN	 
	PSTRES(3)=(2.0*SQRT(VARJ2)/SQRT(3.0))*SIN(THETA)+VARI1/3.          !MID

C	...COMPUTE YIELD CRITERION 

	FS = PSTRES(1)-PSTRES(2)+(PSTRES(1)+PSTRES(2))*SIN(PHI)-Y

	IF(ABS(FS).LE.TOL) GOTO 200
	IF(K.EQ.ITERATION) GOTO 200

C	...DEFINE PREVIOUS STRESSES

	SX   = SIGMA(1)
	SY   = SIGMA(2)
	SZ   = SIGMA(3)
	TXY  = SIGMA(4)
	TXZ  = SIGMA(5)
	TYZ  = SIGMA(6)

C	...DEFINE INCREMENTAL STRESSES

	DSX  = DSIGMAI(1)
	DSY  = DSIGMAI(2)
	DSZ  = DSIGMAI(3)
	DTXY = DSIGMAI(4)
	DTXZ = DSIGMAI(5)
	DTYZ = DSIGMAI(6)

C	...DEFINE STRESSES COEFFICIENT

	TERMA1 = DTXY*DTXY-DSY*DSZ-DSX*DSZ+DTXZ*DTXZ-DSX*DSY+DTYZ*DTYZ
	TERMA2 = (1./3.)*(DSX+DSY+DSZ)*(DSX+DSY+DSZ)

	A = TERMA1 + TERMA2

	TERMB1 = -1.*(SX*DSZ+SZ*DSX+SX*DSY+SY*DSX+SY*DSZ+SZ*DSY)
	TERMB2 = (2./3.)*(SX+SY+SZ)*(DSX+DSY+DSZ)
	TERMB3 = 2.*TYZ*DTYZ+2.*TXZ*DTXZ+2.*TXY*DTXY

	B = TERMB1 + TERMB2 + TERMB3

	TERMC1 = (1./3.)*(SX+SY+SZ)*(SX+SY+SZ)
	TERMC2 = TYZ*TYZ+TXZ*TXZ+TXY*TXY-(SX*SY)-(SY*SZ)-(SX*SZ)

	C = TERMC1 + TERMC2

	TERMD1 = DSX+DSY+DSZ

	D = TERMD1

	TERME1 = SX+SY+SZ

	E = TERME1

C	...DEFINE PROPERTIES COEFFICIENT

	CHI = (-1.+SIN(PHI))/(1.+SIN(PHI))
	UPSILON = -(Y/(1.+SIN(PHI)))

C	...DEFINE PRICIPLE DIRECTION FOR MAXIMUM AND MINIMUM

	PI   = 3.141592654

	ANG1   = THETA + 2.*PI/3.
	THETA1 = SIN(ANG1)

	ANG2   = THETA + 4.*PI/3.
	THETA2 = SIN(ANG2)

C	...DEFINE COUPLE COEFFICIENT

	BETA  = (2./3.)*SQRT(3.)*(THETA1 + CHI*THETA2)
	ALPHA = (1./3.)*D + (1./3.)*CHI*D
	GMMA  = (1./3.)*E + (1./3.)*CHI*E + UPSILON 

C	...DEFINE STATE COEFFICIENT

	XI  = ALPHA/BETA
	ETA = GMMA/BETA

C	...DEFINE QUADRATIC COEFFICIENT

	FL = A - (XI*XI)
	FM = B - (2*XI*ETA)
	FN = C - (ETA*ETA)

C	...COMPUTE R FACTOR POSITIVE ROOT

	SONGVL = FM*FM-4.*FL*FN
	IF(SONGVL.LT.0.0) THEN
c	WRITE(*,*) 'SONGSAK KAK'
	R = 1.0
	GO TO 200
	ENDIF


	R = (1.0/(2.*FL))*(-FM+SQRT(FM*FM-4.*FL*FN))

C	...BACK SUBSTITUE FOR NEXT ITERATION

	SX = SIGMA(1)+ R*DSIGMAI(1)
	SY = SIGMA(2)+ R*DSIGMAI(2)
	SZ = SIGMA(3)+ R*DSIGMAI(3)

	TXY = SIGMA(4)+ R*DSIGMAI(4)
	TXZ = SIGMA(5)+ R*DSIGMAI(5)
	TYZ = SIGMA(6)+ R*DSIGMAI(6)

	GOTO 100

200	CONTINUE

	RFACT = R

C	WRITE(*,10)RFACT,K,FS
C10	FORMAT(1E14.5,I5,4E14.5)

	RETURN
	END

C	*********************************************************

	SUBROUTINE YIEDSF3D(VARI1,VARJ2,PHI,THETA,YIELD)    
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	*************************************
C	COMPUTE YIELD SURFACE AT GAUSS POINT
C	*************************************

	TERM1 = (1./3.)*VARI1*SIN(PHI)
	TERM2 = SQRT(VARJ2)*COS(THETA)
	TERM3 = SQRT(VARJ2)*(1./SQRT(3.))*SIN(THETA)*SIN(PHI)

	YIELD = TERM1+TERM2-TERM3

	RETURN
	END

C	========================================================
c	SUBROUTINE HARDMOHR(HD,EPSTN)
c	IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT INTEGER*4 (I-N)
C	========================================================
c
c	COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
c	COMMON /CRIP/  OCR,VO,DEN,VGL,SWL,SM,OK
c	COMMON /MOHR/  COH,ANG,BETA
c
c	HD = 0.0
c
c	PI  = 3.141592654 
c	PHI = ANG*PI/180.
c
c	IF(EPSTN.NE.0.0) HD = (1.-COS(PHI))*YM/(BETA*EPSTN)
c
C	IF(EPSTN.NE.0.0) HD = YM/(BETA*EPSTN)
c
c
c	WRITE(*,*) HD
c
c
c	RETURN
c	END
c
C	========================================================

C	=======================================================================
C	================== START OF VON-MISE NEW MODULE =======================
C	=======================================================================
C	=======================================================================
      SUBROUTINE DMMISE(SIGP,EPSP,IPEL,EPS,SIG,DP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     MAIN ROUTINE FOR ELASTO-PLASTIC MATERIAL BEHAVIOUR OF
C     THE 3-D SOLID ELEMENTS. THE THEORY IS BASED ON THE VON MISES
C     YIELD CONDITION AND USES NORMALITY RULE AND PRANDTL-REUSS FLOW
C     RULE TO CONVERT ELASTO-PLASTIC STRAIN INCREMENTS INTO STRESS INC
C	----------------------------------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     SIGP(6)    = STRESSES AT THE END OF THE PREVIOUS UPDATE
C     EPSP(6)    = STRAINS AT THE END OF THE PREVIOUS UPDATE
C     YIELD      = YIELD STRESS AT THE END OF THE PREVIOUS UPDATE
C     IPEL       = PLASTICITY FLAG (1=ELASTIC,2=PLASTIC)
C     EPS(6)     = CURRENT TOTAL STRAINS
C     SIG(6)     = RETURNED TOTAL STRESSES
C     DP(6,6)    = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DELEPS (6) = STRAIN INCREMENTS
C     DELSIG(6)  = STRESS INCREMENTS
C     DEPS(6)    = STRAIN SUBINCREMENTS (EQUIVALENCED WITH DELEPS)
C     DEPSE(6)   = ELASTIC PART OF STRAIN SUBINCREMENTS
C     DF(6)      = TOTAL DEVIATORIC STRESSES
C     DS(6)      = INCREMENTAL DEVIATORIC STRESSES
C     SX,SY...   = EQUIVALENCED WITH DF(6)
C     DX,DY...   = EQUIVALENCED WITH DS(6)
C     RATIO      = PART OF STRESS TAKEN ELASTICLY
C     PLAMDA     = PLASTIC STRAIN RATE MULTIPLIER
C     FT         = VON MISES YIELD FUNCTION
C     NINC,IINC  = NUMBER (COUNTER) OF SUBINCREMENTS
C	---------------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/,/HOOK/
C	---------------------------------------
C     NSINC      = VARIABLE CONTROLLING NUMBER OF SUBINCREMENTS
C     ITOLEY     = TOLERANCE ON EXCESS YIELD FUNCTION
C     YLD        = CURRENT UPDATED YIELD STRESS
C     ----------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP
	COMMON /HOOK_S/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION SIGP(6),EPSP(6),EPS(6),SIG(6),DP(6,6),STIFG(6,6),
	1		  DEP(6,6)
      DIMENSION DELEPS(6),DELSIG(6),DEPS(6),DEPSE(6),DF(6),DS(6)
C
      EQUIVALENCE (DELEPS,DEPS)
      EQUIVALENCE (SX,DF(1)),(SY,DF(2)),(SZ,DF(3)),(SXY,DF(4))
      EQUIVALENCE (SXZ,DF(5)),(SYZ,DF(6))
      EQUIVALENCE (DX,DS(1)),(DY,DS(2)),(DZ,DS(3)),(DXY,DS(4))
      EQUIVALENCE (DXZ,DS(5)),(DYZ,DS(6))
C     ------------------------------------------------------------
C     1. CALCULATE INCREMENTAL STRAINS
C     2. CALCULATE INCREMENTAL STRESSES ASSUMING ELASTIC BEHAVIOUR
C     3. CALCULATE TOTAL STRESSES ASSUMING ELASTIC BEHAVIOUR
C     ------------------------------------------------------------
	YOUNG = YM
	POISN = PR

      YLD2 = YLD*YLD
      DO 110  I=1,6
 110  DELEPS(I) = EPS(I) - EPSP(I)

      CALL ELTIFF(YOUNG,POISN,DP)
	DO 120  I = 1,6
	DELSIG(I) = 0.0
	DO 120  J = 1,6
 120	DELSIG(I) = DELSIG(I) + DP(I,J)*DELEPS(J)

      DO 130  I=1,6
 130  SIG(I) = SIGP(I) + DELSIG(I)
C     --------------------------------------------------
C     4. CALCULATE DEVIATORIC STRESSES AND CHECK WHETHER
C        STATE OF STRESS FALLS OUTSIDE THE YIELD SURFACE
C     --------------------------------------------------
      CALL DMDEVI (SIG,DF)
      F = 1.5*(SX*SX + SY*SY + SZ*SZ)+3.0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
      CONTINUE
	IF (F-YLD2)  410,410,450
C     --------------------------------------------------------
C     STATE OF STRESS WITHIN YIELD SURFACE - ELASTIC BEHAVIOUR
C     --------------------------------------------------------
 410  IPEL = 1

      GOTO 900
C     ---------------------------------------------------------
C     STATE OF STRESS OUTSIDE YIELD SURFACE - PLASTIC BEHAVIOUR
C     ---------------------------------------------------------
 450  IPEL = 2

      DO 460  I = 1,6
 460  SIG(I) = SIGP(I)
C     -------------------------------------------------------
C     5. DETERMINE PART OF STRAIN TAKEN ELASTICLY AND ADD
C        STRESSES DUE TO THE ELASTIC STRAINS TO PREVIOUS STRESSES
C     -----------------------------------------------------------
      CALL DMDEVI (SIG,DF)
      CALL DMDEVI (DELSIG,DS)
      A = DX*DX + DY*DY + DZ*DZ + 2.*(DXY*DXY + DXZ*DXZ + DYZ*DYZ)
      B = SX*DX + SY*DY + SZ*DZ + 2.*(SXY*DXY + SXZ*DXZ + SYZ*DYZ)
      E = SX*SX + SY*SY + SZ*SZ + 2.*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
     1  - 2.*YLD2/3.
      RATIO = (-B + SQRT(B*B-A*E))/A
      DO 480  I=1,6
 480  SIG(I) = SIG(I) + RATIO*DELSIG(I)
C     --------------------------------------------------
C     6. DETERMINE NUMBER AND MAGNITUDE OF SUBINCREMENTS
C     --------------------------------------------------
 600  NINC = 30*(SQRT(F)/YLD) + 1 ! NSINC*(SQRT(F)-YLD)/YLD + 1
      IF (NINC.GT.2*NSINC)  NINC = 2*NSINC
	ANINC = NINC
      FACT = (1.0 - RATIO) / ANINC
      DO 610  I=1,6
 610  DEPS(I) = FACT*DELEPS(I)
C     ----------------------------------------------------
C     7. CALCULATE PLASTIC STRAIN RATE MULTIPLIER AND
C        CONVERT ELASTIC PART OF STRAIN SUBINC. TO STRESSES
C     -----------------------------------------------------
      DO 890  IINC=1,NINC
      CALL DMDEVI (SIG,DF)
      PLAMDA = (SX*DEPS(1) + SY*DEPS(2) + SZ*DEPS(3) + SXY*DEPS(4) +
     1          SXZ*DEPS(5) + SYZ*DEPS(6)) * 0.5/YLD2
      IF (PLAMDA.LT.0.0)  PLAMDA = 0.0
      DO 710  I=1,6
 710  DEPSE(I) = DEPS(I) - 3.0*PLAMDA*DEPS(I)
 
	DO 715  I = 1,6
	DELSIG(I) = 0.0
	DO 715  J = 1,6
 715	DELSIG(I) = DELSIG(I) + DP(I,J)*DEPSE(J)


      DO 720  I=1,6
 720  SIG(I) = SIG(I) + DELSIG(I)
C     ----------------------------------------------
C     8. STRESS CORRECTION BACK TO THE YIELD SURFACE
C     ----------------------------------------------
      CALL DMDEVI (SIG,DF)
      F = 1.5*(SX*SX + SY*SY + SZ*SZ)+3.0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
C	ABCD=1./(DFLOAT(ITOLEY))
	ABCD = 1./1000000.
      IF (F-YLD2.LE.ABCD)  GOTO 890
c     COEF = -1.0 + SQRT(YLD2/F)
	COEF = SQRT(YLD2)/SQRT(F)
      DO 810  I = 1,6
 810  SIG(I) = COEF*SIG(I)
C
 890  CONTINUE
C     ----------------------------------------------
C     9. UPDATE PREVIOUS STRESSES, STRAINS AND YIELD
C     ----------------------------------------------
900	CONTINUE

	DO I = 1,6
      SIGP(I) = SIG(I)
	EPSP(I) = EPS(I)
	ENDDO

C     ------------------------------------
C     FORM THE ELASTO-PLASTIC MATERIAL LAW
C     ------------------------------------
      CALL DMDELP(SIG,DEPS,DEP,DP,IPEL,YOUNG,POISN,YLD)
	

	DO I = 1,6
	DO J = 1,6
	DP(I,J) = DEP(I,J)
	ENDDO
	ENDDO

C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DMDEVI (SIG,S)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CALCULATES DEVIATORIC STRESSES  S(IJ) = SIG(IJ)-1/3SIG(KK)D(IJ)
C	---------------------------------------------------------------
C     SIG(6)    = STRESSES
C     S(6)      = DEVIATORIC STRESSES
C     ---------------------------------------------------------------
      DIMENSION SIG(6),S(6)
C
      SM = (SIG(1) + SIG(2) + SIG(3)) / 3.0
      DO 100  I=1,3
      S(I)   = SIG(I) - SM
 100  S(I+3) = SIG(I+3)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DMDELP(SIG,DEPS,DEP,DP,IPEL,YOUNG,POISN,YLD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     FORMS THE ELASTO-PLASTIC STRESS-STRAIN MATRIX FOR SOLID ELEMENTS
C	----------------------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     SIG(6)     = CURRENT TOTAL STRESSES
C     DEPS(6)    = SUBINCREMENTS IN STRAIN
C     DP(6,6)    = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     CP(36)     = DP(6,6)
C     IPEL       = PLASTICITY FLAG (1=ELASTIC,2=PLASTIC)
C	---------------
C     LOCAL VARIABLES
C	---------------
C     S(6)       = DEVIATORIC STRESSES
C     WP         = INCREMENTAL ELASTO-PLASTIC WORK
C
C     FOR VARIABLES IN COMMON BLOCK /HOOK/ SEE ROUTINE HOKLAW
C     --------------------------------------------------------------
C
      DIMENSION SIG(6),DEPS(6),DP(6,6),CP(36),S(6),STIFG(6,6),DEP(6,6)
C     --------------------
C     ELASTIC CONTRIBUTION
C     --------------------

      IF (IPEL.EQ.1)  GOTO 300
C     -----------------------------------------------------
C     CALCULATE DEVIATORIC STRESSES AND CHECK FOR UNLOADING
C     -----------------------------------------------------
      CALL DMDEVI (SIG,S)
      WP = 0.
      DO 100  I=1,6
 100  WP = WP + S(I)*DEPS(I)
      IF (WP.LT.0.)  GOTO 300
C     ---------------------------------------------
C     DEDUCT PLASTIC CONTRIBUTION IN LOWER TRIANGLE
C     ---------------------------------------------
C      C2BETA = C2*1.5/YLD/YLD
      DO 200  I=1,6
      DO 200  J=1,I
 200  DP(I,J) = DP(I,J)  ! - C2BETA*S(J)*S(I)
C     --------------------------------------------
C     FILL IN SYMMETRIC ELEMENTS IN UPPER TRIANGLE
C     --------------------------------------------
 300  DO 390  I=1,6
      DO 390  J=I,6
 390  DP(I,J) = DP(J,I)

	DO I = 1,6
	DO J = 1,6
	DEP(I,J) = DP(I,J)
	ENDDO
	ENDDO
C
      RETURN
      END
C
C	=====================================================================
C	=======================================================================
C	=======================================================================
	SUBROUTINE	ELTIFF(YOUNG,POISN,STIFG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=======================================================================
	DIMENSION STIFG(6,6)
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
	

	RETURN
	END	


C	=======================================================================
C	=======================================================================

C	=======================================================================
C	==================== END OF VON-MISE NEW MODULE =======================
C	=======================================================================

C	*******************************************************************

	SUBROUTINE RFMOHR3D_NEW(RFAC,SG,DSG,PHI,YIELD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C	*****************************************

C	SUBROUTINE COMPUTE MOHR COULOMB R-FACTOR

C	SIGI    = PREVIOUS STRESS INPUT
C	DELSIGI = DELTA STRESS INPUT
C	PHI	    = PRICTION ANGLE(RAD)
C	Y		= YIELD SURFACE
 
C	*****************************************


	DIMENSION SG(6),DSG(6)
	DIMENSION SGR(6)

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

	CALL	SINVAR3D(SGR,VARI1,VARI2,VARI3,VARJ2,VARJ3,THETA)  
	CALL	YIEDSF3D(VARI1,VARJ2,PHI,THETA,FR)

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

C	==============================================================
	SUBROUTINE KDELTA_S (K,I,J)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	K = 0
	IF(I.EQ.J) K = 1

	RETURN

	END
C	==============================================================
C	==============================================================

