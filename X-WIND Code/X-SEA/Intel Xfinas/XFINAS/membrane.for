C========================================================================
C	START 4 NODE MEMBRANE ELEMENT (PLANE STRESS & STRAIN)
C========================================================================
C
C======================================================================
      SUBROUTINE MEMBRA2(PROPM,PROPG,NODEX,WA,S,COORD,
	1                   EDIS,EDISI,RE,MWG,FIN)
C	FIN - ADDED TO PREVIOUS LINE BY GILSON - JUL2003 (INT FORCE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     COMPUTES ELEMENT FORCE VECTOR BALANCED IN CURRENT DISPLACEMENT
C     FIELD. IN CASE IFREF=0 THE STIFFNESS MATRIX IS UPDATED FOR
C	--------------------------------------
C     MATERIAL NONLINEARITIES    (NOLIN=1,3)
C     GEOMETRIC NONLINEARITIES   (NOLIN=2,3)
C     FINDING LINEAR STIFFNESS   (NOLIN=0)
C	---------------
C     INPUT VARIABLES
C	---------------
C     COORD(3,8)     = NODAL COORDINATES X,Y Z
C     PROP(5)        = MEMBER PROPERTIES (SEE ROUTINE CONSTI)
C     NODEX(NNO-4)   = POSITIONS OF MIDSIDE NODES
C     EDIS(16)       = TOTAL DISPLACEMENTS AT ELEMENT NODES
C     EDISI(16)      = NODAL ELEMENT DISPLACEMENT INCREMENTS
C     IE             = ELEMENT NUMBER
C     NI             = NUMBER OF NODES FOR THIS ELEMENT
C     MW             = DIMENSION OF WORKING ARRAY STORING 4 STRESSES
C                      (PLUS 4 STRAINS AND A PLASTICITY FLAG NOLIN=1,3)
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     RE(16)         = ELEMENT FORCE VECTOR (TWO EQUIL. LOADS PER NODE)
C     S(136)         = NEW TANGENTIAL STIFFNESS MATRIX (UPPER TRIANG)
C     WA(MW,NGR*NGS) = WORKING ARRAY STORING STRESS,(STRAIN) COMPONENTS
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/
C	--------------------------------
C     NAME(2)      = NAME OF ELEMENT MODULE
C     ITYPE        = CODE NUMBER FOR ELEMENT MODULE (9)
C     ISTYP        = ELEMENT SUB-TYPE
C     NLOPT        = CODE FOR NONLINEAR OPTION
C     NLOPT = 0      LINEAR ANALYSIS
C     NLOPT = 1      MATERIALLY NONLINEAR ONLY
C     NLOPT = 3      UPDATED LAGRANGIAN
C     MTMOD        = CODE FOR MATERIAL TYPE
C     MTMOD = 1      LINEAR ELASTIC (HOOKEAN)
C     MTMOD = 4      ELASTO-PLASTIC (VON-MISES)
C     NELE         = NUMBER OF ELEMENTS IN THIS GROUP
C     NMPS         = NUMBER OF MATERIAL PROPERTY SETS
C     NGPS         = NUMBER OF GEOMETRIC PROPERTY SETS
C     NMP          = NUMBER OF MATERIAL PROPERTIES PER SET
C     NGP          = NUMBER OF GEOMETRIC PROPERTIES PER SET
C     NNM          = MAXIMUM NUMBER OF NODES FOR ANY ONE ELEMENT
C     NEX          = MAXIMUM NUMBER OF EXCESS NODES
C     NCO          = NUMBER OF NODAL COORDINATES
C     NNF          = NUMBER OF NODAL DEGREES OF FREEDOM
C     NEF          = MAXIMUM NUMBER OF ELEMENT DEGREES OF FREEDOM
C     NWG          = NUMBER OF STORAGE LOCATIONS AT EACH GAUSS POINT
C     NPT          = NUMBER OF GAUSS POINTS
C     NWA          = SIZE OF WORKING ARRAY
C     NWS          = SIZE OF ELEMENT STIFFNESS MATRIX
C     MEL          = CURRENT ELEMENT NUMBER
C     NNO          = NUMBER OF NODES FOR THIS ELEMENT
C     NEF          = NUMBER OF DEGREES OF FREEDOM FOR THIS ELEMENT
C     NELTOT       = TOTAL NUMBER OF ELEMENTS (ALL GROUPS)
C	MTYP		 = FLAG FOR MEMBRANE TYPE
C					(0) - ISOPARAMETRIC
C					(1) - QUASI-CONFORMING
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /GAUS/
C	--------------------------------
C     GLOC(10,10)  = NATURAL GAUSS POINT COORDINATES (1*1 TO 10*10)
C     GWT (10,10)  = GAUSS POINT WEIGHTS
C     NGR          = NUMBER OF GAUSS POINTS ALONG BEAM AXIS (= NPT)
C     NGT          = NUMBER OF CROSS-SECTIONAL INTEGRATION STATIONS
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /FLAG/
C	--------------------------------
C     IFREF        = FLAG FOR REFORMATION OF STIFFNESS (IFREF=0)
C     IFEIG        = FLAG FOR EIGENVALUE SOLUTION (IFEIG=0)
C     ITASK = 1      FIRST ENTRY INTO ELEMENT MODULE
C     ITASK = 2      ENTRY DURING EQUILIBRIUM ITERATIONS
C     ITASK = 3      ENTRY TO WORK OUT STRESSES (LAST STEP ONLY)
C     ITASK = 4      ENTRY TO DETERMINE GEOMETRIC STIFF. MATRIX ONLY
C	------------------------------
C     VARIABLES IN /ELEM/ AND /HOOK/
C	------------------------------
C     YLD       = CURRENT UPDATED YIELD STRESS
C     ISR       = NUMBER OF STRAIN COMPONENTS
C     IST       = NUMBER OF STRESS COMPONENTS
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DISD(5)        = DISPLACEMENT DERIVATIVES
C     STRAIN(4)      = LINEAR OR NONLINEAR STRAIN VECTOR
C     STRESS(4)      = LINEAR OR NONLINEAR STRESS VECTOR
C     D(6,6)         = STRESS - STRAIN MATRIX
C     ----------------------------------------------------------------
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/ A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
c      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
C
C	ADDED BY DE SILVA
	COMMON /HARD/  HP,DEN

	COMMON /DYNA/ CDEN,IMASS  !ADD THIS LINE SONGSAK MAY2006
C
      DIMENSION PROPM(5),PROPG(*),NODEX(*),WA(MWG,*),S(78),COORD(2,8)
      DIMENSION EDIS(12),EDISI(1),RE(12),REDIS(48),COORDI(2,8)
      DIMENSION DR(64),H(8),HD(2,8),XJI(4),HR(8),HS(8)
      DIMENSION DISD(5),EPS(8),EPSQ(8),SIGR(8),FF1(4)
      DIMENSION VR(2),VS(2),EC(3),SL(4,2),FF(6),NOD(4)
C      DIMENSION VR(3),VS(3),VT(3),EC(3),SL(4,2),FF(6),NOD(4)
      DIMENSION COVR(3),COVS(3)
	DIMENSION BA(4,120),FJ(4),RRN(4),SSN(4)
	DIMENSION D(36),P(2,8),B(4,16),STRAIN(4),STRESS(4),TAU(4)
	DIMENSION ST(78)

C	COMMON TO ALL STRAIN TERMS
	DIMENSION RS(2,4),SLN(4),RN(4),SN(4),BL(4),TRS(4)
	DIMENSION FLR(4),FLS(4),FLSS(4),FLRR(4),ANT(4)
	DIMENSION DWR(12),DWRS(12),DWS(12),DWSR(12)
	DIMENSION T(12,12),TT(12,12)
C	FOR BENDING STRAINS
	DIMENSION AB1(4,4),AB2(3,3),CB(11,24),ACB(11,24),PDPB(11,11)
	DIMENSION AKBS(300)
C	FOR MEMBRANE STRAINS
	DIMENSION CM(9,12),ACM(9,12),AM(9,9),PSP(9,9)
	DIMENSION APM(9)
	DIMENSION AMM2(9,9),TEMM1(9,12),ACM2(9,12),TEMM2(12,9)
	DIMENSION AMMI2(9,9),AKM(12,12)
C	FOR GEOMETRIC STIFFNESS
	DIMENSION CG(4,12),ACG(4,12),PFP(4,4),SIGA(8)
C	INTERNAL FORCE VECTOR
	DIMENSION FM(12),FB(24),FS(24),APSH(2)
	DIMENSION PAPM(9),PAPB(11),PAPSH(2)
C	LINEAR STIFFNESS - INELASTIC CASE
	DIMENSION PMPM(9,9),PMPB(9,11),PBPM(11,9),PBPB(11,11),PQPQ(2,2)
C	FOR PLASTICITY
	DIMENSION ANTR(4),ANTS(4),ANTRS(4)
	DIMENSION ANTRR(4),ANTSS(4),ANTRRS(4),ANTSSR(4),ANTRRSS(4)
	DIMENSION ANTRRR(4),ANTSSS(4),ANTRRRR(4),ANTSSSS(4)
	DIMENSION TEDIS(12)

	DIMENSION AX(5),AR(5),AS(5),ARR(5),ASS(5),ARS(5)
	DIMENSION ARRR(5),ASSS(5),ARRS(5),ASSR(5)
	DIMENSION ARRRR(5),ASSSS(5),ARRSS(5)
C
	DIMENSION CRD(3,8),CRDI(3,8)

C	DISPLACEMENT DERIVATIVES WRT REFERENCE COORDINATES
	DIMENSION DIRC(5)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)

      EQUIVALENCE (APEL,IPEL)

C     --------------
C     INITIALIZATION
C     --------------
      PI=3.1415926535898
	DVOL=1.0
	IPEL = 1
	KST = 0
C     -----------------------------------
C     OBTAIN LINEAR STRESS - STRAIN LAW
C     MTMOD=1, IPEL=1 FOR LINEAR ANALYSIS
C	HOKLAW - FLAG=2 FOR PLANE STRESS
C	ISTYP = 1 - AXISYMMETRIC
C	ISTYP = 1 - PLANE STRAIN
C	ISTYP = 2 - PLANE STRESS
C	-----------------------------------
10	ISR = 3
      IST = 3
      IF (ISTYP.LE.1)  IST = 4
      IF (ISTYP.EQ.0)  ISR = 4
      CALL HOKLAW (PROPM,PROPG,ISTYP)
      CALL DELA3D (D,ISTYP)
	IF (MTMOD.LT.3) D = D*TH

C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C     ------------------------------------------------------
      CALL MEMINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)

C	-----------------------------------------------
C	DETERMINE (1)LOCAL CORDINATE VECTORS - VR,VS,VT
C	          (2)LOCAL NODAL COORDINATES - RS
C	-----------------------------------------------
	CALL VRS(COORD,VR,VS,RS,NNO)

	IF (ITASK.NE.5) GOTO 1333

C	----------------
C	FORM MASS MATRIX
C	----------------
	IF (IMASS.EQ.1) THEN
C	  LUMPED MASS MATRIX
	  AREA = (-RS(2,4)/2+RS(2,2)/2)*RS(1,1)+(-RS(2,1)/2+RS(2,3)/2)*R
     #          S(1,2)+(RS(2,4)/2-RS(2,2)/2)*RS(1,3)+(RS(2,1)/2-RS(2,3)/
     #          2)*RS(1,4)
        CALL MASSQC (S,H,PROPM(5),AREA*TH,NNO,NEF,IPT)
	ELSE
C	  CONSISTENT MASS MATRIX
C       LOOP OVER GAUSS POINTS
        DO 911  IGR=1,NGR
          Rg = GLOC(IGR,NGR)
          DO 911  IGS=1,NGS
            Sg = GLOC(IGS,NGS)
            WT = GWT(IGR,NGR)*GWT(IGS,NGS)
C           SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (P)
            CALL SHAP2D (Rg,Sg,H,P,NODEX,NNO)
C           INVERSE JACOBIAN (XJI) , JACOBIAN DETERMINANT (DET)
C           AND STRAIN-DISPLACEMENT MATRIX (B)
            CALL JACO2D (MEL,NNO,COORD,P,XJI,DET,NCO)
            CALL MEBMAT (COORD,H,P,XJI,B,ISTYP,XBAR,NNO)
            IF (ISTYP.NE.0)  XBAR = TH
            DVOL = WT*DET*XBAR
            CALL MASSQC (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
911	  CONTINUE
	ENDIF

	RETURN

C	----------------------------
C	FORM TRANSFORMATION MATRIX T
C	----------------------------
1333	DO 30 I = 1,2
	  TT(1,I) = VR(I)
	  TT(2,I) = VS(I)
30	CONTINUE

C	------------------------------------------
C	COMPUTE FOR ELEMENT BOUNDARY LENGTHS - SLN
C	------------------------------------------
	SLN(1) = DSQRT((RS(1,2)-RS(1,1))**2.0+(RS(2,2)-RS(2,1))**2.0)
	SLN(2) = DSQRT((RS(1,3)-RS(1,2))**2.0+(RS(2,3)-RS(2,2))**2.0)
	SLN(3) = DSQRT((RS(1,4)-RS(1,3))**2.0+(RS(2,4)-RS(2,3))**2.0)
	SLN(4) = DSQRT((RS(1,1)-RS(1,4))**2.0+(RS(2,1)-RS(2,4))**2.0)

C	-----------------------------------------
C	COMPUTE FOR ELEMENT BOUNDARY NORMALS - RN
C	-----------------------------------------
	RN(1) = (RS(2,2)-RS(2,1))/SLN(1)
	RN(2) = (RS(2,3)-RS(2,2))/SLN(2)
	RN(3) = (RS(2,4)-RS(2,3))/SLN(3)
	RN(4) = (RS(2,1)-RS(2,4))/SLN(4)

C	------------------------------------------
C	COMPUTE FOR ELEMENT BOUNDARY TANGENTS - SN
C	------------------------------------------
	SN(1) = (RS(1,1)-RS(1,2))/SLN(1)
	SN(2) = (RS(1,2)-RS(1,3))/SLN(2)
	SN(3) = (RS(1,3)-RS(1,4))/SLN(3)
	SN(4) = (RS(1,4)-RS(1,1))/SLN(4)

C	------------------------------
C	AREA INTEGRATION FACTORS - ANT
C	------------------------------
	ANT(1) = ((RS(2,2)/6-RS(2,4)/6)*RS(1,1)+(-RS(2,1)/6+RS(2,4)/12+
	1         RS(2,3)/12)*RS(1,2)+(-RS(2,2)/12+RS(2,4)/12)*RS(1,3)+
     2         (-RS(2,2)/12-RS(2,3)/12+RS(2,1)/6)*RS(1,4))
	ANT(2) = ((RS(2,2)/6-RS(2,4)/12-RS(2,3)/12)*RS(1,1)+(-RS(2,1)/6+
	1         RS(2,3)/6)*RS(1,2)+(RS(2,1)/12+RS(2,4)/12-RS(2,2)/6)*
     2         RS(1,3)+(-RS(2,3)/12+RS(2,1)/12)*RS(1,4))
	ANT(3) = ((-RS(2,4)/12+RS(2,2)/12)*RS(1,1)+(-RS(2,1)/12-RS(2,4)/12
     1         +RS(2,3)/6)*RS(1,2)+(RS(2,4)/6-RS(2,2)/6)*RS(1,3)+
     2         (-RS(2,3)/6+RS(2,1)/12+RS(2,2)/12)*RS(1,4))
	ANT(4) = ((-RS(2,4)/6+RS(2,3)/12+RS(2,2)/12)*RS(1,1)+(-RS(2,1)/12+
	1        RS(2,3)/12)*RS(1,2)+(-RS(2,2)/12+RS(2,4)/6-RS(2,1)/12)*
     2        RS(1,3)+(RS(2,1)/6-RS(2,3)/6)*RS(1,4))

	IF (MTMOD.GT.2.AND.(NLOPT+ITASK).NE.1) THEN
	  CALL COMF(RS,AX,AR,AS,ARR,ASS,ARS,ARRR,ASSS,ARRS,ASSR,
	1                ARRRR,ASSSS,ARRSS)
	ENDIF

C	-------------------
C	COMPUTE FOR AB1,AB2
C	-------------------
	CALL ABEND(RS,SLN,RN,SN,AB1,AB2)

C	------------------------------------
C	COMP FOR THE MEMBRANE AND DRILLING
C	CONTRIBUTION TO THE LINEAR STIFFNESS
C	------------------------------------
C	Cm FOR MEMBRANE AND DRILLING STRAINS - CM
	CALL CMEM(RS,RN,SN,SLN,CM)
C	MULTIPLY (Am-1)(Cm)  - ACM
	CALL ABCM(CM,RS,AM,ACM)

C     -------------------------------------------------
C     REMOVE RIGID BODY TRANSLATIONS AND ROTATIONS FROM
C     TOTAL DISPLACEMENT VECTOR (NLOPT>1)
C     -------------------------------------------------
	IF (NLOPT.LE.1) GOTO 1800
	CALL MEMDSP(COORD,COORDI,EDIS,VR,VS,REDIS,NNO)

C	--------------------------------
C	DETERMINE STRAIN PARAMETER ALPHA
C	--------------------------------
1800	IF (NLOPT+ITASK.NE.1) THEN
C		TRANSFORM REDIS TO LOCAL - TEDIS
		K = 1
	    DO 1810 I = 1,4
	      TEDIS(K)   = VR(1)*REDIS(K)+VR(2)*REDIS(K+1)
	      TEDIS(K+1) = VS(1)*REDIS(K)+VS(2)*REDIS(K+1)
		  TEDIS(K+2) = REDIS(K+2)
		  K = K+3
1810		CONTINUE
		APM = MATMUL(ACM,TEDIS)
	ENDIF

C	-----------------------------------------
C	LOOP TO DETERMINE NODAL STRESS RESULTANTS
C	-----------------------------------------
	IPT=0

	DO 1000 II=1,5

	IPT=IPT+1
      IF (NLOPT+ITASK.EQ.1) GOTO 6767

C	-------------------------------------------------
C     DISPLACEMENT DERIVATIVES WRT MATERIAL COORDINATES
C	-------------------------------------------------
 200  CALL CLEARA (DISD,5)

	IF (II.NE.5) THEN
	  DISD(1) =APM(1)+RS(2,II)*APM(2)+RS(1,II)*APM(3)+
	1                                            RS(2,II)**2*APM(4)
	  DISD(2) =APM(5)+RS(1,II)*APM(6)+RS(2,II)*APM(7)+
	2                                            RS(1,II)**2*APM(8)
	  DISD(3) =-RS(2,II)*APM(3)/2.0-RS(1,II)*APM(7)/2.0+APM(9)/2.0
	  DISD(4) =-RS(2,II)*APM(3)/2.0-RS(1,II)*APM(7)/2.0+APM(9)/2.0
	ELSE
	  DISD(1) =APM(1)
	  DISD(2) =APM(5)
	  DISD(3) =APM(9)/2.0
	  DISD(4) =APM(9)/2.0
	ENDIF

C	--------------------------------------------------
C	DISPLACEMENT DERIVATIVES WRT REFERENCE COORDINATES
C	--------------------------------------------------
	IF (MTMOD.EQ.6) THEN
	  DETER = (1-DISD(1))*(1-DISD(2))-DISD(2)*DISD(3)
	  DIRC(1) = DISD(1)*(1-DISD(2))/DETER
	  DIRC(2) = DISD(2)*(1-DISD(1))/DETER
	  DIRC(3) = DISD(3)*(1-DISD(1))/DETER
	  DIRC(4) = DISD(4)*(1-DISD(2))/DETER
	  DIRC(5) = 0.0
	ENDIF

      IF (ISTYP.GE.1)  GOTO 300
      DO 290  IEF=1,NEF,2
 290  DISD(5) = DISD(5) + B(4,IEF)*EDIS(IEF)

C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
300	IF (MTMOD.NE.6) THEN
C	  ALMANSI STRAINS
C	  ---------------
	  STRAIN(1) = DISD(1)
        STRAIN(2) = DISD(2)
        STRAIN(3) = DISD(3)+DISD(4)
        STRAIN(4) = 0.
	ELSE
C	  GREEN-LAGRANGE STRAINS
C	  ----------------------
	  STRAIN(1) = DIRC(1)
        STRAIN(2) = DIRC(2)
        STRAIN(3) = DIRC(3)+DIRC(4)
        STRAIN(4) = 0.
	ENDIF

      IF (ISTYP.EQ.0)  STRAIN(4) = DISD(5)

C     --------------------------------------------------------------
C     FOR NLOPT>1 SUBSTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     --------------------------------------------------------------
C      IF (NLOPT.LE.1)  GOTO 400
      IF (NLOPT.LE.1 .AND. MTMOD.NE.6)  GOTO 400

123	IF (MTMOD.NE.6) THEN
C	  ALMANSI STRAINS
C	  ---------------
        STRAIN(1) = STRAIN(1) - 0.5*(DISD(1)*DISD(1) + DISD(4)*DISD(4))
        STRAIN(2) = STRAIN(2) - 0.5*(DISD(2)*DISD(2) + DISD(3)*DISD(3))
        STRAIN(3) = STRAIN(3) -     (DISD(3)*DISD(1) + DISD(4)*DISD(2))
	ELSE
C	  GREEN-LAGRANGE STRAINS
C	  ----------------------
        STRAIN(1) = STRAIN(1) + 0.5*(DIRC(1)*DIRC(1) + DIRC(4)*DIRC(4))
        STRAIN(2) = STRAIN(2) + 0.5*(DIRC(2)*DIRC(2) + DIRC(3)*DIRC(3))
        STRAIN(3) = STRAIN(3) +     (DIRC(3)*DIRC(1) + DIRC(4)*DIRC(2))
	ENDIF

      IF (ISTYP.EQ.0)  STRAIN(4) = STRAIN(4) - 0.5*(DISD(5)*DISD(5))

C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------
 400  GOTO (410,410,430,430,450,460), MTMOD
C	LINEAR ELASTIC
C	--------------
 410  CALL MEMSIG (STRAIN,STRESS,ISTYP)
      DO 415  I=1,4
 415  WA(I,IPT) = STRESS(I)
      GOTO 480

C	ELASTIC-PLASTIC
C	---------------
C430	CALL VMISES (WA(1,IPT),WA(5,IPT),IPL,STRAIN,STRESS,D,ISTYP)
C	WA(10,IPT) = IPL
C      IPEL=2


430	IF (HP.EQ.0.) THEN
	CALL VMISES (WA(1,IPT),WA(5,IPT),IPL,STRAIN,STRESS,D,ISTYP)
	ELSE
	CALL VMISESH (WA(1,IPT),WA(5,IPT),STRAIN,STRESS,D,ISTYP,
     #              WA(9,IPT))
	END IF
	WA(10,IPT) = IPL
      IPEL=2

	GOTO 480

460	CALL MOONEY (DIRC,D,PROPM,STRESS)
	DO 1710  I=1,4
	  WA(I,IPT) = STRESS(I)
	  WA(I+4,IPT) = STRAIN(I)
1710	CONTINUE
	WA(10,IPT) = IPL
      IPEL=2

C     ----------------------------------------
C     ADJUST THICKNESS (NLOPT>1 .AND. ISTYP=2)
C     ----------------------------------------
 480  IF (ITASK.EQ.3)  GOTO 5000
      THICK = TH
C      IF (NLOPT.EQ.1.OR.ISTYP.LT.2)  GOTO 485
      IF ((NLOPT.EQ.1.OR.ISTYP.LT.2).AND.MTMOD.NE.6)  GOTO 485
      STRAIN(4) = -(D(6)*STRAIN(1)+D(12)*STRAIN(2)+D(18)*STRAIN(3))
     1            /D(36)
      EXT  = 1. - 2.*STRAIN(4)
      EXT  = 1. - SQRT(EXT)
C
      THICK = THICK*(1.-EXT)

C	-----------------------------------
C	PRE-INTEGRATE THROUGH THE THICKNESS
C	-----------------------------------
485	DO 487 I=1,4
	  STRESS(I) = STRESS(I)*THICK
487	CONTINUE
	DO 488 I=1,36
	  D(I) = D(I)*THICK
488	CONTINUE

470	CALL PSMEM(IPT,STRESS,AX(II),AR(II),AS(II),ARR(II),ASS(II),PAPM)
	CALL PDMEM(IPT,D,AX(II),AR(II),AS(II),ARS(II),ARR(II),ASS(II),
	1          ARRS(II),ASSR(II),ARRSS(II),ARRR(II),ASSS(II),
	2          ARRRR(II),ASSSS(II),PSP)

 450  CONTINUE

5000	CONTINUE

1000	CONTINUE

C	NEXT LINE CHANGED BY GILSON - JUL2003 (INT FORCE)
C      IF (ITASK.LE.2) GOTO 6767
      IF (ITASK.LE.3.AND.ISOLOP.NE.4) GOTO 6767
C	FOR LINEAR BUCKING ANALYSIS
      IF (IFEIG.EQ.0) GOTO 8000

C	---------------------------------
C	DETEREMINE TRANSPOSE(P)*MODULUS*P
C	MEMBRANE - Nm*S*Pm
C	---------------------------------
6767	IF (IPEL.EQ.1) THEN
C	  ELASTIC CASE
	  CALL NMSPM(AB1,AM,D,RS,PSP)
	ELSE
C	  ELASTIC-PLASTIC CASE - FILL LOWER TRIANGLE
	  DO 700 I = 1,9
	   DO 705 J = I,9
	    PSP(J,I) = PSP(I,J)
705	   CONTINUE
700	  CONTINUE
C	  PSP = PSP*TH
	ENDIF
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
	TEMM2 = TRANSPOSE(ACM)	

      IF (NLOPT+ITASK.EQ.1) GOTO 7000

	IF (IPEL.EQ.1) PAPM = MATMUL(PSP,APM)
c	IF (IPEL.EQ.2) PAPM =PAPM*TH
	IF (IPEL.EQ.2) PAPM =PAPM
	FM = MATMUL(TEMM2,PAPM)

C	TRANSFORM FM TO GLOBAL - RE
	K = 1
	DO 333 I = 1,4
	  RE(K)   = VR(1)*FM(K)+VS(1)*FM(K+1)
	  RE(K+1) = VR(2)*FM(K)+VS(2)*FM(K+1)
	  RE(K+2) = FM(K+2)
	  K = K+3
333	CONTINUE

C     ------------------------------------------------------
C     FIND LINEAR CONTRIBUTION TO STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------------------
	IF (IFREF ) 9000,7000,9000
7000	ST = 0.0
	KST = 1
	CALL KMEM(ACM,PSP,ST)

C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 8100
8000	IF (KST.EQ.0) ST = 0.0
C	- -- -- -- -- -- -- -- -- -- -- --
C	COMPUTE FOR INVERSE(Ag)*(Cg) - ACG
C	- -- -- -- -- -- -- -- -- -- -- --
	CALL BCGEO(AM,CM,ACG)

C	-- -- -- -- -- -- -- -- -- -- --
C	DETERMINE INT(Ng*F*Pg)drds - PFP
C	-- -- -- -- -- -- -- -- -- -- --
5555	IF (IPEL.EQ.1) THEN
C	  ELASTIC CASE
	  CALL NGFPG(AM,APM,D,PFP)
	ELSE
C	  PLASTIC CASE
	  CALL NFPPLAS(PAPM,PFP)
	ENDIF

C	-- -- -- -- -- -- -- -- -- --
C	DETERMINE GEOMETRIC STIFFNESS
C	-- -- -- -- -- -- -- -- -- --
	CALL KGEO(ACG,PFP,ST)

8100	CALL KGMEM(ST,TT,S)
C	TIM(12)=TIM(12)+TIM2
9000	CONTINUE

C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF

      RETURN
      END
C
C=====================================================================
      SUBROUTINE PSMEM(IPT,SIGR,A,AR,AS,ARR,ASS,PAPM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTES FOR INTEGRAL(TRANS(P)*SIGR)dA
C
C	INPUT VARIABLES
C	SIGR(8)  = NODAL STRESS RESULTANTS
C	A   = FACTOR FOR INT(SIGR)dA
C     AR  = FACTOR FOR INT(R*SIGR)dA
C     AS  = FACTOR FOR INT(S*SIGR)dA
C     ARR = FACTOR FOR INT(R*R*SIGR)dA
C     ASS = FACTOR FOR INT(S*S*SIGR)dA
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C	PAPM(5)  = INTEGRAL(TRANS(Pm)*SIGR)
C     ----------------------------------------------------------------
	DIMENSION SIGR(4)
	DIMENSION PAPM(9)

	IF (IPT.EQ.1) THEN
C	  ----------------------------
C	  COMPUTE PAPM - MEMBRANE PART
C	  ----------------------------
	  PAPM(1) = SIGR(1)*A
	  PAPM(2) = SIGR(1)*AS
	  PAPM(3) = SIGR(1)*AR-SIGR(3)*AS
	  PAPM(4) = SIGR(1)*ASS
	  PAPM(5) = SIGR(2)*A
	  PAPM(6) = SIGR(2)*AR
	  PAPM(7) = SIGR(2)*AS-SIGR(3)*AR
	  PAPM(8) = SIGR(2)*ARR
	  PAPM(9) = SIGR(3)*A
	ELSE
C	  ----------------------------
C	  COMPUTE PAPM - MEMBRANE PART
C	  ----------------------------
	  PAPM(1) = PAPM(1)+SIGR(1)*A
	  PAPM(2) = PAPM(2)+SIGR(1)*AS
	  PAPM(3) = PAPM(3)+SIGR(1)*AR-SIGR(3)*AS
	  PAPM(4) = PAPM(4)+SIGR(1)*ASS
	  PAPM(5) = PAPM(5)+SIGR(2)*A
	  PAPM(6) = PAPM(6)+SIGR(2)*AR
	  PAPM(7) = PAPM(7)+SIGR(2)*AS-SIGR(3)*AR
	  PAPM(8) = PAPM(8)+SIGR(2)*ARR
	  PAPM(9) = PAPM(9)+SIGR(3)*A
	ENDIF

      RETURN
      END
C
C=====================================================================
      SUBROUTINE PDMEM(IPT,DR,A,AR,AS,ARS,ARR,ASS,ARRS,ASSR,ARRSS,
	2                ARRR,ASSS,ARRRR,ASSSS,PMPM)
	IMPLICIT REAL*8 (A-H,O-Z)
C     ----------------------------------------------------------------
C     PURPOSE:   COMPUTES FOR INTEGRAL(TRANS(P)*DR*P)dA
C
C	INPUT VARIABLES
C	DR(64)  = NODAL STRESS RESULTANTS
C	A       = FACTOR FOR INT(DR)dA
C     AR      = FACTOR FOR INT(R*DR)dA
C     AS      = FACTOR FOR INT(S*DR)dA
C     ARS     = FACTOR FOR INT(R*S*DR)dA
C     ARR     = FACTOR FOR INT(R*R*DR)dA
C     ASS     = FACTOR FOR INT(S*S*DR)dA
C     ARRS    = FACTOR FOR INT(R*R*S*DR)dA
C     ASSR    = FACTOR FOR INT(S*S*R*DR)dA
C     ARRSS   = FACTOR FOR INT(R*R*S*S*DR)dA
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C	PMPM(9,9)   = INTEGRAL(TRANS(Pm)*DR*Pm) - MEMBRANE PART
C     ----------------------------------------------------------------
	DIMENSION DR(36)
	DIMENSION PMPM(9,9)

	IF (IPT.EQ.1) THEN
C	  ----------------------------
C	  COMPUTE PMPM - MEMBRANE PART
C	  UPPER TRIANGULAR ONLY
C	  ----------------------------
        PMPM(1,1) = A*DR(1)
        PMPM(1,2) = As*DR(1)
        PMPM(1,3) = Ar*DR(1)-As*DR(3)
        PMPM(1,4) = AsS*DR(1)
        PMPM(1,5) = A*DR(2)
        PMPM(1,6) = Ar*DR(2)
        PMPM(1,7) = As*DR(2)-Ar*DR(3)
        PMPM(1,8) = ArR*DR(2)
        PMPM(1,9) = A*DR(3)

        PMPM(2,2) = ASs*DR(1)
        PMPM(2,3) = ARS*DR(1)-ASs*DR(3)
        PMPM(2,4) = ASSs*DR(1)
        PMPM(2,5) = As*DR(2)
        PMPM(2,6) = ARs*DR(2)
        PMPM(2,7) = ASs*DR(2)-ARs*DR(3)
        PMPM(2,8) = ARRs*DR(2)
        PMPM(2,9) = As*DR(3)

        PMPM(3,3) = ARr*DR(1)-2*ARs*DR(3)+ASs*DR(15)
        PMPM(3,4) = ASsR*DR(1)-ASSs*DR(3)
        PMPM(3,5) = Ar*DR(2)-As*DR(9)
        PMPM(3,6) = ARr*DR(2)-ARS*DR(9)
        PMPM(3,7) = ARS*DR(2)-ASs*DR(9)-ARr*DR(3)+Ars*DR(15)
        PMPM(3,8) = ARRr*DR(2)-ARrS*DR(9)
        PMPM(3,9) = Ar*DR(3)-As*DR(15)

        PMPM(4,4) = ASSSs*DR(1)
        PMPM(4,5) = ASs*DR(2)
        PMPM(4,6) = ASsR*DR(2)
        PMPM(4,7) = ASSs*DR(2)-ASsR*DR(3)
        PMPM(4,8) = ARRSS*DR(2)
        PMPM(4,9) = ASs*DR(3)

        PMPM(5,5) = A*DR(8)
        PMPM(5,6) = Ar*DR(8)
        PMPM(5,7) = As*DR(8)-Ar*DR(9)
        PMPM(5,8) = ARr*DR(8)
        PMPM(5,9) = A*DR(9)

        PMPM(6,6) = ARr*DR(8)
        PMPM(6,7) = ARS*DR(8)-ARr*DR(9)
        PMPM(6,8) = ARRr*DR(8)
        PMPM(6,9) = Ar*DR(9)

        PMPM(7,7) = ASs*DR(8)-2*ArS*DR(9)+ARr*DR(15)
        PMPM(7,8) = ARrS*DR(8)-ARRr*DR(9)
        PMPM(7,9) = As*DR(9)-Ar*DR(15)

        PMPM(8,8) = ARRRr*DR(8)
        PMPM(8,9) = ARr*DR(9)

        PMPM(9,9) = A*DR(15)
	ELSE
C	  ----------------------------
C	  COMPUTE PMPM - MEMBRANE PART
C	  UPPER TRIANGULAR ONLY
C	  ----------------------------
        PMPM(1,1) = PMPM(1,1)+A*DR(1)
        PMPM(1,2) = PMPM(1,2)+As*DR(1)
        PMPM(1,3) = PMPM(1,3)+Ar*DR(1)-As*DR(3)
        PMPM(1,4) = PMPM(1,4)+AsS*DR(1)
        PMPM(1,5) = PMPM(1,5)+A*DR(2)
        PMPM(1,6) = PMPM(1,6)+Ar*DR(2)
        PMPM(1,7) = PMPM(1,7)+As*DR(2)-Ar*DR(3)
        PMPM(1,8) = PMPM(1,8)+ArR*DR(2)
        PMPM(1,9) = PMPM(1,9)+A*DR(3)

        PMPM(2,2) = PMPM(2,2)+ASs*DR(1)
        PMPM(2,3) = PMPM(2,3)+ARS*DR(1)-ASs*DR(3)
        PMPM(2,4) = PMPM(2,4)+ASSs*DR(1)
        PMPM(2,5) = PMPM(2,5)+As*DR(2)
        PMPM(2,6) = PMPM(2,6)+ARs*DR(2)
        PMPM(2,7) = PMPM(2,7)+ASs*DR(2)-ARs*DR(3)
        PMPM(2,8) = PMPM(2,8)+ARRs*DR(2)
        PMPM(2,9) = PMPM(2,9)+As*DR(3)

        PMPM(3,3) = PMPM(3,3)+ARr*DR(1)-2*ARs*DR(3)+ASs*DR(15)
        PMPM(3,4) = PMPM(3,4)+ASsR*DR(1)-ASSs*DR(3)
        PMPM(3,5) = PMPM(3,5)+Ar*DR(2)-As*DR(9)
        PMPM(3,6) = PMPM(3,6)+ARr*DR(2)-ARS*DR(9)
        PMPM(3,7) = PMPM(3,7)+ARS*DR(2)-ASs*DR(9)-ARr*DR(3)+Ars*DR(15)
        PMPM(3,8) = PMPM(3,8)+ARRr*DR(2)-ARrS*DR(9)
        PMPM(3,9) = PMPM(3,9)+Ar*DR(3)-As*DR(15)

        PMPM(4,4) = PMPM(4,4)+ASSSs*DR(1)
        PMPM(4,5) = PMPM(4,5)+ASs*DR(2)
        PMPM(4,6) = PMPM(4,6)+ASsR*DR(2)
        PMPM(4,7) = PMPM(4,7)+ASSs*DR(2)-ASsR*DR(3)
        PMPM(4,8) = PMPM(4,8)+ARRSS*DR(2)
        PMPM(4,9) = PMPM(4,9)+ASs*DR(3)

        PMPM(5,5) = PMPM(5,5)+A*DR(8)
        PMPM(5,6) = PMPM(5,6)+Ar*DR(8)
        PMPM(5,7) = PMPM(5,7)+As*DR(8)-Ar*DR(9)
        PMPM(5,8) = PMPM(5,8)+ARr*DR(8)
        PMPM(5,9) = PMPM(5,9)+A*DR(9)

        PMPM(6,6) = PMPM(6,6)+ARr*DR(8)
        PMPM(6,7) = PMPM(6,7)+ARS*DR(8)-ARr*DR(9)
        PMPM(6,8) = PMPM(6,8)+ARRr*DR(8)
        PMPM(6,9) = PMPM(6,9)+Ar*DR(9)

        PMPM(7,7) = PMPM(7,7)+ASs*DR(8)-2*ArS*DR(9)+ARr*DR(15)
        PMPM(7,8) = PMPM(7,8)+ARrS*DR(8)-ARRr*DR(9)
        PMPM(7,9) = PMPM(7,9)+As*DR(9)-Ar*DR(15)

        PMPM(8,8) = PMPM(8,8)+ARRRr*DR(8)
        PMPM(8,9) = PMPM(8,9)+ARr*DR(9)

        PMPM(9,9) = PMPM(9,9)+A*DR(15)
	ENDIF

      RETURN
      END
C
C=====================================================================
      SUBROUTINE COMFACT4(RS,ANTRRR,ANTSSS,ANTRRRR,ANTSSSS)
	IMPLICIT REAL*8 (A-H,O-Z)
C     ----------------------------------------------------------
C     PURPOSE:	COMPUTES FOR INTEGRATION FACTORS FOR INTEGRALS
C				USING ISOPARAMETRIC MAPPING
C
C	INPUT VARIABLES
C	RS(2,4)  = NODAL COORDINATES
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C     ANTRR(4)  = FACTORS FOR INT(R*R*PHI)
C     ANTSS(4)  = FACTORS FOR INT(S*S*PHI)
C     ANTRRS(4) = FACTORS FOR INT(R*R*S*PHI)
C     ANTSSR(4) = FACTORS FOR INT(S*S*R*PHI)
C     ANTRRSS(4) = FACTORS FOR INT(R*R*S*S*PHI)
C     ----------------------------------------------------------
	DIMENSION RS(2,4)
	DIMENSION ANTRRR(4),ANTSSS(4),ANTRRRR(4),ANTSSSS(4)

C	-------------------------
C	COMPUTE FOR FACTOR ANTRRR
C	-------------------------
      s5 = (-RS(2,4)/30+RS(2,2)/30)*RS(1,1)**4+((RS(2,4)/150+RS(2,2)/40+
     #RS(2,3)/600-RS(2,1)/30)*RS(1,2)+(RS(2,4)/600-RS(2,2)/600)*RS(1,3)+
     #(-RS(2,2)/150-RS(2,3)/600+RS(2,1)/30-RS(2,4)/40)*RS(1,4))*RS(1,1)*
     #*3+((RS(2,3)/300+RS(2,4)/200+RS(2,2)/60-RS(2,1)/40)*RS(1,2)**2+((R
     #S(2,3)/1200-RS(2,2)/300+RS(2,4)/400)*RS(1,3)+(RS(2,4)/200-RS(2,2)/
     #200)*RS(1,4))*RS(1,2)+(-RS(2,2)/1200+RS(2,4)/1200)*RS(1,3)**2+(RS(
     #2,4)/300-RS(2,3)/1200-RS(2,2)/400)*RS(1,4)*RS(1,3)+(-RS(2,4)/60-RS
     #(2,2)/200-RS(2,3)/300+RS(2,1)/40)*RS(1,4)**2)*RS(1,1)**2
      s6 = s5+((RS(2,4)/300-RS(2,1)/60+RS(2,2)/120+RS(2,3)/200)*RS(1,2)*
     #*3+((-RS(2,2)/200+RS(2,3)/400+RS(2,4)/400)*RS(1,3)+(RS(2,3)/1200-R
     #S(2,2)/300+RS(2,4)/400)*RS(1,4))*RS(1,2)**2+((-RS(2,2)/400+RS(2,3)
     #/1200+RS(2,4)/600)*RS(1,3)**2+(-RS(2,2)/300+RS(2,4)/300)*RS(1,4)*R
     #S(1,3)+(RS(2,4)/300-RS(2,3)/1200-RS(2,2)/400)*RS(1,4)**2)*RS(1,2)+
     #(-RS(2,2)/1200+RS(2,4)/1200)*RS(1,3)**3+(-RS(2,2)/600-RS(2,3)/1200
     #+RS(2,4)/400)*RS(1,4)*RS(1,3)**2+(RS(2,4)/200-RS(2,2)/400-RS(2,3)/
     #400)*RS(1,4)**2*RS(1,3)+(RS(2,1)/60-RS(2,3)/200-RS(2,4)/120-RS(2,2
     #)/300)*RS(1,4)**3)*RS(1,1)
      s4 = s6+(RS(2,4)/600+RS(2,3)/150-RS(2,1)/120)*RS(1,2)**4+((-RS(2,2
     #)/150+RS(2,3)/200+RS(2,4)/600)*RS(1,3)+(RS(2,3)/1200+RS(2,4)/1200-
     #RS(2,2)/600)*RS(1,4))*RS(1,2)**3
      s5 = s4+((RS(2,3)/300-RS(2,2)/200+RS(2,4)/600)*RS(1,3)**2+(-RS(2,2
     #)/400+RS(2,3)/1200+RS(2,4)/600)*RS(1,4)*RS(1,3)+(-RS(2,2)/1200+RS(
     #2,4)/1200)*RS(1,4)**2)*RS(1,2)**2+((RS(2,4)/600-RS(2,2)/300+RS(2,3
     #)/600)*RS(1,3)**3+(RS(2,4)/400-RS(2,2)/400)*RS(1,4)*RS(1,3)**2+(-R
     #S(2,2)/600-RS(2,3)/1200+RS(2,4)/400)*RS(1,4)**2*RS(1,3)+(RS(2,4)/6
     #00-RS(2,2)/1200-RS(2,3)/1200)*RS(1,4)**3)*RS(1,2)+(RS(2,4)/600-RS(
     #2,2)/600)*RS(1,3)**4
	ANTRRR(1) =
     #s5+(-RS(2,3)/600-RS(2,2)/600+RS(2,4)/300)*RS(1,4)*RS(1,3)**3+
     #(RS(2,4)/200-RS(2,2)/600-RS(2,3)/300)*RS(1,4)**2*RS(1,3)**2+(-RS(2
     #,3)/200-RS(2,2)/600+RS(2,4)/150)*RS(1,4)**3*RS(1,3)+(-RS(2,2)/600-
     #RS(2,3)/150+RS(2,1)/120)*RS(1,4)**4

      s6 = (-RS(2,3)/600+RS(2,2)/120-RS(2,4)/150)*RS(1,1)**4+((-RS(2,3)/
     #300-RS(2,4)/200+RS(2,2)/60-RS(2,1)/120)*RS(1,2)+(RS(2,1)/600-RS(2,
     #4)/1200-RS(2,3)/1200)*RS(1,3)+(-RS(2,4)/200+RS(2,1)/150-RS(2,3)/60
     #0)*RS(1,4))*RS(1,1)**3+((-RS(2,1)/60-RS(2,3)/200+RS(2,2)/40-RS(2,4
     #)/300)*RS(1,2)**2+((RS(2,1)/300-RS(2,3)/400-RS(2,4)/1200)*RS(1,3)+
     #(RS(2,1)/200-RS(2,3)/400-RS(2,4)/400)*RS(1,4))*RS(1,2)+(RS(2,1)/12
     #00-RS(2,3)/1200)*RS(1,3)**2+(-RS(2,3)/600-RS(2,4)/1200+RS(2,1)/400
     #)*RS(1,4)*RS(1,3)+(RS(2,1)/200-RS(2,4)/300-RS(2,3)/600)*RS(1,4)**2
     #)*RS(1,1)**2
      s7 = s6+((RS(2,2)/30-RS(2,1)/40-RS(2,4)/600-RS(2,3)/150)*RS(1,2)**
     #3+((-RS(2,3)/200+RS(2,1)/200)*RS(1,3)+(RS(2,1)/300-RS(2,3)/400-RS(
     #2,4)/1200)*RS(1,4))*RS(1,2)**2+((-RS(2,3)/300+RS(2,4)/1200+RS(2,1)
     #/400)*RS(1,3)**2+(-RS(2,3)/300+RS(2,1)/300)*RS(1,4)*RS(1,3)+(-RS(2
     #,3)/600-RS(2,4)/1200+RS(2,1)/400)*RS(1,4)**2)*RS(1,2)+(-RS(2,3)/60
     #0+RS(2,4)/1200+RS(2,1)/1200)*RS(1,3)**3+(RS(2,1)/600-RS(2,3)/400+R
     #S(2,4)/1200)*RS(1,4)*RS(1,3)**2+(RS(2,1)/400-RS(2,3)/400)*RS(1,4)*
     #*2*RS(1,3)+(RS(2,1)/300-RS(2,3)/600-RS(2,4)/600)*RS(1,4)**3)*RS(1,
     #1)
      s5 = s7+(-RS(2,1)/30+RS(2,3)/30)*RS(1,2)**4+((-RS(2,2)/30+RS(2,4)/
     #600+RS(2,3)/40+RS(2,1)/150)*RS(1,3)+(RS(2,1)/600-RS(2,3)/600)*RS(1
     #,4))*RS(1,2)**3
      s6 = s5+((RS(2,1)/200+RS(2,3)/60-RS(2,2)/40+RS(2,4)/300)*RS(1,3)**
     #2+(-RS(2,3)/300+RS(2,4)/1200+RS(2,1)/400)*RS(1,4)*RS(1,3)+(RS(2,1)
     #/1200-RS(2,3)/1200)*RS(1,4)**2)*RS(1,2)**2+((-RS(2,2)/60+RS(2,3)/1
     #20+RS(2,1)/300+RS(2,4)/200)*RS(1,3)**3+(RS(2,1)/400+RS(2,4)/400-RS
     #(2,3)/200)*RS(1,4)*RS(1,3)**2+(RS(2,1)/600-RS(2,3)/400+RS(2,4)/120
     #0)*RS(1,4)**2*RS(1,3)+(RS(2,1)/1200-RS(2,3)/1200)*RS(1,4)**3)*RS(1
     #,2)+(RS(2,1)/600-RS(2,2)/120+RS(2,4)/150)*RS(1,3)**4
	ANTRRR(2) =
     #s6+(RS(2,4)/200-RS(2,3)/150+RS(2,1)/600)*RS(1,4)*RS(1,3)**3+(
     #RS(2,4)/300+RS(2,1)/600-RS(2,3)/200)*RS(1,4)**2*RS(1,3)**2+(RS(2,4
     #)/600-RS(2,3)/300+RS(2,1)/600)*RS(1,4)**3*RS(1,3)+(RS(2,1)/600-RS(
     #2,3)/600)*RS(1,4)**4

      s7 = (-RS(2,4)/600+RS(2,2)/600)*RS(1,1)**4+((-RS(2,4)/600-RS(2,1)/
     #600+RS(2,2)/300)*RS(1,2)+(RS(2,2)/1200-RS(2,4)/1200)*RS(1,3)+(-RS(
     #2,4)/300+RS(2,1)/600+RS(2,2)/600)*RS(1,4))*RS(1,1)**3+((RS(2,2)/20
     #0-RS(2,1)/300-RS(2,4)/600)*RS(1,2)**2+((RS(2,2)/400-RS(2,1)/1200-R
     #S(2,4)/600)*RS(1,3)+(RS(2,2)/400-RS(2,4)/400)*RS(1,4))*RS(1,2)+(RS
     #(2,2)/1200-RS(2,4)/1200)*RS(1,3)**2+(RS(2,1)/1200+RS(2,2)/600-RS(2
     #,4)/400)*RS(1,4)*RS(1,3)+(RS(2,2)/600+RS(2,1)/300-RS(2,4)/200)*RS(
     #1,4)**2)*RS(1,1)**2
      s8 = s7+((-RS(2,1)/200-RS(2,4)/600+RS(2,2)/150)*RS(1,2)**3+((-RS(2
     #,4)/400+RS(2,2)/200-RS(2,1)/400)*RS(1,3)+(RS(2,2)/400-RS(2,1)/1200
     #-RS(2,4)/600)*RS(1,4))*RS(1,2)**2+((-RS(2,4)/400+RS(2,2)/300-RS(2,
     #1)/1200)*RS(1,3)**2+(RS(2,2)/300-RS(2,4)/300)*RS(1,4)*RS(1,3)+(RS(
     #2,1)/1200+RS(2,2)/600-RS(2,4)/400)*RS(1,4)**2)*RS(1,2)+(-RS(2,4)/6
     #00+RS(2,2)/600)*RS(1,3)**3+(RS(2,2)/400+RS(2,1)/1200-RS(2,4)/300)*
     #RS(1,4)*RS(1,3)**2+(-RS(2,4)/200+RS(2,2)/400+RS(2,1)/400)*RS(1,4)*
     #*2*RS(1,3)+(RS(2,1)/200-RS(2,4)/150+RS(2,2)/600)*RS(1,4)**3)*RS(1,
     #1)
      s6 = s8+(-RS(2,1)/150-RS(2,4)/600+RS(2,3)/120)*RS(1,2)**4+((-RS(2,
     #2)/120+RS(2,3)/60-RS(2,1)/200-RS(2,4)/300)*RS(1,3)+(-RS(2,4)/1200+
     #RS(2,2)/600-RS(2,1)/1200)*RS(1,4))*RS(1,2)**3
      s7 = s6+((-RS(2,2)/60-RS(2,1)/300-RS(2,4)/200+RS(2,3)/40)*RS(1,3)*
     #*2+(-RS(2,4)/400+RS(2,2)/300-RS(2,1)/1200)*RS(1,4)*RS(1,3)+(RS(2,2
     #)/1200-RS(2,4)/1200)*RS(1,4)**2)*RS(1,2)**2+((RS(2,3)/30-RS(2,4)/1
     #50-RS(2,2)/40-RS(2,1)/600)*RS(1,3)**3+(RS(2,2)/200-RS(2,4)/200)*RS
     #(1,4)*RS(1,3)**2+(RS(2,2)/400+RS(2,1)/1200-RS(2,4)/300)*RS(1,4)**2
     #*RS(1,3)+(-RS(2,4)/600+RS(2,1)/1200+RS(2,2)/1200)*RS(1,4)**3)*RS(1
     #,2)+(RS(2,4)/30-RS(2,2)/30)*RS(1,3)**4
	ANTRRR(3) =
     #s7+(-RS(2,3)/30+RS(2,4)/40+RS(2,2)/150+RS(2,1)/600)*RS(1,4)*R
     #S(1,3)**3+(RS(2,4)/60+RS(2,1)/300+RS(2,2)/200-RS(2,3)/40)*RS(1,4)*
     #*2*RS(1,3)**2+(RS(2,2)/300+RS(2,4)/120+RS(2,1)/200-RS(2,3)/60)*RS(
     #1,4)**3*RS(1,3)+(RS(2,2)/600-RS(2,3)/120+RS(2,1)/150)*RS(1,4)**4

      s8 = (RS(2,3)/600+RS(2,2)/150-RS(2,4)/120)*RS(1,1)**4+((RS(2,2)/20
     #0+RS(2,3)/600-RS(2,1)/150)*RS(1,2)+(RS(2,3)/1200-RS(2,1)/600+RS(2,
     #2)/1200)*RS(1,3)+(RS(2,3)/300-RS(2,4)/60+RS(2,1)/120+RS(2,2)/200)*
     #RS(1,4))*RS(1,1)**3+((RS(2,3)/600-RS(2,1)/200+RS(2,2)/300)*RS(1,2)
     #**2+((RS(2,3)/600-RS(2,1)/400+RS(2,2)/1200)*RS(1,3)+(RS(2,3)/400-R
     #S(2,1)/200+RS(2,2)/400)*RS(1,4))*RS(1,2)+(RS(2,3)/1200-RS(2,1)/120
     #0)*RS(1,3)**2+(RS(2,2)/1200+RS(2,3)/400-RS(2,1)/300)*RS(1,4)*RS(1,
     #3)+(RS(2,1)/60+RS(2,3)/200+RS(2,2)/300-RS(2,4)/40)*RS(1,4)**2)*RS(
     #1,1)**2
      s9 = s8+((RS(2,2)/600-RS(2,1)/300+RS(2,3)/600)*RS(1,2)**3+((RS(2,3
     #)/400-RS(2,1)/400)*RS(1,3)+(RS(2,3)/600-RS(2,1)/400+RS(2,2)/1200)*
     #RS(1,4))*RS(1,2)**2+((-RS(2,1)/600-RS(2,2)/1200+RS(2,3)/400)*RS(1,
     #3)**2+(RS(2,3)/300-RS(2,1)/300)*RS(1,4)*RS(1,3)+(RS(2,2)/1200+RS(2
     #,3)/400-RS(2,1)/300)*RS(1,4)**2)*RS(1,2)+(-RS(2,2)/1200-RS(2,1)/12
     #00+RS(2,3)/600)*RS(1,3)**3+(-RS(2,1)/400-RS(2,2)/1200+RS(2,3)/300)
     #*RS(1,4)*RS(1,3)**2+(-RS(2,1)/200+RS(2,3)/200)*RS(1,4)**2*RS(1,3)+
     #(RS(2,1)/40+RS(2,3)/150-RS(2,4)/30+RS(2,2)/600)*RS(1,4)**3)*RS(1,1
     #)
      s7 = s9+(-RS(2,1)/600+RS(2,3)/600)*RS(1,2)**4+((RS(2,3)/300-RS(2,2
     #)/600-RS(2,1)/600)*RS(1,3)+(RS(2,3)/1200-RS(2,1)/1200)*RS(1,4))*RS
     #(1,2)**3
      s8 = s7+((-RS(2,2)/300+RS(2,3)/200-RS(2,1)/600)*RS(1,3)**2+(-RS(2,
     #1)/600-RS(2,2)/1200+RS(2,3)/400)*RS(1,4)*RS(1,3)+(RS(2,3)/1200-RS(
     #2,1)/1200)*RS(1,4)**2)*RS(1,2)**2+((-RS(2,2)/200+RS(2,3)/150-RS(2,
     #1)/600)*RS(1,3)**3+(-RS(2,1)/400+RS(2,3)/200-RS(2,2)/400)*RS(1,4)*
     #RS(1,3)**2+(-RS(2,1)/400-RS(2,2)/1200+RS(2,3)/300)*RS(1,4)**2*RS(1
     #,3)+(-RS(2,1)/600+RS(2,3)/600)*RS(1,4)**3)*RS(1,2)+(RS(2,4)/120-RS
     #(2,1)/600-RS(2,2)/150)*RS(1,3)**4
	ANTRRR(4) =
     #s8+(-RS(2,3)/120-RS(2,2)/200+RS(2,4)/60-RS(2,1)/300)*RS(1,4)*
     #RS(1,3)**3+(-RS(2,2)/300+RS(2,4)/40-RS(2,1)/200-RS(2,3)/60)*RS(1,4
     #)**2*RS(1,3)**2+(RS(2,4)/30-RS(2,3)/40-RS(2,1)/150-RS(2,2)/600)*RS
     #(1,4)**3*RS(1,3)+(-RS(2,3)/30+RS(2,1)/30)*RS(1,4)**4

C	-------------------------
C	COMPUTE FOR FACTOR ANTSSS
C	-------------------------
      s5 = (-RS(2,1)**2*RS(2,4)**2/40-RS(2,1)**3*RS(2,4)/30+RS(2,1)**3*R
     #S(2,2)/30+RS(2,2)**4/120+RS(2,1)*RS(2,2)**3/60+RS(2,1)**2*RS(2,2)*
     #*2/40-RS(2,1)*RS(2,4)**3/60-RS(2,4)**4/120)*RS(1,1)
      s8 = RS(2,2)*RS(2,3)**2*RS(2,4)/400+RS(2,1)*RS(2,2)*RS(2,3)**2/400
     #+RS(2,1)**2*RS(2,2)*RS(2,4)/200+RS(2,1)*RS(2,3)**2*RS(2,4)/600+RS(
     #2,1)**2*RS(2,3)*RS(2,4)/400+RS(2,1)**2*RS(2,2)*RS(2,3)/300+RS(2,1)
     #*RS(2,2)**2*RS(2,3)/200+RS(2,1)*RS(2,3)*RS(2,4)**2/400+RS(2,2)*RS(
     #2,4)**2*RS(2,1)/400+RS(2,2)**2*RS(2,4)*RS(2,3)/400+RS(2,2)**2*RS(2
     #,4)*RS(2,1)/300+RS(2,2)*RS(2,4)**2*RS(2,3)/600+RS(2,1)*RS(2,2)*RS(
     #2,3)*RS(2,4)/300+RS(2,3)**4/600-RS(2,1)**4/30+RS(2,1)*RS(2,3)**3/1
     #200+RS(2,3)**2*RS(2,4)**2/600
      s7 = s8+RS(2,1)**2*RS(2,3)**2/1200+RS(2,2)*RS(2,4)**3/1200+RS(2,2)
     #**3*RS(2,4)/600+RS(2,2)**2*RS(2,4)**2/1200+RS(2,2)*RS(2,3)**3/300+
     #RS(2,2)**3*RS(2,3)/150+RS(2,4)**3*RS(2,3)/600+RS(2,3)**3*RS(2,4)/6
     #00+RS(2,2)**2*RS(2,3)**2/200+RS(2,1)**3*RS(2,3)/600+RS(2,1)**2*RS(
     #2,4)**2/200+RS(2,1)**3*RS(2,4)/150-RS(2,1)**3*RS(2,2)/40+RS(2,4)**
     #4/600+RS(2,1)*RS(2,4)**3/300-RS(2,1)**2*RS(2,2)**2/60-RS(2,1)*RS(2
     #,2)**3/120
      s8 = RS(1,2)
      s6 = s7*s8
      s4 = s5+s6
      s5 = s4
      s9 = -RS(2,1)*RS(2,2)*RS(2,3)**2/1200+RS(2,1)*RS(2,3)**2*RS(2,4)/1
     #200+RS(2,1)**2*RS(2,3)*RS(2,4)/1200-RS(2,1)**2*RS(2,2)*RS(2,3)/120
     #0-RS(2,1)*RS(2,2)**2*RS(2,3)/400+RS(2,1)*RS(2,3)*RS(2,4)**2/400+RS
     #(2,2)*RS(2,4)**2*RS(2,1)/1200-RS(2,2)**2*RS(2,4)*RS(2,3)/1200-RS(2
     #,2)**2*RS(2,4)*RS(2,1)/1200+RS(2,2)*RS(2,4)**2*RS(2,3)/1200+RS(2,3
     #)**2*RS(2,4)**2/300+RS(2,2)*RS(2,4)**3/1200-RS(2,2)**3*RS(2,4)/120
     #0
      s8 = s9-RS(2,2)*RS(2,3)**3/600-RS(2,2)**3*RS(2,3)/200+RS(2,4)**3*R
     #S(2,3)/200+RS(2,3)**3*RS(2,4)/600-RS(2,2)**2*RS(2,3)**2/300+RS(2,1
     #)**2*RS(2,4)**2/300+RS(2,1)**3*RS(2,4)/600-RS(2,1)**3*RS(2,2)/600-
     #RS(2,2)**4/150+RS(2,4)**4/150+RS(2,1)*RS(2,4)**3/200-RS(2,1)**2*RS
     #(2,2)**2/300-RS(2,1)*RS(2,2)**3/200
      s9 = RS(1,3)
      s7 = s8*s9
      s10 = -RS(2,2)*RS(2,3)**2*RS(2,4)/400-RS(2,1)*RS(2,2)*RS(2,3)**2/6
     #00-RS(2,1)**2*RS(2,2)*RS(2,4)/200-RS(2,1)*RS(2,3)**2*RS(2,4)/400-R
     #S(2,1)**2*RS(2,3)*RS(2,4)/300-RS(2,1)**2*RS(2,2)*RS(2,3)/400-RS(2,
     #1)*RS(2,2)**2*RS(2,3)/400-RS(2,1)*RS(2,3)*RS(2,4)**2/200-RS(2,2)*R
     #S(2,4)**2*RS(2,1)/300-RS(2,2)**2*RS(2,4)*RS(2,3)/600-RS(2,2)**2*RS
     #(2,4)*RS(2,1)/400-RS(2,2)*RS(2,4)**2*RS(2,3)/400-RS(2,1)*RS(2,2)*R
     #S(2,3)*RS(2,4)/300-RS(2,3)**4/600+RS(2,1)**4/30-RS(2,1)*RS(2,3)**3
     #/1200-RS(2,3)**2*RS(2,4)**2/200
      s9 = s10-RS(2,1)**2*RS(2,3)**2/1200-RS(2,2)*RS(2,4)**3/600-RS(2,2)
     #**3*RS(2,4)/1200-RS(2,2)**2*RS(2,4)**2/1200-RS(2,2)*RS(2,3)**3/600
     #-RS(2,2)**3*RS(2,3)/600-RS(2,4)**3*RS(2,3)/150-RS(2,3)**3*RS(2,4)/
     #300-RS(2,2)**2*RS(2,3)**2/600-RS(2,1)**3*RS(2,3)/600+RS(2,1)**2*RS
     #(2,4)**2/60+RS(2,1)**3*RS(2,4)/40-RS(2,1)**3*RS(2,2)/150-RS(2,2)**
     #4/600+RS(2,1)*RS(2,4)**3/120-RS(2,1)**2*RS(2,2)**2/200-RS(2,1)*RS(
     #2,2)**3/300
      s10 = RS(1,4)
      s8 = s9*s10
      s6 = s7+s8
      ANTSSS(1) = s5+s6

      s8 = -RS(2,2)*RS(2,3)**2*RS(2,4)/400-RS(2,1)*RS(2,2)*RS(2,3)**2/40
     #0-RS(2,1)**2*RS(2,2)*RS(2,4)/200-RS(2,1)*RS(2,3)**2*RS(2,4)/600-RS
     #(2,1)**2*RS(2,3)*RS(2,4)/400-RS(2,1)**2*RS(2,2)*RS(2,3)/300-RS(2,1
     #)*RS(2,2)**2*RS(2,3)/200-RS(2,1)*RS(2,3)*RS(2,4)**2/400-RS(2,2)*RS
     #(2,4)**2*RS(2,1)/400-RS(2,2)**2*RS(2,4)*RS(2,3)/400-RS(2,2)**2*RS(
     #2,4)*RS(2,1)/300-RS(2,2)*RS(2,4)**2*RS(2,3)/600-RS(2,1)*RS(2,2)*RS
     #(2,3)*RS(2,4)/300-RS(2,3)**4/600-RS(2,1)*RS(2,3)**3/1200-RS(2,3)**
     #2*RS(2,4)**2/600-RS(2,1)**2*RS(2,3)**2/1200
      s7 = s8-RS(2,2)*RS(2,4)**3/1200-RS(2,2)**3*RS(2,4)/600-RS(2,2)**2*
     #RS(2,4)**2/1200-RS(2,2)*RS(2,3)**3/300-RS(2,2)**3*RS(2,3)/150-RS(2
     #,4)**3*RS(2,3)/600-RS(2,3)**3*RS(2,4)/600-RS(2,2)**2*RS(2,3)**2/20
     #0-RS(2,1)**3*RS(2,3)/600-RS(2,1)**2*RS(2,4)**2/200-RS(2,1)**3*RS(2
     #,4)/150+RS(2,1)**3*RS(2,2)/120+RS(2,2)**4/30-RS(2,4)**4/600-RS(2,1
     #)*RS(2,4)**3/300+RS(2,1)**2*RS(2,2)**2/60+RS(2,1)*RS(2,2)**3/40
      s8 = RS(1,1)
      s6 = s7*s8
      s7 = (RS(2,2)**2*RS(2,3)**2/40+RS(2,2)**3*RS(2,3)/30+RS(2,3)**4/12
     #0-RS(2,1)*RS(2,2)**3/30+RS(2,2)*RS(2,3)**3/60-RS(2,1)**2*RS(2,2)**
     #2/40-RS(2,1)**4/120-RS(2,1)**3*RS(2,2)/60)*RS(1,2)
      s5 = s6+s7
      s6 = s5
      s10 = RS(2,2)*RS(2,3)**2*RS(2,4)/200+RS(2,1)*RS(2,2)*RS(2,3)**2/30
     #0+RS(2,1)**2*RS(2,2)*RS(2,4)/400+RS(2,1)*RS(2,3)**2*RS(2,4)/400+RS
     #(2,1)**2*RS(2,3)*RS(2,4)/600+RS(2,1)**2*RS(2,2)*RS(2,3)/400+RS(2,1
     #)*RS(2,2)**2*RS(2,3)/200+RS(2,1)*RS(2,3)*RS(2,4)**2/400+RS(2,2)*RS
     #(2,4)**2*RS(2,1)/600+RS(2,2)**2*RS(2,4)*RS(2,3)/300+RS(2,2)**2*RS(
     #2,4)*RS(2,1)/400+RS(2,2)*RS(2,4)**2*RS(2,3)/400+RS(2,1)*RS(2,2)*RS
     #(2,3)*RS(2,4)/300+RS(2,1)**4/600+RS(2,1)*RS(2,3)**3/600+RS(2,3)**2
     #*RS(2,4)**2/200+RS(2,1)**2*RS(2,3)**2/1200
      s9 = s10+RS(2,2)*RS(2,4)**3/1200+RS(2,2)**3*RS(2,4)/600+RS(2,2)**2
     #*RS(2,4)**2/1200-RS(2,2)*RS(2,3)**3/120-RS(2,2)**3*RS(2,3)/40+RS(2
     #,4)**3*RS(2,3)/300+RS(2,3)**3*RS(2,4)/150-RS(2,2)**2*RS(2,3)**2/60
     #+RS(2,1)**3*RS(2,3)/1200+RS(2,1)**2*RS(2,4)**2/600+RS(2,1)**3*RS(2
     #,4)/600+RS(2,1)**3*RS(2,2)/300-RS(2,2)**4/30+RS(2,4)**4/600+RS(2,1
     #)*RS(2,4)**3/600+RS(2,1)**2*RS(2,2)**2/200+RS(2,1)*RS(2,2)**3/150
      s10 = RS(1,3)
      s8 = s9*s10
      s11 = -RS(2,2)*RS(2,3)**2*RS(2,4)/400-RS(2,1)*RS(2,2)*RS(2,3)**2/1
     #200+RS(2,1)**2*RS(2,2)*RS(2,4)/400-RS(2,1)*RS(2,3)**2*RS(2,4)/1200
     #+RS(2,1)**2*RS(2,3)*RS(2,4)/1200+RS(2,1)**2*RS(2,2)*RS(2,3)/1200+R
     #S(2,2)*RS(2,4)**2*RS(2,1)/1200-RS(2,2)**2*RS(2,4)*RS(2,3)/1200+RS(
     #2,2)**2*RS(2,4)*RS(2,1)/1200-RS(2,2)*RS(2,4)**2*RS(2,3)/1200-RS(2,
     #3)**4/150+RS(2,1)**4/150-RS(2,1)*RS(2,3)**3/1200
      s10 = s11-RS(2,3)**2*RS(2,4)**2/300-RS(2,2)*RS(2,3)**3/200-RS(2,2)
     #**3*RS(2,3)/600-RS(2,4)**3*RS(2,3)/600-RS(2,3)**3*RS(2,4)/200-RS(2
     #,2)**2*RS(2,3)**2/300+RS(2,1)**3*RS(2,3)/1200+RS(2,1)**2*RS(2,4)**
     #2/300+RS(2,1)**3*RS(2,4)/200+RS(2,1)**3*RS(2,2)/200+RS(2,1)*RS(2,4
     #)**3/600+RS(2,1)**2*RS(2,2)**2/300+RS(2,1)*RS(2,2)**3/600
      s11 = RS(1,4)
      s9 = s10*s11
      s7 = s8+s9
      ANTSSS(2) = s6+s7

      s9 = RS(2,1)*RS(2,2)*RS(2,3)**2/1200-RS(2,1)*RS(2,3)**2*RS(2,4)/12
     #00-RS(2,1)**2*RS(2,3)*RS(2,4)/1200+RS(2,1)**2*RS(2,2)*RS(2,3)/1200
     #+RS(2,1)*RS(2,2)**2*RS(2,3)/400-RS(2,1)*RS(2,3)*RS(2,4)**2/400-RS(
     #2,2)*RS(2,4)**2*RS(2,1)/1200+RS(2,2)**2*RS(2,4)*RS(2,3)/1200+RS(2,
     #2)**2*RS(2,4)*RS(2,1)/1200-RS(2,2)*RS(2,4)**2*RS(2,3)/1200-RS(2,3)
     #**2*RS(2,4)**2/300-RS(2,2)*RS(2,4)**3/1200+RS(2,2)**3*RS(2,4)/1200
      s8 = s9+RS(2,2)*RS(2,3)**3/600+RS(2,2)**3*RS(2,3)/200-RS(2,4)**3*R
     #S(2,3)/200-RS(2,3)**3*RS(2,4)/600+RS(2,2)**2*RS(2,3)**2/300-RS(2,1
     #)**2*RS(2,4)**2/300-RS(2,1)**3*RS(2,4)/600+RS(2,1)**3*RS(2,2)/600+
     #RS(2,2)**4/150-RS(2,4)**4/150-RS(2,1)*RS(2,4)**3/200+RS(2,1)**2*RS
     #(2,2)**2/300+RS(2,1)*RS(2,2)**3/200
      s9 = RS(1,1)
      s7 = s8*s9
      s10 = -RS(2,2)*RS(2,3)**2*RS(2,4)/200-RS(2,1)*RS(2,2)*RS(2,3)**2/3
     #00-RS(2,1)**2*RS(2,2)*RS(2,4)/400-RS(2,1)*RS(2,3)**2*RS(2,4)/400-R
     #S(2,1)**2*RS(2,3)*RS(2,4)/600-RS(2,1)**2*RS(2,2)*RS(2,3)/400-RS(2,
     #1)*RS(2,2)**2*RS(2,3)/200-RS(2,1)*RS(2,3)*RS(2,4)**2/400-RS(2,2)*R
     #S(2,4)**2*RS(2,1)/600-RS(2,2)**2*RS(2,4)*RS(2,3)/300-RS(2,2)**2*RS
     #(2,4)*RS(2,1)/400-RS(2,2)*RS(2,4)**2*RS(2,3)/400-RS(2,1)*RS(2,2)*R
     #S(2,3)*RS(2,4)/300+RS(2,3)**4/30-RS(2,1)**4/600-RS(2,1)*RS(2,3)**3
     #/600-RS(2,3)**2*RS(2,4)**2/200
      s9 = s10-RS(2,1)**2*RS(2,3)**2/1200-RS(2,2)*RS(2,4)**3/1200-RS(2,2
     #)**3*RS(2,4)/600-RS(2,2)**2*RS(2,4)**2/1200+RS(2,2)*RS(2,3)**3/40+
     #RS(2,2)**3*RS(2,3)/120-RS(2,4)**3*RS(2,3)/300-RS(2,3)**3*RS(2,4)/1
     #50+RS(2,2)**2*RS(2,3)**2/60-RS(2,1)**3*RS(2,3)/1200-RS(2,1)**2*RS(
     #2,4)**2/600-RS(2,1)**3*RS(2,4)/600-RS(2,1)**3*RS(2,2)/300-RS(2,4)*
     #*4/600-RS(2,1)*RS(2,4)**3/600-RS(2,1)**2*RS(2,2)**2/200-RS(2,1)*RS
     #(2,2)**3/150
      s10 = RS(1,2)
      s8 = s9*s10
      s6 = s7+s8
      s7 = s6
      s9 = (RS(2,3)**2*RS(2,4)**2/40-RS(2,2)**3*RS(2,3)/60+RS(2,4)**3*RS
     #(2,3)/60-RS(2,2)*RS(2,3)**3/30-RS(2,2)**2*RS(2,3)**2/40-RS(2,2)**4
     #/120+RS(2,3)**3*RS(2,4)/30+RS(2,4)**4/120)*RS(1,3)
      s12 = RS(2,2)*RS(2,3)**2*RS(2,4)/200+RS(2,1)*RS(2,2)*RS(2,3)**2/40
     #0+RS(2,1)**2*RS(2,2)*RS(2,4)/400+RS(2,1)*RS(2,3)**2*RS(2,4)/300+RS
     #(2,1)**2*RS(2,3)*RS(2,4)/400+RS(2,1)**2*RS(2,2)*RS(2,3)/600+RS(2,1
     #)*RS(2,2)**2*RS(2,3)/400+RS(2,1)*RS(2,3)*RS(2,4)**2/200+RS(2,2)*RS
     #(2,4)**2*RS(2,1)/400+RS(2,2)**2*RS(2,4)*RS(2,3)/400+RS(2,2)**2*RS(
     #2,4)*RS(2,1)/600+RS(2,2)*RS(2,4)**2*RS(2,3)/300+RS(2,1)*RS(2,2)*RS
     #(2,3)*RS(2,4)/300-RS(2,3)**4/30+RS(2,1)**4/600+RS(2,1)*RS(2,3)**3/
     #600-RS(2,3)**2*RS(2,4)**2/60
      s11 = s12+RS(2,1)**2*RS(2,3)**2/1200+RS(2,2)*RS(2,4)**3/600+RS(2,2
     #)**3*RS(2,4)/1200+RS(2,2)**2*RS(2,4)**2/1200+RS(2,2)*RS(2,3)**3/15
     #0+RS(2,2)**3*RS(2,3)/300-RS(2,4)**3*RS(2,3)/120-RS(2,3)**3*RS(2,4)
     #/40+RS(2,2)**2*RS(2,3)**2/200+RS(2,1)**3*RS(2,3)/1200+RS(2,1)**2*R
     #S(2,4)**2/200+RS(2,1)**3*RS(2,4)/300+RS(2,1)**3*RS(2,2)/600+RS(2,2
     #)**4/600+RS(2,1)*RS(2,4)**3/150+RS(2,1)**2*RS(2,2)**2/600+RS(2,1)*
     #RS(2,2)**3/600
      s12 = RS(1,4)
      s10 = s11*s12
      s8 = s9+s10
      ANTSSS(3) = s7+s8

      s10 = RS(2,2)*RS(2,3)**2*RS(2,4)/400+RS(2,1)*RS(2,2)*RS(2,3)**2/60
     #0+RS(2,1)**2*RS(2,2)*RS(2,4)/200+RS(2,1)*RS(2,3)**2*RS(2,4)/400+RS
     #(2,1)**2*RS(2,3)*RS(2,4)/300+RS(2,1)**2*RS(2,2)*RS(2,3)/400+RS(2,1
     #)*RS(2,2)**2*RS(2,3)/400+RS(2,1)*RS(2,3)*RS(2,4)**2/200+RS(2,2)*RS
     #(2,4)**2*RS(2,1)/300+RS(2,2)**2*RS(2,4)*RS(2,3)/600+RS(2,2)**2*RS(
     #2,4)*RS(2,1)/400+RS(2,2)*RS(2,4)**2*RS(2,3)/400+RS(2,1)*RS(2,2)*RS
     #(2,3)*RS(2,4)/300+RS(2,3)**4/600+RS(2,1)*RS(2,3)**3/1200+RS(2,3)**
     #2*RS(2,4)**2/200+RS(2,1)**2*RS(2,3)**2/1200
      s9 = s10+RS(2,2)*RS(2,4)**3/600+RS(2,2)**3*RS(2,4)/1200+RS(2,2)**2
     #*RS(2,4)**2/1200+RS(2,2)*RS(2,3)**3/600+RS(2,2)**3*RS(2,3)/600+RS(
     #2,4)**3*RS(2,3)/150+RS(2,3)**3*RS(2,4)/300+RS(2,2)**2*RS(2,3)**2/6
     #00+RS(2,1)**3*RS(2,3)/600-RS(2,1)**2*RS(2,4)**2/60-RS(2,1)**3*RS(2
     #,4)/120+RS(2,1)**3*RS(2,2)/150+RS(2,2)**4/600-RS(2,4)**4/30-RS(2,1
     #)*RS(2,4)**3/40+RS(2,1)**2*RS(2,2)**2/200+RS(2,1)*RS(2,2)**3/300
      s10 = RS(1,1)
      s8 = s9*s10
      s11 = RS(2,2)*RS(2,3)**2*RS(2,4)/400+RS(2,1)*RS(2,2)*RS(2,3)**2/12
     #00-RS(2,1)**2*RS(2,2)*RS(2,4)/400+RS(2,1)*RS(2,3)**2*RS(2,4)/1200-
     #RS(2,1)**2*RS(2,3)*RS(2,4)/1200-RS(2,1)**2*RS(2,2)*RS(2,3)/1200-RS
     #(2,2)*RS(2,4)**2*RS(2,1)/1200+RS(2,2)**2*RS(2,4)*RS(2,3)/1200-RS(2
     #,2)**2*RS(2,4)*RS(2,1)/1200+RS(2,2)*RS(2,4)**2*RS(2,3)/1200+RS(2,3
     #)**4/150-RS(2,1)**4/150+RS(2,1)*RS(2,3)**3/1200
      s10 = s11+RS(2,3)**2*RS(2,4)**2/300+RS(2,2)*RS(2,3)**3/200+RS(2,2)
     #**3*RS(2,3)/600+RS(2,4)**3*RS(2,3)/600+RS(2,3)**3*RS(2,4)/200+RS(2
     #,2)**2*RS(2,3)**2/300-RS(2,1)**3*RS(2,3)/1200-RS(2,1)**2*RS(2,4)**
     #2/300-RS(2,1)**3*RS(2,4)/200-RS(2,1)**3*RS(2,2)/200-RS(2,1)*RS(2,4
     #)**3/600-RS(2,1)**2*RS(2,2)**2/300-RS(2,1)*RS(2,2)**3/600
      s11 = RS(1,2)
      s9 = s10*s11
      s7 = s8+s9
      s8 = s7
      s12 = -RS(2,2)*RS(2,3)**2*RS(2,4)/200-RS(2,1)*RS(2,2)*RS(2,3)**2/4
     #00-RS(2,1)**2*RS(2,2)*RS(2,4)/400-RS(2,1)*RS(2,3)**2*RS(2,4)/300-R
     #S(2,1)**2*RS(2,3)*RS(2,4)/400-RS(2,1)**2*RS(2,2)*RS(2,3)/600-RS(2,
     #1)*RS(2,2)**2*RS(2,3)/400-RS(2,1)*RS(2,3)*RS(2,4)**2/200-RS(2,2)*R
     #S(2,4)**2*RS(2,1)/400-RS(2,2)**2*RS(2,4)*RS(2,3)/400-RS(2,2)**2*RS
     #(2,4)*RS(2,1)/600-RS(2,2)*RS(2,4)**2*RS(2,3)/300-RS(2,1)*RS(2,2)*R
     #S(2,3)*RS(2,4)/300-RS(2,1)**4/600-RS(2,1)*RS(2,3)**3/600+RS(2,3)**
     #2*RS(2,4)**2/60-RS(2,1)**2*RS(2,3)**2/1200
      s11 = s12-RS(2,2)*RS(2,4)**3/600-RS(2,2)**3*RS(2,4)/1200-RS(2,2)**
     #2*RS(2,4)**2/1200-RS(2,2)*RS(2,3)**3/150-RS(2,2)**3*RS(2,3)/300+RS
     #(2,4)**3*RS(2,3)/40+RS(2,3)**3*RS(2,4)/120-RS(2,2)**2*RS(2,3)**2/2
     #00-RS(2,1)**3*RS(2,3)/1200-RS(2,1)**2*RS(2,4)**2/200-RS(2,1)**3*RS
     #(2,4)/300-RS(2,1)**3*RS(2,2)/600-RS(2,2)**4/600+RS(2,4)**4/30-RS(2
     #,1)*RS(2,4)**3/150-RS(2,1)**2*RS(2,2)**2/600-RS(2,1)*RS(2,2)**3/60
     #0
      s12 = RS(1,3)
      s10 = s11*s12
      s11 = (-RS(2,3)**3*RS(2,4)/60+RS(2,1)**2*RS(2,4)**2/40-RS(2,3)**2*
     #RS(2,4)**2/40+RS(2,1)*RS(2,4)**3/30+RS(2,1)**3*RS(2,4)/60-RS(2,4)*
     #*3*RS(2,3)/30-RS(2,3)**4/120+RS(2,1)**4/120)*RS(1,4)
      s9 = s10+s11
      ANTSSS(4) = s8+s9

C	--------------------------
C	COMPUTE FOR FACTOR ANTRRRR
C	--------------------------
      s5 = (-RS(1,1)**2*RS(1,2)**3/70+2.E0/105.E0*RS(1,4)**2*RS(1,1)**3-
     #RS(1,2)**5/210+RS(1,4)**5/210-2.E0/105.E0*RS(1,1)**3*RS(1,2)**2+RS
     #(1,4)**3*RS(1,1)**2/70+RS(1,1)**4*RS(1,4)/42+RS(1,4)**4*RS(1,1)/10
     #5-RS(1,1)*RS(1,2)**4/105-RS(1,1)**4*RS(1,2)/42)*RS(2,1)
      s9 = -RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/700-RS(1,4)**2*RS(1,1)**3
     #/315+RS(1,1)**2*RS(1,2)**3/105+2.E0/105.E0*RS(1,1)**4*RS(1,2)+RS(1
     #,1)*RS(1,2)**4/210-RS(1,4)**4*RS(1,1)/630-RS(1,1)**4*RS(1,4)/252-R
     #S(1,4)**3*RS(1,1)**2/420+RS(1,1)**3*RS(1,2)**2/70-RS(1,3)**4*RS(1,
     #2)/630-RS(1,3)**4*RS(1,1)/3150-RS(1,1)**2*RS(1,3)**3/4200-RS(1,3)*
     #*2*RS(1,2)**3/315
      s8 = s9-RS(1,2)**4*RS(1,3)/252-RS(1,3)**4*RS(1,4)/1260-RS(1,4)**3*
     #RS(1,3)**2/1260-RS(1,4)**4*RS(1,3)/1260-RS(1,4)**4*RS(1,2)/3150-RS
     #(1,1)**3*RS(1,3)**2/3150-RS(1,1)**4*RS(1,3)/1260-RS(1,3)**3*RS(1,2
     #)**2/420-RS(1,4)**2*RS(1,2)**3/3150-RS(1,4)**3*RS(1,2)**2/4200-RS(
     #1,4)*RS(1,2)**4/1260-RS(1,3)**3*RS(1,4)**2/1260-2.E0/1575.E0*RS(1,
     #4)*RS(1,2)**3*RS(1,3)-RS(1,3)**3*RS(1,4)*RS(1,1)/1575
      s9 = -RS(1,1)**2*RS(1,3)*RS(1,2)**2/420-RS(1,3)**2*RS(1,4)*RS(1,2)
     #**2/700-RS(1,1)*RS(1,4)**3*RS(1,2)/1050-RS(1,1)**3*RS(1,3)*RS(1,2)
     #/630-RS(1,1)*RS(1,3)**2*RS(1,2)**2/525-RS(1,4)**2*RS(1,2)**2*RS(1,
     #3)/1400-RS(1,4)**2*RS(1,3)**2*RS(1,1)/1050-2.E0/1575.E0*RS(1,3)**3
     #*RS(1,4)*RS(1,2)-RS(1,1)*RS(1,4)**2*RS(1,2)**2/1050-RS(1,1)**2*RS(
     #1,4)**2*RS(1,3)/700-RS(1,1)**2*RS(1,4)*RS(1,2)**2/420-RS(1,4)**3*R
     #S(1,2)*RS(1,3)/1575-RS(1,1)*RS(1,4)*RS(1,2)**3/630-RS(1,4)*RS(1,3)
     #**2*RS(1,1)**2/1400
      s7 = s9-RS(1,1)**3*RS(1,4)*RS(1,2)/315-RS(1,1)*RS(1,3)**3*RS(1,2)/
     #1050-2.E0/1575.E0*RS(1,4)*RS(1,1)**3*RS(1,3)-RS(1,1)**2*RS(1,4)**2
     #*RS(1,2)/525-RS(1,4)**2*RS(1,3)**2*RS(1,2)/1050-RS(1,3)**2*RS(1,1)
     #**2*RS(1,2)/1050-2.E0/1575.E0*RS(1,4)**3*RS(1,3)*RS(1,1)-RS(1,1)*R
     #S(1,3)*RS(1,2)**3/315-RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/525-RS(1,
     #1)*RS(1,3)*RS(1,4)*RS(1,2)**2/525-RS(1,1)*RS(1,3)**2*RS(1,4)*RS(1,
     #2)/700-RS(1,4)**5/1260-RS(1,3)**5/1260+RS(1,1)**5/42+s8
      s8 = RS(2,2)
      s6 = s7*s8
      s4 = s5+s6
      s5 = s4
      s9 = -RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/2100-RS(1,4)**2*RS(1,1)**
     #3/630+RS(1,1)**2*RS(1,2)**3/420+RS(1,1)**4*RS(1,2)/1260+RS(1,1)*RS
     #(1,2)**4/315-RS(1,4)**4*RS(1,1)/315-RS(1,1)**4*RS(1,4)/1260-RS(1,4
     #)**3*RS(1,1)**2/420+RS(1,1)**3*RS(1,2)**2/630+RS(1,3)**4*RS(1,2)/1
     #260+RS(1,3)**2*RS(1,2)**3/420+RS(1,2)**4*RS(1,3)/315-RS(1,3)**4*RS
     #(1,4)/1260-RS(1,4)**3*RS(1,3)**2/420-RS(1,4)**4*RS(1,3)/315-RS(1,4
     #)**4*RS(1,2)/2100+RS(1,3)**3*RS(1,2)**2/630+RS(1,4)**2*RS(1,2)**3/
     #12600-RS(1,4)**3*RS(1,2)**2/12600+RS(1,4)*RS(1,2)**4/2100-RS(1,3)*
     #*3*RS(1,4)**2/630+RS(1,4)*RS(1,2)**3*RS(1,3)/1575
      s10 = s9-RS(1,3)**3*RS(1,4)*RS(1,1)/3150+RS(1,1)**2*RS(1,3)*RS(1,2
     #)**2/1050+RS(1,3)**2*RS(1,4)*RS(1,2)**2/2100-RS(1,1)*RS(1,4)**3*RS
     #(1,2)/1575+RS(1,1)**3*RS(1,3)*RS(1,2)/3150+RS(1,1)*RS(1,3)**2*RS(1
     #,2)**2/1050-RS(1,4)**2*RS(1,3)**2*RS(1,1)/1050-RS(1,1)**2*RS(1,4)*
     #*2*RS(1,3)/1050+RS(1,1)**2*RS(1,4)*RS(1,2)**2/2100-RS(1,4)**3*RS(1
     #,2)*RS(1,3)/1575
      s8 = s10+RS(1,1)*RS(1,4)*RS(1,2)**3/1575-RS(1,4)*RS(1,3)**2*RS(1,1
     #)**2/4200+RS(1,1)*RS(1,3)**3*RS(1,2)/3150-RS(1,4)*RS(1,1)**3*RS(1,
     #3)/3150-RS(1,1)**2*RS(1,4)**2*RS(1,2)/2100-RS(1,4)**2*RS(1,3)**2*R
     #S(1,2)/2100+RS(1,3)**2*RS(1,1)**2*RS(1,2)/4200-RS(1,4)**3*RS(1,3)*
     #RS(1,1)/525+RS(1,1)*RS(1,3)*RS(1,2)**3/525+RS(1,1)*RS(1,3)*RS(1,4)
     #*RS(1,2)**2/2100-RS(1,4)**5/252+RS(1,2)**5/252
      s9 = RS(2,3)
      s7 = s8*s9
      s11 = RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/525-RS(1,4)**2*RS(1,1)**3
     #/70+RS(1,1)**2*RS(1,2)**3/420+RS(1,1)**4*RS(1,2)/252+RS(1,1)*RS(1,
     #2)**4/630-RS(1,4)**4*RS(1,1)/210-2.E0/105.E0*RS(1,1)**4*RS(1,4)-RS
     #(1,4)**3*RS(1,1)**2/105+RS(1,1)**3*RS(1,2)**2/315+RS(1,3)**4*RS(1,
     #2)/1260+RS(1,3)**4*RS(1,1)/3150+RS(1,1)**2*RS(1,3)**3/4200+RS(1,3)
     #**2*RS(1,2)**3/1260
      s10 = s11+RS(1,2)**4*RS(1,3)/1260+RS(1,3)**4*RS(1,4)/630+RS(1,4)**
     #3*RS(1,3)**2/315+RS(1,4)**4*RS(1,3)/252+RS(1,4)**4*RS(1,2)/1260+RS
     #(1,1)**3*RS(1,3)**2/3150+RS(1,1)**4*RS(1,3)/1260+RS(1,3)**3*RS(1,2
     #)**2/1260+RS(1,4)**2*RS(1,2)**3/4200+RS(1,4)**3*RS(1,2)**2/3150+RS
     #(1,4)*RS(1,2)**4/3150+RS(1,3)**3*RS(1,4)**2/420+RS(1,4)*RS(1,2)**3
     #*RS(1,3)/1575+RS(1,3)**3*RS(1,4)*RS(1,1)/1050
      s11 = RS(1,1)**2*RS(1,3)*RS(1,2)**2/700+RS(1,3)**2*RS(1,4)*RS(1,2)
     #**2/1050+RS(1,1)*RS(1,4)**3*RS(1,2)/630+2.E0/1575.E0*RS(1,1)**3*RS
     #(1,3)*RS(1,2)+RS(1,1)*RS(1,3)**2*RS(1,2)**2/1050+RS(1,4)**2*RS(1,2
     #)**2*RS(1,3)/1400+RS(1,4)**2*RS(1,3)**2*RS(1,1)/525+2.E0/1575.E0*R
     #S(1,3)**3*RS(1,4)*RS(1,2)+RS(1,1)*RS(1,4)**2*RS(1,2)**2/1050+RS(1,
     #1)**2*RS(1,4)**2*RS(1,3)/420+RS(1,1)**2*RS(1,4)*RS(1,2)**2/525+2.E
     #0/1575.E0*RS(1,4)**3*RS(1,2)*RS(1,3)+RS(1,1)*RS(1,4)*RS(1,2)**3/10
     #50+RS(1,4)*RS(1,3)**2*RS(1,1)**2/1050
      s9 = s11+RS(1,1)**3*RS(1,4)*RS(1,2)/315+RS(1,1)*RS(1,3)**3*RS(1,2)
     #/1575+RS(1,4)*RS(1,1)**3*RS(1,3)/630+RS(1,1)**2*RS(1,4)**2*RS(1,2)
     #/420+RS(1,4)**2*RS(1,3)**2*RS(1,2)/700+RS(1,3)**2*RS(1,1)**2*RS(1,
     #2)/1400+RS(1,4)**3*RS(1,3)*RS(1,1)/315+2.E0/1575.E0*RS(1,1)*RS(1,3
     #)*RS(1,2)**3+RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/525+RS(1,1)*RS(1,3
     #)*RS(1,4)*RS(1,2)**2/700+RS(1,1)*RS(1,3)**2*RS(1,4)*RS(1,2)/700+RS
     #(1,2)**5/1260+RS(1,3)**5/1260-RS(1,1)**5/42+s10
      s10 = RS(2,4)
      s8 = s9*s10
      s6 = s7+s8
      ANTRRRR(1) = s5+s6

      s9 = RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/700+RS(1,4)**2*RS(1,1)**3/
     #315-RS(1,1)**2*RS(1,2)**3/70-RS(1,1)**4*RS(1,2)/210-2.E0/105.E0*RS
     #(1,1)*RS(1,2)**4+RS(1,4)**4*RS(1,1)/630+RS(1,1)**4*RS(1,4)/252+RS(
     #1,4)**3*RS(1,1)**2/420-RS(1,1)**3*RS(1,2)**2/105+RS(1,3)**4*RS(1,2
     #)/630+RS(1,3)**4*RS(1,1)/3150+RS(1,1)**2*RS(1,3)**3/4200+RS(1,3)**
     #2*RS(1,2)**3/315
      s8 = s9+RS(1,2)**4*RS(1,3)/252+RS(1,3)**4*RS(1,4)/1260+RS(1,4)**3*
     #RS(1,3)**2/1260+RS(1,4)**4*RS(1,3)/1260+RS(1,4)**4*RS(1,2)/3150+RS
     #(1,1)**3*RS(1,3)**2/3150+RS(1,1)**4*RS(1,3)/1260+RS(1,3)**3*RS(1,2
     #)**2/420+RS(1,4)**2*RS(1,2)**3/3150+RS(1,4)**3*RS(1,2)**2/4200+RS(
     #1,4)*RS(1,2)**4/1260+RS(1,3)**3*RS(1,4)**2/1260+2.E0/1575.E0*RS(1,
     #4)*RS(1,2)**3*RS(1,3)+RS(1,3)**3*RS(1,4)*RS(1,1)/1575
      s9 = RS(1,1)**2*RS(1,3)*RS(1,2)**2/420+RS(1,3)**2*RS(1,4)*RS(1,2)*
     #*2/700+RS(1,1)*RS(1,4)**3*RS(1,2)/1050+RS(1,1)**3*RS(1,3)*RS(1,2)/
     #630+RS(1,1)*RS(1,3)**2*RS(1,2)**2/525+RS(1,4)**2*RS(1,2)**2*RS(1,3
     #)/1400+RS(1,4)**2*RS(1,3)**2*RS(1,1)/1050+2.E0/1575.E0*RS(1,3)**3*
     #RS(1,4)*RS(1,2)+RS(1,1)*RS(1,4)**2*RS(1,2)**2/1050+RS(1,1)**2*RS(1
     #,4)**2*RS(1,3)/700+RS(1,1)**2*RS(1,4)*RS(1,2)**2/420+RS(1,4)**3*RS
     #(1,2)*RS(1,3)/1575+RS(1,1)*RS(1,4)*RS(1,2)**3/630+RS(1,4)*RS(1,3)*
     #*2*RS(1,1)**2/1400
      s7 = s9+RS(1,1)**3*RS(1,4)*RS(1,2)/315+RS(1,1)*RS(1,3)**3*RS(1,2)/
     #1050+2.E0/1575.E0*RS(1,4)*RS(1,1)**3*RS(1,3)+RS(1,1)**2*RS(1,4)**2
     #*RS(1,2)/525+RS(1,4)**2*RS(1,3)**2*RS(1,2)/1050+RS(1,3)**2*RS(1,1)
     #**2*RS(1,2)/1050+2.E0/1575.E0*RS(1,4)**3*RS(1,3)*RS(1,1)+RS(1,1)*R
     #S(1,3)*RS(1,2)**3/315+RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/525+RS(1,
     #1)*RS(1,3)*RS(1,4)*RS(1,2)**2/525+RS(1,1)*RS(1,3)**2*RS(1,4)*RS(1,
     #2)/700+RS(1,4)**5/1260-RS(1,2)**5/42+RS(1,3)**5/1260+s8
      s8 = RS(2,1)
      s6 = s7*s8
      s7 = (-RS(1,3)**3*RS(1,2)**2/70-RS(1,3)**5/210-RS(1,3)**4*RS(1,2)/
     #105+RS(1,1)*RS(1,2)**4/42-2.E0/105.E0*RS(1,3)**2*RS(1,2)**3+RS(1,1
     #)**3*RS(1,2)**2/70+RS(1,1)**4*RS(1,2)/105-RS(1,2)**4*RS(1,3)/42+2.
     #E0/105.E0*RS(1,1)**2*RS(1,2)**3+RS(1,1)**5/210)*RS(2,2)
      s5 = s6+s7
      s6 = s5
      s11 = -RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/700-RS(1,4)**2*RS(1,1)**
     #3/1260-RS(1,1)**2*RS(1,2)**3/315-RS(1,1)**4*RS(1,2)/630-RS(1,1)*RS
     #(1,2)**4/252-RS(1,4)**4*RS(1,1)/1260-RS(1,1)**4*RS(1,4)/1260-RS(1,
     #4)**3*RS(1,1)**2/1260-RS(1,1)**3*RS(1,2)**2/420+RS(1,3)**4*RS(1,2)
     #/210-RS(1,3)**4*RS(1,1)/1260-RS(1,1)**2*RS(1,3)**3/3150+RS(1,3)**2
     #*RS(1,2)**3/70
      s10 = s11+2.E0/105.E0*RS(1,2)**4*RS(1,3)-RS(1,3)**4*RS(1,4)/252-RS
     #(1,4)**3*RS(1,3)**2/420-RS(1,4)**4*RS(1,3)/630-RS(1,4)**4*RS(1,2)/
     #3150-RS(1,1)**3*RS(1,3)**2/4200-RS(1,1)**4*RS(1,3)/3150+RS(1,3)**3
     #*RS(1,2)**2/105-RS(1,4)**2*RS(1,2)**3/3150-RS(1,4)**3*RS(1,2)**2/4
     #200-RS(1,4)*RS(1,2)**4/1260-RS(1,3)**3*RS(1,4)**2/315-RS(1,4)*RS(1
     #,2)**3*RS(1,3)/630-2.E0/1575.E0*RS(1,3)**3*RS(1,4)*RS(1,1)
      s11 = -RS(1,1)**2*RS(1,3)*RS(1,2)**2/525-RS(1,3)**2*RS(1,4)*RS(1,2
     #)**2/420-RS(1,1)*RS(1,4)**3*RS(1,2)/1575-RS(1,1)**3*RS(1,3)*RS(1,2
     #)/1050-RS(1,1)*RS(1,3)**2*RS(1,2)**2/420-RS(1,4)**2*RS(1,2)**2*RS(
     #1,3)/1050-RS(1,4)**2*RS(1,3)**2*RS(1,1)/700-RS(1,3)**3*RS(1,4)*RS(
     #1,2)/315-RS(1,1)*RS(1,4)**2*RS(1,2)**2/1400-RS(1,1)**2*RS(1,4)**2*
     #RS(1,3)/1050-RS(1,1)**2*RS(1,4)*RS(1,2)**2/700-RS(1,4)**3*RS(1,2)*
     #RS(1,3)/1050-2.E0/1575.E0*RS(1,1)*RS(1,4)*RS(1,2)**3-RS(1,4)*RS(1,
     #3)**2*RS(1,1)**2/1400
      s9 = s11-2.E0/1575.E0*RS(1,1)**3*RS(1,4)*RS(1,2)-RS(1,1)*RS(1,3)**
     #3*RS(1,2)/630-RS(1,4)*RS(1,1)**3*RS(1,3)/1575-RS(1,1)**2*RS(1,4)**
     #2*RS(1,2)/1050-RS(1,4)**2*RS(1,3)**2*RS(1,2)/525-RS(1,3)**2*RS(1,1
     #)**2*RS(1,2)/1050-2.E0/1575.E0*RS(1,4)**3*RS(1,3)*RS(1,1)-RS(1,1)*
     #RS(1,3)*RS(1,2)**3/315-RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/700-RS(1
     #,1)*RS(1,3)*RS(1,4)*RS(1,2)**2/525-RS(1,1)*RS(1,3)**2*RS(1,4)*RS(1
     #,2)/525-RS(1,4)**5/1260+RS(1,2)**5/42-RS(1,1)**5/1260+s10
      s10 = RS(2,3)
      s8 = s9*s10
      s11 = -RS(1,4)**2*RS(1,1)**3/420-RS(1,1)**2*RS(1,2)**3/630-RS(1,1)
     #**4*RS(1,2)/315-RS(1,1)*RS(1,2)**4/1260-RS(1,4)**4*RS(1,1)/1260-RS
     #(1,1)**4*RS(1,4)/315-RS(1,4)**3*RS(1,1)**2/630-RS(1,1)**3*RS(1,2)*
     #*2/420+RS(1,3)**4*RS(1,2)/315+RS(1,3)**4*RS(1,1)/2100+RS(1,1)**2*R
     #S(1,3)**3/12600+RS(1,3)**2*RS(1,2)**3/630+RS(1,2)**4*RS(1,3)/1260+
     #RS(1,3)**4*RS(1,4)/315+RS(1,4)**3*RS(1,3)**2/630+RS(1,4)**4*RS(1,3
     #)/1260-RS(1,1)**3*RS(1,3)**2/12600-RS(1,1)**4*RS(1,3)/2100+RS(1,3)
     #**3*RS(1,2)**2/420+RS(1,3)**3*RS(1,4)**2/420+RS(1,4)*RS(1,2)**3*RS
     #(1,3)/3150+RS(1,3)**3*RS(1,4)*RS(1,1)/1575
      s12 = s11-RS(1,1)**2*RS(1,3)*RS(1,2)**2/2100+RS(1,3)**2*RS(1,4)*RS
     #(1,2)**2/1050-RS(1,1)*RS(1,4)**3*RS(1,2)/3150-RS(1,1)**3*RS(1,3)*R
     #S(1,2)/1575+RS(1,1)*RS(1,3)**2*RS(1,2)**2/2100+RS(1,4)**2*RS(1,2)*
     #*2*RS(1,3)/4200+RS(1,4)**2*RS(1,3)**2*RS(1,1)/2100+RS(1,3)**3*RS(1
     #,4)*RS(1,2)/525-RS(1,1)*RS(1,4)**2*RS(1,2)**2/4200-RS(1,1)**2*RS(1
     #,4)**2*RS(1,3)/2100
      s10 = s12-RS(1,1)**2*RS(1,4)*RS(1,2)**2/1050+RS(1,4)**3*RS(1,2)*RS
     #(1,3)/3150-RS(1,1)*RS(1,4)*RS(1,2)**3/3150-RS(1,1)**3*RS(1,4)*RS(1
     #,2)/525+RS(1,1)*RS(1,3)**3*RS(1,2)/1575-RS(1,4)*RS(1,1)**3*RS(1,3)
     #/1575-RS(1,1)**2*RS(1,4)**2*RS(1,2)/1050+RS(1,4)**2*RS(1,3)**2*RS(
     #1,2)/1050-RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/2100+RS(1,1)*RS(1,3)*
     #*2*RS(1,4)*RS(1,2)/2100+RS(1,3)**5/252-RS(1,1)**5/252
      s11 = RS(2,4)
      s9 = s10*s11
      s7 = s8+s9
      ANTRRRR(2) = s6+s7

      s9 = RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/2100+RS(1,4)**2*RS(1,1)**3
     #/630-RS(1,1)**2*RS(1,2)**3/420-RS(1,1)**4*RS(1,2)/1260-RS(1,1)*RS(
     #1,2)**4/315+RS(1,4)**4*RS(1,1)/315+RS(1,1)**4*RS(1,4)/1260+RS(1,4)
     #**3*RS(1,1)**2/420-RS(1,1)**3*RS(1,2)**2/630-RS(1,3)**4*RS(1,2)/12
     #60-RS(1,3)**2*RS(1,2)**3/420-RS(1,2)**4*RS(1,3)/315+RS(1,3)**4*RS(
     #1,4)/1260+RS(1,4)**3*RS(1,3)**2/420+RS(1,4)**4*RS(1,3)/315+RS(1,4)
     #**4*RS(1,2)/2100-RS(1,3)**3*RS(1,2)**2/630-RS(1,4)**2*RS(1,2)**3/1
     #2600+RS(1,4)**3*RS(1,2)**2/12600-RS(1,4)*RS(1,2)**4/2100+RS(1,3)**
     #3*RS(1,4)**2/630-RS(1,4)*RS(1,2)**3*RS(1,3)/1575
      s10 = s9+RS(1,3)**3*RS(1,4)*RS(1,1)/3150-RS(1,1)**2*RS(1,3)*RS(1,2
     #)**2/1050-RS(1,3)**2*RS(1,4)*RS(1,2)**2/2100+RS(1,1)*RS(1,4)**3*RS
     #(1,2)/1575-RS(1,1)**3*RS(1,3)*RS(1,2)/3150-RS(1,1)*RS(1,3)**2*RS(1
     #,2)**2/1050+RS(1,4)**2*RS(1,3)**2*RS(1,1)/1050+RS(1,1)**2*RS(1,4)*
     #*2*RS(1,3)/1050-RS(1,1)**2*RS(1,4)*RS(1,2)**2/2100+RS(1,4)**3*RS(1
     #,2)*RS(1,3)/1575
      s8 = s10-RS(1,1)*RS(1,4)*RS(1,2)**3/1575+RS(1,4)*RS(1,3)**2*RS(1,1
     #)**2/4200-RS(1,1)*RS(1,3)**3*RS(1,2)/3150+RS(1,4)*RS(1,1)**3*RS(1,
     #3)/3150+RS(1,1)**2*RS(1,4)**2*RS(1,2)/2100+RS(1,4)**2*RS(1,3)**2*R
     #S(1,2)/2100-RS(1,3)**2*RS(1,1)**2*RS(1,2)/4200+RS(1,4)**3*RS(1,3)*
     #RS(1,1)/525-RS(1,1)*RS(1,3)*RS(1,2)**3/525-RS(1,1)*RS(1,3)*RS(1,4)
     #*RS(1,2)**2/2100+RS(1,4)**5/252-RS(1,2)**5/252
      s9 = RS(2,1)
      s7 = s8*s9
      s11 = RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/700+RS(1,4)**2*RS(1,1)**3
     #/1260+RS(1,1)**2*RS(1,2)**3/315+RS(1,1)**4*RS(1,2)/630+RS(1,1)*RS(
     #1,2)**4/252+RS(1,4)**4*RS(1,1)/1260+RS(1,1)**4*RS(1,4)/1260+RS(1,4
     #)**3*RS(1,1)**2/1260+RS(1,1)**3*RS(1,2)**2/420-2.E0/105.E0*RS(1,3)
     #**4*RS(1,2)+RS(1,3)**4*RS(1,1)/1260+RS(1,1)**2*RS(1,3)**3/3150-RS(
     #1,3)**2*RS(1,2)**3/105
      s10 = s11-RS(1,2)**4*RS(1,3)/210+RS(1,3)**4*RS(1,4)/252+RS(1,4)**3
     #*RS(1,3)**2/420+RS(1,4)**4*RS(1,3)/630+RS(1,4)**4*RS(1,2)/3150+RS(
     #1,1)**3*RS(1,3)**2/4200+RS(1,1)**4*RS(1,3)/3150-RS(1,3)**3*RS(1,2)
     #**2/70+RS(1,4)**2*RS(1,2)**3/3150+RS(1,4)**3*RS(1,2)**2/4200+RS(1,
     #4)*RS(1,2)**4/1260+RS(1,3)**3*RS(1,4)**2/315+RS(1,4)*RS(1,2)**3*RS
     #(1,3)/630+2.E0/1575.E0*RS(1,3)**3*RS(1,4)*RS(1,1)
      s11 = RS(1,1)**2*RS(1,3)*RS(1,2)**2/525+RS(1,3)**2*RS(1,4)*RS(1,2)
     #**2/420+RS(1,1)*RS(1,4)**3*RS(1,2)/1575+RS(1,1)**3*RS(1,3)*RS(1,2)
     #/1050+RS(1,1)*RS(1,3)**2*RS(1,2)**2/420+RS(1,4)**2*RS(1,2)**2*RS(1
     #,3)/1050+RS(1,4)**2*RS(1,3)**2*RS(1,1)/700+RS(1,3)**3*RS(1,4)*RS(1
     #,2)/315+RS(1,1)*RS(1,4)**2*RS(1,2)**2/1400+RS(1,1)**2*RS(1,4)**2*R
     #S(1,3)/1050+RS(1,1)**2*RS(1,4)*RS(1,2)**2/700+RS(1,4)**3*RS(1,2)*R
     #S(1,3)/1050+2.E0/1575.E0*RS(1,1)*RS(1,4)*RS(1,2)**3+RS(1,4)*RS(1,3
     #)**2*RS(1,1)**2/1400
      s9 = s11+2.E0/1575.E0*RS(1,1)**3*RS(1,4)*RS(1,2)+RS(1,1)*RS(1,3)**
     #3*RS(1,2)/630+RS(1,4)*RS(1,1)**3*RS(1,3)/1575+RS(1,1)**2*RS(1,4)**
     #2*RS(1,2)/1050+RS(1,4)**2*RS(1,3)**2*RS(1,2)/525+RS(1,3)**2*RS(1,1
     #)**2*RS(1,2)/1050+2.E0/1575.E0*RS(1,4)**3*RS(1,3)*RS(1,1)+RS(1,1)*
     #RS(1,3)*RS(1,2)**3/315+RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/700+RS(1
     #,1)*RS(1,3)*RS(1,4)*RS(1,2)**2/525+RS(1,1)*RS(1,3)**2*RS(1,4)*RS(1
     #,2)/525+RS(1,4)**5/1260-RS(1,3)**5/42+RS(1,1)**5/1260+s10
      s10 = RS(2,2)
      s8 = s9*s10
      s6 = s7+s8
      s7 = s6
      s9 = (RS(1,3)**2*RS(1,2)**3/70-RS(1,3)**4*RS(1,4)/42+RS(1,2)**4*RS
     #(1,3)/105-2.E0/105.E0*RS(1,3)**3*RS(1,4)**2-RS(1,4)**3*RS(1,3)**2/
     #70+2.E0/105.E0*RS(1,3)**3*RS(1,2)**2+RS(1,3)**4*RS(1,2)/42+RS(1,2)
     #**5/210-RS(1,4)**4*RS(1,3)/105-RS(1,4)**5/210)*RS(2,3)
      s13 = -RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/525-RS(1,4)**2*RS(1,1)**
     #3/420-RS(1,1)**2*RS(1,2)**3/1260-RS(1,1)**4*RS(1,2)/1260-RS(1,1)*R
     #S(1,2)**4/1260-RS(1,4)**4*RS(1,1)/252-RS(1,1)**4*RS(1,4)/630-RS(1,
     #4)**3*RS(1,1)**2/315-RS(1,1)**3*RS(1,2)**2/1260-RS(1,3)**4*RS(1,2)
     #/252-RS(1,3)**4*RS(1,1)/1260-RS(1,1)**2*RS(1,3)**3/3150-RS(1,3)**2
     #*RS(1,2)**3/420
      s12 = s13-RS(1,2)**4*RS(1,3)/630+2.E0/105.E0*RS(1,3)**4*RS(1,4)+RS
     #(1,4)**3*RS(1,3)**2/105+RS(1,4)**4*RS(1,3)/210-RS(1,4)**4*RS(1,2)/
     #1260-RS(1,1)**3*RS(1,3)**2/4200-RS(1,1)**4*RS(1,3)/3150-RS(1,3)**3
     #*RS(1,2)**2/315-RS(1,4)**2*RS(1,2)**3/4200-RS(1,4)**3*RS(1,2)**2/3
     #150-RS(1,4)*RS(1,2)**4/3150+RS(1,3)**3*RS(1,4)**2/70-RS(1,4)*RS(1,
     #2)**3*RS(1,3)/1050-RS(1,3)**3*RS(1,4)*RS(1,1)/630
      s13 = -RS(1,1)**2*RS(1,3)*RS(1,2)**2/1050-RS(1,3)**2*RS(1,4)*RS(1,
     #2)**2/525-2.E0/1575.E0*RS(1,1)*RS(1,4)**3*RS(1,2)-RS(1,1)**3*RS(1,
     #3)*RS(1,2)/1575-RS(1,1)*RS(1,3)**2*RS(1,2)**2/700-RS(1,4)**2*RS(1,
     #2)**2*RS(1,3)/1050-RS(1,4)**2*RS(1,3)**2*RS(1,1)/420-RS(1,3)**3*RS
     #(1,4)*RS(1,2)/315-RS(1,1)*RS(1,4)**2*RS(1,2)**2/1400-RS(1,1)**2*RS
     #(1,4)**2*RS(1,3)/525-RS(1,1)**2*RS(1,4)*RS(1,2)**2/1050-RS(1,4)**3
     #*RS(1,2)*RS(1,3)/630-RS(1,1)*RS(1,4)*RS(1,2)**3/1575-RS(1,4)*RS(1,
     #3)**2*RS(1,1)**2/1050
      s11 = s13-2.E0/1575.E0*RS(1,1)**3*RS(1,4)*RS(1,2)-2.E0/1575.E0*RS(
     #1,1)*RS(1,3)**3*RS(1,2)-RS(1,4)*RS(1,1)**3*RS(1,3)/1050-RS(1,1)**2
     #*RS(1,4)**2*RS(1,2)/700-RS(1,4)**2*RS(1,3)**2*RS(1,2)/420-RS(1,3)*
     #*2*RS(1,1)**2*RS(1,2)/1400-RS(1,4)**3*RS(1,3)*RS(1,1)/315-2.E0/157
     #5.E0*RS(1,1)*RS(1,3)*RS(1,2)**3+s12-RS(1,1)**2*RS(1,3)*RS(1,4)*RS(
     #1,2)/700-RS(1,1)*RS(1,3)*RS(1,4)*RS(1,2)**2/700-RS(1,1)*RS(1,3)**2
     #*RS(1,4)*RS(1,2)/525-RS(1,2)**5/1260+RS(1,3)**5/42-RS(1,1)**5/1260
      s12 = RS(2,4)
      s10 = s11*s12
      s8 = s9+s10
      ANTRRRR(3) = s7+s8

      s11 = -RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/525+RS(1,4)**2*RS(1,1)**
     #3/105-RS(1,1)**2*RS(1,2)**3/420-RS(1,1)**4*RS(1,2)/252-RS(1,1)*RS(
     #1,2)**4/630+2.E0/105.E0*RS(1,4)**4*RS(1,1)+RS(1,1)**4*RS(1,4)/210+
     #RS(1,4)**3*RS(1,1)**2/70-RS(1,1)**3*RS(1,2)**2/315-RS(1,3)**4*RS(1
     #,2)/1260-RS(1,3)**4*RS(1,1)/3150-RS(1,1)**2*RS(1,3)**3/4200-RS(1,3
     #)**2*RS(1,2)**3/1260
      s10 = s11-RS(1,2)**4*RS(1,3)/1260-RS(1,3)**4*RS(1,4)/630-RS(1,4)**
     #3*RS(1,3)**2/315-RS(1,4)**4*RS(1,3)/252-RS(1,4)**4*RS(1,2)/1260-RS
     #(1,1)**3*RS(1,3)**2/3150-RS(1,1)**4*RS(1,3)/1260-RS(1,3)**3*RS(1,2
     #)**2/1260-RS(1,4)**2*RS(1,2)**3/4200-RS(1,4)**3*RS(1,2)**2/3150-RS
     #(1,4)*RS(1,2)**4/3150-RS(1,3)**3*RS(1,4)**2/420-RS(1,4)*RS(1,2)**3
     #*RS(1,3)/1575-RS(1,3)**3*RS(1,4)*RS(1,1)/1050
      s11 = -RS(1,1)**2*RS(1,3)*RS(1,2)**2/700-RS(1,3)**2*RS(1,4)*RS(1,2
     #)**2/1050-RS(1,1)*RS(1,4)**3*RS(1,2)/630-2.E0/1575.E0*RS(1,1)**3*R
     #S(1,3)*RS(1,2)-RS(1,1)*RS(1,3)**2*RS(1,2)**2/1050-RS(1,4)**2*RS(1,
     #2)**2*RS(1,3)/1400-RS(1,4)**2*RS(1,3)**2*RS(1,1)/525-2.E0/1575.E0*
     #RS(1,3)**3*RS(1,4)*RS(1,2)-RS(1,1)*RS(1,4)**2*RS(1,2)**2/1050-RS(1
     #,1)**2*RS(1,4)**2*RS(1,3)/420-RS(1,1)**2*RS(1,4)*RS(1,2)**2/525-2.
     #E0/1575.E0*RS(1,4)**3*RS(1,2)*RS(1,3)-RS(1,1)*RS(1,4)*RS(1,2)**3/1
     #050-RS(1,4)*RS(1,3)**2*RS(1,1)**2/1050
      s9 = s11-RS(1,1)**3*RS(1,4)*RS(1,2)/315-RS(1,1)*RS(1,3)**3*RS(1,2)
     #/1575-RS(1,4)*RS(1,1)**3*RS(1,3)/630-RS(1,1)**2*RS(1,4)**2*RS(1,2)
     #/420-RS(1,4)**2*RS(1,3)**2*RS(1,2)/700-RS(1,3)**2*RS(1,1)**2*RS(1,
     #2)/1400-RS(1,4)**3*RS(1,3)*RS(1,1)/315-2.E0/1575.E0*RS(1,1)*RS(1,3
     #)*RS(1,2)**3-RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/525-RS(1,1)*RS(1,3
     #)*RS(1,4)*RS(1,2)**2/700-RS(1,1)*RS(1,3)**2*RS(1,4)*RS(1,2)/700+RS
     #(1,4)**5/42-RS(1,2)**5/1260-RS(1,3)**5/1260+s10
      s10 = RS(2,1)
      s8 = s9*s10
      s11 = RS(1,4)**2*RS(1,1)**3/420+RS(1,1)**2*RS(1,2)**3/630+RS(1,1)*
     #*4*RS(1,2)/315+RS(1,1)*RS(1,2)**4/1260+RS(1,4)**4*RS(1,1)/1260+RS(
     #1,1)**4*RS(1,4)/315+RS(1,4)**3*RS(1,1)**2/630+RS(1,1)**3*RS(1,2)**
     #2/420-RS(1,3)**4*RS(1,2)/315-RS(1,3)**4*RS(1,1)/2100-RS(1,1)**2*RS
     #(1,3)**3/12600-RS(1,3)**2*RS(1,2)**3/630-RS(1,2)**4*RS(1,3)/1260-R
     #S(1,3)**4*RS(1,4)/315-RS(1,4)**3*RS(1,3)**2/630-RS(1,4)**4*RS(1,3)
     #/1260+RS(1,1)**3*RS(1,3)**2/12600+RS(1,1)**4*RS(1,3)/2100-RS(1,3)*
     #*3*RS(1,2)**2/420-RS(1,3)**3*RS(1,4)**2/420-RS(1,4)*RS(1,2)**3*RS(
     #1,3)/3150-RS(1,3)**3*RS(1,4)*RS(1,1)/1575
      s12 = s11+RS(1,1)**2*RS(1,3)*RS(1,2)**2/2100-RS(1,3)**2*RS(1,4)*RS
     #(1,2)**2/1050+RS(1,1)*RS(1,4)**3*RS(1,2)/3150+RS(1,1)**3*RS(1,3)*R
     #S(1,2)/1575-RS(1,1)*RS(1,3)**2*RS(1,2)**2/2100-RS(1,4)**2*RS(1,2)*
     #*2*RS(1,3)/4200-RS(1,4)**2*RS(1,3)**2*RS(1,1)/2100-RS(1,3)**3*RS(1
     #,4)*RS(1,2)/525+RS(1,1)*RS(1,4)**2*RS(1,2)**2/4200+RS(1,1)**2*RS(1
     #,4)**2*RS(1,3)/2100
      s10 = s12+RS(1,1)**2*RS(1,4)*RS(1,2)**2/1050-RS(1,4)**3*RS(1,2)*RS
     #(1,3)/3150+RS(1,1)*RS(1,4)*RS(1,2)**3/3150+RS(1,1)**3*RS(1,4)*RS(1
     #,2)/525-RS(1,1)*RS(1,3)**3*RS(1,2)/1575+RS(1,4)*RS(1,1)**3*RS(1,3)
     #/1575+RS(1,1)**2*RS(1,4)**2*RS(1,2)/1050-RS(1,4)**2*RS(1,3)**2*RS(
     #1,2)/1050+RS(1,1)**2*RS(1,3)*RS(1,4)*RS(1,2)/2100-RS(1,1)*RS(1,3)*
     #*2*RS(1,4)*RS(1,2)/2100-RS(1,3)**5/252+RS(1,1)**5/252
      s11 = RS(2,2)
      s9 = s10*s11
      s7 = s8+s9
      s8 = s7
      s13 = RS(1,4)**2*RS(1,3)*RS(1,1)*RS(1,2)/525+RS(1,4)**2*RS(1,1)**3
     #/420+RS(1,1)**2*RS(1,2)**3/1260+RS(1,1)**4*RS(1,2)/1260+RS(1,1)*RS
     #(1,2)**4/1260+RS(1,4)**4*RS(1,1)/252+RS(1,1)**4*RS(1,4)/630+RS(1,4
     #)**3*RS(1,1)**2/315+RS(1,1)**3*RS(1,2)**2/1260+RS(1,3)**4*RS(1,2)/
     #252+RS(1,3)**4*RS(1,1)/1260+RS(1,1)**2*RS(1,3)**3/3150+RS(1,3)**2*
     #RS(1,2)**3/420
      s12 = s13+RS(1,2)**4*RS(1,3)/630-RS(1,3)**4*RS(1,4)/210-RS(1,4)**3
     #*RS(1,3)**2/70-2.E0/105.E0*RS(1,4)**4*RS(1,3)+RS(1,4)**4*RS(1,2)/1
     #260+RS(1,1)**3*RS(1,3)**2/4200+RS(1,1)**4*RS(1,3)/3150+RS(1,3)**3*
     #RS(1,2)**2/315+RS(1,4)**2*RS(1,2)**3/4200+RS(1,4)**3*RS(1,2)**2/31
     #50+RS(1,4)*RS(1,2)**4/3150-RS(1,3)**3*RS(1,4)**2/105+RS(1,4)*RS(1,
     #2)**3*RS(1,3)/1050+RS(1,3)**3*RS(1,4)*RS(1,1)/630
      s13 = RS(1,1)**2*RS(1,3)*RS(1,2)**2/1050+RS(1,3)**2*RS(1,4)*RS(1,2
     #)**2/525+2.E0/1575.E0*RS(1,1)*RS(1,4)**3*RS(1,2)+RS(1,1)**3*RS(1,3
     #)*RS(1,2)/1575+RS(1,1)*RS(1,3)**2*RS(1,2)**2/700+RS(1,4)**2*RS(1,2
     #)**2*RS(1,3)/1050+RS(1,4)**2*RS(1,3)**2*RS(1,1)/420+RS(1,3)**3*RS(
     #1,4)*RS(1,2)/315+RS(1,1)*RS(1,4)**2*RS(1,2)**2/1400+RS(1,1)**2*RS(
     #1,4)**2*RS(1,3)/525+RS(1,1)**2*RS(1,4)*RS(1,2)**2/1050+RS(1,4)**3*
     #RS(1,2)*RS(1,3)/630+RS(1,1)*RS(1,4)*RS(1,2)**3/1575+RS(1,4)*RS(1,3
     #)**2*RS(1,1)**2/1050
      s11 = s13+2.E0/1575.E0*RS(1,1)**3*RS(1,4)*RS(1,2)+2.E0/1575.E0*RS(
     #1,1)*RS(1,3)**3*RS(1,2)+RS(1,4)*RS(1,1)**3*RS(1,3)/1050+RS(1,1)**2
     #*RS(1,4)**2*RS(1,2)/700+RS(1,4)**2*RS(1,3)**2*RS(1,2)/420+RS(1,3)*
     #*2*RS(1,1)**2*RS(1,2)/1400+RS(1,4)**3*RS(1,3)*RS(1,1)/315+2.E0/157
     #5.E0*RS(1,1)*RS(1,3)*RS(1,2)**3+s12+RS(1,1)**2*RS(1,3)*RS(1,4)*RS(
     #1,2)/700+RS(1,1)*RS(1,3)*RS(1,4)*RS(1,2)**2/700+RS(1,1)*RS(1,3)**2
     #*RS(1,4)*RS(1,2)/525-RS(1,4)**5/42+RS(1,2)**5/1260+RS(1,1)**5/1260
      s12 = RS(2,3)
      s10 = s11*s12
      s11 = (2.E0/105.E0*RS(1,4)**3*RS(1,3)**2-RS(1,4)**2*RS(1,1)**3/70-
     #RS(1,1)**4*RS(1,4)/105+RS(1,3)**5/210-2.E0/105.E0*RS(1,4)**3*RS(1,
     #1)**2-RS(1,4)**4*RS(1,1)/42-RS(1,1)**5/210+RS(1,3)**4*RS(1,4)/105+
     #RS(1,3)**3*RS(1,4)**2/70+RS(1,4)**4*RS(1,3)/42)*RS(2,4)
      s9 = s10+s11
      ANTRRRR(4) = s8+s9

C	--------------------------
C	COMPUTE FOR FACTOR ANTSSSS
C	--------------------------
      s5 = ((RS(2,2)/42-RS(2,4)/42)*RS(2,1)**4+(-2.E0/105.E0*RS(2,4)**2+
     #2.E0/105.E0*RS(2,2)**2)*RS(2,1)**3+(RS(2,2)**3/70-RS(2,4)**3/70)*R
     #S(2,1)**2+(RS(2,2)**4/105-RS(2,4)**4/105)*RS(2,1)-RS(2,4)**5/210+R
     #S(2,2)**5/210)*RS(1,1)
      s9 = -RS(2,1)**5/42+(-2.E0/105.E0*RS(2,2)+RS(2,4)/252+RS(2,3)/1260
     #)*RS(2,1)**4+(-RS(2,2)**2/70+(RS(2,3)/630+RS(2,4)/315)*RS(2,2)+RS(
     #2,4)**2/315+RS(2,3)**2/3150+2.E0/1575.E0*RS(2,3)*RS(2,4))*RS(2,1)*
     #*3
      s10 = s9+(-RS(2,2)**3/105+(RS(2,3)/420+RS(2,4)/420)*RS(2,2)**2+(RS
     #(2,4)**2/525+RS(2,3)*RS(2,4)/525+RS(2,3)**2/1050)*RS(2,2)+RS(2,4)*
     #*3/420+RS(2,4)**2*RS(2,3)/700+RS(2,3)**3/4200+RS(2,3)**2*RS(2,4)/1
     #400)*RS(2,1)**2
      s8 = s10+(-RS(2,2)**4/210+(RS(2,4)/630+RS(2,3)/315)*RS(2,2)**3+(RS
     #(2,3)**2/525+RS(2,3)*RS(2,4)/525+RS(2,4)**2/1050)*RS(2,2)**2+(RS(2
     #,3)**3/1050+RS(2,4)**2*RS(2,3)/700+RS(2,4)**3/1050+RS(2,3)**2*RS(2
     #,4)/700)*RS(2,2)+2.E0/1575.E0*RS(2,4)**3*RS(2,3)+RS(2,3)**3*RS(2,4
     #)/1575+RS(2,3)**2*RS(2,4)**2/1050+RS(2,3)**4/3150+RS(2,4)**4/630)*
     #RS(2,1)+(RS(2,3)/252+RS(2,4)/1260)*RS(2,2)**4+(RS(2,3)**2/315+2.E0
     #/1575.E0*RS(2,3)*RS(2,4)+RS(2,4)**2/3150)*RS(2,2)**3
      s7 = s8+(RS(2,3)**3/420+RS(2,4)**3/4200+RS(2,4)**2*RS(2,3)/1400+RS
     #(2,3)**2*RS(2,4)/700)*RS(2,2)**2+(RS(2,3)**4/630+RS(2,4)**3*RS(2,3
     #)/1575+RS(2,3)**2*RS(2,4)**2/1050+2.E0/1575.E0*RS(2,3)**3*RS(2,4)+
     #RS(2,4)**4/3150)*RS(2,2)+RS(2,4)**3*RS(2,3)**2/1260+RS(2,4)**4*RS(
     #2,3)/1260+RS(2,4)**5/1260+RS(2,4)**2*RS(2,3)**3/1260+RS(2,3)**5/12
     #60+RS(2,3)**4*RS(2,4)/1260
      s8 = RS(1,2)
      s6 = s7*s8
      s4 = s5+s6
      s5 = s4
      s10 = (-RS(2,2)/1260+RS(2,4)/1260)*RS(2,1)**4+(-RS(2,2)*RS(2,3)/31
     #50+RS(2,3)*RS(2,4)/3150+RS(2,4)**2/630-RS(2,2)**2/630)*RS(2,1)**3+
     #(-RS(2,2)**3/420+(-RS(2,4)/2100-RS(2,3)/1050)*RS(2,2)**2+(RS(2,4)*
     #*2/2100-RS(2,3)**2/4200)*RS(2,2)+RS(2,4)**2*RS(2,3)/1050+RS(2,3)**
     #2*RS(2,4)/4200+RS(2,4)**3/420)*RS(2,1)**2
      s9 = s10+(-RS(2,2)**4/315+(-RS(2,4)/1575-RS(2,3)/525)*RS(2,2)**3+(
     #-RS(2,3)**2/1050-RS(2,3)*RS(2,4)/2100)*RS(2,2)**2+(RS(2,4)**3/1575
     #-RS(2,3)**3/3150+RS(2,4)**2*RS(2,3)/2100)*RS(2,2)+RS(2,3)**3*RS(2,
     #4)/3150+RS(2,3)**2*RS(2,4)**2/1050+RS(2,4)**3*RS(2,3)/525+RS(2,4)*
     #*4/315)*RS(2,1)-RS(2,2)**5/252+(-RS(2,4)/2100-RS(2,3)/315)*RS(2,2)
     #**4+(-RS(2,3)**2/420-RS(2,3)*RS(2,4)/1575-RS(2,4)**2/12600)*RS(2,2
     #)**3
      s8 = s9+(-RS(2,3)**2*RS(2,4)/2100-RS(2,3)**3/630+RS(2,4)**3/12600)
     #*RS(2,2)**2+(-RS(2,3)**4/1260+RS(2,4)**3*RS(2,3)/1575+RS(2,4)**4/2
     #100+RS(2,3)**2*RS(2,4)**2/2100)*RS(2,2)+RS(2,4)**4*RS(2,3)/315+RS(
     #2,4)**3*RS(2,3)**2/420+RS(2,4)**2*RS(2,3)**3/630+RS(2,3)**4*RS(2,4
     #)/1260+RS(2,4)**5/252
      s9 = RS(1,3)
      s7 = s8*s9
      s11 = RS(2,1)**5/42+(-RS(2,2)/252+2.E0/105.E0*RS(2,4)-RS(2,3)/1260
     #)*RS(2,1)**4+(-RS(2,2)**2/315+(-2.E0/1575.E0*RS(2,3)-RS(2,4)/315)*
     #RS(2,2)+RS(2,4)**2/70-RS(2,3)**2/3150-RS(2,3)*RS(2,4)/630)*RS(2,1)
     #**3
      s10 = s11+(-RS(2,2)**3/420+(-RS(2,3)/700-RS(2,4)/525)*RS(2,2)**2+(
     #-RS(2,3)*RS(2,4)/525-RS(2,4)**2/420-RS(2,3)**2/1400)*RS(2,2)+RS(2,
     #4)**3/105-RS(2,4)**2*RS(2,3)/420-RS(2,3)**2*RS(2,4)/1050-RS(2,3)**
     #3/4200)*RS(2,1)**2+(-RS(2,2)**4/630+(-RS(2,4)/1050-2.E0/1575.E0*RS
     #(2,3))*RS(2,2)**3+(-RS(2,3)**2/1050-RS(2,4)**2/1050-RS(2,3)*RS(2,4
     #)/700)*RS(2,2)**2+(-RS(2,3)**3/1575-RS(2,3)**2*RS(2,4)/700-RS(2,4)
     #**2*RS(2,3)/525-RS(2,4)**3/630)*RS(2,2)-RS(2,3)**2*RS(2,4)**2/525-
     #RS(2,3)**3*RS(2,4)/1050-RS(2,3)**4/3150+RS(2,4)**4/210-RS(2,4)**3*
     #RS(2,3)/315)*RS(2,1)-RS(2,2)**5/1260+(-RS(2,3)/1260-RS(2,4)/3150)*
     #RS(2,2)**4
      s9 = s10+(-RS(2,3)**2/1260-RS(2,3)*RS(2,4)/1575-RS(2,4)**2/4200)*R
     #S(2,2)**3+(-RS(2,3)**3/1260-RS(2,4)**3/3150-RS(2,4)**2*RS(2,3)/140
     #0-RS(2,3)**2*RS(2,4)/1050)*RS(2,2)**2+(-RS(2,3)**4/1260-2.E0/1575.
     #E0*RS(2,3)**3*RS(2,4)-2.E0/1575.E0*RS(2,4)**3*RS(2,3)-RS(2,3)**2*R
     #S(2,4)**2/700-RS(2,4)**4/1260)*RS(2,2)-RS(2,4)**3*RS(2,3)**2/315-R
     #S(2,4)**4*RS(2,3)/252-RS(2,3)**4*RS(2,4)/630-RS(2,4)**2*RS(2,3)**3
     #/420-RS(2,3)**5/1260
      s10 = RS(1,4)
      s8 = s9*s10
      s6 = s7+s8
      ANTSSSS(1) = s5+s6

      s9 = (RS(2,2)/210-RS(2,4)/252-RS(2,3)/1260)*RS(2,1)**4+(RS(2,2)**2
     #/105+(-RS(2,4)/315-RS(2,3)/630)*RS(2,2)-RS(2,4)**2/315-RS(2,3)**2/
     #3150-2.E0/1575.E0*RS(2,3)*RS(2,4))*RS(2,1)**3+(RS(2,2)**3/70+(-RS(
     #2,4)/420-RS(2,3)/420)*RS(2,2)**2+(-RS(2,3)*RS(2,4)/525-RS(2,4)**2/
     #525-RS(2,3)**2/1050)*RS(2,2)-RS(2,4)**3/420-RS(2,4)**2*RS(2,3)/700
     #-RS(2,3)**3/4200-RS(2,3)**2*RS(2,4)/1400)*RS(2,1)**2
      s8 = s9+(2.E0/105.E0*RS(2,2)**4+(-RS(2,4)/630-RS(2,3)/315)*RS(2,2)
     #**3+(-RS(2,3)**2/525-RS(2,3)*RS(2,4)/525-RS(2,4)**2/1050)*RS(2,2)*
     #*2+(-RS(2,3)**3/1050-RS(2,4)**3/1050-RS(2,4)**2*RS(2,3)/700-RS(2,3
     #)**2*RS(2,4)/700)*RS(2,2)-RS(2,3)**3*RS(2,4)/1575-2.E0/1575.E0*RS(
     #2,4)**3*RS(2,3)-RS(2,3)**4/3150-RS(2,3)**2*RS(2,4)**2/1050-RS(2,4)
     #**4/630)*RS(2,1)+RS(2,2)**5/42+(-RS(2,4)/1260-RS(2,3)/252)*RS(2,2)
     #**4+(-RS(2,3)**2/315-RS(2,4)**2/3150-2.E0/1575.E0*RS(2,3)*RS(2,4))
     #*RS(2,2)**3
      s7 = s8+(-RS(2,3)**3/420-RS(2,3)**2*RS(2,4)/700-RS(2,4)**3/4200-RS
     #(2,4)**2*RS(2,3)/1400)*RS(2,2)**2+(-RS(2,3)**4/630-RS(2,3)**2*RS(2
     #,4)**2/1050-2.E0/1575.E0*RS(2,3)**3*RS(2,4)-RS(2,4)**3*RS(2,3)/157
     #5-RS(2,4)**4/3150)*RS(2,2)-RS(2,4)**3*RS(2,3)**2/1260-RS(2,4)**4*R
     #S(2,3)/1260-RS(2,4)**5/1260-RS(2,4)**2*RS(2,3)**3/1260-RS(2,3)**4*
     #RS(2,4)/1260-RS(2,3)**5/1260
      s8 = RS(1,1)
      s6 = s7*s8
      s7 = (-RS(2,2)*RS(2,1)**4/105-RS(2,2)**4*RS(2,1)/42+RS(2,3)**4*RS(
     #2,2)/105-RS(2,1)**5/210+RS(2,3)**3*RS(2,2)**2/70+2.E0/105.E0*RS(2,
     #3)**2*RS(2,2)**3+RS(2,3)**5/210-RS(2,2)**2*RS(2,1)**3/70-2.E0/105.
     #E0*RS(2,2)**3*RS(2,1)**2+RS(2,2)**4*RS(2,3)/42)*RS(1,2)
      s5 = s6+s7
      s6 = s5
      s11 = RS(2,1)**5/1260+(RS(2,2)/630+RS(2,4)/1260+RS(2,3)/3150)*RS(2
     #,1)**4+(RS(2,2)**2/420+(RS(2,3)/1050+2.E0/1575.E0*RS(2,4))*RS(2,2)
     #+RS(2,4)**2/1260+RS(2,3)**2/4200+RS(2,3)*RS(2,4)/1575)*RS(2,1)**3
      s10 = s11+(RS(2,2)**3/315+(RS(2,3)/525+RS(2,4)/700)*RS(2,2)**2+(RS
     #(2,3)*RS(2,4)/700+RS(2,4)**2/1050+RS(2,3)**2/1050)*RS(2,2)+RS(2,4)
     #**3/1260+RS(2,4)**2*RS(2,3)/1050+RS(2,3)**2*RS(2,4)/1400+RS(2,3)**
     #3/3150)*RS(2,1)**2+(RS(2,2)**4/252+(2.E0/1575.E0*RS(2,4)+RS(2,3)/3
     #15)*RS(2,2)**3+(RS(2,3)**2/420+RS(2,4)**2/1400+RS(2,3)*RS(2,4)/525
     #)*RS(2,2)**2+(RS(2,3)**3/630+RS(2,3)**2*RS(2,4)/525+RS(2,4)**2*RS(
     #2,3)/700+RS(2,4)**3/1575)*RS(2,2)+RS(2,3)**2*RS(2,4)**2/700+2.E0/1
     #575.E0*RS(2,3)**3*RS(2,4)+RS(2,3)**4/1260+RS(2,4)**4/1260+2.E0/157
     #5.E0*RS(2,4)**3*RS(2,3))*RS(2,1)-RS(2,2)**5/42+(-2.E0/105.E0*RS(2,
     #3)+RS(2,4)/1260)*RS(2,2)**4
      s9 = s10+(-RS(2,3)**2/70+RS(2,3)*RS(2,4)/630+RS(2,4)**2/3150)*RS(2
     #,2)**3+(-RS(2,3)**3/105+RS(2,4)**3/4200+RS(2,4)**2*RS(2,3)/1050+RS
     #(2,3)**2*RS(2,4)/420)*RS(2,2)**2+(-RS(2,3)**4/210+RS(2,3)**3*RS(2,
     #4)/315+RS(2,4)**3*RS(2,3)/1050+RS(2,3)**2*RS(2,4)**2/525+RS(2,4)**
     #4/3150)*RS(2,2)+RS(2,4)**3*RS(2,3)**2/420+RS(2,4)**4*RS(2,3)/630+R
     #S(2,3)**4*RS(2,4)/252+RS(2,4)**2*RS(2,3)**3/315+RS(2,4)**5/1260
      s10 = RS(1,3)
      s8 = s9*s10
      s12 = RS(2,1)**5/252+(RS(2,2)/315+RS(2,3)/2100+RS(2,4)/315)*RS(2,1
     #)**4+(RS(2,2)**2/420+(RS(2,4)/525+RS(2,3)/1575)*RS(2,2)+RS(2,3)*RS
     #(2,4)/1575+RS(2,3)**2/12600+RS(2,4)**2/420)*RS(2,1)**3
      s11 = s12+(RS(2,2)**3/630+(RS(2,3)/2100+RS(2,4)/1050)*RS(2,2)**2+(
     #RS(2,3)*RS(2,4)/2100+RS(2,4)**2/1050)*RS(2,2)+RS(2,4)**2*RS(2,3)/2
     #100-RS(2,3)**3/12600+RS(2,4)**3/630)*RS(2,1)**2+(RS(2,2)**4/1260+R
     #S(2,2)**3*RS(2,4)/3150+(-RS(2,3)**2/2100+RS(2,4)**2/4200)*RS(2,2)*
     #*2+(-RS(2,3)**3/1575-RS(2,3)**2*RS(2,4)/2100+RS(2,4)**3/3150)*RS(2
     #,2)-RS(2,3)**4/2100+RS(2,4)**4/1260-RS(2,3)**2*RS(2,4)**2/2100-RS(
     #2,3)**3*RS(2,4)/1575)*RS(2,1)-RS(2,2)**4*RS(2,3)/1260+(-RS(2,3)*RS
     #(2,4)/3150-RS(2,3)**2/630)*RS(2,2)**3
      s10 = s11+(-RS(2,4)**2*RS(2,3)/4200-RS(2,3)**2*RS(2,4)/1050-RS(2,3
     #)**3/420)*RS(2,2)**2+(-RS(2,3)**2*RS(2,4)**2/1050-RS(2,3)**4/315-R
     #S(2,4)**3*RS(2,3)/3150-RS(2,3)**3*RS(2,4)/525)*RS(2,2)-RS(2,3)**5/
     #252-RS(2,4)**4*RS(2,3)/1260-RS(2,4)**3*RS(2,3)**2/630-RS(2,3)**4*R
     #S(2,4)/315-RS(2,4)**2*RS(2,3)**3/420
      s11 = RS(1,4)
      s9 = s10*s11
      s7 = s8+s9
      ANTSSSS(2) = s6+s7

      s10 = (RS(2,2)/1260-RS(2,4)/1260)*RS(2,1)**4+(RS(2,2)*RS(2,3)/3150
     #-RS(2,3)*RS(2,4)/3150-RS(2,4)**2/630+RS(2,2)**2/630)*RS(2,1)**3+(R
     #S(2,2)**3/420+(RS(2,4)/2100+RS(2,3)/1050)*RS(2,2)**2+(-RS(2,4)**2/
     #2100+RS(2,3)**2/4200)*RS(2,2)-RS(2,4)**2*RS(2,3)/1050-RS(2,3)**2*R
     #S(2,4)/4200-RS(2,4)**3/420)*RS(2,1)**2
      s9 = s10+(RS(2,2)**4/315+(RS(2,4)/1575+RS(2,3)/525)*RS(2,2)**3+(RS
     #(2,3)**2/1050+RS(2,3)*RS(2,4)/2100)*RS(2,2)**2+(-RS(2,4)**3/1575+R
     #S(2,3)**3/3150-RS(2,4)**2*RS(2,3)/2100)*RS(2,2)-RS(2,3)**3*RS(2,4)
     #/3150-RS(2,3)**2*RS(2,4)**2/1050-RS(2,4)**3*RS(2,3)/525-RS(2,4)**4
     #/315)*RS(2,1)+RS(2,2)**5/252+(RS(2,4)/2100+RS(2,3)/315)*RS(2,2)**4
     #+(RS(2,3)**2/420+RS(2,3)*RS(2,4)/1575+RS(2,4)**2/12600)*RS(2,2)**3
      s8 = s9+(RS(2,3)**2*RS(2,4)/2100+RS(2,3)**3/630-RS(2,4)**3/12600)*
     #RS(2,2)**2+(RS(2,3)**4/1260-RS(2,4)**3*RS(2,3)/1575-RS(2,4)**4/210
     #0-RS(2,3)**2*RS(2,4)**2/2100)*RS(2,2)-RS(2,4)**4*RS(2,3)/315-RS(2,
     #4)**3*RS(2,3)**2/420-RS(2,4)**2*RS(2,3)**3/630-RS(2,3)**4*RS(2,4)/
     #1260-RS(2,4)**5/252
      s9 = RS(1,1)
      s7 = s8*s9
      s11 = -RS(2,1)**5/1260+(-RS(2,2)/630-RS(2,4)/1260-RS(2,3)/3150)*RS
     #(2,1)**4+(-RS(2,2)**2/420+(-RS(2,3)/1050-2.E0/1575.E0*RS(2,4))*RS(
     #2,2)-RS(2,4)**2/1260-RS(2,3)**2/4200-RS(2,3)*RS(2,4)/1575)*RS(2,1)
     #**3
      s12 = s11+(-RS(2,2)**3/315+(-RS(2,3)/525-RS(2,4)/700)*RS(2,2)**2+(
     #-RS(2,3)**2/1050-RS(2,4)**2/1050-RS(2,3)*RS(2,4)/700)*RS(2,2)-RS(2
     #,4)**3/1260-RS(2,4)**2*RS(2,3)/1050-RS(2,3)**3/3150-RS(2,3)**2*RS(
     #2,4)/1400)*RS(2,1)**2
      s10 = s12+(-RS(2,2)**4/252+(-2.E0/1575.E0*RS(2,4)-RS(2,3)/315)*RS(
     #2,2)**3+(-RS(2,3)**2/420-RS(2,3)*RS(2,4)/525-RS(2,4)**2/1400)*RS(2
     #,2)**2+(-RS(2,3)**3/630-RS(2,4)**2*RS(2,3)/700-RS(2,4)**3/1575-RS(
     #2,3)**2*RS(2,4)/525)*RS(2,2)-RS(2,3)**4/1260-2.E0/1575.E0*RS(2,3)*
     #*3*RS(2,4)-2.E0/1575.E0*RS(2,4)**3*RS(2,3)-RS(2,3)**2*RS(2,4)**2/7
     #00-RS(2,4)**4/1260)*RS(2,1)+(RS(2,3)/210-RS(2,4)/1260)*RS(2,2)**4+
     #(RS(2,3)**2/105-RS(2,3)*RS(2,4)/630-RS(2,4)**2/3150)*RS(2,2)**3
      s9 = s10+(RS(2,3)**3/70-RS(2,4)**3/4200-RS(2,4)**2*RS(2,3)/1050-RS
     #(2,3)**2*RS(2,4)/420)*RS(2,2)**2+(2.E0/105.E0*RS(2,3)**4-RS(2,4)**
     #3*RS(2,3)/1050-RS(2,3)**2*RS(2,4)**2/525-RS(2,3)**3*RS(2,4)/315-RS
     #(2,4)**4/3150)*RS(2,2)-RS(2,4)**3*RS(2,3)**2/420-RS(2,4)**4*RS(2,3
     #)/630-RS(2,4)**5/1260-RS(2,4)**2*RS(2,3)**3/315+RS(2,3)**5/42-RS(2
     #,3)**4*RS(2,4)/252
      s10 = RS(1,2)
      s8 = s9*s10
      s6 = s7+s8
      s7 = s6
      s9 = (RS(2,4)**3*RS(2,3)**2/70+RS(2,4)**5/210+RS(2,4)**4*RS(2,3)/1
     #05+RS(2,3)**4*RS(2,4)/42-2.E0/105.E0*RS(2,3)**3*RS(2,2)**2-RS(2,3)
     #**2*RS(2,2)**3/70+2.E0/105.E0*RS(2,4)**2*RS(2,3)**3-RS(2,2)**5/210
     #-RS(2,3)**4*RS(2,2)/42-RS(2,2)**4*RS(2,3)/105)*RS(1,3)
      s13 = RS(2,1)**5/1260+(RS(2,2)/1260+RS(2,4)/630+RS(2,3)/3150)*RS(2
     #,1)**4+(RS(2,2)**2/1260+(RS(2,3)/1575+2.E0/1575.E0*RS(2,4))*RS(2,2
     #)+RS(2,4)**2/420+RS(2,3)**2/4200+RS(2,3)*RS(2,4)/1050)*RS(2,1)**3
      s12 = s13+(RS(2,2)**3/1260+(RS(2,3)/1050+RS(2,4)/1050)*RS(2,2)**2+
     #(RS(2,3)*RS(2,4)/700+RS(2,4)**2/700+RS(2,3)**2/1400)*RS(2,2)+RS(2,
     #4)**3/315+RS(2,4)**2*RS(2,3)/525+RS(2,3)**2*RS(2,4)/1050+RS(2,3)**
     #3/3150)*RS(2,1)**2+(RS(2,2)**4/1260+(RS(2,4)/1575+2.E0/1575.E0*RS(
     #2,3))*RS(2,2)**3+(RS(2,3)**2/700+RS(2,4)**2/1400+RS(2,3)*RS(2,4)/7
     #00)*RS(2,2)**2+(2.E0/1575.E0*RS(2,3)**3+RS(2,3)**2*RS(2,4)/525+RS(
     #2,4)**2*RS(2,3)/525+2.E0/1575.E0*RS(2,4)**3)*RS(2,2)+RS(2,3)**2*RS
     #(2,4)**2/420+RS(2,3)**3*RS(2,4)/630+RS(2,3)**4/1260+RS(2,4)**4/252
     #+RS(2,4)**3*RS(2,3)/315)*RS(2,1)+RS(2,2)**5/1260+(RS(2,3)/630+RS(2
     #,4)/3150)*RS(2,2)**4
      s11 = s12+(RS(2,3)**2/420+RS(2,3)*RS(2,4)/1050+RS(2,4)**2/4200)*RS
     #(2,2)**3+(RS(2,3)**3/315+RS(2,4)**3/3150+RS(2,4)**2*RS(2,3)/1050+R
     #S(2,3)**2*RS(2,4)/525)*RS(2,2)**2+(RS(2,3)**4/252+RS(2,3)**3*RS(2,
     #4)/315+RS(2,4)**3*RS(2,3)/630+RS(2,3)**2*RS(2,4)**2/420+RS(2,4)**4
     #/1260)*RS(2,2)-RS(2,4)**3*RS(2,3)**2/105-RS(2,4)**4*RS(2,3)/210-2.
     #E0/105.E0*RS(2,3)**4*RS(2,4)-RS(2,4)**2*RS(2,3)**3/70-RS(2,3)**5/4
     #2
      s12 = RS(1,4)
      s10 = s11*s12
      s8 = s9+s10
      ANTSSSS(3) = s7+s8

      s11 = (RS(2,2)/252-RS(2,4)/210+RS(2,3)/1260)*RS(2,1)**4+(RS(2,2)**
     #2/315+(RS(2,4)/315+2.E0/1575.E0*RS(2,3))*RS(2,2)-RS(2,4)**2/105+RS
     #(2,3)**2/3150+RS(2,3)*RS(2,4)/630)*RS(2,1)**3+(RS(2,2)**3/420+(RS(
     #2,4)/525+RS(2,3)/700)*RS(2,2)**2+(RS(2,3)*RS(2,4)/525+RS(2,4)**2/4
     #20+RS(2,3)**2/1400)*RS(2,2)-RS(2,4)**3/70+RS(2,4)**2*RS(2,3)/420+R
     #S(2,3)**3/4200+RS(2,3)**2*RS(2,4)/1050)*RS(2,1)**2
      s10 = s11+(RS(2,2)**4/630+(RS(2,4)/1050+2.E0/1575.E0*RS(2,3))*RS(2
     #,2)**3+(RS(2,3)*RS(2,4)/700+RS(2,4)**2/1050+RS(2,3)**2/1050)*RS(2,
     #2)**2+(RS(2,3)**3/1575+RS(2,4)**3/630+RS(2,4)**2*RS(2,3)/525+RS(2,
     #3)**2*RS(2,4)/700)*RS(2,2)+RS(2,3)**3*RS(2,4)/1050+RS(2,4)**3*RS(2
     #,3)/315+RS(2,3)**4/3150+RS(2,3)**2*RS(2,4)**2/525-2.E0/105.E0*RS(2
     #,4)**4)*RS(2,1)+RS(2,2)**5/1260+(RS(2,4)/3150+RS(2,3)/1260)*RS(2,2
     #)**4+(RS(2,3)**2/1260+RS(2,4)**2/4200+RS(2,3)*RS(2,4)/1575)*RS(2,2
     #)**3
      s9 = s10+(RS(2,3)**3/1260+RS(2,3)**2*RS(2,4)/1050+RS(2,4)**3/3150+
     #RS(2,4)**2*RS(2,3)/1400)*RS(2,2)**2+(RS(2,3)**4/1260+RS(2,3)**2*RS
     #(2,4)**2/700+2.E0/1575.E0*RS(2,3)**3*RS(2,4)+2.E0/1575.E0*RS(2,4)*
     #*3*RS(2,3)+RS(2,4)**4/1260)*RS(2,2)+RS(2,4)**3*RS(2,3)**2/315+RS(2
     #,4)**4*RS(2,3)/252-RS(2,4)**5/42+RS(2,4)**2*RS(2,3)**3/420+RS(2,3)
     #**4*RS(2,4)/630+RS(2,3)**5/1260
      s10 = RS(1,1)
      s8 = s9*s10
      s12 = -RS(2,1)**5/252+(-RS(2,2)/315-RS(2,3)/2100-RS(2,4)/315)*RS(2
     #,1)**4+(-RS(2,2)**2/420+(-RS(2,4)/525-RS(2,3)/1575)*RS(2,2)-RS(2,3
     #)*RS(2,4)/1575-RS(2,3)**2/12600-RS(2,4)**2/420)*RS(2,1)**3
      s11 = s12+(-RS(2,2)**3/630+(-RS(2,3)/2100-RS(2,4)/1050)*RS(2,2)**2
     #+(-RS(2,3)*RS(2,4)/2100-RS(2,4)**2/1050)*RS(2,2)-RS(2,4)**2*RS(2,3
     #)/2100+RS(2,3)**3/12600-RS(2,4)**3/630)*RS(2,1)**2+(-RS(2,2)**4/12
     #60-RS(2,2)**3*RS(2,4)/3150+(RS(2,3)**2/2100-RS(2,4)**2/4200)*RS(2,
     #2)**2+(RS(2,3)**3/1575+RS(2,3)**2*RS(2,4)/2100-RS(2,4)**3/3150)*RS
     #(2,2)+RS(2,3)**4/2100-RS(2,4)**4/1260+RS(2,3)**2*RS(2,4)**2/2100+R
     #S(2,3)**3*RS(2,4)/1575)*RS(2,1)+RS(2,2)**4*RS(2,3)/1260+(RS(2,3)*R
     #S(2,4)/3150+RS(2,3)**2/630)*RS(2,2)**3
      s10 = s11+(RS(2,4)**2*RS(2,3)/4200+RS(2,3)**2*RS(2,4)/1050+RS(2,3)
     #**3/420)*RS(2,2)**2+(RS(2,3)**2*RS(2,4)**2/1050+RS(2,3)**4/315+RS(
     #2,4)**3*RS(2,3)/3150+RS(2,3)**3*RS(2,4)/525)*RS(2,2)+RS(2,3)**5/25
     #2+RS(2,4)**4*RS(2,3)/1260+RS(2,4)**3*RS(2,3)**2/630+RS(2,3)**4*RS(
     #2,4)/315+RS(2,4)**2*RS(2,3)**3/420
      s11 = RS(1,2)
      s9 = s10*s11
      s7 = s8+s9
      s8 = s7
      s13 = -RS(2,1)**5/1260+(-RS(2,2)/1260-RS(2,4)/630-RS(2,3)/3150)*RS
     #(2,1)**4+(-RS(2,2)**2/1260+(-RS(2,3)/1575-2.E0/1575.E0*RS(2,4))*RS
     #(2,2)-RS(2,4)**2/420-RS(2,3)**2/4200-RS(2,3)*RS(2,4)/1050)*RS(2,1)
     #**3
      s12 = s13+(-RS(2,2)**3/1260+(-RS(2,3)/1050-RS(2,4)/1050)*RS(2,2)**
     #2+(-RS(2,3)*RS(2,4)/700-RS(2,4)**2/700-RS(2,3)**2/1400)*RS(2,2)-RS
     #(2,4)**3/315-RS(2,4)**2*RS(2,3)/525-RS(2,3)**2*RS(2,4)/1050-RS(2,3
     #)**3/3150)*RS(2,1)**2+(-RS(2,2)**4/1260+(-RS(2,4)/1575-2.E0/1575.E
     #0*RS(2,3))*RS(2,2)**3+(-RS(2,3)**2/700-RS(2,4)**2/1400-RS(2,3)*RS(
     #2,4)/700)*RS(2,2)**2+(-2.E0/1575.E0*RS(2,3)**3-RS(2,3)**2*RS(2,4)/
     #525-RS(2,4)**2*RS(2,3)/525-2.E0/1575.E0*RS(2,4)**3)*RS(2,2)-RS(2,3
     #)**2*RS(2,4)**2/420-RS(2,3)**3*RS(2,4)/630-RS(2,3)**4/1260-RS(2,4)
     #**4/252-RS(2,4)**3*RS(2,3)/315)*RS(2,1)-RS(2,2)**5/1260+(-RS(2,3)/
     #630-RS(2,4)/3150)*RS(2,2)**4
      s11 = s12+(-RS(2,3)**2/420-RS(2,3)*RS(2,4)/1050-RS(2,4)**2/4200)*R
     #S(2,2)**3+(-RS(2,3)**3/315-RS(2,4)**3/3150-RS(2,4)**2*RS(2,3)/1050
     #-RS(2,3)**2*RS(2,4)/525)*RS(2,2)**2+(-RS(2,3)**4/252-RS(2,3)**3*RS
     #(2,4)/315-RS(2,4)**3*RS(2,3)/630-RS(2,3)**2*RS(2,4)**2/420-RS(2,4)
     #**4/1260)*RS(2,2)+RS(2,4)**3*RS(2,3)**2/70+2.E0/105.E0*RS(2,4)**4*
     #RS(2,3)+RS(2,3)**4*RS(2,4)/210+RS(2,4)**2*RS(2,3)**3/105+RS(2,4)**
     #5/42
      s12 = RS(1,3)
      s10 = s11*s12
      s11 = (-2.E0/105.E0*RS(2,4)**3*RS(2,3)**2-RS(2,4)**2*RS(2,3)**3/70
     #+RS(2,4)**4*RS(2,1)/42+RS(2,4)**2*RS(2,1)**3/70+RS(2,1)**5/210+RS(
     #2,1)**4*RS(2,4)/105-RS(2,3)**4*RS(2,4)/105+2.E0/105.E0*RS(2,4)**3*
     #RS(2,1)**2-RS(2,4)**4*RS(2,3)/42-RS(2,3)**5/210)*RS(1,4)
      s9 = s10+s11
      ANTSSSS(4) = s8+s9

      RETURN
      END
C
C=====================================================================
      SUBROUTINE MEMINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C	------------------------------------------------------
C     COORD(2,NNO)  = CURRENT NODALCOORDINATES
C     COORDI(2,NNO) = INITIAL NODAL COORDINATES
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     REDIS(48)     = COROTATIONAL FORM OF EDIS
C     NNO           = NUMBER OF NODES
C     NEF           = NUMBER OF ELEMENT FREEDOMS
C     ------------------------------------------------------
C
      DIMENSION COORD(2,8),COORDI(2,8),EDIS(1),REDIS(48)
C
      DO 20 I=1,NEF
   20 REDIS(I)=EDIS(I)

      K=1
      DO 40 I=1,NNO
       DO 30 J=1,2
         COORDI(J,I)=COORD(J,I)-EDIS(K)
30	 K=K+1
40	K=K+1

      RETURN
      END
C
C=====================================================================
      SUBROUTINE MEMDSP(COORD,COORDI,EDIS,VR,VS,REDIS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     PURPOSE:	MODIFIES TOTAL DISPLACEMENT VECTOR BY DEDUCTING
C				RIGID BODY TRANSLATIONS AND ROTATIONS
C
C	INPUT VARIABLES
C     COORD(2,NNO)      = CURRENT NODAL COORDINATES
C     COORDI(2,NNO)     = INITIAL NODAL COORDINATES
C     EDIS(NEF)         = CURRENT NODAL DISPLACEMENTS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     VR(3),VS(3),VT(3) = CURRENT DIRECTION COSINE VECTORS
C
C	LOCAL VARIABLES
C	ARC				  = 2*PI
C	DUM(2,4)		  = DUMMY VARIABLE FOR INITIAL LOCAL COORD.
C
C	OUTPUT VARIABLES
C     REDIS(48)         = COROTATIONAL FORM OF EDIS
C     ----------------------------------------------------
      DIMENSION COORD(2,8),COORDI(2,8),EDIS(*),VR(2),VS(2)
      DIMENSION TM(3,3),XJI(4),CDO(2,8),DUM(2,4),XYZ(2)
      DIMENSION VRO(2),VSO(2)
	DIMENSION REDIS(48)
C
      ARC=6.2831853071796

C     ---------------------------------------------------------------
C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
C     ---------------------------------------------------------------
C	DETERMINE ROTATIONAL COMPONENTS A THE ELEMENT CENTER/REF POINT
      ROT=(EDIS(3)+EDIS(6)+EDIS(9)+EDIS(12))/4
C     FIND ROTATIONAL PSEUDOVECTOR (ROV)
	ROV = 1.0
C	RESULTANT ROTATION (RO)
 	RO  = ROT
      IF (RO.EQ.0.) RETURN

C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
C	ELEMENT INITIAL CENTER COORD/INITIAL REF PT COORD - XYZ
	XYZ(1)=(COORDI(1,1)+COORDI(1,2)+COORDI(1,3)+COORDI(1,4))/4
	XYZ(2)=(COORDI(2,1)+COORDI(2,2)+COORDI(2,3)+COORDI(2,4))/4

C	INITIAL NODAL DISTANCE RELATIVE TO INITIAL REF POINT - CDO
      DO 150 I=1,NNO
      DO 150 J=1,2
150	CDO(J,I)=COORDI(J,I)-XYZ(J)

C	FIND INITIAL COORDINATE DIRECTION VECTORS - VRO,VSO,VTO
	CALL VRS(COORDI,VRO,VSO,DUM,NNO)

C	DETERMINE ORTHOGONAL ROTATION MATRIX - TM
      DO 160 I=1,2
       DO 160 J=1,2
160	  TM(I,J) = VR(I)*VRO(J)+VS(I)*VSO(J)

	VV = ROV

C	DETERMINE ANGLE OF ROTATION - RR
      ABCD=(.5*(TM(1,1)+TM(2,2)))
	IF (ABCD .GT.  1.00000000000000000 ) ABCD =  1.000000000000000000
	IF (ABCD .LT. -1.00000000000000000 ) ABCD = -1.000000000000000000
      RR = ACOS(ABCD)
      SN = SIN(RR)
      IF (SN.EQ.0.0) GOTO 190

C	DETERMINE ROTATION AXIS VECTOR - VV
      F=.5/SN
      VV=F*(TM(2,1)-TM(1,2))
	CS = ROV*VV
      IF (CS.GE.0.) GOTO 190
      RR  = -RR
	VV3 = -VV3

C	???
190	RD=RO-RR
      N=(RD+SIGN(ARC/2.01,RD))/ARC
      RR=RR+N*ARC

C	DETERMINE CO-ROTATIONAL DISPLACEMENTS
      K=1
      DO 220 I=1,NNO
        DO 210 J=1,2
          TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)
          REDIS(K)=COORD(J,I)-TCDO
210	    K=K+1
        REDIS(K)=EDIS(K)-RR*VV
220	K=K+1


      RETURN
      END
C
C=====================================================================
      SUBROUTINE KGMEM(ST,TT,S)

	IMPLICIT REAL*8 (A-H,O-Z)
C	------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR THE GLOBAL STIFFNESS MATRIX T*KL*TT
C
C	INPUT VARIABLES
C	ST(300) = ELEMENT STIFFNESS MATRIX IN LOCAL COORDINATES
C	TT(24,24) = TRANSPOSE OF THE TRANSFORMATION MATRIX
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIAVBLES
C	S(*)		= GLOBAL STIFFNESS MATRIX IN ROW FORMAT
C	------------------------------------------------------------
	DIMENSION ST(78),TT(12,12)
	DIMENSION S(*)

      t3 = TT(1,1)*ST(1)+TT(2,1)*ST(2)
      t7 = TT(1,1)*ST(2)+TT(2,1)*ST(13)
      t18 = TT(1,1)*ST(4)+TT(2,1)*ST(15)
      t22 = TT(1,1)*ST(5)+TT(2,1)*ST(16)
      t33 = TT(1,1)*ST(7)+TT(2,1)*ST(18)
      t37 = TT(1,1)*ST(8)+TT(2,1)*ST(19)
      t48 = TT(1,1)*ST(10)+TT(2,1)*ST(21)
      t52 = TT(1,1)*ST(11)+TT(2,1)*ST(22)
      t75 = TT(1,2)*ST(4)+TT(2,2)*ST(15)
      t79 = TT(1,2)*ST(5)+TT(2,2)*ST(16)
      t90 = TT(1,2)*ST(7)+TT(2,2)*ST(18)
      t94 = TT(1,2)*ST(8)+TT(2,2)*ST(19)
      t105 = TT(1,2)*ST(10)+TT(2,2)*ST(21)
      t109 = TT(1,2)*ST(11)+TT(2,2)*ST(22)
      t138 = TT(1,1)*ST(34)+TT(2,1)*ST(35)
      t142 = TT(1,1)*ST(35)+TT(2,1)*ST(43)
      t153 = TT(1,1)*ST(37)+TT(2,1)*ST(45)
      t157 = TT(1,1)*ST(38)+TT(2,1)*ST(46)
      t168 = TT(1,1)*ST(40)+TT(2,1)*ST(48)
      t172 = TT(1,1)*ST(41)+TT(2,1)*ST(49)
      t195 = TT(1,2)*ST(37)+TT(2,2)*ST(45)
      t199 = TT(1,2)*ST(38)+TT(2,2)*ST(46)
      t210 = TT(1,2)*ST(40)+TT(2,2)*ST(48)
      t214 = TT(1,2)*ST(41)+TT(2,2)*ST(49)
      t237 = TT(1,1)*ST(58)+TT(2,1)*ST(59)
      t241 = TT(1,1)*ST(59)+TT(2,1)*ST(64)
      t252 = TT(1,1)*ST(61)+TT(2,1)*ST(66)
      t256 = TT(1,1)*ST(62)+TT(2,1)*ST(67)
      t279 = TT(1,2)*ST(61)+TT(2,2)*ST(66)
      t283 = TT(1,2)*ST(62)+TT(2,2)*ST(67)
      t300 = TT(1,1)*ST(73)+TT(2,1)*ST(74)
      t304 = TT(1,1)*ST(74)+TT(2,1)*ST(76)
      S(1) = t3*TT(1,1)+t7*TT(2,1)
      S(2) = t3*TT(1,2)+t7*TT(2,2)
      S(3) = TT(1,1)*ST(3)+TT(2,1)*ST(14)
      S(4) = t18*TT(1,1)+t22*TT(2,1)
      S(5) = t18*TT(1,2)+t22*TT(2,2)
      S(6) = TT(1,1)*ST(6)+TT(2,1)*ST(17)
      S(7) = t33*TT(1,1)+t37*TT(2,1)
      S(8) = t33*TT(1,2)+t37*TT(2,2)
      S(9) = TT(1,1)*ST(9)+TT(2,1)*ST(20)
      S(10) = t48*TT(1,1)+t52*TT(2,1)
      S(11) = t48*TT(1,2)+t52*TT(2,2)
      S(12) = TT(1,1)*ST(12)+TT(2,1)*ST(23)
      S(13) = (TT(1,2)*ST(1)+TT(2,2)*ST(2))*TT(1,2)+(TT(1,2)*ST(2)+TT(2,
     #2)*ST(13))*TT(2,2)
      S(14) = TT(1,2)*ST(3)+TT(2,2)*ST(14)
      S(15) = t75*TT(1,1)+t79*TT(2,1)
      S(16) = t75*TT(1,2)+t79*TT(2,2)
      S(17) = TT(1,2)*ST(6)+TT(2,2)*ST(17)
      S(18) = t90*TT(1,1)+t94*TT(2,1)
      S(19) = t90*TT(1,2)+t94*TT(2,2)
      S(20) = TT(1,2)*ST(9)+TT(2,2)*ST(20)
      S(21) = t105*TT(1,1)+t109*TT(2,1)
      S(22) = t105*TT(1,2)+t109*TT(2,2)
      S(23) = TT(1,2)*ST(12)+TT(2,2)*ST(23)
      S(24) = ST(24)
      S(25) = TT(1,1)*ST(25)+TT(2,1)*ST(26)
      S(26) = TT(1,2)*ST(25)+TT(2,2)*ST(26)
      S(27) = ST(27)
      S(28) = TT(1,1)*ST(28)+TT(2,1)*ST(29)
      S(29) = TT(1,2)*ST(28)+TT(2,2)*ST(29)
      S(30) = ST(30)
      S(31) = TT(1,1)*ST(31)+TT(2,1)*ST(32)
      S(32) = TT(1,2)*ST(31)+TT(2,2)*ST(32)
      S(33) = ST(33)
      S(34) = t138*TT(1,1)+t142*TT(2,1)
      S(35) = t138*TT(1,2)+t142*TT(2,2)
      S(36) = TT(1,1)*ST(36)+TT(2,1)*ST(44)
      S(37) = t153*TT(1,1)+t157*TT(2,1)
      S(38) = t153*TT(1,2)+t157*TT(2,2)
      S(39) = TT(1,1)*ST(39)+TT(2,1)*ST(47)
      S(40) = t168*TT(1,1)+t172*TT(2,1)
      S(41) = t168*TT(1,2)+t172*TT(2,2)
      S(42) = TT(1,1)*ST(42)+TT(2,1)*ST(50)
      S(43) = (TT(1,2)*ST(34)+TT(2,2)*ST(35))*TT(1,2)+(TT(1,2)*ST(35)+TT
     #(2,2)*ST(43))*TT(2,2)
      S(44) = TT(1,2)*ST(36)+TT(2,2)*ST(44)
      S(45) = t195*TT(1,1)+t199*TT(2,1)
      S(46) = t195*TT(1,2)+t199*TT(2,2)
      S(47) = TT(1,2)*ST(39)+TT(2,2)*ST(47)
      S(48) = t210*TT(1,1)+t214*TT(2,1)
      S(49) = t210*TT(1,2)+t214*TT(2,2)
      S(50) = TT(1,2)*ST(42)+TT(2,2)*ST(50)
      S(51) = ST(51)
      S(52) = TT(1,1)*ST(52)+TT(2,1)*ST(53)
      S(53) = TT(1,2)*ST(52)+TT(2,2)*ST(53)
      S(54) = ST(54)
      S(55) = TT(1,1)*ST(55)+TT(2,1)*ST(56)
      S(56) = TT(1,2)*ST(55)+TT(2,2)*ST(56)
      S(57) = ST(57)
      S(58) = t237*TT(1,1)+t241*TT(2,1)
      S(59) = t237*TT(1,2)+t241*TT(2,2)
      S(60) = TT(1,1)*ST(60)+TT(2,1)*ST(65)
      S(61) = t252*TT(1,1)+t256*TT(2,1)
      S(62) = t252*TT(1,2)+t256*TT(2,2)
      S(63) = TT(1,1)*ST(63)+TT(2,1)*ST(68)
      S(64) = (TT(1,2)*ST(58)+TT(2,2)*ST(59))*TT(1,2)+(TT(1,2)*ST(59)+TT
     #(2,2)*ST(64))*TT(2,2)
      S(65) = TT(1,2)*ST(60)+TT(2,2)*ST(65)
      S(66) = t279*TT(1,1)+t283*TT(2,1)
      S(67) = t279*TT(1,2)+t283*TT(2,2)
      S(68) = TT(1,2)*ST(63)+TT(2,2)*ST(68)
      S(69) = ST(69)
      S(70) = TT(1,1)*ST(70)+TT(2,1)*ST(71)
      S(71) = TT(1,2)*ST(70)+TT(2,2)*ST(71)
      S(72) = ST(72)
      S(73) = t300*TT(1,1)+t304*TT(2,1)
      S(74) = t300*TT(1,2)+t304*TT(2,2)
      S(75) = TT(1,1)*ST(75)+TT(2,1)*ST(77)
      S(76) = (TT(1,2)*ST(73)+TT(2,2)*ST(74))*TT(1,2)+(TT(1,2)*ST(74)+TT
     #(2,2)*ST(76))*TT(2,2)
      S(77) = TT(1,2)*ST(75)+TT(2,2)*ST(77)
      S(78) = ST(78)

      RETURN
      END
C
C=====================================================================
      SUBROUTINE COMF(RS,AX,AR,AS,ARR,ASS,ARS,ARRR,ASSS,ARRS,ASSR,
	1                ARRRR,ASSSS,ARRSS)
	IMPLICIT REAL*8 (A-H,O-Z)
C     ----------------------------------------------------------
C     PURPOSE:	COMPUTES FOR INTEGRATION FACTORS FOR INTEGRALS
C				USING ISOPARAMETRIC MAPPING
C
C	INPUT VARIABLES
C	RS(2,4)  = NODAL COORDINATES
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLES
C     AX(4)	= FACTORS FOR INT(PHI)
C     AR(4)	= FACTORS FOR INT(R*PHI)
C     AS(4)	= FACTORS FOR INT(S*PHI)
C     ARR(4)	= FACTORS FOR INT(R*R*PHI)
C     ASS(4)	= FACTORS FOR INT(S*S*PHI)
C     ARS(4)	= FACTORS FOR INT(R*S*PHI)
C     ARRR(4)	= FACTORS FOR INT(R*R*R*PHI)
C     ASSS(4)	= FACTORS FOR INT(S*S*S*PHI)
C     ARRS(4)	= FACTORS FOR INT(R*R*S*PHI)
C     ASSR(4)	= FACTORS FOR INT(S*S*R*PHI)
C     ARRRR(4)= FACTORS FOR INT(R*R*R*R*PHI)
C     ASSSS(4)= FACTORS FOR INT(S*S*S*S*PHI)
C     ARRSS(4)= FACTORS FOR INT(R*R*S*S*PHI)
C     ----------------------------------------------------------
	DIMENSION RS(2,4)
	DIMENSION AX(5),AR(5),AS(5),ARR(5),ASS(5),ARS(5)
	DIMENSION ARRR(5),ASSS(5),ARRS(5),ASSR(5)
	DIMENSION ARRRR(5),ASSSS(5),ARRSS(5)

C	---------------------
C	COMPUTE FOR FACTOR AX
C	---------------------
      t0 = (RS(2,2)/6-RS(2,4)/6)*RS(1,1)+(-RS(2,1)/6-RS(2,3)/36-RS(2,4)/
     #36)*RS(1,2)+(RS(2,2)/36-RS(2,4)/36)*RS(1,3)+(RS(2,1)/6+RS(2,3)/36+
     #RS(2,2)/36)*RS(1,4)
	AX(1) = T0

      t0 = (RS(2,2)/6+RS(2,3)/36+RS(2,4)/36)*RS(1,1)+(RS(2,3)/6-RS(2,1)/
     #6)*RS(1,2)+(-RS(2,4)/36-RS(2,1)/36-RS(2,2)/6)*RS(1,3)+(-RS(2,1)/36
     #+RS(2,3)/36)*RS(1,4)
	AX(2) = T0

      t0 = (-RS(2,2)/36+RS(2,4)/36)*RS(1,1)+(RS(2,4)/36+RS(2,1)/36+RS(2,
     #3)/6)*RS(1,2)+(RS(2,4)/6-RS(2,2)/6)*RS(1,3)+(-RS(2,2)/36-RS(2,1)/3
     #6-RS(2,3)/6)*RS(1,4)
	AX(3) = T0

      t0 = (-RS(2,2)/36-RS(2,3)/36-RS(2,4)/6)*RS(1,1)+(-RS(2,3)/36+RS(2,
     #1)/36)*RS(1,2)+(RS(2,4)/6+RS(2,1)/36+RS(2,2)/36)*RS(1,3)+(RS(2,1)/
     #6-RS(2,3)/6)*RS(1,4)
	AX(4) = T0

      t0 = (2.E0/9.E0*RS(2,2)-2.E0/9.E0*RS(2,4))*RS(1,1)+(-2.E0/9.E0*RS(
     #2,1)+2.E0/9.E0*RS(2,3))*RS(1,2)+(-2.E0/9.E0*RS(2,2)+2.E0/9.E0*RS(2
     #,4))*RS(1,3)+(2.E0/9.E0*RS(2,1)-2.E0/9.E0*RS(2,3))*RS(1,4)
	AX(5) = T0

C	---------------------
C	COMPUTE FOR FACTOR AR
C	---------------------
      t0 = (RS(2,2)/12-RS(2,4)/12)*RS(1,1)**2+((-RS(2,1)/12+RS(2,2)/24-R
     #S(2,3)/200-RS(2,4)/450)*RS(1,2)+(RS(2,2)/200-RS(2,4)/200)*RS(1,3)+
     #(RS(2,1)/12-RS(2,4)/24+RS(2,3)/200+RS(2,2)/450)*RS(1,4))*RS(1,1)+(
     #-RS(2,3)/450-RS(2,1)/24-RS(2,4)/200)*RS(1,2)**2+((-RS(2,3)/200-RS(
     #2,4)/600+RS(2,2)/450)*RS(1,3)+(RS(2,2)/200-RS(2,4)/200)*RS(1,4))*R
     #S(1,2)+(RS(2,2)/200-RS(2,4)/200)*RS(1,3)**2+(RS(2,3)/200-RS(2,4)/4
     #50+RS(2,2)/600)*RS(1,4)*RS(1,3)+(RS(2,2)/200+RS(2,3)/450+RS(2,1)/2
     #4)*RS(1,4)**2
	AR(1) = T0

      t0 = (RS(2,4)/450+RS(2,2)/24+RS(2,3)/200)*RS(1,1)**2+((RS(2,4)/200
     #+RS(2,3)/450-RS(2,1)/24+RS(2,2)/12)*RS(1,2)+(RS(2,3)/200-RS(2,1)/2
     #00)*RS(1,3)+(RS(2,3)/600-RS(2,1)/450+RS(2,4)/200)*RS(1,4))*RS(1,1)
     #+(RS(2,3)/12-RS(2,1)/12)*RS(1,2)**2+((-RS(2,1)/450-RS(2,2)/12-RS(2
     #,4)/200+RS(2,3)/24)*RS(1,3)+(RS(2,3)/200-RS(2,1)/200)*RS(1,4))*RS(
     #1,2)+(-RS(2,1)/200-RS(2,2)/24-RS(2,4)/450)*RS(1,3)**2+(-RS(2,4)/20
     #0-RS(2,1)/600+RS(2,3)/450)*RS(1,4)*RS(1,3)+(RS(2,3)/200-RS(2,1)/20
     #0)*RS(1,4)**2
	AR(2) = T0

      t0 = (-RS(2,2)/200+RS(2,4)/200)*RS(1,1)**2+((RS(2,4)/600-RS(2,2)/4
     #50+RS(2,1)/200)*RS(1,2)+(-RS(2,2)/200+RS(2,4)/200)*RS(1,3)+(-RS(2,
     #1)/200+RS(2,4)/450-RS(2,2)/600)*RS(1,4))*RS(1,1)+(RS(2,1)/450+RS(2
     #,3)/24+RS(2,4)/200)*RS(1,2)**2+((RS(2,4)/450-RS(2,2)/24+RS(2,3)/12
     #+RS(2,1)/200)*RS(1,3)+(-RS(2,2)/200+RS(2,4)/200)*RS(1,4))*RS(1,2)+
     #(-RS(2,2)/12+RS(2,4)/12)*RS(1,3)**2+(-RS(2,2)/450-RS(2,1)/200-RS(2
     #,3)/12+RS(2,4)/24)*RS(1,4)*RS(1,3)+(-RS(2,1)/450-RS(2,2)/200-RS(2,
     #3)/24)*RS(1,4)**2
	AR(3) = T0

      t0 = (-RS(2,4)/24-RS(2,3)/200-RS(2,2)/450)*RS(1,1)**2+((-RS(2,2)/2
     #00-RS(2,3)/600+RS(2,1)/450)*RS(1,2)+(RS(2,1)/200-RS(2,3)/200)*RS(1
     #,3)+(-RS(2,3)/450-RS(2,2)/200-RS(2,4)/12+RS(2,1)/24)*RS(1,4))*RS(1
     #,1)+(RS(2,1)/200-RS(2,3)/200)*RS(1,2)**2+((-RS(2,3)/450+RS(2,2)/20
     #0+RS(2,1)/600)*RS(1,3)+(RS(2,1)/200-RS(2,3)/200)*RS(1,4))*RS(1,2)+
     #(RS(2,1)/200+RS(2,2)/450+RS(2,4)/24)*RS(1,3)**2+(-RS(2,3)/24+RS(2,
     #1)/450+RS(2,4)/12+RS(2,2)/200)*RS(1,4)*RS(1,3)+(RS(2,1)/12-RS(2,3)
     #/12)*RS(1,4)**2
	AR(4) = T0

      t0 = (-11.E0/225.E0*RS(2,4)+11.E0/225.E0*RS(2,2))*RS(1,1)**2+((RS(
     #2,3)/225-RS(2,4)/225-11.E0/225.E0*RS(2,1)+11.E0/225.E0*RS(2,2))*RS
     #(1,2)+(11.E0/225.E0*RS(2,1)-11.E0/225.E0*RS(2,4)+RS(2,2)/225-RS(2,
     #3)/225)*RS(1,4))*RS(1,1)+(11.E0/225.E0*RS(2,3)-11.E0/225.E0*RS(2,1
     #))*RS(1,2)**2+(-RS(2,1)/225+11.E0/225.E0*RS(2,3)-11.E0/225.E0*RS(2
     #,2)+RS(2,4)/225)*RS(1,3)*RS(1,2)+(11.E0/225.E0*RS(2,4)-11.E0/225.E
     #0*RS(2,2))*RS(1,3)**2+(11.E0/225.E0*RS(2,4)+RS(2,1)/225-11.E0/225.
     #E0*RS(2,3)-RS(2,2)/225)*RS(1,4)*RS(1,3)+(11.E0/225.E0*RS(2,1)-11.E
     #0/225.E0*RS(2,3))*RS(1,4)**2
	AR(5) = T0

C	---------------------
C	COMPUTE FOR FACTOR AS
C	---------------------
      t0 = (-RS(1,2)/12+RS(1,4)/12)*RS(2,1)**2+((RS(1,4)/450-RS(1,2)/24+
     #RS(1,1)/12+RS(1,3)/200)*RS(2,2)+(RS(1,4)/200-RS(1,2)/200)*RS(2,3)+
     #(-RS(1,1)/12+RS(1,4)/24-RS(1,2)/450-RS(1,3)/200)*RS(2,4))*RS(2,1)+
     #(RS(1,4)/200+RS(1,3)/450+RS(1,1)/24)*RS(2,2)**2+((RS(1,4)/600+RS(1
     #,3)/200-RS(1,2)/450)*RS(2,3)+(RS(1,4)/200-RS(1,2)/200)*RS(2,4))*RS
     #(2,2)+(RS(1,4)/200-RS(1,2)/200)*RS(2,3)**2+(-RS(1,3)/200-RS(1,2)/6
     #00+RS(1,4)/450)*RS(2,4)*RS(2,3)+(-RS(1,3)/450-RS(1,1)/24-RS(1,2)/2
     #00)*RS(2,4)**2
	AS(1) = T0

      t0 = (-RS(1,2)/24-RS(1,4)/450-RS(1,3)/200)*RS(2,1)**2+((-RS(1,2)/1
     #2+RS(1,1)/24-RS(1,4)/200-RS(1,3)/450)*RS(2,2)+(-RS(1,3)/200+RS(1,1
     #)/200)*RS(2,3)+(-RS(1,3)/600+RS(1,1)/450-RS(1,4)/200)*RS(2,4))*RS(
     #2,1)+(-RS(1,3)/12+RS(1,1)/12)*RS(2,2)**2+((-RS(1,3)/24+RS(1,4)/200
     #+RS(1,1)/450+RS(1,2)/12)*RS(2,3)+(-RS(1,3)/200+RS(1,1)/200)*RS(2,4
     #))*RS(2,2)+(RS(1,2)/24+RS(1,4)/450+RS(1,1)/200)*RS(2,3)**2+(RS(1,1
     #)/600-RS(1,3)/450+RS(1,4)/200)*RS(2,4)*RS(2,3)+(-RS(1,3)/200+RS(1,
     #1)/200)*RS(2,4)**2
	AS(2) = T0

      t0 = (RS(1,2)/200-RS(1,4)/200)*RS(2,1)**2+((RS(1,2)/450-RS(1,4)/60
     #0-RS(1,1)/200)*RS(2,2)+(RS(1,2)/200-RS(1,4)/200)*RS(2,3)+(RS(1,2)/
     #600+RS(1,1)/200-RS(1,4)/450)*RS(2,4))*RS(2,1)+(-RS(1,3)/24-RS(1,4)
     #/200-RS(1,1)/450)*RS(2,2)**2+((-RS(1,3)/12-RS(1,4)/450-RS(1,1)/200
     #+RS(1,2)/24)*RS(2,3)+(RS(1,2)/200-RS(1,4)/200)*RS(2,4))*RS(2,2)+(R
     #S(1,2)/12-RS(1,4)/12)*RS(2,3)**2+(RS(1,2)/450+RS(1,3)/12+RS(1,1)/2
     #00-RS(1,4)/24)*RS(2,4)*RS(2,3)+(RS(1,1)/450+RS(1,3)/24+RS(1,2)/200
     #)*RS(2,4)**2
	AS(3) = T0

      t0 = (RS(1,4)/24+RS(1,2)/450+RS(1,3)/200)*RS(2,1)**2+((RS(1,2)/200
     #-RS(1,1)/450+RS(1,3)/600)*RS(2,2)+(RS(1,3)/200-RS(1,1)/200)*RS(2,3
     #)+(-RS(1,1)/24+RS(1,2)/200+RS(1,3)/450+RS(1,4)/12)*RS(2,4))*RS(2,1
     #)+(RS(1,3)/200-RS(1,1)/200)*RS(2,2)**2+((RS(1,3)/450-RS(1,2)/200-R
     #S(1,1)/600)*RS(2,3)+(RS(1,3)/200-RS(1,1)/200)*RS(2,4))*RS(2,2)+(-R
     #S(1,1)/200-RS(1,2)/450-RS(1,4)/24)*RS(2,3)**2+(-RS(1,4)/12+RS(1,3)
     #/24-RS(1,2)/200-RS(1,1)/450)*RS(2,4)*RS(2,3)+(RS(1,3)/12-RS(1,1)/1
     #2)*RS(2,4)**2
	AS(4) = T0

      t0 = (-11.E0/225.E0*RS(1,2)+11.E0/225.E0*RS(1,4))*RS(2,1)**2+((RS(
     #1,4)/225-RS(1,3)/225+11.E0/225.E0*RS(1,1)-11.E0/225.E0*RS(1,2))*RS
     #(2,2)+(11.E0/225.E0*RS(1,4)-11.E0/225.E0*RS(1,1)-RS(1,2)/225+RS(1,
     #3)/225)*RS(2,4))*RS(2,1)+(-11.E0/225.E0*RS(1,3)+11.E0/225.E0*RS(1,
     #1))*RS(2,2)**2+(-11.E0/225.E0*RS(1,3)-RS(1,4)/225+11.E0/225.E0*RS(
     #1,2)+RS(1,1)/225)*RS(2,3)*RS(2,2)+(-11.E0/225.E0*RS(1,4)+11.E0/225
     #.E0*RS(1,2))*RS(2,3)**2+(-RS(1,1)/225+RS(1,2)/225-11.E0/225.E0*RS(
     #1,4)+11.E0/225.E0*RS(1,3))*RS(2,4)*RS(2,3)+(-11.E0/225.E0*RS(1,1)+
     #11.E0/225.E0*RS(1,3))*RS(2,4)**2
	AS(5) = T0

C	----------------------
C	COMPUTE FOR FACTOR ARR
C	----------------------
      s2 = (RS(2,2)/20-RS(2,4)/20)*RS(1,1)**3+((-437.E0/176400.E0*RS(2,3
     #)-RS(2,1)/20+RS(2,2)/30+9.E0/19600.E0*RS(2,4))*RS(1,2)+(-437.E0/17
     #6400.E0*RS(2,4)+437.E0/176400.E0*RS(2,2))*RS(1,3)+(-RS(2,4)/30+437
     #.E0/176400.E0*RS(2,3)-9.E0/19600.E0*RS(2,2)+RS(2,1)/20)*RS(1,4))*R
     #S(1,1)**2
      s1 = s2+((RS(2,2)/60-43.E0/29400.E0*RS(2,4)-RS(2,1)/30-43.E0/29400
     #.E0*RS(2,3))*RS(1,2)**2+((53.E0/88200.E0*RS(2,3)+43.E0/29400.E0*RS
     #(2,2)-RS(2,4)/14700)*RS(1,3)+(43.E0/29400.E0*RS(2,2)-43.E0/29400.E
     #0*RS(2,4))*RS(1,4))*RS(1,2)+(-53.E0/88200.E0*RS(2,2)+53.E0/88200.E
     #0*RS(2,4))*RS(1,3)**2+(-53.E0/88200.E0*RS(2,3)+RS(2,2)/14700-43.E0
     #/29400.E0*RS(2,4))*RS(1,4)*RS(1,3)+(-RS(2,4)/60+43.E0/29400.E0*RS(
     #2,2)+43.E0/29400.E0*RS(2,3)+RS(2,1)/30)*RS(1,4)**2)*RS(1,1)+(-437.
     #E0/176400.E0*RS(2,4)-RS(2,1)/60+9.E0/19600.E0*RS(2,3))*RS(1,2)**3+
     #((-43.E0/29400.E0*RS(2,3)-31.E0/58800.E0*RS(2,4)-9.E0/19600.E0*RS(
     #2,2))*RS(1,3)+(53.E0/88200.E0*RS(2,4)+437.E0/176400.E0*RS(2,2)+9.E
     #0/19600.E0*RS(2,3))*RS(1,4))*RS(1,2)**2
      t0 = s1+((-31.E0/58800.E0*RS(2,4)+43.E0/29400.E0*RS(2,2)-437.E0/17
     #6400.E0*RS(2,3))*RS(1,3)**2+(-RS(2,4)/14700+RS(2,2)/14700)*RS(1,4)
     #*RS(1,3)+(-9.E0/19600.E0*RS(2,3)-53.E0/88200.E0*RS(2,2)-437.E0/176
     #400.E0*RS(2,4))*RS(1,4)**2)*RS(1,2)+(-437.E0/176400.E0*RS(2,4)+437
     #.E0/176400.E0*RS(2,2))*RS(1,3)**3+(-43.E0/29400.E0*RS(2,4)+437.E0/
     #176400.E0*RS(2,3)+31.E0/58800.E0*RS(2,2))*RS(1,4)*RS(1,3)**2+(31.E
     #0/58800.E0*RS(2,2)+9.E0/19600.E0*RS(2,4)+43.E0/29400.E0*RS(2,3))*R
     #S(1,4)**2*RS(1,3)+(RS(2,1)/60+437.E0/176400.E0*RS(2,2)-9.E0/19600.
     #E0*RS(2,3))*RS(1,4)**3
	ARR(1) = T0

      s2 = (437.E0/176400.E0*RS(2,3)+RS(2,2)/60-9.E0/19600.E0*RS(2,4))*R
     #S(1,1)**3+((43.E0/29400.E0*RS(2,4)+RS(2,2)/30-RS(2,1)/60+43.E0/294
     #00.E0*RS(2,3))*RS(1,2)+(-9.E0/19600.E0*RS(2,4)-437.E0/176400.E0*RS
     #(2,1)-53.E0/88200.E0*RS(2,3))*RS(1,3)+(9.E0/19600.E0*RS(2,1)+43.E0
     #/29400.E0*RS(2,4)+31.E0/58800.E0*RS(2,3))*RS(1,4))*RS(1,1)**2
      s1 = s2+((RS(2,2)/20-9.E0/19600.E0*RS(2,3)-RS(2,1)/30+437.E0/17640
     #0.E0*RS(2,4))*RS(1,2)**2+((-43.E0/29400.E0*RS(2,1)+43.E0/29400.E0*
     #RS(2,3))*RS(1,3)+(-43.E0/29400.E0*RS(2,1)-53.E0/88200.E0*RS(2,4)+R
     #S(2,3)/14700)*RS(1,4))*RS(1,2)+(437.E0/176400.E0*RS(2,3)+53.E0/882
     #00.E0*RS(2,1)+9.E0/19600.E0*RS(2,4))*RS(1,3)**2+(RS(2,3)/14700-RS(
     #2,1)/14700)*RS(1,4)*RS(1,3)+(437.E0/176400.E0*RS(2,4)-43.E0/29400.
     #E0*RS(2,1)+31.E0/58800.E0*RS(2,3))*RS(1,4)**2)*RS(1,1)+(-RS(2,1)/2
     #0+RS(2,3)/20)*RS(1,2)**3+((RS(2,3)/30+9.E0/19600.E0*RS(2,1)-437.E0
     #/176400.E0*RS(2,4)-RS(2,2)/20)*RS(1,3)+(437.E0/176400.E0*RS(2,3)-4
     #37.E0/176400.E0*RS(2,1))*RS(1,4))*RS(1,2)**2
      t0 = s1+((-RS(2,2)/30+RS(2,3)/60-43.E0/29400.E0*RS(2,4)-43.E0/2940
     #0.E0*RS(2,1))*RS(1,3)**2+(43.E0/29400.E0*RS(2,3)-RS(2,1)/14700+53.
     #E0/88200.E0*RS(2,4))*RS(1,4)*RS(1,3)+(-53.E0/88200.E0*RS(2,3)+53.E
     #0/88200.E0*RS(2,1))*RS(1,4)**2)*RS(1,2)+(9.E0/19600.E0*RS(2,4)-RS(
     #2,2)/60-437.E0/176400.E0*RS(2,1))*RS(1,3)**3+(-43.E0/29400.E0*RS(2
     #,4)-9.E0/19600.E0*RS(2,3)-31.E0/58800.E0*RS(2,1))*RS(1,4)*RS(1,3)*
     #*2+(43.E0/29400.E0*RS(2,3)-31.E0/58800.E0*RS(2,1)-437.E0/176400.E0
     #*RS(2,4))*RS(1,4)**2*RS(1,3)+(437.E0/176400.E0*RS(2,3)-437.E0/1764
     #00.E0*RS(2,1))*RS(1,4)**3
	ARR(2) = T0

      s2 = (-437.E0/176400.E0*RS(2,2)+437.E0/176400.E0*RS(2,4))*RS(1,1)*
     #*3+((437.E0/176400.E0*RS(2,1)-43.E0/29400.E0*RS(2,2)+31.E0/58800.E
     #0*RS(2,4))*RS(1,2)+(-53.E0/88200.E0*RS(2,4)+53.E0/88200.E0*RS(2,2)
     #)*RS(1,3)+(-31.E0/58800.E0*RS(2,2)+43.E0/29400.E0*RS(2,4)-437.E0/1
     #76400.E0*RS(2,1))*RS(1,4))*RS(1,1)**2
      s1 = s2+((31.E0/58800.E0*RS(2,4)+9.E0/19600.E0*RS(2,2)+43.E0/29400
     #.E0*RS(2,1))*RS(1,2)**2+((-53.E0/88200.E0*RS(2,1)-43.E0/29400.E0*R
     #S(2,2)+RS(2,4)/14700)*RS(1,3)+(RS(2,4)/14700-RS(2,2)/14700)*RS(1,4
     #))*RS(1,2)+(-437.E0/176400.E0*RS(2,2)+437.E0/176400.E0*RS(2,4))*RS
     #(1,3)**2+(43.E0/29400.E0*RS(2,4)-RS(2,2)/14700+53.E0/88200.E0*RS(2
     #,1))*RS(1,4)*RS(1,3)+(-9.E0/19600.E0*RS(2,4)-31.E0/58800.E0*RS(2,2
     #)-43.E0/29400.E0*RS(2,1))*RS(1,4)**2)*RS(1,1)+(RS(2,3)/60-9.E0/196
     #00.E0*RS(2,1)+437.E0/176400.E0*RS(2,4))*RS(1,2)**3+((RS(2,3)/30+43
     #.E0/29400.E0*RS(2,4)+43.E0/29400.E0*RS(2,1)-RS(2,2)/60)*RS(1,3)+(-
     #9.E0/19600.E0*RS(2,1)-53.E0/88200.E0*RS(2,4)-437.E0/176400.E0*RS(2
     #,2))*RS(1,4))*RS(1,2)**2
      t0 = s1+((RS(2,3)/20-RS(2,2)/30+437.E0/176400.E0*RS(2,1)-9.E0/1960
     #0.E0*RS(2,4))*RS(1,3)**2+(43.E0/29400.E0*RS(2,4)-43.E0/29400.E0*RS
     #(2,2))*RS(1,4)*RS(1,3)+(53.E0/88200.E0*RS(2,2)+437.E0/176400.E0*RS
     #(2,4)+9.E0/19600.E0*RS(2,1))*RS(1,4)**2)*RS(1,2)+(-RS(2,2)/20+RS(2
     #,4)/20)*RS(1,3)**3+(-RS(2,3)/20-437.E0/176400.E0*RS(2,1)+RS(2,4)/3
     #0+9.E0/19600.E0*RS(2,2))*RS(1,4)*RS(1,3)**2+(-RS(2,3)/30-43.E0/294
     #00.E0*RS(2,1)-43.E0/29400.E0*RS(2,2)+RS(2,4)/60)*RS(1,4)**2*RS(1,3
     #)+(-437.E0/176400.E0*RS(2,2)+9.E0/19600.E0*RS(2,1)-RS(2,3)/60)*RS(
     #1,4)**3
	ARR(3) = T0

      s2 = (-437.E0/176400.E0*RS(2,3)+9.E0/19600.E0*RS(2,2)-RS(2,4)/60)*
     #RS(1,1)**3+((-43.E0/29400.E0*RS(2,2)-9.E0/19600.E0*RS(2,1)-31.E0/5
     #8800.E0*RS(2,3))*RS(1,2)+(437.E0/176400.E0*RS(2,1)+9.E0/19600.E0*R
     #S(2,2)+53.E0/88200.E0*RS(2,3))*RS(1,3)+(-43.E0/29400.E0*RS(2,2)+RS
     #(2,1)/60-RS(2,4)/30-43.E0/29400.E0*RS(2,3))*RS(1,4))*RS(1,1)**2
      s1 = s2+((-31.E0/58800.E0*RS(2,3)+43.E0/29400.E0*RS(2,1)-437.E0/17
     #6400.E0*RS(2,2))*RS(1,2)**2+((RS(2,1)/14700-RS(2,3)/14700)*RS(1,3)
     #+(43.E0/29400.E0*RS(2,1)+53.E0/88200.E0*RS(2,2)-RS(2,3)/14700)*RS(
     #1,4))*RS(1,2)+(-9.E0/19600.E0*RS(2,2)-53.E0/88200.E0*RS(2,1)-437.E
     #0/176400.E0*RS(2,3))*RS(1,3)**2+(43.E0/29400.E0*RS(2,1)-43.E0/2940
     #0.E0*RS(2,3))*RS(1,4)*RS(1,3)+(-437.E0/176400.E0*RS(2,2)+9.E0/1960
     #0.E0*RS(2,3)+RS(2,1)/30-RS(2,4)/20)*RS(1,4)**2)*RS(1,1)+(-437.E0/1
     #76400.E0*RS(2,3)+437.E0/176400.E0*RS(2,1))*RS(1,2)**3+((-43.E0/294
     #00.E0*RS(2,3)+437.E0/176400.E0*RS(2,2)+31.E0/58800.E0*RS(2,1))*RS(
     #1,3)+(-53.E0/88200.E0*RS(2,1)+53.E0/88200.E0*RS(2,3))*RS(1,4))*RS(
     #1,2)**2
      t0 = s1+((9.E0/19600.E0*RS(2,3)+43.E0/29400.E0*RS(2,2)+31.E0/58800
     #.E0*RS(2,1))*RS(1,3)**2+(RS(2,1)/14700-43.E0/29400.E0*RS(2,3)-53.E
     #0/88200.E0*RS(2,2))*RS(1,4)*RS(1,3)+(-437.E0/176400.E0*RS(2,3)+437
     #.E0/176400.E0*RS(2,1))*RS(1,4)**2)*RS(1,2)+(RS(2,4)/60-9.E0/19600.
     #E0*RS(2,2)+437.E0/176400.E0*RS(2,1))*RS(1,3)**3+(-RS(2,3)/60+43.E0
     #/29400.E0*RS(2,1)+43.E0/29400.E0*RS(2,2)+RS(2,4)/30)*RS(1,4)*RS(1,
     #3)**2+(-RS(2,3)/30-9.E0/19600.E0*RS(2,1)+RS(2,4)/20+437.E0/176400.
     #E0*RS(2,2))*RS(1,4)**2*RS(1,3)+(-RS(2,3)/20+RS(2,1)/20)*RS(1,4)**3
	ARR(4) = T0

      s2 = (-206.E0/11025.E0*RS(2,4)+206.E0/11025.E0*RS(2,2))*RS(1,1)**3
     #+((24.E0/1225.E0*RS(2,2)-206.E0/11025.E0*RS(2,1)+17.E0/11025.E0*RS
     #(2,3)-3.E0/1225.E0*RS(2,4))*RS(1,2)+(13.E0/3675.E0*RS(2,4)-13.E0/3
     #675.E0*RS(2,2))*RS(1,3)+(206.E0/11025.E0*RS(2,1)-17.E0/11025.E0*RS
     #(2,3)-24.E0/1225.E0*RS(2,4)+3.E0/1225.E0*RS(2,2))*RS(1,4))*RS(1,1)
     #**2
      s1 = s2+((-17.E0/11025.E0*RS(2,4)+206.E0/11025.E0*RS(2,2)+3.E0/122
     #5.E0*RS(2,3)-24.E0/1225.E0*RS(2,1))*RS(1,2)**2+((-22.E0/11025.E0*R
     #S(2,3)+22.E0/11025.E0*RS(2,1))*RS(1,3)+(22.E0/11025.E0*RS(2,4)-22.
     #E0/11025.E0*RS(2,2))*RS(1,4))*RS(1,2)+(-13.E0/3675.E0*RS(2,4)+13.E
     #0/3675.E0*RS(2,2))*RS(1,3)**2+(22.E0/11025.E0*RS(2,3)-22.E0/11025.
     #E0*RS(2,1))*RS(1,4)*RS(1,3)+(-206.E0/11025.E0*RS(2,4)+17.E0/11025.
     #E0*RS(2,2)-3.E0/1225.E0*RS(2,3)+24.E0/1225.E0*RS(2,1))*RS(1,4)**2)
     #*RS(1,1)+(-206.E0/11025.E0*RS(2,1)+206.E0/11025.E0*RS(2,3))*RS(1,2
     #)**3+((24.E0/1225.E0*RS(2,3)-206.E0/11025.E0*RS(2,2)+17.E0/11025.E
     #0*RS(2,4)-3.E0/1225.E0*RS(2,1))*RS(1,3)+(13.E0/3675.E0*RS(2,1)-13.
     #E0/3675.E0*RS(2,3))*RS(1,4))*RS(1,2)**2
      t0 = s1+((206.E0/11025.E0*RS(2,3)-17.E0/11025.E0*RS(2,1)-24.E0/122
     #5.E0*RS(2,2)+3.E0/1225.E0*RS(2,4))*RS(1,3)**2+(-22.E0/11025.E0*RS(
     #2,4)+22.E0/11025.E0*RS(2,2))*RS(1,4)*RS(1,3)+(-13.E0/3675.E0*RS(2,
     #1)+13.E0/3675.E0*RS(2,3))*RS(1,4)**2)*RS(1,2)+(-206.E0/11025.E0*RS
     #(2,2)+206.E0/11025.E0*RS(2,4))*RS(1,3)**3+(24.E0/1225.E0*RS(2,4)-3
     #.E0/1225.E0*RS(2,2)-206.E0/11025.E0*RS(2,3)+17.E0/11025.E0*RS(2,1)
     #)*RS(1,4)*RS(1,3)**2+(-24.E0/1225.E0*RS(2,3)-17.E0/11025.E0*RS(2,2
     #)+206.E0/11025.E0*RS(2,4)+3.E0/1225.E0*RS(2,1))*RS(1,4)**2*RS(1,3)
     #+(-206.E0/11025.E0*RS(2,3)+206.E0/11025.E0*RS(2,1))*RS(1,4)**3
	ARR(5) = T0

C	----------------------
C	COMPUTE FOR FACTOR ASS
C	----------------------
      s2 = (-RS(1,2)/20+RS(1,4)/20)*RS(2,1)**3+((RS(1,1)/20+437.E0/17640
     #0.E0*RS(1,3)-RS(1,2)/30-9.E0/19600.E0*RS(1,4))*RS(2,2)+(-437.E0/17
     #6400.E0*RS(1,2)+437.E0/176400.E0*RS(1,4))*RS(2,3)+(-RS(1,1)/20+9.E
     #0/19600.E0*RS(1,2)-437.E0/176400.E0*RS(1,3)+RS(1,4)/30)*RS(2,4))*R
     #S(2,1)**2
      s1 = s2+((RS(1,1)/30+43.E0/29400.E0*RS(1,4)+43.E0/29400.E0*RS(1,3)
     #-RS(1,2)/60)*RS(2,2)**2+((RS(1,4)/14700-53.E0/88200.E0*RS(1,3)-43.
     #E0/29400.E0*RS(1,2))*RS(2,3)+(43.E0/29400.E0*RS(1,4)-43.E0/29400.E
     #0*RS(1,2))*RS(2,4))*RS(2,2)+(53.E0/88200.E0*RS(1,2)-53.E0/88200.E0
     #*RS(1,4))*RS(2,3)**2+(-RS(1,2)/14700+53.E0/88200.E0*RS(1,3)+43.E0/
     #29400.E0*RS(1,4))*RS(2,4)*RS(2,3)+(-RS(1,1)/30+RS(1,4)/60-43.E0/29
     #400.E0*RS(1,2)-43.E0/29400.E0*RS(1,3))*RS(2,4)**2)*RS(2,1)+(437.E0
     #/176400.E0*RS(1,4)-9.E0/19600.E0*RS(1,3)+RS(1,1)/60)*RS(2,2)**3+((
     #9.E0/19600.E0*RS(1,2)+43.E0/29400.E0*RS(1,3)+31.E0/58800.E0*RS(1,4
     #))*RS(2,3)+(-9.E0/19600.E0*RS(1,3)-53.E0/88200.E0*RS(1,4)-437.E0/1
     #76400.E0*RS(1,2))*RS(2,4))*RS(2,2)**2
      t0 = s1+((-43.E0/29400.E0*RS(1,2)+31.E0/58800.E0*RS(1,4)+437.E0/17
     #6400.E0*RS(1,3))*RS(2,3)**2+(RS(1,4)/14700-RS(1,2)/14700)*RS(2,4)*
     #RS(2,3)+(9.E0/19600.E0*RS(1,3)+437.E0/176400.E0*RS(1,4)+53.E0/8820
     #0.E0*RS(1,2))*RS(2,4)**2)*RS(2,2)+(-437.E0/176400.E0*RS(1,2)+437.E
     #0/176400.E0*RS(1,4))*RS(2,3)**3+(-437.E0/176400.E0*RS(1,3)-31.E0/5
     #8800.E0*RS(1,2)+43.E0/29400.E0*RS(1,4))*RS(2,4)*RS(2,3)**2+(-43.E0
     #/29400.E0*RS(1,3)-9.E0/19600.E0*RS(1,4)-31.E0/58800.E0*RS(1,2))*RS
     #(2,4)**2*RS(2,3)+(9.E0/19600.E0*RS(1,3)-437.E0/176400.E0*RS(1,2)-R
     #S(1,1)/60)*RS(2,4)**3
	ASS(1) = T0

      s2 = (-RS(1,2)/60+9.E0/19600.E0*RS(1,4)-437.E0/176400.E0*RS(1,3))*
     #RS(2,1)**3+((-RS(1,2)/30-43.E0/29400.E0*RS(1,4)+RS(1,1)/60-43.E0/2
     #9400.E0*RS(1,3))*RS(2,2)+(9.E0/19600.E0*RS(1,4)+437.E0/176400.E0*R
     #S(1,1)+53.E0/88200.E0*RS(1,3))*RS(2,3)+(-9.E0/19600.E0*RS(1,1)-31.
     #E0/58800.E0*RS(1,3)-43.E0/29400.E0*RS(1,4))*RS(2,4))*RS(2,1)**2
      s1 = s2+((RS(1,1)/30+9.E0/19600.E0*RS(1,3)-RS(1,2)/20-437.E0/17640
     #0.E0*RS(1,4))*RS(2,2)**2+((43.E0/29400.E0*RS(1,1)-43.E0/29400.E0*R
     #S(1,3))*RS(2,3)+(53.E0/88200.E0*RS(1,4)-RS(1,3)/14700+43.E0/29400.
     #E0*RS(1,1))*RS(2,4))*RS(2,2)+(-53.E0/88200.E0*RS(1,1)-9.E0/19600.E
     #0*RS(1,4)-437.E0/176400.E0*RS(1,3))*RS(2,3)**2+(-RS(1,3)/14700+RS(
     #1,1)/14700)*RS(2,4)*RS(2,3)+(43.E0/29400.E0*RS(1,1)-31.E0/58800.E0
     #*RS(1,3)-437.E0/176400.E0*RS(1,4))*RS(2,4)**2)*RS(2,1)+(RS(1,1)/20
     #-RS(1,3)/20)*RS(2,2)**3+((-9.E0/19600.E0*RS(1,1)-RS(1,3)/30+RS(1,2
     #)/20+437.E0/176400.E0*RS(1,4))*RS(2,3)+(-437.E0/176400.E0*RS(1,3)+
     #437.E0/176400.E0*RS(1,1))*RS(2,4))*RS(2,2)**2
      t0 = s1+((-RS(1,3)/60+RS(1,2)/30+43.E0/29400.E0*RS(1,4)+43.E0/2940
     #0.E0*RS(1,1))*RS(2,3)**2+(RS(1,1)/14700-43.E0/29400.E0*RS(1,3)-53.
     #E0/88200.E0*RS(1,4))*RS(2,4)*RS(2,3)+(-53.E0/88200.E0*RS(1,1)+53.E
     #0/88200.E0*RS(1,3))*RS(2,4)**2)*RS(2,2)+(RS(1,2)/60+437.E0/176400.
     #E0*RS(1,1)-9.E0/19600.E0*RS(1,4))*RS(2,3)**3+(31.E0/58800.E0*RS(1,
     #1)+9.E0/19600.E0*RS(1,3)+43.E0/29400.E0*RS(1,4))*RS(2,4)*RS(2,3)**
     #2+(-43.E0/29400.E0*RS(1,3)+31.E0/58800.E0*RS(1,1)+437.E0/176400.E0
     #*RS(1,4))*RS(2,4)**2*RS(2,3)+(-437.E0/176400.E0*RS(1,3)+437.E0/176
     #400.E0*RS(1,1))*RS(2,4)**3
	ASS(2) = T0

      s2 = (-437.E0/176400.E0*RS(1,4)+437.E0/176400.E0*RS(1,2))*RS(2,1)*
     #*3+((-437.E0/176400.E0*RS(1,1)+43.E0/29400.E0*RS(1,2)-31.E0/58800.
     #E0*RS(1,4))*RS(2,2)+(53.E0/88200.E0*RS(1,4)-53.E0/88200.E0*RS(1,2)
     #)*RS(2,3)+(31.E0/58800.E0*RS(1,2)+437.E0/176400.E0*RS(1,1)-43.E0/2
     #9400.E0*RS(1,4))*RS(2,4))*RS(2,1)**2
      s1 = s2+((-9.E0/19600.E0*RS(1,2)-43.E0/29400.E0*RS(1,1)-31.E0/5880
     #0.E0*RS(1,4))*RS(2,2)**2+((-RS(1,4)/14700+53.E0/88200.E0*RS(1,1)+4
     #3.E0/29400.E0*RS(1,2))*RS(2,3)+(-RS(1,4)/14700+RS(1,2)/14700)*RS(2
     #,4))*RS(2,2)+(-437.E0/176400.E0*RS(1,4)+437.E0/176400.E0*RS(1,2))*
     #RS(2,3)**2+(-43.E0/29400.E0*RS(1,4)-53.E0/88200.E0*RS(1,1)+RS(1,2)
     #/14700)*RS(2,4)*RS(2,3)+(9.E0/19600.E0*RS(1,4)+43.E0/29400.E0*RS(1
     #,1)+31.E0/58800.E0*RS(1,2))*RS(2,4)**2)*RS(2,1)+(-RS(1,3)/60-437.E
     #0/176400.E0*RS(1,4)+9.E0/19600.E0*RS(1,1))*RS(2,2)**3+((-43.E0/294
     #00.E0*RS(1,4)-RS(1,3)/30+RS(1,2)/60-43.E0/29400.E0*RS(1,1))*RS(2,3
     #)+(9.E0/19600.E0*RS(1,1)+437.E0/176400.E0*RS(1,2)+53.E0/88200.E0*R
     #S(1,4))*RS(2,4))*RS(2,2)**2
      t0 = s1+((-RS(1,3)/20+9.E0/19600.E0*RS(1,4)+RS(1,2)/30-437.E0/1764
     #00.E0*RS(1,1))*RS(2,3)**2+(43.E0/29400.E0*RS(1,2)-43.E0/29400.E0*R
     #S(1,4))*RS(2,4)*RS(2,3)+(-9.E0/19600.E0*RS(1,1)-437.E0/176400.E0*R
     #S(1,4)-53.E0/88200.E0*RS(1,2))*RS(2,4)**2)*RS(2,2)+(RS(1,2)/20-RS(
     #1,4)/20)*RS(2,3)**3+(437.E0/176400.E0*RS(1,1)+RS(1,3)/20-9.E0/1960
     #0.E0*RS(1,2)-RS(1,4)/30)*RS(2,4)*RS(2,3)**2+(43.E0/29400.E0*RS(1,1
     #)-RS(1,4)/60+RS(1,3)/30+43.E0/29400.E0*RS(1,2))*RS(2,4)**2*RS(2,3)
     #+(437.E0/176400.E0*RS(1,2)-9.E0/19600.E0*RS(1,1)+RS(1,3)/60)*RS(2,
     #4)**3
	ASS(3) = T0

      s2 = (437.E0/176400.E0*RS(1,3)+RS(1,4)/60-9.E0/19600.E0*RS(1,2))*R
     #S(2,1)**3+((43.E0/29400.E0*RS(1,2)+9.E0/19600.E0*RS(1,1)+31.E0/588
     #00.E0*RS(1,3))*RS(2,2)+(-437.E0/176400.E0*RS(1,1)-53.E0/88200.E0*R
     #S(1,3)-9.E0/19600.E0*RS(1,2))*RS(2,3)+(-RS(1,1)/60+43.E0/29400.E0*
     #RS(1,3)+RS(1,4)/30+43.E0/29400.E0*RS(1,2))*RS(2,4))*RS(2,1)**2
      s1 = s2+((437.E0/176400.E0*RS(1,2)-43.E0/29400.E0*RS(1,1)+31.E0/58
     #800.E0*RS(1,3))*RS(2,2)**2+((-RS(1,1)/14700+RS(1,3)/14700)*RS(2,3)
     #+(RS(1,3)/14700-43.E0/29400.E0*RS(1,1)-53.E0/88200.E0*RS(1,2))*RS(
     #2,4))*RS(2,2)+(437.E0/176400.E0*RS(1,3)+53.E0/88200.E0*RS(1,1)+9.E
     #0/19600.E0*RS(1,2))*RS(2,3)**2+(-43.E0/29400.E0*RS(1,1)+43.E0/2940
     #0.E0*RS(1,3))*RS(2,4)*RS(2,3)+(-RS(1,1)/30-9.E0/19600.E0*RS(1,3)+R
     #S(1,4)/20+437.E0/176400.E0*RS(1,2))*RS(2,4)**2)*RS(2,1)+(-437.E0/1
     #76400.E0*RS(1,1)+437.E0/176400.E0*RS(1,3))*RS(2,2)**3+((-437.E0/17
     #6400.E0*RS(1,2)+43.E0/29400.E0*RS(1,3)-31.E0/58800.E0*RS(1,1))*RS(
     #2,3)+(-53.E0/88200.E0*RS(1,3)+53.E0/88200.E0*RS(1,1))*RS(2,4))*RS(
     #2,2)**2
      t0 = s1+((-43.E0/29400.E0*RS(1,2)-9.E0/19600.E0*RS(1,3)-31.E0/5880
     #0.E0*RS(1,1))*RS(2,3)**2+(53.E0/88200.E0*RS(1,2)-RS(1,1)/14700+43.
     #E0/29400.E0*RS(1,3))*RS(2,4)*RS(2,3)+(-437.E0/176400.E0*RS(1,1)+43
     #7.E0/176400.E0*RS(1,3))*RS(2,4)**2)*RS(2,2)+(-437.E0/176400.E0*RS(
     #1,1)-RS(1,4)/60+9.E0/19600.E0*RS(1,2))*RS(2,3)**3+(-RS(1,4)/30+RS(
     #1,3)/60-43.E0/29400.E0*RS(1,1)-43.E0/29400.E0*RS(1,2))*RS(2,4)*RS(
     #2,3)**2+(9.E0/19600.E0*RS(1,1)+RS(1,3)/30-437.E0/176400.E0*RS(1,2)
     #-RS(1,4)/20)*RS(2,4)**2*RS(2,3)+(RS(1,3)/20-RS(1,1)/20)*RS(2,4)**3
	ASS(4) = T0

      s2 = (-206.E0/11025.E0*RS(1,2)+206.E0/11025.E0*RS(1,4))*RS(2,1)**3
     #+((-17.E0/11025.E0*RS(1,3)+206.E0/11025.E0*RS(1,1)+3.E0/1225.E0*RS
     #(1,4)-24.E0/1225.E0*RS(1,2))*RS(2,2)+(13.E0/3675.E0*RS(1,2)-13.E0/
     #3675.E0*RS(1,4))*RS(2,3)+(-3.E0/1225.E0*RS(1,2)+17.E0/11025.E0*RS(
     #1,3)+24.E0/1225.E0*RS(1,4)-206.E0/11025.E0*RS(1,1))*RS(2,4))*RS(2,
     #1)**2
      s1 = s2+((-3.E0/1225.E0*RS(1,3)+24.E0/1225.E0*RS(1,1)+17.E0/11025.
     #E0*RS(1,4)-206.E0/11025.E0*RS(1,2))*RS(2,2)**2+((22.E0/11025.E0*RS
     #(1,3)-22.E0/11025.E0*RS(1,1))*RS(2,3)+(-22.E0/11025.E0*RS(1,4)+22.
     #E0/11025.E0*RS(1,2))*RS(2,4))*RS(2,2)+(-13.E0/3675.E0*RS(1,2)+13.E
     #0/3675.E0*RS(1,4))*RS(2,3)**2+(-22.E0/11025.E0*RS(1,3)+22.E0/11025
     #.E0*RS(1,1))*RS(2,4)*RS(2,3)+(-24.E0/1225.E0*RS(1,1)-17.E0/11025.E
     #0*RS(1,2)+206.E0/11025.E0*RS(1,4)+3.E0/1225.E0*RS(1,3))*RS(2,4)**2
     #)*RS(2,1)+(206.E0/11025.E0*RS(1,1)-206.E0/11025.E0*RS(1,3))*RS(2,2
     #)**3+((3.E0/1225.E0*RS(1,1)+206.E0/11025.E0*RS(1,2)-24.E0/1225.E0*
     #RS(1,3)-17.E0/11025.E0*RS(1,4))*RS(2,3)+(-13.E0/3675.E0*RS(1,1)+13
     #.E0/3675.E0*RS(1,3))*RS(2,4))*RS(2,2)**2
      t0 = s1+((24.E0/1225.E0*RS(1,2)+17.E0/11025.E0*RS(1,1)-206.E0/1102
     #5.E0*RS(1,3)-3.E0/1225.E0*RS(1,4))*RS(2,3)**2+(-22.E0/11025.E0*RS(
     #1,2)+22.E0/11025.E0*RS(1,4))*RS(2,4)*RS(2,3)+(13.E0/3675.E0*RS(1,1
     #)-13.E0/3675.E0*RS(1,3))*RS(2,4)**2)*RS(2,2)+(206.E0/11025.E0*RS(1
     #,2)-206.E0/11025.E0*RS(1,4))*RS(2,3)**3+(206.E0/11025.E0*RS(1,3)-1
     #7.E0/11025.E0*RS(1,1)-24.E0/1225.E0*RS(1,4)+3.E0/1225.E0*RS(1,2))*
     #RS(2,4)*RS(2,3)**2+(24.E0/1225.E0*RS(1,3)+17.E0/11025.E0*RS(1,2)-3
     #.E0/1225.E0*RS(1,1)-206.E0/11025.E0*RS(1,4))*RS(2,4)**2*RS(2,3)+(2
     #06.E0/11025.E0*RS(1,3)-206.E0/11025.E0*RS(1,1))*RS(2,4)**3
	ASS(5) = T0

C	----------------------
C	COMPUTE FOR FACTOR ARS
C	----------------------
      s2 = (-RS(2,1)*RS(2,4)/20-RS(2,4)**2/60+RS(2,2)**2/60+RS(2,1)*RS(2
     #,2)/20)*RS(1,1)**2
      s4 = ((-437.E0/176400.E0*RS(2,1)*RS(2,3)-43.E0/58800.E0*RS(2,4)**2
     #-43.E0/58800.E0*RS(2,2)*RS(2,3)-RS(2,3)*RS(2,4)/29400+RS(2,2)**2/6
     #0-43.E0/58800.E0*RS(2,4)*RS(2,2)+53.E0/176400.E0*RS(2,3)**2+9.E0/1
     #9600.E0*RS(2,1)*RS(2,4)-RS(2,1)**2/20)*RS(1,2)+(53.E0/176400.E0*RS
     #(2,3)*RS(2,4)+437.E0/176400.E0*RS(2,1)*RS(2,2)-437.E0/176400.E0*RS
     #(2,1)*RS(2,4)+43.E0/58800.E0*RS(2,2)**2-53.E0/176400.E0*RS(2,2)*RS
     #(2,3)-43.E0/58800.E0*RS(2,4)**2)*RS(1,3)+(43.E0/58800.E0*RS(2,4)*R
     #S(2,2)+RS(2,1)**2/20+RS(2,2)*RS(2,3)/29400-9.E0/19600.E0*RS(2,1)*R
     #S(2,2)+43.E0/58800.E0*RS(2,3)*RS(2,4)-RS(2,4)**2/60+43.E0/58800.E0
     #*RS(2,2)**2-53.E0/176400.E0*RS(2,3)**2+437.E0/176400.E0*RS(2,1)*RS
     #(2,3))*RS(1,4))*RS(1,1)
      s5 = (-43.E0/58800.E0*RS(2,3)**2-437.E0/176400.E0*RS(2,4)*RS(2,2)-
     #RS(2,3)*RS(2,4)/29400-43.E0/58800.E0*RS(2,1)*RS(2,4)-RS(2,1)*RS(2,
     #2)/60+9.E0/19600.E0*RS(2,2)*RS(2,3)-RS(2,1)**2/60+53.E0/176400.E0*
     #RS(2,4)**2-43.E0/58800.E0*RS(2,1)*RS(2,3))*RS(1,2)**2
      s3 = s4+s5
      s1 = s2+s3
      s2 = s1+((-31.E0/58800.E0*RS(2,3)*RS(2,4)-29.E0/58800.E0*RS(2,4)*R
     #S(2,2)-RS(2,4)**2/29400+43.E0/58800.E0*RS(2,1)*RS(2,2)+53.E0/17640
     #0.E0*RS(2,1)*RS(2,3)-RS(2,1)*RS(2,4)/29400-9.E0/19600.E0*RS(2,2)**
     #2-437.E0/176400.E0*RS(2,3)**2)*RS(1,3)+(437.E0/176400.E0*RS(2,2)**
     #2+43.E0/58800.E0*RS(2,1)*RS(2,2)-43.E0/58800.E0*RS(2,1)*RS(2,4)+29
     #.E0/58800.E0*RS(2,2)*RS(2,3)-29.E0/58800.E0*RS(2,3)*RS(2,4)-437.E0
     #/176400.E0*RS(2,4)**2)*RS(1,4))*RS(1,2)
      t0 = s2+(-437.E0/176400.E0*RS(2,3)*RS(2,4)-43.E0/58800.E0*RS(2,4)*
     #*2+437.E0/176400.E0*RS(2,2)*RS(2,3)+43.E0/58800.E0*RS(2,2)**2+53.E
     #0/176400.E0*RS(2,1)*RS(2,4)-53.E0/176400.E0*RS(2,1)*RS(2,2))*RS(1,
     #3)**2+(-43.E0/58800.E0*RS(2,1)*RS(2,4)+9.E0/19600.E0*RS(2,4)**2+29
     #.E0/58800.E0*RS(2,4)*RS(2,2)-53.E0/176400.E0*RS(2,1)*RS(2,3)+RS(2,
     #1)*RS(2,2)/29400+31.E0/58800.E0*RS(2,2)*RS(2,3)+RS(2,2)**2/29400+4
     #37.E0/176400.E0*RS(2,3)**2)*RS(1,4)*RS(1,3)+(43.E0/58800.E0*RS(2,1
     #)*RS(2,2)+43.E0/58800.E0*RS(2,3)**2+437.E0/176400.E0*RS(2,4)*RS(2,
     #2)+43.E0/58800.E0*RS(2,1)*RS(2,3)+RS(2,1)**2/60+RS(2,2)*RS(2,3)/29
     #400+RS(2,1)*RS(2,4)/60-53.E0/176400.E0*RS(2,2)**2-9.E0/19600.E0*RS
     #(2,3)*RS(2,4))*RS(1,4)**2
	ARS(1) = T0

      s2 = (43.E0/58800.E0*RS(2,4)**2+RS(2,1)*RS(2,2)/60+RS(2,3)*RS(2,4)
     #/29400+43.E0/58800.E0*RS(2,3)*RS(2,2)-53.E0/176400.E0*RS(2,3)**2+R
     #S(2,2)**2/60+437.E0/176400.E0*RS(2,1)*RS(2,3)-9.E0/19600.E0*RS(2,1
     #)*RS(2,4)+43.E0/58800.E0*RS(2,4)*RS(2,2))*RS(1,1)**2
      s4 = ((43.E0/58800.E0*RS(2,1)*RS(2,4)-53.E0/176400.E0*RS(2,4)**2-9
     #.E0/19600.E0*RS(2,3)*RS(2,2)+RS(2,2)**2/20+RS(2,3)*RS(2,4)/29400-R
     #S(2,1)**2/60+43.E0/58800.E0*RS(2,3)**2+43.E0/58800.E0*RS(2,1)*RS(2
     #,3)+437.E0/176400.E0*RS(2,4)*RS(2,2))*RS(1,2)+(29.E0/58800.E0*RS(2
     #,3)*RS(2,4)+43.E0/58800.E0*RS(2,3)*RS(2,2)-437.E0/176400.E0*RS(2,1
     #)**2-29.E0/58800.E0*RS(2,1)*RS(2,4)+437.E0/176400.E0*RS(2,3)**2-43
     #.E0/58800.E0*RS(2,1)*RS(2,2))*RS(1,3)+(-53.E0/176400.E0*RS(2,4)*RS
     #(2,2)+9.E0/19600.E0*RS(2,1)**2+31.E0/58800.E0*RS(2,3)*RS(2,4)+RS(2
     #,3)*RS(2,2)/29400-43.E0/58800.E0*RS(2,1)*RS(2,2)+437.E0/176400.E0*
     #RS(2,4)**2+RS(2,3)**2/29400+29.E0/58800.E0*RS(2,1)*RS(2,3))*RS(1,4
     #))*RS(1,1)
      s5 = (RS(2,3)*RS(2,2)/20-RS(2,1)*RS(2,2)/20-RS(2,1)**2/60+RS(2,3)*
     #*2/60)*RS(1,2)**2
      s3 = s4+s5
      s1 = s2+s3
      s2 = s1+((RS(2,3)**2/60-RS(2,1)*RS(2,4)/29400-437.E0/176400.E0*RS(
     #2,4)*RS(2,2)+9.E0/19600.E0*RS(2,1)*RS(2,2)-RS(2,2)**2/20-43.E0/588
     #00.E0*RS(2,1)*RS(2,3)-43.E0/58800.E0*RS(2,3)*RS(2,4)+53.E0/176400.
     #E0*RS(2,4)**2-43.E0/58800.E0*RS(2,1)**2)*RS(1,3)+(437.E0/176400.E0
     #*RS(2,3)*RS(2,2)+53.E0/176400.E0*RS(2,1)*RS(2,4)-43.E0/58800.E0*RS
     #(2,1)**2-53.E0/176400.E0*RS(2,3)*RS(2,4)-437.E0/176400.E0*RS(2,1)*
     #RS(2,2)+43.E0/58800.E0*RS(2,3)**2)*RS(1,4))*RS(1,2)
      t0 = s2+(53.E0/176400.E0*RS(2,1)**2-RS(2,1)*RS(2,4)/29400-43.E0/58
     #800.E0*RS(2,4)**2-RS(2,3)*RS(2,2)/60-437.E0/176400.E0*RS(2,1)*RS(2
     #,3)-43.E0/58800.E0*RS(2,4)*RS(2,2)+9.E0/19600.E0*RS(2,3)*RS(2,4)-R
     #S(2,2)**2/60-43.E0/58800.E0*RS(2,1)*RS(2,2))*RS(1,3)**2+(-RS(2,1)*
     #*2/29400-31.E0/58800.E0*RS(2,1)*RS(2,4)+53.E0/176400.E0*RS(2,4)*RS
     #(2,2)-29.E0/58800.E0*RS(2,1)*RS(2,3)-RS(2,1)*RS(2,2)/29400-437.E0/
     #176400.E0*RS(2,4)**2-9.E0/19600.E0*RS(2,3)**2+43.E0/58800.E0*RS(2,
     #3)*RS(2,2))*RS(1,4)*RS(1,3)+(-53.E0/176400.E0*RS(2,3)*RS(2,2)+437.
     #E0/176400.E0*RS(2,3)*RS(2,4)-43.E0/58800.E0*RS(2,1)**2+53.E0/17640
     #0.E0*RS(2,1)*RS(2,2)+43.E0/58800.E0*RS(2,3)**2-437.E0/176400.E0*RS
     #(2,1)*RS(2,4))*RS(1,4)**2
	ARS(2) = T0

      s2 = (-43.E0/58800.E0*RS(2,2)**2+43.E0/58800.E0*RS(2,4)**2+53.E0/1
     #76400.E0*RS(2,3)*RS(2,2)+437.E0/176400.E0*RS(2,1)*RS(2,4)-53.E0/17
     #6400.E0*RS(2,3)*RS(2,4)-437.E0/176400.E0*RS(2,1)*RS(2,2))*RS(1,1)*
     #*2
      s4 = ((-53.E0/176400.E0*RS(2,1)*RS(2,3)+437.E0/176400.E0*RS(2,1)**
     #2+29.E0/58800.E0*RS(2,2)*RS(2,4)-43.E0/58800.E0*RS(2,3)*RS(2,2)+RS
     #(2,3)*RS(2,4)/29400+RS(2,4)**2/29400+31.E0/58800.E0*RS(2,1)*RS(2,4
     #)+9.E0/19600.E0*RS(2,2)**2)*RS(1,2)+(43.E0/58800.E0*RS(2,4)**2+53.
     #E0/176400.E0*RS(2,1)*RS(2,2)-437.E0/176400.E0*RS(2,3)*RS(2,2)-53.E
     #0/176400.E0*RS(2,1)*RS(2,4)-43.E0/58800.E0*RS(2,2)**2+437.E0/17640
     #0.E0*RS(2,3)*RS(2,4))*RS(1,3)+(53.E0/176400.E0*RS(2,1)*RS(2,3)-RS(
     #2,3)*RS(2,2)/29400-RS(2,2)**2/29400-29.E0/58800.E0*RS(2,2)*RS(2,4)
     #+43.E0/58800.E0*RS(2,3)*RS(2,4)-9.E0/19600.E0*RS(2,4)**2-31.E0/588
     #00.E0*RS(2,1)*RS(2,2)-437.E0/176400.E0*RS(2,1)**2)*RS(1,4))*RS(1,1
     #)
      s5 = (43.E0/58800.E0*RS(2,1)*RS(2,3)-53.E0/176400.E0*RS(2,4)**2+43
     #7.E0/176400.E0*RS(2,2)*RS(2,4)+43.E0/58800.E0*RS(2,3)*RS(2,4)+43.E
     #0/58800.E0*RS(2,1)**2+RS(2,1)*RS(2,4)/29400+RS(2,3)**2/60-9.E0/196
     #00.E0*RS(2,1)*RS(2,2)+RS(2,3)*RS(2,2)/60)*RS(1,2)**2
      s3 = s4+s5
      s1 = s2+s3
      s2 = s1+((43.E0/58800.E0*RS(2,1)*RS(2,2)+RS(2,1)*RS(2,4)/29400+RS(
     #2,3)**2/20+43.E0/58800.E0*RS(2,2)*RS(2,4)+43.E0/58800.E0*RS(2,4)**
     #2-9.E0/19600.E0*RS(2,3)*RS(2,4)+437.E0/176400.E0*RS(2,1)*RS(2,3)-R
     #S(2,2)**2/60-53.E0/176400.E0*RS(2,1)**2)*RS(1,3)+(29.E0/58800.E0*R
     #S(2,1)*RS(2,4)-437.E0/176400.E0*RS(2,2)**2-43.E0/58800.E0*RS(2,3)*
     #RS(2,2)+437.E0/176400.E0*RS(2,4)**2-29.E0/58800.E0*RS(2,1)*RS(2,2)
     #+43.E0/58800.E0*RS(2,3)*RS(2,4))*RS(1,4))*RS(1,2)
      t0 = s2+(-RS(2,2)**2/60-RS(2,3)*RS(2,2)/20+RS(2,4)**2/60+RS(2,3)*R
     #S(2,4)/20)*RS(1,3)**2+(-RS(2,1)*RS(2,2)/29400-43.E0/58800.E0*RS(2,
     #1)*RS(2,4)+RS(2,4)**2/60+9.E0/19600.E0*RS(2,3)*RS(2,2)-437.E0/1764
     #00.E0*RS(2,1)*RS(2,3)-43.E0/58800.E0*RS(2,2)*RS(2,4)-RS(2,3)**2/20
     #-43.E0/58800.E0*RS(2,2)**2+53.E0/176400.E0*RS(2,1)**2)*RS(1,4)*RS(
     #1,3)+(-43.E0/58800.E0*RS(2,1)*RS(2,3)-43.E0/58800.E0*RS(2,3)*RS(2,
     #2)-43.E0/58800.E0*RS(2,1)**2+53.E0/176400.E0*RS(2,2)**2-RS(2,3)**2
     #/60-RS(2,3)*RS(2,4)/60-437.E0/176400.E0*RS(2,2)*RS(2,4)+9.E0/19600
     #.E0*RS(2,1)*RS(2,4)-RS(2,1)*RS(2,2)/29400)*RS(1,4)**2
	ARS(3) = T0

      s2 = (-RS(2,4)**2/60-437.E0/176400.E0*RS(2,1)*RS(2,3)-43.E0/58800.
     #E0*RS(2,2)**2+53.E0/176400.E0*RS(2,3)**2-RS(2,3)*RS(2,2)/29400-43.
     #E0/58800.E0*RS(2,3)*RS(2,4)-RS(2,1)*RS(2,4)/60+9.E0/19600.E0*RS(2,
     #1)*RS(2,2)-43.E0/58800.E0*RS(2,2)*RS(2,4))*RS(1,1)**2
      s4 = ((-437.E0/176400.E0*RS(2,2)**2-9.E0/19600.E0*RS(2,1)**2-31.E0
     #/58800.E0*RS(2,3)*RS(2,2)+43.E0/58800.E0*RS(2,1)*RS(2,4)-29.E0/588
     #00.E0*RS(2,1)*RS(2,3)-RS(2,3)**2/29400+53.E0/176400.E0*RS(2,2)*RS(
     #2,4)-RS(2,3)*RS(2,4)/29400)*RS(1,2)+(29.E0/58800.E0*RS(2,1)*RS(2,2
     #)+437.E0/176400.E0*RS(2,1)**2-29.E0/58800.E0*RS(2,3)*RS(2,2)+43.E0
     #/58800.E0*RS(2,1)*RS(2,4)-437.E0/176400.E0*RS(2,3)**2-43.E0/58800.
     #E0*RS(2,3)*RS(2,4))*RS(1,3)+(-43.E0/58800.E0*RS(2,3)**2+53.E0/1764
     #00.E0*RS(2,2)**2-RS(2,4)**2/20-43.E0/58800.E0*RS(2,1)*RS(2,3)+RS(2
     #,1)**2/60-43.E0/58800.E0*RS(2,1)*RS(2,2)-437.E0/176400.E0*RS(2,2)*
     #RS(2,4)-RS(2,3)*RS(2,2)/29400+9.E0/19600.E0*RS(2,3)*RS(2,4))*RS(1,
     #4))*RS(1,1)
      s5 = (-43.E0/58800.E0*RS(2,3)**2+437.E0/176400.E0*RS(2,1)*RS(2,2)-
     #437.E0/176400.E0*RS(2,3)*RS(2,2)-53.E0/176400.E0*RS(2,1)*RS(2,4)+5
     #3.E0/176400.E0*RS(2,3)*RS(2,4)+43.E0/58800.E0*RS(2,1)**2)*RS(1,2)*
     #*2
      s3 = s4+s5
      s1 = s2+s3
      s2 = s1+((RS(2,1)**2/29400+9.E0/19600.E0*RS(2,3)**2+29.E0/58800.E0
     #*RS(2,1)*RS(2,3)-53.E0/176400.E0*RS(2,2)*RS(2,4)+RS(2,1)*RS(2,4)/2
     #9400+31.E0/58800.E0*RS(2,1)*RS(2,2)+437.E0/176400.E0*RS(2,2)**2-43
     #.E0/58800.E0*RS(2,3)*RS(2,4))*RS(1,3)+(43.E0/58800.E0*RS(2,1)**2+5
     #3.E0/176400.E0*RS(2,3)*RS(2,2)-53.E0/176400.E0*RS(2,1)*RS(2,2)+437
     #.E0/176400.E0*RS(2,1)*RS(2,4)-437.E0/176400.E0*RS(2,3)*RS(2,4)-43.
     #E0/58800.E0*RS(2,3)**2)*RS(1,4))*RS(1,2)
      t0 = s2+(RS(2,1)*RS(2,2)/29400+43.E0/58800.E0*RS(2,2)*RS(2,4)+43.E
     #0/58800.E0*RS(2,2)**2+43.E0/58800.E0*RS(2,1)*RS(2,4)+RS(2,3)*RS(2,
     #4)/60+RS(2,4)**2/60-9.E0/19600.E0*RS(2,3)*RS(2,2)+437.E0/176400.E0
     #*RS(2,1)*RS(2,3)-53.E0/176400.E0*RS(2,1)**2)*RS(1,3)**2+(437.E0/17
     #6400.E0*RS(2,2)*RS(2,4)+RS(2,1)*RS(2,2)/29400+43.E0/58800.E0*RS(2,
     #3)*RS(2,2)+43.E0/58800.E0*RS(2,1)*RS(2,3)-RS(2,3)**2/60+RS(2,4)**2
     #/20-53.E0/176400.E0*RS(2,2)**2-9.E0/19600.E0*RS(2,1)*RS(2,4)+43.E0
     #/58800.E0*RS(2,1)**2)*RS(1,4)*RS(1,3)+(-RS(2,3)*RS(2,4)/20+RS(2,1)
     #*RS(2,4)/20-RS(2,3)**2/60+RS(2,1)**2/60)*RS(1,4)**2
	ARS(4) = T0

      s2 = (-12.E0/1225.E0*RS(2,4)**2+12.E0/1225.E0*RS(2,2)**2-11.E0/110
     #25.E0*RS(2,3)*RS(2,2)+11.E0/11025.E0*RS(2,3)*RS(2,4)+206.E0/11025.
     #E0*RS(2,1)*RS(2,2)-206.E0/11025.E0*RS(2,1)*RS(2,4))*RS(1,1)**2
      s3 = ((3.E0/1225.E0*RS(2,3)*RS(2,2)-11.E0/11025.E0*RS(2,3)**2+11.E
     #0/11025.E0*RS(2,4)**2+206.E0/11025.E0*RS(2,2)**2-206.E0/11025.E0*R
     #S(2,1)**2-3.E0/1225.E0*RS(2,1)*RS(2,4)+4.E0/1575.E0*RS(2,1)*RS(2,3
     #)-4.E0/1575.E0*RS(2,4)*RS(2,2))*RS(1,2)+(-4.E0/1575.E0*RS(2,1)*RS(
     #2,2)+4.E0/1575.E0*RS(2,3)*RS(2,2)-4.E0/1575.E0*RS(2,3)*RS(2,4)+4.E
     #0/1575.E0*RS(2,1)*RS(2,4))*RS(1,3)+(-4.E0/1575.E0*RS(2,1)*RS(2,3)-
     #11.E0/11025.E0*RS(2,2)**2+3.E0/1225.E0*RS(2,1)*RS(2,2)-206.E0/1102
     #5.E0*RS(2,4)**2-3.E0/1225.E0*RS(2,3)*RS(2,4)+11.E0/11025.E0*RS(2,3
     #)**2+4.E0/1575.E0*RS(2,4)*RS(2,2)+206.E0/11025.E0*RS(2,1)**2)*RS(1
     #,4))*RS(1,1)+(206.E0/11025.E0*RS(2,3)*RS(2,2)+12.E0/1225.E0*RS(2,3
     #)**2+11.E0/11025.E0*RS(2,1)*RS(2,4)-11.E0/11025.E0*RS(2,3)*RS(2,4)
     #-12.E0/1225.E0*RS(2,1)**2-206.E0/11025.E0*RS(2,1)*RS(2,2))*RS(1,2)
     #**2
      s1 = s2+s3
      s2 = s1+((-4.E0/1575.E0*RS(2,1)*RS(2,3)+4.E0/1575.E0*RS(2,4)*RS(2,
     #2)-3.E0/1225.E0*RS(2,1)*RS(2,2)+3.E0/1225.E0*RS(2,3)*RS(2,4)-206.E
     #0/11025.E0*RS(2,2)**2+11.E0/11025.E0*RS(2,1)**2+206.E0/11025.E0*RS
     #(2,3)**2-11.E0/11025.E0*RS(2,4)**2)*RS(1,3)+(-4.E0/1575.E0*RS(2,3)
     #*RS(2,2)-4.E0/1575.E0*RS(2,1)*RS(2,4)+4.E0/1575.E0*RS(2,3)*RS(2,4)
     #+4.E0/1575.E0*RS(2,1)*RS(2,2))*RS(1,4))*RS(1,2)
      t0 = s2+(-11.E0/11025.E0*RS(2,1)*RS(2,4)-12.E0/1225.E0*RS(2,2)**2+
     #11.E0/11025.E0*RS(2,1)*RS(2,2)+206.E0/11025.E0*RS(2,3)*RS(2,4)+12.
     #E0/1225.E0*RS(2,4)**2-206.E0/11025.E0*RS(2,3)*RS(2,2))*RS(1,3)**2+
     #(-11.E0/11025.E0*RS(2,1)**2+11.E0/11025.E0*RS(2,2)**2-4.E0/1575.E0
     #*RS(2,4)*RS(2,2)+3.E0/1225.E0*RS(2,1)*RS(2,4)+206.E0/11025.E0*RS(2
     #,4)**2+4.E0/1575.E0*RS(2,1)*RS(2,3)-206.E0/11025.E0*RS(2,3)**2-3.E
     #0/1225.E0*RS(2,3)*RS(2,2))*RS(1,4)*RS(1,3)+(12.E0/1225.E0*RS(2,1)*
     #*2-12.E0/1225.E0*RS(2,3)**2+11.E0/11025.E0*RS(2,3)*RS(2,2)-206.E0/
     #11025.E0*RS(2,3)*RS(2,4)-11.E0/11025.E0*RS(2,1)*RS(2,2)+206.E0/110
     #25.E0*RS(2,1)*RS(2,4))*RS(1,4)**2
	ARS(5) = T0

C	-----------------------
C	COMPUTE FOR FACTOR ARRR
C	-----------------------
      s2 = (-RS(2,4)/30+RS(2,2)/30)*RS(1,1)**4+((167.E0/198450.E0*RS(2,4
     #)+RS(2,2)/40-997.E0/793800.E0*RS(2,3)-RS(2,1)/30)*RS(1,2)+(-997.E0
     #/793800.E0*RS(2,4)+997.E0/793800.E0*RS(2,2))*RS(1,3)+(-RS(2,4)/40-
     #167.E0/198450.E0*RS(2,2)+RS(2,1)/30+997.E0/793800.E0*RS(2,3))*RS(1
     #,4))*RS(1,1)**3+((RS(2,2)/60-121.E0/264600.E0*RS(2,4)-149.E0/13230
     #0.E0*RS(2,3)-RS(2,1)/40)*RS(1,2)**2+((-73.E0/529200.E0*RS(2,4)+149
     #.E0/132300.E0*RS(2,2)+17.E0/105840.E0*RS(2,3))*RS(1,3)+(-121.E0/26
     #4600.E0*RS(2,4)+121.E0/264600.E0*RS(2,2))*RS(1,4))*RS(1,2)+(17.E0/
     #105840.E0*RS(2,4)-17.E0/105840.E0*RS(2,2))*RS(1,3)**2+(-17.E0/1058
     #40.E0*RS(2,3)+73.E0/529200.E0*RS(2,2)-149.E0/132300.E0*RS(2,4))*RS
     #(1,4)*RS(1,3)+(RS(2,1)/40+121.E0/264600.E0*RS(2,2)-RS(2,4)/60+149.
     #E0/132300.E0*RS(2,3))*RS(1,4)**2)*RS(1,1)**2
      s4 = s2
      s5 = ((RS(2,2)/120-149.E0/132300.E0*RS(2,4)-121.E0/264600.E0*RS(2,
     #3)-RS(2,1)/60)*RS(1,2)**3+((131.E0/529200.E0*RS(2,3)+121.E0/264600
     #.E0*RS(2,2)-RS(2,4)/25200)*RS(1,3)+(73.E0/529200.E0*RS(2,3)+131.E0
     #/529200.E0*RS(2,4)+149.E0/132300.E0*RS(2,2))*RS(1,4))*RS(1,2)**2+(
     #(17.E0/105840.E0*RS(2,3)+RS(2,4)/17640-131.E0/529200.E0*RS(2,2))*R
     #S(1,3)**2+(-13.E0/132300.E0*RS(2,2)+13.E0/132300.E0*RS(2,4))*RS(1,
     #4)*RS(1,3)+(-73.E0/529200.E0*RS(2,3)-149.E0/132300.E0*RS(2,4)-131.
     #E0/529200.E0*RS(2,2))*RS(1,4)**2)*RS(1,2)+(17.E0/105840.E0*RS(2,4)
     #-17.E0/105840.E0*RS(2,2))*RS(1,3)**3+(-RS(2,2)/17640+131.E0/529200
     #.E0*RS(2,4)-17.E0/105840.E0*RS(2,3))*RS(1,4)*RS(1,3)**2+(-131.E0/5
     #29200.E0*RS(2,3)-121.E0/264600.E0*RS(2,4)+RS(2,2)/25200)*RS(1,4)**
     #2*RS(1,3)+(RS(2,1)/60-RS(2,4)/120+149.E0/132300.E0*RS(2,2)+121.E0/
     #264600.E0*RS(2,3))*RS(1,4)**3)*RS(1,1)
      s3 = s4+s5
      s1 = s3+(167.E0/198450.E0*RS(2,3)-RS(2,1)/120-997.E0/793800.E0*RS(
     #2,4))*RS(1,2)**4+((-167.E0/198450.E0*RS(2,2)-121.E0/264600.E0*RS(2
     #,3)-31.E0/113400.E0*RS(2,4))*RS(1,3)+(43.E0/317520.E0*RS(2,3)+997.
     #E0/793800.E0*RS(2,2)+17.E0/105840.E0*RS(2,4))*RS(1,4))*RS(1,2)**3
      s2 = s1+((121.E0/264600.E0*RS(2,2)-149.E0/132300.E0*RS(2,3)-47.E0/
     #264600.E0*RS(2,4))*RS(1,3)**2+(73.E0/529200.E0*RS(2,3)+RS(2,4)/176
     #40+73.E0/529200.E0*RS(2,2))*RS(1,4)*RS(1,3)+(17.E0/105840.E0*RS(2,
     #4)-17.E0/105840.E0*RS(2,2))*RS(1,4)**2)*RS(1,2)**2+((-997.E0/79380
     #0.E0*RS(2,3)-31.E0/113400.E0*RS(2,4)+149.E0/132300.E0*RS(2,2))*RS(
     #1,3)**3+(RS(2,2)/25200-RS(2,4)/25200)*RS(1,4)*RS(1,3)**2+(-RS(2,2)
     #/17640-73.E0/529200.E0*RS(2,3)-73.E0/529200.E0*RS(2,4))*RS(1,4)**2
     #*RS(1,3)+(-17.E0/105840.E0*RS(2,2)-997.E0/793800.E0*RS(2,4)-43.E0/
     #317520.E0*RS(2,3))*RS(1,4)**3)*RS(1,2)+(-997.E0/793800.E0*RS(2,4)+
     #997.E0/793800.E0*RS(2,2))*RS(1,3)**4
      t0 = s2+(997.E0/793800.E0*RS(2,3)-149.E0/132300.E0*RS(2,4)+31.E0/1
     #13400.E0*RS(2,2))*RS(1,4)*RS(1,3)**3+(47.E0/264600.E0*RS(2,2)+149.
     #E0/132300.E0*RS(2,3)-121.E0/264600.E0*RS(2,4))*RS(1,4)**2*RS(1,3)*
     #*2+(31.E0/113400.E0*RS(2,2)+121.E0/264600.E0*RS(2,3)+167.E0/198450
     #.E0*RS(2,4))*RS(1,4)**3*RS(1,3)+(RS(2,1)/120+997.E0/793800.E0*RS(2
     #,2)-167.E0/198450.E0*RS(2,3))*RS(1,4)**4
	ARRR(1) = T0

      s2 = (997.E0/793800.E0*RS(2,3)+RS(2,2)/120-167.E0/198450.E0*RS(2,4
     #))*RS(1,1)**4+((-RS(2,1)/120+149.E0/132300.E0*RS(2,3)+121.E0/26460
     #0.E0*RS(2,4)+RS(2,2)/60)*RS(1,2)+(-17.E0/105840.E0*RS(2,3)-997.E0/
     #793800.E0*RS(2,1)-43.E0/317520.E0*RS(2,4))*RS(1,3)+(31.E0/113400.E
     #0*RS(2,3)+121.E0/264600.E0*RS(2,4)+167.E0/198450.E0*RS(2,1))*RS(1,
     #4))*RS(1,1)**3+((121.E0/264600.E0*RS(2,3)+149.E0/132300.E0*RS(2,4)
     #+RS(2,2)/40-RS(2,1)/60)*RS(1,2)**2+((-73.E0/529200.E0*RS(2,4)-149.
     #E0/132300.E0*RS(2,1)-131.E0/529200.E0*RS(2,3))*RS(1,3)+(-121.E0/26
     #4600.E0*RS(2,1)-131.E0/529200.E0*RS(2,4)+RS(2,3)/25200)*RS(1,4))*R
     #S(1,2)+(-17.E0/105840.E0*RS(2,3)+17.E0/105840.E0*RS(2,1))*RS(1,3)*
     #*2+(-RS(2,3)/17640-73.E0/529200.E0*RS(2,1)-73.E0/529200.E0*RS(2,4)
     #)*RS(1,4)*RS(1,3)+(-121.E0/264600.E0*RS(2,1)+149.E0/132300.E0*RS(2
     #,4)+47.E0/264600.E0*RS(2,3))*RS(1,4)**2)*RS(1,1)**2
      s3 = s2+((997.E0/793800.E0*RS(2,4)-RS(2,1)/40-167.E0/198450.E0*RS(
     #2,3)+RS(2,2)/30)*RS(1,2)**3+((-121.E0/264600.E0*RS(2,1)+121.E0/264
     #600.E0*RS(2,3))*RS(1,3)+(73.E0/529200.E0*RS(2,3)-149.E0/132300.E0*
     #RS(2,1)-17.E0/105840.E0*RS(2,4))*RS(1,4))*RS(1,2)**2+((149.E0/1323
     #00.E0*RS(2,3)+131.E0/529200.E0*RS(2,1)+73.E0/529200.E0*RS(2,4))*RS
     #(1,3)**2+(13.E0/132300.E0*RS(2,1)-13.E0/132300.E0*RS(2,3))*RS(1,4)
     #*RS(1,3)+(-RS(2,3)/17640-17.E0/105840.E0*RS(2,4)+131.E0/529200.E0*
     #RS(2,1))*RS(1,4)**2)*RS(1,2)+(17.E0/105840.E0*RS(2,1)+997.E0/79380
     #0.E0*RS(2,3)+43.E0/317520.E0*RS(2,4))*RS(1,3)**3+(RS(2,1)/17640+73
     #.E0/529200.E0*RS(2,4)+73.E0/529200.E0*RS(2,3))*RS(1,4)*RS(1,3)**2+
     #(-RS(2,1)/25200+RS(2,3)/25200)*RS(1,4)**2*RS(1,3)+(-149.E0/132300.
     #E0*RS(2,1)+31.E0/113400.E0*RS(2,3)+997.E0/793800.E0*RS(2,4))*RS(1,
     #4)**3)*RS(1,1)
      s1 = s3+(RS(2,3)/30-RS(2,1)/30)*RS(1,2)**4+((167.E0/198450.E0*RS(2
     #,1)-RS(2,2)/30+RS(2,3)/40-997.E0/793800.E0*RS(2,4))*RS(1,3)+(-997.
     #E0/793800.E0*RS(2,1)+997.E0/793800.E0*RS(2,3))*RS(1,4))*RS(1,2)**3
      s2 = s1+((RS(2,3)/60-RS(2,2)/40-121.E0/264600.E0*RS(2,1)-149.E0/13
     #2300.E0*RS(2,4))*RS(1,3)**2+(17.E0/105840.E0*RS(2,4)+149.E0/132300
     #.E0*RS(2,3)-73.E0/529200.E0*RS(2,1))*RS(1,4)*RS(1,3)+(-17.E0/10584
     #0.E0*RS(2,3)+17.E0/105840.E0*RS(2,1))*RS(1,4)**2)*RS(1,2)**2+((-RS
     #(2,2)/60-149.E0/132300.E0*RS(2,1)+RS(2,3)/120-121.E0/264600.E0*RS(
     #2,4))*RS(1,3)**3+(-RS(2,1)/25200+121.E0/264600.E0*RS(2,3)+131.E0/5
     #29200.E0*RS(2,4))*RS(1,4)*RS(1,3)**2+(RS(2,1)/17640+17.E0/105840.E
     #0*RS(2,4)-131.E0/529200.E0*RS(2,3))*RS(1,4)**2*RS(1,3)+(-17.E0/105
     #840.E0*RS(2,3)+17.E0/105840.E0*RS(2,1))*RS(1,4)**3)*RS(1,2)+(167.E
     #0/198450.E0*RS(2,4)-RS(2,2)/120-997.E0/793800.E0*RS(2,1))*RS(1,3)*
     #*4
      t0 = s2+(-167.E0/198450.E0*RS(2,3)-121.E0/264600.E0*RS(2,4)-31.E0/
     #113400.E0*RS(2,1))*RS(1,4)*RS(1,3)**3+(-47.E0/264600.E0*RS(2,1)+12
     #1.E0/264600.E0*RS(2,3)-149.E0/132300.E0*RS(2,4))*RS(1,4)**2*RS(1,3
     #)**2+(149.E0/132300.E0*RS(2,3)-31.E0/113400.E0*RS(2,1)-997.E0/7938
     #00.E0*RS(2,4))*RS(1,4)**3*RS(1,3)+(-997.E0/793800.E0*RS(2,1)+997.E
     #0/793800.E0*RS(2,3))*RS(1,4)**4
	ARRR(2) = T0

      s2 = (997.E0/793800.E0*RS(2,4)-997.E0/793800.E0*RS(2,2))*RS(1,1)**
     #4+((997.E0/793800.E0*RS(2,1)+31.E0/113400.E0*RS(2,4)-149.E0/132300
     #.E0*RS(2,2))*RS(1,2)+(-17.E0/105840.E0*RS(2,4)+17.E0/105840.E0*RS(
     #2,2))*RS(1,3)+(149.E0/132300.E0*RS(2,4)-31.E0/113400.E0*RS(2,2)-99
     #7.E0/793800.E0*RS(2,1))*RS(1,4))*RS(1,1)**3+((149.E0/132300.E0*RS(
     #2,1)-121.E0/264600.E0*RS(2,2)+47.E0/264600.E0*RS(2,4))*RS(1,2)**2+
     #((131.E0/529200.E0*RS(2,2)-RS(2,4)/17640-17.E0/105840.E0*RS(2,1))*
     #RS(1,3)+(-RS(2,2)/25200+RS(2,4)/25200)*RS(1,4))*RS(1,2)+(-17.E0/10
     #5840.E0*RS(2,4)+17.E0/105840.E0*RS(2,2))*RS(1,3)**2+(17.E0/105840.
     #E0*RS(2,1)-131.E0/529200.E0*RS(2,4)+RS(2,2)/17640)*RS(1,4)*RS(1,3)
     #+(121.E0/264600.E0*RS(2,4)-47.E0/264600.E0*RS(2,2)-149.E0/132300.E
     #0*RS(2,1))*RS(1,4)**2)*RS(1,1)**2
      s3 = s2+((31.E0/113400.E0*RS(2,4)+121.E0/264600.E0*RS(2,1)+167.E0/
     #198450.E0*RS(2,2))*RS(1,2)**3+((-121.E0/264600.E0*RS(2,2)-131.E0/5
     #29200.E0*RS(2,1)+RS(2,4)/25200)*RS(1,3)+(-RS(2,4)/17640-73.E0/5292
     #00.E0*RS(2,1)-73.E0/529200.E0*RS(2,2))*RS(1,4))*RS(1,2)**2+((-149.
     #E0/132300.E0*RS(2,2)-17.E0/105840.E0*RS(2,1)+73.E0/529200.E0*RS(2,
     #4))*RS(1,3)**2+(13.E0/132300.E0*RS(2,2)-13.E0/132300.E0*RS(2,4))*R
     #S(1,4)*RS(1,3)+(73.E0/529200.E0*RS(2,4)+73.E0/529200.E0*RS(2,1)+RS
     #(2,2)/17640)*RS(1,4)**2)*RS(1,2)+(997.E0/793800.E0*RS(2,4)-997.E0/
     #793800.E0*RS(2,2))*RS(1,3)**3+(17.E0/105840.E0*RS(2,1)-73.E0/52920
     #0.E0*RS(2,2)+149.E0/132300.E0*RS(2,4))*RS(1,4)*RS(1,3)**2+(-RS(2,2
     #)/25200+131.E0/529200.E0*RS(2,1)+121.E0/264600.E0*RS(2,4))*RS(1,4)
     #**2*RS(1,3)+(-31.E0/113400.E0*RS(2,2)-167.E0/198450.E0*RS(2,4)-121
     #.E0/264600.E0*RS(2,1))*RS(1,4)**3)*RS(1,1)
      s1 = s3+(RS(2,3)/120-167.E0/198450.E0*RS(2,1)+997.E0/793800.E0*RS(
     #2,4))*RS(1,2)**4+((RS(2,3)/60-RS(2,2)/120+121.E0/264600.E0*RS(2,1)
     #+149.E0/132300.E0*RS(2,4))*RS(1,3)+(-997.E0/793800.E0*RS(2,2)-17.E
     #0/105840.E0*RS(2,4)-43.E0/317520.E0*RS(2,1))*RS(1,4))*RS(1,2)**3
      s2 = s1+((-RS(2,2)/60+149.E0/132300.E0*RS(2,1)+RS(2,3)/40+121.E0/2
     #64600.E0*RS(2,4))*RS(1,3)**2+(-73.E0/529200.E0*RS(2,1)-149.E0/1323
     #00.E0*RS(2,2)-131.E0/529200.E0*RS(2,4))*RS(1,4)*RS(1,3)+(-17.E0/10
     #5840.E0*RS(2,4)+17.E0/105840.E0*RS(2,2))*RS(1,4)**2)*RS(1,2)**2+((
     #RS(2,3)/30-RS(2,2)/40-167.E0/198450.E0*RS(2,4)+997.E0/793800.E0*RS
     #(2,1))*RS(1,3)**3+(-121.E0/264600.E0*RS(2,2)+121.E0/264600.E0*RS(2
     #,4))*RS(1,4)*RS(1,3)**2+(73.E0/529200.E0*RS(2,1)+131.E0/529200.E0*
     #RS(2,2)+149.E0/132300.E0*RS(2,4))*RS(1,4)**2*RS(1,3)+(43.E0/317520
     #.E0*RS(2,1)+17.E0/105840.E0*RS(2,2)+997.E0/793800.E0*RS(2,4))*RS(1
     #,4)**3)*RS(1,2)+(-RS(2,2)/30+RS(2,4)/30)*RS(1,3)**4
      t0 = s2+(-997.E0/793800.E0*RS(2,1)-RS(2,3)/30+RS(2,4)/40+167.E0/19
     #8450.E0*RS(2,2))*RS(1,4)*RS(1,3)**3+(-149.E0/132300.E0*RS(2,1)-121
     #.E0/264600.E0*RS(2,2)+RS(2,4)/60-RS(2,3)/40)*RS(1,4)**2*RS(1,3)**2
     #+(-149.E0/132300.E0*RS(2,2)-RS(2,3)/60-121.E0/264600.E0*RS(2,1)+RS
     #(2,4)/120)*RS(1,4)**3*RS(1,3)+(-997.E0/793800.E0*RS(2,2)+167.E0/19
     #8450.E0*RS(2,1)-RS(2,3)/120)*RS(1,4)**4
	ARRR(3) = T0

      s2 = (-RS(2,4)/120+167.E0/198450.E0*RS(2,2)-997.E0/793800.E0*RS(2,
     #3))*RS(1,1)**4+((-167.E0/198450.E0*RS(2,1)-31.E0/113400.E0*RS(2,3)
     #-121.E0/264600.E0*RS(2,2))*RS(1,2)+(997.E0/793800.E0*RS(2,1)+43.E0
     #/317520.E0*RS(2,2)+17.E0/105840.E0*RS(2,3))*RS(1,3)+(-121.E0/26460
     #0.E0*RS(2,2)-149.E0/132300.E0*RS(2,3)-RS(2,4)/60+RS(2,1)/120)*RS(1
     #,4))*RS(1,1)**3+((-47.E0/264600.E0*RS(2,3)-149.E0/132300.E0*RS(2,2
     #)+121.E0/264600.E0*RS(2,1))*RS(1,2)**2+((73.E0/529200.E0*RS(2,2)+R
     #S(2,3)/17640+73.E0/529200.E0*RS(2,1))*RS(1,3)+(121.E0/264600.E0*RS
     #(2,1)-RS(2,3)/25200+131.E0/529200.E0*RS(2,2))*RS(1,4))*RS(1,2)+(17
     #.E0/105840.E0*RS(2,3)-17.E0/105840.E0*RS(2,1))*RS(1,3)**2+(73.E0/5
     #29200.E0*RS(2,2)+149.E0/132300.E0*RS(2,1)+131.E0/529200.E0*RS(2,3)
     #)*RS(1,4)*RS(1,3)+(-121.E0/264600.E0*RS(2,3)-149.E0/132300.E0*RS(2
     #,2)+RS(2,1)/60-RS(2,4)/40)*RS(1,4)**2)*RS(1,1)**2
      s3 = s2+((-31.E0/113400.E0*RS(2,3)-997.E0/793800.E0*RS(2,2)+149.E0
     #/132300.E0*RS(2,1))*RS(1,2)**3+((RS(2,1)/25200-RS(2,3)/25200)*RS(1
     #,3)+(RS(2,3)/17640+17.E0/105840.E0*RS(2,2)-131.E0/529200.E0*RS(2,1
     #))*RS(1,4))*RS(1,2)**2+((-RS(2,1)/17640-73.E0/529200.E0*RS(2,3)-73
     #.E0/529200.E0*RS(2,2))*RS(1,3)**2+(13.E0/132300.E0*RS(2,3)-13.E0/1
     #32300.E0*RS(2,1))*RS(1,4)*RS(1,3)+(-73.E0/529200.E0*RS(2,3)+17.E0/
     #105840.E0*RS(2,2)+149.E0/132300.E0*RS(2,1))*RS(1,4)**2)*RS(1,2)+(-
     #17.E0/105840.E0*RS(2,1)-43.E0/317520.E0*RS(2,2)-997.E0/793800.E0*R
     #S(2,3))*RS(1,3)**3+(-149.E0/132300.E0*RS(2,3)-73.E0/529200.E0*RS(2
     #,2)-131.E0/529200.E0*RS(2,1))*RS(1,4)*RS(1,3)**2+(121.E0/264600.E0
     #*RS(2,1)-121.E0/264600.E0*RS(2,3))*RS(1,4)**2*RS(1,3)+(-RS(2,4)/30
     #-997.E0/793800.E0*RS(2,2)+167.E0/198450.E0*RS(2,3)+RS(2,1)/40)*RS(
     #1,4)**3)*RS(1,1)
      s1 = s3+(997.E0/793800.E0*RS(2,1)-997.E0/793800.E0*RS(2,3))*RS(1,2
     #)**4+((-149.E0/132300.E0*RS(2,3)+31.E0/113400.E0*RS(2,1)+997.E0/79
     #3800.E0*RS(2,2))*RS(1,3)+(17.E0/105840.E0*RS(2,3)-17.E0/105840.E0*
     #RS(2,1))*RS(1,4))*RS(1,2)**3
      s2 = s1+((-121.E0/264600.E0*RS(2,3)+149.E0/132300.E0*RS(2,2)+47.E0
     #/264600.E0*RS(2,1))*RS(1,3)**2+(-RS(2,1)/17640+131.E0/529200.E0*RS
     #(2,3)-17.E0/105840.E0*RS(2,2))*RS(1,4)*RS(1,3)+(17.E0/105840.E0*RS
     #(2,3)-17.E0/105840.E0*RS(2,1))*RS(1,4)**2)*RS(1,2)**2+((167.E0/198
     #450.E0*RS(2,3)+31.E0/113400.E0*RS(2,1)+121.E0/264600.E0*RS(2,2))*R
     #S(1,3)**3+(-131.E0/529200.E0*RS(2,2)+RS(2,1)/25200-121.E0/264600.E
     #0*RS(2,3))*RS(1,4)*RS(1,3)**2+(-17.E0/105840.E0*RS(2,2)-149.E0/132
     #300.E0*RS(2,3)+73.E0/529200.E0*RS(2,1))*RS(1,4)**2*RS(1,3)+(997.E0
     #/793800.E0*RS(2,1)-997.E0/793800.E0*RS(2,3))*RS(1,4)**3)*RS(1,2)+(
     #-167.E0/198450.E0*RS(2,2)+997.E0/793800.E0*RS(2,1)+RS(2,4)/120)*RS
     #(1,3)**4
      t0 = s2+(149.E0/132300.E0*RS(2,1)+RS(2,4)/60+121.E0/264600.E0*RS(2
     #,2)-RS(2,3)/120)*RS(1,4)*RS(1,3)**3+(149.E0/132300.E0*RS(2,2)+RS(2
     #,4)/40-RS(2,3)/60+121.E0/264600.E0*RS(2,1))*RS(1,4)**2*RS(1,3)**2+
     #(-167.E0/198450.E0*RS(2,1)+997.E0/793800.E0*RS(2,2)+RS(2,4)/30-RS(
     #2,3)/40)*RS(1,4)**3*RS(1,3)+(-RS(2,3)/30+RS(2,1)/30)*RS(1,4)**4
	ARRR(4) = T0

      s2 = (124.E0/14175.E0*RS(2,2)-124.E0/14175.E0*RS(2,4))*RS(1,1)**4+
     #((-124.E0/14175.E0*RS(2,1)+328.E0/33075.E0*RS(2,2)+8.E0/19845.E0*R
     #S(2,3)-52.E0/33075.E0*RS(2,4))*RS(1,2)+(-22.E0/14175.E0*RS(2,2)+22
     #.E0/14175.E0*RS(2,4))*RS(1,3)+(-8.E0/19845.E0*RS(2,3)+124.E0/14175
     #.E0*RS(2,1)+52.E0/33075.E0*RS(2,2)-328.E0/33075.E0*RS(2,4))*RS(1,4
     #))*RS(1,1)**3+((-328.E0/33075.E0*RS(2,1)+4.E0/4725.E0*RS(2,3)-4.E0
     #/4725.E0*RS(2,4)+328.E0/33075.E0*RS(2,2))*RS(1,2)**2+((38.E0/33075
     #.E0*RS(2,1)+RS(2,3)/33075-2.E0/1323.E0*RS(2,2)+11.E0/33075.E0*RS(2
     #,4))*RS(1,3)+(22.E0/33075.E0*RS(2,4)-22.E0/33075.E0*RS(2,2))*RS(1,
     #4))*RS(1,2)+(2.E0/1323.E0*RS(2,4)-11.E0/33075.E0*RS(2,2)-RS(2,3)/3
     #3075-38.E0/33075.E0*RS(2,1))*RS(1,4)*RS(1,3)+(-328.E0/33075.E0*RS(
     #2,4)+328.E0/33075.E0*RS(2,1)-4.E0/4725.E0*RS(2,3)+4.E0/4725.E0*RS(
     #2,2))*RS(1,4)**2)*RS(1,1)**2
      s3 = s2+((-328.E0/33075.E0*RS(2,1)+124.E0/14175.E0*RS(2,2)+52.E0/3
     #3075.E0*RS(2,3)-8.E0/19845.E0*RS(2,4))*RS(1,2)**3+((-22.E0/33075.E
     #0*RS(2,3)+22.E0/33075.E0*RS(2,1))*RS(1,3)+(-11.E0/33075.E0*RS(2,3)
     #+2.E0/1323.E0*RS(2,1)-RS(2,4)/33075-38.E0/33075.E0*RS(2,2))*RS(1,4
     #))*RS(1,2)**2+((-RS(2,1)/33075-11.E0/33075.E0*RS(2,4)+2.E0/1323.E0
     #*RS(2,2)-38.E0/33075.E0*RS(2,3))*RS(1,3)**2+(11.E0/33075.E0*RS(2,3
     #)+38.E0/33075.E0*RS(2,4)-2.E0/1323.E0*RS(2,1)+RS(2,2)/33075)*RS(1,
     #4)**2)*RS(1,2)+(22.E0/14175.E0*RS(2,2)-22.E0/14175.E0*RS(2,4))*RS(
     #1,3)**3+(-2.E0/1323.E0*RS(2,4)+38.E0/33075.E0*RS(2,3)+RS(2,1)/3307
     #5+11.E0/33075.E0*RS(2,2))*RS(1,4)*RS(1,3)**2+(-22.E0/33075.E0*RS(2
     #,1)+22.E0/33075.E0*RS(2,3))*RS(1,4)**2*RS(1,3)+(-52.E0/33075.E0*RS
     #(2,3)+328.E0/33075.E0*RS(2,1)-124.E0/14175.E0*RS(2,4)+8.E0/19845.E
     #0*RS(2,2))*RS(1,4)**3)*RS(1,1)
      s1 = s3+(124.E0/14175.E0*RS(2,3)-124.E0/14175.E0*RS(2,1))*RS(1,2)*
     #*4+((-124.E0/14175.E0*RS(2,2)+328.E0/33075.E0*RS(2,3)+8.E0/19845.E
     #0*RS(2,4)-52.E0/33075.E0*RS(2,1))*RS(1,3)+(22.E0/14175.E0*RS(2,1)-
     #22.E0/14175.E0*RS(2,3))*RS(1,4))*RS(1,2)**3
      s2 = s1+((-328.E0/33075.E0*RS(2,2)+328.E0/33075.E0*RS(2,3)-4.E0/47
     #25.E0*RS(2,1)+4.E0/4725.E0*RS(2,4))*RS(1,3)**2+(38.E0/33075.E0*RS(
     #2,2)-2.E0/1323.E0*RS(2,3)+RS(2,4)/33075+11.E0/33075.E0*RS(2,1))*RS
     #(1,4)*RS(1,3))*RS(1,2)**2+((-8.E0/19845.E0*RS(2,1)-328.E0/33075.E0
     #*RS(2,2)+124.E0/14175.E0*RS(2,3)+52.E0/33075.E0*RS(2,4))*RS(1,3)**
     #3+(22.E0/33075.E0*RS(2,2)-22.E0/33075.E0*RS(2,4))*RS(1,4)*RS(1,3)*
     #*2+(-11.E0/33075.E0*RS(2,1)+2.E0/1323.E0*RS(2,3)-38.E0/33075.E0*RS
     #(2,4)-RS(2,2)/33075)*RS(1,4)**2*RS(1,3)+(22.E0/14175.E0*RS(2,3)-22
     #.E0/14175.E0*RS(2,1))*RS(1,4)**3)*RS(1,2)+(-124.E0/14175.E0*RS(2,2
     #)+124.E0/14175.E0*RS(2,4))*RS(1,3)**4
      t0 = s2+(328.E0/33075.E0*RS(2,4)-124.E0/14175.E0*RS(2,3)+8.E0/1984
     #5.E0*RS(2,1)-52.E0/33075.E0*RS(2,2))*RS(1,4)*RS(1,3)**3+(328.E0/33
     #075.E0*RS(2,4)+4.E0/4725.E0*RS(2,1)-328.E0/33075.E0*RS(2,3)-4.E0/4
     #725.E0*RS(2,2))*RS(1,4)**2*RS(1,3)**2+(52.E0/33075.E0*RS(2,1)-328.
     #E0/33075.E0*RS(2,3)+124.E0/14175.E0*RS(2,4)-8.E0/19845.E0*RS(2,2))
     #*RS(1,4)**3*RS(1,3)+(124.E0/14175.E0*RS(2,1)-124.E0/14175.E0*RS(2,
     #3))*RS(1,4)**4
	ARRR(5) = T0

C	-----------------------
C	COMPUTE FOR FACTOR ASSS
C	-----------------------
      s2 = (RS(1,4)/30-RS(1,2)/30)*RS(2,1)**4+((997.E0/793800.E0*RS(1,3)
     #-RS(1,2)/40+RS(1,1)/30-167.E0/198450.E0*RS(1,4))*RS(2,2)+(997.E0/7
     #93800.E0*RS(1,4)-997.E0/793800.E0*RS(1,2))*RS(2,3)+(167.E0/198450.
     #E0*RS(1,2)+RS(1,4)/40-997.E0/793800.E0*RS(1,3)-RS(1,1)/30)*RS(2,4)
     #)*RS(2,1)**3+((121.E0/264600.E0*RS(1,4)+RS(1,1)/40+149.E0/132300.E
     #0*RS(1,3)-RS(1,2)/60)*RS(2,2)**2+((73.E0/529200.E0*RS(1,4)-17.E0/1
     #05840.E0*RS(1,3)-149.E0/132300.E0*RS(1,2))*RS(2,3)+(121.E0/264600.
     #E0*RS(1,4)-121.E0/264600.E0*RS(1,2))*RS(2,4))*RS(2,2)+(-17.E0/1058
     #40.E0*RS(1,4)+17.E0/105840.E0*RS(1,2))*RS(2,3)**2+(-73.E0/529200.E
     #0*RS(1,2)+149.E0/132300.E0*RS(1,4)+17.E0/105840.E0*RS(1,3))*RS(2,4
     #)*RS(2,3)+(-RS(1,1)/40+RS(1,4)/60-149.E0/132300.E0*RS(1,3)-121.E0/
     #264600.E0*RS(1,2))*RS(2,4)**2)*RS(2,1)**2
      s4 = s2
      s5 = ((RS(1,1)/60+121.E0/264600.E0*RS(1,3)-RS(1,2)/120+149.E0/1323
     #00.E0*RS(1,4))*RS(2,2)**3+((-121.E0/264600.E0*RS(1,2)+RS(1,4)/2520
     #0-131.E0/529200.E0*RS(1,3))*RS(2,3)+(-149.E0/132300.E0*RS(1,2)-131
     #.E0/529200.E0*RS(1,4)-73.E0/529200.E0*RS(1,3))*RS(2,4))*RS(2,2)**2
     #+((131.E0/529200.E0*RS(1,2)-17.E0/105840.E0*RS(1,3)-RS(1,4)/17640)
     #*RS(2,3)**2+(-13.E0/132300.E0*RS(1,4)+13.E0/132300.E0*RS(1,2))*RS(
     #2,4)*RS(2,3)+(73.E0/529200.E0*RS(1,3)+149.E0/132300.E0*RS(1,4)+131
     #.E0/529200.E0*RS(1,2))*RS(2,4)**2)*RS(2,2)+(-17.E0/105840.E0*RS(1,
     #4)+17.E0/105840.E0*RS(1,2))*RS(2,3)**3+(-131.E0/529200.E0*RS(1,4)+
     #RS(1,2)/17640+17.E0/105840.E0*RS(1,3))*RS(2,4)*RS(2,3)**2+(-RS(1,2
     #)/25200+131.E0/529200.E0*RS(1,3)+121.E0/264600.E0*RS(1,4))*RS(2,4)
     #**2*RS(2,3)+(-RS(1,1)/60-149.E0/132300.E0*RS(1,2)-121.E0/264600.E0
     #*RS(1,3)+RS(1,4)/120)*RS(2,4)**3)*RS(2,1)
      s3 = s4+s5
      s1 = s3+(-167.E0/198450.E0*RS(1,3)+RS(1,1)/120+997.E0/793800.E0*RS
     #(1,4))*RS(2,2)**4+((121.E0/264600.E0*RS(1,3)+167.E0/198450.E0*RS(1
     #,2)+31.E0/113400.E0*RS(1,4))*RS(2,3)+(-43.E0/317520.E0*RS(1,3)-997
     #.E0/793800.E0*RS(1,2)-17.E0/105840.E0*RS(1,4))*RS(2,4))*RS(2,2)**3
      s2 = s1+((149.E0/132300.E0*RS(1,3)+47.E0/264600.E0*RS(1,4)-121.E0/
     #264600.E0*RS(1,2))*RS(2,3)**2+(-73.E0/529200.E0*RS(1,2)-73.E0/5292
     #00.E0*RS(1,3)-RS(1,4)/17640)*RS(2,4)*RS(2,3)+(-17.E0/105840.E0*RS(
     #1,4)+17.E0/105840.E0*RS(1,2))*RS(2,4)**2)*RS(2,2)**2+((-149.E0/132
     #300.E0*RS(1,2)+31.E0/113400.E0*RS(1,4)+997.E0/793800.E0*RS(1,3))*R
     #S(2,3)**3+(RS(1,4)/25200-RS(1,2)/25200)*RS(2,4)*RS(2,3)**2+(73.E0/
     #529200.E0*RS(1,4)+RS(1,2)/17640+73.E0/529200.E0*RS(1,3))*RS(2,4)**
     #2*RS(2,3)+(997.E0/793800.E0*RS(1,4)+17.E0/105840.E0*RS(1,2)+43.E0/
     #317520.E0*RS(1,3))*RS(2,4)**3)*RS(2,2)+(997.E0/793800.E0*RS(1,4)-9
     #97.E0/793800.E0*RS(1,2))*RS(2,3)**4
      t0 = s2+(-31.E0/113400.E0*RS(1,2)+149.E0/132300.E0*RS(1,4)-997.E0/
     #793800.E0*RS(1,3))*RS(2,4)*RS(2,3)**3+(-149.E0/132300.E0*RS(1,3)+1
     #21.E0/264600.E0*RS(1,4)-47.E0/264600.E0*RS(1,2))*RS(2,4)**2*RS(2,3
     #)**2+(-31.E0/113400.E0*RS(1,2)-167.E0/198450.E0*RS(1,4)-121.E0/264
     #600.E0*RS(1,3))*RS(2,4)**3*RS(2,3)+(-RS(1,1)/120+167.E0/198450.E0*
     #RS(1,3)-997.E0/793800.E0*RS(1,2))*RS(2,4)**4
	ASSS(1) = T0

      s2 = (-RS(1,2)/120+167.E0/198450.E0*RS(1,4)-997.E0/793800.E0*RS(1,
     #3))*RS(2,1)**4+((-121.E0/264600.E0*RS(1,4)-RS(1,2)/60+RS(1,1)/120-
     #149.E0/132300.E0*RS(1,3))*RS(2,2)+(997.E0/793800.E0*RS(1,1)+17.E0/
     #105840.E0*RS(1,3)+43.E0/317520.E0*RS(1,4))*RS(2,3)+(-167.E0/198450
     #.E0*RS(1,1)-31.E0/113400.E0*RS(1,3)-121.E0/264600.E0*RS(1,4))*RS(2
     #,4))*RS(2,1)**3+((RS(1,1)/60-RS(1,2)/40-149.E0/132300.E0*RS(1,4)-1
     #21.E0/264600.E0*RS(1,3))*RS(2,2)**2+((131.E0/529200.E0*RS(1,3)+73.
     #E0/529200.E0*RS(1,4)+149.E0/132300.E0*RS(1,1))*RS(2,3)+(121.E0/264
     #600.E0*RS(1,1)+131.E0/529200.E0*RS(1,4)-RS(1,3)/25200)*RS(2,4))*RS
     #(2,2)+(-17.E0/105840.E0*RS(1,1)+17.E0/105840.E0*RS(1,3))*RS(2,3)**
     #2+(73.E0/529200.E0*RS(1,4)+RS(1,3)/17640+73.E0/529200.E0*RS(1,1))*
     #RS(2,4)*RS(2,3)+(-47.E0/264600.E0*RS(1,3)+121.E0/264600.E0*RS(1,1)
     #-149.E0/132300.E0*RS(1,4))*RS(2,4)**2)*RS(2,1)**2
      s3 = s2+((-997.E0/793800.E0*RS(1,4)-RS(1,2)/30+RS(1,1)/40+167.E0/1
     #98450.E0*RS(1,3))*RS(2,2)**3+((121.E0/264600.E0*RS(1,1)-121.E0/264
     #600.E0*RS(1,3))*RS(2,3)+(149.E0/132300.E0*RS(1,1)-73.E0/529200.E0*
     #RS(1,3)+17.E0/105840.E0*RS(1,4))*RS(2,4))*RS(2,2)**2+((-73.E0/5292
     #00.E0*RS(1,4)-131.E0/529200.E0*RS(1,1)-149.E0/132300.E0*RS(1,3))*R
     #S(2,3)**2+(13.E0/132300.E0*RS(1,3)-13.E0/132300.E0*RS(1,1))*RS(2,4
     #)*RS(2,3)+(RS(1,3)/17640-131.E0/529200.E0*RS(1,1)+17.E0/105840.E0*
     #RS(1,4))*RS(2,4)**2)*RS(2,2)+(-997.E0/793800.E0*RS(1,3)-43.E0/3175
     #20.E0*RS(1,4)-17.E0/105840.E0*RS(1,1))*RS(2,3)**3+(-RS(1,1)/17640-
     #73.E0/529200.E0*RS(1,4)-73.E0/529200.E0*RS(1,3))*RS(2,4)*RS(2,3)**
     #2+(RS(1,1)/25200-RS(1,3)/25200)*RS(2,4)**2*RS(2,3)+(149.E0/132300.
     #E0*RS(1,1)-31.E0/113400.E0*RS(1,3)-997.E0/793800.E0*RS(1,4))*RS(2,
     #4)**3)*RS(2,1)
      s1 = s3+(-RS(1,3)/30+RS(1,1)/30)*RS(2,2)**4+((-RS(1,3)/40-167.E0/1
     #98450.E0*RS(1,1)+RS(1,2)/30+997.E0/793800.E0*RS(1,4))*RS(2,3)+(-99
     #7.E0/793800.E0*RS(1,3)+997.E0/793800.E0*RS(1,1))*RS(2,4))*RS(2,2)*
     #*3
      s2 = s1+((121.E0/264600.E0*RS(1,1)+149.E0/132300.E0*RS(1,4)+RS(1,2
     #)/40-RS(1,3)/60)*RS(2,3)**2+(-17.E0/105840.E0*RS(1,4)-149.E0/13230
     #0.E0*RS(1,3)+73.E0/529200.E0*RS(1,1))*RS(2,4)*RS(2,3)+(-17.E0/1058
     #40.E0*RS(1,1)+17.E0/105840.E0*RS(1,3))*RS(2,4)**2)*RS(2,2)**2+((12
     #1.E0/264600.E0*RS(1,4)+RS(1,2)/60-RS(1,3)/120+149.E0/132300.E0*RS(
     #1,1))*RS(2,3)**3+(-131.E0/529200.E0*RS(1,4)-121.E0/264600.E0*RS(1,
     #3)+RS(1,1)/25200)*RS(2,4)*RS(2,3)**2+(-RS(1,1)/17640+131.E0/529200
     #.E0*RS(1,3)-17.E0/105840.E0*RS(1,4))*RS(2,4)**2*RS(2,3)+(-17.E0/10
     #5840.E0*RS(1,1)+17.E0/105840.E0*RS(1,3))*RS(2,4)**3)*RS(2,2)+(997.
     #E0/793800.E0*RS(1,1)+RS(1,2)/120-167.E0/198450.E0*RS(1,4))*RS(2,3)
     #**4
      t0 = s2+(167.E0/198450.E0*RS(1,3)+31.E0/113400.E0*RS(1,1)+121.E0/2
     #64600.E0*RS(1,4))*RS(2,4)*RS(2,3)**3+(47.E0/264600.E0*RS(1,1)-121.
     #E0/264600.E0*RS(1,3)+149.E0/132300.E0*RS(1,4))*RS(2,4)**2*RS(2,3)*
     #*2+(-149.E0/132300.E0*RS(1,3)+997.E0/793800.E0*RS(1,4)+31.E0/11340
     #0.E0*RS(1,1))*RS(2,4)**3*RS(2,3)+(-997.E0/793800.E0*RS(1,3)+997.E0
     #/793800.E0*RS(1,1))*RS(2,4)**4
	ASSS(2) = T0

      s2 = (997.E0/793800.E0*RS(1,2)-997.E0/793800.E0*RS(1,4))*RS(2,1)**
     #4+((149.E0/132300.E0*RS(1,2)-31.E0/113400.E0*RS(1,4)-997.E0/793800
     #.E0*RS(1,1))*RS(2,2)+(-17.E0/105840.E0*RS(1,2)+17.E0/105840.E0*RS(
     #1,4))*RS(2,3)+(-149.E0/132300.E0*RS(1,4)+31.E0/113400.E0*RS(1,2)+9
     #97.E0/793800.E0*RS(1,1))*RS(2,4))*RS(2,1)**3+((-149.E0/132300.E0*R
     #S(1,1)+121.E0/264600.E0*RS(1,2)-47.E0/264600.E0*RS(1,4))*RS(2,2)**
     #2+((17.E0/105840.E0*RS(1,1)+RS(1,4)/17640-131.E0/529200.E0*RS(1,2)
     #)*RS(2,3)+(RS(1,2)/25200-RS(1,4)/25200)*RS(2,4))*RS(2,2)+(-17.E0/1
     #05840.E0*RS(1,2)+17.E0/105840.E0*RS(1,4))*RS(2,3)**2+(-RS(1,2)/176
     #40+131.E0/529200.E0*RS(1,4)-17.E0/105840.E0*RS(1,1))*RS(2,4)*RS(2,
     #3)+(47.E0/264600.E0*RS(1,2)-121.E0/264600.E0*RS(1,4)+149.E0/132300
     #.E0*RS(1,1))*RS(2,4)**2)*RS(2,1)**2
      s3 = s2+((-121.E0/264600.E0*RS(1,1)-31.E0/113400.E0*RS(1,4)-167.E0
     #/198450.E0*RS(1,2))*RS(2,2)**3+((131.E0/529200.E0*RS(1,1)+121.E0/2
     #64600.E0*RS(1,2)-RS(1,4)/25200)*RS(2,3)+(73.E0/529200.E0*RS(1,2)+R
     #S(1,4)/17640+73.E0/529200.E0*RS(1,1))*RS(2,4))*RS(2,2)**2+((-73.E0
     #/529200.E0*RS(1,4)+149.E0/132300.E0*RS(1,2)+17.E0/105840.E0*RS(1,1
     #))*RS(2,3)**2+(13.E0/132300.E0*RS(1,4)-13.E0/132300.E0*RS(1,2))*RS
     #(2,4)*RS(2,3)+(-73.E0/529200.E0*RS(1,4)-73.E0/529200.E0*RS(1,1)-RS
     #(1,2)/17640)*RS(2,4)**2)*RS(2,2)+(997.E0/793800.E0*RS(1,2)-997.E0/
     #793800.E0*RS(1,4))*RS(2,3)**3+(73.E0/529200.E0*RS(1,2)-149.E0/1323
     #00.E0*RS(1,4)-17.E0/105840.E0*RS(1,1))*RS(2,4)*RS(2,3)**2+(-121.E0
     #/264600.E0*RS(1,4)+RS(1,2)/25200-131.E0/529200.E0*RS(1,1))*RS(2,4)
     #**2*RS(2,3)+(167.E0/198450.E0*RS(1,4)+121.E0/264600.E0*RS(1,1)+31.
     #E0/113400.E0*RS(1,2))*RS(2,4)**3)*RS(2,1)
      s1 = s3+(-997.E0/793800.E0*RS(1,4)-RS(1,3)/120+167.E0/198450.E0*RS
     #(1,1))*RS(2,2)**4+((-149.E0/132300.E0*RS(1,4)-121.E0/264600.E0*RS(
     #1,1)-RS(1,3)/60+RS(1,2)/120)*RS(2,3)+(997.E0/793800.E0*RS(1,2)+43.
     #E0/317520.E0*RS(1,1)+17.E0/105840.E0*RS(1,4))*RS(2,4))*RS(2,2)**3
      s2 = s1+((-149.E0/132300.E0*RS(1,1)-RS(1,3)/40+RS(1,2)/60-121.E0/2
     #64600.E0*RS(1,4))*RS(2,3)**2+(131.E0/529200.E0*RS(1,4)+149.E0/1323
     #00.E0*RS(1,2)+73.E0/529200.E0*RS(1,1))*RS(2,4)*RS(2,3)+(-17.E0/105
     #840.E0*RS(1,2)+17.E0/105840.E0*RS(1,4))*RS(2,4)**2)*RS(2,2)**2+((1
     #67.E0/198450.E0*RS(1,4)-997.E0/793800.E0*RS(1,1)+RS(1,2)/40-RS(1,3
     #)/30)*RS(2,3)**3+(121.E0/264600.E0*RS(1,2)-121.E0/264600.E0*RS(1,4
     #))*RS(2,4)*RS(2,3)**2+(-149.E0/132300.E0*RS(1,4)-73.E0/529200.E0*R
     #S(1,1)-131.E0/529200.E0*RS(1,2))*RS(2,4)**2*RS(2,3)+(-43.E0/317520
     #.E0*RS(1,1)-17.E0/105840.E0*RS(1,2)-997.E0/793800.E0*RS(1,4))*RS(2
     #,4)**3)*RS(2,2)+(-RS(1,4)/30+RS(1,2)/30)*RS(2,3)**4
      t0 = s2+(RS(1,3)/30+997.E0/793800.E0*RS(1,1)-RS(1,4)/40-167.E0/198
     #450.E0*RS(1,2))*RS(2,4)*RS(2,3)**3+(121.E0/264600.E0*RS(1,2)+149.E
     #0/132300.E0*RS(1,1)-RS(1,4)/60+RS(1,3)/40)*RS(2,4)**2*RS(2,3)**2+(
     #RS(1,3)/60-RS(1,4)/120+121.E0/264600.E0*RS(1,1)+149.E0/132300.E0*R
     #S(1,2))*RS(2,4)**3*RS(2,3)+(RS(1,3)/120+997.E0/793800.E0*RS(1,2)-1
     #67.E0/198450.E0*RS(1,1))*RS(2,4)**4
	ASSS(3) = T0

      s2 = (RS(1,4)/120+997.E0/793800.E0*RS(1,3)-167.E0/198450.E0*RS(1,2
     #))*RS(2,1)**4+((167.E0/198450.E0*RS(1,1)+121.E0/264600.E0*RS(1,2)+
     #31.E0/113400.E0*RS(1,3))*RS(2,2)+(-43.E0/317520.E0*RS(1,2)-17.E0/1
     #05840.E0*RS(1,3)-997.E0/793800.E0*RS(1,1))*RS(2,3)+(121.E0/264600.
     #E0*RS(1,2)+RS(1,4)/60+149.E0/132300.E0*RS(1,3)-RS(1,1)/120)*RS(2,4
     #))*RS(2,1)**3+((47.E0/264600.E0*RS(1,3)-121.E0/264600.E0*RS(1,1)+1
     #49.E0/132300.E0*RS(1,2))*RS(2,2)**2+((-73.E0/529200.E0*RS(1,1)-73.
     #E0/529200.E0*RS(1,2)-RS(1,3)/17640)*RS(2,3)+(RS(1,3)/25200-131.E0/
     #529200.E0*RS(1,2)-121.E0/264600.E0*RS(1,1))*RS(2,4))*RS(2,2)+(17.E
     #0/105840.E0*RS(1,1)-17.E0/105840.E0*RS(1,3))*RS(2,3)**2+(-131.E0/5
     #29200.E0*RS(1,3)-73.E0/529200.E0*RS(1,2)-149.E0/132300.E0*RS(1,1))
     #*RS(2,4)*RS(2,3)+(-RS(1,1)/60+121.E0/264600.E0*RS(1,3)+149.E0/1323
     #00.E0*RS(1,2)+RS(1,4)/40)*RS(2,4)**2)*RS(2,1)**2
      s3 = s2+((-149.E0/132300.E0*RS(1,1)+997.E0/793800.E0*RS(1,2)+31.E0
     #/113400.E0*RS(1,3))*RS(2,2)**3+((-RS(1,1)/25200+RS(1,3)/25200)*RS(
     #2,3)+(131.E0/529200.E0*RS(1,1)-RS(1,3)/17640-17.E0/105840.E0*RS(1,
     #2))*RS(2,4))*RS(2,2)**2+((RS(1,1)/17640+73.E0/529200.E0*RS(1,2)+73
     #.E0/529200.E0*RS(1,3))*RS(2,3)**2+(-13.E0/132300.E0*RS(1,3)+13.E0/
     #132300.E0*RS(1,1))*RS(2,4)*RS(2,3)+(73.E0/529200.E0*RS(1,3)-17.E0/
     #105840.E0*RS(1,2)-149.E0/132300.E0*RS(1,1))*RS(2,4)**2)*RS(2,2)+(4
     #3.E0/317520.E0*RS(1,2)+997.E0/793800.E0*RS(1,3)+17.E0/105840.E0*RS
     #(1,1))*RS(2,3)**3+(73.E0/529200.E0*RS(1,2)+149.E0/132300.E0*RS(1,3
     #)+131.E0/529200.E0*RS(1,1))*RS(2,4)*RS(2,3)**2+(-121.E0/264600.E0*
     #RS(1,1)+121.E0/264600.E0*RS(1,3))*RS(2,4)**2*RS(2,3)+(-167.E0/1984
     #50.E0*RS(1,3)+997.E0/793800.E0*RS(1,2)-RS(1,1)/40+RS(1,4)/30)*RS(2
     #,4)**3)*RS(2,1)
      s1 = s3+(-997.E0/793800.E0*RS(1,1)+997.E0/793800.E0*RS(1,3))*RS(2,
     #2)**4+((149.E0/132300.E0*RS(1,3)-997.E0/793800.E0*RS(1,2)-31.E0/11
     #3400.E0*RS(1,1))*RS(2,3)+(17.E0/105840.E0*RS(1,1)-17.E0/105840.E0*
     #RS(1,3))*RS(2,4))*RS(2,2)**3
      s2 = s1+((121.E0/264600.E0*RS(1,3)-47.E0/264600.E0*RS(1,1)-149.E0/
     #132300.E0*RS(1,2))*RS(2,3)**2+(-131.E0/529200.E0*RS(1,3)+17.E0/105
     #840.E0*RS(1,2)+RS(1,1)/17640)*RS(2,4)*RS(2,3)+(17.E0/105840.E0*RS(
     #1,1)-17.E0/105840.E0*RS(1,3))*RS(2,4)**2)*RS(2,2)**2+((-31.E0/1134
     #00.E0*RS(1,1)-121.E0/264600.E0*RS(1,2)-167.E0/198450.E0*RS(1,3))*R
     #S(2,3)**3+(-RS(1,1)/25200+131.E0/529200.E0*RS(1,2)+121.E0/264600.E
     #0*RS(1,3))*RS(2,4)*RS(2,3)**2+(149.E0/132300.E0*RS(1,3)+17.E0/1058
     #40.E0*RS(1,2)-73.E0/529200.E0*RS(1,1))*RS(2,4)**2*RS(2,3)+(-997.E0
     #/793800.E0*RS(1,1)+997.E0/793800.E0*RS(1,3))*RS(2,4)**3)*RS(2,2)+(
     #167.E0/198450.E0*RS(1,2)-997.E0/793800.E0*RS(1,1)-RS(1,4)/120)*RS(
     #2,3)**4
      t0 = s2+(RS(1,3)/120-149.E0/132300.E0*RS(1,1)-121.E0/264600.E0*RS(
     #1,2)-RS(1,4)/60)*RS(2,4)*RS(2,3)**3+(RS(1,3)/60-149.E0/132300.E0*R
     #S(1,2)-RS(1,4)/40-121.E0/264600.E0*RS(1,1))*RS(2,4)**2*RS(2,3)**2+
     #(-997.E0/793800.E0*RS(1,2)+RS(1,3)/40-RS(1,4)/30+167.E0/198450.E0*
     #RS(1,1))*RS(2,4)**3*RS(2,3)+(RS(1,3)/30-RS(1,1)/30)*RS(2,4)**4
	ASSS(4) = T0

      s2 = (124.E0/14175.E0*RS(1,4)-124.E0/14175.E0*RS(1,2))*RS(2,1)**4+
     #((52.E0/33075.E0*RS(1,4)-328.E0/33075.E0*RS(1,2)-8.E0/19845.E0*RS(
     #1,3)+124.E0/14175.E0*RS(1,1))*RS(2,2)+(22.E0/14175.E0*RS(1,2)-22.E
     #0/14175.E0*RS(1,4))*RS(2,3)+(-52.E0/33075.E0*RS(1,2)-124.E0/14175.
     #E0*RS(1,1)+328.E0/33075.E0*RS(1,4)+8.E0/19845.E0*RS(1,3))*RS(2,4))
     #*RS(2,1)**3+((328.E0/33075.E0*RS(1,1)+4.E0/4725.E0*RS(1,4)-328.E0/
     #33075.E0*RS(1,2)-4.E0/4725.E0*RS(1,3))*RS(2,2)**2+((-11.E0/33075.E
     #0*RS(1,4)+2.E0/1323.E0*RS(1,2)-RS(1,3)/33075-38.E0/33075.E0*RS(1,1
     #))*RS(2,3)+(-22.E0/33075.E0*RS(1,4)+22.E0/33075.E0*RS(1,2))*RS(2,4
     #))*RS(2,2)+(RS(1,3)/33075+11.E0/33075.E0*RS(1,2)-2.E0/1323.E0*RS(1
     #,4)+38.E0/33075.E0*RS(1,1))*RS(2,4)*RS(2,3)+(-4.E0/4725.E0*RS(1,2)
     #+328.E0/33075.E0*RS(1,4)+4.E0/4725.E0*RS(1,3)-328.E0/33075.E0*RS(1
     #,1))*RS(2,4)**2)*RS(2,1)**2
      s3 = s2+((-124.E0/14175.E0*RS(1,2)+8.E0/19845.E0*RS(1,4)+328.E0/33
     #075.E0*RS(1,1)-52.E0/33075.E0*RS(1,3))*RS(2,2)**3+((-22.E0/33075.E
     #0*RS(1,1)+22.E0/33075.E0*RS(1,3))*RS(2,3)+(-2.E0/1323.E0*RS(1,1)+R
     #S(1,4)/33075+11.E0/33075.E0*RS(1,3)+38.E0/33075.E0*RS(1,2))*RS(2,4
     #))*RS(2,2)**2+((38.E0/33075.E0*RS(1,3)+11.E0/33075.E0*RS(1,4)-2.E0
     #/1323.E0*RS(1,2)+RS(1,1)/33075)*RS(2,3)**2+(-RS(1,2)/33075-38.E0/3
     #3075.E0*RS(1,4)-11.E0/33075.E0*RS(1,3)+2.E0/1323.E0*RS(1,1))*RS(2,
     #4)**2)*RS(2,2)+(22.E0/14175.E0*RS(1,4)-22.E0/14175.E0*RS(1,2))*RS(
     #2,3)**3+(-38.E0/33075.E0*RS(1,3)-RS(1,1)/33075-11.E0/33075.E0*RS(1
     #,2)+2.E0/1323.E0*RS(1,4))*RS(2,4)*RS(2,3)**2+(-22.E0/33075.E0*RS(1
     #,3)+22.E0/33075.E0*RS(1,1))*RS(2,4)**2*RS(2,3)+(-328.E0/33075.E0*R
     #S(1,1)+52.E0/33075.E0*RS(1,3)+124.E0/14175.E0*RS(1,4)-8.E0/19845.E
     #0*RS(1,2))*RS(2,4)**3)*RS(2,1)
      s1 = s3+(-124.E0/14175.E0*RS(1,3)+124.E0/14175.E0*RS(1,1))*RS(2,2)
     #**4+((-328.E0/33075.E0*RS(1,3)-8.E0/19845.E0*RS(1,4)+124.E0/14175.
     #E0*RS(1,2)+52.E0/33075.E0*RS(1,1))*RS(2,3)+(-22.E0/14175.E0*RS(1,1
     #)+22.E0/14175.E0*RS(1,3))*RS(2,4))*RS(2,2)**3
      s2 = s1+((-4.E0/4725.E0*RS(1,4)+4.E0/4725.E0*RS(1,1)-328.E0/33075.
     #E0*RS(1,3)+328.E0/33075.E0*RS(1,2))*RS(2,3)**2+(2.E0/1323.E0*RS(1,
     #3)-11.E0/33075.E0*RS(1,1)-RS(1,4)/33075-38.E0/33075.E0*RS(1,2))*RS
     #(2,4)*RS(2,3))*RS(2,2)**2+((328.E0/33075.E0*RS(1,2)-124.E0/14175.E
     #0*RS(1,3)-52.E0/33075.E0*RS(1,4)+8.E0/19845.E0*RS(1,1))*RS(2,3)**3
     #+(-22.E0/33075.E0*RS(1,2)+22.E0/33075.E0*RS(1,4))*RS(2,4)*RS(2,3)*
     #*2+(11.E0/33075.E0*RS(1,1)+RS(1,2)/33075-2.E0/1323.E0*RS(1,3)+38.E
     #0/33075.E0*RS(1,4))*RS(2,4)**2*RS(2,3)+(-22.E0/14175.E0*RS(1,3)+22
     #.E0/14175.E0*RS(1,1))*RS(2,4)**3)*RS(2,2)+(-124.E0/14175.E0*RS(1,4
     #)+124.E0/14175.E0*RS(1,2))*RS(2,3)**4
      t0 = s2+(52.E0/33075.E0*RS(1,2)-8.E0/19845.E0*RS(1,1)-328.E0/33075
     #.E0*RS(1,4)+124.E0/14175.E0*RS(1,3))*RS(2,4)*RS(2,3)**3+(328.E0/33
     #075.E0*RS(1,3)+4.E0/4725.E0*RS(1,2)-328.E0/33075.E0*RS(1,4)-4.E0/4
     #725.E0*RS(1,1))*RS(2,4)**2*RS(2,3)**2+(328.E0/33075.E0*RS(1,3)-52.
     #E0/33075.E0*RS(1,1)+8.E0/19845.E0*RS(1,2)-124.E0/14175.E0*RS(1,4))
     #*RS(2,4)**3*RS(2,3)+(124.E0/14175.E0*RS(1,3)-124.E0/14175.E0*RS(1,
     #1))*RS(2,4)**4
	ASSS(5) = T0

C	-----------------------
C	COMPUTE FOR FACTOR ARRS
C	-----------------------
      s3 = (RS(2,2)*RS(2,1)/30-RS(2,1)*RS(2,4)/30-RS(2,4)**2/120+RS(2,2)
     #**2/120)*RS(1,1)**3
      s4 = ((167.E0/198450.E0*RS(2,1)*RS(2,4)-997.E0/793800.E0*RS(2,1)*R
     #S(2,3)-73.E0/0.15876E7*RS(2,3)*RS(2,4)-149.E0/396900.E0*RS(2,2)*RS
     #(2,3)-RS(2,1)**2/30-121.E0/793800.E0*RS(2,4)**2-121.E0/793800.E0*R
     #S(2,2)*RS(2,4)+RS(2,2)**2/90+RS(2,2)*RS(2,1)/120+17.E0/317520.E0*R
     #S(2,3)**2)*RS(1,2)+(997.E0/793800.E0*RS(2,2)*RS(2,1)+17.E0/317520.
     #E0*RS(2,3)*RS(2,4)+149.E0/396900.E0*RS(2,2)**2-997.E0/793800.E0*RS
     #(2,1)*RS(2,4)-17.E0/317520.E0*RS(2,2)*RS(2,3)-149.E0/396900.E0*RS(
     #2,4)**2)*RS(1,3)+(-RS(2,4)**2/90+121.E0/793800.E0*RS(2,2)*RS(2,4)-
     #167.E0/198450.E0*RS(2,2)*RS(2,1)-RS(2,1)*RS(2,4)/120+149.E0/396900
     #.E0*RS(2,3)*RS(2,4)+RS(2,1)**2/30+73.E0/0.15876E7*RS(2,2)*RS(2,3)-
     #17.E0/317520.E0*RS(2,3)**2+997.E0/793800.E0*RS(2,1)*RS(2,3)+121.E0
     #/793800.E0*RS(2,2)**2)*RS(1,4))*RS(1,1)**2
      s2 = s3+s4
      s4 = s2
      s7 = (-RS(2,2)*RS(2,1)/180-121.E0/396900.E0*RS(2,2)*RS(2,3)-RS(2,1
     #)**2/60-149.E0/198450.E0*RS(2,2)*RS(2,4)+RS(2,2)**2/120+131.E0/0.1
     #5876E7*RS(2,4)**2+13.E0/396900.E0*RS(2,3)*RS(2,4)-121.E0/396900.E0
     #*RS(2,1)*RS(2,4)+131.E0/0.15876E7*RS(2,3)**2-149.E0/198450.E0*RS(2
     #,1)*RS(2,3))*RS(1,2)**2+((-47.E0/793800.E0*RS(2,2)*RS(2,4)+149.E0/
     #198450.E0*RS(2,2)*RS(2,1)+13.E0/396900.E0*RS(2,4)**2+121.E0/396900
     #.E0*RS(2,2)**2+RS(2,3)*RS(2,4)/26460+17.E0/158760.E0*RS(2,1)*RS(2,
     #3)-73.E0/793800.E0*RS(2,1)*RS(2,4)+17.E0/158760.E0*RS(2,3)**2)*RS(
     #1,3)+(47.E0/793800.E0*RS(2,2)*RS(2,3)+121.E0/396900.E0*RS(2,2)*RS(
     #2,1)-121.E0/396900.E0*RS(2,1)*RS(2,4)-149.E0/198450.E0*RS(2,4)**2+
     #149.E0/198450.E0*RS(2,2)**2-47.E0/793800.E0*RS(2,3)*RS(2,4))*RS(1,
     #4))*RS(1,2)
      s6 = s7+(17.E0/158760.E0*RS(2,3)*RS(2,4)-131.E0/0.15876E7*RS(2,2)*
     #*2-17.E0/158760.E0*RS(2,2)*RS(2,3)+17.E0/158760.E0*RS(2,1)*RS(2,4)
     #+131.E0/0.15876E7*RS(2,4)**2-17.E0/158760.E0*RS(2,2)*RS(2,1))*RS(1
     #,3)**2+(-13.E0/396900.E0*RS(2,2)**2-17.E0/158760.E0*RS(2,1)*RS(2,3
     #)-RS(2,2)*RS(2,3)/26460-17.E0/158760.E0*RS(2,3)**2-149.E0/198450.E
     #0*RS(2,1)*RS(2,4)-121.E0/396900.E0*RS(2,4)**2+73.E0/793800.E0*RS(2
     #,2)*RS(2,1)+47.E0/793800.E0*RS(2,2)*RS(2,4))*RS(1,4)*RS(1,3)+(149.
     #E0/198450.E0*RS(2,1)*RS(2,3)-131.E0/0.15876E7*RS(2,2)**2+RS(2,1)*R
     #S(2,4)/180+121.E0/396900.E0*RS(2,3)*RS(2,4)+121.E0/396900.E0*RS(2,
     #2)*RS(2,1)+149.E0/198450.E0*RS(2,2)*RS(2,4)-13.E0/396900.E0*RS(2,2
     #)*RS(2,3)-131.E0/0.15876E7*RS(2,3)**2+RS(2,1)**2/60-RS(2,4)**2/120
     #)*RS(1,4)**2
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s4 = s3
      s6 = (17.E0/317520.E0*RS(2,4)**2-121.E0/793800.E0*RS(2,1)*RS(2,3)-
     #997.E0/793800.E0*RS(2,2)*RS(2,4)-RS(2,2)*RS(2,1)/120-149.E0/396900
     #.E0*RS(2,1)*RS(2,4)+167.E0/198450.E0*RS(2,2)*RS(2,3)-RS(2,1)**2/18
     #0-73.E0/0.15876E7*RS(2,3)*RS(2,4)-121.E0/793800.E0*RS(2,3)**2)*RS(
     #1,2)**3
      s7 = ((-149.E0/198450.E0*RS(2,3)**2-23.E0/317520.E0*RS(2,3)*RS(2,4
     #)+121.E0/793800.E0*RS(2,2)*RS(2,1)-RS(2,1)*RS(2,4)/75600+131.E0/0.
     #15876E7*RS(2,1)*RS(2,3)-167.E0/198450.E0*RS(2,2)**2-361.E0/0.15876
     #E7*RS(2,2)*RS(2,4)+RS(2,4)**2/52920-121.E0/793800.E0*RS(2,2)*RS(2,
     #3))*RS(1,3)+(149.E0/396900.E0*RS(2,2)*RS(2,1)+131.E0/0.15876E7*RS(
     #2,1)*RS(2,4)+73.E0/0.15876E7*RS(2,1)*RS(2,3)+17.E0/317520.E0*RS(2,
     #2)*RS(2,4)+17.E0/158760.E0*RS(2,4)**2+2.E0/11025.E0*RS(2,2)*RS(2,3
     #)+997.E0/793800.E0*RS(2,2)**2+73.E0/0.15876E7*RS(2,3)**2+RS(2,3)*R
     #S(2,4)/52920)*RS(1,4))*RS(1,2)**2
      s5 = s6+s7
      s1 = s4+s5
      s3 = s1
      s5 = ((17.E0/317520.E0*RS(2,1)*RS(2,3)-131.E0/0.15876E7*RS(2,2)*RS
     #(2,1)+RS(2,1)*RS(2,4)/52920-997.E0/793800.E0*RS(2,3)**2+121.E0/396
     #900.E0*RS(2,2)**2-31.E0/113400.E0*RS(2,3)*RS(2,4)+149.E0/396900.E0
     #*RS(2,2)*RS(2,3)-RS(2,4)**2/75600-167.E0/0.15876E7*RS(2,2)*RS(2,4)
     #)*RS(1,3)**2+(13.E0/396900.E0*RS(2,1)*RS(2,4)-13.E0/396900.E0*RS(2
     #,2)*RS(2,1)-47.E0/396900.E0*RS(2,3)*RS(2,4)-73.E0/793800.E0*RS(2,4
     #)**2+73.E0/793800.E0*RS(2,2)**2+47.E0/396900.E0*RS(2,2)*RS(2,3))*R
     #S(1,4)*RS(1,3)+(-RS(2,2)*RS(2,3)/52920-73.E0/0.15876E7*RS(2,1)*RS(
     #2,3)-73.E0/0.15876E7*RS(2,3)**2-997.E0/793800.E0*RS(2,4)**2-2.E0/1
     #1025.E0*RS(2,3)*RS(2,4)-131.E0/0.15876E7*RS(2,2)*RS(2,1)-149.E0/39
     #6900.E0*RS(2,1)*RS(2,4)-17.E0/158760.E0*RS(2,2)**2-17.E0/317520.E0
     #*RS(2,2)*RS(2,4))*RS(1,4)**2)*RS(1,2)
      s6 = (997.E0/793800.E0*RS(2,2)*RS(2,3)-997.E0/793800.E0*RS(2,3)*RS
     #(2,4)-149.E0/396900.E0*RS(2,4)**2+17.E0/317520.E0*RS(2,1)*RS(2,4)+
     #149.E0/396900.E0*RS(2,2)**2-17.E0/317520.E0*RS(2,2)*RS(2,1))*RS(1,
     #3)**3
      s4 = s5+s6
      s2 = s3+s4
      s3 = s2+(RS(2,2)**2/75600-17.E0/317520.E0*RS(2,1)*RS(2,3)-RS(2,2)*
     #RS(2,1)/52920+31.E0/113400.E0*RS(2,2)*RS(2,3)-149.E0/396900.E0*RS(
     #2,3)*RS(2,4)+997.E0/793800.E0*RS(2,3)**2+167.E0/0.15876E7*RS(2,2)*
     #RS(2,4)-121.E0/396900.E0*RS(2,4)**2+131.E0/0.15876E7*RS(2,1)*RS(2,
     #4))*RS(1,4)*RS(1,3)**2
      t0 = s3+(-131.E0/0.15876E7*RS(2,1)*RS(2,3)-121.E0/793800.E0*RS(2,1
     #)*RS(2,4)+23.E0/317520.E0*RS(2,2)*RS(2,3)+361.E0/0.15876E7*RS(2,2)
     #*RS(2,4)+RS(2,2)*RS(2,1)/75600-RS(2,2)**2/52920+121.E0/793800.E0*R
     #S(2,3)*RS(2,4)+149.E0/198450.E0*RS(2,3)**2+167.E0/198450.E0*RS(2,4
     #)**2)*RS(1,4)**2*RS(1,3)+(121.E0/793800.E0*RS(2,3)**2+RS(2,1)**2/1
     #80-167.E0/198450.E0*RS(2,3)*RS(2,4)-17.E0/317520.E0*RS(2,2)**2+RS(
     #2,1)*RS(2,4)/120+121.E0/793800.E0*RS(2,1)*RS(2,3)+73.E0/0.15876E7*
     #RS(2,2)*RS(2,3)+149.E0/396900.E0*RS(2,2)*RS(2,1)+997.E0/793800.E0*
     #RS(2,2)*RS(2,4))*RS(1,4)**3
	ARRS(1) = T0

      s3 = (121.E0/793800.E0*RS(2,4)*RS(2,2)+73.E0/0.15876E7*RS(2,4)*RS(
     #2,3)+RS(2,2)**2/180+997.E0/793800.E0*RS(2,1)*RS(2,3)-167.E0/198450
     #.E0*RS(2,4)*RS(2,1)+121.E0/793800.E0*RS(2,4)**2-17.E0/317520.E0*RS
     #(2,3)**2+RS(2,1)*RS(2,2)/120+149.E0/396900.E0*RS(2,3)*RS(2,2))*RS(
     #1,1)**3
      s6 = (121.E0/396900.E0*RS(2,4)*RS(2,1)+RS(2,1)*RS(2,2)/180-131.E0/
     #0.15876E7*RS(2,3)**2+121.E0/396900.E0*RS(2,3)*RS(2,2)+RS(2,2)**2/6
     #0-13.E0/396900.E0*RS(2,4)*RS(2,3)+149.E0/198450.E0*RS(2,4)*RS(2,2)
     #-RS(2,1)**2/120-131.E0/0.15876E7*RS(2,4)**2+149.E0/198450.E0*RS(2,
     #1)*RS(2,3))*RS(1,2)
      s7 = (-17.E0/317520.E0*RS(2,1)*RS(2,3)-17.E0/158760.E0*RS(2,3)**2-
     #73.E0/0.15876E7*RS(2,4)*RS(2,2)-131.E0/0.15876E7*RS(2,3)*RS(2,2)-R
     #S(2,4)*RS(2,3)/52920-2.E0/11025.E0*RS(2,4)*RS(2,1)-73.E0/0.15876E7
     #*RS(2,4)**2-149.E0/396900.E0*RS(2,1)*RS(2,2)-997.E0/793800.E0*RS(2
     #,1)**2)*RS(1,3)+(361.E0/0.15876E7*RS(2,1)*RS(2,3)-RS(2,3)**2/52920
     #-131.E0/0.15876E7*RS(2,4)*RS(2,2)+23.E0/317520.E0*RS(2,4)*RS(2,3)+
     #RS(2,3)*RS(2,2)/75600-121.E0/793800.E0*RS(2,1)*RS(2,2)+149.E0/1984
     #50.E0*RS(2,4)**2+167.E0/198450.E0*RS(2,1)**2+121.E0/793800.E0*RS(2
     #,4)*RS(2,1))*RS(1,4)
      s5 = s6+s7
      s6 = RS(1,1)**2
      s4 = s5*s6
      s2 = s3+s4
      s4 = s2
      s7 = (121.E0/793800.E0*RS(2,3)**2-167.E0/198450.E0*RS(2,3)*RS(2,2)
     #-17.E0/317520.E0*RS(2,4)**2-RS(2,1)*RS(2,2)/120+121.E0/793800.E0*R
     #S(2,1)*RS(2,3)-RS(2,1)**2/90+149.E0/396900.E0*RS(2,4)*RS(2,1)+997.
     #E0/793800.E0*RS(2,4)*RS(2,2)+73.E0/0.15876E7*RS(2,4)*RS(2,3)+RS(2,
     #2)**2/30)*RS(1,2)**2+((-121.E0/396900.E0*RS(2,1)*RS(2,2)+47.E0/793
     #800.E0*RS(2,4)*RS(2,3)+149.E0/198450.E0*RS(2,3)**2-149.E0/198450.E
     #0*RS(2,1)**2+121.E0/396900.E0*RS(2,3)*RS(2,2)-47.E0/793800.E0*RS(2
     #,4)*RS(2,1))*RS(1,3)+(-13.E0/396900.E0*RS(2,3)**2+47.E0/793800.E0*
     #RS(2,1)*RS(2,3)-121.E0/396900.E0*RS(2,1)**2-17.E0/158760.E0*RS(2,4
     #)**2-149.E0/198450.E0*RS(2,1)*RS(2,2)-RS(2,4)*RS(2,3)/26460-17.E0/
     #158760.E0*RS(2,4)*RS(2,2)+73.E0/793800.E0*RS(2,3)*RS(2,2))*RS(1,4)
     #)*RS(1,2)
      s6 = s7+(149.E0/396900.E0*RS(2,3)*RS(2,2)+997.E0/793800.E0*RS(2,3)
     #**2+73.E0/0.15876E7*RS(2,4)*RS(2,2)+17.E0/317520.E0*RS(2,1)*RS(2,3
     #)+73.E0/0.15876E7*RS(2,4)**2+17.E0/158760.E0*RS(2,1)**2+131.E0/0.1
     #5876E7*RS(2,1)*RS(2,2)+RS(2,4)*RS(2,1)/52920+2.E0/11025.E0*RS(2,4)
     #*RS(2,3))*RS(1,3)**2+(47.E0/396900.E0*RS(2,4)*RS(2,3)-13.E0/396900
     #.E0*RS(2,3)*RS(2,2)-47.E0/396900.E0*RS(2,4)*RS(2,1)+73.E0/793800.E
     #0*RS(2,3)**2-73.E0/793800.E0*RS(2,1)**2+13.E0/396900.E0*RS(2,1)*RS
     #(2,2))*RS(1,4)*RS(1,3)+(-RS(2,3)*RS(2,2)/52920+31.E0/113400.E0*RS(
     #2,4)*RS(2,3)-149.E0/396900.E0*RS(2,4)*RS(2,1)+997.E0/793800.E0*RS(
     #2,4)**2+167.E0/0.15876E7*RS(2,1)*RS(2,3)-121.E0/396900.E0*RS(2,1)*
     #*2-17.E0/317520.E0*RS(2,4)*RS(2,2)+131.E0/0.15876E7*RS(2,1)*RS(2,2
     #)+RS(2,3)**2/75600)*RS(1,4)**2
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(-RS(2,1)**2/120+RS(2,3)**2/120+RS(2,3)*RS(2,2)/30-RS(2,1)
     #*RS(2,2)/30)*RS(1,2)**3+((-121.E0/793800.E0*RS(2,1)**2-121.E0/7938
     #00.E0*RS(2,1)*RS(2,3)+RS(2,3)*RS(2,2)/120-73.E0/0.15876E7*RS(2,4)*
     #RS(2,1)-997.E0/793800.E0*RS(2,4)*RS(2,2)-RS(2,2)**2/30+RS(2,3)**2/
     #90-149.E0/396900.E0*RS(2,4)*RS(2,3)+167.E0/198450.E0*RS(2,1)*RS(2,
     #2)+17.E0/317520.E0*RS(2,4)**2)*RS(1,3)+(-149.E0/396900.E0*RS(2,1)*
     #*2+17.E0/317520.E0*RS(2,4)*RS(2,1)+997.E0/793800.E0*RS(2,3)*RS(2,2
     #)-17.E0/317520.E0*RS(2,4)*RS(2,3)-997.E0/793800.E0*RS(2,1)*RS(2,2)
     #+149.E0/396900.E0*RS(2,3)**2)*RS(1,4))*RS(1,2)**2
      s3 = s1
      s5 = ((-RS(2,3)*RS(2,2)/180+RS(2,3)**2/120-121.E0/396900.E0*RS(2,1
     #)*RS(2,2)-149.E0/198450.E0*RS(2,4)*RS(2,2)-121.E0/396900.E0*RS(2,4
     #)*RS(2,3)-RS(2,2)**2/60+13.E0/396900.E0*RS(2,4)*RS(2,1)-149.E0/198
     #450.E0*RS(2,1)*RS(2,3)+131.E0/0.15876E7*RS(2,4)**2+131.E0/0.15876E
     #7*RS(2,1)**2)*RS(1,3)**2+(17.E0/158760.E0*RS(2,4)**2+149.E0/198450
     #.E0*RS(2,3)*RS(2,2)+17.E0/158760.E0*RS(2,4)*RS(2,2)-73.E0/793800.E
     #0*RS(2,1)*RS(2,2)+13.E0/396900.E0*RS(2,1)**2-47.E0/793800.E0*RS(2,
     #1)*RS(2,3)+RS(2,4)*RS(2,1)/26460+121.E0/396900.E0*RS(2,3)**2)*RS(1
     #,4)*RS(1,3)+(17.E0/158760.E0*RS(2,4)*RS(2,1)-17.E0/158760.E0*RS(2,
     #4)*RS(2,3)-17.E0/158760.E0*RS(2,3)*RS(2,2)+17.E0/158760.E0*RS(2,1)
     #*RS(2,2)+131.E0/0.15876E7*RS(2,1)**2-131.E0/0.15876E7*RS(2,3)**2)*
     #RS(1,4)**2)*RS(1,2)
      s6 = (17.E0/317520.E0*RS(2,1)**2-149.E0/396900.E0*RS(2,1)*RS(2,2)-
     #73.E0/0.15876E7*RS(2,4)*RS(2,1)-RS(2,3)*RS(2,2)/120+167.E0/198450.
     #E0*RS(2,4)*RS(2,3)-121.E0/793800.E0*RS(2,4)*RS(2,2)-RS(2,2)**2/180
     #-997.E0/793800.E0*RS(2,1)*RS(2,3)-121.E0/793800.E0*RS(2,4)**2)*RS(
     #1,3)**3
      s4 = s5+s6
      s2 = s3+s4
      t0 = s2+(-RS(2,1)*RS(2,2)/75600+131.E0/0.15876E7*RS(2,4)*RS(2,2)+R
     #S(2,1)**2/52920-361.E0/0.15876E7*RS(2,1)*RS(2,3)+121.E0/793800.E0*
     #RS(2,3)*RS(2,2)-23.E0/317520.E0*RS(2,4)*RS(2,1)-149.E0/198450.E0*R
     #S(2,4)**2-167.E0/198450.E0*RS(2,3)**2-121.E0/793800.E0*RS(2,4)*RS(
     #2,3))*RS(1,4)*RS(1,3)**2+(RS(2,1)*RS(2,2)/52920-167.E0/0.15876E7*R
     #S(2,1)*RS(2,3)-RS(2,1)**2/75600+17.E0/317520.E0*RS(2,4)*RS(2,2)-99
     #7.E0/793800.E0*RS(2,4)**2-31.E0/113400.E0*RS(2,4)*RS(2,1)+149.E0/3
     #96900.E0*RS(2,4)*RS(2,3)+121.E0/396900.E0*RS(2,3)**2-131.E0/0.1587
     #6E7*RS(2,3)*RS(2,2))*RS(1,4)**2*RS(1,3)+(149.E0/396900.E0*RS(2,3)*
     #*2-997.E0/793800.E0*RS(2,4)*RS(2,1)-17.E0/317520.E0*RS(2,3)*RS(2,2
     #)-149.E0/396900.E0*RS(2,1)**2+997.E0/793800.E0*RS(2,4)*RS(2,3)+17.
     #E0/317520.E0*RS(2,1)*RS(2,2))*RS(1,4)**3
	ARRS(2) = T0

      s3 = (17.E0/317520.E0*RS(2,3)*RS(2,2)+997.E0/793800.E0*RS(2,4)*RS(
     #2,1)+149.E0/396900.E0*RS(2,4)**2-997.E0/793800.E0*RS(2,1)*RS(2,2)-
     #17.E0/317520.E0*RS(2,4)*RS(2,3)-149.E0/396900.E0*RS(2,2)**2)*RS(1,
     #1)**3
      s4 = ((131.E0/0.15876E7*RS(2,3)*RS(2,2)+997.E0/793800.E0*RS(2,1)**
     #2+RS(2,4)**2/75600-RS(2,4)*RS(2,3)/52920+167.E0/0.15876E7*RS(2,4)*
     #RS(2,2)-17.E0/317520.E0*RS(2,3)*RS(2,1)-121.E0/396900.E0*RS(2,2)**
     #2-149.E0/396900.E0*RS(2,1)*RS(2,2)+31.E0/113400.E0*RS(2,4)*RS(2,1)
     #)*RS(1,2)+(17.E0/158760.E0*RS(2,1)*RS(2,2)+17.E0/158760.E0*RS(2,3)
     #*RS(2,2)-131.E0/0.15876E7*RS(2,4)**2-17.E0/158760.E0*RS(2,4)*RS(2,
     #1)-17.E0/158760.E0*RS(2,4)*RS(2,3)+131.E0/0.15876E7*RS(2,2)**2)*RS
     #(1,3)+(149.E0/396900.E0*RS(2,4)*RS(2,1)+RS(2,3)*RS(2,2)/52920+17.E
     #0/317520.E0*RS(2,3)*RS(2,1)+121.E0/396900.E0*RS(2,4)**2-31.E0/1134
     #00.E0*RS(2,1)*RS(2,2)-131.E0/0.15876E7*RS(2,4)*RS(2,3)-RS(2,2)**2/
     #75600-167.E0/0.15876E7*RS(2,4)*RS(2,2)-997.E0/793800.E0*RS(2,1)**2
     #)*RS(1,4))*RS(1,1)**2
      s2 = s3+s4
      s4 = s2
      s7 = (RS(2,4)*RS(2,3)/75600+361.E0/0.15876E7*RS(2,4)*RS(2,2)-121.E
     #0/793800.E0*RS(2,3)*RS(2,2)+149.E0/198450.E0*RS(2,1)**2+167.E0/198
     #450.E0*RS(2,2)**2+121.E0/793800.E0*RS(2,1)*RS(2,2)-131.E0/0.15876E
     #7*RS(2,3)*RS(2,1)-RS(2,4)**2/52920+23.E0/317520.E0*RS(2,4)*RS(2,1)
     #)*RS(1,2)**2+((-RS(2,4)*RS(2,1)/26460-13.E0/396900.E0*RS(2,4)**2-1
     #7.E0/158760.E0*RS(2,1)**2-121.E0/396900.E0*RS(2,2)**2+47.E0/793800
     #.E0*RS(2,4)*RS(2,2)-149.E0/198450.E0*RS(2,3)*RS(2,2)+73.E0/793800.
     #E0*RS(2,4)*RS(2,3)-17.E0/158760.E0*RS(2,3)*RS(2,1))*RS(1,3)+(-13.E
     #0/396900.E0*RS(2,4)*RS(2,3)-73.E0/793800.E0*RS(2,2)**2+73.E0/79380
     #0.E0*RS(2,4)**2+47.E0/396900.E0*RS(2,4)*RS(2,1)+13.E0/396900.E0*RS
     #(2,3)*RS(2,2)-47.E0/396900.E0*RS(2,1)*RS(2,2))*RS(1,4))*RS(1,2)
      s6 = s7+(149.E0/396900.E0*RS(2,4)**2+997.E0/793800.E0*RS(2,4)*RS(2
     #,3)-149.E0/396900.E0*RS(2,2)**2-997.E0/793800.E0*RS(2,3)*RS(2,2)+1
     #7.E0/317520.E0*RS(2,1)*RS(2,2)-17.E0/317520.E0*RS(2,4)*RS(2,1))*RS
     #(1,3)**2+(149.E0/198450.E0*RS(2,4)*RS(2,3)-73.E0/793800.E0*RS(2,3)
     #*RS(2,2)+121.E0/396900.E0*RS(2,4)**2+17.E0/158760.E0*RS(2,1)**2+13
     #.E0/396900.E0*RS(2,2)**2-47.E0/793800.E0*RS(2,4)*RS(2,2)+RS(2,1)*R
     #S(2,2)/26460+17.E0/158760.E0*RS(2,3)*RS(2,1))*RS(1,4)*RS(1,3)+(131
     #.E0/0.15876E7*RS(2,3)*RS(2,1)-361.E0/0.15876E7*RS(2,4)*RS(2,2)-RS(
     #2,3)*RS(2,2)/75600-149.E0/198450.E0*RS(2,1)**2+RS(2,2)**2/52920-23
     #.E0/317520.E0*RS(2,1)*RS(2,2)-121.E0/793800.E0*RS(2,4)*RS(2,1)-167
     #.E0/198450.E0*RS(2,4)**2+121.E0/793800.E0*RS(2,4)*RS(2,3))*RS(1,4)
     #**2
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s4 = s3
      s6 = (-17.E0/317520.E0*RS(2,4)**2+121.E0/793800.E0*RS(2,3)*RS(2,1)
     #+121.E0/793800.E0*RS(2,1)**2+997.E0/793800.E0*RS(2,4)*RS(2,2)+73.E
     #0/0.15876E7*RS(2,4)*RS(2,1)+149.E0/396900.E0*RS(2,4)*RS(2,3)+RS(2,
     #3)**2/180-167.E0/198450.E0*RS(2,1)*RS(2,2)+RS(2,3)*RS(2,2)/120)*RS
     #(1,2)**3
      s7 = ((149.E0/198450.E0*RS(2,4)*RS(2,2)+149.E0/198450.E0*RS(2,3)*R
     #S(2,1)-13.E0/396900.E0*RS(2,4)*RS(2,1)-131.E0/0.15876E7*RS(2,4)**2
     #+RS(2,3)*RS(2,2)/180+121.E0/396900.E0*RS(2,4)*RS(2,3)-RS(2,2)**2/1
     #20-131.E0/0.15876E7*RS(2,1)**2+RS(2,3)**2/60+121.E0/396900.E0*RS(2
     #,1)*RS(2,2))*RS(1,3)+(-149.E0/396900.E0*RS(2,3)*RS(2,2)-17.E0/1587
     #60.E0*RS(2,4)**2-131.E0/0.15876E7*RS(2,4)*RS(2,3)-2.E0/11025.E0*RS
     #(2,1)*RS(2,2)-73.E0/0.15876E7*RS(2,3)*RS(2,1)-73.E0/0.15876E7*RS(2
     #,1)**2-997.E0/793800.E0*RS(2,2)**2-RS(2,4)*RS(2,1)/52920-17.E0/317
     #520.E0*RS(2,4)*RS(2,2))*RS(1,4))*RS(1,2)**2
      s5 = s6+s7
      s1 = s4+s5
      s3 = s1
      s5 = ((-RS(2,3)*RS(2,2)/120+997.E0/793800.E0*RS(2,3)*RS(2,1)+121.E
     #0/793800.E0*RS(2,4)*RS(2,2)-167.E0/198450.E0*RS(2,4)*RS(2,3)-17.E0
     #/317520.E0*RS(2,1)**2-RS(2,2)**2/90+121.E0/793800.E0*RS(2,4)**2+RS
     #(2,3)**2/30+149.E0/396900.E0*RS(2,1)*RS(2,2)+73.E0/0.15876E7*RS(2,
     #4)*RS(2,1))*RS(1,3)**2+(149.E0/198450.E0*RS(2,4)**2-149.E0/198450.
     #E0*RS(2,2)**2-121.E0/396900.E0*RS(2,3)*RS(2,2)+47.E0/793800.E0*RS(
     #2,4)*RS(2,1)-47.E0/793800.E0*RS(2,1)*RS(2,2)+121.E0/396900.E0*RS(2
     #,4)*RS(2,3))*RS(1,4)*RS(1,3)+(149.E0/396900.E0*RS(2,4)*RS(2,3)+131
     #.E0/0.15876E7*RS(2,3)*RS(2,2)+17.E0/158760.E0*RS(2,2)**2+RS(2,1)*R
     #S(2,2)/52920+997.E0/793800.E0*RS(2,4)**2+73.E0/0.15876E7*RS(2,1)**
     #2+2.E0/11025.E0*RS(2,4)*RS(2,1)+17.E0/317520.E0*RS(2,4)*RS(2,2)+73
     #.E0/0.15876E7*RS(2,3)*RS(2,1))*RS(1,4)**2)*RS(1,2)
      s6 = (-RS(2,3)*RS(2,2)/30+RS(2,4)**2/120-RS(2,2)**2/120+RS(2,4)*RS
     #(2,3)/30)*RS(1,3)**3
      s4 = s5+s6
      s2 = s3+s4
      s3 = s2+(-121.E0/793800.E0*RS(2,4)*RS(2,2)-149.E0/396900.E0*RS(2,4
     #)*RS(2,1)+RS(2,4)**2/90-RS(2,3)**2/30+17.E0/317520.E0*RS(2,1)**2+R
     #S(2,4)*RS(2,3)/120+167.E0/198450.E0*RS(2,3)*RS(2,2)-121.E0/793800.
     #E0*RS(2,2)**2-997.E0/793800.E0*RS(2,3)*RS(2,1)-73.E0/0.15876E7*RS(
     #2,1)*RS(2,2))*RS(1,4)*RS(1,3)**2
      t0 = s3+(RS(2,4)**2/120-149.E0/198450.E0*RS(2,3)*RS(2,1)-RS(2,4)*R
     #S(2,3)/180-121.E0/396900.E0*RS(2,4)*RS(2,1)-149.E0/198450.E0*RS(2,
     #4)*RS(2,2)+13.E0/396900.E0*RS(2,1)*RS(2,2)-121.E0/396900.E0*RS(2,3
     #)*RS(2,2)+131.E0/0.15876E7*RS(2,2)**2-RS(2,3)**2/60+131.E0/0.15876
     #E7*RS(2,1)**2)*RS(1,4)**2*RS(1,3)+(167.E0/198450.E0*RS(2,4)*RS(2,1
     #)-997.E0/793800.E0*RS(2,4)*RS(2,2)-121.E0/793800.E0*RS(2,1)**2-RS(
     #2,3)**2/180-121.E0/793800.E0*RS(2,3)*RS(2,1)-73.E0/0.15876E7*RS(2,
     #1)*RS(2,2)+17.E0/317520.E0*RS(2,2)**2-149.E0/396900.E0*RS(2,3)*RS(
     #2,2)-RS(2,4)*RS(2,3)/120)*RS(1,4)**3
	ARRS(3) = T0

      s3 = (17.E0/317520.E0*RS(2,3)**2-997.E0/793800.E0*RS(2,3)*RS(2,1)-
     #121.E0/793800.E0*RS(2,2)*RS(2,4)+167.E0/198450.E0*RS(2,2)*RS(2,1)-
     #RS(2,1)*RS(2,4)/120-121.E0/793800.E0*RS(2,2)**2-149.E0/396900.E0*R
     #S(2,3)*RS(2,4)-RS(2,4)**2/180-73.E0/0.15876E7*RS(2,2)*RS(2,3))*RS(
     #1,1)**3
      s6 = (-RS(2,3)*RS(2,4)/75600-121.E0/793800.E0*RS(2,2)*RS(2,1)+RS(2
     #,3)**2/52920-361.E0/0.15876E7*RS(2,3)*RS(2,1)-23.E0/317520.E0*RS(2
     #,2)*RS(2,3)-167.E0/198450.E0*RS(2,1)**2+131.E0/0.15876E7*RS(2,2)*R
     #S(2,4)-149.E0/198450.E0*RS(2,2)**2+121.E0/793800.E0*RS(2,1)*RS(2,4
     #))*RS(1,2)
      s7 = (131.E0/0.15876E7*RS(2,3)*RS(2,4)+17.E0/158760.E0*RS(2,3)**2+
     #RS(2,2)*RS(2,3)/52920+2.E0/11025.E0*RS(2,2)*RS(2,1)+17.E0/317520.E
     #0*RS(2,3)*RS(2,1)+73.E0/0.15876E7*RS(2,2)*RS(2,4)+73.E0/0.15876E7*
     #RS(2,2)**2+149.E0/396900.E0*RS(2,1)*RS(2,4)+997.E0/793800.E0*RS(2,
     #1)**2)*RS(1,3)+(-121.E0/396900.E0*RS(2,2)*RS(2,1)+131.E0/0.15876E7
     #*RS(2,2)**2-149.E0/198450.E0*RS(2,3)*RS(2,1)-RS(2,1)*RS(2,4)/180-1
     #21.E0/396900.E0*RS(2,3)*RS(2,4)-149.E0/198450.E0*RS(2,2)*RS(2,4)+1
     #3.E0/396900.E0*RS(2,2)*RS(2,3)-RS(2,4)**2/60+131.E0/0.15876E7*RS(2
     #,3)**2+RS(2,1)**2/120)*RS(1,4)
      s5 = s6+s7
      s6 = RS(1,1)**2
      s4 = s5*s6
      s2 = s3+s4
      s4 = s2
      s7 = (-31.E0/113400.E0*RS(2,2)*RS(2,3)-RS(2,3)**2/75600+121.E0/396
     #900.E0*RS(2,1)**2+17.E0/317520.E0*RS(2,2)*RS(2,4)-997.E0/793800.E0
     #*RS(2,2)**2+149.E0/396900.E0*RS(2,2)*RS(2,1)-167.E0/0.15876E7*RS(2
     #,3)*RS(2,1)+RS(2,3)*RS(2,4)/52920-131.E0/0.15876E7*RS(2,1)*RS(2,4)
     #)*RS(1,2)**2+((-73.E0/793800.E0*RS(2,3)**2+73.E0/793800.E0*RS(2,1)
     #**2-47.E0/396900.E0*RS(2,2)*RS(2,3)+47.E0/396900.E0*RS(2,2)*RS(2,1
     #)-13.E0/396900.E0*RS(2,1)*RS(2,4)+13.E0/396900.E0*RS(2,3)*RS(2,4))
     #*RS(1,3)+(-47.E0/793800.E0*RS(2,3)*RS(2,1)+17.E0/158760.E0*RS(2,2)
     #*RS(2,4)+13.E0/396900.E0*RS(2,3)**2-73.E0/793800.E0*RS(2,3)*RS(2,4
     #)+149.E0/198450.E0*RS(2,1)*RS(2,4)+121.E0/396900.E0*RS(2,1)**2+17.
     #E0/158760.E0*RS(2,2)**2+RS(2,2)*RS(2,3)/26460)*RS(1,4))*RS(1,2)
      s6 = s7+(-17.E0/158760.E0*RS(2,1)**2-73.E0/0.15876E7*RS(2,2)*RS(2,
     #4)-73.E0/0.15876E7*RS(2,2)**2-997.E0/793800.E0*RS(2,3)**2-17.E0/31
     #7520.E0*RS(2,3)*RS(2,1)-131.E0/0.15876E7*RS(2,1)*RS(2,4)-149.E0/39
     #6900.E0*RS(2,3)*RS(2,4)-RS(2,2)*RS(2,1)/52920-2.E0/11025.E0*RS(2,2
     #)*RS(2,3))*RS(1,3)**2+(-149.E0/198450.E0*RS(2,3)**2+121.E0/396900.
     #E0*RS(2,1)*RS(2,4)-47.E0/793800.E0*RS(2,2)*RS(2,3)+149.E0/198450.E
     #0*RS(2,1)**2-121.E0/396900.E0*RS(2,3)*RS(2,4)+47.E0/793800.E0*RS(2
     #,2)*RS(2,1))*RS(1,4)*RS(1,3)+(RS(2,1)*RS(2,4)/120-73.E0/0.15876E7*
     #RS(2,2)*RS(2,3)+167.E0/198450.E0*RS(2,3)*RS(2,4)-149.E0/396900.E0*
     #RS(2,2)*RS(2,1)-121.E0/793800.E0*RS(2,3)*RS(2,1)+17.E0/317520.E0*R
     #S(2,2)**2-121.E0/793800.E0*RS(2,3)**2-RS(2,4)**2/30-997.E0/793800.
     #E0*RS(2,2)*RS(2,4)+RS(2,1)**2/90)*RS(1,4)**2
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(17.E0/317520.E0*RS(2,3)*RS(2,4)-17.E0/317520.E0*RS(2,1)*R
     #S(2,4)-149.E0/396900.E0*RS(2,3)**2-997.E0/793800.E0*RS(2,2)*RS(2,3
     #)+997.E0/793800.E0*RS(2,2)*RS(2,1)+149.E0/396900.E0*RS(2,1)**2)*RS
     #(1,2)**3+((-17.E0/317520.E0*RS(2,2)*RS(2,4)+31.E0/113400.E0*RS(2,2
     #)*RS(2,1)+167.E0/0.15876E7*RS(2,3)*RS(2,1)+131.E0/0.15876E7*RS(2,3
     #)*RS(2,4)+997.E0/793800.E0*RS(2,2)**2+RS(2,1)**2/75600-RS(2,1)*RS(
     #2,4)/52920-149.E0/396900.E0*RS(2,2)*RS(2,3)-121.E0/396900.E0*RS(2,
     #3)**2)*RS(1,3)+(-17.E0/158760.E0*RS(2,1)*RS(2,4)+131.E0/0.15876E7*
     #RS(2,3)**2+17.E0/158760.E0*RS(2,2)*RS(2,3)-131.E0/0.15876E7*RS(2,1
     #)**2+17.E0/158760.E0*RS(2,3)*RS(2,4)-17.E0/158760.E0*RS(2,2)*RS(2,
     #1))*RS(1,4))*RS(1,2)**2
      s3 = s1
      s5 = ((121.E0/793800.E0*RS(2,2)*RS(2,3)+361.E0/0.15876E7*RS(2,3)*R
     #S(2,1)+RS(2,1)*RS(2,4)/75600+167.E0/198450.E0*RS(2,3)**2+23.E0/317
     #520.E0*RS(2,2)*RS(2,1)-131.E0/0.15876E7*RS(2,2)*RS(2,4)-121.E0/793
     #800.E0*RS(2,3)*RS(2,4)-RS(2,1)**2/52920+149.E0/198450.E0*RS(2,2)**
     #2)*RS(1,3)**2+(-149.E0/198450.E0*RS(2,3)*RS(2,4)-RS(2,2)*RS(2,1)/2
     #6460-121.E0/396900.E0*RS(2,3)**2+47.E0/793800.E0*RS(2,3)*RS(2,1)+7
     #3.E0/793800.E0*RS(2,1)*RS(2,4)-13.E0/396900.E0*RS(2,1)**2-17.E0/15
     #8760.E0*RS(2,2)**2-17.E0/158760.E0*RS(2,2)*RS(2,4))*RS(1,4)*RS(1,3
     #)+(149.E0/396900.E0*RS(2,1)**2-149.E0/396900.E0*RS(2,3)**2-997.E0/
     #793800.E0*RS(2,3)*RS(2,4)+17.E0/317520.E0*RS(2,2)*RS(2,3)+997.E0/7
     #93800.E0*RS(2,1)*RS(2,4)-17.E0/317520.E0*RS(2,2)*RS(2,1))*RS(1,4)*
     #*2)*RS(1,2)
      s6 = (73.E0/0.15876E7*RS(2,2)*RS(2,1)+121.E0/793800.E0*RS(2,2)*RS(
     #2,4)-17.E0/317520.E0*RS(2,1)**2+997.E0/793800.E0*RS(2,3)*RS(2,1)+1
     #21.E0/793800.E0*RS(2,2)**2+RS(2,3)*RS(2,4)/120+149.E0/396900.E0*RS
     #(2,1)*RS(2,4)+RS(2,4)**2/180-167.E0/198450.E0*RS(2,2)*RS(2,3))*RS(
     #1,3)**3
      s4 = s5+s6
      s2 = s3+s4
      t0 = s2+(-13.E0/396900.E0*RS(2,2)*RS(2,1)-131.E0/0.15876E7*RS(2,2)
     #**2+149.E0/198450.E0*RS(2,3)*RS(2,1)-131.E0/0.15876E7*RS(2,1)**2+1
     #21.E0/396900.E0*RS(2,1)*RS(2,4)+149.E0/198450.E0*RS(2,2)*RS(2,4)+R
     #S(2,3)*RS(2,4)/180+121.E0/396900.E0*RS(2,2)*RS(2,3)-RS(2,3)**2/120
     #+RS(2,4)**2/60)*RS(1,4)*RS(1,3)**2+(-RS(2,3)*RS(2,4)/120-167.E0/19
     #8450.E0*RS(2,1)*RS(2,4)+121.E0/793800.E0*RS(2,3)*RS(2,1)+149.E0/39
     #6900.E0*RS(2,2)*RS(2,3)+RS(2,4)**2/30+73.E0/0.15876E7*RS(2,2)*RS(2
     #,1)+997.E0/793800.E0*RS(2,2)*RS(2,4)+121.E0/793800.E0*RS(2,1)**2-R
     #S(2,3)**2/90-17.E0/317520.E0*RS(2,2)**2)*RS(1,4)**2*RS(1,3)+(-RS(2
     #,3)**2/120+RS(2,1)*RS(2,4)/30-RS(2,3)*RS(2,4)/30+RS(2,1)**2/120)*R
     #S(1,4)**3
	ARRS(4) = T0

      s3 = (-328.E0/99225.E0*RS(2,4)**2+328.E0/99225.E0*RS(2,2)**2+124.E
     #0/14175.E0*RS(2,1)*RS(2,2)-38.E0/99225.E0*RS(2,3)*RS(2,2)+38.E0/99
     #225.E0*RS(2,3)*RS(2,4)-124.E0/14175.E0*RS(2,1)*RS(2,4))*RS(1,1)**3
      s4 = ((328.E0/99225.E0*RS(2,1)*RS(2,2)+26.E0/33075.E0*RS(2,3)*RS(2
     #,1)+2.E0/33075.E0*RS(2,3)*RS(2,2)+656.E0/99225.E0*RS(2,2)**2+11.E0
     #/99225.E0*RS(2,3)*RS(2,4)-26.E0/33075.E0*RS(2,4)*RS(2,2)+22.E0/992
     #25.E0*RS(2,4)**2-124.E0/14175.E0*RS(2,1)**2-52.E0/33075.E0*RS(2,1)
     #*RS(2,4)+RS(2,3)**2/99225)*RS(1,2)+(2.E0/3969.E0*RS(2,4)**2-116.E0
     #/99225.E0*RS(2,1)*RS(2,2)+116.E0/99225.E0*RS(2,1)*RS(2,4)+RS(2,3)*
     #RS(2,2)/99225-RS(2,3)*RS(2,4)/99225-2.E0/3969.E0*RS(2,2)**2)*RS(1,
     #3)+(-RS(2,3)**2/99225-2.E0/33075.E0*RS(2,3)*RS(2,4)+26.E0/33075.E0
     #*RS(2,4)*RS(2,2)-22.E0/99225.E0*RS(2,2)**2-656.E0/99225.E0*RS(2,4)
     #**2-328.E0/99225.E0*RS(2,1)*RS(2,4)+52.E0/33075.E0*RS(2,1)*RS(2,2)
     #-11.E0/99225.E0*RS(2,3)*RS(2,2)+124.E0/14175.E0*RS(2,1)**2-26.E0/3
     #3075.E0*RS(2,3)*RS(2,1))*RS(1,4))*RS(1,1)**2
      s2 = s3+s4
      s4 = s2
      s7 = (-2.E0/33075.E0*RS(2,1)*RS(2,4)-656.E0/99225.E0*RS(2,1)**2-32
     #8.E0/99225.E0*RS(2,1)*RS(2,2)-11.E0/99225.E0*RS(2,3)*RS(2,4)+52.E0
     #/33075.E0*RS(2,3)*RS(2,2)+124.E0/14175.E0*RS(2,2)**2-RS(2,4)**2/99
     #225-26.E0/33075.E0*RS(2,4)*RS(2,2)+26.E0/33075.E0*RS(2,3)*RS(2,1)-
     #22.E0/99225.E0*RS(2,3)**2)*RS(1,2)**2+((76.E0/99225.E0*RS(2,1)**2-
     #76.E0/99225.E0*RS(2,3)**2-8.E0/14175.E0*RS(2,1)*RS(2,2)+22.E0/9922
     #5.E0*RS(2,1)*RS(2,4)+8.E0/14175.E0*RS(2,3)*RS(2,2)-22.E0/99225.E0*
     #RS(2,3)*RS(2,4))*RS(1,3)+(76.E0/99225.E0*RS(2,4)**2+22.E0/99225.E0
     #*RS(2,3)*RS(2,4)-22.E0/99225.E0*RS(2,3)*RS(2,2)-76.E0/99225.E0*RS(
     #2,2)**2+8.E0/14175.E0*RS(2,1)*RS(2,2)-8.E0/14175.E0*RS(2,1)*RS(2,4
     #))*RS(1,4))*RS(1,2)
      s6 = s7+(2.E0/3969.E0*RS(2,2)**2+RS(2,1)*RS(2,4)/99225-116.E0/9922
     #5.E0*RS(2,3)*RS(2,4)+116.E0/99225.E0*RS(2,3)*RS(2,2)-RS(2,1)*RS(2,
     #2)/99225-2.E0/3969.E0*RS(2,4)**2)*RS(1,3)**2+(-76.E0/99225.E0*RS(2
     #,1)**2+76.E0/99225.E0*RS(2,3)**2-22.E0/99225.E0*RS(2,1)*RS(2,2)-8.
     #E0/14175.E0*RS(2,3)*RS(2,4)+8.E0/14175.E0*RS(2,1)*RS(2,4)+22.E0/99
     #225.E0*RS(2,3)*RS(2,2))*RS(1,4)*RS(1,3)+(11.E0/99225.E0*RS(2,3)*RS
     #(2,2)+26.E0/33075.E0*RS(2,4)*RS(2,2)+2.E0/33075.E0*RS(2,1)*RS(2,2)
     #-26.E0/33075.E0*RS(2,3)*RS(2,1)+22.E0/99225.E0*RS(2,3)**2+328.E0/9
     #9225.E0*RS(2,1)*RS(2,4)-52.E0/33075.E0*RS(2,3)*RS(2,4)+RS(2,2)**2/
     #99225-124.E0/14175.E0*RS(2,4)**2+656.E0/99225.E0*RS(2,1)**2)*RS(1,
     #4)**2
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(38.E0/99225.E0*RS(2,1)*RS(2,4)+124.E0/14175.E0*RS(2,3)*RS
     #(2,2)-38.E0/99225.E0*RS(2,3)*RS(2,4)-328.E0/99225.E0*RS(2,1)**2+32
     #8.E0/99225.E0*RS(2,3)**2-124.E0/14175.E0*RS(2,1)*RS(2,2))*RS(1,2)*
     #*3+((26.E0/33075.E0*RS(2,4)*RS(2,2)+22.E0/99225.E0*RS(2,1)**2+RS(2
     #,4)**2/99225+11.E0/99225.E0*RS(2,1)*RS(2,4)+2.E0/33075.E0*RS(2,3)*
     #RS(2,4)-26.E0/33075.E0*RS(2,3)*RS(2,1)+328.E0/99225.E0*RS(2,3)*RS(
     #2,2)-52.E0/33075.E0*RS(2,1)*RS(2,2)+656.E0/99225.E0*RS(2,3)**2-124
     #.E0/14175.E0*RS(2,2)**2)*RS(1,3)+(116.E0/99225.E0*RS(2,1)*RS(2,2)-
     #RS(2,1)*RS(2,4)/99225+RS(2,3)*RS(2,4)/99225-2.E0/3969.E0*RS(2,3)**
     #2-116.E0/99225.E0*RS(2,3)*RS(2,2)+2.E0/3969.E0*RS(2,1)**2)*RS(1,4)
     #)*RS(1,2)**2
      s3 = s1
      s5 = ((-26.E0/33075.E0*RS(2,3)*RS(2,1)-11.E0/99225.E0*RS(2,1)*RS(2
     #,4)+26.E0/33075.E0*RS(2,4)*RS(2,2)-328.E0/99225.E0*RS(2,3)*RS(2,2)
     #-22.E0/99225.E0*RS(2,4)**2-2.E0/33075.E0*RS(2,1)*RS(2,2)+52.E0/330
     #75.E0*RS(2,3)*RS(2,4)-656.E0/99225.E0*RS(2,2)**2+124.E0/14175.E0*R
     #S(2,3)**2-RS(2,1)**2/99225)*RS(1,3)**2+(76.E0/99225.E0*RS(2,2)**2-
     #22.E0/99225.E0*RS(2,1)*RS(2,4)-8.E0/14175.E0*RS(2,3)*RS(2,2)-76.E0
     #/99225.E0*RS(2,4)**2+22.E0/99225.E0*RS(2,1)*RS(2,2)+8.E0/14175.E0*
     #RS(2,3)*RS(2,4))*RS(1,4)*RS(1,3)+(-116.E0/99225.E0*RS(2,1)*RS(2,4)
     #+2.E0/3969.E0*RS(2,3)**2-2.E0/3969.E0*RS(2,1)**2+RS(2,1)*RS(2,2)/9
     #9225+116.E0/99225.E0*RS(2,3)*RS(2,4)-RS(2,3)*RS(2,2)/99225)*RS(1,4
     #)**2)*RS(1,2)
      s6 = (38.E0/99225.E0*RS(2,1)*RS(2,2)-38.E0/99225.E0*RS(2,1)*RS(2,4
     #)+328.E0/99225.E0*RS(2,4)**2-124.E0/14175.E0*RS(2,3)*RS(2,2)-328.E
     #0/99225.E0*RS(2,2)**2+124.E0/14175.E0*RS(2,3)*RS(2,4))*RS(1,3)**3
      s4 = s5+s6
      s2 = s3+s4
      t0 = s2+(11.E0/99225.E0*RS(2,1)*RS(2,2)+26.E0/33075.E0*RS(2,3)*RS(
     #2,1)+22.E0/99225.E0*RS(2,2)**2+328.E0/99225.E0*RS(2,3)*RS(2,4)+2.E
     #0/33075.E0*RS(2,1)*RS(2,4)+656.E0/99225.E0*RS(2,4)**2+RS(2,1)**2/9
     #9225-52.E0/33075.E0*RS(2,3)*RS(2,2)-26.E0/33075.E0*RS(2,4)*RS(2,2)
     #-124.E0/14175.E0*RS(2,3)**2)*RS(1,4)*RS(1,3)**2+(-2.E0/33075.E0*RS
     #(2,3)*RS(2,2)+26.E0/33075.E0*RS(2,3)*RS(2,1)-328.E0/99225.E0*RS(2,
     #3)*RS(2,4)-11.E0/99225.E0*RS(2,1)*RS(2,2)-656.E0/99225.E0*RS(2,3)*
     #*2+52.E0/33075.E0*RS(2,1)*RS(2,4)-22.E0/99225.E0*RS(2,1)**2+124.E0
     #/14175.E0*RS(2,4)**2-RS(2,2)**2/99225-26.E0/33075.E0*RS(2,4)*RS(2,
     #2))*RS(1,4)**2*RS(1,3)+(-124.E0/14175.E0*RS(2,3)*RS(2,4)+124.E0/14
     #175.E0*RS(2,1)*RS(2,4)-328.E0/99225.E0*RS(2,3)**2+328.E0/99225.E0*
     #RS(2,1)**2+38.E0/99225.E0*RS(2,3)*RS(2,2)-38.E0/99225.E0*RS(2,1)*R
     #S(2,2))*RS(1,4)**3
	ARRS(5) = T0

C	-----------------------
C	COMPUTE FOR FACTOR ASSR
C	-----------------------
      s3 = (-RS(1,1)*RS(1,2)/30+RS(1,1)*RS(1,4)/30-RS(1,2)**2/120+RS(1,4
     #)**2/120)*RS(2,1)**3
      s4 = ((-RS(1,1)*RS(1,2)/120+997.E0/793800.E0*RS(1,1)*RS(1,3)+121.E
     #0/793800.E0*RS(1,4)*RS(1,2)-RS(1,2)**2/90-167.E0/198450.E0*RS(1,1)
     #*RS(1,4)+121.E0/793800.E0*RS(1,4)**2+149.E0/396900.E0*RS(1,2)*RS(1
     #,3)+RS(1,1)**2/30+73.E0/0.15876E7*RS(1,4)*RS(1,3)-17.E0/317520.E0*
     #RS(1,3)**2)*RS(2,2)+(149.E0/396900.E0*RS(1,4)**2+997.E0/793800.E0*
     #RS(1,1)*RS(1,4)-17.E0/317520.E0*RS(1,4)*RS(1,3)+17.E0/317520.E0*RS
     #(1,2)*RS(1,3)-997.E0/793800.E0*RS(1,1)*RS(1,2)-149.E0/396900.E0*RS
     #(1,2)**2)*RS(2,3)+(RS(1,1)*RS(1,4)/120-149.E0/396900.E0*RS(1,4)*RS
     #(1,3)+167.E0/198450.E0*RS(1,1)*RS(1,2)-997.E0/793800.E0*RS(1,1)*RS
     #(1,3)-73.E0/0.15876E7*RS(1,2)*RS(1,3)-RS(1,1)**2/30+17.E0/317520.E
     #0*RS(1,3)**2-121.E0/793800.E0*RS(1,2)**2-121.E0/793800.E0*RS(1,4)*
     #RS(1,2)+RS(1,4)**2/90)*RS(2,4))*RS(2,1)**2
      s2 = s3+s4
      s4 = s2
      s7 = (RS(1,1)*RS(1,2)/180-RS(1,2)**2/120+RS(1,1)**2/60+149.E0/1984
     #50.E0*RS(1,4)*RS(1,2)-131.E0/0.15876E7*RS(1,4)**2+121.E0/396900.E0
     #*RS(1,1)*RS(1,4)-131.E0/0.15876E7*RS(1,3)**2+121.E0/396900.E0*RS(1
     #,2)*RS(1,3)+149.E0/198450.E0*RS(1,1)*RS(1,3)-13.E0/396900.E0*RS(1,
     #4)*RS(1,3))*RS(2,2)**2+((-RS(1,4)*RS(1,3)/26460-121.E0/396900.E0*R
     #S(1,2)**2+73.E0/793800.E0*RS(1,1)*RS(1,4)-13.E0/396900.E0*RS(1,4)*
     #*2-17.E0/158760.E0*RS(1,1)*RS(1,3)+47.E0/793800.E0*RS(1,4)*RS(1,2)
     #-17.E0/158760.E0*RS(1,3)**2-149.E0/198450.E0*RS(1,1)*RS(1,2))*RS(2
     #,3)+(-47.E0/793800.E0*RS(1,2)*RS(1,3)+121.E0/396900.E0*RS(1,1)*RS(
     #1,4)+149.E0/198450.E0*RS(1,4)**2+47.E0/793800.E0*RS(1,4)*RS(1,3)-1
     #21.E0/396900.E0*RS(1,1)*RS(1,2)-149.E0/198450.E0*RS(1,2)**2)*RS(2,
     #4))*RS(2,2)
      s6 = s7+(-131.E0/0.15876E7*RS(1,4)**2+17.E0/158760.E0*RS(1,1)*RS(1
     #,2)-17.E0/158760.E0*RS(1,4)*RS(1,3)+131.E0/0.15876E7*RS(1,2)**2-17
     #.E0/158760.E0*RS(1,1)*RS(1,4)+17.E0/158760.E0*RS(1,2)*RS(1,3))*RS(
     #2,3)**2+(149.E0/198450.E0*RS(1,1)*RS(1,4)-73.E0/793800.E0*RS(1,1)*
     #RS(1,2)+RS(1,2)*RS(1,3)/26460+17.E0/158760.E0*RS(1,3)**2+13.E0/396
     #900.E0*RS(1,2)**2+17.E0/158760.E0*RS(1,1)*RS(1,3)-47.E0/793800.E0*
     #RS(1,4)*RS(1,2)+121.E0/396900.E0*RS(1,4)**2)*RS(2,4)*RS(2,3)+(-121
     #.E0/396900.E0*RS(1,1)*RS(1,2)-121.E0/396900.E0*RS(1,4)*RS(1,3)+131
     #.E0/0.15876E7*RS(1,3)**2+13.E0/396900.E0*RS(1,2)*RS(1,3)-149.E0/19
     #8450.E0*RS(1,1)*RS(1,3)+RS(1,4)**2/120-RS(1,1)*RS(1,4)/180-149.E0/
     #198450.E0*RS(1,4)*RS(1,2)+131.E0/0.15876E7*RS(1,2)**2-RS(1,1)**2/6
     #0)*RS(2,4)**2
      s7 = RS(2,1)
      s5 = s6*s7
      s3 = s4+s5
      s4 = s3
      s6 = (121.E0/793800.E0*RS(1,3)**2+73.E0/0.15876E7*RS(1,4)*RS(1,3)+
     #121.E0/793800.E0*RS(1,1)*RS(1,3)+997.E0/793800.E0*RS(1,4)*RS(1,2)+
     #149.E0/396900.E0*RS(1,1)*RS(1,4)+RS(1,1)**2/180-167.E0/198450.E0*R
     #S(1,2)*RS(1,3)-17.E0/317520.E0*RS(1,4)**2+RS(1,1)*RS(1,2)/120)*RS(
     #2,2)**3
      s7 = ((23.E0/317520.E0*RS(1,4)*RS(1,3)+361.E0/0.15876E7*RS(1,4)*RS
     #(1,2)-RS(1,4)**2/52920-131.E0/0.15876E7*RS(1,1)*RS(1,3)+149.E0/198
     #450.E0*RS(1,3)**2+167.E0/198450.E0*RS(1,2)**2+RS(1,1)*RS(1,4)/7560
     #0-121.E0/793800.E0*RS(1,1)*RS(1,2)+121.E0/793800.E0*RS(1,2)*RS(1,3
     #))*RS(2,3)+(-73.E0/0.15876E7*RS(1,3)**2-149.E0/396900.E0*RS(1,1)*R
     #S(1,2)-2.E0/11025.E0*RS(1,2)*RS(1,3)-131.E0/0.15876E7*RS(1,1)*RS(1
     #,4)-73.E0/0.15876E7*RS(1,1)*RS(1,3)-17.E0/317520.E0*RS(1,4)*RS(1,2
     #)-RS(1,4)*RS(1,3)/52920-997.E0/793800.E0*RS(1,2)**2-17.E0/158760.E
     #0*RS(1,4)**2)*RS(2,4))*RS(2,2)**2
      s5 = s6+s7
      s1 = s4+s5
      s3 = s1
      s5 = ((-RS(1,1)*RS(1,4)/52920-149.E0/396900.E0*RS(1,2)*RS(1,3)+131
     #.E0/0.15876E7*RS(1,1)*RS(1,2)+997.E0/793800.E0*RS(1,3)**2-17.E0/31
     #7520.E0*RS(1,1)*RS(1,3)+167.E0/0.15876E7*RS(1,4)*RS(1,2)+RS(1,4)**
     #2/75600+31.E0/113400.E0*RS(1,4)*RS(1,3)-121.E0/396900.E0*RS(1,2)**
     #2)*RS(2,3)**2+(47.E0/396900.E0*RS(1,4)*RS(1,3)-13.E0/396900.E0*RS(
     #1,1)*RS(1,4)+13.E0/396900.E0*RS(1,1)*RS(1,2)-47.E0/396900.E0*RS(1,
     #2)*RS(1,3)-73.E0/793800.E0*RS(1,2)**2+73.E0/793800.E0*RS(1,4)**2)*
     #RS(2,4)*RS(2,3)+(149.E0/396900.E0*RS(1,1)*RS(1,4)+131.E0/0.15876E7
     #*RS(1,1)*RS(1,2)+17.E0/158760.E0*RS(1,2)**2+73.E0/0.15876E7*RS(1,3
     #)**2+RS(1,2)*RS(1,3)/52920+17.E0/317520.E0*RS(1,4)*RS(1,2)+997.E0/
     #793800.E0*RS(1,4)**2+2.E0/11025.E0*RS(1,4)*RS(1,3)+73.E0/0.15876E7
     #*RS(1,1)*RS(1,3))*RS(2,4)**2)*RS(2,2)
      s6 = (-997.E0/793800.E0*RS(1,2)*RS(1,3)+997.E0/793800.E0*RS(1,4)*R
     #S(1,3)-17.E0/317520.E0*RS(1,1)*RS(1,4)+149.E0/396900.E0*RS(1,4)**2
     #+17.E0/317520.E0*RS(1,1)*RS(1,2)-149.E0/396900.E0*RS(1,2)**2)*RS(2
     #,3)**3
      s4 = s5+s6
      s2 = s3+s4
      s3 = s2+(121.E0/396900.E0*RS(1,4)**2-31.E0/113400.E0*RS(1,2)*RS(1,
     #3)-997.E0/793800.E0*RS(1,3)**2-RS(1,2)**2/75600-167.E0/0.15876E7*R
     #S(1,4)*RS(1,2)+17.E0/317520.E0*RS(1,1)*RS(1,3)+149.E0/396900.E0*RS
     #(1,4)*RS(1,3)+RS(1,1)*RS(1,2)/52920-131.E0/0.15876E7*RS(1,1)*RS(1,
     #4))*RS(2,4)*RS(2,3)**2
      t0 = s3+(-361.E0/0.15876E7*RS(1,4)*RS(1,2)-RS(1,1)*RS(1,2)/75600-2
     #3.E0/317520.E0*RS(1,2)*RS(1,3)-167.E0/198450.E0*RS(1,4)**2-121.E0/
     #793800.E0*RS(1,4)*RS(1,3)-149.E0/198450.E0*RS(1,3)**2+131.E0/0.158
     #76E7*RS(1,1)*RS(1,3)+121.E0/793800.E0*RS(1,1)*RS(1,4)+RS(1,2)**2/5
     #2920)*RS(2,4)**2*RS(2,3)+(-121.E0/793800.E0*RS(1,3)**2-121.E0/7938
     #00.E0*RS(1,1)*RS(1,3)+167.E0/198450.E0*RS(1,4)*RS(1,3)-RS(1,1)**2/
     #180-73.E0/0.15876E7*RS(1,2)*RS(1,3)+17.E0/317520.E0*RS(1,2)**2-RS(
     #1,1)*RS(1,4)/120-149.E0/396900.E0*RS(1,1)*RS(1,2)-997.E0/793800.E0
     #*RS(1,4)*RS(1,2))*RS(2,4)**3
	ASSR(1) = T0

      s3 = (-149.E0/396900.E0*RS(1,2)*RS(1,3)-997.E0/793800.E0*RS(1,1)*R
     #S(1,3)-RS(1,2)**2/180-121.E0/793800.E0*RS(1,4)**2-121.E0/793800.E0
     #*RS(1,2)*RS(1,4)-73.E0/0.15876E7*RS(1,3)*RS(1,4)-RS(1,2)*RS(1,1)/1
     #20+167.E0/198450.E0*RS(1,1)*RS(1,4)+17.E0/317520.E0*RS(1,3)**2)*RS
     #(2,1)**3
      s6 = (-RS(1,2)*RS(1,1)/180+13.E0/396900.E0*RS(1,3)*RS(1,4)-149.E0/
     #198450.E0*RS(1,2)*RS(1,4)-149.E0/198450.E0*RS(1,1)*RS(1,3)+RS(1,1)
     #**2/120-RS(1,2)**2/60-121.E0/396900.E0*RS(1,1)*RS(1,4)-121.E0/3969
     #00.E0*RS(1,2)*RS(1,3)+131.E0/0.15876E7*RS(1,3)**2+131.E0/0.15876E7
     #*RS(1,4)**2)*RS(2,2)
      s7 = (997.E0/793800.E0*RS(1,1)**2+RS(1,3)*RS(1,4)/52920+73.E0/0.15
     #876E7*RS(1,2)*RS(1,4)+131.E0/0.15876E7*RS(1,2)*RS(1,3)+2.E0/11025.
     #E0*RS(1,1)*RS(1,4)+149.E0/396900.E0*RS(1,2)*RS(1,1)+17.E0/317520.E
     #0*RS(1,1)*RS(1,3)+73.E0/0.15876E7*RS(1,4)**2+17.E0/158760.E0*RS(1,
     #3)**2)*RS(2,3)+(131.E0/0.15876E7*RS(1,2)*RS(1,4)-149.E0/198450.E0*
     #RS(1,4)**2-23.E0/317520.E0*RS(1,3)*RS(1,4)-RS(1,2)*RS(1,3)/75600-1
     #21.E0/793800.E0*RS(1,1)*RS(1,4)-361.E0/0.15876E7*RS(1,1)*RS(1,3)+1
     #21.E0/793800.E0*RS(1,2)*RS(1,1)+RS(1,3)**2/52920-167.E0/198450.E0*
     #RS(1,1)**2)*RS(2,4)
      s5 = s6+s7
      s6 = RS(2,1)**2
      s4 = s5*s6
      s2 = s3+s4
      s4 = s2
      s7 = (-73.E0/0.15876E7*RS(1,3)*RS(1,4)-121.E0/793800.E0*RS(1,1)*RS
     #(1,3)+RS(1,1)**2/90+RS(1,2)*RS(1,1)/120-121.E0/793800.E0*RS(1,3)**
     #2-RS(1,2)**2/30+17.E0/317520.E0*RS(1,4)**2+167.E0/198450.E0*RS(1,2
     #)*RS(1,3)-997.E0/793800.E0*RS(1,2)*RS(1,4)-149.E0/396900.E0*RS(1,1
     #)*RS(1,4))*RS(2,2)**2+((121.E0/396900.E0*RS(1,2)*RS(1,1)+149.E0/19
     #8450.E0*RS(1,1)**2-47.E0/793800.E0*RS(1,3)*RS(1,4)-121.E0/396900.E
     #0*RS(1,2)*RS(1,3)+47.E0/793800.E0*RS(1,1)*RS(1,4)-149.E0/198450.E0
     #*RS(1,3)**2)*RS(2,3)+(13.E0/396900.E0*RS(1,3)**2+121.E0/396900.E0*
     #RS(1,1)**2+17.E0/158760.E0*RS(1,4)**2+149.E0/198450.E0*RS(1,2)*RS(
     #1,1)+RS(1,3)*RS(1,4)/26460-73.E0/793800.E0*RS(1,2)*RS(1,3)+17.E0/1
     #58760.E0*RS(1,2)*RS(1,4)-47.E0/793800.E0*RS(1,1)*RS(1,3))*RS(2,4))
     #*RS(2,2)
      s6 = s7+(-131.E0/0.15876E7*RS(1,2)*RS(1,1)-RS(1,1)*RS(1,4)/52920-7
     #3.E0/0.15876E7*RS(1,4)**2-2.E0/11025.E0*RS(1,3)*RS(1,4)-149.E0/396
     #900.E0*RS(1,2)*RS(1,3)-17.E0/158760.E0*RS(1,1)**2-17.E0/317520.E0*
     #RS(1,1)*RS(1,3)-73.E0/0.15876E7*RS(1,2)*RS(1,4)-997.E0/793800.E0*R
     #S(1,3)**2)*RS(2,3)**2+(73.E0/793800.E0*RS(1,1)**2-13.E0/396900.E0*
     #RS(1,2)*RS(1,1)-47.E0/396900.E0*RS(1,3)*RS(1,4)+47.E0/396900.E0*RS
     #(1,1)*RS(1,4)+13.E0/396900.E0*RS(1,2)*RS(1,3)-73.E0/793800.E0*RS(1
     #,3)**2)*RS(2,4)*RS(2,3)+(-131.E0/0.15876E7*RS(1,2)*RS(1,1)+RS(1,2)
     #*RS(1,3)/52920-167.E0/0.15876E7*RS(1,1)*RS(1,3)+17.E0/317520.E0*RS
     #(1,2)*RS(1,4)-997.E0/793800.E0*RS(1,4)**2+121.E0/396900.E0*RS(1,1)
     #**2+149.E0/396900.E0*RS(1,1)*RS(1,4)-31.E0/113400.E0*RS(1,3)*RS(1,
     #4)-RS(1,3)**2/75600)*RS(2,4)**2
      s7 = RS(2,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(-RS(1,2)*RS(1,3)/30-RS(1,3)**2/120+RS(1,1)**2/120+RS(1,2)
     #*RS(1,1)/30)*RS(2,2)**3+((997.E0/793800.E0*RS(1,2)*RS(1,4)+121.E0/
     #793800.E0*RS(1,1)*RS(1,3)+149.E0/396900.E0*RS(1,3)*RS(1,4)-RS(1,3)
     #**2/90+RS(1,2)**2/30-17.E0/317520.E0*RS(1,4)**2+121.E0/793800.E0*R
     #S(1,1)**2+73.E0/0.15876E7*RS(1,1)*RS(1,4)-RS(1,2)*RS(1,3)/120-167.
     #E0/198450.E0*RS(1,2)*RS(1,1))*RS(2,3)+(-997.E0/793800.E0*RS(1,2)*R
     #S(1,3)+149.E0/396900.E0*RS(1,1)**2-17.E0/317520.E0*RS(1,1)*RS(1,4)
     #+17.E0/317520.E0*RS(1,3)*RS(1,4)-149.E0/396900.E0*RS(1,3)**2+997.E
     #0/793800.E0*RS(1,2)*RS(1,1))*RS(2,4))*RS(2,2)**2
      s3 = s1
      s5 = ((RS(1,2)*RS(1,3)/180-13.E0/396900.E0*RS(1,1)*RS(1,4)-RS(1,3)
     #**2/120+RS(1,2)**2/60+121.E0/396900.E0*RS(1,3)*RS(1,4)-131.E0/0.15
     #876E7*RS(1,4)**2+149.E0/198450.E0*RS(1,2)*RS(1,4)+121.E0/396900.E0
     #*RS(1,2)*RS(1,1)-131.E0/0.15876E7*RS(1,1)**2+149.E0/198450.E0*RS(1
     #,1)*RS(1,3))*RS(2,3)**2+(73.E0/793800.E0*RS(1,2)*RS(1,1)-RS(1,1)*R
     #S(1,4)/26460-17.E0/158760.E0*RS(1,4)**2-17.E0/158760.E0*RS(1,2)*RS
     #(1,4)+47.E0/793800.E0*RS(1,1)*RS(1,3)-13.E0/396900.E0*RS(1,1)**2-1
     #49.E0/198450.E0*RS(1,2)*RS(1,3)-121.E0/396900.E0*RS(1,3)**2)*RS(2,
     #4)*RS(2,3)+(-17.E0/158760.E0*RS(1,1)*RS(1,4)-17.E0/158760.E0*RS(1,
     #2)*RS(1,1)+131.E0/0.15876E7*RS(1,3)**2+17.E0/158760.E0*RS(1,2)*RS(
     #1,3)-131.E0/0.15876E7*RS(1,1)**2+17.E0/158760.E0*RS(1,3)*RS(1,4))*
     #RS(2,4)**2)*RS(2,2)
      s6 = (-17.E0/317520.E0*RS(1,1)**2+149.E0/396900.E0*RS(1,2)*RS(1,1)
     #+997.E0/793800.E0*RS(1,1)*RS(1,3)+121.E0/793800.E0*RS(1,4)**2+121.
     #E0/793800.E0*RS(1,2)*RS(1,4)+RS(1,2)**2/180+73.E0/0.15876E7*RS(1,1
     #)*RS(1,4)+RS(1,2)*RS(1,3)/120-167.E0/198450.E0*RS(1,3)*RS(1,4))*RS
     #(2,3)**3
      s4 = s5+s6
      s2 = s3+s4
      t0 = s2+(RS(1,2)*RS(1,1)/75600+23.E0/317520.E0*RS(1,1)*RS(1,4)+121
     #.E0/793800.E0*RS(1,3)*RS(1,4)+361.E0/0.15876E7*RS(1,1)*RS(1,3)+167
     #.E0/198450.E0*RS(1,3)**2-121.E0/793800.E0*RS(1,2)*RS(1,3)+149.E0/1
     #98450.E0*RS(1,4)**2-131.E0/0.15876E7*RS(1,2)*RS(1,4)-RS(1,1)**2/52
     #920)*RS(2,4)*RS(2,3)**2+(-121.E0/396900.E0*RS(1,3)**2+RS(1,1)**2/7
     #5600+31.E0/113400.E0*RS(1,1)*RS(1,4)+131.E0/0.15876E7*RS(1,2)*RS(1
     #,3)-149.E0/396900.E0*RS(1,3)*RS(1,4)-17.E0/317520.E0*RS(1,2)*RS(1,
     #4)-RS(1,2)*RS(1,1)/52920+167.E0/0.15876E7*RS(1,1)*RS(1,3)+997.E0/7
     #93800.E0*RS(1,4)**2)*RS(2,4)**2*RS(2,3)+(149.E0/396900.E0*RS(1,1)*
     #*2-17.E0/317520.E0*RS(1,2)*RS(1,1)-149.E0/396900.E0*RS(1,3)**2-997
     #.E0/793800.E0*RS(1,3)*RS(1,4)+997.E0/793800.E0*RS(1,1)*RS(1,4)+17.
     #E0/317520.E0*RS(1,2)*RS(1,3))*RS(2,4)**3
	ASSR(2) = T0

      s3 = (-149.E0/396900.E0*RS(1,4)**2-997.E0/793800.E0*RS(1,1)*RS(1,4
     #)-17.E0/317520.E0*RS(1,3)*RS(1,2)+17.E0/317520.E0*RS(1,3)*RS(1,4)+
     #149.E0/396900.E0*RS(1,2)**2+997.E0/793800.E0*RS(1,1)*RS(1,2))*RS(2
     #,1)**3
      s4 = ((-31.E0/113400.E0*RS(1,1)*RS(1,4)+149.E0/396900.E0*RS(1,1)*R
     #S(1,2)-RS(1,4)**2/75600-167.E0/0.15876E7*RS(1,4)*RS(1,2)+121.E0/39
     #6900.E0*RS(1,2)**2+RS(1,3)*RS(1,4)/52920-997.E0/793800.E0*RS(1,1)*
     #*2+17.E0/317520.E0*RS(1,1)*RS(1,3)-131.E0/0.15876E7*RS(1,3)*RS(1,2
     #))*RS(2,2)+(17.E0/158760.E0*RS(1,3)*RS(1,4)-131.E0/0.15876E7*RS(1,
     #2)**2-17.E0/158760.E0*RS(1,3)*RS(1,2)+131.E0/0.15876E7*RS(1,4)**2-
     #17.E0/158760.E0*RS(1,1)*RS(1,2)+17.E0/158760.E0*RS(1,1)*RS(1,4))*R
     #S(2,3)+(-149.E0/396900.E0*RS(1,1)*RS(1,4)-17.E0/317520.E0*RS(1,1)*
     #RS(1,3)-121.E0/396900.E0*RS(1,4)**2+997.E0/793800.E0*RS(1,1)**2+31
     #.E0/113400.E0*RS(1,1)*RS(1,2)+167.E0/0.15876E7*RS(1,4)*RS(1,2)+RS(
     #1,2)**2/75600-RS(1,3)*RS(1,2)/52920+131.E0/0.15876E7*RS(1,3)*RS(1,
     #4))*RS(2,4))*RS(2,1)**2
      s2 = s3+s4
      s4 = s2
      s7 = (-361.E0/0.15876E7*RS(1,4)*RS(1,2)+121.E0/793800.E0*RS(1,3)*R
     #S(1,2)-23.E0/317520.E0*RS(1,1)*RS(1,4)-167.E0/198450.E0*RS(1,2)**2
     #-RS(1,3)*RS(1,4)/75600+RS(1,4)**2/52920-149.E0/198450.E0*RS(1,1)**
     #2-121.E0/793800.E0*RS(1,1)*RS(1,2)+131.E0/0.15876E7*RS(1,1)*RS(1,3
     #))*RS(2,2)**2+((-47.E0/793800.E0*RS(1,4)*RS(1,2)+121.E0/396900.E0*
     #RS(1,2)**2-73.E0/793800.E0*RS(1,3)*RS(1,4)+RS(1,1)*RS(1,4)/26460+1
     #7.E0/158760.E0*RS(1,1)*RS(1,3)+17.E0/158760.E0*RS(1,1)**2+149.E0/1
     #98450.E0*RS(1,3)*RS(1,2)+13.E0/396900.E0*RS(1,4)**2)*RS(2,3)+(73.E
     #0/793800.E0*RS(1,2)**2-73.E0/793800.E0*RS(1,4)**2+13.E0/396900.E0*
     #RS(1,3)*RS(1,4)-13.E0/396900.E0*RS(1,3)*RS(1,2)+47.E0/396900.E0*RS
     #(1,1)*RS(1,2)-47.E0/396900.E0*RS(1,1)*RS(1,4))*RS(2,4))*RS(2,2)
      s6 = s7+(-997.E0/793800.E0*RS(1,3)*RS(1,4)+149.E0/396900.E0*RS(1,2
     #)**2-149.E0/396900.E0*RS(1,4)**2-17.E0/317520.E0*RS(1,1)*RS(1,2)+1
     #7.E0/317520.E0*RS(1,1)*RS(1,4)+997.E0/793800.E0*RS(1,3)*RS(1,2))*R
     #S(2,3)**2+(-17.E0/158760.E0*RS(1,1)*RS(1,3)-13.E0/396900.E0*RS(1,2
     #)**2+73.E0/793800.E0*RS(1,3)*RS(1,2)-149.E0/198450.E0*RS(1,3)*RS(1
     #,4)-RS(1,1)*RS(1,2)/26460-17.E0/158760.E0*RS(1,1)**2+47.E0/793800.
     #E0*RS(1,4)*RS(1,2)-121.E0/396900.E0*RS(1,4)**2)*RS(2,4)*RS(2,3)+(1
     #67.E0/198450.E0*RS(1,4)**2+23.E0/317520.E0*RS(1,1)*RS(1,2)-121.E0/
     #793800.E0*RS(1,3)*RS(1,4)-131.E0/0.15876E7*RS(1,1)*RS(1,3)+RS(1,3)
     #*RS(1,2)/75600+121.E0/793800.E0*RS(1,1)*RS(1,4)+361.E0/0.15876E7*R
     #S(1,4)*RS(1,2)+149.E0/198450.E0*RS(1,1)**2-RS(1,2)**2/52920)*RS(2,
     #4)**2
      s7 = RS(2,1)
      s5 = s6*s7
      s3 = s4+s5
      s4 = s3
      s6 = (17.E0/317520.E0*RS(1,4)**2+167.E0/198450.E0*RS(1,1)*RS(1,2)-
     #73.E0/0.15876E7*RS(1,1)*RS(1,4)-997.E0/793800.E0*RS(1,4)*RS(1,2)-1
     #21.E0/793800.E0*RS(1,1)*RS(1,3)-121.E0/793800.E0*RS(1,1)**2-149.E0
     #/396900.E0*RS(1,3)*RS(1,4)-RS(1,3)**2/180-RS(1,3)*RS(1,2)/120)*RS(
     #2,2)**3
      s7 = ((-121.E0/396900.E0*RS(1,1)*RS(1,2)-RS(1,3)*RS(1,2)/180-149.E
     #0/198450.E0*RS(1,4)*RS(1,2)+131.E0/0.15876E7*RS(1,4)**2-RS(1,3)**2
     #/60-121.E0/396900.E0*RS(1,3)*RS(1,4)+RS(1,2)**2/120-149.E0/198450.
     #E0*RS(1,1)*RS(1,3)+131.E0/0.15876E7*RS(1,1)**2+13.E0/396900.E0*RS(
     #1,1)*RS(1,4))*RS(2,3)+(73.E0/0.15876E7*RS(1,1)*RS(1,3)+149.E0/3969
     #00.E0*RS(1,3)*RS(1,2)+2.E0/11025.E0*RS(1,1)*RS(1,2)+17.E0/317520.E
     #0*RS(1,4)*RS(1,2)+RS(1,1)*RS(1,4)/52920+17.E0/158760.E0*RS(1,4)**2
     #+131.E0/0.15876E7*RS(1,3)*RS(1,4)+997.E0/793800.E0*RS(1,2)**2+73.E
     #0/0.15876E7*RS(1,1)**2)*RS(2,4))*RS(2,2)**2
      s5 = s6+s7
      s1 = s4+s5
      s3 = s1
      s5 = ((167.E0/198450.E0*RS(1,3)*RS(1,4)-997.E0/793800.E0*RS(1,1)*R
     #S(1,3)-RS(1,3)**2/30+17.E0/317520.E0*RS(1,1)**2-121.E0/793800.E0*R
     #S(1,4)*RS(1,2)-149.E0/396900.E0*RS(1,1)*RS(1,2)-73.E0/0.15876E7*RS
     #(1,1)*RS(1,4)+RS(1,3)*RS(1,2)/120-121.E0/793800.E0*RS(1,4)**2+RS(1
     #,2)**2/90)*RS(2,3)**2+(149.E0/198450.E0*RS(1,2)**2-121.E0/396900.E
     #0*RS(1,3)*RS(1,4)-149.E0/198450.E0*RS(1,4)**2-47.E0/793800.E0*RS(1
     #,1)*RS(1,4)+121.E0/396900.E0*RS(1,3)*RS(1,2)+47.E0/793800.E0*RS(1,
     #1)*RS(1,2))*RS(2,4)*RS(2,3)+(-RS(1,1)*RS(1,2)/52920-149.E0/396900.
     #E0*RS(1,3)*RS(1,4)-73.E0/0.15876E7*RS(1,1)*RS(1,3)-73.E0/0.15876E7
     #*RS(1,1)**2-997.E0/793800.E0*RS(1,4)**2-17.E0/158760.E0*RS(1,2)**2
     #-131.E0/0.15876E7*RS(1,3)*RS(1,2)-2.E0/11025.E0*RS(1,1)*RS(1,4)-17
     #.E0/317520.E0*RS(1,4)*RS(1,2))*RS(2,4)**2)*RS(2,2)
      s6 = (-RS(1,4)**2/120+RS(1,2)**2/120+RS(1,3)*RS(1,2)/30-RS(1,3)*RS
     #(1,4)/30)*RS(2,3)**3
      s4 = s5+s6
      s2 = s3+s4
      s3 = s2+(997.E0/793800.E0*RS(1,1)*RS(1,3)+149.E0/396900.E0*RS(1,1)
     #*RS(1,4)-167.E0/198450.E0*RS(1,3)*RS(1,2)+121.E0/793800.E0*RS(1,2)
     #**2-RS(1,3)*RS(1,4)/120-RS(1,4)**2/90+RS(1,3)**2/30-17.E0/317520.E
     #0*RS(1,1)**2+73.E0/0.15876E7*RS(1,1)*RS(1,2)+121.E0/793800.E0*RS(1
     #,4)*RS(1,2))*RS(2,4)*RS(2,3)**2
      t0 = s3+(121.E0/396900.E0*RS(1,3)*RS(1,2)-RS(1,4)**2/120+RS(1,3)**
     #2/60+149.E0/198450.E0*RS(1,1)*RS(1,3)+149.E0/198450.E0*RS(1,4)*RS(
     #1,2)+121.E0/396900.E0*RS(1,1)*RS(1,4)-131.E0/0.15876E7*RS(1,2)**2-
     #13.E0/396900.E0*RS(1,1)*RS(1,2)-131.E0/0.15876E7*RS(1,1)**2+RS(1,3
     #)*RS(1,4)/180)*RS(2,4)**2*RS(2,3)+(121.E0/793800.E0*RS(1,1)**2+RS(
     #1,3)*RS(1,4)/120+73.E0/0.15876E7*RS(1,1)*RS(1,2)+RS(1,3)**2/180-17
     #.E0/317520.E0*RS(1,2)**2+149.E0/396900.E0*RS(1,3)*RS(1,2)-167.E0/1
     #98450.E0*RS(1,1)*RS(1,4)+121.E0/793800.E0*RS(1,1)*RS(1,3)+997.E0/7
     #93800.E0*RS(1,4)*RS(1,2))*RS(2,4)**3
	ASSR(3) = T0

      s3 = (121.E0/793800.E0*RS(1,4)*RS(1,2)+149.E0/396900.E0*RS(1,3)*RS
     #(1,4)-17.E0/317520.E0*RS(1,3)**2+RS(1,4)**2/180+73.E0/0.15876E7*RS
     #(1,3)*RS(1,2)+121.E0/793800.E0*RS(1,2)**2+RS(1,1)*RS(1,4)/120+997.
     #E0/793800.E0*RS(1,1)*RS(1,3)-167.E0/198450.E0*RS(1,1)*RS(1,2))*RS(
     #2,1)**3
      s6 = (RS(1,3)*RS(1,4)/75600-121.E0/793800.E0*RS(1,1)*RS(1,4)-131.E
     #0/0.15876E7*RS(1,4)*RS(1,2)+121.E0/793800.E0*RS(1,1)*RS(1,2)+361.E
     #0/0.15876E7*RS(1,1)*RS(1,3)+167.E0/198450.E0*RS(1,1)**2+149.E0/198
     #450.E0*RS(1,2)**2+23.E0/317520.E0*RS(1,3)*RS(1,2)-RS(1,3)**2/52920
     #)*RS(2,2)
      s7 = (-17.E0/317520.E0*RS(1,1)*RS(1,3)-149.E0/396900.E0*RS(1,1)*RS
     #(1,4)-131.E0/0.15876E7*RS(1,3)*RS(1,4)-73.E0/0.15876E7*RS(1,2)**2-
     #997.E0/793800.E0*RS(1,1)**2-17.E0/158760.E0*RS(1,3)**2-RS(1,3)*RS(
     #1,2)/52920-2.E0/11025.E0*RS(1,1)*RS(1,2)-73.E0/0.15876E7*RS(1,4)*R
     #S(1,2))*RS(2,3)+(RS(1,4)**2/60-13.E0/396900.E0*RS(1,3)*RS(1,2)+121
     #.E0/396900.E0*RS(1,3)*RS(1,4)-RS(1,1)**2/120+149.E0/198450.E0*RS(1
     #,4)*RS(1,2)+149.E0/198450.E0*RS(1,1)*RS(1,3)+121.E0/396900.E0*RS(1
     #,1)*RS(1,2)-131.E0/0.15876E7*RS(1,2)**2+RS(1,1)*RS(1,4)/180-131.E0
     #/0.15876E7*RS(1,3)**2)*RS(2,4)
      s5 = s6+s7
      s6 = RS(2,1)**2
      s4 = s5*s6
      s2 = s3+s4
      s4 = s2
      s7 = (-RS(1,3)*RS(1,4)/52920+RS(1,3)**2/75600-149.E0/396900.E0*RS(
     #1,1)*RS(1,2)+997.E0/793800.E0*RS(1,2)**2+167.E0/0.15876E7*RS(1,1)*
     #RS(1,3)-17.E0/317520.E0*RS(1,4)*RS(1,2)-121.E0/396900.E0*RS(1,1)**
     #2+31.E0/113400.E0*RS(1,3)*RS(1,2)+131.E0/0.15876E7*RS(1,1)*RS(1,4)
     #)*RS(2,2)**2+((13.E0/396900.E0*RS(1,1)*RS(1,4)-47.E0/396900.E0*RS(
     #1,1)*RS(1,2)-73.E0/793800.E0*RS(1,1)**2+47.E0/396900.E0*RS(1,3)*RS
     #(1,2)+73.E0/793800.E0*RS(1,3)**2-13.E0/396900.E0*RS(1,3)*RS(1,4))*
     #RS(2,3)+(-RS(1,3)*RS(1,2)/26460+47.E0/793800.E0*RS(1,1)*RS(1,3)-13
     #.E0/396900.E0*RS(1,3)**2+73.E0/793800.E0*RS(1,3)*RS(1,4)-17.E0/158
     #760.E0*RS(1,2)**2-121.E0/396900.E0*RS(1,1)**2-17.E0/158760.E0*RS(1
     #,4)*RS(1,2)-149.E0/198450.E0*RS(1,1)*RS(1,4))*RS(2,4))*RS(2,2)
      s6 = s7+(131.E0/0.15876E7*RS(1,1)*RS(1,4)+2.E0/11025.E0*RS(1,3)*RS
     #(1,2)+17.E0/317520.E0*RS(1,1)*RS(1,3)+149.E0/396900.E0*RS(1,3)*RS(
     #1,4)+73.E0/0.15876E7*RS(1,2)**2+997.E0/793800.E0*RS(1,3)**2+73.E0/
     #0.15876E7*RS(1,4)*RS(1,2)+RS(1,1)*RS(1,2)/52920+17.E0/158760.E0*RS
     #(1,1)**2)*RS(2,3)**2+(47.E0/793800.E0*RS(1,3)*RS(1,2)-121.E0/39690
     #0.E0*RS(1,1)*RS(1,4)-149.E0/198450.E0*RS(1,1)**2+121.E0/396900.E0*
     #RS(1,3)*RS(1,4)+149.E0/198450.E0*RS(1,3)**2-47.E0/793800.E0*RS(1,1
     #)*RS(1,2))*RS(2,4)*RS(2,3)+(149.E0/396900.E0*RS(1,1)*RS(1,2)-RS(1,
     #1)**2/90-17.E0/317520.E0*RS(1,2)**2+997.E0/793800.E0*RS(1,4)*RS(1,
     #2)+73.E0/0.15876E7*RS(1,3)*RS(1,2)+RS(1,4)**2/30+121.E0/793800.E0*
     #RS(1,1)*RS(1,3)-RS(1,1)*RS(1,4)/120-167.E0/198450.E0*RS(1,3)*RS(1,
     #4)+121.E0/793800.E0*RS(1,3)**2)*RS(2,4)**2
      s7 = RS(2,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(-149.E0/396900.E0*RS(1,1)**2+17.E0/317520.E0*RS(1,1)*RS(1
     #,4)+997.E0/793800.E0*RS(1,3)*RS(1,2)+149.E0/396900.E0*RS(1,3)**2-9
     #97.E0/793800.E0*RS(1,1)*RS(1,2)-17.E0/317520.E0*RS(1,3)*RS(1,4))*R
     #S(2,2)**3+((-31.E0/113400.E0*RS(1,1)*RS(1,2)-RS(1,1)**2/75600+17.E
     #0/317520.E0*RS(1,4)*RS(1,2)+RS(1,1)*RS(1,4)/52920+121.E0/396900.E0
     #*RS(1,3)**2-167.E0/0.15876E7*RS(1,1)*RS(1,3)-997.E0/793800.E0*RS(1
     #,2)**2-131.E0/0.15876E7*RS(1,3)*RS(1,4)+149.E0/396900.E0*RS(1,3)*R
     #S(1,2))*RS(2,3)+(-17.E0/158760.E0*RS(1,3)*RS(1,4)+17.E0/158760.E0*
     #RS(1,1)*RS(1,4)-131.E0/0.15876E7*RS(1,3)**2+17.E0/158760.E0*RS(1,1
     #)*RS(1,2)+131.E0/0.15876E7*RS(1,1)**2-17.E0/158760.E0*RS(1,3)*RS(1
     #,2))*RS(2,4))*RS(2,2)**2
      s3 = s1
      s5 = ((-RS(1,1)*RS(1,4)/75600+131.E0/0.15876E7*RS(1,4)*RS(1,2)-121
     #.E0/793800.E0*RS(1,3)*RS(1,2)-149.E0/198450.E0*RS(1,2)**2-361.E0/0
     #.15876E7*RS(1,1)*RS(1,3)+121.E0/793800.E0*RS(1,3)*RS(1,4)-167.E0/1
     #98450.E0*RS(1,3)**2+RS(1,1)**2/52920-23.E0/317520.E0*RS(1,1)*RS(1,
     #2))*RS(2,3)**2+(17.E0/158760.E0*RS(1,2)**2-47.E0/793800.E0*RS(1,1)
     #*RS(1,3)+RS(1,1)*RS(1,2)/26460+121.E0/396900.E0*RS(1,3)**2+17.E0/1
     #58760.E0*RS(1,4)*RS(1,2)+149.E0/198450.E0*RS(1,3)*RS(1,4)-73.E0/79
     #3800.E0*RS(1,1)*RS(1,4)+13.E0/396900.E0*RS(1,1)**2)*RS(2,4)*RS(2,3
     #)+(-149.E0/396900.E0*RS(1,1)**2-17.E0/317520.E0*RS(1,3)*RS(1,2)+17
     #.E0/317520.E0*RS(1,1)*RS(1,2)-997.E0/793800.E0*RS(1,1)*RS(1,4)+149
     #.E0/396900.E0*RS(1,3)**2+997.E0/793800.E0*RS(1,3)*RS(1,4))*RS(2,4)
     #**2)*RS(2,2)
      s6 = (167.E0/198450.E0*RS(1,3)*RS(1,2)-73.E0/0.15876E7*RS(1,1)*RS(
     #1,2)+17.E0/317520.E0*RS(1,1)**2-997.E0/793800.E0*RS(1,1)*RS(1,3)-R
     #S(1,4)**2/180-121.E0/793800.E0*RS(1,2)**2-149.E0/396900.E0*RS(1,1)
     #*RS(1,4)-121.E0/793800.E0*RS(1,4)*RS(1,2)-RS(1,3)*RS(1,4)/120)*RS(
     #2,3)**3
      s4 = s5+s6
      s2 = s3+s4
      t0 = s2+(-149.E0/198450.E0*RS(1,4)*RS(1,2)-121.E0/396900.E0*RS(1,3
     #)*RS(1,2)-149.E0/198450.E0*RS(1,1)*RS(1,3)+131.E0/0.15876E7*RS(1,2
     #)**2+RS(1,3)**2/120-121.E0/396900.E0*RS(1,1)*RS(1,4)-RS(1,3)*RS(1,
     #4)/180+13.E0/396900.E0*RS(1,1)*RS(1,2)-RS(1,4)**2/60+131.E0/0.1587
     #6E7*RS(1,1)**2)*RS(2,4)*RS(2,3)**2+(-149.E0/396900.E0*RS(1,3)*RS(1
     #,2)+167.E0/198450.E0*RS(1,1)*RS(1,4)+RS(1,3)*RS(1,4)/120-997.E0/79
     #3800.E0*RS(1,4)*RS(1,2)+RS(1,3)**2/90-RS(1,4)**2/30-121.E0/793800.
     #E0*RS(1,1)*RS(1,3)-73.E0/0.15876E7*RS(1,1)*RS(1,2)+17.E0/317520.E0
     #*RS(1,2)**2-121.E0/793800.E0*RS(1,1)**2)*RS(2,4)**2*RS(2,3)+(RS(1,
     #3)**2/120-RS(1,1)**2/120+RS(1,3)*RS(1,4)/30-RS(1,1)*RS(1,4)/30)*RS
     #(2,4)**3
	ASSR(4) = T0

      s3 = (328.E0/99225.E0*RS(1,4)**2-38.E0/99225.E0*RS(1,3)*RS(1,4)+12
     #4.E0/14175.E0*RS(1,1)*RS(1,4)+38.E0/99225.E0*RS(1,3)*RS(1,2)-124.E
     #0/14175.E0*RS(1,1)*RS(1,2)-328.E0/99225.E0*RS(1,2)**2)*RS(2,1)**3
      s4 = ((-2.E0/33075.E0*RS(1,3)*RS(1,2)+52.E0/33075.E0*RS(1,1)*RS(1,
     #4)-22.E0/99225.E0*RS(1,4)**2-328.E0/99225.E0*RS(1,1)*RS(1,2)-RS(1,
     #3)**2/99225-656.E0/99225.E0*RS(1,2)**2+26.E0/33075.E0*RS(1,4)*RS(1
     #,2)+124.E0/14175.E0*RS(1,1)**2-26.E0/33075.E0*RS(1,3)*RS(1,1)-11.E
     #0/99225.E0*RS(1,3)*RS(1,4))*RS(2,2)+(-RS(1,3)*RS(1,2)/99225+116.E0
     #/99225.E0*RS(1,1)*RS(1,2)-116.E0/99225.E0*RS(1,1)*RS(1,4)-2.E0/396
     #9.E0*RS(1,4)**2+2.E0/3969.E0*RS(1,2)**2+RS(1,3)*RS(1,4)/99225)*RS(
     #2,3)+(11.E0/99225.E0*RS(1,3)*RS(1,2)+26.E0/33075.E0*RS(1,3)*RS(1,1
     #)-124.E0/14175.E0*RS(1,1)**2-26.E0/33075.E0*RS(1,4)*RS(1,2)+328.E0
     #/99225.E0*RS(1,1)*RS(1,4)+2.E0/33075.E0*RS(1,3)*RS(1,4)-52.E0/3307
     #5.E0*RS(1,1)*RS(1,2)+22.E0/99225.E0*RS(1,2)**2+RS(1,3)**2/99225+65
     #6.E0/99225.E0*RS(1,4)**2)*RS(2,4))*RS(2,1)**2
      s2 = s3+s4
      s4 = s2
      s7 = (11.E0/99225.E0*RS(1,3)*RS(1,4)+328.E0/99225.E0*RS(1,1)*RS(1,
     #2)+RS(1,4)**2/99225+656.E0/99225.E0*RS(1,1)**2+26.E0/33075.E0*RS(1
     #,4)*RS(1,2)-124.E0/14175.E0*RS(1,2)**2-26.E0/33075.E0*RS(1,3)*RS(1
     #,1)+22.E0/99225.E0*RS(1,3)**2+2.E0/33075.E0*RS(1,1)*RS(1,4)-52.E0/
     #33075.E0*RS(1,3)*RS(1,2))*RS(2,2)**2+((22.E0/99225.E0*RS(1,3)*RS(1
     #,4)+76.E0/99225.E0*RS(1,3)**2+8.E0/14175.E0*RS(1,1)*RS(1,2)-76.E0/
     #99225.E0*RS(1,1)**2-8.E0/14175.E0*RS(1,3)*RS(1,2)-22.E0/99225.E0*R
     #S(1,1)*RS(1,4))*RS(2,3)+(22.E0/99225.E0*RS(1,3)*RS(1,2)-76.E0/9922
     #5.E0*RS(1,4)**2-22.E0/99225.E0*RS(1,3)*RS(1,4)-8.E0/14175.E0*RS(1,
     #1)*RS(1,2)+76.E0/99225.E0*RS(1,2)**2+8.E0/14175.E0*RS(1,1)*RS(1,4)
     #)*RS(2,4))*RS(2,2)
      s6 = s7+(116.E0/99225.E0*RS(1,3)*RS(1,4)+2.E0/3969.E0*RS(1,4)**2-R
     #S(1,1)*RS(1,4)/99225+RS(1,1)*RS(1,2)/99225-2.E0/3969.E0*RS(1,2)**2
     #-116.E0/99225.E0*RS(1,3)*RS(1,2))*RS(2,3)**2+(76.E0/99225.E0*RS(1,
     #1)**2+22.E0/99225.E0*RS(1,1)*RS(1,2)-76.E0/99225.E0*RS(1,3)**2-22.
     #E0/99225.E0*RS(1,3)*RS(1,2)+8.E0/14175.E0*RS(1,3)*RS(1,4)-8.E0/141
     #75.E0*RS(1,1)*RS(1,4))*RS(2,4)*RS(2,3)+(-656.E0/99225.E0*RS(1,1)**
     #2-2.E0/33075.E0*RS(1,1)*RS(1,2)-328.E0/99225.E0*RS(1,1)*RS(1,4)+12
     #4.E0/14175.E0*RS(1,4)**2-26.E0/33075.E0*RS(1,4)*RS(1,2)+52.E0/3307
     #5.E0*RS(1,3)*RS(1,4)+26.E0/33075.E0*RS(1,3)*RS(1,1)-11.E0/99225.E0
     #*RS(1,3)*RS(1,2)-RS(1,2)**2/99225-22.E0/99225.E0*RS(1,3)**2)*RS(2,
     #4)**2
      s7 = RS(2,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(-328.E0/99225.E0*RS(1,3)**2+38.E0/99225.E0*RS(1,3)*RS(1,4
     #)+124.E0/14175.E0*RS(1,1)*RS(1,2)+328.E0/99225.E0*RS(1,1)**2-124.E
     #0/14175.E0*RS(1,3)*RS(1,2)-38.E0/99225.E0*RS(1,1)*RS(1,4))*RS(2,2)
     #**3+((124.E0/14175.E0*RS(1,2)**2-328.E0/99225.E0*RS(1,3)*RS(1,2)-R
     #S(1,4)**2/99225-22.E0/99225.E0*RS(1,1)**2-26.E0/33075.E0*RS(1,4)*R
     #S(1,2)-656.E0/99225.E0*RS(1,3)**2-2.E0/33075.E0*RS(1,3)*RS(1,4)+26
     #.E0/33075.E0*RS(1,3)*RS(1,1)-11.E0/99225.E0*RS(1,1)*RS(1,4)+52.E0/
     #33075.E0*RS(1,1)*RS(1,2))*RS(2,3)+(RS(1,1)*RS(1,4)/99225-116.E0/99
     #225.E0*RS(1,1)*RS(1,2)-RS(1,3)*RS(1,4)/99225+2.E0/3969.E0*RS(1,3)*
     #*2-2.E0/3969.E0*RS(1,1)**2+116.E0/99225.E0*RS(1,3)*RS(1,2))*RS(2,4
     #))*RS(2,2)**2
      s3 = s1
      s5 = ((22.E0/99225.E0*RS(1,4)**2+328.E0/99225.E0*RS(1,3)*RS(1,2)+6
     #56.E0/99225.E0*RS(1,2)**2+11.E0/99225.E0*RS(1,1)*RS(1,4)+26.E0/330
     #75.E0*RS(1,3)*RS(1,1)+2.E0/33075.E0*RS(1,1)*RS(1,2)-26.E0/33075.E0
     #*RS(1,4)*RS(1,2)-124.E0/14175.E0*RS(1,3)**2+RS(1,1)**2/99225-52.E0
     #/33075.E0*RS(1,3)*RS(1,4))*RS(2,3)**2+(22.E0/99225.E0*RS(1,1)*RS(1
     #,4)+76.E0/99225.E0*RS(1,4)**2-8.E0/14175.E0*RS(1,3)*RS(1,4)-22.E0/
     #99225.E0*RS(1,1)*RS(1,2)-76.E0/99225.E0*RS(1,2)**2+8.E0/14175.E0*R
     #S(1,3)*RS(1,2))*RS(2,4)*RS(2,3)+(RS(1,3)*RS(1,2)/99225+116.E0/9922
     #5.E0*RS(1,1)*RS(1,4)+2.E0/3969.E0*RS(1,1)**2-2.E0/3969.E0*RS(1,3)*
     #*2-RS(1,1)*RS(1,2)/99225-116.E0/99225.E0*RS(1,3)*RS(1,4))*RS(2,4)*
     #*2)*RS(2,2)
      s6 = (-124.E0/14175.E0*RS(1,3)*RS(1,4)+328.E0/99225.E0*RS(1,2)**2+
     #124.E0/14175.E0*RS(1,3)*RS(1,2)-328.E0/99225.E0*RS(1,4)**2-38.E0/9
     #9225.E0*RS(1,1)*RS(1,2)+38.E0/99225.E0*RS(1,1)*RS(1,4))*RS(2,3)**3
      s4 = s5+s6
      s2 = s3+s4
      t0 = s2+(-22.E0/99225.E0*RS(1,2)**2-2.E0/33075.E0*RS(1,1)*RS(1,4)+
     #52.E0/33075.E0*RS(1,3)*RS(1,2)-328.E0/99225.E0*RS(1,3)*RS(1,4)-RS(
     #1,1)**2/99225+124.E0/14175.E0*RS(1,3)**2-11.E0/99225.E0*RS(1,1)*RS
     #(1,2)+26.E0/33075.E0*RS(1,4)*RS(1,2)-656.E0/99225.E0*RS(1,4)**2-26
     #.E0/33075.E0*RS(1,3)*RS(1,1))*RS(2,4)*RS(2,3)**2+(26.E0/33075.E0*R
     #S(1,4)*RS(1,2)+328.E0/99225.E0*RS(1,3)*RS(1,4)-124.E0/14175.E0*RS(
     #1,4)**2+22.E0/99225.E0*RS(1,1)**2+11.E0/99225.E0*RS(1,1)*RS(1,2)+R
     #S(1,2)**2/99225-52.E0/33075.E0*RS(1,1)*RS(1,4)-26.E0/33075.E0*RS(1
     #,3)*RS(1,1)+2.E0/33075.E0*RS(1,3)*RS(1,2)+656.E0/99225.E0*RS(1,3)*
     #*2)*RS(2,4)**2*RS(2,3)+(124.E0/14175.E0*RS(1,3)*RS(1,4)-124.E0/141
     #75.E0*RS(1,1)*RS(1,4)-38.E0/99225.E0*RS(1,3)*RS(1,2)+328.E0/99225.
     #E0*RS(1,3)**2-328.E0/99225.E0*RS(1,1)**2+38.E0/99225.E0*RS(1,1)*RS
     #(1,2))*RS(2,4)**3
	ASSR(5) = T0

C	------------------------
C	COMPUTE FOR FACTOR ARRRR
C	------------------------
      s3 = (RS(2,2)/42-RS(2,4)/42)*RS(1,1)**5+((37517.E0/0.480249E8*RS(2
     #,4)+2.E0/105.E0*RS(2,2)-RS(2,1)/42-11357.E0/0.160083E8*RS(2,3))*RS
     #(1,2)+(11357.E0/0.160083E8*RS(2,2)-11357.E0/0.160083E8*RS(2,4))*RS
     #(1,3)+(RS(2,1)/42-37517.E0/0.480249E8*RS(2,2)-2.E0/105.E0*RS(2,4)+
     #11357.E0/0.160083E8*RS(2,3))*RS(1,4))*RS(1,1)**4
      s4 = s3
      s6 = ((-17.E0/190575.E0*RS(2,4)-499.E0/686070.E0*RS(2,3)+RS(2,2)/7
     #0-2.E0/105.E0*RS(2,1))*RS(1,2)**2+((881.E0/0.800415E7*RS(2,3)+499.
     #E0/686070.E0*RS(2,2)-494.E0/0.4002075E7*RS(2,4))*RS(1,3)+(17.E0/19
     #0575.E0*RS(2,2)-17.E0/190575.E0*RS(2,4))*RS(1,4))*RS(1,2)+(-881.E0
     #/0.800415E7*RS(2,2)+881.E0/0.800415E7*RS(2,4))*RS(1,3)**2+(494.E0/
     #0.4002075E7*RS(2,2)-881.E0/0.800415E7*RS(2,3)-499.E0/686070.E0*RS(
     #2,4))*RS(1,4)*RS(1,3)+(2.E0/105.E0*RS(2,1)+499.E0/686070.E0*RS(2,3
     #)+17.E0/190575.E0*RS(2,2)-RS(2,4)/70)*RS(1,4)**2)*RS(1,1)**3
      s8 = (-13.E0/23716.E0*RS(2,3)+RS(2,2)/105-RS(2,1)/70-13.E0/23716.E
     #0*RS(2,4))*RS(1,2)**3+((13.E0/23716.E0*RS(2,2)+137.E0/0.160083E7*R
     #S(2,3)-541.E0/0.160083E8*RS(2,4))*RS(1,3)+(19.E0/148225.E0*RS(2,4)
     #+463.E0/0.53361E7*RS(2,3)+13.E0/23716.E0*RS(2,2))*RS(1,4))*RS(1,2)
     #**2+((-59.E0/0.213444E7*RS(2,3)+337.E0/0.320166E8*RS(2,4)-137.E0/0
     #.160083E7*RS(2,2))*RS(1,3)**2+(212.E0/0.4002075E7*RS(2,4)-212.E0/0
     #.4002075E7*RS(2,2))*RS(1,4)*RS(1,3)+(-463.E0/0.53361E7*RS(2,3)-13.
     #E0/23716.E0*RS(2,4)-19.E0/148225.E0*RS(2,2))*RS(1,4)**2)*RS(1,2)+(
     #-59.E0/0.213444E7*RS(2,4)+59.E0/0.213444E7*RS(2,2))*RS(1,3)**3+(-3
     #37.E0/0.320166E8*RS(2,2)+59.E0/0.213444E7*RS(2,3)+137.E0/0.160083E
     #7*RS(2,4))*RS(1,4)*RS(1,3)**2+(541.E0/0.160083E8*RS(2,2)-137.E0/0.
     #160083E7*RS(2,3)-13.E0/23716.E0*RS(2,4))*RS(1,4)**2*RS(1,3)+(-RS(2
     #,4)/105+13.E0/23716.E0*RS(2,2)+13.E0/23716.E0*RS(2,3)+RS(2,1)/70)*
     #RS(1,4)**3
      s9 = RS(1,1)**2
      s7 = s8*s9
      s5 = s6+s7
      s2 = s4+s5
      s4 = s2
      s8 = (-RS(2,1)/105-499.E0/686070.E0*RS(2,4)-17.E0/190575.E0*RS(2,3
     #)+RS(2,2)/210)*RS(1,2)**4+((19.E0/148225.E0*RS(2,3)+17.E0/190575.E
     #0*RS(2,2)-2.E0/190575.E0*RS(2,4))*RS(1,3)+(254.E0/0.4002075E7*RS(2
     #,3)+499.E0/686070.E0*RS(2,2)+137.E0/0.160083E7*RS(2,4))*RS(1,4))*R
     #S(1,2)**3
      s7 = s8+((149.E0/0.800415E7*RS(2,4)-19.E0/148225.E0*RS(2,2)+137.E0
     #/0.160083E7*RS(2,3))*RS(1,3)**2+(-212.E0/0.4002075E7*RS(2,2)-191.E
     #0/0.160083E8*RS(2,4)-163.E0/0.53361E7*RS(2,3))*RS(1,4)*RS(1,3)+(-1
     #37.E0/0.160083E7*RS(2,2)+137.E0/0.160083E7*RS(2,4))*RS(1,4)**2)*RS
     #(1,2)**2+((76.E0/0.2401245E7*RS(2,4)-137.E0/0.160083E7*RS(2,2)+881
     #.E0/0.800415E7*RS(2,3))*RS(1,3)**3+(191.E0/0.160083E8*RS(2,2)-191.
     #E0/0.160083E8*RS(2,4))*RS(1,4)*RS(1,3)**2+(212.E0/0.4002075E7*RS(2
     #,4)+163.E0/0.53361E7*RS(2,3)+191.E0/0.160083E8*RS(2,2))*RS(1,4)**2
     #*RS(1,3)+(-254.E0/0.4002075E7*RS(2,3)-499.E0/686070.E0*RS(2,4)-137
     #.E0/0.160083E7*RS(2,2))*RS(1,4)**3)*RS(1,2)
      s6 = s7+(-881.E0/0.800415E7*RS(2,2)+881.E0/0.800415E7*RS(2,4))*RS(
     #1,3)**4+(137.E0/0.160083E7*RS(2,4)-881.E0/0.800415E7*RS(2,3)-76.E0
     #/0.2401245E7*RS(2,2))*RS(1,4)*RS(1,3)**3+(-137.E0/0.160083E7*RS(2,
     #3)+19.E0/148225.E0*RS(2,4)-149.E0/0.800415E7*RS(2,2))*RS(1,4)**2*R
     #S(1,3)**2+(-17.E0/190575.E0*RS(2,4)-19.E0/148225.E0*RS(2,3)+2.E0/1
     #90575.E0*RS(2,2))*RS(1,4)**3*RS(1,3)+(17.E0/190575.E0*RS(2,3)-RS(2
     #,4)/210+RS(2,1)/105+499.E0/686070.E0*RS(2,2))*RS(1,4)**4
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(-11357.E0/0.160083E8*RS(2,4)+37517.E0/0.480249E8*RS(2,3)-
     #RS(2,1)/210)*RS(1,2)**5+((-17.E0/190575.E0*RS(2,3)-37517.E0/0.4802
     #49E8*RS(2,2)-43.E0/254100.E0*RS(2,4))*RS(1,3)+(733.E0/0.160083E8*R
     #S(2,3)+11357.E0/0.160083E8*RS(2,2)+881.E0/0.800415E7*RS(2,4))*RS(1
     #,4))*RS(1,2)**4+((17.E0/190575.E0*RS(2,2)-173.E0/0.17787E7*RS(2,4)
     #-13.E0/23716.E0*RS(2,3))*RS(1,3)**2+(76.E0/0.2401245E7*RS(2,4)+254
     #.E0/0.4002075E7*RS(2,3)+494.E0/0.4002075E7*RS(2,2))*RS(1,4)*RS(1,3
     #)+(-881.E0/0.800415E7*RS(2,2)-2029.E0/0.960498E8*RS(2,3)-59.E0/0.2
     #13444E7*RS(2,4))*RS(1,4)**2)*RS(1,2)**3
      s3 = s1+((-173.E0/0.17787E7*RS(2,4)+13.E0/23716.E0*RS(2,2)-499.E0/
     #686070.E0*RS(2,3))*RS(1,3)**3+(541.E0/0.160083E8*RS(2,2)+149.E0/0.
     #800415E7*RS(2,4)+463.E0/0.53361E7*RS(2,3))*RS(1,4)*RS(1,3)**2+(-33
     #7.E0/0.320166E8*RS(2,2)+337.E0/0.320166E8*RS(2,4))*RS(1,4)**2*RS(1
     #,3)+(2029.E0/0.960498E8*RS(2,3)+59.E0/0.213444E7*RS(2,2)+881.E0/0.
     #800415E7*RS(2,4))*RS(1,4)**3)*RS(1,2)**2
      s2 = s3+((499.E0/686070.E0*RS(2,2)-43.E0/254100.E0*RS(2,4)-11357.E
     #0/0.160083E8*RS(2,3))*RS(1,3)**4+(-2.E0/190575.E0*RS(2,4)+2.E0/190
     #575.E0*RS(2,2))*RS(1,4)*RS(1,3)**3+(-149.E0/0.800415E7*RS(2,2)-541
     #.E0/0.160083E8*RS(2,4)-463.E0/0.53361E7*RS(2,3))*RS(1,4)**2*RS(1,3
     #)**2+(-494.E0/0.4002075E7*RS(2,4)-254.E0/0.4002075E7*RS(2,3)-76.E0
     #/0.2401245E7*RS(2,2))*RS(1,4)**3*RS(1,3)+(-11357.E0/0.160083E8*RS(
     #2,4)-733.E0/0.160083E8*RS(2,3)-881.E0/0.800415E7*RS(2,2))*RS(1,4)*
     #*4)*RS(1,2)+(11357.E0/0.160083E8*RS(2,2)-11357.E0/0.160083E8*RS(2,
     #4))*RS(1,3)**5
      t0 = s2+(43.E0/254100.E0*RS(2,2)-499.E0/686070.E0*RS(2,4)+11357.E0
     #/0.160083E8*RS(2,3))*RS(1,4)*RS(1,3)**4+(499.E0/686070.E0*RS(2,3)-
     #13.E0/23716.E0*RS(2,4)+173.E0/0.17787E7*RS(2,2))*RS(1,4)**2*RS(1,3
     #)**3+(13.E0/23716.E0*RS(2,3)-17.E0/190575.E0*RS(2,4)+173.E0/0.1778
     #7E7*RS(2,2))*RS(1,4)**3*RS(1,3)**2+(17.E0/190575.E0*RS(2,3)+37517.
     #E0/0.480249E8*RS(2,4)+43.E0/254100.E0*RS(2,2))*RS(1,4)**4*RS(1,3)+
     #(-37517.E0/0.480249E8*RS(2,3)+RS(2,1)/210+11357.E0/0.160083E8*RS(2
     #,2))*RS(1,4)**5
	ARRRR(1) = T0

      s3 = (11357.E0/0.160083E8*RS(2,3)-37517.E0/0.480249E8*RS(2,4)+RS(2
     #,2)/210)*RS(1,1)**5+((17.E0/190575.E0*RS(2,4)-RS(2,1)/210+RS(2,2)/
     #105+499.E0/686070.E0*RS(2,3))*RS(1,2)+(-733.E0/0.160083E8*RS(2,4)-
     #881.E0/0.800415E7*RS(2,3)-11357.E0/0.160083E8*RS(2,1))*RS(1,3)+(37
     #517.E0/0.480249E8*RS(2,1)+17.E0/190575.E0*RS(2,4)+43.E0/254100.E0*
     #RS(2,3))*RS(1,4))*RS(1,1)**4
      s4 = s3
      s6 = ((RS(2,2)/70+13.E0/23716.E0*RS(2,4)-RS(2,1)/105+13.E0/23716.E
     #0*RS(2,3))*RS(1,2)**2+((-137.E0/0.160083E7*RS(2,3)-254.E0/0.400207
     #5E7*RS(2,4)-499.E0/686070.E0*RS(2,1))*RS(1,3)+(-19.E0/148225.E0*RS
     #(2,4)-17.E0/190575.E0*RS(2,1)+2.E0/190575.E0*RS(2,3))*RS(1,4))*RS(
     #1,2)+(881.E0/0.800415E7*RS(2,1)+59.E0/0.213444E7*RS(2,3)+2029.E0/0
     #.960498E8*RS(2,4))*RS(1,3)**2+(-254.E0/0.4002075E7*RS(2,4)-494.E0/
     #0.4002075E7*RS(2,1)-76.E0/0.2401245E7*RS(2,3))*RS(1,4)*RS(1,3)+(13
     #.E0/23716.E0*RS(2,4)+173.E0/0.17787E7*RS(2,3)-17.E0/190575.E0*RS(2
     #,1))*RS(1,4)**2)*RS(1,1)**3
      s8 = (2.E0/105.E0*RS(2,2)+17.E0/190575.E0*RS(2,3)-RS(2,1)/70+499.E
     #0/686070.E0*RS(2,4))*RS(1,2)**3+((-463.E0/0.53361E7*RS(2,4)-13.E0/
     #23716.E0*RS(2,1)-19.E0/148225.E0*RS(2,3))*RS(1,3)+(541.E0/0.160083
     #E8*RS(2,3)-13.E0/23716.E0*RS(2,1)-137.E0/0.160083E7*RS(2,4))*RS(1,
     #4))*RS(1,2)**2+((137.E0/0.160083E7*RS(2,1)-137.E0/0.160083E7*RS(2,
     #3))*RS(1,3)**2+(212.E0/0.4002075E7*RS(2,1)+163.E0/0.53361E7*RS(2,4
     #)+191.E0/0.160083E8*RS(2,3))*RS(1,4)*RS(1,3)+(-149.E0/0.800415E7*R
     #S(2,3)+19.E0/148225.E0*RS(2,1)-137.E0/0.160083E7*RS(2,4))*RS(1,4)*
     #*2)*RS(1,2)+(-2029.E0/0.960498E8*RS(2,4)-881.E0/0.800415E7*RS(2,3)
     #-59.E0/0.213444E7*RS(2,1))*RS(1,3)**3+(-337.E0/0.320166E8*RS(2,3)+
     #337.E0/0.320166E8*RS(2,1))*RS(1,4)*RS(1,3)**2+(-541.E0/0.160083E8*
     #RS(2,1)-149.E0/0.800415E7*RS(2,3)-463.E0/0.53361E7*RS(2,4))*RS(1,4
     #)**2*RS(1,3)+(-13.E0/23716.E0*RS(2,1)+173.E0/0.17787E7*RS(2,3)+499
     #.E0/686070.E0*RS(2,4))*RS(1,4)**3
      s9 = RS(1,1)**2
      s7 = s8*s9
      s5 = s6+s7
      s2 = s4+s5
      s4 = s2
      s8 = (RS(2,2)/42+11357.E0/0.160083E8*RS(2,4)-2.E0/105.E0*RS(2,1)-3
     #7517.E0/0.480249E8*RS(2,3))*RS(1,2)**4+((17.E0/190575.E0*RS(2,3)-1
     #7.E0/190575.E0*RS(2,1))*RS(1,3)+(-499.E0/686070.E0*RS(2,1)-881.E0/
     #0.800415E7*RS(2,4)+494.E0/0.4002075E7*RS(2,3))*RS(1,4))*RS(1,2)**3
      s7 = s8+((463.E0/0.53361E7*RS(2,4)+19.E0/148225.E0*RS(2,1)+13.E0/2
     #3716.E0*RS(2,3))*RS(1,3)**2+(-212.E0/0.4002075E7*RS(2,3)+212.E0/0.
     #4002075E7*RS(2,1))*RS(1,4)*RS(1,3)+(-337.E0/0.320166E8*RS(2,3)+137
     #.E0/0.160083E7*RS(2,1)+59.E0/0.213444E7*RS(2,4))*RS(1,4)**2)*RS(1,
     #2)**2+((254.E0/0.4002075E7*RS(2,4)+499.E0/686070.E0*RS(2,3)+137.E0
     #/0.160083E7*RS(2,1))*RS(1,3)**3+(-163.E0/0.53361E7*RS(2,4)-212.E0/
     #0.4002075E7*RS(2,3)-191.E0/0.160083E8*RS(2,1))*RS(1,4)*RS(1,3)**2+
     #(-191.E0/0.160083E8*RS(2,1)+191.E0/0.160083E8*RS(2,3))*RS(1,4)**2*
     #RS(1,3)+(137.E0/0.160083E7*RS(2,1)-76.E0/0.2401245E7*RS(2,3)-881.E
     #0/0.800415E7*RS(2,4))*RS(1,4)**3)*RS(1,2)
      s6 = s7+(881.E0/0.800415E7*RS(2,1)+11357.E0/0.160083E8*RS(2,3)+733
     #.E0/0.160083E8*RS(2,4))*RS(1,3)**4+(494.E0/0.4002075E7*RS(2,3)+254
     #.E0/0.4002075E7*RS(2,4)+76.E0/0.2401245E7*RS(2,1))*RS(1,4)*RS(1,3)
     #**3+(463.E0/0.53361E7*RS(2,4)+541.E0/0.160083E8*RS(2,3)+149.E0/0.8
     #00415E7*RS(2,1))*RS(1,4)**2*RS(1,3)**2+(2.E0/190575.E0*RS(2,3)-2.E
     #0/190575.E0*RS(2,1))*RS(1,4)**3*RS(1,3)+(43.E0/254100.E0*RS(2,3)-4
     #99.E0/686070.E0*RS(2,1)+11357.E0/0.160083E8*RS(2,4))*RS(1,4)**4
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(RS(2,3)/42-RS(2,1)/42)*RS(1,2)**5+((-11357.E0/0.160083E8*
     #RS(2,4)+37517.E0/0.480249E8*RS(2,1)+2.E0/105.E0*RS(2,3)-RS(2,2)/42
     #)*RS(1,3)+(-11357.E0/0.160083E8*RS(2,1)+11357.E0/0.160083E8*RS(2,3
     #))*RS(1,4))*RS(1,2)**4+((-17.E0/190575.E0*RS(2,1)-2.E0/105.E0*RS(2
     #,2)+RS(2,3)/70-499.E0/686070.E0*RS(2,4))*RS(1,3)**2+(-494.E0/0.400
     #2075E7*RS(2,1)+881.E0/0.800415E7*RS(2,4)+499.E0/686070.E0*RS(2,3))
     #*RS(1,4)*RS(1,3)+(881.E0/0.800415E7*RS(2,1)-881.E0/0.800415E7*RS(2
     #,3))*RS(1,4)**2)*RS(1,2)**3
      s3 = s1+((-13.E0/23716.E0*RS(2,1)-13.E0/23716.E0*RS(2,4)+RS(2,3)/1
     #05-RS(2,2)/70)*RS(1,3)**3+(137.E0/0.160083E7*RS(2,4)+13.E0/23716.E
     #0*RS(2,3)-541.E0/0.160083E8*RS(2,1))*RS(1,4)*RS(1,3)**2+(337.E0/0.
     #320166E8*RS(2,1)-59.E0/0.213444E7*RS(2,4)-137.E0/0.160083E7*RS(2,3
     #))*RS(1,4)**2*RS(1,3)+(59.E0/0.213444E7*RS(2,3)-59.E0/0.213444E7*R
     #S(2,1))*RS(1,4)**3)*RS(1,2)**2
      s2 = s3+((-499.E0/686070.E0*RS(2,1)+RS(2,3)/210-RS(2,2)/105-17.E0/
     #190575.E0*RS(2,4))*RS(1,3)**4+(-2.E0/190575.E0*RS(2,1)+17.E0/19057
     #5.E0*RS(2,3)+19.E0/148225.E0*RS(2,4))*RS(1,4)*RS(1,3)**3+(-19.E0/1
     #48225.E0*RS(2,3)+137.E0/0.160083E7*RS(2,4)+149.E0/0.800415E7*RS(2,
     #1))*RS(1,4)**2*RS(1,3)**2+(76.E0/0.2401245E7*RS(2,1)-137.E0/0.1600
     #83E7*RS(2,3)+881.E0/0.800415E7*RS(2,4))*RS(1,4)**3*RS(1,3)+(881.E0
     #/0.800415E7*RS(2,1)-881.E0/0.800415E7*RS(2,3))*RS(1,4)**4)*RS(1,2)
     #+(-11357.E0/0.160083E8*RS(2,1)+37517.E0/0.480249E8*RS(2,4)-RS(2,2)
     #/210)*RS(1,3)**5
      t0 = s2+(-17.E0/190575.E0*RS(2,4)-43.E0/254100.E0*RS(2,1)-37517.E0
     #/0.480249E8*RS(2,3))*RS(1,4)*RS(1,3)**4+(-173.E0/0.17787E7*RS(2,1)
     #-13.E0/23716.E0*RS(2,4)+17.E0/190575.E0*RS(2,3))*RS(1,4)**2*RS(1,3
     #)**3+(-173.E0/0.17787E7*RS(2,1)+13.E0/23716.E0*RS(2,3)-499.E0/6860
     #70.E0*RS(2,4))*RS(1,4)**3*RS(1,3)**2+(-43.E0/254100.E0*RS(2,1)-113
     #57.E0/0.160083E8*RS(2,4)+499.E0/686070.E0*RS(2,3))*RS(1,4)**4*RS(1
     #,3)+(-11357.E0/0.160083E8*RS(2,1)+11357.E0/0.160083E8*RS(2,3))*RS(
     #1,4)**5
	ARRRR(2) = T0

      s3 = (-11357.E0/0.160083E8*RS(2,2)+11357.E0/0.160083E8*RS(2,4))*RS
     #(1,1)**5+((43.E0/254100.E0*RS(2,4)-499.E0/686070.E0*RS(2,2)+11357.
     #E0/0.160083E8*RS(2,1))*RS(1,2)+(-881.E0/0.800415E7*RS(2,4)+881.E0/
     #0.800415E7*RS(2,2))*RS(1,3)+(499.E0/686070.E0*RS(2,4)-43.E0/254100
     #.E0*RS(2,2)-11357.E0/0.160083E8*RS(2,1))*RS(1,4))*RS(1,1)**4
      s4 = s3
      s6 = ((173.E0/0.17787E7*RS(2,4)-13.E0/23716.E0*RS(2,2)+499.E0/6860
     #70.E0*RS(2,1))*RS(1,2)**2+((-881.E0/0.800415E7*RS(2,1)+137.E0/0.16
     #0083E7*RS(2,2)-76.E0/0.2401245E7*RS(2,4))*RS(1,3)+(2.E0/190575.E0*
     #RS(2,4)-2.E0/190575.E0*RS(2,2))*RS(1,4))*RS(1,2)+(-59.E0/0.213444E
     #7*RS(2,2)+59.E0/0.213444E7*RS(2,4))*RS(1,3)**2+(76.E0/0.2401245E7*
     #RS(2,2)-137.E0/0.160083E7*RS(2,4)+881.E0/0.800415E7*RS(2,1))*RS(1,
     #4)*RS(1,3)+(13.E0/23716.E0*RS(2,4)-173.E0/0.17787E7*RS(2,2)-499.E0
     #/686070.E0*RS(2,1))*RS(1,4)**2)*RS(1,1)**3
      s8 = (173.E0/0.17787E7*RS(2,4)+13.E0/23716.E0*RS(2,1)-17.E0/190575
     #.E0*RS(2,2))*RS(1,2)**3+((19.E0/148225.E0*RS(2,2)-137.E0/0.160083E
     #7*RS(2,1)-149.E0/0.800415E7*RS(2,4))*RS(1,3)+(-541.E0/0.160083E8*R
     #S(2,2)-149.E0/0.800415E7*RS(2,4)-463.E0/0.53361E7*RS(2,1))*RS(1,4)
     #)*RS(1,2)**2+((-337.E0/0.320166E8*RS(2,4)+137.E0/0.160083E7*RS(2,2
     #)+59.E0/0.213444E7*RS(2,1))*RS(1,3)**2+(191.E0/0.160083E8*RS(2,4)-
     #191.E0/0.160083E8*RS(2,2))*RS(1,4)*RS(1,3)+(149.E0/0.800415E7*RS(2
     #,2)+463.E0/0.53361E7*RS(2,1)+541.E0/0.160083E8*RS(2,4))*RS(1,4)**2
     #)*RS(1,2)+(-881.E0/0.800415E7*RS(2,4)+881.E0/0.800415E7*RS(2,2))*R
     #S(1,3)**3+(-59.E0/0.213444E7*RS(2,1)-137.E0/0.160083E7*RS(2,4)+337
     #.E0/0.320166E8*RS(2,2))*RS(1,4)*RS(1,3)**2+(-19.E0/148225.E0*RS(2,
     #4)+149.E0/0.800415E7*RS(2,2)+137.E0/0.160083E7*RS(2,1))*RS(1,4)**2
     #*RS(1,3)+(-173.E0/0.17787E7*RS(2,2)-13.E0/23716.E0*RS(2,1)+17.E0/1
     #90575.E0*RS(2,4))*RS(1,4)**3
      s9 = RS(1,1)**2
      s7 = s8*s9
      s5 = s6+s7
      s2 = s4+s5
      s4 = s2
      s8 = (37517.E0/0.480249E8*RS(2,2)+43.E0/254100.E0*RS(2,4)+17.E0/19
     #0575.E0*RS(2,1))*RS(1,2)**4+((-17.E0/190575.E0*RS(2,2)+2.E0/190575
     #.E0*RS(2,4)-19.E0/148225.E0*RS(2,1))*RS(1,3)+(-254.E0/0.4002075E7*
     #RS(2,1)-76.E0/0.2401245E7*RS(2,4)-494.E0/0.4002075E7*RS(2,2))*RS(1
     #,4))*RS(1,2)**3
      s7 = s8+((541.E0/0.160083E8*RS(2,4)-137.E0/0.160083E7*RS(2,1)-13.E
     #0/23716.E0*RS(2,2))*RS(1,3)**2+(212.E0/0.4002075E7*RS(2,2)+163.E0/
     #0.53361E7*RS(2,1)+191.E0/0.160083E8*RS(2,4))*RS(1,4)*RS(1,3)+(337.
     #E0/0.320166E8*RS(2,2)-337.E0/0.320166E8*RS(2,4))*RS(1,4)**2)*RS(1,
     #2)**2+((494.E0/0.4002075E7*RS(2,4)-499.E0/686070.E0*RS(2,2)-881.E0
     #/0.800415E7*RS(2,1))*RS(1,3)**3+(212.E0/0.4002075E7*RS(2,2)-212.E0
     #/0.4002075E7*RS(2,4))*RS(1,4)*RS(1,3)**2+(-212.E0/0.4002075E7*RS(2
     #,4)-191.E0/0.160083E8*RS(2,2)-163.E0/0.53361E7*RS(2,1))*RS(1,4)**2
     #*RS(1,3)+(76.E0/0.2401245E7*RS(2,2)+494.E0/0.4002075E7*RS(2,4)+254
     #.E0/0.4002075E7*RS(2,1))*RS(1,4)**3)*RS(1,2)
      s6 = s7+(-11357.E0/0.160083E8*RS(2,2)+11357.E0/0.160083E8*RS(2,4))
     #*RS(1,3)**4+(-494.E0/0.4002075E7*RS(2,2)+499.E0/686070.E0*RS(2,4)+
     #881.E0/0.800415E7*RS(2,1))*RS(1,4)*RS(1,3)**3+(137.E0/0.160083E7*R
     #S(2,1)+13.E0/23716.E0*RS(2,4)-541.E0/0.160083E8*RS(2,2))*RS(1,4)**
     #2*RS(1,3)**2+(-2.E0/190575.E0*RS(2,2)+19.E0/148225.E0*RS(2,1)+17.E
     #0/190575.E0*RS(2,4))*RS(1,4)**3*RS(1,3)+(-43.E0/254100.E0*RS(2,2)-
     #37517.E0/0.480249E8*RS(2,4)-17.E0/190575.E0*RS(2,1))*RS(1,4)**4
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(11357.E0/0.160083E8*RS(2,4)+RS(2,3)/210-37517.E0/0.480249
     #E8*RS(2,1))*RS(1,2)**5+((499.E0/686070.E0*RS(2,4)+RS(2,3)/105+17.E
     #0/190575.E0*RS(2,1)-RS(2,2)/210)*RS(1,3)+(-733.E0/0.160083E8*RS(2,
     #1)-881.E0/0.800415E7*RS(2,4)-11357.E0/0.160083E8*RS(2,2))*RS(1,4))
     #*RS(1,2)**4+((RS(2,3)/70-RS(2,2)/105+13.E0/23716.E0*RS(2,1)+13.E0/
     #23716.E0*RS(2,4))*RS(1,3)**2+(-254.E0/0.4002075E7*RS(2,1)-499.E0/6
     #86070.E0*RS(2,2)-137.E0/0.160083E7*RS(2,4))*RS(1,4)*RS(1,3)+(881.E
     #0/0.800415E7*RS(2,2)+2029.E0/0.960498E8*RS(2,1)+59.E0/0.213444E7*R
     #S(2,4))*RS(1,4)**2)*RS(1,2)**3
      s3 = s1+((17.E0/190575.E0*RS(2,4)-RS(2,2)/70+2.E0/105.E0*RS(2,3)+4
     #99.E0/686070.E0*RS(2,1))*RS(1,3)**3+(-463.E0/0.53361E7*RS(2,1)-19.
     #E0/148225.E0*RS(2,4)-13.E0/23716.E0*RS(2,2))*RS(1,4)*RS(1,3)**2+(-
     #137.E0/0.160083E7*RS(2,4)+137.E0/0.160083E7*RS(2,2))*RS(1,4)**2*RS
     #(1,3)+(-59.E0/0.213444E7*RS(2,2)-2029.E0/0.960498E8*RS(2,1)-881.E0
     #/0.800415E7*RS(2,4))*RS(1,4)**3)*RS(1,2)**2
      s2 = s3+((-2.E0/105.E0*RS(2,2)-37517.E0/0.480249E8*RS(2,4)+11357.E
     #0/0.160083E8*RS(2,1)+RS(2,3)/42)*RS(1,3)**4+(17.E0/190575.E0*RS(2,
     #4)-17.E0/190575.E0*RS(2,2))*RS(1,4)*RS(1,3)**3+(13.E0/23716.E0*RS(
     #2,4)+19.E0/148225.E0*RS(2,2)+463.E0/0.53361E7*RS(2,1))*RS(1,4)**2*
     #RS(1,3)**2+(499.E0/686070.E0*RS(2,4)+254.E0/0.4002075E7*RS(2,1)+13
     #7.E0/0.160083E7*RS(2,2))*RS(1,4)**3*RS(1,3)+(881.E0/0.800415E7*RS(
     #2,2)+11357.E0/0.160083E8*RS(2,4)+733.E0/0.160083E8*RS(2,1))*RS(1,4
     #)**4)*RS(1,2)+(RS(2,4)/42-RS(2,2)/42)*RS(1,3)**5
      t0 = s2+(-11357.E0/0.160083E8*RS(2,1)+2.E0/105.E0*RS(2,4)-RS(2,3)/
     #42+37517.E0/0.480249E8*RS(2,2))*RS(1,4)*RS(1,3)**4+(-17.E0/190575.
     #E0*RS(2,2)-2.E0/105.E0*RS(2,3)-499.E0/686070.E0*RS(2,1)+RS(2,4)/70
     #)*RS(1,4)**2*RS(1,3)**3+(RS(2,4)/105-13.E0/23716.E0*RS(2,2)-13.E0/
     #23716.E0*RS(2,1)-RS(2,3)/70)*RS(1,4)**3*RS(1,3)**2+(-499.E0/686070
     #.E0*RS(2,2)-RS(2,3)/105-17.E0/190575.E0*RS(2,1)+RS(2,4)/210)*RS(1,
     #4)**4*RS(1,3)+(-11357.E0/0.160083E8*RS(2,2)-RS(2,3)/210+37517.E0/0
     #.480249E8*RS(2,1))*RS(1,4)**5
	ARRRR(3) = T0

      s3 = (37517.E0/0.480249E8*RS(2,2)-RS(2,4)/210-11357.E0/0.160083E8*
     #RS(2,3))*RS(1,1)**5+((-17.E0/190575.E0*RS(2,2)-43.E0/254100.E0*RS(
     #2,3)-37517.E0/0.480249E8*RS(2,1))*RS(1,2)+(733.E0/0.160083E8*RS(2,
     #2)+11357.E0/0.160083E8*RS(2,1)+881.E0/0.800415E7*RS(2,3))*RS(1,3)+
     #(-17.E0/190575.E0*RS(2,2)+RS(2,1)/210-499.E0/686070.E0*RS(2,3)-RS(
     #2,4)/105)*RS(1,4))*RS(1,1)**4
      s4 = s3
      s6 = ((-173.E0/0.17787E7*RS(2,3)+17.E0/190575.E0*RS(2,1)-13.E0/237
     #16.E0*RS(2,2))*RS(1,2)**2+((76.E0/0.2401245E7*RS(2,3)+494.E0/0.400
     #2075E7*RS(2,1)+254.E0/0.4002075E7*RS(2,2))*RS(1,3)+(-2.E0/190575.E
     #0*RS(2,3)+19.E0/148225.E0*RS(2,2)+17.E0/190575.E0*RS(2,1))*RS(1,4)
     #)*RS(1,2)+(-881.E0/0.800415E7*RS(2,1)-2029.E0/0.960498E8*RS(2,2)-5
     #9.E0/0.213444E7*RS(2,3))*RS(1,3)**2+(137.E0/0.160083E7*RS(2,3)+254
     #.E0/0.4002075E7*RS(2,2)+499.E0/686070.E0*RS(2,1))*RS(1,4)*RS(1,3)+
     #(-13.E0/23716.E0*RS(2,3)+RS(2,1)/105-13.E0/23716.E0*RS(2,2)-RS(2,4
     #)/70)*RS(1,4)**2)*RS(1,1)**3
      s8 = (-499.E0/686070.E0*RS(2,2)-173.E0/0.17787E7*RS(2,3)+13.E0/237
     #16.E0*RS(2,1))*RS(1,2)**3+((149.E0/0.800415E7*RS(2,3)+463.E0/0.533
     #61E7*RS(2,2)+541.E0/0.160083E8*RS(2,1))*RS(1,3)+(137.E0/0.160083E7
     #*RS(2,2)+149.E0/0.800415E7*RS(2,3)-19.E0/148225.E0*RS(2,1))*RS(1,4
     #))*RS(1,2)**2+((-337.E0/0.320166E8*RS(2,1)+337.E0/0.320166E8*RS(2,
     #3))*RS(1,3)**2+(-191.E0/0.160083E8*RS(2,3)-212.E0/0.4002075E7*RS(2
     #,1)-163.E0/0.53361E7*RS(2,2))*RS(1,4)*RS(1,3)+(137.E0/0.160083E7*R
     #S(2,2)-541.E0/0.160083E8*RS(2,3)+13.E0/23716.E0*RS(2,1))*RS(1,4)**
     #2)*RS(1,2)+(881.E0/0.800415E7*RS(2,3)+2029.E0/0.960498E8*RS(2,2)+5
     #9.E0/0.213444E7*RS(2,1))*RS(1,3)**3+(-137.E0/0.160083E7*RS(2,1)+13
     #7.E0/0.160083E7*RS(2,3))*RS(1,4)*RS(1,3)**2+(19.E0/148225.E0*RS(2,
     #3)+13.E0/23716.E0*RS(2,1)+463.E0/0.53361E7*RS(2,2))*RS(1,4)**2*RS(
     #1,3)+(RS(2,1)/70-499.E0/686070.E0*RS(2,2)-17.E0/190575.E0*RS(2,3)-
     #2.E0/105.E0*RS(2,4))*RS(1,4)**3
      s9 = RS(1,1)**2
      s7 = s8*s9
      s5 = s6+s7
      s2 = s4+s5
      s4 = s2
      s7 = (-11357.E0/0.160083E8*RS(2,2)+499.E0/686070.E0*RS(2,1)-43.E0/
     #254100.E0*RS(2,3))*RS(1,2)**4+((2.E0/190575.E0*RS(2,1)-2.E0/190575
     #.E0*RS(2,3))*RS(1,3)+(76.E0/0.2401245E7*RS(2,3)-137.E0/0.160083E7*
     #RS(2,1)+881.E0/0.800415E7*RS(2,2))*RS(1,4))*RS(1,2)**3+((-463.E0/0
     #.53361E7*RS(2,2)-541.E0/0.160083E8*RS(2,3)-149.E0/0.800415E7*RS(2,
     #1))*RS(1,3)**2+(-191.E0/0.160083E8*RS(2,3)+191.E0/0.160083E8*RS(2,
     #1))*RS(1,4)*RS(1,3)+(-137.E0/0.160083E7*RS(2,1)+337.E0/0.320166E8*
     #RS(2,3)-59.E0/0.213444E7*RS(2,2))*RS(1,4)**2)*RS(1,2)**2+((-254.E0
     #/0.4002075E7*RS(2,2)-76.E0/0.2401245E7*RS(2,1)-494.E0/0.4002075E7*
     #RS(2,3))*RS(1,3)**3+(212.E0/0.4002075E7*RS(2,3)+191.E0/0.160083E8*
     #RS(2,1)+163.E0/0.53361E7*RS(2,2))*RS(1,4)*RS(1,3)**2+(212.E0/0.400
     #2075E7*RS(2,3)-212.E0/0.4002075E7*RS(2,1))*RS(1,4)**2*RS(1,3)+(499
     #.E0/686070.E0*RS(2,1)+881.E0/0.800415E7*RS(2,2)-494.E0/0.4002075E7
     #*RS(2,3))*RS(1,4)**3)*RS(1,2)
      s6 = s7+(-11357.E0/0.160083E8*RS(2,3)-733.E0/0.160083E8*RS(2,2)-88
     #1.E0/0.800415E7*RS(2,1))*RS(1,3)**4+(-137.E0/0.160083E7*RS(2,1)-49
     #9.E0/686070.E0*RS(2,3)-254.E0/0.4002075E7*RS(2,2))*RS(1,4)*RS(1,3)
     #**3+(-19.E0/148225.E0*RS(2,1)-13.E0/23716.E0*RS(2,3)-463.E0/0.5336
     #1E7*RS(2,2))*RS(1,4)**2*RS(1,3)**2+(-17.E0/190575.E0*RS(2,3)+17.E0
     #/190575.E0*RS(2,1))*RS(1,4)**3*RS(1,3)+(37517.E0/0.480249E8*RS(2,3
     #)-11357.E0/0.160083E8*RS(2,2)+2.E0/105.E0*RS(2,1)-RS(2,4)/42)*RS(1
     #,4)**4
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(-11357.E0/0.160083E8*RS(2,3)+11357.E0/0.160083E8*RS(2,1))
     #*RS(1,2)**5+((-499.E0/686070.E0*RS(2,3)+43.E0/254100.E0*RS(2,1)+11
     #357.E0/0.160083E8*RS(2,2))*RS(1,3)+(-881.E0/0.800415E7*RS(2,1)+881
     #.E0/0.800415E7*RS(2,3))*RS(1,4))*RS(1,2)**4+((499.E0/686070.E0*RS(
     #2,2)+173.E0/0.17787E7*RS(2,1)-13.E0/23716.E0*RS(2,3))*RS(1,3)**2+(
     #-881.E0/0.800415E7*RS(2,2)+137.E0/0.160083E7*RS(2,3)-76.E0/0.24012
     #45E7*RS(2,1))*RS(1,4)*RS(1,3)+(59.E0/0.213444E7*RS(2,1)-59.E0/0.21
     #3444E7*RS(2,3))*RS(1,4)**2)*RS(1,2)**3
      s3 = s1+((173.E0/0.17787E7*RS(2,1)+13.E0/23716.E0*RS(2,2)-17.E0/19
     #0575.E0*RS(2,3))*RS(1,3)**3+(-137.E0/0.160083E7*RS(2,2)+19.E0/1482
     #25.E0*RS(2,3)-149.E0/0.800415E7*RS(2,1))*RS(1,4)*RS(1,3)**2+(59.E0
     #/0.213444E7*RS(2,2)-337.E0/0.320166E8*RS(2,1)+137.E0/0.160083E7*RS
     #(2,3))*RS(1,4)**2*RS(1,3)+(-881.E0/0.800415E7*RS(2,1)+881.E0/0.800
     #415E7*RS(2,3))*RS(1,4)**3)*RS(1,2)**2
      s2 = s3+((37517.E0/0.480249E8*RS(2,3)+43.E0/254100.E0*RS(2,1)+17.E
     #0/190575.E0*RS(2,2))*RS(1,3)**4+(-17.E0/190575.E0*RS(2,3)-19.E0/14
     #8225.E0*RS(2,2)+2.E0/190575.E0*RS(2,1))*RS(1,4)*RS(1,3)**3+(541.E0
     #/0.160083E8*RS(2,1)-137.E0/0.160083E7*RS(2,2)-13.E0/23716.E0*RS(2,
     #3))*RS(1,4)**2*RS(1,3)**2+(-881.E0/0.800415E7*RS(2,2)+494.E0/0.400
     #2075E7*RS(2,1)-499.E0/686070.E0*RS(2,3))*RS(1,4)**3*RS(1,3)+(-1135
     #7.E0/0.160083E8*RS(2,3)+11357.E0/0.160083E8*RS(2,1))*RS(1,4)**4)*R
     #S(1,2)+(RS(2,4)/210+11357.E0/0.160083E8*RS(2,1)-37517.E0/0.480249E
     #8*RS(2,2))*RS(1,3)**5
      t0 = s2+(RS(2,4)/105+17.E0/190575.E0*RS(2,2)-RS(2,3)/210+499.E0/68
     #6070.E0*RS(2,1))*RS(1,4)*RS(1,3)**4+(13.E0/23716.E0*RS(2,1)-RS(2,3
     #)/105+RS(2,4)/70+13.E0/23716.E0*RS(2,2))*RS(1,4)**2*RS(1,3)**3+(49
     #9.E0/686070.E0*RS(2,2)+2.E0/105.E0*RS(2,4)+17.E0/190575.E0*RS(2,1)
     #-RS(2,3)/70)*RS(1,4)**3*RS(1,3)**2+(-2.E0/105.E0*RS(2,3)+RS(2,4)/4
     #2-37517.E0/0.480249E8*RS(2,1)+11357.E0/0.160083E8*RS(2,2))*RS(1,4)
     #**4*RS(1,3)+(RS(2,1)/42-RS(2,3)/42)*RS(1,4)**5
	ARRRR(4) = T0

      s3 = (-56311.E0/0.1200622E8*RS(2,4)+56311.E0/0.1200622E8*RS(2,2))*
     #RS(1,1)**5+((1817.E0/0.1200622E8*RS(2,3)-56311.E0/0.1200622E8*RS(2
     #,1)-12482.E0/0.1200622E8*RS(2,4)+9568.E0/0.1715175E7*RS(2,2))*RS(1
     #,2)+(-3463.E0/0.4002075E7*RS(2,2)+3463.E0/0.4002075E7*RS(2,4))*RS(
     #1,3)+(12482.E0/0.1200622E8*RS(2,2)-9568.E0/0.1715175E7*RS(2,4)-181
     #7.E0/0.1200622E8*RS(2,3)+56311.E0/0.1200622E8*RS(2,1))*RS(1,4))*RS
     #(1,1)**4
      s4 = s3
      s6 = ((521.E0/88935.E0*RS(2,2)-106.E0/190575.E0*RS(2,4)-9568.E0/0.
     #1715175E7*RS(2,1)+3319.E0/0.1200622E8*RS(2,3))*RS(1,2)**2+((2624.E
     #0/0.1200622E8*RS(2,4)+8572.E0/0.1200622E8*RS(2,1)-10522.E0/0.12006
     #22E8*RS(2,2)-674.E0/0.1200622E8*RS(2,3))*RS(1,3)+(-92.E0/444675.E0
     #*RS(2,2)+92.E0/444675.E0*RS(2,4))*RS(1,4))*RS(1,2)+(1907.E0/0.1200
     #622E8*RS(2,2)-1907.E0/0.1200622E8*RS(2,4))*RS(1,3)**2+(10522.E0/0.
     #1200622E8*RS(2,4)-2624.E0/0.1200622E8*RS(2,2)-8572.E0/0.1200622E8*
     #RS(2,1)+674.E0/0.1200622E8*RS(2,3))*RS(1,4)*RS(1,3)+(9568.E0/0.171
     #5175E7*RS(2,1)-3319.E0/0.1200622E8*RS(2,3)-521.E0/88935.E0*RS(2,4)
     #+106.E0/190575.E0*RS(2,2))*RS(1,4)**2)*RS(1,1)**3
      s9 = (106.E0/190575.E0*RS(2,3)-3319.E0/0.1200622E8*RS(2,4)-521.E0/
     #88935.E0*RS(2,1)+9568.E0/0.1715175E7*RS(2,2))*RS(1,2)**3+((557.E0/
     #0.4002075E7*RS(2,4)+49.E0/81675.E0*RS(2,1)+32.E0/0.1334025E7*RS(2,
     #3)-1018.E0/0.1334025E7*RS(2,2))*RS(1,3)+(-32.E0/0.1334025E7*RS(2,4
     #)-49.E0/81675.E0*RS(2,2)-557.E0/0.4002075E7*RS(2,3)+1018.E0/0.1334
     #025E7*RS(2,1))*RS(1,4))*RS(1,2)**2+((-137.E0/0.1334025E7*RS(2,1)+1
     #37.E0/0.1334025E7*RS(2,3))*RS(1,3)**2+(382.E0/0.4002075E7*RS(2,2)-
     #382.E0/0.4002075E7*RS(2,4))*RS(1,4)*RS(1,3)+(-1018.E0/0.1334025E7*
     #RS(2,1)+32.E0/0.1334025E7*RS(2,2)+557.E0/0.4002075E7*RS(2,3)+49.E0
     #/81675.E0*RS(2,4))*RS(1,4)**2)*RS(1,2)
      s8 = s9+(1907.E0/0.1200622E8*RS(2,4)-1907.E0/0.1200622E8*RS(2,2))*
     #RS(1,3)**3+(-137.E0/0.1334025E7*RS(2,3)+137.E0/0.1334025E7*RS(2,1)
     #)*RS(1,4)*RS(1,3)**2+(-557.E0/0.4002075E7*RS(2,2)-49.E0/81675.E0*R
     #S(2,1)+1018.E0/0.1334025E7*RS(2,4)-32.E0/0.1334025E7*RS(2,3))*RS(1
     #,4)**2*RS(1,3)+(-106.E0/190575.E0*RS(2,3)+521.E0/88935.E0*RS(2,1)-
     #9568.E0/0.1715175E7*RS(2,4)+3319.E0/0.1200622E8*RS(2,2))*RS(1,4)**
     #3
      s9 = RS(1,1)**2
      s7 = s8*s9
      s5 = s6+s7
      s2 = s4+s5
      s4 = s2
      s8 = (-1817.E0/0.1200622E8*RS(2,4)-9568.E0/0.1715175E7*RS(2,1)+563
     #11.E0/0.1200622E8*RS(2,2)+12482.E0/0.1200622E8*RS(2,3))*RS(1,2)**4
     #+((-92.E0/444675.E0*RS(2,3)+92.E0/444675.E0*RS(2,1))*RS(1,3)+(674.
     #E0/0.1200622E8*RS(2,4)+10522.E0/0.1200622E8*RS(2,1)-2624.E0/0.1200
     #622E8*RS(2,3)-8572.E0/0.1200622E8*RS(2,2))*RS(1,4))*RS(1,2)**3
      s7 = s8+((-32.E0/0.1334025E7*RS(2,1)-557.E0/0.4002075E7*RS(2,4)+10
     #18.E0/0.1334025E7*RS(2,2)-49.E0/81675.E0*RS(2,3))*RS(1,3)**2+(382.
     #E0/0.4002075E7*RS(2,3)-382.E0/0.4002075E7*RS(2,1))*RS(1,4)*RS(1,3)
     #+(-137.E0/0.1334025E7*RS(2,4)+137.E0/0.1334025E7*RS(2,2))*RS(1,4)*
     #*2)*RS(1,2)**2+((674.E0/0.1200622E8*RS(2,1)-2624.E0/0.1200622E8*RS
     #(2,4)-8572.E0/0.1200622E8*RS(2,3)+10522.E0/0.1200622E8*RS(2,2))*RS
     #(1,3)**3+(382.E0/0.4002075E7*RS(2,4)-382.E0/0.4002075E7*RS(2,2))*R
     #S(1,4)*RS(1,3)**2+(382.E0/0.4002075E7*RS(2,1)-382.E0/0.4002075E7*R
     #S(2,3))*RS(1,4)**2*RS(1,3)+(-10522.E0/0.1200622E8*RS(2,1)+2624.E0/
     #0.1200622E8*RS(2,3)-674.E0/0.1200622E8*RS(2,2)+8572.E0/0.1200622E8
     #*RS(2,4))*RS(1,4)**3)*RS(1,2)
      s6 = s7+(3463.E0/0.4002075E7*RS(2,2)-3463.E0/0.4002075E7*RS(2,4))*
     #RS(1,3)**4+(-10522.E0/0.1200622E8*RS(2,4)-674.E0/0.1200622E8*RS(2,
     #1)+8572.E0/0.1200622E8*RS(2,3)+2624.E0/0.1200622E8*RS(2,2))*RS(1,4
     #)*RS(1,3)**3+(32.E0/0.1334025E7*RS(2,1)+49.E0/81675.E0*RS(2,3)-101
     #8.E0/0.1334025E7*RS(2,4)+557.E0/0.4002075E7*RS(2,2))*RS(1,4)**2*RS
     #(1,3)**2+(92.E0/444675.E0*RS(2,3)-92.E0/444675.E0*RS(2,1))*RS(1,4)
     #**3*RS(1,3)+(1817.E0/0.1200622E8*RS(2,2)-12482.E0/0.1200622E8*RS(2
     #,3)-56311.E0/0.1200622E8*RS(2,4)+9568.E0/0.1715175E7*RS(2,1))*RS(1
     #,4)**4
      s7 = RS(1,1)
      s5 = s6*s7
      s3 = s4+s5
      s1 = s3+(56311.E0/0.1200622E8*RS(2,3)-56311.E0/0.1200622E8*RS(2,1)
     #)*RS(1,2)**5+((-12482.E0/0.1200622E8*RS(2,1)+9568.E0/0.1715175E7*R
     #S(2,3)+1817.E0/0.1200622E8*RS(2,4)-56311.E0/0.1200622E8*RS(2,2))*R
     #S(1,3)+(3463.E0/0.4002075E7*RS(2,1)-3463.E0/0.4002075E7*RS(2,3))*R
     #S(1,4))*RS(1,2)**4+((-106.E0/190575.E0*RS(2,1)+521.E0/88935.E0*RS(
     #2,3)+3319.E0/0.1200622E8*RS(2,4)-9568.E0/0.1715175E7*RS(2,2))*RS(1
     #,3)**2+(2624.E0/0.1200622E8*RS(2,1)-674.E0/0.1200622E8*RS(2,4)+857
     #2.E0/0.1200622E8*RS(2,2)-10522.E0/0.1200622E8*RS(2,3))*RS(1,4)*RS(
     #1,3)+(-1907.E0/0.1200622E8*RS(2,1)+1907.E0/0.1200622E8*RS(2,3))*RS
     #(1,4)**2)*RS(1,2)**3
      s3 = s1+((9568.E0/0.1715175E7*RS(2,3)+106.E0/190575.E0*RS(2,4)-331
     #9.E0/0.1200622E8*RS(2,1)-521.E0/88935.E0*RS(2,2))*RS(1,3)**3+(32.E
     #0/0.1334025E7*RS(2,4)+557.E0/0.4002075E7*RS(2,1)+49.E0/81675.E0*RS
     #(2,2)-1018.E0/0.1334025E7*RS(2,3))*RS(1,4)*RS(1,3)**2+(137.E0/0.13
     #34025E7*RS(2,4)-137.E0/0.1334025E7*RS(2,2))*RS(1,4)**2*RS(1,3)+(-1
     #907.E0/0.1200622E8*RS(2,3)+1907.E0/0.1200622E8*RS(2,1))*RS(1,4)**3
     #)*RS(1,2)**2
      s2 = s3+((56311.E0/0.1200622E8*RS(2,3)-1817.E0/0.1200622E8*RS(2,1)
     #+12482.E0/0.1200622E8*RS(2,4)-9568.E0/0.1715175E7*RS(2,2))*RS(1,3)
     #**4+(-92.E0/444675.E0*RS(2,4)+92.E0/444675.E0*RS(2,2))*RS(1,4)*RS(
     #1,3)**3+(-32.E0/0.1334025E7*RS(2,2)+1018.E0/0.1334025E7*RS(2,3)-49
     #.E0/81675.E0*RS(2,4)-557.E0/0.4002075E7*RS(2,1))*RS(1,4)**2*RS(1,3
     #)**2+(-8572.E0/0.1200622E8*RS(2,4)+674.E0/0.1200622E8*RS(2,2)-2624
     #.E0/0.1200622E8*RS(2,1)+10522.E0/0.1200622E8*RS(2,3))*RS(1,4)**3*R
     #S(1,3)+(3463.E0/0.4002075E7*RS(2,3)-3463.E0/0.4002075E7*RS(2,1))*R
     #S(1,4)**4)*RS(1,2)+(56311.E0/0.1200622E8*RS(2,4)-56311.E0/0.120062
     #2E8*RS(2,2))*RS(1,3)**5
      t0 = s2+(1817.E0/0.1200622E8*RS(2,1)-12482.E0/0.1200622E8*RS(2,2)-
     #56311.E0/0.1200622E8*RS(2,3)+9568.E0/0.1715175E7*RS(2,4))*RS(1,4)*
     #RS(1,3)**4+(521.E0/88935.E0*RS(2,4)-9568.E0/0.1715175E7*RS(2,3)+33
     #19.E0/0.1200622E8*RS(2,1)-106.E0/190575.E0*RS(2,2))*RS(1,4)**2*RS(
     #1,3)**3+(-521.E0/88935.E0*RS(2,3)+9568.E0/0.1715175E7*RS(2,4)+106.
     #E0/190575.E0*RS(2,1)-3319.E0/0.1200622E8*RS(2,2))*RS(1,4)**3*RS(1,
     #3)**2+(56311.E0/0.1200622E8*RS(2,4)+12482.E0/0.1200622E8*RS(2,1)-1
     #817.E0/0.1200622E8*RS(2,2)-9568.E0/0.1715175E7*RS(2,3))*RS(1,4)**4
     #*RS(1,3)+(-56311.E0/0.1200622E8*RS(2,3)+56311.E0/0.1200622E8*RS(2,
     #1))*RS(1,4)**5
	ARRRR(5) = T0

C	------------------------
C	COMPUTE FOR FACTOR ASSSS
C	------------------------
      s6 = (-RS(1,2)/42+RS(1,4)/42)*RS(2,1)**5+((11357.E0/0.160083E8*RS(
     #1,3)-2.E0/105.E0*RS(1,2)+RS(1,1)/42-37517.E0/0.480249E8*RS(1,4))*R
     #S(2,2)+(-11357.E0/0.160083E8*RS(1,2)+11357.E0/0.160083E8*RS(1,4))*
     #RS(2,3)+(-11357.E0/0.160083E8*RS(1,3)+2.E0/105.E0*RS(1,4)+37517.E0
     #/0.480249E8*RS(1,2)-RS(1,1)/42)*RS(2,4))*RS(2,1)**4
      s7 = s6
      s9 = ((-RS(1,2)/70+499.E0/686070.E0*RS(1,3)+2.E0/105.E0*RS(1,1)+17
     #.E0/190575.E0*RS(1,4))*RS(2,2)**2+((-499.E0/686070.E0*RS(1,2)+494.
     #E0/0.4002075E7*RS(1,4)-881.E0/0.800415E7*RS(1,3))*RS(2,3)+(17.E0/1
     #90575.E0*RS(1,4)-17.E0/190575.E0*RS(1,2))*RS(2,4))*RS(2,2)+(881.E0
     #/0.800415E7*RS(1,2)-881.E0/0.800415E7*RS(1,4))*RS(2,3)**2+(-494.E0
     #/0.4002075E7*RS(1,2)+499.E0/686070.E0*RS(1,4)+881.E0/0.800415E7*RS
     #(1,3))*RS(2,4)*RS(2,3)+(-2.E0/105.E0*RS(1,1)-17.E0/190575.E0*RS(1,
     #2)+RS(1,4)/70-499.E0/686070.E0*RS(1,3))*RS(2,4)**2)*RS(2,1)**3
      s10 = ((43.E0/254100.E0*RS(1,4)-499.E0/686070.E0*RS(1,2)+11357.E0/
     #0.160083E8*RS(1,3))*RS(2,3)**4+(2.E0/190575.E0*RS(1,4)-2.E0/190575
     #.E0*RS(1,2))*RS(2,4)*RS(2,3)**3+(463.E0/0.53361E7*RS(1,3)+149.E0/0
     #.800415E7*RS(1,2)+541.E0/0.160083E8*RS(1,4))*RS(2,4)**2*RS(2,3)**2
     #+(76.E0/0.2401245E7*RS(1,2)+494.E0/0.4002075E7*RS(1,4)+254.E0/0.40
     #02075E7*RS(1,3))*RS(2,4)**3*RS(2,3)+(881.E0/0.800415E7*RS(1,2)+113
     #57.E0/0.160083E8*RS(1,4)+733.E0/0.160083E8*RS(1,3))*RS(2,4)**4)*RS
     #(2,2)
      s8 = s9+s10
      s5 = s7+s8
      s7 = s5
      s9 = (RS(1,1)/70+13.E0/23716.E0*RS(1,4)+13.E0/23716.E0*RS(1,3)-RS(
     #1,2)/105)*RS(2,2)**3+((541.E0/0.160083E8*RS(1,4)-13.E0/23716.E0*RS
     #(1,2)-137.E0/0.160083E7*RS(1,3))*RS(2,3)+(-463.E0/0.53361E7*RS(1,3
     #)-13.E0/23716.E0*RS(1,2)-19.E0/148225.E0*RS(1,4))*RS(2,4))*RS(2,2)
     #**2+((137.E0/0.160083E7*RS(1,2)+59.E0/0.213444E7*RS(1,3)-337.E0/0.
     #320166E8*RS(1,4))*RS(2,3)**2+(212.E0/0.4002075E7*RS(1,2)-212.E0/0.
     #4002075E7*RS(1,4))*RS(2,4)*RS(2,3)+(19.E0/148225.E0*RS(1,2)+13.E0/
     #23716.E0*RS(1,4)+463.E0/0.53361E7*RS(1,3))*RS(2,4)**2)*RS(2,2)+(-5
     #9.E0/0.213444E7*RS(1,2)+59.E0/0.213444E7*RS(1,4))*RS(2,3)**3+(337.
     #E0/0.320166E8*RS(1,2)-137.E0/0.160083E7*RS(1,4)-59.E0/0.213444E7*R
     #S(1,3))*RS(2,4)*RS(2,3)**2+(-541.E0/0.160083E8*RS(1,2)+137.E0/0.16
     #0083E7*RS(1,3)+13.E0/23716.E0*RS(1,4))*RS(2,4)**2*RS(2,3)+(-RS(1,1
     #)/70-13.E0/23716.E0*RS(1,2)+RS(1,4)/105-13.E0/23716.E0*RS(1,3))*RS
     #(2,4)**3
      s10 = RS(2,1)**2
      s8 = s9*s10
      s6 = s7+s8
      s7 = s6+(-11357.E0/0.160083E8*RS(1,2)+11357.E0/0.160083E8*RS(1,4))
     #*RS(2,3)**5
      s8 = s7
      s12 = ((-137.E0/0.160083E7*RS(1,3)+19.E0/148225.E0*RS(1,2)-149.E0/
     #0.800415E7*RS(1,4))*RS(2,3)**2+(212.E0/0.4002075E7*RS(1,2)+191.E0/
     #0.160083E8*RS(1,4)+163.E0/0.53361E7*RS(1,3))*RS(2,4)*RS(2,3)+(137.
     #E0/0.160083E7*RS(1,2)-137.E0/0.160083E7*RS(1,4))*RS(2,4)**2)*RS(2,
     #2)**2+((137.E0/0.160083E7*RS(1,2)-76.E0/0.2401245E7*RS(1,4)-881.E0
     #/0.800415E7*RS(1,3))*RS(2,3)**3+(-191.E0/0.160083E8*RS(1,2)+191.E0
     #/0.160083E8*RS(1,4))*RS(2,4)*RS(2,3)**2+(-191.E0/0.160083E8*RS(1,2
     #)-163.E0/0.53361E7*RS(1,3)-212.E0/0.4002075E7*RS(1,4))*RS(2,4)**2*
     #RS(2,3)+(254.E0/0.4002075E7*RS(1,3)+137.E0/0.160083E7*RS(1,2)+499.
     #E0/686070.E0*RS(1,4))*RS(2,4)**3)*RS(2,2)+(881.E0/0.800415E7*RS(1,
     #2)-881.E0/0.800415E7*RS(1,4))*RS(2,3)**4+(-137.E0/0.160083E7*RS(1,
     #4)+881.E0/0.800415E7*RS(1,3)+76.E0/0.2401245E7*RS(1,2))*RS(2,4)*RS
     #(2,3)**3
      s11 = s12+(149.E0/0.800415E7*RS(1,2)+137.E0/0.160083E7*RS(1,3)-19.
     #E0/148225.E0*RS(1,4))*RS(2,4)**2*RS(2,3)**2+(17.E0/190575.E0*RS(1,
     #3)-RS(1,2)/210+499.E0/686070.E0*RS(1,4)+RS(1,1)/105)*RS(2,2)**4+(1
     #9.E0/148225.E0*RS(1,3)+17.E0/190575.E0*RS(1,4)-2.E0/190575.E0*RS(1
     #,2))*RS(2,4)**3*RS(2,3)+((2.E0/190575.E0*RS(1,4)-17.E0/190575.E0*R
     #S(1,2)-19.E0/148225.E0*RS(1,3))*RS(2,3)+(-137.E0/0.160083E7*RS(1,4
     #)-499.E0/686070.E0*RS(1,2)-254.E0/0.4002075E7*RS(1,3))*RS(2,4))*RS
     #(2,2)**3+(-RS(1,1)/105-17.E0/190575.E0*RS(1,3)-499.E0/686070.E0*RS
     #(1,2)+RS(1,4)/210)*RS(2,4)**4
      s12 = RS(2,1)
      s10 = s11*s12
      s11 = (-43.E0/254100.E0*RS(1,2)-11357.E0/0.160083E8*RS(1,3)+499.E0
     #/686070.E0*RS(1,4))*RS(2,4)*RS(2,3)**4
      s9 = s10+s11
      s4 = s8+s9
      s5 = s4+(11357.E0/0.160083E8*RS(1,4)-37517.E0/0.480249E8*RS(1,3)+R
     #S(1,1)/210)*RS(2,2)**5+(13.E0/23716.E0*RS(1,4)-173.E0/0.17787E7*RS
     #(1,2)-499.E0/686070.E0*RS(1,3))*RS(2,4)**2*RS(2,3)**3+((37517.E0/0
     #.480249E8*RS(1,2)+17.E0/190575.E0*RS(1,3)+43.E0/254100.E0*RS(1,4))
     #*RS(2,3)+(-881.E0/0.800415E7*RS(1,4)-11357.E0/0.160083E8*RS(1,2)-7
     #33.E0/0.160083E8*RS(1,3))*RS(2,4))*RS(2,2)**4
      s6 = s5+(-173.E0/0.17787E7*RS(1,2)+17.E0/190575.E0*RS(1,4)-13.E0/2
     #3716.E0*RS(1,3))*RS(2,4)**3*RS(2,3)**2+((-17.E0/190575.E0*RS(1,2)+
     #13.E0/23716.E0*RS(1,3)+173.E0/0.17787E7*RS(1,4))*RS(2,3)**2+(-76.E
     #0/0.2401245E7*RS(1,4)-254.E0/0.4002075E7*RS(1,3)-494.E0/0.4002075E
     #7*RS(1,2))*RS(2,4)*RS(2,3)+(2029.E0/0.960498E8*RS(1,3)+881.E0/0.80
     #0415E7*RS(1,2)+59.E0/0.213444E7*RS(1,4))*RS(2,4)**2)*RS(2,2)**3
      s3 = s6+(-17.E0/190575.E0*RS(1,3)-37517.E0/0.480249E8*RS(1,4)-43.E
     #0/254100.E0*RS(1,2))*RS(2,4)**4*RS(2,3)+((173.E0/0.17787E7*RS(1,4)
     #+499.E0/686070.E0*RS(1,3)-13.E0/23716.E0*RS(1,2))*RS(2,3)**3+(-541
     #.E0/0.160083E8*RS(1,2)-463.E0/0.53361E7*RS(1,3)-149.E0/0.800415E7*
     #RS(1,4))*RS(2,4)*RS(2,3)**2+(-337.E0/0.320166E8*RS(1,4)+337.E0/0.3
     #20166E8*RS(1,2))*RS(2,4)**2*RS(2,3)+(-59.E0/0.213444E7*RS(1,2)-202
     #9.E0/0.960498E8*RS(1,3)-881.E0/0.800415E7*RS(1,4))*RS(2,4)**3)*RS(
     #2,2)**2+(37517.E0/0.480249E8*RS(1,3)-RS(1,1)/210-11357.E0/0.160083
     #E8*RS(1,2))*RS(2,4)**5
	ASSSS(1) = S3
C      s4 = phi(1)
C      s2 = s3*s4

      s6 = (37517.E0/0.480249E8*RS(1,4)-RS(1,2)/210-11357.E0/0.160083E8*
     #RS(1,3))*RS(2,1)**5+((499.E0/686070.E0*RS(1,4)-RS(1,3)/70+17.E0/19
     #0575.E0*RS(1,1)+2.E0/105.E0*RS(1,2))*RS(2,3)**2+(-499.E0/686070.E0
     #*RS(1,3)+494.E0/0.4002075E7*RS(1,1)-881.E0/0.800415E7*RS(1,4))*RS(
     #2,4)*RS(2,3)+(-881.E0/0.800415E7*RS(1,1)+881.E0/0.800415E7*RS(1,3)
     #)*RS(2,4)**2)*RS(2,2)**3+((-RS(1,2)/105-499.E0/686070.E0*RS(1,3)-1
     #7.E0/190575.E0*RS(1,4)+RS(1,1)/210)*RS(2,2)+(881.E0/0.800415E7*RS(
     #1,3)+11357.E0/0.160083E8*RS(1,1)+733.E0/0.160083E8*RS(1,4))*RS(2,3
     #)+(-37517.E0/0.480249E8*RS(1,1)-17.E0/190575.E0*RS(1,4)-43.E0/2541
     #00.E0*RS(1,3))*RS(2,4))*RS(2,1)**4+(43.E0/254100.E0*RS(1,1)+17.E0/
     #190575.E0*RS(1,4)+37517.E0/0.480249E8*RS(1,3))*RS(2,4)*RS(2,3)**4
      s7 = s6+((RS(1,2)/70+13.E0/23716.E0*RS(1,4)+13.E0/23716.E0*RS(1,1)
     #-RS(1,3)/105)*RS(2,3)**3+(-13.E0/23716.E0*RS(1,3)-137.E0/0.160083E
     #7*RS(1,4)+541.E0/0.160083E8*RS(1,1))*RS(2,4)*RS(2,3)**2+(-337.E0/0
     #.320166E8*RS(1,1)+59.E0/0.213444E7*RS(1,4)+137.E0/0.160083E7*RS(1,
     #3))*RS(2,4)**2*RS(2,3)+(59.E0/0.213444E7*RS(1,1)-59.E0/0.213444E7*
     #RS(1,3))*RS(2,4)**3)*RS(2,2)**2
      s5 = s7+(RS(1,1)/42-RS(1,3)/42)*RS(2,2)**5+((-RS(1,2)/70+RS(1,1)/1
     #05-13.E0/23716.E0*RS(1,4)-13.E0/23716.E0*RS(1,3))*RS(2,2)**2+((254
     #.E0/0.4002075E7*RS(1,4)+499.E0/686070.E0*RS(1,1)+137.E0/0.160083E7
     #*RS(1,3))*RS(2,3)+(17.E0/190575.E0*RS(1,1)+19.E0/148225.E0*RS(1,4)
     #-2.E0/190575.E0*RS(1,3))*RS(2,4))*RS(2,2)+(-881.E0/0.800415E7*RS(1
     #,1)-59.E0/0.213444E7*RS(1,3)-2029.E0/0.960498E8*RS(1,4))*RS(2,3)**
     #2+(494.E0/0.4002075E7*RS(1,1)+254.E0/0.4002075E7*RS(1,4)+76.E0/0.2
     #401245E7*RS(1,3))*RS(2,4)*RS(2,3)+(17.E0/190575.E0*RS(1,1)-173.E0/
     #0.17787E7*RS(1,3)-13.E0/23716.E0*RS(1,4))*RS(2,4)**2)*RS(2,1)**3+(
     #13.E0/23716.E0*RS(1,4)-17.E0/190575.E0*RS(1,3)+173.E0/0.17787E7*RS
     #(1,1))*RS(2,4)**2*RS(2,3)**3
      s7 = s5+((499.E0/686070.E0*RS(1,1)+17.E0/190575.E0*RS(1,4)-RS(1,3)
     #/210+RS(1,2)/105)*RS(2,3)**4+(-17.E0/190575.E0*RS(1,3)-19.E0/14822
     #5.E0*RS(1,4)+2.E0/190575.E0*RS(1,1))*RS(2,4)*RS(2,3)**3+(-149.E0/0
     #.800415E7*RS(1,1)-137.E0/0.160083E7*RS(1,4)+19.E0/148225.E0*RS(1,3
     #))*RS(2,4)**2*RS(2,3)**2+(137.E0/0.160083E7*RS(1,3)-881.E0/0.80041
     #5E7*RS(1,4)-76.E0/0.2401245E7*RS(1,1))*RS(2,4)**3*RS(2,3)+(-881.E0
     #/0.800415E7*RS(1,1)+881.E0/0.800415E7*RS(1,3))*RS(2,4)**4)*RS(2,2)
      s8 = s7
      s10 = ((-2.E0/105.E0*RS(1,3)+11357.E0/0.160083E8*RS(1,4)-37517.E0/
     #0.480249E8*RS(1,1)+RS(1,2)/42)*RS(2,3)+(11357.E0/0.160083E8*RS(1,1
     #)-11357.E0/0.160083E8*RS(1,3))*RS(2,4))*RS(2,2)**4
      s12 = (-17.E0/190575.E0*RS(1,3)+RS(1,1)/70-499.E0/686070.E0*RS(1,4
     #)-2.E0/105.E0*RS(1,2))*RS(2,2)**3+((13.E0/23716.E0*RS(1,1)+463.E0/
     #0.53361E7*RS(1,4)+19.E0/148225.E0*RS(1,3))*RS(2,3)+(-541.E0/0.1600
     #83E8*RS(1,3)+13.E0/23716.E0*RS(1,1)+137.E0/0.160083E7*RS(1,4))*RS(
     #2,4))*RS(2,2)**2+((-137.E0/0.160083E7*RS(1,1)+137.E0/0.160083E7*RS
     #(1,3))*RS(2,3)**2+(-191.E0/0.160083E8*RS(1,3)-212.E0/0.4002075E7*R
     #S(1,1)-163.E0/0.53361E7*RS(1,4))*RS(2,4)*RS(2,3)+(-19.E0/148225.E0
     #*RS(1,1)+137.E0/0.160083E7*RS(1,4)+149.E0/0.800415E7*RS(1,3))*RS(2
     #,4)**2)*RS(2,2)+(2029.E0/0.960498E8*RS(1,4)+881.E0/0.800415E7*RS(1
     #,3)+59.E0/0.213444E7*RS(1,1))*RS(2,3)**3+(-337.E0/0.320166E8*RS(1,
     #1)+337.E0/0.320166E8*RS(1,3))*RS(2,4)*RS(2,3)**2+(463.E0/0.53361E7
     #*RS(1,4)+541.E0/0.160083E8*RS(1,1)+149.E0/0.800415E7*RS(1,3))*RS(2
     #,4)**2*RS(2,3)+(-499.E0/686070.E0*RS(1,4)-173.E0/0.17787E7*RS(1,3)
     #+13.E0/23716.E0*RS(1,1))*RS(2,4)**3
      s13 = RS(2,1)**2
      s11 = s12*s13
      s9 = s10+s11
      s6 = s8+s9
      s7 = s6+(499.E0/686070.E0*RS(1,4)+173.E0/0.17787E7*RS(1,1)-13.E0/2
     #3716.E0*RS(1,3))*RS(2,4)**3*RS(2,3)**2+(11357.E0/0.160083E8*RS(1,1
     #)+RS(1,2)/210-37517.E0/0.480249E8*RS(1,4))*RS(2,3)**5
      s8 = s7+(43.E0/254100.E0*RS(1,1)+11357.E0/0.160083E8*RS(1,4)-499.E
     #0/686070.E0*RS(1,3))*RS(2,4)**4*RS(2,3)
      s9 = s8
      s13 = ((17.E0/190575.E0*RS(1,1)-17.E0/190575.E0*RS(1,3))*RS(2,3)+(
     #-494.E0/0.4002075E7*RS(1,3)+499.E0/686070.E0*RS(1,1)+881.E0/0.8004
     #15E7*RS(1,4))*RS(2,4))*RS(2,2)**3+(499.E0/686070.E0*RS(1,1)-43.E0/
     #254100.E0*RS(1,3)-11357.E0/0.160083E8*RS(1,4))*RS(2,4)**4+((-463.E
     #0/0.53361E7*RS(1,4)-13.E0/23716.E0*RS(1,3)-19.E0/148225.E0*RS(1,1)
     #)*RS(2,3)**2+(-212.E0/0.4002075E7*RS(1,1)+212.E0/0.4002075E7*RS(1,
     #3))*RS(2,4)*RS(2,3)+(337.E0/0.320166E8*RS(1,3)-59.E0/0.213444E7*RS
     #(1,4)-137.E0/0.160083E7*RS(1,1))*RS(2,4)**2)*RS(2,2)**2+((-254.E0/
     #0.4002075E7*RS(1,4)-137.E0/0.160083E7*RS(1,1)-499.E0/686070.E0*RS(
     #1,3))*RS(2,3)**3+(212.E0/0.4002075E7*RS(1,3)+191.E0/0.160083E8*RS(
     #1,1)+163.E0/0.53361E7*RS(1,4))*RS(2,4)*RS(2,3)**2+(-191.E0/0.16008
     #3E8*RS(1,3)+191.E0/0.160083E8*RS(1,1))*RS(2,4)**2*RS(2,3)+(881.E0/
     #0.800415E7*RS(1,4)-137.E0/0.160083E7*RS(1,1)+76.E0/0.2401245E7*RS(
     #1,3))*RS(2,4)**3)*RS(2,2)
      s12 = s13+(-11357.E0/0.160083E8*RS(1,3)-733.E0/0.160083E8*RS(1,4)-
     #881.E0/0.800415E7*RS(1,1))*RS(2,3)**4+(-RS(1,2)/42+2.E0/105.E0*RS(
     #1,1)+37517.E0/0.480249E8*RS(1,3)-11357.E0/0.160083E8*RS(1,4))*RS(2
     #,2)**4+(-254.E0/0.4002075E7*RS(1,4)-76.E0/0.2401245E7*RS(1,1)-494.
     #E0/0.4002075E7*RS(1,3))*RS(2,4)*RS(2,3)**3+(-2.E0/190575.E0*RS(1,3
     #)+2.E0/190575.E0*RS(1,1))*RS(2,4)**3*RS(2,3)+(-149.E0/0.800415E7*R
     #S(1,1)-541.E0/0.160083E8*RS(1,3)-463.E0/0.53361E7*RS(1,4))*RS(2,4)
     #**2*RS(2,3)**2
      s13 = RS(2,1)
      s11 = s12*s13
      s12 = (11357.E0/0.160083E8*RS(1,1)-11357.E0/0.160083E8*RS(1,3))*RS
     #(2,4)**5
      s10 = s11+s12
      s4 = s9+s10
	ASSSS(2) = S4
C      s5 = phi(2)
C      s3 = s4*s5
C      s1 = s2+s3
C      s3 = s1

      s7 = (-11357.E0/0.160083E8*RS(1,4)+11357.E0/0.160083E8*RS(1,2))*RS
     #(2,1)**5+(-RS(1,4)/42+RS(1,2)/42)*RS(2,3)**5+((-13.E0/23716.E0*RS(
     #1,4)-13.E0/23716.E0*RS(1,1)-RS(1,3)/70+RS(1,2)/105)*RS(2,3)**2+(13
     #7.E0/0.160083E7*RS(1,4)+254.E0/0.4002075E7*RS(1,1)+499.E0/686070.E
     #0*RS(1,2))*RS(2,4)*RS(2,3)+(-2029.E0/0.960498E8*RS(1,1)-881.E0/0.8
     #00415E7*RS(1,2)-59.E0/0.213444E7*RS(1,4))*RS(2,4)**2)*RS(2,2)**3+(
     #17.E0/190575.E0*RS(1,2)+499.E0/686070.E0*RS(1,1)+2.E0/105.E0*RS(1,
     #3)-RS(1,4)/70)*RS(2,4)**2*RS(2,3)**3
      s9 = s7
      s12 = (-17.E0/190575.E0*RS(1,4)+2.E0/190575.E0*RS(1,2)-19.E0/14822
     #5.E0*RS(1,1))*RS(2,4)**3*RS(2,3)+(-137.E0/0.160083E7*RS(1,1)+541.E
     #0/0.160083E8*RS(1,2)-13.E0/23716.E0*RS(1,4))*RS(2,4)**2*RS(2,3)**2
     #+((-2.E0/190575.E0*RS(1,4)+17.E0/190575.E0*RS(1,2)+19.E0/148225.E0
     #*RS(1,1))*RS(2,3)+(254.E0/0.4002075E7*RS(1,1)+494.E0/0.4002075E7*R
     #S(1,2)+76.E0/0.2401245E7*RS(1,4))*RS(2,4))*RS(2,2)**3+(37517.E0/0.
     #480249E8*RS(1,4)+43.E0/254100.E0*RS(1,2)+17.E0/190575.E0*RS(1,1))*
     #RS(2,4)**4
      s11 = s12+((-541.E0/0.160083E8*RS(1,4)+137.E0/0.160083E7*RS(1,1)+1
     #3.E0/23716.E0*RS(1,2))*RS(2,3)**2+(-163.E0/0.53361E7*RS(1,1)-212.E
     #0/0.4002075E7*RS(1,2)-191.E0/0.160083E8*RS(1,4))*RS(2,4)*RS(2,3)+(
     #-337.E0/0.320166E8*RS(1,2)+337.E0/0.320166E8*RS(1,4))*RS(2,4)**2)*
     #RS(2,2)**2+(-43.E0/254100.E0*RS(1,4)-17.E0/190575.E0*RS(1,1)-37517
     #.E0/0.480249E8*RS(1,2))*RS(2,2)**4+((499.E0/686070.E0*RS(1,2)-494.
     #E0/0.4002075E7*RS(1,4)+881.E0/0.800415E7*RS(1,1))*RS(2,3)**3+(-212
     #.E0/0.4002075E7*RS(1,2)+212.E0/0.4002075E7*RS(1,4))*RS(2,4)*RS(2,3
     #)**2+(191.E0/0.160083E8*RS(1,2)+163.E0/0.53361E7*RS(1,1)+212.E0/0.
     #4002075E7*RS(1,4))*RS(2,4)**2*RS(2,3)+(-494.E0/0.4002075E7*RS(1,4)
     #-254.E0/0.4002075E7*RS(1,1)-76.E0/0.2401245E7*RS(1,2))*RS(2,4)**3)
     #*RS(2,2)+(-881.E0/0.800415E7*RS(1,1)-499.E0/686070.E0*RS(1,4)+494.
     #E0/0.4002075E7*RS(1,2))*RS(2,4)*RS(2,3)**3+(-11357.E0/0.160083E8*R
     #S(1,4)+11357.E0/0.160083E8*RS(1,2))*RS(2,3)**4
      s12 = RS(2,1)
      s10 = s11*s12
      s8 = s9+s10
      s9 = s8+((-RS(1,3)/42-11357.E0/0.160083E8*RS(1,1)+2.E0/105.E0*RS(1
     #,2)+37517.E0/0.480249E8*RS(1,4))*RS(2,3)**4+(-17.E0/190575.E0*RS(1
     #,4)+17.E0/190575.E0*RS(1,2))*RS(2,4)*RS(2,3)**3+(-19.E0/148225.E0*
     #RS(1,2)-13.E0/23716.E0*RS(1,4)-463.E0/0.53361E7*RS(1,1))*RS(2,4)**
     #2*RS(2,3)**2+(-499.E0/686070.E0*RS(1,4)-254.E0/0.4002075E7*RS(1,1)
     #-137.E0/0.160083E7*RS(1,2))*RS(2,4)**3*RS(2,3)+(-11357.E0/0.160083
     #E8*RS(1,4)-881.E0/0.800415E7*RS(1,2)-733.E0/0.160083E8*RS(1,1))*RS
     #(2,4)**4)*RS(2,2)
      s6 = s9+((499.E0/686070.E0*RS(1,2)-11357.E0/0.160083E8*RS(1,1)-43.
     #E0/254100.E0*RS(1,4))*RS(2,2)+(-881.E0/0.800415E7*RS(1,2)+881.E0/0
     #.800415E7*RS(1,4))*RS(2,3)+(11357.E0/0.160083E8*RS(1,1)-499.E0/686
     #070.E0*RS(1,4)+43.E0/254100.E0*RS(1,2))*RS(2,4))*RS(2,1)**4+((RS(1
     #,2)/210-17.E0/190575.E0*RS(1,1)-RS(1,3)/105-499.E0/686070.E0*RS(1,
     #4))*RS(2,3)+(881.E0/0.800415E7*RS(1,4)+11357.E0/0.160083E8*RS(1,2)
     #+733.E0/0.160083E8*RS(1,1))*RS(2,4))*RS(2,2)**4
      s8 = s6+(11357.E0/0.160083E8*RS(1,1)-37517.E0/0.480249E8*RS(1,2)-2
     #.E0/105.E0*RS(1,4)+RS(1,3)/42)*RS(2,4)*RS(2,3)**4
      s9 = s8
      s12 = (-13.E0/23716.E0*RS(1,1)+17.E0/190575.E0*RS(1,2)-173.E0/0.17
     #787E7*RS(1,4))*RS(2,2)**3+((-19.E0/148225.E0*RS(1,2)+137.E0/0.1600
     #83E7*RS(1,1)+149.E0/0.800415E7*RS(1,4))*RS(2,3)+(463.E0/0.53361E7*
     #RS(1,1)+541.E0/0.160083E8*RS(1,2)+149.E0/0.800415E7*RS(1,4))*RS(2,
     #4))*RS(2,2)**2+((337.E0/0.320166E8*RS(1,4)-59.E0/0.213444E7*RS(1,1
     #)-137.E0/0.160083E7*RS(1,2))*RS(2,3)**2+(191.E0/0.160083E8*RS(1,2)
     #-191.E0/0.160083E8*RS(1,4))*RS(2,4)*RS(2,3)+(-541.E0/0.160083E8*RS
     #(1,4)-149.E0/0.800415E7*RS(1,2)-463.E0/0.53361E7*RS(1,1))*RS(2,4)*
     #*2)*RS(2,2)+(-881.E0/0.800415E7*RS(1,2)+881.E0/0.800415E7*RS(1,4))
     #*RS(2,3)**3+(137.E0/0.160083E7*RS(1,4)-337.E0/0.320166E8*RS(1,2)+5
     #9.E0/0.213444E7*RS(1,1))*RS(2,4)*RS(2,3)**2+(19.E0/148225.E0*RS(1,
     #4)-149.E0/0.800415E7*RS(1,2)-137.E0/0.160083E7*RS(1,1))*RS(2,4)**2
     #*RS(2,3)+(13.E0/23716.E0*RS(1,1)-17.E0/190575.E0*RS(1,4)+173.E0/0.
     #17787E7*RS(1,2))*RS(2,4)**3
      s13 = RS(2,1)**2
      s11 = s12*s13
      s12 = ((-17.E0/190575.E0*RS(1,4)+RS(1,2)/70-499.E0/686070.E0*RS(1,
     #1)-2.E0/105.E0*RS(1,3))*RS(2,3)**3+(463.E0/0.53361E7*RS(1,1)+13.E0
     #/23716.E0*RS(1,2)+19.E0/148225.E0*RS(1,4))*RS(2,4)*RS(2,3)**2+(-13
     #7.E0/0.160083E7*RS(1,2)+137.E0/0.160083E7*RS(1,4))*RS(2,4)**2*RS(2
     #,3)+(2029.E0/0.960498E8*RS(1,1)+59.E0/0.213444E7*RS(1,2)+881.E0/0.
     #800415E7*RS(1,4))*RS(2,4)**3)*RS(2,2)**2
      s10 = s11+s12
      s7 = s9+s10
      s8 = s7+(13.E0/23716.E0*RS(1,1)+13.E0/23716.E0*RS(1,2)-RS(1,4)/105
     #+RS(1,3)/70)*RS(2,4)**3*RS(2,3)**2+(-RS(1,3)/210+37517.E0/0.480249
     #E8*RS(1,1)-11357.E0/0.160083E8*RS(1,4))*RS(2,2)**5
      s5 = s8+(RS(1,3)/105+17.E0/190575.E0*RS(1,1)+499.E0/686070.E0*RS(1
     #,2)-RS(1,4)/210)*RS(2,4)**4*RS(2,3)+((-173.E0/0.17787E7*RS(1,4)+13
     #.E0/23716.E0*RS(1,2)-499.E0/686070.E0*RS(1,1))*RS(2,2)**2+((881.E0
     #/0.800415E7*RS(1,1)+76.E0/0.2401245E7*RS(1,4)-137.E0/0.160083E7*RS
     #(1,2))*RS(2,3)+(-2.E0/190575.E0*RS(1,4)+2.E0/190575.E0*RS(1,2))*RS
     #(2,4))*RS(2,2)+(59.E0/0.213444E7*RS(1,2)-59.E0/0.213444E7*RS(1,4))
     #*RS(2,3)**2+(137.E0/0.160083E7*RS(1,4)-881.E0/0.800415E7*RS(1,1)-7
     #6.E0/0.2401245E7*RS(1,2))*RS(2,4)*RS(2,3)+(173.E0/0.17787E7*RS(1,2
     #)+499.E0/686070.E0*RS(1,1)-13.E0/23716.E0*RS(1,4))*RS(2,4)**2)*RS(
     #2,1)**3+(-37517.E0/0.480249E8*RS(1,1)+11357.E0/0.160083E8*RS(1,2)+
     #RS(1,3)/210)*RS(2,4)**5
	ASSSS(3) = S5
C      s6 = phi(3)
C      s4 = s5*s6
C      s2 = s3+s4
C      s3 = s2

      s8 = (RS(1,4)/210-37517.E0/0.480249E8*RS(1,2)+11357.E0/0.160083E8*
     #RS(1,3))*RS(2,1)**5+(11357.E0/0.160083E8*RS(1,3)-11357.E0/0.160083
     #E8*RS(1,1))*RS(2,2)**5+(37517.E0/0.480249E8*RS(1,2)-RS(1,4)/210-11
     #357.E0/0.160083E8*RS(1,1))*RS(2,3)**5+((-11357.E0/0.160083E8*RS(1,
     #2)+499.E0/686070.E0*RS(1,3)-43.E0/254100.E0*RS(1,1))*RS(2,3)+(881.
     #E0/0.800415E7*RS(1,1)-881.E0/0.800415E7*RS(1,3))*RS(2,4))*RS(2,2)*
     #*4
      s9 = s8+((13.E0/23716.E0*RS(1,2)-17.E0/190575.E0*RS(1,1)+173.E0/0.
     #17787E7*RS(1,3))*RS(2,2)**2+((-76.E0/0.2401245E7*RS(1,3)-254.E0/0.
     #4002075E7*RS(1,2)-494.E0/0.4002075E7*RS(1,1))*RS(2,3)+(-17.E0/1905
     #75.E0*RS(1,1)+2.E0/190575.E0*RS(1,3)-19.E0/148225.E0*RS(1,2))*RS(2
     #,4))*RS(2,2)+(2029.E0/0.960498E8*RS(1,2)+59.E0/0.213444E7*RS(1,3)+
     #881.E0/0.800415E7*RS(1,1))*RS(2,3)**2+(-499.E0/686070.E0*RS(1,1)-1
     #37.E0/0.160083E7*RS(1,3)-254.E0/0.4002075E7*RS(1,2))*RS(2,4)*RS(2,
     #3)+(13.E0/23716.E0*RS(1,2)+RS(1,4)/70-RS(1,1)/105+13.E0/23716.E0*R
     #S(1,3))*RS(2,4)**2)*RS(2,1)**3
      s10 = s9+(RS(1,3)/210-RS(1,4)/105-17.E0/190575.E0*RS(1,2)-499.E0/6
     #86070.E0*RS(1,1))*RS(2,4)*RS(2,3)**4
      s11 = s10
      s13 = ((13.E0/23716.E0*RS(1,3)-499.E0/686070.E0*RS(1,2)-173.E0/0.1
     #7787E7*RS(1,1))*RS(2,3)**2+(881.E0/0.800415E7*RS(1,2)-137.E0/0.160
     #083E7*RS(1,3)+76.E0/0.2401245E7*RS(1,1))*RS(2,4)*RS(2,3)+(59.E0/0.
     #213444E7*RS(1,3)-59.E0/0.213444E7*RS(1,1))*RS(2,4)**2)*RS(2,2)**3
      s15 = (-13.E0/23716.E0*RS(1,1)+499.E0/686070.E0*RS(1,2)+173.E0/0.1
     #7787E7*RS(1,3))*RS(2,2)**3+((-463.E0/0.53361E7*RS(1,2)-149.E0/0.80
     #0415E7*RS(1,3)-541.E0/0.160083E8*RS(1,1))*RS(2,3)+(-137.E0/0.16008
     #3E7*RS(1,2)+19.E0/148225.E0*RS(1,1)-149.E0/0.800415E7*RS(1,3))*RS(
     #2,4))*RS(2,2)**2+((337.E0/0.320166E8*RS(1,1)-337.E0/0.320166E8*RS(
     #1,3))*RS(2,3)**2+(163.E0/0.53361E7*RS(1,2)+212.E0/0.4002075E7*RS(1
     #,1)+191.E0/0.160083E8*RS(1,3))*RS(2,4)*RS(2,3)+(-137.E0/0.160083E7
     #*RS(1,2)-13.E0/23716.E0*RS(1,1)+541.E0/0.160083E8*RS(1,3))*RS(2,4)
     #**2)*RS(2,2)+(-881.E0/0.800415E7*RS(1,3)-59.E0/0.213444E7*RS(1,1)-
     #2029.E0/0.960498E8*RS(1,2))*RS(2,3)**3+(-137.E0/0.160083E7*RS(1,3)
     #+137.E0/0.160083E7*RS(1,1))*RS(2,4)*RS(2,3)**2+(-463.E0/0.53361E7*
     #RS(1,2)-19.E0/148225.E0*RS(1,3)-13.E0/23716.E0*RS(1,1))*RS(2,4)**2
     #*RS(2,3)+(17.E0/190575.E0*RS(1,3)-RS(1,1)/70+499.E0/686070.E0*RS(1
     #,2)+2.E0/105.E0*RS(1,4))*RS(2,4)**3
      s16 = RS(2,1)**2
      s14 = s15*s16
      s12 = s13+s14
      s7 = s11+s12
      s9 = s7+(-13.E0/23716.E0*RS(1,1)-13.E0/23716.E0*RS(1,2)-RS(1,4)/70
     #+RS(1,3)/105)*RS(2,4)**2*RS(2,3)**3
      s10 = s9
      s12 = ((-13.E0/23716.E0*RS(1,2)-173.E0/0.17787E7*RS(1,1)+17.E0/190
     #575.E0*RS(1,3))*RS(2,3)**3+(137.E0/0.160083E7*RS(1,2)+149.E0/0.800
     #415E7*RS(1,1)-19.E0/148225.E0*RS(1,3))*RS(2,4)*RS(2,3)**2+(-137.E0
     #/0.160083E7*RS(1,3)+337.E0/0.320166E8*RS(1,1)-59.E0/0.213444E7*RS(
     #1,2))*RS(2,4)**2*RS(2,3)+(881.E0/0.800415E7*RS(1,1)-881.E0/0.80041
     #5E7*RS(1,3))*RS(2,4)**3)*RS(2,2)**2
      s15 = (499.E0/686070.E0*RS(1,3)+137.E0/0.160083E7*RS(1,1)+254.E0/0
     #.4002075E7*RS(1,2))*RS(2,4)*RS(2,3)**3+(881.E0/0.800415E7*RS(1,1)+
     #733.E0/0.160083E8*RS(1,2)+11357.E0/0.160083E8*RS(1,3))*RS(2,3)**4+
     #(-17.E0/190575.E0*RS(1,1)+17.E0/190575.E0*RS(1,3))*RS(2,4)**3*RS(2
     #,3)+(463.E0/0.53361E7*RS(1,2)+13.E0/23716.E0*RS(1,3)+19.E0/148225.
     #E0*RS(1,1))*RS(2,4)**2*RS(2,3)**2
      s16 = s15+((2.E0/190575.E0*RS(1,3)-2.E0/190575.E0*RS(1,1))*RS(2,3)
     #+(-881.E0/0.800415E7*RS(1,2)+137.E0/0.160083E7*RS(1,1)-76.E0/0.240
     #1245E7*RS(1,3))*RS(2,4))*RS(2,2)**3+(11357.E0/0.160083E8*RS(1,2)+4
     #3.E0/254100.E0*RS(1,3)-499.E0/686070.E0*RS(1,1))*RS(2,2)**4
      s14 = s16+(11357.E0/0.160083E8*RS(1,2)-37517.E0/0.480249E8*RS(1,3)
     #+RS(1,4)/42-2.E0/105.E0*RS(1,1))*RS(2,4)**4+((76.E0/0.2401245E7*RS
     #(1,1)+494.E0/0.4002075E7*RS(1,3)+254.E0/0.4002075E7*RS(1,2))*RS(2,
     #3)**3+(-191.E0/0.160083E8*RS(1,1)-212.E0/0.4002075E7*RS(1,3)-163.E
     #0/0.53361E7*RS(1,2))*RS(2,4)*RS(2,3)**2+(212.E0/0.4002075E7*RS(1,1
     #)-212.E0/0.4002075E7*RS(1,3))*RS(2,4)**2*RS(2,3)+(-499.E0/686070.E
     #0*RS(1,1)-881.E0/0.800415E7*RS(1,2)+494.E0/0.4002075E7*RS(1,3))*RS
     #(2,4)**3)*RS(2,2)+((149.E0/0.800415E7*RS(1,1)+463.E0/0.53361E7*RS(
     #1,2)+541.E0/0.160083E8*RS(1,3))*RS(2,3)**2+(191.E0/0.160083E8*RS(1
     #,3)-191.E0/0.160083E8*RS(1,1))*RS(2,4)*RS(2,3)+(59.E0/0.213444E7*R
     #S(1,2)-337.E0/0.320166E8*RS(1,3)+137.E0/0.160083E7*RS(1,1))*RS(2,4
     #)**2)*RS(2,2)**2
      s15 = RS(2,1)
      s13 = s14*s15
      s11 = s12+s13
      s8 = s10+s11
      s9 = s8+(RS(1,3)/70-2.E0/105.E0*RS(1,4)-17.E0/190575.E0*RS(1,1)-49
     #9.E0/686070.E0*RS(1,2))*RS(2,4)**3*RS(2,3)**2+((-43.E0/254100.E0*R
     #S(1,1)-37517.E0/0.480249E8*RS(1,3)-17.E0/190575.E0*RS(1,2))*RS(2,3
     #)**4+(-2.E0/190575.E0*RS(1,1)+17.E0/190575.E0*RS(1,3)+19.E0/148225
     #.E0*RS(1,2))*RS(2,4)*RS(2,3)**3+(137.E0/0.160083E7*RS(1,2)-541.E0/
     #0.160083E8*RS(1,1)+13.E0/23716.E0*RS(1,3))*RS(2,4)**2*RS(2,3)**2+(
     #-494.E0/0.4002075E7*RS(1,1)+499.E0/686070.E0*RS(1,3)+881.E0/0.8004
     #15E7*RS(1,2))*RS(2,4)**3*RS(2,3)+(11357.E0/0.160083E8*RS(1,3)-1135
     #7.E0/0.160083E8*RS(1,1))*RS(2,4)**4)*RS(2,2)
      s6 = s9+(37517.E0/0.480249E8*RS(1,1)-RS(1,4)/42-11357.E0/0.160083E
     #8*RS(1,2)+2.E0/105.E0*RS(1,3))*RS(2,4)**4*RS(2,3)+((43.E0/254100.E
     #0*RS(1,3)+17.E0/190575.E0*RS(1,2)+37517.E0/0.480249E8*RS(1,1))*RS(
     #2,2)+(-733.E0/0.160083E8*RS(1,2)-881.E0/0.800415E7*RS(1,3)-11357.E
     #0/0.160083E8*RS(1,1))*RS(2,3)+(-RS(1,1)/210+499.E0/686070.E0*RS(1,
     #3)+RS(1,4)/105+17.E0/190575.E0*RS(1,2))*RS(2,4))*RS(2,1)**4+(-RS(1
     #,1)/42+RS(1,3)/42)*RS(2,4)**5
	ASSSS(4) = S6
C      s7 = phi(4)
C      s5 = s6*s7

      s10 = (56311.E0/0.1200622E8*RS(1,4)-56311.E0/0.1200622E8*RS(1,2))*
     #RS(2,1)**5+((9568.E0/0.1715175E7*RS(1,2)-56311.E0/0.1200622E8*RS(1
     #,3)-12482.E0/0.1200622E8*RS(1,4)+1817.E0/0.1200622E8*RS(1,1))*RS(2
     #,3)**4+(-92.E0/444675.E0*RS(1,2)+92.E0/444675.E0*RS(1,4))*RS(2,4)*
     #RS(2,3)**3+(49.E0/81675.E0*RS(1,4)+557.E0/0.4002075E7*RS(1,1)+32.E
     #0/0.1334025E7*RS(1,2)-1018.E0/0.1334025E7*RS(1,3))*RS(2,4)**2*RS(2
     #,3)**2+(-674.E0/0.1200622E8*RS(1,2)-10522.E0/0.1200622E8*RS(1,3)+8
     #572.E0/0.1200622E8*RS(1,4)+2624.E0/0.1200622E8*RS(1,1))*RS(2,4)**3
     #*RS(2,3)+(3463.E0/0.4002075E7*RS(1,1)-3463.E0/0.4002075E7*RS(1,3))
     #*RS(2,4)**4)*RS(2,2)
      s9 = s10+(-56311.E0/0.1200622E8*RS(1,3)+56311.E0/0.1200622E8*RS(1,
     #1))*RS(2,2)**5+((-9568.E0/0.1715175E7*RS(1,2)+12482.E0/0.1200622E8
     #*RS(1,4)+56311.E0/0.1200622E8*RS(1,1)-1817.E0/0.1200622E8*RS(1,3))
     #*RS(2,2)+(-3463.E0/0.4002075E7*RS(1,4)+3463.E0/0.4002075E7*RS(1,2)
     #)*RS(2,3)+(-12482.E0/0.1200622E8*RS(1,2)+1817.E0/0.1200622E8*RS(1,
     #3)-56311.E0/0.1200622E8*RS(1,1)+9568.E0/0.1715175E7*RS(1,4))*RS(2,
     #4))*RS(2,1)**4
      s10 = s9+(56311.E0/0.1200622E8*RS(1,2)-56311.E0/0.1200622E8*RS(1,4
     #))*RS(2,3)**5
      s11 = s10+((12482.E0/0.1200622E8*RS(1,1)-1817.E0/0.1200622E8*RS(1,
     #4)-9568.E0/0.1715175E7*RS(1,3)+56311.E0/0.1200622E8*RS(1,2))*RS(2,
     #3)+(-3463.E0/0.4002075E7*RS(1,1)+3463.E0/0.4002075E7*RS(1,3))*RS(2
     #,4))*RS(2,2)**4
      s8 = s11+((106.E0/190575.E0*RS(1,4)-3319.E0/0.1200622E8*RS(1,3)-52
     #1.E0/88935.E0*RS(1,2)+9568.E0/0.1715175E7*RS(1,1))*RS(2,2)**2+((-2
     #624.E0/0.1200622E8*RS(1,4)+10522.E0/0.1200622E8*RS(1,2)+674.E0/0.1
     #200622E8*RS(1,3)-8572.E0/0.1200622E8*RS(1,1))*RS(2,3)+(92.E0/44467
     #5.E0*RS(1,2)-92.E0/444675.E0*RS(1,4))*RS(2,4))*RS(2,2)+(-1907.E0/0
     #.1200622E8*RS(1,2)+1907.E0/0.1200622E8*RS(1,4))*RS(2,3)**2+(-10522
     #.E0/0.1200622E8*RS(1,4)+2624.E0/0.1200622E8*RS(1,2)+8572.E0/0.1200
     #622E8*RS(1,1)-674.E0/0.1200622E8*RS(1,3))*RS(2,4)*RS(2,3)+(521.E0/
     #88935.E0*RS(1,4)-9568.E0/0.1715175E7*RS(1,1)+3319.E0/0.1200622E8*R
     #S(1,3)-106.E0/190575.E0*RS(1,2))*RS(2,4)**2)*RS(2,1)**3+(12482.E0/
     #0.1200622E8*RS(1,2)+56311.E0/0.1200622E8*RS(1,3)-1817.E0/0.1200622
     #E8*RS(1,1)-9568.E0/0.1715175E7*RS(1,4))*RS(2,4)*RS(2,3)**4
      s10 = s8+((106.E0/190575.E0*RS(1,1)-521.E0/88935.E0*RS(1,3)-3319.E
     #0/0.1200622E8*RS(1,4)+9568.E0/0.1715175E7*RS(1,2))*RS(2,3)**2+(-85
     #72.E0/0.1200622E8*RS(1,2)+10522.E0/0.1200622E8*RS(1,3)+674.E0/0.12
     #00622E8*RS(1,4)-2624.E0/0.1200622E8*RS(1,1))*RS(2,4)*RS(2,3)+(1907
     #.E0/0.1200622E8*RS(1,1)-1907.E0/0.1200622E8*RS(1,3))*RS(2,4)**2)*R
     #S(2,2)**3
      s11 = s10
      s15 = (-9568.E0/0.1715175E7*RS(1,2)+3319.E0/0.1200622E8*RS(1,4)-10
     #6.E0/190575.E0*RS(1,3)+521.E0/88935.E0*RS(1,1))*RS(2,2)**3+((1018.
     #E0/0.1334025E7*RS(1,2)-32.E0/0.1334025E7*RS(1,3)-557.E0/0.4002075E
     #7*RS(1,4)-49.E0/81675.E0*RS(1,1))*RS(2,3)+(49.E0/81675.E0*RS(1,2)+
     #557.E0/0.4002075E7*RS(1,3)-1018.E0/0.1334025E7*RS(1,1)+32.E0/0.133
     #4025E7*RS(1,4))*RS(2,4))*RS(2,2)**2+((137.E0/0.1334025E7*RS(1,1)-1
     #37.E0/0.1334025E7*RS(1,3))*RS(2,3)**2+(-382.E0/0.4002075E7*RS(1,2)
     #+382.E0/0.4002075E7*RS(1,4))*RS(2,4)*RS(2,3)+(-49.E0/81675.E0*RS(1
     #,4)-557.E0/0.4002075E7*RS(1,3)+1018.E0/0.1334025E7*RS(1,1)-32.E0/0
     #.1334025E7*RS(1,2))*RS(2,4)**2)*RS(2,2)
      s14 = s15+(1907.E0/0.1200622E8*RS(1,2)-1907.E0/0.1200622E8*RS(1,4)
     #)*RS(2,3)**3+(-137.E0/0.1334025E7*RS(1,1)+137.E0/0.1334025E7*RS(1,
     #3))*RS(2,4)*RS(2,3)**2+(-1018.E0/0.1334025E7*RS(1,4)+557.E0/0.4002
     #075E7*RS(1,2)+32.E0/0.1334025E7*RS(1,3)+49.E0/81675.E0*RS(1,1))*RS
     #(2,4)**2*RS(2,3)+(-3319.E0/0.1200622E8*RS(1,2)+9568.E0/0.1715175E7
     #*RS(1,4)+106.E0/190575.E0*RS(1,3)-521.E0/88935.E0*RS(1,1))*RS(2,4)
     #**3
      s15 = RS(2,1)**2
      s13 = s14*s15
      s14 = (-521.E0/88935.E0*RS(1,4)+106.E0/190575.E0*RS(1,2)+9568.E0/0
     #.1715175E7*RS(1,3)-3319.E0/0.1200622E8*RS(1,1))*RS(2,4)**2*RS(2,3)
     #**3
      s12 = s13+s14
      s9 = s11+s12
      s11 = s9
      s13 = ((-9568.E0/0.1715175E7*RS(1,3)+521.E0/88935.E0*RS(1,2)+3319.
     #E0/0.1200622E8*RS(1,1)-106.E0/190575.E0*RS(1,4))*RS(2,3)**3+(-32.E
     #0/0.1334025E7*RS(1,4)-557.E0/0.4002075E7*RS(1,1)-49.E0/81675.E0*RS
     #(1,2)+1018.E0/0.1334025E7*RS(1,3))*RS(2,4)*RS(2,3)**2+(137.E0/0.13
     #34025E7*RS(1,2)-137.E0/0.1334025E7*RS(1,4))*RS(2,4)**2*RS(2,3)+(19
     #07.E0/0.1200622E8*RS(1,3)-1907.E0/0.1200622E8*RS(1,1))*RS(2,4)**3)
     #*RS(2,2)**2
      s17 = (-12482.E0/0.1200622E8*RS(1,3)+1817.E0/0.1200622E8*RS(1,4)+9
     #568.E0/0.1715175E7*RS(1,1)-56311.E0/0.1200622E8*RS(1,2))*RS(2,2)**
     #4+(-1817.E0/0.1200622E8*RS(1,2)+12482.E0/0.1200622E8*RS(1,3)+56311
     #.E0/0.1200622E8*RS(1,4)-9568.E0/0.1715175E7*RS(1,1))*RS(2,4)**4
      s16 = s17+((-674.E0/0.1200622E8*RS(1,1)+2624.E0/0.1200622E8*RS(1,4
     #)+8572.E0/0.1200622E8*RS(1,3)-10522.E0/0.1200622E8*RS(1,2))*RS(2,3
     #)**3+(382.E0/0.4002075E7*RS(1,2)-382.E0/0.4002075E7*RS(1,4))*RS(2,
     #4)*RS(2,3)**2+(382.E0/0.4002075E7*RS(1,3)-382.E0/0.4002075E7*RS(1,
     #1))*RS(2,4)**2*RS(2,3)+(674.E0/0.1200622E8*RS(1,2)-2624.E0/0.12006
     #22E8*RS(1,3)+10522.E0/0.1200622E8*RS(1,1)-8572.E0/0.1200622E8*RS(1
     #,4))*RS(2,4)**3)*RS(2,2)+((-1018.E0/0.1334025E7*RS(1,2)+557.E0/0.4
     #002075E7*RS(1,4)+49.E0/81675.E0*RS(1,3)+32.E0/0.1334025E7*RS(1,1))
     #*RS(2,3)**2+(382.E0/0.4002075E7*RS(1,1)-382.E0/0.4002075E7*RS(1,3)
     #)*RS(2,4)*RS(2,3)+(137.E0/0.1334025E7*RS(1,4)-137.E0/0.1334025E7*R
     #S(1,2))*RS(2,4)**2)*RS(2,2)**2
      s15 = s16+(674.E0/0.1200622E8*RS(1,1)-2624.E0/0.1200622E8*RS(1,2)-
     #8572.E0/0.1200622E8*RS(1,3)+10522.E0/0.1200622E8*RS(1,4))*RS(2,4)*
     #RS(2,3)**3+(3463.E0/0.4002075E7*RS(1,4)-3463.E0/0.4002075E7*RS(1,2
     #))*RS(2,3)**4+(-92.E0/444675.E0*RS(1,3)+92.E0/444675.E0*RS(1,1))*R
     #S(2,4)**3*RS(2,3)+(-32.E0/0.1334025E7*RS(1,1)-557.E0/0.4002075E7*R
     #S(1,2)+1018.E0/0.1334025E7*RS(1,4)-49.E0/81675.E0*RS(1,3))*RS(2,4)
     #**2*RS(2,3)**2+((92.E0/444675.E0*RS(1,3)-92.E0/444675.E0*RS(1,1))*
     #RS(2,3)+(8572.E0/0.1200622E8*RS(1,2)+2624.E0/0.1200622E8*RS(1,3)-1
     #0522.E0/0.1200622E8*RS(1,1)-674.E0/0.1200622E8*RS(1,4))*RS(2,4))*R
     #S(2,2)**3
      s16 = RS(2,1)
      s14 = s15*s16
      s12 = s13+s14
      s10 = s11+s12
      s7 = s10+(521.E0/88935.E0*RS(1,3)+3319.E0/0.1200622E8*RS(1,2)-9568
     #.E0/0.1715175E7*RS(1,4)-106.E0/190575.E0*RS(1,1))*RS(2,4)**3*RS(2,
     #3)**2+(9568.E0/0.1715175E7*RS(1,3)+1817.E0/0.1200622E8*RS(1,2)-563
     #11.E0/0.1200622E8*RS(1,4)-12482.E0/0.1200622E8*RS(1,1))*RS(2,4)**4
     #*RS(2,3)+(-56311.E0/0.1200622E8*RS(1,1)+56311.E0/0.1200622E8*RS(1,
     #3))*RS(2,4)**5
	ASSSS(5) = S7
C      s8 = phi(5)
C      s6 = s7*s8
C      s4 = s5+s6
C      t0 = s3+s4

C	------------------------
C	COMPUTE FOR FACTOR ARRSS
C	------------------------
      s6 = (-RS(2,1)**2*RS(2,4)/42+RS(2,1)**2*RS(2,2)/42-RS(2,4)**3/420-
     #RS(2,1)*RS(2,4)**2/105+RS(2,1)*RS(2,2)**2/105+RS(2,2)**3/420)*RS(1
     #,1)**3
      s9 = (-13.E0/142296.E0*RS(2,2)**2*RS(2,3)+881.E0/0.160083E8*RS(2,1
     #)*RS(2,3)**2+RS(2,2)**3/210-541.E0/0.960498E8*RS(2,3)*RS(2,4)**2+R
     #S(2,1)*RS(2,2)**2/140-13.E0/142296.E0*RS(2,2)**2*RS(2,4)-11357.E0/
     #0.160083E8*RS(2,1)**2*RS(2,3)-59.E0/0.1280664E8*RS(2,3)**3+337.E0/
     #0.1920996E9*RS(2,3)**2*RS(2,4)-247.E0/0.4002075E7*RS(2,1)*RS(2,3)*
     #RS(2,4)+37517.E0/0.480249E8*RS(2,1)**2*RS(2,4)-13.E0/142296.E0*RS(
     #2,4)**3-499.E0/0.137214E7*RS(2,1)*RS(2,3)*RS(2,2)-17.E0/381150.E0*
     #RS(2,1)*RS(2,4)**2-17.E0/381150.E0*RS(2,1)*RS(2,4)*RS(2,2)-RS(2,1)
     #**3/42+106.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)+137.E0/0.960498E
     #7*RS(2,3)**2*RS(2,2)+19.E0/889350.E0*RS(2,4)**2*RS(2,2))*RS(1,2)
      s11 = (-463.E0/0.320166E8*RS(2,2)**2*RS(2,4)+13.E0/142296.E0*RS(2,
     #2)**3+499.E0/0.137214E7*RS(2,1)*RS(2,2)**2-13.E0/142296.E0*RS(2,4)
     #**3-59.E0/0.1280664E8*RS(2,3)**2*RS(2,4)+59.E0/0.1280664E8*RS(2,3)
     #**2*RS(2,2)+11357.E0/0.160083E8*RS(2,1)**2*RS(2,2)+463.E0/0.320166
     #E8*RS(2,4)**2*RS(2,2)-499.E0/0.137214E7*RS(2,1)*RS(2,4)**2-881.E0/
     #0.160083E8*RS(2,1)*RS(2,3)*RS(2,2)+137.E0/0.960498E7*RS(2,3)*RS(2,
     #4)**2+881.E0/0.160083E8*RS(2,1)*RS(2,3)*RS(2,4)-11357.E0/0.160083E
     #8*RS(2,1)**2*RS(2,4)-137.E0/0.960498E7*RS(2,2)**2*RS(2,3))*RS(1,3)
      s12 = (13.E0/142296.E0*RS(2,2)**3-106.E0/0.1200622E8*RS(2,3)*RS(2,
     #4)*RS(2,2)+13.E0/142296.E0*RS(2,3)*RS(2,4)**2-RS(2,4)**3/210+RS(2,
     #1)**3/42+59.E0/0.1280664E8*RS(2,3)**3+541.E0/0.960498E8*RS(2,2)**2
     #*RS(2,3)-19.E0/889350.E0*RS(2,2)**2*RS(2,4)+13.E0/142296.E0*RS(2,4
     #)**2*RS(2,2)-337.E0/0.1920996E9*RS(2,3)**2*RS(2,2)+499.E0/0.137214
     #E7*RS(2,1)*RS(2,3)*RS(2,4)+17.E0/381150.E0*RS(2,1)*RS(2,2)**2-RS(2
     #,1)*RS(2,4)**2/140-137.E0/0.960498E7*RS(2,3)**2*RS(2,4)-881.E0/0.1
     #60083E8*RS(2,1)*RS(2,3)**2+17.E0/381150.E0*RS(2,1)*RS(2,4)*RS(2,2)
     #+247.E0/0.4002075E7*RS(2,1)*RS(2,3)*RS(2,2)-37517.E0/0.480249E8*RS
     #(2,1)**2*RS(2,2)+11357.E0/0.160083E8*RS(2,1)**2*RS(2,3))*RS(1,4)
      s10 = s11+s12
      s8 = s9+s10
      s9 = RS(1,1)**2
      s7 = s8*s9
      s5 = s6+s7
      s7 = s5
      s11 = (RS(2,2)**3/210-RS(2,1)**2*RS(2,2)/140-17.E0/381150.E0*RS(2,
     #1)**2*RS(2,4)+137.E0/0.960498E7*RS(2,4)**3-499.E0/0.137214E7*RS(2,
     #1)**2*RS(2,3)-RS(2,1)**3/105+137.E0/0.480249E7*RS(2,1)*RS(2,3)**2-
     #13.E0/35574.E0*RS(2,1)*RS(2,4)*RS(2,2)+212.E0/0.1200622E8*RS(2,1)*
     #RS(2,3)*RS(2,4)+137.E0/0.480249E7*RS(2,4)**2*RS(2,2)-17.E0/381150.
     #E0*RS(2,2)**2*RS(2,3)-191.E0/0.960498E8*RS(2,3)*RS(2,4)**2+137.E0/
     #0.960498E7*RS(2,3)**3+19.E0/444675.E0*RS(2,1)*RS(2,4)**2-13.E0/355
     #74.E0*RS(2,1)*RS(2,3)*RS(2,2)+212.E0/0.1200622E8*RS(2,3)*RS(2,4)*R
     #S(2,2)-191.E0/0.960498E8*RS(2,3)**2*RS(2,4)+19.E0/444675.E0*RS(2,3
     #)**2*RS(2,2)-499.E0/0.137214E7*RS(2,2)**2*RS(2,4))*RS(1,2)**2
      s14 = (881.E0/0.160083E8*RS(2,3)**3+212.E0/0.1200622E8*RS(2,1)*RS(
     #2,4)**2+13.E0/35574.E0*RS(2,1)*RS(2,2)**2-247.E0/0.4002075E7*RS(2,
     #1)**2*RS(2,4)+17.E0/381150.E0*RS(2,2)**3+107.E0/0.960498E8*RS(2,3)
     #*RS(2,4)**2+38.E0/0.2401245E7*RS(2,3)**2*RS(2,4)+137.E0/0.960498E7
     #*RS(2,3)**2*RS(2,2)-RS(2,2)**2*RS(2,4)/43659+499.E0/0.137214E7*RS(
     #2,1)**2*RS(2,2)-59.E0/0.320166E7*RS(2,1)*RS(2,3)**2-191.E0/0.96049
     #8E8*RS(2,4)**2*RS(2,2)+881.E0/0.160083E8*RS(2,1)**2*RS(2,3)+106.E0
     #/0.1200622E8*RS(2,4)**3-193.E0/0.480249E7*RS(2,1)*RS(2,4)*RS(2,2)+
     #149.E0/0.2401245E8*RS(2,3)*RS(2,4)*RS(2,2)+337.E0/0.480249E8*RS(2,
     #1)*RS(2,3)*RS(2,4)-19.E0/889350.E0*RS(2,2)**2*RS(2,3))*RS(1,3)
      s15 = (499.E0/0.137214E7*RS(2,2)**3-169.E0/0.1200622E8*RS(2,3)*RS(
     #2,4)**2-137.E0/0.960498E7*RS(2,2)**2*RS(2,4)-499.E0/0.137214E7*RS(
     #2,4)**3-787.E0/0.960498E8*RS(2,3)**2*RS(2,2)+169.E0/0.1200622E8*RS
     #(2,2)**2*RS(2,3)+13.E0/35574.E0*RS(2,1)*RS(2,2)**2+137.E0/0.960498
     #E7*RS(2,4)**2*RS(2,2)-17.E0/381150.E0*RS(2,1)**2*RS(2,4)+17.E0/381
     #150.E0*RS(2,1)**2*RS(2,2)+193.E0/0.480249E7*RS(2,1)*RS(2,3)*RS(2,2
     #)-13.E0/35574.E0*RS(2,1)*RS(2,4)**2-193.E0/0.480249E7*RS(2,1)*RS(2
     #,3)*RS(2,4)+787.E0/0.960498E8*RS(2,3)**2*RS(2,4))*RS(1,4)
      s13 = s14+s15
      s14 = RS(1,2)
      s12 = s13*s14
      s10 = s11+s12
      s11 = s10+(137.E0/0.480249E7*RS(2,1)*RS(2,4)**2-881.E0/0.160083E8*
     #RS(2,3)**2*RS(2,2)+163.E0/0.320166E8*RS(2,2)**2*RS(2,4)-137.E0/0.4
     #80249E7*RS(2,1)*RS(2,2)**2+881.E0/0.160083E8*RS(2,1)**2*RS(2,4)-19
     #.E0/889350.E0*RS(2,2)**3+59.E0/0.320166E7*RS(2,1)*RS(2,3)*RS(2,2)-
     #881.E0/0.160083E8*RS(2,1)**2*RS(2,2)+881.E0/0.160083E8*RS(2,3)**2*
     #RS(2,4)+137.E0/0.480249E7*RS(2,3)*RS(2,4)**2-163.E0/0.320166E8*RS(
     #2,4)**2*RS(2,2)+19.E0/889350.E0*RS(2,4)**3-137.E0/0.480249E7*RS(2,
     #2)**2*RS(2,3)-59.E0/0.320166E7*RS(2,1)*RS(2,3)*RS(2,4))*RS(1,3)**2
      s12 = s11
      s14 = (-106.E0/0.1200622E8*RS(2,2)**3-881.E0/0.160083E8*RS(2,3)**3
     #-881.E0/0.160083E8*RS(2,1)**2*RS(2,3)-17.E0/381150.E0*RS(2,4)**3-4
     #99.E0/0.137214E7*RS(2,1)**2*RS(2,4)+59.E0/0.320166E7*RS(2,1)*RS(2,
     #3)**2+RS(2,4)**2*RS(2,2)/43659+193.E0/0.480249E7*RS(2,1)*RS(2,4)*R
     #S(2,2)+191.E0/0.960498E8*RS(2,2)**2*RS(2,4)-137.E0/0.960498E7*RS(2
     #,3)**2*RS(2,4)-38.E0/0.2401245E7*RS(2,3)**2*RS(2,2)-212.E0/0.12006
     #22E8*RS(2,1)*RS(2,2)**2+247.E0/0.4002075E7*RS(2,1)**2*RS(2,2)-337.
     #E0/0.480249E8*RS(2,1)*RS(2,3)*RS(2,2)-13.E0/35574.E0*RS(2,1)*RS(2,
     #4)**2-149.E0/0.2401245E8*RS(2,3)*RS(2,4)*RS(2,2)-107.E0/0.960498E8
     #*RS(2,2)**2*RS(2,3)+19.E0/889350.E0*RS(2,3)*RS(2,4)**2)*RS(1,4)*RS
     #(1,3)
      s15 = (-137.E0/0.960498E7*RS(2,3)**3-137.E0/0.960498E7*RS(2,2)**3-
     #RS(2,4)**3/210+RS(2,1)**3/105+13.E0/35574.E0*RS(2,1)*RS(2,4)*RS(2,
     #2)-212.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-212.E0/0.1200622E8*R
     #S(2,3)*RS(2,4)*RS(2,2)-137.E0/0.480249E7*RS(2,1)*RS(2,3)**2+191.E0
     #/0.960498E8*RS(2,2)**2*RS(2,3)+191.E0/0.960498E8*RS(2,3)**2*RS(2,2
     #)-19.E0/444675.E0*RS(2,3)**2*RS(2,4)+17.E0/381150.E0*RS(2,3)*RS(2,
     #4)**2+RS(2,1)**2*RS(2,4)/140-19.E0/444675.E0*RS(2,1)*RS(2,2)**2+49
     #9.E0/0.137214E7*RS(2,4)**2*RS(2,2)+499.E0/0.137214E7*RS(2,1)**2*RS
     #(2,3)-137.E0/0.480249E7*RS(2,2)**2*RS(2,4)+17.E0/381150.E0*RS(2,1)
     #**2*RS(2,2)+13.E0/35574.E0*RS(2,1)*RS(2,3)*RS(2,4))*RS(1,4)**2
      s13 = s14+s15
      s9 = s12+s13
      s10 = RS(1,1)
      s8 = s9*s10
      s6 = s7+s8
      s7 = s6
      s9 = (-13.E0/142296.E0*RS(2,1)**2*RS(2,3)+19.E0/889350.E0*RS(2,1)*
     #RS(2,3)**2-541.E0/0.960498E8*RS(2,3)**2*RS(2,4)+337.E0/0.1920996E9
     #*RS(2,3)*RS(2,4)**2-17.E0/381150.E0*RS(2,1)*RS(2,3)*RS(2,2)-499.E0
     #/0.137214E7*RS(2,1)*RS(2,4)*RS(2,2)+106.E0/0.1200622E8*RS(2,1)*RS(
     #2,3)*RS(2,4)-13.E0/142296.E0*RS(2,1)**2*RS(2,4)-RS(2,1)**2*RS(2,2)
     #/210-247.E0/0.4002075E7*RS(2,3)*RS(2,4)*RS(2,2)+37517.E0/0.480249E
     #8*RS(2,2)**2*RS(2,3)-RS(2,1)*RS(2,2)**2/210+881.E0/0.160083E8*RS(2
     #,4)**2*RS(2,2)-RS(2,1)**3/420-17.E0/381150.E0*RS(2,3)**2*RS(2,2)-1
     #1357.E0/0.160083E8*RS(2,2)**2*RS(2,4)-13.E0/142296.E0*RS(2,3)**3-5
     #9.E0/0.1280664E8*RS(2,4)**3+137.E0/0.960498E7*RS(2,1)*RS(2,4)**2)*
     #RS(1,2)**3
      s12 = (337.E0/0.1920996E9*RS(2,4)**3-37517.E0/0.480249E8*RS(2,2)**
     #3+901.E0/0.640332E8*RS(2,4)**2*RS(2,2)-499.E0/0.137214E7*RS(2,3)**
     #3-1721.E0/0.160083E8*RS(2,2)**2*RS(2,4)+13.E0/142296.E0*RS(2,1)**2
     #*RS(2,2)+17.E0/381150.E0*RS(2,1)*RS(2,2)**2+137.E0/0.480249E7*RS(2
     #,1)*RS(2,3)**2-191.E0/0.960498E8*RS(2,1)*RS(2,4)**2+137.E0/0.96049
     #8E7*RS(2,1)**2*RS(2,3)+149.E0/0.2401245E8*RS(2,3)*RS(2,4)**2+19.E0
     #/889350.E0*RS(2,1)*RS(2,3)*RS(2,2)-541.E0/0.960498E8*RS(2,1)**2*RS
     #(2,4)-169.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2)-1303.E0/0.2401245
     #E8*RS(2,3)*RS(2,4)*RS(2,2)+107.E0/0.960498E8*RS(2,1)*RS(2,3)*RS(2,
     #4)-631.E0/0.320166E8*RS(2,3)**2*RS(2,4)-13.E0/47432.E0*RS(2,3)**2*
     #RS(2,2))*RS(1,3)
      s13 = (11357.E0/0.160083E8*RS(2,2)**3+901.E0/0.640332E8*RS(2,3)*RS
     #(2,4)**2-47.E0/0.53361E7*RS(2,3)*RS(2,4)*RS(2,2)+881.E0/0.160083E8
     #*RS(2,4)**3+RS(2,1)*RS(2,3)*RS(2,2)/43659+463.E0/0.320166E8*RS(2,3
     #)**3+1721.E0/0.160083E8*RS(2,2)**2*RS(2,3)-59.E0/0.426888E7*RS(2,4
     #)**2*RS(2,2)+3589.E0/0.960498E8*RS(2,3)**2*RS(2,2)+137.E0/0.960498
     #E7*RS(2,1)*RS(2,4)*RS(2,2)+499.E0/0.137214E7*RS(2,1)*RS(2,2)**2+13
     #7.E0/0.480249E7*RS(2,1)*RS(2,4)**2-163.E0/0.320166E8*RS(2,1)*RS(2,
     #3)**2+19.E0/889350.E0*RS(2,1)**2*RS(2,4)+463.E0/0.320166E8*RS(2,1)
     #**2*RS(2,3)+13.E0/142296.E0*RS(2,1)**2*RS(2,2)+149.E0/0.480249E8*R
     #S(2,3)**2*RS(2,4)-191.E0/0.960498E8*RS(2,1)*RS(2,3)*RS(2,4))*RS(1,
     #4)
      s11 = s12+s13
      s12 = RS(1,2)**2
      s10 = s11*s12
      s8 = s9+s10
      s4 = s7+s8
      s6 = s4
      s10 = (149.E0/0.480249E8*RS(2,4)**2*RS(2,2)-137.E0/0.960498E7*RS(2
     #,1)**2*RS(2,2)+337.E0/0.1920996E9*RS(2,1)**2*RS(2,4)-631.E0/0.3201
     #66E8*RS(2,3)*RS(2,4)**2-191.E0/0.960498E8*RS(2,1)*RS(2,4)**2-3589.
     #E0/0.960498E8*RS(2,2)**2*RS(2,4)-11357.E0/0.160083E8*RS(2,3)**3-19
     #.E0/444675.E0*RS(2,1)*RS(2,2)**2+881.E0/0.160083E8*RS(2,1)*RS(2,3)
     #**2-137.E0/0.960498E7*RS(2,1)*RS(2,3)*RS(2,2)-59.E0/0.1280664E8*RS
     #(2,1)**2*RS(2,3)+17.E0/381150.E0*RS(2,2)**3-43.E0/254100.E0*RS(2,3
     #)**2*RS(2,4)-101.E0/0.160083E7*RS(2,3)*RS(2,4)*RS(2,2)+787.E0/0.96
     #0498E8*RS(2,1)*RS(2,4)*RS(2,2)-541.E0/0.960498E8*RS(2,4)**3+13.E0/
     #47432.E0*RS(2,2)**2*RS(2,3)+38.E0/0.2401245E7*RS(2,1)*RS(2,3)*RS(2
     #,4))*RS(1,3)**2
      s12 = (247.E0/0.4002075E7*RS(2,2)**3+149.E0/0.2401245E8*RS(2,1)*RS
     #(2,3)*RS(2,4)+212.E0/0.1200622E8*RS(2,1)*RS(2,4)**2-101.E0/0.16008
     #3E7*RS(2,3)**2*RS(2,4)+106.E0/0.1200622E8*RS(2,1)**2*RS(2,4)-47.E0
     #/0.53361E7*RS(2,4)**2*RS(2,2)+47.E0/0.53361E7*RS(2,2)**2*RS(2,4)-2
     #12.E0/0.1200622E8*RS(2,1)*RS(2,2)**2-149.E0/0.2401245E8*RS(2,1)*RS
     #(2,3)*RS(2,2)+101.E0/0.160083E7*RS(2,3)**2*RS(2,2)-247.E0/0.400207
     #5E7*RS(2,4)**3+1303.E0/0.2401245E8*RS(2,2)**2*RS(2,3)-106.E0/0.120
     #0622E8*RS(2,1)**2*RS(2,2)-1303.E0/0.2401245E8*RS(2,3)*RS(2,4)**2)*
     #RS(1,4)*RS(1,3)
      s13 = (-463.E0/0.320166E8*RS(2,3)**3-881.E0/0.160083E8*RS(2,2)**3-
     #11357.E0/0.160083E8*RS(2,4)**3+191.E0/0.960498E8*RS(2,1)*RS(2,3)*R
     #S(2,2)+47.E0/0.53361E7*RS(2,3)*RS(2,4)*RS(2,2)+163.E0/0.320166E8*R
     #S(2,1)*RS(2,3)**2+59.E0/0.426888E7*RS(2,2)**2*RS(2,4)-901.E0/0.640
     #332E8*RS(2,2)**2*RS(2,3)-3589.E0/0.960498E8*RS(2,3)**2*RS(2,4)-172
     #1.E0/0.160083E8*RS(2,3)*RS(2,4)**2-13.E0/142296.E0*RS(2,1)**2*RS(2
     #,4)-137.E0/0.480249E7*RS(2,1)*RS(2,2)**2-137.E0/0.960498E7*RS(2,1)
     #*RS(2,4)*RS(2,2)-463.E0/0.320166E8*RS(2,1)**2*RS(2,3)-499.E0/0.137
     #214E7*RS(2,1)*RS(2,4)**2-149.E0/0.480249E8*RS(2,3)**2*RS(2,2)-RS(2
     #,1)*RS(2,3)*RS(2,4)/43659-19.E0/889350.E0*RS(2,1)**2*RS(2,2))*RS(1
     #,4)**2
      s11 = s12+s13
      s9 = s10+s11
      s10 = RS(1,2)
      s8 = s9*s10
      s9 = (13.E0/142296.E0*RS(2,2)**3+11357.E0/0.160083E8*RS(2,3)**2*RS
     #(2,2)-59.E0/0.1280664E8*RS(2,1)**2*RS(2,4)-881.E0/0.160083E8*RS(2,
     #1)*RS(2,3)*RS(2,2)+59.E0/0.1280664E8*RS(2,1)**2*RS(2,2)+499.E0/0.1
     #37214E7*RS(2,2)**2*RS(2,3)-463.E0/0.320166E8*RS(2,2)**2*RS(2,4)-13
     #7.E0/0.960498E7*RS(2,1)*RS(2,2)**2+463.E0/0.320166E8*RS(2,4)**2*RS
     #(2,2)-499.E0/0.137214E7*RS(2,3)*RS(2,4)**2+881.E0/0.160083E8*RS(2,
     #1)*RS(2,3)*RS(2,4)+137.E0/0.960498E7*RS(2,1)*RS(2,4)**2-13.E0/1422
     #96.E0*RS(2,4)**3-11357.E0/0.160083E8*RS(2,3)**2*RS(2,4))*RS(1,3)**
     #3
      s7 = s8+s9
      s5 = s6+s7
      s6 = s5+(3589.E0/0.960498E8*RS(2,4)**2*RS(2,2)+19.E0/444675.E0*RS(
     #2,1)*RS(2,4)**2-17.E0/381150.E0*RS(2,4)**3-337.E0/0.1920996E9*RS(2
     #,1)**2*RS(2,2)+137.E0/0.960498E7*RS(2,1)**2*RS(2,4)-38.E0/0.240124
     #5E7*RS(2,1)*RS(2,3)*RS(2,2)+631.E0/0.320166E8*RS(2,2)**2*RS(2,3)+1
     #1357.E0/0.160083E8*RS(2,3)**3+541.E0/0.960498E8*RS(2,2)**3-881.E0/
     #0.160083E8*RS(2,1)*RS(2,3)**2+101.E0/0.160083E7*RS(2,3)*RS(2,4)*RS
     #(2,2)-13.E0/47432.E0*RS(2,3)*RS(2,4)**2+191.E0/0.960498E8*RS(2,1)*
     #RS(2,2)**2+137.E0/0.960498E7*RS(2,1)*RS(2,3)*RS(2,4)+43.E0/254100.
     #E0*RS(2,3)**2*RS(2,2)-149.E0/0.480249E8*RS(2,2)**2*RS(2,4)+59.E0/0
     #.1280664E8*RS(2,1)**2*RS(2,3)-787.E0/0.960498E8*RS(2,1)*RS(2,4)*RS
     #(2,2))*RS(1,4)*RS(1,3)**2
      s7 = s6
      s9 = (-149.E0/0.2401245E8*RS(2,2)**2*RS(2,3)+541.E0/0.960498E8*RS(
     #2,1)**2*RS(2,2)-107.E0/0.960498E8*RS(2,1)*RS(2,3)*RS(2,2)-337.E0/0
     #.1920996E9*RS(2,2)**3+169.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2)+1
     #721.E0/0.160083E8*RS(2,4)**2*RS(2,2)+631.E0/0.320166E8*RS(2,3)**2*
     #RS(2,2)+1303.E0/0.2401245E8*RS(2,3)*RS(2,4)*RS(2,2)+37517.E0/0.480
     #249E8*RS(2,4)**3+499.E0/0.137214E7*RS(2,3)**3-901.E0/0.640332E8*RS
     #(2,2)**2*RS(2,4)-137.E0/0.960498E7*RS(2,1)**2*RS(2,3)-13.E0/142296
     #.E0*RS(2,1)**2*RS(2,4)+13.E0/47432.E0*RS(2,3)**2*RS(2,4)-17.E0/381
     #150.E0*RS(2,1)*RS(2,4)**2-19.E0/889350.E0*RS(2,1)*RS(2,3)*RS(2,4)+
     #191.E0/0.960498E8*RS(2,1)*RS(2,2)**2-137.E0/0.480249E7*RS(2,1)*RS(
     #2,3)**2)*RS(1,4)**2*RS(1,3)
      s10 = (-881.E0/0.160083E8*RS(2,2)**2*RS(2,4)+RS(2,1)**3/420+17.E0/
     #381150.E0*RS(2,1)*RS(2,3)*RS(2,4)+13.E0/142296.E0*RS(2,1)**2*RS(2,
     #2)+13.E0/142296.E0*RS(2,1)**2*RS(2,3)+17.E0/381150.E0*RS(2,3)**2*R
     #S(2,4)+11357.E0/0.160083E8*RS(2,4)**2*RS(2,2)+499.E0/0.137214E7*RS
     #(2,1)*RS(2,4)*RS(2,2)-137.E0/0.960498E7*RS(2,1)*RS(2,2)**2-37517.E
     #0/0.480249E8*RS(2,3)*RS(2,4)**2-337.E0/0.1920996E9*RS(2,2)**2*RS(2
     #,3)+59.E0/0.1280664E8*RS(2,2)**3+247.E0/0.4002075E7*RS(2,3)*RS(2,4
     #)*RS(2,2)+RS(2,1)**2*RS(2,4)/210-19.E0/889350.E0*RS(2,1)*RS(2,3)**
     #2+541.E0/0.960498E8*RS(2,3)**2*RS(2,2)-106.E0/0.1200622E8*RS(2,1)*
     #RS(2,3)*RS(2,2)+13.E0/142296.E0*RS(2,3)**3+RS(2,1)*RS(2,4)**2/210)
     #*RS(1,4)**3
      s8 = s9+s10
      s3 = s7+s8
	ARRSS(1) = S3
C      s4 = phi(1)
C      s2 = s3*s4

      s7 = (-337.E0/0.1920996E9*RS(2,3)**2*RS(2,4)+13.E0/142296.E0*RS(2,
     #4)**3-106.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)+RS(2,2)**3/420+59
     #.E0/0.1280664E8*RS(2,3)**3+13.E0/142296.E0*RS(2,2)**2*RS(2,3)-3751
     #7.E0/0.480249E8*RS(2,1)**2*RS(2,4)+RS(2,1)*RS(2,2)**2/210+11357.E0
     #/0.160083E8*RS(2,1)**2*RS(2,3)+541.E0/0.960498E8*RS(2,3)*RS(2,4)**
     #2-881.E0/0.160083E8*RS(2,1)*RS(2,3)**2+13.E0/142296.E0*RS(2,2)**2*
     #RS(2,4)-19.E0/889350.E0*RS(2,4)**2*RS(2,2)-137.E0/0.960498E7*RS(2,
     #3)**2*RS(2,2)+17.E0/381150.E0*RS(2,1)*RS(2,4)**2+RS(2,1)**2*RS(2,2
     #)/210+499.E0/0.137214E7*RS(2,1)*RS(2,3)*RS(2,2)+247.E0/0.4002075E7
     #*RS(2,1)*RS(2,3)*RS(2,4)+17.E0/381150.E0*RS(2,1)*RS(2,4)*RS(2,2))*
     #RS(1,1)**3
      s10 = (-137.E0/0.480249E7*RS(2,1)*RS(2,3)**2+RS(2,2)**3/105-137.E0
     #/0.960498E7*RS(2,3)**3+191.E0/0.960498E8*RS(2,3)*RS(2,4)**2+191.E0
     #/0.960498E8*RS(2,3)**2*RS(2,4)+499.E0/0.137214E7*RS(2,2)**2*RS(2,4
     #)+499.E0/0.137214E7*RS(2,1)**2*RS(2,3)-212.E0/0.1200622E8*RS(2,1)*
     #RS(2,3)*RS(2,4)+13.E0/35574.E0*RS(2,1)*RS(2,4)*RS(2,2)-19.E0/44467
     #5.E0*RS(2,1)*RS(2,4)**2-137.E0/0.960498E7*RS(2,4)**3-RS(2,1)**3/21
     #0+RS(2,1)*RS(2,2)**2/140-19.E0/444675.E0*RS(2,3)**2*RS(2,2)+13.E0/
     #35574.E0*RS(2,1)*RS(2,3)*RS(2,2)+17.E0/381150.E0*RS(2,1)**2*RS(2,4
     #)-212.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)+17.E0/381150.E0*RS(2,
     #2)**2*RS(2,3)-137.E0/0.480249E7*RS(2,4)**2*RS(2,2))*RS(1,2)
      s12 = (-463.E0/0.320166E8*RS(2,4)**3-881.E0/0.160083E8*RS(2,3)**3-
     #13.E0/142296.E0*RS(2,1)*RS(2,2)**2-3589.E0/0.960498E8*RS(2,1)*RS(2
     #,4)**2-19.E0/889350.E0*RS(2,2)**2*RS(2,3)-137.E0/0.480249E7*RS(2,3
     #)**2*RS(2,2)-499.E0/0.137214E7*RS(2,1)**2*RS(2,2)+163.E0/0.320166E
     #8*RS(2,4)**2*RS(2,2)-463.E0/0.320166E8*RS(2,2)**2*RS(2,4)-11357.E0
     #/0.160083E8*RS(2,1)**3-149.E0/0.480249E8*RS(2,3)*RS(2,4)**2-137.E0
     #/0.960498E7*RS(2,1)*RS(2,3)*RS(2,2)+59.E0/0.426888E7*RS(2,1)*RS(2,
     #3)**2-RS(2,1)*RS(2,4)*RS(2,2)/43659+191.E0/0.960498E8*RS(2,3)*RS(2
     #,4)*RS(2,2)+47.E0/0.53361E7*RS(2,1)*RS(2,3)*RS(2,4)-1721.E0/0.1600
     #83E8*RS(2,1)**2*RS(2,4)-901.E0/0.640332E8*RS(2,3)**2*RS(2,4))*RS(1
     #,3)
      s13 = (-337.E0/0.1920996E9*RS(2,3)**3+499.E0/0.137214E7*RS(2,4)**3
     #+37517.E0/0.480249E8*RS(2,1)**3-137.E0/0.480249E7*RS(2,4)**2*RS(2,
     #2)+191.E0/0.960498E8*RS(2,3)**2*RS(2,2)+541.E0/0.960498E8*RS(2,2)*
     #*2*RS(2,3)-137.E0/0.960498E7*RS(2,2)**2*RS(2,4)+13.E0/47432.E0*RS(
     #2,1)*RS(2,4)**2-149.E0/0.2401245E8*RS(2,3)**2*RS(2,4)-19.E0/889350
     #.E0*RS(2,1)*RS(2,4)*RS(2,2)-13.E0/142296.E0*RS(2,1)*RS(2,2)**2-107
     #.E0/0.960498E8*RS(2,3)*RS(2,4)*RS(2,2)-17.E0/381150.E0*RS(2,1)**2*
     #RS(2,2)+169.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)+631.E0/0.320166
     #E8*RS(2,3)*RS(2,4)**2+1721.E0/0.160083E8*RS(2,1)**2*RS(2,3)-901.E0
     #/0.640332E8*RS(2,1)*RS(2,3)**2+1303.E0/0.2401245E8*RS(2,1)*RS(2,3)
     #*RS(2,4))*RS(1,4)
      s11 = s12+s13
      s9 = s10+s11
      s10 = RS(1,1)**2
      s8 = s9*s10
      s6 = s7+s8
      s8 = s6
      s12 = (-19.E0/889350.E0*RS(2,1)*RS(2,3)**2-37517.E0/0.480249E8*RS(
     #2,2)**2*RS(2,3)-RS(2,1)**3/210+59.E0/0.1280664E8*RS(2,4)**3+17.E0/
     #381150.E0*RS(2,3)**2*RS(2,2)-RS(2,1)**2*RS(2,2)/140+13.E0/142296.E
     #0*RS(2,1)**2*RS(2,4)+13.E0/142296.E0*RS(2,3)**3+13.E0/142296.E0*RS
     #(2,1)**2*RS(2,3)-881.E0/0.160083E8*RS(2,4)**2*RS(2,2)+17.E0/381150
     #.E0*RS(2,1)*RS(2,3)*RS(2,2)-337.E0/0.1920996E9*RS(2,3)*RS(2,4)**2+
     #541.E0/0.960498E8*RS(2,3)**2*RS(2,4)-137.E0/0.960498E7*RS(2,1)*RS(
     #2,4)**2+11357.E0/0.160083E8*RS(2,2)**2*RS(2,4)+247.E0/0.4002075E7*
     #RS(2,3)*RS(2,4)*RS(2,2)+RS(2,2)**3/42-106.E0/0.1200622E8*RS(2,1)*R
     #S(2,3)*RS(2,4)+499.E0/0.137214E7*RS(2,1)*RS(2,4)*RS(2,2))*RS(1,2)*
     #*2
      s15 = (-17.E0/381150.E0*RS(2,1)*RS(2,2)**2+499.E0/0.137214E7*RS(2,
     #3)**3-787.E0/0.960498E8*RS(2,3)*RS(2,4)**2-13.E0/35574.E0*RS(2,1)*
     #*2*RS(2,2)+787.E0/0.960498E8*RS(2,1)*RS(2,4)**2+169.E0/0.1200622E8
     #*RS(2,3)**2*RS(2,4)+13.E0/35574.E0*RS(2,3)**2*RS(2,2)+137.E0/0.960
     #498E7*RS(2,1)**2*RS(2,3)-499.E0/0.137214E7*RS(2,1)**3-137.E0/0.960
     #498E7*RS(2,1)*RS(2,3)**2+17.E0/381150.E0*RS(2,2)**2*RS(2,3)-193.E0
     #/0.480249E7*RS(2,1)*RS(2,4)*RS(2,2)+193.E0/0.480249E7*RS(2,3)*RS(2
     #,4)*RS(2,2)-169.E0/0.1200622E8*RS(2,1)**2*RS(2,4))*RS(1,3)
      s16 = (-881.E0/0.160083E8*RS(2,4)**3-17.E0/381150.E0*RS(2,1)**3+19
     #3.E0/0.480249E7*RS(2,1)*RS(2,3)*RS(2,2)-106.E0/0.1200622E8*RS(2,3)
     #**3+247.E0/0.4002075E7*RS(2,2)**2*RS(2,3)-881.E0/0.160083E8*RS(2,2
     #)**2*RS(2,4)+59.E0/0.320166E7*RS(2,4)**2*RS(2,2)-212.E0/0.1200622E
     #8*RS(2,3)**2*RS(2,2)-499.E0/0.137214E7*RS(2,1)*RS(2,2)**2-137.E0/0
     #.960498E7*RS(2,1)*RS(2,4)**2-337.E0/0.480249E8*RS(2,3)*RS(2,4)*RS(
     #2,2)-38.E0/0.2401245E7*RS(2,3)*RS(2,4)**2+RS(2,1)**2*RS(2,3)/43659
     #-13.E0/35574.E0*RS(2,1)**2*RS(2,2)+191.E0/0.960498E8*RS(2,1)*RS(2,
     #3)**2-149.E0/0.2401245E8*RS(2,1)*RS(2,3)*RS(2,4)+19.E0/889350.E0*R
     #S(2,1)**2*RS(2,4)-107.E0/0.960498E8*RS(2,3)**2*RS(2,4))*RS(1,4)
      s14 = s15+s16
      s15 = RS(1,2)
      s13 = s14*s15
      s11 = s12+s13
      s12 = s11+(-59.E0/0.426888E7*RS(2,1)**2*RS(2,3)-163.E0/0.320166E8*
     #RS(2,4)**2*RS(2,2)+901.E0/0.640332E8*RS(2,1)**2*RS(2,4)+3589.E0/0.
     #960498E8*RS(2,3)*RS(2,4)**2+149.E0/0.480249E8*RS(2,1)*RS(2,4)**2+4
     #99.E0/0.137214E7*RS(2,3)**2*RS(2,2)+137.E0/0.480249E7*RS(2,1)**2*R
     #S(2,2)+19.E0/889350.E0*RS(2,1)*RS(2,2)**2+137.E0/0.960498E7*RS(2,1
     #)*RS(2,3)*RS(2,2)+463.E0/0.320166E8*RS(2,2)**2*RS(2,4)+881.E0/0.16
     #0083E8*RS(2,1)**3+11357.E0/0.160083E8*RS(2,3)**3+1721.E0/0.160083E
     #8*RS(2,3)**2*RS(2,4)+RS(2,3)*RS(2,4)*RS(2,2)/43659-191.E0/0.960498
     #E8*RS(2,1)*RS(2,4)*RS(2,2)+463.E0/0.320166E8*RS(2,4)**3+13.E0/1422
     #96.E0*RS(2,2)**2*RS(2,3)-47.E0/0.53361E7*RS(2,1)*RS(2,3)*RS(2,4))*
     #RS(1,3)**2
      s13 = s12
      s15 = (149.E0/0.2401245E8*RS(2,1)*RS(2,4)*RS(2,2)-1303.E0/0.240124
     #5E8*RS(2,1)**2*RS(2,4)-101.E0/0.160083E7*RS(2,1)*RS(2,4)**2-212.E0
     #/0.1200622E8*RS(2,3)**2*RS(2,2)-149.E0/0.2401245E8*RS(2,3)*RS(2,4)
     #*RS(2,2)+47.E0/0.53361E7*RS(2,1)*RS(2,3)**2-106.E0/0.1200622E8*RS(
     #2,2)**2*RS(2,3)-47.E0/0.53361E7*RS(2,1)**2*RS(2,3)+212.E0/0.120062
     #2E8*RS(2,1)**2*RS(2,2)+101.E0/0.160083E7*RS(2,3)*RS(2,4)**2-247.E0
     #/0.4002075E7*RS(2,1)**3+106.E0/0.1200622E8*RS(2,1)*RS(2,2)**2+1303
     #.E0/0.2401245E8*RS(2,3)**2*RS(2,4)+247.E0/0.4002075E7*RS(2,3)**3)*
     #RS(1,4)*RS(1,3)
      s16 = (541.E0/0.960498E8*RS(2,3)**3+11357.E0/0.160083E8*RS(2,4)**3
     #-17.E0/381150.E0*RS(2,1)**3-38.E0/0.2401245E7*RS(2,3)*RS(2,4)*RS(2
     #,2)-149.E0/0.480249E8*RS(2,1)*RS(2,3)**2+59.E0/0.1280664E8*RS(2,2)
     #**2*RS(2,4)+191.E0/0.960498E8*RS(2,3)**2*RS(2,2)-337.E0/0.1920996E
     #9*RS(2,2)**2*RS(2,3)+43.E0/254100.E0*RS(2,3)*RS(2,4)**2-13.E0/4743
     #2.E0*RS(2,1)**2*RS(2,4)+19.E0/444675.E0*RS(2,1)**2*RS(2,2)-881.E0/
     #0.160083E8*RS(2,4)**2*RS(2,2)+137.E0/0.960498E7*RS(2,1)*RS(2,2)**2
     #-787.E0/0.960498E8*RS(2,1)*RS(2,3)*RS(2,2)+137.E0/0.960498E7*RS(2,
     #1)*RS(2,4)*RS(2,2)+3589.E0/0.960498E8*RS(2,1)**2*RS(2,3)+101.E0/0.
     #160083E7*RS(2,1)*RS(2,3)*RS(2,4)+631.E0/0.320166E8*RS(2,3)**2*RS(2
     #,4))*RS(1,4)**2
      s14 = s15+s16
      s10 = s13+s14
      s11 = RS(1,1)
      s9 = s10*s11
      s7 = s8+s9
      s8 = s7
      s10 = (-RS(2,1)*RS(2,2)**2/42-RS(2,1)**2*RS(2,2)/105+RS(2,2)**2*RS
     #(2,3)/42+RS(2,3)**3/420-RS(2,1)**3/420+RS(2,3)**2*RS(2,2)/105)*RS(
     #1,2)**3
      s13 = (-59.E0/0.1280664E8*RS(2,4)**3-11357.E0/0.160083E8*RS(2,2)**
     #2*RS(2,4)+RS(2,3)**3/210+137.E0/0.960498E7*RS(2,3)*RS(2,4)**2+881.
     #E0/0.160083E8*RS(2,4)**2*RS(2,2)-541.E0/0.960498E8*RS(2,1)**2*RS(2
     #,4)+RS(2,3)**2*RS(2,2)/140-17.E0/381150.E0*RS(2,1)**2*RS(2,2)-RS(2
     #,2)**3/42-13.E0/142296.E0*RS(2,1)*RS(2,3)**2+337.E0/0.1920996E9*RS
     #(2,1)*RS(2,4)**2+19.E0/889350.E0*RS(2,1)**2*RS(2,3)-13.E0/142296.E
     #0*RS(2,1)**3+37517.E0/0.480249E8*RS(2,1)*RS(2,2)**2-17.E0/381150.E
     #0*RS(2,1)*RS(2,3)*RS(2,2)-247.E0/0.4002075E7*RS(2,1)*RS(2,4)*RS(2,
     #2)-499.E0/0.137214E7*RS(2,3)*RS(2,4)*RS(2,2)+106.E0/0.1200622E8*RS
     #(2,1)*RS(2,3)*RS(2,4)-13.E0/142296.E0*RS(2,3)**2*RS(2,4))*RS(1,3)
      s14 = (13.E0/142296.E0*RS(2,3)**3-881.E0/0.160083E8*RS(2,3)*RS(2,4
     #)*RS(2,2)-13.E0/142296.E0*RS(2,1)**3+463.E0/0.320166E8*RS(2,1)**2*
     #RS(2,3)-11357.E0/0.160083E8*RS(2,1)*RS(2,2)**2+11357.E0/0.160083E8
     #*RS(2,2)**2*RS(2,3)+499.E0/0.137214E7*RS(2,3)**2*RS(2,2)+137.E0/0.
     #960498E7*RS(2,1)**2*RS(2,4)+59.E0/0.1280664E8*RS(2,3)*RS(2,4)**2-5
     #9.E0/0.1280664E8*RS(2,1)*RS(2,4)**2-463.E0/0.320166E8*RS(2,1)*RS(2
     #,3)**2+881.E0/0.160083E8*RS(2,1)*RS(2,4)*RS(2,2)-499.E0/0.137214E7
     #*RS(2,1)**2*RS(2,2)-137.E0/0.960498E7*RS(2,3)**2*RS(2,4))*RS(1,4)
      s12 = s13+s14
      s13 = RS(1,2)**2
      s11 = s12*s13
      s9 = s10+s11
      s5 = s8+s9
      s7 = s5
      s11 = (-191.E0/0.960498E8*RS(2,1)**2*RS(2,4)+137.E0/0.960498E7*RS(
     #2,1)**3+19.E0/444675.E0*RS(2,1)**2*RS(2,2)-RS(2,2)**2*RS(2,3)/140+
     #RS(2,3)**3/210-RS(2,2)**3/105-499.E0/0.137214E7*RS(2,1)*RS(2,3)**2
     #-499.E0/0.137214E7*RS(2,2)**2*RS(2,4)+137.E0/0.480249E7*RS(2,1)**2
     #*RS(2,3)+137.E0/0.960498E7*RS(2,4)**3-13.E0/35574.E0*RS(2,1)*RS(2,
     #3)*RS(2,2)-17.E0/381150.E0*RS(2,3)**2*RS(2,4)+19.E0/444675.E0*RS(2
     #,3)*RS(2,4)**2+137.E0/0.480249E7*RS(2,4)**2*RS(2,2)-13.E0/35574.E0
     #*RS(2,3)*RS(2,4)*RS(2,2)+212.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2
     #)-17.E0/381150.E0*RS(2,1)*RS(2,2)**2-191.E0/0.960498E8*RS(2,1)*RS(
     #2,4)**2+212.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4))*RS(1,3)**2
      s13 = (-19.E0/889350.E0*RS(2,3)**2*RS(2,4)-59.E0/0.320166E7*RS(2,4
     #)**2*RS(2,2)+499.E0/0.137214E7*RS(2,2)**2*RS(2,3)-247.E0/0.4002075
     #E7*RS(2,1)*RS(2,2)**2-191.E0/0.960498E8*RS(2,1)**2*RS(2,3)+106.E0/
     #0.1200622E8*RS(2,1)**3+38.E0/0.2401245E7*RS(2,1)*RS(2,4)**2+137.E0
     #/0.960498E7*RS(2,3)*RS(2,4)**2+17.E0/381150.E0*RS(2,3)**3-RS(2,1)*
     #RS(2,3)**2/43659+13.E0/35574.E0*RS(2,3)**2*RS(2,2)+881.E0/0.160083
     #E8*RS(2,4)**3-193.E0/0.480249E7*RS(2,1)*RS(2,3)*RS(2,2)+149.E0/0.2
     #401245E8*RS(2,1)*RS(2,3)*RS(2,4)+337.E0/0.480249E8*RS(2,1)*RS(2,4)
     #*RS(2,2)+107.E0/0.960498E8*RS(2,1)**2*RS(2,4)+212.E0/0.1200622E8*R
     #S(2,1)**2*RS(2,2)+881.E0/0.160083E8*RS(2,2)**2*RS(2,4))*RS(1,4)*RS
     #(1,3)
      s14 = (-19.E0/889350.E0*RS(2,3)**3+19.E0/889350.E0*RS(2,1)**3-163.
     #E0/0.320166E8*RS(2,1)**2*RS(2,3)-59.E0/0.320166E7*RS(2,1)*RS(2,4)*
     #RS(2,2)+59.E0/0.320166E7*RS(2,3)*RS(2,4)*RS(2,2)-137.E0/0.480249E7
     #*RS(2,3)**2*RS(2,4)-881.E0/0.160083E8*RS(2,2)**2*RS(2,3)+137.E0/0.
     #480249E7*RS(2,1)**2*RS(2,4)+137.E0/0.480249E7*RS(2,1)**2*RS(2,2)+8
     #81.E0/0.160083E8*RS(2,1)*RS(2,4)**2+881.E0/0.160083E8*RS(2,1)*RS(2
     #,2)**2-881.E0/0.160083E8*RS(2,3)*RS(2,4)**2+163.E0/0.320166E8*RS(2
     #,1)*RS(2,3)**2-137.E0/0.480249E7*RS(2,3)**2*RS(2,2))*RS(1,4)**2
      s12 = s13+s14
      s10 = s11+s12
      s11 = RS(1,2)
      s9 = s10*s11
      s10 = (-17.E0/381150.E0*RS(2,3)*RS(2,4)**2+881.E0/0.160083E8*RS(2,
     #1)**2*RS(2,3)-RS(2,2)**2*RS(2,3)/210-11357.E0/0.160083E8*RS(2,1)*R
     #S(2,3)**2-RS(2,2)**3/420-541.E0/0.960498E8*RS(2,1)*RS(2,4)**2-13.E
     #0/142296.E0*RS(2,2)**2*RS(2,4)-59.E0/0.1280664E8*RS(2,1)**3-RS(2,3
     #)**2*RS(2,2)/210+37517.E0/0.480249E8*RS(2,3)**2*RS(2,4)+19.E0/8893
     #50.E0*RS(2,4)**2*RS(2,2)+337.E0/0.1920996E9*RS(2,1)**2*RS(2,4)-247
     #.E0/0.4002075E7*RS(2,1)*RS(2,3)*RS(2,4)-13.E0/142296.E0*RS(2,1)*RS
     #(2,2)**2+137.E0/0.960498E7*RS(2,1)**2*RS(2,2)-13.E0/142296.E0*RS(2
     #,4)**3-499.E0/0.137214E7*RS(2,1)*RS(2,3)*RS(2,2)-17.E0/381150.E0*R
     #S(2,3)*RS(2,4)*RS(2,2)+106.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2))
     #*RS(1,3)**3
      s8 = s9+s10
      s6 = s7+s8
      s7 = s6+(137.E0/0.480249E7*RS(2,4)**2*RS(2,2)-541.E0/0.960498E8*RS
     #(2,1)*RS(2,2)**2-169.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-1303.E
     #0/0.2401245E8*RS(2,1)*RS(2,3)*RS(2,4)+19.E0/889350.E0*RS(2,3)*RS(2
     #,4)*RS(2,2)+107.E0/0.960498E8*RS(2,1)*RS(2,4)*RS(2,2)-37517.E0/0.4
     #80249E8*RS(2,3)**3+149.E0/0.2401245E8*RS(2,1)**2*RS(2,4)+13.E0/142
     #296.E0*RS(2,2)**2*RS(2,3)-13.E0/47432.E0*RS(2,3)*RS(2,4)**2-191.E0
     #/0.960498E8*RS(2,1)**2*RS(2,2)+337.E0/0.1920996E9*RS(2,1)**3-499.E
     #0/0.137214E7*RS(2,4)**3-1721.E0/0.160083E8*RS(2,1)*RS(2,3)**2+137.
     #E0/0.960498E7*RS(2,2)**2*RS(2,4)+901.E0/0.640332E8*RS(2,1)**2*RS(2
     #,3)+17.E0/381150.E0*RS(2,3)**2*RS(2,2)-631.E0/0.320166E8*RS(2,1)*R
     #S(2,4)**2)*RS(1,4)*RS(1,3)**2
      s8 = s7
      s10 = (149.E0/0.480249E8*RS(2,1)**2*RS(2,3)-191.E0/0.960498E8*RS(2
     #,1)**2*RS(2,2)-3589.E0/0.960498E8*RS(2,1)*RS(2,3)**2+17.E0/381150.
     #E0*RS(2,3)**3-137.E0/0.960498E7*RS(2,2)**2*RS(2,3)-19.E0/444675.E0
     #*RS(2,3)**2*RS(2,2)-11357.E0/0.160083E8*RS(2,4)**3-59.E0/0.1280664
     #E8*RS(2,2)**2*RS(2,4)+13.E0/47432.E0*RS(2,3)**2*RS(2,4)+881.E0/0.1
     #60083E8*RS(2,4)**2*RS(2,2)-541.E0/0.960498E8*RS(2,1)**3+787.E0/0.9
     #60498E8*RS(2,1)*RS(2,3)*RS(2,2)-101.E0/0.160083E7*RS(2,1)*RS(2,3)*
     #RS(2,4)+38.E0/0.2401245E7*RS(2,1)*RS(2,4)*RS(2,2)-631.E0/0.320166E
     #8*RS(2,1)**2*RS(2,4)-137.E0/0.960498E7*RS(2,3)*RS(2,4)*RS(2,2)+337
     #.E0/0.1920996E9*RS(2,1)*RS(2,2)**2-43.E0/254100.E0*RS(2,1)*RS(2,4)
     #**2)*RS(1,4)**2*RS(1,3)
      s11 = (-137.E0/0.960498E7*RS(2,3)**2*RS(2,2)+137.E0/0.960498E7*RS(
     #2,1)**2*RS(2,2)-463.E0/0.320166E8*RS(2,1)*RS(2,3)**2-881.E0/0.1600
     #83E8*RS(2,3)*RS(2,4)*RS(2,2)+881.E0/0.160083E8*RS(2,1)*RS(2,4)*RS(
     #2,2)+499.E0/0.137214E7*RS(2,3)**2*RS(2,4)-59.E0/0.1280664E8*RS(2,1
     #)*RS(2,2)**2+463.E0/0.320166E8*RS(2,1)**2*RS(2,3)+11357.E0/0.16008
     #3E8*RS(2,3)*RS(2,4)**2+59.E0/0.1280664E8*RS(2,2)**2*RS(2,3)+13.E0/
     #142296.E0*RS(2,3)**3-13.E0/142296.E0*RS(2,1)**3-11357.E0/0.160083E
     #8*RS(2,1)*RS(2,4)**2-499.E0/0.137214E7*RS(2,1)**2*RS(2,4))*RS(1,4)
     #**3
      s9 = s10+s11
      s4 = s8+s9
	ARRSS(2) = S4
C      s5 = phi(2)
C      s3 = s4*s5
C      s1 = s2+s3
C      s3 = s1

      s8 = (499.E0/0.137214E7*RS(2,1)*RS(2,4)**2+463.E0/0.320166E8*RS(2,
     #2)**2*RS(2,4)-11357.E0/0.160083E8*RS(2,1)**2*RS(2,2)-13.E0/142296.
     #E0*RS(2,2)**3-463.E0/0.320166E8*RS(2,4)**2*RS(2,2)-137.E0/0.960498
     #E7*RS(2,3)*RS(2,4)**2+13.E0/142296.E0*RS(2,4)**3-59.E0/0.1280664E8
     #*RS(2,3)**2*RS(2,2)-499.E0/0.137214E7*RS(2,1)*RS(2,2)**2+137.E0/0.
     #960498E7*RS(2,2)**2*RS(2,3)+11357.E0/0.160083E8*RS(2,1)**2*RS(2,4)
     #+59.E0/0.1280664E8*RS(2,3)**2*RS(2,4)+881.E0/0.160083E8*RS(2,1)*RS
     #(2,3)*RS(2,2)-881.E0/0.160083E8*RS(2,1)*RS(2,3)*RS(2,4))*RS(1,1)**
     #3
      s11 = (59.E0/0.1280664E8*RS(2,1)*RS(2,3)**2+19.E0/444675.E0*RS(2,2
     #)**2*RS(2,3)-13.E0/47432.E0*RS(2,1)*RS(2,2)**2+631.E0/0.320166E8*R
     #S(2,1)*RS(2,4)**2+3589.E0/0.960498E8*RS(2,2)**2*RS(2,4)+541.E0/0.9
     #60498E8*RS(2,4)**3+191.E0/0.960498E8*RS(2,3)*RS(2,4)**2-337.E0/0.1
     #920996E9*RS(2,3)**2*RS(2,4)-38.E0/0.2401245E7*RS(2,1)*RS(2,3)*RS(2
     #,4)+43.E0/254100.E0*RS(2,1)**2*RS(2,4)-881.E0/0.160083E8*RS(2,1)**
     #2*RS(2,3)+137.E0/0.960498E7*RS(2,1)*RS(2,3)*RS(2,2)-17.E0/381150.E
     #0*RS(2,2)**3+101.E0/0.160083E7*RS(2,1)*RS(2,4)*RS(2,2)+11357.E0/0.
     #160083E8*RS(2,1)**3-787.E0/0.960498E8*RS(2,3)*RS(2,4)*RS(2,2)+137.
     #E0/0.960498E7*RS(2,3)**2*RS(2,2)-149.E0/0.480249E8*RS(2,4)**2*RS(2
     #,2))*RS(1,2)
      s13 = (-19.E0/889350.E0*RS(2,4)**3+19.E0/889350.E0*RS(2,2)**3+137.
     #E0/0.480249E7*RS(2,1)*RS(2,2)**2-163.E0/0.320166E8*RS(2,2)**2*RS(2
     #,4)-881.E0/0.160083E8*RS(2,3)**2*RS(2,4)+881.E0/0.160083E8*RS(2,3)
     #**2*RS(2,2)+881.E0/0.160083E8*RS(2,1)**2*RS(2,2)+163.E0/0.320166E8
     #*RS(2,4)**2*RS(2,2)-137.E0/0.480249E7*RS(2,1)*RS(2,4)**2-59.E0/0.3
     #20166E7*RS(2,1)*RS(2,3)*RS(2,2)-137.E0/0.480249E7*RS(2,3)*RS(2,4)*
     #*2+59.E0/0.320166E7*RS(2,1)*RS(2,3)*RS(2,4)-881.E0/0.160083E8*RS(2
     #,1)**2*RS(2,4)+137.E0/0.480249E7*RS(2,2)**2*RS(2,3))*RS(1,3)
      s14 = (-541.E0/0.960498E8*RS(2,2)**3-19.E0/444675.E0*RS(2,3)*RS(2,
     #4)**2+38.E0/0.2401245E7*RS(2,1)*RS(2,3)*RS(2,2)+17.E0/381150.E0*RS
     #(2,4)**3-11357.E0/0.160083E8*RS(2,1)**3-191.E0/0.960498E8*RS(2,2)*
     #*2*RS(2,3)+149.E0/0.480249E8*RS(2,2)**2*RS(2,4)-3589.E0/0.960498E8
     #*RS(2,4)**2*RS(2,2)+337.E0/0.1920996E9*RS(2,3)**2*RS(2,2)+787.E0/0
     #.960498E8*RS(2,3)*RS(2,4)*RS(2,2)-631.E0/0.320166E8*RS(2,1)*RS(2,2
     #)**2+13.E0/47432.E0*RS(2,1)*RS(2,4)**2-137.E0/0.960498E7*RS(2,3)**
     #2*RS(2,4)-137.E0/0.960498E7*RS(2,1)*RS(2,3)*RS(2,4)-101.E0/0.16008
     #3E7*RS(2,1)*RS(2,4)*RS(2,2)+881.E0/0.160083E8*RS(2,1)**2*RS(2,3)-4
     #3.E0/254100.E0*RS(2,1)**2*RS(2,2)-59.E0/0.1280664E8*RS(2,1)*RS(2,3
     #)**2)*RS(1,4)
      s12 = s13+s14
      s10 = s11+s12
      s11 = RS(1,1)**2
      s9 = s10*s11
      s7 = s8+s9
      s9 = s7
      s13 = (541.E0/0.960498E8*RS(2,3)**2*RS(2,4)+37517.E0/0.480249E8*RS
     #(2,2)**3+1303.E0/0.2401245E8*RS(2,1)*RS(2,4)*RS(2,2)-137.E0/0.9604
     #98E7*RS(2,1)*RS(2,3)**2+13.E0/47432.E0*RS(2,1)**2*RS(2,2)-13.E0/14
     #2296.E0*RS(2,3)**2*RS(2,2)-337.E0/0.1920996E9*RS(2,4)**3-107.E0/0.
     #960498E8*RS(2,1)*RS(2,3)*RS(2,4)+499.E0/0.137214E7*RS(2,1)**3+191.
     #E0/0.960498E8*RS(2,3)*RS(2,4)**2+631.E0/0.320166E8*RS(2,1)**2*RS(2
     #,4)-137.E0/0.480249E7*RS(2,1)**2*RS(2,3)-901.E0/0.640332E8*RS(2,4)
     #**2*RS(2,2)-19.E0/889350.E0*RS(2,1)*RS(2,3)*RS(2,2)+169.E0/0.12006
     #22E8*RS(2,3)*RS(2,4)*RS(2,2)+1721.E0/0.160083E8*RS(2,2)**2*RS(2,4)
     #-149.E0/0.2401245E8*RS(2,1)*RS(2,4)**2-17.E0/381150.E0*RS(2,2)**2*
     #RS(2,3))*RS(1,2)**2
      s16 = (-107.E0/0.960498E8*RS(2,1)*RS(2,4)**2+19.E0/889350.E0*RS(2,
     #1)*RS(2,2)**2-38.E0/0.2401245E7*RS(2,1)**2*RS(2,4)-13.E0/35574.E0*
     #RS(2,2)**2*RS(2,3)-17.E0/381150.E0*RS(2,2)**3+247.E0/0.4002075E7*R
     #S(2,3)**2*RS(2,4)-499.E0/0.137214E7*RS(2,3)**2*RS(2,2)+RS(2,2)**2*
     #RS(2,4)/43659-106.E0/0.1200622E8*RS(2,4)**3-881.E0/0.160083E8*RS(2
     #,1)*RS(2,3)**2+191.E0/0.960498E8*RS(2,4)**2*RS(2,2)+59.E0/0.320166
     #E7*RS(2,1)**2*RS(2,3)-881.E0/0.160083E8*RS(2,1)**3-212.E0/0.120062
     #2E8*RS(2,3)*RS(2,4)**2-149.E0/0.2401245E8*RS(2,1)*RS(2,4)*RS(2,2)+
     #193.E0/0.480249E7*RS(2,3)*RS(2,4)*RS(2,2)-337.E0/0.480249E8*RS(2,1
     #)*RS(2,3)*RS(2,4)-137.E0/0.960498E7*RS(2,1)**2*RS(2,2))*RS(1,3)
      s17 = (-247.E0/0.4002075E7*RS(2,2)**3+149.E0/0.2401245E8*RS(2,1)*R
     #S(2,3)*RS(2,2)+47.E0/0.53361E7*RS(2,4)**2*RS(2,2)+247.E0/0.4002075
     #E7*RS(2,4)**3+212.E0/0.1200622E8*RS(2,2)**2*RS(2,3)-47.E0/0.53361E
     #7*RS(2,2)**2*RS(2,4)+1303.E0/0.2401245E8*RS(2,1)*RS(2,4)**2+106.E0
     #/0.1200622E8*RS(2,3)**2*RS(2,2)-101.E0/0.160083E7*RS(2,1)**2*RS(2,
     #2)-1303.E0/0.2401245E8*RS(2,1)*RS(2,2)**2-149.E0/0.2401245E8*RS(2,
     #1)*RS(2,3)*RS(2,4)+101.E0/0.160083E7*RS(2,1)**2*RS(2,4)-212.E0/0.1
     #200622E8*RS(2,3)*RS(2,4)**2-106.E0/0.1200622E8*RS(2,3)**2*RS(2,4))
     #*RS(1,4)
      s15 = s16+s17
      s16 = RS(1,2)
      s14 = s15*s16
      s12 = s13+s14
      s13 = s12+(59.E0/0.1280664E8*RS(2,1)**2*RS(2,4)-137.E0/0.960498E7*
     #RS(2,1)*RS(2,4)**2+463.E0/0.320166E8*RS(2,2)**2*RS(2,4)+137.E0/0.9
     #60498E7*RS(2,1)*RS(2,2)**2-13.E0/142296.E0*RS(2,2)**3-11357.E0/0.1
     #60083E8*RS(2,3)**2*RS(2,2)+881.E0/0.160083E8*RS(2,1)*RS(2,3)*RS(2,
     #2)-59.E0/0.1280664E8*RS(2,1)**2*RS(2,2)+11357.E0/0.160083E8*RS(2,3
     #)**2*RS(2,4)+499.E0/0.137214E7*RS(2,3)*RS(2,4)**2-463.E0/0.320166E
     #8*RS(2,4)**2*RS(2,2)+13.E0/142296.E0*RS(2,4)**3-499.E0/0.137214E7*
     #RS(2,2)**2*RS(2,3)-881.E0/0.160083E8*RS(2,1)*RS(2,3)*RS(2,4))*RS(1
     #,3)**2
      s14 = s13
      s16 = (106.E0/0.1200622E8*RS(2,2)**3+499.E0/0.137214E7*RS(2,3)**2*
     #RS(2,4)-59.E0/0.320166E7*RS(2,1)**2*RS(2,3)+17.E0/381150.E0*RS(2,4
     #)**3+137.E0/0.960498E7*RS(2,1)**2*RS(2,4)+881.E0/0.160083E8*RS(2,1
     #)*RS(2,3)**2-RS(2,4)**2*RS(2,2)/43659+149.E0/0.2401245E8*RS(2,1)*R
     #S(2,4)*RS(2,2)-191.E0/0.960498E8*RS(2,2)**2*RS(2,4)+337.E0/0.48024
     #9E8*RS(2,1)*RS(2,3)*RS(2,2)+881.E0/0.160083E8*RS(2,1)**3+107.E0/0.
     #960498E8*RS(2,1)*RS(2,2)**2+38.E0/0.2401245E7*RS(2,1)**2*RS(2,2)+1
     #3.E0/35574.E0*RS(2,3)*RS(2,4)**2-19.E0/889350.E0*RS(2,1)*RS(2,4)**
     #2-193.E0/0.480249E7*RS(2,3)*RS(2,4)*RS(2,2)+212.E0/0.1200622E8*RS(
     #2,2)**2*RS(2,3)-247.E0/0.4002075E7*RS(2,3)**2*RS(2,2))*RS(1,4)*RS(
     #1,3)
      s17 = (337.E0/0.1920996E9*RS(2,2)**3-37517.E0/0.480249E8*RS(2,4)**
     #3-499.E0/0.137214E7*RS(2,1)**3-1303.E0/0.2401245E8*RS(2,1)*RS(2,4)
     #*RS(2,2)+107.E0/0.960498E8*RS(2,1)*RS(2,3)*RS(2,2)-169.E0/0.120062
     #2E8*RS(2,3)*RS(2,4)*RS(2,2)+137.E0/0.960498E7*RS(2,1)*RS(2,3)**2-1
     #91.E0/0.960498E8*RS(2,2)**2*RS(2,3)-541.E0/0.960498E8*RS(2,3)**2*R
     #S(2,2)+13.E0/142296.E0*RS(2,3)**2*RS(2,4)+17.E0/381150.E0*RS(2,3)*
     #RS(2,4)**2-13.E0/47432.E0*RS(2,1)**2*RS(2,4)+149.E0/0.2401245E8*RS
     #(2,1)*RS(2,2)**2-1721.E0/0.160083E8*RS(2,4)**2*RS(2,2)+137.E0/0.48
     #0249E7*RS(2,1)**2*RS(2,3)+901.E0/0.640332E8*RS(2,2)**2*RS(2,4)-631
     #.E0/0.320166E8*RS(2,1)**2*RS(2,2)+19.E0/889350.E0*RS(2,1)*RS(2,3)*
     #RS(2,4))*RS(1,4)**2
      s15 = s16+s17
      s11 = s14+s15
      s12 = RS(1,1)
      s10 = s11*s12
      s8 = s9+s10
      s9 = s8
      s11 = (499.E0/0.137214E7*RS(2,3)*RS(2,4)*RS(2,2)+59.E0/0.1280664E8
     #*RS(2,4)**3+247.E0/0.4002075E7*RS(2,1)*RS(2,4)*RS(2,2)+11357.E0/0.
     #160083E8*RS(2,2)**2*RS(2,4)+RS(2,3)**3/420+541.E0/0.960498E8*RS(2,
     #1)**2*RS(2,4)+17.E0/381150.E0*RS(2,1)**2*RS(2,2)+13.E0/142296.E0*R
     #S(2,1)*RS(2,3)**2+17.E0/381150.E0*RS(2,1)*RS(2,3)*RS(2,2)+13.E0/14
     #2296.E0*RS(2,3)**2*RS(2,4)+13.E0/142296.E0*RS(2,1)**3-106.E0/0.120
     #0622E8*RS(2,1)*RS(2,3)*RS(2,4)-137.E0/0.960498E7*RS(2,3)*RS(2,4)**
     #2+RS(2,2)**2*RS(2,3)/210-37517.E0/0.480249E8*RS(2,1)*RS(2,2)**2-88
     #1.E0/0.160083E8*RS(2,4)**2*RS(2,2)-337.E0/0.1920996E9*RS(2,1)*RS(2
     #,4)**2+RS(2,3)**2*RS(2,2)/210-19.E0/889350.E0*RS(2,1)**2*RS(2,3))*
     #RS(1,2)**3
      s14 = (-RS(2,2)**3/210+191.E0/0.960498E8*RS(2,1)**2*RS(2,4)-137.E0
     #/0.960498E7*RS(2,4)**3-19.E0/444675.E0*RS(2,3)*RS(2,4)**2+17.E0/38
     #1150.E0*RS(2,3)**2*RS(2,4)+499.E0/0.137214E7*RS(2,2)**2*RS(2,4)-19
     #.E0/444675.E0*RS(2,1)**2*RS(2,2)+499.E0/0.137214E7*RS(2,1)*RS(2,3)
     #**2+191.E0/0.960498E8*RS(2,1)*RS(2,4)**2-137.E0/0.480249E7*RS(2,1)
     #**2*RS(2,3)-137.E0/0.480249E7*RS(2,4)**2*RS(2,2)+17.E0/381150.E0*R
     #S(2,1)*RS(2,2)**2+13.E0/35574.E0*RS(2,1)*RS(2,3)*RS(2,2)+RS(2,3)**
     #3/105-137.E0/0.960498E7*RS(2,1)**3+13.E0/35574.E0*RS(2,3)*RS(2,4)*
     #RS(2,2)-212.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4)+RS(2,3)**2*RS(2
     #,2)/140-212.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2))*RS(1,3)
      s15 = (-11357.E0/0.160083E8*RS(2,2)**3-137.E0/0.480249E7*RS(2,3)*R
     #S(2,4)**2-499.E0/0.137214E7*RS(2,2)**2*RS(2,3)-881.E0/0.160083E8*R
     #S(2,4)**3-463.E0/0.320166E8*RS(2,1)**3+163.E0/0.320166E8*RS(2,1)**
     #2*RS(2,3)-1721.E0/0.160083E8*RS(2,1)*RS(2,2)**2+59.E0/0.426888E7*R
     #S(2,4)**2*RS(2,2)-13.E0/142296.E0*RS(2,3)**2*RS(2,2)+191.E0/0.9604
     #98E8*RS(2,1)*RS(2,3)*RS(2,4)-RS(2,1)*RS(2,3)*RS(2,2)/43659-901.E0/
     #0.640332E8*RS(2,1)*RS(2,4)**2-137.E0/0.960498E7*RS(2,3)*RS(2,4)*RS
     #(2,2)-149.E0/0.480249E8*RS(2,1)**2*RS(2,4)-3589.E0/0.960498E8*RS(2
     #,1)**2*RS(2,2)-463.E0/0.320166E8*RS(2,1)*RS(2,3)**2-19.E0/889350.E
     #0*RS(2,3)**2*RS(2,4)+47.E0/0.53361E7*RS(2,1)*RS(2,4)*RS(2,2))*RS(1
     #,4)
      s13 = s14+s15
      s14 = RS(1,2)**2
      s12 = s13*s14
      s10 = s11+s12
      s6 = s9+s10
      s8 = s6
      s12 = (13.E0/142296.E0*RS(2,2)**2*RS(2,4)-337.E0/0.1920996E9*RS(2,
     #1)**2*RS(2,4)-19.E0/889350.E0*RS(2,4)**2*RS(2,2)+59.E0/0.1280664E8
     #*RS(2,1)**3+541.E0/0.960498E8*RS(2,1)*RS(2,4)**2-137.E0/0.960498E7
     #*RS(2,1)**2*RS(2,2)+RS(2,3)**3/42-RS(2,2)**3/210+11357.E0/0.160083
     #E8*RS(2,1)*RS(2,3)**2-881.E0/0.160083E8*RS(2,1)**2*RS(2,3)+499.E0/
     #0.137214E7*RS(2,1)*RS(2,3)*RS(2,2)+13.E0/142296.E0*RS(2,4)**3-3751
     #7.E0/0.480249E8*RS(2,3)**2*RS(2,4)+17.E0/381150.E0*RS(2,3)*RS(2,4)
     #**2+17.E0/381150.E0*RS(2,3)*RS(2,4)*RS(2,2)-106.E0/0.1200622E8*RS(
     #2,1)*RS(2,4)*RS(2,2)-RS(2,2)**2*RS(2,3)/140+13.E0/142296.E0*RS(2,1
     #)*RS(2,2)**2+247.E0/0.4002075E7*RS(2,1)*RS(2,3)*RS(2,4))*RS(1,3)**
     #2
      s14 = (-499.E0/0.137214E7*RS(2,2)**3+193.E0/0.480249E7*RS(2,1)*RS(
     #2,3)*RS(2,4)-13.E0/35574.E0*RS(2,2)**2*RS(2,3)+137.E0/0.960498E7*R
     #S(2,2)**2*RS(2,4)-169.E0/0.1200622E8*RS(2,1)*RS(2,2)**2-787.E0/0.9
     #60498E8*RS(2,1)**2*RS(2,4)+499.E0/0.137214E7*RS(2,4)**3+169.E0/0.1
     #200622E8*RS(2,1)*RS(2,4)**2+13.E0/35574.E0*RS(2,3)*RS(2,4)**2+17.E
     #0/381150.E0*RS(2,3)**2*RS(2,4)-137.E0/0.960498E7*RS(2,4)**2*RS(2,2
     #)-193.E0/0.480249E7*RS(2,1)*RS(2,3)*RS(2,2)-17.E0/381150.E0*RS(2,3
     #)**2*RS(2,2)+787.E0/0.960498E8*RS(2,1)**2*RS(2,2))*RS(1,4)*RS(1,3)
      s15 = (881.E0/0.160083E8*RS(2,2)**3+11357.E0/0.160083E8*RS(2,4)**3
     #+463.E0/0.320166E8*RS(2,1)**3-47.E0/0.53361E7*RS(2,1)*RS(2,4)*RS(2
     #,2)-191.E0/0.960498E8*RS(2,1)*RS(2,3)*RS(2,2)+137.E0/0.960498E7*RS
     #(2,3)*RS(2,4)*RS(2,2)-59.E0/0.426888E7*RS(2,2)**2*RS(2,4)+137.E0/0
     #.480249E7*RS(2,2)**2*RS(2,3)+13.E0/142296.E0*RS(2,3)**2*RS(2,4)+49
     #9.E0/0.137214E7*RS(2,3)*RS(2,4)**2+3589.E0/0.960498E8*RS(2,1)**2*R
     #S(2,4)+901.E0/0.640332E8*RS(2,1)*RS(2,2)**2+19.E0/889350.E0*RS(2,3
     #)**2*RS(2,2)-163.E0/0.320166E8*RS(2,1)**2*RS(2,3)+463.E0/0.320166E
     #8*RS(2,1)*RS(2,3)**2+1721.E0/0.160083E8*RS(2,1)*RS(2,4)**2+149.E0/
     #0.480249E8*RS(2,1)**2*RS(2,2)+RS(2,1)*RS(2,3)*RS(2,4)/43659)*RS(1,
     #4)**2
      s13 = s14+s15
      s11 = s12+s13
      s12 = RS(1,2)
      s10 = s11*s12
      s11 = (RS(2,4)**3/420-RS(2,2)**2*RS(2,3)/105-RS(2,2)**3/420+RS(2,3
     #)*RS(2,4)**2/105-RS(2,3)**2*RS(2,2)/42+RS(2,3)**2*RS(2,4)/42)*RS(1
     #,3)**3
      s9 = s10+s11
      s7 = s8+s9
      s8 = s7+(-RS(2,3)**3/42-13.E0/142296.E0*RS(2,1)*RS(2,4)**2-247.E0/
     #0.4002075E7*RS(2,1)*RS(2,3)*RS(2,2)+106.E0/0.1200622E8*RS(2,1)*RS(
     #2,4)*RS(2,2)-17.E0/381150.E0*RS(2,3)*RS(2,4)*RS(2,2)-541.E0/0.9604
     #98E8*RS(2,1)*RS(2,2)**2-499.E0/0.137214E7*RS(2,1)*RS(2,3)*RS(2,4)+
     #RS(2,4)**3/210-59.E0/0.1280664E8*RS(2,1)**3+881.E0/0.160083E8*RS(2
     #,1)**2*RS(2,3)+337.E0/0.1920996E9*RS(2,1)**2*RS(2,2)+RS(2,3)*RS(2,
     #4)**2/140-17.E0/381150.E0*RS(2,2)**2*RS(2,3)+137.E0/0.960498E7*RS(
     #2,1)**2*RS(2,4)-13.E0/142296.E0*RS(2,2)**3+19.E0/889350.E0*RS(2,2)
     #**2*RS(2,4)-13.E0/142296.E0*RS(2,4)**2*RS(2,2)+37517.E0/0.480249E8
     #*RS(2,3)**2*RS(2,2)-11357.E0/0.160083E8*RS(2,1)*RS(2,3)**2)*RS(1,4
     #)*RS(1,3)**2
      s9 = s8
      s11 = (19.E0/444675.E0*RS(2,2)**2*RS(2,3)-13.E0/35574.E0*RS(2,1)*R
     #S(2,3)*RS(2,4)+137.E0/0.960498E7*RS(2,2)**3-RS(2,3)**3/105+212.E0/
     #0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-17.E0/381150.E0*RS(2,3)**2*RS(
     #2,2)+19.E0/444675.E0*RS(2,1)**2*RS(2,4)+RS(2,4)**3/210+137.E0/0.96
     #0498E7*RS(2,1)**3-RS(2,3)**2*RS(2,4)/140-499.E0/0.137214E7*RS(2,4)
     #**2*RS(2,2)+212.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2)+137.E0/0.48
     #0249E7*RS(2,1)**2*RS(2,3)-13.E0/35574.E0*RS(2,3)*RS(2,4)*RS(2,2)+1
     #37.E0/0.480249E7*RS(2,2)**2*RS(2,4)-17.E0/381150.E0*RS(2,1)*RS(2,4
     #)**2-499.E0/0.137214E7*RS(2,1)*RS(2,3)**2-191.E0/0.960498E8*RS(2,1
     #)*RS(2,2)**2-191.E0/0.960498E8*RS(2,1)**2*RS(2,2))*RS(1,4)**2*RS(1
     #,3)
      s12 = (106.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-499.E0/0.137214E
     #7*RS(2,3)*RS(2,4)*RS(2,2)-247.E0/0.4002075E7*RS(2,1)*RS(2,4)*RS(2,
     #2)+19.E0/889350.E0*RS(2,1)**2*RS(2,3)-RS(2,3)**2*RS(2,4)/210-11357
     #.E0/0.160083E8*RS(2,4)**2*RS(2,2)+337.E0/0.1920996E9*RS(2,1)*RS(2,
     #2)**2-RS(2,3)*RS(2,4)**2/210+137.E0/0.960498E7*RS(2,2)**2*RS(2,3)-
     #59.E0/0.1280664E8*RS(2,2)**3+881.E0/0.160083E8*RS(2,2)**2*RS(2,4)-
     #13.E0/142296.E0*RS(2,1)*RS(2,3)**2-13.E0/142296.E0*RS(2,3)**2*RS(2
     #,2)-541.E0/0.960498E8*RS(2,1)**2*RS(2,2)-13.E0/142296.E0*RS(2,1)**
     #3-17.E0/381150.E0*RS(2,1)**2*RS(2,4)+37517.E0/0.480249E8*RS(2,1)*R
     #S(2,4)**2-RS(2,3)**3/420-17.E0/381150.E0*RS(2,1)*RS(2,3)*RS(2,4))*
     #RS(1,4)**3
      s10 = s11+s12
      s5 = s9+s10
	ARRSS(3) = S5
C      s6 = phi(3)
C      s4 = s5*s6
C      s2 = s3+s4
C      s3 = s2

      s9 = (-17.E0/381150.E0*RS(2,1)*RS(2,2)**2-RS(2,1)**2*RS(2,4)/210-5
     #41.E0/0.960498E8*RS(2,2)**2*RS(2,3)+137.E0/0.960498E7*RS(2,3)**2*R
     #S(2,4)-59.E0/0.1280664E8*RS(2,3)**3-13.E0/142296.E0*RS(2,3)*RS(2,4
     #)**2-13.E0/142296.E0*RS(2,2)**3-RS(2,1)*RS(2,4)**2/210-11357.E0/0.
     #160083E8*RS(2,1)**2*RS(2,3)+337.E0/0.1920996E9*RS(2,3)**2*RS(2,2)+
     #881.E0/0.160083E8*RS(2,1)*RS(2,3)**2+19.E0/889350.E0*RS(2,2)**2*RS
     #(2,4)-13.E0/142296.E0*RS(2,4)**2*RS(2,2)+37517.E0/0.480249E8*RS(2,
     #1)**2*RS(2,2)-RS(2,4)**3/420+106.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS
     #(2,2)-247.E0/0.4002075E7*RS(2,1)*RS(2,3)*RS(2,2)-499.E0/0.137214E7
     #*RS(2,1)*RS(2,3)*RS(2,4)-17.E0/381150.E0*RS(2,1)*RS(2,4)*RS(2,2))*
     #RS(1,1)**3
      s12 = (901.E0/0.640332E8*RS(2,1)*RS(2,3)**2+337.E0/0.1920996E9*RS(
     #2,3)**3-541.E0/0.960498E8*RS(2,3)*RS(2,4)**2-191.E0/0.960498E8*RS(
     #2,3)**2*RS(2,4)+137.E0/0.480249E7*RS(2,2)**2*RS(2,4)-1721.E0/0.160
     #083E8*RS(2,1)**2*RS(2,3)-169.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4
     #)+19.E0/889350.E0*RS(2,1)*RS(2,4)*RS(2,2)+13.E0/142296.E0*RS(2,1)*
     #RS(2,4)**2-499.E0/0.137214E7*RS(2,2)**3-37517.E0/0.480249E8*RS(2,1
     #)**3-13.E0/47432.E0*RS(2,1)*RS(2,2)**2+149.E0/0.2401245E8*RS(2,3)*
     #*2*RS(2,2)-1303.E0/0.2401245E8*RS(2,1)*RS(2,3)*RS(2,2)+17.E0/38115
     #0.E0*RS(2,1)**2*RS(2,4)+107.E0/0.960498E8*RS(2,3)*RS(2,4)*RS(2,2)-
     #631.E0/0.320166E8*RS(2,2)**2*RS(2,3)+137.E0/0.960498E7*RS(2,4)**2*
     #RS(2,2))*RS(1,2)
      s14 = (463.E0/0.320166E8*RS(2,2)**3+881.E0/0.160083E8*RS(2,3)**3+3
     #589.E0/0.960498E8*RS(2,1)*RS(2,2)**2+13.E0/142296.E0*RS(2,1)*RS(2,
     #4)**2+149.E0/0.480249E8*RS(2,2)**2*RS(2,3)+901.E0/0.640332E8*RS(2,
     #3)**2*RS(2,2)+1721.E0/0.160083E8*RS(2,1)**2*RS(2,2)+463.E0/0.32016
     #6E8*RS(2,4)**2*RS(2,2)-163.E0/0.320166E8*RS(2,2)**2*RS(2,4)+11357.
     #E0/0.160083E8*RS(2,1)**3+19.E0/889350.E0*RS(2,3)*RS(2,4)**2-47.E0/
     #0.53361E7*RS(2,1)*RS(2,3)*RS(2,2)-59.E0/0.426888E7*RS(2,1)*RS(2,3)
     #**2+RS(2,1)*RS(2,4)*RS(2,2)/43659-191.E0/0.960498E8*RS(2,3)*RS(2,4
     #)*RS(2,2)+137.E0/0.960498E7*RS(2,1)*RS(2,3)*RS(2,4)+499.E0/0.13721
     #4E7*RS(2,1)**2*RS(2,4)+137.E0/0.480249E7*RS(2,3)**2*RS(2,4))*RS(1,
     #3)
      s15 = (137.E0/0.960498E7*RS(2,2)**3+137.E0/0.960498E7*RS(2,3)**3-R
     #S(2,4)**3/105+RS(2,1)**3/210-499.E0/0.137214E7*RS(2,4)**2*RS(2,2)-
     #191.E0/0.960498E8*RS(2,3)**2*RS(2,2)-191.E0/0.960498E8*RS(2,2)**2*
     #RS(2,3)+137.E0/0.480249E7*RS(2,2)**2*RS(2,4)-RS(2,1)*RS(2,4)**2/14
     #0+19.E0/444675.E0*RS(2,3)**2*RS(2,4)-13.E0/35574.E0*RS(2,1)*RS(2,4
     #)*RS(2,2)+19.E0/444675.E0*RS(2,1)*RS(2,2)**2+212.E0/0.1200622E8*RS
     #(2,3)*RS(2,4)*RS(2,2)-17.E0/381150.E0*RS(2,1)**2*RS(2,2)+212.E0/0.
     #1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-17.E0/381150.E0*RS(2,3)*RS(2,4)*
     #*2-499.E0/0.137214E7*RS(2,1)**2*RS(2,3)+137.E0/0.480249E7*RS(2,1)*
     #RS(2,3)**2-13.E0/35574.E0*RS(2,1)*RS(2,3)*RS(2,4))*RS(1,4)
      s13 = s14+s15
      s11 = s12+s13
      s12 = RS(1,1)**2
      s10 = s11*s12
      s8 = s9+s10
      s10 = s8
      s14 = (-19.E0/444675.E0*RS(2,1)**2*RS(2,4)+337.E0/0.1920996E9*RS(2
     #,3)*RS(2,4)**2-631.E0/0.320166E8*RS(2,3)**2*RS(2,2)-43.E0/254100.E
     #0*RS(2,2)**2*RS(2,3)-137.E0/0.960498E7*RS(2,1)*RS(2,4)*RS(2,2)+17.
     #E0/381150.E0*RS(2,1)**3-3589.E0/0.960498E8*RS(2,1)**2*RS(2,3)-59.E
     #0/0.1280664E8*RS(2,4)**2*RS(2,2)+149.E0/0.480249E8*RS(2,1)*RS(2,3)
     #**2-101.E0/0.160083E7*RS(2,1)*RS(2,3)*RS(2,2)-191.E0/0.960498E8*RS
     #(2,3)**2*RS(2,4)-137.E0/0.960498E7*RS(2,1)*RS(2,4)**2+881.E0/0.160
     #083E8*RS(2,2)**2*RS(2,4)-541.E0/0.960498E8*RS(2,3)**3+38.E0/0.2401
     #245E7*RS(2,3)*RS(2,4)*RS(2,2)-11357.E0/0.160083E8*RS(2,2)**3+787.E
     #0/0.960498E8*RS(2,1)*RS(2,3)*RS(2,4)+13.E0/47432.E0*RS(2,1)**2*RS(
     #2,2))*RS(1,2)**2
      s17 = (106.E0/0.1200622E8*RS(2,3)*RS(2,4)**2+101.E0/0.160083E7*RS(
     #2,1)*RS(2,2)**2-247.E0/0.4002075E7*RS(2,3)**3-1303.E0/0.2401245E8*
     #RS(2,3)**2*RS(2,2)+1303.E0/0.2401245E8*RS(2,1)**2*RS(2,2)+212.E0/0
     #.1200622E8*RS(2,3)**2*RS(2,4)+47.E0/0.53361E7*RS(2,1)**2*RS(2,3)-4
     #7.E0/0.53361E7*RS(2,1)*RS(2,3)**2+247.E0/0.4002075E7*RS(2,1)**3-10
     #6.E0/0.1200622E8*RS(2,1)*RS(2,4)**2-101.E0/0.160083E7*RS(2,2)**2*R
     #S(2,3)-149.E0/0.2401245E8*RS(2,1)*RS(2,4)*RS(2,2)+149.E0/0.2401245
     #E8*RS(2,3)*RS(2,4)*RS(2,2)-212.E0/0.1200622E8*RS(2,1)**2*RS(2,4))*
     #RS(1,3)
      s18 = (881.E0/0.160083E8*RS(2,2)**3+17.E0/381150.E0*RS(2,1)**3+149
     #.E0/0.2401245E8*RS(2,1)*RS(2,3)*RS(2,2)+106.E0/0.1200622E8*RS(2,3)
     #**3+38.E0/0.2401245E7*RS(2,2)**2*RS(2,3)-59.E0/0.320166E7*RS(2,2)*
     #*2*RS(2,4)+881.E0/0.160083E8*RS(2,4)**2*RS(2,2)+107.E0/0.960498E8*
     #RS(2,3)**2*RS(2,2)+137.E0/0.960498E7*RS(2,1)*RS(2,2)**2+499.E0/0.1
     #37214E7*RS(2,1)*RS(2,4)**2+13.E0/35574.E0*RS(2,1)**2*RS(2,4)-247.E
     #0/0.4002075E7*RS(2,3)*RS(2,4)**2-RS(2,1)**2*RS(2,3)/43659-19.E0/88
     #9350.E0*RS(2,1)**2*RS(2,2)+337.E0/0.480249E8*RS(2,3)*RS(2,4)*RS(2,
     #2)-193.E0/0.480249E7*RS(2,1)*RS(2,3)*RS(2,4)-191.E0/0.960498E8*RS(
     #2,1)*RS(2,3)**2+212.E0/0.1200622E8*RS(2,3)**2*RS(2,4))*RS(1,4)
      s16 = s17+s18
      s17 = RS(1,2)
      s15 = s16*s17
      s13 = s14+s15
      s14 = s13+(59.E0/0.426888E7*RS(2,1)**2*RS(2,3)+163.E0/0.320166E8*R
     #S(2,2)**2*RS(2,4)-463.E0/0.320166E8*RS(2,4)**2*RS(2,2)-137.E0/0.48
     #0249E7*RS(2,1)**2*RS(2,4)-19.E0/889350.E0*RS(2,1)*RS(2,4)**2-1721.
     #E0/0.160083E8*RS(2,3)**2*RS(2,2)-11357.E0/0.160083E8*RS(2,3)**3-46
     #3.E0/0.320166E8*RS(2,2)**3+47.E0/0.53361E7*RS(2,1)*RS(2,3)*RS(2,2)
     #-901.E0/0.640332E8*RS(2,1)**2*RS(2,2)-881.E0/0.160083E8*RS(2,1)**3
     #-499.E0/0.137214E7*RS(2,3)**2*RS(2,4)-13.E0/142296.E0*RS(2,3)*RS(2
     #,4)**2-RS(2,3)*RS(2,4)*RS(2,2)/43659+191.E0/0.960498E8*RS(2,1)*RS(
     #2,4)*RS(2,2)-3589.E0/0.960498E8*RS(2,2)**2*RS(2,3)-149.E0/0.480249
     #E8*RS(2,1)*RS(2,2)**2-137.E0/0.960498E7*RS(2,1)*RS(2,3)*RS(2,4))*R
     #S(1,3)**2
      s15 = s14
      s17 = (193.E0/0.480249E7*RS(2,1)*RS(2,4)*RS(2,2)-169.E0/0.1200622E
     #8*RS(2,3)**2*RS(2,2)+13.E0/35574.E0*RS(2,1)**2*RS(2,4)-137.E0/0.96
     #0498E7*RS(2,1)**2*RS(2,3)-193.E0/0.480249E7*RS(2,3)*RS(2,4)*RS(2,2
     #)+137.E0/0.960498E7*RS(2,1)*RS(2,3)**2-13.E0/35574.E0*RS(2,3)**2*R
     #S(2,4)-787.E0/0.960498E8*RS(2,1)*RS(2,2)**2+169.E0/0.1200622E8*RS(
     #2,1)**2*RS(2,2)-499.E0/0.137214E7*RS(2,3)**3-17.E0/381150.E0*RS(2,
     #3)*RS(2,4)**2+17.E0/381150.E0*RS(2,1)*RS(2,4)**2+787.E0/0.960498E8
     #*RS(2,2)**2*RS(2,3)+499.E0/0.137214E7*RS(2,1)**3)*RS(1,4)*RS(1,3)
      s18 = (-13.E0/142296.E0*RS(2,3)**3-59.E0/0.1280664E8*RS(2,2)**3-RS
     #(2,4)**3/42+RS(2,1)**3/210-247.E0/0.4002075E7*RS(2,3)*RS(2,4)*RS(2
     #,2)+881.E0/0.160083E8*RS(2,2)**2*RS(2,4)-541.E0/0.960498E8*RS(2,3)
     #**2*RS(2,2)+337.E0/0.1920996E9*RS(2,2)**2*RS(2,3)+37517.E0/0.48024
     #9E8*RS(2,3)*RS(2,4)**2+RS(2,1)**2*RS(2,4)/140-13.E0/142296.E0*RS(2
     #,1)**2*RS(2,2)-11357.E0/0.160083E8*RS(2,4)**2*RS(2,2)+137.E0/0.960
     #498E7*RS(2,1)*RS(2,2)**2+19.E0/889350.E0*RS(2,1)*RS(2,3)**2+106.E0
     #/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-13.E0/142296.E0*RS(2,1)**2*RS
     #(2,3)-17.E0/381150.E0*RS(2,1)*RS(2,3)*RS(2,4)-17.E0/381150.E0*RS(2
     #,3)**2*RS(2,4)-499.E0/0.137214E7*RS(2,1)*RS(2,4)*RS(2,2))*RS(1,4)*
     #*2
      s16 = s17+s18
      s12 = s15+s16
      s13 = RS(1,1)
      s11 = s12*s13
      s9 = s10+s11
      s10 = s9
      s12 = (137.E0/0.960498E7*RS(2,3)**2*RS(2,4)-11357.E0/0.160083E8*RS
     #(2,2)**2*RS(2,3)-59.E0/0.1280664E8*RS(2,3)*RS(2,4)**2+499.E0/0.137
     #214E7*RS(2,1)**2*RS(2,2)+881.E0/0.160083E8*RS(2,3)*RS(2,4)*RS(2,2)
     #-137.E0/0.960498E7*RS(2,1)**2*RS(2,4)+13.E0/142296.E0*RS(2,1)**3+1
     #1357.E0/0.160083E8*RS(2,1)*RS(2,2)**2+463.E0/0.320166E8*RS(2,1)*RS
     #(2,3)**2-499.E0/0.137214E7*RS(2,3)**2*RS(2,2)-881.E0/0.160083E8*RS
     #(2,1)*RS(2,4)*RS(2,2)-13.E0/142296.E0*RS(2,3)**3+59.E0/0.1280664E8
     #*RS(2,1)*RS(2,4)**2-463.E0/0.320166E8*RS(2,1)**2*RS(2,3))*RS(1,2)*
     #*3
      s15 = (-881.E0/0.160083E8*RS(2,2)**2*RS(2,4)-337.E0/0.1920996E9*RS
     #(2,1)*RS(2,4)**2-17.E0/381150.E0*RS(2,3)**3+137.E0/0.960498E7*RS(2
     #,3)*RS(2,4)**2+191.E0/0.960498E8*RS(2,1)**2*RS(2,4)-13.E0/47432.E0
     #*RS(2,3)**2*RS(2,2)+631.E0/0.320166E8*RS(2,1)**2*RS(2,2)+11357.E0/
     #0.160083E8*RS(2,2)**3+59.E0/0.1280664E8*RS(2,4)**2*RS(2,2)-149.E0/
     #0.480249E8*RS(2,1)**2*RS(2,3)+3589.E0/0.960498E8*RS(2,1)*RS(2,3)**
     #2+541.E0/0.960498E8*RS(2,1)**3+43.E0/254100.E0*RS(2,1)*RS(2,2)**2+
     #101.E0/0.160083E7*RS(2,1)*RS(2,3)*RS(2,2)-38.E0/0.2401245E7*RS(2,1
     #)*RS(2,4)*RS(2,2)+137.E0/0.960498E7*RS(2,3)*RS(2,4)*RS(2,2)-787.E0
     #/0.960498E8*RS(2,1)*RS(2,3)*RS(2,4)+19.E0/444675.E0*RS(2,3)**2*RS(
     #2,4))*RS(1,3)
      s16 = (-59.E0/0.320166E7*RS(2,3)*RS(2,4)*RS(2,2)-19.E0/889350.E0*R
     #S(2,1)**3+163.E0/0.320166E8*RS(2,1)**2*RS(2,3)+19.E0/889350.E0*RS(
     #2,3)**3+881.E0/0.160083E8*RS(2,2)**2*RS(2,3)+137.E0/0.480249E7*RS(
     #2,3)**2*RS(2,2)+59.E0/0.320166E7*RS(2,1)*RS(2,4)*RS(2,2)-881.E0/0.
     #160083E8*RS(2,1)*RS(2,2)**2-881.E0/0.160083E8*RS(2,1)*RS(2,4)**2-1
     #37.E0/0.480249E7*RS(2,1)**2*RS(2,4)+881.E0/0.160083E8*RS(2,3)*RS(2
     #,4)**2-137.E0/0.480249E7*RS(2,1)**2*RS(2,2)-163.E0/0.320166E8*RS(2
     #,1)*RS(2,3)**2+137.E0/0.480249E7*RS(2,3)**2*RS(2,4))*RS(1,4)
      s14 = s15+s16
      s15 = RS(1,2)**2
      s13 = s14*s15
      s11 = s12+s13
      s7 = s10+s11
      s9 = s7
      s13 = (541.E0/0.960498E8*RS(2,1)*RS(2,4)**2-137.E0/0.480249E7*RS(2
     #,2)**2*RS(2,4)+1721.E0/0.160083E8*RS(2,1)*RS(2,3)**2-13.E0/142296.
     #E0*RS(2,3)*RS(2,4)**2-337.E0/0.1920996E9*RS(2,1)**3+37517.E0/0.480
     #249E8*RS(2,3)**3-149.E0/0.2401245E8*RS(2,1)**2*RS(2,2)+191.E0/0.96
     #0498E8*RS(2,1)**2*RS(2,4)+1303.E0/0.2401245E8*RS(2,1)*RS(2,3)*RS(2
     #,2)+499.E0/0.137214E7*RS(2,2)**3-901.E0/0.640332E8*RS(2,1)**2*RS(2
     #,3)-17.E0/381150.E0*RS(2,3)**2*RS(2,4)-19.E0/889350.E0*RS(2,3)*RS(
     #2,4)*RS(2,2)-107.E0/0.960498E8*RS(2,1)*RS(2,4)*RS(2,2)-137.E0/0.96
     #0498E7*RS(2,4)**2*RS(2,2)+13.E0/47432.E0*RS(2,2)**2*RS(2,3)+169.E0
     #/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4)+631.E0/0.320166E8*RS(2,1)*RS(
     #2,2)**2)*RS(1,3)**2
      s15 = (-881.E0/0.160083E8*RS(2,4)**2*RS(2,2)+RS(2,1)*RS(2,3)**2/43
     #659+193.E0/0.480249E7*RS(2,1)*RS(2,3)*RS(2,4)-17.E0/381150.E0*RS(2
     #,3)**3-106.E0/0.1200622E8*RS(2,1)**3+19.E0/889350.E0*RS(2,3)**2*RS
     #(2,2)-499.E0/0.137214E7*RS(2,3)*RS(2,4)**2-107.E0/0.960498E8*RS(2,
     #1)**2*RS(2,2)+191.E0/0.960498E8*RS(2,1)**2*RS(2,3)-212.E0/0.120062
     #2E8*RS(2,1)**2*RS(2,4)-881.E0/0.160083E8*RS(2,2)**3-13.E0/35574.E0
     #*RS(2,3)**2*RS(2,4)-38.E0/0.2401245E7*RS(2,1)*RS(2,2)**2+59.E0/0.3
     #20166E7*RS(2,2)**2*RS(2,4)-149.E0/0.2401245E8*RS(2,1)*RS(2,3)*RS(2
     #,2)-337.E0/0.480249E8*RS(2,1)*RS(2,4)*RS(2,2)+247.E0/0.4002075E7*R
     #S(2,1)*RS(2,4)**2-137.E0/0.960498E7*RS(2,2)**2*RS(2,3))*RS(1,4)*RS
     #(1,3)
      s16 = (11357.E0/0.160083E8*RS(2,1)*RS(2,4)**2-13.E0/142296.E0*RS(2
     #,3)**3+881.E0/0.160083E8*RS(2,3)*RS(2,4)*RS(2,2)+13.E0/142296.E0*R
     #S(2,1)**3+137.E0/0.960498E7*RS(2,3)**2*RS(2,2)-499.E0/0.137214E7*R
     #S(2,3)**2*RS(2,4)-11357.E0/0.160083E8*RS(2,3)*RS(2,4)**2-59.E0/0.1
     #280664E8*RS(2,2)**2*RS(2,3)-137.E0/0.960498E7*RS(2,1)**2*RS(2,2)+5
     #9.E0/0.1280664E8*RS(2,1)*RS(2,2)**2+463.E0/0.320166E8*RS(2,1)*RS(2
     #,3)**2-463.E0/0.320166E8*RS(2,1)**2*RS(2,3)+499.E0/0.137214E7*RS(2
     #,1)**2*RS(2,4)-881.E0/0.160083E8*RS(2,1)*RS(2,4)*RS(2,2))*RS(1,4)*
     #*2
      s14 = s15+s16
      s12 = s13+s14
      s13 = RS(1,2)
      s11 = s12*s13
      s12 = (RS(2,3)*RS(2,4)**2/210+247.E0/0.4002075E7*RS(2,1)*RS(2,3)*R
     #S(2,2)-881.E0/0.160083E8*RS(2,1)**2*RS(2,3)-19.E0/889350.E0*RS(2,2
     #)**2*RS(2,4)-137.E0/0.960498E7*RS(2,1)**2*RS(2,4)+13.E0/142296.E0*
     #RS(2,1)*RS(2,4)**2+11357.E0/0.160083E8*RS(2,1)*RS(2,3)**2+59.E0/0.
     #1280664E8*RS(2,1)**3+499.E0/0.137214E7*RS(2,1)*RS(2,3)*RS(2,4)+RS(
     #2,4)**3/420+13.E0/142296.E0*RS(2,4)**2*RS(2,2)+RS(2,3)**2*RS(2,4)/
     #210+17.E0/381150.E0*RS(2,2)**2*RS(2,3)+541.E0/0.960498E8*RS(2,1)*R
     #S(2,2)**2-37517.E0/0.480249E8*RS(2,3)**2*RS(2,2)+13.E0/142296.E0*R
     #S(2,2)**3-337.E0/0.1920996E9*RS(2,1)**2*RS(2,2)+17.E0/381150.E0*RS
     #(2,3)*RS(2,4)*RS(2,2)-106.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2))*
     #RS(1,3)**3
      s10 = s11+s12
      s8 = s9+s10
      s9 = s8+(RS(2,3)*RS(2,4)**2/140+499.E0/0.137214E7*RS(2,1)*RS(2,3)*
     #*2-RS(2,3)**3/210+13.E0/35574.E0*RS(2,1)*RS(2,3)*RS(2,4)-212.E0/0.
     #1200622E8*RS(2,1)*RS(2,3)*RS(2,2)+RS(2,4)**3/105-19.E0/444675.E0*R
     #S(2,2)**2*RS(2,3)+191.E0/0.960498E8*RS(2,1)**2*RS(2,2)+13.E0/35574
     #.E0*RS(2,3)*RS(2,4)*RS(2,2)+191.E0/0.960498E8*RS(2,1)*RS(2,2)**2-1
     #37.E0/0.960498E7*RS(2,1)**3+17.E0/381150.E0*RS(2,1)*RS(2,4)**2-19.
     #E0/444675.E0*RS(2,1)**2*RS(2,4)-137.E0/0.480249E7*RS(2,2)**2*RS(2,
     #4)-212.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2)+17.E0/381150.E0*RS(2
     #,3)**2*RS(2,2)-137.E0/0.480249E7*RS(2,1)**2*RS(2,3)-137.E0/0.96049
     #8E7*RS(2,2)**3+499.E0/0.137214E7*RS(2,4)**2*RS(2,2))*RS(1,4)*RS(1,
     #3)**2
      s10 = s9
      s12 = (-137.E0/0.960498E7*RS(2,2)**2*RS(2,3)+11357.E0/0.160083E8*R
     #S(2,4)**2*RS(2,2)+59.E0/0.1280664E8*RS(2,2)**3-337.E0/0.1920996E9*
     #RS(2,1)*RS(2,2)**2-19.E0/889350.E0*RS(2,1)**2*RS(2,3)-RS(2,3)**3/2
     #10-RS(2,3)**2*RS(2,4)/140+13.E0/142296.E0*RS(2,3)**2*RS(2,2)-881.E
     #0/0.160083E8*RS(2,2)**2*RS(2,4)+17.E0/381150.E0*RS(2,1)**2*RS(2,4)
     #+541.E0/0.960498E8*RS(2,1)**2*RS(2,2)+247.E0/0.4002075E7*RS(2,1)*R
     #S(2,4)*RS(2,2)+RS(2,4)**3/42+13.E0/142296.E0*RS(2,1)*RS(2,3)**2+49
     #9.E0/0.137214E7*RS(2,3)*RS(2,4)*RS(2,2)-37517.E0/0.480249E8*RS(2,1
     #)*RS(2,4)**2+17.E0/381150.E0*RS(2,1)*RS(2,3)*RS(2,4)+13.E0/142296.
     #E0*RS(2,1)**3-106.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2))*RS(1,4)*
     #*2*RS(1,3)
      s13 = (-RS(2,3)*RS(2,4)**2/42-RS(2,3)**2*RS(2,4)/105-RS(2,3)**3/42
     #0+RS(2,1)**2*RS(2,4)/105+RS(2,1)*RS(2,4)**2/42+RS(2,1)**3/420)*RS(
     #1,4)**3
      s11 = s12+s13
      s6 = s10+s11
	ARRSS(4) = S6
C      s7 = phi(4)
C      s5 = s6*s7

      s10 = (56311.E0/0.1200622E8*RS(2,1)**2*RS(2,2)+521.E0/533610.E0*RS
     #(2,2)**3+49.E0/490050.E0*RS(2,3)*RS(2,4)**2+137.E0/0.800415E7*RS(2
     #,3)**2*RS(2,2)-49.E0/490050.E0*RS(2,2)**2*RS(2,3)-56311.E0/0.12006
     #22E8*RS(2,1)**2*RS(2,4)+4784.E0/0.1715175E7*RS(2,1)*RS(2,2)**2-137
     #.E0/0.800415E7*RS(2,3)**2*RS(2,4)-509.E0/0.4002075E7*RS(2,2)**2*RS
     #(2,4)+509.E0/0.4002075E7*RS(2,4)**2*RS(2,2)-4784.E0/0.1715175E7*RS
     #(2,1)*RS(2,4)**2-521.E0/533610.E0*RS(2,4)**3-4286.E0/0.1200622E8*R
     #S(2,1)*RS(2,3)*RS(2,2)+4286.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4)
     #)*RS(1,1)**3
      s13 = (19.E0/800415.E0*RS(2,2)**2*RS(2,3)-31.E0/686070.E0*RS(2,1)*
     #RS(2,3)**2-2707.E0/0.800415E7*RS(2,2)**2*RS(2,4)+521.E0/177870.E0*
     #RS(2,1)*RS(2,2)**2-19.E0/800415.E0*RS(2,1)*RS(2,4)**2+137.E0/0.800
     #415E7*RS(2,3)**3+RS(2,3)*RS(2,4)**2/137214+49.E0/490050.E0*RS(2,4)
     #**3-16.E0/0.4002075E7*RS(2,4)**2*RS(2,2)-12482.E0/0.1200622E8*RS(2
     #,1)**2*RS(2,4)+6103.E0/0.1200622E8*RS(2,1)**2*RS(2,3)+17.E0/444675
     #.E0*RS(2,1)*RS(2,3)*RS(2,2)+1312.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS
     #(2,4)-1622.E0/0.4002075E7*RS(2,1)*RS(2,4)*RS(2,2)-56311.E0/0.12006
     #22E8*RS(2,1)**3+191.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)+4784.E0
     #/0.1715175E7*RS(2,2)**3+32.E0/0.4002075E7*RS(2,3)**2*RS(2,2))*RS(1
     #,2)
      s15 = (509.E0/0.4002075E7*RS(2,4)**3-509.E0/0.4002075E7*RS(2,2)**3
     #-2707.E0/0.800415E7*RS(2,1)*RS(2,2)**2+313.E0/0.800415E7*RS(2,2)**
     #2*RS(2,4)+31.E0/686070.E0*RS(2,3)**2*RS(2,4)-31.E0/686070.E0*RS(2,
     #3)**2*RS(2,2)-6103.E0/0.1200622E8*RS(2,1)**2*RS(2,2)-313.E0/0.8004
     #15E7*RS(2,4)**2*RS(2,2)+2707.E0/0.800415E7*RS(2,1)*RS(2,4)**2+1159
     #.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-16.E0/0.4002075E7*RS(2,3)*
     #RS(2,4)**2-1159.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4)+6103.E0/0.1
     #200622E8*RS(2,1)**2*RS(2,4)+16.E0/0.4002075E7*RS(2,2)**2*RS(2,3))*
     #RS(1,3)
      s16 = (-49.E0/490050.E0*RS(2,2)**3-19.E0/800415.E0*RS(2,3)*RS(2,4)
     #**2-1312.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-4784.E0/0.1715175E
     #7*RS(2,4)**3+56311.E0/0.1200622E8*RS(2,1)**3-137.E0/0.800415E7*RS(
     #2,3)**3-RS(2,2)**2*RS(2,3)/137214+16.E0/0.4002075E7*RS(2,2)**2*RS(
     #2,4)+2707.E0/0.800415E7*RS(2,4)**2*RS(2,2)-191.E0/0.1200622E8*RS(2
     #,3)*RS(2,4)*RS(2,2)+19.E0/800415.E0*RS(2,1)*RS(2,2)**2-521.E0/1778
     #70.E0*RS(2,1)*RS(2,4)**2-32.E0/0.4002075E7*RS(2,3)**2*RS(2,4)-17.E
     #0/444675.E0*RS(2,1)*RS(2,3)*RS(2,4)+1622.E0/0.4002075E7*RS(2,1)*RS
     #(2,4)*RS(2,2)-6103.E0/0.1200622E8*RS(2,1)**2*RS(2,3)+12482.E0/0.12
     #00622E8*RS(2,1)**2*RS(2,2)+31.E0/686070.E0*RS(2,1)*RS(2,3)**2)*RS(
     #1,4)
      s14 = s15+s16
      s12 = s13+s14
      s13 = RS(1,1)**2
      s11 = s12*s13
      s9 = s10+s11
      s11 = s9
      s15 = (2707.E0/0.800415E7*RS(2,1)**2*RS(2,3)-32.E0/0.4002075E7*RS(
     #2,1)*RS(2,4)**2+56311.E0/0.1200622E8*RS(2,2)**3-RS(2,3)**2*RS(2,4)
     #/137214-1312.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)-6103.E0/0.1200
     #622E8*RS(2,2)**2*RS(2,4)+12482.E0/0.1200622E8*RS(2,2)**2*RS(2,3)-1
     #91.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4)-17.E0/444675.E0*RS(2,1)*
     #RS(2,4)*RS(2,2)-49.E0/490050.E0*RS(2,3)**3-521.E0/177870.E0*RS(2,1
     #)**2*RS(2,2)-137.E0/0.800415E7*RS(2,4)**3-19.E0/800415.E0*RS(2,1)*
     #*2*RS(2,4)+1622.E0/0.4002075E7*RS(2,1)*RS(2,3)*RS(2,2)-4784.E0/0.1
     #715175E7*RS(2,1)**3+16.E0/0.4002075E7*RS(2,1)*RS(2,3)**2+31.E0/686
     #070.E0*RS(2,4)**2*RS(2,2)+19.E0/800415.E0*RS(2,3)**2*RS(2,2))*RS(1
     #,2)**2
      s18 = (-191.E0/0.1200622E8*RS(2,1)*RS(2,4)**2+1312.E0/0.1200622E8*
     #RS(2,1)**2*RS(2,4)-1622.E0/0.4002075E7*RS(2,1)*RS(2,2)**2+17.E0/44
     #4675.E0*RS(2,3)**2*RS(2,2)-4286.E0/0.1200622E8*RS(2,3)**3-17.E0/44
     #4675.E0*RS(2,1)**2*RS(2,2)-1312.E0/0.1200622E8*RS(2,3)**2*RS(2,4)-
     #1159.E0/0.1200622E8*RS(2,1)**2*RS(2,3)+1159.E0/0.1200622E8*RS(2,1)
     #*RS(2,3)**2+4286.E0/0.1200622E8*RS(2,1)**3+191.E0/0.1200622E8*RS(2
     #,3)*RS(2,4)**2+1114.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2)-1114.E0
     #/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)+1622.E0/0.4002075E7*RS(2,2)**
     #2*RS(2,3))*RS(1,3)
      s19 = (-4286.E0/0.1200622E8*RS(2,2)**3+4286.E0/0.1200622E8*RS(2,4)
     #**3+191.E0/0.1200622E8*RS(2,3)**2*RS(2,2)-1312.E0/0.1200622E8*RS(2
     #,2)**2*RS(2,3)+1159.E0/0.1200622E8*RS(2,2)**2*RS(2,4)-1159.E0/0.12
     #00622E8*RS(2,4)**2*RS(2,2)+1312.E0/0.1200622E8*RS(2,3)*RS(2,4)**2+
     #17.E0/444675.E0*RS(2,1)*RS(2,2)**2-17.E0/444675.E0*RS(2,1)*RS(2,4)
     #**2+1114.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,4)-1114.E0/0.1200622E
     #8*RS(2,1)*RS(2,3)*RS(2,2)-1622.E0/0.4002075E7*RS(2,1)**2*RS(2,4)+1
     #622.E0/0.4002075E7*RS(2,1)**2*RS(2,2)-191.E0/0.1200622E8*RS(2,3)**
     #2*RS(2,4))*RS(1,4)
      s17 = s18+s19
      s18 = RS(1,2)
      s16 = s17*s18
      s14 = s15+s16
      s15 = s14+(31.E0/686070.E0*RS(2,1)**2*RS(2,2)-313.E0/0.800415E7*RS
     #(2,2)**2*RS(2,4)+313.E0/0.800415E7*RS(2,4)**2*RS(2,2)-2707.E0/0.80
     #0415E7*RS(2,3)*RS(2,4)**2+6103.E0/0.1200622E8*RS(2,3)**2*RS(2,2)+5
     #09.E0/0.4002075E7*RS(2,2)**3-16.E0/0.4002075E7*RS(2,1)*RS(2,2)**2+
     #16.E0/0.4002075E7*RS(2,1)*RS(2,4)**2-509.E0/0.4002075E7*RS(2,4)**3
     #+2707.E0/0.800415E7*RS(2,2)**2*RS(2,3)-1159.E0/0.1200622E8*RS(2,1)
     #*RS(2,3)*RS(2,2)-31.E0/686070.E0*RS(2,1)**2*RS(2,4)-6103.E0/0.1200
     #622E8*RS(2,3)**2*RS(2,4)+1159.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,
     #4))*RS(1,3)**2
      s16 = s15
      s18 = (191.E0/0.1200622E8*RS(2,1)*RS(2,2)**2-1312.E0/0.1200622E8*R
     #S(2,1)**2*RS(2,2)+17.E0/444675.E0*RS(2,1)**2*RS(2,4)-1159.E0/0.120
     #0622E8*RS(2,1)*RS(2,3)**2-191.E0/0.1200622E8*RS(2,2)**2*RS(2,3)-17
     #.E0/444675.E0*RS(2,3)**2*RS(2,4)+4286.E0/0.1200622E8*RS(2,3)**3-11
     #14.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2)-4286.E0/0.1200622E8*RS(2
     #,1)**3-1622.E0/0.4002075E7*RS(2,3)*RS(2,4)**2+1312.E0/0.1200622E8*
     #RS(2,3)**2*RS(2,2)+1622.E0/0.4002075E7*RS(2,1)*RS(2,4)**2+1114.E0/
     #0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)+1159.E0/0.1200622E8*RS(2,1)**2
     #*RS(2,3))*RS(1,4)*RS(1,3)
      s19 = (49.E0/490050.E0*RS(2,3)**3+137.E0/0.800415E7*RS(2,2)**3-563
     #11.E0/0.1200622E8*RS(2,4)**3+4784.E0/0.1715175E7*RS(2,1)**3-31.E0/
     #686070.E0*RS(2,2)**2*RS(2,4)+RS(2,3)**2*RS(2,2)/137214-19.E0/80041
     #5.E0*RS(2,3)**2*RS(2,4)-12482.E0/0.1200622E8*RS(2,3)*RS(2,4)**2+52
     #1.E0/177870.E0*RS(2,1)**2*RS(2,4)+19.E0/800415.E0*RS(2,1)**2*RS(2,
     #2)+6103.E0/0.1200622E8*RS(2,4)**2*RS(2,2)+32.E0/0.4002075E7*RS(2,1
     #)*RS(2,2)**2-16.E0/0.4002075E7*RS(2,1)*RS(2,3)**2+17.E0/444675.E0*
     #RS(2,1)*RS(2,4)*RS(2,2)+191.E0/0.1200622E8*RS(2,1)*RS(2,3)*RS(2,2)
     #-2707.E0/0.800415E7*RS(2,1)**2*RS(2,3)+1312.E0/0.1200622E8*RS(2,3)
     #*RS(2,4)*RS(2,2)-1622.E0/0.4002075E7*RS(2,1)*RS(2,3)*RS(2,4))*RS(1
     #,4)**2
      s17 = s18+s19
      s13 = s16+s17
      s14 = RS(1,1)
      s12 = s13*s14
      s10 = s11+s12
      s11 = s10
      s13 = (-521.E0/533610.E0*RS(2,1)**3-509.E0/0.4002075E7*RS(2,1)*RS(
     #2,3)**2+509.E0/0.4002075E7*RS(2,1)**2*RS(2,3)+56311.E0/0.1200622E8
     #*RS(2,2)**2*RS(2,3)-49.E0/490050.E0*RS(2,3)**2*RS(2,4)-56311.E0/0.
     #1200622E8*RS(2,1)*RS(2,2)**2+137.E0/0.800415E7*RS(2,3)*RS(2,4)**2+
     #4784.E0/0.1715175E7*RS(2,3)**2*RS(2,2)-4286.E0/0.1200622E8*RS(2,3)
     #*RS(2,4)*RS(2,2)+521.E0/533610.E0*RS(2,3)**3+4286.E0/0.1200622E8*R
     #S(2,1)*RS(2,4)*RS(2,2)+49.E0/490050.E0*RS(2,1)**2*RS(2,4)-4784.E0/
     #0.1715175E7*RS(2,1)**2*RS(2,2)-137.E0/0.800415E7*RS(2,1)*RS(2,4)**
     #2)*RS(1,2)**3
      s16 = (32.E0/0.4002075E7*RS(2,3)*RS(2,4)**2+4784.E0/0.1715175E7*RS
     #(2,3)**3-19.E0/800415.E0*RS(2,1)**2*RS(2,2)-56311.E0/0.1200622E8*R
     #S(2,2)**3-31.E0/686070.E0*RS(2,4)**2*RS(2,2)+521.E0/177870.E0*RS(2
     #,3)**2*RS(2,2)+6103.E0/0.1200622E8*RS(2,2)**2*RS(2,4)+137.E0/0.800
     #415E7*RS(2,4)**3+RS(2,1)**2*RS(2,4)/137214-16.E0/0.4002075E7*RS(2,
     #1)**2*RS(2,3)-2707.E0/0.800415E7*RS(2,1)*RS(2,3)**2+49.E0/490050.E
     #0*RS(2,1)**3-12482.E0/0.1200622E8*RS(2,1)*RS(2,2)**2-1622.E0/0.400
     #2075E7*RS(2,1)*RS(2,3)*RS(2,2)+1312.E0/0.1200622E8*RS(2,1)*RS(2,4)
     #*RS(2,2)+17.E0/444675.E0*RS(2,3)*RS(2,4)*RS(2,2)+191.E0/0.1200622E
     #8*RS(2,1)*RS(2,3)*RS(2,4)+19.E0/800415.E0*RS(2,3)**2*RS(2,4))*RS(1
     #,3)
      s17 = (1159.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)-6103.E0/0.12006
     #22E8*RS(2,2)**2*RS(2,3)+509.E0/0.4002075E7*RS(2,1)**3+2707.E0/0.80
     #0415E7*RS(2,1)**2*RS(2,2)-509.E0/0.4002075E7*RS(2,3)**3+31.E0/6860
     #70.E0*RS(2,1)*RS(2,4)**2-2707.E0/0.800415E7*RS(2,3)**2*RS(2,2)-31.
     #E0/686070.E0*RS(2,3)*RS(2,4)**2+6103.E0/0.1200622E8*RS(2,1)*RS(2,2
     #)**2+313.E0/0.800415E7*RS(2,1)*RS(2,3)**2-1159.E0/0.1200622E8*RS(2
     #,1)*RS(2,4)*RS(2,2)-313.E0/0.800415E7*RS(2,1)**2*RS(2,3)+16.E0/0.4
     #002075E7*RS(2,3)**2*RS(2,4)-16.E0/0.4002075E7*RS(2,1)**2*RS(2,4))*
     #RS(1,4)
      s15 = s16+s17
      s16 = RS(1,2)**2
      s14 = s15*s16
      s12 = s13+s14
      s8 = s11+s12
      s10 = s8
      s14 = (2707.E0/0.800415E7*RS(2,2)**2*RS(2,4)-521.E0/177870.E0*RS(2
     #,2)**2*RS(2,3)+12482.E0/0.1200622E8*RS(2,3)**2*RS(2,4)+31.E0/68607
     #0.E0*RS(2,1)**2*RS(2,3)+16.E0/0.4002075E7*RS(2,4)**2*RS(2,2)-32.E0
     #/0.4002075E7*RS(2,1)**2*RS(2,2)+56311.E0/0.1200622E8*RS(2,3)**3-19
     #.E0/800415.E0*RS(2,1)*RS(2,2)**2-RS(2,1)*RS(2,4)**2/137214-17.E0/4
     #44675.E0*RS(2,1)*RS(2,3)*RS(2,2)-6103.E0/0.1200622E8*RS(2,1)*RS(2,
     #3)**2-4784.E0/0.1715175E7*RS(2,2)**3+1622.E0/0.4002075E7*RS(2,3)*R
     #S(2,4)*RS(2,2)-191.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(2,2)-137.E0/0
     #.800415E7*RS(2,1)**3-49.E0/490050.E0*RS(2,4)**3-1312.E0/0.1200622E
     #8*RS(2,1)*RS(2,3)*RS(2,4)+19.E0/800415.E0*RS(2,3)*RS(2,4)**2)*RS(1
     #,3)**2
      s16 = (-1159.E0/0.1200622E8*RS(2,2)**2*RS(2,4)-1114.E0/0.1200622E8
     #*RS(2,1)*RS(2,3)*RS(2,4)+4286.E0/0.1200622E8*RS(2,2)**3+1114.E0/0.
     #1200622E8*RS(2,1)*RS(2,3)*RS(2,2)+1312.E0/0.1200622E8*RS(2,1)*RS(2
     #,2)**2+1159.E0/0.1200622E8*RS(2,4)**2*RS(2,2)-191.E0/0.1200622E8*R
     #S(2,1)**2*RS(2,2)+17.E0/444675.E0*RS(2,3)*RS(2,4)**2+1622.E0/0.400
     #2075E7*RS(2,3)**2*RS(2,4)-1622.E0/0.4002075E7*RS(2,3)**2*RS(2,2)-1
     #7.E0/444675.E0*RS(2,2)**2*RS(2,3)+191.E0/0.1200622E8*RS(2,1)**2*RS
     #(2,4)-1312.E0/0.1200622E8*RS(2,1)*RS(2,4)**2-4286.E0/0.1200622E8*R
     #S(2,4)**3)*RS(1,4)*RS(1,3)
      s17 = (-1159.E0/0.1200622E8*RS(2,3)*RS(2,4)*RS(2,2)+509.E0/0.40020
     #75E7*RS(2,3)**3-313.E0/0.800415E7*RS(2,1)*RS(2,3)**2-509.E0/0.4002
     #075E7*RS(2,1)**3+2707.E0/0.800415E7*RS(2,3)**2*RS(2,4)+6103.E0/0.1
     #200622E8*RS(2,3)*RS(2,4)**2-2707.E0/0.800415E7*RS(2,1)**2*RS(2,4)+
     #31.E0/686070.E0*RS(2,2)**2*RS(2,3)+1159.E0/0.1200622E8*RS(2,1)*RS(
     #2,4)*RS(2,2)+313.E0/0.800415E7*RS(2,1)**2*RS(2,3)-31.E0/686070.E0*
     #RS(2,1)*RS(2,2)**2-16.E0/0.4002075E7*RS(2,3)**2*RS(2,2)+16.E0/0.40
     #02075E7*RS(2,1)**2*RS(2,2)-6103.E0/0.1200622E8*RS(2,1)*RS(2,4)**2)
     #*RS(1,4)**2
      s15 = s16+s17
      s13 = s14+s15
      s14 = RS(1,2)
      s12 = s13*s14
      s13 = (-137.E0/0.800415E7*RS(2,1)**2*RS(2,2)+509.E0/0.4002075E7*RS
     #(2,2)**2*RS(2,4)-4784.E0/0.1715175E7*RS(2,2)**2*RS(2,3)+4286.E0/0.
     #1200622E8*RS(2,1)*RS(2,3)*RS(2,2)-521.E0/533610.E0*RS(2,2)**3-49.E
     #0/490050.E0*RS(2,1)*RS(2,4)**2+521.E0/533610.E0*RS(2,4)**3+4784.E0
     #/0.1715175E7*RS(2,3)*RS(2,4)**2-56311.E0/0.1200622E8*RS(2,3)**2*RS
     #(2,2)+56311.E0/0.1200622E8*RS(2,3)**2*RS(2,4)-509.E0/0.4002075E7*R
     #S(2,4)**2*RS(2,2)+137.E0/0.800415E7*RS(2,1)**2*RS(2,4)-4286.E0/0.1
     #200622E8*RS(2,1)*RS(2,3)*RS(2,4)+49.E0/490050.E0*RS(2,1)*RS(2,2)**
     #2)*RS(1,3)**3
      s11 = s12+s13
      s9 = s10+s11
      s10 = s9+(-56311.E0/0.1200622E8*RS(2,3)**3-16.E0/0.4002075E7*RS(2,
     #2)**2*RS(2,4)+RS(2,1)*RS(2,2)**2/137214+19.E0/800415.E0*RS(2,1)*RS
     #(2,4)**2+137.E0/0.800415E7*RS(2,1)**3+1312.E0/0.1200622E8*RS(2,1)*
     #RS(2,3)*RS(2,2)+49.E0/490050.E0*RS(2,2)**3+521.E0/177870.E0*RS(2,3
     #)*RS(2,4)**2+17.E0/444675.E0*RS(2,1)*RS(2,3)*RS(2,4)-19.E0/800415.
     #E0*RS(2,2)**2*RS(2,3)-1622.E0/0.4002075E7*RS(2,3)*RS(2,4)*RS(2,2)+
     #6103.E0/0.1200622E8*RS(2,1)*RS(2,3)**2+4784.E0/0.1715175E7*RS(2,4)
     #**3+32.E0/0.4002075E7*RS(2,1)**2*RS(2,4)+191.E0/0.1200622E8*RS(2,1
     #)*RS(2,4)*RS(2,2)-2707.E0/0.800415E7*RS(2,4)**2*RS(2,2)-12482.E0/0
     #.1200622E8*RS(2,3)**2*RS(2,2)-31.E0/686070.E0*RS(2,1)**2*RS(2,3))*
     #RS(1,4)*RS(1,3)**2
      s11 = s10
      s13 = (-4784.E0/0.1715175E7*RS(2,3)**3-191.E0/0.1200622E8*RS(2,1)*
     #RS(2,3)*RS(2,2)-17.E0/444675.E0*RS(2,3)*RS(2,4)*RS(2,2)-49.E0/4900
     #50.E0*RS(2,1)**3-32.E0/0.4002075E7*RS(2,2)**2*RS(2,3)-137.E0/0.800
     #415E7*RS(2,2)**3+56311.E0/0.1200622E8*RS(2,4)**3-1312.E0/0.1200622
     #E8*RS(2,1)*RS(2,4)*RS(2,2)+16.E0/0.4002075E7*RS(2,1)**2*RS(2,3)-61
     #03.E0/0.1200622E8*RS(2,4)**2*RS(2,2)+31.E0/686070.E0*RS(2,2)**2*RS
     #(2,4)-19.E0/800415.E0*RS(2,3)**2*RS(2,2)-521.E0/177870.E0*RS(2,3)*
     #*2*RS(2,4)+1622.E0/0.4002075E7*RS(2,1)*RS(2,3)*RS(2,4)+2707.E0/0.8
     #00415E7*RS(2,1)*RS(2,3)**2+19.E0/800415.E0*RS(2,1)**2*RS(2,4)-RS(2
     #,1)**2*RS(2,2)/137214+12482.E0/0.1200622E8*RS(2,1)*RS(2,4)**2)*RS(
     #1,4)**2*RS(1,3)
      s14 = (49.E0/490050.E0*RS(2,3)**2*RS(2,2)-49.E0/490050.E0*RS(2,1)*
     #*2*RS(2,2)+509.E0/0.4002075E7*RS(2,1)*RS(2,3)**2+4286.E0/0.1200622
     #E8*RS(2,3)*RS(2,4)*RS(2,2)-4286.E0/0.1200622E8*RS(2,1)*RS(2,4)*RS(
     #2,2)-4784.E0/0.1715175E7*RS(2,3)**2*RS(2,4)+137.E0/0.800415E7*RS(2
     #,1)*RS(2,2)**2-509.E0/0.4002075E7*RS(2,1)**2*RS(2,3)-521.E0/533610
     #.E0*RS(2,3)**3-137.E0/0.800415E7*RS(2,2)**2*RS(2,3)-56311.E0/0.120
     #0622E8*RS(2,3)*RS(2,4)**2+521.E0/533610.E0*RS(2,1)**3+4784.E0/0.17
     #15175E7*RS(2,1)**2*RS(2,4)+56311.E0/0.1200622E8*RS(2,1)*RS(2,4)**2
     #)*RS(1,4)**3
      s12 = s13+s14
      s7 = s11+s12
	ARRSS(5) = S7
C      s8 = phi(5)
C      s6 = s7*s8
C      s4 = s5+s6
C      t0 = s3+s4


      RETURN
      END
C
C=====================================================================
      SUBROUTINE MOONEY (DISD,DR,PROP,STRESS)                          
      IMPLICIT REAL*8(A-H,O-Z)                                                  
C     ----------------------------------------------------------------
C                                                                               
C     INCOMPRESSIBLE NONLINEAR ELASTIC MATERIAL                                 
C         ( MOONEY-RIVLIN DESCRIPTION IN STATE OF PLANE STRESS )                
C                                                                               
C     ----------------------------------------------------------------
      DIMENSION DISD(5),DR(6,6),PROP(5),STRESS(4)
      DIMENSION STRAIN(4)
	DIMENSION C(3,3)                              

C	-------------------------
C	DEFINE MATERIAL CONSTANTS
C	-------------------------
      C1 = PROP(3)                                                                
      C2 = PROP(4)                                                                

C	------------------------------
C	CALCULATE DEFORMATION GRADIENT
C	------------------------------
	X11 = 1.0+DISD(1)
	X22 = 1.0+DISD(2)
	X12 = DISD(3)
	X21 = DISD(4)

C	-----------------------------------------------
C	CALCULATE RIGHT CAUCHY-GREEN DEFORMATION TENSOR
C	-----------------------------------------------
	C11 = X11**2+X21**2
	C22 = X22**2+X12**2
	C12 = X11*X12+X22*X21
	C21 = C12

C	-------------------------
C	CALCULATE EXTENSION RATIO
C	-------------------------
      DENOM = C11*C22 - C12*C12                                                   
      EX2   = 1.0/DENOM
      EX4   = EX2*EX2                                                               
      C1122 = C11 + C22                                                           

C	--------------------------------------
C	CALCULATE 2ND PIOLA-KIRCHHOFF STRESSES
C	--------------------------------------
      T1 = 2.0*(C1 + C2*EX2)
      T2 = 2.0*(C1*EX4 - C2*(1.D0 - EX4*C1122))
      S11       = T1 - T2*C22                                                     
      S22       = T1 - T2*C11                                                     
      S12       = T2*C12                                                          

C	----------------------------------------------
C	CALCULATE DENSITY RATIO
C       DET    = (INITIAL DENSITY / CURRENT DENSITY)
C       INCOMPRESSIBILITY CONDITION YIELDS  DET=1.
C	----------------------------------------------
      DET=1.0

C	-------------------------
C	CALCULATE CAUCHY STRESSES
C	-------------------------
      STRESS(1) = DET*(S11*X11*X11+2.0*S12*X11*X12+S22*X12*X12)
      STRESS(2) = DET*(S11*X21*X21+2.0*S12*X21*X22+S22*X22*X22)
      STRESS(3) = DET*(S11*X11*X21+S12*(X11*X22+X12*X21)+S22*X12*X22)            
      STRESS(4) = 0.0

C	------------------------------
C	CALCULATE STRESS-STRAIN MATRIX
C	------------------------------
	C1F  = C1*EX4*4.0
	C2F  = C2*EX4*4.0
	EC11 = EX2*C11
	EC22 = EX2*C22
	EC1  = EC11*C1122
	EC2  = EC22*C1122
	EC12 = EC11*C22
C
	C(1,1) = 2.0*C22*(C1F*EC22 + C2F*(-1.0 + EC2))
	C(1,2) = C1F*(-1.0 + 2.0*EC12) + C2F*(1.0/EX4 + 2.0*C1122*(-1.0
     *         + EC12))
	C(1,3) = C12*(-2.0*C1F*EC22 + C2F*(1.0 - 2.0*EC2))
	C(2,2) = 2.0*C11*(C1F*EC11 + C2F*(-1.0 + EC1))
	C(2,3) = C12*(-2.0*C1F*EC11 + C2F*(1.0 - 2.0*EC1))
	C(3,3) = (2.0*C12*C12*EX2 + 0.50)*(C1F + C2F*C1122) - 2.0*C2

	C(2,1) = C(1,2)
	C(3,1) = C(1,3)
	C(3,2 )= C(2,3)
C	-------------------------
C	TRANSFORM TO CURRENT TIME
C	-------------------------
	DR(1,1) = X11*X11*C(1,1)*X11*X11 + X11*X11*C(1,2)*X12*X12+
	1                                           X11*X11*C(1,3)*X11*X12+
	2          X12*X12*C(2,1)*X11*X11 + X12*X12*C(2,2)*X12*X12+
	3	                                       X12*X12*C(2,3)*X11*X12+
	4          X11*X12*C(3,1)*X11*X11 + X11*X12*C(3,2)*X12*X12+
	5                                           X11*X12*C(3,3)*X11*X12

	DR(1,2) = X11*X11*C(1,1)*X21*X21 + X11*X11*C(1,2)*X22*X22+
	1                                           X11*X11*C(1,3)*X21*X22+
	2          X12*X12*C(2,1)*X21*X21 + X12*X12*C(2,2)*X22*X22+
	3	                                       X12*X12*C(2,3)*X21*X22+
	4          X11*X12*C(3,1)*X22*X22 + X11*X12*C(3,2)*X22*X22+
	5                                           X11*X12*C(3,3)*X21*X22

	DR(1,3) = X11*X11*C(1,1)*X11*X21 + X11*X11*C(1,2)*X12*X22+
	1                                           X11*X11*C(1,3)*X11*X22+
	2          X12*X12*C(2,1)*X11*X21 + X12*X12*C(2,2)*X12*X22+
	3	                                       X12*X12*C(2,3)*X11*X22+
	4          X11*X12*C(3,1)*X11*X21 + X11*X12*C(3,2)*X12*X22+
	5                                           X11*X12*C(3,3)*X11*X22

	DR(2,1) = DR(1,2)

	DR(2,2) = X21*X21*C(1,1)*X21*X21 + X21*X21*C(1,2)*X22*X22+
	1                                           X21*X21*C(1,3)*X21*X22+
	2          X22*X22*C(2,1)*X21*X21 + X22*X22*C(2,2)*X22*X22+
	3	                                       X22*X22*C(2,3)*X22*X21+
	4          X21*X21*C(3,1)*X21*X22 + X12*X22*C(3,2)*X22*X22+
	5                                           X21*X21*C(3,3)*X22*X22

	DR(2,3) = X21*X21*C(1,1)*X21*X11 + X21*X21*C(1,2)*X12*X22+
	1                                           X21*X21*C(1,3)*X11*X22+
	2          X22*X22*C(2,1)*X11*X21 + X22*X22*C(2,2)*X22*X12+
	3	                                       X22*X22*C(2,3)*X22*X11+
	4          X21*X21*C(3,1)*X22*X11 + X21*X22*C(3,2)*X22*X12+
	5                                           X21*X22*C(3,3)*X22*X11

	DR(3,1) = DR(1,3)

	DR(3,2) = DR(2,3)

	DR(3,3) = X11*X11*C(1,1)*X12*X12 + X11*X21*C(1,2)*X12*X22+
	1                                           X11*X11*C(1,3)*X21*X22+
	2          X12*X22*C(2,1)*X11*X21 + X12*X12*C(2,2)*X22*X22+
	3	                                       X22*X22*C(2,3)*X12*X11+
	4          X11*X11*C(3,1)*X22*X21 + X11*X22*C(3,2)*X22*X12+
	5                                           X11*X11*C(3,3)*X22*X22

      DR(1,4) = 0.E0
      DR(1,5) = 0.E0
      DR(1,6) = 0.E0
      DR(2,4) = 0.E0
      DR(2,5) = 0.E0
      DR(2,6) = 0.E0
      DR(3,4) = 0.E0
      DR(3,5) = 0.E0
      DR(3,6) = 0.E0
      DR(4,1) = 0.E0
      DR(4,2) = 0.E0
      DR(4,3) = 0.E0
      DR(4,4) = 0.E0
      DR(4,5) = 0.E0
      DR(4,6) = 0.E0
      DR(5,1) = 0.E0
      DR(5,2) = 0.E0
      DR(5,3) = 0.E0
      DR(5,4) = 0.E0
      DR(5,5) = 0.E0
      DR(5,6) = 0.E0
      DR(6,1) = 0.E0
      DR(6,2) = 0.E0
      DR(6,3) = 0.E0
      DR(6,4) = 0.E0
      DR(6,5) = 0.E0
      DR(6,6) = 0.E0


      RETURN                                                                    
      END                                                                       
C
C=====================================================================
      SUBROUTINE CMEM(RS,RN,SN,SLN,CM)
	IMPLICIT REAL*8 (A-H,O-Z)
C	-------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Cm FOR MEMBRANE STRAINS
C				USING THE ALLMAN STRING FUNCTIONS
C
C	INPUT VARIABLES
C	RS(2,4)		= ELEMENT COORDINATES IN LOCAL SPACE
C	RN(4),SN(4)	= ELEMENET BOUNDARY NORMALS AND TANGENTIALS
C	SLN(4)		= BOUNDARY LENGTHS
C
C	OUTPUT VARIAVBLES
C	CM(9,12)	= MEMBRANE AND DRILLING STRAIN PARAMETER
C	-------------------------------------------------------
      DIMENSION RS(2,4),RN(4),SN(4),SLN(4)
	DIMENSION CM(9,12)

C	------------
C	SOLVE FOR Cm
C	------------
	CM(1,1) =
     #((RN(1)*SN(1)**2/2+RN(1)**3/2)*SLN(1)+(RN(4)**3/2+RN(4)*SN(4)
     #**2/2)*SLN(4))
	CM(1,4) =((RN(1)*SN(1)**2/2+RN(1)**3/2)*SLN(1)+(RN(2)*
     #*3/2+RN(2)*SN(2)**2/2)*SLN(2))
	CM(1,7)=((RN(2)**3/2+RN(2)*SN(2)**2/2
     #)*SLN(2)+(RN(3)*SN(3)**2/2+RN(3)**3/2)*SLN(3))
	CM(1,10)=((RN(3)*SN(3)
     #**2/2+RN(3)**3/2)*SLN(3)+(RN(4)**3/2+RN(4)*SN(4)**2/2)*SLN(4))
	CM(1,3) =(-RN(1)**2*SLN(1)**2/12+RN(4)**2*SLN(4)**2/12)
	CM(1,6)=(-RN(2)**2*SLN(2)**2/12+RN(1)**2*SLN(1)**2/12)
	CM(1,9)=(RN(2)**2*SLN(2)**2/12-RN(3)**2*SLN(3)**2/12)
	CM(1,12)=(-RN(4)**2*SLN(4)**2/12+RN(3)**2*SLN(3)**2/12)

	CM(2,1) =
     #(((RS(2,1)/3+RS(2,2)/6)*RN(1)**3+(SN(1)**2*RS(2,2)/6+SN(1)**2
     #*RS(2,1)/3)*RN(1))*SLN(1)+((RS(2,4)/6+RS(2,1)/3)*RN(4)**3+(SN(4)**
     #2*RS(2,4)/6+SN(4)**2*RS(2,1)/3)*RN(4))*SLN(4))
	CM(2,4) =(((RS(2,2)/3+
     #RS(2,1)/6)*RN(1)**3+(SN(1)**2*RS(2,1)/6+SN(1)**2*RS(2,2)/3)*RN(1))
     #*SLN(1)+((RS(2,2)/3+RS(2,3)/6)*RN(2)**3+(SN(2)**2*RS(2,3)/6+SN(2)*
     #*2*RS(2,2)/3)*RN(2))*SLN(2))
	CM(2,7)=
     #(((RS(2,3)/3+RS(2,2)/6)*RN(2)**3+(SN(2)**2*RS(2,2)/6+SN(2)
     #**2*RS(2,3)/3)*RN(2))*SLN(2)+((RS(2,3)/3+RS(2,4)/6)*RN(3)**3+(SN(3
     #)**2*RS(2,4)/6+SN(3)**2*RS(2,3)/3)*RN(3))*SLN(3))
	CM(2,10)=(((RS(2,3)
     #/6+RS(2,4)/3)*RN(3)**3+(SN(3)**2*RS(2,3)/6+SN(3)**2*RS(2,4)/3)*RN(
     #3))*SLN(3)+((RS(2,4)/3+RS(2,1)/6)*RN(4)**3+(SN(4)**2*RS(2,1)/6+SN(
     #4)**2*RS(2,4)/3)*RN(4))*SLN(4))
	CM(2,3) =
     #((-RS(2,2)/24-RS(2,1)/24)*RN(1)**2*SLN(1)**2+(RS(2,1)/24+R
     #S(2,4)/24)*RN(4)**2*SLN(4)**2)
	CM(2,6)=((RS(2,1)/24+RS(2,2)/24)*RN(
     #1)**2*SLN(1)**2+(-RS(2,3)/24-RS(2,2)/24)*RN(2)**2*SLN(2)**2)
	CM(2,9)=
     #((RS(2,2)/24+RS(2,3)/24)*RN(2)**2*SLN(2)**2+(-RS(2,3)/24-RS(2,4)
     #/24)*RN(3)**2*SLN(3)**2)
	CM(2,12)=((RS(2,3)/24+RS(2,4)/24)*RN(3)**2*
     #SLN(3)**2+(-RS(2,4)/24-RS(2,1)/24)*RN(4)**2*SLN(4)**2)

	CM(3,1) =
     #(((RS(1,1)/3+RS(1,2)/6)*RN(1)**3+(-RS(2,2)/6-RS(2,1)/3)*SN(1)
     #*RN(1)**2+(RS(1,1)/3+RS(1,2)/6)*SN(1)**2*RN(1)+(-RS(2,2)/6-RS(2,1)
     #/3)*SN(1)**3)*SLN(1)+((RS(1,1)/3+RS(1,4)/6)*RN(4)**3+(-RS(2,4)/6-R
     #S(2,1)/3)*SN(4)*RN(4)**2+(RS(1,1)/3+RS(1,4)/6)*SN(4)**2*RN(4)+(-RS
     #(2,4)/6-RS(2,1)/3)*SN(4)**3)*SLN(4))
	CM(3,4) =
     #(((RS(1,1)/6+RS(1,2)/3)*RN(1)**3+(-RS(2,2)/3-RS(2,1)/6)*SN(1)
     #*RN(1)**2+(RS(1,1)/6+RS(1,2)/3)*SN(1)**2*RN(1)+(-RS(2,2)/3-RS(2,1)
     #/6)*SN(1)**3)*SLN(1)+((RS(1,3)/6+RS(1,2)/3)*RN(2)**3+(-RS(2,3)/6-R
     #S(2,2)/3)*SN(2)*RN(2)**2+(RS(1,3)/6+RS(1,2)/3)*SN(2)**2*RN(2)+(-RS
     #(2,3)/6-RS(2,2)/3)*SN(2)**3)*SLN(2))
	CM(3,7)=
     #(((RS(1,3)/3+RS(1,2)/6)*RN(2)**3+(-RS(2,3)/3-RS(2,2)/6)*SN(2)
     #*RN(2)**2+(RS(1,3)/3+RS(1,2)/6)*SN(2)**2*RN(2)+(-RS(2,3)/3-RS(2,2)
     #/6)*SN(2)**3)*SLN(2)+((RS(1,4)/6+RS(1,3)/3)*RN(3)**3+(-RS(2,3)/3-R
     #S(2,4)/6)*SN(3)*RN(3)**2+(RS(1,4)/6+RS(1,3)/3)*SN(3)**2*RN(3)+(-RS
     #(2,3)/3-RS(2,4)/6)*SN(3)**3)*SLN(3))
	CM(3,10)=
     #(((RS(1,3)/6+RS(1,4)/3)*RN(3)**3+(-RS(2,4)/3-RS(2,3)/6)*SN
     #(3)*RN(3)**2+(RS(1,3)/6+RS(1,4)/3)*SN(3)**2*RN(3)+(-RS(2,4)/3-RS(2
     #,3)/6)*SN(3)**3)*SLN(3)+((RS(1,4)/3+RS(1,1)/6)*RN(4)**3+(-RS(2,1)/
     #6-RS(2,4)/3)*SN(4)*RN(4)**2+(RS(1,4)/3+RS(1,1)/6)*SN(4)**2*RN(4)+(
     #-RS(2,1)/6-RS(2,4)/3)*SN(4)**3)*SLN(4))
	CM(3,2) =
     #(((-RS(2,2)/6-RS(2,1)/3)*RN(1)**3+(-RS(2,2)/6-RS(2,1)/3)*S
     #N(1)**2*RN(1))*SLN(1)+((-RS(2,4)/6-RS(2,1)/3)*RN(4)**3+(-RS(2,4)/6
     #-RS(2,1)/3)*SN(4)**2*RN(4))*SLN(4))
	CM(3,5) =(((-RS(2,2)/3-RS(2,1)/6)
     #*RN(1)**3+(-RS(2,2)/3-RS(2,1)/6)*SN(1)**2*RN(1))*SLN(1)+((-RS(2,3)
     #/6-RS(2,2)/3)*RN(2)**3+(-RS(2,3)/6-RS(2,2)/3)*SN(2)**2*RN(2))*SLN(
     #2))
	CM(3,8)=
     #(((-RS(2,3)/3-RS(2,2)/6)*RN(2)**3+(-RS(2,3)/3-RS(2,2)/6)*S
     #N(2)**2*RN(2))*SLN(2)+((-RS(2,3)/3-RS(2,4)/6)*RN(3)**3+(-RS(2,3)/3
     #-RS(2,4)/6)*SN(3)**2*RN(3))*SLN(3))
	CM(3,11)=(((-RS(2,4)/3-RS(2,3)/6)
     #*RN(3)**3+(-RS(2,4)/3-RS(2,3)/6)*SN(3)**2*RN(3))*SLN(3)+((-RS(2,1)
     #/6-RS(2,4)/3)*RN(4)**3+(-RS(2,1)/6-RS(2,4)/3)*SN(4)**2*RN(4))*SLN(
     #4))
	CM(3,3) =
     #(((-RS(1,1)/24-RS(1,2)/24)*RN(1)**2+(RS(2,1)/12+RS(2,2)/12
     #)*SN(1)*RN(1))*SLN(1)**2+((RS(1,4)/24+RS(1,1)/24)*RN(4)**2+(-RS(2,
     #4)/12-RS(2,1)/12)*SN(4)*RN(4))*SLN(4)**2)
	CM(3,6)=
     #(((RS(1,1)/24+RS(1,2)/24)*RN(1)**2+(-RS(2,1)/12-RS(2,2)/12
     #)*SN(1)*RN(1))*SLN(1)**2+((-RS(1,3)/24-RS(1,2)/24)*RN(2)**2+(RS(2,
     #2)/12+RS(2,3)/12)*SN(2)*RN(2))*SLN(2)**2)
	CM(3,9)=(((RS(1,2)/24+RS(
     #1,3)/24)*RN(2)**2+(-RS(2,3)/12-RS(2,2)/12)*SN(2)*RN(2))*SLN(2)**2+
     #((-RS(1,4)/24-RS(1,3)/24)*RN(3)**2+(RS(2,3)/12+RS(2,4)/12)*SN(3)*R
     #N(3))*SLN(3)**2)
	CM(3,12)=(((RS(1,3)/24+RS(1,4)/24)*RN(3)**2+(-RS(2,
     #3)/12-RS(2,4)/12)*SN(3)*RN(3))*SLN(3)**2+((-RS(1,4)/24-RS(1,1)/24)
     #*RN(4)**2+(RS(2,1)/12+RS(2,4)/12)*SN(4)*RN(4))*SLN(4)**2)

	CM(4,1) =
     #(((RS(2,1)**2/4+RS(2,1)*RS(2,2)/6+RS(2,2)**2/12)*RN(1)**3+(SN
     #(1)**2*RS(2,2)**2/12+SN(1)**2*RS(2,1)**2/4+SN(1)**2*RS(2,1)*RS(2,2
     #)/6)*RN(1))*SLN(1)+((RS(2,1)**2/4+RS(2,4)*RS(2,1)/6+RS(2,4)**2/12)
     #*RN(4)**3+(SN(4)**2*RS(2,4)*RS(2,1)/6+SN(4)**2*RS(2,1)**2/4+SN(4)*
     #*2*RS(2,4)**2/12)*RN(4))*SLN(4))
	CM(4,4) =
     #(((RS(2,1)**2/12+RS(2,1)*RS(2,2)/6+RS(2,2)**2/4)*RN(1)**3+(SN
     #(1)**2*RS(2,1)**2/12+SN(1)**2*RS(2,1)*RS(2,2)/6+SN(1)**2*RS(2,2)**
     #2/4)*RN(1))*SLN(1)+((RS(2,3)**2/12+RS(2,2)*RS(2,3)/6+RS(2,2)**2/4)
     #*RN(2)**3+(SN(2)**2*RS(2,2)*RS(2,3)/6+SN(2)**2*RS(2,3)**2/12+SN(2)
     #**2*RS(2,2)**2/4)*RN(2))*SLN(2))
	CM(4,7)=
     #(((RS(2,3)**2/4+RS(2,2)**2/12+RS(2,2)*RS(2,3)/6)*RN(2)**3+(SN
     #(2)**2*RS(2,2)*RS(2,3)/6+SN(2)**2*RS(2,3)**2/4+SN(2)**2*RS(2,2)**2
     #/12)*RN(2))*SLN(2)+((RS(2,3)**2/4+RS(2,4)**2/12+RS(2,3)*RS(2,4)/6)
     #*RN(3)**3+(SN(3)**2*RS(2,4)**2/12+SN(3)**2*RS(2,3)**2/4+SN(3)**2*R
     #S(2,3)*RS(2,4)/6)*RN(3))*SLN(3))
	CM(4,10)=
     #(((RS(2,4)**2/4+RS(2,3)**2/12+RS(2,3)*RS(2,4)/6)*RN(3)**3+(SN
     #(3)**2*RS(2,3)*RS(2,4)/6+SN(3)**2*RS(2,4)**2/4+SN(3)**2*RS(2,3)**2
     #/12)*RN(3))*SLN(3)+((RS(2,1)**2/12+RS(2,4)**2/4+RS(2,4)*RS(2,1)/6)
     #*RN(4)**3+(SN(4)**2*RS(2,1)**2/12+SN(4)**2*RS(2,4)**2/4+SN(4)**2*R
     #S(2,4)*RS(2,1)/6)*RN(4))*SLN(4))
	CM(4,3) =
     #((-RS(2,1)**2/40-RS(2,1)*RS(2,2)/30-RS(2,2)**2/40)*RN(1)**
     #2*SLN(1)**2+(RS(2,4)**2/40+RS(2,1)**2/40+RS(2,4)*RS(2,1)/30)*RN(4)
     #**2*SLN(4)**2)
	CM(4,6)=((RS(2,1)*RS(2,2)/30+RS(2,1)**2/40+RS(2,2)**
     #2/40)*RN(1)**2*SLN(1)**2+(-RS(2,2)**2/40-RS(2,3)**2/40-RS(2,2)*RS(
     #2,3)/30)*RN(2)**2*SLN(2)**2)
	CM(4,9)=((RS(2,2)**2/40+RS(2,2)*RS(2,3
     #)/30+RS(2,3)**2/40)*RN(2)**2*SLN(2)**2+(-RS(2,3)**2/40-RS(2,3)*RS(
     #2,4)/30-RS(2,4)**2/40)*RN(3)**2*SLN(3)**2)
	CM(4,12)=((RS(2,4)**2/40+
     #RS(2,3)*RS(2,4)/30+RS(2,3)**2/40)*RN(3)**2*SLN(3)**2+(-RS(2,1)**2/
     #40-RS(2,4)*RS(2,1)/30-RS(2,4)**2/40)*RN(4)**2*SLN(4)**2)

	CM(5,2) =
     #((RN(1)**2*SN(1)/2+SN(1)**3/2)*SLN(1)+(SN(4)**3/2+RN(4)**2*SN
     #(4)/2)*SLN(4))
	CM(5,5) =((RN(1)**2*SN(1)/2+SN(1)**3/2)*SLN(1)+(SN(2)*
     #*3/2+RN(2)**2*SN(2)/2)*SLN(2))
	CM(5,8)=((SN(2)**3/2+RN(2)**2*SN(2)/2
     #)*SLN(2)+(RN(3)**2*SN(3)/2+SN(3)**3/2)*SLN(3))
	CM(5,11)=((RN(3)**2*SN
     #(3)/2+SN(3)**3/2)*SLN(3)+(SN(4)**3/2+RN(4)**2*SN(4)/2)*SLN(4))
	CM(5,3) =(-SN(1)**2*SLN(1)**2/12+SN(4)**2*SLN(4)**2/12)
	CM(5,6)=(-SN(2)**2*SLN(2)**2/12+SN(1)**2*SLN(1)**2/12)
	CM(5,9)=(SN(2)**2*SLN(2)**2/12-SN(3)**2*SLN(3)**2/12)
	CM(5,12)=(-SN(4)**2*SLN(4)**2/12+SN(3)**2*SLN(3)**2/12)

	CM(6,2) =
     #(((RS(1,1)/3+RS(1,2)/6)*SN(1)**3+(RN(1)**2*RS(1,2)/6+RN(1)**2
     #*RS(1,1)/3)*SN(1))*SLN(1)+((RS(1,1)/3+RS(1,4)/6)*SN(4)**3+(RN(4)**
     #2*RS(1,1)/3+RN(4)**2*RS(1,4)/6)*SN(4))*SLN(4))
	CM(6,5) =(((RS(1,1)/6+
     #RS(1,2)/3)*SN(1)**3+(RN(1)**2*RS(1,1)/6+RN(1)**2*RS(1,2)/3)*SN(1))
     #*SLN(1)+((RS(1,3)/6+RS(1,2)/3)*SN(2)**3+(RN(2)**2*RS(1,2)/3+RN(2)*
     #*2*RS(1,3)/6)*SN(2))*SLN(2))
	CM(6,8)=
     #(((RS(1,3)/3+RS(1,2)/6)*SN(2)**3+(RN(2)**2*RS(1,2)/6+RN(2)
     #**2*RS(1,3)/3)*SN(2))*SLN(2)+((RS(1,4)/6+RS(1,3)/3)*SN(3)**3+(RN(3
     #)**2*RS(1,4)/6+RN(3)**2*RS(1,3)/3)*SN(3))*SLN(3))
	CM(6,11)=(((RS(1,3)
     #/6+RS(1,4)/3)*SN(3)**3+(RN(3)**2*RS(1,3)/6+RN(3)**2*RS(1,4)/3)*SN(
     #3))*SLN(3)+((RS(1,4)/3+RS(1,1)/6)*SN(4)**3+(RN(4)**2*RS(1,1)/6+RN(
     #4)**2*RS(1,4)/3)*SN(4))*SLN(4))
	CM(6,3) =
     #((-RS(1,1)/24-RS(1,2)/24)*SN(1)**2*SLN(1)**2+(RS(1,4)/24+R
     #S(1,1)/24)*SN(4)**2*SLN(4)**2)
	CM(6,6)=((RS(1,1)/24+RS(1,2)/24)*SN(
     #1)**2*SLN(1)**2+(-RS(1,3)/24-RS(1,2)/24)*SN(2)**2*SLN(2)**2)
	CM(6,9)=
     #((RS(1,2)/24+RS(1,3)/24)*SN(2)**2*SLN(2)**2+(-RS(1,4)/24-RS(1,3)
     #/24)*SN(3)**2*SLN(3)**2)
	CM(6,12)=((RS(1,3)/24+RS(1,4)/24)*SN(3)**2*
     #SLN(3)**2+(-RS(1,4)/24-RS(1,1)/24)*SN(4)**2*SLN(4)**2)

	CM(7,1) =
     #(((-RS(1,1)/3-RS(1,2)/6)*SN(1)*RN(1)**2+(-RS(1,1)/3-RS(1,2)/6
     #)*SN(1)**3)*SLN(1)+((-RS(1,1)/3-RS(1,4)/6)*SN(4)*RN(4)**2+(-RS(1,1
     #)/3-RS(1,4)/6)*SN(4)**3)*SLN(4))
	CM(7,4) =(((-RS(1,1)/6-RS(1,2)/3)*SN
     #(1)*RN(1)**2+(-RS(1,1)/6-RS(1,2)/3)*SN(1)**3)*SLN(1)+((-RS(1,2)/3-
     #RS(1,3)/6)*SN(2)*RN(2)**2+(-RS(1,2)/3-RS(1,3)/6)*SN(2)**3)*SLN(2))
	CM(7,7)=
     #(((-RS(1,3)/3-RS(1,2)/6)*SN(2)*RN(2)**2+(-RS(1,3)/3-RS(1,2)
     #/6)*SN(2)**3)*SLN(2)+((-RS(1,3)/3-RS(1,4)/6)*SN(3)*RN(3)**2+(-RS(1
     #,3)/3-RS(1,4)/6)*SN(3)**3)*SLN(3))
	CM(7,10)=
     #(((-RS(1,4)/3-RS(1,3)/6)*SN(3)*RN(3)**2+(-RS(1,4)/3-RS(1,3
     #)/6)*SN(3)**3)*SLN(3)+((-RS(1,4)/3-RS(1,1)/6)*SN(4)*RN(4)**2+(-RS(
     #1,4)/3-RS(1,1)/6)*SN(4)**3)*SLN(4))
	CM(7,2) =
     #(((-RS(1,1)/3-RS(1,2)/6)*RN(1)**3+(RS(2,1)/3+RS(2,2)/6)*SN(1)
     #*RN(1)**2+(-RS(1,1)/3-RS(1,2)/6)*SN(1)**2*RN(1)+(RS(2,1)/3+RS(2,2)
     #/6)*SN(1)**3)*SLN(1)+((-RS(1,1)/3-RS(1,4)/6)*RN(4)**3+(RS(2,4)/6+R
     #S(2,1)/3)*SN(4)*RN(4)**2+(-RS(1,1)/3-RS(1,4)/6)*SN(4)**2*RN(4)+(RS
     #(2,4)/6+RS(2,1)/3)*SN(4)**3)*SLN(4))
	CM(7,5) =
     #(((-RS(1,1)/6-RS(1,2)/3)*RN(1)**3+(RS(2,2)/3+RS(2,1)/6)*SN(1)
     #*RN(1)**2+(-RS(1,1)/6-RS(1,2)/3)*SN(1)**2*RN(1)+(RS(2,2)/3+RS(2,1)
     #/6)*SN(1)**3)*SLN(1)+((-RS(1,2)/3-RS(1,3)/6)*RN(2)**3+(RS(2,2)/3+R
     #S(2,3)/6)*SN(2)*RN(2)**2+(-RS(1,2)/3-RS(1,3)/6)*SN(2)**2*RN(2)+(RS
     #(2,2)/3+RS(2,3)/6)*SN(2)**3)*SLN(2))
	CM(7,8)=
     #(((-RS(1,3)/3-RS(1,2)/6)*RN(2)**3+(RS(2,3)/3+RS(2,2)/6)*SN(2)
     #*RN(2)**2+(-RS(1,3)/3-RS(1,2)/6)*SN(2)**2*RN(2)+(RS(2,3)/3+RS(2,2)
     #/6)*SN(2)**3)*SLN(2)+((-RS(1,3)/3-RS(1,4)/6)*RN(3)**3+(RS(2,3)/3+R
     #S(2,4)/6)*SN(3)*RN(3)**2+(-RS(1,3)/3-RS(1,4)/6)*SN(3)**2*RN(3)+(RS
     #(2,3)/3+RS(2,4)/6)*SN(3)**3)*SLN(3))
	CM(7,11)=
     #(((-RS(1,4)/3-RS(1,3)/6)*RN(3)**3+(RS(2,3)/6+RS(2,4)/3)*SN(3)
     #*RN(3)**2+(-RS(1,4)/3-RS(1,3)/6)*SN(3)**2*RN(3)+(RS(2,3)/6+RS(2,4)
     #/3)*SN(3)**3)*SLN(3)+((-RS(1,4)/3-RS(1,1)/6)*RN(4)**3+(RS(2,4)/3+R
     #S(2,1)/6)*SN(4)*RN(4)**2+(-RS(1,4)/3-RS(1,1)/6)*SN(4)**2*RN(4)+(RS
     #(2,4)/3+RS(2,1)/6)*SN(4)**3)*SLN(4))
	CM(7,3) =
     #(((RS(1,2)/12+RS(1,1)/12)*SN(1)*RN(1)+(-RS(2,2)/24-RS(2,1)
     #/24)*SN(1)**2)*SLN(1)**2+((-RS(1,4)/12-RS(1,1)/12)*SN(4)*RN(4)+(RS
     #(2,1)/24+RS(2,4)/24)*SN(4)**2)*SLN(4)**2)
	CM(7,6)=
     #(((-RS(1,2)/12-RS(1,1)/12)*SN(1)*RN(1)+(RS(2,1)/24+RS(2,2)
     #/24)*SN(1)**2)*SLN(1)**2+((RS(1,3)/12+RS(1,2)/12)*SN(2)*RN(2)+(-RS
     #(2,3)/24-RS(2,2)/24)*SN(2)**2)*SLN(2)**2)
	CM(7,9)=(((-RS(1,2)/12-RS
     #(1,3)/12)*SN(2)*RN(2)+(RS(2,2)/24+RS(2,3)/24)*SN(2)**2)*SLN(2)**2+
     #((RS(1,4)/12+RS(1,3)/12)*SN(3)*RN(3)+(-RS(2,3)/24-RS(2,4)/24)*SN(3
     #)**2)*SLN(3)**2)
	CM(7,12)=(((-RS(1,3)/12-RS(1,4)/12)*SN(3)*RN(3)+(RS
     #(2,3)/24+RS(2,4)/24)*SN(3)**2)*SLN(3)**2+((RS(1,1)/12+RS(1,4)/12)*
     #SN(4)*RN(4)+(-RS(2,4)/24-RS(2,1)/24)*SN(4)**2)*SLN(4)**2)

	CM(8,2) =
     #(((RS(1,2)**2/12+RS(1,1)*RS(1,2)/6+RS(1,1)**2/4)*SN(1)**3+(RN
     #(1)**2*RS(1,1)*RS(1,2)/6+RN(1)**2*RS(1,1)**2/4+RN(1)**2*RS(1,2)**2
     #/12)*SN(1))*SLN(1)+((RS(1,4)*RS(1,1)/6+RS(1,1)**2/4+RS(1,4)**2/12)
     #*SN(4)**3+(RN(4)**2*RS(1,1)**2/4+RN(4)**2*RS(1,4)**2/12+RN(4)**2*R
     #S(1,4)*RS(1,1)/6)*SN(4))*SLN(4))
	CM(8,5) =
     #(((RS(1,1)**2/12+RS(1,1)*RS(1,2)/6+RS(1,2)**2/4)*SN(1)**3+(RN
     #(1)**2*RS(1,1)*RS(1,2)/6+RN(1)**2*RS(1,2)**2/4+RN(1)**2*RS(1,1)**2
     #/12)*SN(1))*SLN(1)+((RS(1,2)**2/4+RS(1,2)*RS(1,3)/6+RS(1,3)**2/12)
     #*SN(2)**3+(RN(2)**2*RS(1,2)**2/4+RN(2)**2*RS(1,3)**2/12+RN(2)**2*R
     #S(1,2)*RS(1,3)/6)*SN(2))*SLN(2))
	CM(8,8)=
     #(((RS(1,2)*RS(1,3)/6+RS(1,2)**2/12+RS(1,3)**2/4)*SN(2)**3+(RN
     #(2)**2*RS(1,2)**2/12+RN(2)**2*RS(1,2)*RS(1,3)/6+RN(2)**2*RS(1,3)**
     #2/4)*SN(2))*SLN(2)+((RS(1,4)**2/12+RS(1,3)**2/4+RS(1,3)*RS(1,4)/6)
     #*SN(3)**3+(RN(3)**2*RS(1,3)**2/4+RN(3)**2*RS(1,4)**2/12+RN(3)**2*R
     #S(1,3)*RS(1,4)/6)*SN(3))*SLN(3))
	CM(8,11)=
     #(((RS(1,3)**2/12+RS(1,4)**2/4+RS(1,3)*RS(1,4)/6)*SN(3)**3+(RN
     #(3)**2*RS(1,4)**2/4+RN(3)**2*RS(1,3)*RS(1,4)/6+RN(3)**2*RS(1,3)**2
     #/12)*SN(3))*SLN(3)+((RS(1,4)**2/4+RS(1,4)*RS(1,1)/6+RS(1,1)**2/12)
     #*SN(4)**3+(RN(4)**2*RS(1,4)**2/4+RN(4)**2*RS(1,1)**2/12+RN(4)**2*R
     #S(1,4)*RS(1,1)/6)*SN(4))*SLN(4))
	CM(8,3) =
     #((-RS(1,1)*RS(1,2)/30-RS(1,1)**2/40-RS(1,2)**2/40)*SN(1)**
     #2*SLN(1)**2+(RS(1,4)**2/40+RS(1,4)*RS(1,1)/30+RS(1,1)**2/40)*SN(4)
     #**2*SLN(4)**2)
	CM(8,6)=
     #((RS(1,1)*RS(1,2)/30+RS(1,2)**2/40+RS(1,1)**
     #2/40)*SN(1)**2*SLN(1)**2+(-RS(1,2)**2/40-RS(1,3)**2/40-RS(1,2)*RS(
     #1,3)/30)*SN(2)**2*SLN(2)**2)
	CM(8,9)=((RS(1,3)**2/40+RS(1,2)**2/40+
     #RS(1,2)*RS(1,3)/30)*SN(2)**2*SLN(2)**2+(-RS(1,4)**2/40-RS(1,3)*RS(
     #1,4)/30-RS(1,3)**2/40)*SN(3)**2*SLN(3)**2)
	CM(8,12)=((RS(1,4)**2/40+
     #RS(1,3)**2/40+RS(1,3)*RS(1,4)/30)*SN(3)**2*SLN(3)**2+(-RS(1,4)**2/
     #40-RS(1,4)*RS(1,1)/30-RS(1,1)**2/40)*SN(4)**2*SLN(4)**2)

	CM(9,1) =
     #((RN(1)**2*SN(1)/2+SN(1)**3/2)*SLN(1)+(SN(4)**3/2+RN(4)**2*SN
     #(4)/2)*SLN(4))
	CM(9,4) =((RN(1)**2*SN(1)/2+SN(1)**3/2)*SLN(1)+(SN(2)*
     #*3/2+RN(2)**2*SN(2)/2)*SLN(2))
	CM(9,7)=((SN(2)**3/2+RN(2)**2*SN(2)/2
     #)*SLN(2)+(RN(3)**2*SN(3)/2+SN(3)**3/2)*SLN(3))
	CM(9,10)=((RN(3)**2*SN(3)/2+SN(3)**3/2)*SLN(3)+
     #(SN(4)**3/2+RN(4)**2*SN(4)/2)*SLN(4))
	CM(9,2) =
     #((RN(1)*SN(1)**2/2+RN(1)**3/2)*SLN(1)+(RN(4)**3/2+RN(4)*SN(4)**
     #2/2)*SLN(4))
	CM(9,5) =((RN(1)*SN(1)**2/2+RN(1)**3/2)*SLN(1)+(RN(2)**3
     #/2+RN(2)*SN(2)**2/2)*SLN(2))
	CM(9,8)=
     #((RN(2)**3/2+RN(2)*SN(2)**2/2)*SLN(2)+(RN(3)*SN(3)**2/2+RN
     #(3)**3/2)*SLN(3))
	CM(9,11)=((RN(3)*SN(3)**2/2+RN(3)**3/2)*SLN(3)+(RN(
     #4)**3/2+RN(4)*SN(4)**2/2)*SLN(4))
	CM(9,3) =(-RN(1)*SN(1)*SLN(1)**2/6+RN(4)*SN(4)*SLN(4)**2/6)
	CM(9,6)=(-RN(2)*SN(2)*SLN(2)**2/6+RN(1)*SN(1)*SLN(1)**2/6)
	CM(9,9)=(-RN(3)*SN(3)*SLN(3)**2/6+RN(2)*SN(2)*SLN(2)**2/6)
	CM(9,12)=(RN(3)*SN(3)*SLN(3)**2/6-RN(4)*SN(4)*SLN(4)**2/6)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE ABCM(CM,RS,AM,ACM)
	IMPLICIT REAL*8 (A-H,O-Z)
C	------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR INVERSE(Am)*(Cm)=ACM
C
C	INPUT VARIABLES
C	CM(9,12)	= MEMBRANE AND DRILLING STRAIN PARAMETER
C	RS(2,4)		= ELEMENT NODAL COORDINATES IN LOCAL SPACE
C
C	LOCAL VARIABLE
C	AMI(9,9)	= INVERSE OF Am
C
C	OUTPUT VARIAVBLES
C	AM(9,9)		= MATRIX Am
C	ACM(9,12)	= INVERSE(Am)*(Cm)
C	------------------------------------------------------------
	DIMENSION AB1(4,4),CM(9,12),RS(2,4)
	DIMENSION AM(9,9),AMI(9,9)
	DIMENSION ACM(9,12)

C	------------------
C	ASSEMBLE MATRIX Am
C	------------------
	AM(1,1)=(-RS(2,4)/2+RS(2,2)/2)*RS(1,1)+(-RS(2,1)/2+RS(2,3)/2)*R
     #S(1,2)+(RS(2,4)/2-RS(2,2)/2)*RS(1,3)+(RS(2,1)/2-RS(2,3)/2)*RS(1,4)
	AM(1,2)=((-RS(2,4)/6+RS(2,2)/6)*RS(2,1)+RS(2,2)**2/6-RS(2,4)**2
     #/6)*RS(1,1)+(-RS(2,1)**2/6-RS(2,1)*RS(2,2)/6+RS(2,2)*RS(2,3)/6+RS(
     #2,3)**2/6)*RS(1,2)+(-RS(2,2)**2/6-RS(2,2)*RS(2,3)/6+RS(2,3)*RS(2,4
     #)/6+RS(2,4)**2/6)*RS(1,3)+(RS(2,1)**2/6+RS(2,1)*RS(2,4)/6-RS(2,3)*
     #*2/6-RS(2,3)*RS(2,4)/6)*RS(1,4)
	AM(1,3)=(-RS(2,4)/6+RS(2,2)/6)*RS(1,1)**2+((RS(2,2)/6-RS(2,1)/6
     #)*RS(1,2)+(-RS(2,4)/6+RS(2,1)/6)*RS(1,4))*RS(1,1)+(-RS(2,1)/6+RS(2
     #,3)/6)*RS(1,2)**2+(RS(2,3)/6-RS(2,2)/6)*RS(1,3)*RS(1,2)+(-RS(2,2)/
     #6+RS(2,4)/6)*RS(1,3)**2+(RS(2,4)/6-RS(2,3)/6)*RS(1,4)*RS(1,3)+(-RS
     #(2,3)/6+RS(2,1)/6)*RS(1,4)**2
	AM(1,4)=((-RS(2,4)/12+RS(2,2)/12)*RS(2,1)**2+(RS(2,2)**2/12-RS(
     #2,4)**2/12)*RS(2,1)+RS(2,2)**3/12-RS(2,4)**3/12)*RS(1,1)+(-RS(2,1)
     #**3/12-RS(2,2)*RS(2,1)**2/12-RS(2,2)**2*RS(2,1)/12+RS(2,3)*RS(2,2)
     #**2/12+RS(2,3)**2*RS(2,2)/12+RS(2,3)**3/12)*RS(1,2)+(-RS(2,2)**3/1
     #2-RS(2,3)*RS(2,2)**2/12-RS(2,3)**2*RS(2,2)/12+RS(2,4)*RS(2,3)**2/1
     #2+RS(2,4)**2*RS(2,3)/12+RS(2,4)**3/12)*RS(1,3)+(RS(2,1)**3/12+RS(2
     #,4)*RS(2,1)**2/12+RS(2,4)**2*RS(2,1)/12-RS(2,3)**3/12-RS(2,4)*RS(2
     #,3)**2/12-RS(2,4)**2*RS(2,3)/12)*RS(1,4)

	AM(2,2)= AM(1,4)
	AM(2,3)=((-RS(2,4)/12+RS(2,2)/12)*RS(2,1)+RS(2,2)**2/24-RS(2,4)
     #**2/24)*RS(1,1)**2+((-RS(2,1)**2/12+RS(2,2)**2/12)*RS(1,2)+(RS(2,1
     #)**2/12-RS(2,4)**2/12)*RS(1,4))*RS(1,1)+(-RS(2,1)**2/24-RS(2,1)*RS
     #(2,2)/12+RS(2,2)*RS(2,3)/12+RS(2,3)**2/24)*RS(1,2)**2+(-RS(2,2)**2
     #/12+RS(2,3)**2/12)*RS(1,3)*RS(1,2)+(-RS(2,2)**2/24-RS(2,2)*RS(2,3)
     #/12+RS(2,3)*RS(2,4)/12+RS(2,4)**2/24)*RS(1,3)**2+(-RS(2,3)**2/12+R
     #S(2,4)**2/12)*RS(1,4)*RS(1,3)+(RS(2,1)**2/24+RS(2,1)*RS(2,4)/12-RS
     #(2,3)**2/24-RS(2,3)*RS(2,4)/12)*RS(1,4)**2
      s1 = ((RS(2,2)/20-RS(2,4)/20)*RS(2,1)**3+(RS(2,2)**2/20-RS(2,4)**2
     #/20)*RS(2,1)**2+(RS(2,2)**3/20-RS(2,4)**3/20)*RS(2,1)+RS(2,2)**4/2
     #0-RS(2,4)**4/20)*RS(1,1)+(-RS(2,1)**4/20-RS(2,1)**3*RS(2,2)/20-RS(
     #2,1)**2*RS(2,2)**2/20-RS(2,1)*RS(2,2)**3/20+RS(2,3)*RS(2,2)**3/20+
     #RS(2,3)**2*RS(2,2)**2/20+RS(2,3)**3*RS(2,2)/20+RS(2,3)**4/20)*RS(1
     #,2)
      AM(2,4) = s1+(-RS(2,2)**4/20-RS(2,3)*RS(2,2)**3/20-RS(2,3)**2*RS(
     #2,2)**2/20-RS(2,3)**3*RS(2,2)/20+RS(2,3)**3*RS(2,4)/20+RS(2,3)**2*
     #RS(2,4)**2/20+RS(2,3)*RS(2,4)**3/20+RS(2,4)**4/20)*RS(1,3)+(RS(2,1
     #)**4/20+RS(2,4)*RS(2,1)**3/20+RS(2,4)**2*RS(2,1)**2/20+RS(2,4)**3*
     #RS(2,1)/20-RS(2,3)**4/20-RS(2,3)**3*RS(2,4)/20-RS(2,3)**2*RS(2,4)*
     #*2/20-RS(2,3)*RS(2,4)**3/20)*RS(1,4)

	AM(5,8)=(-RS(2,4)/12+RS(2,2)/12)*RS(1,1)**3+((-RS(2,1)/12+RS(2,
     #2)/12)*RS(1,2)+(-RS(2,4)/12+RS(2,1)/12)*RS(1,4))*RS(1,1)**2+((-RS(
     #2,1)/12+RS(2,2)/12)*RS(1,2)**2+(-RS(2,4)/12+RS(2,1)/12)*RS(1,4)**2
     #)*RS(1,1)+(RS(2,3)/12-RS(2,1)/12)*RS(1,2)**3+(RS(2,3)/12-RS(2,2)/1
     #2)*RS(1,3)*RS(1,2)**2+(RS(2,3)/12-RS(2,2)/12)*RS(1,3)**2*RS(1,2)+(
     #RS(2,4)/12-RS(2,2)/12)*RS(1,3)**3+(RS(2,4)/12-RS(2,3)/12)*RS(1,4)*
     #RS(1,3)**2+(RS(2,4)/12-RS(2,3)/12)*RS(1,4)**2*RS(1,3)+(RS(2,1)/12-
     #RS(2,3)/12)*RS(1,4)**3

	AM(3,3)=AM(1,4)+AM(5,8)
      s1 = ((RS(2,2)/20-RS(2,4)/20)*RS(2,1)**2+(RS(2,2)**2/30-RS(2,4)**2
     #/30)*RS(2,1)+RS(2,2)**3/60-RS(2,4)**3/60)*RS(1,1)**2+((-RS(2,1)**3
     #/20-RS(2,1)**2*RS(2,2)/60+RS(2,1)*RS(2,2)**2/60+RS(2,2)**3/20)*RS(
     #1,2)+(RS(2,1)**3/20+RS(2,1)**2*RS(2,4)/60-RS(2,1)*RS(2,4)**2/60-RS
     #(2,4)**3/20)*RS(1,4))*RS(1,1)+(-RS(2,1)**3/60-RS(2,1)**2*RS(2,2)/3
     #0-RS(2,1)*RS(2,2)**2/20+RS(2,2)**2*RS(2,3)/20+RS(2,2)*RS(2,3)**2/3
     #0+RS(2,3)**3/60)*RS(1,2)**2
      AM(3,4) = s1+(-RS(2,2)**3/20-RS(2,2)**2*RS(2,3)/60+RS(2,2)*RS(2,3
     #)**2/60+RS(2,3)**3/20)*RS(1,3)*RS(1,2)+(-RS(2,2)**3/60-RS(2,2)**2*
     #RS(2,3)/30-RS(2,2)*RS(2,3)**2/20+RS(2,3)**2*RS(2,4)/20+RS(2,3)*RS(
     #2,4)**2/30+RS(2,4)**3/60)*RS(1,3)**2+(-RS(2,3)**3/20-RS(2,3)**2*RS
     #(2,4)/60+RS(2,3)*RS(2,4)**2/60+RS(2,4)**3/20)*RS(1,4)*RS(1,3)+(RS(
     #2,1)**3/60+RS(2,1)**2*RS(2,4)/30+RS(2,1)*RS(2,4)**2/20-RS(2,3)**3/
     #60-RS(2,3)**2*RS(2,4)/30-RS(2,3)*RS(2,4)**2/20)*RS(1,4)**2
	AM(3,7)=AM(2,3)
	AM(3,9)=-AM(1,2)

      s1 = ((RS(2,2)/30-RS(2,4)/30)*RS(2,1)**4+(RS(2,2)**2/30-RS(2,4)**2
     #/30)*RS(2,1)**3+(RS(2,2)**3/30-RS(2,4)**3/30)*RS(2,1)**2+(RS(2,2)*
     #*4/30-RS(2,4)**4/30)*RS(2,1)+RS(2,2)**5/30-RS(2,4)**5/30)*RS(1,1)+
     #(-RS(2,1)**5/30-RS(2,2)*RS(2,1)**4/30-RS(2,2)**2*RS(2,1)**3/30-RS(
     #2,2)**3*RS(2,1)**2/30-RS(2,2)**4*RS(2,1)/30+RS(2,3)*RS(2,2)**4/30+
     #RS(2,3)**2*RS(2,2)**3/30+RS(2,3)**3*RS(2,2)**2/30+RS(2,3)**4*RS(2,
     #2)/30+RS(2,3)**5/30)*RS(1,2)
      AM(4,4) = s1+(-RS(2,2)**5/30-RS(2,3)*RS(2,2)**4/30-RS(2,3)**2*RS(
     #2,2)**3/30-RS(2,3)**3*RS(2,2)**2/30-RS(2,3)**4*RS(2,2)/30+RS(2,4)*
     #RS(2,3)**4/30+RS(2,4)**2*RS(2,3)**3/30+RS(2,4)**3*RS(2,3)**2/30+RS
     #(2,4)**4*RS(2,3)/30+RS(2,4)**5/30)*RS(1,3)+(RS(2,1)**5/30+RS(2,4)*
     #RS(2,1)**4/30+RS(2,4)**2*RS(2,1)**3/30+RS(2,4)**3*RS(2,1)**2/30+RS
     #(2,4)**4*RS(2,1)/30-RS(2,3)**5/30-RS(2,4)*RS(2,3)**4/30-RS(2,4)**2
     #*RS(2,3)**3/30-RS(2,4)**3*RS(2,3)**2/30-RS(2,4)**4*RS(2,3)/30)*RS(
     #1,4)

	AM(5,5)=AM(1,1)
	AM(5,6)=AM(1,3)
	AM(5,7)=AM(1,2)

	AM(6,6)=AM(5,8)
	AM(6,7)=AM(2,3)
      s1 = (RS(2,2)/20-RS(2,4)/20)*RS(1,1)**4+((-RS(2,1)/20+RS(2,2)/20)*
     #RS(1,2)+(-RS(2,4)/20+RS(2,1)/20)*RS(1,4))*RS(1,1)**3+((-RS(2,1)/20
     #+RS(2,2)/20)*RS(1,2)**2+(-RS(2,4)/20+RS(2,1)/20)*RS(1,4)**2)*RS(1,
     #1)**2+((-RS(2,1)/20+RS(2,2)/20)*RS(1,2)**3+(-RS(2,4)/20+RS(2,1)/20
     #)*RS(1,4)**3)*RS(1,1)+(RS(2,3)/20-RS(2,1)/20)*RS(1,2)**4+(-RS(2,2)
     #/20+RS(2,3)/20)*RS(1,3)*RS(1,2)**3
      AM(6,8) = s1+(-RS(2,2)/20+RS(2,3)/20)*RS(1,3)**2*RS(1,2)**2+(-RS(
     #2,2)/20+RS(2,3)/20)*RS(1,3)**3*RS(1,2)+(-RS(2,2)/20+RS(2,4)/20)*RS
     #(1,3)**4+(RS(2,4)/20-RS(2,3)/20)*RS(1,4)*RS(1,3)**3+(RS(2,4)/20-RS
     #(2,3)/20)*RS(1,4)**2*RS(1,3)**2+(RS(2,4)/20-RS(2,3)/20)*RS(1,4)**3
     #*RS(1,3)+(-RS(2,3)/20+RS(2,1)/20)*RS(1,4)**4

	AM(7,7)=AM(3,3)
      s1 = ((RS(2,2)/20-RS(2,4)/20)*RS(2,1)+RS(2,2)**2/60-RS(2,4)**2/60)
     #*RS(1,1)**3+((-RS(2,1)**2/20+RS(2,1)*RS(2,2)/60+RS(2,2)**2/30)*RS(
     #1,2)+(RS(2,1)**2/20-RS(2,1)*RS(2,4)/60-RS(2,4)**2/30)*RS(1,4))*RS(
     #1,1)**2+((-RS(2,1)**2/30-RS(2,1)*RS(2,2)/60+RS(2,2)**2/20)*RS(1,2)
     #**2+(RS(2,1)**2/30+RS(2,1)*RS(2,4)/60-RS(2,4)**2/20)*RS(1,4)**2)*R
     #S(1,1)+(-RS(2,1)**2/60-RS(2,1)*RS(2,2)/20+RS(2,2)*RS(2,3)/20+RS(2,
     #3)**2/60)*RS(1,2)**3+(-RS(2,2)**2/20+RS(2,2)*RS(2,3)/60+RS(2,3)**2
     #/30)*RS(1,3)*RS(1,2)**2
      AM(7,8) = s1+(-RS(2,2)**2/30-RS(2,2)*RS(2,3)/60+RS(2,3)**2/20)*RS
     #(1,3)**2*RS(1,2)+(-RS(2,2)**2/60-RS(2,2)*RS(2,3)/20+RS(2,3)*RS(2,4
     #)/20+RS(2,4)**2/60)*RS(1,3)**3+(-RS(2,3)**2/20+RS(2,3)*RS(2,4)/60+
     #RS(2,4)**2/30)*RS(1,4)*RS(1,3)**2+(-RS(2,3)**2/30-RS(2,3)*RS(2,4)/
     #60+RS(2,4)**2/20)*RS(1,4)**2*RS(1,3)+(RS(2,1)**2/60+RS(2,1)*RS(2,4
     #)/20-RS(2,3)**2/60-RS(2,3)*RS(2,4)/20)*RS(1,4)**3
	AM(7,9)=-AM(1,3)

      s1 = (RS(2,2)/30-RS(2,4)/30)*RS(1,1)**5+((-RS(2,1)/30+RS(2,2)/30)*
     #RS(1,2)+(-RS(2,4)/30+RS(2,1)/30)*RS(1,4))*RS(1,1)**4+((-RS(2,1)/30
     #+RS(2,2)/30)*RS(1,2)**2+(-RS(2,4)/30+RS(2,1)/30)*RS(1,4)**2)*RS(1,
     #1)**3+((-RS(2,1)/30+RS(2,2)/30)*RS(1,2)**3+(-RS(2,4)/30+RS(2,1)/30
     #)*RS(1,4)**3)*RS(1,1)**2+((-RS(2,1)/30+RS(2,2)/30)*RS(1,2)**4+(-RS
     #(2,4)/30+RS(2,1)/30)*RS(1,4)**4)*RS(1,1)+(-RS(2,1)/30+RS(2,3)/30)*
     #RS(1,2)**5+(-RS(2,2)/30+RS(2,3)/30)*RS(1,3)*RS(1,2)**4+(-RS(2,2)/3
     #0+RS(2,3)/30)*RS(1,3)**2*RS(1,2)**3
      AM(8,8) = s1+(-RS(2,2)/30+RS(2,3)/30)*RS(1,3)**3*RS(1,2)**2+(-RS(
     #2,2)/30+RS(2,3)/30)*RS(1,3)**4*RS(1,2)+(-RS(2,2)/30+RS(2,4)/30)*RS
     #(1,3)**5+(RS(2,4)/30-RS(2,3)/30)*RS(1,4)*RS(1,3)**4+(RS(2,4)/30-RS
     #(2,3)/30)*RS(1,4)**2*RS(1,3)**3+(RS(2,4)/30-RS(2,3)/30)*RS(1,4)**3
     #*RS(1,3)**2+(RS(2,4)/30-RS(2,3)/30)*RS(1,4)**4*RS(1,3)+(RS(2,1)/30
     #-RS(2,3)/30)*RS(1,4)**5

	AM(9,9)=AM(1,1)

C	-----------------
C	FILTER ZERO TERMS
C	-----------------
	AMAX=0.0
	DO 10 I=1,9
	  DO 10 J=I,9
	    IF (DABS(AM(I,J)).GE.AMAX) AMAX=DABS(AM(I,J))
10	CONTINUE
	DO 20 I=1,9
	  DO 20 J=I,9
	    RES=DABS(AM(I,J)/AMAX)
	    IF (RES.LE.1.0E-8) AM(I,J)=0.0
20	CONTINUE

C	----------------------
C	FILL IN LOWER TRIANGLE
C	----------------------
	DO 30 I=1,9
	  DO 30 J=I,9
	    AM(J,I)=AM(I,J)
30	CONTINUE

C	--------------------
C	SOLVE FOR Am INVERSE
C	--------------------
	CALL INVMAT(AM,AMI,9)

C	--------------------------------
C	SOLVE FOR INVERSE(Am)*(Cm) - ACM
C	--------------------------------
      ACM(1,1) = AMI(1,1)*CM(1,1)+AMI(1,2)*CM(2,1)+AMI(1,3)*CM(3,1)+AMI(
     #1,4)*CM(4,1)+AMI(1,7)*CM(7,1)+AMI(1,9)*CM(9,1)
      ACM(1,2) = AMI(1,3)*CM(3,2)+AMI(1,5)*CM(5,2)+AMI(1,6)*CM(6,2)+AMI(
     #1,7)*CM(7,2)+AMI(1,8)*CM(8,2)+AMI(1,9)*CM(9,2)
      ACM(1,3) = AMI(1,1)*CM(1,3)+AMI(1,2)*CM(2,3)+AMI(1,3)*CM(3,3)+AMI(
     #1,4)*CM(4,3)+AMI(1,5)*CM(5,3)+AMI(1,6)*CM(6,3)+AMI(1,7)*CM(7,3)+AM
     #I(1,8)*CM(8,3)+AMI(1,9)*CM(9,3)
      ACM(1,4) = AMI(1,1)*CM(1,4)+AMI(1,2)*CM(2,4)+AMI(1,3)*CM(3,4)+AMI(
     #1,4)*CM(4,4)+AMI(1,7)*CM(7,4)+AMI(1,9)*CM(9,4)
      ACM(1,5) = AMI(1,3)*CM(3,5)+AMI(1,5)*CM(5,5)+AMI(1,6)*CM(6,5)+AMI(
     #1,7)*CM(7,5)+AMI(1,8)*CM(8,5)+AMI(1,9)*CM(9,5)
      ACM(1,6) = AMI(1,1)*CM(1,6)+AMI(1,2)*CM(2,6)+AMI(1,3)*CM(3,6)+AMI(
     #1,4)*CM(4,6)+AMI(1,5)*CM(5,6)+AMI(1,6)*CM(6,6)+AMI(1,7)*CM(7,6)+AM
     #I(1,8)*CM(8,6)+AMI(1,9)*CM(9,6)
      ACM(1,7) = AMI(1,1)*CM(1,7)+AMI(1,2)*CM(2,7)+AMI(1,3)*CM(3,7)+AMI(
     #1,4)*CM(4,7)+AMI(1,7)*CM(7,7)+AMI(1,9)*CM(9,7)
      ACM(1,8) = AMI(1,3)*CM(3,8)+AMI(1,5)*CM(5,8)+AMI(1,6)*CM(6,8)+AMI(
     #1,7)*CM(7,8)+AMI(1,8)*CM(8,8)+AMI(1,9)*CM(9,8)
      ACM(1,9) = AMI(1,1)*CM(1,9)+AMI(1,2)*CM(2,9)+AMI(1,3)*CM(3,9)+AMI(
     #1,4)*CM(4,9)+AMI(1,5)*CM(5,9)+AMI(1,6)*CM(6,9)+AMI(1,7)*CM(7,9)+AM
     #I(1,8)*CM(8,9)+AMI(1,9)*CM(9,9)
      ACM(1,10) = AMI(1,1)*CM(1,10)+AMI(1,2)*CM(2,10)+AMI(1,3)*CM(3,10)+
     #AMI(1,4)*CM(4,10)+AMI(1,7)*CM(7,10)+AMI(1,9)*CM(9,10)
      ACM(1,11) = AMI(1,3)*CM(3,11)+AMI(1,5)*CM(5,11)+AMI(1,6)*CM(6,11)+
     #AMI(1,7)*CM(7,11)+AMI(1,8)*CM(8,11)+AMI(1,9)*CM(9,11)
      ACM(1,12) = AMI(1,1)*CM(1,12)+AMI(1,2)*CM(2,12)+AMI(1,3)*CM(3,12)+
     #AMI(1,4)*CM(4,12)+AMI(1,5)*CM(5,12)+AMI(1,6)*CM(6,12)+AMI(1,7)*CM(
     #7,12)+AMI(1,8)*CM(8,12)+AMI(1,9)*CM(9,12)
      ACM(2,1) = AMI(2,1)*CM(1,1)+AMI(2,2)*CM(2,1)+AMI(2,3)*CM(3,1)+AMI(
     #2,4)*CM(4,1)+AMI(2,7)*CM(7,1)+AMI(2,9)*CM(9,1)
      ACM(2,2) = AMI(2,3)*CM(3,2)+AMI(2,5)*CM(5,2)+AMI(2,6)*CM(6,2)+AMI(
     #2,7)*CM(7,2)+AMI(2,8)*CM(8,2)+AMI(2,9)*CM(9,2)
      ACM(2,3) = AMI(2,1)*CM(1,3)+AMI(2,2)*CM(2,3)+AMI(2,3)*CM(3,3)+AMI(
     #2,4)*CM(4,3)+AMI(2,5)*CM(5,3)+AMI(2,6)*CM(6,3)+AMI(2,7)*CM(7,3)+AM
     #I(2,8)*CM(8,3)+AMI(2,9)*CM(9,3)
      ACM(2,4) = AMI(2,1)*CM(1,4)+AMI(2,2)*CM(2,4)+AMI(2,3)*CM(3,4)+AMI(
     #2,4)*CM(4,4)+AMI(2,7)*CM(7,4)+AMI(2,9)*CM(9,4)
      ACM(2,5) = AMI(2,3)*CM(3,5)+AMI(2,5)*CM(5,5)+AMI(2,6)*CM(6,5)+AMI(
     #2,7)*CM(7,5)+AMI(2,8)*CM(8,5)+AMI(2,9)*CM(9,5)
      ACM(2,6) = AMI(2,1)*CM(1,6)+AMI(2,2)*CM(2,6)+AMI(2,3)*CM(3,6)+AMI(
     #2,4)*CM(4,6)+AMI(2,5)*CM(5,6)+AMI(2,6)*CM(6,6)+AMI(2,7)*CM(7,6)+AM
     #I(2,8)*CM(8,6)+AMI(2,9)*CM(9,6)
      ACM(2,7) = AMI(2,1)*CM(1,7)+AMI(2,2)*CM(2,7)+AMI(2,3)*CM(3,7)+AMI(
     #2,4)*CM(4,7)+AMI(2,7)*CM(7,7)+AMI(2,9)*CM(9,7)
      ACM(2,8) = AMI(2,3)*CM(3,8)+AMI(2,5)*CM(5,8)+AMI(2,6)*CM(6,8)+AMI(
     #2,7)*CM(7,8)+AMI(2,8)*CM(8,8)+AMI(2,9)*CM(9,8)
      ACM(2,9) = AMI(2,1)*CM(1,9)+AMI(2,2)*CM(2,9)+AMI(2,3)*CM(3,9)+AMI(
     #2,4)*CM(4,9)+AMI(2,5)*CM(5,9)+AMI(2,6)*CM(6,9)+AMI(2,7)*CM(7,9)+AM
     #I(2,8)*CM(8,9)+AMI(2,9)*CM(9,9)
      ACM(2,10) = AMI(2,1)*CM(1,10)+AMI(2,2)*CM(2,10)+AMI(2,3)*CM(3,10)+
     #AMI(2,4)*CM(4,10)+AMI(2,7)*CM(7,10)+AMI(2,9)*CM(9,10)
      ACM(2,11) = AMI(2,3)*CM(3,11)+AMI(2,5)*CM(5,11)+AMI(2,6)*CM(6,11)+
     #AMI(2,7)*CM(7,11)+AMI(2,8)*CM(8,11)+AMI(2,9)*CM(9,11)
      ACM(2,12) = AMI(2,1)*CM(1,12)+AMI(2,2)*CM(2,12)+AMI(2,3)*CM(3,12)+
     #AMI(2,4)*CM(4,12)+AMI(2,5)*CM(5,12)+AMI(2,6)*CM(6,12)+AMI(2,7)*CM(
     #7,12)+AMI(2,8)*CM(8,12)+AMI(2,9)*CM(9,12)
      ACM(3,1) = AMI(3,1)*CM(1,1)+AMI(3,2)*CM(2,1)+AMI(3,3)*CM(3,1)+AMI(
     #3,4)*CM(4,1)+AMI(3,7)*CM(7,1)+AMI(3,9)*CM(9,1)
      ACM(3,2) = AMI(3,3)*CM(3,2)+AMI(3,5)*CM(5,2)+AMI(3,6)*CM(6,2)+AMI(
     #3,7)*CM(7,2)+AMI(3,8)*CM(8,2)+AMI(3,9)*CM(9,2)
      ACM(3,3) = AMI(3,1)*CM(1,3)+AMI(3,2)*CM(2,3)+AMI(3,3)*CM(3,3)+AMI(
     #3,4)*CM(4,3)+AMI(3,5)*CM(5,3)+AMI(3,6)*CM(6,3)+AMI(3,7)*CM(7,3)+AM
     #I(3,8)*CM(8,3)+AMI(3,9)*CM(9,3)
      ACM(3,4) = AMI(3,1)*CM(1,4)+AMI(3,2)*CM(2,4)+AMI(3,3)*CM(3,4)+AMI(
     #3,4)*CM(4,4)+AMI(3,7)*CM(7,4)+AMI(3,9)*CM(9,4)
      ACM(3,5) = AMI(3,3)*CM(3,5)+AMI(3,5)*CM(5,5)+AMI(3,6)*CM(6,5)+AMI(
     #3,7)*CM(7,5)+AMI(3,8)*CM(8,5)+AMI(3,9)*CM(9,5)
      ACM(3,6) = AMI(3,1)*CM(1,6)+AMI(3,2)*CM(2,6)+AMI(3,3)*CM(3,6)+AMI(
     #3,4)*CM(4,6)+AMI(3,5)*CM(5,6)+AMI(3,6)*CM(6,6)+AMI(3,7)*CM(7,6)+AM
     #I(3,8)*CM(8,6)+AMI(3,9)*CM(9,6)
      ACM(3,7) = AMI(3,1)*CM(1,7)+AMI(3,2)*CM(2,7)+AMI(3,3)*CM(3,7)+AMI(
     #3,4)*CM(4,7)+AMI(3,7)*CM(7,7)+AMI(3,9)*CM(9,7)
      ACM(3,8) = AMI(3,3)*CM(3,8)+AMI(3,5)*CM(5,8)+AMI(3,6)*CM(6,8)+AMI(
     #3,7)*CM(7,8)+AMI(3,8)*CM(8,8)+AMI(3,9)*CM(9,8)
      ACM(3,9) = AMI(3,1)*CM(1,9)+AMI(3,2)*CM(2,9)+AMI(3,3)*CM(3,9)+AMI(
     #3,4)*CM(4,9)+AMI(3,5)*CM(5,9)+AMI(3,6)*CM(6,9)+AMI(3,7)*CM(7,9)+AM
     #I(3,8)*CM(8,9)+AMI(3,9)*CM(9,9)
      ACM(3,10) = AMI(3,1)*CM(1,10)+AMI(3,2)*CM(2,10)+AMI(3,3)*CM(3,10)+
     #AMI(3,4)*CM(4,10)+AMI(3,7)*CM(7,10)+AMI(3,9)*CM(9,10)
      ACM(3,11) = AMI(3,3)*CM(3,11)+AMI(3,5)*CM(5,11)+AMI(3,6)*CM(6,11)+
     #AMI(3,7)*CM(7,11)+AMI(3,8)*CM(8,11)+AMI(3,9)*CM(9,11)
      ACM(3,12) = AMI(3,1)*CM(1,12)+AMI(3,2)*CM(2,12)+AMI(3,3)*CM(3,12)+
     #AMI(3,4)*CM(4,12)+AMI(3,5)*CM(5,12)+AMI(3,6)*CM(6,12)+AMI(3,7)*CM(
     #7,12)+AMI(3,8)*CM(8,12)+AMI(3,9)*CM(9,12)
      ACM(4,1) = AMI(4,1)*CM(1,1)+AMI(4,2)*CM(2,1)+AMI(4,3)*CM(3,1)+AMI(
     #4,4)*CM(4,1)+AMI(4,7)*CM(7,1)+AMI(4,9)*CM(9,1)
      ACM(4,2) = AMI(4,3)*CM(3,2)+AMI(4,5)*CM(5,2)+AMI(4,6)*CM(6,2)+AMI(
     #4,7)*CM(7,2)+AMI(4,8)*CM(8,2)+AMI(4,9)*CM(9,2)
      ACM(4,3) = AMI(4,1)*CM(1,3)+AMI(4,2)*CM(2,3)+AMI(4,3)*CM(3,3)+AMI(
     #4,4)*CM(4,3)+AMI(4,5)*CM(5,3)+AMI(4,6)*CM(6,3)+AMI(4,7)*CM(7,3)+AM
     #I(4,8)*CM(8,3)+AMI(4,9)*CM(9,3)
      ACM(4,4) = AMI(4,1)*CM(1,4)+AMI(4,2)*CM(2,4)+AMI(4,3)*CM(3,4)+AMI(
     #4,4)*CM(4,4)+AMI(4,7)*CM(7,4)+AMI(4,9)*CM(9,4)
      ACM(4,5) = AMI(4,3)*CM(3,5)+AMI(4,5)*CM(5,5)+AMI(4,6)*CM(6,5)+AMI(
     #4,7)*CM(7,5)+AMI(4,8)*CM(8,5)+AMI(4,9)*CM(9,5)
      ACM(4,6) = AMI(4,1)*CM(1,6)+AMI(4,2)*CM(2,6)+AMI(4,3)*CM(3,6)+AMI(
     #4,4)*CM(4,6)+AMI(4,5)*CM(5,6)+AMI(4,6)*CM(6,6)+AMI(4,7)*CM(7,6)+AM
     #I(4,8)*CM(8,6)+AMI(4,9)*CM(9,6)
      ACM(4,7) = AMI(4,1)*CM(1,7)+AMI(4,2)*CM(2,7)+AMI(4,3)*CM(3,7)+AMI(
     #4,4)*CM(4,7)+AMI(4,7)*CM(7,7)+AMI(4,9)*CM(9,7)
      ACM(4,8) = AMI(4,3)*CM(3,8)+AMI(4,5)*CM(5,8)+AMI(4,6)*CM(6,8)+AMI(
     #4,7)*CM(7,8)+AMI(4,8)*CM(8,8)+AMI(4,9)*CM(9,8)
      ACM(4,9) = AMI(4,1)*CM(1,9)+AMI(4,2)*CM(2,9)+AMI(4,3)*CM(3,9)+AMI(
     #4,4)*CM(4,9)+AMI(4,5)*CM(5,9)+AMI(4,6)*CM(6,9)+AMI(4,7)*CM(7,9)+AM
     #I(4,8)*CM(8,9)+AMI(4,9)*CM(9,9)
      ACM(4,10) = AMI(4,1)*CM(1,10)+AMI(4,2)*CM(2,10)+AMI(4,3)*CM(3,10)+
     #AMI(4,4)*CM(4,10)+AMI(4,7)*CM(7,10)+AMI(4,9)*CM(9,10)
      ACM(4,11) = AMI(4,3)*CM(3,11)+AMI(4,5)*CM(5,11)+AMI(4,6)*CM(6,11)+
     #AMI(4,7)*CM(7,11)+AMI(4,8)*CM(8,11)+AMI(4,9)*CM(9,11)
      ACM(4,12) = AMI(4,1)*CM(1,12)+AMI(4,2)*CM(2,12)+AMI(4,3)*CM(3,12)+
     #AMI(4,4)*CM(4,12)+AMI(4,5)*CM(5,12)+AMI(4,6)*CM(6,12)+AMI(4,7)*CM(
     #7,12)+AMI(4,8)*CM(8,12)+AMI(4,9)*CM(9,12)
      ACM(5,1) = AMI(5,1)*CM(1,1)+AMI(5,2)*CM(2,1)+AMI(5,3)*CM(3,1)+AMI(
     #5,4)*CM(4,1)+AMI(5,7)*CM(7,1)+AMI(5,9)*CM(9,1)
      ACM(5,2) = AMI(5,3)*CM(3,2)+AMI(5,5)*CM(5,2)+AMI(5,6)*CM(6,2)+AMI(
     #5,7)*CM(7,2)+AMI(5,8)*CM(8,2)+AMI(5,9)*CM(9,2)
      ACM(5,3) = AMI(5,1)*CM(1,3)+AMI(5,2)*CM(2,3)+AMI(5,3)*CM(3,3)+AMI(
     #5,4)*CM(4,3)+AMI(5,5)*CM(5,3)+AMI(5,6)*CM(6,3)+AMI(5,7)*CM(7,3)+AM
     #I(5,8)*CM(8,3)+AMI(5,9)*CM(9,3)
      ACM(5,4) = AMI(5,1)*CM(1,4)+AMI(5,2)*CM(2,4)+AMI(5,3)*CM(3,4)+AMI(
     #5,4)*CM(4,4)+AMI(5,7)*CM(7,4)+AMI(5,9)*CM(9,4)
      ACM(5,5) = AMI(5,3)*CM(3,5)+AMI(5,5)*CM(5,5)+AMI(5,6)*CM(6,5)+AMI(
     #5,7)*CM(7,5)+AMI(5,8)*CM(8,5)+AMI(5,9)*CM(9,5)
      ACM(5,6) = AMI(5,1)*CM(1,6)+AMI(5,2)*CM(2,6)+AMI(5,3)*CM(3,6)+AMI(
     #5,4)*CM(4,6)+AMI(5,5)*CM(5,6)+AMI(5,6)*CM(6,6)+AMI(5,7)*CM(7,6)+AM
     #I(5,8)*CM(8,6)+AMI(5,9)*CM(9,6)
      ACM(5,7) = AMI(5,1)*CM(1,7)+AMI(5,2)*CM(2,7)+AMI(5,3)*CM(3,7)+AMI(
     #5,4)*CM(4,7)+AMI(5,7)*CM(7,7)+AMI(5,9)*CM(9,7)
      ACM(5,8) = AMI(5,3)*CM(3,8)+AMI(5,5)*CM(5,8)+AMI(5,6)*CM(6,8)+AMI(
     #5,7)*CM(7,8)+AMI(5,8)*CM(8,8)+AMI(5,9)*CM(9,8)
      ACM(5,9) = AMI(5,1)*CM(1,9)+AMI(5,2)*CM(2,9)+AMI(5,3)*CM(3,9)+AMI(
     #5,4)*CM(4,9)+AMI(5,5)*CM(5,9)+AMI(5,6)*CM(6,9)+AMI(5,7)*CM(7,9)+AM
     #I(5,8)*CM(8,9)+AMI(5,9)*CM(9,9)
      ACM(5,10) = AMI(5,1)*CM(1,10)+AMI(5,2)*CM(2,10)+AMI(5,3)*CM(3,10)+
     #AMI(5,4)*CM(4,10)+AMI(5,7)*CM(7,10)+AMI(5,9)*CM(9,10)
      ACM(5,11) = AMI(5,3)*CM(3,11)+AMI(5,5)*CM(5,11)+AMI(5,6)*CM(6,11)+
     #AMI(5,7)*CM(7,11)+AMI(5,8)*CM(8,11)+AMI(5,9)*CM(9,11)
      ACM(5,12) = AMI(5,1)*CM(1,12)+AMI(5,2)*CM(2,12)+AMI(5,3)*CM(3,12)+
     #AMI(5,4)*CM(4,12)+AMI(5,5)*CM(5,12)+AMI(5,6)*CM(6,12)+AMI(5,7)*CM(
     #7,12)+AMI(5,8)*CM(8,12)+AMI(5,9)*CM(9,12)
      ACM(6,1) = AMI(6,1)*CM(1,1)+AMI(6,2)*CM(2,1)+AMI(6,3)*CM(3,1)+AMI(
     #6,4)*CM(4,1)+AMI(6,7)*CM(7,1)+AMI(6,9)*CM(9,1)
      ACM(6,2) = AMI(6,3)*CM(3,2)+AMI(6,5)*CM(5,2)+AMI(6,6)*CM(6,2)+AMI(
     #6,7)*CM(7,2)+AMI(6,8)*CM(8,2)+AMI(6,9)*CM(9,2)
      ACM(6,3) = AMI(6,1)*CM(1,3)+AMI(6,2)*CM(2,3)+AMI(6,3)*CM(3,3)+AMI(
     #6,4)*CM(4,3)+AMI(6,5)*CM(5,3)+AMI(6,6)*CM(6,3)+AMI(6,7)*CM(7,3)+AM
     #I(6,8)*CM(8,3)+AMI(6,9)*CM(9,3)
      ACM(6,4) = AMI(6,1)*CM(1,4)+AMI(6,2)*CM(2,4)+AMI(6,3)*CM(3,4)+AMI(
     #6,4)*CM(4,4)+AMI(6,7)*CM(7,4)+AMI(6,9)*CM(9,4)
      ACM(6,5) = AMI(6,3)*CM(3,5)+AMI(6,5)*CM(5,5)+AMI(6,6)*CM(6,5)+AMI(
     #6,7)*CM(7,5)+AMI(6,8)*CM(8,5)+AMI(6,9)*CM(9,5)
      ACM(6,6) = AMI(6,1)*CM(1,6)+AMI(6,2)*CM(2,6)+AMI(6,3)*CM(3,6)+AMI(
     #6,4)*CM(4,6)+AMI(6,5)*CM(5,6)+AMI(6,6)*CM(6,6)+AMI(6,7)*CM(7,6)+AM
     #I(6,8)*CM(8,6)+AMI(6,9)*CM(9,6)
      ACM(6,7) = AMI(6,1)*CM(1,7)+AMI(6,2)*CM(2,7)+AMI(6,3)*CM(3,7)+AMI(
     #6,4)*CM(4,7)+AMI(6,7)*CM(7,7)+AMI(6,9)*CM(9,7)
      ACM(6,8) = AMI(6,3)*CM(3,8)+AMI(6,5)*CM(5,8)+AMI(6,6)*CM(6,8)+AMI(
     #6,7)*CM(7,8)+AMI(6,8)*CM(8,8)+AMI(6,9)*CM(9,8)
      ACM(6,9) = AMI(6,1)*CM(1,9)+AMI(6,2)*CM(2,9)+AMI(6,3)*CM(3,9)+AMI(
     #6,4)*CM(4,9)+AMI(6,5)*CM(5,9)+AMI(6,6)*CM(6,9)+AMI(6,7)*CM(7,9)+AM
     #I(6,8)*CM(8,9)+AMI(6,9)*CM(9,9)
      ACM(6,10) = AMI(6,1)*CM(1,10)+AMI(6,2)*CM(2,10)+AMI(6,3)*CM(3,10)+
     #AMI(6,4)*CM(4,10)+AMI(6,7)*CM(7,10)+AMI(6,9)*CM(9,10)
      ACM(6,11) = AMI(6,3)*CM(3,11)+AMI(6,5)*CM(5,11)+AMI(6,6)*CM(6,11)+
     #AMI(6,7)*CM(7,11)+AMI(6,8)*CM(8,11)+AMI(6,9)*CM(9,11)
      ACM(6,12) = AMI(6,1)*CM(1,12)+AMI(6,2)*CM(2,12)+AMI(6,3)*CM(3,12)+
     #AMI(6,4)*CM(4,12)+AMI(6,5)*CM(5,12)+AMI(6,6)*CM(6,12)+AMI(6,7)*CM(
     #7,12)+AMI(6,8)*CM(8,12)+AMI(6,9)*CM(9,12)
      ACM(7,1) = AMI(7,1)*CM(1,1)+AMI(7,2)*CM(2,1)+AMI(7,3)*CM(3,1)+AMI(
     #7,4)*CM(4,1)+AMI(7,7)*CM(7,1)+AMI(7,9)*CM(9,1)
      ACM(7,2) = AMI(7,3)*CM(3,2)+AMI(7,5)*CM(5,2)+AMI(7,6)*CM(6,2)+AMI(
     #7,7)*CM(7,2)+AMI(7,8)*CM(8,2)+AMI(7,9)*CM(9,2)
      ACM(7,3) = AMI(7,1)*CM(1,3)+AMI(7,2)*CM(2,3)+AMI(7,3)*CM(3,3)+AMI(
     #7,4)*CM(4,3)+AMI(7,5)*CM(5,3)+AMI(7,6)*CM(6,3)+AMI(7,7)*CM(7,3)+AM
     #I(7,8)*CM(8,3)+AMI(7,9)*CM(9,3)
      ACM(7,4) = AMI(7,1)*CM(1,4)+AMI(7,2)*CM(2,4)+AMI(7,3)*CM(3,4)+AMI(
     #7,4)*CM(4,4)+AMI(7,7)*CM(7,4)+AMI(7,9)*CM(9,4)
      ACM(7,5) = AMI(7,3)*CM(3,5)+AMI(7,5)*CM(5,5)+AMI(7,6)*CM(6,5)+AMI(
     #7,7)*CM(7,5)+AMI(7,8)*CM(8,5)+AMI(7,9)*CM(9,5)
      ACM(7,6) = AMI(7,1)*CM(1,6)+AMI(7,2)*CM(2,6)+AMI(7,3)*CM(3,6)+AMI(
     #7,4)*CM(4,6)+AMI(7,5)*CM(5,6)+AMI(7,6)*CM(6,6)+AMI(7,7)*CM(7,6)+AM
     #I(7,8)*CM(8,6)+AMI(7,9)*CM(9,6)
      ACM(7,7) = AMI(7,1)*CM(1,7)+AMI(7,2)*CM(2,7)+AMI(7,3)*CM(3,7)+AMI(
     #7,4)*CM(4,7)+AMI(7,7)*CM(7,7)+AMI(7,9)*CM(9,7)
      ACM(7,8) = AMI(7,3)*CM(3,8)+AMI(7,5)*CM(5,8)+AMI(7,6)*CM(6,8)+AMI(
     #7,7)*CM(7,8)+AMI(7,8)*CM(8,8)+AMI(7,9)*CM(9,8)
      ACM(7,9) = AMI(7,1)*CM(1,9)+AMI(7,2)*CM(2,9)+AMI(7,3)*CM(3,9)+AMI(
     #7,4)*CM(4,9)+AMI(7,5)*CM(5,9)+AMI(7,6)*CM(6,9)+AMI(7,7)*CM(7,9)+AM
     #I(7,8)*CM(8,9)+AMI(7,9)*CM(9,9)
      ACM(7,10) = AMI(7,1)*CM(1,10)+AMI(7,2)*CM(2,10)+AMI(7,3)*CM(3,10)+
     #AMI(7,4)*CM(4,10)+AMI(7,7)*CM(7,10)+AMI(7,9)*CM(9,10)
      ACM(7,11) = AMI(7,3)*CM(3,11)+AMI(7,5)*CM(5,11)+AMI(7,6)*CM(6,11)+
     #AMI(7,7)*CM(7,11)+AMI(7,8)*CM(8,11)+AMI(7,9)*CM(9,11)
      ACM(7,12) = AMI(7,1)*CM(1,12)+AMI(7,2)*CM(2,12)+AMI(7,3)*CM(3,12)+
     #AMI(7,4)*CM(4,12)+AMI(7,5)*CM(5,12)+AMI(7,6)*CM(6,12)+AMI(7,7)*CM(
     #7,12)+AMI(7,8)*CM(8,12)+AMI(7,9)*CM(9,12)
      ACM(8,1) = AMI(8,1)*CM(1,1)+AMI(8,2)*CM(2,1)+AMI(8,3)*CM(3,1)+AMI(
     #8,4)*CM(4,1)+AMI(8,7)*CM(7,1)+AMI(8,9)*CM(9,1)
      ACM(8,2) = AMI(8,3)*CM(3,2)+AMI(8,5)*CM(5,2)+AMI(8,6)*CM(6,2)+AMI(
     #8,7)*CM(7,2)+AMI(8,8)*CM(8,2)+AMI(8,9)*CM(9,2)
      ACM(8,3) = AMI(8,1)*CM(1,3)+AMI(8,2)*CM(2,3)+AMI(8,3)*CM(3,3)+AMI(
     #8,4)*CM(4,3)+AMI(8,5)*CM(5,3)+AMI(8,6)*CM(6,3)+AMI(8,7)*CM(7,3)+AM
     #I(8,8)*CM(8,3)+AMI(8,9)*CM(9,3)
      ACM(8,4) = AMI(8,1)*CM(1,4)+AMI(8,2)*CM(2,4)+AMI(8,3)*CM(3,4)+AMI(
     #8,4)*CM(4,4)+AMI(8,7)*CM(7,4)+AMI(8,9)*CM(9,4)
      ACM(8,5) = AMI(8,3)*CM(3,5)+AMI(8,5)*CM(5,5)+AMI(8,6)*CM(6,5)+AMI(
     #8,7)*CM(7,5)+AMI(8,8)*CM(8,5)+AMI(8,9)*CM(9,5)
      ACM(8,6) = AMI(8,1)*CM(1,6)+AMI(8,2)*CM(2,6)+AMI(8,3)*CM(3,6)+AMI(
     #8,4)*CM(4,6)+AMI(8,5)*CM(5,6)+AMI(8,6)*CM(6,6)+AMI(8,7)*CM(7,6)+AM
     #I(8,8)*CM(8,6)+AMI(8,9)*CM(9,6)
      ACM(8,7) = AMI(8,1)*CM(1,7)+AMI(8,2)*CM(2,7)+AMI(8,3)*CM(3,7)+AMI(
     #8,4)*CM(4,7)+AMI(8,7)*CM(7,7)+AMI(8,9)*CM(9,7)
      ACM(8,8) = AMI(8,3)*CM(3,8)+AMI(8,5)*CM(5,8)+AMI(8,6)*CM(6,8)+AMI(
     #8,7)*CM(7,8)+AMI(8,8)*CM(8,8)+AMI(8,9)*CM(9,8)
      ACM(8,9) = AMI(8,1)*CM(1,9)+AMI(8,2)*CM(2,9)+AMI(8,3)*CM(3,9)+AMI(
     #8,4)*CM(4,9)+AMI(8,5)*CM(5,9)+AMI(8,6)*CM(6,9)+AMI(8,7)*CM(7,9)+AM
     #I(8,8)*CM(8,9)+AMI(8,9)*CM(9,9)
      ACM(8,10) = AMI(8,1)*CM(1,10)+AMI(8,2)*CM(2,10)+AMI(8,3)*CM(3,10)+
     #AMI(8,4)*CM(4,10)+AMI(8,7)*CM(7,10)+AMI(8,9)*CM(9,10)
      ACM(8,11) = AMI(8,3)*CM(3,11)+AMI(8,5)*CM(5,11)+AMI(8,6)*CM(6,11)+
     #AMI(8,7)*CM(7,11)+AMI(8,8)*CM(8,11)+AMI(8,9)*CM(9,11)
      ACM(8,12) = AMI(8,1)*CM(1,12)+AMI(8,2)*CM(2,12)+AMI(8,3)*CM(3,12)+
     #AMI(8,4)*CM(4,12)+AMI(8,5)*CM(5,12)+AMI(8,6)*CM(6,12)+AMI(8,7)*CM(
     #7,12)+AMI(8,8)*CM(8,12)+AMI(8,9)*CM(9,12)
      ACM(9,1) = AMI(9,1)*CM(1,1)+AMI(9,2)*CM(2,1)+AMI(9,3)*CM(3,1)+AMI(
     #9,4)*CM(4,1)+AMI(9,7)*CM(7,1)+AMI(9,9)*CM(9,1)
      ACM(9,2) = AMI(9,3)*CM(3,2)+AMI(9,5)*CM(5,2)+AMI(9,6)*CM(6,2)+AMI(
     #9,7)*CM(7,2)+AMI(9,8)*CM(8,2)+AMI(9,9)*CM(9,2)
      ACM(9,3) = AMI(9,1)*CM(1,3)+AMI(9,2)*CM(2,3)+AMI(9,3)*CM(3,3)+AMI(
     #9,4)*CM(4,3)+AMI(9,5)*CM(5,3)+AMI(9,6)*CM(6,3)+AMI(9,7)*CM(7,3)+AM
     #I(9,8)*CM(8,3)+AMI(9,9)*CM(9,3)
      ACM(9,4) = AMI(9,1)*CM(1,4)+AMI(9,2)*CM(2,4)+AMI(9,3)*CM(3,4)+AMI(
     #9,4)*CM(4,4)+AMI(9,7)*CM(7,4)+AMI(9,9)*CM(9,4)
      ACM(9,5) = AMI(9,3)*CM(3,5)+AMI(9,5)*CM(5,5)+AMI(9,6)*CM(6,5)+AMI(
     #9,7)*CM(7,5)+AMI(9,8)*CM(8,5)+AMI(9,9)*CM(9,5)
      ACM(9,6) = AMI(9,1)*CM(1,6)+AMI(9,2)*CM(2,6)+AMI(9,3)*CM(3,6)+AMI(
     #9,4)*CM(4,6)+AMI(9,5)*CM(5,6)+AMI(9,6)*CM(6,6)+AMI(9,7)*CM(7,6)+AM
     #I(9,8)*CM(8,6)+AMI(9,9)*CM(9,6)
      ACM(9,7) = AMI(9,1)*CM(1,7)+AMI(9,2)*CM(2,7)+AMI(9,3)*CM(3,7)+AMI(
     #9,4)*CM(4,7)+AMI(9,7)*CM(7,7)+AMI(9,9)*CM(9,7)
      ACM(9,8) = AMI(9,3)*CM(3,8)+AMI(9,5)*CM(5,8)+AMI(9,6)*CM(6,8)+AMI(
     #9,7)*CM(7,8)+AMI(9,8)*CM(8,8)+AMI(9,9)*CM(9,8)
      ACM(9,9) = AMI(9,1)*CM(1,9)+AMI(9,2)*CM(2,9)+AMI(9,3)*CM(3,9)+AMI(
     #9,4)*CM(4,9)+AMI(9,5)*CM(5,9)+AMI(9,6)*CM(6,9)+AMI(9,7)*CM(7,9)+AM
     #I(9,8)*CM(8,9)+AMI(9,9)*CM(9,9)
      ACM(9,10) = AMI(9,1)*CM(1,10)+AMI(9,2)*CM(2,10)+AMI(9,3)*CM(3,10)+
     #AMI(9,4)*CM(4,10)+AMI(9,7)*CM(7,10)+AMI(9,9)*CM(9,10)
      ACM(9,11) = AMI(9,3)*CM(3,11)+AMI(9,5)*CM(5,11)+AMI(9,6)*CM(6,11)+
     #AMI(9,7)*CM(7,11)+AMI(9,8)*CM(8,11)+AMI(9,9)*CM(9,11)
      ACM(9,12) = AMI(9,1)*CM(1,12)+AMI(9,2)*CM(2,12)+AMI(9,3)*CM(3,12)+
     #AMI(9,4)*CM(4,12)+AMI(9,5)*CM(5,12)+AMI(9,6)*CM(6,12)+AMI(9,7)*CM(
     #7,12)+AMI(9,8)*CM(8,12)+AMI(9,9)*CM(9,12)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE NMSPM(AB1,AM,DR,RS,PSP)

	IMPLICIT REAL*8 (A-H,O-Z)
C	------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Nm*S*Pm - PSP
C
C	INPUT VARIABLES
C	AB1(4,4)	= STRAIN PARAMETERS
C	AM(9,9)		= MEMBRANE STRAIN PARAMETERS
C	DR(36)		= MEMBERANE RIGIDITY
C	RS(2,4)		= ELEMENT NODAL COORDINATES IN LOCAL SPACE
C
C	OUTPUT VARIAVBLES
C	PSP(5,5)	= Nm*S*Pm
C	------------------------------------------------------
      DIMENSION AB1(4,4),AM(9,9),DR(36),RS(2,4)
	DIMENSION PSP(9,9)

C	SOLVE FOR INT(R**2*S**2,drds)
      s1 = ((RS(2,2)/30-RS(2,4)/30)*RS(2,1)**2+(RS(2,2)**2/60-RS(2,4)**2
     #/60)*RS(2,1)+RS(2,2)**3/180-RS(2,4)**3/180)*RS(1,1)**3+((-RS(2,1)*
     #*3/30+RS(2,1)*RS(2,2)**2/60+RS(2,2)**3/60)*RS(1,2)+(RS(2,1)**3/30-
     #RS(2,1)*RS(2,4)**2/60-RS(2,4)**3/60)*RS(1,4))*RS(1,1)**2+((-RS(2,1
     #)**3/60-RS(2,1)**2*RS(2,2)/60+RS(2,2)**3/30)*RS(1,2)**2+(RS(2,1)**
     #3/60+RS(2,1)**2*RS(2,4)/60-RS(2,4)**3/30)*RS(1,4)**2)*RS(1,1)+(-RS
     #(2,1)**3/180-RS(2,1)**2*RS(2,2)/60-RS(2,1)*RS(2,2)**2/30+RS(2,2)**
     #2*RS(2,3)/30+RS(2,2)*RS(2,3)**2/60+RS(2,3)**3/180)*RS(1,2)**3+(-RS
     #(2,2)**3/30+RS(2,2)*RS(2,3)**2/60+RS(2,3)**3/60)*RS(1,3)*RS(1,2)**
     #2
      FACT = s1+(-RS(2,2)**3/60-RS(2,2)**2*RS(2,3)/60+RS(2,3)**3/30)
     #*RS(1,3)**2*RS(1,2)+(-RS(2,2)**3/180-RS(2,2)**2*RS(2,3)/60-RS(2,2)
     #*RS(2,3)**2/30+RS(2,3)**2*RS(2,4)/30+RS(2,3)*RS(2,4)**2/60+RS(2,4)
     #**3/180)*RS(1,3)**3+(-RS(2,3)**3/30+RS(2,3)*RS(2,4)**2/60+RS(2,4)*
     #*3/60)*RS(1,4)*RS(1,3)**2+(-RS(2,3)**3/60-RS(2,3)**2*RS(2,4)/60+RS
     #(2,4)**3/30)*RS(1,4)**2*RS(1,3)+(RS(2,1)**3/180+RS(2,1)**2*RS(2,4)
     #/60+RS(2,1)*RS(2,4)**2/30-RS(2,3)**3/180-RS(2,3)**2*RS(2,4)/60-RS(
     #2,3)*RS(2,4)**2/30)*RS(1,4)**3

      PSP(1,1) = AB1(1,1)*DR(1)
      PSP(1,2) = AB1(1,3)*DR(1)
      PSP(1,3) = AB1(1,2)*DR(1)
      PSP(1,4) = AB1(3,3)*DR(1)
      PSP(1,5) = AB1(1,1)*DR(2)
      PSP(1,6) = AB1(1,2)*DR(2)
      PSP(1,7) = AB1(1,3)*DR(2)
      PSP(1,8) = AB1(2,2)*DR(2)

      PSP(2,1) = PSP(1,2)
      PSP(2,2) = AB1(3,3)*DR(1)
      PSP(2,3) = AB1(1,4)*DR(1)
      PSP(2,4) = AM(2,4)*DR(1)
      PSP(2,5) = AB1(1,3)*DR(2)
      PSP(2,6) = AB1(1,4)*DR(2)
      PSP(2,7) = AB1(3,3)*DR(2)
      PSP(2,8) = AM(7,8)*DR(2)

      PSP(3,1) = PSP(1,3)
      PSP(3,2) = PSP(2,3)
      PSP(3,3) = AB1(2,2)*DR(1)+AB1(3,3)*DR(15)
      PSP(3,4) = AM(3,4)*DR(1)
      PSP(3,5) = AB1(1,2)*DR(2)
      PSP(3,6) = AB1(2,2)*DR(2)
      PSP(3,7) = AB1(1,4)*(DR(2)+DR(15))
      PSP(3,8) = AM(6,8)*DR(2)
      PSP(3,9) = -AB1(1,3)*DR(15)

      PSP(4,1) = PSP(1,4)
      PSP(4,2) = PSP(2,4)
      PSP(4,3) = PSP(3,4)
      PSP(4,4) = AM(4,4)*DR(1)
      PSP(4,5) = AB1(3,3)*DR(2)
      PSP(4,6) = AM(3,4)*DR(2)
      PSP(4,7) = AM(2,4)*DR(2)
      PSP(4,8) = FACT*DR(2)

      PSP(5,1) = PSP(1,5)
      PSP(5,2) = PSP(2,5)
      PSP(5,3) = PSP(3,5)
      PSP(5,4) = PSP(4,5)
      PSP(5,5) = AB1(1,1)*DR(8)
      PSP(5,6) = AB1(1,2)*DR(8)
      PSP(5,7) = AB1(1,3)*DR(8)
      PSP(5,8) = AB1(2,2)*DR(8)

      PSP(6,1) = PSP(1,6)
      PSP(6,2) = PSP(2,6)
      PSP(6,3) = PSP(3,6)
      PSP(6,4) = PSP(4,6)
      PSP(6,5) = PSP(5,6)
      PSP(6,6) = AB1(2,2)*DR(8)
      PSP(6,7) = AB1(1,4)*DR(8)
      PSP(6,8) = AM(6,8)*DR(8)

      PSP(7,1) = PSP(1,7)
      PSP(7,2) = PSP(2,7)
      PSP(7,3) = PSP(3,7)
      PSP(7,4) = PSP(4,7)
      PSP(7,5) = PSP(5,7)
      PSP(7,6) = PSP(6,7)
      PSP(7,7) = AB1(3,3)*DR(8)+AB1(2,2)*DR(15)
      PSP(7,8) = AM(7,8)*DR(8)
      PSP(7,9) = -AB1(1,2)*DR(15)

      PSP(8,1) = PSP(1,8)
      PSP(8,2) = PSP(2,8)
      PSP(8,3) = PSP(3,8)
      PSP(8,4) = PSP(4,8)
      PSP(8,5) = PSP(5,8)
      PSP(8,6) = PSP(6,8)
      PSP(8,7) = PSP(7,8)
      PSP(8,8) = AM(8,8)*DR(8)

      PSP(9,3) = PSP(3,9)
      PSP(9,7) = PSP(7,9)
      PSP(9,9) = AB1(1,1)*DR(15)

	RETURN
	END
C
C=====================================================================
	SUBROUTINE BCGEO(AM,CM,ACG)
	IMPLICIT REAL*8 (A-H,O-Z)
C	----------------------------------------------------
C	PURPOSE:	TO SOLVE FOR INV(Ag)*(Cg) - ACG
C
C	INPUT VARIABLES
C	AM(9,9)				= STRAIN PARAMETER
C	CG(4,12)			= NONLINEAR STRAIN PARAMETER
C
C	LOCAL VARIABLE
C	AGI(4,4)	= INVERSE OF A SUB-MATRIX OF Ag
C
C	OUTPUT VARIAVBLE
C	ACG(4,12)	= INVERSE(Ag)*(Cg)
C	----------------------------------------------------
      DIMENSION AM(9,9),CM(9,12)
	DIMENSION AGI(4,4)
	DIMENSION ACG(4,12)

C	------------------------------------------------
C	SOLVE FOR THE INVERSE(Ag)*(Cg) - ACG
C	------------------------------------------------
C	--------------------
C	COMPONENTS FOR DU/DR
C	--------------------
	ACG(1,1) = CM(1,1)/AM(1,1)
	ACG(1,4) = CM(1,4)/AM(1,1)
	ACG(1,7) = CM(1,7)/AM(1,1)
	ACG(1,10)= CM(1,10)/AM(1,1)

	ACG(1,3) = CM(1,3)/AM(1,1)
	ACG(1,6) = CM(1,6)/AM(1,1)
	ACG(1,9) = CM(1,9)/AM(1,1)
	ACG(1,12)= CM(1,12)/AM(1,1)

C	--------------------
C	COMPONENTS FOR DU/DS
C	--------------------
	ACG(2,1) = CM(5,2)/AM(1,1)
	ACG(2,4) = CM(5,5)/AM(1,1)
	ACG(2,7) = CM(5,8)/AM(1,1)
	ACG(2,10)= CM(5,11)/AM(1,1)

	ACG(2,3) = CM(5,3)/AM(1,1)
	ACG(2,6) = CM(5,6)/AM(1,1)
	ACG(2,9) = CM(5,9)/AM(1,1)
	ACG(2,12)= CM(5,12)/AM(1,1)

C	--------------------
C	COMPONENTS FOR DV/DR
C	--------------------
	ACG(3,2) = CM(1,1)/AM(1,1)
	ACG(3,5) = CM(1,4)/AM(1,1)
	ACG(3,8) = CM(1,7)/AM(1,1)
	ACG(3,11)= CM(1,10)/AM(1,1)

	ACG(3,3) = CM(1,3)/AM(1,1)
	ACG(3,6) = CM(1,6)/AM(1,1)
	ACG(3,9) = CM(1,9)/AM(1,1)
	ACG(3,12)= CM(1,12)/AM(1,1)

C	--------------------
C	COMPONENTS FOR DV/DS
C	--------------------
	ACG(4,2) = CM(5,2)/AM(1,1)
	ACG(4,5) = CM(5,5)/AM(1,1)
	ACG(4,8) = CM(5,8)/AM(1,1)
	ACG(4,11)= CM(5,11)/AM(1,1)

	ACG(4,3) = CM(5,3)/AM(1,1)
	ACG(4,6) = CM(5,6)/AM(1,1)
	ACG(4,9) = CM(5,9)/AM(1,1)
	ACG(4,12)= CM(5,12)/AM(1,1)

	RETURN
	END
C
C=====================================================================
	SUBROUTINE NGFPG(AM,APM,DR,PFP)
	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------
C	PURPOSE:	TO SOLVE FOR Ng*F*Pg - PFP
C
C	INPUT VARIABLES
C	AM(9,9)	= SUB-MATRIX OF INTEGRAL{TRANS(Pb)*Pb}
C	DR(36)		= MODULI MATRIX
C
C
C	OUTPUT VARIAVBLE
C	PFP(51,51)	= Ng*F*Pg
C	--------------------------------------------------

      DIMENSION AM(9,9),APM(9),DR(36)
	DIMENSION PFP(4,4)

C	-----------------------
C	SOLVE FOR Ng*F*Pg - PFP
C	-----------------------
      t3 = AM(1,4)
      t5 = AM(1,1)*APM(1)+AM(1,2)*APM(2)+AM(1,3)*APM(3)+t3*APM(4)
      t9 = AM(5,8)
      t11 = AM(1,1)*APM(5)+AM(1,3)*APM(6)+AM(1,2)*APM(7)+t9*APM(8)
      t13 = t5*DR(1)+t11*DR(2)
      t17 = (-AM(1,2)*APM(3)-AM(1,3)*APM(7)+AM(1,1)*APM(9))*DR(15)
      t38 = t11*DR(1)+t5*DR(2)

      PFP(1,1) = t13
      PFP(1,2) = t17
	PFP(1,3) = 0.0
	PFP(1,4) = 0.0

      PFP(2,1) = t17
      PFP(2,2) = t38
      PFP(2,3) = 0.0
	PFP(2,4) = 0.0

	PFP(3,1) = 0.0
      PFP(3,2) = 0.0
      PFP(3,3) = t13
      PFP(3,4) = t17

      PFP(4,1) = 0.0
      PFP(4,2) = 0.0
      PFP(4,3) = t17
      PFP(4,4) = t38

	RETURN
	END
C
C=====================================================================
      SUBROUTINE NFPPLAS(PAPM,PFP)
	IMPLICIT REAL*8 (A-H,O-Z)
C=============================================================
C	PURPOSE: TO COMPUTE FOR THE INTEGRAL(TRANSPOSE(Pg)*F*Pg)
C
C	INPUT VARIABLES
C		PAPM(9)   = INTEGRAL(TRANS(Pm)*N)
C
C	LOCAL VARIABLES
C
C	OUTPUT VARIABLE
C		PFP(4,4) = INTEGRAL(TRANSPOSE(Pg)*F*Pg)
C=============================================================
	DIMENSION PAPM(9),PFP(4,4)

C	-------------------------------------
C	ASSEMBLE INTEGRAL(TRANSPOSE(Pg)*F*Pg)
C	-------------------------------------
      PFP(1,1) = PAPM(1)
      PFP(1,2) = PAPM(5)
      PFP(1,3) = 0.0
      PFP(1,4) = 0.0

      PFP(2,1) = PAPM(5)
      PFP(2,2) = PAPM(9)
      PFP(2,3) = 0.0
      PFP(2,4) = 0.0

      PFP(3,1) = 0.0
      PFP(3,2) = 0.0
      PFP(3,3) = PAPM(1)
      PFP(3,4) = PAPM(5)

      PFP(4,1) = 0.0
      PFP(4,2) = 0.0
      PFP(4,3) = PAPM(5)
      PFP(4,4) = PAPM(9)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE KGEO(ACG,PFP,ST)
	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR THE GEOMETRIC STIFFNESS - S
C
C	INPUT VARIABLES
C	ACG(4,12)	= INVERSE(Ag)*(Cg)*TRANSPOSE(T)
C	PFP(4,4)	= Ng*F*Pg
C
C	LOACAL VARIABLES
C	AKG(12,12)		= GEOMETRIC STIFFNESS IN FULL MATRIX FORM
C	TEMP1(4,12)		= PFP*ACG
C	TEMP2(12,4)		= TRANSPOSE(ACG)
C
C	OUTPUT VARIAVBLES
C	S(*)		= GEOMETRIC STIFFNESS - UPPER TRIANG - VECTOR FORM
C	--------------------------------------------------------------
      DIMENSION ACG(4,12),PFP(4,4)
	DIMENSION ST(78)

      t3 = ACG(1,1)*PFP(1,1)+ACG(2,1)*PFP(1,2)
      t7 = ACG(1,1)*PFP(1,2)+ACG(2,1)*PFP(2,2)
      t33 = ACG(3,2)*PFP(1,1)+ACG(4,2)*PFP(1,2)
      t37 = ACG(3,2)*PFP(1,2)+ACG(4,2)*PFP(2,2)
      t63 = ACG(1,3)*PFP(1,1)+ACG(2,3)*PFP(1,2)
      t67 = ACG(1,3)*PFP(1,2)+ACG(2,3)*PFP(2,2)
      t71 = ACG(3,3)*PFP(1,1)+ACG(4,3)*PFP(1,2)
      t75 = ACG(3,3)*PFP(1,2)+ACG(4,3)*PFP(2,2)
      t113 = ACG(1,4)*PFP(1,1)+ACG(2,4)*PFP(1,2)
      t117 = ACG(1,4)*PFP(1,2)+ACG(2,4)*PFP(2,2)
      t137 = ACG(3,5)*PFP(1,1)+ACG(4,5)*PFP(1,2)
      t141 = ACG(3,5)*PFP(1,2)+ACG(4,5)*PFP(2,2)
      t161 = ACG(1,6)*PFP(1,1)+ACG(2,6)*PFP(1,2)
      t165 = ACG(1,6)*PFP(1,2)+ACG(2,6)*PFP(2,2)
      t169 = ACG(3,6)*PFP(1,1)+ACG(4,6)*PFP(1,2)
      t173 = ACG(3,6)*PFP(1,2)+ACG(4,6)*PFP(2,2)
      t200 = ACG(1,7)*PFP(1,1)+ACG(2,7)*PFP(1,2)
      t204 = ACG(1,7)*PFP(1,2)+ACG(2,7)*PFP(2,2)
      t218 = ACG(3,8)*PFP(1,1)+ACG(4,8)*PFP(1,2)
      t222 = ACG(3,8)*PFP(1,2)+ACG(4,8)*PFP(2,2)
      t236 = ACG(1,9)*PFP(1,1)+ACG(2,9)*PFP(1,2)
      t240 = ACG(1,9)*PFP(1,2)+ACG(2,9)*PFP(2,2)
      t244 = ACG(3,9)*PFP(1,1)+ACG(4,9)*PFP(1,2)
      t248 = ACG(3,9)*PFP(1,2)+ACG(4,9)*PFP(2,2)
      t264 = ACG(1,10)*PFP(1,1)+ACG(2,10)*PFP(1,2)
      t268 = ACG(1,10)*PFP(1,2)+ACG(2,10)*PFP(2,2)
      t276 = ACG(3,11)*PFP(1,1)+ACG(4,11)*PFP(1,2)
      t280 = ACG(3,11)*PFP(1,2)+ACG(4,11)*PFP(2,2)

      ST(1) = ST(1)+t3*ACG(1,1)+t7*ACG(2,1)
C      ST(2) = ST(2)
      ST(3) = ST(3)+t3*ACG(1,3)+t7*ACG(2,3)
      ST(4) = ST(4)+t3*ACG(1,4)+t7*ACG(2,4)
C      ST(5) = ST(5)
      ST(6) = ST(6)+t3*ACG(1,6)+t7*ACG(2,6)
      ST(7) = ST(7)+t3*ACG(1,7)+t7*ACG(2,7)
C      ST(8) = ST(8)
      ST(9) = ST(9)+t3*ACG(1,9)+t7*ACG(2,9)
      ST(10) = ST(10)+t3*ACG(1,10)+t7*ACG(2,10)
C      ST(11) = ST(11)
      ST(12) = ST(12)+t3*ACG(1,12)+t7*ACG(2,12)
      ST(13) = ST(13)+t33*ACG(3,2)+t37*ACG(4,2)
      ST(14) = ST(14)+t33*ACG(3,3)+t37*ACG(4,3)
C      ST(15) = ST(15)
      ST(16) = ST(16)+t33*ACG(3,5)+t37*ACG(4,5)
      ST(17) = ST(17)+t33*ACG(3,6)+t37*ACG(4,6)
C      ST(18) = ST(18)
      ST(19) = ST(19)+t33*ACG(3,8)+t37*ACG(4,8)
      ST(20) = ST(20)+t33*ACG(3,9)+t37*ACG(4,9)
C      ST(21) = ST(21)
      ST(22) = ST(22)+t33*ACG(3,11)+t37*ACG(4,11)
      ST(23) = ST(23)+t33*ACG(3,12)+t37*ACG(4,12)
      ST(24) =ST(24)+t63*ACG(1,3)+t67*ACG(2,3)+t71*ACG(3,3)+t75*ACG(4,3)
      ST(25) = ST(25)+t63*ACG(1,4)+t67*ACG(2,4)
      ST(26) = ST(26)+t71*ACG(3,5)+t75*ACG(4,5)
      ST(27) =ST(27)+t63*ACG(1,6)+t67*ACG(2,6)+t71*ACG(3,6)+t75*ACG(4,6)
      ST(28) = ST(28)+t63*ACG(1,7)+t67*ACG(2,7)
      ST(29) = ST(29)+t71*ACG(3,8)+t75*ACG(4,8)
      ST(30) =ST(30)+t63*ACG(1,9)+t67*ACG(2,9)+t71*ACG(3,9)+t75*ACG(4,9)
      ST(31) = ST(31)+t63*ACG(1,10)+t67*ACG(2,10)
      ST(32) = ST(32)+t71*ACG(3,11)+t75*ACG(4,11)
      ST(33) =ST(33)+t63*ACG(1,12)+t67*ACG(2,12)+t71*ACG(3,12)+t75*ACG(4
     #,12)
      ST(34) = ST(34)+t113*ACG(1,4)+t117*ACG(2,4)
C      ST(35) = ST(35)
      ST(36) = ST(36)+t113*ACG(1,6)+t117*ACG(2,6)
      ST(37) = ST(37)+t113*ACG(1,7)+t117*ACG(2,7)
C      ST(38) = ST(38)
      ST(39) = ST(39)+t113*ACG(1,9)+t117*ACG(2,9)
      ST(40) = ST(40)+t113*ACG(1,10)+t117*ACG(2,10)
C      ST(41) = ST(41)
      ST(42) = ST(42)+t113*ACG(1,12)+t117*ACG(2,12)
      ST(43) = ST(43)+t137*ACG(3,5)+t141*ACG(4,5)
      ST(44) = ST(44)+t137*ACG(3,6)+t141*ACG(4,6)
C      ST(45) = ST(45)
      ST(46) = ST(46)+t137*ACG(3,8)+t141*ACG(4,8)
      ST(47) = ST(47)+t137*ACG(3,9)+t141*ACG(4,9)
C      ST(48) = ST(48)
      ST(49) = ST(49)+t137*ACG(3,11)+t141*ACG(4,11)
      ST(50) = ST(50)+t137*ACG(3,12)+t141*ACG(4,12)
      ST(51) =ST(51)+t161*ACG(1,6)+t165*ACG(2,6)+t169*ACG(3,6)+t173*ACG(
     #4,6)
      ST(52) = ST(52)+t161*ACG(1,7)+t165*ACG(2,7)
      ST(53) = ST(53)+t169*ACG(3,8)+t173*ACG(4,8)
      ST(54) =ST(54)+t161*ACG(1,9)+t165*ACG(2,9)+t169*ACG(3,9)+t173*ACG(
     #4,9)
      ST(55) = ST(55)+t161*ACG(1,10)+t165*ACG(2,10)
      ST(56) = ST(56)+t169*ACG(3,11)+t173*ACG(4,11)
      ST(57) =ST(57)+t161*ACG(1,12)+t165*ACG(2,12)+t169*ACG(3,12)+t173*A
     #CG(4,12)
      ST(58) = ST(58)+t200*ACG(1,7)+t204*ACG(2,7)
C      ST(59) = ST(59)
      ST(60) = ST(60)+t200*ACG(1,9)+t204*ACG(2,9)
      ST(61) = ST(61)+t200*ACG(1,10)+t204*ACG(2,10)
C      ST(62) = ST(62)
      ST(63) = ST(63)+t200*ACG(1,12)+t204*ACG(2,12)
      ST(64) = ST(64)+t218*ACG(3,8)+t222*ACG(4,8)
      ST(65) = ST(65)+t218*ACG(3,9)+t222*ACG(4,9)
C      ST(66) = ST(66)
      ST(67) = ST(67)+t218*ACG(3,11)+t222*ACG(4,11)
      ST(68) = ST(68)+t218*ACG(3,12)+t222*ACG(4,12)
      ST(69) =ST(69)+t236*ACG(1,9)+t240*ACG(2,9)+t244*ACG(3,9)+t248*ACG(
     #4,9)
      ST(70) = ST(70)+t236*ACG(1,10)+t240*ACG(2,10)
      ST(71) = ST(71)+t244*ACG(3,11)+t248*ACG(4,11)
      ST(72) =ST(72)+t236*ACG(1,12)+t240*ACG(2,12)+t244*ACG(3,12)+t248*A
     #CG(4,12)
      ST(73) = ST(73)+t264*ACG(1,10)+t268*ACG(2,10)
C      ST(74) = ST(74)
      ST(75) = ST(75)+t264*ACG(1,12)+t268*ACG(2,12)
      ST(76) = ST(76)+t276*ACG(3,11)+t280*ACG(4,11)
      ST(77) = ST(77)+t276*ACG(3,12)+t280*ACG(4,12)
      ST(78) =ST(78)+(ACG(1,12)*PFP(1,1)+ACG(2,12)*PFP(1,2))*ACG(1,12)+(
     #ACG(1,12)*PFP(1,2)+ACG(2,12)*PFP(2,2))*ACG(2,12)+(ACG(3,12)*PFP(1,
     #1)+ACG(4,12)*PFP(1,2))*ACG(3,12)+(ACG(3,12)*PFP(1,2)+ACG(4,12)*PFP
     #(2,2))*ACG(4,12)

	RETURN
	END
C
C=====================================================================
      SUBROUTINE MASSQC (CM,H,DE,DVOL,NNO,NEF,IPT)
      IMPLICIT REAL*8(A-H,O-Z)
C	------------------------------------------------------------------
C	PROGRAMED BY GILSON - DEC2002
C	PURPOSE: CALCULATE THE MASS MATRIX OF A QUADRILATERAL ELEMENT
C
C	DESCRIPTION OF VARIABLES
C		CM		= ELEMENT MASS MATRIX
C		DE		= DENSITY
C		DVOL	= VOLUME FACTOR
C		FAC		= MASS FACTOR
C		H,D		= SHAPE FUNCTION
C		IPT		= INTEGRATION POINT NUMBER
C		NNO		= NUMBER OF ELEMENT NODES
C		IMASS	= FLAG FOR TYPE OF MASS MATRIX
C					1,2 - CONSISTENTLY LUMPED MASS MATRIX
C					3   - CONSISTENT MASS MATRIX
C		ISTYP	= TYPE OF MEMBRANE ELEMENT
C					0 - AXISYMMETRIC
C					1 - PLANE STRAIN
C					2 - PLANE STRESS
C		KL,KS	= POINTERS
C		NEF		= NUMBER OF DEGREES OF FREEDOM (2 * NUMBER OF NODES)
C	------------------------------------------------------------------
      COMMON /DYNA/ CDEN,IMASS
C
      DIMENSION CM(1),D(12),H(8)
C
	FAC=DVOL*DE
C
      IF (IMASS.LT.3) GO TO 320
C
C	---------------
C     CONSISTENT MASS
C	---------------
      DO 200 I = 1,NNO
        D(3*I-2) = H(I)
        D(3*I-1) = H(I)
        D(3*I  ) = H(I)
200   CONTINUE

      KL = 1
      DO 300 I = 1,NEF,3
	  DO 301 J = I,NEF,3
		CM(KL) = CM(KL) + D(I)*D(J)*FAC
301     KL = KL+3
300   KL = KL+(NEF-I)*2-1
C
      KL = 1
      DO 401 I = 1,NEF,3
        KS = KL+NEF-I+1
        DO 400 J = I,NEF,3
		CM(KS) = CM(KL)
		KS = KS+3
400	  KL = KL+3
401	KL = KL+(NEF-I)*2-1
C
      RETURN
C
C	-----------
C     LUMPED MASS
C	-----------
320	CONTINUE

	DO 500 I = 1,NEF,3
	  FACM = FAC/NNO
	  CM(I) = CM(I)+FACM
	  CM(I+1) = CM(I)
500	CONTINUE

	RETURN
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE VRS(COORD,VR,VS,RS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
C     -------------------------------------------------------
C     PROGRAM TO DETERMINE LOCAL CORDINATE VECTORS - VR,VS,VT
C     -------------------------------------------------------
      DIMENSION COORD(2,8)
	DIMENSION XB(2,8),EC(2),VL1(2),VL2(2),DUMMY(2)
	DIMENSION VR(2),VS(2),RS(2,4)

C	--------------------------------------
C	COMPUTE FOR MID-POINT COORDINATES - XB
C	--------------------------------------
	XB(1,1)=(COORD(1,1)+COORD(1,4))/2.0
	XB(2,1)=(COORD(2,1)+COORD(2,4))/2.0
	DO 10 I=2,NNO
	 XB(1,I)=(COORD(1,I)+COORD(1,I-1))/2.0
	 XB(2,I)=(COORD(2,I)+COORD(2,I-1))/2.0
10	CONTINUE
C	-------------------------------------------
C	COMPUTE FOR ELEMENT CENTER COORDINATES - EC
C	-------------------------------------------
	EC(1)=(XB(1,3)+XB(1,1))/2.0
	EC(2)=(XB(2,3)+XB(2,1))/2.0
C	----------------------------------------------
C	COMPUTE FOR OPPOSITE MID-POINT VECTORS - L1,L2
C	----------------------------------------------
	VL1(1)=XB(1,1)-XB(1,3)
	VL1(2)=XB(2,1)-XB(2,3)
	CALL SCALEN(VL1,VL1,DUM,2)

	VL2(1)=XB(1,2)-XB(1,4)
	VL2(2)=XB(2,2)-XB(2,4)
	CALL SCALEN(VL2,VL2,DUM,2)
C	-----------------------------------------------
C	COMPUTE FOR LOCAL COORDINATE VECTORS - VR,VS,VT
C	-----------------------------------------------
	VSMAG = SQRT((-VL1(2)+VL2(1))**2+(VL1(1)+VL2(2))**2)
	VS(1) = (-VL1(2)+VL2(1))/VSMAG
	VS(2) = (VL1(1)+VL2(2))/VSMAG

	VR(1) = VS(2)
	VR(2) = -VS(1)
C	-----------------------------------------------
C	DETERMINE NODAL COORDINATES IN LOCAL SPACE - RS
C	-----------------------------------------------
	DO 20 I=1,4
	 DO 30 J=1,2
	  DUMMY(J)=COORD(J,I)-EC(J)
30	 CONTINUE
	 RS(1,I)=DUMMY(1)*VR(1)+DUMMY(2)*VR(2)
	 RS(2,I)=DUMMY(1)*VS(1)+DUMMY(2)*VS(2)
20	CONTINUE

      RETURN
      END
C
C=====================================================================
      SUBROUTINE KMEM(ACM,PSP,ST)
	IMPLICIT REAL*8 (A-H,O-Z)
C	--------------------------------------------------------------
C	PURPOSE:	TO SOLVE FOR THE LINEAR STIFFNESS - ST
C
C	INPUT VARIABLES
C	ACM(9,12)	= INVERSE(Am)*(Cm)
C	PFP(9,9)	= Nm*Dm*Pm
C
C
C	OUTPUT VARIAVBLES
C	ST(78)		= LINEAR STIFFNESS - UPPER TRIANG - VECTOR FORM
C	--------------------------------------------------------------
      DIMENSION ACM(9,12),PSP(9,9)
	DIMENSION ST(78)

      t10 = ACM(1,1)*PSP(1,1)+ACM(2,1)*PSP(1,2)+ACM(3,1)*PSP(1,3)+ACM(4,
     #1)*PSP(1,4)+ACM(5,1)*PSP(1,5)+ACM(6,1)*PSP(1,6)+ACM(7,1)*PSP(1,7)+
     #ACM(8,1)*PSP(1,8)+ACM(9,1)*PSP(1,9)
      t21 = ACM(1,1)*PSP(1,2)+ACM(2,1)*PSP(2,2)+ACM(3,1)*PSP(2,3)+ACM(4,
     #1)*PSP(2,4)+ACM(5,1)*PSP(2,5)+ACM(6,1)*PSP(2,6)+ACM(7,1)*PSP(2,7)+
     #ACM(8,1)*PSP(2,8)+ACM(9,1)*PSP(2,9)
      t32 = ACM(1,1)*PSP(1,3)+ACM(2,1)*PSP(2,3)+ACM(3,1)*PSP(3,3)+ACM(4,
     #1)*PSP(3,4)+ACM(5,1)*PSP(3,5)+ACM(6,1)*PSP(3,6)+ACM(7,1)*PSP(3,7)+
     #ACM(8,1)*PSP(3,8)+ACM(9,1)*PSP(3,9)
      t43 = ACM(1,1)*PSP(1,4)+ACM(2,1)*PSP(2,4)+ACM(3,1)*PSP(3,4)+ACM(4,
     #1)*PSP(4,4)+ACM(5,1)*PSP(4,5)+ACM(6,1)*PSP(4,6)+ACM(7,1)*PSP(4,7)+
     #ACM(8,1)*PSP(4,8)+ACM(9,1)*PSP(4,9)
      t54 = ACM(1,1)*PSP(1,5)+ACM(2,1)*PSP(2,5)+ACM(3,1)*PSP(3,5)+ACM(4,
     #1)*PSP(4,5)+ACM(5,1)*PSP(5,5)+ACM(6,1)*PSP(5,6)+ACM(7,1)*PSP(5,7)+
     #ACM(8,1)*PSP(5,8)+ACM(9,1)*PSP(5,9)
      t65 = ACM(1,1)*PSP(1,6)+ACM(2,1)*PSP(2,6)+ACM(3,1)*PSP(3,6)+ACM(4,
     #1)*PSP(4,6)+ACM(5,1)*PSP(5,6)+ACM(6,1)*PSP(6,6)+ACM(7,1)*PSP(6,7)+
     #ACM(8,1)*PSP(6,8)+ACM(9,1)*PSP(6,9)
      t76 = ACM(1,1)*PSP(1,7)+ACM(2,1)*PSP(2,7)+ACM(3,1)*PSP(3,7)+ACM(4,
     #1)*PSP(4,7)+ACM(5,1)*PSP(5,7)+ACM(6,1)*PSP(6,7)+ACM(7,1)*PSP(7,7)+
     #ACM(8,1)*PSP(7,8)+ACM(9,1)*PSP(7,9)
      t87 = ACM(1,1)*PSP(1,8)+ACM(2,1)*PSP(2,8)+ACM(3,1)*PSP(3,8)+ACM(4,
     #1)*PSP(4,8)+ACM(5,1)*PSP(5,8)+ACM(6,1)*PSP(6,8)+ACM(7,1)*PSP(7,8)+
     #ACM(8,1)*PSP(8,8)+ACM(9,1)*PSP(8,9)
      t98 = ACM(1,1)*PSP(1,9)+ACM(2,1)*PSP(2,9)+ACM(3,1)*PSP(3,9)+ACM(4,
     #1)*PSP(4,9)+ACM(5,1)*PSP(5,9)+ACM(6,1)*PSP(6,9)+ACM(7,1)*PSP(7,9)+
     #ACM(8,1)*PSP(8,9)+ACM(9,1)*PSP(9,9)
      t220 = ACM(1,2)*PSP(1,1)+ACM(2,2)*PSP(1,2)+ACM(3,2)*PSP(1,3)+ACM(4
     #,2)*PSP(1,4)+ACM(5,2)*PSP(1,5)+ACM(6,2)*PSP(1,6)+ACM(7,2)*PSP(1,7)
     #+ACM(8,2)*PSP(1,8)+ACM(9,2)*PSP(1,9)
      t231 = ACM(1,2)*PSP(1,2)+ACM(2,2)*PSP(2,2)+ACM(3,2)*PSP(2,3)+ACM(4
     #,2)*PSP(2,4)+ACM(5,2)*PSP(2,5)+ACM(6,2)*PSP(2,6)+ACM(7,2)*PSP(2,7)
     #+ACM(8,2)*PSP(2,8)+ACM(9,2)*PSP(2,9)
      t242 = ACM(1,2)*PSP(1,3)+ACM(2,2)*PSP(2,3)+ACM(3,2)*PSP(3,3)+ACM(4
     #,2)*PSP(3,4)+ACM(5,2)*PSP(3,5)+ACM(6,2)*PSP(3,6)+ACM(7,2)*PSP(3,7)
     #+ACM(8,2)*PSP(3,8)+ACM(9,2)*PSP(3,9)
      t253 = ACM(1,2)*PSP(1,4)+ACM(2,2)*PSP(2,4)+ACM(3,2)*PSP(3,4)+ACM(4
     #,2)*PSP(4,4)+ACM(5,2)*PSP(4,5)+ACM(6,2)*PSP(4,6)+ACM(7,2)*PSP(4,7)
     #+ACM(8,2)*PSP(4,8)+ACM(9,2)*PSP(4,9)
      t264 = ACM(1,2)*PSP(1,5)+ACM(2,2)*PSP(2,5)+ACM(3,2)*PSP(3,5)+ACM(4
     #,2)*PSP(4,5)+ACM(5,2)*PSP(5,5)+ACM(6,2)*PSP(5,6)+ACM(7,2)*PSP(5,7)
     #+ACM(8,2)*PSP(5,8)+ACM(9,2)*PSP(5,9)
      t275 = ACM(1,2)*PSP(1,6)+ACM(2,2)*PSP(2,6)+ACM(3,2)*PSP(3,6)+ACM(4
     #,2)*PSP(4,6)+ACM(5,2)*PSP(5,6)+ACM(6,2)*PSP(6,6)+ACM(7,2)*PSP(6,7)
     #+ACM(8,2)*PSP(6,8)+ACM(9,2)*PSP(6,9)
      t286 = ACM(1,2)*PSP(1,7)+ACM(2,2)*PSP(2,7)+ACM(3,2)*PSP(3,7)+ACM(4
     #,2)*PSP(4,7)+ACM(5,2)*PSP(5,7)+ACM(6,2)*PSP(6,7)+ACM(7,2)*PSP(7,7)
     #+ACM(8,2)*PSP(7,8)+ACM(9,2)*PSP(7,9)
      t297 = ACM(1,2)*PSP(1,8)+ACM(2,2)*PSP(2,8)+ACM(3,2)*PSP(3,8)+ACM(4
     #,2)*PSP(4,8)+ACM(5,2)*PSP(5,8)+ACM(6,2)*PSP(6,8)+ACM(7,2)*PSP(7,8)
     #+ACM(8,2)*PSP(8,8)+ACM(9,2)*PSP(8,9)
      t308 = ACM(1,2)*PSP(1,9)+ACM(2,2)*PSP(2,9)+ACM(3,2)*PSP(3,9)+ACM(4
     #,2)*PSP(4,9)+ACM(5,2)*PSP(5,9)+ACM(6,2)*PSP(6,9)+ACM(7,2)*PSP(7,9)
     #+ACM(8,2)*PSP(8,9)+ACM(9,2)*PSP(9,9)
      t420 = ACM(1,3)*PSP(1,1)+ACM(2,3)*PSP(1,2)+ACM(3,3)*PSP(1,3)+ACM(4
     #,3)*PSP(1,4)+ACM(5,3)*PSP(1,5)+ACM(6,3)*PSP(1,6)+ACM(7,3)*PSP(1,7)
     #+ACM(8,3)*PSP(1,8)+ACM(9,3)*PSP(1,9)
      t431 = ACM(1,3)*PSP(1,2)+ACM(2,3)*PSP(2,2)+ACM(3,3)*PSP(2,3)+ACM(4
     #,3)*PSP(2,4)+ACM(5,3)*PSP(2,5)+ACM(6,3)*PSP(2,6)+ACM(7,3)*PSP(2,7)
     #+ACM(8,3)*PSP(2,8)+ACM(9,3)*PSP(2,9)
      t442 = ACM(1,3)*PSP(1,3)+ACM(2,3)*PSP(2,3)+ACM(3,3)*PSP(3,3)+ACM(4
     #,3)*PSP(3,4)+ACM(5,3)*PSP(3,5)+ACM(6,3)*PSP(3,6)+ACM(7,3)*PSP(3,7)
     #+ACM(8,3)*PSP(3,8)+ACM(9,3)*PSP(3,9)
      t453 = ACM(1,3)*PSP(1,4)+ACM(2,3)*PSP(2,4)+ACM(3,3)*PSP(3,4)+ACM(4
     #,3)*PSP(4,4)+ACM(5,3)*PSP(4,5)+ACM(6,3)*PSP(4,6)+ACM(7,3)*PSP(4,7)
     #+ACM(8,3)*PSP(4,8)+ACM(9,3)*PSP(4,9)
      t464 = ACM(1,3)*PSP(1,5)+ACM(2,3)*PSP(2,5)+ACM(3,3)*PSP(3,5)+ACM(4
     #,3)*PSP(4,5)+ACM(5,3)*PSP(5,5)+ACM(6,3)*PSP(5,6)+ACM(7,3)*PSP(5,7)
     #+ACM(8,3)*PSP(5,8)+ACM(9,3)*PSP(5,9)
      t475 = ACM(1,3)*PSP(1,6)+ACM(2,3)*PSP(2,6)+ACM(3,3)*PSP(3,6)+ACM(4
     #,3)*PSP(4,6)+ACM(5,3)*PSP(5,6)+ACM(6,3)*PSP(6,6)+ACM(7,3)*PSP(6,7)
     #+ACM(8,3)*PSP(6,8)+ACM(9,3)*PSP(6,9)
      t486 = ACM(1,3)*PSP(1,7)+ACM(2,3)*PSP(2,7)+ACM(3,3)*PSP(3,7)+ACM(4
     #,3)*PSP(4,7)+ACM(5,3)*PSP(5,7)+ACM(6,3)*PSP(6,7)+ACM(7,3)*PSP(7,7)
     #+ACM(8,3)*PSP(7,8)+ACM(9,3)*PSP(7,9)
      t497 = ACM(1,3)*PSP(1,8)+ACM(2,3)*PSP(2,8)+ACM(3,3)*PSP(3,8)+ACM(4
     #,3)*PSP(4,8)+ACM(5,3)*PSP(5,8)+ACM(6,3)*PSP(6,8)+ACM(7,3)*PSP(7,8)
     #+ACM(8,3)*PSP(8,8)+ACM(9,3)*PSP(8,9)
      t508 = ACM(1,3)*PSP(1,9)+ACM(2,3)*PSP(2,9)+ACM(3,3)*PSP(3,9)+ACM(4
     #,3)*PSP(4,9)+ACM(5,3)*PSP(5,9)+ACM(6,3)*PSP(6,9)+ACM(7,3)*PSP(7,9)
     #+ACM(8,3)*PSP(8,9)+ACM(9,3)*PSP(9,9)
      t610 = ACM(1,4)*PSP(1,1)+ACM(2,4)*PSP(1,2)+ACM(3,4)*PSP(1,3)+ACM(4
     #,4)*PSP(1,4)+ACM(5,4)*PSP(1,5)+ACM(6,4)*PSP(1,6)+ACM(7,4)*PSP(1,7)
     #+ACM(8,4)*PSP(1,8)+ACM(9,4)*PSP(1,9)
      t621 = ACM(1,4)*PSP(1,2)+ACM(2,4)*PSP(2,2)+ACM(3,4)*PSP(2,3)+ACM(4
     #,4)*PSP(2,4)+ACM(5,4)*PSP(2,5)+ACM(6,4)*PSP(2,6)+ACM(7,4)*PSP(2,7)
     #+ACM(8,4)*PSP(2,8)+ACM(9,4)*PSP(2,9)
      t632 = ACM(1,4)*PSP(1,3)+ACM(2,4)*PSP(2,3)+ACM(3,4)*PSP(3,3)+ACM(4
     #,4)*PSP(3,4)+ACM(5,4)*PSP(3,5)+ACM(6,4)*PSP(3,6)+ACM(7,4)*PSP(3,7)
     #+ACM(8,4)*PSP(3,8)+ACM(9,4)*PSP(3,9)
      t643 = ACM(1,4)*PSP(1,4)+ACM(2,4)*PSP(2,4)+ACM(3,4)*PSP(3,4)+ACM(4
     #,4)*PSP(4,4)+ACM(5,4)*PSP(4,5)+ACM(6,4)*PSP(4,6)+ACM(7,4)*PSP(4,7)
     #+ACM(8,4)*PSP(4,8)+ACM(9,4)*PSP(4,9)
      t654 = ACM(1,4)*PSP(1,5)+ACM(2,4)*PSP(2,5)+ACM(3,4)*PSP(3,5)+ACM(4
     #,4)*PSP(4,5)+ACM(5,4)*PSP(5,5)+ACM(6,4)*PSP(5,6)+ACM(7,4)*PSP(5,7)
     #+ACM(8,4)*PSP(5,8)+ACM(9,4)*PSP(5,9)
      t665 = ACM(1,4)*PSP(1,6)+ACM(2,4)*PSP(2,6)+ACM(3,4)*PSP(3,6)+ACM(4
     #,4)*PSP(4,6)+ACM(5,4)*PSP(5,6)+ACM(6,4)*PSP(6,6)+ACM(7,4)*PSP(6,7)
     #+ACM(8,4)*PSP(6,8)+ACM(9,4)*PSP(6,9)
      t676 = ACM(1,4)*PSP(1,7)+ACM(2,4)*PSP(2,7)+ACM(3,4)*PSP(3,7)+ACM(4
     #,4)*PSP(4,7)+ACM(5,4)*PSP(5,7)+ACM(6,4)*PSP(6,7)+ACM(7,4)*PSP(7,7)
     #+ACM(8,4)*PSP(7,8)+ACM(9,4)*PSP(7,9)
      t687 = ACM(1,4)*PSP(1,8)+ACM(2,4)*PSP(2,8)+ACM(3,4)*PSP(3,8)+ACM(4
     #,4)*PSP(4,8)+ACM(5,4)*PSP(5,8)+ACM(6,4)*PSP(6,8)+ACM(7,4)*PSP(7,8)
     #+ACM(8,4)*PSP(8,8)+ACM(9,4)*PSP(8,9)
      t698 = ACM(1,4)*PSP(1,9)+ACM(2,4)*PSP(2,9)+ACM(3,4)*PSP(3,9)+ACM(4
     #,4)*PSP(4,9)+ACM(5,4)*PSP(5,9)+ACM(6,4)*PSP(6,9)+ACM(7,4)*PSP(7,9)
     #+ACM(8,4)*PSP(8,9)+ACM(9,4)*PSP(9,9)
      t790 = ACM(1,5)*PSP(1,1)+ACM(2,5)*PSP(1,2)+ACM(3,5)*PSP(1,3)+ACM(4
     #,5)*PSP(1,4)+ACM(5,5)*PSP(1,5)+ACM(6,5)*PSP(1,6)+ACM(7,5)*PSP(1,7)
     #+ACM(8,5)*PSP(1,8)+ACM(9,5)*PSP(1,9)
      t801 = ACM(1,5)*PSP(1,2)+ACM(2,5)*PSP(2,2)+ACM(3,5)*PSP(2,3)+ACM(4
     #,5)*PSP(2,4)+ACM(5,5)*PSP(2,5)+ACM(6,5)*PSP(2,6)+ACM(7,5)*PSP(2,7)
     #+ACM(8,5)*PSP(2,8)+ACM(9,5)*PSP(2,9)
      t812 = ACM(1,5)*PSP(1,3)+ACM(2,5)*PSP(2,3)+ACM(3,5)*PSP(3,3)+ACM(4
     #,5)*PSP(3,4)+ACM(5,5)*PSP(3,5)+ACM(6,5)*PSP(3,6)+ACM(7,5)*PSP(3,7)
     #+ACM(8,5)*PSP(3,8)+ACM(9,5)*PSP(3,9)
      t823 = ACM(1,5)*PSP(1,4)+ACM(2,5)*PSP(2,4)+ACM(3,5)*PSP(3,4)+ACM(4
     #,5)*PSP(4,4)+ACM(5,5)*PSP(4,5)+ACM(6,5)*PSP(4,6)+ACM(7,5)*PSP(4,7)
     #+ACM(8,5)*PSP(4,8)+ACM(9,5)*PSP(4,9)
      t834 = ACM(1,5)*PSP(1,5)+ACM(2,5)*PSP(2,5)+ACM(3,5)*PSP(3,5)+ACM(4
     #,5)*PSP(4,5)+ACM(5,5)*PSP(5,5)+ACM(6,5)*PSP(5,6)+ACM(7,5)*PSP(5,7)
     #+ACM(8,5)*PSP(5,8)+ACM(9,5)*PSP(5,9)
      t845 = ACM(1,5)*PSP(1,6)+ACM(2,5)*PSP(2,6)+ACM(3,5)*PSP(3,6)+ACM(4
     #,5)*PSP(4,6)+ACM(5,5)*PSP(5,6)+ACM(6,5)*PSP(6,6)+ACM(7,5)*PSP(6,7)
     #+ACM(8,5)*PSP(6,8)+ACM(9,5)*PSP(6,9)
      t856 = ACM(1,5)*PSP(1,7)+ACM(2,5)*PSP(2,7)+ACM(3,5)*PSP(3,7)+ACM(4
     #,5)*PSP(4,7)+ACM(5,5)*PSP(5,7)+ACM(6,5)*PSP(6,7)+ACM(7,5)*PSP(7,7)
     #+ACM(8,5)*PSP(7,8)+ACM(9,5)*PSP(7,9)
      t867 = ACM(1,5)*PSP(1,8)+ACM(2,5)*PSP(2,8)+ACM(3,5)*PSP(3,8)+ACM(4
     #,5)*PSP(4,8)+ACM(5,5)*PSP(5,8)+ACM(6,5)*PSP(6,8)+ACM(7,5)*PSP(7,8)
     #+ACM(8,5)*PSP(8,8)+ACM(9,5)*PSP(8,9)
      t878 = ACM(1,5)*PSP(1,9)+ACM(2,5)*PSP(2,9)+ACM(3,5)*PSP(3,9)+ACM(4
     #,5)*PSP(4,9)+ACM(5,5)*PSP(5,9)+ACM(6,5)*PSP(6,9)+ACM(7,5)*PSP(7,9)
     #+ACM(8,5)*PSP(8,9)+ACM(9,5)*PSP(9,9)
      t960 = ACM(1,6)*PSP(1,1)+ACM(2,6)*PSP(1,2)+ACM(3,6)*PSP(1,3)+ACM(4
     #,6)*PSP(1,4)+ACM(5,6)*PSP(1,5)+ACM(6,6)*PSP(1,6)+ACM(7,6)*PSP(1,7)
     #+ACM(8,6)*PSP(1,8)+ACM(9,6)*PSP(1,9)
      t971 = ACM(1,6)*PSP(1,2)+ACM(2,6)*PSP(2,2)+ACM(3,6)*PSP(2,3)+ACM(4
     #,6)*PSP(2,4)+ACM(5,6)*PSP(2,5)+ACM(6,6)*PSP(2,6)+ACM(7,6)*PSP(2,7)
     #+ACM(8,6)*PSP(2,8)+ACM(9,6)*PSP(2,9)
      t982 = ACM(1,6)*PSP(1,3)+ACM(2,6)*PSP(2,3)+ACM(3,6)*PSP(3,3)+ACM(4
     #,6)*PSP(3,4)+ACM(5,6)*PSP(3,5)+ACM(6,6)*PSP(3,6)+ACM(7,6)*PSP(3,7)
     #+ACM(8,6)*PSP(3,8)+ACM(9,6)*PSP(3,9)
      t993 = ACM(1,6)*PSP(1,4)+ACM(2,6)*PSP(2,4)+ACM(3,6)*PSP(3,4)+ACM(4
     #,6)*PSP(4,4)+ACM(5,6)*PSP(4,5)+ACM(6,6)*PSP(4,6)+ACM(7,6)*PSP(4,7)
     #+ACM(8,6)*PSP(4,8)+ACM(9,6)*PSP(4,9)
      t1004 = ACM(1,6)*PSP(1,5)+ACM(2,6)*PSP(2,5)+ACM(3,6)*PSP(3,5)+ACM(
     #4,6)*PSP(4,5)+ACM(5,6)*PSP(5,5)+ACM(6,6)*PSP(5,6)+ACM(7,6)*PSP(5,7
     #)+ACM(8,6)*PSP(5,8)+ACM(9,6)*PSP(5,9)
      t1015 = ACM(1,6)*PSP(1,6)+ACM(2,6)*PSP(2,6)+ACM(3,6)*PSP(3,6)+ACM(
     #4,6)*PSP(4,6)+ACM(5,6)*PSP(5,6)+ACM(6,6)*PSP(6,6)+ACM(7,6)*PSP(6,7
     #)+ACM(8,6)*PSP(6,8)+ACM(9,6)*PSP(6,9)
      t1026 = ACM(1,6)*PSP(1,7)+ACM(2,6)*PSP(2,7)+ACM(3,6)*PSP(3,7)+ACM(
     #4,6)*PSP(4,7)+ACM(5,6)*PSP(5,7)+ACM(6,6)*PSP(6,7)+ACM(7,6)*PSP(7,7
     #)+ACM(8,6)*PSP(7,8)+ACM(9,6)*PSP(7,9)
      t1037 = ACM(1,6)*PSP(1,8)+ACM(2,6)*PSP(2,8)+ACM(3,6)*PSP(3,8)+ACM(
     #4,6)*PSP(4,8)+ACM(5,6)*PSP(5,8)+ACM(6,6)*PSP(6,8)+ACM(7,6)*PSP(7,8
     #)+ACM(8,6)*PSP(8,8)+ACM(9,6)*PSP(8,9)
      t1048 = ACM(1,6)*PSP(1,9)+ACM(2,6)*PSP(2,9)+ACM(3,6)*PSP(3,9)+ACM(
     #4,6)*PSP(4,9)+ACM(5,6)*PSP(5,9)+ACM(6,6)*PSP(6,9)+ACM(7,6)*PSP(7,9
     #)+ACM(8,6)*PSP(8,9)+ACM(9,6)*PSP(9,9)
      t1120 = ACM(1,7)*PSP(1,1)+ACM(2,7)*PSP(1,2)+ACM(3,7)*PSP(1,3)+ACM(
     #4,7)*PSP(1,4)+ACM(5,7)*PSP(1,5)+ACM(6,7)*PSP(1,6)+ACM(7,7)*PSP(1,7
     #)+ACM(8,7)*PSP(1,8)+ACM(9,7)*PSP(1,9)
      t1131 = ACM(1,7)*PSP(1,2)+ACM(2,7)*PSP(2,2)+ACM(3,7)*PSP(2,3)+ACM(
     #4,7)*PSP(2,4)+ACM(5,7)*PSP(2,5)+ACM(6,7)*PSP(2,6)+ACM(7,7)*PSP(2,7
     #)+ACM(8,7)*PSP(2,8)+ACM(9,7)*PSP(2,9)
      t1142 = ACM(1,7)*PSP(1,3)+ACM(2,7)*PSP(2,3)+ACM(3,7)*PSP(3,3)+ACM(
     #4,7)*PSP(3,4)+ACM(5,7)*PSP(3,5)+ACM(6,7)*PSP(3,6)+ACM(7,7)*PSP(3,7
     #)+ACM(8,7)*PSP(3,8)+ACM(9,7)*PSP(3,9)
      t1153 = ACM(1,7)*PSP(1,4)+ACM(2,7)*PSP(2,4)+ACM(3,7)*PSP(3,4)+ACM(
     #4,7)*PSP(4,4)+ACM(5,7)*PSP(4,5)+ACM(6,7)*PSP(4,6)+ACM(7,7)*PSP(4,7
     #)+ACM(8,7)*PSP(4,8)+ACM(9,7)*PSP(4,9)
      t1164 = ACM(1,7)*PSP(1,5)+ACM(2,7)*PSP(2,5)+ACM(3,7)*PSP(3,5)+ACM(
     #4,7)*PSP(4,5)+ACM(5,7)*PSP(5,5)+ACM(6,7)*PSP(5,6)+ACM(7,7)*PSP(5,7
     #)+ACM(8,7)*PSP(5,8)+ACM(9,7)*PSP(5,9)
      t1175 = ACM(1,7)*PSP(1,6)+ACM(2,7)*PSP(2,6)+ACM(3,7)*PSP(3,6)+ACM(
     #4,7)*PSP(4,6)+ACM(5,7)*PSP(5,6)+ACM(6,7)*PSP(6,6)+ACM(7,7)*PSP(6,7
     #)+ACM(8,7)*PSP(6,8)+ACM(9,7)*PSP(6,9)
      t1186 = ACM(1,7)*PSP(1,7)+ACM(2,7)*PSP(2,7)+ACM(3,7)*PSP(3,7)+ACM(
     #4,7)*PSP(4,7)+ACM(5,7)*PSP(5,7)+ACM(6,7)*PSP(6,7)+ACM(7,7)*PSP(7,7
     #)+ACM(8,7)*PSP(7,8)+ACM(9,7)*PSP(7,9)
      t1197 = ACM(1,7)*PSP(1,8)+ACM(2,7)*PSP(2,8)+ACM(3,7)*PSP(3,8)+ACM(
     #4,7)*PSP(4,8)+ACM(5,7)*PSP(5,8)+ACM(6,7)*PSP(6,8)+ACM(7,7)*PSP(7,8
     #)+ACM(8,7)*PSP(8,8)+ACM(9,7)*PSP(8,9)
      t1208 = ACM(1,7)*PSP(1,9)+ACM(2,7)*PSP(2,9)+ACM(3,7)*PSP(3,9)+ACM(
     #4,7)*PSP(4,9)+ACM(5,7)*PSP(5,9)+ACM(6,7)*PSP(6,9)+ACM(7,7)*PSP(7,9
     #)+ACM(8,7)*PSP(8,9)+ACM(9,7)*PSP(9,9)
      t1270 = ACM(1,8)*PSP(1,1)+ACM(2,8)*PSP(1,2)+ACM(3,8)*PSP(1,3)+ACM(
     #4,8)*PSP(1,4)+ACM(5,8)*PSP(1,5)+ACM(6,8)*PSP(1,6)+ACM(7,8)*PSP(1,7
     #)+ACM(8,8)*PSP(1,8)+ACM(9,8)*PSP(1,9)
      t1281 = ACM(1,8)*PSP(1,2)+ACM(2,8)*PSP(2,2)+ACM(3,8)*PSP(2,3)+ACM(
     #4,8)*PSP(2,4)+ACM(5,8)*PSP(2,5)+ACM(6,8)*PSP(2,6)+ACM(7,8)*PSP(2,7
     #)+ACM(8,8)*PSP(2,8)+ACM(9,8)*PSP(2,9)
      t1292 = ACM(1,8)*PSP(1,3)+ACM(2,8)*PSP(2,3)+ACM(3,8)*PSP(3,3)+ACM(
     #4,8)*PSP(3,4)+ACM(5,8)*PSP(3,5)+ACM(6,8)*PSP(3,6)+ACM(7,8)*PSP(3,7
     #)+ACM(8,8)*PSP(3,8)+ACM(9,8)*PSP(3,9)
      t1303 = ACM(1,8)*PSP(1,4)+ACM(2,8)*PSP(2,4)+ACM(3,8)*PSP(3,4)+ACM(
     #4,8)*PSP(4,4)+ACM(5,8)*PSP(4,5)+ACM(6,8)*PSP(4,6)+ACM(7,8)*PSP(4,7
     #)+ACM(8,8)*PSP(4,8)+ACM(9,8)*PSP(4,9)
      t1314 = ACM(1,8)*PSP(1,5)+ACM(2,8)*PSP(2,5)+ACM(3,8)*PSP(3,5)+ACM(
     #4,8)*PSP(4,5)+ACM(5,8)*PSP(5,5)+ACM(6,8)*PSP(5,6)+ACM(7,8)*PSP(5,7
     #)+ACM(8,8)*PSP(5,8)+ACM(9,8)*PSP(5,9)
      t1325 = ACM(1,8)*PSP(1,6)+ACM(2,8)*PSP(2,6)+ACM(3,8)*PSP(3,6)+ACM(
     #4,8)*PSP(4,6)+ACM(5,8)*PSP(5,6)+ACM(6,8)*PSP(6,6)+ACM(7,8)*PSP(6,7
     #)+ACM(8,8)*PSP(6,8)+ACM(9,8)*PSP(6,9)
      t1336 = ACM(1,8)*PSP(1,7)+ACM(2,8)*PSP(2,7)+ACM(3,8)*PSP(3,7)+ACM(
     #4,8)*PSP(4,7)+ACM(5,8)*PSP(5,7)+ACM(6,8)*PSP(6,7)+ACM(7,8)*PSP(7,7
     #)+ACM(8,8)*PSP(7,8)+ACM(9,8)*PSP(7,9)
      t1347 = ACM(1,8)*PSP(1,8)+ACM(2,8)*PSP(2,8)+ACM(3,8)*PSP(3,8)+ACM(
     #4,8)*PSP(4,8)+ACM(5,8)*PSP(5,8)+ACM(6,8)*PSP(6,8)+ACM(7,8)*PSP(7,8
     #)+ACM(8,8)*PSP(8,8)+ACM(9,8)*PSP(8,9)
      t1358 = ACM(1,8)*PSP(1,9)+ACM(2,8)*PSP(2,9)+ACM(3,8)*PSP(3,9)+ACM(
     #4,8)*PSP(4,9)+ACM(5,8)*PSP(5,9)+ACM(6,8)*PSP(6,9)+ACM(7,8)*PSP(7,9
     #)+ACM(8,8)*PSP(8,9)+ACM(9,8)*PSP(9,9)
      t1410 = ACM(1,9)*PSP(1,1)+ACM(2,9)*PSP(1,2)+ACM(3,9)*PSP(1,3)+ACM(
     #4,9)*PSP(1,4)+ACM(5,9)*PSP(1,5)+ACM(6,9)*PSP(1,6)+ACM(7,9)*PSP(1,7
     #)+ACM(8,9)*PSP(1,8)+ACM(9,9)*PSP(1,9)
      t1421 = ACM(1,9)*PSP(1,2)+ACM(2,9)*PSP(2,2)+ACM(3,9)*PSP(2,3)+ACM(
     #4,9)*PSP(2,4)+ACM(5,9)*PSP(2,5)+ACM(6,9)*PSP(2,6)+ACM(7,9)*PSP(2,7
     #)+ACM(8,9)*PSP(2,8)+ACM(9,9)*PSP(2,9)
      t1432 = ACM(1,9)*PSP(1,3)+ACM(2,9)*PSP(2,3)+ACM(3,9)*PSP(3,3)+ACM(
     #4,9)*PSP(3,4)+ACM(5,9)*PSP(3,5)+ACM(6,9)*PSP(3,6)+ACM(7,9)*PSP(3,7
     #)+ACM(8,9)*PSP(3,8)+ACM(9,9)*PSP(3,9)
      t1443 = ACM(1,9)*PSP(1,4)+ACM(2,9)*PSP(2,4)+ACM(3,9)*PSP(3,4)+ACM(
     #4,9)*PSP(4,4)+ACM(5,9)*PSP(4,5)+ACM(6,9)*PSP(4,6)+ACM(7,9)*PSP(4,7
     #)+ACM(8,9)*PSP(4,8)+ACM(9,9)*PSP(4,9)
      t1454 = ACM(1,9)*PSP(1,5)+ACM(2,9)*PSP(2,5)+ACM(3,9)*PSP(3,5)+ACM(
     #4,9)*PSP(4,5)+ACM(5,9)*PSP(5,5)+ACM(6,9)*PSP(5,6)+ACM(7,9)*PSP(5,7
     #)+ACM(8,9)*PSP(5,8)+ACM(9,9)*PSP(5,9)
      t1465 = ACM(1,9)*PSP(1,6)+ACM(2,9)*PSP(2,6)+ACM(3,9)*PSP(3,6)+ACM(
     #4,9)*PSP(4,6)+ACM(5,9)*PSP(5,6)+ACM(6,9)*PSP(6,6)+ACM(7,9)*PSP(6,7
     #)+ACM(8,9)*PSP(6,8)+ACM(9,9)*PSP(6,9)
      t1476 = ACM(1,9)*PSP(1,7)+ACM(2,9)*PSP(2,7)+ACM(3,9)*PSP(3,7)+ACM(
     #4,9)*PSP(4,7)+ACM(5,9)*PSP(5,7)+ACM(6,9)*PSP(6,7)+ACM(7,9)*PSP(7,7
     #)+ACM(8,9)*PSP(7,8)+ACM(9,9)*PSP(7,9)
      t1487 = ACM(1,9)*PSP(1,8)+ACM(2,9)*PSP(2,8)+ACM(3,9)*PSP(3,8)+ACM(
     #4,9)*PSP(4,8)+ACM(5,9)*PSP(5,8)+ACM(6,9)*PSP(6,8)+ACM(7,9)*PSP(7,8
     #)+ACM(8,9)*PSP(8,8)+ACM(9,9)*PSP(8,9)
      t1498 = ACM(1,9)*PSP(1,9)+ACM(2,9)*PSP(2,9)+ACM(3,9)*PSP(3,9)+ACM(
     #4,9)*PSP(4,9)+ACM(5,9)*PSP(5,9)+ACM(6,9)*PSP(6,9)+ACM(7,9)*PSP(7,9
     #)+ACM(8,9)*PSP(8,9)+ACM(9,9)*PSP(9,9)
      t1540 = ACM(1,10)*PSP(1,1)+ACM(2,10)*PSP(1,2)+ACM(3,10)*PSP(1,3)+A
     #CM(4,10)*PSP(1,4)+ACM(5,10)*PSP(1,5)+ACM(6,10)*PSP(1,6)+ACM(7,10)*
     #PSP(1,7)+ACM(8,10)*PSP(1,8)+ACM(9,10)*PSP(1,9)
      t1551 = ACM(1,10)*PSP(1,2)+ACM(2,10)*PSP(2,2)+ACM(3,10)*PSP(2,3)+A
     #CM(4,10)*PSP(2,4)+ACM(5,10)*PSP(2,5)+ACM(6,10)*PSP(2,6)+ACM(7,10)*
     #PSP(2,7)+ACM(8,10)*PSP(2,8)+ACM(9,10)*PSP(2,9)
      t1562 = ACM(1,10)*PSP(1,3)+ACM(2,10)*PSP(2,3)+ACM(3,10)*PSP(3,3)+A
     #CM(4,10)*PSP(3,4)+ACM(5,10)*PSP(3,5)+ACM(6,10)*PSP(3,6)+ACM(7,10)*
     #PSP(3,7)+ACM(8,10)*PSP(3,8)+ACM(9,10)*PSP(3,9)
      t1573 = ACM(1,10)*PSP(1,4)+ACM(2,10)*PSP(2,4)+ACM(3,10)*PSP(3,4)+A
     #CM(4,10)*PSP(4,4)+ACM(5,10)*PSP(4,5)+ACM(6,10)*PSP(4,6)+ACM(7,10)*
     #PSP(4,7)+ACM(8,10)*PSP(4,8)+ACM(9,10)*PSP(4,9)
      t1584 = ACM(1,10)*PSP(1,5)+ACM(2,10)*PSP(2,5)+ACM(3,10)*PSP(3,5)+A
     #CM(4,10)*PSP(4,5)+ACM(5,10)*PSP(5,5)+ACM(6,10)*PSP(5,6)+ACM(7,10)*
     #PSP(5,7)+ACM(8,10)*PSP(5,8)+ACM(9,10)*PSP(5,9)
      t1595 = ACM(1,10)*PSP(1,6)+ACM(2,10)*PSP(2,6)+ACM(3,10)*PSP(3,6)+A
     #CM(4,10)*PSP(4,6)+ACM(5,10)*PSP(5,6)+ACM(6,10)*PSP(6,6)+ACM(7,10)*
     #PSP(6,7)+ACM(8,10)*PSP(6,8)+ACM(9,10)*PSP(6,9)
      t1606 = ACM(1,10)*PSP(1,7)+ACM(2,10)*PSP(2,7)+ACM(3,10)*PSP(3,7)+A
     #CM(4,10)*PSP(4,7)+ACM(5,10)*PSP(5,7)+ACM(6,10)*PSP(6,7)+ACM(7,10)*
     #PSP(7,7)+ACM(8,10)*PSP(7,8)+ACM(9,10)*PSP(7,9)
      t1617 = ACM(1,10)*PSP(1,8)+ACM(2,10)*PSP(2,8)+ACM(3,10)*PSP(3,8)+A
     #CM(4,10)*PSP(4,8)+ACM(5,10)*PSP(5,8)+ACM(6,10)*PSP(6,8)+ACM(7,10)*
     #PSP(7,8)+ACM(8,10)*PSP(8,8)+ACM(9,10)*PSP(8,9)
      t1628 = ACM(1,10)*PSP(1,9)+ACM(2,10)*PSP(2,9)+ACM(3,10)*PSP(3,9)+A
     #CM(4,10)*PSP(4,9)+ACM(5,10)*PSP(5,9)+ACM(6,10)*PSP(6,9)+ACM(7,10)*
     #PSP(7,9)+ACM(8,10)*PSP(8,9)+ACM(9,10)*PSP(9,9)
      t1660 = ACM(1,11)*PSP(1,1)+ACM(2,11)*PSP(1,2)+ACM(3,11)*PSP(1,3)+A
     #CM(4,11)*PSP(1,4)+ACM(5,11)*PSP(1,5)+ACM(6,11)*PSP(1,6)+ACM(7,11)*
     #PSP(1,7)+ACM(8,11)*PSP(1,8)+ACM(9,11)*PSP(1,9)
      t1671 = ACM(1,11)*PSP(1,2)+ACM(2,11)*PSP(2,2)+ACM(3,11)*PSP(2,3)+A
     #CM(4,11)*PSP(2,4)+ACM(5,11)*PSP(2,5)+ACM(6,11)*PSP(2,6)+ACM(7,11)*
     #PSP(2,7)+ACM(8,11)*PSP(2,8)+ACM(9,11)*PSP(2,9)
      t1682 = ACM(1,11)*PSP(1,3)+ACM(2,11)*PSP(2,3)+ACM(3,11)*PSP(3,3)+A
     #CM(4,11)*PSP(3,4)+ACM(5,11)*PSP(3,5)+ACM(6,11)*PSP(3,6)+ACM(7,11)*
     #PSP(3,7)+ACM(8,11)*PSP(3,8)+ACM(9,11)*PSP(3,9)
      t1693 = ACM(1,11)*PSP(1,4)+ACM(2,11)*PSP(2,4)+ACM(3,11)*PSP(3,4)+A
     #CM(4,11)*PSP(4,4)+ACM(5,11)*PSP(4,5)+ACM(6,11)*PSP(4,6)+ACM(7,11)*
     #PSP(4,7)+ACM(8,11)*PSP(4,8)+ACM(9,11)*PSP(4,9)
      t1704 = ACM(1,11)*PSP(1,5)+ACM(2,11)*PSP(2,5)+ACM(3,11)*PSP(3,5)+A
     #CM(4,11)*PSP(4,5)+ACM(5,11)*PSP(5,5)+ACM(6,11)*PSP(5,6)+ACM(7,11)*
     #PSP(5,7)+ACM(8,11)*PSP(5,8)+ACM(9,11)*PSP(5,9)
      t1715 = ACM(1,11)*PSP(1,6)+ACM(2,11)*PSP(2,6)+ACM(3,11)*PSP(3,6)+A
     #CM(4,11)*PSP(4,6)+ACM(5,11)*PSP(5,6)+ACM(6,11)*PSP(6,6)+ACM(7,11)*
     #PSP(6,7)+ACM(8,11)*PSP(6,8)+ACM(9,11)*PSP(6,9)
      t1726 = ACM(1,11)*PSP(1,7)+ACM(2,11)*PSP(2,7)+ACM(3,11)*PSP(3,7)+A
     #CM(4,11)*PSP(4,7)+ACM(5,11)*PSP(5,7)+ACM(6,11)*PSP(6,7)+ACM(7,11)*
     #PSP(7,7)+ACM(8,11)*PSP(7,8)+ACM(9,11)*PSP(7,9)
      t1737 = ACM(1,11)*PSP(1,8)+ACM(2,11)*PSP(2,8)+ACM(3,11)*PSP(3,8)+A
     #CM(4,11)*PSP(4,8)+ACM(5,11)*PSP(5,8)+ACM(6,11)*PSP(6,8)+ACM(7,11)*
     #PSP(7,8)+ACM(8,11)*PSP(8,8)+ACM(9,11)*PSP(8,9)
      t1748 = ACM(1,11)*PSP(1,9)+ACM(2,11)*PSP(2,9)+ACM(3,11)*PSP(3,9)+A
     #CM(4,11)*PSP(4,9)+ACM(5,11)*PSP(5,9)+ACM(6,11)*PSP(6,9)+ACM(7,11)*
     #PSP(7,9)+ACM(8,11)*PSP(8,9)+ACM(9,11)*PSP(9,9)
      ST(1) = t10*ACM(1,1)+t21*ACM(2,1)+t32*ACM(3,1)+t43*ACM(4,1)+t54*AC
     #M(5,1)+t65*ACM(6,1)+t76*ACM(7,1)+t87*ACM(8,1)+t98*ACM(9,1)
      ST(2) = t10*ACM(1,2)+t21*ACM(2,2)+t32*ACM(3,2)+t43*ACM(4,2)+t54*AC
     #M(5,2)+t65*ACM(6,2)+t76*ACM(7,2)+t87*ACM(8,2)+t98*ACM(9,2)
      ST(3) = t10*ACM(1,3)+t21*ACM(2,3)+t32*ACM(3,3)+t43*ACM(4,3)+t54*AC
     #M(5,3)+t65*ACM(6,3)+t76*ACM(7,3)+t87*ACM(8,3)+t98*ACM(9,3)
      ST(4) = t10*ACM(1,4)+t21*ACM(2,4)+t32*ACM(3,4)+t43*ACM(4,4)+t54*AC
     #M(5,4)+t65*ACM(6,4)+t76*ACM(7,4)+t87*ACM(8,4)+t98*ACM(9,4)
      ST(5) = t10*ACM(1,5)+t21*ACM(2,5)+t32*ACM(3,5)+t43*ACM(4,5)+t54*AC
     #M(5,5)+t65*ACM(6,5)+t76*ACM(7,5)+t87*ACM(8,5)+t98*ACM(9,5)
      ST(6) = t10*ACM(1,6)+t21*ACM(2,6)+t32*ACM(3,6)+t43*ACM(4,6)+t54*AC
     #M(5,6)+t65*ACM(6,6)+t76*ACM(7,6)+t87*ACM(8,6)+t98*ACM(9,6)
      ST(7) = t10*ACM(1,7)+t21*ACM(2,7)+t32*ACM(3,7)+t43*ACM(4,7)+t54*AC
     #M(5,7)+t65*ACM(6,7)+t76*ACM(7,7)+t87*ACM(8,7)+t98*ACM(9,7)
      ST(8) = t10*ACM(1,8)+t21*ACM(2,8)+t32*ACM(3,8)+t43*ACM(4,8)+t54*AC
     #M(5,8)+t65*ACM(6,8)+t76*ACM(7,8)+t87*ACM(8,8)+t98*ACM(9,8)
      ST(9) = t10*ACM(1,9)+t21*ACM(2,9)+t32*ACM(3,9)+t43*ACM(4,9)+t54*AC
     #M(5,9)+t65*ACM(6,9)+t76*ACM(7,9)+t87*ACM(8,9)+t98*ACM(9,9)
      ST(10) = t10*ACM(1,10)+t21*ACM(2,10)+t32*ACM(3,10)+t43*ACM(4,10)+t
     #54*ACM(5,10)+t65*ACM(6,10)+t76*ACM(7,10)+t87*ACM(8,10)+t98*ACM(9,1
     #0)
      ST(11) = t10*ACM(1,11)+t21*ACM(2,11)+t32*ACM(3,11)+t43*ACM(4,11)+t
     #54*ACM(5,11)+t65*ACM(6,11)+t76*ACM(7,11)+t87*ACM(8,11)+t98*ACM(9,1
     #1)
      ST(12) = t10*ACM(1,12)+t21*ACM(2,12)+t32*ACM(3,12)+t43*ACM(4,12)+t
     #54*ACM(5,12)+t65*ACM(6,12)+t76*ACM(7,12)+t87*ACM(8,12)+t98*ACM(9,1
     #2)
      ST(13) = t220*ACM(1,2)+t231*ACM(2,2)+t242*ACM(3,2)+t253*ACM(4,2)+t
     #264*ACM(5,2)+t275*ACM(6,2)+t286*ACM(7,2)+t297*ACM(8,2)+t308*ACM(9,
     #2)
      ST(14) = t220*ACM(1,3)+t231*ACM(2,3)+t242*ACM(3,3)+t253*ACM(4,3)+t
     #264*ACM(5,3)+t275*ACM(6,3)+t286*ACM(7,3)+t297*ACM(8,3)+t308*ACM(9,
     #3)
      ST(15) = t220*ACM(1,4)+t231*ACM(2,4)+t242*ACM(3,4)+t253*ACM(4,4)+t
     #264*ACM(5,4)+t275*ACM(6,4)+t286*ACM(7,4)+t297*ACM(8,4)+t308*ACM(9,
     #4)
      ST(16) = t220*ACM(1,5)+t231*ACM(2,5)+t242*ACM(3,5)+t253*ACM(4,5)+t
     #264*ACM(5,5)+t275*ACM(6,5)+t286*ACM(7,5)+t297*ACM(8,5)+t308*ACM(9,
     #5)
      ST(17) = t220*ACM(1,6)+t231*ACM(2,6)+t242*ACM(3,6)+t253*ACM(4,6)+t
     #264*ACM(5,6)+t275*ACM(6,6)+t286*ACM(7,6)+t297*ACM(8,6)+t308*ACM(9,
     #6)
      ST(18) = t220*ACM(1,7)+t231*ACM(2,7)+t242*ACM(3,7)+t253*ACM(4,7)+t
     #264*ACM(5,7)+t275*ACM(6,7)+t286*ACM(7,7)+t297*ACM(8,7)+t308*ACM(9,
     #7)
      ST(19) = t220*ACM(1,8)+t231*ACM(2,8)+t242*ACM(3,8)+t253*ACM(4,8)+t
     #264*ACM(5,8)+t275*ACM(6,8)+t286*ACM(7,8)+t297*ACM(8,8)+t308*ACM(9,
     #8)
      ST(20) = t220*ACM(1,9)+t231*ACM(2,9)+t242*ACM(3,9)+t253*ACM(4,9)+t
     #264*ACM(5,9)+t275*ACM(6,9)+t286*ACM(7,9)+t297*ACM(8,9)+t308*ACM(9,
     #9)
      ST(21) = t220*ACM(1,10)+t231*ACM(2,10)+t242*ACM(3,10)+t253*ACM(4,1
     #0)+t264*ACM(5,10)+t275*ACM(6,10)+t286*ACM(7,10)+t297*ACM(8,10)+t30
     #8*ACM(9,10)
      ST(22) = t220*ACM(1,11)+t231*ACM(2,11)+t242*ACM(3,11)+t253*ACM(4,1
     #1)+t264*ACM(5,11)+t275*ACM(6,11)+t286*ACM(7,11)+t297*ACM(8,11)+t30
     #8*ACM(9,11)
      ST(23) = t220*ACM(1,12)+t231*ACM(2,12)+t242*ACM(3,12)+t253*ACM(4,1
     #2)+t264*ACM(5,12)+t275*ACM(6,12)+t286*ACM(7,12)+t297*ACM(8,12)+t30
     #8*ACM(9,12)
      ST(24) = t420*ACM(1,3)+t431*ACM(2,3)+t442*ACM(3,3)+t453*ACM(4,3)+t
     #464*ACM(5,3)+t475*ACM(6,3)+t486*ACM(7,3)+t497*ACM(8,3)+t508*ACM(9,
     #3)
      ST(25) = t420*ACM(1,4)+t431*ACM(2,4)+t442*ACM(3,4)+t453*ACM(4,4)+t
     #464*ACM(5,4)+t475*ACM(6,4)+t486*ACM(7,4)+t497*ACM(8,4)+t508*ACM(9,
     #4)
      ST(26) = t420*ACM(1,5)+t431*ACM(2,5)+t442*ACM(3,5)+t453*ACM(4,5)+t
     #464*ACM(5,5)+t475*ACM(6,5)+t486*ACM(7,5)+t497*ACM(8,5)+t508*ACM(9,
     #5)
      ST(27) = t420*ACM(1,6)+t431*ACM(2,6)+t442*ACM(3,6)+t453*ACM(4,6)+t
     #464*ACM(5,6)+t475*ACM(6,6)+t486*ACM(7,6)+t497*ACM(8,6)+t508*ACM(9,
     #6)
      ST(28) = t420*ACM(1,7)+t431*ACM(2,7)+t442*ACM(3,7)+t453*ACM(4,7)+t
     #464*ACM(5,7)+t475*ACM(6,7)+t486*ACM(7,7)+t497*ACM(8,7)+t508*ACM(9,
     #7)
      ST(29) = t420*ACM(1,8)+t431*ACM(2,8)+t442*ACM(3,8)+t453*ACM(4,8)+t
     #464*ACM(5,8)+t475*ACM(6,8)+t486*ACM(7,8)+t497*ACM(8,8)+t508*ACM(9,
     #8)
      ST(30) = t420*ACM(1,9)+t431*ACM(2,9)+t442*ACM(3,9)+t453*ACM(4,9)+t
     #464*ACM(5,9)+t475*ACM(6,9)+t486*ACM(7,9)+t497*ACM(8,9)+t508*ACM(9,
     #9)
      ST(31) = t420*ACM(1,10)+t431*ACM(2,10)+t442*ACM(3,10)+t453*ACM(4,1
     #0)+t464*ACM(5,10)+t475*ACM(6,10)+t486*ACM(7,10)+t497*ACM(8,10)+t50
     #8*ACM(9,10)
      ST(32) = t420*ACM(1,11)+t431*ACM(2,11)+t442*ACM(3,11)+t453*ACM(4,1
     #1)+t464*ACM(5,11)+t475*ACM(6,11)+t486*ACM(7,11)+t497*ACM(8,11)+t50
     #8*ACM(9,11)
      ST(33) = t420*ACM(1,12)+t431*ACM(2,12)+t442*ACM(3,12)+t453*ACM(4,1
     #2)+t464*ACM(5,12)+t475*ACM(6,12)+t486*ACM(7,12)+t497*ACM(8,12)+t50
     #8*ACM(9,12)
      ST(34) = t610*ACM(1,4)+t621*ACM(2,4)+t632*ACM(3,4)+t643*ACM(4,4)+t
     #654*ACM(5,4)+t665*ACM(6,4)+t676*ACM(7,4)+t687*ACM(8,4)+t698*ACM(9,
     #4)
      ST(35) = t610*ACM(1,5)+t621*ACM(2,5)+t632*ACM(3,5)+t643*ACM(4,5)+t
     #654*ACM(5,5)+t665*ACM(6,5)+t676*ACM(7,5)+t687*ACM(8,5)+t698*ACM(9,
     #5)
      ST(36) = t610*ACM(1,6)+t621*ACM(2,6)+t632*ACM(3,6)+t643*ACM(4,6)+t
     #654*ACM(5,6)+t665*ACM(6,6)+t676*ACM(7,6)+t687*ACM(8,6)+t698*ACM(9,
     #6)
      ST(37) = t610*ACM(1,7)+t621*ACM(2,7)+t632*ACM(3,7)+t643*ACM(4,7)+t
     #654*ACM(5,7)+t665*ACM(6,7)+t676*ACM(7,7)+t687*ACM(8,7)+t698*ACM(9,
     #7)
      ST(38) = t610*ACM(1,8)+t621*ACM(2,8)+t632*ACM(3,8)+t643*ACM(4,8)+t
     #654*ACM(5,8)+t665*ACM(6,8)+t676*ACM(7,8)+t687*ACM(8,8)+t698*ACM(9,
     #8)
      ST(39) = t610*ACM(1,9)+t621*ACM(2,9)+t632*ACM(3,9)+t643*ACM(4,9)+t
     #654*ACM(5,9)+t665*ACM(6,9)+t676*ACM(7,9)+t687*ACM(8,9)+t698*ACM(9,
     #9)
      ST(40) = t610*ACM(1,10)+t621*ACM(2,10)+t632*ACM(3,10)+t643*ACM(4,1
     #0)+t654*ACM(5,10)+t665*ACM(6,10)+t676*ACM(7,10)+t687*ACM(8,10)+t69
     #8*ACM(9,10)
      ST(41) = t610*ACM(1,11)+t621*ACM(2,11)+t632*ACM(3,11)+t643*ACM(4,1
     #1)+t654*ACM(5,11)+t665*ACM(6,11)+t676*ACM(7,11)+t687*ACM(8,11)+t69
     #8*ACM(9,11)
      ST(42) = t610*ACM(1,12)+t621*ACM(2,12)+t632*ACM(3,12)+t643*ACM(4,1
     #2)+t654*ACM(5,12)+t665*ACM(6,12)+t676*ACM(7,12)+t687*ACM(8,12)+t69
     #8*ACM(9,12)
      ST(43) = t790*ACM(1,5)+t801*ACM(2,5)+t812*ACM(3,5)+t823*ACM(4,5)+t
     #834*ACM(5,5)+t845*ACM(6,5)+t856*ACM(7,5)+t867*ACM(8,5)+t878*ACM(9,
     #5)
      ST(44) = t790*ACM(1,6)+t801*ACM(2,6)+t812*ACM(3,6)+t823*ACM(4,6)+t
     #834*ACM(5,6)+t845*ACM(6,6)+t856*ACM(7,6)+t867*ACM(8,6)+t878*ACM(9,
     #6)
      ST(45) = t790*ACM(1,7)+t801*ACM(2,7)+t812*ACM(3,7)+t823*ACM(4,7)+t
     #834*ACM(5,7)+t845*ACM(6,7)+t856*ACM(7,7)+t867*ACM(8,7)+t878*ACM(9,
     #7)
      ST(46) = t790*ACM(1,8)+t801*ACM(2,8)+t812*ACM(3,8)+t823*ACM(4,8)+t
     #834*ACM(5,8)+t845*ACM(6,8)+t856*ACM(7,8)+t867*ACM(8,8)+t878*ACM(9,
     #8)
      ST(47) = t790*ACM(1,9)+t801*ACM(2,9)+t812*ACM(3,9)+t823*ACM(4,9)+t
     #834*ACM(5,9)+t845*ACM(6,9)+t856*ACM(7,9)+t867*ACM(8,9)+t878*ACM(9,
     #9)
      ST(48) = t790*ACM(1,10)+t801*ACM(2,10)+t812*ACM(3,10)+t823*ACM(4,1
     #0)+t834*ACM(5,10)+t845*ACM(6,10)+t856*ACM(7,10)+t867*ACM(8,10)+t87
     #8*ACM(9,10)
      ST(49) = t790*ACM(1,11)+t801*ACM(2,11)+t812*ACM(3,11)+t823*ACM(4,1
     #1)+t834*ACM(5,11)+t845*ACM(6,11)+t856*ACM(7,11)+t867*ACM(8,11)+t87
     #8*ACM(9,11)
      ST(50) = t790*ACM(1,12)+t801*ACM(2,12)+t812*ACM(3,12)+t823*ACM(4,1
     #2)+t834*ACM(5,12)+t845*ACM(6,12)+t856*ACM(7,12)+t867*ACM(8,12)+t87
     #8*ACM(9,12)
      ST(51) = t960*ACM(1,6)+t971*ACM(2,6)+t982*ACM(3,6)+t993*ACM(4,6)+t
     #1004*ACM(5,6)+t1015*ACM(6,6)+t1026*ACM(7,6)+t1037*ACM(8,6)+t1048*A
     #CM(9,6)
      ST(52) = t960*ACM(1,7)+t971*ACM(2,7)+t982*ACM(3,7)+t993*ACM(4,7)+t
     #1004*ACM(5,7)+t1015*ACM(6,7)+t1026*ACM(7,7)+t1037*ACM(8,7)+t1048*A
     #CM(9,7)
      ST(53) = t960*ACM(1,8)+t971*ACM(2,8)+t982*ACM(3,8)+t993*ACM(4,8)+t
     #1004*ACM(5,8)+t1015*ACM(6,8)+t1026*ACM(7,8)+t1037*ACM(8,8)+t1048*A
     #CM(9,8)
      ST(54) = t960*ACM(1,9)+t971*ACM(2,9)+t982*ACM(3,9)+t993*ACM(4,9)+t
     #1004*ACM(5,9)+t1015*ACM(6,9)+t1026*ACM(7,9)+t1037*ACM(8,9)+t1048*A
     #CM(9,9)
      ST(55) = t960*ACM(1,10)+t971*ACM(2,10)+t982*ACM(3,10)+t993*ACM(4,1
     #0)+t1004*ACM(5,10)+t1015*ACM(6,10)+t1026*ACM(7,10)+t1037*ACM(8,10)
     #+t1048*ACM(9,10)
      ST(56) = t960*ACM(1,11)+t971*ACM(2,11)+t982*ACM(3,11)+t993*ACM(4,1
     #1)+t1004*ACM(5,11)+t1015*ACM(6,11)+t1026*ACM(7,11)+t1037*ACM(8,11)
     #+t1048*ACM(9,11)
      ST(57) = t960*ACM(1,12)+t971*ACM(2,12)+t982*ACM(3,12)+t993*ACM(4,1
     #2)+t1004*ACM(5,12)+t1015*ACM(6,12)+t1026*ACM(7,12)+t1037*ACM(8,12)
     #+t1048*ACM(9,12)
      ST(58) = t1120*ACM(1,7)+t1131*ACM(2,7)+t1142*ACM(3,7)+t1153*ACM(4,
     #7)+t1164*ACM(5,7)+t1175*ACM(6,7)+t1186*ACM(7,7)+t1197*ACM(8,7)+t12
     #08*ACM(9,7)
      ST(59) = t1120*ACM(1,8)+t1131*ACM(2,8)+t1142*ACM(3,8)+t1153*ACM(4,
     #8)+t1164*ACM(5,8)+t1175*ACM(6,8)+t1186*ACM(7,8)+t1197*ACM(8,8)+t12
     #08*ACM(9,8)
      ST(60) = t1120*ACM(1,9)+t1131*ACM(2,9)+t1142*ACM(3,9)+t1153*ACM(4,
     #9)+t1164*ACM(5,9)+t1175*ACM(6,9)+t1186*ACM(7,9)+t1197*ACM(8,9)+t12
     #08*ACM(9,9)
      ST(61) = t1120*ACM(1,10)+t1131*ACM(2,10)+t1142*ACM(3,10)+t1153*ACM
     #(4,10)+t1164*ACM(5,10)+t1175*ACM(6,10)+t1186*ACM(7,10)+t1197*ACM(8
     #,10)+t1208*ACM(9,10)
      ST(62) = t1120*ACM(1,11)+t1131*ACM(2,11)+t1142*ACM(3,11)+t1153*ACM
     #(4,11)+t1164*ACM(5,11)+t1175*ACM(6,11)+t1186*ACM(7,11)+t1197*ACM(8
     #,11)+t1208*ACM(9,11)
      ST(63) = t1120*ACM(1,12)+t1131*ACM(2,12)+t1142*ACM(3,12)+t1153*ACM
     #(4,12)+t1164*ACM(5,12)+t1175*ACM(6,12)+t1186*ACM(7,12)+t1197*ACM(8
     #,12)+t1208*ACM(9,12)
      ST(64) = t1270*ACM(1,8)+t1281*ACM(2,8)+t1292*ACM(3,8)+t1303*ACM(4,
     #8)+t1314*ACM(5,8)+t1325*ACM(6,8)+t1336*ACM(7,8)+t1347*ACM(8,8)+t13
     #58*ACM(9,8)
      ST(65) = t1270*ACM(1,9)+t1281*ACM(2,9)+t1292*ACM(3,9)+t1303*ACM(4,
     #9)+t1314*ACM(5,9)+t1325*ACM(6,9)+t1336*ACM(7,9)+t1347*ACM(8,9)+t13
     #58*ACM(9,9)
      ST(66) = t1270*ACM(1,10)+t1281*ACM(2,10)+t1292*ACM(3,10)+t1303*ACM
     #(4,10)+t1314*ACM(5,10)+t1325*ACM(6,10)+t1336*ACM(7,10)+t1347*ACM(8
     #,10)+t1358*ACM(9,10)
      ST(67) = t1270*ACM(1,11)+t1281*ACM(2,11)+t1292*ACM(3,11)+t1303*ACM
     #(4,11)+t1314*ACM(5,11)+t1325*ACM(6,11)+t1336*ACM(7,11)+t1347*ACM(8
     #,11)+t1358*ACM(9,11)
      ST(68) = t1270*ACM(1,12)+t1281*ACM(2,12)+t1292*ACM(3,12)+t1303*ACM
     #(4,12)+t1314*ACM(5,12)+t1325*ACM(6,12)+t1336*ACM(7,12)+t1347*ACM(8
     #,12)+t1358*ACM(9,12)
      ST(69) = t1410*ACM(1,9)+t1421*ACM(2,9)+t1432*ACM(3,9)+t1443*ACM(4,
     #9)+t1454*ACM(5,9)+t1465*ACM(6,9)+t1476*ACM(7,9)+t1487*ACM(8,9)+t14
     #98*ACM(9,9)
      ST(70) = t1410*ACM(1,10)+t1421*ACM(2,10)+t1432*ACM(3,10)+t1443*ACM
     #(4,10)+t1454*ACM(5,10)+t1465*ACM(6,10)+t1476*ACM(7,10)+t1487*ACM(8
     #,10)+t1498*ACM(9,10)
      ST(71) = t1410*ACM(1,11)+t1421*ACM(2,11)+t1432*ACM(3,11)+t1443*ACM
     #(4,11)+t1454*ACM(5,11)+t1465*ACM(6,11)+t1476*ACM(7,11)+t1487*ACM(8
     #,11)+t1498*ACM(9,11)
      ST(72) = t1410*ACM(1,12)+t1421*ACM(2,12)+t1432*ACM(3,12)+t1443*ACM
     #(4,12)+t1454*ACM(5,12)+t1465*ACM(6,12)+t1476*ACM(7,12)+t1487*ACM(8
     #,12)+t1498*ACM(9,12)
      ST(73) = t1540*ACM(1,10)+t1551*ACM(2,10)+t1562*ACM(3,10)+t1573*ACM
     #(4,10)+t1584*ACM(5,10)+t1595*ACM(6,10)+t1606*ACM(7,10)+t1617*ACM(8
     #,10)+t1628*ACM(9,10)
      ST(74) = t1540*ACM(1,11)+t1551*ACM(2,11)+t1562*ACM(3,11)+t1573*ACM
     #(4,11)+t1584*ACM(5,11)+t1595*ACM(6,11)+t1606*ACM(7,11)+t1617*ACM(8
     #,11)+t1628*ACM(9,11)
      ST(75) = t1540*ACM(1,12)+t1551*ACM(2,12)+t1562*ACM(3,12)+t1573*ACM
     #(4,12)+t1584*ACM(5,12)+t1595*ACM(6,12)+t1606*ACM(7,12)+t1617*ACM(8
     #,12)+t1628*ACM(9,12)
      ST(76) = t1660*ACM(1,11)+t1671*ACM(2,11)+t1682*ACM(3,11)+t1693*ACM
     #(4,11)+t1704*ACM(5,11)+t1715*ACM(6,11)+t1726*ACM(7,11)+t1737*ACM(8
     #,11)+t1748*ACM(9,11)
      ST(77) = t1660*ACM(1,12)+t1671*ACM(2,12)+t1682*ACM(3,12)+t1693*ACM
     #(4,12)+t1704*ACM(5,12)+t1715*ACM(6,12)+t1726*ACM(7,12)+t1737*ACM(8
     #,12)+t1748*ACM(9,12)
      s2 = (ACM(1,12)*PSP(1,1)+ACM(2,12)*PSP(1,2)+ACM(3,12)*PSP(1,3)+ACM
     #(4,12)*PSP(1,4)+ACM(5,12)*PSP(1,5)+ACM(6,12)*PSP(1,6)+ACM(7,12)*PS
     #P(1,7)+ACM(8,12)*PSP(1,8)+ACM(9,12)*PSP(1,9))*ACM(1,12)+(ACM(1,12)
     #*PSP(1,2)+ACM(2,12)*PSP(2,2)+ACM(3,12)*PSP(2,3)+ACM(4,12)*PSP(2,4)
     #+ACM(5,12)*PSP(2,5)+ACM(6,12)*PSP(2,6)+ACM(7,12)*PSP(2,7)+ACM(8,12
     #)*PSP(2,8)+ACM(9,12)*PSP(2,9))*ACM(2,12)
      s1 = s2+(ACM(1,12)*PSP(1,3)+ACM(2,12)*PSP(2,3)+ACM(3,12)*PSP(3,3)+
     #ACM(4,12)*PSP(3,4)+ACM(5,12)*PSP(3,5)+ACM(6,12)*PSP(3,6)+ACM(7,12)
     #*PSP(3,7)+ACM(8,12)*PSP(3,8)+ACM(9,12)*PSP(3,9))*ACM(3,12)+(ACM(1,
     #12)*PSP(1,4)+ACM(2,12)*PSP(2,4)+ACM(3,12)*PSP(3,4)+ACM(4,12)*PSP(4
     #,4)+ACM(5,12)*PSP(4,5)+ACM(6,12)*PSP(4,6)+ACM(7,12)*PSP(4,7)+ACM(8
     #,12)*PSP(4,8)+ACM(9,12)*PSP(4,9))*ACM(4,12)
      s2 = s1+(ACM(1,12)*PSP(1,5)+ACM(2,12)*PSP(2,5)+ACM(3,12)*PSP(3,5)+
     #ACM(4,12)*PSP(4,5)+ACM(5,12)*PSP(5,5)+ACM(6,12)*PSP(5,6)+ACM(7,12)
     #*PSP(5,7)+ACM(8,12)*PSP(5,8)+ACM(9,12)*PSP(5,9))*ACM(5,12)+(ACM(1,
     #12)*PSP(1,6)+ACM(2,12)*PSP(2,6)+ACM(3,12)*PSP(3,6)+ACM(4,12)*PSP(4
     #,6)+ACM(5,12)*PSP(5,6)+ACM(6,12)*PSP(6,6)+ACM(7,12)*PSP(6,7)+ACM(8
     #,12)*PSP(6,8)+ACM(9,12)*PSP(6,9))*ACM(6,12)
      ST(78) = s2+(ACM(1,12)*PSP(1,7)+ACM(2,12)*PSP(2,7)+ACM(3,12)*PSP(3
     #,7)+ACM(4,12)*PSP(4,7)+ACM(5,12)*PSP(5,7)+ACM(6,12)*PSP(6,7)+ACM(7
     #,12)*PSP(7,7)+ACM(8,12)*PSP(7,8)+ACM(9,12)*PSP(7,9))*ACM(7,12)+(AC
     #M(1,12)*PSP(1,8)+ACM(2,12)*PSP(2,8)+ACM(3,12)*PSP(3,8)+ACM(4,12)*P
     #SP(4,8)+ACM(5,12)*PSP(5,8)+ACM(6,12)*PSP(6,8)+ACM(7,12)*PSP(7,8)+A
     #CM(8,12)*PSP(8,8)+ACM(9,12)*PSP(8,9))*ACM(8,12)+(ACM(1,12)*PSP(1,9
     #)+ACM(2,12)*PSP(2,9)+ACM(3,12)*PSP(3,9)+ACM(4,12)*PSP(4,9)+ACM(5,1
     #2)*PSP(5,9)+ACM(6,12)*PSP(6,9)+ACM(7,12)*PSP(7,9)+ACM(8,12)*PSP(8,
     #9)+ACM(9,12)*PSP(9,9))*ACM(9,12)


      RETURN
      END
C
C=====================================================================
