C=====================================================================
C	START 4 NODE SHELL ELEMENT ROUTINES
C=====================================================================
C
C=====================================================================
      SUBROUTINE SHELL4(PROPM,PROPG,NODEX,WA,AMV,S,COORD,
	1					EDIS,EDISI,RE,MWG,FIN,MSET)
C	FIN - ADDED TO PREVIOUS LINE BY GILSON - JUL2003 (INT FORCE)

	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	modified by Hari December 2000
C     ----------------------------------------------------------------
C     EVALUATES THE TANGENTIAL STIFFNESS, GAUSS POINT
C     STRAINS, GAUSS POINT STRESS-RESULTANTS AND NODAL STRESS-
C     RESULTANTS FOR 4 TO 8 NODED MINDLIN TYPE ISOPAROMETRIC
C     SHELL ELEMENTS WITH 6 DOF AT EACH NODE
C	--------------------------------------
C     INPUT VARIABLES
C	---------------
C     PROPM(NMP)    = MATERIAL PROPERTIES (YM,PR,YLD,HP,DEN)
C     PROPG(NGP)    = GEOMETRIC PROPERTIES (NNO)
C     NODEX(NEX)    = LOCATIONS OF EXCESS NODES (MIDSIDE NODES)
C     WA(MWG,NPT)   = WORKING ARRAY (8 STRESSES + (8 STRAINS,YLD,IPEL))
C     COORD(3,NNO)  = CURRENT NODAL COORDINATES X,Y,Z
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     EDISI(NEF)    = CURRENT NODAL DISPLACEMENT INCREMENTS
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     S(NWS) 1176   = ELEMENT STIFFNESS MATRIX (UPPER TRIANG.ROW-WISE)
C     RE(NEF)       = EQUILIBRIUM LOADS AT ELEMENT NODES
C	--------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/
C	--------------------------------
C     NAME(2)       = NAME OF ELEMENT MODULE
C     ITYPE         = CODE NUMBER FOR ELEMENT MODULE
C     NLOPT         = CODE FOR NONLINEAR OPTION
C     NLOPT=0         LINEAR ANALYSIS
C     NLOPT=1         MATERIALLY NONLINEAR ONLY
C     NLOPT=2,3       TOTAL LAGRANGIAN,UPDATED LAGRANGIAN
C     MTMOD         = CODE FOR MATERIAL MODULE
C     MTMOD=1         LINEAR ELASTIC,ISOTROPIC
C     MTMOD=2         LINEAR ELASTIC,ORTHOTROPIC
C     MTMOD=3         ELASTO-PLASTIC (IVANOV)
C     MTMOD=4         ELASTO-PLASTIC (MULTI-LAYER)
C     MTMOD=5         CONCRETE WITH CRACKING
C     NSINC         = FACTOR CONTROLLING NUMBER OF SUBINCREMENTS
C     ITOLEY         = TOLERANCE ON YIELD FUNCTION
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
C     GLOC(10,10)   = NATURAL GAUSS POINT COORDINATES (1*1 TO 10*10)
C     GWT (10,10)   = GAUSS POINT WEIGHTS
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
C     COORDI(3,8)   = INITIAL NODAL COORDINATES
C     REDIS(48)     = COROTATIONAL FORM OF EDIS
C     DISD(12)      = DISPLACEMENT DERIVATIVES
C     EPS(8)        = GAUSS POINT STRAINS
C     EPSQ(8)       = QUADRATIC PART OF GAUSS POINT STRAINS
C     SIGR(8)       = GAUSS POINT STRESS-RESULTANTS
C     DR(64)        = ELASTO-PLASTIC RIGIDITY MATRIX STORED COLUMN-WISE
C     PR            = POISSON'S RATIO
C     TH            = GAUSS POINT THICKNESS
C     RN,SN         = NATURAL (NON-DIMENSIONAL) COORDINATES
C     SLR,SLS       = SHEAR LOCKING CONSTRAINT FACTORS ALONG RN,SN
C     NGR,NGS       = NUMBER OF GAUSS PONTS ALONG RN,SN
C     RWT(4),SWT(4) = WEIGHTING FACTORS FOR GAUSS POINTS ALONG RN,SN
C     H(8)          = SHAPE FUNCTIONS
C     HD(2,8)       = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)   = SHAPE FUNCTION DERIVATIVES W.R.T R,S
C     XJI(4)        = INVERSE JACOBIAN STORED COLUMN-WISE
C     DET           = JACOBIAN DETERMINANT
C     DVOL          = INTEGRATION FACTOR (DR*DS=DVOL*DRN*DSN)
C     VT(3)         = DIRECTION COSINE VECTOR ALONG OUTWARD NORMAL TO
C                     LOCAL RN/SN SURFACE (TRUE GEOMETRIC NORMAL)
C     VR(3),VS(3)   = DIRECTION COSINE VECTORS WITH VR TANGENTIAL TO RN
C                     AND VS NORMAL TO VR/VT PLANE (VR,VS,VT FORM A
C                     RIGHT-HANDED ORTHOGONAL SET)
C     RR,SS         = SQUARED BASE VECTOR LENGTHS
C     SNA           = SIN OF ANGLE SUBTENDED BY BASE VECTORS
C     IPEL          = SECTION PLASTICITY INDICATOR (1=EL,2=EL-PL)
C     ------------------------------------------------------------------
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /DYNA/  CDEN,IMASS
C	NEXT COMMON ADDED BY GILSON - SEPT2002
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
C	ADDED BY DE SILVA
	COMMON /HARD/  HP,DEN

C
      DIMENSION PROPM(5),PROPG(*),NODEX(*),WA(MWG,*),S(1176),COORD(3,8)
      DIMENSION EDIS(48),EDISI(48),RE(48),REDIS(48),COORDI(3,8)
      DIMENSION DR(64),H(8),HD(2,8),XJI(4),HR(8),HS(8),B(240)
      DIMENSION BDRL(48),DISD(12),EPS(8),EPSQ(8),SIGR(8)
      DIMENSION VR(3),VS(3),VT(3),SL(4,2),FF(6),NOD(4)
      DIMENSION COVR(3),COVS(3)
	DIMENSION BA(4,120),FJ(4),RRN(4),SSN(4)
C	NEXT LINE ADDED BY GILSON - SEPT2002
	DIMENSION AMV(3)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
C
      EQUIVALENCE (APEL,IPEL)
C
      DATA SL /2.5E5,.1334,.1334,.1334,1.E20,19.2,4.00,4.00/
C
      PI=3.1415926535898
      CALL CLEARA (BDRL,48)
C     ---------------------------------
C     OBTAIN LINEAR STRESS - STRAIN LAW
C     mtmode=1, ipel=1 for linear analysis
	IPEL=1
      If (MTMOD.NE.2) CALL HOKLAW (PROPM,PROPG,2)
C	here 2 indicates plane stress parameters

C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C     ------------------------------------------------------
      CALL SHINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)
C	this subroutine updated coordinate adding displacement for nonlinear 
C	analysis Loop over sampling points, four points for transverse strains
      RRN(1)=1.0
	SSN(1)=0.0
	RRN(2)=-1.0
	SSN(2)=0.0
	RRN(3)=0.0
	SSN(3)=1.0
	RRN(4)=0.0
	SSN(4)=-1.0
C	make all assumed strain Bs values zero at the start
      DO 5 I=1,4
	FJ(I)=0.00
	XJI(I)=0.0
	DO 5 J=1,120
 5	BA(I,J)=0.00
 	NSP=4
      DO 6 I=1,NSP
	CALL SHAP2D (RRN(I),SSN(I),H,HD,NODEX,NNO) 
C	nodex IS NOT USED here for four node element, 
      CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FJ)
C	unit vectors and jacobians are calculated at sampling points
	CALL SHBMATS (NNO,H,HD,VR,VS,VT,XJI,HR,HS,BA,I,FJ)
 6    CONTINUE
C	BA matrix at 4 nodes are completed
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
      IPT=0
      DO 900  IGR=1,NGR
      RN = GLOC(IGR,NGR)
      DO 900  IGS=1,NGS
      SN = GLOC(IGS,NGS)
      WT = GWT(IGR,NGR)*GWT(IGS,NGS)
      IPT= IPT+1

C     -----------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C     -----------------------------------------------------
      CALL SHAP2D (RN,SN,H,HD,NODEX,NNO)
C     ----------------------------------------------------
C     INVERSE JACOBIAN (XJI) , JACOBIAN coefficient F, 
C	DETERMINANT (DET) AND STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------------------
      CALL SHJACO (NNO,COORD,HD,VR,VS,VT,XJI,DET,RR,SS,SNA,1,FJ)
      DVOL=WT*DET
	
	CALL SHSTORL1 (COORD,PROPG,WA(1,IPT),NWG,VR,VS,VT)	!STORE SHELL LAX FOR STRESS CALCULATION JAN09

      CALL RELFILL('TWGH',DVOL*TH*PROPM(5),1,KEG,2)
      
	GOTO 120      
C      IF (NGP.LE.2) GOTO 120
C     ----------------------------------------
C     INTERPOLATE NODAL THICKNESSES (NGP.GT.2)
C     ----------------------------------------
      TH=0.
      DO 100 I=1,NNO
  100 TH=TH+H(I)*PROPG(I+1)   
  
	GOTO 121
	  
  120 TH=PROPG(2)

121	IF (ITASK.NE.5) GOTO 140

C     ---------------------
C     MASS MATRIX (ITASK=5)
C     ---------------------
	IF (MTMOD.EQ.2) PROPM(5) = CDEN
C	CONMSS ADDED BY SONGSAK FOR RC SHELL
	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6) THEN
	CALL CONMSS(TH,MSET,RHORC)
	CALL SHMASS (S,H,VR,VS,DVOL,RHORC,TH,IMASS,NNO,NEF)
	ELSE
      CALL SHMASS (S,H,VR,VS,DVOL,PROPM(5),TH,IMASS,NNO,NEF)
	ENDIF
      GOTO 900


  140 CALL SHBMAT1 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,RN,SN,BA,NEF)
      CALL SHBMAT2 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B)
C     -------------------------------------------
C     SHEAR LOCKING CONSTRAINT FACTORS (SLR,SLS), 
C	MODIFIED WITH ASSUME STRAIN
C     -------------------------------------------
      SLR=1.
	SLS=1.
C	--------------------------------
C	DETERMINE INITIAL MATERIAL ANGLE
C	added by gislon - sept2002
C	--------------------------------
      SELECTCASE(MTMOD)
      CASE(2,5,6)
	CALL CANGLE(VR,VS,AMV,ANG)
	ENDSELECT

C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES (IFREF.EQ.0)
C     -----------------------------------------
      IF (MTMOD.NE.2) CALL SHDELA (DR,DVOL,SLR,SLS)
      IF (NLOPT+ITASK.EQ.1) GOTO 700
C     -------------------------------------------------
C     REMOVE RIGID BODY TRANSLATIONS AND ROTATIONS FROM
C     TOTAL DISPLACEMENT VECTOR (NLOPT>1)
C     -------------------------------------------------
      IF (NLOPT.LE.1) GOTO 180
      CALL SHMDSP (COORD,COORDI,EDIS,REDIS,H,HD,VR,VS,VT,NNO)
C     ------------------------
C     DISPLACEMENT DERIVATIVES
C     ------------------------
  180 CALL CLEARA (DISD,12)
      K=1
      M=1
      DO 200 I=1,NNO
      DO 190 J=1,3
      L=K+3
      DISD(1)=DISD(1)+B(M)*REDIS(K)	
      DISD(2)=DISD(2)+B(M+1)*REDIS(K)
      DISD(3)=DISD(3)+(B(M+2)-HR(I)*VS(J))*REDIS(K)
      DISD(4)=DISD(4)+(B(M+2)-HS(I)*VR(J))*REDIS(K)
      DISD(5)=DISD(5)+B(M+15)*REDIS(L)
      DISD(6)=DISD(6)+B(M+16)*REDIS(L)
      DISD(7)=DISD(7)+B(M+17)*REDIS(L)
      DISD(9)=DISD(9)+B(M+3)*REDIS(K)
      DISD(10)=DISD(10)+B(M+18)*REDIS(L)
      DISD(11)=DISD(11)+B(M+4)*REDIS(K)
      DISD(12)=DISD(12)+B(M+19)*REDIS(L)
      K=K+1
 190  M=M+5
      K=K+3
 200  M=M+15
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      EPS(1) = DISD(1)
      EPS(2) = DISD(2)
      EPS(3) = DISD(3)+DISD(4)
      EPS(4) = DISD(5)
      EPS(5) = DISD(6)
      EPS(6) = DISD(7)
      EPS(7) = DISD(9) +DISD(10)
      EPS(8) = DISD(11)+DISD(12)
C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400
      EPSQ(1)=0.5*(DISD(1)*DISD(1)+DISD(4)*DISD(4)+DISD(9)*DISD(9))
      EPSQ(2)=0.5*(DISD(2)*DISD(2)+DISD(3)*DISD(3)+DISD(11)*DISD(11))
      EPSQ(3)=     DISD(3)*DISD(1)+DISD(2)*DISD(4)+DISD(11)*DISD(9)
      EPSQ(7)=     DISD(10)*DISD(1)+DISD(12)*DISD(4)
      EPSQ(8)=     DISD(10)*DISD(3)+DISD(12)*DISD(2)
      DO 300  I=1,3
 300  EPS(I)= EPS(I)-EPSQ(I)
      EPS(7)= EPS(7)-EPSQ(7)
      EPS(8)= EPS(8)-EPSQ(8)
 400  EPS(7)= EPS(7)
      EPS(8)= EPS(8)
C     ---------------------------------------
C     GAUSS POINT STRESS - RESULTANTS SIGR(8)
C     ---------------------------------------
      GO TO (410,465,450,460,467,468), MTMOD
 410  CALL MPSIGA (EPS,SIGR)

	CALL SHSTRS(VR,VS,VT,EPS,WA(9,IPT),PROPM,TH)

       DO 420 I=1,8
 420  WA(I,IPT)=SIGR(I)
      GO TO 500

c-----COMPOSITE 
 465  CALL COMRGD(EPS,SIGR,DR,PROPM,TH,ANG,DVOL,'SHELL')
      DO 426 I=1,8
	WA(I  ,IPT) =  EPS(I)
 426	WA(I+8,IPT) = SIGR(I)
      GO TO 500

 450	IF (HP.EQ.0.) THEN
	CALL IVANOV (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,SIGR,DR,DVOL)
      ELSE
	EPCA = 0.
	CALL IVANOVH (WA(1,IPT),WA(9,IPT),WA(17,IPT),EPS,
	1					SIGR,DR,DVOL,WA(18,IPT),EPCA)
	END IF
	  APEL = WA(17,IPT)
	  IPEL = 2
      GOTO 500

 460	IF (HP.EQ.0.) THEN
	CALL MLAYER (WA(1,IPT),EPS,SIGR,DR,DVOL)
      IPEL=2
      DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GOTO 500
	ELSE
	CALL MLAYERH (WA(1,IPT),EPS,SIGR,DR,DVOL)
      IPEL=2
      DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GO TO 500
	END IF

 467	CALL RCLAYR (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GO TO 500

 468	CALL RCLAYRE (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510 I=1,8
 510  SIGR(I)=SIGR(I)*DVOL

C	NEXT LINE CHANGED BY GILSON - JUL2003 (INT FORCE)
C      IF (ITASK.LE.2) GOTO 520
      IF (ITASK.LE.3.AND.ISOLOP.NE.4) GOTO 520 

c      CALL SECOND(T1,TIM1)
      IF (IFEIG.EQ.0) GOTO 800

      GOTO 900
 520  K=1
      M=1
      DO 560 I=1,NNO
      DO 540 J=1,3
      L=K+3
      RE(K)=RE(K)+B(M)*SIGR(1)+B(M+1)*SIGR(2)+B(M+2)*SIGR(3)
     1                      +B(M+3)*SIGR(7)+B(M+4)*SIGR(8)
      RE(L)=RE(L)+B(M+15)*SIGR(4)+B(M+16)*SIGR(5)+B(M+17)*SIGR(6)
     1                         +B(M+18)*SIGR(7)+B(M+19)*SIGR(8)
      K=K+1
 540  M=M+5
      K=K+3
 560  M=M+15

C	---------------------------------------
C	FORCES FROM DRILLING DOF SONGSAK JAN09
      IF(NLOPT.EQ.0) THEN
      CALL SHBDRL (NNO,H,HR,HS,VR,VS,VT,BDRL)
	FDRL=10.*DR(19)
	DO IEF = 1,NEF
	DO JEF = 1,NEF
	RE(IEF) = RE(IEF) + BDRL(IEF)*FDRL*BDRL(JEF)*REDIS(JEF)
	ENDDO
	ENDDO
	ENDIF
C	---------------------------------------

C     ------------------------------------------------------
C     FIND LINEAR CONTRIBUTION TO STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------------------
      IF (IFREF ) 900,700,900
c 700  CALL SECOND(T1,TIM1)
C     -----------------------------------------------------------
C     REMOVE SINGULARITY FROM STIFFNESS OF EXACTLY PLANE ELEMENTS
C     -----------------------------------------------------------
 700     IJ=1
      N=3*NEF-2
      LROW=NEF-3
      FAC=DMIN1(DR(19),DR(28)/DET)*1.E-6
      DO 720 I=1,3
      FI=FAC*VT(I)
      DO 720 J=I,3
      FF(IJ)=FI*VT(J)
 720  IJ=IJ+1
      DO 740 I=1,NNO
      HH=H(I)*H(I)
      S(N)=S(N)+HH*FF(1)
      S(N+1)=S(N+1)+HH*FF(2)
      S(N+2)=S(N+2)+HH*FF(3)
      N=N+LROW
      S(N)=S(N)+HH*FF(4)
      S(N+1)=S(N+1)+HH*FF(5)
      N=N+LROW-1
      S(N)=S(N)+HH*FF(6)
      N=N+4*LROW-14
 740  LROW=LROW-6
      CALL SHBDRL (NNO,H,HR,HS,VR,VS,VT,BDRL)
      CALL SHKLIN (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD)
C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 810
C 800  IF (ISTYP.EQ.1) GOTO 805
 800  CALL SHGEO1 (S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)
C      GOTO 810
C 805  CALL SHGEO2 (S,SIGR,HR,HS,NNO)
c 810  CALL SECOND(T1,TIM2)
 810  CONTINUE !TIM(12)=TIM(12)+TIM2
 900  CONTINUE


C	KEF = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	ENDDO
C	ENDDO
C	MEF = 0
C	DO IEF = 1  ,NEF
C	DO JEF = IEF,NEF
C	KEF = KEF + 1
C	MEF = MEF + 1
C	S(KEF) = S(MEF) 
C	ENDDO
C	ENDDO



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
      SUBROUTINE SHINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C	------------------------------------------------------
C     COORD(3,NNO)  = CURRENT NODAL COORDINATES
C     COORDI(3,NNO) = INITIAL NODAL COORDINATES
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     REDIS(48)     = COROTATIONAL FORM OF EDIS
C     NNO           = NUMBER OF NODES
C     NEF           = NUMBER OF ELEMENT FREEDOMS
C     ------------------------------------------------------
C
      DIMENSION COORD(3,8),COORDI(3,8),EDIS(1),REDIS(1)
C

CC      DO 20 I=1,NEF
CC   20 REDIS(I)=EDIS(I)
      REDIS(1:NEF) = EDIS(1:NEF)
      
      K=1
      DO 40 I=1,NNO
      DO 30 J=1,3
      COORDI(J,I)=COORD(J,I)-EDIS(K)
   30 K=K+1
   40 K=K+3
      
      RETURN
      END
C

C=====================================================================
	SUBROUTINE SHBMATS (NNO,H,HD,VR,VS,VT,XJI,HR,HS,BA,ISP,FJ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------
C	THIS SUBROUTINE IS REQUIRED ONLY FOR ASSUMED STRAIN, 4NODE SHELL
c     ----------------------------------------------------------------
c     EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
C	------------------------------------------
c     FJ               = JACOBIAN COEFFICIENT F1,F2,F3,F4 COLUMNWISE
C     isp              = ID OF SAMPLING POINT (1,4)
c     NNO              = NUMBER OF NODES FOR ELEMENT
C     H(4)             = SHAPE FUNCTIONS
C     HD(2,4)          = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(4),HS(4)      = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3)= LOCAL DIRECTION COSINE VECTORS
C     XJI(4)           = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(120)           = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     -----------------------------------------------------------------
c	modified B matrix at gauss point Hari
C	-------------------------------------  
      DIMENSION H(4),HD(2,4),VR(3),VS(3),VT(3),XJI(4),HR(4),HS(4)
      DIMENSION B(120), BA(4,120),FJ(4)
      M=1
	N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      B(M+3)=A1*VT(J)
      B(M+4)=A2*VT(J)
      B(N+3)=A3*VS(J)
      B(N+4)=-A3*VR(J)
      M=M+5
      N=N+5
  10  CONTINUE
      M=M+15
      N=N+15
  20  CONTINUE
C	----------------------------------------------------------------
c	now transfer shear terms from Local TO natural coordinate system
C	----------------------------------------------------------------
	M=4
	DO 30 I=1,4
	DO 30 J=1,6
	BA(ISP,M) = B(M)*FJ(1) + B(M+1)*FJ(3)
	BA(ISP,M+1) = B(M)*FJ(2) + B(M+1)*FJ(4)
	M=M+5
 30	CONTINUE
	RETURN
      END
C
C=====================================================================
      SUBROUTINE SHBMAT1 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B,RN,SN,BA,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	----------------------------------------------------------------
C	B matrix containing transverse shear terms SHBMAT1, Hari
c     ----------------------------------------------------------------
c      EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
c	-------------------------------------------
c     NNO              = NUMBER OF NODES FOR ELEMENT
C     H(4)             = SHAPE FUNCTIONS
C     HD(2,4)          = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(4),HS(4)      = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3)= LOCAL DIRECTION COSINE VECTORS
C     XJI(4)           = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(120)           = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     BA(4,120)        = B VALUES AT SAMPLI NG POINTS, ASSUME STRAIN
C     SN1,SN2,SN3,SN4  = strain interpolation function in natural 
C						coordinate system
C    -----------------------------------------------------------------
      DIMENSION H(4),HD(2,4),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION B(120),BA(4,120)

C	---------------------------------------------------------------
C	COMPUTE INTERPOLATION VALUES AT GAUSS POINTS, 4 SAMPLING POINTS
C	---------------------------------------------------------------
	SN1 = 0.5*(1+RN)
	SN2 = 0.5*(1-RN)
	SN3 = 0.5*(1+SN)
	SN4 = 0.5*(1-SN)
	M = 4
	DO 20 J=1,4
	DO 20 I=1,6
	B(M)=0.0
	B(M+1)=0.0
	B(M)  = BA(3,M)*SN3  +BA(4,M)*SN4
 	B(M+1)= BA(1,M+1)*SN1+BA(2,M+1)*SN2
	M=M+5
  20	CONTINUE
	
C	--------------------------------------------------------------
C	NOW TRANSFER BS FROM NATURAL COORDINATE SYSTEM TO LOCAL SYSTEM
C	NOTE ALL BA(120) VALUES ARE NOT USED
C	------------------------------------
      M=4
	DO 30 I=1,NEF
	B(M)   = XJI(1)*B(M) +XJI(3)*B(M+1)
	B(M+1) = XJI(2)*B(M) +XJI(4)*B(M+1)
  	M=M+5
  30	CONTINUE
C	----------------------------------
C	OTHER TERMS SEE ANOTHER SUBROUTINE
C	----------------------------------
      RETURN
      END 
C
C=====================================================================
      SUBROUTINE SHBMAT2 (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
C	------------------------------------------
C     NNO              = NUMBER OF NODES FOR ELEMENT
C     H(8)             = SHAPE FUNCTIONS
C     HD(2,8)          = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)      = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3)= LOCAL DIRECTION COSINE VECTORS
C     XJI(4)           = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)           = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     ----------------------------------------------------------------
      DIMENSION H(8),HD(2,8),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION B(240)

      M=1
      N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      F1=A2*VR(J)
      F2=A1*VS(J)
      B(M)=A1*VR(J)
      B(M+1)=A2*VS(J)
      B(M+2)=F1+F2
C THIS TEARM IS TRANSVERSE SHEAR INTERPOLATED IN SHBMAT1
C      B(M+3)=A1*VT(J)
C      B(M+4)=A2*VT(J)
      B(N)=F2
      B(N+1)=-F1
      B(N+2)=B(M+1)-B(M)
C THIS TEARM IS TRANSVERSE SHEAR INTERPOLATED IN SHBMAT1
C      B(N+3)=A3*VS(J)
C      B(N+4)=-A3*VR(J)
      M=M+5
      N=N+5
   10 CONTINUE
      M=M+15
      N=N+15
   20 CONTINUE
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHMDSP (COORD,COORDI,EDIS,REDIS,H,HD,VR,VS,VT,NNO)
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
      DIMENSION COORD(3,1),COORDI(3,1),EDIS(1),REDIS(1)
      DIMENSION XYZ(3),XYZO(3),CD(3,9),CDO(3,9)
      DIMENSION TM(3,3),XJI(4),HD(2,1),H(1),VR(3),VS(3),VT(3)
      DIMENSION VV(3),VRO(3),VSO(3),VTO(3),ROT(3),ROV(3),FA(4)
C
      ARC=6.2831853071796
      K=4
      CALL CLEARA (ROT,3)
      DO 50 I=1,NNO
      DO 40 J=1,3
      ROT(J)=ROT(J)+H(I)*EDIS(K)
   40 K=K+1
   50 K=K+3
C     ---------------------------------------------------------------
C     FIND ROTATIONAL PSEUDOVECTOR (ROV) PLUS RESULTANT ROTATION (RO)
C     ---------------------------------------------------------------
      CALL SCALEN (ROT,ROV,RO,3)
      IF (RO.EQ.0.) RETURN
C     -------------------------------------------------------------
C     SET UP CO-ROTATIONAL DISPLACEMENT VECTOR (REDIS) BY DEDUCTING
C     RIGID BODY ROTATIONS FROM EDIS
C     -------------------------------------------------------------
      CALL CLEARA (XYZO,3)
	CALL CLEARA (XYZ ,3)
      DO 140 I=1,NNO
      DO 140 J=1,3
	XYZO(J)=XYZO(J)+H(I)*COORDI(J,I)
  140 XYZ (J)=XYZ (J)+H(I)*COORD (J,I)
      DO 150 I=1,NNO
      DO 150 J=1,3
	CDO(J,I)=COORDI(J,I)-XYZO(J)
  150 CD (J,I)=COORD (J,I)-XYZ (J)
      CALL SHJACO (NNO,COORDI,HD,VRO,VSO,VTO,XJI,DET,RR,SS,SNA,2,FA)
      DO 160 I=1,3
      VV(I)=ROV(I)
      DO 160 J=1,3
 160  TM(I,J)=VR(I)*VRO(J)+VS(I)*VSO(J)+VT(I)*VTO(J)
      ABCD=(.5*(TM(1,1)+TM(2,2)+TM(3,3)-1.))
	IF (ABCD .GT. 1.00000000000000000 )ABCD=1.000000000000000000
	IF (ABCD .LT. -1.00000000000000000 )ABCD=-1.000000000000000000
      RR=ACOS(ABCD)
      SN=SIN(RR)
      IF (SN.EQ.0.) GOTO 190
      F=.5/SN
      VV(1)=F*(TM(3,2)-TM(2,3))
      VV(2)=F*(TM(1,3)-TM(3,1))
      VV(3)=F*(TM(2,1)-TM(1,2))
      CALL SCAPRD (ROV,VV,CS,3)
      IF (CS.GE.0.) GOTO 190
      RR=-RR
      DO 180 I=1,3
  180 VV(I)=-VV(I)
  190 RD=RO-RR
      N=(RD+SIGN(ARC/2.01,RD))/ARC
      RR=RR+N*ARC
      K=1
      DO 220 I=1,NNO
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)
      REDIS(K)=CD(J,I)-TCDO !COORD(J,I)-TCDO
      REDIS(K+3)=EDIS(K+3)-RR*VV(J)
  210 K=K+1
  220 K=K+3
      RETURN
      END
C
C	=====================================================================
      SUBROUTINE MPSIGA (STRAIN,STRESS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------
C     CONVERTS GLOBAL STRAINS TO GLOBAL STRESS-RESULTANTS
C	---------------------------------------------------
C     STRAIN(8) = STRAINS (ER,ES,ERS,KR,KS,KRS,ERT,EST)
C     STRESS(8) = STRESS-RESULTANTS (NR,NS,NRS,MR,MS,MRS,QRT,QST)
C     FOR VARIABLES IN COMMON BLOCK /HOOK/ SEE ROUTINE HOKLAW
C     -----------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION STRAIN(8),STRESS(8)
C
      TH3 = TH*TH*TH/12.
      STRESS(1) = (A1*STRAIN(1) + B1*STRAIN(2)) * TH
      STRESS(2) = (B1*STRAIN(1) + A1*STRAIN(2)) * TH
      STRESS(3) = C1*STRAIN(3)*TH
      STRESS(4) = (A1*STRAIN(4) + B1*STRAIN(5)) * TH3
      STRESS(5) = (B1*STRAIN(4) + A1*STRAIN(5)) * TH3
      STRESS(6) = C1*STRAIN(6)*TH3
      STRESS(7) = D1*STRAIN(7)*TH
      STRESS(8) = D1*STRAIN(8)*TH
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE IVANOV (SIG,EPS,IPEL,STRAIN,STRESS,DP,DVOL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CONVERTS ELASTO-PLASTIC STRAINS INTO STRESSES USING FLOW RULE
C     AND NORMALITY RULE (BASED ON IVANOVS YIELD CRITERION)
C	-----------------------------------------------------
C     INPUT VARIABLES
C	---------------
C     SIG(8)    = STRESS RESULTANTS AT THE END OF THE PREVIOUS UPDATE
C     EPS(8)    = STRAINS AT THE END OF THE PREVIOUS UPDATE
C     IPEL      = PLASTICITY FLAG (1=ELASTIC,2=PLASTIC)
C     STRAIN(8) = CURRENT TOTAL STRAINS
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     STRESS(8) = CURRENT TOTAL STRESS RESULTANTS
C     DP(8,8) = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C	---------------
C     LOCAL VARIABLES
C	---------------
C     DELEPS(8) = INCREMENT IN STRAINS
C     DELSIG(8) = INCREMENT IN STRESS RESULTANTS
C     DEPS(8)   = SUBINCREMENT IN STRAINS (EQUIVALENCED WITH DELEPS)
C     DEPSE(8)  = ELASTIC PART OF STRAIN INCREMENTS
C     DF(6)     = YIELD FUNCTION DERIVATIVES DF/DNX,DF/DNY,DF/DNXY,..
C     RATIO     = PART OF STRESS TAKEN ELASTICLY
C     NSINC     = NUMBER OF SUBINCREMENTS
C     PLAMDA    = PLASTIC STRAIN RATE MULTIPLIER
C     /HOOK/    = FOR VARIABLES IN COMMON BLOCK HOOK SEE CONSTI
C     ----------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
      DIMENSION SIG(8),EPS(8),STRAIN(8),STRESS(8),DP(8,8)
      DIMENSION DELEPS(8),DELSIG(8),DEPS(8),DEPSE(8),DF(6)
      EQUIVALENCE (DELEPS(1),DEPS(1))
C     --------------------------------------------------------------
C     1. CALCULATE INCREMENTAL STRAINS
C     2. CALCULATE INCREMENTAL STRESS RES. ASSUMING ELASTIC BEHAVIOR
C     3. CALCULATE TOTAL STRESS RESULTANTS,ASSUMING ELASTIC BEHAVIOR
C     --------------------------------------------------------------
      DO 110  I=1,8
 110  DELEPS(I) = STRAIN(I)-EPS(I)
      CALL MPSIGA (DELEPS,DELSIG)
      DO 300  I=1,8
 300  STRESS(I) = SIG(I) + DELSIG(I)
C     ----------------------------------------------------------------
C     4. CHECK WHETHER STATE OF STRESS FALLS OUTSIDE THE YIELD SURFACE
C     ----------------------------------------------------------------
      CALL YIELDF (STRESS(1),STRESS(2),STRESS(3),STRESS(4),STRESS(5),
     1             STRESS(6),UYN,UYM,UYN2,UYM2,QT,QM,QTM,Q,R,S,FT,2)
      IF (FT)  410,410,450
C     -------------------------------------------------------
C     STATE OF STRESS WITHIN YIELD SURFACE - ELASTIC BEHAVIOR
C     -------------------------------------------------------
 410  IPEL = 1
      GOTO 700
C     ------------------------------------------------------------
C     STATE OF STRESS OUTSIDE YIELD SURFACE - PLASTIC BEHAVIOR
C     DETERMINE PART OF STRAIN TAKEN ELASTICLY AND ADD STRESSES
C     DUE TO THE ELASTIC STRAIN INCREMENT TO THE PREVIOUS STRESSES
C     ------------------------------------------------------------
 450  IPEL = 2
      CALL YIELDF (SIG(1),SIG(2),SIG(3),SIG(4),SIG(5),SIG(6),
     1             UYN,UYM,UYN2,UYM2,QT,QM,QTM,Q,R,S,YFT,2)
      IF (YFT.LT.0)  GOTO 470
      RATIO = 0.
      DO 460  I=1,8
 460  STRESS(I) = SIG(I)
      GOTO 500
 470  CALL TRANSI (SIG(1),SIG(2),SIG(3),SIG(4),SIG(5),SIG(6),
     1             DELSIG(1),DELSIG(2),DELSIG(3),DELSIG(4),DELSIG(5),
     2             DELSIG(6),UYN,UYM,UYN2,UYM2,RATIO)
      DO 480  I=1,8
 480  STRESS(I) = SIG(I) + RATIO*DELSIG(I)
C     --------------------------------------------------
C     5. DETERMINE NUMBER AND MAGNITUDE OF SUBINCREMENTS
C     --------------------------------------------------
 500  NSINC = 20.*SQRT(FT) + 1
      IF (NSINC.GT.30)  NSINC = 30
      FACT = (1.0-RATIO) / NSINC
      DO 510  I=1,8
 510  DEPS(I) = FACT*DELEPS(I)
C     ----------------------------------------
C     6. CALCULATION OF ELASTOPLASTIC STRESSES
C     ----------------------------------------
      DO 690  ISINC=1,NSINC
      CALL LAMDAP (STRESS(1),STRESS(2),STRESS(3),STRESS(4),STRESS(5),
     1             STRESS(6),DEPS,PLAMDA,DF,DENOM)
      IF (PLAMDA.LT.0)  PLAMDA = 0.
C
      DO 610  I=1,6
 610  DEPSE(I) = DEPS(I) - PLAMDA*DF(I)
      DEPSE(7) = DEPS(7)
      DEPSE(8) = DEPS(8)
      CALL MPSIGA (DEPSE,DELSIG)
      DO 620  I=1,8
 620  STRESS(I) = STRESS(I) + DELSIG(I)
C     FORCE STRESS STATE BACK ON YIELD SURFACE
C
 690  CONTINUE
C     -----------------------------------------
C     7. UPDATING STRESS RESULTANTS AND STRAINS
C     -----------------------------------------
 700  DO 710  I=1,8
      SIG(I) = STRESS(I)
 710  EPS(I) = STRAIN(I)
C     -------------------------------------------
C     8. FORM THE MATERIAL LAW (MATRIX DP) IF THE
C        STIFNESS IS TO BE REFORMED (IFREF=0)
C     -------------------------------------------
      IF (IFREF.NE.0) RETURN
      IF (IPEL.EQ.1) RETURN
      CALL LAMDAP (STRESS(1),STRESS(2),STRESS(3),STRESS(4),STRESS(5),
     1             STRESS(6),DEPS,PLAMDA,DF,DENOM)
      IF (PLAMDA.GT.0.) CALL MPDELP (DF,DVOL/DENOM,DP,DP)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE YIELDF (NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     1                   QT,QM,QTM,Q,R,S,FT,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     RETURNS CURRENT VALUE OF IVANOVS YIELD FUNCTION
C	-----------------------------------------------
C     NX,NY,NXY,MX,MY,MXY= STRESS RESULTANTS
C     UYN,UYM,UYN2,UYM2  = UNIAXIAL YIELD FORCE/MOMENT PER UNIT WIDTH
C     QT,QM,QTM          = QUADRATIC NON-DIMENSIONAL STRESS INTENSITIE
C     Q,R,S              = COMPONENTS IN IVANOVS YIELD FUNCTION
C     FT                 = VALUE OF IVANOV'S YIELD FUNCTION
C     IND                = INDICATES WHETHER FT IS REQUIRED (IND=2)
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      REAL*8 NX,NY,NXY,MX,MY,MXY
C
      UYN  = YLD*TH
      UYM  = YLD*TH*TH/4.
      UYN2 = UYN*UYN
      UYM2 = UYM*UYM
C
      QT =(NX*NX + NY*NY - NX*NY + 3.*NXY*NXY) / UYN2
      QM =(MX*MX + MY*MY - MX*MY + 3.*MXY*MXY) / UYM2
      QTM=(MX*NX + MY*NY - MX*NY/2. - MY*NX/2. + 3.*MXY*NXY)/(UYN*UYM)
C
      ENTRY FTVALU(NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     +             QT,QM,QTM,Q,R,S,FT,IND)
      Q   = QT + .48*QM
      R   = SQRT(QM*QM/4. + QTM*QTM)
      S   = QT*QM - QTM*QTM
C
      IF (IND.EQ.1)  RETURN
      FT = -1.
      IF (Q.EQ.0.)  RETURN
      FT = QT + QM/2. + R - .25*S/Q - 1.
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE TRANSI (NX,NY,NXY,MX,MY,MXY,DNX,DNY,DNXY,
     1                   DMX,DMY,DMXY,UYN,UYM,UYN2,UYM2,RATIO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     RETURNS RATIO OF STRAIN INCREMENTS TAKEN ELASTICLY
C	--------------------------------------------------
C     NX,NY,NXY,MX,MY,MXY      = STRESS RESULTANTS AT PREVIOUS UPDATE
C     DNX,DNY,DNXY,DMX,DMY,DMXY= INCREMENTS IN STRESS RESULTANTS
C     UYN,UYM,UYN2,UYM2        = UNIAXIAL YIELD FORCE/UNIT WIDTH
C     RATIO                    = PART OF STRAIN INCR. TAKEN ELASTICLY
C     ----------------------------------------------------------------
	
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
C
      REAL*8 NX,NY,NXY,MX,MY,MXY
C     -------------------------------------
C     QUADRATIC IN PLANE STRESS INTENSITIES
C     N = AM*RATIO**2 + BM*RATIO + CM
C     -------------------------------------
	IF (ITOLEY .EQ.0) ITOLEY = 100000
      AM = DNX*DNX + DNY*DNY - DNX*DNY + 3.*DNXY*DNXY
      BM = 2.*(NX*DNX + NY*DNY) - NX*DNY - DNX*NY + 6.*NXY*DNXY
      CM = NX*NX + NY*NY - NX*NY + 3.*NXY*NXY
C     ------------------------------------
C     QUADRATIC BENDING STRESS INTENSITIES
C     M = AB*RATIO**2 + BB*RATIO + CB
C     ------------------------------------
      AB = DMX*DMX + DMY*DMY - DMX*DMY + 3.*DMXY*DMXY
      BB = 2.*(MX*DMX + MY*DMY) - MX*DMY - DMX*MY + 6.*MXY*DMXY
      CB = MX*MX + MY*MY - MX*MY + 3.*MXY*MXY
C     --------------------------------------------------------
C     COMPUTE RATIO IN CASE OF PURE STRETCHING OR PURE BENDING
C     --------------------------------------------------------
      IF (AM+CM.NE.0. .AND. AB+CB.NE.0.)  GOTO 200
      RATIO = -BM-BB+SQRT(BM*BM-4.*AM*(CM-UYN2)+BB*BB-4.*AB*(CB-UYM2))
      RATIO = RATIO/2./(AM+AB)
      RETURN
C     -----------------------------------
C     MIXED QUADRATIC STRESS INTENSITIES
C     MN = AMB*RATIO**2 + BMB*RATIO + CMB
C     -----------------------------------
 200  AMB = DMX*DNX + DMY*DNY - .5*(DMX*DNY+DMY*DNX) + 3.*DMXY*DNXY
      BMB = MX*DNX + DMX*NX + MY*DNY + DMY*NY - .5*(MX*DNY+DMX*NY+
     1      MY*DNX+DMY*NX) + 3.*(MXY*DNXY+DMXY*NXY)
      CMB = MX*NX + MY*NY - .5*(MX*NY+MY*NX) + 3.*MXY*NXY
C     ---------------------------------------
C     NUMERICAL LOOP TO FIND RATIO FOR FT = 0
C     ---------------------------------------
      GAP = 1.0
      FAC = 1.0
      RATIO = 0.0
C
 300  GAP    = GAP/2.
      RATIO  = RATIO + FAC*GAP
      RATIO2 = RATIO*RATIO
      FAC    = 1.
C
      QT  = (AM *RATIO2 + BM *RATIO + CM) / UYN2
      QM  = (AB *RATIO2 + BB *RATIO + CB) / UYM2
      QTM = (AMB*RATIO2 + BMB*RATIO + CMB) / (UYN*UYM)
      CALL FTVALU (NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     1             QT,QM,QTM,Q,R,S,FT,2)
C
      IF (FT.GT.0.)  FAC = -1.
C	ABCD=1./(DFLOAT(ITOLEY))
	ABCD=1./1000000.
      IF (ABS(FT)-ABCD)  400,400,300
C
  400 RETURN
      END
C
C=====================================================================
      SUBROUTINE LAMDAP (NX,NY,NXY,MX,MY,MXY,DEPS,PLAMDA,DF,B)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     RETURNS THE PLASTIC STRAIN RATE MULTIPLIER AND THE YIELD
C     FUNCTION DERIVATIVES
C	--------------------
C     NX,NY,NXY,MX,MY,MXY  = STRESS RESULTANTS
C     DEPS(8)              = SUBINCREMENT OF STRAIN
C     DF(6)                = YIELD FUNCTION DERIVATIVES
C     PLAMDA               = PLASTIC STRAIN RATE MULTIPLIER
C     B                    = DENOMINATOR OF PLAMDA
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION DEPS(8),DQ(6),DF(6),DFD(6)
C
      REAL*8 NX,NY,NXY,MX,MY,MXY
C     -----------------------------------------------------
C     UNIAXIAL YIELD FORCE AND MOMENT, QUADRATIC STRESS
C     INTENSITIES N,M,NM, CONSTANTS Q,R,S OF YIELD FUNCTION
C     -----------------------------------------------------
      CALL YIELDF (NX,NY,NXY,MX,MY,MXY,UYN,UYM,UYN2,UYM2,
     1             QT,QM,QTM,Q,R,S,FT,1)
C     -----------------------------------------------------
C     CHECK FOR DISCONTINUITY OF YIELD FUNCTION DERIVATIVES
C     -----------------------------------------------------
      HR = 0.
      IF (R.GT.1.E-04)  HR = 1./R
C
      Q2 = Q*Q
      Q4 = 4.*Q
      C = 1. - QM/Q4 + S/(4.*Q2)
      D = .5 + HR*QM/4. - QT/Q4 + .12*S/Q2
      UYNM = UYM/UYN
      E = .5*QTM * (HR+.5/Q)
      F = E/UYNM
      E = E*UYNM
C     -------------------------------------------
C     VECTOR OF DERIVATIVES OF STRESS INTENSITIES
C     -------------------------------------------
      DQ(1) = (2.*NX-NY) / UYN2
      DQ(2) = (2.*NY-NX) / UYN2
      DQ(3) = 6.*NXY/UYN2
      DQ(4) = (2.*MX-MY) / UYM2
      DQ(5) = (2.*MY-MX) / UYM2
      DQ(6) = 6.*MXY/UYM2
C     ------------------------------------
C     VECTOR OF YIELD FUNCTION DERIVATIVES
C     ------------------------------------
      DO 100  I=1,3
      DF(I)   = C*DQ(I) + E*DQ(I+3)
 100  DF(I+3) = F*DQ(I) + D*DQ(I+3)
C     ------------------------------
C     PLASTIC STRAIN RATE MULTIPLIER
C     ------------------------------
      TH3 = TH*TH*TH/12.
      DFD(1) = (A1*DF(1) + B1*DF(2)) * TH
      DFD(2) = (B1*DF(1) + A1*DF(2)) * TH
      DFD(3) = C1*DF(3)*TH
      DFD(4) = (A1*DF(4) + B1*DF(5)) * TH3
      DFD(5) = (B1*DF(4) + A1*DF(5)) * TH3
      DFD(6) = C1*DF(6)*TH3
C
      A = 0.
      B = 0.
      DO 200  I=1,6
      A = A + DFD(I)*DEPS(I)
 200  B = B + DFD(I)*DF(I)
      PLAMDA = A/B
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MPDELP (DF,FAC,DP,CP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     FORMS THE ELASTO-PLASTIC STRESS - STRAIN MATRIX
C	-----------------------------------------------
C     DF(6)    = YIELD FUNCTION DERIVATIVES
C     FAC      = DVOL/DENOMINATOR OF STRAIN RATE MULTIPLIER
C     DP(64)   = STRESS - STRAIN MATRIX
C     CP(8,8)  = DP(64)
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION DF(6),DP(64),CP(8,8),FFT(6,6),FD(6)
C     -----------------------------------------------
C     1. PRODUCT MATRIX OF YIELD FUNCTION DERIVATIVES
C        FFT(I,J) = DF(I)*DF(J)
C     -----------------------------------------------
      DO 50  I=1,6
      DO 50  J=I,6
      FFT(I,J) = DF(I)*DF(J)
 50   FFT(J,I) = FFT(I,J)
C     ---------------------------------------------------------
C     2. MATRIX MULTIPLICATION [DP] = [DP] - [D]*[FFT]*[D]*FACT
C     ---------------------------------------------------------
      TH3   = TH*TH*TH/12.
      FACM  = TH*TH*FAC
      FACB  = TH3*TH3*FAC
      FACMB = TH*TH3*FAC
C
      DO 100  I=1,6
 100  FD(I)  = FFT(I,1)*A1 + FFT(I,2)*B1
      DP(1)  = DP(1)  - (A1*FD(1) + B1*FD(2)) * FACM
      DP(2)  = DP(2)  - (B1*FD(1) + A1*FD(2)) * FACM
      DP(3)  =        -  C1*FD(3)*FACM
      DP(4)  =        - (A1*FD(4) + B1*FD(5)) * FACMB
      DP(5)  =        - (B1*FD(4) + A1*FD(5)) * FACMB
      DP(6)  =        -  C1*FD(6)*FACMB
C
      DO 200  I=1,6
 200  FD(I)  = FFT(I,1)*B1 + FFT(I,2)*A1
      DP(10) = DP(10) - (B1*FD(1) + A1*FD(2)) * FACM
      DP(11) =        -  C1*FD(3)*FACM
      DP(12) =        - (A1*FD(4) + B1*FD(5)) * FACMB
      DP(13) =        - (B1*FD(4) + A1*FD(5)) * FACMB
      DP(14) =        -  C1*FD(6)*FACMB
C
      DO 300  I=3,6
 300  FD(I)  = FFT(I,3)*C1
      DP(19) = DP(19) -  C1*FD(3)*FACM
      DP(20) =        - (A1*FD(4) + B1*FD(5)) * FACMB
      DP(21) =        - (B1*FD(4) + A1*FD(5)) * FACMB
      DP(22) =        -  C1*FD(6)*FACMB
C
      DO 400  I=4,6
 400  FD(I)  = FFT(I,4)*A1 + FFT(I,5)*B1
      DP(28) = DP(28) - (A1*FD(4) + B1*FD(5)) * FACB
      DP(29) = DP(29) - (B1*FD(4) + A1*FD(5)) * FACB
      DP(30) =        -  C1*FD(6)*FACB
C
      DO 500  I=4,6
 500  FD(I)  = FFT(I,4)*B1 + FFT(I,5)*A1
      DP(37) = DP(37) - (B1*FD(4) + A1*FD(5)) * FACB
      DP(38) =        -  C1*FD(6)*FACB
C
      FD(6)  = FFT(6,6)*C1
      DP(46) = DP(46) -  C1*FD(6)*FACB
C     -----------------------------------------------
C     3. FILL IN SYMMETRIC ELEMENTS OF UPPER TRIANGLE
C     -----------------------------------------------
 900  DO 950  I=1,6
      DO 950  J=I,6
 950  CP(I,J) = CP(J,I)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MLAYER (WA,STRAIN,STRESS,DP,DVOL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	integer apel
C     --------------------------------------------------------------
C     SET UP LOOP OVER NUMBER OF LAYERS
C     CALL VMISES TO COMPUTE ELASTO-PLASTIC STRESSES AND STRESS RES.
C     INTEGRATES ELESTO-PLASTIC MATERIAL MATRIX DP OVER THE DEPTH
C	-----------------------------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     WA(8,NGT+1) = WORKING ARRAY STORING STRESSES AND STRAINS
C     STRAIN(8)   = CURRENT TOTAL STRAINS
C     STRESS(8)   = RETURNED STRESS RESULTANTS
C     DP(8,8)     = MATRIX CONTAINING ELASTO-PLASTIC RIGIDITIES
C     DVOL        = WT*DET
C	---------------
C     LOCAL VARIABLES
C	---------------
C     EPSZ(4)     = STRAINS AT LEVEL Z
C     SIGZ(4)     = STRESSES AT LEVEL Z
C     CP(6,6)     = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     DEPTH(26)   = Z-COORDINATES OF LAYER BOUNDARIES
C     --------------------------------------------------------------
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
      DIMENSION WA(8,1),STRAIN(8),STRESS(8),DP(8,8)
      DIMENSION EPSZ(4),SIGZ(4),CP(6,6),DEPTH(26)
      EQUIVALENCE (IPE,APEL)
C
      IND = 3
      ISR = 3
      IST = 3
      IPEL = 1
      CALL CLEARA (STRESS,8)
      IF (IFREF.EQ.0) CALL CLEARA (DP,64)
      IF (NGT.GT.1)  GOTO 400
C     --------------------------
C     PURE MEMBRANE STRESS STATE
C     --------------------------
      DO 110  I=1,3
 110  EPSZ(I) = STRAIN(I)
c      CALL VMISES (WA(1,1),WA(4,1),WA(8,1),EPSZ,SIGZ,CP,IND)
      CALL VMISES (WA(1,1),WA(4,1),apel,EPSZ,SIGZ,CP,IND)
      WA(8,IGT) = APEL !added by gilson - aug2004
      DO 120  I=1,3
      STRESS(I) = SIGZ(I)*TH
 120  WA(8+I,1) = STRESS(I)
      IF (IFREF.NE.0)  RETURN
C     -------------------------
C     ELASTO-PLASTIC RIGIDITIES
C     -------------------------
      TH3 = TH*TH*TH/12.
      FACM = DVOL*TH
      FACB = DVOL*TH3
      DO 150  I=1,3
      DO 150  J=I,3
      K = I+3
      L = J+3
      DP(I,J) = CP(I,J)*FACM
      DP(K,L) = CP(I,J)*FACB
      DP(J,I) = DP(I,J)
 150  DP(L,K) = DP(K,L)
      DP(7,7) = D1*FACM
      DP(8,8) = D1*FACM
      RETURN
C     ------------------------------------
C	MULTI-LAYER SOLUTION, INITIALISATION
C     ------------------------------------
 400  DZ = TH/(NGT-1)
      Z  = -TH/2-DZ
C     -------------------------------------------------------
C     LOOP OVER NUMBER OF LAYERS,CALCULATE STRAINS AT LEVEL Z
C     -------------------------------------------------------
      DO 800  IGT=1,NGT
      Z = Z+DZ
      DO 550  I=1,3
 550  EPSZ(I) = STRAIN(I) + Z*STRAIN(I+3)
C     ----------------------
C     COMPUTE LAYER STRESSES
C     ----------------------
c      CALL VMISES (WA(1,IGT),WA(4,IGT),WA(8,IGT),EPSZ,SIGZ,CP,IND)
      CALL VMISES (WA(1,IGT),WA(4,IGT),apel,EPSZ,SIGZ,CP,IND)
      WA(8,IGT) = APEL
      IF (IPE.EQ.2)  IPEL = 2
C     -------------------------------------------
C     ADD LAYER CONTRIBUTION TO STRESS RESULTANTS
C     -------------------------------------------
      WTMEM=DZ
      ZZ = Z
      IF (IGT.NE.1)  GOTO 620
      WTMEM = DZ/2.
      ZZ = Z+DZ/3.
 620  IF (IGT.NE.NGT)  GOTO 650
      WTMEM = DZ/2.
      ZZ = Z-DZ/3.
 650  WTBEN = ZZ*WTMEM
      DO 690  I=1,3
      J = I+3
      STRESS(I) = STRESS(I) + WTMEM*SIGZ(I)
 690  STRESS(J) = STRESS(J) + WTBEN*SIGZ(I)
      IF (IFREF.NE.0)  GOTO 800
C     -------------------------
C     INTEGRATE MATERIAL MATRIX
C     -------------------------
      FACM = WTMEM*DVOL
      FACC = Z*FACM
      FACB = Z*FACC
      DO 750  I=1,3
      DO 750  J=I,3
      K = I+3
      L = J+3
      DP(I,J) = DP(I,J) + CP(I,J)*FACM
 750  DP(K,L) = DP(K,L) + CP(I,J)*FACB
      DO 760  I=1,3
      DO 760  J=1,3
      L = J+3
 760  DP(I,L) = DP(I,L) + CP(I,J)*FACC
C
 800  CONTINUE
C     ---------------------------------------------------
C     FILL IN LOWER TRIANGLE OF RIGIDITY MATRIX (IFREF=0)
C     ---------------------------------------------------
      IF (IFREF.NE.0)  GOTO 900
      IF (IPEL.EQ.1)   GOTO 850
      DO 830  I=1,6
      DO 830  J=I,6
 830  DP(J,I) = DP(I,J)
      DP(7,7) = D1*TH*DVOL
      DP(8,8) = D1*TH*DVOL
      GOTO 900
C
 850  	SBB=1.
	SAA=1.
	CALL SHDELA (DP,DVOL,SBB,SAA)
C     ------------------------------------
C     STORE STRESS RESULTANTS AT END OF WA
C     ------------------------------------
 900  STRESS(7) = D1*STRAIN(7)*TH
      STRESS(8) = D1*STRAIN(8)*TH
      DO 990  I=1,8
 990  WA(I,NGT+1) = STRESS(I)
C
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHDELA (DR,DVOL,SLA,SLB)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATES ELASTIC RIGIDITIES FOR MEMBRANE-PLATE OR SHELL ELEM.
C	---------------------------------------------------------------
C     DR(8,8) = ELASTIC MEMBRANE,BENDING AND SHEAR RIGIDITIES
C     DVOL    = WT*DET
C     SLA,SLB = SHEAR LOCKING CONSTRAINT FACTORS
C	----------------------------------------------------------------
C	PLEASE NOT SLA,SLB ARE 1.0 MODIFIED ASSUMED STRAIN METHOD, 
C	NO FACTOR REQUIRED
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION DR(64)
C
      CALL CLEARA (DR,64)

C	WRITE(*,*) PR
C	PAUSE

      TH3  = TH*TH*TH/12.
      FACM = DVOL*TH
      FACB = DVOL*TH3
C
      DR(1)  = FACM*A1
      DR(10) = FACM*A1
      DR(2)  = FACM*B1
      DR(9)  = FACM*B1
      DR(19) = FACM*C1
      DR(28) = FACB*A1
      DR(37) = FACB*A1
      DR(29) = FACB*B1
      DR(36) = FACB*B1
      DR(46) = FACB*C1
      DR(55) = FACM*SLA*D1
      DR(64) = FACM*SLB*D1

C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHBDRL (NNO,H,HR,HS,VR,VS,VT,BDRL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     EVALUATES DERIVATIVE OPERATORS ASSOCIATED WITH THE DRILLING
C     CONSTRAINT THETA(T)=0.5*(DV/DR-DU/DS) (PENALTY FUNCTION METHOD)
C	---------------------------------------------------------------
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     H(8)              = SHAPE FUNCTIONS
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     BDRL(48)          = DERIVATIVE OPERATORS FOR DRILLING STRAIN
C     ---------------------------------------------------------------
      DIMENSION H(8),HR(8),HS(8),VR(3),VS(3),VT(3),BDRL(48)
C
      N=1
      DO 20 I=1,NNO
      A1=.5*HR(I)
      A2=.5*HS(I)
      A3=H(I)
      DO 10 J=1,3
      BDRL(N)=A2*VR(J)-A1*VS(J)
      BDRL(N+3)=A3*VT(J)
   10 N=N+1
   20 N=N+3
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHKLIN (S,DR,B,BDRL,NNO,NEF,IPEL,MTMOD)
C
C     ------------------------------------------------------------
C     EVALUATES LINEAR CONTRIBUTION TO TANGENTIAL STIFFNESS MATRIX
C
C     S(1176)  = STIFFNESS MATRIX STORED UPPER-TRIANGULAR ROW-WISE
C     DR(64)   = ELASTO-PLASTIC RIGIDITY MATRIX STORED COLUMN-WISE
C     B(240)   = STRAIN-DISPLACEMENT MATRIX
C     BDRL(48) = DERIVATIVE OPERATORS FOR DRILLING STRAIN
C     NNO      = NUMBER OF NODES FOR ELEMENT
C     NEF      = NUMBER OF D.O.F FOR ELEMENT
C     IPEL     = SECTION PLASTICITY INDICATOR (1=EL,2=EL-PL)
C     ------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(*),DR(64),B(240),BDRL(48)
C
      N=1
      I=4
      FDRL=10.*DR(19)
C     -----------------------------
C     TRANSVERSE SHEAR CONTRIBUTION
C     -----------------------------
      DO 30 NROW=1,NEF
      A1=DR(55)*B(I) + DR(56)*B(I+1)
      A2=DR(63)*B(I) + DR(64)*B(I+1)
      A3=FDRL*BDRL(NROW)
      J=I
      DO 20 NCOL=NROW,NEF
      S(N)=S(N)+A1*B(J)+A2*B(J+1)+A3*BDRL(NCOL)
      N=N+1
   20 J=J+5
   30 I=I+5
      LROW=NEF
      I=1
      N1=1
      N2=3*NEF-2
      IF (IPEL.EQ.2.OR.MTMOD.EQ.2) GO TO 120
C     -----------------------
C     ELASTIC RIGIDITY MATRIX
C     -----------------------
      DO 100 IR=1,NNO
      DO 90 IRR=1,3
      A1=DR(1)*B(I)+DR(2)*B(I+1)
      A2=DR(2)*B(I)+DR(10)*B(I+1)
      A3=DR(19)*B(I+2)
      A4=DR(28)*B(I+15)+DR(29)*B(I+16)
      A5=DR(29)*B(I+15)+DR(37)*B(I+16)
      A6=DR(46)*B(I+17)
      J=I
      IF=4-IRR
C     ------------------------------------------------
C     UPPER TRIANGULAR PART OF DIAGONAL 3*3 PARTITIONS
C     ------------------------------------------------
      DO 40 JR=1,IF
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
   40 J=J+5
C     ----------------------------------------
C     INTERMEDIATE OFF-DIAGONAL 3*3 PARTITIONS
C     ----------------------------------------
      N1=N1+3
      IF (IR.EQ.NNO) GO TO 70
      NB=NNO-IR
      N2=N2+3
      J=J+15
      DO 60 JB=1,NB
      DO 50 JR=1,3
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
   50 J=J+5
      N1=N1+3
      N2=N2+3
   60 J=J+15
      N2=N2-3
   70 I=I+5
   90 LROW=LROW-1
      I=I+15
      N1=N2
      LROW=LROW-3
  100 N2=N2+3*LROW-3
      RETURN
C     --------------------------------------
C     MULTILAYERED ANISOTROPIC COMPOSITE AND
C     ELASTO-PLASTIC RIGIDITY MATRIX
C     --------------------------------------
  120 DO 200 IR=1,NNO
      DO 190 IRR=1,3
      A1=DR(1)*B(I)+DR(2)*B(I+1)+DR(3)*B(I+2)
      A2=DR(2)*B(I)+DR(10)*B(I+1)+DR(11)*B(I+2)
      A3=DR(3)*B(I)+DR(11)*B(I+1)+DR(19)*B(I+2)
      A4=DR(28)*B(I+15)+DR(29)*B(I+16)+DR(30)*B(I+17)
      A5=DR(29)*B(I+15)+DR(37)*B(I+16)+DR(38)*B(I+17)
      A6=DR(30)*B(I+15)+DR(38)*B(I+16)+DR(46)*B(I+17)
      A7=DR(4)*B(I)+DR(12)*B(I+1)+DR(20)*B(I+2)
      A8=DR(5)*B(I)+DR(13)*B(I+1)+DR(21)*B(I+2)
      A9=DR(6)*B(I)+DR(14)*B(I+1)+DR(22)*B(I+2)
      J=I
      IF=4-IRR
C     ------------------------------------------------
C     UPPER TRIANGULAR PART OF DIAGONAL 6*6 PARTITIONS
C     ------------------------------------------------
      DO 140 JR=1,IF
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      N1=N1+1
      N2=N2+1
  140 J=J+5
      S(N1)=S(N1)+A7*B(J)+A8*B(J+1)+A9*B(J+2)
      S(N1+1)=S(N1+1)+A7*B(J+5)+A8*B(J+6)+A9*B(J+7)
      S(N1+2)=S(N1+2)+A7*B(J+10)+A8*B(J+11)+A9*B(J+12)
C     ---------------------------
C     OFF-DIAGONAL 6*6 PARTITIONS
C     ---------------------------
      N1=N1+3
      IF (IR.EQ.NNO) GO TO 170
      NB=NNO-IR
      N2=N2+3
      J=J+15
      A10=DR( 4)*B(I+15)+DR( 5)*B(I+16)+DR( 6)*B(I+17)
      A11=DR(12)*B(I+15)+DR(13)*B(I+16)+DR(14)*B(I+17)
      A12=DR(20)*B(I+15)+DR(21)*B(I+16)+DR(22)*B(I+17)
      DO 160 JB=1,NB
      DO 150 JR=1,3
      S(N1)=S(N1)+A1*B(J)+A2*B(J+1)+A3*B(J+2)
      S(N2)=S(N2)+A4*B(J+15)+A5*B(J+16)+A6*B(J+17)
      S(N1+3)=S(N1+3)+A7*B(J+15)+A8*B(J+16)+A9*B(J+17)
      S(N2-3)=S(N2-3)+A10*B(J)+A11*B(J+1)+A12*B(J+2)
      N1=N1+1
      N2=N2+1
  150 J=J+5
      N1=N1+3
      N2=N2+3
  160 J=J+15
      N2=N2-3
  170 I=I+5
  190 LROW=LROW-1
      I=I+15
      N1=N2
      LROW=LROW-3
  200 N2=N2+3*LROW-3
      RETURN
      END
C=====================================================================
      SUBROUTINE SHGEO1 (S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     EVALUATES GEOMETRIC CONTRIBUTION TO TANGENTIAL STIFFNESS MATRIX
C     THIS VERSION IS BASED ON FULL GREEN'S STRAIN EXPANSION
C	------------------------------------------------------
C     S(1176)           = GEOMETRIC STIFFNESS STORED UPPER TRIANGULAR
C                         ROW-WISE
C     SIGR(8)           = WEIGHTED STRESS-RESULTANT VECTOR
C     H(8)              = SHAPE FUNCTIONS
C     HR(8),HS(8)       = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3) = LOCAL DIRECTION COSINE VECTORS
C     TH                = GAUSS POINT THICKNESS
C     NNO               = NUMBER OF NODES FOR ELEMENT
C     ---------------------------------------------------------------
      DIMENSION S(*),SIGR(8),H(8),HR(8),HS(8),VR(3),VS(3),VT(3),F(6)
C
      F1=VR(1)*VR(1)+VS(1)*VS(1)
      F2=VR(1)*VR(2)+VS(1)*VS(2)
      F3=VR(1)*VR(3)+VS(1)*VS(3)
      F4=VR(2)*VR(2)+VS(2)*VS(2)
      F5=VR(2)*VR(3)+VS(2)*VS(3)
      F6=VR(3)*VR(3)+VS(3)*VS(3)
      FB=TH*TH/12.
      NS=1
      LROW=6*NNO
      DO 30 I=1,NNO
      A1=SIGR(1)*HR(I)+SIGR(3)*HS(I)
      A2=SIGR(3)*HR(I)+SIGR(2)*HS(I)
      A3=SIGR(4)*HR(I)+SIGR(6)*HS(I)
      A4=SIGR(6)*HR(I)+SIGR(5)*HS(I)
      A5=SIGR(7)*HR(I)+SIGR(8)*HS(I)
      HH=.5*H(I)
      HHR=.5*HR(I)
      HHS=.5*HS(I)
      B1=SIGR(7)*HH+SIGR(4)*HHR+SIGR(6)*HHS
      B2=SIGR(8)*HH+SIGR(6)*HHR+SIGR(5)*HHS
      B3=SIGR(4)*HH
      B4=SIGR(5)*HH
      B5=SIGR(6)*HH
      DO 20 J=I,NNO
      N=NS
      A=HR(J)*A1+HS(J)*A2
      E=HR(J)*A3+HS(J)*A4
      B=H(I)*(SIGR(7)*HR(J)+SIGR(8)*HS(J))+E
      C=H(J)*A5+E
      D=A*FB
      S1=H(J)*B1+HR(J)*B3+HS(J)*B5
      S2=H(J)*B2+HR(J)*B5+HS(J)*B4
      M=1
      DO 10 II=1,3
      DO 10 JJ=II,3
      F(M)=S1*(VR(II)*VT(JJ)+VT(II)*VR(JJ))
     1    +S2*(VS(II)*VT(JJ)+VT(II)*VS(JJ))
   10 M=M+1
C     ----------------------------------------------
C     UPPER TRIANGULAR PORTION OF EACH 6*6 PARTITION
C     ----------------------------------------------
      S(N)=S(N)+A
      S(N+4)=S(N+4)+C*VT(3)
      S(N+5)=S(N+5)-C*VT(2)
      N=N+LROW
      S(N)=S(N)+A
      S(N+2)=S(N+2)-C*VT(3)
      S(N+4)=S(N+4)+C*VT(1)
      N=N+LROW-1
      S(N)=S(N)+A
      S(N+1)=S(N+1)+C*VT(2)
      S(N+2)=S(N+2)-C*VT(1)
      N=N+LROW-2
      S(N)=S(N)+D*F1+F(1)
      S(N+1)=S(N+1)+D*F2+F(2)
      S(N+2)=S(N+2)+D*F3+F(3)
      N=N+LROW-3
      S(N)=S(N)+D*F4+F(4)
      S(N+1)=S(N+1)+D*F5+F(5)
      N=N+LROW-4
      S(N)=S(N)+D*F6+F(6)
      IF (I.EQ.J) GOTO 20
C     --------------------------------------------------------------
C     FILL IN REMAINING COEFFICIENTS FOR OFF-DIAGONAL 6*6 PARTITIONS
C     --------------------------------------------------------------
      N=NS+3*(LROW-2)
      S(N+1)=S(N+1)-B*VT(3)
      S(N+2)=S(N+2)+B*VT(2)
      N=N+LROW-4
      S(N)=S(N)+B*VT(3)
      S(N+2)=S(N+2)-B*VT(1)
      S(N+3)=S(N+3)+D*F2+F(2)
      N=N+LROW-5
      S(N)=S(N)-B*VT(2)
      S(N+1)=S(N+1)+B*VT(1)
      S(N+3)=S(N+3)+D*F3+F(3)
      S(N+4)=S(N+4)+D*F5+F(5)
   20 NS=NS+6
      NS=N+6
   30 LROW=LROW-6
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHGEO2 (S,SIGR,HR,HS,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     EVALUATES GEOMETRIC CONTRIBUTION TO TANGENTIAL STIFFNESS MATRIX
C     THIS VERSION IS BASED ON VON-KARMAN'S STRAIN EXPANSION
C	------------------------------------------------------
C     S(1176)       = GEOMETRIC STIFFNESS STORED UPPER TRIANGULAR
C                     ROW-WISE
C     SIGR(8)       = WEIGHTED STRESS-RESULTANT VECTOR
C     HR(8),HS(8)   = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     NNO           = NUMBER OF NODES FOR ELEMENT
C     ---------------------------------------------------------------
      DIMENSION S(*),SIGR(8),HR(8),HS(8)
C
      NS=1
      LROW=6*NNO
      DO 20 I=1,NNO
      A1=SIGR(1)*HR(I)+SIGR(3)*HS(I)
      A2=SIGR(3)*HR(I)+SIGR(2)*HS(I)
      DO 10 J=I,NNO
      N=NS
      A=HR(J)*A1+HS(J)*A2
C     ----------------------------------------------
C     UPPER TRIANGULAR PORTION OF EACH 6*6 PARTITION
C     ----------------------------------------------
      S(N)=S(N)+A
      N=N+LROW
      S(N)=S(N)+A
      N=N+LROW-1
      S(N)=S(N)+A
      IF (I.EQ.J) GOTO 10
      N=NS+2*LROW-3
   10 NS=NS+6
      NS=N+3*LROW-6
   20 LROW=LROW-6
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SHMASS (S,H,VR,VS,DVOL,DEN,TH,IMASS,NNO,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------
C     LUMPED, CONSISTENT DIAGONAL AND CONSISTENT MASS MATRICES
C	--------------------------------------------------------
C     S(NWS)      = STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     H(8)        = SHAPE FUNCTIONS
C     VR(3),VS(3) = LOCAL DIRECTION COSINE VECTORS
C     DVOL        = INTEGRATION FACTOR
C     DEN         = MASS DENSITY
C     IMASS       = 1   LUMPED MASS
C     IMASS       = 2   CONSISTENT DIAGONAL MASS
C     IMASS       = 3   CONSISTENT MASS
C     NNO         = NUMBER OF NODES
C     NEF         = NUMBER OF DEGREES OF FREEDOM
C     --------------------------------------------------------
      DIMENSION S(1),H(8),VR(3),VS(3),F(6),SM(6)
C
      FA=DVOL*DEN*TH
      FB=FA*TH*TH/12.
      N=1
      DO 10 I=1,3
      DO 10 J=I,3
      F(N)=FB*(VR(I)*VR(J)+VS(I)*VS(J))
   10 N=N+1
      N=1
      LROW=NEF
      GOTO (220,120,20),IMASS
C     CALL GOTOER
C     ----------------------
C     CONSISTENT MASS MATRIX
C     ----------------------
   20 NS=1
      DO 100 I=1,NNO
      DO  90 J=I,NNO
      N=NS
      HH=H(I)*H(J)
      SA=HH*FA
      DO 30 K=1,6
   30 SM(K)=HH*F(K)
      DO 40 K=1,3
      S(N)=S(N)+SA
   40 N=N+LROW+1-K
      S(N)=S(N)+SM(1)
      S(N+1)=S(N+1)+SM(2)
      S(N+2)=S(N+2)+SM(3)
      N=N+LROW-3
      S(N)=S(N)+SM(4)
      S(N+1)=S(N+1)+SM(5)
      N=N+LROW-4
      S(N)=S(N)+SM(6)
      IF (I.EQ.J) GOTO 90
      N=NS+4*LROW-7
      S(N)=S(N)+SM(2)
      N=N+LROW-5
      S(N)=S(N)+SM(3)
      S(N+1)=S(N+1)+SM(5)
   90 NS=NS+6
      NS=N+3
  100 LROW=LROW-6
      RETURN
C     -------------------------------
C     CONSISTENT DIAGONAL MASS MATRIX
C     -------------------------------
  120 DO 200 I=1,NNO
      HH=H(I)*H(J)
      SA=FA*HH
      DO 140 K=1,3
      S(N)=S(N)+SA
  140 N=N+LROW+1-K
      S(N)=S(N)+HH*F(1)
      N=N+LROW-3
      S(N)=S(N)+HH*F(4)
      N=N+LROW-4
      S(N)=S(N)+HH*F(6)
      N=N+LROW-5
  200 LROW=LROW-6
      RETURN
C     ------------------
C     LUMPED MASS MATRIX
C     ------------------
  220 FA=FA/(2*NNO-4)
      SA=FA
      DO 300 I=1,NNO
      IF (I.EQ.5) SA=2.*SA
      DO 240 K=1,3
      S(N)=S(N)+SA
  240 N=N+LROW+1-K
      N=N+3*LROW-12
  300 LROW=LROW-6
      RETURN
      END
C
C=====================================================================
C=====================================================================
C	END 4 NODE SHELL ELEMENT ROUTINES
C=====================================================================


C=====================================================================
      SUBROUTINE SHBMAT (NNO,H,HD,VR,VS,VT,XJI,HR,HS,B)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

C     ----------------------------------------------------------------
C     EVALUATES LOCAL STRAIN-DISPLACEMENT MATRIX
C	------------------------------------------
c     NNO              = NUMBER OF NODES FOR ELEMENT
C     H(8)             = SHAPE FUNCTIONS
C     HD(2,8)          = SHAPE FUNCTION DERIVATIVES W.R.T RN,SN
C     HR(8),HS(8)      = SHAPE FUNCTION DERIVATIVES W.R.T R,S RESP.
C     VR(3),VS(3),VT(3)= LOCAL DIRECTION COSINE VECTORS
C     XJI(4)           = INVERSE JACOBIAN MATRIX STORED COLUMN-WISE
C     B(240)           = STRAIN-DISPLACEMENT MATRIX STORED COLUMN-WISE
C     ----------------------------------------------------------------
      DIMENSION H(8),HD(2,8),VR(3),VS(3),VT(3),XJI(4),HR(8),HS(8)
      DIMENSION B(240)

      M=1
      N=16
      DO 20 I=1,NNO
      HR(I)=XJI(1)*HD(1,I)+XJI(3)*HD(2,I)
      HS(I)=XJI(2)*HD(1,I)+XJI(4)*HD(2,I)
      A1=HR(I)
      A2=HS(I)
      A3=H(I)
      DO 10 J=1,3
      F1=A2*VR(J)
      F2=A1*VS(J)
      B(M)=A1*VR(J)
      B(M+1)=A2*VS(J)
      B(M+2)=F1+F2
      B(M+3)=A1*VT(J)
      B(M+4)=A2*VT(J)
      B(N)=F2
      B(N+1)=-F1
      B(N+2)=B(M+1)-B(M)
      B(N+3)=A3*VS(J)
      B(N+4)=-A3*VR(J)
      M=M+5
      N=N+5
   10 CONTINUE
      M=M+15
      N=N+15
   20 CONTINUE
C
      RETURN
      END
C
C=====================================================================
C=====================================================================
      SUBROUTINE SHDELA8 (DR,DVOL,SLA,SLB)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------
C     ELASTIC RIGIDITIES FOR MEMBRANE-PLATE OR SHELL ELEMENT
C	------------------------------------------------------
C     DR(8,8) = ELASTIC MEMBRANE,BENDING AND SHEAR RIGIDITIES
C     DVOL    = WT*DET
C     SLA,SLB = SHEAR LOCKING CONSTRAINT FACTORS
C     -------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION DR(64)
C
      CALL CLEARA (DR,64)
      TH3  = TH*TH*TH/12.
      FACM = DVOL*TH
      FACB = DVOL*TH3
C
      DR(1)  = FACM*A1
      DR(10) = FACM*A1
      DR(2)  = FACM*B1
      DR(9)  = FACM*B1
      DR(19) = FACM*C1
      DR(28) = FACB*A1
      DR(37) = FACB*A1
      DR(29) = FACB*B1
      DR(36) = FACB*B1
      DR(46) = FACB*C1
      DR(55) = FACM*SLA*D1
      DR(64) = FACM*SLB*D1
C
      RETURN
      END
C


C=====================================================================
C	END 8 NODE SHELL ELEMENT ROUTINES
C=====================================================================












