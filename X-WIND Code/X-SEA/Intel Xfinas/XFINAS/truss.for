C=====================================================================
C	START 3D TRUSS ELEMENT ROUTINES
C=====================================================================

	SUBROUTINE TRUS3DM (PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,
	1					FIN,MSET,CABFF,DISLI)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      DIMENSION PROPM(1),PROPG(1),EDIS(1),EDISI(1)
      DIMENSION COORD(1),ELOD(1),S(1)
	DIMENSION WA(1) !ADDED JAN2004 - CABLE

	DIMENSION CABFF(1)

      DIMENSION DISLI(1)
      
	GOTO (110,120,130,140,150,160,170,180) ISTYP
110	CONTINUE
	RETURN

120	CONTINUE
	RETURN

C	TRUSS WITH GEOMETRIC STIFFNESS ELEMENT
130	CALL STRUSS2(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG)
	RETURN

C	INEXTENSIBLE CABLE - ELAXTIC TRUSS ELEMENT
140	CALL SPCABLE(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG,FIN,
	1			 CABFF)
	RETURN

150	CALL NONGAP(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG,FIN)
	RETURN


C	BONDLINK ELEMENT SONGSAK FEB2006
160	CALL BONDLNK(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG,FIN
	1		    ,MSET)
	RETURN

C	CATENARY - ELASTIC CABLE
170   CALL CMCABLE(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG,FIN)
      !CALL CNCABLE(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG,FIN)
      RETURN
    
C     MOORING      
180   CALL CNCABLE_NEW_MOORING(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG,FIN)
      !IF(NMOOR.EQ.1) THEN !TESTING MOORING LINE BY BJ 
      !!ORIGINAL => CALL CNCABLE_MOORING (COORD,EDIS,PROPM,PROPG,RHO,GAMMA,WX,WY,WZ,FHOR,NSEG,NPTS,S,RE,CCOOR)
      !CALL CNCABLE_MOORING (COORD,EDIS,PROPM,PROPG,RHO,GAMMA,WX,WY,WZ,FHOR,NSEG,NPTS,S,ELOD,FIN,DISLI,WA)
      !CALL CNCABLE_NEW_MOORING(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG,FIN)
      !ENDIF
      RETURN

	END
C=====================================================================
C==================================================================
      SUBROUTINE MATRIX3MUL0(A,B,C,D,N1,N2,N3)
C     ---------  
C     COMPUTE THE PRODUCT OF THREE MATRICES
C     TO COMPUTE  TRANSPOSE(A)*B*C=D
C     AT(N1,N2),B(N2,N2),C(N2,N3)=D(N1,N3)
C     ---------  
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N2,N1),B(N2,N2),C(N2,N3),D(N1,N3)
	DIMENSION AT(N1,N2),AB(N1,N2)
	DO 10 I=1,N1
	DO 10 J=1,N3
	D(I,J)=0.0
  10  CONTINUE
      DO 15 I=1,N1
	DO 15 J=1,N2
	AB(I,J)=0.0
  15  CONTINUE

C     TRANSPOSE A(N2,N1) TO AT(N1,N2)
      DO 20 I=1,N1
	DO 20 J=1,N2
	AT(I,J)=A(J,I)
  20  CONTINUE
      
C     COMPUTE AT*B=AB
      DO 30 I=1,N1
	DO 30 J=1,N2
	DO 30 K=1,N2
	AB(I,J)=AB(I,J)+AT(I,K)*B(K,J)
  30  CONTINUE

C     COMPUTE AB*C=D
      DO 40 I=1,N1
	DO 40 J=1,N3
	DO 40 K=1,N2
	D(I,J)=D(I,J)+AB(I,K)*C(K,J)
  40  CONTINUE
      RETURN
      END
C=====================================================================
C	START 3D TRUSS ELEMENT ROUTINES
C=====================================================================
C	CREATED BY WANG, SEPT 10,2002
C	RE-ARRANGED BY GILSON - APRIL2003
C=====================================================================
      SUBROUTINE STRUSS2(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,
	1				   MWG)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     ELEMENT PROPERTIES
C     ------------------
C	HN(2,6)		TRANSFORMATION MATRIX
C	DISP(2,1)	LOCAL DISPLACEMENT
C	QSTR		NONLINEAR STRAIN
C	RL(2,1)		LOCAL INTERNAL RESISTING FORCE 
C	RE(6)		GLOBAL INTERNAL RESISTING FORCE
C	SL(2,2)		LOCAL LINEAR STIFFNESS MATRIX
C	S1(6,6)		GLOBAL LINEAR STIFFNESS MATRIX
C	S2(6,6)		GLOBAL NONLINEAR STIFFNESS MATRIX
C	S(21)		TOTAL STIFFNESS (UPPER-TRIANGLE)
C     ------------------------------------------------------------- 
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

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
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
	COMMON /TRUSS/ AREA

C 
      DIMENSION PROPM(5),PROPG(5),EDIS(6),EDISI(6)
      DIMENSION COORD(6)
C
	DIMENSION HN(6,6),DISP(2,1),RE(6)
	DIMENSION SL(6,6),S1(6,6),S2(6,6),S(21)

	DIMENSION WA(1)

C     FOR DGEMM FUNCTION	
	DIMENSION COUT1(6,6)

C	--------------------------------
C	INITIALIZATION OF SOME VARIABLES
C	--------------------------------
	DO 10 I=1,6
	 DO 10 J=1,6
	  S1(I,J)=0.0
	  S2(I,J)=0.0
10	CONTINUE
	AREA = PROPG(2)

	IFOPT = INT(PROPG(4))  !IFOPT    1=INITIAL FORCE WILL ADD TO LOAD VECTOR    0=INITIAL FORCE WILL NOT ADD TO LOAD VECTOR FEB09
	PRET = PROPG(3) !PRETENSION FORCE
	STSO = PRET/AREA !PRETENSION STRESS
C	---------------------------------------
C	SET VALUES FOR LINEAR STRESS-STRAIN LAW 
C	INITIALISATION OF INTEGRATION RULE
C	---------------------------------------
	CALL HOKLAW (PROPM,PROPG,1)

C	---------------------------------------
C	TO COMPUTE THE LENGTH OF EACH TRUSS BAR
C	---------------------------------------
	A1=COORD(4)-COORD(1)
	A2=COORD(5)-COORD(2)
	A3=COORD(6)-COORD(3)
	WL=SQRT(A1**2+A2**2+A3**2)

C	--------------------------------------------------
C	TRANSFORMATION MATRIX HN(2*6) FROM GLOBAL TO LOCAL
C	--------------------------------------------------
	CALL TRANMATR(COORD,HN)

C	-----------------------------------
C	MASS MATRIX FOR FREQUENCY ANALYSIS
C	-----------------------------------
	IF(ITASK.NE.5) GOTO 50
	CALL SPTRUMAS(S,COORD,PROPM,PROPG(2))
	GOTO 900
	
50	CONTINUE
C	---------------------------------------
C	GLOBAL NODAL DISPS TO LOCAL NODAL DISPS
C	---------------------------------------
	DISP(1,1)=HN(1,1)*EDIS(1)+HN(1,2)*EDIS(2)+HN(1,3)*EDIS(3)
	DISP(2,1)=HN(1,1)*EDIS(4)+HN(1,2)*EDIS(5)+HN(1,3)*EDIS(6)
C	---------------
C	ORIGINAL LENGTH
C	---------------
	A10 = (COORD(4)-EDIS(4)) - (COORD(1)-EDIS(1))
	A20 = (COORD(5)-EDIS(5)) - (COORD(2)-EDIS(2))
	A30 = (COORD(6)-EDIS(6)) - (COORD(3)-EDIS(3))
	WL0=SQRT(A10**2+A20**2+A30**2)

C	-------------------
C	LINEAR LOCAL STRAIN
C	-------------------
100	STRAIN=(DISP(2,1)-DISP(1,1))/WL

C	------------------------------------------------------------
C	FOR NLOPT>1 SUBTRACT QUADRATIC STRAIN TERMS (ALMANSI STRAIN)
C	------------------------------------------------------------
	IF (NLOPT.LE.1) GOTO 400
	QSTR=((EDIS(4)-EDIS(1))**2+(EDIS(5)-EDIS(2))**2
     +	   +(EDIS(6)-EDIS(3))**2)*0.5/WL0**2

	STRAIN = STRAIN-QSTR

	STRAIN = (WL-WL0)/WL0
C	--------------------------
C	COMPUTE LOCAL AXIAL STRESS 
C	--------------------------
400	STRESS=PROPM(1)*STRAIN + STSO !ADD ALSO PRETENSION STRESS

C	ELASTO PLASTIC SONGSAK FEB2006
	IF(MTMOD.EQ.3) THEN

	STRN  = WA(2)
	STRS  = WA(3) + STSO !ADD ALSO PRETENSION STRESS
	EPSTN = WA(4)

	DSTN = STRAIN - STRN
	DSTS = PROPM(1)*DSTN

	CALL	TRUSPLS(STRS,DSTS,PROPM(1),PROPM(3),PROPM(4),
	1				EPSTN)

	STRESS = STRS

	WA(2) = STRAIN
	WA(3) = STRESS - STSO
	WA(4) = EPSTN

	ENDIF
C	==============================

C	WA=STRESS
C	--------------------------------
C	COMPUTE INTERNAL RESISTING FORCE
C	--------------------------------
	P=STRESS*PROPG(2)

      PP = P
      IF(IFOPT.EQ.0) PP = PP - PRET !NOT INCLUDE PRETENSION IF IFOPT = 0
      
	WA(1) = P  !CHANGE HERE SEE ALSO ELPRIN

      
	IF (IFEIG.EQ.0.AND.ISOLOP.EQ.4) GOTO 800

C	WRITE(*,*) MEL,STRAIN,WL
C	P=STRESS*AREA

500	RE(1) = -PP*HN(1,1)
	RE(2) = -PP*HN(1,2)
	RE(3) = -PP*HN(1,3)
	RE(4) =  PP*HN(1,1)
	RE(5) =  PP*HN(1,2)
	RE(6) =  PP*HN(1,3)
C	-----------------------------------------------------
C	COMPUTE LINEAR STIFFNESS SL(2*2) IN LOCAL COORDIANTES 
C	-----------------------------------------------------
1000	IF (IFREF) 900,750,900

750	SL(1,1) =  1.0*PROPG(2)*PROPM(1)/WL
	SL(1,4) = -1.0*PROPG(2)*PROPM(1)/WL
	SL(4,1) = -1.0*PROPG(2)*PROPM(1)/WL
	SL(4,4) =  1.0*PROPG(2)*PROPM(1)/WL

C	----------------------------------
C	TRANSFER LOCAL STIFFNESS TO GLOBAL
C	----------------------------------

C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M) 
      ALPHA = 1.0D0
      BETA  = 1.0D0
      COUT1 = 0.0D0 
      CALL DGEMM('N','N',6,6,6,ALPHA,TRANSPOSE(HN),6,SL,6,BETA,COUT1,6) 
      S1    = 0.0D0
      CALL DGEMM('N','N',6,6,6,ALPHA,COUT1,6,HN,6,BETA,S1,6)  
C	S1 = MATMUL(MATMUL(TRANSPOSE(HN),SL),HN)

	IF (ISOLOP.EQ.4) GOTO 1100
C	----------------------------------
C	COMPUTE NONLINEAR GLOBAL STIFFNESS
C	----------------------------------
800	W=P/WL

	S2(1,1) = W
	S2(2,2) = W
	S2(3,3) = W
	S2(4,4) = W
	S2(5,5) = W
	S2(6,6) = W
	S2(1,4) =-W
	S2(2,5) =-W
	S2(3,6) =-W
	S2(4,1) =-W
	S2(5,2) =-W
	S2(6,3) =-W

C	--------------------------------------------------
C	COMPUTE FULL STIFFNESS MATRIX S(21) UPPER-TRIANGLE
C	--------------------------------------------------
1100	KK=1
	DO 700 I=1,6
	 DO 700 J=I,6
	  S(KK)=S1(I,J)+S2(I,J)
	  KK=KK+1
700	CONTINUE


900	RETURN
	END
C
C=====================================================================
C=====================================================================
      SUBROUTINE INISIGCABL(IEG,IGIDM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------------------
C     READS, GENERATES AND PRESETS INITIAL STRESSES AT GAUSS POINTS
C	-------------------------------------------------------------
C     WA(NWA,NELE)  = WORKING ARRAY STORING STRESSES AND STRAINS
C     -------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
C
      DIMENSION IGIDM(1)
	ALLOCATABLE WA(:)

	CALL MINTFIL('GWOK',MWA,IEG,3,0)
      ALLOCATE(WA(MWA))

C
	KEL = 0
	READ (ITI,*)
 100  READ (ITI,*)  MEMBA,CFORE

	CALL ELEREODER(IGIDM,NELE,MEMBA)
	CALL ADREWT(IEG,MEMBA,WA,'RED')
      WA(3) = CFORE
	CALL ADREWT(IEG,MEMBA,WA,'WRT')

	KEL = KEL + 1
      IF (KEL.LT.NELE)  GOTO 100
C
      DEALLOCATE(WA)


      RETURN
      END
C
C	==========================================================
	SUBROUTINE SPTRUMAS(S,COORD,PROPM,AREA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION S(21),S1(6,6),COORD(6),PROPM(5),SL(6,6),HN(6,6)
C	--------------------
C	TRANSFOMATION MATRIX
C	--------------------
	A1 = COORD(4)-COORD(1)
	A2 = COORD(5)-COORD(2)
	A3 = COORD(6)-COORD(3)
	WL = SQRT(A1**2+A2**2+A3**2)
	Ax = A1/WL
	Ay = A2/WL
	Az = A3/WL

C	FORMULATE TRANSFORMATION MATRIX
C	-------------------------------
	CALL TRANMATR(COORD,HN)

C	FORMULATE LOCAL MASS MATRIX
C	---------------------------
	SMAS = PROPM(5)*AREA*WL/6.0
	SL = 0.0
	SL(1,1) = 2.0
	SL(2,2) = 2.0
	SL(3,3) = 2.0
	SL(4,4) = 2.0
	SL(5,5) = 2.0
	SL(6,6) = 2.0

	SL(1,4) = 1.0
	SL(2,5) = 1.0
	SL(3,6) = 1.0
	SL(4,1) = 1.0
	SL(5,2) = 1.0
	SL(6,3) = 1.0
C	-------------------------------------
C	TRANSFER LOCAL MASS MATRIX TO GLOBAL
C	-------------------------------------
	S1 = MATMUL(MATMUL(TRANSPOSE(HN),SL),HN)

	KK=1
	DO 700 I=1,6
	 DO 700 J=I,6
	  S(KK)=SMAS*S1(I,J)
	  KK=KK+1
700	CONTINUE

      RETURN
      END
C	==========================================================
	SUBROUTINE TRANMATR(COORD,HN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION COORD(6),HN(6,6)
C	--------------------
C	TRANSFOMATION MATRIX
C	--------------------
	A1 = COORD(4)-COORD(1)
	A2 = COORD(5)-COORD(2)
	A3 = COORD(6)-COORD(3)
	WL = SQRT(A1**2+A2**2+A3**2)
	Ax = A1/WL
	Ay = A2/WL
	Az = A3/WL

C	LET DEFINE LOCAL AXIS SYSTEM
C	DEFINE A POINT ON THE LOCAL XY PLANE
C	P1 = COORD(4)+10.0
C	P2 = COORD(2)+5.0
C	P3 = COORD(3)+2.0
	IF(Az.EQ.0.0) THEN
	P1 = 0.0
	P2 = 0.0
	P3 = 1.0
	ELSEIF(Ay.EQ.0.0) THEN
	P1 = 0.0
	P2 = 1.0
	P3 = 0.0
	ELSEIF(Ax.EQ.0.0) THEN
	P1 = 1.0
	P2 = 0.0
	P3 = 0.0
	ELSE
	P1 = 0.0
	P2 = 0.0
	P3 = 1.0
	END IF

	Q1 = P1-COORD(1)
	Q2 = P2-COORD(2)
	Q3 = P3-COORD(3)	
	
	C1 = A2*Q3-A3*Q2
	C2 = A3*Q1-A1*Q3
	C3 = A1*Q2-A2*Q1
	WLZ = SQRT(C1**2+C2**2+C3**2)
	Cx = C1/WLZ
	Cy = C2/WLZ
	Cz = C3/WLZ

	B1 = Ay*Cz - Az*Cy
	B2 = Az*Cx - Ax*Cz
	B3 = Ax*Cy - Ay*Cx
	WLY = SQRT(B1**2+B2**2+B3**2)
	Bx = B1/WLY
	By = B2/WLY
	Bz = B3/WLY
	
C	TRANSFOMATION MATRIX
C	----------------------
	HN = 0.0
	HN(1,1)=Ax
	HN(1,2)=Ay
	HN(1,3)=Az
	HN(2,1)=Bx
	HN(2,2)=By
	HN(2,3)=Bz
	HN(3,1)=Cx
	HN(3,2)=Cy
	HN(3,3)=Cz	

	HN(4,4)=Ax
	HN(4,5)=Ay
	HN(4,6)=Az
	HN(5,4)=Bx
	HN(5,5)=By
	HN(5,6)=Bz
	HN(6,4)=Cx
	HN(6,5)=Cy
	HN(6,6)=Cz	

      RETURN
      END
C	===================================================
	SUBROUTINE TRVRST(COORD,VR,VS,VT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION COORD(6),VR(3),VS(3),VT(3)
C	--------------------
C	TRANSFOMATION MATRIX
C	--------------------
	A1 = COORD(4)-COORD(1)
	A2 = COORD(5)-COORD(2)
	A3 = COORD(6)-COORD(3)
	WL = SQRT(A1**2+A2**2+A3**2)
	Ax = A1/WL
	Ay = A2/WL
	Az = A3/WL

C	LET DEFINE LOCAL AXIS SYSTEM
C	DEFINE A POINT ON THE LOCAL XY PLANE
C	P1 = COORD(4)+10.0
C	P2 = COORD(2)+5.0
C	P3 = COORD(3)+2.0
	IF(Az.EQ.0.0) THEN
	P1 = 0.0
	P2 = 0.0
	P3 = 1.0
	ELSEIF(Ay.EQ.0.0) THEN
	P1 = 0.0
	P2 = 1.0
	P3 = 0.0
	ELSEIF(Ax.EQ.0.0) THEN
	P1 = 1.0
	P2 = 0.0
	P3 = 0.0
	ELSE
	P1 = 0.0
	P2 = 0.0
	P3 = 1.0
	END IF

	Q1 = P1-COORD(1)
	Q2 = P2-COORD(2)
	Q3 = P3-COORD(3)	
	
	C1 = A2*Q3-A3*Q2
	C2 = A3*Q1-A1*Q3
	C3 = A1*Q2-A2*Q1
	WLZ = SQRT(C1**2+C2**2+C3**2)
	Cx = C1/WLZ
	Cy = C2/WLZ
	Cz = C3/WLZ

	B1 = Ay*Cz - Az*Cy
	B2 = Az*Cx - Ax*Cz
	B3 = Ax*Cy - Ay*Cx
	WLY = SQRT(B1**2+B2**2+B3**2)
	Bx = B1/WLY
	By = B2/WLY
	Bz = B3/WLY
	
C	BASE VECTORS
C	------------
	VR(1) = Ax
	VR(2) = Ay
	VR(3) = Az
	VS(1) = Bx
	VS(2) = By
	VS(3) = Bz
	VT(1) = Cx
	VT(2) = Cy
	VT(3) = Cz	

      RETURN
      END
C	===================================================





