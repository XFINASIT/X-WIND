C	=====================================================================
      SUBROUTINE SHELL4J(PROPM,PROPG,NODEX,WA,AMV,S,COORD,
	1					EDIS,EDISI,RE,MWG,FIN,MSET)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
C
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /DYNA/  CDEN,IMASS

      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET

	COMMON /HARD/  HP,DEN
C
      DIMENSION PROPM(5),PROPG(*),NODEX(*),WA(MWG,*),S(1176),COORD(3,8)
      DIMENSION EDIS(48),EDISI(48),RE(48),REDIS(48),COORDI(3,8)
      DIMENSION DR(64),H(8),HD(2,8),XJI(4),HR(8),HS(8),B(240)
      DIMENSION BDRL(24),DISD(12),EPS(8),EPSQ(8),SIGR(8)
      DIMENSION VR(3),VS(3),VT(3),SL(4,2),FF(6),NOD(4)
      DIMENSION COVR(3),COVS(3)
	DIMENSION BA(4,120),FJ(4),RRN(4),SSN(4)
C
	DIMENSION AMV(3)
C
	DIMENSION FIN(NEF)

	DIMENSION BMX(12,24),BMB(8,24)
C
      EQUIVALENCE (APEL,IPEL)
C
      DATA SL /2.5E5,.1334,.1334,.1334,1.E20,19.2,4.00,4.00/
C
      PI=3.1415926535898
      CALL CLEARA (BDRL,24)
C     ---------------------------------

      CALL HOKLAW (PROPM,PROPG,2)


C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C     ------------------------------------------------------
      CALL SHINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)


	RE = 0.0
	KKK = 0
	DO I = 1,24
	DO J = I,24
	KKK = KKK + 1
	S(KKK) = 0.0
	ENDDO
	ENDDO

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


	GOTO 120      
C	IF (NGP.LE.2) GOTO 120  !ALL NODE HAS AN EQUAL THICKNESS
C     ----------------------------------------
C     INTERPOLATE NODAL THICKNESSES (NGP.GT.2) EACH NODE HAVE DIFFERENT THICKNESS
C     ----------------------------------------
      TH=0.0
      DO 100 I=1,NNO
  100 TH=TH+H(I)*PROPG(I+1) 
	GO TO 130    

  120 TH=PROPG(2)
      
130	IF (ITASK.NE.5) GOTO 140

C     ---------------------
C     MASS MATRIX (ITASK=5)
C     ---------------------
C	CONMSS ADDED BY SONGSAK FOR RC SHELL
	IF(MTMOD.EQ.5.OR.MTMOD.EQ.6) THEN
	CALL CONMSS(TH,MSET,RHORC)
	CALL SHMASS (S,H,VR,VS,DVOL,RHORC,TH,IMASS,NNO,NEF)
	ELSE
      CALL SHMASS (S,H,VR,VS,DVOL,PROPM(5),TH,IMASS,NNO,NEF)
	ENDIF
      GOTO 900


  140 CALL SHBJ1(BMX,BMB,XJI,H,HD,VR,VS,VT)

	SLR = 1.0
	SLS = 1.0

C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES (IFREF.EQ.0)
C     -----------------------------------------
      IF (IFREF.EQ.0) CALL SHDELA (DR,DVOL,SLR,SLS)

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
      

	DO IMX = 1,12
	DO JMX = 1,24
	DISD(IMX) = DISD(IMX) + BMX(IMX,JMX)*REDIS(JMX)
	ENDDO
	ENDDO
	DISD(7) = DISD(7) + DISD(8)


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

C	LINEAR
 410  CALL MPSIGA (EPS,SIGR)

      DO I=1,8
	WA(I,IPT) = SIGR(I)
	ENDDO

      GO TO 500

C	COMPOSITE 
 465  CALL COMSIG (EPS,SIGR)

      DO I=1,8
	WA(I  ,IPT) =  EPS(I)
	WA(I+8,IPT) = SIGR(I)
	ENDDO

      GO TO 500

C	IVANOV ELASTO PLASTICITY
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

C	MULTILAYER VON-MISES ELASTO PLASTICITY
 460	IF (HP.EQ.0.) THEN
	CALL MLAYER (WA(1,IPT),EPS,SIGR,DR,DVOL)
      IPEL=2
      DR(55)=DR(55)
      DR(64)=DR(64)
	GOTO 500

	ELSE
	CALL MLAYERH (WA(1,IPT),EPS,SIGR,DR,DVOL)
	END IF
      IPEL=2
      DR(55)=DR(55)
      DR(64)=DR(64)

	GO TO 500

C	MULTILAYER REINFORCED CONCRETE (ELASTO PLASTICITY)
 467	CALL RCLAYR (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET)
	IPEL=2
	DR(55)=DR(55)
      DR(64)=DR(64)

	GO TO 500

C	MULTILAYER REINFORCED CONCRETE (EPF MODEL)
 468	CALL RCLAYRE (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET)
	IPEL=2
	DR(55)=DR(55)
      DR(64)=DR(64)

	GOTO 500

C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510 I=1,8
 510  SIGR(I)=SIGR(I)*DVOL


      IF (ITASK.LE.3.AND.ISOLOP.NE.4) GOTO 520 


      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900

520	CONTINUE

	DO IMB = 1,24
	DO JMB = 1,8
	RE(IMB) = RE(IMB) + BMB(JMB,IMB)*SIGR(JMB)
	ENDDO
	ENDDO

C     ------------------------------------------------------
C     FIND LINEAR CONTRIBUTION TO STIFFNESS MATRIX (IFREF=0)
C     ------------------------------------------------------
      IF (IFREF ) 900,700,900

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
      CALL SHKJ1(S,DR,BMB,BDRL)
C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 810

 800  CALL SHGEO1 (S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)

 810  TIM(12)=TIM(12)+TIM2


 900  CONTINUE

c	KKK = 0
c	DO I = 1,24
c	DO J = I,24
c	KKK = KKK + 1
c	if(I.EQ.J) write(*,*) MEL,DR(55),DR(64) !S(KKK),I
c	ENDDO
c	ENDDO
c	stop

	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF

      RETURN

      END
C
C	=====================================================================

C	=====================================================================
      SUBROUTINE SHBJ1(BM,BMB,XJI,H,HD,VR,VS,VT)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION XJI(2,2),HD(2,4),HR(4),HS(4),BM(12,24),H(4),TM(24,24)
	DIMENSION VR(3),VS(3),VT(3),BMB(8,24)

	
      DO I = 1,4
      HR(I)=XJI(1,1)*HD(1,I)+XJI(1,2)*HD(2,I)
      HS(I)=XJI(2,1)*HD(1,I)+XJI(2,2)*HD(2,I)
	ENDDO


	BM = 0.0

	DO I = 1,4

	BM(1,6*I-5)  =  HR(I)
	BM(2,6*I-4)  =  HS(I)
	BM(3,6*I-5)  =  HS(I)
	BM(4,6*I-4)  =  HR(I)
	BM(5,6*I-1)  =  HR(I)
	BM(6,6*I-2)  = -HS(I)
	BM(7,6*I-2)  = -HR(I)
	BM(8,6*I-1)  =  HS(I)
	BM(9,6*I-3)  =  HR(I)
	BM(10,6*I-1) =  H(I)
	BM(11,6*I-3) =  HS(I)
	BM(12,6*I-2) = -H(I)

	ENDDO

	BMB = 0.0 

	DO I = 1,4

	BMB(1,6*I-5)  =  HR(I)
	BMB(2,6*I-4)  =  HS(I)
	BMB(3,6*I-5)  =  HS(I)
	BMB(3,6*I-4)  =  HR(I)
	BMB(4,6*I-1)  =  HR(I)
	BMB(5,6*I-2)  = -HS(I)
	BMB(6,6*I-2)  = -HR(I)
	BMB(6,6*I-1)  =  HS(I)
	BMB(7,6*I-3)  =  HR(I)
	BMB(7,6*I-1)  =  H(I)
	BMB(8,6*I-3)  =  HS(I)
	BMB(8,6*I-2)  = -H(I)

	ENDDO

	DO I = 1,4

	TM(6*I-5,6*I-5) = VR(1)
	TM(6*I-4,6*I-5) = VR(2)
	TM(6*I-3,6*I-5) = VR(3)

	TM(6*I-5,6*I-4) = VS(1)
	TM(6*I-4,6*I-4) = VS(2)
	TM(6*I-3,6*I-4) = VS(3)

	TM(6*I-5,6*I-3) = VT(1)
	TM(6*I-4,6*I-3) = VT(2)
	TM(6*I-3,6*I-3) = VT(3)

	TM(6*I-2,6*I-2) = VR(1)
	TM(6*I-1,6*I-2) = VR(2)
	TM(6*I-0,6*I-2) = VR(3)

	TM(6*I-2,6*I-1) = VS(1)
	TM(6*I-1,6*I-1) = VS(2)
	TM(6*I-0,6*I-1) = VS(3)

	TM(6*I-2,6*I-0) = VT(1)
	TM(6*I-1,6*I-0) = VT(2)
	TM(6*I-0,6*I-0) = VT(3)

	ENDDO

	BM = MATMUL(BM,TRANSPOSE(TM))

	BMB = MATMUL(BMB,TRANSPOSE(TM))


	RETURN 

	END

C	=====================================================================
C	=====================================================================


C	=====================================================================
      SUBROUTINE SHKJ1(S,DR,BMB,BDRL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=====================================================================

      DIMENSION SK(24,24),DR(8,8),BMB(8,24),BDRL(24,1),S(*)
C
      SK = 0.0

      FDRL=10.*DR(3,3)


	SK = SK + MATMUL(BDRL,TRANSPOSE(BDRL))*FDRL


	SK = SK + MATMUL(TRANSPOSE(BMB),MATMUL(DR,BMB))

	K = 0
	DO I = 1,24
	DO J = I,24
	K = K+1
	S(K) = S(K) + SK(I,J)
	ENDDO
	ENDDO


      RETURN

      END
C	=====================================================================
C	=====================================================================






	
	
