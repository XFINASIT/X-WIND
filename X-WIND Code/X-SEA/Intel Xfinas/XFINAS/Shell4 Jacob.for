C	=====================================================================
      SUBROUTINE SHELL4J(PROPM,PROPG,S,COORD,WA,
	1				   EDIS,EDISI,RE,FIN,MSET,AMV)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2               NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /DYNA/  CDEN,IMASS
	COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST

C
      DIMENSION PROPM(NMP),PROPG(NGP),S(1176),COORD(3,8)
      DIMENSION EDIS(48),EDISI(48),RE(48),REDIS(48),COORDI(3,8)
      DIMENSION DR(64),H(8),HD(2,8),XJI(4),HR(8),HS(8)
      DIMENSION DISD(12),EPS(8),EPSQ(8),SIGR(8)
      DIMENSION VR(3),VS(3),VT(3),FJ(4),FF(6)
	DIMENSION FIN(NEF),WA(NWG,NPT)

	DIMENSION BMX(12,24),BMB(8,24),BA(2,24),BANS(8,24)
	DIMENSION BDID1(8,24),BDID2(8,24),DID1(2,24),DID2(2,24)
	DIMENSION TM(24,24),BDRL(48),BDL(8,24)
	DIMENSION REE(24)

	DIMENSION RRN(4),SSN(4),AMV(3)


C     ---------------------------------
      CALL CLEARA (BDRL,48)
C     ---------------------------------

      CALL HOKLAW (PROPM,PROPG,2)

	
C     ------------------------------------------------------
C     INITIAL COORDS (COORDI) AND COROTATIONAL DISPS (REDIS)
C     ------------------------------------------------------
      CALL SHINIT (COORD,COORDI,EDIS,REDIS,NNO,NEF)

C     ------------------------------------------------------
	FJ     = 0.0
	XJI    = 0.0
	BANS   = 0.0
	RRN(1) = 1.0
	SSN(1) = 0.0
	RRN(2) =-1.0
	SSN(2) = 0.0
	RRN(3) = 0.0
	SSN(3) = 1.0
	RRN(4) = 0.0
	SSN(4) =-1.0
	DO ISP = 1,4
	RN = RRN(ISP)
	SN = SSN(ISP)
	CALL SHAP4D(RN,SN,H,HD)
	CALL SHJACB(NNO,COORD,HD,VR,VS,VT,DET,FJ,XJI)
	CALL TRANM1(VR,VS,VT,TM)
	CALL SHSJ1 (BA,FJ,XJI,H,HD,TM)
	CALL SDID1 (DID1,DID2,FJ,XJI,H,HD,TM)
	DO I = 1,24
	BANS(2*ISP-1,I) = BA(1,I)
	BANS(2*ISP-0,I) = BA(2,I)
	BDID1(2*ISP-1,I) = DID1(1,I)
	BDID1(2*ISP-0,I) = DID1(2,I)
	BDID2(2*ISP-1,I) = DID2(1,I)
	BDID2(2*ISP-0,I) = DID2(2,I)
	ENDDO
	ENDDO

C     ------------------------------------------------------
	DO I = 1,NEF
	RE(I) = 0.0
	ENDDO
	DO I = 1,(NEF*NEF-NEF)/2+NEF
	S(I) = 0.0
	ENDDO
	DO I = 1,8
	EPS(I)  = 0.0
	EPSQ(I) = 0.0
	ENDDO
C     ------------------------------------------------------

	IPT = 0
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
      DO 900  IGR=1,NGR
      RN = GLOC(IGR,NGR)
      DO 900  IGS=1,NGS
      SN = GLOC(IGS,NGS)
      WT = GWT(IGR,NGR)*GWT(IGS,NGS)
      IPT= IPT+1

C     -----------------------------------------------------
C     SHAPE FUNCTIONS (H) , SHAPE FUNCTION DERIVATIVES (HD)
C     -----------------------------------------------------
      CALL SHAP4D(RN,SN,H,HD)
C     ----------------------------------------------------
C     INVERSE JACOBIAN (XJI) , JACOBIAN coefficient F,
C	DETERMINANT (DET) AND STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------------------
      CALL SHJACB(NNO,COORD,HD,VR,VS,VT,DET,FJ,XJI)
	CALL TRANM1(VR,VS,VT,TM)

      DVOL=WT*DET

	
	CALL SHSTORL1 (COORD,PROPG,WA(1,IPT),NWG,VR,VS,VT)	!STORE SHELL LAX FOR STRESS CALCULATION JAN09

      SELECTCASE(MTMOD)
      CASE(2,5,6)
	CALL CANGLE(VR,VS,AMV,ANG)
	ENDSELECT
	      
C     ----------------------------------------------------
	TH=PROPG(2)
C     ----------------------------------------------------


	IF(ITASK.NE.5) GOTO 140
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
C     ----------------------------------------------------
C     ----------------------------------------------------
  140 CONTINUE 
C	CALL SHBJ1(BMX,XJI,H,HD,HR,HS,TM)
	CALL SDID2(BMX,BDID1,BDID2,XJI,H,HD,HR,HS,RN,SN,TM)
	CALL SHSJ2(BA,BANS,XJI,RN,SN)
	CALL SHBJ2(BMB,BA,XJI,H,HD,TM)

C     ----------------------------------------------------

C     -----------------------------------------
C     DETERMINE ELASTIC RIGIDITIES (IFREF.EQ.0)
C     -----------------------------------------
	SLR = 1.0
	SLS = 1.0
      CALL SHDELA(DR,DVOL,SLR,SLS)


C     -------------------------------------------------
C     REMOVE RIGID BODY TRANSLATIONS AND ROTATIONS FROM
C     TOTAL DISPLACEMENT VECTOR (NLOPT>1)
C     -------------------------------------------------
      IF (NLOPT.LE.1) GOTO 180

      CALL SHMDSP (COORD,COORDI,EDIS,REDIS,H,HD,VR,VS,VT,NNO)

C     ------------------------
C     DISPLACEMENT DERIVATIVES
C     ------------------------
	
  180 DISD = 0.0
	DO IMX = 1,12
	DO JMX = 1,24
	DISD(IMX) = DISD(IMX) + BMX(IMX,JMX)*REDIS(JMX)
	ENDDO
	ENDDO


C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
	DO IEP = 1,8
	EPS(IEP) = 0.0
	DO JEP = 1,24
	EPS(IEP) = EPS(IEP) + BMB(IEP,JEP)*REDIS(JEP)
	ENDDO
	ENDDO
C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400

      EPSQ(1) = 0.5*(DISD(1)*DISD(1)+DISD(3)*DISD(3)+DISD(5)*DISD(5))
      EPSQ(2) = 0.5*(DISD(2)*DISD(2)+DISD(4)*DISD(4)+DISD(6)*DISD(6))
      EPSQ(3) =      DISD(1)*DISD(2)+DISD(3)*DISD(4)+DISD(5)*DISD(6)

	EPSQ(4) =      DISD(1)*DISD(11) - DISD(3)*DISD(8)
	EPSQ(5) =      DISD(2)*DISD(12) - DISD(4)*DISD(9)
	EPSQ(6) =      DISD(1)*DISD(12) + DISD(2)*DISD(11)
	1			 - DISD(4)*DISD(8)  - DISD(3)*DISD(9)

      EPSQ(7) =      DISD(1)*DISD(10) - DISD(3)*DISD(7) 
      EPSQ(8) =      DISD(2)*DISD(10) - DISD(4)*DISD(7)


      DO 300  I=1,8
300   EPS(I)= EPS(I)-EPSQ(I)

400	CONTINUE

C     -----------------------------------------------
C	SELECT MATERIAL MODULE
C     -----------------------------------------------
	GO TO (410,420,430,440,450,460), MTMOD

C	LINEAR
410   CALL MPSIGA(EPS,SIGR)
	DO IWG = 1,NWG
	WA(IWG,IPT) = SIGR(IWG)
	ENDDO
	GO TO 500

420	CONTINUE
	GO TO 500

430	CONTINUE
	GO TO 500

440	CONTINUE
	GO TO 500

450	CALL RCLAYR (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GO TO 500


460	CALL RCLAYRE (WA(1,IPT),EPS,SIGR,MEL,TH,DR,DVOL,MSET,ANG)
	IPEL=2
	DR(55)=DR(55)*SLR
      DR(64)=DR(64)*SLS
	GO TO 500

C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
500   DO I=1,8
	SIGR(I) = SIGR(I)*DVOL
	ENDDO

	REE = MATMUL(TRANSPOSE(BMB),SIGR)

	DO IMB = 1,24
	RE(IMB) = RE(IMB) + REE(IMB)
	ENDDO


	IF (IFREF ) 900,700,900

C700   CALL DRILL(RN,SN,BDL,XJI,TM,HD)   
C      CALL SHKJ2(S,DR,BMB,BDL)  


700	CALL SHBDRL(NNO,H,HR,HS,VR,VS,VT,BDRL)
	CALL SHKJ1(S,DR,BMB,BDRL)



	     
C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 810

800   CALL SHGEO1(S,SIGR,H,HR,HS,VR,VS,VT,TH,NNO)

810   CONTINUE !TIM(12)=TIM(12)+TIM2


900   CONTINUE

	IF (ITASK.EQ.3) THEN
	DO  950 I = 1,NEF
	FIN(I) = RE(I)
950   CONTINUE
	ENDIF

      RETURN

      END
C
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHBJ1(BR,XJI,H,HD,HR,HS,TM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION XJI(2,2),HD(2,8),HR(8),HS(8),BM(12,24),H(8),TM(24,24)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION BR(12,24)
	
	HR = 0.0
	HS = 0.0
      DO I = 1,4
      HR(I)=XJI(1,1)*HD(1,I)+XJI(1,2)*HD(2,I)
      HS(I)=XJI(2,1)*HD(1,I)+XJI(2,2)*HD(2,I)
	ENDDO

	BM = 0.0

	DO I = 1,4

	BM(1,6*I-5)  =  HR(I)
	BM(2,6*I-5)  =  HS(I)
	BM(3,6*I-4)  =  HR(I)
	BM(4,6*I-4)  =  HS(I)
	BM(5,6*I-3)  =  HR(I)
	BM(6,6*I-3)  =  HS(I)
	BM(7,6*I-2)  =  H(I)
	BM(8,6*I-2)  =  HR(I)
	BM(9,6*I-2)  =  HS(I)
	BM(10,6*I-1) =  H(I)
	BM(11,6*I-1) =  HR(I)
	BM(12,6*I-1) =  HS(I)

	ENDDO



	BR = MATMUL(BM,TRANSPOSE(TM))


	RETURN

	END

C	=====================================================================
C	=====================================================================
      SUBROUTINE SHBJ2(BR,BA,XJI,H,HD,TM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION XJI(2,2),HD(2,8),HR(8),HS(8),H(8),BANS(8,24)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION BM(12,24),BR(8,24),TM(24,24),TEPS(8,12)
	DIMENSION BA(2,24)

	HR = 0.0
	HS = 0.0
      DO I = 1,4
      HR(I) = XJI(1,1)*HD(1,I) + XJI(1,2)*HD(2,I)
      HS(I) = XJI(2,1)*HD(1,I) + XJI(2,2)*HD(2,I)
	ENDDO
	

	TEPS = 0.0

	TEPS(1,1)  = 1.0
	TEPS(2,4)  = 1.0
	TEPS(3,2)  = 1.0
	TEPS(3,3)  = 1.0
	TEPS(4,11) = 1.0
	TEPS(5,9)  =-1.0
	TEPS(6,12) = 1.0
	TEPS(6,8)  =-1.0


	BM = 0.0

	DO I = 1,4

	BM(1 ,6*I-5)  =  HR(I)
	BM(2 ,6*I-5)  =  HS(I)
	BM(3 ,6*I-4)  =  HR(I)
	BM(4 ,6*I-4)  =  HS(I)
	BM(5 ,6*I-3)  =  HR(I)
	BM(6 ,6*I-3)  =  HS(I)
	BM(7 ,6*I-2)  =  H(I)
	BM(8 ,6*I-2)  =  HR(I)
	BM(9 ,6*I-2)  =  HS(I)
	BM(10,6*I-1)  =  H(I)
	BM(11,6*I-1)  =  HR(I)
	BM(12,6*I-1)  =  HS(I)

	ENDDO


	BR = MATMUL(TEPS,MATMUL(BM,TRANSPOSE(TM)))

	DO J = 1,24
	BR(7,J) = BA(1,J)
	BR(8,J) = BA(2,J)
	ENDDO


	RETURN

	END

C	=====================================================================
C	=====================================================================
      SUBROUTINE SHKJ1(S,DR,BMB,BDRL)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=====================================================================

      DIMENSION SK(24,24),DR(8,8),BMB(8,24),BDRL(48),S(1176),
	1		  SK1(24,24),BBRL(24,1),SKK(24,24)
C
      SK  = 0.0
	SK1 = 0.0

      FDRL = 0.01*DR(3,3) !10.0*DR(3,3)

	DO I = 1,24
	BBRL(I,1) = BDRL(I)
	ENDDO


	SK1 =  MATMUL(BBRL,TRANSPOSE(BBRL))*FDRL

C	WRITE(*,*) SK1
C	pause


	SK  =  MATMUL(TRANSPOSE(BMB),MATMUL(DR,BMB))

	DO I = 1,24
	DO J = 1,24
	SKK(I,J) = SK(I,J)  + SK1(I,J)
	ENDDO
	ENDDO


	K = 0
	DO I = 1,24
	DO J = I,24
	K = K+1
	S(K) = S(K) + SKK(I,J)
	ENDDO
	ENDDO


      RETURN

      END
C	================================================================
C	================================================================
C	================================================================
	SUBROUTINE SHSJ1(BA,FJ,XJI,H,HD,TM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION FJ(2,2),XJI(2,2),HD(2,8),HR(8),HS(8),H(8)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION BM(12,24),TM(24,24),TEPS(2,12)
	DIMENSION BA(2,24),BR(2,24)

	HR = 0.0
	HS = 0.0
      DO I = 1,4
      HR(I) = XJI(1,1)*HD(1,I) + XJI(1,2)*HD(2,I)
      HS(I) = XJI(2,1)*HD(1,I) + XJI(2,2)*HD(2,I)
	ENDDO
	

	TEPS = 0.0
	TEPS(1,5)  = 1.0
	TEPS(1,10) = 1.0
	TEPS(2,6)  = 1.0
	TEPS(2,7)  =-1.0


	BM = 0.0
	DO I = 1,4
	BM(1 ,6*I-5)  =  HR(I)
	BM(2 ,6*I-5)  =  HS(I)
	BM(3 ,6*I-4)  =  HR(I)
	BM(4 ,6*I-4)  =  HS(I)
	BM(5 ,6*I-3)  =  HR(I)
	BM(6 ,6*I-3)  =  HS(I)
	BM(7 ,6*I-2)  =  H(I)
	BM(8 ,6*I-2)  =  HR(I)
	BM(9 ,6*I-2)  =  HS(I)
	BM(10,6*I-1)  =  H(I)
	BM(11,6*I-1)  =  HR(I)
	BM(12,6*I-1)  =  HS(I)
	ENDDO


	BR = MATMUL(TEPS,MATMUL(BM,TRANSPOSE(TM)))
	BA = MATMUL(FJ,BR)


	RETURN

      END

C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SHSJ2(BA,BANS,XJI,RN,SN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION BANS(8,24),BA(2,24),BAIN(2,8)
	DIMENSION XJI(2,2)
	

	SN1 = 0.5*(1+RN)
	SN2 = 0.5*(1-RN)
	SN3 = 0.5*(1+SN)
	SN4 = 0.5*(1-SN)

	BAIN = 0.0

	BAIN(1,5) = SN3
	BAIN(1,7) = SN4

C	BAIN(1,1) = SN1 !!
C	BAIN(1,3) = SN2 !!


	BAIN(2,2) = SN1
	BAIN(2,4) = SN2

C	BAIN(2,6) = SN3 !!
C	BAIN(2,8) = SN4 !!


	BA = MATMUL(XJI,MATMUL(BAIN,BANS))

	RETURN

	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHAP4D (R,S,H,P)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
      DIMENSION  H(8),P(2,8)
C
	H = 0.0
	P = 0.0

      RP  = 1.0+R
      SP  = 1.0+S
      RM  = 1.0-R
      SM  = 1.0-S
      R2  = 1.0-R*R
      S2  = 1.0-S*S

      H(1)   = 0.25*RP*SP
      H(2)   = 0.25*RM*SP
      H(3)   = 0.25*RM*SM
      H(4)   = 0.25*RP*SM
      P(1,1) = 0.25*SP
      P(1,2) = -P(1,1)
      P(1,3) = -0.25*SM
      P(1,4) = -P(1,3)
      P(2,1) = 0.25*RP
      P(2,2) = 0.25*RM
      P(2,3) = -P(2,2)
      P(2,4) = -P(2,1)

	RETURN

	END
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE SDID1(BA1,BA2,FJ,XJI,H,HD,TM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION FJ(2,2),XJI(2,2),HD(2,8),HR(8),HS(8),H(8)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION BM(12,24),TM(24,24),TEPS(2,12)
	DIMENSION BA1(2,24),BA2(2,24),BR(2,24)

	HR = 0.0
	HS = 0.0
      DO I = 1,4
      HR(I) = XJI(1,1)*HD(1,I) + XJI(1,2)*HD(2,I)
      HS(I) = XJI(2,1)*HD(1,I) + XJI(2,2)*HD(2,I)
	ENDDO
	


	BM = 0.0
	DO I = 1,4
	BM(1 ,6*I-5)  =  HR(I)
	BM(2 ,6*I-5)  =  HS(I)
	BM(3 ,6*I-4)  =  HR(I)
	BM(4 ,6*I-4)  =  HS(I)
	BM(5 ,6*I-3)  =  HR(I)
	BM(6 ,6*I-3)  =  HS(I)
	BM(7 ,6*I-2)  =  H(I)
	BM(8 ,6*I-2)  =  HR(I)
	BM(9 ,6*I-2)  =  HS(I)
	BM(10,6*I-1)  =  H(I)
	BM(11,6*I-1)  =  HR(I)
	BM(12,6*I-1)  =  HS(I)
	ENDDO


	TEPS = 0.0
	TEPS(1,5)  = 1.0
	TEPS(1,10) = 0.0
	TEPS(2,6)  = 1.0
	TEPS(2,7)  = 0.0

	BR  = MATMUL(TEPS,MATMUL(BM,TRANSPOSE(TM)))
	BA1 = MATMUL(FJ,BR)


	TEPS = 0.0
	TEPS(1,5)  = 0.0
	TEPS(1,10) = 1.0
	TEPS(2,6)  = 0.0
	TEPS(2,7)  =-1.0

	BR  = MATMUL(TEPS,MATMUL(BM,TRANSPOSE(TM)))
	BA2 = MATMUL(FJ,BR)


	RETURN


      END

C	=====================================================================
C	=====================================================================
      SUBROUTINE SDID2(BR,BDID1,BDID2,XJI,H,HD,HR,HS,RN,SN,TM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION XJI(2,2),HD(2,8),HR(8),HS(8),BM(12,24),H(8),TM(24,24)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION BR(12,24),BDID1(8,24),BDID2(8,24)

	DIMENSION BA(2,24),BAIN(2,8)
	
	HR = 0.0
	HS = 0.0
      DO I = 1,4
      HR(I)=XJI(1,1)*HD(1,I)+XJI(1,2)*HD(2,I)
      HS(I)=XJI(2,1)*HD(1,I)+XJI(2,2)*HD(2,I)
	ENDDO

	BM = 0.0
	DO I = 1,4
	BM(1,6*I-5)  =  HR(I)
	BM(2,6*I-5)  =  HS(I)
	BM(3,6*I-4)  =  HR(I)
	BM(4,6*I-4)  =  HS(I)
	BM(5,6*I-3)  =  HR(I)
	BM(6,6*I-3)  =  HS(I)
	BM(7,6*I-2)  =  H(I)
	BM(8,6*I-2)  =  HR(I)
	BM(9,6*I-2)  =  HS(I)
	BM(10,6*I-1) =  H(I)
	BM(11,6*I-1) =  HR(I)
	BM(12,6*I-1) =  HS(I)
	ENDDO


	BR = MATMUL(BM,TRANSPOSE(TM))


	

	SN1 = 0.5*(1+RN)
	SN2 = 0.5*(1-RN)
	SN3 = 0.5*(1+SN)
	SN4 = 0.5*(1-SN)

	BAIN = 0.0

	BAIN(1,5) = SN3
	BAIN(1,7) = SN4

C	BAIN(1,1) = SN1 !!
C	BAIN(1,3) = SN2 !!


	BAIN(2,2) = SN1
	BAIN(2,4) = SN2

C	BAIN(2,6) = SN3 !!
C	BAIN(2,8) = SN4 !!


	BA = MATMUL(XJI,MATMUL(BAIN,BDID1))

	DO J = 1,24
	BR(5,J) = BA(1,J)
	BR(6,J) = BA(2,J)
	ENDDO


	BA = MATMUL(XJI,MATMUL(BAIN,BDID2))

	DO J = 1,24
	BR(10,J) = BA(1,J)
	BR(7,J)  = BA(2,J)
	ENDDO


	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE TRANM1 (VR,VS,VT,TM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION TM(24,24)
	DIMENSION VR(3),VS(3),VT(3)

	TM = 0.0

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


	RETURN

	END

C	=====================================================================
C	=====================================================================
      SUBROUTINE SHJACB (NNO,COORD,HD,VR,VS,VT,DET,FJ,XJI)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
      DIMENSION COORD(3,NNO),HD(2,8),VR(3),VS(3),VT(3)
      DIMENSION COVR(3),COVS(3),RV(3),SV(3),FJ(4),XJI(4)
C
      CALL CLEARA (COVR,3)
      CALL CLEARA (COVS,3)

      DO 20 I=1,NNO
      DO 20 J=1,3
      COVR(J)=COVR(J)+HD(1,I)*(COORD(J,I))
   20 COVS(J)=COVS(J)+HD(2,I)*(COORD(J,I))

      CALL VECPRD (COVR,COVS,VT)
      CALL SCALEN (VT,VT,DET,3)
      CALL SCALEN (COVR,RV,RL,3)
      CALL SCALEN (COVS,SV,SL,3)
      CALL VECPRD (VT,RV,VS)
      CALL ADDVEC (SV,VS,VS)
      CALL SCALEN (VS,VS,DM,3)
      CALL VECPRD (VS,VT,VR)


      CALL SCAPRD (COVR,VR,F1,3)
      CALL SCAPRD (COVS,VR,F2,3)
      CALL SCAPRD (COVR,VS,F3,3)
      CALL SCAPRD (COVS,VS,F4,3)

C	---------------
C	 [J] = [F1  F3]
C		   [F2  F4]
C	---------------
	DET  = F1*F4 - F2*F3

      FJ(1) = F1
	FJ(2) = F2
	FJ(3) = F3
	FJ(4) = F4

	XJI(1) = F4/DET
      XJI(2) =-F2/DET
      XJI(3) =-F3/DET
      XJI(4) = F1/DET


      RETURN

      END

C	=====================================================================
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE TNEAS (FJO,TTI)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION FJO(2,2),TTO(8,8),TTI(8,8)

      !GOTO 10
      
C	FROM Z-SOIL BOOK
	TTO = 0.0

	DO I = 1,8
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
10    TTO = 0.0

	DO I = 1,8
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


	CALL INVMATRIX(TTO,TTI,8)


	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE MTEAS(RN,SN,TTO,DETO,DET,GE,MM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION TTO(8,8),GE(8,MM),TTI(8,8),GEO(8,MM),C(8,MM)

	CALL EASGE(RN,SN,GEO,MM)

CC	DO I = 1,3
CC	DO J = 1,3
CC	TTI(I,J) = DETO*TTO(I,J)/DET
CC	ENDDO
CC    ENDDO
      FACT = DETO/DET
      TTI(1:3,1:3) = FACT*TTO(1:3,1:3)
      
      GE = MATMUL(TTI,GEO)

	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE EASGE(RN,SN,GE,MM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION GE(8,MM)

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
C	=====================================================================
C	=====================================================================
	SUBROUTINE DRILL(RN,SN,BDL,XJI,TM,HD)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION BDL(8,24),HU(4),HV(4),PU(2,4),PV(2,4),XJI(2,2)
	DIMENSION TM(24,24),PR(2,4),PS(2,4),BDLT(8,24)
	DIMENSION HD(2,8),HR(8),HS(8)
	
	HR = 0.0
	HS = 0.0
      DO I = 1,4
      HR(I)=XJI(1,1)*HD(1,I)+XJI(1,2)*HD(2,I)
      HS(I)=XJI(2,1)*HD(1,I)+XJI(2,2)*HD(2,I)
	ENDDO



	BDLT = 0.0

	HU(1) =  - 1.0 - RN - SN - RN*SN + RN*RN + RN*RN*RN
	1		 + RN*RN*SN + RN*RN*RN*SN

	HU(2) =  - 1.0 - RN + SN + RN*SN + RN*RN + RN*RN*RN
	1		 - RN*RN*SN - RN*RN*RN*SN

	HU(3) =  + 1.0 - RN - SN + RN*SN - RN*RN + RN*RN*RN
	1		 + RN*RN*SN - RN*RN*RN*SN

	HU(4) =  + 1.0 - RN + SN - RN*SN - RN*RN + RN*RN*RN
	1		 - RN*RN*SN + RN*RN*RN*SN

	PU(1,1) =  - 1.0 - SN + 2.0*RN + 3.0*RN*RN + 2.0*RN*SN 
	1		   + 3.0*RN*RN*SN

	PU(1,2) =  - 1.0 + SN + 2.0*RN + 3.0*RN*RN - 2.0*RN*SN 
	1		   - 3.0*RN*RN*SN

	PU(1,3) =  + 1.0 + SN - 2.0*RN + 3.0*RN*RN + 2.0*RN*SN 
	1		   - 3.0*RN*RN*SN

	PU(1,4) =  + 1.0 - SN - 2.0*RN + 3.0*RN*RN - 2.0*RN*SN 
	1		   + 3.0*RN*RN*SN



	PU(2,1) = - 1.0 - RN + RN*RN + RN*RN*RN

	PU(2,2) = + 1.0 + RN - RN*RN - RN*RN*RN

	PU(2,3) = - 1.0 + RN + RN*RN - RN*RN*RN

	PU(2,4) = + 1.0 - RN - RN*RN + RN*RN*RN



	HV(1) =  + 1.0 + RN + SN + RN*SN - SN*SN - RN*SN*SN
	1		 - SN*SN*SN - RN*SN*SN*SN

	HV(2) =  - 1.0 - RN + SN + RN*SN + SN*SN + RN*SN*SN
	1		 - SN*SN*SN - RN*SN*SN*SN

	HV(3) =  - 1.0 + RN + SN + RN*SN - SN*SN + RN*SN*SN
	1		 - SN*SN*SN + RN*SN*SN*SN

	HV(4) =  + 1.0 - RN + SN - RN*SN - SN*SN + RN*SN*SN
	1		 - SN*SN*SN + RN*SN*SN*SN



	PV(1,1) = + 1.0 + SN - SN*SN - SN*SN*SN

	PV(1,2) = - 1.0 + SN + SN*SN - SN*SN*SN

	PV(1,3) = + 1.0 - SN - SN*SN + SN*SN*SN

	PV(1,4) = - 1.0 - SN + SN*SN + SN*SN*SN


	PV(2,1) =  + 1.0 + RN - 2.0*SN - 2.0*RN*SN - 3.0*RN*RN 
	1		   - 3.0*RN*SN*SN

	PV(2,2) =  + 1.0 + RN + 2.0*SN + 2.0*RN*SN - 3.0*RN*RN 
	1		   - 3.0*RN*SN*SN

	PV(2,3) =  + 1.0 - RN + 2.0*SN - 2.0*RN*SN - 3.0*RN*RN 
	1		   + 3.0*RN*SN*SN

	PV(2,4) =  + 1.0 - RN - 2.0*SN + 2.0*RN*SN - 3.0*RN*RN 
	1		   + 3.0*RN*SN*SN


	DO I = 1,4 
	HU(I) = HU(I)/8.0
	HV(I) = HV(I)/8.0
	PU(1,I) = PU(1,I)/8.0
	PU(2,I) = PU(2,I)/8.0
	PV(1,I) = PV(1,I)/8.0
	PV(2,I) = PV(2,I)/8.0
	ENDDO


	PR = MATMUL(XJI,PU)
	PS = MATMUL(XJI,PV)


	DO I = 1,4
	BDLT(1,6*I-0) = PR(1,I)
	BDLT(2,6*I-0) = PS(2,I)
	BDLT(3,6*I-0) = PR(2,I) + PS(1,I)
	ENDDO


	BDL = MATMUL(BDLT,TRANSPOSE(TM))

	


	RETURN

	END

C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE SHKJ2(S,DR,BMB,BDL)
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	=====================================================================

      DIMENSION SK(24,24),DR(8,8),BMB(8,24),BDL(8,24),S(1176),
	1		  SK1(24,24),BBRL(24,1),SKK(24,24)
C
      SK  = 0.0
	SK1 = 0.0


	SK1 =  MATMUL(TRANSPOSE(BDL),MATMUL(DR,BDL))
	SK  =  MATMUL(TRANSPOSE(BMB),MATMUL(DR,BMB))

	DO I = 1,24
	DO J = 1,24
	SKK(I,J) = SK(I,J)  + SK1(I,J)
	ENDDO
	ENDDO


	K = 0
	DO I = 1,24
	DO J = I,24
	K = K+1
	S(K) = S(K) + SKK(I,J)
	ENDDO
	ENDDO


      RETURN

      END
C	================================================================
C	================================================================
C	================================================================

