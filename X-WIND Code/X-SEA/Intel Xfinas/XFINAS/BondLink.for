C	=====================================================================
C	=====================================================================
C	START 3D BOND LINKAGE ELEMENT
C	=====================================================================
C	CREATED BY SARA, DEC 20,2004
C	=====================================================================
      SUBROUTINE BONDLNK(PROPM,PROPG,NODEX,WABND,S,COORD,
	1				   EDIS,EDISI,RE,MWG,FIN,MSET)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     ELEMENT PROPERTIES
C     ------------------
C	HN(3,6)		TRANSFORMATION MATRIX
C	DSLIP(3,1)	LOCAL SLIP 
C	RE(6)		GLOBAL INTERNAL RESISTING FORCE
C	SL(3,3)		LOCAL LINEAR STIFFNESS MATRIX
C	S1(6,6)		GLOBAL STIFFNESS MATRIX
C     ------------------------------------------------------------- 
      
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,KONEQ,NIREF,ITOPT,
     2			  ICONV,NOLIN,KSTEP,LIMEQ(2),ITEMAX,NUMREF,
     3			  NUMITE,ITETOT

C	=======================================================================================
C	BONDST READS 4 COLS IN INPUT DATA: FC, D, SMAX & FALG
C	STEND READS 5 COLS IN INPUT DATA: ELEMENT NO, FIRST NODE, MIDDLE NODE, LAST NODE & FALG
C	=======================================================================================

      DIMENSION EDIS(6),S1(6,6),S(21),COODEX(9),RE(6)
      DIMENSION COORD(6),SLIP(3),TSLIP(3),FIN(NEF),SGASH(3,6)
	DIMENSION HN(3,6),DSLIP(3),SIGR(3),KXYZ(3)
	DIMENSION SL(3,3),SIG(3),TSIG(3),SLPCR(3)
	DIMENSION PROPM(*),PROPG(*),SK(3)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION VRi(3),VSi(3),VTi(3)
	DIMENSION VRj(3),VSj(3),VTj(3)
	DIMENSION WABND(1)


	PII = 3.141569


c	WRITE(*,*) PROPG(2),PROPG(3),PROPG(4)	
C	WRITE(*,*) PROPG(5),PROPG(6),PROPG(7)
C	WRITE(*,*) PROPG(8),PROPG(9),PROPG(10)
c	PAUSE

	FGCOR = PROPG(11)
	IF(FGCOR.EQ.0.0) THEN
	COODEX(1) = PROPG(2)
	COODEX(2) = PROPG(3)
	COODEX(3) = PROPG(4)
	COODEX(4) = PROPG(5)
	COODEX(5) = PROPG(6)	
	COODEX(6) = PROPG(7)
	COODEX(7) = PROPG(8)
	COODEX(8) = PROPG(9)
	COODEX(9) = PROPG(10)
	ELSE
	COODEX(1) = PROPG(2)
	COODEX(2) = PROPG(3)
	COODEX(3) = PROPG(4)
	COODEX(4) = PROPG(5)
	COODEX(5) = PROPG(6)	
	COODEX(6) = PROPG(7)
	COODEX(7) = 0.0
	COODEX(8) = 0.0
	COODEX(9) = 0.0
	ENDIF

C	--------------------------------
C	INITIALIZATION OF SOME VARIABLES
C	--------------------------------
	DO 10 I = 1,3
	DO 10 J = 1,3
	SL(I,J) = 0.0
10	CONTINUE

C	-------------------------------------------
C	TO COMPUTE THE LENGTH AND DIRECTION COSINES
C	-------------------------------------------
C	YI=COODEX(1), YJ=COODEX(4), YK=COODEX(7)
C	XI=COODEX(2), XJ=COODEX(5), XK=COODEX(8)
C	ZI=COODEX(3), ZJ=COODEX(6), ZK=COODEX(9)
C	-------------------------------------------
C	A2 = (COODEX(4) - COODEX(1))
C	A1 = (COODEX(5) - COODEX(2))
C	A3 = (COODEX(6) - COODEX(3))
	
C	A5 = (COODEX(7) - COODEX(4))
C	A4 = (COODEX(8) - COODEX(5))
C	A6 = (COODEX(9) - COODEX(6))


c 	COODEX(1) = 10.0
c	COODEX(2) = 0.0
c	COODEX(3) = 0.0
c	COODEX(4) = 10.0
c	COODEX(5) = 0.05	
c	COODEX(6) = 0.0
c	COODEX(7) = 12.0
c	COODEX(8) = 0.0
c	COODEX(9) = 0.0
	     
	A1 = (COODEX(4) - COODEX(1))
	A2 = (COODEX(5) - COODEX(2))
	A3 = (COODEX(6) - COODEX(3))
	
	A4 = (COODEX(7) - COODEX(4))
	A5 = (COODEX(8) - COODEX(5))
	A6 = (COODEX(9) - COODEX(6))

	WL1 = SQRT(A1**2 + A2**2 + A3**2)
	WL2 = SQRT(A4**2 + A5**2 + A6**2)

	VRi(1) = A1/WL1
	VRi(2) = A2/WL1
	VRi(3) = A3/WL1

	VRj(1) = A4/WL2
	VRj(2) = A5/WL2
	VRj(3) = A6/WL2

	CALL FMVEVR (VRi,VSi,VTi)
	CALL FMVEVR (VRj,VSj,VTj)


C	--------------------------------------------------
C	TRANSFORMATION MATRIX HN(3*6) FROM GLOBAL TO LOCAL
C	--------------------------------------------------
	IF(FGCOR.EQ.0.0) THEN 

	WL = (WL1+WL2)/2.0

	DO I = 1,3
	VR(I) = 0.5*(VRi(I) + VRj(I))
	VS(I) = 0.5*(VSi(I) + VSj(I))
	VT(I) = 0.5*(VTi(I) + VTj(I))
	ENDDO

	ELSE

	WL = WL1/2.0

	DO I = 1,3
	VR(I) = VRi(I)
	VS(I) = VSi(I)
	VT(I) = VTi(I)
	ENDDO

	ENDIF

	HN = 0.0
	HN(1,1) = -VR(1)
	HN(1,2) = -VR(2)
	HN(1,3) = -VR(3)
	HN(1,4) =  VR(1)
	HN(1,5) =  VR(2)
	HN(1,6) =  VR(3)
	HN(2,1) = -VS(1)
	HN(2,2) = -VS(2)
	HN(2,3) = -VS(3)
	HN(2,4) =  VS(1)
	HN(2,5) =  VS(2)
	HN(2,6) =  VS(3)
	HN(3,1) = -VT(1)
	HN(3,2) = -VT(2)
	HN(3,3) = -VT(3)
	HN(3,4) =  VT(1)
	HN(3,5) =  VT(2)
	HN(3,6) =  VT(3)


C	DO I =1,3
C	WRITE(*,111) (HN(I,J),J=1,6)
C	WRITE(*,*) MEL,CU1,CU2
C	ENDDO
C111	FORMAT(6E17.5)
C	PAUSE

C	-----------------------
C	INITIALIZE SOME VALUES
C	-----------------------
	SLIP(1) = WABND(1)
	SLIP(2) = WABND(2)
	SLIP(3) = WABND(3)
	SIG(1)  = WABND(4)
	SIG(2)  = WABND(5)
	SIG(3)  = WABND(6)

C	-----------------------------------------------------------
C	CALCULATE THE SLIP FROM THE NODAL DISPLACEMENTS
C	-----------------------------------------------------------
C	{Si}=[A]*{Ui}
C	-----------------------------------------------------------

	SLPCR(1)=HN(1,1)*EDIS(1)+HN(1,2)*EDIS(2)+HN(1,3)*EDIS(3)
	1	+HN(1,4)*EDIS(4)+HN(1,5)*EDIS(5)+HN(1,6)*EDIS(6)
	SLPCR(2)=HN(2,1)*EDIS(1)+HN(2,2)*EDIS(2)+HN(2,3)*EDIS(3)
	1	+HN(2,4)*EDIS(4)+HN(2,5)*EDIS(5)+HN(2,6)*EDIS(6)
	SLPCR(3)=HN(3,1)*EDIS(1)+HN(3,2)*EDIS(2)+HN(3,3)*EDIS(3)
	1	+HN(3,4)*EDIS(4)+HN(3,5)*EDIS(5)+HN(3,6)*EDIS(6)

	DO 30 I = 1,3
	DSLIP(I) = 0.0
30	DSLIP(I) = SLPCR(I) - SLIP(I)

C	--------------------------------------
C	CALCULATE TOTAL SLIP AT ITERATION I
C	--------------------------------------
	DO 40 I = 1,3
	TSLIP(I) = 0.0
40	TSLIP(I) = SLIP(I) + DSLIP(I)

C	------------------------
C	INITIALIZING FC',D ,SMAX
C	------------------------

	SK(1)   = PROPM(6)
	SK(2)   = PROPM(7)
	SK(3)   = PROPM(8) 
	
	FC    = PROPM(9)
	D     = PROPM(10)
	SMAX  = PROPM(11)

C	----------------------------------------------
C	COMPUTE BOND STIFFNESS AND BOND STRESSES
C	----------------------------------------------
	DO 60 I = 1,3

	KXYZ(I) = 0.0D0
	TSIG(I) = 0.0D0

	SELECTCASE(MTMOD)

	CASE(1) 
C	----------------------------------
C	LINEAR FORMULATION 
C	----------------------------------
	KXYZ(I) = SK(I)
	TSIG(I) = SK(I)*TSLIP(I)

	DVOL = 1.0

	CASE(3)
C	----------------------------------
C	FORMULATION FROM BOOK OF MAEKAWA
C	----------------------------------
	IF(ABS(TSLIP(I)).GT.0.0000000583) THEN

	DUMM1 = ABS(TSLIP(I))
	DUMM2 = FC**(2.0/3.0)
	DUMM3 = DUMM1/D
	DUMM4 = -40.0*(DUMM3)**0.6
	DUMM5 = D*(DUMM3**0.4)
	DUMM6 = EXP(DUMM4)
	KXYZ(I) = 21.6*DUMM2*DUMM6/DUMM5               
	IF (ABS(TSLIP(I)).GE.SMAX) THEN
	TSIG(I)=0.0001 
	ELSE
	TSIG(I) = (DUMM1/TSLIP(I))*0.9*DUMM2*(1.0-DUMM6)
	ENDIF
	TSIG(I) = TSIG(I)

	ENDIF

	DVOL =  PII*D*WL

	CASE(4)
C	-----------------------------------------
C	FORMULATION OF NILSON (BECAREFUL FOR UNIT)
C	-----------------------------------------
	DUMM1 = DUMM1/25.4
	DUMM2 = (3606.0E+3)*DUMM1
	DUMM3 = (-5356.0E+6)*DUMM1*DUMM1
	DUMM4 = (1986.0E+9)*DUMM1*DUMM1*DUMM1
	DUMM5 = 3606.0E+3
	DUMM6 = (-2.0*5356.0E+6)*DUMM1
	DUMM7 = (3.0*1986.0E+9)*DUMM1*DUMM1
	KXYZ(I) = (DUMM5 + DUMM6 + DUMM7)*2.7154E-4
	IF(DUMM1.GE.450.0E-6) KXYZ(I) = 0.0

	IF (ABS(TSLIP(I)).GE.SMAX) THEN
	TSIG(I)=0.0001 
	ELSE
	DUMM8 = (DUMM2 + DUMM3 + DUMM4)
	IF(DUMM1.GE.450.0E-6) DUMM8 = 650.0
	TSIG(I) = (DUMM1/(TSLIP(I)/25.4))*DUMM8
	TSIG(I) = TSIG(I)/145.0
	ENDIF

	DVOL =  PII*D*WL

	ENDSELECT

60	CONTINUE


C	--------------------------------------
C	DEFINE ELEMENT LOCAL STIFFNESS SL(3*3)
C	--------------------------------------
	DO 65 I = 1,3
	IF(KXYZ(I).LE.1.0) THEN
	KXYZ(I) = 1.0
	ENDIF
65	CONTINUE

	IF (IFREF) 760,750,760
750	SL(1,1) = KXYZ(1)
	SL(1,2) = 0.0
	SL(1,3) = 0.0
	SL(2,1) = 0.0
	SL(2,2) = KXYZ(2)
	SL(2,3) = 0.0
	SL(3,1) = 0.0
	SL(3,2) = 0.0
	SL(3,3) = KXYZ(3)

C	----------------------------------
C	TRANSFER LOCAL STIFFNESS TO GLOBAL
C	----------------------------------
	DO 110 I = 1,3
	DO 110 J = 1,6
	SGASH(I,J) = 0.0
	DO 110 K = 1,3
110	SGASH(I,J) = SGASH(I,J) + SL(I,K)*HN(K,J)

	DO 120 I = 1,6
	DO 120 J = 1,6
	S1(I,J) = 0.0
	DO 120 K = 1,3
120	S1(I,J) = S1(I,J) + HN(K,I)*SGASH(K,J)

C	--------------------------------------------------
C	COMPUTE FULL STIFFNESS MATRIX S(21) UPPER-TRIANGLE
C	--------------------------------------------------
	KK = 1
	DO 700 I = 1,6
	DO 700 J = I,6
	S(KK) = DVOL*S1(I,J)
	KK = KK+1
700	CONTINUE

C	----------------------
C	COMPUTE FORCES RE(6*1)
C	----------------------
760	DO 100 I = 1,6
	RE(I) = 0.0
	DO 100 J = 1,3
100	RE(I) = RE(I) + DVOL*HN(J,I)*TSIG(J)

	WABND(1) = SLPCR(1) 
	WABND(2) = SLPCR(2)
	WABND(3) = SLPCR(3)
	WABND(4) = TSIG(1)
	WABND(5) = TSIG(2) 
	WABND(6) = TSIG(3)
	

	IF (ITASK.EQ.3) THEN
	DO 800 I = 1,NEF
	FIN(I) = RE(I)
800	CONTINUE
	ENDIF


70	RETURN


	END
C
C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE	TRUSPLS(STSPS,DSTSZ,YOUNG,YIELD,YOUN2,EPSTN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	ELASTIC & PLASTIC BEHAVIOR OF REINFORCING STEEL
C	TRUSS MODIFIED FOR HARDENING-SARA MARCH2005
C	==================================================================


	HARDS = YOUNG*YOUN2/(YOUNG-YOUN2)
	STREP = DSTSZ
	STRPP = STSPS
	STOTP = STRPP + STREP

	FUNSO = STRPP - HARDS*EPSTN
	FUNSN = STOTP - HARDS*EPSTN

	PREYS = YIELD
	IF(ABS(FUNSO).GE.PREYS) GO TO 20
	ESCUR = ABS(FUNSN) - PREYS
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(STREP)
	GO TO 30

20	IF(ABS(FUNSN).LT.ABS(FUNSO)) GO TO 40

	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	STRPP = STRPP + REDUC*STREP + 
	1		RFACT*STREP*(1.0-YOUNG/(YOUNG+HARDS))
	EPSTN = EPSTN + RFACT*STREP/(YOUNG+HARDS)
	EPSTN = EPSTN
	GO TO 50

40	STRPP = STRPP + STREP
50	STSPS = STRPP

	RETURN
	END


C	==================================================================
C	==================================================================
C	==================================================================
	SUBROUTINE	BONGEN(COORD,MONCT,VECTR,MATVC,PROPM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     =================================================================


	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

	COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT


	DIMENSION COORD(NSN,NSC),MONCT(NEF,NELE)
 	DIMENSION VECTR(NGP,NGPS),MATVC(NELE)

	DIMENSION NMDIC(NELE,5),NSTOR(NELE,NSN+1)
	DIMENSION DIREC(NELE,NSN),RELEN(NELE,NSN)
	DIMENSION PROPM(NMP,NMPS)


	DO 600 IELE = 1,NELE

	IGPS = MATVC(IELE)
	VR1  = VECTR(2,IGPS)
	VR2  = VECTR(3,IGPS)
	VR3  = VECTR(4,IGPS)

	
	NNOD1 = MONCT(1,IELE)
	NNOD  = MONCT(2,IELE)

	X1 = COORD(NNOD,1)
	Y1 = COORD(NNOD,2)
	Z1 = COORD(NNOD,3)

	ICOUT = 0

	DO 500 INNO = 1,NSN

	X2 = COORD(INNO,1)
	Y2 = COORD(INNO,2)
	Z2 = COORD(INNO,3)

	VEC1 = (X2-X1)
	VEC2 = (Y2-Y1)
	VEC3 = (Z2-Z1)

	TEST = ABS(VEC1) + ABS(VEC2) + ABS(VEC3)   	

	IF(INNO.EQ.NNOD1) GO TO 500

	IF(TEST.NE.0.0) GO TO 400  !SAME POINT

	ICOUT = ICOUT + 1
	NSTOR(IELE,ICOUT) = INNO
	RELEN(IELE,ICOUT) = 0.0
	DIREC(IELE,ICOUT) = 0.0

	GO TO 500

400	CONTINUE


	CALL SCABED(VR1,VR2,VR3,VRL)
	CALL SCABED(VEC1,VEC2,VEC3,VECL)

C	DOT PRODUCT TO FINE THE NODES LAYING ON THE LINE OF VECTOR
	DOTP = VR1*VEC1 + VR2*VEC2 + VR3*VEC3   
	ADOT = ABS(DOTP)

	
	IF(ADOT.GE.0.99) THEN   !DIFFERENT POINT 
	ICOUT = ICOUT + 1
	NSTOR(IELE,ICOUT) = INNO
	RELEN(IELE,ICOUT) = VECL
	DIREC(IELE,ICOUT) = DOTP
	ENDIF


500	CONTINUE


	NSTOR(IELE,NSN+1) = ICOUT

600	CONTINUE


	DO 700 IELE = 1,NELE

	NCOUT = NSTOR(IELE,NSN+1)



	VLN = 0.0
	DO ICOUT = 1,NCOUT
	IF(RELEN(IELE,ICOUT).GE.VLN) THEN
	VLN = RELEN(IELE,ICOUT)
	NUM = NSTOR(IELE,ICOUT)
	ENDIF
	ENDDO

	VL1 = VLN
	DO ICOUT = 1,NCOUT
	IF(RELEN(IELE,ICOUT).LE.VL1) THEN
	VL1  = RELEN(IELE,ICOUT)
	NUM1 = NSTOR(IELE,ICOUT)
	ENDIF
	ENDDO


	VL2 = VLN
	ICN2 = 0
	DO ICOUT = 1,NCOUT
	IF(DIREC(IELE,ICOUT).LT.0.0.AND.RELEN(IELE,ICOUT).LE.VL2) THEN
	VL2  = RELEN(IELE,ICOUT)
	NUM2 = NSTOR(IELE,ICOUT)
	ICN2 = 1   !!!
	ENDIF
	ENDDO
	



	VL3 = VLN
	ICN3 = 0
	DO ICOUT = 1,NCOUT
	IF(DIREC(IELE,ICOUT).GT.0.0.AND.RELEN(IELE,ICOUT).LE.VL3) THEN
	VL3  = RELEN(IELE,ICOUT)
	NUM3 = NSTOR(IELE,ICOUT)
	ICN3 = 1   !!!
	ENDIF
	ENDDO


	IF(ICN2.NE.0.AND.ICN3.NE.0) THEN
	NMDIC(IELE,1) = IELE
	NMDIC(IELE,2) = NUM2
	NMDIC(IELE,3) = NUM1
	NMDIC(IELE,4) = NUM3
	NMDIC(IELE,5) = 0
	ELSEIF(ICN2.EQ.0.AND.ICN3.NE.0) THEN
	NMDIC(IELE,1) = IELE
	NMDIC(IELE,2) = NUM3
	NMDIC(IELE,3) = NUM1
	NMDIC(IELE,4) = 0
	NMDIC(IELE,5) = 1
	ELSEIF(ICN2.NE.0.AND.ICN3.EQ.0) THEN
	NMDIC(IELE,1) = IELE
	NMDIC(IELE,2) = NUM2
	NMDIC(IELE,3) = NUM1
	NMDIC(IELE,4) = 0
	NMDIC(IELE,5) = 1
	ENDIF


700	CONTINUE


	DO IELE = 1,NELE

	NM1 = NMDIC(IELE,2)
	NM2 = NMDIC(IELE,3)
	NM3 = NMDIC(IELE,4)

	IF(NM3.NE.0) THEN
	VECTR(2,IELE)  = COORD(NM1,1)
	VECTR(3,IELE)  = COORD(NM1,2)
	VECTR(4,IELE)  = COORD(NM1,3)
	VECTR(5,IELE)  = COORD(NM2,1)
	VECTR(6,IELE)  = COORD(NM2,2)
	VECTR(7,IELE)  = COORD(NM2,3)
	VECTR(8,IELE)  = COORD(NM3,1)
	VECTR(9,IELE)  = COORD(NM3,2)
	VECTR(10,IELE) = COORD(NM3,3)
	VECTR(11,IELE) = 0.0
	ELSE
	VECTR(2,IELE)  = COORD(NM1,1)
	VECTR(3,IELE)  = COORD(NM1,2)
	VECTR(4,IELE)  = COORD(NM1,3)
	VECTR(5,IELE)  = COORD(NM2,1)
	VECTR(6,IELE)  = COORD(NM2,2)
	VECTR(7,IELE)  = COORD(NM2,3)
	VECTR(8,IELE)  = 0.0
	VECTR(9,IELE)  = 0.0
	VECTR(10,IELE) = 0.0
	VECTR(11,IELE) = 1.0

	ENDIF

	ENDDO


110	FORMAT(5I6)

	RETURN 
	END


C	===========================================================
C	===========================================================
C	===========================================================

	SUBROUTINE	SCABED(V1,V2,V3,VL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	===========================================================
C	===========================================================

	VL = SQRT(V1*V1 + V2*V2 + V3*V3)
	V1 = V1/VL
	V2 = V2/VL
	V3 = V3/VL

	RETURN
	END

C	===========================================================
C	===========================================================
C	===========================================================
	SUBROUTINE	TRUSPLS_ISOTRO(STSPS,DSTSZ,YOUNG,YIELD,YOUN2,EPSTN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C	==================================================================
C	ELASTIC & PLASTIC BEHAVIOR OF REINFORCING STEEL
C	TRUSS MODIFIED FOR HARDENING-SARA MARCH2005
C	==================================================================

	HARDS = YOUNG*YOUN2/(YOUNG-YOUN2)
	STREP = DSTSZ
	STRPP = STSPS
	STOTP = STRPP + STREP

	PREYS = YIELD + HARDS*EPSTN
	IF(ABS(STRPP).GE.PREYS) GO TO 20
	ESCUR = ABS(STOTP) - PREYS
	IF(ESCUR.LE.0.0) GO TO 40
	RFACT = ESCUR/ABS(STREP)
	GO TO 30

20	IF(STOTP.GT.0.0.AND.STREP.LT.0.0) GO TO 40
	IF(STOTP.LT.0.0.AND.STREP.GT.0.0) GO TO 40
	RFACT = 1.0

30	REDUC = 1.0 - RFACT
	STRPP = STRPP + REDUC*STREP + 
	1		RFACT*STREP*(1.0-YOUNG/(YOUNG+HARDS))
	EPSTN = EPSTN + RFACT*STREP/(YOUNG+HARDS)
	EPSTN = ABS(EPSTN)
	GO TO 50

40	STRPP = STRPP + STREP
50	STSPS = STRPP

	RETURN
	END


C	==================================================================
C	==================================================================
C	==================================================================
	