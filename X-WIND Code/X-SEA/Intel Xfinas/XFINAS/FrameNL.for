C	-------------------------------------------------------------------
      SUBROUTINE BFRASPC (PROPM,PROPG,S,COORD,EDIS,EDISI,RE,PROPO,
	1					FIN,ISET,WA,FIXEN,FIXLR,FIXEO,FIXLO,LREAS,
     2					PMATRL)
C	-------------------------------------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------------
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	-------------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG

      COMMON /INOU/  ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1               IFPR(10),IFPL(10)
      COMMON /DYNA/  CDEN,IMASS


C	-------------------------------------------------------------------    
C	-------------------------------------------------------------------
      DIMENSION PROPM(*),PROPG(100),S(105),COORD(3,2),EDIS(14)
      DIMENSION EDISI(14),RE(14),DR(49)
      DIMENSION DISD(9),REDIS(14),COORDI(3,2),SIGR(8)
      DIMENSION EPS(8),VRO(3),VSO(3),VTO(3),VR(3),VS(3),VT(3),EC(3)
     	DIMENSION R(1,2)
	DIMENSION FIN(NEF)
	DIMENSION AKG(14,14),AKL(14,14),BKG1(14,14),BKG(14,14)
	DIMENSION BMGG(14),BWG(10),BPG(10),BMSG(7,14),DRMAT(7,7)
	DIMENSION BBX(14),TT(14,14),EDISL(14)
	DIMENSION WA(1),LSEH(14),FIXEN(NEF),FIXLR(NEF)
	DIMENSION FIXEO(NEF),FIXLO(NEF)
	DIMENSION TRANH(14,14),LREAS(14),TRANO(14,14),PROPO(6)
	DIMENSION COORO(3,2),TRANOL(14,14)

	DIMENSION DSIGM(5,5),QMAT(5,5),FMAT(5,5),TRANM(5,5),TRNFLX(14,5)
	DIMENSION IJ(4)

	DIMENSION PMATRL(NMP,1)


C	CALLING ELEMENt GAUSS POINT NUMBER
	CALL INTFILL('OGRP',NGG ,5 ,KEG,0) !

C	WRITE(*,*) MEL,IFEIG,ISOLOP,ITASK
C	----------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	CALL TRNOOF(TRANO,PROPO)

	IF(NLOPT.EQ.3) THEN
	CALL TRNOFT2(COORD,COORO,EDIS,PROPO)
	ELSE
	CALL TRNOFT1(COORD,COORO,EDIS,PROPO,TRANO)
	ENDIF

	

C	---------------------
C     GEOMETRIC STIFFNESS
C     INTERNAL FORCE VECTOR
C     ---------------------------------------
C     DEFINE GEOMETRIC AND MATERIAL CONSTANTS
C     ---------------------------------------
      CALL FMINIT (EDIS,COORO,COORDI,VSO,NLOPT)


C	CALL SECTION PROPERTIES
	CALL XFSECTION(KEG,MEL,1)

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,ROTAG,1,5 ,0) !
C	ROTAG   = PROPG(14)


      DO 10 I=1,14
   10 REDIS(I)=EDIS(I)

C     --------------------------------------------------------
C     LOCAL DIRECTION COSINES AND NODAL COORDINATES (VR,VS,VT)
C     REMOVES RIGID BODY TRANSLATIONS AND ROTATIONS FROM
C     TOTAL DISPLACEMENT VECTOR (NLOPT>1) 
C     --------------------------------------------------------
	

      CALL BMRIGD (COORO,COORDI,EDIS,REDIS,VSO,VR,VS,VT,NNO,NLOPT,EC,
     1             R,ELN,ROTAG)

C	Vr,Vs,Vt.....ARE THE UNIT VECTORS ALONG THE LOCAL DIRECTIONS r,s,t
C     --------------------------------------------

C	CALLULATE INITIAL LENGTH OF THE BEAM
	ELN = SQRT((COORO(1,2)-COORO(1,1))**2+(COORO(2,2)-
	1		    COORO(2,1))**2+(COORO(3,2)-COORO(3,1))**2)


C	--------------------------------------------------------------
C	OBTAIN RELEASE TRANSFORMATION MATRIX
	CALL TRNHIG(TRANH,ELN,LREAS)
C	--------------------------------------------------------------
C	--------------------------------------------------------------
C	OBTAIN LOCAL OFFSET TRANSFORMATION MATRIX
	CALL TRNOOFL(TRANOL,PROPO,VR,VS,VT)
C	--------------------------------------------------------------

C	CALLING ELEMENT GAUSS DATA
	CALL FRGAUS (NGG,ELN,BPG,BWG)

      
	DO I = 1,105
	S(I)  = 0.0
	ENDDO
	DO I = 1,49
	DR(I)  = 0.0
	ENDDO
	DO I = 1,8
	SIGR(I)  = 0.0
	ENDDO
	DO I = 1,14
	RE(I) = 0.0
	DO J = 1,14
	AKG(I,J) = 0.0
	ENDDO
	ENDDO

C	----------------------------------------------------------
C	TRANSFORMATION TO LOCAL COORDINATE
C	----------------------------------------------------------
	CALL TT1A (VR,VS,VT,TT)
	EDISL = MATMUL(TRANSPOSE(TT),REDIS)

C	TRANSFORM CORRESPONDING RELEASE CONDITION {Ur} = [T]{U}
	CALL TRNMUL(TRANH,EDISL,1)

	TRANM = 0.0D0
	QMAT  = 0.0D0
	FMAT  = 0.0D0
C     ----------------------------------------------------------
C     LOOP OVER GAUSS TO DET STRESS CONT TO ELEMENT FORCE VECTOR
C     ----------------------------------------------------------
      DO 300 II = 1,NGG

C	CALL SECTION PROPERTIES
	CALL XFSECTION(KEG,MEL,II)
	
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,ASS ,1,58,0) !
	CALL MRELFIL(NAMEI,ATT ,1,59,0) !	

C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
	BXD = BPG(II)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
	CALL BBXFRM(BBX,ELN,BXD)

	
C     ------------------------------------------------------
C     ELASTIC RIGIDITY MATRIX (DR)
C     ------------------------------------------------------
      IF(MTMOD.GT.2) GOTO 220
	
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,AE   ,1,51,0) !
	CALL MRELFIL(NAMEI,QSE  ,1,52,0) !
	CALL MRELFIL(NAMEI,QTE  ,1,53,0) !
	CALL MRELFIL(NAMEI,SIE  ,1,54,0) !
	CALL MRELFIL(NAMEI,TIE  ,1,55,0) !
	CALL MRELFIL(NAMEI,SITE ,1,56,0) !
	CALL MRELFIL(NAMEI,PJL  ,1,57,0) !

C	RIGIDITY MATRIX 
	DR(1) = +AE     !AE
	DR(4) = +QSE    !QsE
	DR(5) = -QTE    !QtE
	DR(25)= +SIE    !IsE
	DR(26)= -SITE   !IstE
	DR(33)= +TIE    !ItE
	DR(22)= +QSE    !QsE
	DR(29)= -QTE    !QtE
	DR(32)= -SITE   !IstE
	DR(41)= +PJL    !GJ
	DR(49)= +AE     !DUMMY FOR WARPING STABILITY
	
		
  220 IF (ITASK.NE.5) GOTO 245
C     -------------------------------------------
C     CONSISTENT AND LUMPED MASS MATRIX (ITASK=5)
C     -------------------------------------------

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,AREAD  ,1,41,0) !
	CALL MRELFIL(NAMEI,PJLD   ,1,47,0) !
	STE = PJLD/AREAD   !POLAR/AREA FACTOR
      CALL FRMASS(S,IMASS,AREAD,ELN,VR,VS,VT,STE,GAR)
      GOTO 900

C     -------------------------------------------------------------
C     FIND LINEAR GLOBAL DAMPING MATRIX OF FRAME ELEMENT (IFREF=0)
C     -------------------------------------------------------------
  245 IF (ITASK.NE.6) GOTO 250

	GOTO 900

  250 CONTINUE


C     ------------------------------------------------------------
C     STRAIN TERMS MEMBRANE - BENDING - TORSION
C     ------------------------------------------------------------
	EPS(1)= BBX(1)*EDISL(1)+BBX(8)*EDISL(8)                         !Er
      EPS(4)= BBX(3)*EDISL(3)+BBX(4)*EDISL(5)+BBX(10)*EDISL(10)+
     +		BBX(11)*EDISL(12)                                       !Ps
      EPS(5)= BBX(2)*EDISL(2)+BBX(5)*EDISL(6)+BBX(9)*EDISL(9)+
     +		BBX(12)*EDISL(13)                                       !Pt
	EPS(6)= BBX(6)*EDISL(4)+BBX(13)*EDISL(11)                       !Tor

C     ------------------------------------------------------------
C     IF NLOPT>1 SUBTRACT NONLINEAR STRAIN TERMS (ALMANSI STRAINS)
C     ------------------------------------------------------------
      IF (NLOPT.LE.1) GOTO 400
      EPS(1)=EPS(1) -.5*EPS(1)*EPS(1)
      EPS(2)=EPS(2)
      EPS(3)=EPS(3)


400	CONTINUE

C     ---------------------------------
C     NODAL STRESS - RESULTANTS SIGR(8)
C     ---------------------------------

C     LINEAR ELASTIC
	IF (MTMOD.EQ.1) THEN

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,AE   ,1,51,0) !
	CALL MRELFIL(NAMEI,QSE  ,1,52,0) !
	CALL MRELFIL(NAMEI,QTE  ,1,53,0) !
	CALL MRELFIL(NAMEI,SIE  ,1,54,0) !
	CALL MRELFIL(NAMEI,TIE  ,1,55,0) !
	CALL MRELFIL(NAMEI,SITE ,1,56,0) !
	CALL MRELFIL(NAMEI,PJL  ,1,57,0) !

	SIGR(1) = AE*EPS(1) + QSE*EPS(4) - QTE*EPS(5)   !AXIAL FORCE
	SIGR(4) = QSE*EPS(1)+ SIE*EPS(4) - SITE*EPS(5)  !MOMENT ABOUT S-AXIS (MY FOR XY PLANE)
	SIGR(5) =-QTE*EPS(1)- SITE*EPS(4)+ TIE*EPS(5)   !MOMENT ABOUT T-AXIS (MZ FOR XY PLANE)
	SIGR(6) = PJL*EPS(6)                            !TORSIONAL MOMENT

	ENDIF



C     ELASTO PLASTIC
	IF (MTMOD.EQ.3) THEN
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FNFIB,1,4,0) !
	NFIB = INT(FNFIB)
	CALL FMPLAS(PMATRL,EPS,DR,SIGR,WA,II,NFIB,NMP,ISET,KEG) 
	ENDIF


C     CONCRETE
	IF (MTMOD.EQ.5) THEN
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FNFIB,1,4,0) !
	NFIB = INT(FNFIB)
	CALL SPCONT(EPS,II,DR,SIGR,NMP,PMATRL,
	1			WA(1),WA(8*NGG+1),NLOPT,NFIB,ISET,KEG,NGG)
	ENDIF


	D11 = AE
	D12 = QTE
	D13 = QSE
	D21 = QTE 
	D22 = TIE
	D23 = SITE
	D31 = QSE
	D32 = SITE
	D33 = SIE

	DSIGM = 0.0D0

	DETS = D11*D22*D33+D21*D13*D32+D31*D12*D23
	1	  -D11*D23*D32-D21*D12*D33-D31*D13*D22
	DSIGM(1,1) = (D22*D33-D23*D32)/DETS
	DSIGM(1,2) = (D13*D32-D12*D33)/DETS
	DSIGM(1,3) = (D12*D23-D13*D22)/DETS

	DSIGM(2,1) = (D23*D31-D21*D33)/DETS
	DSIGM(2,2) = (D11*D33-D13*D31)/DETS
	DSIGM(2,3) = (D13*D21-D11*D23)/DETS

	DSIGM(3,1) = (D21*D32-D22*D31)/DETS
	DSIGM(3,2) = (D12*D31-D11*D32)/DETS
	DSIGM(3,3) = (D11*D22-D12*D21)/DETS

	IF(ASS.GT.0.0D0) DSIGM(4,4) = 1.0/ASS
	IF(ATT.GT.0.0D0) DSIGM(5,5) = 1.0/ATT

	DSIGM = DSIGM*BWG(II)
 
	TRANM(1,1) = 1.0D0
	TRANM(2,2) = ELN-BXD
	TRANM(3,3) =-1.0D0*(ELN-BXD)
	TRANM(3,4) = 1.0D0
	TRANM(2,5) = 1.0D0
	TRANM(4,2) = 1.0D0
	TRANM(5,3) = 1.0D0

C	QMAT(1:5,1:5) = 0.0D0
C	QMAT(1:3,1)  = [D11,D21,D31]
C	QMAT(1:3,2)  = [D12,D22,D32]
C	QMAT(1:3,3)  = [D13,D23,D33]
C	QMAT(4,4)    = 1.0/ASS
C	QMAT(5,5)    = 1.0/ATT
C	QMAT = QMAT*BWG(II)

	FMAT = FMAT + MATMUL(TRANSPOSE(TRANM),MATMUL(DSIGM,TRANM))


C	INTERNAL FORCE

C	MEMBRANE
	RE(1)  = RE(1)  + BWG(II)*SIGR(1)*BBX(1) 
	RE(8)  = RE(8)  + BWG(II)*SIGR(1)*BBX(8) 
C	BENDING
      RE(2)  = RE(2)  + BWG(II)*SIGR(5)*BBX(2) 
      RE(3)  = RE(3)  + BWG(II)*SIGR(4)*BBX(3) 
      RE(5)  = RE(5)  + BWG(II)*SIGR(4)*BBX(4) 
      RE(6)  = RE(6)  + BWG(II)*SIGR(5)*BBX(5) 
      RE(9)  = RE(9)  + BWG(II)*SIGR(5)*BBX(9) 
      RE(10) = RE(10) + BWG(II)*SIGR(4)*BBX(10) 
      RE(12) = RE(12) + BWG(II)*SIGR(4)*BBX(11)
      RE(13) = RE(13) + BWG(II)*SIGR(5)*BBX(12)
C	TORSION
	RE(4)  = RE(4)  + BWG(II)*SIGR(6)*BBX(6)
	RE(11) = RE(11) + BWG(II)*SIGR(6)*BBX(13) 


	IDT = 0
	DO I = 1,7
	DO J = 1,7
	IDT = IDT + 1
	DRMAT(I,J) = BWG(II)*DR(IDT)
	ENDDO
	ENDDO
	DO I = 1,7
	DO J = 1,14
	BMSG(I,J) = 0.0
	ENDDO
	ENDDO 


	BMSG(1,1)  =  BBX(1)
	BMSG(1,8)  =  BBX(8)
	BMSG(4,3)  =  BBX(3)
	BMSG(4,5)  =  BBX(4)
	BMSG(4,10) =  BBX(10)
	BMSG(4,12) =  BBX(11)
	BMSG(5,2)  =  BBX(2)
	BMSG(5,6)  =  BBX(5)
	BMSG(5,9)  =  BBX(9)
	BMSG(5,13) =  BBX(12)
	BMSG(6,4)  =  BBX(6)
	BMSG(6,11) =  BBX(13)
	BMSG(7,7)  = -1.0/ELN
	BMSG(7,14) =  1.0/ELN

	AKL = MATMUL(TRANSPOSE(BMSG),MATMUL(DRMAT,BMSG))


	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = AKG(I,J) + AKL(I,J)
	ENDDO
	ENDDO


  300 CONTINUE !END LOOP GAUSS POINT

	
	IF(MTMOD.EQ.1) THEN

	QMAT = 0.0
C	CALL INVMATRIX(FMAT,QMAT,5)
	CALL INVMATF(FMAT,QMAT,5,IB)
	IF(IB.EQ.1) GOTO 305

	TRNFLX = 0.0D0
	TRNFLX(1,1) = -1.0D0
	TRNFLX(2,2) = -1.0D0
	TRNFLX(3,3) = -1.0D0
	TRNFLX(5,4) = -1.0D0
	TRNFLX(6,5) = -1.0D0

	TRNFLX(5,3) =  ELN
	TRNFLX(6,2) = -ELN

	TRNFLX(8 ,1) =  1.0D0
	TRNFLX(9 ,2) =  1.0D0
	TRNFLX(10,3) =  1.0D0
	TRNFLX(12,4) =  1.0D0
	TRNFLX(13,5) =  1.0D0

	AKL = MATMUL(TRNFLX,MATMUL(QMAT,TRANSPOSE(TRNFLX)))

	IJ(1:4) = [4,7,11,14]
	DO II = 1,4
	I = IJ(II)
	DO J = 1,14
	AKL(I,J) = AKG(I,J)
	AKL(J,I) = AKG(J,I)
	ENDDO
	ENDDO
	AKG = AKL


	RE = MATMUL(AKG,EDISL)

305	CONTINUE

	ENDIF


888	CONTINUE
C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
	CALL TRNMUL(TRANH,RE,2)
	CALL TRNMUM(TRANH,AKG)
C	------------------------------------------------------------

	PLN = RE(8)

	IF (IFEIG.EQ.0.AND.ISOLOP.EQ.4) GOTO 800

c	IF (ITASK.LE.3) GOTO 500
c	IF (IFEIG.NE.0) GOTO 900
c	GOTO 800
c500	CONTINUE

	DO IEF = 1,NEF
	FIN(IEF) = RE(IEF)
	ENDDO
C	----------------------------------------------------------
C	TRANSFORMATION DUE TO LOCAL OFFSET
C	----------------------------------------------------------
	CALL TRNMUL(TRANOL,FIN,2)
C	----------------------------------------------------------


C	TRANSFORMATION
	CALL TT1A (VR,VS,VT,TT)
	RE = MATMUL(TT,RE)

C	----------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	CALL TRNMUL(TRANO,RE,2)
C	----------------------------------------------------------

	DO IEF = 1,NEF
	FIXLR(IEF) = 0.0
	FIXLO(IEF) = 0.0
	DO JEF = 1,NEF
	FIXLR(IEF) = FIXLR(IEF) + TT(IEF,JEF)*FIXEN(JEF)  !VARY FIXEND
	FIXLO(IEF) = FIXLO(IEF) + TT(IEF,JEF)*FIXEO(JEF)  !CONT FIXEND
	ENDDO
	ENDDO
C	CALL TRNMUL(TRANOL,FIXLR,1)


	CALL TT1A (VR,VS,VT,TT)
	BKG1 = MATMUL(TT,AKG)
	BKG  = MATMUL(BKG1,TRANSPOSE(TT))

C	----------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	CALL TRNMUM(TRANO,BKG)
C	------------------------------------------------------------

	IDT = 0
	DO I = 1,14
	DO J = I,14
	IDT = IDT + 1
	S(IDT) = S(IDT) + BKG(I,J)
	ENDDO
	ENDDO

C     --------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>1)
C     --------------------------------------------------------
	IF (NLOPT.LE.1) GOTO 820

800	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = 0.0
	ENDDO
	ENDDO

	DO II = 1,NGG

	BXD = BPG(II)
		
	DO I = 1,14
	BMGG(I) = 0.0
	ENDDO

	BMGG(1)  =  -1.0/ELN
	BMGG(8)  =   1.0/ELN
C	BMGG(4)  =  -1.0/ELN
C	BMGG(11) =   1.0/ELN
C	BMGG(7)  =  -1.0/ELN
C	BMGG(14) =   1.0/ELN

	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = AKG(I,J) + BWG(II)*PLN*BMGG(I)*BMGG(J)
	ENDDO
	ENDDO

	DO I = 1,14
	BMGG(I) = 0.0
	ENDDO

	BMGG(3)  =   (6.0*BXD/ELN/ELN) - (6.0*BXD*BXD/ELN/ELN/ELN)
	BMGG(5)  =   1.0-(4.0*BXD/ELN) + (3.0*BXD*BXD/ELN/ELN)
	BMGG(10) =  -(6.0*BXD/ELN/ELN) + (6.0*BXD*BXD/ELN/ELN/ELN)
	BMGG(12) =  -(2.0*BXD/ELN)     + (3.0*BXD*BXD/ELN/ELN)

	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = AKG(I,J) + BWG(II)*PLN*BMGG(I)*BMGG(J)
	ENDDO
	ENDDO

	DO I = 1,14
	BMGG(I) = 0.0
	ENDDO

	BMGG(2)  =  -(6.0*BXD/ELN/ELN) + (6.0*BXD*BXD/ELN/ELN/ELN)
	BMGG(6)  =   1.0-(4.0*BXD/ELN) + (3.0*BXD*BXD/ELN/ELN)
	BMGG(9)  =   (6.0*BXD/ELN/ELN) - (6.0*BXD*BXD/ELN/ELN/ELN)
	BMGG(13) =  -(2.0*BXD/ELN)     + (3.0*BXD*BXD/ELN/ELN)

	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = AKG(I,J) + BWG(II)*PLN*BMGG(I)*BMGG(J)
	ENDDO
	ENDDO


	ENDDO !END II


C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
	CALL TRNMUM(TRANH,AKG)
C	------------------------------------------------------------


	CALL TT1A (VR,VS,VT,TT)
	BKG1 = MATMUL(TT,AKG)
	BKG  = MATMUL(BKG1,TRANSPOSE(TT))

C	----------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	CALL TRNMUM(TRANO,BKG)
C	------------------------------------------------------------

	IDT = 0
	DO I = 1,14
	DO J = I,14
	IDT = IDT + 1
	S(IDT) = S(IDT) + BKG(I,J)
	ENDDO
	ENDDO



820   CONTINUE !TIM(12)=TIM(12)+TIM2


900   CONTINUE


C	WRITE(*,*) 'KAK',MEL,AE,QSE,QTE,SIE,TIE,SITE,PJL 


	RETURN

      END

C     =================================================================
C     =================================================================
C     =================================================================
	SUBROUTINE BBXFRM(BBX,ELN,BXD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION BBX(14)


	F1  =  1.0/ELN                             
      F2  = -6.0*(1.0-2.0*BXD/ELN)/ELN**2.0    
	F3  = -6.0*(1.0-2.0*BXD/ELN)/ELN**2.0
	F4  =  1.0/ELN 
	F5  = -2.0*(2.0-3.0*BXD/ELN)/ELN
	F6  = -2.0*(2.0-3.0*BXD/ELN)/ELN
	F7  =  1.0/ELN
	F8  =  6.0*(1.0-2.0*BXD/ELN)/ELN**2.0 
	F9  =  6.0*(1.0-2.0*BXD/ELN)/ELN**2.0 
	F10 =  1.0/ELN
	F11 = -2.0*(1.0-3.0*BXD/ELN)/ELN
	F12 = -2.0*(1.0-3.0*BXD/ELN)/ELN
		            

	BBX(1) = -F1                               
      BBX(2) = +F2                
      BBX(3) = -F3        
      BBX(6) = -F4      !TORSION                                                          
      BBX(4) = +F5                   
      BBX(5) = +F6               
      BBX(8) = +F7                                
      BBX(9) = +F8               
      BBX(10)= -F9              
      BBX(13)= +F10     !TORSION                                                                   
      BBX(11)= +F11
      BBX(12)= +F12

	RETURN

	END


C	==========================================================================
C	==========================================================================
C	==========================================================================
	SUBROUTINE FMTRES(PST,IPT,FBSTR,AE,QSE,QTE,SIE,TIE,SITE,YNG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	
	DIMENSION PST(2,4)
	DIMENSION FBSTR(28)

	DO IQR = 1,4
	SCOR = PST(1,IQR)
	TCOR = PST(2,IQR)
C	FEPS  = EPS(1) + TCOR*EPS(4)  - SCOR*EPS(5)
	FBSTR(IQR+14*(IPT-1))   = SCOR !YOUNG*FEPS
	FBSTR(IQR+14*(IPT-1)+4) = TCOR 
C	WRITE(*,*) IQR,SCOR,TCOR
C	PAUSE
	ENDDO

	FBSTR(9 +14*(IPT-1)) = AE/YNG
	FBSTR(10+14*(IPT-1)) = QSE/YNG
	FBSTR(11+14*(IPT-1)) = QTE/YNG
	FBSTR(12+14*(IPT-1)) = SIE/YNG
	FBSTR(13+14*(IPT-1)) = TIE/YNG
	FBSTR(14+14*(IPT-1)) = SITE/YNG


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE TRNHIG(TRANH,ELN,LREAS)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION TRANH(14,14),LREAS(14)


	DO I = 1,14
	DO J = 1,14
	TRANH(I,J) = 0.0
	ENDDO
	ENDDO	

	DO I = 1,14
	TRANH(I,I) = 1.0
	ENDDO

	CALL RESCON(LREAS,ELN,TRANH)


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE TRNMUL(TRANH,A,ID)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION TRANH(14,14),A(14),B(14)


	IF(ID.EQ.1) THEN

	DO I = 1,14
	B(I) = 0.0
	DO J = 1,14
	B(I) = B(I) + TRANH(I,J)*A(J)
	ENDDO
	ENDDO	
	

	ELSEIF(ID.EQ.2) THEN

	DO I = 1,14
	B(I) = 0.0
	DO J = 1,14
	B(I) = B(I) + TRANH(J,I)*A(J)
	ENDDO
	ENDDO	
	
	ENDIF

	DO I = 1,14
	A(I) = B(I)
	ENDDO


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE TRNMUM(TRANH,A)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION TRANH(14,14),A(14,14),B(14,14)
	
	B = MATMUL(TRANSPOSE(TRANH),MATMUL(A,TRANH))

	A(1:14,1:14) = B(1:14,1:14)



	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE RESCON(LREAS,EL,B)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION LREAS(14),A(14,14),B(14,14),C(14,14),D(14,14)


	DO I = 1,14	
	DO J = 1,14
	B(I,J) = 0.0
	C(I,J) = 0.0
	D(I,J) = 0.0
	A(I,J) = 0.0
	ENDDO
	ENDDO


	A(1,1)   = 1.0/EL
	A(2,2)   = 12.0/EL/EL/EL
	A(3,3)   = 12.0/EL/EL/EL
	A(4,4)   = 6.0/5.0/EL
	A(5,5)   = 4.0/EL
	A(6,6)   = 4.0/EL
	A(7,7)   = 2.0*EL/15.0

	A(8,8)   = 1.0/EL
	A(9,9)   = 12.0/EL/EL/EL
	A(10,10) = 12.0/EL/EL/EL
	A(11,11) = 6.0/5.0/EL
	A(12,12) = 4.0/EL
	A(13,13) = 4.0/EL
	A(14,14) = 2.0*EL/15.0

	A(1,8)   = -1.0/EL
	A(2,6)   =  6.0/EL/EL
	A(2,9)   = -12.0/EL/EL/EL
	A(2,13)  =  6.0/EL/EL
	A(3,5)   = -6.0/EL/EL
	A(3,10)  = -12.0/EL/EL/EL
	A(3,12)  = -6.0/EL/EL
	A(4,7)   =  1.0/10.0
	A(4,11)  = -6.0/5.0/EL
	A(4,14)  =  1.0/10.0
	A(5,10)  =  6.0/EL/EL
	A(5,12)  =  2.0/EL
	A(6,9)   = -6.0/EL/EL
	A(6,13)  =  2.0/EL
	A(7,14)  = -EL/30.0

	A(9,13)  = -6.0/EL/EL
	A(10,12) =  6.0/EL/EL
	A(11,14) = -1.0/10.0

	DO I = 1,14
	DO J = 1,I
	A(I,J) = A(J,I)
	ENDDO
	ENDDO


	DO I = 1,14	
	K1 = LREAS(I)
	DO J = 1,14
	K2 = LREAS(J)
	IF(K1*K2.NE.0)           B(I,J) = A(I,J)
	IF(K1.NE.0.AND.K2.EQ.0)  C(I,J) =-A(I,J)
	ENDDO
	ENDDO

	DO I = 1,14	
	K1 = LREAS(I)
	IF(K1.EQ.0) B(I,I) = 1.0D0
	IF(K1.EQ.0) C(I,I) = 1.0D0
	ENDDO
	
	N = 14
	CALL INVMATRIX(B,D,N)

      B = MATMUL(D,C)

	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE TRNOOF(TRANO,PROPO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION TRANO(14,14),PROPO(6)

	TRANO = 0.0
	DO I = 1,14
	TRANO(I,I) = 1.0D0
	ENDDO

	TRANO(1,5) = PROPO(3)
	TRANO(1,6) =-PROPO(2)
	TRANO(2,4) =-PROPO(3)
	TRANO(2,6) = PROPO(1)
	TRANO(3,4) = PROPO(2)
	TRANO(3,5) =-PROPO(1)
	

	TRANO(8 ,12) = PROPO(6)
	TRANO(8 ,13) =-PROPO(5)
	TRANO(9 ,11) =-PROPO(6)
	TRANO(9 ,13) = PROPO(4)
	TRANO(10,11) = PROPO(5)
	TRANO(10,12) =-PROPO(4)

C	WRITE(*,*) (PROPO(I),I=1,6)
C	PAUSE


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE TRNOFT1(COORD,COORO,EDIS,PROPO,TRANO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION COORD(3,2),COORO(3,2),EDIS(14),PROPO(6)
	DIMENSION VE(3)
	DIMENSION TRANO(14,14)
	
C     FOR DGMM FUNCTION
	DIMENSION C(14)
	
	EDIS = MATMUL(TRANO,EDIS)
      
	VE(1) = PROPO(1)
	VE(2) = PROPO(2)
	VE(3) = PROPO(3)


	DO I = 1,3
	COORO(I,1) = COORD(I,1) + VE(I)
	ENDDO

C	=================================================================


	VE(1) = PROPO(4)
	VE(2) = PROPO(5)
	VE(3) = PROPO(6)


	DO I = 1,3
	COORO(I,2) = COORD(I,2) + VE(I) 
	ENDDO


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE TRNOFT2(COORD,COORO,EDIS,PROPO)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION COORD(3,2),COORO(3,2),EDIS(14),PROPO(6)
	DIMENSION ROV(3),VE(3),RV(3),V1(3),V2(3),V3(3)

	ROV(1) = EDIS(4)
	ROV(2) = EDIS(5)
	ROV(3) = EDIS(6)

	VE(1) = PROPO(1)
	VE(2) = PROPO(2)
	VE(3) = PROPO(3)


	V1 = 0.0
	V2 = 0.0
	V3 = 0.0


	CALL SCALEN(ROV,RV,ALP,3)

	IF(ALP.EQ.0.0D0) GO TO 100

	ALP2 = ALP*ALP
	SS = SIN(ALP)/ALP
	CC = COS(ALP)
	C1 = (1.0-CC)/ALP2

	CALL VECPRD(ROV,VE,V1)
	CALL VECPRD(ROV,V1,V2)

	DO I = 1,3
	V3(I) = SS*V1(I) + C1*V2(I)
	ENDDO

100	CONTINUE


	DO I = 1,3
	EDIS(I) = EDIS(I) + V3(I)
	ENDDO


	DO I = 1,3
	COORO(I,1) = COORD(I,1) + VE(I) + V3(I)
	ENDDO

C	=================================================================

	ROV(1) = EDIS(11)
	ROV(2) = EDIS(12)
	ROV(3) = EDIS(13)

	VE(1) = PROPO(4)
	VE(2) = PROPO(5)
	VE(3) = PROPO(6)


	CALL SCALEN(ROV,RV,ALP,3)

	IF(ALP.EQ.0.0D0) GO TO 200

	ALP2 = ALP*ALP
	SS = SIN(ALP)/ALP
	CC = COS(ALP)
	C1 = (1.0-CC)/ALP2

	CALL VECPRD(ROV,VE,V1)
	CALL VECPRD(ROV,V1,V2)

	DO I = 1,3
	V3(I) = SS*V1(I) + C1*V2(I)
	ENDDO

200	CONTINUE

	DO I = 1,3
	EDIS(I+7) = EDIS(I+7) + V3(I)
	ENDDO


	DO I = 1,3
	COORO(I,2) = COORD(I,2) + VE(I) + V3(I)
	ENDDO


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE TRNOOFL(TRANO,PROPO,VR,VS,VT)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION TRANO(14,14),PROPO(6),VR(3),VS(3),VT(3)
	DIMENSION VF(3)

	VF(1) = VR(1)*PROPO(1)+VR(2)*PROPO(2)+VR(3)*PROPO(3)
	VF(2) = VS(1)*PROPO(1)+VS(2)*PROPO(2)+VS(3)*PROPO(3)
	VF(3) = VT(1)*PROPO(1)+VT(2)*PROPO(2)+VT(3)*PROPO(3)

	TRANO = 0.0
	DO I = 1,14
	TRANO(I,I) = 1.0D0
	ENDDO

	TRANO(1,5) = VF(3)
	TRANO(1,6) =-VF(2)
	TRANO(2,4) =-VF(3)
	TRANO(2,6) = VF(1)
	TRANO(3,4) = VF(2)
	TRANO(3,5) =-VF(1)
	
	VF(1) = VR(1)*PROPO(4)+VR(2)*PROPO(5)+VR(3)*PROPO(6)
	VF(2) = VS(1)*PROPO(4)+VS(2)*PROPO(5)+VS(3)*PROPO(6)
	VF(3) = VT(1)*PROPO(4)+VT(2)*PROPO(5)+VT(3)*PROPO(6)

	TRANO(8 ,12) = VF(3)
	TRANO(8 ,13) =-VF(2)
	TRANO(9 ,11) =-VF(3)
	TRANO(9 ,13) = VF(1)
	TRANO(10,11) = VF(2)
	TRANO(10,12) =-VF(1)

C	WRITE(*,*) (PROPO(I),I=1,6)
C	PAUSE


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE INVMATF(A,B,N,IND)
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

	EPS=0.0D0
	IND = 0

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
   33 CONTINUE
	IND = 1
	RETURN
	WRITE(*,*) 'SINGULAR MATRIX EPS- INSIDE SUB INVMATRIX'
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


C	=====================================================================
C	=====================================================================
C	=====================================================================
      SUBROUTINE FMINIT (EDIS,COORD,COORDI,VSO,NLOPT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------------------------
C     FINDS NODAL COORDINATES AND DIRECTION COSINES OF S-AXIS FOR
C     ORIGINAL CONFIGURATION AND TRANSFERS GEOMETRIC DATA FROM
C     PROPG TO LOCAL VARIABLES
C	------------------------------------------------------------
C     PROPG(NGP)    = GEOMETRIC PROPERTIES FOR CURRENT SET
C     EDIS(NEF)     = CURRENT NODAL DISPLACEMENTS
C     COORD(3,NNO)  = CURRENT NODAL COORDINATES (X,Y,Z)
C     COORDI(3,NNO) = ORIGINAL NODAL COORDINATES (X,Y,Z)
C     VSO(3)        = DIRECTION COSINES OF ORIGINAL S-AXIS
C                     (EXCLUDING PRETWIST)
C     SVF           = ASPECT RATIO FACTOR FOR ST. VENANT TORSION
C     NSG           = NUMBER OF SEGMENTS (MAX 20)
C     NNO           = NUMBER OF NODES
C     NLOPT         = NONLINEAR OPTION CODE
C     IWRP          = WARPING OPTION CODE
C
C     V12(3)        = POSITION VECTOR OF NODE 2 RELATIVE TO NODE 1
C     V13(3)        = POSITION VECTOR OF NODE 3 RELATIVE TO NODE 1
C     V14(3)        = POSITION VECTOR OF ATTITUDE NODE RELATIVE TO
C                     NODE 1
C     ------------------------------------------------------------
      DIMENSION EDIS(1),COORD(3,1),COORDI(3,3)
      DIMENSION VSO(3),V12(3)
	DIMENSION VR(3),VT(3)

C
      IF (NLOPT.EQ.3) GOTO 20

      DO 10 I=1,2
      DO 10 J=1,3
   10 COORDI(J,I)=COORD(J,I)

      GOTO 50

   20 CONTINUE

	K=1
      DO 40 I=1,2
      DO 30 J=1,3
      COORDI(J,I)=COORD(J,I)-EDIS(K)
   30 K=K+1
   40 K=K+4

   50 CONTINUE

	DO 60 I=1,3
      V12(I)=COORDI(I,2)-COORDI(I,1)
   60 CONTINUE


	X1 = COORDI(1,1)
	Y1 = COORDI(2,1)
	Z1 = COORDI(3,1)

	X2 = COORDI(1,2)
	Y2 = COORDI(2,2)
	Z2 = COORDI(3,2)

C	NEW BASE VECTOR SONGSAK JUN2006 (SAME AS SPC FRAME)
      CALL FMVECB (X1,Y1,Z1,X2,Y2,Z2,VR,VSO,VT,DUM)

      CALL SCALEN (VSO,VSO,DUM,3)


      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
	SUBROUTINE FRMSTRS(EDATA,FIN,KEG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)

	DIMENSION FIN(7),EDATA(1)
	DIMENSION STCON(4,4),SSCON(4,4)
	DIMENSION ABC(4),PMT(4)

	ALLOCATABLE STC(:,:),STS(:),STN(:)


	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)

	CALL MRELFIL(NAMEI,FNSTS  ,1,2 ,0) !NUMBER OF STRESS POINT 

	CALL MRELFIL(NAMEI,YOUNG  ,1,11,0) !REFERENCE YOUNG MODULUS

	NSTS = INT(FNSTS)

	ALLOCATE(STC(4,NSTS),STS(NSTS),STN(NSTS))
C	--------------------------------------
	DO I = 1,NSTS
	CALL MRELFIL(NAMEI,SN ,1,1+4*(I-1)+120,0) !
	CALL MRELFIL(NAMEI,SS ,1,2+4*(I-1)+120,0) !
	CALL MRELFIL(NAMEI,TT ,1,3+4*(I-1)+120,0) !
	CALL MRELFIL(NAMEI,WP ,1,4+4*(I-1)+120,0) !

	STN(I) = SN

	STC(1,I) = 1.0
	STC(2,I) = SS
	STC(3,I) = TT
	STC(4,I) = WP

	ENDDO
C	--------------------------------------

	CALL MRELFIL(NAMEI,A    ,1,91,0) !
	CALL MRELFIL(NAMEI,QS   ,1,92,0) !
	CALL MRELFIL(NAMEI,QT   ,1,93,0) !
	CALL MRELFIL(NAMEI,SIS  ,1,94,0) !
	CALL MRELFIL(NAMEI,SIT  ,1,95,0) !
	CALL MRELFIL(NAMEI,SIST ,1,96,0) !

	CALL MRELFIL(NAMEI,FWWP ,1,100,0) !
	CALL MRELFIL(NAMEI,FWSP ,1,101,0) !
	CALL MRELFIL(NAMEI,FWTP ,1,102,0) !
	CALL MRELFIL(NAMEI,FWP  ,1,103,0) !


C	A
	STCON(1,1) = A
C	Qs
	STCON(2,1) = QT
C	Qt
	STCON(3,1) = QS
C	Qw
	STCON(4,1) = FWP

C	Qt
	STCON(1,2) = QT
C	It
	STCON(2,2) = SIT
C	Ist
	STCON(3,2) = SIST
C	Isw
	STCON(4,2) = FWSP

C	Qs
	STCON(1,3) = QS
C	Ist
	STCON(2,3) = SIST
C	Is
	STCON(3,3) = SIS
C	Itw
	STCON(4,3) = FWTP

C	Qw
	STCON(1,4) = FWP
C	Isw
	STCON(2,4) = FWSP
C	Itw
	STCON(3,4) = FWTP
C	Iw
	STCON(4,4) = FWWP

C	--------------------------------------
	PMT(1) = FIN(1)
	PMT(2) =-FIN(6)
	PMT(3) = FIN(5)
	PMT(4) = FIN(7)

C	--------------------------------------

	STCON(1:4,1:4) = STCON(1:4,1:4)/YOUNG
	IF(ABS(STCON(4,4)).LT.1.0E-6) THEN
	STCON(4,1:4) = 0.0D0
	STCON(1:4,4) = 0.0D0
	STCON(4,4)   = 1.0D0
	PMT(4) = 0.0D0
	ENDIF
C	--------------------------------------
	N = 4
	CALL INVMATF(STCON,SSCON,N,IB)

C	--------------------------------------
C     =========================================================      
C     MULTIPLYING MATRICES USING 'DGEMM'
C     M, N, KIntegers indicating the size of the matrices:
C     A: M rows by K columns
C     B: K rows by N columns
C     C: M rows by N columns
C     =========================================================   
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
C     SSCON(4,4),PMT(4)
      ALPHA = 1.0D0
      BETA  = 1.0D0
      ABC   = 0.0D0
      CALL DGEMM('N','N',4,1,4,ALPHA,SSCON,4,PMT,4,BETA,ABC,4) 
C	ABC = MATMUL(SSCON,PMT)
C	ABC(1,4),STC(4,NSTS)
      STS   = 0.0D0
	CALL DGEMM('N','N',1,NSTS,4,ALPHA,ABC,1,STC,4,BETA,STS,1)
C	STS = MATMUL(ABC,STC)
	
C	--------------------------------------

	EDATA(1:12) = 0.0D0

	DO I = 1,NSTS
	NS = INT(STN(I))
	IF(I.LE.12.AND.NS.NE.0) EDATA(I) = STS(I)
	ENDDO


C	WRITE(*,'(20E15.4)') STS(1:3) !STC(2,I),ABC(2),STS(I)  !,FIN(6),A,SIT,QT


	DEALLOCATE(STC,STS,STN)

	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE STORFORC(FIN,FCN)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)

	DIMENSION FIN(1)
	DIMENSION FCN(7)

	DO I = 1,7
	FCN(I) = FIN(I)
	ENDDO

	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
	SUBROUTINE CENFORCE(FIN,KEG)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 ITRAN
	
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
      
	DIMENSION FCN(7),FIN(7),ECCMT(3,3)

      
C     -------------------------------------------------	

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)

	CALL MRELFIL(NAMEI,A    ,1,91,0) !
	CALL MRELFIL(NAMEI,QS   ,1,92,0) !
	CALL MRELFIL(NAMEI,QT   ,1,93,0) !
	
	TC = QS/A
	SC = QT/A
	
	ECS = -SC
	ECT = -TC
		    
C     -------------------------------------------------	
      
	ECCMT = 0.0
	ECCMT(1,1) =  0.0
	ECCMT(1,2) = -ECT
	ECCMT(1,3) =  ECS
	ECCMT(2,1) =  ECT
	ECCMT(3,1) = -ECS


	FCN(1:7) = FIN(1:7)

	DO I = 1,3
	DO J = 1,3
	FCN(I+3) = FCN(I+3) + ECCMT(I,J)*FIN(J)
	ENDDO
	ENDDO
	
	FIN(1:7) = FCN(1:7)


	RETURN
	END
C	==============================================================================
C	==============================================================================
C	==============================================================================
      SUBROUTINE FRMASS (S,IMASS,A,AL,VR,VS,VT,STE,AREA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------
C     LUMPED, CONSISTENT DIAGONAL AND CONSISTENT MASS MATRICES
C	--------------------------------------------------------
C     S(NWS)    = STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     DMASS     = MASS-DENSITY 
C     IMASS = 1   LUMPED MASS
C     IMASS = 2   CONSISTENT DIAGONAL MASS
C     IMASS = 3   CONSISTENT MASS
c      COMMON / OPTIONMASS/ IMASS
C     --------------------------------------------------------
      DIMENSION S(*),AM(14,14),AMG(14,14),AMG1(14,14),TT(14,14)
	DIMENSION VR(3),VS(3),VT(3)


	SIR   = STE
	ACON  = A*AL
      !ACON  = 4854.52947501958D0
	BCON  = AL*AL
	AM    = 0.
      IF (IMASS.EQ.0) IMASS = 3
      IF (IMASS.EQ.3)THEN
	AM(1,1)   = ACON/3.
	AM(2,2)   = ACON*13./35.
	AM(3,3)   = ACON*13./35.
	AM(4,4)   = ACON*SIR/3.
	AM(5,5)   = ACON*BCON/105.
	AM(6,6)   = ACON*BCON/105.
	AM(7,7)   = 1.
	AM(8,8)   = ACON/3.
	AM(9,9)   = ACON*13./35.
	AM(10,10) = ACON*13./35.
	AM(11,11) = ACON*SIR/3.
	AM(12,12) = ACON*BCON/105.
	AM(13,13) = ACON*BCON/105.
	AM(14,14) = 1.
      
	AM(1,8)   = ACON/6.
	AM(2,6)   = ACON*AL*11./210.	
	AM(2,9)   = ACON*9./70.
	AM(2,13)  = -ACON*AL*13./420.
	AM(3,5)   = -ACON*AL*11./210.	
	AM(3,10)  = ACON*9./70.
	AM(3,12)  = ACON*AL*13./420.
	AM(4,11)  = ACON*SIR/6.
	AM(5,10)  = -ACON*AL*13./420.
	AM(5,12)  = -ACON*AL*AL/140.
	AM(6,9)   = ACON*AL*13./420.
	AM(6,13)  = -ACON*AL*AL/140.
	AM(9,13)  = -ACON*AL*11./210.
	AM(10,12) = ACON*AL*11./210.

	AM(8,1)   = ACON/6.
	AM(6,2)   = ACON*AL*11./210.	
	AM(9,2)   = ACON*9./70.
	AM(13,2)  = -ACON*AL*13./420.
	AM(5,3)   = -ACON*AL*11./210.	
	AM(10,3)  = ACON*9./70.
	AM(12,3)  = ACON*AL*13./420.
	AM(11,4)  = ACON*SIR/6.
	AM(10,5)  = -ACON*AL*13./420.
	AM(12,5)  = -ACON*AL*AL/140.
	AM(9,6)   = ACON*AL*13./420.
	AM(13,6)  = -ACON*AL*AL/140.
	AM(13,9)  = -ACON*AL*11./210.
	AM(12,10) = ACON*AL*11./210.
      ELSEIF (IMASS.EQ.1)THEN
      AM(1,1)   = ACON*1D0/2D0
	AM(2,2)   = ACON*1D0/2D0
	AM(3,3)   = ACON*1D0/2D0
	AM(4,4)   = ACON*(SIR/AREA)/2D0
	AM(5,5)   = ACON*0D0
	AM(6,6)   = ACON*0D0
	AM(7,7)   = 1.0D0
      AM(8,8)   = ACON*1D0/2D0
	AM(9,9)   = ACON*1D0/2D0
	AM(10,10) = ACON*1D0/2D0
	AM(11,11) = ACON*(SIR/AREA)/2D0
	AM(12,12) = ACON*0D0
	AM(13,13) = ACON*0D0
	AM(14,14) = 1.
      
      ENDIF


C	TRANSFORM LOCAL MASS INTO GLOBAL MASS
	CALL TT1A (VR,VS,VT,TT)
	AMG1 = MATMUL(TT,AM)
	AMG  = MATMUL(AMG1,TRANSPOSE(TT))

C	ARRANGE GLOBAL STIFFNESS MATRIX UPPER TRAGULE ROW WISE
	IT = 0.
	DO I=1,14
	DO J=I,14
	IT = IT + 1
	S(IT) = S(IT)+AMG(I,J)
	END DO 
	END DO
      RETURN

		
      RETURN
      END
C
C	==============================================================================
C	==============================================================================
C	==============================================================================
      SUBROUTINE FRGAUS (NGG,ELN,BPG,BWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SET UP FRAME ELEMENT GAUSS POINT
C	FIRST AND LAST IS FOR BOTTH END WHICH HAVE WEIGH = 0.0
C	=============================================================
	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
C	=============================================================
C	BPG = IS GAUSS POSITION ALONG THE LENGTH OF ELEMENT
C	BWG = IS GAUSS WEIGHT   ALONG THE LENGTH OF ELEMENT
	DIMENSION BPG(1),BWG(1)
	ALLOCATABLE GPL(:),GPW(:)

	ALLOCATE(GPL(NGG),GPW(NGG))

	DO IGG = 1,NGG
	IF(IGG.EQ.1  ) GPL(IGG) = -1.0D0
	IF(IGG.EQ.NGG) GPL(IGG) =  1.0D0 
	IF(IGG.NE.1.AND.IGG.NE.NGG) GPL(IGG) =  GAUSP(IGG-1,NGG-2)
	IF(IGG.EQ.1  ) GPW(IGG) =  0.0D0
	IF(IGG.EQ.NGG) GPW(IGG) =  0.0D0 
	IF(IGG.NE.1.AND.IGG.NE.NGG) GPW(IGG) =  GAUSW(IGG-1,NGG-2)
	ENDDO

	DO IGG = 1,NGG
	RS = GPL(IGG) !GAUSP(IGG,NGG)
	RW = GPW(IGG) !GAUSW(IGG,NGG)
	BPG(IGG) = 0.5*ELN*(1.0 + RS)
	BWG(IGG) = 0.5*ELN*RW
	ENDDO

	DEALLOCATE(GPL,GPW)
		
      RETURN
      END
C
C	================================================================
C
C	END OF FRAME ELEMENT ROUTINE
C
C	================================================================
C	ELASTIC FOUNDATION STIFFNESS BY SONGSAK FEB2006
c	DO IELS = 1,14
c	DO JELS = 1,14
c	SELAS(IELS,JELS) = 0.0
c	ENDDO
c	ENDDO
c
c	STIEL  = 300000.0E3
c	SLAMDA = STIEL/4.0/E/SIZ
c
c	STES1 = STIEL*AL/3.0
c	STES2 = SLAMDA*AL*AL*AL*AL*AL/945.0
c	STES3 = SLAMDA*AL*AL*AL/45.0
c
c	SELAS(2,2) = STES1
c	SELAS(9,9) = STES1
c	SELAS(2,9) = 0.5*STES1
c	SELAS(9,2) = 0.5*STES1
c
c	SELAS(6,6)   = 8.0*STES3
c	SELAS(13,13) = 8.0*STES3
c	SELAS(6,13)  = 31.0*STES3/4.0
c	SELAS(13,6)  = 31.0*STES3/4.0
c
c	SELAS(2,6)   = 4.0*STES2
c	SELAS(2,13)  = 7.0*STES2/2.0
c	SELAS(9,6)   = 7.0*STES2/2.0
c	SELAS(9,13)  = 4.0*STES2
c
c	KPLUS = 0
c	DO IELS = 1,14
c	DO JELS = IELS,14
c	KPLUS = KPLUS+1
c	S(KPLUS) = S(KPLUS) + SELAS(IELS,JELS)
c	ENDDO
c	ENDDO
C	WRITE(*,*) STIEL,AL,E,SIZ







