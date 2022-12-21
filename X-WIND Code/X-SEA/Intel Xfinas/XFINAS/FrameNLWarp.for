C	-------------------------------------------------------------------
      SUBROUTINE BFRAMEW (PROPM,PROPG,S,COORD,EDIS,EDISI,RE,PROPO,
	1					FIN,ISET,WA,FIXEN,FIXLR,FIXEO,FIXLO,LREAS,
     2					PMATRL,FIXEN_OFF,FIXLR_OFF,FIXEO_OFF,FIXLO_OFF)
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

C     SAVE ELEMENT NUMBER      
      COMMON /NELEM/ IEL  

C	-------------------------------------------------------------------    
C	-------------------------------------------------------------------
      DIMENSION PROPM(1),PROPG(1),S(105),COORD(3,2),EDIS(14)
      DIMENSION EDISI(14),DR(49)
      DIMENSION DISD(9),REDIS(14),COORDI(3,2)
      DIMENSION VRO(3),VSO(3),VTO(3),VR(3),VS(3),VT(3),EC(3)
     	DIMENSION R(1,2)
	DIMENSION FIN(NEF)
	DIMENSION AKG(14,14),AKL(14,14),BKG1(14,14),BKG(14,14)
	DIMENSION BMGG(14),BWG(10),BPG(10),BMSG(7,14),DRMAT(7,7)
	DIMENSION BBX(14),TT(14,14),EDISL(14)
	DIMENSION WA(1),LSEH(14),FIXEN(NEF),FIXLR(NEF)
      DIMENSION FIXEN_OFF(NEF),FIXLR_OFF(NEF),FIXEO_OFF(NEF),FIXLO_OFF(NEF)
	DIMENSION FIXEO(NEF),FIXLO(NEF)
	DIMENSION TRANH(14,14),LREAS(14),TRANO(14,14),PROPO(6)
	DIMENSION COORO(3,2),TRANOL(14,14)

	DIMENSION DSIGM(5,5),QMAT(5,5),FMAT(5,5),TRANM(5,5),TRNFLX(14,5)
	DIMENSION IJ(4)

	DIMENSION PMATRL(NMP,1)
	DIMENSION IGSET(NELE),NELEMENT(NELE)


	DIMENSION EPS4(4)
	DIMENSION BMATX(7,14),DM(7,7),DG(7,7),EPS(7),SIGR(7),RE(14)
      DIMENSION SVAL(30)
	
C     FOR DGMM FUNCTION
	DIMENSION COUT1(5,5),COUT2(5,5),COUT3(14),COUT4(7,14),COUT5(5,14),COUT6(14,14),COUT7(14)
	
	DIMENSION EDISLT(14)
      
      DIMENSION WFLOCAL(14)
      

	ALLOCATABLE FGVEC(:,:),GEOPR(:,:),GEOVC(:,:)



C	CALLING ELEMENT GAUSS POINT NUMBER
	CALL INTFILL('OGRP',NGG ,5 ,KEG,0) !


	ALLOCATE(FGVEC(11,NGG),GEOPR(8,NGG),GEOVC(5,NGG))

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
	
C	TT(14,14),REDIS(14)	
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
      
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
	BXD = BPG(II)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
	CALL BBXFRMW(BMATX,ELN,BXD) 

C     ------------------------------------------------------
C     ELASTIC RIGIDITY MATRIX (DR)
C     ------------------------------------------------------
      IF(MTMOD.GT.2) GOTO 220
	
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFILA(NAMEI,SVAL,1,91,103,0) !
      EA   = SVAL(1)
      EQS  = SVAL(2)
      EQt  = SVAL(3)
      EIs  = SVAL(4)
      EIt  = SVAL(5)
      EIst = SVAL(6)
      GJr  = SVAL(7)
      ASS  = SVAL(8)
      ATT  = SVAL(9)
      EIw  = SVAL(10)
      EIsw = SVAL(11)
      EItw = SVAL(12)
      EQw  = SVAL(13)

CC	CALL MRELFIL(NAMEI,EA   ,1,91 ,0) !
CC	CALL MRELFIL(NAMEI,EQS  ,1,92 ,0) !
CC	CALL MRELFIL(NAMEI,EQt  ,1,93 ,0) !
CC	CALL MRELFIL(NAMEI,EIs  ,1,94 ,0) !
CC	CALL MRELFIL(NAMEI,EIt  ,1,95 ,0) !
CC	CALL MRELFIL(NAMEI,EIst ,1,96 ,0) !
      !EIst = -2.142786979675293D-005
CC	CALL MRELFIL(NAMEI,GJr  ,1,97 ,0) !
CC	CALL MRELFIL(NAMEI,EIw  ,1,100,0) !
CC	CALL MRELFIL(NAMEI,EIsw ,1,101,0) !
CC	CALL MRELFIL(NAMEI,EItw ,1,102,0) !
CC	CALL MRELFIL(NAMEI,EQw  ,1,103,0) !

CC	CALL MRELFIL(NAMEI,ASS  ,1,98 ,0) !
CC	CALL MRELFIL(NAMEI,ATT  ,1,99 ,0) !


	CALL MRELFIL(NAMEI,ASS1 ,1,28 ,0) !
	CALL MRELFIL(NAMEI,ATT1 ,1,29 ,0) !


	DM(1:7,1:7) = 0.0D0

	DM(1,1) = EA
	DM(4,4) = GJr

	DM(1,5) = EQs
	DM(5,5) = EIs
	
	DM(1,6) =-EQt
	DM(5,6) =-EIst
	DM(6,6) = EIt

	DM(1,7) = EQw
	DM(5,7) = EItw
	DM(6,7) =-EIsw
	DM(7,7) = EIw


	DO I = 1,7
	DO J = I,7
	DM(J,I) = DM(I,J)
	ENDDO
	ENDDO


220	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFILA(NAMEI,SVAL,1,21,34,0) !
      GAr = SVAL(1)
      GIs = SVAL(4)
      GIt = SVAL(5)
      GIst= SVAL(6)
      GIw = SVAL(10)
      GIsw= SVAL(11)
      GItw= SVAL(12)
      GItp= SVAL(14)
CC	CALL MRELFIL(NAMEI,GAr  ,1,21,0) !
CC	CALL MRELFIL(NAMEI,GIs  ,1,24,0) !
CC	CALL MRELFIL(NAMEI,GIt  ,1,25,0) !
CC	CALL MRELFIL(NAMEI,GIst ,1,26,0) !
CC	CALL MRELFIL(NAMEI,GIw  ,1,30,0) !
CC	CALL MRELFIL(NAMEI,GIsw ,1,31,0) !
CC	CALL MRELFIL(NAMEI,GItw ,1,32,0) !
CC	CALL MRELFIL(NAMEI,GIp  ,1,34,0) !

	GEOPR(1,II) = GIs
	GEOPR(2,II) = GIt
	GEOPR(3,II) = GIst
	GEOPR(4,II) = GIsw
	GEOPR(5,II) = GItw
	GEOPR(6,II) = GIw
	GEOPR(7,II) = GIp
	GEOPR(8,II) = GAr

	CALL MRELFIL(NAMEI,AREAD  ,1,61,0) !
      CALL RELFILL('TWGH',AREAD*BWG(II),1,KEG,2)
	
	IF (ITASK.NE.5) GOTO 245
C     -------------------------------------------
C     CONSISTENT AND LUMPED MASS MATRIX (ITASK=5)
C     -------------------------------------------

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,AREAD  ,1,61,0) !
	CALL MRELFIL(NAMEI,PJLD   ,1,74,0) !
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
C     STRAIN TERMS 
C     ------------------------------------------------------------
	EPS(1:7) = 0.0D0
	DO IEF = 1,14
	EPS(1) = EPS(1) + BMATX(1,IEF)*EDISL(IEF)  !dU /dX
	EPS(2) = EPS(2) + BMATX(2,IEF)*EDISL(IEF)  !dV /dX
	EPS(3) = EPS(3) + BMATX(3,IEF)*EDISL(IEF)  !dW /dX
	EPS(4) = EPS(4) + BMATX(4,IEF)*EDISL(IEF)  !dPr/dX
	EPS(5) = EPS(5) + BMATX(5,IEF)*EDISL(IEF)  !dPs/dX
	EPS(6) = EPS(6) + BMATX(6,IEF)*EDISL(IEF)  !dPt/dX
	EPS(7) = EPS(7) + BMATX(7,IEF)*EDISL(IEF)  !ddPr/ddX
	ENDDO

C     ---------------------------------
C     RESULTANTS SIGR(8)
C     ---------------------------------

C	---------------------------
C     LINEAR ELASTIC
	IF (MTMOD.EQ.1) THEN

	SIGR(1:7) = 0.0D0
	DO ISR = 1,7
	SIGR(1) = SIGR(1) + DM(1,ISR)*EPS(ISR)  !P
	SIGR(4) = SIGR(4) + DM(4,ISR)*EPS(ISR)  !Mr
	SIGR(5) = SIGR(5) + DM(5,ISR)*EPS(ISR)  !Ms
	SIGR(6) = SIGR(6) + DM(6,ISR)*EPS(ISR)  !Mt
	SIGR(7) = SIGR(7) + DM(7,ISR)*EPS(ISR)  !Mw
	ENDDO

	EPS4(1:4) = [EPS(1),EPS(5),EPS(6),EPS(7)]
	CALL FORCGEO(FGVEC(1,II),EPS4,KEG)

	GEOVC(1:5,II) = [SIGR(1),SIGR(4),SIGR(5),-SIGR(6),SIGR(7)]

	ENDIF
C	---------------------------

C     ELASTO PLASTIC
	IF (MTMOD.EQ.3) THEN
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FNFIB,1,4,0) !
	NFIB = INT(FNFIB)
	CALL FMPLAS(PMATRL,EPS,DM,SIGR,WA,II,NFIB,NMP,ISET,KEG) 

	DO I = 1,7
	DO J = I,7
	DM(J,I) = DM(I,J)
	ENDDO
	ENDDO

	FGVEC(1:11,II) = 0.0D0
	FGVEC(1,II) = SIGR(1)
	FGVEC(2,II) = SIGR(5)
	FGVEC(3,II) = SIGR(6)
	FGVEC(7,II) = SIGR(7)

	
	GEOVC(1:5,II) = [SIGR(1),SIGR(4),SIGR(5),-SIGR(6),SIGR(7)]

	ENDIF
C	---------------------------
C     CONCRETE
	IF (MTMOD.EQ.5) THEN
	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFIL(NAMEI,FNFIB,1,4,0) !
	NFIB = INT(FNFIB)
	CALL SPCONT(EPS,II,DM,SIGR,NMP,PMATRL,
	1			WA(1),WA(8*NGG+1),NLOPT,NFIB,ISET,KEG,NGG)

	DO I = 1,7
	DO J = I,7
	DM(J,I) = DM(I,J)
	ENDDO
	ENDDO

	FGVEC(1:11,II) = 0.0D0
	FGVEC(1,II) = SIGR(1)
	FGVEC(2,II) = SIGR(5)
	FGVEC(3,II) = SIGR(6)
	FGVEC(7,II) = SIGR(7)

	
	GEOVC(1:5,II) = [SIGR(1),SIGR(4),SIGR(5),-SIGR(6),SIGR(7)]


	ENDIF
C	---------------------------


C	D11 = EA 
C	D12 = EQt
C	D13 = EQs
C	D21 = EQt 
C	D22 = EIt
C	D23 = EIst
C	D31 = EQs
C	D32 = EIst
C	D33 = EIs


	D11 = DM(1,1) 
	D12 = DM(1,6)  
	D13 = DM(1,5)
	D21 = D12 
	D22 = DM(6,6)
	D23 =-DM(5,6)
	D31 = D13
	D32 = D23
	D33 = DM(5,5)

	DSIGM = 0.0D0

	DETS = D11*D22*D33+D21*D13*D32+D31*D12*D23
	1	  -D11*D23*D32-D21*D12*D33-D31*D13*D22
	IF(DETS.NE.0.0D0) THEN
	    DSIGM(1,1) = (D22*D33-D23*D32)/DETS
	    DSIGM(1,2) = (D13*D32-D12*D33)/DETS
	    DSIGM(1,3) = (D12*D23-D13*D22)/DETS

	    DSIGM(2,1) = (D23*D31-D21*D33)/DETS
	    DSIGM(2,2) = (D11*D33-D13*D31)/DETS
	    DSIGM(2,3) = (D13*D21-D11*D23)/DETS

	    DSIGM(3,1) = (D21*D32-D22*D31)/DETS
	    DSIGM(3,2) = (D12*D31-D11*D32)/DETS
	    DSIGM(3,3) = (D11*D22-D12*D21)/DETS
      ENDIF
      
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
	RE = RE + BWG(II)*MATMUL(TRANSPOSE(BMATX),SIGR)

C     BMATX(7,14),DM(7,7)
	AKL = MATMUL(TRANSPOSE(BMATX),MATMUL(DM,BMATX))	

	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = AKG(I,J) + BWG(II)*AKL(I,J)
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
      
      chana = 3

305	CONTINUE

	ENDIF

C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
	CALL TRNMUL(TRANH,RE,2)
	CALL TRNMUM(TRANH,AKG)
C	------------------------------------------------------------

	IF (IFEIG.EQ.0.AND.ISOLOP.EQ.4) GOTO 800

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
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
	RE = MATMUL(TT,RE)
C	----------------------------------------------------------
C	TRANSFORMATION DUE TO OFFSET
C	----------------------------------------------------------
	CALL TRNMUL(TRANO,RE,2)
C	----------------------------------------------------------


      CALL TRAPIZOIDAL_TRANSFORMATION(VR,VS,VT,'WRT',WFLOCAL,ILC)
      
	DO IEF = 1,NEF
	FIXLR(IEF) = 0.0
	FIXLO(IEF) = 0.0
	DO JEF = 1,NEF
	FIXLR(IEF) = FIXLR(IEF) + TT(IEF,JEF)*FIXEN(JEF)  !VARY FIXEND
	FIXLO(IEF) = FIXLO(IEF) + TT(IEF,JEF)*FIXEO(JEF)  !CONT FIXEND

	ENDDO
      ENDDO
      
      DO IEF = 1,NEF
	FIXLR_OFF(IEF) = 0.0
	FIXLO_OFF(IEF) = 0.0
	DO JEF = 1,NEF
	FIXLR_OFF(IEF) = FIXLR_OFF(IEF) + TT(IEF,JEF)*FIXEN_OFF(JEF)  !VARY FIXEND (OFFSHORE LOAD CASE)
	FIXLO_OFF(IEF) = FIXLO_OFF(IEF) + TT(IEF,JEF)*FIXEO_OFF(JEF)  !CONT FIXEND (OFFSHORE LOAD CASE)
	ENDDO
      ENDDO
C	CALL TRNMUL(TRANOL,FIXLR,1)
      
      chana =3

	CALL TT1A (VR,VS,VT,TT)
	
      
C	TT(14,14),AKG(14,14)
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

800	AKG(1:14,1:14) = 0.0


	DO II = 1,NGG

	BXD = BPG(II)

C     OBTAIN LINEAR B MATRIX AT THE REFERENCE AXIS AND STRAIN
	CALL BBXFRMW(BMATX,ELN,BXD)

C	CALL DGEOMAT_OLD(FGVEC(1,II),DG)
	CALL DGEOMAT(ELN,BXD,GEOPR(1,II),GEOVC(1,II),AKL)

C	AKL = MATMUL(TRANSPOSE(BMATX),MATMUL(DG,BMATX))

	DO I = 1,14
	DO J = 1,14
	AKG(I,J) = AKG(I,J) + BWG(II)*AKL(I,J)
	ENDDO
	ENDDO

	ENDDO !END II


C	------------------------------------------------------------
C	TRANSFORM CORRESPONDING RELEASE CONDITION
C	------------------------------------------------------------
	CALL TRNMUM(TRANH,AKG)
C	------------------------------------------------------------


	CALL TT1A (VR,VS,VT,TT)

C	TT(14,14),AKG(14,14)
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

820   CONTINUE 
C     --------------------------------------------------------


900   CONTINUE

	DEALLOCATE(FGVEC,GEOPR,GEOVC)



	RETURN

      END

C     =================================================================
C     =================================================================
C     =================================================================
	SUBROUTINE BBXFRMW(BBX,ELN,BXD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /NUMB2/ IDOFB(9)

	DIMENSION BBX(7,14)

C	CHECK WHETHER WARPING IS INCLUDE OR NOT
	LWARP = 0
	IF(IDOFB(7).EQ.0) LWARP = 1 !INCLUDE WARPING

	dUr1  = -1.0/ELN
	dUr8  =  1.0/ELN

	dVs2  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
	dVs6  = +1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
	dVs9  = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
	dVs13 = +1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN

	 
	dWt3  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
	dWt5  = -1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
	dWt10 = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
	dWt12 = -1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN


	dPr4  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
	dPr7  = +1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
	dPr11 = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
	dPr14 = +1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN


	dPs3  = -6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0      
	dPs5  = +2.0*(-2.0+3.0*BXD/ELN)/ELN
	dPs10 = -6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0  
	dPs12 = +2.0*(-1.0+3.0*BXD/ELN)/ELN	 
 
	dPt2  = +6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0   
	dPt6  = +2.0*(-2.0+3.0*BXD/ELN)/ELN
	dPt9  = +6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0 
	dPt13 = +2.0*(-1.0+3.0*BXD/ELN)/ELN
	
	ddPr4 = +6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0   
	ddPr7 = +2.0*(-2.0+3.0*BXD/ELN)/ELN
	ddPr11= +6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0 
	ddPr14= +2.0*(-1.0+3.0*BXD/ELN)/ELN
	
	
	IF(LWARP.EQ.0) THEN !IF WARPING DO NOT INCLUDE
	    dPr4  = -1.0/ELN
	    dPr7  = 0.0D0
	    dPr11 = +1.0/ELN
	    dPr14 = 0.0D0
    	
	    ddPr4 = 0.0D0   
	    ddPr7 = 0.0D0  
	    ddPr11= 0.0D0   
	    ddPr14= 0.0D0  
	ENDIF
	
	
	BBX(1:7,1:14) = 0.0D0
	
C	dU/dX
	BBX(1,1 ) = dUr1
	BBX(1,8 ) = dUr8

C	dV/dX
	BBX(2,2 ) = dVs2
	BBX(2,6 ) = dVs6
	BBX(2,9 ) = dVs9
	BBX(2,13) = dVs13
	
C	dW/dX
	BBX(3,3 ) = dWt3
	BBX(3,5 ) = dWt5
	BBX(3,10) = dWt10
	BBX(3,12) = dWt12
	

C	dPr/dX
	BBX(4,4 ) = dPr4
	BBX(4,7 ) = dPr7
	BBX(4,11) = dPr11
	BBX(4,14) = dPr14

C	dPs/dX
	BBX(5,3 ) = dPs3
	BBX(5,5 ) = dPs5
	BBX(5,10) = dPs10
	BBX(5,12) = dPs12
	
C	dPt/dX
	BBX(6,2 ) = dPt2
	BBX(6,6 ) = dPt6
	BBX(6,9 ) = dPt9
	BBX(6,13) = dPt13


C	ddPr/ddX
	BBX(7,4 ) = ddPr4
	BBX(7,7 ) = ddPr7
	BBX(7,11) = ddPr11
	BBX(7,14) = ddPr14



	RETURN

	END


C	==========================================================================
C	==========================================================================
C	==========================================================================
	SUBROUTINE FORCGEO(FGVEC,EPS4,KEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	-------------------------------------------------------------------
      CHARACTER*1 NAMEI(4)
      DIMENSION   INAME(4)
C	-------------------------------------------------------------------

	DIMENSION FGVEC(11),FGMAT(11,4),EPS4(4),SVAL(20)

	INAME(1:4) = [5,0,1,KEG] !XSEC
	CALL ICONC(INAME,NAMEI)
	CALL MRELFILA(NAMEI,SVAL,1,91,103,0) !
      EA   = SVAL(1)
      EQS  = SVAL(2)
      EQt  = SVAL(3)
      EIs  = SVAL(4)
      EIt  = SVAL(5)
      EIst = SVAL(6)
      GJr  = SVAL(7)
      EIw  = SVAL(10)
      EIsw = SVAL(11)
      EItw = SVAL(12)
      EQw  = SVAL(13)   
CC	CALL MRELFIL(NAMEI,EA   ,1,91,0) !
CC	CALL MRELFIL(NAMEI,EQs  ,1,92,0) !
CC	CALL MRELFIL(NAMEI,EQt  ,1,93,0) !
CC	CALL MRELFIL(NAMEI,EIs  ,1,94,0) !
CC	CALL MRELFIL(NAMEI,EIt  ,1,95,0) !
CC	CALL MRELFIL(NAMEI,EIst ,1,96,0) !
CC	CALL MRELFIL(NAMEI,GJr  ,1,97,0) !
CC	CALL MRELFIL(NAMEI,EIw  ,1,100,0) !
CC	CALL MRELFIL(NAMEI,EIsw ,1,101,0) !
CC	CALL MRELFIL(NAMEI,EItw ,1,102,0) !
CC	CALL MRELFIL(NAMEI,EQw  ,1,103,0) !

	EQw  = 0.0D0
	EIsw = 0.0D0
	EItw = 0.0D0
	EIw  = 0.0D0

	EIttt= 0.0D0
	EIsss= 0.0D0
	EIwww= 0.0D0
	EIstt= 0.0D0
	EIsst= 0.0D0
	EIttw= 0.0D0
	EIssw= 0.0D0
	EIsww= 0.0D0
	EItww= 0.0D0
	EIstw= 0.0D0

	FGMAT(1:11,1:4) = 0.0D0

	FGMAT(1 ,1:4) = [EA   ,EQs  ,-EQt  ,EQw  ]	!P
	FGMAT(2 ,1:4) = [EQs  ,EIs  ,-EIst ,EItw ]	!Ms
	FGMAT(3 ,1:4) = [EQt  ,EIst ,-EIt  ,EIsw ]	!Mt
	FGMAT(4 ,1:4) = [EIs  ,EIttt,-EIstt,EIttw]	!Fs
	FGMAT(5 ,1:4) = [EIt  ,EIsst,-EIsss,EIssw]	!Ft
	FGMAT(6 ,1:4) = [EIst ,EIstt,-EIsst,EIstw]	!Fst
	FGMAT(7 ,1:4) = [EQw  ,EItw ,-EIsw ,EIw  ]	!Mw
	FGMAT(8 ,1:4) = [EItw ,EIttw,-EIstw,EItww]	!Fsw
	FGMAT(9 ,1:4) = [EIsw ,EIstw,-EIssw,EIsww]	!Ftw
	FGMAT(10,1:4) = [EIw  ,EItww,-EIsww,EIwww]	!Fw
	FGMAT(11,1:4) = [EIt+EIs,EIsst+EIttt,-EIsss-EIstt,EIssw+EIttw]	!Fp

C     FGMAT(11,4),EPS4(4)
	FGVEC = MATMUL(FGMAT,EPS4)


	RETURN

	END


C	==========================================================================
C	==========================================================================
C	==========================================================================
	SUBROUTINE DGEOMAT_OLD(FGVEC,DG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION FGVEC(11),DG(7,7)

	DG(1:7,1:7) = 0.0D0

	DG(1,1) = FGVEC(1)		! P
	DG(2,2) = FGVEC(1)		! P
	DG(3,3) = FGVEC(1)		! P 

	DG(2,4) =-FGVEC(2)		!-Ms
	DG(3,4) = FGVEC(3)		! Mt
	DG(4,4) = FGVEC(11)		! Fp

	DG(1,5) = FGVEC(2)		! Ms
	DG(5,5) = FGVEC(4)		! Fs

	DG(1,6) =-FGVEC(3)		!-Mt
	DG(5,6) =-FGVEC(6)		!-Fst
	DG(6,6) = FGVEC(5)		! Ft


	DG(1,7) = FGVEC(7)		! Mw
	DG(5,7) = FGVEC(8)		! Fsw
	DG(6,7) =-FGVEC(9)		!-Ftw
	DG(7,7) = FGVEC(10)		! Fw

	DO I = 1,7
	DO J = I,7
	DG(J,I) = DG(I,J)
	ENDDO
	ENDDO


	RETURN

	END


C	==========================================================================
C	==========================================================================
C	==========================================================================
c	SUBROUTINE DGEOMAT_OLD(ELN,BXD,GEOPR,GEOVC,AKL)
c	IMPLICIT REAL*8 (A-H,O-Z)
c      IMPLICIT INTEGER*4 (I-N)
c
c      COMMON /NUMB2/ IDOFB(9)
c
c	DIMENSION BBX(8,14),H(4),GEOPR(1),GEOVC(1),DG(8,8),AKL(14,14)
c	
cC	CHECK WHETHER WARPING IS INCLUDE OR NOT
c	LWARP = 0
c	IF(IDOFB(7).EQ.0) LWARP = 1 !INCLUDE WARPING
c	
c
c	dUr1  = -1.0/ELN
c	dUr8  =  1.0/ELN
c
c	dVs2  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
c	dVs6  = +1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
c	dVs9  = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
c	dVs13 = +1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN
c
c	 
c	dWt3  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
c	dWt5  = -1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
c	dWt10 = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
c	dWt12 = -1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN
c
c
c	dPr4  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
c	dPr7  = +1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
c	dPr11 = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
c	dPr14 = +1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN
c
c
c	dPs3  = -6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0      
c	dPs5  = +2.0*(-2.0+3.0*BXD/ELN)/ELN
c	dPs10 = -6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0  
c	dPs12 = +2.0*(-1.0+3.0*BXD/ELN)/ELN	 
c 
c	dPt2  = +6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0   
c	dPt6  = +2.0*(-2.0+3.0*BXD/ELN)/ELN
c	dPt9  = +6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0 
c	dPt13 = +2.0*(-1.0+3.0*BXD/ELN)/ELN
c	
c	ddPr4 = +6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0   
c	ddPr7 = +2.0*(-2.0+3.0*BXD/ELN)/ELN
c	ddPr11= +6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0 
c	ddPr14= +2.0*(-1.0+3.0*BXD/ELN)/ELN
c
c
c	H(1) = 1.0 - 3.0*BXD*BXD/ELN**2.0 + 2.0*BXD*BXD*BXD/ELN**3.0
c	H(2) = BXD - 2.0*BXD*BXD/ELN      +     BXD*BXD*BXD/ELN**2.0
c	H(3) =       3.0*BXD*BXD/ELN**2.0 - 2.0*BXD*BXD*BXD/ELN**3.0
c	H(4) =          -BXD*BXD/ELN      +     BXD*BXD*BXD/ELN**2.0
c
c	IF(LWARP.EQ.0) THEN !IF WARPING DO NOT INCLUDE
c	    dPr4  = -1.0/ELN
c	    dPr7  = 0.0D0
c	    dPr11 = +1.0/ELN
c	    dPr14 = 0.0D0
c    	
c	    ddPr4 = 0.0D0   
c	    ddPr7 = 0.0D0  
c	    ddPr11= 0.0D0   
c	    ddPr14= 0.0D0  
c	    
c	    H(1) = 1.0 - BXD/ELN
c	    H(2) = 0.0D0  
c	    H(3) = BXD/ELN
c	    H(4) = 0.0D0  
c	ENDIF
c
c	BBX(1:8,1:14) = 0.0D0
c	
cC	dU/dX
c	BBX(1,1 ) = dUr1
c	BBX(1,8 ) = dUr8
c
cC	dV/dX
c	BBX(2,2 ) = dVs2
c	BBX(2,6 ) = dVs6
c	BBX(2,9 ) = dVs9
c	BBX(2,13) = dVs13
c	
cC	dW/dX
c	BBX(3,3 ) = dWt3
c	BBX(3,5 ) = dWt5
c	BBX(3,10) = dWt10
c	BBX(3,12) = dWt12
c
cC	Pr
c	BBX(4,4 ) = H(1)
c	BBX(4,7 ) = H(2)
c	BBX(4,11) = H(3)
c	BBX(4,14) = H(4)
c	
c
cC	dPr/dX
c	BBX(5,4 ) = dPr4
c	BBX(5,7 ) = dPr7
c	BBX(5,11) = dPr11
c	BBX(5,14) = dPr14
c
cC	dPs/dX
c	BBX(6,3 ) = dPs3
c	BBX(6,5 ) = dPs5
c	BBX(6,10) = dPs10
c	BBX(6,12) = dPs12
c	
cC	dPt/dX
c	BBX(7,2 ) = dPt2
c	BBX(7,6 ) = dPt6
c	BBX(7,9 ) = dPt9
c	BBX(7,13) = dPt13
c
c
cC	ddPr/ddX
c	BBX(8,4 ) = ddPr4
c	BBX(8,7 ) = ddPr7
c	BBX(8,11) = ddPr11
c	BBX(8,14) = ddPr14
c
c	
c	DG(1:8,1:8) = 0.0D0
c	ALPHA = 0.5
c
c
c	GIs = GEOPR(1) 
c	GIt = GEOPR(2) 
c	GIst= GEOPR(3) 
c	GIsw= GEOPR(4) 
c	GItw= GEOPR(5) 
c	GIw = GEOPR(6)
c	GIp = GEOPR(7) 
c	GAr = GEOPR(8) 
c
c
c	Fr = GEOVC(1) 
c	FMr= GEOVC(2) 
c	FMs= GEOVC(3) 
c	FMt= GEOVC(4) 
c	FMw= GEOVC(5) 
c
cC	-------------------------------------
c	DG(1,1) = Fr
c	DG(6,1) =-FMs
c	DG(7,1) =-FMt
c	DG(8,1) = FMw
c
c	DG(2,2) = Fr
cC	DG(5,2) = FMs/2.0D0
c	DG(6,2) = FMr*(1.0D0-ALPHA)/2.0D0
c
c	DG(3,3) = Fr
cC	DG(5,3) = FMt/2.0D0
c	DG(7,3) = FMr*ALPHA/2.0D0
c
c	DG(6,4) = FMt !/2.0D0
c	DG(7,4) =-FMs !/2.0D0
c
c	DG(5,5) = Fr*GIp/GAr
c
c	DG(6,6) = Fr*GIs/GAr
c	DG(7,6) =-Fr*GIst/GAr
c	DG(8,6) = Fr*GItw/GAr
c
c	DG(7,7) = Fr*GIt/GAr
c	DG(8,7) =-Fr*GIsw/GAr
c
c	DG(8,8) = Fr*GIw/GAr
cC	-------------------------------------
c
c
c	DO I = 1,8
c	DO J = I,8
c	DG(I,J) = DG(J,I)
c	ENDDO
c	ENDDO
c
c
cC	DO I = 1,8
cC	WRITE(110,111) DG(I,1:8)
cC	ENDDO
cC	PAUSE
c111	FORMAT(8E13.4)
c
c	AKL = MATMUL(TRANSPOSE(BBX),MATMUL(DG,BBX))
c
c
c
c	RETURN
c
c	END
c
c
cC	==========================================================================
C	==========================================================================
C	==========================================================================

	SUBROUTINE DGEOMAT(ELN,BXD,GEOPR,GEOVC,AKL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      
      COMMON /NUMB2/ IDOFB(9)

	DIMENSION BBX(8,14),H(4),GEOPR(1),GEOVC(1),DG(8,8),AKL(14,14)
	
C	CHECK WHETHER WARPING IS INCLUDE OR NOT
	LWARP = 0
	IF(IDOFB(7).EQ.0) LWARP = 1 !INCLUDE WARPING

	dUr1  = -1.0/ELN
	dUr8  =  1.0/ELN

	dVs2  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
	dVs6  = +1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
	dVs9  = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
	dVs13 = +1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN

	 
	dWt3  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
	dWt5  = -1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
	dWt10 = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
	dWt12 = -1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN
	
c	dWt3  = -1.0D0*dWt3 
c	dWt5  = -1.0D0*dWt5
c	dWt10 = -1.0D0*dWt10
c	dWt12 = -1.0D0*dWt12 

	dPr4  = +6.0*(       -BXD+    BXD*BXD/ELN)/ELN**2.0
	dPr7  = +1.0*(ELN-4.0*BXD+3.0*BXD*BXD/ELN)/ELN
	dPr11 = +6.0*(        BXD-    BXD*BXD/ELN)/ELN**2.0
	dPr14 = +1.0*(   -2.0*BXD+3.0*BXD*BXD/ELN)/ELN


	dPs3  = -6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0      
	dPs5  = +2.0*(-2.0+3.0*BXD/ELN)/ELN
	dPs10 = -6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0  
	dPs12 = +2.0*(-1.0+3.0*BXD/ELN)/ELN	 
	
	dPs3  = -1.0D0*dPs3 
	dPs5  = -1.0D0*dPs5
	dPs10 = -1.0D0*dPs10
	dPs12 = -1.0D0*dPs12 
 
	dPt2  = +6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0   
	dPt6  = +2.0*(-2.0+3.0*BXD/ELN)/ELN
	dPt9  = +6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0 
	dPt13 = +2.0*(-1.0+3.0*BXD/ELN)/ELN
	
	ddPr4 = +6.0*(-1.0+2.0*BXD/ELN)/ELN**2.0   
	ddPr7 = +2.0*(-2.0+3.0*BXD/ELN)/ELN
	ddPr11= +6.0*( 1.0-2.0*BXD/ELN)/ELN**2.0 
	ddPr14= +2.0*(-1.0+3.0*BXD/ELN)/ELN


	H(1) = 1.0 - 3.0*BXD*BXD/ELN**2.0 + 2.0*BXD*BXD*BXD/ELN**3.0
	H(2) = BXD - 2.0*BXD*BXD/ELN      +     BXD*BXD*BXD/ELN**2.0
	H(3) =       3.0*BXD*BXD/ELN**2.0 - 2.0*BXD*BXD*BXD/ELN**3.0
	H(4) =          -BXD*BXD/ELN      +     BXD*BXD*BXD/ELN**2.0

	IF(LWARP.EQ.0) THEN !IF WARPING DO NOT INCLUDE
	    dPr4  = -1.0/ELN
	    dPr7  = 0.0D0
	    dPr11 = +1.0/ELN
	    dPr14 = 0.0D0
    	
	    ddPr4 = 0.0D0   
	    ddPr7 = 0.0D0  
	    ddPr11= 0.0D0   
	    ddPr14= 0.0D0  
	    
	    H(1) = 1.0 - BXD/ELN
	    H(2) = 0.0D0  
	    H(3) = BXD/ELN
	    H(4) = 0.0D0  
	ENDIF

	BBX(1:8,1:14) = 0.0D0
	
C	dU/dX
	BBX(1,1 ) = dUr1
	BBX(1,8 ) = dUr8

C	dV/dX
	BBX(2,2 ) = dVs2
	BBX(2,6 ) = dVs6
	BBX(2,9 ) = dVs9
	BBX(2,13) = dVs13
	
C	dW/dX
	BBX(3,3 ) = dWt3
	BBX(3,5 ) = dWt5
	BBX(3,10) = dWt10
	BBX(3,12) = dWt12

C	Pr
	BBX(4,4 ) = H(1)
	BBX(4,7 ) = H(2)
	BBX(4,11) = H(3)
	BBX(4,14) = H(4)
	

C	dPr/dX
	BBX(5,4 ) = dPr4
	BBX(5,7 ) = dPr7
	BBX(5,11) = dPr11
	BBX(5,14) = dPr14

C	dPs/dX
	BBX(6,3 ) = dPs3
	BBX(6,5 ) = dPs5
	BBX(6,10) = dPs10
	BBX(6,12) = dPs12
	
C	dPt/dX
	BBX(7,2 ) = dPt2
	BBX(7,6 ) = dPt6
	BBX(7,9 ) = dPt9
	BBX(7,13) = dPt13


C	ddPr/ddX
	BBX(8,4 ) = ddPr4
	BBX(8,7 ) = ddPr7
	BBX(8,11) = ddPr11
	BBX(8,14) = ddPr14

	
	DG(1:8,1:8) = 0.0D0
	ALPHA = 0.5


	GIs = GEOPR(1) 
	GIt = GEOPR(2) 
	GIst= GEOPR(3) 
	GIsw= GEOPR(4) 
	GItw= GEOPR(5) 
	GIw = GEOPR(6)
	GIp = GEOPR(7) 
	GAr = GEOPR(8) 


	Fr = GEOVC(1) 
	FMr= GEOVC(2) 
	FMs= GEOVC(3) 
	FMt= GEOVC(4) 
	FMw= GEOVC(5) 

C	-------------------------------------
C	DG(1,1) = Fr
	DG(6,1) = FMs
	DG(7,1) = FMt
	DG(8,1) = FMw

	DG(2,2) = Fr
 	DG(5,2) =-FMs !/2.0D0
	DG(6,2) = FMr*(1.0D0-ALPHA)/2.0D0

	DG(3,3) = Fr
	DG(5,3) = FMt !/2.0D0
	DG(7,3) = FMr*ALPHA/2.0D0

C	DG(6,4) =-FMt/2.0D0
C	DG(7,4) = FMs/2.0D0

	DG(5,5) = Fr*GIp/GAr

C	DG(6,6) = Fr*GIs/GAr
C	DG(7,6) =-Fr*GIst/GAr
C	DG(8,6) = Fr*GItw/GAr
C
C	DG(7,7) = Fr*GIt/GAr
C	DG(8,7) =-Fr*GIsw/GAr
C
C	DG(8,8) = Fr*GIw/GAr
C	-------------------------------------


	DO I = 1,8
	DO J = I,8
	DG(I,J) = DG(J,I)
	ENDDO
	ENDDO


C	DO I = 1,8
C	WRITE(110,111) DG(I,1:8)
C	ENDDO
C	PAUSE
111	FORMAT(8E13.4)

	AKL = MATMUL(TRANSPOSE(BBX),MATMUL(DG,BBX))



	RETURN

	END


C	==========================================================================
C	==========================================================================
C	==========================================================================
C      AKG(4,5) =  AKG(4,5) - REL(6)
C      AKG(4,6) =  AKG(4,6) + REL(5)
C      AKG(11,12) =  AKG(11,12) - REL(13)
C      AKG(11,13) =  AKG(11,13) + REL(12)
      
      
C      AKG(5,4) =  AKG(4,5) 
C      AKG(6,4) =  AKG(4,6) 
C      AKG(12,11) =  AKG(11,12) 
C      AKG(13,11) =  AKG(11,13) 