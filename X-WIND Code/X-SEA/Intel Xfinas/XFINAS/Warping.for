C	====================================================
C	====================================================
C	====================================================
	SUBROUTINE WARPCAL(NNODE,NELEM,COORDO,LCONT,MMNN,MONST,KONST,
	1				   FKT,FKS,FKST,WARP,WARPG,TORCG,IGAS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /NUMB2/ IDOFB(9)

	COMMON /GASEC/  GAUSP(10,10),GAUSW(10,10)
	
	COMMON /SECTIONPROP/ SECTIONPROP(5000,32)
		
	DIMENSION COORDO(2,1),LCONT(8,1),LSUPP(1),WARPG(1),TORCG(1)
	DIMENSION MMNN(2,1)
	DIMENSION XY(2,8),BPG(3),BWG(3)
      DIMENSION H(8),P(2,8)
	DIMENSION XJ(2,2),XJI(2,2),STIFE(8,8)
	DIMENSION STIFG(NNODE,NNODE),STIFI(NNODE*NNODE),STDIA(NNODE)
	DIMENSION RHE(8),RHS(NNODE),TS(2),WARP(NNODE),WARS(NNODE),RHES(8)
	DIMENSION RFS(NNODE),RFT(NNODE),RFSE(8),RFTE(8)
	DIMENSION QMAT(3,3),QMAI(3,3),QRHS(3)
	DIMENSION RESUT(40)
	DIMENSION WARO(NNODE),RHEO(8)

	DIMENSION MONST(1),MCONT(8,NELEM)

	DIMENSION KKONT(NNODE),MCHEK(NNODE),COORD(2,NNODE)
	
C     FOR DGEMM FUNCTION	
	DIMENSION COUT1(8,8),COUT2(8),COUT3(3)

	ALLOCATABLE MG(:)

C	CHECK WHETHER WARPING IS INCLUDE OR NOT
	LWARP = 0
	IF(IDOFB(7).EQ.0) LWARP = 1

	DO IOD = 1,NNODE
	COORD(1:2,IOD) = COORDO(1:2,IOD)
	ENDDO

C	-------------------------------------------------------------------------------------------------
C	SCAN FOR NODE WHICH VERY CLOSE TO EACHOTHER
C
C	GET THE DATA FROM PREVIOUS CONSTRAINT VECTOR
	KKK = 0
	DO I = 1,KONST
	KKK = KKK + 1
	NG = MONST(KKK)
	DO J = 1,NG
	KKK = KKK + 1
	ENDDO
	ENDDO

C	NODE FOUND FLAG
	MCHEK(1:NNODE) = 0

	DO 112 IOD = 1,NNODE
	TI = COORD(1,IOD)
	SI = COORD(2,IOD)

	MCUNT = 1
	KKONT(MCUNT) = IOD

	DO JOD = IOD,NNODE
	TJ = COORD(1,JOD)
	SJ = COORD(2,JOD)
	SLEN = SQRT((SJ-SI)**2.0+(TJ-TI)**2.0)
	KOD = MCHEK(JOD)
	IF(SLEN.LT.1.0E-10.AND.IOD.NE.JOD.AND.KOD.NE.1) THEN
	MCUNT = MCUNT + 1    !COUNTING NUMBER OF NODE WHICH OVERLAP
	KKONT(MCUNT) = JOD   !ADD NODE TO FOUNDED NODE VECTOR
	ENDIF
	ENDDO

C	IF OVERLAP FOUND -- APPEND IT INTO PREVIOUS CONSTRAIN VECTOR
	IF(MCUNT.NE.1) THEN
	KONST = KONST + 1
	KKK = KKK + 1
	MONST(KKK) = MCUNT   !APPEND NUMBER OF NODE TO CONSTRAINT INTO CONSTRAINT VECTOR
	DO ICUNT = 1,MCUNT
	KKK = KKK + 1
	KOD = KKONT(ICUNT)   !GET NODE FROM FOUNDED NODE VECTOR
	MCHEK(KOD) = 1       !SET NODE FOUND FLAG TO BE ONE-- SO THAT IT WILL NOT LOOP AGAIN DURING SCANING FASE
	MONST(KKK) = KOD     !APPEND NODE NUMBER TO CONSTRAINT VECTOR
	ENDDO
	ENDIF

112	CONTINUE
C	-------------------------------------------------------------------------------------------------


C	================================================
C	PUT THE NODAL CONSTRAINT
C	================================================
	DO IEL = 1,NELEM
	MCONT(1:8,IEL) = LCONT(1:8,IEL)
	ENDDO

	K = 0
	DO I = 1,KONST

		K = K + 1
		NG = MONST(K)
		ALLOCATE(MG(NG))

		DO J = 1,NG
			K = K + 1
			MG(J) = MONST(K)
		ENDDO

		NPOLE = MG(1)	 
		DO IEL = 1,NELEM
			DO ICN = 1,8
			ND = MCONT(ICN,IEL)
				DO J = 1,NG
					IF(MG(J).EQ.ND) THEN
						MCONT(ICN,IEL) = NPOLE
						EXIT
					ENDIF
				ENDDO
			ENDDO
		ENDDO

		DEALLOCATE(MG)

	ENDDO
C	================================================
C	================================================


C	REFERENCE POISSON RATIO
	POISN = 0.0

C	ONE SUPPORT AT NODE ONE
	NSUPP    = 1
	LSUPP(1) = 1

C	---------------------------------------------------------
C	SETUP GAUSS POINT DATA
C	---------------------------------------------------------
	BPG(1)=  -SQRT(0.6)
	BPG(2)=   0.0D0
	BPG(3)=   SQRT(0.6)
	BWG(1)= 5.0/9.0
	BWG(2)= 8.0/9.0
	BWG(3)= 5.0/9.0



C	---------------------------------------------------------
C	LOOP TO CALCULATE SECTION PROPERTIES WRT 0,0 POINT
C	---------------------------------------------------------
	ARA = 0.0
	QS  = 0.0
	QT  = 0.0
	SIT = 0.0
	SIS = 0.0
	SST = 0.0
	DO IEL = 1,NELEM

	DO ICN = 1,8
	ND = LCONT(ICN,IEL)
	DO IXY = 1,2
	VALV = COORD(IXY,ND)
	XY(IXY,ICN) = VALV
	ENDDO
	ENDDO


	DO IGR = 1,3
	RI = BPG(IGR)
	DO IGS = 1,3
	SI = BPG(IGS)
	WG = BWG(IGR)*BWG(IGS)
	CALL WRP2D(RI,SI,H,P)
	CALL JAWRP2H(XY,P,XJ,XJI,DET)
	DVOL = WG*DET
	
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M) 
C     H(1,8),TRANSPOSE(XY)(8,2)
      ALPHA = 1.0
      BETA  = 1.0
      TS    = 0.0D0
      CALL DGEMM('N','N',1,2,8,ALPHA,H,1,TRANSPOSE(XY),8,BETA,TS,1)  	
C	TS = MATMUL(H,TRANSPOSE(XY))
	
	TT = TS(1)
	SS = TS(2)
	ARA = ARA + DVOL
	QS  = QS  + TT*DVOL
	QT  = QT  + SS*DVOL
	SIT = SIT + SS*SS*DVOL
	SIS = SIS + TT*TT*DVOL
	SST = SST + SS*TT*DVOL
	ENDDO
	ENDDO
	
	ENDDO


	QS0 = QS
	QT0 = QT
	SIT0= SIT
	SIS0= SIS
	SST0= SST
	POL0= SIT+SIS

C	---------------------------------------------------------
C	LOOP TO CALCULATE NEW COORDINATE WRT CENTROID
C	---------------------------------------------------------
	TC  = QS/ARA
	SC  = QT/ARA
	DO INODE = 1,NNODE
	COORD(1,INODE) = COORD(1,INODE) - TC 
	COORD(2,INODE) = COORD(2,INODE) - SC 
	ENDDO

	SSC = -SC
	TTC = -TC 
C	--------------------------------------------------------------------------
C	LOOP TO ASSEMBLY STIFFNESS AND RHS FORCE
C	ALAO CALCULATE SECTION PROPERTIES WRT CENTROID
C	--------------------------------------------------------------------------
	STIFG = 0.0D0
	RHS   = 0.0D0

	SIT = 0.0
	SIS = 0.0
	SST = 0.0
	ARA = 0.0
	QS  = 0.0
	QT  = 0.0
	DO IEL = 1,NELEM
	DO ICN = 1,8
	ND = LCONT(ICN,IEL)
	DO IXY = 1,2
	VALV = COORD(IXY,ND)
	XY(IXY,ICN) = VALV
	ENDDO
	ENDDO

	STIFE = 0.0D0
	RHE   = 0.0D0
C	------------------------------------------		

	DO IGR = 1,3
	RI = BPG(IGR)
	DO IGS = 1,3
	SI = BPG(IGS)
	WG = BWG(IGR)*BWG(IGS)
	CALL WRP2D(RI,SI,H,P)
	CALL JAWRP2H(XY,P,XJ,XJI,DET)
	DVOL = WG*DET
	
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M) 
C     TRANSPOSE(P)(8,2),P(2,8)
      ALPHA = 1.0
      BETA  = 1.0
      COUT1 = 0.0D0
      CALL DGEMM('N','N',8,8,2,ALPHA,TRANSPOSE(P),8,P,2,BETA,COUT1,8) 
      STIFE = STIFE + COUT1*DVOL
C	STIFE = STIFE + MATMUL(TRANSPOSE(P),P)*DVOL

C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M) 
C     H(1,8),TRANSPOSE(XY)(8,2)
      TS = 0.0D0
	CALL DGEMM('N','N',1,2,8,ALPHA,H,1,TRANSPOSE(XY),8,BETA,TS,1) 
C	TS = MATMUL(H,TRANSPOSE(XY))
	
	CALL CPU_TIME(TOEY2)
	TT = TS(1)
	SS = TS(2)
	ARA = ARA + DVOL
	SIT = SIT + SS*SS*DVOL
	SIS = SIS + TT*TT*DVOL
	SST = SST + SS*TT*DVOL
	QS  = QS  + TT*DVOL
	QT  = QT  + SS*DVOL
	TS(1) =   SS
	TS(2) =  -TT
	
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M) 	
C	TRANSPOSE(P)(8,2),TS(2)
      COUT2 = 0.0D0
      CALL DGEMM('N','N',8,1,2,ALPHA,TRANSPOSE(P),8,TS,2,BETA,COUT2,8) 
      RHE = RHE +  COUT2*DVOL
C	RHE = RHE +  MATMUL(TRANSPOSE(P),TS)*DVOL

	ENDDO
	ENDDO

C	------------------------------------------	
	DO IN = 1,8
	NI = MCONT(IN,IEL)
	RHS(NI) = RHS(NI) + RHE(IN)
	DO JN = 1,8
	NJ = MCONT(JN,IEL)
	STIFG(NI,NJ) = STIFG(NI,NJ) + STIFE(IN,JN)
	ENDDO
	ENDDO


	ENDDO


C	--------------------------------------------------------------------------
C	CALCULATE FOR SECTION PROPERTIES (PRINCIPAL AXES, DIRECTION OF BENDING)
C	--------------------------------------------------------------------------
	SICA = 0.5*(SIS + SIT)
	SSTM = SQRT(SST*SST + 0.25*(SIS - SIT)*(SIS - SIT))
	SMAX = SICA + SSTM
	SMIN = SICA - SSTM

	IF(SSTM.EQ.0.0D0) THEN
	ANGR = 0.0D0
	GOTO 100
	ENDIF
	ANGR = DASIN(SST/SSTM)

100	PII  = 3.141592654
	ANGD = ANGR*180.0D0/PII
C	--------------------------------------------------------------------------
C	LOOP TO ASSEMBLY FORCE FOR CALCULATE SHEAR COEFFICIENT
C	--------------------------------------------------------------------------
	RFS   = 0.0D0
	RFT   = 0.0D0
	DO IEL = 1,NELEM

	DO ICN = 1,8
	ND = LCONT(ICN,IEL)
	DO IXY = 1,2
	VALV = COORD(IXY,ND)
	XY(IXY,ICN) = VALV
	ENDDO
	ENDDO

	RFSE   = 0.0D0
	RFTE   = 0.0D0
C	------------------------------------------	
	DO IGR = 1,3
	RI = BPG(IGR)
	DO IGS = 1,3
	SI = BPG(IGS)
	WG = BWG(IGR)*BWG(IGS)
	CALL WRP2D(RI,SI,H,P)
	CALL JAWRP2H(XY,P,XJ,XJI,DET)
	DVOL = WG*DET
	
C     H(1,8),TRANSPOSE(XY)(8,2)
      TS = 0.0D0
	CALL DGEMM('N','N',1,2,8,ALPHA,H,1,TRANSPOSE(XY),8,BETA,TS,1) 
C	TS = MATMUL(H,TRANSPOSE(XY))

      CALL SHEARFS(RFSE,H,P,TS,DVOL,POISN,SIS,SIT,SST,8)
      CALL SHEARFT(RFTE,H,P,TS,DVOL,POISN,SIS,SIT,SST,8)


	ENDDO
	ENDDO
C	------------------------------------------	

	DO IN = 1,8
	NI = MCONT(IN,IEL)
	RFS(NI) = RFS(NI) + RFSE(IN)
	RFT(NI) = RFT(NI) + RFTE(IN)
	ENDDO


	ENDDO


C	---------------------------------------------------------
C	ASSIGN SUPPORT CONDITION
C	---------------------------------------------------------
	DO ISUPP = 1,NSUPP
	LCOL = LSUPP(ISUPP)
	RHS(LCOL) = 0.0D0
	RFS(LCOL) = 0.0D0
	RFT(LCOL) = 0.0D0
	DO I = 1,NNODE
	STIFG(I,LCOL) = 0.0D0
	STIFG(LCOL,I) = 0.0D0
	ENDDO
	STIFG(LCOL,LCOL) = 1.0D0
	ENDDO


C	---------------------------------------------------------
C	DETECT ZERO DIAGONAL
C	---------------------------------------------------------
	DO INODE = 1,NNODE
	IF(STIFG(INODE,INODE).EQ.0.0D0) THEN
	STIFG(INODE,INODE) = 1.0D0
	RHS(INODE) = 0.0D0
	RFS(INODE) = 0.0D0
	RFT(INODE) = 0.0D0
	ENDIF
	ENDDO

C	---------------------------------------------------------
C	INVERSE TO GET NODAL VALUE OF 
C	WARPING FUNCTION -- WARP
C	S SHEAR FUNCTION -- RFS
C	T SHEAR FUNCTION -- RFT
C	---------------------------------------------------------

C	CALL INVMATRIX(STIFG,STIFI,NNODE)
C	WARP = MATMUL(STIFI,RHS)
C	RFS  = MATMUL(STIFI,RFS)
C	RFT  = MATMUL(STIFI,RFT)

	K = 0	
	DO I = 1,NNODE
	JJ = I
	DO J = 1,I
	K = K + 1
	STIFI(K) = STIFG(I,JJ)
	JJ = JJ - 1
	ENDDO
	STDIA(I) = 0.0D0
	ENDDO

	CALL COLSYM (STIFI,STDIA,WARP,1,NNODE,2)

	WARP(1:NNODE) = RHS(1:NNODE)
	CALL COLSYM (STIFI,STDIA,WARP,2,NNODE,2)
	CALL COLSYM (STIFI,STDIA,RFS,2,NNODE,2)
	CALL COLSYM (STIFI,STDIA,RFT,2,NNODE,2)


C	================================================
C	PUT THE NODAL CONSTRAINT
C	================================================

	K = 0
	DO I = 1,KONST

		K = K + 1
		NG = MONST(K)
		ALLOCATE(MG(NG))

		DO J = 1,NG
			K = K + 1
			MG(J) = MONST(K)
		ENDDO

		NPOLE = MG(1)	 
		DO INODE = 1,NNODE
			DO J = 1,NG
				IF(MG(J).EQ.INODE) THEN
					WARP(INODE) = WARP(NPOLE)
					 RFS(INODE) =  RFS(NPOLE)
					 RFT(INODE) =  RFT(NPOLE)
					EXIT
				ENDIF
			ENDDO
		ENDDO

		DEALLOCATE(MG)

	ENDDO
C	================================================
C	================================================

C	---------------------------------------------------------
C	CALCULATE TORSIONAL CONSTANT  J = Is + It - WT*R
C	---------------------------------------------------------
	TORC = SIT+SIS
	DO INODE = 1,NNODE
	TORC =TORC - WARP(INODE)*RHS(INODE)
	ENDDO

C	---------------------------------------------------------
C	INTEGRAL THE WARPING TERM OVER THE CROSS-SECTION
C	QW  = INTG(W   )dA
C	SWS = INTG(W*TT)dA
C	SWT = INTG(W*SS)dA
C	INTEGRAL THE SHEAR FUNCTION OVER THE CROSS-SECTION
C	---------------------------------------------------------
	QW  = 0.0D0 ; SWS = 0.0D0 ; SWT = 0.0D0		!WARPING
	QFS = 0.0D0 ; QFT = 0.0D0	 				!SHEAR FUNCTION WEIGHTING
	FKT = 0.0D0 ; FKS = 0.0D0 ; FKST= 0.0D0		!SHEAR DEFORM COEF.
	TSHE= 0.0D0 ; SSHE= 0.0D0  					!SHEAR CENTER
	DO IEL = 1,NELEM

	DO ICN = 1,8
	ND = LCONT(ICN,IEL)
	DO IXY = 1,2
	VALV = COORD(IXY,ND)
	XY(IXY,ICN) = VALV
	ENDDO
	ND = MCONT(ICN,IEL)
	 RHE(ICN) = WARP(ND)
	RFSE(ICN) =  RFS(ND)
	RFTE(ICN) =  RFT(ND)
	ENDDO
	
C	------------------------------------------
	
	DO IGR = 1,3
	RI = BPG(IGR)
	DO IGS = 1,3
	SI = BPG(IGS)
	WG = BWG(IGR)*BWG(IGS)
	CALL WRP2D(RI,SI,H,P)
	CALL JAWRP2H(XY,P,XJ,XJI,DET)
	
C     H(1,8),TRANSPOSE(XY)(8,2)
      TS = 0.0D0
	CALL DGEMM('N','N',1,2,8,ALPHA,H,1,TRANSPOSE(XY),8,BETA,TS,1) 
C	TS = MATMUL(H,TRANSPOSE(XY))

	TT = TS(1)
	SS = TS(2)
	DVOL = WG*DET

	WW = 0.0D0
	FSE= 0.0D0
	FTE= 0.0D0
	DO II = 1,8
	WW = WW + RHE(II)*H(II)
	FSE= FSE+RFSE(II)*H(II)
	FTE= FTE+RFTE(II)*H(II)
	ENDDO
	QW  =  QW +    WW*DVOL
	SWS = SWS + TT*WW*DVOL
	SWT = SWT + SS*WW*DVOL

	QFS  =QFS +   FSE*DVOL
	QFT  =QFT +   FTE*DVOL

	CALL SHFACT(FKT,FKS,FKST,TSHE,SSHE,RFTE,RFSE,H,P,
	1			TS,DVOL,ARA,POISN,SIS,SIT,SST,8)
	
	ENDDO
	ENDDO


	ENDDO


C	---------------------------------------------------------
C	SOLVE 3x3 MATRIX TO GET THE THREE CONSTANT
C	SSH = S SHEAR CENTER WRT CENTROID
C	TSH = T SHEAR CENTER WRT CENTROID
C	 CW = CONSTANT TO TRANSFORM WARPING FUNCTION FROM CENTROID TO SHEAR CENTER
C	---------------------------------------------------------
	QMAT(1,1:3) = [-QS , QT , ARA]
	QMAT(2,1:3) = [-SIS, SST, QS ]
	QMAT(3,1:3) = [-SST, SIT, QT ]
	QRHS(1:3)   = [-QW ,-SWS,-SWT]
	CALL INVMATRIX(QMAT,QMAI,3)
	
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M) 	
C	QMAI(3,3),QRHS(3,1)
      COUT3 = 0.0D0
      CALL DGEMM('N','N',3,1,3,ALPHA,QMAI,3,QRHS,3,BETA,COUT3,3)
      QRHS = COUT3
C	QRHS = MATMUL(QMAI,QRHS)

	SSH = QRHS(1)
	TSH = QRHS(2)
	CW  = QRHS(3)

C	---------------------------------------------------------
C	COMPUTE WARPING FUNCTION AT SHEAR CENTER WARS
C	Ws = Wc - SSH*TT + TSH*SS + CW
C	AND UPDATE THE WARPING FUNTION AT CENTROID TO MAKE A ZERO INTERGRAL Wc = Wc+CW
C	---------------------------------------------------------
	DO I = 1,NNODE
	TT = COORD(1,I)
	SS = COORD(2,I)
	WARS(I) = WARP(I) - SSH*TT + TSH*SS + CW
	WARO(I) = WARP(I) - SSC*TT + TTC*SS + CW
	WARP(I) = WARP(I) + CW

	 RFS(I) =  RFS(I) - QFS/ARA
	 RFT(I) =  RFT(I) - QFT/ARA
	ENDDO

C	---------------------------------------------------------
C	CALCULATE WARPING CONSTANT Iw = INTG(W*W)dA
C	---------------------------------------------------------
	SIWC = 0.0D0
	SIWS = 0.0D0
	SIWO = 0.0D0
	DO IEL = 1,NELEM

	DO ICN = 1,8
	ND = LCONT(ICN,IEL)
	DO IXY = 1,2
	VALV = COORD(IXY,ND)
	XY(IXY,ICN) = VALV
	ENDDO
	ND = MCONT(ICN,IEL)
	 RHE(ICN) = WARP(ND)
	RHES(ICN) = WARS(ND)
	RHEO(ICN) = WARO(ND)
	ENDDO
	
C	------------------------------------------
	
	DO IGR = 1,3
	RI = BPG(IGR)
	DO IGS = 1,3
	SI = BPG(IGS)
	WG = BWG(IGR)*BWG(IGS)
	CALL WRP2D(RI,SI,H,P)
	CALL JAWRP2H(XY,P,XJ,XJI,DET)
	DVOL = WG*DET

	WW = 0.0D0
	WWS= 0.0D0
	WWO= 0.0D0
	DO II = 1,8
	WW = WW +  RHE(II)*H(II)
	WWS= WWS+ RHES(II)*H(II)
	WWO= WWO+ RHEO(II)*H(II)
	ENDDO
	SIWC = SIWC + WW *WW *DVOL
	SIWS = SIWS + WWS*WWS*DVOL
	SIWO = SIWO + WWO*WWO*DVOL

	ENDDO
	ENDDO
	

	ENDDO

C	---------------------------------------------------------
C	PRINT TO OUTPUT
C	---------------------------------------------------------
	RESUT(1 :5 ) = [ARA,QS0,QT0,SC,TC]
	RESUT(6 :9 ) = [SSHE,TSHE,SSH,TSH]
	RESUT(10:13) = [SIS0,SIT0,SST0,SIS0+SIT0]
	RESUT(14:17) = [SIS,SIT,SST,SIS+SIT]
	RESUT(18:22) = [ANGR,ANGD,SMAX,SMIN,SSTM]
	RESUT(23:24) = [0.0,0.0]
	RESUT(25   ) = POISN

	FKTI = 0.0D0
	FKSI = 0.0D0
	FKSTI = 0.0D0

	IF(FKT.GT.0.0D0) FKTI = 1.0/FKT
	IF(FKS.GT.0.0D0) FKSI = 1.0/FKS
	IF(FKST.GT.0.0D0) FKSTI = 1.0/FKST

	RESUT(26:28) = [FKSI,FKTI,FKSTI]
	RESUT(29:32) = [TORC,SIWS,SIWC,SIWO]

	DO I = 1,32
	IF(ABS(RESUT(I)).LT.1.0E-12) RESUT(I) = 0.0D0
      ENDDO
      
      IF (IGAS.EQ.8)THEN
      TOEY = 1
      ENDIF

      SECTIONPROP(IGAS,1:32) = RESUT(1:32) 
      
	WRITE(110,500) RESUT(1:32)
C	WRITE(ITO,500) RESUT(1:32)


500	FORMAT(/X,'Section Properties',//
	1		X,'Cross-Sectional Area. . . . . . . . . ',E13.6/		!1
	2		X,'Moment of Area Qs . . . (Origin). . . ',E13.6/
	3		X,'Moment of Area Qt . . . (Origin). . . ',E13.6//
	4		X,'S Centroid. . . . . . . (Origin). . . ',E13.6/
	5		X,'T Centroid. . . . . . . (Origin). . . ',E13.6//		!5
	6		X,'S Shear Center. . . . . (Centroid). . ',E13.6/		!6
	7		X,'T Shear Center. . . . . (Centroid). . ',E13.6/
	8		X,'S Shear Center-Trefftz. (Centroid). . ',E13.6/
	9		X,'T Shear Center-Trefftz. (Centroid). . ',E13.6//		!9
	1		X,'Moment  of Inertia Is . (Origin). . . ',E13.6/		!10
	2		X,'Moment  of Inertia It . (Origin). . . ',E13.6/
	3		X,'Product of Inertia Ist. (Origin). . . ',E13.6/
	4		X,'Polar Moment of Inertia (Origin). . . ',E13.6//		!13
	5		X,'Moment  of Inertia Is . (Centroid). . ',E13.6/		!14
	6		X,'Moment  of Inertia It . (Centroid). . ',E13.6/
	7		X,'Product of Inertia Ist. (Centroid). . ',E13.6/
	8		X,'Polar Moment of Inertia.(Centroid). . ',E13.6//		!17
	9		X,'Princ Bending Angle-rad (Centroid). . ',E13.6/		!18
	1		X,'Princ Bending Angle-deg (Centroid). . ',E13.6/
	2		X,'Princ Moment of Int Max (Centroid). . ',E13.6/
	3		X,'Princ Moment of Int Min (Centroid). . ',E13.6/		
	4		X,'Princ Prd of Inetia Max (Centroid). . ',E13.6//		!22
	5		X,'S Norminal Size . . . . . . . . . . . ',E13.6/		!23
	6		X,'T Norminal Size . . . . . . . . . . . ',E13.6//		!24
	7		X,'Reference Poisson Ratio . . . . . . . ',E13.6//		!25
	8		X,'S  Shear Coefficient. (Shear Center). ',E13.6/		!26
	9		X,'T  Shear Coefficient. (Shear Center). ',E13.6/
	1		X,'ST Shear Coefficient. (Shear Center). ',E13.6//		!28
	2		X,'Torsional Constant J. . . . . . . . . ',E13.6/		!29
	3		X,'Warping Constant . . .(Shear Center). ',E13.6/
	4		X,'Warping Constant . . . .(Centroid). . ',E13.6/  		!31
	4		X,'Warping Constant . . . . (Offset) . . ',E13.6//)		!32


C	IF(NNODE.GT.2000) THEN
C	---------------------------------------------------------
	
C	WRITE(110,700 ) 
C	WRITE(110,800 ) 
C	WRITE(110,900 ) 
C	WRITE(110,1000) 

C	DO I = 1,NNODE
C	WRITE(110,600) I,WARS(I),WARP(I),CW
C	ENDDO
C	WRITE(110,1100) 

C	---------------------------------------------------------

C	WRITE(110,1200) 
C	WRITE(110,1300) 
C	WRITE(110,1000) 

C	DO I = 1,NNODE
C	WRITE(110,600) I,RFS(I),RFT(I),0.0
C	ENDDO
C	WRITE(110,1100) 

C	---------------------------------------------------------
C	STOP
C	ENDIF

C	---------------------------------------------------------
C	ELEMENT LOOP
C	---------------------------------------------------------
	IFIB = 0
	DO IEL = 1,NELEM

C	LOOP OVER ELEMENt NODES
	DO ICN = 1,8
	ND = LCONT(ICN,IEL)
	DO IXY = 1,2
	VALV = COORD(IXY,ND)
	XY(IXY,ICN) = VALV
	ENDDO
	ND = MCONT(ICN,IEL)
	 RHE(ICN) = WARO(ND)        !USE WARPING AT OFFSET POINt !SHEAR CENTER
	ENDDO
	
	MM = MMNN(1,IEL)
	NN = MMNN(2,IEL)

C	GAUSS LOOP	
	DO IGR = 1,MM
	RI = GAUSP(IGR,MM)
	DO IGS = 1,NN
	SI = GAUSP(IGS,NN)
	WG = GAUSW(IGR,MM)*GAUSW(IGS,NN)

	CALL WRP2D(RI,SI,H,P)
	CALL JAWRP2H(XY,P,XJ,XJI,DET)	

C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M) 
C     H(1,8),TRANSPOSE(XY)(8,2)
      TS = 0.0D0
	CALL DGEMM('N','N',1,2,8,ALPHA,H,1,TRANSPOSE(XY),8,BETA,TS,1) 
C	TS = MATMUL(H,TRANSPOSE(XY))
C	TS = MATMUL(H,TRANSPOSE(XY))
	
	TT = TS(1)
	SS = TS(2)
	TS(1) =   SS
	TS(2) =  -TT

	IFIB = IFIB + 1
	WARPG(IFIB) = 0.0D0
	TORCG(IFIB) = 0.0D0

	WPT = 0.0D0
	WPS = 0.0D0
	DO II = 1,8
	IF(LWARP.NE.0) WARPG(IFIB) = WARPG(IFIB) +  RHE(II)*H(II)
	WPT         = WPT + RHE(II)*P(1,II)
	WPS         = WPS + RHE(II)*P(2,II)
	ENDDO
	TORCG(IFIB) = TS(1)*WPT + TS(2)*WPS

	ENDDO
	ENDDO
C	END GAUSS LOOP

	ENDDO
C	---------------------------------------------------------
C	END ELEMENT LOOP
C	---------------------------------------------------------

C	RETURN THE VALUE OF WARPING AT CONTOUR POINT
	WARP(1:NNODE) = 0.0D0
	IF(LWARP.NE.0) THEN
	DO IOD = 1,NNODE
	WARP(IOD) = WARO(IOD)
	ENDDO
	ENDIF



600	FORMAT(I5,X,100E15.6)

700	FORMAT('GiD Post Results File 1.0'//)

800	FORMAT('Result "Warping" "Section Analysis" 1 Vector OnNodes')

900	FORMAT('ComponentNames "Warping-Shear-Center" ', 
	1					  '"Warping-Centroid" ', 
     2					  '"Warping-Centroid-Constant"')

1000	FORMAT('Values')

1100	FORMAT('End Values'//)


1200	FORMAT('Result "Shear" "Section Analysis" 1 Vector OnNodes')

1300	FORMAT('ComponentNames "S-Shear-Func" ', 
	1					  '"T-Shear-Func" ', 
     2					  '"------------"')


	RETURN

	END
C	====================================================
C	====================================================
C	====================================================
      SUBROUTINE SHEARFT(RFT,H,P,TS,DVOL,POISN,SIS,SIT,SST,NNM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	
	DIMENSION RFT(1)
      DIMENSION H(1),P(2,1),TS(1)

	TT = TS(1)
	SS = TS(2)
	R  = TT*TT - SS*SS
	Q  = 2.0D0*TT*SS

	D1 = 0.5*POISN*(SIT*R - SST*Q)
	D2 = 0.5*POISN*(SST*R + SIT*Q)
	D3 = 2.0*(1.0+POISN)*(SIT*TT-SST*SS)

	DO INM = 1,NNM
	RFT(INM) = RFT(INM) + P(1,INM)*D1*DVOL
	RFT(INM) = RFT(INM) + P(2,INM)*D2*DVOL
	RFT(INM) = RFT(INM) + H(  INM)*D3*DVOL
	ENDDO
	

      RETURN
      END
C
C	====================================================
C	====================================================
C	====================================================
      SUBROUTINE SHEARFS(RFS,H,P,TS,DVOL,POISN,SIS,SIT,SST,NNM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	
	DIMENSION RFS(1)
      DIMENSION H(1),P(2,1),TS(1)

	TT = TS(1)
	SS = TS(2)
	R  = TT*TT - SS*SS
	Q  = 2.0D0*TT*SS

	H1 = 0.5*POISN*(-SST*R + SIS*Q)
	H2 = 0.5*POISN*(-SIS*R - SST*Q)
	H3 = 2.0*(1.0+POISN)*(SIS*SS-SST*TT)

	DO INM = 1,NNM
	RFS(INM) = RFS(INM) + P(1,INM)*H1*DVOL
	RFS(INM) = RFS(INM) + P(2,INM)*H2*DVOL
	RFS(INM) = RFS(INM) + H(  INM)*H3*DVOL
	ENDDO
	

      RETURN
      END
C
C	====================================================
C	====================================================
C	====================================================
      SUBROUTINE SHFACT(FKT,FKS,FKST,TSHE,SSHE,RFT,RFS,H,P,
	1				  TS,DVOL,ARA,POISN,SIS,SIT,SST,NNM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	
	DIMENSION RFT(1),RFS(1)
      DIMENSION H(1),P(2,1),TS(1)

	TT = TS(1)
	SS = TS(2)

	DELTA  = 2.0*(1.0+POISN)*(SIS*SIT-SST*SST)
	DELTA2 = DELTA*DELTA 


	R  = TT*TT - SS*SS
	Q  = 2.0D0*TT*SS
	D1 = 0.5*POISN*(SIT*R - SST*Q)
	D2 = 0.5*POISN*(SST*R + SIT*Q)

	R  = TT*TT - SS*SS
	Q  = 2.0D0*TT*SS
	H1 = 0.5*POISN*(-SST*R + SIS*Q)
	H2 = 0.5*POISN*(-SIS*R - SST*Q)

	PT1 = -D1
	PT2 = -D2
	PS1 = -H1
	PS2 = -H2
	DO INM = 1,NNM
	PT1 = PT1 + P(1,INM)*RFT(INM)
	PT2 = PT2 + P(2,INM)*RFT(INM)
	PS1 = PS1 + P(1,INM)*RFS(INM)
	PS2 = PS2 + P(2,INM)*RFS(INM)
	ENDDO

C	SHEAR DEFORMATION COEF.
	FKT  = FKT + (PT1*PT1 + PT2*PT2)*DVOL*ARA/DELTA2	!ALONG T
	FKS  = FKS + (PS1*PS1 + PS2*PS2)*DVOL*ARA/DELTA2	!ALONG S
	FKST = FKST+ (PT1*PS1 + PT2*PS2)*DVOL*ARA/DELTA2	!INTERACTION ST
	
C	SHEAR CENTER ALONG T
	TSHE = TSHE+ (SIS*TT+SST*SS)*(TT*TT+SS*SS)*DVOL*0.5*POISN/DELTA
	TSHE = TSHE- (SS*PS1 - TT*PS2)*DVOL/DELTA

C	SHEAR CENTER ALONG S
	SSHE = SSHE+ (SIT*SS+SST*TT)*(TT*TT+SS*SS)*DVOL*0.5*POISN/DELTA
	SSHE = SSHE+ (SS*PT1 - TT*PT2)*DVOL/DELTA

      RETURN
      END
C
C	====================================================
C	====================================================
C	====================================================
      SUBROUTINE WRP2D (R,S,HM,PM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO FIND INTERPOLATION FUNCTIONS AND THERE DERIVATIVES
C     AT THE NODAL POINTS OF A 4 TO 8 NODE,ISOPARAMETRIC QUADRILATERAL
C	----------------------------------------------------------------
C                         NODE NUMBERING CONVENTION
C                         -------------------------
C
C                   2                 5                 1
C
C                     0 . . . . . . . 0 . . . . . . . 0
C                     .                               .
C                     .                               .
C                     .               S               .
C                     .               .               .
C                     .               .               .
C                   6 0               . . . R         0 8
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     .                               .
C                     0 . . . . . . . 0 . . . . . . . 0
C
C                   3                 7                 4
C
C     R,S        = NATURAL COORDINATES OF POINT TO BE INTERPOLATED
C     H(8)       = INTERPOLATION (SHAPE) FUNCTIONS
C     P(2,8)     = FUNCTION DERIVATIVES WITH RESPECT TO R,S RESP.
C     NODEX(NEX) = POSITION OF MIDSIDE (EXCESSIVE) NODES
C     NNO        = NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ----------------------------------------------------------------
      DIMENSION  H(8),P(2,8),IPERM(4),NODEX(4),HM(8),PM(2,8),MAPP(8)
      DATA (IPERM(I), I=1,4)  /2,3,4,1/

	NNO = 8
	NODEX(1:4) = [5,6,7,8]
C
      NMI = NNO-4
      RP  = 1.0+R
      SP  = 1.0+S
      RM  = 1.0-R
      SM  = 1.0-S
      R2  = 1.0-R*R
      S2  = 1.0-S*S
C     ---------------------------------------------
C     INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
C     FOR A FOUR NODE ELEMENT
C     ---------------------------------------------
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
      IF (NNO.EQ.4)  RETURN
C     -------------------------------------
C     ADD DEGREES OF FREEDOM IN EXCESS OF 4
C     -------------------------------------
      DO 100  IMI=1,NMI
      NN = NODEX(IMI)-4
      GOTO (50,60,70,80), NN
C      CALL GOTOER!NON-EXIST SUBROUTINE
 50   H(5)   = 0.5*R2*SP
      P(1,5) = -R*SP
      P(2,5) = 0.5*R2
      GOTO 100
 60   H(6)   = 0.5*RM*S2
      P(1,6) = -0.5*S2
      P(2,6) = -RM*S
      GOTO 100
 70   H(7)   = 0.5*R2*SM
      P(1,7) = -R*SM
      P(2,7) = -0.5*R2
      GOTO 100
 80   H(8)   = 0.5*RP*S2
      P(1,8) = 0.5*S2
      P(2,8) = -RP*S
 100  CONTINUE
C     ----------------------------------------------
C     CORRECT FUNCTIONS AND DERIVATIVES IF 5 OR MORE
C     NODES ARE USED TO DESCRIBE THE ELENENT
C     ----------------------------------------------
      DO 200  IMI=1,NMI
      IN = NODEX(IMI)
      I1 = IN-4
      I2 = IPERM(I1)
      H(I1)      = H(I1)-0.5*H(IN)
      H(I2)      = H(I2)-0.5*H(IN)
      H(IMI+4)   = H(IN)
      DO 200  J=1,2
      P(J,I1)    = P(J,I1)-0.5*P(J,IN)
      P(J,I2)    = P(J,I2)-0.5*P(J,IN)
 200  P(J,IMI+4) = P(J,IN)
C

	MAPP(1:8) = [1,3,5,7,2,4,6,8]

	DO I = 1,8
	NN= MAPP(I)
	HM(NN) = H(I)
	PM(1,NN) = P(1,I)
	PM(2,NN) = P(2,I)
	ENDDO


      RETURN
      END

C	====================================================
C	====================================================
C	====================================================
	SUBROUTINE JAWRP2H(XY,P,XJ,XJI,DET)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------------
C     FINDS JACOBIAN (XJ), ITS DETERMINANT (DET)
C     AND THE INVERSE (XJI) OF THE JACOBIAN
C     ------------------------------------------
      DIMENSION XY(2,1),P(2,8),XJ(2,2),XJI(2,2)
      
C     FOR DGEMM FUNCTION	
      DIMENSION COUT1(2,8)
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
	NNO = 8
      DO 100  I=1,2
      DO 100  J=1,2
      DUM = 0.0
      DO 90   K=1,NNO
 90   DUM = DUM + P(I,K)*XY(J,K)
 100  XJ(I,J) = DUM
C     ---------------------------------
C     DETERMINANT (DET) OF THE JACOBIAN
C     ---------------------------------
      DET = XJ(1,1)*XJ(2,2) - XJ(2,1)*XJ(1,2)

	IF(ABS(DET).LE.1.0E-12) THEN
	P = 0.0D0
	DET = 0.0D0
	XJ = 0.0D0
	XJI = 0.0D0
	RETURN
	ENDIF

C     -----------------------------
C     INVERSE (XJI) OF THE JACOBIAN
C     -----------------------------
      DUM = 1.0/DET
      XJI(1,1) =  XJ(2,2)*DUM
      XJI(2,1) = -XJ(2,1)*DUM
      XJI(1,2) = -XJ(1,2)*DUM
      XJI(2,2) =  XJ(1,1)*DUM
      
C     CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)  
C     XJI(2,2),P(2,8)
      ALPHA = 1.0D0
      BETA  = 1.0D0
      COUT1 = 0.0D0
      CALL DGEMM('N','N',2,8,2,ALPHA,XJI,2,P,2,BETA,COUT1,2) 
      P = COUT1
C	P = MATMUL(XJI,P)
C
      RETURN
      END
C
C	====================================================
C	====================================================
C	====================================================
 