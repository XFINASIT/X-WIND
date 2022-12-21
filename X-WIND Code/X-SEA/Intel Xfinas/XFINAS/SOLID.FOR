C	=====================================================
C	6 KINDS OF SOLID ELEMENTS (ITYPE=10 AND ISTYP=1---6)
C	CREADTED BY WANG, OCTOBER 6,2002
C	=====================================================
      SUBROUTINE SOLID3 (PROPM,PROPG,NODEX,WA,WA2,S,COORD,EDIS,EDISI,RE,MWG,
     +           ALPHA,SEL,SEDI,FIN,HINFC,MSET,AMV)



	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ------------------------------------
C     CALLS THE APPROPRIATE SOLID ELEMENT MODULE 
C	ISTYP=1:	ASSUMED NATURAL 8-NODE SOLID ELEMENT (HEXAHEDRON)
C	ISTYP=2:	ENHANCED 8-NODE SOLID ELEMENT   (HEXAHEDRON-LIU)
C	ISTYP=3:	ENHANCED 8-NODE SOLID ELEMENT   (HEXAHEDRON-WANG)
C	ISTYP=4:	ENHANCED 6-NODE SOLID ELEMENT (PENTRAHEDRON)
C	ISTYP=5:	HEAT TRANSFER ELEMENT
C	ISTYP=6:	ENHANCED 4-NODE SOLID ELEMENT	(TETRAHEDRON)
C	ISTYP=7:	COMPATIBLE 8-NODE SOLID ELEMENT (HEXAHEDRON)
C     ------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT

C	EAS ARRAY SONGSAK MAR2006
	COMMON /MMENH/ MM,MM1,MM2,NDIMC


      DIMENSION PROPM(1),PROPG(1),NODEX(1),WA(1),S(1),COORD(1),EDIS(1)
      DIMENSION EDISI(1),ELOD(1),FIN(1),AMV(3)
      
      DIMENSION WA2(1) !---------NEW WORKING ARRAY BY BJ
	

	DIMENSION ALPHA(MM,1),SEL(MM,24),SEDI(MM,MM),RH(MM,1),
	1		  HINFC(MM)

C 
	GOTO (11,12,13,14,15,16,17,18,19,20,21,22,23), ISTYP

 11	CALL SOLICO(PROPM,PROPG,NODEX,WA,WA2,S,COORD,EDIS,EDISI,RE
     1           ,MWG,ALPHA,SEL,SEDI,FIN,HINFC) ! Pramin Solid-Shell
	RETURN	

 12	CALL SOLID82(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,MWG
     +				,ALPHA,SEL,SEDI,FIN) !LIU EAS + ANS 5 TERM
	RETURN

 13	CALL SOLI632(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,
     +			 MWG,ALPHA,SEL,SEDI,FIN,HINFC) !wang
      RETURN

 14	CALL SOLINEW(PROPM,PROPG,NODEX,WA,WA2,S,COORD,EDIS,EDISI,RE,MWG,
	1			 FIN,MSET,SEDI,SEL,ALPHA,HINFC,AMV) !LIU SOLID SHELL
 	RETURN

C	HEAT TRANSFER ELEMENT
C	---------------------
 15	CONTINUE
	RETURN
	 
 16	CALL SOLIDEAS(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE
     +			 ,MWG,ALPHA,SEL,SEDI,FIN,HINFC) ! Sucharuck AIT

	RETURN !SOLID 85

 17	CONTINUE
	RETURN 

 18   CALL SOLI81(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,MWG)
	RETURN !STANDARD COMPATIBLE SOLID

 19	CONTINUE
C	CALL SOLITHIN(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE
C     +			 ,MWG,FIN) ! Thin Solid Joint Element Songsak FEB2006
	CALL SOLIDINF(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE
     +			 ,MWG,FIN) ! Contact Zero Thickness Interface Element Songsak FEB2006
	RETURN 

20	CONTINUE
	RETURN 

21	CONTINUE
	RETURN 

22	CONTINUE
	RETURN

23	CALL SOLIDEL3T(PROPM,PROPG,WA,S,COORD,EDIS,EDISI,RE,MWG,FIN)
	RETURN !TETRAHEDRA 



	RETURN
	END

C=====================================================================
C	START 3D 8-NODE COMPATIBLE SOLID ELEMENT 
C=====================================================================
      SUBROUTINE SOLI81(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,
	1				  MWG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     MAIN PROGRAM FOR THE 3-D SOLID
C     EVALUATES THE TANGENTIAL STIFFNESS MATRIX,STRAINS AND STRESSES
C     AND EQUILIBRIUM LOADS FOR THE CURVILINEAR ISOPARAMETRIC
C     HEXAHEDRON (8 TO 21 NODES)
C     ----------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
	COMMON /THERM/ STRAINO(6)

C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(*),COORD(*)
      DIMENSION EDIS(*),EDISI(*),RE(*),XJ(9)
      DIMENSION DP(6,6),H(21),P(3,21),XJI(9),B(63),DISD(9)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6)

C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------
      CALL HOKLAW (PROPM,PROPG,1)
      MGR = NGR
      MGS = NGS
      MGT = NGT
      IF (ITASK.NE.5) GOTO 10
      MGR = 3
      MGS = 3
      MGT = 3
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
 10   IPT = 0
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
      CALL SHAP3D (RI,SI,TI,H,P,NODEX,NNO)

	
      CALL JACO3D (COORD,P,XJ,XJI,DET,MEL,NNO)
      DVOL = WT*DET
C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
      IF (ITASK.NE.5)  GOTO 50
      CALL SOMASS (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
      GOTO 900
C     ----------------------------------------
C     COMPACTED STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------
 50   CALL SOBMAT (P,XJI,B,NNO)
      IF (NLOPT+ITASK.EQ.1)  GOTO 700
C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------
      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,6
 100  TAU(I) = WA(I,IPT)*DVOL
      GOTO 800
C     -------------------------------
C     DISPLACEMENT DERIVATIVES (DISD)
C     -------------------------------
 200  CALL CLEARA (DISD,9)
      DO 290  IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDIS(IEF)
      DISD(2) = DISD(2) + B(JEF)*EDIS(JEF)
      DISD(3) = DISD(3) + B(KEF)*EDIS(KEF)
      DISD(4) = DISD(4) + B(JEF)*EDIS(IEF)
      DISD(5) = DISD(5) + B(KEF)*EDIS(IEF)
      DISD(6) = DISD(6) + B(IEF)*EDIS(JEF)
      DISD(7) = DISD(7) + B(KEF)*EDIS(JEF)
      DISD(8) = DISD(8) + B(IEF)*EDIS(KEF)
 290  DISD(9) = DISD(9) + B(JEF)*EDIS(KEF)
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)
      STRAIN(4) = DISD(4) + DISD(6)
      STRAIN(5) = DISD(5) + DISD(8)
      STRAIN(6) = DISD(7) + DISD(9)
C	---------------
C	FOR THERMAL STRAIN (INITIAL STRAIN)
C	---------------
	DO 292 I=1,6
 292	STRAIN(I)=STRAIN(I)-STRAINO(I)
C     -------------------------------------------------------------
C     FOR NLOPT>1 SUBTRACT QUADRATIC STRAIN TERMS (ALMANSI STRAINS)
C     -------------------------------------------------------------
	IF (NLOPT.LE.1) GOTO 395
      QSTRAI(1) = .5*(DISD(1)*DISD(1)+DISD(6)*DISD(6)+DISD(8)*DISD(8))
      QSTRAI(2) = .5*(DISD(4)*DISD(4)+DISD(2)*DISD(2)+DISD(9)*DISD(9))
      QSTRAI(3) = .5*(DISD(5)*DISD(5)+DISD(7)*DISD(7)+DISD(3)*DISD(3))
      QSTRAI(4) =     DISD(1)*DISD(4)+DISD(6)*DISD(2)+DISD(8)*DISD(9)
      QSTRAI(5) =     DISD(1)*DISD(5)+DISD(6)*DISD(7)+DISD(8)*DISD(3)
      QSTRAI(6) =     DISD(4)*DISD(5)+DISD(2)*DISD(7)+DISD(9)*DISD(3)
      DO 390  I=1,6
 390  STRAIN(I) = STRAIN(I) - QSTRAI(I)
C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------
 395  IF(MTMOD.EQ.6) GO TO 455
 400  IF (MTMOD-3)  405,450,460
 405  CALL SOLSIG (STRAIN,STRESS)
      DO 410  I=1,6
 410  WA(I,IPT) = STRESS(I)
      GOTO 500
C
 450  CALL MISE3D (WA(1,IPT),WA(7,IPT),WA(13,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DP)
      GOTO 500

 455	CALL DAMAGE(STRAIN,STRESS,DP,PROPM,WA(1,IPT))

	DO II=1,6
	WA(II,IPT)   = STRESS(II)
	WA(II+6,IPT) = STRAIN(II)
	ENDDO
	GOTO 500


C460  CALL CONC3D -- DELETED BECAUSE OF SUBROUTINE NOT EXIST
 460  CONTINUE

	
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510  I=1,6
 510  TAU(I) = STRESS(I)*DVOL
      IF (ITASK.LE.2) GOTO 520
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900
C     ---------
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     ---------
 520	DO 550  IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      RE(IEF)=RE(IEF) + B(IEF)*TAU(1) + B(JEF)*TAU(4) + B(KEF)*TAU(5)
      RE(JEF)=RE(JEF) + B(JEF)*TAU(2) + B(IEF)*TAU(4) + B(KEF)*TAU(6)
      RE(KEF)=RE(KEF) + B(KEF)*TAU(3) + B(IEF)*TAU(5) + B(JEF)*TAU(6)
 550  CONTINUE
C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------
      IF (IFREF) 900,700,900
 700  IF (NLOPT.NE.0) GOTO 750
	FAC=SQRT(DVOL)
	DO 720 IEF=1,NEF
 720	B(IEF)=FAC*B(IEF)
	KS=0
	DO 730 IEF=1,NEF
	DO 730 JEF=IEF,NEF
	KS=KS+1
 730	S(KS)=S(KS)+B(IEF)*B(JEF)
	GOTO 790
C
 750  IF (MTMOD.LE.2) CALL SOK0NL (S,B,NEF,DVOL)
c      IF (MTMOD.GT.2) CALL SOK0PL (S,DP,B,NEF,DVOL)
	IF (MTMOD.EQ.3) CALL SOK0PL (S,DP,B,NEF,DVOL)

	IF (MTMOD.EQ.6) THEN
	CALL SOK0PL (S,DP,B,NEF,DVOL)
	END IF
C     -------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------
 790  IF (NLOPT.LE.1) GOTO 810
 800  CALL SOKSIG (S,TAU,B,NEF)
 810  CONTINUE !TIM(12) = TIM(12) + TIM2-TIM1
C
 900  CONTINUE

C	------------
C	MULTIPLY [D] BY THE INTEGRATED [BT]*[B]
C	------------
	IF (ITASK.GT.1) RETURN
	IF (IFREF.NE.0 .OR. NLOPT.NE.0) RETURN
	CALL SOK0EL(S,NEF,NNO)

      RETURN
      END


C============END OF COMPATIBLE 8-NODE SOLID ELEMENT===========





C=====================================================================
      SUBROUTINE SHAP3D (R,S,T,H,P,NODEX,NNO)
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
 9    H(9)    = .25*RR*SP*TP
      P(1,9)  =-.50*R*SP*TP
      P(2,9)  = .25*RR*TP
      P(3,9)  = .25*RR*SP
      GOTO 100
C
 10   H(10)   = .25*RM*SS*TP
      P(1,10) =-.25*SS*TP
      P(2,10) =-.50*RM*S*TP
      P(3,10) = .25*RM*SS
      GOTO 100
C
 11   H(11)   = .25*RR*SM*TP
      P(1,11) =-.50*R*SM*TP
      P(2,11) =-.25*RR*TP
      P(3,11) = .25*RR*SM
      GOTO 100
C
 12   H(12)   = .25*RP*SS*TP
      P(1,12) = .25*SS*TP
      P(2,12) =-.50*RP*S*TP
      P(3,12) = .25*RP*SS
      GOTO 100
C
 13   H(13)   = .25*RR*SP*TM
      P(1,13) =-.50*R*SP*TM
      P(2,13) = .25*RR*TM
      P(3,13) =-.25*RR*SP
      GOTO 100
C
 14   H(14)   = .25*RM*SS*TM
      P(1,14) =-.25*SS*TM
      P(2,14) =-.50*RM*S*TM
      P(3,14) =-.25*RM*SS
      GOTO 100
C
 15   H(15)   = .25*RR*SM*TM
      P(1,15) =-.50*R*SM*TM
      P(2,15) =-.25*RR*TM
      P(3,15) =-.25*RR*SM
      GOTO 100
C
 16   H(16)   = .25*RP*SS*TM
      P(1,16) = .25*SS*TM
      P(2,16) =-.50*RP*S*TM
      P(3,16) =-.25*RP*SS
      GOTO 100
C
 17   H(17)   = .25*RP*SP*TT
      P(1,17) = .25*SP*TT
      P(2,17) = .25*RP*TT
      P(3,17) =-.50*RP*SP*T
      GOTO 100
C
 18   H(18)   = .25*RM*SP*TT
      P(1,18) =-.25*SP*TT
      P(2,18) = .25*RM*TT
      P(3,18) =-.50*RM*SP*T
      GOTO 100
C
 19   H(19)   = .25*RM*SM*TT
      P(1,19) =-.25*SM*TT
      P(2,19) =-.25*RM*TT
      P(3,19) =-.50*RM*SM*T
      GOTO 100
C
 20   H(20)   = .25*RP*SM*TT
      P(1,20) = .25*SM*TT
      P(2,20) =-.25*RP*TT
      P(3,20) =-.50*RP*SM*T
      GOTO 100
C
 21   H(21)   = RR*SS*TT
      P(1,21) =-2.*R*SS*TT
      P(2,21) =-2.*S*RR*TT
      P(3,21) =-2.*T*RR*SS
 100  CONTINUE
C     ------------------------------------------------------
C     MODIFY FIRS 8 FUNCTIONS IF 9 OR MORE NODES ARE PRESENT
C     ------------------------------------------------------
      DO 290  IMI=1,NMI
      IN = NODEX(IMI)
      IF (IN.GT.16)  GOTO 210
      I1 = IN-8
      I2 = IPERM(I1)
      GOTO 250
 210  IF (IN.EQ.21)  GOTO 300
      I1 = IN-16
      I2 = I1+4
C
 250  H(I1) = H(I1) - .5*H(IN)
      H(I2) = H(I2) - .5*H(IN)
      H(IMI+8) = H(IN)
      DO 290  J=1,3
      P(J,I1) = P(J,I1) - .5*P(J,IN)
      P(J,I2) = P(J,I2) - .5*P(J,IN)
 290  P(J,IMI+8) = P(J,IN)
      RETURN
C     -----------------------------------------------
C     MODIFY FIRST 20 FUNCTIONS IF NODE 21 IS PRESENT
C     -----------------------------------------------
 300  DO 390  IMI=1,NMI
      IN = NODEX(IMI)
      IF (IN.GT.16)  GOTO 310
      I1 = IN-8
      I2 = IPERM(I1)
      GOTO 350
 310  IF (IN.EQ.21)  GOTO 400
      I1 = IN-16
      I2 = I1+4
C
 350  H(I1) = H(I1) + .125*H(21)
      H(I2) = H(I2) + .125*H(21)
      DO 390  J=1,3
      P(J,I1) = P(J,I1) + .125*P(J,21)
 390  P(J,I2) = P(J,I2) + .125*P(J,21)
C
 400  DO 410  I=1,8
      H(I) = H(I) - .125*H(21)
      DO 410  J=1,3
 410  P(J,I) = P(J,I) - .125*P(J,21)
      NN = NMI+7
      IF (NN.EQ.8)  RETURN
      DO 420  I=9,NN
      H(I) = H(I) - .25*H(21)
      DO 420  J=1,3
 420  P(J,I) = P(J,I) - .25*P(J,21)
      H(NMI+8) = H(21)
      DO 430  J=1,3
 430  P(J,NMI+8) = P(J,21)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE JACO3D (XY,P,XJ,XJI,DET,MEL,NNO)
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
      DIMENSION XY(3,NNO),P(3,21),XJ(3,3),AJ(9),XJI(9)
C
C     --------------------
C     JACOBIAN MATRIX (XJ)
C     --------------------
      DO 100  I=1,3
      DO 100  J=1,3
      DUM = 0.
      DO 90  INO=1,NNO
 90   DUM = DUM + P(I,INO)*XY(J,INO)
 100  XJ(I,J) = DUM
      K=1
      DO 110 J=1,3
      DO 110 I=1,3
	AJ(K)=XJ(I,J)
	K=K+1
 110  CONTINUE
C     ----------------------------------------
C     DETERMINANT OF THE JACOBIAN MATRIX (DET)
C     ----------------------------------------
      DET = AJ(1)*AJ(5)*AJ(9) + AJ(4)*AJ(8)*AJ(3) + AJ(7)*AJ(2)*AJ(6)
     1    - AJ(7)*AJ(5)*AJ(3) - AJ(4)*AJ(2)*AJ(9) - AJ(1)*AJ(8)*AJ(6)

C	IF (DET.LT.0.0D0) THEN
C	WRITE(110,*) 'NNO',NNO,MEL
C	WRITE(110,*) 'DET',DET
C	DO INO = 1,NNO 
C	WRITE(110,*) INO,XY(1:3,INO)
C	ENDDO
C	STOP
C	ENDIF

      IF (DET.LE.0.0D0) CALL ERRORS (15,H,MEL,'JACOB.DET.')

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
      SUBROUTINE SOBMAT (P,XJI,B,NNO)
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
      DIMENSION P(3,21),XJI(3,3),B(63)
C
      CALL CLEARA (B,63)
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
      RETURN
      END
C
  

C===========================================================
C	END OF ENHACED 8-NODE SOLID ELEMNT
C===========================================================






C=========================================================
C	COMMON SUBROUTINES FOR SOLID ELEMENT
C=========================================================
C=====================================================================
      SUBROUTINE SOLSIG (STRAIN,STRESS)
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
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION STRAIN(6),STRESS(6)
C
      STRESS(1) = A1*STRAIN(1) + B1*(STRAIN(2)+STRAIN(3))
      STRESS(2) = A1*STRAIN(2) + B1*(STRAIN(1)+STRAIN(3))
      STRESS(3) = A1*STRAIN(3) + B1*(STRAIN(1)+STRAIN(2))
      DO 100  I=4,6
 100  STRESS(I) = C1*STRAIN(I)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SOK0EL (S,NEF,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     CALCULATES ELEMENT STIFFNESS MATRIX FOR THE 3-D SOLID
C	-----------------------------------------------------
C     VARIABLES IN ARGUMENT LIST
C	--------------------------
C     S(NWS)    = ELEMENT STIFFNESS MATRIX (UPPER TRIANGULAR ROW-WISE)
C     NEF       = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     NNO       = NUMBER OF NODES USED TO DESCRIBE ELEMENT
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION S(1),BB(9)
C
	KL = 1
      DO 900  INO=1,NNO
      I0 = 3*(INO-1)
      DO 800  JNO=INO,NNO
C     -----------------------------------
C     PRESERVE VALUES OF [BT]*[B] IN [BB]
C     -----------------------------------
      KS = KL
      KB = 0
      DO 190  I=1,3
      DO 150  J=1,3
      KB = KB+1
      BB(KB) = S(KS)
 150  KS = KS+1
 190  KS = KS+NEF-I0-I-3
C     ---------------------------------------
C     MULTIPLY [D] BY THE INTEGRATED [BT]*[B]
C     ---------------------------------------
      KS1 = KL
      KS2 = KS1+NEF-I0-1
      KS3 = KS2+NEF-I0-2
      S(KS1)   = BB(1)*A1 + (BB(5)+BB(9))*C1
      S(KS2+1) = BB(5)*A1 + (BB(1)+BB(9))*C1
      S(KS3+2) = BB(9)*A1 + (BB(1)+BB(5))*C1
      IF (INO.EQ.JNO)  GOTO 300
C
      S(KS1+1) = BB(2)*B1 + BB(4)*C1
      S(KS1+2) = BB(3)*B1 + BB(7)*C1
      S(KS2)   = BB(4)*B1 + BB(2)*C1
      S(KS2+2) = BB(6)*B1 + BB(8)*C1
      S(KS3)   = BB(7)*B1 + BB(3)*C1
      S(KS3+1) = BB(8)*B1 + BB(6)*C1
      GOTO 800
C
 300  S(KS1+1) = BB(2)*(B1+C1)
      S(KS1+2) = BB(3)*(B1+C1)
      S(KS2+2) = BB(6)*(B1+C1)
C
 800  KL = KL+3
 900  KL = KL + 2*(NEF-I0) - 3
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SOK0NL (S,B,NEF,DVOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ADDS LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX
C	----------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     S(NWS)    = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     B(3*NNO)  = COMPACTED STRAIN-DISPLACEMENT MATRIX
C     NEF       = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     DVOL      = INTEGRATION CONSTANT =WT*DET
C     ----------------------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION S(1),B(63)
C
      D11 = A1*DVOL
      D12 = B1*DVOL
      D44 = C1*DVOL
      KS  = 1
C
      DO 190  JEF=1,NEF,3
      DB1 = D11*B(JEF)
      DB2 = D12*B(JEF)
      DB3 = DB2
      DB4 = D44*B(JEF+1)
      DB5 = D44*B(JEF+2)
      DO 150  IEF=JEF,NEF,3
      S(KS)   = S(KS)   + B(IEF)  *DB1 + B(IEF+1)*DB4 + B(IEF+2)*DB5
      S(KS+1) = S(KS+1) + B(IEF+1)*DB2 + B(IEF)  *DB4
      S(KS+2) = S(KS+2) + B(IEF+2)*DB3 + B(IEF)  *DB5
 150  KS = KS+3
 190  KS = KS + 2*(NEF-JEF) - 1
C
      KS = NEF+1
      DO 290  JEF=2,NEF,3
      DB1 = D12*B(JEF)
      DB2 = D11*B(JEF)
      DB3 = DB1
      DB4 = D44*B(JEF-1)
      DB6 = D44*B(JEF+1)
      S(KS)   = S(KS)   + B(JEF)  *DB2 + B(JEF-1)*DB4 + B(JEF+1)*DB6
      S(KS+1) = S(KS+1) + B(JEF+1)*DB3 + B(JEF)  *DB6
      KS   = KS+2
      JEF2 = JEF+2
      IF (NEF-JEF2)  290,210,210
 210  DO 250  IEF=JEF2,NEF,3
      S(KS)   = S(KS)   + B(IEF)  *DB1 + B(IEF+1)*DB4
      S(KS+1) = S(KS+1) + B(IEF+1)*DB2 + B(IEF)  *DB4 + B(IEF+2)*DB6
      S(KS+2) = S(KS+2) + B(IEF+2)*DB3 + B(IEF+1)*DB6
 250  KS = KS+3
 290  KS = KS + 2*(NEF-JEF) - 1
C
      KS = 2*NEF
      DO 390  JEF=3,NEF,3
      DB1 = D12*B(JEF)
      DB2 = DB1
      DB3 = D11*B(JEF)
      DB5 = D44*B(JEF-2)
      DB6 = D44*B(JEF-1)
      S(KS) = S(KS) + B(JEF)*DB3 + B(JEF-2)*DB5 + B(JEF-1)*DB6
      KS   = KS+1
      JEF1 = JEF+1
      IF (NEF-JEF1)  390,310,310
 310  DO 350  IEF=JEF1,NEF,3
      S(KS)   = S(KS)   + B(IEF)  *DB1 + B(IEF+2)*DB5
      S(KS+1) = S(KS+1) + B(IEF+1)*DB2 + B(IEF+2)*DB6
      S(KS+2) = S(KS+2) + B(IEF+2)*DB3 + B(IEF)  *DB5 + B(IEF+1)*DB6
 350  KS = KS+3
 390  KS = KS + 2*(NEF-JEF) - 1
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SOK0PL (S,DP,B,NEF,DVOL)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ADDS LINEAR CONTRIBUTION TO STIFFNESS MATRIX OF THE 3-D ELASTO-
C     PLASTIC SOLID ELEMENT
C	---------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     S(NWS)    = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     DP(6,6)   = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     B(3,NEF)  = COMPACTED STRAIN-DISPLACEMENT MATRIX
C     NEF       = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     DVOL      = INTEGRATION CONSTANT (WT*DET)
C     ----------------------------------------------------------------
      DIMENSION S(1),DP(6,6),B(1),DB(6)
C
      DO 90  I=1,6
      DO 90  J=1,I
      DP(I,J) = DP(I,J)*DVOL
 90   DP(J,I) = DP(I,J)
C
      KS = 1
      DO 190  J=1,NEF,3
      DO 110  K=1,6
 110  DB(K) = DP(K,1)*B(J) + DP(K,4)*B(J+1) + DP(K,5)*B(J+2)
      DO 150  I=J,NEF,3
      S(KS)   = S(KS)   + B(I)  *DB(1) + B(I+1)*DB(4) + B(I+2)*DB(5)
      S(KS+1) = S(KS+1) + B(I+1)*DB(2) + B(I)  *DB(4) + B(I+2)*DB(6)
      S(KS+2) = S(KS+2) + B(I+2)*DB(3) + B(I)  *DB(5) + B(I+1)*DB(6)
 150  KS = KS + 3
 190  KS = KS + 2*(NEF-J) - 1
C
      KS = NEF+1
      DO 290  J=2,NEF,3
      DO 210  K=1,6
 210  DB(K) = DP(K,2)*B(J) + DP(K,4)*B(J-1) + DP(K,6)*B(J+1)
      S(KS)   = S(KS)   + B(J)  *DB(2) + B(J-1)*DB(4) + B(J+1)*DB(6)
      S(KS+1) = S(KS+1) + B(J+1)*DB(3) + B(J-1)*DB(5) + B(J)  *DB(6)
      K  = J+2
      KS = KS+2
      IF (NEF-K)  300,250,250
 250  DO 280  I=K,NEF,3
      S(KS)   = S(KS)   + B(I)  *DB(1) + B(I+1)*DB(4) + B(I+2)*DB(5)
      S(KS+1) = S(KS+1) + B(I+1)*DB(2) + B(I)  *DB(4) + B(I+2)*DB(6)
      S(KS+2) = S(KS+2) + B(I+2)*DB(3) + B(I)  *DB(5) + B(I+1)*DB(6)
 280  KS = KS+3
 290  KS = KS + 2*(NEF-J) - 1
C
 300  KS = 2*NEF
      DO 390  J=3,NEF,3
      DO 310  K=1,6
 310  DB(K) = DP(K,3)*B(J) + DP(K,5)*B(J-2) + DP(K,6)*B(J-1)
      S(KS) = S(KS) + B(J)*DB(3) + B(J-2)*DB(5) + B(J-1)*DB(6)
      K  = J+1
      KS = KS+1
      IF (NEF-K)  400,350,350
 350  DO 380  I=K,NEF,3
      S(KS)   = S(KS)   + B(I)  *DB(1) + B(I+1)*DB(4) + B(I+2)*DB(5)
      S(KS+1) = S(KS+1) + B(I+1)*DB(2) + B(I)  *DB(4) + B(I+2)*DB(6)
      S(KS+2) = S(KS+2) + B(I+2)*DB(3) + B(I)  *DB(5) + B(I+1)*DB(6)
 380  KS = KS + 3
 390  KS = KS + 2*(NEF-J) - 1
C
 400  RETURN
      END
C
C=====================================================================
      SUBROUTINE SOKSIG (S,SIG,B,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ADDS INITIAL STRESS CONTRIBUTION TO ELEMENT STIFFNESS MATRIX
C	------------------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     S(NWS)    = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     SIG(6)    = GAUSS POINT STRESSES (MULTIPLIED BY DVOL=WT*DET)
C     B(3*NNO)  = COMPACTED STRAIN DISPLACEMENT MATRIX
C     NEF       = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ----------------------------------------------------------------
      DIMENSION S(1),SIG(6),B(63)
C
      KL = 1
      DO 200  JEF=1,NEF,3
      SB1 = SIG(1)*B(JEF) + SIG(4)*B(JEF+1) + SIG(5)*B(JEF+2)
      SB2 = SIG(4)*B(JEF) + SIG(2)*B(JEF+1) + SIG(6)*B(JEF+2)
      SB3 = SIG(5)*B(JEF) + SIG(6)*B(JEF+1) + SIG(3)*B(JEF+2)
C
      KS1 = KL
      DO 100  IEF=JEF,NEF,3
      KS2 = KS1+NEF-JEF+1
      KS3 = KS2+NEF-JEF
      DUM = B(IEF)*SB1 + B(IEF+1)*SB2 + B(IEF+2)*SB3
      S(KS1) = S(KS1) + DUM
      S(KS2) = S(KS2) + DUM
      S(KS3) = S(KS3) + DUM
 100  KS1 = KS1 + 3
C
 200  KL  = KL + 3*(NEF-JEF)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE SOMASS (CM,H,DEN,DVOL,NNO,NEF,IPT)
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
      COMMON /DYNA/  IMASS
      DIMENSION CM(1),H(21)
C     --------------
C     INITIALISATION
C     --------------
      IF (IPT.GT.1) GOTO 50
      VOL = 0.
 50   DVOL = DVOL*DEN
      IF (IMASS-2) 100,100,300
C     ----------------------------------------
C     INTEGRATE LUMPED OR DIAGONAL MASS MATRIX
C     ----------------------------------------
 100  VOL = VOL + DVOL
      IF (IMASS.EQ.1) GOTO 500
C
      DO 250  INO=1,NNO
      KM = 3*(INO-1) + 1
 250  CM(KM) = CM(KM) + H(INO)*H(INO)
      GOTO 500
C     --------------------------------
C     INTEGRATE CONSISTENT MASS MATRIX
C     --------------------------------
 300  KM = 1
      DO 390  INO=1,NNO
      DO 350  JNO=INO,NNO
      CM(KM) = CM(KM) + H(INO)*H(JNO)*DVOL
 350  KM = KM+3
 390  KM = KM + 6*(NNO-INO) + 3
C     --------------------------------------------
C     FILL IN REMAINING TERMS IN Y AND Z DIRECTION
C     --------------------------------------------
 500  IF (IPT.LT.27) RETURN
      IF (IMASS-2) 600,700,800
C     ------------------------------
C     LUMPED OR DIAGONAL MASS MATRIX
C     ------------------------------
 600  FACM = VOL/NNO
      DO 650  IEF=1,NEF
 650  CM(IEF) = FACM
      RETURN
C
 700  SUM = 0.
      DO 710  IEF=1,NEF,3
 710  SUM = SUM + CM(IEF)
      FACM = VOL/SUM
      DO 750  IEF=1,NEF,3
      CM(IEF) = CM(IEF)*FACM
      CM(IEF+1) = CM(IEF)
 750  CM(IEF+2) = CM(IEF)
      RETURN
C     ----------------------
C     CONSISTENT MASS MATRIX
C     ----------------------
 800  KM1 = 1
      DO 890  IEF=1,NEF,3
      KM2 = KM1 + NEF-IEF+1
      KM3 = KM2 + NEF-IEF
      DO 850  JEF=IEF,NEF,3
      CM(KM2) = CM(KM1)
      CM(KM3) = CM(KM1)
      KM1 = KM1+3
      KM2 = KM2+3
 850  KM3 = KM3+3
 890  KM1 = KM1 + 2*(NEF-IEF) - 1
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE MISE3D (SIGP,EPSP,YIELD,IPEL,EPS,SIG,DP)
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
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
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
 110  DELEPS(I) = EPS(I) - EPSP(I)
      CALL SOLSIG (DELEPS,DELSIG)
      DO 130  I=1,6
 130  SIG(I) = SIGP(I) + DELSIG(I)
C     --------------------------------------------------
C     4. CALCULATE DEVIATORIC STRESSES AND CHECK WHETHER
C        STATE OF STRESS FALLS OUTSIDE THE YIELD SURFACE
C     --------------------------------------------------
      CALL DEVI3D (SIG,DF)
      F = 1.5*(SX*SX + SY*SY + SZ*SZ)+3.*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
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
      DO 460  I=1,6
 460  SIG(I) = SIGP(I)
C     -------------------------------------------------------
C     5. DETERMINE PART OF STRAIN TAKEN ELASTICLY AND ADD
C        STRESSES DUE TO THE ELASTIC STRAINS TO PREVIOUS STRESSES
C     -----------------------------------------------------------
      CALL DEVI3D (SIG,DF)
      CALL DEVI3D (DELSIG,DS)
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
 600  NINC = NSINC*(SQRT(F)-YLD)/YLD + 1
      IF (NINC.GT.2*NSINC)  NINC = 2*NSINC
      FACT = (1.-RATIO) / NINC
      DO 610  I=1,6
 610  DEPS(I) = FACT*DELEPS(I)
C     ----------------------------------------------------
C     7. CALCULATE PLASTIC STRAIN RATE MULTIPLIER AND
C        CONVERT ELASTIC PART OF STRAIN SUBINC. TO STRESSES
C     -----------------------------------------------------
      DO 890  IINC=1,NINC
      CALL DEVI3D (SIG,DF)
      PLAMDA = (SX*DEPS(1) + SY*DEPS(2) + SZ*DEPS(3) + SXY*DEPS(4) +
     1          SXZ*DEPS(5) + SYZ*DEPS(6)) * 0.5/YLD2
      IF (PLAMDA.LT.0.)  PLAMDA = 0.
      DO 710  I=1,6
 710  DEPSE(I) = DEPS(I) - 3.*PLAMDA*DF(I)
      CALL SOLSIG (DEPSE,DELSIG)
      DO 720  I=1,6
 720  SIG(I) = SIG(I) + DELSIG(I)
C     ----------------------------------------------
C     8. STRESS CORRECTION BACK TO THE YIELD SURFACE
C     ----------------------------------------------
      CALL DEVI3D (SIG,DF)
      F = 1.5*(SX*SX + SY*SY + SZ*SZ)+3.*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
C	ABCD=1./(DFLOAT(ITOLEY))
	ABCD=1./1000000.
      IF (F-YLD2.LE.ABCD)  GOTO 890
      COEF = -1. + SQRT(YLD2/F)
      DO 810  I=1,6
 810  SIG(I) = SIG(I) + COEF*DF(I)
C
 890  CONTINUE
C     ----------------------------------------------
C     9. UPDATE PREVIOUS STRESSES, STRAINS AND YIELD
C     ----------------------------------------------
 900  DO 910  I=1,6
      SIGP(I) = SIG(I)
 910  EPSP(I) = EPS(I)
C     ------------------------------------
C     FORM THE ELASTO-PLASTIC MATERIAL LAW
C     ------------------------------------
      CALL DELP3D (SIG,DEPS,DP,DP,IPEL)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DEVI3D (SIG,S)
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
 100  S(I+3) = SIG(I+3)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DELP3D (SIG,DEPS,DP,CP,IPEL)
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
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
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
      CALL DEVI3D (SIG,S)
      WP = 0.
      DO 100  I=1,6
 100  WP = WP + S(I)*DEPS(I)
      IF (WP.LT.0.)  GOTO 300
C     ---------------------------------------------
C     DEDUCT PLASTIC CONTRIBUTION IN LOWER TRIANGLE
C     ---------------------------------------------
      C2BETA = C2*1.5/YLD/YLD
      DO 200  I=1,6
      DO 200  J=1,I
 200  DP(I,J) = DP(I,J) - C2BETA*S(J)*S(I)
C     --------------------------------------------
C     FILL IN SYMMETRIC ELEMENTS IN UPPER TRIANGLE
C     --------------------------------------------
 300  DO 390  I=1,6
      DO 390  J=I,6
 390  DP(I,J) = DP(J,I)
C
      RETURN
      END
C
C=====================================================================
C	END 3D SOLID ELEMENT
C=====================================================================
C=======================================================================
      SUBROUTINE SOLID82(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,MWG
     +,ALPHA,SEL,SEDI,FIN)
c 200  CALL       SOLI63(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,ELOD,NWG
c     +,ALPHA,SEL,SEDI)
C
 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------------------
C     MAIN PROGRAM FOR THE 3-D SOLID-solidliu
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
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
C
      PARAMETER (MM=5)
C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(*),COORD(*)
      DIMENSION EDIS(24),EDISI(*),RE(*)
      DIMENSION DP(6,6),H(21),P(3,21),XJI(9),B(63),DISD(9)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
C     ------------
C     EAS METHOD
C     ------------    
	DIMENSION XJ(3,3),TTO(6,6),TT(6,6),XJO(3,3),TM(6,MM)
	DIMENSION EAS(6,1),ALPHA(*),TMT(MM,6),RE1(24),RE2(24,1)
      DIMENSION SED(MM,MM),SEL(MM,1),SEDI(*),RHH(MM,1),RH(MM,1)
c    ---------------------------
C           ADDED BY LIU
C    ---------------------------
      DIMENSION BB(6,24)
C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------
      CALL HOKLAW (PROPM,PROPG,1)
C     --------------
C     COMPUTE JACOBIAN MATRIX AT THE ORIGIN (R=0.0,S=0.0,T=0.0)
C     --------------
      CALL SHAP3D (0.0D0,0.0D0,0.0D0,H,P,NODEX,NNO)
	CALL JACO3D(COORD,P,XJO,XJI,DETO,MEL,NNO)
C     ----------------
C     COMPUTE T MATRIX (FOR EAS METHOD)
C     ---------------
       CALL MATRIXT(XJI,TT)
       
 	DO 5 I=1,6
 	DO 5 J=1,6
  	TTO(I,J)=TT(J,I)
  5   CONTINUE

      MGR = NGR
      MGS = NGS
      MGT = NGT
      IF (ITASK.NE.5) GOTO 10
      MGR = 3
      MGS = 3
      MGT = 3
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
 10   IPT = 0
C     ---------------
C     INITIALIZATION OF SOME VARIABLE MATRICES
C     ---------------
      DO 11 I=1,MM
	DO 12 J=1,MM
	SED(I,J)=0.0
 12   CONTINUE
      DO 13 K=1,24
	SEL(I,K)=0.0
 13   CONTINUE
      RH(I,1)=0.0
 11   CONTINUE
      DO 17 I=1,24
	RE1(I)=0.0
 17   CONTINUE

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
      CALL SHAP3D (RI,SI,TI,H,P,NODEX,NNO)
      CALL JACO3D (COORD,P,XJ,XJI,DET,MEL,NNO)
      DVOL = WT*DET
C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
      IF (ITASK.NE.5)  GOTO 50
      CALL SOMASS (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
      GOTO 900
C     ----------------------------------------
C     COMPACTED STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------
 50   CALL SOBMAT (P,XJI,B,NNO)
      WEIGHTL=DETO/DET
      CALL BMAT(B,XJI,XJ,BB,RI,SI,TTO,P,WEIGHTL) !added by Liu
C     ----------------     
C     FOR ENHANCED STRAIN METHOD 
C     COMPUTE MATRIX M AND ITS TRANSPOSE MATRIX
C     ----------------
      CALL MATRIXTMliu(RI,SI,TI,TTO,TM,MM)
      DO 51 I=1,6
  	DO 51 J=1,MM
  	TM(I,J)=TM(I,J)*DETO/DET
  	TMT(J,I)=TM(I,J)
 51   CONTINUE
C     ------------------
C     D=TRANSPOSE(M)*E*M
C     L=TRANSPOSE(M)*E*B
C     ------------------
      CALL MATRIXD(TM,SED,DVOL,MM)
      CALL MATRIXLIU(TM,SEL,BB,DVOL)
      IF (NLOPT+ITASK.EQ.1)  GOTO 700
C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------
      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,6
 100  TAU(I) = WA(I,IPT)*DVOL
      GOTO 800
C     -------------------------------
C     DISPLACEMENT DERIVATIVES (DISD)
C     -------------------------------
 200  IF (ITASK.GT.1) THEN
      CALL CLEARA (DISD,9)
      DO 290  IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDIS(IEF)
      DISD(2) = DISD(2) + B(JEF)*EDIS(JEF)
      DISD(3) = DISD(3) + B(KEF)*EDIS(KEF)
      DISD(4) = DISD(4) + B(JEF)*EDIS(IEF)
      DISD(5) = DISD(5) + B(KEF)*EDIS(IEF)
      DISD(6) = DISD(6) + B(IEF)*EDIS(JEF)
      DISD(7) = DISD(7) + B(KEF)*EDIS(JEF)
      DISD(8) = DISD(8) + B(IEF)*EDIS(KEF)
 290  DISD(9) = DISD(9) + B(JEF)*EDIS(KEF)
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
 !     STRAIN(1) = DISD(1)
 !     STRAIN(2) = DISD(2)
 !     STRAIN(3) = DISD(3)
 !     STRAIN(4) = DISD(4) + DISD(6)
  !    STRAIN(5) = DISD(5) + DISD(8)
  !    STRAIN(6) = DISD(7) + DISD(9)
      strain=matmul(bb,edis)
C     --------
C     FOR EAS TERM (EPS(EAS)=M*ALPHA)
C     --------
       CALL MATRIX2MUL(TM,ALPHA,EAS,6,MM,1)
 	DO 395 I=1,6
         STRAIN(I)=STRAIN(I)-EAS(I,1)
 395  CONTINUE
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
       DO 390  I=1,6
 390  STRAIN(I) = STRAIN(I) - QSTRAI(I)
C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------
 400  IF (MTMOD-3)  405,450,460
 405  CALL SOLSIG (STRAIN,STRESS)
      DO 410  I=1,6
 410  WA(I,IPT) = STRESS(I)
      GOTO 500
C
 450  CALL MISE3D (WA(1,IPT),WA(7,IPT),WA(13,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DP)
      GOTO 500
C
C460  CALL CONC3D -- DELETED BECAUSE OF SUBROUTINE NOT EXIST
 460  CONTINUE
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510  I=1,6
      TAU(I) = STRESS(I)*DVOL
 510  CONTINUE
      ENDIF
C      IF (ITASK.LE.2) GOTO 520
      IF (ITASK.LE.3) GOTO 520
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900
 520  IF (ITASK.GT.1) THEN
C     ---------
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     ---------
!	DO 550  IEF=1,NEF,3
  !    JEF = IEF+1
  !    KEF = IEF+2
  !    RE1(IEF)=RE1(IEF) + B(IEF)*TAU(1) + B(JEF)*TAU(4) + B(KEF)*TAU(5)
   !   RE1(JEF)=RE1(JEF) + B(JEF)*TAU(2) + B(IEF)*TAU(4) + B(KEF)*TAU(6)
   !   RE1(KEF)=RE1(KEF) + B(KEF)*TAU(3) + B(IEF)*TAU(5) + B(JEF)*TAU(6)
      re1=matmul(transpose(bb),tau)+re1
 550  CONTINUE
C     ---------
C     COMPUTE EQUILIBRIUM FORCE FOR EAS METHOD
C     ---------
        CALL MATRIX2MUL(TMT,TAU,RHH,MM,6,1)
 	DO 552 I=1,MM
  	RH(I,1)=RH(I,1)+RHH(I,1)
 552  CONTINUE
      ENDIF
C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------
      IF (IFREF) 900,700,900
 700  IF (NLOPT.NE.0) GOTO 750
      CALL SOK0NLIU(S,BB,DVOL)
!	CALL SOK0NL (S,B,NEF,DVOL)
	GOTO 790
C
 750  IF (MTMOD.LE.2) CALL SOK0NLIU(S,BB,DVOL)
      IF (MTMOD.GT.2) CALL SOK0PLIU (S,DP,BB,NEF,DVOL)
C     -------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------
 790  IF (NLOPT.LE.1) GOTO 810
 800  CALL SOKSIG (S,TAU,B,NEF)
 810  CONTINUE !TIM(12) = TIM(12) + TIM2-TIM1
      
C
 900  CONTINUE
C     --------
C     FOR NONLINEAR ANALYSIS OF EAS METHOD
C     --------
     	CALL INVMAT1(SED,SEDI,MM)
        CALL MATRIX3MUL(SEL,SEDI,RH,RE2,24,MM,1)
        IF (ITASK.GT.1) THEN
 	DO 910 I=1,24
   	RE(I)=RE1(I)-RE2(I,1)
 910  CONTINUE
         ENDIF

C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF

C     ------------------
C     FOR ENHANCED STRAIN METHOD
C     K=TRANSPOSE(L)*INVERSE(D)*L
C     K(FINAL)=K(PURE DISP.)-K(EAS)
C     ------------------ 
       IF (ITASK.GT.1) RETURN
        CALL MATRIXES(SEL,SEDI,S,MM)
      RETURN
      END
C
C===========================================================
      SUBROUTINE MATRIX3MUL(A,B,C,D,N1,N2,N3)
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
C
C=============================================================
	SUBROUTINE INVMAT1(A,B,N)
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
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N,N),B(N,N),C(N,N)

C     CONSTRUCT IDENTITY MATRIX B(I,J)=I
      EPS=1.0E-12
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
   33 STOP
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
C

C
C	==============================================================
C	SUBROUTINES FOR 3D SOLID ELEMENT
C     SUBROUTINES FOR ENHANCED ASSUMED STRAIN (EAS) METHOD
C	==============================================================
      SUBROUTINE MATRIXD(TM,DD,DVOL,MM)
C     ---------
C     TO COMPUTE DD=MT*E*M
C     ---------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
	DIMENSION TM(6,MM),DD(MM,MM),DE(MM,6)
	DO 10 I=1,MM
	DE(I,1)=TM(1,I)*A1+TM(2,I)*B1+TM(3,I)*B1
	DE(I,2)=TM(1,I)*B1+TM(2,I)*A1+TM(3,I)*B1
	DE(I,3)=TM(1,I)*B1+TM(2,I)*B1+TM(3,I)*A1
	DE(I,4)=TM(4,I)*C1
	DE(I,5)=TM(5,I)*C1
	DE(I,6)=TM(6,I)*C1
 10   CONTINUE
	DO 20 I=1,MM
	DO 20 J=1,MM
	DO 20 K=1,6
      DD(I,J)=DD(I,J)+DE(I,K)*TM(K,J)*DVOL
  20  CONTINUE
      RETURN
	END

C	==============================================================
      SUBROUTINE MATRIXL(TM,DD,B,DVOL,MM)
C     ---------
C     TO COMPUTE L=MT*E*B
C     ---------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
	DIMENSION TM(6,MM),DD(MM,24),DE(MM,6),B(63)
	DO 10 I=1,MM
	DE(I,1)=TM(1,I)*A1+TM(2,I)*B1+TM(3,I)*B1
	DE(I,2)=TM(1,I)*B1+TM(2,I)*A1+TM(3,I)*B1
	DE(I,3)=TM(1,I)*B1+TM(2,I)*B1+TM(3,I)*A1
	DE(I,4)=TM(4,I)*C1
	DE(I,5)=TM(5,I)*C1
	DE(I,6)=TM(6,I)*C1
 10   CONTINUE
      DO 20 I=1,MM
	DO 30 K1=1,24,3
	K2=K1+1
	K3=K1+2
	DD(I,K1)=DD(I,K1)+(DE(I,1)*B(K1)+DE(I,4)*B(K2)+DE(I,5)*B(K3))*DVOL
	DD(I,K2)=DD(I,K2)+(DE(I,2)*B(K2)+DE(I,4)*B(K1)+DE(I,6)*B(K3))*DVOL
	DD(I,K3)=DD(I,K3)+(DE(I,3)*B(K3)+DE(I,5)*B(K1)+DE(I,6)*B(K2))*DVOL
 30   CONTINUE
 20   CONTINUE     
	RETURN
	END

C	==============================================================
      SUBROUTINE MATRIXES(A,B,S,MM)
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
 10   CONTINUE
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
 50   CONTINUE
      KK=KK+1
 40   CONTINUE
      RETURN
      END
	

C	==============================================================
      SUBROUTINE MATRIX2MUL(A,B,C,N1,N2,N3)
C	--------     
C     COMPUTE THE PRODUCT OF TWO MATRICES
C     COMPUTE A(N1,N2)*B(N2,N3)=C(N1,N3)
C	--------     
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION A(N1,N2),B(N2,N3),C(N1,N3)
	DO 5 I=1,N1
	DO 5 J=1,N3
	C(I,J)=0.0
  5   CONTINUE
	DO 10 I=1,N1
	DO 10 J=1,N3
	DO 10 K=1,N2
      C(I,J)=C(I,J)+A(I,K)*B(K,J)
  10  CONTINUE
      RETURN
      END

C	==============================================================
	
C	==============================================================


C	==============================================================
      SUBROUTINE MATRIXT(XJO,TTO)
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
 5    CONTINUE
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
      SUBROUTINE LOADVE3 (PROPG,IGSET,XYZ,NODEX,LM,RT,RP,R,
     1                   MGP,MXY,MEX,MEF,IND)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     CONVERTES UNIFORMLY DISTRIBUTED GLOBAL TRACTIONS OR NORMAL
C     PRESSURES TO EQUIVALENT NODAL LOADS
C	FOR 3D SOLID ELEMENT (NNO=8)
C	-----------------------------------
C     PROPG(NGP,NGPS)  = GEOMETRIC PROPERTIES
C     IGSET(NELE)      = GEOMETRIC SET NUMBER
C     XYZ(MXY,NELE)    = NODAL COORDINATES FOR ELEMENTS
C     NODEX(NEX,NELE)  = LOCATIONS OF EXCESSIVE NODES
C     LM(NEF,NELE)     = EQUATION NUMBERS FOR ELEMENT D.O.F.
C     RT(3,NELE)       = GLOBAL TRACTIONS IN X,Y,Z
C     RP(NELE)         = NORMAL PRESSURES
C     R(NEQ)           = LOAD VECTOR (EQUIVALENT NODAL LOADS)
C     IND              = FLAG FOR TRACTIONS (IND=2),PRESSURES (IND=3)
C     ---------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/ GLOC(10,10),GWT(10,10),NGR,NGS,NGT
	COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK
      COMMON A(9000000),IA(9000000)

C	ELEMENT LINK POSITION AND DOF SONGSAK MAR2007
	COMMON /EFLINK/ NFLINK(30,30)
C
      DIMENSION PROPG(MGP,1),IGSET(1),XYZ(MXY,1),
     +          NODEX(MEX,1),LM(MEF,1)
      DIMENSION RT(3,1),RP(1),R(1),H(8),P(3,8),XJ(3,3),XJI(9)
      DIMENSION VR(3),VS(3),VT(3),RL(3,8),IGPOS(9)
C	SUNIL 18/01/01 NEXT LINE REMOVED
C	EQUIVALENCE (ANO,NNO)

      LINK = NFLINK(ITYPE,ISTYP)
      CALL EDOFLINK(LINK,NNF,IGPOS)
C     -----------------------------------
C     ELEMENT LOOP,SKIP UNLOADED ELEMENTS
C     -----------------------------------
      DO 900  IELE=1,NELE
      INDEL = 0
      IF (IND.NE.2) GOTO 110
      DO 100  I=1,3
 100  IF (RT(I,IELE).NE.0.) INDEL = 1
      GOTO 120
 110  IF (RP(IELE).NE.0.) INDEL = 1
 120  IF (INDEL.EQ.0) GOTO 900
C
      CALL CLEARA (RL,24)
	MGR=NGR
	MGS=NGS
	MGT=NGT
	IF (ITASK.NE.5) GOTO 121
	MGR=3
	MGS=3
	MGT=3
C     ----------------
C     GAUSS POINT LOOP
C     ----------------
 121  DO 800  IGR=1,MGR
      RI = GLOC(IGR,MGR)
      DO 800  IGS=1,MGS
      SI = GLOC(IGS,MGS)
	DO 800  IGT=1,MGT
	TI = GLOC(IGT,MGT)
      WT = GWT(IGR,MGR)*GWT(IGS,MGS)*GWT(IGT,MGT)
C     ----------------------------------------------
C     SHAPE FUNCTIONS (H),JACOBIAN DETERMINANT (DET)
C     AND DIRECTION COSINES (VR,VS,VT)
C     ----------------------------------------------
      CALL SHAP3D (RI,SI,TI,H,P,NODEX(1,IELE),NNO)
      CALL JACO3D (XYZ(1,IELE),P,XJ,XJI,DET,MEL,NNO)
C     --------------------------------------------------------
C     UNIFORMLY DISTRIBUTED GLOBAL TRACTION (RT); CONSERVATIVE
C     --------------------------------------------------------
 300  IF (IND.NE.2) GOTO 400
      DO 390  INO=1,NNO
      FAC = H(INO)*WT*DET
      DO 380  I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 390
 380  RL(I,INO) = RL(I,INO) + RT(ICO,IELE)*FAC
 390  CONTINUE
      GOTO 800
C     ------------------------------------------------------------
C     UNIFORMLY DISTRIBUTED NORMAL PRESSURE (RP); NON-CONSERVATIVE
C     ------------------------------------------------------------
 400  DO 500  INO=1,NNO
      FAC = RP(IELE)*H(INO)*WT*DET
      DO 480  I=1,3
      ICO = IGPOS(I)
      IF (ICO.GT.3) GOTO 500
 480  RL(I,INO) = RL(I,INO) + RP(ICO)*FAC
 500  CONTINUE
C
 800  CONTINUE
C     ------------------------------------------------------
C     TRANSFORM INTO LOCAL COORDINATES AT SKEW NODES, IF ANY
C     ------------------------------------------------------
      IF (NLS.EQ.0) GOTO 700
      CALL LOCRES (IA(LID),IA(LDS),A(LDC),LM(1,IELE),A(LES),A(LED),
     1            A(LEI),RL,NSF,NNF,4)
C
 700  DO 600  INO=1,NNO
      II = (INO-1)*NNF
      DO 600  INF=1,3
      IEQ = LM(II+INF,IELE)
      IF (IEQ.NE.0) R(IEQ) = R(IEQ) + RL(INF,INO)
 600  CONTINUE
 900  CONTINUE
C
      RETURN
      END

C	========================================================
C	SUBROUTINE TO CHANGE THE SEQUENCE OF NODE NUMBERING
C	AFTER CHANGE, THE NODE NUMBERING OF 2D (4 NODES) IS
C	 1(+,+),2(-,+),3(-,-),4(+,-)
C	FOR SURFACE PRESSURE LOAD 
C	AND ONLY FOR 3D SOLID ELEMENT (8-NODE)
C	--------------------------
	SUBROUTINE CHANGE(A,F,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	DIMENSION A(2*NNO,1),B(2*NNO,1)
	INTEGER F(4),N(4)
	X1=A(1,1)
	Y1=A(2,1)

	X2=X1
	Y2=Y1
	
	X3=X1
	Y3=Y1
	
	X4=X1
	Y4=Y1
	
	DO 1 I=1,2*NNO,2
	IF (A(I,1).GE.X1 .AND. A(I+1,1).GE.Y1) THEN
	B(1,1)=A(I,1)
	B(2,1)=A(I+1,1)
	X1=A(I,1)
	Y1=A(I+1,1)
	N(1)=F((I+1)/2)
	ENDIF
	IF (A(I,1).LE.X2 .AND. A(I+1,1).GE.Y2) THEN
	B(3,1)=A(I,1)
	B(4,1)=A(I+1,1)
	X2=A(I,1)
	Y2=A(I+1,1)
	N(2)=F((I+1)/2)
	ENDIF
	IF (A(I,1).LE.X3 .AND. A(I+1,1).LE.Y3) THEN
	B(5,1)=A(I,1)
	B(6,1)=A(I+1,1)
	X3=A(I,1)
	Y3=A(I+1,1)
	N(3)=F((I+1)/2)
	ENDIF
	IF (A(I,1).GE.X4 .AND. A(I+1,1).LE.Y4) THEN
	B(7,1)=A(I,1)
	B(8,1)=A(I+1,1)
	X4=A(I,1)
	Y4=A(I+1,1)
	N(4)=F((I+1)/2)
	ENDIF
 1    CONTINUE

	DO 2 I=1,NNO*2
	A(I,1)=B(I,1)
 2    CONTINUE
	
	DO 3 I=1,NNO
	F(I)=N(I)
 3    CONTINUE
 	
	RETURN
	END

C==================================================	                
	SUBROUTINE MATRIXLIU(TM,DD,B,DVOL)
C     ---------
C     TO COMPUTE L=MT*E*B
C     ---------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
	DIMENSION TM(6,5),DD(5,24),DE(5,6),B(6,24)
	DIMENSION D(6,6)
C
      D=0.
	D(1,1)=A2
	D(1,2)=B2
	D(1,3)=B2
      D(2,2)=A2
	D(2,3)=B2
	D(2,1)=B2
	D(3,1)=B2
	D(3,2)=B2
	D(3,3)=A2
	D(4,4)=0.5
	D(5,5)=0.5
	D(6,6)=0.5
	D=C2*D
C
      DD=MATMUL(MATMUL(TRANSPOSE(TM),D),B)*dvol+DD
C
      RETURN
	END
c
c

c
c
      SUBROUTINE BMAT(B,XJI,XJ,BB,RI,SI,TTO,P,WEIGHTL)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------
C        SUBROUTINE TO COMPUTE THE ASSUMED STRAIN
C        BMATRIX IN FULL FORM (6,24)
C     -------------------------------------------------
      DIMENSION B(24),BB(6,24),XJI(3,3),XJ(3,3)
	DIMENSION PTA(8),PSA(8),PTB(8),PRB(8)
	DIMENSION PTC(8),PSC(8),PTD(8),PRD(8)
	DIMENSION PTA1(8),PTA2(8),PTA3(8),PTA4(8)
	DIMENSION RL(4),SL(4)
	DIMENSION TT(6),TTO(6,6),P(3,21),XJJ(24,24)
C
      CALL SHDSAMP(PSA,PTA,PRB,PTB,PSC,PTC,PRD,PTD,
     @                   PTA1,PTA2,PTA3,PTA4,RL,SL)
C
C     ------------------------------------------------
C      MODIFY TTO
C     ------------------------------------------------
      TTO(1,5)=0.
  	TTO(1,6)=0.
  	TTO(2,5)=0.
  	TTO(2,6)=0.
 	TTO(3,5)=0.
	TTO(3,6)=0.

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
C     ----DEFINE THE TRANSFORM MATRIX TT
 !     TT(1)=XJI(1,1)*XJI(3,3)+XJI(3,1)*XJI(1,3)
 ! 	TT(2)=XJI(1,2)*XJI(3,3)+XJI(3,2)*XJI(1,3)
!	TT(3)=XJI(2,1)*XJI(3,3)+XJI(3,1)*XJI(2,3)
!	TT(4)=XJI(2,2)*XJI(3,3)+XJI(3,2)*XJI(2,3)
!	TT(5)=XJI(1,1)*XJI(2,3)+XJI(2,1)*XJI(1,3)
!	TT(6)=XJI(1,2)*XJI(2,3)+XJI(2,2)*XJI(1,3) 
C
      BB=0.
	DO I=1,8
      BB(1,3*I-2)=P(1,I)
	BB(2,3*I-1)=P(2,I)
!	BB(3,3*I)=B(3*I)
      BB(4,3*I-2)=P(2,I)
	BB(4,3*I-1)=P(1,I)
!	BB(5,3*I-2)=P(3,I)
!	BB(5,3*I)=P(1,I)
!	BB(6,3*I-1)=P(3,I)
!	BB(6,3*I)=P(2,I)
	ENDDO
	DO I=1,8
C     ----MODIFY BB DUE TO ASSUMED SHEAR STRAIN--------
      BDT=(PTB(I)*(1-SI)+(1+SI)*PTD(I))*.5
      BDR=(PRB(I)*(1-SI)+(1+SI)*PRD(I))*.5
      ACT=(PTA(I)*(1-RI)+(1+RI)*PTC(I))*.5
      ACS=(PSA(I)*(1-RI)+(1+RI)*PSC(I))*.5
      BB(5,3*I-2) =BDT 

	BB(5,3*I)   = BDR
 
	BB(6,3*I-1) = ACT
      BB(6,3*I)   = ACS
!	BB(4,3*I-2)=BB(4,3*I-2)+TT(5)*(BDT*XJ(1,1)+BDR*XJ(3,1))
  !   @	+TT(6)*(ACT*XJ(2,1)+ACS*XJ(3,1))
!	BB(4,3*I-1)=BB(4,3*I-1)+TT(6)*(ACT*XJ(2,2)+ACS*XJ(3,2))+
 !    @     TT(5)*(BDT*XJ(1,2)+BDR*XJ(3,2))
	
C     ----MODIFY BB DUE TO ASSUMED THICKNESS STRAIN----     
	TH=((1+RL(1)*RI)*(1+SL(1)*SI)*PTA1(I)
     @	+(1+RL(2)*RI)*(1+SL(2)*SI)*PTA2(I)
     @	+(1+RL(3)*RI)*(1+SL(3)*SI)*PTA3(I)
     @	+(1+RL(4)*RI)*(1+SL(4)*SI)*PTA4(I))*0.25
C			    
      BB(3,3*I)=TH
!	BB(1,3*I-2)=BB(1,3*I-2)+XJI(1,3)*XJI(1,3)*TH*XJ(3,1)
!	BB(2,3*I-1)=BB(2,3*I-1)+XJI(2,3)*XJI(2,3)*TH*XJ(3,2)
      ENDDO
C
      BB=MATMUL(MATMUL(TTO,BB),XJJ)*WEIGHTL
	RETURN
	END
c
c

c=======================================================
      SUBROUTINE SOK0PLIU(S,DP,BB,NEF,DVOL)
	IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     ADDS LINEAR CONTRIBUTION TO STIFFNESS MATRIX OF THE 3-D ELASTO-
C     PLASTIC SOLID ELEMENT
C	---------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     S(NWS)    = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     DP(6,6)   = ELASTO-PLASTIC STRESS-STRAIN MATRIX
C     B(3,NEF)  = COMPACTED STRAIN-DISPLACEMENT MATRIX
C     NEF       = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     DVOL      = INTEGRATION CONSTANT (WT*DET)
C     ----------------------------------------------------------------
      DIMENSION S(1),DP(6,6),BB(6,24),DB(6)
	DIMENSION STIFF(24,24)
C
      STIFF=MATMUL(MATMUL(TRANSPOSE(BB),DP),BB)*DVOL
	KS=0                                 
      DO I=1,24                    
	   DO J=I,24                           
	      KS=KS+1 		                     
              S(KS)=STIFF(I,J)+S(KS)                 
         END DO
      END DO 
      RETURN
      END
c
C========================================================
      SUBROUTINE MATRIXTMliu(R,S,T,TT,EM,MM)
C	--------
C     FOR ENHANCED STRAIN TERM OF ENHANCED STRAIN METHOD
C     COMPUTE COEFFICIENT MATRIX EM(6,30)=T*M
C     IN MASTER OR LOCAL COORDINATE SYSTEM
C	--------
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION EM(6,MM),TT(6,6)
      DO 100 I=1,6
	EM(I,1)=TT(I,1)*R
	EM(I,2)=TT(I,2)*S
	EM(I,3)=TT(I,3)*T
	EM(I,4)=TT(I,4)*R
	EM(I,5)=TT(I,4)*S
 100	CONTINUE
      RETURN
      END
C=====================================================================
C	END 3D SOLID ELEMENT
C=====================================================================
	
	SUBROUTINE SOLICO(PROPM,PROPG,NODEX,WA,WA2,S,COORDI,EDIS,EDISI,RE
     1           ,MWG,ALPHA,SEL,SEDI,FIN,HINFC)
	  
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
C     XY(3,8)       = LOCAL NODAL COORDINATES
C     ----------------------------------------------------------------

      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /GAUS/  GLOC(10,10),GWT(10,10),NGR,NGS,NGT
      COMMON /FTIM/  TIM(20),IDATE,ITIME
      COMMON /FLAG/  IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

C	EAS
	COMMON /MMENH/ MM,MM1,MM2,NDIMC                                !***************
	   
	   
      COMMON A(9000000),IA(9000000)
      
      COMMON /IELEM/ IEL
      
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      

C      COMMON /NEWWA/WA2(6,8) !---------NEW WORKING ARRAY BY BJ
C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(300),WAA(MWG,1)
      DIMENSION COORDI(3,8),COORD(3,8)
      DIMENSION EDIS(24),EDISO(24),EDISI(24),RE(24),XJ(9)
      DIMENSION DP(6,6),H(21),P(3,21),XJI(9),B(63),DISD(9)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6),GSTRESS(6),FSTRESS(6)
      DIMENSION WA2(MWG,1) !---------NEW WORKING ARRAY BY BJ
C
      DIMENSION VR(3),VS(3),VT(3)
	DIMENSION XY(3,8),BB(6,24)   
	DIMENSION EDISLO(24),SS(300),FE(24)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)
	DIMENSION REP(24),SEP1(24,24),SEP2(24,24),EDISP(24),MPNOD(8),MPDIS(24),MPGUS(8)
C	EAS
	DIMENSION TTO(6,6),GE(6,MM),SG(300)
	DIMENSION REAS(24),SEAS(24,24)
	DIMENSION EAS(6),ALPHA(MM)
      DIMENSION SED(MM,MM),SEL(MM,24),SEDI(MM,MM),RH(MM),D(6,6)
	DIMENSION HINFC(MM),SEL1(MM,24)	

C     REMOVE RIGID BODY MOTION
      DIMENSION COORDO(3,8),REDIS(24)
      
      DIMENSION DIRCOS(9,1),VR1(3),VS1(3),VT1(3),CSTRESS(6)
 
CC 	ALLOCATABLE WA2(:,:)     
      
C     ------------------------------------------------------
C	DETECTED THE THIN DIRECTION
C     ------------------------------------------------------
	XYZ1 = SQRT( (COORDI(1,2)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,2)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,2)-COORDI(3,1))**2.0 )
	XYZ2 = SQRT( (COORDI(1,4)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,4)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,4)-COORDI(3,1))**2.0 )
	XYZ3 = SQRT( (COORDI(1,5)-COORDI(1,1))**2.0 +
	1	         (COORDI(2,5)-COORDI(2,1))**2.0 +
	2	         (COORDI(3,5)-COORDI(3,1))**2.0 )

	NMX = 0
	IF(XYZ1.LT.XYZ2.AND.XYZ1.LT.XYZ3) NMX = 1
	IF(XYZ2.LT.XYZ1.AND.XYZ2.LT.XYZ3) NMX = 2
	IF(XYZ3.LT.XYZ1.AND.XYZ3.LT.XYZ2) NMX = 3

	IF(NMX.EQ.0) THEN
	IF(XYZ1.EQ.XYZ2) NMX = 1
	IF(XYZ1.EQ.XYZ3) NMX = 1
	IF(XYZ2.EQ.XYZ3) NMX = 2
	ENDIF

	CALL SLCOFACE(NMX,MPNOD,MPDIS,MPGUS)
	
C     FLIPPING COORDINATE
      DO I = 1,NNO
      NOD = MPNOD(I)
      COORD(1:3,I) = COORDI(1:3,NOD)
      ENDDO
      
C     FLIPPING DISPLACEMENT
      EDISP(1:24) = EDIS(1:24)
      DO I = 1,24
      IEF = MPDIS(I)
      EDISO(IEF) = EDISP(I)
      ENDDO
C     ------------------------------------------------------------
C	FIND THE LOCAL VECTORS,LOCAL COORDINATES AND DISPLACEMENTS
C     -----------------------------------------------------------      
      CALL SHAP3D (0.0,0.0,0.0,H,P,NODEX,NNO)
      CALL SOLIRST(COORD,P,VR,VS,VT)   ! LOCAL VECTORS
          	
C     REMOVE RIGID BODY MOTION
  	REDIS(1:24) = EDISO(1:24)
      IF(NLOPT.EQ.3) THEN
          K = 0
          DO I = 1,NNO
              DO J = 1,3
                  K = K + 1
                  COORDO(J,I) = COORD(J,I) - EDISO(K)
              ENDDO
          ENDDO
          CALL SOMDSP(COORD,COORDO,EDISO,REDIS,H,P,VR,VS,VT,NNO)
      ENDIF
C     REDIS = DISPLACEMENT REMOVING RIGID BOBY      
      CALL SCOOR(COORD,REDIS,XY,EDISLO,VR,VS,VT) !REDIS INSTEAD OF EDIS HERE 
C     --------------------------------------------
C     EAS TRANSFORMATION MATRIX
C     --------------------------------------------	

	CALL JACO3D (XY,P,XJ,XJI,DETO,MEL,NNO)	! *****************
	CALL TNEAS6(XJ,TTO)                       ! *****************	 

C     SET ZEROS STIFFNESS FOR EAS	
	SED  = 0.0                                 ! ***************
	SEDI = 0.0                                 ! ***************
	SEL  = 0.0                                 ! ***************
	RH   = 0.0                                 ! ***************	
	SG = 0.0 
	HINFC = 0.0 
C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------	
	CALL HOKLAW (PROPM,PROPG,1)
      MGR = NGR
      MGS = NGS
      MGT = NGT
      IF (ITASK.NE.5) GOTO 10
      MGR = 3
      MGS = 3
      MGT = 3
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
 10   IPT = 0
      FE =0.0
	S =0.0
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
C     ---------------------------------------------------
      CALL SHAP3D (RI,SI,TI,H,P,NODEX,NNO)
C     -----------------------------------------------
C     SETUP THE COROTATIONAL COORDINATE (VR,VS,VT)
C     -----------------------------------------------
 !	CALL SOLIRST(COORD,P,VR,VS,VT)
C     -----------------------------------------------------------------
C      TRANSLATE THE GLOBAL NODAL COORDINATES/DISPLACEMENTS INTO LOCAL 
C     -----------------------------------------------------------------
! 	CALL SCOOR(COORD,EDISO,XY,EDISLO,VR,VS,VT)
C     -------------------------------------------------------------------
C	JACOBIAN (XJ), INVERSE OF THE JACOBIAN (XJI) AND DETERMINANT (DET)
C     -------------------------------------------------------------------	    
  	    
      CALL JACO3D (XY,P,XJ,XJI,DET,MEL,NNO)	     
      DVOL = WT*DET	
      
C     EAS TERMS    
      CALL MTEAS7(RI,SI,TI,TTO,DETO,DET,GE,MM)        !***************
      
C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
 !     IF (ITASK.NE.5)  GOTO 50
 !     CALL SOMASS (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
  !    GOTO 900
C     ----------------------------------------
C     COMPACTED STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------
 50   CALL SOBMAT (P,XJI,B,NNO)
C     ----------------------------------------
C      ASSUMED STRAIN-DISPLACEMENT MATRIX (BB)
C     ----------------------------------------       
C      CALL SOBMATA(XY,BB,B,RI,SI,NNO,MEL,NODEX)       ! OLD ANS METHOD
      CALL SOBMATANS(XY,BB,B,RI,SI,TI,NNO,MEL,NODEX)   ! NEW ANS METHOD	
      IF (NLOPT+ITASK.EQ.1)  GOTO 700
C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------
      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,6
 100  TAU(I) = WA(I,IPT)*DVOL
      GOTO 800
C     -------------------------------
C     DISPLACEMENT DERIVATIVES (DISD)
C     -------------------------------
 200  DISD=0.0
      DO 290  IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDISLO(IEF)
      DISD(2) = DISD(2) + B(JEF)*EDISLO(JEF)
      DISD(3) = DISD(3) + B(KEF)*EDISLO(KEF)
      DISD(4) = DISD(4) + B(JEF)*EDISLO(IEF)
      DISD(5) = DISD(5) + B(KEF)*EDISLO(IEF)
      DISD(6) = DISD(6) + B(IEF)*EDISLO(JEF)
      DISD(7) = DISD(7) + B(KEF)*EDISLO(JEF)
      DISD(8) = DISD(8) + B(IEF)*EDISLO(KEF)
 290  DISD(9) = DISD(9) + B(JEF)*EDISLO(KEF)
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)
      STRAIN(4) = DISD(4) + DISD(6)
C      STRAIN(5) = DISD(5) + DISD(8)
C      STRAIN(6) = DISD(7) + DISD(9)
      STRAIN(5)=0.0
	STRAIN(6)=0.0
	DO I=1,24     
      STRAIN(5)=STRAIN(5)+BB(5,I)*EDISLO(I)
	STRAIN(6)=STRAIN(6)+BB(6,I)*EDISLO(I)
      ENDDO
      
C     -------------------------------
C     FOR EAS TERM  EAS=[GE]*{ALPHA}
C     -------------------------------
	EAS = MATMUL(GE,ALPHA)                  ! -------------------------
	
	DO I=1,6
       STRAIN(I)=STRAIN(I) - EAS(I)            ! -------------------------
	ENDDO    
	
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
      DO 390  I=1,6
 390  STRAIN(I) = STRAIN(I) - QSTRAI(I)
C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------
 400  IF (MTMOD-3)  405,450,460
 405  CONTINUE
   !   CALL SOLSIG (STRAIN,STRESS)
      CALL SOLSIGCO(STRAIN,STRESS)
      
CC    --------------------------------------------- NEW WAKING ARRAY WA2 BY BJ

C	GROUP FILE

CC      NFL1 = 808 
cc      NFIL = 803
      
C	OPEN FILE FOR STORE GROUP DATA
cc	CALL INTFILL('NOUT',NFIL,1,3,0)   !ELEMENT GROUP STORAGE FILE
	
cc      DO 408 IEG = 1,NEG

CC      IEG = NEG
      
cc      NFL5 = 5 + 4*(IEG-1) + NFIL
      
cc	CALL INTFILL('OGRF',NFL ,11,IEG,1) !

cc408   CONTINUE    

      
	
cc	CALL INTFILL('OGRF',NFL5 ,11,KEG,0) !
      
cc	ALLOCATE(WA2(MWG,8))
cc	REWIND(NFL5)
cc	READ(NFL5) WA2 
  
  
      ! ------------------- STORE STRESSES --------------------------------------
      CALL STRESSTRANS(VR,VS,VT,STRESS,GSTRESS)
cc      CALL STRESSLOCALTRANS(VR1,VS1,VT1,GSTRESS,CSTRESS) 
      NOD = MPGUS(IPT)    
      DO 410  I=1,6
      WA(I,NOD)  = GSTRESS(I) !GSTRESS = GLOBAL STRESS
 410  WA2(I,NOD)  = STRESS(I) !SOLID-SHELL(P) FOR LOCAL STRESS
 
!---------------------------------------------------------------FOR WA2 BY BJ
cc      REWIND(NFL5)
cc	WRITE(NFL5) WA2
	
cc	DEALLOCATE(WA2)
!---------------------------------------------------------------FOR WA2 BY BJ
      GOTO 500
      
      CALL SOLSIG_S (STRAIN,STRESS)
      DO 415  I=1,6
 415  WA(I,IPT) = STRESS(I)
      GOTO 500      
            
C     -----------------------------------------------------------------------------
 450  CALL MISE3D (WA(1,IPT),WA(7,IPT),WA(13,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DP)
      GOTO 500
C
C460  CALL CONC3D -- DELETED BECAUSE OF SUBROUTINE NOT EXIST
 460  CONTINUE
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510  I=1,6
 510  TAU(I) = STRESS(I)*DVOL
C      IF (ITASK.LE.2) GOTO 520
      IF (ITASK.LE.3) GOTO 520
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900

 520  CONTINUE
   
      FE = FE + MATMUL(TRANSPOSE(BB),TAU)
      RH = RH + MATMUL(TRANSPOSE(GE),TAU)           !***************

C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------
   !   IF (IFREF) 900,700,900
 700  CONTINUE  
C
      CALL SOK0NLCO(S,BB,D,DVOL)  !      
    !  IF (MTMOD.GT.2) CALL SOK0PL (S,DP,B,NEF,DVOL)
C     -------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------
 790  IF (NLOPT.LE.1) GOTO 810
 800  CONTINUE
C      CALL SOKSIG (S,TAU,B,NEF)
      CALL KNLSTIFF (S,TAU,B)
 810  CONTINUE !TIM(12) = TIM(12) + TIM2-TIM1
 
	SED = SED + MATMUL(TRANSPOSE(GE),MATMUL(D,GE))*DVOL           !***************
	SEL = SEL + MATMUL(TRANSPOSE(GE),MATMUL(D,BB))*DVOL           !*************** 	
      
 900  CONTINUE
C     ----------------------------------------------------------------     
	DO I = 1,MM
	  HINFC(I) = RH(I)                             !***************
	ENDDO
	    
	CALL INVMATRIX(SED,SEDI,MM)                                   !***************

	REAS = MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,RH))
	SEAS = MATMUL(TRANSPOSE(SEL),MATMUL(SEDI,SEL))
	
	K = 0
	DO I = 1,24		                      
	  DO J = I,24
	    K = K + 1
C	    S(K) = S(K) + SG(K) - SEAS(I,J)                            !*************** 
	    S(K) = S(K) - SEAS(I,J)                            !*************** 
	  ENDDO
	ENDDO 	
	  			
	DO I = 1,24
	  RE(I) = FE(I) + REAS(I)                      !***************
	ENDDO 
	
	CALL STIFFTRANS(VR,VS,VT,S,RE)	
	CALL STIFFTRANE(VR,VS,VT,MM,SEL) !TRANSFORM EAS MATRIX TO GLOBAL
C     ----------------------------------------------------------------       
      
C     ---------------------------------------
C     FLIPPING DUE TO THICKNESS DETECTION
C     ---------------------------------------
      REP(1:24)  = RE(1:24)
      DO IEF = 1,24
          II = MPDIS(IEF)
          RE(IEF) = REP(II)
      ENDDO
      
      SEP1(1:24,1:24) =  0.0D0
      SEP2(1:24,1:24) =  0.0D0
      KEF = 0
      DO IEF = 1,24
          DO JEF = IEF,24
              KEF = KEF+1
              SEP1(IEF,JEF) = S(KEF)
              SEP1(JEF,IEF) = S(KEF)
              SEP2(IEF,JEF) = S(KEF)
              SEP2(JEF,IEF) = S(KEF)
          ENDDO
      ENDDO

      DO IEF = 1,24
          II = MPDIS(IEF)
          DO JEF = 1,24
              JJ = MPDIS(JEF)
              SEP1(IEF,JEF) = SEP2(II,JJ) 
          ENDDO
      ENDDO
      
C     REARRANGE EAS MATRIX      
      SEL1 = SEL
      DO JEF = 1,24
          JJ = MPDIS(JEF)
          SEL(1:MM,JEF) = SEL1(1:MM,JJ)
      ENDDO
      
      KEF = 0
      DO IEF = 1,24
          DO JEF = IEF,24
              KEF = KEF+1
              S(KEF) = SEP1(IEF,JEF)
          ENDDO
      ENDDO      
      
  
          
C     ---------------------------------------
C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF

      RETURN
      END
C	=====================================================================
	SUBROUTINE TNEAS6(XJO,TTO)
C	=====================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION XJO(3,3),TT(6,6),TTO(6,6),TTI(6,6)

C     The Transformation Matrix of JACOBIAN     
      TT(1,1) = XJO(1,1)*XJO(1,1)
      TT(1,2) = XJO(2,1)*XJO(2,1)
      TT(1,3) = XJO(3,1)*XJO(3,1)
      TT(1,4) = 2*XJO(1,1)*XJO(2,1)
      TT(1,5) = 2*XJO(1,1)*XJO(3,1)
      TT(1,6) = 2*XJO(2,1)*XJO(3,1)
    
      TT(2,1) = XJO(1,2)*XJO(1,2) 
      TT(2,2) = XJO(2,2)*XJO(2,2)
      TT(2,3) = XJO(3,2)*XJO(3,2)
      TT(2,4) = 2*XJO(1,2)*XJO(2,2)
      TT(2,5) = 2*XJO(1,2)*XJO(3,2)
      TT(2,6) = 2*XJO(2,2)*XJO(3,2)
      
      TT(3,1) = XJO(1,3)*XJO(1,3)  
      TT(3,2) = XJO(2,3)*XJO(2,3)
      TT(3,3) = XJO(3,3)*XJO(3,3)
      TT(3,4) = 2*XJO(1,3)*XJO(2,3)
      TT(3,5) = 2*XJO(1,3)*XJO(3,3)
      TT(3,6) = 2*XJO(2,3)*XJO(3,3) 
      
      TT(4,1) = XJO(1,1)*XJO(1,2)  
      TT(4,2) = XJO(2,1)*XJO(2,2)
      TT(4,3) = XJO(3,1)*XJO(3,2)
      TT(4,4) = XJO(1,2)*XJO(2,1)+XJO(1,1)*XJO(2,2)
      TT(4,5) = XJO(1,2)*XJO(3,1)+XJO(3,2)*XJO(1,1)
      TT(4,6) = XJO(2,2)*XJO(3,1)+XJO(3,2)*XJO(2,1)
      
      TT(5,1) = XJO(1,1)*XJO(1,3) 
      TT(5,2) = XJO(2,1)*XJO(2,3)
      TT(5,3) = XJO(3,1)*XJO(3,3)
      TT(5,4) = XJO(1,3)*XJO(2,1)+XJO(2,3)*XJO(1,1)
      TT(5,5) = XJO(1,3)*XJO(3,1)+XJO(3,3)*XJO(1,1)
      TT(5,6) = XJO(2,3)*XJO(3,1)+XJO(3,3)*XJO(2,1)
      
      TT(6,1) = XJO(1,2)*XJO(1,3) 
      TT(6,2) = XJO(2,2)*XJO(2,3)
      TT(6,3) = XJO(3,2)*XJO(3,3)
      TT(6,4) = XJO(1,3)*XJO(2,2)+XJO(2,3)*XJO(1,2)
      TT(6,5) = XJO(1,3)*XJO(3,2)+XJO(3,3)*XJO(1,2)
      TT(6,6) = XJO(2,3)*XJO(3,2)+XJO(3,3)*XJO(2,2)              
      
      ! Inverse The Transformation Matrix 
      CALL INVMATRIX(TT,TTI,6) 
      ! Transpose of The Inverse of The Transformation Matrix                       
      TTO = TRANSPOSE(TTI)   
      
	RETURN
	END
C	=====================================================================
C	=====================================================================
	SUBROUTINE MTEAS7(RN,SN,TN,TTO,DETO,DET,GE,MM)          
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION TTO(6,6),GE(6,MM),TITO(6,6),GEO(6,MM)	

	CALL EASGE7(RN,SN,TN,GEO,MM)

	DO I = 1,6
	    DO J = 1,6
	      TITO(I,J) = (DETO/DET)*TTO(I,J)
	    ENDDO
	ENDDO

	GE = MATMUL(TITO,GEO)

	RETURN
	END

C	=====================================================================
C	=====================================================================
	SUBROUTINE EASGE7(RN,SN,TN,GEO,MM)
C	=====================================================================
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=====================================================================
	DIMENSION GEO(6,MM) 
	
	IF (MM == 5) THEN   

	GEO = 0.0D0	
	GEO(1,1) = RN
	GEO(2,2) = SN	
	GEO(3,3) = TN
	GEO(4,4) = RN
	GEO(4,5) = SN
			
	ELSEIF (MM == 7) THEN   

	GEO = 0.0D0	
	GEO(1,1) = RN
	GEO(2,2) = SN	
	GEO(4,3) = RN
	GEO(4,4) = SN
	GEO(1,5) = RN*SN
	GEO(2,6) = RN*SN
	GEO(4,7) = RN*SN
	
	ELSEIF (MM == 11) THEN  ! GOOD
	 
	GEO = 0.0D0
	GEO(1,1) =  RN
	GEO(2,2) =  SN	
	! ---------------------------------------
	GEO(4,3) =  RN
	GEO(4,4) =  SN
	GEO(5,5) =  RN
	GEO(6,6) =  SN
	! ---------------------------------------	
	GEO(5,7) =  RN*SN
	GEO(6,8) =  RN*SN	
      ! ---------------------------------------
	GEO(1,9)  = RN*SN
	GEO(2,10) = RN*SN
	GEO(4,11) = RN*SN
	
	ELSEIF (MM == 15) THEN  ! GOOD
	 
	GEO = 0.0D0
	GEO(1,1) =  RN
	GEO(2,2) =  SN	
	! ---------------------------------------
	GEO(4,3) =  RN
	GEO(4,4) =  SN
	GEO(5,5) =  RN
	GEO(5,6) =  TN	
	GEO(6,7) =  SN
	GEO(6,8) =  TN	
      ! --------------------------------------- 	
	GEO(5,9) = RN*SN 
	GEO(5,10) = SN*TN
	GEO(6,11) = RN*SN 	
	GEO(6,12) = RN*TN
	! ---------------------------------------
	GEO(1,13) = RN*SN
	GEO(2,14) = RN*SN
	GEO(4,15) = RN*SN
	
	ELSEIF (MM == 18) THEN  
	 
	GEO = 0.0D0
	GEO(1,1) =  RN
	GEO(2,2) =  SN
	GEO(3,3) =  TN
	! ---------------------------------------
	GEO(4,4) =  RN
	GEO(4,5) =  SN
	GEO(5,6) =  RN
	GEO(5,7) =  TN	
	GEO(6,8) =  SN
	GEO(6,9) =  TN	
      ! --------------------------------------- 	
	GEO(4,10) = SN*TN
	GEO(4,11) = RN*TN	
	GEO(5,12) = RN*SN 
	GEO(5,13) = SN*TN
	GEO(6,14) = RN*SN 	
	GEO(6,15) = RN*TN
	
	GEO(1,16) = RN*SN
	GEO(2,17) = RN*SN
	GEO(4,18) = RN*SN
      
	
	ELSEIF (MM == 24) THEN  ! BACK UP
	 
	GEO = 0.0D0
	
	GEO(1,1) =  RN
	GEO(2,2) =  SN
	GEO(3,3) =  TN
	GEO(4,4) =  RN
	GEO(4,5) =  SN

	GEO(5,6) =  RN
	GEO(5,7) =  TN	
	GEO(6,8) =  SN
	GEO(6,9) =  TN	

	GEO(4,10) = RN*SN 
	GEO(4,11) = SN*TN
	GEO(4,12) = RN*TN
	
	GEO(5,13) = RN*SN 
	GEO(5,14) = SN*TN
	GEO(5,15) = RN*TN

	GEO(6,16) = RN*SN 
	GEO(6,17) = SN*TN
	GEO(6,18) = RN*TN

	GEO(1,19) = RN*TN 
	GEO(1,20) = RN*SN
	 
	GEO(2,21) = RN*SN 
	GEO(2,22) = SN*TN 

	GEO(3,23) = SN*TN 
	GEO(3,24) = RN*TN
		 
	ENDIF
      
	RETURN
	END
C	=====================================================================      
C****************************************************************
      SUBROUTINE SLCOFACE(NMX,MPNOD,MPDIS,MPGUS)                 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ============================================
C     DEFINE THE FACES FOR THE 8 NODE SOLID ELEMET
C     ============================================
      DIMENSION MPNOD(8),MPDIS(24),MPGUS(8)

C     FLIPING NODES
	MPNOD = 0
	IF(NMX.EQ.1) THEN !14852376
	MPNOD(1 ) = 1  
	MPNOD(2 ) = 4
	MPNOD(3 ) = 8
	MPNOD(4 ) = 5  
	MPNOD(5 ) = 2
	MPNOD(6 ) = 3
	MPNOD(7 ) = 7  
	MPNOD(8 ) = 6
	ELSEIF(NMX.EQ.2) THEN !21563487
	MPNOD(1 ) = 2 
	MPNOD(2 ) = 1
	MPNOD(3 ) = 5
	MPNOD(4 ) = 6  
	MPNOD(5 ) = 3
	MPNOD(6 ) = 4
	MPNOD(7 ) = 8  
	MPNOD(8 ) = 7
	ELSEIF(NMX.EQ.3) THEN !12345678
	MPNOD(1 ) = 1  
	MPNOD(2 ) = 2
	MPNOD(3 ) = 3
	MPNOD(4 ) = 4  
	MPNOD(5 ) = 5
	MPNOD(6 ) = 6
	MPNOD(7 ) = 7  
	MPNOD(8 ) = 8
	ENDIF
		
C     FLIPING DISPLACEMENTS		
	MPDIS = 0
	IF(NMX.EQ.1) THEN !14852376
	MPDIS(1 ) = 1  !NODE 1
	MPDIS(2 ) = 2
	MPDIS(3 ) = 3
	MPDIS(4 ) = 13 !NODE 2
	MPDIS(5 ) = 14
	MPDIS(6 ) = 15
	MPDIS(7 ) = 16 !NODE 3
	MPDIS(8 ) = 17
	MPDIS(9 ) = 18
	MPDIS(10) = 4  !NODE 4
	MPDIS(11) = 5
	MPDIS(12) = 6
	MPDIS(13) = 10 !NODE 5
	MPDIS(14) = 11
	MPDIS(15) = 12
	MPDIS(16) = 22 !NODE 6
	MPDIS(17) = 23
	MPDIS(18) = 24
	MPDIS(19) = 19 !NODE 7
	MPDIS(20) = 20
	MPDIS(21) = 21
	MPDIS(22) = 7  !NODE 8
	MPDIS(23) = 8
	MPDIS(24) = 9
	ELSEIF(NMX.EQ.2) THEN !21563487
	MPDIS(1 ) = 4  !NODE 1
	MPDIS(2 ) = 5
	MPDIS(3 ) = 6
	MPDIS(4 ) = 1  !NODE 2
	MPDIS(5 ) = 2
	MPDIS(6 ) = 3
	MPDIS(7 ) = 13 !NODE 3
	MPDIS(8 ) = 14
	MPDIS(9 ) = 15
	MPDIS(10) = 16 !NODE 4
	MPDIS(11) = 17
	MPDIS(12) = 18
	MPDIS(13) = 7  !NODE 5
	MPDIS(14) = 8
	MPDIS(15) = 9
	MPDIS(16) = 10 !NODE 6
	MPDIS(17) = 11
	MPDIS(18) = 12
	MPDIS(19) = 22 !NODE 7
	MPDIS(20) = 23
	MPDIS(21) = 24
	MPDIS(22) = 19 !NODE 8
	MPDIS(23) = 20
	MPDIS(24) = 21
	ELSEIF(NMX.EQ.3) THEN !12345678
	MPDIS(1 ) = 1  !NODE 1
	MPDIS(2 ) = 2
	MPDIS(3 ) = 3
	MPDIS(4 ) = 4  !NODE 2
	MPDIS(5 ) = 5
	MPDIS(6 ) = 6
	MPDIS(7 ) = 7  !NODE 3
	MPDIS(8 ) = 8
	MPDIS(9 ) = 9
	MPDIS(10) = 10 !NODE 4
	MPDIS(11) = 11
	MPDIS(12) = 12
	MPDIS(13) = 13 !NODE 5
	MPDIS(14) = 14
	MPDIS(15) = 15
	MPDIS(16) = 16 !NODE 6
	MPDIS(17) = 17
	MPDIS(18) = 18
	MPDIS(19) = 19 !NODE 7
	MPDIS(20) = 20
	MPDIS(21) = 21
	MPDIS(22) = 22 !NODE 8
	MPDIS(23) = 23
	MPDIS(24) = 24
	ENDIF

C     FLIPING GUASS POINTS	
	MPGUS = 0
	IF(NMX.EQ.1) THEN !14852376
	MPGUS(1) = 1 
	MPGUS(3) = 2
	MPGUS(5) = 3 
	MPGUS(7) = 4  
	MPGUS(2) = 5
	MPGUS(4) = 6
	MPGUS(6) = 7  
	MPGUS(8) = 8
	ELSEIF(NMX.EQ.2) THEN !21563487
	MPGUS(5) = 1 
	MPGUS(7) = 2
	MPGUS(6) = 3 
	MPGUS(8) = 4  
	MPGUS(1) = 5
	MPGUS(3) = 6
	MPGUS(2) = 7  
	MPGUS(4) = 8
	ELSEIF(NMX.EQ.3) THEN !12345678
	MPGUS(1) = 1  
	MPGUS(2) = 2
	MPGUS(3) = 3
	MPGUS(4) = 4  
	MPGUS(5) = 5
	MPGUS(6) = 6
	MPGUS(7) = 7  
	MPGUS(8) = 8
	ENDIF

	RETURN
	END
C****************************************************************      
      
C****************************************************************
C     ---------------------------------------------------------      
      SUBROUTINE SOLIRST(COORD,P,VR,VS,VT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------
C     SUBROUTINE TO FIND THE COROTATIONAL VECTORS
C     --------------------------------------------
C      
C      VT(3) = UNIT VECTOR NORMAL TO THE MIDSURFACE
C      VR(3),VS(3)=UNIT VECTORS TANGENTIAL TO MIDSURFACE
C     COVR(3),COVS(3)   = COVARIENT BASE VECTORS ALONG RN,SN
C     ---------------------------------------------------
      DIMENSION COORD(3,8),P(3,8)
	DIMENSION VR(3),VS(3),VT(3)
	DIMENSION COVR(3),COVS(3)
C
      COVR=0.
	COVS=0.
      DO I=1,8
	DO J=1,3
      COVR(J)=COVR(J)+P(1,I)*COORD(J,I)
      COVS(J)=COVS(J)+P(2,I)*COORD(J,I)
	ENDDO
	ENDDO
C
      CALL VECPRD (COVR,COVS,VT)
	CALL SCALEN (VT,VT,DD,3)
	CALL VECPRD (VT,COVR,VS)
	CALL SCALEN (VS,VS,DD,3)
	CALL SCALEN (COVR,VR,DD,3)
C	
	RETURN
	END
C****************************************************************
      SUBROUTINE SCOOR(COORD,EDIS,XY,EDISLO,VR,VS,VT)  
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------
C      CALCULATE COROTATIONAL FORM OF NODAL COORDINATES
C      AND NODAL DISPLACEMENTS
C     --------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      DIMENSION COORD(3,8),XY(3,8),VR(3),VS(3),VT(3)
	DIMENSION EDIS(24),EDISLO(24)
C
      XY=0.0
	EDISLO=0.0
	DO I=1,8
	DO J=1,3
      XY(1,I)=XY(1,I)+COORD(J,I)*VR(J)
      XY(2,I)=XY(2,I)+COORD(J,I)*VS(J)
	XY(3,I)=XY(3,I)+COORD(J,I)*VT(J)
C	
    !  IF(NLOPT.LE.1)GOTO 100
	EDISLO(3*I-2)=EDISLO(3*I-2)+EDIS(3*I+J-3)*VR(J)
	EDISLO(3*I-1)=EDISLO(3*I-1)+EDIS(3*I+J-3)*VS(J)
	EDISLO(3*I)  =EDISLO(3*I)  +EDIS(3*I+J-3)*VT(J)
  100 CONTINUE
	ENDDO
	ENDDO
C
      RETURN
	END
C*****************************************************************
      SUBROUTINE  SOBMATA(XY,BB,B,RI,SI,NNO,MEL,NODEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C    -----------------------------------------------------
C     COMPUTE THE STABLE B MATRICES IN TAYLOR SERIES FORM
C      USING ASSUMED STRAIN METHOD
C    -----------------------------------------------------
C     B1,B2,B3,B4,B5,B6 = B MATRIX DERIVATIVES WITH 
C                         RESPECT TO R,S,T,RS,RT,ST
C     X(8),Y(8),Z(8) = NODAL COORDINATES
C     COORD(3,NNO)   = NODAL COORDINATES
C     PP(24)         = UNIFORM B MATRIX IN ONE DIMENSION
C     GAMMA1-4       = THE GAMMA VECTORS
C     HD(4,3)        = THE 4 H DERIVATIVES WITH RESPECT TO X,Y,Z
C    -----------------------------------------------------
C
      DIMENSION COORD(3,NNO),XJI(3,3),XJ(3,3)
      DIMENSION BB(6,24),B(24)  
      DIMENSION PTA(8),PSA(8),PTB(8),PRB(8)
	DIMENSION PTC(8),PSC(8),PTD(8),PRD(8)
	DIMENSION PTA1(8),PTA2(8),PTA3(8),PTA4(8)
	DIMENSION RL(4),SL(4)
	DIMENSION TT(4)
C
	CALL SHDSAMP(PSA,PTA,PRB,PTB,PSC,PTC,PRD,PTD,
     @                   PTA1,PTA2,PTA3,PTA4,RL,SL)
      CALL SHAP3D(RI,SI,0.0,H,P,NODEX,NNO)      
	CALL JACO3D (XY,P,XJ,XJI,DET,MEL,NNO)
C     ----DEFINE THE TRANSFORM MATRIX TT
      TT(1)=XJI(1,1)*XJI(3,3)+XJI(3,1)*XJI(1,3)
  	TT(2)=XJI(1,2)*XJI(3,3)+XJI(3,2)*XJI(1,3)
	TT(3)=XJI(2,1)*XJI(3,3)+XJI(3,1)*XJI(2,3)
	TT(4)=XJI(2,2)*XJI(3,3)+XJI(3,2)*XJI(2,3)
C
      BB=0.
C
	DO I=1,8
      BB(1,3*I-2)=B(3*I-2)
	BB(2,3*I-1)=B(3*I-1)
	BB(3,3*I)  =B(3*I)
	BB(4,3*I-2)=B(3*I-1)
	BB(4,3*I-1)=B(3*I-2)
	BB(5,3*I-2)=B(3*I)
	BB(5,3*I)  =B(3*I-2)
	BB(6,3*I-1)=B(3*I)
	BB(6,3*I)=B(3*I-1)
C     ----MODIFY DUE TO ASSUMED SHEAR STRAIN--------
      BDT=(PTB(I)*(1-SI)+(1+SI)*PTD(I))*.5
      BDR=(PRB(I)*(1-SI)+(1+SI)*PRD(I))*.5
      ACT=(PTA(I)*(1-RI)+(1+RI)*PTC(I))*.5
      ACS=(PSA(I)*(1-RI)+(1+RI)*PSC(I))*.5
      !BDT=(PTB(I)*(1-RI)+(1+RI)*PTD(I))*.5
      !BDR=(PRB(I)*(1-RI)+(1+RI)*PRD(I))*.5
      !ACT=(PTA(I)*(1-SI)+(1+SI)*PTC(I))*.5
      !ACS=(PSA(I)*(1-SI)+(1+SI)*PSC(I))*.5
      BB(5,3*I-2) = TT(1)*BDT*XJ(1,1)+TT(2)*ACT*XJ(2,1)+
     @ 	XJ(3,1)*TT(1)*BDR+XJ(3,1)*TT(2)*ACS
	BB(5,3*I)   = TT(1)*BDT*XJ(1,3)+TT(2)*ACT*XJ(2,3)+
     @ 	XJ(3,3)*TT(1)*BDR+XJ(3,3)*TT(2)*ACS
      BB(5,3*I-1)=TT(1)*BDt*XJ(1,2)+TT(2)*ACt*XJ(2,2)
     @ 	+(TT(1)*BDr+TT(2)*ACs)*XJ(3,2) 
	BB(6,3*I-1) = TT(3)*BDT*XJ(1,2)+TT(4)*ACT*XJ(2,2)+
     @	XJ(3,2)*TT(3)*BDR+XJ(3,2)*TT(4)*ACS
      BB(6,3*I)   = TT(3)*BDT*XJ(1,3)+TT(4)*ACT*XJ(2,3)+
     @	XJ(3,3)*TT(3)*BDR+XJ(3,3)*TT(4)*ACS	
	BB(6,3*I-2)=TT(3)*BDt*XJ(1,1)+TT(4)*ACt*XJ(2,1)
     @	+(TT(3)*BDr+TT(4)*ACs)*XJ(3,1) 	
	 
C     ----MODIFY DUE TO ASSUMED THICKNESS STRAIN----
      TH=(PTA1(I)*(1+RL(1)*RI)*(1+SL(1)*SI)	
     @ +PTA2(I)*(1+RL(2)*RI)*(1+SL(2)*SI)
     @ +PTA3(I)*(1+RL(3)*RI)*(1+SL(3)*SI)
     @ +PTA4(I)*(1+RL(4)*RI)*(1+SL(4)*SI))*0.25      
      BB(3,3*I)=XJI(3,3)*XJI(3,3)*TH*XJ(3,3)
 	BB(3,3*I-2)=XJI(3,3)*XJI(3,3)*TH*XJ(3,1)
 	BB(3,3*I-1)=XJI(3,3)*XJI(3,3)*TH*XJ(3,2)     
      ENDDO
C
      RETURN 
	END
C*************************************************************************
C*****************************************************************
      SUBROUTINE  SOBMATANS(XY,BB,B,RI,SI,TI,NNO,MEL,NODEX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C    -----------------------------------------------------
C     COMPUTE THE STABLE B MATRICES IN TAYLOR SERIES FORM
C      USING ASSUMED STRAIN METHOD
C    -----------------------------------------------------
C     B1,B2,B3,B4,B5,B6 = B MATRIX DERIVATIVES WITH 
C                         RESPECT TO R,S,T,RS,RT,ST
C     X(8),Y(8),Z(8) = NODAL COORDINATES
C     COORD(3,NNO)   = NODAL COORDINATES
C     PP(24)         = UNIFORM B MATRIX IN ONE DIMENSION
C     GAMMA1-4       = THE GAMMA VECTORS
C     HD(4,3)        = THE 4 H DERIVATIVES WITH RESPECT TO X,Y,Z
C    -----------------------------------------------------
C
      DIMENSION COORD(3,NNO),XJI(3,3),XJ(3,3),XY(3,NNO),NODEX(1)
      DIMENSION BB(6,24),B(24)  
      DIMENSION PTA(8),PSA(8),PTB(8),PRB(8)
	DIMENSION PTC(8),PSC(8),PTD(8),PRD(8)
	DIMENSION PTA1(8),PTA2(8),PTA3(8),PTA4(8)
	DIMENSION RL(4),SL(4)
	DIMENSION TT(4)
C	ANS
	DIMENSION PAT(3,8),PAB(3,8),PBT(3,8),PBB(3,8),PCT(3,8),PCB(3,8),PDT(3,8),PDB(3,8)
	DIMENSION PE(3,8),PF(3,8),PG(3,8),PH(3,8)
C
	CALL SIMPANS(PAT,PAB,PBT,PBB,PCT,PCB,PDT,PDB,PE,PF,PG,PH,NODEX,NNO)
      CALL SHAP3D(RI,SI,TI,H,P,NODEX,NNO)      
	CALL JACO3D (XY,P,XJ,XJI,DET,MEL,NNO)
C     ----DEFINE THE TRANSFORM MATRIX TT
      TT(1)=XJI(1,1)*XJI(3,3)+XJI(3,1)*XJI(1,3)
  	TT(3)=XJI(1,2)*XJI(3,3)+XJI(3,2)*XJI(1,3)
	TT(2)=XJI(2,1)*XJI(3,3)+XJI(3,1)*XJI(2,3)
	TT(4)=XJI(2,2)*XJI(3,3)+XJI(3,2)*XJI(2,3)
C
      BB=0.
C
	DO I=1,8
      BB(1,3*I-2) = B(3*I-2)
	BB(2,3*I-1) = B(3*I-1)
	BB(3,3*I-0) = B(3*I)
	BB(4,3*I-2) = B(3*I-1)
	BB(4,3*I-1) = B(3*I-2)
	BB(5,3*I-2) = B(3*I)
	BB(5,3*I-0) = B(3*I-2)
	BB(6,3*I-1) = B(3*I)
	BB(6,3*I-0) = B(3*I-1)
C     ----MODIFY DUE TO ASSUMED SHEAR STRAIN--------   
      
      BS11 = (1+SI)*((1-TI)*PAB(3,I)+(1+TI)*PAT(3,I)) + (1-SI)*((1-TI)*PCB(3,I)+(1+TI)*PCT(3,I))
      BS13 = (1+SI)*((1-TI)*PAB(1,I)+(1+TI)*PAT(1,I)) + (1-SI)*((1-TI)*PCB(1,I)+(1+TI)*PCT(1,I))
      BS22 = (1+RI)*((1-TI)*PDB(3,I)+(1+TI)*PDT(3,I)) + (1-RI)*((1-TI)*PBB(3,I)+(1+TI)*PBT(3,I))
      BS23 = (1+RI)*((1-TI)*PDB(2,I)+(1+TI)*PDT(2,I)) + (1-RI)*((1-TI)*PBB(2,I)+(1+TI)*PBT(2,I))     
      
      GAM11 = 0.25*(XJ(1,1)*BS11 + XJ(3,1)*BS13)
      GAM12 = 0.25*(XJ(1,2)*BS11 + XJ(3,2)*BS13)
      GAM13 = 0.25*(XJ(1,3)*BS11 + XJ(3,3)*BS13)
      GAM21 = 0.25*(XJ(2,1)*BS22 + XJ(3,1)*BS23)
      GAM22 = 0.25*(XJ(2,2)*BS22 + XJ(3,2)*BS23)
      GAM23 = 0.25*(XJ(2,3)*BS22 + XJ(3,3)*BS23)
      
      BB(5,3*I-2) = TT(1)*GAM11+TT(2)*GAM21
      BB(5,3*I-1) = TT(1)*GAM12+TT(2)*GAM22     
	BB(5,3*I-0) = TT(1)*GAM13+TT(2)*GAM23
      
      BB(6,3*I-2) = TT(3)*GAM11+TT(4)*GAM21    
	BB(6,3*I-1) = TT(3)*GAM12+TT(4)*GAM22   
      BB(6,3*I-0) = TT(3)*GAM13+TT(4)*GAM23	      
	 
C     ----MODIFY DUE TO ASSUMED THICKNESS STRAIN----
      BT = (1-RI)*(1-SI)*PE(3,I) 
     @    +(1+RI)*(1-SI)*PF(3,I)  
     @    +(1+RI)*(1+SI)*PG(3,I)  
     @    +(1-RI)*(1+SI)*PH(3,I)  
     
      BBT11 = 0.25*XJ(3,1)*BT
      BBT12 = 0.25*XJ(3,2)*BT
      BBT13 = 0.25*XJ(3,3)*BT      
      
      BB(3,3*I-2) = XJI(3,3)*XJI(3,3)*BBT11  
      BB(3,3*I-1) = XJI(3,3)*XJI(3,3)*BBT12  
      BB(3,3*I-0) = XJI(3,3)*XJI(3,3)*BBT13 	
 	
      ENDDO
C
      RETURN 
	END
C*************************************************************************
C*************************************************************************
      SUBROUTINE SIMPANS(PAT,PAB,PBT,PBB,PCT,PCB,PDT,PDB,PE,PF,PG,PH,NODEX,NNO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -------------------------------------------------- 
C	ANS
	DIMENSION PAT(3,8),PAB(3,8),PBT(3,8),PBB(3,8),PCT(3,8),PCB(3,8),PDT(3,8),PDB(3,8)
	DIMENSION PE(3,8),PF(3,8),PG(3,8),PH(3,8),NODEX(1)
C     -------------------------------------------------- 
C     (1) TRANSVERSE SHEAR STRAINS
C     --------------------------------------------------       
C     AT - POINT
      RI =  0.0
      SI =  1.0
      TI =  1.0
      CALL SHAP3D (RI,SI,TI,H,PAT,NODEX,NNO)
C     AB - POINT
      RI =  0.0
      SI =  1.0
      TI = -1.0
      CALL SHAP3D (RI,SI,TI,H,PAB,NODEX,NNO)
C     BT - POINT
      RI = -1.0
      SI =  0.0
      TI =  1.0
      CALL SHAP3D (RI,SI,TI,H,PBT,NODEX,NNO)
C     BB - POINT
      RI = -1.0
      SI =  0.0
      TI = -1.0
      CALL SHAP3D (RI,SI,TI,H,PBB,NODEX,NNO)
C     CT - POINT
      RI =  0.0
      SI = -1.0
      TI =  1.0
      CALL SHAP3D (RI,SI,TI,H,PCT,NODEX,NNO)
C     CB - POINT
      RI =  0.0
      SI = -1.0
      TI = -1.0
      CALL SHAP3D (RI,SI,TI,H,PCB,NODEX,NNO)
C     DT - POINT
      RI =  1.0
      SI =  0.0
      TI =  1.0
      CALL SHAP3D (RI,SI,TI,H,PDT,NODEX,NNO)  
C     DB - POINT
      RI =  1.0
      SI =  0.0
      TI = -1.0
      CALL SHAP3D (RI,SI,TI,H,PDB,NODEX,NNO)   
      
C     -------------------------------------------------- 
C     (2) NORMAL TRANSVERSE SHEAR STRAINS
C     --------------------------------------------------  
C     E - POINT
      RI = -1.0
      SI = -1.0
      TI =  0.0
      CALL SHAP3D (RI,SI,TI,H,PE,NODEX,NNO)
C     F - POINT
      RI =  1.0
      SI = -1.0
      TI =  0.0
      CALL SHAP3D (RI,SI,TI,H,PF,NODEX,NNO)
C     G - POINT
      RI =  1.0
      SI =  1.0
      TI =  0.0
      CALL SHAP3D (RI,SI,TI,H,PG,NODEX,NNO)  
C     H - POINT
      RI = -1.0
      SI =  1.0
      TI =  0.0
      CALL SHAP3D (RI,SI,TI,H,PH,NODEX,NNO)               
     
	RETURN
	END

C*************************************************************************
      SUBROUTINE SOLSIGCO(STRAIN,STRESS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------------
C      FIND LOCAL STRESS USING PLANE STRAIN MATERIAL LAW
C     --------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
      DIMENSION STRAIN(6),STRESS(6)
      SMU=YM/(2*(1+PR))
      SLA=YM*PR/(1-PR*PR)
	STRESS(1)=(SLA+2*SMU)*STRAIN(1)+SLA*STRAIN(2)
      STRESS(2)=SLA*STRAIN(1)+(SLA+2*SMU)*STRAIN(2)
	STRESS(3)=YM*STRAIN(3)
	STRESS(4)=SMU*STRAIN(4)
	STRESS(5)=SMU*STRAIN(5)
	STRESS(6)=SMU*STRAIN(6)
	RETURN
	END


C*************************************************************************
      SUBROUTINE SOK0NLCO(S,BB,D,DVOL) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------
C     ADDS LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX
C      PLANE STRESS/STRAIN MATERIAL LAW
C	----------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     STIFF(24,24)    = ELEMENT STIFFNESS MATRIX 
C     B(6,24)  =  STRAIN-DISPLACEMENT MATRIX     
C     DVOL      = INTEGRATION CONSTANT =WT*DET
C     T(24,24): UL={T}UG
C     ----------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION STIFF(24,24),BB(6,24),S(300)
	DIMENSION D(6,6)
	
C     PLANE STRESS      
      NMAT = 1
      IF(NMAT.EQ.1) THEN
      D=0.
      SMU=YM/(2*(1+PR))
      SLA=YM*PR/(1-PR*PR)

	D(1,1)=SLA+2*SMU
	D(1,2)=SLA
	
      D(2,2)=SLA+2*SMU
	
	D(2,1)=SLA
	
	D(3,3)=YM
	D(4,4)=SMU
	D(5,5)=SMU
	D(6,6)=SMU
	
C     3D MATERIAL	
      ELSEIF(NMAT.EQ.2) THEN 
      D=0.
      SMU=YM*(1-PR)/((1+PR)*(1-2*PR))      

	D(1,1)=SMU*1
	D(1,2)=SMU*PR/(1-PR)
	D(1,3)=SMU*PR/(1-PR)
	
	D(2,1)=SMU*PR/(1-PR)
      D(2,2)=SMU*1
      D(2,3)=SMU*PR/(1-PR)
	
	D(3,1)=SMU*PR/(1-PR)
      D(3,2)=SMU*PR/(1-PR)
      D(3,3)=SMU*1	
	
	D(4,4)=SMU*(1-2*PR)/(2*(1-PR))
	D(5,5)=SMU*(1-2*PR)/(2*(1-PR))
	D(6,6)=SMU*(1-2*PR)/(2*(1-PR))
	ENDIF

C
      STIFF=MATMUL(MATMUL(TRANSPOSE(BB),D),BB)*DVOL
	KS=0                                 
      DO I=1,24                    
	   DO J=I,24                           
	      KS=KS+1 		                     
              S(KS)=STIFF(I,J)+S(KS)                 
         END DO
      END DO 
C
      RETURN
	END
c***************************************************************
        SUBROUTINE SHDSAMP(PSA,PTA,PRB,PTB,PSC,PTC,PRD,PTD,
     @                   PTA1,PTA2,PTA3,PTA4,RL,SL)
C     --------------------------------------
C	COMPUTE THE SHAPE FUNCTION DERIVATIVES
C     AT THE SAMPLING POINTS
C     --------------------------------------
C     PTA(8),PSA(8) = SHAPE FUNCTION DERIVATIVES AT A(-1,0,0)
C                     WITH RESPECT TO T,S RESPECTIVELY  
C     PTB(8),PRB(8) = SHAPE FUNCTION DERIVATIVES AT B(0,-1,0)
C                     WITH RESPECT TO T,R RESPECTIVELY
C     PTC(8),PSC(8) = SHAPE FUNCTION DERIVATIVES AT C(1,0,0)
C                     WITH RESPECT TO T,S RESPECTIVELY
C     PTD(8),PRD(8) = SHAPE FUNCTION DERIVATIVE AT D(0,1,0)
C                     WITH RESPECT TO T,R RESPECTIVELY
C     -------------------------------------------------------
C     PTA1(8),PTA2(8),PTA3(8),PTA4(8)
C            =SHAPE FUNCTION DERIVATIVES WITH RESPECT TO T AT                                      
C             A1(-1,-1,0),A2(1,-1,0),A3(1,1,0),A4(-1,1,0) RESPECTIVELY
C     RL(4),SL(4) = NODAL COORDINATES AT A1--A4
C     ----------------------------------------------------------------
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      DIMENSION PTA(8),PSA(8),PTB(8),PRB(8)
	DIMENSION PTC(8),PSC(8),PTD(8),PRD(8)
	DIMENSION PTA1(8),PTA2(8),PTA3(8),PTA4(8)
	DIMENSION RL(4),SL(4)
C
      PTA(1)=0.
	PTA(2)=.25
	PTA(3)=.25
	PTA(4)=0.
	PTA(5)=0.
	PTA(6)=-.25
	PTA(7)=-.25
	PTA(8)=0.
C
      PSA(1)=0.
	PSA(2)=.25
	PSA(3)=-.25
	PSA(4)=0.
	PSA(5)=0.
	PSA(6)=.25
	PSA(7)=-.25
	PSA(8)=0.
C
      PTB(1)=0.
	PTB(2)=0.
	PTB(3)=.25
	PTB(4)=.25
	PTB(5)=0.
	PTB(6)=0.
	PTB(7)=-.25
	PTB(8)=-.25
C
      PRB(1)=0.
	PRB(2)=0.
	PRB(3)=-.25
	PRB(4)=.25
	PRB(5)=0.
	PRB(6)=0.
	PRB(7)=-.25
	PRB(8)=.25
C
      PTC(1)=.25
	PTC(2)=0.
	PTC(3)=0.
	PTC(4)=.25
	PTC(5)=-.25
	PTC(6)=0.
	PTC(7)=0.
	PTC(8)=-.25
C
      PSC(1)=.25
	PSC(2)=0.
	PSC(3)=0.
	PSC(4)=-.25
	PSC(5)=.25
	PSC(6)=0.
	PSC(7)=0.
	PSC(8)=-.25
C
      PTD(1)=.25
	PTD(2)=.25
	PTD(3)=0.
	PTD(4)=0.
	PTD(5)=-.25
	PTD(6)=-.25
	PTD(7)=0.
	PTD(8)=0.
C
      PRD(1)=.25
	PRD(2)=-.25
	PRD(3)=0.
	PRD(4)=0.
	PRD(5)=.25
	PRD(6)=-.25
	PRD(7)=0.
	PRD(8)=0.
C    
      RL(1)=-1.
	RL(2)=1.
	RL(3)=1.
	RL(4)=-1.
C
      SL(1)=-1.
	SL(2)=-1.
	SL(3)=1.
	SL(4)=1.
C
      PTA1(1)=0.
      PTA1(2)=0.
	PTA1(3)=.5
	PTA1(4)=0.
	PTA1(5)=0.
	PTA1(6)=0.
	PTA1(7)=-.5
	PTA1(8)=0.
C
      PTA2(1)=0.
      PTA2(2)=0.
	PTA2(3)=0.
	PTA2(4)=.5
	PTA2(5)=0.
	PTA2(6)=0.
	PTA2(7)=0.
	PTA2(8)=-.5
C
      PTA3(1)=.5
      PTA3(2)=0.
	PTA3(3)=0.
	PTA3(4)=0.
	PTA3(5)=-.5
	PTA3(6)=0.
	PTA3(7)=0.
	PTA3(8)=0.
C
      PTA4(1)=0.
      PTA4(2)=.5
	PTA4(3)=0.
	PTA4(4)=0.
	PTA4(5)=0.
	PTA4(6)=-.5
	PTA4(7)=0.
	PTA4(8)=0.
C
      RETURN
	END 
	
C     ----------------------------------------------------
      SUBROUTINE SOK0NLIU(S,BB,DVOL) 
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------
C     ADDS LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX
C     FOR ONE-POINT INTEGRATION SCHEME
C	----------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     STIFF(24,24)    = ELEMENT STIFFNESS MATRIX 
C     B(6,24)  =  STRAIN-DISPLACEMENT MATRIX     
C     DVOL      = INTEGRATION CONSTANT =WT*DET
C     ----------------------------------------------------
      COMMON /HOOK/  A1,B1,C1,D1,A2,B2,C2,D2,BM,YM,PR,TH,YLD,ISR,IST
C
      DIMENSION STIFF(24,24),BB(6,24),S(300)
	DIMENSION D(6,6)
C
      D=0
	D(1,1)=A2
	D(1,2)=B2
	D(1,3)=B2
      D(2,2)=A2
	D(2,3)=B2
	D(2,1)=B2
	D(3,1)=B2
	D(3,2)=B2
	D(3,3)=A2
	D(4,4)=0.5
	D(5,5)=0.5
	D(6,6)=0.5
	D=C2*D

C
      STIFF=MATMUL(MATMUL(TRANSPOSE(BB),D),BB)*DVOL
	KS=0                                 
      DO I=1,24                    
	   DO J=I,24                           
	      KS=KS+1 		                     
              S(KS)=STIFF(I,J)+S(KS)                 
         END DO
      END DO 
C   
      RETURN
      END	
C*******************************************************************************
      SUBROUTINE STIFFTRANS(VR,VS,VT,S,RE)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C      TRANSFORM THE STIFFNESS MATRIX IN LOCAL COORDINATE INTO GLOBAL
C      S=TRANSPOSE(T)*S*T
C     -----------------------------------------------------------------
      COMMON /ELEM/  NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1               NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      DIMENSION VR(3),VS(3),VT(3),S(300)
      DIMENSION T(3,3),SK(3,3),RE(24)
	DIMENSION STIFF(24,24),TG(24,24),R(24),KUA(24,24)
C
  !    DO I=1,3
!	T(I,1)=VR(I)
!	T(I,2)=VS(I)
!	T(I,3)=VT(I)
!	ENDDO
	TG = 0.0D0
	DO I=1,8
      DO J=0,2
      TG(3*I-2,3*I-J)=VR(3-J)
	TG(3*I-1,3*I-J)=VS(3-J)
	TG(3*I,3*I-J)  =VT(3-J)
      ENDDO
      ENDDO
      
C     K = TGT.K.TG  ------------------------------------------       
      KS=0
	DO I=1,24
	DO J=I,24
	KS=KS+1
	STIFF(I,J)=S(KS)
	ENDDO 
	ENDDO
	DO I=1,24
	DO J=1,I-1
	STIFF(I,J)=STIFF(J,I)
	ENDDO
	ENDDO

       STIFF=MATMUL(MATMUL(TRANSPOSE(TG),STIFF),TG)
	KS=0
	DO I=1,24
	DO J=I,24
	KS=KS+1
	S(KS)=STIFF(I,J)
	ENDDO
	ENDDO

C     RE = TGT.RE  ------------------------------------------  
    !  IF(NLOPT.LE.1)RETURN
       RE=MATMUL(TRANSPOSE(TG),RE)
	RETURN
C
      RETURN
	END
C=====================================================================
C*******************************************************************************
      SUBROUTINE STRESSTRANS(VR,VS,VT,STRESS,GSTRESS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C      TRANSFORM THE STIFFNESS MATRIX IN LOCAL COORDINATE INTO GLOBAL
C      S=TRANSPOSE(T)*S*T
C     -----------------------------------------------------------------     
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION TT(3,3),STRESS(6),GSTRESS(6),STL(3,3),STG(3,3)	
C
      DO I=1,3
	TT(I,1) = VR(I)
	TT(I,2) = VS(I)
	TT(I,3) = VT(I)
	ENDDO	
C     
      STL(1,1) = STRESS(1)
      STL(1,2) = STRESS(4)
      STL(1,3) = STRESS(5)
      STL(2,1) = STRESS(4)
      STL(2,2) = STRESS(2)
      STL(2,3) = STRESS(6)
      STL(3,1) = STRESS(5)
      STL(3,2) = STRESS(6)
      STL(3,3) = STRESS(3)
      
      STG = MATMUL(MATMUL(TT,STL),TRANSPOSE(TT))   
    
      GSTRESS(1) = STG(1,1)
      GSTRESS(2) = STG(2,2) 
      GSTRESS(3) = STG(3,3)    
      GSTRESS(4) = STG(1,2)
      GSTRESS(5) = STG(1,3)      
      GSTRESS(6) = STG(2,3)           

      RETURN
      END
C	=====================================================================

C*******************************************************************************
      SUBROUTINE STRESSTRANS_LIU(VR,VS,VT,STRESS,GSTRESS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C      TRANSFORM THE STIFFNESS MATRIX IN LOCAL COORDINATE INTO GLOBAL
C      S=TRANSPOSE(T)*S*T
C     -----------------------------------------------------------------     
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION TT(3,3),TT2(3,3),STRESS(9),GSTRESS(9),STL(3,3),STL2(3,3),STG(3,3),STG2(3,3)      
C
      TT(1:3,1:3) = 0.0     
      TT2(1:3,1:3) = 0.0     

      DO I=1,3
	TT(I,1) = VR(I)
	TT(I,2) = VS(I)
	TT(I,3) = VT(I)
     
      ENDDO	
      
C     
      STL(1,1) = STRESS(1)
      STL(1,2) = STRESS(4)
      STL(1,3) = STRESS(5)
      STL(2,1) = STRESS(4)
      STL(2,2) = STRESS(2)
      STL(2,3) = STRESS(6)
      STL(3,1) = STRESS(5)
      STL(3,2) = STRESS(6)
      STL(3,3) = STRESS(3)
      
      STL2(1,1) = STRESS(7)
      STL2(2,2) = STRESS(8)
      STL2(3,3) = STRESS(9)
      
      STG = MATMUL(MATMUL(TT,STL),TRANSPOSE(TT))   
      STG2 = MATMUL(MATMUL(TT,STL2),TRANSPOSE(TT))
    
      GSTRESS(1) = STG(1,1)
      GSTRESS(2) = STG(2,2) 
      GSTRESS(3) = STG(3,3)    
      GSTRESS(4) = STG(1,2)
      GSTRESS(5) = STG(1,3)      
      GSTRESS(6) = STG(2,3)    
      GSTRESS(7) = STG2(1,1)
      GSTRESS(8) = STG2(2,2)      
      GSTRESS(9) = STG2(3,3)         

      RETURN
      END
C	=====================================================================
      SUBROUTINE STRESSLOCALTRANS(VR,VS,VT,STRESS,STRESSLO)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C      TRANSFORM THE STIFFNESS MATRIX IN LOCAL COORDINATE INTO GLOBAL
C      S=TRANSPOSE(T)*S*T
C     -----------------------------------------------------------------     
      DIMENSION VR(3),VS(3),VT(3)
      DIMENSION TT(3,3),STRESS(6),STRESSLO(6),STL(3,3),STG(3,3)	
C
      
      DO I=1,3
	TT(I,1) = VR(I)
	TT(I,2) = VS(I)
	TT(I,3) = VT(I)
	ENDDO	
C     
   
      STL(1,1) = STRESS(1)
      STL(1,2) = STRESS(4)
      STL(1,3) = STRESS(5)
      STL(2,1) = STRESS(4)
      STL(2,2) = STRESS(2)
      STL(2,3) = STRESS(6)
      STL(3,1) = STRESS(5)
      STL(3,2) = STRESS(6)
      STL(3,3) = STRESS(3)
      
      
      STG = MATMUL(MATMUL(TRANSPOSE(TT),STL),TT)   

    
      STRESSLO(1) = STG(1,1)
      STRESSLO(2) = STG(2,2) 
      STRESSLO(3) = STG(3,3)    
      STRESSLO(4) = STG(1,2)
      STRESSLO(5) = STG(1,3)      
      STRESSLO(6) = STG(2,3)           

      RETURN
	END
C	=====================================================================


      SUBROUTINE SOMDSP(COORD,COORDI,EDIS,REDIS,H,P,VR,VS,VT,NNO)
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
      DIMENSION XYZ(3),XYZO(3),CD(3,NNO),CDO(3,NNO)
      DIMENSION H(1),P(3,1),VR(3),VS(3),VT(3)
      DIMENSION VRO(3),VSO(3),VTO(3),TM(3,3)
C
      ARC=6.2831853071796
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
  
  	CALL SOLIRST(COORDI,P,VRO,VSO,VTO)
      
      DO 160 I=1,3
      DO 160 J=1,3
 160  TM(I,J)=VR(I)*VRO(J)+VS(I)*VSO(J)+VT(I)*VTO(J)
 
      ABCD=(.5*(TM(1,1)+TM(2,2)+TM(3,3)-1.))
	IF (ABCD .GT. 1.00000000000000000 )ABCD=1.000000000000000000
	IF (ABCD .LT. -1.00000000000000000 )ABCD=-1.000000000000000000
      RR=ACOS(ABCD)
      SN=SIN(RR)
 
      K=1
      DO 210 I=1,NNO
      DO 210 J=1,3
      TCDO=TM(J,1)*CDO(1,I)+TM(J,2)*CDO(2,I)+TM(J,3)*CDO(3,I)
      REDIS(K)=CD(J,I)-TCDO
  210 K=K+1
  
      RETURN
      END
C
C	=====================================================================
      SUBROUTINE STIFFTRANE(VR,VS,VT,MM,SEL)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     -----------------------------------------------------------------
C      TRANSFORM THE STIFFNESS MATRIX IN LOCAL COORDINATE INTO GLOBAL
C      SEL=SEL*T
C     -----------------------------------------------------------------

      DIMENSION VR(3),VS(3),VT(3),SEL(MM,24),SELO(MM,24)
	DIMENSION TG(24,24)

	TG = 0.0D0
	DO I=1,8
      DO J=0,2
      TG(3*I-2,3*I-J)=VR(3-J)
	TG(3*I-1,3*I-J)=VS(3-J)
	TG(3*I,3*I-J)  =VT(3-J)
      ENDDO
      ENDDO
      
      SELO = SEL
      SEL=MATMUL(SELO,TG)

	RETURN
C
      RETURN
	END
C=====================================================================

	
C=====================================================================
      SUBROUTINE SOLI632(PROPM,PROPG,NODEX,WA,S,COORD,EDIS,EDISI,RE,MWG
     +                  ,ALPHA,SEL,SEDI,FIN,HINFC)
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
C
      PARAMETER (MM=21)
C
      DIMENSION PROPM(*),PROPG(*),NODEX(*),WA(MWG,1),S(*),COORD(*)
      DIMENSION EDIS(*),EDISI(*),RE(*)
      DIMENSION DP(6,6),H(21),P(3,21),XJI(9),B(63),DISD(9)
      DIMENSION STRAIN(6),QSTRAI(6),STRESS(6),TAU(6)
C     ------------
C     EAS METHOD
C     ------------    
	DIMENSION XJ(3,3),TTO(6,6),TT(6,6),XJO(3,3),TM(6,MM)
	DIMENSION EAS(6,1),ALPHA(MM),TMT(MM,6),RE1(24),RE2(24,1)
      DIMENSION SED(MM,MM),SEL(MM,1),SEDI(*),RHH(MM,1),RH(MM,1)
C	NEXT ADDED LINE BY GILSON - JUL2003 (INT FORCE)
	DIMENSION FIN(NEF)

	DIMENSION HINFC(MM)
C     ------------------------------------------------------------
C     SET VALUES FOR LINEAR STRESS-STRAIN LAW (COMMON BLOCK /HOOK/
C     INITIALISATION OF INTEGRATION RULE
C     ------------------------------------------------------------
      CALL HOKLAW (PROPM,PROPG,1)
C     --------------
C     COMPUTE JACOBIAN MATRIX AT THE ORIGIN (R=0.0,S=0.0,T=0.0)
C     --------------
      CALL SHAP3D (0.0D0,0.0D0,0.0D0,H,P,NODEX,NNO)
	CALL JACO3D(COORD,P,XJO,XJI,DETO,MEL,NNO)
C     ----------------
C     COMPUTE T MATRIX (FOR EAS METHOD)
C     ---------------
	CALL MATRIXT(XJO,TTO)
      CALL INVMAT1(TTO,TT,6)
	DO 5 I=1,6
	DO 5 J=1,6
	TTO(I,J)=TT(J,I)
  5   CONTINUE

      MGR = NGR
      MGS = NGS
      MGT = NGT
      IF (ITASK.NE.5) GOTO 10
      MGR = 3
      MGS = 3
      MGT = 3
C     ----------------------
C     LOOP OVER GAUSS POINTS
C     ----------------------
 10   IPT = 0
C     ---------------
C     INITIALIZATION OF SOME VARIABLE MATRICES
C     ---------------
      DO 11 I=1,MM
	DO 12 J=1,MM
	SED(I,J)=0.0
 12   CONTINUE
      DO 13 K=1,24
	SEL(I,K)=0.0
 13   CONTINUE
      RH(I,1)=0.0
 11   CONTINUE
      DO 17 I=1,24
	RE1(I)=0.0
 17   CONTINUE

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
      CALL SHAP3D (RI,SI,TI,H,P,NODEX,NNO)
      CALL JACO3D (COORD,P,XJ,XJI,DET,MEL,NNO)
      DVOL = WT*DET
C     -----------------------------------------
C     ADD CONTRIBUTION TO MASS MATRIX (ITASK=5)
C     -----------------------------------------
      IF (ITASK.NE.5)  GOTO 50
      CALL SOMASS (S,H,PROPM(5),DVOL,NNO,NEF,IPT)
      GOTO 900
C     ----------------------------------------
C     COMPACTED STRAIN-DISPLACEMENT MATRIX (B)
C     ----------------------------------------
 50   CALL SOBMAT (P,XJI,B,NNO)
C     ----------------     
C     FOR ENHANCED STRAIN METHOD 
C     COMPUTE MATRIX M AND ITS TRANSPOSE MATRIX
C     ----------------
      CALL MATRIXTM(RI,SI,TI,TTO,TM,MM)
      DO 51 I=1,6
	DO 51 J=1,MM
	TM(I,J)=TM(I,J)*DETO/DET
	TMT(J,I)=TM(I,J)
 51   CONTINUE
C     ------------------
C     D=TRANSPOSE(M)*E*M
C     L=TRANSPOSE(M)*E*B
C     ------------------
      CALL MATRIXD(TM,SED,DVOL,MM)
      CALL MATRIXL(TM,SEL,B,DVOL,MM)
      IF (NLOPT+ITASK.EQ.1)  GOTO 700
C     ----------------------------------------------------------------
C     FIND STRESSES AND CALCULATE GEOMETRIC STIFFNESS MATRIX (ITASK=4)
C     ----------------------------------------------------------------
      IF (ITASK.NE.4)  GOTO 200
      DO 100  I=1,6
 100  TAU(I) = WA(I,IPT)*DVOL
      GOTO 800
C     -------------------------------
C     DISPLACEMENT DERIVATIVES (DISD)
C     -------------------------------
 200  IF (ITASK.GT.1) THEN
      CALL CLEARA (DISD,9)
      DO 290  IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      DISD(1) = DISD(1) + B(IEF)*EDIS(IEF)
      DISD(2) = DISD(2) + B(JEF)*EDIS(JEF)
      DISD(3) = DISD(3) + B(KEF)*EDIS(KEF)
      DISD(4) = DISD(4) + B(JEF)*EDIS(IEF)
      DISD(5) = DISD(5) + B(KEF)*EDIS(IEF)
      DISD(6) = DISD(6) + B(IEF)*EDIS(JEF)
      DISD(7) = DISD(7) + B(KEF)*EDIS(JEF)
      DISD(8) = DISD(8) + B(IEF)*EDIS(KEF)
 290  DISD(9) = DISD(9) + B(JEF)*EDIS(KEF)
C     -------------------
C     LINEAR STRAIN TERMS
C     -------------------
      STRAIN(1) = DISD(1)
      STRAIN(2) = DISD(2)
      STRAIN(3) = DISD(3)
      STRAIN(4) = DISD(4) + DISD(6)
      STRAIN(5) = DISD(5) + DISD(8)
      STRAIN(6) = DISD(7) + DISD(9)
C     --------
C     FOR EAS TERM (EPS(EAS)=M*ALPHA)
C     --------
      CALL MATRIX2MUL(TM,ALPHA,EAS,6,MM,1)
	DO 395 I=1,6
      STRAIN(I)=STRAIN(I) - EAS(I,1)
 395  CONTINUE
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
      DO 390  I=1,6
 390  STRAIN(I) = STRAIN(I) - QSTRAI(I)
C     ------------------------------------
C     COMPUTE AND STORE NONLINEAR STRESSES
C     ------------------------------------
 400  IF (MTMOD-3)  405,450,460
 405  CALL SOLSIG (STRAIN,STRESS)
      DO 410  I=1,6
 410  WA(I,IPT) = STRESS(I)
      GOTO 500
C
 450  CALL MISE3D (WA(1,IPT),WA(7,IPT),WA(13,IPT),WA(13,IPT),
     1             STRAIN,STRESS,DP)
      GOTO 500
C
C460  CALL CONC3D -- DELETED BECAUSE OF SUBROUTINE NOT EXIST
 460  CONTINUE
C     -----------------------------------------------
C     ADD STRESS CONTRIBUTION TO ELEMENT FORCE VECTOR
C     -----------------------------------------------
 500  DO 510  I=1,6
      TAU(I) = STRESS(I)*DVOL
 510  CONTINUE
      ENDIF
C      IF (ITASK.LE.2) GOTO 520
      IF (ITASK.LE.3) GOTO 520
      IF (IFEIG.EQ.0) GOTO 800
      GOTO 900
c 520  IF (ITASK.GT.1) THEN !blocked - jan2005
C     ---------
C     EQUILIBRIUM FORCE OF COMPATIBLE ELEMENT
C     ---------
520	DO 550  IEF=1,NEF,3
      JEF = IEF+1
      KEF = IEF+2
      RE1(IEF)=RE1(IEF) + B(IEF)*TAU(1) + B(JEF)*TAU(4) + B(KEF)*TAU(5)
      RE1(JEF)=RE1(JEF) + B(JEF)*TAU(2) + B(IEF)*TAU(4) + B(KEF)*TAU(6)
      RE1(KEF)=RE1(KEF) + B(KEF)*TAU(3) + B(IEF)*TAU(5) + B(JEF)*TAU(6)
 550  CONTINUE
C     ---------
C     COMPUTE EQUILIBRIUM FORCE FOR EAS METHOD
C     ---------
      CALL MATRIX2MUL(TMT,TAU,RHH,MM,6,1)
	DO 552 I=1,MM
	RH(I,1)=RH(I,1)+RHH(I,1)
 552  CONTINUE
c      ENDIF
C     -------------------------------------------------------------
C     FOR STIFFNESS REFORMATION ONLY (IFREF=0)
C     ADD CONTRIBUTIONS OF INTEGRATED [B]T*[B] INTO [S]   (MTMOD<2)
C     ADD LINEAR CONTRIBUTION TO ELEMENT STIFFNESS MATRIX (MTMOD>2)
C     -------------------------------------------------------------
      IF (IFREF) 900,700,900
 700  IF (NLOPT.NE.0) GOTO 750
	CALL SOK0NL (S,B,NEF,DVOL)
	GOTO 790
C
 750  IF (MTMOD.LE.2) CALL SOK0NL (S,B,NEF,DVOL)
      IF (MTMOD.GT.2) CALL SOK0PL (S,DP,B,NEF,DVOL)
C     -------------------------------------------------------
C     ADD NONLINEAR CONTRIBUTION TO STIFFNESS MATRIX (NLOPT>2)
C     --------------------------------------------------------
 790  IF (NLOPT.LE.1) GOTO 810
 800  CALL SOKSIG (S,TAU,B,NEF)
 810  CONTINUE !TIM(12) = TIM(12) + TIM2-TIM1
C
 900  CONTINUE
C     --------
C     FOR NONLINEAR ANALYSIS OF EAS METHOD
C     --------

C	STORE THE EAS FORCE SONGSAK MAR2006
	DO I = 1,MM
	HINFC(I) = RH(I,1)
	ENDDO


  	CALL INVMAT1(SED,SEDI,MM)
      CALL MATRIX3MUL(SEL,SEDI,RH,RE2,24,MM,1)
      IF (ITASK.GT.1) THEN
	DO 910 I=1,24
	RE(I)=RE1(I) + RE2(I,1)
 910  CONTINUE
      ENDIF

C	NEXT BLOCK ADDED BY GILSON - JUL2003 (INT FORCE)
	IF (ITASK.EQ.3) THEN
	  DO 2000 I = 1,NEF
	    FIN(I) = RE(I)
2000	  CONTINUE
	ENDIF

C     ------------------
C     FOR ENHANCED STRAIN METHOD
C     K=TRANSPOSE(L)*INVERSE(D)*L
C     K(FINAL)=K(PURE DISP.)-K(EAS)
C     ------------------ 
      IF (ITASK.GT.1) RETURN
      CALL MATRIXES(SEL,SEDI,S,MM)
      RETURN
      END
C	==============================================================
      SUBROUTINE MATRIXTM(R,S,T,TT,EM,MM)
C	--------
C     FOR ENHANCED STRAIN TERM OF ENHANCED STRAIN METHOD
C     COMPUTE COEFFICIENT MATRIX EM(6,30)=T*M
C     IN MASTER OR LOCAL COORDINATE SYSTEM
C	--------
      IMPLICIT REAL*8 (A-H,O-Z)
	IMPLICIT INTEGER*4 (I-N)
	DIMENSION EM(6,MM),TT(6,6)
      DO 100 I=1,6
	EM(I,1)=TT(I,1)*R
	EM(I,2)=TT(I,2)*S
	EM(I,3)=TT(I,3)*T
	EM(I,4)=TT(I,4)*R
	EM(I,5)=TT(I,4)*S
	EM(I,6)=TT(I,5)*R
	EM(I,7)=TT(I,5)*T
	EM(I,8)=TT(I,6)*S
	EM(I,9)=TT(I,6)*T

	IF (MM.EQ.9) GOTO 100
	EM(I,10)=TT(I,4)*R*T
	EM(I,11)=TT(I,4)*S*T
	EM(I,12)=TT(I,5)*R*S
	EM(I,13)=TT(I,5)*S*T
	EM(I,14)=TT(I,6)*R*S
	EM(I,15)=TT(I,6)*R*T
	
	IF (MM.EQ.15) GOTO 100
	EM(I,16)=TT(I,1)*R*S
	EM(I,17)=TT(I,1)*R*T
	EM(I,18)=TT(I,2)*R*S
	EM(I,19)=TT(I,2)*S*T
	EM(I,20)=TT(I,3)*R*T
	EM(I,21)=TT(I,3)*S*T
	
      IF (MM.EQ.21) GOTO 100
	EM(I,22)=TT(I,4)*R*S
	EM(I,23)=TT(I,5)*R*T
	EM(I,24)=TT(I,6)*S*T
	EM(I,25)=TT(I,1)*R*S*T
	EM(I,26)=TT(I,2)*R*S*T
	EM(I,27)=TT(I,3)*R*S*T
	EM(I,28)=TT(I,4)*R*S*T
	EM(I,29)=TT(I,5)*R*S*T
	EM(I,30)=TT(I,6)*R*S*T

 100	CONTINUE
      RETURN
      END
