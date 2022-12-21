C=====================================================================
C	START OF SUBROUTINE FOR LANCZOS EIGEN SOLVER BY SUNIL
C=====================================================================
C
C=====================================================================
      SUBROUTINE LANC (W,ID,MAXA,N11,N10,AA,BB,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP1,TYP2
C	09/02/15 SUNIL 
C     ----------------------------------------------------------------
C     SOLUTION OF STANDARD AND GENERALIZED EIGENVALUE PROBLEMS BY
C	SUB-SPACE INTERATIONS INCLUDING SHIFTING OF NON +VE DEFINIT [K]
c	-----------
C	The pointer N8,N10 added Jan18,2003 by NguyenDV for Mode super
C	(The pointer N8 changed to N11, Oct11,2003 by NguyenDV)
C     ----------------------------------------------------------------

      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C
      COMMON A(9000000),IA(9000000)
C	--------------------------------------------------------------------
      DIMENSION W(1),ID(NSF,1),MAXA(1),
     +            RVEC(NEQ,NROOT),REIGV(NROOT)
	DIMENSION AA(1),BB(1)
C	
C	--------------
C     INITIALISATION
C     --------------
      IF (NROOT.LE.0) NROOT = 1
      IF (NROOT.GT.NEQ) NROOT = NEQ
	IF (IVPRT.LT.0) IVPRT = 0
	IF (IVPRT.GT.NROOT) IVPRT = NROOT
      NEIG = NEIG+1
      IF (KSTEP-1.GT.NSEIG) GOTO 100
      IF (IEIG.EQ.1) WRITE (ISO,3100)
      IF (IEIG.EQ.2) WRITE (ISO,3200)
      IF (IEIG.EQ.3) WRITE (ISO,3300)
      WRITE (ISO,3500) NSEIG,NROOT,IVPRT,NC,NITEM,IFSS,EPS,SHIFT0
	GOTO 150

 100	WRITE (ISO,9999)
	STOP
C     -------------------------------------------------------------
C     READ TANGENTIAL STIF MATRIX (2) AND GEO. STIF/MASS MATRIX (3)
C     -------------------------------------------------------------
 150	IF (IEIG.NE.1) GOTO 500

C	DO 160  IEQ=1,NEQ
C 160	BB(IEQ) = 1.0
	FAC = 1.0D0
	CALL MDMAKE_DIA(TYP2,FAC,IA(LMA))  !MAKE DIAGONAL B MATRIX

C     -------------------------------------------------------------
C	FLIP STIFFNESS AND MASS COLUMN UPSIDE DOWN
C     -------------------------------------------------------------
 500	CALL CPU_TIME (TIM1)
C     ---------------------------------------------------------
C     SOLVE CURRENT EIGENVALUE PROBLEM BY LANCZOS VECTOR METHOD
C     ---------------------------------------------------------
	IF (ISOLV.EQ.2)
     1	CALL LANCZOS(ID,MAXA,W,NEQ,NROOT,IVPRT,NSF,N11,N10,
     2				 AA,BB,TYP1,TYP2)

C     ------------------------------------------------------------
C     SOLVE CURRENT EIGENVALUE PROBLEM BY INVERSE ITERATION METHOD
C     ------------------------------------------------------------
C	IF (ISOLV.EQ.3)
C     +	CALL INVERSIT(ID,JDIAG,AA,BB,W,NEQ,NROOT,IVPRT)
	CALL CPU_TIME (TIM2)

	TIM(18) = TIM(18) + (TIM2-TIM1)
      ITASK = 1
C
 3100 FORMAT(//1X,79(1HX)//14X,15HS T A N D A R D
     1    ,40H   E I G E N V A L U E   A N A L Y S I S/14X,55(1H-)/)
 3200 FORMAT(//1X,79(1HX)//14X,15HB U C K L I N G
     1    ,40H   E I G E N V A L U E   A N A L Y S I S/14X,55(1H-)/)
 3300 FORMAT(//1X,79(1HX)//14X,17HV I B R A T I O N
     1    ,38H   F R E Q U E N C Y   A N A L Y S I S/14X,55(1H-)/)
 3500 FORMAT(
     +    14X,'INTERVAL BETWEEN EIGENVALUE ANALYSIS . . NSEIG =',I7/
     +    14X,'NUMBER OF EIGENVALUES REQUIRED . . . . . NROOT =',I7/
     +    14X,'NUMBER OF EIGENVECTORS REQUIRED. . . . . IVPRT =',I7/
     +    14X,'NUMBER OF ITERATION VECTORS USED . . . . .  NC =',I7/
     +    14X,'NUMBER OF SUBSPACE ITERATIONS PERMITTED  NITEM =',I7/
     +    14X,'FLAG FOR STURM SEQUENCE CHECK  . . . . .  IFSS =',I7/
     +    14X,'CONVERGENCE TOLERANCE ON EIGENVALUES . . . EPS =',E12.4/
     +    14X,'INITIAL SHIFT TO BE APPLIED . . . . . . SHIFT0 =',E12.4/
     +    /1X,79(1HX))
9999	FORMAT(//1X,'THIS OPTION IS FOR NONLINEAR EIGHE VALUE ANALYSIS'/
     +         1X,'THIS IS OPTION IS NOT AVAILABLE NOW')
C
      RETURN
      END
C
C=====================================================================
C *****************************************************************
C
C     E I G E N S O L V E R   B Y   L A N C Z O S   M E T H O D
C    	---------------------------------------------------------
C          SUBROUTINE LANCZOS(JDIAG,A,B,NEQ,M,NR,DFL)
C		 --------------------------------------------
C                 INSTRUCTIONS OF PROGRAM LANCZOS
C     		    -------------------------------
C  1. GENERAL FEATURE OF PROGRAM LANCZOS
C
C  THIS SUBROUTINE IS DESIGNED TO SOLVE EIGENVALUES AND
C  EIGENVECTORS OF PROBLEMS IN THE FOLLOWING FORM:
C
C     (K)(X) = W**2 (M)(X)  ;WHERE (K) IS POSITIVE DEFINITE
C
C  A GIVEN NUMBER OF THE LOWEST EIGENSOLUTIONS CAN BE COMPUTED BY
C  A TECHNIQUE CALLED LANCZOS METHOD. THIS SUBROUTINE CAN HANDLE
C  BOTH LUMPED MASS MATRIX AND CONSISTEN MASS MATRIX. THE ALGORITHM
C  IMPLEMENTED HERE CLOSELY FOLLOWS THE ONE SUGGESTED BY GOLUB ET.
C  AL. (1972) AND FURTHER DEVELOPED BY CHOWDHURY (1975). (SEE
C  REFERENCES). THE PROGRAM LANCZOS WAS ORIGINALLY DEVELOPED BY
C  CHOWDHURY AND HAS BEEN MODIFIED BY KENJI KIMURA IN 1988 FOR 
C  THE M.ENG. THESIS AT ASIAN INSTITUTE OF TECHNOLOGY.
C
C  2. USERS MANUAL
C
C  LANCZOS CAN BE USED AS AN INDEPENDENT PROGRAM FOR SOLVING EIGEN
C  VALUE PROBLEMS BY BEING SUPPLIED WITH THE FOLLOWING SETS OF DATA:
C     A(I)     :SYMMETRIC POSITIVE DEFINITE STIFFNESS MATRIX STORED
C               IN A PROFILE FORM
C     JDIAG(I) :DIAGONAL POINTERS OF SYMMETRIC STIFFNESS MATRIX
C     B(I)     :LUMPED MASS MATRIX STORED IN A VECTOR FORM OR
C               CONSISTENT MASS MATRIX STORED IN A PROFILE FORM
C     NEQ      :ORDER OF ORIGINAL EIGENSYSTEM
C     M        :ORDER OF TRUNCATED (REDUCED) EIGENSYSTEM
C     NR       :NUMBER OF LOWEST EIGENVALUES TO BE SOLVED
C     DFL      :LOGICAL FLAG TO DISTINGUISH THE TYPE OF MASS MATRIX
C               LUMPED MASS     : DFL=.TRUE.
C               CONSISTENT MASS : DFL=.FALSE.
C
C  3. ORGANIZATION OF PROGRAM LANCZOS
C
C  THE PROGRAM LANCZOS CONSISTS OF 11 SUBROUTINES. THE FUNCTIONS OF
C  EACH SUBROUTINE AND BRIEF ALGORITHM ARE AS FOLLOWS:
C
C     LANCZOS
C     DETERMINATION OF PROBLEM CONFIGURATION
C     ALLOCATION OF WORK SPACE
C                             ]
C     CHOLS
C     CHOLESKY FACTORIZATION OF STIFFNESS MATRIX
C                             ]
C     GOLUB
C     CALCULATION OF LANCZOS VECTOR AND TRIDIAGONAL
C     MATRIX ELEMENTS OF REDUCED SYSTEM
C     (NLNVC AND MULPRO ARE CALLED IN GOLUB)
C                             ]
C     TDBS
C     CALCULATION OF EIGENVALUES OF REDUCED SYSTEM USING
C     BISECTION METHOD WITH STRUM SEQUENCE CHECK
C     (STURM IS CALLED IN TDBS)
C                             ]
C     TDII & ORTH
C     CALCULATION OF EIGENVECTORS OF REDUCED SYSTEM USING
C     INVERSE ITERATION THROUGH ORTHOGONALIZATION TECHNIQUE
C     TO ISOLATE EIGENVECTORS CORRESPONDING TO CLOSE SETS
C     OF EIGENVALUES (ORTH IS CALLED FOR ORTHOGONALIZATION)
C                             ]
C     BKTFM & BKSUB
C     CALCULATION OF EIGENVECTORS OF ORIGINAL SYSTEM
C     USING EIGENVECTORS OBTAINED IN TDII AND ORTH
C
C  4. SIZE OF SUBSPACE
C
C  ACCORDING TO THE NUMBER OF REQUIRED EIGENVALUES, THE SUBSPACE
C  SIZE M IS DETERMINED AS FOLLOWS (UNLESS USER DEFINES
C  THE SUBSPACE SIZE)
C      M = 4*NR  FOR NR >= 3
C      M = 10    FOR NR = 1 OR 2
C  USER CAN HAVE AN OPTION TO SET THE SUBSPACE SIZE ARBITRARILY.
C  IF NR AND M ARE INADEQUATE, THE PROGRAM WILL AUTOMATICALLY
C  CORRECT THE VALUES.     
C      IF NR =< 0   THEN NR = 5
C      IF NR >  N   THEN NR = N
C      IF M  < NR*4 THEN M  = 4*NR
C      IF M  >  N   THEN M  = N
C
C  5. MEMORY ESTIMATE
C
C  ALL THE ARRAYS USED IN THE PROGRAM EXCEPT UER-SUPPLIED-DATA
C  RESIDE IN A SINGLE ARRAY W(I). THE POINTERS INDICATING THE
C  LOCATION OF EACH INDIVIDUAL ARRAY ARE AS FOLLOWS:
C      N1  = 1+N        N7  = N6+M
C      N2  = N1+N       N8  = N7+M*(2*N-M+1)/2
C      N3  = N2+N       N9  = N8+NR
C      N4  = N3+N+2     N10 = N9+M+((NR-1)*NN+1)/2
C      N5  = N4+M       NSIZE = N10+NR*N
C      N6  = N5+M
C
C  WHERE     N  : ORDER OF ORIGINAL SYSTEM (N=NEQ)
C            M  : ORDER OF REDUCED SYSTEM
C            NR : NUMBER OF EIGENVALUES TO BE SOLVED
C            NN : MAX (N,2*M)
C            NSIZE : REQUIRED MEMORY SIZE
C
C  PRINCIPAL VARIABLES ARE STORED IN THE ARRAY W(I) FROM THE
C  FOLLOWING STARTING POINTS:
C     W(N4)  :ALFA(I) DIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX
C     W(N5)  :BETA(I) SEMIDIAGONAL ELEMENTS OF TRIDIAGONAL MATRIX
C     W(N7)  :U(I) VECTORS TO FORM HOUSEHOLDER REFLECTION MATRIX
C              DESCRIBED IN SUBROUTINE GOLUB
C     W(N8)  :REQUIRED EIGENVALUES
C     W(N9)  :TVEC(I) EIGENVECTORS OF REDUCED SYSTEM
C     W(N10) :EIGENVECTORS CORRESPONDING TO REQUIRED EIGENVALUES
C              FOR ORIGINAL SYSTEM
C
C  6. OTHER IMPORTANT PARAMETERS
C
C  OTHER IMPORTANT PARAMETERS USED IN THE PROGRAM ARE DESCRIBED AS
C  FOLLOWS:
C     MEM     : SIZE OF MEMORY CAPACITY
C               (IN CASE OF MEMORY SHORTAGE, INCREASE THIS VALUE)
C     RL      : LOWER LIMIT OF EIGENVALUE RANGE SET TO 1.0E-50
C               (SPECIFY HIGHEST PERMISSIBLE VALUE TO MINIMIZE
C                NUMBER OF BISECTIONS)
C     RU      : UPPER LIMIT OF EIGENVALUE RANGE SET TO 1.0E+70
C               (IN THE ABSENCE OF A RELIABLE UPPER-BOUND AT THE
C                FIRST STAGE OF COMPUTATION, SUPPLY A GROSS OVER
C                ESTIMATE WHICH WILL BE APPROPRIATELY REDUCED)
C     EPS     : ACCURACY TO WHICH EIGENVALUES ARE REQUIRED
C
C *****************************************************************
C
C      SUBROUTINE LANCZOS(ID,JDIAG,A,B,W,NEQ,NR,DFL,IVPRT,NF)
C      SUBROUTINE LANCZOS(ID,JDIAG,A,B,W,NEQ,NR,DFL,IVPRT,NF,N8,N10)
C
	SUBROUTINE LANCZOS(ID,MAXA,W,NEQ,NR,IVPRT,NF,N11,N10,
	1				   AA,BB,TYP1,TYP2)
      IMPLICIT REAL*8(A-H,O-Z)
	CHARACTER*4 TYP1,TYP2


	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI

      DIMENSION ID(NF,1),JDIAG(1),W(1),MAXA(1)
	DIMENSION AA(1),BB(1),DD(NEQ)

      DATA RL/-1.0D70/,RU/1.0D70/,EPS/1.0D-09/

C.... DETERMINATION OF M (SIZE OF SUBSPACE)

C	SUNIL
C	M = 3*NR
C      IF (NR.GT.25) M = 75+(NR-25)*2
C      IF (NR.LE.10) M = 40
C      IF (NR.LE.05) M = 25

C	SUNIL (SIZE OF SUBSPACE IS DETERMINED BY LINEAR INTERPOLATION
C	BETWEEN  M=50 FOR NR=5 AND M=100 FOR NR=50
	M = 2*NR
      IF (NR.LE.50) M = IFIX(50.0+(FLOAT(NR)-5.0)*50.0/45.0)+1
      IF (NR.LE.05) M = 50
      IF (M.GT.NEQ) M = NEQ

C.... ALLOCATION OF THE WORK SPACE
      N = NEQ
      NN = MAX0(N,2*M)

      N1 = 1 +N   
      N2 = N1+N   
      N3 = N2+N   
      N4 = N3+N+2   
      N5 = N4+M   
      N6 = N5+M   
      N7 = N6+M   
      N8 = N7+M*(2*N-M+1)/2   
      N9 = N8+NR   
      N10 = N9+M+((NR-1)*NN+1)/2   
      N11 = N10+NR*N   
C
	N12 = N11+NR
      NSIZE = N12

10	CALL CLEARA(W,NSIZE)


C	USING STANDARD FACTORIZE  K = L*D*LT  INSTEAD OF CHOLESKY FACTORIZATION  K = L*LT 
	INDPD = 0
	CALL COLSOL (MAXA,AA,DD,VV,1,INDPD,TYP1,'TEMP')  !FACTORIZE STIFFNESS MATRIX

      CALL GOLUB(AA,BB,N,M,W,W(N1),W(N2),W(N3),W(N4),W(N5),W(N7),
     1           MAXA,DD,TYP1,TYP2,'TEMP')

	CALL TDBS(W(N8),W(N4),W(N5),W(N6),RL,RU,NR,M,EPS,ACC,W(N11))


      IF(IVPRT.NE.0)CALL TDII(W(N9),W(N8),W(N4),W(N5),W(N6),W,
     +			    		W(N1),W(N2),W(N3),IVPRT,M,NN,EPS,ACC)
      IF(IVPRT.NE.0)CALL BKTFM(ID,N,M,IVPRT,W(N9),W(N10),W,W(N1),
     + 					     W(N7),NN,MAXA,AA,DD,'TEMP')


      RETURN
1000  FORMAT(/,5X,43H*******  WARNING : MEMORY SHORTAGE  *******,/,
     *       /,14X,18HMEMORY CAPACITY = ,I6,/,14X,
     *       18HMEMORY REQUIRED = ,I6)
      END
C
C=====================================================================
      SUBROUTINE GOLUB(AA,BB,N,M,Q,Q0,Q1,R,ALFA,BETA,U,MAXA,DD,
	1				 TYP1,TYP2,TYP3)
      IMPLICIT REAL*8(A-H,O-Z)
	CHARACTER*4 TYP1,TYP2,TYP3
      DIMENSION Q(N),Q0(N),Q1(N),R(1),ALFA(M),BETA(M)
      DIMENSION MAXA(1),U(1),DD(1)
C
C ****************************************************************
C * THIS SUBROUTINE CALCULATES THE LANCZOS VECTORS Q(K) AND THE **
C *  COMPONENTS OF THE TRIDIAGONAL MATRIX ALFA(K) AND BETA(K)   **
C *            BRIEF ALGORITHM IS AS FOLLOWS                    **
C *                                                             **
C *                 K TH LANCZOS VECTOR                         **
C *                          ]                                  **
C *                CALCULATION OF A*Q(K)                        **
C *                          ]                                  **
C *         CALCULATION OF ALFA(K) AND BETA(K-1)                **
C *                          ]                                  **
C *    CALCULATION OF K+1 TH HOUSEHOLDER REFLECTION MATRIX      **
C *                          ]                                  **
C *                       K = K+1                               **
C *                                                             **
C ****************************************************************
C
C     VARIABLES USED IN THIS SUBROUTINE
C
C     N    :ORDER OF THE ORIGINAL EIGENSYSTEM CONSIDERED
C     M    :ORDER OF THE REDUCED SYSTEM 
C           (I.E. OF DERIVED TRIANGULAR MATRIX)
C     Q    :K+1 TH LANCZOS VECTOR
C     Q0   :K-1 TH LANCZOS VECTOR
C     Q1   :K TH LANCZOS VECTOR
C     R    :R(K)=A*Q(K) CALCULATED BY SUBROUTINE NLNVC
C     ALFA :DIAGONAL COMPONENTS OF TRIDIAGONAL MATRIX
C     BETA :SEMIDIAGONAL COMPONENTS OF TRIDIAGONAL MATRIX
C     U    :VECTORS TO FORM HOUSEHOLDER REFLECTION MATRIX P
C           FOR REORTHOGONALIZATION
C     JDIAG:DIAGONAL POINTERS OF THE STIFFNESS MATRIX
C     A    :SYMMETRIC STIFFNESS MATRIX IN A PROFILE FORM
C     B    :LUMPED MASS MATRIX IN A VECTOR OR
C           CONSISTENT MASS MATRIX IN A PROFILE FORM
C     DFL  :LOGICAL FLAG TO DISTINGUISH THE TYPE OF MASS MATRIX
C
C     SUBROUTINE TO BE CALLED
C
C     NLNVC  :CALCULATION OF R(K)=L(-1)*M*L(-T)*Q(K)
C
C *****************************************************************
C
	DIMENSION AA(1),BB(1)

C.... INITIAL LANCZOS VECTOR
      Q(1) = 1.0D0
      DO 100 I=2,N
		Q(I) = 0.0
100	CONTINUE


      BETA(M) = 0.0D0
C.... CALCULATION OF ALFA(K), BETA(K-1) AND K+1 LANCZOS VECTOR
      DO 500 K=1,M
		DO 230 J=1,N
			Q0(J) = Q1(J)
			Q1(J) = Q(J)
230		CONTINUE
C.... CALCULATE R(K)=A*Q(K)
		CALL NLNVC(AA,BB,DD,Q,N,R,MAXA,TYP1,TYP2,TYP3)
C.... CALCULATE ALFA(K) USING ORTHOGONALITY
		ZNORM = 0.0D0
		DO 310 J=1,N
			ZNORM = ZNORM+Q1(J)*R(J)
310		CONTINUE 
		ALFA(K) = ZNORM
C.... CALCULATE Q(K+1)=R(K)-ALFA(K)*Q(K)
		DO 320 J=1,N
			Q(J) = R(J)-ALFA(K)*Q1(J)
320		CONTINUE 
		IF(K.EQ.1)GOTO 330
C.... CALCULATE BETA(K-1) USING ORTHOGONALITY
		ZNORM = 0.0D0
		DO 340 J=1,N
C SUNIL
C			ZNORM = ZNORM+Q0(J)*Q(J)
			ZNORM = ZNORM+Q0(J)*R(J)
340		CONTINUE
		BETA(K-1) = ZNORM
C.... CALCULATE Q(K+1)=R(K)-ALFA(K)*Q(K)-BETA(K-1)*Q(K-1)
		IF(K.EQ.M)GOTO 500
		DO 350 J=1,N
			Q(J) = Q(J)-BETA(K-1)*Q0(J)
350		CONTINUE
C.... CALCULATE U(K+1)=P(K)*P(K-1)*.....*P(1)*R(K+1)
C.... AND U(K+1) STORED INTO Q(I) (P(K): OPERATOR MATRIX)
C.... U(K+1)=(0,0,...,0,Q(K+1),Q(K+2),....,Q(N))
		DO 360 I=2,K
			II=((I-1)*(2*N-I))/2
			ZNORM = 0.0D0
			DO 370 J=I,N
				ZNORM = ZNORM+U(II+J)*Q(J)
370			CONTINUE	
			ZNORM = ZNORM/U(I)
			DO 380 J=I,N
				Q(J) = Q(J)-ZNORM*U(II+J)
380			CONTINUE
360		CONTINUE
C.... P(K+1) DETERMINED BY USING U(K+1) AND STORED IN DURRENT Q(K)
330		ZNORM = 0.0D0
		II = K+1
		DO 390 I=II,N
			ZNORM = ZNORM+Q(I)*Q(I)
390		CONTINUE
		IQ = (K*(2*N-II))/2
C SUNIL
C		IF(ZNORM.LT.1.0D-10)GOTO 400
C
		ZNORM = DSQRT(ZNORM)
		IF(Q(II).LT.0.0D0)ZNORM =-ZNORM
		U(IQ+II) = Q(II)+ZNORM
		U(II) = U(IQ+II)*ZNORM
		JJ = II+1
		IF(JJ.GT.N)GOTO 420
		DO 410 J=JJ,N
			U(IQ+J) = Q(J)
410		CONTINUE
		GOTO 420
C SUNIL
C400		U(II) = 1.0
C		DO 430 I=II,N
C			U(IQ+I)=0.0D0
C430		CONTINUE


C.... CALCULATE REOTHOGONALIZED K+1 LANCZOS VECTOR
C.... Q(K+1)=P(1)*P(2)*....*P(K)*P(K+1)*E(K+1) (P: OPERATOR MATRIX)
C.... Q(K+1) WILL BE STORED IN Q()
420		DO 440 J=1,N
			Q(J) = 0.0D0
440		CONTINUE
		Q(II) = 1.0D0
		DO 450 IP=1,K
			I = II+1-IP
			JJ = ((I-1)*(2*N-I))/2
			ZNORM = 0.0D0
			DO 460 J=I,N
				ZNORM = ZNORM+U(JJ+J)*Q(J)
460			CONTINUE
			ZNORM = ZNORM/U(I)
			DO 470 J=I,N
				Q(J) = Q(J)-ZNORM*U(JJ+J)
470			CONTINUE 
450		CONTINUE
500   CONTINUE
      RETURN
      END
C
C=====================================================================
      SUBROUTINE NLNVC(A,B,D,Q,NEQ,R,MAXA,TYP1,TYP2,TYP3)
      IMPLICIT REAL*8(A-H,O-Z)
	CHARACTER*4 TYP1,TYP2,TYP3
      DIMENSION MAXA(1)
      DIMENSION Q(1),R(1)
C
C *****************************************************************
C **  CALCULATE L(-1)*M*L(-T)*Q(I)                               **
C **  ----------------------------                               **
C **  THIS IS A SUBROUTINE TO CALCULATE A*Q(I)=L(-1)*M*L(-T)Q(I) **
C **  THE COMPUTATION OF L(-1) IS AVOIDED IN THIS SUBROUTINE BY  **
C **  USING THE FOLLOWING TECHNIQUE.                             **
C **                                                             **
C **  R(I)=A*Q(I)=L(-1)*M*L(-T)*Q(I), WHERE A=LL(T)              **
C **                                                             **
C **  THE S(I) IS EVALUATED BY FIRST SETTING                     **
C **                                                             **
C **  L(T)*S(I)=Q(I)                                             **
C **                                                             **
C **  SOLVE FOR S(I) TO GET L(-T)*Q(I). A*Q(I) WILL BE OBTAINED  **
C **  BY SOVING THE MATRIX EQUATION                              **
C **                                                             **
C **  L*(A*Q(I))=M*S(I)                                          **
C **                                                             **
C *****************************************************************
C
C     VARIABLES USED IN THIS SUBROUTINE
C
C     Q    :INPUT VECTOR Q(I)
C     NEQ  :ORDER OF THE ORIGINAL EIGENSYSTEM
C     R    :OUTPUT VECTOR A*Q(I)
C     JDIAG:DIAGONAL POINTERS OF THE TRIANGULAR MATRIX
C     A    :TRIANGULAR MATRIX L
C     B    :LUMPED MASS MATRIX IN A VECTOR FORM OR
C           CONSISTENT MASS MATRIX IN A PROFILE FORM
C     DFL  :LOGICAL FLAG TO DISTINGUISH THE TYPE OF MASS MATRIX
C
C     SUBROUTINE TO BE CALLED
C
C     MULPRO :MULTIPLICATION OF M*S(I) IN THE INSTRUCTION ABOVE
C
C *****************************************************************
C
	DIMENSION A(1),B(1),D(1)

C.... CALCULATE S(I) IN L(T)*S(I)=Q(I) BY BACK SUBSTITUTION
	CALL BKCHOL(MAXA,A,D,Q,NEQ,2,TYP3)

C.... MULTIPLICATION OF M*S(I) (M :CONSISTENT MASS)
	R(1:NEQ) = 0.0D0
	CALL MAMULT(MAXA,B,Q,R,TYP2,'STD')  

C.... CALCULATE A*Q(I) IN L*(A*Q(I))=M*S(I) BY FORWARD ELIMINATION
	CALL BKCHOL(MAXA,A,D,R,NEQ,1,TYP3)

      RETURN

      END
C
C=====================================================================
      SUBROUTINE TDBS(E,ALFA,BETA,S,RL,RU,NR,M,EPS,ACC,EVAL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C	Output EVAL added Oct11,2003 by NguyenDV to store Eigenvalues of original system
C
C *****************************************************************
C **  BISECTION METHOD THROUGH STURM SEQUENCE                    **
C **  ---------------------------------------                    **
C **  THIS IS A SUBROUTINE FOR SOLVING EIGENVALUES OF THE TRI    **
C **  DIAGONAL MATRIX DERIVED IN SUBROUTINE GOLUB BY BISECTION   **
C **  METHOD.                                                    **
C **  BRIFE ALGORITHM                                            **
C **  ---------------                                            **
C **  DETERMINATION OF UPPER AND LOWER BOUNDS OF EIGENVALUES     **
C **                            ]                                **
C **  STURM SEQUENCE CHECK                                       **
C **                            ]                                **
C **  DETERMINATION OF ACTUAL NUMBER OF EIGENVALUES              **
C **                            ]                                **
C **  CALCULATION OF EIGENVALUES USING BISECTION                 **
C **  METHOD WITH STURM SEQUENCE CHECK                           **
C **                                                             **
C **                                                             **
C **  AT FIRST THE LARGEST EIGENVALUES WANTED (UPPER BOUND LIMIT)**
C **  AND LOWER BOUND IS DETERMINED. SUITABLE VALUES FOR UPPER   **
C **  BOUND AND LOWER BOUND ARE +-|C|, WHERE |C|=MAX{|BETA(I)|+  **
C **  |ALFA(I)|+|BETA(I-1)|,|ALFA(1)|+|BETA(1)|}.                **
C **                                                             **
C *****************************************************************
C
C     VARIABLES USED IN THIS SUBROUTINE
C
C     E    :EIGENVALUES OF REDUCED SYSTEM
C     ALFA :ALFA IN TRIDIAGONAL MATRIX
C     BETA :BETA IN TRIDIAGONAL MATRIX
C     S    :BETA(I-1)*BETA(I-1)
C     RL   :LOWER BOUND OF EIGENVALUE
C     RU   :UPPER BOUND OF EIGENVALUE
C     NR   :REQUIRED NUMBER OF EIGENVALUES
C     M    :ORDER OF REDUCED SYSTEM
C     EPS  :ACCURACY OF EIGENVALUE
C     GL   :REDETERMINED LOWER BOUND OF EIGENVALUE
C     GU   :REDETERMINED UPPER BOUND OF EIGENVALUE
C     EPS  :CONVERGENCY TOLERANCE FOR BISECTION METHOD ADJUSTED TO
C           THE VALUES OF UPPER AND LOWER BOUNDS
C     EIGEN:EIGENVALUES OF THE ORIGINAL EIGENSYSTEM
C     FREQU:NATURAL CIRCULAR FREQUENCIES OF THE ORIGINAL SYSTEM
C
C     SUBROUTINE TO BE CALLED
C
C     STURM :STURM SEQUENCE CHECK FOR BISECTION METHOD
C
C *****************************************************************

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION E(1),ALFA(1),BETA(1),S(1),EVAL(1)
C	Array EVAL added Oct11,2003 by NguyenDV to store Eigenvalues of original system
C	Next line addded Oct11,2003 by NguyenDV
      EVAL(1) = 0.0D0
c
C     DETERMINE UPPER AND LOWER BOUNDS OF EIGENVALUES
      GL = RL
      GU = RU
	RNORM = DABS(ALFA(1))+DABS(BETA(1))
      DO 10 I=2,M
		X = DABS(BETA(I-1))+DABS(ALFA(I))+DABS(BETA(I))
		IF(X.GT.RNORM)RNORM=X
10	CONTINUE
C SUNIL
      RNORM = RNORM*1.1D0
      IF(GL.LT.-RNORM)GL = -RNORM
      IF(GU.GT.RNORM)GU = RNORM
      ACC = EPS*(GU-GL)
C.... DETERMINE THE ACTUAL NUMBER OF EIGENVALUES BETWEEN LOWER
C.... AND UPPER BOUND BY STURM SEQUENCE PROPERTY CHECK 
      S(1) = 0.0D0
C SUNIL
C	EPSS=ACC*ACC
	DO 12 I=2,M
		S(I) = BETA(I-1)*BETA(I-1)
C		IF(S(I).LT.EPSS)S(I) = EPSS
12	CONTINUE

C.... COUNT NO OF EIGEN VALUES GREATER THAN LOWER BOUND GL
      CALL STURM(ALFA,S,GL,ML,M)
C.... COUNT NO OF EIGEN VALUES GREATER THAN UPPER BOUND GU
      CALL STURM(ALFA,S,GU,MU,M)
      NNR = ML-MU
      IF(NR.GT.NNR)NR = NNR
C.... FIND THE EIGENVALUES BY BISECTION METHOD
      IF(NR.EQ.0)RETURN
	ML     = MU+NR
	RLAMDA = GL
      DO 30 K=1,NR
		I      = NR-K+1
		ML     = ML-1
C SUNIL
C		X      = RLAMDA-ACC
		X      = RLAMDA-DABS(RLAMDA)*1.0D-09
		Y      = GU
15		RLAMDA = 0.5D0*(X+Y)
C SUNIL
C		IF(Y-X.LE.ACC)GOTO 29
		IF((DABS(DABS(Y)-DABS(X)))/DABS(RLAMDA).LE.1.0D-09)GOTO 29
		CALL STURM(ALFA,S,RLAMDA,IS,M)
		IF(IS-ML)25,25,20
20		X = RLAMDA
		GOTO 15
25		Y = RLAMDA
		GOTO 15
29		E(I) = RLAMDA
30	CONTINUE


C.... PRINT EIGENVALUES OF ORIGINAL SYSTEM

      DO 40 I=1,NR
	EIGEN = 1.0D0/E(I)

	CALL RELFILL('EIVV',EIGEN,1,I,1) !STORE EIGEN VALUE

C     Jan 17, 2003 added next line to make {E} or W(N8) contain egenvalues
C	Changed	to the second line, Oct11,2003, NguyenDV
C	E(I) = EIGEN
	EVAL(I) = EIGEN
40	CONTINUE

      RETURN

 900  FORMAT (//32X,16(1H*)/32X,1H*,14X,1H*/
     1        32X,16H* JOB PROGRESS */32X,1H*,14X,1H*/32X,16(1H*)//)
 1100 FORMAT (1X,A2,I3,4I7,4X,A6,E11.3,2X,24H * DIVERGENCE OCCURRED *)
1000  FORMAT(//8X,'E I G E N S O L U T I O N   B Y   L A N C Z O S',
     *  3X,'M E T H O D'/8X,61(1H-)//
c     +  04X,70(1H-)//
     *  19X,'O U T P U T   OF   E I G E N V A L U E S'/19X,40(1H-)//
     *  14X,'REQUIRED NUMBER OF EIGEN VALUES. . . . . . =',I5/
     +  14X,'SIZE OF SUBSPACE . . . . . . . . . . . . . =',I5//
     *  14X,'MODE NUMBER',27X,'EIGENVALUE',/14X,11(1H-),27X,10(1H-))
1001  FORMAT(//8X,'E I G E N S O L U T I O N   B Y   L A N C Z O S',
     *  3X,'M E T H O D'/8X,61(1H-)//
c     +  04X,70(1H-)//
     *  19X,'O U T P U T   OF   E I G E N V A L U E S'/19X,40(1H-)//
     *  14X,'REQUIRED NUMBER OF EIGEN VALUES. . . . . . =',I5/
     *  14X,'MODE NUMBER',27X,'EIGENVALUE',/14X,11(1H-),27X,10(1H-)/)
2000  FORMAT(15X,I5,24X,E19.6)
2100  FORMAT(I5,7X,'WARNING: NEGATIVE EIGEN VALUE IS DETECTED')

      END
C
C=====================================================================
      SUBROUTINE STURM(ALFA,S,RLAMDA,IS,M)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ALFA(1),S(1)
C
C ******************************************************************
C **  STURM SEQUENCE PROPERTY                                     **
C **  -----------------------                                     **
C **  THIS SUBROUTINE IS TO CHECK THE STURM SEQUENCE PROPERTY     **
C **  OF A TRIDIAGONAL MATRIX. STURM SEQUENCE IS DETERMINED AS    **
C **  FOLLOWS.                                                    **
C **  P_0(X) = 1                                                  **
C **  P_1(X) = ALFA(1)-X                                          **
C **  P_I(X) = (ALFA(I)-X)*P_(I-1)(X)-BETA(I-1)**2*P_(I-2)(X)     **
C **           (I = 2,3,....,M)                                   **
C **                                                              **
C **  S IS THE NUMBER OF EIGENVALUES OF THE TRI DIAGONAL MATRIX   **
C **  WHICH ARE GREATER THAN X, IF S IS EQUAL TO THE NUMBER OF    **
C **  AGREEMENTS IN SIGNS OF SUCCESSIVE MEMBERS OF SEQUENCE P_I(X)**
C ******************************************************************
C
C     VARIABLES USED IN THIS SUBROUTINE
C
C     S(I)  :BETA(I-1)**2
C     RLAMDA:THE VALUE X IN THE DESCRIPTION ABOVE
C     IS    :THE NUMBER OF EIGENVALUES GREATER THAN THE VALUE OF RLAMDA
C     M     :ORDER OF REDUCED EIGENSYSTEM
C
C ******************************************************************
C
C SUNIL
	IS = 0
      X  = 0.0D0
      Y  = 1.0D0
      DO 35 I=1,M
		Z  = (ALFA(I)-RLAMDA)*Y-S(I)*X
		AZ = DABS(Z)
		IF(AZ.LT.1.0D30)GOTO 5
		Z  = Z *1.0D-30
		Y  = Y *1.0D-30
		GOTO 10
5		IF(AZ.GT.1.0D-30)GOTO 10
		Z  = Z *1.0D30
		Y  = Y *1.0D30
10		IF(Z)15,15,20
15		IF(Y)25,25,30
20		IF(Y)30,25,25
25		IS = IS+1
30		X  = Y
		Y  = Z
35	CONTINUE
      RETURN
      END
C
C=====================================================================
      SUBROUTINE TDII(Z,E,ALFA,BETA,S,Q,Q0,Q1,R,NR,M,NN,EPS,ACC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Z(1),E(NR),ALFA(M),BETA(M),S(M),Q(M),Q0(M),Q1(M),R(M+2)
C
C *****************************************************************
C     INVERSE ITERATION METHOD
C     ------------------------
C     INVERSE ITERATION ALGORITHM IS USED TO EVALUATE THE 
C     CORRESPONDING EIGENVECTORS OF TRIDIAGONAL MATRIX, AN 
C     ORTHOGONALIZATION TECHNIQUE IS USED (IN SUBROUTINE ORTH)
C     TO ISOLATE EIGENVECTORS CORRESPONDING TO CLOSE  SETS OF 
C     EIGENVALUES.
C *****************************************************************
C
C     VARIABLES USED IN THIS SUBROUTINE
C
C     Z    :EIGENVECTORS OF REDUCED SYSTEM
C     E    :EIGENVALUES OF REDUCED SYSTEM
C     ALFA :ALFA IN TRIDIAGONAL MATRIX
C     BETA :BETA IN TRIDIAGLNAL MATRIX
C     NR   :NUMBER OF REQUIRED EIGENVALUES
C     M    :ORDER OF REDUCED EIGENSYSTEM
C     EPS  :ACCURACY OF EIGENVECTORS
C     ACC  :EPS*(GU-GL)
C     NN   :DESCRIBED IN LNCZS
C
C     SUBROUTINE TO BE CALLED
C
C     ORTH :ORTHOGONALIZATION PROCESS TO ISOLATE EIGENVECTORS
C           CORRESPONDING TO CLOSE SETS OF EIGENVALUES
C
C ******************************************************************
C
      IF(M.EQ.1)RETURN
      MM1    = M-1
      MP1    = M+1
      MP2    = M+2
      ID     = 0
      TLAM   = 1.0D70
      R(MP1) = 0.0D0
      R(MP2) = 0.0D0
      Q0(M)  = 0.0D0
	ACCD   = ACC
      EPSD = 1000.0D0*ACCD
C SUNIL
C      DO 7 J = 1,MM1
C		IF(DABS(BETA(J)).LT.ACCD)BETA(J) = ACCD
C7	CONTINUE
      DO 100 I = 1,NR
		IN = ((I-1)*NN+1)/2
		RLAMDA = E(I)
C		IF(RLAMDA.GT.TLAM-ACCD)RLAMDA = TLAM-ACCD
		S(1) = ALFA(1)-RLAMDA
		IF(S(1).EQ.0.0D0) S(1) = ACCD
		Q(1) = BETA(1)
		R(M) = 1.0D0
C ...	LU DECOMPOSITION WITH PIVOTING
		DO 20 J=1,MM1
			R(J) =1.0D0
			JP1 = J+1
			IF(DABS(BETA(J))-DABS(S(J)))10,10,15
10			Q1(J) = BETA(J)/S(J)
			S(JP1) = ALFA(JP1)-RLAMDA-Q1(J)*Q(J)
			Q(JP1) = BETA(JP1)
			Q0(J) = 0.0D0
C SUNIL
			IF(S(JP1).EQ.0.0D0)S(JP1) = ACCD
			GOTO 20
15			Q1(J) = S(J)/BETA(J)
			S(JP1) = Q(J)-Q1(J)*(ALFA(JP1)-RLAMDA)
			Q(JP1) = -Q1(J)*BETA(JP1)
			S(J) = BETA(J)
			Q(J) = ALFA(JP1)-RLAMDA
			Q0(J) = BETA(JP1)
			IF(S(JP1).EQ.0.0D0)S(JP1) = ACCD
20		CONTINUE
C....	INVERSE ITERATION IS PERFORM ONLY FOR 3 STEPS. IN THE FIRST STEP 
C.... FORWARD ELIMINATION IS OMITTED
		IT = 0
30		IT = IT+1
C....	BACK SUBSTITUTION
		DO 40 J = M,1,-1
			R(J) = (R(J)-Q(J)*R(J+1)-Q0(J)*R(J+2))/S(J)
			IF(DABS(R(J)).LT.1.0D30)GOTO 40
			DO 35 K = 1,M
				R(K) = R(K)*1.0D-30
35			CONTINUE
40		CONTINUE
C SUNIL
		IF(IT.GE.3)GOTO 60
C.... FORWARD ELIMINATION
		DO 55 J = 1,MM1
			IF(Q0(J))50,45,50
45			R(J+1) = R(J+1)-Q1(J)*R(J)
			GOTO 55
50			TEMP = R(J)
			R(J) = R(J+1)
			R(J+1) = TEMP-Q1(J)*R(J)
55		CONTINUE
		GOTO 30
C SUNIL
60		IF(TLAM-RLAMDA.GT.EPSD)GOTO 87
C60		IF(TLAM-RLAMDA.GT.DABS(RLABDA)*1.0D-06)GOTO 87
		ID = ID+1
		CALL ORTH(Z,R,Q1,I,ID,M,NN)
		GOTO 92
87		ID = 0
92		TLAM = RLAMDA
		SUM = 0.0D0
		DO 90 J = 1,M
			SUM = SUM+R(J)*R(J)
90		CONTINUE
		TEMP = 1.0D0/DSQRT(SUM)
		DO 95 J = 1,M
			Z(J+IN) = TEMP*R(J)
95		CONTINUE
100   CONTINUE
      RETURN
      END
C
C=====================================================================
      SUBROUTINE ORTH(Z,R,Q1,K,ID,M,NN)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Z(1),R(1),Q1(1)
C
      JN = K-ID-2
      DO 10 I = 1,ID
		IN = ((JN+I)*NN+1)/2
		SUM = 0.D0
		DO 5 J = 1,M
			SUM = SUM+Z(J+IN)*R(J)
5		CONTINUE
		Q1(I) = SUM
10	CONTINUE
      DO 20 I = 1,M
		SUM = 0.D0
		DO 15 J = 1,ID
			IJ =  ((JN+J)*NN+1)/2
			SUM = SUM+Q1(J)*Z(IJ)
15		CONTINUE
		R(I) = R(I)-SUM
20	CONTINUE
      RETURN
      END
C
C=====================================================================
      SUBROUTINE BKTFM(ID,N,M,NR,TVEC,EVEC,R,Y,U,NN,MAXA,A,D,TYP)
      IMPLICIT REAL*8(A-H,O-Z)
	CHARACTER*4 TYP
C
C *****************************************************************
C **  BACK TRANSFORMATION OF EIGENVECTORS                        **
C **  -----------------------------------                        **
C **  THIS SUBROUTINE CALCULATES THE EIGENVECTORS OF A           **
C **  (=L(-1)*M*L(-T)) USING THE EIGENVECTORS OBTAINED IN THE    **
C **  SUB TDII AND THE HOUSEHOLDER REFLECTION MATRICES OBTAINED  **
C **  IN THE SUB GOLUB. EIGENVECTOR Y IS OBTAINED AS FOLLOWS:    **
C **                                                             **
C **  Y = P(1)*P(2)*...*P(M)*{Z}                                 **
C **                                                             **
C **  WHERE Z IS EIGENVECTOR OF THE REDUCED SYSTEM               **
C **        P(I) IS I-TH HOUSEHOLDER REFLECTION MATRIX           **
C *****************************************************************
C
C     VARIABLES USED IN THIS SUBROUTINE
C
C     N    :THE SIZE OF ORIGINAL SYSTEM
C     M    :THE SIZE OF REDUCED SYSTEM
C     NR   :THE NUMBER OF REQUIRED EIGENFREQUENCIES
C     TVEC :EIGENVECTORS OF REDUCED TRIDIAGONAL MATRIX 
C           OBTAINED IN TDII
C     EVEC :EIGENVECTORS OF ORIGINAL SYSTEM
C     R    :THE VECTOR AFTER THE OPERATION OF P(I)*V WHICH BECOMES
C           FINALLY Y (=EIGENVECTOR OF J-TH MODE).
C     U    :VECTORS TO FORM OPERATOR MATRIX
C
C     SUBROUTINE TO BE CALLED
C
C     BKSUB  : BACK SUBSTITUTION OF EIGENVECTORS OF STANDARD
C              EIGENSYSTEM INTO THOSE OF GENERALIZED EIGENSYSTEM
C
C ******************************************************************
C

      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)

      DIMENSION TVEC(1),EVEC(1),R(N),U(1),Y(1)
      DIMENSION ID(1),MAXA(1),A(1),D(1)

      MM1 = M-1
      MP1 = M+1
      JN = -N
      J1 = (MM1*(2*N-M))/2
	WRITE (ISO,1000)

      DO 22 J = 1,NR
		IN = ((J-1)*NN+1)/2
		I1 = J1
C.... FIRST OPERATION OF P(M)*Z (WHERE Z: TVEC)
		SUM = U(I1+M)*TVEC(IN+M)/U(M)
		DO 33 I=1,MM1
			R(I) = TVEC(IN+I)
33		CONTINUE
		R(M) = TVEC(IN+M)-SUM*U(I1+M)
		IF(MP1.GT.N)GOTO 45
		DO 44 JJ = MP1,N
			R(JJ) = -SUM*U(I1+JJ)
44		CONTINUE
C.... OPERATION OF P(I)*R(I=M-1,1)
45		DO 55 II = 2,MM1
			I = MP1-II
			I1 = I1-N+I
			IF(DABS(U(I)-1.0D0).LT.1.0D-16)GOTO 55
			SUM = 0.D0
			DO 66 JJ = I,N
				SUM = SUM+U(I1+JJ)*R(JJ)
66			CONTINUE
			SUM = SUM/U(I)
			DO 77 JJ = I,N
				R(JJ) = R(JJ)-SUM*U(I1+JJ)
77			CONTINUE
55		CONTINUE
C.... TRANSFORMATION OF EIGENVECTORS OF A TO THOSE OF THE ORIGINAL
C.... SYSTEM WITH THE INVERSE OF X = L(-1)*Y. 
C....	SUBROUTINE BKSUB PERFORM THE BACK SUBSTITUTION
		CALL BKCHOL(MAXA,A,D,R,N,2,TYP)
		JN = JN+N

C	SUNIL
		CALL NORM(Y,R,N)
		CALL PRIDIS(ID,Y,J)

     		DO 111 I = 1,N
			EVEC(JN+I) = Y(I)
111		CONTINUE
22    CONTINUE

      RETURN

 1000 FORMAT (//17X,'O U T P U T   O F   E I G E N V E C T O R S'/
     +		  17X,'-------------------------------------------')
      END
C
C=====================================================================
      SUBROUTINE NORM(X,Y,N)
      IMPLICIT REAL*8(A-H,O-Z)
C.... NORMALIZE VECTOR Y TO UNIT VECTOR X
      DIMENSION X(1),Y(1)
      SCALE = DSQRT(DOT(Y,Y,N))
      IF(SCALE.GT.0.D0) SCALE = 1.D0/SCALE
      DO 100 I = 1,N
100   X(I) = Y(I)*SCALE
      RETURN
      END
C
C=====================================================================
      FUNCTION DOT(A,B,N)
      IMPLICIT REAL*8(A-H,O-Z)
C.... VECTOR DOT PRODUCT
      DIMENSION A(1),B(1)
      DOT = 0.D0
      DO 100 I = 1,N
 100  DOT = DOT+A(I)*B(I)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE PRIDIS(ID,R,J)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	SUNIL 16/02/01 PROGRAMMED BY SUNIL
C	----------------------------------------
C	PRINT DISPLACEMENTS AT EACH NODAL POINTS
C	----------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),ISOLOP
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON A(9000000),IA(9000000)
C	-------------------------------------------------------------------
      DIMENSION ID(NSF,1),R(1)


	CALL MODOUT(ID,R,J)          !STORE MODE SHAPE TO FILE


	RETURN

 3000 FORMAT (//19X,16HEIGENVALUE LAMDA,I3,3H  =,E20.12/19X,42(1H-))
      WRITE (ISO,3100)J
 3100 FORMAT (//26X,'EIGENVECTOR OF MODE ='I5/26X,26(1H-)//
     1          1X,5HJOINT,7X,6HX-DISP,6X,6HY-DISP,6X,6HZ-DISP,
     2          6X,6HX-ROTA,6X,6HY-ROTA,6X,6HZ-ROTA/1X,79(1H-))
 3200 FORMAT (I5,3X,6E12.4)

C	NEXT FIVE LINES ADDED BY GILSON - NOV2003 - FOR GIDOUT
7000	FORMAT (I5,3X,7E12.4)
8000	FORMAT ('Displacement     4     'I4'.0     Estatic     2   1   1')
8100	FORMAT ('X-Displacement')
8200	FORMAT ('Y-Displacement')
8300	FORMAT ('Z-Displacement')

8010	FORMAT ('Result "ModeShapes"  "ModeShapes"  'I4' Vector OnNodes')
8110	FORMAT ('ComponentNames  "X-Component", "Y-Component",',
     #' "Z-Component"')
8210	FORMAT ('Values')
8310	FORMAT ('End Values',/)


8011	FORMAT ('Result "ModeShapesMat" "ModeShapes" 'I4' Matrix OnNodes')
8111	FORMAT ('ComponentNames  "X", "Y", "Z", "RX", "RY", "RZ" ')


	END
C
C=====================================================================
C	END OF SUBROUTINE FOR LANCZOS EIGEN SOLVER BY SUNIL
C=====================================================================
C
C=====================================================================
C	START OF SUBROUTINE FOR INVERSE ITERATION EIGEN SOLVER BY SUNIL
C=====================================================================
C
C=====================================================================
      SUBROUTINE INVERSIT(ID,JDIAG,A,B,W,NEQ,NR,MM)
      IMPLICIT REAL*8(A-H,O-Z)
C	------------------------------------------------------------
C	ALL THE SUBROUTINES FOR EIGEN SOLUTOION BY INVERSE ITERATION
C	WERE IMPLEMENTED BY SUNIL
C	FIND EIGEN VALUES AND VECTORS BY INVERSE ITERATION
C	------------------------------------------------------------

	COMMON /MEMO/ MEMA,MEMI,LASTA,LASTI,NELEA,NELEI
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
	COMMON /WLAS/ MERW,MEIW,NLASI,NLASW
	
      DIMENSION ID(1),JDIAG(1),A(1),B(1),W(1)

C.... ALLOCATION OF THE WORK SPACE

	NHALF = (NEQ*(NEQ+1))/2
      NAU   = NEQ*NEQ
      NTT   = (2*NHALF+NAU)

C	SUNIL 20/02/01 
      IF(NR.LE.0  )NR = 5
      IF(NR.GT.NEQ)NR = NEQ

      IF(MM.LT.0  )MM = 0
      IF(MM.GT.NR )MM = NR

      NR11  = 1
      NR12  = NR11+NHALF
      NR13  = NR12+NHALF
      NR14  = NR13+NAU
      NSIZE = NR14+NEQ

C.... MEMORY REQUIREMENT CHECK
      IF(NSIZE.LT.MERW)GOTO 10
      WRITE(ISO,1000)MERW,NSIZE
      RETURN

10	CALL CLEARA(W,NSIZE)
      CALL CONVT(A,W(NR11),NEQ,JDIAG,0)
      JJ = 0
      CALL CONVT(B,W(NR12),NEQ,JDIAG,JJ)
	CALL DMTS1(W(NR11),W(NR12),NEQ,W(NR13),IER)
      IF(IER.NE.0)WRITE(ISO,4006)IER
      CALL CLEARA(W(NR13),NAU)
      W(NR13) = NR
      CALL DMEVJ(W(NR11),NEQ,W(NR13),NEQ,IER)
      IF(IER.NE.0)WRITE(ISO,4007)IER

      DO 4008 JJ=1,NEQ
		NEZ = NEQ*(JJ-1)+NR13
		CALL BSUB(W(NR12),W(NEZ),JDIAG,NEQ)
4008  CONTINUE

      CALL PSORTE(W(NR11),W(NR13),NEQ,EIGVEC)
      WRITE(100,1320)NR,MM
      WRITE(ISO,1320)NR,MM
      NUMMOD=NR
C	IF(NUMMOD.LE.0)NUMMOD=5
C	IF(NUMMOD.GE.NEQ)NUMMOD=NEQ
      DO 1305 I=1,NUMMOD
		EIGEN = W(NR11-1+I)
		WRITE(100,1330)I,EIGEN
			WRITE(ISO,1330)I,EIGEN
1305	CONTINUE
      IF(MM.EQ.0)GOTO 130
	NMM1=NR13
	DO 1310 K = 1,MM
c		WRITE(100,1300)K
			WRITE(ISO,1300)K
C	SUNIL
		CALL NORM(W(NR14),W(NMM1),NEQ)
C		CALL PRTDIS(ID,X,W(NMM1),F,NDM,NEF,1.0D0,1.0D0)
C		CALL PRTDIS(ID,X,W(NR14),F,NDM,NEF,1.0D0,1.0D0)
C	SUNIL
		CALL PRIDIS(ID,(NR14),K)
		NMM1 = NMM1+NEQ
1310	CONTINUE

130	RETURN
1000  FORMAT(/,5X,43H*******  WARNING : MEMORY SHORTAGE  *******,/,
     *       /,14X,18HMEMORY CAPACITY = ,I6,/,14X,
     *       18HMEMORY REQUIRED = ,I6)
1300  FORMAT (//15X,'O U T P U T   O F   E I G E N V E C T O R S'
     * //'MODE NUMBER  -------> ',I5/)
1320  FORMAT(//6X,'E I G E N S O L U T I O N   B Y   I N V E R S E'
     *       ,3X,'I T E R A T I O N'//
     +     12X,51(1H-)//
     *       //15X,'O U T P U T   O F   E I G E N V A L U E S'//
     *         16X,'REQUIRED NUMBER OF EIGENVALUES  ='I5/
     *         16X,'REQUIRED NUMVER OF EIGENVECTORS ='I5//
     *         20X,'MODE NUMBER',10X,'EIGENVALUE'/)
1330  FORMAT(20X,I5,7X,E19.6)

4006  FORMAT(5X,27H**ERROR MESSAGE FROM DMTS1 ,I5)
4007  FORMAT(5X,27H**ERROR MESSAGE FROM DMEVJ ,I5)
	END
C
C=====================================================================
      SUBROUTINE CONVT(A,B,NEQ,JDIAG,KODE)
      IMPLICIT REAL*8 (A-H,O-Z)
C.... SUBROUTINE TO CONVERT PROFILE MATRIX (KODE=0) OR 
C.... DIAGONAL MATRIX (KODE.NE.0) TO TRIANGULAR FORM
      DIMENSION A(1),B(1),JDIAG(1)
      NTOT = NEQ*(NEQ+1)/2
      DO 10 I=1,NTOT
10    B(I) = 0.D0
      LA = 1
      LB = 1
      DO 200 I=1,NEQ
      I2 = I*(I+1)/2
      IF(KODE.GT.0) B(I2) = A(I)
      IF(KODE.GT.0) GOTO 200
      J1 = 0
      IF(I.GT.1) J1 = JDIAG(I-1)
      J2 = JDIAG(I)
      I1 = I*(I-1)/2
      IB = I2-I1
      IA = J2-J1
      IVOID = IB-IA
      DO 100 K=1,IB
      IF(K.LE.IVOID) GOTO 100
      B(LB) = A(LA)
      LA = LA+1
100   LB = LB+1
200   CONTINUE
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DMTS1(A,B,N,AUX,IER)
C	-----------------------------------------------------------
C	REDUCTION OF THE REAL SYMMETRIC EIGENPROBLEM A*X=B*X*LAMBDA
C	TO STANDARD FORM
C	-----------------------------------------------------------
      DIMENSION A(1),B(1),AUX(1)
      DOUBLE PRECISION A,B,AUX,S

      ISW  =IER
      IER  =1000
      IF (N) 18,18,1
    1 IER  =-12345
      CALL DMFS(B,N,IER)
      IF (IER-10) 2,2,18
    2 LS   =1
      NN   =N*N
      IIDP =(NN+N)/2
      K1   =0
      DO 13 K=1,N
         LK   =LS
         LI   =1
         L1   =0
         IF (K1) 8,8,3
    3    DO 7 L=1,K1
            S    =A(LK)
            IDP  =IIDP
            IF (L1) 6,6,4
    4       DO 5 I=1,L1
               IDP  =IDP+1
               S    =S-AUX(IDP)*B(LI)
    5          LI   =LI+1
    6       L1   =L
            LK   =LK+1
            IDP  =IDP+1
            AUX(IDP)=S/B(LI)
    7       LI   =LI+1
    8    LS   =LK+1
         DO 12 L=K,N
            S    =A(LK)
            IDP  =IIDP
            IF (L1) 11,11,9
    9       DO 10 I=1,L1
               IDP  =IDP+1
               S    =S-AUX(IDP)*B(LI)
   10          LI   =LI+1
   11       L1   =L
            IDP  =IDP+1
            AUX(IDP)=S/B(LI)
            AUX(LK)=AUX(IDP)
            LK   =LK+L
   12       LI   =LI+1
   13    K1   =K
      LK   =1
      DO 17 K=1,N
         L1   =0
         LI   =1
         DO 17 L=1,K
            S    =AUX(LK)
            IF (L1) 16,16,14
   14       DO 15 I=1,L1
               IDP  =IIDP+I
               S    =S-AUX(IDP)*B(LI)
   15          LI   =LI+1
   16       IDP  =IIDP+L
            AUX(IDP)=S/B(LI)
            A(LK)=AUX(IDP)
            LK   =LK+1
            LI   =LI+1
   17       L1   =L
      IF (IER) 20,20,18
   18 IF (ISW+12345) 19,20,19
   19 CALL WIER(IER,10341)
   20 RETURN
      END
C
C=====================================================================
      SUBROUTINE DMFS(A,N,IER)

      DOUBLE PRECISION SUM,F1,F2,A,EPS
      DIMENSION A(1)

      INDER=IER
      IER  =0
      IF (N) 1,1,2
    1 IER  =1000
      GOTO 14
    2 EPS  =1D-15
      IB   =1
      IN   =0
      DO 13 K=1,N
         KL   =0
    3    SUM  =0.D0
         L    =IB
    4    IF (L-IN) 5,5,6
    5    KL   =KL+1
         F1   =A(L)
         F2   =A(KL)
         SUM  =SUM+F1*F2
         L    =L+1
         GOTO 4
    6    KL   =KL+1
         IN   =IN+1
         SUM  =A(IN)-SUM
         IF (IN-KL) 8,8,7
    7    A(IN)=SUM/A(KL)
         GOTO 3
    8    IF (SUM) 9,9,10
    9    IER  =IER+2000
         GOTO 14
   10    IF (SUM-DABS(EPS*A(IN))) 11,11,12
   11    IER  =10
   12    A(IN)=DSQRT(SUM)
   13    IB   =IB+K
      IF (IER) 14,16,14
   14 IF (INDER+12345) 15,16,15
   15 CALL WIER(IER,10221)
   16 RETURN
      END
C
C=====================================================================
      SUBROUTINE WIER(IER,NO)
      WRITE (100,1) NO,IER
    1 FORMAT (//' ***** DKO',I5,' RAISED ERROR INDICATOR TO ',I4,
     1     ' *****'///)
      RETURN
      END
C
C=====================================================================
      SUBROUTINE DMEVJ(A,N,EV,IEV,IER)
C	-------------------------------------------------------
C	EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC MATRIX
C	JACOBI'S METHOD WITH THRESHOLD
C	-------------------------------------------------------
      DIMENSION A(1),EV(1)

      DOUBLE PRECISION A,EV,ANORM,ANRMX,THR,U,X,Y,SINX,SINX2,
	1                 COSX,COSX2,E1,SINCS,FN,Z,TEMP,AUX
      EQUIVALENCE (SEV,ISEV)

      ISW  =IER
      IER  =0
      IF (N) 1,1,2
    1 IER  =1000
      GOTO 48
    2 SEV  =EV(1)
      IF (ISEV) 3,12,3
    3 IF (IABS(IEV)-N) 1,4,4
    4 INEV =1
      IF (IEV) 5,6,6
    5 IREV =-IEV
      ICEV =INEV
      IER  =1
      GOTO 7
    6 IREV =INEV
      ICEV =IEV
    7 IV1  =1
      DO 11 I=1,N
         IV   =IV1
         DO 10 J=1,N
            IF (I-J) 8,9,8
    8       EV(IV)=0.D0
            GOTO 10
    9       EV(IV)=1.D0
   10       IV   =IV+IREV
   11    IV1  =IV1+ICEV
   12 E1   =1.D-16
      FN   =N
      ANORM=0.D0
      J    =0
      I1   =N-1
      IF (I1) 50,50,13
   13 DO 14 I=1,I1
         J    =J+I
         J1   =J+1
         J2   =J+I
         DO 14 K=J1,J2
            AUX  =A(K)
   14       ANORM=ANORM+AUX*AUX
      IF (ANORM) 46,46,15
   15 ANORM=DSQRT(ANORM+ANORM)
      ANRMX=ANORM*E1/FN
      IND  =0
      THR  =ANORM
C.... INITIAL THRESHOLD
   16 THR  =THR/FN
   17 LQ   =0
      IF (ISEV) 18,19,18
   18 IV1  =1
   19 DO 43 L=1,I1
         MQ   =LQ+L
         L1   =L+1
         IF (ISEV) 20,21,20
   20    IV2  =IV1+ICEV
   21    DO 41 M=L1,N
            LM   =L+MQ
            Y    =A(LM)
            IF (DABS(Y)-THR) 39,22,22
   22       IND  =1
            LL   =L+LQ
            MM   =M+MQ
            Z    =A(LL)-A(MM)
            X    =Z*.5D0
            IF (DABS(Y)-DABS(X)) 26,26,23
   23       X    =-X/Y
            U    =1.D0/DSQRT(X*X+1.D0)
            IF (X) 24,25,25
   24       U    =-U
   25       X    =1.D0+X*U
            GOTO 27
   26       X    =-Y/X
            U    =1.D0/DSQRT(X*X+1.D0)
            TEMP =X*U
            X    =U+1.D0
            U    =TEMP
   27       X    =X+X
            SINX =U/DSQRT(X)
            SINX2=SINX*SINX
            COSX2=1.D0-SINX2
            COSX =DSQRT(COSX2)
            SINCS=SINX*COSX
            IQ   =0
            IF (ISEV) 28,29,28
   28       IVL  =IV1
            IVM  =IV2
   29       DO 38 I=1,N
               IF (I-L) 30,36,30
   30          IF (I-M) 31,36,32
   31          IM   =I+MQ
               GOTO 33
   32          IM   =M+IQ
   33          IL   =L+IQ
               IF (I-L) 34,36,35
   34          IL   =I+LQ
   35          X    =A(IL)*COSX-A(IM)*SINX
               A(IM)=A(IL)*SINX+A(IM)*COSX
               A(IL)=X
   36          IF (ISEV) 37,38,37
   37          X    =EV(IVL)*COSX-EV(IVM)*SINX
               EV(IVM)=EV(IVL)*SINX+EV(IVM)*COSX
               EV(IVL)=X
               IVL  =IVL+IREV
               IVM  =IVM+IREV
   38          IQ   =IQ+I
            X    =(Y+Y)*SINCS
            U    =A(LL)*COSX2+A(MM)*SINX2-X
            X    =A(LL)*SINX2+A(MM)*COSX2+X
            A(LM)=Z*SINCS+Y*(COSX2-SINX2)
            A(LL)=U
            A(MM)=X
   39       IF (ISEV) 40,41,40
   40       IV2  =IV2+ICEV
   41       MQ   =MQ+M
         IF (ISEV) 42,43,42
   42    IV1  =IV1+ICEV
   43    LQ   =LQ+L
      IF (IND) 45,45,44
   44 IND  =0
      GOTO 17
   45 IF (ANRMX-THR) 16,46,46
C.... RESTORE EIGENVALUES
   46 J    =1
      DO 47 I=2,N
         J    =J+I
   47    A(I) =A(J)
      IF (IER) 48,50,48
   48 IF (ISW+12345) 49,50,49
   49 CALL WIER(IER,10337)
   50 RETURN
      END
C
C=====================================================================
      SUBROUTINE BSUB(A,B,JDIAG,NEQ)
C	-------------------------------------------
C	BACK SUBSTITUTION OF FULL TRIANGULAR MATRIX
C	-------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1),B(1),JDIAG(1)
      J = NEQ
      NHA=NEQ*(NEQ+1)/2
      JD = NHA
      JR=NHA
800   B(J)=B(J)/A(JR)
      D = B(J)
      J = J - 1
      IF(J.LE.0) RETURN
      JR = (J*(J+1))/2
      IF(JD-JR.LE.1) GOTO 1000
      IS = J - JD + JR + 2
      K = JR - IS + 1
      DO 900 I = IS,J
 900  B(I) = B(I) - A(I+K)*D
1000  JD = JR
      GOTO 800
      END
C
C=====================================================================
      SUBROUTINE PSORTE(D,VEC,N,EIVEC)
      IMPLICIT REAL*8(A-H,O-Z)
C	------------------------------------
C	SORT EIGENVALUES TO ASSCENDING ORDER
C	------------------------------------
C	SUNIL
      LOGICAL EIVEC
      DIMENSION D(1),VEC(N,1)

      N1 = N - 1
      DO 200 I = 1,N1
		E = D(I)
		II = I
		I1 = I + 1
		DO 100 J = I1,N
C SUNIL
C			IF(DABS(E).GT.DABS(D(J)))GOTO 100
			IF(DABS(E).LE.DABS(D(J)))GOTO 100
			E = D(J)
			II = J
100		CONTINUE
		D(II) = D(I)
		D(I) = E
		IF(.NOT.EIVEC)GOTO 200
C.... SORT EIGENVECTORS IF REQUIRED
		DO 150 J = 1,N
			E = VEC(J,I)
			VEC(J,I) = VEC(J,II)
			VEC(J,II) = E
150		CONTINUE
200   CONTINUE
      RETURN
      END
C
C=====================================================================
C	END OF SUBROUTINE FOR INVERSE ITERATION EIGEN SOLVER BY SUNIL
C=====================================================================
C
C=====================================================================
C	START OF SUBROUTINE FOR SUBSPACE EIGEN SOLVER
C=====================================================================
C
C=====================================================================
      SUBROUTINE STABIL (AA,BB,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C	09/02/01 SUNIL IMPLEMENTED IN-CORE SUB-SPACE ITERATION METHOD
C	MODIFIED AND CHANGED ALL THE RELEVENT SUBROUTINES
C     ----------------------------------------------------------------
C     SOLUTION OF STANDARD AND GENERALIZED EIGENVALUE PROBLEMS BY
C	SUB-SPACE INTERATIONS INCLUDING SHIFTING OF NON +VE DEFINIT [K]
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C
      COMMON A(9000000),IA(9000000)
C
      DIMENSION AA(1),BB(1)
     
C	--------------
C     INITIALISATION
C     --------------
      CALL CPU_TIME (TIM1)
C
      IF (NROOT.LE.0) NROOT = 1
      IF (NROOT.GT.NEQ) NROOT = NEQ
	IF (IVPRT.LT.0) IVPRT = 0
	IF (IVPRT.GT.NROOT) IVPRT = NROOT
      NEIG = NEIG+1

      IF (IEIG.EQ.1) WRITE (ISO,3100)
      IF (IEIG.EQ.2) WRITE (ISO,3200)
      IF (IEIG.EQ.3) WRITE (ISO,3300)
      WRITE (ISO,3500) NSEIG,NROOT,IVPRT,NC,NITEM,IFSS,EPS,SHIFT0
	GOTO 150

 100	WRITE (ISO,9999)
	STOP
C     -------------------------------------------------------------
C     READ TANGENTIAL STIF MATRIX (2) AND GEO. STIF/MASS MATRIX (3)
C     -------------------------------------------------------------
 150	IF (IEIG.NE.1) GOTO 500

C	THIS BLOCK FOR STANDARD EIGENVALUE IEIG=1
C	DO 160  IEQ=1,NEQ
C160	BB(IEQ) = 1.0
	FAC = 1.0D0
	CALL MDMAKE_DIA(TYP2,FAC,IA(LMA))  !MAKE DIAGONAL B MATRIX

C     -------------------------------------------------------------
C     SOLVE CURRENT EIGENVALUE PROBLEM BY SUBSPACE ITERATION METHOD
C     -------------------------------------------------------------
 500	CALL SSPACE (IA(LMA),A(LER),A(LEV),A(LTT),A(LWV),
     1             A(LAR),A(LBR),A(LVE),A(LDD),A(LRT),A(LBU),
     2             A(LBC),NEQ,NC,IA(LID),NSF,NSN,AA,BB,TYP1,TYP2)

C
 900  CALL CPU_TIME (TIM2)
	TIM(18) = TIM(18) + (TIM2-TIM1)
      ITASK = 1
C
 3100 FORMAT(//1X,79(1HX)//14X,15HS T A N D A R D
     1    ,40H   E I G E N V A L U E   A N A L Y S I S/14X,55(1H-)/)
 3200 FORMAT(//1X,79(1HX)//14X,15HB U C K L I N G
     1    ,40H   E I G E N V A L U E   A N A L Y S I S/14X,55(1H-)/)
 3300 FORMAT(//1X,79(1HX)//14X,17HV I B R A T I O N
     1    ,38H   F R E Q U E N C Y   A N A L Y S I S/14X,55(1H-)/)
 3500 FORMAT(
     +    14X,'INTERVAL BETWEEN EIGENVALUE ANALYSIS . . NSEIG =',I7/
     +    14X,'NUMBER OF EIGENVALUES REQUIRED . . . . . NROOT =',I7/
     +    14X,'NUMBER OF EIGENVECTORS REQUIRED. . . . . IVPRT =',I7/
     +    14X,'NUMBER OF ITERATION VECTORS USED . . . . .  NC =',I7/
     +    14X,'NUMBER OF SUBSPACE ITERATIONS PERMITTED  NITEM =',I7/
     +    14X,'FLAG FOR STURM SEQUENCE CHECK  . . . . .  IFSS =',I7/
     +    14X,'CONVERGENCE TOLERANCE ON EIGENVALUES . . . EPS =',E12.4/
     +    14X,'INITIAL SHIFT TO BE APPLIED . . . . . . SHIFT0 =',E12.4/
     +    /1X,79(1HX))
9999	FORMAT(//1X,'THIS OPTION IS FOR NONLINEAR EIGHE VALUE ANALYSIS'/
     +         1X,'THIS IS OPTION IS NOT AVAILABLE NOW')
C
 9990 RETURN
      END
C
C=====================================================================
      SUBROUTINE SSPACE (MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,BUP,
     1                   BUPC,MEQ,MC,ID,MF,MJ,AA,BB,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C	09/02/01 SUNIL IMPLEMENTED IN-CORE SUB-SPACE ITERATION METHOD
C	MODIFIED AND CHANGED THIS SUBROUTINE
C	----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
      COMMON A(9000000),IA(9000000)
C	-----------------------------------------------------------------
      DIMENSION MAXA(1),R(MEQ,1),TT(MEQ),W(MEQ),EIGV(1)
      DIMENSION D(MC),VEC(MC,MC),AR(*),BR(*),RTOLV(MC),BUP(MC)
      DIMENSION BUPC(MC),DISP(6),ID(MF,1)

	DIMENSION AA(1),BB(1)

C     --------------
C     INITIALISATION
C     --------------
      ICONV  = 0
      INDPD  = KPOSD
      IFEIG  = 1
      TOLJ   = 1.0E-12
      NSMAX  = 12
      N1     = NC+1
      NC1    = NC-1
      CALL CLEARA (D,NC)
      CALL CLEARA (R,NC*NEQ)

C     -------------------------------------------
C     APPLY SHIFT IF AA IS NOT POSITIVE DEFINITE
C     OR IN ORDER TO INCREASE RATE OF CONVERGENCE
C     -------------------------------------------
 50   IF (SHIFT0.EQ.0.0) GOTO 60
      CALL SHIFTA (MAXA,AA,BB,A(LDK),SHIFT0,NEI,0,TYP1,TYP2)

C     ----------------------------
C     STARTING  ITERATION  VECTORS
C     ----------------------------
 60   CALL INVECT (MAXA,R,TT,W,NEQ,AA,BB,TYP1,TYP2)

C     --------------------------
C     L*D*LT DECOMPOSITION OF AA
C     --------------------------
	CALL COLSOL (MAXA,AA,A(LDK),TT,1,INDPD,TYP1,'TEMP')

C     -----------------------
C     START OF ITERATION LOOP
C     -----------------------
      NITE = 0
 100  NITE = NITE+1

      IF (IFPR(10).EQ.1) WRITE (ISO,7000) NITE
C	----------------------------------------
C     CALCULATE THE PROJECTIONS OF AA  AND  BB
C	----------------------------------------
 110  IJ = 0
      DO 160  J=1,NC
      DO 120  K=1,NEQ
 120  TT(K) = R(K,J)
	CALL COLSOL (MAXA,AA,A(LDK),TT,2,INDPD,'TEMP','TEMP')
      DO 140  I=J,NC
      ART = 0.0
      DO 130  K=1,NEQ
 130  ART = ART + R(K,I)*TT(K)
      IJ  = IJ+1
 140  AR(IJ) = ART
      DO 150  K=1,NEQ
 150  R(K,J) = TT(K)
 160  CONTINUE

      IJ = 0
      DO 200  J=1,NC
	CALL MAMULT (MAXA,BB,R(1,J),TT,TYP2,'STD')
      DO 180  I=J,NC
      BRT = 0.0
      DO 170  K=1,NEQ
 170  BRT = BRT + R(K,I)*TT(K)
      IJ  = IJ+1
 180  BR(IJ) = BRT
      IF (ICONV.GT.0) GOTO 200
      DO 190  K=1,NEQ
 190  R(K,J) = TT(K)
 200  CONTINUE

C
      IF (IFPR(10).NE.1) GOTO 220
      CALL MATOUT (AR,NC,NC,6,11,'E',4,4,3,'A-PROJ(AR)')
      CALL MATOUT (BR,NC,NC,6,11,'E',4,4,3,'B-PROJ(BR)')
C     -------------------------------------------
C     SOLVE FOR EIGENSYSTEM OF SUBSPACE OPERATORS
C     -------------------------------------------
 220  CALL JACOBI (AR,BR,VEC,EIGV,W,NC,NNC,TOLJ,NSMAX)
      IF (IFPR(10).NE.1) GOTO 250
      WRITE (ISO,8200)
      CALL MATOUT (AR,NC,NC,6,11,'E',4,4,3,'A-PROJ(AR)')
      CALL MATOUT (BR,NC,NC,6,11,'E',4,4,3,'B-PROJ(BR)')
C     --------------------------------------
C     ARRANGE EIGENVALUES IN ASCENDING ORDER
C     --------------------------------------
 250  IS = 0
      II = 1
      DO 280  I=1,NC1
      ITEMP = II+N1-I
      IF (DABS(EIGV(I+1)).GE.DABS(EIGV(I))) GOTO 280
      IS = IS+1
      SWAP      = EIGV(I+1)
      EIGV(I+1) = EIGV(I)
      EIGV(I)   = SWAP
      SWAP      = BR(ITEMP)
      BR(ITEMP) = BR(II)
      BR(II)    = SWAP
      DO 270  K=1,NC
      SWAP      = VEC(K,I+1)
      VEC(K,I+1)= VEC(K,I)
 270  VEC(K,I)  = SWAP
 280  II = ITEMP
      IF (IS.GT.0) GOTO 250
C
      IF (IFPR(10).NE.1) GOTO 300
      WRITE (ISO,8500)
      WRITE (ISO,8600) (EIGV(I), I=1,NC)
C     -------------------------------------------------------
C     CALCULATE B TIMES APPROXIMATE EIGENVECTORS (ICONV.EQ.0)
C     OR FINAL EIGENVECTOR APPROXIMATIONS        (ICONV.GT.0)
C     -------------------------------------------------------
 300  DO 390  I=1,NEQ
      DO 310  J=1,NC
 310  TT(J) = R(I,J)
      DO 350  K=1,NC
      RT = 0.0
      DO 330  L=1,NC
 330  RT = RT + TT(L)*VEC(L,K)
 350  R(I,K) = RT
 390  CONTINUE
      IF (ICONV.GT.0) GOTO 500
C     ------------------------------------
C     CHECK FOR CONVERGENCE OF EIGENVALUES
C     ------------------------------------
      DO 400 I=1,NC
      DIF      = EIGV(I)-D(I)
 400  RTOLV(I) = DABS(DIF/EIGV(I))
      IF (IFPR(10).NE.1) GOTO 420
      WRITE (ISO,9000)
      WRITE (ISO,9100) (I,I,RTOLV(I), I=1,NC)
C
 420  DO 430  I=1,NROOT
 430  IF (RTOLV(I).GT.EPS) GOTO 450
      ICONV = 1
      GOTO 100
C
 450  IF (NITE.LT.NITEM) GOTO 460
      ICONV = 2
      IFSS  = 0
      GOTO 100
C
 460  DO 470  I=1,NC
 470  D(I) = EIGV(I)
      GOTO 100
C     ---------------------
C     END OF ITERATION LOOP
C     ---------------------
 500  CONTINUE
 
 
C     =================================================================
      DO 540  J=1,NROOT ! NUMBER OF MODE SHAPE
      EIGVAL = EIGV(J)-SHIFT0

	CALL RELFILL('EIVV',EIGVAL,1,J,1) !STORE EIGEN VALUE

	IF (IVPRT.EQ.0 .OR. (J-IVPRT).GT.0) GOTO 540

	CALL NORM (R(1,J),R(1,J),NEQ)

	CALL MODOUT(ID,R(1,J),J)          !STORE MODE SHAPE TO FILE

540	CONTINUE

C     =================================================================

C     -----------------------------------
C     APPLY STURM SEQUENCE CHECK (IFSS=0)
C     -----------------------------------
	IF (IFSS.EQ.0) GOTO 700
      CALL SCHECK (EIGV,RTOLV,BUPC,D,BUP,NC,NEI,EPS,SHIFT)
 700	IF (SHIFT0.NE.0.0) CALL SHIFTA (MAXA,AA,BB,A(LDK),-SHIFT,NEI,1,
	1								TYP1,TYP2)
C
	IFEIG = 0
C
 1000 FORMAT (//24X,32(1H*)/24X,1H*,30X,1H*/
     1        24X,9H* LOWEST ,I2,21H EIGENVALUE(S) AND  */
     2        24X,32H* CORRESPONDING EIGENVECTOR(S) */
     3        24X,1H*,30X,1H*/24X,32(1H*))
1002  FORMAT(//8X,'E I G E N S O L U T I O N   B Y  SUBSPACE ITERATION',
     *  3X,'M E T H O D'/8X,61(1H-)//
     *  19X,'O U T P U T   OF   E I G E N V A L U E S'/19X,40(1H-)//
     *  14X,'REQUIRED NUMBER OF EIGEN VALUES. . . . . . =',I5/
     +  14X,'SIZE OF SUBSPACE . . . . . . . . . . . . . =',I5//
     *  14X,'MODE NUMBER',27X,'EIGENVALUE',/14X,11(1H-),27X,10(1H-))
1222  FORMAT(//8X,'E I G E N S O L U T I O N   B Y  SUBSPACE ITERATION',
     *  3X,'M E T H O D'/8X,61(1H-)//
     *  19X,'O U T P U T   OF   E I G E N V A L U E S'/19X,40(1H-)//
     *  14X,'REQUIRED NUMBER OF EIGEN VALUES. . . . . . =',I5/
     +  14X,'SIZE OF SUBSPACE . . . . . . . . . . . . . =',I5//
     *  4X,'MODE NUMBER',10X,'EIGENVALUE',10X,'NATURAL FREQUENZY (Hz)',
     *  8X,'NATURAL PERIOD (Sec)',/
     *  4X,11(1H-),10X,10(1H-),10X,22(1H-),8X,20(1H-)) !SONGSAK NOV2005 NATURAL FREQUENZY
1001  FORMAT(//8X,'E I G E N S O L U T I O N   B Y  SUBSPACE ITERATION',
     *  3X,'M E T H O D'/8X,61(1H-)//)
3002  format(//14X,'MODE NUMBER',27X,
     *     'EIGENVALUE',/14X,11(1H-),27X,10(1H-))
 2000 FORMAT (//19X,'CONVERGENCE REACHED FOR EPS  =',E12.4/
     1          19X,'NUMBER OF ITERATIONS USED    =',I7)
 2100 FORMAT (//19X,18HNO CONVERGENCE FOR,I3,21H PERMITTED ITERATIONS//
     1          19X,34HWE ACCEPT CURRENT ITERATION VALUES//
     2          19X,41HTHE STURM SEQUENCE CHECK IS NOT PERFORMED)
 3000 FORMAT (//21X,16HEIGENVALUE LAMDA,I3,3H  =,E14.6/21X,36(1H-))
 3001  FORMAT(15X,I5,24X,E19.6)
 3100 FORMAT (//26X,'EIGENVECTOR OF MODE ='I5/26X,26(1H-)//
     1          1X,5HJOINT,7X,6HX-DISP,6X,6HY-DISP,6X,6HZ-DISP,
     2          6X,6HX-ROTA,6X,6HY-ROTA,6X,6HZ-ROTA/1X,79(1H-))
 3111  FORMAT(5X,I5,7X,E19.6,6X,E20.6,9X,E20.6)
 3200 FORMAT (I5,3X,6E12.4)
 4000 FORMAT (//25X,30HERROR NORMS ON THE EIGENVALUES/25X,30(1H-)/)
 4100 FORMAT (21X,9HFOR LAMDA,I3,4X,2HD(,I3,3H) =,E14.6)
 7000 FORMAT (//22X,32HI T E R A T I O N    N U M B E R,I4//
     1        22X,36(1H*)//)
 8200 FORMAT (//21X,38HAR AND BR AFTER JACOBI DIAGONALISATION/
     1        21X,38(1H-))
 8500 FORMAT (//,27X,26HEIGENVALUES OF AR-LAMDA*BR/27X,26(1H-)//)
 8600 FORMAT (7X,3E22.14)
 9000 FORMAT (//19X,41HRELATIVE TOLERANCE REACHED ON EIGENVALUES/
     1        19X,41(1H-)//)
 9100 FORMAT (19X,9HFOR LAMDA,I3,7X,6HRTOLV(,I2,4H) = ,E12.4/)
C
C	NEXT FIVE LINES ADDED BY GILSON - NOV2003 - FOR GIDOUT
c8010	FORMAT ('Displacement     4     'I4'.0     Estatic     2   1   1')
8010	FORMAT ('Result "ModeShapes"  "ModeShapes"  'I4' Vector  OnNodes')
8110	FORMAT ('ComponentNames  "X-Component", "Y-Component",',
     #' "Z-Component"')
8210	FORMAT ('Values')
8310	FORMAT ('End Values',/)

8016	FORMAT ('Result "ModeShapesMat" "ModeShapes" 'I4' Matrix OnNodes')
8116	FORMAT ('ComponentNames  "X", "Y", "Z", "RX", "RY", "RZ" ')
c
	RETURN
      END
C
C=====================================================================
      SUBROUTINE SHIFTA (MAXA,AA,BB,DD,SHIFT,NEI,IFACT,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C	09/02/01 SUNIL IMPLEMENTED IN-CORE SUB-SPACE ITERATION METHOD
C	MODIFIED AND CHANGED THIS SUBROUTINE
C     ----------------------------------------------------------------
C     APPLIES A SHIFT TO MATRIX AA, [AA] = [AA] + SHIFT*[BB]
C     FOR IFACT>0 INVERTS SHIFTED MATRIX AND COUNTS NO.OF NEG.PIVOTS
C	--------------------------------------------------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     MAXA(NEQ1)    = DIAGONAL ADDRESSES OF GLOBAL STIFFNESS
C     DD(NEQ)       = DIAGONAL COEF. OF INVERTED MATRIX [AA]
C     AA(ISTOR)     = GLOBAL COMPACTED STIFFNESS MATRIX
C     BB(ISTOR)     = GLOBAL COMPACTED GEOM/MASS MATRIX
C     SHIFT         = SHIFT TO BE APPLIED
C     NEI           = NUMBER OF CALCULATED EIGENVALUES
C     IFACT (0,1)   = FLAG TO SHIFT ONLY OR SHIFT AND INVERT
C	--------------------------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/,/INOU/ AND /SOLU/
C	--------------------------------------------------
C     NEQ           = NUMBER OF EQUATIONS
C     ISTOR         = MAXIMUM LENGTH OF ANY ONE BLOCK
C	---------------
C     LOCAL VARIABLES
C	---------------
C     IEQ           = EQUATION COUNTER
C     NPIV,NMIS     = NO.OF NEG.PIVOTS, NO.OF MISSING EIGENVALUES
C     ----------------------------------------------------------------
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /LOCA/ LID,LDS,LEL,LDC,LXY,LCH,LNU,LMP,LGP,LMS,LGS,
     1              LCO,LEX,LLM,LES,LEC,LED,LEI,LEE,LMA,LLF,LLV,
     2              LRE,LDI,LDL,LDT,LDK,LER,LEV,LTT,LWV,LAR,LBR,
     3              LVE,LDD,LRT,LBU,LBC,LVL,LAL,LEF,LDU,LPR,LLO,
	4              LRV,LRT1,LRET,LRET1,LDM,LDPT,LVL1,LMV,LXI,LCM,LCC,
	5			    LCN,LDIM,LFRE,LSFC,LLOF
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C	-----------------------------------------------------------------
      DIMENSION MAXA(1),DD(1),AA(1),BB(1)

C     -------------------------------------------------
C     APPLY SHIFT TO [A] FOR CONSISTENT MASS MATRIX [B]
C     -------------------------------------------------

C	DO 290  I=1,ISTOR
C 290  AA(I) = AA(I) + SHIFT*BB(I)
	CALL MDOPER(TYP1,TYP1,TYP2,1.0D0,SHIFT,'ADD')

C     -----------------------------------------------------
C     INVERT SHIFTED MATRIX AND COUNT NUMBER OF NEG. PIVOTS
C     -----------------------------------------------------
400	IF (IFACT.LE.0) RETURN

	CALL COLSOL (MAXA,AA,DD,DD,1,2,TYP1,'TEMP')

      NPIV = 0
      DO 600  IEQ=1,NEQ
 600  IF (DD(IEQ).LT.0.) NPIV = NPIV+1
      NMIS = NPIV-NEI
      WRITE (ISO,5000) SHIFT
      IF (NPIV.EQ.NEI) WRITE (ISO,5200) NPIV
      IF (NPIV.NE.NEI) WRITE (ISO,5100) NMIS
C
 5000 FORMAT (/////22X,35HCHECK IF NO EIGENVALUES ARE MISSING/
     1        22X,35(1H-)///19X,22HCHECK APPLIED AT SHIFT,E20.12/)
 5100 FORMAT (24X,9HTHERE ARE,I3,20H EIGENVALUES MISSING)
 5200 FORMAT (23X,19HWE FOUND THE LOWEST,I3,12H EIGENVALUES)
C
      RETURN
      END
C
C=====================================================================
      SUBROUTINE INVECT (MAXA,R,TT,W,MEQ,AA,BB,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C	09/02/01 SUNIL IMPLEMENTED IN-CORE SUB-SPACE ITERATION METHOD
C	MODIFIED AND CHANGED THIS SUBROUTINE
C     ----------------------------------------------------------------
C     FINDS STARTING ITERATION VECTORS FROM DIAGONAL TERMS OF
C     STIFFNESS MATRIX
C	----------------------
C     INPUT,OUTPUT VARIABLES
C	----------------------
C     MAXA(NEQ1)    = DIAGONAL ADDRESSES OF COMPACTED GLOBAL STIFF/MASS
C     AA(ISTOR)     = GLOBAL COMPACTED STIFFNESS BLOCK
C     BB(ISTOR)     = GLOBAL COMPACTED GEOM/MASS BLOCK
C     R(NEQ,1)      = STARTING ITERATION VECTORS
C     TT(NEQ)       = EXCITED DEGREES OF FREEDOM BY UNIT STARTING VEC.
C     W(NEQ)        = RATIOS OF DIAGONAL TERMS BB(KDIA)/AA(KDIA)
C	--------------------------------------------------
C     VARIABLES IN COMMON BLOCK /ELEM/,/INOU/ AND /SOLU/
C	--------------------------------------------------
C     NEQ           = NUMBER OF EQUATIONS
C     NWM           = NUMBER OF COEFFICIENTS IN GEOM/MASS MATRIX
C     ISTOR         = MAXIMUM LENGTH OF ANY ONE BLOCK
C     NROOT         = NUMBER OF REQUIRED EIGENVALUES
C     NC            = NUMBER OF VECTORS USED IN SUBSPACE ITER.
C	---------------
C     LOCAL VARIABLES
C	---------------
C     NMASS         = NUMBER OF APPLIED NODAL MASSES
C     IEQ           = EQUATION COUNTER
C     KDIA          = ADDRESS OF DIAGONAL TERM IN COMPACTED BLOCK
C     WMAX          = ABSOLUTE MAXIMUM VALUE IN W
C     ----------------------------------------------------------------
      COMMON /ELEM/ NAME(2),ITYPE,ISTYP,NLOPT,MTMOD,NSINC,ITOLEY,
     1              NELE,NMPS,NGPS,NMP,NGP,NNM,NEX,NCO,NNF,NWG,NEFC,
     2              NPT,NWA,NWS,KEG,MEL,NNO,NEF,NELTOT,NMV,MTYP,ISECT
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C	-----------------------------------------------------------------
C
      DIMENSION MAXA(1),R(MEQ,1),TT(1),W(1),AA(1),BB(1)
C     --------------
C     INITIALISATION
C     --------------
      NMASS = 0
C     ---------------------------------------------
C     CONSISTENT MASS OR GEOMETRIC STIFFNESS MATRIX
C     ---------------------------------------------

C	DO 290  IEQ=1,NEQ
C	KDIA = MAXA(IEQ)
C	R(IEQ,1) = BB(KDIA)
C290  W(IEQ)   = BB(KDIA)/AA(KDIA)
	CALL MDCALL_DIA(TYP2,R(1,1),'RED',MAXA)
	CALL MDCALL_DIA(TYP1,W     ,'RED',MAXA)
	DO 290  IEQ=1,NEQ
 290  W(IEQ) = R(IEQ,1)/W(IEQ)

C     -----------------------------------------
C     EXCITE STARTING ITERATION VECTORS 2 TO NC
C     -----------------------------------------
 300	ND = NEQ/NC
      LEQ = NEQ-ND
      DO 400  J=2,NC
      WMAX = 0.
      DO 330 IEQ=1,LEQ
      IF (DABS(W(IEQ)).LT.WMAX) GOTO 330
      WMAX = DABS(W(IEQ))
      IPOS = IEQ
 330  CONTINUE
      DO 350  IEQ=LEQ,NEQ
      IF (DABS(W(IEQ)).LE.WMAX) GOTO 350
      WMAX = DABS(W(IEQ))
      IPOS = IEQ
 350  CONTINUE
      TT(J) = FLOAT(IPOS)
      W(IPOS) = 0.
      LEQ     = LEQ-ND
 400  R(IPOS,J) = 1.0
C
      IF (IFPR(10).NE.1) RETURN
      WRITE (ISO,1000) NROOT
      WRITE (ISO,6000)
      WRITE (ISO,6100) (TT(J), J=2,NC)
C
 1000 FORMAT (///24X,32(1H*)/24X,1H*,30X,1H*/
     1        24X,9H* LOWEST ,I2,21H EIGENVALUE(S) AND  */
     2         24X,32H* CORRESPONDING EIGENVECTOR(S) */
     3         24X,1H*,30X,1H*/24X,32(1H*)////)
 6000 FORMAT (10X,30HDEGREES OF FREEDOM EXCITED BY
     1        ,30HUNIT STARTING ITERATION VECTOR//)
 6100 FORMAT (8X,10F6.0,/)
C
      RETURN
      END
C
C=====================================================================

	SUBROUTINE JACOBI (A,B,X,EIGV,D,N,NWA,RTOL,NSMAX)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ----------------------------------------------------------------
C     PROGRAM TO SOLVE THE GENERALIZED EIGENPROBLEM USING
C     THE GENERALIZED JACOBI ITERATION
C	--------------------------------
C     INPUT VARIABLES
C	---------------
C     A(NWA)    = STIFFNESS MATRIX (ASSUMED POSITIVE DEFINITE)
C     B(NWM)    = MASS MATRIX (ASSUMED POSITIVE DEFINITE)
C     X(N,N)    = MATRIX STORING EIGENVECTOR APPROXIMATIONS
C     EIGV(N)   = VECTOR STORING EIGENVALUES (CURRENT APPROXIMATION)
C     D(N)      = WORKING VECTOR TO STORE PREVIOUS EIGENV. APPROX.
C     N         = ORDER OF MATRICES A AND B
C     NWM       = NUMBER OF ELEMENTS IN A AND B (UPPER TRIANGULAR FORM)
C     RTOL      = CONVERGENCE TOLERANCE (USUALY SET TO 10.**-12)
C     NSMAX     = MAXIMUM NUMBER OF SWEEPS ALLOWED (USUALLY 15)
C     IFPR(10)  = FLAG FOR PRINTING DURING ITERATION (1 FOR PRINTING)
C	----------------
C     OUTPUT VARIABLES
C	----------------
C     A(NWM)    = DIAGONALIZED STIFFNESS MATRIX
C     B(NWM)    = DIAGONALIZED MASS MATRIX
C     X(N,N)    = EIGENVECTORS STORED COLUMNWISE
C     EIGV(N)   = EIGENVALUES
C     ----------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
C      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,STOL,NFAC,
C     +              NRED,KPOSD,DETK,DET1,DAVR
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
C	----------------------------------------------------------------
      DIMENSION A(NWA),B(NWA),X(N,N),EIGV(N),D(N)

C     ----------------------------------------------
C     INITIALIZE EIGENVALUE AND EIGENVECTOR MATRICES
C     ----------------------------------------------
      N1 = N+1
      II = 1
      DO 10  I=1,N
	IF (A(II).LE.0.0 .AND. KPOSD.EQ.0)
	1	CALL ERRORS (12,II,A(II),'JACOBI ITE')
	IF (DABS(B(II)).EQ.0.) CALL ERRORS (12,II,B(II),'JACOBI ITE')
      D(I) = A(II)/B(II)
      EIGV(I) = D(I)
 10   II = II+N1-I
      DO 30  I=1,N
      DO 20  J=1,N
 20   X(I,J) = 0.0
 30   X(I,I) = 1.0
      IF (N .EQ. 1) RETURN
C     --------------------------------------------
C     INITIALIZE SWEEP COUNTER AND BEGIN ITERATION
C     --------------------------------------------
      NSWEEP = 0
      NR     = N-1
 40   NSWEEP = NSWEEP+1
      IF (IFPR(10).EQ.1) WRITE (ISO,1000) NSWEEP
C     -----------------------------------------
C     CHECK IF PRESENT OFF-DIAGONAL  ELEMENT IS
C     LARGE ENOUGH TO REQUIRE ZEROING
C     -----------------------------------------
      EPS = (0.01**NSWEEP)**2
      DO 210  J=1,NR
      JP1 = J+1
      JM1 = J-1
      LJK = JM1*N - JM1*J/2
      JJ  = LJK+J
      DO 210  K=JP1,N
      KP1 = K+1
      KM1 = K-1
      JK  = LJK+K
      KK  = KM1*N - KM1*K/2 + K
      EPTOLA = (A(JK)*A(JK))/(A(JJ)*A(KK))
      EPTOLB = (B(JK)*B(JK))/(B(JJ)*B(KK))
      IF ((EPTOLA.LT.EPS) .AND. (EPTOLB.LT.EPS)) GOTO 210
C     ----------------------------------------------
C     IF ZEROING IS REQUIRED, CALCULATE THE ROTATION
C     MATRIX ELEMENTS CA AND CG
C     ----------------------------------------------
      AKK = A(KK)*B(JK)-B(KK)*A(JK)
      AJJ = A(JJ)*B(JK)-B(JJ)*A(JK)
      AB  = A(JJ)*B(KK)-A(KK)*B(JJ)
      CHECK = (AB*AB+4.*AKK*AJJ)/4.
      IF (CHECK)  50,60,60
 50   CALL ERRORS (12,J,CHECK,'JACOBI ITE')
 60   SQCH = DSQRT(CHECK)
      D1 = AB/2.+SQCH
      D2 = AB/2.-SQCH
      DEN = D1
      IF (DABS(D2) .GT. DABS(D1)) DEN = D2
      IF (DEN) 80,70,80
 70   CA = 0.0
      CG = -A(JK)/A(KK)
      GOTO 90
 80   CA = AKK/DEN
      CG = -AJJ/DEN
C     ----------------------------------------
C     PERFORM THE GENERALIZED ROTATION
C     TO ZERO THE PRESENT OFF-DIAGONAL ELEMENT
C     ----------------------------------------
 90   IF (N-2)   100,190,100
 100  IF (JM1-1) 130,110,110
 110  DO 120  I=1,JM1
      IM1 = I-1
      IJ  = IM1*N - IM1*I/2 + J
      IK  = IM1*N - IM1*I/2 + K
      AJ  = A(IJ)
      BJ  = B(IJ)
      AK  = A(IK)
      BK  = B(IK)
      A(IJ) = AJ+CG*AK
      B(IJ) = BJ+CG*BK
      A(IK) = AK+CA*AJ
 120  B(IK) = BK+CA*BJ
 130  IF (KP1-N) 140,140,160
 140  LJI = JM1*N - JM1*J/2
      LKI = KM1*N - KM1*K/2
      DO 150  I=KP1,N
      JI  = LJI+I
      KI  = LKI+I
      AJ  = A(JI)
      BJ  = B(JI)
      AK  = A(KI)
      BK  = B(KI)
      A(JI) = AJ+CG*AK
      B(JI) = BJ+CG*BK
      A(KI) = AK+CA*AJ
 150  B(KI) = BK+CA*BJ
 160  IF (JP1-KM1)  170,170,190
 170  LJI = JM1*N - JM1*J/2
      DO 180  I=JP1,KM1
      JI  = LJI+I
      IM1 = I-1
      IK  = IM1*N - IM1*I/2 + K
      AJ  = A(JI)
      BJ  = B(JI)
      AK  = A(IK)
      BK  = B(IK)
      A(JI) = AJ+CG*AK
      B(JI) = BJ+CG*BK
      A(IK) = AK+CA*AJ
 180  B(IK) = BK+CA*BJ
 190  AK  = A(KK)
      BK  = B(KK)
      A(KK) = AK + 2.0*CA*A(JK) + CA*CA*A(JJ)
      B(KK) = BK + 2.0*CA*B(JK) + CA*CA*B(JJ)
      A(JJ) = A(JJ) + 2.0*CG*A(JK) + CG*CG*AK
      B(JJ) = B(JJ) + 2.0*CG*B(JK) + CG*CG*BK
      A(JK) = 0.0
      B(JK) = 0.0
C     -------------------------------------------------
C     UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION
C     -------------------------------------------------
      DO 200  I=1,N
      XJ = X(I,J)
      XK = X(I,K)
      X(I,J) = XJ+CG*XK
 200  X(I,K) = XK+CA*XJ
 210  CONTINUE
C     ---------------------------------------
C     UPDATE THE EIGENVALUES AFTER EACH SWEEP
C     ---------------------------------------
      II = 1
      DO 220  I=1,N
      IF (A(II).LE.0.0 .AND. KPOSD.EQ.0)
     1      CALL ERRORS (12,II,A(II),'JACOBI ITE')
      IF (DABS(B(II)).EQ.0.) CALL ERRORS (12,II,B(II),'JACOBI ITE')
      EIGV(I) = A(II)/B(II)
 220  II = II+N1-I
      IF (IFPR(10).EQ.1) WRITE (ISO,1100) (EIGV(I),I=1,N)
C     -----------------------
C     CHECK  FOR  CONVERGENCE
C     -----------------------
      DO 240  I=1,N
      DIF = EIGV(I)-D(I)
 240  IF (DABS(DIF/D(I)).GT.RTOL) GOTO 280
C     --------------------------------------
C     CHECK ALL OFF-DIAGONAL ELEMENTS TO SEE
C     IF ANATHOR SWEEP IS REQUIRED
C     --------------------------------------
      EPS = RTOL*RTOL
      DO 250  J=1,NR
      JM1 = J-1
      JP1 = J+1
      LJK = JM1*N - JM1*J/2
      JJ  = LJK+J
      DO 250  K=JP1,N
      KM1 = K-1
      JK  = LJK+K
      KK  = KM1*N - KM1*K/2 + K
      EPSA = (A(JK)*A(JK))/(A(JJ)*A(KK))
      EPSB = (B(JK)*B(JK))/(B(JJ)*B(KK))
      IF ((EPSA.LT.EPS) .AND. (EPSB.LT.EPS)) GOTO 250
      GOTO 280
 250  CONTINUE
C     ------------------
C     SCALE EIGENVECTORS
C     ------------------
 255  II = 1
      DO 275  I=1,N
      BB = DSQRT(DABS(B(II)))
      DO 270  K=1,N
 270  X(K,I) = X(K,I)/BB
 275  II = II+N1-I
      RETURN
C     -----------------------------------------------
C     UPDATE D MATRIX AND START NEW SWEEP, IF ALLOWED
C     -----------------------------------------------
 280  DO 290  I=1,N
 290  D(I) = EIGV(I)
      IF (NSWEEP .LT. NSMAX) GOTO 40
      GOTO 255
C
 1000 FORMAT (///,10X,36HCURRENT EIGENVALUES FOR SWEEP NUMBER,I4,
     1        10H IN JACOBI/10X,50(1H-)/)
 1100 FORMAT (8X,3E20.12)
      END
C
C=====================================================================
      SUBROUTINE SCHECK (EIGV,RTOLV,BUPC,NEIV,NEIL,MC,NEI,RTOL,SHIFT)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     --------------------------------------------
C     EVALUATION OF SHIFT FOR STURM SEQUENCE CHECK
C     --------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /EIGN/ NSEIG,NROOT,NC,NNC,NITEM,IFSS,SHIFT0,EPS,IEIG,NEIG,
     +              ISOLV,IVPRT
C	---------------------------------------------
      DIMENSION EIGV(MC),RTOLV(MC),BUPC(MC),NEIV(MC),NEIL(MC)
C
      MROOT = 0
      DO 100  I=1,NC
  100 IF (RTOLV(I) .LT. RTOL) MROOT = MROOT+1
C     ----------------------------------------
C     FIND UPPER BOUNDS ON EIGENVALUE CLUSTERS
C     ----------------------------------------
      DO 240  I=1,MROOT
 240  NEIV(I) = 1
      L       = 1
      I       = 2
      IF (MROOT .NE. 1) GOTO 270
      BUPC(1) = 1.01*EIGV(1)
      GOTO 300
 270  IF (1.01*EIGV(I-1) .LE. 0.99*EIGV(I)) GOTO 280
      NEIV(L) = NEIV(L)+1
      GOTO 285
 280  BUPC(L) = 1.01*EIGV(I-1)
      L = L+1
 285  I = I+1
      IF (I .LE. MROOT) GOTO 270
      BUPC(L) = 1.01*EIGV(I-1)
C     -----------
C     FIND  SHIFT
C     -----------
 300  NEIL(1) = NEIV(1)
      IF (L .EQ. 1) GOTO 400
      DO 320  IC=2,L
 320  NEIL(IC) = NEIV(IC)+NEIL(IC-1)
C
 400  SHIFT = BUPC(1)
      NEI   = NEIL(1)
      DO 490  IC=1,L
      IF (NEIL(IC).GE.NROOT) GOTO 900
      SHIFT = BUPC(IC+1)
 490  NEI   = NEIL(IC+1)
C
C     IF (IFPR(10).NE.1)  RETURN
 900  WRITE (ISO,1000)
      WRITE (ISO,1100) (I,BUPC(I),NEIV(I),NEIL(I), I=1,L)
C
 1000 FORMAT (//29X,19HEIGENVALUE CLUSTERS/29X,19(1H-)//
     1        3X,7HCLUSTER,5X,11HUPPER BOUND,5X,17HNO OF EIGENVALUES
     2        ,5X,24HNO LESS THAN UPPER BOUND/3X,76(1H-)//)
 1100 FORMAT (6X,I2,3X,E19.12,12X,I2,20X,I2/)
C
      RETURN
      END
C
C=====================================================================
C	END OF SUBROUTINE FOR SUBSPACE EIGEN SOLVER
C=====================================================================

C=====================================================================
C	SUBROUTINE MASS PARTICPATION FACTOR
C=====================================================================
      SUBROUTINE MPFCAL(MAXA,ID,RVEC,REIG,DIMR,RFRE,NMOD,NITEM,AM,TYP,KSC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP
C     ---------------------------------------------------------------
C     PROGRAM TO 
C		+ COMPUTE TOTAL STRUCTURE'S MASS PARTICIPATED IN EACH
C			   TRANSLATIONAL DIRECTIONS, 
C		+ CALCULATE MODAL MASS PARTICIPATION FACTORS
C	added by Nguyen, December 2009 
C	---------------------------------------------------------------
C	INPUT ARGUMENT:
C	---------------
C		MAXA = ARRAY CONTAIN GLOBAL DIAGONAL ADDRESSES
C		AM   = GLOBAL MASS MATRIX
C		ID(NSF,NSN) = IDENTITY ARRAY (EQUATION NUMBERS)
C		NMOD = NUMBER OF MODES CONSIDERED
C		RVEC(NEQ,NMOD)= A(LER) = REDUCED MATRIX OF EIGENVECTOR TO SUPERPOSE
C		REIG(NMOD)    = REDUCED EIGENVALUE ARRAY,CORRESPONDING TO [RVEC]
C	----------------
C	OUTPUT ARGUMENTS:
C	----------------
C		XMP(NMOD) = MODAL MASS PARTICIPATION FACTORS IN X-DIRECTION
C		YMP(NMOD) = MODAL MASS PARTICIPATION FACTORS IN Y-DIRECTION
C		ZMP(NMOD) = MODAL MASS PARTICIPATION FACTORS IN Z-DIRECTION
C		DIMR(NMOD)= Diagonal elements array of MODAL MASS MATRIX
C		RFRE(NMOD)= CIRCULAR FREQUENCY ARRAY
C     ---------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC  
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
     	COMMON /LINEAT/ KTRAF,KEATH,KCSAL,KOFFL,KSPEC,KDESIGN,KFATM,KFATJ,KFATL,KFAST,KOREV !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1 !SONGSAK AUG2007 RESPONSE SPECTRUM FOR ISOLOP 1
C	------------------------------------------------------------
	DIMENSION ID(NSF,NSN),MAXA(1),AM(NWM),RVEC(NEQ,NMOD),REIG(NMOD)
	DIMENSION COX(NEQ),COY(NEQ),COZ(NEQ),CRX(NEQ),CRY(NEQ),CRZ(NEQ)	
	DIMENSION AX(NEQ),AY(NEQ),AZ(NEQ),ARX(NEQ),ARY(NEQ),ARZ(NEQ)
	DIMENSION XMP(NMOD),YMP(NMOD),ZMP(NMOD),XRP(NMOD),YRP(NMOD)
	DIMENSION ZRP(NMOD),DIMR(NMOD),RFRE(NMOD)
C	------------------------------------------------------------
C	Initilization
	CALL CLEARA(COX,NEQ)
	CALL CLEARA(COY,NEQ)
	CALL CLEARA(COZ,NEQ)
	CALL CLEARA(CRX,NEQ)
	CALL CLEARA(CRY,NEQ)
	CALL CLEARA(CRZ,NEQ)

	CALL CLEARA(AX,NEQ)
	CALL CLEARA(AY,NEQ)
	CALL CLEARA(AZ,NEQ)
	CALL CLEARA(ARX,NEQ)
	CALL CLEARA(ARY,NEQ)
	CALL CLEARA(ARZ,NEQ)

	CALL CLEARA(XMP,NMOD)
	CALL CLEARA(YMP,NMOD)
	CALL CLEARA(ZMP,NMOD)
	CALL CLEARA(XRP,NMOD)
	CALL CLEARA(YRP,NMOD)
	CALL CLEARA(ZRP,NMOD)

	CALL CLEARA(DIMR,NMOD)

	DO 100 J = 1,NSN
	  IXE = ID(1,J)  !Equation number of X-dof, node J
	  IF(IXE.EQ.0) GOTO 150
		 COX(IXE) = 1.0

  150	  IYE = ID(2,J)	 !Equation number of Y-dof, node J
	  IF(IYE.EQ.0) GOTO 200
		COY(IYE) = 1.0

  200	  IZE = ID(3,J)	 !Equation number of Z-dof, node J
	  IF(IZE.EQ.0) GOTO 250
		COZ(IZE) = 1.0	
		
  250	  IRX = ID(4,J)	 !Equation number of XROT-dof, node J
	  IF(IRX.EQ.0) GOTO 270
		CRX(IRX) = 1.0	 

  270	  IRY = ID(5,J)	 !Equation number of YROT-dof, node J
	  IF(IRY.EQ.0) GOTO 280
		CRY(IRY) = 1.0	
		
  280	  IRZ = ID(6,J)	 !Equation number of ZROT-dof, node J
	  IF(IRZ.EQ.0) GOTO 100
		CRZ(IRZ) = 1.0					
				   
  100 CONTINUE

c	Form vectors {AX} = [AM]*{COX};{AY} = [AM]*{COY};{AZ} = [AM]*{COZ}
c  	CALL MAMULT(MAXA,AM,COX,AX,NWM)  !changed next 08Mar2010 by Nguyen
      CALL MAMULT(MAXA,AM,COX,AX,TYP,'STD')

c  	CALL MAMULT(MAXA,AM,COY,AY,NWM)
      CALL MAMULT(MAXA,AM,COY,AY,TYP,'STD')

c  	CALL MAMULT(MAXA,AM,COZ,AZ,NWM)
      CALL MAMULT(MAXA,AM,COZ,AZ,TYP,'STD')

c  	CALL MAMULT(MAXA,AM,CRX,ARX,NWM)
  	CALL MAMULT(MAXA,AM,CRX,ARX,TYP,'STD')

c  	CALL MAMULT(MAXA,AM,CRY,ARY,NWM)
      CALL MAMULT(MAXA,AM,CRY,ARY,TYP,'STD')

c  	CALL MAMULT(MAXA,AM,CRZ,ARZ,NWM)
      CALL MAMULT(MAXA,AM,CRZ,ARZ,TYP,'STD')

c	Form total mass participated in X,Y,Z & RX,RY,RZ directions
C	NORM PRODUCT OF TWO VECTORS {COX_T}*{AX}
C	 VNORM = SUM(XT(I)*AX(I),I=1,NEQ) = (XT(1)*AX(1)+...+(XT(N)*AX(N)
      XMAS = 0.
	YMAS = 0.
	ZMAS = 0.

      XRMAS = 0.
	YRMAS = 0.
	ZRMAS = 0

      DO 300  I=1,NEQ
 		XMAS = XMAS + COX(I)*AX(I)
 		YMAS = YMAS + COY(I)*AY(I)
 		ZMAS = ZMAS + COZ(I)*AZ(I)

 		XRMAS = XRMAS + CRX(I)*ARX(I)
 		YRMAS = YRMAS + CRY(I)*ARY(I)
 		ZRMAS = ZRMAS + CRZ(I)*ARZ(I)
  300 CONTINUE

C	FORM MODAL MASS MATRIX [RMAS(NMOD,NMOD)] = [RVEC]T *[MASS]*[RVEC]
C	-->DIMR(NMOD) = Diagonal elements array of product matrix [RMAS]
c	CALL VFVMULD (MAXA,AM,RVEC,DIMR,NEQ,NMOD,NWM) !!changed next 08Mar2010 by Nguyen
      CALL VFVMULD (MAXA,AM,RVEC,DIMR,NEQ,NMOD,'MASS')
      
C     ================================================
C                   MODIFILED BY TOEY
C     ================================================
      IF (KSPEC.EQ.1.AND.KSC.EQ.1)THEN
      CALL STOREMODESHAPE (RVEC,NMOD)
C      CALL MUTIFORCE (RVEC,NMOD)
      ENDIF
      
C	CALCULATE MODAL MASS PARTICIPATION FACTORS
	DO 400 IM = 1,NMOD
        XMOD = 0.
	  YMOD = 0.
	  ZMOD = 0.

        XRMOD = 0.
	  YRMOD = 0.
	  ZRMOD = 0.

        DO 450  IE=1,NEQ
 		XMOD = XMOD + RVEC(IE,IM)*AX(IE)
 		YMOD = YMOD + RVEC(IE,IM)*AY(IE)
 		ZMOD = ZMOD + RVEC(IE,IM)*AZ(IE)

 		XRMOD = XRMOD + RVEC(IE,IM)*ARX(IE)
 		YRMOD = YRMOD + RVEC(IE,IM)*ARY(IE)
 		ZRMOD = ZRMOD + RVEC(IE,IM)*ARZ(IE)
  450   CONTINUE
	  
C	  XMP(IM) = (XMOD**2)/XMAS
C	  YMP(IM) = (YMOD**2)/YMAS
C	  ZMP(IM) = (ZMOD**2)/ZMAS

C	  XMP(IM) = ABS(XMOD/XMAS)
C	  YMP(IM) = ABS(YMOD/YMAS)
C	  ZMP(IM) = ABS(ZMOD/ZMAS)

c	  XMP(IM) = 100*(XMOD**2)/(XMAS*DIMR(IM))
c	  YMP(IM) = 100*(YMOD**2)/(YMAS*DIMR(IM))
c	  ZMP(IM) = 100*(ZMOD**2)/(ZMAS*DIMR(IM))

c	  XRP(IM) = 100*(XRMOD**2)/(XRMAS*DIMR(IM))
c	  YRP(IM) = 100*(YRMOD**2)/(YRMAS*DIMR(IM))
c	  ZRP(IM) = 100*(ZRMOD**2)/(ZRMAS*DIMR(IM))

c	added IF 08Jan2010 when one direction DOFs are fixed (ex. plane struc)

	  IF(XMAS.NE.0.0) THEN 
		XMP(IM) = 100*(XMOD**2)/(XMAS*DIMR(IM))
	  ELSE
          XMP(IM) = 0.
	  ENDIF

	  IF(YMAS.NE.0.0) THEN 
		YMP(IM) = 100*(YMOD**2)/(YMAS*DIMR(IM))
	  ELSE
		YMP(IM) = 0.
	  ENDIF
	  
	  IF(ZMAS.NE.0.0) THEN 	  	  
		ZMP(IM) = 100*(ZMOD**2)/(ZMAS*DIMR(IM))
	  ELSE
		ZMP(IM) = 0.
	  ENDIF

	  IF(XRMAS.NE.0.0) THEN 	
		XRP(IM) = 100*(XRMOD**2)/(XRMAS*DIMR(IM))
	  ELSE
		XRP(IM) = 0.0
	  ENDIF

	  IF(YRMAS.NE.0.0) THEN 
		YRP(IM) = 100*(YRMOD**2)/(YRMAS*DIMR(IM))
	  ELSE
		YRP(IM) = 0.
	  ENDIF

	  IF(ZRMAS.NE.0.0) THEN 
	    ZRP(IM) = 100*(ZRMOD**2)/(ZRMAS*DIMR(IM))
	  ELSE
		ZRP(IM) = 0.
	  ENDIF

  400 CONTINUE

C	CALCULATE CUMULATIVE MASS PARTICIPATION OF NMOD-MODES CONSIDERED
	CALL SUMVEC(XMP,NMOD,SX)
	CALL SUMVEC(YMP,NMOD,SY)
	CALL SUMVEC(ZMP,NMOD,SZ)

	CALL SUMVEC(XRP,NMOD,SRX)
	CALL SUMVEC(YRP,NMOD,SRY)
	CALL SUMVEC(ZRP,NMOD,SRZ)
c
	CALL FREPRT(REIG,XMP,YMP,ZMP,XRP,YRP,ZRP,RFRE,NMOD,NITEM,KSC)

      IF (KSC.NE.1)THEN
	WRITE(100,1000) SX,SY,SZ,SRX,SRY,SRZ
	WRITE(ITO,1000) SX,SY,SZ,SRX,SRY,SRZ
      WRITE(10,1000) SX,SY,SZ,SRX,SRY,SRZ
c 1000 FORMAT (/,6x,'Cumulative Sum of MPF:',6(2X,E10.4))
 1000 FORMAT (/,3x,'Cumulative Sum of MPF(%):',6(2X,f8.4))
      ENDIF
	 
      RETURN
      END
C

C=====================================================================
C	SUBROUTINE MASS PARTICPATION FACTOR
C=====================================================================
      SUBROUTINE FREPRT(EIGV,XMP,YMP,ZMP,XRP,YRP,ZRP,ROOT,NMOD,NITEM,KSC)     
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	--------------------------------------------------------------------------
C     PROGRAM TO 
C		+ CALCULATE NATURAL FREQUENCIES AND PERIODS FROM EIGENVALUES 
C		+ PRINT DYNAMIC PROPERTIES STRUCTURE (FREQ,PERIODS)
C	added December 2009 by Nguyen, used for any Eigen Solver
C	---------------------------------------------------------------------------
C	INPUT:
C	EIGV    = EIGENVALUE ARRAY
C     NMOD    = NUMBER OF EIGENVALUES REQUIRED
C     NITEM   = NUMBER OF SUBSPACE ITERATIONS PERMITTED
C				IF NITEM.LE.0, TWO EIGENVALUES SHOULD BE INPUTTED
C	OUTPUT:
C	ROOT(NMOD)= CIRCULAR FREQUENCY ARRAY
C	-------------------------------------------------------------------------
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON / EIGVPED / EIGVFREQ(1000),EIGVPERI(1000)
      DIMENSION EIGV(NMOD),ROOT(NMOD),XMP(NMOD),YMP(NMOD),ZMP(NMOD)
	DIMENSION XRP(NMOD),YRP(NMOD),ZRP(NMOD)  
C	--------------------------------------------------------------------------
      PI=3.141592654
      IF (KSC.NE.1)THEN
      WRITE(100,2010)
      WRITE(ITO,2010)
      WRITE(10,2010)
      ENDIF

      IF(NITEM.LE.0) THEN
C	READ INPUT CIRCULAR FREQUENCIES OF BRIDGE
	  READ(ITI,*)
        READ(ITI,1000) (ROOT(I),I=1,NMOD) 

	  DO 20 I=1,NROOT
		EIGV=ROOT(I)**2
		FREQ=ROOT(I)/(2.*PI)
		PERI=1./FREQ
		WRITE(ISO,2020) I,PERI,FREQ,ROOT(I),EIGV
   20	  CONTINUE
        GOTO 500
      ENDIF

C	CALCULATE NATURAL FREQUENCIES AND PERIODS
      DO 40 I=1,NMOD
        ROOT(I)=DSQRT(EIGV(I))
        FREQ=ROOT(I)/(2.*PI)
        PERI=1./FREQ
        EIGVFREQ(I) = FREQ
        EIGVPERI(I) = PERI
        
        IF (KSC.NE.1)THEN
c       WRITE(100,2020) I,EIGV(I),FREQ,PERI,XMP(I),YMP(I),ZMP(I),XRP(I),
c	1				  YRP(I),ZRP(I) !changed next, not to print eigenvalues
        WRITE(100,2020) I,FREQ,PERI,XMP(I),YMP(I),ZMP(I),XRP(I),YRP(I),
	1				  ZRP(I)
        WRITE(ITO,2020) I,FREQ,PERI,XMP(I),YMP(I),ZMP(I),XRP(I),YRP(I),
	1				  ZRP(I)
        WRITE(10,2020) I,FREQ,PERI,XMP(I),YMP(I),ZMP(I),XRP(I),YRP(I),
	1				  ZRP(I)
	  ENDIF
   40 CONTINUE

 1000 FORMAT(2F10.0)
 2010 FORMAT(//,1X,'MODE',2X,'FREQUENCY',4X,'PERIOD',16X    
     1,'MASS PARTICIPATION FACTORS (MPF,%)',/
     2,8X,'(Hz)',8X,'(sec)',8X,'X-DIR',4X,'Y-DIR',
     3 5X,'Z-DIR',5X,'X-ROT',5X,'Y-ROT',5X,'Z-ROT',/
     4,1X,3(1H-),3X,9(1H-),3X,9(1H-),4X,56(1H-))

c 2020 FORMAT(I5,3(2X,E13.6),6(2X,E10.4)) !changed next, not to print eigenvalues
 2020 FORMAT(I3,1X,2(2X,E11.4),6(2X,f8.4))
       
  500	RETURN
      END
C   

C=====================================================================
C	SUBROUTINE MASS PARTICPATION FACTOR
C=====================================================================
C	NEXT PART FOR MASS PARTICIPATION FACTOR
C	added by Nguyen DV, December 2009 & March 2010
C	----------------------------------------------------------------
      SUBROUTINE SUMVEC(VEC,N,SUM)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C     ---------------------------------------------------------------
C     PROGRAM TO SUM ALL COMPONENTS OF VECTOR
C	---------------------------------------------------------------
C	VEC(N) =  VECTOR
C	
C	OUTPUT:
C	-------
C	SUM
C	---------------------------------------------------------------
	DIMENSION VEC(1)

	SUM = 0.0
      DO 100  I=1,N
  100		SUM = SUM + VEC(I)

      RETURN
      END
C
C=====================================================================
C	SUBROUTINE MASS PARTICPATION FACTOR
C=====================================================================


