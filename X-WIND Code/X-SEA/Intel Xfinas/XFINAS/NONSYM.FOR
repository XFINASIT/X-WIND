C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE COLSOLN (MAXA,A,D,V,IND,INDPD,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2
C     ----------------------------------------------------------------
C     SOLVE FINITE ELEMENT STATIC EQUILIBRIUM EQUATIONS USING IN CORE
C	SOLVER BASED ON COMPACTED STORAGE AND COLUMN REDUCTION SCHEME
C	-------------------------------------------------------------
C     INPUT VARIABLES IN ARGUMENT LIST
C	--------------------------------
C     MAXA(NEQ1)    = ADDRESSES OF DIAGONAL ELEMENTS IN A
C     A(ISTOR)      = STIFFNESS BLOCK TO BE REDUCED
C     B(ISTOR)      = COUPLING BLOCK
C     D(NEQ)        = DIAGONAL ELEMENTS OF GLOBAL STIFFNESS
C     V(NEQ)        = RIGHT-HAND-SIDE LOAD VECTOR
C     IND=1         = FLAG FOR FACTORISATION OF STIFFNESS MATRIX OR
C     IND=2         = REDUCTION AND BACK-SUBSTITUTION OF LOAD VECTOR
C     INDPD         = FLAG FOR FACTORISATION OF NON-POSITIVE DEFINITE A
C	---------------------------------
C     OUTPUT VARIABLES IN ARGUMENT LIST
C	---------------------------------
C     A(ISTOR)      = L AND FACTORS OF STIFFNESS BLOCK
C     D(NEQ)        = D FACTORS OF GLOBAL STIFFNESS MATRIX
C     V(NEQ)        = VECTOR OF GENERALIZED DISPLACEMENTS
C	-------------------------------------------
C     VARIABLES IN COMMON BLOCK /SOLU/ AND /INOU/
C	-------------------------------------------
C     NEQ           = NUMBER OF EQUATIONS (NEQ1=NEQ+1)
C     NBLOCK        = NUMBER OF STIFFNESS BLOCKS
C     ISTOR         = MAXIMUM LENGTH OF ANY ONE BLOCK
C     DETK          = DETERMINANT OF GLOBAL STIFFNESS MATRIX
C     NKFAC         = UNIT NUMBER OF RANDOM ACCESS FILE STORING
C                     FACTORIZED STIFFNESS BLOCKS
C     LSIGN         = STABILITY INDICATER (+1=STABEL, -1=INSTABEL)
C	---------------
C     LOCAL VARIABLES
C	---------------
C     NCOLUM        = NUMBER OF COLUMNS FOR THIS BLOCK
C     IBLO,JBLO     = LOOP COUNTER FOR REDUCTION/COUPLING BLOCK
C     ICOL,JCOL     = LOOP COUNTER FOR REDUCTION/COUPLING COLUMN
C     KBLO,KCOL     = BLOCK/COLUMN COUNTER
C     JBF,JBL       = FIRST,LAST COUPLED BLOCK
C     JCF,JCL       = FIRST,LAST COUPLED COLUMN
C     KDIA,KLOW,KUPP= ADDRESSES OF DIAGONAL,FIRST AND LAST ELEMENT IN
C                     COLUMN TO BE REDUCED
C     KHEI,KMOD,LCOP= TOTAL,EFFECTIVE HEIGHT AND NUM.OF COUPLED ELEM.
C     KTOP          = TOP ELEMENT TO BE REDUCED
C     LDIA,LHEI     = DIAGONAL ELEMENT AND HEIGHT OF COUPLING COLUMN
C     IEQ           = EQUATION NUMBER
C     -----------------------------------------------------------------
      COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     +              NRED,KPOSD,DETK,DET1,DAVR,STOL
      COMMON /ITER/ RHO,RHOP,RHOPREV,RTOL,ETOL,DLMAX,ALP,
	1              NSTEP,NPRIN,NDRAW,
	2			  KONEQ,NIREF,ITOPT,ICONV,NOLIN,KSTEP,
     3              LIMEQ(2),ITEMAX,NUMREF,NUMITE,ITETOT,LIMET
      COMMON /INOU/ ITI,ITO,ISO,NDATI,NPLOT,NKFAC,NELEM,
     1              IFPR(10),IFPL(10)
      COMMON /FTIM/ TIM(20),IDATE,ITIME
      COMMON /FLAG/ IFPRI,ISPRI,IFPLO,IFREF,IFEIG,ITASK,IFFLAG
C
	DIMENSION MAXA(1),A(1),D(1),V(1)
	ALLOCATABLE NCLM(:),ICOP(:),IHIG(:),B(:),AL(:),BL(:)

C	-------------------------------------------
C	IMOPT = 0 SKYLINE + COLSOL , IMOPT = 1 SPARSE  + ITERATION
	CALL MINTFIL('SOLV',IMOPT,1,1,0)		!STIFFNESS PROFILE AND SOLVER OPTION
	CALL MINTFIL('SOLV',ITPRE,1,2,0)		!PRECONDITIONER SCHEME
	IF(IMOPT.EQ.1) THEN
	IF(ITPRE.EQ.0) CALL CONJUGS(MAXA,A,V,NEQ,IND,TYP1,TYP2)		!STANDARD
	IF(ITPRE.EQ.1) CALL CONJUGJ(MAXA,A,V,NEQ,IND,TYP1,TYP2)		!JACOBI
	IF(ITPRE.EQ.2) CALL CONJUGI(MAXA,A,V,NEQ,IND,TYP1,TYP2,D)	!IMCOMPLETE CHOLESKI
	RETURN
	ENDIF
C	-------------------------------------------

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)


	ALLOCATE(NCLM(NBLOCK),ICOP(NBLOCK),IHIG(NBLOCK))

	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCLM(IBLO),1,IBLO,0)
	CALL MINTFIL('BLOC',ICOP(IBLO),2,IBLO,0)
	CALL MINTFIL('BLOC',IHIG(IBLO),3,IBLO,0)
	ENDDO


	CALL CPU_TIME (TIM1)

      IF (IND-2) 100,800,800

100	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)
	CALL MCALFIL(KTREC,TYP2)
	REWIND(KTREC)

	ALLOCATE(B(MSTOR),AL(MSTOR),BL(MSTOR))
	
	CALL CPU_TIME (TIME1)
C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	NCOLB = 0

	DO 1000 IBLO=1,NBLOCK

	NCOL = NCLM(IBLO)
      NPRE = MAXA(NCOLB+1)-1
	NHIG = IHIG(IBLO)

	READ(KSREC)  A(1:NHIG)
	READ(KSREC) AL(1:NHIG)

C     --------------------------------------------------------
C     FACTORIZE STIFFNESS MATRIX A (A = L*D*L'T DECOMPOSISION)
C     --------------------------------------------------------
      IF (IBLO.EQ.ICOP(IBLO))  GOTO 400

C     ------------------------------------------------
C     FIND SUM (MCOLB) OF PRECEEDING UNCOUPLED COLUMNS
C     FIRST (JBF) AND LAST (JBL) COUPLED BLOCK
C     ------------------------------------------------
      JBF = ICOP(IBLO)-1
      MCOLB = 0
      IF (JBF)  400,150,120
 120  DO 130  JBLO=1,JBF
 130  MCOLB = MCOLB+NCLM(JBLO)
 150  MCOL  = NCOLB-MCOLB
      JBF = JBF+1
      JBL = IBLO-1
C
C     ----------------------------------------------
C     REDUCE BLOCK BY THE PRECEEDING COUPLING BLOCKS
C     LOOP OVER NUMBER OF COUPLING BLOCKS
C     READ BLOCK FROM RANDOM ACESS TAPE
C     ----------------------------------------------
	REWIND(KTREC)
	DO II = 1,JBF-1
	READ(KTREC) 
	READ(KTREC) 
	ENDDO

      DO 390  JBLO=JBF,JBL
	NHIG = IHIG(JBLO)

	READ(KTREC)  B(1:NHIG)
	READ(KTREC) BL(1:NHIG)

      MCOL = MCOL-NCLM(JBLO)
      MPRE = MAXA(MCOLB+1)-1
C     ------------------------------------------
C     LOOP OVER NUMBER OF COLUMNS OF BLOCK TO BE
C     REDUCED,ADDRESSES OF COLUMN ELEMENTS
C     ------------------------------------------
      DO 295  ICOL=1,NCOL
      KDIA = MAXA(NCOLB+ICOL)-NPRE
      KLOW = KDIA+1
      KUPP = MAXA(NCOLB+ICOL+1)-1-NPRE
      KHEI = KUPP-KLOW-ICOL+1
      KMOD = KHEI-MCOL
      IF (KMOD.LE.0)  GOTO 295
C
C     -------------------------------------------------
C     LOOP OVER NUMBER OF COUPLING COLUMNS (JCF TO JCL)
C     ADDRESSES OF ELEMENTS IN COUPLING COLUMN
C     -------------------------------------------------
      LCOP = 0
      JCF  = NCLM(JBLO)-KMOD+1
      JCL  = NCLM(JBLO)
      IF (JCF.GT.0)  GOTO 210
      LCOP = KMOD-NCLM(JBLO)
      JCF  = 1
 210  KTOP = KUPP-LCOP
C
      DO 290  JCOL=JCF,JCL
      LCOP = LCOP+1
      KTOP = KTOP-1
      LDIA = MAXA(MCOLB+JCOL)-MPRE
      LHEI = MAXA(MCOLB+JCOL+1)-1-LDIA-MPRE
      IF (LHEI)  290,290,230
 230  NL = MIN0(LCOP,LHEI)

      DO 250  IL=1,NL
	AL(KTOP) = AL(KTOP) -  B(LDIA+IL)*AL(KTOP+IL)
 250   A(KTOP) =  A(KTOP) - BL(LDIA+IL)* A(KTOP+IL)

 290  CONTINUE
 295  CONTINUE
C
 390  MCOLB = MCOLB+NCLM(JBLO)

C     -----------------------------------------
C     REDUCE BLOCK, LOOP OVER NUMBER OF COLUMNS
C     -----------------------------------------
 400  DO 690  ICOL=1,NCOL
      KDIA = MAXA(NCOLB+ICOL)-NPRE
      KLOW = KDIA+1
      KUPP = MAXA(NCOLB+ICOL+1)-1-NPRE
      KCOL = ICOL + NCOLB
      KHEI = KUPP-KLOW
      KMOD = MIN0(KHEI,ICOL-1)

      IF (KMOD) 600,500,410
C     ---------------------------------------------
C     LOOP OVER NUMBER OF COUPLING COLUMNS JCF, JCL
C     ---------------------------------------------
 410  LCOP = 0
      JCF  = ICOL-KMOD
      JCL  = ICOL-1
      KTOP = KLOW+KMOD
      IF (ICOL-1.LT.KHEI) LCOP = KHEI-ICOL+1	
      DO 490 JCOL=JCF,JCL
      LCOP = LCOP+1
      KTOP = KTOP-1
      LDIA = MAXA(NCOLB+JCOL)-NPRE
      LHEI = MAXA(NCOLB+JCOL+1)-1-LDIA-NPRE
      IF (LHEI) 490,490,430
 430  NL = MIN0(LCOP,LHEI)

      DO 450  IL=1,NL
	AL(KTOP) = AL(KTOP) -  A(LDIA+IL)*AL(KTOP+IL)
 450	 A(KTOP) =  A(KTOP) - AL(LDIA+IL)* A(KTOP+IL)

 490  CONTINUE

C     -------------------------------------------
C     FINAL COLUMN TERMS LIJ = GIJ/DII  AND
C     DIAGONAL TERMS     DJJ = KJJ - SUM(IRJ*GRJ)
C     -------------------------------------------
 500	IEQ = KCOL
      SUM = 0.0D0
      DO 590  KA=KLOW,KUPP
      IEQ = IEQ-1
	SUM = SUM + A(KA)*AL(KA)/D(IEQ)
	AL(KA) = AL(KA)/D(IEQ)
 590   A(KA) =  A(KA)/D(IEQ)
	AL(KDIA) = AL(KDIA) - SUM
	 A(KDIA) =  A(KDIA) - SUM
C     ----------------------------------------------------------
C     SET DIAGONAL TERMS IN D AND TEST WHETHER POSITIVE DEFINITE
C     ----------------------------------------------------------
 600  D(KCOL) = A(KDIA)

C	KK = MAXA(KCOL+1) - MAXA(KCOL)
C	WRITE(*,*) KCOL,A(KDIA)

      PIV = A(KDIA)
      IF (PIV-STOL)       610,610,690
 610  IF (DABS(PIV)-STOL) 650,650,620
 620  IF (INDPD-1)        630,690,690
 630  CALL ERRORS (12,KCOL,PIV,' SOLUTION ')
 650  IF (INDPD-1)        630,630,660
 660  PIV = STOL
      IF (PIV.EQ.0.0) PIV = -1.0D-16
      D(KCOL) = PIV
      A(KDIA) = PIV
      AL(KDIA) = PIV
 690  CONTINUE
C     ------------------------------------------------
C     END OF LOOP
C     CALCULATE DETERMINANT OF STIFFNESS MATRIX (DETK)
C     ------------------------------------------------
	NHIG = IHIG(IBLO)
	WRITE(KTREC)  A(1:NHIG)
	WRITE(KTREC) AL(1:NHIG)

C	DO IH = 1,NHIG
C	WRITE(*,*) IBLO,A(IH),AL(IH)
C	ENDDO

C     ------------------------------------------------
1000  NCOLB = NCOLB + NCOL  !END BLOCK LOOP
C     ------------------------------------------------
C	PAUSE

	DEALLOCATE(NCLM,ICOP,IHIG)

	DEALLOCATE(B,AL,BL)

	CALL CPU_TIME (TIME2)
C	WRITE(*,6000) NFAC+1,TIME2-TIME1

 700	IF (ITASK-2) 750,750,795
C
 750  DETK=1.0
      POWER = 1.0/DFLOAT(NEQ)
      DO 720 IEQ=1,NEQ
      DETK=DETK*DABS(D(IEQ))**POWER
 720  IF (D(IEQ).LT. 0.0) DETK = -DETK

 790  NFAC = NFAC+1
C
	CALL CPU_TIME (TIM2)
      TIM(14) = TIM(14) + (TIM2-TIM1)
 795  RETURN




 800	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)


	ALLOCATE(AL(MSTOR))
C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	NCOLB = 0

	DO 2000 IBLO=1,NBLOCK

	NCOL = NCLM(IBLO)
      NPRE = MAXA(NCOLB+1)-1

	NHIG = IHIG(IBLO)
	READ(KSREC)  A(1:NHIG)
	READ(KSREC) AL(1:NHIG)
C     --------------------------------
C     FORWARD REDUCTION OF LOAD VECTOR
C     --------------------------------


C     ------------------------------------------------
C     LOOP OVER NUMBER OF COLUMNS, VI = RI-SUM(IRI*VR)
C     ------------------------------------------------
	DO 880  ICOL=1,NCOL
      KLOW = MAXA(NCOLB+ICOL)+1-NPRE
      KUPP = MAXA(NCOLB+ICOL+1)-1-NPRE
      IF (KUPP-KLOW) 880,830,830
 830  IEQ = NCOLB+ICOL
      KV = IEQ
      C = 0.0
      DO 850  KA=KLOW,KUPP
      KV = KV-1
 850  C = C + AL(KA)*V(KV)
      V(IEQ) = V(IEQ) - C
 880  CONTINUE

C     ------------------------------------------------
2000  NCOLB = NCOLB + NCOL  !END BLOCK LOOP
C     ------------------------------------------------

C     -----------------------------
C     BACK SUBSTITUTION , V=[D]-1*V
C     -----------------------------
 890	DO 900  IEQ=1,NEQ
 900  V(IEQ) = V(IEQ)/D(IEQ)


 905	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)
C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	KBLO = NBLOCK

	DO 3000 IBLO=1,NBLOCK
	REWIND(KSREC)

	NCOL  =  NCLM(KBLO)
	NCOLB = NCOLB - NCOL
      NPRE  = MAXA(NCOLB+1)-1

	DO II = 1,KBLO-1
	NHIG = IHIG(II)
	READ(KSREC) ! A(1:NHIG)
	READ(KSREC) !AL(1:NHIG)
	ENDDO
	NHIG = IHIG(KBLO)
	READ(KSREC)  A(1:NHIG)
	READ(KSREC) AL(1:NHIG)

C     ---------------------------
C     LOOP OVER NUMBER OF COLUMNS
C     ---------------------------
 910	KCOL = NCOL
      DO 980  ICOL=1,NCOL
      KLOW = MAXA(NCOLB+KCOL)+1-NPRE
      KUPP = MAXA(NCOLB+KCOL+1)-1-NPRE
      IF (KUPP-KLOW) 980,930,930
 930  IEQ = NCOLB+KCOL
      KV = IEQ
      DO 950  KA=KLOW,KUPP
      KV = KV-1
 950  V(KV) = V(KV) - A(KA)*V(IEQ)
 980  KCOL = KCOL-1

	KBLO = KBLO-1
C     ------------------------------------------------
3000  CONTINUE  !END BLOCK LOOP
C     ------------------------------------------------


	DEALLOCATE(AL)

	DEALLOCATE(NCLM,ICOP,IHIG)

C
 990  IF (ITASK.GT.2) RETURN
      NRED = NRED+1
C
	CALL CPU_TIME (TIM2)
      TIM(15) = TIM(15) + (TIM2-TIM1)


 6000	FORMAT (X,'FACTORIZATION NO. .',I6, '  TIME TO PERFORM. .',E15.6)


      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================








