C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE CSABLOK (MAXA,NEQ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C     ----------------------------------------------------------------
C	CALCULATE NUMBER OF BLOCK
C     ----------------------------------------------------------------
C
      DIMENSION MAXA(1)
	ALLOCATABLE MCOL(:),ICOP(:),IHIG(:)

	MXCOL = 0
	DO IEQ = 1,NEQ
	NCOL = MAXA(IEQ+1) - MAXA(IEQ)
	IF(NCOL.GT.MXCOL) MXCOL = NCOL
	ENDDO

C	NWK = MAXA(NEQ+1) -1
	MSTOR = 5*NEQ   !BLOCK SIZE HERE!!!!!!!!!!!!!!
	IF(MXCOL.GT.MSTOR) MSTOR = MXCOL
C	------------------------------------
	NCO    = 1
	NBLOCK = 0
	DO IEQ = 1,NEQ
	
	CALL BLOCKC(MAXA,NEQ,NCO,NCN,MSTOR)

	NBLOCK = NBLOCK + 1
	IF(NCN.EQ.NEQ+1) GOTO 1000

	NCO = NCN

	ENDDO
C	------------------------------------
1000	CONTINUE


	ALLOCATE(MCOL(NBLOCK),ICOP(NBLOCK),IHIG(NBLOCK))


C	-------------------------------------
C	NCOL
C	-------------------------------------
	NCO = 1
	DO IBLO=1,NBLOCK
	CALL BLOCKC(MAXA,NEQ,NCO,NCN,MSTOR)
	NCOL = NCN - NCO
	MCOL(IBLO) = NCOL
	IHIG(IBLO) = MAXA(NCN) - MAXA(NCO)
	NCO = NCN
	ENDDO
C	-------------------------------------
C	COUPLING ICOL
C	-------------------------------------
      ICOP(1) = 1
      NCOLB   = MCOL(1)
      DO 390  IBLO=2,NBLOCK
C
C     LOOP OVER NUMBER OF COLUMNS (NCOLUM) FOR THIS BLOCK (IBLO)
C     TO FIND MAXIMUM EXCESS COLUMN LENGTH (ICLMAX)
C
      ICLMAX = 0
      NCOLUM = MCOL(IBLO)
      DO 320  ICOL=1,NCOLUM
      ICL = MAXA(NCOLB+ICOL+1) - MAXA(NCOLB+ICOL) - ICOL - 1
 320  IF (ICL.GT.ICLMAX)  ICLMAX = ICL
C
      IBLOP = IBLO
 350  IF (ICLMAX.LE.0)  GOTO 380
      IBLOP = IBLOP-1
      ICLMAX = ICLMAX - MCOL(IBLOP)
      GOTO 350
 380  ICOP(IBLO) = IBLOP

 390  NCOLB = NCOLB + MCOL(IBLO)
C	-------------------------------------
C	-------------------------------------
C	CLEAR ALL MEMORY
	CALL DELMINT('BLOK')
	CALL  DEMINT('BLOK',KBLOK,1,10)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,1)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,1)
	CALL MINTFIL('BLOK',0     ,1,3 ,1)  !NUMBER OF LINE TO WRITE FOR STIFFNESS ASSEMBLY

	ICPRO = 0   !STIFNESS SIZE REDUCTION PROCESS IS NOT YET DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1 SEE ALSO SUB. CONTRO)
C	ICPRO = 1   !STIFNESS SIZE REDUCTION PROCESS IS ALREADY DONE SEE ALSO SUB. MCSPARS  (USE WITH IMOPT=1 SEE ALSO SUB. CONTRO)
	CALL MINTFIL('BLOK',ICPRO ,1,4 ,1)  !INDEX FOR STIFFNESS SIZE REDUCTION PROCESS

	IPFAC = 0   !IPFAC=IPFAG DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT   IN PARDISO
C	IPFAC = 1   !IPFAC=IPFAG DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT   IN PARDISO
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,1)  

C     PIVOT FLAG  0=POSITIVE DEF  1=NEGATIVE DEF  ONLY APPEARED IN PARDISO,NEWSOLVER,LOADUP
      IPIVT = 0 
      CALL MINTFIL('BLOK',IPIVT ,1,6 ,1)  
      
C	CLEAR ALL MEMORY
	CALL DELMINT('BLOC')
	CALL  DEMINT('BLOC',KBLOK,3,NBLOCK)


	DO IBLO=1,NBLOCK
	CALL MINTFIL('BLOC',MCOL(IBLO) ,1,IBLO,1)
	CALL MINTFIL('BLOC',ICOP(IBLO) ,2,IBLO,1)
	CALL MINTFIL('BLOC',IHIG(IBLO) ,3,IBLO,1)
	ENDDO	



	DEALLOCATE(MCOL,ICOP,IHIG)


	RETURN

	END
C	=============================================================
C	=============================================================
C	=============================================================

      SUBROUTINE XCOLSOL (MAXA,A,D,V,NEQ,IND,INDPD,STOL,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP1,TYP2

C	NEED ONLY LSYMM
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM
C     ----------------------------------------------------------------
C     SOLVE FINITE ELEMENT STATIC EQUILIBRIUM EQUATIONS USING IN CORE
C	SOLVER BASED ON COMPACTED STORAGE AND COLUMN REDUCTION SCHEME
C	----------------------------------------------------------------
	DIMENSION MAXA(1),A(1),D(1),V(1)
	ALLOCATABLE NCLM(:),ICOP(:),IHIG(:),B(:)

C	IF NON-SYMMETRIC
	IF(LSYMM.EQ.1) THEN
	CALL COLSOLN (MAXA,A,D,V,IND,INDPD,TYP1,TYP2)
	RETURN
	ENDIF

C	-------------------------------------------
C	IMOPT = 0 SKYLINE + COLSOL , IMOPT = 1 SPARSE  + ITERATION
	CALL MINTFIL('SOLV',IMOPT,1,1,0)		!STIFFNESS PROFILE AND SOLVER OPTION
	CALL MINTFIL('SOLV',ITPRE,1,2,0)		!PRECONDITIONER SCHEME
	IF(IMOPT.EQ.1) THEN
	IF(ITPRE.EQ.0) CALL CONJUGS(MAXA,A,V,NEQ,IND,TYP1,TYP2)		!STANDARD
	IF(ITPRE.EQ.1) CALL CONJUGJ(MAXA,A,V,NEQ,IND,TYP1,TYP2)		!JACOBI
	IF(ITPRE.EQ.2) CALL CONJUGI(MAXA,A,V,NEQ,IND,TYP1,TYP2,D)	    !IMCOMPLETE CHOLESKI
	IF(ITPRE.EQ.3) CALL PARDISI(MAXA,V,NEQ,IND,TYP1,TYP2)             !PARDISO
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


      IF (IND-2) 100,800,800

100	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)
	CALL MCALFIL(KTREC,TYP2)
	REWIND(KTREC)

	ALLOCATE(B(MSTOR))
	
	CALL CPU_TIME (TIME1)
C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	NCOLB = 0

	DO 1000 IBLO=1,NBLOCK

	NCOL = NCLM(IBLO)
      NPRE = MAXA(NCOLB+1)-1
	
	NHIG = IHIG(IBLO)
	READ(KSREC) A(1:NHIG)

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
	NHIG = IHIG(II)
	READ(KTREC) !B(1:NHIG)
	ENDDO

      DO 390  JBLO=JBF,JBL
	NHIG = IHIG(JBLO)
	READ(KTREC) B(1:NHIG)
      MCOL = MCOL-NCLM(JBLO)
      MPRE = MAXA(MCOLB+1)-1
C
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
      C = 0.
      DO 250  IL=1,NL
 250  C = C + B(LDIA+IL)*A(KTOP+IL)
      A(KTOP) = A(KTOP)-C
 290  CONTINUE
 295  CONTINUE
C
 390  MCOLB = MCOLB+NCLM(JBLO)
C

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
      C = 0.0D0
      DO 450  IL=1,NL
 450  C = C + A(LDIA+IL)*A(KTOP+IL)
      A(KTOP) = A(KTOP) - C
 490  CONTINUE

C     -------------------------------------------
C     FINAL COLUMN TERMS LIJ = GIJ/DII  AND
C     DIAGONAL TERMS     DJJ = KJJ - SUM(IRJ*GRJ)
C     -------------------------------------------
 500  IEQ = KCOL
      SUM = 0.0
      DO 590  KA=KLOW,KUPP
      IEQ = IEQ-1
      C = A(KA)/D(IEQ)
      SUM = SUM + C*A(KA)
 590  A(KA) = C
      A(KDIA) = A(KDIA) - SUM
C     ----------------------------------------------------------
C     SET DIAGONAL TERMS IN D AND TEST WHETHER POSITIVE DEFINITE
C     ----------------------------------------------------------
 600  D(KCOL) = A(KDIA)

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
 690  CONTINUE
C     ------------------------------------------------
C     END OF LOOP
C     CALCULATE DETERMINANT OF STIFFNESS MATRIX (DETK)
C     ------------------------------------------------

	NHIG = IHIG(IBLO)
	WRITE(KTREC) A(1:NHIG)

C     ------------------------------------------------
1000  NCOLB = NCOLB + NCOL  !END BLOCK LOOP
C     ------------------------------------------------
C	STOP

	DEALLOCATE(NCLM,ICOP,IHIG)
	DEALLOCATE(B)

	CALL CPU_TIME (TIME2)
C	WRITE(*,6000) NFAC+1,TIME2-TIME1

 790  RETURN



 800	CALL MCALFIL(KSREC,TYP1)
	REWIND(KSREC)

C     ----------------------------
C     INITIALISATION FOR BLOCK LOOP
C     ----------------------------
	NCOLB = 0

	DO 2000 IBLO=1,NBLOCK

	NCOL = NCLM(IBLO)
      NPRE = MAXA(NCOLB+1)-1

	NHIG = IHIG(IBLO)
	READ(KSREC) A(1:NHIG)

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
 850  C = C + A(KA)*V(KV)
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
	READ(KSREC) !A(1:NHIG)
	ENDDO
	NHIG = IHIG(KBLO)
	READ(KSREC) A(1:NHIG)
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

	DEALLOCATE(NCLM,ICOP,IHIG)
C
 990  RETURN


 6000	FORMAT (X,'FACTORIZATION NO. .',I6, '  TIME TO PERFORM. .',E15.6)


      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================