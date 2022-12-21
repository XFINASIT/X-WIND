C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE ADMEMY
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C     Define & Open sequential files      
C	=============================================================

	CALL DEMINT('MDSA',KMDSA,1,20)

	NFL  = 301
	CALL MINTFIL('MDSA',NFL,1,1 ,1)
	CALL SEQOPEN(NFL)  

	NFL  = 302
	CALL MINTFIL('MDSA',NFL,1,2 ,1)
	CALL SEQOPEN(NFL) 

	NFL  = 303
	CALL MINTFIL('MDSA',NFL,1,3 ,1)
	CALL SEQOPEN(NFL) 

	NFL = 304
	CALL MINTFIL('MDSA',NFL,1,4 ,1)
	CALL SEQOPEN(NFL) 

	NFL = 305
	CALL MINTFIL('MDSA',NFL,1,5 ,1)
	CALL SEQOPEN(NFL) 

	NFL = 306
	CALL MINTFIL('MDSA',NFL,1,6 ,1)
	CALL SEQOPEN(NFL) 

	NFL = 307
	CALL MINTFIL('MDSA',NFL,1,7 ,1)
	CALL SEQOPEN(NFL) 

	NFL = 308
	CALL MINTFIL('MDSA',NFL,1,8 ,1)
	CALL SEQOPEN(NFL) 

	NFL = 309
	CALL MINTFIL('MDSA',NFL,1,9 ,1)
	CALL SEQOPEN(NFL) 

	NFL = 310
	CALL MINTFIL('MDSA',NFL,1,10,1)
	CALL SEQOPEN(NFL) 

	NFL = 311
	CALL MINTFIL('MDSA',NFL,1,11,1)
	CALL SEQOPEN(NFL) 

	NFL = 312
	CALL MINTFIL('MDSA',NFL,1,12,1)
	CALL SEQOPEN(NFL) 

	NFL = 313
	CALL MINTFIL('MDSA',NFL,1,13,1)
	CALL SEQOPEN(NFL) 

	NFL = 314
	CALL MINTFIL('MDSA',NFL,1,14,1)
	CALL SEQOPEN(NFL) 

	NFL = 315
	CALL MINTFIL('MDSA',NFL,1,15,1)
	CALL SEQOPEN(NFL) 

	NFL = 316
	CALL MINTFIL('MDSA',NFL,1,16,1)
	CALL SEQOPEN(NFL) 


	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE MBLOCK (MAXA,NEQ)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C     ----------------------------------------------------------------
C	CALCULATE NUMBER OF BLOCK
C     ----------------------------------------------------------------
C
      DIMENSION MAXA(1)
	ALLOCATABLE MCOL(:),ICOP(:),IHIG(:)

C     Get Max column height      
	MXCOL = 0
	DO IEQ = 1,NEQ
	NCOL = MAXA(IEQ+1) - MAXA(IEQ)
	IF(NCOL.GT.MXCOL) MXCOL = NCOL
      ENDDO

C     Define max "block" storage size      
C	NWK = MAXA(NEQ+1) -1
	MSTOR = 5*NEQ   !BLOCK SIZE HERE!!!!!!!!!!!!!!
	IF(MXCOL.GT.MSTOR) MSTOR = MXCOL
      
      
C	------------------------------------
C     Calculate number of block to be used      
	NCO    = 1
	NBLOCK = 0
	DO IEQ = 1,NEQ

C     Determine how much memory to be stored into each block          
	CALL BLOCKC(MAXA,NEQ,NCO,NCN,MSTOR) 

	NBLOCK = NBLOCK + 1
	IF(NCN.EQ.NEQ+1) GOTO 1000

	NCO = NCN

	ENDDO
C	------------------------------------
1000	CONTINUE


	ALLOCATE(MCOL(NBLOCK),ICOP(NBLOCK),IHIG(NBLOCK))


C	-------------------------------------
C	NCOL / Check how many column which stored in each block
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


C     Create small memory for block control data (only 10 data)    
	CALL  DEMINT('BLOK',KBLOK,1,10)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,1)  !Number of block
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,1)  !Maximum storage size for one block
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

C     NEFM = MAXIMUM NUMBER OF ELEMENT NDOF (MAX FROM OVERALL GROUP)  .... UPDATE IN ELLOOP
      NEFM = 0 !INITIALIZE
      CALL MINTFIL('BLOK',NEFM  ,1,7 ,1)        
      
C     Allocate memory of properties for each block (3 data per block)      
	CALL DEMINT('BLOC',KBLOK,3,NBLOCK)

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
	SUBROUTINE ADIWOK(IEG,NELE,NWA,NSTH,LCS)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C     PREPARE FILE FOR WORKING ARRAY     
C	=============================================================
	ALLOCATABLE A(:),B(:)
C	=============================================================

	NBIT  = 2*NWA
	IF(NSTH.GT.NWA) NBIT  = 2*NSTH

	NFLB  = 3200          !BACK UP WORKING ARRAY FILE
	KWBAK = NFLB + IEG
	CALL DIROPEN(KWBAK,NBIT) 

	NFL0  = 3300          !NORMAL  WORKING ARRAY FILE
	KWREC = NFL0 + IEG
	CALL DIROPEN(KWREC,NBIT) 
      
	NFLN  = 3600          !NORMAL  WORKING ARRAY FILE BY BJ FOR SOLID-SHELL(P)
	KWRECN = NFLN + IEG
	CALL DIROPEN(KWRECN,NBIT)       


	CALL MINTFIL('GWOK',NELE     ,IEG  ,1,1) !Number of element in this group
	MELE = 0
      
	IF(IEG.GT.1) 
	1CALL MINTFIL('GWOK',MELE     ,IEG-1,2,0) !Read accumulate number of element from prvious group
      
	CALL MINTFIL('GWOK',NELE+MELE,IEG  ,2,1) !Accumulate number of element (combine with previous group)
	CALL MINTFIL('GWOK',NWA      ,IEG  ,3,1) !Size of working array
	CALL MINTFIL('GWOK',NSTH     ,IEG  ,4,1) !Size of thermal load / initial stress

C	INITIALIZED THE POINTER FOR ELEMENT GROUP FOR TEMPERATURE LOAD (0 = NO TEMP LOAD ON THIS ELEMENT,1 = TEMP LOAD ON THIS ELEMENT)
C	SEE ALSO THERMALLOAD.FOR 
	NSTHG = 0
	CALL MINTFIL('GWOK',NSTHG    ,IEG  ,5,1)   


	ALLOCATE(A(NWA),B(NSTH))

	A(1:NWA ) = 0.0D0
	B(1:NSTH) = 0.0D0
	DO IELE = 1,NELE

C	FRIST NELE BLOCK FOR WORKING ARRAY
	IRC = IELE 
	WRITE(KWREC,REC=IRC)  A(1:NWA)
      
C	FIRST NELE BLOCK FOR NEW WORKING ARRAY BY BJ FOR SOLID-SHELL(P) 
	IRC = IELE 
	WRITE(KWRECN,REC=IRC)  A(1:NWA)      
      
      GOTO 10
      
C	FOLLOWING BLOCK FOR THERMAL LOAD INITIAL STRESS FOR EACH LOADCASE
	DO ILC = 1,LCS                    !LOOP OVER LOAD CASE
	IRC0= IELE + NELE + 2*NELE*(ILC-1)
	WRITE(KWREC,REC=IRC0) B(1:NSTH)   !FOR INITIAL STRESS FOR VARY LOAD
	IRC0= IELE + NELE + 2*NELE*(ILC-1) + NELE
	WRITE(KWREC,REC=IRC0) B(1:NSTH)   !FOR INITIAL STRESS FOR CONT LOAD
	ENDDO

10    CONTINUE

	ENDDO 


	DEALLOCATE(A,B)


	RETURN
	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE ADWBAK(IEG,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 OPER
	ALLOCATABLE WA(:)
C	=============================================================
C     BACKUP PROCESS OF WORKING ARRAY      
C	=============================================================
	IF(OPER.EQ.'CAL') IND = 1  !CALLING WORKING ARRAY FM BACK UP FILE 3200+IEG
	IF(OPER.EQ.'BAK') IND = 2  !BACK UP WORKING ARRAY TO BACK UP FILE 3200+IEG
	IF(OPER.EQ.'BAD') IND = 3  !BACK UP WORKING ARRAY TO BACK UP FILE 3200+IEG DEEP BACKUP

	NFLB = 3200
	KWBAK = NFLB + IEG

	NFL0 = 3300
	KWREC = NFL0 + IEG

	CALL MINTFIL('GWOK',NELE,IEG,1,0) !CALL NUMBER OF ELEMENT IN THIS GROUP
	CALL MINTFIL('GWOK',MELE,IEG,2,0) !CALL ACCUMULATE NUMBER OF ELEMENT IN THIS GROUP
	CALL MINTFIL('GWOK',NWA ,IEG,3,0) !CALL NUMBER OF WORKING ARRAY OF EACH ELEMENT IN THIS GROUP

	ALLOCATE(WA(NWA))

	DO IRC = 1,NELE

	SELECTCASE(IND)

	CASE(1)
	 READ(KWBAK,REC=IRC) WA(1:NWA)
	WRITE(KWREC,REC=IRC) WA(1:NWA)
	CASE(2)
	 READ(KWREC,REC=IRC) WA(1:NWA)
	WRITE(KWBAK,REC=IRC) WA(1:NWA)
	CASE(3)
	 READ(KWREC,REC=IRC) WA(1:NWA)
	WRITE(KWBAK,REC=IRC) WA(1:NWA)
	IRD = IRC + NELE
	WRITE(KWBAK,REC=IRD) WA(1:NWA)

	ENDSELECT

	ENDDO


	DEALLOCATE(WA)

	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE ADWINIT(NEG)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*3 OPER
	ALLOCATABLE WA(:)
C	=============================================================
C	REFRESH THE WORKING ARRAY BY CALLING FROM DEEP RECORD
C	=============================================================

	DO 100 IEG = 1,NEG

	NFLB = 3200
	KWBAK = NFLB + IEG

	NFL0 = 3300
	KWREC = NFL0 + IEG

	CALL MINTFIL('GWOK',NELE,IEG,1,0)
	CALL MINTFIL('GWOK',MELE,IEG,2,0)
	CALL MINTFIL('GWOK',NWA ,IEG,3,0)

	ALLOCATE(WA(NWA))


	DO IRC = 1,NELE

	WA(1:NWA) = 0.0D0
	WRITE(KWREC,REC=IRC) WA(1:NWA) !STORE DATA TO NORMAL WORKING ARRAY FILE FIRST
	IRD = IRC + NELE
	 READ(KWBAK,REC=IRD) WA(1:NWA) !CALL DATA_READ FROM DEEP BACKUP
	WRITE(KWBAK,REC=IRC) WA(1:NWA) !WRITE IT TO NORMAL BACKUP RECORD (READY TO USE NOW)

	ENDDO


	DEALLOCATE(WA)

100	CONTINUE


	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE ADREWT(IEG,IELE,WA,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	CHARACTER*3 OPER
	DIMENSION WA(1)
C	=============================================================
	IF(OPER.EQ.'RED') IND = 1
	IF(OPER.EQ.'WRT') IND = 2

	NFL0  = 3300
	KWREC = NFL0 + IEG

	CALL MINTFIL('GWOK',NELE,IEG,1,0)
	CALL MINTFIL('GWOK',MELE,IEG,2,0)
	CALL MINTFIL('GWOK',NWA ,IEG,3,0)

C	MELE0 = MELE - NELE
	IRC   = IELE 

	SELECTCASE(IND)

	CASE(1)
	 READ(KWREC,REC=IRC) WA(1:NWA) !CALL WORKING ARRAY
	CASE(2)
	WRITE(KWREC,REC=IRC) WA(1:NWA) !WRITE WORKING ARRAY TO NORMAL STORAGE FILE


	ENDSELECT


	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================

C	=============================================================
C	NEW WORKING ARRAY SAVING BY BJ
C	=============================================================
	SUBROUTINE ADREWTN(IEG,IELE,WA,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM

	CHARACTER*3 OPER
	DIMENSION WA(1)
C	=============================================================
	IF(OPER.EQ.'RED') IND = 1
	IF(OPER.EQ.'WRT') IND = 2

	NFLN  = 3600
	KWRECN = NFLN + IEG

	CALL MINTFIL('GWOK',NELE,IEG,1,0)
	CALL MINTFIL('GWOK',MELE,IEG,2,0)
	CALL MINTFIL('GWOK',NWA ,IEG,3,0)
	
C	MELE0 = MELE - NELE
	IRC   = IELE 

	SELECTCASE(IND)

	CASE(1)
	 READ(KWRECN,REC=IRC) WA(1:NWA)
	CASE(2)
	WRITE(KWRECN,REC=IRC) WA(1:NWA)


	ENDSELECT


	RETURN
	END	

C	=============================================================
C	=============================================================
C	=============================================================

	SUBROUTINE ADRINI(IEG,IELE,WA,ILC,OPER,TYP)   !READ AND WRITE INITIAL STRESS
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	CHARACTER*3 OPER
	CHARACTER*4 TYP
	DIMENSION WA(1)
	ALLOCATABLE A(:)
C	=============================================================
	IF(ILC.LE.0) RETURN

	IF(OPER.EQ.'RED') IND = 1
	IF(OPER.EQ.'WRT') IND = 2
	IF(OPER.EQ.'ADD') IND = 3
	IF(IND.EQ.0) RETURN

	NFL0  = 3300
	KWREC = NFL0 + IEG

	CALL MINTFIL('GWOK',NELE,IEG,1,0)
	CALL MINTFIL('GWOK',NSTH,IEG,4,0)

	ALLOCATE(A(NSTH))

C	REMIND THAT FIRST NELE BLOCK OF 3300 FILE IS FOR STORING WORKING ARRAY

	IRC = 0
	IRO = IELE + NELE*(ILC-1)
	IF(TYP.EQ.'VARY') IRC = (2*IRO-1) + NELE         !LINE NUMBER 2*N-1
	IF(TYP.EQ.'CONT') IRC = (2*IRO-0) + NELE	     !LINE NUMBER 2*N
	IF(IRC.EQ.0) RETURN

	

	SELECTCASE(IND)

	CASE(1) !Read (set to zero if not exist before)
	WA(1:NSTH) = 0.0D0
	 READ(KWREC,REC=IRC,ERR=10) WA(1:NSTH)
10    CONTINUE	
 
	CASE(2) !Write to file
	WRITE(KWREC,REC=IRC) WA(1:NSTH)
	
	CASE(3) !Add to previous (add to zero if not exist before)
	A(1:NSTH) = 0.0D0
	 READ(KWREC,REC=IRC,ERR=20)  A(1:NSTH)
20    CONTINUE	
 	WRITE(KWREC,REC=IRC)  A(1:NSTH) + WA(1:NSTH)

	ENDSELECT


	DEALLOCATE(A)

	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE COLSOL (MAXA,A,D,V,IND,INDPD,TYP1,TYP2)
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

	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
	1              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
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
	ALLOCATABLE NCLM(:),ICOP(:),IHIG(:),B(:)

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


	CALL CPU_TIME (TIM1)

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

C	WRITE(*,*) KCOL,D(KCOL)

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
      SUBROUTINE BKCHOL (MAXA,A,D,V,NEQ,IND,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP
C	BACK SUBSTITUTION OF CHOLESKY TRIANGULAR MATRIX
C
	DIMENSION MAXA(1),A(1),D(1),V(1)
	ALLOCATABLE NCLM(:),ICOP(:),IHIG(:)


	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

	ALLOCATE(NCLM(NBLOCK),ICOP(NBLOCK),IHIG(NBLOCK))

	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCLM(IBLO),1,IBLO,0)
	CALL MINTFIL('BLOC',ICOP(IBLO),2,IBLO,0)
	CALL MINTFIL('BLOC',IHIG(IBLO),3,IBLO,0)
	ENDDO
	

      GOTO(800,900) IND


 800	NCOLB = 0
	CALL MCALFIL(KSREC,TYP)
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
C     IND = 1 FORWARD REDUCTION OF LOAD VECTOR    L*SQRT(D)*V = V
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

	DO IEQ=1,NEQ
	V(IEQ) = V(IEQ)/SQRT(D(IEQ))
	ENDDO


	DEALLOCATE(NCLM,ICOP,IHIG)

	RETURN
C     ------------------------------------------------
C     ------------------------------------------------



 900	NCOLB = NEQ
	CALL MCALFIL(KSREC,TYP)
	REWIND(KSREC)

	DO IEQ=1,NEQ
	V(IEQ) = V(IEQ)/SQRT(D(IEQ))
	ENDDO

C     ----------------------------
C     IND = 2 BACKWARD REDUCTION OF LOAD VECTOR    SQRT(D)*LT*V = V
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

      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MESTIF(LM,S,NWK,NEF,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C     MANIPULATION OF ELEMENT STIFF MATRIX      
C	=============================================================
      CHARACTER*3 OPER
C	=============================================================
	DIMENSION LM(1),S(1)


	CALL MCALFIL(NFL,'ETIF')

	IND = -1
	IF(OPER.EQ.'RED') IND = 0 
	IF(OPER.EQ.'WRT') IND = 1
	IF(OPER.EQ.'REL') IND = 2 !READ ONLY LM 
	IF(IND.EQ.-1) RETURN
	

	NK = NWK
	SELECT CASE(IND)

C	---------------------------------------
	CASE(0)
	 READ(NFL) NEF,NWK,LM(1:NEF),S(1:NWK)
C	---------------------------------------
	CASE(1)
	IF(NK.EQ.0) NK = (NEF*NEF+NEF)/2 

	WRITE(NFL) NEF,NK,LM(1:NEF),S(1:NK)

	CALL MINTFIL('BLOK',1,1,3,2)  !ACCUMULATE NUMBER OF LINE
C	---------------------------------------
	CASE(2)
	 READ(NFL) NEF,NWK,LM(1:NEF)

C	---------------------------------------

	ENDSELECT	 


	RETURN

      END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE ASSEMB (MAXA,NEQ,TYP,KSC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP


      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
      COMMON /OPT_ADDED_STIFFNESS_MATRIX/ NMSTIFF

C     ----------------------------------------------------------------
C     ASSEMBLES GLOBAL COMPACTED STIFFNESS, INITIAL STRESS OR MASS
C     BLOCKS AND WRITES THESE BLOCKS ONTO RANDOM ACCESS TAPE NKFAC
C
C     MAXA(NEQ1)   = ADRESSES OF DIAGONAL ELEM.IN COMPACTED STIFFNESS
C     AA(ISTOR)    = STIFFNESS BLOCK TO BE ASSEMBLED
C	TYP = MATRIX TYPE 'STIF' 'MASS' 'DAMP'
C     ----------------------------------------------------------------
C
      DIMENSION MAXA(1)
      DIMENSION DD(MAXA(NEQ+1) - 1)
	ALLOCATABLE AA(:),LM(:),S(:) ,MAXC(:),BB(:)

	ALLOCATE(LM(200),S(21000))

	CALL CPU_TIME (TIME1)

C	-------------------------------------------
C	IMOPT = 0 SKYLINE + COLSOL , IMOPT = 1 SPARSE  + ITERATION
	CALL MINTFIL('SOLV',IMOPT,1,1,0)		!STIFFNESS PROFILE AND SOLVER OPTION
C	-------------------------------------------

C	-----------------------------------------------------------------------------
C	FOR SPARSE SOLVER - THE STIFFNASS SIZE AND INDEX WILL BE UPDATED BY NEXT SUB.
C	-----------------------------------------------------------------------------
	CALL MCSPARS (MAXA,NEQ,KSC)
      
	CALL CPU_TIME (TIME2)
C	-----------------------------------------------------------------------------

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	CALL MINTFIL('BLOK',MLINE ,1,3, 0)  !NUMBER OF LINE
      
C     -----------------------------------------------

	IF(IMOPT.EQ.1) CALL MCALFIL(NFLCH,'MAXC')	!SPARSE
	IF(IMOPT.EQ.1) REWIND(NFLCH)				!SPARSE

C     -----------------------------------------------
      NWKM = MAXA(NEQ+1) - 1
	ALLOCATE(AA(NWKM),MAXC(NWKM))
      IF(LSYMM.EQ.1) ALLOCATE(BB(NWKM))
      
      AA = 0.0D0
      IF(LSYMM.EQ.1) BB = 0.0D0
      
C     --------------------------------------------- 
	ICC = 1
	DO IBLO = 1,NBLOCK
	    CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	    READ(NFLCH) MAXC(ICC:ICC+NHIG-1)
	    ICC = ICC + NHIG
      ENDDO
C     ---------------------------------------------  

	CALL MCALFIL(KEREC,'ETIF')
	REWIND(KEREC)
      
      NPRE = 0
      NEQF = 1
      NEQL = NEQ
      IF(NEQ.EQ.0) GOTO 250
      
	DO 200 MELM = 1,MLINE

      CALL MESTIF(LM,S,NWK,NEF,'RED')
      
	IF(IMOPT.EQ.0) CALL ASSEME (LM,MAXA,AA,S,NEQF,NEQL,NPRE,NEF		)	!SKYLINE
	IF(IMOPT.EQ.1) CALL ASSEMS (LM,MAXA,AA,S,NEQF,NEQL,NPRE,NEF,MAXC)	!SPARSE
      
	IF(LSYMM.EQ.1) THEN
	NWK2= NWK/2 + 1
	IF(IMOPT.EQ.0) 
	1	CALL ASSEME (LM,MAXA,BB,S(NWK2),NEQF,NEQL,NPRE,NEF	   )	!SKYLINE
	IF(IMOPT.EQ.1) 
	1	CALL ASSEMS (LM,MAXA,BB,S(NWK2),NEQF,NEQL,NPRE,NEF,MAXC)	!SPARSE
      ENDIF

 200  CONTINUE
      
      IF (NMSTIFF.NE.0.AND.TYP.EQ."STIF")THEN
      DD(1:NWKM) = AA(1:NWKM)
      CALL ADDSTIFFNESS (DD,NWKM,MAXA)
      AA(1:NWKM) = DD(1:NWKM)
      ENDIF
C     ---------------------------------------------      

 250  CALL MCALFIL(KSREC,TYP)
	REWIND(KSREC)
      

	IAA = 1
	DO IBLO = 1,NBLOCK
	    CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	    WRITE(KSREC) AA(IAA:IAA+NHIG-1)
	    IF(LSYMM.EQ.1) WRITE(KSREC) BB(IAA:IAA+NHIG-1)
	    IAA = IAA + NHIG
      ENDDO      
C     ---------------------------------------------  
      
      
	CALL CPU_TIME (TIME3)


	DEALLOCATE(LM,S)
      DEALLOCATE(AA,MAXC)
      IF(LSYMM.EQ.1) DEALLOCATE(BB)

C      WRITE(*,3000) TIME2-TIME1
C      WRITE(*,3100) TIME3-TIME2
      
 3000 FORMAT (X,'TIME FOR SPARSE REDUCTION . . .',E10.3)
 3100 FORMAT (X,'TIME FOR STIFNESS ASSEMBLE. . .',E10.3)


      RETURN
      END

C	=============================================================
C	=============================================================
C	=============================================================     
      SUBROUTINE ASSEMB_WORK_WELL_2 (MAXA,NEQ,TYP,KSC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP


      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

C     ----------------------------------------------------------------
C     ASSEMBLES GLOBAL COMPACTED STIFFNESS, INITIAL STRESS OR MASS
C     BLOCKS AND WRITES THESE BLOCKS ONTO RANDOM ACCESS TAPE NKFAC
C
C     MAXA(NEQ1)   = ADRESSES OF DIAGONAL ELEM.IN COMPACTED STIFFNESS
C     AA(ISTOR)    = STIFFNESS BLOCK TO BE ASSEMBLED
C	TYP = MATRIX TYPE 'STIF' 'MASS' 'DAMP'
C     ----------------------------------------------------------------
C
      DIMENSION MAXA(1)
	ALLOCATABLE AA(:),LM(:),S(:) ,MAXC(:),BB(:)
	ALLOCATABLE NBEQ(:,:),MEBL(:),MEBM(:,:)
      ALLOCATABLE LMEMB(:,:)

	ALLOCATE(LM(200),S(21000))

	CALL CPU_TIME (TIME1)

C	-------------------------------------------
C	IMOPT = 0 SKYLINE + COLSOL , IMOPT = 1 SPARSE  + ITERATION
	CALL MINTFIL('SOLV',IMOPT,1,1,0)		!STIFFNESS PROFILE AND SOLVER OPTION
C	-------------------------------------------

C	-----------------------------------------------------------------------------
C	FOR SPARSE SOLVER - THE STIFFNASS SIZE AND INDEX WILL BE UPDATED BY NEXT SUB.
C	-----------------------------------------------------------------------------
	CALL MCSPARS (MAXA,NEQ,KSC)
C	-----------------------------------------------------------------------------

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	CALL MINTFIL('BLOK',MLINE ,1,3, 0)  !NUMBER OF LINE

      
	ALLOCATE(NBEQ(NBLOCK,2),MEBL(NBLOCK))

      ALLOCATE(LMEMB(MLINE,24+2))
      
	CALL MCALFIL(KSREC,TYP)
	REWIND(KSREC)

C     FLAG TO CHECK WHETHER ELEMENT HAS BEEN CALL AT LEAST ONCE OR NOT      
      LMEMB = 0
      
C     -----------------------------------------------
      NEQF = 1
      NEQL = 0
      DO IBLO=1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
      NEQL = NEQL+NCOL
	NBEQ(IBLO,1:2) = [NEQF,NEQL]
	NEQF = NEQF + NCOL
	ENDDO


	CALL MCALFIL(KEREC,'ETIF')
	REWIND(KEREC)
C     -----------------------------------------------

 
C     -----------------------------------------------
C     FIND MAX SIZE OF MEBM ...MAXBM
C     LOOP OVER NUMBER OF LINE RECORDS
C     SCAN ELEMENT EQUATION NUMBERS FOR CONTRIBUTIONS
C     -----------------------------------------------
      MAXBM = 0

      MEBL(1:NBLOCK) = 0
	DO 5 MELM = 1,MLINE
          
C	CALL MESTIF(LM,S,NWK,NEF,'REL')
      ICOUNT = LMEMB(MELM,1)
      IF(ICOUNT.EQ.0) THEN  
	    CALL MESTIF(LM,S,NWK,NEF,'REL')
          IF(NEF.LE.24) THEN
              LMEMB(MELM,1) = 1
              LMEMB(MELM,2) = NEF
          
              LMEMB(MELM,2+1:2+NEF) = LM(1:NEF)
          ENDIF
      ELSE
          NEF = LMEMB(MELM,2)
          LM(1:NEF) = LMEMB(MELM,2+1:2+NEF)
      ENDIF

	DO 4 IBLO = 1,NBLOCK
	NEQF = NBEQ(IBLO,1)
	NEQL = NBEQ(IBLO,2)
	DO IEF=1,NEF
      IEQ = LM(IEF)
      IF(IEQ.GE.NEQF .AND. IEQ.LE.NEQL)  THEN
	MEBL(IBLO) = MEBL(IBLO) + 1
	IF(MEBL(IBLO).GT.MAXBM) MAXBM = MEBL(IBLO)
	GOTO 4
	ENDIF
	ENDDO
 4	CONTINUE

C     -----------------------------------------------
 5    CONTINUE
C     -----------------------------------------------           

	ALLOCATE(MEBM(NBLOCK,MAXBM))
      
	REWIND(KEREC)
C     -----------------------------------------------
C     LOOP OVER NUMBER OF LINE RECORDS
C     SCAN ELEMENT EQUATION NUMBERS FOR CONTRIBUTIONS
C     -----------------------------------------------
      MEBL(1:NBLOCK) = 0
      
	DO 10 MELM = 1,MLINE

C	CALL MESTIF(LM,S,NWK,NEF,'REL')
      ICOUNT = LMEMB(MELM,1)
      IF(ICOUNT.EQ.0) THEN  
	    CALL MESTIF(LM,S,NWK,NEF,'REL')
      ELSE
          NEF = LMEMB(MELM,2)
          LM(1:NEF) = LMEMB(MELM,2+1:2+NEF)
      ENDIF

	DO 8 IBLO = 1,NBLOCK
	NEQF = NBEQ(IBLO,1)
	NEQL = NBEQ(IBLO,2)
	DO IEF=1,NEF
      IEQ = LM(IEF)
      IF(IEQ.GE.NEQF .AND. IEQ.LE.NEQL)  THEN
	MEBL(IBLO) = MEBL(IBLO) + 1
	MEBM(IBLO,MEBL(IBLO)) = MELM
	GOTO 8
	ENDIF
	ENDDO
 8	CONTINUE

C     -----------------------------------------------
 10	CONTINUE
C     -----------------------------------------------

	IF(IMOPT.EQ.1) CALL MCALFIL(NFLCH,'MAXC')	!SPARSE
	IF(IMOPT.EQ.1) REWIND(NFLCH)				!SPARSE
C     --------------
C     INITIALISATION
C     --------------
      NEQF = 1
      NEQL = 0
C     --------------------------
C     LOOP OVER NUMBER OF BLOCKS
C     --------------------------
      DO 900  IBLO=1,NBLOCK

	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	ALLOCATE(AA(NHIG))
	IF(LSYMM.EQ.1) ALLOCATE(BB(NHIG))

	IF(IMOPT.EQ.1) ALLOCATE(MAXC(NHIG))	!SPARSE

      CALL CLEARA (AA,NHIG)
      IF(LSYMM.EQ.1) CALL CLEARA (BB,NHIG)

	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
      NEQL   = NEQL+NCOL
      NPRE   = MAXA(NEQF)-1

	CALL MCALFIL(KEREC,'ETIF')
	REWIND(KEREC)

	IF(IMOPT.EQ.1) READ(NFLCH) MAXC(1:NHIG)		!SPARSE
C     ---------------------------------------------
      
      NUME = 0
	DO 600 MELM = 1,MLINE

      ICOUNT = 0
      DO IEBL = 1,MEBL(IBLO)
          KELM = MEBM(IBLO,IEBL)
          IF(KELM.EQ.MELM) THEN
              ICOUNT = 1
              EXIT
          ENDIF
      ENDDO
      
      IF(ICOUNT.EQ.1) THEN  
	    CALL MESTIF(LM,S,NWK,NEF,'RED')
          NUME = NUME + 1
      ELSE	
          READ(KEREC)
          GOTO 600
      ENDIF
      
      
	IF(IMOPT.EQ.0) CALL ASSEME (LM,MAXA,AA,S,NEQF,NEQL,NPRE,NEF		)	!SKYLINE
	IF(IMOPT.EQ.1) CALL ASSEMS (LM,MAXA,AA,S,NEQF,NEQL,NPRE,NEF,MAXC)	!SPARSE
      
	IF(LSYMM.EQ.1) THEN
	NWK2= NWK/2 + 1
	IF(IMOPT.EQ.0) 
	1	CALL ASSEME (LM,MAXA,BB,S(NWK2),NEQF,NEQL,NPRE,NEF	   )	!SKYLINE
	IF(IMOPT.EQ.1) 
	1	CALL ASSEMS (LM,MAXA,BB,S(NWK2),NEQF,NEQL,NPRE,NEF,MAXC)	!SPARSE
      ENDIF

C     IF WE ALREADY READ ALL ELEMENT IN CURRENT BLOCK, THEN WE JUMP TO NEXT BLOCK       
      IF(NUME.EQ.MEBL(IBLO)) EXIT

 600	CONTINUE
C     ---------------------------------------------

C     ------------------------------------------------------
C     END OF BLOCK LOOP, WRITE ONTO RANDOM ACCESS TAPE 
C     ------------------------------------------------------
      IF (IBLO.EQ.1)THEN
C          AA(1) = AA(1)-10D0+10000D0
      ENDIF
	WRITE(KSREC) AA(1:NHIG)
C	TOEY WORK
	!WRITE(300,1) AA(1:NHIG)
1     FORMAT (E12.5)	
	IF(LSYMM.EQ.1) WRITE(KSREC) BB(1:NHIG)
	DEALLOCATE(AA)
	IF(LSYMM.EQ.1) DEALLOCATE(BB)
	IF(IMOPT.EQ.1) DEALLOCATE(MAXC)	!SPARSE
C     ------------------------------------------------------
 900  NEQF = NEQF + NCOL
C     ------------------------------------------------------

	CALL CPU_TIME (TIME2)

C	WRITE(*,3000) TIME2-TIME1


	DEALLOCATE(LM,S)
	DEALLOCATE(NBEQ,MEBL,MEBM)
      DEALLOCATE(LMEMB)


 3000 FORMAT (X,'TIME TO PERFORM STIFNESS ASSEMBLE. . .',E15.6)


      RETURN
      END

C	=============================================================
C	=============================================================
C	=============================================================      
      SUBROUTINE ASSEMB_WORK_WELL (MAXA,NEQ,TYP,KSC)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP


      COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC

C     ----------------------------------------------------------------
C     ASSEMBLES GLOBAL COMPACTED STIFFNESS, INITIAL STRESS OR MASS
C     BLOCKS AND WRITES THESE BLOCKS ONTO RANDOM ACCESS TAPE NKFAC
C
C     MAXA(NEQ1)   = ADRESSES OF DIAGONAL ELEM.IN COMPACTED STIFFNESS
C     AA(ISTOR)    = STIFFNESS BLOCK TO BE ASSEMBLED
C	TYP = MATRIX TYPE 'STIF' 'MASS' 'DAMP'
C     ----------------------------------------------------------------
C
      DIMENSION MAXA(1)
	ALLOCATABLE AA(:),LM(:),S(:) ,MAXC(:),BB(:)
	ALLOCATABLE NBEQ(:,:),MEBL(:),MEBM(:,:)

	ALLOCATE(LM(200),S(21000))

	CALL CPU_TIME (TIME1)

C	-------------------------------------------
C	IMOPT = 0 SKYLINE + COLSOL , IMOPT = 1 SPARSE  + ITERATION
	CALL MINTFIL('SOLV',IMOPT,1,1,0)		!STIFFNESS PROFILE AND SOLVER OPTION
C	-------------------------------------------

C	-----------------------------------------------------------------------------
C	FOR SPARSE SOLVER - THE STIFFNASS SIZE AND INDEX WILL BE UPDATED BY NEXT SUB.
C	-----------------------------------------------------------------------------
	CALL MCSPARS (MAXA,NEQ,KSC)
C	-----------------------------------------------------------------------------

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	CALL MINTFIL('BLOK',MLINE ,1,3, 0)  !NUMBER OF LINE

      
	ALLOCATE(NBEQ(NBLOCK,2),MEBL(NBLOCK))

	CALL MCALFIL(KSREC,TYP)
	REWIND(KSREC)

C     -----------------------------------------------
      NEQF = 1
      NEQL = 0
      DO IBLO=1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
      NEQL = NEQL+NCOL
	NBEQ(IBLO,1:2) = [NEQF,NEQL]
	NEQF = NEQF + NCOL
	ENDDO


	CALL MCALFIL(KEREC,'ETIF')
	REWIND(KEREC)
C     -----------------------------------------------

 
C     -----------------------------------------------
C     FIND MAX SIZE OF MEBM ...MAXBM
C     LOOP OVER NUMBER OF LINE RECORDS
C     SCAN ELEMENT EQUATION NUMBERS FOR CONTRIBUTIONS
C     -----------------------------------------------
      MAXBM = 0

      MEBL(1:NBLOCK) = 0
	DO 5 MELM = 1,MLINE
	CALL MESTIF(LM,S,NWK,NEF,'REL')

	DO 4 IBLO = 1,NBLOCK
	NEQF = NBEQ(IBLO,1)
	NEQL = NBEQ(IBLO,2)
	DO IEF=1,NEF
      IEQ = LM(IEF)
      IF(IEQ.GE.NEQF .AND. IEQ.LE.NEQL)  THEN
	MEBL(IBLO) = MEBL(IBLO) + 1
	IF(MEBL(IBLO).GT.MAXBM) MAXBM = MEBL(IBLO)
	GOTO 4
	ENDIF
	ENDDO
 4	CONTINUE

C     -----------------------------------------------
 5    CONTINUE
C     -----------------------------------------------           

	ALLOCATE(MEBM(NBLOCK,MAXBM))
      
	REWIND(KEREC)
C     -----------------------------------------------
C     LOOP OVER NUMBER OF LINE RECORDS
C     SCAN ELEMENT EQUATION NUMBERS FOR CONTRIBUTIONS
C     -----------------------------------------------
      MEBL(1:NBLOCK) = 0
      
	DO 10 MELM = 1,MLINE
	CALL MESTIF(LM,S,NWK,NEF,'REL')

	DO 8 IBLO = 1,NBLOCK
	NEQF = NBEQ(IBLO,1)
	NEQL = NBEQ(IBLO,2)
	DO IEF=1,NEF
      IEQ = LM(IEF)
      IF(IEQ.GE.NEQF .AND. IEQ.LE.NEQL)  THEN
	MEBL(IBLO) = MEBL(IBLO) + 1
	MEBM(IBLO,MEBL(IBLO)) = MELM
	GOTO 8
	ENDIF
	ENDDO
 8	CONTINUE

C     -----------------------------------------------
 10	CONTINUE
C     -----------------------------------------------

	IF(IMOPT.EQ.1) CALL MCALFIL(NFLCH,'MAXC')	!SPARSE
	IF(IMOPT.EQ.1) REWIND(NFLCH)				!SPARSE
C     --------------
C     INITIALISATION
C     --------------
      NEQF = 1
      NEQL = 0
C     --------------------------
C     LOOP OVER NUMBER OF BLOCKS
C     --------------------------
      DO 900  IBLO=1,NBLOCK

	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	ALLOCATE(AA(NHIG))
	IF(LSYMM.EQ.1) ALLOCATE(BB(NHIG))

	IF(IMOPT.EQ.1) ALLOCATE(MAXC(NHIG))	!SPARSE

      CALL CLEARA (AA,NHIG)
      IF(LSYMM.EQ.1) CALL CLEARA (BB,NHIG)

	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
      NEQL   = NEQL+NCOL
      NPRE   = MAXA(NEQF)-1

	CALL MCALFIL(KEREC,'ETIF')
	REWIND(KEREC)

	IF(IMOPT.EQ.1) READ(NFLCH) MAXC(1:NHIG)		!SPARSE
C     ---------------------------------------------
	IEL0 = 1
	DO 600 IEBL = 1,MEBL(IBLO)
	IELM = MEBM(IBLO,IEBL)
	DO II = IEL0,IELM-1
	READ(KEREC)
	ENDDO
	IEL0 = IELM+1
	CALL MESTIF(LM,S,NWK,NEF,'RED')
      
	IF(IMOPT.EQ.0) CALL ASSEME (LM,MAXA,AA,S,NEQF,NEQL,NPRE,NEF		)	!SKYLINE
	IF(IMOPT.EQ.1) CALL ASSEMS (LM,MAXA,AA,S,NEQF,NEQL,NPRE,NEF,MAXC)	!SPARSE
      
	IF(LSYMM.EQ.1) THEN
	NWK2= NWK/2 + 1
	IF(IMOPT.EQ.0) 
	1	CALL ASSEME (LM,MAXA,BB,S(NWK2),NEQF,NEQL,NPRE,NEF	   )	!SKYLINE
	IF(IMOPT.EQ.1) 
	1	CALL ASSEMS (LM,MAXA,BB,S(NWK2),NEQF,NEQL,NPRE,NEF,MAXC)	!SPARSE
	ENDIF


 600	CONTINUE
C     ---------------------------------------------

C     ------------------------------------------------------
C     END OF BLOCK LOOP, WRITE ONTO RANDOM ACCESS TAPE 
C     ------------------------------------------------------
      IF (IBLO.EQ.1)THEN
C          AA(1) = AA(1)-10D0+10000D0
      ENDIF
	WRITE(KSREC) AA(1:NHIG)
C	TOEY WORK
	!WRITE(300,1) AA(1:NHIG)
1     FORMAT (E12.5)	
	IF(LSYMM.EQ.1) WRITE(KSREC) BB(1:NHIG)
	DEALLOCATE(AA)
	IF(LSYMM.EQ.1) DEALLOCATE(BB)
	IF(IMOPT.EQ.1) DEALLOCATE(MAXC)	!SPARSE
C     ------------------------------------------------------
 900  NEQF = NEQF + NCOL
C     ------------------------------------------------------

	CALL CPU_TIME (TIME2)

C	WRITE(*,3000) TIME2-TIME1


	DEALLOCATE(LM,S)
	DEALLOCATE(NBEQ,MEBL,MEBM)


 3000 FORMAT (X,'TIME TO PERFORM STIFNESS ASSEMBLE. . .',E15.6)


      RETURN
      END

C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE BLOCKC(MAXA,NEQ,NCO,NCN,MXBLOK)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
	DIMENSION MAXA(1)
C	=============================================================
	
	NC1   = NCO

	DO IEQ = 1,NEQ

	NC1   = NC1 + 1
	NHIGH = MAXA(NC1) - MAXA(NCO)
	
	IF(NHIGH.GT.MXBLOK) GOTO 100

	NCN   = NC1	

	IF(NC1.EQ.NEQ+1) GOTO 100

	ENDDO

100	CONTINUE


	RETURN
	END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE ASSEME (LM,MAXA,A,S,NEQF,NEQL,NPRE,NEF)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C     ----------------------------------------------------------------
C     ASSEMBLES UPPER TRIANGULAR ELEMENT STIFFNESS INTO COMPACTED
C     GLOBAL STIFFNESS BLOCK (CONTRIBUTIONS BETWEEN NEQF AND NEQL)
C
C     LM(NEF)       = EQUATION NUMBERS FOR ELEMENT DEGREE OF FREEDOMS
C     MAXA(NEQ1)    = ADDRESSES OF DIAGONAL ELEMENTS IN A
C     A(ISTOR)      = GLOBAL COMPACTED STIFFNESS BLOCK
C     S(NWS)        = ELEMENT STIFFNESS MATRIX (UPPER TRIANG. ROW-WISE)
C     NEQF,NEQL     = FIRST AND LAST EQUATION CONTAINED IN BLOCK
C     NPRE          = NUMBER OF PREVIOUS ELEMENTS IN A
C     NEF           = NUMBER OF DEGREES OF FREEDOM FOR ELEMENT
C     ----------------------------------------------------------------
C
      DIMENSION LM(1),MAXA(1),A(1),S(1)
C
      NDI = 0
      DO 200  IEF=1,NEF
      II = LM(IEF)
      IF (II.LT.NEQF .OR. II.GT.NEQL)  GOTO 200
      MI = MAXA(II) - NPRE
      IEF1 = IEF-1
      NIE = 1 + IEF1*NEF - (IEF1*IEF1-IEF1)/2
      DO 220  JEF=1,NEF
      JJ = LM(JEF)
      IF (JJ)  220,220,110
 110  IJ = JJ-II
      IF (IJ)  220,210,210
 210  KK = MI+IJ
      KSS = NIE+JEF-IEF
      IF (JEF.LT.IEF) THEN
      JEF1 = JEF-1
      NJE = 1 + JEF1*NEF - (JEF1*JEF1-JEF1)/2
      KSS = NJE+IEF-JEF
      ENDIF
      A(KK) = A(KK)+S(KSS)
 220  CONTINUE
 200  CONTINUE
C
C
      RETURN

      END
C	=============================================================
C	=============================================================
C	=============================================================

	SUBROUTINE MDSALL(TYP,VALV)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*4 TYP
C	=============================================================
	ALLOCATABLE AA(:)

	CALL MCALFIL(KSREC,TYP)
	REWIND(KSREC)
	CALL MCALFIL(KTEMP,'TEMP')
	REWIND(KTEMP)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

	ALLOCATE(AA(MSTOR))


	DO 1000	IBLO = 1,NBLOCK

	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)

	 READ(KSREC) AA(1:NHIG)
	WRITE(KTEMP) AA(1:NHIG)*VALV

1000	CONTINUE

C	MOVE BACK TO STIF
	CALL MDMOVE('TEMP',TYP)


	DEALLOCATE(AA)


	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MCALFIL(NFL,TYP)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*4 TYP
C	=============================================================

	IF(TYP.EQ.'STIF') NUM = 1
	IF(TYP.EQ.'MASS') NUM = 2
	IF(TYP.EQ.'DAMP') NUM = 3
	IF(TYP.EQ.'STIG') NUM = 4
	IF(TYP.EQ.'STI0') NUM = 5
	IF(TYP.EQ.'MAS0') NUM = 6
	IF(TYP.EQ.'DAM0') NUM = 7
	IF(TYP.EQ.'STG0') NUM = 8
	IF(TYP.EQ.'ETIF') NUM = 9
	IF(TYP.EQ.'TEMP') NUM = 10

	IF(TYP.EQ.'STFR') NUM = 11
	IF(TYP.EQ.'EFTF') NUM = 12
	IF(TYP.EQ.'MAXC') NUM = 13
	IF(TYP.EQ.'ICCG') NUM = 14
	IF(TYP.EQ.'TEM1') NUM = 15
	IF(TYP.EQ.'TEM2') NUM = 16

	CALL MINTFIL('MDSA',NFL,1,NUM,0)


	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MDMOVE(TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*4 TYP1,TYP2
C	=============================================================
	ALLOCATABLE AA(:)
		
	CALL MCALFIL(NFL1,TYP1)
	CALL MCALFIL(NFL2,TYP2)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	ALLOCATE(AA(MSTOR))
	
	REWIND(NFL1)
	REWIND(NFL2)

	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	 READ(NFL1) AA(1:NHIG)
	WRITE(NFL2) AA(1:NHIG)
	ENDDO

	REWIND(NFL1)
	REWIND(NFL2)

	DEALLOCATE(AA)

	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MDMOVI(TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*4 TYP1,TYP2
C	=============================================================
	ALLOCATABLE IA(:)
		
	CALL MCALFIL(NFL1,TYP1)
	CALL MCALFIL(NFL2,TYP2)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	ALLOCATE(IA(MSTOR))
	
	REWIND(NFL1)
	REWIND(NFL2)

	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	 READ(NFL1) IA(1:NHIG)
	WRITE(NFL2) IA(1:NHIG)
	ENDDO

	REWIND(NFL1)
	REWIND(NFL2)

	DEALLOCATE(IA)

	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE SEQOPEN(LUN)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C-----OPEN THE UNFORMATTED DIRECT ACCESS FILE "FIN.FXT"-----
C        LRECL = NUMBER OF BYTES IN A RECORD

      OPEN (LUN,STATUS='SCRATCH',FORM='UNFORMATTED')
C      OPEN (LUN,STATUS='SCRATCH',FORM='UNFORMATTED',BUFFERED='YES')
C
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MDOPER(TYP1,TYP2,TYP3,FAC1,FAC2,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*4 TYP1,TYP2,TYP3
      CHARACTER*3 OPER
C	=============================================================
	ALLOCATABLE AA(:),BB(:)
		

	IND = 0
	IF(OPER.EQ.'ADD') IND = 1
	IF(OPER.EQ.'SUB') IND = 2
	IF(OPER.EQ.'MOV') IND = 3
	IF(OPER.EQ.'DIV') IND = 4
	IF(IND.EQ.0) RETURN

	CALL MCALFIL(NFL1,TYP1)
	CALL MCALFIL(NFL2,TYP2)
	CALL MCALFIL(NFL3,TYP3)

	REWIND(NFL1)
	REWIND(NFL2)
	REWIND(NFL3)

	IF(TYP1.EQ.TYP2) THEN
	CALL MCALFIL(NFLT,'TEM1')
	CALL  MDMOVE(TYP2,'TEM1')
	REWIND(NFLT)

	ELSEIF(TYP1.EQ.TYP3) THEN
	CALL MCALFIL(NFLT,'TEM1')
	CALL  MDMOVE(TYP3,'TEM1')
	REWIND(NFLT)

	ELSEIF(TYP2.EQ.TYP3) THEN
	CALL MCALFIL(NFLT,'TEM1')
	CALL  MDMOVE(TYP2,'TEM1')
	REWIND(NFLT)
	ENDIF


	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	ALLOCATE(AA(MSTOR),BB(MSTOR))
	

	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)

	IF(TYP1.EQ.TYP2) THEN
	 READ(NFLT) AA(1:NHIG)
	 READ(NFL3) BB(1:NHIG)
	ELSEIF(TYP1.EQ.TYP3) THEN
	 READ(NFL2) AA(1:NHIG)
	 READ(NFLT) BB(1:NHIG)
	ELSEIF(TYP2.EQ.TYP3) THEN
	 READ(NFLT) AA(1:NHIG)
	 READ(NFL3) BB(1:NHIG)
	ELSE
	 READ(NFL2) AA(1:NHIG)
	 READ(NFL3) BB(1:NHIG)
	ENDIF
	
		
	SELECTCASE(IND)
	CASE(1)
	AA(1:NHIG) = FAC1*AA(1:NHIG) + FAC2*BB(1:NHIG)
	CASE(2)
	AA(1:NHIG) = FAC1*AA(1:NHIG) - FAC2*BB(1:NHIG)
	CASE(3)
	AA(1:NHIG) = FAC1*AA(1:NHIG) * FAC2*BB(1:NHIG)
	CASE(4)
	AA(1:NHIG) = FAC1*AA(1:NHIG) / FAC2*BB(1:NHIG)
	ENDSELECT

	WRITE(NFL1) AA(1:NHIG)
	ENDDO


	DEALLOCATE(AA,BB)

	RETURN

	END


C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE MDCALL_DIA(TYP1,DIA,OPER,MAXA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*4 TYP1,TYP2
      CHARACTER*3 OPER
C	=============================================================
	DIMENSION MAXA(1),DIA(1)
	ALLOCATABLE AA(:)
		

	IND = 0
	IF(OPER.EQ.'RED') IND = 1
	IF(OPER.EQ.'WRT') IND = 2
	IF(OPER.EQ.'ADD') IND = 3
	IF(OPER.EQ.'SUB') IND = 4
	IF(OPER.EQ.'MUL') IND = 5
	IF(OPER.EQ.'DIV') IND = 6
	IF(IND.EQ.0) RETURN

	CALL MCALFIL(NFL1,TYP1)
	REWIND(NFL1)
	CALL MCALFIL(NFLT,'TEM1')
	CALL  MDMOVE(TYP1,'TEM1')
	REWIND(NFLT)


	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	ALLOCATE(AA(MSTOR))
	

      NEQF = 1
      NEQL = 0
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
      NEQL = NEQL+NCOL
	NPRE = MAXA(NEQF) - 1
	 READ(NFLT) AA(1:NHIG)
	
	DO IEQ = NEQF,NEQL
	IDIA = MAXA(IEQ) - NPRE
	SELECTCASE(IND)
	CASE(1)
	DIA(IEQ) = AA(IDIA)
	CASE(2)
	AA(IDIA) = DIA(IEQ)
	CASE(3)
	AA(IDIA) = AA(IDIA) + DIA(IEQ)
	CASE(4)
	AA(IDIA) = AA(IDIA) - DIA(IEQ)
	CASE(5)
	AA(IDIA) = AA(IDIA) * DIA(IEQ)
	CASE(6)
	AA(IDIA) = AA(IDIA) / DIA(IEQ)
	ENDSELECT
	ENDDO

	IF(IND.NE.1) WRITE(NFL1) AA(1:NHIG)
	NEQF = NEQF + NCOL
	ENDDO


	DEALLOCATE(AA)


	RETURN

	END



C	=============================================================
C	=============================================================
C	=============================================================

	SUBROUTINE MDMAKE_DIA(TYP,FAC,MAXA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
C	MAKE DIAGONAL MATRIX WITH 'FAC' VALUE
      CHARACTER*4 TYP
C	=============================================================
	DIMENSION MAXA(1)
	ALLOCATABLE AA(:)
		
	CALL MCALFIL(NFL,TYP)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	ALLOCATE(AA(MSTOR))
	
	REWIND(NFL)

      NEQF = 1
      NEQL = 0
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
      NEQL = NEQL+NCOL
	NPRE = MAXA(NEQF) - 1
	AA(1:NHIG) = 0.0D0

	DO IEQ = NEQF,NEQL
	IDIA = MAXA(IEQ) - NPRE
	AA(IDIA) = FAC
	ENDDO

	WRITE(NFL) AA(1:NHIG)
	NEQF = NEQF + NCOL
	ENDDO


	DEALLOCATE(AA)


	RETURN

	END



C	=============================================================
C	=============================================================
C	=============================================================	
      SUBROUTINE MAMULT (MAXA,B,RR,TT,TYP,OPER)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 TYP
      CHARACTER*3 OPER

C     ----------------------------------------------------------------
C     EVALUATES MATRIX PRODUCT OF {TT} =        [B]*{RR} IF OPER='STD' 
C     EVALUATES MATRIX PRODUCT OF {TT} = {TT} + [B]*{RR} IF OPER='ADD' 
C	----------------------------------------------------------------
      DIMENSION MAXA(1),B(1),RR(1),TT(1)

C	-------------------------------------------
C	IMOPT = 0 SKYLINE + COLSOL , IMOPT = 1 SPARSE  + ITERATION
	CALL MINTFIL('SOLV',IMOPT,1,1,0)		!STIFFNESS PROFILE AND SOLVER OPTION
	IF(IMOPT.EQ.1) THEN
	CALL MASULT (MAXA,B,RR,TT,TYP,OPER)			!PERFORM SPARSE MULTIPLICATION
	RETURN
	ENDIF
C	-------------------------------------------

	CALL MCALFIL(NFL,TYP)
	REWIND(NFL)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

      NEQF = 1
      NEQL = 0
	DO 1000 IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
      NEQL = NEQL+NCOL
	NPRE = MAXA(NEQF) - 1

	READ(NFL) B(1:NHIG)
C     --------------------------------------
C     LOOP OVER NUMBER OF COLUMNS IN BLOCK
C     MULTIPLY COLUMN COEFFICIENTS WITH (RR)
C     --------------------------------------
      DO 500  IEQ = NEQF,NEQL
      KLOW  = MAXA(IEQ)       - NPRE
      KUPP  = MAXA(IEQ+1) - 1 - NPRE
      KT    = IEQ + 1
      R     = RR(IEQ)
	IF(OPER.EQ.'STD') TT(IEQ) = 0.0D0  !IF(OPER.EQ.'ADD') ---> NO INITIALIZATION 
      DO 300  KB = KLOW,KUPP
      KT    = KT - 1
 300  TT(KT)= TT(KT) + B(KB)*R
C     ------------------------------------------------
C     MULTIPLY SYMMETRIC COEFFICIENTS OF ROW WITH (RR)
C     ------------------------------------------------
      KLOW = KLOW + 1
      IF (KUPP-KLOW) 500,410,410
 410  KR = IEQ
      CC = 0.0
      DO 490  KB = KLOW,KUPP
      KR = KR - 1
 490  CC = CC + B(KB)*RR(KR)
      TT(IEQ) = TT(IEQ) + CC
C
 500  CONTINUE


	NEQF = NEQF + NCOL
1000	CONTINUE
C     --------------------------------------


C
      RETURN
      END
C
C	=============================================================
C	=============================================================
C	=============================================================
	SUBROUTINE AZFILL(ID,MAXA,VEC,ILC,AA,TYP1,TYP2)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
	CHARACTER*4 TYP1,TYP2
C	-----------------------------------------------------------------------
C	-----------------------------------------------------------------------
	COMMON /NUMB/ HED(20),MODEX,NRE,NSN,NEG,NBS,NLS,NLA,
     +              NSC,NSF,IDOF(9),LCS,ISOLOP,LSYMM,ICONTROLSPEC
	COMMON /SOLU/ NEQ,NEQ1,NBLOCK,MK,BM,NWK,NWM,ISTOR,NFAC,
     1              NRED,KPOSD,DETK,DET1,DAVR,STOL
C	-----------------------------------------------------------------------
	COMMON /SETTM/ NSETT,LSTL             !SETTLEMENT MODULE
C	-----------------------------------------------------------------------
	DIMENSION ID(NSF,NSN),MAXA(NEQ+1),VEC(NEQ),AA(1)
C	-----------------------------------------------------------------------

C	-------------------------------------------
	CALL MCALFIL(NFL1,TYP1)
	CALL MCALFIL(NFL2,TYP2)

	CALL MDMOVE(TYP1,TYP2)

	REWIND(NFL1)
	REWIND(NFL2)

	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

      NEQF = 1
      NEQL = 0
	DO 1000 IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
      NEQL = NEQL+NCOL
	NPRE = MAXA(NEQF) - 1

	READ(NFL2) AA(1:NHIG)
C	-----------------------------------------------------------------------
	DO II = 1,NSETT

	CALL MRELFIL('SETL',FNODS,II,1,0) !CALLING THE SETTLEMENT DATA
	CALL MRELFIL('SETL',FIDIR,II,2,0)
	CALL MRELFIL('SETL',VALU ,II,3,0)
	CALL MRELFIL('SETL',FILCN,II,4,0)
	CALL MRELFIL('SETL',FILCC,II,5,0)
	CALL MRELFIL('SETL',FIFG ,II,6,0)
	NODS = INT(FNODS)
	IDIR = INT(FIDIR)
	ILN  = INT(FILCN)
	ILO  = INT(FILCC)
	IFG  = INT(FIFG )

	IDIR = IDOFCALL(IDOF,IDIR)
	IEQ  = ID(IDIR,NODS)

	IF(IEQ.GE.NEQF.AND.IEQ.LE.NEQL) THEN
	IF(ILN.EQ.ILC.OR.ILO.EQ.ILC) THEN
	IF(IFG.EQ.0.AND.IEQ.NE.0) THEN
	L1 = MAXA(IEQ  ) - NPRE
	AA(L1)     = 1.0D0
	VEC(IEQ)   = 0.0D0
	ENDIF
	ENDIF
	ENDIF

	ENDDO
C	-----------------------------------------------------------------------
	NEQF = NEQF + NCOL

	WRITE(NFL1) AA(1:NHIG)
1000	CONTINUE
C     --------------------------------------

	RETURN
	END

C	=============================================================
C	=============================================================
C	=============================================================



	SUBROUTINE MDCALL_WRT(TYP1,MAXA)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C	=============================================================
      CHARACTER*4 TYP1,TYP2
      CHARACTER*3 OPER
C	=============================================================
	DIMENSION MAXA(1),DIA(1)
	ALLOCATABLE AA(:)
		

	IND = 1
	IF(IND.EQ.0) RETURN

	CALL MCALFIL(NFL1,TYP1)
	REWIND(NFL1)
	CALL MCALFIL(NFLT,'TEM1')
	CALL  MDMOVE(TYP1,'TEM1')
	REWIND(NFLT)


	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)
	ALLOCATE(AA(MSTOR))
	

      NEQF = 1
      NEQL = 0
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NCOL,1,IBLO,0)
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
      NEQL = NEQL+NCOL
	NPRE = MAXA(NEQF) - 1
	 READ(NFLT) AA(1:NHIG)
	DO II = 1,NHIG
	WRITE(110,*) AA(II)
	ENDDO


	IF(IND.NE.1) WRITE(NFL1) AA(1:NHIG)
	NEQF = NEQF + NCOL
	ENDDO


	DEALLOCATE(AA)


	RETURN

	END



C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE COLSYM (A,D,V,IND,NEQ,INDPD)
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

	DIMENSION MAXA(NEQ+1),A(1),D(1),V(1)
C
	MAXA(1) = 1
	DO IEQ = 1,NEQ
	MAXA(IEQ+1) = MAXA(IEQ) + IEQ
	ENDDO

	STOL = 1.0E-10

	NCOLUM = NEQ
      IF (IND-2) 400,800,800
C     --------------------------------------------------------
C     FACTORIZE STIFFNESS MATRIX A (A = L*D*L'T DECOMPOSISION)
C     --------------------------------------------------------
C
C     -----------------------------------------
C     REDUCE BLOCK, LOOP OVER NUMBER OF COLUMNS
C     -----------------------------------------
 400  DO 690  ICOL=1,NCOLUM
      KDIA = MAXA(ICOL)
      KLOW = KDIA+1
      KUPP = MAXA(ICOL+1)-1
      KCOL = ICOL
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
      DO 490  JCOL=JCF,JCL
      LCOP = LCOP+1
      KTOP = KTOP-1
      LDIA = MAXA(JCOL)
      LHEI = MAXA(JCOL+1)-1-LDIA
      IF (LHEI) 490,490,430
 430  NL = MIN0(LCOP,LHEI)
      C = 0.0
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
 650  IF (INDPD-1)       630,630,660
 660  PIV = STOL
      IF (PIV.EQ.0.0) PIV = -1.0D-16
      D(KCOL) = PIV
      A(KDIA) = PIV
 690  CONTINUE
C     ------------------------------------------------
C     END OF LOOP
C     CALCULATE DETERMINANT OF STIFFNESS MATRIX (DETK)
C     ------------------------------------------------


C     --------------------------------
C     FORWARD REDUCTION OF LOAD VECTOR
C     --------------------------------
C
C     ------------------------------------------------
C     LOOP OVER NUMBER OF COLUMNS, VI = RI-SUM(IRI*VR)
C     ------------------------------------------------
 800	DO 880  ICOL=1,NCOLUM
      KLOW = MAXA(ICOL)+1
      KUPP = MAXA(ICOL+1)-1
      IF (KUPP-KLOW) 880,830,830
 830  IEQ = ICOL
      KV = IEQ
      C = 0.0
      DO 850  KA=KLOW,KUPP
      KV = KV-1
 850  C = C + A(KA)*V(KV)
      V(IEQ) = V(IEQ) - C
 880  CONTINUE
C     -----------------------------
C     BACK SUBSTITUTION , V=[D]-1*V
C     -----------------------------
 890	DO 900  IEQ=1,NEQ
 900  V(IEQ) = V(IEQ)/D(IEQ)
C     ---------------------------
C     LOOP OVER NUMBER OF COLUMNS
C     ---------------------------
 910	KCOL = NCOLUM
      DO 980  ICOL=1,NCOLUM
      KLOW = MAXA(KCOL)+1
      KUPP = MAXA(KCOL+1)-1
      IF (KUPP-KLOW) 980,930,930
 930  IEQ = KCOL
      KV = IEQ
      DO 950  KA=KLOW,KUPP
      KV = KV-1
 950  V(KV) = V(KV) - A(KA)*V(IEQ)
 980  KCOL = KCOL-1
C


      RETURN
      END
C
C	=====================================================================
C	=====================================================================
C	=====================================================================
