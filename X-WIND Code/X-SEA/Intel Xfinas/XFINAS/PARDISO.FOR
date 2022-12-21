C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE CLRPARD
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION V(1),MAXA(2)
      NEQ = 1
      MAXA(1:2) = [1,2]
      V(1) = 0.0D0
      IND = 3  !CLEAR ALL PARDISO MEMORY
      CALL PARDISI(MAXA,V,NEQ,IND,'TEMP','TEMP')  
      END
C	=============================================================
C	=============================================================
C	=============================================================
      SUBROUTINE PARDISI(MAXA,V,NEQ,IND,TYP1,TYP2)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C      INCLUDE 'MKL_PARDISO.F77'
      CHARACTER*4 TYP1,TYP2  
      
C     FOR PARDISO x64
      INTEGER*8 LPT,LPTO
      
C     =============================================================        
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C     STORE PARDISO POINTER NUMBER 
      COMMON /PARDIS/ LPTO(64) 
C     =============================================================   
C     .. Internal solver memory pointer for 64-bit architectures
C     .. INTEGER*8 LPT(64)
C     .. Internal solver memory pointer for 32-bit architectures
C     .. INTEGER*4 LPT(64)
C     .. This is OK in both cases
      ALLOCATABLE LPT(:)
C     ============================================================= 
      POINTER ::  AA(:),MAXC(:),IPERM(:),X(:) 
C      ALLOCATABLE AA(:),MAXC(:),IPERM(:),X(:)   
      DIMENSION MAXA(NEQ),V(NEQ),IPARM(64),IIPARM(64)
C     =============================================================   
      
      ALLOCATE(IPERM(NEQ),X(NEQ))
      
C     IPFAC DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,0)  

      !NWK = MAXA(NEQ)
      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

	ALLOCATE(AA(NWK),MAXC(NWK))
      
C     FOR CLEARING MEMORY PHASE NO NEED TO READ THE STIFFNESS MATRIX 
	IF(IND.EQ.3) GOTO 50
	
      IF(IND.EQ.1) THEN
	CALL MDMOVE(TYP1,TYP2)
	ENDIF

	CALL MCALFIL(KSREC,TYP2)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
	
	IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	READ(NFLCH) MAXC(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
CC      WRITE (300,"(E12.5)") AA !!
      
C..
C.. Set up PARDISO control parameter
C..
50    NRHS = 1
      MAXFCT = 1
      MNUM = 1
      
      DO I = 1,64
         IPARM(I) = 0
      ENDDO
      IPARM(1)  =  1 ! no solver default
      IPARM(2)  =  2 !3 ! 2 = fill-in reordering from METIS , 3 = parallel (OpenMP) nested dissection algorithm 
      IPARM(3)  =  1 !MKL_GET_MAX_THREADS() !1 ! numbers of processors
      IPARM(4)  =  0 !32 ! no iterative-direct algorithm
      IPARM(5)  =  0 ! no user fill-in reducing permutation
      IPARM(6)  =  0 ! =0 solution on the first n compoments of x
      IPARM(7)  =  0 ! not in use
      IPARM(8)  =  9 ! numbers of iterative refinement steps
      IPARM(9)  =  0 ! not in use
      IPARM(10) =  13 ! perturbe the pivot elements with 1E-13
      IPARM(11) =  1 ! use nonsymmetric permutation and scaling MPS
      IPARM(12) =  0 ! not in use
      IPARM(13) =  0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try IPARM(13) = 1 in case of inappropriate accuracy
      IPARM(14) =  0 ! Output: number of perturbed pivots
      IPARM(15) =  0 ! not in use
      IPARM(16) =  0 ! not in use
      IPARM(17) =  0 ! not in use
      IPARM(18) = -1 ! Output: number of nonzeros in the factor LU
      IPARM(19) = -1 ! Output: Mflops for LU factorization
      IPARM(20) =  0 ! Output: Numbers of CG Iterations
      IPARM(27) =  0 ! Output: Numbers of CG Iterations
      IPARM(60) =  1 ! 0 IN-CORE 1 IN&OUT 2 OUT-CORE
      IERROR    =  0 ! initialize error flag
      MSGLVL    =  0 ! print statistical information
      MTYPE     = -2 ! symmetric, indefinite

      ALLOCATE( LPT(64) )
      
      GOTO (100,200,300),IND
      
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
100   IF(IPFAC.EQ.1) GOTO 300

120   DO I = 1, 64
      LPT(I) = 0
      ENDDO

      CALL CPU_TIME (TIME1)
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
c      IPHASE = 11 ! only reordering 
      IPHASE = 12 ! reordering and factorization in the same time
      CALL PARDISO(LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
      
      IF (IERROR.EQ.-2) THEN
          WRITE(*,*) '**PARDISO RUNOUT OF MEMORY**'
          WRITE(*,*) '**NO RESULT HAS BEEN SOLVED**'
          STOP
      ENDIF
      IF (IERROR.EQ.-9) THEN
          MEMX  = IPARM(15)
          MEMX2 = IPARM(16)+IPARM(63)
          IF(MEMX2.GT.MEMX) MEMX = MEMX2
          MEMX = MEMX/1000          
          WRITE(*,*) '**PARDISO RUNOUT OF MEMORY FOR OOC**'
          WRITE(*,*) '**    NO RESULT HAS BEEN SOLVED   **'
          WRITE(*,*) '** CHECK MKL_PARDISO_MAX_CORE_SIZE**'
          WRITE(*,*) '** SET MKL_PARDISO_MAX_CORE_SIZE TO BE GREATER THAN: ',MEMX,' MB'
          STOP
      ENDIF
      
C.. Factorization.      
c      IPHASE = 22 ! only factorization
c      CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
c     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
      CALL CPU_TIME (TIME2)
      TIM(14) =  TIM(14) + (TIME2-TIME1)
      
      DO I = 1, 64
      LPTO(I) = LPT(I)
      ENDDO

C     UPDATE IPFAC
      IPFAC = 1 !STATE THAT THE MEMORY IS ALREADY ALLOCATE
      CALL MINTFIL('BLOK',IPFAC ,1,5 ,1)  

C     PIVOT FLAG  0=POSITIVE DEF  1=NEGATIVE DEF  ONLY APPEARED IN PARDISO,NEWSOLVER,LOADUP    
C     THIS WILL BE USE IN NONLINEAR ANALYIS   SUBROUTINE--LOADUP  	
      IPIVT = 0
      IF(IPARM(23).GT.0) IPIVT = 1 !OPARM(23) CONTAIN NUMBER OF NEGATIVE EIGENVALUE
	CALL MINTFIL('BLOK',IPIVT ,1,6 ,1)  

      GOTO 900
C     --------------------------------------------------------------
C     --------------------------------------------------------------  
200   CALL CPU_TIME (TIME1)
      DO I = 1, 64
      LPT(I) = LPTO(I)
      ENDDO
C.. Back substitution and iterative refinement
      IPARM(8) = 2 ! max numbers of iterative refinement steps
      IPHASE = 33 ! back substitution and iterative refinement
      CALL PARDISO(LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
      CALL CPU_TIME (TIME2)
      TIM(15) =  TIM(15) + (TIME2-TIME1)
      
      V(1:NEQ) = X(1:NEQ)

      GOTO 900
C     --------------------------------------------------------------
C     --------------------------------------------------------------  
C.. Termination and release of memory
300   IF(IPFAC.EQ.0) GOTO 900
      DO I = 1, 64
      LPT(I) = LPTO(I)
      ENDDO
      IPHASE = -1 ! release internal memory
      CALL PARDISO(LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)

C     UPDATE IPFAC
      IPFAC = 0 !STATE THAT THE MEMORY IS ALREADY RELEASE
      CALL MINTFIL('BLOK',IPFAC ,1,5 ,1)  
      
      IF(IND.EQ.1) GOTO 120 !BACK TO FACTORIZATION AFTER RELEASE MEMORY
      GOTO 900
C     --------------------------------------------------------------
C     --------------------------------------------------------------  
     
900	DEALLOCATE(AA,MAXC,IPERM,X,LPT)
	
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================
      
      SUBROUTINE PARDISI_PCG(MAXA,V,NEQ,IND,TYP1,TYP2)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)   
C      INCLUDE 'mkl_rci.fi'
      CHARACTER*4 TYP1,TYP2  

!---------------------------------------------------------------------------
! Define arrays for the upper triangle of the coefficient matrix and rhs vector
! Compressed sparse row storage is used for sparse representation
!---------------------------------------------------------------------------
        INTEGER RCI_request

      
      
C     =============================================================        
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C     =============================================================   
      
C     ---------------------------------------------------------------------------
C     Allocate storage for the solver ?par and the initial solution vector
C     ---------------------------------------------------------------------------
        PARAMETER (length=128, method=1, nRhs=1)
        DIMENSION solution(NEQ)
C     ---------------------------------------------------------------------------
C     Some additional variables to use with the RCI (P)CG solver
C     ---------------------------------------------------------------------------
        DIMENSION ipar(length+2*nRhs),itercount(1)
        DIMENSION dpar(length+2*nRhs),TMP(NEQ,3 + 2*nRhs)
 
        EXTERNAL DNRM2            
      
C     ============================================================= 
      POINTER ::  AA(:),MAXC(:),IPERM(:),X(:),EX(:) 
      DIMENSION MAXA(NEQ),V(NEQ),IPARM(64),IIPARM(64)
      DIMENSION IDUM(1)
C     =============================================================   

      
      ALLOCATE(X(NEQ),EX(NEQ))
      
      X(1:NEQ) = 0.0D0
      
      
C     IPFAC DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,0)  

      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

	ALLOCATE(AA(NWK),MAXC(NWK))
      
      nRows = NEQ
      nCols = NEQ

      

C     FOR CLEARING MEMORY PHASE NO NEED TO READ THE STIFFNESS MATRIX 
	IF(IND.EQ.3) GOTO 50
	
      IF(IND.EQ.1) THEN
	CALL MDMOVE(TYP1,TYP2)
	ENDIF

	CALL MCALFIL(KSREC,TYP2)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
	
	IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	READ(NFLCH) MAXC(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
!      WRITE (300,"(E12.5)") AA
C..
C.. DSS Direct Sparse Solver
C..
50    CONTINUE     

      IF (IND.EQ.1) THEN
C     UPDATE IPFAC
          IPFAC = 1 !STATE THAT THE MEMORY IS ALREADY ALLOCATE
          CALL MINTFIL('BLOK',IPFAC ,1,5 ,1)  
      ENDIF

      
      if (IND.NE.2) GOTO 999
      
C     ---------------------------------------------------------------------------
C      Initialize the right hand side through matrix-vector product
C     ---------------------------------------------------------------------------
C         CALL MKL_DCSRSYMV('U', NEQ, AA, MAXA, MAXC, expected_sol, V)
         
100   CONTINUE          
C     ---------------------------------------------------------------------------
C     Initialize the initial guess
C     ---------------------------------------------------------------------------
      solution(1:NEQ) = 0.0D0
          
200   CONTINUE
      
Cc        CALL dcg_init(NEQ, solution,V, RCI_request,ipar,dpar,TMP)
        
        DO I = 1, length + 2 * nRhs
            ipar(i) = 0
            dpar(i) = 0
        END DO

        CALL dcgmrhs_init(NEQ,solution,nRhs,V,method,RCI_request,ipar,dpar,TMP)
        
        
        IF (RCI_request .NE. 0 ) GO TO 999

300   CONTINUE
C---------------------------------------------------------------------------
C Set the desired parameters:
C LOGICAL parameters:
C do residual stopping test
C do not request for the user defined stopping test
C DOUBLE PRECISION parameters
C set the relative tolerance to 1.0D-5 instead of default value 1.0D-6
C---------------------------------------------------------------------------
        ipar(9) = 1
        ipar(10) = 0
        dpar(1) = 1.D-5    
        
        ipar(5)  = 100
        ipar(11) = 1
        ipar(9)  = 1
        ipar(10) = 0

        dpar(1) = 1.D-5        
        
C     ---------------------------------------------------------------------------
C     Check the correctness and consistency of the newly set parameters
C     ---------------------------------------------------------------------------
        CALL dcg_check(NEQ,solution,V,RCI_request,ipar,dpar,TMP)
        IF (RCI_request .NE. 0 ) GO TO 999        
400   CONTINUE
C     ---------------------------------------------------------------------------
C     Compute the solution by RCI (P)CG solver without preconditioning
C     Reverse Communications starts here
C     ---------------------------------------------------------------------------
450   CALL dcgmrhs(NEQ,solution,nRhs,V,RCI_request,ipar,dpar,TMP)
   
cc450     CALL dcg(NEQ,solution,V,RCI_request,ipar,dpar,TMP)
   
C---------------------------------------------------------------------------
C If RCI_request=0, then the solution was found with the required precision
C---------------------------------------------------------------------------
        IF (RCI_request .EQ. 0) THEN
            GO TO 500
C---------------------------------------------------------------------------
C If RCI_request=1, then compute the vector A*TMP(:,1)
C and put the result in vector TMP(:,2)
C---------------------------------------------------------------------------
        ELSE IF (RCI_request .EQ. 1) THEN
            CALL MKL_DCSRSYMV('U', NEQ, AA, MAXA, MAXC, TMP, TMP(1,2))
            GO TO 450
       ELSE IF (RCI_request .EQ. 3) THEN
            DO i = 1, NEQ
                tmp(i, 3) = tmp (i, 3+ipar(3))
            END DO
            GO TO 450    
        ELSE
C---------------------------------------------------------------------------
C If RCI_request=anything else, then dcg subroutine failed
C to compute the solution vector: solution(N)
C---------------------------------------------------------------------------
            GO TO 999
      END IF

500   CONTINUE
      
Cc      CALL dcg_get(NEQ,solution,V,RCI_request,ipar,dpar,TMP,itercount) 
      CALL dcgmrhs_get(NEQ,solution,nRhs,V,RCI_request,ipar,dpar,TMP,itercount)
      
      if (IND.EQ.2) V(1:NEQ) = SOLUTION(1:NEQ)
      
600   CONTINUE
        
999   CONTINUE     

      CALL MKL_FREE_BUFFERS
      
1000  DEALLOCATE(AA,MAXC,X,EX)
        
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================
      
      SUBROUTINE PARDISI_DSS(MAXA,V,NEQ,IND,TYP1,TYP2)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)   
      INCLUDE 'mkl_dss.f77'
      CHARACTER*4 TYP1,TYP2  
      INTEGER*8 handle
        
C     =============================================================        
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C     =============================================================   
C     ============================================================= 
      POINTER ::  AA(:),MAXC(:),IPERM(:),X(:)  
      DIMENSION MAXA(NEQ),V(NEQ),IPARM(64),IIPARM(64)
      DIMENSION IDUM(1)
C     =============================================================   
      
      ALLOCATE(X(NEQ))
      
      X(1:NEQ) = 0.0D0
      
C     IPFAC DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,0)  

      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

	ALLOCATE(AA(NWK),MAXC(NWK))
      
      nRows = NEQ
      nCols = NEQ
      nRhs = 1

      

C     FOR CLEARING MEMORY PHASE NO NEED TO READ THE STIFFNESS MATRIX 
	IF(IND.EQ.3) GOTO 50
	
      IF(IND.EQ.1) THEN
	CALL MDMOVE(TYP1,TYP2)
	ENDIF

	CALL MCALFIL(KSREC,TYP2)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
	
	IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	READ(NFLCH) MAXC(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
!      WRITE (300,"(E12.5)") AA
C..
C.. DSS Direct Sparse Solver
C..
50    CONTINUE     

      CALL MKL_FREE_BUFFERS
C      if (IND.NE.1) GOTO 400
C---------------------------------------------------------------------------
C Initialize the solver
C---------------------------------------------------------------------------
cc          ipt(1) = MKL_DSS_DEFAULTS
cc          ipt(2) =  MKL_DSS_OOC_STRONG
          ierror = dss_create(handle,MKL_DSS_OOC_STRONG)
          IF (ierror .NE. MKL_DSS_SUCCESS ) GO TO 999 
100   CONTINUE          
C---------------------------------------------------------------------------
C Define the non-zero structure of the matrix
C---------------------------------------------------------------------------
        ierror = dss_define_structure( handle, MKL_DSS_SYMMETRIC,
     &  MAXA, nRows, nCols, MAXC, NWK )
        IF (ierror .NE. MKL_DSS_SUCCESS ) GO TO 999          
          
200     CONTINUE

C---------------------------------------------------------------------------
C Reorder the matrix
C---------------------------------------------------------------------------
        ierror = dss_reorder( handle, MKL_DSS_DEFAULTS, idum)
        IF (ierror .NE. MKL_DSS_SUCCESS ) GO TO 999
300     CONTINUE
C---------------------------------------------------------------------------
C Factor the matrix
C---------------------------------------------------------------------------
        ierror = dss_factor_real( handle, MKL_DSS_DEFAULTS, AA)
        IF (ierror .NE. MKL_DSS_SUCCESS ) GO TO 999        
        
400     CONTINUE
C---------------------------------------------------------------------------
C Get the solution vector
C---------------------------------------------------------------------------
      IF (IND.EQ.1) THEN
C     UPDATE IPFAC
          IPFAC = 1 !STATE THAT THE MEMORY IS ALREADY ALLOCATE
          CALL MINTFIL('BLOK',IPFAC ,1,5 ,1)  
      ENDIF
        
        if (IND.NE.2) GOTO 500
        
        ierror = dss_solve_real( handle, MKL_DSS_DEFAULTS, V, nRhs,
     &  X)
        V(1:NEQ) = X(1:NEQ)
        IF (ierror .NE. MKL_DSS_SUCCESS ) GO TO 999
        
500     CONTINUE
        
C---------------------------------------------------------------------------
C Deallocate solver storage
C---------------------------------------------------------------------------
      IF (IND.EQ.3) THEN
        ierror = dss_delete( handle, MKL_DSS_DEFAULTS )
        IF (ierror .NE. MKL_DSS_SUCCESS ) GO TO 999
      ENDIF
        
999     CONTINUE     
        
        
1000     DEALLOCATE(AA,MAXC,X)
        
      RETURN
      END
C	=============================================================
C	=============================================================
C	=============================================================
      
      
      
      SUBROUTINE PARDISI_WORK_WELL(MAXA,V,NEQ,IND,TYP1,TYP2)  
	IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C      INCLUDE 'MKL_PARDISO.F77'
      CHARACTER*4 TYP1,TYP2  

C     =============================================================        
      COMMON /FTIM/ TIM(20),IDATE,ITIME
C     STORE PARDISO POINTER NUMBER 
      COMMON /PARDIS/ LPTO(64) 
C     =============================================================   
C     .. Internal solver memory pointer for 64-bit architectures
C     .. INTEGER*8 LPT(64)
C     .. Internal solver memory pointer for 32-bit architectures
C     .. INTEGER*4 LPT(64)
C     .. This is OK in both cases
      ALLOCATABLE LPT(:)
C     ============================================================= 
      POINTER ::  AA(:),MAXC(:),IPERM(:),X(:)  
      DIMENSION MAXA(NEQ),V(NEQ),IPARM(64),IIPARM(64)
C     =============================================================   
      
      ALLOCATE(IPERM(NEQ),X(NEQ))
      
C     IPFAC DEFINE WHETHER THE MEMORY HAS BEEN ALLOCATED OR NOT
	CALL MINTFIL('BLOK',IPFAC ,1,5 ,0)  

      NWK = MAXA(NEQ+1) - 1
	CALL MINTFIL('BLOK',NBLOCK,1,1 ,0)
	CALL MINTFIL('BLOK',MSTOR ,1,2 ,0)

	ALLOCATE(AA(NWK),MAXC(NWK))
      
C     FOR CLEARING MEMORY PHASE NO NEED TO READ THE STIFFNESS MATRIX 
	IF(IND.EQ.3) GOTO 50
	
      IF(IND.EQ.1) THEN
	CALL MDMOVE(TYP1,TYP2)
	ENDIF

	CALL MCALFIL(KSREC,TYP2)
	REWIND(KSREC)
	CALL MCALFIL(NFLCH,'MAXC')
	REWIND(NFLCH)
	
	IAA = 1
	ICC = 1
	DO IBLO = 1,NBLOCK
	CALL MINTFIL('BLOC',NHIG,3,IBLO,0)
	READ(KSREC) AA(IAA:IAA+NHIG-1)
	READ(NFLCH) MAXC(ICC:ICC+NHIG-1)
	IAA = IAA + NHIG
	ICC = ICC + NHIG
      ENDDO
!      WRITE (300,"(E12.5)") AA
C..
C.. Set up PARDISO control parameter
C..
50    NRHS = 1
      MAXFCT = 1
      MNUM = 1
      
      DO I = 1,64
         IPARM(I) = 0
      ENDDO
      IPARM(1)  =  1 ! no solver default
      IPARM(2)  =  2 !3 ! 2 = fill-in reordering from METIS , 3 = parallel (OpenMP) nested dissection algorithm 
      IPARM(3)  =  1 !MKL_GET_MAX_THREADS() !1 ! numbers of processors
      IPARM(4)  =  0 !32 ! no iterative-direct algorithm
      IPARM(5)  =  0 ! no user fill-in reducing permutation
      IPARM(6)  =  0 ! =0 solution on the first n compoments of x
      IPARM(7)  =  0 ! not in use
      IPARM(8)  =  9 ! numbers of iterative refinement steps
      IPARM(9)  =  0 ! not in use
      IPARM(10) =  13 ! perturbe the pivot elements with 1E-13
      IPARM(11) =  1 ! use nonsymmetric permutation and scaling MPS
      IPARM(12) =  0 ! not in use
      IPARM(13) =  0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try IPARM(13) = 1 in case of inappropriate accuracy
      IPARM(14) =  0 ! Output: number of perturbed pivots
      IPARM(15) =  0 ! not in use
      IPARM(16) =  0 ! not in use
      IPARM(17) =  0 ! not in use
      IPARM(18) = -1 ! Output: number of nonzeros in the factor LU
      IPARM(19) = -1 ! Output: Mflops for LU factorization
      IPARM(20) =  0 ! Output: Numbers of CG Iterations
      IPARM(27) =  0 ! Output: Numbers of CG Iterations
      IPARM(60) =  1 ! 0 IN-CORE 1 IN&OUT 2 OUT-CORE
      IERROR    =  0 ! initialize error flag
      MSGLVL    =  0 ! print statistical information
      MTYPE     = -2 ! symmetric, indefinite

      ALLOCATE( LPT(64) )
      
      GOTO (100,200,300),IND
      
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
100   IF(IPFAC.EQ.1) GOTO 300

120   DO I = 1, 64
      LPT(I) = 0
      ENDDO

      CALL CPU_TIME (TIME1)
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
c      IPHASE = 11 ! only reordering 
      IPHASE = 12 ! reordering and factorization in the same time
      CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)

C.. Factorization.      
c      IPHASE = 22 ! only factorization
c      CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
c     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
      CALL CPU_TIME (TIME2)
      TIM(14) =  TIM(14) + (TIME2-TIME1)
      
      DO I = 1, 64
      LPTO(I) = LPT(I)
      ENDDO

C     UPDATE IPFAC
      IPFAC = 1 !STATE THAT THE MEMORY IS ALREADY ALLOCATE
      CALL MINTFIL('BLOK',IPFAC ,1,5 ,1)  

C     PIVOT FLAG  0=POSITIVE DEF  1=NEGATIVE DEF  ONLY APPEARED IN PARDISO,NEWSOLVER,LOADUP    
C     THIS WILL BE USE IN NONLINEAR ANALYIS   SUBROUTINE--LOADUP  	
      IPIVT = 0
      IF(IPARM(23).GT.0) IPIVT = 1 !OPARM(23) CONTAIN NUMBER OF NEGATIVE EIGENVALUE
	CALL MINTFIL('BLOK',IPIVT ,1,6 ,1)  

      GOTO 900
C     --------------------------------------------------------------
C     --------------------------------------------------------------  
200   CALL CPU_TIME (TIME1)
      DO I = 1, 64
      LPT(I) = LPTO(I)
      ENDDO
C.. Back substitution and iterative refinement
      IPARM(8) = 2 ! max numbers of iterative refinement steps
      IPHASE = 33 ! back substitution and iterative refinement
      CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
      CALL CPU_TIME (TIME2)
      TIM(15) =  TIM(15) + (TIME2-TIME1)
      
      V(1:NEQ) = X(1:NEQ)

      GOTO 900
C     --------------------------------------------------------------
C     --------------------------------------------------------------  
C.. Termination and release of memory
300   IF(IPFAC.EQ.0) GOTO 900
      DO I = 1, 64
      LPT(I) = LPTO(I)
      ENDDO
      IPHASE = -1 ! release internal memory
      CALL PARDISO (LPT,MAXFCT,MNUM,MTYPE,IPHASE,NEQ,AA,MAXA,MAXC,
     1 IPERM,NRHS,IPARM,MSGLVL,V,X,IERROR)
     
C     UPDATE IPFAC
      IPFAC = 0 !STATE THAT THE MEMORY IS ALREADY RELEASE
      CALL MINTFIL('BLOK',IPFAC ,1,5 ,1)  
      
      IF(IND.EQ.1) GOTO 120 !BACK TO FACTORIZATION AFTER RELEASE MEMORY
      GOTO 900
C     --------------------------------------------------------------
C     --------------------------------------------------------------  
     
900	DEALLOCATE(AA,MAXC,IPERM,X,LPT)
	
      RETURN
      END








